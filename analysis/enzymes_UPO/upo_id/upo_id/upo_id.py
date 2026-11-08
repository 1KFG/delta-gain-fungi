#!/usr/bin/env python3
"""upo_id: classify heme-thiolate peroxidase (PF01328) proteins into UPO classes.

Primary scheme: Kinner et al. 2021 Family I (short) / Family II (long).
Secondary scheme: UPObase (Faiza et al. 2019a,b) motif subfamilies I-V and Pog-a..g.

Subcommands
  build     build alignment, HMMs and position map from a reference FASTA
  classify  classify proteins in a FASTA (plain, .gz or .zst)

External tools (must be on PATH, e.g. `module load hmmer/3.4 mafft/7.505`):
  hmmsearch, hmmbuild, hmmalign, phmmer, mafft (build only)
Python: 3.11+, Biopython. pyarrow is optional (Parquet output, Parquet SP table).

All paths come from arguments. The script does not resolve its own location.
"""
import argparse
import csv
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict

VERSION = "0.1.0"

# ---------------------------------------------------------------------------
# Key sites, AaeUPO precursor numbering (UniProt B9W4V6).
# Mature numbering (literature, PDB 2YOR) = precursor - 43 (2YOR DBREF 4-328 = UNP 47-371).
# Sources: 2YOR LINK records (Cys36-Fe; Glu122, Gly123 O, Ser126 -> Mg), 2YOR SSBOND 278-319,
# Faiza 2019a text (Arg189-Glu196 acid-base pair; EGD = Glu122-Gly123-Asp124).
# ---------------------------------------------------------------------------
AAE_OFFSET = 43
SITES = OrderedDict([
    ("pcp_P1", 78), ("prox_cys", 79), ("pcp_P2", 80),
    ("sig_S", 137), ("sig_I", 138), ("sig_G", 139),
    ("mg_glu", 165), ("exd_mid", 166), ("exd_asp", 167), ("mg_ser", 169),
    ("sxxrxd_R", 172), ("sxxrxd_D", 174),
    ("cs_arg", 232), ("acid_base_glu", 239),
    ("cterm_cys1", 321), ("cterm_cys2", 362),
])
CORE = (60, 280)            # AaeUPO precursor span used for core coverage
RXE_WINDOW = (180, 262)     # AaeUPO precursor span searched for R-X(6,7)-E
LONG_MIN_LEN = 330          # calibrated on reference set (short refs <= 308, long refs >= 352)

# UPObase Table 1 (Faiza et al. 2019b, Database baz122). '-' separates ordered motifs; X = any.
UPOBASE_MOTIFS = OrderedDict([
    ("Ia",   "PCP-[NS]HG-SIG-HXXF-EGD-SXXRXD-RXXXXXXE-FXD-C-C"),
    ("Ib",   "PCP-[NS]HG-GVARPD-SIG-HXXF-EGD-SXXRXD-G[AVFY]NG-RXXXXXXE-FXD-RQP-C-RV[IV]P-C"),
    ("II",   "PCP-NHG-RGN-S[IL]G-VPPLPG-IDG-HGRF-EGD-SMTRXD-RXXXXXXE-TXXXXXXR"),
    ("IIIa", "PCP-NH[NG]-G[ML]G-SIG-E[GA]D-SXTRXD-GPXTG-RXXXXXXE-TGG-CXXXE"),
    ("IIIb", "PCP-NH[NG]-G[ML]G-SIG-E[GT]D-SXTRXD-RXXXXXXE-TXG-CXXXQ"),
    ("IV",   "PCP-N[HY][NG]-FXXXD-S[IL]G-CDA-HXXF-EGD-SLTRXD-RXXXXXXE-GAAXXXYE"),
    ("V",    "EDXXH-PCP-NHG-SIG-GXG-EGD-SVTRXD-RXXXXXXE"),
    ("Pog-a", "RGPCP-NTL[AT]N-PXXG-NXT-HXXL-EHD-RXD-PXXXFG"),
    ("Pog-b", "RXPCP-PRXG-[EQ]HD-S[FMV]T-RXD"),
    ("Pog-c", "RXPCP-NTLXN-PXXGR-EHD-S[ML]S-RXD-GWXP"),
    ("Pog-d", "RXPCP-E[IHF]D-GSLS-RXD-RIPY"),
    ("Pog-e", "RXPCP-NSLAN-PRXG-LIXGM-GLNL-HXLI-EHD-SLS-RXD"),
    ("Pog-f", "RXPCP-[EQ]HD-S[LM]S-RXD-DXXXFN-RXXR"),
])
UPOBASE_UPO_SUBFAMILIES = {"Ia", "Ib", "II", "IIIa", "IIIb", "IV", "V"}

FAMILY_HMMS = [("FamII", "famII.hmm"), ("FamI_UPO", "famI_upo.hmm"), ("CPO", "cpo.hmm")]


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
def open_text(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz") or path.endswith(".bgz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    if path.endswith(".zst"):
        try:
            import zstandard  # type: ignore
            fh = open(path, "rb")
            return io.TextIOWrapper(zstandard.ZstdDecompressor().stream_reader(fh))
        except ImportError:
            if not shutil.which("zstd"):
                sys.exit("ERROR: .zst input needs the zstandard module or the zstd binary")
            p = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE)
            return io.TextIOWrapper(p.stdout)
    return open(path)


def read_fasta(path):
    name, seq = None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0] if line[1:].split() else "", []
            elif line:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def clean_seq(s):
    s = s.upper().rstrip("*")
    return re.sub(r"[^A-Z]", "X", s.replace("*", "X"))


def write_fasta(recs, path):
    with open(path, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


def run(cmd, **kw):
    sys.stderr.write("[upo_id] " + " ".join(cmd) + "\n")
    r = subprocess.run(cmd, **kw)
    if r.returncode != 0:
        sys.exit(f"ERROR: command failed ({r.returncode}): {' '.join(cmd)}")
    return r


def need(tools):
    miss = [t for t in tools if not shutil.which(t)]
    if miss:
        sys.exit("ERROR: missing tools on PATH: " + ", ".join(miss) +
                 " (try: module load hmmer/3.4 mafft/7.505)")


def read_afa(path):
    d, n = OrderedDict(), None
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            n = line[1:].split()[0]
            d[n] = []
        elif n:
            d[n].append(line.strip())
    return OrderedDict((k, "".join(v).upper()) for k, v in d.items())


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------
def cmd_build(a):
    need(["mafft", "hmmbuild"])
    os.makedirs(a.out, exist_ok=True)
    train = {r["name"]: r["train_set"] for r in csv.DictReader(open(a.train), delimiter="\t")}
    exclude = set(a.exclude or [])
    refs = OrderedDict((n, clean_seq(s)) for n, s in read_fasta(a.refs))
    pos = [(n, s) for n, s in refs.items()
           if train.get(n) in ("FamII", "FamI_UPO", "CPO", "holdout") and n not in exclude]
    if a.anchor not in dict(pos):
        sys.exit(f"ERROR: anchor {a.anchor} not in positive references")
    tmp = os.path.join(a.out, "_build_tmp")
    os.makedirs(tmp, exist_ok=True)
    write_fasta(pos, os.path.join(tmp, "pos.fa"))
    aln_path = os.path.join(a.out, "ref_alignment.afa")
    with open(aln_path, "w") as fh:
        run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", str(a.cpu), "--quiet",
             os.path.join(tmp, "pos.fa")], stdout=fh)
    aln = read_afa(aln_path)
    anchor = aln[a.anchor]
    ncol = len(anchor)
    # anchor HMM: match states = anchor (AaeUPO) residues, so state k == AaeUPO precursor residue k
    rf = "".join("x" if c != "-" else "." for c in anchor)
    sto = os.path.join(tmp, "anchor.sto")
    with open(sto, "w") as fh:
        fh.write("# STOCKHOLM 1.0\n")
        for n, s in aln.items():
            fh.write(f"{n} {s}\n")
        fh.write(f"#=GC RF {rf}\n//\n")
    run(["hmmbuild", "--hand", "--amino", "-n", "UPO_anchor", os.path.join(a.out, "anchor.hmm"), sto],
        stdout=subprocess.DEVNULL)
    # position map: AaeUPO precursor position -> residue/position in other characterized refs
    colpos = {}
    for n, s in aln.items():
        i, m = 0, {}
        for c, ch in enumerate(s):
            if ch != "-":
                i += 1
                m[c] = (i, ch)
        colpos[n] = m
    mapped = [n for n in a.map_refs if n in aln]
    with open(os.path.join(a.out, "posmap.tsv"), "w") as fh:
        fh.write("aae_precursor\taae_mature\taae_res\t" + "\t".join(f"{n}_pos\t{n}_res" for n in mapped) + "\n")
        k = 0
        for c in range(ncol):
            if anchor[c] == "-":
                continue
            k += 1
            row = [str(k), str(k - AAE_OFFSET), anchor[c]]
            for n in mapped:
                p = colpos[n].get(c)
                row += [str(p[0]), p[1]] if p else ["", "-"]
            fh.write("\t".join(row) + "\n")
    # family HMMs from subsets of the reference alignment
    for fam, fn in FAMILY_HMMS:
        names = [n for n in aln if train.get(n) == fam]
        if len(names) < 2:
            sys.exit(f"ERROR: fewer than 2 sequences for {fam}")
        cols = [c for c in range(ncol) if any(aln[n][c] != "-" for n in names)]
        sub = os.path.join(tmp, f"{fam}.afa")
        with open(sub, "w") as fh:
            for n in names:
                fh.write(f">{n}\n{''.join(aln[n][c] for c in cols)}\n")
        run(["hmmbuild", "--amino", "-n", fam, os.path.join(a.out, fn), sub], stdout=subprocess.DEVNULL)
        shutil.copy(sub, os.path.join(a.out, f"{fam}_train.afa"))
    shutil.copy(a.pfam_hmm, os.path.join(a.out, "PF01328.hmm"))
    # characterized references (nearest-reference step) and phmmer seeds
    roles = {r["name"]: r for r in csv.DictReader(open(a.manifest), delimiter="\t")}
    char = [(n, refs[n]) for n in refs if roles.get(n, {}).get("role") == "characterized" and n not in exclude]
    seeds = [(n, refs[n]) for n in refs if n not in exclude and
             (roles.get(n, {}).get("role") == "upobase_seed_UPO" or n in ("AaeUPO", "MroUPO", "LfuCPO"))]
    write_fasta(char, os.path.join(a.out, "characterized.fa"))
    write_fasta(seeds, os.path.join(a.out, "phmmer_seeds.fa"))
    meta = {"version": VERSION, "anchor": a.anchor, "excluded": sorted(exclude),
            "n_aligned": len(aln), "train_counts": {f: sum(1 for n in aln if train.get(n) == f)
                                                    for f, _ in FAMILY_HMMS},
            "characterized": [n for n, _ in char], "n_seeds": len(seeds)}
    json.dump(meta, open(os.path.join(a.out, "models.json"), "w"), indent=1)
    shutil.rmtree(tmp)
    sys.stderr.write(json.dumps(meta) + "\n")


# ---------------------------------------------------------------------------
# classify helpers
# ---------------------------------------------------------------------------
def parse_tblout(path):
    """full-sequence E-value and score per target from hmmsearch/phmmer --tblout."""
    best = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        t, q, ev, sc = f[0], f[2], float(f[4]), float(f[5])
        if t not in best or sc > best[t][2]:
            best[t] = (q, ev, sc)
    return best


def parse_domtbl(path):
    out = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        t, ev_full, sc_full = f[0], float(f[6]), float(f[7])
        hf, ht, af, at = int(f[15]), int(f[16]), int(f[17]), int(f[18])
        d = out.setdefault(t, {"evalue": ev_full, "score": sc_full, "ndom": 0, "hmm_cov": set()})
        d["ndom"] += 1
        d["hmm_cov"].update(range(hf, ht + 1))
    return out


def parse_a2m(path):
    """hmmalign A2M: upper/'-' = match state, lower = insert. Returns name -> list of
    (state_k or None, residue or '-', query_pos or None)."""
    out = OrderedDict()
    for n, s in read_afa_raw(path).items():
        k, q, ev = 0, 0, []
        for ch in s:
            if ch == "-":
                k += 1
                ev.append((k, "-", None))
            elif ch.isupper():
                k += 1
                q += 1
                ev.append((k, ch, q))
            elif ch.islower():
                q += 1
                ev.append((None, ch.upper(), q))
            # '.' (A2M pads) ignored
        out[n] = ev
    return out


def read_afa_raw(path):
    d, n = OrderedDict(), None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            n = line[1:].split()[0]
            d[n] = []
        elif n:
            d[n].append(line.strip())
    return OrderedDict((k, "".join(v)) for k, v in d.items())


def motif_regex(m):
    return re.compile(m.replace("X", "."))


UPOBASE_RX = OrderedDict((k, [motif_regex(x) for x in v.split("-")]) for k, v in UPOBASE_MOTIFS.items())


def ordered_motif_match(seq, rxs):
    """Earliest-match chaining of ordered motifs (no gap limits). Returns positions or None."""
    pos, start = [], 0
    for rx in rxs:
        m = rx.search(seq, start)
        if not m:
            return None
        pos.append(m.start() + 1)
        start = m.end()
    return pos


def upobase_classify(seq, has_cys):
    hits = []
    for name, rxs in UPOBASE_RX.items():
        p = ordered_motif_match(seq, rxs)
        if p:
            hits.append((name, len(rxs)))
    if hits:
        best = max(hits, key=lambda x: x[1])[0]   # most motif elements wins; ties -> table order
    elif has_cys:
        best = "Pog-g"
    else:
        best = "none"
    group = "UPO_subfamily" if best in UPOBASE_UPO_SUBFAMILIES else ("Pog" if best.startswith("Pog") else "none")
    return best, group, ",".join(h[0] for h in hits)


def load_sp_ids(path, wanted):
    if path is None:
        return None
    ids = set()
    if path.endswith(".parquet"):
        import pyarrow.dataset as ds  # type: ignore
        import pyarrow.compute as pc  # type: ignore
        d = ds.dataset(path)
        t = d.to_table(columns=["protein_id"], filter=pc.field("protein_id").isin(list(wanted)))
        ids = set(t.column("protein_id").to_pylist())
    else:
        with open_text(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f and f[0] and f[0] != "protein_id":
                    ids.add(f[0])
    return ids


def nearest_reference(seq, refs, aligner):
    best = ("", 0.0, 0.0, 0.0)
    for rn, rs in refs:
        try:
            aln = aligner.align(seq, rs)[0]
        except Exception:
            continue
        c = aln.counts()
        span = c.identities + c.mismatches + c.gaps
        if span == 0:
            continue
        pid = 100.0 * c.identities / span
        qcov = 100.0 * (c.identities + c.mismatches) / len(seq)
        if aln.score > best[3]:
            best = (rn, pid, qcov, aln.score)
    return best


def decide(r):
    """Decision table. Returns (final_class, kinner_family, confidence, reason)."""
    if not r["htp_hit"]:
        return "no_HTP_hit", "", "", "no PF01328 (GA) or family-HMM hit"
    fails = []
    if r["core_cov"] < 0.5:
        fails.append(f"core coverage {r['core_cov']:.2f} < 0.5")
    if r["prox_cys"] != "C":
        fails.append(f"no proximal Cys at AaeUPO C36 (found {r['prox_cys']})")
    if r["mg_glu"] != "E":
        fails.append(f"no Mg-site Glu at AaeUPO E122 (found {r['mg_glu']})")
    if not r["acid_base_ok"]:
        fails.append("no acid-base Glu (AaeUPO E196 column or R-X(6,7)-E window)")
    if fails:
        return "HTP_fails_motif", "", "", "; ".join(fails)
    fam = r["best_family_hmm"]
    cls = {"FamII": "UPO_FamilyII", "FamI_UPO": "UPO_FamilyI", "CPO": "CPO_like_FamilyI"}[fam]
    kin = "II" if fam == "FamII" else "I"
    # supporting evidence
    motif_fam = "II" if r["exd_mid"] in "GATS" and r["exd_mid"] != "" else ("I" if r["exd_mid"] == "H" else "?")
    len_fam = "II" if r["length"] >= LONG_MIN_LEN else "I"
    ss_fam = "II" if r["cterm_ss_pair"] else "I"
    disc = []
    if motif_fam != kin:
        disc.append(f"metal-site E{r['exd_mid']}D suggests Family {motif_fam}")
    if len_fam != kin:
        disc.append(f"length {r['length']} suggests Family {len_fam}")
    if ss_fam != kin:
        disc.append(("C-terminal Cys pair present" if r["cterm_ss_pair"] else "no C-terminal Cys pair") +
                    f" suggests Family {ss_fam}")
    if r["family_margin"] < 20:
        disc.append(f"family HMM margin {r['family_margin']:.1f} bits < 20")
    conf = "high" if not disc else ("medium" if len(disc) == 1 else "low")
    reason = f"best family HMM {fam}; Cys/Mg-Glu/acid-base Glu present"
    if disc:
        reason += "; discordant: " + "; ".join(disc)
    return cls, kin, conf, reason


# ---------------------------------------------------------------------------
# classify
# ---------------------------------------------------------------------------
def cmd_classify(a):
    need(["hmmsearch", "hmmalign"] + ([] if a.no_phmmer else ["phmmer"]))
    from Bio import Align
    from Bio.Align import substitution_matrices
    md = a.models
    tmpbase = a.tmpdir or os.environ.get("SCRATCH") or os.environ.get("TMPDIR") or None
    tmp = tempfile.mkdtemp(prefix="upo_id_", dir=tmpbase)
    seqs = OrderedDict()
    for n, s in read_fasta(a.input):
        if n in seqs:
            sys.exit(f"ERROR: duplicate sequence ID {n}")
        seqs[n] = clean_seq(s)
    if not seqs:
        sys.exit("ERROR: no sequences read")
    allfa = os.path.join(tmp, "in.fa")
    write_fasta(seqs.items(), allfa)
    cpu = str(a.cpu)
    # (a) HMM detection
    pf_dom = os.path.join(tmp, "pf.dom")
    run(["hmmsearch", "--cut_ga", "--cpu", cpu, "--noali", "--domtblout", pf_dom, "-o", os.devnull,
         os.path.join(md, "PF01328.hmm"), allfa])
    pf = parse_domtbl(pf_dom)
    fam = {}
    for f, fn in FAMILY_HMMS:
        tb = os.path.join(tmp, f + ".tbl")
        run(["hmmsearch", "-E", "1", "--cpu", cpu, "--noali", "--tblout", tb, "-o", os.devnull,
             os.path.join(md, fn), allfa])
        fam[f] = parse_tblout(tb)
    cand = [n for n in seqs if n in pf or any(n in fam[f] and fam[f][n][1] <= a.evalue for f, _ in FAMILY_HMMS)]
    sys.stderr.write(f"[upo_id] {len(seqs)} proteins; {len(cand)} HTP candidates\n")
    # (b) alignment to anchor model
    aligned = {}
    if cand:
        cfa = os.path.join(tmp, "cand.fa")
        write_fasta([(n, seqs[n]) for n in cand], cfa)
        a2m = os.path.join(tmp, "cand.a2m")
        run(["hmmalign", "--outformat", "A2M", "-o", a2m, os.path.join(md, "anchor.hmm"), cfa])
        aligned = parse_a2m(a2m)
    # optional UPObase-style phmmer detection with the published seeds
    seedhit = {}
    if cand and not a.no_phmmer:
        tb = os.path.join(tmp, "phmmer.tbl")
        run(["phmmer", "-E", "10", "--incE", "0.01", "--cpu", cpu, "--noali", "--tblout", tb, "-o", os.devnull,
             os.path.join(md, "phmmer_seeds.fa"), os.path.join(tmp, "cand.fa")])
        # phmmer tblout: target = candidate, query = seed
        for line in open(tb):
            if line.startswith("#"):
                continue
            f = line.split()
            t, q, ev, sc = f[0], f[2], float(f[4]), float(f[5])
            if t not in seedhit or ev < seedhit[t][1]:
                seedhit[t] = (q, ev, sc)
    # nearest characterized reference
    refs = list(read_fasta(os.path.join(md, "characterized.fa")))
    aligner = Align.PairwiseAligner(mode="local", open_gap_score=-11, extend_gap_score=-1)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    sp_ids = load_sp_ids(a.sp_table, set(cand) if not a.all else set(seqs))
    rows = []
    ids = list(seqs) if a.all else cand
    for n in ids:
        s = seqs[n]
        r = OrderedDict()
        r["protein_id"] = n
        r["length"] = len(s)
        r["mw_kda"] = round(len(s) * 0.110, 1)
        p = pf.get(n)
        r["pf01328_ga_hit"] = p is not None
        r["pf01328_evalue"] = p["evalue"] if p else None
        r["pf01328_score"] = p["score"] if p else None
        r["pf01328_ndom"] = p["ndom"] if p else 0
        r["pf01328_hmm_cov"] = round(len(p["hmm_cov"]) / 186, 3) if p else 0.0
        scores = {}
        for f, _ in FAMILY_HMMS:
            h = fam[f].get(n)
            r[f"{f}_evalue"] = h[1] if h else None
            r[f"{f}_score"] = h[2] if h else None
            scores[f] = h[2] if h else -1e9
        order = sorted(scores, key=lambda f: -scores[f])
        r["htp_hit"] = n in cand
        r["best_family_hmm"] = order[0] if scores[order[0]] > -1e9 else ""
        sec = scores[order[1]] if scores[order[1]] > -1e9 else 0.0
        r["family_margin"] = round(scores[order[0]] - sec, 1) if r["best_family_hmm"] else 0.0
        ev = aligned.get(n)
        state = {}
        if ev:
            for k, res, q in ev:
                if k is not None:
                    state[k] = (res, q)
        for key, k in SITES.items():
            r[key] = state.get(k, ("-", None))[0] if ev else ""
        r["pcp_triad"] = (r["pcp_P1"] + r["prox_cys"] + r["pcp_P2"]) if ev else ""
        r["prox_cys_qpos"] = state.get(SITES["prox_cys"], ("", None))[1] if ev else None
        r["exd_motif"] = (r["mg_glu"] + r["exd_mid"] + r["exd_asp"]) if ev else ""
        r["sig_motif"] = (r["sig_S"] + r["sig_I"] + r["sig_G"]) if ev else ""
        r["acid_base_qpos"] = state.get(SITES["acid_base_glu"], ("", None))[1] if ev else None
        # R-X(6,7)-E in the query segment aligned to the AaeUPO acid-base window
        rxe = ""
        if ev:
            # query span from the first to the last residue matched to window states (inserts included)
            win_q = [q for k, res, q in ev if k is not None and q is not None and RXE_WINDOW[0] <= k <= RXE_WINDOW[1]]
            if win_q:
                lo, hi = min(win_q), max(win_q)
                seg = s[lo - 1:hi]
                m = re.search(r"R.{6,7}E", seg)
                if m:
                    rxe = f"{m.group()}@{lo + m.start()}"
        r["rxe_motif"] = rxe
        strict = r["acid_base_glu"] == "E" if ev else False
        r["acid_base_ok"] = bool(strict or rxe)
        r["acid_base_basis"] = "E196_column" if strict else ("RxE_window" if rxe else "")
        r["cterm_ss_pair"] = (r["cterm_cys1"] == "C" and r["cterm_cys2"] == "C") if ev else False
        core = [k for k in range(CORE[0], CORE[1] + 1)]
        r["core_cov"] = round(sum(1 for k in core if state.get(k, ("-",))[0] != "-") / len(core), 3) if ev else 0.0
        r["length_class"] = "long" if len(s) >= LONG_MIN_LEN else "short"
        r["signal_peptide"] = ("yes" if n in sp_ids else "no") if sp_ids is not None else "NA"
        # UPObase motif label from the alignment (Faiza 2019b intro: PCP-EGD-R-E vs PCP-EHD-E)
        if ev and r["prox_cys"] == "C" and r["exd_mid"] == "G" and ((r["cs_arg"] == "R" and strict) or rxe):
            lab = "PCP-EGD-R-E"
        elif ev and r["prox_cys"] == "C" and r["exd_mid"] == "H" and r["acid_base_ok"]:
            lab = "PCP-EHD-E"
        else:
            lab = "neither"
        r["upobase_catalytic_label"] = lab
        sub, grp, allm = upobase_classify(s, "PCP" in s)
        r["upobase_subfamily"] = sub
        r["upobase_group"] = grp
        r["upobase_all_matches"] = allm
        sh = seedhit.get(n)
        r["phmmer_best_seed"] = sh[0] if sh else ""
        r["phmmer_seed_evalue"] = sh[1] if sh else None
        r["phmmer_included"] = bool(sh and sh[1] <= 0.01)
        if r["htp_hit"]:
            rn, pid, qcov, _ = nearest_reference(s, refs, aligner)
        else:
            rn, pid, qcov = "", 0.0, 0.0
        r["nearest_characterized"] = rn
        r["nearest_pct_identity"] = round(pid, 1)
        r["nearest_query_cov"] = round(qcov, 1)
        cls, kin, conf, reason = decide(r)
        r["kinner_family"] = kin
        r["final_class"] = cls
        r["confidence"] = conf
        r["reason"] = reason
        rows.append(r)
    write_rows(rows, a.out, a.parquet)
    if not a.keep_tmp:
        shutil.rmtree(tmp)
    from collections import Counter
    c = Counter(r["final_class"] for r in rows)
    sys.stderr.write("[upo_id] classes: " + json.dumps(dict(c)) + "\n")


def fmt(v):
    if v is None:
        return "NA"
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, float):
        return f"{v:.3g}" if (v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e6)) else (f"{v:.4g}")
    return str(v)


def write_rows(rows, out, parquet):
    if not rows:
        sys.stderr.write("[upo_id] no rows to write\n")
    cols = list(rows[0].keys()) if rows else ["protein_id"]
    text = "\t".join(cols) + "\n" + "".join("\t".join(fmt(r[c]) for c in cols) + "\n" for r in rows)
    if out.endswith(".zst"):
        subprocess.run(["zstd", "-q", "-f", "-o", out], input=text.encode(), check=True)
    elif out.endswith(".gz"):
        with gzip.open(out, "wt") as fh:
            fh.write(text)
    elif out == "-":
        sys.stdout.write(text)
    else:
        with open(out, "w") as fh:
            fh.write(text)
    if parquet:
        import pyarrow as pa  # type: ignore
        import pyarrow.parquet as pq  # type: ignore
        pq.write_table(pa.Table.from_pylist(rows), parquet, compression="zstd")


def main():
    ap = argparse.ArgumentParser(prog="upo_id", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=VERSION)
    sp = ap.add_subparsers(dest="cmd", required=True)
    b = sp.add_parser("build", help="build models from references")
    b.add_argument("--refs", required=True, help="reference FASTA (ref/refs_all.fa)")
    b.add_argument("--train", required=True, help="train_sets.tsv (name, train_set)")
    b.add_argument("--manifest", required=True, help="ref/manifest.tsv (name, role)")
    b.add_argument("--pfam-hmm", required=True, help="PF01328 HMM (hmmfetch from Pfam-A.hmm)")
    b.add_argument("--out", required=True, help="model directory")
    b.add_argument("--anchor", default="AaeUPO")
    b.add_argument("--map-refs", nargs="*", default=["LfuCPO", "MroUPO", "HspUPO", "CviUPO"])
    b.add_argument("--exclude", nargs="*", help="reference names to leave out (leave-one-out tests)")
    b.add_argument("--cpu", type=int, default=4)
    c = sp.add_parser("classify", help="classify proteins")
    c.add_argument("--input", "-i", required=True, help="protein FASTA (.fa, .gz, .zst)")
    c.add_argument("--models", "-m", required=True, help="model directory from 'build'")
    c.add_argument("--out", "-o", required=True, help="output TSV (.tsv, .tsv.gz, .tsv.zst or -)")
    c.add_argument("--parquet", help="also write Parquet here (needs pyarrow)")
    c.add_argument("--sp-table", help="signal peptide positives: TSV/list (first column protein_id) "
                                      "or Parquet with a protein_id column (e.g. BFD signalp table)")
    c.add_argument("--evalue", type=float, default=1e-5, help="family-HMM full-sequence E-value for candidates")
    c.add_argument("--all", action="store_true", help="write a row for every input protein")
    c.add_argument("--no-phmmer", action="store_true", help="skip phmmer seed search")
    c.add_argument("--cpu", type=int, default=2)
    c.add_argument("--tmpdir", help="temporary directory (default $SCRATCH, then $TMPDIR)")
    c.add_argument("--keep-tmp", action="store_true")
    a = ap.parse_args()
    {"build": cmd_build, "classify": cmd_classify}[a.cmd](a)


if __name__ == "__main__":
    main()
