#!/usr/bin/env python3.12
"""
upo_sites -- find fungal unspecific peroxygenases (PF01328), annotate their active-site
residues in HspUPO numbering, assign them to Clades I-V, and score the aliphatic pocket.

WHAT THIS TOOL DOES AND DOES NOT CLAIM
--------------------------------------
It reports MEASURED SEQUENCE FEATURES and a CLADE ASSIGNMENT. It does not predict activity.

The calibration set behind it is 19 enzymes (Barros-Reguera et al. 2025,
DOI 10.1186/s13068-025-02675-w). In that set:

  * the aliphatic pocket does NOT predict aromatic ring hydroxylation (Pearson r = 0.18);
  * the aliphatic pocket DOES correlate with peroxygenase-vs-peroxidase propensity
    (r = 0.75), but it is nearly collinear with clade, and the partial correlation given
    clade is only r = 0.13.

So clade carries the signal, not the pocket. This tool therefore emits clade (with a
confidence), plus the raw site residues as evidence. It deliberately does NOT emit a
predicted activity number. See CALIBRATION_MATRIX.md sections 4, 5 and 8.

METHOD
------
1. hmmsearch the PF01328 profile against the input proteome(s) at the gathering threshold.
2. hmmalign the hits, in chunks, together with the 19-enzyme calibration backbone.
3. Read residues at fixed HSPUPO positions. Position -> HMM match-state index is computed
   once from the profile, so it is identical in every chunk and at any scale.
4. Assign clade by nearest neighbour to the backbone over match-state columns, reporting
   identity and the margin to the best member of any other clade.

Only HMM MATCH STATES are comparable across sequences. Three of the paper's eight aliphatic
pocket positions (223, 228, 232) are not usable: 228 falls in an insert column of PF01328,
and 223/232 are the least reliably aligned columns in the calibration set. The pocket score
is therefore over the five positions that are reliable: 73, 77, 81, 183, 187.
"""

from __future__ import annotations
import argparse, gzip, os, shutil, subprocess, sys, tempfile
from pathlib import Path

# ---------------------------------------------------------------- site definitions
# HspUPO (UniProt A0A1Y2TH07 / NCBI OTA57433.1, PDB 7O1R), 261 aa.
# Locations are the paper's own labels from Figure S2, defined on Clade V.
SITES: dict[int, str] = {
    39: "Cys_proximal",
    73: "AP", 77: "AP", 80: "E", 81: "AP", 84: "E",
    108: "CP", 109: "E_of_EHD", 110: "Base_EHD", 111: "D_of_EHD",
    172: "E", 173: "base_long", 175: "E", 176: "CP", 179: "CP",
    180: "Acid", 183: "AP", 187: "AP",
    223: "AP/E", 228: "AP", 232: "AP",
}
POCKET5 = (73, 77, 81, 183, 187)      # alignment-reliable aliphatic pocket
ALIPHATIC = frozenset("ILVAM")
CATALYTIC = (39, 109, 110, 111, 180)

CLADE_T1 = {  # in-sample T1 (peroxygenase %) from Figure S8; n=19; NOT a prediction
    "I": "0-53 (bimodal, n=5)", "II": "0 (n=2)", "III": "0-3 (n=2)",
    "IV": "0-5 (n=1)", "V": "33-82 (n=9)",
}


def log(msg: str) -> None:
    print(f"[upo_sites] {msg}", file=sys.stderr, flush=True)


def need(tool: str) -> str:
    p = shutil.which(tool)
    if not p:
        sys.exit(f"ERROR: {tool} not on PATH. On UCR HPCC: module load hmmer/3.4")
    return p


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    if path.endswith(".zst"):
        return subprocess.Popen(["zstd", "-dcq", path], stdout=subprocess.PIPE,
                                text=True).stdout
    return open(path, "rt")


def read_fasta(path: str):
    name, buf = None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            elif line:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


# ---------------------------------------------------------------- stockholm
def parse_stockholm(path: str) -> tuple[dict[str, str], str]:
    """Return {seqname: aligned}, RF line. Handles interleaved blocks."""
    seqs: dict[str, list[str]] = {}
    rf: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line == "//":
                continue
            if line.startswith("#=GC RF"):
                rf.append(line.split()[-1])
                continue
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            seqs.setdefault(parts[0], []).append(parts[-1])
    return {k: "".join(v) for k, v in seqs.items()}, "".join(rf)


def match_column_indices(rf: str) -> list[int]:
    """Columns of the alignment that are HMM match states (RF marks them non-'.')."""
    return [i for i, c in enumerate(rf) if c != "." and c != "-"]


def run_hmmalign(hmm: str, fasta: str, out_sto: str, trim: bool = False) -> None:
    cmd = [need("hmmalign"), "--amino", "-o", out_sto, hmm, fasta]
    if trim:
        cmd.insert(2, "--trim")
    subprocess.run(cmd, check=True)


# ---------------------------------------------------------------- position map
def build_position_map(hmm: str, anchor_fa: str, anchor_name: str,
                       positions) -> dict[int, int]:
    """
    Map anchor (HspUPO) 1-based residue number -> match-state ordinal (0-based).
    Computed once from the profile alone, so every chunk uses the same mapping.
    Positions that land in insert columns are omitted and reported.
    """
    with tempfile.TemporaryDirectory() as td:
        sto = os.path.join(td, "anchor.sto")
        run_hmmalign(hmm, anchor_fa, sto)
        seqs, rf = parse_stockholm(sto)
        key = next((k for k in seqs if k.split("|")[0] == anchor_name), None)
        if key is None:
            sys.exit(f"ERROR: anchor {anchor_name} not found in {anchor_fa}")
        aln = seqs[key]
        if not rf:
            sys.exit("ERROR: no #=GC RF line in hmmalign output; cannot identify match states")
        mcols = match_column_indices(rf)
        col_to_m = {c: i for i, c in enumerate(mcols)}
        pmap, residue, skipped = {}, 0, []
        for col, ch in enumerate(aln):
            if ch.isalpha():
                residue += 1
                if residue in positions:
                    if col in col_to_m:
                        pmap[residue] = col_to_m[col]
                    else:
                        skipped.append(residue)
        if skipped:
            log(f"positions in INSERT columns, not scoreable across sequences: {sorted(skipped)}")
        return pmap


# ---------------------------------------------------------------- search
def hmmsearch_ids(hmm: str, proteome: str, workdir: str, cpu: int,
                  use_ga: bool, evalue: float) -> dict[str, float]:
    tbl = os.path.join(workdir, "hits.tbl")
    cmd = [need("hmmsearch"), "--cpu", str(cpu), "--tblout", tbl, "-o", os.devnull]
    cmd += ["--cut_ga"] if use_ga else ["-E", str(evalue)]
    src = proteome
    tmp_plain = None
    if proteome.endswith((".gz", ".zst")):
        tmp_plain = os.path.join(workdir, "proteome.fa")
        with open(tmp_plain, "w") as out:
            for n, s in read_fasta(proteome):
                out.write(f">{n}\n{s}\n")
        src = tmp_plain
    subprocess.run(cmd + [hmm, src], check=True)
    best: dict[str, float] = {}
    with open(tbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) > 4:
                tid, ev = f[0], float(f[4])
                if tid not in best or ev < best[tid]:
                    best[tid] = ev
    if tmp_plain:
        os.unlink(tmp_plain)
    return best


# ---------------------------------------------------------------- scoring
def load_backbone(path: str):
    """Backbone headers look like Name|cladeV|Short|ACC."""
    out = []
    for name, seq in read_fasta(path):
        f = name.split("|")
        clade = f[1].replace("clade", "") if len(f) > 1 else "?"
        # hmmalign takes UNALIGNED sequence. Accept either the raw backbone or an aligned one
        # (.afa): degapping here is always correct and is the single point where backbone
        # sequences enter, so both the anchor map and every alignment chunk stay valid.
        out.append((f[0], clade, seq.replace("-", "").replace(".", "")))
    return out


def identity(a: str, b: str) -> float:
    n = same = 0
    for x, y in zip(a, b):
        if x == "-" and y == "-":
            continue
        n += 1
        if x == y and x != "-":
            same += 1
    return 100.0 * same / n if n else 0.0


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Annotate PF01328 active sites and assign UPO clades.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--proteome", required=True, nargs="+",
                    help="protein FASTA (.fa/.gz/.zst), one or more")
    ap.add_argument("--hmm", required=True, help="PF01328.hmm")
    ap.add_argument("--backbone", required=True,
                    help="19-enzyme calibration FASTA (data/clade_backbone/raw.fa)")
    ap.add_argument("--anchor", default="HspUPO", help="anchor sequence name in backbone")
    ap.add_argument("--out", required=True, help="output TSV (.gz/.zst ok)")
    ap.add_argument("--chunk", type=int, default=4000, help="sequences per hmmalign chunk")
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--evalue", type=float, default=1e-5,
                    help="used only when --no-ga is given")
    ap.add_argument("--no-ga", action="store_true",
                    help="use an E-value cutoff instead of the Pfam gathering threshold")
    ap.add_argument("--validated-identity", type=float, default=52.7,
                    help="lowest identity at which the clade caller was actually validated "
                         "(leave-one-out minimum over the 19 calibration enzymes). Calls at or "
                         "above this are flagged 'validated_range'.")
    ap.add_argument("--unassign-below", type=float, default=30.0,
                    help="below this identity to the nearest backbone member the clade is "
                         "reported as 'unassigned' rather than forced to a call")
    ap.add_argument("--min-margin", type=float, default=3.0,
                    help="identity margin over the best member of any OTHER clade. NOTE: in the "
                         "n=19 calibration set margin did NOT separate correct from incorrect "
                         "calls, so this is a heuristic, not a calibrated threshold.")
    ap.add_argument("--pre-selected", action="store_true",
                    help="input sequences are already known PF01328 members (e.g. pulled from "
                         "the BFD pfam table); skip hmmsearch and annotate all of them")
    ap.add_argument("--keep-temp", action="store_true")
    args = ap.parse_args()

    need("hmmalign")
    if not args.pre_selected:
        need("hmmsearch")
    backbone = load_backbone(args.backbone)
    bb_names = {n for n, _, _ in backbone}
    log(f"backbone: {len(backbone)} calibration sequences, "
        f"clades {sorted({c for _, c, _ in backbone})}")

    workdir = tempfile.mkdtemp(prefix="upo_sites.")
    try:
        # anchor-only fasta for the position map
        anchor_fa = os.path.join(workdir, "anchor.fa")
        with open(anchor_fa, "w") as fh:
            for n, c, s in backbone:
                if n == args.anchor:
                    fh.write(f">{n}\n{s}\n")   # already degapped by load_backbone
                    break
            else:
                sys.exit(f"ERROR: anchor {args.anchor} not in backbone")
        pmap = build_position_map(args.hmm, anchor_fa, args.anchor, set(SITES))
        usable = sorted(pmap)
        log(f"scoreable positions ({len(usable)}/{len(SITES)}): {usable}")
        pocket_use = [p for p in POCKET5 if p in pmap]
        if len(pocket_use) != len(POCKET5):
            log(f"WARNING: pocket positions missing from map: "
                f"{sorted(set(POCKET5) - set(pocket_use))}")

        # search every proteome
        hits: list[tuple[str, str, str, float]] = []   # (sample, id, seq, evalue)
        for prot in args.proteome:
            sample = Path(prot).name
            for suf in (".gz", ".zst", ".fa", ".faa", ".fasta", ".proteins"):
                sample = sample[: -len(suf)] if sample.endswith(suf) else sample
            if args.pre_selected:
                got = 0
                for n, s in read_fasta(prot):
                    hits.append((sample, n.split()[0], s, float("nan"))); got += 1
                log(f"{sample}: {got} pre-selected sequences")
                continue
            ev = hmmsearch_ids(args.hmm, prot, workdir, args.cpu,
                               not args.no_ga, args.evalue)
            if not ev:
                log(f"{sample}: 0 PF01328 hits")
                continue
            got = 0
            for n, s in read_fasta(prot):
                sid = n.split()[0]
                if sid in ev:
                    hits.append((sample, sid, s, ev[sid])); got += 1
            log(f"{sample}: {got} PF01328 hits")
        if not hits:
            sys.exit("No PF01328 hits in any input. Nothing to do.")
        log(f"total hits: {len(hits)}")

        # output
        opener = (lambda p: gzip.open(p, "wt")) if args.out.endswith(".gz") else open
        zstd_proc = None
        if args.out.endswith(".zst"):
            zstd_proc = subprocess.Popen(["zstd", "-qf", "-o", args.out],
                                         stdin=subprocess.PIPE, text=True)
            out = zstd_proc.stdin
        else:
            out = opener(args.out, "w") if args.out.endswith(".gz") else open(args.out, "w")

        cols = (["sample", "protein_id", "hmm_evalue", "length",
                 "clade", "clade_identity", "clade_margin", "clade_confidence",
                 "nearest_ref", "pocket_score", "pocket_max", "pocket_residues",
                 "catalytic_ok", "catalytic_detail", "sites_mapped", "clade_T1_insample"]
                + [f"p{p}_{SITES[p]}" for p in usable])
        out.write("# upo_sites: measured sequence features and clade assignment.\n")
        out.write("# NOT an activity prediction. See CALIBRATION_MATRIX.md sections 4, 5, 8.\n")
        out.write("# Residue numbering is HspUPO (PDB 7O1R). pocket_score is over "
                  f"{pocket_use} only.\n")
        out.write("# clade_T1_insample is the OBSERVED range in the 19-enzyme calibration "
                  "set, not a prediction for this protein.\n")
        out.write(f"# clade_confidence tiers: validated_range = identity >= "
                  f"{args.validated_identity} (the leave-one-out validated range, 18/19 = 94.7% "
                  f"accurate); extrapolated = {args.unassign_below}-{args.validated_identity} "
                  f"identity, BELOW the validated range, treat as a hypothesis; "
                  f"below_identity_floor = < {args.unassign_below}, reported as unassigned. "
                  f"Suffix _low_margin: < {args.min_margin} identity margin over the next clade "
                  f"(heuristic; margin did not separate correct from incorrect at n=19). "
                  f"Suffix +provisional_cladeIV: Clade IV has one characterized member.\n")
        out.write("\t".join(cols) + "\n")

        n_written = 0
        for start in range(0, len(hits), args.chunk):
            chunk = hits[start:start + args.chunk]
            cf = os.path.join(workdir, "chunk.fa")
            with open(cf, "w") as fh:
                for n, c, s in backbone:
                    fh.write(f">__BB__{n}\n{s}\n")
                for i, (sample, sid, seq, ev) in enumerate(chunk):
                    fh.write(f">__Q__{start + i}\n{seq}\n")
            sto = os.path.join(workdir, "chunk.sto")
            run_hmmalign(args.hmm, cf, sto)
            seqs, rf = parse_stockholm(sto)
            mcols = match_column_indices(rf)

            def at(aligned: str, pos: int) -> str:
                return aligned[mcols[pmap[pos]]].upper()

            bb_aln = {}
            for n, c, _ in backbone:
                k = next((k for k in seqs if k.split("|")[0] == f"__BB__{n}"), None)
                if k:
                    bb_aln[n] = (c, "".join(seqs[k][i] for i in mcols).upper())

            for i, (sample, sid, seq, ev) in enumerate(chunk):
                k = next((k for k in seqs if k.split("|")[0] == f"__Q__{start + i}"), None)
                if k is None:
                    continue
                aligned = seqs[k]
                qm = "".join(aligned[j] for j in mcols).upper()

                best_by_clade: dict[str, tuple[float, str]] = {}
                for rn, (rc, rm) in bb_aln.items():
                    idv = identity(qm, rm)
                    if rc not in best_by_clade or idv > best_by_clade[rc][0]:
                        best_by_clade[rc] = (idv, rn)
                ranked = sorted(best_by_clade.items(), key=lambda kv: -kv[1][0])
                clade, (cid, cref) = ranked[0]
                margin = cid - ranked[1][1][0] if len(ranked) > 1 else cid
                # Three evidence-derived tiers. The caller was validated by leave-one-out only
                # over 52.7-96.2 % identity; essentially all of BFD falls below that, so a call
                # there is an extrapolation and is labelled as one rather than hidden.
                if cid < args.unassign_below:
                    clade, conf = "unassigned", "below_identity_floor"
                elif cid < args.validated_identity:
                    conf = "extrapolated" if margin >= args.min_margin else "extrapolated_low_margin"
                else:
                    conf = "validated_range" if margin >= args.min_margin else "validated_low_margin"
                if clade == "IV":
                    # Clade IV has exactly one characterized member (MphUPO) and failed
                    # leave-one-out by construction. Never report it as settled.
                    conf += "+provisional_cladeIV"

                res = {p: at(aligned, p) for p in usable}
                pres = "".join(res.get(p, "-") for p in pocket_use)
                pscore = sum(1 for p in pocket_use if res.get(p, "-") in ALIPHATIC)
                cat_detail, cat_ok = [], True
                for p, want in ((39, "C"), (109, "E"), (111, "D"), (180, "E")):
                    if p in res:
                        good = res[p] == want
                        cat_detail.append(f"{p}{res[p]}{'' if good else '!'}")
                        cat_ok &= good
                if 110 in res:
                    b = res[110]
                    cat_detail.append(f"110{b}")
                mapped = sum(1 for p in usable if res.get(p, "-") != "-")

                row = [sample, sid, "NA" if ev != ev else f"{ev:.2g}", str(len(seq)),
                       clade, f"{cid:.1f}", f"{margin:.1f}", conf, cref,
                       str(pscore), str(len(pocket_use)), pres,
                       "yes" if cat_ok else "no", ",".join(cat_detail),
                       f"{mapped}/{len(usable)}", CLADE_T1.get(clade, "n/a")]
                row += [res.get(p, "-") for p in usable]
                out.write("\t".join(row) + "\n")
                n_written += 1
            log(f"  scored {min(start + args.chunk, len(hits))}/{len(hits)}")

        if zstd_proc:
            out.close(); zstd_proc.wait()
        else:
            out.close()
        log(f"wrote {n_written} rows to {args.out}")
    finally:
        if args.keep_temp:
            log(f"temp kept: {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
