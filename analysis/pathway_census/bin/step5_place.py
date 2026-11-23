#!/usr/bin/env python3.12
"""Census step 5 (reference tree, placement, step assignment) for ONE step. Step-generic.

Subcommands (the step id is always an argument):
  extract      candidates table -> cut domain FASTA + coordinates TSV (BFD protein files)
  select-reps  mmseqs cluster TSV -> representatives (stratified cap by phylum) + report
  clades       reference tree + labels -> clade_members.tsv + rooted tree + clade report
  placements   upo_place `assign` output -> placements TSV for bin/assign_steps.py + detail table
  summarize    gather MEASURED numbers into SUMMARY.tsv
  make-refset  write the reference-alignment FASTA (refs, decoys, members, outgroups, representatives)
  make-labels  table of labelled leaves (references and named members) with their taxonomic class
  subclade-report  candidates per lineage subclade and nearest labelled leaf of each representative
Named members (role named_member) are tree leaves only: never references, decoys or cutoff input.

The clade rule and the tree code are reused from analysis/enzymes_UPO (imported, not edited):
  upo_place/bin/upo_place.py  (assign, edge_bipartitions, to_newick)
  refset/bin/tree_outline.py  (parse_newick, reroot)
All text input is read as plain, .gz or .zst. Run with /usr/bin/python3.12.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import random
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UPO_BIN = ROOT.parent / "enzymes_UPO" / "upo_place" / "bin" / "upo_place.py"
PAD_DEFAULT = 20  # PROPOSED: residues added on each side of the aligned domain


# ---------- generic text input ----------

def open_text(path):
    """Open a plain, .gz or .zst text file for reading (zst through `zstd -dc`)."""
    path = str(path)
    if path.endswith(".zst"):
        out = subprocess.run(["zstd", "-dc", path], check=True, capture_output=True).stdout
        return io.StringIO(out.decode())
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, newline="")


def read_tsv(path, comment="#"):
    """Read a TSV into a list of dicts; lines starting with `comment` are skipped."""
    with open_text(path) as fh:
        lines = [ln for ln in fh if not ln.startswith(comment)]
    return list(csv.DictReader(lines, delimiter="\t"))


def read_fasta_text(text):
    """FASTA text -> list of (id, sequence); id is the first word of the header."""
    out, name, buf = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def write_fasta(fh, name, seq):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), 60):
        fh.write(seq[i:i + 60] + "\n")


def load_upo():
    spec = importlib.util.spec_from_file_location("upo_place", UPO_BIN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------- extract ----------

def locustag_of(protein_id: str) -> str:
    """BFD protein ids start with the LOCUSTAG before the first underscore."""
    return protein_id.split("_")[0]


def protein_path(seqdb: Path, locustag: str) -> Path:
    h = hashlib.sha1(locustag.encode()).hexdigest()[:2]
    return Path(seqdb) / "proteins" / h / f"{locustag}.proteins.fa.gz"


def best_rows(rows, step_id):
    """Keep rows of one step; one row per protein (highest score). Returns (rows, n_dup)."""
    best, n_dup = {}, 0
    for r in rows:
        if r["step_id"] != step_id:
            continue
        s = float(r["score"])
        if r["protein_id"] in best:
            n_dup += 1
            if s <= float(best[r["protein_id"]]["score"]):
                continue
        best[r["protein_id"]] = r
    return list(best.values()), n_dup


def cut_domain(seq: str, ali_from: int, ali_to: int, pad: int = PAD_DEFAULT):
    """Return (start, end, subsequence) with 1-based inclusive coordinates.
    Pads `pad` residues each side and clamps to the protein ends."""
    if ali_from < 1 or ali_to < ali_from:
        raise ValueError(f"bad alignment coordinates {ali_from}-{ali_to}")
    start = max(1, ali_from - pad)
    end = min(len(seq), ali_to + pad)
    if start > end:
        raise ValueError(f"alignment {ali_from}-{ali_to} is outside a protein of length {len(seq)}")
    return start, end, seq[start - 1:end]


def fetch_tag(seqdb, tag, ids):
    """samtools faidx for the ids of one genome. Returns (dict id->seq, error string or '')."""
    p = protein_path(seqdb, tag)
    if not p.exists():
        return {}, "protein_file_missing"
    if not Path(str(p) + ".fai").exists():
        return {}, "fai_missing"
    r = subprocess.run(["samtools", "faidx", str(p), *ids], capture_output=True, text=True)
    return dict(read_fasta_text(r.stdout)), ""


def extract(cand_paths, step_id, seqdb, out_fasta, out_coords, out_fail, pad=PAD_DEFAULT, workers=16):
    rows = []
    for p in cand_paths:
        rows.extend(read_tsv(p))
    sel, n_dup = best_rows(rows, step_id)
    by_tag = defaultdict(list)
    for r in sel:
        by_tag[locustag_of(r["protein_id"])].append(r)
    with ThreadPoolExecutor(workers) as ex:
        res = dict(zip(by_tag, ex.map(lambda t: fetch_tag(seqdb, t, [r["protein_id"] for r in by_tag[t]]),
                                      by_tag)))
    n_ok, fails = 0, []
    with open(out_fasta, "w") as fa, open(out_coords, "w", newline="") as fc:
        w = csv.writer(fc, delimiter="\t", lineterminator="\n")
        w.writerow(["protein_id", "locustag", "protein_length", "ali_from", "ali_to", "cut_from",
                    "cut_to", "cut_length", "score"])
        for tag, rs in by_tag.items():
            seqs, err = res[tag]
            for r in rs:
                pid = r["protein_id"]
                if err or pid not in seqs:
                    fails.append((pid, err or "id_not_in_protein_file"))
                    continue
                try:
                    s, e, sub = cut_domain(seqs[pid], int(r["ali_from"]), int(r["ali_to"]), pad)
                except ValueError as exc:
                    fails.append((pid, f"bad_coordinates:{exc}"))
                    continue
                write_fasta(fa, pid, sub)
                w.writerow([pid, tag, len(seqs[pid]), r["ali_from"], r["ali_to"], s, e, len(sub), r["score"]])
                n_ok += 1
    with open(out_fail, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["protein_id", "reason"])
        w.writerows(fails)
    return {"rows_in_step": len(sel), "duplicate_rows_dropped": n_dup, "genomes": len(by_tag),
            "domains_extracted": n_ok, "failures": len(fails)}


# ---------- representatives ----------

def parse_clusters(path):
    """mmseqs easy-cluster `_cluster.tsv` (rep<TAB>member) -> {rep: [members]} (rep included)."""
    out = defaultdict(list)
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            rep, mem = line.rstrip("\n").split("\t")[:2]
            out[rep].append(mem)
    return dict(out)


def stratified_quota(counts: dict, n: int) -> dict:
    """Allocate n picks to strata in proportion to counts, at least 1 per stratum, never more
    than the stratum size. Largest-remainder method; ties broken by stratum name. If there are
    more strata than n, every stratum still gets 1 (the total then exceeds n; the caller reports it)."""
    total = sum(counts.values())
    if total <= n:
        return dict(counts)
    quota = {k: min(c, max(1, int(n * c / total))) for k, c in counts.items()}
    def remainder(k):
        return (n * counts[k] / total) - int(n * counts[k] / total)
    while sum(quota.values()) < n:
        cand = [k for k in counts if quota[k] < counts[k]]
        if not cand:
            break
        k = max(cand, key=lambda x: (remainder(x), x))
        quota[k] += 1
    while sum(quota.values()) > n:
        cand = [k for k in counts if quota[k] > 1]
        if not cand:
            break
        k = max(cand, key=lambda x: (quota[x] - n * counts[x] / total, x))
        quota[k] -= 1
    return quota


def stratified_draw(reps_by_phylum: dict, n: int, seed: int = 1):
    """reps_by_phylum: {phylum: [rep ids]}. Returns (chosen ids sorted, quota)."""
    quota = stratified_quota({k: len(v) for k, v in reps_by_phylum.items()}, n)
    rng = random.Random(seed)
    chosen = []
    for k in sorted(reps_by_phylum):
        ids = sorted(reps_by_phylum[k])
        chosen.extend(rng.sample(ids, quota[k]))
    return sorted(chosen), quota


def select_reps(cluster_tsv, genomes_tsv, max_reps, seed, out_ids, out_report):
    clusters = parse_clusters(cluster_tsv)
    phylum = {r["locustag"]: (r["phylum"] or "unknown") for r in read_tsv(genomes_tsv)}
    by_phy = defaultdict(list)
    for rep in clusters:
        by_phy[phylum.get(locustag_of(rep), "unknown")].append(rep)
    chosen, quota = stratified_draw(by_phy, max_reps, seed)
    Path(out_ids).write_text("".join(c + "\n" for c in chosen))
    with open(out_report, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["phylum", "n_clusters", "n_domains_in_clusters", "n_selected"])
        for k in sorted(by_phy):
            n_dom = sum(len(clusters[r]) for r in by_phy[k])
            w.writerow([k, len(by_phy[k]), n_dom, quota[k]])
    return {"clusters": len(clusters), "representatives_total": len(clusters),
            "representatives_selected": len(chosen), "phyla": len(by_phy),
            "capped": len(clusters) > max_reps}


# ---------- clades ----------

def all_sides(bip):
    """Every bipartition side of the tree as frozensets (both sides of every edge)."""
    out = set()
    for a, b in bip.values():
        out.add(a)
        out.add(b)
    return out


def smallest_side(sides, must_have, must_not_have):
    """Smallest side that holds all of `must_have` and none of `must_not_have`; None if absent."""
    must_have, must_not_have = frozenset(must_have), frozenset(must_not_have)
    ok = [s for s in sides if must_have <= s and not (s & must_not_have)]
    # ties in size (possible in an unrooted tree) go to the side whose sorted names sort first
    return min(ok, key=lambda x: (len(x), sorted(x))) if ok else None


def define_clades(bip, refs, decoys, outgroups, members=()):
    """Clade membership from the reference tree (root-independent: bipartitions).

    target  = smallest side with all `refs` and no decoy/outgroup. None -> the refs do not form
              a clade free of decoys and outgroups (then decoy = the decoy leaves only).
    decoy   = every leaf that is neither in `target` nor an outgroup: the decoys and the
              representatives outside the target clade. A placement counts for `decoy` only when
              its edge lies inside a subtree made only of such leaves (assign rule: proper subset).
    outgroup= the outgroup leaves (a one-leaf outgroup is a leaf, so nothing can place inside it).
    decoys_share_side_without_refs: some side holds all decoys and no ref (a sanity report, not used
              for the assignment).
    members: named tree members; they are kept out of `decoy` (see below).
    Each leaf gets one label (see write_clade_members)."""
    sides = all_sides(bip)
    everything = frozenset().union(*sides) if sides else frozenset()
    refs, decoys, outgroups = frozenset(refs), frozenset(decoys), frozenset(outgroups)
    members = frozenset(members)
    target = smallest_side(sides, refs, decoys | outgroups)
    if target is None:
        decoy = decoys - outgroups
    else:
        # named members are never decoys: a member outside `target` gets no clade label
        decoy = everything - target - outgroups - members
    shared = bool(decoys) and smallest_side(sides, decoys, refs) is not None
    return {"target": target, "decoy": frozenset(decoy), "outgroup": outgroups,
            "decoys_share_side_without_refs": shared}


def leaf_roles(all_leaves, refs, decoys, outgroups, reps, members=(), outgroup_role="outgroup"):
    """Role per leaf. Priority: outgroup, decoy, reference, named_member, representative, other."""
    role = {}
    for l in all_leaves:
        role[l] = (outgroup_role if l in outgroups else "decoy" if l in decoys else
                   "reference" if l in refs else "named_member" if l in members else
                   "representative" if l in reps else "other")
    return role


def write_clade_members(path, clades, roles, subclades=None):
    """clade_members.tsv: label, our_clade, fig1_clade (read by upo_place `assign`), plus `role`.
    our_clade = target / decoy / outgroup. fig1_clade = the lineage subclade name (second clade
    system of the UPO tool) or empty; without subclades it stays empty."""
    lab = {}
    for name in ("target", "decoy", "outgroup"):
        for l in clades[name] or ():
            lab[l] = name
    sub = {}
    for name, members in (subclades or {}).items():
        for l in members:
            sub[l] = name
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["label", "our_clade", "fig1_clade", "role"])
        for l in sorted(roles):
            w.writerow([l, lab.get(l, ""), sub.get(l, ""), roles[l]])


def ids_in_fasta(path):
    with open(path) as fh:
        return [n for n, _ in read_fasta_text(fh.read())]


def mrca_side(sides, leaves):
    """Smallest bipartition side that holds all `leaves` (the MRCA clade once the tree is rooted
    outside it). None if `leaves` is empty."""
    return smallest_side(sides, leaves, ()) if leaves else None


def describe_outgroups(sides, outgroups, decoys, roles):
    """Where the outgroups and the decoys sit relative to each other.
    outgroups_monophyletic: some side is exactly the outgroup set.
    The MRCA side of the outgroups is reported with its role counts; so is the MRCA side of the decoys."""
    outgroups, decoys = frozenset(outgroups), frozenset(decoys)
    out = {}
    if outgroups:
        o = mrca_side(sides, outgroups)
        out["outgroups_monophyletic"] = outgroups in sides
        out["outgroup_mrca_size"] = len(o)
        out["outgroup_mrca_roles"] = dict(Counter(roles[l] for l in o))
        out["decoys_inside_outgroup_mrca"] = len(decoys & o)
    if decoys:
        d = mrca_side(sides, decoys)
        out["decoys_monophyletic"] = decoys in sides
        out["decoy_mrca_size"] = len(d)
        out["decoy_mrca_roles"] = dict(Counter(roles[l] for l in d))
    return out


def lineage_subclades(sides, labels, excluded=frozenset()):
    """PROPOSED rule. labels: {leaf: class} for the labelled leaves (references and named members).
    For each class with at least 2 labelled leaves: take the smallest side that holds them all. It is
    a named subclade `<class>` only if that side holds no labelled leaf of another class and no leaf
    of `excluded` (decoys and outgroups). Otherwise no subclade is made and the reason is reported.
    Returns (subclades {class: frozenset}, report rows)."""
    by_class = defaultdict(set)
    for leaf, cls in labels.items():
        by_class[cls].add(leaf)
    subclades, rows = {}, []
    for cls in sorted(by_class):
        leaves = frozenset(by_class[cls])
        row = {"class": cls, "n_labelled": len(leaves), "subclade": False, "size": "", "reason": ""}
        if len(leaves) < 2:
            row["reason"] = "fewer than 2 labelled leaves"
        else:
            side = mrca_side(sides, leaves)
            row["size"] = len(side)
            other = sorted(l for l in side if l in labels and labels[l] != cls)
            bad = sorted(side & frozenset(excluded))
            if other:
                row["reason"] = "holds labelled leaves of other classes: " + ",".join(
                    f"{l}({labels[l]})" for l in other)
            elif bad:
                row["reason"] = "holds decoy/outgroup leaves: " + ",".join(bad)
            else:
                row["subclade"] = True
                subclades[cls] = side
        rows.append(row)
    return subclades, rows


def read_labels(path):
    """labelled_leaves.tsv (leaf, role, class, ...) -> {leaf: class}; empty classes are skipped."""
    return {r["leaf"]: r["class"] for r in read_tsv(path) if r.get("class")}


def make_labels(refs_fa, ref_classes_tsv, members_fa, provenance_tsv, step_id, out):
    """Write the table of labelled leaves. Class of a reference comes from `ref_classes_tsv`
    (accession, class); class of a named member from the curated provenance table.
    A leaf without a class raises ValueError (nothing is guessed)."""
    ref_class = {r["accession"]: r["class"] for r in read_tsv(ref_classes_tsv)}
    prov = {r["accession"]: r for r in read_tsv(provenance_tsv)
            if r["step_id"] == step_id and r["role"] == "named_member"}
    rows = []
    for acc in ids_in_fasta(refs_fa):
        if acc not in ref_class:
            raise ValueError(f"reference {acc} has no class in {ref_classes_tsv}")
        rows.append({"leaf": acc, "role": "reference", "class": ref_class[acc], "evidence": ""})
    for acc in (ids_in_fasta(members_fa) if members_fa else []):
        if acc not in prov:
            raise ValueError(f"named member {acc} has no named_member row for step {step_id} in {provenance_tsv}")
        rows.append({"leaf": acc, "role": "named_member", "class": prov[acc]["class"],
                     "evidence": prov[acc]["evidence_codes"]})
    _write_dicts(out, rows, ["leaf", "role", "class", "evidence"])
    return len(rows)


def make_refset(out, refs, decoys, members, outgroups, reps_fa):
    """Write the reference-alignment FASTA: references, decoys, named members, outgroups, then the
    representatives. Headers are cut to the first word. A duplicate id raises ValueError. Nothing is
    dropped: the representative cap applies to candidates only."""
    seen, n = set(), 0
    with open(out, "w") as fo:
        for path in (refs, decoys, members, outgroups, reps_fa):
            if not path:
                continue
            with open(path) as fh:
                for name, seq in read_fasta_text(fh.read()):
                    if name in seen:
                        raise ValueError(f"duplicate id in the reference set: {name}")
                    seen.add(name)
                    write_fasta(fo, name, seq)
                    n += 1
    return n


def run_clades(tree, refs_fa, decoys_fa, outgroup_fa, reps_ids, out_dir, step_id, root_on=None,
               members_fa=None, labels_tsv=None, outgroup_role="outgroup"):
    upo = load_upo()
    out_dir = Path(out_dir)
    refs = set(ids_in_fasta(refs_fa))
    decoys = set(ids_in_fasta(decoys_fa))
    outgroups = set(ids_in_fasta(outgroup_fa)) if outgroup_fa else set()
    members = set(ids_in_fasta(members_fa)) if members_fa else set()
    if members & (refs | decoys | outgroups):
        raise SystemExit("ERROR: a named member is also a reference, decoy or outgroup")
    reps = {ln.strip() for ln in open(reps_ids) if ln.strip()}
    top = upo.parse_newick(open(tree).read())
    leaves = {l.name for l in top.leaves()}
    missing = (refs | decoys | outgroups | members) - leaves
    if missing:
        raise SystemExit(f"ERROR: leaves missing from the tree: {sorted(missing)}")
    # the first edge numbering is irrelevant here: bipartitions are computed from the topology
    for i, n in enumerate(_nodes(top)[1:]):   # the root has no edge
        n.edge = i
    bip = upo.edge_bipartitions(top)
    sides = all_sides(bip)
    rooting = outgroups if outgroups else decoys
    if root_on:
        rooting = set(root_on)
    clades = define_clades(bip, refs, decoys, outgroups, members)
    roles = leaf_roles(leaves, refs, decoys, outgroups, reps, members, outgroup_role)
    rooted = _reroot(upo, top, rooting)
    (out_dir / "ref.rooted.nwk").write_text(rooted + "\n")
    rep = {"step_id": step_id, "n_leaves": len(leaves), "rooting_leaves": sorted(rooting),
           "rooting_source": "outgroup file" if outgroups else "decoys (closest paralogs)",
           "refs_form_clade_free_of_decoys": clades["target"] is not None,
           "decoys_share_side_without_refs": clades["decoys_share_side_without_refs"]}
    if outgroups or members:
        rep.update(describe_outgroups(sides, outgroups, decoys, roles))
    subclades = {}
    if clades["target"] is not None:
        t = clades["target"]
        rep["target_size"] = len(t)
        rep["target_refs"] = len(refs & t)
        rep["target_representatives"] = len(reps & t)
        if members:
            rep["target_named_members"] = len(members & t)
            rep["named_members_outside_target"] = sorted(members - t)
        rep["decoy_size"] = len(clades["decoy"])
        rep["decoy_representatives"] = len(reps & clades["decoy"])
        if labels_tsv:
            labels = read_labels(labels_tsv)
            subclades, rows = lineage_subclades(sides, labels, decoys | outgroups)
            _write_dicts(out_dir / "lineage_subclades.tsv",
                         [dict(r, n_representatives=len(reps & subclades[r["class"]]) if r["subclade"] else "")
                          for r in rows],
                         ["class", "n_labelled", "subclade", "size", "n_representatives", "reason"])
            rep["subclades"] = {c: {"size": len(m), "labelled": len([l for l in m if l in labels]),
                                    "representatives": len(reps & m)} for c, m in subclades.items()}
            rep["classes_without_subclade"] = [r["class"] for r in rows if not r["subclade"]]
        write_clade_members(out_dir / "clade_members.tsv", clades, roles, subclades)
    (out_dir / "clade_report.json").write_text(json.dumps(rep, indent=1) + "\n")
    return rep


def _nodes(root):
    out = [root]
    for c in root.children:
        out.extend(_nodes(c))
    return out


def _reroot(upo, top, rooting_leaves):
    spec = importlib.util.spec_from_file_location("tree_outline", UPO_BIN.parents[2] / "refset" / "bin" / "tree_outline.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("tree_outline", mod)
    spec.loader.exec_module(mod)
    return upo.to_newick(mod.reroot(top, set(rooting_leaves)))


# ---------- placements ----------

def to_placements(assign_rows, step_id, target="target", system="our"):
    """upo_place `assign` rows -> (placements rows, detail rows).

    system "our" (default): placements: protein_id, step_id, clade, support. `clade` is the target
    name only when the assign rule called the target clade; a decoy, outgroup or unassigned call gives
    an EMPTY clade, because bin/assign_steps.py treats any non-empty clade as in-step. support = LWR
    sum of the best clade (kept for every row).
    system "fig1": the clade is the lineage subclade call (its name), empty when there is none;
    support = the subclade LWR sum. `target` is ignored. The target stays the parent: the detail
    table keeps both calls.
    detail keeps both calls, the supports, the stem LWR, the pendant length and `long_pendant`."""
    plc, det = [], []
    for r in assign_rows:
        call = r["our_clade"]
        sub = r.get("fig1_clade", "") or ""
        if system == "our":
            clade, support = (target if call == target else ""), r["our_lwr"]
        elif system == "fig1":
            clade, support = ("" if sub in ("", "unassigned") else sub), r.get("fig1_lwr", "0.000")
        else:
            raise ValueError(f"unknown clade system: {system}")
        plc.append({"protein_id": r["query"], "step_id": step_id, "clade": clade, "support": support})
        det.append({"protein_id": r["query"], "step_id": step_id, "clade_call": call,
                    "best_clade": r["our_best"], "support": r["our_lwr"],
                    "stem_lwr": r["our_stem_lwr"], "pendant_length": r["pendant_length"],
                    "best_edge_lwr": r["best_edge_lwr"],
                    "long_pendant": "long_pendant" in (r.get("flags") or ""),
                    "subclade_call": sub, "subclade_support": r.get("fig1_lwr", "")})
    return plc, det


def run_placements(assign_tsv, step_id, out_plc, out_detail, target="target", system="our"):
    rows = read_tsv(assign_tsv)
    plc, det = to_placements(rows, step_id, target, system)
    _write_dicts(out_plc, plc, ["protein_id", "step_id", "clade", "support"])
    _write_dicts(out_detail, det, list(det[0]) if det else ["protein_id"])
    return {"queries": len(rows)}


# ---------- nearest labelled leaf, subclade counts ----------

def leaf_depths(root):
    """{leaf: list of (node id, depth)} from the leaf up to the root; depth = distance from the root."""
    depth = {}

    def walk(n, d):
        depth[id(n)] = d
        for c in n.children:
            walk(c, d + c.length)
    walk(root, 0.0)
    out = {}
    for l in root.leaves():
        path, n = [], l
        while n is not None:
            path.append((id(n), depth[id(n)]))
            n = n.parent
        out[l.name] = path
    return out


def tree_distance(paths, a, b):
    """Patristic distance between leaves a and b: depth(a)+depth(b)-2*depth(LCA)."""
    pa, pb = paths[a], paths[b]
    ids_b = {i for i, _ in pb}
    lca_depth = next(d for i, d in pa if i in ids_b)
    return pa[0][1] + pb[0][1] - 2 * lca_depth


def nearest_labelled(root, queries, labelled):
    """For each query leaf, the nearest labelled leaf by tree distance.
    Returns [(query, nearest, distance)]; ties go to the leaf name that sorts first."""
    paths = leaf_depths(root)
    labelled = sorted(labelled)
    out = []
    for q in sorted(queries):
        best = min(((tree_distance(paths, q, l), l) for l in labelled), key=lambda x: (x[0], x[1]))
        out.append((q, best[1], best[0]))
    return out


def subclade_counts(detail_rows, subclade_info, target_name="target"):
    """Counts of placed candidates per subclade. subclade_info: {name: {...}} from clade_report.
    A candidate counts for a subclade when its subclade_call is that name (LWR sum >= the
    --min-lwr used in assign). Also: target calls in no subclade, subclade calls outside target."""
    c = Counter(r["subclade_call"] for r in detail_rows if r["subclade_call"] not in ("", "unassigned"))
    out = {k: c.get(k, 0) for k in subclade_info}
    tgt = [r for r in detail_rows if r["clade_call"] == target_name]
    out["_target_total"] = len(tgt)
    out["_target_in_no_subclade"] = sum(r["subclade_call"] in ("", "unassigned") for r in tgt)
    out["_subclade_call_not_target"] = sum(
        r["subclade_call"] not in ("", "unassigned") and r["clade_call"] != target_name for r in detail_rows)
    return out


def run_subclade_report(tree, clade_members, labels_tsv, detail_tsv, clusters, out_dir):
    upo = load_upo()
    out_dir = Path(out_dir)
    top = upo.parse_newick(open(tree).read())
    labels = read_labels(labels_tsv)
    members = read_tsv(clade_members)
    reps = [r["label"] for r in members if r["role"] == "representative"]
    size = {rep: len(m) for rep, m in parse_clusters(clusters).items()} if clusters else {}
    rows = [{"representative": q, "nearest_labelled": n, "class": labels[n],
             "distance": f"{d:.4f}", "n_candidates_in_cluster": size.get(q, "")}
            for q, n, d in nearest_labelled(top, reps, labels)]
    _write_dicts(out_dir / "nearest_labelled.tsv", rows,
                 ["representative", "nearest_labelled", "class", "distance", "n_candidates_in_cluster"])
    by = defaultdict(lambda: [0, 0, []])
    for r in rows:
        e = by[r["class"]]
        e[0] += 1
        e[1] += int(r["n_candidates_in_cluster"] or 0)
        e[2].append(float(r["distance"]))
    summ = [{"class": k, "n_representatives": v[0], "n_candidates_in_their_clusters": v[1],
             "median_distance": f"{sorted(v[2])[len(v[2]) // 2]:.4f}"} for k, v in sorted(by.items())]
    _write_dicts(out_dir / "nearest_labelled_by_class.tsv", summ,
                 ["class", "n_representatives", "n_candidates_in_their_clusters", "median_distance"])
    rep = json.loads((out_dir / "clade_report.json").read_text())
    counts = subclade_counts(read_tsv(detail_tsv), rep.get("subclades", {}))
    srows = [{"subclade": k, "size": v["size"], "n_labelled": v["labelled"],
              "n_representatives": v["representatives"], "n_candidates_called": counts[k]}
             for k, v in rep.get("subclades", {}).items()]
    _write_dicts(out_dir / "subclade_counts.tsv", srows,
                 ["subclade", "size", "n_labelled", "n_representatives", "n_candidates_called"])
    (out_dir / "subclade_summary.json").write_text(json.dumps(
        {k: v for k, v in counts.items() if k.startswith("_")}, indent=1) + "\n")
    return counts


def _write_dicts(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


# ---------- summary ----------

def support_distribution(values, edges=(0.5, 0.7, 0.9, 0.95)):
    """Counts of values >= each edge, plus min/median/max. Empty input gives zeros."""
    v = sorted(float(x) for x in values)
    out = {f"n_ge_{e}": sum(x >= e for x in v) for e in edges}
    out["n"] = len(v)
    out["median"] = v[len(v) // 2] if v else 0.0
    out["min"] = v[0] if v else 0.0
    out["max"] = v[-1] if v else 0.0
    return out


def genomes_per_assigned(assigned_rows):
    """Distinct genomes (LOCUSTAG prefix) holding assigned proteins, and the copy distribution."""
    per = Counter(locustag_of(r["protein_id"]) for r in assigned_rows if r["label"] == "assigned")
    return {"assigned_proteins": sum(per.values()), "genomes_with_assigned": len(per),
            "copy_distribution": dict(sorted(Counter(per.values()).items()))}


def run_summarize(out_dir, assigned_tsv, summary_path):
    out_dir = Path(out_dir)
    items = []  # (key, value, tag)

    def add(k, v, tag="MEASURED"):
        items.append((k, v, tag))
    for fn, label in (("extract_summary.json", "extract"), ("reps_summary.json", "reps"),
                      ("clade_report.json", "clades"), ("subclade_summary.json", "subclades")):
        p = out_dir / fn
        if p.exists():
            for k, v in json.loads(p.read_text()).items():
                add(f"{label}.{k}", json.dumps(v) if isinstance(v, (list, dict)) else v)
    sc = out_dir / "subclade_counts.tsv"
    if sc.exists():
        for r in read_tsv(sc):
            add(f"subclade_candidates.{r['subclade']}", r["n_candidates_called"])
    det = out_dir / "placements_detail.tsv"
    if det.exists():
        d = read_tsv(det)
        add("placement.queries", len(d))
        for k, n in sorted(Counter(r["clade_call"] for r in d).items()):
            add(f"placement.call.{k}", n)
        add("placement.long_pendant", sum(r["long_pendant"] == "True" for r in d))
        tgt = [r["support"] for r in d if r["clade_call"] == "target"]
        for k, v in support_distribution(tgt).items():
            add(f"support_target.{k}", v)
        for k, v in support_distribution([r["support"] for r in d], (0.5, 0.9)).items():
            add(f"support_all_queries_best_clade.{k}", v)
    if assigned_tsv and Path(assigned_tsv).exists():
        a = read_tsv(assigned_tsv)
        for k, n in sorted(Counter(r["label"] for r in a).items()):
            add(f"assign_steps.label.{k}", n)
        add("assign_steps.fusion_candidate_true", sum(r["fusion_candidate"] == "True" for r in a))
        g = genomes_per_assigned(a)
        for k, v in g.items():
            add(f"assign_steps.{k}", json.dumps(v) if isinstance(v, dict) else v)
    with open(summary_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["key", "value", "status"])
        w.writerows(items)
    return len(items)


# ---------- CLI ----------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract")
    e.add_argument("--candidates", nargs="+", required=True)
    e.add_argument("--step", required=True)
    e.add_argument("--seqdb", type=Path, required=True)
    e.add_argument("--pad", type=int, default=PAD_DEFAULT)
    e.add_argument("--workers", type=int, default=16)
    e.add_argument("--out-fasta", required=True)
    e.add_argument("--out-coords", required=True)
    e.add_argument("--out-fail", required=True)
    e.add_argument("--out-json", required=True)
    s = sub.add_parser("select-reps")
    s.add_argument("--clusters", required=True)
    s.add_argument("--genomes", required=True)
    s.add_argument("--max-reps", type=int, default=300)
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--out-ids", required=True)
    s.add_argument("--out-report", required=True)
    s.add_argument("--out-json", required=True)
    c = sub.add_parser("clades")
    c.add_argument("--tree", required=True)
    c.add_argument("--refs", required=True)
    c.add_argument("--decoys", required=True)
    c.add_argument("--outgroup", "--outgroup-fa", dest="outgroup",
                   help="FASTA of one or more outgroup leaves (rooted on their MRCA); omit to root on the decoys")
    c.add_argument("--members-fa", help="FASTA of named tree members (role named_member)")
    c.add_argument("--labels", help="labelled_leaves.tsv from make-labels; turns on lineage subclades")
    c.add_argument("--outgroup-role", default="outgroup", help="role written for outgroup leaves")
    c.add_argument("--reps", required=True)
    c.add_argument("--step", required=True)
    c.add_argument("--out-dir", required=True)
    p = sub.add_parser("placements")
    p.add_argument("--assign", required=True)
    p.add_argument("--step", required=True)
    p.add_argument("--target", default="target")
    p.add_argument("--system", choices=["our", "fig1"], default="our",
                   help="our = target call (default); fig1 = lineage subclade call")
    p.add_argument("--out", required=True)
    p.add_argument("--out-detail", required=True)
    r_ = sub.add_parser("make-refset")
    r_.add_argument("--refs", required=True)
    r_.add_argument("--decoys", required=True)
    r_.add_argument("--members")
    r_.add_argument("--outgroups")
    r_.add_argument("--reps", required=True)
    r_.add_argument("--out", required=True)
    lb = sub.add_parser("make-labels")
    lb.add_argument("--refs", required=True)
    lb.add_argument("--ref-classes", required=True, help="TSV: accession, class")
    lb.add_argument("--members")
    lb.add_argument("--provenance", required=True)
    lb.add_argument("--step", required=True)
    lb.add_argument("--out", required=True)
    sr = sub.add_parser("subclade-report")
    sr.add_argument("--tree", required=True)
    sr.add_argument("--clade-members", required=True)
    sr.add_argument("--labels", required=True)
    sr.add_argument("--detail", required=True)
    sr.add_argument("--clusters")
    sr.add_argument("--out-dir", required=True)
    m = sub.add_parser("summarize")
    m.add_argument("--out-dir", required=True)
    m.add_argument("--assigned")
    m.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "extract":
        r = extract(a.candidates, a.step, a.seqdb, a.out_fasta, a.out_coords, a.out_fail,
                    a.pad, a.workers)
        Path(a.out_json).write_text(json.dumps(r, indent=1) + "\n")
    elif a.cmd == "select-reps":
        r = select_reps(a.clusters, a.genomes, a.max_reps, a.seed, a.out_ids, a.out_report)
        Path(a.out_json).write_text(json.dumps(r, indent=1) + "\n")
    elif a.cmd == "clades":
        r = run_clades(a.tree, a.refs, a.decoys, a.outgroup, a.reps, a.out_dir, a.step,
                       members_fa=a.members_fa, labels_tsv=a.labels, outgroup_role=a.outgroup_role)
        if not r["refs_form_clade_free_of_decoys"]:
            print("RESULT: the reference sequences do not form a clade free of decoys/outgroup",
                  file=sys.stderr)
            print(json.dumps(r), file=sys.stderr)
            sys.exit(3)
    elif a.cmd == "placements":
        r = run_placements(a.assign, a.step, a.out, a.out_detail, a.target, a.system)
    elif a.cmd == "make-refset":
        r = make_refset(a.out, a.refs, a.decoys, a.members, a.outgroups, a.reps)
    elif a.cmd == "make-labels":
        r = make_labels(a.refs, a.ref_classes, a.members, a.provenance, a.step, a.out)
    elif a.cmd == "subclade-report":
        r = run_subclade_report(a.tree, a.clade_members, a.labels, a.detail, a.clusters, a.out_dir)
    else:
        r = run_summarize(a.out_dir, a.assigned, a.out)
    print(json.dumps(r), file=sys.stderr)


if __name__ == "__main__":
    main()
