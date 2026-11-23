#!/usr/bin/env python3
"""Sensitivity audit of the step HMMs and cutoffs against reviewed fungal UniProt entries of the same EC number.

For every step that has its own EC number in the knowledge-base tables, fetch the reviewed fungal entries of that EC
(UniProt REST: `ec:<EC> AND taxonomy_id:4751 AND reviewed:true`), remove the step's own references and decoys, score the rest
with the step HMM (hmmsearch, E <= 1000, best domain score) and report how many reach the cutoff.

What this measures: recall of the HMM and cutoff against ANNOTATED entries. What it does not measure: truth. UniProt EC
assignments of reviewed entries are often by similarity, and one EC can sit on multi-domain proteins (for example a terpene
cyclase with a prenyltransferase domain), so a miss is a lead to read, not proof of a wrong cutoff. Entries are listed with
their names so each miss can be read.

A step is auditable only when each of its reactions with an EC belongs to that step alone (a reaction shared by two steps, such
as the astaxanthin P450 and its reductase, says nothing about either step's own EC).
"""
from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pathway_census import hits  # noqa: E402

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ec_by_step(kb_dir):
    """step -> sorted list of EC numbers, only from reactions that belong to that step alone."""
    kb = Path(kb_dir)
    ec = {r["reaction_id"]: (r["ec"] or "").strip() for r in read_tsv(kb / "reaction.tsv")}
    steps_of = {}
    for r in read_tsv(kb / "step_reaction.tsv"):
        steps_of.setdefault(r["reaction_id"], set()).add(r["step_id"])
    out = {}
    for rid, steps in steps_of.items():
        if len(steps) == 1 and ec.get(rid):
            out.setdefault(next(iter(steps)), set()).add(ec[rid])
    return {s: sorted(v) for s, v in out.items()}


def read_fasta(text):
    name, hdr, seq = None, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, hdr, "".join(seq)
            hdr = line[1:]
            parts = hdr.split("|")
            name = parts[1] if len(parts) > 2 and parts[0] in ("sp", "tr") else hdr.split()[0]
            seq = []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        yield name, hdr, "".join(seq)


def fetch_fasta(ec, taxon=4751, timeout=120):
    query = f"(ec:{ec}) AND (taxonomy_id:{taxon}) AND (reviewed:true)"
    r = subprocess.run(["curl", "-sf", "-m", str(timeout), "-G", UNIPROT_STREAM, "--data-urlencode", f"query={query}",
                        "--data-urlencode", "format=fasta"], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"UniProt request failed for EC {ec} (curl exit {r.returncode})")
    return r.stdout


def own_accessions(refs_dir, step):
    out = set()
    for f in (Path(refs_dir) / f"{step}.faa", Path(refs_dir) / "decoys" / f"{step}.faa"):
        if f.exists():
            out |= {n for n, _, _ in read_fasta(f.read_text())}
    return out


def best_scores(hmm, fasta_text, cpus=2):
    with tempfile.TemporaryDirectory() as td:
        fa, dt = Path(td) / "q.faa", Path(td) / "d.domtbl"
        fa.write_text(fasta_text)
        hits.run_tool(["hmmsearch", "--cpu", str(cpus), "--noali", "-o", "/dev/null", "--domtblout", str(dt),
                       "-E", "1000", "--domE", "1000", str(hmm), str(fa)])
        df = hits.parse_domtbl(dt)
    if df.empty:
        return {}
    return df.groupby("protein_id").score.max().to_dict()


def audit_step(entries, scores, cutoff, rule_cutoff=None):
    """entries: list of (accession, header). scores: accession -> best score. Return (per-entry rows, summary dict)."""
    rows = []
    for acc, hdr in entries:
        sc = scores.get(acc)
        rows.append(dict(accession=acc, name=hdr.split(" OS=")[0].split(" ", 1)[-1], organism=(hdr.split(" OS=")[1].split(" OX=")[0] if " OS=" in hdr else ""),
                         score=sc, at_cutoff=(sc is not None and sc >= cutoff),
                         at_rule=(None if rule_cutoff is None else (sc is not None and sc >= rule_cutoff))))
    n = len(rows)
    found = [r["score"] for r in rows if r["score"] is not None]
    summ = dict(n_entries=n, n_no_hit=sum(r["score"] is None for r in rows),
                n_at_cutoff=sum(r["at_cutoff"] for r in rows),
                recall_at_cutoff=(sum(r["at_cutoff"] for r in rows) / n if n else None),
                median_score=(statistics.median(found) if found else None))
    if rule_cutoff is not None:
        summ["n_at_rule"] = sum(bool(r["at_rule"]) for r in rows)
        summ["recall_at_rule"] = summ["n_at_rule"] / n if n else None
    return rows, summ


def main(argv=None, fetch=fetch_fasta):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--kb-dir", required=True, type=Path)
    ap.add_argument("--refs-dir", required=True, type=Path)
    ap.add_argument("--hmm-dir", required=True, type=Path)
    ap.add_argument("--min-scores", required=True, type=Path, help="min_scores TSV used in the run (has min_score and rule_min_score)")
    ap.add_argument("--tag", default="PROVISIONAL")
    ap.add_argument("--out-dir", required=True, type=Path)
    a = ap.parse_args(argv)
    ms = {r["step"]: r for r in read_tsv(a.min_scores)}
    ecs = ec_by_step(a.kb_dir)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    summary, per_entry, skipped = [], [], []
    for step in sorted(ms):
        if step not in ecs:
            skipped.append((step, "no EC of its own in the tables (shared reaction or none)"))
            continue
        seen, entries = set(), []
        for ec in ecs[step]:
            for acc, hdr, seq in read_fasta(fetch(ec)):
                if acc not in seen:
                    seen.add(acc)
                    entries.append((acc, hdr, seq))
        excluded = own_accessions(a.refs_dir, step)
        kept = [(acc, hdr, seq) for acc, hdr, seq in entries if acc not in excluded]
        fasta = "".join(f">{acc}\n{seq}\n" for acc, _, seq in kept)
        hmm = a.hmm_dir / f"{step}.{a.tag}.hmm" if a.tag else a.hmm_dir / f"{step}.hmm"
        scores = best_scores(hmm, fasta) if kept else {}
        cutoff = float(ms[step]["min_score"])
        rule = float(ms[step]["rule_min_score"]) if ms[step].get("rule_min_score") else None
        rows, summ = audit_step([(acc, hdr) for acc, hdr, _ in kept], scores, cutoff, rule)
        summ.update(step=step, ec=",".join(ecs[step]), cutoff=cutoff, rule_cutoff=rule, n_fetched=len(entries), n_own_excluded=len(entries) - len(kept))
        summary.append(summ)
        for r in rows:
            per_entry.append(dict(step=step, **r))
    cols = ["step", "ec", "cutoff", "rule_cutoff", "n_fetched", "n_own_excluded", "n_entries", "n_no_hit", "n_at_cutoff", "recall_at_cutoff", "n_at_rule", "recall_at_rule", "median_score"]
    with open(a.out_dir / "sensitivity_summary.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(summary)
    with open(a.out_dir / "sensitivity_entries.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=["step", "accession", "name", "organism", "score", "at_cutoff", "at_rule"], delimiter="\t")
        w.writeheader()
        w.writerows(per_entry)
    with open(a.out_dir / "not_auditable.tsv", "w") as fh:
        fh.write("step\treason\n")
        for s, why in skipped:
            fh.write(f"{s}\t{why}\n")
    print(f"audited {len(summary)} steps, {len(skipped)} not auditable")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
