#!/usr/bin/env python3
"""List what in a pathway definition is not yet verified, as Markdown.

Usage: /usr/bin/python3.12 bin/list_unverified_refs.py pathways/<pathway>.yaml

Prints: the pathway-level verified flag, every ref with verified: false,
refs missing an accession or a PMID, and steps that have no verified ref.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Make the package importable when run from any directory (Python __file__,
# resolved at run time; no shell path tricks).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pathway_census.defs import load_pathway, validate  # noqa: E402


def report(path: Path) -> str:
    p = load_pathway(path)
    out: list[str] = [f"# Unverified items: pathway `{p.id}`", ""]
    out.append(f"- Source: `{path}`")
    out.append(f"- Pathway-level `verified`: **{str(p.verified).lower()}**")
    errs = validate(p)
    out.append(f"- `validate()` errors: {len(errs)}")
    out += [f"  - {e}" for e in errs]
    n_refs = sum(len(s.refs) for s in p.steps)
    n_ver = sum(r.verified for s in p.steps for r in s.refs)
    out.append(f"- Refs: {n_refs} total, {n_ver} verified, {n_refs - n_ver} unverified.")
    out += ["", "## Refs with verified: false", ""]
    rows = [(s.id, r) for s in p.steps for r in s.refs if not r.verified]
    if rows:
        out.append("| step | accession | species | pmid |")
        out.append("|---|---|---|---|")
        for sid, r in rows:
            out.append(f"| {sid} | {r.accession or '-'} | {r.species} | {r.pmid or '-'} |")
    else:
        out.append("None.")
    out += ["", "## Refs missing an accession or a PMID", ""]
    miss = [(s.id, r) for s in p.steps for r in s.refs if not (r.accession and r.pmid)]
    out += [f"- {sid}: accession={r.accession or '-'}, pmid={r.pmid or '-'}"
            for sid, r in miss] or ["None."]
    out += ["", "## Steps with no verified ref", ""]
    bare = [s for s in p.steps if not any(r.verified for r in s.refs)]
    out += [f"- {s.id} (optional={str(s.optional).lower()})" for s in bare] or ["None."]
    return "\n".join(out) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("pathway_yaml", type=Path)
    args = ap.parse_args(argv)
    sys.stdout.write(report(args.pathway_yaml))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
