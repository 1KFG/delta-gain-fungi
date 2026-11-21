#!/usr/bin/env python3.12
"""Assign triage bins to assigned pathway copies. Rules: pathway_census/bins.py docstring.

Input --copies: TSV (plain, .gz or .zst) of already-joined copies; columns are listed in the
bins.py docstring. Boolean columns must be exactly True/False (empty = NA).
Output: TSV with a first line "# params: ..."; read it with pandas comment="#".
Refuses a verified: false pathway unless --provisional; the output name then gets .PROVISIONAL.
Bins are triage, not a ranking and not an activity prediction.
"""
import argparse
import gzip
import io
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pathway_census import hits  # noqa: E402
from pathway_census.bins import DEFAULT_REQUIRE, triage  # noqa: E402
from pathway_census.defs import load_pathway  # noqa: E402
from pathway_census.matrix import parse_bool_col  # noqa: E402

BOOL_COLS = ("truncated", "contam_suspect", "domain_ok", "motif_ok", "hap_dup_suspect")


def read_tsv(path: Path) -> pd.DataFrame:
    kw = dict(sep="\t", comment="#", dtype=str, keep_default_na=False)
    if str(path).endswith(".zst"):
        data = subprocess.run(["zstd", "-dc", str(path)], check=True, capture_output=True).stdout
        return pd.read_csv(io.BytesIO(data), **kw)
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            return pd.read_csv(fh, **kw)
    return pd.read_csv(path, **kw)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pathway", type=Path, required=True)
    ap.add_argument("--copies", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True, help="output TSV (.zst suffix compresses)")
    ap.add_argument("--window-bp", type=int, default=20000, help="PROPOSED; cluster span UNVERIFIED")
    ap.add_argument("--tier1-min", type=int, default=3)
    ap.add_argument("--tier2-min", type=int, default=2)
    ap.add_argument("--divergent-max-identity", type=float, default=0.5)
    ap.add_argument("--expansion-min", type=int, default=2)
    ap.add_argument("--min-length", type=int, default=100)
    ap.add_argument("--max-length", type=int, default=1500)
    ap.add_argument("--tm-max", type=int, default=1)
    ap.add_argument("--soluble-steps", nargs="*", default=[])
    ap.add_argument("--require", nargs="+", default=list(DEFAULT_REQUIRE))
    ap.add_argument("--provisional", action="store_true")
    a = ap.parse_args()
    pw = load_pathway(a.pathway)
    try:
        tag = hits.provisional_tag(pw, a.provisional)
    except hits.UnverifiedPathwayError as e:
        sys.exit(f"ERROR: {e}")
    try:
        c = read_tsv(a.copies)
        for col in BOOL_COLS + ("signal_peptide",):
            if col in c.columns:
                c[col] = parse_bool_col(c[col], col, allow_na=(col != "signal_peptide"))
        for col in ("ref_identity", "copy_filtered", "start", "end", "length", "tm_helices"):
            if col in c.columns:
                c[col] = pd.to_numeric(c[col].replace("", pd.NA), errors="raise")
        out = triage(c, a.window_bp, tier1_min=a.tier1_min, tier2_min=a.tier2_min,
                     divergent_max_identity=a.divergent_max_identity, expansion_min=a.expansion_min,
                     min_length=a.min_length, max_length=a.max_length, tm_max=a.tm_max,
                     soluble_steps=tuple(a.soluble_steps), require=tuple(a.require))
    except ValueError as e:
        sys.exit(f"ERROR: {e}")
    if not a.soluble_steps:
        print("NOTE: --soluble-steps is empty; signal/TM check OFF", file=sys.stderr)
    zst = str(a.out).endswith(".zst")
    plain = Path(str(a.out)[:-4]) if zst else a.out
    plain = plain.with_name(hits.tagged_name(plain.stem, plain.suffix.lstrip(".") or "tsv", tag))
    params = (f"window_bp={a.window_bp} tier1_min={a.tier1_min} tier2_min={a.tier2_min} "
              f"divergent_max_identity={a.divergent_max_identity} expansion_min={a.expansion_min} "
              f"min_length={a.min_length} max_length={a.max_length} tm_max={a.tm_max} "
              f"soluble_steps={','.join(a.soluble_steps) or 'NONE (signal/TM check OFF)'} require={','.join(a.require)} "
              f"pathway={pw.id} (all PROPOSED)")
    plain.parent.mkdir(parents=True, exist_ok=True)
    with open(plain, "w") as fh:
        fh.write(f"# params: {params}\n")
        if tag:
            fh.write(f"# {tag}: pathway '{pw.id}' verified: false; not user-reviewed\n")
        out.to_csv(fh, sep="\t", index=False)
    if zst:
        subprocess.run(["zstd", "-q", "-f", "-T0", "--rm", str(plain), "-o", str(plain) + ".zst"], check=True)
        plain = Path(str(plain) + ".zst")
    print(f"{len(out)} rows -> {plain}; bins: {out.bin.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
