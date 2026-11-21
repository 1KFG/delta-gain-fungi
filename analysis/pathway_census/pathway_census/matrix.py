"""Copy matrix and phylum summary (spec 5.4 step 5, anti-goals).

build_matrix: long table, one row per (genome, step): the FULL cross product of genome table x
pathway steps, zeros included.

Columns
-------
locustag, step_id
copy_raw        number of `assigned` rows for (genome, step). `unassigned` rows never count (a bare
                HMM hit is not a copy). A protein with a fused-step pair has one row per step and is
                counted once in each step.
copy_filtered   `assigned` rows where clean_mask(flags, require) is True, i.e. every required flag
                is known False. Default require = ("truncated", "contam_suspect").
copy_unknown    `assigned` rows with no required flag known True and at least one required flag NA
                (not evaluated). A missing flags row counts as not evaluated (all flags NA).
copy_flagged    `assigned` rows with ANY required flag known True (excluded for a known reason),
                even if another required flag is NA.
                The three buckets are mutually exclusive (priority: flagged, then unknown, then
                filtered), so copy_filtered + copy_unknown + copy_flagged == copy_raw for every row.
n_hap_suspect / n_hap_unknown   count of assigned rows with hap_dup_suspect True / NA. These enter
                copy_filtered only if the caller puts "hap_dup_suspect" in `require`.
n_proteins_distinct   distinct proteins with any assigned row in the genome (repeated on each of the
                genome's rows). A fused protein counts once here.
has_proteins    from the genome table. False means no proteome was searched: counts are 0 but this
                is NOT evidence of absence. phylum_summary excludes such genomes from n_genomes.
                Parsed with parse_bool_col(allow_na=False) in build_matrix and phylum_summary: real
                booleans or the exact strings "True"/"False"; anything else raises ValueError.
genome_flags    ";"-joined tokens: no_proteins; unknown_flags (copy_unknown > 0).

Errors (ValueError): assigned protein whose 8-hex LOCUSTAG prefix is not in the genome table;
step_id not in the pathway; duplicate (protein_id, step_id) among assigned rows; assigned rows in a
genome with has_proteins False.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .quality import clean_mask

DEFAULT_REQUIRE = ("truncated", "contam_suspect")
MATRIX_COLUMNS = ["locustag", "step_id", "copy_raw", "copy_filtered", "copy_unknown", "copy_flagged",
                  "n_hap_suspect", "n_hap_unknown", "n_proteins_distinct", "has_proteins", "genome_flags"]
SUMMARY_COLUMNS = ["phylum", "step_id", "n_genomes", "n_with_step", "median_copy", "max_copy",
                   "frac_with_step", "n_unknown", "copy_raw_total", "copy_filtered_total",
                   "copy_unknown_total", "copy_flagged_total"]


def parse_bool_col(series: pd.Series, name: str, allow_na: bool) -> pd.Series:
    """Strict boolean parse. Accepts bool values and the strings "True" / "False" (exact).

    Empty string, None and NaN become NA if allow_na, else ValueError. Any other value
    (for example "1", "yes", "nan", "true", "FALSE ") raises ValueError naming the column and value.
    Returns the pandas nullable boolean dtype.
    """
    out = []
    for v in series:
        if isinstance(v, (bool, np.bool_)):
            out.append(bool(v))
        elif v is None or v is pd.NA or (isinstance(v, float) and np.isnan(v)) or (isinstance(v, str) and v == ""):
            if not allow_na:
                raise ValueError(f"column '{name}': missing value not allowed")
            out.append(pd.NA)
        elif isinstance(v, str) and v in ("True", "False"):
            out.append(v == "True")
        else:
            raise ValueError(f"column '{name}': invalid boolean value {v!r}")
    return pd.Series(out, index=series.index, dtype="boolean")


def build_matrix(assigned: pd.DataFrame, flags: pd.DataFrame, genomes: pd.DataFrame, pathway,
                 require: tuple[str, ...] = DEFAULT_REQUIRE) -> pd.DataFrame:
    steps = [s.id for s in pathway.steps]
    g = genomes[["locustag", "has_proteins"]].copy()
    if g.locustag.duplicated().any():
        raise ValueError("duplicate locustag in genome table")
    g["has_proteins"] = parse_bool_col(g.has_proteins, "has_proteins", allow_na=False).astype(bool)

    a = assigned[assigned.label == "assigned"][["protein_id", "step_id"]].copy()
    a["locustag"] = a.protein_id.str.split("_").str[0]
    bad = sorted(set(a.locustag) - set(g.locustag))
    if bad:
        raise ValueError(f"protein LOCUSTAG(s) not in genomes table: {bad[:5]}")
    bad = sorted(set(a.step_id) - set(steps))
    if bad:
        raise ValueError(f"step_id not in pathway: {bad[:5]}")
    if a.duplicated(["protein_id", "step_id"]).any():
        raise ValueError("duplicate (protein_id, step_id) among assigned rows")
    noprot = set(g.locustag[~g.has_proteins]) & set(a.locustag)
    if noprot:
        raise ValueError(f"assigned rows in genome(s) with has_proteins False: {sorted(noprot)[:5]}")

    # attach flags; a missing flags row means not evaluated (NA), never clean
    cols = list(dict.fromkeys(list(require) + ["hap_dup_suspect"]))
    f = flags[["protein_id", "step_id"] + cols]
    if f.duplicated(["protein_id", "step_id"]).any():
        raise ValueError("duplicate (protein_id, step_id) in flags")
    a = a.merge(f, on=["protein_id", "step_id"], how="left")
    for c in cols:
        a[c] = a[c].astype("boolean")
    a["is_clean"] = clean_mask(a, require=require)
    any_true = a[list(require)].eq(True).any(axis=1).fillna(False).astype(bool) if require else \
        pd.Series(False, index=a.index)
    any_na = a[list(require)].isna().any(axis=1) if require else pd.Series(False, index=a.index)
    a["is_flagged"] = any_true
    a["is_unknown"] = any_na & ~any_true
    a["hap_true"] = a.hap_dup_suspect.eq(True).fillna(False).astype(bool)
    a["hap_na"] = a.hap_dup_suspect.isna()

    per = a.groupby(["locustag", "step_id"]).agg(
        copy_raw=("protein_id", "size"), copy_filtered=("is_clean", "sum"),
        copy_unknown=("is_unknown", "sum"), copy_flagged=("is_flagged", "sum"), n_hap_suspect=("hap_true", "sum"),
        n_hap_unknown=("hap_na", "sum")).reset_index()
    distinct = a.groupby("locustag").protein_id.nunique().rename("n_proteins_distinct").reset_index()

    grid = g.merge(pd.DataFrame({"step_id": steps}), how="cross")
    m = grid.merge(per, on=["locustag", "step_id"], how="left").merge(distinct, on="locustag", how="left")
    count_cols = ["copy_raw", "copy_filtered", "copy_unknown", "copy_flagged", "n_hap_suspect", "n_hap_unknown",
                  "n_proteins_distinct"]
    m[count_cols] = m[count_cols].fillna(0).astype(int)  # fills genome x step pairs with no assigned row
    m["genome_flags"] = [
        ";".join(t for t, on in (("no_proteins", not hp), ("unknown_flags", cu > 0)) if on)
        for hp, cu in zip(m.has_proteins, m.copy_unknown)]
    order = {s: i for i, s in enumerate(steps)}
    m = m.sort_values(["locustag", "step_id"], key=lambda c: c.map(order) if c.name == "step_id" else c,
                      kind="stable").reset_index(drop=True)
    return m[MATRIX_COLUMNS]


def phylum_summary(matrix: pd.DataFrame, genomes: pd.DataFrame, count_col: str = "copy_filtered") -> pd.DataFrame:
    """Per (phylum, step). median_copy, max_copy and frac_with_step use `count_col`
    (default copy_filtered; pass "copy_raw" for the unfiltered summary).

    n_genomes counts genomes with has_proteins True only. n_with_step = genomes with count > 0.
    n_unknown = genomes with copy_unknown > 0 for the step (copies with unevaluated flags).
    copy_raw_total, copy_filtered_total, copy_unknown_total, copy_flagged_total = sums of the matrix
    buckets over the included genomes (has_proteins True). They do not depend on `count_col`, and
    filtered + unknown + flagged == raw holds for the totals too.
    With no contamination table, contam_suspect is NA for every copy, so every copy is copy_unknown and
    copy_filtered (and frac_with_step with the default count_col) is 0 everywhere.
    A phylum with no included genome gets n_genomes 0, NaN median/frac and totals 0.
    """
    ph = genomes[["locustag", "phylum"]].copy()
    ph["phylum"] = ph.phylum.fillna("(missing)")
    m = matrix.merge(ph, on="locustag", how="left", validate="many_to_one")
    steps = list(dict.fromkeys(matrix.step_id))
    inc = m[parse_bool_col(m.has_proteins, "has_proteins", allow_na=False).astype(bool)]
    rows = []
    for p in sorted(set(ph.phylum)):
        for s in steps:
            x = inc[(inc.phylum == p) & (inc.step_id == s)]
            n = len(x)
            w = int((x[count_col] > 0).sum())
            rows.append(dict(phylum=p, step_id=s, n_genomes=n, n_with_step=w,
                             median_copy=float(x[count_col].median()) if n else float("nan"),
                             max_copy=int(x[count_col].max()) if n else 0,
                             frac_with_step=w / n if n else float("nan"),
                             n_unknown=int((x.copy_unknown > 0).sum()),
                             copy_raw_total=int(x.copy_raw.sum()),
                             copy_filtered_total=int(x.copy_filtered.sum()),
                             copy_unknown_total=int(x.copy_unknown.sum()),
                             copy_flagged_total=int(x.copy_flagged.sum())))
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)
