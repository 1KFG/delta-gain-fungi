"""Paired deposited-versus-BFD comparison (spec 5.6, G6).

Design. Genomes that have BOTH a deposited annotation and a BFD annotation of the same assembly are
compared. The same pipeline (same HMMs, same placement) runs once per protein source and gives two
copy matrices (`matrix.build_matrix`). This module compares the two.

What the test shows: the effect of the annotation source on identical assemblies.
What it does NOT show: anything about the genomes with no prior annotation (BFD-scoped, MEASURED,
Task 2 genome table: 18,055 no_prior_annotation_present and 673 no_prior_annotation_not_downloaded;
spec 5.6 gives the TSV-wide 18,145). Those genomes are outside the pairing.

Keys. A genome is identified by the caller-supplied `locustag` in BOTH matrices. Deposited protein IDs
do not start with the LOCUSTAG prefix before the first underscore (8 hex characters in all but 3 BFD genomes). The caller builds the deposited matrix with
the genome's locustag, never from a protein-ID prefix. This module takes counts only. It never joins
or deduplicates by protein ID, so a protein ID that occurs in both sources stays counted once per
source.

Pair rule (controller ruling). A genome is a pair only if prior_status == "baseline_ok", the BFD
proteins exist (genome table has_proteins True, and the BFD matrix has the genome with
has_proteins True), the deposited proteins exist (genome table has_deposited_proteins True) and the
deposited matrix has the genome with has_proteins True (the deposited file was found and searched).
A genome flag (has_proteins, has_deposited_proteins) that is NA is treated as not True. NA and False
share the same reason code. The extra key `unknown_flag` in `n_unpaired` counts how many of those
non-pairs were caused by an NA flag. It is a subset of the reason counts, not an additional reason, so
do not add it to their sum.
Each non-pair is counted once in `n_unpaired` under the FIRST failing reason, in this order:
not_baseline, no_bfd_proteins, no_deposited_proteins, deposited_file_missing.

paired_table columns: locustag, step_id, copy_deposited, copy_bfd, delta (= bfd - deposited, copy_raw),
annotation_provider, busco_complete, copy_filtered_deposited, copy_filtered_bfd, unknown_deposited,
unknown_bfd, delta_filtered. delta_filtered = filtered_bfd - filtered_deposited where BOTH unknown
counts are 0, else NA (a quality-flag unknown makes the filtered count incomplete; no fillna).

summarize_paired: bootstrap over GENOMES. For each (step, provider) group the per-genome deltas are
resampled with replacement n_boot times (numpy default_rng(seed), integer indices, one draw of
shape (n_boot, n_genomes)); the CI is the 2.5 and 97.5 percentiles of the resampled means.
PROPOSED defaults (no measured basis): min_pairs=20, n_boot=2000, seed=20261120.
Resampling is per step: each step is resampled on its own; steps are NOT resampled jointly, so the
intervals of different steps are not independent draws of the same genomes.
Reserved provider names: "ALL" (overall row) and "(missing)" (null provider). A real provider with
either name raises ValueError.
`dropped_providers` uses the column name `provider`; the summary uses `annotation_provider`.
A provider with fewer than min_pairs pairs gets no per-provider row. It is returned in the
`dropped_providers` frame. It stays in the "ALL" rows. A null provider is labelled "(missing)".
The CI of a group with all deltas equal collapses to that value (for example 0, 0); that is correct
for a resample of identical values and is not a sign of an error.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

UNPAIRED_REASONS = ("not_baseline", "no_bfd_proteins", "no_deposited_proteins", "deposited_file_missing")
UNKNOWN_FLAG_KEY = "unknown_flag"  # subset count, NOT an extra reason: see module docstring
RESERVED_PROVIDERS = ("ALL", "(missing)")
TABLE_COLUMNS = ["locustag", "step_id", "copy_deposited", "copy_bfd", "delta", "annotation_provider",
                 "busco_complete", "copy_filtered_deposited", "copy_filtered_bfd", "unknown_deposited",
                 "unknown_bfd", "delta_filtered"]
SUMMARY_COLUMNS = ["step_id", "annotation_provider", "n_pairs", "mean_delta", "ci_low", "ci_high",
                   "frac_bfd_more", "frac_bfd_less"]
DEFAULT_MIN_PAIRS = 20   # PROPOSED
DEFAULT_N_BOOT = 2000    # PROPOSED
DEFAULT_SEED = 20261120  # PROPOSED
_ACC = re.compile(r"^(GC[AF]_[0-9]+\.[0-9]+)")


def accession_from_asmid(asmid: str) -> str:
    m = _ACC.match(str(asmid))
    if not m:
        raise ValueError(f"cannot read accession from ASMID {asmid!r}")
    return m.group(1)


def deposited_protein_path(deposited_root, folder: str) -> Path:
    """Deposited protein file of one genome, from the prior-annotation TSV `folder` column.

    Order: (1) <deposited_root>/<folder>/<folder>_protein.faa.gz (the MEASURED NCBI_ASM layout);
    (2) <deposited_root>/<folder>/<folder>_protein.faa, used only if (1) does not exist and (2) does.
    If neither exists, path (1) is returned (the caller's existence check then fails).
    A null or empty folder raises ValueError.
    """
    if folder is None or pd.isna(folder) or str(folder) == "":
        raise ValueError("deposited_protein_path needs a non-empty folder (prior-annotation TSV `folder`)")
    d = Path(deposited_root) / str(folder)
    gz, plain = d / f"{folder}_protein.faa.gz", d / f"{folder}_protein.faa"
    if not gz.exists() and plain.exists():
        return plain
    return gz


def _scalar_true(v) -> bool:
    """True only for a known True. NA, None and NaN are not True."""
    return (not pd.isna(v)) and bool(v)


def _is_true(series: pd.Series) -> pd.Series:
    return series.astype("boolean").eq(True).fillna(False).astype(bool)


def count_pairs(genomes: pd.DataFrame, deposited_root) -> dict:
    """Count pairs from the genome table and file existence. Read-only.

    `genomes` needs the load_genome_table columns, including `folder` (prior-annotation TSV). The
    deposited file is found with `deposited_protein_path(deposited_root, folder)`. A null folder
    counts as deposited_file_missing. Returns n_genomes, n_pairs, n_unpaired (reason -> count) and
    pair_locustags.
    """
    if "folder" not in genomes.columns:
        raise ValueError("genome table needs a 'folder' column to locate "
                         "<deposited_root>/<folder>/<folder>_protein.faa.gz")
    base = genomes.prior_status.eq("baseline_ok")
    bfd = _is_true(genomes.has_proteins)
    dep = _is_true(genomes.has_deposited_proteins)

    def _found(fo) -> bool:
        if fo is None or pd.isna(fo) or str(fo) == "":
            return False
        return os.path.exists(deposited_protein_path(deposited_root, fo))

    exists = pd.Series([bool(b and bf and d and _found(fo))
                        for b, bf, d, fo in zip(base, bfd, dep, genomes.folder)], index=genomes.index)
    reason = np.select([~base, ~bfd, ~dep, ~exists], list(UNPAIRED_REASONS[:4]), default="")
    n_unpaired = {r: int((reason == r).sum()) for r in UNPAIRED_REASONS}
    na_bfd = genomes.has_proteins.isna() & (reason == "no_bfd_proteins")
    na_dep = genomes.has_deposited_proteins.isna() & (reason == "no_deposited_proteins")
    n_unpaired[UNKNOWN_FLAG_KEY] = int((na_bfd | na_dep).sum())
    return dict(n_genomes=len(genomes), n_pairs=int((reason == "").sum()), n_unpaired=n_unpaired,
                pair_locustags=list(genomes.locustag[reason == ""]))


def _per_genome_ok(matrix: pd.DataFrame) -> dict:
    return matrix.groupby("locustag").has_proteins.apply(lambda s: bool(_is_true(s).all())).to_dict()


def paired_table(matrix_deposited: pd.DataFrame, matrix_bfd: pd.DataFrame, genomes: pd.DataFrame,
                 return_unpaired: bool = False):
    """Join the two matrices per (locustag, step_id) for paired genomes only. See module docstring.

    `n_unpaired` is stored in the result's `.attrs`, and returned as a second value if
    return_unpaired is True.
    """
    known = set(genomes.locustag)
    for name, m in (("deposited", matrix_deposited), ("bfd", matrix_bfd)):
        bad = sorted(set(m.locustag) - known)
        if bad:
            raise ValueError(f"{name} matrix locustag(s) not in genome table: {bad[:5]}")
    if genomes.locustag.duplicated().any():
        raise ValueError("duplicate locustag in genome table")
    dep_ok, bfd_ok = _per_genome_ok(matrix_deposited), _per_genome_ok(matrix_bfd)

    n_unpaired = {r: 0 for r in UNPAIRED_REASONS}
    n_unpaired[UNKNOWN_FLAG_KEY] = 0
    pairs = []
    for row in genomes.itertuples():
        lt = row.locustag
        if row.prior_status != "baseline_ok":
            n_unpaired["not_baseline"] += 1
        elif not (_scalar_true(row.has_proteins) and bfd_ok.get(lt, False)):
            n_unpaired["no_bfd_proteins"] += 1
            n_unpaired[UNKNOWN_FLAG_KEY] += int(pd.isna(row.has_proteins))
        elif not _scalar_true(row.has_deposited_proteins):
            n_unpaired["no_deposited_proteins"] += 1
            n_unpaired[UNKNOWN_FLAG_KEY] += int(pd.isna(row.has_deposited_proteins))
        elif not dep_ok.get(lt, False):
            n_unpaired["deposited_file_missing"] += 1
        else:
            pairs.append(lt)

    keep = ["locustag", "step_id", "copy_raw", "copy_filtered", "copy_unknown"]
    d = matrix_deposited[matrix_deposited.locustag.isin(pairs)][keep].rename(columns={
        "copy_raw": "copy_deposited", "copy_filtered": "copy_filtered_deposited",
        "copy_unknown": "unknown_deposited"})
    b = matrix_bfd[matrix_bfd.locustag.isin(pairs)][keep].rename(columns={
        "copy_raw": "copy_bfd", "copy_filtered": "copy_filtered_bfd", "copy_unknown": "unknown_bfd"})
    for nm, x in (("deposited", d), ("bfd", b)):
        if x.duplicated(["locustag", "step_id"]).any():
            raise ValueError(f"{nm} matrix has duplicate (locustag, step_id)")
    t = d.merge(b, on=["locustag", "step_id"], how="outer", validate="one_to_one", indicator=True)
    if (t["_merge"] != "both").any():
        raise ValueError("deposited and BFD matrices cover different (locustag, step_id) sets for paired genomes")
    t = t.drop(columns="_merge")
    t = t.merge(genomes[["locustag", "annotation_provider", "busco_complete"]], on="locustag",
                how="left", validate="many_to_one")
    t["delta"] = t.copy_bfd - t.copy_deposited
    both_known = (t.unknown_deposited == 0) & (t.unknown_bfd == 0)
    df = (t.copy_filtered_bfd - t.copy_filtered_deposited).astype("Int64")
    t["delta_filtered"] = df.where(both_known, pd.NA)
    t = t[TABLE_COLUMNS].reset_index(drop=True)
    t.attrs["n_unpaired"] = n_unpaired
    return (t, n_unpaired) if return_unpaired else t


def _boot_ci(delta: np.ndarray, n_boot: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(delta)
    means = delta[rng.integers(0, n, size=(n_boot, n))].mean(axis=1)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return float(lo), float(hi)


def _row(step, prov, delta, n_boot, seed):
    lo, hi = _boot_ci(delta, n_boot, seed)
    return dict(step_id=step, annotation_provider=prov, n_pairs=len(delta), mean_delta=float(delta.mean()),
                ci_low=lo, ci_high=hi, frac_bfd_more=float((delta > 0).mean()),
                frac_bfd_less=float((delta < 0).mean()))


def summarize_paired(table: pd.DataFrame, min_pairs: int = DEFAULT_MIN_PAIRS,
                     n_boot: int = DEFAULT_N_BOOT, seed: int = DEFAULT_SEED):
    """Return (summary, dropped_providers). See module docstring.

    Rows: one per (step, provider) with n_pairs >= min_pairs, plus one per step with provider "ALL".
    dropped_providers columns: provider, n_pairs (n_pairs = genomes of that provider in the table).
    """
    if min_pairs < 1 or n_boot < 1:
        raise ValueError("min_pairs and n_boot must be >= 1")
    t = table.copy()
    bad = sorted(set(t.annotation_provider.dropna()) & set(RESERVED_PROVIDERS))
    if bad:
        raise ValueError(f"annotation_provider value(s) {bad} are reserved (overall / null-provider labels)")
    t["annotation_provider"] = t.annotation_provider.fillna("(missing)").replace("", "(missing)")
    steps = list(dict.fromkeys(t.step_id))
    n_by_prov = t.groupby("annotation_provider").locustag.nunique()
    dropped = pd.DataFrame({"provider": [p for p, n in n_by_prov.items() if n < min_pairs],
                            "n_pairs": [int(n) for p, n in n_by_prov.items() if n < min_pairs]},
                           columns=["provider", "n_pairs"])
    rows = []
    for s in steps:
        x = t[t.step_id == s]
        for p in sorted(n_by_prov.index):
            if n_by_prov[p] >= min_pairs:
                rows.append(_row(s, p, x.delta[x.annotation_provider == p].to_numpy(dtype=float), n_boot, seed))
        rows.append(_row(s, "ALL", x.delta.to_numpy(dtype=float), n_boot, seed))
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS), dropped


def summarize_paired_by_busco(table: pd.DataFrame, bins=(0, 80, 95, 100)) -> pd.DataFrame:
    """Mean delta per (step, BUSCO bin). BUSCO control for spec 5.6.

    Bins are [b0,b1), [b1,b2), ... and the last is closed: [b_{n-1},b_n]. Null busco_complete goes to an
    explicit "unknown" bin. A value outside [bins[0], bins[-1]] raises ValueError. Empty bins are kept
    (n_pairs 0, mean_delta NaN).
    """
    edges = list(bins)
    labels = [f"[{edges[i]},{edges[i + 1]})" for i in range(len(edges) - 2)] + [f"[{edges[-2]},{edges[-1]}]"]
    b = table.busco_complete
    if (b.notna() & ((b < edges[0]) | (b > edges[-1]))).any():
        raise ValueError(f"busco_complete outside [{edges[0]}, {edges[-1]}]")
    idx = np.searchsorted(edges, b.to_numpy(dtype=float), side="right") - 1
    idx = np.minimum(idx, len(labels) - 1)  # the top edge belongs to the last bin
    t = table.assign(busco_bin=[labels[i] if not np.isnan(v) else "unknown" for i, v in zip(idx, b.to_numpy(dtype=float))])
    rows = []
    for s in dict.fromkeys(t.step_id):
        for lab in labels + ["unknown"]:
            x = t[(t.step_id == s) & (t.busco_bin == lab)]
            rows.append(dict(step_id=s, busco_bin=lab, n_pairs=len(x),
                             mean_delta=float(x.delta.mean()) if len(x) else float("nan")))
    return pd.DataFrame(rows, columns=["step_id", "busco_bin", "n_pairs", "mean_delta"])
