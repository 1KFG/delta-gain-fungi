"""Quality, contamination and haplotype-duplication flags (spec section 2 anti-goals, section 5.4 step 4).

Output: one row per input (protein_id, step_id) row of `assigned`, with columns
protein_id, step_id, domain_ok, motif_ok, truncated, contam_suspect, hap_dup_suspect.
Flags are only flags. No row is dropped here. All flag columns use the nullable boolean dtype.
<NA> means "not evaluated because the needed data was not supplied". <NA> is never False.

All thresholds are function parameters. Every value below is PROPOSED. None is MEASURED. No census
data exists yet to calibrate them.

Rules (defined before coding)
-----------------------------
Gene id of a protein: protein_id without the transcript suffix `-T<n>` (for example
`F77FFF1D_000001-T1` -> `F77FFF1D_000001`). The genome of a protein is
the LOCUSTAG prefix before the first underscore (8 hex characters in all but 3 BFD genomes). gene_info is keyed by that gene id. Its columns are
gene_id, LOCUSTAG, chrom, start, end, plus `length` (protein length in aa, from BFD
`gene_proteins.parquet`, joined by the caller).

domain_ok
    True if every expected domain of the step (`expected_domains[step_id]`, hmm names such as Pfam
    ids) occurs in the `domains` frame for that protein with evalue <= domain_evalue
    (PROPOSED 1e-5). Counting is at domain level, one hit table row per domain. A fused gene is
    judged per step: it needs the domains of that step only, not the domains of its fusion partner.
    <NA> if `expected_domains` has no entry for the step.

motif_ok
    The reference definition gives no catalytic motif yet. The value comes from `motif_hits`
    (protein_id, step_id, motif_present) when the caller supplies it. Otherwise <NA>.

truncated
    True if protein length < length_floor * refs_len[step_id][0] (lower bound of the reference
    length range; length_floor PROPOSED 0.7). Equality is not truncated. <NA> if length or the
    reference range is missing.

contam_suspect
    Data source: BFD has no per-gene taxonomy table (checked in tables/ on 2026-11-20). The flag uses
    a precomputed `best_hits` table (protein_id, taxon_class in {"fungal", "non_fungal"}). The
    caller makes it by searching proteins against a reference set labelled fungal vs
    bacterial/algal. This module does not run that search.
    Rule: take the genes on the same scaffold (same LOCUSTAG and chrom) within `neighbour_window`
    genes of the focal gene on each side, ordered by start (PROPOSED 10), excluding the focal gene.
    Let n be the neighbours that have a best hit. If n >= min_neighbours_with_hits (PROPOSED 3) and
    the fraction of those n that are "non_fungal" >= non_fungal_fraction (PROPOSED 0.5), the flag is
    True. Fewer than min_neighbours_with_hits: <NA> (insufficient evidence, unknown, NOT clean).
    No `best_hits` table, or focal gene not in gene_info: <NA>. The focal gene's own best hit is NOT used. A true carotenoid gene from a lichen
    or algal partner would hit the foreign class itself, so the scaffold context decides.

hap_dup_suspect
    Data source: `pair_identity` table (protein_a, protein_b, pident as a fraction 0-1, qcov as a
    fraction 0-1), made by the caller from an all-vs-all search inside one genome. It must cover the
    neighbour genes, not only the pathway genes.
    Rule: two proteins of the same step and the same LOCUSTAG, on DIFFERENT scaffolds, with
    pident >= hap_identity (PROPOSED 0.98) and qcov >= hap_coverage (PROPOSED 0.9), are a
    candidate pair. The pair is suspect (both members flagged) if the scaffolds also duplicate each
    other around the pair: at least min_shared_neighbours (PROPOSED 3) genes within
    `neighbour_window` genes of one member have a near-identical partner (same two thresholds)
    within `neighbour_window` genes of the other member. Same-scaffold near-identical pairs are
    tandem copies and are not flagged. No `pair_identity` table, or gene not in gene_info: <NA>.
"""
from __future__ import annotations

import re

import pandas as pd

OUT_COLUMNS = ["protein_id", "step_id", "domain_ok", "motif_ok", "truncated", "contam_suspect",
               "hap_dup_suspect"]
_TSUFFIX = re.compile(r"-T\d+$")


def _gid(protein_id: str) -> str:
    return _TSUFFIX.sub("", protein_id)


def _neighbours(gene_info: pd.DataFrame, window: int) -> dict[str, list[str]]:
    """gene_id -> gene ids within `window` genes on each side, same LOCUSTAG and chrom, by start."""
    out: dict[str, list[str]] = {}
    for _, g in gene_info.groupby(["LOCUSTAG", "chrom"], sort=False):
        ids = g.sort_values("start").gene_id.tolist()
        for i, gid in enumerate(ids):
            out[gid] = ids[max(0, i - window):i] + ids[i + 1:i + 1 + window]
    return out


def _bool_series(values: list) -> pd.Series:
    return pd.Series(values, dtype="boolean")


def quality_flags(assigned: pd.DataFrame, domains: pd.DataFrame, gene_info: pd.DataFrame,
                  refs_len: dict[str, tuple[int, int]], genomes: pd.DataFrame, *,
                  expected_domains: dict[str, list[str]] | None = None,
                  motif_hits: pd.DataFrame | None = None,
                  best_hits: pd.DataFrame | None = None,
                  pair_identity: pd.DataFrame | None = None,
                  domain_evalue: float = 1e-5,
                  length_floor: float = 0.7,
                  neighbour_window: int = 10,
                  min_neighbours_with_hits: int = 3,
                  non_fungal_fraction: float = 0.5,
                  hap_identity: float = 0.98,
                  hap_coverage: float = 0.9,
                  min_shared_neighbours: int = 3) -> pd.DataFrame:
    """Return the flag table described in the module docstring. All thresholds are PROPOSED."""
    known = set(genomes.locustag)
    unknown = sorted({p.split("_")[0] for p in assigned.protein_id} - known)
    if unknown:
        raise ValueError(f"protein LOCUSTAG(s) not in genomes table: {unknown[:5]}")

    info = gene_info.set_index("gene_id")
    nbrs = _neighbours(gene_info, neighbour_window)

    # domain level evidence: (protein, hmm_name) with evalue <= threshold
    good = domains[domains.evalue <= domain_evalue] if len(domains) else domains
    dom_have = set(zip(good.protein_id, good.hmm_name)) if len(good) else set()
    motif = {}
    if motif_hits is not None:
        motif = {(r.protein_id, r.step_id): bool(r.motif_present) for r in motif_hits.itertuples()}
    cls = dict(zip(best_hits.protein_id.map(_gid), best_hits.taxon_class)) if best_hits is not None else None

    # haplotype pair structures
    partners: dict[str, set[str]] = {}
    if pair_identity is not None:
        ok = pair_identity[(pair_identity.pident >= hap_identity) & (pair_identity.qcov >= hap_coverage)]
        for a, b in zip(ok.protein_a.map(_gid), ok.protein_b.map(_gid)):
            partners.setdefault(a, set()).add(b)
            partners.setdefault(b, set()).add(a)
    by_step: dict[str, set[str]] = {}
    for r in assigned.itertuples():
        by_step.setdefault(r.step_id, set()).add(_gid(r.protein_id))

    def hap(gid: str, step: str) -> bool:
        for other in partners.get(gid, ()):
            if other not in by_step[step] or other not in info.index or gid not in info.index:
                continue
            if info.at[other, "LOCUSTAG"] != info.at[gid, "LOCUSTAG"]:
                continue  # different genome
            if info.at[other, "chrom"] == info.at[gid, "chrom"]:
                continue  # same scaffold: tandem copy
            other_nb = set(nbrs.get(other, ()))
            shared = sum(1 for n in nbrs.get(gid, ()) if partners.get(n, set()) & other_nb)
            if shared >= min_shared_neighbours:
                return True
        return False

    dom_ok, mot_ok, trunc, contam, hapf = [], [], [], [], []
    for r in assigned.itertuples():
        pid, step, gid = r.protein_id, r.step_id, _gid(r.protein_id)
        exp = (expected_domains or {}).get(step)
        dom_ok.append(pd.NA if exp is None else all((pid, d) in dom_have for d in exp))
        mot_ok.append(motif.get((pid, step), pd.NA))

        length = info.at[gid, "length"] if gid in info.index and "length" in info.columns else pd.NA
        rng = refs_len.get(step)
        trunc.append(pd.NA if pd.isna(length) or rng is None else bool(length < length_floor * rng[0]))

        if cls is None or gid not in nbrs:
            contam.append(pd.NA)
        else:
            c = [cls[n] for n in nbrs[gid] if n in cls]
            if len(c) < min_neighbours_with_hits:
                contam.append(pd.NA)  # insufficient evidence: unknown, not clean
            else:
                contam.append(bool(sum(x == "non_fungal" for x in c) / len(c) >= non_fungal_fraction))

        hapf.append(pd.NA if pair_identity is None or gid not in info.index else hap(gid, step))

    return pd.DataFrame({
        "protein_id": assigned.protein_id.values, "step_id": assigned.step_id.values,
        "domain_ok": _bool_series(dom_ok), "motif_ok": _bool_series(mot_ok),
        "truncated": _bool_series(trunc), "contam_suspect": _bool_series(contam),
        "hap_dup_suspect": _bool_series(hapf)})[OUT_COLUMNS]


FLAG_COLUMNS = ["domain_ok", "motif_ok", "truncated", "contam_suspect", "hap_dup_suspect"]


def clean_mask(flags: pd.DataFrame,
               require: tuple[str, ...] = ("truncated", "contam_suspect", "hap_dup_suspect")) -> pd.Series:
    """True only where every required flag is known False. NA or True gives False.

    Use this to filter. Never use fillna(False): that treats unknown as clean.
    """
    ok = pd.Series(True, index=flags.index)
    for col in require:
        ok &= flags[col].eq(False).fillna(False).astype(bool)
    return ok


def unknown_counts(flags: pd.DataFrame) -> pd.Series:
    """Per flag column, the number of NA (not evaluated) rows."""
    cols = [c for c in FLAG_COLUMNS if c in flags.columns]
    return flags[cols].isna().sum().astype(int)
