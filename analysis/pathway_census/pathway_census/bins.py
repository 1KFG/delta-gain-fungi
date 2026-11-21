"""Triage bins for assigned pathway copies (spec section 5.5).

THE BINS ARE TRIAGE, NOT A RANKING AND NOT AN ACTIVITY PREDICTION. They say which copies to look
at first. They say nothing about enzyme activity. No weights are fitted. There is no additive
score. These rules were written down before any real data was binned (Task 10, brief step 1).
Every number below is PROPOSED: it has no data basis. Nothing here is MEASURED.

Input (`triage(copies, ...)`): one row per (protein_id, step_id), already joined by the caller.
Only rows with label == "assigned" are used. Rows with another label are dropped and never appear
in the output. Columns:
  protein_id, step_id, label
  truncated, contam_suspect, domain_ok, motif_ok, hap_dup_suspect
                  pandas nullable booleans from quality.quality_flags. NA = not evaluated.
  ref_identity    best identity to a reference protein of the step, as a fraction 0..1.
                  NA = not computed.
  copy_filtered   the genome's clean copy count for the step, from matrix.build_matrix.
  chrom, start, end   gene coordinates from gene_info (same scaffold name = same chrom).
  length          protein length in aa.
  signal_peptide  boolean. signalp is POSITIVES ONLY: a protein absent from the signalp table is
                  passed as False ("no signal peptide predicted"). ASSUMPTION: the table lists
                  every protein that was run, and absence means a negative call. Not verified.
  tm_helices      integer count of predicted helices from tmhmm. NA = not computed.
The genome is the LOCUSTAG prefix of protein_id: the text before the first "_" (8 hex characters
in all but 3 BFD genomes).

Bin order of decision (first match wins):
1. excluded. Any REQUIRED flag is known bad: truncated True, contam_suspect True, domain_ok False,
   or motif_ok False. Required flags default to (truncated, contam_suspect, domain_ok, motif_ok)
   and are a parameter (`require`). A known bad flag excludes even if another flag is NA.
2. unknown. No required flag is known bad and at least one required flag is NA. These copies are
   never promoted or demoted.
3. Otherwise count the four promotion terms. Each term is a boolean output column.
   - divergent          ref_identity < divergent_max_identity. NA does not promote.
   - lineage_expansion  copy_filtered >= expansion_min.
   - in_cluster         another assigned gene (a different protein_id AND a different step_id),
                        same genome and same chrom (a gene never neighbours itself, so a fused
                        protein assigned to two steps is not its own neighbour; an empty or NA
                        chrom never neighbours),
                        whose interval lies within window_bp of this copy (gap between the two
                        intervals <= window_bp, inclusive; overlap counts as gap 0). Copies of
                        the same step never count. The neighbour may itself be excluded or unknown: it
                        is still a pathway gene on the scaffold.
   - expressible        min_length <= length <= max_length; and, only for steps in
                        `soluble_steps`, no signal peptide and tm_helices <= tm_max. NA tm_helices
                        for a soluble step does not promote. `soluble_steps` is a parameter
                        (default empty): the pathway YAML has no soluble flag yet.
   The term columns are plain booleans for counting: NA (not evaluated) counts as not true.
   A term is "unevaluated" when its input is NA: divergent (ref_identity), lineage_expansion
   (copy_filtered), in_cluster (chrom empty/NA, or start/end NA), expressible (length NA, or for a
   soluble step tm_helices or signal_peptide NA).
   tier1: at least tier1_min terms true. tier2: at least tier2_min true, fewer than tier1_min.
   tier3: fewer than tier2_min true (including zero).
4. Demotion: a tier1 copy with hap_dup_suspect True becomes tier2, reason
   "hap_dup_suspect_demoted". hap_dup_suspect NA changes nothing; reason "hap_dup_unevaluated"
   is shown. hap_dup_suspect is not a required flag by default.

Output columns: protein_id, step_id, bin, reasons, divergent, lineage_expansion, in_cluster,
expressible. `bin` is one of tier1, tier2, tier3, excluded, unknown. `reasons` is ";"-joined in
this fixed order: exclusion reasons (truncated, contam_suspect, domain_not_ok, motif_not_ok),
unevaluated flags ("unevaluated:<flag>", in `require` order), fired terms (divergent,
lineage_expansion, in_cluster, expressible), unevaluated terms ("unevaluated:<term>", same term
order), then hap_dup_suspect_demoted / hap_dup_unevaluated.
Term columns are filled for every row, but `reasons` lists terms only for tiered rows.
`n_terms_unevaluated` (int) counts unevaluated terms for tiered rows and is 0 for excluded and
unknown rows.
Output is sorted by (protein_id, step_id), so row order does not depend on input order. Duplicate
(protein_id, step_id) among assigned rows raise ValueError. Required flag columns must be pandas
boolean, or object holding only True/False/pd.NA/None; anything else (np.nan in an object column,
strings, numbers) raises ValueError.

PROPOSED defaults (no data basis):
  window_bp=20000. The cluster span of carotenoid genes in the reference species (Fusarium,
      Neurospora) is UNVERIFIED. Check it before any real run.
  tier1_min=3, tier2_min=2: tier1 needs most terms, tier2 needs a half. Chosen for simplicity.
  divergent_max_identity=0.5, expansion_min=2, min_length=100, max_length=1500, tm_max=1.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

DEFAULT_REQUIRE = ("truncated", "contam_suspect", "domain_ok", "motif_ok")


TERMS = ("divergent", "lineage_expansion", "in_cluster", "expressible")
_BAD = {"truncated": (True, "truncated"), "contam_suspect": (True, "contam_suspect"),
        "domain_ok": (False, "domain_not_ok"), "motif_ok": (False, "motif_not_ok")}
_EXCL_ORDER = ("truncated", "contam_suspect", "domain_ok", "motif_ok")
_NEEDED = ("protein_id", "step_id", "label", "ref_identity", "copy_filtered", "chrom", "start", "end",
           "length", "signal_peptide", "tm_helices", "hap_dup_suspect")


def _na(v) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and v != v)


def _check_flag_col(df: pd.DataFrame, col: str) -> None:
    s = df[col]
    if str(s.dtype) in ("boolean", "bool"):
        return
    if s.dtype == object:
        for v in s:
            if v is None or v is pd.NA or isinstance(v, (bool, np.bool_)):
                continue
            raise ValueError(f"column '{col}': invalid value {v!r} (need True/False/pd.NA/None)")
        return
    raise ValueError(f"column '{col}': dtype {s.dtype} is not boolean or object of True/False/NA")


def _known(v):
    """None if NA, else bool."""
    return None if _na(v) else bool(v)


def _in_cluster(df: pd.DataFrame, window_bp: int) -> tuple[list[bool], list[bool]]:
    """Return (in_cluster, evaluated). A gene never neighbours itself (same protein_id)."""
    genome = df.protein_id.str.split("_").str[0].tolist()
    pid, step = df.protein_id.tolist(), df.step_id.tolist()
    chrom, start, end = df.chrom.tolist(), df.start.tolist(), df.end.tolist()
    ok = [not (_na(chrom[i]) or chrom[i] == "" or _na(start[i]) or _na(end[i])) for i in range(len(df))]
    by_scaf: dict = {}
    for i in range(len(df)):
        if ok[i]:
            by_scaf.setdefault((genome[i], chrom[i]), []).append(i)
    out = [False] * len(df)
    for idx in by_scaf.values():
        for i in idx:
            for j in idx:
                if step[j] == step[i] or pid[j] == pid[i]:
                    continue
                if max(start[i], start[j]) - min(end[i], end[j]) <= window_bp:
                    out[i] = True
                    break
    return out, ok


def triage(copies: pd.DataFrame, window_bp: int = 20000, *, tier1_min: int = 3, tier2_min: int = 2,
           divergent_max_identity: float = 0.5, expansion_min: int = 2, min_length: int = 100,
           max_length: int = 1500, tm_max: int = 1, soluble_steps=(),
           require=DEFAULT_REQUIRE) -> pd.DataFrame:
    """Bin assigned copies. See the module docstring for the rules. All parameters are PROPOSED."""
    missing = [c for c in _NEEDED + tuple(require) if c not in copies.columns]
    if missing:
        raise ValueError(f"triage input is missing columns: {missing}")
    if tier2_min > tier1_min:
        raise ValueError("tier2_min must be <= tier1_min")
    for c in dict.fromkeys(tuple(require) + ("hap_dup_suspect",)):
        _check_flag_col(copies, c)
    df = copies[copies.label == "assigned"]
    dup = df[df.duplicated(["protein_id", "step_id"], keep=False)]
    if len(dup):
        raise ValueError(f"duplicate (protein_id, step_id) rows: {dup[['protein_id', 'step_id']].drop_duplicates().values.tolist()[:5]}")
    df = df.sort_values(["protein_id", "step_id"], kind="mergesort").reset_index(drop=True)
    n = len(df)
    soluble = set(soluble_steps)
    cluster, cl_ok = _in_cluster(df, window_bp)
    rows = []
    for i in range(n):
        r = df.iloc[i]
        flags = {c: _known(r[c]) for c in require}
        bad = [_BAD[c][1] for c in _EXCL_ORDER if c in flags and flags[c] is _BAD[c][0]]
        unev = [f"unevaluated:{c}" for c in require if flags[c] is None]
        ident, ncopy, tm, ln = r.ref_identity, r.copy_filtered, r.tm_helices, r.length
        div = bool(not _na(ident) and ident < divergent_max_identity)
        exp = bool(not _na(ncopy) and ncopy >= expansion_min)
        expr = bool(not _na(ln) and min_length <= ln <= max_length)
        expr_ok = not _na(ln)
        if r.step_id in soluble:
            if _na(tm) or _na(r.signal_peptide):
                expr_ok = False
                expr = False
            elif expr:
                expr = (not bool(r.signal_peptide)) and tm <= tm_max
        terms = dict(divergent=div, lineage_expansion=exp, in_cluster=cluster[i], expressible=expr)
        evald = dict(divergent=not _na(ident), lineage_expansion=not _na(ncopy), in_cluster=cl_ok[i],
                     expressible=expr_ok)
        n_unev = 0
        if bad:
            bin_, reasons = "excluded", bad
        elif unev:
            bin_, reasons = "unknown", unev
        else:
            k = sum(terms.values())
            bin_ = "tier1" if k >= tier1_min else "tier2" if k >= tier2_min else "tier3"
            reasons = [t for t in TERMS if terms[t]]
            un = [t for t in TERMS if not evald[t]]
            n_unev = len(un)
            reasons += [f"unevaluated:{t}" for t in un]
            hap = _known(r.hap_dup_suspect)
            if hap is None:
                reasons.append("hap_dup_unevaluated")
            elif hap and bin_ == "tier1":
                bin_ = "tier2"
                reasons.append("hap_dup_suspect_demoted")
        rows.append(dict(protein_id=r.protein_id, step_id=r.step_id, bin=bin_, reasons=";".join(reasons),
                         **terms, n_terms_unevaluated=n_unev))
    cols = ["protein_id", "step_id", "bin", "reasons", *TERMS, "n_terms_unevaluated"]
    return pd.DataFrame(rows, columns=cols)
