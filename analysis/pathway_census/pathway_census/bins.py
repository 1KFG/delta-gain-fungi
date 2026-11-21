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
The genome is the 8-hex LOCUSTAG prefix of protein_id (text before the first "_").

Bin order of decision (first match wins):
1. excluded. Any REQUIRED flag is known bad: truncated True, contam_suspect True, domain_ok False,
   or motif_ok False. Required flags default to (truncated, contam_suspect, domain_ok, motif_ok)
   and are a parameter (`require`). A known bad flag excludes even if another flag is NA.
2. unknown. No required flag is known bad and at least one required flag is NA. These copies are
   never promoted or demoted.
3. Otherwise count the four promotion terms. Each term is a boolean output column.
   - divergent          ref_identity < divergent_max_identity. NA does not promote.
   - lineage_expansion  copy_filtered >= expansion_min.
   - in_cluster         another assigned copy of a DIFFERENT step, same genome and same chrom,
                        whose interval lies within window_bp of this copy (gap between the two
                        intervals <= window_bp; overlap counts as gap 0). Copies of the same
                        step never count. The neighbour may itself be excluded or unknown: it
                        is still a pathway gene on the scaffold.
   - expressible        min_length <= length <= max_length; and, only for steps in
                        `soluble_steps`, no signal peptide and tm_helices <= tm_max. NA tm_helices
                        for a soluble step does not promote. `soluble_steps` is a parameter
                        (default empty): the pathway YAML has no soluble flag yet.
   tier1: at least tier1_min terms true. tier2: at least tier2_min true, fewer than tier1_min.
   tier3: fewer than tier2_min true (including zero).
4. Demotion: a tier1 copy with hap_dup_suspect True becomes tier2, reason
   "hap_dup_suspect_demoted". hap_dup_suspect NA changes nothing; reason "hap_dup_unevaluated"
   is shown. hap_dup_suspect is not a required flag by default.

Output columns: protein_id, step_id, bin, reasons, divergent, lineage_expansion, in_cluster,
expressible. `bin` is one of tier1, tier2, tier3, excluded, unknown. `reasons` is ";"-joined in
this fixed order: exclusion reasons (truncated, contam_suspect, domain_not_ok, motif_not_ok),
unevaluated flags ("unevaluated:<flag>", in `require` order), fired terms (divergent,
lineage_expansion, in_cluster, expressible), then hap_dup_suspect_demoted / hap_dup_unevaluated.
Term columns are filled for every row, but `reasons` lists terms only for tiered rows.
Output order = input order of the assigned rows (stable). Same input gives the same output.

PROPOSED defaults (no data basis):
  window_bp=20000. The cluster span of carotenoid genes in the reference species (Fusarium,
      Neurospora) is UNVERIFIED. Check it before any real run.
  tier1_min=3, tier2_min=2: tier1 needs most terms, tier2 needs a half. Chosen for simplicity.
  divergent_max_identity=0.5, expansion_min=2, min_length=100, max_length=1500, tm_max=1.
"""
from __future__ import annotations

import pandas as pd

DEFAULT_REQUIRE = ("truncated", "contam_suspect", "domain_ok", "motif_ok")


def triage(copies: pd.DataFrame, window_bp: int = 20000, **params) -> pd.DataFrame:
    raise NotImplementedError("rules fixed first; implementation follows in a later commit")
