# Open side-questions (return to these)

## Q1. Has lycopene cyclase diversified in fungi, and in which groups? (owner, 2026-11-22)
- Context: the fungal lycopene cyclase is the cyclase domain of a fused phytoene synthase / cyclase protein (CrtYB type, domain rule any-of TIGR03462, PF18916). Five references exist in `refs/carotenoid/lycopene_cyclase.faa`.
- Question: do copies of the cyclase domain differ in domain architecture (fused, stand-alone, duplicated) or in clade across fungal groups? Which groups show expansion or loss?
- Needs: the carotenoid census (steps 4 to 8) and the step 5 tree for `lycopene_cyclase`. Not testable before that.
- Status: not started. No data.
- Related: ortholog naming (separate Stage 2 task).

## Q2. A Saccharomycotina-specific BUSCO marker set (owner, 2026-11-23)
- Context: BFD ran `fungi_odb12` (1,122 markers) for 23,582 genomes. Saccharomycetes has median completeness 75.4% (10th to 90th percentile 72.2 to 77.9, 2,589 genomes) against 99.3% for Sordariomycetes. A narrow band at a low value fits markers missing from the lineage, not assembly quality (UNVERIFIED).
- Effect on the census: a BUSCO filter at 80 or higher removes the budding yeasts. A floor of 65 to 70 was tried as a sensitivity check (it changes counts by 13% or less and cannot catch an incomplete yeast assembly). Owner decision, 2026-11-23: no BUSCO filter in the main gain/loss results.
- Decision (owner): do not run a yeast BUSCO set now. Add a Saccharomycotina-specific run (about 6,035 genomes) to the next BFD BUSCO run, because the earlier run was expensive and should not be doubled.
- Note written in the BFD run folder: `Fungi_BFD_runs/docs/NOTES_TO_FOLLOW.md`, section "BUSCO marker set for Saccharomycotina". It is not committed in that checkout.
- Status: note only; nothing run; dataset name and cost not checked.
