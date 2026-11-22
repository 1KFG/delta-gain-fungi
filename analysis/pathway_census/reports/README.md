# Pathway census reports

One Quarto report per pathway. Format rules: `WRITING.md` at the repository root (twelve rules, callout, caption and figure idioms).
Layout, as in `analysis/enzymes_UPO`:

```
analysis/pathway_census/reports/<pathway>/REPORT.md        source of truth (prose, tables, figure links)
analysis/pathway_census/reports/<pathway>/figures/bin/     committed, tested figure code
analysis/pathway_census/reports/<pathway>/figures/out/     PNG (page) and PDF (linked from the caption)
reports/pathway_<pathway>/index.qmd                        thin wrapper: front matter + {{< include >}}
reports/pathway_<pathway>/figures -> ../../analysis/pathway_census/reports/<pathway>/figures
```

The H1 of `REPORT.md` and the `title:` of the wrapper must be identical (the site filter drops the duplicate).

## Section order used by the carotenoid report (the order may change)

1. Opening callout: what the report is, what it is not, where it stops being reliable.
2. Findings so far: bold the claim, give the number with its denominator and its bound.
3. The pathway, references, decoys and hits per step (table written by `bin/make_step_table.py`).
4. Evidence sections, each with a figure whose caption says what it is evidence for and names the misreading:
   score separation (step 2), genome-wide search (step 4), cutoff side evaluation, and, when done, placement (step 5),
   copy numbers by phylum (step 8), deposited versus BFD (step 10), absence check (step 11), triage (step 12).
5. Negative and unresolved results, in their own section.
6. What would change this, and what was not done.
7. Methods, tools, versions, job ids, commits and files.

## Making a report for another pathway

1. Finish census steps 1 to 4 for the pathway (`analysis/pathway_census/README.md`, "Batched BFD search").
2. Copy `reports/carotenoid/` to `reports/<pathway>/`. `figures/bin/fig_census_search.py` and `bin/make_step_table.py` take the census folders as arguments and need no code change.
3. Run the figure script and the table script (commands in the carotenoid `REPORT.md`, section 8), then rewrite the prose from the new numbers. Do not reuse a sentence without checking its number.
4. Add the wrapper and the `figures` symlink. Render in a scratch copy, serve it over HTTP and look at it before touching `reports/_quarto.yml` or `docs/`.
5. Colours: first three slots of the validated categorical palette (blue `#2a78d6`, orange `#eb6834`, aqua `#1baf7a`). Aqua is below 3:1 contrast on the light surface, so every series has a direct label and its own marker shape.

## Status

| Pathway | Report | State |
|---|---|---|
| carotenoid | `reports/carotenoid/REPORT.md` | draft, steps 1 to 4 (2026-11-22) |
| melanin_dhn, melanin_ldopa, mycosporine | none | steps 1 to 4 not run |
