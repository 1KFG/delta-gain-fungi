# How reports in this repository are written

For `analysis/<name>/REPORT.md`, the pages built from them, and the pull requests that
publish them.

## The look is stock. Do not go looking for the stylesheet.

The site is Quarto's `cosmo` theme with `toc: true` and **no custom CSS or SCSS at all**.
Nothing here is a design system and nothing needs one. What makes these pages read the way
they do is the structure below, which is editorial, costs nothing, and works just as well in
a plain Markdown file, an email or a PR body as it does in HTML.

If a page looks wrong, the fix is almost always a structural one from this document, not a
styling one.

---

## The rules that actually do the work

### 1. One source of truth, included — never copied

The prose lives in `analysis/<name>/REPORT.md`, beside the code that produced it. The page is
a thin wrapper that pulls it in:

```
{{< include ../../analysis/<name>/REPORT.md >}}
```

A copy is a second thing to keep in step, and it will drift. If a document belongs to another
repository, reproduce it at a pinned commit with a pinned checksum rather than copying it in
(`reports/bin/fetch_external.sh`).

### 2. Say what the thing is — and is not — before the reader can over-read it

Every page opens with a callout that states, in order: that it is a working analysis, that
nothing has been tested at the bench, and **where the work stops being reliable**. That box is
the single highest-value element on the page. It goes above the content, not at the end,
because a limitation discovered after the reader has formed a conclusion is a limitation that
did not do its job.

### 3. Lead with the number, and put the limit in the same breath

> Ranking clade-V-like pockets is reliable — LOO AUC 0.956. Predicting whether an arbitrary
> protein is an active peroxygenase is not: out-of-clade AUC **0.556**.

Not "performs well in some contexts." A claim without its number is an opinion; a number
without its bound is a trap.

### 4. Negative results get sections, not footnotes

"The negative results" is a numbered section with the four failed attempts in it, each with
the statistic that killed it. A result that came out negative cost exactly as much to obtain
as one that came out positive, and burying it means someone pays for it again. The same goes
for anything discarded by a pre-declared kill switch: say that it was declared in advance.

### 5. State what would change the conclusion, and what was not done

Every report carries a "what would change this" and an explicit "not done here." Where a test
would settle an open question, specify the test — sample sizes, outgroups, what each outcome
would mean — even when you are not running it. An unresolved question that names its own
decisive experiment is a contribution. One that trails off is a liability.

### 6. Resolve the reading traps instead of passing them on

If two tables in the same document can be read as contradicting each other, say so in the
document and explain why they differ. Two examples from `enzymes_UPO/REPORT.md`:

- Table T1 and the trust-tier table count differently, so Clade I reads 41,679 in one and
  2,565 in the other.
- Clade I and Zamocky's "Basidiomycetous 3rd" are the same group under two names.

Nobody reading cold will work either out, and both look like errors. This is the single thing
most worth copying from these reports.

### 7. A caption says what the figure is evidence for

Not what it depicts — the reader can see that. State what the figure demonstrates, and warn
about the available misreading:

> **This figure counts the way Table T1 does, not the way the tier table does**: "unassigned"
> here means *placed outside Clades I–V* (41,001), not *below the 0.5 likelihood-weight
> threshold* (76,913).

Where a figure's n is small, the caption says what that buys and what it does not: *"the
intervals are wide because n = 19; this shows the caller abstains rather than guesses, not
that it is accurate to within a few percent."*

### 8. Put a correction where it cannot falsify the thing it corrects

When a reproduced document's result is later qualified, the qualification goes in a callout
on the wrapper — never into the body, because the page claims to reproduce that body
verbatim and an edit would make the claim false. The callout names the affected section,
gives the corrected magnitude, says explicitly what the new test does *not* establish, and
says which parts of the document are untouched.

### 9. Numbers carry their denominator, their n, and their direction

`16/16`, `95/95`, `2,265 proteins in 2,190 genomes`, `1,033 genomes = 5.8 % of 17,827`.
A bare percentage hides its sample size. Give the absolute and the fraction together when the
absolute is small or the denominator is surprising.

### 10. Bold marks the claim, not the emphasis

One bolded span per paragraph at most, on the sentence a skimming reader must not miss.
Bolding for emphasis spends the budget and nothing stands out.

### 11. Link to the exact thing

`G1_FINDINGS.md §5b`, `PLAN.md:74`, `scripts/dg_rebin.py` line 48, a commit SHA. "See the
documentation" is not a reference.

### 12. Report outcomes faithfully

If a job failed, say so with the output. If a step was skipped, say that. If something is
done and verified, state it plainly without hedging. Do not describe predicted results in
language that implies measurement — "candidates for assay, not a prediction of activity."

---

## The index-page contract

The reports index carries two sections and both are load-bearing:

- **Available** — each report with a few sentences that include its main limitation, not a
  teaser.
- **Not published here** — what exists in the repository with no page, **and why**. A reader
  who can see what was withheld and the reason can calibrate everything else. "Its main
  result was fixed upstream the same day, so its tables describe a snapshot being rebuilt"
  tells them more about the project's standards than any published page does.

---

## Page mechanics

| need | mechanism |
|---|---|
| page from a report | thin `.qmd` wrapper + `{{< include >}}` |
| report with figures | own subdirectory + scoped `figures` symlink — see `reports/README.md` |
| document from another repo | `reports/bin/fetch_external.sh`, commit **and** sha256 pinned |
| title rendering twice | `reports/_filters/drop-duplicate-title.lua` |
| checking the result | serve `docs/` over HTTP and look at it — `file://` lies |

---

## Before publishing

- [ ] Opening callout states what it is, what it is not, and where it stops being reliable
- [ ] Every headline number carries its bound, its n and its denominator
- [ ] Negative results are in a section of their own
- [ ] Any two tables that could be read as contradicting each other are reconciled in the text
- [ ] Each figure caption says what the figure is evidence for, and names the misreading
- [ ] "Not done here" is explicit, and names the test that would settle it
- [ ] Every reference points at a file, section or line
- [ ] Rendered, served over HTTP, and **looked at** — images loaded, no dead links, no console errors
- [ ] The index's "Not published here" is still true
