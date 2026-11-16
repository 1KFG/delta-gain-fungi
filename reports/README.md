# `reports/` — the website source

This folder is **Quarto source for the public site only**. It renders to `../docs/reports/`,
which GitHub Pages serves at <https://dgf.fungalgenomes.org/reports/>.

## What belongs here

- `_quarto.yml` — the site config (navbar, theme, output directory).
- `index.qmd` — the report index; the navbar's landing page.
- One `<analysis>.qmd` per published report.

A report `.qmd` is a **thin wrapper, not a copy**. It carries the front matter and any
framing callout, then pulls in the report body with

```
{{< include ../analysis/<analysis>/REPORT.md >}}
```

so the published page and the `REPORT.md` that lives beside the code cannot drift. Write
the analysis in `analysis/<name>/REPORT.md`; write only the wrapper here.

## What does not belong here

- **Data files and analysis output.** `reports/` is not an output directory. Results stay
  under `analysis/<name>/out/` (gitignored) or on lab storage. A 2,000-row TSV lived here
  for a while because `scripts/dg_rebin.py` writes to a `reports/` directory *inside the
  shared DeltaGain folder on /bigdata* — a different `reports/` that happens to share the
  name. Nothing in the repository ever read the copy that landed here.
- **Rendered HTML.** That is `docs/reports/`, written by `quarto render` and committed.
  Do not hand-edit it.
- **The prose of a report.** See above: it goes in `analysis/<name>/REPORT.md`.

## Rendering

```sh
cd reports && quarto render      # -> ../docs/reports/
```

Quarto warns that the output directory is not a subdirectory of the project directory, and
that it is therefore refusing to clean `docs/reports/site_libs`. Both are expected: the
output is deliberately outside this folder so that Pages can serve it from `docs/`
alongside the hand-written Jekyll landing page at `docs/index.md`. Nothing is lost —
`site_libs` is simply left in place and overwritten rather than deleted first.

Commit the re-rendered `docs/reports/` together with the source change. The live site only
updates when the result reaches `main`.
