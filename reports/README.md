# `reports/` — the website source

This folder is **Quarto source for the public site only**. It renders to `../docs/reports/`,
which GitHub Pages serves at <https://dgf.fungalgenomes.org/reports/>.

## What belongs here

- `_quarto.yml` — the site config (navbar, theme, output directory).
- `index.qmd` — the report index; the navbar's landing page.
- One wrapper per published report: either `<analysis>.qmd`, or `<analysis>/index.qmd`
  if the report has figures (see below).
- `_filters/` — Lua filters applied to every page. The leading `_` keeps Quarto from
  treating them as pages.

A report `.qmd` is a **thin wrapper, not a copy**. It carries the front matter and any
framing callout, then pulls in the report body with

```
{{< include ../analysis/<analysis>/REPORT.md >}}
```

so the published page and the `REPORT.md` that lives beside the code cannot drift. Write
the analysis in `analysis/<name>/REPORT.md`; write only the wrapper here.

### Reports with figures

A figure has to resolve from two places at once: from `analysis/<name>/REPORT.md` when
GitHub renders it, and from `reports/` when Quarto renders it. `{{< include >}}` is a
textual include and does not rewrite paths, so a path that works in one breaks in the
other — Quarto passes an unresolvable `src` straight through and the image 404s on the
site with no warning at render time.

The fix is to make the same relative path valid from both ends. `REPORT.md` keeps the
GitHub-correct path (`figures/out/fig1.png`), and the report gets **its own subdirectory
here** holding a symlink that makes that path resolve:

```
reports/upo_clades/index.qmd            # wrapper, includes ../../analysis/enzymes_UPO/REPORT.md
reports/upo_clades/figures -> ../../analysis/enzymes_UPO/figures
```

Quarto follows the symlink and copies the referenced files into
`docs/reports/upo_clades/figures/out/` — PDFs linked from captions as well as the PNGs
themselves. The subdirectory is what keeps this safe: the symlink is scoped to the one
report that owns it, so two reports can both have a `figures/` without colliding. The
page is then served at `/reports/<name>/`.

A report with no figures does not need any of this and can stay a flat `<name>.qmd`.

### Reports whose source is another repository

`reports/delta_gain/` publishes `docs/FINDINGS.md` from
[1KFG/delta-gain](https://github.com/1KFG/delta-gain), which is where the funnel pipeline and
its narrative live. Copying that file here would create a second copy to keep in step with
the first, so the page reproduces it at a **pinned commit** instead, fetched by

```sh
reports/bin/fetch_external.sh      # -> reports/delta_gain/{_external,figures}/, gitignored
```

Both the commit and each file's sha256 are pinned in that script: if the content behind the
pinned URL ever differs, the fetch fails and writes nothing, because the page states that it
reproduces the document verbatim. Re-fetch before re-rendering. The rendered output under
`docs/reports/` is committed, so the live site never depends on the fetch.

Because the reproduction is verbatim, anything we need to say about that document goes in
the wrapper's callouts, not into the text — see the qualification box on that page.

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

A local render is not a check that the page is right. Serve `docs/` over HTTP and look at
it — `python3 -m http.server` from `docs/`, then open `/reports/<name>/`. Opening the HTML
as a `file://` URL reports CORS errors for `site_libs/*.js` that do not happen on the real
site, and will not tell you whether a figure actually loaded.
