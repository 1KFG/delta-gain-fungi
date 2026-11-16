# delta-gain-fungi

## Writing

Reports, published pages and PR bodies follow [`WRITING.md`](WRITING.md). Read it before
writing or editing an `analysis/<name>/REPORT.md`, a page under `reports/`, or a PR
description. The short version: the visual style is stock Quarto `cosmo` with no custom CSS,
and everything that matters is structural — limitations stated up front, every number with
its bound and denominator, negative results in their own section, reading traps resolved
rather than passed on, and "not done here" made explicit.

## Site

`reports/` is Quarto source for <https://dgf.fungalgenomes.org>; `docs/` is the published
output and is committed. See [`reports/README.md`](reports/README.md) for the folder
contract, the figures-symlink convention and the pinned external fetch. Render with
`cd reports && quarto render`, then serve `docs/` over HTTP and look at the page —
`file://` reports CORS errors that do not occur on the live site.

## Decisions that are the owner's

Merging pull requests, DNS and HTTPS settings, and deletions in `/bigdata/stajichlab/shared/`
are Jason's calls, not Claude's.
