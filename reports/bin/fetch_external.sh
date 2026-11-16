#!/bin/bash
# Fetch the delta-gain funnel findings for the report page at reports/delta_gain/.
#
#   reports/bin/fetch_external.sh
#
# The funnel narrative is not in this repository: it belongs to the pipeline
# repository 1KFG/delta-gain, which produced it. Copying it here would create a
# second copy to keep in step with the first, so instead the page reproduces it
# from that repository at a PINNED COMMIT, and the pin is part of what the page
# claims. The commit is named on the page itself.
#
# Both the commit and the content are pinned: a raw.githubusercontent URL can
# in principle be re-served, and a silent change in the text behind a page that
# says "reproduced verbatim" is exactly the failure worth catching. If a
# checksum does not match, nothing is written and the script exits non-zero.
#
# Output is gitignored. The published HTML under docs/reports/ is committed, so
# the live site does not depend on this script or on network access; only
# re-rendering does.
set -euo pipefail

REPO=1KFG/delta-gain
# Findings synthesis (docs/FINDINGS.md) with E4 accumulation figure, 2026-11-08.
SHA=91eb5fe6be888707fcf9c6953ef68a99bd91d414

cd "$(dirname "$0")/.."
DEST=delta_gain
mkdir -p "$DEST/_external" "$DEST/figures"

# published path : source path in the pinned repo : expected sha256
FILES=(
  "$DEST/_external/FINDINGS.md:docs/FINDINGS.md:2ba2b1029fd1fe547a1e90c90b4511d0ce24062f427b5c36d7c433e973989f05"
  "$DEST/figures/E4_accumulation.svg:docs/figures/E4_accumulation.svg:ab887c0d3262fb45a0f6a90849789b1f03e057215387b888eb0e2b898ad05652"
)

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

for entry in "${FILES[@]}"; do
    out=${entry%%:*}
    rest=${entry#*:}
    src=${rest%%:*}
    want=${rest#*:}

    url="https://raw.githubusercontent.com/$REPO/$SHA/$src"
    curl -fsSL -o "$tmp/f" "$url"

    got=$(shasum -a 256 "$tmp/f" | cut -d' ' -f1)
    if [[ "$got" != "$want" ]]; then
        echo "CHECKSUM MISMATCH for $src" >&2
        echo "  url      $url" >&2
        echo "  expected $want" >&2
        echo "  got      $got" >&2
        echo "The pinned content changed under the pinned commit. Do not render." >&2
        exit 1
    fi

    mv "$tmp/f" "$out"
    echo "ok  $out  ($src @ ${SHA:0:7})"
done

echo "fetched from $REPO @ $SHA"
