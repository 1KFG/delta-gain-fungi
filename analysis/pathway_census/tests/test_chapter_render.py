"""Render test for chapters/A_carotenoid.qmd (Task 11).

The chapter is a draft scaffold. With no results files it must render (exit 0), show a
"NOT YET RUN" placeholder for every results-driven section, and contain no results table.

How "no fake tables" is checked (kept simple on purpose): the only tables allowed with empty
results are the two definition tables (pathway steps, reference provenance). They come before the
`sec-phylum` heading. Everything from the `sec-phylum` heading onward is results-driven, so that
part of the HTML must contain no `<table` element and no `<img`. The number of placeholders must
equal the number of results files the chapter reads (7).

A second test renders against synthetic result files in a tmp dir (numbers invented for the test
only) to check that the loaders read the `# params:` first line correctly and that placeholders go
away when files exist. Those numbers are never written to the repository.

Both tests skip, with a reason, if quarto cannot be loaded.
"""
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
QMD = ROOT / "chapters" / "A_carotenoid.qmd"
QUARTO_MODULE = "quarto/1.7.31"
N_RESULT_FILES = 7  # phylum summary, 3 paired TSVs, paired JSON, absence summary, triage


def _quarto_cmd():
    if shutil.which("quarto"):
        return "quarto"
    r = subprocess.run(["bash", "-lc", f"module load {QUARTO_MODULE} >/dev/null 2>&1 && command -v quarto"],
                       capture_output=True, text=True)
    return f"module load {QUARTO_MODULE} >/dev/null 2>&1 && quarto" if r.returncode == 0 else None


QUARTO = _quarto_cmd()
pytestmark = pytest.mark.skipif(QUARTO is None, reason=f"quarto not available (PATH or module {QUARTO_MODULE})")


def _render(tmp_path, results_dir):
    out = tmp_path / "out"
    cmd = f"{QUARTO} render {QMD} --output-dir {out}"
    r = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True, cwd=ROOT,
                       env={**os.environ, "PC_RESULTS_DIR": str(results_dir)})
    return r, out / "A_carotenoid.html"


def test_renders_with_no_results_and_shows_placeholders(tmp_path):
    empty = tmp_path / "no_results"
    empty.mkdir()
    r, html = _render(tmp_path, empty)
    assert r.returncode == 0, r.stderr[-2000:]
    text = html.read_text()
    assert "NOT YET RUN" in text
    assert text.count('class="not-yet-run"') == N_RESULT_FILES
    results_part = text[text.index('id="sec-phylum"'):]
    assert "<table" not in results_part
    assert "<img" not in results_part
    # Definition tables still render (pathway steps and provenance).
    assert text.count("<table") == 2


def test_renders_with_synthetic_results_and_drops_placeholders(tmp_path):
    res = tmp_path / "res"
    res.mkdir()
    (res / "carotenoid_phylum_summary.tsv").write_text(
        "# params: FIXTUREPARAMS\nphylum\tstep_id\tn_genomes\tn_with_step\tmedian_copy\tmax_copy\tfrac_with_step\tn_unknown\n"
        "PhylA\tggpp_synthase\t4\t2\t1\t2\t0.5\t0\nPhylA\tcarx\t4\t1\t0\t1\t0.25\t0\n"
        "PhylB\tggpp_synthase\t3\t3\t1\t1\t1.0\t0\nPhylB\tcarx\t3\t0\t0\t0\t0.0\t0\n")
    (res / "carotenoid_triage.tsv").write_text(
        "# params: test\n# PROVISIONAL x\nprotein_id\tstep_id\tbin\treasons\nP1\tcarx\ttier2\tdivergent\n")
    r, html = _render(tmp_path, res)
    assert r.returncode == 0, r.stderr[-2000:]
    text = html.read_text()
    assert text.count('class="not-yet-run"') == N_RESULT_FILES - 2
    assert "PhylA" in text and "P1" in text
    assert "FIXTUREPARAMS" not in text  # comment lines were skipped on reading
    assert re.search(r'id="phylum-figure"', text)


def _sec(text, name):
    """HTML of one level-2 section (from its <section id=...> to the next section)."""
    start = text.index(f'<section id="sec-{name}"')
    nxt = text.find('<section id="sec-', start + 10)
    return text[start:nxt if nxt > 0 else len(text)]


def _text_for(tmp_path, files):
    """Render with the given {name: text} results files; return the HTML text."""
    res = tmp_path / "res"
    res.mkdir()
    for name, body in files.items():
        (res / name).write_text(body)
    r, html = _render(tmp_path, res)
    assert r.returncode == 0, r.stderr[-2000:]
    return html.read_text()


def test_paired_loader_reads_params_line_and_columns(tmp_path):
    # Fixture mimics bin/paired_compare.py: a "# params:" first line, then SUMMARY_COLUMNS.
    # Values are invented for this test only.
    text = _text_for(tmp_path, {"carotenoid_paired_summary.tsv":
        "# params: min_pairs=20 n_boot=2000 seed=1 busco_bins=0,80,95,100\n"
        "step_id\tannotation_provider\tn_pairs\tmean_delta\tci_low\tci_high\tfrac_bfd_more\tfrac_bfd_less\n"
        "carx\tProvX\t30\t0.5\t0.1\t0.9\t0.6\t0.2\n"})
    assert "ProvX" in text and "ci_low" in text and "frac_bfd_less" in text
    assert "min_pairs=20" not in _sec(text, "paired").split("<table")[1]
    assert text.count('class="not-yet-run"') == N_RESULT_FILES - 1


def test_absence_loader_skips_hash_lines_and_counts_outcomes(tmp_path):
    text = _text_for(tmp_path, {"carotenoid_absence_summary.tsv":
        "# min_bitscore=50 min_query_cover=0.7\nlocustag\tstep_id\toutcome\treason\n"
        "AAAA0001\tcarx\tabsent\t\nAAAA0002\tcarx\tunconfirmed\tno_bfd_fasta\n"})
    sec = _sec(text, "absence")
    assert "absent" in sec and "unconfirmed" in sec
    assert "min_bitscore=50" not in sec
    assert "NOT YET RUN" not in sec


def test_provisional_only_file_is_not_shown_as_result(tmp_path):
    text = _text_for(tmp_path, {"carotenoid_absence_summary.PROVISIONAL.tsv":
        "# PROVISIONAL: pathway 'carotenoid' verified: false\n"
        "locustag\tstep_id\toutcome\treason\nPROVROW\tcarx\tabsent\t\n"})
    sec = _sec(text, "absence")
    assert "only a PROVISIONAL file" in sec
    assert "PROVROW" not in text


@pytest.mark.skipif(shutil.which("zstd") is None, reason="zstd not on PATH")
def test_triage_loader_accepts_zst(tmp_path):
    res = tmp_path / "res"
    res.mkdir()
    plain = res / "carotenoid_triage.tsv"
    plain.write_text("# params: x\nprotein_id\tstep_id\tbin\treasons\nZSTROW\tcarx\ttier3\t\n")
    subprocess.run(["zstd", "-q", "--rm", str(plain), "-o", str(plain) + ".zst"], check=True)
    r, html = _render(tmp_path, res)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "ZSTROW" in html.read_text()


def test_readme_merge_recipe_uses_chapter_absence_name():
    assert "carotenoid_absence_summary.tsv" in (ROOT / "README.md").read_text()
    assert "carotenoid_absence_summary.tsv" in QMD.read_text()
