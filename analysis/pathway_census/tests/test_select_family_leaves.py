import importlib.util
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location("select_family_leaves", Path(__file__).resolve().parents[1] / "bin" / "select_family_leaves.py")
sfl = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(sfl)


def srow(prot, score, cov=0.9, ok=True):
    return {"hmm": "h", "protein": prot, "entry": prot, "best_score": str(score), "best_cov": str(cov), "union_cov": str(cov),
            "pass_best": str(ok), "pass_union": str(ok)}


def mrow(acc, org, fungal=True, ec="", length=100):
    return {"Entry": acc, "Organism": org, "Taxonomic lineage": "Eukaryota, Fungi (kingdom)" if fungal else "Eukaryota, Metazoa",
            "Gene Names": "g", "Protein names": "p", "Length": str(length), "EC number": ec}


def test_accession_of():
    assert sfl.accession_of("sp|P12345|X_HUMAN") == "P12345"
    assert sfl.accession_of("P12345") == "P12345"


def test_organism_short():
    assert sfl.organism_short("Aspergillus fumigatus (strain Af293) (Neosartorya)") == "Aspergillus fumigatus"


def test_eligible_excludes_reference_and_same_ec_and_failed_screen():
    rows = [srow("sp|A1|N", 300), srow("sp|A2|N", 200), srow("sp|A3|N", 150), srow("sp|A4|N", 100, ok=False)]
    meta = {"A1": mrow("A1", "x"), "A2": mrow("A2", "y", ec="1.1.1.1"), "A3": mrow("A3", "z")}
    elig, excl = sfl.eligible_rows(rows, meta, "best", {"A1"}, ["1.1.1.1"])
    assert [d["accession"] for d in elig] == ["A3"]
    reasons = dict(excl)
    assert "reference" in reasons["A1"] and "same EC" in reasons["A2"]
    assert "A4" not in reasons   # failed screen is not an exclusion, it was never eligible


def test_eligible_sorted_by_score_and_missing_meta():
    rows = [srow("sp|B1|N", 50), srow("sp|B2|N", 90), srow("sp|B3|N", 70)]
    meta = {"B1": mrow("B1", "a"), "B2": mrow("B2", "b")}
    elig, excl = sfl.eligible_rows(rows, meta, "best", set(), [])
    assert [d["accession"] for d in elig] == ["B2", "B1"]
    assert excl == [("B3", "no UniProt table row")]


def _d(acc, score, org, fungal):
    return dict(accession=acc, best_score=score, organism=org, fungal=fungal)


def test_choose_one_per_organism_and_caps():
    elig = [_d("F1", 300, "Aspergillus niger (s)", True), _d("F2", 290, "Aspergillus niger (t)", True),
            _d("F3", 280, "Fusarium oxysporum", True), _d("N1", 100, "Homo sapiens", False),
            _d("N2", 90, "Mus musculus", False), _d("N3", 80, "Danio rerio", False), _d("N4", 70, "Bos taurus", False)]
    named, outg = sfl.choose(elig, n_named=2, n_outgroup=2)
    assert [d["accession"] for d in named] == ["F1", "F3"]
    assert [d["accession"] for d in outg] == ["N1", "N4"]   # first and last rank of the non-fungal pool


def test_choose_fungal_leftover_goes_to_outgroups_when_no_other():
    elig = [_d("F1", 300, "A a", True), _d("F2", 290, "B b", True), _d("F3", 280, "C c", True)]
    named, outg = sfl.choose(elig, n_named=1, n_outgroup=5)
    assert [d["accession"] for d in named] == ["F1"]
    assert [d["accession"] for d in outg] == ["F2", "F3"]


def test_choose_zero_outgroups():
    named, outg = sfl.choose([_d("F1", 1, "A a", True)], 3, 0)
    assert len(named) == 1 and outg == []


def test_main_writes_files(tmp_path):
    fa = tmp_path / "f.faa"
    fa.write_text(">sp|A1|N_X Protein one OS=X\nMKT\n>sp|A2|N_Y Protein two OS=Y\nMKTA\n")
    meta = tmp_path / "m.tsv"
    meta.write_text("Entry\tOrganism\tTaxonomic lineage\tGene Names\tProtein names\tLength\tEC number\n"
                    "A1\tX x\tFungi\tg1\tone\t3\t\nA2\tY y\tMetazoa\tg2\ttwo\t4\t\n")
    scr = tmp_path / "s.tsv"
    scr.write_text("hmm\tprotein\tentry\tbest_score\tbest_cov\tn_domains_ge_min\tunion_cov\tpass_best\tpass_union\n"
                   "h\tsp|A1|N_X\te\t50\t0.9\t1\t0.9\tTrue\tTrue\nh\tsp|A2|N_Y\te\t40\t0.9\t1\t0.9\tTrue\tTrue\n")
    out = {k: tmp_path / k for k in ("n.faa", "o.faa", "p.tsv", "e.tsv")}
    rc = sfl.main(["--screen", str(scr), "--meta", str(meta), "--fasta", str(fa), "--hmm", "h", "--out-named", str(out["n.faa"]),
                   "--out-outgroups", str(out["o.faa"]), "--out-provenance", str(out["p.tsv"]), "--out-exclusions", str(out["e.tsv"])])
    assert rc == 0
    assert out["n.faa"].read_text().startswith(">A1 ")
    assert out["o.faa"].read_text().startswith(">A2 ")
    assert "named_member" in out["p.tsv"].read_text() and "tree_outgroup" in out["p.tsv"].read_text()


def test_eligible_excludes_scores_at_or_above_floor():
    rows = [srow("sp|C1|N", 400), srow("sp|C2|N", 300), srow("sp|C3|N", 299.9)]
    meta = {k: mrow(k, k) for k in ("C1", "C2", "C3")}
    elig, excl = sfl.eligible_rows(rows, meta, "best", set(), [], max_score=300)
    assert [d["accession"] for d in elig] == ["C3"]
    assert {a for a, _ in excl} == {"C1", "C2"}
    assert all("floor" in r for _, r in excl)


def test_eligible_excludes_reference_gene_names_case_insensitive():
    rows = [srow("sp|G1|N", 100), srow("sp|G2|N", 90), srow("sp|G3|N", 80)]
    meta = {"G1": dict(mrow("G1", "a"), **{"Gene Names": "Arp1 MAN_00844"}), "G2": dict(mrow("G2", "b"), **{"Gene Names": "other"}),
            "G3": dict(mrow("G3", "c"), **{"Gene Names": ""})}
    elig, excl = sfl.eligible_rows(rows, meta, "best", set(), [], exclude_genes=["arp1", "scd1"])
    assert [d["accession"] for d in elig] == ["G2", "G3"]
    assert excl[0][0] == "G1" and "gene name" in excl[0][1]
