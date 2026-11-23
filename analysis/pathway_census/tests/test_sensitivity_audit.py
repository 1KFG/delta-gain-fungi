"""Tests for bin/sensitivity_audit.py (network access is replaced by an injected fetch function)."""
import importlib.util
from pathlib import Path

import pytest

from pathway_census import hits

SPEC = importlib.util.spec_from_file_location("sensitivity_audit", Path(__file__).parents[1] / "bin" / "sensitivity_audit.py")
sa = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sa)
DATA = Path(__file__).parent / "data"
needs_hmmer = pytest.mark.skipif(not all(hits.tool_available(t) for t in ["hmmbuild", "hmmsearch", "mafft"]), reason="hmmer/mafft not loadable")


def write_kb(d):
    d.mkdir(parents=True, exist_ok=True)
    (d / "reaction.tsv").write_text("reaction_id\tec\tkegg_rn\tdirection\tcurated\tnote\trhea_id\n"
                                    "R1\t2.5.1.29\t\tf\tf\t\t\nR2\t1.14.99.63\t\tf\tf\t\t\nR3\t\t\tf\tf\t\t\n")
    (d / "step_reaction.tsv").write_text("step_id\treaction_id\torganism\tgene_name\tprotein_accession\ttier\tassay\tpmid\tquote\n"
                                         "s1\tR1\to\tg\tP\tT3\t\t\t\ns2\tR2\to\tg\tP\tT3\t\t\t\ns3\tR2\to\tg\tP\tT3\t\t\t\ns4\tR3\to\tg\tP\tT3\t\t\t\n")


def test_ec_by_step_only_keeps_reactions_owned_by_one_step(tmp_path):
    write_kb(tmp_path / "kb")
    assert sa.ec_by_step(tmp_path / "kb") == {"s1": ["2.5.1.29"]}   # R2 is shared by s2 and s3, R3 has no EC


def test_read_fasta_uses_uniprot_accession_from_sp_and_tr_headers():
    txt = ">sp|P12345|GGPS_X Geranylgeranyl synthase OS=Some fungus OX=1 GN=g PE=1 SV=1\nMKV\n>tr|Q9ZZZ9|Q9ZZZ9_Y Foo OS=Other OX=2\nMM\n>plain desc\nMA\n"
    assert [n for n, _, _ in sa.read_fasta(txt)] == ["P12345", "Q9ZZZ9", "plain"]


def test_audit_step_recall_at_cutoff_and_rule_and_no_hit():
    entries = [("A", "A n1 OS=Org One OX=1"), ("B", "B n2 OS=Org Two OX=2"), ("C", "C n3 OS=Org Three OX=3")]
    rows, s = sa.audit_step(entries, {"A": 500.0, "B": 300.0}, cutoff=400.0, rule_cutoff=250.0)
    assert s["n_entries"] == 3 and s["n_no_hit"] == 1 and s["n_at_cutoff"] == 1 and s["recall_at_cutoff"] == pytest.approx(1 / 3)
    assert s["n_at_rule"] == 2 and s["median_score"] == 400.0
    assert rows[2]["score"] is None and rows[2]["at_cutoff"] is False and rows[0]["organism"] == "Org One"


def test_audit_step_with_no_entries_gives_none_recall():
    rows, s = sa.audit_step([], {}, 10.0)
    assert rows == [] and s["recall_at_cutoff"] is None and s["n_entries"] == 0


def test_own_accessions_reads_refs_and_decoys(tmp_path):
    (tmp_path / "decoys").mkdir()
    (tmp_path / "s1.faa").write_text(">P1 ref\nMK\n")
    (tmp_path / "decoys" / "s1.faa").write_text(">P2 decoy\nMK\n")
    assert sa.own_accessions(tmp_path, "s1") == {"P1", "P2"}
    assert sa.own_accessions(tmp_path, "none") == set()


@needs_hmmer
def test_main_end_to_end_with_fake_fetch(tmp_path):
    write_kb(tmp_path / "kb")
    refs = tmp_path / "refs"
    (refs / "decoys").mkdir(parents=True)
    tiny = (DATA / "tiny_refs.faa").read_text()
    (refs / "s1.faa").write_text(tiny)
    hmms = tmp_path / "hmms"
    hmms.mkdir()
    hits.build_hmm(refs / "s1.faa", hmms / "s1.PROVISIONAL.hmm")
    ms = tmp_path / "ms.tsv"
    ms.write_text("step\tmin_score\tbasis\tn_refs\trule_min_score\ns1\t20\tx\t1\t30\ns4\t5\tx\t1\t5\n")
    proteome = (DATA / "tiny_proteome.faa").read_text()
    fake = lambda ec: ">sp|HOMOLOGUE|N_X Homologue OS=Fungus A OX=1\n" + "".join(l for l in proteome.split(">homologue\n")[1].split(">")[0].splitlines(True))
    out = tmp_path / "out"
    assert sa.main(["--kb-dir", str(tmp_path / "kb"), "--refs-dir", str(refs), "--hmm-dir", str(hmms), "--min-scores", str(ms), "--out-dir", str(out)], fetch=fake) == 0
    summ = sa.read_tsv(out / "sensitivity_summary.tsv")
    assert summ[0]["step"] == "s1" and summ[0]["n_entries"] == "1" and summ[0]["n_at_cutoff"] == "1"
    assert "s4" in (out / "not_auditable.tsv").read_text()
