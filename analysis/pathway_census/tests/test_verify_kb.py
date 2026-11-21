"""Tests for bin/verify_kb.py. No network: a fake fetch serves stored fixtures (tests/data/kb)."""
import importlib.util
import json
import re
from pathlib import Path

import pandas as pd
import pytest
import yaml

from pathway_census.kb import COLUMNS
from pathway_census.kbsources import KbFetchError

ROOT = Path(__file__).resolve().parents[1]
KB = Path(__file__).parent / "data" / "kb"
spec = importlib.util.spec_from_file_location("verify_kb", ROOT / "bin" / "verify_kb.py")
vk = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vk)

ABSTRACT = ("A gene has been cloned from Xanthophyllomyces dendrorhous. "
            "The Ast protein is a beta-carotene 3-hydroxylase/4-ketolase.")
QUOTE = "The Ast protein is a beta-carotene 3-hydroxylase/4-ketolase."


def pubmed_xml(title="A title", abstract=ABSTRACT):
    ab = f"<Abstract><AbstractText>{abstract}</AbstractText></Abstract>" if abstract else ""
    return (f"<PubmedArticleSet><PubmedArticle><MedlineCitation><Article>"
            f"<Journal><JournalIssue><PubDate><Year>2006</Year></PubDate></JournalIssue></Journal>"
            f"<ArticleTitle>{title}</ArticleTitle>{ab}</Article></MedlineCitation></PubmedArticle>"
            f"</PubmedArticleSet>")


class Fake:
    def __init__(self, pubmed=None, errors=None):
        self.pubmed = pubmed or {}      # pmid -> xml
        self.errors = errors or {}      # url substring -> exception
        self.urls = []

    def __call__(self, url):
        self.urls.append(url)
        for key, exc in self.errors.items():
            if key in url:
                raise exc
        m = re.search(r"uniprotkb/(\w+)\.json", url)
        if m:
            f = KB / f"uniprot_{m.group(1)}.json"
            if f.exists():
                return f.read_bytes()
            raise KbFetchError(f"request failed for {url}: HTTP Error 404: Not Found")
        if "db=pubmed" in url:
            pmid = re.search(r"id=(\d+)", url).group(1)
            return self.pubmed.get(pmid, "<PubmedArticleSet></PubmedArticleSet>")
        if "db=protein" in url:
            return (KB / "ncbi_protein_WEV88679.gp").read_bytes()
        if "compound/40968" in url:
            return (KB / "chebi_40968.json").read_bytes()
        if "pubchem" in url:
            return (KB / "pubchem_MQZIGYBFDRPAKN.json").read_bytes()
        raise KbFetchError(f"request failed for {url}: HTTP Error 404: Not Found")


def seq_of(acc):
    return json.loads((KB / f"uniprot_{acc}.json").read_text())["sequence"]["value"]


def gp_seq():
    t = (KB / "ncbi_protein_WEV88679.gp").read_text()
    return re.sub(r"[^a-z]", "", t.split("ORIGIN")[1].split("//")[0]).upper()


def world(tmp_path, ref=None, link=None, decoys=None, fasta=None, compound=None):
    """Write a one-step pathway 'tp' (yaml, refs, kb tables) and return (yaml, kb_dir, cache)."""
    ref = {"accession": "Q3HR17", "species": "Phaffia rhodozyma", "gene_name": "Ast",
           "pmid": "16416328", "verified": True, **(ref or {})}
    step = {"id": "s1", "name": "step one", "ec": None, "refs": [ref], "copy_expect": "single",
            "role": "catalytic", "reactions": ["R1"], "domains_any": []}
    if decoys:
        step["decoys"] = decoys
    (tmp_path / "pathways").mkdir()
    ypath = tmp_path / "pathways" / "tp.yaml"
    ypath.write_text(yaml.safe_dump({"id": "tp", "verified": False, "steps": [step]}))
    refs = tmp_path / "refs" / "tp"
    (refs / "decoys").mkdir(parents=True)
    acc = ref["accession"]
    seq = fasta if fasta is not None else (gp_seq() if acc.startswith("WEV") else seq_of(acc))
    (refs / "s1.faa").write_text(f">{acc} test\n{seq}\n")
    for d in decoys or []:
        (refs / "decoys" / "s1.faa").write_text(f">{d['accession']} x\n{seq_of(d['accession'])}\n")
    kbd = tmp_path / "kb" / "tp"
    kbd.mkdir(parents=True)
    link = {"step_id": "s1", "reaction_id": "R1", "organism": ref["species"],
            "gene_name": ref["gene_name"], "protein_accession": acc, "tier": "T1a",
            "assay": "expression", "pmid": ref["pmid"], "quote": QUOTE, **(link or {})}
    comp = {"chebi_id": "CHEBI:40968", "inchikey": "MQZIGYBFDRPAKN-UWFIBFSHSA-N",
            "inchikey_block1": "MQZIGYBFDRPAKN", "smiles": "C", "formula": "C40H52O4",
            "neutral_exact_mass": "596.3865", "kegg_id": "", "np_class": "", **(compound or {})}
    rows = {"step_reaction": [link], "compound": [comp], "reaction": [], "reaction_participant": [],
            "step": []}
    for name, cols in COLUMNS.items():
        pd.DataFrame(rows[name], columns=cols).to_csv(kbd / f"{name}.tsv", sep="\t", index=False)
    return ypath, kbd, tmp_path / "cache"


def run(tmp_path, fake=None, **kw):
    y, k, c = world(tmp_path, **kw)
    fake = fake or Fake({"16416328": pubmed_xml()})
    return vk.verify_pathway(y, k, fetch=fake, cache_dir=c)


def row(df, check, subject_part=""):
    r = df[(df.check == check) & df.subject.str.contains(subject_part, regex=False)]
    assert len(r) >= 1, (check, df[["check", "subject", "status"]])
    return r.iloc[0]


def test_columns_and_matching_reference_passes(tmp_path):
    df = run(tmp_path)
    assert list(df.columns) == ["check", "subject", "expected", "found", "status"]
    assert set(df.status) <= {"pass", "fail", "unchecked"}
    for c in ("ref_exists", "ref_organism", "ref_gene", "ref_sequence", "ref_pmid_exists"):
        assert row(df, c).status == "pass", df[df.check == c]
    assert row(df, "ref_gene").found.startswith("gene name")
    assert row(df, "link_quote_in_abstract").status == "pass"
    assert row(df, "link_quote_one_sentence").status == "pass"
    assert row(df, "link_quote_names_protein").status == "pass"
    assert row(df, "chebi_formula").status == "pass"
    assert row(df, "chebi_inchikey").status == "pass"
    assert row(df, "pubchem_block1").status == "pass"
    assert (df.status == "fail").sum() == 0, df[df.status == "fail"]


def test_swapped_organism_fails(tmp_path):
    df = run(tmp_path, ref={"species": "Neurospora crassa"})
    assert row(df, "ref_organism").status == "fail"
    assert "Phaffia rhodozyma" in row(df, "ref_organism").found


def test_missing_pmid_fails(tmp_path):
    df = run(tmp_path, ref={"pmid": "99999999"}, link={"pmid": "99999999"})
    assert row(df, "ref_pmid_exists").status == "fail"


def test_quote_not_in_abstract_fails(tmp_path):
    df = run(tmp_path, link={"quote": "This sentence is not in the abstract."})
    assert row(df, "link_quote_in_abstract").status == "fail"


def test_quote_whitespace_normalised(tmp_path):
    df = run(tmp_path, link={"quote": "The Ast  protein is a beta-carotene\n3-hydroxylase/4-ketolase."})
    assert row(df, "link_quote_in_abstract").status == "pass"


def test_unavailable_abstract_is_unchecked(tmp_path):
    fake = Fake({"16416328": pubmed_xml(abstract=None)})
    df = run(tmp_path, fake=fake)
    r = row(df, "link_quote_in_abstract")
    assert r.status == "unchecked" and "abstract" in r.found
    assert row(df, "ref_pmid_exists").status == "pass"


def test_pubmed_network_error_is_unchecked(tmp_path):
    fake = Fake({}, errors={"db=pubmed": KbFetchError("request failed for u: timed out")})
    df = run(tmp_path, fake=fake)
    assert row(df, "ref_pmid_exists").status == "unchecked"
    assert row(df, "link_quote_in_abstract").status == "unchecked"


def test_ordered_locus_name_reported(tmp_path):
    df = run(tmp_path, ref={"accession": "B2J6X9", "species": "Nostoc punctiforme (strain ATCC 29133 / PCC 73102)",
                            "gene_name": "Npun_R5600"}, link={"tier": "T2", "quote": QUOTE})
    r = row(df, "ref_gene")
    assert r.status == "pass" and "ordered locus name" in r.found


def test_orf_name_reported(tmp_path):
    df = run(tmp_path, ref={"accession": "P21334", "species": "Neurospora crassa (strain ATCC 24698 / 74-OR23-1A)",
                            "gene_name": "NCU00552"}, link={"tier": "T2"})
    r = row(df, "ref_gene")
    assert r.status == "pass" and "ORF name" in r.found


def test_gene_name_not_in_record_fails(tmp_path):
    df = run(tmp_path, ref={"gene_name": "nope1"})
    assert row(df, "ref_gene").status == "fail"


def test_genbank_reference(tmp_path):
    title = "Unraveling the Molecular Basis of Mycosporine Biosynthesis in Fungi"
    fake = Fake({"36983003": pubmed_xml(title=title)})
    df = run(tmp_path, fake=fake,
             ref={"accession": "WEV88679.1", "species": "Phaffia rhodozyma (strain UCD 67-385)",
                  "gene_name": "DDGS", "pmid": "36983003"},
             link={"tier": "T2", "quote": ""})
    assert row(df, "ref_exists").status == "pass"
    assert row(df, "ref_organism").status == "pass"
    assert row(df, "ref_gene").status == "pass"
    assert row(df, "ref_sequence").status == "pass"
    assert "477" in row(df, "ref_sequence").found
    r = row(df, "ref_record_cites_pmid")
    assert r.status == "pass" and "title" in r.found


def test_genbank_not_citing_pmid_fails(tmp_path):
    fake = Fake({"36983003": pubmed_xml(title="Another paper entirely")})
    df = run(tmp_path, fake=fake,
             ref={"accession": "WEV88679.1", "species": "Phaffia rhodozyma", "gene_name": "DDGS",
                  "pmid": "36983003"}, link={"tier": "T2", "quote": ""})
    assert row(df, "ref_record_cites_pmid").status == "fail"


def test_genbank_unfetchable_is_unchecked(tmp_path):
    fake = Fake({"36983003": pubmed_xml()}, errors={"db=protein": KbFetchError("request failed for u: timed out")})
    df = run(tmp_path, fake=fake,
             ref={"accession": "WEV88679.1", "species": "Phaffia rhodozyma", "gene_name": "DDGS",
                  "pmid": "36983003"}, link={"tier": "T2", "quote": ""})
    r = row(df, "ref_exists")
    assert r.status == "unchecked" and "timed out" in r.found


def test_sequence_mismatch_fails(tmp_path):
    s = seq_of("Q3HR17")
    df = run(tmp_path, fasta=s[:-1] + ("A" if s[-1] != "A" else "G"))
    assert row(df, "ref_sequence").status == "fail"


def test_quote_must_be_one_sentence(tmp_path):
    two = "A gene has been cloned from Xanthophyllomyces dendrorhous. The Ast protein is a beta-carotene 3-hydroxylase/4-ketolase."
    df = run(tmp_path, link={"quote": two})
    assert row(df, "link_quote_in_abstract").status == "pass"
    assert row(df, "link_quote_one_sentence").status == "fail"


def test_one_sentence_ignores_abbreviations(tmp_path):
    assert vk.is_one_sentence("The enzyme of N. crassa and M. circinelloides f. lusitanicus was 3.5 fold higher (Fig. 2).")
    assert not vk.is_one_sentence("It was cloned. It was expressed.")


def test_quote_not_naming_protein_fails_for_t1a_only(tmp_path):
    q = "A gene has been cloned from Xanthophyllomyces dendrorhous."
    df = run(tmp_path, link={"quote": q})
    r = row(df, "link_quote_names_protein")
    assert r.status == "fail" and "quote does not name protein" in r.found
    df = run(tmp_path / "b", link={"quote": q, "tier": "T1b"}) if (tmp_path / "b").mkdir() is None else None
    r = row(df, "link_quote_names_protein")
    assert r.status == "unchecked" and "quote does not name protein" in r.found


def test_decoy_missing_fails(tmp_path):
    y, k, c = world(tmp_path, decoys=[{"accession": "P21334", "species": "Neurospora crassa", "gene_name": "al-1"}])
    # rewrite the decoy to an accession that does not exist
    d = yaml.safe_load(y.read_text())
    d["steps"][0]["decoys"][0]["accession"] = "Q00000"
    y.write_text(yaml.safe_dump(d))
    df = vk.verify_pathway(y, k, fetch=Fake({"16416328": pubmed_xml()}), cache_dir=c)
    assert row(df, "decoy_exists", "Q00000").status == "fail"


def test_decoy_present_passes(tmp_path):
    df = run(tmp_path, decoys=[{"accession": "P21334", "species": "Neurospora crassa", "gene_name": "al-1"}])
    assert row(df, "decoy_exists", "P21334").status == "pass"
    assert row(df, "decoy_organism", "P21334").status == "pass"
    assert row(df, "decoy_sequence", "P21334").status == "pass"


def test_chebi_formula_mismatch_fails(tmp_path):
    df = run(tmp_path, compound={"formula": "C40H52O5"})
    assert row(df, "chebi_formula").status == "fail"


def test_pubchem_no_compound_is_unchecked(tmp_path):
    fake = Fake({"16416328": pubmed_xml()}, errors={"pubchem": KbFetchError("request failed for u: HTTP Error 404: Not Found")})
    df = run(tmp_path, fake=fake)
    r = row(df, "pubchem_block1")
    assert r.status == "unchecked" and "no compound" in r.found


def test_link_organism_must_match_record(tmp_path):
    df = run(tmp_path, link={"organism": "Neurospora crassa"})
    assert row(df, "link_organism").status == "fail"


def test_link_gene_mismatch_fails(tmp_path):
    df = run(tmp_path, link={"gene_name": "zzz"})
    assert row(df, "link_gene").status == "fail"


def test_decoy_gene_mismatch_fails(tmp_path):
    df = run(tmp_path, decoys=[{"accession": "P21334", "species": "Neurospora crassa", "gene_name": "zzz"}])
    assert row(df, "decoy_gene", "P21334").status == "fail"


def test_chebi_inchikey_mismatch_fails(tmp_path):
    df = run(tmp_path, compound={"inchikey": "AAAAAAAAAAAAAA-UWFIBFSHSA-N", "inchikey_block1": "AAAAAAAAAAAAAA"})
    assert row(df, "chebi_inchikey").status == "fail"


def test_pubchem_block_mismatch_fails(tmp_path):
    df = run(tmp_path, compound={"inchikey": "AAAAAAAAAAAAAA-UWFIBFSHSA-N", "inchikey_block1": "AAAAAAAAAAAAAA"})
    assert row(df, "pubchem_block1").status == "fail"


def test_pubchem_formula_mismatch_fails(tmp_path):
    df = run(tmp_path, compound={"formula": "C40H52O5"})
    assert row(df, "pubchem_formula").status == "fail"
    assert row(df, "pubchem_block1").status == "pass"


def test_pubchem_formula_passes(tmp_path):
    assert row(run(tmp_path), "pubchem_formula").status == "pass"


def test_t1a_without_quote_fails(tmp_path):
    df = run(tmp_path, link={"quote": ""})
    r = row(df, "link_evidence_present")
    assert r.status == "fail" and "quote" in r.found


def test_t1a_without_pmid_fails(tmp_path):
    df = run(tmp_path, link={"pmid": ""})
    assert row(df, "link_evidence_present").status == "fail"


def test_t3_without_pmid_and_quote_passes_silently(tmp_path):
    df = run(tmp_path, link={"tier": "T3", "pmid": "", "quote": ""})
    assert (df.check == "link_evidence_present").sum() == 0
    assert not (df[df.check.str.startswith("link_")].status == "fail").any()


def test_region_slice_passes(tmp_path):
    s = seq_of("Q3HR17")
    df = run(tmp_path, ref={"region": [11, 60]}, fasta=s[10:60])
    assert row(df, "ref_sequence").status == "pass"


def test_region_wrong_position_fails(tmp_path):
    s = seq_of("Q3HR17")
    df = run(tmp_path, ref={"region": [1, 50]}, fasta=s[10:60])
    r = row(df, "ref_sequence")
    assert r.status == "fail" and "substring at 11" in r.found


def test_substring_without_region_fails(tmp_path):
    s = seq_of("Q3HR17")
    assert row(run(tmp_path, fasta=s[10:60]), "ref_sequence").status == "fail"


def test_short_names_match_case_sensitively():
    v = vk.Verifier()
    rec = {"names": [], "protein_names": []}
    assert v.names_quote(rec, "T", "the t was low") is None
    assert v.names_quote(rec, "T", "the T gene was low") == "T"
    assert v.names_quote(rec, "tdh", "the TDH gene") == "tdh"
    assert v.names_quote(rec, "tdh", "the tdhx gene") is None
    assert v.names_quote(rec, "al-3", "the al3 mutant") == "al-3"


def test_unusable_response_is_unchecked(tmp_path):
    fake = Fake({}, errors={"db=pubmed": KbFetchError("pubmed 1: unusable response from u (ValueError('x'))")})
    assert row(run(tmp_path, fake=fake), "ref_pmid_exists").status == "unchecked"


def test_genbank_version_mismatch_fails(tmp_path):
    fake = Fake({"36983003": pubmed_xml()})
    df = run(tmp_path, fake=fake, ref={"accession": "WEV88679.2", "species": "Phaffia rhodozyma",
                                       "gene_name": "DDGS", "pmid": "36983003"}, link={"tier": "T3", "quote": ""})
    assert row(df, "ref_exists").status == "fail"


def test_report_counts_exact():
    df = pd.DataFrame([
        ["ref_gene", "a", "x", "ORF name 'n'", "pass"], ["ref_gene", "b", "x", "gene name 'g'", "pass"],
        ["ref_organism", "c", "x", "y", "fail"], ["chebi_formula", "d", "x", "reason one", "unchecked"]], columns=vk.COLS)
    md = vk.render_report({"tp": df})
    assert "| tp | 4 | 2 | 1 | 1 |" in md
    assert "| all | 4 | 2 | 1 | 1 |" in md
    assert "## Failures (1)" in md and "## Unchecked (1)" in md
    assert "| tp | 1 | 0 | 1 | 0 |" in md          # gene-name passes by kind
    assert "| ref_organism | 0 | 1 | 0 |" in md
    assert "| tp | ref_organism | c | x | y |" in md and "| tp | chebi_formula | d | reason one |" in md
