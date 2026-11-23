"""Build the mycosporine / mycosporine-like amino acid curation files from checked records.

Writes pathways/mycosporine.yaml, refs/mycosporine/ (FASTA, provenance, decoys,
fungal_lineages.tsv), kb/mycosporine/*.tsv and tests/data/kb/mycosporine_quotes.tsv. It does not
write UNVERIFIED.md. Every curated fact below is checked against a fetched record by an assert.
Run in this directory: /usr/bin/python3.12 build.py
"""
from records import *
from fetch_fulltext import fulltext
from find_compounds import pubchem_name, chebi_by_text
import csv, re, textwrap
from pathway_census.kb import formula_mass, COLUMNS, species_key

ROOT = PKG
PW = "mycosporine"
REFS = ROOT / "refs" / PW; DEC = REFS / "decoys"; KB = ROOT / "kb" / PW
TODAY = "2026-11-21"
FUNGAL_PMID = "36983003"   # Phaffia rhodozyma UCD 67-385: DDGS, OMT, ATPG
MYS_PMID = "21890703"      # Nostoc punctiforme ATCC 29133: NpR5600, NpR5599, NpR5598, NpF5597

def wrap(seq):
    return "\n".join(seq[i:i+60] for i in range(0, len(seq), 60))

def norm(s):
    return " ".join(s.split())

def ab(pmid):
    return pubmed_abstract(pmid)["abstract"]

def quote(pmid, start, end):
    """Cut one sentence out of the fetched abstract. The text is never typed by hand."""
    a = ab(pmid)
    i = a.index(start); j = a.index(end, i) + len(end)
    q = a[i:j]
    assert q.endswith(".") and len(q) <= 400, (pmid, q)
    assert norm(q) in norm(a)
    assert not re.search(r"[.!?] [A-Z]", q[:-1]), (pmid, q)
    return q

def tsv(path, cols, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar=None)
        w.writerow(cols)
        for r in rows:
            assert len(r) == len(cols), r
            assert not any("\t" in str(x) or "\n" in str(x) for x in r), r
            w.writerow(r)

# One sentence per link. A paper can give different sentences to different links.
Q = {
 "del": quote(FUNGAL_PMID, "Homozygous deletion mutants of the entire cluster", "did not produce mycosporines."),
 "het": quote(FUNGAL_PMID, "Heterologous expression of the DDGS and OMT", "MG production, respectively."),
 "mys": quote(MYS_PMID, "Heterologous expression in Escherichia coli of three contiguous", "an oxomycosporine."),
 "cmp": quote("22741921", "Comparative bioinformatics and biochemical analyses", "plausible catalytic mechanisms."),
}

# ---------- the fungal paper and its proteins (each fact is an assert) ----------
a = ab(FUNGAL_PMID)
assert a.startswith("The Phaffia rhodozyma UCD 67-385 genome harbors a 7873 bp cluster containing DDGS, OMT, and ATPG")
assert "However, Δatpg-/- accumulated the intermediate 4-deoxygadusol." in a
assert "non-mycosporine-producing CBS 6938 wild-type strain" in a and "produced MG and mycosporine glutaminol glucoside" in a
summ = raw("pubmed_summary", FUNGAL_PMID, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=" + FUNGAL_PMID)["result"][FUNGAL_PMID]
DOI = [x["value"] for x in summ["articleids"] if x["idtype"] == "doi"][0]
assert DOI == "10.3390/ijms24065930"
for acc, gene, length in (("WEV88679.1", "DDGS", 477), ("WEV88680.1", "OMT", 284), ("WEV88681.1", "ATPG", 500)):
    r = rec(acc)
    assert r["gene_names"] == [gene] and r["length"] == length and r["taxid"] == "264483", acc
    assert r["species"] == "Phaffia rhodozyma (strain UCD 67-385)", acc
    # the GenBank record cites the paper by title, journal and DOI (it has no PUBMED line)
    assert (pubmed_abstract(FUNGAL_PMID)["title"].rstrip("."), "Int J Mol Sci 24 (6), 5930 (2023)", "DOI: " + DOI) in r["refs"], acc
    # UniParc knows the sequence and gives it no UniProtKB accession
    up = raw("uniparc_search", acc.split(".")[0], f"https://rest.uniprot.org/uniparc/search?query={acc.split('.')[0]}&fields=upi,accession,length,organism&format=tsv", as_json=False)
    row = up.strip().splitlines()[1].split("\t")
    assert len(up.strip().splitlines()) == 2 and row[1] == "" and row[2] == str(length), (acc, up)
# full text (Europe PMC, PMC10057719; read on 2026-11-21): accession of the cluster, method, glucosylation, AROM
pmc, ft = fulltext(FUNGAL_PMID)
assert pmc == "PMC10057719"
for frag in ("the MYC cluster (Genbank OQ547787 )",
             "Mycosporines were identified according to their retention time and absorption spectra.",
             "catalyzed by the product of a P. rhodozyma -specific gene that is not part of the MYC cluster and whose function is not present in S. cerevisiae",
             "there is no relationship between the Dhqs activity of the AROM gene product and MGG production",
             "The enzyme ATP-grasp ligase subsequently incorporates glutaminol into this intermediate molecule, leading to the formation of MG",
             "filamentous fungi, such as Ascochyta fabae, that produce mycosporine glutaminol glucoside (MGG)"):
    assert frag in ft, frag
assert "mass spectrometry" not in ft.split("3.3.")[1].split("3.4.")[0]  # the mycosporine analysis section reports HPLC only
# ATPG of strain VOH (UniProt A0A9Y1GGR7): another strain, 10 differences in 500 residues, no cited paper
voh, ucd = rec("A0A9Y1GGR7"), rec("WEV88681.1")
assert voh["gene_names"] == ["ATPG"] and voh["organism"] == "Phaffia rhodozyma" and voh["pmids"] == [] and voh["length"] == 500
assert sum(x != y for x, y in zip(voh["sequence"], ucd["sequence"])) == 10
assert "OP021064" in [x["id"] for x in full("A0A9Y1GGR7")["uniProtKBCrossReferences"] if x["database"] == "EMBL"]
assert "PS50975" in voh["xrefs"] and not {x for x in voh["xrefs"] if x.startswith("PF")}
# the bacterial papers
assert "Nostoc punctiforme ATCC 29133" in ab(MYS_PMID)
assert "Additional expression of gene NpF5597 (mysD) led to the conversion of mycosporine-glycine into iminomycosporines (preferentially shinorine but also others like mycosporine-2-glycine and porphyra-334)." in ab(MYS_PMID)
assert "These enzymes catalyze the conversion of sedoheptulose 7-phosphate to 2-epi-valiolone" in ab("22741921") and "Actinosynnema mirum" in ab("22741921")
x = raw("ncbi_pmc", "3405187", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=3405187&retmode=xml", as_json=False)
t = norm(re.sub(r"<[^>]+>", " ", x))
assert "The catalytic function of Ava_3858 has not previously been characterized, however, its high sequence identity to Npun_5600 suggested that it is a DDGS. This activity was confirmed in the present study" in t
x = raw("ncbi_pmc", "3116657", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=3116657&retmode=xml", as_json=False)
t = norm(re.sub(r"<[^>]+>", " ", x))
assert "treatment of SH 7-P with NpR5600 and NpR5599 in the presence of SAM, NAD + , and Co 2+ produced a single product, 4-deoxygadusol" in t
assert "ATP-grasp homolog Ava_3856 and NRPS-like enzyme Ava_3855" in t
assert "characterized all four biosynthetic enzymes in vitro" in ab("20813918")
assert "Pseudonocardia sp. strain P1" in ab("24907338") and "produced mainly shinorine" in ab("24907338")
assert "A biosynthetic gene cluster was identified in the N. friedmannii genome" in ab("40036870")

# ---------- references: step -> [(accession, gene_name, reference pmid, verified)] ----------
REF = {
 "ddgs": [("WEV88679.1", "DDGS", FUNGAL_PMID, False), ("B2J6X9", "Npun_R5600", "20813918", True), ("Q3M6C3", "Ava_3858", "22741921", True)],
 "omt": [("WEV88680.1", "OMT", FUNGAL_PMID, False), ("B2J6X8", "Npun_R5599", MYS_PMID, False)],
 "atp_grasp_ligase": [("WEV88681.1", "ATPG", FUNGAL_PMID, False), ("A0A9Y1GGR7", "ATPG", FUNGAL_PMID, False), ("B2J6X7", "Npun_R5598", MYS_PMID, False)],
}
REF_NOTE = {
 ("ddgs", "WEV88679.1"): [
   "GenBank protein of the Phaffia rhodozyma UCD 67-385 cluster OQ547787 (the record cites PMID 36983003 by title,",
   "journal and DOI). UniProtKB has no entry for this sequence, so the UniProt rule cannot be met: verified: false.",
   "species is the GenBank organism and strain. USER DECISION 1 in UNVERIFIED.md."],
 ("ddgs", "B2J6X9"): [
   "Cyanobacterium; never fungal evidence (T2 at most). gene_name is the ordered locus name of the UniProt entry;",
   "the papers write NpR5600 and Npun_5600 (mysA). pmid 20813918 is cited by the entry for the catalytic activity."],
 ("ddgs", "Q3M6C3"): [
   "Cyanobacterium; never fungal evidence. UniProt organism: Trichormus variabilis; the papers say Anabaena variabilis",
   "ATCC 29413 (one NCBI taxid, 264691). gene_name is the ordered locus name. Link held at T3."],
 ("omt", "WEV88680.1"): ["GenBank protein of cluster OQ547787; no UniProtKB entry: verified: false (USER DECISION 1 in UNVERIFIED.md)."],
 ("omt", "B2J6X8"): [
   "Cyanobacterium; never fungal evidence. NpR5599 (mysB) of PMID 21890703 (pmid below; the entry does not cite it).",
   "Unreviewed entry: its only citation is PMID 23463784, a sugar-transporter paper that UniProt uses for the genome",
   "sequence; it is not a paper about this step. EC 2.1.1.104 is automatic annotation: verified: false. 17% identical",
   "to the Phaffia OMT (MEASURED): the same enzyme family is not shown. Link held at T3 (UNVERIFIED.md, B3)."],
 ("atp_grasp_ligase", "WEV88681.1"): ["GenBank protein of cluster OQ547787; no UniProtKB entry: verified: false (USER DECISION 1 in UNVERIFIED.md)."],
 ("atp_grasp_ligase", "A0A9Y1GGR7"): [
   "ATPG of another strain (VOH, mRNA OP021064): 490 of 500 residues equal to WEV88681.1. The entry cites no paper",
   "and the strain is not the assayed one: verified: false, link T3. pmid is the paper of the UCD 67-385 gene."],
 ("atp_grasp_ligase", "B2J6X7"): [
   "Cyanobacterium; never fungal evidence. NpR5598 (mysC) of PMID 21890703: it ligates glycine, not glutaminol.",
   "pmid is that paper; the entry does not cite it. Unreviewed entry with no EC and no reaction; its only citation is",
   "PMID 23463784 (a sugar-transporter paper used for the genome sequence, not a paper about this step): verified: false."],
}
REF_EC = {"ddgs": "4.2.3.154"}
REF_RHEA = {"ddgs": {"RHEA:49560"}}
# family ids that every reference of the step carries (UniProt cross-references; for the GenBank
# proteins: hmmscan against Pfam 38.2 with --cut_ga, run on 2026-11-21, see README)
HMMSCAN_GA = {"WEV88679.1": {"PF01761", "PF24621", "PF13685"}, "WEV88680.1": {"PF13578"}, "WEV88681.1": set()}
DOMAINS = {"ddgs": ["PF01761", "PF24621"], "omt": ["PF13578", "PF01596"], "atp_grasp_ligase": []}

def family_ids(acc):
    return HMMSCAN_GA[acc] if is_genbank(acc) else rec(acc)["xrefs"]

REFS.mkdir(parents=True, exist_ok=True); DEC.mkdir(exist_ok=True); KB.mkdir(parents=True, exist_ok=True)
prov = []
for step, refs in REF.items():
    recs = []
    for acc, gene, pmid, ver in refs:
        r = rec(acc)
        assert gene in r["gene_names"], (acc, gene, r["gene_names"])
        assert pubmed_abstract(pmid)["title"]
        if r["source"] == "UniProtKB" and acc not in ("A0A9Y1GGR7", "B2J6X8", "B2J6X7"):
            assert pmid in r["pmids"], (acc, pmid, r["pmids"])
        # verified: true only for a UniProt entry with a gene name that cites the paper and carries the EC and Rhea id of the step
        ok = (r["source"] == "UniProtKB" and pmid in r["pmids"] and step in REF_EC and REF_EC[step] in r["ec"] and REF_RHEA[step] <= r["rhea"])
        assert ver == ok, (step, acc, ok)
        assert set(DOMAINS[step]) & family_ids(acc) or not DOMAINS[step], (step, acc)
        recs.append(f">{r['header']}\n{wrap(r['sequence'])}\n")
        src = r["srcdb"] + ("; PubMed (paper cited by the GenBank record)" if is_genbank(acc) else "; PubMed" + (" (cited by the entry)" if pmid in r["pmids"] else " (not cited by the entry)"))
        if acc in ("B2J6X8", "B2J6X7"):
            src += "; the only citation of the entry, PMID 23463784, is a sugar-transporter paper used for the genome sequence"
        prov.append([step, acc, r["species"], gene, "", pmid, src, TODAY, str(ver).lower()])
    (REFS / f"{step}.faa").write_text("".join(recs), encoding="utf-8")
tsv(REFS / "provenance.tsv", ["step_id", "accession", "species", "gene_name", "region", "pmid", "source_database", "retrieved_on", "verified"], prov)
# ddgs: all three references carry both Pfam families; omt: the fungal and the bacterial reference carry different ones
assert all({"PF01761", "PF24621"} <= family_ids(a) for a, *_ in REF["ddgs"])
assert family_ids("WEV88680.1") == {"PF13578"} and {x for x in rec("B2J6X8")["xrefs"] if x.startswith("PF")} == {"PF01596"}
assert {"PF02655", "PS50975"} <= rec("B2J6X7")["xrefs"]   # MysC: ATP-grasp_3 and the PROSITE ATP-grasp profile
# the bacterial references: ordered locus names only (kbsources.uniprot_entry returns no gene name for them)
for acc in ("B2J6X9", "Q3M6C3", "B2J6X8", "B2J6X7"):
    assert uniprot_entry(acc)["gene_names"] == [] and len(ordered_locus_names(acc)) == 1, acc
assert rec("B2J6X9")["reviewed"] and rec("Q3M6C3")["reviewed"] and "22741921" in rec("B2J6X9")["pmids"]
assert MYS_PMID not in rec("B2J6X9")["pmids"] + rec("B2J6X8")["pmids"] + rec("B2J6X7")["pmids"]
# the only citation of the two unreviewed Nostoc entries is a sugar-transporter paper (genome sequence), not a paper about the step
assert rec("B2J6X8")["pmids"] == ["23463784"] == rec("B2J6X7")["pmids"] and "sugar transporter" in pubmed_abstract("23463784")["title"]
def ca_evidence(acc):
    return [(e["evidenceCode"], e.get("id")) for c in full(acc).get("comments", []) if c["commentType"] == "CATALYTIC ACTIVITY" for e in (c["reaction"].get("evidences") or [])]
# B2J6X9: experimental evidence (ECO:0000269) for RHEA:49560 from two papers; B2J6X8: no catalytic activity record, existence by homology
assert ca_evidence("B2J6X9") == [("ECO:0000269", "20813918"), ("ECO:0000269", "22741921")]
assert ca_evidence("B2J6X8") == [] and full("B2J6X8")["proteinExistence"].startswith("3: Inferred from homology")
# O42898: the catalytic activity has no evidence code; 'Probable' in the protein name
assert ca_evidence("O42898") == [] and "RHEA:17877" in rec("O42898")["rhea"]
assert full("O42898")["proteinDescription"]["recommendedName"]["fullName"]["value"] == "Probable catechol O-methyltransferase 1"
# Ascochyta fabae: a mycosporin glucosyltransferase activity (abstract read; no sequence)
GT_PMID = "24258137"
ag = ab(GT_PMID)
assert pubmed_abstract(GT_PMID)["title"].startswith("Soluble uridine diphospho-D-glucose: mycosporin glucosyltransferase from spores of Ascochyta fabae")
for frag in ("uridine 5'-diphosphate (UDP) glucose: mycosporin-2 glucosyltransferase from spores of Ascochyta fabae", "(the amide form being the best acceptor)", "pH of about 8.5", "divalent cations"):
    assert frag in ag, frag
assert "glutaminol" not in ag   # the abstract does not say that mycosporin-2 is mycosporine-glutaminol
assert rec("B2J6X8")["ec"] == ["2.1.1.104"] and not rec("B2J6X8")["rhea"] and rec("B2J6X7")["ec"] == [] and not rec("B2J6X7")["rhea"]
assert "(PubMed:22741921)" in function_text("Q3M6C3")
# bacterial genes that are neither references nor decoys
for acc, locus in (("Q3M6C4", "Ava_3857"), ("Q3M6C5", "Ava_3856"), ("Q3M6C6", "Ava_3855")):
    assert locus in rec(acc)["gene_names"] and not rec(acc)["reviewed"], acc
assert {"PF00501", "PF00550", "PF00975"} <= rec("Q3M6C6")["xrefs"]   # the NRPS-like enzyme: no fungal step has this fold

# ---------- decoys: (accession, gene_name, pmid or None, text, ECs on the entry, Rhea ids on the entry, verified) ----------
SPC = "sugar phosphate cyclase of the 3-dehydroquinate synthase fold (Pfam PF01761, PF24621); UniProt reviewed entry: "
OMTF = "S-adenosylmethionine-dependent O-methyltransferase (Pfam PF01596); UniProt reviewed entry: "
ATPF = "ATP-grasp fold protein (PROSITE profile PS50975); "
DECOY = {
 "ddgs": [
  ("H2K887", "valA", None, SPC + "2-epi-5-epi-valiolone synthase of validamycin biosynthesis (EC 4.2.3.152, RHEA:44184); same substrate (sedoheptulose 7-phosphate), another product. annotation only: no abstract was read for this entry", ["4.2.3.152"], ["RHEA:44184"], True),
  ("Q9ZAE9", "acbC", None, SPC + "2-epi-5-epi-valiolone synthase of acarbose biosynthesis (EC 4.2.3.152, RHEA:44184); same substrate, another product. annotation only: no abstract was read for this entry", ["4.2.3.152"], ["RHEA:44184"], True),
  ("C6WFL3", "Amir_2000", "22741921", SPC + "2-epi-valiolone synthase of Actinosynnema mirum (EC 4.2.3.155, RHEA:49564); same substrate, another product. assay read: the abstract of PMID 22741921 states the conversion of sedoheptulose 7-phosphate to 2-epi-valiolone by the Actinosynnema mirum enzyme. gene_name is the ordered locus name", ["4.2.3.155"], ["RHEA:49564"], True),
  ("P07639", "aroB", None, SPC + "3-dehydroquinate synthase of the shikimate pathway (EC 4.2.3.4, RHEA:21968); another substrate. annotation only: no abstract was read for this entry", ["4.2.3.4"], ["RHEA:21968"], True),
  ("P07547", "aromA", None, SPC + "pentafunctional AROM polypeptide of Aspergillus nidulans; its first domain is the 3-dehydroquinate synthase (EC 4.2.3.4, RHEA:21968; the entry also carries RHEA:21096, RHEA:17737, RHEA:13121 and RHEA:21256 for the other four domains, and no EC on the protein name). This is the fungal protein that a PF01761 search finds in every genome. annotation only for this entry; in Phaffia rhodozyma the deletion of the AROM gene does not change mycosporine production (full text of PMID 36983003)", [], ["RHEA:21968", "RHEA:21096", "RHEA:17737", "RHEA:13121", "RHEA:21256"], True),
  ("Q9S5E2", "btrC", None, SPC + "2-deoxy-scyllo-inosose synthase of butirosin biosynthesis (EC 4.2.3.124, RHEA:33071); another substrate. annotation only: no abstract was read for this entry", ["4.2.3.124"], ["RHEA:33071"], True),
 ],
 "omt": [
  ("O42898", "SPBC119.03", None, OMTF + "probable catechol O-methyltransferase 1 of Schizosaccharomyces pombe (EC 2.1.1.6, RHEA:17877). annotation only: the name says 'Probable' and the catalytic activity has no evidence code in UniProt; not assayed with demethyl-4-deoxygadusol. It is not excluded as a true member of the step: it is 30% identical to the Phaffia OMT (MEASURED), closer than MysB (17%), and no assay of this protein was read. The step needs tree placement; do not count a hit to this protein as a negative (USER DECISION 7)", ["2.1.1.6"], ["RHEA:17877"], True),
  ("B8NI24", "imqG", None, OMTF + "O-methyltransferase of imizoquin biosynthesis in Aspergillus flavus (EC 2.1.1.-, no Rhea reaction on the entry). annotation only; not assayed with demethyl-4-deoxygadusol", ["2.1.1.-"], [], True),
  ("A0A397HK53", "ankF", None, OMTF + "O-methyltransferase of NK13650 biosynthesis in Aspergillus thermomutatus (EC 2.1.1.-, RHEA:80259). annotation only: no abstract was read for this entry; not assayed with demethyl-4-deoxygadusol", ["2.1.1.-"], ["RHEA:80259"], True),
  ("P21964", "COMT", None, OMTF + "human catechol O-methyltransferase (EC 2.1.1.6, RHEA:17877, RHEA:53108, RHEA:53104, RHEA:53100, RHEA:53096, RHEA:53092, RHEA:53088). annotation only: no abstract was read for this entry; not assayed with demethyl-4-deoxygadusol. 29% identical to the Phaffia OMT (MEASURED)", ["2.1.1.6"], ["RHEA:17877", "RHEA:53108", "RHEA:53104", "RHEA:53100", "RHEA:53096", "RHEA:53092", "RHEA:53088"], True),
 ],
 "atp_grasp_ligase": [
  ("Q4WD45", "fsqD", None, ATPF + "UniProt reviewed entry: ATP-grasp enzyme of fumisoquin biosynthesis in Aspergillus fumigatus (EC 6.3.2.-, no Rhea reaction on the entry). annotation only; not assayed with 4-deoxygadusol", ["6.3.2.-"], [], True),
  ("C4JF81", "uncJ", None, ATPF + "UniProt reviewed entry: L-Ala-L-Gln ligase of isoindolamide biosynthesis in Uncinocarpus reesii (EC 6.3.2.-, RHEA:85751). annotation only: no abstract was read for this entry; not assayed with 4-deoxygadusol", ["6.3.2.-"], ["RHEA:85751"], True),
  ("A0A5N6TQK6", "aveB", None, ATPF + "UniProt reviewed entry: amide-bond-forming ATP-grasp enzyme of isoindolamide A biosynthesis in Aspergillus avenaceus (EC 6.3.2.-, RHEA:85807, RHEA:85783). annotation only: no abstract was read for this entry; not assayed with 4-deoxygadusol", ["6.3.2.-"], ["RHEA:85807", "RHEA:85783"], True),
  ("P07862", "ddlB", None, ATPF + "UniProt reviewed entry: D-alanine--D-alanine ligase B of Escherichia coli (EC 6.3.2.4, RHEA:11224). annotation only: no abstract was read for this entry; not assayed with 4-deoxygadusol", ["6.3.2.4"], ["RHEA:11224"], True),
  ("P04425", "gshB", None, ATPF + "UniProt reviewed entry: glutathione synthetase of Escherichia coli (EC 6.3.2.3, RHEA:13557). annotation only: no abstract was read for this entry; not assayed with 4-deoxygadusol", ["6.3.2.3"], ["RHEA:13557"], True),
  ("B2J6X6", "Npun_F5597", MYS_PMID, ATPF + "UniProt unreviewed entry named 'D-alanine--D-alanine ligase domain protein' (EC 6.3.2.4 from automatic annotation, no Rhea reaction on the entry): NpF5597 (mysD) of Nostoc punctiforme. assay read: the abstract of PMID 21890703 states that expression of NpF5597 converts mycosporine-glycine into shinorine and other iminomycosporines, which is another reaction (the second amino acid); it was not assayed with 4-deoxygadusol alone. The entry does not cite the paper and the label NpF5597 is the short form of the ordered locus name Npun_F5597, so the decoy is verified: false", ["6.3.2.4"], [], False),
 ],
}
DECOY_FAMILY = {"ddgs": {"PF01761", "PF24621"}, "omt": {"PF01596"}, "atp_grasp_ligase": {"PS50975"}}
FORBIDDEN_EC = {"ddgs": "4.2.3.154", "omt": None, "atp_grasp_ligase": None}
FORBIDDEN_RHEA = {"ddgs": {"RHEA:49560"}, "omt": set(), "atp_grasp_ligase": set()}
SUFFIX = {"ddgs": "; no EC 4.2.3.154 and no RHEA:49560 on the entry",
          "omt": "; the step has no EC and no Rhea reaction, so the check is the other annotation of the entry",
          "atp_grasp_ligase": "; the step has no EC and no Rhea reaction, so the check is the other annotation of the entry"}
ALL_REF = {a for refs in REF.values() for a, *_ in refs}
dprov = []
for step, ds in DECOY.items():
    recs = []
    for acc, gene, pmid, func, ecs, rhs, ver in ds:
        r = rec(acc)
        assert acc not in ALL_REF and gene in r["gene_names"], (acc, gene, r["gene_names"])
        assert r["reviewed"] == (acc != "B2J6X6"), acc
        assert DECOY_FAMILY[step] <= r["xrefs"], (step, acc, r["xrefs"])
        assert not (r["rhea"] & FORBIDDEN_RHEA[step]) and (FORBIDDEN_EC[step] is None or FORBIDDEN_EC[step] not in r["ec"]), acc
        assert sorted(ecs) == sorted(r["ec"]) and set(rhs) == r["rhea"], (acc, r["ec"], r["rhea"])
        assert set(re.findall(r"RHEA:\d+", func)) == set(rhs), acc
        assert set(re.findall(r"EC (\d+\.[\d-]+\.[\d-]+\.[\d-]+)", func)) == set(ecs) | ({"4.2.3.4"} if acc == "P07547" else set()), acc
        assert ("annotation only" in func) != ("assay read" in func), acc
        if pmid:
            assert pubmed_abstract(pmid)["title"] and ((pmid in r["pmids"]) == ver), acc
        recs.append(f">{r['header']}\n{wrap(r['sequence'])}\n")
        dprov.append([step, acc, r["species"], gene, pmid or "", func + SUFFIX[step], r["srcdb"] + ("; PubMed" if pmid else ""), TODAY, str(ver).lower()])
    (DEC / f"{step}.faa").write_text("".join(recs), encoding="utf-8")
import hashlib
tsv(ROOT / "tests" / "data" / "kb" / f"{PW}_sequences.tsv", ["accession", "source", "length", "sha256"],
    sorted({(a, rec(a)["source"], str(rec(a)["length"]), hashlib.sha256(rec(a)["sequence"].encode()).hexdigest())
            for a in ALL_REF | {d[0] for ds in DECOY.values() for d in ds}}))
tsv(DEC / "provenance.tsv", ["step_id", "accession", "species", "gene_name", "pmid", "reason", "source_database", "retrieved_on", "verified"], dprov)
assert "4.2.3.4" in [e["value"] for c in full("P07547")["comments"] if c["commentType"] == "CATALYTIC ACTIVITY" for e in [{"value": c["reaction"].get("ecNumber")}]]
assert ordered_locus_names("C6WFL3") == ["Amir_2000"] and "22741921" in rec("C6WFL3")["pmids"]
assert ordered_locus_names("B2J6X6") == ["Npun_F5597"]

# ---------- compounds: ChEBI number -> (KEGG id, how the KEGG id is tied) ----------
def kegg(cid):
    return raw("kegg", cid, f"https://rest.kegg.jp/get/{cid}", as_json=False)
CPD = {"15721": ("C05382", "chebi"), "132137": ("C21281", "name"), "144275": ("", ""), "144274": ("", ""), "7039": ("C10607", "chebi"),
       "172445": ("C22846", "name"), "131712": ("", ""), "35671": ("", ""), "182637": ("", "")}
PUBCHEM_NAME = {"131712": ("shinorine", False), "35671": ("porphyra-334", False), "182637": ("palythine", True)}  # same first block as ChEBI?
crows = []
for num, (kid, how) in CPD.items():
    c = chebi_compound(num)
    assert c["charge"] in (0, "0") and not re.search(r"[+-]\]", c["smiles"]), (num, c["charge"], c["smiles"])
    m = formula_mass(c["formula"])
    assert abs(m - c["monoisotopic_mass"]) < 0.0005, num
    if num in PUBCHEM_NAME:   # PubChem has no compound with the ChEBI InChIKey; compare by name
        try:
            pubchem_by_inchikey(c["inchikey"]); raise AssertionError(num)
        except KbFetchError:
            pass
        name, same = PUBCHEM_NAME[num]
        props = pubchem_name(name)
        assert all(p["MolecularFormula"] == c["formula"] for p in props), num
        assert any(p["InChIKey"][:14] == c["inchikey"][:14] for p in props) == same, num
    else:
        p = pubchem_by_inchikey(c["inchikey"])
        assert p["inchikey"][:14] == c["inchikey"][:14] and abs(m - p["exact_mass"]) < 0.0005, num
    r = raw("chebi_raw", num, f"https://www.ebi.ac.uk/chebi/backend/api/public/compound/{num}/")
    kx = {str(x.get("accession_number")) for lst in (r.get("database_accessions") or {}).values() for x in lst if re.fullmatch(r"C\d{5}", str(x.get("accession_number")))}
    if how == "chebi":
        assert kx == {kid}, (num, kx)
    else:
        assert kx == set(), (num, kx)
    if how == "name":   # neither database links the other: KEGG name and formula
        k = kegg(kid)
        assert f"FORMULA     {c['formula']}" in k and "ChEBI:" not in k, (num, kid)
    crows.append([f"CHEBI:{num}", c["inchikey"], c["inchikey"][:14], c["smiles"], c["formula"], f"{m:.4f}", kid, ""])
tsv(KB / "compound.tsv", COLUMNS["compound"], crows)
assert chebi_compound("7039")["name"] == "Mycosporine" and pubchem_name("mycosporine serinol")[0]["InChIKey"] == chebi_compound("7039")["inchikey"]
assert "NAME        Demethyl-4-deoxygadusol" in kegg("C21281") and "NAME        Gadusol" in kegg("C22846")
assert pubchem_name("gadusol")[0]["InChIKey"] == chebi_compound("172445")["inchikey"]
assert "NAME        Demethylgadusol" in kegg("C22843") and "ChEBI: 7039" in kegg("C10607")
# Rhea participants -> neutral parents (same atoms but H, same first block)
for charged, neutral, nh in [("57483", "15721", 2), ("132139", "132137", 1)]:
    x, y = chebi_compound(charged), chebi_compound(neutral)
    hx = int(re.search(r"H(\d+)", x["formula"]).group(1)); hy = int(re.search(r"H(\d+)", y["formula"]).group(1))
    assert hy - hx == nh and int(x["charge"]) == -nh and x["inchikey"][:14] == y["inchikey"][:14], charged
# compounds with no ChEBI entry (USER DECISIONS): PubChem knows the name, ChEBI has neither the key nor its first block
NO_CHEBI = {"4-deoxygadusol": "C8H12O5", "mycosporine-glycine": "C10H15NO6", "mycosporine-glutaminol-glucoside": "C19H32N2O11",
            "mycosporine glutamicol glucoside": "C19H31NO12", "mycosporine-alanine": "C11H17NO6"}
for name, formula in NO_CHEBI.items():
    props = pubchem_name(name)
    assert props[0]["MolecularFormula"] == formula, name
    assert not chebi_by_text(props[0]["InChIKey"]) and not chebi_by_text(props[0]["InChIKey"][:14]), name
for term in ("deoxygadusol", "gadusol", "demethylgadusol"):
    hits = {h["_source"]["chebi_accession"] for h in chebi_by_text(term)}
    assert hits == {"deoxygadusol": {"CHEBI:131711", "CHEBI:132139", "CHEBI:132052", "CHEBI:132137"}, "gadusol": {"CHEBI:228248"}, "demethylgadusol": {"CHEBI:228249"}}[term], (term, hits)
assert int(chebi_compound("228248")["charge"]) == -1 and int(chebi_compound("228249")["charge"]) == -1
for term in ("asterina-330", "mycosporine glycine"):
    assert chebi_by_text(term) == [], term

# ---------- reactions ----------
C = "CUR:mycosporine:"
SM = "Small molecules (phosphate, H2O, H+, S-adenosylmethionine, ATP, ADP) are not listed as participants."
NODG = "4-deoxygadusol has no ChEBI entry (search on 2026-11-21; PubChem CID 85926716), so it is not listed as a participant (USER DECISION 6, UNVERIFIED.md, D1)."
RX = [  # id, ec, note, substrates, products
 ("RHEA:49560", "4.2.3.154", f"Demethyl-4-deoxygadusol synthase: D-sedoheptulose 7-phosphate = (R)-demethyl-4-deoxygadusol + phosphate + H2O. Rhea lists the charged species CHEBI:57483 (2-) and CHEBI:132139 (1-); mapped to the neutral parents. KEGG KO K21342 is a cross-reference of the step. {SM}", [("15721", 1)], [("132137", 1)]),
 (C + "ddg_o_methylation", "", f"Curated: O-methylation of demethyl-4-deoxygadusol to 4-deoxygadusol with S-adenosylmethionine. Rhea has no reaction for it (searches for gadusol, deoxygadusol and mycosporine on 2026-11-21 found RHEA:49560, RHEA:50220 and the vertebrate gadusol reactions RHEA:78171 and RHEA:78167 only). {NODG} So the reaction has a substrate and no product in the table. {SM}", [("132137", 1)], []),
 (C + "mycosporine_glutaminol_formation", "", f"Curated: ATP-dependent ligation of 4-deoxygadusol with glutaminol to mycosporine-glutaminol (the fungal ATP-grasp ligase, PMID 36983003; the full text names glutaminol as the amine). {NODG} Glutaminol has no ChEBI entry that was found. So the reaction has a product and no substrate in the table, and path_possible reaches it only as the first step of a call. {SM}", [], [("144275", 1)]),
 (C + "mycosporine_glycine_formation", "", f"Curated, bacterial: ATP-dependent ligation of 4-deoxygadusol with glycine to mycosporine-glycine (MysC of cyanobacteria, PMID 21890703). No fungal protein is linked to this reaction. {NODG} Mycosporine-glycine has no ChEBI entry (PubChem CID 14444485). The reaction has no participant in the table.", [], []),
 (C + "mg_glucosylation", "", "Curated: glucosylation of mycosporine-glutaminol to mycosporine-glutaminol-glucoside. The gene is not identified in any fungus (T4). An enzyme activity is reported for Ascochyta fabae spores: a soluble glucosyltransferase that moves glucose from UDP-glucose to mycosporin-2, the amide form being the best acceptor (PMID 24258137, abstract). Whether mycosporin-2 (amide form) is mycosporine-glutaminol is UNVERIFIED. The compound is made by Phaffia rhodozyma and by strain CBS 6938 after integration of the cluster (PMID 36983003), and not by Saccharomyces cerevisiae with the three genes; the gene is not identified (T4). Mycosporine-glutaminol-glucoside has no ChEBI entry (PubChem CID 171990409), so no product is listed.", [("144275", 1)], []),
 (C + "mycosporine_glutamicol_formation", "", "Curated, hypothesis: formation of mycosporine-glutamicol, the compound with a carboxyl group where mycosporine-glutaminol has the amide. No paper that was read gives the enzyme or the precursor (ligation of another amine to 4-deoxygadusol, or hydrolysis of the amide of mycosporine-glutaminol). T4. No substrate is listed. The glucoside (mycosporine-glutamicol-glucoside, PMID 12717761, 30182254) has no ChEBI entry (PubChem CID 171990350).", [], [("144274", 1)]),
 (C + "shinorine_formation", "", "Curated, bacterial context, no step: attachment of serine to mycosporine-glycine gives the iminomycosporine shinorine (NpF5597 / MysD of Nostoc punctiforme, PMID 21890703; an NRPS-like enzyme in Anabaena variabilis, full text of PMID 20813918; actinomycete clusters, PMID 24907338). There is no fungal enzyme and no fungal report of the compound in the papers that were read. Mycosporine-glycine has no ChEBI entry, so no substrate is listed.", [], [("131712", 1)]),
 (C + "porphyra334_formation", "", "Curated, bacterial context, no step: attachment of threonine to mycosporine-glycine gives porphyra-334 (a minor product of NpF5597 / MysD of Nostoc punctiforme, PMID 21890703). There is no fungal enzyme and no fungal report of the compound in the papers that were read. Mycosporine-glycine has no ChEBI entry, so no substrate is listed.", [], [("35671", 1)]),
]
def rhea_search(q):
    return raw("rhea", q, "https://www.rhea-db.org/rhea/?query=" + q.replace(":", "%3A") + "&columns=rhea-id,equation,chebi-id,ec&format=tsv&limit=10", as_json=False)
assert "RHEA:49560\tD-sedoheptulose 7-phosphate = (R)-demethyl-4-deoxygadusol + phosphate + H2O + H(+)\tCHEBI:57483;CHEBI:132139;CHEBI:43474;CHEBI:15377;CHEBI:15378\tEC:4.2.3.154" in rhea_search("deoxygadusol")
assert {l.split("\t")[0] for l in rhea_search("deoxygadusol").strip().splitlines()[1:]} == {"RHEA:49560", "RHEA:50220"}
assert {l.split("\t")[0] for l in rhea_search("gadusol").strip().splitlines()[1:]} == {"RHEA:78171"}
assert len(rhea_search("mycosporine").strip().splitlines()) == 1 and len(rhea_search("shinorine").strip().splitlines()) == 1
assert "RHEA:78167" in rhea_search("valiolone")
assert "K21342" in raw("kegg_list", "ko_deoxygadusol", "https://rest.kegg.jp/find/ko/deoxygadusol", as_json=False)
tsv(KB / "reaction.tsv", COLUMNS["reaction"],
    [[rid, ec, "", "forward", "false" if rid.startswith("RHEA:") else "true", note, rid if rid.startswith("RHEA:") else ""] for rid, ec, note, s, p in RX])
tsv(KB / "reaction_participant.tsv", COLUMNS["reaction_participant"],
    [[rid, f"CHEBI:{c}", role, str(n)] for rid, ec, note, s, p in RX for role, lst in (("substrate", s), ("product", p)) for c, n in lst])

# ---------- step_reaction links ----------
# Tier rule (binding rule 1): T1a only when the quoted sentence names the protein of the link and states the assay or the
# product of that reaction. A bacterial protein is T2 at most, and only when the sentence names its gene. Else T3.
ID = ("Products were identified by HPLC retention time and absorption spectrum, with co-elution with the Phaffia rhodozyma compounds (full text); "
      "the full text reports no mass spectrometry. Protein: GenBank {} of cluster OQ547787, which cites the paper.")
SHORT = "The abstract writes {} for the ordered locus name {} of UniProt {} (short form); the entry does not cite the paper."
L = []  # step, reaction, accession, tier, assay, pmid, quote key or None
def link(*x):
    L.append(x)
link("ddgs", "RHEA:49560", "WEV88679.1", "T1c",
     "homozygous deletion in the fungus (Phaffia rhodozyma UCD 67-385): the ddgs omt double mutant and the single-gene mutants make no mycosporines. Phenotype only: no intermediate is reported for this gene, "
     "and the sentence names ddgs in the double mutant only. The product of this reaction, demethyl-4-deoxygadusol is not named in the abstract, so the link is not T1a. The same abstract reports 4-deoxygadusol "
     "after expression of DDGS together with OMT in Saccharomyces cerevisiae. USER DECISION: T1a if that co-expression product is accepted as the assay of DDGS. " + ID.format("WEV88679.1"), FUNGAL_PMID, "del")
link("ddgs", "RHEA:49560", "B2J6X9", "T2",
     "another organism (cyanobacterium Nostoc punctiforme ATCC 29133): expression of NpR5600, NpR5599 and NpR5598 in Escherichia coli gives mycosporine-glycine. " + SHORT.format("NpR5600 (mysA)", "Npun_R5600", "B2J6X9") +
     " The sentence names the gene; the product in it is that of the three genes together. The activity of this protein alone rests on a record: UniProt B2J6X9 (reviewed) gives EC 4.2.3.154 and RHEA:49560 "
     "with experimental evidence (ECO:0000269) from PMID 20813918 and 22741921.", MYS_PMID, "mys")
link("ddgs", "RHEA:49560", "Q3M6C3", "T3",
     "Held at T3. Another organism (cyanobacterium Anabaena variabilis ATCC 29413, UniProt name Trichormus variabilis). The quoted sentence names the enzyme class and does not name the protein Ava_3858. "
     "UniProt Q3M6C3 (reviewed) cites the paper for the catalytic activity (EC 4.2.3.154, RHEA:49560), and the full text (read) reports the activity of Ava_3858. Raise to T2 if the user accepts the UniProt record as the tie.", "22741921", "cmp")
link("omt", C + "ddg_o_methylation", "WEV88680.1", "T1a",
     "heterologous co-expression of the DDGS and OMT cDNAs of Phaffia rhodozyma in Saccharomyces cerevisiae gives 4-deoxygadusol, the product of this reaction. It is a co-expression: the abstract reports no expression of OMT alone, "
     "and the full text has no strain with DDGS alone, so the order of the two enzymes rests on the homology of DDGS. The DDGS link rests on the same experiment and is held at T1c for that reason (USER DECISION 3). " + ID.format("WEV88680.1"), FUNGAL_PMID, "het")
link("omt", C + "ddg_o_methylation", "B2J6X8", "T3",
     "Held at T3. Another organism (cyanobacterium Nostoc punctiforme ATCC 29133): expression of NpR5600, NpR5599 and NpR5598 in Escherichia coli gives mycosporine-glycine. " + SHORT.format("NpR5599 (mysB)", "Npun_R5599", "B2J6X8") +
     " The sentence names the gene and does not state the product of the methylation, and no database record supports the activity: UniProt B2J6X8 is unreviewed, has no catalytic activity record with evidence "
     "(EC 2.1.1.104 is automatic annotation) and exists 'by homology'. A sentence of the full text of PMID 20813918 (found by word search) reports 4-deoxygadusol from sedoheptulose 7-phosphate with NpR5600 and NpR5599 together in vitro; "
     "it is not an abstract sentence and it does not test NpR5599 alone. The same enzyme family is not shown: this protein is 17% identical to the Phaffia OMT (MEASURED) "
     "and has Pfam PF01596 where the Phaffia OMT has PF13578 at the gathering cutoff.", MYS_PMID, "mys")
link("atp_grasp_ligase", C + "mycosporine_glutaminol_formation", "WEV88681.1", "T1a",
     "heterologous co-expression of the DDGS, OMT and ATPG cDNAs of Phaffia rhodozyma in Saccharomyces cerevisiae gives mycosporine-glutaminol (MG); without ATPG the product is 4-deoxygadusol. "
     "In the fungus the atpg deletion accumulates 4-deoxygadusol (same abstract, another sentence). " + ID.format("WEV88681.1"), FUNGAL_PMID, "het")
link("atp_grasp_ligase", C + "mycosporine_glutaminol_formation", "A0A9Y1GGR7", "T3",
     "Inferred. ATPG of Phaffia rhodozyma strain VOH (UniProt A0A9Y1GGR7, mRNA OP021064): 490 of 500 residues equal to the assayed protein WEV88681.1 of strain UCD 67-385 (MEASURED). No assay of this strain was read, and the entry cites no paper.", "", None)
link("atp_grasp_ligase", C + "mycosporine_glycine_formation", "B2J6X7", "T2",
     "another organism (cyanobacterium Nostoc punctiforme ATCC 29133): expression of NpR5600, NpR5599 and NpR5598 in Escherichia coli gives mycosporine-glycine. " + SHORT.format("NpR5598 (mysC)", "Npun_R5598", "B2J6X7") +
     " The amine is glycine, not glutaminol: this is T2 for the enzyme family and not for the fungal reaction.", MYS_PMID, "mys")
acc_gene = {(s, a): g for s, refs in REF.items() for a, g, *_ in refs}
FUNGAL_ACC = {"WEV88679.1", "WEV88680.1", "WEV88681.1", "A0A9Y1GGR7"}
srows = []
for step, rid, acc, tier, assay, pmid, q in L:
    assert rid in {r[0] for r in RX} and (step, acc) in acc_gene, (step, rid, acc)
    qt = Q[q] if q else ""
    assert qt or tier not in ("T1a", "T1b", "T1c", "T2"), (step, acc, tier)
    assert not qt or (pubmed_abstract(pmid)["title"] and norm(qt) in norm(ab(pmid))), (step, acc, pmid)
    assert (tier == "T3") == assay.startswith(("Held at T3.", "Inferred.")), (step, acc)
    # a T1 tier needs a fungal protein and the fungal paper; a bacterial protein is T2 or T3
    assert (tier.startswith("T1")) <= (acc in FUNGAL_ACC and pmid == FUNGAL_PMID), (step, acc)
    assert (acc in FUNGAL_ACC) or tier in ("T2", "T3"), (step, acc)
    assert (acc in FUNGAL_ACC) == (species_key(rec(acc)["species"]) == "xanthophyllomyces dendrorhous"), acc
    srows.append([step, rid, rec(acc)["species"], acc_gene[(step, acc)], acc, tier, assay, pmid, qt])
# the gene named in each T1 / T2 sentence
for acc, word in (("WEV88679.1", "Δddgs"), ("WEV88680.1", "OMT"), ("WEV88681.1", "ATPG"), ("B2J6X9", "NpR5600"), ("B2J6X7", "NpR5598")):
    row = [r for r in srows if r[4] == acc][0]
    assert word in row[8] and row[5] != "T3", acc
assert "4-deoxygadusol or MG production" in Q["het"] and "Ava_3858" not in Q["cmp"]
assert set(Q.values()) == {r[8] for r in srows if r[8]}
# T4: no gene identified. The steps are in the YAML with refs: [] so that the route variants can be declared.
NO_ENZYME = "not assigned (no enzyme identified)"
T4 = [("mg_glucosyltransferase", C + "mg_glucosylation",
       "hypothesis for the gene only; no sequence is known, so no protein link is possible. An enzyme activity is published for one fungus: a soluble UDP-glucose:mycosporin-2 glucosyltransferase from spores of Ascochyta fabae "
       "(PMID 24258137, abstract read: the amide form is the best acceptor; pH optimum about 8.5; divalent cations needed). UNVERIFIED: whether 'mycosporin-2, amide form' is mycosporine-glutaminol, the aglycone of mycosporine-glutaminol-glucoside. "
       "The full text of PMID 36983003 assigns the glucosylation in Phaffia rhodozyma to a gene that is not part of the cluster and is absent from Saccharomyces cerevisiae; the gene is not identified"),
      ("glutamicol_formation", C + "mycosporine_glutamicol_formation",
       "hypothesis; no paper that was read gives an enzyme or a precursor for mycosporine-glutamicol in fungi")]
for sid, rid, assay in T4:
    srows.append([sid, rid, NO_ENZYME, "unidentified", "", "T4", assay, "", ""])
tsv(KB / "step_reaction.tsv", COLUMNS["step_reaction"], srows)
tsv(ROOT / "tests" / "data" / "kb" / f"{PW}_quotes.tsv", ["step_id", "reaction_id", "protein_accession", "pmid", "quote", "tier", "assay"],
    [[r[0], r[1], r[4], r[7], r[8], r[5], r[6]] for r in srows])

# ---------- fungal lineage table: organism, text of the abstract for the group, compound, pmid, evidence, note, fragments ----------
MGG, MGL, MGOH, MGOHG = "mycosporine-glutaminol-glucoside", "mycosporine-glutaminol", "mycosporine-glutamicol", "mycosporine-glutamicol-glucoside"
LIN = [
 ("Phaffia rhodozyma", "basidiomycetous yeast (Cystofilobasidiales)", MGG, "20955516", "abstract", "48 wild and collection strains were positive; the paper says Xanthophyllomyces dendrorhous", ["Xanthophyllomyces dendrorhous (Phaffia rhodozyma)", "The 48 wild and collection strains screened were positive", "mycosporine-glutaminol-glucoside (MGG)"]),
 ("Phaffia rhodozyma", "basidiomycetous yeast", MGL, FUNGAL_PMID, "full text", "wild-type strain UCD 67-385 (peak 1 of Figure 3B in the full text); the abstract shows MG only in Saccharomyces cerevisiae with DDGS, OMT and ATPG and in strain CBS 6938 with the cluster", ["with peak 1 corresponding to MG and peak 2 corresponding to MGG"]),
 ("Phaffia rhodozyma", "basidiomycetous yeast", MGG, FUNGAL_PMID, "abstract", "strain CBS 6938 with the integrated cluster makes MG and MGG", ["produced MG and mycosporine glutaminol glucoside"]),
 ("Phaffia rhodozyma (strain CBS 6938)", "basidiomycetous yeast", "none detected", FUNGAL_PMID, "abstract", "negative report: the wild-type strain makes no mycosporine and lacks the cluster (full text)", ["non-mycosporine-producing CBS 6938 wild-type strain"]),
 ("Phaffia rhodozyma", "basidiomycetous yeast", MGG, "38329645", "abstract", "", ["Phaffia rhodozyma is a basidiomycetous yeast", "mycosporine-glutaminol-glucoside (MGG)"]),
 ("Naganishia friedmannii", "yeast", MGG, "40036870", "abstract", "isolate FBU002, after PAR-UVR exposure; a gene cluster is reported, no gene is named in the abstract; the full text was not available", ["an isolate of Naganishia friedmannii", "mycosporine-glutaminol-glucoside (MGG)"]),
 ("Rhodotorula minuta", "basidiomycetous yeast", MGG, "15484286", "abstract", "", ["Rhodotorula minuta and R. slooffiae", "Mycosporine-glutaminol-glucoside"]),
 ("Rhodotorula slooffiae", "basidiomycetous yeast", MGG, "15484286", "abstract", "", ["Rhodotorula minuta and R. slooffiae", "Mycosporine-glutaminol-glucoside"]),
 ("Rhodotorula minuta", "red yeast", MGG, "14993945", "abstract", "after PAR or PAR + UVR exposure", ["Rhodotorula minuta, Rh. pinicola, and Rhodotorula sp.", "mycosporine-glutaminol-glucoside"]),
 ("Rhodotorula pinicola", "red yeast", MGG, "14993945", "abstract", "after PAR or PAR + UVR exposure", ["Rh. pinicola", "mycosporine-glutaminol-glucoside"]),
 ("Cryptococcus laurentii", "colourless yeast", MGG, "14993945", "abstract", "after PAR or PAR + UVR exposure", ["the colourless Cryptococcus laurentii", "mycosporine-glutaminol-glucoside"]),
 ("Cryptococcus steppossus", "basidiomycetous yeast", MGG, "21175346", "abstract", "nine isolates", ["nine isolates of the basidiomycetous yeast Cryptococcus steppossus", "mycosporine-glutaminol glucoside (MGG)"]),
 ("Pucciniomycotina (94 species tested)", "dimorphic basidiomycetes", MGG, "21744380", "abstract", "almost half of the species were positive; five of seven classes; species are not named in the abstract", ["Of the 94 fungal species (377 strains and 33 genera) tested, almost half were MYC-positive and MGG was the main compound produced"]),
 ("Cystobasidium keelungensis", "basidiomycetous yeast", MGG, "30182254", "abstract", "", ["Cystobasidium keelungensis", "mycosporine-glutaminol-glucoside and mycosporine-glutamicol-glucoside"]),
 ("Cystobasidium keelungensis", "basidiomycetous yeast", MGOHG, "30182254", "abstract", "", ["mycosporine-glutaminol-glucoside and mycosporine-glutamicol-glucoside"]),
 ("rock-inhabiting microcolonial fungi (species not named in the abstract)", "microcolonial ascomycetes", MGOHG, "12717761", "abstract", "", ["Microcolonial ascomycetes", "both mycosporine-glutamicol-glucoside and mycosporine-glutaminol-glucoside"]),
 ("rock-inhabiting microcolonial fungi (species not named in the abstract)", "microcolonial ascomycetes", MGG, "12717761", "abstract", "", ["both mycosporine-glutamicol-glucoside and mycosporine-glutaminol-glucoside"]),
 ("Hortaea werneckii", "melanized ascomycetous fungus", MGG, "18048939", "abstract", "stationary-phase cells", ["Hortaea werneckii", "stationary-phase cells also contained mycosporine-glutaminol-glucoside"]),
 ("Colletotrichum graminicola", "filamentous fungus", "mycosporine-glutamine", "31550526", "abstract", "germination self-inhibitor of falcate conidia; absent from oval conidia cultures", ["Colletotrichum graminicola", "self-inhibitor mycosporine-glutamine"]),
 ("Ascochyta fabae", "filamentous fungus", MGG, FUNGAL_PMID, "full text", "stated in the introduction of PMID 36983003 with two citations; one is REPORTED by the review to be PMID 24258137 (next row), the other was not read", []),
 ("Ascochyta fabae", "Fungi imperfecti (spores)", "mycosporin-2 (acceptor of a spore glucosyltransferase)", GT_PMID, "abstract", "enzyme activity, not a compound report: UDP-glucose:mycosporin-2 glucosyltransferase, amide form the best acceptor; UNVERIFIED whether mycosporin-2 (amide form) is mycosporine-glutaminol", ["mycosporin-2 glucosyltransferase from spores of Ascochyta fabae Speg. (Fungi imperfecti)", "the amide form being the best acceptor"]),
 ("Dermatocarpon luridum var. luridum", "chlorolichen (Verrucariaceae, Ascomycota)", MGOH, "30893758", "abstract", "lichen thallus: the producing partner is not stated in the abstract", ["Dermatocarpon luridum var. luridum", "(mycosporine glutamicol, mycosporine glutaminol)"]),
 ("Dermatocarpon luridum var. luridum", "chlorolichen (Verrucariaceae, Ascomycota)", MGL, "30893758", "abstract", "lichen thallus: the producing partner is not stated in the abstract", ["(mycosporine glutamicol, mycosporine glutaminol)"]),
 ("Cladonia arbuscula subsp. squarrosa", "lichen", MGOH, "35888756", "abstract", "lichen thallus: the abstract says that the partner that makes the compound is not identified", ["Cladonia arbuscula", "Mycosporine-glutamicol", "identification of a lichen partner actively involved"]),
 ("Lichina pygmaea", "cyanolichen", "mycosporine-serinol", "25800181", "abstract", "lichen thallus: the producing partner is not stated in the abstract", ["Lichina pygmaea", "mycosporine serinol"]),
 ("Degelia plumbea", "cyanolichen", MGOH, "21550091", "abstract", "lichen thallus: the producing partner is not stated in the abstract", ["mycosporine glutamicol in Degelia plumbea"]),
 ("Nephroma laevigatum", "cyanolichen", "mycosporine-hydroxyglutamicol", "21550091", "abstract", "lichen thallus: the producing partner is not stated in the abstract", ["From Nephroma laevigatum, the mycosporine hydroxyglutamicol was preparatively isolated"]),
 ("Collema cristatum", "lichenized ascomycete", "a mycosporine (not named in the abstract)", "14764094", "abstract", "lichen thallus: the producing partner is not stated in the abstract", ["A novel photo protective mycosporine was isolated from the lichenized ascomycete Collema cristatum"]),
 ("Metarhizium anisopliae var. anisopliae", "entomopathogenic fungus", "none detected", "20943158", "abstract", "negative report: conidia of strain ESALQ-1037, LC-MS", ["Metarhizium anisopliae", "search for mycosporines (the initial targets), but no compounds of this class were detected"]),
]
lrows = []
for org, group, cpd, pmid, ev, note, frags in LIN:
    for f in frags:
        assert f in (ab(pmid) if ev == "abstract" else ft), (pmid, f)
    assert ev == "abstract" or (pmid == FUNGAL_PMID and org.split(" (")[0] in ft), org
    lrows.append([org, group, cpd, pmid, ev, note])
tsv(REFS / "fungal_lineages.tsv", ["organism", "group", "compound", "pmid", "evidence", "note"], lrows)

# ---------- steps ----------
ORG = "A gene name is valid only with its organism."
STEPS = [  # id, name, ec, copy, role, optional, route, ko, lineage_note, family
 ("ddgs", "demethyl-4-deoxygadusol synthase (DDGS of Phaffia rhodozyma; sedoheptulose 7-phosphate cyclase)", "4.2.3.154", "single", "catalytic", True, None, "K21342",
  "First step of the fungal route: D-sedoheptulose 7-phosphate to demethyl-4-deoxygadusol (RHEA:49560). Fungal evidence: DDGS of Phaffia rhodozyma UCD 67-385 (GenBank WEV88679.1, PMID 36983003), T1c: the deletion strains make no mycosporine; "
  "the abstract names no product of this enzyme alone. The cyanobacterial enzymes (Npun_R5600 of Nostoc punctiforme, Ava_3858 of Anabaena variabilis) are T2 at most for a fungal ortholog and are never fungal evidence; they are 43 to 44% identical to the Phaffia protein, "
  "and 73% to each other (MEASURED 2026-11-21, global alignment). The enzyme is a sugar phosphate cyclase of the 3-dehydroquinate synthase fold (Pfam PF01761, PF24621; hmmscan of the Phaffia protein, Pfam 38.2). The family is broader than the step: "
  "the 3-dehydroquinate synthase domain of the fungal AROM protein, 2-epi-5-epi-valiolone synthase and 2-epi-valiolone synthase have the same Pfam ids and other products (decoys; the five single-domain decoys are 20 to 30% identical to the Phaffia protein, MEASURED; the five-domain AROM protein was not aligned). "
  "In Phaffia rhodozyma the AROM deletion does not change mycosporine production (full text of PMID 36983003). A step call needs tree placement with the decoys. The cluster is not in every strain: Phaffia rhodozyma CBS 6938 lacks it and makes no mycosporine, so absence is not a defect. "
  "Gene names: DDGS (Phaffia rhodozyma); mysA, NpR5600, Npun_R5600 (Nostoc punctiforme); Ava_3858 (Anabaena variabilis). " + ORG,
  "sugar phosphate cyclase, 3-dehydroquinate synthase fold"),
 ("omt", "O-methyltransferase (OMT of Phaffia rhodozyma; demethyl-4-deoxygadusol to 4-deoxygadusol)", None, "single", "catalytic", True, None, "",
  "Second step: O-methylation of demethyl-4-deoxygadusol to 4-deoxygadusol. Rhea has no reaction and no EC is assigned; 4-deoxygadusol has no ChEBI entry. Fungal evidence: OMT of Phaffia rhodozyma UCD 67-385 (GenBank WEV88680.1, PMID 36983003), T1a: "
  "co-expression with DDGS in Saccharomyces cerevisiae gives 4-deoxygadusol. The fungal enzyme is not shown to be an ortholog of the cyanobacterial one: the Phaffia OMT is 17% identical to Npun_R5599 (mysB) of Nostoc punctiforme and 19% to Ava_3857, "
  "and 29 to 30% identical to catechol O-methyltransferases of fission yeast and human, which are decoys (MEASURED 2026-11-21, global alignment). The Phaffia protein has Pfam PF13578 at the gathering cutoff and PF01596 below it (hmmscan, Pfam 38.2); "
  "the cyanobacterial proteins have PF01596. So the cyanobacterial link is held at T3, not T2 (no record supports the methylation by Npun_R5599 alone, and the same family is not shown), it is never fungal evidence, and the decoys are closer to the fungal reference than the bacterial reference is. The decoy O42898 (probable catechol O-methyltransferase of fission yeast, no evidence code on its activity) is not excluded as a true member. "
  "A step call needs the Phaffia protein as the query and tree placement with the decoys. Gene names: OMT (Phaffia rhodozyma); mysB, NpR5599, Npun_R5599 (Nostoc punctiforme). " + ORG,
  "class I S-adenosylmethionine-dependent O-methyltransferase"),
 ("atp_grasp_ligase", "ATP-grasp ligase (ATPG of Phaffia rhodozyma; 4-deoxygadusol and glutaminol to mycosporine-glutaminol)", None, "single", "catalytic", True, None, "",
  "Third step: ligation of an amine to 4-deoxygadusol. In the fungus the amine is glutaminol and the product is mycosporine-glutaminol: ATPG of Phaffia rhodozyma UCD 67-385 (GenBank WEV88681.1, PMID 36983003), T1a; the atpg deletion accumulates 4-deoxygadusol. "
  "In cyanobacteria the homologous ATP-grasp enzyme MysC (Npun_R5598 of Nostoc punctiforme, 28% identical to the Phaffia protein, MEASURED) ligates glycine and makes mycosporine-glycine; that link is T2 for the enzyme family, on another reaction, and is never fungal evidence. "
  "The cyanobacterial route goes on to a second amino acid with MysD (a D-Ala-D-Ala ligase-like protein, NpF5597) or with an NRPS-like enzyme (Ava_3855); that is not a step here: no paper that was read shows a second-amine step or an iminomycosporine in a fungus. "
  "MysD is a decoy of this step; the NRPS-like enzyme has another fold and is neither reference nor decoy. The Phaffia protein has no Pfam family at the gathering cutoff (hmmscan, Pfam 38.2; best hits PF02786 and PF02655 with E-values 0.05 and 0.13), "
  "so domains_any is empty; the UniProt entry of the ATPG of strain VOH (A0A9Y1GGR7, 98% identical) carries the PROSITE ATP-grasp profile PS50975. The decoys are ATP-grasp proteins with PS50975 and another reaction (15 to 17% identical, MEASURED). "
  "A step call needs the Phaffia protein as the query and tree placement. Gene names: ATPG (Phaffia rhodozyma); mysC, NpR5598, Npun_R5598 (Nostoc punctiforme). " + ORG,
  "ATP-grasp ligase"),
 ("mg_glucosyltransferase", "glucosylation of mycosporine-glutaminol (gene not identified)", None, "single", "catalytic", True, "mycosporine_glutaminol_glucoside", "",
  "T4: the gene is not identified. Mycosporine-glutaminol-glucoside is the main mycosporine of the yeasts in the lineage table. The enzyme activity is shown in one fungus and no sequence is known: a soluble UDP-glucose:mycosporin-2 glucosyltransferase from spores of Ascochyta fabae, "
  "with the amide form as the best acceptor, a pH optimum of about 8.5 and a need for divalent cations (PMID 24258137, abstract read). UNVERIFIED: whether 'mycosporin-2, amide form' is mycosporine-glutaminol, the aglycone of mycosporine-glutaminol-glucoside. "
  "UDP-glucose as the donor is a lead for the enzyme family (a UDP-glycosyltransferase), not a reference. The full text of PMID 36983003 states that the glucosylation is done by a gene that is not part of the cluster "
  "and whose function is absent from Saccharomyces cerevisiae; strain CBS 6938 with the integrated cluster makes the glucoside. No reference protein, no decoy, no domain. The cyanobacterial glycosyltransferase MysG1 (PMID 39247679) acts on a palythine and is not evidence for this step.",
  "unknown (no fungal enzyme identified)"),
 ("glutamicol_formation", "formation of mycosporine-glutamicol (enzyme and precursor not identified)", None, "single", "catalytic", True, "mycosporine_glutamicol", "",
  "T4: the enzyme is not identified. Mycosporine-glutamicol and its glucoside are reported in rock-inhabiting microcolonial fungi (PMID 12717761), in Cystobasidium keelungensis (PMID 30182254) and in lichen thalli (PMID 30893758, 35888756, 21550091). "
  "No paper that was read gives a gene, an enzyme or the precursor of the glutamicol side chain. No reference protein, no decoy, no domain.",
  "unknown (no fungal enzyme identified)"),
]
NOPROD = " In the tables the route gives no product: 4-deoxygadusol has no ChEBI entry, so the chain stops after the O-methyltransferase step (USER DECISION 6, UNVERIFIED.md, D1)."
ROUTES = [
 ("mycosporine_glutaminol", "direct route to mycosporine-glutaminol (Phaffia rhodozyma DDGS, OMT, ATPG)", ["ddgs", "omt", "atp_grasp_ligase"],
  "fungal_evidence: true. Three genes of one cluster in Phaffia rhodozyma UCD 67-385 (PMID 36983003): deletion of each gene stops mycosporine production, and the three cDNAs in Saccharomyces cerevisiae give mycosporine-glutaminol. "
  "Sedoheptulose 7-phosphate to demethyl-4-deoxygadusol to 4-deoxygadusol to mycosporine-glutaminol." + NOPROD),
 ("mycosporine_glutaminol_glucoside", "glucoside route to mycosporine-glutaminol-glucoside (Phaffia rhodozyma; glucosylation gene not identified)", ["ddgs", "omt", "atp_grasp_ligase", "mg_glucosyltransferase"],
  "fungal_evidence: partial. The first three steps are shown in Phaffia rhodozyma (PMID 36983003). The last step is T4: the glucosylation gene is not identified, so the route is never reported as possible." + NOPROD +
  " Mycosporine-glutaminol-glucoside has no ChEBI entry either."),
 ("mycosporine_glutamicol", "glutamicol route to mycosporine-glutamicol and its glucoside (enzymes not identified)", ["ddgs", "omt", "glutamicol_formation"],
  "fungal_evidence: partial. The compounds are reported in fungi and lichens (PMID 12717761, 30182254, 30893758); no gene is known for the glutamicol side chain or for the glucoside. The last step is T4. "
  "That the route shares DDGS and OMT with the glutaminol route is an assumption from the shared 4-deoxygadusol core; it is not shown in a glutamicol-producing fungus." + NOPROD),
]
step_rx = {}
for r in srows:
    step_rx.setdefault(r[0], [])
    if r[1] not in step_rx[r[0]]:
        step_rx[r[0]].append(r[1])
def y(s):  # YAML double-quoted scalar
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'
out = [textwrap.dedent('''\
    # Mycosporines and mycosporine-like amino acids, step set v1 (stage 1, task S6, 2026-11-21).
    #
    # Status: DRAFT for user review. Pathway-level verified: false.
    # Do not run a census against this file until the user has reviewed the step list.
    #
    # The step set starts from the fungal literature: one paper shows three genes of one cluster in the
    # yeast Phaffia rhodozyma UCD 67-385 (PMID 36983003; abstract and full text read). It is not a copy
    # of the cyanobacterial gene list. Cyanobacterial proteins are references for tree placement only:
    # they are T2 at most and never fungal evidence. No other fungal gene of the pathway was found in
    # the papers that were read; the glucosylation and the glutamicol steps are T4 (refs: []).
    #
    # Reference checks: refs/mycosporine/provenance.tsv, refs/mycosporine/decoys/provenance.tsv and
    # refs/mycosporine/UNVERIFIED.md. Evidence tiers, quotes and reactions: kb/mycosporine/*.tsv.
    # Reported producers: refs/mycosporine/fungal_lineages.tsv.
    #
    # A ref has verified: true only when, on 2026-11-21, (1) the UniProt entry was fetched and its
    # organism and gene name matched, (2) the UniProt entry carries the EC and Rhea reaction of the
    # step, and (3) the PMID exists in PubMed and is cited by the entry. It never means a person
    # reviewed it. The three Phaffia proteins are GenBank records with no UniProtKB entry, so they
    # are verified: false although the GenBank record cites the paper (USER DECISION).
    # A decoy has verified: true when the UniProt entry was fetched on 2026-11-21, its organism and
    # gene name matched, it carries the family ids of the step, and it carries neither the EC nor a
    # Rhea reaction of the step. A decoy pmid is given only when the abstract was read.
    #
    # Field notes:
    # - Routes: each route note starts with "fungal_evidence: true" or "fungal_evidence: partial"
    #   (text, not a schema field). In the tables no declared route gives a product, because
    #   4-deoxygadusol has no ChEBI entry.
    # - species: the organism name of the UniProt entry, copied exactly. For the GenBank proteins:
    #   the GenBank organism with the strain, "Phaffia rhodozyma (strain UCD 67-385)".
    # - gene_name: a gene name of the record. For the cyanobacterial entries it is the ordered locus
    #   name (Npun_R5600), the only gene identifier that UniProt gives. It is valid only with species.
    # - pmid of a ref: for a UniProt entry a paper cited by the entry; for a GenBank protein the paper
    #   cited by the GenBank record. The papers that show each reaction are in step_reaction.tsv.
    #   Exception: B2J6X8 and B2J6X7 carry PMID 21890703, which names their genes and which the
    #   entries do not cite (their only citation, PMID 23463784, is a sugar-transporter paper).
    # - domains_any: Pfam ids. For the Phaffia proteins they come from hmmscan (Pfam 38.2, gathering
    #   cutoff), because the proteins have no UniProt cross-references. atp_grasp_ligase has none.
    # - reactions: ids of kb/mycosporine/reaction.tsv. The tier of each (protein, reaction) link is
    #   in kb/mycosporine/step_reaction.tsv, not here.
    # - copy_expect, optional and role are PROPOSED values. No census data exists yet. Every step is
    #   optional: one Phaffia rhodozyma strain lacks the whole cluster, so absence is not a defect.
    # - KEGG: K21342 is a cross-reference in kb/mycosporine/step.tsv only.
    id: mycosporine
    verified: false
    routes:
    ''')]
for rid, name, sids, note in ROUTES:
    out.append(f"  - id: {rid}\n    name: {y(name)}\n    step_ids: [{', '.join(sids)}]\n    note: {y(note)}\n")
out.append("steps:\n")
krows = []
for sid, name, ec, copy, role, opt, route, ko, note, family in STEPS:
    doms = DOMAINS.get(sid, [])
    o = [f"  - id: {sid}", f"    name: {y(name)}", f"    ec: {y(ec) if ec else 'null'}"]
    if route:
        o.append(f"    route: {route}")
    if sid in REF:
        o.append("    refs:")
        for acc, gene, pmid, ver in REF[sid]:
            o += [f"      # {c}" for c in REF_NOTE.get((sid, acc), [])]
            o += [f"      - accession: {acc}", f"        species: {y(rec(acc)['species'])}", f"        gene_name: {y(gene)}",
                  f"        pmid: {y(pmid)}", f"        verified: {str(ver).lower()}"]
    else:
        o.append("    refs: []")
    if DECOY.get(sid):
        o.append("    decoys:")
        for acc, gene, pmid, func, ecs, rhs, ver in DECOY[sid]:
            o += [f"      - accession: {acc}", f"        species: {y(rec(acc)['species'])}", f"        gene_name: {y(gene)}",
                  f"        pmid: {y(pmid) if pmid else 'null'}", f"        verified: {str(ver).lower()}"]
    o += [f"    domains_any: [{', '.join(doms)}]", f"    copy_expect: {copy}", f"    role: {role}",
          "    fused_with: []", f"    optional: {str(opt).lower()}",
          "    reactions:"] + [f"      - {y(r)}" for r in step_rx[sid]] + [f"    lineage_note: {y(note)}", ""]
    out.append("\n".join(o) + "\n")
    krows.append([sid, route or "", ko, ";".join(doms), family, str(opt).lower()])
(ROOT / "pathways" / f"{PW}.yaml").write_text("".join(out).rstrip("\n") + "\n", encoding="utf-8")
tsv(KB / "step.tsv", COLUMNS["step"], krows)
n = {t: sum(1 for r in srows if r[5].startswith(t)) for t in ("T1", "T2", "T3", "T4")}
print("refs", len(prov), "decoy rows", len(dprov), "compounds", len(crows), "reactions", len(RX), "links", len(srows), n, "lineage rows", len(lrows))
