"""Build the L-DOPA melanin / other melanin curation files from checked database records.

Writes pathways/melanin_ldopa.yaml, refs/melanin_ldopa/ (FASTA, provenance, decoys),
kb/melanin_ldopa/*.tsv and tests/data/kb/melanin_ldopa_quotes.tsv. It does not write
UNVERIFIED.md. Every curated fact below is checked against the fetched record by an assert.
"""
from h import *
import csv, re, textwrap, urllib.parse
from pathway_census.kb import formula_mass, COLUMNS

ROOT = PKG
PW = "melanin_ldopa"
REFS = ROOT / "refs" / PW; DEC = REFS / "decoys"; KB = ROOT / "kb" / PW
TODAY = "2026-11-21"

def full(acc):
    return raw("uniprot_full", acc, f"https://rest.uniprot.org/uniprotkb/{acc}.json")

def header(acc):
    d = full(acc)
    db = "sp" if d["entryType"].startswith("UniProtKB reviewed") else "tr"
    pdn = d["proteinDescription"]
    nm = (pdn.get("recommendedName") or (pdn.get("submissionNames") or [{}])[0])["fullName"]["value"]
    gn = (d.get("genes") or [{}])[0]
    g = (gn.get("geneName") or (gn.get("orfNames") or [{}])[0] or {}).get("value")
    h = f"{acc} {db}|{acc}|{d['uniProtkbId']} {nm} OS={d['organism']['scientificName']} OX={d['organism']['taxonId']}"
    return h + (f" GN={g}" if g else "")

def srcdb(acc):
    d = full(acc); au = d["entryAudit"]
    kind = "Swiss-Prot" if d["entryType"].startswith("UniProtKB reviewed") else "TrEMBL"
    return f"UniProtKB/{kind} (entry v{au['entryVersion']}, sequence v{au['sequenceVersion']})"

def rheas(acc):
    out = set()
    for c in full(acc).get("comments", []):
        if c["commentType"] == "CATALYTIC ACTIVITY":
            out |= {x["id"] for x in c["reaction"].get("reactionCrossReferences", []) if x["database"] == "Rhea" and x["id"].startswith("RHEA:")}
    return out

def function_text(acc):
    return " ".join(t["value"] for c in full(acc).get("comments", []) if c["commentType"] == "FUNCTION" for t in c.get("texts", []))

def xref_ids(acc):
    return {x["id"] for x in full(acc).get("uniProtKBCrossReferences", []) if x["database"] in ("Pfam", "NCBIfam")}

def wrap(seq):
    return "\n".join(seq[i:i+60] for i in range(0, len(seq), 60))

def norm(s):
    return " ".join(s.split())

def ab(pmid):
    return pubmed_abstract(pmid)["abstract"]

def summary(pmids):
    return raw("pubmed_summary", ",".join(pmids), "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(pmids))["result"]

def quote(pmid, start, end):
    """Cut one sentence out of the fetched abstract. The text is never typed by hand."""
    a = ab(pmid)
    i = a.index(start); j = a.index(end, i) + len(end)
    q = a[i:j]
    assert q.endswith(".") and len(q) <= 400, (pmid, q)
    assert norm(q) in norm(a)
    assert not re.search(r"[.!?] [A-Z]", q[:-1]), (pmid, q)
    return q

# One sentence per link. A paper can give different sentences to different links (key pmid:name).
Q = {
 "2147580": quote("2147580", "A transformant with high constitutive tyrosinase levels", "1 mg/ml L-tyrosine."),
 "23053534": quote("23053534", "To compare the mechanism of hydroxylation by the acid-treated tyrosinase", "Clark-type oxygen electrode."),
 "22213164": quote("22213164", "The pro form tyrosinase itself showed no catalytic activity", "to induce tyrosinase activity."),
 "7893753": quote("7893753", "The predicted protein product derived from the melO gene", "activity determination."),
 "28500345": quote("28500345", "The proteolytically activated mushroom tyrosinase shows", "flavonols and chalcones."),
 "21598903": quote("21598903", "The H subunit originates from the ppo3 gene", "residues larger."),
 "17617709": quote("17617709", "Tyrosinase (monophenol, 3,4-dihydroxy L-phenylalanine", "purified to homogeneity."),
 "16151802": quote("16151802", "The best transformant obtained, A. niger D15#26-e", "ca. 20 mg l(-1)."),
 "8760791": quote("8760791", "To assess the effect of CNLAC1 deletion on virulence", "5' end of the gene."),
 "8300520": quote("8300520", "Dopamine was shown by mass spectroscopy", "an intermediate of melanin synthesis."),
 "17098915": quote("17098915", "HGA melanization was laccase dependent", "concentrations of glucose."),
 "15643074:lac1": quote("15643074", "LAC2 overexpression suppresses the melanin defects", "of these strains."),
 "15643074:lac2": quote("15643074", "Accordingly, a lac2 mutation results", "melanin formation."),
 "19028908:hppD": quote("19028908", "In contrast, homogentisic acid and pyomelanin were not observed", "hppD deletion mutant."),
 "19028908:hmgA": quote("19028908", "In the hmgA deletion strain, HmgA activity was abolished", "increased pigment formation."),
 "7673153": quote("7673153", "The HmgA protein catalyzes an essential step", "containing phenylalanine."),
 "27133313": quote("27133313", "We demonstrate that MelA produces aspulvinone E", "polymerization by TyrP."),
 "28791090": quote("28791090", "In our study, heterologous expression of all the A. terreus NRPS-like genes", "aspulvinone E."),
}

# ---------- references: step -> [(accession, gene_name, reference pmid, verified)] ----------
# The reference pmid is a paper cited by the UniProt entry. Papers that show the reaction are in
# the step_reaction links.
NO_GENE = {"Q2TL94": "tyrosinase"}  # UniProt entry without any gene name: protein name of the entry
REF = {
 "tyrosinase": [("P00440", "T", "2529259", True), ("A0A1S9DK56", "OAory_01107480", "32356371", True), ("Q00234", "melO", "7893753", True),
                ("C7FF05", "PPO4", "20676921", True), ("C7FF04", "PPO3", "21598903", True), ("A7BHQ9", "tyr1", "17617709", True),
                ("Q2TL94", "tyrosinase", "16151802", False)],
 "laccase_diphenol": [("Q55P57", "LAC1", "8300520", True), ("J9VY90", "LAC1", "12581355", True), ("J9VQZ4", "LAC2", "24743168", True)],
 "hppd": [("Q4WHU1", "hppD", "19028908", True)],
 "hmga": [("Q4WHT9", "hmgA", "19028908", True), ("Q00667", "hmgA", "7673153", True)],
 "aspulvinone_e_synthase": [("Q0CRX1", "melA", "27133313", True)],
 "asp_melanin_tyrosinase": [("Q0CRX0", "tyrP", "27133313", True)],
}
# What the UniProt entry of every reference of a step must carry (the entry supports the reaction).
REF_RHEA = {"tyrosinase": {"RHEA:18117", "RHEA:34287"}, "laccase_diphenol": {"RHEA:11276"}, "hppd": {"RHEA:16189"},
            "hmga": {"RHEA:15449"}, "aspulvinone_e_synthase": {"RHEA:63824"}, "asp_melanin_tyrosinase": {"RHEA:74195", "RHEA:74199"}}
REF_EC = {"tyrosinase": "1.14.18.1", "laccase_diphenol": "1.10.3.2", "hppd": "1.13.11.27", "hmga": "1.13.11.5"}
REF_NOTE = {
 ("tyrosinase", "P00440"): [
   "pmid 2529259 is the gene paper cited by the entry (no assay in the abstract). The T1c link rests on",
   "PMID 2147580 (same first and last author), which the entry does not cite."],
 ("tyrosinase", "A0A1S9DK56"): [
   "melB of Aspergillus oryzae (PMID 16233650, 22213164, 23053534). The UniProt entry has the ORF name",
   "OAory_01107480 only; the name melB is not in the entry. The tie is the sequence: PDB 3W6W (title: melB",
   "holo-protyrosinase) has this sequence at 615 of 616 residues, and the entry has the 616 residues that",
   "PMID 16233650 gives for melB. Q2UP46 (RIB40 gene model of the locus, 638 residues, 91% identical) is not",
   "the assayed sequence and is not a reference. Links held at T3. See UNVERIFIED.md, A2."],
 ("tyrosinase", "Q00234"): ["melO of PMID 7893753 (539 residues in the abstract and in the entry). Links held at T3."],
 ("tyrosinase", "C7FF05"): [
   "PPO4 ('Tyrosinase 4' in UniProt) = tyrosinase 4 / AbPPO4 of PMID 28500345, which the entry does not cite. Links held at T3.",
   "pmid 20676921 is the cDNA paper cited by the entry (the Escherichia coli product there had no activity)."],
 ("tyrosinase", "C7FF04"): ["PPO3: the H subunit of the crystallised tyrosinase complex (PMID 21598903). Links held at T3."],
 ("tyrosinase", "A7BHQ9"): ["UniProt organism: Pholiota microspora. The paper says Pholiota nameko (one NCBI taxid, 1538424). Links held at T3."],
 ("tyrosinase", "Q2TL94"): [
   "verified: false, for one reason only: the UniProt entry has NO gene name, so the gene-name rule cannot",
   "be met. gene_name below is the protein name of the entry; it is not a UniProt gene name. The entry",
   "cites PMID 16151802. UniProt organism: Trametes sanguinea (taxid 158606); the paper says Pycnoporus",
   "sanguineus, a name with no NCBI taxid hit. The user decides at the review gate (UNVERIFIED.md)."],
 ("laccase_diphenol", "Q55P57"): [
   "CNLAC1 of PMID 8300520 (624 residues in the abstract and in the entry). UniProt organism: Cryptococcus",
   "deneoformans (serotype D); the papers say Cryptococcus neoformans."],
 ("laccase_diphenol", "J9VY90"): [
   "LAC1 of strain H99 (serotype A). UniProt gives its function 'by similarity'. pmid 12581355 is cited by",
   "the entry for catalytic activity. Link held at T3: the strain of PMID 15643074 is not confirmed."],
 ("laccase_diphenol", "J9VQZ4"): [
   "LAC2 of strain H99 (CNAG_03464, next to LAC1 CNAG_03465). UniProt gives its function 'by similarity'.",
   "pmid 24743168 is the genome paper cited by the entry. Link held at T3 (strain of PMID 15643074 not confirmed)."],
 ("hmga", "Q00667"): ["Aspergillus nidulans hmgA: disruption gives accumulation of homogentisate (PMID 7673153)."],
 ("aspulvinone_e_synthase", "Q0CRX1"): ["apvA (Q0CWD0) of the same fungus does the same reaction and is not a reference of this step."],
 ("asp_melanin_tyrosinase", "Q0CRX0"): ["UniProt gene names: tyrP, melB (synonym). Links held at T3 (UNVERIFIED.md, B)."],
}

# ---------- decoys: (accession, gene_name, pmid or None, function text, ECs on the entry, Rhea ids on the entry) ----------
T3CU = "type-3 copper protein (Pfam PF00264) member; "
D_IVOB = ("Q5BGU9", "ivoB", "2126551", T3CU + "UniProt reviewed entry: N-acetyl-6-hydroxytryptophan oxidase of the conidiophore pigment of Aspergillus nidulans "
          "(EC 1.14.-.-, no Rhea reaction on the entry). The purified enzyme oxidises p-cresol but not tyrosine or 3,4-dihydroxyphenylalanine (abstract of PMID 2126551)",
          ["1.14.-.-"], [])
D_ORSC = ("Q5AUW8", "orsC", None, T3CU + "UniProt reviewed entry: tyrosinase-like protein of the orsellinic acid / F9775 gene cluster (EC 1.14.-.-, no Rhea reaction on the entry); "
          "function read from the UniProt record", ["1.14.-.-"], [])
D_PHOM = ("A0A8J9R8H5", "phomQ1", None, T3CU + "UniProt reviewed entry: tyrosinase-like protein of phomopsin biosynthesis (EC 1.14.18.-, no Rhea reaction on the entry); "
          "function read from the UniProt record", ["1.14.18.-"], [])
D_AOCO4 = ("Q2UNF9", "AO090001000383", "24043469", T3CU + "UniProt unreviewed entry with no EC and no Rhea reaction: catechol oxidase AoCO4 of Aspergillus oryzae. The entry cites "
           "PMID 24043469 (crystal structure; the abstract names AoCO4). The abstract of PMID 19798497 (AoCO4, 408 residues; the entry has 408 residues) reports "
           "no activity on L-tyrosine and on L-dopa", [], [])
D_GRIF = ("B1VTI5", "griF", "16282322", T3CU + "UniProt reviewed entry: o-aminophenol oxidase GriF of grixazone biosynthesis in Streptomyces griseus (EC 1.10.3.15, EC 1.10.3.4; RHEA:41420, RHEA:41424, RHEA:40963). "
          "The abstract of PMID 16282322 reports no monophenolase activity, so GriF is excluded from the monophenolase reaction only. It oxidises o-aminophenols in preference to catechol-type substrates: "
          "it is not shown to lack diphenol oxidation, and the abstract does not name L-DOPA",
          ["1.10.3.15", "1.10.3.4"], ["RHEA:41420", "RHEA:41424", "RHEA:40963"])
MCO = "multicopper oxidase (Pfam PF07732, PF00394, PF07731) member; UniProt reviewed entry: "
D_MCO = [
 ("P38993", "FET3", None, MCO + "iron transport ferroxidase (EC 1.16.3.1, RHEA:11148; EC 1.16.3.4, RHEA:30083)", ["1.16.3.1", "1.16.3.4"], ["RHEA:11148", "RHEA:30083"]),
 ("Q0CCX6", "gedJ", "7665560", MCO + "dihydrogeodin oxidase of geodin biosynthesis in Aspergillus terreus (EC 1.10.3.-, RHEA:64316; phenol oxidative coupling of dihydrogeodin)", ["1.10.3.-"], ["RHEA:64316"]),
 ("P0DYC3", "phlC", None, MCO + "multicopper oxidase of phlegmacin biosynthesis (EC 1.10.3.18, RHEA:86387; dimerisation of torosachrysone)", ["1.10.3.18"], ["RHEA:86387"]),
 ("Q8LPL3", "AAO", None, MCO + "L-ascorbate oxidase of Arabidopsis thaliana (EC 1.10.3.3, RHEA:30243)", ["1.10.3.3"], ["RHEA:30243"]),
 ("Q4WZB4", "abr1", None, MCO + "multicopper oxidase of the DHN-melanin cluster of Aspergillus fumigatus (EC 1.-.-.-, no Rhea reaction on the entry); a reference of step mco_abr1 of pathway melanin_dhn. The entry cites PMID 14970241, which reports a laccase-type activity on L-DOPA in conidial extracts and names no protein: abr1 is not shown to lack diphenol oxidation", ["1.-.-.-"], []),
]
DECOY = {
 "tyrosinase": [D_IVOB, D_ORSC, D_PHOM, D_AOCO4, D_GRIF],
 "laccase_diphenol": D_MCO,
 "asp_melanin_tyrosinase": [D_IVOB, D_ORSC, D_PHOM, D_AOCO4],
}
DECOY_XREF = {"tyrosinase": {"PF00264"}, "laccase_diphenol": {"PF07732", "PF00394", "PF07731"}, "asp_melanin_tyrosinase": {"PF00264"}}
# What a decoy of the step must not carry: the EC and the Rhea reactions of the step.
FORBIDDEN = {"tyrosinase": {"RHEA:18117", "RHEA:34287"}, "laccase_diphenol": {"RHEA:11276"}, "asp_melanin_tyrosinase": {"RHEA:74195", "RHEA:74199"}}
FORBIDDEN_EC = {"tyrosinase": "1.14.18.1", "laccase_diphenol": "1.10.3.2", "asp_melanin_tyrosinase": None}
SUFFIX = {"tyrosinase": "; no EC 1.14.18.1, no RHEA:18117 and no RHEA:34287 on the entry",
          "laccase_diphenol": "; no EC 1.10.3.2 and no RHEA:11276 on the entry",
          "asp_melanin_tyrosinase": "; no RHEA:74195 and no RHEA:74199 on the entry (the reactions of TyrP with aspulvinone E)"}
UNREVIEWED_DECOY = {"Q2UNF9"}
# Entries that carry the EC of a step and are neither references nor decoys (checked below).
NOT_DECOY = {"B8NM74": ("ustQ", "1.14.18.1"), "P17489": ("yA", "1.10.3.2"), "P06811": ("lacc", "1.10.3.2"), "Q70KY3": ("LAC1", "1.10.3.2")}

def species(acc):
    return uniprot_entry(acc)["organism"]

def check_entry(acc, gene):
    e = uniprot_entry(acc)
    assert e["accession"] == acc, acc
    if acc in NO_GENE:
        assert e["gene_names"] == [] and gene == NO_GENE[acc], (acc, e["gene_names"])
    else:
        assert gene in e["gene_names"], (acc, gene, e["gene_names"])
    return e

def tsv(path, cols, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar=None)
        w.writerow(cols)
        for r in rows:
            assert len(r) == len(cols), r
            assert not any("\t" in str(x) or "\n" in str(x) for x in r), r
            w.writerow(r)

# ---------- write reference FASTA and provenance ----------
REFS.mkdir(parents=True, exist_ok=True); DEC.mkdir(exist_ok=True); KB.mkdir(parents=True, exist_ok=True)
prov = []
for step, refs in REF.items():
    recs = []
    for acc, gene, pmid, ver in refs:
        e = check_entry(acc, gene)
        assert pubmed_abstract(pmid)["title"]
        assert pmid in e["pmids"], (step, acc, pmid, e["pmids"])
        assert REF_RHEA[step] <= rheas(acc), (step, acc, rheas(acc))       # the entry supports the reaction
        if step in REF_EC:
            assert REF_EC[step] in e["ec"], (step, acc, e["ec"])
        assert ver == (acc not in NO_GENE), acc
        recs.append(f">{header(acc)}\n{wrap(e['sequence'])}\n")
        prov.append([step, acc, e["organism"], gene, "", pmid, srcdb(acc) + "; PubMed (cited by the entry)", TODAY, str(ver).lower()])
    (REFS / f"{step}.faa").write_text("".join(recs), encoding="utf-8")
tsv(REFS / "provenance.tsv", ["step_id", "accession", "species", "gene_name", "region", "pmid", "source_database", "retrieved_on", "verified"], prov)

# ---------- ties between references and the papers of their links (each is an assert) ----------
# P00440: same first and last author in the gene paper (cited by the entry) and the overexpression paper
s = summary(["2147580", "2529259"])
a1, a2 = ([x["name"] for x in s[p]["authors"]] for p in ("2147580", "2529259"))
assert a1[0] == a2[0] == "Kupper U" and a1[-1] == a2[-1] == "Lerch K"
assert "2147580" not in uniprot_entry("P00440")["pmids"] and "Neurospora crassa" in ab("2147580") and "protyrosinase" in ab("2147580")
# A0A1S9DK56 = melB: PDB 3W6W (melB holo-protyrosinase) has its sequence; Q2UP46 is another gene model of the locus
pdb = raw("rcsb_entry", "3W6W", "https://data.rcsb.org/rest/v1/core/entry/3W6W")
assert "melB holo-protyrosinase" in pdb["struct"]["title"]
pseq = raw("rcsb_polymer", "3W6W_1", "https://data.rcsb.org/rest/v1/core/polymer_entity/3W6W/1")["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", "")
mel, rib = uniprot_entry("A0A1S9DK56"), uniprot_entry("Q2UP46")
assert len(pseq) == 620 and mel["length"] == 616 and rib["length"] == 638
assert sum(x == y for x, y in zip(pseq[-616:], mel["sequence"])) == 615 and pseq[-615:] == mel["sequence"][1:]
assert sum(x == y for x, y in zip(pseq[-616:], rib["sequence"])) < 600
# residues named in the abstracts of PMID 23749993 (Cys92, His94, Phe513, Cys522, Cys525) and 22213164 (Glu458, Ala616)
for pos, aa in ((92, "C"), (94, "H"), (458, "E"), (513, "F"), (522, "C"), (525, "C"), (616, "A")):
    assert mel["sequence"][pos - 1] == aa, pos
assert [rib["sequence"][q - 1] for q in (458, 513, 522, 525, 616)] == list("DLAEE")
for frag in ("Phe(513)", "Cys(522)", "Cys(525)"):
    assert frag in ab("23749993"), frag
assert "Glu458-Ala616" in ab("22213164")
assert "3W6W" in [x["id"] for x in full("Q2UP46")["uniProtKBCrossReferences"] if x["database"] == "PDB"]
assert mel["gene_names"] == ["OAory_01107480"] and "32356371" in mel["pmids"] and "with substrate l-tyrosine" in ab("32356371")
for q in ("22213164", "23053534", "16233650"):
    assert "melB" in ab(q) and "Aspergillus oryzae" in ab(q), q
assert "616 amino acid residues" in ab("16233650")
# Q00234 = melO; C7FF05 = tyrosinase 4; C7FF04 = ppo3
assert "539 amino acid residues" in ab("7893753") and uniprot_entry("Q00234")["length"] == 539
assert "tyrosinase 4 from Agaricus bisporus" in ab("28500345") and "AbPPO4" in pubmed_abstract("28500345")["title"]
assert "Tyrosinase 4" in [n["fullName"]["value"] for n in full("C7FF05")["proteinDescription"]["alternativeNames"]]
assert "28500345" not in uniprot_entry("C7FF05")["pmids"] and "without PPO activity" in ab("20676921")
assert "ppo3 gene" in ab("21598903") and "21598903" in uniprot_entry("C7FF04")["pmids"]
# Cryptococcus
e = uniprot_entry("Q55P57")
assert {"8300520", "8760791"} <= set(e["pmids"]) and "17098915" not in e["pmids"] and e["length"] == 624 and "624 amino acids" in ab("8300520")
assert "but not monophenolic groups such as in phenol, tyramine, or tyrosine (PubMed:8300520)" in function_text("Q55P57")
assert "By similarity" in function_text("J9VY90") and "By similarity" in function_text("J9VQZ4")
for acc in ("J9VY90", "J9VQZ4"):
    assert "15643074" not in uniprot_entry(acc)["pmids"]
assert "CNAG_03465" in uniprot_entry("J9VY90")["gene_names"] and "CNAG_03464" in uniprot_entry("J9VQZ4")["gene_names"]
assert "(LAC2) is adjacent to the LAC1 gene" in ab("15643074")
assert "produces melanin only from exogenous substrates" in ab("17098915") and "requires the availability of diphenolic precursors" in ab("17098915")
assert "including L-dopa, to quinones, which then polymerize to melanin" in ab("8760791")
# Asp-melanin
assert "28791090" not in uniprot_entry("Q0CRX1")["pmids"] and "atmelA" in ab("28791090")
assert "28791090" in uniprot_entry("Q0CRX0")["pmids"] and "28791090" in uniprot_entry("Q0CWD0")["pmids"]
assert "RHEA:63824" in rheas("Q0CWD0") and "melB" in uniprot_entry("Q0CRX0")["gene_names"]
assert "lacks genes for synthesis of DHN-melanin" in ab("27133313")
# pyomelanin
assert "crucial for the conversion of p-hydroxyphenylpyruvate to homogentisic acid" in ab("22046314")
assert "spontaneous oxidation and polymerisation of HGA" in ab("9570125")
assert "evolve spontaneously" in function_text("C7FF04")
assert "hmgA" in uniprot_entry("Q12649")["gene_names"] and "1.1.1.34" in uniprot_entry("Q12649")["ec"]
assert "LAC1" in uniprot_entry("G2WRK6")["gene_names"] and uniprot_entry("G2WRK6")["organism"].startswith("Verticillium dahliae")

# ---------- decoys ----------
ALL_REF = {a for refs in REF.values() for a, *_ in refs}
dprov = []
for step, ds in DECOY.items():
    recs = []
    for acc, gene, pmid, func, ecs, rhs in ds:
        assert acc not in ALL_REF, acc
        e = check_entry(acc, gene)
        assert e["reviewed"] != (acc in UNREVIEWED_DECOY), acc
        assert DECOY_XREF[step] <= xref_ids(acc), (step, acc, xref_ids(acc))
        assert not (rheas(acc) & FORBIDDEN[step]), (step, acc, rheas(acc))
        assert FORBIDDEN_EC[step] is None or FORBIDDEN_EC[step] not in e["ec"], (step, acc, e["ec"])
        assert sorted(ecs) == sorted(e["ec"]), (acc, ecs, e["ec"])
        assert set(rhs) == rheas(acc), (acc, rhs, rheas(acc))
        assert set(re.findall(r"RHEA:\d+", func)) == set(rhs) and set(re.findall(r"EC (\d+\.[\d-]+\.[\d-]+\.[\d-]+)", func)) == set(ecs), acc
        if "no Rhea reaction" in func:
            assert not rheas(acc), acc
        if pmid:
            assert pmid in e["pmids"] and pubmed_abstract(pmid)["title"]
        recs.append(f">{header(acc)}\n{wrap(e['sequence'])}\n")
        dprov.append([step, acc, e["organism"], gene, pmid or "", func + SUFFIX[step], srcdb(acc) + ("; PubMed" if pmid else ""), TODAY, "true"])
    (DEC / f"{step}.faa").write_text("".join(recs), encoding="utf-8")
tsv(DEC / "provenance.tsv", ["step_id", "accession", "species", "gene_name", "pmid", "reason", "source_database", "retrieved_on", "verified"], dprov)
# what the decoy reasons quote from abstracts
assert "but not tyrosine, 3,4-dihydroxyphenylalanine" in ab("2126551") and "ivoB" in ab("2126551")
assert "AoCO4" in ab("24043469") and "408 aa, AoCO4" in ab("19798497") and uniprot_entry("Q2UNF9")["length"] == 408
assert "no activity was detected on L-tyrosine and on L-dopa" in ab("19798497")
assert "GriF showed no monophenolase activity" in ab("16282322") and "as preferable substrates rather than catechol-type substrates" in ab("16282322")
assert "14970241" in uniprot_entry("Q4WZB4")["pmids"] and "laccase-type activity in conidial extracts" in ab("14970241")
assert "dihydrogeodin to (+)- geodin" in ab("7665560")
for acc, (gene, ec) in NOT_DECOY.items():  # carry the EC of a step: not decoys
    e = check_entry(acc, gene)
    assert ec in e["ec"] and acc not in ALL_REF, acc
assert {"RHEA:18117", "RHEA:34287"} <= rheas("B8NM74")

# ---------- compounds (ChEBI number of the neutral parent -> KEGG cross-reference found in the ChEBI record) ----------
KEGG = {"17895": "C00082", "15765": "C00355", "16852": "C00822", "15772": "C01693", "2003": "C04185", "18243": "C03758", "189002": "",
        "15999": "C01179", "44747": "C00544", "47904": "C01036", "17704": "C02006"}
crows = []
for num, kegg in KEGG.items():
    c = chebi_compound(num)
    assert c["charge"] in (0, "0"), (num, c["charge"])
    assert not re.search(r"[+-]\]", c["smiles"]), (num, c["smiles"])  # not a zwitterion
    m = formula_mass(c["formula"])
    p = pubchem_by_inchikey(c["inchikey"])
    assert p["inchikey"][:14] == c["inchikey"][:14], num
    assert abs(m - p["exact_mass"]) < 0.0005, num
    assert abs(m - c["monoisotopic_mass"]) < 0.0005, num
    r = raw("chebi_raw", num, f"https://www.ebi.ac.uk/chebi/backend/api/public/compound/{num}/")
    xs = {x.get("accession_number") for lst in (r.get("database_accessions") or {}).values() for x in lst}
    kx = {x for x in xs if re.fullmatch(r"C\d{5}", str(x))}
    assert kx == ({kegg} if kegg else set()), (num, kx)
    crows.append([f"CHEBI:{num}", c["inchikey"], c["inchikey"][:14], c["smiles"], c["formula"], f"{m:.4f}", kegg, ""])
tsv(KB / "compound.tsv", COLUMNS["compound"], crows)
# Rhea participants -> neutral parents. nh = protons added (0 for a zwitterion). Same atoms, same first block.
for charged, neutral, nh in [("58315", "17895", 0), ("57504", "15765", 0), ("57924", "16852", 0), ("57509", "15772", 1), ("16875", "2003", 1),
                             ("36242", "15999", 1), ("16169", "44747", 1), ("17105", "47904", 2), ("58240", "17704", 1)]:
    a, b = chebi_compound(charged), chebi_compound(neutral)
    ha = int(re.search(r"H(\d+)", a["formula"]).group(1)); hb = int(re.search(r"H(\d+)", b["formula"]).group(1))
    assert hb - ha == nh and re.sub(r"H\d+", "", a["formula"]) == re.sub(r"H\d+", "", b["formula"]), (charged, neutral)
    assert int(a["charge"]) == -nh and a["inchikey"][:14] == b["inchikey"][:14], (charged, a["charge"])
    if nh == 0:
        assert re.search(r"[+-]\]", a["smiles"]) and "zwitterion" in a["name"], charged
# the TyrP products: anions only, no neutral parent entry in ChEBI
for num in ("193114", "193115"):
    c = chebi_compound(num)
    assert int(c["charge"]) == -1, num
    r = raw("chebi_raw", num, f"https://www.ebi.ac.uk/chebi/backend/api/public/compound/{num}/")
    rel = r["ontology_relations"]["incoming_relations"] + r["ontology_relations"]["outgoing_relations"]
    assert not [x for x in rel if "conjugate" in x["relation_type"]], num
    hits = raw("chebi_search", c["inchikey"][:14], "https://www.ebi.ac.uk/chebi/backend/api/public/es_search/?term=" + c["inchikey"][:14] + "&page=1&size=8").get("results", [])
    assert {h["_source"]["chebi_accession"] for h in hits} <= {f"CHEBI:{num}"}, (num, hits)
# 'decarboxy dopachrome' (PMID 8300520) -> dopaminechrome: L-dopachrome minus CO2
assert abs(formula_mass(chebi_compound("15772")["formula"]) - formula_mass("CO2") - formula_mass(chebi_compound("189002")["formula"])) < 1e-6
assert chebi_compound("189002")["name"] == "dopaminechrome (keto form)"

# ---------- reactions ----------
C = "CUR:melanin_ldopa:"
SM = "Small molecules (O2, H2O, CO2, H+, ATP, AMP, diphosphate, acceptor) are not listed as participants."
GEN = "Rhea has only the generic laccase reaction RHEA:11276 (4 hydroquinone + O2 = 4 benzosemiquinone + 2 H2O, EC 1.10.3.2), whose participants are compound classes."
RX = [  # id, ec, note, substrates, products
 ("RHEA:18117", "1.14.18.1", f"Monophenolase reaction of tyrosinase. Rhea writes L-tyrosine + O2 = L-dopaquinone + H2O as one reaction: free L-DOPA is not a product (pathway review, O7). Rhea lists the zwitterions CHEBI:58315 and CHEBI:57924; mapped to the neutral parents. KEGG KO K00505 is a cross-reference of the step. {SM}", [("17895", 1)], [("16852", 1)]),
 ("RHEA:34287", "1.14.18.1", f"Diphenolase reaction of tyrosinase (2 L-dopa + O2 = 2 L-dopaquinone + 2 H2O). L-DOPA is a substrate here; no reaction of the table makes it. Rhea lists the zwitterions CHEBI:57504 and CHEBI:57924; mapped to the neutral parents. {SM}", [("15765", 2)], [("16852", 2)]),
 (C + "ldopa_oxidation_laccase", "1.10.3.2", f"Curated: oxidation of exogenous L-DOPA to L-dopaquinone by a laccase. {GEN} The product is written as the quinone, as in the abstract of PMID 8760791 (L-dopa 'to quinones, which then polymerize to melanin'); the semiquinone step is not shown. RHEA:34287 has the same net equation and is the tyrosinase reaction (EC 1.14.18.1); the two are kept apart so that no laccase is linked to a tyrosinase reaction. In the reference organism the substrate is exogenous: Cryptococcus neoformans produces melanin only from exogenous substrates (PMID 17098915). KEGG KO K00421 and K05909 are cross-references of the step. {SM}", [("15765", 1)], [("16852", 1)]),
 (C + "dopamine_oxidation_laccase", "1.10.3.2", f"Curated: oxidation of exogenous dopamine to decarboxy dopachrome by the Cryptococcus laccase (mass spectrometry, PMID 8300520). The name 'decarboxy dopachrome' of the abstract is mapped by the curator to dopaminechrome (keto form), CHEBI:189002 (C8H7NO2, the formula of L-dopachrome minus CO2); PROPOSED. The steps between dopamine and dopaminechrome are not shown in the abstract. {GEN} {SM}", [("18243", 1)], [("189002", 1)]),
 (C + "hga_oxidation_laccase", "1.10.3.2", f"Curated: laccase-dependent pigment formation from exogenous homogentisic acid in Cryptococcus neoformans (PMID 17098915). The pigment has no formula, so no product is listed. The abstract does not name the laccase gene. {GEN}", [("44747", 1)], []),
 (C + "dopaquinone_cyclisation", "", "Curated, non-enzymatic in this table: L-dopaquinone to L-dopachrome. No enzyme link above T4. UniProt states for the Agaricus bisporus tyrosinases (C7FF04, function text) that the o-quinones 'evolve spontaneously' to pigment intermediates. The single line from L-dopaquinone to L-dopachrome is a simplification: the intermediate leucodopachrome (CHEBI:60872) is not listed. The chemistry is from the design and the pathway review (RECALLED there); no abstract was read for it.", [("16852", 1)], [("15772", 1)]),
 ("RHEA:13041", "5.3.3.12", "L-dopachrome to 5,6-dihydroxyindole-2-carboxylate (dopachrome tautomerase; an animal enzyme, for example human DCT, UniProt P40126). For fungi the link is T4 only: UniProt has no reviewed fungal entry with EC 5.3.3.12 (search on 2026-11-21); the unreviewed fungal entries with this EC have the automatic protein name 'L-dopachrome isomerase'. One PubMed query for a fungal dopachrome tautomerase found none. Rhea lists the anions CHEBI:57509 and CHEBI:16875; mapped to the neutral parents.", [("15772", 1)], [("2003", 1)]),
 ("RHEA:16189", "1.13.11.27", f"4-hydroxyphenylpyruvate dioxygenase. The reaction is a step of tyrosine catabolism and the source of homogentisic acid for pyomelanin. Rhea lists the anions CHEBI:36242 and CHEBI:16169; mapped to the neutral parents. KEGG KO K00457 is a cross-reference of the step. {SM}", [("15999", 1)], [("44747", 1)]),
 ("RHEA:15449", "1.13.11.5", f"Homogentisate 1,2-dioxygenase: removes homogentisic acid (tyrosine catabolism). Not a pigment-forming reaction: in Aspergillus fumigatus the hmgA deletion raises pigment formation (PMID 19028908). Rhea lists CHEBI:16169 (1-) and CHEBI:17105 (2-); mapped to the neutral parents. KEGG KO K00451 is a cross-reference of the step. {SM}", [("44747", 1)], [("47904", 1)]),
 (C + "hga_autooxidation", "", "Curated, no enzyme link: oxidation and polymerisation of homogentisic acid to pyomelanin. The pigment has no formula, so no product is listed. The abstract of PMID 9570125 calls it spontaneous (homogentisic acid made by a fungal HPPD in Escherichia coli); in Aspergillus fumigatus the accumulation of homogentisic acid raises pigment formation (PMID 19028908). Whether a fungal enzyme takes part was not checked.", [("44747", 1)], []),
 ("RHEA:63824", "", f"Aspulvinone E synthesis from two molecules of 4-hydroxyphenylpyruvate by an NRPS-like enzyme. Rhea has no EC for this reaction. Rhea lists the anions CHEBI:36242 and CHEBI:58240; mapped to the neutral parents. {SM}", [("15999", 2)], [("17704", 1)]),
 ("RHEA:74195", "", f"TyrP reaction of UniProt Q0CRX0: aspulvinone E + O2 = a product hydroxylated on both rings (CHEBI:193114). ChEBI has this product as an anion only (charge -1) and has no neutral parent entry (no conjugate-acid relation in the record; an InChIKey first-block search on 2026-11-21 found the anion only), so no product is listed. Rhea has no EC for this reaction. Rhea lists aspulvinone E as CHEBI:58240 (1-); mapped to the neutral parent. {SM}", [("17704", 1)], []),
 ("RHEA:74199", "", f"TyrP reaction of UniProt Q0CRX0: aspulvinone E + O2 = an ortho-quinone product (CHEBI:193115) + H2O. ChEBI has this product as an anion only (charge -1) and has no neutral parent entry (no conjugate-acid relation in the record; an InChIKey first-block search on 2026-11-21 found the anion only), so no product is listed. Rhea has no EC for this reaction. Rhea lists aspulvinone E as CHEBI:58240 (1-); mapped to the neutral parent. {SM}", [("17704", 1)], []),
]
for rid, ec, *_ in RX:  # re-check each Rhea id and EC against the Rhea response
    if rid.startswith("RHEA:"):
        t = raw("rhea", "id_" + rid[5:], "https://www.rhea-db.org/rhea/?query=" + f"rhea%3A{rid[5:]}" + "&columns=rhea-id,equation,chebi-id,ec&format=tsv&limit=3", as_json=False)
        line = [l for l in t.splitlines() if l.startswith(rid + "\t")]
        assert line and (line[0].split("\t") + [""])[3] == (f"EC:{ec}" if ec else ""), (rid, t)
t = raw("rhea", "id_11276", "https://www.rhea-db.org/rhea/?query=rhea%3A11276&columns=rhea-id,equation,chebi-id,ec&format=tsv&limit=3", as_json=False)
assert "RHEA:11276\t4 hydroquinone + O2 = 4 benzosemiquinone + 2 H2O\t" in t and "EC:1.10.3.2" in t
# no reviewed fungal dopachrome tautomerase in UniProt; the human enzyme carries the EC
q = "(ec:5.3.3.12) AND (taxonomy_id:4751) AND (reviewed:true)"
t = raw("uniprot_search", q + "|10", "https://rest.uniprot.org/uniprotkb/search?query=" + urllib.parse.quote(q) +
        "&fields=accession,reviewed,id,protein_name,gene_names,organism_name,length,ec,lit_pubmed_id&format=tsv&size=10", as_json=False)
assert len(t.strip().splitlines()) == 1, t
assert "5.3.3.12" in uniprot_entry("P40126")["ec"]
tsv(KB / "reaction.tsv", COLUMNS["reaction"],
    [[rid, ec, "", "forward", "false" if rid.startswith("RHEA:") else "true", note, rid if rid.startswith("RHEA:") else ""] for rid, ec, note, s, p in RX])
tsv(KB / "reaction_participant.tsv", COLUMNS["reaction_participant"],
    [[rid, f"CHEBI:{c}", role, str(n)] for rid, ec, note, s, p in RX for role, lst in (("substrate", s), ("product", p)) for c, n in lst])

# ---------- step_reaction links ----------
# Tier rule (binding rule 1, applied as written in fix round 1): a link is T1a only when the quoted sentence names the
# protein of the link and states the assay or the product for the reaction; T1b deletion with intermediate analysis;
# T1c colour phenotype in a sentence that names the protein; else T3 (the sentence stays as context).
L = []  # step, reaction, accession, tier, assay, pmid, quote key (True = Q[pmid], str = Q[key], False = none)
def link(step, rid, acc, tier, assay, pmid, q=True):
    L.append((step, rid, acc, tier, assay, pmid, q))
T, D = "RHEA:18117", "RHEA:34287"
link("tyrosinase", T, "P00440", "T1c", "protyrosinase gene expressed from a metallothionein promoter in Neurospora crassa; a transformant with high constitutive tyrosinase levels forms melanin on medium with L-tyrosine. Colour phenotype of an overexpression strain (not a deletion); no product analysis. The abstract also reports activity assays without a substrate. UniProt P00440 does not cite the paper; the tie is the organism and the protyrosinase gene of PMID 2529259 (same first and last author)", "2147580")
link("tyrosinase", D, "P00440", "T3", "Inferred. No abstract that was read states an assay of the Neurospora crassa enzyme with L-DOPA or a diphenol. UniProt P00440 carries RHEA:34287 without a cited assay paper.", "", False)
link("tyrosinase", T, "A0A1S9DK56", "T3", "Held at T3. Recombinant melB pro-tyrosinase activated by acid treatment; steady-state kinetic analysis of the phenolase activity (phenol hydroxylation) by O2 consumption. The quoted sentence does not name melB and names no substrate. A0A1S9DK56 is tied to melB through PDB 3W6W (UNVERIFIED.md, A2).", "23053534")
link("tyrosinase", D, "A0A1S9DK56", "T3", "Held at T3. The abstract reports tyrosinase activity of the trypsin-activated recombinant melB enzyme; the quoted sentence names no substrate.", "22213164")
MELO = "Held at T3. melO cDNA expressed in yeast; the product is identified by Western blotting and activity determination. The abstract names no substrate and no substrate class."
link("tyrosinase", T, "Q00234", "T3", MELO, "7893753")
link("tyrosinase", D, "Q00234", "T3", MELO, "7893753")
PPO4 = ("Held at T3. Tyrosinase 4 (AbPPO4) produced in Escherichia coli and activated by proteinase K; the activated enzyme accepts mono- and diphenols. The quoted sentence does not name tyrosinase 4 "
        "(it says 'the proteolytically activated mushroom tyrosinase') and names neither L-tyrosine nor L-DOPA. UniProt C7FF05 (PPO4, 'Tyrosinase 4') does not cite the paper.")
link("tyrosinase", T, "C7FF05", "T3", PPO4, "28500345")
link("tyrosinase", D, "C7FF05", "T3", PPO4, "28500345")
PPO3 = "Held at T3. Crystal structure of the Agaricus bisporus tyrosinase complex; the H subunit is the ppo3 product. The abstract reports no assay of PPO3. Kinetic papers on commercial mushroom tyrosinase (for example PMID 10552675, with L-tyrosine and L-DOPA) do not name the isoform."
link("tyrosinase", T, "C7FF04", "T3", PPO3, "21598903")
link("tyrosinase", D, "C7FF04", "T3", PPO3, "21598903")
PHO = "Held at T3. Enzyme purified to homogeneity from fruit bodies of Pholiota nameko and named by its systematic name; the abstract reports no assay result and no product. UniProt A7BHQ9 cites the paper for catalytic activity."
link("tyrosinase", T, "A7BHQ9", "T3", PHO, "17617709")
link("tyrosinase", D, "A7BHQ9", "T3", PHO, "17617709")
TRA = ("Held at T3. cDNA over-expressed in Aspergillus niger; the recombinant enzyme has monophenolase and diphenolase activities (534 and 1,668 U per litre). The quoted sentence names the activity classes "
       "and neither L-tyrosine nor L-DOPA, and the UniProt entry has no gene name. The paper says Pycnoporus sanguineus; UniProt Q2TL94 (Trametes sanguinea) cites the paper.")
link("tyrosinase", T, "Q2TL94", "T3", TRA, "16151802")
link("tyrosinase", D, "Q2TL94", "T3", TRA, "16151802")
LD, DO, HG = C + "ldopa_oxidation_laccase", C + "dopamine_oxidation_laccase", C + "hga_oxidation_laccase"
HGA = "Held at T3. Pigment formation from homogentisic acid is laccase dependent in Cryptococcus neoformans strains of all four serotypes; the abstract does not name the laccase gene or the strain."
link("laccase_diphenol", LD, "Q55P57", "T1c", "disruption of the 5' end of CNLAC1 gives a Mel- mutant, and complementation of another Mel- mutant with CNLAC1 restores the Mel+ phenotype. Colour phenotype only. The abstract states as background that the laccase oxidises L-dopa to quinones; it reports no product analysis", "8760791")
link("laccase_diphenol", DO, "Q55P57", "T3", "Held at T3. Purified enzyme (diphenol oxidase, identified as a laccase; gene CNLAC1 cloned); the oxidation product of dopamine was identified by mass spectrometry. The quoted sentence does not name the enzyme. UniProt Q55P57 cites the paper for the sequence, the function and the catalytic activity.", "8300520")
link("laccase_diphenol", HG, "Q55P57", "T3", HGA, "17098915")
STRAIN = " The abstract names no strain and no serotype, and UniProt {} does not cite the paper, so the protein of the paper is not confirmed. A fetch of the full text (Europe PMC, PMC544166) failed on 2026-11-21."
link("laccase_diphenol", LD, "J9VY90", "T3", "Held at T3. lac1 mutants have melanin defects, which LAC2 overexpression suppresses (colour phenotype)." + STRAIN.format("J9VY90 (LAC1 of strain H99)"), "15643074", "15643074:lac1")
link("laccase_diphenol", HG, "J9VY90", "T3", HGA, "17098915")
link("laccase_diphenol", LD, "J9VQZ4", "T3", "Held at T3. A lac2 mutation gives a modest delay in melanin formation (colour phenotype). The abstract places LAC2 adjacent to LAC1 (ORF names CNAG_03464 and CNAG_03465 in H99)." + STRAIN.format("J9VQZ4 (LAC2 of strain H99)"), "15643074", "15643074:lac2")
link("hppd", "RHEA:16189", "Q4WHU1", "T1b", "gene deletion in the fungus; no homogentisic acid and no pyomelanin in the hppD deletion mutant (intermediate analysis)", "19028908", "19028908:hppD")
link("hmga", "RHEA:15449", "Q4WHT9", "T1b", "gene deletion in the fungus; HmgA activity abolished and homogentisic acid accumulates (intermediate analysis)", "19028908", "19028908:hmgA")
link("hmga", "RHEA:15449", "Q00667", "T1b", "gene disruption in the fungus; homogentisate accumulates in broths with phenylalanine (intermediate analysis). The same abstract reports homogentisate dioxygenase activity of the purified fusion protein made in Escherichia coli, in a sentence that does not name HmgA", "7673153")
link("aspulvinone_e_synthase", "RHEA:63824", "Q0CRX1", "T1a", "the abstract states the product of MelA (aspulvinone E); it does not state the method. UniProt Q0CRX1 cites the paper for catalytic activity", "27133313")
link("aspulvinone_e_synthase", "RHEA:63824", "Q0CRX1", "T1a", "heterologous expression of the Aspergillus terreus NRPS-like genes; atmelA and apvA both release aspulvinone E. UniProt Q0CRX1 (melA) does not cite the paper; the entries of tyrP (Q0CRX0) and apvA (Q0CWD0) do", "28791090")
TYRP = "Held at T3. The abstract states that TyrP activates aspulvinone E for polymerisation; it states no assay and no product of TyrP. UniProt Q0CRX0 cites the paper for catalytic activity and biophysicochemical properties (RHEA:74195, RHEA:74199). Raise to T1a if the user accepts the UniProt record as the evidence."
link("asp_melanin_tyrosinase", "RHEA:74195", "Q0CRX0", "T3", TYRP, "27133313")
link("asp_melanin_tyrosinase", "RHEA:74199", "Q0CRX0", "T3", TYRP, "27133313")

acc_gene = {(s, a): g for s, refs in REF.items() for a, g, *_ in refs}
srows = []
for step, rid, acc, tier, assay, pmid, q in L:
    assert rid in {r[0] for r in RX} and (step, acc) in acc_gene, (step, rid, acc)
    if pmid:
        assert pubmed_abstract(pmid)["title"]
    qt = Q[q] if isinstance(q, str) else (Q[pmid] if q else "")
    assert qt or tier not in ("T1a", "T1b", "T1c", "T2"), (step, acc, tier)
    assert not qt or norm(qt) in norm(ab(pmid)), (step, acc, pmid)
    assert (tier == "T3") == assay.startswith(("Held at T3.", "Inferred.")), (step, acc)
    srows.append([step, rid, species(acc), acc_gene[(step, acc)], acc, tier, assay, pmid, qt])
assert set(Q.values()) == {r[8] for r in srows if r[8]}
# T4: dopaquinone onward. No reference protein; the step is in step.tsv only (carotenoid convention).
T4_STEP, NO_ENZYME = "dopachrome_onward", "not assigned (no enzyme identified)"
T4_ASSAY = "hypothesis; non-enzymatic in this table. No fungal enzyme was found (UniProt has no reviewed fungal entry with EC 5.3.3.12; one PubMed query)"
for rid in (C + "dopaquinone_cyclisation", "RHEA:13041"):
    srows.append([T4_STEP, rid, NO_ENZYME, "unidentified", "", "T4", T4_ASSAY, "", ""])
tsv(KB / "step_reaction.tsv", COLUMNS["step_reaction"], srows)
# change-detector fixture for the tests (all links)
tsv(ROOT / "tests" / "data" / "kb" / f"{PW}_quotes.tsv", ["step_id", "reaction_id", "protein_accession", "pmid", "quote", "tier", "assay"],
    [[r[0], r[1], r[4], r[7], r[8], r[5], r[6]] for r in srows])

# ---------- steps ----------
FLAG = "exogenous_substrate: true"
SEP = ("Tyrosinase and laccase both oxidise L-DOPA to the quinone, and they are not related in sequence. The annotation that separates them: a tyrosinase is EC 1.14.18.1 with "
       "RHEA:18117 (L-tyrosine to L-dopaquinone) and RHEA:34287 (L-DOPA to L-dopaquinone), in the type-3 copper family (Pfam PF00264); a laccase is EC 1.10.3.2 with the generic "
       "reaction RHEA:11276 only, in the multicopper oxidase family (Pfam PF07732, PF00394, PF07731). A laccase does not hydroxylate tyrosine: UniProt states for LAC1 of Cryptococcus "
       "(Q55P57, from PMID 8300520) that it does not oxidise monophenols such as tyrosine. ")
ORG = "A gene name is valid only with its organism."
STEPS = [  # id, name, ec, domains, copy, role, optional, route, ko, lineage_note, family
 ("tyrosinase", "tyrosinase (L-tyrosine and L-DOPA to L-dopaquinone)", "1.14.18.1", ["PF00264"], "multi", "catalytic", True, "tyrosinase_dopa", "K00505",
  SEP + "A PF00264 hit is not a tyrosinase call. MEASURED on 2026-11-21 (global alignment): the references are 20 to 55% identical to each other; the decoys (the short tyrosinase-like proteins ivoB, orsC and phomQ1, "
  "the catechol oxidase AoCO4, the bacterial GriF) and TyrP of Asp-melanin are 12 to 18% identical to the references and 20 to 45% to each other. UniProt gives EC 1.14.18.1 also to tyrosinases of other pathways, "
  "for example ustQ of ustiloxin biosynthesis (B8NM74); such entries are neither decoys nor references. Rhea writes L-tyrosine to L-dopaquinone as one reaction, so free L-DOPA is not a product of the table. "
  "From L-dopaquinone onward (L-dopachrome and later) the table is non-enzymatic; any fungal enzyme there is T4. No tyrosinase link is above T1c: no abstract that was read names the protein and L-tyrosine or L-DOPA in one sentence (assays of substrate classes are held at T3). "
  "Gene names: T (Neurospora crassa), melO and melB (Aspergillus oryzae; UniProt A0A1S9DK56 has the ORF name OAory_01107480 only), PPO3 and PPO4 (Agaricus bisporus), tyr1 (Pholiota microspora). "
  "melB is also a UniProt synonym of tyrP of Aspergillus terreus (Q0CRX0, Asp-melanin). " + ORG,
  "type-3 copper protein family (tyrosinase)"),
 ("laccase_diphenol", "laccase, diphenol oxidation only (LAC1 and LAC2 of Cryptococcus)", "1.10.3.2", ["PF07732", "PF00394", "PF07731"], "multi", "catalytic", True, "laccase_exogenous", "K00421;K05909",
  SEP + FLAG + ". Cryptococcus neoformans produces melanin only from exogenous substrates and needs diphenolic precursors (PMID 17098915), so the substrate of this step comes from outside and the fungus is not shown to make L-DOPA. "
  "The compound that enters the table depends on the substrate supplied: L-DOPA, dopamine (oxidised to decarboxy dopachrome, PMID 8300520) or homogentisic acid (PMID 17098915). "
  "The step is not specific to melanin: every laccase (EC 1.10.3.2) has the same generic reaction, so laccases with other roles (lignin laccases, yA of Aspergillus nidulans) are not decoys, and a laccase hit does not show melanin formation. "
  "The decoys are multicopper oxidases with another annotated reaction (22 to 29% identical to the references, MEASURED 2026-11-21). LAC1 and LAC2 of Cryptococcus neoformans H99 are 72% identical; LAC1 of serotype D "
  "(Cryptococcus deneoformans) and LAC1 of serotype A are 90% identical. Gene names: LAC1 names the diphenol laccase in Cryptococcus and a DHN-melanin cluster oxidase in Verticillium dahliae (G2WRK6). " + ORG,
  "multicopper oxidase family (laccase)"),
 ("hppd", "4-hydroxyphenylpyruvate dioxygenase (homogentisic acid for pyomelanin; hppD of Aspergillus fumigatus)", "1.13.11.27", ["TIGR01263"], "single", "catalytic", True, None, "K00457",
  "HppD makes homogentisic acid. Pyomelanin forms when homogentisic acid accumulates and oxidises: in Aspergillus fumigatus the hppD deletion gives no homogentisic acid and no pyomelanin, and the hmgA deletion gives more pigment (PMID 19028908). "
  "hppD and hmgA are enzymes of tyrosine catabolism, so their presence shows catabolism, not pigment production; the census cannot tell the two apart. How widely the two genes occur in fungi was not measured here. "
  "In Aspergillus fumigatus HmgX is needed for the HppD reaction (PMID 22046314); HmgX is not a step here. The step before (L-tyrosine to 4-hydroxyphenylpyruvate) is not curated. "
  "The references also carry Pfam PF00903 (glyoxalase family), which is broader and is not listed.",
  "4-hydroxyphenylpyruvate dioxygenase family"),
 ("hmga", "homogentisate 1,2-dioxygenase (removes homogentisic acid; hmgA of Aspergillus fumigatus)", "1.13.11.5", ["TIGR01015"], "single", "catalytic", True, "tyrosine_degradation", "K00451",
  "HmgA removes homogentisic acid (tyrosine catabolism). It is not a pigment-forming step: pyomelanin rises when it is lost (hmgA deletion of Aspergillus fumigatus, PMID 19028908). "
  "In a capability table its presence argues against accumulation of the pigment precursor, and its absence is not a defect. hmgA of Aspergillus fumigatus and hmgA of Aspergillus nidulans are 89% identical (MEASURED 2026-11-21). "
  "The gene name hmgA names HMG-CoA reductase in Phycomyces blakesleeanus (UniProt Q12649). " + ORG,
  "homogentisate 1,2-dioxygenase family"),
 ("aspulvinone_e_synthase", "aspulvinone E synthetase, NRPS-like (melA of Aspergillus terreus)", None, ["PF00975"], "single", "catalytic", True, "asp_melanin", "",
  "NRPS-like enzyme (adenylation, thiolation and thioesterase domains) that makes aspulvinone E from 4-hydroxyphenylpyruvate (RHEA:63824). Asp-melanin is described for Aspergillus terreus, which lacks the genes for DHN-melanin (PMID 27133313); "
  "the range in other lineages is not checked, so absence is not a defect. apvA of Aspergillus terreus (Q0CWD0, 65% identical to melA, MEASURED 2026-11-21) does the same reaction in hyphae, where aspulvinone E goes to aspulvinones and not to melanin (PMID 28791090). "
  "So a hit to this step does not show Asp-melanin without TyrP, and apvA cannot be a decoy. The references carry Pfam PF00501, PF00550 and PF00975; all three are broad families (PF00975 is listed), so the call needs tree placement.",
  "NRPS-like enzyme (A-T-TE domains)"),
 ("asp_melanin_tyrosinase", "tyrosinase of Asp-melanin (tyrP of Aspergillus terreus; acts on aspulvinone E)", None, ["PF00264"], "single", "catalytic", True, "asp_melanin", "",
  "Tyrosinase of the Asp-melanin cluster. UniProt Q0CRX0 gives it EC 1.14.18.- (not EC 1.14.18.1) and two reactions with aspulvinone E (RHEA:74195, RHEA:74199). The abstract of PMID 27133313 states that TyrP activates aspulvinone E for polymerisation "
  "and states no assay, so the links are held at T3. The products have charged ChEBI entries only, so the route gives no product in the tables. TyrP is a short tyrosinase-like protein (356 residues): 45% identical to ivoB, 36% to orsC and 35% to the "
  "catechol oxidase AoCO4 (decoys), and 13 to 17% to the references of step tyrosinase (MEASURED 2026-11-21). Whether TyrP oxidises L-tyrosine or L-DOPA is not stated in the records that were read. "
  "Gene names: UniProt gives tyrP the synonym melB; melB of Aspergillus oryzae is another, long tyrosinase (step tyrosinase). " + ORG + " Described for Aspergillus terreus only; absence elsewhere is not a defect.",
  "type-3 copper protein family (short tyrosinase-like)"),
]
NOFLAG = "exogenous_substrate: false. "
ROUTES = [
 ("tyrosinase_dopa", "tyrosinase route to L-dopaquinone (Neurospora crassa T, Aspergillus oryzae melB, Agaricus bisporus PPO4)", ["tyrosinase"],
  NOFLAG + "One step. Tyrosinase makes L-dopaquinone from L-tyrosine (RHEA:18117) or from L-DOPA (RHEA:34287). From L-dopaquinone onward (L-dopachrome, 5,6-dihydroxyindole-2-carboxylic acid, the polymer) the table is non-enzymatic; "
  "a fungal enzyme there is T4 (step dopachrome_onward, in kb/melanin_ldopa/step.tsv only). The polymer is not a compound of the table."),
 ("laccase_exogenous", "laccase route from exogenous diphenols (Cryptococcus neoformans and Cryptococcus deneoformans LAC1, LAC2)", ["laccase_diphenol"],
  FLAG + ". One step, and its substrate comes from outside the fungus: Cryptococcus neoformans produces melanin only from exogenous substrates and needs diphenolic precursors (PMID 17098915). "
  "The laccase oxidises the diphenol that is supplied: L-DOPA to L-dopaquinone, dopamine to decarboxy dopachrome (PMID 8300520), homogentisic acid to a pigment (PMID 17098915). There is no tyrosine hydroxylation. "
  "A capability call for this route needs the substrate in the environment; gene presence alone shows nothing about the pigment."),
 ("pyomelanin", "pyomelanin route through homogentisic acid (Aspergillus fumigatus hppD)", ["hppd"],
  NOFLAG + "One enzyme step: HppD makes homogentisic acid. The pigment forms by oxidation of homogentisic acid with no enzyme link in the table (curated reaction hga_autooxidation). "
  "HmgA is not in this route: it removes homogentisic acid (route tyrosine_degradation). Reference organism: Aspergillus fumigatus (PMID 19028908)."),
 ("tyrosine_degradation", "removal of homogentisic acid, competing with pyomelanin (Aspergillus fumigatus hppD, hmgA)", ["hppd", "hmga"],
  NOFLAG + "This is not a pigment route. HppD then HmgA: homogentisic acid is removed as 4-maleylacetoacetic acid. It competes with route pyomelanin; the hmgA deletion of Aspergillus fumigatus raises pigment formation (PMID 19028908). "
  "The route is declared so that hmgA counts against pyomelanin and not for it."),
 ("asp_melanin", "Asp-melanin route (Aspergillus terreus melA, tyrP)", ["aspulvinone_e_synthase", "asp_melanin_tyrosinase"],
  NOFLAG + "MelA makes aspulvinone E; TyrP activates it for polymerisation (Aspergillus terreus, PMID 27133313). The TyrP products have charged ChEBI entries only (CHEBI:193114, CHEBI:193115), "
  "so the tables give no product for this route; the first step alone gives aspulvinone E."),
]
step_rx = {}
for r in srows:
    step_rx.setdefault(r[0], [])
    if r[1] not in step_rx[r[0]]:
        step_rx[r[0]].append(r[1])
def y(s):  # YAML double-quoted scalar
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'
out = [textwrap.dedent('''\
    # L-DOPA melanin and other melanin routes, step set v2 (stage 1, task S5, 2026-11-21).
    #
    # Status: DRAFT for user review. Pathway-level verified: false.
    # Do not run a census against this file until the user has reviewed the step list.
    #
    # Reference checks: refs/melanin_ldopa/provenance.tsv, refs/melanin_ldopa/decoys/provenance.tsv and
    # refs/melanin_ldopa/UNVERIFIED.md. Evidence tiers, quotes and reactions: kb/melanin_ldopa/*.tsv.
    #
    # A ref has verified: true only when, on 2026-11-21, (1) the UniProt entry was fetched and its
    # organism and gene name matched, (2) the UniProt entry carries the EC and Rhea reaction of the
    # step, and (3) the PMID exists in PubMed and is cited by the entry. It never means a person
    # reviewed it. Q2TL94 is verified: false because UniProt gives the entry no gene name.
    # A decoy has verified: true when the UniProt entry was fetched on 2026-11-21, its organism and
    # gene name matched, it carries the family ids of the step, and it carries neither the EC nor a
    # Rhea reaction of the step. The decoys are reviewed entries, except Q2UNF9 (AoCO4), whose
    # function is in two abstracts. A decoy pmid is given only when the abstract was read.
    #
    # Field notes:
    # - Routes: five routes. They are alternatives; no route joins two of them. Tyrosinase and
    #   laccase are separate steps in separate routes. A laccase is never linked to tyrosine
    #   hydroxylation.
    # - exogenous_substrate: the schema has no field for it. Each route note starts with
    #   "exogenous_substrate: true" or "exogenous_substrate: false". Only laccase_exogenous is true:
    #   Cryptococcus neoformans makes melanin only from exogenous substrates (PMID 17098915). The
    #   lineage_note of step laccase_diphenol carries the same text. Tests check both.
    # - Dopaquinone to dopachrome and onward is non-enzymatic in the tables. The T4 step
    #   dopachrome_onward is in kb/melanin_ldopa/step.tsv only (it has no reference protein).
    # - species: the organism name of the UniProt entry, copied exactly. UniProt names serotype D
    #   Cryptococcus neoformans "Cryptococcus deneoformans", Pholiota nameko "Pholiota microspora"
    #   and Pycnoporus sanguineus "Trametes sanguinea".
    # - gene_name: a gene name of the UniProt entry (an ORF name where the entry has no other
    #   name). It is valid only together with species.
    # - pmid of a ref: a paper cited by the UniProt entry. The papers that show each reaction are
    #   in step_reaction.tsv; some of them are not cited by the entry (UNVERIFIED.md, section A).
    # - domains_any: Pfam / NCBIfam ids from the UniProt cross-references of the references. The
    #   families are broader than the steps. A step call needs HMM plus tree placement with decoys.
    # - reactions: ids of kb/melanin_ldopa/reaction.tsv. The tier of each (protein, reaction) link
    #   is in kb/melanin_ldopa/step_reaction.tsv, not here.
    # - copy_expect, optional and role are PROPOSED values. No census data exists yet. Every step
    #   is optional: the routes are alternatives, and absence of one is not a defect.
    # - KEGG: K00505, K00421, K05909, K00457 and K00451 are cross-references in
    #   kb/melanin_ldopa/step.tsv only (names read from KEGG list on 2026-11-21).
    id: melanin_ldopa
    verified: false
    routes:
    ''')]
for rid, name, sids, note in ROUTES:
    out.append(f"  - id: {rid}\n    name: {y(name)}\n    step_ids: [{', '.join(sids)}]\n" + (f"    note: {y(note)}\n" if note else ""))
out.append("steps:\n")
krows = []
for sid, name, ec, doms, copy, role, opt, route, ko, note, family in STEPS:
    for acc, *_ in REF[sid]:
        assert set(doms) <= xref_ids(acc), (sid, acc, xref_ids(acc))
    o = [f"  - id: {sid}", f"    name: {y(name)}", f"    ec: {y(ec) if ec else 'null'}"]
    if route:
        o.append(f"    route: {route}")
    o.append("    refs:")
    for acc, gene, pmid, ver in REF[sid]:
        o += [f"      # {c}" for c in REF_NOTE.get((sid, acc), [])]
        o += [f"      - accession: {acc}", f"        species: {y(species(acc))}", f"        gene_name: {y(gene)}",
              f"        pmid: {y(pmid)}", f"        verified: {str(ver).lower()}"]
    if DECOY.get(sid):
        o.append("    decoys:")
        for acc, gene, pmid, *_ in DECOY[sid]:
            o += [f"      - accession: {acc}", f"        species: {y(species(acc))}", f"        gene_name: {y(gene)}",
                  f"        pmid: {y(pmid) if pmid else 'null'}", "        verified: true"]
    o += [f"    domains_any: [{', '.join(doms)}]", f"    copy_expect: {copy}", f"    role: {role}",
          "    fused_with: []", f"    optional: {str(opt).lower()}",
          "    reactions:"] + [f"      - {y(r)}" for r in step_rx[sid]] + [f"    lineage_note: {y(note)}", ""]
    out.append("\n".join(o) + "\n")
    krows.append([sid, route or "", ko, ";".join(doms), family, str(opt).lower()])
krows.append([T4_STEP, "", "", "", "unknown (no fungal enzyme identified)", "true"])
(ROOT / "pathways" / f"{PW}.yaml").write_text("".join(out).rstrip("\n") + "\n", encoding="utf-8")
tsv(KB / "step.tsv", COLUMNS["step"], krows)
n = {t: sum(1 for r in srows if r[5].startswith(t)) for t in ("T1", "T2", "T3", "T4")}
print("refs", len(prov), "decoy rows", len(dprov), "distinct decoys", len({r[1] for r in dprov}),
      "compounds", len(crows), "reactions", len(RX), "links", len(srows), n)
