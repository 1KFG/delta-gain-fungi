"""Build the DHN-melanin curation files from checked database records (cached fetches only).

Writes pathways/melanin_dhn.yaml, refs/melanin_dhn/ (FASTA, provenance, decoys),
kb/melanin_dhn/*.tsv and tests/data/kb/melanin_dhn_quotes.tsv. It does not write UNVERIFIED.md.
Every curated fact below is checked against the fetched record by an assert.
"""
from h import *
import csv, re, textwrap
from pathway_census.kb import formula_mass, COLUMNS

ROOT = PKG
PW = "melanin_dhn"
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

def xref_ids(acc):
    return {x["id"] for x in full(acc).get("uniProtKBCrossReferences", []) if x["database"] in ("Pfam", "NCBIfam")}

def wrap(seq):
    return "\n".join(seq[i:i+60] for i in range(0, len(seq), 60))

def norm(s):
    return " ".join(s.split())

def quote(pmid, start, end):
    """Cut one sentence out of the fetched abstract. The text is never typed by hand."""
    ab = pubmed_abstract(pmid)["abstract"]
    i = ab.index(start); j = ab.index(end, i) + len(end)
    q = ab[i:j]
    assert q.endswith(".") and len(q) <= 400, (pmid, q)
    assert norm(q) in norm(ab)
    assert not re.search(r"[.!?] [A-Z]", q[:-1]), (pmid, q)
    return q

Q = {
 "10501004": quote("10501004", "The production of T4HN by the PKS1 transformant", "pentaketide T4HN."),
 "23261597": quote("23261597", "Here, in vitro characterization of Pks1", "hexaketide substrate."),
 "11179356": quote("11179356", "The resulting wdpks1Delta disruptants", "wild-type parent."),
 "11040426": quote("11040426", "The Alb1p PKS product in this heterologous expression system", "pentaketide T4HN."),
 "21277986": quote("21277986", "The fawn phenotype was complemented by a PKS protein", "in A. fumigatus."),
 "15310761": quote("15310761", "The Ayg1p converted the naphthopyrone YWA1", "acetoacetic acid."),
 "11350964": quote("11350964", "Cell-free enzyme assays, using the crude Ayg1p protein extract", "to 1,3,6,8-THN."),
 "18676950": quote("18676950", "Chemical analyses of mutant cultures demonstrated", "in the culture media."),
 "10956664": quote("10956664", "Substrate competition experiments indicate", "by a factor of 4.2."),
 "10515939": quote("10515939", "The presence of a DHN-melanin pathway in A. fumigatus was supported", "in the arp2 deletants."),
 # second sentence of PMID 10515939, for the abr1 link (one sentence per link, fix round 1)
 "10515939:abr1": quote("10515939", "Involvement of the six genes in conidial pigmentation", "disruption of each gene in A. fumigatus."),
 "28517364": quote("28517364", "Combining genetic deletion and combinatorial gene expression", "synthesis of scytalone."),
 "8672814": quote("8672814", "A melanin deficient mutant 9141 (Thr-) has a defect", "melanin biosynthetic pathway."),
 "12967213": quote("12967213", "However melanin production by the double mutant was restored by THR1", "in vivo."),
 "8953707": quote("8953707", "A plasmid containing the SCD1 gene transformed", "C. lagenarium 8015."),
 "9466792": quote("9466792", "Scytalone dehydratase was eight-fold more specific", "than it was for scytalone."),
 "10320327": quote("10320327", "Relative substrate specificities for the two physiological substrates", "within the active site."),
 "16988817": quote("16988817", "Deletion of the abr2 gene changed the gray-green conidial pigment", "wild-type conidia."),
}

# ---------- references: step -> [(accession, gene_name, reference pmid, verified)] ----------
# The reference pmid is a paper tied to the entry by UniProt (cited by the entry), or, where
# marked in REF_NOTE, by the abstract (locus tag or protein length). Papers that show the
# reaction are in the step_reaction links.
NO_GENE = {"Q9HFV6": "4HNR"}  # UniProt entry without any gene name: protein abbreviation of the paper
REF = {
 "thn_synthase": [("P79068", "PKS1", "7500937", True)],
 "hexaketide_pks": [("Q9Y7A7", "PKS1", "11179356", True)],
 "heptaketide_pks": [("Q4WZA8", "alb1", "11040426", True), ("A2QUI2", "albA", "21277986", True)],
 "ayg1": [("Q4WZB3", "ayg1", "15310761", True), ("Q6BCB5", "yg1", "18676950", True), ("A2R3S8", "An14g05350", "21277986", True)],
 "t4hn_reductase": [("Q9HFV6", "4HNR", "10956664", False), ("E9QUT3", "arp2", "10515939", True), ("W3XAA0", "PfmaG", "28517364", True),
                    ("F8WTI4", "T4HR1", "12967213", True), ("A7DWM4", "teh", "17681008", True)],
 "t3hn_reductase": [("Q12634", "MGG_02252", "8112349", True), ("P87025", "THR1", "8672814", True), ("A7DWM6", "tih", "17681008", True)],
 "scytalone_dehydratase": [("Q00455", "SCD1", "8953707", True), ("O14434", "arp1", "9383199", True), ("P56221", "SDH1", "9571787", True),
                           ("A7DWM5", "sdh", "17681008", True)],
 "laccase_abr2": [("E9RBR0", "abr2", "16988817", True)],
 "mco_abr1": [("Q4WZB4", "abr1", "10515939", True), ("A2R3T0", "An14g05370", "21277986", True)],
}
# references whose pmid is NOT cited by the UniProt entry: how the paper is tied to the entry
ABSTRACT_TIE = {
 ("heptaketide_pks", "A2QUI2"): ("An09g05730", "FwnA"),   # locus tag in the abstract and in UniProt
 ("ayg1", "A2R3S8"): ("An14g05350", "OlvA"),
 ("mco_abr1", "A2R3T0"): ("An14g05370", "BrnA"),
}
REF_NOTE = {
 ("thn_synthase", "P79068"): [
   "UniProt organism name: Glomerella lagenarium (UniProt taxid 5462; NCBI scientific name Colletotrichum",
   "lagenaria). The papers say Colletotrichum lagenarium. pmid 7500937 is the paper cited by the entry",
   "(PKS1, 2187 residues; the entry has 2187 residues). Product papers: PMID 10501004, 23261597 (links)."],
 ("hexaketide_pks", "Q9Y7A7"): [
   "UniProt organism name: Exophiala dermatitidis. The papers say Wangiella dermatitidis (WdPKS1)."],
 ("heptaketide_pks", "A2QUI2"): [
   "FwnA of PMID 21277986: the abstract gives the locus tag An09g05730, a gene name of this entry.",
   "UniProt names the gene albA and does not cite PMID 21277986."],
 ("ayg1", "Q6BCB5"): ["WdYG1 of PMID 18676950 (UniProt gene name: yg1)."],
 ("ayg1", "A2R3S8"): [
   "OlvA of PMID 21277986: the abstract gives the locus tag An14g05350. The entry has this ORF name only",
   "and does not cite the paper."],
 ("t4hn_reductase", "Q9HFV6"): [
   "verified: false, for one reason only: the UniProt entry has NO gene name, so the gene-name rule cannot",
   "be met. gene_name below is NOT a UniProt gene name; it is the protein abbreviation 4HNR of the paper",
   "(PMID 11467929). The identity is not in doubt: the entry cites PMID 10956664 and 11467929 (PDB 1JA9)",
   "and is 45.3% identical to the 3HNR Q12634. UniProt organism: Pyricularia grisea (taxid 148305); the",
   "papers say Magnaporthe grisea. The user decides at the review gate (UNVERIFIED.md)."],
 ("t4hn_reductase", "F8WTI4"): [
   "The entry cites PMID 12967213 for the sequence. The abstract deduces a T4HN-specific reductase and",
   "does not name T4HR1. Link tier T3."],
 ("t3hn_reductase", "Q12634"): [
   "The UniProt entry has an ORF name only (MGG_02252). UniProt names it tetrahydroxynaphthalene",
   "reductase (EC 1.1.1.252). It is the trihydroxynaphthalene reductase (3HNR) of PMID 10956664: it cites",
   "PMID 8112349 and the 3HNR crystal structures (PMID 8939741, 11342131), and it is 45% identical to",
   "Q9HFV6 (the abstract of PMID 10956664 gives 46% between 4HNR and 3HNR)."],
 ("scytalone_dehydratase", "P56221"): [
   "Enzyme assay papers cited by the entry: PMID 9466792 (scytalone), 10320327 (vermelone)."],
 ("mco_abr1", "A2R3T0"): [
   "BrnA of PMID 21277986: the abstract gives the locus tag An14g05370. The entry has this ORF name only,",
   "does not cite the paper, and carries an automatic EC 1.16.3.1 (ferroxidase)."],
}

# ---------- decoys: step -> [(accession, gene_name, pmid or None, reason)] ----------
PKS_FAM = "NR-PKS (NCBIfam TIGR04532, Pfam PF16073) member; UniProt reviewed entry: "
D_WA = ("Q03149", "wA", "11251292", PKS_FAM + "heptaketide naphthopyrone YWA1 synthase (RHEA:62652) of the conidial yellow pigment; "
        "not a T4HN or 2-acetyl-T4HN synthase. It is NOT a decoy of heptaketide_pks: it does that reference reaction")
D_PKS = [
 ("Q12053", "aflC", None, PKS_FAM + "norsolorinic acid synthase of aflatoxin biosynthesis (EC 2.3.1.221, RHEA:35179)"),
 ("Q12397", "stcA", None, PKS_FAM + "norsolorinic acid synthase of sterigmatocystin biosynthesis (EC 2.3.1.221, RHEA:35179)"),
 ("Q6DQW3", "CTB1", None, PKS_FAM + "nor-toralactone synthase of cercosporin biosynthesis (RHEA:62892)"),
 ("S0DZM7", "bik1", None, PKS_FAM + "bikaverin polyketide synthase (pre-bikaverin; no Rhea reaction on the entry)"),
]
SDR = "Pfam PF13561 (short-chain dehydrogenase) member; UniProt reviewed entry: "
D_RED = [
 ("P50161", "aflM", None, SDR + "versicolorin reductase of aflatoxin biosynthesis (EC 1.3.1.-, RHEA:64296); no EC 1.1.1.252, no RHEA:21908"),
 ("Q00791", "stcU", None, SDR + "versicolorin reductase of sterigmatocystin biosynthesis (EC 1.3.1.-, RHEA:64296); no EC 1.1.1.252, no RHEA:21908"),
 ("Q8TFD5", "ver1", None, SDR + "versicolorin reductase of dothistromin biosynthesis (EC 1.1.-.-; no Rhea reaction on the entry)"),
 ("P87219", "SOU1", None, SDR + "sorbose reductase (EC 1.1.1.289, RHEA:14609)"),
]
MCO = "multicopper oxidase (Pfam PF07732, PF00394, PF07731) member; UniProt reviewed entry: "
D_YA = ("P17489", "yA", "2659435", MCO + "conidial laccase (EC 1.10.3.2, RHEA:11276) that converts the yellow pigment made by wA to the green conidial pigment; not in a DHN-melanin cluster")
D_LAC = [
 ("Q55P57", "LAC1", None, MCO + "laccase (EC 1.10.3.2, RHEA:11276) that oxidises L-dopa and other diphenols to melanin precursors (laccase melanin, not DHN-melanin)"),
 ("P06811", "lacc", None, MCO + "laccase (EC 1.10.3.2, RHEA:11276); lignin degradation and detoxification"),
]
D_FET3 = ("P38993", "FET3", None, MCO + "iron transport ferroxidase (EC 1.16.3.1, RHEA:11148; EC 1.16.3.4, RHEA:30083)")
DECOY = {
 "thn_synthase": [D_WA] + D_PKS,
 "hexaketide_pks": [D_WA] + D_PKS,
 "heptaketide_pks": D_PKS,
 "t4hn_reductase": D_RED,
 "t3hn_reductase": D_RED,
 "laccase_abr2": [D_YA] + D_LAC + [("Q70KY3", "LAC1", None, MCO + "laccase (EC 1.10.3.2, RHEA:11276); lignin degradation and detoxification"), D_FET3],
 "mco_abr1": [D_FET3, D_YA] + D_LAC,
}
DECOY_XREF = {"thn_synthase": {"TIGR04532", "PF16073"}, "hexaketide_pks": {"TIGR04532", "PF16073"}, "heptaketide_pks": {"TIGR04532", "PF16073"},
              "t4hn_reductase": {"PF13561"}, "t3hn_reductase": {"PF13561"},
              "laccase_abr2": {"PF07732", "PF00394", "PF07731"}, "mco_abr1": {"PF07732", "PF00394", "PF07731"}}
# Rhea reactions a decoy of the step must not carry (the reference reactions; for the two curated
# PKS reactions also the bacterial T4HN synthase reaction RHEA:41524)
FORBIDDEN = {"thn_synthase": {"RHEA:41524"}, "hexaketide_pks": {"RHEA:41524"}, "heptaketide_pks": {"RHEA:62652"},
             "t4hn_reductase": {"RHEA:21908"}, "t3hn_reductase": {"RHEA:21908"},
             "laccase_abr2": set(), "mco_abr1": set()}

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
        ab = pubmed_abstract(pmid)
        assert ab["title"]
        if (step, acc) in ABSTRACT_TIE:
            locus, name = ABSTRACT_TIE[(step, acc)]
            assert pmid not in e["pmids"] and locus in e["gene_names"] and f"{locus} ({name})" in ab["abstract"].replace("(FwnA, An09g05730)", "An09g05730 (FwnA)"), (step, acc)
            src = srcdb(acc) + f"; PubMed (abstract gives the locus tag {locus} for {name}; the entry does not cite the paper)"
        else:
            assert pmid in e["pmids"], (step, acc, pmid, e["pmids"])
            src = srcdb(acc) + "; PubMed (cited by the entry)"
        recs.append(f">{header(acc)}\n{wrap(e['sequence'])}\n")
        prov.append([step, acc, e["organism"], gene, "", pmid, src, TODAY, str(ver).lower()])
    (REFS / f"{step}.faa").write_text("".join(recs), encoding="utf-8")
assert uniprot_entry("P79068")["length"] == 2187 and "2187 amino acids" in pubmed_abstract("7500937")["abstract"]
tsv(REFS / "provenance.tsv", ["step_id", "accession", "species", "gene_name", "region", "pmid", "source_database", "retrieved_on", "verified"], prov)

# ---------- decoys ----------
ALL_REF = {a for refs in REF.values() for a, *_ in refs}
dprov = []
for step, ds in DECOY.items():
    recs = []
    for acc, gene, pmid, reason in ds:
        assert acc not in ALL_REF, acc
        e = check_entry(acc, gene)
        assert e["reviewed"], acc
        assert DECOY_XREF[step] <= xref_ids(acc), (step, acc, xref_ids(acc))
        assert not (rheas(acc) & FORBIDDEN[step]), (step, acc, rheas(acc))
        for rid in re.findall(r"RHEA:\d+", reason.split("It is NOT")[0].split("; no ")[0]):
            assert rid in rheas(acc), (acc, rid)
        if "no Rhea reaction" in reason:
            assert not rheas(acc), acc
        for ec in re.findall(r"EC (\d+\.\d+\.[\d-]+\.[\d-]+)", reason.split("; no EC")[0]):
            assert ec in e["ec"], (acc, ec, e["ec"])
        assert "1.1.1.252" not in e["ec"], acc
        if pmid:
            assert pmid in e["pmids"] and pubmed_abstract(pmid)["title"]
        recs.append(f">{header(acc)}\n{wrap(e['sequence'])}\n")
        dprov.append([step, acc, e["organism"], gene, pmid or "", reason, srcdb(acc) + ("; PubMed" if pmid else ""), TODAY, "true"])
    (DEC / f"{step}.faa").write_text("".join(recs), encoding="utf-8")
tsv(DEC / "provenance.tsv", ["step_id", "accession", "species", "gene_name", "pmid", "reason", "source_database", "retrieved_on", "verified"], dprov)

# ---------- compounds (ChEBI id -> KEGG cross-reference found in the ChEBI record) ----------
KEGG = {"15351": "C00024", "15531": "C00083", "133763": "", "15344": "C00164", "18365": "C04033", "16945": "C00779",
        "18393": "C01173", "9957": "C01624", "167604": "", "42646": "C18012"}
# ChEBI compound whose InChIKey is not in PubChem: ChEBI number -> (PubChem name, PubChem InChIKey)
PUBCHEM_MISMATCH = {"42646": ("flaviolin", "XNPCAGMCQDGQKK-UHFFFAOYSA-N")}
crows = []
for num, kegg in KEGG.items():
    c = chebi_compound(num)
    assert c["charge"] in (0, "0"), (num, c["charge"])
    m = formula_mass(c["formula"])
    if num in PUBCHEM_MISMATCH:
        # PubChem has no compound with the ChEBI InChIKey; the PubChem compound of that name is
        # a tautomer with another first block. Recorded in UNVERIFIED.md. Formula and mass agree.
        try:
            pubchem_by_inchikey(c["inchikey"]); raise AssertionError(f"{num}: PubChem now has the ChEBI key")
        except KbFetchError as e:
            assert "404" in str(e), e
        name, key = PUBCHEM_MISMATCH[num]
        pn = raw("pubchem_name", name, "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + name +
                 "/property/InChIKey,MolecularFormula,ExactMass,IUPACName,ConnectivitySMILES/JSON")["PropertyTable"]["Properties"][0]
        assert pn["InChIKey"] == key and key[:14] != c["inchikey"][:14] and pn["MolecularFormula"] == c["formula"], pn
        assert abs(m - float(pn["ExactMass"])) < 0.0005
    else:
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
# charged Rhea participants -> neutral parents (formula check: neutral = charged + n H)
for charged, neutral, nh in [("57288", "15351", 4), ("57384", "15531", 5), ("13705", "15344", 1)]:
    a, b = chebi_compound(charged), chebi_compound(neutral)
    ha = int(re.search(r"H(\d+)", a["formula"]).group(1)); hb = int(re.search(r"H(\d+)", b["formula"]).group(1))
    assert hb - ha == nh and re.sub(r"H\d+", "", a["formula"]) == re.sub(r"H\d+", "", b["formula"]), (charged, neutral)
    assert int(a["charge"]) == -nh and a["inchikey"][:14] == b["inchikey"][:14], (charged, a["charge"])

# ---------- reactions ----------
C = "CUR:melanin_dhn:"
SM = "Small molecules (H2O, CO2, CoA, NADP(H), H+) are not listed as participants."
NO_AT4HN = ("ChEBI has no entry for 2-acetyl-1,3,6,8-tetrahydroxynaphthalene (name searches and an InChIKey search on 2026-11-21 found none; "
            "PubChem CID 86018477, C12H10O5, InChIKey FYPSLERIHMOYRN-UHFFFAOYSA-N), so it is not a participant.")
RX = [  # id, ec, note, substrates, products
 (C + "t4hn_synthesis", "", "Curated: a true THN synthase releases 1,3,6,8-tetrahydroxynaphthalene (T4HN) directly (PMID 10501004, 23261597). The starter and extender units are not stated in the abstracts that were read, so no substrate is listed. The bacterial type III PKS reaction of Rhea (EC 2.3.1.233, KEGG rppA) is a different enzyme class and is not used.", [], [("18365", 1)]),
 (C + "at4hn_synthesis", "", f"Curated: a hexaketide PKS releases 2-acetyl-1,3,6,8-tetrahydroxynaphthalene (AT4HN) (PMID 18676950). {NO_AT4HN} The abstract shows AT4HN in the WdYG1 deletion mutant; it does not show the product of the PKS directly.", [], []),
 ("RHEA:62652", "", f"Rhea has no EC for this reaction. Rhea lists the charged species CHEBI:57288 and CHEBI:57384; mapped to the neutral parents. {SM}", [("15351", 1), ("15531", 6)], [("133763", 1)]),
 ("RHEA:62700", "", f"Chain shortening of the heptaketide YWA1 by the hydrolase Ayg1. Rhea has no EC for this reaction. Rhea lists acetoacetate CHEBI:13705 (1-); mapped to the neutral parent CHEBI:15344. {SM}", [("133763", 1)], [("18365", 1), ("15344", 1)]),
 (C + "at4hn_deacetylation", "", f"Curated: deacetylation of the hexaketide AT4HN to T4HN by the Ayg1 orthologue Yg1 (PMID 18676950). {NO_AT4HN} With no substrate in the table, path_possible cannot reach this reaction after another step.", [], [("18365", 1)]),
 ("RHEA:21908", "1.1.1.252", f"Rhea writes the reaction in the oxidation direction (scytalone + NADP+ = T4HN + NADPH + H+). The roles are stored in the reduction direction that the step uses. KEGG KO K17739 is a cross-reference of the step. {SM}", [("18365", 1)], [("16945", 1)]),
 (C + "t3hn_reduction", "", f"Curated: NADPH-dependent reduction of 1,3,8-trihydroxynaphthalene to vermelone. A Rhea search for 'vermelone' on 2026-11-21 returned no reaction. ChEBI:9957 is vermelone without stereochemistry (2 stars); ChEBI also has (+)-vermelone CHEBI:198683. {SM}", [("18393", 1)], [("9957", 1)]),
 ("RHEA:24396", "4.2.1.94", f"KEGG KO K17740 is a cross-reference of the step. {SM}", [("16945", 1)], [("18393", 1)]),
 (C + "vermelone_dehydration", "", f"Curated: dehydration of vermelone to 1,8-dihydroxynaphthalene. Rhea has no reaction for vermelone (search on 2026-11-21). This reaction has its own evidence records; they are not copies of the scytalone edge. Which enzyme acts is not settled for every fungus: the abstract of PMID 10320327 names scytalone and vermelone as the two physiological substrates of scytalone dehydratase (Pyricularia oryzae enzyme, UniProt P56221); for the Aspergillus fumigatus cluster the abstracts that were read do not show which enzyme dehydrates vermelone (pathway review, O4). {SM}", [("9957", 1)], [("167604", 1)]),
 (C + "dhn_polymerisation", "", "Curated: oxidative polymerisation of 1,8-dihydroxynaphthalene to DHN-melanin. The polymer has no formula (ChEBI:167605, 1,8-dihydroxynaphthalene-melanin, is a class entry), so no product is listed. The abstracts show colour phenotypes only: the substrate of each oxidase and the order of the laccase and the multicopper oxidase are not shown.", [("167604", 1)], []),
 (C + "t4hn_autooxidation", "", "Curated shunt, no enzyme link: oxidation of T4HN to flaviolin. Flaviolin accumulates when the reduction of T4HN is blocked (arp2 deletants, PMID 10515939). That the oxidation is non-enzymatic is from the pathway review and was not checked in an abstract.", [("18365", 1)], [("42646", 1)]),
]
for rid, ec, *_ in RX:  # re-check each Rhea id and EC against the cached Rhea response
    if rid.startswith("RHEA:"):
        t = raw("rhea", "id_" + rid[5:], "https://www.rhea-db.org/rhea/?query=" + f"rhea%3A{rid[5:]}" + "&columns=rhea-id,equation,chebi-id,ec&format=tsv&limit=3", as_json=False)
        line = [l for l in t.splitlines() if l.startswith(rid + "\t")]
        assert line and line[0].split("\t")[3] == (f"EC:{ec}" if ec else ""), (rid, t)
assert raw("rhea", "vermelone", "https://www.rhea-db.org/rhea/?query=vermelone&columns=rhea-id,equation,chebi-id,ec&format=tsv&limit=10", as_json=False).strip().count("\n") == 0
tsv(KB / "reaction.tsv", COLUMNS["reaction"],
    [[rid, ec, "", "forward", "false" if rid.startswith("RHEA:") else "true", note, rid if rid.startswith("RHEA:") else ""] for rid, ec, note, s, p in RX])
tsv(KB / "reaction_participant.tsv", COLUMNS["reaction_participant"],
    [[rid, f"CHEBI:{c}", role, str(n)] for rid, ec, note, s, p in RX for role, lst in (("substrate", s), ("product", p)) for c, n in lst])

# ---------- step_reaction links ----------
L = []  # step, reaction, accession, tier, assay, pmid, with_quote
def link(step, rid, acc, tier, assay, pmid, q=True):
    L.append((step, rid, acc, tier, assay, pmid, q))
NIGER = "complementation of the colour mutant and targeted disruption in the fungus; colour phenotype. The abstract does not name the product of the enzyme"
SORD = "homology only: the abstract of PMID 17681008 reports isolation of the gene and its homology to DHN-melanin genes"
link("thn_synthase", C + "t4hn_synthesis", "P79068", "T1a", "heterologous expression in Aspergillus oryzae (host named in the first sentence of the abstract); the PKS1 transformant produces T4HN, isolated and identified as its tetraacetate", "10501004")
link("thn_synthase", C + "t4hn_synthesis", "P79068", "T1a", "in vitro characterisation of Pks1; the thioesterase domain cyclises and deacetylates an enzyme-bound hexaketide", "23261597")
link("hexaketide_pks", C + "at4hn_synthesis", "Q9Y7A7", "T1c", "gene disruption in the fungus; albino phenotype. The abstract does not show the product of WdPks1p; the hexaketide AT4HN is shown in the WdYG1 deletion mutant (PMID 18676950)", "11179356")
link("heptaketide_pks", "RHEA:62652", "Q4WZA8", "T1a", "heterologous expression in Aspergillus oryzae; product identified as the heptaketide naphthopyrone", "11040426")
link("heptaketide_pks", "RHEA:62652", "A2QUI2", "T1c", NIGER + " (FwnA, An09g05730). The same abstract reports that naphtho-gamma-pyrones depend on fwnA; YWA1 is not named", "21277986")
link("ayg1", "RHEA:62700", "Q4WZB3", "T1a", "enzyme purified from an Aspergillus oryzae overexpression strain; YWA1 converted to T4HN with release of acetoacetic acid", "15310761")
link("ayg1", "RHEA:62700", "Q4WZB3", "T1a", "cell-free assay with crude Ayg1p extract; heptaketide product of Alb1p shortened to T4HN", "11350964")
link("ayg1", "RHEA:62700", "A2R3S8", "T1c", NIGER + " (OlvA, An14g05350)", "21277986")
link("ayg1", C + "at4hn_deacetylation", "Q6BCB5", "T1b", "gene replacement in the fungus; the mutants accumulate AT4HN and 3-acetylflaviolin (chemical analysis of cultures)", "18676950")
COMP = "enzyme produced in Escherichia coli and purified; substrate competition between T4HN and T3HN"
link("t4hn_reductase", "RHEA:21908", "Q9HFV6", "T1a", COMP + ". The second reductase (4HNR) prefers T4HN by a factor of 310", "10956664")
link("t4hn_reductase", C + "t3hn_reduction", "Q9HFV6", "T1a", COMP + ". Minor activity: T3HN is accepted with a 310-fold lower preference than T4HN", "10956664")
link("t4hn_reductase", "RHEA:21908", "E9QUT3", "T1b", "gene deletion in the fungus; arp2 deletants accumulate flaviolin only, the arp1 deletant accumulates scytalone and flaviolin", "10515939")
link("t4hn_reductase", C + "t3hn_reduction", "E9QUT3", "T3", "Inferred. The abstracts that were read (PMID 10515939, 11350964) do not show the reduction of T3HN by Arp2. The six-gene cluster has one reductase; a second reductase outside the cluster was not checked.", "", False)
link("t4hn_reductase", "RHEA:21908", "W3XAA0", "T1a", "heterologous expression in Aspergillus nidulans; the PKS PfmaE and PfmaG are sufficient for scytalone, whose structure was elucidated", "28517364")
link("t4hn_reductase", "RHEA:21908", "F8WTI4", "T3", "The abstract of PMID 12967213 deduces a T4HN-specific reductase from a double mutant; it does not name T4HR1 and shows no assay of this protein. The UniProt entry cites the paper for the sequence.", "12967213", False)
link("t4hn_reductase", "RHEA:21908", "A7DWM4", "T3", SORD + " (teh).", "17681008", False)
link("t3hn_reductase", C + "t3hn_reduction", "Q12634", "T1a", "substrate competition between T3HN and T4HN; trihydroxynaphthalene reductase (3HNR) prefers T3HN by a factor of 4.2. The abstract does not name an accession; Q12634 is tied to 3HNR by the papers it cites (PMID 8112349, 8939741, 11342131) and by 45% identity to Q9HFV6", "10956664")
link("t3hn_reductase", "RHEA:21908", "Q12634", "T1a", "substrate competition between T3HN and T4HN; 3HNR accepts T4HN (preference for T3HN only 4.2-fold). Side activity of the 3HNR type", "10956664")
link("t3hn_reductase", C + "t3hn_reduction", "P87025", "T1b", "the Thr- mutant 9141 has a defect in the conversion of T3HN to vermelone and is complemented by the cosmid that carries THR1; gene disruption gives a light brown phenotype", "8672814")
link("t3hn_reductase", "RHEA:21908", "P87025", "T1b", "THR1 restores melanin production in a double mutant that metabolises neither T4HN nor T3HN. Side activity of the 3HNR type", "12967213")
link("t3hn_reductase", C + "t3hn_reduction", "A7DWM6", "T3", SORD + " (tih).", "17681008", False)
link("scytalone_dehydratase", "RHEA:24396", "Q00455", "T1c", "complementation of the melanin-deficient Scd- mutant 9201Y to the wild phenotype. The abstract reports no intermediate analysis", "8953707")
link("scytalone_dehydratase", "RHEA:24396", "O14434", "T1b", "gene deletion in the fungus; the arp1 deletant accumulates scytalone and flaviolin", "10515939")
link("scytalone_dehydratase", "RHEA:24396", "P56221", "T1a", "continuous multiwavelength enzyme assay with scytalone and a synthetic substrate (enzyme named as scytalone dehydratase of Magnaporthe grisea; UniProt P56221 cites the paper for catalytic activity)", "9466792")
link("scytalone_dehydratase", "RHEA:24396", "A7DWM5", "T3", SORD + " (sdh). The same abstract reports RNA interference of sdh; it states no phenotype and no intermediate.", "17681008", False)
link("scytalone_dehydratase", C + "vermelone_dehydration", "Q00455", "T3", "Inferred. UniProt Q00455 states dehydration of vermelone and cites PMID 8953707; the abstract of that paper does not show it.", "8953707", False)
link("scytalone_dehydratase", C + "vermelone_dehydration", "O14434", "T3", "Inferred. Which enzyme dehydrates vermelone in Aspergillus fumigatus is not settled by the abstracts that were read (pathway review, O4).", "", False)
link("scytalone_dehydratase", C + "vermelone_dehydration", "P56221", "T3", "Held at T3 (rule 1, strict reading, harmonisation round): the quoted sentence states steady-state kinetics with the two physiological substrates scytalone and vermelone (wild type and Ser129 mutant); it names no enzyme, no organism and no product. Only the first sentence of the abstract says scytalone dehydratase. UniProt P56221 cites this paper for catalytic activity and for the mutagenesis of the residues named in the abstract. Under the context reading the link is T1a.", "10320327")
link("scytalone_dehydratase", C + "vermelone_dehydration", "A7DWM5", "T3", SORD + " (sdh).", "17681008", False)
link("laccase_abr2", C + "dhn_polymerisation", "E9RBR0", "T1c", "gene deletion in the fungus; conidial colour changes from gray-green to brown. The abstract does not show the substrate of the laccase", "16988817")
link("mco_abr1", C + "dhn_polymerisation", "Q4WZB4", "T3", "Held at T3. The quoted sentence reports altered conidial colour after disruption of each of the six cluster genes; it does not name abr1 (the abstract names abr1 as one of the six genes in other sentences). Raise to T1c if the user accepts the sentence for abr1. The abstract does not show the substrate of the oxidase.", "10515939", "10515939:abr1")
link("mco_abr1", C + "dhn_polymerisation", "A2R3T0", "T1c", NIGER + " (BrnA, An14g05370)", "21277986")

acc_gene = {(s, a): g for s, refs in REF.items() for a, g, *_ in refs}
p56221_pm = uniprot_entry("P56221")["pmids"]
assert "9466792" in p56221_pm and "10320327" in p56221_pm
assert "disruption of each gene in A. fumigatus" in pubmed_abstract("10515939")["abstract"]
srows = []
for step, rid, acc, tier, assay, pmid, q in L:
    assert rid in {r[0] for r in RX} and (step, acc) in acc_gene, (step, rid, acc)
    if pmid:
        assert pubmed_abstract(pmid)["title"]
    qt = Q[q] if isinstance(q, str) else (Q[pmid] if q else "")  # q: True = sentence of the paper, str = key of another sentence
    assert qt or tier not in ("T1a", "T1b", "T1c", "T2"), (step, acc, tier)
    assert not qt or norm(qt) in norm(pubmed_abstract(pmid)["abstract"]), (step, acc, pmid)
    srows.append([step, rid, species(acc), acc_gene[(step, acc)], acc, tier, assay, pmid, qt])
assert set(Q.values()) == {r[8] for r in srows if r[8]}
tsv(KB / "step_reaction.tsv", COLUMNS["step_reaction"], srows)
# change-detector fixture for the tests (all links)
tsv(ROOT / "tests" / "data" / "kb" / "melanin_dhn_quotes.tsv", ["step_id", "reaction_id", "protein_accession", "tier", "assay", "pmid", "quote"],
    [[r[0], r[1], r[4], r[5], r[6], r[7], r[8]] for r in srows])

# ---------- steps ----------
PKS_NOTE = ("Three PKS routes give 1,3,6,8-tetrahydroxynaphthalene (T4HN): a true THN synthase (this route needs no Ayg1), a hexaketide PKS plus Yg1, and a "
            "heptaketide PKS plus Ayg1. All are non-reducing iterative type I PKS; the product chain length cannot be read from a domain hit, so the call needs "
            "tree placement with the references and decoys of all three PKS steps. ")
SAME = ("The same gene name PKS1 is used for two enzymes: PKS1 of Colletotrichum lagenarium is a THN synthase, PKS1 (WdPKS1) of Exophiala (Wangiella) "
        "dermatitidis is the hexaketide synthase. A gene name is valid only with its organism. ")
LINEAGE = ("Lineage restrictions: Aspergillus terreus lacks the genes for DHN-melanin and makes Asp-melanin (MelA, TyrP; PMID 27133313), so absence there is not a defect. "
           "Botrytis cinerea and other Leotiomycetes carry two melanin PKS (bcpks12 for sclerotia, bcpks13 for conidia; PMID 26514268); BcYGH1 converts the BcPKS13 product.")
PARA = ("The two naphthol reductase types are paralogues of one short-chain dehydrogenase family with real substrate overlap (PMID 10956664, 12967213). "
        "One paralogue present does not mean the other reaction is absent, and one paralogue absent does not mean the step is absent: ")
STEPS = [  # id, name, ec, domains, copy, role, optional, route, ko, lineage_note, family
 ("thn_synthase", "1,3,6,8-tetrahydroxynaphthalene synthase, non-reducing PKS (PKS1 of Colletotrichum lagenarium)", None, ["TIGR04532", "PF16073"], "single", "catalytic", True, "pks_thn", "",
  PKS_NOTE + SAME + "The UniProt organism name of P79068 is Glomerella lagenarium (UniProt taxid 5462); the other Colletotrichum references are Colletotrichum orbiculare (NCBI taxid 5465). The two names are not merged. " + LINEAGE,
  "non-reducing iterative type I polyketide synthase"),
 ("hexaketide_pks", "hexaketide synthase, non-reducing PKS (WdPKS1 of Exophiala dermatitidis)", None, ["TIGR04532", "PF16073"], "single", "catalytic", True, "pks_hexaketide", "",
  PKS_NOTE + SAME + "The product 2-acetyl-1,3,6,8-tetrahydroxynaphthalene is shown in the WdYG1 deletion mutant (PMID 18676950), not for the PKS itself. " + LINEAGE,
  "non-reducing iterative type I polyketide synthase"),
 ("heptaketide_pks", "heptaketide naphthopyrone YWA1 synthase, non-reducing PKS (alb1/pksP of Aspergillus fumigatus)", None, ["TIGR04532", "PF16073"], "single", "catalytic", True, "pks_heptaketide", "",
  PKS_NOTE + "Aspergillus nidulans wA (Q03149) does the same reaction (YWA1) for a pigment that is not DHN-melanin, so a hit in this clade does not show DHN-melanin without Ayg1 and the downstream steps. "
  "The gene is named alb1 or pksP in Aspergillus fumigatus and albA or fwnA in Aspergillus niger; UniProt ties albA of Aspergillus niger to aurasperone (naphthopyrone) biosynthesis. " + LINEAGE,
  "non-reducing iterative type I polyketide synthase"),
 ("ayg1", "polyketide chain-shortening hydrolase (ayg1 of Aspergillus fumigatus, yg1 of Exophiala dermatitidis)", None, ["PF06500"], "single", "catalytic", True, None, "",
  "Needed only where the PKS releases a heptaketide (YWA1, Aspergillus fumigatus) or a hexaketide (2-acetyl-T4HN, Exophiala dermatitidis). Not needed with a true THN synthase (Colletotrichum lagenarium PKS1), so absence on route pks_thn is not a defect. "
  "Botrytis cinerea BcYGH1 converts the BcPKS13 product (PMID 26514268; no UniProt entry was tied to it). The two reference reactions differ: Ayg1 releases acetoacetic acid from YWA1, Yg1 removes an acetyl group.",
  "FrsA-like alpha/beta hydrolase family"),
 ("t4hn_reductase", "tetrahydroxynaphthalene reductase, 4HNR/Arp2 type (T4HN to scytalone)", "1.1.1.252", ["PF13561"], "single", "catalytic", False, None, "K17739",
  PARA + "the 3HNR type also reduces T4HN (Pyricularia oryzae 3HNR, Colletotrichum THR1), and the Aspergillus fumigatus cluster has one reductase (arp2). "
  "Versicolorin reductases of aflatoxin-type pathways are in the same family (decoys). MEASURED on 2026-11-21 (global alignment): 4HNR Q9HFV6 is 45% identical to 3HNR Q12634; Arp2 E9QUT3 is 52% identical to Q9HFV6 and 44% to Q12634. "
  "KEGG has one KO (K17739) for tetrahydroxynaphthalene reductase; which paralogue type its members belong to was not checked.",
  "short-chain dehydrogenase/reductase family (naphthol reductase, 4HNR type)"),
 ("t3hn_reductase", "trihydroxynaphthalene reductase, 3HNR/THR1/BUF1/BRM2 type (T3HN to vermelone)", None, ["PF13561"], "single", "catalytic", False, None, "",
  PARA + "the 3HNR type reduces T4HN as well (4.2-fold preference for T3HN in Pyricularia oryzae), and the 4HNR type accepts T3HN weakly (310-fold preference for T4HN). "
  "UniProt names Q12634 tetrahydroxynaphthalene reductase with EC 1.1.1.252; it is the 3HNR of the papers. The names BUF1 and BRM2 are not gene names of a reference here (no UniProt entry was tied to them). "
  "MEASURED on 2026-11-21: the versicolorin reductase decoys are 59 to 63% identical to the 3HNR references, more than the 4HNR references are (43 to 50%).",
  "short-chain dehydrogenase/reductase family (naphthol reductase, 3HNR type)"),
 ("scytalone_dehydratase", "scytalone dehydratase (SCD1, arp1, SDH1)", "4.2.1.94", ["PF02982"], "single", "catalytic", False, None, "K17740",
  "Two reactions with separate evidence: scytalone to 1,3,8-trihydroxynaphthalene (RHEA:24396) and vermelone to 1,8-dihydroxynaphthalene (curated). The step is listed twice in each PKS route for this reason. "
  "The abstract of PMID 10320327 reports kinetics of scytalone dehydratase with scytalone and vermelone, in a sentence that names no enzyme, so the vermelone link of the Pyricularia oryzae enzyme is held at T3 under the strict reading of rule 1 (T1a under the context reading); for the other references the reaction is inferred (T3), and for Aspergillus fumigatus the acting enzyme is not settled. "
  "The gene names differ by organism: SCD1 (Colletotrichum), arp1 (Aspergillus fumigatus), SDH1 (Pyricularia oryzae), sdh (Sordaria macrospora). sdh and SDH1 also name succinate dehydrogenase genes in other organisms (not checked here).",
  "scytalone dehydratase family"),
 ("laccase_abr2", "laccase of the DHN-melanin cluster, Abr2 type (abr2 of Aspergillus fumigatus)", None, ["PF07732", "PF00394", "PF07731"], "single", "catalytic", True, "polymerisation", "",
  "Evidence is a colour phenotype only (T1c); the substrate and the order of the two oxidases are not shown. Laccases are a large family with many functions; the nearest characterised relative of Abr2 is the Aspergillus nidulans conidial laccase yA "
  "(39% identity, MEASURED 2026-11-21), which acts on the wA pigment and is a decoy. One reference only. Described for the Aspergillus fumigatus cluster; the range in other lineages is not checked, so absence is not a defect.",
  "multicopper oxidase family (laccase, Abr2 type)"),
 ("mco_abr1", "multicopper oxidase of the DHN-melanin cluster, Abr1 type (abr1 of Aspergillus fumigatus, brnA of Aspergillus niger)", None, ["PF07732", "PF00394", "PF07731"], "single", "catalytic", True, "polymerisation", "",
  "Evidence is a colour phenotype only: T1c for BrnA of Aspergillus niger; the Abr1 link is held at T3 because the quoted sentence of PMID 10515939 does not name abr1. The substrate and the order of the two oxidases are not shown. Abr1 and BrnA are 43% identical (MEASURED 2026-11-21) and 28 to 31% identical to the iron transport ferroxidase FET3 (decoy); "
  "the unreviewed Aspergillus niger entry carries an automatic ferroxidase EC. Described for Aspergillus fumigatus and Aspergillus niger; the range in other lineages is not checked, so absence is not a defect.",
  "multicopper oxidase family (Fet3-like, Abr1 type)"),
]
TAIL = ["t4hn_reductase", "scytalone_dehydratase", "t3hn_reductase", "scytalone_dehydratase"]
TAILNOTE = " Then T4HN to scytalone, 1,3,8-trihydroxynaphthalene, vermelone and 1,8-dihydroxynaphthalene. scytalone_dehydratase is listed twice (two dehydrations). The route ends at the monomer; the oxidases are in route polymerisation."
ROUTES = [
 ("pks_thn", "(a) THN synthase route (Colletotrichum lagenarium PKS1) to 1,8-dihydroxynaphthalene", ["thn_synthase"] + TAIL,
  "The PKS releases T4HN directly. No Ayg1/Yg1 step." + TAILNOTE),
 ("pks_hexaketide", "(b) hexaketide PKS plus Yg1 route (Exophiala dermatitidis WdPKS1, WdYg1) to 1,8-dihydroxynaphthalene", ["hexaketide_pks", "ayg1"] + TAIL,
  "The PKS releases 2-acetyl-T4HN; Yg1 deacetylates it to T4HN." + TAILNOTE + " The hexaketide has no ChEBI entry, so the tables cannot follow this route past the PKS step."),
 ("pks_heptaketide", "(c) heptaketide PKS plus Ayg1 route (Aspergillus fumigatus Alb1, Ayg1) to 1,8-dihydroxynaphthalene", ["heptaketide_pks", "ayg1"] + TAIL,
  "The PKS releases the naphthopyrone YWA1; Ayg1 shortens it to T4HN." + TAILNOTE),
 ("polymerisation", "oxidases of the polymerisation of 1,8-dihydroxynaphthalene (Aspergillus fumigatus Abr2 and Abr1)", ["laccase_abr2", "mco_abr1"],
  "Not a chain: two oxidase families with colour-phenotype evidence only (T1c; the Abr1 link is held at T3). The order is not settled and the polymer is not a compound of the table. It follows each of the three PKS routes."),
]
step_rx = {}
for r in srows:
    step_rx.setdefault(r[0], [])
    if r[1] not in step_rx[r[0]]:
        step_rx[r[0]].append(r[1])
def y(s):  # YAML double-quoted scalar
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'
out = [textwrap.dedent('''\
    # DHN-melanin pathway definition, step set v2 (stage 1, task S4, 2026-11-21).
    #
    # Status: DRAFT for user review. Pathway-level verified: false.
    # Do not run a census against this file until the user has reviewed the step list.
    #
    # Reference checks: refs/melanin_dhn/provenance.tsv, refs/melanin_dhn/decoys/provenance.tsv and
    # refs/melanin_dhn/UNVERIFIED.md. Evidence tiers, quotes and reactions: kb/melanin_dhn/*.tsv.
    #
    # A ref has verified: true only when, on 2026-11-21, (1) the UniProt entry was fetched and its
    # organism and gene name matched, (2) the UniProt entry or the PubMed abstract supports the
    # reaction of the step, and (3) the PMID exists in PubMed. It never means a person reviewed it.
    # Q9HFV6 is verified: false because UniProt gives the entry no gene name (UNVERIFIED.md).
    # A decoy has verified: true when the reviewed UniProt entry was fetched on 2026-11-21, its
    # organism and gene name matched, it carries the family ids of the step, it does not carry the
    # Rhea reaction of the step, and its record states a function that is not the function of the
    # step. A decoy pmid is given only when the abstract was read.
    #
    # Field notes:
    # - Routes: the three PKS routes (a), (b), (c) differ in the first one or two steps and share
    #   the steps from T4HN to 1,8-dihydroxynaphthalene. The Ayg1/Yg1 step is in routes (b) and
    #   (c) only. The PKS steps are alternatives; they are never consecutive steps of one route.
    # - species: the organism name of the UniProt entry, copied exactly. UniProt names
    #   Colletotrichum lagenarium PKS1 "Glomerella lagenarium", the other Colletotrichum
    #   references "Colletotrichum orbiculare", Wangiella dermatitidis "Exophiala dermatitidis",
    #   and the rice blast fungus "Pyricularia oryzae" (70-15) or "Pyricularia grisea".
    # - gene_name: a gene name of the UniProt entry (an ORF name where the entry has no other
    #   name). It is valid only together with species.
    # - pmid of a ref: a paper tied to the entry (cited by the UniProt entry, or the abstract
    #   gives the locus tag). The papers that show each reaction are in step_reaction.tsv.
    # - domains_any: Pfam / NCBIfam ids from the UniProt cross-references of the references. The
    #   families are broader than the steps. A step call needs HMM plus tree placement with decoys.
    # - reactions: ids of kb/melanin_dhn/reaction.tsv. The tier of each (protein, reaction) link is
    #   in kb/melanin_dhn/step_reaction.tsv, not here.
    # - copy_expect, optional and role are PROPOSED values. No census data exists yet.
    # - KEGG: K17739 and K17740 are cross-references in kb/melanin_dhn/step.tsv only. K19580
    #   (rppA) is a bacterial enzyme and is not used.
    id: melanin_dhn
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
        for acc, gene, pmid, reason in DECOY[sid]:
            o += [f"      - accession: {acc}", f"        species: {y(species(acc))}", f"        gene_name: {y(gene)}",
                  f"        pmid: {y(pmid) if pmid else 'null'}", "        verified: true"]
    o += [f"    domains_any: [{', '.join(doms)}]", f"    copy_expect: {copy}", f"    role: {role}",
          "    fused_with: []", f"    optional: {str(opt).lower()}",
          "    reactions:"] + [f"      - {y(r)}" for r in step_rx[sid]] + [f"    lineage_note: {y(note)}", ""]
    out.append("\n".join(o) + "\n")
    krows.append([sid, route or "", ko, ";".join(doms), family, str(opt).lower()])
(ROOT / "pathways" / f"{PW}.yaml").write_text("".join(out).rstrip("\n") + "\n", encoding="utf-8")
tsv(KB / "step.tsv", COLUMNS["step"], krows)
n = {t: sum(1 for r in srows if r[5].startswith(t)) for t in ("T1", "T2", "T3", "T4")}
print("refs", len(prov), "decoy rows", len(dprov), "distinct decoys", len({r[1] for r in dprov}),
      "compounds", len(crows), "reactions", len(RX), "links", len(srows), n)
