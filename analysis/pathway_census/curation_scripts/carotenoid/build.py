"""Build the carotenoid curation files from checked database records (cached fetches only)."""
from h import *
import csv, io, re, textwrap
from pathlib import Path
from pathway_census.kb import formula_mass, COLUMNS

ROOT = PKG
REFS = ROOT / "refs" / "carotenoid"; DEC = REFS / "decoys"; KB = ROOT / "kb" / "carotenoid"
TODAY = "2026-11-21"

def full(acc):
    return raw("uniprot_full", acc, f"https://rest.uniprot.org/uniprotkb/{acc}.json")

def header(acc):
    d = full(acc)
    db = "sp" if d["entryType"].startswith("UniProtKB reviewed") else "tr"
    pdn = d["proteinDescription"]
    nm = (pdn.get("recommendedName") or (pdn.get("submissionNames") or [{}])[0])["fullName"]["value"]
    gn = d.get("genes", [{}])[0]
    g = (gn.get("geneName") or (gn.get("orfNames") or [{}])[0] or {}).get("value")
    h = f"{acc} {db}|{acc}|{d['uniProtkbId']} {nm} OS={d['organism']['scientificName']} OX={d['organism']['taxonId']}"
    return h + (f" GN={g}" if g else "")

def srcdb(acc):
    d = full(acc); au = d["entryAudit"]
    kind = "Swiss-Prot" if d["entryType"].startswith("UniProtKB reviewed") else "TrEMBL"
    return f"UniProtKB/{kind} (entry v{au['entryVersion']}, sequence v{au['sequenceVersion']})"

def wrap(seq):
    return "\n".join(seq[i:i+60] for i in range(0, len(seq), 60))

def norm(s):
    return " ".join(s.split())

def quote(pmid, start, end):
    ab = pubmed_abstract(pmid)["abstract"]
    i = ab.index(start); j = ab.index(end, i) + len(end)
    q = ab[i:j]
    assert q.endswith("."), (pmid, q)
    assert norm(q) in norm(ab)
    # one sentence: no sentence break inside (". " followed by a capital), except known abbreviations
    inner = q[:-1].replace("E. coli", "E_coli").replace("M. circinelloides", "M_c").replace("X. dendrorhous", "X_d").replace("3, 4-", "3,4-")
    assert not re.search(r"\. [A-Z]", inner), (pmid, q)
    return q

Q = {
 "12695851": quote("12695851", "Functional analyses in Escherichia coli showed", "synthase activity."),
 "16928467": quote("16928467", "The enzyme comprised the two catalytic activities", "lycopene cyclase."),
 "12172798": quote("12172798", "The resulting mutants were albino", "carotene cyclase."),
 "11172012": quote("11172012", "The sequences of the wild type and six mutants", "mutant phenotypes."),
 "10951210": quote("10951210", "Using functional analyses in Escherichia coli and Mucor circinelloides", "other than fungi."),
 "10589832": quote("10589832", "Expression studies in different carotenogenic", "to beta-carotene."),
 "11862485": quote("11862485", "The mutants lack neurosporaxanthin", "and HPLC."),
 "11017770": quote("11017770", "The major reaction products were", "and lycopene."),
 "9079885": quote("9079885", "Heterologous expression of the Phycomyces gene", "phytoene to lycopene."),
 "10397832": quote("10397832", "We conclude that in X. dendrorhous", "single gene product."),
 "19645721": quote("19645721", "HPLC analyses indicated a specific impairment", "four-step desaturation."),
 "17293483": quote("17293483", "The expression of the enzyme in beta-carotene-accumulating", "chromophore retinal."),
 "19584000": quote("19584000", "Heterologous overexpression and in vitro analyses", "two molecules of retinal."),
 "17493127": quote("17493127", "The purified CarT enzyme efficiently cleaved torulene", "synthetic substrates."),
 "17610084": quote("17610084", "The biological function was confirmed by in vitro assays", "aldehyde of neurosporaxanthin."),
 "18627463": quote("18627463", "In vitro incubation of purified YLO-1 protein", "carboxyl group."),
 "21749649": quote("21749649", "To confirm this function in vivo", "expression of carD."),
 "25504221": quote("25504221", "Introduction of the crtS on autonomously replicating plasmid", "was verified."),
 "16416328": quote("16416328", "Functional characterization of the enzyme by genetic modification", "carbon 3 and 4."),
 "18837978": quote("18837978", "CBSTr forms yellow colonies", "accumulating beta-carotene."),
 "21854466": quote("21854466", "According to these and other results", "by AcaA."),
 "22445959": quote("22445959", "Phycomyces mutants in gene carS accumulate", "sexual process."),
 "17822929": quote("17822929", "The functionality of the tsp3 gene product", "tsp3 expression."),
 "18931040": quote("18931040", "The NADP-dependent 4-dihydrotrisporin-dehydrogenase is", "to trisporic acid."),
 "8828234": quote("8828234", "We have purified the NADP-dependent", "Mucor mucedo."),
}

# ---------- references: step -> [(accession, gene_name, pmid, verified)] ----------
FUSED = [("P37295", "al-2"), ("Q8X0Z1", "carRA"), ("Q9P854", "carRA"), ("Q9UUQ6", "CARRP"), ("Q7Z859", "crtYB")]
PSY_PMID = {"P37295": "8163509", "Q8X0Z1": "12172798", "Q9P854": "11172012", "Q9UUQ6": "10951210", "Q7Z859": "10589832"}
CYC_PMID = dict(PSY_PMID, P37295="11862485")
REF = {
 "ggpp_synthase": [("P24322", "al-3", "1826006", True), ("Q9P885", "carG", "12695851", True)],
 "phytoene_synthase": [(a, g, PSY_PMID[a], True) for a, g in FUSED],
 "lycopene_cyclase": [(a, g, CYC_PMID[a], True) for a, g in FUSED],
 "phytoene_desaturase": [("P21334", "al-1", "11017770", True), ("Q8X0Z0", "carB", "12172798", True),
                         ("P54982", "carB", "9079885", True), ("O13506", "crtI", "10397832", True)],
 "carx": [("S0ENQ2", "carX", "17293483", True), (None, "cco1", "19584000", False)],
 "cart": [("A1KQY3", "carT", "17493127", True), ("A7UXI1", "cao-2", "17610084", True)],
 "ylo1": [("Q1K615", "ylo-1", "18627463", True), ("F6IBC7", "carD", "21749649", True)],
 "astaxanthin_synthase": [("Q3HR17", "crtS", "25504221", True), ("B5LEC8", "crtS", "18769767", True)],
 "crtr": [("B6D970", "crtR", "18837978", True)],
 "cars": [("G8Z9M3", "carS", "22445959", True)],
 "acaa": [("A0ABR3AZC2", "acaA", "21854466", False)],
 "tsp3": [("A5HXP5", "tsp3", "17822929", True)],
 "tsp2": [("B9DR41", "tsp2", "18931040", True)],
 "tdh": [("Q01213", "tdh", "8828234", True)],
}
REF_NOTE = {
 ("astaxanthin_synthase", "Q3HR17"): [
   "pmid 25504221: cited by the UniProt entry; expression in Mucor circinelloides (T1a for ketolation).",
   "Functional paper in X. dendrorhous: PMID 16416328. UniProt does not cite it for Q3HR17. The link rests",
   "on sequence identity (557 of 557) with NCBI CAC16261.1, the patent sequence AX034665.1, and on shared",
   "authors (Hoshino, Ojima, Setoguchi on the patent record)."],
 ("cars", "G8Z9M3"): [
   "pmid 22445959: cited by GenBank ADU04395.1, which is identical to G8Z9M3 (628 of 628). UniProt cites no",
   "paper for G8Z9M3. The product (beta-apo-12'-carotenal) is shown in PMID 21854466."],
 ("acaa", "A0ABR3AZC2"): [
   "verified: false. The gene name acaA on this entry comes from a 2024 genome annotation (EMBL",
   "JBCLYO010000008, KAL0086242.1). No record ties the sequence to the protein assayed in PMID 21854466.",
   "GenBank OAD75744.1 (PHYBLDRAFT_180852, strain NRRL 1555(-)) has the identical sequence (602 of 602).",
   "The link in step_reaction.tsv is held at T3 until the protein id in the paper is compared."],
 ("carx", None): [
   "verified: false. No UniProt accession could be tied to Ustilago maydis cco1 (see UNVERIFIED.md A1)."],
}
NULL_SPECIES = {"cco1": "Ustilago maydis"}

CCO = [("Q7S860", "cao-1", "23893079", "PF03055 member; did not convert carotenoids; cleaves resveratrol and piceatannol (abstract of PMID 23893079)"),
       ("A0A0D1E6L2", "RCO1", None, "PF03055 member; UniProt reviewed entry: resveratrol cleavage oxygenase 1 (EC 1.13.11.-)"),
       ("A0A384JIP2", "rco1", None, "PF03055 member; UniProt reviewed entry: resveratrol cleavage oxygenase 1 (EC 1.13.11.-)"),
       ("Q4W9H1", "rco1", None, "PF03055 member; UniProt reviewed entry: resveratrol cleavage oxygenase 1 (EC 1.13.11.-)"),
       ("G2QIL8", "MYCTH_2309496", None, "PF03055 member; UniProt reviewed entry: aromatic dioxygenase (EC 1.13.11.88, 1.13.11.94); the gene has an ordered locus name only"),
       ("Q52008", "lsdB", None, "PF03055 member (bacterial); UniProt reviewed entry: lignostilbene-alpha,beta-dioxygenase isozyme III (EC 1.13.11.43)")]
P450FOXY = ("Q9Y8G7", "CYP505", None, "PF00067 member; UniProt reviewed entry: bifunctional cytochrome P450/NADPH--P450 reductase (fatty acid hydroxylase P450foxy); the P450 domain is fused to a reductase domain (PF00667, PF00258, PF00175). Decoy of astaxanthin_synthase only: its reductase domain does the reaction of crtr, so it is not a crtr decoy")
DECOY = {
 "ggpp_synthase": [("P08524", "ERG20", None, "PF00348 member; UniProt reviewed entry: farnesyl pyrophosphate synthase (EC 2.5.1.10, 2.5.1.1), no EC 2.5.1.29"),
                   ("O14230", "fps1", None, "PF00348 member; UniProt reviewed entry: farnesyl pyrophosphate synthase (EC 2.5.1.10, 2.5.1.1), no EC 2.5.1.29"),
                   ("A0A1D8PH78", "ERG20", None, "PF00348 member; UniProt reviewed entry: farnesyl pyrophosphate synthase (EC 2.5.1.10, 2.5.1.1), no EC 2.5.1.29"),
                   ("P18900", "COQ1", None, "PF00348 member; UniProt reviewed entry: hexaprenyl pyrophosphate synthase, mitochondrial (EC 2.5.1.-)")],
 "phytoene_synthase": [(a, g, None, "PF00494 member; UniProt reviewed entry: squalene synthase (EC 2.5.1.21)") for a, g in
                       [("P29704", "ERG9"), ("P36596", "erg9"), ("Q92459", "ERG9"), ("Q4WAG4", "erg9")]],
 "phytoene_desaturase": [("P50264", "FMS1", None, "PF01593 member; UniProt reviewed entry: polyamine oxidase (EC 1.5.3.17)"),
                         ("P46882", "maoN", None, "PF01593 member; UniProt reviewed entry: monoamine oxidase N (EC 1.4.3.4)"),
                         ("Q10062", "hem14", None, "PF01593 member; UniProt reviewed entry: protoporphyrinogen oxidase (EC 1.3.3.4)"),
                         ("P23623", "lox", None, "PF01593 member; UniProt reviewed entry: L-amino-acid oxidase (EC 1.4.3.2)"),
                         ("B8NI10", "cpaO", None, "PF01593 member; UniProt reviewed entry: beta-cyclopiazonate dehydrogenase (EC 1.21.99.1)")],
 "carx": CCO, "cart": CCO, "cars": CCO, "acaa": CCO, "tsp3": CCO,
 "ylo1": [("P54115", "ALD6", None, "PF00171 member; UniProt reviewed entry: magnesium-activated aldehyde dehydrogenase, cytosolic (EC 1.2.1.4)"),
          ("P46367", "ALD4", None, "PF00171 member; UniProt reviewed entry: potassium-activated aldehyde dehydrogenase, mitochondrial (EC 1.2.1.5)"),
          ("P08157", "aldA", None, "PF00171 member; UniProt reviewed entry: aldehyde dehydrogenase (EC 1.2.1.3)"),
          ("Q04458", "HFD1", None, "PF00171 member; UniProt reviewed entry: fatty aldehyde dehydrogenase (EC 1.2.1.3, 1.2.1.64)"),
          ("P38067", "UGA2", None, "PF00171 member; UniProt reviewed entry: succinate-semialdehyde dehydrogenase [NADP(+)] (EC 1.2.1.16)")],
 "astaxanthin_synthase": [("P10614", "ERG11", None, "PF00067 member; UniProt reviewed entry: lanosterol 14-alpha demethylase CYP51 (EC 1.14.14.154)"),
                          ("P54781", "ERG5", None, "PF00067 member; UniProt reviewed entry: C-22 sterol desaturase CYP61 (EC 1.14.19.41)"),
                          ("P23295", "CYP55A1", None, "PF00067 member; UniProt reviewed entry: nitric oxide reductase P450nor (EC 1.7.1.14)"),
                          ("P21595", "DIT2", None, "PF00067 member; UniProt reviewed entry: cytochrome P450-DIT2, CYP56 (EC 1.14.14.-)"),
                          P450FOXY],
 "crtr": [("Q12181", "TAH18", None, "PF00667/PF00258/PF00175 member; UniProt reviewed entry: NADPH-dependent diflavin oxidoreductase 1 (EC 1.18.1.-), not a cytochrome P450 reductase"),
          ("O94613", "tah18", None, "PF00667/PF00258/PF00175 member; UniProt reviewed entry: NADPH-dependent diflavin oxidoreductase 1 (EC 1.18.1.-), not a cytochrome P450 reductase"),
          ("P39692", "MET10", None, "PF00667/PF00175 member; UniProt reviewed entry: sulfite reductase [NADPH] flavoprotein component (EC 1.8.1.2)")],
 "tsp2": [("P21158", "csgA", None, "PF00106 member (bacterial); the abstract of PMID 18931040 names the Myxococcus xanthus CsgA gene product as the closest relative of TSP2 outside the zygomycetes; UniProt reviewed entry: C-signal")],
 "tdh": [("P38715", "GRE3", None, "PF00248 member; UniProt reviewed entry: NADPH-dependent aldose reductase (EC 1.1.1.21)"),
         ("Q12458", "YPR1", None, "PF00248 member; UniProt reviewed entry: aldehyde reductase (EC 1.1.1.2)"),
         ("P14065", "GCY1", None, "PF00248 member; UniProt reviewed entry: glycerol 2-dehydrogenase (NADP(+)) (EC 1.1.1.156)")],
}
DECOY_PFAM = {"ggpp_synthase": ["PF00348"], "phytoene_synthase": ["PF00494"], "phytoene_desaturase": ["PF01593"],
              "carx": ["PF03055"], "cart": ["PF03055"], "cars": ["PF03055"], "acaa": ["PF03055"], "tsp3": ["PF03055"],
              "ylo1": ["PF00171"], "astaxanthin_synthase": ["PF00067"], "crtr": ["PF00667"], "tsp2": ["PF00106"], "tdh": ["PF00248"]}

def species(acc):
    return uniprot_entry(acc)["organism"]

def check_entry(acc, gene):
    e = uniprot_entry(acc)
    assert e["accession"] == acc, acc
    assert gene in e["gene_names"], (acc, gene, e["gene_names"])
    return e

# ---------- carS: GenBank record that ties G8Z9M3 to PMID 22445959 ----------
_gp = raw("ncbi_protein", "ADU04395_gp", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=ADU04395&rettype=gp&retmode=text", as_json=False)
assert re.sub(r"[^a-z]", "", _gp.split("ORIGIN")[1].split("//")[0]).upper() == uniprot_entry("G8Z9M3")["sequence"]
assert "PUBMED   22445959" in _gp
_gp = raw("ncbi_protein", "CAC16261_gp", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=CAC16261&rettype=gp&retmode=text", as_json=False)
assert re.sub(r"[^a-z]", "", _gp.split("ORIGIN")[1].split("//")[0]).upper() == uniprot_entry("Q3HR17")["sequence"]
assert "AX034665.1" in _gp
_gp = raw("ncbi_protein", "OAD75744_gp", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=OAD75744&rettype=gp&retmode=text", as_json=False)
assert re.sub(r"[^a-z]", "", _gp.split("ORIGIN")[1].split("//")[0]).upper() == uniprot_entry("A0ABR3AZC2")["sequence"]
assert "PHYBLDRAFT_180852" in _gp
assert "activity is found exclusively in the (-) mating type" in pubmed_abstract("18931040")["abstract"]

# ---------- regions of the fused proteins ----------
REGION = {"phytoene_synthase": {}, "lycopene_cyclase": {}}
SUPPORT = {"phytoene_synthase": {}, "lycopene_cyclase": {}}
for acc, gene in FUSED:
    d = full(acc)
    reg = {f["description"]: (f["location"]["start"]["value"], f["location"]["end"]["value"])
           for f in d["features"] if f["type"] == "Region"}
    ip = raw("interpro", acc, f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{acc}?page_size=100")
    loc = {}
    for r in ip["results"]:
        for p in r["proteins"]:
            for l in p["entry_protein_locations"]:
                loc.setdefault(r["metadata"]["accession"], []).extend((x["start"], x["end"]) for x in l["fragments"])
    psy, cyc = reg["Phytoene synthase"], reg["Lycopene beta-cyclase"]
    pf = sorted(loc["PF00494"]); tg = sorted(loc["TIGR03462"])
    assert all(psy[0] <= s and e <= psy[1] for s, e in pf), (acc, psy, pf)
    assert all(cyc[0] <= s and e <= cyc[1] for s, e in tg), (acc, cyc, tg)
    assert cyc[1] < psy[0]
    REGION["phytoene_synthase"][acc] = psy; REGION["lycopene_cyclase"][acc] = cyc
    fmt = lambda xs: ",".join(f"{s}-{e}" for s, e in xs)
    SUPPORT["phytoene_synthase"][acc] = f"region = UniProt feature Region 'Phytoene synthase'; it contains Pfam PF00494 at {fmt(pf)} (InterPro)"
    extra = f"; Pfam PF18916 at {fmt(sorted(loc['PF18916']))}" if "PF18916" in loc else "; no Pfam PF18916 match in InterPro"
    SUPPORT["lycopene_cyclase"][acc] = f"region = UniProt feature Region 'Lycopene beta-cyclase'; it contains NCBIfam TIGR03462 at {fmt(tg)} (InterPro){extra}"

# ---------- write reference FASTA and provenance ----------
REFS.mkdir(parents=True, exist_ok=True); DEC.mkdir(exist_ok=True); KB.mkdir(parents=True, exist_ok=True)
(REFS / "full_length").mkdir(exist_ok=True)
prov = []
for step, refs in REF.items():
    recs = []
    for acc, gene, pmid, ver in refs:
        if acc is None:
            prov_sp = NULL_SPECIES[gene]
            continue
        e = check_entry(acc, gene)
        assert pubmed_abstract(pmid)["title"]
        region = REGION.get(step, {}).get(acc)
        seq = e["sequence"]
        h = header(acc)
        src = srcdb(acc) + "; PubMed"
        if region:
            seq = seq[region[0]-1:region[1]]
            h = h.replace(f"{acc} ", f"{acc} region={region[0]}-{region[1]} ", 1)
            src += "; " + SUPPORT[step][acc]
        recs.append(f">{h}\n{wrap(seq)}\n")
        prov.append([step, acc, e["organism"], gene, f"{region[0]}-{region[1]}" if region else "", pmid, src, TODAY, str(ver).lower()])
    (REFS / f"{step}.faa").write_text("".join(recs), encoding="utf-8")
(REFS / "full_length" / "fused.faa").write_text(
    "".join(f">{header(a)}\n{wrap(uniprot_entry(a)['sequence'])}\n" for a, _ in FUSED), encoding="utf-8")
def tsv(path, cols, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar=None)
        w.writerow(cols)
        for r in rows:
            assert len(r) == len(cols), r
            assert not any("\t" in str(x) or "\n" in str(x) for x in r), r
            w.writerow(r)
tsv(REFS / "provenance.tsv", ["step_id", "accession", "species", "gene_name", "region", "pmid", "source_database", "retrieved_on", "verified"], prov)

# ---------- decoys ----------
dprov = []
for step, ds in DECOY.items():
    recs = []
    ref_acc = {r[0] for r in REF[step]}
    for acc, gene, pmid, reason in ds:
        assert acc not in ref_acc
        e = check_entry(acc, gene)
        d = full(acc)
        pf = {x["id"] for x in d.get("uniProtKBCrossReferences", []) if x["database"] == "Pfam"}
        assert set(DECOY_PFAM[step]) <= pf, (step, acc, pf)
        if pmid:
            assert pmid in e["pmids"] and pubmed_abstract(pmid)["title"]
        recs.append(f">{header(acc)}\n{wrap(e['sequence'])}\n")
        dprov.append([step, acc, e["organism"], gene, pmid or "", reason, srcdb(acc) + ("; PubMed" if pmid else ""), TODAY, "true"])
    (DEC / f"{step}.faa").write_text("".join(recs), encoding="utf-8")
tsv(DEC / "provenance.tsv", ["step_id", "accession", "species", "gene_name", "pmid", "reason", "source_database", "retrieved_on", "verified"], dprov)

# ---------- compounds ----------
KEGG = {"16584": "C00129", "17407": "C00448", "48861": "C00353", "27787": "C05421", "15948": "C05432", "62474": "C15867",
        "27740": "C05435", "17579": "C02094", "9638": "C08613", "17898": "C00376", "53157": "", "15825": "C07330",
        "7540": "C08607", "4746": "C08592", "3362": "C08583", "80215": "C15966", "80216": "C15967", "40968": "C08580",
        "53160": "", "53175": "", "211238": ""}
crows = []
for num, kegg in KEGG.items():
    c = chebi_compound(num)
    assert c["charge"] in (0, "0"), (num, c["charge"])
    p = pubchem_by_inchikey(c["inchikey"])
    assert p["inchikey"][:14] == c["inchikey"][:14], num
    m = formula_mass(c["formula"])
    assert abs(m - c["monoisotopic_mass"]) < 0.0005 and abs(m - p["exact_mass"]) < 0.0005, num
    if kegg:
        r = raw("chebi_raw", num, f"https://www.ebi.ac.uk/chebi/backend/api/public/compound/{num}/")
        xs = {x.get("accession_number") for lst in (r.get("database_accessions") or {}).values() for x in lst}
        assert kegg in xs, (num, kegg)
    crows.append([f"CHEBI:{num}", c["inchikey"], c["inchikey"][:14], c["smiles"], c["formula"], f"{m:.4f}", kegg, ""])
tsv(KB / "compound.tsv", COLUMNS["compound"], crows)
# charged Rhea participants -> neutral parents (formula check: neutral = charged + n H)
for charged, neutral, nh in [("128769", "16584", 3), ("175763", "17407", 3), ("58756", "48861", 3), ("63069", "7540", 1)]:
    a, b = chebi_compound(charged), chebi_compound(neutral)
    ha = int(re.search(r"H(\d+)", a["formula"]).group(1)); hb = int(re.search(r"H(\d+)", b["formula"]).group(1))
    assert hb - ha == nh and re.sub(r"H\d+", "", a["formula"]) == re.sub(r"H\d+", "", b["formula"]), (charged, neutral)
    assert int(a["charge"]) == -nh and a["inchikey"][:14] == b["inchikey"][:14], (charged, a["charge"], a["inchikey"], b["inchikey"])

# ---------- reactions ----------
C = "CUR:carotenoid:"
SM = "Small molecules (O2, H2O, acceptor A/AH2, NAD(P), diphosphate) are not listed as participants."
RX = [  # id, ec, note, substrates, products
 ("RHEA:17653", "2.5.1.29", f"Rhea lists the charged species CHEBI:128769, CHEBI:175763 and CHEBI:58756; mapped to the neutral parents. {SM}", [("16584", 1), ("17407", 1)], [("48861", 1)]),
 ("RHEA:34475", "2.5.1.32", f"Rhea lists CHEBI:58756 (3-); mapped to the neutral parent CHEBI:48861. {SM}", [("48861", 2)], [("27787", 1)]),
 ("RHEA:15585", "1.3.99.31", f"Four desaturations, lycopene-forming. {SM}", [("27787", 1)], [("15948", 1)]),
 ("RHEA:30975", "1.3.99.30", f"Five desaturations, 3,4-didehydrolycopene-forming. {SM}", [("27787", 1)], [("62474", 1)]),
 ("RHEA:32219", "5.5.1.19", "First cyclisation.", [("15948", 1)], [("27740", 1)]),
 ("RHEA:32239", "5.5.1.19", "Second cyclisation.", [("27740", 1)], [("17579", 1)]),
 (C + "didehydrolycopene_to_torulene", "", "Curated: cyclisation of one end of 3,4-didehydrolycopene to torulene. A Rhea search for 'torulene' on 2026-11-21 returned no such reaction. Order from the abstract of PMID 11017770 (gamma-carotene is not a substrate of AL-1).", [("62474", 1)], [("9638", 1)]),
 ("RHEA:32887", "1.13.11.63", SM, [("17579", 1)], [("17898", 2)]),
 ("RHEA:31519", "1.13.11.59", SM, [("9638", 1)], [("53157", 1), ("15825", 1)]),
 ("RHEA:31515", "1.2.1.82", f"Rhea lists neurosporaxanthin(1-) CHEBI:63069; mapped to the neutral parent CHEBI:7540. {SM}", [("53157", 1)], [("7540", 1)]),
 ("RHEA:55660", "1.14.99.63", f"4-ketolation of beta-carotene. The EC is the Rhea annotation and is a cross-reference only; the step keeps ec null. {SM}", [("17579", 1)], [("4746", 1)]),
 ("RHEA:55664", "1.14.99.63", f"4'-ketolation of echinenone. The EC is the Rhea annotation and is a cross-reference only; the step keeps ec null. {SM}", [("4746", 1)], [("3362", 1)]),
 (C + "canthaxanthin_3_hydroxylation", "", f"Curated: 3-hydroxylation of canthaxanthin to phoenicoxanthin. Rhea searches for 'phoenicoxanthin' and 'adonirubin' on 2026-11-21 returned no reaction. {SM}", [("3362", 1)], [("80216", 1)]),
 (C + "phoenicoxanthin_3p_hydroxylation", "", f"Curated: 3'-hydroxylation of phoenicoxanthin to astaxanthin. Not RHEA:56772 (adonixanthin is not on this route). {SM}", [("80216", 1)], [("40968", 1)]),
 (C + "echinenone_3_hydroxylation", "", f"Curated: hydroxylation-first branch of the grid, echinenone to 3-hydroxyechinenone. The branch comes from the pathway review (full text of PMID 16416328); the abstract does not name it. ChEBI has two isomers (CHEBI:80215 3-hydroxyechinenone, CHEBI:80214 3'-hydroxyechinenone); CHEBI:80215 is used and the isomer is not checked. {SM}", [("4746", 1)], [("80215", 1)]),
 (C + "3_hydroxyechinenone_4p_ketolation", "", f"Curated: hydroxylation-first branch of the grid, 3-hydroxyechinenone to phoenicoxanthin. Not in the abstract of PMID 16416328. {SM}", [("80215", 1)], [("80216", 1)]),
 (C + "torulene_4_ketolation", "", f"Curated: torulene to 4-ketotorulene. From the pathway review (full text of PMID 16416328); not in the abstract. ChEBI has no entry found for 4-ketotorulene on 2026-11-21, so the product is not listed. {SM}", [("9638", 1)], []),
 (C + "4_ketotorulene_3_hydroxylation", "", f"Curated: 4-ketotorulene to 3-hydroxy-4-ketotorulene. From the pathway review; not in the abstract. ChEBI has no entry found for either compound on 2026-11-21, so no participant is listed. {SM}", [], []),
 (C + "beta_carotene_11p12p_cleavage", "", f"Curated: cleavage of beta-carotene at C11'-C12' to beta-apo-12'-carotenal (C25). The second fragment is not named in the abstract of PMID 21854466 and is not listed. A Rhea search for 12'-apo-beta-carotenal on 2026-11-21 returned no reaction. {SM}", [("17579", 1)], [("53160", 1)]),
 (C + "apo12p_carotenal_13_14_cleavage", "", f"Curated: cleavage of beta-apo-12'-carotenal at C13-C14 to beta-apo-13-carotenone (C18). The second fragment is not named in the abstract and is not listed. Rhea has 13-apo-beta-carotenone only from 10'- and 8'-apo-beta-carotenal (RHEA:26401, RHEA:36043). {SM}", [("53160", 1)], [("53175", 1)]),
 (C + "4_dihydrotrisporin_oxidation", "", "Curated: NADP-dependent oxidation of 4-dihydrotrisporin (to trisporin; the product is not named in the abstract of PMID 18931040). ChEBI searches for 'trisporin' and 'dihydrotrisporin' on 2026-11-21 found no entry, so no participant is listed. Rhea search for 'trisporin' returned no reaction.", [], []),
 (C + "4_dihydromethyltrisporate_oxidation", "", "Curated: NADP-dependent oxidation of 4-dihydromethyltrisporate (to methyltrisporate; the product is not named in the abstract of PMID 8828234). A ChEBI search for 'methyl trisporate' on 2026-11-21 found no entry, so no participant is listed. Rhea search for 'trisporate' returned no reaction.", [], []),
 (C + "torulene_to_torularhodin", "", f"Curated hypothesis (T4): oxidation of torulene to torularhodin. No enzyme is identified. Rhea search for 'torularhodin' on 2026-11-21 returned no reaction. {SM}", [("9638", 1)], [("211238", 1)]),
]
for rid, ec, *_ in RX:  # re-check each Rhea id and EC against the cached Rhea response
    if rid.startswith("RHEA:"):
        t = raw("rhea", "id_" + rid[5:], "https://www.rhea-db.org/rhea/?query=" + f"rhea%3A{rid[5:]}" + "&columns=rhea-id,equation,chebi-id,ec&format=tsv&limit=3", as_json=False)
        line = [l for l in t.splitlines() if l.startswith(rid + "\t")]
        assert line and line[0].split("\t")[3] == f"EC:{ec}", (rid, t)
tsv(KB / "reaction.tsv", COLUMNS["reaction"],
    [[rid, ec, "", "forward", "false" if rid.startswith("RHEA:") else "true", note, rid if rid.startswith("RHEA:") else ""] for rid, ec, note, s, p in RX])
tsv(KB / "reaction_participant.tsv", COLUMNS["reaction_participant"],
    [[rid, f"CHEBI:{c}", role, str(n)] for rid, ec, note, s, p in RX for role, lst in (("substrate", s), ("product", p)) for c, n in lst])

# ---------- step_reaction links ----------
E_COLI = "heterologous expression in Escherichia coli; product identified"
L = []  # step, reaction, accession, gene, tier, assay, pmid, with_quote
def link(step, rid, acc, tier, assay, pmid, q=True):
    L.append((step, rid, acc, tier, assay, pmid, q))
link("ggpp_synthase", "RHEA:17653", "P24322", "T3", "The abstract of PMID 7929398 reports GGPP synthase activity kept in an al-3 deletion mutant; it does not say how the activity was measured. Tier held at T3 until the full text is read.", "7929398", False)
link("ggpp_synthase", "RHEA:17653", "Q9P885", "T1a", "functional analysis in Escherichia coli: GGPP synthase activity", "12695851")
link("phytoene_synthase", "RHEA:34475", "P37295", "T1a", "al-2 cDNA cloned, expressed and functionally characterised; the sentence names no expression host", "16928467")
link("phytoene_synthase", "RHEA:34475", "Q8X0Z1", "T1b", "gene disruption in the fungus; mutants no longer produce phytoene", "12172798")
link("phytoene_synthase", "RHEA:34475", "Q9P854", "T1c", "mutant sequences correlated with colour phenotypes (white, red)", "11172012")
link("phytoene_synthase", "RHEA:34475", "Q9UUQ6", "T1a", "functional analyses in Escherichia coli and in Mucor circinelloides", "10951210")
link("phytoene_synthase", "RHEA:34475", "Q7Z859", "T1a", "expression in carotenogenic Escherichia coli strains", "10589832")
link("lycopene_cyclase", "RHEA:32219", "P37295", "T1a", "al-2 cDNA cloned, expressed and functionally characterised; the sentence names no expression host; the same abstract states that single cyclisations were preferentially catalysed", "16928467")
link("lycopene_cyclase", "RHEA:32239", "P37295", "T3", "The abstract of PMID 16928467 states that AL-2 preferentially does single cyclisations; formation of beta-carotene by AL-2 is not stated.", "16928467", False)
link("lycopene_cyclase", C + "didehydrolycopene_to_torulene", "P37295", "T1b", "mutants in the cyclase domain of al-2, complemented by wild-type al-2; carotenoids analysed by column chromatography and HPLC; the accumulated carotenoid was not identified", "11862485")
link("lycopene_cyclase", "RHEA:32219", "Q8X0Z1", "T3", "carRA disruptants are albino (no phytoene), so the cyclase function is not shown separately in the abstract", "12172798", False)
link("lycopene_cyclase", "RHEA:32239", "Q8X0Z1", "T3", "carRA disruptants are albino (no phytoene), so the cyclase function is not shown separately in the abstract", "12172798", False)
link("lycopene_cyclase", C + "didehydrolycopene_to_torulene", "Q8X0Z1", "T3", "inferred; not shown in the abstract", "12172798", False)
link("lycopene_cyclase", "RHEA:32219", "Q9P854", "T1c", "mutant sequences correlated with colour phenotypes (white, red)", "11172012")
link("lycopene_cyclase", "RHEA:32239", "Q9P854", "T1c", "mutant sequences correlated with colour phenotypes (white, red)", "11172012")
link("lycopene_cyclase", "RHEA:32219", "Q9UUQ6", "T1a", "functional analyses in Escherichia coli and in Mucor circinelloides: lycopene cyclase activity", "10951210")
link("lycopene_cyclase", "RHEA:32239", "Q9UUQ6", "T3", "The abstract of PMID 10951210 states lycopene cyclase activity; it does not name the product (gamma-carotene or beta-carotene).", "10951210", False)
link("lycopene_cyclase", "RHEA:32219", "Q7Z859", "T1a", "expression in carotenogenic Escherichia coli strains: lycopene to beta-carotene", "10589832")
link("lycopene_cyclase", "RHEA:32239", "Q7Z859", "T1a", "expression in carotenogenic Escherichia coli strains: lycopene to beta-carotene", "10589832")
link("lycopene_cyclase", C + "didehydrolycopene_to_torulene", "Q7Z859", "T3", "inferred; not shown in the abstract", "10589832", False)
link("phytoene_desaturase", "RHEA:30975", "P21334", "T1a", "AL-1 expressed in Escherichia coli; active enzyme isolated; products identified in vitro", "11017770")
link("phytoene_desaturase", "RHEA:15585", "P21334", "T1a", "AL-1 expressed in Escherichia coli; active enzyme isolated; lycopene is one of the two major products", "11017770")
link("phytoene_desaturase", "RHEA:15585", "P54982", "T1a", E_COLI + " (lycopene)", "9079885")
link("phytoene_desaturase", "RHEA:15585", "O13506", "T1a", "crtI cDNA expressed in Escherichia coli with Erwinia uredovora GGPP synthase and phytoene synthase genes; lycopene accumulated", "10397832")
link("phytoene_desaturase", "RHEA:30975", "Q8X0Z0", "T1b", "mutant carB allele (SF21) in the fungus with HPLC analysis; targeted allele replacement; expression of the mutant allele in Escherichia coli", "19645721")
link("carx", "RHEA:32887", "S0ENQ2", "T1a", E_COLI + " (retinal, HPLC and GC-MS); in vitro assay with purified enzyme", "17293483")
link("carx", "RHEA:32887", None, "T1a", "heterologous overexpression; in vitro assay with purified enzyme; deletion and over-expression strains", "19584000")
link("cart", "RHEA:31519", "A1KQY3", "T1a", "in vitro assay with purified enzyme", "17493127")
link("cart", "RHEA:31519", "A7UXI1", "T1a", "in vitro assay with purified enzyme; targeted disruption", "17610084")
link("ylo1", "RHEA:31515", "Q1K615", "T1a", "in vitro assay with purified enzyme", "18627463")
link("ylo1", "RHEA:31515", "F6IBC7", "T1a", "in vitro assay with heterologously expressed protein; expression in an Escherichia coli strain that makes beta-apo-4'-carotenal; targeted disruption", "21749649")
HET = "crtS expressed in Mucor circinelloides (heterologous fungal host); echinenone and canthaxanthin detected"
GEN = "complementation of a beta-carotene accumulating mutant and genetic modification of X. dendrorhous; the abstract does not name the intermediates"
for rid in ("RHEA:55660", "RHEA:55664"):
    link("astaxanthin_synthase", rid, "Q3HR17", "T1a", HET, "25504221")
    link("astaxanthin_synthase", rid, "Q3HR17", "T1b", GEN, "16416328")
for r in ("canthaxanthin_3_hydroxylation", "phoenicoxanthin_3p_hydroxylation"):
    link("astaxanthin_synthase", C + r, "Q3HR17", "T1b", GEN + "; the Mucor host of PMID 25504221 has its own beta-carotene hydroxylase, so that paper does not raise the tier", "16416328")
NOTABS = "Reported by the pathway review from the full text of PMID 16416328. Not in the abstract; the abstract states the order as 4-ketolation first, then 3-hydroxylation. Not checked by the curator; tier held at T3."
for r in ("echinenone_3_hydroxylation", "3_hydroxyechinenone_4p_ketolation", "torulene_4_ketolation", "4_ketotorulene_3_hydroxylation"):
    link("astaxanthin_synthase", C + r, "Q3HR17", "T3", NOTABS, "16416328", False)
link("crtr", "RHEA:55660", "B6D970", "T1b", "crtR deletion in the fungus; beta-carotene accumulates and astaxanthin is absent; CrtR is the electron donor of CrtS, not the catalyst", "18837978")
link("cars", C + "beta_carotene_11p12p_cleavage", "G8Z9M3", "T1a", "carS expressed in beta-carotene-producing Escherichia coli; beta-apo-12'-carotenal formed", "21854466")
link("cars", C + "beta_carotene_11p12p_cleavage", "G8Z9M3", "T1b", "carS mutants of the fungus accumulate beta-carotene", "22445959")
link("acaa", C + "apo12p_carotenal_13_14_cleavage", "A0ABR3AZC2", "T3", "The abstract shows the reaction for Phycomyces AcaA (in vitro assay; joint expression with carS in Escherichia coli gave beta-apo-13-carotenone). No record ties the accession A0ABR3AZC2 to the assayed protein (reference verified: false), so the link is held at T3. Raise to T1a on a full-text check of PMID 21854466: compare the protein id given there with PHYBLDRAFT_180852 (identical sequence to A0ABR3AZC2, record OAD75744).", "21854466")
link("tsp3", C + "beta_carotene_11p12p_cleavage", "A5HXP5", "T3", "Loss of carotene in a carotene-producing Escherichia coli strain; the abstract names no product and does not say which species' tsp3 was expressed. The product is inferred from Phycomyces CarS.", "17822929")
link("tsp2", C + "4_dihydrotrisporin_oxidation", "B9DR41", "T3", "The abstract reports purification of the native enzyme from the (-) mating type, peptide sequencing and cloning of the gene. It also states that enzyme activity is found only in the (-) mating type. It names no assay and no product. Tier held at T3; raise to T1 on a full-text check.", "18931040")
link("tdh", C + "4_dihydromethyltrisporate_oxidation", "Q01213", "T3", "The abstract reports purification of the native enzyme from the (-) mating type, peptide sequencing and cloning of the gene. It names no assay and no product. Tier held at T3; raise to T1 on a full-text check.", "8828234")
ref_of = {(s, a): (g, pm) for s, rs in REF.items() for a, g, pm, v in rs}
srows = []
for step, rid, acc, tier, assay, pmid, q in L:
    gene, _ = ref_of[(step, acc)]
    org = species(acc) if acc else NULL_SPECIES[gene]
    assert pubmed_abstract(pmid)["title"]
    if acc and pmid != _ :
        pass
    qt = Q[pmid] if q else ""
    assert tier in ("T3",) or qt, (step, rid)
    srows.append([step, rid, org, gene, acc or "", tier, assay, pmid, qt])
srows.append(["torularhodin_formation", C + "torulene_to_torularhodin", "not assigned (no enzyme identified)", "unidentified", "", "T4",
              "hypothesis; a Rhea search and the pathway review found no enzyme", "", ""])
tsv(KB / "step_reaction.tsv", COLUMNS["step_reaction"], srows)
# change-detector fixture for the tests (quotes of all T1-T3 links)
tsv(ROOT / "tests" / "data" / "kb" / "carotenoid_quotes.tsv", ["step_id", "reaction_id", "protein_accession", "tier", "assay", "pmid", "quote"],
    [[r[0], r[1], r[4], r[5], r[6], r[7], r[8]] for r in srows if r[5] != "T4"])
for r in srows:
    assert len(r[8]) <= 400, (r[7], len(r[8]))
assert all(k in {r[7] for r in srows} for k in Q), [k for k in Q if k not in {r[7] for r in srows}]

# ---------- steps ----------
STEPS = [  # id, name, ec, domains, copy, role, fused, optional, route, lineage_note, comment
 ("ggpp_synthase", "geranylgeranyl diphosphate synthase (al-3, carG, crtE)", "2.5.1.29", ["PF00348"], "multi", "catalytic", [], False, None,
  "GGPP synthases also serve protein prenylation and other terpenoid pathways. Yeasts without carotenoids carry one (the pathway review names S. cerevisiae BTS1, Q12051). Presence of this step alone says nothing about carotenoids.", "polyprenyl synthase family"),
 ("phytoene_synthase", "phytoene synthase (C-terminal domain of al-2 / carRA / carRP / crtYB)", "2.5.1.32", ["PF00494"], "single", "catalytic", ["lycopene_cyclase"], False, None,
  "One protein with lycopene_cyclase in the reference fungi. References are cut to the phytoene synthase region. PF00494 also holds squalene synthase (ERG9), present in nearly all fungi (decoys).", "squalene/phytoene synthase family"),
 ("lycopene_cyclase", "lycopene beta-cyclase (N-terminal domain of al-2 / carRA / carRP / crtYB)", "5.5.1.19", ["TIGR03462", "PF18916"], "single", "catalytic", ["phytoene_synthase"], False, None,
  "References are cut to the cyclase region. InterPro places Pfam PF18916 on P37295 only among the five references; NCBIfam TIGR03462 is on all five. Neurospora AL-2 preferentially does single cyclisations (PMID 16928467).", "CarR-type lycopene cyclase domain"),
 ("phytoene_desaturase", "phytoene desaturase (al-1, carB, crtI)", None, ["PF01593"], "single", "catalytic", [], False, None,
  "Two reactions. The number of desaturations cannot be read from the sequence family: four steps to lycopene are shown for Phycomyces carB (PMID 9079885) and X. dendrorhous crtI (PMID 10397832); five steps to 3,4-didehydrolycopene are shown for Neurospora AL-1 (PMID 11017770, with lycopene as the second major product) and Fusarium carB (PMID 19645721). Only the five-step enzyme gives torulene and neurosporaxanthin. A tier holds for the named organism and protein only.", "FAD-dependent amine oxidase family (CrtI type)"),
 ("carx", "beta-carotene cleaving oxygenase, retinal-forming (CarX, Cco1)", "1.13.11.63", ["PF03055"], "single", "catalytic", [], True, "retinal",
  "Linked to opsins. Not needed for neurosporaxanthin (PMID 17293483). carX and carT are close relatives with different substrates.", "carotenoid cleavage oxygenase family"),
 ("cart", "torulene cleaving oxygenase (CarT; cao-2 in Neurospora)", "1.13.11.59", ["PF03055"], "single", "catalytic", [], True, "neurosporaxanthin",
  "Shown in Neurospora and Fusarium. The order of cyclisation, cleavage and oxidation differs between the two (abstract of PMID 21749649).", "carotenoid cleavage oxygenase family"),
 ("ylo1", "beta-apo-4'-carotenal dehydrogenase (ylo-1, carD)", "1.2.1.82", ["PF00171"], "single", "catalytic", [], True, "neurosporaxanthin",
  "Shown in Neurospora and Fusarium.", "aldehyde dehydrogenase family"),
 ("astaxanthin_synthase", "astaxanthin synthase, cytochrome P450 (crtS / asy)", None, ["PF00067"], "single", "catalytic", [], True, "astaxanthin",
  "Shown in Xanthophyllomyces dendrorhous (UniProt organism name: Phaffia rhodozyma). Needs the P450 reductase step crtr. The route is a reaction grid, not a line. Functional paper in the fungus: PMID 16416328; its link to Q3HR17 rests on sequence identity with the patent sequence CAC16261.1 (AX034665.1) and shared authors.", "cytochrome P450 family"),
 ("crtr", "cytochrome P450 reductase, electron donor of astaxanthin synthase (crtR)", "1.6.2.4", ["PF00667", "PF00258", "PF00175"], "single", "auxiliary", [], True, "astaxanthin",
  "A general cytochrome P450 reductase. Most fungi carry one, so it is a weak marker; it is needed for the reactions of astaxanthin_synthase.", "diflavin reductase family"),
 ("cars", "beta-carotene C11'-C12' cleaving oxygenase (carS of Phycomyces)", None, ["PF03055"], "single", "catalytic", [], True, "trisporoid",
  "Mucorales, sexual signalling. Gene: PMID 22445959; product beta-apo-12'-carotenal: PMID 21854466. The name carS belongs to Phycomyces blakesleeanus here. In Fusarium fujikuroi carS names a different gene (a regulator, per the pathway review); that gene is not a reference of this step.", "carotenoid cleavage oxygenase family"),
 ("acaa", "beta-apo-12'-carotenal C13-C14 cleaving oxygenase (acaA of Phycomyces)", None, ["PF03055"], "single", "catalytic", [], True, "trisporoid",
  "Mucorales, sexual signalling.", "carotenoid cleavage oxygenase family"),
 ("tsp3", "carotene oxygenase induced by trisporic acid (tsp3 of Blakeslea trispora)", None, ["PF03055"], "single", "catalytic", [], True, "trisporoid_tsp3",
  "Mucorales, sexual signalling. The product is not identified in the abstract; the link to the CarS reaction is T3.", "carotenoid cleavage oxygenase family"),
 ("tsp2", "4-dihydrotrisporin dehydrogenase (tsp2 of Mucor mucedo)", None, ["PF00106"], "single", "catalytic", [], True, "trisporoid_late",
  "Mucorales, sexual signalling. Enzyme activity is found only in the (-) mating type although the protein is in both (PMID 18931040).", "short-chain dehydrogenase family"),
 ("tdh", "4-dihydromethyltrisporate dehydrogenase (tdh of Mucor mucedo)", None, ["PF00248"], "single", "catalytic", [], True, "trisporoid_late",
  "Mucorales, sexual signalling. Protein purified from the (-) mating type.", "aldo/keto reductase family"),
]
CORE = ["ggpp_synthase", "phytoene_synthase", "phytoene_desaturase", "lycopene_cyclase"]
ROUTES = [
 ("core", "GGPP to beta-carotene and torulene", CORE, "Shared by all routes. Torulene needs the five-step desaturase."),
 ("neurosporaxanthin", "torulene to neurosporaxanthin (Neurospora, Fusarium)", CORE + ["cart", "ylo1"], "The order of cyclisation, cleavage and oxidation differs between Neurospora and Fusarium (PMID 21749649)."),
 ("retinal", "beta-carotene to retinal (Fusarium CarX, Ustilago Cco1)", CORE + ["carx"], None),
 ("astaxanthin", "beta-carotene to astaxanthin (Xanthophyllomyces dendrorhous)", CORE + ["crtr", "astaxanthin_synthase"], "A grid of ketolation and hydroxylation reactions by one P450 with its reductase. crtr comes first in the step order because its link is the first ketolation (beta-carotene to echinenone), which needs the reductase."),
 ("trisporoid", "beta-carotene to beta-apo-13-carotenone by CarS then AcaA (Phycomyces, mating)", CORE + ["cars", "acaa"], "Two cleavage steps."),
 ("trisporoid_tsp3", "beta-carotene cleavage by tsp3 (Blakeslea, Rhizopus; mating)", CORE + ["tsp3"], "tsp3 is an alternative of cars, not a successor of acaa. Its product is not identified in the abstract (T3)."),
 ("trisporoid_late", "trisporoid dehydrogenases tsp2 and tdh (Mucor mucedo, mating)", ["tsp2", "tdh"], "Not a chain. The substrates and products are not in ChEBI, so the tables cannot join these steps to the cleavage steps or to each other. Later oxidations are excluded (T4)."),
]
step_rx = {}
for r in srows:
    step_rx.setdefault(r[0], [])
    if r[1] not in step_rx[r[0]]:
        step_rx[r[0]].append(r[1])
def y(s):  # YAML double-quoted scalar
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'
out = [textwrap.dedent('''\
    # Carotenoid pathway definition, step set v2 (stage 1, task S3, 2026-11-21).
    #
    # Status: DRAFT for user review. Pathway-level verified: false.
    # Do not run a census against this file until the user has reviewed the step list.
    #
    # Reference checks: refs/carotenoid/provenance.tsv, refs/carotenoid/decoys/provenance.tsv and
    # refs/carotenoid/UNVERIFIED.md. Evidence tiers, quotes and reactions: kb/carotenoid/*.tsv.
    #
    # A ref has verified: true only when, on 2026-11-21, (1) the UniProt entry was fetched and its
    # organism and gene name matched, (2) the UniProt entry or the PubMed abstract supports the
    # reaction of the step, and (3) the PMID exists in PubMed. It never means a person reviewed it.
    # A decoy has verified: true when the UniProt entry was fetched on 2026-11-21, its organism and
    # gene name matched, it carries the Pfam family of the step, and its reviewed UniProt record
    # states a function that is not the function of the step. A decoy pmid is given only when the
    # abstract was read.
    #
    # Field notes:
    # - species: the organism name of the UniProt entry, copied exactly. UniProt names
    #   Xanthophyllomyces dendrorhous "Phaffia rhodozyma", Fusarium fujikuroi IMI 58289
    #   "Gibberella fujikuroi (strain CBS 195.34 / IMI 58289 / NRRL A-6831)" and Ustilago maydis
    #   "Mycosarcoma maydis".
    # - gene_name: a gene name of the UniProt entry. It is valid only together with species.
    # - region: [start, end], 1-based inclusive, in coordinates of the FASTA record of the step.
    #   For the fused phytoene synthase / lycopene cyclase proteins the FASTA record is already cut
    #   to the domain, so region covers the whole record. The coordinates on the full UniProt
    #   sequence are in the comment above each region and in the `region` column of
    #   refs/carotenoid/provenance.tsv. Reason: defs.validate checks region against the record in
    #   <step>.faa, so full-protein coordinates of a cut C-terminal domain do not validate.
    #   Full sequences: refs/carotenoid/full_length/fused.faa.
    # - domains_any: Pfam / NCBIfam ids from the UniProt cross-references of the references. The
    #   families are broader than the steps. A step call needs HMM plus tree placement with decoys.
    # - reactions: ids of kb/carotenoid/reaction.tsv. The tier of each (protein, reaction) link is
    #   in kb/carotenoid/step_reaction.tsv, not here.
    # - copy_expect, optional and role are PROPOSED values. No census data exists yet.
    # - Hypothesis-level steps (T4) with no reference protein are not steps here. They are listed in
    #   refs/carotenoid/UNVERIFIED.md; torulene to torularhodin also has a T4 row in the tables.
    id: carotenoid
    verified: false
    routes:
    ''')]
for rid, name, sids, note in ROUTES:
    out.append(f"  - id: {rid}\n    name: {y(name)}\n    step_ids: [{', '.join(sids)}]\n" + (f"    note: {y(note)}\n" if note else ""))
out.append("steps:\n")
krows = []
for sid, name, ec, doms, copy, role, fused, opt, route, note, family in STEPS:
    o = [f"  - id: {sid}", f"    name: {y(name)}", f"    ec: {y(ec) if ec else 'null'}"]
    if route:
        o.append(f"    route: {route}")
    o.append("    refs:")
    for acc, gene, pmid, ver in REF[sid]:
        sp = species(acc) if acc else NULL_SPECIES[gene]
        o += [f"      # {c}" for c in REF_NOTE.get((sid, acc), [])]
        o += [f"      - accession: {acc if acc else 'null'}", f"        species: {y(sp)}", f"        gene_name: {y(gene)}",
              f"        pmid: {y(pmid)}", f"        verified: {str(ver).lower()}"]
        reg = REGION.get(sid, {}).get(acc)
        if reg:
            # defs.validate checks region against the record in <step>.faa, which is the cut
            # domain. So region is in record coordinates; full-protein coordinates are a comment
            # here and the `region` column of provenance.tsv.
            o.append(f"        # FASTA record = residues {reg[0]}-{reg[1]} of the full UniProt sequence (UniProt feature Region)")
            o.append(f"        region: [1, {reg[1] - reg[0] + 1}]")
    if DECOY.get(sid):
        o.append("    decoys:")
        for acc, gene, pmid, reason in DECOY[sid]:
            o += [f"      - accession: {acc}", f"        species: {y(species(acc))}", f"        gene_name: {y(gene)}",
                  f"        pmid: {y(pmid) if pmid else 'null'}", "        verified: true"]
    o += [f"    domains_any: [{', '.join(doms)}]", f"    copy_expect: {copy}", f"    role: {role}",
          f"    fused_with: [{', '.join(fused)}]", f"    optional: {str(opt).lower()}",
          "    reactions:"] + [f"      - {y(r)}" for r in step_rx[sid]] + [f"    lineage_note: {y(note)}", ""]
    out.append("\n".join(o) + "\n")
    krows.append([sid, route or "", "", ";".join(doms), family, str(opt).lower()])
krows.append(["torularhodin_formation", "", "", "", "unknown (no enzyme identified)", "true"])
(ROOT / "pathways" / "carotenoid.yaml").write_text("".join(out).rstrip("\n") + "\n", encoding="utf-8")
tsv(KB / "step.tsv", COLUMNS["step"], krows)
print("refs", len(prov), "decoys", len(dprov), "compounds", len(crows), "reactions", len(RX), "links", len(srows))
