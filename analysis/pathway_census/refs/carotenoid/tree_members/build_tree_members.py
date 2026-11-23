"""Write the tree-member files of this folder from fetched UniProt and InterPro records.

Run: /usr/bin/python3.12 analysis/pathway_census/refs/carotenoid/tree_members/build_tree_members.py
The fetches go through pathway_census.kbsources (rate limit, user agent, cache in .cache/).
The script stops on any record that does not match the expected organism, name or EC.
It writes only into its own folder. It is not run by the tests.
"""
import csv, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PKG = HERE.parents[2]                       # analysis/pathway_census
sys.path.insert(0, str(PKG / "curation_scripts" / "carotenoid"))
from h import raw                            # noqa: E402

RETRIEVED = "2026-11-22"
REFS = PKG / "refs" / "carotenoid"

# accession: (expected start of organism name, expected EC list)
NAME = "Bifunctional lycopene cyclase/phytoene synthase"
MEMBERS = {
    "Q67GH9": "Blakeslea trispora",
    "B2ATB0": "Podospora anserina",
    "B2WAQ3": "Pyrenophora tritici-repentis",
    "E4ZUB5": "Leptosphaeria maculans",
    "Q0V6M5": "Parastagonospora nodorum",
    "D4D802": "Trichophyton verrucosum",
    "Q4P0S6": "Mycosarcoma maydis",
    "D1Z4K7": "Sordaria macrospora",
    "E3Q717": "Colletotrichum graminicola",
    "E4UPP6": "Nannizzia gypsea",
    "Q2U4X9": "Aspergillus oryzae",
    "A2QM49": "Aspergillus niger",
    "E6ZZ11": "Sporisorium reilianum",
}
# Fetched and not written: no UniProt Region feature, so the cut rule gives no cut.
EXCLUDED = {"S0EQ07": "Gibberella fujikuroi"}

# accession: (organism start, gene, protein name, EC, type, reason)
OUTGROUPS = {
    "P21683": ("Pantoea ananas", "crtB", "15-cis-phytoene synthase", "2.5.1.32", "bacterial CrtB (Gammaproteobacteria)"),
    "P54975": ("Paracoccus sp.", "crtB", "15-cis-phytoene synthase", "2.5.1.32", "bacterial CrtB (Alphaproteobacteria)"),
    "P37269": ("Synechococcus elongatus", "crtB", "15-cis-phytoene synthase", "2.5.1.32", "cyanobacterial CrtB"),
    "P37294": ("Synechocystis sp.", "crtB", "15-cis-phytoene synthase", "2.5.1.32", "cyanobacterial CrtB"),
    "P37271": ("Arabidopsis thaliana", "PSY1", "Phytoene synthase 1, chloroplastic", "2.5.1.32", "plant PSY"),
    "P37272": ("Capsicum annuum", "PSY1", "Bifunctional 15-cis-phytoene synthase, chromoplastic", "2.5.1.32", "plant PSY"),
    "Q2FV59": ("Staphylococcus aureus", "crtM", "4,4'-diapophytoene synthase", "2.5.1.96", "bacterial CrtM (dehydrosqualene synthase, C30)"),
}
STEP_EC = {"phytoene_synthase": "2.5.1.32", "lycopene_cyclase": "5.5.1.19"}
REGION_NAME = {"phytoene_synthase": "Phytoene synthase", "lycopene_cyclase": "Lycopene beta-cyclase"}
FAMILY = {"phytoene_synthase": "PF00494", "lycopene_cyclase": "TIGR03462"}
ECO = {"ECO:0000269": "experimental", "ECO:0000250": "by similarity", "ECO:0000305": "curator inference",
       "ECO:0000303": "non-traceable author statement", "ECO:0000255": "sequence model", "ECO:0000256": "sequence model"}


def full(acc):
    return raw("uniprot_full", acc, f"https://rest.uniprot.org/uniprotkb/{acc}.json")


def interpro(acc):
    ip = raw("interpro", acc, f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{acc}?page_size=100")
    loc = {}
    for r in ip["results"]:
        for p in r["proteins"]:
            for l in p["entry_protein_locations"]:
                loc.setdefault(r["metadata"]["accession"], []).extend((x["start"], x["end"]) for x in l["fragments"])
    return {k: sorted(v) for k, v in loc.items()}


def gene(d):
    g = d.get("genes", [{}])[0]
    if "geneName" in g:
        return g["geneName"]["value"]
    for key in ("orderedLocusNames", "orfNames"):
        if g.get(key):
            return g[key][0]["value"]
    return ""


def header(d):
    acc = d["primaryAccession"]
    assert d["entryType"].startswith("UniProtKB reviewed"), acc
    nm = d["proteinDescription"]["recommendedName"]["fullName"]["value"]
    h = f"{acc} sp|{acc}|{d['uniProtkbId']} {nm} OS={d['organism']['scientificName']} OX={d['organism']['taxonId']}"
    g = gene(d)
    return h + (f" GN={g}" if g else "")


def wrap(seq):
    return "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))


def ec_evidence(d, ec):
    """Evidence codes on the EC number and on the catalytic-activity lines with this EC."""
    ev = []
    pd = d["proteinDescription"]
    for p in [pd] + pd.get("includes", []):
        for e in (p.get("recommendedName") or {}).get("ecNumbers", []):
            if e["value"] == ec:
                ev += e.get("evidences", [])
    for c in d.get("comments", []):
        if c["commentType"] == "CATALYTIC ACTIVITY" and c["reaction"].get("ecNumber") == ec:
            ev += c["reaction"].get("evidences", [])
    assert ev, (d["primaryAccession"], ec, "no evidence on the EC number")
    out = sorted({f"{e['evidenceCode']} ({ECO[e['evidenceCode']]}{', ' + e['source'] + ':' + e['id'] if e.get('id') else ''})" for e in ev})
    return "; ".join(out)


def function_pmids(d):
    out = []
    for r in d["references"]:
        pos = r.get("referencePositions", [])
        if "FUNCTION" in pos or "CATALYTIC ACTIVITY" in pos:
            out += [x["id"] for x in r["citation"].get("citationCrossReferences", []) if x["database"] == "PubMed"]
    return out


def lineage(d, fungal):
    lin = d["organism"]["lineage"]
    if fungal:
        assert "Fungi" in lin, d["primaryAccession"]
        phy = [x for x in lin if x.endswith("mycota")]
        cls = [x for x in lin if x.endswith("mycetes")]
        assert len(phy) == 1 and len(cls) == 1, (d["primaryAccession"], lin)
        return cls[0], phy[0]
    assert "Fungi" not in lin, d["primaryAccession"]
    if lin[0] == "Bacteria":                 # Bacteria > kingdom > phylum > class
        return lin[3], lin[2]
    assert "Streptophyta" in lin and "Magnoliopsida" in lin, lin
    return "Magnoliopsida", "Streptophyta"


def srcdb(d):
    au = d["entryAudit"]
    return f"UniProtKB/Swiss-Prot (entry v{au['entryVersion']}, sequence v{au['sequenceVersion']})"


def all_ecs(d):
    pd = d["proteinDescription"]
    return {e["value"] for p in [pd] + pd.get("includes", []) for e in (p.get("recommendedName") or {}).get("ecNumbers", [])}


def ids_of(path):
    return {l[1:].split()[0] for l in path.read_text().splitlines() if l.startswith(">")}


def main():
    taken = {s: ids_of(REFS / f"{s}.faa") | ids_of(REFS / "decoys" / f"{s}.faa") for s in STEP_EC}
    taken["lycopene_cyclase"] |= ids_of(REFS / "outgroups" / "lycopene_cyclase.faa")
    rows, fa = [], {"phytoene_synthase": [], "lycopene_cyclase": []}
    full_fa = []
    for acc, org in MEMBERS.items():
        d = full(acc)
        assert d["primaryAccession"] == acc and d["organism"]["scientificName"].startswith(org), acc
        assert d["proteinDescription"]["recommendedName"]["fullName"]["value"] == NAME, acc
        assert set(STEP_EC.values()) <= all_ecs(d), acc
        reg = {f["description"]: (f["location"]["start"]["value"], f["location"]["end"]["value"])
               for f in d["features"] if f["type"] == "Region"}
        psy, cyc = reg["Phytoene synthase"], reg["Lycopene beta-cyclase"]
        assert cyc[1] < psy[0], acc
        loc = interpro(acc)
        cls, phy = lineage(d, True)
        seq = d["sequence"]["value"]
        full_fa.append(f">{header(d)}\n{wrap(seq)}\n")
        fp = function_pmids(d)
        for step in STEP_EC:
            assert acc not in taken[step], (step, acc)
            s, e = reg[REGION_NAME[step]]
            fam = loc[FAMILY[step]]
            assert all(s <= a and b <= e for a, b in fam), (acc, step, (s, e), fam)   # same check as build.py
            h = header(d).replace(f"{acc} ", f"{acc} region={s}-{e} ", 1)
            fa[step].append(f">{h}\n{wrap(seq[s - 1:e])}\n")
            ev = ec_evidence(d, STEP_EC[step])
            exp = "ECO:0000269" in ev
            reason = ("Named fungal member for the tree. " +
                      ("The UniProt record gives experimental evidence for the EC number." if exp else
                       "The UniProt record gives the EC number by similarity to P37295 only. No experiment on this protein is cited.") +
                      " Not a reference, not a decoy, not evidence for a step call.")
            src = (f"{srcdb(d)}; region = UniProt feature Region '{REGION_NAME[step]}'; it contains "
                   f"{'Pfam' if step == 'phytoene_synthase' else 'NCBIfam'} {FAMILY[step]} at "
                   + ",".join(f"{a}-{b}" for a, b in fam) + " (InterPro)")
            rows.append([step, acc, d["organism"]["scientificName"], cls, phy, gene(d), ",".join(fp), "named_member",
                         d["proteinExistence"], d["annotationScore"], ev, STEP_EC[step], f"{s}-{e}", reason, src,
                         RETRIEVED, "true"])
    for acc, org in EXCLUDED.items():           # keep the reason for the exclusion checkable
        d = full(acc)
        assert d["organism"]["scientificName"].startswith(org)
        assert not [f for f in d["features"] if f["type"] == "Region"], acc
    og = []
    for acc, (org, g, name, ec, typ) in OUTGROUPS.items():
        d = full(acc)
        assert d["primaryAccession"] == acc and d["organism"]["scientificName"].startswith(org), acc
        assert gene(d) == g and d["proteinDescription"]["recommendedName"]["fullName"]["value"] == name, acc
        assert all_ecs(d) == {ec}, (acc, all_ecs(d))
        assert acc not in taken["phytoene_synthase"], acc
        assert "PF00494" in {x["id"] for x in d["uniProtKBCrossReferences"] if x["database"] == "Pfam"}, acc
        cls, phy = lineage(d, False)
        og.append(f">{header(d)}\n{wrap(d['sequence']['value'])}\n")
        reason = (f"Tree outgroup for rooting: {typ}, EC {ec}, Pfam PF00494 in the UniProt record, not fungal. "
                  "Full-length sequence (the record has no Region feature for the domain). "
                  "Not a reference, not a decoy, not evidence for a step call.")
        rows.append(["phytoene_synthase", acc, d["organism"]["scientificName"], cls, phy, g, ",".join(function_pmids(d)),
                     "tree_outgroup", d["proteinExistence"], d["annotationScore"], ec_evidence(d, ec), ec, "", reason,
                     srcdb(d), RETRIEVED, "true"])
    for step, recs in fa.items():
        (HERE / f"{step}.faa").write_text("".join(recs), encoding="utf-8")
    (HERE / "full_length").mkdir(exist_ok=True)
    (HERE / "full_length" / "named_members.faa").write_text("".join(full_fa), encoding="utf-8")
    (HERE / "outgroups_phytoene_synthase.faa").write_text("".join(og), encoding="utf-8")
    cols = ["step_id", "accession", "species", "class", "phylum", "gene_name", "pmid", "role", "protein_existence",
            "annotation_score", "evidence_codes", "ec", "region", "reason", "source_database", "retrieved_on", "verified"]
    with open(HERE / "provenance.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar=None)
        w.writerow(cols)
        for r in rows:
            assert len(r) == len(cols) and not any("\t" in str(x) or "\n" in str(x) for x in r), r
            w.writerow(r)
    print(f"wrote {len(rows)} provenance rows")


if __name__ == "__main__":
    main()
