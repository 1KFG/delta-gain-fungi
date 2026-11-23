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


# ---------------------------------------------------------------------------------------------
# Second set (2026-11-22): phytoene_desaturase and cart. These proteins have one domain, so the
# sequences are FULL LENGTH and nothing is cut. Tree members only, as above.
# ---------------------------------------------------------------------------------------------
HMM_DIR = PKG / "results" / "hmms" / "carotenoid"
MIN_SCORE, MIN_COVERAGE = 10.0, 0.5          # the rule of bin/filter_decoys_for_tree.py (PROPOSED there)
FAMILY2 = {"phytoene_desaturase": "Belongs to the carotenoid/retinoid oxidoreductase family",
           "cart": "Belongs to the carotenoid oxygenase family"}
# step: {accession: (organism start, gene, protein name, EC, expected length)}
MEMBERS2 = {"phytoene_desaturase": {"P48537": ("Cercospora nicotianae", "PDH1", "Phytoene desaturase", "1.3.99.30", 621)},
            "cart": {}}
# Fetched and not written: the same gene as a reference, from another strain record.
# accession: (step, reference accession, organism start, number of different residues at equal length)
SAME_AS_REFERENCE = {"S0EPU6": ("phytoene_desaturase", "Q8X0Z0", "Gibberella fujikuroi", 0),
                     "S0E9A6": ("cart", "A1KQY3", "Gibberella fujikuroi", 7)}
# step: {accession: (organism start, gene, protein name, EC, type)}
OUTGROUPS2 = {
    "phytoene_desaturase": {
        "P21685": ("Pantoea ananas", "crtI", "Phytoene desaturase (lycopene-forming)", "1.3.99.31", "bacterial CrtI (Gammaproteobacteria)"),
        "P17054": ("Rhodobacter capsulatus", "crtI", "Phytoene desaturase (neurosporene-forming)", "1.3.99.28", "bacterial CrtI (Alphaproteobacteria)"),
        "P54979": ("Myxococcus xanthus", "carA2", "zeta-carotene-forming phytoene desaturase", "1.3.99.29", "bacterial phytoene desaturase (Myxococcota)"),
        "Q02861": ("Myxococcus xanthus", "carC", "All-trans-zeta-carotene desaturase", "1.3.99.26", "bacterial zeta-carotene desaturase of the CrtI family (Myxococcota)"),
        "Q7WT72": ("Flavobacterium sp.", "crtD", "1-hydroxycarotenoid 3,4-desaturase", "1.3.99.27", "bacterial CrtD (other function: 3,4-desaturase)"),
        "Q64FW2": ("Mus musculus", "Retsat", "All-trans-retinol 13,14-reductase", "1.3.99.23", "animal RetSat (other function, CrtISO subfamily)"),
        "Q8S4R4": ("Solanum lycopersicum", "CRTISO", "Prolycopene isomerase, chloroplastic", "5.2.1.13", "plant CRTISO (other function, CrtISO subfamily)"),
    },
    "cart": {
        "Q9HAY6": ("Homo sapiens", "BCO1", "Beta,beta-carotene 15,15'-dioxygenase", "1.13.11.63", "animal BCO1 (retinal-forming)"),
        "Q9JJS6": ("Mus musculus", "Bco1", "Beta,beta-carotene 15,15'-dioxygenase", "1.13.11.63", "animal BCO1 (retinal-forming)"),
        "Q9I993": ("Gallus gallus", "BCO1", "Beta,beta-carotene 15,15'-dioxygenase", "1.13.11.63", "animal BCO1 (retinal-forming)"),
        "Q16518": ("Homo sapiens", "RPE65", "Retinoid isomerohydrolase", "3.1.1.64", "animal RPE65 (other function: retinoid isomerohydrolase)"),
        "Q28175": ("Bos taurus", "RPE65", "Retinoid isomerohydrolase", "3.1.1.64", "animal RPE65 (other function: retinoid isomerohydrolase)"),
    },
}
# MEASURED 2026-11-22: best domain score (bits) and fraction of the HMM that this one domain covers.
# hmmsearch -E 1000 --domE 1000 (HMMER 3.4) with results/hmms/carotenoid/<step>.PROVISIONAL.hmm
# (phytoene_desaturase: 579 match states, 4 references; cart: 621 match states, 2 references).
# main() measures again when hmmsearch and the HMM are present and stops on a difference.
ALIGNABILITY = {
    "phytoene_desaturase": {"P48537": (776.2, 0.801), "P21685": (288.7, 0.853), "P17054": (317.0, 0.874),
                            "P54979": (494.4, 0.862), "Q02861": (302.7, 0.865), "Q7WT72": (323.1, 0.851),
                            "Q64FW2": (45.5, 0.544), "Q8S4R4": (41.4, 0.501)},
    "cart": {"Q9HAY6": (132.0, 0.628), "Q9JJS6": (142.9, 0.633), "Q9I993": (133.7, 0.628),
             "Q16518": (118.2, 0.625), "Q28175": (116.2, 0.625)},
}


def lineage2(d):
    """(class, phylum) of a non-fungal record: Bacteria, Metazoa (Chordata) or seed plants."""
    lin = d["organism"]["lineage"]
    assert "Fungi" not in lin, d["primaryAccession"]
    if "Metazoa" in lin:
        cls = [x for x in lin if x in ("Mammalia", "Aves")]
        assert "Chordata" in lin and len(cls) == 1, lin
        return cls[0], "Chordata"
    return lineage(d, False)


def seqs_of(path):
    out, k = {}, None
    for l in path.read_text().splitlines():
        if l.startswith(">"):
            k = l[1:].split()[0]
            out[k] = ""
        elif l.strip():
            out[k] += l.strip()
    return out


def measure_alignability(step, fasta_text):
    """Best domain score and HMM coverage with the code of bin/filter_decoys_for_tree.py. None when not possible here."""
    import importlib.util, shutil, tempfile
    hmm = HMM_DIR / f"{step}.PROVISIONAL.hmm"
    if not shutil.which("hmmsearch") or not hmm.exists():
        return None
    spec = importlib.util.spec_from_file_location("filter_decoys_for_tree", PKG / "bin" / "filter_decoys_for_tree.py")
    fd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fd)
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "q.faa"
        fa.write_text(fasta_text, encoding="utf-8")
        best, n = fd.best_scores(hmm, fa), fd.hmm_length(hmm)
    return {k: (round(sc, 1), round((b - a + 1) / n, 3)) for k, (sc, a, b) in best.items()}


def second_set(full_fa):
    """Rows and files for phytoene_desaturase and cart. Appends the full-length named members to full_fa."""
    rows, files = [], {}
    taken = {s: ids_of(REFS / f"{s}.faa") | ids_of(REFS / "decoys" / f"{s}.faa") for s in MEMBERS2}
    for acc, (step, ref, org, ndiff) in SAME_AS_REFERENCE.items():    # keep the reason for the exclusion checkable
        d = full(acc)
        assert d["organism"]["scientificName"].startswith(org) and ref in taken[step], acc
        a, b = d["sequence"]["value"], seqs_of(REFS / f"{step}.faa")[ref]
        assert len(a) == len(b) and sum(x != y for x, y in zip(a, b)) == ndiff, acc
    for step, members in MEMBERS2.items():
        recs = []
        for acc, (org, g, name, ec, length) in members.items():
            d = full(acc)
            assert d["primaryAccession"] == acc and d["organism"]["scientificName"].startswith(org), acc
            assert gene(d) == g and d["proteinDescription"]["recommendedName"]["fullName"]["value"] == name, acc
            assert all_ecs(d) == {ec} and acc not in taken[step], acc
            sim = [c["texts"][0]["value"] for c in d["comments"] if c["commentType"] == "SIMILARITY"]
            assert sim and sim[0].startswith(FAMILY2[step]), (acc, sim)
            seq = d["sequence"]["value"]
            assert len(seq) == length, acc
            # This record has no evidence code on the EC number or on the catalytic activity line.
            # The only evidence code is on the FUNCTION text. The script stops if that changes.
            fev = {e["evidenceCode"] for c in d["comments"] if c["commentType"] == "FUNCTION"
                   for t in c["texts"] for e in t.get("evidences", [])}
            ecev = [e for e in d["proteinDescription"]["recommendedName"]["ecNumbers"] for e in e.get("evidences", [])]
            ecev += [e for c in d["comments"] if c["commentType"] == "CATALYTIC ACTIVITY" for e in c["reaction"].get("evidences", [])]
            assert fev == {"ECO:0000250"} and not ecev and not function_pmids(d), (acc, fev, ecev)
            cls, phy = lineage(d, True)
            sc, cov = ALIGNABILITY[step][acc]
            full_fa.append(f">{header(d)}\n{wrap(seq)}\n")
            recs.append(f">{header(d).replace(f'{acc} ', f'{acc} region=1-{length} ', 1)}\n{wrap(seq)}\n")
            reason = ("Named fungal member for the tree. The UniProt record gives the function by similarity only "
                      "(ECO:0000250 on the FUNCTION text, no source protein named). The EC number and the catalytic "
                      "activity line carry no evidence code. No experiment on this protein is cited. "
                      f"Full-length sequence, nothing cut. Step HMM best domain score {sc} bits, HMM coverage {cov}. "
                      "Not a reference, not a decoy, not evidence for a step call.")
            rows.append([step, acc, d["organism"]["scientificName"], cls, phy, g, "", "named_member",
                         d["proteinExistence"], d["annotationScore"], "ECO:0000250 (by similarity)", ec, f"1-{length}",
                         reason, f"{srcdb(d)}; full length; UniProt family comment '{sim[0]}'", RETRIEVED, "true"])
        if recs:                                   # no file for a step with no named member
            files[f"{step}.faa"] = "".join(recs)
            got = measure_alignability(step, "".join(recs))
            assert got is None or all(got[a] == ALIGNABILITY[step][a] for a in members), (step, got)
    for step, ogs in OUTGROUPS2.items():
        recs = []
        for acc, (org, g, name, ec, typ) in ogs.items():
            d = full(acc)
            assert d["primaryAccession"] == acc and d["organism"]["scientificName"].startswith(org), acc
            assert gene(d) == g and d["proteinDescription"]["recommendedName"]["fullName"]["value"] == name, acc
            assert all_ecs(d) == {ec} and acc not in taken[step] and acc not in MEMBERS2[step], (acc, all_ecs(d))
            sim = [c["texts"][0]["value"] for c in d["comments"] if c["commentType"] == "SIMILARITY"]
            assert sim and sim[0].startswith(FAMILY2[step]), (acc, sim)
            ev = ec_evidence(d, ec)
            assert "ECO:0000269" in ev, (acc, ev)   # every outgroup has experimental evidence on its EC number
            sc, cov = ALIGNABILITY[step][acc]
            assert sc >= MIN_SCORE and cov >= MIN_COVERAGE, acc
            cls, phy = lineage2(d)
            recs.append(f">{header(d)}\n{wrap(d['sequence']['value'])}\n")
            reason = (f"Tree outgroup: {typ}, EC {ec}, same UniProt family as the step, not fungal. "
                      f"Alignability: step HMM best domain score {sc} bits, HMM coverage {cov} "
                      f"(rule: score >= {MIN_SCORE:g} bits and coverage >= {MIN_COVERAGE}). Full-length sequence. "
                      "Not a reference, not a decoy, not evidence for a step call.")
            rows.append([step, acc, d["organism"]["scientificName"], cls, phy, g, ",".join(function_pmids(d)),
                         "tree_outgroup", d["proteinExistence"], d["annotationScore"], ev, ec, "", reason,
                         f"{srcdb(d)}; UniProt family comment '{sim[0]}'", RETRIEVED, "true"])
        files[f"outgroups_{step}.faa"] = "".join(recs)
        got = measure_alignability(step, "".join(recs))
        if got is None:
            print(f"{step}: hmmsearch or the HMM is not present; alignability was not measured again")
        else:
            assert all(got[a] == ALIGNABILITY[step][a] for a in ogs), (step, got)
    return rows, files


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
    rows2, files2 = second_set(full_fa)          # appended after the rows and records of the first set
    rows += rows2
    for name, text in files2.items():
        (HERE / name).write_text(text, encoding="utf-8")
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
