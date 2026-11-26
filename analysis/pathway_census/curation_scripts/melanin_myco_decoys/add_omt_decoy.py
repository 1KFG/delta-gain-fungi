#!/usr/bin/env python3.12
"""Add UniProt Q8NKC1 (S. pombe 'Probable catechol O-methyltransferase 2') as a decoy of the mycosporine step omt (2026-11-23).

One-off build helper (not tested code). Needs network access to rest.uniprot.org.
It fetches the UniProt JSON record and the FASTA, checks that the record is reviewed, is not a reference of
any pathway, and that accession, organism and length agree between the two fetches. It then appends the
sequence to refs/mycosporine/decoys/omt.faa, a row to refs/mycosporine/decoys/provenance.tsv and an entry
to the omt `decoys:` block of pathways/mycosporine.yaml. Every field comes from the fetched record.
Run from analysis/pathway_census:  /usr/bin/python3.12 curation_scripts/melanin_myco_decoys/add_omt_decoy.py
"""
import csv, json, sys, urllib.request
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ACC = "Q8NKC1"
UA = {"User-Agent": "pathway-census-curation"}
COLS = ["step_id", "accession", "species", "gene_name", "pmid", "reason", "source_database", "retrieved_on", "verified"]


def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
        return r.read().decode()


def main():
    refs = set()
    for p in (ROOT / "refs").glob("*/provenance.tsv"):
        refs |= {r["accession"] for r in csv.DictReader(open(p), delimiter="\t")}
    if ACC in refs:
        sys.exit(f"{ACC} is a reference")
    pv = ROOT / "refs/mycosporine/decoys/provenance.tsv"
    if any(r["accession"] == ACC and r["step_id"] == "omt" for r in csv.DictReader(open(pv), delimiter="\t")):
        sys.exit(f"{ACC} is already an omt decoy")
    j = json.loads(get(f"https://rest.uniprot.org/uniprotkb/{ACC}.json"))
    fasta = get(f"https://rest.uniprot.org/uniprotkb/{ACC}.fasta")
    seq = "".join(fasta.splitlines()[1:])
    assert j["primaryAccession"] == ACC
    assert j["entryType"].startswith("UniProtKB reviewed")
    assert j["sequence"]["length"] == len(seq) == len(j["sequence"]["value"])
    org = j["organism"]["scientificName"]
    assert f"OS={org}" in fasta.splitlines()[0]
    g = j["genes"][0]
    gene = g["geneName"]["value"] if "geneName" in g else g["orfNames"][0]["value"]  # ORF name, as for O42898 (SPBC119.03)
    rn = j["proteinDescription"]["recommendedName"]
    name, ecs = rn["fullName"]["value"], [e["value"] for e in rn.get("ecNumbers", [])]
    pf = sorted({x["id"] for x in j["uniProtKBCrossReferences"] if x["database"] == "Pfam"})
    rhea = sorted({x["id"] for c in j["comments"] if c["commentType"] == "CATALYTIC ACTIVITY"
                   for x in c["reaction"].get("reactionCrossReferences", []) if x["id"].startswith("RHEA:")})
    a = j["entryAudit"]
    print(ACC, org, gene, name, len(seq), ecs, rhea, pf, a["entryVersion"], a["sequenceVersion"], sep="\t")
    reason = (f"S-adenosylmethionine-dependent O-methyltransferase (Pfam {', '.join(pf)}); UniProt reviewed entry: {name[0].lower() + name[1:]} of "
              f"{org} (EC {', '.join(ecs)}, {', '.join(rhea)}). annotation only: the name says 'Probable'; the entry is not assayed with "
              "demethyl-4-deoxygadusol; the step has no EC and no Rhea reaction, so the check is the other annotation of the entry. "
              "Scores 120.7 on the omt step HMM, the highest of the omt decoys (MEASURED, step 2 side check, omt_check). "
              "gene_name is the ordered locus name (UniProt gives no gene name). Sister of the decoy O42898 (COMT1). "
              "It is not excluded as a true member of the step; the step needs tree placement; do not count a hit to this protein as a negative (USER DECISION 7)")
    row = dict(step_id="omt", accession=ACC, species=org, gene_name=gene, pmid="", reason=reason,
               source_database=f"UniProtKB/Swiss-Prot (entry v{a['entryVersion']}, sequence v{a['sequenceVersion']})",
               retrieved_on=date.today().isoformat(), verified="true")
    hdr, body = fasta.split("\n", 1)
    with open(ROOT / "refs/mycosporine/decoys/omt.faa", "a") as fh:
        fh.write(f">{ACC} {hdr[1:]}\n{body}" + ("" if body.endswith("\n") else "\n"))
    with open(pv, "a", newline="") as fh:
        csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n").writerow(row)
    path = ROOT / "pathways/mycosporine.yaml"
    lines = path.read_text().split("\n")
    start = next(i for i, l in enumerate(lines) if l == "  - id: omt")
    end = next(i for i in range(start + 1, len(lines)) if lines[i].startswith("  - id: "))
    at = start + next(i for i, l in enumerate(lines[start:end]) if l.startswith("    domains_any:"))
    lines[at:at] = [f"      - accession: {ACC}", f"        species: {json.dumps(org)}", f"        gene_name: {json.dumps(gene)}",
                    "        pmid: null", "        verified: true"]
    path.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
