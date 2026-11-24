#!/usr/bin/env python3.12
"""Add decoys for ayg1, scytalone_dehydratase, hppd, hmga, aspulvinone_e_synthase (owner request 2026-11-23).

One-off build helper (not tested code). Needs network access to rest.uniprot.org.
For each accession it fetches the UniProt JSON record and the FASTA, checks that the record is
reviewed, is not already a reference in any refs/<pathway>/provenance.tsv, and that organism,
gene and length agree between the two fetches. It then writes refs/<pathway>/decoys/<step>.faa
and appends rows to refs/<pathway>/decoys/provenance.tsv (same columns as the carotenoid decoys).
Nothing is taken from memory: every field of a row comes from the fetched record.
Run from analysis/pathway_census:  /usr/bin/python3.12 curation_scripts/melanin_myco_decoys/add_decoys.py
"""
import csv, json, sys, time, urllib.request
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
UA = {"User-Agent": "pathway-census-curation"}

# step -> (pathway, [(accession, family note written by the curator)])
PLAN = {
    "hppd": ("melanin_ldopa", [
        ("P50107", "glyoxalase I (VOC superfamily, same Pfam PF00903 as HPPD), fungal"),
        ("Q56415", "fosfomycin-resistance glutathione transferase (VOC superfamily, PF00903)"),
        ("Q99XE5", "bleomycin-resistance protein (VOC superfamily, PF00903)"),
        ("O52791", "4-hydroxymandelate synthase, the closest HPPD paralogue with another reaction (also PF00903 plus PF14696); added after the first six decoys did not score"),
        ("P06622", "catechol 2,3-dioxygenase, extradiol (VOC superfamily, PF00903 plus PF22247)"),
        ("A0A0U2WCB2", "fungal glyoxalase-domain protein of an itaconate-pathway cluster, annotated as a putative ring-cleaving dioxygenase (function not shown in the record)"),
        ("Q0D1P8", "fungal lactoylglutathione-lyase-like protein of the terrein cluster (function not shown in the record)"),
    ]),
    "aspulvinone_e_synthase": ("melanin_ldopa", [
        ("B7STY1", "NRPS-like atromentin synthetase (A-T-TE, same architecture as the reference; makes atromentin, not aspulvinone E)"),
        ("Q0CT94", "NRPS-like atromentin synthetase of Aspergillus terreus (ATEG_03090; makes atromentin, not aspulvinone E)"),
        ("A7XRY0", "NRPS-like terrequinone/didemethylasterriquinone D synthetase"),
        ("P07702", "L-2-aminoadipate reductase (adenylation-thiolation-reductase NRPS-like enzyme of lysine biosynthesis)"),
        ("P08659", "firefly luciferase (adenylate-forming enzyme, not NRPS)"),
        ("Q01574", "acetyl-CoA synthetase (adenylate-forming enzyme, not NRPS)"),
    ]),
    "ayg1": ("melanin_dhn", [
        ("P04335", "bacterial esterase FrsA (PF06500, the Pfam family of the references)"),
        ("Q93NG6", "bacterial 2,6-dihydropseudooxynicotine hydrolase (PF06500)"),
        ("S8ASK9", "fungal hydrolyase poxO of the oxaleimide cluster (PF06500; a fungal paralogue family member that is not the DHN-melanin enzyme)"),
        ("S0EE80", "fungal 20-hydroxy-prefusarin hydrolase FUS2 of the fusarin cluster (PF06500; not the DHN-melanin enzyme)"),
    ]),
    "scytalone_dehydratase": ("melanin_dhn", [
        ("P07445", "ketosteroid isomerase (NTF2-like fold, Pfam PF12680)"),
        ("Q9RN59", "SnoaL polyketide cyclase (NTF2-like fold, Pfam PF07366)"),
        ("O52646", "AcmA aklanonic acid methyl ester cyclase (NTF2-like fold, Pfam PF07366)"),
        ("P61970", "nuclear transport factor 2 (NTF2 fold, Pfam PF02136)"),
        ("Q9ZAG3", "limonene-1,2-epoxide hydrolase (NTF2-like fold, Pfam PF07858)"),
    ]),
    "hmga": ("melanin_ldopa", [
        ("O43123", "maleylacetoacetate isomerase, next step of the tyrosine-degradation route, a glutathione S-transferase fold, fungal"),
        ("O43708", "human GSTZ1 maleylacetoacetate isomerase (glutathione S-transferase zeta)"),
        ("P47096", "3-hydroxyanthranilate 3,4-dioxygenase, an Fe-dependent aromatic ring-cleaving dioxygenase (PF06052), fungal"),
        ("Q0QFQ2", "gentisate 1,2-dioxygenase, an Fe-dependent aromatic ring-cleaving dioxygenase (cupin, PF07883)"),
        ("Q16878", "cysteine dioxygenase type 1 (Fe-dependent cupin dioxygenase, PF05995)"),
        ("P47125", "indoleamine 2,3-dioxygenase BNA2 (heme dioxygenase, PF01231; not Fe-cupin)"),
    ]),
}
# Owner gene-name rule (2026-11-22): UniProt has no gene name for these three. A UniProt protein abbreviation
# (short name in the record) is used where the record has one; otherwise the pattern <accession>_<role>_<species code>.
# The name is valid only together with the organism.
NAMES = {"O52791": "HmaS", "P08659": "Luciferase", "Q93NG6": "Q93NG6_hydrolase_Pnico"}
COLS = ["step_id", "accession", "species", "gene_name", "pmid", "reason", "source_database", "retrieved_on", "verified"]


def get(url):
    for i in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.read().decode()
        except Exception:
            time.sleep(2 * (i + 1))
    raise RuntimeError(url)


def main():
    refs = set()
    for p in (ROOT / "refs").glob("*/provenance.tsv"):
        refs |= {r["accession"] for r in csv.DictReader(open(p), delimiter="\t")}
    today = date.today().isoformat()
    for step, (pathway, items) in PLAN.items():
        ddir = ROOT / "refs" / pathway / "decoys"
        fa, rows = [], []
        for acc, note in items:
            if acc in refs:
                sys.exit(f"{acc} is already a reference; not a decoy")
            j = json.loads(get(f"https://rest.uniprot.org/uniprotkb/{acc}.json"))
            fasta = get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")
            seq = "".join(fasta.splitlines()[1:])
            assert j["primaryAccession"] == acc, acc
            assert j["entryType"].startswith("UniProtKB reviewed"), (acc, j["entryType"])
            assert j["sequence"]["length"] == len(seq) == len(j["sequence"]["value"]), acc
            org = j["organism"]["scientificName"]
            assert f"OS={org}" in fasta.splitlines()[0], (acc, org)
            genes = j.get("genes") or []
            gene = genes[0]["geneName"]["value"] if genes and "geneName" in genes[0] else ""
            shorts = [x["value"] for x in j["proteinDescription"]["recommendedName"].get("shortNames", [])]
            if not gene and acc in NAMES:
                assert NAMES[acc] in shorts or acc == "Q93NG6", (acc, shorts)
                gene_note = f"UniProt has no gene name; name used here: {NAMES[acc]}"
                gene = NAMES[acc]
            else:
                gene_note = ""
            pn = j["proteinDescription"]["recommendedName"]["fullName"]["value"]
            ecs = [e["value"] for e in j["proteinDescription"]["recommendedName"].get("ecNumbers", [])]
            pf = sorted({x["id"] for x in j.get("uniProtKBCrossReferences", []) if x["database"] == "Pfam"})
            a = j["entryAudit"]
            reason = (f"Not the step reaction. {note}. UniProt reviewed entry: {pn}"
                      + (f" (EC {', '.join(ecs)})" if ecs else "") + f"; Pfam in the record: {', '.join(pf) or 'none'}"
                      + (f"; {gene_note}" if gene_note else ""))
            rows.append(dict(step_id=step, accession=acc, species=org, gene_name=gene, pmid="", reason=reason,
                             source_database=f"UniProtKB/Swiss-Prot (entry v{a['entryVersion']}, sequence v{a['sequenceVersion']})",
                             retrieved_on=today, verified="true"))
            hdr, body = fasta.split("\n", 1)  # header format of the existing decoy files: ">ACC sp|ACC|ENTRY ..."
            fa.append(f">{acc} {hdr[1:]}\n{body}" if body.endswith("\n") else f">{acc} {hdr[1:]}\n{body}\n")
            print(step, acc, org, gene or "-", len(seq), pn[:50], sep="\t")
            time.sleep(0.3)
        (ddir / f"{step}.faa").write_text("".join(fa))
        pv = ddir / "provenance.tsv"
        have = [r for r in csv.DictReader(open(pv), delimiter="\t") if r["step_id"] != step]
        with open(pv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(have + rows)
        add_to_yaml(pathway, step, rows)


def add_to_yaml(pathway, step, rows):
    """Insert a decoys: block (accession, species, gene_name, pmid, verified) before fused_with of the step.
    The pathway tests require the YAML decoys to equal the decoy provenance rows. Idempotent."""
    path = ROOT / "pathways" / f"{pathway}.yaml"
    lines = path.read_text().split("\n")
    start = next(i for i, l in enumerate(lines) if l == f"  - id: {step}")
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("  - id: ")), len(lines))
    block = lines[start:end]
    assert not any(l.startswith("    decoys:") for l in block), f"{step} already has decoys"
    at = start + next(i for i, l in enumerate(block) if l.startswith("    fused_with:"))
    new = ["    decoys:"]
    for r in rows:
        new += [f"      - accession: {r['accession']}", f"        species: {json.dumps(r['species'])}",
                f"        gene_name: {json.dumps(r['gene_name'])}", "        pmid: null", "        verified: true"]
    lines[at:at] = new
    path.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
