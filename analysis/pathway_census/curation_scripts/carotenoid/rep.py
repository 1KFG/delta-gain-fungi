from h import *
import csv
from pathlib import Path
from pathway_census.defs import load_pathway
from pathway_census.kb import load_kb
R = PKG
p = load_pathway(R / "pathways/carotenoid.yaml"); kb = load_kb(R / "kb/carotenoid")
norm = lambda s: " ".join(s.split())
prov = {(r["step_id"], r["accession"]): r for r in csv.DictReader(open(R / "refs/carotenoid/provenance.tsv"), delimiter="\t")}
o = ["### References (every row re-fetched from UniProt and PubMed on 2026-11-21)", "",
     "| step | accession | UniProt organism | gene name in entry | reviewed | length | region on full protein | ref PMID | PMID exists | entry cites PMID | verified |",
     "|---|---|---|---|---|---|---|---|---|---|---|"]
for s in p.steps:
    for r in s.refs:
        pm = "yes" if pubmed_abstract(r.pmid)["title"] else "no"
        if not r.accession:
            o.append(f"| {s.id} | none | {r.species} (name from the paper) | {r.gene_name} (not found in UniProt) | - | - | - | {r.pmid} | {pm} | - | false |"); continue
        e = uniprot_entry(r.accession)
        o.append(f"| {s.id} | {r.accession} | {e['organism']} | {r.gene_name} ({'match' if r.gene_name in e['gene_names'] else 'NO MATCH'}) | {'yes' if e['reviewed'] else 'no'} | {e['length']} | {prov[(s.id, r.accession)]['region'] or '-'} | {r.pmid} | {pm} | {'yes' if r.pmid in e['pmids'] else 'no'} | {str(r.verified).lower()} |")
o += ["", "### Protein-reaction links (kb/carotenoid/step_reaction.tsv)", "",
      "| step | reaction | organism | gene | accession | tier | PMID | quote present | quote found in fetched abstract |", "|---|---|---|---|---|---|---|---|---|"]
for _, r in kb["step_reaction"].iterrows():
    found = "-"
    if r["quote"]:
        found = "yes" if norm(r["quote"]) in norm(pubmed_abstract(r["pmid"])["abstract"]) else "NO"
    o.append(f"| {r['step_id']} | {r['reaction_id']} | {r['organism'].split(' (')[0]} | {r['gene_name']} | {r['protein_accession'] or 'none'} | {r['tier']} | {r['pmid'] or '-'} | {'yes' if r['quote'] else 'no'} | {found} |")
o += ["", "### Quoted sentences (one per paper)", ""]
seen = {}
for _, r in kb["step_reaction"].iterrows():
    if r["quote"] and r["pmid"] not in seen:
        seen[r["pmid"]] = r["quote"]; o.append(f"- PMID {r['pmid']}: \"{r['quote']}\"")
o += ["", "### Decoys (every accession re-fetched from UniProt on 2026-11-21)", "",
      "| accession | UniProt organism | gene name (match) | reviewed | EC in entry | steps | PMID read | reason |", "|---|---|---|---|---|---|---|---|"]
drows = list(csv.DictReader(open(R / "refs/carotenoid/decoys/provenance.tsv"), delimiter="\t"))
by = {}
for r in drows:
    by.setdefault(r["accession"], {"r": r, "steps": []})["steps"].append(r["step_id"])
for acc, d in by.items():
    e = uniprot_entry(acc); r = d["r"]
    o.append(f"| {acc} | {e['organism']} | {r['gene_name']} ({'match' if r['gene_name'] in e['gene_names'] else 'NO MATCH'}) | {'yes' if e['reviewed'] else 'no'} | {', '.join(e['ec']) or '-'} | {', '.join(d['steps'])} | {r['pmid'] or 'none'} | {r['reason'] if len(d['steps'])==1 or acc!='Q9Y8G7' else r['reason']} |")
o += ["", "### Compounds (ChEBI and PubChem fetched on 2026-11-21)", "", "| ChEBI | ChEBI name | formula | computed mass | InChIKey | PubChem CID | first block agrees | KEGG xref (from ChEBI) |", "|---|---|---|---|---|---|---|---|"]
for _, c in kb["compound"].iterrows():
    ch = chebi_compound(c["chebi_id"]); pc = pubchem_by_inchikey(c["inchikey"])
    o.append(f"| {c['chebi_id']} | {ch['name']} | {c['formula']} | {c['neutral_exact_mass']} | {c['inchikey']} | {pc['cid']} | {'yes' if pc['inchikey'][:14]==c['inchikey_block1'] else 'NO'} | {c['kegg_id'] or '-'} |")
Path("tables.md").write_text("\n".join(o) + "\n")
print(len(o))
