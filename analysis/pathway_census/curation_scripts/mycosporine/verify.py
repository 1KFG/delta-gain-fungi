"""Second pass from the written files (not from build.py's data): for every reference and decoy
re-read the record (UniProt, or GenBank for the three Phaffia proteins) and compare organism,
gene name and sequence; check that each PMID exists and whether the UniProt entry cites it; check
that each quote is in the abstract of its PMID; compare each compound with ChEBI; check each
lineage row has a PubMed record. Prints Markdown tables (used in the task report).
Exit code 1 on any mismatch. Run in this directory: /usr/bin/python3.12 verify.py"""
from records import *
import csv
from pathway_census.defs import load_pathway, validate
from pathway_census.kb import load_kb, validate_kb, formula_mass, path_possible, species_key

PW = "mycosporine"
bad = []
def norm(s):
    return " ".join(s.split())
def fasta(path):
    out, name = {}, None
    for ln in path.read_text(encoding="utf-8").splitlines():
        if ln.startswith(">"):
            name = ln[1:].split()[0]; out[name] = ""
        elif name:
            out[name] += ln.strip()
    return out
def yn(x):
    return "yes" if x else "NO"

p = load_pathway(PKG / "pathways" / f"{PW}.yaml")
assert validate(p) == [], validate(p)
kb = load_kb(PKG / "kb" / PW)
assert validate_kb(kb) == [], validate_kb(kb)
sr = kb["step_reaction"]

print("### References\n\n| step | accession | source | organism = species | gene name in record | sequence = FASTA | ref PMID exists | record cites ref PMID | reviewed | verified |\n|---|---|---|---|---|---|---|---|---|---|")
for st in p.steps:
    if not st.refs:
        print(f"| {st.id} | none (T4 step) | - | - | - | - | - | - | - | - |")
        continue
    seqs = fasta(PKG / "refs" / PW / f"{st.id}.faa")
    for r in st.refs:
        e = rec(r.accession)
        org, gene, seq = e["species"] == r.species, r.gene_name in e["gene_names"], seqs.get(r.accession) == e["sequence"]
        pm = bool(pubmed_abstract(r.pmid)["title"])
        if e["source"] == "GenBank":   # the GenBank record cites the paper by title (no PUBMED line)
            cites = "by title and DOI" if any(t.rstrip(".") == pubmed_abstract(r.pmid)["title"].rstrip(".") for t, *_ in e["refs"]) else "NO"
        else:
            cites = "yes" if r.pmid in e["pmids"] else "no"
        if not (org and gene and seq and pm) or cites == "NO":
            bad.append(("ref", st.id, r.accession))
        print(f"| {st.id} | {r.accession} | {e['source']} | {yn(org)} ({e['species']}) | {yn(gene)} ({r.gene_name}) | {yn(seq)} ({e['length']} aa) | {yn(pm)} ({r.pmid}) | {cites} | {'yes' if e['reviewed'] else 'no'} | {str(r.verified).lower()} |")

print("\n### Decoys\n\n| step | accession | organism = species | gene name in record | sequence = FASTA | reviewed | EC on the entry | Rhea on the entry | PMID (read) | verified |\n|---|---|---|---|---|---|---|---|---|---|")
for st in p.steps:
    if not st.decoys:
        continue
    seqs = fasta(PKG / "refs" / PW / "decoys" / f"{st.id}.faa")
    for d in st.decoys:
        e = rec(d.accession)
        ok = e["species"] == d.species and d.gene_name in e["gene_names"] and seqs.get(d.accession) == e["sequence"]
        if d.pmid:
            ok = ok and bool(pubmed_abstract(d.pmid)["title"])
        if st.ec and st.ec in e["ec"]:
            ok = False
        if set(st.reactions) & e["rhea"]:
            ok = False
        if not ok:
            bad.append(("decoy", st.id, d.accession))
        print(f"| {st.id} | {d.accession} | {yn(e['species'] == d.species)} ({e['species']}) | {yn(d.gene_name in e['gene_names'])} ({d.gene_name}) | {yn(seqs.get(d.accession) == e['sequence'])} | {'yes' if e['reviewed'] else 'no'} | {', '.join(e['ec']) or '-'} | {', '.join(sorted(e['rhea'])) or '-'} | {d.pmid or '-'} | {str(d.verified).lower()} |")

print("\n### Links (step_reaction.tsv)\n\n| step | reaction | accession | organism | gene | tier | PMID | PMID exists | quote found in the abstract | record cites the PMID |\n|---|---|---|---|---|---|---|---|---|---|")
for _, r in sr.iterrows():
    ex = found = "-"
    if r["pmid"]:
        a = pubmed_abstract(r["pmid"])
        ex = yn(a["title"])
        if r["quote"]:
            found = yn(norm(r["quote"]) in norm(a["abstract"]))
    if "NO" in (ex, found):
        bad.append(("link", r["step_id"], r["protein_accession"], r["pmid"]))
    cites = "-"
    if r["pmid"] and r["protein_accession"]:
        e = rec(r["protein_accession"])
        cites = "by title and DOI (GenBank)" if e["source"] == "GenBank" else ("yes" if r["pmid"] in e["pmids"] else "no")
    print(f"| {r['step_id']} | {r['reaction_id']} | {r['protein_accession'] or '-'} | {r['organism'].split(' (')[0]} | {r['gene_name']} | {r['tier']} | {r['pmid'] or '-'} | {ex} | {found} | {cites} |")

print("\n### Compounds\n\n| ChEBI | name | formula | neutral exact mass | InChIKey = ChEBI | PubChem CID (same InChIKey) | KEGG |\n|---|---|---|---|---|---|---|")
for _, r in kb["compound"].iterrows():
    c = chebi_compound(r["chebi_id"])
    ok = c["inchikey"] == r["inchikey"] and c["formula"] == r["formula"] and abs(formula_mass(r["formula"]) - c["monoisotopic_mass"]) < 0.0005
    try:
        pc = pubchem_by_inchikey(r["inchikey"])
        cid = str(pc["cid"]) if pc["inchikey"][:14] == r["inchikey_block1"] else "NO"
    except KbFetchError:
        cid = "none: PubChem has no compound with this InChIKey"
    if not ok or cid == "NO":
        bad.append(("compound", r["chebi_id"]))
    print(f"| {r['chebi_id']} | {c['name']} | {r['formula']} | {r['neutral_exact_mass']} | {yn(ok)} | {cid} | {r['kegg_id'] or '-'} |")

print("\n### Routes (path_possible, all organisms)\n\n| route or step list | steps | product | tier | organisms | single_organism | n_equal_paths |\n|---|---|---|---|---|---|---|")
calls = [(rt.id, rt.step_ids) for rt in p.routes] + [("prefix: ddgs", ["ddgs"]), ("prefix: ddgs, omt", ["ddgs", "omt"]), ("one step: atp_grasp_ligase", ["atp_grasp_ligase"])]
for name, steps in calls:
    got = path_possible(kb, steps)
    if got.empty:
        print(f"| {name} | {', '.join(steps)} | none | - | - | - | - |")
    for _, g in got.iterrows():
        print(f"| {name} | {', '.join(steps)} | {g['compound']} | {g['min_tier']} | {g['organisms']} | {g['single_organism']} | {g['n_equal_paths']} |")
print("\n### One organism at a time\n\n| step list | organism | product (tier) |\n|---|---|---|")
orgs = sorted({species_key(o) for o in sr[sr["tier"] != "T4"]["organism"]})
for name, steps in calls[3:]:
    for org in orgs:
        got = path_possible(kb, steps, organism=org)
        print(f"| {', '.join(steps)} | {org} | {'; '.join(f'{c} ({t})' for c, t in zip(got['compound'], got['min_tier'])) or 'none'} |")

with open(PKG / "refs" / PW / "fungal_lineages.tsv", encoding="utf-8", newline="") as fh:
    lin = list(csv.DictReader(fh, delimiter="\t"))
for r in lin:
    if not pubmed_abstract(r["pmid"])["title"]:
        bad.append(("lineage", r["organism"], r["pmid"]))
print(f"\nlineage rows: {len(lin)}, distinct PMIDs: {len({r['pmid'] for r in lin})}")
print("\nmismatches:", bad)
sys.exit(1 if bad else 0)
