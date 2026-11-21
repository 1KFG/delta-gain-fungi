"""Second pass from the written files (not from build.py's data): for every reference and decoy
re-read the UniProt record (cache) and compare organism, gene name and sequence; check that each
PMID exists and whether the UniProt entry cites it; check that each quote is in the abstract of
its PMID; compare each compound with ChEBI. Prints Markdown tables (used in the task report).
Exit code 1 on any mismatch."""
from h import *
import csv
from pathway_census.defs import load_pathway, validate
from pathway_census.kb import load_kb, validate_kb, formula_mass, path_possible, species_key

PW = "melanin_ldopa"
NO_GENE = {"Q2TL94"}  # UniProt entry without a gene name (UNVERIFIED.md, A1)
UNREVIEWED_DECOY = {"Q2UNF9"}  # AoCO4: unreviewed entry, function in two abstracts (UNVERIFIED.md, F)
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

p = load_pathway(PKG / "pathways" / f"{PW}.yaml")
assert validate(p) == [], validate(p)
kb = load_kb(PKG / "kb" / PW)
assert validate_kb(kb) == [], validate_kb(kb)
sr = kb["step_reaction"]

print("### References\n\n| step | accession | UniProt organism = species | gene name in entry | sequence = FASTA | ref PMID exists | entry cites ref PMID | reviewed |\n|---|---|---|---|---|---|---|---|")
for st in p.steps:
    seqs = fasta(PKG / "refs" / PW / f"{st.id}.faa")
    for r in st.refs:
        e = uniprot_entry(r.accession)
        org = e["organism"] == r.species
        gene = (r.gene_name in e["gene_names"]) if r.accession not in NO_GENE else None
        seq = seqs.get(r.accession) == e["sequence"]
        pm = bool(pubmed_abstract(r.pmid)["title"])
        if not (org and seq and pm and gene is not False):
            bad.append(("ref", st.id, r.accession))
        g = "no gene name in the entry" if gene is None else ("yes" if gene else "NO")
        print(f"| {st.id} | {r.accession} | {'yes' if org else 'NO'} ({e['organism']}) | {g} ({r.gene_name}) | {'yes' if seq else 'NO'} ({e['length']} aa) | {'yes' if pm else 'NO'} ({r.pmid}) | {'yes' if r.pmid in e['pmids'] else 'no'} | {'yes' if e['reviewed'] else 'no'} |")

print("\n### Decoys\n\n| step | accession | UniProt organism = species | gene name in entry | sequence = FASTA | reviewed | EC on the entry | PMID (read) |\n|---|---|---|---|---|---|---|---|")
for st in p.steps:
    if not st.decoys:
        continue
    seqs = fasta(PKG / "refs" / PW / "decoys" / f"{st.id}.faa")
    for d in st.decoys:
        e = uniprot_entry(d.accession)
        ok = e["organism"] == d.species and d.gene_name in e["gene_names"] and seqs.get(d.accession) == e["sequence"] and (e["reviewed"] or d.accession in UNREVIEWED_DECOY)
        if d.pmid:
            ok = ok and d.pmid in e["pmids"] and bool(pubmed_abstract(d.pmid)["title"])
        if not ok:
            bad.append(("decoy", st.id, d.accession))
        print(f"| {st.id} | {d.accession} | {'yes' if e['organism'] == d.species else 'NO'} ({e['organism']}) | {'yes' if d.gene_name in e['gene_names'] else 'NO'} ({d.gene_name}) | {'yes' if seqs.get(d.accession) == e['sequence'] else 'NO'} | {'yes' if e['reviewed'] else 'no'} | {', '.join(e['ec']) or '-'} | {d.pmid or '-'} |")

print("\n### Links (step_reaction.tsv)\n\n| step | reaction | accession | organism | gene | tier | PMID | PMID exists | quote found in the abstract | entry cites the PMID |\n|---|---|---|---|---|---|---|---|---|---|")
for _, r in sr.iterrows():
    ex = found = "-"
    if r["pmid"]:
        ab = pubmed_abstract(r["pmid"])
        ex = "yes" if ab["title"] else "NO"
        if r["quote"]:
            found = "yes" if norm(r["quote"]) in norm(ab["abstract"]) else "NO"
    if "NO" in (ex, found):
        bad.append(("link", r["step_id"], r["protein_accession"], r["pmid"]))
    cites = "-"
    if r["pmid"] and r["protein_accession"]:
        cites = "yes" if r["pmid"] in uniprot_entry(r["protein_accession"])["pmids"] else "no"
    print(f"| {r['step_id']} | {r['reaction_id']} | {r['protein_accession'] or '-'} | {r['organism'].split(' (')[0]} | {r['gene_name']} | {r['tier']} | {r['pmid'] or '-'} | {ex} | {found} | {cites} |")

print("\n### Compounds\n\n| ChEBI | name | formula | neutral exact mass | InChIKey = ChEBI | PubChem CID (same first block) | KEGG |\n|---|---|---|---|---|---|---|")
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
    print(f"| {r['chebi_id']} | {c['name']} | {r['formula']} | {r['neutral_exact_mass']} | {'yes' if ok else 'NO'} | {cid} | {r['kegg_id'] or '-'} |")

print("\n### Routes (path_possible, all organisms)\n\n| route | steps | product | tier | organisms | single_organism | n_equal_paths |\n|---|---|---|---|---|---|---|")
for rt in p.routes:
    got = path_possible(kb, rt.step_ids)
    if got.empty:
        print(f"| {rt.id} | {', '.join(rt.step_ids)} | none | - | - | - | - |")
    for _, g in got.iterrows():
        print(f"| {rt.id} | {', '.join(rt.step_ids)} | {g['compound']} | {g['min_tier']} | {g['organisms']} | {g['single_organism']} | {g['n_equal_paths']} |")
print("\n### Routes, one organism at a time (longest prefix of the route that the organism completes)\n\n| route | organism | steps completed | products of the last completed step (tier) |\n|---|---|---|---|")
orgs = sorted({species_key(o) for o in sr[sr["tier"] != "T4"]["organism"]})
for rt in p.routes:
    for org in orgs:
        best = None
        for n in range(1, len(rt.step_ids) + 1):
            got = path_possible(kb, rt.step_ids[:n], organism=org)
            if got.empty:
                break
            best = (n, "; ".join(f"{c} ({t})" for c, t in zip(got["compound"], got["min_tier"])))
        if best:
            print(f"| {rt.id} | {org} | {best[0]} of {len(rt.step_ids)} | {best[1]} |")
print("\nmismatches:", bad)
sys.exit(1 if bad else 0)
