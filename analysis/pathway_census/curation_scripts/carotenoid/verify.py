"""Second pass from the written files only."""
from h import *
from pathlib import Path
from pathway_census.defs import load_pathway, validate
from pathway_census.kb import load_kb, validate_kb, formula_mass
R = PKG
p = load_pathway(R / "pathways/carotenoid.yaml"); assert validate(p) == []
kb = load_kb(R / "kb/carotenoid"); assert validate_kb(kb) == []
def fasta(path):
    out, n = {}, None
    for ln in Path(path).read_text().splitlines():
        if ln.startswith(">"): n = ln[1:].split()[0]; out[n] = ""
        else: out[n] += ln.strip()
    return out
norm = lambda s: " ".join(s.split())
nref = ndec = 0
import csv
prov = {(r["step_id"], r["accession"]): r for r in csv.DictReader(open(R / "refs/carotenoid/provenance.tsv"), delimiter="\t")}
for s in p.steps:
    fa = fasta(R / f"refs/carotenoid/{s.id}.faa")
    for r in s.refs:
        assert pubmed_abstract(r.pmid)["title"]
        if not r.accession:
            assert not r.verified; print("  no accession:", s.id, r.species, r.gene_name, r.pmid); continue
        e = uniprot_entry(r.accession)
        assert e["organism"] == r.species and r.gene_name in e["gene_names"], (s.id, r.accession)
        seq = e["sequence"]
        reg = prov[(s.id, r.accession)]["region"]
        if reg:
            a, b = map(int, reg.split("-")); seq = seq[a-1:b]
        assert fa[r.accession] == seq, (s.id, r.accession)
        cited = r.pmid in e["pmids"]
        if not cited: print("  UniProt does not cite the ref PMID:", s.id, r.accession, r.pmid, "entry cites", e["pmids"])
        nref += 1
    if s.decoys:
        da = fasta(R / f"refs/carotenoid/decoys/{s.id}.faa")
        for d in s.decoys:
            e = uniprot_entry(d.accession)
            assert e["organism"] == d.species and d.gene_name in e["gene_names"] and e["reviewed"], (s.id, d.accession)
            assert da[d.accession] == e["sequence"]
            ndec += 1
sr = kb["step_reaction"]
nq = 0
for _, r in sr.iterrows():
    if r["pmid"]:
        ab = pubmed_abstract(r["pmid"])
        if r["quote"]:
            assert norm(r["quote"]) in norm(ab["abstract"]), (r["step_id"], r["pmid"]); nq += 1
    if r["protein_accession"]:
        e = uniprot_entry(r["protein_accession"])
        assert e["organism"] == r["organism"] and r["gene_name"] in e["gene_names"]
        if r["pmid"] not in e["pmids"]:
            print("  link PMID not cited by the UniProt entry:", r["step_id"], r["protein_accession"], r["pmid"], r["tier"])
for _, c in kb["compound"].iterrows():
    ch = chebi_compound(c["chebi_id"]); pc = pubchem_by_inchikey(c["inchikey"])
    assert ch["inchikey"] == c["inchikey"] and ch["formula"] == c["formula"] and ch["smiles"] == c["smiles"] and int(ch["charge"]) == 0
    assert pc["inchikey"][:14] == c["inchikey_block1"] and pc["formula"] == c["formula"]
    assert abs(float(c["neutral_exact_mass"]) - pc["exact_mass"]) < 0.0005
print("refs with accession checked", nref, "decoy rows checked", ndec, "quotes found in abstract", nq, "compounds", len(kb["compound"]))
