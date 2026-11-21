"""Protein records for the mycosporine curation: one interface for UniProtKB entries and for the
three GenBank proteins of the Phaffia rhodozyma cluster (no UniProtKB entry).

rec(acc) -> dict(accession, source ('UniProtKB' or 'GenBank'), organism, species (the text used
in the tables), gene_names, sequence, length, ec, rhea, xrefs (Pfam, NCBIfam, PROSITE ids),
pmids, reviewed, header (FASTA header), srcdb (provenance text)).

Gene names: kbsources.uniprot_entry reads geneName, synonyms and orfNames. It does not read
orderedLocusNames. Cyanobacterial entries have ordered locus names only (Npun_R5600, Ava_3858),
so they are added here from the full UniProt record.
"""
from helpers import *
from genbank_proteins import gp

def is_genbank(acc):
    return acc.startswith("WEV")

def full(acc):
    return raw("uniprot_full", acc, f"https://rest.uniprot.org/uniprotkb/{acc}.json")

def ordered_locus_names(acc):
    return [o["value"] for g in full(acc).get("genes", []) for o in g.get("orderedLocusNames", [])]

def rec(acc):
    if is_genbank(acc):
        d = gp(acc)
        nuc = "OQ547787.1"
        assert "OQ547787.1" in d["coded_by"], acc
        return {"accession": acc, "source": "GenBank", "organism": d["organism"],
                "species": f"{d['organism']} (strain {d['strain']})", "gene_names": [d["gene"]],
                "sequence": d["sequence"], "length": d["length"], "ec": [], "rhea": set(), "xrefs": set(),
                "pmids": [], "reviewed": False, "taxid": d["taxid"], "refs": d["refs"],
                "header": f"{acc} gb|{acc}| {d['product']} OS={d['organism']} OX={d['taxid']} GN={d['gene']} strain={d['strain']} coded_by={nuc}",
                "srcdb": f"GenBank protein {acc} (CDS of nucleotide record OQ547787.1; no UniProtKB entry)"}
    e, d = uniprot_entry(acc), full(acc)
    assert e["accession"] == acc, acc
    rhea = set()
    for c in d.get("comments", []):
        if c["commentType"] == "CATALYTIC ACTIVITY":
            rhea |= {x["id"] for x in c["reaction"].get("reactionCrossReferences", []) if x["database"] == "Rhea" and x["id"].startswith("RHEA:")}
    xrefs = {x["id"] for x in d.get("uniProtKBCrossReferences", []) if x["database"] in ("Pfam", "NCBIfam", "PROSITE")}
    db = "sp" if e["reviewed"] else "tr"
    pdn = d["proteinDescription"]
    nm = (pdn.get("recommendedName") or (pdn.get("submissionNames") or [{}])[0])["fullName"]["value"]
    genes = e["gene_names"] + [g for g in ordered_locus_names(acc) if g not in e["gene_names"]]
    au = d["entryAudit"]
    kind = "Swiss-Prot" if e["reviewed"] else "TrEMBL"
    return {"accession": acc, "source": "UniProtKB", "organism": e["organism"], "species": e["organism"],
            "gene_names": genes, "sequence": e["sequence"], "length": e["length"], "ec": e["ec"], "rhea": rhea,
            "xrefs": xrefs, "pmids": e["pmids"], "reviewed": e["reviewed"], "taxid": str(d["organism"]["taxonId"]),
            "header": f"{acc} {db}|{acc}|{d['uniProtkbId']} {nm} OS={e['organism']} OX={d['organism']['taxonId']}" + (f" GN={genes[0]}" if genes else ""),
            "srcdb": f"UniProtKB/{kind} (entry v{au['entryVersion']}, sequence v{au['sequenceVersion']})"}

def function_text(acc):
    return " ".join(t["value"] for c in full(acc).get("comments", []) if c["commentType"] == "FUNCTION" for t in c.get("texts", []))
