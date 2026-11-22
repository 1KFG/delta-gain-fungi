#!/usr/bin/env python3
"""Build the region table for the EBC comparison from the per-genome antiSMASH subset tables.

One row per antiSMASH region: location, class, contig-edge flag, the rank-1 MIBiG
known-cluster hit (accession, compound names, similarity), a status flag, MIBiG compound
identifiers (PubChem CID, formula, mass) and, for hits at or above --inchikey-min, the InChIKey
looked up from PubChem by CID. Biosynthetic gene locus tags are listed per region.

Status (similarity of the rank-1 MIBiG hit, the antiSMASH 0-100 score):
  known >= 75, weak 1-74, none (no hit). The cutoff is a PROPOSAL; the similarity is kept in
  the table so the reader can choose another.

Usage: build_ebc_table.py RESULTS_DIR MIBIG_JSON_DIR OUT.tsv [--inchikey-min 75] [--no-pubchem]
Reads RESULTS_DIR/<id>.{regions,cds,known}.tsv.gz.
"""
import argparse, csv, glob, gzip, json, os, re, sys, urllib.request


def rows(path):
    with gzip.open(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def mibig_compounds(jdir, acc):
    base = acc.split(".")[0]
    p = os.path.join(jdir, base + ".json")
    if not os.path.exists(p):
        return []
    d = json.load(open(p))
    out = []
    for c in d.get("compounds", []):
        cid = ""
        for x in c.get("databaseIds", []):
            if x.startswith("pubchem:"):
                cid = x.split(":", 1)[1]
        out.append({"name": c.get("name", ""), "cid": cid, "formula": c.get("formula", ""),
                    "mass": c.get("mass", "")})
    return out


def inchikeys(cids):
    res = {}
    cids = sorted(set(c for c in cids if re.fullmatch(r"\d+", c)))
    for i in range(0, len(cids), 100):
        chunk = cids[i:i + 100]
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/InChIKey/JSON" % ",".join(chunk)
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                for p in json.load(r)["PropertyTable"]["Properties"]:
                    res[str(p["CID"])] = p["InChIKey"]
        except Exception as e:  # leave blank rather than guess
            print("PubChem lookup failed for a chunk:", e, file=sys.stderr)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results"); ap.add_argument("mibig"); ap.add_argument("out")
    ap.add_argument("--inchikey-min", type=int, default=75)
    ap.add_argument("--no-pubchem", action="store_true")
    a = ap.parse_args()
    table = []
    for rf in sorted(glob.glob(os.path.join(a.results, "*.regions.tsv.gz"))):
        gid = os.path.basename(rf).split(".")[0]
        base = os.path.join(a.results, gid)
        reg = rows(rf)
        cds = rows(base + ".cds.tsv.gz")
        known = rows(base + ".known.tsv.gz")
        genes = {}
        for c in cds:
            if c["gene_kind"] == "biosynthetic":
                genes.setdefault((c["record"], c["region"]), []).append(c["locus_tag"])
        top = {(k["record"], k["region"]): k for k in known if k["rank"] == "1"}
        for r in reg:
            key = (r["record"], r["region"])
            k = top.get(key)
            sim = int(float(k["similarity"])) if k and k["similarity"] not in ("", "None") else None
            status = "none" if sim is None else ("known" if sim >= 75 else "weak")
            comps = mibig_compounds(a.mibig, k["mibig_accession"]) if k else []
            table.append({
                "genome": gid, "record": r["record"], "region": r["region"], "start": r["start"], "end": r["end"],
                "products": r["products"], "contig_edge": r["contig_edge"], "n_cds": r["n_cds"],
                "status": status, "mibig_accession": k["mibig_accession"] if k else "",
                "mibig_compounds": ";".join(c["name"] for c in comps) or (k["hit_description"] if k else ""),
                "mibig_class": k["hit_cluster_type"] if k else "",
                "similarity": "" if sim is None else sim, "core_gene_hits": k["core_gene_hits"] if k else "",
                "pubchem_cids": ";".join(c["cid"] for c in comps if c["cid"]),
                "formulas": ";".join(c["formula"] for c in comps if c["formula"]),
                "inchikeys": "",
                "biosynthetic_genes": ",".join(genes.get(key, [])),
            })
    if not a.no_pubchem:
        want = [c for t in table if t["similarity"] != "" and t["similarity"] >= a.inchikey_min
                for c in t["pubchem_cids"].split(";") if c]
        ik = inchikeys(want)
        for t in table:
            if t["similarity"] != "" and t["similarity"] >= a.inchikey_min:
                t["inchikeys"] = ";".join(ik.get(c, "") for c in t["pubchem_cids"].split(";") if c)
    cols = list(table[0].keys())
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(table)
    n = len(table)
    print(n, "regions;", {s: sum(1 for t in table if t["status"] == s) for s in ("known", "weak", "none")})


if __name__ == "__main__":
    main()
