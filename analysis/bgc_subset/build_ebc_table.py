#!/usr/bin/env python3
"""Build the region table for the EBC comparison from the per-genome antiSMASH subset tables.

One row per antiSMASH region: location, class, contig-edge flag, the rank-1 MIBiG
known-cluster hit (accession, compound names, similarity), a status flag, MIBiG compound
identifiers (PubChem CID, formula, mass) and, for hits at or above --inchikey-min, the InChIKey
looked up from PubChem by CID. Biosynthetic gene locus tags are listed per region.

Status (similarity of the best of the top-3 MIBiG hits, the antiSMASH 0-100 score):
  known >= 75, weak 1-74, none (no hit). The cutoff is a PROPOSAL; the similarity is kept in
  the table so the reader can choose another.

Adds entry size (n_entry_proteins), fraction of the entry's genes hit, synteny score, and MIBiG quality/status.
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


def mibig_meta(jdir, acc):
    p = os.path.join(jdir, acc.split(".")[0] + ".json")
    if not os.path.exists(p):
        return "", "", ""
    d = json.load(open(p))
    return d.get("quality", ""), d.get("status", ""), d.get("completeness", "")


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
        # antiSMASH ranks hits by core-gene hits and score, not by the similarity number, so the
        # rank-1 hit can have a lower similarity than rank 2 or 3. Use the highest similarity of
        # the (up to) 3 kept hits as the region's best hit; keep all three in top3_hits.
        hits = {}
        for k in known:
            hits.setdefault((k["record"], k["region"]), []).append(k)
        def simval(k):
            try:
                return int(float(k["similarity"]))
            except (ValueError, TypeError):
                return -1
        top = {key: max(v, key=lambda k: (simval(k), -int(k["rank"]))) for key, v in hits.items()}
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
                "best_hit_rank": k["rank"] if k else "",
                "top3_hits": ";".join("%s|%s|%s" % (h["mibig_accession"], h["hit_description"][:40], h["similarity"]) for h in sorted(hits.get(key, []), key=lambda h: int(h["rank"]))),
                "similarity": "" if sim is None else sim, "core_gene_hits": k["core_gene_hits"] if k else "",
                "gene_hits": k["gene_hits"] if k else "",
                "n_entry_proteins": k.get("n_entry_proteins", "") if k else "",
                "entry_fraction_hit": (round(int(k["gene_hits"]) / int(k["n_entry_proteins"]), 2) if k and k.get("n_entry_proteins") not in (None, "", "0") and k["gene_hits"] not in ("", "None") else ""),
                "synteny_score": k.get("synteny_score", "") if k else "",
                "mibig_quality": mibig_meta(a.mibig, k["mibig_accession"])[0] if k else "",
                "mibig_status": mibig_meta(a.mibig, k["mibig_accession"])[1] if k else "",
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
