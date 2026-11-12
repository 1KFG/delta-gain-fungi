#!/usr/bin/env python3
"""Tables behind bfd_findings/REPORT.md: three observations made while building the UPO shortlist.

  report_tables.py --placed upo_place_bfd_all.tsv.zst --models upo_models.tsv \
      --genome-qc genome_qc.tsv --tax genome_tax.tsv --fasta all_pf01328.faa \
      --genomes-dir Fungi_BFD_runs/input_clean_genomes --threads 8 --out bfd_findings/out

1. Gene models. Every PF01328 protein gets a completeness category from BFD's own flags
   (`has_start_codon`, `has_stop_codon`): complete, no_start, no_stop, neither, flags_null.
   Rates are tabulated by phylum and set against the genome-wide rate of the same genomes, so an
   annotation problem specific to UPOs can be told apart from one affecting whole genomes. For
   Mucoromycota and the Mucorales basal clade, each model's distance to the nearest contig end is
   computed from the input assembly: a model truncated at a contig end is an assembly break, one
   truncated mid-contig is an annotation (prediction) error.
2. Glomeromycete clade: every protein placed in it, with phylum and the `long_pendant` flag.
3. Ascomycete proteins placed in Basidiomycetous 1st: every one, with taxonomy, copies per genome,
   placement statistics and the sequence; class-level carrier rates and a pendant-length comparison
   with the clade's basidiomycete members.

Standard library only. Writes TSVs, one FASTA and summary.md (markdown tables for the report).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

CATS = ["complete", "no_start", "no_stop", "neither", "flags_null"]
MUC_CLADE = "Mucoromycetous basal"
B1 = "Basidiomycetous 1st"
END_NEAR = 1000  # bp: a model within this distance of a contig end counts as "at contig end"


# ---------------------------------------------------------------- io
def text(path):
    if path.endswith(".zst"):
        return io.StringIO(subprocess.run(["zstd", "-dcq", path], capture_output=True, text=True,
                                          check=True).stdout)
    return open(path, newline="")


def tsv(path):
    return list(csv.DictReader((l for l in text(path) if not l.startswith("#")), delimiter="\t"))


def write_tsv(path, rows, cols=None):
    cols = cols or (list(rows[0]) if rows else ["empty"])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def pid(x):
    return x.split()[0]


# ---------------------------------------------------------------- logic (unit-tested)
def category(start, stop):
    """BFD flags -> completeness category. Flags are the strings TRUE/FALSE, or empty (NULL)."""
    s, t = (start or "").upper(), (stop or "").upper()
    if s not in ("TRUE", "FALSE") or t not in ("TRUE", "FALSE"):
        return "flags_null"
    return {("TRUE", "TRUE"): "complete", ("FALSE", "TRUE"): "no_start",
            ("TRUE", "FALSE"): "no_stop", ("FALSE", "FALSE"): "neither"}[(s, t)]


def end_distance(start, end, contig_len):
    """Distance (bp) from the model to the nearer contig end; None if the contig is unknown.
    Coordinates are 1-based inclusive as in the BFD tables."""
    if not contig_len:
        return None
    return max(0, min(int(start) - 1, int(contig_len) - int(end)))


def pct(k, n):
    return f"{100 * k / n:.1f}" if n else "NA"


def quartiles(vals):
    v = sorted(vals)
    if not v:
        return ("NA", "NA", "NA")
    if len(v) < 4:
        m = statistics.median(v)
        return (f"{v[0]:.3f}", f"{m:.3f}", f"{v[-1]:.3f}")
    q = statistics.quantiles(v, n=4)
    return tuple(f"{x:.3f}" for x in q)


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


# ---------------------------------------------------------------- contig lengths
def contig_lengths(job):
    """(asm_fasta, wanted contigs) -> {contig: length} for the wanted contigs only."""
    path, wanted = job
    out, name, n = {}, None, 0
    try:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if name in wanted:
                        out[name] = n
                    name, n = line[1:].split()[0], 0
                else:
                    n += len(line.strip())
        if name in wanted:
            out[name] = n
    except OSError as e:
        print(f"warn: {path}: {e}", file=sys.stderr)
    return path, out


# ---------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for k in ("placed", "models", "genome-qc", "tax", "fasta", "genomes-dir", "out"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    md = ["<!-- written by bfd_findings/bin/report_tables.py; do not edit by hand -->"]

    placed = {pid(r["protein_id"]): r for r in tsv(a.placed)}
    models = {r["protein_id"]: r for r in tsv(a.models)}
    gqc = {r["LOCUSTAG"]: r for r in tsv(a.genome_qc)}
    tax = {r["LOCUSTAG"]: r for r in tsv(a.tax)}
    for p, m in models.items():
        m["category"] = category(m["has_start_codon"], m["has_stop_codon"])

    # ---------- 1. gene-model completeness: UPO vs genome-wide, by phylum
    upo_by = defaultdict(Counter)
    genomes_by = defaultdict(set)
    for p, m in models.items():
        ph = m.get("PHYLUM") or "?"
        upo_by[ph][m["category"]] += 1
        genomes_by[ph].add(p.split("_")[0])
    rows, mdrows = [], []
    for ph in sorted(upo_by, key=lambda k: -sum(upo_by[k].values())):
        c, n = upo_by[ph], sum(upo_by[ph].values())
        gw = Counter()
        for lt in genomes_by[ph]:
            g = gqc.get(lt, {})
            for k in ("n_models", "n_complete", "n_no_start", "n_no_stop", "n_neither",
                      "n_flags_null"):
                gw[k] += int(g.get(k) or 0)
        rows.append({"phylum": ph, "genomes": len(genomes_by[ph]), "upo_models": n,
                     **{f"upo_{k}": c[k] for k in CATS}, "upo_incomplete_pct": pct(n - c["complete"], n),
                     "genomewide_models": gw["n_models"],
                     "genomewide_incomplete_pct": pct(gw["n_models"] - gw["n_complete"], gw["n_models"]),
                     "genomewide_neither_pct": pct(gw["n_neither"], gw["n_models"])})
        mdrows.append([ph, len(genomes_by[ph]), n, pct(n - c["complete"], n), pct(c["neither"], n),
                       pct(gw["n_models"] - gw["n_complete"], gw["n_models"]),
                       pct(gw["n_neither"], gw["n_models"])])
    write_tsv(os.path.join(a.out, "gene_model_qc_by_phylum.tsv"), rows)
    md += ["## 1. Gene-model completeness by phylum",
           "UPO = the PF01328 proteins; genome-wide = every gene model in the same genomes.",
           md_table(["phylum", "genomes", "UPO models", "UPO incomplete %", "UPO neither %",
                     "genome-wide incomplete %", "genome-wide neither %"], mdrows)]

    # Mucoromycota + Mucorales basal clade, model by model, with contig-end distance
    muc = [p for p, m in models.items()
           if m.get("PHYLUM") == "Mucoromycota"
           or placed.get(p, {}).get("place_fig1_clade") == MUC_CLADE]
    want = defaultdict(set)
    for p in muc:
        want[models[p]["ASMID"]].add(models[p]["chrom"])
    jobs = [(os.path.join(a.genomes_dir, f"{asm}.fa.gz"), cs) for asm, cs in want.items() if asm]
    clen = {}
    with ProcessPoolExecutor(max_workers=a.threads) as ex:
        for path, d in ex.map(contig_lengths, jobs):
            asm = os.path.basename(path)[:-len(".fa.gz")]
            for c, n in d.items():
                clen[(asm, c)] = n
    muc_rows = []
    for p in sorted(muc):
        m, pl = models[p], placed.get(p, {})
        L = clen.get((m["ASMID"], m["chrom"]))
        d = end_distance(m["start"], m["end"], L)
        lt = p.split("_")[0]
        g = gqc.get(lt, {})
        muc_rows.append({
            "protein_id": p, "gene_id": m["gene_id"], "ASMID": m["ASMID"], "LOCUSTAG": lt,
            "species": m["SPECIES"], "phylum": m["PHYLUM"], "order": m["ORDER_"],
            "family": m["FAMILY"], "category": m["category"], "is_partial": m["is_partial"],
            "has_start_codon": m["has_start_codon"], "has_stop_codon": m["has_stop_codon"],
            "CDS_length": m["CDS_length"], "cds_mult3": "yes" if int(m["CDS_length"] or 0) % 3 == 0 else "no",
            "protein_length": pl.get("length", ""), "n_exons": m["n_exons"],
            "chrom": m["chrom"], "start": m["start"], "end": m["end"], "strand": m["strand"],
            "contig_length": L if L is not None else "",
            "dist_to_contig_end": d if d is not None else "",
            "at_contig_end": "" if d is None else ("yes" if d <= END_NEAR else "no"),
            "genome_contigs": m["contig_count"], "genome_N50": m["N50_bp"],
            "genome_models": g.get("n_models", ""),
            "genome_incomplete_pct": pct(int(g.get("n_models") or 0) - int(g.get("n_complete") or 0),
                                         int(g.get("n_models") or 0)),
            "fig1_clade": pl.get("place_fig1_clade", ""), "our_clade": pl.get("place_our_clade", ""),
        })
    write_tsv(os.path.join(a.out, "mucoromycota_upo_gene_models.tsv"), muc_rows)

    def tab_muc(sel, label):
        n = len(sel)
        c = Counter(r["category"] for r in sel)
        inc = [r for r in sel if r["category"] != "complete"]
        known = [r for r in inc if r["at_contig_end"]]
        at_end = sum(r["at_contig_end"] == "yes" for r in known)
        return [label, n, pct(n - c["complete"], n), c["no_start"], c["no_stop"], c["neither"],
                c["flags_null"], f"{at_end}/{len(known)}",
                sum(r["cds_mult3"] == "no" for r in inc)]
    clade_rows = [r for r in muc_rows if r["fig1_clade"] == MUC_CLADE]
    md += ["### Mucoromycota and the Mucorales basal clade",
           f"`at contig end`: incomplete models within {END_NEAR} bp of a contig end, of those whose "
           "contig was found in the input assembly.",
           md_table(["set", "models", "incomplete %", "no start", "no stop", "neither",
                     "flags null", "incomplete at contig end", "incomplete, CDS not x3"],
                    [tab_muc(muc_rows, "all Mucoromycota PF01328 (+ basal clade members)"),
                     tab_muc(clade_rows, "Mucorales basal clade members"),
                     tab_muc([r for r in muc_rows if r["phylum"] == "Mucoromycota"],
                             "Mucoromycota proteins")])]
    found = sum(1 for r in muc_rows if r["contig_length"] != "")
    md.append(f"Contig found for {found} of {len(muc_rows)} models.")

    # per-genome table for Mucoromycota
    per = defaultdict(list)
    for r in muc_rows:
        per[(r["ASMID"], r["LOCUSTAG"], r["species"], r["order"])].append(r)
    gen_rows = []
    for (asm, lt, sp, order), rs in sorted(per.items()):
        g = gqc.get(lt, {})
        n_m = int(g.get("n_models") or 0)
        gen_rows.append({"ASMID": asm, "LOCUSTAG": lt, "species": sp, "order": order,
                         "upo_models": len(rs),
                         "upo_incomplete": sum(r["category"] != "complete" for r in rs),
                         "upo_incomplete_at_contig_end": sum(r["at_contig_end"] == "yes" for r in rs
                                                             if r["category"] != "complete"),
                         "genome_models": n_m,
                         "genome_incomplete_pct": pct(n_m - int(g.get("n_complete") or 0), n_m),
                         "genome_neither_pct": pct(int(g.get("n_neither") or 0), n_m),
                         "genome_cds_not_mult3": g.get("n_cds_not_mult3", ""),
                         "genome_contigs": rs[0]["genome_contigs"], "genome_N50": rs[0]["genome_N50"]})
    write_tsv(os.path.join(a.out, "mucoromycota_genome_qc.tsv"), gen_rows)
    inc_g = [float(r["genome_incomplete_pct"]) for r in gen_rows if r["genome_incomplete_pct"] != "NA"]
    md.append(f"Mucoromycota genomes with a PF01328 protein: {len(gen_rows)}; genome-wide incomplete "
              f"gene models per genome, quartiles {quartiles(inc_g)} %.")

    # ---------- 2. glomeromycete clade
    glo = [r for r in placed.values() if r["place_fig1_clade"].startswith("Glomeromycetous")]
    glo_rows = [{"protein_id": pid(r["protein_id"]), "species": r["SPECIES"], "phylum": r["PHYLUM"],
                 "class": r["CLASS"], "lineage": r["LINEAGE_GROUP"],
                 "fig1_clade": r["place_fig1_clade"], "fig1_lwr": r["place_fig1_lwr"],
                 "pendant_length": r["place_pendant_length"], "flags": r["place_flags"],
                 "our_clade": r["place_our_clade"], "length": r["length"], "p180": r["p180_Acid"],
                 "gene_model": models.get(pid(r["protein_id"]), {}).get("category", "")}
                for r in glo]
    write_tsv(os.path.join(a.out, "glomeromycete_clade_members.tsv"), glo_rows)
    by = defaultdict(Counter)
    pend = defaultdict(list)
    for r in glo_rows:
        k = r["phylum"] if r["phylum"] != "Ascomycota" else f"Ascomycota: {r['class']}"
        by[k]["long_pendant" if r["flags"] else "short"] += 1
        pend[k].append(float(r["pendant_length"] or 0))
    md += ["## 2. Glomeromycete clade members",
           md_table(["phylum / class", "members", "long_pendant", "short pendant",
                     "pendant median"],
                    [[k, sum(c.values()), c["long_pendant"], c["short"],
                      f"{statistics.median(pend[k]):.2f}"]
                     for k, c in sorted(by.items(), key=lambda x: -sum(x[1].values()))])]

    # ---------- 3. ascomycetes in Basidiomycetous 1st
    b1 = [r for r in placed.values() if r["place_fig1_clade"] == B1]
    asco = [r for r in b1 if r["PHYLUM"] == "Ascomycota"]
    copies = Counter(pid(r["protein_id"]).split("_")[0] for r in asco)
    a_rows = []
    for r in asco:
        p = pid(r["protein_id"])
        lt = p.split("_")[0]
        t = tax.get(lt, {})
        a_rows.append({"protein_id": p, "species": r["SPECIES"], "class": t.get("CLASS", r["CLASS"]),
                       "order": t.get("ORDER_", ""), "family": t.get("FAMILY", ""),
                       "genus": t.get("GENUS", ""), "ASMID": t.get("ASMID", ""), "LOCUSTAG": lt,
                       "copies_in_genome": copies[lt], "fig1_lwr": r["place_fig1_lwr"],
                       "best_edge_lwr": r["place_best_edge_lwr"],
                       "pendant_length": r["place_pendant_length"], "flags": r["place_flags"],
                       "our_clade": r["place_our_clade"], "length": r["length"],
                       "p39": r["p39_Cys_proximal"], "p109": r["p109_E_of_EHD"],
                       "p111": r["p111_D_of_EHD"], "p180": r["p180_Acid"],
                       "old_nn_clade": r["clade"], "old_nearest_ref": r["nearest_ref"],
                       "gene_model": models.get(p, {}).get("category", "")})
    a_rows.sort(key=lambda r: (r["class"], r["order"], r["species"], r["protein_id"]))
    write_tsv(os.path.join(a.out, "asco_in_basidiomycetous_1st.tsv"), a_rows)
    want = {r["protein_id"] for r in a_rows}
    with open(os.path.join(a.out, "asco_in_basidiomycetous_1st.faa"), "w") as fo:
        keep = False
        for line in text(a.fasta):
            if line.startswith(">"):
                keep = pid(line[1:]) in want
            if keep:
                fo.write(line)

    genomes_in_class = Counter(t.get("CLASS") or "?" for t in tax.values() if t.get("PHYLUM") == "Ascomycota")
    carriers = defaultdict(set)
    for r in a_rows:
        carriers[r["class"] or "?"].add(r["LOCUSTAG"])
    cls_n = Counter(r["class"] or "?" for r in a_rows)
    md += ["## 3. Ascomycete proteins placed in Basidiomycetous 1st",
           f"{len(asco)} of {len(b1)} members are Ascomycota, in {len(copies)} genomes.",
           md_table(["class", "proteins", "genomes carrying", "BFD genomes in class", "carrier %"],
                    [[c, cls_n[c], len(carriers[c]), genomes_in_class[c],
                      pct(len(carriers[c]), genomes_in_class[c])]
                     for c in sorted(cls_n, key=lambda k: -cls_n[k])])]
    ords = Counter((r["class"], r["order"]) for r in a_rows)
    md += ["Top orders:", md_table(["class", "order", "proteins"],
                                   [[c, o, n] for (c, o), n in ords.most_common(15)])]
    cpg = Counter(copies.values())
    md.append("Copies per carrier genome: " + ", ".join(f"{k}: {cpg[k]}" for k in sorted(cpg)))
    comp = []
    for ph in ("Basidiomycota", "Ascomycota"):
        rs = [r for r in b1 if r["PHYLUM"] == ph]
        q = quartiles([float(r["place_pendant_length"] or 0) for r in rs])
        lw = quartiles([float(r["place_fig1_lwr"] or 0) for r in rs])
        e180 = Counter(r["p180_Acid"] for r in rs).most_common(4)
        comp.append([ph, len(rs), " / ".join(q), " / ".join(lw),
                     ", ".join(f"{k or '-'} {pct(v, len(rs))}%" for k, v in e180),
                     sum(1 for r in rs if r["place_flags"])])
    md += ["Members of the clade by phylum (quartiles Q1 / median / Q3):",
           md_table(["phylum", "members", "pendant length", "clade LWR", "residue at 180",
                     "long_pendant"], comp)]
    gm = Counter(r["gene_model"] for r in a_rows)
    md.append("Gene models of the ascomycete members: " + ", ".join(f"{k} {v}" for k, v in gm.most_common()))

    with open(os.path.join(a.out, "summary.md"), "w") as fh:
        fh.write("\n\n".join(md) + "\n")
    print(f"report tables -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
