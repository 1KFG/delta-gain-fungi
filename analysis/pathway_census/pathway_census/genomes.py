"""Genome table: taxonomy, quality and prior-annotation status for every BFD genome.

Real schemas (DESCRIBE, 2026-11-20):
  species.parquet (24,363 rows): ASMID, LOCUSTAG, PHYLUM, SUBPHYLUM, CLASS, ORDER,
      FAMILY, GENUS, SPECIES, BUSCO_LINEAGE, SUBCLASS, NCBI_TAXONID, STRAIN
  busco_genome.parquet (23,582 rows): ASMID, complete_pct, single_pct, duplicated_pct,
      fragmented_pct, missing_pct, n_markers, lineage
  asm_stats.parquet (23,584 rows): ASMID, SPECIES, STRAIN, contig_count, total_length_bp,
      ..., N50_bp, ..., annotation_status
  samples.parquet (24,363 rows): ASMID, SPECIES_IN, STRAIN, BIOPROJECT, NCBI_TAXONID, ... LOCUSTAG
  prior_annotation_status.tsv: accession, folder, annotation_provider, annotation_release,
      protein_coding, protein_file_ok, status, assembly_release, source_database

Join key: ASMID in the BFD tables (e.g. "GCF_010015735.1_Aaoar1") joins on ASMID between
BFD tables. To the prior-annotation TSV it joins on the accession prefix of ASMID,
regexp '^(GC[AF]_[0-9]+\\.[0-9]+)' = TSV `accession`. (The TSV `folder` column equals ASMID
for only some genomes, so `folder` is not used.)
The TSV has more rows (27,401) than species.parquet; unmatched TSV rows are dropped.

Column semantics (decided here):
  BFD file layout: seqdb/genomes.tsv (LOCUSTAG, out-name, proteins path, cds path) has 23,412
  rows; each genome has genome_annotation/<out>/predict_results/<out>.{proteins,scaffolds}.fa.gz.
  has_proteins = LOCUSTAG in seqdb/genomes.tsv and the BFD proteins.fa.gz exists on disk.
  downloaded   = BFD genome FASTA (<out>.scaffolds.fa.gz) exists on disk, i.e. searchable by tblastn.
  has_deposited_proteins = TSV protein_file_ok == 1 (prior-annotation protein file; <NA> if no TSV row).
  prior_folder_state = present | not_downloaded | unknown, from the TSV status
                 (not_downloaded: no_prior_annotation_not_downloaded, annotated_NOT_downloaded).
  unmatched TSV status -> prior_status "unknown".
  genus_source: genus_field | parsed_species | unresolved. Unresolved rows get a non-null label
                 "<first token> sp." (e.g. "Ascomycota sp.") or "unresolved_<phylum>".
"""
import logging
import os
import re
from pathlib import Path

import duckdb
import pandas as pd

log = logging.getLogger(__name__)

NOT_DOWNLOADED = {"no_prior_annotation_not_downloaded", "annotated_NOT_downloaded"}
_SQL = Path(__file__).resolve().parent.parent / "sql" / "genomes.sql"
# Higher-rank suffixes: "Xaceae sp." is a family, not a genus.
_RANK = re.compile(r"(mycota|mycotina|mycetes|mycetidae|ales|aceae|ineae|idia|ae)$")


def parse_genus(species):
    """Genus from a SPECIES string when GENUS is empty, or None."""
    if not species:
        return None
    toks = [t for t in species.split() if t not in ("aff.", "cf.")]
    if not toks or not re.fullmatch(r"[A-Z][a-z]{2,}", toks[0]):
        return None
    if "sp." in toks:
        i = toks.index("sp.")
        if i + 1 < len(toks) and re.fullmatch(r"[A-Z][a-z]{2,}", toks[i + 1]):
            return toks[i + 1]
        first = toks[0]
        # "Candida sp." -> Candida; "Ascomycota sp." / "Xaceae sp." -> None
        return None if _RANK.search(first) else first
    first = toks[0]
    return first if first[:1].isupper() and not _RANK.search(first) else None


def unresolved_label(species, phylum):
    """Non-null genus label for rows whose genus cannot be resolved."""
    toks = [t for t in (species or "").split() if t not in ("aff.", "cf.")]
    if toks:
        return f"{toks[0]} sp."
    return f"unresolved_{phylum or 'unknown'}"


def _exists_many(paths, workers=32):
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(workers) as ex:
        return list(ex.map(lambda q: isinstance(q, str) and os.path.exists(q), paths))


def load_genome_table(tables_dir, prior_status_tsv, seqdb_genomes_tsv=None):
    """One row per species.parquet row.

    has_proteins / downloaded come from BFD files, not from the prior-annotation TSV:
      has_proteins = LOCUSTAG is in seqdb/genomes.tsv and its BFD proteins.fa.gz exists on disk.
      downloaded   = the BFD genome FASTA (<out>.scaffolds.fa.gz, same predict_results
                     folder as the proteins file) exists on disk.
    has_deposited_proteins and prior_folder_state come from the prior-annotation TSV.
    """
    tables_dir = Path(tables_dir)
    if seqdb_genomes_tsv is None:
        seqdb_genomes_tsv = tables_dir.parent / "seqdb" / "genomes.tsv"
    sql = _SQL.read_text()
    df = duckdb.connect().execute(sql, {
        "species": str(tables_dir / "species.parquet"),
        "busco": str(tables_dir / "busco_genome.parquet"),
        "asm": str(tables_dir / "asm_stats.parquet"),
        "prior": str(prior_status_tsv),
        "seqdb": str(seqdb_genomes_tsv),
    }).df()
    n_species = duckdb.connect().execute(
        "SELECT count(*) FROM read_parquet($p)", {"p": str(tables_dir / "species.parquet")}
    ).fetchone()[0]
    assert len(df) == n_species, f"row count {len(df)} != species rows {n_species}"

    parsed = df["genus_field"].isna()
    df["genus"] = df["genus_field"].astype(object)
    df.loc[parsed, "genus"] = df.loc[parsed, "species"].map(parse_genus)
    df["genus_source"] = "genus_field"
    df.loc[parsed & df["genus"].notna(), "genus_source"] = "parsed_species"
    unres = df["genus"].isna()
    df.loc[unres, "genus"] = [unresolved_label(sp, ph) for sp, ph in
                              zip(df.loc[unres, "species"], df.loc[unres, "phylum"])]
    df.loc[unres, "genus_source"] = "unresolved"
    log.info("genus: %s", df["genus_source"].value_counts().to_dict())

    df["prior_status"] = df["prior_status"].fillna("unknown")
    df["prior_folder_state"] = [
        "unknown" if st == "unknown" else ("not_downloaded" if st in NOT_DOWNLOADED else "present")
        for st in df["prior_status"]]
    df["has_deposited_proteins"] = pd.array(
        [None if pd.isna(v) else (v == "1") for v in df["protein_file_ok"]], dtype="boolean")
    prot = df["bfd_proteins_path"].where(df["bfd_proteins_path"].notna(), None)
    scaf = [q.replace(".proteins.fa.gz", ".scaffolds.fa.gz") if isinstance(q, str) else None
            for q in prot]
    df["has_proteins"] = _exists_many(list(prot))
    df["downloaded"] = _exists_many(scaf)
    cols = ["locustag", "species", "phylum", "class_", "genus", "genus_source", "prior_status",
            "prior_folder_state", "annotation_provider", "busco_complete", "n50", "downloaded",
            "has_proteins", "has_deposited_proteins"]
    return df[cols]


def main(argv=None):
    """Write the genome table as TSV compressed with the zstd CLI."""
    import argparse
    import subprocess

    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", default="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables")
    ap.add_argument("--prior", default="/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/prior_annotation_status.tsv")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent.parent / "results" / "genomes.tsv.zst"))
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    df = load_genome_table(a.tables, a.prior)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    tsv = df.to_csv(sep="\t", index=False).encode()
    with open(a.out, "wb") as fh:
        subprocess.run(["zstd", "-q", "-T0", "-c"], input=tsv, stdout=fh, check=True)
    print(df["prior_status"].value_counts().to_string())


if __name__ == "__main__":
    main()
