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
  downloaded   = prior_status not in {no_prior_annotation_not_downloaded, annotated_NOT_downloaded}
  has_proteins = TSV protein_file_ok == 1 (a prior-annotation protein file was found).
                 It does NOT say whether BFD itself has re-annotated proteins.
  unmatched TSV status -> prior_status "unknown", downloaded/has_proteins <NA>.
"""
import logging
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


def load_genome_table(tables_dir, prior_status_tsv):
    tables_dir = Path(tables_dir)
    sql = _SQL.read_text().format(
        species=tables_dir / "species.parquet",
        busco=tables_dir / "busco_genome.parquet",
        asm=tables_dir / "asm_stats.parquet",
        prior=Path(prior_status_tsv),
    )
    df = duckdb.connect().sql(sql).df()
    parsed = df["genus_field"].isna()
    df["genus"] = df["genus_field"]
    df.loc[parsed, "genus"] = df.loc[parsed, "species"].map(parse_genus)
    df["genus"] = df["genus"].astype(object).where(df["genus"].notna(), None)
    df["genus_source"] = None
    df.loc[~parsed, "genus_source"] = "genus_field"
    df.loc[parsed & df["genus"].notna(), "genus_source"] = "parsed_species"
    df["genus_source"] = df["genus_source"].astype(object)
    log.info("genus parsed from species for %d of %d rows (%d unresolved)",
             int((df["genus_source"] == "parsed_species").sum()), len(df),
             int(df["genus"].isna().sum()))
    df["prior_status"] = df["prior_status"].fillna("unknown")
    known = df["prior_status"] != "unknown"
    df["downloaded"] = pd.array(
        [None if not k else (s not in NOT_DOWNLOADED) for k, s in zip(known, df["prior_status"])],
        dtype="boolean")
    df["has_proteins"] = pd.array(
        [None if pd.isna(v) else (v == "1") for v in df["protein_file_ok"]], dtype="boolean")
    cols = ["locustag", "species", "phylum", "class_", "genus", "genus_source", "prior_status",
            "annotation_provider", "busco_complete", "n50", "downloaded", "has_proteins"]
    return df[cols]


def main(argv=None):
    """Write the genome table as TSV compressed with the zstd CLI."""
    import argparse
    import subprocess

    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", default="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables")
    ap.add_argument("--prior", default="/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/prior_annotation_status.tsv")
    ap.add_argument("--out", default="results/genomes.tsv.zst")
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
