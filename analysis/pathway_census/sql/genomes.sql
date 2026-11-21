-- One row per species.parquet row (left joins only). Parameters substituted by genomes.py.
SELECT
  s.LOCUSTAG  AS locustag,
  s.ASMID     AS asmid,
  s.SPECIES   AS species,
  s.PHYLUM    AS phylum,
  s."CLASS"   AS class_,
  NULLIF(trim(s.GENUS), '') AS genus_field,
  p.status    AS prior_status,
  NULLIF(p.annotation_provider, '') AS annotation_provider,
  b.complete_pct AS busco_complete,
  a.N50_bp    AS n50,
  p.protein_file_ok AS protein_file_ok
FROM read_parquet('{species}') s
LEFT JOIN read_parquet('{busco}') b ON b.ASMID = s.ASMID
LEFT JOIN read_parquet('{asm}') a ON a.ASMID = s.ASMID
LEFT JOIN read_csv('{prior}', delim='\t', header=true, all_varchar=true) p
  ON p.accession = regexp_extract(s.ASMID, '^(GC[AF]_[0-9]+\.[0-9]+)', 1)
