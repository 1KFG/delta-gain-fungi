-- One row per species.parquet row (left joins only). Bound parameters: $species $busco $asm $prior $seqdb.
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
  p.protein_file_ok AS protein_file_ok,
  sd.proteins_path  AS bfd_proteins_path
FROM read_parquet($species) s
LEFT JOIN read_parquet($busco) b ON b.ASMID = s.ASMID
LEFT JOIN read_parquet($asm) a ON a.ASMID = s.ASMID
LEFT JOIN read_csv($prior, delim='\t', header=true, all_varchar=true) p
  ON p.accession = regexp_extract(s.ASMID, '^(GC[AF]_[0-9]+\.[0-9]+)', 1)
LEFT JOIN read_csv($seqdb, delim='\t', header=false, all_varchar=true,
                   columns={'locustag': 'VARCHAR', 'out_name': 'VARCHAR',
                            'proteins_path': 'VARCHAR', 'cds_path': 'VARCHAR'}) sd
  ON sd.locustag = s.LOCUSTAG
