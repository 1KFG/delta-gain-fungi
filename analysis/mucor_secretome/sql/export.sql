-- Annotation tables for the Mucoromycota secretome build (bin/build_secretome.py).
-- Run from Fungi_BFD_runs/tables with duckdb 1.1.3; ${OUT} is substituted by the driver.
-- Scope: every Mucoromycota genome in BFD (Mucorales, Umbelopsidales, Endogonales).
-- LOCUSTAG is the protein-ID prefix before the first '_'; signalp/tmhmm/cazy carry it as
-- species_prefix, pfam does not.
.mode csv
.separator "\t"

create temp table g as
  select LOCUSTAG, ASMID, SPECIES, STRAIN, "ORDER" as ORD, FAMILY, GENUS
  from read_parquet('species.parquet') where PHYLUM = 'Mucoromycota';
copy (select * from g) to '${OUT}/genomes.tsv' (header, delimiter '\t');

-- Candidate secreted proteins: any signal peptide at p >= 0.5. The 0.7 cut and the
-- TM/GPI/mitochondrial filters are applied in build_secretome.py, not here, so the thresholds
-- stay visible in one place.
create temp table sp as
  select s.protein_id, s.species_prefix, s.peptide_start, s.peptide_end, s.probability
  from read_parquet('signalp.signal_peptide.parquet') s
  where s.species_prefix in (select LOCUSTAG from g) and s.probability >= 0.5;
copy (select * from sp) to '${OUT}/signalp.tsv' (header, delimiter '\t');

copy (select protein_id, len, ExpAA, First60, PredHel, Topology from read_parquet('tmhmm.parquet')
      where species_prefix in (select LOCUSTAG from g))
  to '${OUT}/tmhmm.tsv' (header, delimiter '\t');

copy (select distinct protein_id, feature, start, "end" from read_parquet('predgpi.parquet')
      where species_prefix in (select LOCUSTAG from g))
  to '${OUT}/gpi.tsv' (header, delimiter '\t');

copy (select protein_id, prediction, probability from read_parquet('targetP.parquet')
      where species_prefix in (select LOCUSTAG from g))
  to '${OUT}/targetp.tsv' (header, delimiter '\t');

copy (select protein_id, pfam_acc, pfam_id, domain_i_evalue, ali_from, ali_to
      from read_parquet('pfam.parquet')
      where split_part(protein_id, '_', 1) in (select LOCUSTAG from g) and domain_i_evalue < 1e-5)
  to '${OUT}/pfam.tsv' (header, delimiter '\t');

copy (select protein_id, EC, cazyme_fam, sub_fam, substrate, toolcount
      from read_parquet('cazy.overview.parquet')
      where species_prefix in (select LOCUSTAG from g))
  to '${OUT}/cazy.tsv' (header, delimiter '\t');

copy (select p.protein_id, p.transcript_id, p.length
      from read_parquet('gene_proteins.parquet') p
      where split_part(p.protein_id, '_', 1) in (select LOCUSTAG from g))
  to '${OUT}/lengths.tsv' (header, delimiter '\t');

-- Gene-model flags, and the CDS segments the sequence-level completeness check needs.
-- The flags are known to be wrong before Fungi_BFD fb1d7cd (enzymes_UPO/bfd_findings/REPORT.md);
-- they are exported only as the input check_codons.py compares against.
copy (select 'mucoromycota_secretome' as "set", t.transcript_id, g.ASMID,
             t.has_start_codon, t.has_stop_codon, t.is_partial
      from read_parquet('gene_transcripts.parquet') t
      join g on g.LOCUSTAG = split_part(t.transcript_id, '_', 1)
      where t.transcript_id in (select protein_id from sp))
  to '${OUT}/codon_models.tsv' (header, delimiter '\t');

copy (select c.transcript_id, c.chrom, c.start, c."end", c.strand
      from read_parquet('gene_CDS.parquet') c
      where c.transcript_id in (select protein_id from sp))
  to '${OUT}/codon_cds.tsv' (header, delimiter '\t');
