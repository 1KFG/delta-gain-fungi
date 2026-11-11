# BFD coverage check (2026-11-08)

Verification that the site + score layer covers every PF01328 protein in BFD, rather than
assuming it does.

## Method

BFD's own annotation is the independent reference. Protein IDs were taken straight from
`Fungi_BFD_runs/tables/pfam.parquet` where `pfam_acc LIKE 'PF01328%'` and compared as sets
against `data/bfd_full/upo_sites_bfd_scored.tsv.zst`.

## Result

| quantity | value |
|---|---|
| PF01328 rows in pfam.parquet | 98,159 |
| **distinct** PF01328 proteins | **97,520** |
| proteins in our scored table | **97,520** |
| in BFD but not ours | **0** |
| in ours but not BFD | **0** |
| genomes with PF01328 | 15,137 (both) |

The 98,159 figure is a row count: a protein with several PF01328 domain hits contributes several
rows. Distinct proteins match exactly, in both directions.

Total genomes in BFD's pfam table are 23,410. The 8,273 genomes absent here carry no PF01328
protein at all, so they are correctly absent rather than skipped.

## Reproduce

```bash
module load duckdb
duckdb -c "copy (select distinct protein_id from read_parquet(
  '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/pfam.parquet')
  where pfam_acc like 'PF01328%') to '/tmp/bfd_pf.txt' (header false)"
sort /tmp/bfd_pf.txt > /tmp/bfd_pf.s
zstd -dcq data/bfd_full/upo_sites_bfd_scored.tsv.zst | grep -v '^#' | tail -n +2 \
  | cut -f1 | sort > /tmp/ours.s
comm -23 /tmp/bfd_pf.s /tmp/ours.s | wc -l   # expect 0
comm -13 /tmp/bfd_pf.s /tmp/ours.s | wc -l   # expect 0
```

Note the pfam column is `pfam_acc`, not `pfam_id`; filtering on `pfam_id` returns zero rows.
