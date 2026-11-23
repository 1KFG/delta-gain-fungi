# Small result tables kept in git (the full per-genome outputs stay in `results/`, which is git-ignored)

- `ebc_regions_v3_35genomes.tsv.gz`: one row per antiSMASH region, 35 genomes, 2,278 regions, with the best of the top-3 MIBiG hits (the table sent to the EBC agent as message 014).
- `timings.tsv`: antiSMASH wall seconds (id, genome folder, seconds, CPUs). A genome that was run more than once has more than one row.
- `variance_by_species_32genomes.csv`: regions and `known` regions per genome within a species (32 genomes, before the three Penicillium genomes were added).
- `contig_filter_by_genome.csv`: contig and base fractions below 1, 5, 10 and 20 kb per genome, with wall time.

Method and caveats: see `../FINDINGS.md`. Per-genome gz tables (regions, cds, known, known_pairs, run.json) are on HPCC in `results/` (about 85 MB) and are not archived here.
