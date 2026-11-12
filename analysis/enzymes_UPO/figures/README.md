# UPO clade summary: figures and Table T1

Built by `slurm/figures.sh` (on the HPCC): `bin/prep_figures.py` writes `data/*.tsv`,
`bin/plot_figures.R` (ggtree, ggalluvial, ggplot2 via pinned pixi environment) writes `out/*.{pdf,png}`.
Tests: `python3 -m unittest discover -s figures/test`.

| file | shows |
|---|---|
| `out/fig1_reference_tree` | The 187-leaf reference tree, rooted on Mucorales; Zamocky 2022 lineages as sectors (11 clades folded into 8 colour groups); our Clades I-V as arcs with SH-aLRT/UFBoot; the 19 characterized enzymes named |
| `out/fig3_old_vs_placement` | 97,520 BFD proteins: old nearest-neighbour call -> placement call, flows coloured by the old call's confidence tier |
| `out/fig4_matched_identity` | Matched-identity validation of both callers on the 19 characterized enzymes: precision with exact 95 % CI, and abstention rate |
| `data/table_T1.md` | Per clade: references, support, leave-one-out, BFD members, catalytic core, Glu180, length, measured T1 activity range |

## Read Table T1 with these caveats

- **Glomeromycete clade:** 784 BFD members but only 121 are Glomeromycota; 585 carry `long_pendant`.
  Treat its count as an upper bound inflated by divergent proteins attaching to a small isolated
  clade (`../CANDIDATES.md`).
- **Basidiomycetous 1st:** about 2,000 of 6,378 members are ascomycetes, *not* long-pendant. Unresolved.
- **Family I (short):** its node support is low (SH-aLRT/UFBoot 79/65); Clades III-V inside it are
  each better supported than Family I as a whole.
- **Clade III:** 3 reference leaves; only 60 % of its BFD calls reach LWR >= 0.95.
- T1 activity ranges are from the 19 characterized enzymes only (n per clade shown); Clade IV is n = 1.
