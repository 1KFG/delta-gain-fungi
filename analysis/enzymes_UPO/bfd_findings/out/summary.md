<!-- written by bfd_findings/bin/report_tables.py; do not edit by hand -->

## 1. Gene-model completeness by phylum

UPO = the PF01328 proteins; genome-wide = every gene model in the same genomes.

| phylum | genomes | UPO models | UPO incomplete % | UPO neither % | genome-wide incomplete % | genome-wide neither % |
|---|---|---|---|---|---|---|
| Ascomycota | 12061 | 78382 | 9.4 | 0.5 | 12.7 | 0.5 |
| Basidiomycota | 2677 | 18327 | 7.9 | 0.2 | 14.8 | 1.6 |
| Mucoromycota | 252 | 370 | 83.0 | 44.3 | 48.9 | 10.9 |
| Chytridiomycota | 56 | 188 | 22.9 | 1.6 | 27.7 | 3.4 |
| Glomeromycota | 57 | 124 | 75.0 | 44.4 | 60.8 | 17.7 |
| Basidiobolomycota | 18 | 86 | 11.6 | 0.0 | 20.6 | 5.6 |
| Zoopagomycota | 10 | 30 | 30.0 | 0.0 | 27.7 | 2.4 |
| Monoblepharomycota | 2 | 4 | 0.0 | 0.0 | 13.4 | 0.7 |
| ? | 1 | 4 | 0.0 | 0.0 | 21.1 | 1.0 |
| Aphelidiomycota | 2 | 3 | 33.3 | 0.0 | 20.7 | 1.2 |
| Calcarisporiellomycota | 1 | 2 | 0.0 | 0.0 | 23.1 | 1.2 |

### Mucoromycota and the Mucorales basal clade

`at contig end`: incomplete models within 1000 bp of a contig end, of those whose contig was found in the input assembly.

| set | models | incomplete % | no start | no stop | neither | flags null | incomplete at contig end | incomplete, CDS not x3 |
|---|---|---|---|---|---|---|---|---|
| all Mucoromycota PF01328 (+ basal clade members) | 371 | 83.0 | 47 | 97 | 164 | 0 | 19/308 | 0 |
| Mucorales basal clade members | 315 | 92.7 | 39 | 94 | 159 | 0 | 19/292 | 0 |
| Mucoromycota proteins | 370 | 83.0 | 46 | 97 | 164 | 0 | 19/307 | 0 |

Contig found for 371 of 371 models.

Mucoromycota genomes with a PF01328 protein: 253; genome-wide incomplete gene models per genome, quartiles ('30.150', '54.000', '56.750') %.

## 2. Glomeromycete clade members

| phylum / class | members | long_pendant | short pendant | pendant median |
|---|---|---|---|---|
| Ascomycota: Dothideomycetes | 360 | 341 | 19 | 1.17 |
| Glomeromycota | 121 | 0 | 121 | 0.04 |
| Ascomycota: Sordariomycetes | 86 | 71 | 15 | 1.19 |
| Ascomycota: Eurotiomycetes | 66 | 65 | 1 | 1.22 |
| Basidiomycota | 62 | 61 | 1 | 1.36 |
| Basidiobolomycota | 29 | 12 | 17 | 1.00 |
| Ascomycota: Orbiliomycetes | 22 | 3 | 19 | 0.97 |
| Ascomycota: Lecanoromycetes | 20 | 20 | 0 | 1.17 |
| Ascomycota: Leotiomycetes | 17 | 12 | 5 | 1.13 |
| Ascomycota: Geoglossomycetes | 1 | 0 | 1 | 0.85 |

## 3. Ascomycete proteins placed in Basidiomycetous 1st

2265 of 6378 members are Ascomycota, in 2190 genomes.

| class | proteins | genomes carrying | BFD genomes in class | carrier % |
|---|---|---|---|---|
| Dothideomycetes | 1364 | 1298 | 2731 | 47.5 |
| Eurotiomycetes | 680 | 675 | 3220 | 21.0 |
| Leotiomycetes | 77 | 75 | 435 | 17.2 |
| Orbiliomycetes | 69 | 67 | 67 | 100.0 |
| Sordariomycetes | 65 | 65 | 6016 | 1.1 |
| ? | 3 | 3 | 48 | 6.2 |
| Arthoniomycetes | 3 | 3 | 13 | 23.1 |
| Lecanoromycetes | 2 | 2 | 479 | 0.4 |
| Geoglossomycetes | 1 | 1 | 12 | 8.3 |
| Lipomycetes | 1 | 1 | 59 | 1.7 |

Top orders:

| class | order | proteins |
|---|---|---|
| Eurotiomycetes | Eurotiales | 670 |
| Dothideomycetes | Pleosporales | 446 |
| Dothideomycetes | Mycosphaerellales | 401 |
| Dothideomycetes | Botryosphaeriales | 310 |
| Dothideomycetes | Venturiales | 111 |
| Dothideomycetes | Cladosporiales | 82 |
| Orbiliomycetes | Orbiliales | 69 |
| Sordariomycetes | Hypocreales | 56 |
| Leotiomycetes | Helotiales | 40 |
| Leotiomycetes | Thelebolales | 35 |
| Eurotiomycetes | Chaetothyriales | 10 |
| Dothideomycetes |  | 6 |
| Sordariomycetes | Amphisphaeriales | 4 |
|  |  | 3 |
| Dothideomycetes | Myriangiales | 3 |

Copies per carrier genome: 1: 2123, 2: 60, 3: 6, 4: 1

Members of the clade by phylum (quartiles Q1 / median / Q3):

| phylum | members | pendant length | clade LWR | residue at 180 | long_pendant |
|---|---|---|---|---|---|
| Basidiomycota | 4097 | 0.171 / 0.279 / 0.377 | 0.995 / 1.000 / 1.000 | N 87.4%, - 4.7%, K 3.4%, S 3.4% | 14 |
| Ascomycota | 2265 | 0.170 / 0.229 / 0.401 | 0.694 / 1.000 / 1.000 | N 87.4%, T 6.2%, - 5.2%, G 0.8% | 0 |

Gene models of the ascomycete members: complete 2186, no_stop 52, no_start 24, neither 3
