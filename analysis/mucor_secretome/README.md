# *Mucor hiemalis* secretome: input set for a structure-based plant-interaction screen

**Rules written 2026-11-11, before the run.** The purpose is to produce the fungal side of an
AlphaFold3-style screen of secreted proteins against candidate plant genes, with the proteins
classified and ranked so the screen does not have to fold everything.

Code: `bin/extract_proteins.py`, `bin/build_secretome.py`; driver `slurm/secretome.sh`; tests
`test/` . The 5' completeness check reuses `../enzymes_UPO/bfd_findings/bin/check_codons.py`.

## Scope

Focal genome **GCA_977110905.1_gzMucHiem1** (*Mucor hiemalis*, LOCUSTAG `F00CDB58`, 10,434
proteins) — the only *M. hiemalis* genome in BFD. The candidate pool for protein families is every
**Mucoromycota** genome in BFD: 275 Mucorales, 14 Umbelopsidales, 5 Endogonales.

There is no within-species replication, so nothing here can separate *hiemalis*-specific biology
from assembly- or annotation-specific artefacts. Conservation is measured across species, and a
protein found only in this genome is reported as `genome_specific` rather than as novel.

## Why this framing

*M. hiemalis* is a soil saprotroph, not a plant pathogen, and Mucoromycota lack the expanded
effector families that the plant-pathogen screens were built around. Reading its secretome as an
effector arsenal would assume the answer. The tiers below put first the two questions that stand on
their own:

1. **What the plant perceives.** Mucorales cell walls are chitosan-rich — chitin deacetylated by
   CE4 enzymes — and plant LysM receptors (CERK1, LYK5, CEBiP) bind chitosan poorly. Perception and
   masking is a receptor–ligand question with a defined ligand, which is what a co-folding method
   is actually good at.
2. **What acts on plant material.** Secreted enzymes with plant polysaccharide, cutin or lipid
   substrates, including the peroxygenases this repository works on (PF01328).

Proteins with no domain at all are kept, but ranked below both, and labelled as a hypothesis.

## Secretion call

From BFD's own per-protein predictions (`signalp.signal_peptide`, `tmhmm`, `predgpi`, `targetP`):

| call | rule |
|---|---|
| signal peptide | SignalP probability >= 0.7 (the pool is exported at >= 0.5) |
| `membrane` | a TMHMM helix starting more than 5 aa after the cleavage site |
| `gpi_anchored` | a PredGPI feature; kept in the screen set, labelled |
| `mitochondrial` | TargetP mTP at probability >= 0.7; dropped |
| `secreted` | signal peptide, none of the above |

A helix lying inside the signal peptide is the signal peptide, not a transmembrane domain, so it
does not disqualify a protein.

## 5' completeness

Taken from the **sequence**, never from BFD's `has_start_codon` flag. Those flags are wrong in
these genomes: the Mucorales rate of models flagged incomplete is about 54 %, and 99.3 % of
flagged-missing starts are in fact an ATG in the assembly (`../enzymes_UPO/bfd_findings/REPORT.md`
§1; fixed in Fungi_BFD fb1d7cd, but `tables/` predates the fix). This matters more here than
anywhere else in the repository, because a model that really has lost its 5' end has lost its
signal peptide and should not be in a secretome at all.

Only proteins whose CDS begins with ATG in the assembly enter the screen set.

## Tiers

A protein takes the first tier it qualifies for.

| tier | what it is | rule |
|---|---|---|
| T1 | chitin / chitosan | CAZy CE4, GH18, GH20, GH46, GH75, or a chitin-binding module (CBM12/14/18/50); or Pfam PF01522 (polysaccharide deacetylase), PF01476 (LysM) |
| T2 | plant-facing enzyme | dbCAN substrate naming a plant polysaccharide (cellulose, xylan, pectin, xyloglucan, mannan, arabinan, galactan, lignin, cutin, starch, ...); or cutinase, lipase, peroxygenase/peroxidase Pfam |
| T3 | effector-like unknown | no Pfam domain, <= 300 aa, >= 4 cysteines and >= 2 % cysteine |
| T4 | other secreted | everything else |

Ranking inside a tier is by the number of Mucoromycota species sharing the protein's family, most
widely shared first, then by protein ID. Wide conservation is treated as evidence that the protein
is real and expressed, not as evidence that it interacts with a plant.

## Protein families

MMseqs2 `easy-cluster` over the whole Mucoromycota candidate secretome at 50 % identity and 80 %
coverage of the shorter sequence. Each focal protein reports its cluster's size, and the species,
genera and orders holding a member. `family_presence_absence.tsv` is the per-species matrix for
the screen-set clusters.

These clusters are sequence-similarity families, not orthogroups. A cluster spanning several orders
means the family is old and widespread; it does not establish one-to-one orthology.

## Outputs (`out/`)

| file | contents |
|---|---|
| `focal_secretome.tsv` | every *M. hiemalis* protein with a signal peptide at p >= 0.5, classified, tiered, with family conservation and the sequence-level start/stop check |
| `screen_set.tsv`, `screen_set.faa` | the ranked proteins to fold: secreted or GPI-anchored, real ATG start |
| `family_presence_absence.tsv` | each screen-set cluster across the Mucoromycota species |
| `codon_check.tsv` | BFD's start/stop flags against the assembly, for the whole candidate pool |
| `mucoromycota_genomes.tsv` | the 294 genomes, with taxonomy |
| `SUMMARY.md` | counts, written by the script |

## What this does not do

- It does not choose the plant side. That set should be a defined receptor/immunity panel, and the
  pipeline should first be benchmarked on pairs with solved structures (flg22/FLS2–BAK1,
  chitin/CERK1–LYK5, AVR4/Cf-4, Avr2/Rcr3). If those are not recovered at the chosen interface
  threshold, a *Mucor* screen cannot be interpreted.
- It does not predict interactions, and a tier is a reason to fold a protein, not a claim about it.
- It has no expression evidence. A protein secreted in culture is not necessarily made in contact
  with a root.
