# Calibration set: published UPO activity on aromatic rings and PAHs

Status: literature survey, 2026-11-07. Source: PubMed. Nothing computed here. Every number below is
quoted from a paper; none is measured by us. DOIs are given for each claim.

Purpose: decide whether published enzyme-substrate data is sufficient to calibrate a computational screen
that predicts the substrate band (fatty acid / 1 ring / 2-3 ring / 4-5 ring / none) of a PF01328 protein.

## 1. Headline answer

A calibration set exists and it is better than expected, but it is **small**. Two enzymes have a full PAH
panel measured side by side. Two or three more have partial data. That is enough to fix the *direction* of
the structure-activity relationship. It is not enough to fit a quantitative predictor.

The relationship we proposed to test is already published, with numbers:
**heme access channel width controls which PAHs are accepted, and the discriminating substrate dimension is
the minimum width, not the maximum.**

## 2. The core dataset — Karich et al. 2017

Karich A, Ullrich R, Scheibner K, Hofrichter M. "Fungal Unspecific Peroxygenases Oxidize the Majority of
Organic EPA Priority Pollutants." Front Microbiol 8:1463.
[DOI](https://doi.org/10.3389/fmicb.2017.01463) · PMID 28848501 · PMC5552789 (open access)

44 compounds tested, 40 of them EPA priority pollutants. Two basidiomycete UPOs:
- **AaeUPO** (*Agrocybe*/*Cyclocybe aegerita*), long type / Family II, 44 kDa
- **MroUPO** (*Marasmius rotula*), short type / Family I, 29 kDa

Their Table 2, quoted verbatim. Conversion scale: `0` none, `t` traces, `+` <20 % substrate loss,
`++` <50 %, `+++` >50 %.

| Substrate | C atoms | Min/max width (Å) | AaeUPO | MroUPO | Water sol. (mg/L) | IP (eV) |
|---|---|---|---|---|---|---|
| Acenaphthylene | 12 | 4.0 / 4.6 | +++ | +++ | 16.1 | 8.22 |
| Acenaphthene | 12 | 4.0 / 4.6 | ++ | ++ | 3.93 | 7.86 |
| Fluorene | 13 | 3.0 / 7.4 | +++ * | n.t. | 1.98 | 7.89 |
| Phenanthrene | 14 | 4.0 / 7.6 | ++ * | n.t. | 1.29 | 7.91 |
| Anthracene | 14 | 3.0 / 7.9 | ++ * | n.t. | 0.073 | 7.43 |
| Pyrene | 16 | 5.2 / 7.6 | ++ * | n.t. | 0.135 | 7.43 |
| Benzo[a]anthracene | 18 | 4.8 / 8.1 | ++ | + | 0.014 | 7.44 |
| Benzo[a]pyrene | 20 | 4.8 / 8.1 | + | ++ | 0.0038 | 7.12 |
| Perylene | 20 | 4.6 / 6.6 | t | + | 0.0004 | 6.97 |
| Indeno[1,2,3-cd]pyrene | 22 | 5.8 / 8.4 | t | t | 0.062 | – |
| Benzo[b]fluoranthene | 20 | 4.8 / 9.6 | 0 | t | 0.0012 | 7.70 |
| Benzo[k]fluoranthene | 20 | 4.6 / 10.1 | t | t | 0.00055 | 7.48 |
| Dibenz[a,h]anthracene | 22 | 6.9 / 10.6 | **0** | **0** | 0.0005 | 7.38 |
| Benzo[g,h,i]perylene | 22 | 6.7 / 6.9 | **0** | **0** | 0.00026 | 7.16 |

`*` these four values are cited by Karich et al. from Aranda et al. 2010, not re-measured.
`n.t.` not tested.

Smaller arenes, from the same Table 1 (AaeUPO, cited from earlier papers): benzene `+`, naphthalene `++`,
phenol `+++`, anisole `+++`, toluene `+++`, ethylbenzene `+++`.

### What the authors conclude, in their words

- MroUPO "was capable of hydroxylating more bulky PAHs than" AaeUPO.
- This is attributed to "the wider heme channel of [MroUPO] (**11 Å**) compared to the relatively narrow
  channel of [AaeUPO] (**7 Å**)", which limits AaeUPO to PAHs no larger than about **6 Å in diameter**.
- Three stated limits on conversion overall: (i) steric hindrance / bulkiness, (ii) ring deactivation by
  electron-withdrawing groups, (iii) low water solubility.
- **UPOs favour low-molecular-mass PAHs.** Ligninolytic peroxidases (MnP, LiP) and laccase do the opposite —
  they oxidise the *large* PAHs faster and cannot touch phenanthrene or fluoranthene. Substrate ionization
  potential, which governs the laccase/peroxidase route, "seems to be of minor relevance" for UPOs.
- Products were mono- and dihydroxylated, not quinones, so these are true oxygen transfers rather than
  one-electron oxidation followed by water addition.

### The pattern that matters for a screen

Minimum width separates the data better than maximum width or carbon count:
- Every substrate with min width ≤ 5.8 Å is converted by at least one enzyme.
- Both substrates with min width ≥ 6.7 Å give zero for both enzymes.
- Benzo[b]fluoranthene has a max width of 9.6 Å and is still converted (traces, MroUPO) because its min
  width is 4.8 Å.

This supports the "slot, not tube" design point: a probe-sphere radius would misrank planar PAHs.

## 3. Supporting evidence for the channel-width relationship

Three independent lines, all consistent.

- **Olmedo et al. 2017**, Chemistry 23:16985. MroUPO shortens fatty acid chains one carbon at a time, which
  AaeUPO cannot do. Attributed to "a wider heme access channel, enabling fatty acid positioning with the
  carboxylic end near the heme cofactor". [DOI](https://doi.org/10.1002/chem.201704773) · PMID 29083064
- **Ullrich et al. 2018**, J Inorg Biochem 183:84. Models on the MroUPO crystal structure show "the bulky
  cortisone molecule suitably fits into the enzyme's access channel".
  [DOI](https://doi.org/10.1016/j.jinorgbio.2018.03.011) · PMID 29604496
- **Ramirez-Escudero et al. 2018**, ACS Chem Biol 13:3259. Crystal structure of a laboratory-evolved AaeUPO.
  Directed evolution produced "a broader heme access channel", and the authors name "dynamic trafficking
  through the heme channel as the main driving force for the exceptional substrate promiscuity".
  [DOI](https://doi.org/10.1021/acschembio.8b00500) · PMID 30376293

So channel width is not an inference from one paper. It is the accepted explanation, it has been changed
experimentally by directed evolution, and the change went in the predicted direction.

## 4. Single-residue evidence for a sequence-level screen

**Linde et al. 2020**, Appl Environ Microbiol 86:e02899-19. Two new ascomycete UPOs, CviUPO
(*Collariella virescens*) and DcaUPO (*Daldinia caldariorum*), expressed in *E. coli*.
[DOI](https://doi.org/10.1128/AEM.02899-19) · PMID 31980430 · PMC7082571

- The two differ 10-20 fold in catalytic efficiency on simple aromatics (veratryl alcohol, **naphthalene**,
  benzyl alcohol).
- Homology models showed "three conserved and two differing residues in the distal side of the heme", the
  difference being **two different positions of a phenylalanine**.
- Replacing **Phe88** with the homologous residue of the other enzyme gave an **F88L variant with 5- to
  21-fold higher efficiency** on those aromatics.

This is the strongest argument that a Tier 1 sequence-only readout of channel-lining residues can carry
real information: one phenylalanine position moved aromatic activity by up to 21-fold.

## 5. A second discriminating enzyme pair

**Aranda, Ullrich, Hofrichter 2010**, Biodegradation 21:267-81. AaeUPO ("AaP") vs *Coprinellus radians*
peroxygenase ("CrP") on 2- and 1-methylnaphthalene, dibenzofuran, fluorene, phenanthrene, anthracene,
pyrene. [DOI](https://doi.org/10.1007/s10532-009-9299-2) · PMID 19823936

- AaeUPO hydroxylated the aromatic ring of **all** substrates, at varying positions.
- CrP "showed a limited capacity for aromatic ring-hydroxylation and **did not hydroxylate phenanthrene**",
  preferring the non-aromatic C9 of fluorene and the side chain of methylnaphthalenes.

A second enzyme with a published inability to attack a specific 3-ring PAH. Useful as a negative in the
calibration set. Not read in full; the above is the abstract.

## 6. Within-organism contrast

**Song et al. 2024**, J Hazard Mater 471:134437. *Moesziomyces aphidis* XM01 (Ustilaginomycetes,
Basidiomycota) degrades crude oil using two extracellular UPOs with **opposite preferences**: MaUPO.1 prefers
PAHs, MaUPO.2 prefers n-alkanes. [DOI](https://doi.org/10.1016/j.jhazmat.2024.134437) · PMID 38691934

Two paralogs in one genome split the substrate space. Directly relevant to the lineage-expansion question in
goal 1: expansions may partition substrate classes rather than duplicate one function.

## 7. Outgroup / family boundary

**Lopez-Tavera et al. 2026**, Biochemistry 65:2283. Bacterial heme-thiolate peroxygenases (BUPOs) found by
structure database mining. Clear structural homology to **Family I (short) fungal UPOs** at <20 % sequence
identity. Active on **naphthalene**, indole, 3-phenyl-1-propanol, 16-hydroxypalmitic acid.
[DOI](https://doi.org/10.1021/acs.biochem.6c00230) · PMID 42374716 · PMC13394420

Two consequences. Structure-based search finds real family members that sequence identity misses, which
supports the Tier 2 structural approach. And the family is not fungal-only, so any "first non-Dikarya UPO"
style claim needs care.

## 8. Limits of this calibration set — read before designing the screen

1. **n is small.** Only AaeUPO and MroUPO have the full PAH panel measured side by side. Two calibrated
   points cannot fit a curve. We can anchor the direction and a rough cutoff, nothing more.
2. **The readout is ordinal, not kinetic.** `0 / t / + / ++ / +++` is semi-quantitative substrate depletion
   at one condition. There are no kcat or KM values for the PAH panel. So the screen's target variable is an
   ordered category, and a regression against it would be overfitting a 5-level scale.
3. **Size and solubility are confounded, and the authors say so.** Larger PAHs are simultaneously bulkier and
   far less water soluble — across the panel, solubility falls from 16.1 to 0.00026 mg/L. Conversion falls
   with both. This dataset cannot separate a steric effect from a bioavailability effect. Any claim that a
   predicted wide channel will convert a large PAH inherits this confound.
4. **The biggest PAHs may be the wrong target for this enzyme family.** Both UPOs gave zero on
   dibenz[a,h]anthracene and benzo[g,h,i]perylene, and the authors state that UPOs favour low-mass PAHs while
   laccase and ligninolytic peroxidases favour high-mass ones. A screen for 5-ring activity may be looking
   for something the family does not do. Finding a copy that does would be a genuine result — but the prior
   is against it, and the project should say so rather than assume the band is populated.
5. **Channel widths of 7 Å and 11 Å are quoted by Karich et al. from secondary sources**
   (Poraj-Kobielska thesis; Piontek et al.). We have not read the primary measurement, nor confirmed how the
   width was defined or at what point along the channel. Before these numbers anchor a screen, read the
   primary source and reproduce the measurement on 2YOR and 5FUJ with whatever tool the screen will use.
6. **The Phe88 result is from one enzyme pair** and was measured on simple aromatics (naphthalene, benzyl
   alcohol, veratryl alcohol), not on 3- to 5-ring PAHs. Do not extrapolate it to the large-PAH end.

## 9. What this means for the screen design

Supported by the data:
- Rank by **minimum channel cross-section**, and compare against substrate **minimum width**. Both the
  published table and the physical argument agree this is the discriminating dimension.
- Use AaeUPO (narrow, ~7 Å) and MroUPO (wide, ~11 Å) as the two anchor points, with their full panel as the
  validation target. A screen that cannot reproduce "MroUPO accepts bulkier PAHs than AaeUPO" is not working.
- Include a fatty acid in the ladder. MroUPO's wide channel gives it chain-shortening activity AaeUPO lacks,
  so the fatty-acid and large-PAH predictions are coupled through the same geometry.
- Read out channel-lining phenylalanine positions in Tier 1. Phe88 alone moved activity 5-21x.

Not supported by the data:
- Any continuous or quantitative prediction of activity.
- Confident calls in the 4-5 ring band. Report those as extrapolation beyond the calibrated range, and state
  that both calibration enzymes failed on the two largest substrates tested.

## 10. Still to check

- Primary source for the 7 Å / 11 Å channel widths (Piontek et al.; Poraj-Kobielska thesis). Not read.
- Aranda et al. 2010 full text, for the per-substrate detail behind the AaeUPO vs CrP contrast. Abstract only.
- Whether any UPO has published activity on a 5-ring PAH above trace level. Not found in this search.
- Whether MthUPO, HspUPO, DcaUPO, CviUPO or the newer PDB entries (9TM7, 9WA6) have PAH panel data. Not
  searched.
- Published Pichia / *Komagataella* expression outcomes per enzyme, for the expressibility filter. Not
  searched.

---

# ADDENDUM (same day): a 19-enzyme calibration set, and a change of target variable

## 11. Why channel width fails in the 1-3 ring range

Checked against the Karich Table 2 numbers in section 2. Minimum widths of the 1-3 ring substrates:
fluorene 3.0 Å, anthracene 3.0 Å, phenanthrene 4.0 Å, acenaphthylene and acenaphthene 4.0 Å.
The discriminating cutoff is about 6 Å. **Every 1-3 ring substrate clears both the narrow (7 Å) and the wide
(11 Å) channel**, and AaeUPO converts all of them at `++` or `+++`.

Conclusion: channel width separates 4-5 ring PAHs from everything else. Inside the 1-3 ring band it is flat
and carries no signal. A model keyed on width would be fitting noise in this range.

## 12. The better calibration set — Biotechnol Biofuels 2025

"Structure-function relationships in unspecific peroxygenases revealed by a comparative study of their action
on the phenolic lignin monomer 4-propylguaiacol." Biotechnol Biofuels Bioprod.
[DOI](https://doi.org/10.1186/s13068-025-02675-w) · PMID 40722037 · PMC12306115 (open access)

This is close to the study we were proposing to do computationally, already done experimentally.

- **19 UPOs**, chosen across **5 monophyletic clades (I-V)**, biased toward enzymes with solved crystal
  structures or demonstrated recombinant expression, plus one uncharacterized UPO from a distinct clade.
- **All expressed in one host and assayed under identical conditions.** This is what the Karich comparison
  lacked. Cross-study unit mismatch is not an issue within this dataset.
- Substrates: ABTS and 2,6-DMP (peroxidase), NBD, **naphthalene**, 4-nitroanisole (peroxygenase), and
  4-propylguaiacol (4-PG) in depth.
- Naphthalene assay, fixed: 0.1 µM UPO, 1 mM naphthalene in 10 % acetone, 1 mM H2O2, 100 mM citrate pH 4.5,
  A324, initial rates from the linear region.

### Why this fits the recentred 1-3 ring scope

4-PG is a single aromatic ring. Naphthalene is two rings. Both sit inside the band we chose.

### Measured outcomes we can use

- **Naphthalene negatives:** no activity detected for four of the 19 (AbrUPO-II, CviUPO, PanUPO, TteUPO).
  A classifier needs negatives; these are negatives measured under the same assay as the positives.
- **Clade does not determine function.** "Some UPOs from the same phylogenetic clade exhibited distinct
  regioselectivities." But two UPOs at 96 % identity gave the same regioselectivity. So the signal lives
  below clade level, at individual residues.
- **Peroxygenase-to-peroxidase ratio varies widely** across the 19 and is the headline functional axis.

### The structural determinants they report

- **Clade V** (Ascomycota): narrow substrate channel **plus an "aliphatic pocket"** near the heme, conserved
  across the clade, residues **I73, L77, L81, A183, I187, L228**. These show marked peroxygenase activity.
- **Clades I-IV:** funnel-shaped access channel with a narrow base, **no aliphatic pocket**, mainly
  peroxidase activity toward 4-PG.
- **AbrUPO-II** breaks its clade: constrained entrance and narrow channel give it high peroxygenase activity
  despite sitting in Clade I. Against PaDa-I at 73 % identity, the difference is a more restricted entrance,
  largely **Phe302 (Pro277 in PaDa-I)**, plus bulky residues narrowing the channel.
- Their mechanistic reading: a constrained channel plus an aliphatic pocket holds the substrate near the iron
  long enough for hydrogen abstraction and oxygen rebound (peroxygenase). A broad funnel lets substrate in
  and product out quickly, which favours one-electron oxidation (peroxidase).
- Geometric criterion: about **4 Å** between the reactive hydrogen and the heme iron is needed for efficient
  rebound to Compound II.

### Methods they used — the same stack we specified

- **CAVER (Caver Analyst 2.0)** for access channel analysis.
- **AlphaFold3** to model the one UPO without a crystal structure.
- **Docking**, with solutions filtered by the oxidized carbon sitting **3.5-4 Å from the heme** and binding
  energy at least 20.72 kJ/mol.

That the published study converged on CAVER plus AF3 plus a heme-distance docking filter is evidence the
Tier 2 and Tier 3 design is workable. It is not evidence that our implementation will reproduce their result.

## 13. Revised target variable for the model

Not "predicted ring-size band". In the 1-3 ring range the useful prediction is:

1. **Peroxygenase vs peroxidase propensity.** This is the industrially decisive axis and the structurally
   predictable one. It is also exactly the overoxidation failure mode flagged earlier: peroxidase activity on
   a phenolic product is what drives radical coupling and polymerization.
2. **Naphthalene activity, present or absent.** 19 enzymes, one assay, four measured negatives.

Candidate predictive features, all from the paper above:
- presence of the aliphatic pocket (the six named residues, mapped by structural alignment)
- channel shape: constrained versus funnel with a narrow base
- entrance residues, with Phe302/Pro277 as the worked example
- channel-lining phenylalanines, with Phe88 from Linde et al. 2020 as a second worked example

## 14. Honest limits of the revised plan

1. **n = 19 is enough to rank and to sanity-check, not to fit a many-parameter model.** With four naphthalene
   negatives, any classifier has four negative examples. Keep the feature count very small or it will overfit.
2. **The detailed substrate is 4-PG, a lignin monomer, not a PAH.** Naphthalene data from this study appears
   to be presence/absence plus relative percentages, not kinetic constants. Check the supplementary tables
   before assuming otherwise.
3. **Direct conflict to resolve: CviUPO.** No naphthalene activity in this study; Linde et al. 2020
   ([DOI](https://doi.org/10.1128/AEM.02899-19)) report a catalytic efficiency for naphthalene oxidation by
   CviUPO. Different expression hosts and assay conditions are the likely cause. Resolve before CviUPO is
   used as either a positive or a negative.
4. **The aliphatic pocket is defined on Clade V Ascomycota.** Whether those six positions are recognisable in
   early-diverging or lichen sequences at 30 % identity is untested, and is the main risk to transferring
   this to the BFD set.
5. Supplementary tables S1 (the 19 enzymes), S2 (docking), S4 (colorimetric activities) have **not** been
   retrieved. All per-enzyme numbers above are from the main text or the abstract.

## 15. Next actions for the model

1. Retrieve supplementary tables S1, S2, S4 from PMC12306115 and build the enzyme x assay matrix.
2. Add the Karich 1-3 ring rows and the Aranda AaeUPO/CrP contrast as independent held-out cases.
3. Resolve the CviUPO conflict.
4. Map the aliphatic pocket positions and the channel-entrance positions onto a structural alignment, and
   test whether they are callable at BFD-scale sequence identity before any full run.
