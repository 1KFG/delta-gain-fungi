# Unverified items: pathway `carotenoid`

- Source: `pathways/carotenoid.yaml`
- Pathway-level `verified`: **false**
- `validate()` errors: 0
- Refs: 23 total, 22 verified, 1 unverified.

## Refs with verified: false

| step | accession | species | pmid |
|---|---|---|---|
| astaxanthin_synthase | - | Xanthophyllomyces dendrorhous | 16416328 |

## Refs missing an accession or a PMID

- astaxanthin_synthase: accession=-, pmid=16416328

## Steps with no verified ref

None.

---

## Items for a person to check (hand-written, 2026-11-20)

The section above is the output of `bin/list_unverified_refs.py pathways/carotenoid.yaml`.
The items below are not reference checks. They are decisions in the draft that no source settled.
Labels: MEASURED = read from a database record on 2026-11-20; PROPOSED = builder's choice; UNVERIFIED = not checked.

1. **Step list (user review gate, spec 5.8).** Review the 8 steps: `ggpp_synthase`, `phytoene_synthase`, `lycopene_cyclase`, `phytoene_desaturase`, `carx`, `cart`, `ylo1`, `astaxanthin_synthase`. `ylo1` was added (verified refs). `cao2` was not added as a step: cao-2 is a ref of `cart` (PROPOSED; PMID 17610084 shows CAO-2 does the torulene cleavage).
2. **Steps not included.** Mucorales beta-carotene cleavage for trisporic acid (carS and others; Fable review item 6), crtR (P450 reductase partner of astaxanthin synthase; PMID 18837978), Rhodotorula torulene/torularhodin steps, and carO (opsin, non-catalytic). None was curated. Decide if any belongs in v1.
3. **Lycopene cyclase domain model.** MEASURED from UniProt cross-references: NCBIfam TIGR03462 (CarR_dom_SF, InterPro IPR017825) is on all 5 fused refs; Pfam PF18916 is cross-referenced only on P37295 (al-2). Whether PF18916 detects the other fungal cyclase domains is UNVERIFIED. The draft lists `[TIGR03462, PF18916]`. Test both models on the 5 refs before the prefilter is built.
4. **copy_expect, optional, role.** All PROPOSED. `ggpp_synthase` is `multi` because PF00348 also covers non-carotenoid GGPP and FPP synthases; the carotenoid-specific clade is not defined. `carx`, `cart`, `ylo1`, `astaxanthin_synthase` are `optional: true` (lineage-specific or off the main route).
5. **Domain families are broad.** PF00348, PF00494, PF01593, PF03055, PF00171 and PF00067 are much broader than the steps. Decoy clades (squalene synthase for PF00494; non-carotenoid CCOs for PF03055; general P450s) must be in each reference tree. No decoys are curated yet.
6. **Fused proteins as full-length refs.** The 5 fused accessions are full-length in both `phytoene_synthase.faa` and `lycopene_cyclase.faa`. UniProt region coordinates (MEASURED) put the cyclase in roughly residues 1-238..1-268 and the synthase after it. Splitting into per-domain FASTA is left for the HMM/tree task.
7. **Reference coverage by lineage.** Refs come from Ascomycota (Neurospora, Fusarium), Mucoromycota (Phycomyces, Mucor) and Basidiomycota (Xanthophyllomyces only). No Rhodotorula/Pucciniomycotina, no chytrid and no Dothideomycete refs. Whether this spread is enough for tree placement is UNVERIFIED.
8. **Astaxanthin synthase accession.** B5LEC8 (TrEMBL) is cited by PMID 18769767 (cloning of crtS and mutants). The functional paper (PMID 16416328, Ojima 2006) has no linked accession. A person may wish to find the Ojima sequence and confirm it matches B5LEC8.
9. **Strain mixing in Fusarium refs.** Q8X0Z1, Q8X0Z0, A1KQY3 are "Fusarium fujikuroi" (strain not given); S0ENQ2 and F6IBC7 are strain IMI 58289. PMID 12172798 states carB codes for 541 aa; UniProt Q8X0Z0 is 570 aa at sequence version 2. The difference was not investigated.
