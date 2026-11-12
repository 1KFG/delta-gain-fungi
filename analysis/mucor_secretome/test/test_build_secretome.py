"""Offline tests for build_secretome.py and extract_proteins.py.
Run: python3 -m unittest discover -s mucor_secretome/test"""
import gzip
import io
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import build_secretome as bs  # noqa: E402
import extract_proteins as ep  # noqa: E402
import annotate_to_tables as at  # noqa: E402
import gff3_to_cds as gc  # noqa: E402
import unite_tax as ut  # noqa: E402


class Classify(unittest.TestCase):
    def test_topology_parsing(self):
        self.assertEqual(bs.parse_topology("i7-29o44-66i"), [(7, 29), (44, 66)])
        self.assertEqual(bs.parse_topology(""), [])

    def test_signal_peptide_helix_is_not_a_membrane_helix(self):
        # the only helix sits inside the signal peptide (cleaved at 22)
        self.assertEqual(bs.membrane_helices("n5-20o", 1, 22), 0)
        # a helix well past the cleavage site is a real TM
        self.assertEqual(bs.membrane_helices("i60-82o", 1, 22), 1)

    def test_membrane_helices_without_topology_falls_back_to_count(self):
        self.assertEqual(bs.membrane_helices("", 3, 22), 2)
        self.assertEqual(bs.membrane_helices(None, 1, 0), 1)

    def test_classify(self):
        self.assertEqual(bs.classify(0.99, "n5-20o", 1, 22, False, 0.0), "secreted")
        self.assertEqual(bs.classify(0.6, "", 0, 20, False, 0.0), "weak_signal")
        self.assertEqual(bs.classify(0.99, "i60-82o", 1, 22, False, 0.0), "membrane")
        self.assertEqual(bs.classify(0.99, "", 0, 22, True, 0.0), "gpi_anchored")
        self.assertEqual(bs.classify(0.99, "", 0, 22, False, 0.9), "mitochondrial")
        # a membrane helix outranks a GPI anchor
        self.assertEqual(bs.classify(0.99, "i60-82o", 1, 22, True, 0.0), "membrane")


class Tiers(unittest.TestCase):
    def test_chitin_beats_plant(self):
        t, why = bs.tier_of(["CE4(145-269)"], ["chitin"], [], 300, 5, 0.02)
        self.assertEqual(t, "T1")
        self.assertIn("CE4", why)

    def test_lysm_pfam_is_t1(self):
        self.assertEqual(bs.tier_of([], [], ["PF01476"], 200, 0, 0.0)[0], "T1")

    def test_plant_substrate_is_t2(self):
        t, why = bs.tier_of(["GH28(20-300)"], ["pectin, polygalacturonan"], [], 350, 2, 0.01)
        self.assertEqual(t, "T2")
        self.assertIn("pectin", why)

    def test_upo_pfam_is_t2(self):
        self.assertEqual(bs.tier_of([], [], ["PF01328"], 280, 6, 0.03)[0], "T2")

    def test_effector_like_needs_no_domain(self):
        self.assertEqual(bs.tier_of([], [], [], 180, 8, 0.045)[0], "T3")
        # same protein with a Pfam domain is not "unknown"
        self.assertEqual(bs.tier_of([], [], ["PF99999"], 180, 8, 0.045)[0], "T4")
        # too long
        self.assertEqual(bs.tier_of([], [], [], 400, 8, 0.02)[0], "T4")
        # too few cysteines
        self.assertEqual(bs.tier_of([], [], [], 180, 3, 0.05)[0], "T4")

    def test_cazy_root(self):
        self.assertEqual(bs.cazy_root("CE4(145-269)"), "CE4")
        self.assertEqual(bs.cazy_root("GH18_3"), "GH18")
        self.assertEqual(bs.cazy_root(""), "")

    def test_cys_stats_ignores_trailing_stop(self):
        n, frac, length = bs.cys_stats("MCCAA*")
        self.assertEqual((n, length), (2, 5))
        self.assertAlmostEqual(frac, 0.4)


class Extract(unittest.TestCase):
    def test_read_fasta_filters_and_keeps_last_record(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "p.fa.gz")
            with gzip.open(p, "wt") as fh:
                fh.write(">A_1-T1 gene1\nMKV\nAA\n>A_2-T1 gene2\nMMM\n")
            self.assertEqual(dict(ep.read_fasta(p)), {"A_1-T1": "MKVAA", "A_2-T1": "MMM"})
            self.assertEqual(dict(ep.read_fasta(p, {"A_2-T1"})), {"A_2-T1": "MMM"})

    def test_read_index(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "genomes.tsv")
            with open(p, "w") as fh:
                fh.write("LT1\tout1\t/x/p.fa.gz\t/x/c.fa.gz\n\nLT2\tout2\t/y/p.fa.gz\n")
            self.assertEqual(ep.read_index(p), {"LT1": "/x/p.fa.gz", "LT2": "/y/p.fa.gz"})



class Gff3(unittest.TestCase):
    GFF = (
        "##gff-version 3\n"
        "s1\tfunannotate\tmRNA\t100\t900\t.\t-\t.\tID=T1;Parent=G1;\n"
        "s1\tfunannotate\tCDS\t800\t900\t.\t-\t0\tID=T1.cds;Parent=T1;\n"
        "s1\tfunannotate\tCDS\t100\t300\t.\t-\t2\tID=T1.cds;Parent=T1;\n"
        "s1\tfunannotate\texon\t100\t300\t.\t-\t.\tID=T1.exon;Parent=T1;\n"
        "s2\tfunannotate\tCDS\t10\t40\t.\t+\t0\tID=T2.cds;Parent=T2;\n"
    )

    def test_parse(self):
        got = gc.parse_gff3(io.StringIO(self.GFF))
        self.assertEqual(got["T1"], ("s1", -1, [(100, 300), (800, 900)]))
        self.assertEqual(got["T2"], ("s2", 1, [(10, 40)]))
        self.assertNotIn("G1", got)   # only CDS features are used

    def test_attrs(self):
        self.assertEqual(gc.attrs("ID=a;Parent=b;"), {"ID": "a", "Parent": "b"})

class Annot(unittest.TestCase):
    SP6 = (
        "# SignalP-6.0\tOrganism: Eukarya\n"
        "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
        "P1 desc\tSP\t0.0021\t0.9979\tCS pos: 22-23. Pr: 0.98\n"
        "P2\tOTHER\t0.9900\t0.0100\t\n"
    )
    TMHMM = ("P1\tlen=255\tExpAA=20.11\tFirst60=19.90\tPredHel=1\tTopology=n5-20o\n"
             "# a comment\n")
    DOMTBL = (
        "# target name        accession   tlen query name   accession  qlen ...\n"
        "P1 - 255 LysM PF01476.25 44 1e-30 100.0 0.1 1 1 1e-12 2.5e-9 40.0 0.1 "
        "1 44 31 74 30 75 0.95 desc here\n"
    )

    def test_signalp6_probability_is_one_minus_other(self):
        got = at.parse_signalp6(io.StringIO(self.SP6))
        self.assertEqual(got[0], ("P1", 0.9979, 22, 23))
        self.assertEqual(got[1], ("P2", 0.01, 0, 0))

    def test_tmhmm_short(self):
        got = at.parse_tmhmm_short(io.StringIO(self.TMHMM))
        self.assertEqual(got, [("P1", 255, 20.11, 19.9, 1, "n5-20o")])

    def test_domtbl_target_is_protein_query_is_hmm(self):
        got = at.parse_domtbl(io.StringIO(self.DOMTBL))
        self.assertEqual(got[0][:3], ("P1", "LysM", "PF01476.25"))
        self.assertAlmostEqual(got[0][3], 2.5e-9)
        self.assertEqual(got[0][4:], (31, 74))

    def test_fam_from_hmm(self):
        self.assertEqual(at.fam_from_hmm("GH18.hmm"), "GH18")
        self.assertEqual(at.fam_from_hmm("CBM50"), "CBM50")


class BfdParity(unittest.TestCase):
    """Each parser against the merge_*.py it mirrors in Fungi_BFD/nextflow/bin."""

    def test_signalp_gff3_is_sp_positive_only_and_scores_from_column_6(self):
        gff = ("##gff-version 3\n"
               "# comment\n"
               "P1 some description\tSignalP-6.0\tsignal_peptide\t1\t22\t0.9979\t.\t.\t.\n"
               "P3\tSignalP-6.0\tlipoprotein_signal_peptide\t1\t19\t0.8810\t.\t.\t.\n")
        got = at.parse_signalp_gff3(io.StringIO(gff))
        # merge_signalp.py takes the id up to the first space and columns 4, 5, 6
        self.assertEqual(got, [("P1", "0.9979", "1", "22"), ("P3", "0.8810", "1", "19")])

    def test_predgpi_keeps_only_anchor_rows(self):
        gff = ("##gff-version 3\n"
               "P1\tPredGPI\tChain\t1\t255\t1.0\t.\t.\tevidence=ECO:0000256\n"
               "P1\tPredGPI\tGPI-anchor\t230\t230\t0.99\t.\t.\t.\n"
               "P2\tPredGPI\tChain\t1\t100\t1.0\t.\t.\t.\n")
        self.assertEqual(at.parse_predgpi_gff3(io.StringIO(gff)),
                         [("P1", "GPI-anchor", "230", "230")])

    def test_targetp_drops_notp_and_scores_the_predicted_class(self):
        tp = ("# TargetP-2.0\n"
              "P1\tSP\t0.0012\t0.9980\t0.0008\tCS pos: 22-23. ANA-AT. Pr: 0.76\n"
              "P2\tmTP\t0.0100\t0.0300\t0.9600\tCS pos: 31-32. VRA-SS. Pr: 0.51\n"
              "P3\tnoTP\t0.9900\t0.0050\t0.0050\t\n")
        got = at.parse_targetp(io.StringIO(tp))
        self.assertEqual([r[0] for r in got], ["P1", "P2"])
        # the SP protein's probability belongs to SP, not to mTP - build_secretome.py must look
        # at `prediction` before reading it as a mitochondrial score
        self.assertEqual(got[0][1:3], ("SP", "0.9980"))
        self.assertEqual(got[1][1:], ("mTP", "0.9600", "31", "32"))

    def test_cazy_overview_carries_substrates_one_row_per_family(self):
        ov = ("Gene_ID\tEC\tcazyme_fam\tsub_fam\tdiamond_fam\tSubstrate\t#ofTools\n"
              "P1\t3.2.1.14\tGH18+CBM18\tGH18_e1+-\tGH18\tchitin\t3\n"
              "P2\t-\tGH28\t-\tGH28\tpectin, polygalacturonan\t2\n")
        got = at.parse_cazy_overview(io.StringIO(ov))
        self.assertEqual([(r[0], r[2], r[4]) for r in got],
                         [("P1", "GH18", "chitin"), ("P1", "CBM18", "chitin"),
                          ("P2", "GH28", "pectin, polygalacturonan")])
        self.assertEqual(got[0][1], "3.2.1.14")
        self.assertEqual(got[1][3], "")          # '-' sub_fam becomes empty
        self.assertEqual(got[2][1], "")          # '-' EC becomes empty

    def test_cazy_overview_rows_feed_the_plant_tier(self):
        rows = at.parse_cazy_overview(io.StringIO(
            "Gene_ID\tEC\tcazyme_fam\tsub_fam\tdiamond_fam\tSubstrate\t#ofTools\n"
            "P2\t-\tGH28\t-\tGH28\tpectin, polygalacturonan\t2\n"))
        fams = [f"{r[2]}" for r in rows]
        subs = [r[4] for r in rows]
        self.assertEqual(bs.tier_of(fams, subs, [], 350, 2, 0.01)[0], "T2")


class Unite(unittest.TestCase):
    TAX = (
        "Feature ID\tTaxon\n"
        "SH1.10FU_X1_refs\tk__Fungi;p__Ascomycota;c__D;o__A;f__A;g__Abrothallus;"
        "s__Abrothallus_subhalei;sh__SH1\n"
        "SH2.10FU_X2_refs\tk__Fungi;p__Mucoromycota;c__Mucoromycetes;o__Mucorales;"
        "f__Mucoraceae;g__Mucor;s__Mucor_hiemalis;sh__SH2\n"
        "SH3.10FU_X3_refs\tk__Fungi;p__Mucoromycota;c__Mucoromycetes;o__Mucorales;"
        "f__Mucoraceae;g__Mucor;s__Mucor_circinelloides;sh__SH3\n"
    )
    FASTA = ">SH1.10FU_X1_refs\nAAAA\n>SH2.10FU_X2_refs desc\nCCCC\n>SH3.10FU_X3_refs\nGGGG\n"

    def tax(self):
        return ut.parse_taxonomy(io.StringIO(self.TAX))

    def test_header_row_is_not_a_record(self):
        self.assertNotIn("Feature ID", self.tax())
        self.assertEqual(len(self.tax()), 3)

    def test_rank_of(self):
        lin = self.tax()["SH2.10FU_X2_refs"]
        self.assertEqual(ut.rank_of(lin, "p"), "Mucoromycota")
        self.assertEqual(ut.rank_of(lin, "s"), "Mucor_hiemalis")
        self.assertEqual(ut.rank_of(lin, "x"), "")

    def test_phylum_filter_uses_the_taxonomy_not_the_header(self):
        # the headers carry no lineage at all, which is why a grep on the FASTA finds nothing
        self.assertNotIn("Mucoromycota", self.FASTA)
        ids = ut.ids_in_phylum(self.tax(), "Mucoromycota")
        out = io.StringIO()
        self.assertEqual(ut.filter_fasta(io.StringIO(self.FASTA), ids, out), 2)
        self.assertEqual(out.getvalue().count(">"), 2)
        self.assertIn("CCCC", out.getvalue())
        self.assertNotIn("AAAA", out.getvalue())

    def test_filter_fasta_keeps_only_the_accession_part_of_the_header(self):
        out = io.StringIO()
        ut.filter_fasta(io.StringIO(self.FASTA), {"SH2.10FU_X2_refs"}, out)
        self.assertEqual(out.getvalue(), ">SH2.10FU_X2_refs desc\nCCCC\n")

    def test_rank_hits_keeps_the_best_hit_per_species(self):
        def row(acc, pid, ln, bits, qlen=500):
            return dict(zip(ut.BLAST_COLS,
                            [acc, "ctg1", pid, ln, "0", "0", "1", ln, "1", ln, "0.0", bits, qlen]))
        ranked = ut.rank_hits([row("SH2.10FU_X2_refs", "99.0", "450", "800"),
                               row("SH2.10FU_X2_refs", "95.0", "400", "700"),
                               row("SH3.10FU_X3_refs", "88.0", "300", "500")], self.tax())
        self.assertEqual([r[0] for r in ranked], ["Mucor_hiemalis", "Mucor_circinelloides"])
        self.assertEqual(ranked[0][2:5], (99.0, 450, 0.9))   # the 800-bit hit, not the 700

    def test_rank_hits_falls_back_to_the_accession_when_unknown(self):
        row = dict(zip(ut.BLAST_COLS,
                       ["NOPE", "c", "90", "100", "0", "0", "1", "100", "1", "100", "0", "90",
                        "100"]))
        self.assertEqual(ut.rank_hits([row], self.tax())[0][0], "NOPE")


class Gpi(unittest.TestCase):
    def test_only_anchor_rows_count(self):
        rows = [{"protein_id": "A", "feature": "Chain"},
                {"protein_id": "B", "feature": "GPI-anchor"},
                {"protein_id": "C", "feature": "Chain"},
                {"protein_id": "C", "feature": "GPI-anchor"}]
        self.assertEqual(bs.gpi_set(rows), {"B", "C"})

    def test_missing_feature_column_is_treated_as_anchor(self):
        self.assertEqual(bs.gpi_set([{"protein_id": "A"}]), {"A"})


if __name__ == "__main__":
    unittest.main()


class UniteRanks(unittest.TestCase):
    def test_subset_at_any_rank(self):
        tax = ut.parse_taxonomy(io.StringIO(Unite.TAX))
        self.assertEqual(ut.ids_at_rank(tax, "g", "Mucor"),
                         {"SH2.10FU_X2_refs", "SH3.10FU_X3_refs"})
        self.assertEqual(ut.ids_at_rank(tax, "p", "Ascomycota"), {"SH1.10FU_X1_refs"})
        # the phylum helper is the same call, so both subsets stay in step
        self.assertEqual(ut.ids_in_phylum(tax, "Mucoromycota"),
                         ut.ids_at_rank(tax, "p", "Mucoromycota"))


import call_effectors as ce  # noqa: E402
import fold_report as fr  # noqa: E402


class Effectors(unittest.TestCase):
    def test_fold_reason_picks_unknown_small_and_chitin(self):
        # no Pfam, short -> folded even though it is cysteine-poor (T3 would have dropped it)
        self.assertEqual(ce.fold_reason("secreted", "T4", [], 150, 300), "unknown_small")
        # a chitin-active enzyme is folded for its binding site, Pfam or not
        self.assertEqual(ce.fold_reason("secreted", "T1", ["PF01522"], 400, 300), "chitin_active")
        # a domain-bearing protein that is not chitin-active is not folded
        self.assertEqual(ce.fold_reason("secreted", "T2", ["PF01328"], 280, 300), "")
        # too long to be an unknown effector
        self.assertEqual(ce.fold_reason("secreted", "T4", [], 600, 300), "")
        # not secreted at all
        self.assertEqual(ce.fold_reason("membrane", "T4", [], 150, 300), "")
        self.assertEqual(ce.fold_reason("gpi_anchored", "T4", [], 150, 300), "unknown_small")

    def test_mito_prob_only_counts_mtp(self):
        self.assertEqual(ce.mito_prob({"prediction": "mTP", "probability": "0.96"}), 0.96)
        # an SP-predicted protein's probability belongs to SP, not to mTP
        self.assertEqual(ce.mito_prob({"prediction": "SP", "probability": "0.99"}), 0.0)
        self.assertEqual(ce.mito_prob({}), 0.0)


class FoldCalls(unittest.TestCase):
    HITS = (
        "P1\t1abc_A\t1abc_A LysM domain\t0.11\t60\t0.9\t0.8\t1e-3\t55\t0.72\t0.95\n"
        "P1\t2xyz_B\t2xyz_B something else\t0.09\t55\t0.8\t0.7\t1e-1\t40\t0.41\t0.60\n"
        "P2\t3def_A\t3def_A unrelated\t0.08\t40\t0.5\t0.4\t2.0\t20\t0.31\t0.20\n"
    )

    def test_only_hits_at_or_above_tm_threshold_count(self):
        best = fr.best_per_query(fr.read_hits_from_string(self.HITS), 0.5)
        self.assertEqual(set(best), {"P1"})           # P2's best is TM 0.31
        self.assertEqual(best["P1"]["target"], "1abc_A")

    def test_low_plddt_model_is_held_back(self):
        best = {"pdb": fr.best_per_query(fr.read_hits_from_string(self.HITS), 0.5)}
        rows = fr.call_folds(best, {"P1": {"plddt": "55.0", "ptm": "0.3", "length": "120"}})
        self.assertEqual(rows["P1"]["fold_call"], "low_confidence")
        rows = fr.call_folds(best, {"P1": {"plddt": "88.0", "ptm": "0.8", "length": "120"}})
        self.assertEqual(rows["P1"]["fold_call"], "assigned")
        self.assertEqual(rows["P1"]["best_tm"], 0.72)

    def test_folded_but_unmatched_is_a_result_not_a_failure(self):
        rows = fr.call_folds({"pdb": {}}, {"P9": {"plddt": "90.0", "ptm": "0.8", "length": "99"}})
        self.assertEqual(rows["P9"]["fold_call"], "no_match")


import plant_panel as pp  # noqa: E402


class PlantPanel(unittest.TestCase):
    def test_class_assignment_prefers_the_specific_class(self):
        # a LysM receptor kinase carries both LysM and a kinase domain; LysM wins
        self.assertEqual(pp.classify_protein(["PF01476", "PF00069"]), "lysm_receptor")
        self.assertEqual(pp.classify_protein(["PF00704"]), "chitinase")
        self.assertEqual(pp.classify_protein(["PF00112"]), "cys_protease")
        self.assertEqual(pp.classify_protein(["PF99999"]), "")

    def test_receptor_gives_an_ectodomain_not_the_whole_chain(self):
        # 600 aa receptor, SP cleaved at 25, single TM helix at 300-322
        start, end, kind = pp.construct_range(600, 25, "o26-299i300-322o", 1)
        self.assertEqual((start, end, kind), (26, 299, "ectodomain"))
        # the kinase domain (323-600) is excluded
        self.assertLess(end, 323)

    def test_secreted_protein_with_no_helix_is_the_mature_chain(self):
        self.assertEqual(pp.construct_range(300, 22, "", 0), (23, 300, "mature"))

    def test_helix_inside_the_signal_peptide_is_not_a_tm(self):
        # the only helix sits within the cleavage site + slack, so the chain stays "mature"
        self.assertEqual(pp.construct_range(300, 22, "n5-20o", 1)[2], "mature")

    def test_no_signal_peptide_is_not_a_surface_target(self):
        self.assertEqual(pp.construct_range(300, 0, "i10-32o", 1)[2], "none")

    def test_ectodomain_too_short_to_model_is_dropped(self):
        # helix starts at 40, leaving 15 residues outside the membrane
        self.assertEqual(pp.construct_range(400, 24, "o25-39i40-62o", 1)[2], "none")
