"""Tests for refs/carotenoid/tree_members/: named tree members and tree outgroups.

These sequences are tree members only. They are not references, not decoys, and no step call,
absence check or decoy-calibrated cutoff may read them.
"""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REFS = ROOT / "refs" / "carotenoid"
TM = REFS / "tree_members"
PROV = TM / "provenance.tsv"
COLS = ["step_id", "accession", "species", "class", "phylum", "gene_name", "pmid", "role",
        "protein_existence", "annotation_score", "evidence_codes", "ec", "region", "reason",
        "source_database", "retrieved_on", "verified"]
ROLES = {"named_member", "tree_outgroup"}
STEPS = {"phytoene_synthase", "lycopene_cyclase"}
MEMBER_FA = {s: TM / f"{s}.faa" for s in STEPS}
OUTGROUP_FA = {"phytoene_synthase": TM / "outgroups_phytoene_synthase.faa"}
# Second set (phytoene_desaturase, cart): one-domain proteins, full length, no cut.
# `cart` has no named member, so it has no member file.
NEW_STEPS = {"phytoene_desaturase", "cart"}
MEMBER_FA["phytoene_desaturase"] = TM / "phytoene_desaturase.faa"
OUTGROUP_FA["phytoene_desaturase"] = TM / "outgroups_phytoene_desaturase.faa"
OUTGROUP_FA["cart"] = TM / "outgroups_cart.faa"
FULL_FA = TM / "full_length" / "named_members.faa"
ALL_FA = list(MEMBER_FA.values()) + list(OUTGROUP_FA.values()) + [FULL_FA]
AA = set("ACDEFGHIKLMNPQRSTVWY")
# Code that makes step calls, absence calls or decoy-calibrated cutoffs.
MUST_NOT_READ = ["bin/score_refs_decoys.py", "bin/absence_check.py", "bin/make_min_scores.py",
                 "bin/merge_absence.py", "bin/assign_steps.py", "bin/build_step_hmms.py",
                 "pathway_census/absence.py", "pathway_census/defs.py",
                 "slurm/score_refs_decoys.sbatch"]


def _fasta(path):
    """Return [(header, sequence)] in file order. Fails on text before the first header."""
    recs, head, seq = [], None, []
    for ln in Path(path).read_text(encoding="utf-8").splitlines():
        if ln.startswith(">"):
            if head is not None:
                recs.append((head, "".join(seq)))
            head, seq = ln[1:], []
        elif ln.strip():
            assert head is not None, f"{path}: sequence before the first header"
            seq.append(ln.strip())
    if head is not None:
        recs.append((head, "".join(seq)))
    return recs


def _ids(path):
    return [h.split()[0] for h, _ in _fasta(path)]


def _rows():
    with open(PROV, encoding="utf-8", newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        assert r.fieldnames == COLS
        return list(r)


def test_files_parse():
    for p in ALL_FA:
        recs = _fasta(p)
        assert recs, p
        for h, s in recs:
            assert re.fullmatch(r"[A-Z0-9]{6,10}", h.split()[0]), (p, h)
            assert f"|{h.split()[0]}|" in h, (p, h)
            assert s and set(s) <= AA, (p, h)
    assert _rows()


def test_fasta_headers_are_unique():
    for p in ALL_FA:
        ids = _ids(p)
        assert len(ids) == len(set(ids)), p
        heads = [h for h, _ in _fasta(p)]
        assert len(heads) == len(set(heads)), p


def test_provenance_columns_and_roles():
    rows = _rows()
    keys = [(r["step_id"], r["accession"]) for r in rows]
    assert len(keys) == len(set(keys))
    for r in rows:
        assert r["role"] in ROLES, r
        assert r["step_id"] in STEPS | NEW_STEPS, r
        assert r["verified"] in ("true", "false"), r
        assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", r["retrieved_on"]), r
        assert re.match(r"[1-5]: ", r["protein_existence"]), r
        assert 1.0 <= float(r["annotation_score"]) <= 5.0, r
        assert "ECO:" in r["evidence_codes"], r
        assert re.fullmatch(r"\d+\.\d+\.\d+\.\d+", r["ec"]), r
        assert r["species"] and r["class"] and r["phylum"] and r["reason"], r
        assert all(re.fullmatch(r"\d+", p) for p in r["pmid"].split(",") if p), r
        assert (r["role"] == "named_member") == bool(r["region"]), r


def test_every_provenance_row_has_a_sequence_in_the_right_file():
    rows = _rows()
    for r in rows:
        fa = (MEMBER_FA if r["role"] == "named_member" else OUTGROUP_FA)[r["step_id"]]
        assert r["accession"] in _ids(fa), r
    # and the reverse: no sequence without a provenance row
    for role, files in (("named_member", MEMBER_FA), ("tree_outgroup", OUTGROUP_FA)):
        for step, fa in files.items():
            want = {r["accession"] for r in rows if r["role"] == role and r["step_id"] == step}
            assert set(_ids(fa)) == want, fa
    members = {r["accession"] for r in rows if r["role"] == "named_member"}
    assert set(_ids(FULL_FA)) == members


def test_cut_matches_region_and_full_length():
    full = {h.split()[0]: s for h, s in _fasta(FULL_FA)}
    region = {(r["step_id"], r["accession"]): r["region"] for r in _rows() if r["role"] == "named_member"}
    for step, fa in MEMBER_FA.items():
        for h, s in _fasta(fa):
            acc = h.split()[0]
            m = re.fullmatch(r"region=(\d+)-(\d+)", h.split()[1])
            assert m, h
            a, b = int(m.group(1)), int(m.group(2))
            assert region[(step, acc)] == f"{a}-{b}", h
            assert 1 <= a < b <= len(full[acc]), h
            assert s == full[acc][a - 1:b], h
    # the cyclase region ends before the synthase region starts (fused protein, two cuts)
    for (step, acc), reg in region.items():
        if step == "lycopene_cyclase":
            assert int(reg.split("-")[1]) < int(region[("phytoene_synthase", acc)].split("-")[0]), acc


def test_no_tree_member_is_a_reference_decoy_or_outgroup_of_the_same_step():
    for r in _rows():
        step, acc = r["step_id"], r["accession"]
        for other in (REFS / f"{step}.faa", REFS / "decoys" / f"{step}.faa", REFS / "outgroups" / f"{step}.faa"):
            if other.exists():
                assert acc not in _ids(other), (acc, other)


def test_roles_do_not_overlap():
    rows = _rows()
    for step in STEPS:
        mem = {r["accession"] for r in rows if r["step_id"] == step and r["role"] == "named_member"}
        og = {r["accession"] for r in rows if r["step_id"] == step and r["role"] == "tree_outgroup"}
        assert not (mem & og), step


def test_tree_members_are_not_in_the_pathway_definition():
    """The YAML holds references and decoys only. No tree member may be in it."""
    text = (ROOT / "pathways" / "carotenoid.yaml").read_text(encoding="utf-8")
    for acc in {r["accession"] for r in _rows()}:
        assert acc not in text, acc


def test_scoring_and_absence_code_does_not_name_the_folder():
    for rel in MUST_NOT_READ:
        p = ROOT / rel
        assert p.exists(), rel
        assert "tree_members" not in p.read_text(encoding="utf-8"), rel


def test_reference_dir_glob_does_not_reach_the_folder():
    """bin/score_refs_decoys.py lists steps with a non-recursive glob of refs/<pathway>/*.faa."""
    assert not [p for p in REFS.glob("*.faa") if TM in p.parents]
    src = (ROOT / "bin" / "score_refs_decoys.py").read_text(encoding="utf-8")
    assert "rglob" not in src and "**" not in src


# ---- second set: phytoene_desaturase and cart ------------------------------------------------

def _new_rows():
    return [r for r in _rows() if r["step_id"] in NEW_STEPS]


def test_new_steps_have_rows_and_first_set_is_unchanged():
    rows = _rows()
    old = [r for r in rows if r["step_id"] in STEPS]
    assert len(old) == 33 and rows[:33] == old           # the first set comes first and is complete
    assert {r["step_id"] for r in _new_rows()} == NEW_STEPS
    for step in NEW_STEPS:
        og = [r for r in _new_rows() if r["step_id"] == step and r["role"] == "tree_outgroup"]
        assert len(og) >= 3, step


def test_new_steps_roles_do_not_overlap():
    rows = _new_rows()
    for step in NEW_STEPS:
        mem = {r["accession"] for r in rows if r["step_id"] == step and r["role"] == "named_member"}
        og = {r["accession"] for r in rows if r["step_id"] == step and r["role"] == "tree_outgroup"}
        assert not (mem & og), step


def test_cart_has_no_named_member_and_no_member_file():
    assert not [r for r in _rows() if r["step_id"] == "cart" and r["role"] == "named_member"]
    assert not (TM / "cart.faa").exists()               # an empty file with a step name is a trap


def test_new_named_members_are_full_length_and_fungal():
    full = {h.split()[0]: s for h, s in _fasta(FULL_FA)}
    for r in _new_rows():
        if r["role"] != "named_member":
            continue
        assert r["region"] == f"1-{len(full[r['accession']])}", r
        assert r["phylum"].endswith("mycota") and r["class"].endswith("mycetes"), r
        assert "Full-length sequence" in r["reason"], r


def test_new_named_member_says_by_similarity_when_no_experiment():
    for r in _new_rows():
        if r["role"] == "named_member" and "ECO:0000269" not in r["evidence_codes"]:
            assert "by similarity" in r["evidence_codes"] and "No experiment" in r["reason"], r
            assert r["pmid"] == "", r


def test_new_outgroups_are_not_fungal_and_have_experimental_evidence():
    for r in _new_rows():
        if r["role"] != "tree_outgroup":
            continue
        assert not r["phylum"].endswith("mycota") and not r["class"].endswith("mycetes"), r
        assert "ECO:0000269" in r["evidence_codes"] and r["pmid"], r
        assert r["verified"] == "true", r


def test_new_outgroup_rows_give_a_passing_alignability():
    """The reason text gives the measured score and coverage. Both pass the rule of filter_decoys_for_tree.py."""
    pat = re.compile(r"best domain score (\d+\.\d) bits, HMM coverage (0\.\d+|1\.0) ")
    for r in _new_rows():
        if r["role"] == "tree_outgroup":
            m = pat.search(r["reason"])
            assert m, r
            assert float(m.group(1)) >= 10.0 and float(m.group(2)) >= 0.5, r


def test_new_outgroup_alignability_measured_again():
    """Measure with hmmsearch and the code of bin/filter_decoys_for_tree.py. Skipped without the tool or the HMM
    (results/ is not in git)."""
    import importlib.util
    import pytest
    from pathway_census import hits
    hmms = {s: ROOT / "results" / "hmms" / "carotenoid" / f"{s}.PROVISIONAL.hmm" for s in NEW_STEPS}
    if not hits.tool_available("hmmsearch") or not all(h.exists() for h in hmms.values()):
        pytest.skip("hmmsearch or the PROVISIONAL step HMMs are not present")
    spec = importlib.util.spec_from_file_location("filter_decoys_for_tree", ROOT / "bin" / "filter_decoys_for_tree.py")
    fd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fd)
    for step in sorted(NEW_STEPS):
        recs = list(fd.read_fasta(OUTGROUP_FA[step]))
        kept, table = fd.select(recs, fd.best_scores(hmms[step], OUTGROUP_FA[step]), 10.0, 0.5, fd.hmm_length(hmms[step]))
        assert len(kept) == len(recs), (step, table)
        want = {r["accession"]: r["reason"] for r in _new_rows() if r["step_id"] == step}
        for name, sc, cov, keep in table:
            assert f"best domain score {round(sc, 1)} bits, HMM coverage {round(cov, 3)} " in want[name], (name, sc, cov)


def test_strain_copies_of_references_are_not_tree_members():
    """S0EPU6 (same sequence as Q8X0Z0) and S0E9A6 (7 residues from A1KQY3) were fetched and left out."""
    for p in ALL_FA:
        assert not ({"S0EPU6", "S0E9A6"} & set(_ids(p))), p
    assert not ({"S0EPU6", "S0E9A6"} & {r["accession"] for r in _rows()})


def test_no_tree_member_sequence_equals_a_reference_or_decoy_sequence():
    for step in NEW_STEPS:
        other = set()
        for f in (REFS / f"{step}.faa", REFS / "decoys" / f"{step}.faa"):
            other |= {s for _, s in _fasta(f)}
        for f in (MEMBER_FA.get(step), OUTGROUP_FA[step]):
            if f is not None:
                assert not ({s for _, s in _fasta(f)} & other), f


def test_cutoff_and_absence_code_does_not_name_the_folder_any_file():
    """Wider than MUST_NOT_READ: every bin/, pathway_census/ and slurm/ file. Only the step 5 tree code may name it."""
    allowed = {"slurm/step5_prepare.sbatch", "bin/step5_place.py"}
    for sub in ("bin", "pathway_census", "slurm"):
        for p in sorted((ROOT / sub).glob("*")):
            rel = f"{sub}/{p.name}"
            if p.is_file() and p.suffix in (".py", ".sbatch", ".sh") and rel not in allowed:
                assert "tree_members" not in p.read_text(encoding="utf-8"), rel
