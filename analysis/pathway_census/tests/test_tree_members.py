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
        assert r["step_id"] in STEPS, r
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
