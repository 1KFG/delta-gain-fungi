"""Pathway definition loader and validator (spec section 5.2)."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import yaml

COPY_EXPECT = ("single", "multi")
ROLES = ("catalytic", "auxiliary")


@dataclass
class Ref:
    accession: str | None
    species: str
    pmid: str | None
    verified: bool
    gene_name: str | None = None
    region: tuple[int, int] | None = None  # 1-based inclusive, inside the sequence


@dataclass
class Step:
    id: str
    name: str
    ec: str | None
    refs: list[Ref]
    domains_any: list[str]  # any-of: one hit among these is enough
    copy_expect: Literal["single", "multi"]
    role: Literal["catalytic", "auxiliary"]
    fused_with: list[str] = field(default_factory=list)
    optional: bool = False
    route: str | None = None
    reactions: list[str] = field(default_factory=list)
    decoys: list[Ref] = field(default_factory=list)
    lineage_note: str | None = None


@dataclass
class Route:
    id: str
    name: str
    step_ids: list[str]
    note: str | None = None


@dataclass
class Pathway:
    id: str
    steps: list[Step]
    verified: bool
    routes: list[Route] = field(default_factory=list)
    # Directory with <step id>.faa reference files. Set by load_pathway. None: no region check.
    refs_dir: Path | None = None


def _check_keys(d: dict, allowed: set[str], where: str) -> None:
    if not isinstance(d, dict):
        raise ValueError(f"{where}: expected a mapping, got {type(d).__name__}")
    unknown = sorted(set(d) - allowed)
    if unknown:
        raise ValueError(f"{where}: unknown key(s) {unknown}")


def _str_or_none(v):
    return None if v is None else str(v)


def _bool(d: dict, key: str, default: bool, where: str) -> bool:
    """Return d[key] as a bool. Only a real YAML bool is accepted.

    A string such as "false" or "no" would be truthy under bool(); reject it.
    """
    v = d.get(key, default)
    if not isinstance(v, bool):
        raise ValueError(f"{where}: '{key}' must be a YAML bool (true/false), "
                         f"got {type(v).__name__} {v!r}")
    return v


def _str_list(d: dict, key: str, where: str) -> list[str]:
    """Return d[key] as a list of strings. A scalar is an error, not split into chars."""
    v = d.get(key)
    if v is None:
        return []
    if not isinstance(v, list):
        raise ValueError(f"{where}: '{key}' must be a list, got {type(v).__name__} {v!r}")
    return [str(x) for x in v]


def _load_region(d: dict, where: str) -> tuple[int, int] | None:
    v = d.get("region")
    if v is None:
        return None
    ok = (isinstance(v, list) and len(v) == 2
          and all(isinstance(x, int) and not isinstance(x, bool) for x in v))
    if not ok:
        raise ValueError(f"{where}: 'region' must be [start, end] of two integers, got {v!r}")
    if v[0] < 1 or v[1] < v[0]:
        raise ValueError(f"{where}: 'region' {v!r} must satisfy 1 <= start <= end")
    return (v[0], v[1])


def _fasta_lengths(path: Path) -> dict[str, int]:
    """Map every whitespace/pipe separated header token to the sequence length."""
    out: dict[str, int] = {}
    toks: list[str] = []
    n = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                for t in toks:
                    out[t] = n
                toks = [t for t in line[1:].replace("|", " ").split() if t]
                n = 0
            else:
                n += len(line)
    for t in toks:
        out[t] = n
    return out


def _load_ref(d: dict, where: str) -> Ref:
    _check_keys(d, {"accession", "species", "pmid", "verified", "gene_name", "region"}, where)
    return Ref(
        accession=_str_or_none(d.get("accession")),
        species=str(d.get("species", "")),
        pmid=_str_or_none(d.get("pmid")),
        verified=_bool(d, "verified", False, where),
        gene_name=_str_or_none(d.get("gene_name")),
        region=_load_region(d, where),
    )


def _load_step(d: dict, i: int) -> Step:
    where = f"step[{i}]"
    if isinstance(d, dict) and "domains" in d:
        raise ValueError(f"{where}: key 'domains' was renamed, use domains_any")
    _check_keys(d, {"id", "name", "ec", "refs", "domains_any", "copy_expect",
                    "role", "fused_with", "optional", "route", "reactions",
                    "decoys", "lineage_note"}, where)
    for req in ("id", "name", "copy_expect", "role"):
        if req not in d:
            raise ValueError(f"{where}: missing required key '{req}'")
    where = f"step '{d['id']}'"
    refs = [_load_ref(r, f"{where} ref[{j}]") for j, r in enumerate(d.get("refs") or [])]
    decoys = [_load_ref(r, f"{where} decoy[{j}]") for j, r in enumerate(d.get("decoys") or [])]
    return Step(
        id=str(d["id"]),
        name=str(d["name"]),
        ec=_str_or_none(d.get("ec")),
        refs=refs,
        domains_any=_str_list(d, "domains_any", where),
        copy_expect=d["copy_expect"],
        role=d["role"],
        fused_with=_str_list(d, "fused_with", where),
        optional=_bool(d, "optional", False, where),
        route=_str_or_none(d.get("route")),
        reactions=_str_list(d, "reactions", where),
        decoys=decoys,
        lineage_note=_str_or_none(d.get("lineage_note")),
    )


def _load_route(d: dict, i: int) -> Route:
    where = f"route[{i}]"
    _check_keys(d, {"id", "name", "step_ids", "note"}, where)
    for req in ("id", "name"):
        if req not in d:
            raise ValueError(f"{where}: missing required key '{req}'")
    return Route(id=str(d["id"]), name=str(d["name"]),
                 step_ids=_str_list(d, "step_ids", f"route '{d['id']}'"),
                 note=_str_or_none(d.get("note")))


def _region_errors(p: Pathway) -> list[str]:
    """Check ref regions against the reference FASTA files.

    No refs directory (None or absent): skip. Directory present: a ref with a region needs an
    accession, a FASTA file for the step and the accession in it, and the region must end
    inside the sequence.
    """
    if p.refs_dir is None or not Path(p.refs_dir).is_dir():
        return []
    errs: list[str] = []
    for s in p.steps:
        lengths: dict[str, int] | None = None
        for r in s.refs:
            if r.region is None:
                continue
            if not r.accession:
                errs.append(f"step '{s.id}': region needs an accession")
                continue
            if lengths is None:
                faa = Path(p.refs_dir) / f"{s.id}.faa"
                lengths = _fasta_lengths(faa) if faa.exists() else {}
                have_file = faa.exists()
            if not have_file:
                errs.append(f"step '{s.id}' ref '{r.accession}': region cannot be checked, "
                            f"no FASTA file {s.id}.faa in the refs directory")
            elif r.accession not in lengths:
                errs.append(f"step '{s.id}' ref '{r.accession}': region cannot be checked, "
                            f"accession not found in {s.id}.faa")
            elif r.region[1] > lengths[r.accession]:
                errs.append(f"step '{s.id}' ref '{r.accession}': region {list(r.region)} "
                            f"lies outside the sequence (length {lengths[r.accession]})")
    return errs


def load_pathway(path: Path, refs_dir: Path | None = None) -> Pathway:
    """Load a pathway YAML.

    `refs_dir` defaults to <yaml dir>/../refs/<pathway id>. This assumes the YAML sits in
    `<root>/pathways/` and the references in `<root>/refs/<id>/`. Pass `refs_dir` for any
    other layout. The directory is stored on the Pathway and used by `validate` for the
    region check; if it does not exist the region check is skipped.
    """
    path = Path(path)
    with open(path, encoding="utf-8") as fh:
        d = yaml.safe_load(fh)
    _check_keys(d, {"id", "verified", "steps", "routes"}, "pathway")
    if "id" not in d:
        raise ValueError("pathway: missing required key 'id'")
    steps = [_load_step(s, i) for i, s in enumerate(d.get("steps") or [])]
    routes = [_load_route(r, i) for i, r in enumerate(d.get("routes") or [])]
    p = Pathway(id=str(d["id"]), steps=steps,
                verified=_bool(d, "verified", False, "pathway"), routes=routes)
    p.refs_dir = Path(refs_dir) if refs_dir is not None else path.parent.parent / "refs" / p.id
    return p


def validate(p: Pathway) -> list[str]:
    errs: list[str] = []
    ids = [s.id for s in p.steps]
    for sid in sorted({i for i in ids if ids.count(i) > 1}):
        errs.append(f"duplicate step id '{sid}'")
    known = set(ids)
    for s in p.steps:
        if s.copy_expect not in COPY_EXPECT:
            errs.append(f"step '{s.id}': invalid copy_expect '{s.copy_expect}'")
        if s.role not in ROLES:
            errs.append(f"step '{s.id}': invalid role '{s.role}'")
        if p.verified and not s.optional and not s.refs:
            errs.append(f"step '{s.id}': pathway verified=true but non-optional "
                        f"step has no refs")
        for r in s.refs:
            if not r.accession and not r.pmid:
                errs.append(f"step '{s.id}': ref has neither accession nor pmid")
            if p.verified and not r.verified:
                errs.append(f"step '{s.id}': pathway verified=true but ref "
                            f"'{r.accession or r.pmid}' is not verified")
        for t in s.fused_with:
            if t not in known:
                errs.append(f"step '{s.id}': fused_with target '{t}' is not a step id")
    route_ids = [r.id for r in p.routes]
    for rid in sorted({i for i in route_ids if route_ids.count(i) > 1}):
        errs.append(f"duplicate route id '{rid}'")
    for r in p.routes:
        for sid in r.step_ids:
            if sid not in known:
                errs.append(f"route '{r.id}': step '{sid}' is not a step id")
    for s in p.steps:
        if s.route is not None and s.route not in route_ids:
            errs.append(f"step '{s.id}': route '{s.route}' is not declared")
        elif s.route is not None and not any(
                r.id == s.route and s.id in r.step_ids for r in p.routes):
            errs.append(f"step '{s.id}': route '{s.route}' does not list the step "
                        f"in its step_ids")
        for rx in sorted({x for x in s.reactions if s.reactions.count(x) > 1}):
            errs.append(f"step '{s.id}': duplicate reaction id '{rx}'")
        ref_acc = {r.accession for r in s.refs if r.accession}
        dacc = [dc.accession for dc in s.decoys if dc.accession]
        for a in sorted({x for x in dacc if dacc.count(x) > 1}):
            errs.append(f"step '{s.id}': duplicate decoy accession '{a}'")
        for dc in s.decoys:
            if dc.accession and dc.accession in ref_acc:
                errs.append(f"step '{s.id}': decoy accession '{dc.accession}' is also "
                            f"a reference of the step")
    errs.extend(_region_errors(p))
    return errs
