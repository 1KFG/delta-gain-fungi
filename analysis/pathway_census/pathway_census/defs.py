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


@dataclass
class Step:
    id: str
    name: str
    ec: str | None
    refs: list[Ref]
    domains: list[str]
    copy_expect: Literal["single", "multi"]
    role: Literal["catalytic", "auxiliary"]
    fused_with: list[str] = field(default_factory=list)
    optional: bool = False


@dataclass
class Pathway:
    id: str
    steps: list[Step]
    verified: bool


def _check_keys(d: dict, allowed: set[str], where: str) -> None:
    if not isinstance(d, dict):
        raise ValueError(f"{where}: expected a mapping, got {type(d).__name__}")
    unknown = sorted(set(d) - allowed)
    if unknown:
        raise ValueError(f"{where}: unknown key(s) {unknown}")


def _str_or_none(v):
    return None if v is None else str(v)


def _load_ref(d: dict, where: str) -> Ref:
    _check_keys(d, {"accession", "species", "pmid", "verified"}, where)
    return Ref(
        accession=_str_or_none(d.get("accession")),
        species=str(d.get("species", "")),
        pmid=_str_or_none(d.get("pmid")),
        verified=bool(d.get("verified", False)),
    )


def _load_step(d: dict, i: int) -> Step:
    where = f"step[{i}]"
    _check_keys(d, {"id", "name", "ec", "refs", "domains", "copy_expect",
                    "role", "fused_with", "optional"}, where)
    for req in ("id", "name"):
        if req not in d:
            raise ValueError(f"{where}: missing required key '{req}'")
    where = f"step '{d['id']}'"
    refs = [_load_ref(r, f"{where} ref[{j}]") for j, r in enumerate(d.get("refs") or [])]
    return Step(
        id=str(d["id"]),
        name=str(d["name"]),
        ec=_str_or_none(d.get("ec")),
        refs=refs,
        domains=[str(x) for x in d.get("domains") or []],
        copy_expect=d.get("copy_expect", "single"),
        role=d.get("role", "catalytic"),
        fused_with=[str(x) for x in d.get("fused_with") or []],
        optional=bool(d.get("optional", False)),
    )


def load_pathway(path: Path) -> Pathway:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    _check_keys(d, {"id", "verified", "steps"}, "pathway")
    if "id" not in d:
        raise ValueError("pathway: missing required key 'id'")
    steps = [_load_step(s, i) for i, s in enumerate(d.get("steps") or [])]
    return Pathway(id=str(d["id"]), steps=steps, verified=bool(d.get("verified", False)))


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
        for r in s.refs:
            if not r.accession and not r.pmid:
                errs.append(f"step '{s.id}': ref has neither accession nor pmid")
            if p.verified and not r.verified:
                errs.append(f"step '{s.id}': pathway verified=true but ref "
                            f"'{r.accession or r.pmid}' is not verified")
        for t in s.fused_with:
            if t not in known:
                errs.append(f"step '{s.id}': fused_with target '{t}' is not a step id")
    return errs
