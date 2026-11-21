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


def _load_ref(d: dict, where: str) -> Ref:
    _check_keys(d, {"accession", "species", "pmid", "verified"}, where)
    return Ref(
        accession=_str_or_none(d.get("accession")),
        species=str(d.get("species", "")),
        pmid=_str_or_none(d.get("pmid")),
        verified=_bool(d, "verified", False, where),
    )


def _load_step(d: dict, i: int) -> Step:
    where = f"step[{i}]"
    _check_keys(d, {"id", "name", "ec", "refs", "domains", "copy_expect",
                    "role", "fused_with", "optional"}, where)
    for req in ("id", "name", "copy_expect", "role"):
        if req not in d:
            raise ValueError(f"{where}: missing required key '{req}'")
    where = f"step '{d['id']}'"
    refs = [_load_ref(r, f"{where} ref[{j}]") for j, r in enumerate(d.get("refs") or [])]
    return Step(
        id=str(d["id"]),
        name=str(d["name"]),
        ec=_str_or_none(d.get("ec")),
        refs=refs,
        domains=_str_list(d, "domains", where),
        copy_expect=d["copy_expect"],
        role=d["role"],
        fused_with=_str_list(d, "fused_with", where),
        optional=_bool(d, "optional", False, where),
    )


def load_pathway(path: Path) -> Pathway:
    with open(path, encoding="utf-8") as fh:
        d = yaml.safe_load(fh)
    _check_keys(d, {"id", "verified", "steps"}, "pathway")
    if "id" not in d:
        raise ValueError("pathway: missing required key 'id'")
    steps = [_load_step(s, i) for i, s in enumerate(d.get("steps") or [])]
    return Pathway(id=str(d["id"]), steps=steps,
                   verified=_bool(d, "verified", False, "pathway"))


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
    return errs
