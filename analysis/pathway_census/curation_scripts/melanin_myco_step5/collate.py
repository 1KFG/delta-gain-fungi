#!/usr/bin/env python3.12
"""Collate the step 5 runs of the other melanin and mycosporine steps into one table (one-off helper, not tested code).

For each run folder: runs bin/step5_run_report.py (writes <run>/report/), reads summary.json, the reps summary and the extract
summary, and writes results/step5/other_steps_summary.tsv. Run from analysis/pathway_census:
  /usr/bin/python3.12 curation_scripts/melanin_myco_step5/collate.py results/step5/other_steps_summary.tsv RUN_DIR [RUN_DIR ...]
RUN_DIR may be <pathway>/<step>[_idNN]; a folder without SUMMARY.tsv is listed as `not finished`.
"""
import csv, json, subprocess, sys
from pathlib import Path

out = Path(sys.argv[1])
cols = ["run", "step", "identity", "candidates", "domains_extracted", "extract_failures", "clusters", "representatives", "leaves",
        "n_refs", "n_decoys", "n_named", "n_outgroups", "refs_clade_free", "smallest_rule_target", "target_size", "target_reps",
        "call_target", "call_decoy", "call_unassigned", "call_other", "support_min", "assigned_proteins", "genomes_assigned",
        "genomes_searched", "copies_1", "copies_2plus", "candidates_equal_placements"]
rows = []
for d in map(Path, sys.argv[2:]):
    d = Path("results/step5") / d if not d.exists() else d
    if not (d / "SUMMARY.tsv").exists():
        rows.append({"run": d.name, "step": "not finished"})
        continue
    subprocess.run(["/usr/bin/python3.12", "bin/step5_run_report.py", "--run", str(d), "--genomes", "results/genomes.tsv.zst"],
                   check=True, capture_output=True)
    s = json.loads((d / "report" / "summary.json").read_text())
    cr = json.loads((d / "clade_report.json").read_text())
    if (d / "extract_summary.json").exists():
        ex = json.loads((d / "extract_summary.json").read_text())
    else:   # joint SDR step folders: domains.faa holds this step's candidates; extract ran in sdr_joint
        n_dom = sum(1 for l in open(d / "domains.faa") if l.startswith(">"))
        ex = {"domains_extracted": n_dom, "failures": s["n_candidates_step_rows"] - n_dom}
    rp = json.loads((d / "reps_summary.json").read_text())
    c = s["calls"]
    cp = s["copies_per_genome"]
    name = d.name
    ident = name.split("_id")[-1] if "_id" in name else ""
    rows.append({"run": name, "step": s["step_id"], "identity": ident, "candidates": s["n_candidates_step_rows"],
                 "domains_extracted": ex["domains_extracted"], "extract_failures": ex["failures"], "clusters": rp["clusters"],
                 "representatives": rp["representatives_selected"], "leaves": cr["n_leaves"], "n_refs": cr["n_references"],
                 "n_decoys": cr["n_decoys"], "n_named": cr["n_named_members"], "n_outgroups": cr["n_outgroups"],
                 "refs_clade_free": s["refs_form_clade_free_of_labelled_leaves"], "smallest_rule_target": s["smallest_rule_target_size"],
                 "target_size": s["target_size"], "target_reps": s["target_representatives"], "call_target": c.get("target", 0),
                 "call_decoy": c.get("decoy", 0), "call_unassigned": c.get("unassigned", 0),
                 "call_other": sum(v for k, v in c.items() if k not in ("target", "decoy", "unassigned")),
                 "support_min": s["support_target_min"], "assigned_proteins": s["assigned_proteins"],
                 "genomes_assigned": s["genomes_with_assigned"], "genomes_searched": s["genomes_searched"],
                 "copies_1": cp.get("1", 0), "copies_2plus": sum(v for k, v in cp.items() if int(k) >= 2),
                 "candidates_equal_placements": s["candidates_equal_placements"]})
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print(out, len(rows), "rows")
