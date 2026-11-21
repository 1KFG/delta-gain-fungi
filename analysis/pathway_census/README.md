# pathway_census

Shared pathway census package (carotenoid pilot). Task 1: definition loader (`pathway_census.defs`).

Run tests with the system Python 3.12 (default `python3` is 3.9):

    /usr/bin/python3.12 -m pytest analysis/pathway_census/tests -v

`pytest` and `pyyaml` were already present under `/usr/bin/python3.12`, so no venv was needed.
If they go missing, create `analysis/pathway_census/.venv` (git-ignored) and install `pyyaml pytest` there.
