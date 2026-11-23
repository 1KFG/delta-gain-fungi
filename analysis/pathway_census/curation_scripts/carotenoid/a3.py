from h import *
import time, sys
pm = sys.argv[2:]
with open(sys.argv[1],"w") as fh:
    for p in pm:
        for attempt in range(3):
            try:
                d = pubmed_abstract(p)
                fh.write(f"=== {p} ({d['year']}) {d['title']}\n{d['abstract']}\n\n"); break
            except KbFetchError as x:
                err = x; time.sleep(5)
        else:
            fh.write(f"=== {p} ERROR {err}\n\n"); print(p, "ERROR", err)
        time.sleep(1.0)
