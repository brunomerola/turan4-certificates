#!/usr/bin/env python3
"""Negative controls for sk_check_keys.py: corrupted copies of one key file (written to the scratch dir given as argv[1],
together with links to the cert files) must FAIL: two flags swapped, one |Aut| changed, one forbidden graph changed,
one factor-row count changed."""
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, "..", "certificates", "sigma_catalogue", "K5lt")
SRC = os.path.join(CERT, "cosig3_s3_K5lt_f.keys.json")
out = sys.argv[1]
os.makedirs(out, exist_ok=True)
for f in os.listdir(CERT):
    if f.startswith("cosig3_s3_K5lt_f.cosig3"):
        shutil.copy2(os.path.join(CERT, f), out)
cases = {}
K = json.load(open(SRC))
k1 = json.loads(json.dumps(K)); fl = k1["keys"][40]["flags"]; fl[3], fl[4] = fl[4], fl[3]; cases["swap"] = k1
k2 = json.loads(json.dumps(K)); k2["keys"][20]["aut_size"] += 1; cases["aut"] = k2
k3 = json.loads(json.dumps(K)); k3["predicate"]["forbidden"] = ["123 124 125 134 135 145 234 245"]; cases["forb"] = k3
k4 = json.loads(json.dumps(K)); k4["keys"][7]["factor_rows"] += 1; cases["rows"] = k4
ok = True
for name, data in cases.items():
    p = os.path.join(out, f"neg_{name}.keys.json")
    data["certificate"] = "cosig3_s3_K5lt_f"
    json.dump(data, open(p, "w"))
    rc = subprocess.run([sys.executable, "-B", os.path.join(HERE, "sk_check_keys.py"), p], capture_output=True, text=True)
    err = json.loads(rc.stdout.splitlines()[0])["errors"]
    print(f"{name}: exit {rc.returncode} ({'FAIL as required' if rc.returncode else 'PASSED -- CONTROL BROKEN'}) {err[:1]}")
    ok &= rc.returncode != 0
print("NEGATIVE CONTROLS OK" if ok else "NEGATIVE CONTROLS BROKEN")
sys.exit(0 if ok else 1)
