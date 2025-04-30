#!/usr/bin/env python3
import os
import subprocess
import re
import matplotlib.pyplot as plt

# ──────────────────────────────────────────────────────────────────────────────
# Configuration
apps    = ["app3" ,"app4","app5"]
exe     = "./L1simulate"
outdir  = "logs"
# Default cache params:
default_s = 6   # 2^6 = 64 sets
default_E = 2   # 2-way
default_b = 5   # 2^5 = 32-byte blocks
# Parameter sweeps (exponents)
s_vals = [1,2,3,4, 5, 6, 7]
E_vals = [1, 2,3, 4,5, 8]
b_vals = [1,2,3,4, 5, 6, 7]
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs(outdir, exist_ok=True)

def run_sim(trace, s, E, b, tag):
    out = f"{outdir}/{trace}_{tag}.txt"
    cmd = [exe, "-t", trace, "-s", str(s), "-E", str(E), "-b", str(b), "-o", out]
    print(f"→ {trace}: running", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return out

def parse_max(fname):
    patt = re.compile(r"Total Execution Cycles:\s*(\d+)")
    mx = 0
    with open(fname) as f:
        for L in f:
            m = patt.search(L)
            if m:
                mx = max(mx, int(m.group(1)))
    if mx == 0:
        raise RuntimeError(f"No cycles found in {fname}")
    return mx

def sweep(param, values):
    """Return dict app→(xlabels, yvalues)."""
    results = {}
    for app in apps:
        xs, ys = [], []
        for v in values:
            if param == "s":
                s, E, b = v, default_E, default_b
            elif param == "E":
                s, E, b = default_s, v, default_b
            else:  # "b"
                s, E, b = default_s, default_E, v

            tag = f"{param}{v}"
            fname = run_sim(app, s, E, b, tag)
            mx    = parse_max(fname)

            xs.append(f"2^{v}")
            ys.append(mx)
            print(f"   {app} @ {param}={v} → {mx}")
        results[app] = (xs, ys)
    return results

if __name__ == "__main__":
    sweeps = [("s", s_vals), ("E", E_vals), ("b", b_vals)]
    for param, vals in sweeps:
        data = sweep(param, vals)

        plt.figure()
        for app in apps:
            xs, ys = data[app]
            plt.plot(xs, ys, marker='o', label=app)
        plt.title(f"Max Exec Time vs {param.upper()} (default s={default_s}, E={default_E}, b={default_b})")
        plt.xlabel(f"{param.upper()} (2^x)")
        plt.ylabel("Max Execution Time (cycles)")
        plt.grid(True)
        plt.legend()
        out_png = f"plot_{param}.png"
        plt.savefig(out_png)
        print(f"→ Saved {out_png}\n")
