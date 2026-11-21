"""Run every check in this folder. Exit status is nonzero if any assertion fails.
Requires NumPy and SciPy. Usage: python3 hamiltonian/check_all.py"""
import os, subprocess, sys, time
here = os.path.dirname(os.path.abspath(__file__))
scripts = ["check_su2.py", "plaquette.py", "validate_ladder.py", "validate_ladder2.py", "ladder_gap.py"]
failed = []
for s in scripts:
    t = time.time()
    r = subprocess.run([sys.executable, "-W", "error", os.path.join(here, s)], capture_output=True, text=True)
    status = "ok" if r.returncode == 0 else "FAILED"
    print(f"{status:6s} {s:22s} {time.time() - t:5.1f} s")
    if r.returncode:
        failed.append(s)
        print(r.stdout[-2000:], r.stderr[-2000:])
print("all checks passed" if not failed else f"failed: {failed}")
sys.exit(1 if failed else 0)
