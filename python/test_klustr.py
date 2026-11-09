"""Test of the Python version: python test_klustr.py
1. the 50 reference cases of test_reference.txt (N = 4, 8 and 16, open and closed systems, both countings, eta = 1
   and 3) to a relative difference of 1e-12;
2. eta scales out: J with eta = 3 equals 3 J with eta = 1 and CS_A/3;
3. with a barrier, C_A C_B < c_sat^2, open and closed systems give the same J;
4. with clusters=True, the same J and N cluster concentrations, the pair at most at c_sat;
5. the documented condition of klustr.py gives J = 2.235716e+02 cm-3 s-1.
"""
import os, sys, math, random
import klustr

here = os.path.dirname(os.path.abspath(__file__))
fails = 0

def check(name, ok):
    global fails
    print(f"{'ok  ' if ok else 'FAIL'} {name}")
    fails += not ok

# 1. reference values
worst = 0.0
for line in open(os.path.join(here, 'test_reference.txt')):
    if line.startswith('#') or not line.strip():
        continue
    v = line.split()
    N, closed, first, eta = int(v[0]), v[1] == '1', v[2] == '1', float(v[3])
    a, b, cs, T, g, Jref = map(float, v[4:])
    c = klustr.setup(N, -1.7, 98.08, 45.08, 1830.0, 680.0, eta, closed=closed, first_collision=first)
    J = klustr.rate(a, b, cs, T, g, c)
    worst = max(worst, abs(J/Jref - 1) if Jref > 0 else abs(J))
check(f'reference values (largest relative difference {worst:.1e})', worst < 1e-12)

# 2. eta scales out, 3. open and closed systems agree with a barrier, 4. cluster concentrations
rng = random.Random(1)
c1 = klustr.setup(8, -1.7, 98.08, 45.08, 1830.0, 680.0, 1.0)
c3 = klustr.setup(8, -1.7, 98.08, 45.08, 1830.0, 680.0, 3.0)
c1open = klustr.setup(8, -1.7, 98.08, 45.08, 1830.0, 680.0, 1.0, closed=False)
worst_eta = worst_open = 0.0
ok_clusters = True
for _ in range(1000):
    T = rng.uniform(220, 320); g = -rng.uniform(15, 35)
    sat = klustr.P0/(klustr.KB*T)*1e-6*math.exp(g)
    a, b, cs = sat*10**rng.uniform(-4, 4), sat*10**rng.uniform(-4, 6), 10**rng.uniform(-5, 0)
    J3, J1 = klustr.rate(a, b, cs, T, g, c3), klustr.rate(a, b, cs/3, T, g, c1)
    if J1 > 0:
        worst_eta = max(worst_eta, abs(J3/(3*J1) - 1))
    J, c = klustr.rate(a, b, cs, T, g, c1, clusters=True)
    ok_clusters &= J == klustr.rate(a, b, cs, T, g, c1) and len(c) == 8 and c[0] <= sat*(1 + 1e-15)
    if a*b < sat*sat:
        Jc, Jo = klustr.rate(a, b, cs, T, g, c1), klustr.rate(a, b, cs, T, g, c1open)
        if Jc > 0:
            worst_open = max(worst_open, abs(Jo/Jc - 1))
check(f'eta scales out (largest relative difference {worst_eta:.1e})', worst_eta < 1e-12)
check(f'open = closed with a barrier (largest relative difference {worst_open:.1e})', worst_open < 1e-12)
check('cluster concentrations returned with the same J', ok_clusters)

# 5. the documented condition (klustr.py __main__)
cd = klustr.setup(N=8, gamma=-1.7, mA=98.08, mB=45.08, densA=1830.0, densB=680.0, eta=3.0)
Jd = klustr.rate(1e6, 1e7, 1e-3, 278.0, (-24.6465441592 - 278.0*(-31.008553e-3))/(1.987204259e-3*278.0), cd)
check(f'documented condition J = {Jd:.6e} cm-3 s-1 (expected 2.235716e+02)', f'{Jd:.6e}' == '2.235716e+02')

sys.exit(1 if fails else 0)
