"""KLUSTR v1.0 verification against ACDC (Sect. 5 of the paper): the conditions of the comparison and the KLUSTR
values for them. Python 3.6 or newer, standard library only.

    python3 klustr_verification.py          every acdc/acdc_<free energies>_4x4.csv (and _band8.csv) found

Reads the output of acdc/acdc_verification.m and writes, next to it, one row per accepted condition:

  klustr_<free energies>_4x4.csv     T, CS, CA, CB, free acid of ACDC (1A), J_ACDC, J_KLUSTR, J_KLUSTR counted as in
                                     Algorithm 1, c1 of ACDC (all clusters with one A and at least one D), c1 of KLUSTR,
                                     share of 2A2D formed by collisions of two pairs in KLUSTR, sqrt(C_A C_B)/c_sat
  klustr_<free energies>_band8.csv   T, CS, CA, CB, J past 8A8D of ACDC and KLUSTR, J past 4A4D of ACDC and KLUSTR,
                                     J past 8A8D of KLUSTR with the collisions with larger clusters by iteration

A condition is accepted when, in the 4x4 run, ACDC converged with J > 1e-10 cm-3 s-1, the largest cluster 4A4D
collides more often than it evaporates (column stable), the pair has a barrier, C_A C_B < c_sat^2, and the pair is
the critical cluster, Eq. (2) for 2A2D, 3A3D and 4A4D, with C_A the free acid of ACDC; for the band to 8A8D also
when ACDC converged there with J > 1e-10. All ratios of acid to base are kept. KLUSTR as in ACDC: eta = 1,
gamma = -1.7, open system, a cluster counted at its first collision with a molecule past N pairs (unless stated).
Units: CS in s-1, concentrations in cm-3, J in cm-3 s-1, T in K.
"""
import csv, glob, os, sys
from math import exp, log, log10, sqrt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'python'))
import klustr

R = 8.314462618/4184                                   # kcal mol-1 K-1
SA_DMA = dict(gamma=-1.7, mA=98.08, mB=45.08, densA=1830.0, densB=680.0, eta=1.0, closed=False)


def free_energies(fe, T):
    """Formation free energies (kcal/mol, 1 atm) of 1A1D ... 4A4D at T: H and S (olenius2013) or G per temperature."""
    names = ['%dA%dD' % (n, n) for n in range(1, 5)]
    if fe == 'olenius2013':
        path, hs = os.path.join(HERE, 'acdc', 'SA-DMA_olenius2013_HS.txt'), True
    else:
        path, hs = os.path.join(HERE, 'acdc', fe, 'G%dK.txt' % round(T)), False
    g = {}
    for line in open(path):
        f = line.split()
        if f and f[0] in names:
            g[f[0]] = float(f[1]) - T*float(f[2])*1e-3 if hs else float(f[1])
    return [g[n] for n in names]


def rows(path):
    with open(path) as f:
        for r in csv.DictReader(f):
            yield r


def key(r):
    return (round(float(r['T'])), float(r['CS']), round(log10(float(r['CA'])), 3), round(log10(float(r['CB'])), 3))


def one_acid(names):
    """Clusters with one A and at least one D, e.g. 1A1D ... 1A4D."""
    return [c for c in names if c.startswith('1A') and c.endswith('D') and c != '1A']


def iterated(a, b, s, const, passes=12):
    """Reduced J past N pairs with the collisions of each cluster with larger clusters, taken from the previous pass
    (mass conserving), and the exact pair concentration (bisection); counted at the first collision with a molecule."""
    rA, rB, rho, sig = const['rhoA'], const['rhoB'], const['rho'], const['sig']
    N = len(sig) - 1
    kap = lambda a1, n: rA[n]*a1*rB[n]*b/(rA[n]*a1 + rB[n]*b + s*sig[n])
    big = [0.0]*(N + 1)
    for _ in range(passes):
        def pair(a1):
            L = 1 + kap(a1, 1) + s*sig[1] + big[1]
            return 2*a1*b/(L + sqrt(L*L + 4*rho[1][1]*a1*b))
        lo, hi = 0.0, a                                  # y1 = pair(a - y1), y1 - pair(a - y1) increases with y1
        for _ in range(100):
            mid = 0.5*(lo + hi)
            if mid - pair(a - mid) > 0:
                hi = mid
            else:
                lo = mid
        y1 = 0.5*(lo + hi); a1 = a - y1
        k = [0.0] + [kap(a1, n) for n in range(1, N + 1)]
        k[N] = rA[N]*a1 + rB[N]*b
        y = [0.0, y1]
        for n in range(2, N + 1):
            S = k[n-1]*y[n-1] + 0.5*sum(rho[i][n-i]*y[i]*y[n-i] for i in range(1, n))
            D = k[n] + sum(rho[i][n]*y[i] for i in range(1, n)) + s*sig[n] + big[n]
            y.append(2*S/(D + sqrt(D*D + 4*rho[n][n]*S)))
        big = [0.0] + [sum(rho[i][j]*y[j] for j in range(i + 1, N + 1)) for i in range(1, N + 1)]
    return k[N]*y[N] + sum((0.5 if i == j else 1.0)*rho[i][j]*y[i]*y[j]
                           for i in range(1, N + 1) for j in range(i, N + 1) if i + j > N)


def write(path, header, out):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in out:
            w.writerow([r[0], '%g' % r[1]] + ['%.6e' % v for v in r[2:]])


def report(name, n, bad, out, ab):
    print('%s: %d conditions, at most %d not converged at one temperature; accepted %d, of them %d with CA/CB <= 10 '
          '(%d with J_ACDC >= 1e-3)' % (name, n, max(bad.values()) if bad else 0, len(out),
                                       sum(1 for r in out if ab(r) <= 10), sum(1 for r in out if ab(r) <= 10 and r[5] >= 1e-3)))


def verify(fe):
    c4 = klustr.setup(4, first_collision=True, **SA_DMA)
    c4a = klustr.setup(4, first_collision=False, **SA_DMA)
    c8 = klustr.setup(8, first_collision=True, **SA_DMA)
    acc, out4, bad, n = {}, [], {}, 0
    path = os.path.join(HERE, 'acdc', 'acdc_%s_4x4.csv' % fe)
    with open(path) as f:
        if 'stable' not in next(csv.reader(f)):
            sys.exit('%s has no column stable: run acdc/acdc_verification.m of KLUSTR v1.0 again' % path)
    for r in rows(path):
        n += 1
        T, CS, CA, CB, J = (float(r[c]) for c in ('T', 'CS', 'CA', 'CB', 'J'))
        if r['converged'] != '1':
            bad[T] = bad.get(T, 0) + 1
            continue
        if not (J > 1e-10 and r['stable'] == '1'):
            continue
        g = [v/(R*T) for v in free_energies(fe, T)]       # dG_n/kB T
        sat = klustr.P0/(klustr.KB*T)*1e-6*exp(g[0])
        A1 = float(r['1A'])
        x = sqrt(A1*CB)/sat
        if x >= 1 or any(g[m-1] - (2*m - 1)*g[0] >= (2*m - 2)*log(x) for m in (2, 3, 4)):
            continue                                       # no barrier, or a larger cluster critical
        Jk, c = klustr.rate(CA, CB, CS, T, g[0], c4, clusters=True)
        Ja = klustr.rate(CA, CB, CS, T, g[0], c4a)
        c1a = sum(float(r[k]) for k in one_acid(r))
        bAB = c4['bABf']*sqrt(T)
        a, b, s, y1 = CA/sat, CB/sat, CS/(bAB*sat), c[0]/sat
        a1 = a - y1
        kap1 = c4['rhoA'][1]*a1*c4['rhoB'][1]*b/(c4['rhoA'][1]*a1 + c4['rhoB'][1]*b + s*c4['sig'][1])
        pp = 0.5*c4['rho'][1][1]*y1*y1
        out4.append((round(T), CS, CA, CB, A1, J, Jk, Ja, c1a, c[0], pp/(pp + kap1*y1), x))
        acc[key(r)] = (J, Jk, g[0], sat, bAB)
    write(os.path.join(HERE, 'acdc', 'klustr_%s_4x4.csv' % fe),
          ['T', 'CS', 'CA', 'CB', 'A1_ACDC', 'J_ACDC', 'J_KLUSTR', 'J_KLUSTR_alg1', 'c1_ACDC', 'c1_KLUSTR',
           'share_2A2D_pair_collisions', 'sqrt_CACB_over_csat'], out4)
    report('%s 4x4' % fe, n, bad, out4, lambda r: r[2]/r[3])

    path = os.path.join(HERE, 'acdc', 'acdc_%s_band8.csv' % fe)
    if not os.path.isfile(path):
        return
    out8, bad, n = [], {}, 0
    for r in rows(path):
        n += 1
        T, CS, CA, CB, J = (float(r[c]) for c in ('T', 'CS', 'CA', 'CB', 'J'))
        if r['converged'] != '1':
            bad[T] = bad.get(T, 0) + 1
            continue
        if J <= 1e-10 or key(r) not in acc:
            continue
        J4, J4k, g1, sat, bAB = acc[key(r)]
        J8k = klustr.rate(CA, CB, CS, T, g1, c8)
        J8i = bAB*sat*sat*iterated(CA/sat, CB/sat, CS/(bAB*sat), c8)
        out8.append((round(T), CS, CA, CB, J, J8k, J4, J4k, J8i))
    write(os.path.join(HERE, 'acdc', 'klustr_%s_band8.csv' % fe),
          ['T', 'CS', 'CA', 'CB', 'J8_ACDC', 'J8_KLUSTR', 'J4_ACDC', 'J4_KLUSTR', 'J8_KLUSTR_iterated'], out8)
    report('%s band8' % fe, n, bad, [(r[0], r[1], r[2], r[3], None, r[4]) for r in out8], lambda r: r[2]/r[3])


if __name__ == '__main__':
    sets = sorted(os.path.basename(p)[5:-8] for p in glob.glob(os.path.join(HERE, 'acdc', 'acdc_*_4x4.csv')))
    if not sets:
        sys.exit('no acdc/acdc_<free energies>_4x4.csv: run acdc/acdc_verification.m in MATLAB first')
    for fe in sets:
        verify(fe)
