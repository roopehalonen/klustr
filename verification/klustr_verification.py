"""KLUSTR v1.0 against ACDC (Sect. 5 of the paper). Python 3.6 or newer, standard library only.

Reads acdc/acdc_<free energies>_<set>.csv of acdc_verification.m and writes acdc/klustr_<free energies>_<set>.csv
for the conditions of the comparison: ACDC converged with J > 1e-10 cm-3 s-1, 4A4D stable, and the pair the critical
cluster with a barrier (Eq. 2 for 2A2D to 4A4D, with the free acid 1A of ACDC). The band to 8A8D uses the same
conditions. KLUSTR as ACDC: eta = 1, gamma = -1.7, open system, counted at the first collision past N pairs.
"""
import csv, os, sys
from math import exp, log, sqrt

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', 'python'))
import klustr

R = 1.98720425864e-3                                    # kcal mol-1 K-1
SA_DMA = dict(gamma=-1.7, mA=98.08, mB=45.08, densA=1830.0, densB=680.0, eta=1.0, closed=False)


def free_energies(fe, T):
    """dG/kB T of 1A1D ... 4A4D at 1 atm, from H and S (olenius2013) or from G at each temperature."""
    f, hs = ('SA-DMA_olenius2013_HS.txt', True) if fe == 'olenius2013' else ('%s/G%dK.txt' % (fe, T), False)
    g = {}
    for line in open(os.path.join(here, 'acdc', f)):
        w = line.split()
        if len(w) > 1 and w[0] in ('1A1D', '2A2D', '3A3D', '4A4D'):
            g[w[0]] = (float(w[1]) - T*float(w[2])*1e-3 if hs else float(w[1]))/(R*T)
    return [g['%dA%dD' % (n, n)] for n in range(1, 5)]


def iterated(a, b, s, c):
    """Reduced J past N pairs with the collisions of each cluster with larger clusters from the previous pass,
    the exact pair concentration and counting at the first collision (Table 1)."""
    rA, rB, rho, sig = c['rhoA'], c['rhoB'], c['rho'], c['sig']
    N = len(sig) - 1
    kap = lambda a1, n: rA[n]*a1*rB[n]*b/(rA[n]*a1 + rB[n]*b + s*sig[n])
    big = [0.0]*(N + 1)
    for _ in range(12):
        lo, hi = 0.0, a                                 # pair balance with a1 = a - y1, by bisection
        for _ in range(100):
            y1 = 0.5*(lo + hi); a1 = a - y1
            M = 1 + kap(a1, 1) + s*sig[1] + big[1]
            if y1 > 2*a1*b/(M + sqrt(M*M + 4*rho[1][1]*a1*b)):
                hi = y1
            else:
                lo = y1
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


c4 = klustr.setup(4, first_collision=True, **SA_DMA)
c4p = klustr.setup(4, **SA_DMA)                         # counted as in Algorithm 1
c8 = klustr.setup(8, first_collision=True, **SA_DMA)
for fe in ('olenius2013', 'myllys19_anharmonic'):
    out4, keep = [], {}
    for r in csv.DictReader(open(os.path.join(here, 'acdc', 'acdc_%s_4x4.csv' % fe))):
        T, CS, CA, CB, J, A1 = (float(r[k]) for k in ('T', 'CS', 'CA', 'CB', 'J', '1A'))
        if r['converged'] != '1' or r['stable'] != '1' or not J > 1e-10:
            continue
        g = free_energies(fe, T)
        sat = klustr.P0/(klustr.KB*T)*1e-6*exp(g[0])
        x = sqrt(A1*CB)/sat                             # below 1: a barrier at the pair
        if x >= 1 or any(g[n-1] - (2*n - 1)*g[0] >= (2*n - 2)*log(x) for n in (2, 3, 4)):
            continue
        Jk, c = klustr.rate(CA, CB, CS, T, g[0], c4, clusters=True)
        # share of 2A2D from collisions of two pairs, 1/2 beta_11 c1^2 against k_1 c1
        bAB = c4['bABf']*sqrt(T)
        a1, b, s, y1 = (CA - c[0])/sat, CB/sat, CS/(bAB*sat), c[0]/sat
        k1 = c4['rhoA'][1]*a1*c4['rhoB'][1]*b/(c4['rhoA'][1]*a1 + c4['rhoB'][1]*b + s*c4['sig'][1])
        pp = 0.5*c4['rho'][1][1]*y1*y1
        c1 = sum(float(r['1A%dD' % n]) for n in range(1, 5))
        out4.append([r['T'], r['CS'], r['CA'], r['CB'], J, Jk, klustr.rate(CA, CB, CS, T, g[0], c4p),
                     c1, c[0], pp/(pp + k1*y1), x])
        keep[r['T'], r['CS'], r['CA'], r['CB']] = (J, Jk, g[0])
    out8 = []
    for r in csv.DictReader(open(os.path.join(here, 'acdc', 'acdc_%s_band8.csv' % fe))):
        key = r['T'], r['CS'], r['CA'], r['CB']
        if key not in keep or r['converged'] != '1' or not float(r['J']) > 1e-10:
            continue
        T, CS, CA, CB = (float(v) for v in key)
        J4, J4k, g1 = keep[key]
        sat, bAB = klustr.P0/(klustr.KB*T)*1e-6*exp(g1), c8['bABf']*sqrt(T)
        out8.append(list(key) + [float(r['J']), klustr.rate(CA, CB, CS, T, g1, c8), J4, J4k,
                                 bAB*sat*sat*iterated(CA/sat, CB/sat, CS/(bAB*sat), c8)])
    for name, head, out in (
            ('4x4', ['J_ACDC', 'J_KLUSTR', 'J_KLUSTR_alg1', 'c1_ACDC', 'c1_KLUSTR', 'share_pair_collisions',
                     'sqrt_CACB_csat'], out4),
            ('band8', ['J8_ACDC', 'J8_KLUSTR', 'J4_ACDC', 'J4_KLUSTR', 'J8_KLUSTR_iterated'], out8)):
        with open(os.path.join(here, 'acdc', 'klustr_%s_%s.csv' % (fe, name)), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['T', 'CS', 'CA', 'CB'] + head)
            w.writerows(out)
        print('klustr_%s_%s.csv: %d conditions' % (fe, name, len(out)))
