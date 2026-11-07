"""KLUSTR v1.0: formation rate of new particles from an acid and a base whose pair is the critical cluster or has no
barrier (Algorithm 1 of Halonen, R., Geosci. Model Dev., submitted, 2026). Standard library only.

    const = setup(N=8, gamma=-1.7, mA=98.08, mB=45.08, densA=1830.0, densB=680.0, eta=3.0)    # once
    J = rate(a_mono, b_free, CS_A, T, dG_kT, const)                                          # cm-3 s-1

Copyright (c) 2026 Roope Halonen. MIT licence.
"""
from math import sqrt, exp, pi

__version__ = '1.0'

KB, P0, AMU, NA = 1.380649e-23, 101325.0, 1.66053907e-27, 6.02214076e23


def setup(N, gamma, mA, mB, densA, densB, eta, closed=True, first_collision=False):
    """Constants of one chemistry. N: pairs past which a cluster is formed (8: about 1.7 nm for sulfuric acid and
    dimethylamine); gamma: exponent of the loss rate coefficient in the diameter; mA, mB: masses (amu); densA, densB:
    bulk densities (kg m-3); eta: factor on the collision coefficients (1: hard spheres); closed: pair at most at its
    saturation concentration; first_collision: a cluster counted at its first collision past N pairs, as in ACDC."""
    if N < 1:
        raise ValueError('N must be at least 1')
    if min(mA, mB, densA, densB) <= 0:
        raise ValueError('masses and densities must be positive')
    vA, vB = mA*1e-3/(densA*NA), mB*1e-3/(densB*NA)                    # m3
    rad = lambda v: (3*v/(4*pi))**(1/3)
    sp = {'A': (mA*AMU, rad(vA)), 'B': (mB*AMU, rad(vB))}
    sp.update({n: (n*(mA + mB)*AMU, rad(n*(vA + vB))) for n in range(1, N + 1)})
    def g(p, q):                                                        # hard-sphere beta/sqrt(T) (m3 s-1 K-1/2)
        (m1, r1), (m2, r2) = sp[p], sp[q]
        return sqrt(8*KB*(m1 + m2)/(pi*m1*m2))*pi*(r1 + r2)**2
    gAB = g('A', 'B')
    rhoA = [0.0] + [g('A', n)/gAB for n in range(1, N + 1)]
    rhoB = [0.0] + [g('B', n)/gAB for n in range(1, N + 1)]
    rho = [[0.0]*(N + 1)] + [[0.0] + [g(i, j)/gAB for j in range(1, N + 1)] for i in range(1, N + 1)]
    sig = [0.0] + [(n*(vA + vB)/vA)**(gamma/3) for n in range(1, N + 1)]
    return {'rhoA': rhoA, 'rhoB': rhoB, 'rho': rho, 'sig': sig, 'bABf': eta*1e6*gAB, 'closed': closed,
            'first': first_collision}


def reduced(a, b, s, const, clusters=False):
    """Reduced formation rate F = J/(beta_AB sat^2) past N pairs for a = a_mono/sat, b = b_free/sat and
    s = CS_A/(beta_AB sat); with clusters=True also the reduced concentrations y_n = c_n/sat, n = 1 ... N."""
    rA, rB, rho, sig = const['rhoA'], const['rhoB'], const['rho'], const['sig']
    N = len(sig) - 1
    if a <= 0 or b <= 0:
        return (0.0, [0.0]*N) if clusters else 0.0
    kap = lambda a1, n: rA[n]*a1*rB[n]*b/(rA[n]*a1 + rB[n]*b + s*sig[n])
    a1 = a                                                   # free A; the pair balance with a1 = a - y1 is
    for _ in range(2):                                       # a quadratic in y1, kappa_1 taken from the previous pass
        M = 1 + kap(a1, 1) + s*sig[1] + b
        y1 = 2*a*b/(M + sqrt(M*M + 4*rho[1][1]*a*b))
        if const['closed']:
            y1 = min(y1, 1.0)                                # closed system: the pair at most saturated
        a1 = a - y1
    k = [0.0] + [kap(a1, n) for n in range(1, N + 1)]
    if const['first']:                                       # counted at the first collision with a molecule
        k[N] = rA[N]*a1 + rB[N]*b
    y = [0.0, y1]
    for n in range(2, N + 1):                                # larger clusters in order of size
        S = k[n-1]*y[n-1] + 0.5*sum(rho[i][n-i]*y[i]*y[n-i] for i in range(1, n))
        D = k[n] + sum(rho[i][n]*y[i] for i in range(1, n)) + s*sig[n]
        y.append(2*S/(D + sqrt(D*D + 4*rho[n][n]*S)))
    F = k[N]*y[N] + sum((0.5 if i == j else 1.0)*rho[i][j]*y[i]*y[j]
                        for i in range(1, N + 1) for j in range(i, N + 1) if i + j > N)
    return (F, y[1:]) if clusters else F


def rate(a_mono, b_free, CS_A, T, dG_kT, const, clusters=False):
    """Formation rate J (cm-3 s-1) past N pairs. a_mono: acid monomer concentration (all clusters with one acid),
    b_free: free base concentration (cm-3); CS_A: loss rate coefficient of the acid (s-1); T (K); dG_kT: standard
    formation free energy of the pair at 1 atm over kB T. clusters=True also returns c_1 ... c_N (cm-3)."""
    sat = P0/(KB*T)*1e-6*exp(dG_kT)                          # saturation concentration of the pair (cm-3)
    bAB = const['bABf']*sqrt(T)                              # collision coefficient of A and B (cm3 s-1)
    if clusters:
        F, y = reduced(a_mono/sat, b_free/sat, CS_A/(bAB*sat), const, clusters=True)
        return bAB*sat*sat*F, [yn*sat for yn in y]
    return bAB*sat*sat*reduced(a_mono/sat, b_free/sat, CS_A/(bAB*sat), const)


if __name__ == '__main__':
    # One documented condition, the first check to run: sulfuric acid and dimethylamine past 8 pairs, eta = 3,
    # 278 K, acid monomer 1e6 cm-3, free base 1e7 cm-3, CS_A = 1e-3 s-1, pair free energy of Olenius et al. (2013)
    const = setup(N=8, gamma=-1.7, mA=98.08, mB=45.08, densA=1830.0, densB=680.0, eta=3.0)
    T = 278.0
    dG_kT = (-24.6465441592 - T*(-31.008553e-3))/(1.987204259e-3*T)
    J = rate(1e6, 1e7, 1e-3, T, dG_kT, const)
    print(f'KLUSTR {__version__}: J = {J:.6e} cm-3 s-1 (expected 2.235716e+02)')
