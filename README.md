# KLUSTR

KLUSTR calculates the formation rate of new particles from an acid and a base when their pair is the critical
cluster or there is no barrier. The only thermodynamic input is the formation free energy of the pair.

The model is described in

> Halonen, R.: KLUSTR v1.0: a closed-form steady-state model of atmospheric new particle formation from the free
> energy of the acid–base pair, Geosci. Model Dev., submitted, 2026.

## Contents

```
python/        KLUSTR in Python and a test
fortran/       KLUSTR as a Fortran module and a test
verification/  the ACDC simulations and KLUSTR rates of Sect. 5 of the paper
```

## Usage

Python 3.6 or newer, no other packages:

```python
import klustr
const = klustr.setup(N=8, gamma=-1.7, mA=98.08, mB=45.08, densA=1830.0, densB=680.0, eta=3.0)
J = klustr.rate(a_mono, b_free, CS_A, T, dG_kT, const)
```

Fortran 2008:

```fortran
use klustr
call klustr_setup(8, -1.7_wp, 98.08_wp, 45.08_wp, 1830.0_wp, 680.0_wp, 3.0_wp)
J = klustr_rate(a_mono, b_free, cs_a, t, dg_kt)
```

The inputs are the monomer concentration of the acid and the free concentration of the base (cm-3), the loss rate
coefficient of the acid (s-1), the temperature (K) and the standard formation free energy of the pair at 1 atm in
units of kB T. J is the formation rate past N pairs (cm-3 s-1). N = 8 is about 1.7 nm in mobility diameter for
sulfuric acid and dimethylamine (Kürten et al., 2018, https://doi.org/10.5194/acp-18-845-2018). The arguments are
explained in `klustr.py` and `klustr.F90`.

Tests: `python test_klustr.py` in `python/` and `make test` in `fortran/`.

## Verification

In MATLAB, run `acdc_verification.m` in `verification/acdc/` (ACDC 2024_02_12 is included). Then
`verification/klustr_rates.sh` calculates the KLUSTR rates for the ACDC output.

## Licence

MIT. `verification/acdc/ACDC/` is an unchanged copy of ACDC 2024_02_12 (https://github.com/tolenius/ACDC, commit
870b82a), under its own GPL-3.0 licence.
