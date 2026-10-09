# Validation record

Checks in Python on 7 and 9 October 2026, checks in Stata on 8 and 9 October 2026.

## Run in Python

`python/verify_math.py` passes 35 mathematical and numerical checks, plus 3009
pointwise comparisons between the analytic AR inversion and the F test. The
checks cover inner gradients and Hessians, envelope gradients, independent
primal problems solved with SLSQP, weighted moments, CUE Newton against the
closed form, invariance to invertible transformations of the instruments,
equality under exact identification, the ET Pythagorean identity and block
accumulation. The simulation also checks the size of K at zero concentration
and of the conditional CLR at four values of TT. Thresholds and errors are in
`output/math_validation.json`.

`python/audit_ak_design.py` was run on the original minimal MIT extract. The
checksum of the downloaded archive is

    0b6c505e52fd6c9f6194c356c91688acf2f8937b5023a84e754e942c9a2f2dcc

The sample has 329509 observations before and after the declared filters, it is
common to the three specifications, and the controls are 59 plus the intercept.

| Excluded instruments | Total moments | Instrument rank | 2SLS education coefficient | Homoskedastic first-stage F |
|---:|---:|---:|---:|---:|
| 3 | 63 | 63 | 0.107693705 | 36.036354 |
| 30 | 90 | 90 | 0.086930200 | 5.363572 |
| 180 | 240 | 240 | 0.092818056 | 2.582341 |

The OLS coefficient with the same controls is 0.067338971. These figures are
Python estimates from sparse cross products, produced before any Stata, Mata or
GEL run. The 2SLS trajectory is not monotone, which already rules out a
predetermined empirical narrative, and the distance from OLS measures no bias
without knowledge of the true structural coefficient.

`python/verify_cells.py`, run on 9 October 2026, implements the row-level and
cell-level computations independently on 1200 simulated observations grouped in
60 cells in shuffled order. Scaling vectors, 2SLS and OLS starting values, the
dual value, gradient and Hessian, the envelope gradient for ET, EL and CUE, the
covariance and the reordering of implied weights agree within 1.5e-14, against a
tolerance of 1e-10. The report is `output/cells_validation.json`.

## Run in Stata

On 8 October 2026, on Stata/BE 19.0 for macOS, the two Mata modules and the
`geliv` command of version 0.1.3 compiled without errors. The synthetic
self-test of `01_selftest.do`, on 800 observations with two excluded
instruments, passes every check, including the comparison between Newton and
the closed form of the CUE, the finite-difference check of the outer gradients
of ET, EL and CUE, invariance to the block size, coincidence with 2SLS under
exact identification and agreement between the analytic AR inversion and a grid
of 400 points. The same file compares the overidentified CUE with
`ivreg2, cue robust`, and the two coefficients differ by 1.1e-10.

On Card (1995), with 3010 observations, `10_card.do` reproduces the 2SLS
coefficient 0.131504 and verifies that ET, EL and CUE coincide with it under
exact identification, with uniform implied weights. With `nearc4` and `nearc2`
as excluded instruments, the CUE of `geliv` and that of `ivreg2` agree on all 17
coefficients within 1e-6, with an education coefficient of 0.1623756, while ET
and EL converge to 0.1622988 and 0.1622307 with moment errors below 1e-15 and
no negative weight.

`11_weakiv_benchmark.do` compares the canonical homoskedastic statistics at the
null value 0.1 with those of `weakiv`. The module computes AR and K with the
residual variance of the full reduced form, whereas `weakiv` uses the Lagrange
multiplier form built on the residual from the controls alone, and the two
versions are linked by the identity AR_LM = AR/(1 + AR/n), which holds in the
same way for K. After the conversion, AR, K and the CLR statistic agree with
`weakiv` to the seventh significant digit, with values 2.8329641, 1.4888246 and
1.6032545, and the CLR p-value differs by 0.002, within the Monte Carlo
tolerance of the conditional simulation.

On Angrist and Krueger (1991), `20_ak_prepare.do` rebuilds in Stata the same
sample as the Python audit, with 329509 observations and sets of 3, 30 and 180
excluded instruments. With 3 excluded instruments the row-level CUE of `geliv`
estimates the education coefficient at 0.108965 and ET at 0.1089056, both
converged, with a J statistic of 3.10 on 2 degrees of freedom and implied
weights between 0.96 and 1.07 times 1/n. The CUE estimate took 362 seconds on
Stata/BE.

## Pending

The cell-level engine of version 0.2.0 has not yet been compiled in Stata. Its
acceptance requires the second block of `01_selftest.do`, which compares the two
engines on shuffled cell data for ET, EL and CUE, the unchanged Card gates, and
the reproduction of the row-level estimates on Angrist and Krueger with 3
excluded instruments, 0.108965 for CUE and 0.1089056 for ET, through
`cells(gel_cell)`.

The estimates with 30 and 180 excluded instruments, the full-sample comparison
in `22_ak_weakiv.do`, the analysis of implied weights in `40_geometry.do` and
the calibrated simulation in `30_mc.do` remain to be run. The gates already
passed certify the engine on the cases described above and extend to the
Angrist and Krueger results only once each script has run its own checks.
