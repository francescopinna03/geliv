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

On 9 October 2026 the cell-level engine of version 0.2.0 compiled without errors
on the same installation. In the second block of `01_selftest.do`, on 1200
observations grouped in 60 cells in shuffled order, all 25 checks pass for ET,
EL and CUE. The dual value, gradient and Hessian of the two engines differ by at
most 2.2e-16, the envelope gradients by at most 2.4e-17, the coefficients by at
most 2.8e-17, and covariances, criteria and implied weights in the original order
by at most 5.6e-17. The row-level block, the CUE benchmark against `ivreg2` and
all the Card gates still pass. On Angrist and Krueger with 3 excluded
instruments, `cells(gel_cell)` reproduces the row-level CUE and ET estimates on
all 61 coefficients and their standard errors to the displayed precision, with
education coefficients of 0.108965 and 0.1089056, identical diagnostics for ET,
including J=3.0924534 and an effective sample size of 329505.9, and running times
of 52 seconds for CUE and 80 seconds for ET, against 362 seconds for the
row-level CUE.

## Angrist and Krueger estimates

`21_ak_estimate.do` ran on 9 October 2026 with the cell-level engine and
completed its 18 estimates in 1097 seconds on Stata/BE 19.0, with every
estimator converging. The table reports the coefficient on education, the
conventional standard error in parentheses, and for the GEL estimators the
overidentification statistic `Jgel` with its chi-squared p-value. The sample has
329509 observations in every cell of the table.

| Excluded instruments | OLS | 2SLS | LIML | CUE | ET | EL |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 0.0673 (0.0004) | 0.1077 (0.0196) | 0.1089 (0.0202) | 0.1090 (0.0196) | 0.1089 (0.0196) | 0.1088 (0.0196) |
| 30 | 0.0673 (0.0004) | 0.0869 (0.0159) | 0.0905 (0.0188) | 0.0925 (0.0159) | 0.0922 (0.0159) | 0.0919 (0.0159) |
| 180 | 0.0673 (0.0004) | 0.0928 (0.0097) | 0.1064 (0.0150) | 0.1068 (0.0096) | 0.1063 (0.0096) | 0.1056 (0.0096) |

| Excluded instruments | Degrees of freedom | J for CUE | J for ET | J for EL | Smallest effective sample size |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 3.10 (0.21) | 3.09 (0.21) | 3.09 (0.21) | 329505.9 |
| 30 | 29 | 23.80 (0.74) | 23.95 (0.73) | 24.05 (0.73) | 329485.1 |
| 180 | 179 | 161.26 (0.83) | 162.71 (0.80) | 162.53 (0.81) | 329341.5 |

The 2SLS coefficients reproduce the Python audit to every reported digit, at
0.1076937, 0.0869302 and 0.0928181, and the Cragg-Donald statistics reported by
`ivreg2`, 36.04, 5.36 and 2.58, coincide with the audit as well, which confirms
that the Stata construction of the 30 and 180 instruments matches the
independent one. The overidentification statistics reject none of the three
specifications at conventional levels, and the effective sample sizes never
fall below 329341, so the tilting required to satisfy the moments remains mild
throughout.

With 180 instruments the three GEL estimates, between 0.1056 and 0.1068, lie
close to LIML at 0.1064 and well above 2SLS at 0.0928. This pattern is
consistent with the many-instrument bias of 2SLS toward OLS and with the
higher-order results of Newey and Smith (2004), although on real data the
distance between estimators identifies no bias without knowledge of the true
coefficient, and with heterogeneous effects different instrument sets may weight
different subpopulations. The calibrated simulation is the step that measures
bias. Within the GEL family the ordering CUE, ET, EL holds in all three
specifications, with differences of at most 0.12 conventional standard errors.
With 180 instruments the conventional GEL standard errors, about 0.0096, are as
small as those of 2SLS and omit the variance inflation that many weak
instruments produce, so inference on these specifications should rest on the
robust tests of `22_ak_weakiv.do`.

## Pending

The full-sample comparison in `22_ak_weakiv.do`, the analysis of implied weights
in `40_geometry.do` and the calibrated simulation in `30_mc.do` remain to be
run. The gates already passed certify the engine on the cases described above,
and the remaining Angrist and Krueger results acquire the same status only once
each script has run its own checks.
