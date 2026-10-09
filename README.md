# GELIV, a Mata laboratory for weak instruments

Version 0.2.1, 9 October 2026. Research code for Stata 16 or later.

GELIV implements exponential tilting (ET), empirical likelihood (EL) and the
continuously updated estimator (CUE) from first principles in Mata for the
linear instrumental variables model. The inner dual problem is solved by damped
Newton with analytic gradient and Hessian, and the outer problem by BFGS through
`optimize()`, with the envelope gradient and several starting points. Controls
and the intercept remain in the model and in the moment conditions, and the
moments are accumulated in blocks without ever storing the full n by m matrix.

When instruments and exogenous regressors are constant within the cells of a
discrete grouping, as in designs built on year, quarter and state of birth, the
option `cells()` switches to an engine that aggregates individual weights within
each cell before forming outer products of the moments. The cost of one
evaluation then falls from the order of n m² to the order of n + C m², where C
is the number of cells, and every estimator, covariance and implied weight is
the same as in the row-level engine.

The formulas and algorithms were checked numerically in Python, with the report
in `output/math_validation.json`. Since 8 October 2026 the row-level engine and
the `geliv` command have been compiled and run in Stata/BE 19.0, where they pass
the synthetic self-test, the comparison of the overidentified CUE with `ivreg2`
and the acceptance gates on Card, all documented in `docs/VALIDATION.md`. The
cell-level engine introduced in version 0.2.0 was verified against the row-level
formulas in Python, with the report in `output/cells_validation.json`, and on 9
October 2026 it compiled in Stata, matched the row-level engine to machine
precision in the self-test and reproduced the row-level Angrist and Krueger
estimates with 3 excluded instruments exactly. Full estimates with 30 and 180
instruments and the calibrated simulation remain to be run, and the package
contains no empirical table that anticipates them.

## Getting started

Extract the package and set its folder as the current directory in Stata. The
benchmarks rely on external packages, which `install_deps.do` installs once.

```stata
cd "/path/to/geliv_toolkit"
do install_deps.do
do run_selftest.do
do do/run_card_checks.do
```

`run_selftest.do` loads the Mata modules through `do/00_setup.do`, runs
`do/01_selftest.do` and writes `run_selftest.log`. The Card runner downloads the
data through `bcuse` and writes `output/acceptance.log`.

For Angrist and Krueger, the official extract is downloaded from the terminal
and the pipeline then runs in Stata. The downloader compares the size and the
SHA-256 hash of the archive with the tracked record in
`data/asciiqob.zip.source.json` and stops on any difference.

```bash
python3 python/download_data.py ak-minimal
```

```stata
do do/00_setup.do
do do/20_ak_prepare.do
do do/21_ak_estimate.do
do do/22_ak_weakiv.do
do do/40_geometry.do
do do/30_mc.do
```

`30_mc.do` holds a pilot configuration with 20 replications, which should be
raised before a final simulation, and the large simulation should only follow
the synthetic and Card checks. The Mata engine performs no network operation.
On Stata/BE 19.0 a single row-level CUE estimate with 3 excluded instruments
took 362 seconds on the full sample, against 52 seconds for the same estimate
with `cells()` and 80 seconds for ET. Since the row-level cost grows with the
square of the number of moments, the pipeline passes the cell index built by
`20_ak_prepare.do` to every call of `geliv`.

## Interface

```stata
geliv lwage, endog(educ) controls(exper expersq black smsa south smsa66 ///
    reg662 reg663 reg664 reg665 reg666 reg667 reg668 reg669) ///
    instruments(nearc4 nearc2) method(et) weights(p_et)
matrix list e(diagnostics)
matrix list e(starts)
```

`method()` accepts `et`, `el` and `cue`, and `block(4096)` sets the size of the
temporary blocks in the row-level engine. `noconstant` removes the intercept.
`from(M)` supplies a matrix of starting values, one per row, ordered as
endogenous regressors, controls and intercept. Coefficients are reported in the
original units, and the GEL engine accepts several endogenous regressors.
Variables must be expanded beforehand, without factor-variable syntax, and
collinear instruments are never dropped silently, so the rank of the instrument
matrix is the responsibility of the user.

`cells(varname)` names a numeric variable that indexes the cells within which
instruments and controls are constant. The command checks this constancy
exactly, observation by observation, and stops with an error when it fails. The
endogenous regressors and the outcome may vary freely within cells, at least one
control or the intercept must be present, and implied weights are returned in
the original order of the observations. `e(cells)` stores the name of the cell
variable.

The Mata interface uses `gel_data_init(y,X,Z,gamma,block)` for the row-level
engine, `gel_data_cells(y,E,W,Z,cell,gamma,block)` for the cell-level engine and
`gel_fit(d,starts)` for both. `X` contains the endogenous regressors, the
controls and the constant, while `Z` contains the excluded instruments, the
controls and the constant. In the cell-level call, `E` holds the endogenous
regressors and `W` the controls and the constant. The three values of `gamma`
are 0 for ET, -1 for EL and 1 for CUE. The engine is specific to the linear IV
moment function.

The Mata result contains coefficients, multipliers in the original units,
implied weights, the conventional covariance, numerical residuals and the log of
starting points. `e(diagnostics)` holds the following quantities in this order.

| Field | Meaning |
|---|---|
| converged | outer convergence and numerical certification of the dual |
| q | mean criterion, equal to the KL divergence for ET |
| Jgel | normalized GEL criterion multiplied by 2n |
| Jdistance | 2n times the criterion used by the optimizer |
| inner_error | maximum norm of the inner gradient, also normalized for EL |
| outer_error | maximum norm of the outer gradient in scaled coordinates |
| moment_error | maximum norm of the weighted moments, scaled instruments |
| ess | 1/sum(weights²), for nonnegative weights only |
| min_nw, max_nw | extremes of the ratio between weight and 1/n |
| negative | number of negative weights, possible for CUE |
| N, m, p | observations, total moments, total coefficients |

`e(starts)` records the return code, the convergence flag, the criterion and
the gradient norm for each starting point. The selected minimum is the best
among the solutions found, which carries no certificate of global optimality.
`e(V)` is the first-order HC0 covariance under strong identification and valid
moments, and Wald intervals built from it remain invalid under weak
identification.

## What is tested

`01_selftest.do` checks the CUE Newton iterations against the closed form, the
outer gradients against finite differences, invariance to the block size, the
weighted moments, exact identification and the AR inversion. When `ivreg2` is
installed it compares the overidentified CUE with `ivreg2, cue robust`, which is
the relevant benchmark because the homoskedastic CUE minimizes a different
criterion, with a tolerance of 1e-6 on the coefficient of interest. A second
block of the same file builds data with 60 cells in shuffled order and requires
the row-level and cell-level engines to agree on the dual value, gradient and
Hessian, on the envelope gradient, on the 2SLS and OLS starting values and, after
estimation, on coefficients, covariance, criterion and implied weights in the
original order, for ET, EL and CUE alike.

`10_card.do` checks N=3010 and the 2SLS coefficient 0.131504 with the declared
controls, and requires ET, EL and CUE to coincide with it in the exactly
identified case. The second instrument `nearc2` creates an overidentified
numerical test, which by itself offers no new economic defence of its validity.

`11_weakiv_benchmark.do` compares the canonical homoskedastic statistics with
`weakiv`, after converting them to the Lagrange multiplier form used by that
command, AR/(1 + AR/n) and likewise for K. The CLR p-value relies on conditional
simulation and therefore carries a Monte Carlo tolerance. A failed assertion
keeps the gate open until the discrepancy is understood.

## Inference included and its boundaries

| Function | Object | Assumptions and limits |
|---|---|---|
| `gel_full_tests(d,theta0)` | GMM-type AR/S and K on the full structural vector | heteroskedasticity, independent observations, fixed m, regular rank of the corrected Jacobian, row-level data only |
| `gel_null_grid(d,grid)` | GEL criterion at the null, without subtracting its minimum | p=1, no nuisance parameters, chi-squared with m degrees of freedom |
| `gel_rf_test(r,beta,draws)` | canonical AR, K and CLR for one endogenous regressor | homoskedastic branch, estimated reduced-form covariance, asymptotic CLR |
| `gel_ar_set(r,alpha)` | AR F inversion through a quadratic inequality | infinite tails determined analytically, homoskedastic Gaussian errors |
| `gel_rf_grid(r,grid,draws)` | numerical inversion of K, CLR and AR | grid edges certify no infinite tail |
| `weakiv` in the scripts | benchmark and robust scalar inference with controls | external implementation, options to be kept consistent |

The CUE-based CLR is never transferred to ET or EL by substituting their
criteria. A self-contained heteroskedastic or cluster-robust CLR is absent from
the Mata code, and so is profiled GEL inversion over the many coefficients of
the controls. For Card and Angrist and Krueger, robust scalar inference with
controls is delegated to `weakiv, robust` after a robust IV estimate, and the
canonical branch serves the homoskedastic comparison and the calibrated Gaussian
simulation.

A failure of the dual is an unresolved numerical outcome and carries no
automatic statistical meaning. The code requires accurate roots and contains no
linear-programming certificate of membership in the interior of the convex hull.

## Empirical choices

The main comparison uses the same sample and the same controls, including year
and state of birth, with 3, 30 and 180 excluded instruments. These are the three
quarter dummies, their 30 interactions with year, and the latter together with
150 interactions between quarter and state. The fourth quarter and Michigan are
the base categories. With 59 controls and the intercept, the total number of
moments is 63, 90 and 240.

The design is a controlled extension of the classical one and differs both from
a literal replication of Table V and from a numerical reproduction of the
results of Bound, Jaeger and Baker. Adding instruments may also change the
weighting of the return to education when effects are heterogeneous, so the
distance between estimates on real data identifies no bias by itself.

The simulation keeps the instrument design fixed, calibrates a first-stage
direction and the error covariance on the data, and then imposes an explicit
true coefficient and a total concentration of 10, 30 or 100. The same replication
is used for every specification. It reports bias, RMSE, median error, Wald
coverage, AR, K and CLR coverage and failures. Bias and RMSE are conditional on
success, while coverage also reports lower and upper bounds that account for
failures. The coverage of the robust tests is a property of the test procedure
for the model as a whole.

`MANIFEST.sha256` lists the SHA-256 hash of every tracked file and is rebuilt
with `python3 python/write_manifest.py` after any change, so that
`sha256sum -c MANIFEST.sha256` verifies a checkout. The only tracked files that
scripts rewrite are the manifest itself and the validation reports in `output/`,
which record the run that produced them.

`docs/MATHEMATICS.md` contains formulas, normalizations and interpretation,
`docs/DATA.md` the sources and the full extract, and `docs/VALIDATION.md` the
record of the checks already run and of those still pending, together with the
audit on the MIT extract.
