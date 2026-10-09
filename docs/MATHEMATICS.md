# Definitions and mathematical decisions

## 1. Criterion and signs

The vector of structural coefficients, of dimension p, is denoted theta. The
vector of included and excluded instruments for observation i, of dimension m,
is z_i, and the vector of regressors, of dimension p, is x_i. The moments are

    g_i(theta) = z_i (y_i - x_i' theta).

The multipliers, of dimension m, are lambda. With the concave functions
rho_ET(v)=-exp(v), rho_EL(v)=log(1-v) and rho_CUE(v)=-(1+v)^2/2, the normalized
mean criterion is

    q_gamma(theta) = sup_lambda mean[rho_gamma(lambda'g_i)-rho_gamma(0)].

The criterion q is minimized over theta, and EL requires 1-lambda'g_i>0 for every
observation. The limit gamma=-1 of the Cressie-Read parametrization is understood
net of the divergent constant that depends on neither parameters nor multipliers.

The software minimizes the negative of the concave inner objective. For ET it
uses the monotone logarithmic transformation of the exponential. At fixed theta,

    L(theta,lambda) = log(mean(exp(lambda'g_i(theta)))),
    q_ET_log(theta) = - min_lambda L(theta,lambda).

The inner gradient is the mean of the moments under the exponential weights and
the Hessian is their weighted covariance, which allows stable Newton steps with
the log-sum-exp device. The original GEL criterion equals 1-exp(-q_ET_log), a
strictly increasing transformation with the same minimizers, and both criteria
are stored.

For EL the convex inner objective is -mean(log(1-lambda'g_i)), and for CUE it is
mean(lambda'g_i + (lambda'g_i)^2/2). Both Hessians are positive semidefinite,
although convexity alone guarantees no finite optimum.

The sample mean of the moments is gbar, and the matrix of uncentered second
moments is Omega2=mean(g_i g_i'). When Omega2 is invertible,

    lambda_CUE = -solve(Omega2,gbar),
    q_CUE = (gbar' solve(Omega2,gbar))/2.

The centered covariance differs from Omega2 away from the null, so the two cannot
be exchanged on the strength of a claimed numerical equality between
implementations. The CUE estimator itself is the same under either choice,
because the centered criterion is a strictly increasing function of the
uncentered one.

## 2. Outer gradient

The derivative of the residual with respect to theta is -x_i. Define the scalar
t_i=z_i'lambda. The derivatives of the three criteria actually minimized are

    ET log     sum_i pi_i x_i t_i,
    EL         mean_i x_i t_i / (1-u_i t_i),
    CUE        mean_i (1+u_i t_i) x_i t_i.

The code optimizes over coefficients scaled by the root mean square of the
regressors and scales the instruments by their own root mean squares. These
transformations are invertible and fixed with respect to theta, so they leave
every estimator unchanged, and results are returned in the original units. No
Frisch-Waugh-Lovell partialling is applied to GEL.

Newton regularization, when a step requires it, affects only the Hessian used to
propose a direction. The criterion receives no penalty and the moment system is
left intact, and convergence is always assessed on the gradient of the
unpenalized problem. Singular systems in the closed-form CUE or in the
covariance are never resolved through a silent pseudoinverse. Newton for EL also
monitors the normalized moment, which prevents false convergence when all the
unnormalized weights tend to zero.

## 3. Probabilities and divergences

The uniform empirical distribution assigns 1/n to each observation. The implied
weights are

    ET      pi_i = exp(v_i) / sum_j exp(v_j),
    EL      pi_i = [1/(1-v_i)] / sum_j [1/(1-v_j)],
    CUE     pi_i = (1+v_i) / sum_j (1+v_j).

For ET and EL, under interior feasibility, these are probabilities. CUE imposes
no condition 1+v_i>=0 and may therefore produce signed weights, which then admit
neither an interpretation as a distribution nor an effective sample size.

ET minimizes KL(pi||uniform) and EL minimizes KL(uniform||pi), subject to weights
summing to one and to weighted moments equal to zero. The ET dual without an
explicit multiplier for the normalization yields the same solution through the
logarithmic transformation above. At its inner optimum

    KL(pi_ET||uniform) = -log(mean(exp(v_i))).

For EL, at the inner optimum, the unnormalized weights sum to n, and

    KL(uniform||pi_EL) = mean(log(1-v_i)).

The quadratic case of GEL is equivalent to the CUE with uncentered second
moments. A quadratic primal with an explicit positivity constraint differs from
it whenever the unconstrained CUE produces negative weights.

## 4. Pythagorean identity and its domain

Fix theta. The set of probability vectors that satisfy the moments at that theta
is an affine family intersected with the simplex. Let pi* be the ET projection of
the uniform distribution, assumed interior to the simplex. For every other
feasible distribution r,

    KL(r||uniform) = KL(r||pi*) + KL(pi*||uniform).

The identity follows from log(pi*_i/(1/n))=lambda'g_i-log(mean exp(lambda'g_i)),
whose first term has mean zero under r as under pi*, by subtracting the two
divergences. It compares no projections computed at different values of theta.
The subsequent minimization over theta moves across different affine families,
whose union need not be convex.

With nested moments and at fixed theta, adding groups of constraints yields an
informational decomposition. Once the estimates of theta differ, calling such a
difference an exact Pythagorean decomposition would be incorrect. Likewise,
extreme weights flag the observations that the reconciliation of the moments
relies on, and they say nothing by themselves about the validity of an
individual instrument.

## 5. Conventional covariance and overidentification

The mean Jacobian, an m by p matrix, is D=-mean(z_i x_i'). Under strong
identification and correct specification,

    V(theta_hat) = inverse(D' inverse(Omega2) D) / n.

This is the first-order efficient covariance shared by the three estimators. It
includes neither higher-order corrections nor robustness to misspecification,
and under weak instruments the usual Wald intervals may have incorrect coverage
even when the matrix is invertible.

At the estimated parameter and under regularity with fixed m, the GEL criterion
multiplied by 2n has a chi-squared(m-p) reference. For ET the package returns
both 2n*(1-exp(-KL)) and 2n*KL, which are asymptotically equivalent under the
null and numerically distinct in finite samples. A small distance statistic is
no proof of the economic validity of the instruments.

## 6. Structural null and weak instruments

For a null on the full vector theta0, a system of m restrictions on the mean of
the moments is evaluated. The reference for the AR/S-type statistic is
chi-squared(m), which requires no strong identification of theta, with m held
fixed in the approximation, a nondegenerate variance of the moments and a
central limit theorem. The GEL criterion at the null, without subtracting its
minimum over theta, can be used for this purpose. The profiled ratio that
subtracts the minimum does not become chi-squared(p) under weak instruments.

For the robust K statistic of the module, Omega denotes the centered covariance
of the moments at the null. For each coefficient j, the covariance between
column j of the individual Jacobian and the moments is an m by m matrix denoted
C_j. The corrected Jacobian has columns

    Dtilde_j = D_j - C_j inverse(Omega) gbar.

The two statistics are

    AR = n*gbar' inverse(Omega) gbar,
    K = n*gbar' inverse(Omega) Dtilde
          inverse(Dtilde' inverse(Omega) Dtilde)
          Dtilde' inverse(Omega) gbar.

The reference for K is chi-squared(p), under the appropriate rank conditions and
joint limit. The covariance correction distinguishes K from the square of the
plain GEL gradient weighted by the inverse observed Hessian, and it is essential.
The function rejects a singular projected matrix. Both are tests on the full
vector, and plugging in estimates of nuisance parameters does not automatically
preserve the justification just given.

## 7. Canonical AR, K and CLR with one endogenous regressor

For this branch only, the outcome, the endogenous regressor and the excluded
instruments are projected off the exogenous controls W. Let Y be the matrix of
the two residualized outcomes, in the order y and x, and Z the matrix of
residualized excluded instruments. The two 2 by 2 matrices of quadratic
products are

    R = Y' Z inverse(Z'Z) Z' Y,
    E = Y'Y - R.

The number of excluded instruments is k and the number of controls including
the intercept is l. The module uses Sigma=E/(n-l-k). At the null beta, the 2 by 1
vectors are b=(1,-beta)' and a=(beta,1)', and the scalar sufficient statistics
are

    SS = b' R b / (b' Sigma b),
    TT = a' inverse(Sigma) R inverse(Sigma) a / (a' inverse(Sigma) a),
    ST = b' R inverse(Sigma) a / sqrt((b' Sigma b)*(a' inverse(Sigma) a)).

K equals ST²/TT, and the conditional likelihood ratio statistic is

    CLR = (SS-TT + sqrt((SS-TT)²+4*ST²))/2.

Its conditional distribution is simulated holding TT fixed. Two independent
variables U1, distributed chi-squared(1), and U2, distributed chi-squared(k-1),
are drawn and substituted as SS=U1+U2 and ST²=TT*U1. The p-value is reported
together with its Monte Carlo error, and a single matrix of draws is reused across
the whole grid. With known covariance and normality the canonical construction
is exact, and once the covariance is replaced by its reduced-form estimate the
implemented CLR is asymptotic.

The `weakiv` command reports AR and K in Lagrange multiplier form, with the
residual variance taken from the projection on the controls alone. The two forms
are linked by AR_LM = AR/(1 + AR/n), and the same relation holds for K, which is
how the benchmark in `11_weakiv_benchmark.do` converts the statistics of the
module before comparing them.

The classical AR test uses SS/k with the F(k,n-l-k) reference, exact under
homoskedastic Gaussian errors. Its inversion is equivalent to

    (1,-beta) [R - k*Fcrit*E/(n-l-k)] (1,-beta)' <= 0.

The code solves this quadratic inequality and distinguishes the empty set, a
bounded interval, a half-line, the union of two half-lines and the whole line,
so that infinite tails are certified without mistaking the edge of a grid for
infinity. K and CLR remain finite grids, and their topology between nodes or
outside the domain is not certified. An empty AR set signals that the moments
may be incompatible, and it should never be read as extremely precise
identification.

## 8. Scope of the Newey and Smith results

Newey and Smith (2004) show that the GEL family removes a component of the
asymptotic bias due to the correlation between moments and Jacobian, and that EL
also removes a component linked to the estimation of the weighting matrix. Their
results establish neither that every bias of ET, EL and CUE is independent of the
number of moments nor any uniform guarantee under weak-instrument or
many-instrument asymptotics. The project therefore measures performance and
failures directly and makes no prior assumption about the ranking of the
estimators.

The simulation is an experiment on finite-sample performance inspired by these
expansions, and it tests no theorem outside its assumptions. Fixed m, growing m,
constant overall strength and constant strength per instrument are distinct
regimes. The program fixes the overall concentration in the largest model, so
the concentration projected on the subsets of instruments may be lower.

## 9. Aggregation by cells

Suppose the observations are partitioned into C cells and that the instruments
and the exogenous regressors take a single value within each cell, while the
outcome and the endogenous regressors vary across individuals. Denote by c(i)
the cell of observation i, by n_c the number of observations in cell c, by z_c
the vector of instruments of cell c, of dimension m, and by w_c the vector of
exogenous regressors of cell c, of dimension p-e, where e is the number of
endogenous regressors. The regressor vector splits as x_i=(x_i^E, w_c(i)), with
x_i^E of dimension e. The residual and the scalar argument of rho become

    u_i = y_i - x_i^E' theta_E - w_c(i)' theta_W,
    t_c = z_c' lambda,
    v_i = u_i t_c(i).

Every sum over observations that multiplies a cell-level matrix can be computed
by first summing scalar or e-dimensional quantities within cells. With
psi_i the derivative weight of the method at observation i, namely exp(v_i-shift)
for ET, 1/(1-v_i) for EL and 1+v_i for CUE, and with omega_i the corresponding
Hessian weight, exp(v_i-shift) for ET, 1/(1-v_i)^2 for EL and 1 for CUE, the
inner gradient and Hessian are

    sum_i psi_i u_i z_c(i)          = sum_c z_c a_c,         a_c = sum_{i in c} psi_i u_i,
    sum_i omega_i u_i^2 z_c(i) z_c(i)' = sum_c z_c z_c' h_c,  h_c = sum_{i in c} omega_i u_i^2,

up to the normalizations of Section 1. Denote by kappa_i the scalar factor that
multiplies x_i t_i in the expression of Section 2 for the method in use. The
outer gradient then splits into a part for the endogenous regressors,
sum_i kappa_i t_c(i) x_i^E, which remains an individual sum of dimension e, and a
part for the exogenous regressors, sum_c w_c t_c (sum_{i in c} kappa_i).
The Jacobian, the matrix Omega2, the 2SLS and OLS starting values and the
weighted moments decompose in the same way.

The cost of one evaluation thus falls from the order of n m² to the order of
n + C m², and no n by m matrix is ever formed. Every quantity is algebraically
identical to its row-level counterpart, so estimators, covariances and implied
weights are unchanged and only the order of summation differs. The engine sorts
the observations by cell in Mata through `order()`, `panelsetup()` and
`panelsum()`, checks the constancy of instruments and exogenous regressors
within each cell exactly, and returns the implied weights in the original order
of the observations. In the Angrist and Krueger design the cells are defined by
year, quarter and state of birth, with at most 2040 of them for 329509
observations.

## Primary sources

[1] Newey, W. K. and Smith, R. J. (2004), Higher order properties of GMM and
generalized empirical likelihood estimators, *Econometrica* 72, 219-255.
https://doi.org/10.1111/j.1468-0262.2004.00482.x

[2] Kleibergen, F. (2005), Testing parameters in GMM without assuming that they
are identified, *Econometrica* 73, 1103-1123.
https://doi.org/10.1111/j.1468-0262.2005.00610.x

[3] Moreira, M. J. (2003), A conditional likelihood ratio test for structural
models, *Econometrica* 71, 1027-1048. https://doi.org/10.1111/1468-0262.00438

[4] Kitamura, Y. and Stutzer, M. (1997), An information-theoretic alternative to
generalized method of moments estimation, *Econometrica* 65, 861-874.

[5] Finlay, K., Magnusson, L. M. and Schaffer, M. E., code and documentation of
`weakiv`. https://github.com/kfinlay/weakiv

[6] StataCorp, Mata reference manual, `optimize()` and `panelsetup()`.
https://www.stata.com/manuals/m-5optimize.pdf

[7] Angrist, J. D., data archive. https://economics.mit.edu/people/faculty/josh-angrist/angrist-data-archive
