version 16.0
* Run 10_card.do OR 20_ak_prepare.do first, keeping its dataset in memory.
* Calibrated fixed-design Gaussian experiment. No empirical result is imposed.
* 20 reps is a PILOT. Set 1000+ for final results and report MC uncertainty.
local reps 20
local beta0 .10
local mu2grid 10 30 100
set seed 77102026
tempvar cons ysim xsim
generate double `cons'=1
generate double `ysim'=.
generate double `xsim'=.
capture mata: mata drop gel_mc_rf()
mata:
real rowvector gel_mc_rf(real colvector ys, real colvector xs, real matrix W,
    string scalar zz, real scalar beta0)
{
    struct gel_rf scalar rr
    rr=gel_rf_init(ys,xs,W,st_data(.,tokens(zz)))
    return(gel_rf_test(rr,beta0,gel_clr_draws(10000,rr.k)))
}
W=st_data(.,tokens(st_global("GEL_W")+" `cons'"))
Z=st_data(.,tokens(st_global("GEL_ZMAX")))
x0=st_data(.,st_global("GEL_X")); y0=st_data(.,st_global("GEL_Y"))
ZR=Z-W*qrsolve(W,Z)
baseX=W*qrsolve(W,x0)
signal=ZR*qrsolve(ZR,x0)
v0=x0-baseX-signal
baseY=W*qrsolve(W,y0-`beta0'*x0)
u0=y0-`beta0'*x0-baseY
u0=u0-ZR*qrsolve(ZR,u0)
Sigma=quadcross((u0,v0),(u0,v0))/(rows(W)-cols(W)-cols(ZR))
L=cholesky(Sigma)
assert(!hasmissing(L))
assert(quadcross(signal,signal)>0)
/* Freeze estimated direction; independently set total concentration mu^2.
   Thus adding instruments does not mechanically increase DGP signal. */
unitSignal=signal*sqrt(Sigma[2,2]/quadcross(signal,signal))
end
tempname P T
postfile `P' double mu2 int rep excluded str8 estimator double beta se byte ok cover double ess using output/mc_estimates.dta, replace
postfile `T' double mu2 int rep excluded double ar_p k_p clr_p clr_mcse using output/mc_tests.dta, replace
foreach strength of local mu2grid {
    forvalues rep=1/`reps' {
        mata: UV=rnormal(rows(W),2,0,1)*L'
        mata: xs=baseX+sqrt(`strength')*unitSignal+UV[.,2]
        mata: ys=`beta0'*xs+baseY+UV[.,1]
        mata: st_store(.,"`xsim'",xs); st_store(.,"`ysim'",ys)
        foreach k in $GEL_SPECS {
            local zz ${GEL_Z`k'}
            foreach method in ols tsls liml cue et el {
                local ess=.
                if "`method'"=="ols" capture quietly regress `ysim' `xsim' $GEL_W, vce(robust)
                else if "`method'"=="tsls" capture quietly ivregress 2sls `ysim' $GEL_W (`xsim'=`zz'), vce(robust)
                else if "`method'"=="liml" capture quietly ivregress liml `ysim' $GEL_W (`xsim'=`zz'), vce(robust)
                else capture quietly geliv `ysim', endog(`xsim') controls($GEL_W) instruments(`zz') method(`method')
                local rc=_rc
                if `rc'==0 {
                    if inlist("`method'","cue","et","el") local ess=e(ess)
                    local cover=cond(missing(_b[`xsim'],_se[`xsim']),.,abs(_b[`xsim']-`beta0')<=invnormal(.975)*_se[`xsim'])
                    post `P' (`strength') (`rep') (`k') ("`method'") (_b[`xsim']) (_se[`xsim']) (1) (`cover') (`ess')
                }
                else post `P' (`strength') (`rep') (`k') ("`method'") (.) (.) (0) (.) (.)
            }
            mata: st_matrix("MCT",gel_mc_rf(ys,xs,W,"`zz'",`beta0'))
            post `T' (`strength') (`rep') (`k') (MCT[1,2]) (MCT[1,4]) (MCT[1,6]) (MCT[1,7])
        }
        display "mu2=`strength', replication `rep'/`reps' completed"
    }
}
postclose `P'
postclose `T'
preserve
use output/mc_estimates.dta, clear
generate double error=beta-`beta0'
generate double squared_error=error^2
generate byte failure=1-ok
* Failure-aware coverage bounds: failed runs count first as all uncovered,
* then as all covered. Bias/RMSE are explicitly conditional on success.
generate byte cover_lower=cond(missing(cover),0,cover)
generate byte cover_upper=cond(missing(cover),1,cover)
collapse (mean) bias=error mse=squared_error coverage_conditional=cover failure_rate=failure coverage_lower=cover_lower coverage_upper=cover_upper (p50) median_error=error (sd) sd_error=error (count) successful=beta attempted=ok, by(mu2 excluded estimator)
generate double rmse=sqrt(mse)
generate double bias_mcse=sd_error/sqrt(successful)
generate double coverage_mcse=sqrt(coverage_conditional*(1-coverage_conditional)/successful)
export delimited using output/mc_summary.csv, replace
use output/mc_tests.dta, clear
generate byte cover_ar=ar_p>=.05 if ar_p<.
generate byte cover_k=k_p>=.05 if k_p<.
generate byte cover_clr=clr_p>=.05 if clr_p<.
collapse (mean) cover_ar cover_k cover_clr (count) reps=cover_ar, by(mu2 excluded)
export delimited using output/mc_weak_coverage.csv, replace
restore
drop `cons' `ysim' `xsim'
