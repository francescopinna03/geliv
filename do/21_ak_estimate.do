version 16.0
* Run 20_ak_prepare.do in this session first.
confirm variable $GEL_Y $GEL_X
capture log close ak
log using output/ak.log, text replace name(ak)
tempname P
postfile `P' int excluded str8 estimator double beta se N ess J int rc using output/ak_results.dta, replace
foreach k in $GEL_SPECS {
    local zz ${GEL_Z`k'}
    foreach method in ols tsls liml cue et el {
        local rc=0
        local ess=.
        local J=.
        if "`method'"=="ols" capture noisily regress $GEL_Y $GEL_X $GEL_W, vce(robust)
        else if "`method'"=="tsls" capture noisily ivreg2 $GEL_Y $GEL_W ($GEL_X=`zz'), robust first
        else if "`method'"=="liml" capture noisily ivreg2 $GEL_Y $GEL_W ($GEL_X=`zz'), liml robust
        else {
            capture noisily geliv $GEL_Y, endog($GEL_X) controls($GEL_W) instruments(`zz') method(`method') weights(p`k'_`method') cells($GEL_CELL)
        }
        local rc=_rc
        if `rc'==0 {
            if inlist("`method'","cue","et","el") {
                local ess=e(ess)
                local J=e(Jgel)
            }
            post `P' (`k') ("`method'") (_b[$GEL_X]) (_se[$GEL_X]) (e(N)) (`ess') (`J') (0)
        }
        else post `P' (`k') ("`method'") (.) (.) (_N) (.) (.) (`rc')
    }
}
postclose `P'
save output/ak_weights.dta, replace
preserve
use output/ak_results.dta, clear
export delimited using output/ak_results.csv, replace
list, separator(6)
restore
log close ak
