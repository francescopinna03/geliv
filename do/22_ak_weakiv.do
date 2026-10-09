version 16.0
* Heteroskedastic scalar confidence sets with controls, using external weakiv.
* Run 20_ak_prepare.do first. The scalar null does not fix nuisance coefficients.
set seed 7102026
capture log close weakak
log using output/ak_weakiv.log, text replace name(weakak)
foreach k in $GEL_SPECS {
    local zz ${GEL_Z`k'}
    ivreg2 $GEL_Y $GEL_W ($GEL_X=`zz'), robust
    weakiv, null(.1) usegrid gridmin(-.5) gridmax(.5) gridpoints(401) clrsims(50000)
    display "Excluded instruments: `k'"
    display "AR set: `e(ar_cset)'"
    display "K set: `e(k_cset)'"
    display "CLR set: `e(clr_cset)'"
    estimates save output/ak_weakiv_`k', replace
}
log close weakak
display "These are finite-grid sets. Extend/refine the grid where necessary."
