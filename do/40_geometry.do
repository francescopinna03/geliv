version 16.0
* Current design globals required. No partialling out inside GEL.
local zz $GEL_ZMAX
capture drop gel_probability gel_nweight gel_logratio gel_klcell
geliv $GEL_Y, endog($GEL_X) controls($GEL_W) instruments(`zz') method(et) weights(gel_probability)
generate double gel_nweight=e(N)*gel_probability if e(sample)
generate double gel_logratio=ln(gel_nweight) if gel_probability>0
generate double gel_klcell=gel_probability*gel_logratio if gel_probability>0
summarize gel_nweight, detail
scalar KL=e(Jdistance)/(2*e(N))
display "KL(p||uniform)=" scalar(KL)
histogram gel_nweight, fraction name(gel_weights, replace)
graph export output/et_weights.png, replace
* Observation ranking is descriptive: large weights are not causal influence.
preserve
generate double departure=abs(gel_logratio)
gsort -departure
list $GEL_Y $GEL_X gel_nweight gel_logratio in 1/20
keep in 1/100
export delimited using output/et_top_departures.csv, replace
restore
* Canonical homoskedastic inversion, separate from empirical HC0 GEL fits.
tempvar cons
generate byte `cons'=1
set seed 27102026
capture mata: mata drop gel_geometry_grid()
mata:
real matrix gel_geometry_grid(string scalar cons, string scalar zz)
{
    struct gel_rf scalar rf
    real matrix draws,components
    rf=gel_rf_init(st_data(.,st_global("GEL_Y")),st_data(.,st_global("GEL_X")),
        st_data(.,tokens(st_global("GEL_W")+" "+cons)),st_data(.,tokens(zz)))
    draws=gel_clr_draws(50000,rf.k)
    components=gel_ar_set(rf,.05)
    st_numscalar("AR_EMPTY",rows(components)==0)
    if(rows(components)) st_matrix("AR_COMPONENTS",components)
    return(gel_rf_grid(rf,rangen(-1,1,801),draws))
}
out=gel_geometry_grid("`cons'","`zz'")
end
if scalar(AR_EMPTY) display "AR confidence set is empty."
else {
    matrix colnames AR_COMPONENTS=lower upper left_infinite right_infinite
    matrix list AR_COMPONENTS
}
preserve
clear
getmata (beta AR_chi AR_F_p K K_p CLR CLR_p CLR_mcse TT)=out
export delimited using output/weak_grid_homoskedastic.csv, replace
twoway line AR_F_p K_p CLR_p beta, yline(.05) name(gel_tests, replace)
graph export output/weak_grid_homoskedastic.png, replace
restore
drop `cons'
display "K/CLR endpoints touching the finite grid are unresolved, not certified infinite."
