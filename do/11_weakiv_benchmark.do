version 16.0
* Designed for the current Card/AK design globals, one endogenous regressor.
* The custom canonical branch is HOMOSKEDASTIC only.
* This script records differences and refuses to certify a failed comparison.
set seed 7102026
local zz $GEL_ZMAX
tempvar one
generate byte `one'=1
capture mata: mata drop gel_bench_weakiv()
mata:
void gel_bench_weakiv(string scalar one, string scalar zz)
{
    struct gel_rf scalar rf
    real matrix draws
    rf=gel_rf_init(st_data(.,st_global("GEL_Y")),st_data(.,st_global("GEL_X")),
        st_data(.,tokens(st_global("GEL_W")+" "+one)),st_data(.,tokens(zz)))
    /* weakiv's large-sample iid normalization uses n rather than residual df.
       Set this explicitly for the benchmark, keeping the exact F branch separate. */
    rf.Sigma=rf.E/rf.n
    draws=gel_clr_draws(200000,rf.k)
    st_matrix("GELTEST",gel_rf_test(rf,.1,draws))
    st_numscalar("GELN",rf.n)
}
gel_bench_weakiv("`one'","`zz'")
end
quietly ivreg2 $GEL_Y $GEL_W ($GEL_X=`zz')
weakiv, null(.1)
scalar lm_factor=1+GELTEST[1,1]/scalar(GELN)
matrix BENCH=(GELTEST[1,1]/scalar(lm_factor),e(ar_chi2) \ GELTEST[1,3]/scalar(lm_factor),e(k_chi2) \ GELTEST[1,5],e(clr_stat) \ GELTEST[1,6],e(clr_p))
matrix rownames BENCH=AR_chi K CLR CLR_p
matrix colnames BENCH=custom weakiv
matrix list BENCH
* Failure means investigate normalization/version before interpreting estimates.
assert abs(BENCH[1,1]-BENCH[1,2])<1e-6
assert abs(BENCH[2,1]-BENCH[2,2])<1e-6
assert abs(BENCH[3,1]-BENCH[3,2])<1e-6
assert abs(BENCH[4,1]-BENCH[4,2])<max(.001,6*GELTEST[1,7])
display "Canonical weak-IV benchmark passed (CLR p-value compared at MC precision)."
