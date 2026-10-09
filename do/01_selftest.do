version 16.0
clear
set obs 800
set seed 7102026
generate double z1=rnormal()
generate double z2=rnormal()
generate double v=rnormal()
generate double x=z1+.3*z2+v
generate double y=.7*x+.2+.5*v+rnormal()
generate double one=1
capture mata: mata drop gel_check()
capture mata: mata drop gel_show_starts()
capture mata: mata drop gel_selftest()
mata:
real scalar gel_check(string scalar label, real scalar pass, real scalar value)
{
    printf("%s %-52s %14.6e\n", pass ? "[ OK ]" : "[FAIL]", label, value)
    return(!pass)
}
void gel_show_starts(struct gel_result scalar r)
{
    real scalar h
    for(h=1;h<=rows(r.starts);h++) {
        printf("    start %g: rc %g  converged %g  q %g  outer_error %g\n",
            h, r.starts[h,1], r.starts[h,2], r.starts[h,3], r.starts[h,4])
    }
    printf("    b %g %g  inner_error %g  moment_error %g  ess %g  negative %g\n",
        r.b[1], r.b[2], r.inner_error, r.moment_error, r.ess, r.negative)
}
void gel_selftest()
{
    struct gel_data scalar d,dj
    struct gel_dual scalar s,t
    struct gel_result scalar r,r2
    struct gel_rf scalar rf
    fails=0
    y=st_data(.,"y"); X=st_data(.,("x","one")); Z=st_data(.,("z1","z2","one"))
    d=gel_data_init(y,X,Z,1,73)
    b=(.71,.21):*d.sx
    s=gel_inner(d,b,0); t=gel_inner(d,b,1)
    fails=fails+gel_check("cue closed form ok (value: inner error)", s.ok, s.error)
    fails=fails+gel_check(sprintf("cue newton ok, %g iterations", t.iter), t.ok, t.error)
    fails=fails+gel_check("cue closed form vs newton lambda", max(abs(s.lambda-t.lambda))<1e-8, max(abs(s.lambda-t.lambda)))
    for(method=-1;method<=1;method++) {
        d.method=method
        s=gel_inner(d,b,0)
        fails=fails+gel_check(sprintf("method %g inner ok, %g iterations, f %g", method, s.iter, s.f), s.ok, s.error)
        gel_eval(1,b,d,q,g,H)
        fails=fails+gel_check(sprintf("method %g eval q nonmissing", method), !hasmissing(q), q)
        for(j=1;j<=d.p;j++) {
            bp=b; bm=b; bp[j]=bp[j]+1e-5; bm[j]=bm[j]-1e-5
            gel_eval(0,bp,d,qp,tmp,H); gel_eval(0,bm,d,qm,tmp,H)
            fd=(qp-qm)/2e-5
            fails=fails+gel_check(sprintf("method %g gradient %g: fd %g analytic %g", method, j, fd, g[j]), abs(fd-g[j])<1e-6, abs(fd-g[j]))
        }
        starts=gel_2sls(d)\qrsolve(X,y)'
        printf("method %g starts: 2sls %g %g, ols %g %g\n", method, starts[1,1], starts[1,2], starts[2,1], starts[2,2])
        r=gel_fit(d,starts)
        gel_show_starts(r)
        fails=fails+gel_check(sprintf("method %g fit converged (value: q)", method), r.converged, r.q)
        fails=fails+gel_check(sprintf("method %g moment error", method), r.moment_error<1e-7, r.moment_error)
        fails=fails+gel_check(sprintf("method %g weights sum to one", method), abs(sum(r.weights)-1)<1e-10, abs(sum(r.weights)-1))
        d.block=800; r2=gel_fit(d,starts); d.block=73
        gel_show_starts(r2)
        fails=fails+gel_check(sprintf("method %g block 800 converged (value: q)", method), r2.converged, r2.q)
        fails=fails+gel_check(sprintf("method %g block invariance", method), max(abs(r.b-r2.b))<1e-6, max(abs(r.b-r2.b)))
    }
    dj=gel_data_init(y,X,Z[.,(1,3)],0,73)
    b0=gel_2sls(dj)
    printf("just identified 2sls: %g %g\n", b0[1], b0[2])
    for(method=-1;method<=1;method++) {
        dj.method=method; r=gel_fit(dj,b0)
        gel_show_starts(r)
        fails=fails+gel_check(sprintf("just-id method %g converged (value: q)", method), r.converged, r.q)
        fails=fails+gel_check(sprintf("just-id method %g equals 2sls", method), max(abs(r.b-b0))<1e-6, max(abs(r.b-b0)))
        fails=fails+gel_check(sprintf("just-id method %g uniform weights", method), max(abs(r.weights:-1/rows(y)))<1e-8, max(abs(r.weights:-1/rows(y))))
    }
    rf=gel_rf_init(y,X[.,1],X[.,2],Z[.,(1,2)])
    draws=gel_clr_draws(20000,rf.k)
    ci=gel_ar_set(rf,.05)
    printf("AR set components (lower upper left_inf right_inf), %g rows\n", rows(ci))
    for(j=1;j<=rows(ci);j++) printf("    %g %g %g %g\n", ci[j,1], ci[j,2], ci[j,3], ci[j,4])
    mismatch=0; npts=0
    for(beta=-2;beta<=2;beta=beta+.01) {
        tt=gel_rf_test(rf,beta,draws); inside=0
        for(j=1;j<=rows(ci);j++) {
            inside=inside | ((ci[j,3] | beta>=ci[j,1]) & (ci[j,4] | beta<=ci[j,2]))
        }
        if(inside!=(tt[2]>=.05)) {
            if(mismatch<5) printf("    mismatch at beta %g: inside %g, AR_F_p %g\n", beta, inside, tt[2])
            mismatch++
        }
        npts++
    }
    fails=fails+gel_check(sprintf("AR inversion vs grid over %g points", npts), mismatch==0, mismatch)
    printf("\nTOTAL FAILED CHECKS: %g\n", fails)
    st_numscalar("GEL_FAILS",fails)
}
gel_selftest()
end
capture which ivreg2
if _rc == 0 {
    capture noisily ivreg2 y (x=z1 z2), cue robust
    scalar b_reference=_b[x]
    capture noisily geliv y, endog(x) instruments(z1 z2) method(cue)
    if _rc {
        display as error "geliv cue failed with rc " _rc
        scalar GEL_FAILS=scalar(GEL_FAILS)+1
    }
    else {
        display "ivreg2 cue: " scalar(b_reference) "   geliv cue: " _b[x] "   diff: " abs(_b[x]-scalar(b_reference))
        if abs(_b[x]-scalar(b_reference))<1e-6 display as result "Overidentified CUE benchmark passed."
        else scalar GEL_FAILS=scalar(GEL_FAILS)+1
    }
}
else display "ivreg2 benchmark skipped: package not installed."
assert scalar(GEL_FAILS)==0
display as result "Mata synthetic checks passed."
