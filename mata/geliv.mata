version 16.0
mata:
/* GELIV 0.1. Linear IV GEL, independent observations, HC0 moments.
   No n x m moment matrix is formed. X and Z may be Stata views.
   Internal b = diag(sx)*theta; internal lambda = diag(sz)*lambda_raw.
   All functions use the gel_ prefix. */
struct gel_data {
    real colvector y
    real matrix X, Z
    real rowvector sx, sz
    real scalar n, p, m, block, maxinner, tol, method
    real scalar cells, e
    real colvector cid, ord, cnt
    real matrix info, E, Xc, Zc
}
struct gel_dual {
    real colvector lambda, gradient
    real matrix H
    real scalar f, ok, iter, error, ridge, margin, mass
}
struct gel_result {
    real rowvector b
    real colvector lambda, weights
    real matrix V, starts
    real scalar q, Jgel, Jdistance, converged, inner_error, outer_error
    real scalar ess, max_nw, min_nw, negative, moment_error
}

real matrix gel_solve(real matrix A, real matrix B)
{
    real matrix L
    L=cholesky((A+A')/2)
    if (hasmissing(L)) return(J(cols(A),cols(B),.))
    return(cholsolve((A+A')/2,B))
}

struct gel_data scalar gel_data_init(real colvector y, real matrix X,
    real matrix Z, real scalar method, real scalar block)
{
    struct gel_data scalar d
    real scalar a,c
    if (rows(y)!=rows(X) | rows(X)!=rows(Z) | cols(Z)<cols(X)) _error(3200)
    if (rows(y)<=cols(Z) | !anyof((-1,0,1),method) | block<1) _error(3498)
    if (hasmissing(y) | hasmissing(X) | hasmissing(Z)) _error(3351)
    d.y=y; d.X=X; d.Z=Z
    d.n=rows(y); d.p=cols(X); d.m=cols(Z)
    d.method=method; d.block=block; d.tol=1e-10; d.maxinner=150
    d.cells=0; d.e=0
    d.sx=J(1,d.p,0); d.sz=J(1,d.m,0)
    for(a=1;a<=d.n;a=a+block) {
        c=min((a+block-1,d.n))
        d.sx=d.sx+colsum(X[|a,1\c,d.p|]:^2)
        d.sz=d.sz+colsum(Z[|a,1\c,d.m|]:^2)
    }
    d.sx=sqrt(d.sx/d.n); d.sz=sqrt(d.sz/d.n)
    if (min(d.sx)<=0 | min(d.sz)<=0) _error(3498)
    return(d)
}

struct gel_data scalar gel_data_cells(real colvector y, real matrix E,
    real matrix W, real matrix Z, real colvector cell, real scalar method,
    real scalar block)
{
    struct gel_data scalar d
    real colvector cs,first,co
    real scalar a,c,k
    if (rows(y)!=rows(E) | rows(E)!=rows(W) | rows(W)!=rows(Z) | rows(Z)!=rows(cell)) _error(3200)
    if (cols(E)<1 | cols(W)<1 | cols(Z)<cols(E)+cols(W)) _error(3200)
    if (!anyof((-1,0,1),method) | block<1) _error(3498)
    if (hasmissing(y) | hasmissing(E) | hasmissing(W) | hasmissing(Z) | hasmissing(cell)) _error(3351)
    d.n=rows(y); d.e=cols(E); d.p=cols(E)+cols(W); d.m=cols(Z)
    if (d.n<=d.m) _error(3498)
    d.method=method; d.block=block; d.tol=1e-10; d.maxinner=150
    d.ord=order(cell,1)
    cs=cell[d.ord]
    d.info=panelsetup(cs,1)
    d.cells=rows(d.info)
    d.cnt=(d.info[.,2]-d.info[.,1]):+1
    d.cid=J(d.n,1,.)
    for(k=1;k<=d.cells;k++) {
        d.cid[|d.info[k,1]\d.info[k,2]|]=J(d.cnt[k],1,k)
    }
    first=d.ord[d.info[.,1]]
    d.Zc=Z[first,.]
    d.Xc=W[first,.]
    co=J(d.n,1,.)
    co[d.ord]=d.cid
    for(a=1;a<=d.n;a=a+block) {
        c=min((a+block-1,d.n))
        if (max(abs(Z[|a,1\c,d.m|]-d.Zc[co[|a\c|],.]))>0) _error(3498)
        if (max(abs(W[|a,1\c,cols(W)|]-d.Xc[co[|a\c|],.]))>0) _error(3498)
    }
    d.y=y[d.ord]
    d.E=E[d.ord,.]
    d.sx=(sqrt(quadcolsum(d.E:^2)/d.n),sqrt(quadcolsum(d.cnt:*(d.Xc:^2))/d.n))
    d.sz=sqrt(quadcolsum(d.cnt:*(d.Zc:^2))/d.n)
    if (min(d.sx)<=0 | min(d.sz)<=0) _error(3498)
    return(d)
}

real colvector gel_resid_cells(struct gel_data scalar d, real rowvector theta)
{
    real colvector fc
    fc=d.Xc*theta[|d.e+1\d.p|]'
    return(d.y-d.E*theta[|1\d.e|]'-fc[d.cid])
}

struct gel_dual scalar gel_stats_cells(struct gel_data scalar d, real rowvector b,
    real colvector lam, real scalar derivatives)
{
    struct gel_dual scalar s
    real matrix Zs
    real colvector u,t,v,w
    real scalar shift,sw
    s.lambda=lam; s.ok=0; s.iter=0; s.ridge=0; s.error=.
    s.f=0; s.gradient=J(d.m,1,0); s.H=J(d.m,d.m,0); s.margin=.; s.mass=0
    Zs=d.Zc:/d.sz
    u=gel_resid_cells(d,b:/d.sx)
    t=Zs*lam
    v=u:*t[d.cid]
    if(hasmissing(v)) {
        s.f=.
        return(s)
    }
    if(d.method==0) {
        shift=max(v)
        w=exp(v:-shift)
        sw=quadsum(w)
        if(sw<=0 | sw>=.) {
            s.f=.
            return(s)
        }
        s.f=shift+ln(sw/d.n)
        if(derivatives) {
            s.gradient=quadcross(Zs,panelsum(w:*u,d.info))/sw
            s.H=quadcross(Zs,panelsum(w:*(u:^2),d.info),Zs)/sw
            s.H=s.H-s.gradient*s.gradient'
        }
    }
    else if(d.method==-1) {
        if(min(1:-v)<=1e-14) {
            s.f=.
            return(s)
        }
        s.margin=min(1:-v)
        s.f=-quadsum(ln(1:-v))/d.n
        if(derivatives) {
            w=1:/(1:-v)
            s.mass=quadsum(w)/d.n
            s.gradient=quadcross(Zs,panelsum(w:*u,d.info))/d.n
            s.H=quadcross(Zs,panelsum((w:^2):*(u:^2),d.info),Zs)/d.n
        }
    }
    else {
        s.f=quadsum(v:+(0.5*(v:^2)))/d.n
        if(derivatives) {
            s.gradient=quadcross(Zs,panelsum((1:+v):*u,d.info))/d.n
            s.H=quadcross(Zs,panelsum(u:^2,d.info),Zs)/d.n
        }
    }
    s.ok=!hasmissing(s.f)
    if(derivatives) {
        s.error=max(abs(s.gradient))
        if(d.method==-1) s.error=max((s.error,s.error/s.mass))
    }
    return(s)
}

real rowvector gel_grad_cells(struct gel_data scalar d, real rowvector b,
    struct gel_dual scalar s)
{
    real colvector u,t,tt,v,w
    u=gel_resid_cells(d,b:/d.sx)
    t=(d.Zc:/d.sz)*s.lambda
    tt=t[d.cid]
    v=u:*tt
    if(d.method==0) w=exp(v:-s.f)/d.n
    else if(d.method==-1) w=(1:/(1:-v))/d.n
    else w=(1:+v)/d.n
    return((quadcross(d.E:/d.sx[|1\d.e|],w:*tt)',quadcross(d.Xc:/d.sx[|d.e+1\d.p|],t:*panelsum(w,d.info))'))
}

real matrix gel_vcov_cells(struct gel_data scalar d, real rowvector b)
{
    real matrix Zs,D,O,A,V
    real colvector u
    Zs=d.Zc:/d.sz
    u=gel_resid_cells(d,b:/d.sx)
    D=-(quadcross(Zs,panelsum(d.E:/d.sx[|1\d.e|],d.info)),quadcross(Zs,d.cnt,d.Xc:/d.sx[|d.e+1\d.p|]))/d.n
    O=quadcross(Zs,panelsum(u:^2,d.info),Zs)/d.n
    A=gel_solve(O,D)
    if(hasmissing(A)) return(J(d.p,d.p,.))
    V=gel_solve(D'*A,I(d.p))/d.n
    return(V:/(d.sx'*d.sx))
}

void gel_finish_cells(struct gel_data scalar d, struct gel_result scalar r,
    struct gel_dual scalar s)
{
    real colvector u,t,v,w,mg,rw
    real scalar den
    u=gel_resid_cells(d,r.b)
    t=d.Zc*r.lambda
    v=u:*t[d.cid]
    if(d.method==0) w=exp(v:-s.f)/d.n
    else if(d.method==-1) w=1:/(1:-v)
    else w=1:+v
    den=quadsum(w)
    if(abs(den)<1e-12 | hasmissing(den)) {
        r.converged=0
        return
    }
    w=w/den
    mg=quadcross(d.Zc:/d.sz,panelsum(w:*u,d.info))
    r.moment_error=max(abs(mg))
    r.negative=sum(w:<0)
    r.max_nw=max(d.n*w)
    r.min_nw=min(d.n*w)
    if(r.negative==0) r.ess=1/quadsum(w:^2)
    if(r.moment_error>1e-7 | hasmissing(r.moment_error)) r.converged=0
    rw=J(d.n,1,.)
    rw[d.ord]=w
    r.weights=rw
}

real rowvector gel_2sls_cells(struct gel_data scalar d)
{
    real matrix Zs,ZZ,ZX
    real colvector Zy
    Zs=d.Zc:/d.sz
    ZZ=quadcross(Zs,d.cnt,Zs)
    ZX=(quadcross(Zs,panelsum(d.E:/d.sx[|1\d.e|],d.info)),quadcross(Zs,d.cnt,d.Xc:/d.sx[|d.e+1\d.p|]))
    Zy=quadcross(Zs,panelsum(d.y,d.info))
    return(gel_solve(ZX'*gel_solve(ZZ,ZX),ZX'*gel_solve(ZZ,Zy))':/d.sx)
}

real rowvector gel_ols_cells(struct gel_data scalar d)
{
    real matrix EW,XX
    real colvector Xy
    EW=quadcross(panelsum(d.E,d.info),d.Xc)
    XX=(quadcross(d.E,d.E),EW\EW',quadcross(d.Xc,d.cnt,d.Xc))
    Xy=quadcross(d.E,d.y)\quadcross(d.Xc,panelsum(d.y,d.info))
    return(gel_solve(XX,Xy)')
}

/* Convex internal objective: ET log mean exp(v), EL -mean log(1-v),
   CUE mean(v+v^2/2). ET uses a global shift over streaming blocks. */
struct gel_dual scalar gel_stats(struct gel_data scalar d, real rowvector b,
    real colvector lam, real scalar derivatives)
{
    struct gel_dual scalar s
    real scalar a,c,shift,sw
    real matrix Z,G
    real colvector u,v,w
    if(d.cells) return(gel_stats_cells(d,b,lam,derivatives))
    s.lambda=lam; s.ok=0; s.iter=0; s.ridge=0; s.error=.
    s.f=0; s.gradient=J(d.m,1,0); s.H=J(d.m,d.m,0); s.margin=.; s.mass=0
    shift=0; sw=0
    if(d.method==0) {
        shift=-1e300
        for(a=1;a<=d.n;a=a+d.block) {
            c=min((a+d.block-1,d.n))
            u=d.y[|a\c|]-d.X[|a,1\c,d.p|]*(b:/d.sx)'
            v=u:*(d.Z[|a,1\c,d.m|]*(lam:/d.sz'))
            shift=max((shift,max(v)))
        }
    }
    for(a=1;a<=d.n;a=a+d.block) {
        c=min((a+d.block-1,d.n))
        Z=d.Z[|a,1\c,d.m|]:/d.sz
        u=d.y[|a\c|]-d.X[|a,1\c,d.p|]*(b:/d.sx)'
        G=Z:*u; v=G*lam
        if(hasmissing(v)) {
            s.f=.; return(s)
        }
        if(d.method==0) {
            w=exp(v:-shift); sw=sw+sum(w)
            if(derivatives) {
                s.gradient=s.gradient+quadcross(G,w)
                s.H=s.H+quadcross(G,w,G)
            }
        }
        else if(d.method==-1) {
            if(min(1:-v)<=1e-14) {
                s.f=.; return(s)
            }
            s.margin=min((s.margin,min(1:-v)))
            s.f=s.f-sum(ln(1:-v))/d.n
            if(derivatives) {
                w=1:/(1:-v)
                s.mass=s.mass+sum(w)/d.n
                s.gradient=s.gradient+quadcross(G,w)/d.n
                s.H=s.H+quadcross(G,w:^2,G)/d.n
            }
        }
        else {
            s.f=s.f+sum(v:+(0.5*(v:^2)))/d.n
            if(derivatives) {
                s.gradient=s.gradient+quadcross(G,1:+v)/d.n
                s.H=s.H+quadcross(G,G)/d.n
            }
        }
    }
    if(d.method==0) {
        if(sw<=0 | sw>=.) {
            s.f=.; return(s)
        }
        s.f=shift+ln(sw/d.n)
        if(derivatives) {
            s.gradient=s.gradient/sw
            s.H=s.H/sw-s.gradient*s.gradient'
        }
    }
    s.ok=!hasmissing(s.f)
    if(derivatives) {
        s.error=max(abs(s.gradient))
        /* EL raw gradients can approach zero along an escaping dual ray.
           The normalized moment root must converge as well. */
        if(d.method==-1) s.error=max((s.error,s.error/s.mass))
    }
    return(s)
}

struct gel_dual scalar gel_inner(struct gel_data scalar d, real rowvector b,
    real scalar force_newton)
{
    struct gel_dual scalar s,t
    real colvector lam,step
    real scalar it,j,rate,descent,ridge,totalridge
    lam=J(d.m,1,0); totalridge=0
    s=gel_stats(d,b,lam,1)
    if(d.method==1 & !force_newton) {
        lam=-gel_solve(s.H,s.gradient)
        if(hasmissing(lam)) {
            s.ok=0; return(s)
        }
        s=gel_stats(d,b,lam,1)
        s.ok=s.ok & s.error<d.tol; return(s)
    }
    for(it=1;it<=d.maxinner;it++) {
        s=gel_stats(d,b,lam,1); s.iter=it
        if(!s.ok) return(s)
        if(s.error<d.tol) {
            s.ridge=totalridge; return(s)
        }
        ridge=0
        for(j=1;j<=9;j++) {
            step=-gel_solve(s.H+ridge*I(d.m),s.gradient)
            if(!hasmissing(step)) break
            ridge=max((1e-12,10*ridge,1e-12*trace(s.H)/d.m))
        }
        if(hasmissing(step)) {
            s.ok=0; return(s)
        }
        totalridge=totalridge+(ridge>0)
        descent=quadcross(s.gradient,step); rate=1
        if(descent>=0) {
            s.ok=0; return(s)
        }
        for(j=1;j<=50;j++) {
            t=gel_stats(d,b,lam+rate*step,0)
            if(t.ok) {
                if(t.f<=s.f+1e-4*rate*descent+1e-15) break
            }
            rate=rate/2
        }
        if(j>50) {
            s.ok=0; return(s)
        }
        lam=lam+rate*step
        if(max(abs(lam))>1e12) {
            s.ok=0; return(s)
        }
    }
    s=gel_stats(d,b,lam,1); s.ok=0; s.iter=d.maxinner
    return(s)
}

void gel_eval(real scalar todo, real rowvector b, struct gel_data scalar d,
    q, grad, H)
{
    struct gel_dual scalar s
    real scalar a,c
    real colvector u,t,v,w
    real matrix X,Z
    s=gel_inner(d,b,0); q=.; grad=J(1,d.p,.)
    if(!s.ok) return
    q=-s.f
    if(!todo) return
    if(d.cells) {
        grad=gel_grad_cells(d,b,s)
        return
    }
    grad=J(1,d.p,0)
    for(a=1;a<=d.n;a=a+d.block) {
        c=min((a+d.block-1,d.n))
        X=d.X[|a,1\c,d.p|]:/d.sx; Z=d.Z[|a,1\c,d.m|]:/d.sz
        u=d.y[|a\c|]-X*b'; t=Z*s.lambda; v=u:*t
        if(d.method==0) w=exp(v:-s.f)/d.n
        else if(d.method==-1) w=(1:/(1:-v))/d.n
        else w=(1:+v)/d.n
        grad=grad+quadcross(X,w:*t)'
    }
}

real matrix gel_vcov(struct gel_data scalar d, real rowvector b)
{
    real matrix D,O,Z,X,G,A,V
    real colvector u
    real scalar a,c
    if(d.cells) return(gel_vcov_cells(d,b))
    D=J(d.m,d.p,0); O=J(d.m,d.m,0)
    for(a=1;a<=d.n;a=a+d.block) {
        c=min((a+d.block-1,d.n))
        Z=d.Z[|a,1\c,d.m|]:/d.sz; X=d.X[|a,1\c,d.p|]:/d.sx
        u=d.y[|a\c|]-X*b'; G=Z:*u
        D=D-quadcross(Z,X)/d.n; O=O+quadcross(G,G)/d.n
    }
    A=gel_solve(O,D)
    if(hasmissing(A)) return(J(d.p,d.p,.))
    V=gel_solve(D'*A,I(d.p))/d.n
    return(V:/(d.sx'*d.sx))
}

struct gel_result scalar gel_fit(struct gel_data scalar d, real matrix starts)
{
    struct gel_result scalar r
    struct gel_dual scalar s
    transmorphic S
    real rowvector b,g,best
    real matrix H
    real colvector u,v,w,mg
    real scalar h,rc,q,a,c,den,oe
    r.q=.; r.converged=0; r.b=J(1,d.p,.)
    r.V=J(d.p,d.p,.); r.weights=J(d.n,1,.); r.lambda=J(d.m,1,.)
    r.Jgel=.; r.Jdistance=.; r.inner_error=.; r.outer_error=.
    r.ess=.; r.max_nw=.; r.min_nw=.; r.negative=.; r.moment_error=.
    r.starts=J(rows(starts),4,.)
    for(h=1;h<=rows(starts);h++) {
        S=optimize_init()
        optimize_init_evaluator(S,&gel_eval())
        optimize_init_evaluatortype(S,"d1")
        optimize_init_argument(S,1,d)
        optimize_init_which(S,"min")
        optimize_init_technique(S,"bfgs")
        optimize_init_params(S,starts[h,.]:*d.sx)
        optimize_init_conv_maxiter(S,400)
        optimize_init_conv_ptol(S,1e-10)
        optimize_init_conv_vtol(S,1e-12)
        optimize_init_conv_nrtol(S,1e-9)
        optimize_init_tracelevel(S,"none")
        rc=_optimize(S)
        if(rc) {
            r.starts[h,.]=(rc,0,.,.); continue
        }
        b=optimize_result_params(S)
        gel_eval(1,b,d,q,g,H); oe=max(abs(g))
        r.starts[h,.]=(rc,optimize_result_converged(S),q,oe)
        if(!optimize_result_converged(S) | oe>1e-6 | hasmissing(q)) continue
        if(q<r.q) {
            r.q=q; best=b; r.outer_error=oe; r.converged=1
        }
    }
    if(!r.converged) return(r)
    r.b=best:/d.sx; s=gel_inner(d,best,0)
    r.lambda=s.lambda:/d.sz'; r.inner_error=s.error
    r.V=gel_vcov(d,best); r.Jdistance=2*d.n*r.q
    if(d.method==0) r.Jgel=2*d.n*(1-exp(-r.q))
    else r.Jgel=r.Jdistance
    if(d.cells) {
        gel_finish_cells(d,r,s)
        return(r)
    }
    den=0
    for(a=1;a<=d.n;a=a+d.block) {
        c=min((a+d.block-1,d.n))
        u=d.y[|a\c|]-d.X[|a,1\c,d.p|]*r.b'
        v=u:*(d.Z[|a,1\c,d.m|]*r.lambda)
        if(d.method==0) w=exp(v:-s.f)/d.n
        else if(d.method==-1) w=1:/(1:-v)
        else w=1:+v
        r.weights[|a\c|]=w; den=den+sum(w)
    }
    if(abs(den)<1e-12 | hasmissing(den)) {
        r.converged=0; return(r)
    }
    r.weights=r.weights/den; mg=J(d.m,1,0)
    for(a=1;a<=d.n;a=a+d.block) {
        c=min((a+d.block-1,d.n))
        u=d.y[|a\c|]-d.X[|a,1\c,d.p|]*r.b'
        mg=mg+quadcross(d.Z[|a,1\c,d.m|]:/d.sz,r.weights[|a\c|]:*u)
    }
    r.moment_error=max(abs(mg)); r.negative=sum(r.weights:<0)
    r.max_nw=max(d.n*r.weights); r.min_nw=min(d.n*r.weights)
    if(r.negative==0) r.ess=1/sum(r.weights:^2)
    if(r.moment_error>1e-7 | hasmissing(r.moment_error)) r.converged=0
    return(r)
}

/* HC0 2SLS starting values, all controls included in X and Z. */
real rowvector gel_2sls(struct gel_data scalar d)
{
    real matrix ZZ,ZX
    real colvector Zy
    if(d.cells) return(gel_2sls_cells(d))
    ZZ=quadcross(d.Z,d.Z):/(d.sz'*d.sz)
    ZX=quadcross(d.Z,d.X):/(d.sz'*d.sx)
    Zy=quadcross(d.Z,d.y):/d.sz'
    return(gel_solve(ZX'*gel_solve(ZZ,ZX),ZX'*gel_solve(ZZ,Zy))':/d.sx)
}

void gel_stata(string scalar yname,string scalar xnames,string scalar znames,
    string scalar sample,real scalar method,real scalar block,
    string scalar startname,string scalar bname,string scalar vname,
    string scalar diagname,string scalar weightsname,string scalar startsname,
    string scalar cellname,real scalar ne)
{
    real colvector y,cell
    real matrix X,Z,starts,E,W
    string rowvector xn
    struct gel_data scalar d
    struct gel_result scalar r
    st_view(y,.,yname,sample); st_view(X,.,tokens(xnames),sample)
    st_view(Z,.,tokens(znames),sample)
    if(cellname!="") {
        xn=tokens(xnames)
        if(ne<1 | ne>=cols(xn)) _error(3498)
        st_view(E,.,xn[|1,1\1,ne|],sample)
        st_view(W,.,xn[|1,ne+1\1,cols(xn)|],sample)
        cell=st_data(.,cellname,sample)
        d=gel_data_cells(y,E,W,Z,cell,method,block)
    }
    else d=gel_data_init(y,X,Z,method,block)
    if(startname!="") starts=st_matrix(startname)
    else if(d.cells) starts=gel_2sls(d)\gel_ols_cells(d)
    else starts=gel_2sls(d)\qrsolve(X,y)'
    if(cols(starts)!=d.p | hasmissing(starts)) _error(3498)
    r=gel_fit(d,starts)
    st_matrix(bname,r.b); st_matrix(vname,r.V); st_matrix(startsname,r.starts)
    st_matrix(diagname,(r.converged,r.q,r.Jgel,r.Jdistance,r.inner_error,
        r.outer_error,r.moment_error,r.ess,r.min_nw,r.max_nw,r.negative,d.n,d.m,d.p))
    if(weightsname!="" & r.converged) st_store(.,weightsname,sample,r.weights)
}
end
