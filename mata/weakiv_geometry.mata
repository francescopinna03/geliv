version 16.0
mata:
/* Independent observations only. Full-vector robust tests and a separate
   scalar homoskedastic IV module. No generic ET/EL-CLR claim is made. */

/* Full null theta=theta0. Output AR, p_AR, K, p_K, J=AR-K, p_J.
   K uses the covariance-orthogonalized Jacobian, not raw GEL gradient.
   Fixed m,p asymptotics, centered HC0 covariance, full-rank adjusted D. */
real rowvector gel_full_tests(struct gel_data scalar d,real rowvector theta)
{
    real matrix O,D,G,Z,X,C,A,B
    real colvector g,u,h,t
    real scalar a,c,j,AR,K
    if(d.cells) _error(3498)
    g=J(d.m,1,0); O=J(d.m,d.m,0); D=J(d.m,d.p,0)
    for(a=1;a<=d.n;a=a+d.block) {
        c=min((a+d.block-1,d.n))
        Z=d.Z[|a,1\c,d.m|]:/d.sz; X=d.X[|a,1\c,d.p|]:/d.sx
        u=d.y[|a\c|]-d.X[|a,1\c,d.p|]*theta'
        G=Z:*u; g=g+colsum(G)'/d.n
        O=O+quadcross(G,G)/d.n; D=D-quadcross(Z,X)/d.n
    }
    O=O-g*g'; h=gel_solve(O,g)
    if(hasmissing(h)) return(J(1,6,.))
    A=D
    for(j=1;j<=d.p;j++) {
        C=J(d.m,d.m,0)
        for(a=1;a<=d.n;a=a+d.block) {
            c=min((a+d.block-1,d.n))
            Z=d.Z[|a,1\c,d.m|]:/d.sz
            u=d.y[|a\c|]-d.X[|a,1\c,d.p|]*theta'
            C=C-quadcross(Z,(d.X[|a,j\c,j|]/d.sx[j]):*u,Z)/d.n
        }
        C=C-D[.,j]*g'; A[.,j]=D[.,j]-C*h
    }
    B=gel_solve(O,A); t=A'*h
    AR=d.n*quadcross(g,h)
    if(hasmissing(B)) return((AR,chi2tail(d.m,AR),.,.,.,.))
    K=d.n*quadcross(t,gel_solve(A'*B,t))
    if(hasmissing(K)) return((AR,chi2tail(d.m,AR),.,.,.,.))
    K=max((0,min((K,AR))))
    if(d.m==d.p) return((AR,chi2tail(d.m,AR),K,chi2tail(d.p,K),AR-K,.))
    return((AR,chi2tail(d.m,AR),K,chi2tail(d.p,K),AR-K,chi2tail(d.m-d.p,AR-K)))
}

struct gel_rf {
    real matrix R,E,Sigma
    real scalar n,k,l,df
}

/* Residualizing W is valid HERE for linear reduced-form homoskedastic tests.
   It is deliberately not used in the GEL estimator. W must contain a constant
   when one is wanted. Z contains excluded instruments only. */
struct gel_rf scalar gel_rf_init(real colvector y,real colvector x,
    real matrix W,real matrix Z)
{
    struct gel_rf scalar r
    real matrix Y,Zr,C,E,B
    Y=(y,x); Zr=Z
    if(cols(W)>0) {
        Y=Y-W*qrsolve(W,Y); Zr=Z-W*qrsolve(W,Z)
    }
    B=quadcross(Zr,Zr); C=quadcross(Zr,Y)
    r.R=C'*gel_solve(B,C); E=quadcross(Y,Y)-r.R
    r.n=rows(y); r.k=cols(Z); r.l=cols(W); r.df=r.n-r.l-r.k
    if(r.df<=0 | hasmissing(r.R)) _error(3498)
    r.E=(E+E')/2; r.Sigma=r.E/r.df
    if(hasmissing(gel_solve(r.Sigma,I(2)))) _error(3498)
    return(r)
}

/* Common random numbers are supplied by caller: col 1 ~ chi2(1),
   col 2 ~ independent chi2(k-1). k=1 uses zeros in col 2. */
real matrix gel_clr_draws(real scalar reps,real scalar k)
{
    real matrix A
    A=rnormal(reps,1,0,1):^2
    if(k>1) return((A,rowsum(rnormal(reps,k-1,0,1):^2)))
    return((A,J(reps,1,0)))
}

/* Output AR_chi, AR_F_p, K, K_p, CLR, CLR_p_MC, MCse, TT.
   Canonical LR conditioning is exact if Sigma is KNOWN Gaussian covariance;
   using unrestricted reduced-form Sigma_hat gives asymptotic CLR validity.
   AR_F is finite-sample exact under normal homoskedastic errors. */
real rowvector gel_rf_test(struct gel_rf scalar r,real scalar beta,
    real matrix draws)
{
    real colvector b,a,t,Ldraw,U
    real scalar SS,TT,ST,K,LR,pc,denb,dent
    b=(1\-beta); a=(beta\1); t=gel_solve(r.Sigma,a)
    denb=quadcross(b,r.Sigma*b); dent=quadcross(a,t)
    SS=quadcross(b,r.R*b)/denb; TT=quadcross(t,r.R*t)/dent
    ST=quadcross(b,r.R*t)/sqrt(denb*dent)
    K=0
    if(TT>1e-14) K=ST^2/TT
    K=min((SS,max((0,K))))
    LR=(SS-TT+sqrt((SS-TT)^2+4*ST^2))/2
    U=rowsum(draws)
    Ldraw=((U:-TT):+sqrt(((U:-TT):^2):+(4*TT*draws[.,1])))/2
    pc=(1+sum(Ldraw:>=LR))/(rows(draws)+1)
    return((SS,Ftail(r.k,r.df,SS/r.k),K,chi2tail(1,K),LR,pc,
        sqrt(pc*(1-pc)/(rows(draws)+1)),TT))
}

/* Exact algebraic inversion of classical AR F test. Rows are components.
   Columns lower,upper,left_infinite,right_infinite. Infinity is represented
   by . plus an explicit flag. 0 rows means empty set. */
real matrix gel_ar_set(struct gel_rf scalar r,real scalar alpha)
{
    real matrix M
    real scalar h,A,B,C,disc,lo,hi,tol
    h=r.k*invFtail(r.k,r.df,alpha)/r.df
    M=r.R-h*r.E; A=M[2,2]; B=-2*M[1,2]; C=M[1,1]
    tol=1e-12*max((1,max(abs(M))))
    if(abs(A)<tol) {
        if(abs(B)<tol) {
            if(C<=0) return((.,.,1,1))
            return(J(0,4,.))
        }
        if(B>0) return((.,-C/B,1,0))
        return((-C/B,.,0,1))
    }
    disc=B^2-4*A*C
    if(disc<0) {
        if(A<0) return((.,.,1,1))
        return(J(0,4,.))
    }
    /* Stable quadratic roots; avoid cancellation of -B + sqrt(disc). */
    if(B>=0) h=-0.5*(B+sqrt(disc))
    else h=-0.5*(B-sqrt(disc))
    if(abs(h)<tol) {
        lo=-B/(2*A); hi=lo
    }
    else {
        lo=min((h/A,C/h)); hi=max((h/A,C/h))
    }
    if(A>0) return((lo,hi,0,0))
    return((.,lo,1,0)\(hi,.,0,1))
}

/* Grid is a set of evaluated points, NOT a certificate about its tails.
   Use gel_ar_set for certified AR tails. K/CLR grid edges remain unresolved. */
real matrix gel_rf_grid(struct gel_rf scalar r,real colvector grid,real matrix draws)
{
    real matrix out
    real scalar j
    out=J(rows(grid),9,.)
    for(j=1;j<=rows(grid);j++) out[j,.]=(grid[j],gel_rf_test(r,grid[j],draws))
    return(out)
}

/* Scalar p=1 GEL null inversion without nuisance parameters.
   The full-null reference is chi2(m), not chi2(m-p), and the minimized
   unrestricted criterion must NOT be subtracted under weak identification. */
real matrix gel_null_grid(struct gel_data scalar d,real colvector grid)
{
    struct gel_dual scalar s
    real matrix out
    real scalar j,T
    if(d.p!=1) _error(3498)
    out=J(rows(grid),4,.)
    for(j=1;j<=rows(grid);j++) {
        s=gel_inner(d,grid[j]*d.sx,0)
        out[j,1]=grid[j]; out[j,4]=s.ok
        if(!s.ok) continue
        T=-2*d.n*s.f
        if(d.method==0) T=2*d.n*(1-exp(s.f))
        out[j,2]=T; out[j,3]=chi2tail(d.m,max((0,T)))
    }
    return(out)
}
end
