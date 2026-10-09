"""Independent numerical checks, NOT a substitute for executing Mata/Stata.
Run: OPENBLAS_NUM_THREADS=1 python3 python/verify_math.py
Dependencies: numpy, scipy. Output output/math_validation.json.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import chi2, f

RNG = np.random.default_rng(7102026)


def dual_stats(G, lam, method):
    v = G @ lam
    n = len(G)
    if method == 'et':
        val = logsumexp(v) - np.log(n)
        w = np.exp(v-logsumexp(v))
        g = G.T @ w
        H = G.T @ (w[:, None]*G)-np.outer(g, g)
    elif method == 'el':
        if np.min(1-v) <= 0:
            return np.inf, None, None
        w = 1/(1-v)
        val = -np.log(1-v).mean()
        g = G.T @ w/n
        H = G.T @ (w[:, None]**2*G)/n
    else:
        val = np.mean(v+v*v/2)
        g = G.T @ (1+v)/n
        H = G.T @ G/n
    return val, g, H


def solve_dual(G, method):
    lam = np.zeros(G.shape[1])
    for _ in range(150):
        val, grad, H = dual_stats(G, lam, method)
        if np.max(np.abs(grad)) < 1e-11:
            return lam, val
        step = np.linalg.solve(H, -grad)
        rate = 1.
        for _ in range(60):
            new = dual_stats(G, lam+rate*step, method)[0]
            if np.isfinite(new) and new <= val+1e-4*rate*(grad@step)+1e-15:
                break
            rate /= 2
        else:
            raise RuntimeError('line search')
        lam += rate*step
    raise RuntimeError('dual did not converge')


def criterion(y, X, Z, b, method):
    G = Z*(y-X@b)[:, None]
    lam, val = solve_dual(G, method)
    v = G@lam
    if method == 'et':
        w = np.exp(v-logsumexp(v))
    elif method == 'el':
        w = 1/(1-v)/len(y)
    else:
        w = (1+v)/len(y)
    return -val, X.T@(w*(Z@lam)), lam, w/w.sum()


def fd(fun, x, h=1e-5):
    E = np.eye(len(x))*h
    return np.array([(fun(x+e)-fun(x-e))/(2*h) for e in E]).T


def rf(y, x, W, Z):
    Y = np.column_stack((y, x))
    Y = Y-W@np.linalg.lstsq(W, Y, rcond=None)[0]
    Z = Z-W@np.linalg.lstsq(W, Z, rcond=None)[0]
    C = Z.T@Y
    R = C.T@np.linalg.solve(Z.T@Z, C)
    E = Y.T@Y-R
    df = len(y)-W.shape[1]-Z.shape[1]
    return R, E, df


def ar_components(R, E, df, k, alpha=.05):
    M = R-k*f.isf(alpha, k, df)/df*E
    A, B, C = M[1, 1], -2*M[0, 1], M[0, 0]
    if abs(A) < 1e-12:
        if abs(B) < 1e-12:
            return [(-np.inf, np.inf)] if C <= 0 else []
        return [(-np.inf, -C/B)] if B > 0 else [(-C/B, np.inf)]
    disc = B*B-4*A*C
    if disc < 0:
        return [(-np.inf, np.inf)] if A < 0 else []
    roots = np.sort(np.roots([A, B, C]))
    return [tuple(roots)] if A > 0 else [(-np.inf, roots[0]), (roots[1], np.inf)]


def run():
    out = {'scope': 'Python mathematical checks only; Stata/Mata NOT executed', 'checks': []}

    def check(name, error, tol):
        error = float(error)
        out['checks'].append(dict(name=name, error=error, tolerance=tol, passed=error <= tol))
        assert error <= tol, (name, error, tol)

    # Moderate nonzero means, interior feasible moment hull.
    G = RNG.normal(size=(240, 4)) + np.array([.1, -.07, .04, .08])
    for method in ['et', 'el', 'cue']:
        lam, val = solve_dual(G, method)
        _, grad, H = dual_stats(G, lam, method)
        probe = lam + .01
        check(f'{method}: inner gradient', np.max(np.abs(
            fd(lambda t: dual_stats(G, t, method)[0], probe)-dual_stats(G, probe, method)[1])), 1e-7)
        check(f'{method}: inner Hessian', np.max(np.abs(
            fd(lambda t: dual_stats(G, t, method)[1], probe)-dual_stats(G, probe, method)[2])), 1e-7)
        check(f'{method}: moment root', np.max(np.abs(grad)), 1e-9)
        if method == 'cue':
            exact = -np.linalg.solve(G.T@G/len(G), G.mean(axis=0))
            check('CUE Newton versus closed form', np.max(np.abs(lam-exact)), 1e-10)
        else:
            n = len(G)
            # Independent primal constrained optimizer, 240 probabilities.
            if method == 'et':
                obj = lambda p: np.sum(p*np.log(n*p))
                jac = lambda p: np.log(n*p)+1
                pdual = np.exp(G@lam-logsumexp(G@lam))
            else:
                obj = lambda p: -np.log(n*p).mean()
                jac = lambda p: -1/(n*p)
                pdual = 1/(1-G@lam)/n
            A = np.vstack((np.ones(n), G.T))
            target = np.r_[1., np.zeros(G.shape[1])]
            res = minimize(obj, np.full(n, 1/n), jac=jac, method='SLSQP',
                           constraints={'type': 'eq', 'fun': lambda p: A@p-target,
                                        'jac': lambda p: A},
                           bounds=[(1e-12, 1)]*n, options={'ftol': 1e-12, 'maxiter': 1000})
            assert res.success, res.message
            check(f'{method}: independent primal versus dual weights', np.max(abs(res.x-pdual)), 2e-7)
            check(f'{method}: primal divergence equals profiled criterion', abs(obj(pdual)+val), 1e-9)

    n=450
    Z = np.c_[RNG.normal(size=(n, 3)), np.ones(n)]
    v = RNG.normal(size=n)
    X = np.c_[Z[:, 0]+.3*Z[:, 1]+v, np.ones(n)]
    y = X@np.array([.7, .2])+.5*v+RNG.normal(size=n)
    b = np.array([.72, .22])
    for method in ['et', 'el', 'cue']:
        q, grad, lam, w = criterion(y, X, Z, b, method)
        check(f'{method}: envelope gradient', np.max(abs(
            fd(lambda t: criterion(y, X, Z, t, method)[0], b)-grad)), 1e-7)
        A = RNG.normal(size=(4, 4))+3*np.eye(4)
        check(f'{method}: invertible instrument transformation',
              abs(criterion(y, X, Z@A, b, method)[0]-q), 1e-9)

    # Just identification, uniform probabilities and common coefficient.
    Zj=Z[:, [0, 3]]
    biv=np.linalg.solve(Zj.T@X, Zj.T@y)
    for method in ['et', 'el', 'cue']:
        q, grad, lam, w=criterion(y, X, Zj, biv, method)
        check(f'{method}: just-identified criterion', abs(q), 1e-10)
        check(f'{method}: just-identified uniform weights', max(abs(w-1/n)), 1e-10)

    # Pythagorean identity for another feasible distribution at FIXED theta.
    G=Z*(y-X@b)[:, None]
    lam, val=solve_dual(G, 'et')
    p=np.exp(G@lam-logsumexp(G@lam))
    A=np.vstack([np.ones(n), G.T])
    h=RNG.normal(size=n); h-=A.T@np.linalg.solve(A@A.T,A@h)
    h*=.2*np.min(p)/max(abs(h)); q=p+h
    kl=lambda a,b: np.sum(a*np.log(a/b))
    unif=np.ones(n)/n
    check('ET Pythagorean identity', abs(kl(q,unif)-kl(q,p)-kl(p,unif)), 1e-10)

    # Streaming sums with a single global shift agree with dense calculations.
    shift=max(G@lam); sw=0.; s1=np.zeros(4); s2=np.zeros((4,4))
    for chunk in np.array_split(G, 17):
        w=np.exp(chunk@lam-shift); sw+=sum(w)
        s1+=chunk.T@w; s2+=chunk.T@(w[:,None]*chunk)
    _,g,H=dual_stats(G,lam,'et')
    check('ET block gradient', max(abs(s1/sw-g)), 1e-12)
    check('ET block Hessian', np.max(abs(s2/sw-np.outer(s1/sw,s1/sw)-H)), 1e-12)

    # Analytic AR inversion against pointwise F tests, including weak designs.
    for strength in [0., .05, .4]:
        z=RNG.normal(size=(300,3)); w=np.ones((300,1)); v=RNG.normal(size=300)
        x=strength*z[:,0]+v; y=.7*x+.6*v+RNG.normal(size=300)
        R,E,df=rf(y,x,w,z); components=ar_components(R,E,df,3)
        for beta in np.r_[np.linspace(-100,100,1001),-1e8,1e8]:
            b0=np.array([1,-beta]); stat=(b0@R@b0)/(b0@E@b0/df)/3
            inside=any(lo<=beta<=hi for lo,hi in components)
            check(f'AR inversion strength={strength} beta={beta:g}',
                  float(inside != (f.sf(stat,3,df)>=.05)), 0)
    # Controlled zero-strength null simulation: K and CLR finite known-Sigma.
    B=10000; k=5
    S=RNG.normal(size=(B,k)); T=RNG.normal(size=(B,k))
    ss=np.sum(S*S,axis=1); tt=np.sum(T*T,axis=1); st=np.sum(S*T,axis=1)
    kval=st*st/tt
    check('K zero-strength size versus 5% (MC tolerance)',
          abs(np.mean(chi2.sf(kval,1)<.05)-.05), .012)
    # Conditional critical values for several fixed T norms.
    draws1=RNG.chisquare(1, size=200000); drawsr=RNG.chisquare(k-1,size=200000)
    for t2 in [0.,1.,10.,1000.]:
        u=draws1+drawsr
        sim=(u-t2+np.sqrt((u-t2)**2+4*t2*draws1))/2
        crit=np.quantile(sim,.95)
        lr=(ss-t2+np.sqrt((ss-t2)**2+4*t2*S[:,0]**2))/2
        check(f'CLR conditional size T2={t2}', abs(np.mean(lr>crit)-.05), .012)
    # A dense algebra check of the adjusted full-vector K covariance identity.
    u=y-.7*x; gg=z*u[:,None]; g=gg.mean(axis=0); O=np.cov(gg,rowvar=False,bias=True)
    gi=-z*x[:,None]; D=gi.mean(axis=0)
    C=(gi-D).T@(gg-g)/len(y); dh=D-C@np.linalg.solve(O,g)
    K=len(y)*(dh@np.linalg.solve(O,g))**2/(dh@np.linalg.solve(O,dh))
    AR=len(y)*g@np.linalg.solve(O,g)
    check('K is a projection bounded by AR', max(0.,K-AR,-K), 1e-10)

    out['passed']=all(c['passed'] for c in out['checks'])
    out['count']=len(out['checks'])
    path=Path(__file__).resolve().parents[1]/'output/math_validation.json'
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(out,indent=2))
    print(json.dumps({'passed':out['passed'],'count':out['count'],'report':str(path)}))


if __name__=='__main__':
    run()
