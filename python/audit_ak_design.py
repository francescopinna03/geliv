"""Independent sparse audit of the actual AK design; no GEL estimation here."""
from pathlib import Path
import json
import numpy as np
from scipy.sparse import csr_matrix, hstack

root=Path(__file__).resolve().parents[1]
a=np.loadtxt(root/'data/asciiqob.txt')
raw_n=len(a)
y,x,yob,qob,pob=a.T
yob=np.where(yob<100,yob+1900,yob)
states=np.array([1,2,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,53,54,55,56])
mask=np.isfinite(a).all(axis=1)&(yob>=1930)&(yob<=1939)&np.isin(pob,states)
y,x,yob,qob,pob=[v[mask] for v in [y,x,yob,qob,pob]]
n=len(y); other=states[states!=26]
W=csr_matrix(np.column_stack([*(yob==k for k in range(1931,1940)),*(pob==s for s in other),np.ones(n)]).astype(float))
Z3=csr_matrix(np.column_stack([qob==q for q in [1,2,3]]).astype(float))
Z30=csr_matrix(np.column_stack([(qob==q)&(yob==k) for q in [1,2,3] for k in range(1930,1940)]).astype(float))
Z150=csr_matrix(np.column_stack([(qob==q)&(pob==s) for q in [1,2,3] for s in other]).astype(float))
X=hstack([csr_matrix(x[:,None]),W],format='csr')
XX=(X.T@X).toarray(); Xy=X.T@y
ols=np.linalg.solve(XX,Xy)[0]
result={'raw_n':raw_n,'common_n':n,'controls_including_constant':W.shape[1],
        'ols_education':float(ols),'specifications':[],'scope':'Python sparse crossproduct audit, not Stata or GEL validation'}
for k,Z in [(3,Z3),(30,Z30),(180,hstack([Z30,Z150],format='csr'))]:
    A=hstack([Z,W],format='csr')
    AA=(A.T@A).toarray(); AX=(A.T@X).toarray(); Ay=A.T@y
    rank=int(np.linalg.matrix_rank(AA))
    assert rank==A.shape[1],(k,rank,A.shape[1])
    b=np.linalg.solve(AX.T@np.linalg.solve(AA,AX),AX.T@np.linalg.solve(AA,Ay))
    vx=x-A@np.linalg.solve(AA,A.T@x)
    WW=(W.T@W).toarray()
    vw=x-W@np.linalg.solve(WW,W.T@x)
    ssr=float(vx@vx); restricted=float(vw@vw)
    F=((restricted-ssr)/k)/(ssr/(n-A.shape[1]))
    result['specifications'].append({'excluded':k,'total_moments':A.shape[1],
        'instrument_rank':rank,'tsls_education':float(b[0]),'first_stage_F_homoskedastic':F})
(root/'output/ak_design_audit.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
