"""Reconstruct the archived transfer artifact from explicit planar skein states."""
from __future__ import annotations
import argparse, hashlib, json, time
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from boundary_algebra import Y, raw_realization, BASIS_WORDS, STATES


def payload(m):
    return [[str(x) for x in row] for row in m.to_list()]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=Path("_build/fixed_closure_proof"))
    args=parser.parse_args()
    if args.out.exists():
        raise FileExistsError("Use a new output directory to preserve earlier evidence")
    start=time.perf_counter()
    L, A, B, R=raw_realization()
    matrices=[DomainMatrix.from_Matrix(m).to_field() for m in (L,A,B,R)]
    F=matrices[1].domain
    L,A,B,R=[m.convert_to(F).to_dense() for m in matrices]
    print('Explicit realization ready',F,time.perf_counter()-start,flush=True)
    lefts={};rights={}
    for w in BASIS_WORDS:
        left=L;right=R
        for ch in w: left=left.matmul(A if ch=='A' else B)
        for ch in reversed(w): right=(A if ch=='A' else B).matmul(right)
        lefts[w]=left;rights[w]=right
    def make(mid):
        rows=[]
        for u in BASIS_WORDS:
            row=[]
            for v in BASIS_WORDS:
                middle=rights[v] if mid is None else mid.matmul(rights[v])
                row.append(lefts[u].matmul(middle).to_list()[0][0])
            rows.append(row)
        return DomainMatrix.from_list(rows,F)
    H=make(None);HA=make(A);HB=make(B)
    print('Hankel blocks assembled',time.perf_counter()-start,flush=True)
    C=DomainMatrix.from_list([[rights[w].to_list()[i][0] for w in BASIS_WORDS] for i in range(15)],F)
    O=DomainMatrix.from_list([lefts[w].to_list()[0] for w in BASIS_WORDS],F)
    assert O.matmul(C).sub(H).is_zero_matrix
    assert O.matmul(A.matmul(C)).sub(HA).is_zero_matrix
    assert O.matmul(B.matmul(C)).sub(HB).is_zero_matrix
    raw_blocks={'variable':'Y','basis':BASIS_WORDS,
                **{name:[[str(sp.expand(value)) for value in row] for row in dm.to_Matrix().tolist()]
                   for name,dm in [('H',H),('H_A',HA),('H_B',HB)]}}
    block_bytes=(json.dumps(raw_blocks,indent=2,sort_keys=True)+'\n').encode()
    git_blob=hashlib.sha1(b'blob '+str(len(block_bytes)).encode()+b'\0'+block_bytes).hexdigest()
    expected='ed8d20c2c97d7cdcf62d4e87cae153e658d71224'
    det=H.to_Matrix().subs(Y,2).det(method='domain-ge')
    assert det!=0
    identity=DomainMatrix.eye((15,15),F).to_dense()
    for op in (A,B):
        cubic=identity
        for exponent in (2,-2,-4):
            cubic=cubic.matmul(op.sub(identity.scalarmul(F.convert(Y**exponent))))
        assert cubic.is_zero_matrix
    assert not A.matmul(B).sub(B.matmul(A)).is_zero_matrix
    result={'independent_basis_dimension':15,'basis_states':STATES,
            'H_equals_O_C':True,'H_A_equals_O_MA_C':True,'H_B_equals_O_MB_C':True,
            'determinant_H_at_2':str(det),'rank_15_verified':True,
            'cubic_A_identity':True,'cubic_B_identity':True,'operators_noncommute':True,
            'generated_hankel_blob_sha':git_blob,'expected_archive_hankel_blob_sha':expected,
            'all_three_hankel_blocks_byte_identical_to_archive':git_blob==expected,
            'elapsed_seconds':time.perf_counter()-start,'input_training_coefficients_used':False,
            'scope':'Fixed blackboard-framed three-strand closure only; not a volume/spine certificate.',
            'proof':'The explicit 15-state skein realization and invertible H imply the all-word Hankel realization by change of basis.'}
    out=args.out;out.mkdir(parents=True,exist_ok=False)
    (out/'certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'independent_hankel_input_matrices.json').write_bytes(block_bytes)
    (out/'explicit_planar_realization.json').write_text(json.dumps({'field':str(F),'variable':'Y','basis_states':STATES,
        'M_A':payload(A),'M_B':payload(B),'ell':payload(L),'r':payload(R),'change_of_basis_C':payload(C),
        'observability_O':payload(O)},indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    if git_blob!=expected:
        raise AssertionError('Hankel block identity mismatch')

if __name__=='__main__':main()
