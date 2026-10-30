from types import SimpleNamespace
import numpy as np
from scipy import sparse
from src.gate_core import identifier_keep,normalize_counts,scores,zero_shift

def mock(center=True):
    return SimpleNamespace(features=np.array(['g1','g2']),cell_types=np.array(['a','b','c']),scaler=SimpleNamespace(with_mean=center,mean_=np.array([2.,1.]),scale_=np.array([2.,1.])),classifier=SimpleNamespace(coef_=np.array([[2.,-1.],[1.,1.],[-3.,2.]]),intercept_=np.array([0.,2.,-1.])))
def test_zero_deletion_and_all_classes():
    m=mock();x=sparse.csr_matrix([[0.,3.],[0.,1.]])
    a=scores(x,m.features,m);b=scores(x[:,1:],m.features[1:],m)
    np.testing.assert_allclose(a+zero_shift(['g1'],m),b,atol=1e-12)
def test_uncentered_zero():
    m=mock(False);np.testing.assert_equal(zero_shift(['g1'],m),[0,0,0])
def test_duplicate_denominator_order_and_real_hyphen():
    names=np.array(['dup','dup','HLA-DRA','g']); x=sparse.csr_matrix([[0,5,3,2]])
    keep=identifier_keep(names);np.testing.assert_equal(keep,[False,False,True,True])
    y,t=normalize_counts(x);assert t[0]==10
    np.testing.assert_allclose(np.expm1(y[:,keep].toarray()).sum(),5000)
    perm=np.array([2,1,3,0])
    assert set(names[keep])==set(names[perm][identifier_keep(names[perm])])
def test_common_shift_cancellation_and_nonrunner():
    a=np.array([[3.,2.,0.]])
    assert (a+np.ones(3)*9).argmax()==a.argmax()
    np.testing.assert_equal(np.array([2.,-1.,0.])+np.array([-2.,1.,0.]),[0,0,0])
    assert (a+np.array([0.,0.,4.])).argmax()==2
    margin=1.;spread=np.ptp([0.,0.,4.]);assert margin<spread

def test_D_and_clipping():
    m=mock(); x=sparse.csr_matrix([[40.,3.],[1.,0.]])
    a=scores(x,m.features,m);b=scores(x[:,1:],m.features[1:],m)
    scaled=np.minimum((x[:,0].toarray().ravel()-2)/2,10)
    np.testing.assert_allclose(a-scaled[:,None]*m.classifier.coef_[:,0],b)


def test_explicit_stored_zero_not_detected():
    x=sparse.csr_matrix(([0.,2.],([0,0],[0,1])),shape=(1,2))
    y,_=normalize_counts(x)
    assert y.getnnz(axis=0).tolist()==[0,1]
