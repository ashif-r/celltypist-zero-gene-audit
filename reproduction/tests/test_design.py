import pytest,numpy as np
from scipy import sparse
from src.design import retained_mask,detection_counts,partitions,paired_summary,stability_certificate,material_effect

@pytest.mark.parametrize('matrix',[np.array([[0,1,0],[0,0,2]]),sparse.csr_matrix([[0,1,0],[0,0,2]])])
def test_zero_mask(matrix):np.testing.assert_equal(retained_mask(matrix),[False,True,True])
@pytest.mark.parametrize('minimum,expected',[(1,[False,True,True]),(2,[False,True,False]),(3,[False,False,False])])
def test_min_cells(minimum,expected):np.testing.assert_equal(retained_mask(sparse.csr_matrix([[0,1,0],[0,2,3]]),minimum),expected)
def test_explicit_zeros():
 x=sparse.csr_matrix(([0.,1.],([0,0],[0,1])),shape=(1,2));np.testing.assert_equal(retained_mask(x),[False,True])
def test_pool_mask_is_union():
 x=sparse.csr_matrix([[1,0,0],[0,3,0]]);np.testing.assert_equal(retained_mask(x),retained_mask(x[:1])|retained_mask(x[1:]))
def test_computation_chunks_cannot_redefine_mask():
 x=sparse.csr_matrix([[1,0],[0,1]]);mask=retained_mask(x);assert mask.tolist()==[True,True];assert retained_mask(x[:1]).tolist()!=mask.tolist()
@pytest.mark.parametrize('x',[np.array([[np.nan]]),np.array([[np.inf]]),np.array([[-1.]])])
def test_bad_values(x):
 with pytest.raises(ValueError):retained_mask(x)
@pytest.mark.parametrize('minimum',[0,-1,1.5])
def test_bad_threshold(minimum):
 with pytest.raises(ValueError):retained_mask(np.zeros((1,1)),minimum)
def test_empty_query():np.testing.assert_equal(retained_mask(np.zeros((0,3))),[False,False,False])
@pytest.mark.parametrize('size',[1,2,8,20])
def test_partition_membership_and_remainder(size):
 ids=[str(i) for i in range(11)];p=partitions(ids,size,0);flat=np.concatenate(p);assert set(flat)==set(ids) and len(flat)==len(ids);assert all(0<len(v)<=size for v in p);assert len(p[-1])==11%size or len(p[-1])==size
@pytest.mark.parametrize('seed',[0,1,2,3,4])
def test_partition_order_invariant_and_repeated(seed):
 a=partitions(['d','b','c','a'],3,seed);b=partitions(['a','b','c','d'],3,seed)
 assert [v.tolist() for v in a]==[v.tolist() for v in b]
def test_partition_duplicate():
 with pytest.raises(ValueError):partitions(['a','a'],1,0)
def test_invalid_size():
 with pytest.raises(ValueError):partitions(['a'],0,0)
def test_pairs_join_by_id():assert paired_summary(['a','b'],{'b':'T','a':'B'},{'a':'B','b':'B'})=={'eligible':2,'valid':2,'changed':1,'failed':0,'fraction':.5,'complete':True}
@pytest.mark.parametrize('bad',[None,np.nan,''])
def test_failure_not_unchanged(bad):
 s=paired_summary(['a','b'],{'a':'T','b':bad},{'a':'T','b':'T'});assert s['valid']==1 and s['failed']==1 and s['changed']==0 and not s['complete']
def test_unmatched_is_failure():
 s=paired_summary(['a','b'],{'a':'T','b':'T'},{'a':'B'});assert s['changed']==1 and s['fraction']==1 and s['failed']==1
def test_no_pairs_is_missing_not_zero():assert paired_summary(['a'],{}, {})['fraction'] is None
def test_unrelated_extra_ids_not_counted():assert paired_summary(['a'],{'a':'T','x':'B'},{'a':'T','x':'T'})['eligible']==1
def test_duplicate_pair_ids():
 with pytest.raises(ValueError):paired_summary(['a','a'],{}, {})
def test_certificate_strict_boundary():
 a=np.array([[3.,2.,0.],[4.,2.,0.],[2.,2.,0.]])
 np.testing.assert_equal(stability_certificate(a,[0,1,0]),[False,True,False])
def test_certificate_agrees_with_exhaustive_perturbations():
 rng=np.random.default_rng(0);a=rng.normal(size=(100,5));d=rng.normal(size=5);safe=stability_certificate(a,d);assert np.array_equal(a.argmax(1)[safe],(a+d).argmax(1)[safe])
def test_material_boundary_not_null_publication_gate():
 assert material_effect([{'valid':100,'fraction':.01}]);assert not material_effect([{'valid':100,'fraction':.0099}]);assert not material_effect([{'valid':0,'fraction':None}])
