import pytest,numpy as np
from scipy import sparse
from examples.preserve_zero_schema import filter_and_record,restore_recorded_zeros

def test_restores_recorded_zero_values_only():
 x=sparse.csr_matrix([[1.,0.,2.,0.],[0.,0.,3.,0.]])
 f,g,r=filter_and_record(x,['a','b','c','d'],['x','y']);assert g==('a','c');z=restore_recorded_zeros(f,g,['x','y'],r);assert (z!=x).nnz==0
@pytest.mark.parametrize('change',['cells','genes','values'])
def test_refuses_changed_alignment_or_values(change):
 x=sparse.csr_matrix([[1.,0.,2.],[0.,0.,3.]]);f,g,r=filter_and_record(x,['a','b','c'],['x','y'])
 ids=['y','x'] if change=='cells' else ['x','y'];gn=g[::-1] if change=='genes' else g
 if change=='values':f=f*2
 with pytest.raises(ValueError):restore_recorded_zeros(f,gn,ids,r)
def test_all_zero_schema_can_restore_empty_columns():
 x=sparse.csr_matrix((2,3));f,g,r=filter_and_record(x,['a','b','c'],['x','y']);assert f.shape==(2,0);assert restore_recorded_zeros(f,g,['x','y'],r).shape==(2,3)
