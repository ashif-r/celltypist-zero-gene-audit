import numpy as np,pytest
from scipy import sparse
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from celltypist.models import Model
from src.inference import api_scores,direct_scores,query,restore_zeros,compare

def tiny_model():
 c=LogisticRegression();c.classes_=np.array(['a','b','c']);c.coef_=np.array([[2.,1.,0.],[-1.,0.,1.],[0.,-1.,-1.]]);c.intercept_=np.array([0.,1.,0.]);c.n_features_in_=3;c.features=np.array(['g0','g1','g2'])
 s=StandardScaler().fit([[0.,0.,1.],[2.,2.,3.],[1.,1.,2.]])
 return Model(c,s,{})

def test_native_query_and_physical_restore(tmp_path):
 m=tiny_model();x=sparse.csr_matrix([[0,1,1],[0,2,0],[0,3,1],[0,4,0]],dtype=float);genes=np.array(['g0','g1','g2']);ids=np.array(['c1','c2','c3','c4'])
 A=api_scores(x,genes,ids,m,chunk=2);compare(direct_scores(x,genes,m,chunk=3),A,'tiny')
 before=m.classifier.coef_.copy();v=query(x,genes,ids,m,A,tmp_path/'query',do_D=True)
 assert v['removed_genes']==1 and v['D_status']=='passed';assert v['B_max_score_residual']<1e-12
 np.testing.assert_equal(before,m.classifier.coef_)
 saved=np.load(tmp_path/'query/observed.npz');np.testing.assert_allclose(saved['C'],A,atol=1e-12)

def test_zero_restore_keeps_order():
 x=sparse.csr_matrix([[1.,0.,2.],[0.,0.,4.]]);mask=np.array([True,False,True]);assert (restore_zeros(x[:,mask],mask,x.shape)!=x).nnz==0

def test_broken_control_raises():
 with pytest.raises(AssertionError):compare(np.ones((2,3)),np.zeros((2,3)),'bad')

def test_api_duplicate_genes_refused():
 with pytest.raises(ValueError,match='Nonunique'):api_scores(sparse.csr_matrix([[1.,1.]]),np.array(['g1','g1']),np.array(['c']),tiny_model())

def test_model_prefilter_equals_full_native_api():
 from anndata import AnnData
 from celltypist import annotate
 import copy,pandas as pd
 m=tiny_model();genes=np.array(['ignored','g2','g0','g1']);ids=np.array(['c1','c2'])
 x=sparse.csr_matrix([[9.,1.,0.,2.],[3.,0.,0.,4.]])
 full=annotate(AnnData(x,obs=pd.DataFrame(index=ids),var=pd.DataFrame(index=genes)),model=copy.deepcopy(m),majority_voting=False)
 np.testing.assert_array_equal(api_scores(x,genes,ids,m),full.decision_matrix.to_numpy())
 np.testing.assert_array_equal(full.predicted_labels.iloc[:,0].astype(str),m.cell_types[api_scores(x,genes,ids,m).argmax(1)])

def test_invalid_ignored_gene_cannot_bypass_contract():
 with pytest.raises(ValueError,match='before model-feature'):
  api_scores(sparse.csr_matrix([[10.,1.]]),np.array(['ignored','g1']),np.array(['c']),tiny_model())
