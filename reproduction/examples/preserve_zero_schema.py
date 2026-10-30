"""Keep evidence when deleting observed-zero columns, then restore that exact schema.

This utility records zeros in a supplied nonnegative measured matrix. It does not
infer missing/unmeasured genes from a model or turn an unknown absence into zero.
Normalize first if required; neither function renormalizes or modifies values.
"""
from dataclasses import dataclass
import hashlib
import numpy as np
from scipy import sparse

@dataclass(frozen=True)
class RecordedZeroColumns:
    genes:tuple
    cell_ids:tuple
    kept_columns:tuple
    retained_values_sha256:str

def _digest(x):
    x=sparse.csr_matrix(x).copy();x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
    h=hashlib.sha256()
    for a in [np.asarray(x.shape,dtype=np.int64),x.indptr.astype(np.int64),x.indices.astype(np.int64),x.data.astype(np.float64)]:h.update(a.tobytes())
    return h.hexdigest()

def filter_and_record(x,genes,cell_ids):
    x=sparse.csr_matrix(x);genes=tuple(genes);ids=tuple(cell_ids)
    if x.shape!=(len(ids),len(genes)) or len(set(genes))!=len(genes) or len(set(ids))!=len(ids):raise ValueError('Shape or identifier ambiguity')
    if not np.isfinite(x.data).all() or (x.data<0).any():raise ValueError('Nonnegative finite measured values required')
    keep=np.flatnonzero(np.asarray((x>0).sum(axis=0)).ravel()>0);filtered=x[:,keep].tocsr()
    return filtered,tuple(genes[i] for i in keep),RecordedZeroColumns(genes,ids,tuple(keep.tolist()),_digest(filtered))

def restore_recorded_zeros(filtered,genes,cell_ids,record):
    expected=tuple(record.genes[i] for i in record.kept_columns)
    if tuple(genes)!=expected or tuple(cell_ids)!=record.cell_ids or _digest(filtered)!=record.retained_values_sha256:raise ValueError('Cell/feature alignment or expression changed')
    c=sparse.coo_matrix(filtered);keep=np.asarray(record.kept_columns,dtype=int)
    return sparse.csr_matrix((c.data,(c.row,keep[c.col])),shape=(len(record.cell_ids),len(record.genes)))
