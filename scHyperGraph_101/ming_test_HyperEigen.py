import scipy.sparse as ss
import scipy.sparse.linalg as ssalg
from scipy.sparse import csr_matrix
import numpy as np

# generate a random matrix
#mat=np.random.random((3, 3))
mat=[[0.1,0.9],[1,0]]
inc_mat=csr_matrix(mat)
inc_mat.todense()

# get degree matrix
sqr_mat = inc_mat.power(2)
sqr_sum = ss.csr_matrix.sum(sqr_mat, 1)
deg_mat = ss.csr_matrix.multiply(ss.identity(sqr_sum.shape[0]), sqr_sum)
deg_mat = csr_matrix(deg_mat)
deg_mat.todense()
#matrix([[0.82, 0.  ],
#        [0.  , 1.  ]])

# get laplacian matrix
#lap_mat = get_hyperedge_normalised_laplacian_py(inc_mat, deg_mat)
    
for i in range(deg_mat.shape[0]):
    if deg_mat[i,i] != 0:
        deg_mat[i,i] = deg_mat[i,i]**-1
    
deg_mat.todense()
inc_T = inc_mat.T
part_lap = inc_T.dot(deg_mat)
lap_mat = part_lap.dot(inc_mat)
    
lap_mat.todense()
#matrix([[1.01219512, 0.1097561 ],
#        [0.1097561 , 0.98780488]])

# eigenvalues        
eigen_val = ssalg.eigs(lap_mat, k = 1)
eigen_val = eigen_val[0] / inc_mat.shape[0]


