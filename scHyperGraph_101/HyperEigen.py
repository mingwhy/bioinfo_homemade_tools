import scipy.sparse as ss
import scipy.sparse.linalg as ssalg
from scipy.sparse import csr_matrix
import numpy as np

#Function to return the hypergraph degree matrix from its incidence matrix - sparse implementation
def get_degree_matrix_py(inc_mat):
    '''
    Parameters
    ----------
    inc_mat : csr_matrix
        Incidence matrix of hypergraph

    Returns
    -------
    deg_mat : csr_matrix
        Degree matrix of hypergraph

    '''
    if ss.isspmatrix(inc_mat) == False:
        inc_mat = csr_matrix(inc_mat)
    sqr_mat = inc_mat.power(2)
    sqr_sum = ss.csr_matrix.sum(sqr_mat, 1)
    deg_mat = ss.csr_matrix.multiply(ss.identity(sqr_sum.shape[0]), sqr_sum)
    deg_mat = csr_matrix(deg_mat)
    return deg_mat

#Function to return the hypergraph hyperedge normalise laplacian matrix from its incidence matrix - sparse implementation
def get_hyperedge_normalised_laplacian_py(inc_mat, deg_mat = None):
    '''
    Parameters
    ----------
    inc_mat : csr_matrix
        Incidence matrix of hypergraph
    deg_mat : csr_matrix, optional
        Degree matrix of hypergraph. If None then it is calculated from the incidence matrix. The default is None.

    Returns
    -------
    lap_mat : csr_matrix
        The hyperedge normalised laplacian matrix of the hypergraph

    '''
    if not ss.isspmatrix(inc_mat):
        inc_mat = csr_matrix(inc_mat)
    if deg_mat is None:
        deg_mat = get_degree_matrix_py(inc_mat)
    elif not ss.isspmatrix_csr(deg_mat):
        deg_mat = csr_matrix(deg_mat)
    else:
        deg_mat = deg_mat
    for i in range(deg_mat.shape[0]):
        if deg_mat[i,i] != 0:
            deg_mat[i,i] = deg_mat[i,i]**-1
    inc_T = inc_mat.T
    part_lap = inc_T.dot(deg_mat)
    lap_mat = part_lap.dot(inc_mat)
    return lap_mat

#Function to return the largest eigenvalue of the hyperedge normalise laplacian matrix from its incidence matrix - sparse implementation
def get_largest_eigenvalue_py(inc_mat, lap_mat = None, deg_mat = None):
    '''
    Parameters
    ----------
    inc_mat : csr_matrix
        Incidence matrix of hypergraph
    lap_mat : csr_matrix optional
        The hyperedge normalised laplacian of the hypergraph. If None then it is calculated from incidence matrix. The default is None.
    deg_mat : csr_matrix, optional
        The degree matrix of the hypergraph. If needed, is used to speed up finding the hyperedge normalised laplacian matrix. The default is None.

    Returns
    -------
    eigen_val : complex
        The largest eigenvalue of the hyperedge normalised laplacian matrix.

    '''
    if not ss.isspmatrix(inc_mat):
        inc_mat = csr_matrix(inc_mat)
    if lap_mat is None:
        lap_mat = get_hyperedge_normalised_laplacian_py(inc_mat, deg_mat)
    elif ss.isspmatrix_csr(lap_mat):
        lap_mat = csr_matrix(lap_mat)
    lap_mat = get_hyperedge_normalised_laplacian_py(inc_mat)
    eigen_val = ssalg.eigs(lap_mat, k = 1)
    eigen_val = eigen_val[0] / inc_mat.shape[0]
    return eigen_val


count_mat = r.CountsMatrix
large_val = get_largest_eigenvalue_py(inc_mat = count_mat)
