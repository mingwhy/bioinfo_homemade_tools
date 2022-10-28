
import numpy as np
import nimfa
import datetime
import scipy.sparse as spr

data=np.load('hvg2k_mat.npy',allow_pickle =True)
Vm= data.all() #extract sparse matrix from array[sparse matrix ... ]

start= datetime.datetime.now()
nmf = nimfa.Nmf(Vm, rank=50, max_iter=100, update='euclidean', objective='fro')
nmf_fit = nmf()
end= datetime.datetime.now()
print(end-start)

start= datetime.datetime.now()
lsnmf = nimfa.Nmf(Vm, seed='random_vcol', rank=50, max_iter=100,update='euclidean', objective='fro')
lsnmf_fit = lsnmf()
end= datetime.datetime.now()
print(start)
print(end)
print(end-start) # 0:37:48.012543

print('Rss: %5.4f' % lsnmf_fit.fit.rss())
print('Evar: %5.4f' % lsnmf_fit.fit.evar())
print('K-L divergence: %5.4f' % lsnmf_fit.distance(metric='kl'))
print('Sparseness, W: %5.4f, H: %5.4f' % lsnmf_fit.fit.sparseness())

print('Iterations: %d' % lsnmf_fit.n_iter)
print('Target estimate:\n%s' % np.dot(W, H))

W = lsnmf_fit.basis()
print('Basis matrix:\n%s' % W)

H = lsnmf_fit.coef()
print('Mixture matrix:\n%s' % H)

