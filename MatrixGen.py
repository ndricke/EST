import numpy as np
import scipy as sc
#from scipy.stats import ortho_group


def hubbH(n,U):
    h = np.diag(np.ones(n-1),1)
    h[0,-1] = -1
    h += h.T
    V = np.zeros((n,n,n,n))
    for i in range(n):
        V[i,i,i,i] = U
	return h, V


def antisymMat(n,scale=1.):
    rand_mat = np.random.rand(n,n)*scale
    lhalf = np.tril(rand_mat)
    return lhalf.T - lhalf

def randUnitary(n,method='ortho'):
#    if method == 'ortho':
#        return sc.stats.ortho_group.rvs(dim=n)
#    elif method == 'exp':
        D = antisymMat(n)
        return sc.linalg.expm(D)





