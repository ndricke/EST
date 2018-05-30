import numpy as np
import itertools
import copy
import re


def chem2phys(C):
    n = C.shape[0] #assuming that C is np.array
    P = np.zeros((n,n,n,n))
    nz = range(n)
    for i,j,k,l in itertools.product(nz,nz,nz,nz):
        P[i,j,k,l] = C[i,k,j,l]
    return P

def AssertPhysSym(C): #Doesn't care about numerical accuracy; this should be true if strictly enforced
    n = C.shape[0]
    nz = range(n)
    try:
        for i,j,k,l in itertools.product(nz,nz,nz,nz):
            assert C[i,j,k,l]==C[j,i,l,k]
            assert C[i,j,k,l]==C[k,l,i,j]
            assert C[i,j,k,l]==C[l,k,j,i]
            assert C[i,j,k,l]==C[k,j,i,l]
            assert C[i,j,k,l]==C[l,i,j,k]
            assert C[i,j,k,l]==C[i,l,k,j]
            assert C[i,j,k,l]==C[j,k,l,i]
        return True
    except:
        return False
    
def AssertChemSym(C): #Doesn't care about numerical accuracy; this should be true if strictly enforced
    n = C.shape[0]
    nz = range(n)
    try:
        for i,j,k,l in itertools.product(nz,nz,nz,nz):
            assert C[i,j,k,l]==C[j,i,l,k]
            assert C[i,j,k,l]==C[k,l,i,j]
            assert C[i,j,k,l]==C[l,k,j,i]
            assert C[i,j,k,l]==C[j,i,k,l]
            assert C[i,j,k,l]==C[l,k,i,j]
            assert C[i,j,k,l]==C[i,j,l,k]
            assert C[i,j,k,l]==C[k,l,j,i]
        return True
    except:
        return False

def EnforceChemSym(V,p_zero=10.0**-8): #Blithely enforce symmetry of real 2ei to match sym of chemical notation
    n = V.shape[0]
    C = copy.deepcopy(V) #so that we aren't modifying V in-place
    nz = range(n)
    for i,j,k,l in itertools.product(nz,nz,nz,nz):
        if np.abs(C[i,j,k,l]) > p_zero: #if the entry isn't empty, copy it to all the other ones
            C[j,i,l,k]=C[i,j,k,l]
            C[k,l,i,j]=C[i,j,k,l]
            C[l,k,j,i]=C[i,j,k,l]
            C[j,i,k,l]=C[i,j,k,l]
            C[l,k,i,j]=C[i,j,k,l]
            C[i,j,l,k]=C[i,j,k,l]
            C[k,l,j,i]=C[i,j,k,l]
    return C

def EnforcePhysSym(V,p_zero=10.0**-8): #Blithely enforce symmetry of real 2ei to match sym of chemical notation
    n = V.shape[0]
    C = copy.deepcopy(V) #so that we aren't modifying V in-place
    nz = range(n)
    for i,j,k,l in itertools.product(nz,nz,nz,nz):
        if np.abs(C[i,j,k,l]) > p_zero: #if the entry isn't empty, copy it to all the other ones
            C[j,i,l,k]=C[i,j,k,l]
            C[k,l,i,j]=C[i,j,k,l]
            C[l,k,j,i]=C[i,j,k,l]
            C[k,j,i,l]=C[i,j,k,l]
            C[l,i,j,k]=C[i,j,k,l]
            C[i,l,k,j]=C[i,j,k,l]
            C[j,k,l,i]=C[i,j,k,l]
    return C

def CheckPhysSym(C):
    n = C.shape[0]
    nz = range(n)
    tracker = 0.0
    for i,j,k,l in itertools.product(nz,nz,nz,nz):
        tracker += np.abs(C[i,j,k,l]-C[j,i,l,k])
        tracker += np.abs(C[i,j,k,l]-C[k,l,i,j])
        tracker += np.abs(C[i,j,k,l]-C[l,k,j,i])
        tracker += np.abs(C[i,j,k,l]-C[k,j,i,l])
        tracker += np.abs(C[i,j,k,l]-C[l,i,j,k])
        tracker += np.abs(C[i,j,k,l]-C[i,l,k,j])
        tracker += np.abs(C[i,j,k,l]-C[j,k,l,i])
    return tracker
    
def CheckChemSym(C):
    n = C.shape[0]
    nz = range(n)
    tracker = 0.0
    for i,j,k,l in itertools.product(nz,nz,nz,nz):
        tracker += np.abs(C[i,j,k,l]-C[j,i,l,k])
        tracker += np.abs(C[i,j,k,l]-C[k,l,i,j])
        tracker += np.abs(C[i,j,k,l]-C[l,k,j,i])
        tracker += np.abs(C[i,j,k,l]-C[j,i,k,l])
        tracker += np.abs(C[i,j,k,l]-C[l,k,i,j])
        tracker += np.abs(C[i,j,k,l]-C[i,j,l,k])
        tracker += np.abs(C[i,j,k,l]-C[k,l,j,i])
    return tracker

def printSymInd(V,rep=5):
    n = V.shape[0]
    for r in range(rep):
        i = np.random.randint(n);j = np.random.randint(n);k = np.random.randint(n);l = np.random.randint(n);
        print "Random sample",r
        print i,j,k,l
        print V[i,j,k,l],V[j,i,l,k],V[k,l,i,j],V[l,k,j,i]
        print V[j,i,k,l],V[l,k,i,j],V[i,j,l,k],V[k,l,j,i]

def hfEnergy(P,h,F):
    Emat = np.dot(P,h+F)
    print np.shape(Emat)
    return np.trace(Emat)

#I think this assumes Chemist's Notation (re. the exchange terms)
def genFock(P,h,V):
    n = h.shape[0]
    F = np.array(h)
    nz = range(n)
    for i,j,k,l in itertools.product(nz,nz,nz,nz):
        F[i,j] += P[k,l]*(2.0*V[i,j,k,l] - V[i,k,j,l])
#        F[i,j] += P[k,l]*(V[i,j,k,l] - 0.5*V[i,k,j,l])
    return F

def pdmEnergy(P,h,V):
    n = h.shape[0]
    nz = range(n)
    e1 = 2.0*np.trace(np.dot(h,P))
    e2 = 0.0
    for i,j,k,l in itertools.product(nz,nz,nz,nz):
        e2 += (2.0*V[i,j,k,l]-V[i,j,l,k])*P[k,i]*P[j,l]
    return e1, e2

def matrixGrab(f,n):
    M = np.zeros((n,n))
    for line in f: 
        if not re.match("^[0-9 .-]+$", line): break #means we've hit text, so matrix complete
        spline = line.split()
        if re.match("^[0-9 ]+$", line): #no signs or decimals, so it's an index line
            col_inds = [int(i)-1 for i in spline]
        else: #it's got numbers, and either negatives or decimals, so it's float entries for the matrix
            row_indx = int(spline[0])-1
            row_elements = [float(i) for i in spline[1:]]
            M[row_indx,col_inds] = row_elements
    return M

def tensor4grab(f,n):
    f.next() #skip the blank line that comes before the 2ei
    V = np.zeros((n,n,n,n))
    for line in f:
        if line.isspace(): break
        spline = line.split()
        ix = tuple([int(i)-1 for i in spline[:4]])
        V[ix] = float(spline[-1])
    return V

def toDense(indices,ei2,notation):
    n = np.max(indices)
    V = np.zeros((n,n,n,n))
    for i in range(indices.shape[0]):
        ix = tuple(j-1 for j in indices[i,:])
        V[ix] = ei2[i]
    if notation == "chem":
        V = EnforceChemSym(V)
    if notation == "phys":
        V = EnforcePhysSym(V)
    return V

def saveHV(outfile,n,h,V,tol):
    with open(outfile,'w') as f:
        f.write(str(n)+"\n")
        rn = range(n)
        for i,j in itertools.product(rn,rn):
#            if abs(h[i,j]) > tol:
                txt = ' '.join([str(a) for a in [i,j]])
                txt += ' %.17f' % h[i,j]
                f.write(txt+"\n")
        for i,j,k,l in itertools.product(rn,rn,rn,rn):
#            if abs(V[i,j,k,l]) > tol:
                txt = ' '.join([str(a) for a in [i,j,k,l,]])
                txt += ' %.17f' % V[i,j,k,l]
                f.write(txt+"\n")




