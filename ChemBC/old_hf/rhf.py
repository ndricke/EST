import numpy as np
import re
import hfWel

def genMatrix(filename):
  data = np.genfromtxt(filename)
  #print(data)
  entry = np.size(data,0)
  iter_start = data[0,0]
  nrow = data[-1,0]
  mat = np.zeros((nrow,nrow))
  for i in range(entry):
    [a,b] = data[i,0:2]
    mat[a-iter_start,b-iter_start] = data[i,2]
  return mat

def kListEst(klist):
 enuc = klist[0][0]
 

def genSymMatrix(filename):
  mat = genMatrix(filename)
  return mat + mat.T - np.diag(mat.diagonal())

def superIndex(i,j):
  if i >= j:
    ij = i*(i+1)/2+j
  elif j > i:
    ij = j*(j+1)/2+i
  return ij

def teiIndex(i,j,k,l):
  ij = superIndex(i,j)
  kl = superIndex(k,l)
  ijkl = superIndex(ij,kl)
  return ijkl
"""
def densityMat(C,P,dim,n_occ):
 P_old = np.matrix(P)
 P = np.zeros((dim,dim))
 Carr = np.array(C)
 for u in range(dim):
  for n in range(dim):
   P[u,n] = np.sum(2*Carr[u,0:n_occ]*Carr[n,0:n_occ])
 return P, P_old
"""

def densityMat(C,P,dim,n_occ):
 P_old = np.matrix(P)
 P = np.zeros((dim,dim))
 Carr = C[:,0:n_occ]
 P = 2.0*Carr*Carr.T
 return P, P_old

def hfEnergy(P,H_core,F):
 return np.sum(0.5*np.array(P)*(np.array(H_core)+np.array(F)))

#F = genFock(H,P_mat,teiD)
def genFock(H,P,teiD,dim):
 F = np.array(H)
 for i in range(dim):
  for j in range(dim):
   for k in range(dim):
    for l in range(dim):
     F[i,j] = F[i,j] + P[k,l]*(tei_dense[i,j,k,l] 
              -0.5*tei_dense[i,k,j,l])
 return F

Nelec = 10
nucatt = genSymMatrix('nucatt.dat')
KE = genSymMatrix('ke.dat')
S_orb_overlap = genSymMatrix('s.dat')
S_size = np.size(S_orb_overlap,1)
n_occ = 5
H = KE + nucatt

#reading 2 electron integrals
tei_file = open('eri.dat')
"""
teiD = {}
for line in tei_file:
  lspl = line.split()
  lspl[0:4] = [float(n) for n in lspl[0:4]]
  uvyo = teiIndex(lspl[0],lspl[1],lspl[2],lspl[3])
  teiD[uvyo] = lspl[4]
"""

tei_dense = np.zeros((S_size,S_size,S_size,S_size))
for line in tei_file:
  lspl = line.split()
  lspl[0:4] = [(float(n)-1) for n in lspl[0:4]]
  tei_dense[lspl[0],lspl[1],lspl[2],lspl[3]] = lspl[4]
  tei_dense[lspl[1],lspl[0],lspl[2],lspl[3]] = lspl[4]
  tei_dense[lspl[0],lspl[1],lspl[3],lspl[2]] = lspl[4]
  tei_dense[lspl[1],lspl[0],lspl[3],lspl[2]] = lspl[4]
  tei_dense[lspl[2],lspl[3],lspl[0],lspl[1]] = lspl[4]
  tei_dense[lspl[3],lspl[2],lspl[0],lspl[1]] = lspl[4]
  tei_dense[lspl[2],lspl[3],lspl[1],lspl[0]] = lspl[4]
  tei_dense[lspl[3],lspl[2],lspl[1],lspl[0]] = lspl[4]

eigval, eigvec = np.linalg.eig(S_orb_overlap)
eigvec = np.matrix(eigvec)
inroot_val = eigval**-0.5
ireig_mat = np.matrix(np.diag(inroot_val))
eigvecT = eigvec.T
#print(eigvecT*eigvec)
S_ir = eigvec*ireig_mat*eigvecT
S_irT = S_ir.T
#print(S_ir)

Fock_prime = (S_irT)*H*S_ir
#Fock_matT2 = np.dot(np.transpose(S_ir),np.dot(H,S_ir))
eval_Fock, evec_Fock = np.linalg.eigh(Fock_prime)
evec_FockAO = np.dot(S_ir,evec_Fock)

P_mat = np.zeros((S_size,S_size))
P_mat, P_old = densityMat(evec_FockAO,P_mat,S_size,n_occ)
E_elec = hfEnergy(P_mat,H,H)
#print P_mat
F = genFock(H,P_mat,tei_dense,S_size)
print F
'''
for i in range(100):
 print(E_elec)
 F = genFock(H,P_mat,tei_dense,S_size)
 Fprime = S_irT*F*S_ir
 e, Cprime = np.linalg.eigh(Fprime)
 C = S_ir*Cprime
 P_mat, P_old = densityMat(C,P_mat,S_size,n_occ)
 E_elec = hfEnergy(P_mat,H,F)
 #E = E_elec + E_nuc
 #conv = convTest(E,E_old,P_mat,P_old)
 #if conv == True: break

#E_elec = hfEnergy(P_mat,H,F)
print(E_elec)
'''
hfw = hfWel.HF(H,tei_dense,n_occ,Sov=S_orb_overlap,diis=True)
print hfw.get_energy()
"""
#hint files
sir_hint = np.genfromtxt('sir.dat')
fock_hint = np.genfromtxt('fockhint.dat')
mo_hint = np.genfromtxt('mohint.dat')
d_hint = np.genfromtxt('dhint.dat')
sir_test = sir_hint - S_ir
fock_test = fock_hint - Fock_prime
mo_test = mo_hint - evec_FockAO
#print(np.max(sir_test))
#print(np.max(fock_test))
#print(np.max(mo_test))
mat_hint = np.matrix(mo_hint)
D_test, D_old = densityMat(mo_hint,np.zeros((7,7)),7,5)
D2, D2o = makedensity(mo_hint,np.zeros((7,7)),7,10)
"""











