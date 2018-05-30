import numpy as np
import re
import pickle

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
 data2 = klist[1].split()
 dim = data2[0]
 nelec = data2[1]
 sfile = open('s.dat','w')
 hcfile = open('Hc.dat','w')
 teifile = open('eri.dat','w')
 pickle.dump(klist[2],sfile)
 pickle.dump(klist[3],hcfile)
 pickle.dump(klist[4],teifile)
 return enuc, dim, nelec

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

def densityMat(C,P,dim,n_occ):
 P_old = np.matrix(P)
 P = np.zeros((dim,dim))
 Carr = np.array(C)
 for u in range(dim):
  for n in range(dim):
   P[u,n] = np.sum(2*Carr[u,0:n_occ]*Carr[n,0:n_occ])
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

def fileSplit(filename,key_list,end_list):
 file = open(filename)
 k = 0
 klist = []
 for item in key_list: klist.append([])

 parsing = False
 for line in file:
  if line in end_list:
   parsing == False
   k += 1
  elif line in key_list:
   parsing == True
  if parsing == True:
   if line[0] == '##':
    continue
   else:
    klist[k].append(line.strip())
 return klist

"""
def ESTread(filename):
 file = open(filename)
 for line in file:
  if line[0:1] == '##':
   enuc = float(line.next())
  if line[0:1] == '##':
   line.next()
 return enuc

Nelec = 10
nucatt = genSymMatrix('nucatt.dat')
n_occ = 5
H = KE + nucatt
"""
filename = "H2O.Integrals"
key_list = ['##Nuclear Repulsion Energy','##Number of Basis Functions and Occupied Orbitals','##Overlap Matrix','##One Electron Integrals','##Two Electron Integrals (ij|kl)']
end_list = ['##Number of Basis Functions and Occupied Orbitals','##Overlap Matrix','##One Electron Integrals','##Two Electron Integrals (ij|kl)']
k_list = fileSplit(filename,key_list,end_list)
enuc, S_size, Nelec = kListEst(k_list)
S_orb_overlap = genSymMatrix('s.dat')
H = genSymMatrix('Hc.dat')
S_size = np.size(S_orb_overlap,1)
n_occ = Nelec/2

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

F = genFock(H,P_mat,tei_dense,S_size)


for i in range(100):
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











