import numpy as np
import re
import pickle

def genMatrix(filename):
  data = np.genfromtxt(filename)
  entry = np.size(data,0)
  iter_start = data[0,0]
  nrow = data[-1,0]
  mat = np.zeros((nrow,nrow))
  for i in range(entry):
    [a,b] = data[i,0:2]
    mat[a-iter_start,b-iter_start] = data[i,2]
  return mat

def genSymMatrix(filename):
  mat = genMatrix(filename)
  return mat + mat.T - np.diag(mat.diagonal())

def densityMat(C,P,dim,n_occ):
 P_old = np.matrix(P)
 P = np.zeros((dim,dim))
 Carr = C[:,0:n_occ]
 P = 2.0*Carr*Carr.T
 return P, P_old

def hfEnergy(P,H_core,F):
 return np.sum(0.5*np.array(P)*(np.array(H_core)+np.array(F)))

def genFock(H,P,teiD,dim):
 F = np.array(H)
 for i in range(dim):
  for j in range(dim):
   for k in range(dim):
    for l in range(dim):
     F[i,j] = F[i,j] + P[k,l]*(tei_dense[i,j,k,l] 
              -0.5*tei_dense[i,k,j,l])
 return F

def fileSplit(filename,delim_list,name_list):
 outfile = None
 file_num = 0
 nameL = len(name_list)
 file = open(filename)
 for line in file:
  if file_num < nameL and delim_list[file_num] in line:
   if outfile != None:
    outfile.close()
    outfile = None
   outfile = open(name_list[file_num],'w')
   file_num += 1
  else:
   if line[1] != '#':
    outfile.write(line)
 outfile.close()
 file.close()

def convTest(E_elec,E_old,P_mat,P_old,del1,del2):
 if abs(E_elec-E_old) < del1 and np.sum((P_mat-P_old)**2)**0.5 < del2:
  return True
 else:
  return False


filename = "LiF.Integrals"
key_list = ['##Nuclear','##Number','##Overlap','##One','##Two']
name_list = ['enuc.dat','neo.dat','s.dat','Hc.dat','eri.dat']
k_list = fileSplit(filename,key_list,name_list)
S_orb_overlap = genMatrix('s.dat')
H = genMatrix('Hc.dat')
S_size = np.size(S_orb_overlap,1)
neo = open('neo.dat')
neo_sp = neo.readline().split()
S_size = int(neo_sp[0])
n_occ = int(neo_sp[1])
Nelec = n_occ*2
nuc_file = open('enuc.dat')
enuc = float(nuc_file.readline())
nuc_file.close()
neo.close()


print('Reading 2 electron integrals')
tei_file = open('eri.dat')

tei_dense = np.zeros((S_size,S_size,S_size,S_size))
for line in tei_file:
  lspl = line.split()
  lspl[0:4] = [(float(n)-1) for n in lspl[0:4]]
  tei_dense[lspl[0],lspl[1],lspl[2],lspl[3]] = lspl[4]

tei_file.close()

print('Calculating eigh of S')
eigval, eigvec = np.linalg.eigh(S_orb_overlap)
eigvec = np.matrix(eigvec)
print('Inv Sq Root')
inroot_val = eigval**-0.5
ireig_mat = np.matrix(np.diag(inroot_val))
eigvecT = eigvec.T
S_ir = eigvec*ireig_mat*eigvecT
S_irT = S_ir.T

Fock_prime = (S_irT)*H*S_ir
eval_Fock, evec_Fock = np.linalg.eigh(Fock_prime)
evec_FockAO = np.dot(S_ir,evec_Fock)

P_mat = np.zeros((S_size,S_size))
P_mat, P_old = densityMat(evec_FockAO,P_mat,S_size,n_occ)
E_elec = hfEnergy(P_mat,H,H)

F = genFock(H,P_mat,tei_dense,S_size)

print('Beginning SCF Procedure')

for i in range(100):
 E_old = E_elec
 F = genFock(H,P_mat,tei_dense,S_size)
 #if i<10:
 # F += np.random.rand(S_size,S_size)*0.00001
 Fprime = S_irT*F*S_ir
 e, Cprime = np.linalg.eigh(Fprime)
 C = S_ir*Cprime
 P_mat, P_old = densityMat(C,P_mat,S_size,n_occ)
 E_elec = hfEnergy(P_mat,H,F)
 print(i,E_elec)
 conv = convTest(E_elec,E_old,P_mat,P_old,0.00001,0.01)
 if conv == True: break

Etot = E_elec + enuc
print("SCF Converged, Total Energy:",Etot)
















