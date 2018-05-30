import numpy as np



class HF:

  def __init__(s,S_data,T_data,V_data,Nuc_data,eri_data):
    s.S = genMatrix(S_data)
    s.T = genMatrix(T_data)
    s.V = genMatrix(V_data)
    s.n_occ = n_elec/2
    s.nuc_rep = genMatrix(Nuc_data)
    s.h = T + V
    s.tei = genTei(eri_data)
    s.S_ir = np.zeros((s.S_size,s.S_size))
    s.S_irT = np.zeros((s.S_size,s.S_size))
    s.tolerance = [10^-7,10^-5]
    
  def genTei(s,eri_data):
    tei_file = open('eri.dat')
    tei_dense = np.zeros((s.S_size,s.S_size,s.S_size,s.S_size))
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
    return tei_dense

  def rHF(s):
    s.S_ir = invRootMat(s.S)
    s.S_irT = s.S_ir.T
    Fock_prime = (s.S_irT)*H*s.S_ir
    eval_Fock, evec_Fock = np.linalg.eigh(Fock_prime)
    evec_FockAO = np.dot(S_ir,evec_Fock)
    s.P_mat = np.zeros((s.S_size,s.S_size))
    s.P_mat, s.P_old = densityMat(evec_FockAO,s.P_mat,s.S_size,s.n_occ)
    
    for i in range(iter_max):
      HFcycle();
      if convTest(E_elec,E_old,s.tolerance[1],s.P_mat,s.P_old,s.tolerance[1]) == True:
        break

  def convTest(s,E,E_old,dE_tol,P,P_old,dP_tol):
   dP = np.linalg.norm(P - P_old)
   if abs(E-E_old) < dE_tol and dP < dP_tol:
    return True
   else:
    return False

  def genMatrix(s, filename):
    data = np.genfromtxt(filename)
    entry = np.size(data,0)
    iter_start = data[0,0]
    nrow = data[-1,0]
    mat = np.zeros((nrow,nrow))
    for i in range(entry):
      [a,b] = data[i,0:2]
      mat[a-iter_start,b-iter_start] = data[i,2]
    return mat

  def genSymMatrix(s, filename):
    mat = genMatrix(filename)
    return mat + mat.T - np.diag(mat.diagonal())

  def densityMat(s, C,P,dim,n_occ):
    P_old = np.matrix(P)
    P = np.zeros((dim,dim))
    Carr = C[:,0:n_occ]
    P = 2.0*Carr*Carr.T
      return P, P_old

  def hfEnergy(s, P,H_core,F):
    return np.sum(0.5*np.array(P)*(np.array(H_core)+np.array(F)))

  def genFock(s, H,P,tei_dense,dim):
   F = np.array(H)
   for i in range(dim):
    for j in range(dim):
     for k in range(dim):
      for l in range(dim):
       F[i,j] = F[i,j] + P[k,l]*(tei_dense[i,j,k,l]
              -0.5*tei_dense[i,k,j,l])
   return F  

  def invRootMat(s, mat):
    eigval, eigvec = np.linalg.eig(mat)
    eigvec = np.matrix(eigvec)
    inroot_val = eigval**-0.5
    ireig_mat = np.matrix(np.diag(inroot_val))
    eigvecT = eigvec.T
    ir_mat = eigvec*ireig_mat*eigvecT
    return ir_mat

  def HFcycle(s):
    F = genFock(s.h,s.P_mat,s.tei,s.S_size)
    Fprime = s.S_irT*F*s.S_ir
    eps, Cprime = np.linalg.eigh(Fprime)
    C = s.S_ir*Cprime
    s.P_mat, s.P_old = densityMat(C,s.P_mat,s.S_size,s.n_occ)
  
  def mp2(s):
    MO_int = np.zeros((s.S_size,s.S_size,s.S_size,s.S_size))
    for p in range(s.S_size):
     for q in range(s.S_size):
      for r in range(s.S_size):
       for s in range(s.S_size):
        for i in range(s.S_size):
         for j in range(s.S_size):
          for k in range(s.S_size):
           for l in range(s.S_size):
            MO_int[p,q,r,s] += s.P_mat[i,p]*s.P_mat[j,q]*s.P_mat[k,r]*s.P_mat[l,s]*s.tei[i,j,k,l]

    E_mp2 = 0.0
    for i in range(s.n_occ):
     for j in range(s.n_occ):
      for a in np.arange(s.n_occ,s.S_size,1):
       for b in np.arange(s.n_occ,s.S_size,1):
        E_mp2 += MO_int[i,a,j,b]*(2*MO_int[i,a,j,b]-MO_int[i,b,j,a])/(s.eps[i]+s.eps[j]-s.eps[a]-s.eps[b])
    return E_mp2


if __name__ == "__main__":
 eri_data = "eri.dat"
 S_data = "s.dat"
 T_data = "t.dat"
 V_data = "v.dat"
 Nuc_data = "enuc.dat"
 gen_hf = HF(S_data,T_data,V_data,Nuc_data,eri_data)
 gen_hf.rHF()
 e_mp2 = gen_hf.mp2()
 print(gen_hf.elec_energy)
 print(e_mp2)
























