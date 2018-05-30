import numpy
import itertools
import logging


class HF:
    def __init__(s,h,V,m,diis=True,tol=1e-8,Sov=None):
        s.h = h
        s.V = V
        s.m = m
        s.Sov = Sov
        s.diis = diis
        s.tol = tol
        s.C = None
        s.e = None

        if s.Sov != None:
          eigval, eigvec = numpy.linalg.eig(s.Sov)
          inroot_val = eigval**-0.5
          ireig_mat = numpy.diag(inroot_val)
          eigvecT = eigvec.T
          s.S_ir = eigvec.dot(ireig_mat).dot(eigvecT)
          s.S_irT = s.S_ir.T


    def _do_hf(s):
        n = s.h.shape[0]
        if s.Sov == None:
          F,Cocc,C = s.hubbHF()
        else:
          F,Cocc,C = s.regHF()
        emat = Cocc.T.dot(s.h+F).dot(Cocc)
        e = 0
        for i in range(s.m):
            e += emat[i,i]
        s.C = C
        s.e = e

    def regHF(s,max_iter=100):
      n = s.h.shape[0]
      rho,Cocc,C = s.genDenMat(s.h)
      F = s.genFock(rho,phys_notation = False)
      converged = False
      errs = []
      Fs = []
      niter = 0
      for i in range(max_iter):
        print niter
        rho_old = rho.copy()
        rho,Cocc,C = s.genDenMat(F)  
        F = s.genFock(rho, phys_notation = False) 
        err = numpy.linalg.norm(rho-rho_old)
        converged = err<s.tol*n
        if converged: break
        if s.diis:
          errs.append(err)
          Fs.append(F)
          F = HF._next_diis(errs, Fs)
      return F,Cocc,C

    def hubbHF(s,max_iter=100):
        n = s.h.shape[0]
        rho = numpy.eye(n)
        if type(s.V) == (int or float): F = s.h + numpy.diag(numpy.diag(rho)*s.V)
        elif len(s.V.shape)==4: F = s.genFock(rho,phys_notation=True)
        converged = False
        errs = []
        Fs = []
        for i in range(max_iter):
            w,C = numpy.linalg.eigh(F)
            
            # Check for HOMO-LUMO degeneracy
            homo = w[s.m-1]
            lumo = w[s.m]
            if abs(lumo-homo) < 1e-4:
                kennyloggins.warning("DANGER ZONE: HOMO and LUMO are degenerate")

            Cocc = C[:,:s.m]
            rho = Cocc.dot(Cocc.T)

            if type(s.V) == (int or float):
              F = s.h + numpy.diag(numpy.diag(rho)*s.V)
            elif len(s.V.shape)==4:
              F = s.genFock(rho,phys_notation=True)
            err = F.dot(rho) - rho.dot(F)
            converged = numpy.linalg.norm(err)<s.tol*n
            if converged: break
            if s.diis:
                errs.append(err)
                Fs.append(F)
                F = HF._next_diis(errs, Fs)
        return F,Cocc,C



    def genFock(s,rho,phys_notation=True):
      F = s.h.copy()
      ri = range(len(s.h))
      for mu,nu,la,si in itertools.product(ri,ri,ri,ri):
        if phys_notation == True:
          F[mu,nu] += (2*s.V[mu,si,nu,la] - s.V[mu,si,la,nu])*rho[la,si]
        elif phys_notation == False:
          F[mu,nu] += (2*s.V[mu,nu,la,si] - s.V[mu,la,nu,si])*rho[la,si] #Chemist's Notation
      return F

    def genDenMat(s,F):
      Fprime = s.S_irT.dot(F).dot(s.S_ir)
      w,C = numpy.linalg.eigh(Fprime)
      Cocc = numpy.dot(s.S_ir,C)[:,:s.m]
      return Cocc.dot(Cocc.T),Cocc,C

    @staticmethod
    def _next_diis(errs, Fs):
        n = len(errs)
        B = numpy.zeros((n,n))
        for i,j in itertools.product(range(n), range(n)):
            B[i,j] = numpy.dot(errs[i].ravel(), errs[j].ravel())
        A = numpy.zeros((n+1, n+1))
        A[:n,:n] = B
        A[n,:] = -1
        A[:,n] = -1
        A[n,n] = 0
        b = numpy.zeros(n+1)
        b[n] = -1
        x = numpy.linalg.solve(A,b)
        w = x[:n]

        F = numpy.zeros(Fs[0].shape)
        for i in range(n):
            F += w[i] * Fs[i]
        return F

    def get_Cocc(s):
        if s.C is None: s._do_hf()
        return s.C[:,:s.m]

    def get_C(s):
        if s.C is None: s._do_hf()
        return s.C

    def get_energy(s):
        if s.e is None: s._do_hf()
        return s.e

    def get_rho(s):
        C = s.get_Cocc()
        return C.dot(C.T)


if __name__ == "__main__":
    numpy.set_printoptions(precision=3, suppress = True)
    n = 480
    m = 311
    U = 8
    h = numpy.diag(numpy.ones(n-1),1)
    h[0,-1] = 1
    h += h.T
    h*=-1
    for i in range(0,n,2):
        h[i,i] = 1

    hf = HF(h,U,m)
    psi = hf.get_Cocc()
    print psi.dot(psi.T)*2
    print hf.get_energy()
