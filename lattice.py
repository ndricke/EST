import sys
import numpy
import os.path
import logging
kennyloggins = logging.getLogger('lattice')

class Lattice:
    @staticmethod
    def fromfile(fname, m):
        kennyloggins.debug("Loading lattice from file " + os.path.abspath(fname))
        f = open(fname)
        #Some lattice files will have the node they ran on printed at the top; bypass 1st line by trying again
        try: n = int(f.readline())
        except ValueError:
            n = int(f.readline())
        h = numpy.zeros((n,n))
        V = None #wait until later to allocate, depending on whether 2ei are local (Hubbard)
        for line in f:
            s = line.split()
            if len(s) == 3:
                i = int(s[0])
                j = int(s[1])
                v = float(s[2])
                h[i,j] = v
#Due to memory constraints, define 2 options, depending on the length of V
            elif s[0] == 'U':
                spline = line.split()
                V = float(spline[1])
            elif len(s) == 5: #Nonlocal
                if V == None: V = numpy.zeros((n,n,n,n))
                i = int(s[0])
                j = int(s[1])
                k = int(s[2])
                l = int(s[3])
                v = float(s[4])
                V[i,j,k,l] = v
            else:
                kennyloggins.warning("Unknown line length in lattice file")
                sys.exit(-1)
        f.close()
        return Lattice(n,h,V,m)

    def __init__(s,n,h,V,m):
        s.m = m
        s.n = n
        s.h = h
        s.V = V
        #I'm eliminating the check for V because sometimes we want it to just be a float
        if not (s.n == len(s.h)):
            kennyloggins.error("Number of sites and integrals misaligned")
            import sys; sys.exit(1)

    def filling(s):
        return float(s.m) / s.n



        
#if __name__ == "__main__":
