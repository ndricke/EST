import numpy as np
from scipy.special import comb
from itertools import dropwhile
from copy import copy
import sys

from embed import c1c2
"""
1. need h_a, h_b, V_aa, V_ab, V_bb for systems with different numbers of spin
   for the moment, I will assume that N[0] == N[1]
"""
class FCI:

    def __init__(self, K, N, h, V):
        self.N = N
        self.K = K
        self.h = h
        self.V = V

        self.KcN = [0,0]
        for spin in range(2):
            self.KcN[spin] = comb(K, N[spin], exact=True)

        KcN2 = self.KcN[0] * self.KcN[1]
        H = np.zeros((KcN2, KcN2))

        #self.back_string = reversed(range(K)) #right now need to keep recreating iterator in exciteString

        string_list = [[],[]] #a list of 2 empty lists
        for spin in range(2):
            next_excitation = np.zeros(K, dtype=int)
            next_excitation[:N[spin]] = 1
            for i in range(self.KcN[spin]):
                print("iteration: ", i)
                string_list[spin].append(next_excitation) #append before excitation to include ground state
                next_excitation = self.exciteString(next_excitation, self.N[spin])

        num_str_a = len(string_list[0]) #number of alpha strings
        num_str_b = len(string_list[1]) #number of beta strings
        self.string_ab = np.zeros((num_str_a*num_str_b, self.K*2)) #array for combined alpha and beta strings
        for i in range(num_str_a):
            for j in range(num_str_b):
                row = j + i*num_str_b
                self.string_ab[row,:self.K] = string_list[0][i]
                self.string_ab[row,self.K:] = string_list[1][j]

        for i in range(KcN2):
            for j in range(i+1):
                H[i,j] = self.compareString(i,j)
                H[j,i] = H[i,j]

        print(H)

    def exciteString(self, string, N):
        next_string = copy(string)
        rightmost_1 = next(dropwhile(lambda x: string[x] != 1, reversed(range(self.K))))

        #if the highest 1 is not at the top, move it up 1
        if rightmost_1 != self.K-1:
            next_string[rightmost_1] = 0
            next_string[rightmost_1+1] = 1
            return next_string

        else:
            #if the highest 1 is at the top, look for the highest 0
            rightmost_0 = next(dropwhile(lambda x: string[x] != 0, reversed(range(self.K))))
            #if all the ones are at the top, return error (cannot excite further)
            if rightmost_0 == self.K - N - 1: #string is already fully excited
                return ValueError("String already fully excited")

            #if some (but not all) are the top, move the highest non-top 1 up by 1, and put the former top 1s onto it
            next_1 = next(dropwhile(lambda x: string[x] != 1, reversed(range(rightmost_0))))
            print("Next to excite: ", next_1)
            print("Rightmost 0: ", rightmost_0)
            next_string[next_1] = 0
            next_string[rightmost_0:] = 0
            next_string[next_1+1: next_1+1+self.K-rightmost_0] = 1
            return next_string

    ## This was from when I had the alpha and beta strings as separate lists. I may go back to it
    #def calcStringdex(self, i):
    #    return i//self.KcN[spin], i % self.KcN[spin]

    def excitOperators(self, str_1, str_2):
        """
        params:
            str_1, str_2 : (np.array) occupations for excitation strings
        ----
        returns:
             (list of tuples) tuples containing indices for (creation, destruction) operators
        """
        str_dif = str_1 - str_2
        lowering_dex = np.where(str_dif == 1)[0]
        raising_dex = np.where(str_dif == -1)[0]
        return zip(raising_dex, lowering_dex)

    def compareString(self, stringdex_1, stringdex_2):
        """
        params:
            stringdex_1, stringdex_2 : (int) index for excitation strings in string_list
        returns:
            (int) FCI Hamiltonian coefficient
        """
        credes_op = self.excitOperators(self.string_ab[stringdex_1], self.string_ab[stringdex_2])

        excit_count = len(credes_op) #number of excitations between alpha and beta strings
        if excit_count > 2:
            return 0.0
        elif excit_count == 2:
            return self.doubleExcitation(credes_op)
        elif excit_count == 1:
            return self.singleExcitation(stringdex_1, stringdex_2, credes_op)
        elif excit_count == 0:
            return self.evalEnergy(stringdex_1)
        else:
            raise ValueError("Excitation count is not a whole number")

    def evalEnergy(self, stringdex): #if there are no excitations, only 1 string need be passed
        """
        params:
            stringdex : (list) excitation string index
        returns:
            (double) : E, the energy of this string
        """
        ex_string = self.string_ab[stringdex,:]
        E = 0.0
        occ_ind = np.where(ex_string == 1)[0]
        for a in occ_ind:
            E += self.h[a,a]
        for a,b in itertools.product(occ_ind,occ_ind):
            E += 0.5*(self.V[a,a,b,b]-self.V[a,b,b,a]) ## check: is this correct for the unrestricted basis?
        return E

    def singleExcitation(self, stringdex_1, stringdex_2, credes_op):
        """
        params:
            str_1 : (list) excitation string
        returns:
            (double) : E, the energy of this string
        """
        occ_ind_1 = np.where(self.string_ab[stringdex_1] == 1)[0]
        occ_ind_2 = np.where(self.string_ab[stringdex_2] == 1)[0]

        E = 0.0
        E += self.h[credes_op[0][0], credes_op[0][1]] #only one 1-particle term

        m = credes_op[0]  #creation operator
        p = credes_op[1]  #destruction operator
        for n in occ_ind_1: ## check: is this correct given each string has different occupations?
            E += self.V[m,p,n,n] - self.V[m,n,n,p]
        return E


    def doubleExcitation(self, credes_op):
        m = credes_op[0][0]
        n = credes_op[0][1]
        p = credes_op[1][0]
        q = credes_op[1][1]
        return V[m,p,n,q] - V[m,q,n,p]





if __name__ == "__main__":
    n = 6
    U = 2.0
    na, nb = 3, 3

    h, V = c1c2.hubb(n, U, pbc=True)
    fci_6 = FCI(n, [na,nb], h=h, V=V)
