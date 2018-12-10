import numpy as np
from scipy.special import comb
from itertools import dropwhile
from copy import copy
import sys
"""
1. Need to implement this to work properly over both spins
2. Create direct product of alpha and beta spin determinants to make full determinant space
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

        print(string_list)
        print(len(string_list))
        print(len(string_list[0]))
        sys.exit(-1)

        for i in range(KcN2):
            for j in range(i+1):
                H[i,j] = self.compareString(i,j)
                H[j,i] = H[i,j]

    def exciteString(self, string, N):
        next_string = copy(string)
        print(string)
        rightmost_1 = next(dropwhile(lambda x: string[x] != 1, reversed(range(self.K))))
        print(rightmost_1)

        #if the highest 1 is not at the top, move it up 1
        if rightmost_1 != self.K-1:
            print("not at top")
            next_string[rightmost_1] = 0
            next_string[rightmost_1+1] = 1
            return next_string

        else:
            print("Hit the top!")
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

    def calcStringdex(self, i):
        return i//self.KcN[spin], i % self.KcN[spin]

    def excitOperators(self, str_1, str_2):
        """

        ----
        returns:
            lowering_dex : (int) index for lowering operator
            raising_dex : (int) index for raising operator
        """
        str_dif_a = str_1a - str_2a
        lowering_dex = np.where(str_dif_a == 1)[0]
        raising_dex = np.where(str_dif_a == -1)[0]
        return lowering_dex, raising_dex

    def compareString(self, stringdex_1, stringdex_2):

        str_ind_1a, str_ind_1b = self.calcStringdex(stringdex_1)
        str_ind_2a, str_ind_2b = self.calcStringdex(stringdex_2)
        str_1a, str_1b = string_list[str_ind_1a], string_list[str_ind_1b]
        str_2a, str_2b = string_list[str_ind_2a], string_list[str_ind_2b]

        al_a, ar_a = self.excitOperators(str_1a, str_2a)
        al_b, ar_b = self.excitOperators(str_1b, str_2b)

        excit_count = len(al_a) + len(al_b) #number of excitations between alpha and beta strings
        if excit_count > 2:
            return 0.0
        elif excit_count == 2:
            return self.doubleExcitation(al_a, al_b, ar_a, ar_b, str_1a, str_1b, str_2a, str_2b)
        elif excit_count == 1:
            return self.singleExcitation(al_a, al_b, ar_a, ar_b, str_1a, str_1b, str_2a, str_2b)
        elif excit_count == 0:
            return self.evalEnergy(str_1a, str_1b)
        else:
            raise ValueError("Excitation count is not a whole number")

    def evalEnergy(self, string_a, string_b):


    def singleExcitation(self, str_1a, str_1b, str_2a, str_2b):

    def doubleExcitation(self):
        pass




        #return next_string


if __name__ == "__main__":
    bob = FCI(6, [3,3], h=1, V=1)
