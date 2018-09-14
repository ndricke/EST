import numpy as np
from scipy.special import comb
from itertools import dropwhile
from copy import copy
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

        KcN = [0,0]
        for spin in range(2):
            KcN[spin] = comb(K, N[spin], exact=True) 

        KcN2 = KcN[0] * KcN[1]
        H = np.zeros((KcN2, KcN2))

        #self.back_string = reversed(range(K)) #right now need to keep recreating iterator in exciteString

        string_list = [[],[]] #a list of 2 empty lists
        for spin in range(2):
            next_excitation = np.zeros(K)
            next_excitation[:N[spin]] = 1
            for i in range(KcN[spin]):
                print("iteration: ", i)
                string_list[spin].append(next_excitation) #append before excitation to include ground state
                next_excitation = self.exciteString(next_excitation, self.N[spin])

        for i in range(KcN2):
            for j in range(i+1):
                H[i,j] = compareString(i,j)
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

    def compareString(self, string1_index, string2_index):
        truth_arr = np.equal(string1, string2)
        ##Compare to see how many excitations are conserved
        excitation_index = [i for i, b in enumerate(truth_arr) if b == False]

        if len(excitation_index) > 2:
            return 0
        elif len(excitation_index) == 2:
            return doubleExcitation(excitation_index)
        elif len(excitation_index) == 1:
            return singleExcitation(excitation_index)
        elif len(excitation_index) == 0:
            return evalEnergy(string1)

    def doubleExcitation(self):
        pass




        #return next_string


if __name__ == "__main__":
    bob = FCI(6, [3,3], h=1, V=1)


















