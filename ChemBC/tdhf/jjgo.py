#!/usr/bin/python
 
####################################
#
#  CI Singles (CIS) and
#
#  Time-Dependent Hartree Fock (TDHF)
#
####################################
 
from __future__ import division
import math
import numpy as np
 
####################################
#
#   FUNCTIONS
#
####################################
 
# Return compound index given four indices
def eint(a,b,c,d):
  if a > b: ab = a*(a+1)/2 + b
  else: ab = b*(b+1)/2 + a
  if c > d: cd = c*(c+1)/2 + d
  else: cd = d*(d+1)/2 + c
  if ab > cd: abcd = ab*(ab+1)/2 + cd
  else: abcd = cd*(cd+1)/2 + ab
  return abcd
 
# Return Value of spatial MO two electron integral
# Example: (12|34) = tei(1,2,3,4)
def teimo(a,b,c,d):
  return ttmo.get(eint(a,b,c,d),0.0e0)
 
####################################
#
#  INITIALIZE ORBITAL ENERGIES
#  AND TRANSFORMED TWO ELECTRON
#  INTEGRALS
#
####################################
 
Nelec = 2 # we have 2 electrons in HeH+
dim = 2 # we have two spatial basis functions in STO-3G
# Spatial orbital eigenvalues
E = [-1.52378656, -0.26763148]
# Two electron integrals, in a python dictionary
# Using compound indices, see definition of eint()
# and teimo() above.
ttmo = {5.0: 0.94542695583037617, 12.0: 0.17535895381500544, 14.0: 0.12682234020148653, 17.0: 0.59855327701641903, 19.0: -0.056821143621433257, 20.0: 0.74715464784363106}
 
####################################################
#
#  CONVERT SPATIAL TO SPIN ORBITAL MO
#
####################################################
 
# This makes the spin basis double bar integral (physicists' notation)
 
ints=np.zeros((dim*2,dim*2,dim*2,dim*2))
for p in range(1,dim*2+1):
  for q in range(1,dim*2+1):
    for r in range(1,dim*2+1):
      for s in range(1,dim*2+1):
        value1 = teimo((p+1)//2,(r+1)//2,(q+1)//2,(s+1)//2) * (p%2 == r%2) * (q%2 == s%2)
        value2 = teimo((p+1)//2,(s+1)//2,(q+1)//2,(r+1)//2) * (p%2 == s%2) * (q%2 == r%2)
        ints[p-1,q-1,r-1,s-1] = value1 - value2
 
#####################################################
#
#  Spin basis fock matrix eigenvalues
#
#####################################################
 
fs = np.zeros((dim*2))
for i in range(0,dim*2):
    fs[i] = E[i//2]
fs = np.diag(fs)
 
######################################################
#
#   CIS & TDHF CALCULATION
#
######################################################
 
A = np.zeros((2*dim,2*dim))
B = np.zeros((2*dim,2*dim))
I = -1
for i in range(0,Nelec):
  for a in range(Nelec,dim*2):
    I = I + 1
    J = -1
    for j in range(0,Nelec):
      for b in range(Nelec,dim*2):
        J = J+1
        A[I,J] = (fs[a,a] - fs[i,i]) * ( i == j ) * (a == b) + ints[a,j,i,b]
        B[I,J] =  ints[a,b,i,j]
 
# Solve CIS matrix equation
ECIS,CCIS = np.linalg.eig(A)
print "E(CIS) = ", np.amax(ECIS)
# Solve TDHF matrix equation
M = np.bmat([[A,B],[-B,-A]])
ETD,CTD = np.linalg.eig(M)
print "E(TDHF) = ", abs(np.amax(ETD))


