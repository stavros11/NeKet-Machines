# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 21:58:33 2018

@author: Stavros
"""

import numpy as np

def calculate_psi(mps, v):
    ## Returns <v|MPS> for a given |v> in computational basis
    prod = np.copy(mps[0][v[0]])
    for i in range(1, len(v)):
        prod = prod.dot(mps[i][v[i]])
    return prod / np.sqrt(calculate_inner_product(mps).real)

def calculate_inner_product(mps):
    ## Returns <PSI|PSI> where |PSI> is in MPS form
    totem = np.einsum('iab,icd->abcd', np.conj(mps[0]), mps[0])
    extended_totem = np.einsum('abcd,ibe->icade', totem, mps[1])
    for i in range(1, len(mps)-1):
        totem = np.einsum('iabcd, ice->abed', extended_totem, np.conj(mps[i]))
        extended_totem = np.einsum('abcd, ide->iabce', totem, mps[i+1])
    totem = np.einsum('iabcd, ice->abed', extended_totem, np.conj(mps[-1]))
    return np.einsum('abab', totem)

def str2complex(x):
    y = x[1:-1].split(',')
    return float(y[0]) + 1j * float(y[1])

def mps_from_file(filename):
    file_str = open(filename).read().split(' \n\n')[:-1]
    
    k = 0
    mps = []
    for i in range(N):
        pars = D[i] * D[i+1] * d
        comp_array = [str2complex(x) for x in file_str[k : k + pars]]
        mps.append(np.array(comp_array).reshape(d, D[i], D[i+1]))
        k += pars
        
    return mps


N, Duser, d = 6, 4, 2
D = [1, d] + [Duser for i in range(N-3)] + [d, 1]

ver = 3
mps_init = mps_from_file('cpp/tests/Winit%d.txt'%ver)
mps_norm = mps_from_file('cpp/tests/Wsvd%d.txt'%ver)
print("Before SVD: " + str(calculate_inner_product(mps_init)))
print("After SVD:" + str(calculate_inner_product(mps_norm)))