# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 14:28:27 2018

@author: Admin
"""

import numpy as np
from numpy.linalg import svd

N, d, D = 6, 2, 4
sigma = 0.1
PBC = False

def normal_complex(sigma, shape):
    return np.random.normal(0, sigma, size=shape) + 1j * np.random.normal(0, sigma, size=shape)

def calculate_psi(mps, v):
    ## Returns <v|MPS> for a given |v> in computational basis
    prod = np.copy(mps[0][v[0]])
    for i in range(1, len(v)):
        prod = prod.dot(mps[i][v[i]])
    return prod / np.sqrt(calculate_inner_product(mps).real)
        

def calculate_inner_product(mps):
    ## Returns <MPS|MPS> where |PSI> is in MPS form
    totem = np.einsum('iab,icd->abcd', np.conj(mps[0]), mps[0])
    extended_totem = np.einsum('abcd,ibe->icade', totem, mps[1])
    for i in range(1, len(mps)-1):
        totem = np.einsum('iabcd, ice->abed', extended_totem, np.conj(mps[i]))
        extended_totem = np.einsum('abcd, ide->iabce', totem, mps[i+1])
    totem = np.einsum('iabcd, ice->abed', extended_totem, np.conj(mps[-1]))
    return np.einsum('abab', totem)

def PBC_2canonical(mps):
    ## Puts states into canonical MPS form
    states = np.copy(mps)
    for i in range(N - 1):
        U, S, V = svd(states[i].reshape([D * d, D]), full_matrices=False)
        states[i] = U.reshape([d, D, D])
        states[i+1] = np.einsum('ab, ibc->iac', np.diag(S).dot(V), states[i+1])
        
    U, S, V = svd(states[-1].reshape([D * d, D]), full_matrices=False)
    states[-1] = U.reshape([d, D, D])
    states[0] = np.einsum('ab, ibc->iac', np.diag(S).dot(V), states[0])
    
    return states
    
def OBC_2canonical(mps):
    ## Puts states into canonical MPS form
    ## Warning: now mps is a list, not an numpy array
    U, S, V = svd(mps[0].reshape([d, d]), full_matrices=False)
    ## Warning states[0] changes shape from (d, 1, D) to (d, 1, d)
    states = [U.reshape([d, 1, d])]

    states.append(np.einsum('ab, ibc->iac', np.diag(S).dot(V), mps[1]))
    U, S, V = svd(states[1].reshape([d * d, D]), full_matrices=False)
    states[1] = U.reshape([d, d, D])
    for i in range(2, N - 2):
        states.append(np.einsum('ab, ibc->iac', np.diag(S).dot(V), mps[i]))
        U, S, V = svd(states[i].reshape([d * D, D]), full_matrices=False)
        states[i] = U.reshape([d, D, D])
    states.append(np.einsum('ab, ibc->iac', np.diag(S).dot(V), mps[N-2]))
    U, S, V = svd(states[N-2].reshape([d * D, d]), full_matrices=False)
    states[N-2] = U.reshape([d, D, d])
    
    last = np.einsum('ab, ibc->iac', np.diag(S).dot(V), mps[N-1])
    #states.append(last / np.sqrt(np.sum(np.conj(last) * last).real))
    states.append(last)
    
    #U, S, V = svd(states[-1][:, :, 0], full_matrices=False)
    #states[-1] = U[:, :, np.newaxis] / np.sqrt(d)
    return states

if PBC:
    states = normal_complex(sigma, shape=(N, d, D, D))
    to_canonical = PBC_2canonical
else:
    states = [normal_complex(sigma, shape=(d, 1, d)), normal_complex(sigma, shape=(d, d, D))]
    states += [normal_complex(sigma, shape=(d, D, D)) for i in range(N-4)]
    states += [normal_complex(sigma, shape=(d, D, d)), normal_complex(sigma, shape=(d, d, 1))]
    to_canonical = OBC_2canonical

#print('Before SVD: ' + str(calculate_inner_product(states)))
#states_can = to_canonical(states)
#print('After SVD: ' + str(calculate_inner_product(states_can)))

## Check whether SVD changes the state
states_up = [x for x in states]

#U, S, V = svd(states_up[0][:, 0, :], full_matrices=False)
#states_up[0] = U[:, np.newaxis, :]
#states_up[1] = np.einsum('ab,ibc->iac', np.diag(S).dot(V), states_up[1])

U, S, V = svd(states_up[1].reshape([d * d, D]), full_matrices=False)
states_up[1] = U.reshape([d, d, D])
states_up[2] = np.einsum('ab,ibc->iac', np.diag(S).dot(V), states_up[2])

#U, S, V = svd(states_up[2].reshape([d * D, D]), full_matrices=False)
#states_up[2] = U.reshape([d, D, D])
#states_up[3] = np.einsum('ab,ibc->iac', np.diag(S).dot(V), states_up[3])
#
#U, S, V = svd(states_up[3].reshape([d * D, D]), full_matrices=False)
#states_up[3] = U.reshape([d, D, D])
#states_up[4] = np.einsum('ab,ibc->iac', np.diag(S).dot(V), states_up[4])
#
#U, S, V = svd(states_up[4].reshape([d * D, d]), full_matrices=False)
#states_up[4] = U.reshape([d, D, d])
#states_up[5] = np.einsum('ab,ibc->iac', np.diag(S).dot(V), states_up[5])

#states_up[2] = states_up[2] / np.sqrt(calculate_inner_product(states_up).real)

for i in range(N):
    print((states[i] - states_up[i]).mean())
    
print('\nNorms:')
print(calculate_inner_product(states))
print(calculate_inner_product(states_up))

print('\nPsi:')
vis = np.random.randint(0, d, N)
print(vis)
print(calculate_psi(states, vis))
print(calculate_psi(states_up, vis))


