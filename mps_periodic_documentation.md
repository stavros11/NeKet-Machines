---
title: MPS
permalink: /docs/mps/
---

<h2 class="bg-primary">Periodic MPS</h2>
The periodic MPS machine is defined by:

$$
\Psi(s_1, \dots , s_N) = \mathrm{Tr}\left [ A^{(s_1)}_1A^{(s_2)}_2\dots A^{(s_N)}_N \right ]
$$

for arbitrary local quantum numbers $$ s_i $$. $$ A^{(s_i)} $$ are complex $$ D\times D $$ matrices 

The total number of variational parameters of this long-range Jastrow operator is $$ \mathcal{O} (N (N-1)/2) $$, and all parameters are taken complex-valued. This wavefunction needs less input parameters compared to the RBM wavefunction, since the hidden layer is not present. In the language of neural networks a Jastrow machine can be seen as a fully visible RBM, with all the intra-layer connections active.



|---
| Parameter | Possible values | Description | Default value |
|-|-|-|-
| `InitFile` | String |  If specified, matrix parameters are loaded from the given file | None |
| `InitRandom` | Boolean |  Whether to initialize the parameters with random gaussian-distributed values | True |
| `SigmaRand` | Float |  If InitRandom is chosen, this is the standard deviation of the gaussian  | 0.1 |
|===

### Example
```python
pars['Machine']={
    'Name'           : 'Jastrow',
}
```

<h2 class="bg-primary">JastrowSymm</h2>
This type of Jastrow wavefunction is constructed to be invariant under specific permutations of the
graph indices. In particular, $$ \Psi(s_1,\dots s_N) = \Psi(s_{P_k(1)} \dots s_{P_k(N)}) $$, where $$ P_k(i) $$ is the $$ k $$ -th permutation of the index $$ i $$,
and assuming a total of $$ N_p $$ distinct permutations.  

Let us consider the case of a N=20 spin chain, under periodic boundary conditions. The total number of parameters is $$N (N-1)1/2 =190$$, corresponding to the full upper triangle (excluding the diagonal, always set to $$0$$) to the $$W$$ matrix. However it is true that, if permutation symmetry holds, the parameter $$W[1,2]$$  controlling the interaction between quantum operator $$ s_1 $$ and $$ s_2 $$ must be equal to $$W[2,3]$$, and so on.
Therefore, in this case, the effective number of independent parameter is no longer $$190$$ but $$10$$.
In the general case, this number is automatically computed inside the class by checking the symmetry table associated with the graph.




|---
| Parameter | Possible values | Description | Default value |
|-|-|-|-
| `InitFile` | String |  If specified, matrix parameters are loaded from the given file | None |
| `InitRandom` | Boolean |  Whether to initialize the parameters with random gaussian-distributed values | True |
| `SigmaRand` | Float |  If InitRandom is chosen, this is the standard deviation of the gaussian  | 0.1 |
|===

### Example
```python
pars['Machine']={
    'Name'           : 'JastrowSymm',
}
```




## References
---------------
1. [W. L. McMillan, Phys. Rev. 138, A442 (1965), Ground State of Liquid
$$He_4$$](https://journals.aps.org/pr/abstract/10.1103/PhysRev.138.A442)
