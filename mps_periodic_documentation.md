---
title: MPS
permalink: /docs/mps_periodic/
---

<h2 class="bg-primary">Periodic MPS</h2>
The periodic MPS machine is defined by:

$$
\Psi(s_1, \dots , s_N) = \mathrm{Tr}\left [ A^{(s_1)}_1A^{(s_2)}_2\dots A^{(s_N)}_N \right ]
$$

where $$ A^{(s_i)}_{i} $$ are complex $$ D\times D $$ matrices and $$ s_i $$ are arbitrary local quantum numbers. $$ D $$ is known as the bond dimension and is the same for all the matrices. The total number of parameters is $$ N D^2 $$. It is possible to reduce the number of parameters to $$ N D $$ using diagonal matrices, or to enforce translational symmetry by using the same matrix in all sites. Generally it is possible to repeat $$ n < N $$ distinct matrices to cover all the sites, where $$ n $$ is a divisor of $$ N $$.

|---
| Parameter | Possible values | Description | Default value |
|-|-|-|-
| `InitFile` | String |  If specified, matrix parameters are loaded from the given file | None |
| `InitRandom` | Boolean |  Whether to initialize the parameters with random gaussian-distributed values | True |
| `SigmaRand` | Float |  If InitRandom is chosen, this is the standard deviation of the gaussian  | 0.1 |
| `BondDim` | Integer |  Bond dimension | None |
| `Diagonal` | Boolean |  Use diagonal matrices | False |
| `SymmetryPeriod` | Integer |  The number $$ n $$ of distinct matrices to repeat (for translational symmetry) | $$ N $$ |
|===

### Example
```python
pars['Machine']={
    'Name'           : 'MPSperiodic',
    'BondDim'        : 10,
}
```

## References
---------------
1. [Cirac J. I. & Verstraete F., Journal of Physics A: Mathematical and Theoretical (2009), Renormalization and tensor product states in spin chains and lattices](http://iopscience.iop.org/article/10.1088/1751-8113/42/50/504004)
2. [Sandvik A. W. & Vidal G., Phys. Rev. Lett. 99, 220602 (2007), Variational Quantum Monte Carlo Simulations with Tensor-Network States](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.220602)
