## MPSTranslation / MPSDiagonal
#### Must specify
- BondDim (`int`): Bond dimension of matrices (constant across the chain)
#### Optional
- Nspins (`int`): Number of spins in the chain. Must agree with Hilbert (default from Hilbert)
- PhysDim (`int`): Physical dimension of system. Must agree with Hilbert (default from Hilbert)
- SymmetryPeriod (`int`): Number of alternating MPS matrices (default is Nspins, resp. no translational symmetry)

## MPSOpen
#### Must specify
- BondDim (`int`/`list`): Bond dimension of matrices.
  - If `int` is given the same bond dimension is used in all matrices.
  - If `list` is given the bond dimension of each matrix is assigned from this. Length must agree with Nspins.
#### Optional
- Nspins (`int`): Number of spins in the chain. Must agree with Hilbert (default from Hilbert)
- PhysDim (`int`): Physical dimension of system. Must agree with Hilbert (default from Hilbert)
- CanonicalForm (`bool`): Transform to canonical form after every update (default is false) (*not working*)

