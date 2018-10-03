## MPSTranslation / MPSDiagonal
#### Must specify
- BondDim (`int`): Bond dimension of matrices (constant across the chain)
#### Optional
- Nspins (`int`): Number of spins in the chain. Must agree with Hilbert (default from Hilbert)
- PhysDim (`int`): Physical dimension of system. Must agree with Hilbert (default from Hilbert)
- SymmetryPeriod (`int`): Number of alternating MPS matrices (default is Nspins, resp. no translational symmetry)

## MPSOpen
#### Must specify
- BondDim (`int`/`vector<int>`): Bond dimension of matrices.
  - If `int` is given the same bond dimension is used in all matrices.
  - If `vector<int>` is given the bond dimension of each matrix is assigned from this. Length must agree with Nspins.
#### Optional
- Nspins (`int`): Number of spins in the chain. Must agree with Hilbert (default from Hilbert)
- PhysDim (`int`): Physical dimension of system. Must agree with Hilbert (default from Hilbert)
- (*not working*) CanonicalForm (`bool`): Transform to canonical form after every update (default is false)

## SBS
#### Must specify
- Nstrings (`int`): Number of strings.
- BondDim (`int`/`vector<int>`): Bond dimension of each string (bond dimension is constant across one string like MPSTranslation)
  - If `int` is given the same bond dimension is used in all strings.
  - If `vector<int>` is given the bond dimension of each string is assigned from this. Length must agree with Nstrings.
- (**TODO**) `sites2strings` and `strings2sites`. This will define `Lstr_` too.
#### Optional
- Nspins (`int`): Number of spins in the chain. Must agree with Hilbert (default from Hilbert)
- PhysDim (`int`): Physical dimension of system. Must agree with Hilbert (default from Hilbert)
- Diagonal (`bool`/`vector<bool>`): Use diagonal MPS (default is false)
- SymmetryPeriod (`int`/`vector<int>`): Symmetry period of each string (default is no symmetry in all strings)
  - If `int` is given the symmetry is used in all strings.
  - If `vector<int>` is given the symmetry of each string is assigned from this. Length must agree with Nstrings and each symmetry must be smaller or equal to the corresponding string length.
