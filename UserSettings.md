## MPSTranslation / MPSDiagonal
#### Must specify
- **BondDim** (`int`): Bond dimension of matrices (constant across the chain)
#### Optional
- **Nspins** (`int`): Number of spins in the chain. Must agree with Hilbert (default from Hilbert)
- **PhysDim** (`int`): Physical dimension of system. Must agree with Hilbert (default from Hilbert)
- **SymmetryPeriod** (`int`): Number of alternating MPS matrices (default is Nspins, resp. no translational symmetry)

## MPSOpen
#### Must specify
- **BondDim** (`int`/`vector<int>`): Bond dimension of matrices.
  - If `int` is given the same bond dimension is used in all matrices.
  - If `vector<int>` is given the bond dimension of each matrix is assigned from this. Length must agree with Nspins.
#### Optional
- **Nspins** (`int`): Number of spins in the chain. Must agree with Hilbert (default from Hilbert)
- **PhysDim** (`int`): Physical dimension of system. Must agree with Hilbert (default from Hilbert)
- (*not working*) **CanonicalForm** (`bool`): Transform to canonical form after every update (default is false)

## SBS
#### Must specify
- One of: 
  - **StringSites** (`vector<vector<int>>`): Distribution of site numbers to strings.
  - **Strings** (`int`): Number of strings.
  If both are given, the sizes must be compatible. If only **Strings** is given all strings cover by default the whole configuration.
- And
  - **BondDim** (`int`/`vector<int>`): Bond dimension of each string (bond dimension is constant across one string like MPSTranslation)
    - If `int` is given the same bond dimension is used in all strings.
    - If `vector<int>` is given the bond dimension of each string is assigned from this. Length must be compatible with number of strings.

#### Optional
- **Nspins** (`int`): Number of spins in the chain. Must agree with Hilbert (default from Hilbert)
- **PhysDim** (`int`): Physical dimension of system. Must agree with Hilbert (default from Hilbert)
- **SymmetryPeriod** (`int`/`vector<int>`): Symmetry period of each string (default is no symmetry in all strings)
  - If `int` is given the symmetry is used in all strings.
  - If `vector<int>` is given the symmetry of each string is assigned from this. Length must agree with number of strings and each symmetry must be smaller or equal to the corresponding string length.
- **Diagonal** (`bool`/`vector<bool>`): Use diagonal MPS (default is false). Vector use is same as above.
