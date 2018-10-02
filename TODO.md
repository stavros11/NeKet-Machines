### Currently working
- Open MPS, versions with constant (*deprecated*) and vector D (`mps_open.hpp`).
- Periodic MPS (`mps_periodic.hpp`) - This is special case of `MPSTranslation` with `symperiod_ = N`.
- Translation invariant periodic MPS with arbitrary symmetry period (`mps_translation.hpp`).
- Diagonal MPS (`mps_diagonal.hpp`) - Same as `MPSTranslation` but matrices are diagonal and we use `cwiseProduct` instead of matrix product (a bit faster).
- SBS (`sbs.hpp`). Uses `MPSTranslation` or `MPSDiagonal` classes for the MPS calculations. These are defined by the `AbstractMPS` class.

All tested in 1D Ising and Heisenberg. Not sure about BoseHubbard1D as RBM does not converge either.

### Implemented but not working
- Canonical form in Open MPS (doesn't work in Python either).

### Additions needed

#### Physics related
- Possible improvement of the lookups in all cases for the calculation of middle contractions. Currently we save only left and right contractions.
- Setting for ordering of spins and how to put the strings into graph. Possible connection to graph object. Currently SBS cover the whole graph by default.

#### Code related
- Currently SBS uses MPS classes for calculations, but each MPS function is defined twice in the class, once for pure MPS and once for SBS (incompatibility with `confindex_` calculation). It might be possible to combine some of these functions.
- Setting to save weights for all cases (it is currently commented).
- Use of lookups in the derivative (is this possible?). Currently we do the contractions from scratch in the `DerLog` functions.
- Special case of diagonal SBS.
- Setting for different dimensions in open MPS and different string lengths and bond dimensions in SBS. Code supports that but currently it is not controlled by user.

