### Currently working
#### Final
- `mps_open.hpp`: Open MPS with vector D.
- `mps_periodic.hpp`: Translation invariant periodic MPS with arbitrary symmetry period.
- `mps_diagonal.hpp`: Diagonal MPS (Same as `MPSTranslation` but matrices are diagonal and we use `cwiseProduct` instead of matrix product - possibly a bit faster).
- `sbs.hpp`: SBS.

All tested in 1D Ising and Heisenberg. Not sure about BoseHubbard1D as RBM does not converge either.

#### Dependencies
- `MPSTranslation` and `MPSDiagonal` are derived from `AbstractMPS`.
- SBS uses `MPSTranslation` or `MPSDiagonal` (user selects) for MPS calculations.

#### Deprecated
- Open MPS with constant D. Keep the vector D version.
- Periodic MPS. Keep the translational symmetric version and set `symmperiod_=N_`.

### Additions needed

#### Physics related
- [ ] Canonical form in Open MPS (doesn't work in the Python code either).
- [ ] Possible improvement of the lookups in all cases for the calculation of middle contractions. Currently we save only left and right contractions. Not sure if saving middle contractions will improve efficiency.
- [x] Setting for ordering of spins and how to put the strings into graph. Possible connection to graph object.
  - Currently implemented with the user **explicitly** giving the site distribution in strings in the json.
  - Easy improvement using a Python wrapper.
  - Would it help to connect to the `graph` object within C++?

#### Code related
- [ ] Test weight loading and enable weight saving (incompatibilities with json format)
- [ ] Currently SBS uses MPS classes for calculations, but each MPS function is defined twice in the class, once for pure MPS and once for SBS . It might be possible to combine some of these functions.
  - Incompatibility with `confindex_` calculation. It might be less efficient to combine, but possibly negligible difference.
- [ ] Use of lookups in the derivative (is this possible?). Currently we do the contractions from scratch in the `DerLog` functions.
  - Easy fix requires editing the VMC part of the code.
- [x] Setting for different string lengths and bond dimensions in SBS.
