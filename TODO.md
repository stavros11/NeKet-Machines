### Currently working
- Open MPS, versions with constant and vector D.
- Periodic MPS.
- Translation invariant periodic MPS with arbitrary symmetry period.
- General SBS without lookups.

All tested only in 1D Ising and Heisenberg. Not sure about BoseHubbard1D as RBM does not converge either.

### Implemented but not working
- Lookups in SBS: Works only for the special case M=1 (which is MPS). Possibly small bug.
- Canonical form in Open MPS (doesn't work in Python either).

### Additions needed
- Setting to save weights for all cases (it is currently commented).
- Possible use of the `MPSCalculator` class in Periodic and Translation MPS for better looking code. Currently it is used only in MPS.
- Possible improvement of the lookups in all cases for the calculation of middle contractions. Currently we save only left and right contractions.
- Use of lookups in the derivative (is this possible?). Currently we do the contractions from scratch in the `DerLog` functions.
- Special case of diagonal SBS.
- Setting for ordering of spins and how to put the strings into graph. Possible connection to graph object. Currently SBS cover the whole graph by default.
