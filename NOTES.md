The uploaded files are the following:
- `mps_periodic.hpp`: Implements a periodic MPS machine where all matrices are square with the same (bond) dimension given by user. It is possible to set translational symmetry (same set of matrices repeated across the chain).
- `mps_diagonal.hpp`: Implements a machine identical to `mps_periodic.hpp` where all matrices are diagonal.
- `abstract_mps.hpp`: Contains virtual methods so that the above classes can be used for SBS calculations.
- `sbs.hpp`: Implements a SBS (String Bond State = Product of MPS strings) machine. All SBS strings are periodic.

**WARNING: Unit test failures**
- The real test of `DerLog` for all machines. The imaginary part has no problems and the machines seem to converge to the correct energies for small (`N<10`) 1D Ising and Heisenberg Hamiltonians. I am not sure where the unit test error comes from - in my computer the Jastrow machines fail the same test too.
- I use exceptions in the SBS machine which I think are not allowed in the Google style. This is only in the `from_json` function and does not affect calculations, so it will be replaced soon (sorry for that).

Some additional notes:
1. The look up tables store all left and right contractions. This is efficient in the one-flip case but might not be the optimal when more than one spins are flipped. Storing middle contractions could assist in computations for large systems but would add time in the look up update. It is not straightforward (to me) what is the best to do. Optimizing this is among my priorities.
2. All contractions are calculated from scratch every time `DerLog` is called. [Using look up tables in `DerLog`](https://github.com/netket/netket/issues/73) would help this and other machines, but this requires small changes in the VMC part and possibly other parts of the code (I haven't forgot about that). 
3. User can explicitly assign which sites go to each string in SBS through the JSON. This allows the greatest flexibility but might be a bit hard/non-automatic from the user's perspective and also allows non-physical cases (eg. strings that jump sites). An easy way to bypass this is by writing a simple [Python wrapper](https://github.com/netket/netket/pull/65) that creates the "sites to strings" lists for the JSON. A more formal way could be to use the graph object directly in C++. I am not yet sure what's the best way to simplify things while still keeping the settings quite general.
4. The "sites to strings" option in SBS can be used to generalize MPS too (1-string SBS = MPS). Hence the SBS machine can be used to get different MPS shapes (eg. snake MPS), which is not supported in my pure MPS implementations.
5. Weight saving/loading is not tested.
6. My machines are possibly slower than all current machines and the number of parameters is typically larger (eg. ~N^3 for MPS vs ~N^2 for RBM). In any case, I am far from C++ expert, so please comment if you find any "newbie" errors in my code. It seems that the diagonal version is faster even for small systems, which might be just because the number of parameters is smaller or because my contractions are inefficient.
7. I have an implementation for open boundary MPS which allows different matrix dimensions across the chain, but it might be better to focus on the submitted ones for now. If there is interest I add that one too, probably after I manage to add a feature to transform it to its canonical form.
