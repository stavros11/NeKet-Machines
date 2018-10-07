**WARNING: The machines fail one unit test**, the real test for `DerLog`. The imaginary part seems to have no problems and they seem to converge to the correct energies for small (`N<10`) 1D Ising and Heisenberg Hamiltonians. I am not sure where the unit test comes from - in my computer also the Jastrow machines fail the same test.

The uploaded files are the following:
- `mps_periodic.hpp`: Implements a periodic MPS machine where all matrices are square with the same (bond) dimension given by user. It is possible to set translational symmetry (same set of matrices repeated across the chain).
- `mps_diagonal.hpp`: Implements a machine identical to `mps_periodic.hpp` where all matrices are diagonal.
- `abstract_mps.hpp`: Contains virtual methods so that the above classes can be used for SBS calculations.
- `sbs.hpp`: Implements a SBS (String Bond State = Product of MPS strings) machine. All SBS strings are periodic.
- `UserSettings.md`: Explanations of the different JSON settings currently implemented.

Some notes:
1. Currently contractions are calculated from scratch every time `DerLog` is called. As mentioned in ..., using look up tables in `DerLog` could speed up some machines, but this would require small changes in the VMC part, which I plan to try. 
2. The look up tables currently used involve storing all left and right contractions. This is efficient in the one-flip case but could be inefficient when more than one spins are flipped. It is not straightforward (to me) whether storing middle contractions can assist or no in the latter case, so that's another possible improvement I will try to look at.
3. Currently user can explicitly assign which sites go to each string in SBS through the JSON. This allows the greatest flexibility but might be a bit hard/non-automatic from the user's perspective and also allows non-physical cases (eg. strings that jump sites). An easy way to bypass this is through a simple Python wrapper that creates the "sites to strings" lists for the JSON. A more formal way could be to use the graph object directly in C++. I am not yet sure what's the best way to simplify things while still keeping it quite general.
4. Following 2, the "sites to strings" option in SBS can be used to generalize MPS too (1-string SBS = MPS). Hence the SBS machine can be used to get different MPS shapes (eg. snake MPS), which is not supported in my MPS implementations.
5. The machines are not very efficient, possibly slower than all current machines, but they also provide generalizations (eg. diagonal SBS = RBM). Under basic assumptions the number of parameters is typically larger: ~N^3 for MPS vs ~N^2 for RBM.
6. I have an implementation for open boundary MPS which allows different matrix dimensions across the chain, but it might be better to focus on the current ones for now. If there is interest I can try adding that one too, after I manage to implement its canonical form.
