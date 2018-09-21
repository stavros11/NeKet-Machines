#### Functions that change from Periodic to Open MPS
- Init
- GetParameters
- SetParameters
- SetParametersIdentity
- InitLookup
- LogVal (because the product changes type from (1, D) to (1, 1))
- LogValDiff full
- DerLog

#### New functions for Open MPS
- InitLookup(Left/Right/Boundary)_check instead of InitLookup_check
- mps_contraction(Left/Right)
- tensor product for derivative (cannot do this in eigen because (D, 1) is stored as (D, D))
