# bdistrib

Parameter: nu_plasma
Type: int
Default: 16
When it matters: always
Description: Number of grid points in the poloidal direction on the plasma surface

Parameter: nv_plasma
Type: int
Default: 32
When it matters: always
Description: Number of grid points in the toroidal direction (per period) on the plasma surface

Parameter: surface_option_current
Type: int
Default: 0
When it matters: always
Description:
 0 - The current surface will be a torus. The major radius will be Rmajor_p from the wout file.
     The minor radius will be a_current.
 1 - Same as option 0, except the major radius will be R_current.
 2 - The current surface will computing by expanding the plasma LCFS uniformly by a distance separation_current.


