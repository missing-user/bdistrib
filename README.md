# bdistrib

There are 3 surfaces, named 'plasma', 'middle', and 'current',
The 'plasma' surface corresponds to the outermost surface of the VMEC equilibrium.
The name of the relevant surface appears as a suffix on most variables.

On each surface we use a poloidal angle u and a toroidal angle v, defined as in NESCOIL.
The coordinate u lies in the range [0,1]. The stellarator has nfp identical toroidal periods,
and an increase in v by 1 corresponds to 1 of these periods. Thus, v increases by nfp
in a complete toroidal revolution. In the output file, there is an array v corresponding
to one toroidal period, as well as an array vl corresponding to all nfp toroidal periods.
Note that v is proportional to the standard cylindrical angle phi: phi = 2*pi*v/nfp.
Generally, field lines are not straight in the (u,v) coordinates.  On the plasma surface,
u is identical to the VMEC poloidal angle divided by 2*pi, while the VMEC toroidal angle
differs by a sign relative to v.

In the various variable names in the code and output file, 'r' refers to the position
vector, not to a radius.  In various arrays with a dimension of length 3, this dimension
always corresponds to Cartesian coordinates (x,y,z).

The 'normal' quantities in the code and output file refer to the surface normal vector
N = (dr/dv) cross (dr/du) as in the NESCOIL paper. Note that this vector does not have unit magnitude.

The current potential on the outermost surface is described by a Fourier expansion with
variable names analogous to VMEC but beginning 'currentPotential_'. Hence, the maximum
poloidal and toroidal mode numbers are currentPotential_mpol and currentPotential_ntor respectively,
the total number of modes is currentPotential_mnmax, and the poloidal and toroidal mode
numbers are stored in arrays currrentPotential_xm and currentPotential_xn respectively.
As in VMEC, the toroidal mode numbers can be positive or negative or zero, the poloidal mode
numbers are non-negative.  For poloidal mode numbers of zero, the toroidal mode number
must be positive, rather than non-negative as in VMEC. This difference arises because
the current potential is only defined up to a constant, so we force the constant mode to be zero.


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


