# bdistrib

This program computes efficiency-ordered distributions of magnetic fields as described in Boozer, Nuclear Fusion 55, 025001 (2015).  The main method is described on page 12, and some key definitions are given on page 9.

There are 3 surfaces, named 'plasma', 'middle', and 'outer', The 'plasma' surface can correspond to the outermost surface of a VMEC equilibrium, or it can be a plain circular toroidal surface. The name of the relevant surface appears as a suffix on
most variables.

On each surface we use a poloidal angle u and a toroidal angle v, defined as in NESCOIL.  The coordinate u lies in the range [0,1]. The stellarator has nfp identical toroidal periods, and an increase in v by 1 corresponds to 1 of these periods. Thus,
v increases by nfp in a complete toroidal revolution. In the output file, there is an array v corresponding to one toroidal period, as well as an array vl corresponding to all nfp toroidal periods.  Note that v is proportional to the standard
cylindrical angle phi: phi = 2*pi*v/nfp.  Generally, field lines are not straight in the (u,v) coordinates.  On the plasma surface, u is identical to the VMEC poloidal angle divided by 2*pi, while the VMEC toroidal angle differs by a sign relative to
v.

In the various variable names in the code and output file, 'r' refers to the position vector, not to a radius.  In various arrays with a dimension of length 3, this dimension always corresponds to Cartesian coordinates (x,y,z).

The 'normal' quantities in the code and output file refer to the surface normal vector N = (dr/dv) cross (dr/du) as in the NESCOIL paper. Note that this vector does not have unit magnitude.

*******************************************************************
Resolution parameters
*******************************************************************

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

Parameter: nu_middle
Type: int
Default: 16
When it matters: always
Description: Number of grid points in the poloidal direction on the middle surface

Parameter: nv_middle
Type: int
Default: 32
When it matters: always
Description: Number of grid points in the toroidal direction (per period) on the middle surface

Parameter: nu_outer
Type: int
Default: 16
When it matters: always
Description: Number of grid points in the poloidal direction on the outer surface

Parameter: nv_outer
Type: int
Default: 32
When it matters: always
Description: Number of grid points in the toroidal direction (per period) on the outer surface

*******************************************************************
Parameters giving the geometry of the 3 surfaces
*******************************************************************

Parameter: surface_option_plasma
Type: int
Default: 0
When it matters: always
Description:
 0 - The plasma surface will be a plain circular torus. The major radius will be R0_plasma.
     The minor radius will be a_plasma.
 1 - Same as option 0.
 2 - The plasma surface will be the last surface in the VMEC file specified by woutFilename.

Parameter: R0_plasma
Type: real
Default: 5.5
When it matter: only if surface_option_plasma = 0 or 1.
Description: Major radius of the plasma surface, in meters.

Parameter: a_plasma
Type: real
Default: 0.5
When it matter: only if surface_option_plasma = 0 or 1.
Description: Minor radius of the plasma surface, in meters.

Parameter: surface_option_middle
Type: int
Default: 0
When it matters: always
Description:
 0 - The middle surface will be a torus. The major radius will be Rmajor_p from the wout file.
     The minor radius will be a_middle.
 1 - Same as option 0, except the major radius will be R_middle.
 2 - The middle surface will computing by expanding the plasma LCFS uniformly by a distance separation_middle.

Parameter: surface_option_outer
Type: int
Default: 0
When it matters: always
Description:
 0 - The outer surface will be a torus. The major radius will be Rmajor_p from the wout file.
     The minor radius will be a_outer.
 1 - Same as option 0, except the major radius will be R_outer.
 2 - The outer surface will computing by expanding the plasma LCFS uniformly by a distance separation_outer.


*******************************************************************
Other parameters
*******************************************************************

Parameter: save_level
Type: int
Default: 1
When it matters:
Description: How much information to save in the .nc output file.
  0: Save everything
  1: Save everything but the interpolation matrices


