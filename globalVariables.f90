module globalVariables

  use stel_kinds

  implicit none

  integer :: nu_plasma=64, nv_plasma=64, nvl_plasma
  integer :: nu_middle=64, nv_middle=64, nvl_middle
  integer :: nu_outer =64, nv_outer =64, nvl_outer

  integer :: geometry_option_plasma = 0
  integer :: geometry_option_middle = 0
  integer :: geometry_option_outer  = 0

  real(dp) :: R0_plasma = 10.0, R0_middle = 10.0, R0_outer = 10.0
  real(dp) :: a_plasma = 0.5, a_middle = 1.0, a_outer = 1.5
  real(dp) :: separation_middle=0.3, separation_outer=0.6

  character(len=200) :: woutFilename=""
  character(len=200) :: outputFilename

  real(dp), dimension(:), allocatable :: u_plasma, v_plasma, vl_plasma
  real(dp), dimension(:,:,:), allocatable :: r_plasma, drdu_plasma, drdv_plasma, normal_plasma

  real(dp), dimension(:), allocatable :: u_middle, v_middle, vl_middle
  real(dp), dimension(:,:,:), allocatable :: r_middle, drdu_middle, drdv_middle, normal_middle

  real(dp), dimension(:), allocatable :: u_outer, v_outer, vl_outer
  real(dp), dimension(:,:,:), allocatable :: r_outer, drdu_outer, drdv_outer, normal_outer

  real(dp), dimension(:,:), allocatable :: norm_normal_plasma, norm_normal_middle, norm_normal_outer

  real(dp) :: du_plasma, dv_plasma, du_middle, dv_middle, du_outer, dv_outer

  integer :: mpol_plasma=8, mpol_middle=8, mpol_outer=8
  integer :: ntor_plasma=8, ntor_middle=8, ntor_outer=8
  integer :: mnmax_plasma, mnmax_middle, mnmax_outer
  integer :: num_basis_functions_plasma, num_basis_functions_middle, num_basis_functions_outer
  integer, dimension(:), allocatable :: xm_plasma, xm_middle, xm_outer
  integer, dimension(:), allocatable :: xn_plasma, xn_middle, xn_outer

  real(dp), dimension(:,:), allocatable :: inductance_plasma, inductance_middle

  integer :: n_singular_values_inductance_plasma, n_singular_values_inductance_middle
  real(dp), dimension(:), allocatable :: svd_s_inductance_plasma, svd_s_inductance_middle

  real(dp), dimension(:,:), allocatable :: svd_uT_inductance_middle, svd_v_inductance_middle

  integer, parameter :: nmax_pseudoinverse_thresholds = 1000
  integer :: n_pseudoinverse_thresholds
  real(dp) :: pseudoinverse_thresholds(nmax_pseudoinverse_thresholds)

  integer :: save_level = 2
  integer :: n_singular_vectors_to_save = 5, n_singular_values_transferMatrix
  integer, dimension(:), allocatable :: n_singular_values_retained
  real(dp), dimension(:,:), allocatable :: svd_s_transferMatrix

  ! EZCDF doesn't allow 4D arrays, so I can't include sin/cos as a 4th dimension.
  real(dp), dimension(:,:,:), allocatable :: svd_u_transferMatrix_sin, svd_v_transferMatrix_sin
  real(dp), dimension(:,:,:), allocatable :: svd_u_transferMatrix_cos, svd_v_transferMatrix_cos

  logical :: allSVDsSucceeded
  integer :: nfp_imposed = 1

  integer :: basis_set_option = 1
  real(dp) :: totalTime

end module globalVariables

