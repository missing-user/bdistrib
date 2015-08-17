module globalVariables

  use stel_kinds

  implicit none

  integer :: nu_plasma=16, nv_plasma=18, nvl_plasma
  integer :: nu_middle=20, nv_middle=22, nvl_middle
  integer :: nu_outer=24, nv_outer=26, nvl_outer

  integer :: surface_option_plasma  = 2
  integer :: surface_option_middle  = 0
  integer :: surface_option_outer = 0

  real(dp) :: R0_plasma = 5.5, R0_middle = 5.5, R0_outer = 5.5
  real(dp) :: a_plasma = 1.0, a_middle = 1.5, a_outer = 2.0
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

  integer :: currentPotential_mnmax, currentPotential_mpol=4, currentPotential_ntor=4
  integer, dimension(:), allocatable :: currentPotential_xm, currentPotential_xn

  real(dp), dimension(:,:), allocatable :: inductance_plasma, inductance_middle

  integer :: n_singular_values_inductance_plasma, n_singular_values_inductance_middle
  real(dp), dimension(:), allocatable :: svd_s_inductance_plasma, svd_s_inductance_middle

  real(dp), dimension(:,:), allocatable :: svd_uT_inductance_middle, svd_v_inductance_middle

  integer, parameter :: nmax_pseudoinverse_thresholds = 1000
  integer :: n_pseudoinverse_thresholds
  real(dp) :: pseudoinverse_thresholds(nmax_pseudoinverse_thresholds)

  integer :: save_level = 1
  integer :: n_singular_vectors_to_save = 5, n_singular_values_transferMatrix
  integer, dimension(:), allocatable :: n_singular_values_retained
  real(dp), dimension(:,:), allocatable :: svd_s_transferMatrix
  real(dp), dimension(:,:,:), allocatable :: svd_u_transferMatrix, svd_v_transferMatrix

end module globalVariables

