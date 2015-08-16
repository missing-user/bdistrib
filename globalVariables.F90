module globalVariables

  use stel_kinds

  implicit none

  integer :: nu_plasma=16, nv_plasma=18, nvl_plasma
  integer :: nu_middle=20, nv_middle=22, nvl_middle
  integer :: nu_current=24, nv_current=26, nvl_current

  integer :: surface_option_plasma  = 2
  integer :: surface_option_middle  = 0
  integer :: surface_option_current = 0

  real(dp) :: R0_middle=5.5, R0_current=5.5
  real(dp) :: a_middle=0.3, a_current=0.6
  real(dp) :: separation_middle=0.3, separation_current=0.6

  character(len=200) :: woutFilename=""
  character(len=200) :: outputFilename

  real(dp), dimension(:), allocatable :: u_plasma, v_plasma, vl_plasma
  real(dp), dimension(:,:,:), allocatable :: r_plasma, drdu_plasma, drdv_plasma, normal_plasma

  real(dp), dimension(:), allocatable :: u_middle, v_middle, vl_middle
  real(dp), dimension(:,:,:), allocatable :: r_middle, drdu_middle, drdv_middle, normal_middle

  real(dp), dimension(:), allocatable :: u_current, v_current, vl_current
  real(dp), dimension(:,:,:), allocatable :: r_current, drdu_current, drdv_current, normal_current

  real(dp), dimension(:,:), allocatable :: norm_normal_plasma, norm_normal_middle, norm_normal_current

  integer :: currentPotential_mnmax, currentPotential_mpol=4, currentPotential_ntor=4
  integer, dimension(:), allocatable :: currentPotential_xm, currentPotential_xn

  real(dp), dimension(:,:), allocatable :: inductance_plasma, inductance_middle

  integer :: n_singular_values_inductance_plasma, n_singular_values_inductance_middle
  real(dp), dimension(:), allocatable :: svd_s_inductance_plasma, svd_s_inductance_middle

end module globalVariables

