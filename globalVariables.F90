module globalVariables

  use stel_kinds

  implicit none

  integer :: nu_plasma=16, nv_plasma=18, nvl_plasma
  integer :: nu_middle=20, nv_middle=22, nvl_middle
  integer :: nu_current=24, nv_current=26, nvl_current

  integer :: surface_option_middle=0
  integer :: surface_option_current=0

  real(dp) :: R_middle=5.5, R_current=5.5
  real(dp) :: a_middle=0.3, a_current=0.6
  real(dp) :: separation_middle=0.3, separation_current=0.6

  character(len=200) :: woutFilename=""
  character(len=200) :: outputFilename

  real(dp) :: norm_normal_plasma, norm_normal_middle, norm_normal_current

  real(dp), dimension(:), allocatable :: u_plasma, v_plasma, vl_plasma
  real(dp), dimension(:,:), allocatable :: x_plasma, y_plasma, z_plasma
  real(dp), dimension(:,:), allocatable :: dxdu_plasma, dydu_plasma, dzdu_plasma
  real(dp), dimension(:,:), allocatable :: dxdv_plasma, dydv_plasma, dzdv_plasma
  real(dp), dimension(:,:), allocatable :: normal_x_plasma, normal_y_plasma, normal_z_plasma

  real(dp), dimension(:), allocatable :: u_middle, v_middle, vl_middle
  real(dp), dimension(:,:), allocatable :: x_middle, y_middle, z_middle
  real(dp), dimension(:,:), allocatable :: dxdu_middle, dydu_middle, dzdu_middle
  real(dp), dimension(:,:), allocatable :: dxdv_middle, dydv_middle, dzdv_middle
  real(dp), dimension(:,:), allocatable :: normal_x_middle, normal_y_middle, normal_z_middle

  real(dp), dimension(:), allocatable :: u_current, v_current, vl_current
  real(dp), dimension(:,:), allocatable :: x_current, y_current, z_current
  real(dp), dimension(:,:), allocatable :: dxdu_current, dydu_current, dzdu_current
  real(dp), dimension(:,:), allocatable :: dxdv_current, dydv_current, dzdv_current
  real(dp), dimension(:,:), allocatable :: normal_x_current, normal_y_current, normal_z_current

end module globalVariables

