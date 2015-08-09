module globalVariables

  use stel_kinds

  implicit none

  integer :: nu_plasma=16, nv_plasma=18
  integer :: nu_middle=20, nv_middle=22
  integer :: nu_current=24, nv_current=26

  integer :: surface_option_middle=0
  integer :: surface_option_current=0

  real(dp) :: R_middle=5.5, R_current=5.5
  real(dp) :: a_middle=0.3, a_current=0.6
  real(dp) :: separation_middle=0.3, separation_current=0.6

  character(len=200) :: woutFilename=""
  character(len=200) :: outputFilename

  real(dp), dimension(:), allocatable :: u_plasma, v_plasma
  real(dp), dimension(:,:), allocatable :: x_plasma, y_plasma, z_plasma

  real(dp), dimension(:), allocatable :: u_middle, v_middle
  real(dp), dimension(:,:), allocatable :: x_middle, y_middle, z_middle

  real(dp), dimension(:), allocatable :: u_current, v_current
  real(dp), dimension(:,:), allocatable :: x_current, y_current, z_current

end module globalVariables

