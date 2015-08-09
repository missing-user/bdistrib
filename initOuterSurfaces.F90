subroutine initOuterSurfaces

  use globalVariables
  use initSurfaceMod

  implicit none

  print *,"AAA"
  call initSurface(nu_middle, nv_middle, u_middle, v_middle, &
       x_middle, y_middle, z_middle, &
       surface_option_middle, R_middle, a_middle, separation_middle)

  print *,"BBB"
  call initSurface(nu_current, nv_current, u_current, v_current, &
       x_current, y_current, z_current, &
       surface_option_current, R_current, a_current, separation_current)
  print *,"CCC"

end subroutine initOuterSurfaces
