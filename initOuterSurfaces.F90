subroutine initOuterSurfaces

  use globalVariables
  use initSurfaceMod

  implicit none

  call initSurface(nu_middle, nv_middle, nvl_middle, u_middle, v_middle, vl_middle, &
       x_middle, y_middle, z_middle, &
       normal_x_middle, normal_y_middle, normal_z_middle, &
       surface_option_middle, R_middle, a_middle, separation_middle)

  call initSurface(nu_current, nv_current, nvl_current, u_current, v_current, vl_current, &
       x_current, y_current, z_current, &
       normal_x_current, normal_y_current, normal_z_current, &
       surface_option_current, R_current, a_current, separation_current)

end subroutine initOuterSurfaces
