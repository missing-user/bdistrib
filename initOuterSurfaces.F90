subroutine initOuterSurfaces

  use globalVariables
  use initSurfaceMod

  implicit none

  call initSurface(nu_middle, nv_middle, nvl_middle, u_middle, v_middle, vl_middle, &
       r_middle, drdu_middle, drdv_middle, normal_middle, norm_normal_middle, &
       surface_option_middle, R0_middle, a_middle, separation_middle)

  call initSurface(nu_current, nv_current, nvl_current, u_current, v_current, vl_current, &
       r_current, drdu_current, drdv_current, normal_current, norm_normal_current, &
       surface_option_current, R0_current, a_current, separation_current)

end subroutine initOuterSurfaces
