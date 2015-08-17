subroutine initOuterSurfaces

  use globalVariables
  use initSurfaceMod

  implicit none

  integer :: tic, toc, countrate

  call system_clock(tic,countrate)
  print *,"Initializing middle surface."
  call initSurface(nu_middle, nv_middle, nvl_middle, u_middle, v_middle, vl_middle, &
       r_middle, drdu_middle, drdv_middle, normal_middle, norm_normal_middle, &
       surface_option_middle, R0_middle, a_middle, separation_middle)
  call system_clock(toc)
  print *,"Done initializing middle surface. Took ",real(toc-tic)/countrate," sec."

  call system_clock(tic)
  print *,"Initializing current surface."
  call initSurface(nu_current, nv_current, nvl_current, u_current, v_current, vl_current, &
       r_current, drdu_current, drdv_current, normal_current, norm_normal_current, &
       surface_option_current, R0_current, a_current, separation_current)
  call system_clock(toc)
  print *,"Done initializing current surface. Took ",real(toc-tic)/countrate," sec."

end subroutine initOuterSurfaces
