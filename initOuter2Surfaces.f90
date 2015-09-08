subroutine initOuter2Surfaces

  use globalVariables
  use initSurfaceMod

  implicit none

  integer :: tic, toc, countrate

  call system_clock(tic,countrate)
  print *,"Initializing middle surface."
  call initSurface(nu_middle, nv_middle, nvl_middle, u_middle, v_middle, vl_middle, &
       r_middle, drdu_middle, drdv_middle, normal_middle, norm_normal_middle, area_middle, &
       geometry_option_middle, R0_middle, a_middle, separation_middle, du_middle, dv_middle, &
       nescin_filename_middle)
  call system_clock(toc)
  print *,"Done initializing middle surface. Took ",real(toc-tic)/countrate," sec."

  call system_clock(tic)
  print *,"Initializing outer surface."
  call initSurface(nu_outer, nv_outer, nvl_outer, u_outer, v_outer, vl_outer, &
       r_outer, drdu_outer, drdv_outer, normal_outer, norm_normal_outer, area_outer, &
       geometry_option_outer, R0_outer, a_outer, separation_outer, du_outer, dv_outer, &
       nescin_filename_outer)
  call system_clock(toc)
  print *,"Done initializing outer surface. Took ",real(toc-tic)/countrate," sec."

end subroutine initOuter2Surfaces
