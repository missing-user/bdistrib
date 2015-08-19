subroutine buildInductanceMatrices()

  use buildInductanceMatrixMod
  use globalVariables

  implicit none

  integer :: tic, toc, countrate

  call system_clock(tic,countrate)
  print *,"Building inductance matrix between plasma surface and outer surface."
  call buildInductanceMatrix(inductance_plasma, r_plasma, normal_plasma, nu_plasma, nv_plasma, &
       mnmax_plasma, xm_plasma, xn_plasma, u_plasma, v_plasma)
  call system_clock(toc)
  print *,"Done building inductance matrix. Took ",real(toc-tic)/countrate," sec."

  call system_clock(tic)
  print *,"Building inductance matrix between middle surface and outer surface."
  call buildInductanceMatrix(inductance_middle, r_middle, normal_middle, nu_middle, nv_middle, &
       mnmax_middle, xm_middle, xn_middle, u_middle, v_middle)
  call system_clock(toc)
  print *,"Done building inductance matrix. Took ",real(toc-tic)/countrate," sec."

end subroutine buildInductanceMatrices
