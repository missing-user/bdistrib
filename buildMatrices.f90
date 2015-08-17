subroutine buildMatrices()

  use buildInductanceMatrixMod
  use globalVariables

  implicit none

  integer :: tic, toc, countrate

  call system_clock(tic,countrate)
  print *,"Building inductance matrix between plasma surface and current surface."
  call buildInductanceMatrix(inductance_plasma, r_plasma, normal_plasma, nu_plasma, nv_plasma)
  call system_clock(toc)
  print *,"Done building inductance matrix. Took ",real(toc-tic)/countrate," sec."

  call system_clock(tic)
  print *,"Building inductance matrix between plasma surface and current surface."
  call buildInductanceMatrix(inductance_middle, r_middle, normal_middle, nu_middle, nv_middle)
  call system_clock(toc)
  print *,"Done building inductance matrix. Took ",real(toc-tic)/countrate," sec."

end subroutine buildMatrices
