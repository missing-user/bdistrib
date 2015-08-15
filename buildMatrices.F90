subroutine buildMatrices()

  use buildInductanceMatrixMod
  use globalVariables

  implicit none

  print *,"Building inductance matrix between plasma surface and current surface."
  call buildInductanceMatrix(inductance_plasma, r_plasma, normal_plasma, nu_plasma, nv_plasma)
  print *,"Done."

  print *,"Building inductance matrix between plasma surface and current surface."
  call buildInductanceMatrix(inductance_middle, r_middle, normal_middle, nu_middle, nv_middle)
  print *,"Done."

end subroutine buildMatrices
