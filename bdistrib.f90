! Main program

program bdistrib

  implicit none

  print *,"This is BDISTRIB."

  call readInput()
  call validateInput()

  ! Define the position vector and normal vector at each grid point for the 3 surfaces:
  call initPlasma()
  call initOuter2Surfaces()

  ! This next subroutine is not used yet.
  call initCurrentPotentialModes()

  ! Compute the mutual inductance matrices, which relate current on the outer surface to B_n on the inner surfaces:
  call buildInductanceMatrices()

  ! Compute SVD of each of the inductance matrices:
  call svdInductanceMatrices()

  call buildTransferMatrix()

  call writeOutput()

  print *,"BDISTRIB complete."

end program bdistrib
