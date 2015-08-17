! Main program

program bdistrib

  use svdMod

  implicit none

  print *,"This is BDISTRIB."

  call readInput()

  ! Define the position vector and normal vector at each grid point for the 3 surfaces:
  call initPlasma()
  call initOuterSurfaces()
  call initCurrentPotentialModes()

  ! Compute the mutual inductance matrices, which relate current on the outer surface to B_n on the inner surfaces:
  call buildMatrices()

  ! Compute SVD of each of the inductance matrices:
  call svd1()

  call writeOutput()

  print *,"BDISTRIB complete."

end program bdistrib
