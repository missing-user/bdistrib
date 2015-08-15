! Main program

program bdistrib

  implicit none

  print *,"This is BDISTRIB."

  call readInput()

  ! Define the position vector and normal vector at each grid point for the 3 surfaces:
  call initPlasma()
  call initOuterSurfaces()
  call initCurrentPotentialModes()

  ! Compute the matrices relating current on the outer surface to B_n on the inner surfaces
  call buildMatrices()

  call writeOutput()

  print *,"BDISTRIB complete."

end program bdistrib
