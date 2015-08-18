! Main program

program bdistrib

  use globalVariables
  use initFourierModesMod

  implicit none

  integer :: tic, toc, countrate

  print *,"This is BDISTRIB."
  call system_clock(tic,countrate)

  call readInput()
  call validateInput()

  ! Define the position vector and normal vector at each grid point for the 3 surfaces:
  call initPlasma()
  call initOuter2Surfaces()

  call initFourierModes(mpol_plasma, ntor_plasma, mnmax_plasma, xm_plasma, xn_plasma)
  call initFourierModes(mpol_middle, ntor_middle, mnmax_middle, xm_middle, xn_middle)
  call initFourierModes(mpol_outer, ntor_outer, mnmax_outer, xm_outer, xn_outer)

  ! Compute the mutual inductance matrices, which relate current on the outer surface to B_n on the inner surfaces:
  call buildInductanceMatrices()

  ! Compute SVD of each of the inductance matrices:
  call svdInductanceMatrices()

  call buildTransferMatrix()

  call writeOutput()

  if (allSVDsSucceeded) then
     print *,"All SVDs succeeded."
  else
     print *,"**************************"
     print *,"At least one SVD failed!!!"
     print *,"**************************"
  end if

  call system_clock(toc)
  print *,"BDISTRIB complete. Total time=",real(toc-tic)/countrate,"sec."

end program bdistrib
