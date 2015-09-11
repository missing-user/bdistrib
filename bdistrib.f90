! Main program

program bdistrib

  use globalVariables
  use init_basis_functions_mod

  implicit none

  integer :: tic, toc, countrate, tic1, toc1

  print *,"This is BDISTRIB."
  call system_clock(tic,countrate)

  call readInput()
  call validateInput()

  ! Define the position vector and normal vector at each grid point for the 3 surfaces:
  call initPlasma()
  call initOuter2Surfaces()

  ! Initialize basis functions on the 3 surfaces:
  call system_clock(tic1)
  print *,"Initializing basis functions on the plasma surface."
  call init_basis_functions(basis_option_plasma, mpol_plasma, ntor_plasma, mnmax_plasma, xm_plasma, xn_plasma, &
       u_plasma, v_plasma, num_basis_functions_plasma, basis_functions_plasma, area_plasma, norm_normal_plasma, &
       should_be_identity_plasma)
  call system_clock(toc1)
  print *,"Done. Took ",real(toc1-tic1)/countrate," sec."

  call system_clock(tic1)
  print *,"Initializing basis functions on the middle surface."
  call init_basis_functions(basis_option_middle, mpol_middle, ntor_middle, mnmax_middle, xm_middle, xn_middle, &
       u_middle, v_middle, num_basis_functions_middle, basis_functions_middle, area_middle, norm_normal_middle, &
       should_be_identity_middle)
  call system_clock(toc1)
  print *,"Done. Took ",real(toc1-tic1)/countrate," sec."

  call system_clock(tic1)
  print *,"Initializing basis functions on the outer surface."
  call init_basis_functions(basis_option_outer, mpol_outer, ntor_outer, mnmax_outer, xm_outer, xn_outer, &
       u_outer, v_outer, num_basis_functions_outer, basis_functions_outer, area_outer, norm_normal_outer, &
       should_be_identity_outer)
  call system_clock(toc1)
  print *,"Done. Took ",real(toc1-tic1)/countrate," sec."


  ! Compute the mutual inductance matrices, which relate current on 1 surface to B_normal on an interior surface:
  call buildInductanceMatrices()

  ! Compute SVD of each of the inductance matrices:
  call svdInductanceMatrices()

  ! Assemble transfer matrix and compute its SVD:
  call transfer_matrix()

  call system_clock(toc)
  totalTime = real(toc-tic)/countrate

  call writeOutput()

  if (allSVDsSucceeded) then
     print *,"All SVDs succeeded."
  else
     print *,"**************************"
     print *,"At least one SVD failed!!!"
     print *,"**************************"
  end if

  print *,"BDISTRIB complete. Total time=",totalTime,"sec."
  print *,"You can run bdistribPlot ",trim(outputFilename)," to plot results."

end program bdistrib
