subroutine buildInductanceMatrices()

  use buildInductanceMatrixMod
  use globalVariables

  implicit none

  integer :: tic, toc, countrate

  select case (basis_set_option)
  case (1,2)
     num_basis_functions_plasma = mnmax_plasma
     num_basis_functions_middle = mnmax_middle
     num_basis_functions_outer  = mnmax_outer
  case (3)
     num_basis_functions_plasma = mnmax_plasma * 2
     num_basis_functions_middle = mnmax_middle * 2
     num_basis_functions_outer  = mnmax_outer  * 2
  case default
     print *,"Error! Invalid setting for basis_set_option:",basis_set_option
     stop
  end select

  print *,"Number of basis functions on plasma surface:",num_basis_functions_plasma
  print *,"Number of basis functions on middle surface:",num_basis_functions_middle
  print *,"Number of basis functions on outer surface: ",num_basis_functions_outer

  n_singular_vectors_to_save = min(n_singular_vectors_to_save, &
       num_basis_functions_plasma, num_basis_functions_middle)

  call system_clock(tic,countrate)
  print *,"Building inductance matrix between plasma surface and outer surface."
  call buildInductanceMatrix(inductance_plasma, r_plasma, normal_plasma, nu_plasma, nv_plasma, &
       mnmax_plasma, num_basis_functions_plasma, xm_plasma, xn_plasma, u_plasma, v_plasma)
  call system_clock(toc)
  print *,"Done building inductance matrix. Took ",real(toc-tic)/countrate," sec."

  call system_clock(tic)
  print *,"Building inductance matrix between middle surface and outer surface."
  call buildInductanceMatrix(inductance_middle, r_middle, normal_middle, nu_middle, nv_middle, &
       mnmax_middle, num_basis_functions_middle, xm_middle, xn_middle, u_middle, v_middle)
  call system_clock(toc)
  print *,"Done building inductance matrix. Took ",real(toc-tic)/countrate," sec."

end subroutine buildInductanceMatrices
