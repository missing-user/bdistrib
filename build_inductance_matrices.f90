subroutine build_inductance_matrices()

  use build_inductance_matrix_mod
  use global_variables

  implicit none

  integer :: tic, toc, countrate

  print *,"Number of basis functions on plasma surface:",num_basis_functions_plasma
  print *,"Number of basis functions on middle surface:",num_basis_functions_middle
  print *,"Number of basis functions on outer surface: ",num_basis_functions_outer

  n_singular_vectors_to_save = min(n_singular_vectors_to_save, &
       num_basis_functions_plasma, num_basis_functions_middle)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  call system_clock(tic,countrate)
  print *,"Building inductance matrix between plasma surface and outer surface."

  call build_inductance_matrix(inductance_plasma_outer, u_plasma, v_plasma, r_plasma, normal_plasma, basis_functions_plasma, &
       u_outer, v_outer, r_outer, normal_outer, basis_functions_outer)

  call system_clock(toc)
  print *,"Done building inductance matrix. Took ",real(toc-tic)/countrate," sec."

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  call system_clock(tic,countrate)
  print *,"Building inductance matrix between middle surface and outer surface."

  call build_inductance_matrix(inductance_middle_outer, u_middle, v_middle, r_middle, normal_middle, basis_functions_middle, &
       u_outer, v_outer, r_outer, normal_outer, basis_functions_outer)

  call system_clock(toc)
  print *,"Done building inductance matrix. Took ",real(toc-tic)/countrate," sec."

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  call system_clock(tic,countrate)
  print *,"Building inductance matrix between plasma surface and middle surface."

  call build_inductance_matrix(inductance_plasma_middle, u_plasma, v_plasma, r_plasma, normal_plasma, basis_functions_plasma, &
       u_middle, v_middle, r_middle, normal_middle, basis_functions_middle)

  call system_clock(toc)
  print *,"Done building inductance matrix. Took ",real(toc-tic)/countrate," sec."

end subroutine build_inductance_matrices
