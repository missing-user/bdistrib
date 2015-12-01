subroutine one_over_R_field

  use global_variables
  use stel_kinds

  implicit none

  integer :: iu, iv, index, iflag
  real(dp) :: R2, temp, factors, du, dv
  integer :: tic, toc, countrate
  real(dp), dimension(:), allocatable :: tempVec

  call system_clock(tic, countrate)
  print *,"Computing quantities related to the 1/R field."
  
  allocate(normal_component_of_1_over_R_field(num_basis_functions_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(normal_component_of_1_over_R_field_uv(nu_plasma,nv_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(tempVec(nu_plasma*nv_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  du = u_plasma(2)-u_plasma(1)
  dv = v_plasma(2)-v_plasma(1)
  factors = nfp * du * dv
  do iv = 1,nv_plasma
     do iu = 1,nu_plasma
        index = (iv-1)*nu_plasma + iu
        R2 = r_plasma(1,iu,iv)*r_plasma(1,iu,iv) + r_plasma(2,iu,iv)*r_plasma(2,iu,iv)
        temp = (r_plasma(1,iu,iv)*normal_plasma(2,iu,iv)-r_plasma(2,iu,iv)*normal_plasma(1,iu,iv)) / R2
        normal_component_of_1_over_R_field_uv(iu,iv) = temp / norm_normal_plasma(iu,iv)
        tempVec(index) = temp*factors
     end do
  end do

  ! This next line could be optimized using BLAS
  normal_component_of_1_over_R_field = matmul(tempVec,basis_functions_plasma)
  !normal_component_of_1_over_R_field = matmul(transpose(basis_functions_plasma),tempVec)

  deallocate(tempVec)
  call system_clock(toc)
  print *,"Done with 1/R field computations. Took ",real(toc-tic)/countrate," sec."

end subroutine one_over_R_field
