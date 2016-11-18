subroutine read_bnorm()

  use global_variables, only: load_bnorm, bnorm_filename, num_basis_functions_plasma, nu_plasma, nv_plasma, &
       Bnormal_from_plasma_current, Bnormal_from_plasma_current_inductance, &
       Bnormal_from_plasma_current_transfer, Bnormal_from_plasma_current_uv, &
       u_plasma, v_plasma, norm_normal_plasma, basis_functions_plasma, nfp, curpol
  use safe_open_mod
  use stel_constants
  use stel_kinds
  
  implicit none

  integer :: iunit, i, mm, nn, iu, iv, index, iflag
  real(dp) :: bf, du, dv
  real(dp), dimension(:), allocatable :: tempVec
  integer :: tic, toc, countrate

  call system_clock(tic,countrate)

  du = u_plasma(2)-u_plasma(1)
  dv = v_plasma(2)-v_plasma(1)
  allocate(Bnormal_from_plasma_current(num_basis_functions_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Bnormal_from_plasma_current_inductance(num_basis_functions_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Bnormal_from_plasma_current_transfer(num_basis_functions_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Bnormal_from_plasma_current_uv(nu_plasma,nv_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  Bnormal_from_plasma_current_uv = 0

  if (.not. load_bnorm) then
     print *,"Not reading a bnorm file, so Bnormal_from_plasma_current arrays will all be 0."
     Bnormal_from_plasma_current = 0
     Bnormal_from_plasma_current_inductance = 0
     Bnormal_from_plasma_current_transfer = 0
     return
  end if

  print *,"Loading B_normal on the plasma surface due to plasma current from file ",trim(bnorm_filename)

  call safe_open(iunit, i, trim(bnorm_filename), 'old', 'formatted')
  if (i .ne. 0 ) then
     stop 'Unable to open bnorm_filename.'
  end if

  do
     read(iunit,*,iostat = i) mm, nn, bf
     if (i .ne. 0) exit
     print *,"  Adding a mode with m=",mm,",  n=",nn
     do iv = 1,nv_plasma
        do iu = 1,nu_plasma
           Bnormal_from_plasma_current_uv(iu,iv) = Bnormal_from_plasma_current_uv(iu,iv) + &
                bf*sin(2*pi*(mm*u_plasma(iu)+nn*v_plasma(iv)))
           ! To see that it should be (mu+nv) rather than (mu-nv) in the above line, you can examine
           ! BNORM/Sources/bn_fouri.f (where the arrays in the bnorm files are computed)
           ! or
           ! either NESCOIL/Sources/bnfld.f (where bnorm files are read)
        end do
     end do
  end do

  close(iunit)

  ! BNORM scales B_n by curpol=(2*pi/nfp)*bsubv(m=0,n=0)
  ! where bsubv is the extrapolation to the last full mesh point of
  ! bsubvmnc.  Let's undo this scaling now.
  Bnormal_from_plasma_current_uv = Bnormal_from_plasma_current_uv * curpol

  ! Convert from (u,v) grid to a vector of fluxes:
 
  allocate(tempVec(nu_plasma*nv_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  do iv = 1,nv_plasma
     do iu = 1,nu_plasma
        index = (iv-1)*nu_plasma + iu
        tempVec(index) = Bnormal_from_plasma_current_uv(iu,iv) * norm_normal_plasma(iu,iv)
     end do
  end do
  tempVec = tempVec * (nfp * du * dv)
  ! This next line could be optimized using BLAS
  Bnormal_from_plasma_current = matmul(tempVec,basis_functions_plasma)

  deallocate(tempVec)

  call system_clock(toc)
  print *,"Done reading B_normal on the plasma surface due to plasma current."
  print *,"Took ",real(toc-tic)/countrate," sec."

end subroutine  read_bnorm
