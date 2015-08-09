subroutine initPlasma

  use globalVariables
  use read_wout_mod

  implicit none

  integer :: i

  allocate(u_plasma(nu_plasma))
  allocate(v_plasma(nv_plasma))

  do i=1,nu_plasma
     u_plasma(i) = (i-1.0_rprec)/nu_plasma
  end do

  do i=1,nv_plasma
     v_plasma(i) = (i-1.0_rprec)/nv_plasma
  end do

  

end subroutine initPlasma
