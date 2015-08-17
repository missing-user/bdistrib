subroutine initCurrentPotentialModes()

  use globalVariables, only: currentPotential_mnmax, currentPotential_xm, currentPotential_xn, currentPotential_mpol, currentPotential_ntor

  implicit none

  integer :: jn, jm, index

  ! xm is nonnegative.
  ! xn can be negative, zero, or positive.
  ! When xm is 0, xn must be positive.
  currentPotential_mnmax = currentPotential_mpol*(currentPotential_ntor*2+1) + currentPotential_ntor

  allocate(currentPotential_xm(currentPotential_mnmax))
  allocate(currentPotential_xn(currentPotential_mnmax))

  ! Handle the xm=0 modes:
  currentPotential_xm=0
  do jn=1,currentPotential_ntor
     currentPotential_xn(jn)=jn
  end do

  ! Handle the xm>0 modes:
  index = currentPotential_ntor
  do jm = 1,currentPotential_mpol
     do jn = -currentPotential_ntor, currentPotential_ntor
        index = index + 1
        currentPotential_xn(index) = jn
        currentPotential_xm(index) = jm
     end do
  end do

  if (index .ne. currentPotential_mnmax) then
     print *,"Error!  index=",index," but currentPotential_mnmax=",currentPotential_mnmax
     stop
  end if

end subroutine initCurrentPotentialModes
