module initFourierModesMod

  implicit none

contains

  subroutine initFourierModes(mpol, ntor, mnmax, xm, xn)

    implicit none

    integer :: mpol, ntor, mnmax
    integer, dimension(:), allocatable :: xm, xn
    
    integer :: jn, jm, index
    
    ! xm is nonnegative.
    ! xn can be negative, zero, or positive.
    ! When xm is 0, xn must be positive.
    mnmax = mpol*(ntor*2+1) + ntor
    
    allocate(xm(mnmax))
    allocate(xn(mnmax))
    
    ! Handle the xm=0 modes:
    xm=0
    do jn=1,ntor
       xn(jn)=jn
    end do
    
    ! Handle the xm>0 modes:
    index = ntor
    do jm = 1,mpol
       do jn = -ntor, ntor
          index = index + 1
          xn(index) = jn
          xm(index) = jm
       end do
    end do
    
    if (index .ne. mnmax) then
       print *,"Error!  index=",index," but mnmax=",mnmax
       stop
    end if
    
  end subroutine initFourierModes
  
end module initFourierModesMod
