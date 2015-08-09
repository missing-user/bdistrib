module initSurfaceMod

  ! Passing un-allocated arrays is valid in modules but not in standalone subroutine
  ! files unless using pointers or an explicit interface.

  implicit none
  
  contains

    subroutine initSurface(nu, nv, u, v, &
         x, y, z, &
         surface_option, R, a, separation)

      use read_wout_mod
      use stel_kinds
      
      implicit none

      integer :: nu, nv, surface_option
      real(rprec) :: R, a, separation
      real(rprec), dimension(:), allocatable :: u, v
      real(rprec), dimension(:,:), allocatable :: x, y, z
      integer :: i

      allocate(u(nu))
      allocate(v(nv))

      do i = 1,nu
         u(i) = (i-1.0_rprec)/nu
      end do

      do i = 1,nv
         v(i) = (i-1.0_rprec)/nv
      end do

    end subroutine initSurface

  end module initSurfaceMod
  
