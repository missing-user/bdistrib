module initSurfaceMod

  ! Passing un-allocated arrays is valid in modules but not in standalone subroutine
  ! files unless using pointers or an explicit interface.

  implicit none

  contains

    subroutine initSurface(nu, nv, nvl, u, v, vl, &
         x, y, z, &
         normal_x, normal_y, normal_z, &
         surface_option, R_specified, a, separation)

      use read_wout_mod
      use stel_kinds
      use stel_constants

      implicit none

      integer :: nu, nv, nvl, surface_option
      real(rprec) :: R_specified, a, separation
      real(rprec), dimension(:), allocatable :: u, v, vl
      real(rprec), dimension(:,:), allocatable :: x, y, z
      real(rprec), dimension(:,:), allocatable :: normal_x, normal_y, normal_z
      real(rprec), dimension(:,:), allocatable :: dxdu, dxdv, dydu, dydv, dzdu, dzdv
      real(rprec) :: R0_to_use, R
      real(rprec) :: angle, sinangle, cosangle, dsinangledu, dcosangledu
      real(rprec) :: angle2, sinangle2, cosangle2, dsinangle2dv, dcosangle2dv
      integer :: i, iu, iv, fzeroFlag
      real(rprec) :: u_rootSolve, rootSolve_abserr, rootSolve_relerr, v_rootSolve_min, v_rootSolve_max
      real(rprec) :: v_rootSolve_target, v_plasma_rootSolveSolution, x_new, y_new, z_new

      allocate(u(nu))
      allocate(v(nv))
      allocate(vl(nvl))

      do i = 1,nu
         u(i) = (i-1.0_rprec)/nu
      end do

      do i = 1,nv
         v(i) = (i-1.0_rprec)/nv
      end do

      do i = 1,nvl
         vl(i) = (i-1.0_rprec)/nv
      end do

      allocate(x(nu,nvl))
      allocate(y(nu,nvl))
      allocate(z(nu,nvl))
      allocate(dxdu(nu,nvl))
      allocate(dydu(nu,nvl))
      allocate(dzdu(nu,nvl))
      allocate(dxdv(nu,nvl))
      allocate(dydv(nu,nvl))
      allocate(dzdv(nu,nvl))
      allocate(normal_x(nu,nvl))
      allocate(normal_y(nu,nvl))
      allocate(normal_z(nu,nvl))

      x = 0
      y = 0
      z = 0
      dxdu = 0
      dydu = 0
      dzdu = 0
      dxdv = 0
      dydv = 0
      dzdv = 0

      select case (surface_option)
      case (0,1)
         ! Torus with circular cross-section

         if (surface_option==0) then
            R0_to_use = Rmajor
         else
            R0_to_use = R_specified
         end if

         do iu = 1,nu
            angle = twopi*u(iu)
            sinangle = sin(angle)
            cosangle = cos(angle)
            dsinangledu = cosangle*twopi
            dcosangledu = -sinangle*twopi
            do iv = 1,nvl
               angle2 = twopi * vl(iv) / nfp
               sinangle2 = sin(angle2)
               cosangle2 = cos(angle2)
               dsinangle2dv = cosangle2*twopi/nfp
               dcosangle2dv = -sinangle2*twopi/nfp

               x(iu,iv) = (R0_to_use + a * cosangle) * cosangle2
               y(iu,iv) = (R0_to_use + a * cosangle) * sinangle2
               z(iu,iv) = a * sinangle

               dxdu(iu,iv) = (a * dcosangledu) * cosangle2
               dydu(iu,iv) = (a * dcosangledu) * sinangle2
               dzdu(iu,iv) = a * dsinangledu

               dxdv(iu,iv) = (R0_to_use + a * cosangle) * dcosangle2dv
               dydv(iu,iv) = (R0_to_use + a * cosangle) * dsinangle2dv
               !dzdv(iu,iv) = 0, so no equation needed for it here.
            end do
         end do

      case (2)

         !rootSolve_abserr = 0
         !rootSolve_relerr = 0
         rootSolve_abserr = 1.0e-10_rprec
         rootSolve_relerr = 1.0e-10_rprec
         do iu = 1,nu
            u_rootSolve = u(iu)
            print *,"u=",u_rootSolve
            do iv = 1,nvl
               v_rootSolve_target = vl(iv)
               v_rootSolve_min = v_rootSolve_target - 0.3
               v_rootSolve_max = v_rootSolve_target + 0.3
               print *,"***",fzero_residual(v_rootSolve_min), fzero_residual(v_rootSolve_target), fzero_residual(v_rootSolve_max)
               call fzero(fzero_residual, v_rootSolve_min, v_rootSolve_max, v_rootSolve_target, &
                    rootSolve_relerr, rootSolve_abserr, fzeroFlag)
               ! Note: fzero returns its answer in v_rootSolve_min
               v_plasma_rootSolveSolution = v_rootSolve_min
               if (fzeroFlag == 4) then
                  stop "ERROR: fzero returned error 4: no sign change in residual"
               else if (fzeroFlag > 1) then
                  print *,"WARNING: fzero returned an error code:",fzeroFlag
               end if

               call expandPlasmaSurface(u_rootSolve, v_plasma_rootSolveSolution, separation, x_new, y_new, z_new)
               x(iu,iv) = x_new
               y(iu,iv) = y_new
               z(iu,iv) = z_new

            end do

         end do

      case default
         print *,"Invalid setting for surface_option: ",surface_option
         stop
      end select

      ! Evaluate cross product:
      normal_x = dydv * dzdu - dydu * dzdv
      normal_y = dzdv * dxdu - dzdu * dxdv
      normal_z = dxdv * dydu - dxdu * dydv

      contains

        function fzero_residual(v_plasma_test)

          implicit none

          real(rprec) :: v_plasma_test, fzero_residual
          real(rprec) :: x_outer, y_outer, z_outer, v_outer_new, v_error

          call expandPlasmaSurface(u_rootSolve, v_plasma_test, separation, x_outer, y_outer, z_outer)
          v_outer_new = atan2(y_outer,x_outer)*nfp/twopi
          v_error = v_outer_new - v_rootSolve_target
          if (v_error < -nfp/2.0_rprec) then
             v_error = v_error + nfp
          end if
          if (v_error > nfp/2.0_rprec) then
             v_error = v_error - nfp
          end if
          fzero_residual = v_error

        end function fzero_residual

    end subroutine initSurface

  end module initSurfaceMod
  
