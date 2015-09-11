module init_surface_mod

  ! Passing un-allocated arrays is valid in modules but not in standalone subroutine
  ! files unless using pointers or an explicit interface.

  implicit none

  contains

    subroutine init_surface(nu, nv, nvl, u, v, vl, &
         r, drdu, drdv, normal, norm_normal, area, &
         geometry_option, R_specified, a, separation, du, dv, nescin_filename)

      use global_variables, only: R0_plasma, nfp
      use stel_kinds
      use stel_constants
      use omp_lib

      implicit none

      character(*) :: nescin_filename
      integer :: nu, nv, nvl, geometry_option, iflag
      real(dp) :: R_specified, a, separation, du, dv, area
      real(dp), dimension(:), allocatable :: u, v, vl
      real(dp), dimension(:,:,:), allocatable :: r, drdu, drdv, normal
      real(dp), dimension(:,:), allocatable :: norm_normal
      real(dp) :: R0_to_use
      real(dp) :: angle, sinangle, cosangle, dsinangledu, dcosangledu
      real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dv, dcosangle2dv
      integer :: i, iu, iv
      real(dp) :: x_new, y_new, z_new, x_old, y_old, z_old, delta_u, delta_v, temp

      allocate(u(nu),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(v(nv),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(vl(nvl),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'

      do i = 1,nu
         u(i) = (i-1.0_dp)/nu
      end do

      do i = 1,nv
         v(i) = (i-1.0_dp)/nv
      end do

      do i = 1,nvl
         vl(i) = (i-1.0_dp)/nv
      end do

      du = u(2)-u(1)
      dv = v(2)-v(1)

      ! First dimension is the Cartesian component x, y, or z.
      allocate(r(3,nu,nvl),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(drdu(3,nu,nvl),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(drdv(3,nu,nvl),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      allocate(normal(3,nu,nvl),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'

      r = 0
      drdu = 0
      drdv = 0

      select case (geometry_option)
      case (0,1)
         ! Torus with circular cross-section

         print *,"  Building a plain circular torus."

         if (geometry_option==0) then
            R0_to_use = R0_plasma
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

               r(1,iu,iv) = (R0_to_use + a * cosangle) * cosangle2
               r(2,iu,iv) = (R0_to_use + a * cosangle) * sinangle2
               r(3,iu,iv) = a * sinangle

               drdu(1,iu,iv) = (a * dcosangledu) * cosangle2
               drdu(2,iu,iv) = (a * dcosangledu) * sinangle2
               drdu(3,iu,iv) = a * dsinangledu

               drdv(1,iu,iv) = (R0_to_use + a * cosangle) * dcosangle2dv
               drdv(2,iu,iv) = (R0_to_use + a * cosangle) * dsinangle2dv
               !drdv(3,iu,iv) = 0, so no equation needed for it here.
            end do
         end do

      case (2)

         print *,"  Constructing a surface offset from the plasma by ",separation

         ! Finite differences to use:
         ! (Numerical Recipes suggests (machine epsilon)^(1/3)
         delta_u = 1e-5;
         delta_v = 1e-5;
         ! Trick from Numerical Recipes for improved accuracy:
         temp = 1.0_dp + delta_u
         delta_u = temp - 1.0_dp
         temp = 1.0_dp + delta_v
         delta_v = temp - 1.0_dp
 

         !$OMP PARALLEL

         !$OMP MASTER
         print *,"  Number of OpenMP threads:",omp_get_num_threads()
         !$OMP END MASTER

         !$OMP DO PRIVATE(x_new,y_new,z_new,x_old,y_old,z_old)
         do iu = 1,nu
            do iv = 1,nvl

               ! Compute r:
               call compute_offset_surface_xyz_of_uv(u(iu),vl(iv),x_new,y_new,z_new,separation)
               r(1,iu,iv) = x_new
               r(2,iu,iv) = y_new
               r(3,iu,iv) = z_new

               ! Compute dr/du:
               call compute_offset_surface_xyz_of_uv(u(iu)-delta_u,vl(iv),x_old,y_old,z_old,separation)
               call compute_offset_surface_xyz_of_uv(u(iu)+delta_u,vl(iv),x_new,y_new,z_new,separation)
               drdu(1,iu,iv) = (x_new-x_old)/(2*delta_u)
               drdu(2,iu,iv) = (y_new-y_old)/(2*delta_u)
               drdu(3,iu,iv) = (z_new-z_old)/(2*delta_u)

               ! Compute dr/dv:
               call compute_offset_surface_xyz_of_uv(u(iu),vl(iv)-delta_v,x_old,y_old,z_old,separation)
               call compute_offset_surface_xyz_of_uv(u(iu),vl(iv)+delta_v,x_new,y_new,z_new,separation)
               drdv(1,iu,iv) = (x_new-x_old)/(2*delta_v)
               drdv(2,iu,iv) = (y_new-y_old)/(2*delta_v)
               drdv(3,iu,iv) = (z_new-z_old)/(2*delta_v)

            end do

         end do
         !$OMP END DO
         !$OMP END PARALLEL

      case (3)

         print *,"  Reading coil surface from nescin file ",trim(nescin_filename)

         call read_nescin(nescin_filename, r, drdu, drdv, nu, nvl, u, vl)

      case default
         print *,"Invalid setting for geometry_option: ",geometry_option
         stop
      end select

      ! Evaluate cross product:
      normal(1,:,:) = drdv(2,:,:) * drdu(3,:,:) - drdu(2,:,:) * drdv(3,:,:)
      normal(2,:,:) = drdv(3,:,:) * drdu(1,:,:) - drdu(3,:,:) * drdv(1,:,:)
      normal(3,:,:) = drdv(1,:,:) * drdu(2,:,:) - drdu(1,:,:) * drdv(2,:,:)

      allocate(norm_normal(nu, nvl),stat=iflag)
      if (iflag .ne. 0) stop 'Allocation error!'
      norm_normal = sqrt(normal(1,:,:)**2 + normal(2,:,:)**2 + normal(3,:,:)**2)

      area = du * dv * sum(norm_normal)

    end subroutine init_surface

  end module init_surface_mod
  

