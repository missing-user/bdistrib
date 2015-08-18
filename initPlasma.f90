subroutine initPlasma

  use globalVariables
  use read_wout_mod, only: nfp, xm, xn, rmnc, zmns, mnmax, ns, Rmajor, read_wout_file
  use stel_constants

  implicit none

  integer :: i, iu, iv, imn, tic, toc, countrate, iflag, ierr, iopen
  real(dp) :: angle, sinangle, cosangle, dsinangledu, dsinangledv, dcosangledu, dcosangledv
  real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dv, dcosangle2dv

  call system_clock(tic, countrate)
  print *,"Initializing plasma surface."

  select case (surface_option_plasma)
  case (0,1)
     ! Plain circular torus
     print *,"  Building a plain circular torus."

     nfp = 1
     mnmax = 2
     Rmajor = R0_plasma
     ns = 1
     allocate(xm(2),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(xn(2),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(rmnc(2,1),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     allocate(zmns(2,1),stat=iflag)
     if (iflag .ne. 0) stop 'Allocation error!'
     xm = (/0,1/)
     xn = (/0,0/)
     rmnc(1,1) = R0_plasma
     rmnc(2,1) = a_plasma
     zmns(1,1) = 0
     zmns(2,1) = a_plasma
     
  case (2)
     call read_wout_file(woutFilename, ierr, iopen)
     if (iopen .ne. 0) stop 'error opening wout in bn_read_vmecf90'
     if (ierr .ne. 0) stop 'error reading wout in bn_read_vmecf90'
     print *,"  Successfully read VMEC data from ",woutFilename

  case default
     print *,"Error! Invalid setting for surface_option_plasma:",surface_option_plasma
     stop
  end select

  nvl_plasma = nv_plasma * nfp
  nvl_middle = nv_middle * nfp
  nvl_outer = nv_outer * nfp

  allocate(u_plasma(nu_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(v_plasma(nv_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(vl_plasma(nvl_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  do i=1,nu_plasma
     u_plasma(i) = (i-1.0_dp)/nu_plasma
  end do

  do i=1,nv_plasma
     v_plasma(i) = (i-1.0_dp)/nv_plasma
  end do

  do i=1,nvl_plasma
     vl_plasma(i) = (i-1.0_dp)/nv_plasma
  end do

  ! Last coordinate is the Cartesian component x, y, or z
  allocate(r_plasma(nu_plasma,nvl_plasma,3),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(drdu_plasma(nu_plasma,nvl_plasma,3),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(drdv_plasma(nu_plasma,nvl_plasma,3),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(normal_plasma(nu_plasma,nvl_plasma,3),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  r_plasma=0
  drdu_plasma=0
  drdv_plasma=0
 
  do iv = 1,nvl_plasma
     angle2 = twopi*vl_plasma(iv)/nfp
     sinangle2 = sin(angle2)
     cosangle2 = cos(angle2)
     dsinangle2dv = cosangle2*twopi/nfp
     dcosangle2dv = -sinangle2*twopi/nfp
     do iu = 1,nu_plasma
        do imn = 1,mnmax
           angle = twopi*(xm(imn)*u_plasma(iu) - xn(imn)*vl_plasma(iv)/nfp)
           sinangle = sin(angle)
           cosangle = cos(angle)
           dsinangledu = cosangle*twopi*xm(imn)
           dcosangledu = -sinangle*twopi*xm(imn)
           dsinangledv = -cosangle*twopi*xn(imn)/nfp
           dcosangledv = sinangle*twopi*xn(imn)/nfp

           r_plasma(iu,iv,1) = r_plasma(iu,iv,1) + rmnc(imn,ns) * cosangle * cosangle2
           r_plasma(iu,iv,2) = r_plasma(iu,iv,2) + rmnc(imn,ns) * cosangle * sinangle2
           r_plasma(iu,iv,3) = r_plasma(iu,iv,3) + zmns(imn,ns) * sinangle

           drdu_plasma(iu,iv,1) = drdu_plasma(iu,iv,1) + rmnc(imn,ns) * dcosangledu * cosangle2
           drdu_plasma(iu,iv,2) = drdu_plasma(iu,iv,2) + rmnc(imn,ns) * dcosangledu * sinangle2
           drdu_plasma(iu,iv,3) = drdu_plasma(iu,iv,3) + zmns(imn,ns) * dsinangledu

           drdv_plasma(iu,iv,1) = drdv_plasma(iu,iv,1) + rmnc(imn,ns) * (dcosangledv * cosangle2 + cosangle * dcosangle2dv)
           drdv_plasma(iu,iv,2) = drdv_plasma(iu,iv,2) + rmnc(imn,ns) * (dcosangledv * sinangle2 + cosangle * dsinangle2dv)
           drdv_plasma(iu,iv,3) = drdv_plasma(iu,iv,3) + zmns(imn,ns) * dsinangledv
        end do
     end do
  end do

  ! Evaluate cross product
  normal_plasma(:,:,1) = drdv_plasma(:,:,2) * drdu_plasma(:,:,3) - drdu_plasma(:,:,2) * drdv_plasma(:,:,3)
  normal_plasma(:,:,2) = drdv_plasma(:,:,3) * drdu_plasma(:,:,1) - drdu_plasma(:,:,3) * drdv_plasma(:,:,1)
  normal_plasma(:,:,3) = drdv_plasma(:,:,1) * drdu_plasma(:,:,2) - drdu_plasma(:,:,1) * drdv_plasma(:,:,2)

  allocate(norm_normal_plasma(nu_plasma, nvl_plasma),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  norm_normal_plasma = sqrt(normal_plasma(:,:,1)**2 + normal_plasma(:,:,2)**2 + normal_plasma(:,:,3)**2)

  du_plasma = u_plasma(2)-u_plasma(1)
  dv_plasma = v_plasma(2)-v_plasma(1)
  
  call system_clock(toc)
  print *,"Done initializing plasma surface. Took ",real(toc-tic)/countrate," sec."

end subroutine initPlasma
