subroutine initPlasma

  use globalVariables
  use read_wout_mod, only: nfp, xm, xn, rmnc, zmns, mnmax, ns
  use stel_constants

  implicit none

  integer :: i, iu, iv, imn
  real(rprec) :: angle, sinangle, cosangle, dsinangledu, dsinangledv, dcosangledu, dcosangledv
  real(rprec) :: angle2, sinangle2, cosangle2, dsinangle2dv, dcosangle2dv

  allocate(u_plasma(nu_plasma))
  allocate(v_plasma(nv_plasma))
  allocate(vl_plasma(nvl_plasma))

  do i=1,nu_plasma
     u_plasma(i) = (i-1.0_rprec)/nu_plasma
  end do

  do i=1,nv_plasma
     v_plasma(i) = (i-1.0_rprec)/nv_plasma
  end do

  do i=1,nvl_plasma
     vl_plasma(i) = (i-1.0_rprec)/nv_plasma
  end do

  ! Last coordinate is the Cartesian component x, y, or z
  allocate(r_plasma(nu_plasma,nvl_plasma,3))
  allocate(drdu_plasma(nu_plasma,nvl_plasma,3))
  allocate(drdv_plasma(nu_plasma,nvl_plasma,3))
  allocate(normal_plasma(nu_plasma,nvl_plasma,3))

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
           angle = twopi*(xm(imn)*u_plasma(iu) + xn(imn)*vl_plasma(iv)/nfp)
           cosangle = cos(angle)
           sinangle = sin(angle)
           dsinangledu = cosangle*twopi*xm(imn)
           dcosangledu = -sinangle*twopi*xm(imn)
           dsinangledv = cosangle*twopi*xn(imn)/nfp
           dcosangledv = -sinangle*twopi*xn(imn)/nfp

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

  allocate(norm_normal_plasma(nu_plasma, nvl_plasma))
  norm_normal_plasma = sqrt(normal_plasma(:,:,1)**2 + normal_plasma(:,:,2)**2 + normal_plasma(:,:,3)**2)


end subroutine initPlasma
