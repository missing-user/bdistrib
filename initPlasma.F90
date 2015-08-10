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

  allocate(x_plasma(nu_plasma,nvl_plasma))
  allocate(y_plasma(nu_plasma,nvl_plasma))
  allocate(z_plasma(nu_plasma,nvl_plasma))
  allocate(dxdu_plasma(nu_plasma,nvl_plasma))
  allocate(dydu_plasma(nu_plasma,nvl_plasma))
  allocate(dzdu_plasma(nu_plasma,nvl_plasma))
  allocate(dxdv_plasma(nu_plasma,nvl_plasma))
  allocate(dydv_plasma(nu_plasma,nvl_plasma))
  allocate(dzdv_plasma(nu_plasma,nvl_plasma))
  allocate(normal_x_plasma(nu_plasma,nvl_plasma))
  allocate(normal_y_plasma(nu_plasma,nvl_plasma))
  allocate(normal_z_plasma(nu_plasma,nvl_plasma))

  x_plasma=0
  y_plasma=0
  z_plasma=0
  dxdu_plasma=0
  dydu_plasma=0
  dzdu_plasma=0
  dxdv_plasma=0
  dydv_plasma=0
  dzdv_plasma=0

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

           x_plasma(iu,iv) = x_plasma(iu,iv) + rmnc(imn,ns) * cosangle * cosangle2
           y_plasma(iu,iv) = y_plasma(iu,iv) + rmnc(imn,ns) * cosangle * sinangle2
           z_plasma(iu,iv) = z_plasma(iu,iv) + zmns(imn,ns) * sinangle

           dxdu_plasma(iu,iv) = dxdu_plasma(iu,iv) + rmnc(imn,ns) * dcosangledu * cosangle2
           dydu_plasma(iu,iv) = dydu_plasma(iu,iv) + rmnc(imn,ns) * dcosangledu * sinangle2
           dzdu_plasma(iu,iv) = dzdu_plasma(iu,iv) + zmns(imn,ns) * dsinangledu

           dxdv_plasma(iu,iv) = dxdv_plasma(iu,iv) + rmnc(imn,ns) * (dcosangledv * cosangle2 + cosangle * dcosangle2dv)
           dydv_plasma(iu,iv) = dydv_plasma(iu,iv) + rmnc(imn,ns) * (dcosangledv * sinangle2 + cosangle * dsinangle2dv)
           dzdv_plasma(iu,iv) = dzdv_plasma(iu,iv) + zmns(imn,ns) * dsinangledv
        end do
     end do
  end do

  ! Evaluate cross product
  normal_x_plasma = dydv_plasma * dzdu_plasma - dydu_plasma * dzdv_plasma
  normal_y_plasma = dzdv_plasma * dxdu_plasma - dzdu_plasma * dxdv_plasma
  normal_z_plasma = dxdv_plasma * dydu_plasma - dxdu_plasma * dydv_plasma


end subroutine initPlasma
