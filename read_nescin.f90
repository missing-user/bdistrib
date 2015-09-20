subroutine read_nescin(nescin_filename, r, drdu, drdv, nu, nvl, u, vl)

  use global_variables, only: nfp, geometry_option_outer, xm, xn, mnmax, rmnc_global => rmnc, zmns_global => zmns
  use safe_open_mod
  use stel_constants
  use stel_kinds
  
  implicit none
  
  character(*) :: nescin_filename
  integer, intent(in) :: nu, nvl
  real(dp), dimension(3,nu,nvl) :: r, drdu, drdv
  real(dp), dimension(nu)  :: u
  real(dp), dimension(nvl) :: vl
  
  integer :: iunit = 7, iu, iv, iflag
  integer :: m, n, ntotal, k, mr, nr, istat
  real(dp) :: rmnc, zmns, rmns, zmnc
  real(dp) :: angle, sinangle, cosangle, dsinangledu, dsinangledv, dcosangledu, dcosangledv
  real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dv, dcosangle2dv

  character(300) :: myline
  character(*), parameter :: matchString = "------ Current Surface"

  r=0
  drdu=0
  drdv=0

  call safe_open(iunit, istat, trim(nescin_filename), 'old', 'formatted')
  if (istat .ne. 0) then
     stop 'Error opening nescin file'
  endif


  ! Skip down to the line in the file that begins with matchString
  do
     read (iunit,"(a)") myline
     if (myline(:len(matchString)) == matchString) then
        exit
     end if
  end do
  read (iunit, *)

  read (iunit, *) ntotal
  print *,"  Reading",ntotal,"Fourier modes from nescin"


  if (geometry_option_outer==4) then
     ! Clear arrays associated with the plasma surface for offsetting,
     ! and replace them with the nescin values.
     deallocate(xm,xn,rmnc_global,zmns_global)
     mnmax = ntotal
     allocate(xm(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error"
     allocate(xn(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error"
     allocate(rmnc_global(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error"
     allocate(zmns_global(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error"
  end if

  read (iunit, *)
  read (iunit, *)
  do k = 1, ntotal
     read (iunit, *) m, n, rmnc, zmns, rmns, zmnc

     if (geometry_option_outer==4) then
        ! Set arrays associated with offsetting surfaces
        xm(k) = m
        xn(k) = -n*nfp
        rmnc_global(k) = rmnc
        zmns_global(k) = zmns
     end if

     do iv = 1,nvl
        angle2 = twopi*vl(iv)/nfp
        sinangle2 = sin(angle2)
        cosangle2 = cos(angle2)
        dsinangle2dv = cosangle2*twopi/nfp
        dcosangle2dv = -sinangle2*twopi/nfp
        do iu = 1,nu
           angle = twopi*(m*u(iu) + n*vl(iv))
           sinangle = sin(angle)
           cosangle = cos(angle)
           dsinangledu = cosangle*twopi*m
           dcosangledu = -sinangle*twopi*m
           dsinangledv = cosangle*twopi*n
           dcosangledv = -sinangle*twopi*n
           
           r(1,iu,iv) = r(1,iu,iv) + rmnc * cosangle * cosangle2
           r(2,iu,iv) = r(2,iu,iv) + rmnc * cosangle * sinangle2
           r(3,iu,iv) = r(3,iu,iv) + zmns * sinangle
           
           drdu(1,iu,iv) = drdu(1,iu,iv) + rmnc * dcosangledu * cosangle2
           drdu(2,iu,iv) = drdu(2,iu,iv) + rmnc * dcosangledu * sinangle2
           drdu(3,iu,iv) = drdu(3,iu,iv) + zmns * dsinangledu
           
           drdv(1,iu,iv) = drdv(1,iu,iv) + rmnc * (dcosangledv * cosangle2 + cosangle * dcosangle2dv)
           drdv(2,iu,iv) = drdv(2,iu,iv) + rmnc * (dcosangledv * sinangle2 + cosangle * dsinangle2dv)
           drdv(3,iu,iv) = drdv(3,iu,iv) + zmns * dsinangledv
        end do
     end do
  end do
  
  
  close(iunit)

end subroutine  read_nescin
