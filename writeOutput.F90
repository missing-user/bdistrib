subroutine writeOutput

  use globalVariables
  use read_wout_mod
  use ezcdf

  implicit none

  integer :: ierr, ncid

  ! Same convention as in VMEC:
  ! Prefix vn_ indicates the variable name used in the .nc file.
  character(len=*), parameter :: &
       vn_nu_plasma = "nu_plasma", &
       vn_nv_plasma = "nv_plasma", &
       vn_nu_middle = "nu_middle", &
       vn_nv_middle = "nv_middle", &
       vn_nu_current = "nu_current", &
       vn_nv_current = "nv_current", &
       vn_u_plasma = "u_plasma", &
       vn_v_plasma = "v_plasma", &
       vn_u_middle = "u_middle", &
       vn_v_middle = "v_middle", &
       vn_u_current = "u_current", &
       vn_v_current = "v_current"

  ! Arrays with dimension 1:
  character(len=*), parameter, dimension(1) :: &
       nu_plasma_dim = (/'nu_plasma'/), &
       nv_plasma_dim = (/'nv_plasma'/), &
       nu_middle_dim = (/'nu_middle'/), &
       nv_middle_dim = (/'nv_middle'/), &
       nu_current_dim = (/'nu_current'/), &
       nv_current_dim = (/'nv_current'/)

  print *,"nfp=",nfp
  print *,"Rmajor=",Rmajor

  call cdf_open(ncid,outputFilename,'w',ierr)
  IF (ierr .ne. 0) then
     print *,"Error opening output file ",outputFilename
     stop
  end IF

  ! Scalars

  call cdf_define(ncid, vn_nu_plasma, nu_plasma)
  call cdf_define(ncid, vn_nv_plasma, nv_plasma)
  call cdf_define(ncid, vn_nu_middle, nu_middle)
  call cdf_define(ncid, vn_nv_middle, nv_middle)
  call cdf_define(ncid, vn_nu_current, nu_current)
  call cdf_define(ncid, vn_nv_current, nv_current)

  ! Arrays with dimension 1

  call cdf_define(ncid, vn_u_plasma, u_plasma, dimname=nu_plasma_dim)
  call cdf_define(ncid, vn_v_plasma, v_plasma, dimname=nv_plasma_dim)
!!$  call cdf_define(ncid, vn_u_middle, u_middle, dimname=nu_middle_dim)
!!$  call cdf_define(ncid, vn_v_middle, v_middle, dimname=nv_middle_dim)
!!$  call cdf_define(ncid, vn_u_current, u_current, dimname=nu_current_dim)
!!$  call cdf_define(ncid, vn_v_current, v_current, dimname=nv_current_dim)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! Done with cdf_define calls. Now write the data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  ! Scalars

  call cdf_write(ncid, vn_nu_plasma, nu_plasma)
  call cdf_write(ncid, vn_nv_plasma, nv_plasma)
  call cdf_write(ncid, vn_nu_middle, nu_middle)
  call cdf_write(ncid, vn_nv_middle, nv_middle)
  call cdf_write(ncid, vn_nu_current, nu_current)
  call cdf_write(ncid, vn_nv_current, nv_current)

  ! Arrays with dimension 1

  call cdf_write(ncid, vn_u_plasma, u_plasma)
  call cdf_write(ncid, vn_v_plasma, v_plasma)
!!$  call cdf_write(ncid, vn_u_middle, u_middle)
!!$  call cdf_write(ncid, vn_v_middle, v_middle)
!!$  call cdf_write(ncid, vn_u_current, u_current)
!!$  call cdf_write(ncid, vn_v_current, v_current)

  call cdf_close(ncid)

end subroutine writeOutput
