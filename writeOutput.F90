subroutine writeOutput

  use globalVariables
  use read_wout_mod, only: nfp
  use ezcdf

  implicit none

  integer :: ierr, ncid

  ! Same convention as in VMEC:
  ! Prefix vn_ indicates the variable name used in the .nc file.
  character(len=*), parameter :: &
       vn_nfp = "nfp", &
       vn_nu_plasma = "nu_plasma", &
       vn_nv_plasma = "nv_plasma", &
       vn_nvl_plasma = "nvl_plasma", &
       vn_nu_middle = "nu_middle", &
       vn_nv_middle = "nv_middle", &
       vn_nvl_middle = "nvl_middle", &
       vn_nu_current = "nu_current", &
       vn_nv_current = "nv_current", &
       vn_nvl_current = "nvl_current", &
       vn_u_plasma = "u_plasma", &
       vn_v_plasma = "v_plasma", &
       vn_vl_plasma = "vl_plasma", &
       vn_u_middle = "u_middle", &
       vn_v_middle = "v_middle", &
       vn_vl_middle = "vl_middle", &
       vn_u_current = "u_current", &
       vn_v_current = "v_current", &
       vn_vl_current = "vl_current", &
       vn_x_plasma = "x_plasma", &
       vn_y_plasma = "y_plasma", &
       vn_z_plasma = "z_plasma", &
       vn_x_middle = "x_middle", &
       vn_y_middle = "y_middle", &
       vn_z_middle = "z_middle", &
       vn_x_current = "x_current", &
       vn_y_current = "y_current", &
       vn_z_current = "z_current"

  ! Arrays with dimension 1:
  character(len=*), parameter, dimension(1) :: &
       nu_plasma_dim = (/'nu_plasma'/), &
       nv_plasma_dim = (/'nv_plasma'/), &
       nvl_plasma_dim = (/'nvl_plasma'/), &
       nu_middle_dim = (/'nu_middle'/), &
       nv_middle_dim = (/'nv_middle'/), &
       nvl_middle_dim = (/'nvl_middle'/), &
       nu_current_dim = (/'nu_current'/), &
       nv_current_dim = (/'nv_current'/), &
       nvl_current_dim = (/'nvl_current'/)

  ! Arrays with dimension 2:
  character(len=*), parameter, dimension(2) :: &
       uvl_plasma_dim = (/'nu_plasma','nvl_plasma'/), &
       uvl_middle_dim = (/'nu_middle','nvl_middle'/), &
       uvl_current_dim = (/'nu_current','nvl_current'/)
!       uvl_plasma_dim = (/'nvl_plasma','nu_plasma'/)

  call cdf_open(ncid,outputFilename,'w',ierr)
  IF (ierr .ne. 0) then
     print *,"Error opening output file ",outputFilename
     stop
  end IF

  ! Scalars

  call cdf_define(ncid, vn_nfp, nfp)
  call cdf_define(ncid, vn_nu_plasma, nu_plasma)
  call cdf_define(ncid, vn_nv_plasma, nv_plasma)
  call cdf_define(ncid, vn_nvl_plasma, nvl_plasma)
  call cdf_define(ncid, vn_nu_middle, nu_middle)
  call cdf_define(ncid, vn_nv_middle, nv_middle)
  call cdf_define(ncid, vn_nvl_middle, nvl_middle)
  call cdf_define(ncid, vn_nu_current, nu_current)
  call cdf_define(ncid, vn_nv_current, nv_current)
  call cdf_define(ncid, vn_nvl_current, nvl_current)

  ! Arrays with dimension 1

  call cdf_define(ncid, vn_u_plasma, u_plasma, dimname=nu_plasma_dim)
  call cdf_define(ncid, vn_v_plasma, v_plasma, dimname=nv_plasma_dim)
  call cdf_define(ncid, vn_vl_plasma, vl_plasma, dimname=nvl_plasma_dim)
  call cdf_define(ncid, vn_u_middle, u_middle, dimname=nu_middle_dim)
  call cdf_define(ncid, vn_v_middle, v_middle, dimname=nv_middle_dim)
  call cdf_define(ncid, vn_vl_middle, vl_middle, dimname=nvl_middle_dim)
  call cdf_define(ncid, vn_u_current, u_current, dimname=nu_current_dim)
  call cdf_define(ncid, vn_v_current, v_current, dimname=nv_current_dim)
  call cdf_define(ncid, vn_vl_current, vl_current, dimname=nvl_current_dim)

  ! Arrays with dimension 2

  call cdf_define(ncid, vn_x_plasma, x_plasma, dimname=uvl_plasma_dim)
  call cdf_define(ncid, vn_y_plasma, y_plasma, dimname=uvl_plasma_dim)
  call cdf_define(ncid, vn_z_plasma, z_plasma, dimname=uvl_plasma_dim)
  call cdf_define(ncid, vn_x_middle, x_middle, dimname=uvl_middle_dim)
  call cdf_define(ncid, vn_y_middle, y_middle, dimname=uvl_middle_dim)
  call cdf_define(ncid, vn_z_middle, z_middle, dimname=uvl_middle_dim)
  call cdf_define(ncid, vn_x_current, x_current, dimname=uvl_current_dim)
  call cdf_define(ncid, vn_y_current, y_current, dimname=uvl_current_dim)
  call cdf_define(ncid, vn_z_current, z_current, dimname=uvl_current_dim)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! Done with cdf_define calls. Now write the data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  ! Scalars

  call cdf_write(ncid, vn_nfp, nfp)
  call cdf_write(ncid, vn_nu_plasma, nu_plasma)
  call cdf_write(ncid, vn_nv_plasma, nv_plasma)
  call cdf_write(ncid, vn_nvl_plasma, nvl_plasma)
  call cdf_write(ncid, vn_nu_middle, nu_middle)
  call cdf_write(ncid, vn_nv_middle, nv_middle)
  call cdf_write(ncid, vn_nvl_middle, nvl_middle)
  call cdf_write(ncid, vn_nu_current, nu_current)
  call cdf_write(ncid, vn_nv_current, nv_current)
  call cdf_write(ncid, vn_nvl_current, nvl_current)

  ! Arrays with dimension 1

  call cdf_write(ncid, vn_u_plasma, u_plasma)
  call cdf_write(ncid, vn_v_plasma, v_plasma)
  call cdf_write(ncid, vn_vl_plasma, vl_plasma)
  call cdf_write(ncid, vn_u_middle, u_middle)
  call cdf_write(ncid, vn_v_middle, v_middle)
  call cdf_write(ncid, vn_vl_middle, vl_middle)
  call cdf_write(ncid, vn_u_current, u_current)
  call cdf_write(ncid, vn_v_current, v_current)
  call cdf_write(ncid, vn_vl_current, vl_current)

  ! Arrays with dimension 2

  call cdf_write(ncid, vn_x_plasma, x_plasma)
  call cdf_write(ncid, vn_y_plasma, y_plasma)
  call cdf_write(ncid, vn_z_plasma, z_plasma)
  call cdf_write(ncid, vn_x_middle, x_middle)
  call cdf_write(ncid, vn_y_middle, y_middle)
  call cdf_write(ncid, vn_z_middle, z_middle)
  call cdf_write(ncid, vn_x_current, x_current)
  call cdf_write(ncid, vn_y_current, y_current)
  call cdf_write(ncid, vn_z_current, z_current)

  call cdf_close(ncid)

end subroutine writeOutput
