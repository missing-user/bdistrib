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
       vn_r_plasma  = "r_plasma", &
       vn_r_middle  = "r_middle", &
       vn_r_current = "r_current", &
       vn_drdu_plasma  = "drdu_plasma", &
       vn_drdu_middle  = "drdu_middle", &
       vn_drdu_current = "drdu_current", &
       vn_drdv_plasma  = "drdv_plasma", &
       vn_drdv_middle  = "drdv_middle", &
       vn_drdv_current = "drdv_current", &
       vn_normal_plasma = "normal_plasma", &
       vn_normal_middle = "normal_middle", &
       vn_normal_current = "normal_current", &
       vn_norm_normal_plasma  = "norm_normal_plasma", &
       vn_norm_normal_middle  = "norm_normal_middle", &
       vn_norm_normal_current = "norm_normal_current", &
       vn_currentPotential_mpol = "currentPotential_mpol", &
       vn_currentPotential_ntor = "currentPotential_ntor", &
       vn_currentPotential_mnmax = "currentPotential_mnmax", &
       vn_currentPotential_xm = "currentPotential_xm", &
       vn_currentPotential_xn = "currentPotential_xn"

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
       nvl_current_dim = (/'nvl_current'/), &
       currentPotential_mnmax_dim = (/'currentPotential_mnmax'/)

  ! Arrays with dimension 2:
  character(len=*), parameter, dimension(2) :: &
       u_vl_plasma_dim = (/'nu_plasma','nvl_plasma'/), &
       u_vl_middle_dim = (/'nu_middle','nvl_middle'/), &
       u_vl_current_dim = (/'nu_current','nvl_current'/)
!       uvl_plasma_dim = (/'nvl_plasma','nu_plasma'/)

  ! Arrays with dimension 3:
  character(len=*), parameter, dimension(3) :: &
       u_vl_xyz_plasma_dim  = (/'nu_plasma' ,'nvl_plasma' ,'xyz'/), &
       u_vl_xyz_middle_dim  = (/'nu_middle' ,'nvl_middle' ,'xyz'/), &
       u_vl_xyz_current_dim = (/'nu_current','nvl_current','xyz'/)

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
  call cdf_define(ncid, vn_currentPotential_mpol, currentPotential_mpol)
  call cdf_define(ncid, vn_currentPotential_ntor, currentPotential_ntor)
  call cdf_define(ncid, vn_currentPotential_mnmax, currentPotential_mnmax)

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
  call cdf_define(ncid, vn_currentPotential_xm, currentPotential_xm, dimname=currentPotential_mnmax_dim)
  call cdf_define(ncid, vn_currentPotential_xn, currentPotential_xn, dimname=currentPotential_mnmax_dim)

  ! Arrays with dimension 2

  call cdf_define(ncid, vn_norm_normal_plasma,  norm_normal_plasma,  dimname=u_vl_plasma_dim)
  call cdf_define(ncid, vn_norm_normal_middle,  norm_normal_middle,  dimname=u_vl_middle_dim)
  call cdf_define(ncid, vn_norm_normal_current,  norm_normal_current,  dimname=u_vl_current_dim)

  ! Arrays with dimension 3

  call cdf_define(ncid, vn_r_plasma,  r_plasma,  dimname=u_vl_xyz_plasma_dim)
  call cdf_define(ncid, vn_r_middle,  r_middle,  dimname=u_vl_xyz_middle_dim)
  call cdf_define(ncid, vn_r_current, r_current, dimname=u_vl_xyz_current_dim)

  call cdf_define(ncid, vn_drdu_plasma,  drdu_plasma,  dimname=u_vl_xyz_plasma_dim)
  call cdf_define(ncid, vn_drdu_middle,  drdu_middle,  dimname=u_vl_xyz_middle_dim)
  call cdf_define(ncid, vn_drdu_current, drdu_current, dimname=u_vl_xyz_current_dim)

  call cdf_define(ncid, vn_drdv_plasma,  drdv_plasma,  dimname=u_vl_xyz_plasma_dim)
  call cdf_define(ncid, vn_drdv_middle,  drdv_middle,  dimname=u_vl_xyz_middle_dim)
  call cdf_define(ncid, vn_drdv_current, drdv_current, dimname=u_vl_xyz_current_dim)

  call cdf_define(ncid, vn_normal_plasma,  normal_plasma,  dimname=u_vl_xyz_plasma_dim)
  call cdf_define(ncid, vn_normal_middle,  normal_middle,  dimname=u_vl_xyz_middle_dim)
  call cdf_define(ncid, vn_normal_current, normal_current, dimname=u_vl_xyz_current_dim)

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
  call cdf_write(ncid, vn_currentPotential_mpol, currentPotential_mpol)
  call cdf_write(ncid, vn_currentPotential_ntor, currentPotential_ntor)
  call cdf_write(ncid, vn_currentPotential_mnmax, currentPotential_mnmax)

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
  call cdf_write(ncid, vn_currentPotential_xm, currentPotential_xm)
  call cdf_write(ncid, vn_currentPotential_xn, currentPotential_xn)

  ! Arrays with dimension 2

  call cdf_write(ncid, vn_norm_normal_plasma,  norm_normal_plasma)
  call cdf_write(ncid, vn_norm_normal_middle,  norm_normal_middle)
  call cdf_write(ncid, vn_norm_normal_current,  norm_normal_current)

  ! Arrays with dimension 3

  call cdf_write(ncid, vn_r_plasma,  r_plasma)
  call cdf_write(ncid, vn_r_middle,  r_middle)
  call cdf_write(ncid, vn_r_current, r_current)

  call cdf_write(ncid, vn_drdu_plasma,  drdu_plasma)
  call cdf_write(ncid, vn_drdu_middle,  drdu_middle)
  call cdf_write(ncid, vn_drdu_current, drdu_current)

  call cdf_write(ncid, vn_drdv_plasma,  drdv_plasma)
  call cdf_write(ncid, vn_drdv_middle,  drdv_middle)
  call cdf_write(ncid, vn_drdv_current, drdv_current)

  call cdf_write(ncid, vn_normal_plasma,  normal_plasma)
  call cdf_write(ncid, vn_normal_middle,  normal_middle)
  call cdf_write(ncid, vn_normal_current, normal_current)

  call cdf_close(ncid)

end subroutine writeOutput
