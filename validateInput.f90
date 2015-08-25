subroutine validateInput

  use globalVariables

  implicit none

  if (nu_plasma < 1) then
     stop "Error! nu_plasma must be >= 1."
  end if

  if (nu_middle < 1) then
     stop "Error! nu_middle must be >= 1."
  end if

  if (nu_outer < 1) then
     stop "Error! nu_outer must be >= 1."
  end if


  if (nv_plasma < 1) then
     stop "Error! nv_plasma must be >= 1."
  end if

  if (nv_middle < 1) then
     stop "Error! nv_middle must be >= 1."
  end if

  if (nv_outer < 1) then
     stop "Error! nv_outer must be >= 1."
  end if


  if (mpol_plasma < 0) then
     stop "Error! mpol_plasma must be >= 0."
  end if

  if (mpol_middle < 0) then
     stop "Error! mpol_middle must be >= 0."
  end if

  if (mpol_outer < 0) then
     stop "Error! mpol_outer must be >= 0."
  end if


  if (ntor_plasma < 0) then
     stop "Error! ntor_plasma must be >= 0."
  end if

  if (ntor_middle < 0) then
     stop "Error! ntor_middle must be >= 0."
  end if

  if (ntor_outer < 0) then
     stop "Error! ntor_outer must be >= 0."
  end if



  if (save_level < 0) then
     stop "Error! save_level must be >= 0."
  end if

  if (save_level > 2) then
     stop "Error! save_level must be <= 2."
  end if


  if (basis_set_option < 1) then
     stop "Error! basis_set_option must be >= 1."
  end if

  if (basis_set_option > 3) then
     stop "Error! basis_set_option must be <= 3."
  end if


  if (geometry_option_plasma < 0) then
     stop "Error! geometry_option_plasma must be >= 0."
  end if

  if (geometry_option_plasma > 2) then
     stop "Error! geometry_option_plasma must be <= 2."
  end if


  if (geometry_option_middle < 0) then
     stop "Error! geometry_option_middle must be >= 0."
  end if

  if (geometry_option_middle > 3) then
     stop "Error! geometry_option_middle must be <= 3."
  end if


  if (geometry_option_outer < 0) then
     stop "Error! geometry_option_outer must be >= 0."
  end if

  if (geometry_option_outer > 3) then
     stop "Error! geometry_option_outer must be <= 3."
  end if


  if (separation_middle < 0) then
     stop "Error! separation_middle must be >= 0."
  end if

  if (separation_outer <= 0) then
     stop "Error! separation_outer must be > 0."
  end if

  if (separation_outer < separation_middle) then
     stop "Error! separation_outer must be >= separation_middle."
  end if



end subroutine validateInput
