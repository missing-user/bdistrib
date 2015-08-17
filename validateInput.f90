subroutine validateInput

  use globalVariables

  implicit none

  if (save_level < 0) then
     stop "Error! save_level must be >= 0."
  end if

  if (save_level > 1) then
     stop "Error! save_level must be <= 1."
  end if

end subroutine validateInput
