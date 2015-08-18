subroutine readInput

  use globalVariables

  implicit none

  integer :: numargs
  character(len=200) :: inputFilename
  integer :: fileUnit, didFileAccessWork, i
  integer, parameter :: uninitialized = -9999

  namelist / bdistrib / nu_plasma, nv_plasma, nu_middle, nv_middle, nu_outer, nv_outer, &
       surface_option_plasma, surface_option_middle, surface_option_outer, &
       R0_plasma, R0_middle, R0_outer, a_plasma, a_middle, a_outer, &
       separation_middle, separation_outer, woutFilename, pseudoinverse_thresholds, &
       save_level, n_singular_vectors_to_save, &
       mpol_plasma, ntor_plasma, mpol_middle, ntor_middle, mpol_outer, ntor_outer

  ! getcarg is in LIBSTELL
  call getcarg(1, inputFilename, numargs)

  if (numargs<1) then
     stop "One argument is required: the input namelist file, which must be named bdistrib_in.XXXXX"
  end if
  if (numargs>1) then
     print *,"WARNING: Arguments after the first will be ignored."
  end if
  if (inputFilename(1:12) .ne. "bdistrib_in.") then
     stop "Input file must be named bdistrib_in.XXX for some extension XXX"
  end if

  outputFilename = "bdistrib_out" // trim(inputFilename(12:)) // ".nc"

  pseudoinverse_thresholds = uninitialized
  ! Default: a single threshold
  pseudoinverse_thresholds(1) = 1e-10

  fileUnit=11
  open(unit=fileUnit, file=inputFilename, action="read", status="old", iostat=didFileAccessWork)
  if (didFileAccessWork /= 0) then
     print *,"Error opening input file ", trim(inputFilename)
     stop
  else
     read(fileUnit, nml=bdistrib, iostat=didFileAccessWork)
     if (didFileAccessWork /= 0) then
        print *,"Error!  I was able to open the file ", trim(inputFilename), &
               " but not read data from the bdistrib namelist in it."
        stop
     end if
     print *,"Successfully read parameters from bdistrib namelist in ", trim(inputFilename), "."
  end if
  close(unit = fileUnit)

  do i=1,nmax_pseudoinverse_thresholds
     if (pseudoinverse_thresholds(i) == uninitialized) then
        n_pseudoinverse_thresholds = i-1
        exit
     end if
  end do

end subroutine readInput
