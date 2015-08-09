! Main program

program bdistrib

  implicit none

  print *,"This is BDISTRIB."

  call readInput()

  call initPlasma()
  call initOuterSurfaces()

  call writeOutput()

  print *,"BDISTRIB complete."

end program bdistrib
