module buildInductanceMatrixMod

  implicit none

contains

  subroutine buildInductanceMatrix(inductance, r, normal, nu, nv)

    use globalVariables, only: r_current, normal_current, u_current, vl_current, nu_current, nv_current
    use read_wout_mod, only: nfp
    use stel_constants
    use stel_kinds

    implicit none

    real(rprec), dimension(:,:), allocatable :: inductance
    real(rprec), dimension(:,:,:), allocatable :: r, normal
    integer, intent(in) :: nu, nv

    integer :: lprime, iu, iv, iuprime, ivprime, index, index_current, ivlprime, iflag
    real(rprec) :: x, y, z, dx, dy, dz, dr2, dr32

    allocate(inductance(nu*nv, nu_current*nv_current),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'

    inductance = 0

    do iu = 1, nu
       do iv = 1, nv
          index = (iv-1)*nu + iu
          x = r(iu,iv,1)
          y = r(iu,iv,2)
          z = r(iu,iv,3)
          do iuprime = 1, nu_current
             do ivprime = 1, nv_current
                index_current = (ivprime-1)*nu_current + iuprime
                do lprime = 0, (nfp-1)
                   ivlprime = ivprime + lprime*nv_current
                   dx = x - r_current(iuprime,ivlprime,1)
                   dy = y - r_current(iuprime,ivlprime,2)
                   dz = z - r_current(iuprime,ivlprime,3)

                   dr2 = dx*dx + dy*dy + dz*dz
                   dr32 = dr2*sqrt(dr2)

                   inductance(index,index_current) = inductance(index,index_current) + &
                        (normal(iu,iv,1)*normal_current(iuprime,ivlprime,1) &
                        +normal(iu,iv,2)*normal_current(iuprime,ivlprime,2) &
                        +normal(iu,iv,3)*normal_current(iuprime,ivlprime,3) &
                        - (3/dr2) * &
                        (normal(iu,iv,1)*dx + normal(iu,iv,2)*dy + normal(iu,iv,3)*dz) * &
                        (normal_current(iuprime,ivlprime,1)*dx &
                        +normal_current(iuprime,ivlprime,2)*dy &
                        +normal_current(iuprime,ivlprime,3)*dz)) / dr32
                end do
             end do
          end do
       end do
    end do

    inductance = inductance * mu0 / (4*pi)
                   
  end subroutine buildInductanceMatrix

end module buildInductanceMatrixMod
