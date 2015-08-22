module buildInductanceMatrixMod

  implicit none

contains

  subroutine buildInductanceMatrix(inductance, r, normal, nu, nv, mnmax, xm, xn, u, v)

    use globalVariables, only: r_outer, normal_outer, u_outer, v_outer, nu_outer, nv_outer, &
         du_outer, dv_outer, mnmax_outer, xn_outer, xm_outer
    use read_wout_mod, only: nfp
    use stel_constants
    use stel_kinds
    use omp_lib

    implicit none

    real(dp), dimension(:,:), allocatable :: inductance
    real(dp), dimension(:,:,:), allocatable, intent(in) :: r, normal
    integer, intent(in) :: nu, nv, mnmax
    integer, dimension(:), allocatable, intent(in) :: xm, xn
    real(dp), dimension(:), allocatable, intent(in) :: u, v

    integer :: l_outer, iu, iv, iu_outer, iv_outer, ivl_outer, iflag
    real(dp) :: x, y, z, dx, dy, dz, dr2, dr32, du, dv
    integer :: imn, imn_outer, index, index_outer
    real(dp), dimension(:,:), allocatable :: inductance_xbasis, xToFourier, xToFourier_outer
    integer :: tic, toc, countrate, omp_num_threads

    ! Variables needed by BLAS DGEMM:
    character :: TRANSA='N', TRANSB='N'
    integer :: M, N, K, LDA, LDB, LDC
    real(dp) :: ALPHA=1, BETA=0
    real(dp), dimension(:,:), allocatable :: tempMatrix


    du = u(2)-u(1)
    dv = v(2)-v(1)

    allocate(inductance(mnmax, mnmax_outer),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(inductance_xbasis(nu*nv, nu_outer*nv_outer),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(xToFourier(mnmax, nu*nv),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(xToFourier_outer(nu_outer*nv_outer, mnmax_outer),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'

    call system_clock(tic,countrate)

    do iu = 1, nu
       do iv = 1, nv
          index = (iv-1)*nu + iu
          do imn = 1, mnmax
             xToFourier(imn, index) = sin(twopi*(xm(imn)*u(iu)+xn(imn)*v(iv)))
          end do
       end do
    end do

    do iu_outer = 1, nu_outer
       do iv_outer = 1, nv_outer
          index_outer = (iv_outer-1)*nu_outer + iu_outer
          do imn_outer = 1, mnmax_outer
             xToFourier_outer(index_outer, imn_outer) = sin(twopi*(xm_outer(imn_outer)*u_outer(iu_outer) &
                  + xn_outer(imn_outer)*v_outer(iv_outer)))
          end do
       end do
    end do

    call system_clock(toc)
    print *,"  xToFourier matrices:",real(toc-tic)/countrate,"sec."
    call system_clock(tic)

    inductance_xbasis = 0

    !$OMP PARALLEL

    !$OMP CRITICAL
    !omp_num_threads = omp_get_thread_num()
    write (*,*) "  Hello from thread ",omp_get_thread_num()
    !write (*,*) "  Hello from thread ",omp_num_threads
    !$OMP END CRITICAL

    !$OMP DO PRIVATE(index_outer,index,x,y,z,ivl_outer,dx,dy,dz,dr2,dr32)
    do iv_outer = 1, nv_outer
       do iu_outer = 1, nu_outer
          index_outer = (iv_outer-1)*nu_outer + iu_outer
          do iv = 1, nv
             do iu = 1, nu
                index = (iv-1)*nu + iu
                x = r(iu,iv,1)
                y = r(iu,iv,2)
                z = r(iu,iv,3)
                do l_outer = 0, (nfp-1)
                   ivl_outer = iv_outer + l_outer*nv_outer
                   dx = x - r_outer(iu_outer,ivl_outer,1)
                   dy = y - r_outer(iu_outer,ivl_outer,2)
                   dz = z - r_outer(iu_outer,ivl_outer,3)

                   dr2 = dx*dx + dy*dy + dz*dz
                   dr32 = dr2*sqrt(dr2)

                   inductance_xbasis(index,index_outer) = inductance_xbasis(index,index_outer) + &
                        (normal(iu,iv,1)*normal_outer(iu_outer,ivl_outer,1) &
                        +normal(iu,iv,2)*normal_outer(iu_outer,ivl_outer,2) &
                        +normal(iu,iv,3)*normal_outer(iu_outer,ivl_outer,3) &
                        - (3/dr2) * &
                        (normal(iu,iv,1)*dx + normal(iu,iv,2)*dy + normal(iu,iv,3)*dz) * &
                        (normal_outer(iu_outer,ivl_outer,1)*dx &
                        +normal_outer(iu_outer,ivl_outer,2)*dy &
                        +normal_outer(iu_outer,ivl_outer,3)*dz)) / dr32
                end do
             end do
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    call system_clock(toc)
    print *,"  inductance_xbasis:",real(toc-tic)/countrate,"sec."

    call system_clock(tic)
    inductance = matmul(xToFourier, matmul(inductance_xbasis, xToFourier_outer))

!!$    !*******************************************************
!!$    ! Call BLAS3 subroutine DGEMM for matrix multiplications:
!!$    !*******************************************************
!!$
!!$
!!$    ! tempMatrix = inductance_xbasis * xToFourier_outer
!!$    ! A = inductance_xbasis
!!$    ! B = xToFourier_outer
!!$    ! C = tempMatrix
!!$    M = nu*nv ! # rows of A
!!$    N = mnmax_outer ! # cols of B
!!$    K = nu_outer*nv_outer ! Common dimension of A and B
!!$    LDA = M
!!$    LDB = K
!!$    LDC = M
!!$    allocate(tempMatrix(M,N))
!!$    print *,"AAA"
!!$
!!$    !$OMP PARALLEL
!!$    !$OMP CRITICAL
!!$    print *, "  Hello from thread ",omp_get_thread_num()
!!$    !$OMP END CRITICAL
!!$    !$OMP END PARALLEL
!!$
!!$    !$OMP BARRIER
!!$    call DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,inductance_xbasis,LDA,xToFourier_outer,LDB,BETA,tempMatrix,LDC)
!!$
!!$
!!$    !$OMP PARALLEL
!!$    !$OMP CRITICAL
!!$    print *, "  Hello from thread ",omp_get_thread_num()
!!$    !$OMP END CRITICAL
!!$    !$OMP END PARALLEL
!!$    
!!$    print *,"BBB"
!!$!    !$OMP END PARALLEL
!!$    ! inductance = xToFourier * tempMatrix
!!$    ! A = inductance_xbasis
!!$    ! B = xToFourier_outer
!!$    ! C = tempMatrix
!!$    M = mnmax ! # rows of A
!!$    N = mnmax_outer ! # cols of B
!!$    K = nu*nv ! Common dimension of A and B
!!$    LDA = M
!!$    LDB = K
!!$    LDC = M
!!$    !$OMP BARRIER
!!$    call DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,xToFourier,LDA,tempMatrix,LDB,BETA,inductance,LDC)
!!$    !$OMP BARRIER
!!$
!!$    !$OMP PARALLEL
!!$    !$OMP CRITICAL
!!$    print *, "  Hello from thread ",omp_get_thread_num()
!!$    !$OMP END CRITICAL
!!$    !$OMP END PARALLEL
!!$    print *,"DDD"
!!$    deallocate(tempMatrix)
    
    call system_clock(toc)
    print *,"  matmul:",real(toc-tic)/countrate,"sec."

    ! Multiply by some overall constants:
    inductance = inductance * (2 * nfp * du * dv * du_outer * dv_outer * mu0 / (4*pi))
                   
    deallocate(inductance_xbasis, xToFourier, xToFourier_outer)

  end subroutine buildInductanceMatrix

end module buildInductanceMatrixMod
