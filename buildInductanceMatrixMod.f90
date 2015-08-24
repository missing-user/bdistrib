module buildInductanceMatrixMod

  implicit none

contains

  subroutine buildInductanceMatrix(inductance, r, normal, nu, nv, mnmax, num_basis_functions, xm, xn, u, v)

    use globalVariables, only: r_outer, normal_outer, u_outer, v_outer, nu_outer, nv_outer, &
         du_outer, dv_outer, mnmax_outer, num_basis_functions_outer, xn_outer, xm_outer, basis_set_option
    use read_wout_mod, only: nfp
    use stel_constants
    use stel_kinds
    use omp_lib

    implicit none

    real(dp), dimension(:,:), allocatable :: inductance
    real(dp), dimension(:,:,:), allocatable, intent(in) :: r, normal
    integer, intent(in) :: nu, nv, mnmax, num_basis_functions
    integer, dimension(:), allocatable, intent(in) :: xm, xn
    real(dp), dimension(:), allocatable, intent(in) :: u, v

    integer :: l_outer, iu, iv, iu_outer, iv_outer, ivl_outer, iflag
    real(dp) :: x, y, z, dx, dy, dz, dr2, dr32, du, dv
    integer :: imn, imn_outer, index, index_outer
    real(dp), dimension(:,:), allocatable :: inductance_xbasis, xToFourier, xToFourier_outer
    integer :: tic, toc, countrate, whichSymmetry, minSymmetry, maxSymmetry, offset

    ! Variables needed by BLAS DGEMM:
    character :: TRANSA='N', TRANSB='N'
    integer :: M, N, K, LDA, LDB, LDC
    real(dp) :: ALPHA=1, BETA=0
    real(dp), dimension(:,:), allocatable :: tempMatrix


    du = u(2)-u(1)
    dv = v(2)-v(1)

    allocate(inductance(num_basis_functions, num_basis_functions_outer),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(inductance_xbasis(nu*nv, nu_outer*nv_outer),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(xToFourier(num_basis_functions, nu*nv),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(xToFourier_outer(nu_outer*nv_outer, num_basis_functions_outer),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'

    select case (basis_set_option)
    case (1)
       minSymmetry = 1
       maxSymmetry = 1
    case (2)
       minSymmetry = 2
       maxSymmetry = 2
    case (3)
       minSymmetry = 1
       maxSymmetry = 2
    end select

    call system_clock(tic,countrate)

    ! These loops to assemble xToFourier and xToFourier_outer could be made faster
    ! by using the sum-angle trig identities and pretabulating the trig functions.
    ! But these loops are not the rate-limiting step, so I'll use the more transparent direct method here.
    do whichSymmetry = minSymmetry, maxSymmetry

       if (whichSymmetry==2 .and. basis_set_option==3) then
          offset = mnmax
       else
          offset = 0
       end if

       do iu = 1, nu
          do iv = 1, nv
             index = (iv-1)*nu + iu
             do imn = 1, mnmax
                if (whichSymmetry==1) then
                   xToFourier(imn, index) = sin(twopi*(xm(imn)*u(iu)+xn(imn)*v(iv)))
                else
                   xToFourier(imn+offset, index) = cos(twopi*(xm(imn)*u(iu)+xn(imn)*v(iv)))
                end if
             end do
          end do
       end do
    end do


    do whichSymmetry = minSymmetry, maxSymmetry

       if (whichSymmetry==2 .and. basis_set_option==3) then
          offset = mnmax_outer
       else
          offset = 0
       end if

       do imn_outer = 1, mnmax_outer
          do iu_outer = 1, nu_outer
             do iv_outer = 1, nv_outer
                index_outer = (iv_outer-1)*nu_outer + iu_outer
                if (whichSymmetry==1) then
                   xToFourier_outer(index_outer, imn_outer) = sin(twopi*(xm_outer(imn_outer)*u_outer(iu_outer) &
                        + xn_outer(imn_outer)*v_outer(iv_outer)))
                else
                   xToFourier_outer(index_outer, imn_outer+offset) = cos(twopi*(xm_outer(imn_outer)*u_outer(iu_outer) &
                        + xn_outer(imn_outer)*v_outer(iv_outer)))
                end if
             end do
          end do
       end do
    end do

    call system_clock(toc)
    print *,"  xToFourier matrices:",real(toc-tic)/countrate,"sec."
    call system_clock(tic)

    inductance_xbasis = 0

    !$OMP PARALLEL

    !$OMP MASTER
    print *,"  Number of OpenMP threads:",omp_get_num_threads()
    !$OMP END MASTER

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

    ! For some reason, the BLAS matrix-matrix multiplication function DGEMM is causing the
    ! program to crash on Edison. So here I'll use a slower but reliable method:
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
