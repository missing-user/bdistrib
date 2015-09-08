module buildInductanceMatrixMod

  implicit none

contains

  subroutine buildInductanceMatrix(inductance, r, normal, norm_normal, &
       nu, nv, mnmax, num_basis_functions, xm, xn, u, v, &
       basis_to_Fourier, area)

    use globalVariables, only: r_outer, normal_outer, u_outer, v_outer, nu_outer, nv_outer, &
         du_outer, dv_outer, mnmax_outer, num_basis_functions_outer, xn_outer, xm_outer, &
         basis_set_option, weight_option, nfp
    use stel_constants
    use stel_kinds
    use omp_lib

    implicit none

    real(dp), dimension(:,:), allocatable :: inductance
    real(dp), dimension(:,:,:), allocatable, intent(in) :: r, normal
    real(dp), dimension(:,:), allocatable, intent(in) :: norm_normal
    integer, intent(in) :: nu, nv, mnmax, num_basis_functions
    integer, dimension(:), allocatable, intent(in) :: xm, xn
    real(dp), dimension(:), allocatable, intent(in) :: u, v
    real(dp), dimension(:,:), allocatable :: basis_to_Fourier
    real(dp) :: area

    integer :: l_outer, iu, iv, iu_outer, iv_outer, ivl_outer, iflag
    real(dp) :: x, y, z, dx, dy, dz, dr2, dr32, du, dv
    integer :: imn, imn_outer, index, index_outer
    real(dp), dimension(:,:), allocatable :: inductance_xbasis, xToFourier, xToFourier_outer
    integer :: tic, toc, countrate, whichSymmetry, minSymmetry, maxSymmetry, offset
    character :: UPLO
    integer, dimension(:), allocatable :: IPIV

    ! Variables needed by BLAS DGEMM:
    character :: TRANSA='N', TRANSB='N'
    integer :: M, N, K, LDA, LDB, LDC, INFO
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

    ! Normalization so \int d^2a w f_i f_j = \delta_{i,j}
    xToFourier       = sqrt2 * xToFourier
    xToFourier_outer = sqrt2 * xToFourier_outer

    call system_clock(toc)
    print *,"  xToFourier matrices:",real(toc-tic)/countrate,"sec."

    ! If needed, compute the transformation between the Fourier functions and the basis functions.
    if (weight_option>1) then
       call system_clock(tic)
       allocate(basis_to_Fourier(num_basis_functions, num_basis_functions),stat=iflag)
       if (iflag .ne. 0) stop 'Allocation error!'
       allocate(tempMatrix(num_basis_functions, nu*nv),stat=iflag)
       if (iflag .ne. 0) stop 'Allocation error!'

       index = 0
       do iv = 1,nv
          do iu = 1,nu
             !index = (iv-1)*nu + iu
             index = index + 1
             ! The 2nd dimension of norm_normal is nvl rather than nv, but we can just ignore all periods after the 1st.
             tempMatrix(:,index) = xToFourier(:,index) * norm_normal(iu,iv)
          end do
       end do
       ! The 'basis_to_Fourier' array is not yet what the name suggests, since we will shortly do an in-place Cholesky decomposition
       basis_to_Fourier = matmul(tempMatrix,transpose(xToFourier))*du*dv*nfp/area

       call system_clock(toc)
       print *,"  Assemble A:",real(toc-tic)/countrate,"sec."
       call system_clock(tic)

       ! Compute Cholesky factorization:
       UPLO = 'L'
       call DPOTRF(UPLO, num_basis_functions, basis_to_Fourier, num_basis_functions, INFO)
       if (INFO < 0) then
          print *,"Error in Cholesky decomposition DPOTRF. The i-th argument had an illegal value. INFO=",INFO
       elseif (INFO > 0) then
          print *,"Error in Cholesky decomposition DPOTRF. The leading minor of order i is not positive definite, and the factorization could not be completed. INFO=",INFO
       end if

       call system_clock(toc)
       print *,"  Cholesky decomp:",real(toc-tic)/countrate,"sec."
       deallocate(tempMatrix)

       ! LAPACK's DPOTRF leaves the upper-triangular part nonzero, so clean it up now.
       do imn = 2,num_basis_functions
          basis_to_Fourier(1:imn-1, imn) = 0
       end do
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Done computing everything related to the basis functions.
    ! Now compute S from Boozer's eq (39)-(40).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    inductance_xbasis = 0
    call system_clock(tic)

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
                x = r(1,iu,iv)
                y = r(2,iu,iv)
                z = r(3,iu,iv)
                do l_outer = 0, (nfp-1)
                   ivl_outer = iv_outer + l_outer*nv_outer
                   dx = x - r_outer(1,iu_outer,ivl_outer)
                   dy = y - r_outer(2,iu_outer,ivl_outer)
                   dz = z - r_outer(3,iu_outer,ivl_outer)

                   dr2 = dx*dx + dy*dy + dz*dz
                   dr32 = dr2*sqrt(dr2)

                   inductance_xbasis(index,index_outer) = inductance_xbasis(index,index_outer) + &
                        (normal(1,iu,iv)*normal_outer(1,iu_outer,ivl_outer) &
                        +normal(2,iu,iv)*normal_outer(2,iu_outer,ivl_outer) &
                        +normal(3,iu,iv)*normal_outer(3,iu_outer,ivl_outer) &
                        - (3/dr2) * &
                        (normal(1,iu,iv)*dx + normal(2,iu,iv)*dy + normal(3,iu,iv)*dz) * &
                        (normal_outer(1,iu_outer,ivl_outer)*dx &
                        +normal_outer(2,iu_outer,ivl_outer)*dy &
                        +normal_outer(3,iu_outer,ivl_outer)*dz)) / dr32
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

    ! For some reason, the BLAS matrix-matrix multiplication function DGEMM sometimes causes the
    ! program to crash on Edison. In this case, you can use the following method which is slower but more reliable:
!    inductance = matmul(xToFourier, matmul(inductance_xbasis, xToFourier_outer))

    !*******************************************************
    ! Call BLAS3 subroutine DGEMM for matrix multiplications:
    !*******************************************************

    ! Here we carry out tempMatrix = inductance_xbasis * xToFourier_outer
    ! A = inductance_xbasis
    ! B = xToFourier_outer
    ! C = tempMatrix
    M = nu*nv ! # rows of A
    N = num_basis_functions_outer ! # cols of B
    K = nu_outer*nv_outer ! Common dimension of A and B
    LDA = M
    LDB = K
    LDC = M
    allocate(tempMatrix(M,N),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    tempMatrix = 0

    call DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,inductance_xbasis,LDA,xToFourier_outer,LDB,BETA,tempMatrix,LDC)

    ! Here we carry out inductance = xToFourier * tempMatrix
    ! A = inductance_xbasis
    ! B = xToFourier_outer
    ! C = tempMatrix
    M = num_basis_functions ! # rows of A
    N = num_basis_functions_outer ! # cols of B
    K = nu*nv ! Common dimension of A and B
    LDA = M
    LDB = K
    LDC = M
    call DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,xToFourier,LDA,tempMatrix,LDB,BETA,inductance,LDC)

    deallocate(tempMatrix)
    
    call system_clock(toc)
    print *,"  matmul:",real(toc-tic)/countrate,"sec."

    ! Multiply by some overall constants:
    inductance = inductance * (nfp * du * dv * du_outer * dv_outer * mu0 / (4*pi))

    ! At this point, the inductance matrix assumes we have used the Fourier basis.
    ! If needed, convert to the 'real' basis:
    if (weight_option > 1) then
       call system_clock(tic)
       allocate(IPIV(num_basis_functions),stat=iflag)
       if (iflag .ne. 0) stop 'Allocation error!'
       allocate(tempMatrix(num_basis_functions, num_basis_functions),stat=iflag)
       if (iflag .ne. 0) stop 'Allocation error!'

       tempMatrix = basis_to_Fourier
       ! DGESV will overwrite tempMatrix with the L & U factors,
       ! and overwrite the old 'inductance' with the solution.
       call DGESV(num_basis_functions, num_basis_functions_outer, tempMatrix, num_basis_functions, &
            IPIV, inductance, num_basis_functions, INFO)

       if (INFO < 0) then
          print *,"Error in DGESV. The i-th argument had an illegal value. INFO=",INFO
       elseif (INFO > 0) then
          print *,"Error in DGESV. Matrix is singular. INFO=",INFO
       end if
       deallocate(IPIV,tempMatrix)
       call system_clock(toc)
       print *,"  Convert basis:",real(toc-tic)/countrate,"sec."
    end if

    deallocate(inductance_xbasis, xToFourier, xToFourier_outer)

  end subroutine buildInductanceMatrix

end module buildInductanceMatrixMod

! Documentation for LAPACK subroutine for Cholesky decomposition:

!!$*       SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
!!$* 
!!$*       .. Scalar Arguments ..
!!$*       CHARACTER          UPLO
!!$*       INTEGER            INFO, LDA, N
!!$*       ..
!!$*       .. Array Arguments ..
!!$*       DOUBLE PRECISION   A( LDA, * )
!!$*       ..
!!$*  
!!$*
!!$*> \par Purpose:
!!$*  =============
!!$*>
!!$*> \verbatim
!!$*>
!!$*> DPOTRF computes the Cholesky factorization of a real symmetric
!!$*> positive definite matrix A.
!!$*>
!!$*> The factorization has the form
!!$*>    A = U**T * U,  if UPLO = 'U', or
!!$*>    A = L  * L**T,  if UPLO = 'L',
!!$*> where U is an upper triangular matrix and L is lower triangular.
!!$*>
!!$*> This is the block version of the algorithm, calling Level 3 BLAS.
!!$*> \endverbatim
!!$*
!!$*  Arguments:
!!$*  ==========
!!$*
!!$*> \param[in] UPLO
!!$*> \verbatim
!!$*>          UPLO is CHARACTER*1
!!$*>          = 'U':  Upper triangle of A is stored;
!!$*>          = 'L':  Lower triangle of A is stored.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] N
!!$*> \verbatim
!!$*>          N is INTEGER
!!$*>          The order of the matrix A.  N >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in,out] A
!!$*> \verbatim
!!$*>          A is DOUBLE PRECISION array, dimension (LDA,N)
!!$*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!!$*>          N-by-N upper triangular part of A contains the upper
!!$*>          triangular part of the matrix A, and the strictly lower
!!$*>          triangular part of A is not referenced.  If UPLO = 'L', the
!!$*>          leading N-by-N lower triangular part of A contains the lower
!!$*>          triangular part of the matrix A, and the strictly upper
!!$*>          triangular part of A is not referenced.
!!$*>
!!$*>          On exit, if INFO = 0, the factor U or L from the Cholesky
!!$*>          factorization A = U**T*U or A = L*L**T.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDA
!!$*> \verbatim
!!$*>          LDA is INTEGER
!!$*>          The leading dimension of the array A.  LDA >= max(1,N).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] INFO
!!$*> \verbatim
!!$*>          INFO is INTEGER
!!$*>          = 0:  successful exit
!!$*>          < 0:  if INFO = -i, the i-th argument had an illegal value
!!$*>          > 0:  if INFO = i, the leading minor of order i is not
!!$*>                positive definite, and the factorization could not be
!!$*>                completed.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Documentation for LAPACK's subroutine for solving a linear system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!$
!!$
!!$*       SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!!$* 
!!$*       .. Scalar Arguments ..
!!$*       INTEGER            INFO, LDA, LDB, N, NRHS
!!$*       ..
!!$*       .. Array Arguments ..
!!$*       INTEGER            IPIV( * )
!!$*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!!$*       ..
!!$*  
!!$*
!!$*> \par Purpose:
!!$*  =============
!!$*>
!!$*> \verbatim
!!$*>
!!$*> DGESV computes the solution to a real system of linear equations
!!$*>    A * X = B,
!!$*> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
!!$*>
!!$*> The LU decomposition with partial pivoting and row interchanges is
!!$*> used to factor A as
!!$*>    A = P * L * U,
!!$*> where P is a permutation matrix, L is unit lower triangular, and U is
!!$*> upper triangular.  The factored form of A is then used to solve the
!!$*> system of equations A * X = B.
!!$*> \endverbatim
!!$*
!!$*  Arguments:
!!$*  ==========
!!$*
!!$*> \param[in] N
!!$*> \verbatim
!!$*>          N is INTEGER
!!$*>          The number of linear equations, i.e., the order of the
!!$*>          matrix A.  N >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] NRHS
!!$*> \verbatim
!!$*>          NRHS is INTEGER
!!$*>          The number of right hand sides, i.e., the number of columns
!!$*>          of the matrix B.  NRHS >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in,out] A
!!$*> \verbatim
!!$*>          A is DOUBLE PRECISION array, dimension (LDA,N)
!!$*>          On entry, the N-by-N coefficient matrix A.
!!$*>          On exit, the factors L and U from the factorization
!!$*>          A = P*L*U; the unit diagonal elements of L are not stored.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDA
!!$*> \verbatim
!!$*>          LDA is INTEGER
!!$*>          The leading dimension of the array A.  LDA >= max(1,N).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] IPIV
!!$*> \verbatim
!!$*>          IPIV is INTEGER array, dimension (N)
!!$*>          The pivot indices that define the permutation matrix P;
!!$*>          row i of the matrix was interchanged with row IPIV(i).
!!$*> \endverbatim
!!$*>
!!$*> \param[in,out] B
!!$*> \verbatim
!!$*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!!$*>          On entry, the N-by-NRHS matrix of right hand side matrix B.
!!$*>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDB
!!$*> \verbatim
!!$*>          LDB is INTEGER
!!$*>          The leading dimension of the array B.  LDB >= max(1,N).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] INFO
!!$*> \verbatim
!!$*>          INFO is INTEGER
!!$*>          = 0:  successful exit
!!$*>          < 0:  if INFO = -i, the i-th argument had an illegal value
!!$*>          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!!$*>                has been completed, but the factor U is exactly
!!$*>                singular, so the solution could not be computed.
