MODULE utilities
  ! use, intrinsic :: iso_fortran_env, only: dp=>real64

  IMPLICIT NONE

  ! TYPE cholesky_cache_item
  !   integer :: ninds = 0
  !   ! Usiamo un allocatable perché la dimensione (ninds x ninds) varia per ogni pattern
  !   double precision, dimension(:,:), allocatable :: LT_phi
  ! END TYPE cholesky_cache_item

  !---------------------------------------------------------------------------
  ! Interfaccia per DSYEV (Autovalori/Autovettori di una matrice simmetrica)
  ! Descrizione: Calcola gli autovalori e (opzionalmente) gli autovettori
  !              di una matrice reale simmetrica.
  ! Documentazione LAPACK: http://www.netlib.org/lapack/explore-html/d4/d86/dsyev_8f.html
  !---------------------------------------------------------------------------
  INTERFACE
    SUBROUTINE DSYEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
      CHARACTER(LEN=*), INTENT(IN) :: JOBZ  ! 'N' (solo autovalori) o 'V' (autovalori e autovettori)
      CHARACTER(LEN=*), INTENT(IN) :: UPLO  ! 'U' (upper triangle) o 'L' (lower triangle)
      INTEGER, INTENT(IN) :: N     ! Ordine della matrice A
      DOUBLE PRECISION, DIMENSION(LDA,*), INTENT(INOUT) :: A ! Matrice simmetrica di input/output (autovettori se JOBZ='V')
      INTEGER, INTENT(IN) :: LDA   ! Leading dimension di A (>=N)
      DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: W   ! Autovalori in ordine crescente
      DOUBLE PRECISION, DIMENSION(*), INTENT(OUT) :: WORK ! Area di lavoro
      INTEGER, INTENT(IN) :: LWORK ! Dimensione di WORK (almeno max(1, 3*N-1) o 2*N-1 per N=1)
      INTEGER, INTENT(OUT) :: INFO  ! =0 per successo, <0 per errore input, >0 per mancata convergenza
    END SUBROUTINE DSYEV
  END INTERFACE

    INTERFACE
    SUBROUTINE DPOTRS(UPLO, N, NRHS, A, LDA, B, LDB, INFO)
      CHARACTER(LEN=*), INTENT(IN) :: UPLO
      INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
      DOUBLE PRECISION, DIMENSION(LDA,*), INTENT(IN) :: A
      DOUBLE PRECISION, DIMENSION(LDB,*), INTENT(INOUT) :: B
      INTEGER, INTENT(OUT) :: INFO
    END SUBROUTINE DPOTRS
  END INTERFACE

  INTERFACE
    SUBROUTINE DCOPY(N, DX, INCX, DY, INCY)
      INTEGER, INTENT(IN) :: N, INCX, INCY
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: DX
      DOUBLE PRECISION, DIMENSION(*), INTENT(OUT) :: DY
    END SUBROUTINE DCOPY
  END INTERFACE


  !---------------------------------------------------------------------------
  ! Interfaccia per DGEMM (Moltiplicazione Generale Matrice)
  ! Descrizione: Esegue una moltiplicazione di matrici C = alpha*op(A)*op(B) + beta*C
  ! Documentazione LAPACK/BLAS: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
  !---------------------------------------------------------------------------
  INTERFACE
    SUBROUTINE DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
      CHARACTER(LEN=*), INTENT(IN) :: TRANSA ! 'N' (no transpose), 'T' (transpose), 'C' (conjugate transpose) for matrix A
      CHARACTER(LEN=*), INTENT(IN) :: TRANSB ! 'N', 'T', 'C' for matrix B
      INTEGER, INTENT(IN) :: M      ! Numero di righe della matrice op(A) e C
      INTEGER, INTENT(IN) :: N      ! Numero di colonne della matrice op(B) e C
      INTEGER, INTENT(IN) :: K      ! Numero di colonne di op(A) e righe di op(B)
      DOUBLE PRECISION, INTENT(IN) :: ALPHA    ! Scalare alpha
      DOUBLE PRECISION, DIMENSION(LDA,*), INTENT(IN) :: A ! Matrice A
      INTEGER, INTENT(IN) :: LDA    ! Leading dimension di A
      DOUBLE PRECISION, DIMENSION(LDB,*), INTENT(IN) :: B ! Matrice B
      INTEGER, INTENT(IN) :: LDB    ! Leading dimension di B
      DOUBLE PRECISION, INTENT(IN) :: BETA     ! Scalare beta
      DOUBLE PRECISION, DIMENSION(LDC,*), INTENT(INOUT) :: C ! Matrice C (output)
      INTEGER, INTENT(IN) :: LDC    ! Leading dimension di C
    END SUBROUTINE DGEMM
  END INTERFACE

  !---------------------------------------------------------------------------
  ! Interfaccia per DPOTRF (Decomposizione di Cholesky)
  ! Descrizione: Calcola la fattorizzazione di Cholesky di una matrice
  !              simmetrica definita positiva A = U^T*U o A = L*L^T.
  ! Documentazione LAPACK: http://www.netlib.org/lapack/explore-html/d8/d53/dpotrf_8f.html
  !---------------------------------------------------------------------------
  INTERFACE
    SUBROUTINE DPOTRF(UPLO, N, A, LDA, INFO)
      CHARACTER(LEN=*), INTENT(IN) :: UPLO  ! 'U' (upper triangle) o 'L' (lower triangle)
      INTEGER, INTENT(IN) :: N     ! Ordine della matrice A
      DOUBLE PRECISION, DIMENSION(LDA,*), INTENT(INOUT) :: A ! Matrice di input (contiene L o U in output)
      INTEGER, INTENT(IN) :: LDA   ! Leading dimension di A
      INTEGER, INTENT(OUT) :: INFO  ! =0 per successo, <0 per errore input, >0 se matrice non definita positiva
    END SUBROUTINE DPOTRF
  END INTERFACE

  !---------------------------------------------------------------------------
  ! Interfaccia per DPOTRI (Inversa di una matrice da Cholesky)
  ! Descrizione: Calcola l'inversa di una matrice simmetrica definita positiva
  !              usando la sua fattorizzazione di Cholesky.
  ! Documentazione LAPACK: http://www.netlib.org/lapack/explore-html/d6/d34/dpotri_8f.html
  !---------------------------------------------------------------------------
  INTERFACE
    SUBROUTINE DPOTRI(UPLO, N, A, LDA, INFO)
      CHARACTER(LEN=*), INTENT(IN) :: UPLO  ! 'U' (upper triangle) o 'L' (lower triangle)
      INTEGER, INTENT(IN) :: N     ! Ordine della matrice A
      DOUBLE PRECISION, DIMENSION(LDA,*), INTENT(INOUT) :: A ! Matrice di input (fattorizzazione di Cholesky)/output (inversa)
      INTEGER, INTENT(IN) :: LDA   ! Leading dimension di A
      INTEGER, INTENT(OUT) :: INFO  ! =0 per successo, <0 per errore input, >0 se singolare
    END SUBROUTINE DPOTRI
  END INTERFACE

  INTERFACE
    SUBROUTINE DGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
      CHARACTER(LEN=*), INTENT(IN) :: TRANS
      INTEGER, INTENT(IN) :: M, N, LDA, INCX, INCY
      DOUBLE PRECISION, INTENT(IN) :: ALPHA, BETA
      DOUBLE PRECISION, DIMENSION(LDA,*), INTENT(IN) :: A
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X
      DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: Y
    END SUBROUTINE DGEMV
  END INTERFACE




  CONTAINS

  ! dll fortran
  function alnorm ( x, upper )
    !*****************************************************************************80
    !
    !! ALNORM computes the cumulative density of the standard normal distribution.
    !
    !  Modified:
    !
    !    13 January 2008
    !
    !  Author:
    !
    !    David Hill
    !    FORTRAN90 version by John Burkardt
    !
    !  Reference:
    !
    !    David Hill,
    !    Algorithm AS 66:
    !    The Normal Integral,
    !    Applied Statistics,
    !    Volume 22, Number 3, 1973, pages 424-427.
    !
    !  Parameters:
    !
    !    Input, double precision X, is one endpoint of the semi-infinite interval
    !    over which the integration takes place.
    !
    !    Input, logical UPPER, determines whether the upper or lower
    !    interval is to be integrated:
    !    .TRUE.  => integrate from X to + Infinity;
    !    .FALSE. => integrate from - Infinity to X.
    !
    !    Output, double precision ALNORM, the integral of the standard normal
    !    distribution over the desired interval.
    !
    implicit none

    double precision, parameter :: a1 = 5.75885480458D+00
    double precision, parameter :: a2 = 2.62433121679D+00
    double precision, parameter :: a3 = 5.92885724438D+00
    double precision alnorm
    double precision, parameter :: b1 = -29.8213557807D+00
    double precision, parameter :: b2 = 48.6959930692D+00
    double precision, parameter :: c1 = -0.000000038052D+00
    double precision, parameter :: c2 = 0.000398064794D+00
    double precision, parameter :: c3 = -0.151679116635D+00
    double precision, parameter :: c4 = 4.8385912808D+00
    double precision, parameter :: c5 = 0.742380924027D+00
    double precision, parameter :: c6 = 3.99019417011D+00
    double precision, parameter :: con = 1.28D+00
    double precision, parameter :: d1 = 1.00000615302D+00
    double precision, parameter :: d2 = 1.98615381364D+00
    double precision, parameter :: d3 = 5.29330324926D+00
    double precision, parameter :: d4 = -15.1508972451D+00
    double precision, parameter :: d5 = 30.789933034D+00
    double precision, parameter :: ltone = 7.0D+00
    double precision, parameter :: p = 0.398942280444D+00
    double precision, parameter :: q = 0.39990348504D+00
    double precision, parameter :: r = 0.398942280385D+00
    logical up
    logical upper
    double precision, parameter :: utzero = 18.66D+00
    double precision x
    double precision y
    double precision z

    up = upper
    z = x

    if ( z < 0.0D+00 ) then
      up = .not. up
      z = - z
    end if

    if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then
      if ( up ) then
          alnorm = 0.0D+00
      else
          alnorm = 1.0D+00
      end if
      return
    end if

    y = 0.5D+00 * z * z

    if ( z <= con ) then
      alnorm = 0.5D+00 - z * ( p - q * y  &
            / ( y + a1 + b1/ ( y + a2 + b2 / ( y + a3 ))))
    else
      alnorm = r * exp ( - y )   &
            / ( z + c1 + d1 / ( z + c2 + d2 / ( z + c3 + d3  &
            / ( z + c4 + d4 / ( z + c5 + d5 / ( z + c6 ))))))
    end if

    if ( .not. up ) then
      alnorm = 1.0D+00 - alnorm
    end if

    return
  end function alnorm



  function alnorm_fast(x, upper) result(p)
    implicit none
    double precision, intent(in) :: x
    logical, intent(in) :: upper
    double precision :: p
    double precision, parameter :: inv_sqrt_2 = 0.70710678118654752440084436210485d0
    if (upper) then
       p = 0.5d0 * erfc( x*inv_sqrt_2 )   ! coda superiore
    else
       p = 0.5d0 * erfc(-x*inv_sqrt_2 )   !coda inferiore
    end if
  end function

  SUBROUTINE reconstruct_covariance_matrix(cov_vector, variances_vector, dimA, cov_matrix)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: dimA
    double precision, INTENT(IN), OPTIONAL    :: cov_vector(dimA*(dimA-1)/2)       ! Vettore dei coefficienti di covarianza (elementi non diagonali)
    double precision, INTENT(IN)    :: variances_vector(dimA) ! Vettore delle varianze (elementi diagonali)

    ! Ora cov_matrix è INTENT(OUT) o INTENT(INOUT) e deve essere allocata PRIMA di chiamare la subroutine
    double precision, INTENT(OUT)   :: cov_matrix(dimA,dimA)     ! Matrice di covarianza da ricostruire


    INTEGER :: k               ! Contatore per cov_vector
    INTEGER :: i, j            ! Indici per la matrice



    cov_matrix = 0.0

    ! 4. Inserisci le varianze sulla diagonale
    DO i = 1, dimA
      cov_matrix(i, i) = variances_vector(i)
    END DO

    IF (PRESENT(cov_vector)) THEN
      k = 1
      DO j = 1, dimA
        DO i = 1, j-1 ! Iterare solo sugli elementi della triangolare superiore
          cov_matrix(i, j) = cov_vector(k)
          cov_matrix(j, i) = cov_vector(k)
          k = k + 1
        END DO
      END DO
    END IF
  END SUBROUTINE reconstruct_covariance_matrix



  SUBROUTINE build_Lphi_ar1_full(n, Q, A, sqrt_sig1, L_phi)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, Q
    DOUBLE PRECISION, INTENT(IN) :: A(Q,Q)
    DOUBLE PRECISION, INTENT(IN) :: sqrt_sig1(Q)          ! sqrt(diag(Sigma1))
    DOUBLE PRECISION, INTENT(OUT) :: L_phi(n*Q, n*Q)

    INTEGER :: j, m, r, c, kpow
    INTEGER :: row0, col0
    DOUBLE PRECISION, ALLOCATABLE :: Apow(:,:,:)

    ALLOCATE(Apow(Q,Q,n))

    ! Apow(:,:,1)=A^0=I
    Apow(:,:,1) = 0.d0
    DO r = 1, Q
      Apow(r,r,1) = 1.d0
    END DO

    ! Apow(:,:,k)=A^(k-1)
    DO kpow = 2, n
      Apow(:,:,kpow) = MATMUL(Apow(:,:,kpow-1), A)
    END DO

    L_phi = 0.d0

    ! Blocchi: per j=1..n, m=1..j
    DO j = 1, n
      row0 = (j-1)*Q

      DO m = 1, j
        col0 = (m-1)*Q

        IF (m == 1) THEN
          ! Block(j,1) = A^(j-1) * sqrt(Sigma1)  (Sigma1 diag -> scala le colonne)
          DO c = 1, Q
            DO r = 1, Q
              L_phi(row0 + r, col0 + c) = Apow(r,c,j) * sqrt_sig1(c)
            END DO
          END DO

        ELSE
          ! Block(j,m) = A^(j-m)
          kpow = j - m + 1   ! perché Apow(:,:,1)=A^0
          DO c = 1, Q
            DO r = 1, Q
              L_phi(row0 + r, col0 + c) = Apow(r,c,kpow)
            END DO
          END DO
        END IF

      END DO
    END DO

    DEALLOCATE(Apow)
  END SUBROUTINE build_Lphi_ar1_full


  SUBROUTINE eta_from_z_ar1(n, Q, kappa, A, sqrt_sig1, z, eta)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, Q
    DOUBLE PRECISION, INTENT(IN) :: kappa(n*Q)
    DOUBLE PRECISION, INTENT(IN) :: A(Q,Q), sqrt_sig1(Q)
    DOUBLE PRECISION, INTENT(IN) :: z(n*Q)
    DOUBLE PRECISION, INTENT(OUT):: eta(n*Q)

    INTEGER :: j, t, qq, idx1, idx0
    DOUBLE PRECISION :: tmp(Q), v(Q)

    ! t = 1
    idx1 = 1
    DO qq = 1, Q
      eta(idx1 + qq - 1) = kappa(idx1 + qq - 1) + sqrt_sig1(qq) * z(idx1 + qq - 1)
    END DO

    ! t = 2..n
    DO t = 2, n
      idx1 = (t-1)*Q + 1        ! start corrente
      idx0 = (t-2)*Q + 1        ! start precedente

      ! tmp = eta_{t-1} - kappa_{t-1}
      DO qq = 1, Q
          tmp(qq) = eta(idx0 + qq - 1) - kappa(idx0 + qq - 1)
      END DO

      ! v = A * tmp
      DO qq = 1, Q
          v(qq) = 0.d0
      END DO
      ! (volendo qqui puoi mettere DGEMV, ma loop va già bene)
      DO qq = 1, Q
          v(qq) = 0.d0
      END DO
      DO qq = 1, Q
          ! v(qq) = sum_j A(qq,j)*tmp(j)
          DO j = 1, Q
            v(qq) = v(qq) + A(qq,j) * tmp(j)
          END DO
      END DO

      ! eta_t = kappa_t + v + z_t   (Sigma_eps = I)
      DO qq = 1, Q
          eta(idx1 + qq - 1) = kappa(idx1 + qq - 1) + v(qq) + z(idx1 + qq - 1)
      END DO
    END DO

  END SUBROUTINE eta_from_z_ar1



 END MODULE utilities

 MODULE shared
   IMPLICIT NONE

   INTEGER, SAVE :: nvar, nsub, nobs, Q, n, ntre, nquad, npm
   INTEGER, SAVE :: p_fix, p_ran

   INTEGER, SAVE :: re_item_flag
   INTEGER, SAVE :: theta_time_flag
   INTEGER, SAVE :: tau_invar_flag = 1
   INTEGER, SAVE :: gamma_invar_flag = 1
   INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: cfa_gamma_time_free
   INTEGER, SAVE :: lambda_invar_flag = 1
   INTEGER, SAVE :: theta_invar_flag = 1

   ! ===== NUOVO: modalità integrazione =====
  ! 0 = integra su eta (attuale)
  ! 1 = integra su b (random effects), epsilon=0
  INTEGER, SAVE :: integ_mode = 0

   ! ===== NUOVO: fattori correlati nel processo =====
   INTEGER, SAVE :: dimB    ! = Q * p_ran

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: L_b     ! (dimB,dimB)
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: Sigma_b ! (dimB,dimB)


   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: Sigma_eps ! (Q,Q)

   ! dati
   INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: nmes      ! (nsub,nvar)
   INTEGER, DIMENSION(:),   ALLOCATABLE, SAVE :: c, cm1    ! (nvar)
   INTEGER, DIMENSION(:),   ALLOCATABLE, SAVE :: Y         ! (nobs)
   INTEGER, DIMENSION(:),   ALLOCATABLE, SAVE :: nitem     ! (nobs) indice tempo 1..n

   ! soglie
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: tau
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: tau_long

   ! misurazione
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: lambda   ! (nvar,Q)
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: lambda_long
   DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, SAVE :: gamma    ! (nvar)
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: gamma_long
   DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, SAVE :: theta_k  ! (nvar)
   DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, SAVE :: invsd_k  ! (nvar)
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: invsd_long

   ! fixed effects del processo
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: beta     ! (Q,p_fix)

   ! random effects sugli item a_{ik}
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: sigmak
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: sqrt_sigmak

   ! design del processo
   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: X_eta ! (nsub,n,Q,p_fix)
   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: Z_eta ! (nsub,n,Q,p_ran)

   ! growth design generalizzato
   INTEGER, SAVE :: growth_delta_x_col = 0
   INTEGER, SAVE :: growth_delta_z_col = 0
   ! p_fix0 < 0 in the C ABI means that X/Z/delta are common to all subjects.
   INTEGER, SAVE :: growth_subject_specific_design = 1
   INTEGER, SAVE :: growth_design_nsub = 0

   ! growth identification/configuration
   INTEGER, SAVE :: growth_ident_flag = 0
   INTEGER, SAVE :: growth_threshold_id_flag = 0
   INTEGER, SAVE :: growth_delta_mode = 0
   INTEGER, SAVE :: growth_rho_fixed_flag = 0
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: growth_rho_fixed
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: growth_zero_time
   LOGICAL, SAVE :: growth_options_ready = .FALSE.
   INTEGER, SAVE :: growth_options_nvar = 0, growth_options_n = 0, growth_options_Q = 0
   INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: growth_gamma_time_free
   INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: growth_gamma_dep
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: growth_x_slope_col
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: growth_marker_item
   INTEGER, SAVE :: growth_options_pfix = 0
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: growth_delta_mask ! (growth_design_nsub,n,Q)
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: kappa_sub           ! (n*Q,nsub)
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: L_phi_sub         ! (n*Q,n*Q,nsub)
   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: growth_tt_mu0_grid

   ! nodi integrazione: dimensione = nquad * dimZ, dimZ = Q*n (+ nvar)
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: nodes
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: weights

     INTEGER, SAVE :: cache_i_eta = -1
  INTEGER, SAVE :: cache_dimEta = -1
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: m_eta_cache(:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: L_V_cache(:,:)

  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Zbig_cache(:,:)   ! (dimEta, dimB)

  ! agh twotier
  LOGICAL :: agh_tt_grid_ready = .FALSE.
  INTEGER :: agh_tt_dim_cached = -1
  INTEGER :: agh_tt_nquad_cached = -1
  INTEGER :: agh_tt_nsub_cached = -1
  INTEGER :: agh_tt_fit_id_cached = -1

  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: agh_tt_z
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: agh_tt_modes
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: agh_tt_logconst


  ! agh per long
  INTEGER :: agh_fit_id_cached = -999999999
  LOGICAL :: agh_grid_ready = .FALSE.
INTEGER :: agh_dimZ_cached = -1
INTEGER :: agh_nquad_cached = -1
INTEGER :: agh_nsub_cached = -1


DOUBLE PRECISION, ALLOCATABLE :: agh_z(:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: agh_logconst(:)
DOUBLE PRECISION, ALLOCATABLE :: agh_modes(:,:)

    ! ===== NUOVO: come integrare i random item a_ik =====
  ! 0 = joint dentro al QMC (vecchio comportamento: dimZ = dimLat + nvar)
  ! 1 = GH univariato per ogni item (QMC solo sulle latenti: dimZ = dimLat)
  INTEGER, SAVE :: re_item_integ_mode = 0

  ! GH 1D per a_ik ~ N(0, sigmak)
  INTEGER, SAVE :: nquad_a = 0
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: nodes_a
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: weights_a

  ! two-tier
  INTEGER :: nquad_eta, nquad_re
  INTEGER :: cache_nquad_eta, cache_nquad_re

  DOUBLE PRECISION, ALLOCATABLE :: nodes_eta(:), weights_eta(:), logweights_eta(:)
  DOUBLE PRECISION, ALLOCATABLE :: nodes_re(:),  weights_re(:),  logweights_re(:)

  ! CFA
  ! misurazione
  ! in MODULE shared
  ! --- CFA longitudinale multidimensionale ---
  INTEGER, SAVE :: dim_eta   ! = n * Q (lo setti nella loglik)

  DOUBLE PRECISION, ALLOCATABLE, SAVE :: kappa(:)      ! (dim_eta)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: L_phi(:,:)    ! (dim_eta, dim_eta)

  INTEGER, ALLOCATABLE, SAVE :: loadmask(:,:)          ! (nvar, Q)

  DOUBLE PRECISION, DIMENSION(:,:),    ALLOCATABLE, SAVE :: theta_long ! opzione per varianze varianti nel tempo

  ! CFA longitudinal identification / covariates
  INTEGER, SAVE :: cfa_p_fix = 0
  INTEGER, SAVE :: cfa_threshold_id_flag = 0
  INTEGER, SAVE :: cfa_mean_ident_flag = 0
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: cfa_marker_item
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: cfa_X
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: cfa_beta


  ! Cholesky dei random effects b_{iq}: L_chol(q, p_ran, p_ran)
  DOUBLE PRECISION, DIMENSION(:,:,:),  ALLOCATABLE, SAVE :: L_chol   ! (Q,p_ran,p_ran)

  INTEGER, SAVE :: phi_struct_flag


  DOUBLE PRECISION, ALLOCATABLE :: A_ar(:,:)      ! Q x Q
  DOUBLE PRECISION, ALLOCATABLE :: sqrt_sig1(:)   ! Q

  LOGICAL :: is_init_long = .FALSE.

  INTEGER, ALLOCATABLE :: subj_start(:)
  INTEGER, ALLOCATABLE :: kk_base(:)
  INTEGER, ALLOCATABLE :: nload_item(:)
  INTEGER, ALLOCATABLE :: load_idx(:,:)

  DOUBLE PRECISION, ALLOCATABLE :: logweights(:)


 LOGICAL :: cache_ready = .FALSE.
  INTEGER :: cache_fit_id = -1
  INTEGER :: cache_nvar = -1, cache_nsub = -1, cache_nobs = -1
  INTEGER :: cache_n = -1, cache_Q = -1, cache_nquad = -1

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: lambda_act


    ! GH two-tier fast precomputation
  LOGICAL, SAVE :: tt_precomp_ready = .FALSE.
  INTEGER, SAVE :: tt_precomp_dimEta = -1
  INTEGER, SAVE :: tt_precomp_nquad_eta = -1
  INTEGER, SAVE :: tt_precomp_nvar = -1
  INTEGER, SAVE :: tt_precomp_n = -1
  INTEGER, SAVE :: tt_precomp_Q = -1
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE, SAVE :: tt_eta_grid
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: tt_mu0_grid

  ! AGH two-tier equivalent to full joint AGH
  LOGICAL, SAVE :: agh_tt_fe_grid_ready = .FALSE.
  INTEGER, SAVE :: agh_tt_fe_dimEta_cached = -1
  INTEGER, SAVE :: agh_tt_fe_nvar_cached = -1
  INTEGER, SAVE :: agh_tt_fe_nquad_eta_cached = -1
  INTEGER, SAVE :: agh_tt_fe_nquad_re_cached = -1
  INTEGER, SAVE :: agh_tt_fe_nsub_cached = -1
  INTEGER, SAVE :: agh_tt_fe_fit_id_cached = -1

  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: agh_tt_fe_zeta_eta
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: agh_tt_fe_alpha
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE, SAVE :: agh_tt_fe_sdev
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE, SAVE :: agh_tt_fe_logconst
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE, SAVE :: agh_tt_fe_modes
END MODULE shared

MODULE datindiv
  IMPLICIT NONE
  INTEGER :: nmescur
END MODULE datindiv


MODULE core
  USE shared
  USE utilities
  USE datindiv
  IMPLICIT NONE

CONTAINS

SUBROUTINE set_cfa_long_options(marker_item0, Q0, nvar0, nsub0, n0, p_fix0, X0, ok0) &
     BIND(C, name="set_cfa_long_options")
  USE, INTRINSIC :: ISO_C_BINDING
  USE shared
  IMPLICIT NONE
  INTEGER(C_INT), INTENT(IN) :: Q0, nvar0, nsub0, n0, p_fix0
  INTEGER(C_INT), DIMENSION(Q0), INTENT(IN) :: marker_item0
  REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: X0
  INTEGER(C_INT), INTENT(OUT) :: ok0
  INTEGER :: nX

  ok0 = 0_C_INT
  IF (Q0 <= 0 .OR. nvar0 <= 0 .OR. nsub0 <= 0 .OR. n0 <= 0 .OR. p_fix0 < 0) RETURN
  IF (ANY(marker_item0 < 1_C_INT) .OR. ANY(marker_item0 > nvar0)) RETURN

  IF (ALLOCATED(cfa_marker_item)) DEALLOCATE(cfa_marker_item)
  ALLOCATE(cfa_marker_item(Q0))
  cfa_marker_item = INT(marker_item0)
  cfa_p_fix = INT(p_fix0)

   IF (ALLOCATED(cfa_gamma_time_free)) DEALLOCATE(cfa_gamma_time_free)
   ALLOCATE(cfa_gamma_time_free(nvar0,n0))
   cfa_gamma_time_free = 0

  IF (ALLOCATED(cfa_X)) DEALLOCATE(cfa_X)
  IF (cfa_p_fix > 0) THEN
     nX = nsub0 * n0 * Q0 * cfa_p_fix
     ALLOCATE(cfa_X(nsub0,n0,Q0,cfa_p_fix))
     cfa_X = RESHAPE(X0(1:nX), (/ nsub0, n0, Q0, cfa_p_fix /))
  END IF

  ok0 = 1_C_INT
END SUBROUTINE set_cfa_long_options

SUBROUTINE set_cfa_long_gamma_time_free(gamma_time_free0, nvar0, n0, ok0) &
     BIND(C, name="set_cfa_long_gamma_time_free")
  USE, INTRINSIC :: ISO_C_BINDING
  USE shared
  IMPLICIT NONE
  INTEGER(C_INT), INTENT(IN) :: nvar0, n0
  INTEGER(C_INT), DIMENSION(nvar0,n0), INTENT(IN) :: gamma_time_free0
  INTEGER(C_INT), INTENT(OUT) :: ok0

  ok0 = 0_C_INT
  IF (nvar0 < 1 .OR. n0 < 1) RETURN
  IF (MINVAL(gamma_time_free0) < 0_C_INT .OR. MAXVAL(gamma_time_free0) > 1_C_INT) RETURN
  IF (ANY(gamma_time_free0(:,1) /= 0_C_INT)) RETURN
  IF (.NOT. ALLOCATED(cfa_gamma_time_free)) THEN
     ALLOCATE(cfa_gamma_time_free(nvar0,n0))
  ELSE IF (SIZE(cfa_gamma_time_free,1) /= nvar0 .OR. SIZE(cfa_gamma_time_free,2) /= n0) THEN
     DEALLOCATE(cfa_gamma_time_free)
     ALLOCATE(cfa_gamma_time_free(nvar0,n0))
  END IF

  IF (ANY(cfa_gamma_time_free /= INT(gamma_time_free0))) THEN
     agh_tt_fe_grid_ready = .FALSE.
  END IF
  cfa_gamma_time_free = INT(gamma_time_free0)
  ok0 = 1_C_INT
END SUBROUTINE set_cfa_long_gamma_time_free

SUBROUTINE set_cfa_long_threshold_id(threshold_id0, ok0) &
     BIND(C, name="set_cfa_long_threshold_id")
  USE, INTRINSIC :: ISO_C_BINDING
  USE shared
  IMPLICIT NONE
  INTEGER(C_INT), INTENT(IN) :: threshold_id0
  INTEGER(C_INT), INTENT(OUT) :: ok0

  ok0 = 0_C_INT
  IF (threshold_id0 /= 0_C_INT .AND. threshold_id0 /= 1_C_INT) RETURN

  cfa_threshold_id_flag = INT(threshold_id0)

  agh_grid_ready = .FALSE.
  agh_tt_grid_ready = .FALSE.
  agh_tt_fe_grid_ready = .FALSE.
  tt_precomp_ready = .FALSE.

  ok0 = 1_C_INT
END SUBROUTINE set_cfa_long_threshold_id

SUBROUTINE set_cfa_long_mean_ident(mean_ident0, ok0) &
     BIND(C, name="set_cfa_long_mean_ident")
  USE, INTRINSIC :: ISO_C_BINDING
  USE shared
  IMPLICIT NONE
  INTEGER(C_INT), INTENT(IN) :: mean_ident0
  INTEGER(C_INT), INTENT(OUT) :: ok0

  ok0 = 0_C_INT
  IF (mean_ident0 /= 0_C_INT .AND. mean_ident0 /= 1_C_INT) RETURN

  cfa_mean_ident_flag = INT(mean_ident0)

  agh_grid_ready = .FALSE.
  agh_tt_grid_ready = .FALSE.
  agh_tt_fe_grid_ready = .FALSE.
  tt_precomp_ready = .FALSE.

  ok0 = 1_C_INT
END SUBROUTINE set_cfa_long_mean_ident

SUBROUTINE cfa_eta_mean(i, mean_eta)
  USE shared, ONLY: n, Q, p_fix, kappa, kappa_sub, cfa_X, cfa_beta
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(n*Q), INTENT(OUT) :: mean_eta
  INTEGER :: t, qidx, r, idx

  mean_eta = kappa
  IF (ALLOCATED(kappa_sub)) RETURN
  IF (p_fix <= 0) RETURN
  IF (.NOT. ALLOCATED(cfa_X)) THEN
     CALL rexit('ERROR (cfa_eta_mean): CFA covariate design not initialized')
  END IF
  IF (.NOT. ALLOCATED(cfa_beta)) THEN
     CALL rexit('ERROR (cfa_eta_mean): CFA beta not allocated')
  END IF

  DO t = 1, n
     DO qidx = 1, Q
        idx = (t-1)*Q + qidx
        DO r = 1, p_fix
           mean_eta(idx) = mean_eta(idx) + cfa_X(i,t,qidx,r) * cfa_beta(qidx,r)
        END DO
     END DO
  END DO
END SUBROUTINE cfa_eta_mean

SUBROUTINE cfa_fill_tau_item(m, u, monotone_tau0, out)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: m, monotone_tau0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: u
  DOUBLE PRECISION, DIMENSION(m), INTENT(OUT) :: out
  INTEGER :: nf, h
  DOUBLE PRECISION :: prev

  IF (m < 3) THEN
     CALL rexit('ERROR (CFA thresholds): identification requires at least 4 response categories')
  END IF

  nf = m - 2
  out = 0.d0

  SELECT CASE (cfa_threshold_id_flag)

  CASE (0)
     out(2) = 0.d0
     out(3) = 1.d0

     IF (monotone_tau0 == 0) THEN
        out(1) = u(1)
        IF (m > 3) out(4:m) = u(2:nf)
     ELSE
        out(1) = -EXP(u(1))
        prev = 1.d0
        DO h = 2, nf
           prev = prev + EXP(u(h))
           out(h+2) = prev
        END DO
     END IF

  CASE (1)
     out(1) = 0.d0
     out(2) = 1.d0

     IF (monotone_tau0 == 0) THEN
        DO h = 1, nf
           out(h+2) = u(h)
        END DO
     ELSE
        prev = 1.d0
        DO h = 1, nf
           prev = prev + EXP(u(h))
           out(h+2) = prev
        END DO
     END IF

  CASE DEFAULT
     CALL rexit('ERROR (CFA thresholds): invalid cfa_threshold_id_flag')
  END SELECT
END SUBROUTINE cfa_fill_tau_item

SUBROUTINE ensure_cfa_beta_storage()
  USE shared, ONLY: Q, p_fix, cfa_beta
  IMPLICIT NONE
  LOGICAL :: bad
  bad = .NOT. ALLOCATED(cfa_beta)
  IF (.NOT. bad) bad = SIZE(cfa_beta,1) /= Q .OR. &
                       SIZE(cfa_beta,2) /= MAX(1,p_fix)
  IF (bad) THEN
     IF (ALLOCATED(cfa_beta)) DEALLOCATE(cfa_beta)
     ALLOCATE(cfa_beta(Q,MAX(1,p_fix)))
  END IF
  cfa_beta = 0.d0
END SUBROUTINE ensure_cfa_beta_storage

DOUBLE PRECISION FUNCTION tau_at(idx, t_idx)
  USE shared, ONLY: tau, tau_long
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: idx, t_idx

  IF (ALLOCATED(tau_long)) THEN
     tau_at = tau_long(idx, t_idx)
  ELSE
     tau_at = tau(idx)
  END IF
END FUNCTION tau_at

DOUBLE PRECISION FUNCTION gamma_at(k_item, t_idx)
  USE shared, ONLY: gamma, gamma_long
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: k_item, t_idx

  IF (ALLOCATED(gamma_long)) THEN
     gamma_at = gamma_long(k_item, t_idx)
  ELSE
     gamma_at = gamma(k_item)
  END IF
END FUNCTION gamma_at

DOUBLE PRECISION FUNCTION lambda_at(k_item, qq, t_idx)
  USE shared, ONLY: lambda, lambda_long
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: k_item, qq, t_idx

  IF (ALLOCATED(lambda_long)) THEN
     lambda_at = lambda_long(k_item, qq, t_idx)
  ELSE
     lambda_at = lambda(k_item, qq)
  END IF
END FUNCTION lambda_at

DOUBLE PRECISION FUNCTION theta_at(k_item, t_idx)
  USE shared, ONLY: theta_k, theta_long
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: k_item, t_idx

  IF (ALLOCATED(theta_long)) THEN
     theta_at = theta_long(k_item, t_idx)
  ELSE
     theta_at = theta_k(k_item)
  END IF
END FUNCTION theta_at

DOUBLE PRECISION FUNCTION vrais_multo_i(i, sm)
  USE shared
  USE utilities
  USE datindiv
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, sm

  DOUBLE PRECISION, PARAMETER :: LIK_FLOOR = 1.0D-300
  DOUBLE PRECISION, PARAMETER :: LOG_LIK_FLOOR = LOG(LIK_FLOOR)

  INTEGER :: l, k_item, j_mes, idx_obs, kk, yijk, tuk
  INTEGER :: qq, r, row, col, t_idx, d
  INTEGER :: dimZ, baseZ, dimEta, dimLat, base_items
  INTEGER :: info

  ! GH1D locals
  INTEGER :: mnode, mobs, idx0, idx_tmp
  DOUBLE PRECISION :: maxlw_a, sumexp_a, lw_a, logI_k, a_shift

  DOUBLE PRECISION :: mu_ijk
  DOUBLE PRECISION :: z_item, tausup, tauinf
  DOUBLE PRECISION :: invsd_lat_k
  DOUBLE PRECISION :: Phi_sup, Phi_inf, p
  DOUBLE PRECISION :: loglik_l, lw, maxlw, sumexp
  LOGICAL :: has_any, has_any_a

  DOUBLE PRECISION, ALLOCATABLE :: a(:)
  DOUBLE PRECISION, ALLOCATABLE :: eta_vec(:), m_eta(:)

  ! mode 0
  DOUBLE PRECISION, ALLOCATABLE :: z_lat_eta(:)
  DOUBLE PRECISION, ALLOCATABLE :: V_eta(:,:), L_V(:,:), tmp(:,:)

  ! mode 1
  DOUBLE PRECISION, ALLOCATABLE :: z_lat_b(:), b_vec(:)

  DOUBLE PRECISION, ALLOCATABLE :: Zbig(:,:)

  ! GH1D buffers per item (max = n osservazioni per item)
  DOUBLE PRECISION, ALLOCATABLE :: mu0(:)
  INTEGER,          ALLOCATABLE :: ybuf(:), tbuf(:), tukbuf(:)

  IF (sm < 0) THEN
    vrais_multo_i = -1.d9
    RETURN
  END IF

  dimEta = Q*n

  IF (integ_mode == 0) THEN
    dimLat = dimEta
  ELSE
    dimLat = dimB
  END IF

  ! dimZ (QMC) dipende anche da come integriamo gli item
  IF (re_item_flag == 1 .AND. re_item_integ_mode == 0) THEN
    dimZ = dimLat + nvar
  ELSE
    dimZ = dimLat
  END IF

  base_items = dimLat

  ALLOCATE(a(nvar))
  ALLOCATE(eta_vec(dimEta), m_eta(dimEta))
  ALLOCATE(Zbig(dimEta, dimB))

  a(:) = 0.d0

  ! m_eta = X beta
  m_eta(:) = 0.d0
  DO t_idx = 1, n
    DO qq = 1, Q
      row = (t_idx-1)*Q + qq
      DO r = 1, p_fix
        m_eta(row) = m_eta(row) + X_eta(i,t_idx,qq,r) * beta(qq,r)
      END DO
    END DO
  END DO

  ! Zbig (sempre utile se integ_mode=1 o per costruire V_eta in mode 0)
  Zbig(:,:) = 0.d0
  DO t_idx = 1, n
    DO qq = 1, Q
      row = (t_idx-1)*Q + qq
      DO r = 1, p_ran
        Zbig(row, (qq-1)*p_ran + r) = Z_eta(i, t_idx, qq, r)
      END DO
    END DO
  END DO

  !======================
  ! Precomputazioni parte latente
  !======================
  IF (integ_mode == 0) THEN
    ALLOCATE(V_eta(dimEta,dimEta), L_V(dimEta,dimEta), tmp(dimEta,dimB))
    V_eta(:,:) = 0.d0
    tmp(:,:)   = 0.d0

    CALL DGEMM('N','N', dimEta, dimB, dimB, 1.d0, Zbig, dimEta, Sigma_b, dimB, 0.d0, tmp, dimEta)
    CALL DGEMM('N','T', dimEta, dimEta, dimB, 1.d0, tmp,  dimEta, Zbig,   dimEta, 0.d0, V_eta, dimEta)

    DO t_idx = 1, n
      DO qq = 1, Q
        row = (t_idx-1)*Q + qq
        V_eta(row,row) = V_eta(row,row) + Sigma_eps(qq,qq)
      END DO
    END DO

    L_V(:,:) = V_eta(:,:)
    CALL DPOTRF('L', dimEta, L_V, dimEta, info)
    IF (info /= 0) THEN
      vrais_multo_i = -1.d9
      DEALLOCATE(a, eta_vec, m_eta, Zbig, V_eta, L_V, tmp)
      RETURN
    END IF

    ALLOCATE(z_lat_eta(dimEta))
  ELSE
    ALLOCATE(z_lat_b(dimB), b_vec(dimB))
  END IF

  ! Buffers GH1D (solo se richiesti)
  IF (re_item_flag == 1 .AND. re_item_integ_mode == 1) THEN
    ALLOCATE(mu0(n), ybuf(n), tbuf(n), tukbuf(n))
  END IF

  !======================
  ! LOG-SUM-EXP su nodi QMC
  !======================
  has_any = .FALSE.
  maxlw   = -HUGE(1.d0)
  sumexp  = 0.d0

  DO l = 1, nquad
    IF (weights(l) <= 0.d0) CYCLE
    baseZ = (l-1)*dimZ

    ! ---- costruisci eta_vec dato il nodo QMC (solo latenti) ----
    IF (integ_mode == 0) THEN
      DO d = 1, dimEta
        z_lat_eta(d) = nodes(baseZ + d)
      END DO

      eta_vec(:) = m_eta(:)
      DO row = 1, dimEta
        DO col = 1, row
          eta_vec(row) = eta_vec(row) + L_V(row,col) * z_lat_eta(col)
        END DO
      END DO

    ELSE
      DO d = 1, dimB
        z_lat_b(d) = nodes(baseZ + d)
      END DO

      ! b = L_b z_b
      b_vec(:) = 0.d0
      DO row = 1, dimB
        DO col = 1, row
          b_vec(row) = b_vec(row) + L_b(row,col) * z_lat_b(col)
        END DO
      END DO

      ! eta = m + Zbig b
      eta_vec(:) = m_eta(:)
      DO row = 1, dimEta
        DO col = 1, dimB
          eta_vec(row) = eta_vec(row) + Zbig(row,col) * b_vec(col)
        END DO
      END DO
    END IF

    ! ---- loglik condizionata al nodo QMC ----
    loglik_l = 0.d0
    kk       = 0
    idx_obs  = nmescur

    DO k_item = 1, nvar

      IF (nmes(i,k_item) == 0) THEN
        kk = kk + (c(k_item) - 1)
        CYCLE
      END IF

      invsd_lat_k = invsd_k(k_item)

      ! ==========================
      ! CASO A: re_item joint (vecchio) oppure re_item assente
      ! ==========================
      IF (re_item_flag == 0 .OR. re_item_integ_mode == 0) THEN

        ! a_k dal nodo QMC solo se re_item_flag=1
        IF (re_item_flag == 1) THEN
          z_item    = nodes(baseZ + base_items + k_item)
          a(k_item) = sqrt_sigmak(k_item) * z_item
        ELSE
          a(k_item) = 0.d0
        END IF

        DO j_mes = 1, nmes(i,k_item)
          idx_obs = idx_obs + 1
          yijk    = Y(idx_obs)
          t_idx   = nitem(idx_obs)

          mu_ijk = gamma(k_item) + a(k_item)
          DO qq = 1, Q
            mu_ijk = mu_ijk + lambda(k_item,qq) * eta_vec((t_idx-1)*Q + qq)
          END DO

          tuk = yijk + kk

          IF (yijk == 1) THEN
            tausup  = (tau(tuk) - mu_ijk) * invsd_lat_k
            p       = alnorm_fast(tausup, .FALSE.)
          ELSE IF (yijk == c(k_item)) THEN
            tauinf  = (tau(tuk-1) - mu_ijk) * invsd_lat_k
            p       = 1.d0 - alnorm_fast(tauinf, .FALSE.)
          ELSE
            tausup  = (tau(tuk)   - mu_ijk) * invsd_lat_k
            tauinf  = (tau(tuk-1) - mu_ijk) * invsd_lat_k
            p       = alnorm_fast(tausup, .FALSE.) - alnorm_fast(tauinf, .FALSE.)
          END IF

          IF (p /= p) THEN
            loglik_l = loglik_l + LOG_LIK_FLOOR
          ELSE IF (p <= LIK_FLOOR) THEN
            loglik_l = loglik_l + LOG_LIK_FLOOR
          ELSE
            loglik_l = loglik_l + LOG(p)
          END IF
        END DO

      ! ==========================
      ! CASO B: re_item integrato via GH1D indipendente (NUOVO)
      ! ==========================
      ELSE
        mobs = nmes(i,k_item)
        idx0 = idx_obs  ! punta all'ultima obs già consumata prima di questo item

        ! precompute mu0(j) = gamma + lambda*eta (senza a)
        DO j_mes = 1, mobs
          idx_tmp     = idx0 + j_mes
          ybuf(j_mes) = Y(idx_tmp)
          tbuf(j_mes) = nitem(idx_tmp)
          tukbuf(j_mes) = ybuf(j_mes) + kk

          mu_ijk = gamma(k_item)
          DO qq = 1, Q
            mu_ijk = mu_ijk + lambda(k_item,qq) * eta_vec((tbuf(j_mes)-1)*Q + qq)
          END DO
          mu0(j_mes) = mu_ijk
        END DO

        ! integra su a con GH: log I_k = log sum_m w_m * exp(loglik_m)
        has_any_a = .FALSE.
        maxlw_a   = -HUGE(1.d0)
        sumexp_a  = 0.d0

        DO mnode = 1, nquad_a
          IF (weights_a(mnode) <= 0.d0) CYCLE

          a_shift = sqrt_sigmak(k_item) * nodes_a(mnode)

          lw_a = 0.d0
          DO j_mes = 1, mobs
            yijk  = ybuf(j_mes)
            t_idx = tbuf(j_mes)
            tuk   = tukbuf(j_mes)

            mu_ijk = mu0(j_mes) + a_shift

            IF (yijk == 1) THEN
              tausup  = (tau(tuk) - mu_ijk) * invsd_lat_k
              p       = alnorm_fast(tausup, .FALSE.)
            ELSE IF (yijk == c(k_item)) THEN
              tauinf  = (tau(tuk-1) - mu_ijk) * invsd_lat_k
              p       = 1.d0 - alnorm_fast(tauinf, .FALSE.)
            ELSE
              tausup  = (tau(tuk)   - mu_ijk) * invsd_lat_k
              tauinf  = (tau(tuk-1) - mu_ijk) * invsd_lat_k
              p       = alnorm_fast(tausup, .FALSE.) - alnorm_fast(tauinf, .FALSE.)
            END IF

            IF (p /= p) THEN
              lw_a = lw_a + LOG_LIK_FLOOR
            ELSE IF (p <= LIK_FLOOR) THEN
              lw_a = lw_a + LOG_LIK_FLOOR
            ELSE
              lw_a = lw_a + LOG(p)
            END IF
          END DO

          lw_a = LOG(weights_a(mnode)) + lw_a

          IF (.NOT. has_any_a) THEN
            maxlw_a   = lw_a
            sumexp_a  = 1.d0
            has_any_a = .TRUE.
          ELSE
            IF (lw_a > maxlw_a) THEN
              sumexp_a = sumexp_a * EXP(maxlw_a - lw_a) + 1.d0
              maxlw_a  = lw_a
            ELSE
              sumexp_a = sumexp_a + EXP(lw_a - maxlw_a)
            END IF
          END IF
        END DO

        IF (.NOT. has_any_a .OR. sumexp_a <= 0.d0) THEN
          logI_k = LOG_LIK_FLOOR
        ELSE
          logI_k = maxlw_a + LOG(sumexp_a)
        END IF

        loglik_l = loglik_l + logI_k
        idx_obs  = idx_obs + mobs
      END IF

      kk = kk + (c(k_item) - 1)
    END DO

    lw = LOG(weights(l)) + loglik_l

    IF (.NOT. has_any) THEN
      maxlw   = lw
      sumexp  = 1.d0
      has_any = .TRUE.
    ELSE
      IF (lw > maxlw) THEN
        sumexp = sumexp * EXP(maxlw - lw) + 1.d0
        maxlw  = lw
      ELSE
        sumexp = sumexp + EXP(lw - maxlw)
      END IF
    END IF

  END DO

  IF (.NOT. has_any) THEN
    vrais_multo_i = -1.d9
  ELSE
    IF (sumexp <= 0.d0) THEN
      vrais_multo_i = -1.d9
    ELSE
      vrais_multo_i = maxlw + LOG(sumexp)
    END IF
  END IF

  IF (re_item_flag == 1 .AND. re_item_integ_mode == 1) THEN
    DEALLOCATE(mu0, ybuf, tbuf, tukbuf)
  END IF

  IF (integ_mode == 0) THEN
    DEALLOCATE(z_lat_eta, V_eta, L_V, tmp)
  ELSE
    DEALLOCATE(z_lat_b, b_vec)
  END IF

  DEALLOCATE(a, eta_vec, m_eta, Zbig)

END FUNCTION vrais_multo_i


  DOUBLE PRECISION FUNCTION vrais_multo()
    IMPLICIT NONE
    INTEGER :: i, sm
    DOUBLE PRECISION :: temp

    nmescur     = 0
    vrais_multo = 0.d0

    DO i = 1, nsub
      sm   = SUM(nmes(i,:))
      temp = vrais_multo_i(i, sm)

      IF (temp == -1.d9) THEN
        vrais_multo = -1.d9
        RETURN
      END IF

      vrais_multo = vrais_multo + temp
      nmescur = nmescur + sm
    END DO

  END FUNCTION vrais_multo



 SUBROUTINE loglik(b0, Y0, nvar0, nsub0, nobs0, Q0, n0, &
                  nmes0, c0, cm10, ntre0, nquad0, &
                  nitem0, nodes0, weights0, npm0, &
                  p_fix0, p_ran0, X0, Z0, fac_item0, &
                  ident_type0, re_item0, integ_mode0, &
                  re_item_integ0, nquad_a0, nodes_a0, weights_a0, &
                  vrais) &
     BIND(C, name="loglik")

  USE shared
  USE utilities
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, Q0, n0
  INTEGER, INTENT(IN) :: ntre0, nquad0, npm0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0
  INTEGER, INTENT(IN) :: ident_type0, re_item0
  INTEGER, INTENT(IN) :: integ_mode0
  INTEGER, INTENT(IN) :: re_item_integ0, nquad_a0

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0),       INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: fac_item0

  DOUBLE PRECISION, DIMENSION(*),      INTENT(IN) :: nodes0
  DOUBLE PRECISION, DIMENSION(nquad0), INTENT(IN) :: weights0
  DOUBLE PRECISION, DIMENSION(npm0),   INTENT(IN) :: b0

  DOUBLE PRECISION, DIMENSION(*),      INTENT(IN) :: nodes_a0
  DOUBLE PRECISION, DIMENSION(*),      INTENT(IN) :: weights_a0

  DOUBLE PRECISION, INTENT(IN), TARGET :: X0(*)
  DOUBLE PRECISION, INTENT(IN), TARGET :: Z0(*)

  DOUBLE PRECISION, INTENT(OUT) :: vrais

  INTEGER :: k, qq, rr, rowi, coli, jj
  INTEGER :: dimZ, dimEta, dimLat, nX, nZ
  INTEGER :: npar_eps, npar_Lb
  INTEGER :: ntheta_here, nload, ngamma
  INTEGER :: pos, pos_tau
  INTEGER :: expected_npm

  INTEGER, ALLOCATABLE :: cm3(:)
  DOUBLE PRECISION, ALLOCATABLE :: ta(:)

  INTEGER, ALLOCATABLE :: anchor(:), Jq(:)
  LOGICAL, ALLOCATABLE :: is_anchor(:)

  ! ------------------ globali ------------------
  nvar  = nvar0
  nsub  = nsub0
  nobs  = nobs0
  Q     = Q0
  n     = n0
  ntre  = ntre0
  nquad = nquad0
  npm   = npm0
  p_fix = p_fix0
  p_ran = p_ran0

  dimEta = Q*n
  dimB   = Q*p_ran

  re_item_flag       = re_item0
  integ_mode         = integ_mode0
  re_item_integ_mode = re_item_integ0

  ! dimLat
  IF (integ_mode == 0) THEN
    dimLat = dimEta
  ELSE
    dimLat = dimB
  END IF

  ! dimZ (QMC) dipende da re_item_integ_mode
  IF (re_item_flag == 1 .AND. re_item_integ_mode == 0) THEN
    dimZ = dimLat + nvar
  ELSE
    dimZ = dimLat
  END IF

  ! Var_eps: se integ_mode=1 => epsilon=0 (nessun parametro Var_eps)
  IF (integ_mode == 1) THEN
    npar_eps = 0
  ELSE
    IF (ident_type0 == 1) THEN
      npar_eps = 0
    ELSE
      npar_eps = Q
    END IF
  END IF

  npar_Lb = dimB*(dimB+1)/2

  ! ------------------ GH1D per a_ik (solo se serve) ------------------
  IF (ALLOCATED(nodes_a))   DEALLOCATE(nodes_a)
  IF (ALLOCATED(weights_a)) DEALLOCATE(weights_a)

  IF (re_item_flag == 1 .AND. re_item_integ_mode == 1) THEN
    nquad_a = nquad_a0
    IF (nquad_a <= 0) THEN
      vrais = 1.d9
      RETURN
    END IF
    ALLOCATE(nodes_a(nquad_a), weights_a(nquad_a))
    nodes_a(1:nquad_a)   = nodes_a0(1:nquad_a)
    weights_a(1:nquad_a) = weights_a0(1:nquad_a)
  ELSE
    nquad_a = 0
    ALLOCATE(nodes_a(1), weights_a(1))
    nodes_a(1)   = 0.d0
    weights_a(1) = 0.d0
  END IF

  ! ------------------ anchor/Jq ------------------
  ALLOCATE(anchor(Q), Jq(Q), is_anchor(nvar))
  anchor = 0
  Jq = 0
  is_anchor = .FALSE.

  DO k = 1, nvar
    qq = fac_item0(k)
    Jq(qq) = Jq(qq) + 1
    IF (anchor(qq) == 0) anchor(qq) = k
  END DO
  DO qq = 1, Q
    is_anchor(anchor(qq)) = .TRUE.
  END DO

  ! ------------------ reset cache ------------------
  cache_i_eta   = -1
  cache_dimEta  = -1
  IF (ALLOCATED(m_eta_cache)) DEALLOCATE(m_eta_cache)
  IF (ALLOCATED(L_V_cache))   DEALLOCATE(L_V_cache)
  IF (ALLOCATED(Zbig_cache))  DEALLOCATE(Zbig_cache)

  ! ------------------ dati ------------------
  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs), c(nvar), cm1(nvar))
  nmes  = nmes0
  Y     = Y0
  nitem = nitem0
  c     = c0
  cm1   = cm10

  ntheta_here = 0
  DO k = 1, nvar
    IF (cm1(k) > 1) ntheta_here = ntheta_here + 1
  END DO

  ! ------------------ nodi/pesi QMC ------------------
  ALLOCATE(nodes(nquad*dimZ))
  nodes(1:nquad*dimZ) = nodes0(1:nquad*dimZ)
  ALLOCATE(weights(nquad))
  weights(1:nquad) = weights0(1:nquad)

  ! ------------------ parametri ------------------
  ALLOCATE(lambda(nvar,Q))
  ALLOCATE(gamma(nvar))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(beta(Q,p_fix))
  ALLOCATE(L_b(dimB,dimB), Sigma_b(dimB,dimB))
  ALLOCATE(Sigma_eps(Q,Q))
  ALLOCATE(sigmak(nvar), sqrt_sigmak(nvar))
  ALLOCATE(X_eta(nsub,n,Q,p_fix))
  ALLOCATE(Z_eta(nsub,n,Q,p_ran))

  nX = nsub*n*Q*p_fix
  nZ = nsub*n*Q*p_ran
  X_eta = RESHAPE( X0(1:nX), (/ nsub, n, Q, p_fix /) )
  Z_eta = RESHAPE( Z0(1:nZ), (/ nsub, n, Q, p_ran /) )

  ! ------------------ thresholds ------------------
  ALLOCATE(cm3(nvar))
  DO k = 1, nvar
    cm3(k) = cm1(k) - 2
    IF (cm3(k) < 0) cm3(k) = 0
  END DO
  ntre = SUM(cm3)

  IF (ntre > 0) THEN
    ALLOCATE(ta(ntre))
    ta = b0(1:ntre)
  ELSE
    ALLOCATE(ta(1))
    ta = 0.d0
  END IF

  ALLOCATE(tau(SUM(cm1)))

  pos_tau = 0
  jj = 1
  DO k = 1, nvar
    IF (cm1(k) == 1) THEN
      tau(pos_tau+1) = 0.d0
    ELSE
      tau(pos_tau+1) = 0.d0
      tau(pos_tau+2) = 1.d0
      IF (cm1(k) > 2) THEN
        tau(pos_tau+3 : pos_tau+cm1(k)) = ta(jj : jj+cm3(k)-1)
        jj = jj + cm3(k)
      END IF
    END IF
    pos_tau = pos_tau + cm1(k)
  END DO

  pos = ntre

  ! ------------------ lambda ------------------
  lambda = 0.d0

  IF (ident_type0 == 1) THEN
    DO k = 1, nvar
      qq = fac_item0(k)
      pos = pos + 1
      lambda(k,qq) = b0(pos)
    END DO

  ELSE IF (ident_type0 == 0) THEN
    DO qq = 1, Q
      lambda(anchor(qq),qq) = 1.d0
    END DO
    DO qq = 1, Q
      DO k = 1, nvar
        IF (fac_item0(k) == qq .AND. k /= anchor(qq)) THEN
          pos = pos + 1
          lambda(k,qq) = b0(pos)
        END IF
      END DO
    END DO

  ELSE IF (ident_type0 == 2 .OR. ident_type0 == 3) THEN
    DO qq = 1, Q
      DO k = 1, nvar
        IF (fac_item0(k) == qq .AND. k /= anchor(qq)) THEN
          pos = pos + 1
          lambda(k,qq) = b0(pos)
        END IF
      END DO

      IF (ident_type0 == 2) THEN
        lambda(anchor(qq),qq) = 1.d0
      ELSE
        lambda(anchor(qq),qq) = DBLE(Jq(qq))  ! EC-1: somma=Jq -> media=1
      END IF

      DO k = 1, nvar
        IF (fac_item0(k) == qq .AND. k /= anchor(qq)) THEN
          lambda(anchor(qq),qq) = lambda(anchor(qq),qq) - lambda(k,qq)
        END IF
      END DO
    END DO
  END IF

  ! ------------------ gamma ------------------
  gamma = 0.d0

  IF (ident_type0 == 2 .OR. ident_type0 == 3) THEN
    DO qq = 1, Q
      DO k = 1, nvar
        IF (fac_item0(k) == qq .AND. k /= anchor(qq)) THEN
          pos = pos + 1
          gamma(k) = b0(pos)
        END IF
      END DO

      gamma(anchor(qq)) = 0.d0
      DO k = 1, nvar
        IF (fac_item0(k) == qq .AND. k /= anchor(qq)) THEN
          gamma(anchor(qq)) = gamma(anchor(qq)) - gamma(k)
        END IF
      END DO
    END DO
  ELSE
    DO k = 1, nvar
      pos = pos + 1
      gamma(k) = b0(pos)
    END DO
  END IF

  ! ------------------ theta ------------------
  DO k = 1, nvar
    IF (cm1(k) == 1) THEN
      theta_k(k) = 1.d0
      invsd_k(k) = 1.d0
    ELSE
      pos = pos + 1
      theta_k(k) = EXP(b0(pos))
      invsd_k(k) = 1.d0 / SQRT(theta_k(k))
    END IF
  END DO

  ! ------------------ beta ------------------
  DO qq = 1, Q
    DO rr = 1, p_fix
      pos = pos + 1
      beta(qq,rr) = b0(pos)
    END DO
  END DO

  ! ------------------ L_b ------------------
  L_b(:,:) = 0.d0
  DO rr = 1, dimB
    pos = pos + 1
    L_b(rr,rr) = EXP(b0(pos))
  END DO
  DO coli = 1, dimB
    DO rowi = coli+1, dimB
      pos = pos + 1
      L_b(rowi,coli) = b0(pos)
    END DO
  END DO

  ! ------------------ Sigma_eps ------------------
  Sigma_eps(:,:) = 0.d0
  IF (integ_mode == 0) THEN
    IF (ident_type0 == 1) THEN
      DO qq = 1, Q
        Sigma_eps(qq,qq) = 1.d0
      END DO
    ELSE
      DO qq = 1, Q
        pos = pos + 1
        Sigma_eps(qq,qq) = EXP(b0(pos))
      END DO
    END IF
  ELSE
    Sigma_eps(:,:) = 0.d0   ! epsilon=0
  END IF

  ! ------------------ Sigma_b = L_b L_b' ------------------
  Sigma_b(:,:) = 0.d0
  CALL DGEMM('N','T', dimB, dimB, dimB, 1.d0, L_b, dimB, L_b, dimB, 0.d0, Sigma_b, dimB)

  ! ------------------ sigmak ------------------
  IF (re_item0 == 1) THEN
    DO k = 1, nvar
      pos = pos + 1
      sigmak(k)      = EXP(b0(pos))
      sqrt_sigmak(k) = SQRT(sigmak(k))
    END DO
  ELSE
    sigmak(:)      = 0.d0
    sqrt_sigmak(:) = 0.d0
  END IF

  ! ------------------ check conteggi (debug) ------------------
  IF (ident_type0 == 1) THEN
    nload  = nvar
    ngamma = nvar
  ELSE IF (ident_type0 == 0) THEN
    nload  = nvar - Q
    ngamma = nvar
  ELSE
    nload  = nvar - Q
    ngamma = nvar - Q
  END IF

  IF (re_item0 == 1) THEN
    expected_npm = ntre + nload + ngamma + ntheta_here + Q*p_fix + npar_Lb + npar_eps + nvar
  ELSE
    expected_npm = ntre + nload + ngamma + ntheta_here + Q*p_fix + npar_Lb + npar_eps
  END IF
  ! (se vuoi) puoi controllare pos vs expected_npm qui

  ! ------------------ loglik QMC ------------------
  vrais = -vrais_multo()
  ! print *, vrais

  ! ------------------ dealloc ------------------
  DEALLOCATE(nmes, Y, nitem, c, cm1)
  DEALLOCATE(nodes, weights)
  DEALLOCATE(lambda, gamma, theta_k, invsd_k, beta)
  DEALLOCATE(L_b, Sigma_b, Sigma_eps)
  DEALLOCATE(sigmak, sqrt_sigmak)
  DEALLOCATE(tau, cm3, ta)
  DEALLOCATE(X_eta, Z_eta)
  DEALLOCATE(anchor, Jq, is_anchor)

  IF (ALLOCATED(nodes_a))   DEALLOCATE(nodes_a)
  IF (ALLOCATED(weights_a)) DEALLOCATE(weights_a)

END SUBROUTINE loglik








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                   !
!                                   !
!              ADAPTIVE             !
!                                   !
!                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




SUBROUTINE ensure_eta_cache(i, info)
  USE shared
  USE utilities
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: i
  INTEGER, INTENT(OUT) :: info

  INTEGER :: dimEta, t_idx, qq, r, idx_eta
  DOUBLE PRECISION, ALLOCATABLE :: tmp(:,:), V_eta(:,:)

  dimEta = Q*n
  info   = 0

  IF (cache_i_eta == i .AND. cache_dimEta == dimEta .AND. ALLOCATED(m_eta_cache)) RETURN

  cache_i_eta  = i
  cache_dimEta = dimEta

  IF (ALLOCATED(m_eta_cache)) DEALLOCATE(m_eta_cache)
  IF (ALLOCATED(L_V_cache))   DEALLOCATE(L_V_cache)
  IF (ALLOCATED(Zbig_cache))  DEALLOCATE(Zbig_cache)

  ALLOCATE(m_eta_cache(dimEta))
  ALLOCATE(Zbig_cache(dimEta, dimB))

  ! ---- m_eta ----
  m_eta_cache(:) = 0.d0
  DO t_idx = 1, n
    DO qq = 1, Q
      idx_eta = (t_idx-1)*Q + qq
      DO r = 1, p_fix
        m_eta_cache(idx_eta) = m_eta_cache(idx_eta) + X_eta(i,t_idx,qq,r) * beta(qq,r)
      END DO
    END DO
  END DO

  ! ---- Zbig_cache ----
  Zbig_cache(:,:) = 0.d0
  DO t_idx = 1, n
    DO qq = 1, Q
      idx_eta = (t_idx-1)*Q + qq
      DO r = 1, p_ran
        Zbig_cache(idx_eta, (qq-1)*p_ran + r) = Z_eta(i,t_idx,qq,r)
      END DO
    END DO
  END DO

  ! ---- mode 0: chol(V_eta) ----
  IF (integ_mode == 0) THEN
    ALLOCATE(tmp(dimEta, dimB), V_eta(dimEta,dimEta))
    tmp(:,:) = 0.d0
    V_eta(:,:) = 0.d0

    CALL DGEMM('N','N', dimEta, dimB, dimB, 1.d0, Zbig_cache, dimEta, Sigma_b, dimB, 0.d0, tmp, dimEta)
    CALL DGEMM('N','T', dimEta, dimEta, dimB, 1.d0, tmp, dimEta, Zbig_cache, dimEta, 0.d0, V_eta, dimEta)

    DO t_idx = 1, n
      DO qq = 1, Q
        idx_eta = (t_idx-1)*Q + qq
        V_eta(idx_eta,idx_eta) = V_eta(idx_eta,idx_eta) + Sigma_eps(qq,qq)
      END DO
    END DO

    ALLOCATE(L_V_cache(dimEta,dimEta))
    L_V_cache(:,:) = V_eta(:,:)
    CALL DPOTRF('L', dimEta, L_V_cache, dimEta, info)

    DEALLOCATE(tmp, V_eta)

  ELSE
    ! mode 1: NON serve L_V_cache
    ALLOCATE(L_V_cache(1,1))
    L_V_cache(1,1) = 0.d0
    info = 0
  END IF

END SUBROUTINE ensure_eta_cache







   !=========================================================
  ! loglik_cond_i:
  !  log-likelihood condizionata a un singolo vettore z (dimZ)
  !  per il soggetto i.
  !  z contiene i random standard:
  !     primi Q*p_ran = blocchi per i fattori
  !     ultimi nvar   = random sugli item (se re_item_flag=1)
  !=========================================================
    SUBROUTINE loglik_cond_i(i, sm, dimZ, z, loglik_i)
  USE shared
  USE utilities
  USE datindiv
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i, sm, dimZ
  DOUBLE PRECISION, DIMENSION(dimZ), INTENT(IN) :: z
  DOUBLE PRECISION, INTENT(OUT) :: loglik_i

  DOUBLE PRECISION, PARAMETER :: LIK_FLOOR = 1.0D-300
  INTEGER :: info
  INTEGER :: dimEta, dimLat, base_items
  INTEGER :: k_item, j_mes, idx_obs, kk, yijk, tuk, t_idx, qq, col, rowEta
  DOUBLE PRECISION :: mu_ijk, tausup, tauinf, p
  DOUBLE PRECISION :: invsd_lat_k

  DOUBLE PRECISION, ALLOCATABLE :: eta_vec(:), a(:)

  ! mode 0
  DOUBLE PRECISION, ALLOCATABLE :: z_lat_eta(:)

  ! mode 1
  DOUBLE PRECISION, ALLOCATABLE :: z_lat_b(:), b_vec(:)

  IF (sm < 0) THEN
    loglik_i = -1.d300
    RETURN
  END IF

  dimEta = Q*n
  IF (integ_mode == 0) THEN
    dimLat = dimEta
  ELSE
    dimLat = dimB
  END IF
  base_items = dimLat

  IF (re_item_flag == 1) THEN
    IF (dimZ /= dimLat + nvar) THEN
      loglik_i = -1.d300
      RETURN
    END IF
  ELSE
    IF (dimZ /= dimLat) THEN
      loglik_i = -1.d300
      RETURN
    END IF
  END IF

  CALL ensure_eta_cache(i, info)
  IF (info /= 0) THEN
    loglik_i = -1.d300
    RETURN
  END IF

  ALLOCATE(eta_vec(dimEta), a(nvar))

  ! =========================
  ! costruisci eta_vec
  ! =========================
  IF (integ_mode == 0) THEN
    ALLOCATE(z_lat_eta(dimEta))

    z_lat_eta(:) = z(1:dimEta)

    eta_vec(:) = m_eta_cache(:)
    DO rowEta = 1, dimEta
      DO col = 1, rowEta
        eta_vec(rowEta) = eta_vec(rowEta) + L_V_cache(rowEta,col) * z_lat_eta(col)
      END DO
    END DO

  ELSE
    ALLOCATE(z_lat_b(dimB), b_vec(dimB))

    z_lat_b(:) = z(1:dimB)

    ! b = L_b z_b
    b_vec(:) = 0.d0
    DO rowEta = 1, dimB
      DO col = 1, rowEta
        b_vec(rowEta) = b_vec(rowEta) + L_b(rowEta,col) * z_lat_b(col)
      END DO
    END DO

    ! eta = m + Zbig b
    eta_vec(:) = m_eta_cache(:)
    DO rowEta = 1, dimEta
      DO col = 1, dimB
        eta_vec(rowEta) = eta_vec(rowEta) + Zbig_cache(rowEta,col) * b_vec(col)
      END DO
    END DO
  END IF

  ! random item a_k
  IF (re_item_flag == 1) THEN
    DO k_item = 1, nvar
      a(k_item) = sqrt_sigmak(k_item) * z(base_items + k_item)
    END DO
  ELSE
    a(:) = 0.d0
  END IF

  ! log-lik condizionata
  loglik_i = 0.d0
  kk       = 0
  idx_obs  = nmescur

  DO k_item = 1, nvar

    IF (nmes(i,k_item) == 0) THEN
      kk = kk + (c(k_item) - 1)
      CYCLE
    END IF

    invsd_lat_k = invsd_k(k_item)

    DO j_mes = 1, nmes(i,k_item)
      idx_obs = idx_obs + 1
      yijk    = Y(idx_obs)
      t_idx   = nitem(idx_obs)

      mu_ijk = gamma(k_item) + a(k_item)
      DO qq = 1, Q
        mu_ijk = mu_ijk + lambda(k_item,qq) * eta_vec((t_idx-1)*Q + qq)
      END DO

      tuk = yijk + kk

      IF (yijk == 1) THEN
        tausup = (tau(tuk) - mu_ijk) * invsd_lat_k
        p = alnorm_fast(tausup, .FALSE.)
      ELSE IF (yijk == c(k_item)) THEN
        tauinf = (tau(tuk-1) - mu_ijk) * invsd_lat_k
        p = 1.d0 - alnorm_fast(tauinf, .FALSE.)
      ELSE
        tausup = (tau(tuk)   - mu_ijk) * invsd_lat_k
        tauinf = (tau(tuk-1) - mu_ijk) * invsd_lat_k
        p = alnorm_fast(tausup, .FALSE.) - alnorm_fast(tauinf, .FALSE.)
      END IF

      IF (p < LIK_FLOOR) p = LIK_FLOOR
      loglik_i = loglik_i + LOG(p)
    END DO

    kk = kk + (c(k_item) - 1)
  END DO

  IF (integ_mode == 0) THEN
    DEALLOCATE(z_lat_eta)
  ELSE
    DEALLOCATE(z_lat_b, b_vec)
  END IF

  DEALLOCATE(eta_vec, a)

END SUBROUTINE loglik_cond_i


  SUBROUTINE logpost_i(i, sm, dimZ, z, f)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, sm, dimZ
    DOUBLE PRECISION, DIMENSION(dimZ), INTENT(IN) :: z
    DOUBLE PRECISION, INTENT(OUT) :: f
    DOUBLE PRECISION :: loglik_i

    CALL loglik_cond_i(i, sm, dimZ, z, loglik_i)
    f = loglik_i - 0.5d0 * DOT_PRODUCT(z, z)
  END SUBROUTINE logpost_i




   SUBROUTINE grad_hess_logpost_i(i, sm, dimZ, z, f, grad, H)
  USE shared
  USE utilities
  USE datindiv
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i, sm, dimZ
  DOUBLE PRECISION, DIMENSION(dimZ), INTENT(IN)  :: z
  DOUBLE PRECISION, INTENT(OUT) :: f
  DOUBLE PRECISION, DIMENSION(dimZ), INTENT(OUT) :: grad
  DOUBLE PRECISION, DIMENSION(dimZ,dimZ), INTENT(OUT) :: H

 DOUBLE PRECISION, PARAMETER :: LIK_FLOOR = 1.0D-300
  DOUBLE PRECISION, PARAMETER :: INV_SQRT_2PI = 0.39894228040143267794d0

  INTEGER :: info, dimEta, dimLat, base_items
  INTEGER :: k_item, j_mes, idx_obs, kk, yijk, tuk, t_idx, qq, col, rowEta
  INTEGER :: d1, d2, d_item, m
  DOUBLE PRECISION :: mu_ijk, invsd_lat_k
  DOUBLE PRECISION :: alpha, beta_thr, p
  DOUBLE PRECISION :: phi_a, phi_b
  DOUBLE PRECISION :: dP_dmu, d2P_dmu2
  DOUBLE PRECISION :: s_mu, h_mu
  DOUBLE PRECISION :: loglik_i, qz
  DOUBLE PRECISION :: lam, tmp

  DOUBLE PRECISION, ALLOCATABLE :: eta_vec(:), a(:), dmu_dz(:)

  ! mode 0
  DOUBLE PRECISION, ALLOCATABLE :: z_lat_eta(:)

  ! mode 1
  DOUBLE PRECISION, ALLOCATABLE :: z_lat_b(:), b_vec(:)

  IF (sm < 0) THEN
    f = -1.d300
    grad = 0.d0
    H = 0.d0
    RETURN
  END IF

  dimEta = Q*n
  IF (integ_mode == 0) THEN
    dimLat = dimEta
  ELSE
    dimLat = dimB
  END IF
  base_items = dimLat

  IF (re_item_flag == 1) THEN
    IF (dimZ /= dimLat + nvar) THEN
      f = -1.d300; grad = 0.d0; H = 0.d0; RETURN
    END IF
  ELSE
    IF (dimZ /= dimLat) THEN
      f = -1.d300; grad = 0.d0; H = 0.d0; RETURN
    END IF
  END IF

  CALL ensure_eta_cache(i, info)
  IF (info /= 0) THEN
    f = -1.d300; grad = 0.d0; H = 0.d0; RETURN
  END IF

  ALLOCATE(eta_vec(dimEta), a(nvar), dmu_dz(dimZ))

  ! =========================
  ! eta_vec
  ! =========================
  IF (integ_mode == 0) THEN
    ALLOCATE(z_lat_eta(dimEta))
    z_lat_eta(:) = z(1:dimEta)

    eta_vec(:) = m_eta_cache(:)
    DO rowEta = 1, dimEta
      DO col = 1, rowEta
        eta_vec(rowEta) = eta_vec(rowEta) + L_V_cache(rowEta,col) * z_lat_eta(col)
      END DO
    END DO
  ELSE
    ALLOCATE(z_lat_b(dimB), b_vec(dimB))
    z_lat_b(:) = z(1:dimB)

    b_vec(:) = 0.d0
    DO rowEta = 1, dimB
      DO col = 1, rowEta
        b_vec(rowEta) = b_vec(rowEta) + L_b(rowEta,col) * z_lat_b(col)
      END DO
    END DO

    eta_vec(:) = m_eta_cache(:)
    DO rowEta = 1, dimEta
      DO col = 1, dimB
        eta_vec(rowEta) = eta_vec(rowEta) + Zbig_cache(rowEta,col) * b_vec(col)
      END DO
    END DO
  END IF

  ! a_k
  IF (re_item_flag == 1) THEN
    DO k_item = 1, nvar
      a(k_item) = sqrt_sigmak(k_item) * z(base_items + k_item)
    END DO
  ELSE
    a(:) = 0.d0
  END IF

  grad     = 0.d0
  H        = 0.d0
  loglik_i = 0.d0
  kk       = 0
  idx_obs  = nmescur

  DO k_item = 1, nvar

    IF (nmes(i,k_item) == 0) THEN
      kk = kk + (c(k_item) - 1)
      CYCLE
    END IF

    invsd_lat_k = invsd_k(k_item)

    DO j_mes = 1, nmes(i,k_item)
      idx_obs = idx_obs + 1
      yijk    = Y(idx_obs)
      t_idx   = nitem(idx_obs)

      ! mu
      mu_ijk = gamma(k_item) + a(k_item)
      DO qq = 1, Q
        mu_ijk = mu_ijk + lambda(k_item,qq) * eta_vec((t_idx-1)*Q + qq)
      END DO

      tuk = yijk + kk

      ! p e derivate wrt mu
      IF (yijk == 1) THEN
        alpha = (tau(tuk) - mu_ijk) * invsd_lat_k
        p = alnorm_fast(alpha, .FALSE.)
        IF (p < LIK_FLOOR) p = LIK_FLOOR
        phi_a    = INV_SQRT_2PI * EXP(-0.5d0*alpha*alpha)
        dP_dmu   = -phi_a * invsd_lat_k
        d2P_dmu2 = -alpha * phi_a * invsd_lat_k * invsd_lat_k
      ELSE IF (yijk == c(k_item)) THEN
        beta_thr = (tau(tuk-1) - mu_ijk) * invsd_lat_k
        p = alnorm_fast(beta_thr, .TRUE.)
        IF (p < LIK_FLOOR) p = LIK_FLOOR
        phi_b    = INV_SQRT_2PI * EXP(-0.5d0*beta_thr*beta_thr)
        dP_dmu   =  phi_b * invsd_lat_k
        d2P_dmu2 =  beta_thr * phi_b * invsd_lat_k * invsd_lat_k
      ELSE
        alpha    = (tau(tuk)   - mu_ijk) * invsd_lat_k
        beta_thr = (tau(tuk-1) - mu_ijk) * invsd_lat_k
        IF (beta_thr > 0.d0) THEN
                  p = alnorm_fast(beta_thr, .TRUE.) - alnorm_fast(alpha, .TRUE.)
               ELSE
                  p = alnorm_fast(alpha, .FALSE.) - alnorm_fast(beta_thr, .FALSE.)
               END IF
        IF (p < LIK_FLOOR) p = LIK_FLOOR
        phi_a = INV_SQRT_2PI * EXP(-0.5d0*alpha*alpha)
        phi_b = INV_SQRT_2PI * EXP(-0.5d0*beta_thr*beta_thr)
        dP_dmu   = (phi_b - phi_a) * invsd_lat_k
        d2P_dmu2 = (beta_thr*phi_b - alpha*phi_a) * invsd_lat_k * invsd_lat_k
      END IF

      s_mu = dP_dmu / p
      h_mu = d2P_dmu2 / p - s_mu*s_mu
      loglik_i = loglik_i + LOG(p)

      ! =========================
      ! dmu/dz
      ! =========================
      dmu_dz(:) = 0.d0

      IF (integ_mode == 0) THEN
        ! z = [z_eta, z_item]
        DO qq = 1, Q
          rowEta = (t_idx-1)*Q + qq
          lam = lambda(k_item,qq)
          IF (lam == 0.d0) CYCLE
          DO col = 1, rowEta
            dmu_dz(col) = dmu_dz(col) + lam * L_V_cache(rowEta,col)
          END DO
        END DO
      ELSE
        ! z = [z_b, z_item],  b = L_b z_b, eta = m + Zbig b
        ! d eta(rowEta)/ d z_b(col) = sum_m Zbig(rowEta,m) * L_b(m,col)
        DO qq = 1, Q
          rowEta = (t_idx-1)*Q + qq
          lam = lambda(k_item,qq)
          IF (lam == 0.d0) CYCLE
          DO col = 1, dimB
            tmp = 0.d0
            DO m = 1, dimB
              tmp = tmp + Zbig_cache(rowEta,m) * L_b(m,col)
            END DO
            dmu_dz(col) = dmu_dz(col) + lam * tmp
          END DO
        END DO
      END IF

      ! item random
      IF (re_item_flag == 1) THEN
        d_item = base_items + k_item
        dmu_dz(d_item) = dmu_dz(d_item) + sqrt_sigmak(k_item)
      END IF

      ! grad/H
      DO d1 = 1, dimZ
        IF (dmu_dz(d1) /= 0.d0) grad(d1) = grad(d1) + s_mu * dmu_dz(d1)
      END DO

      DO d1 = 1, dimZ
        IF (dmu_dz(d1) == 0.d0) CYCLE
        tmp = h_mu * dmu_dz(d1)
        DO d2 = 1, dimZ
          IF (dmu_dz(d2) == 0.d0) CYCLE
          H(d1,d2) = H(d1,d2) + tmp * dmu_dz(d2)
        END DO
      END DO

    END DO

    kk = kk + (c(k_item) - 1)
  END DO

  qz = DOT_PRODUCT(z, z)
  f  = loglik_i - 0.5d0*qz
  DO d1 = 1, dimZ
    grad(d1) = grad(d1) - z(d1)
    H(d1,d1) = H(d1,d1) - 1.d0
  END DO

  IF (integ_mode == 0) THEN
    DEALLOCATE(z_lat_eta)
  ELSE
    DEALLOCATE(z_lat_b, b_vec)
  END IF

  DEALLOCATE(eta_vec, a, dmu_dz)

END SUBROUTINE grad_hess_logpost_i



  SUBROUTINE num_grad_hess_logpost_i(i, sm, dimZ, z, grad, H)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, sm, dimZ
    DOUBLE PRECISION, DIMENSION(dimZ),      INTENT(IN)  :: z
    DOUBLE PRECISION, DIMENSION(dimZ),      INTENT(OUT) :: grad
    DOUBLE PRECISION, DIMENSION(dimZ,dimZ), INTENT(OUT) :: H

    INTEGER :: d, e
    DOUBLE PRECISION :: f0, fp, fm, fpp, fpm, fmp, fmm
    DOUBLE PRECISION :: eps
    DOUBLE PRECISION, DIMENSION(dimZ) :: zp, zm

    eps = 1.d-5

    CALL logpost_i(i, sm, dimZ, z, f0)

    ! grad
    DO d = 1, dimZ
      zp = z
      zm = z
      zp(d) = zp(d) + eps
      zm(d) = zm(d) - eps
      CALL logpost_i(i, sm, dimZ, zp, fp)
      CALL logpost_i(i, sm, dimZ, zm, fm)
      grad(d) = (fp - fm) / (2.d0*eps)
    END DO

    ! Hessiana (centrale)
    DO d = 1, dimZ
      DO e = d, dimZ
        zp = z; zm = z
        zp(d) = zp(d) + eps; zp(e) = zp(e) + eps
        CALL logpost_i(i, sm, dimZ, zp, fpp)

        zp = z
        zp(d) = zp(d) + eps; zp(e) = zp(e) - eps
        CALL logpost_i(i, sm, dimZ, zp, fpm)

        zp = z
        zp(d) = zp(d) - eps; zp(e) = zp(e) + eps
        CALL logpost_i(i, sm, dimZ, zp, fmp)

        zp = z
        zp(d) = zp(d) - eps; zp(e) = zp(e) - eps
        CALL logpost_i(i, sm, dimZ, zp, fmm)

        H(d,e) = (fpp - fpm - fmp + fmm) / (4.d0*eps*eps)
        H(e,d) = H(d,e)
      END DO
    END DO
  END SUBROUTINE num_grad_hess_logpost_i



  SUBROUTINE find_mode_and_R(i, sm, dimZ, z_mode, R, logdetR)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, sm, dimZ
    DOUBLE PRECISION, DIMENSION(dimZ),   INTENT(OUT) :: z_mode
    DOUBLE PRECISION, DIMENSION(dimZ,dimZ), INTENT(OUT) :: R
    DOUBLE PRECISION, INTENT(OUT) :: logdetR

    INTEGER, PARAMETER :: max_iter = 15
    DOUBLE PRECISION, PARAMETER :: tol_grad = 1.d-4
    DOUBLE PRECISION, PARAMETER :: tol_step = 1.d-4
    DOUBLE PRECISION, PARAMETER :: alpha_min = 1.d-4

    INTEGER :: iter, info, j, row, col
    DOUBLE PRECISION :: f0, f_trial, alpha
    DOUBLE PRECISION :: max_step, max_g

    DOUBLE PRECISION, DIMENSION(dimZ) :: z, grad, delta, z_trial, grad_trial
    DOUBLE PRECISION, DIMENSION(dimZ,dimZ) :: H, A, H_trial

    ! ------------------------
    ! start dal centro 0
    ! ------------------------
    z = 0.d0
    CALL grad_hess_logpost_i(i, sm, dimZ, z, f0, grad, H)

    ! Se qualcosa va storto già qui, fallback
    IF (f0 <= -1.d250) THEN
      z_mode  = z
      R       = 0.d0
      DO j = 1, dimZ
        R(j,j) = 1.d0
      END DO
      logdetR = 0.d0
      RETURN
    END IF

    ! ------------------------
    ! Newton ascent con line search
    ! ------------------------
    DO iter = 1, max_iter

      ! (opzionale ma robusto) forza simmetria di H
      DO row = 1, dimZ
        DO col = row+1, dimZ
          H(row,col) = 0.5d0*(H(row,col) + H(col,row))
          H(col,row) = H(row,col)
        END DO
      END DO

      ! A = (-H)^{-1}; nessuna regolarizzazione diagonale.
      A = -H
      CALL DPOTRF('L', dimZ, A, dimZ, info)

      IF (info /= 0) THEN
        ! fallback: niente adattamento
        z_mode  = z
        R       = 0.d0
        DO j = 1, dimZ
          R(j,j) = 1.d0
        END DO
        logdetR = 0.d0
        RETURN
      END IF

      CALL DPOTRI('L', dimZ, A, dimZ, info)
      IF (info /= 0) THEN
        z_mode  = z
        R       = 0.d0
        DO j = 1, dimZ
          R(j,j) = 1.d0
        END DO
        logdetR = 0.d0
        RETURN
      END IF

      ! DPOTRI('L') riempie la triangolare inferiore: simmetrizza A
      DO row = 1, dimZ
        DO col = row+1, dimZ
          A(row,col) = A(col,row)
        END DO
      END DO

      ! delta = A * grad
      CALL DGEMV('N', dimZ, dimZ, 1.d0, A, dimZ, grad, 1, 0.d0, delta, 1)

      ! line search
      alpha = 1.d0
      DO
        z_trial = z + alpha * delta
        CALL grad_hess_logpost_i(i, sm, dimZ, z_trial, f_trial, grad_trial, H_trial)

        IF (f_trial >= f0) EXIT
        alpha = alpha * 0.5d0
        IF (alpha < alpha_min) EXIT
      END DO

      IF (alpha < alpha_min) EXIT

      ! criteri stop
      max_step = MAXVAL(ABS(alpha*delta))
      max_g    = MAXVAL(ABS(grad_trial))

      z    = z_trial
      f0   = f_trial
      grad = grad_trial
      H    = H_trial

      IF (max_step < tol_step .AND. max_g < tol_grad) EXIT

    END DO

    z_mode = z

    ! ------------------------
    ! Sigma = (-H(mode))^{-1}, R = chol(Sigma)
    ! ------------------------
    CALL grad_hess_logpost_i(i, sm, dimZ, z_mode, f0, grad, H)

    ! simmetrizza H
    DO row = 1, dimZ
      DO col = row+1, dimZ
        H(row,col) = 0.5d0*(H(row,col) + H(col,row))
        H(col,row) = H(row,col)
      END DO
    END DO

    A = -H
    CALL DPOTRF('L', dimZ, A, dimZ, info)

    IF (info /= 0) THEN
      R       = 0.d0
      DO j = 1, dimZ
        R(j,j) = 1.d0
      END DO
      logdetR = 0.d0
      RETURN
    END IF

    CALL DPOTRI('L', dimZ, A, dimZ, info)
    IF (info /= 0) THEN
      R       = 0.d0
      DO j = 1, dimZ
        R(j,j) = 1.d0
      END DO
      logdetR = 0.d0
      RETURN
    END IF

    ! simmetrizza Sigma in A
    DO row = 1, dimZ
      DO col = row+1, dimZ
        A(row,col) = A(col,row)
      END DO
    END DO

    ! R = chol(Sigma)  (lower)
    R = A
    CALL DPOTRF('L', dimZ, R, dimZ, info)
    IF (info /= 0) THEN
      R       = 0.d0
      DO j = 1, dimZ
        R(j,j) = 1.d0
      END DO
      logdetR = 0.d0
      RETURN
    END IF

    ! log|R| = sum log diag(R)
    logdetR = 0.d0
    DO j = 1, dimZ
      logdetR = logdetR + LOG(R(j,j))
    END DO

    ! pulisci triangolare superiore
    DO col = 1, dimZ
      DO row = 1, col-1
        R(row,col) = 0.d0
      END DO
    END DO

  END SUBROUTINE find_mode_and_R

  !=============================================================
  !  Versione AGH: nodi/pesi specifici per soggetto
  !  nodes ha dimensione: nsub * nquad * dimZ
  !  weights ha dimensione: nsub * nquad
  !=============================================================
  DOUBLE PRECISION FUNCTION vrais_multo_i_agh(i, sm)
  USE shared
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i, sm

  INTEGER :: l, dimZ, dimEta, dimLat, d, baseZ
  DOUBLE PRECISION :: logdetR, loglik_cond, qz, qzeta
  DOUBLE PRECISION :: max_term, sum_exp, term

  DOUBLE PRECISION, ALLOCATABLE :: z_mode(:)
  DOUBLE PRECISION, ALLOCATABLE :: R(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: z(:), zeta(:)
  DOUBLE PRECISION, ALLOCATABLE :: term_log(:)

  dimEta = Q*n
  IF (integ_mode == 0) THEN
    dimLat = dimEta
  ELSE
    dimLat = dimB
  END IF

  dimZ = dimLat
  IF (re_item_flag == 1) dimZ = dimZ + nvar

  ALLOCATE(z_mode(dimZ), R(dimZ,dimZ), z(dimZ), zeta(dimZ))
  ALLOCATE(term_log(nquad))

  CALL find_mode_and_R(i, sm, dimZ, z_mode, R, logdetR)

  max_term = -1.d300
  DO l = 1, nquad
    baseZ = (l-1)*dimZ
    DO d = 1, dimZ
      zeta(d) = nodes(baseZ + d)
    END DO

    z = z_mode
    CALL DGEMV('N', dimZ, dimZ, 1.d0, R, dimZ, zeta, 1, 1.d0, z, 1)

    CALL loglik_cond_i(i, sm, dimZ, z, loglik_cond)

    qz    = DOT_PRODUCT(z, z)
    qzeta = DOT_PRODUCT(zeta, zeta)

    term = loglik_cond + (-0.5d0*(qz - qzeta)) + LOG(weights(l))
    term_log(l) = term
    IF (term > max_term) max_term = term
  END DO

  sum_exp = 0.d0
  DO l = 1, nquad
    sum_exp = sum_exp + EXP(term_log(l) - max_term)
  END DO

  IF (sum_exp <= 0.d0) THEN
    vrais_multo_i_agh = -1.d9
  ELSE
    vrais_multo_i_agh = logdetR + max_term + LOG(sum_exp)
  END IF

  DEALLOCATE(z_mode, R, z, zeta, term_log)

END FUNCTION vrais_multo_i_agh


  DOUBLE PRECISION FUNCTION vrais_multo_agh()
    IMPLICIT NONE
    INTEGER :: i, sm
    DOUBLE PRECISION :: temp

    nmescur        = 0
    vrais_multo_agh = 0.d0

    DO i = 1, nsub
       sm   = SUM(nmes(i,:))
       temp = vrais_multo_i_agh(i, sm)

       IF (temp == -1.d9) THEN
          vrais_multo_agh = -1.d9
          RETURN
       END IF

       vrais_multo_agh = vrais_multo_agh + temp
       nmescur = nmescur + sm
    END DO

  END FUNCTION vrais_multo_agh


  !************************************************************
  !  loglik_agh:
  !    stessa parametrizzazione di loglik, ma integra usando
  !    nodi/pesi soggetto-specifici (Adaptive GH à la Pinheiro)
  !************************************************************
SUBROUTINE loglik_agh(b0, Y0, nvar0, nsub0, nobs0, Q0, n0, &
                      nmes0, c0, cm10, ntre0, nquad0, &
                      nitem0, nodes0, weights0, npm0, &
                      p_fix0, p_ran0, X0, Z0, fac_item0, &
                      ident_type0, re_item0, integ_mode0, vrais) &
     BIND(C, name="loglik_agh")

  USE shared
  USE utilities
  IMPLICIT NONE

  ! ---- input ----
  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, Q0, n0
  INTEGER, INTENT(IN) :: ntre0, nquad0, npm0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0
  INTEGER, INTENT(IN) :: ident_type0, re_item0
  INTEGER, INTENT(IN) :: integ_mode0

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0),       INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: fac_item0

  DOUBLE PRECISION, DIMENSION(*),      INTENT(IN) :: nodes0
  DOUBLE PRECISION, DIMENSION(nquad0), INTENT(IN) :: weights0
  DOUBLE PRECISION, DIMENSION(npm0),   INTENT(IN) :: b0

  DOUBLE PRECISION, INTENT(IN), TARGET :: X0(*)
  DOUBLE PRECISION, INTENT(IN), TARGET :: Z0(*)

  ! ---- output ----
  DOUBLE PRECISION, INTENT(OUT) :: vrais

  ! ---- locals ----
  INTEGER :: k, q_idx, qq, rr, row, col, jj
  INTEGER :: dimEta, dimLat, dimZ, nX, nZ
  INTEGER :: npar_eps, npar_Lb
  INTEGER :: ntheta_here, nload, ngamma
  INTEGER :: pos, pos_tau, expected_npm

  INTEGER, ALLOCATABLE :: cm3(:)
  DOUBLE PRECISION, ALLOCATABLE :: ta(:)

  INTEGER, ALLOCATABLE :: anchor(:), Jq(:)
  LOGICAL, ALLOCATABLE :: is_anchor(:)

  ! ------------------ set global dims ------------------
  nvar  = nvar0
  nsub  = nsub0
  nobs  = nobs0
  Q     = Q0
  n     = n0
  ntre  = ntre0
  nquad = nquad0
  npm   = npm0
  p_fix = p_fix0
  p_ran = p_ran0

  dimEta = Q*n
  dimB   = Q*p_ran

  re_item_flag = re_item0
  integ_mode   = integ_mode0

  ! dimLat e dimZ
  IF (integ_mode == 0) THEN
    dimLat = dimEta
  ELSE
    dimLat = dimB
  END IF

  dimZ = dimLat
  IF (re_item0 == 1) dimZ = dimZ + nvar

  ! eps variances: SOLO in mode 0 e SOLO se ident_type != unit_variance
  IF (integ_mode == 1) THEN
    npar_eps = 0
  ELSE
    IF (ident_type0 == 1) THEN
      npar_eps = 0
    ELSE
      npar_eps = Q
    END IF
  END IF

  npar_Lb = dimB*(dimB+1)/2

  ! ------------------ anchors & Jq ------------------
  ALLOCATE(anchor(Q), Jq(Q), is_anchor(nvar))
  anchor = 0
  Jq     = 0
  is_anchor = .FALSE.

  DO k = 1, nvar
    q_idx = fac_item0(k)
    Jq(q_idx) = Jq(q_idx) + 1
    IF (anchor(q_idx) == 0) anchor(q_idx) = k
  END DO
  DO q_idx = 1, Q
    is_anchor(anchor(q_idx)) = .TRUE.
  END DO

  ! ------------------ reset cache (AGH) ------------------
  cache_i_eta   = -1
  cache_dimEta  = -1
  IF (ALLOCATED(m_eta_cache)) DEALLOCATE(m_eta_cache)
  IF (ALLOCATED(L_V_cache))   DEALLOCATE(L_V_cache)
  IF (ALLOCATED(Zbig_cache))  DEALLOCATE(Zbig_cache)

  ! ------------------ data ------------------
  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs), c(nvar), cm1(nvar))
  nmes  = nmes0
  Y     = Y0
  nitem = nitem0
  c     = c0
  cm1   = cm10

  ! theta stimata solo per non-binary (cm1>1)
  ntheta_here = 0
  DO k = 1, nvar
    IF (cm1(k) > 1) ntheta_here = ntheta_here + 1
  END DO

  ! ------------------ nodes/weights ------------------
  ALLOCATE(nodes(nquad*dimZ))
  nodes(1:nquad*dimZ) = nodes0(1:nquad*dimZ)

  ALLOCATE(weights(nquad))
  weights(1:nquad) = weights0(1:nquad)

  ! ------------------ allocate global params ------------------
  ALLOCATE(lambda(nvar,Q))
  ALLOCATE(gamma(nvar))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(beta(Q,p_fix))
  ALLOCATE(L_b(dimB,dimB), Sigma_b(dimB,dimB))
  ALLOCATE(Sigma_eps(Q,Q))
  ALLOCATE(sigmak(nvar), sqrt_sigmak(nvar))
  ALLOCATE(X_eta(nsub,n,Q,p_fix))
  ALLOCATE(Z_eta(nsub,n,Q,p_ran))

  nX = nsub*n*Q*p_fix
  nZ = nsub*n*Q*p_ran
  X_eta = RESHAPE( X0(1:nX), (/ nsub, n, Q, p_fix /) )
  Z_eta = RESHAPE( Z0(1:nZ), (/ nsub, n, Q, p_ran /) )

  ! ------------------ thresholds ------------------
  ALLOCATE(cm3(nvar))
  DO k = 1, nvar
    cm3(k) = cm1(k) - 2
    IF (cm3(k) < 0) cm3(k) = 0
  END DO
  ntre = SUM(cm3)

  IF (ntre > 0) THEN
    ALLOCATE(ta(ntre))
    ta = b0(1:ntre)
  ELSE
    ALLOCATE(ta(1))
    ta = 0.d0
  END IF

  ALLOCATE(tau(SUM(cm1)))

  pos_tau = 0
  jj = 1
  DO k = 1, nvar
    IF (cm1(k) == 1) THEN
      tau(pos_tau+1) = 0.d0
    ELSE
      tau(pos_tau+1) = 0.d0
      tau(pos_tau+2) = 1.d0
      IF (cm1(k) > 2) THEN
        tau(pos_tau+3 : pos_tau+cm1(k)) = ta(jj : jj+cm3(k)-1)
        jj = jj + cm3(k)
      END IF
    END IF
    pos_tau = pos_tau + cm1(k)
  END DO

  pos = ntre

  ! ------------------ loadings ------------------
  lambda = 0.d0

  IF (ident_type0 == 1) THEN
    ! unit_variance: tutte le lambda stimate
    DO k = 1, nvar
      q_idx = fac_item0(k)
      pos = pos + 1
      lambda(k,q_idx) = b0(pos)
    END DO

  ELSE IF (ident_type0 == 0) THEN
    ! anchor_loadings
    DO q_idx = 1, Q
      lambda(anchor(q_idx),q_idx) = 1.d0
    END DO
    DO q_idx = 1, Q
      DO k = 1, nvar
        IF (fac_item0(k) == q_idx .AND. k /= anchor(q_idx)) THEN
          pos = pos + 1
          lambda(k,q_idx) = b0(pos)
        END IF
      END DO
    END DO

  ELSE IF (ident_type0 == 2 .OR. ident_type0 == 3) THEN
    ! ECI2 / ECI1
    DO q_idx = 1, Q
      DO k = 1, nvar
        IF (fac_item0(k) == q_idx .AND. k /= anchor(q_idx)) THEN
          pos = pos + 1
          lambda(k,q_idx) = b0(pos)
        END IF
      END DO

      IF (ident_type0 == 2) THEN
        lambda(anchor(q_idx),q_idx) = 1.d0
      ELSE
        lambda(anchor(q_idx),q_idx) = DBLE(Jq(q_idx))  ! EC-1: somma=Jq -> media=1
      END IF

      DO k = 1, nvar
        IF (fac_item0(k) == q_idx .AND. k /= anchor(q_idx)) THEN
          lambda(anchor(q_idx),q_idx) = lambda(anchor(q_idx),q_idx) - lambda(k,q_idx)
        END IF
      END DO
    END DO
  END IF

  ! ------------------ gamma ------------------
  gamma = 0.d0

  IF (ident_type0 == 2 .OR. ident_type0 == 3) THEN
    DO q_idx = 1, Q
      DO k = 1, nvar
        IF (fac_item0(k) == q_idx .AND. k /= anchor(q_idx)) THEN
          pos = pos + 1
          gamma(k) = b0(pos)
        END IF
      END DO

      gamma(anchor(q_idx)) = 0.d0
      DO k = 1, nvar
        IF (fac_item0(k) == q_idx .AND. k /= anchor(q_idx)) THEN
          gamma(anchor(q_idx)) = gamma(anchor(q_idx)) - gamma(k)
        END IF
      END DO
    END DO
  ELSE
    DO k = 1, nvar
      pos = pos + 1
      gamma(k) = b0(pos)
    END DO
  END IF

  ! ------------------ theta ------------------
  DO k = 1, nvar
    IF (cm1(k) == 1) THEN
      theta_k(k) = 1.d0
      invsd_k(k) = 1.d0
    ELSE
      pos = pos + 1
      theta_k(k) = EXP(b0(pos))
      invsd_k(k) = 1.d0 / SQRT(theta_k(k))
    END IF
  END DO

  ! ------------------ beta ------------------
  DO qq = 1, Q
    DO rr = 1, p_fix
      pos = pos + 1
      beta(qq,rr) = b0(pos)
    END DO
  END DO

  ! ------------------ L_b ------------------
  L_b(:,:) = 0.d0
  DO rr = 1, dimB
    pos = pos + 1
    L_b(rr,rr) = EXP(b0(pos))
  END DO
  DO col = 1, dimB
    DO row = col+1, dimB
      pos = pos + 1
      L_b(row,col) = b0(pos)
    END DO
  END DO

  ! ------------------ Sigma_eps ------------------
  Sigma_eps(:,:) = 0.d0
  IF (integ_mode == 0) THEN
    IF (ident_type0 == 1) THEN
      DO qq = 1, Q
        Sigma_eps(qq,qq) = 1.d0
      END DO
    ELSE
      DO qq = 1, Q
        pos = pos + 1
        Sigma_eps(qq,qq) = EXP(b0(pos))
      END DO
    END IF
  ELSE
    ! mode 1: epsilon=0 => nessun parametro letto
    Sigma_eps(:,:) = 0.d0
  END IF

  ! ------------------ Sigma_b = L_b L_b' ------------------
  Sigma_b(:,:) = 0.d0
  CALL DGEMM('N','T', dimB, dimB, dimB, 1.d0, L_b, dimB, L_b, dimB, 0.d0, Sigma_b, dimB)

  ! ------------------ sigmak ------------------
  IF (re_item0 == 1) THEN
    DO k = 1, nvar
      pos = pos + 1
      sigmak(k)      = EXP(b0(pos))
      sqrt_sigmak(k) = SQRT(sigmak(k))
    END DO
  ELSE
    sigmak(:)      = 0.d0
    sqrt_sigmak(:) = 0.d0
  END IF

  ! ------------------ (opzionale) check conteggi ------------------
  IF (ident_type0 == 1) THEN
    nload  = nvar
    ngamma = nvar
  ELSE IF (ident_type0 == 0) THEN
    nload  = nvar - Q
    ngamma = nvar
  ELSE
    nload  = nvar - Q
    ngamma = nvar - Q
  END IF

  IF (re_item0 == 1) THEN
    expected_npm = ntre + nload + ngamma + ntheta_here + Q*p_fix + npar_Lb + npar_eps + nvar
  ELSE
    expected_npm = ntre + nload + ngamma + ntheta_here + Q*p_fix + npar_Lb + npar_eps
  END IF
  ! se vuoi: IF (expected_npm /= npm0) THEN ... END IF

  ! ------------------ likelihood (AGH) ------------------
  vrais = -vrais_multo_agh()

  ! ------------------ cleanup ------------------
  DEALLOCATE(nmes, Y, nitem, c, cm1)
  DEALLOCATE(nodes, weights)
  DEALLOCATE(lambda, gamma, theta_k, invsd_k, beta)
  DEALLOCATE(L_b, Sigma_b, Sigma_eps)
  DEALLOCATE(sigmak, sqrt_sigmak)
  DEALLOCATE(tau, cm3, ta)
  DEALLOCATE(X_eta, Z_eta)
  DEALLOCATE(anchor, Jq, is_anchor)

  cache_i_eta   = -1
  cache_dimEta  = -1
  IF (ALLOCATED(m_eta_cache)) DEALLOCATE(m_eta_cache)
  IF (ALLOCATED(L_V_cache))   DEALLOCATE(L_V_cache)
  IF (ALLOCATED(Zbig_cache))  DEALLOCATE(Zbig_cache)

END SUBROUTINE loglik_agh






  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!                          !!!!
  !!!                            !!!
  !!      CFA longitudinale       !!
  !!!                            !!!
  !!!!                          !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !===========================================================
!  CFA LONGITUDINALE MULTIDIMENSIONALE + (opz.) RE item
!  eta_{t,q}, t=1..n, q=1..Q  -> vettore eta_vec( (t-1)*Q + q )
!  eta_vec = kappa + L_phi * z_eta,  z_eta ~ N(0, I_{nQ})
!  a_item(k) = sqrt_sigmak(k) * z_item(k)  se re_item_flag=1
!===========================================================
DOUBLE PRECISION FUNCTION vrais_long_i(i, obs0)
   USE shared, ONLY: nvar, n, Q, nquad, nodes, logweights, re_item_flag, &
                  nmes, c, Y, nitem, invsd_long, gamma, lambda, &
                  sqrt_sigmak, tau, kappa, L_phi, &
                  phi_struct_flag, A_ar, sqrt_sig1, &
                  kk_base, nload_item, load_idx, lambda_act
  USE utilities, ONLY: alnorm_fast, eta_from_z_ar1
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, obs0
  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9

  INTEGER :: l, k_item, j_mes, idx_obs, yijk, tuk
  INTEGER :: t_idx, dimEta, dimZ, baseZ
  INTEGER :: idim, irow, icol, qq, idx_tq, ll
  INTEGER :: kk0, nmes_k, c_k, base_eta

  DOUBLE PRECISION :: logp
  DOUBLE PRECISION :: mu_ijk, mu0_k, linpred
  DOUBLE PRECISION :: invsd_lat_k
  DOUBLE PRECISION :: tausup, tauinf, p, z_item
  DOUBLE PRECISION :: loglik_l, logterm
  DOUBLE PRECISION :: log_max, sumexp

  DOUBLE PRECISION :: eta_lat(n*Q), mean_eta(n*Q), z_lat(n*Q), a_re(nvar)


  dimEta  = n * Q

  IF (re_item_flag == 1) THEN
     dimZ = dimEta + nvar
  ELSE
     dimZ = dimEta
  END IF

  IF (phi_struct_flag == 2) THEN
     IF (.NOT. ALLOCATED(A_ar) .OR. .NOT. ALLOCATED(sqrt_sig1)) THEN
        CALL rexit("ERROR: A_ar/sqrt_sig1 not allocated")
     END IF
  END IF

  CALL cfa_eta_mean(i, mean_eta)

  ! som = 0.d0
  log_max = -1.d300
  sumexp  = 0.d0


  DO l = 1, nquad
     baseZ = (l-1) * dimZ


     ! z_lat ~ N(0, I)
     DO idim = 1, dimEta
        z_lat(idim) = nodes(baseZ + idim)
     END DO

     ! eta_lat
     IF (phi_struct_flag == 2) THEN
        CALL eta_from_z_ar1(n, Q, mean_eta, A_ar, sqrt_sig1, z_lat, eta_lat)
     ELSE
        DO irow = 1, dimEta
           eta_lat(irow) = mean_eta(irow)
           DO icol = 1, irow
              eta_lat(irow) = eta_lat(irow) + L_phi(irow, icol) * z_lat(icol)
           END DO
        END DO
     END IF

     ! random item effects a_ik
     IF (re_item_flag == 1) THEN
        DO k_item = 1, nvar
           z_item      = nodes(baseZ + dimEta + k_item)
           a_re(k_item) = sqrt_sigmak(k_item) * z_item
        END DO
     ELSE
        a_re = 0.d0
     END IF

     ! Likelihood condizionata
      ! lik_l   = 1.d0
     loglik_l = 0.d0
     idx_obs  = obs0 - 1

     DO k_item = 1, nvar

        nmes_k = nmes(i, k_item)
        IF (nmes_k == 0) CYCLE

        kk0   = kk_base(k_item)
        c_k   = c(k_item)
        mu0_k = gamma(k_item) + a_re(k_item)


        DO j_mes = 1, nmes_k
           idx_obs = idx_obs + 1
           t_idx   = nitem(idx_obs)

           invsd_lat_k = invsd_long(k_item, t_idx)
           base_eta    = (t_idx - 1) * Q

           IF (Q == 1) THEN
            linpred = lambda(k_item,1) * eta_lat(base_eta + 1)
          ELSE
            linpred = 0.d0
            DO ll = 1, nload_item(k_item)
                idx_tq = base_eta + load_idx(ll, k_item)
                linpred = linpred + lambda_act(ll, k_item) * eta_lat(idx_tq)
            END DO
          END IF

           mu_ijk = mu0_k + linpred
           yijk   = Y(idx_obs)
           tuk    = yijk + kk0

           IF (c_k == 2) THEN
            tausup = - mu_ijk * invsd_lat_k
            IF (yijk == 1) THEN
                p = alnorm_fast(tausup, .FALSE.)
            ELSE
                p = 1.d0 - alnorm_fast(tausup, .FALSE.)
            END IF
          ELSE
              IF (yijk == 1) THEN
                 tausup = (tau(tuk) - mu_ijk) * invsd_lat_k
                 p = alnorm_fast(tausup, .FALSE.)
              ELSE IF (yijk == c_k) THEN
                 tauinf = (tau(tuk-1) - mu_ijk) * invsd_lat_k
                 p = 1.d0 - alnorm_fast(tauinf, .FALSE.)
              ELSE
                 tausup = (tau(tuk)   - mu_ijk) * invsd_lat_k
                 tauinf = (tau(tuk-1) - mu_ijk) * invsd_lat_k
                 p = alnorm_fast(tausup, .FALSE.) - alnorm_fast(tauinf, .FALSE.)
              END IF
           END IF

            !  IF (p < LIK_FLOOR) p = LIK_FLOOR
          !  lik_l = lik_l * p

          IF (p > 0.d0) THEN
            logp = LOG(p)
            IF (logp == logp) THEN
                loglik_l = loglik_l + logp
            END IF

          ELSE IF (p == 0.d0) THEN
            loglik_l = LOG_ZERO

          ELSE
            ! p negativo o NaN ignorato
          END IF
            !  ---- DIAGNOSTICA p ----
            ! if (p /= p .or. p <= 0.d0 .or. p > 1.d0) then
            !    write(*,*) "BAD p in vrais_long_i"
            !    write(*,*) "i=", i, " l(node)=", l, " k_item=", k_item, " j_mes=", j_mes, " t_idx=", t_idx
            !    write(*,*) "y=", yijk, " c=", c(k_item), " kk_thr=", kk_thr, " tuk=", tuk
            !    write(*,*) "mu_ijk=", mu_ijk, " invsd_lat_k=", invsd_lat_k
            !    write(*,*) "tausup=", tausup, " tauinf=", tauinf, " p=", p

            !    ! stampa le soglie dell'item (come sono in tau)
            !    write(*,*) "tau for this item:"
            !    write(*,'(*(ES24.16,1X))') ( tau(kk_thr + q), q=1, c(k_item)-1 )

            !    stop
            ! end if

            ! -----------------------
        END DO


     END DO

      ! som = som + weights(l) * lik_l
      logterm = logweights(l) + loglik_l

      IF (logterm == logterm) THEN
        IF (sumexp == 0.d0) THEN
            log_max = logterm
            sumexp  = 1.d0
        ELSE IF (logterm > log_max) THEN
            sumexp  = sumexp * EXP(log_max - logterm) + 1.d0
            log_max = logterm
        ELSE IF (logterm == log_max) THEN
            sumexp = sumexp + 1.d0
        ELSE
            sumexp = sumexp + EXP(logterm - log_max)
        END IF
      END IF

  END DO

  ! IF (som <= 0.d0) THEN
  !    vrais_long_i = -1.d9
  ! ELSE
  !    vrais_long_i = LOG(som)
  ! END IF
  IF (sumexp <= 0.d0) THEN
   vrais_long_i = LOG_ZERO
ELSE
   vrais_long_i = log_max + LOG(sumexp)
END IF


END FUNCTION vrais_long_i






DOUBLE PRECISION FUNCTION vrais_long()
  USE shared, ONLY: nsub, nmes
  IMPLICIT NONE

  INTEGER :: i, sm, obs0
  DOUBLE PRECISION :: temp

  obs0       = 1
  vrais_long = 0.d0

  DO i = 1, nsub
     sm   = SUM(nmes(i,:))
     temp = vrais_long_i(i, obs0)

   vrais_long = vrais_long + temp

     obs0 = obs0 + sm

  END DO
END FUNCTION vrais_long

!===========================================================
! loglik_long: multidim con loadmask (BIND(C))
! Ordine parametri b0:
!  1) soglie libere (ntre0)
!  2) lambda liberi (solo loadmask==1, esclusi Q riferimenti fissati a 1)
!  3) gamma (nvar)
!  4) theta (solo ordinali; per tempo se theta_time_flag=1)
!  5) sigmak (se re_item_flag=1)
!  6) kappa(t>=2,q) -> (n-1)*Q  (kappa(1,q)=0 fissati)
!  7) L_phi su dim_eta=nQ -> dim_eta*(dim_eta+1)/2  (diag in log)
!===========================================================
SUBROUTINE loglik_long(b0, Y0, nvar0, nsub0, nobs0, n0,      &
                       nmes0, c0, cm10, ntre0,               &
                       nquad0, nitem0,                       &
                       nodes0, weights0,                     &
                       npm0, re_item0, theta_time0,          &
                       Q0, loadmask0, phi_struct0,           &
                       monotone_tau0, fit_id0,               &
                       vrais)                                &
        BIND(C, name="loglik_long")

  USE shared
  USE utilities
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: fit_id0
  INTEGER, INTENT(IN) :: monotone_tau0
  INTEGER, INTENT(IN) :: phi_struct0
  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0
  INTEGER, INTENT(IN) :: ntre0, nquad0, npm0
  INTEGER, INTENT(IN) :: re_item0, theta_time0
  INTEGER, INTENT(IN) :: Q0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0),       INTENT(IN) :: Y0, nitem0

  DOUBLE PRECISION, DIMENSION(*),      INTENT(IN) :: nodes0
  DOUBLE PRECISION, DIMENSION(nquad0), INTENT(IN) :: weights0
  DOUBLE PRECISION, DIMENSION(npm0),   INTENT(IN) :: b0

  DOUBLE PRECISION, INTENT(OUT) :: vrais


  INTEGER :: j, k, t, row, col, qq, ll, r, maxload
  INTEGER :: pos, expected_npm
  INTEGER :: dimZ, dimEta
  INTEGER :: n_theta_est, n_lambda_free, ref_item
  INTEGER :: ntre_calc

  LOGICAL :: need_data_init, need_param_alloc

      ! monotonia lcmm-style: incrementi positivi via exp
      DOUBLE PRECISION :: prev
      INTEGER :: h, nf


    ! =========================
  ! DEBUG I/O (R -> Fortran)
  ! =========================
  INTEGER :: jj_end, nn_end, ww_end
  DOUBLE PRECISION :: sw



  INTEGER :: cm3(nvar0)
  DOUBLE PRECISION :: ta(MAX(1, ntre0))
  DOUBLE PRECISION :: A_tmp(Q0,Q0), sqrt1_tmp(Q0)
  INTEGER :: n_phi_par
  INTEGER :: jj, jj2, rr, cc
  !-----------------------------
  ! 1) global (shared)
  !-----------------------------
  nvar  = nvar0
  nsub  = nsub0
  nobs  = nobs0
  n     = n0
  nquad = nquad0
  npm   = npm0

  phi_struct_flag = phi_struct0

  re_item_flag    = re_item0
  theta_time_flag = theta_time0
  Q               = Q0

  p_fix = cfa_p_fix
  p_ran = 0

  dimEta = n * Q
  dim_eta = dimEta

  need_data_init = (.NOT. cache_ready) .OR. &
                 (fit_id0 /= cache_fit_id) .OR. &
                 (nvar0 /= cache_nvar) .OR. &
                 (nsub0 /= cache_nsub) .OR. &
                 (nobs0 /= cache_nobs) .OR. &
                 (n0   /= cache_n) .OR. &
                 (Q0   /= cache_Q) .OR. &
                 (nquad0 /= cache_nquad)


  !-----------------------------
  ! 2) dati (shared) - pulizia
  !-----------------------------
  IF (need_data_init) THEN
    IF (ALLOCATED(nmes))    DEALLOCATE(nmes)
    IF (ALLOCATED(Y))       DEALLOCATE(Y)
    IF (ALLOCATED(nitem))   DEALLOCATE(nitem)
    IF (ALLOCATED(c))       DEALLOCATE(c)
    IF (ALLOCATED(cm1))     DEALLOCATE(cm1)

    ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs), c(nvar), cm1(nvar))
    nmes  = nmes0
    Y     = Y0
    nitem = nitem0
    c     = c0
    cm1   = cm10

    IF (ALLOCATED(kk_base)) DEALLOCATE(kk_base)
    ALLOCATE(kk_base(nvar))

    kk_base(1) = 0
    DO k = 2, nvar
      kk_base(k) = kk_base(k-1) + (c(k-1) - 1)
    END DO

    ! loadmask globale
    IF (ALLOCATED(loadmask)) DEALLOCATE(loadmask)
    ALLOCATE(loadmask(nvar,Q))
    loadmask = loadmask0

    IF (ALLOCATED(nload_item)) DEALLOCATE(nload_item)
    ALLOCATE(nload_item(nvar))

    DO k = 1, nvar
      nload_item(k) = COUNT(loadmask(k,:) == 1)
    END DO

    maxload = MAXVAL(nload_item)
    IF (maxload < 1) maxload = 1

    IF (ALLOCATED(load_idx)) DEALLOCATE(load_idx)
    ALLOCATE(load_idx(maxload, nvar))
    load_idx = 0

    DO k = 1, nvar
      ll = 0
      DO qq = 1, Q
          IF (loadmask(k,qq) == 1) THEN
            ll = ll + 1
            load_idx(ll, k) = qq
          END IF
      END DO
    END DO

    ! dimensione integrazione
    IF (re_item_flag == 1) THEN
      dimZ = dimEta + nvar
    ELSE
      dimZ = dimEta
    END IF

    IF (ALLOCATED(nodes))      DEALLOCATE(nodes)
    IF (ALLOCATED(weights))    DEALLOCATE(weights)
    IF (ALLOCATED(logweights)) DEALLOCATE(logweights)

    ALLOCATE(nodes(nquad*dimZ))
    nodes(1:nquad*dimZ) = nodes0(1:nquad*dimZ)

    ALLOCATE(weights(nquad))
    weights(1:nquad) = weights0(1:nquad)

    ALLOCATE(logweights(nquad))
    DO ll = 1, nquad
      logweights(ll) = LOG(weights(ll))
    END DO

    cache_fit_id = fit_id0
    cache_nvar   = nvar0
    cache_nsub   = nsub0
    cache_nobs   = nobs0
    cache_n      = n0
    cache_Q      = Q0
    cache_nquad  = nquad0
    cache_ready  = .TRUE.
  END IF


  !-----------------------------
  ! 3) parametri (shared) - pulizia
  !-----------------------------
need_param_alloc = (.NOT. ALLOCATED(lambda)) .OR. &
                   (SIZE(lambda,1) /= nvar) .OR. (SIZE(lambda,2) /= Q) .OR. &
                   (.NOT. ALLOCATED(gamma)) .OR. (SIZE(gamma) /= nvar) .OR. &
                   (.NOT. ALLOCATED(theta_k)) .OR. (SIZE(theta_k) /= nvar) .OR. &
                   (.NOT. ALLOCATED(invsd_k)) .OR. (SIZE(invsd_k) /= nvar) .OR. &
                   (.NOT. ALLOCATED(sigmak)) .OR. (SIZE(sigmak) /= nvar) .OR. &
                   (.NOT. ALLOCATED(sqrt_sigmak)) .OR. (SIZE(sqrt_sigmak) /= nvar) .OR. &
                   (.NOT. ALLOCATED(kappa)) .OR. (SIZE(kappa) /= dimEta) .OR. &
                   (.NOT. ALLOCATED(L_phi)) .OR. (SIZE(L_phi,1) /= dimEta) .OR. (SIZE(L_phi,2) /= dimEta) .OR. &
                   (.NOT. ALLOCATED(theta_long)) .OR. (SIZE(theta_long,1) /= nvar) .OR. (SIZE(theta_long,2) /= n) .OR. &
                   (.NOT. ALLOCATED(invsd_long)) .OR. (SIZE(invsd_long,1) /= nvar) .OR. (SIZE(invsd_long,2) /= n) .OR. &
                   (.NOT. ALLOCATED(tau)) .OR. (SIZE(tau) /= SUM(cm1))

IF (need_param_alloc) THEN
   IF (ALLOCATED(lambda))       DEALLOCATE(lambda)
   IF (ALLOCATED(gamma))        DEALLOCATE(gamma)
   IF (ALLOCATED(theta_k))      DEALLOCATE(theta_k)
   IF (ALLOCATED(invsd_k))      DEALLOCATE(invsd_k)
   IF (ALLOCATED(sigmak))       DEALLOCATE(sigmak)
   IF (ALLOCATED(sqrt_sigmak))  DEALLOCATE(sqrt_sigmak)
   IF (ALLOCATED(kappa))        DEALLOCATE(kappa)
   IF (ALLOCATED(L_phi))        DEALLOCATE(L_phi)
   IF (ALLOCATED(tau))          DEALLOCATE(tau)
   IF (ALLOCATED(theta_long))   DEALLOCATE(theta_long)
   IF (ALLOCATED(invsd_long))   DEALLOCATE(invsd_long)

   ALLOCATE(lambda(nvar,Q))
   ALLOCATE(gamma(nvar))
   ALLOCATE(theta_k(nvar), invsd_k(nvar))
   ALLOCATE(sigmak(nvar),  sqrt_sigmak(nvar))
   ALLOCATE(kappa(dimEta))
   ALLOCATE(L_phi(dimEta,dimEta))
   ALLOCATE(tau(SUM(cm1)))
   ALLOCATE(theta_long(nvar,n))
   ALLOCATE(invsd_long(nvar,n))
END IF


    IF (Q > 1) THEN
      maxload = MAXVAL(nload_item)
      IF (maxload < 1) maxload = 1

      IF ((.NOT. ALLOCATED(lambda_act)) .OR. &
          (SIZE(lambda_act,1) /= maxload) .OR. (SIZE(lambda_act,2) /= nvar)) THEN
          IF (ALLOCATED(lambda_act)) DEALLOCATE(lambda_act)
          ALLOCATE(lambda_act(maxload, nvar))
      END IF
    END IF
  ! thresholds: identification controlled by cfa_threshold_id_flag
  DO k = 1, nvar
     cm3(k) = cm1(k) - 2
     IF (cm1(k) < 3) THEN
        CALL rexit('ERROR (CFA): threshold identification requires at least 4 response categories for every item')
     END IF
  END DO

  ntre_calc = SUM(cm3)
  IF (ntre_calc /= ntre0) THEN
     CALL rexit('ERROR (CFA): ntre0 mismatch. ntre0')
  END IF
  ntre = ntre0
  IF (ntre0 > 0) ta(1:ntre0) = b0(1:ntre0)

  pos = 0
  j = 1
  tau = 0.d0
  DO k = 1, nvar
     CALL cfa_fill_tau_item(cm1(k), ta(j), monotone_tau0, tau(pos+1:pos+cm1(k)))
     j = j + cm3(k)
     pos = pos + cm1(k)
  END DO

  !-----------------------------
  ! 4) lambda  !-----------------------------
  ! 4) lambda con loadmask:
  !    per ogni fattore qq: primo item con loadmask=1 fissato a 1
  !-----------------------------
  lambda = 0.d0
  pos = ntre0

  DO qq = 1, Q
     ref_item = 0
     DO k = 1, nvar
        IF (loadmask(k,qq) == 1) THEN
           ref_item = k
           EXIT
        END IF
     END DO
     IF (ALLOCATED(cfa_marker_item)) ref_item = cfa_marker_item(qq)
     IF (ref_item == 0 .OR. loadmask(ref_item,qq) /= 1) THEN
        CALL rexit('ERROR (loglik_long): factor')
     END IF

     lambda(ref_item,qq) = 1.d0

     DO k = 1, nvar
        IF (loadmask(k,qq) == 0) THEN
           lambda(k,qq) = 0.d0
        ELSE
           IF (k == ref_item) CYCLE
           pos = pos + 1
           lambda(k,qq) = b0(pos)
        END IF
     END DO
  END DO

  IF (Q > 1) THEN
   lambda_act = 0.d0
   DO k = 1, nvar
      DO ll = 1, nload_item(k)
         qq = load_idx(ll, k)
         lambda_act(ll, k) = lambda(k, qq)
      END DO
   END DO
END IF

  !-----------------------------
  ! 5) gamma
  !-----------------------------
  DO k = 1, nvar
     pos = pos + 1
     gamma(k) = b0(pos)
  END DO

  !-----------------------------
  ! 6) theta (solo ordinali; binari fissati a 1)
  !-----------------------------
  IF (theta_time_flag == 0) THEN
   DO k = 1, nvar
      IF (cm1(k) == 1) THEN
         theta_k(k) = 1.d0
         invsd_k(k) = 1.d0
         DO t = 1, n
            theta_long(k,t) = 1.d0
            invsd_long(k,t) = 1.d0
         END DO
      ELSE
         pos = pos + 1
         theta_k(k) = EXP(b0(pos))
         invsd_k(k) = 1.d0 / SQRT(theta_k(k))
         DO t = 1, n
            theta_long(k,t) = theta_k(k)
            invsd_long(k,t) = invsd_k(k)
         END DO
      END IF
   END DO
  ELSE
    DO k = 1, nvar
        IF (cm1(k) == 1) THEN
          theta_k(k) = 1.d0
          invsd_k(k) = 1.d0
          DO t = 1, n
              theta_long(k,t) = 1.d0
              invsd_long(k,t) = 1.d0
          END DO
        ELSE
          DO t = 1, n
              pos = pos + 1
              theta_long(k,t) = EXP(b0(pos))
              invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
          END DO
          theta_k(k) = theta_long(k,1)
          invsd_k(k) = invsd_long(k,1)
        END IF
    END DO
  END IF

  !-----------------------------
  ! 7) sigmak
  !-----------------------------
  IF (re_item_flag == 1) THEN
     DO k = 1, nvar
        pos = pos + 1
        sigmak(k)      = EXP(b0(pos))
        sqrt_sigmak(k) = SQRT(sigmak(k))
     END DO
  ELSE
     DO k = 1, nvar
        sigmak(k)      = 0.d0
        sqrt_sigmak(k) = 0.d0
     END DO
  END IF

  !-----------------------------
  ! 8) kappa: kappa(t=1,qq)=0 per ogni qq
  !-----------------------------
  DO qq = 1, Q
     kappa( (1-1)*Q + qq ) = 0.d0
  END DO

  IF (n > 1) THEN
     DO t = 2, n
        DO qq = 1, Q
           pos = pos + 1
           kappa( (t-1)*Q + qq ) = b0(pos)
        END DO
     END DO
  END IF


  ! factor regressions on covariates with time-invariant beta
  CALL ensure_cfa_beta_storage()
  IF (p_fix > 0) THEN
     DO qq = 1, Q
        DO r = 1, p_fix
           pos = pos + 1
           cfa_beta(qq,r) = b0(pos)
        END DO
     END DO
  END IF

  !-----------------------------
  ! 9) L_phi / AR1 params
  !-----------------------------
  L_phi = 0.d0

  IF (phi_struct_flag == 0) THEN
    ! L_phi libero
    DO t = 1, dimEta
      pos = pos + 1
      L_phi(t,t) = EXP(b0(pos))
    END DO
    DO col = 1, dimEta
      DO row = col+1, dimEta
          pos = pos + 1
          L_phi(row,col) = b0(pos)
      END DO
    END DO

  ELSE IF (phi_struct_flag == 1) THEN
    ! AR1: costruisco L_phi esplicito, A piena LIBERA, Sigma1 diag


    DO col = 1, Q
      DO row = 1, Q
        pos = pos + 1
        A_tmp(row,col) = b0(pos)
      END DO
    END DO

    DO row = 1, Q
      pos = pos + 1
      sqrt1_tmp(row) = SQRT(EXP(b0(pos)))
    END DO

    CALL build_Lphi_ar1_full(n, Q, A_tmp, sqrt1_tmp, L_phi)


  ELSE IF (phi_struct_flag == 2) THEN
    ! AR1 via ricorsione: salvo A_ar e sqrt_sig1, NON costruisco L_phi
    IF (ALLOCATED(A_ar)) DEALLOCATE(A_ar)
    IF (ALLOCATED(sqrt_sig1)) DEALLOCATE(sqrt_sig1)
    ALLOCATE(A_ar(Q,Q), sqrt_sig1(Q))

    DO col = 1, Q
      DO row = 1, Q
        pos = pos + 1
        A_ar(row,col) = b0(pos)
      END DO
    END DO

    DO row = 1, Q
      pos = pos + 1
      sqrt_sig1(row) = SQRT(EXP(b0(pos)))
    END DO

    L_phi = 0.d0

  ELSE
    CALL rexit("ERROR: phi_struct_flag non valido")
  END IF

  !-----------------------------
  ! 10) expected npm (warning)
  !-----------------------------
  n_theta_est = 0
  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) n_theta_est = n_theta_est + 1
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) n_theta_est = n_theta_est + n
     END DO
  END IF

  n_lambda_free = 0
  DO qq = 1, Q
     DO k = 1, nvar
        IF (loadmask(k,qq) == 1) n_lambda_free = n_lambda_free + 1
     END DO
  END DO
  n_lambda_free = n_lambda_free - Q

  IF (phi_struct_flag == 0) THEN
     n_phi_par = dim_eta*(dim_eta+1)/2
  ELSE
     n_phi_par = Q*Q + Q
  END IF

  expected_npm = ntre0                            &
              + n_lambda_free                     &
              + nvar                              &
              + n_theta_est                       &
              + MERGE(nvar, 0, re_item_flag == 1) &
              + (n-1)*Q                           &
              + Q*p_fix                          &
              + n_phi_par

  IF (npm0 /= expected_npm) THEN
     CALL rexit("The supplied parameter vector has the wrong length for this CFA model.")
  END IF

  ! ===== TRACE leggibile: vettori/matrici in blocchi (3 decimali) =====

  ! ! ---- b0 (vettore)
  ! WRITE(*,'(A)') "b0:"
  ! DO jj = 1, npm0, 8
  !   jj2 = MIN(jj+7, npm0)
  !   WRITE(*,'(*(F10.3,1X))') b0(jj:jj2)
  ! END DO

  ! ! ---- tau (vettore)
  ! WRITE(*,'(A,I0)') "tau_len = ", SIZE(tau)
  ! WRITE(*,'(A)') "tau:"
  ! DO jj = 1, SIZE(tau), 8
  !   jj2 = MIN(jj+7, SIZE(tau))
  !   WRITE(*,'(*(F10.3,1X))') tau(jj:jj2)
  ! END DO

  ! ! ---- lambda (matrice nvar x Q) : una riga per item
  ! WRITE(*,'(A)') "lambda:"
  ! DO rr = 1, nvar
  !   WRITE(*,'(*(F10.3,1X))') (lambda(rr,cc), cc=1, Q)
  ! END DO

  ! ! ---- gamma (vettore)
  ! WRITE(*,'(A)') "gamma:"
  ! DO jj = 1, nvar, 8
  !   jj2 = MIN(jj+7, nvar)
  !   WRITE(*,'(*(F10.3,1X))') gamma(jj:jj2)
  ! END DO

  ! ! ---- theta_long (matrice nvar x n) : una riga per item
  ! WRITE(*,'(A)') "theta_long:"
  ! DO rr = 1, nvar
  !   WRITE(*,'(*(F10.3,1X))') (theta_long(rr,cc), cc=1, n)
  ! END DO

  ! ! ---- sigmak (vettore)
  ! WRITE(*,'(A)') "sigmak:"
  ! DO jj = 1, nvar, 8
  !   jj2 = MIN(jj+7, nvar)
  !   WRITE(*,'(*(F10.3,1X))') sigmak(jj:jj2)
  ! END DO

  ! ! ---- kappa (vettore dimEta)
  ! WRITE(*,'(A,I0)') "kappa_len = ", dimEta
  ! WRITE(*,'(A)') "kappa:"
  ! DO jj = 1, dimEta, 8
  !   jj2 = MIN(jj+7, dimEta)
  !   WRITE(*,'(*(F10.3,1X))') kappa(jj:jj2)
  ! END DO

  ! ! ---- struttura phi
  ! WRITE(*,'(A,I0)') "phi_struct_flag = ", phi_struct_flag
  ! IF (phi_struct_flag == 0 .OR. phi_struct_flag == 1) THEN
  !   WRITE(*,'(A)') "L_phi:"
  !   DO rr = 1, dimEta
  !     WRITE(*,'(*(F10.3,1X))') (L_phi(rr,cc), cc=1, dimEta)
  !   END DO
  ! ELSE IF (phi_struct_flag == 2) THEN
  !   WRITE(*,'(A)') "A_ar:"
  !   DO rr = 1, Q
  !     WRITE(*,'(*(F10.3,1X))') (A_ar(rr,cc), cc=1, Q)
  !   END DO
  !   WRITE(*,'(A)') "sqrt_sig1:"
  !   DO jj = 1, Q, 8
  !     jj2 = MIN(jj+7, Q)
  !     WRITE(*,'(*(F10.3,1X))') sqrt_sig1(jj:jj2)
  !   END DO
  ! END IF

  ! WRITE(*,'(A)') "----"
  ! ! ===== fine TRACE =====

  ! --- END TRACE ---

  !-----------------------------
  ! 11) -loglik
  !-----------------------------


  vrais = -vrais_long()
  ! print *, vrais



END SUBROUTINE loglik_long




!===========================================================
!======================  AGH PART  =========================
!===========================================================
SUBROUTINE loglik_long_cond_i(i, sm, dimZ, z, loglik_i, ok)
  USE shared,   ONLY: nvar, n, Q, re_item_flag, nmes, c, Y, nitem, theta_long, &
                      lambda, gamma, sqrt_sigmak, tau, kappa, L_phi, loadmask
  USE datindiv, ONLY: nmescur
  USE utilities, ONLY: alnorm_fast
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, sm, dimZ
  DOUBLE PRECISION, DIMENSION(dimZ), INTENT(IN) :: z
  DOUBLE PRECISION, INTENT(OUT) :: loglik_i
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: k_item, j_mes, idx_obs, yijk, tuk
  INTEGER :: t_idx, kk_thr, qq, s_idx, idx_tq
  INTEGER :: dimEta

  DOUBLE PRECISION :: mu_ijk, linpred
  DOUBLE PRECISION :: var_lat_k, invsd_lat_k
  DOUBLE PRECISION :: tausup, tauinf, p
  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9
  DOUBLE PRECISION :: logp

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eta_vec, a_item, z_eta

  ok = .TRUE.
  dimEta = n * Q

  IF (re_item_flag == 1) THEN
     IF (dimZ /= dimEta + nvar) THEN
        CALL rexit('ERROR (loglik_long_cond_i): dimZ must be n*Q + nvar')
     END IF
  ELSE
     IF (dimZ /= dimEta) THEN
        CALL rexit('ERROR (loglik_long_cond_i): dimZ must be n*Q')
     END IF
  END IF

  ALLOCATE(eta_vec(dimEta), a_item(nvar), z_eta(dimEta))

  z_eta(1:dimEta) = z(1:dimEta)

  CALL cfa_eta_mean(i, eta_vec)
  DO s_idx = 1, dimEta
     DO qq = 1, s_idx
        eta_vec(s_idx) = eta_vec(s_idx) + L_phi(s_idx,qq) * z_eta(qq)
     END DO
  END DO

  IF (re_item_flag == 1) THEN
     DO k_item = 1, nvar
        a_item(k_item) = sqrt_sigmak(k_item) * z(dimEta + k_item)
     END DO
  ELSE
     a_item = 0.d0
  END IF

  loglik_i = 0.d0
  kk_thr   = 0
  idx_obs  = nmescur

  DO k_item = 1, nvar

     IF (nmes(i,k_item) == 0) THEN
        kk_thr = kk_thr + (c(k_item) - 1)
        CYCLE
     END IF

     DO j_mes = 1, nmes(i,k_item)

        idx_obs = idx_obs + 1
        t_idx   = nitem(idx_obs)

        var_lat_k = theta_at(k_item, t_idx)


      IF (var_lat_k <= 0.d0 .OR. var_lat_k /= var_lat_k) THEN
        loglik_i = LOG_ZERO
        ok = .TRUE.
        DEALLOCATE(eta_vec, a_item, z_eta)
        RETURN
      END IF

        invsd_lat_k = 1.d0 / SQRT(var_lat_k)

        linpred = 0.d0
        DO qq = 1, Q
           IF (loadmask(k_item,qq) == 0) CYCLE
           idx_tq  = (t_idx-1)*Q + qq
           linpred = linpred + lambda_at(k_item,qq,t_idx) * eta_vec(idx_tq)
        END DO

        mu_ijk = gamma_at(k_item,t_idx) + a_item(k_item) + linpred

        yijk = Y(idx_obs)
        tuk  = yijk + kk_thr

        IF (c(k_item) == 2) THEN
           tausup = (tau_at(kk_thr+1,t_idx) - mu_ijk) * invsd_lat_k

           IF (yijk == 1) THEN
              p = alnorm_fast(tausup, .FALSE.)
           ELSE
              p = 1.d0 - alnorm_fast(tausup, .FALSE.)
           END IF

        ELSE
           IF (yijk == 1) THEN
              tausup = (tau_at(tuk,t_idx) - mu_ijk) * invsd_lat_k
              p = alnorm_fast(tausup, .FALSE.)

           ELSE IF (yijk == c(k_item)) THEN
              tauinf = (tau_at(tuk-1,t_idx) - mu_ijk) * invsd_lat_k
              p = 1.d0 - alnorm_fast(tauinf, .FALSE.)

           ELSE
              tausup = (tau_at(tuk,t_idx)   - mu_ijk) * invsd_lat_k
              tauinf = (tau_at(tuk-1,t_idx) - mu_ijk) * invsd_lat_k
              p = alnorm_fast(tausup, .FALSE.) - alnorm_fast(tauinf, .FALSE.)
           END IF
        END IF

        IF (p > 0.d0) THEN
           logp = LOG(p)
           IF (logp == logp) THEN
              loglik_i = loglik_i + logp
           END IF

        ELSE IF (p == 0.d0) THEN
           loglik_i = LOG_ZERO

        ELSE
           ! p negativo o NaN ignorato
        END IF

     END DO

     kk_thr = kk_thr + (c(k_item) - 1)

  END DO

  DEALLOCATE(eta_vec, a_item, z_eta)

END SUBROUTINE loglik_long_cond_i



SUBROUTINE blockroot_tt_from_precision_long(dimEta, dimZ, H, B, logdetB, logdetEta, ok_root)
  USE shared, ONLY: nvar, re_item_flag
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: dimEta, dimZ
  DOUBLE PRECISION, DIMENSION(dimZ,dimZ), INTENT(IN) :: H
  DOUBLE PRECISION, DIMENSION(dimZ,dimZ), INTENT(OUT) :: B
  DOUBLE PRECISION, INTENT(OUT) :: logdetB, logdetEta
  LOGICAL, INTENT(OUT) :: ok_root

  INTEGER :: info, r, c, k, rr
  DOUBLE PRECISION :: dkk, sd_k
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: S_eta, V_eta, T_eta

  ok_root = .TRUE.
  B = 0.d0
  logdetB = 0.d0
  logdetEta = 0.d0

  IF (dimEta <= 0 .OR. dimZ < dimEta) THEN
     ok_root = .FALSE.
     RETURN
  END IF

  ALLOCATE(S_eta(dimEta,dimEta), V_eta(dimEta,dimEta), T_eta(dimEta,dimEta))

  S_eta = H(1:dimEta,1:dimEta)

  IF (re_item_flag == 1) THEN
     DO k = 1, nvar
        dkk = H(dimEta+k,dimEta+k)
        IF (dkk <= 0.d0 .OR. dkk /= dkk .OR. ABS(dkk) > 1.d300) THEN
           ok_root = .FALSE.
           DEALLOCATE(S_eta, V_eta, T_eta)
           RETURN
        END IF

        DO r = 1, dimEta
           DO c = 1, dimEta
              S_eta(r,c) = S_eta(r,c) - H(r,dimEta+k) * H(dimEta+k,c) / dkk
           END DO
        END DO
     END DO
  END IF

  S_eta = 0.5d0 * (S_eta + TRANSPOSE(S_eta))
  V_eta = S_eta

  CALL DPOTRF('U', dimEta, V_eta, dimEta, info)
  IF (info /= 0) THEN
     ok_root = .FALSE.
     DEALLOCATE(S_eta, V_eta, T_eta)
     RETURN
  END IF

  CALL DPOTRI('U', dimEta, V_eta, dimEta, info)
  IF (info /= 0) THEN
     ok_root = .FALSE.
     DEALLOCATE(S_eta, V_eta, T_eta)
     RETURN
  END IF

  DO r = 1, dimEta
     DO c = r+1, dimEta
        V_eta(c,r) = V_eta(r,c)
     END DO
  END DO
  V_eta = 0.5d0 * (V_eta + TRANSPOSE(V_eta))

  T_eta = V_eta
  CALL DPOTRF('L', dimEta, T_eta, dimEta, info)
  IF (info /= 0) THEN
     ok_root = .FALSE.
     DEALLOCATE(S_eta, V_eta, T_eta)
     RETURN
  END IF

  IF (dimEta > 1) THEN
     DO r = 1, dimEta-1
        T_eta(r,r+1:dimEta) = 0.d0
     END DO
  END IF

  B(1:dimEta,1:dimEta) = T_eta

  DO r = 1, dimEta
     IF (T_eta(r,r) <= 0.d0 .OR. T_eta(r,r) /= T_eta(r,r)) THEN
        ok_root = .FALSE.
        DEALLOCATE(S_eta, V_eta, T_eta)
        RETURN
     END IF
     logdetEta = logdetEta + LOG(T_eta(r,r))
  END DO

  logdetB = logdetEta

  IF (re_item_flag == 1) THEN
     DO k = 1, nvar
        dkk = H(dimEta+k,dimEta+k)
        sd_k = 1.d0 / SQRT(dkk)

        IF (sd_k <= 0.d0 .OR. sd_k /= sd_k .OR. ABS(sd_k) > 1.d300) THEN
           ok_root = .FALSE.
           DEALLOCATE(S_eta, V_eta, T_eta)
           RETURN
        END IF

        DO c = 1, dimEta
           B(dimEta+k,c) = 0.d0
           DO rr = 1, dimEta
              B(dimEta+k,c) = B(dimEta+k,c) - H(dimEta+k,rr) * T_eta(rr,c) / dkk
           END DO
        END DO

        B(dimEta+k,dimEta+k) = sd_k
        logdetB = logdetB + LOG(sd_k)
     END DO
  END IF

  DEALLOCATE(S_eta, V_eta, T_eta)

END SUBROUTINE blockroot_tt_from_precision_long


SUBROUTINE find_mode_and_Btt_long(i, sm, dimZ, z_start, z_mode, B, logdetB, logdetEta, ok_mode)
  USE shared, ONLY: n, Q
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, sm, dimZ
  DOUBLE PRECISION, DIMENSION(dimZ), INTENT(IN) :: z_start
  DOUBLE PRECISION, DIMENSION(dimZ), INTENT(OUT) :: z_mode
  DOUBLE PRECISION, DIMENSION(dimZ,dimZ), INTENT(OUT) :: B
  DOUBLE PRECISION, INTENT(OUT) :: logdetB, logdetEta
  LOGICAL, INTENT(OUT) :: ok_mode

  INTEGER :: dimEta
  DOUBLE PRECISION :: obj, logdet_dummy
  LOGICAL :: ok_hess, ok_root
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: grad
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: H, R_dummy

  ok_mode = .TRUE.
  dimEta = n * Q

  ALLOCATE(grad(dimZ), H(dimZ,dimZ), R_dummy(dimZ,dimZ))

  CALL find_mode_and_R_long(i, sm, dimZ, z_start, z_mode, R_dummy, logdet_dummy, ok_mode)

  IF (.NOT. ok_mode) THEN
     B = 0.d0
     logdetB = 0.d0
     logdetEta = 0.d0
     DEALLOCATE(grad, H, R_dummy)
     RETURN
  END IF

  CALL grad_hess_logpost_long_i(i, sm, dimZ, z_mode, obj, grad, H, ok_hess)

  IF (.NOT. ok_hess .OR. obj /= obj .OR. ABS(obj) > 1.d300) THEN
     ok_mode = .FALSE.
     B = 0.d0
     logdetB = 0.d0
     logdetEta = 0.d0
     DEALLOCATE(grad, H, R_dummy)
     RETURN
  END IF

  H = 0.5d0 * (H + TRANSPOSE(H))

  CALL blockroot_tt_from_precision_long(dimEta, dimZ, H, B, logdetB, logdetEta, ok_root)

  IF (.NOT. ok_root) THEN
     ok_mode = .FALSE.
     B = 0.d0
     logdetB = 0.d0
     logdetEta = 0.d0
  END IF

  DEALLOCATE(grad, H, R_dummy)

END SUBROUTINE find_mode_and_Btt_long

SUBROUTINE build_agh_grid_long(ok_grid)
  USE shared,   ONLY: nsub, nvar, n, Q, nquad, nodes, re_item_flag, nmes, &
                      agh_grid_ready, agh_dimZ_cached, agh_nquad_cached,  &
                      agh_nsub_cached, agh_z, agh_logconst, agh_modes
  USE datindiv, ONLY: nmescur
  IMPLICIT NONE

  INTEGER :: i, sm, l, d, dimZ, dimEta, baseZ
  DOUBLE PRECISION :: logdetR

  LOGICAL, INTENT(OUT) :: ok_grid
  LOGICAL :: ok_mode
  DOUBLE PRECISION, PARAMETER :: BAD_LOG = -1.d9
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: z_mode, z, zeta, z_start
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R

  ok_grid = .TRUE.

  dimEta = n * Q

  IF (re_item_flag == 1) THEN
     dimZ = dimEta + nvar
  ELSE
     dimZ = dimEta
  END IF

  IF (ALLOCATED(agh_z)) THEN
     IF (agh_dimZ_cached /= dimZ .OR. agh_nquad_cached /= nquad .OR. &
         agh_nsub_cached /= nsub) THEN
        DEALLOCATE(agh_z)
        DEALLOCATE(agh_logconst)
        DEALLOCATE(agh_modes)
        agh_grid_ready = .FALSE.
     END IF
  END IF

  IF (.NOT. ALLOCATED(agh_z)) THEN
     ALLOCATE(agh_z(dimZ, nquad, nsub))
     ALLOCATE(agh_logconst(nsub))
     ALLOCATE(agh_modes(dimZ, nsub))

     agh_dimZ_cached  = dimZ
     agh_nquad_cached = nquad
     agh_nsub_cached  = nsub
  END IF

  ALLOCATE(z_mode(dimZ), z(dimZ), zeta(dimZ), z_start(dimZ), R(dimZ, dimZ))

  nmescur = 0

  DO i = 1, nsub

     sm = SUM(nmes(i,:))

     IF (agh_grid_ready) THEN
        z_start = agh_modes(:, i)
     ELSE
        z_start = 0.d0
     END IF

     CALL find_mode_and_R_long(i, sm, dimZ, z_start, z_mode, R, logdetR, ok_mode)

      IF (.NOT. ok_mode) THEN
        ok_grid = .FALSE.
        agh_grid_ready = .FALSE.
        DEALLOCATE(z_mode, z, zeta, z_start, R)
        RETURN
      END IF
      agh_logconst(i) = logdetR
      agh_modes(:, i) = z_mode

     DO l = 1, nquad

        baseZ = (l - 1) * dimZ

        DO d = 1, dimZ
           zeta(d) = nodes(baseZ + d)
        END DO

        z = z_mode

        CALL DGEMV('N', dimZ, dimZ, 1.d0, R, dimZ, zeta, 1, 1.d0, z, 1)

        agh_z(:, l, i) = z

     END DO

     nmescur = nmescur + sm

  END DO

  agh_grid_ready = .TRUE.

  DEALLOCATE(z_mode, z, zeta, z_start, R)

END SUBROUTINE build_agh_grid_long


SUBROUTINE grad_hess_logpost_long_i(i, sm, dimZ, z, f, grad, H, ok)
  USE shared,   ONLY: nvar, n, Q, re_item_flag, nmes, c, Y, nitem, theta_long, &
                      lambda, gamma, sqrt_sigmak, tau, kappa, L_phi, loadmask
  USE datindiv, ONLY: nmescur
  USE utilities, ONLY: alnorm_fast
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, sm, dimZ
  DOUBLE PRECISION, DIMENSION(dimZ),      INTENT(IN)  :: z
  DOUBLE PRECISION,                        INTENT(OUT) :: f
  DOUBLE PRECISION, DIMENSION(dimZ),      INTENT(OUT) :: grad
  DOUBLE PRECISION, DIMENSION(dimZ,dimZ), INTENT(OUT) :: H
  LOGICAL, INTENT(OUT) :: ok

  DOUBLE PRECISION, PARAMETER :: INV_SQRT_2PI = 0.39894228040143267794d0

  INTEGER :: k_item, j_mes, idx_obs, yijk, tuk
  INTEGER :: t_idx, kk_thr
  INTEGER :: d1, d2, s_idx, qq, idx_tq
  INTEGER :: dimEta

  DOUBLE PRECISION :: mu_ijk, linpred
  DOUBLE PRECISION :: var_lat_k, invsd_lat_k
  DOUBLE PRECISION :: alpha, beta_thr, p
  DOUBLE PRECISION :: phi_a, phi_b
  DOUBLE PRECISION :: dP_dmu, d2P_dmu2
  DOUBLE PRECISION :: s_mu, h_mu
  DOUBLE PRECISION :: loglik_i, qz
  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9
  DOUBLE PRECISION :: logp

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eta_vec, a_item, z_eta
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dmu_dz

  ok = .TRUE.

  dimEta = n * Q

  IF (re_item_flag == 1) THEN
     IF (dimZ /= dimEta + nvar) THEN
        CALL rexit('ERROR (grad_hess_logpost_long_i): dimZ must be n*Q + nvar')
     END IF
  ELSE
     IF (dimZ /= dimEta) THEN
        CALL rexit('ERROR (grad_hess_logpost_long_i): dimZ must be n*Q')
     END IF
  END IF

  ALLOCATE(eta_vec(dimEta), a_item(nvar), z_eta(dimEta))
  ALLOCATE(dmu_dz(dimZ))

  z_eta(1:dimEta) = z(1:dimEta)

  CALL cfa_eta_mean(i, eta_vec)
  DO s_idx = 1, dimEta
     DO qq = 1, s_idx
        eta_vec(s_idx) = eta_vec(s_idx) + L_phi(s_idx,qq) * z_eta(qq)
     END DO
  END DO

  IF (re_item_flag == 1) THEN
     DO k_item = 1, nvar
        a_item(k_item) = sqrt_sigmak(k_item) * z(dimEta + k_item)
     END DO
  ELSE
     a_item = 0.d0
  END IF

  grad     = 0.d0
  H        = 0.d0
  loglik_i = 0.d0

  kk_thr  = 0
  idx_obs = nmescur

  DO k_item = 1, nvar

     IF (nmes(i,k_item) == 0) THEN
        kk_thr = kk_thr + (c(k_item) - 1)
        CYCLE
     END IF

     DO j_mes = 1, nmes(i,k_item)

        idx_obs = idx_obs + 1
        t_idx   = nitem(idx_obs)


        var_lat_k = theta_at(k_item, t_idx)

        IF (var_lat_k <= 0.d0 .OR. var_lat_k /= var_lat_k) THEN
          f = 1.d9
          grad = 0.d0
          H = 0.d0

          DO d1 = 1, dimZ
              H(d1,d1) = 1.d0
          END DO

          ok = .TRUE.
          DEALLOCATE(eta_vec, a_item, z_eta, dmu_dz)
          RETURN
        END IF

        invsd_lat_k = 1.d0 / SQRT(var_lat_k)

        linpred = 0.d0
        DO qq = 1, Q
           IF (loadmask(k_item,qq) == 0) CYCLE
           idx_tq  = (t_idx-1)*Q + qq
           linpred = linpred + lambda_at(k_item,qq,t_idx) * eta_vec(idx_tq)
        END DO

        mu_ijk = gamma_at(k_item,t_idx) + a_item(k_item) + linpred

        yijk = Y(idx_obs)
        tuk  = yijk + kk_thr

        IF (c(k_item) == 2) THEN

           alpha = (tau_at(kk_thr+1,t_idx) - mu_ijk) * invsd_lat_k
           phi_a = INV_SQRT_2PI * EXP(-0.5d0 * alpha * alpha)

           IF (yijk == 1) THEN
              p = alnorm_fast(alpha, .FALSE.)

              dP_dmu   = -phi_a * invsd_lat_k
              d2P_dmu2 = -alpha * phi_a * invsd_lat_k * invsd_lat_k

           ELSE
              p = alnorm_fast(alpha, .TRUE.)

              dP_dmu   =  phi_a * invsd_lat_k
              d2P_dmu2 =  alpha * phi_a * invsd_lat_k * invsd_lat_k
           END IF

        ELSE

           IF (yijk == 1) THEN

              alpha = (tau_at(tuk,t_idx) - mu_ijk) * invsd_lat_k
              p = alnorm_fast(alpha, .FALSE.)

              phi_a    = INV_SQRT_2PI * EXP(-0.5d0 * alpha * alpha)
              dP_dmu   = -phi_a * invsd_lat_k
              d2P_dmu2 = -alpha * phi_a * invsd_lat_k * invsd_lat_k

           ELSE IF (yijk == c(k_item)) THEN

              beta_thr = (tau_at(tuk-1,t_idx) - mu_ijk) * invsd_lat_k
              p = alnorm_fast(beta_thr, .TRUE.)

              phi_b    = INV_SQRT_2PI * EXP(-0.5d0 * beta_thr * beta_thr)
              dP_dmu   =  phi_b * invsd_lat_k
              d2P_dmu2 =  beta_thr * phi_b * invsd_lat_k * invsd_lat_k

           ELSE

              alpha    = (tau_at(tuk,t_idx)   - mu_ijk) * invsd_lat_k
              beta_thr = (tau_at(tuk-1,t_idx) - mu_ijk) * invsd_lat_k

               IF (beta_thr > 0.d0) THEN
                  p = alnorm_fast(beta_thr, .TRUE.) - alnorm_fast(alpha, .TRUE.)
               ELSE
                  p = alnorm_fast(alpha, .FALSE.) - alnorm_fast(beta_thr, .FALSE.)
               END IF

              phi_a    = INV_SQRT_2PI * EXP(-0.5d0 * alpha * alpha)
              phi_b    = INV_SQRT_2PI * EXP(-0.5d0 * beta_thr * beta_thr)

              dP_dmu   = (phi_b - phi_a) * invsd_lat_k
              d2P_dmu2 = (beta_thr * phi_b - alpha * phi_a) * invsd_lat_k * invsd_lat_k

           END IF
        END IF

       IF (p > 0.d0) THEN
          logp = LOG(p)

          IF (logp == logp) THEN
              loglik_i = loglik_i + logp
          ELSE
              CYCLE
          END IF

        ELSE IF (p == 0.d0) THEN
          loglik_i = LOG_ZERO
          CYCLE

        ELSE
          ! p negativo o NaN ignorato
          CYCLE
        END IF

        s_mu = dP_dmu / p
        h_mu = d2P_dmu2 / p - s_mu * s_mu

        IF (s_mu /= s_mu .OR. h_mu /= h_mu) THEN
          CYCLE
        END IF


        dmu_dz = 0.d0

        DO qq = 1, Q
           IF (loadmask(k_item,qq) == 0) CYCLE
           idx_tq = (t_idx-1)*Q + qq
           DO s_idx = 1, idx_tq
              dmu_dz(s_idx) = dmu_dz(s_idx) + lambda_at(k_item,qq,t_idx) * L_phi(idx_tq, s_idx)
           END DO
        END DO

        IF (re_item_flag == 1) THEN
           d1 = dimEta + k_item
           dmu_dz(d1) = dmu_dz(d1) + sqrt_sigmak(k_item)
        END IF

        DO d1 = 1, dimZ
           IF (dmu_dz(d1) == 0.d0) CYCLE

           grad(d1) = grad(d1) + s_mu * dmu_dz(d1)

           DO d2 = 1, dimZ
              IF (dmu_dz(d2) == 0.d0) CYCLE
              H(d1,d2) = H(d1,d2) + h_mu * dmu_dz(d1) * dmu_dz(d2)
           END DO
        END DO

     END DO

     kk_thr = kk_thr + (c(k_item) - 1)

  END DO

  qz = DOT_PRODUCT(z, z)

  f = -loglik_i + 0.5d0 * qz

  grad = -grad
  H    = -H

  DO d1 = 1, dimZ
     grad(d1) = grad(d1) + z(d1)
     H(d1,d1) = H(d1,d1) + 1.d0
  END DO

  DEALLOCATE(eta_vec, a_item, z_eta, dmu_dz)

END SUBROUTINE grad_hess_logpost_long_i


SUBROUTINE cd_hessian_score_logpost_long_i(i, sm, dimZ, z, Hcd)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, sm, dimZ
  DOUBLE PRECISION, DIMENSION(dimZ), INTENT(IN) :: z
  DOUBLE PRECISION, DIMENSION(dimZ,dimZ), INTENT(OUT) :: Hcd

  DOUBLE PRECISION, PARAMETER :: eps = 0.001d0

  INTEGER :: d, r, c
  DOUBLE PRECISION :: ex, diffx
  DOUBLE PRECISION :: f1, f2
  LOGICAL :: ok1, ok2

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: z1, z2, g1, g2
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Htmp, res

  ALLOCATE(z1(dimZ), z2(dimZ), g1(dimZ), g2(dimZ))
  ALLOCATE(Htmp(dimZ,dimZ), res(dimZ,dimZ))

  res = 0.d0

  DO d = 1, dimZ

     z1 = z
     z2 = z

     ex = MAX(ABS(z(d)), 1.d0)

     z1(d) = z(d) + eps * ex
     z2(d) = z(d) - eps * ex

     CALL grad_hess_logpost_long_i(i, sm, dimZ, z1, f1, g1, Htmp, ok1)
     CALL grad_hess_logpost_long_i(i, sm, dimZ, z2, f2, g2, Htmp, ok2)

     IF (.NOT. ok1 .OR. .NOT. ok2) THEN
        CALL rexit('ERROR (cd_hessian_score_logpost_long_i): invalid finite difference value for subject')
     END IF

     IF (f1 /= f1 .OR. f2 /= f2) THEN
        CALL rexit('ERROR (cd_hessian_score_logpost_long_i): NaN finite difference value for subject')
     END IF

     diffx = z1(d) - z2(d)

     DO r = 1, dimZ
        res(r,d) = (g1(r) - g2(r)) / diffx
     END DO

  END DO

  DO r = 1, dimZ
     DO c = 1, dimZ
        Hcd(r,c) = 0.5d0 * (res(r,c) + res(c,r))
     END DO
  END DO

  DEALLOCATE(z1, z2, g1, g2, Htmp, res)

END SUBROUTINE cd_hessian_score_logpost_long_i

SUBROUTINE find_mode_and_R_long(i, sm, dimZ, z_start, z_mode, R, logdetR, ok_mode)
  USE shared
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE vmmin_long_mode(i, sm, dimZ, x, fmin, fail, fncount, grcount, &
                                maxit, abstol, reltol, trace) BIND(C, name="vmmin_long_mode")
       USE, INTRINSIC :: ISO_C_BINDING
       IMPLICIT NONE

       INTEGER(C_INT), INTENT(IN) :: i, sm, dimZ
       REAL(C_DOUBLE), DIMENSION(*), INTENT(INOUT) :: x
       REAL(C_DOUBLE), INTENT(OUT) :: fmin
       INTEGER(C_INT), INTENT(OUT) :: fail, fncount, grcount
       INTEGER(C_INT), INTENT(IN) :: maxit, trace
       REAL(C_DOUBLE), INTENT(IN) :: abstol, reltol
     END SUBROUTINE vmmin_long_mode
  END INTERFACE

  INTEGER, INTENT(IN) :: i, sm, dimZ
  DOUBLE PRECISION, DIMENSION(dimZ), INTENT(IN) :: z_start
  DOUBLE PRECISION, DIMENSION(dimZ), INTENT(OUT) :: z_mode
  DOUBLE PRECISION, DIMENSION(dimZ,dimZ), INTENT(OUT) :: R
  DOUBLE PRECISION, INTENT(OUT) :: logdetR
  LOGICAL, INTENT(OUT) :: ok_mode

  INTEGER :: info, j, jj
  LOGICAL :: ok
  DOUBLE PRECISION :: obj
  DOUBLE PRECISION, DIMENSION(dimZ) :: g
  DOUBLE PRECISION, DIMENSION(dimZ,dimZ) :: A

  INTEGER(C_INT) :: i_c, sm_c, dimZ_c
  INTEGER(C_INT) :: fail_vmmin, fncount_vmmin, grcount_vmmin
  INTEGER(C_INT) :: maxit_vmmin, trace_vmmin
  REAL(C_DOUBLE) :: abstol_vmmin, reltol_vmmin
  REAL(C_DOUBLE) :: obj_vmmin

  ok_mode = .TRUE.
  z_mode = z_start
  R = 0.d0
  DO j = 1, dimZ
     R(j,j) = 1.d0
  END DO
  logdetR = 0.d0

  i_c    = INT(i, C_INT)
  sm_c   = INT(sm, C_INT)
  dimZ_c = INT(dimZ, C_INT)

  maxit_vmmin  = 100_C_INT
  trace_vmmin  = 0_C_INT
  abstol_vmmin = -HUGE(1.0D0)
  reltol_vmmin = SQRT(EPSILON(1.0D0))

  CALL vmmin_long_mode(i_c, sm_c, dimZ_c, z_mode, obj_vmmin, fail_vmmin, &
                       fncount_vmmin, grcount_vmmin, maxit_vmmin,        &
                       abstol_vmmin, reltol_vmmin, trace_vmmin)

  IF (fail_vmmin /= 0_C_INT) THEN
     CALL intpr("Adaptive mode search stopped early for subject", -1, (/ i /), 1)
  END IF

  IF (obj_vmmin /= obj_vmmin .OR. ABS(obj_vmmin) > 1.d300) THEN
     ok_mode = .FALSE.
     z_mode = z_start
     R = 0.d0
     DO j = 1, dimZ
        R(j,j) = 1.d0
     END DO
     logdetR = 0.d0
     RETURN
  END IF

  CALL grad_hess_logpost_long_i(i, sm, dimZ, z_mode, obj, g, A, ok)

  IF (.NOT. ok .OR. obj /= obj .OR. ABS(obj) > 1.d300) THEN
     ok_mode = .FALSE.
     z_mode = z_start
     R = 0.d0
     DO j = 1, dimZ
        R(j,j) = 1.d0
     END DO
     logdetR = 0.d0
     RETURN
  END IF

  A = 0.5d0 * (A + TRANSPOSE(A))

  CALL DPOTRF('U', dimZ, A, dimZ, info)

  IF (info /= 0) THEN
     ok_mode = .FALSE.
     z_mode = z_start
     R = 0.d0
     DO j = 1, dimZ
        R(j,j) = 1.d0
     END DO
     logdetR = 0.d0
     RETURN
  END IF

  logdetR = 0.d0

  DO j = 1, dimZ
     IF (A(j,j) <= 0.d0 .OR. A(j,j) /= A(j,j) .OR. ABS(A(j,j)) > 1.d300) THEN
        ok_mode = .FALSE.
        z_mode = z_start
        R = 0.d0
        DO jj = 1, dimZ
           R(jj,jj) = 1.d0
        END DO
        logdetR = 0.d0
        RETURN
     END IF

     logdetR = logdetR - LOG(A(j,j))
  END DO

  R = A

  IF (dimZ > 1) THEN
     DO j = 1, dimZ - 1
        R(j+1:dimZ,j) = 0.d0
     END DO
  END IF

  CALL DTRTRI('U', 'N', dimZ, R, dimZ, info)

  IF (info /= 0) THEN
     ok_mode = .FALSE.
     z_mode = z_start
     R = 0.d0
     DO j = 1, dimZ
        R(j,j) = 1.d0
     END DO
     logdetR = 0.d0
     RETURN
  END IF

END SUBROUTINE find_mode_and_R_long

SUBROUTINE fg_long_for_vmmin(i_c, sm_c, dimZ_c, z, f, g, ok_c) BIND(C, name="fg_long_for_vmmin")
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE

  INTEGER(C_INT), INTENT(IN) :: i_c, sm_c, dimZ_c
  REAL(C_DOUBLE), DIMENSION(dimZ_c), INTENT(IN) :: z
  REAL(C_DOUBLE), INTENT(OUT) :: f
  REAL(C_DOUBLE), DIMENSION(dimZ_c), INTENT(OUT) :: g
  INTEGER(C_INT), INTENT(OUT) :: ok_c

  INTEGER :: i, sm, dimZ
  LOGICAL :: ok
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Htmp

  i    = INT(i_c)
  sm   = INT(sm_c)
  dimZ = INT(dimZ_c)

  ALLOCATE(Htmp(dimZ, dimZ))

  CALL grad_hess_logpost_long_i(i, sm, dimZ, z, f, g, Htmp, ok)

  IF (ok) THEN
     ok_c = 1_C_INT
  ELSE
     ok_c = 0_C_INT
  END IF

  DEALLOCATE(Htmp)

END SUBROUTINE fg_long_for_vmmin


SUBROUTINE vrais_long_i_long_agh_fixed(i, sm, out, ok)
  USE shared,   ONLY: nvar, n, Q, nquad, nodes, weights, re_item_flag, &
                      agh_z, agh_logconst
  USE datindiv, ONLY: nmescur
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, sm
  DOUBLE PRECISION, INTENT(OUT) :: out
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: n_valid
  INTEGER :: l, d, dimZ, dimEta, baseZ
  LOGICAL :: ok_node

  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9

  DOUBLE PRECISION :: loglik_cond
  DOUBLE PRECISION :: qz, qzeta
  DOUBLE PRECISION :: max_term, sum_exp, term

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: z, zeta

  ok = .TRUE.
  out = 0.d0

  dimEta = n * Q

  IF (re_item_flag == 1) THEN
     dimZ = dimEta + nvar
  ELSE
     dimZ = dimEta
  END IF

  ALLOCATE(z(dimZ), zeta(dimZ))

  n_valid = 0
  max_term = 0.d0
  sum_exp  = 0.d0

  DO l = 1, nquad

     baseZ = (l - 1) * dimZ

     DO d = 1, dimZ
        zeta(d) = nodes(baseZ + d)
     END DO

     z = agh_z(:, l, i)

     CALL loglik_long_cond_i(i, sm, dimZ, z, loglik_cond, ok_node)

     IF (.NOT. ok_node) CYCLE
     IF (loglik_cond /= loglik_cond) CYCLE

     qz    = DOT_PRODUCT(z, z)
     qzeta = DOT_PRODUCT(zeta, zeta)

     IF (qz /= qz .OR. qzeta /= qzeta) CYCLE

     term = loglik_cond - 0.5d0 * (qz - qzeta) + LOG(weights(l))

     IF (term /= term) CYCLE

     IF (n_valid == 0) THEN
        max_term = term
        sum_exp  = 1.d0
        n_valid  = 1

     ELSE IF (term > max_term) THEN
        sum_exp  = sum_exp * EXP(max_term - term) + 1.d0
        max_term = term
        n_valid  = n_valid + 1

     ELSE IF (term == max_term) THEN
        sum_exp = sum_exp + 1.d0
        n_valid = n_valid + 1

     ELSE
        sum_exp = sum_exp + EXP(term - max_term)
        n_valid = n_valid + 1
     END IF

  END DO

  IF (n_valid == 0) THEN
     out = LOG_ZERO
     DEALLOCATE(z, zeta)
     RETURN
  END IF

  IF (sum_exp <= 0.d0 .OR. sum_exp /= sum_exp) THEN
     out = LOG_ZERO
     DEALLOCATE(z, zeta)
     RETURN
  END IF

  out = agh_logconst(i) + max_term + LOG(sum_exp)

  IF (out /= out) THEN
     out = LOG_ZERO
  END IF

  DEALLOCATE(z, zeta)

END SUBROUTINE vrais_long_i_long_agh_fixed

SUBROUTINE vrais_long_long_agh_fixed(out)
  USE shared,   ONLY: nsub, nmes, nvar, agh_grid_ready
  USE datindiv, ONLY: nmescur
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(OUT) :: out

  INTEGER :: i, sm
  DOUBLE PRECISION :: temp
  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9
  LOGICAL :: ok_i

  IF (.NOT. agh_grid_ready) THEN
     CALL rexit('ERROR (vrais_long_long_agh_fixed): AGH grid not ready')
  END IF

  nmescur = 0
  out = 0.d0

  DO i = 1, nsub

     sm = SUM(nmes(i,:))

     CALL vrais_long_i_long_agh_fixed(i, sm, temp, ok_i)

     IF (ok_i .AND. temp == temp) THEN
        out = out + temp
     END IF

     nmescur = nmescur + sm

  END DO

  IF (out /= out) THEN
     out = LOG_ZERO
  END IF

END SUBROUTINE vrais_long_long_agh_fixed

SUBROUTINE loglik_long_agh(b0, Y0, nvar0, nsub0, nobs0, n0,      &
                          nmes0, c0, cm10, ntre0,                &
                          nquad0, nitem0,                        &
                          nodes0, weights0,                      &
                          npm0, re_item0, theta_time0,           &
                          Q0, loadmask0, phi_struct0,            &
                          monotone_tau0, update_gh0,fit_id0,     &
                          vrais)                                 &
        BIND(C, name="loglik_long_agh")

  USE shared
  USE utilities
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: monotone_tau0
  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0
  INTEGER, INTENT(IN) :: ntre0, nquad0, npm0
  INTEGER, INTENT(IN) :: re_item0, theta_time0
  INTEGER, INTENT(IN) :: Q0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0
  INTEGER, INTENT(IN) :: update_gh0
  INTEGER, INTENT(IN) :: fit_id0
  INTEGER, INTENT(IN) :: phi_struct0

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.d9
  LOGICAL :: ok_grid
  DOUBLE PRECISION :: ll_agh

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0),       INTENT(IN) :: Y0, nitem0

  DOUBLE PRECISION, DIMENSION(*),      INTENT(IN) :: nodes0
  DOUBLE PRECISION, DIMENSION(nquad0), INTENT(IN) :: weights0
  DOUBLE PRECISION, DIMENSION(npm0),   INTENT(IN) :: b0

  DOUBLE PRECISION, INTENT(OUT) :: vrais

  INTEGER :: j, k, t, row, col, qq, ll, r
  INTEGER :: pos, expected_npm
  INTEGER :: dimZ, dimEta
  INTEGER :: n_theta_est, n_lambda_free, ref_item
  INTEGER :: ntre_calc

          ! monotonia lcmm-style: incrementi positivi via exp
        DOUBLE PRECISION :: prev
        INTEGER :: h, nf
  INTEGER, DIMENSION(:), ALLOCATABLE :: cm3
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ta

  INTEGER :: n_phi_par


  ! global
  nvar  = nvar0
  nsub  = nsub0
  nobs  = nobs0
  n     = n0
  nquad = nquad0
  npm   = npm0

  phi_struct_flag = phi_struct0

  re_item_flag    = re_item0
  theta_time_flag = theta_time0
  Q               = Q0

  p_fix = cfa_p_fix
  p_ran = 0

  dimEta = n * Q
  dim_eta = dimEta

  IF (agh_fit_id_cached /= fit_id0) THEN
   agh_grid_ready = .FALSE.
   agh_fit_id_cached = fit_id0
   END IF

  IF (phi_struct0 == 2) THEN
    CALL rexit("ERROR: phi_struct_flag=2 (ricorsione) non supportato con AGH")
  END IF

  ! dati (cleanup)
  IF (ALLOCATED(nmes))    DEALLOCATE(nmes)
  IF (ALLOCATED(Y))       DEALLOCATE(Y)
  IF (ALLOCATED(nitem))   DEALLOCATE(nitem)
  IF (ALLOCATED(c))       DEALLOCATE(c)
  IF (ALLOCATED(cm1))     DEALLOCATE(cm1)

  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs), c(nvar), cm1(nvar))
  nmes  = nmes0
  Y     = Y0
  nitem = nitem0
  c     = c0
  cm1   = cm10

IF (nobs <= 0) THEN
   CALL rexit('ERROR: nobs <= 0')
END IF

IF (SUM(nmes) <= 0) THEN
   CALL rexit('ERROR: SUM(nmes) <= 0. No observations are being used')
END IF

IF (MAXVAL(Y) <= 0) THEN
   CALL rexit('ERROR: Y seems not correctly coded. MAXVAL(Y)')
END IF

IF (MINVAL(nitem) < 1 .OR. MAXVAL(nitem) > n) THEN
   CALL rexit('ERROR: nitem out of range. MIN')
END IF

  IF (ALLOCATED(loadmask)) DEALLOCATE(loadmask)
  ALLOCATE(loadmask(nvar,Q))
  loadmask = loadmask0

  IF (re_item_flag == 1) THEN
     dimZ = dimEta + nvar
  ELSE
     dimZ = dimEta
  END IF

  IF (ALLOCATED(nodes))   DEALLOCATE(nodes)
  IF (ALLOCATED(weights)) DEALLOCATE(weights)
  ALLOCATE(nodes(nquad*dimZ))
  nodes(1:nquad*dimZ) = nodes0(1:nquad*dimZ)
  ALLOCATE(weights(nquad))
  weights(1:nquad) = weights0(1:nquad)

  ! parametri (cleanup)
  IF (ALLOCATED(lambda))       DEALLOCATE(lambda)
  IF (ALLOCATED(gamma))        DEALLOCATE(gamma)
  IF (ALLOCATED(theta_k))      DEALLOCATE(theta_k)
  IF (ALLOCATED(invsd_k))      DEALLOCATE(invsd_k)
  IF (ALLOCATED(sigmak))       DEALLOCATE(sigmak)
  IF (ALLOCATED(sqrt_sigmak))  DEALLOCATE(sqrt_sigmak)
  IF (ALLOCATED(kappa))        DEALLOCATE(kappa)
  IF (ALLOCATED(L_phi))        DEALLOCATE(L_phi)
  IF (ALLOCATED(tau))          DEALLOCATE(tau)
  IF (ALLOCATED(theta_long))   DEALLOCATE(theta_long)

  ALLOCATE(lambda(nvar,Q))
  ALLOCATE(gamma(nvar))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(sigmak(nvar),  sqrt_sigmak(nvar))
  ALLOCATE(kappa(dimEta))
  ALLOCATE(L_phi(dimEta,dimEta))
  ALLOCATE(theta_long(nvar,n))

  IF (ALLOCATED(tau_long)) DEALLOCATE(tau_long)
  IF (ALLOCATED(gamma_long)) DEALLOCATE(gamma_long)
  IF (ALLOCATED(lambda_long)) DEALLOCATE(lambda_long)
  ALLOCATE(cm3(nvar), ta(MAX(1, ntre0)), tau(SUM(cm1)))

  ! thresholds: identification controlled by cfa_threshold_id_flag
  DO k = 1, nvar
     cm3(k) = cm1(k) - 2
     IF (cm1(k) < 3) THEN
        CALL rexit('ERROR (CFA): threshold identification requires at least 4 response categories for every item')
     END IF
  END DO

  ntre_calc = SUM(cm3)
  IF (ntre_calc /= ntre0) THEN
     CALL rexit('ERROR (CFA): ntre0 mismatch. ntre0')
  END IF
  ntre = ntre0
  IF (ntre0 > 0) ta(1:ntre0) = b0(1:ntre0)

  pos = 0
  j = 1
  tau = 0.d0
  DO k = 1, nvar
     CALL cfa_fill_tau_item(cm1(k), ta(j), monotone_tau0, tau(pos+1:pos+cm1(k)))
     j = j + cm3(k)
     pos = pos + cm1(k)
  END DO

  ! lambda
  ! lambda
  lambda = 0.d0
  pos = ntre0

  DO qq = 1, Q
     ref_item = 0
     DO k = 1, nvar
        IF (loadmask(k,qq) == 1) THEN
           ref_item = k
           EXIT
        END IF
     END DO
     IF (ALLOCATED(cfa_marker_item)) ref_item = cfa_marker_item(qq)
     IF (ref_item == 0 .OR. loadmask(ref_item,qq) /= 1) THEN
        CALL rexit('ERROR (loglik_long_agh): factor')
     END IF

     lambda(ref_item,qq) = 1.d0

     DO k = 1, nvar
        IF (loadmask(k,qq) == 0) THEN
           lambda(k,qq) = 0.d0
        ELSE
           IF (k == ref_item) CYCLE
           pos = pos + 1
           lambda(k,qq) = b0(pos)
        END IF
     END DO
  END DO

  ! gamma
  DO k = 1, nvar
     pos = pos + 1
     gamma(k) = b0(pos)
  END DO

  ! theta
  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) == 1) THEN
           theta_k(k) = 1.d0
           invsd_k(k) = 1.d0
           DO t = 1, n
              theta_long(k,t) = 1.d0
           END DO
        ELSE
           pos = pos + 1
           theta_k(k) = EXP(b0(pos))
           invsd_k(k) = 1.d0 / SQRT(theta_k(k))
           DO t = 1, n
              theta_long(k,t) = theta_k(k)
           END DO
        END IF
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) == 1) THEN
           theta_k(k) = 1.d0
           invsd_k(k) = 1.d0
           DO t = 1, n
              theta_long(k,t) = 1.d0
           END DO
        ELSE
           DO t = 1, n
              pos = pos + 1
              theta_long(k,t) = EXP(b0(pos))
           END DO
           theta_k(k) = theta_long(k,1)
           invsd_k(k) = 1.d0 / SQRT(theta_k(k))
        END IF
     END DO
  END IF

  ! sigmak
  IF (re_item_flag == 1) THEN
     DO k = 1, nvar
        pos = pos + 1
        sigmak(k)      = EXP(b0(pos))
        sqrt_sigmak(k) = SQRT(sigmak(k))
     END DO
  ELSE
     DO k = 1, nvar
        sigmak(k)      = 0.d0
        sqrt_sigmak(k) = 0.d0
     END DO
  END IF

  ! kappa
  DO qq = 1, Q
     kappa( (1-1)*Q + qq ) = 0.d0
  END DO
  IF (n > 1) THEN
     DO t = 2, n
        DO qq = 1, Q
           pos = pos + 1
           kappa( (t-1)*Q + qq ) = b0(pos)
        END DO
     END DO
  END IF

  ! factor regressions on covariates with time-invariant beta
  CALL ensure_cfa_beta_storage()
  IF (p_fix > 0) THEN
     DO qq = 1, Q
        DO r = 1, p_fix
           pos = pos + 1
           cfa_beta(qq,r) = b0(pos)
        END DO
     END DO
  END IF

  ! L_phi
  L_phi = 0.d0

  IF (phi_struct_flag == 0) THEN
     ! L_phi libero
     DO t = 1, dim_eta
        pos = pos + 1
        L_phi(t,t) = EXP(b0(pos))
     END DO
     DO col = 1, dim_eta
        DO row = col+1, dim_eta
           pos = pos + 1
           L_phi(row,col) = b0(pos)
        END DO
     END DO

  ELSE
     ! AR1: A piena LIBERA, Sigma1 diag, costruisco L_phi
     ALLOCATE(A_ar(Q,Q), sqrt_sig1(Q))

     DO col = 1, Q
       DO row = 1, Q
         pos = pos + 1
         A_ar(row,col) = b0(pos)
       END DO
     END DO

     DO row = 1, Q
       pos = pos + 1
       sqrt_sig1(row) = SQRT(EXP(b0(pos)))
     END DO

     CALL build_Lphi_ar1_full(n, Q, A_ar, sqrt_sig1, L_phi)

     DEALLOCATE(A_ar, sqrt_sig1)
  END IF

  ! expected npm warning
  n_theta_est = 0
  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) n_theta_est = n_theta_est + 1
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) n_theta_est = n_theta_est + n
     END DO
  END IF

  n_lambda_free = 0
  DO qq = 1, Q
     DO k = 1, nvar
        IF (loadmask(k,qq) == 1) n_lambda_free = n_lambda_free + 1
     END DO
  END DO
  n_lambda_free = n_lambda_free - Q

  IF (phi_struct_flag == 0) THEN
     n_phi_par = dim_eta*(dim_eta+1)/2
  ELSE
     n_phi_par = Q*Q + Q
  END IF

  expected_npm = ntre0                            &
              + n_lambda_free                     &
              + nvar                              &
              + n_theta_est                       &
              + MERGE(nvar, 0, re_item_flag == 1) &
              + (n-1)*Q                           &
              + Q*p_fix                          &
              + n_phi_par

  IF (npm0 /= expected_npm) THEN
     CALL rexit("The supplied parameter vector has the wrong length for this CFA model.")
  END IF

  ! ! -loglik AGH
  ! vrais = -vrais_long_long_agh()
  ! -loglik AGH fast, stile GLMMadaptive

 IF (update_gh0 == 1 .OR. .NOT. agh_grid_ready) THEN

   CALL build_agh_grid_long(ok_grid)

  IF (.NOT. ok_grid .OR. .NOT. agh_grid_ready) THEN
    vrais = BAD_OBJ
    IF (ALLOCATED(cm3)) DEALLOCATE(cm3)
    IF (ALLOCATED(ta)) DEALLOCATE(ta)
    RETURN
  END IF

  END IF

CALL vrais_long_long_agh_fixed(ll_agh)

IF (ll_agh /= ll_agh .OR. ABS(ll_agh) > 1.d300) THEN
   vrais = BAD_OBJ
ELSE
   vrais = -ll_agh
END IF


  ! print *, vrais

  DEALLOCATE(cm3, ta)

END SUBROUTINE loglik_long_agh















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                           !!!
!!!         AGH nested        !!!
!!!                           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!===========================================================
!======================  AGH PART  =========================
!===========================================================


SUBROUTINE build_agh_grid_eta_nested(ok_grid)
  USE shared, ONLY: nsub, nvar, n, Q, nquad, nodes, nmes,              &
                    agh_tt_grid_ready, agh_tt_dim_cached,              &
                    agh_tt_nquad_cached, agh_tt_nsub_cached,           &
                    agh_tt_z, agh_tt_logconst, agh_tt_modes
  USE datindiv, ONLY: nmescur
  IMPLICIT NONE

  LOGICAL, INTENT(OUT) :: ok_grid

  INTEGER :: i, sm, l, d, dimEta, baseZ
  DOUBLE PRECISION :: logdetR
  LOGICAL :: ok_mode

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: z_mode, z, zeta, z_start
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R

  ok_grid = .TRUE.
  dimEta = n * Q

  IF (ALLOCATED(agh_tt_z)) THEN
     IF (agh_tt_dim_cached /= dimEta .OR. agh_tt_nquad_cached /= nquad .OR. &
         agh_tt_nsub_cached /= nsub) THEN
        DEALLOCATE(agh_tt_z)
        DEALLOCATE(agh_tt_logconst)
        DEALLOCATE(agh_tt_modes)
        agh_tt_grid_ready = .FALSE.
     END IF
  END IF

  IF (.NOT. ALLOCATED(agh_tt_z)) THEN
     ALLOCATE(agh_tt_z(dimEta, nquad, nsub))
     ALLOCATE(agh_tt_logconst(nsub))
     ALLOCATE(agh_tt_modes(dimEta, nsub))

     agh_tt_dim_cached = dimEta
     agh_tt_nquad_cached = nquad
     agh_tt_nsub_cached = nsub
     agh_tt_grid_ready = .FALSE.
  END IF

  ALLOCATE(z_mode(dimEta), z(dimEta), zeta(dimEta), z_start(dimEta), R(dimEta,dimEta))

  nmescur = 0

  DO i = 1, nsub

     CALL activate_growth_subject(i)
     sm = SUM(nmes(i,:))

     IF (agh_tt_grid_ready) THEN
        z_start = agh_tt_modes(:,i)
     ELSE
        z_start = 0.d0
     END IF

     CALL find_mode_and_R_eta_nested(i, dimEta, z_start, z_mode, R, logdetR, ok_mode)

     IF (.NOT. ok_mode) THEN
        ok_grid = .FALSE.
        agh_tt_grid_ready = .FALSE.
        DEALLOCATE(z_mode, z, zeta, z_start, R)
        RETURN
     END IF

     agh_tt_logconst(i) = logdetR
     agh_tt_modes(:,i) = z_mode

     DO l = 1, nquad
        baseZ = (l - 1) * dimEta

        DO d = 1, dimEta
           zeta(d) = nodes(baseZ + d)
        END DO

        z = z_mode
        CALL DGEMV('N', dimEta, dimEta, 1.d0, R, dimEta, zeta, 1, 1.d0, z, 1)

        agh_tt_z(:,l,i) = z
     END DO

     nmescur = nmescur + sm

  END DO

  agh_tt_grid_ready = .TRUE.

  DEALLOCATE(z_mode, z, zeta, z_start, R)

END SUBROUTINE build_agh_grid_eta_nested


SUBROUTINE item_logh_derivs_a_cond_eta(i, k_item, eta_vec, a, logh, score_a, hess_a, ok)
  USE shared, ONLY: nvar, n, Q, nmes, c, Y, nitem, theta_long, lambda, gamma, tau, loadmask
  USE datindiv, ONLY: nmescur
  USE utilities, ONLY: alnorm_fast
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, k_item
  DOUBLE PRECISION, DIMENSION(n*Q), INTENT(IN) :: eta_vec
  DOUBLE PRECISION, INTENT(IN) :: a
  DOUBLE PRECISION, INTENT(OUT) :: logh, score_a, hess_a
  LOGICAL, INTENT(OUT) :: ok

  DOUBLE PRECISION, PARAMETER :: INV_SQRT_2PI = 0.39894228040143267794d0
  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9

  INTEGER :: kk_thr, idx_obs, kk, j_mes, qq, t_idx, idx_tq
  INTEGER :: yijk, tuk
  DOUBLE PRECISION :: mu_ijk, linpred
  DOUBLE PRECISION :: var_lat_k, invsd_lat_k
  DOUBLE PRECISION :: alpha, beta_thr
  DOUBLE PRECISION :: p, phi_a, phi_b
  DOUBLE PRECISION :: dP_dmu, d2P_dmu2
  DOUBLE PRECISION :: s_mu, h_mu

  ok = .TRUE.
  logh = 0.d0
  score_a = 0.d0
  hess_a = 0.d0

  kk_thr = 0
  idx_obs = nmescur

  IF (k_item > 1) THEN
     DO kk = 1, k_item - 1
        kk_thr = kk_thr + (c(kk) - 1)
        idx_obs = idx_obs + nmes(i,kk)
     END DO
  END IF

  IF (nmes(i,k_item) == 0) RETURN

  DO j_mes = 1, nmes(i,k_item)

     idx_obs = idx_obs + 1
     t_idx   = nitem(idx_obs)

     var_lat_k = theta_at(k_item, t_idx)

     IF (var_lat_k <= 0.d0 .OR. var_lat_k /= var_lat_k) THEN
        ok = .FALSE.
        logh = LOG_ZERO
        RETURN
     END IF

     invsd_lat_k = 1.d0 / SQRT(var_lat_k)

     linpred = 0.d0
     DO qq = 1, Q
        IF (loadmask(k_item,qq) == 0) CYCLE
        idx_tq = (t_idx - 1) * Q + qq
        linpred = linpred + lambda_at(k_item,qq,t_idx) * eta_vec(idx_tq)
     END DO

     mu_ijk = gamma_at(k_item,t_idx) + a + linpred

     yijk = Y(idx_obs)
     tuk  = yijk + kk_thr

     IF (c(k_item) == 2) THEN

        alpha = (tau_at(kk_thr+1,t_idx) - mu_ijk) * invsd_lat_k
        phi_a = INV_SQRT_2PI * EXP(-0.5d0 * alpha * alpha)

        IF (yijk == 1) THEN
           p = alnorm_fast(alpha, .FALSE.)
           dP_dmu   = -phi_a * invsd_lat_k
           d2P_dmu2 = -alpha * phi_a * invsd_lat_k * invsd_lat_k
        ELSE
           p = alnorm_fast(alpha, .TRUE.)
           dP_dmu   =  phi_a * invsd_lat_k
           d2P_dmu2 =  alpha * phi_a * invsd_lat_k * invsd_lat_k
        END IF

     ELSE

        IF (yijk == 1) THEN

           alpha = (tau_at(tuk,t_idx) - mu_ijk) * invsd_lat_k
           p = alnorm_fast(alpha, .FALSE.)

           phi_a    = INV_SQRT_2PI * EXP(-0.5d0 * alpha * alpha)
           dP_dmu   = -phi_a * invsd_lat_k
           d2P_dmu2 = -alpha * phi_a * invsd_lat_k * invsd_lat_k

        ELSE IF (yijk == c(k_item)) THEN

           beta_thr = (tau_at(tuk-1,t_idx) - mu_ijk) * invsd_lat_k
           p = alnorm_fast(beta_thr, .TRUE.)

           phi_b    = INV_SQRT_2PI * EXP(-0.5d0 * beta_thr * beta_thr)
           dP_dmu   =  phi_b * invsd_lat_k
           d2P_dmu2 =  beta_thr * phi_b * invsd_lat_k * invsd_lat_k

        ELSE

           alpha    = (tau_at(tuk,t_idx)   - mu_ijk) * invsd_lat_k
           beta_thr = (tau_at(tuk-1,t_idx) - mu_ijk) * invsd_lat_k

           IF (beta_thr > 0.d0) THEN
                  p = alnorm_fast(beta_thr, .TRUE.) - alnorm_fast(alpha, .TRUE.)
               ELSE
                  p = alnorm_fast(alpha, .FALSE.) - alnorm_fast(beta_thr, .FALSE.)
               END IF

           phi_a    = INV_SQRT_2PI * EXP(-0.5d0 * alpha * alpha)
           phi_b    = INV_SQRT_2PI * EXP(-0.5d0 * beta_thr * beta_thr)

           dP_dmu   = (phi_b - phi_a) * invsd_lat_k
           d2P_dmu2 = (beta_thr * phi_b - alpha * phi_a) * invsd_lat_k * invsd_lat_k

        END IF
     END IF

     IF (p <= 0.d0 .OR. p /= p) THEN
        ok = .FALSE.
        logh = LOG_ZERO
        RETURN
     END IF

     s_mu = dP_dmu / p
     h_mu = d2P_dmu2 / p - s_mu * s_mu

     IF (s_mu /= s_mu .OR. h_mu /= h_mu) THEN
        ok = .FALSE.
        logh = LOG_ZERO
        RETURN
     END IF

     logh = logh + LOG(p)
     score_a = score_a + s_mu
     hess_a  = hess_a  + h_mu

  END DO

END SUBROUTINE item_logh_derivs_a_cond_eta


! Likelihood-only version used inside fixed-grid quadrature.
! It deliberately avoids normal-density, score and Hessian calculations.
SUBROUTINE item_logh_only_a_cond_eta(i, k_item, eta_vec, a, logh, ok)
  USE shared, ONLY: n, Q, nmes, c, Y, nitem, loadmask
  USE datindiv, ONLY: nmescur
  USE utilities, ONLY: alnorm_fast
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, k_item
  DOUBLE PRECISION, DIMENSION(n*Q), INTENT(IN) :: eta_vec
  DOUBLE PRECISION, INTENT(IN) :: a
  DOUBLE PRECISION, INTENT(OUT) :: logh
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: kk_thr, idx_obs, kk, j_mes, qq, t_idx, idx_tq
  INTEGER :: yijk, tuk
  DOUBLE PRECISION :: mu_ijk, linpred, var_lat_k, invsd_lat_k
  DOUBLE PRECISION :: alpha, beta_thr, p
  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9

  ok = .TRUE.
  logh = 0.d0
  kk_thr = 0
  idx_obs = nmescur

  IF (k_item > 1) THEN
     DO kk = 1, k_item - 1
        kk_thr = kk_thr + (c(kk) - 1)
        idx_obs = idx_obs + nmes(i,kk)
     END DO
  END IF

  IF (nmes(i,k_item) == 0) RETURN

  DO j_mes = 1, nmes(i,k_item)
     idx_obs = idx_obs + 1
     t_idx = nitem(idx_obs)

     var_lat_k = theta_at(k_item, t_idx)
     IF (var_lat_k <= 0.d0 .OR. var_lat_k /= var_lat_k) THEN
        ok = .FALSE.
        logh = LOG_ZERO
        RETURN
     END IF
     invsd_lat_k = 1.d0 / SQRT(var_lat_k)

     linpred = 0.d0
     DO qq = 1, Q
        IF (loadmask(k_item,qq) == 0) CYCLE
        idx_tq = (t_idx - 1) * Q + qq
        linpred = linpred + lambda_at(k_item,qq,t_idx) * eta_vec(idx_tq)
     END DO
     mu_ijk = gamma_at(k_item,t_idx) + a + linpred

     yijk = Y(idx_obs)
     tuk = yijk + kk_thr

     IF (c(k_item) == 2) THEN
        alpha = (tau_at(kk_thr+1,t_idx) - mu_ijk) * invsd_lat_k
        IF (yijk == 1) THEN
           p = alnorm_fast(alpha, .FALSE.)
        ELSE
           p = alnorm_fast(alpha, .TRUE.)
        END IF
     ELSE IF (yijk == 1) THEN
        alpha = (tau_at(tuk,t_idx) - mu_ijk) * invsd_lat_k
        p = alnorm_fast(alpha, .FALSE.)
     ELSE IF (yijk == c(k_item)) THEN
        beta_thr = (tau_at(tuk-1,t_idx) - mu_ijk) * invsd_lat_k
        p = alnorm_fast(beta_thr, .TRUE.)
     ELSE
        alpha = (tau_at(tuk,t_idx) - mu_ijk) * invsd_lat_k
        beta_thr = (tau_at(tuk-1,t_idx) - mu_ijk) * invsd_lat_k
        IF (beta_thr > 0.d0) THEN
                  p = alnorm_fast(beta_thr, .TRUE.) - alnorm_fast(alpha, .TRUE.)
               ELSE
                  p = alnorm_fast(alpha, .FALSE.) - alnorm_fast(beta_thr, .FALSE.)
               END IF
     END IF

     IF (p <= 0.d0 .OR. p /= p) THEN
        ok = .FALSE.
        logh = LOG_ZERO
        RETURN
     END IF
     logh = logh + LOG(p)
  END DO
END SUBROUTINE item_logh_only_a_cond_eta


SUBROUTINE item_neglogpost_a_cond_eta(i, k_item, eta_vec, a, f, grad, curv, ok)
  USE shared, ONLY: n, Q, sigmak
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, k_item
  DOUBLE PRECISION, DIMENSION(n*Q), INTENT(IN) :: eta_vec
  DOUBLE PRECISION, INTENT(IN) :: a
  DOUBLE PRECISION, INTENT(OUT) :: f, grad, curv
  LOGICAL, INTENT(OUT) :: ok

  DOUBLE PRECISION :: logh, score_a, hess_a, psi
  DOUBLE PRECISION, PARAMETER :: BAD = 1.d12

  psi = sigmak(k_item)

  IF (psi <= 0.d0 .OR. psi /= psi) THEN
     ok = .FALSE.
     f = BAD
     grad = 0.d0
     curv = 1.d0
     RETURN
  END IF

  CALL item_logh_derivs_a_cond_eta(i, k_item, eta_vec, a, logh, score_a, hess_a, ok)

  IF (.NOT. ok .OR. logh /= logh) THEN
     f = BAD
     grad = 0.d0
     curv = 1.d0
     ok = .FALSE.
     RETURN
  END IF

  f    = -logh + 0.5d0 * a * a / psi
  grad = -score_a + a / psi
  curv = -hess_a + 1.d0 / psi

  IF (f /= f .OR. grad /= grad .OR. curv /= curv) THEN
     ok = .FALSE.
     f = BAD
     grad = 0.d0
     curv = 1.d0
  END IF

END SUBROUTINE item_neglogpost_a_cond_eta


SUBROUTINE inner_mode_a_item_cond_eta(i, k_item, eta_vec, alpha, sdev, ok_mode)
  USE shared, ONLY: n, Q, sigmak
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, k_item
  DOUBLE PRECISION, DIMENSION(n*Q), INTENT(IN) :: eta_vec
  DOUBLE PRECISION, INTENT(OUT) :: alpha, sdev
  LOGICAL, INTENT(OUT) :: ok_mode

  INTEGER :: iter, ls
  DOUBLE PRECISION :: a, a_try, f, f_try, grad, curv, curv_try
  DOUBLE PRECISION :: step, fac, psi
  LOGICAL :: ok, ok_try

  ok_mode = .TRUE.
  psi = sigmak(k_item)

  IF (psi <= 0.d0 .OR. psi /= psi) THEN
     alpha = 0.d0
     sdev = 0.d0
     ok_mode = .FALSE.
     RETURN
  END IF

  a = 0.d0

  DO iter = 1, 35

     CALL item_neglogpost_a_cond_eta(i, k_item, eta_vec, a, f, grad, curv, ok)

     IF (.NOT. ok) THEN
        ok_mode = .FALSE.
        alpha = 0.d0
        sdev = SQRT(psi)
        RETURN
     END IF

     IF (curv <= 1.d-8) curv = 1.d-8

     step = -grad / curv

     IF (ABS(step) > 5.d0 * SQRT(psi)) THEN
        step = SIGN(5.d0 * SQRT(psi), step)
     END IF

     fac = 1.d0

     DO ls = 1, 20
        a_try = a + fac * step

        CALL item_neglogpost_a_cond_eta(i, k_item, eta_vec, a_try, f_try, grad, curv_try, ok_try)

        IF (ok_try .AND. f_try < f) THEN
           a = a_try
           EXIT
        END IF

        fac = 0.5d0 * fac
     END DO

     IF (ABS(fac * step) < 1.d-7 * (ABS(a) + 1.d0)) EXIT

  END DO

  CALL item_neglogpost_a_cond_eta(i, k_item, eta_vec, a, f, grad, curv, ok)

  IF (.NOT. ok .OR. curv <= 0.d0 .OR. curv /= curv) THEN
     ok_mode = .FALSE.
     alpha = 0.d0
     sdev = SQRT(psi)
     RETURN
  END IF

  alpha = a
  sdev = 1.d0 / SQRT(curv)

END SUBROUTINE inner_mode_a_item_cond_eta

SUBROUTINE inneragh_item_logint_cond_eta(i, k_item, eta_vec, logint, ok)
  USE shared, ONLY: n, Q, nquad_re, nodes_re, weights_re, sigmak, re_item_flag, nmes
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, k_item
  DOUBLE PRECISION, DIMENSION(n*Q), INTENT(IN) :: eta_vec
  DOUBLE PRECISION, INTENT(OUT) :: logint
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: r, n_valid
  DOUBLE PRECISION :: psi, alpha, sdev, v, a
  DOUBLE PRECISION :: logh
  DOUBLE PRECISION :: term, max_term, sum_exp
  LOGICAL :: ok_h, ok_mode

  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9

  ok = .TRUE.
  logint = 0.d0

  IF (nmes(i,k_item) == 0) RETURN

  IF (re_item_flag == 0) THEN
     CALL item_logh_only_a_cond_eta(i, k_item, eta_vec, 0.d0, logh, ok_h)
     IF (.NOT. ok_h) THEN
        ok = .FALSE.
        logint = LOG_ZERO
     ELSE
        logint = logh
     END IF
     RETURN
  END IF

  psi = sigmak(k_item)

  IF (psi <= 0.d0 .OR. psi /= psi) THEN
     ok = .FALSE.
     logint = LOG_ZERO
     RETURN
  END IF

  CALL inner_mode_a_item_cond_eta(i, k_item, eta_vec, alpha, sdev, ok_mode)

  IF (.NOT. ok_mode .OR. sdev <= 0.d0 .OR. sdev /= sdev) THEN
     ok = .FALSE.
     logint = LOG_ZERO
     RETURN
  END IF

  n_valid = 0
  max_term = 0.d0
  sum_exp = 0.d0

  DO r = 1, nquad_re

     v = nodes_re(r)
     a = alpha + sdev * v

     CALL item_logh_only_a_cond_eta(i, k_item, eta_vec, a, logh, ok_h)

     IF (.NOT. ok_h) CYCLE
     IF (logh /= logh) CYCLE
     IF (weights_re(r) <= 0.d0) CYCLE

     term = LOG(weights_re(r)) + LOG(sdev) + logh              &
          - 0.5d0 * LOG(psi)                                  &
          - 0.5d0 * a * a / psi                               &
          + 0.5d0 * v * v

     IF (term /= term) CYCLE

     IF (n_valid == 0) THEN
        max_term = term
        sum_exp = 1.d0
        n_valid = 1
     ELSE IF (term > max_term) THEN
        sum_exp = sum_exp * EXP(max_term - term) + 1.d0
        max_term = term
        n_valid = n_valid + 1
     ELSE
        sum_exp = sum_exp + EXP(term - max_term)
        n_valid = n_valid + 1
     END IF

  END DO

  IF (n_valid == 0 .OR. sum_exp <= 0.d0 .OR. sum_exp /= sum_exp) THEN
     ok = .FALSE.
     logint = LOG_ZERO
     RETURN
  END IF

  logint = max_term + LOG(sum_exp)

END SUBROUTINE inneragh_item_logint_cond_eta

SUBROUTINE logg_eta_nested_i(i, dimEta, z_eta, logg, ok)
  USE shared, ONLY: nvar, n, Q, kappa, L_phi
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, dimEta
  DOUBLE PRECISION, DIMENSION(dimEta), INTENT(IN) :: z_eta
  DOUBLE PRECISION, INTENT(OUT) :: logg
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: s_idx, qq, k_item
  DOUBLE PRECISION :: logint_k
  LOGICAL :: ok_k
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eta_vec

  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9

  ok = .TRUE.
  logg = 0.d0

  IF (dimEta /= n * Q) THEN
     ok = .FALSE.
     logg = LOG_ZERO
     RETURN
  END IF

  ALLOCATE(eta_vec(dimEta))

  CALL cfa_eta_mean(i, eta_vec)
  DO s_idx = 1, dimEta
     DO qq = 1, s_idx
        eta_vec(s_idx) = eta_vec(s_idx) + L_phi(s_idx,qq) * z_eta(qq)
     END DO
  END DO

  DO k_item = 1, nvar
     CALL inneragh_item_logint_cond_eta(i, k_item, eta_vec, logint_k, ok_k)

     IF (.NOT. ok_k .OR. logint_k /= logint_k) THEN
        ok = .FALSE.
        logg = LOG_ZERO
        DEALLOCATE(eta_vec)
        RETURN
     END IF

     logg = logg + logint_k
  END DO

  DEALLOCATE(eta_vec)

END SUBROUTINE logg_eta_nested_i

SUBROUTINE fd_grad_outer_eta_nested_i(i, dimEta, z, f0, grad, ok)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, dimEta
  DOUBLE PRECISION, DIMENSION(dimEta), INTENT(IN) :: z
  DOUBLE PRECISION, INTENT(OUT) :: f0
  DOUBLE PRECISION, DIMENSION(dimEta), INTENT(OUT) :: grad
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: r
  DOUBLE PRECISION :: hr, fp, fm
  LOGICAL :: ok0, okp, okm
  DOUBLE PRECISION, PARAMETER :: eps = 1.d-5

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: zp, zm

  ok = .TRUE.
  grad = 0.d0

  ALLOCATE(zp(dimEta), zm(dimEta))

  CALL outer_obj_eta_nested_i(i, dimEta, z, f0, ok0)

  IF (.NOT. ok0 .OR. f0 /= f0 .OR. ABS(f0) > 1.d300) THEN
     ok = .FALSE.
     DEALLOCATE(zp, zm)
     RETURN
  END IF

  DO r = 1, dimEta
     hr = eps * MAX(ABS(z(r)), 1.d0)

     zp = z
     zm = z

     zp(r) = z(r) + hr
     zm(r) = z(r) - hr

     CALL outer_obj_eta_nested_i(i, dimEta, zp, fp, okp)
     CALL outer_obj_eta_nested_i(i, dimEta, zm, fm, okm)

     IF (.NOT. okp .OR. .NOT. okm .OR. fp /= fp .OR. fm /= fm) THEN
        ok = .FALSE.
        DEALLOCATE(zp, zm)
        RETURN
     END IF

     grad(r) = (fp - fm) / (2.d0 * hr)
  END DO

  DEALLOCATE(zp, zm)

END SUBROUTINE fd_grad_outer_eta_nested_i

SUBROUTINE fg_eta_nested_for_vmmin(i_c, sm_c, dimZ_c, z, f, g, ok_c) BIND(C, name="fg_eta_nested_for_vmmin")
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE

  INTEGER(C_INT), INTENT(IN) :: i_c, sm_c, dimZ_c
  REAL(C_DOUBLE), DIMENSION(dimZ_c), INTENT(IN) :: z
  REAL(C_DOUBLE), INTENT(OUT) :: f
  REAL(C_DOUBLE), DIMENSION(dimZ_c), INTENT(OUT) :: g
  INTEGER(C_INT), INTENT(OUT) :: ok_c

  INTEGER :: i, dimZ
  LOGICAL :: ok

  i = INT(i_c)
  dimZ = INT(dimZ_c)

  CALL fd_grad_outer_eta_nested_i(i, dimZ, z, f, g, ok)

  IF (ok) THEN
     ok_c = 1_C_INT
  ELSE
     f = 1.d12
     g = 0.d0
     ok_c = 0_C_INT
  END IF

END SUBROUTINE fg_eta_nested_for_vmmin

SUBROUTINE outer_obj_eta_nested_i(i, dimEta, z_eta, f, ok)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, dimEta
  DOUBLE PRECISION, DIMENSION(dimEta), INTENT(IN) :: z_eta
  DOUBLE PRECISION, INTENT(OUT) :: f
  LOGICAL, INTENT(OUT) :: ok

  DOUBLE PRECISION :: logg, qz
  DOUBLE PRECISION, PARAMETER :: BAD = 1.d12

  CALL logg_eta_nested_i(i, dimEta, z_eta, logg, ok)

  IF (.NOT. ok .OR. logg /= logg) THEN
     f = BAD
     ok = .FALSE.
     RETURN
  END IF

  qz = DOT_PRODUCT(z_eta, z_eta)
  f = -logg + 0.5d0 * qz

  IF (f /= f .OR. ABS(f) > 1.d300) THEN
     f = BAD
     ok = .FALSE.
  END IF

END SUBROUTINE outer_obj_eta_nested_i


SUBROUTINE fd_grad_hess_outer_eta_nested_i(i, dimEta, z, f0, grad, H, ok)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, dimEta
  DOUBLE PRECISION, DIMENSION(dimEta), INTENT(IN) :: z
  DOUBLE PRECISION, INTENT(OUT) :: f0
  DOUBLE PRECISION, DIMENSION(dimEta), INTENT(OUT) :: grad
  DOUBLE PRECISION, DIMENSION(dimEta,dimEta), INTENT(OUT) :: H
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: r, c
  DOUBLE PRECISION :: hr, hc
  DOUBLE PRECISION :: fp, fm, fpp, fpm, fmp, fmm
  LOGICAL :: ok0, okp, okm, okpp, okpm, okmp, okmm
  DOUBLE PRECISION, PARAMETER :: eps = 1.d-4

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: zp, zm, zpp, zpm, zmp, zmm

  ok = .TRUE.
  grad = 0.d0
  H = 0.d0

  ALLOCATE(zp(dimEta), zm(dimEta), zpp(dimEta), zpm(dimEta), zmp(dimEta), zmm(dimEta))

  CALL outer_obj_eta_nested_i(i, dimEta, z, f0, ok0)

  IF (.NOT. ok0) THEN
     ok = .FALSE.
     DEALLOCATE(zp, zm, zpp, zpm, zmp, zmm)
     RETURN
  END IF

  DO r = 1, dimEta
     hr = eps * MAX(ABS(z(r)), 1.d0)

     zp = z
     zm = z

     zp(r) = zp(r) + hr
     zm(r) = zm(r) - hr

     CALL outer_obj_eta_nested_i(i, dimEta, zp, fp, okp)
     CALL outer_obj_eta_nested_i(i, dimEta, zm, fm, okm)

     IF (.NOT. okp .OR. .NOT. okm) THEN
        ok = .FALSE.
        DEALLOCATE(zp, zm, zpp, zpm, zmp, zmm)
        RETURN
     END IF

     grad(r) = (fp - fm) / (2.d0 * hr)
     H(r,r)  = (fp - 2.d0 * f0 + fm) / (hr * hr)
  END DO

  DO r = 1, dimEta - 1
     hr = eps * MAX(ABS(z(r)), 1.d0)

     DO c = r + 1, dimEta
        hc = eps * MAX(ABS(z(c)), 1.d0)

        zpp = z
        zpm = z
        zmp = z
        zmm = z

        zpp(r) = zpp(r) + hr
        zpp(c) = zpp(c) + hc

        zpm(r) = zpm(r) + hr
        zpm(c) = zpm(c) - hc

        zmp(r) = zmp(r) - hr
        zmp(c) = zmp(c) + hc

        zmm(r) = zmm(r) - hr
        zmm(c) = zmm(c) - hc

        CALL outer_obj_eta_nested_i(i, dimEta, zpp, fpp, okpp)
        CALL outer_obj_eta_nested_i(i, dimEta, zpm, fpm, okpm)
        CALL outer_obj_eta_nested_i(i, dimEta, zmp, fmp, okmp)
        CALL outer_obj_eta_nested_i(i, dimEta, zmm, fmm, okmm)

        IF (.NOT. okpp .OR. .NOT. okpm .OR. .NOT. okmp .OR. .NOT. okmm) THEN
           ok = .FALSE.
           DEALLOCATE(zp, zm, zpp, zpm, zmp, zmm)
           RETURN
        END IF

        H(r,c) = (fpp - fpm - fmp + fmm) / (4.d0 * hr * hc)
        H(c,r) = H(r,c)

     END DO
  END DO

  H = 0.5d0 * (H + TRANSPOSE(H))

  DEALLOCATE(zp, zm, zpp, zpm, zmp, zmm)

END SUBROUTINE fd_grad_hess_outer_eta_nested_i

SUBROUTINE find_mode_and_R_eta_nested(i, dimEta, z_start, z_mode, R, logdetR, ok_mode)
  USE shared
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE vmmin_eta_nested_mode(i, sm, dimZ, x, fmin, fail, fncount, grcount, &
                                      maxit, abstol, reltol, trace) BIND(C, name="vmmin_eta_nested_mode")
       USE, INTRINSIC :: ISO_C_BINDING
       IMPLICIT NONE

       INTEGER(C_INT), INTENT(IN) :: i, sm, dimZ
       REAL(C_DOUBLE), DIMENSION(*), INTENT(INOUT) :: x
       REAL(C_DOUBLE), INTENT(OUT) :: fmin
       INTEGER(C_INT), INTENT(OUT) :: fail, fncount, grcount
       INTEGER(C_INT), INTENT(IN) :: maxit, trace
       REAL(C_DOUBLE), INTENT(IN) :: abstol, reltol
     END SUBROUTINE vmmin_eta_nested_mode
  END INTERFACE

  INTEGER, INTENT(IN) :: i, dimEta
  DOUBLE PRECISION, DIMENSION(dimEta), INTENT(IN) :: z_start
  DOUBLE PRECISION, DIMENSION(dimEta), INTENT(OUT) :: z_mode
  DOUBLE PRECISION, DIMENSION(dimEta,dimEta), INTENT(OUT) :: R
  DOUBLE PRECISION, INTENT(OUT) :: logdetR
  LOGICAL, INTENT(OUT) :: ok_mode

  INTEGER :: info, j, jj
  LOGICAL :: ok
  DOUBLE PRECISION :: obj
  DOUBLE PRECISION, DIMENSION(dimEta) :: g
  DOUBLE PRECISION, DIMENSION(dimEta,dimEta) :: A

  INTEGER(C_INT) :: i_c, sm_c, dimZ_c
  INTEGER(C_INT) :: fail_vmmin, fncount_vmmin, grcount_vmmin
  INTEGER(C_INT) :: maxit_vmmin, trace_vmmin
  REAL(C_DOUBLE) :: abstol_vmmin, reltol_vmmin
  REAL(C_DOUBLE) :: obj_vmmin

  ok_mode = .TRUE.
  z_mode = z_start

  R = 0.d0
  DO j = 1, dimEta
     R(j,j) = 1.d0
  END DO
  logdetR = 0.d0

  i_c = INT(i, C_INT)
  sm_c = 0_C_INT
  dimZ_c = INT(dimEta, C_INT)

  maxit_vmmin = 60_C_INT
  trace_vmmin = 0_C_INT
  abstol_vmmin = -HUGE(1.0D0)
  reltol_vmmin = SQRT(EPSILON(1.0D0))

  CALL vmmin_eta_nested_mode(i_c, sm_c, dimZ_c, z_mode, obj_vmmin, fail_vmmin, &
                             fncount_vmmin, grcount_vmmin, maxit_vmmin,        &
                             abstol_vmmin, reltol_vmmin, trace_vmmin)

  IF (fail_vmmin /= 0_C_INT) THEN
     CALL intpr("Adaptive mode search stopped early for subject", -1, (/ i /), 1)
  END IF

  IF (obj_vmmin /= obj_vmmin .OR. ABS(obj_vmmin) > 1.d300) THEN
     ok_mode = .FALSE.
     z_mode = z_start
     R = 0.d0
     DO j = 1, dimEta
        R(j,j) = 1.d0
     END DO
     logdetR = 0.d0
     RETURN
  END IF

  CALL fd_grad_hess_outer_eta_nested_i(i, dimEta, z_mode, obj, g, A, ok)

  IF (.NOT. ok .OR. obj /= obj .OR. ABS(obj) > 1.d300) THEN
     ok_mode = .FALSE.
     z_mode = z_start
     R = 0.d0
     DO j = 1, dimEta
        R(j,j) = 1.d0
     END DO
     logdetR = 0.d0
     RETURN
  END IF

  A = 0.5d0 * (A + TRANSPOSE(A))

  CALL DPOTRF('U', dimEta, A, dimEta, info)

  IF (info /= 0 .OR. .NOT. ok) THEN
     ok_mode = .FALSE.
     z_mode = z_start
     R = 0.d0
     DO j = 1, dimEta
        R(j,j) = 1.d0
     END DO
     logdetR = 0.d0
     RETURN
  END IF

  logdetR = 0.d0

  DO j = 1, dimEta
     IF (A(j,j) <= 0.d0 .OR. A(j,j) /= A(j,j) .OR. ABS(A(j,j)) > 1.d300) THEN
        ok_mode = .FALSE.
        z_mode = z_start
        R = 0.d0
        DO jj = 1, dimEta
           R(jj,jj) = 1.d0
        END DO
        logdetR = 0.d0
        RETURN
     END IF

     logdetR = logdetR - LOG(A(j,j))
  END DO

  R = A

  IF (dimEta > 1) THEN
     DO j = 1, dimEta - 1
        R(j+1:dimEta,j) = 0.d0
     END DO
  END IF

  CALL DTRTRI('U', 'N', dimEta, R, dimEta, info)

  IF (info /= 0) THEN
     ok_mode = .FALSE.
     z_mode = z_start
     R = 0.d0
     DO j = 1, dimEta
        R(j,j) = 1.d0
     END DO
     logdetR = 0.d0
     RETURN
  END IF

END SUBROUTINE find_mode_and_R_eta_nested

SUBROUTINE vrais_long_i_tt_agh_fixed(i, sm, out, ok)
  USE shared, ONLY: n, Q, nquad, nodes, weights, agh_tt_z, agh_tt_logconst
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, sm
  DOUBLE PRECISION, INTENT(OUT) :: out
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: l, d, dimEta, baseZ, n_valid
  DOUBLE PRECISION :: logg, qz, qzeta, term, max_term, sum_exp
  LOGICAL :: ok_node

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: z_eta, zeta

  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9

  ok = .TRUE.
  out = 0.d0

  dimEta = n * Q
  CALL activate_growth_subject(i)

  ALLOCATE(z_eta(dimEta), zeta(dimEta))

  n_valid = 0
  max_term = 0.d0
  sum_exp = 0.d0

  DO l = 1, nquad

     baseZ = (l - 1) * dimEta

     DO d = 1, dimEta
        zeta(d) = nodes(baseZ + d)
     END DO

     z_eta = agh_tt_z(:,l,i)

     CALL logg_eta_nested_i(i, dimEta, z_eta, logg, ok_node)

     IF (.NOT. ok_node) CYCLE
     IF (logg /= logg) CYCLE
     IF (weights(l) <= 0.d0) CYCLE

     qz = DOT_PRODUCT(z_eta, z_eta)
     qzeta = DOT_PRODUCT(zeta, zeta)

     term = logg - 0.5d0 * (qz - qzeta) + LOG(weights(l))

     IF (term /= term) CYCLE

     IF (n_valid == 0) THEN
        max_term = term
        sum_exp = 1.d0
        n_valid = 1
     ELSE IF (term > max_term) THEN
        sum_exp = sum_exp * EXP(max_term - term) + 1.d0
        max_term = term
        n_valid = n_valid + 1
     ELSE
        sum_exp = sum_exp + EXP(term - max_term)
        n_valid = n_valid + 1
     END IF

  END DO

  IF (n_valid == 0 .OR. sum_exp <= 0.d0 .OR. sum_exp /= sum_exp) THEN
     out = LOG_ZERO
     ok = .FALSE.
     DEALLOCATE(z_eta, zeta)
     RETURN
  END IF

  out = agh_tt_logconst(i) + max_term + LOG(sum_exp)

  IF (out /= out) THEN
     out = LOG_ZERO
     ok = .FALSE.
  END IF

  DEALLOCATE(z_eta, zeta)

END SUBROUTINE vrais_long_i_tt_agh_fixed


SUBROUTINE vrais_long_tt_agh_fixed(out)
  USE shared, ONLY: nsub, nmes, nvar, agh_tt_grid_ready
  USE datindiv, ONLY: nmescur
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(OUT) :: out

  INTEGER :: i, sm
  DOUBLE PRECISION :: temp
  LOGICAL :: ok_i

  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9

  IF (.NOT. agh_tt_grid_ready) THEN
     CALL rexit('ERROR (vrais_long_tt_agh_fixed): nested AGH grid not ready')
  END IF

  nmescur = 0
  out = 0.d0

  DO i = 1, nsub

     sm = SUM(nmes(i,:))

     CALL vrais_long_i_tt_agh_fixed(i, sm, temp, ok_i)

     IF (ok_i .AND. temp == temp) THEN
        out = out + temp
     ELSE
        out = out + LOG_ZERO
     END IF

     nmescur = nmescur + sm

  END DO

  IF (out /= out) out = LOG_ZERO

END SUBROUTINE vrais_long_tt_agh_fixed



!===========================================================
! Full joint AGH written in two-tier form
! Uses the joint mode and joint Hessian, then factors the
! adaptive Gaussian proposal into q_eta * q_a|eta.
!===========================================================
SUBROUTINE build_agh_grid_full_equiv_tt(ok_grid)
  USE shared
  USE datindiv, ONLY: nmescur
  IMPLICIT NONE

  LOGICAL, INTENT(OUT) :: ok_grid

  INTEGER :: i, sm, l, d, cc, k, dimEta, dimZ, baseZ
  DOUBLE PRECISION :: logdetB, logdetEta, sdev_k
  LOGICAL :: ok_mode

  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: z_mode, z_start
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: u_eta, z_eta_node
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B_root

  ok_grid = .TRUE.

  dimEta = n * Q
  IF (re_item_flag == 1) THEN
     dimZ = dimEta + nvar
  ELSE
     dimZ = dimEta
  END IF

  IF (ALLOCATED(agh_tt_fe_zeta_eta)) THEN
     IF (agh_tt_fe_dimEta_cached /= dimEta .OR.                         &
         agh_tt_fe_nvar_cached /= nvar .OR.                              &
         agh_tt_fe_nquad_eta_cached /= nquad .OR.                        &
         agh_tt_fe_nquad_re_cached /= nquad_re .OR.                      &
         agh_tt_fe_nsub_cached /= nsub .OR.                              &
         SIZE(agh_tt_fe_modes,1) /= dimZ) THEN
        DEALLOCATE(agh_tt_fe_zeta_eta)
        DEALLOCATE(agh_tt_fe_alpha)
        DEALLOCATE(agh_tt_fe_sdev)
        DEALLOCATE(agh_tt_fe_logconst)
        DEALLOCATE(agh_tt_fe_modes)
        agh_tt_fe_grid_ready = .FALSE.
     END IF
  END IF

  IF (.NOT. ALLOCATED(agh_tt_fe_zeta_eta)) THEN
     ALLOCATE(agh_tt_fe_zeta_eta(dimEta, nquad, nsub))
     ALLOCATE(agh_tt_fe_alpha(nvar, nquad, nsub))
     ALLOCATE(agh_tt_fe_sdev(nvar, nsub))
     ALLOCATE(agh_tt_fe_logconst(nsub))
     ALLOCATE(agh_tt_fe_modes(dimZ, nsub))

     agh_tt_fe_dimEta_cached = dimEta
     agh_tt_fe_nvar_cached = nvar
     agh_tt_fe_nquad_eta_cached = nquad
     agh_tt_fe_nquad_re_cached = nquad_re
     agh_tt_fe_nsub_cached = nsub
  END IF

  ALLOCATE(z_mode(dimZ), z_start(dimZ), u_eta(dimEta), z_eta_node(dimEta))
  ALLOCATE(B_root(dimZ,dimZ))

  nmescur = 0

  DO i = 1, nsub

     CALL activate_growth_subject(i)
     sm = SUM(nmes(i,:))

     IF (agh_tt_fe_grid_ready) THEN
        z_start = agh_tt_fe_modes(:, i)
     ELSE
        z_start = 0.d0
     END IF

     CALL find_mode_and_Btt_long(i, sm, dimZ, z_start, z_mode, B_root,       &
          logdetB, logdetEta, ok_mode)

     IF (.NOT. ok_mode) THEN
        ok_grid = .FALSE.
        agh_tt_fe_grid_ready = .FALSE.
        DEALLOCATE(z_mode, z_start, u_eta, z_eta_node, B_root)
        RETURN
     END IF

     agh_tt_fe_logconst(i) = logdetEta
     agh_tt_fe_modes(:,i) = z_mode

     IF (re_item_flag == 1) THEN
        DO k = 1, nvar
           sdev_k = B_root(dimEta+k,dimEta+k)
           IF (sdev_k <= 0.d0 .OR. sdev_k /= sdev_k) THEN
              ok_grid = .FALSE.
              agh_tt_fe_grid_ready = .FALSE.
              DEALLOCATE(z_mode, z_start, u_eta, z_eta_node, B_root)
              RETURN
           END IF
           agh_tt_fe_sdev(k,i) = sdev_k
        END DO
     ELSE
        agh_tt_fe_sdev(:,i) = 1.d0
     END IF

     DO l = 1, nquad

        baseZ = (l - 1) * dimEta
        DO d = 1, dimEta
           u_eta(d) = nodes(baseZ + d)
        END DO

        DO d = 1, dimEta
           z_eta_node(d) = z_mode(d)
           DO cc = 1, d
              z_eta_node(d) = z_eta_node(d) + B_root(d,cc) * u_eta(cc)
           END DO
        END DO

        agh_tt_fe_zeta_eta(:,l,i) = z_eta_node

        IF (re_item_flag == 1) THEN
           DO k = 1, nvar
              agh_tt_fe_alpha(k,l,i) = z_mode(dimEta+k)
              DO d = 1, dimEta
                 agh_tt_fe_alpha(k,l,i) = agh_tt_fe_alpha(k,l,i)             &
                      + B_root(dimEta+k,d) * u_eta(d)
              END DO
           END DO
        ELSE
           agh_tt_fe_alpha(:,l,i) = 0.d0
        END IF

     END DO

     nmescur = nmescur + sm

  END DO

  agh_tt_fe_grid_ready = .TRUE.

  DEALLOCATE(z_mode, z_start, u_eta, z_eta_node, B_root)

END SUBROUTINE build_agh_grid_full_equiv_tt



SUBROUTINE inneragh_item_logint_full_equiv(i, k_item, eta_vec, alpha_z, sdev_z, logint, ok)
  USE shared, ONLY: n, Q, nquad_re, nodes_re, weights_re, sqrt_sigmak, re_item_flag, nmes
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, k_item
  DOUBLE PRECISION, DIMENSION(n*Q), INTENT(IN) :: eta_vec
  DOUBLE PRECISION, INTENT(IN) :: alpha_z, sdev_z
  DOUBLE PRECISION, INTENT(OUT) :: logint
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: r, n_valid
  DOUBLE PRECISION :: v, z_a, a, logh
  DOUBLE PRECISION :: term, max_term, sum_exp
  LOGICAL :: ok_h
  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9

  ok = .TRUE.
  logint = 0.d0

  IF (nmes(i,k_item) == 0) RETURN

  IF (re_item_flag == 0) THEN
     CALL item_logh_only_a_cond_eta(i, k_item, eta_vec, 0.d0, logh, ok_h)
     IF (.NOT. ok_h) THEN
        ok = .FALSE.
        logint = LOG_ZERO
     ELSE
        logint = logh
     END IF
     RETURN
  END IF

  IF (sdev_z <= 0.d0 .OR. sdev_z /= sdev_z) THEN
     ok = .FALSE.
     logint = LOG_ZERO
     RETURN
  END IF

  n_valid = 0
  max_term = 0.d0
  sum_exp = 0.d0

  DO r = 1, nquad_re

     IF (weights_re(r) <= 0.d0) CYCLE

     v = nodes_re(r)
     z_a = alpha_z + sdev_z * v
     a = sqrt_sigmak(k_item) * z_a

     CALL item_logh_only_a_cond_eta(i, k_item, eta_vec, a, logh, ok_h)

     IF (.NOT. ok_h) CYCLE
     IF (logh /= logh) CYCLE

     term = LOG(weights_re(r)) + LOG(sdev_z) + logh - 0.5d0 * (z_a*z_a - v*v)

     IF (term /= term) CYCLE

     IF (n_valid == 0) THEN
        max_term = term
        sum_exp = 1.d0
        n_valid = 1
     ELSE IF (term > max_term) THEN
        sum_exp = sum_exp * EXP(max_term - term) + 1.d0
        max_term = term
        n_valid = n_valid + 1
     ELSE
        sum_exp = sum_exp + EXP(term - max_term)
        n_valid = n_valid + 1
     END IF

  END DO

  IF (n_valid == 0 .OR. sum_exp <= 0.d0 .OR. sum_exp /= sum_exp) THEN
     ok = .FALSE.
     logint = LOG_ZERO
     RETURN
  END IF

  logint = max_term + LOG(sum_exp)

END SUBROUTINE inneragh_item_logint_full_equiv


SUBROUTINE vrais_long_i_tt_agh_full_equiv_fixed(i, sm, out, ok)
  USE shared, ONLY: nvar, n, Q, nquad, nodes, weights, re_item_flag,       &
                    kappa, L_phi, agh_tt_fe_zeta_eta, agh_tt_fe_alpha,    &
                    agh_tt_fe_sdev, agh_tt_fe_logconst
  USE datindiv, ONLY: nmescur
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, sm
  DOUBLE PRECISION, INTENT(OUT) :: out
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: l, d, r, k_item, dimEta, baseZ, n_valid
  DOUBLE PRECISION :: q_eta, q_u, logg, logint
  DOUBLE PRECISION :: term, max_term, sum_exp
  LOGICAL :: ok_item
  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: z_eta, u_eta, eta_vec

  ok = .TRUE.
  out = 0.d0
  dimEta = n * Q
  CALL activate_growth_subject(i)

  ALLOCATE(z_eta(dimEta), u_eta(dimEta), eta_vec(dimEta))

  n_valid = 0
  max_term = 0.d0
  sum_exp = 0.d0

  DO l = 1, nquad

     baseZ = (l - 1) * dimEta
     DO d = 1, dimEta
        u_eta(d) = nodes(baseZ + d)
     END DO

     z_eta = agh_tt_fe_zeta_eta(:,l,i)

     CALL cfa_eta_mean(i, eta_vec)
     DO r = 1, dimEta
        DO d = 1, r
           eta_vec(r) = eta_vec(r) + L_phi(r,d) * z_eta(d)
        END DO
     END DO

     logg = 0.d0

     DO k_item = 1, nvar
        IF (re_item_flag == 1) THEN
           CALL inneragh_item_logint_full_equiv(i, k_item, eta_vec,              &
                agh_tt_fe_alpha(k_item,l,i), agh_tt_fe_sdev(k_item,i),           &
                logint, ok_item)
        ELSE
           CALL inneragh_item_logint_full_equiv(i, k_item, eta_vec, 0.d0, 1.d0,  &
                logint, ok_item)
        END IF

        IF (.NOT. ok_item) THEN
           logg = LOG_ZERO
           EXIT
        END IF

        logg = logg + logint
     END DO

     IF (logg == LOG_ZERO) CYCLE

     q_eta = DOT_PRODUCT(z_eta, z_eta)
     q_u = DOT_PRODUCT(u_eta, u_eta)

     IF (weights(l) <= 0.d0) CYCLE
     IF (q_eta /= q_eta .OR. q_u /= q_u) CYCLE

     term = LOG(weights(l)) + logg - 0.5d0 * (q_eta - q_u)

     IF (term /= term) CYCLE

     IF (n_valid == 0) THEN
        max_term = term
        sum_exp = 1.d0
        n_valid = 1
     ELSE IF (term > max_term) THEN
        sum_exp = sum_exp * EXP(max_term - term) + 1.d0
        max_term = term
        n_valid = n_valid + 1
     ELSE
        sum_exp = sum_exp + EXP(term - max_term)
        n_valid = n_valid + 1
     END IF

  END DO

  IF (n_valid == 0 .OR. sum_exp <= 0.d0 .OR. sum_exp /= sum_exp) THEN
     out = LOG_ZERO
  ELSE
     out = agh_tt_fe_logconst(i) + max_term + LOG(sum_exp)
  END IF

  IF (out /= out) out = LOG_ZERO

  DEALLOCATE(z_eta, u_eta, eta_vec)

END SUBROUTINE vrais_long_i_tt_agh_full_equiv_fixed


SUBROUTINE vrais_long_tt_agh_full_equiv_fixed(out)
  USE shared, ONLY: nsub, nmes, nvar, agh_tt_fe_grid_ready
  USE datindiv, ONLY: nmescur
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(OUT) :: out

  INTEGER :: i, sm
  DOUBLE PRECISION :: temp
  LOGICAL :: ok_i
  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9

  IF (.NOT. agh_tt_fe_grid_ready) THEN
     CALL rexit('ERROR (vrais_long_tt_agh_full_equiv_fixed): AGH full-equivalent grid not ready')
  END IF

  nmescur = 0
  out = 0.d0

  DO i = 1, nsub

     sm = SUM(nmes(i,:))

     CALL vrais_long_i_tt_agh_full_equiv_fixed(i, sm, temp, ok_i)

     IF (ok_i .AND. temp == temp) THEN
        out = out + temp
     ELSE
        out = LOG_ZERO
     END IF

     nmescur = nmescur + sm

  END DO

  IF (out /= out) out = LOG_ZERO

END SUBROUTINE vrais_long_tt_agh_full_equiv_fixed


SUBROUTINE loglik_long_tt_agh(b0, Y0, nvar0, nsub0, nobs0, n0,      &
                              nmes0, c0, cm10, ntre0,              &
                              nquad_eta0, nitem0,                  &
                              nodes_eta0, weights_eta0,            &
                              nquad_re0, nodes_re0, weights_re0,   &
                              npm0, re_item0, theta_time0,         &
                              Q0, loadmask0, phi_struct0,          &
                              monotone_tau0, update_gh0, fit_id0,  &
                              vrais)                               &
        BIND(C, name="loglik_long_tt_agh")

  USE shared
  USE utilities
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: monotone_tau0
  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0
  INTEGER, INTENT(IN) :: ntre0, nquad_eta0, nquad_re0, npm0
  INTEGER, INTENT(IN) :: re_item0, theta_time0
  INTEGER, INTENT(IN) :: Q0
  INTEGER, INTENT(IN) :: update_gh0
  INTEGER, INTENT(IN) :: fit_id0
  INTEGER, INTENT(IN) :: phi_struct0

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0),       INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0),    INTENT(IN) :: loadmask0

  DOUBLE PRECISION, DIMENSION(npm0),       INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*),          INTENT(IN) :: nodes_eta0
  DOUBLE PRECISION, DIMENSION(nquad_eta0), INTENT(IN) :: weights_eta0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: weights_re0
  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.d9
  LOGICAL :: ok_grid
  DOUBLE PRECISION :: ll_agh

  DOUBLE PRECISION, INTENT(OUT) :: vrais

  INTEGER :: j, k, t, row, col, qq, ll, r
  INTEGER :: pos, expected_npm
  INTEGER :: dimZ, dimEta
  INTEGER :: n_theta_est, n_lambda_free, ref_item
  INTEGER :: ntre_calc

          ! monotonia lcmm-style: incrementi positivi via exp
        DOUBLE PRECISION :: prev
        INTEGER :: h, nf
  INTEGER, DIMENSION(:), ALLOCATABLE :: cm3
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ta

  INTEGER :: n_phi_par


  ! global
  nvar  = nvar0
  nsub  = nsub0
  nobs  = nobs0
  n     = n0
  nquad = nquad_eta0
  nquad_re = nquad_re0
  npm   = npm0

  phi_struct_flag = phi_struct0

  re_item_flag    = re_item0
  theta_time_flag = theta_time0
  Q               = Q0

  p_fix = cfa_p_fix
  p_ran = 0

  dimEta = n * Q
  dim_eta = dimEta

  IF (agh_tt_fit_id_cached /= fit_id0) THEN
   agh_tt_grid_ready = .FALSE.
   agh_tt_fit_id_cached = fit_id0
END IF

  IF (phi_struct0 == 2) THEN
    CALL rexit("ERROR: phi_struct_flag=2 (ricorsione) non supportato con AGH")
  END IF

  ! dati (cleanup)
  IF (ALLOCATED(nmes))    DEALLOCATE(nmes)
  IF (ALLOCATED(Y))       DEALLOCATE(Y)
  IF (ALLOCATED(nitem))   DEALLOCATE(nitem)
  IF (ALLOCATED(c))       DEALLOCATE(c)
  IF (ALLOCATED(cm1))     DEALLOCATE(cm1)

  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs), c(nvar), cm1(nvar))
  nmes  = nmes0
  Y     = Y0
  nitem = nitem0
  c     = c0
  cm1   = cm10

IF (nobs <= 0) THEN
   CALL rexit('ERROR: nobs <= 0')
END IF

IF (SUM(nmes) <= 0) THEN
   CALL rexit('ERROR: SUM(nmes) <= 0. No observations are being used')
END IF

IF (MAXVAL(Y) <= 0) THEN
   CALL rexit('ERROR: Y seems not correctly coded. MAXVAL(Y)')
END IF

IF (MINVAL(nitem) < 1 .OR. MAXVAL(nitem) > n) THEN
   CALL rexit('ERROR: nitem out of range. MIN')
END IF

  IF (ALLOCATED(loadmask)) DEALLOCATE(loadmask)
  ALLOCATE(loadmask(nvar,Q))
  loadmask = loadmask0

  IF (re_item_flag == 1) THEN
     dimZ = dimEta + nvar
  ELSE
     dimZ = dimEta
  END IF


  IF (ALLOCATED(nodes)) DEALLOCATE(nodes)
  IF (ALLOCATED(weights)) DEALLOCATE(weights)
  IF (ALLOCATED(nodes_re)) DEALLOCATE(nodes_re)
  IF (ALLOCATED(weights_re)) DEALLOCATE(weights_re)

  ALLOCATE(nodes(nquad * dimEta))
  nodes(1:nquad * dimEta) = nodes_eta0(1:nquad * dimEta)

  ALLOCATE(weights(nquad))
  weights(1:nquad) = weights_eta0(1:nquad)

  ALLOCATE(nodes_re(nquad_re))
  nodes_re(1:nquad_re) = nodes_re0(1:nquad_re)

  ALLOCATE(weights_re(nquad_re))
  weights_re(1:nquad_re) = weights_re0(1:nquad_re)


  ! parametri (cleanup)
  IF (ALLOCATED(lambda))       DEALLOCATE(lambda)
  IF (ALLOCATED(gamma))        DEALLOCATE(gamma)
  IF (ALLOCATED(theta_k))      DEALLOCATE(theta_k)
  IF (ALLOCATED(invsd_k))      DEALLOCATE(invsd_k)
  IF (ALLOCATED(sigmak))       DEALLOCATE(sigmak)
  IF (ALLOCATED(sqrt_sigmak))  DEALLOCATE(sqrt_sigmak)
  IF (ALLOCATED(kappa))        DEALLOCATE(kappa)
  IF (ALLOCATED(L_phi))        DEALLOCATE(L_phi)
  IF (ALLOCATED(tau))          DEALLOCATE(tau)
  IF (ALLOCATED(theta_long))   DEALLOCATE(theta_long)

  ALLOCATE(lambda(nvar,Q))
  ALLOCATE(gamma(nvar))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(sigmak(nvar),  sqrt_sigmak(nvar))
  ALLOCATE(kappa(dimEta))
  ALLOCATE(L_phi(dimEta,dimEta))
  ALLOCATE(theta_long(nvar,n))

  IF (ALLOCATED(tau_long)) DEALLOCATE(tau_long)
  IF (ALLOCATED(gamma_long)) DEALLOCATE(gamma_long)
  IF (ALLOCATED(lambda_long)) DEALLOCATE(lambda_long)
  ALLOCATE(cm3(nvar), ta(MAX(1, ntre0)), tau(SUM(cm1)))

  ! thresholds: identification controlled by cfa_threshold_id_flag
  DO k = 1, nvar
     cm3(k) = cm1(k) - 2
     IF (cm1(k) < 3) THEN
        CALL rexit('ERROR (CFA): threshold identification requires at least 4 response categories for every item')
     END IF
  END DO

  ntre_calc = SUM(cm3)
  IF (ntre_calc /= ntre0) THEN
     CALL rexit('ERROR (CFA): ntre0 mismatch. ntre0')
  END IF
  ntre = ntre0
  IF (ntre0 > 0) ta(1:ntre0) = b0(1:ntre0)

  pos = 0
  j = 1
  tau = 0.d0
  DO k = 1, nvar
     CALL cfa_fill_tau_item(cm1(k), ta(j), monotone_tau0, tau(pos+1:pos+cm1(k)))
     j = j + cm3(k)
     pos = pos + cm1(k)
  END DO

  ! lambda
  ! lambda
  lambda = 0.d0
  pos = ntre0

  DO qq = 1, Q
     ref_item = 0
     DO k = 1, nvar
        IF (loadmask(k,qq) == 1) THEN
           ref_item = k
           EXIT
        END IF
     END DO
     IF (ALLOCATED(cfa_marker_item)) ref_item = cfa_marker_item(qq)
     IF (ref_item == 0 .OR. loadmask(ref_item,qq) /= 1) THEN
        CALL rexit('ERROR (loglik_long_agh): factor')
     END IF

     lambda(ref_item,qq) = 1.d0

     DO k = 1, nvar
        IF (loadmask(k,qq) == 0) THEN
           lambda(k,qq) = 0.d0
        ELSE
           IF (k == ref_item) CYCLE
           pos = pos + 1
           lambda(k,qq) = b0(pos)
        END IF
     END DO
  END DO

  ! gamma
  DO k = 1, nvar
     pos = pos + 1
     gamma(k) = b0(pos)
  END DO

  ! theta
  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) == 1) THEN
           theta_k(k) = 1.d0
           invsd_k(k) = 1.d0
           DO t = 1, n
              theta_long(k,t) = 1.d0
           END DO
        ELSE
           pos = pos + 1
           theta_k(k) = EXP(b0(pos))
           invsd_k(k) = 1.d0 / SQRT(theta_k(k))
           DO t = 1, n
              theta_long(k,t) = theta_k(k)
           END DO
        END IF
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) == 1) THEN
           theta_k(k) = 1.d0
           invsd_k(k) = 1.d0
           DO t = 1, n
              theta_long(k,t) = 1.d0
           END DO
        ELSE
           DO t = 1, n
              pos = pos + 1
              theta_long(k,t) = EXP(b0(pos))
           END DO
           theta_k(k) = theta_long(k,1)
           invsd_k(k) = 1.d0 / SQRT(theta_k(k))
        END IF
     END DO
  END IF

  ! sigmak
  IF (re_item_flag == 1) THEN
     DO k = 1, nvar
        pos = pos + 1
        sigmak(k)      = EXP(b0(pos))
        sqrt_sigmak(k) = SQRT(sigmak(k))
     END DO
  ELSE
     DO k = 1, nvar
        sigmak(k)      = 0.d0
        sqrt_sigmak(k) = 0.d0
     END DO
  END IF

  ! kappa
  DO qq = 1, Q
     kappa( (1-1)*Q + qq ) = 0.d0
  END DO
  IF (n > 1) THEN
     DO t = 2, n
        DO qq = 1, Q
           pos = pos + 1
           kappa( (t-1)*Q + qq ) = b0(pos)
        END DO
     END DO
  END IF

  ! factor regressions on covariates with time-invariant beta
  CALL ensure_cfa_beta_storage()
  IF (p_fix > 0) THEN
     DO qq = 1, Q
        DO r = 1, p_fix
           pos = pos + 1
           cfa_beta(qq,r) = b0(pos)
        END DO
     END DO
  END IF

  ! L_phi
  L_phi = 0.d0

  IF (phi_struct_flag == 0) THEN
     ! L_phi libero
     DO t = 1, dim_eta
        pos = pos + 1
        L_phi(t,t) = EXP(b0(pos))
     END DO
     DO col = 1, dim_eta
        DO row = col+1, dim_eta
           pos = pos + 1
           L_phi(row,col) = b0(pos)
        END DO
     END DO

  ELSE
     ! AR1: A piena LIBERA, Sigma1 diag, costruisco L_phi
     ALLOCATE(A_ar(Q,Q), sqrt_sig1(Q))

     DO col = 1, Q
       DO row = 1, Q
         pos = pos + 1
         A_ar(row,col) = b0(pos)
       END DO
     END DO

     DO row = 1, Q
       pos = pos + 1
       sqrt_sig1(row) = SQRT(EXP(b0(pos)))
     END DO

     CALL build_Lphi_ar1_full(n, Q, A_ar, sqrt_sig1, L_phi)

     DEALLOCATE(A_ar, sqrt_sig1)
  END IF

  ! expected npm warning
  n_theta_est = 0
  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) n_theta_est = n_theta_est + 1
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) n_theta_est = n_theta_est + n
     END DO
  END IF

  n_lambda_free = 0
  DO qq = 1, Q
     DO k = 1, nvar
        IF (loadmask(k,qq) == 1) n_lambda_free = n_lambda_free + 1
     END DO
  END DO
  n_lambda_free = n_lambda_free - Q

  IF (phi_struct_flag == 0) THEN
     n_phi_par = dim_eta*(dim_eta+1)/2
  ELSE
     n_phi_par = Q*Q + Q
  END IF

  expected_npm = ntre0                            &
              + n_lambda_free                     &
              + nvar                              &
              + n_theta_est                       &
              + MERGE(nvar, 0, re_item_flag == 1) &
              + (n-1)*Q                           &
              + Q*p_fix                          &
              + n_phi_par

  IF (npm0 /= expected_npm) THEN
     CALL rexit("The supplied parameter vector has the wrong length for this CFA model.")
  END IF

  ! ! -loglik AGH
  ! vrais = -vrais_long_long_agh()
  ! -loglik AGH fast, stile GLMMadaptive

 IF (update_gh0 == 1 .OR. .NOT. agh_tt_grid_ready) THEN

   CALL build_agh_grid_eta_nested(ok_grid)

   IF (.NOT. ok_grid .OR. .NOT. agh_tt_grid_ready) THEN
      vrais = BAD_OBJ
      IF (ALLOCATED(cm3)) DEALLOCATE(cm3)
      IF (ALLOCATED(ta)) DEALLOCATE(ta)
      RETURN
   END IF

END IF

CALL vrais_long_tt_agh_fixed(ll_agh)


IF (ll_agh /= ll_agh .OR. ABS(ll_agh) > 1.d300) THEN
   vrais = BAD_OBJ
ELSE
   vrais = -ll_agh
END IF


  ! print *, vrais

  DEALLOCATE(cm3, ta)

END SUBROUTINE loglik_long_tt_agh


SUBROUTINE loglik_long_tt_agh_full_equiv(b0, Y0, nvar0, nsub0, nobs0, n0,      &
                              nmes0, c0, cm10, ntre0,              &
                              nquad_eta0, nitem0,                  &
                              nodes_eta0, weights_eta0,            &
                              nquad_re0, nodes_re0, weights_re0,   &
                              npm0, re_item0, theta_time0,         &
                              meas_invar0,                         &
                              Q0, loadmask0, phi_struct0,          &
                              monotone_tau0, update_gh0, fit_id0,  &
                              vrais)                               &
        BIND(C, name="loglik_long_tt_agh_full_equiv")

  USE shared
  USE utilities
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: monotone_tau0
  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0
  INTEGER, INTENT(IN) :: ntre0, nquad_eta0, nquad_re0, npm0
  INTEGER, INTENT(IN) :: re_item0, theta_time0
  INTEGER, DIMENSION(4), INTENT(IN) :: meas_invar0
  INTEGER, INTENT(IN) :: Q0
  INTEGER, INTENT(IN) :: update_gh0
  INTEGER, INTENT(IN) :: fit_id0
  INTEGER, INTENT(IN) :: phi_struct0

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0),       INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0),    INTENT(IN) :: loadmask0

  DOUBLE PRECISION, DIMENSION(npm0),       INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*),          INTENT(IN) :: nodes_eta0
  DOUBLE PRECISION, DIMENSION(nquad_eta0), INTENT(IN) :: weights_eta0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: weights_re0
  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.d9
  LOGICAL :: ok_grid
  DOUBLE PRECISION :: ll_agh

  DOUBLE PRECISION, INTENT(OUT) :: vrais

  INTEGER :: j, k, t, row, col, qq, ll, r
  INTEGER :: pos, expected_npm
  INTEGER :: dimZ, dimEta
  INTEGER :: n_theta_est, n_lambda_free, ref_item
  INTEGER :: ntre_calc
  INTEGER :: ntre_base, n_lambda_base, n_gamma_base, n_gamma_est

          ! monotonia lcmm-style: incrementi positivi via exp
        DOUBLE PRECISION :: prev
        INTEGER :: h, nf
  INTEGER, DIMENSION(:), ALLOCATABLE :: cm3
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ta

  INTEGER :: n_phi_par


  ! global
  nvar  = nvar0
  nsub  = nsub0
  nobs  = nobs0
  n     = n0
  nquad = nquad_eta0
  nquad_re = nquad_re0
  npm   = npm0

  phi_struct_flag = phi_struct0

  re_item_flag    = re_item0
  theta_time_flag = theta_time0
  Q               = Q0

  p_fix = cfa_p_fix
  p_ran = 0

  IF (ALLOCATED(kappa_sub)) DEALLOCATE(kappa_sub)
  IF (ALLOCATED(L_phi_sub)) DEALLOCATE(L_phi_sub)

  dimEta = n * Q
  dim_eta = dimEta

  IF (SIZE(meas_invar0) /= 4) THEN
     CALL rexit('ERROR (loglik_long_tt_agh_full_equiv): meas_invar0 must have length 4')
  END IF

  IF (MINVAL(meas_invar0) < 0 .OR. MAXVAL(meas_invar0) > 1) THEN
     CALL rexit('ERROR (loglik_long_tt_agh_full_equiv): meas_invar0 must contain only 0/1')
  END IF

  IF (tau_invar_flag    /= meas_invar0(1) .OR. &
      gamma_invar_flag  /= meas_invar0(2) .OR. &
      lambda_invar_flag /= meas_invar0(3) .OR. &
      theta_invar_flag  /= meas_invar0(4)) THEN
     agh_tt_fe_grid_ready = .FALSE.
  END IF

  tau_invar_flag    = meas_invar0(1)
  gamma_invar_flag  = meas_invar0(2)
  lambda_invar_flag = meas_invar0(3)
  theta_invar_flag  = meas_invar0(4)

  theta_time_flag = 1 - theta_invar_flag

  IF (theta_time0 /= theta_time_flag) THEN
     CALL rwarn("The residual invariance flag overrides theta_time0.")
  END IF

  IF (agh_tt_fe_fit_id_cached /= fit_id0) THEN
   agh_tt_fe_grid_ready = .FALSE.
   agh_tt_fe_fit_id_cached = fit_id0
END IF

  IF (phi_struct0 == 2) THEN
    CALL rexit("ERROR: phi_struct_flag=2 (ricorsione) non supportato con AGH")
  END IF

  ! dati (cleanup)
  IF (ALLOCATED(nmes))    DEALLOCATE(nmes)
  IF (ALLOCATED(Y))       DEALLOCATE(Y)
  IF (ALLOCATED(nitem))   DEALLOCATE(nitem)
  IF (ALLOCATED(c))       DEALLOCATE(c)
  IF (ALLOCATED(cm1))     DEALLOCATE(cm1)

  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs), c(nvar), cm1(nvar))
  nmes  = nmes0
  Y     = Y0
  nitem = nitem0
  c     = c0
  cm1   = cm10

IF (nobs <= 0) THEN
   CALL rexit('ERROR: nobs <= 0')
END IF

IF (SUM(nmes) <= 0) THEN
   CALL rexit('ERROR: SUM(nmes) <= 0. No observations are being used')
END IF

IF (MAXVAL(Y) <= 0) THEN
   CALL rexit('ERROR: Y seems not correctly coded. MAXVAL(Y)')
END IF

IF (MINVAL(nitem) < 1 .OR. MAXVAL(nitem) > n) THEN
   CALL rexit('ERROR: nitem out of range. MIN')
END IF

  IF (ALLOCATED(loadmask)) DEALLOCATE(loadmask)
  ALLOCATE(loadmask(nvar,Q))
  loadmask = loadmask0

  IF (re_item_flag == 1) THEN
     dimZ = dimEta + nvar
  ELSE
     dimZ = dimEta
  END IF


  IF (ALLOCATED(nodes)) DEALLOCATE(nodes)
  IF (ALLOCATED(weights)) DEALLOCATE(weights)
  IF (ALLOCATED(nodes_re)) DEALLOCATE(nodes_re)
  IF (ALLOCATED(weights_re)) DEALLOCATE(weights_re)

  ALLOCATE(nodes(nquad * dimEta))
  nodes(1:nquad * dimEta) = nodes_eta0(1:nquad * dimEta)

  ALLOCATE(weights(nquad))
  weights(1:nquad) = weights_eta0(1:nquad)

  ALLOCATE(nodes_re(nquad_re))
  nodes_re(1:nquad_re) = nodes_re0(1:nquad_re)

  ALLOCATE(weights_re(nquad_re))
  weights_re(1:nquad_re) = weights_re0(1:nquad_re)


  ! parametri (cleanup)
  IF (ALLOCATED(lambda))       DEALLOCATE(lambda)
  IF (ALLOCATED(lambda_long))  DEALLOCATE(lambda_long)
  IF (ALLOCATED(gamma))        DEALLOCATE(gamma)
  IF (ALLOCATED(gamma_long))   DEALLOCATE(gamma_long)
  IF (ALLOCATED(theta_k))      DEALLOCATE(theta_k)
  IF (ALLOCATED(invsd_k))      DEALLOCATE(invsd_k)
  IF (ALLOCATED(sigmak))       DEALLOCATE(sigmak)
  IF (ALLOCATED(sqrt_sigmak))  DEALLOCATE(sqrt_sigmak)
  IF (ALLOCATED(kappa))        DEALLOCATE(kappa)
  IF (ALLOCATED(L_phi))        DEALLOCATE(L_phi)
  IF (ALLOCATED(tau))          DEALLOCATE(tau)
  IF (ALLOCATED(tau_long))     DEALLOCATE(tau_long)
  IF (ALLOCATED(theta_long))   DEALLOCATE(theta_long)
  IF (ALLOCATED(invsd_long))   DEALLOCATE(invsd_long)

  ALLOCATE(lambda(nvar,Q))
  ALLOCATE(lambda_long(nvar,Q,n))
  ALLOCATE(gamma(nvar))
  ALLOCATE(gamma_long(nvar,n))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(sigmak(nvar),  sqrt_sigmak(nvar))
  ALLOCATE(kappa(dimEta))
  ALLOCATE(L_phi(dimEta,dimEta))
  ALLOCATE(theta_long(nvar,n))
  ALLOCATE(invsd_long(nvar,n))

  ! thresholds: identification controlled by cfa_threshold_id_flag
  ALLOCATE(cm3(nvar))
  DO k = 1, nvar
     cm3(k) = cm1(k) - 2
     IF (cm1(k) < 3) THEN
        CALL rexit('ERROR (CFA): threshold identification requires at least 4 response categories for every item')
     END IF
  END DO

  ntre_base = SUM(cm3)
  IF (tau_invar_flag == 1) THEN
     ntre_calc = ntre_base
  ELSE
     ntre_calc = n * ntre_base
  END IF
  IF (ntre_calc /= ntre0) THEN
     CALL rexit('ERROR (loglik_long_tt_agh_full_equiv): ntre0 mismatch. ntre0')
  END IF
  ntre = ntre0

  ALLOCATE(ta(MAX(1,ntre0)))
  IF (ntre0 > 0) ta(1:ntre0) = b0(1:ntre0)
  ALLOCATE(tau(SUM(cm1)))
  ALLOCATE(tau_long(SUM(cm1),n))
  tau = 0.d0
  tau_long = 0.d0
  j = 1

  IF (tau_invar_flag == 1) THEN
     pos = 0
     DO k = 1, nvar
        CALL cfa_fill_tau_item(cm1(k), ta(j), monotone_tau0, tau_long(pos+1:pos+cm1(k),1))
        DO t = 2, n
           tau_long(pos+1:pos+cm1(k),t) = tau_long(pos+1:pos+cm1(k),1)
        END DO
        j = j + cm3(k)
        pos = pos + cm1(k)
     END DO
  ELSE
     DO t = 1, n
        pos = 0
        DO k = 1, nvar
           CALL cfa_fill_tau_item(cm1(k), ta(j), monotone_tau0, tau_long(pos+1:pos+cm1(k),t))
           j = j + cm3(k)
           pos = pos + cm1(k)
        END DO
     END DO
  END IF
  tau = tau_long(:,1)

  ! lambda
  ! lambda
  lambda = 0.d0
  lambda_long = 0.d0
  pos = ntre0

  IF (lambda_invar_flag == 1) THEN
     DO qq = 1, Q
        ref_item = 0
        DO k = 1, nvar
           IF (loadmask(k,qq) == 1) THEN
              ref_item = k
              EXIT
           END IF
        END DO
        IF (ALLOCATED(cfa_marker_item)) ref_item = cfa_marker_item(qq)
      IF (ref_item == 0 .OR. loadmask(ref_item,qq) /= 1) THEN
           CALL rexit('ERROR (loglik_long_tt_agh_full_equiv): factor')
        END IF

        lambda(ref_item,qq) = 1.d0

        DO k = 1, nvar
           IF (loadmask(k,qq) == 0) THEN
              lambda(k,qq) = 0.d0
           ELSE
              IF (k == ref_item) CYCLE
              pos = pos + 1
              lambda(k,qq) = b0(pos)
           END IF
        END DO
     END DO

     DO t = 1, n
        lambda_long(:,:,t) = lambda
     END DO
  ELSE
     DO t = 1, n
        DO qq = 1, Q
           ref_item = 0
           DO k = 1, nvar
              IF (loadmask(k,qq) == 1) THEN
                 ref_item = k
                 EXIT
              END IF
           END DO
           IF (ALLOCATED(cfa_marker_item)) ref_item = cfa_marker_item(qq)
      IF (ref_item == 0 .OR. loadmask(ref_item,qq) /= 1) THEN
              CALL rexit('ERROR (loglik_long_tt_agh_full_equiv): factor')
           END IF

           lambda_long(ref_item,qq,t) = 1.d0

           DO k = 1, nvar
              IF (loadmask(k,qq) == 0) THEN
                 lambda_long(k,qq,t) = 0.d0
              ELSE
                 IF (k == ref_item) CYCLE
                 pos = pos + 1
                 lambda_long(k,qq,t) = b0(pos)
              END IF
           END DO
        END DO
     END DO

     lambda = lambda_long(:,:,1)
  END IF

  ! gamma
  gamma = 0.d0
  gamma_long = 0.d0

  IF (cfa_mean_ident_flag == 1 .AND. .NOT. ALLOCATED(cfa_marker_item)) THEN
     CALL rexit('ERROR (loglik_long_tt_agh_full_equiv): cfa_marker_item not allocated')
  END IF

  IF (gamma_invar_flag == 1) THEN
     DO k = 1, nvar
        IF (cfa_mean_ident_flag == 1) THEN
           IF (ANY(cfa_marker_item == k)) CYCLE
        END IF
        pos = pos + 1
        gamma(k) = b0(pos)
        gamma_long(k,:) = gamma(k)
     END DO
     IF (ALLOCATED(cfa_gamma_time_free)) THEN
        IF (SIZE(cfa_gamma_time_free,1) == nvar .AND. SIZE(cfa_gamma_time_free,2) == n) THEN
           DO t = 2, n
              DO k = 1, nvar
                 IF (cfa_gamma_time_free(k,t) /= 1) CYCLE
                 IF (cfa_mean_ident_flag == 1) THEN
                    IF (ANY(cfa_marker_item == k)) CYCLE
                 END IF
                 pos = pos + 1
                 gamma_long(k,t) = b0(pos)
              END DO
           END DO
        END IF
     END IF
  ELSE
     DO t = 1, n
        DO k = 1, nvar
           IF (cfa_mean_ident_flag == 1) THEN
              IF (ANY(cfa_marker_item == k)) CYCLE
           END IF
           pos = pos + 1
           gamma_long(k,t) = b0(pos)
        END DO
     END DO
     gamma = gamma_long(:,1)
  END IF

  ! theta
  theta_long = 1.d0
  invsd_long = 1.d0

  IF (theta_invar_flag == 1) THEN
     DO k = 1, nvar
        IF (cm1(k) == 1) THEN
           theta_k(k) = 1.d0
           invsd_k(k) = 1.d0
           theta_long(k,:) = 1.d0
           invsd_long(k,:) = 1.d0
        ELSE
           pos = pos + 1
           theta_k(k) = EXP(b0(pos))
           invsd_k(k) = 1.d0 / SQRT(theta_k(k))
           theta_long(k,:) = theta_k(k)
           invsd_long(k,:) = invsd_k(k)
        END IF
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) == 1) THEN
           theta_k(k) = 1.d0
           invsd_k(k) = 1.d0
           theta_long(k,:) = 1.d0
           invsd_long(k,:) = 1.d0
        ELSE
           DO t = 1, n
              pos = pos + 1
              theta_long(k,t) = EXP(b0(pos))
              invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
           END DO
           theta_k(k) = theta_long(k,1)
           invsd_k(k) = invsd_long(k,1)
        END IF
     END DO
  END IF

  ! sigmak
  IF (re_item_flag == 1) THEN
     DO k = 1, nvar
        pos = pos + 1
        sigmak(k)      = EXP(b0(pos))
        sqrt_sigmak(k) = SQRT(sigmak(k))
     END DO
  ELSE
     DO k = 1, nvar
        sigmak(k)      = 0.d0
        sqrt_sigmak(k) = 0.d0
     END DO
  END IF

  ! kappa
  kappa = 0.d0
  IF (cfa_mean_ident_flag == 1) THEN
     DO t = 1, n
        DO qq = 1, Q
           pos = pos + 1
           kappa( (t-1)*Q + qq ) = b0(pos)
        END DO
     END DO
  ELSE
     IF (n > 1) THEN
        DO t = 2, n
           DO qq = 1, Q
              pos = pos + 1
              kappa( (t-1)*Q + qq ) = b0(pos)
           END DO
        END DO
     END IF
  END IF

  ! factor regressions on covariates with time-invariant beta
  CALL ensure_cfa_beta_storage()
  IF (p_fix > 0) THEN
     DO qq = 1, Q
        DO r = 1, p_fix
           pos = pos + 1
           cfa_beta(qq,r) = b0(pos)
        END DO
     END DO
  END IF

  ! L_phi
  L_phi = 0.d0

  IF (phi_struct_flag == 0) THEN
     ! L_phi libero
     DO t = 1, dim_eta
        pos = pos + 1
        L_phi(t,t) = EXP(b0(pos))
     END DO
     DO col = 1, dim_eta
        DO row = col+1, dim_eta
           pos = pos + 1
           L_phi(row,col) = b0(pos)
        END DO
     END DO

  ELSE
     ! AR1: A piena LIBERA, Sigma1 diag, costruisco L_phi
     ALLOCATE(A_ar(Q,Q), sqrt_sig1(Q))

     DO col = 1, Q
       DO row = 1, Q
         pos = pos + 1
         A_ar(row,col) = b0(pos)
       END DO
     END DO

     DO row = 1, Q
       pos = pos + 1
       sqrt_sig1(row) = SQRT(EXP(b0(pos)))
     END DO

     CALL build_Lphi_ar1_full(n, Q, A_ar, sqrt_sig1, L_phi)

     DEALLOCATE(A_ar, sqrt_sig1)
  END IF

  ! expected npm warning
  n_theta_est = 0
  IF (theta_invar_flag == 1) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) n_theta_est = n_theta_est + 1
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) n_theta_est = n_theta_est + n
     END DO
  END IF

  n_lambda_base = 0
  DO qq = 1, Q
     DO k = 1, nvar
        IF (loadmask(k,qq) == 1) n_lambda_base = n_lambda_base + 1
     END DO
  END DO
  n_lambda_base = n_lambda_base - Q

  IF (lambda_invar_flag == 1) THEN
     n_lambda_free = n_lambda_base
  ELSE
     n_lambda_free = n * n_lambda_base
  END IF

  n_gamma_base = nvar
  IF (cfa_mean_ident_flag == 1) THEN
     n_gamma_base = 0
     DO k = 1, nvar
        IF (.NOT. ANY(cfa_marker_item == k)) n_gamma_base = n_gamma_base + 1
     END DO
  END IF

  IF (gamma_invar_flag == 1) THEN
     n_gamma_est = n_gamma_base
     IF (ALLOCATED(cfa_gamma_time_free) .AND. n > 1) THEN
        IF (SIZE(cfa_gamma_time_free,1) == nvar .AND. SIZE(cfa_gamma_time_free,2) == n) THEN
           DO t = 2, n
              DO k = 1, nvar
                 IF (cfa_gamma_time_free(k,t) /= 1) CYCLE
                 IF (cfa_mean_ident_flag == 1) THEN
                    IF (ANY(cfa_marker_item == k)) CYCLE
                 END IF
                 n_gamma_est = n_gamma_est + 1
              END DO
           END DO
        END IF
     END IF
  ELSE
     n_gamma_est = n * n_gamma_base
  END IF

  IF (phi_struct_flag == 0) THEN
     n_phi_par = dim_eta*(dim_eta+1)/2
  ELSE
     n_phi_par = Q*Q + Q
  END IF

  expected_npm = ntre_calc                         &
              + n_lambda_free                     &
              + n_gamma_est                       &
              + n_theta_est                       &
              + MERGE(nvar, 0, re_item_flag == 1) &
              + MERGE(n*Q, (n-1)*Q, cfa_mean_ident_flag == 1) &
              + Q*p_fix                          &
              + n_phi_par

  IF (npm0 /= expected_npm) THEN
     CALL rexit("The supplied parameter vector has the wrong length for this CFA model.")
  END IF

  ! ! -loglik AGH
  ! vrais = -vrais_long_long_agh()
  ! -loglik AGH fast, stile GLMMadaptive

 IF (update_gh0 == 1 .OR. .NOT. agh_tt_fe_grid_ready) THEN

   CALL build_agh_grid_full_equiv_tt(ok_grid)

   IF (.NOT. ok_grid .OR. .NOT. agh_tt_fe_grid_ready) THEN
      vrais = BAD_OBJ
      IF (ALLOCATED(cm3)) DEALLOCATE(cm3)
      IF (ALLOCATED(ta)) DEALLOCATE(ta)
      RETURN
   END IF

END IF

CALL vrais_long_tt_agh_full_equiv_fixed(ll_agh)


IF (ll_agh /= ll_agh .OR. ABS(ll_agh) > 1.d300) THEN
   vrais = BAD_OBJ
ELSE
   vrais = -ll_agh
END IF


  ! print *, vrais

  DEALLOCATE(cm3, ta)

END SUBROUTINE loglik_long_tt_agh_full_equiv





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                            !!!
!!!          Two-tier          !!!
!!!                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!===========================================================
! two tier longitudinal cfa
! primary dimensions: eta_{t,q}, t=1..n, q=1..Q
! item specific random effects: a_ik, mutually independent
! and independent of eta
!===========================================================
SUBROUTINE build_tt_precomp()
  USE shared, ONLY: nvar, nsub, n, Q, nquad_eta, nodes_eta, kappa, L_phi, &
                    phi_struct_flag, A_ar, sqrt_sig1, gamma, lambda, &
                    nload_item, load_idx, lambda_act, tt_precomp_ready, &
                    tt_precomp_dimEta, tt_precomp_nquad_eta, tt_precomp_nvar, &
                    tt_precomp_n, tt_precomp_Q, tt_eta_grid, tt_mu0_grid
  USE utilities, ONLY: eta_from_z_ar1
  IMPLICIT NONE

  INTEGER :: dimEta, l_eta, baseZ, idim, irow, icol, i
  INTEGER :: k_item, t, ll, idx_tq, base_eta
  LOGICAL :: need_alloc
  DOUBLE PRECISION :: linpred
  DOUBLE PRECISION, ALLOCATABLE :: z_lat(:), mean_eta(:)

  dimEta = n * Q

  need_alloc = (.NOT. tt_precomp_ready) .OR. &
               (.NOT. ALLOCATED(tt_eta_grid)) .OR. &
               (.NOT. ALLOCATED(tt_mu0_grid)) .OR. &
               (tt_precomp_dimEta /= dimEta) .OR. &
               (tt_precomp_nquad_eta /= nquad_eta) .OR. &
               (tt_precomp_nvar /= nvar) .OR. &
               (tt_precomp_n /= n) .OR. &
               (tt_precomp_Q /= Q)

  IF (.NOT. need_alloc) THEN
     need_alloc = SIZE(tt_mu0_grid,4) /= nsub
  END IF

  IF (need_alloc) THEN
     IF (ALLOCATED(tt_eta_grid)) DEALLOCATE(tt_eta_grid)
     IF (ALLOCATED(tt_mu0_grid)) DEALLOCATE(tt_mu0_grid)

     ALLOCATE(tt_eta_grid(dimEta, nquad_eta))
     ALLOCATE(tt_mu0_grid(nvar, n, nquad_eta, nsub))

     tt_precomp_dimEta     = dimEta
     tt_precomp_nquad_eta  = nquad_eta
     tt_precomp_nvar       = nvar
     tt_precomp_n          = n
     tt_precomp_Q          = Q
     tt_precomp_ready      = .TRUE.
  END IF

  ALLOCATE(z_lat(dimEta), mean_eta(dimEta))

  DO l_eta = 1, nquad_eta
     baseZ = (l_eta - 1) * dimEta
     DO idim = 1, dimEta
        z_lat(idim) = nodes_eta(baseZ + idim)
     END DO

     DO i = 1, nsub
        CALL cfa_eta_mean(i, mean_eta)

        IF (phi_struct_flag == 2) THEN
           IF (.NOT. ALLOCATED(A_ar) .OR. .NOT. ALLOCATED(sqrt_sig1)) THEN
              CALL rexit('ERROR (build_tt_precomp): A_ar/sqrt_sig1 not allocated')
           END IF
           CALL eta_from_z_ar1(n, Q, mean_eta, A_ar, sqrt_sig1, z_lat, tt_eta_grid(:, l_eta))
        ELSE
           DO irow = 1, dimEta
              tt_eta_grid(irow, l_eta) = mean_eta(irow)
              DO icol = 1, irow
                 tt_eta_grid(irow, l_eta) = tt_eta_grid(irow, l_eta) + &
                                            L_phi(irow, icol) * z_lat(icol)
              END DO
           END DO
        END IF

        IF (Q == 1) THEN
           DO t = 1, n
              base_eta = t
              DO k_item = 1, nvar
                 tt_mu0_grid(k_item, t, l_eta, i) = gamma(k_item) + &
                                                    lambda(k_item, 1) * &
                                                    tt_eta_grid(base_eta, l_eta)
              END DO
           END DO
        ELSE
           DO t = 1, n
              base_eta = (t - 1) * Q
              DO k_item = 1, nvar
                 linpred = 0.d0
                 DO ll = 1, nload_item(k_item)
                    idx_tq  = base_eta + load_idx(ll, k_item)
                    linpred = linpred + lambda_act(ll, k_item) * &
                                         tt_eta_grid(idx_tq, l_eta)
                 END DO
                 tt_mu0_grid(k_item, t, l_eta, i) = gamma(k_item) + linpred
              END DO
           END DO
        END IF
     END DO
  END DO

  DEALLOCATE(z_lat, mean_eta)
END SUBROUTINE build_tt_precomp


DOUBLE PRECISION FUNCTION log_item_block_tt(i, obs_start, k_item, l_eta)
  USE shared, ONLY: n, re_item_flag, nquad_re, nodes_re, logweights_re, &
                    nmes, c, Y, nitem, invsd_long, sqrt_sigmak, tau, &
                    kk_base, tt_mu0_grid
  USE utilities, ONLY: alnorm_fast
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, obs_start, k_item, l_eta

  INTEGER :: j_mes, idx_obs, yijk, t_idx
  INTEGER :: kk0, nmes_k, c_k, l_re
  INTEGER :: y_vec(n), tuk_vec(n)

  DOUBLE PRECISION :: xmu_vec(n), invsd_vec(n)
  DOUBLE PRECISION :: mu_ijk, a_k
  DOUBLE PRECISION :: loglik_re, logterm, log_max, sumexp

  nmes_k = nmes(i, k_item)

  IF (nmes_k == 0) THEN
     log_item_block_tt = 0.d0
     RETURN
  END IF

  kk0 = kk_base(k_item)
  c_k = c(k_item)

  idx_obs = obs_start - 1

  DO j_mes = 1, nmes_k
     idx_obs = idx_obs + 1

     t_idx = nitem(idx_obs)
     yijk  = Y(idx_obs)

     y_vec(j_mes)     = yijk
     tuk_vec(j_mes)   = yijk + kk0
     invsd_vec(j_mes) = invsd_long(k_item, t_idx)
     xmu_vec(j_mes)   = tt_mu0_grid(k_item, t_idx, l_eta, i)
  END DO

  IF (re_item_flag == 0) THEN
     loglik_re = 0.d0

     DO j_mes = 1, nmes_k
        loglik_re = loglik_re + log_prob_obs_tt(xmu_vec(j_mes),       &
                                                invsd_vec(j_mes),     &
                                                y_vec(j_mes),         &
                                                tuk_vec(j_mes), c_k)
     END DO

     log_item_block_tt = loglik_re
     RETURN
  END IF

  log_max = -HUGE(1.d0)
  sumexp  = 0.d0

  DO l_re = 1, nquad_re
     a_k       = sqrt_sigmak(k_item) * nodes_re(l_re)
     loglik_re = 0.d0

     DO j_mes = 1, nmes_k
        mu_ijk    = xmu_vec(j_mes) + a_k
        loglik_re = loglik_re + log_prob_obs_tt(mu_ijk,               &
                                                invsd_vec(j_mes),     &
                                                y_vec(j_mes),         &
                                                tuk_vec(j_mes), c_k)
     END DO

     logterm = logweights_re(l_re) + loglik_re

     IF (logterm > log_max) THEN
        sumexp  = sumexp * EXP(log_max - logterm) + 1.d0
        log_max = logterm
     ELSE
        sumexp = sumexp + EXP(logterm - log_max)
     END IF
  END DO

  log_item_block_tt = log_max + LOG(sumexp)

CONTAINS

  DOUBLE PRECISION FUNCTION log_prob_obs_tt(mu, invsd, y, tuk, ck)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: mu, invsd
    INTEGER, INTENT(IN) :: y, tuk, ck

    DOUBLE PRECISION :: tausup, tauinf, p

    IF (ck == 2) THEN
       tausup = -mu * invsd

       IF (y == 1) THEN
          p = alnorm_fast(tausup, .FALSE.)
       ELSE
          p = alnorm_fast(tausup, .TRUE.)
       END IF

    ELSE

       IF (y == 1) THEN
          tausup = (tau(tuk) - mu) * invsd
          p      = alnorm_fast(tausup, .FALSE.)

       ELSE IF (y == ck) THEN
          tauinf = (tau(tuk-1) - mu) * invsd
          p      = alnorm_fast(tauinf, .TRUE.)

       ELSE
          tausup = (tau(tuk)   - mu) * invsd
          tauinf = (tau(tuk-1) - mu) * invsd

          IF (tauinf > 0.d0) THEN
             p = alnorm_fast(tauinf, .TRUE.) - alnorm_fast(tausup, .TRUE.)
          ELSE
             p = alnorm_fast(tausup, .FALSE.) - alnorm_fast(tauinf, .FALSE.)
          END IF
       END IF

    END IF

    log_prob_obs_tt = LOG(p)

  END FUNCTION log_prob_obs_tt

END FUNCTION log_item_block_tt

DOUBLE PRECISION FUNCTION vrais_long_tt_i(i, obs0)
  USE shared, ONLY: nvar, nmes, nquad_eta, logweights_eta
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, obs0

  INTEGER :: l_eta, k_item, obs_start
  DOUBLE PRECISION :: loglik_eta, logterm
  DOUBLE PRECISION :: log_max, sumexp

  log_max = -HUGE(1.d0)
  sumexp  = 0.d0

  DO l_eta = 1, nquad_eta
     loglik_eta = 0.d0
     obs_start  = obs0

     DO k_item = 1, nvar
        IF (nmes(i, k_item) > 0) THEN
           loglik_eta = loglik_eta + log_item_block_tt(i, obs_start, &
                                                       k_item, l_eta)
        END IF
        obs_start = obs_start + nmes(i, k_item)
     END DO

     logterm = logweights_eta(l_eta) + loglik_eta

     IF (logterm > log_max) THEN
        sumexp  = sumexp * EXP(log_max - logterm) + 1.d0
        log_max = logterm
     ELSE
        sumexp = sumexp + EXP(logterm - log_max)
     END IF
  END DO

  vrais_long_tt_i = log_max + LOG(sumexp)

END FUNCTION vrais_long_tt_i


DOUBLE PRECISION FUNCTION vrais_long_tt()
  USE shared, ONLY: nsub, nmes
  IMPLICIT NONE

  INTEGER :: i, sm, obs0
  DOUBLE PRECISION :: temp

  obs0          = 1
  vrais_long_tt = 0.d0

  DO i = 1, nsub
     sm   = SUM(nmes(i,:))
     temp = vrais_long_tt_i(i, obs0)

     vrais_long_tt = vrais_long_tt + temp

     obs0 = obs0 + sm
  END DO

END FUNCTION vrais_long_tt


SUBROUTINE loglik_long_tt(b0, Y0, nvar0, nsub0, nobs0, n0,       &
                       nmes0, c0, cm10, ntre0,                  &
                       nquad_eta0, nodes_eta0, weights_eta0,    &
                       nquad_re0, nodes_re0, weights_re0,       &
                       nitem0, npm0, re_item0, theta_time0,     &
                       Q0, loadmask0, phi_struct0,              &
                       monotone_tau0, fit_id0,                  &
                       vrais)                                   &
        BIND(C, name="loglik_long_tt")

  USE shared
  USE utilities
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: fit_id0
  INTEGER, INTENT(IN) :: monotone_tau0
  INTEGER, INTENT(IN) :: phi_struct0
  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0
  INTEGER, INTENT(IN) :: ntre0, npm0
  INTEGER, INTENT(IN) :: nquad_eta0, nquad_re0
  INTEGER, INTENT(IN) :: re_item0, theta_time0
  INTEGER, INTENT(IN) :: Q0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0),       INTENT(IN) :: Y0, nitem0

  DOUBLE PRECISION, DIMENSION(*),           INTENT(IN) :: nodes_eta0
  DOUBLE PRECISION, DIMENSION(nquad_eta0),  INTENT(IN) :: weights_eta0
  DOUBLE PRECISION, DIMENSION(*),           INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0),   INTENT(IN) :: weights_re0
  DOUBLE PRECISION, DIMENSION(npm0),        INTENT(IN) :: b0

  DOUBLE PRECISION, INTENT(OUT) :: vrais

  INTEGER :: j, k, t, row, col, qq, ll, r, maxload
  INTEGER :: pos, expected_npm
  INTEGER :: dimEta
  INTEGER :: n_theta_est, n_lambda_free, ref_item
  INTEGER :: ntre_calc
  INTEGER :: cm3(nvar0)
  INTEGER :: n_phi_par
  INTEGER :: h

  LOGICAL :: need_data_init, need_param_alloc

  DOUBLE PRECISION :: prev, ll_tt
  DOUBLE PRECISION :: ta(MAX(1, ntre0))
  DOUBLE PRECISION :: A_tmp(Q0,Q0), sqrt1_tmp(Q0)

  nvar  = nvar0
  nsub  = nsub0
  nobs  = nobs0
  n     = n0
  npm   = npm0

  nquad_eta = nquad_eta0
  nquad_re  = nquad_re0

  phi_struct_flag = phi_struct0
  re_item_flag    = re_item0
  theta_time_flag = theta_time0
  Q               = Q0

  p_fix = cfa_p_fix
  p_ran = 0

  dimEta  = n * Q
  dim_eta = dimEta

  need_data_init = (.NOT. cache_ready) .OR. &
                   (fit_id0 /= cache_fit_id) .OR. &
                   (nvar0 /= cache_nvar) .OR. &
                   (nsub0 /= cache_nsub) .OR. &
                   (nobs0 /= cache_nobs) .OR. &
                   (n0   /= cache_n) .OR. &
                   (Q0   /= cache_Q) .OR. &
                   (nquad_eta0 /= cache_nquad_eta) .OR. &
                   (MERGE(nquad_re0, 0, re_item0 == 1) /= cache_nquad_re)

  IF (need_data_init) THEN

    IF (ALLOCATED(nmes))    DEALLOCATE(nmes)
    IF (ALLOCATED(Y))       DEALLOCATE(Y)
    IF (ALLOCATED(nitem))   DEALLOCATE(nitem)
    IF (ALLOCATED(c))       DEALLOCATE(c)
    IF (ALLOCATED(cm1))     DEALLOCATE(cm1)

    ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs), c(nvar), cm1(nvar))

    nmes  = nmes0
    Y     = Y0
    nitem = nitem0
    c     = c0
    cm1   = cm10

    IF (ALLOCATED(kk_base)) DEALLOCATE(kk_base)
    ALLOCATE(kk_base(nvar))

    kk_base(1) = 0
    DO k = 2, nvar
      kk_base(k) = kk_base(k-1) + (c(k-1) - 1)
    END DO

    IF (ALLOCATED(subj_start)) DEALLOCATE(subj_start)
    ALLOCATE(subj_start(nsub+1))

    subj_start(1) = 1
    DO k = 1, nsub
      subj_start(k+1) = subj_start(k) + SUM(nmes(k,:))
    END DO

    IF (subj_start(nsub+1) /= nobs + 1) THEN
       CALL rexit('ERROR (loglik_long_tt): SUM(nmes) != nobs')
    END IF

    IF (ALLOCATED(loadmask)) DEALLOCATE(loadmask)
    ALLOCATE(loadmask(nvar,Q))
    loadmask = loadmask0

    IF (ALLOCATED(nload_item)) DEALLOCATE(nload_item)
    ALLOCATE(nload_item(nvar))

    DO k = 1, nvar
      nload_item(k) = COUNT(loadmask(k,:) == 1)
    END DO

    IF (ANY(nload_item < 1)) THEN
       CALL rexit('ERROR (loglik_long_tt): item without loading')
    END IF

    maxload = MAXVAL(nload_item)
    IF (maxload < 1) maxload = 1

    IF (ALLOCATED(load_idx)) DEALLOCATE(load_idx)
    ALLOCATE(load_idx(maxload, nvar))
    load_idx = 0

    DO k = 1, nvar
      ll = 0
      DO qq = 1, Q
        IF (loadmask(k,qq) == 1) THEN
          ll = ll + 1
          load_idx(ll, k) = qq
        END IF
      END DO
    END DO

    IF (ALLOCATED(nodes_eta))      DEALLOCATE(nodes_eta)
    IF (ALLOCATED(weights_eta))    DEALLOCATE(weights_eta)
    IF (ALLOCATED(logweights_eta)) DEALLOCATE(logweights_eta)

    ALLOCATE(nodes_eta(nquad_eta * dimEta))
    nodes_eta(1:nquad_eta*dimEta) = nodes_eta0(1:nquad_eta*dimEta)

    ALLOCATE(weights_eta(nquad_eta))
    weights_eta(1:nquad_eta) = weights_eta0(1:nquad_eta)

    ALLOCATE(logweights_eta(nquad_eta))
    DO ll = 1, nquad_eta
      logweights_eta(ll) = LOG(weights_eta(ll))
    END DO

    IF (re_item_flag == 1) THEN

      IF (ALLOCATED(nodes_re))      DEALLOCATE(nodes_re)
      IF (ALLOCATED(weights_re))    DEALLOCATE(weights_re)
      IF (ALLOCATED(logweights_re)) DEALLOCATE(logweights_re)

      ALLOCATE(nodes_re(nquad_re))
      nodes_re(1:nquad_re) = nodes_re0(1:nquad_re)

      ALLOCATE(weights_re(nquad_re))
      weights_re(1:nquad_re) = weights_re0(1:nquad_re)

      ALLOCATE(logweights_re(nquad_re))
      DO ll = 1, nquad_re
        logweights_re(ll) = LOG(weights_re(ll))
      END DO

    ELSE

      IF (ALLOCATED(nodes_re))      DEALLOCATE(nodes_re)
      IF (ALLOCATED(weights_re))    DEALLOCATE(weights_re)
      IF (ALLOCATED(logweights_re)) DEALLOCATE(logweights_re)

    END IF

    cache_fit_id     = fit_id0
    cache_nvar       = nvar0
    cache_nsub       = nsub0
    cache_nobs       = nobs0
    cache_n          = n0
    cache_Q          = Q0
    cache_nquad_eta  = nquad_eta0
    cache_nquad_re   = MERGE(nquad_re0, 0, re_item0 == 1)
    cache_ready      = .TRUE.

  END IF

  need_param_alloc = (.NOT. ALLOCATED(lambda)) .OR. &
                     (SIZE(lambda,1) /= nvar) .OR. (SIZE(lambda,2) /= Q) .OR. &
                     (.NOT. ALLOCATED(gamma)) .OR. (SIZE(gamma) /= nvar) .OR. &
                     (.NOT. ALLOCATED(theta_k)) .OR. (SIZE(theta_k) /= nvar) .OR. &
                     (.NOT. ALLOCATED(invsd_k)) .OR. (SIZE(invsd_k) /= nvar) .OR. &
                     (.NOT. ALLOCATED(sigmak)) .OR. (SIZE(sigmak) /= nvar) .OR. &
                     (.NOT. ALLOCATED(sqrt_sigmak)) .OR. (SIZE(sqrt_sigmak) /= nvar) .OR. &
                     (.NOT. ALLOCATED(kappa)) .OR. (SIZE(kappa) /= dimEta) .OR. &
                     (.NOT. ALLOCATED(L_phi)) .OR. (SIZE(L_phi,1) /= dimEta) .OR. &
                     (SIZE(L_phi,2) /= dimEta) .OR. &
                     (.NOT. ALLOCATED(theta_long)) .OR. (SIZE(theta_long,1) /= nvar) .OR. &
                     (SIZE(theta_long,2) /= n) .OR. &
                     (.NOT. ALLOCATED(invsd_long)) .OR. (SIZE(invsd_long,1) /= nvar) .OR. &
                     (SIZE(invsd_long,2) /= n) .OR. &
                     (.NOT. ALLOCATED(tau)) .OR. (SIZE(tau) /= SUM(cm1))

  IF (need_param_alloc) THEN

     IF (ALLOCATED(lambda))       DEALLOCATE(lambda)
     IF (ALLOCATED(gamma))        DEALLOCATE(gamma)
     IF (ALLOCATED(theta_k))      DEALLOCATE(theta_k)
     IF (ALLOCATED(invsd_k))      DEALLOCATE(invsd_k)
     IF (ALLOCATED(sigmak))       DEALLOCATE(sigmak)
     IF (ALLOCATED(sqrt_sigmak))  DEALLOCATE(sqrt_sigmak)
     IF (ALLOCATED(kappa))        DEALLOCATE(kappa)
     IF (ALLOCATED(L_phi))        DEALLOCATE(L_phi)
     IF (ALLOCATED(tau))          DEALLOCATE(tau)
     IF (ALLOCATED(theta_long))   DEALLOCATE(theta_long)
     IF (ALLOCATED(invsd_long))   DEALLOCATE(invsd_long)

     ALLOCATE(lambda(nvar,Q))
     ALLOCATE(gamma(nvar))
     ALLOCATE(theta_k(nvar), invsd_k(nvar))
     ALLOCATE(sigmak(nvar), sqrt_sigmak(nvar))
     ALLOCATE(kappa(dimEta))
     ALLOCATE(L_phi(dimEta,dimEta))
     ALLOCATE(tau(SUM(cm1)))
     ALLOCATE(theta_long(nvar,n))
     ALLOCATE(invsd_long(nvar,n))

  END IF

  IF (Q > 1) THEN
    maxload = MAXVAL(nload_item)
    IF (maxload < 1) maxload = 1

    IF ((.NOT. ALLOCATED(lambda_act)) .OR. &
        (SIZE(lambda_act,1) /= maxload) .OR. &
        (SIZE(lambda_act,2) /= nvar)) THEN
      IF (ALLOCATED(lambda_act)) DEALLOCATE(lambda_act)
      ALLOCATE(lambda_act(maxload, nvar))
    END IF
  END IF

  ! thresholds: identification controlled by cfa_threshold_id_flag
  DO k = 1, nvar
     cm3(k) = cm1(k) - 2
     IF (cm1(k) < 3) THEN
        CALL rexit('ERROR (CFA): threshold identification requires at least 4 response categories for every item')
     END IF
  END DO

  ntre_calc = SUM(cm3)
  IF (ntre_calc /= ntre0) THEN
     CALL rexit('ERROR (CFA): ntre0 mismatch. ntre0')
  END IF
  ntre = ntre0
  IF (ntre0 > 0) ta(1:ntre0) = b0(1:ntre0)

  pos = 0
  j = 1
  tau = 0.d0
  DO k = 1, nvar
     CALL cfa_fill_tau_item(cm1(k), ta(j), monotone_tau0, tau(pos+1:pos+cm1(k)))
     j = j + cm3(k)
     pos = pos + cm1(k)
  END DO

  lambda = 0.d0
  lambda = 0.d0
  pos = ntre0

  DO qq = 1, Q
    ref_item = 0

    DO k = 1, nvar
      IF (loadmask(k,qq) == 1) THEN
        ref_item = k
        EXIT
      END IF
    END DO

    IF (ALLOCATED(cfa_marker_item)) ref_item = cfa_marker_item(qq)
    IF (ref_item == 0 .OR. loadmask(ref_item,qq) /= 1) THEN
      CALL rexit('ERROR (loglik_long_tt): factor')
    END IF

    lambda(ref_item,qq) = 1.d0

    DO k = 1, nvar
      IF (loadmask(k,qq) == 0) THEN
        lambda(k,qq) = 0.d0
      ELSE
        IF (k == ref_item) CYCLE
        pos = pos + 1
        lambda(k,qq) = b0(pos)
      END IF
    END DO
  END DO

  IF (Q > 1) THEN
    lambda_act = 0.d0
    DO k = 1, nvar
      DO ll = 1, nload_item(k)
        qq = load_idx(ll, k)
        lambda_act(ll, k) = lambda(k, qq)
      END DO
    END DO
  END IF

  DO k = 1, nvar
    pos = pos + 1
    gamma(k) = b0(pos)
  END DO

  IF (theta_time_flag == 0) THEN

    DO k = 1, nvar
      IF (cm1(k) == 1) THEN
        theta_k(k) = 1.d0
        invsd_k(k) = 1.d0

        DO t = 1, n
          theta_long(k,t) = 1.d0
          invsd_long(k,t) = 1.d0
        END DO

      ELSE
        pos = pos + 1
        theta_k(k) = EXP(b0(pos))
        invsd_k(k) = 1.d0 / SQRT(theta_k(k))

        DO t = 1, n
          theta_long(k,t) = theta_k(k)
          invsd_long(k,t) = invsd_k(k)
        END DO
      END IF
    END DO

  ELSE

    DO k = 1, nvar
      IF (cm1(k) == 1) THEN
        theta_k(k) = 1.d0
        invsd_k(k) = 1.d0

        DO t = 1, n
          theta_long(k,t) = 1.d0
          invsd_long(k,t) = 1.d0
        END DO

      ELSE
        DO t = 1, n
          pos = pos + 1
          theta_long(k,t) = EXP(b0(pos))
          invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
        END DO

        theta_k(k) = theta_long(k,1)
        invsd_k(k) = invsd_long(k,1)
      END IF
    END DO

  END IF

  IF (re_item_flag == 1) THEN
    DO k = 1, nvar
      pos = pos + 1
      sigmak(k)      = EXP(b0(pos))
      sqrt_sigmak(k) = SQRT(sigmak(k))
    END DO
  ELSE
    DO k = 1, nvar
      sigmak(k)      = 0.d0
      sqrt_sigmak(k) = 0.d0
    END DO
  END IF

  DO qq = 1, Q
    kappa(qq) = 0.d0
  END DO

  IF (n > 1) THEN
    DO t = 2, n
      DO qq = 1, Q
        pos = pos + 1
        kappa((t-1)*Q + qq) = b0(pos)
      END DO
    END DO
  END IF


  ! factor regressions on covariates with time-invariant beta
  CALL ensure_cfa_beta_storage()
  IF (p_fix > 0) THEN
     DO qq = 1, Q
        DO r = 1, p_fix
           pos = pos + 1
           cfa_beta(qq,r) = b0(pos)
        END DO
     END DO
  END IF

  L_phi = 0.d0

  IF (phi_struct_flag == 0) THEN

    DO t = 1, dimEta
      pos = pos + 1
      L_phi(t,t) = EXP(b0(pos))
    END DO

    DO col = 1, dimEta
      DO row = col + 1, dimEta
        pos = pos + 1
        L_phi(row,col) = b0(pos)
      END DO
    END DO

  ELSE IF (phi_struct_flag == 1) THEN

    DO col = 1, Q
      DO row = 1, Q
        pos = pos + 1
        A_tmp(row,col) = b0(pos)
      END DO
    END DO

    DO row = 1, Q
      pos = pos + 1
      sqrt1_tmp(row) = SQRT(EXP(b0(pos)))
    END DO

    CALL build_Lphi_ar1_full(n, Q, A_tmp, sqrt1_tmp, L_phi)

  ELSE IF (phi_struct_flag == 2) THEN

    IF (ALLOCATED(A_ar)) DEALLOCATE(A_ar)
    IF (ALLOCATED(sqrt_sig1)) DEALLOCATE(sqrt_sig1)

    ALLOCATE(A_ar(Q,Q), sqrt_sig1(Q))

    DO col = 1, Q
      DO row = 1, Q
        pos = pos + 1
        A_ar(row,col) = b0(pos)
      END DO
    END DO

    DO row = 1, Q
      pos = pos + 1
      sqrt_sig1(row) = SQRT(EXP(b0(pos)))
    END DO

    L_phi = 0.d0

  ELSE
    CALL rexit('ERROR (loglik_long_tt): phi_struct_flag invalid')
  END IF

  n_theta_est = 0

  IF (theta_time_flag == 0) THEN
    DO k = 1, nvar
      IF (cm1(k) > 1) n_theta_est = n_theta_est + 1
    END DO
  ELSE
    DO k = 1, nvar
      IF (cm1(k) > 1) n_theta_est = n_theta_est + n
    END DO
  END IF

  n_lambda_free = 0
  DO qq = 1, Q
    DO k = 1, nvar
      IF (loadmask(k,qq) == 1) n_lambda_free = n_lambda_free + 1
    END DO
  END DO

  n_lambda_free = n_lambda_free - Q

  IF (phi_struct_flag == 0) THEN
    n_phi_par = dim_eta * (dim_eta + 1) / 2
  ELSE
    n_phi_par = Q * Q + Q
  END IF

  expected_npm = ntre0                             &
               + n_lambda_free                     &
               + nvar                              &
               + n_theta_est                       &
               + MERGE(nvar, 0, re_item_flag == 1) &
               + (n - 1) * Q                       &
               + Q * p_fix                      &
               + n_phi_par

  IF (npm0 /= expected_npm) THEN
    CALL rexit('ERROR (loglik_long_tt): npm0')
  END IF

  IF (pos /= npm0) THEN
    CALL rexit('ERROR (loglik_long_tt): decoded pos')
  END IF

  CALL build_tt_precomp()

  ll_tt = vrais_long_tt()

  vrais = -ll_tt
  ! print *, vrais

END SUBROUTINE loglik_long_tt



DOUBLE PRECISION FUNCTION log_item_block_growth_tt(i, obs_start, k_item, l_eta)
  USE shared
  USE utilities, ONLY: alnorm_fast
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, obs_start, k_item, l_eta

  INTEGER :: j_mes, idx_obs, yijk, t_idx
  INTEGER :: kk0, nmes_k, c_k, l_re, design_i
  DOUBLE PRECISION :: mu_ijk, a_k
  DOUBLE PRECISION :: loglik_re, logterm, log_max, sumexp

  design_i = MERGE(i, 1, growth_subject_specific_design == 1)
  nmes_k = nmes(i, k_item)

  IF (nmes_k == 0) THEN
     log_item_block_growth_tt = 0.d0
     RETURN
  END IF

  kk0 = kk_base(k_item)
  c_k = c(k_item)

  IF (re_item_flag == 0) THEN

     loglik_re = 0.d0
     idx_obs = obs_start - 1

     DO j_mes = 1, nmes_k
        idx_obs = idx_obs + 1

        yijk  = Y(idx_obs)
        t_idx = nitem(idx_obs)

        mu_ijk = growth_tt_mu0_grid(k_item, t_idx, l_eta, design_i)

        loglik_re = loglik_re + log_prob_obs_growth_tt( &
                    mu_ijk, invsd_long(k_item,t_idx),   &
                    yijk, yijk + kk0, c_k)
     END DO

     log_item_block_growth_tt = loglik_re
     RETURN

  END IF

  log_max = -HUGE(1.d0)
  sumexp  = 0.d0

  DO l_re = 1, nquad_re

     a_k = sqrt_sigmak(k_item) * nodes_re(l_re)

     loglik_re = 0.d0
     idx_obs = obs_start - 1

     DO j_mes = 1, nmes_k
        idx_obs = idx_obs + 1

        yijk  = Y(idx_obs)
        t_idx = nitem(idx_obs)

        mu_ijk = growth_tt_mu0_grid(k_item, t_idx, l_eta, design_i) + a_k

        loglik_re = loglik_re + log_prob_obs_growth_tt( &
                    mu_ijk, invsd_long(k_item,t_idx),   &
                    yijk, yijk + kk0, c_k)
     END DO

     logterm = logweights_re(l_re) + loglik_re

     IF (logterm > log_max) THEN
        sumexp  = sumexp * EXP(log_max - logterm) + 1.d0
        log_max = logterm
     ELSE
        sumexp = sumexp + EXP(logterm - log_max)
     END IF

  END DO

  log_item_block_growth_tt = log_max + LOG(sumexp)

CONTAINS

  DOUBLE PRECISION FUNCTION log_prob_obs_growth_tt(mu, invsd, y, tuk, ck)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: mu, invsd
    INTEGER, INTENT(IN) :: y, tuk, ck

    DOUBLE PRECISION, PARAMETER :: LIK_FLOOR = 1.0D-300
    DOUBLE PRECISION :: tausup, tauinf, p

    IF (ck == 2) THEN

       tausup = -mu * invsd

       IF (y == 1) THEN
          p = alnorm_fast(tausup, .FALSE.)
       ELSE
          p = alnorm_fast(tausup, .TRUE.)
       END IF

    ELSE

       IF (y == 1) THEN

          tausup = (tau(tuk) - mu) * invsd
          p = alnorm_fast(tausup, .FALSE.)

       ELSE IF (y == ck) THEN

          tauinf = (tau(tuk-1) - mu) * invsd
          p = alnorm_fast(tauinf, .TRUE.)

       ELSE

          tausup = (tau(tuk)   - mu) * invsd
          tauinf = (tau(tuk-1) - mu) * invsd

          IF (tauinf > 0.d0) THEN
             p = alnorm_fast(tauinf, .TRUE.) - alnorm_fast(tausup, .TRUE.)
          ELSE
             p = alnorm_fast(tausup, .FALSE.) - alnorm_fast(tauinf, .FALSE.)
          END IF

       END IF

    END IF

    IF (p /= p .OR. p <= LIK_FLOOR) THEN
       log_prob_obs_growth_tt = LOG(LIK_FLOOR)
    ELSE
       log_prob_obs_growth_tt = LOG(p)
    END IF

  END FUNCTION log_prob_obs_growth_tt

END FUNCTION log_item_block_growth_tt

DOUBLE PRECISION FUNCTION vrais_growth_tt_i(i, obs0)
  USE shared
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, obs0

  INTEGER :: l_eta, k_item, obs_start
  DOUBLE PRECISION :: loglik_eta, logterm
  DOUBLE PRECISION :: log_max, sumexp

  log_max = -HUGE(1.d0)
  sumexp  = 0.d0

  DO l_eta = 1, nquad_eta

     loglik_eta = 0.d0
     obs_start  = obs0

     DO k_item = 1, nvar

        IF (nmes(i,k_item) > 0) THEN
           loglik_eta = loglik_eta + log_item_block_growth_tt( &
                        i, obs_start, k_item, l_eta)
        END IF

        obs_start = obs_start + nmes(i,k_item)

     END DO

     logterm = logweights_eta(l_eta) + loglik_eta

     IF (logterm > log_max) THEN
        sumexp  = sumexp * EXP(log_max - logterm) + 1.d0
        log_max = logterm
     ELSE
        sumexp = sumexp + EXP(logterm - log_max)
     END IF

  END DO

  vrais_growth_tt_i = log_max + LOG(sumexp)

END FUNCTION vrais_growth_tt_i


DOUBLE PRECISION FUNCTION vrais_growth_tt()
  USE shared
  IMPLICIT NONE

  INTEGER :: i, obs0
  DOUBLE PRECISION :: tmp

  obs0 = 1
  vrais_growth_tt = 0.d0

  DO i = 1, nsub

     tmp = vrais_growth_tt_i(i, obs0)

     vrais_growth_tt = vrais_growth_tt + tmp

     obs0 = obs0 + SUM(nmes(i,:))

  END DO

END FUNCTION vrais_growth_tt

!===========================================================
! Flexible identification/options for second-order latent growth
! ident0: 0=effects, 1=baseline_std, 2=mixed_marker, 3=marker
!===========================================================
SUBROUTINE set_growth_tt_options(ident0, threshold_id0, delta_mode0, &
                                 nvar0, n0, Q0, pfix0, gamma_time_free0, gamma_dep0, &
                                 x_slope_col0, marker_item0, rho_fixed_flag0, rho_fixed0, zero_time0, ok0) &
     BIND(C, name="set_growth_tt_options")
  USE, INTRINSIC :: ISO_C_BINDING
  USE shared
  IMPLICIT NONE
  INTEGER(C_INT), INTENT(IN) :: ident0, threshold_id0, delta_mode0
  INTEGER(C_INT), INTENT(IN) :: nvar0, n0, Q0, pfix0
  INTEGER(C_INT), DIMENSION(nvar0,n0), INTENT(IN) :: gamma_time_free0
  INTEGER(C_INT), DIMENSION(Q0,n0), INTENT(IN) :: gamma_dep0
  INTEGER(C_INT), DIMENSION(pfix0), INTENT(IN) :: x_slope_col0
  INTEGER(C_INT), DIMENSION(Q0), INTENT(IN) :: marker_item0
  INTEGER(C_INT), INTENT(IN) :: rho_fixed_flag0
  REAL(C_DOUBLE), DIMENSION(Q0), INTENT(IN) :: rho_fixed0
  INTEGER(C_INT), DIMENSION(Q0), INTENT(IN) :: zero_time0
  INTEGER(C_INT), INTENT(OUT) :: ok0

  ok0 = 0_C_INT
  IF (nvar0 < 1 .OR. n0 < 2 .OR. Q0 < 1 .OR. pfix0 < 1) RETURN
  IF (ident0 < 0 .OR. ident0 > 3) RETURN
  IF (ANY(zero_time0 < 1) .OR. ANY(zero_time0 > n0)) RETURN
  IF (threshold_id0 < 0 .OR. threshold_id0 > 1) RETURN
  IF (delta_mode0 < 0 .OR. delta_mode0 > 1) RETURN
  IF (rho_fixed_flag0 < 0 .OR. rho_fixed_flag0 > 1) RETURN
  IF (rho_fixed_flag0 == 1 .AND. delta_mode0 /= 1) RETURN
  IF (rho_fixed_flag0 == 1) THEN
     IF (ANY(rho_fixed0 /= rho_fixed0)) RETURN
     IF (ANY(ABS(rho_fixed0) > HUGE(1.d0))) RETURN
  END IF
  IF (ident0 == 2 .OR. ident0 == 3) THEN
     IF (ANY(marker_item0 < 1_C_INT) .OR. ANY(marker_item0 > nvar0)) RETURN
  END IF

  IF (ALLOCATED(growth_gamma_time_free)) DEALLOCATE(growth_gamma_time_free)
  IF (ALLOCATED(growth_gamma_dep)) DEALLOCATE(growth_gamma_dep)
  IF (ALLOCATED(growth_x_slope_col)) DEALLOCATE(growth_x_slope_col)
  IF (ALLOCATED(growth_marker_item)) DEALLOCATE(growth_marker_item)
  IF (ALLOCATED(growth_rho_fixed)) DEALLOCATE(growth_rho_fixed)
  IF (ALLOCATED(growth_zero_time)) DEALLOCATE(growth_zero_time)
  ALLOCATE(growth_gamma_time_free(nvar0,n0), growth_gamma_dep(Q0,n0), growth_x_slope_col(pfix0))
  ALLOCATE(growth_marker_item(Q0), growth_rho_fixed(Q0), growth_zero_time(Q0))
  growth_zero_time = INT(zero_time0)
  growth_gamma_time_free = INT(gamma_time_free0)
  growth_gamma_dep = INT(gamma_dep0)
  growth_x_slope_col = INT(x_slope_col0)
  growth_marker_item = INT(marker_item0)
  growth_gamma_time_free(:,1) = 0
  growth_rho_fixed_flag = INT(rho_fixed_flag0)
  IF (growth_rho_fixed_flag == 1) THEN
     growth_rho_fixed = DBLE(rho_fixed0)
  ELSE
     growth_rho_fixed = 0.d0
  END IF

  growth_ident_flag = INT(ident0)
  growth_threshold_id_flag = INT(threshold_id0)
  growth_delta_mode = INT(delta_mode0)
  growth_options_nvar = nvar0
  growth_options_n = n0
  growth_options_Q = Q0
  growth_options_pfix = pfix0
  growth_options_ready = .TRUE.
  ok0 = 1_C_INT
END SUBROUTINE set_growth_tt_options

SUBROUTINE ensure_growth_tt_options()
  USE shared
  IMPLICIT NONE
  INTEGER :: qidx, k
  INTEGER, ALLOCATABLE :: used(:)

  IF (growth_options_ready .AND. growth_options_nvar == nvar .AND. &
      growth_options_n == n .AND. growth_options_Q == Q .AND. growth_options_pfix == p_fix .AND. &
      ALLOCATED(growth_gamma_time_free) .AND. ALLOCATED(growth_gamma_dep) .AND. &
      ALLOCATED(growth_x_slope_col) .AND. ALLOCATED(growth_marker_item) .AND. &
      ALLOCATED(growth_rho_fixed) .AND. ALLOCATED(growth_zero_time)) RETURN

  IF (ALLOCATED(growth_gamma_time_free)) DEALLOCATE(growth_gamma_time_free)
  IF (ALLOCATED(growth_gamma_dep)) DEALLOCATE(growth_gamma_dep)
  IF (ALLOCATED(growth_x_slope_col)) DEALLOCATE(growth_x_slope_col)
  IF (ALLOCATED(growth_marker_item)) DEALLOCATE(growth_marker_item)
  IF (ALLOCATED(growth_rho_fixed)) DEALLOCATE(growth_rho_fixed)
  IF (ALLOCATED(growth_zero_time)) DEALLOCATE(growth_zero_time)
  ALLOCATE(growth_gamma_time_free(nvar,n), growth_gamma_dep(Q,n), growth_x_slope_col(p_fix), used(nvar))
  ALLOCATE(growth_marker_item(Q), growth_rho_fixed(Q), growth_zero_time(Q))
  growth_zero_time = 1
  growth_gamma_time_free = 0
  growth_gamma_dep = 0
  growth_x_slope_col = 0
  growth_marker_item = 0
  growth_rho_fixed = 0.d0
  growth_rho_fixed_flag = 0
  used = 0
  DO qidx = 1, Q
     DO k = nvar, 1, -1
        IF (loadmask(k,qidx) == 1 .AND. used(k) == 0) THEN
           growth_gamma_dep(qidx,1) = k
           growth_marker_item(qidx) = k
           used(k) = 1
           EXIT
        END IF
     END DO
  END DO
  DEALLOCATE(used)

  growth_ident_flag = 0
  growth_threshold_id_flag = 0
  growth_delta_mode = 0
  growth_options_nvar = nvar
  growth_options_n = n
  growth_options_Q = Q
  growth_options_pfix = p_fix
  growth_options_ready = .TRUE.
END SUBROUTINE ensure_growth_tt_options

SUBROUTINE decode_growth_delta(b0, npm0, pos, delta_free0, delta_g, info)
  USE shared
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npm0, delta_free0
  INTEGER, INTENT(INOUT) :: pos
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(Q), INTENT(OUT) :: delta_g
  INTEGER, INTENT(OUT) :: info
  INTEGER :: qidx

  info = 0
  CALL ensure_growth_tt_options()

  IF (delta_free0 == 1) THEN
     IF (growth_rho_fixed_flag == 1) THEN
        info = -1
        RETURN
     END IF
     IF (pos + Q > npm0) THEN
        info = -2
        RETURN
     END IF
     DO qidx = 1, Q
        pos = pos + 1
        delta_g(qidx) = b0(pos)
     END DO
  ELSE
     IF (growth_delta_mode == 1 .AND. growth_rho_fixed_flag == 1) THEN
        IF (.NOT. ALLOCATED(growth_rho_fixed)) THEN
           info = -3
           RETURN
        END IF
        IF (SIZE(growth_rho_fixed) /= Q) THEN
           info = -4
           RETURN
        END IF
        delta_g = growth_rho_fixed
     ELSE
        delta_g = 0.d0
     END IF
  END IF
END SUBROUTINE decode_growth_delta


INTEGER FUNCTION growth_n_lambda_par()
  USE shared
  IMPLICIT NONE
  INTEGER :: qidx
  CALL ensure_growth_tt_options()
  IF (growth_ident_flag == 1) THEN
     growth_n_lambda_par = COUNT(loadmask == 1)
  ELSE
     growth_n_lambda_par = 0
     DO qidx = 1, Q
        growth_n_lambda_par = growth_n_lambda_par + COUNT(loadmask(:,qidx) == 1) - 1
     END DO
  END IF
END FUNCTION growth_n_lambda_par

INTEGER FUNCTION growth_n_gamma_par()
  USE shared
  IMPLICIT NONE
  INTEGER :: t, m
  CALL ensure_growth_tt_options()
  IF (growth_ident_flag == 1 .OR. growth_ident_flag == 2) THEN
     growth_n_gamma_par = nvar
     IF (n > 1) growth_n_gamma_par = growth_n_gamma_par + COUNT(growth_gamma_time_free(:,2:n) == 1)
  ELSE IF (growth_ident_flag == 3) THEN
     growth_n_gamma_par = nvar - Q
     IF (n > 1) growth_n_gamma_par = growth_n_gamma_par + COUNT(growth_gamma_time_free(:,2:n) == 1)
  ELSE
     growth_n_gamma_par = nvar - Q
     DO t = 2, n
        m = COUNT(growth_gamma_time_free(:,t) == 1)
        IF (m > 0) growth_n_gamma_par = growth_n_gamma_par + MAX(0, m - Q)
     END DO
  END IF
END FUNCTION growth_n_gamma_par

INTEGER FUNCTION growth_n_beta_par()
  USE shared
  IMPLICIT NONE
  CALL ensure_growth_tt_options()
  IF (growth_ident_flag == 1 .OR. growth_ident_flag == 2) THEN
     growth_n_beta_par = Q * MAX(0, p_fix - 1)
  ELSE
     growth_n_beta_par = Q * p_fix
  END IF
END FUNCTION growth_n_beta_par

INTEGER FUNCTION growth_n_zeta_par()
  USE shared
  IMPLICIT NONE
  CALL ensure_growth_tt_options()
  IF (growth_ident_flag == 1) THEN
     growth_n_zeta_par = MAX(0, n - 1) * Q
  ELSE
     growth_n_zeta_par = n * Q
  END IF
END FUNCTION growth_n_zeta_par

DOUBLE PRECISION FUNCTION growth_adjust_score(base, delta, mask)
  USE shared, ONLY: growth_delta_mode
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: base, delta, mask
  IF (mask == 0.d0) THEN
     growth_adjust_score = base
  ELSE IF (growth_delta_mode == 1) THEN
     growth_adjust_score = delta
  ELSE
     growth_adjust_score = base + delta * mask
  END IF
END FUNCTION growth_adjust_score

SUBROUTINE decode_growth_thresholds(b0, npm0, ntre0, monotone_tau0, info)
  USE shared, ONLY: nvar, cm1, tau, kk_base, growth_threshold_id_flag
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npm0, ntre0, monotone_tau0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  INTEGER, INTENT(OUT) :: info
  INTEGER :: k, m, idx, pos_tau, h, j0, j1
  DOUBLE PRECISION :: tmp

  info = 0
  idx = 0
  pos_tau = 0
  DO k = 1, nvar
     m = cm1(k)
     kk_base(k) = pos_tau
     IF (m == 1) THEN
        tau(pos_tau+1) = 0.d0
     ELSE IF (growth_threshold_id_flag == 0) THEN
        tau(pos_tau+1) = 0.d0
        tau(pos_tau+2) = 1.d0
        IF (monotone_tau0 == 1) THEN
           tmp = 1.d0
           DO h = 3, m
              idx = idx + 1
              IF (idx > ntre0 .OR. idx > npm0) THEN
                 info = -100
                 RETURN
              END IF
              tmp = tmp + EXP(b0(idx))
              tau(pos_tau+h) = tmp
           END DO
        ELSE
           DO h = 3, m
              idx = idx + 1
              IF (idx > ntre0 .OR. idx > npm0) THEN
                 info = -100
                 RETURN
              END IF
              tau(pos_tau+h) = b0(idx)
           END DO
        END IF
     ELSE
        j0 = (m + 1) / 2
        j1 = j0 + 1
        tau(pos_tau+j0) = 0.d0
        tau(pos_tau+j1) = 1.d0
        IF (monotone_tau0 == 1) THEN
           tmp = 0.d0
           DO h = j0-1, 1, -1
              idx = idx + 1
              IF (idx > ntre0 .OR. idx > npm0) THEN
                 info = -100
                 RETURN
              END IF
              tmp = tmp - EXP(b0(idx))
              tau(pos_tau+h) = tmp
           END DO
           tmp = 1.d0
           DO h = j1+1, m
              idx = idx + 1
              IF (idx > ntre0 .OR. idx > npm0) THEN
                 info = -100
                 RETURN
              END IF
              tmp = tmp + EXP(b0(idx))
              tau(pos_tau+h) = tmp
           END DO
        ELSE
           DO h = j0-1, 1, -1
              idx = idx + 1
              IF (idx > ntre0 .OR. idx > npm0) THEN
                 info = -100
                 RETURN
              END IF
              tau(pos_tau+h) = b0(idx)
           END DO
           DO h = j1+1, m
              idx = idx + 1
              IF (idx > ntre0 .OR. idx > npm0) THEN
                 info = -100
                 RETURN
              END IF
              tau(pos_tau+h) = b0(idx)
           END DO
        END IF
     END IF
     pos_tau = pos_tau + m
  END DO
  IF (idx /= ntre0) info = -101
END SUBROUTINE decode_growth_thresholds

SUBROUTINE decode_growth_beta(b0, npm0, pos, info)
  USE shared, ONLY: Q, p_fix, beta, growth_ident_flag
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npm0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  INTEGER, INTENT(INOUT) :: pos
  INTEGER, INTENT(OUT) :: info
  INTEGER :: qidx, r
  info = 0
  beta = 0.d0
  DO r = 1, p_fix
     DO qidx = 1, Q
        IF ((growth_ident_flag == 1 .OR. growth_ident_flag == 2) .AND. r == 1) THEN
           beta(qidx,r) = 0.d0
        ELSE
           pos = pos + 1
           IF (pos > npm0) THEN
              info = -100
              RETURN
           END IF
           beta(qidx,r) = b0(pos)
        END IF
     END DO
  END DO
END SUBROUTINE decode_growth_beta

SUBROUTINE decode_growth_random_structure(b0, npm0, pos, dimU, Lu, Sigma_u, zeta_var, info)
  USE shared, ONLY: n, Q, growth_ident_flag, growth_zero_time
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npm0, dimU
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  INTEGER, INTENT(INOUT) :: pos
  DOUBLE PRECISION, DIMENSION(dimU,dimU), INTENT(OUT) :: Lu, Sigma_u
  DOUBLE PRECISION, DIMENSION(n,Q), INTENT(OUT) :: zeta_var
  INTEGER, INTENT(OUT) :: info
  INTEGER :: d, row, col, t, qidx
  DOUBLE PRECISION :: x, prod, pc, eps
  DOUBLE PRECISION, ALLOCATABLE :: Lcorr(:,:), sd(:)

  info = 0
  eps = 1.d-10
  Lu = 0.d0
  Sigma_u = 0.d0
  zeta_var = 0.d0
  ALLOCATE(Lcorr(dimU,dimU), sd(dimU))
  Lcorr = 0.d0
  sd = 0.d0

  IF (growth_ident_flag /= 1) THEN
     DO col = 1, dimU
        DO row = col, dimU
           pos = pos + 1
           IF (pos > npm0) THEN
              info = -100
              DEALLOCATE(Lcorr,sd)
              RETURN
           END IF
           IF (row == col) THEN
              Lu(row,col) = EXP(b0(pos))
           ELSE
              Lu(row,col) = b0(pos)
           END IF
        END DO
     END DO
     Sigma_u = MATMUL(Lu, TRANSPOSE(Lu))
     DO t = 1, n
        DO qidx = 1, Q
           pos = pos + 1
           IF (pos > npm0) THEN
              info = -100
              DEALLOCATE(Lcorr,sd)
              RETURN
           END IF
           zeta_var(t,qidx) = EXP(b0(pos))
        END DO
     END DO
  ELSE
     DO d = 1, dimU
        pos = pos + 1
        IF (pos > npm0) THEN
           info = -100
           DEALLOCATE(Lcorr,sd)
           RETURN
        END IF
        IF (d <= Q) THEN
           x = b0(pos)
           IF (x >= 0.d0) THEN
              sd(d) = SQRT(1.d0 / (1.d0 + EXP(-x)))
           ELSE
              sd(d) = SQRT(EXP(x) / (1.d0 + EXP(x)))
           END IF
        ELSE
           sd(d) = EXP(b0(pos))
        END IF
     END DO

     Lcorr(1,1) = 1.d0
     DO row = 2, dimU
        prod = 1.d0
        DO col = 1, row-1
           pos = pos + 1
           IF (pos > npm0) THEN
              info = -100
              DEALLOCATE(Lcorr,sd)
              RETURN
           END IF
           pc = TANH(b0(pos))
           Lcorr(row,col) = pc * prod
           prod = prod * SQRT(MAX(eps,1.d0-pc*pc))
        END DO
        Lcorr(row,row) = prod
     END DO
     DO row = 1, dimU
        DO col = 1, row
           Lu(row,col) = sd(row) * Lcorr(row,col)
        END DO
     END DO
     Sigma_u = MATMUL(Lu, TRANSPOSE(Lu))

     DO qidx = 1, Q
        zeta_var(growth_zero_time(qidx),qidx) = 1.d0 - Sigma_u(qidx,qidx)
        IF (zeta_var(growth_zero_time(qidx),qidx) <= eps) THEN
           info = -200 - qidx
           DEALLOCATE(Lcorr,sd)
           RETURN
        END IF
     END DO
     DO t = 1, n
        DO qidx = 1, Q
           IF (t == growth_zero_time(qidx)) CYCLE
           pos = pos + 1
           IF (pos > npm0) THEN
              info = -100
              DEALLOCATE(Lcorr,sd)
              RETURN
           END IF
           zeta_var(t,qidx) = EXP(b0(pos))
        END DO
     END DO
  END IF
  DEALLOCATE(Lcorr,sd)
END SUBROUTINE decode_growth_random_structure

!===========================================================
! growth two-tier likelihood, general n and general Q
! replace the old growth-only block with this block.
!
! important conventions:
!   eta index: idx = (t - 1) * Q + q
!   growth random effects index: u0_q = q, u1_q = Q + q
!   xscore(t,q) = dble(t - 1), plus delta(q) only at t = n
!   identification:
!       effects      -> sum of loadings fixed and sum of intercepts fixed
!       baseline_std -> latent intercept mean fixed to 0, baseline latent variance fixed to 1
!       mixed_marker -> selected loading fixed to 1, latent intercept mean fixed to 0
!       marker       -> selected loading fixed to 1, selected item intercept fixed to 0
!   for cross-loadings, the Q x Q dependent-intercept matrix must be nonsingular.
!===========================================================
! paste these inside MODULE core CONTAINS, replacing the existing versions

SUBROUTINE decode_lambda_sum_constraints_growth(b0, npm0, pos, info)
  USE shared, ONLY: nvar, Q, loadmask, lambda, growth_ident_flag, growth_marker_item
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npm0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  INTEGER, INTENT(INOUT) :: pos
  INTEGER, INTENT(OUT) :: info
  INTEGER :: qidx, k, ref_item, nload_q

  CALL ensure_growth_tt_options()
  info = 0
  lambda = 0.d0
  IF (growth_ident_flag == 1) THEN
     DO qidx = 1, Q
        DO k = 1, nvar
           IF (loadmask(k,qidx) == 0) CYCLE
           pos = pos + 1
           IF (pos > npm0) THEN
              info = -100
              RETURN
           END IF
           lambda(k,qidx) = b0(pos)
        END DO
     END DO
     RETURN
  END IF

  IF (growth_ident_flag == 2 .OR. growth_ident_flag == 3) THEN
     IF (.NOT. ALLOCATED(growth_marker_item)) THEN
        info = -20
        RETURN
     END IF
     IF (SIZE(growth_marker_item) /= Q) THEN
        info = -21
        RETURN
     END IF
     DO qidx = 1, Q
        ref_item = growth_marker_item(qidx)
        IF (ref_item < 1 .OR. ref_item > nvar) THEN
           info = -30 - qidx
           RETURN
        END IF
        IF (loadmask(ref_item,qidx) /= 1) THEN
           info = -40 - qidx
           RETURN
        END IF
        lambda(ref_item,qidx) = 1.d0
        DO k = 1, nvar
           IF (loadmask(k,qidx) == 0 .OR. k == ref_item) CYCLE
           pos = pos + 1
           IF (pos > npm0) THEN
              info = -100
              RETURN
           END IF
           lambda(k,qidx) = b0(pos)
        END DO
     END DO
     RETURN
  END IF

  DO qidx = 1, Q
     ref_item = 0
     nload_q = 0
     DO k = 1, nvar
        IF (loadmask(k,qidx) == 1) THEN
           nload_q = nload_q + 1
           ref_item = k
        END IF
     END DO
     IF (ref_item == 0 .OR. nload_q <= 0) THEN
        info = -10 - qidx
        RETURN
     END IF
     lambda(ref_item,qidx) = DBLE(nload_q)
     DO k = 1, nvar
        IF (loadmask(k,qidx) == 0 .OR. k == ref_item) CYCLE
        pos = pos + 1
        IF (pos > npm0) THEN
           info = -100
           RETURN
        END IF
        lambda(k,qidx) = b0(pos)
        lambda(ref_item,qidx) = lambda(ref_item,qidx) - lambda(k,qidx)
     END DO
  END DO
END SUBROUTINE decode_lambda_sum_constraints_growth


SUBROUTINE decode_gamma_sum_constraints_growth(b0, npm0, pos, info)
  USE shared, ONLY: nvar, n, Q, loadmask, gamma, gamma_long, growth_ident_flag, &
                    growth_gamma_time_free, growth_gamma_dep, growth_marker_item
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npm0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  INTEGER, INTENT(INOUT) :: pos
  INTEGER, INTENT(OUT) :: info
  INTEGER :: qidx, k, d, t, m, ref_item
  INTEGER, ALLOCATABLE :: dep_idx(:), is_dep(:), ipiv(:)
  DOUBLE PRECISION, ALLOCATABLE :: Adep(:,:), rhs(:,:), delta(:)

  CALL ensure_growth_tt_options()
  info = 0
  gamma = 0.d0
  IF (ALLOCATED(gamma_long)) DEALLOCATE(gamma_long)
  ALLOCATE(gamma_long(nvar,n))
  gamma_long = 0.d0


  IF (growth_ident_flag == 3) THEN
     IF (.NOT. ALLOCATED(growth_marker_item)) THEN
        info = -20
        RETURN
     END IF
     IF (SIZE(growth_marker_item) /= Q) THEN
        info = -21
        RETURN
     END IF
     ALLOCATE(is_dep(nvar))
     is_dep = 0
     DO qidx = 1, Q
        ref_item = growth_marker_item(qidx)
        IF (ref_item < 1 .OR. ref_item > nvar) THEN
           info = -30 - qidx
           DEALLOCATE(is_dep)
           RETURN
        END IF
        is_dep(ref_item) = 1
     END DO
     DO k = 1, nvar
        IF (is_dep(k) == 0) THEN
           pos = pos + 1
           IF (pos > npm0) THEN
              info = -100
              DEALLOCATE(is_dep)
              RETURN
           END IF
           gamma(k) = b0(pos)
        END IF
     END DO
     gamma_long(:,1) = gamma
     DO t = 2, n
        gamma_long(:,t) = gamma
        DO k = 1, nvar
           IF (growth_gamma_time_free(k,t) == 1) THEN
              pos = pos + 1
              IF (pos > npm0) THEN
                 info = -100
                 DEALLOCATE(is_dep)
                 RETURN
              END IF
              gamma_long(k,t) = b0(pos)
           END IF
        END DO
     END DO
     DEALLOCATE(is_dep)
     RETURN
  END IF

  IF (growth_ident_flag == 1 .OR. growth_ident_flag == 2) THEN
     DO k = 1, nvar
        pos = pos + 1
        IF (pos > npm0) THEN
           info = -100
           RETURN
        END IF
        gamma(k) = b0(pos)
     END DO
     gamma_long(:,1) = gamma
     DO t = 2, n
        gamma_long(:,t) = gamma
        DO k = 1, nvar
           IF (growth_gamma_time_free(k,t) == 1) THEN
              pos = pos + 1
              IF (pos > npm0) THEN
                 info = -100
                 RETURN
              END IF
              gamma_long(k,t) = b0(pos)
           END IF
        END DO
     END DO
     RETURN
  END IF

  ALLOCATE(dep_idx(Q), is_dep(nvar), Adep(Q,Q), rhs(Q,1), ipiv(Q))
  dep_idx = growth_gamma_dep(:,1)
  is_dep = 0
  DO d = 1, Q
     IF (dep_idx(d) < 1 .OR. dep_idx(d) > nvar) THEN
        info = -10-d
        DEALLOCATE(dep_idx,is_dep,Adep,rhs,ipiv)
        RETURN
     END IF
     is_dep(dep_idx(d)) = 1
  END DO
  DO k = 1, nvar
     IF (is_dep(k) == 0) THEN
        pos = pos + 1
        IF (pos > npm0) THEN
           info = -100
           DEALLOCATE(dep_idx,is_dep,Adep,rhs,ipiv)
           RETURN
        END IF
        gamma(k) = b0(pos)
     END IF
  END DO
  Adep = 0.d0
  rhs = 0.d0
  DO qidx = 1, Q
     DO d = 1, Q
        Adep(qidx,d) = DBLE(loadmask(dep_idx(d),qidx))
     END DO
     DO k = 1, nvar
        IF (is_dep(k) == 0) rhs(qidx,1) = rhs(qidx,1) - DBLE(loadmask(k,qidx))*gamma(k)
     END DO
  END DO
  CALL DGESV(Q,1,Adep,Q,ipiv,rhs,Q,info)
  IF (info /= 0) THEN
     DEALLOCATE(dep_idx,is_dep,Adep,rhs,ipiv)
     RETURN
  END IF
  DO d = 1, Q
     gamma(dep_idx(d)) = rhs(d,1)
  END DO
  gamma_long(:,1) = gamma

  ALLOCATE(delta(nvar))
  DO t = 2, n
     gamma_long(:,t) = gamma
     m = COUNT(growth_gamma_time_free(:,t) == 1)
     IF (m == 0) CYCLE
     IF (m <= Q) THEN
        info = -200-t
        DEALLOCATE(dep_idx,is_dep,Adep,rhs,ipiv,delta)
        RETURN
     END IF
     dep_idx = growth_gamma_dep(:,t)
     is_dep = 0
     delta = 0.d0
     DO d = 1, Q
        IF (dep_idx(d) < 1 .OR. dep_idx(d) > nvar .OR. &
            growth_gamma_time_free(dep_idx(d),t) /= 1) THEN
           info = -300-t
           DEALLOCATE(dep_idx,is_dep,Adep,rhs,ipiv,delta)
           RETURN
        END IF
        is_dep(dep_idx(d)) = 1
     END DO
     DO k = 1, nvar
        IF (growth_gamma_time_free(k,t) == 1 .AND. is_dep(k) == 0) THEN
           pos = pos + 1
           IF (pos > npm0) THEN
              info = -100
              DEALLOCATE(dep_idx,is_dep,Adep,rhs,ipiv,delta)
              RETURN
           END IF
           delta(k) = b0(pos)
        END IF
     END DO
     Adep = 0.d0
     rhs = 0.d0
     DO qidx = 1, Q
        DO d = 1, Q
           Adep(qidx,d) = DBLE(loadmask(dep_idx(d),qidx))
        END DO
        DO k = 1, nvar
           IF (growth_gamma_time_free(k,t) == 1 .AND. is_dep(k) == 0) THEN
              rhs(qidx,1) = rhs(qidx,1) - DBLE(loadmask(k,qidx))*delta(k)
           END IF
        END DO
     END DO
     CALL DGESV(Q,1,Adep,Q,ipiv,rhs,Q,info)
     IF (info /= 0) THEN
        DEALLOCATE(dep_idx,is_dep,Adep,rhs,ipiv,delta)
        RETURN
     END IF
     DO d = 1, Q
        delta(dep_idx(d)) = rhs(d,1)
     END DO
     gamma_long(:,t) = gamma + delta
  END DO
  DEALLOCATE(dep_idx,is_dep,Adep,rhs,ipiv,delta)
END SUBROUTINE decode_gamma_sum_constraints_growth


SUBROUTINE loglik_growth_tt(b0, Y0, nvar0, nsub0, nobs0, n0,        &
                            nmes0, c0, cm10, ntre0,               &
                            nquad_eta0, nodes_eta0, weights_eta0, &
                            nquad_re0, nodes_re0, weights_re0,    &
                            nitem0, npm0, re_item0, theta_time0,  &
                            delta_free0, monotone_tau0,           &
                            Q0, loadmask0,                        &
                            p_fix0, p_ran0, X0, Z0, delta_mask0,    &
                            delta_x_col0, delta_z_col0,              &
                            vrais)                                &
        BIND(C, name="loglik_growth_tt")

  USE shared
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0
  INTEGER, INTENT(IN) :: ntre0, npm0
  INTEGER, INTENT(IN) :: nquad_eta0, nquad_re0
  INTEGER, INTENT(IN) :: re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0
  INTEGER, INTENT(IN) :: Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0
  INTEGER, INTENT(IN) :: delta_x_col0, delta_z_col0

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0),       INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0),    INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0

  DOUBLE PRECISION, DIMENSION(npm0),       INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*),          INTENT(IN) :: nodes_eta0
  DOUBLE PRECISION, DIMENSION(nquad_eta0), INTENT(IN) :: weights_eta0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: weights_re0

  DOUBLE PRECISION, INTENT(OUT) :: vrais

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.0D9

  INTEGER :: k, t, r, row, col, pos, jj, pos_tau, info
  INTEGER :: qidx, dimEta, dimU, ntheta_est, expected_npm, ntre_calc, nX, nZ
  INTEGER :: n_lambda_free, n_gamma_free, n_Lu_est
  INTEGER :: nf, h
  DOUBLE PRECISION :: tmp_tau, ll_val

  INTEGER, ALLOCATABLE :: cm3(:)
  DOUBLE PRECISION, ALLOCATABLE :: delta_g(:)
  DOUBLE PRECISION, ALLOCATABLE :: Lu(:,:), Sigma_u(:,:), zeta_var(:,:)

  vrais = BAD_OBJ

  IF (n0 < 2 .OR. Q0 < 1 .OR. ABS(p_fix0) < 1 .OR. p_ran0 < 1) RETURN

  nvar = nvar0
  nsub = nsub0
  nobs = nobs0
  n    = n0
  Q    = Q0
  npm  = npm0

  dimEta = n * Q
  dimU   = Q * p_ran0
  dim_eta = dimEta

  nquad_eta = nquad_eta0
  nquad_re  = nquad_re0

  re_item_flag    = re_item0
  theta_time_flag = theta_time0
  phi_struct_flag = 0
  growth_subject_specific_design = MERGE(1, 0, p_fix0 > 0)
  growth_design_nsub = MERGE(nsub0, 1, growth_subject_specific_design == 1)
  p_fix = ABS(p_fix0)
  p_ran = p_ran0
  growth_delta_x_col = delta_x_col0
  growth_delta_z_col = delta_z_col0

  tt_precomp_ready = .FALSE.

  CALL cleanup_growth_tt_allocations()

  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs))
  ALLOCATE(c(nvar), cm1(nvar))
  ALLOCATE(nodes_eta(nquad_eta*dimEta), weights_eta(nquad_eta), logweights_eta(nquad_eta))
  ALLOCATE(nodes_re(nquad_re), weights_re(nquad_re), logweights_re(nquad_re))
  ALLOCATE(loadmask(nvar,Q))
  ALLOCATE(lambda(nvar,Q), gamma(nvar))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(theta_long(nvar,n), invsd_long(nvar,n))
  ALLOCATE(sigmak(nvar), sqrt_sigmak(nvar))
  ALLOCATE(kappa(dimEta), L_phi(dimEta,dimEta))
  ALLOCATE(kk_base(nvar))
  ALLOCATE(tau(SUM(cm10)))
  ALLOCATE(cm3(nvar))
  ALLOCATE(delta_g(Q))
  ALLOCATE(beta(Q,p_fix), X_eta(growth_design_nsub,n,Q,p_fix), &
           Z_eta(growth_design_nsub,n,Q,p_ran))
  ALLOCATE(growth_delta_mask(growth_design_nsub,n,Q))
  ALLOCATE(Lu(dimU,dimU), Sigma_u(dimU,dimU), zeta_var(n,Q))

  nmes     = nmes0
  Y        = Y0
  nitem    = nitem0
  c        = c0
  cm1      = cm10
  loadmask = loadmask0

  nX = growth_design_nsub * n * Q * p_fix
  nZ = growth_design_nsub * n * Q * p_ran
  X_eta = RESHAPE(X0(1:nX), (/ growth_design_nsub, n, Q, p_fix /))
  Z_eta = RESHAPE(Z0(1:nZ), (/ growth_design_nsub, n, Q, p_ran /))
  growth_delta_mask = RESHAPE(delta_mask0(1:growth_design_nsub*n*Q), &
                              (/ growth_design_nsub, n, Q /))

  nodes_eta(1:nquad_eta*dimEta) = nodes_eta0(1:nquad_eta*dimEta)
  weights_eta(1:nquad_eta) = weights_eta0(1:nquad_eta)

  DO r = 1, nquad_eta
     IF (weights_eta(r) > 0.d0) THEN
        logweights_eta(r) = LOG(weights_eta(r))
     ELSE
        logweights_eta(r) = -HUGE(1.d0)
     END IF
  END DO

  nodes_re(1:nquad_re) = nodes_re0(1:nquad_re)
  weights_re(1:nquad_re) = weights_re0(1:nquad_re)

  DO r = 1, nquad_re
     IF (weights_re(r) > 0.d0) THEN
        logweights_re(r) = LOG(weights_re(r))
     ELSE
        logweights_re(r) = -HUGE(1.d0)
     END IF
  END DO

  CALL setup_load_indices_from_loadmask_growth(info)
  IF (info /= 0) GOTO 999

  DO k = 1, nvar
     cm3(k) = cm1(k) - 2
     IF (cm3(k) < 0) cm3(k) = 0
  END DO

  ntre_calc = SUM(cm3)
  IF (ntre_calc /= ntre0) GOTO 999
  CALL decode_growth_thresholds(b0, npm0, ntre0, monotone_tau0, info)
  IF (info /= 0) GOTO 999

  pos = ntre0

  CALL decode_lambda_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999

  CALL refresh_lambda_act_growth()

  CALL decode_gamma_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999

  theta_k(:) = 1.d0
  invsd_k(:) = 1.d0
  theta_long(:,:) = 1.d0
  invsd_long(:,:) = 1.d0

  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) THEN
           pos = pos + 1
           theta_k(k) = EXP(b0(pos))
           invsd_k(k) = 1.d0 / SQRT(theta_k(k))

           DO t = 1, n
              theta_long(k,t) = theta_k(k)
              invsd_long(k,t) = invsd_k(k)
           END DO
        END IF
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) THEN
           DO t = 1, n
              pos = pos + 1
              theta_long(k,t) = EXP(b0(pos))
              invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
           END DO

           theta_k(k) = theta_long(k,1)
           invsd_k(k) = invsd_long(k,1)
        END IF
     END DO
  END IF
  CALL decode_growth_beta(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999
   CALL decode_growth_delta(b0, npm0, pos, delta_free0, delta_g, info)
   IF (info /= 0) GOTO 999
  CALL decode_growth_random_structure(b0, npm0, pos, dimU, Lu, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 999

  IF (re_item_flag == 1) THEN
     DO k = 1, nvar
        pos = pos + 1
        sigmak(k) = EXP(b0(pos))
        sqrt_sigmak(k) = SQRT(sigmak(k))
     END DO
  ELSE
     sigmak(:) = 0.d0
     sqrt_sigmak(:) = 0.d0
  END IF

  ntheta_est = 0

  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) ntheta_est = ntheta_est + 1
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) ntheta_est = ntheta_est + n
     END DO
  END IF
  n_lambda_free = growth_n_lambda_par()
  n_gamma_free = growth_n_gamma_par()
  n_Lu_est = dimU * (dimU + 1) / 2

  expected_npm = ntre0                                      &
               + n_lambda_free                              &
               + n_gamma_free                               &
               + ntheta_est                                 &
               + growth_n_beta_par()                        &
               + MERGE(Q, 0, delta_free0 == 1)              &
               + n_Lu_est                                   &
               + growth_n_zeta_par()                        &
               + MERGE(nvar, 0, re_item_flag == 1)

  IF (npm0 /= expected_npm) GOTO 999
  IF (pos /= npm0) GOTO 999

  CALL build_growth_eta_moments_design(delta_g, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 999

  CALL build_growth_tt_precomp(info)
  IF (info /= 0) GOTO 999

  ll_val = vrais_growth_tt()

  IF (ll_val /= ll_val .OR. ABS(ll_val) > 1.d300) THEN
     vrais = BAD_OBJ
  ELSE
     vrais = -ll_val
  END IF

999 CONTINUE

  CALL cleanup_growth_tt_allocations()

  IF (ALLOCATED(cm3)) DEALLOCATE(cm3)
  IF (ALLOCATED(delta_g)) DEALLOCATE(delta_g)
  IF (ALLOCATED(Lu)) DEALLOCATE(Lu)
  IF (ALLOCATED(Sigma_u)) DEALLOCATE(Sigma_u)
  IF (ALLOCATED(zeta_var)) DEALLOCATE(zeta_var)

END SUBROUTINE loglik_growth_tt


SUBROUTINE loglik_growth_tt_agh(b0, Y0, nvar0, nsub0, nobs0, n0,        &
                                nmes0, c0, cm10, ntre0,               &
                                nquad_eta0, nodes_eta0, weights_eta0, &
                                nquad_re0, nodes_re0, weights_re0,    &
                                nitem0, npm0, re_item0, theta_time0,  &
                                delta_free0, monotone_tau0,           &
                                update_gh0, fit_id0,                  &
                                Q0, loadmask0,                        &
                            p_fix0, p_ran0, X0, Z0, delta_mask0,    &
                            delta_x_col0, delta_z_col0,              &
                                vrais)                                &
        BIND(C, name="loglik_growth_tt_agh")

  USE shared
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0
  INTEGER, INTENT(IN) :: ntre0, npm0
  INTEGER, INTENT(IN) :: nquad_eta0, nquad_re0
  INTEGER, INTENT(IN) :: re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0
  INTEGER, INTENT(IN) :: update_gh0, fit_id0
  INTEGER, INTENT(IN) :: Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0
  INTEGER, INTENT(IN) :: delta_x_col0, delta_z_col0

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0),       INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0),    INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0

  DOUBLE PRECISION, DIMENSION(npm0),       INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*),          INTENT(IN) :: nodes_eta0
  DOUBLE PRECISION, DIMENSION(nquad_eta0), INTENT(IN) :: weights_eta0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: weights_re0

  DOUBLE PRECISION, INTENT(OUT) :: vrais

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.0D9

  INTEGER :: k, t, r, row, col, pos, jj, pos_tau, info
  INTEGER :: qidx, dimEta, dimU, ntheta_est, expected_npm, ntre_calc, nX, nZ
  INTEGER :: n_lambda_free, n_gamma_free, n_Lu_est
  INTEGER :: nf, h
  DOUBLE PRECISION :: tmp_tau, ll_agh
  LOGICAL :: ok_grid

  INTEGER, ALLOCATABLE :: cm3(:)
  DOUBLE PRECISION, ALLOCATABLE :: delta_g(:)
  DOUBLE PRECISION, ALLOCATABLE :: Lu(:,:), Sigma_u(:,:), zeta_var(:,:)

  vrais = BAD_OBJ

  IF (n0 < 2 .OR. Q0 < 1 .OR. ABS(p_fix0) < 1 .OR. p_ran0 < 1) RETURN

  nvar = nvar0
  nsub = nsub0
  nobs = nobs0
  n    = n0
  Q    = Q0
  npm  = npm0

  dimEta = n * Q
  dimU   = Q * p_ran0
  dim_eta = dimEta

  nquad     = nquad_eta0
  nquad_eta = nquad_eta0
  nquad_re  = nquad_re0

  re_item_flag    = re_item0
  theta_time_flag = theta_time0
  phi_struct_flag = 0
  growth_subject_specific_design = MERGE(1, 0, p_fix0 > 0)
  growth_design_nsub = MERGE(nsub0, 1, growth_subject_specific_design == 1)
  p_fix = ABS(p_fix0)
  p_ran = p_ran0
  growth_delta_x_col = delta_x_col0
  growth_delta_z_col = delta_z_col0

  IF (fit_id0 /= agh_tt_fit_id_cached) THEN
     agh_tt_grid_ready = .FALSE.
     agh_tt_fit_id_cached = fit_id0

     IF (ALLOCATED(agh_tt_z)) DEALLOCATE(agh_tt_z)
     IF (ALLOCATED(agh_tt_modes)) DEALLOCATE(agh_tt_modes)
     IF (ALLOCATED(agh_tt_logconst)) DEALLOCATE(agh_tt_logconst)

     agh_tt_dim_cached = -1
     agh_tt_nquad_cached = -1
     agh_tt_nsub_cached = -1
  END IF

  CALL cleanup_growth_tt_allocations()

  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs))
  ALLOCATE(c(nvar), cm1(nvar))
  ALLOCATE(nodes(nquad*dimEta), weights(nquad))
  ALLOCATE(nodes_re(nquad_re), weights_re(nquad_re), logweights_re(nquad_re))
  ALLOCATE(loadmask(nvar,Q))
  ALLOCATE(lambda(nvar,Q), gamma(nvar))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(theta_long(nvar,n), invsd_long(nvar,n))
  ALLOCATE(sigmak(nvar), sqrt_sigmak(nvar))
  ALLOCATE(kappa(dimEta), L_phi(dimEta,dimEta))
  ALLOCATE(kk_base(nvar))
  ALLOCATE(tau(SUM(cm10)))
  ALLOCATE(cm3(nvar))
  ALLOCATE(delta_g(Q))
  ALLOCATE(beta(Q,p_fix), X_eta(growth_design_nsub,n,Q,p_fix), &
           Z_eta(growth_design_nsub,n,Q,p_ran))
  ALLOCATE(growth_delta_mask(growth_design_nsub,n,Q))
  ALLOCATE(Lu(dimU,dimU), Sigma_u(dimU,dimU), zeta_var(n,Q))

  nmes     = nmes0
  Y        = Y0
  nitem    = nitem0
  c        = c0
  cm1      = cm10
  loadmask = loadmask0

  nX = growth_design_nsub * n * Q * p_fix
  nZ = growth_design_nsub * n * Q * p_ran
  X_eta = RESHAPE(X0(1:nX), (/ growth_design_nsub, n, Q, p_fix /))
  Z_eta = RESHAPE(Z0(1:nZ), (/ growth_design_nsub, n, Q, p_ran /))
  growth_delta_mask = RESHAPE(delta_mask0(1:growth_design_nsub*n*Q), &
                              (/ growth_design_nsub, n, Q /))

  nodes(1:nquad*dimEta) = nodes_eta0(1:nquad*dimEta)
  weights(1:nquad) = weights_eta0(1:nquad)

  nodes_re(1:nquad_re) = nodes_re0(1:nquad_re)
  weights_re(1:nquad_re) = weights_re0(1:nquad_re)

  DO r = 1, nquad_re
     IF (weights_re(r) > 0.d0) THEN
        logweights_re(r) = LOG(weights_re(r))
     ELSE
        logweights_re(r) = -HUGE(1.d0)
     END IF
  END DO

  CALL setup_load_indices_from_loadmask_growth(info)
  IF (info /= 0) GOTO 999

  DO k = 1, nvar
     cm3(k) = cm1(k) - 2
     IF (cm3(k) < 0) cm3(k) = 0
  END DO

  ntre_calc = SUM(cm3)
  IF (ntre_calc /= ntre0) GOTO 999
  CALL decode_growth_thresholds(b0, npm0, ntre0, monotone_tau0, info)
  IF (info /= 0) GOTO 999

  pos = ntre0

  CALL decode_lambda_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999

  CALL refresh_lambda_act_growth()

  CALL decode_gamma_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999

  theta_k(:) = 1.d0
  invsd_k(:) = 1.d0
  theta_long(:,:) = 1.d0
  invsd_long(:,:) = 1.d0

  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) THEN
           pos = pos + 1
           theta_k(k) = EXP(b0(pos))
           invsd_k(k) = 1.d0 / SQRT(theta_k(k))

           DO t = 1, n
              theta_long(k,t) = theta_k(k)
              invsd_long(k,t) = invsd_k(k)
           END DO
        END IF
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) THEN
           DO t = 1, n
              pos = pos + 1
              theta_long(k,t) = EXP(b0(pos))
              invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
           END DO

           theta_k(k) = theta_long(k,1)
           invsd_k(k) = invsd_long(k,1)
        END IF
     END DO
  END IF
  CALL decode_growth_beta(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999
   CALL decode_growth_delta(b0, npm0, pos, delta_free0, delta_g, info)
   IF (info /= 0) GOTO 999
  CALL decode_growth_random_structure(b0, npm0, pos, dimU, Lu, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 999

  IF (re_item_flag == 1) THEN
     DO k = 1, nvar
        pos = pos + 1
        sigmak(k) = EXP(b0(pos))
        sqrt_sigmak(k) = SQRT(sigmak(k))
     END DO
  ELSE
     sigmak(:) = 0.d0
     sqrt_sigmak(:) = 0.d0
  END IF

  ntheta_est = 0

  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) ntheta_est = ntheta_est + 1
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) ntheta_est = ntheta_est + n
     END DO
  END IF
  n_lambda_free = growth_n_lambda_par()
  n_gamma_free = growth_n_gamma_par()
  n_Lu_est = dimU * (dimU + 1) / 2

  expected_npm = ntre0                                      &
               + n_lambda_free                              &
               + n_gamma_free                               &
               + ntheta_est                                 &
               + growth_n_beta_par()                        &
               + MERGE(Q, 0, delta_free0 == 1)              &
               + n_Lu_est                                   &
               + growth_n_zeta_par()                        &
               + MERGE(nvar, 0, re_item_flag == 1)

  IF (npm0 /= expected_npm) GOTO 999
  IF (pos /= npm0) GOTO 999

  CALL build_growth_eta_moments_design(delta_g, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 999

  IF (update_gh0 == 1 .OR. .NOT. agh_tt_grid_ready) THEN
     CALL build_agh_grid_eta_nested(ok_grid)

     IF (.NOT. ok_grid .OR. .NOT. agh_tt_grid_ready) THEN
        vrais = BAD_OBJ
        GOTO 999
     END IF
  END IF

  CALL vrais_long_tt_agh_fixed(ll_agh)

  IF (ll_agh /= ll_agh .OR. ABS(ll_agh) > 1.d300) THEN
     vrais = BAD_OBJ
  ELSE
     vrais = -ll_agh
  END IF

999 CONTINUE

  CALL cleanup_growth_tt_allocations()

  IF (ALLOCATED(cm3)) DEALLOCATE(cm3)
  IF (ALLOCATED(delta_g)) DEALLOCATE(delta_g)
  IF (ALLOCATED(Lu)) DEALLOCATE(Lu)
  IF (ALLOCATED(Sigma_u)) DEALLOCATE(Sigma_u)
  IF (ALLOCATED(zeta_var)) DEALLOCATE(zeta_var)

END SUBROUTINE loglik_growth_tt_agh

SUBROUTINE setup_load_indices_from_loadmask_growth(info)
  USE shared, ONLY: nvar, Q, loadmask, nload_item, load_idx, lambda_act
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: info
  INTEGER :: k, qidx, ll, maxload

  info = 0

  IF (.NOT. ALLOCATED(loadmask)) THEN
     info = -1
     RETURN
  END IF

  IF (ALLOCATED(nload_item)) DEALLOCATE(nload_item)
  ALLOCATE(nload_item(nvar))

  DO k = 1, nvar
     nload_item(k) = COUNT(loadmask(k,1:Q) == 1)
  END DO

  IF (ANY(nload_item < 1)) THEN
     info = -2
     RETURN
  END IF

  maxload = MAXVAL(nload_item)

  IF (ALLOCATED(load_idx)) DEALLOCATE(load_idx)
  ALLOCATE(load_idx(maxload, nvar))
  load_idx = 0

  DO k = 1, nvar
     ll = 0
     DO qidx = 1, Q
        IF (loadmask(k,qidx) == 1) THEN
           ll = ll + 1
           load_idx(ll,k) = qidx
        END IF
     END DO
  END DO

  IF (ALLOCATED(lambda_act)) DEALLOCATE(lambda_act)
  ALLOCATE(lambda_act(maxload, nvar))
  lambda_act = 0.d0

END SUBROUTINE setup_load_indices_from_loadmask_growth


SUBROUTINE refresh_lambda_act_growth()
  USE shared, ONLY: nvar, nload_item, load_idx, lambda, lambda_act
  IMPLICIT NONE

  INTEGER :: k, ll, qidx

  lambda_act = 0.d0

  DO k = 1, nvar
     DO ll = 1, nload_item(k)
        qidx = load_idx(ll,k)
        lambda_act(ll,k) = lambda(k,qidx)
     END DO
  END DO

END SUBROUTINE refresh_lambda_act_growth




SUBROUTINE activate_growth_subject(i)
  USE shared
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  INTEGER :: design_i

  IF (.NOT. ALLOCATED(kappa_sub) .OR. .NOT. ALLOCATED(L_phi_sub)) RETURN
  IF (i < 1 .OR. i > nsub) RETURN

  design_i = MERGE(i, 1, growth_subject_specific_design == 1)
  IF (design_i < 1 .OR. design_i > SIZE(kappa_sub,2)) RETURN

  kappa(:) = kappa_sub(:,design_i)
  L_phi(:,:) = L_phi_sub(:,:,design_i)
END SUBROUTINE activate_growth_subject


SUBROUTINE build_growth_eta_moments_design(delta_g, Sigma_u, zeta_var, info)
  USE shared
  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(Q), INTENT(IN) :: delta_g
  DOUBLE PRECISION, DIMENSION(Q*p_ran,Q*p_ran), INTENT(IN) :: Sigma_u
  DOUBLE PRECISION, DIMENSION(n,Q), INTENT(IN) :: zeta_var
  INTEGER, INTENT(OUT) :: info

  INTEGER :: i, t, qidx, r, row, col, dimEta, dimU, n_design
  DOUBLE PRECISION :: xval, zval
  DOUBLE PRECISION, ALLOCATABLE :: Zbig(:,:), Zref(:,:), tmp(:,:)
  LOGICAL :: same_z

  info = 0
  dimEta = n * Q
  dimU = Q * p_ran
  n_design = growth_design_nsub

  IF (nsub < 1 .OR. n < 2 .OR. Q < 1 .OR. p_fix < 1 .OR. p_ran < 1 .OR. &
      n_design < 1) THEN
     info = -1
     RETURN
  END IF

  IF (.NOT. ALLOCATED(beta) .OR. .NOT. ALLOCATED(X_eta) .OR. &
      .NOT. ALLOCATED(Z_eta) .OR. .NOT. ALLOCATED(growth_delta_mask)) THEN
     info = -2
     RETURN
  END IF

  IF (growth_delta_x_col < 0 .OR. growth_delta_x_col > p_fix .OR. &
      growth_delta_z_col < 0 .OR. growth_delta_z_col > p_ran) THEN
     info = -3
     RETURN
  END IF

  IF (ALLOCATED(kappa_sub)) DEALLOCATE(kappa_sub)
  IF (ALLOCATED(L_phi_sub)) DEALLOCATE(L_phi_sub)
  ALLOCATE(kappa_sub(dimEta,n_design), L_phi_sub(dimEta,dimEta,n_design))
  ALLOCATE(Zbig(dimEta,dimU), Zref(dimEta,dimU), tmp(dimEta,dimU))
  Zref = 0.d0

  DO i = 1, n_design
     Zbig = 0.d0

     DO t = 1, n
        DO qidx = 1, Q
           row = (t - 1) * Q + qidx
           kappa_sub(row,i) = 0.d0

           DO r = 1, p_fix
              xval = X_eta(i,t,qidx,r)
              IF (r == growth_delta_x_col) THEN
                 xval = growth_adjust_score(xval, delta_g(qidx), growth_delta_mask(i,t,qidx))
              ELSE IF (ALLOCATED(growth_x_slope_col) .AND. growth_x_slope_col(r) == 1) THEN
                 zval = X_eta(i,t,qidx,growth_delta_x_col)
                 zval = growth_adjust_score(zval, delta_g(qidx), growth_delta_mask(i,t,qidx))
                 xval = xval * zval
              END IF
              kappa_sub(row,i) = kappa_sub(row,i) + xval * beta(qidx,r)
           END DO

           DO r = 1, p_ran
              zval = Z_eta(i,t,qidx,r)
              IF (r == growth_delta_z_col) THEN
                 zval = growth_adjust_score(zval, delta_g(qidx), growth_delta_mask(i,t,qidx))
              END IF
              col = (r - 1) * Q + qidx
              Zbig(row,col) = zval
           END DO
        END DO
     END DO

     same_z = .FALSE.
     IF (i == 1) THEN
        Zref = Zbig
     ELSE
        same_z = ALL(Zbig == Zref)
     END IF

     IF (same_z) THEN
        L_phi_sub(:,:,i) = L_phi_sub(:,:,1)
        CYCLE
     END IF

     CALL DGEMM('N','N', dimEta, dimU, dimU, 1.d0, Zbig, dimEta, &
                Sigma_u, dimU, 0.d0, tmp, dimEta)
     CALL DGEMM('N','T', dimEta, dimEta, dimU, 1.d0, tmp, dimEta, &
                Zbig, dimEta, 0.d0, L_phi_sub(:,:,i), dimEta)

     DO t = 1, n
        DO qidx = 1, Q
           row = (t - 1) * Q + qidx
           IF (zeta_var(t,qidx) <= 0.d0 .OR. zeta_var(t,qidx) /= zeta_var(t,qidx)) THEN
              info = -4
              DEALLOCATE(Zbig, Zref, tmp)
              RETURN
           END IF
           L_phi_sub(row,row,i) = L_phi_sub(row,row,i) + zeta_var(t,qidx)
        END DO
     END DO

     CALL DPOTRF('L', dimEta, L_phi_sub(:,:,i), dimEta, info)
     IF (info /= 0) THEN
        DEALLOCATE(Zbig, Zref, tmp)
        RETURN
     END IF
     DO col = 2, dimEta
        L_phi_sub(1:col-1,col,i) = 0.d0
     END DO
  END DO

  IF (.NOT. ALLOCATED(kappa)) ALLOCATE(kappa(dimEta))
  IF (.NOT. ALLOCATED(L_phi)) ALLOCATE(L_phi(dimEta,dimEta))
  CALL activate_growth_subject(1)

  DEALLOCATE(Zbig, Zref, tmp)
END SUBROUTINE build_growth_eta_moments_design


SUBROUTINE build_growth_eta_moments(n_in, Q_in, alpha0, alpha1, delta_g, Sigma_u, &
                                    zeta_var, kappa_out, L_eta, info)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n_in, Q_in
  DOUBLE PRECISION, DIMENSION(Q_in), INTENT(IN) :: alpha0, alpha1, delta_g
  DOUBLE PRECISION, DIMENSION(2*Q_in,2*Q_in), INTENT(IN) :: Sigma_u
  DOUBLE PRECISION, DIMENSION(n_in,Q_in), INTENT(IN) :: zeta_var
  DOUBLE PRECISION, DIMENSION(n_in*Q_in), INTENT(OUT) :: kappa_out
  DOUBLE PRECISION, DIMENSION(n_in*Q_in,n_in*Q_in), INTENT(OUT) :: L_eta
  INTEGER, INTENT(OUT) :: info

  INTEGER :: t, q, a, b, idx1, idx2, dimEta, dimU
  DOUBLE PRECISION :: xtq, covtr
  DOUBLE PRECISION, ALLOCATABLE :: Bmat(:,:)

  info = 0
  dimEta = n_in * Q_in
  dimU = 2 * Q_in

  IF (n_in < 2 .OR. Q_in < 1) THEN
     info = -1
     RETURN
  END IF

  ALLOCATE(Bmat(dimEta, dimU))
  Bmat = 0.d0

  DO t = 1, n_in
     DO q = 1, Q_in
        idx1 = (t - 1) * Q_in + q
        xtq = DBLE(t - 1)
        IF (t == n_in) xtq = xtq + delta_g(q)

        kappa_out(idx1) = alpha0(q) + xtq * alpha1(q)
        Bmat(idx1, q) = 1.d0
        Bmat(idx1, Q_in + q) = xtq
     END DO
  END DO

  L_eta = 0.d0

  DO idx1 = 1, dimEta
     DO idx2 = 1, dimEta
        covtr = 0.d0
        DO a = 1, dimU
           DO b = 1, dimU
              covtr = covtr + Bmat(idx1,a) * Sigma_u(a,b) * Bmat(idx2,b)
           END DO
        END DO
        L_eta(idx1,idx2) = covtr
     END DO
  END DO

  DO t = 1, n_in
     DO q = 1, Q_in
        idx1 = (t - 1) * Q_in + q
        IF (zeta_var(t,q) <= 0.d0 .OR. zeta_var(t,q) /= zeta_var(t,q)) THEN
           info = -2
           DEALLOCATE(Bmat)
           RETURN
        END IF
        L_eta(idx1,idx1) = L_eta(idx1,idx1) + zeta_var(t,q)
     END DO
  END DO

  CALL DPOTRF('L', dimEta, L_eta, dimEta, info)

  DEALLOCATE(Bmat)

END SUBROUTINE build_growth_eta_moments


SUBROUTINE build_growth_tt_precomp(info)
  USE shared
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: info

  INTEGER :: i, dimEta, l_eta, baseZ, irow, icol
  INTEGER :: k_item, t, ll, idx_tq, base_eta
  LOGICAL :: need_alloc
  DOUBLE PRECISION :: linpred
  DOUBLE PRECISION, ALLOCATABLE :: z_lat(:), eta_lat(:)

  info = 0
  dimEta = n * Q

  IF (dimEta <= 0 .OR. nsub <= 0 .OR. growth_design_nsub <= 0) THEN
     info = -1
     RETURN
  END IF

  IF (.NOT. ALLOCATED(nload_item) .OR. .NOT. ALLOCATED(load_idx) .OR. &
      .NOT. ALLOCATED(lambda_act) .OR. .NOT. ALLOCATED(kappa_sub) .OR. &
      .NOT. ALLOCATED(L_phi_sub)) THEN
     info = -2
     RETURN
  END IF

  need_alloc = (.NOT. tt_precomp_ready) .OR.                         &
               (.NOT. ALLOCATED(growth_tt_mu0_grid))

  IF (.NOT. need_alloc) THEN
     need_alloc = SIZE(growth_tt_mu0_grid,1) /= nvar .OR.             &
                  SIZE(growth_tt_mu0_grid,2) /= n .OR.                &
                  SIZE(growth_tt_mu0_grid,3) /= nquad_eta .OR.        &
                  SIZE(growth_tt_mu0_grid,4) /= growth_design_nsub
  END IF

  IF (need_alloc) THEN
     IF (ALLOCATED(growth_tt_mu0_grid)) DEALLOCATE(growth_tt_mu0_grid)

     ALLOCATE(growth_tt_mu0_grid(nvar,n,nquad_eta,growth_design_nsub))

     tt_precomp_dimEta = dimEta
     tt_precomp_nquad_eta = nquad_eta
     tt_precomp_nvar = nvar
     tt_precomp_n = n
     tt_precomp_Q = Q
     tt_precomp_ready = .TRUE.
  END IF

  ALLOCATE(z_lat(dimEta), eta_lat(dimEta))

  DO i = 1, growth_design_nsub
     CALL activate_growth_subject(i)

     DO l_eta = 1, nquad_eta
        baseZ = (l_eta - 1) * dimEta
        z_lat = nodes_eta(baseZ+1:baseZ+dimEta)

        DO irow = 1, dimEta
           eta_lat(irow) = kappa(irow)
           DO icol = 1, irow
              eta_lat(irow) = eta_lat(irow) + L_phi(irow,icol) * z_lat(icol)
           END DO
        END DO

        DO t = 1, n
           base_eta = (t - 1) * Q
           DO k_item = 1, nvar
              linpred = 0.d0
              DO ll = 1, nload_item(k_item)
                 idx_tq = base_eta + load_idx(ll,k_item)
                 linpred = linpred + lambda_act(ll,k_item) * eta_lat(idx_tq)
              END DO
              growth_tt_mu0_grid(k_item,t,l_eta,i) = gamma_at(k_item,t) + linpred
           END DO
        END DO
     END DO
  END DO

  CALL activate_growth_subject(1)
  DEALLOCATE(z_lat, eta_lat)
END SUBROUTINE build_growth_tt_precomp



SUBROUTINE cleanup_growth_tt_allocations()
  USE shared
  IMPLICIT NONE

  IF (ALLOCATED(nmes)) DEALLOCATE(nmes)
  IF (ALLOCATED(Y)) DEALLOCATE(Y)
  IF (ALLOCATED(nitem)) DEALLOCATE(nitem)
  IF (ALLOCATED(c)) DEALLOCATE(c)
  IF (ALLOCATED(cm1)) DEALLOCATE(cm1)

  IF (ALLOCATED(nodes_eta)) DEALLOCATE(nodes_eta)
  IF (ALLOCATED(weights_eta)) DEALLOCATE(weights_eta)
  IF (ALLOCATED(logweights_eta)) DEALLOCATE(logweights_eta)

  IF (ALLOCATED(nodes)) DEALLOCATE(nodes)
  IF (ALLOCATED(weights)) DEALLOCATE(weights)

  IF (ALLOCATED(nodes_re)) DEALLOCATE(nodes_re)
  IF (ALLOCATED(weights_re)) DEALLOCATE(weights_re)
  IF (ALLOCATED(logweights_re)) DEALLOCATE(logweights_re)

  IF (ALLOCATED(loadmask)) DEALLOCATE(loadmask)
  IF (ALLOCATED(nload_item)) DEALLOCATE(nload_item)
  IF (ALLOCATED(load_idx)) DEALLOCATE(load_idx)
  IF (ALLOCATED(lambda_act)) DEALLOCATE(lambda_act)

  IF (ALLOCATED(lambda)) DEALLOCATE(lambda)
  IF (ALLOCATED(gamma)) DEALLOCATE(gamma)
  IF (ALLOCATED(gamma_long)) DEALLOCATE(gamma_long)

  IF (ALLOCATED(theta_k)) DEALLOCATE(theta_k)
  IF (ALLOCATED(invsd_k)) DEALLOCATE(invsd_k)
  IF (ALLOCATED(theta_long)) DEALLOCATE(theta_long)
  IF (ALLOCATED(invsd_long)) DEALLOCATE(invsd_long)

  IF (ALLOCATED(sigmak)) DEALLOCATE(sigmak)
  IF (ALLOCATED(sqrt_sigmak)) DEALLOCATE(sqrt_sigmak)

  IF (ALLOCATED(kappa)) DEALLOCATE(kappa)
  IF (ALLOCATED(L_phi)) DEALLOCATE(L_phi)
  IF (ALLOCATED(kappa_sub)) DEALLOCATE(kappa_sub)
  IF (ALLOCATED(L_phi_sub)) DEALLOCATE(L_phi_sub)
  IF (ALLOCATED(beta)) DEALLOCATE(beta)
  IF (ALLOCATED(X_eta)) DEALLOCATE(X_eta)
  IF (ALLOCATED(Z_eta)) DEALLOCATE(Z_eta)
  IF (ALLOCATED(growth_delta_mask)) DEALLOCATE(growth_delta_mask)
  IF (ALLOCATED(tau)) DEALLOCATE(tau)
  IF (ALLOCATED(kk_base)) DEALLOCATE(kk_base)

  IF (ALLOCATED(growth_tt_mu0_grid)) DEALLOCATE(growth_tt_mu0_grid)

  tt_precomp_ready = .FALSE.

END SUBROUTINE cleanup_growth_tt_allocations



!===========================================================
! Individual log-likelihood vectors for growth two-tier models.
! Insert this block inside MODULE core, before END MODULE core.
! It returns log-likelihood contributions ell_i, not negative log-likelihood.
!===========================================================

SUBROUTINE loglik_growth_tt_i_vec(b0, Y0, nvar0, nsub0, nobs0, n0,        &
                            nmes0, c0, cm10, ntre0,               &
                            nquad_eta0, nodes_eta0, weights_eta0, &
                            nquad_re0, nodes_re0, weights_re0,    &
                            nitem0, npm0, re_item0, theta_time0,  &
                            delta_free0, monotone_tau0,           &
                            Q0, loadmask0,                        &
                            p_fix0, p_ran0, X0, Z0, delta_mask0,    &
                            delta_x_col0, delta_z_col0,              &
                            loglik_i)                             &
        BIND(C, name="loglik_growth_tt_i_vec")

  USE shared
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0
  INTEGER, INTENT(IN) :: ntre0, npm0
  INTEGER, INTENT(IN) :: nquad_eta0, nquad_re0
  INTEGER, INTENT(IN) :: re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0
  INTEGER, INTENT(IN) :: Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0
  INTEGER, INTENT(IN) :: delta_x_col0, delta_z_col0

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0),       INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0),    INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0

  DOUBLE PRECISION, DIMENSION(npm0),       INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*),          INTENT(IN) :: nodes_eta0
  DOUBLE PRECISION, DIMENSION(nquad_eta0), INTENT(IN) :: weights_eta0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: weights_re0

  DOUBLE PRECISION, DIMENSION(nsub0), INTENT(OUT) :: loglik_i

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.0D9

  INTEGER :: k, t, r, row, col, pos, jj, pos_tau, info
  INTEGER :: qidx, dimEta, dimU, ntheta_est, expected_npm, ntre_calc, nX, nZ
  INTEGER :: n_lambda_free, n_gamma_free, n_Lu_est
  INTEGER :: nf, h
  DOUBLE PRECISION :: tmp_tau, tmp_i
  INTEGER :: i, obs0

  INTEGER, ALLOCATABLE :: cm3(:)
  DOUBLE PRECISION, ALLOCATABLE :: delta_g(:)
  DOUBLE PRECISION, ALLOCATABLE :: Lu(:,:), Sigma_u(:,:), zeta_var(:,:)

  loglik_i(:) = -BAD_OBJ

  IF (n0 < 2 .OR. Q0 < 1 .OR. ABS(p_fix0) < 1 .OR. p_ran0 < 1) RETURN

  nvar = nvar0
  nsub = nsub0
  nobs = nobs0
  n    = n0
  Q    = Q0
  npm  = npm0

  dimEta = n * Q
  dimU   = Q * p_ran0
  dim_eta = dimEta

  nquad_eta = nquad_eta0
  nquad_re  = nquad_re0

  re_item_flag    = re_item0
  theta_time_flag = theta_time0
  phi_struct_flag = 0
  growth_subject_specific_design = MERGE(1, 0, p_fix0 > 0)
  growth_design_nsub = MERGE(nsub0, 1, growth_subject_specific_design == 1)
  p_fix = ABS(p_fix0)
  p_ran = p_ran0
  growth_delta_x_col = delta_x_col0
  growth_delta_z_col = delta_z_col0

  tt_precomp_ready = .FALSE.

  CALL cleanup_growth_tt_allocations()

  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs))
  ALLOCATE(c(nvar), cm1(nvar))
  ALLOCATE(nodes_eta(nquad_eta*dimEta), weights_eta(nquad_eta), logweights_eta(nquad_eta))
  ALLOCATE(nodes_re(nquad_re), weights_re(nquad_re), logweights_re(nquad_re))
  ALLOCATE(loadmask(nvar,Q))
  ALLOCATE(lambda(nvar,Q), gamma(nvar))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(theta_long(nvar,n), invsd_long(nvar,n))
  ALLOCATE(sigmak(nvar), sqrt_sigmak(nvar))
  ALLOCATE(kappa(dimEta), L_phi(dimEta,dimEta))
  ALLOCATE(kk_base(nvar))
  ALLOCATE(tau(SUM(cm10)))
  ALLOCATE(cm3(nvar))
  ALLOCATE(delta_g(Q))
  ALLOCATE(beta(Q,p_fix), X_eta(growth_design_nsub,n,Q,p_fix), &
           Z_eta(growth_design_nsub,n,Q,p_ran))
  ALLOCATE(growth_delta_mask(growth_design_nsub,n,Q))
  ALLOCATE(Lu(dimU,dimU), Sigma_u(dimU,dimU), zeta_var(n,Q))

  nmes     = nmes0
  Y        = Y0
  nitem    = nitem0
  c        = c0
  cm1      = cm10
  loadmask = loadmask0

  nX = growth_design_nsub * n * Q * p_fix
  nZ = growth_design_nsub * n * Q * p_ran
  X_eta = RESHAPE(X0(1:nX), (/ growth_design_nsub, n, Q, p_fix /))
  Z_eta = RESHAPE(Z0(1:nZ), (/ growth_design_nsub, n, Q, p_ran /))
  growth_delta_mask = RESHAPE(delta_mask0(1:growth_design_nsub*n*Q), &
                              (/ growth_design_nsub, n, Q /))

  nodes_eta(1:nquad_eta*dimEta) = nodes_eta0(1:nquad_eta*dimEta)
  weights_eta(1:nquad_eta) = weights_eta0(1:nquad_eta)

  DO r = 1, nquad_eta
     IF (weights_eta(r) > 0.d0) THEN
        logweights_eta(r) = LOG(weights_eta(r))
     ELSE
        logweights_eta(r) = -HUGE(1.d0)
     END IF
  END DO

  nodes_re(1:nquad_re) = nodes_re0(1:nquad_re)
  weights_re(1:nquad_re) = weights_re0(1:nquad_re)

  DO r = 1, nquad_re
     IF (weights_re(r) > 0.d0) THEN
        logweights_re(r) = LOG(weights_re(r))
     ELSE
        logweights_re(r) = -HUGE(1.d0)
     END IF
  END DO

  CALL setup_load_indices_from_loadmask_growth(info)
  IF (info /= 0) GOTO 999

  DO k = 1, nvar
     cm3(k) = cm1(k) - 2
     IF (cm3(k) < 0) cm3(k) = 0
  END DO

  ntre_calc = SUM(cm3)
  IF (ntre_calc /= ntre0) GOTO 999
  CALL decode_growth_thresholds(b0, npm0, ntre0, monotone_tau0, info)
  IF (info /= 0) GOTO 999

  pos = ntre0

  CALL decode_lambda_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999

  CALL refresh_lambda_act_growth()

  CALL decode_gamma_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999

  theta_k(:) = 1.d0
  invsd_k(:) = 1.d0
  theta_long(:,:) = 1.d0
  invsd_long(:,:) = 1.d0

  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) THEN
           pos = pos + 1
           theta_k(k) = EXP(b0(pos))
           invsd_k(k) = 1.d0 / SQRT(theta_k(k))

           DO t = 1, n
              theta_long(k,t) = theta_k(k)
              invsd_long(k,t) = invsd_k(k)
           END DO
        END IF
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) THEN
           DO t = 1, n
              pos = pos + 1
              theta_long(k,t) = EXP(b0(pos))
              invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
           END DO

           theta_k(k) = theta_long(k,1)
           invsd_k(k) = invsd_long(k,1)
        END IF
     END DO
  END IF
  CALL decode_growth_beta(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999
   CALL decode_growth_delta(b0, npm0, pos, delta_free0, delta_g, info)
   IF (info /= 0) GOTO 999
  CALL decode_growth_random_structure(b0, npm0, pos, dimU, Lu, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 999

  IF (re_item_flag == 1) THEN
     DO k = 1, nvar
        pos = pos + 1
        sigmak(k) = EXP(b0(pos))
        sqrt_sigmak(k) = SQRT(sigmak(k))
     END DO
  ELSE
     sigmak(:) = 0.d0
     sqrt_sigmak(:) = 0.d0
  END IF

  ntheta_est = 0

  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) ntheta_est = ntheta_est + 1
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) ntheta_est = ntheta_est + n
     END DO
  END IF
  n_lambda_free = growth_n_lambda_par()
  n_gamma_free = growth_n_gamma_par()
  n_Lu_est = dimU * (dimU + 1) / 2

  expected_npm = ntre0                                      &
               + n_lambda_free                              &
               + n_gamma_free                               &
               + ntheta_est                                 &
               + growth_n_beta_par()                        &
               + MERGE(Q, 0, delta_free0 == 1)              &
               + n_Lu_est                                   &
               + growth_n_zeta_par()                        &
               + MERGE(nvar, 0, re_item_flag == 1)

  IF (npm0 /= expected_npm) GOTO 999
  IF (pos /= npm0) GOTO 999

  CALL build_growth_eta_moments_design(delta_g, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 999

  CALL build_growth_tt_precomp(info)
  IF (info /= 0) GOTO 999

  obs0 = 1

  DO i = 1, nsub
     tmp_i = vrais_growth_tt_i(i, obs0)

     IF (tmp_i == tmp_i .AND. ABS(tmp_i) <= 1.d300) THEN
        loglik_i(i) = tmp_i
     ELSE
        loglik_i(i) = -BAD_OBJ
     END IF

     obs0 = obs0 + SUM(nmes(i,:))
  END DO

999 CONTINUE

  CALL cleanup_growth_tt_allocations()

  IF (ALLOCATED(cm3)) DEALLOCATE(cm3)
  IF (ALLOCATED(delta_g)) DEALLOCATE(delta_g)
  IF (ALLOCATED(Lu)) DEALLOCATE(Lu)
  IF (ALLOCATED(Sigma_u)) DEALLOCATE(Sigma_u)
  IF (ALLOCATED(zeta_var)) DEALLOCATE(zeta_var)

END SUBROUTINE loglik_growth_tt_i_vec

SUBROUTINE loglik_growth_tt_agh_i_vec(b0, Y0, nvar0, nsub0, nobs0, n0,        &
                                nmes0, c0, cm10, ntre0,               &
                                nquad_eta0, nodes_eta0, weights_eta0, &
                                nquad_re0, nodes_re0, weights_re0,    &
                                nitem0, npm0, re_item0, theta_time0,  &
                                delta_free0, monotone_tau0,           &
                                update_gh0, fit_id0,                  &
                                Q0, loadmask0,                        &
                            p_fix0, p_ran0, X0, Z0, delta_mask0,    &
                            delta_x_col0, delta_z_col0,              &
                                loglik_i)                             &
        BIND(C, name="loglik_growth_tt_agh_i_vec")

  USE shared
  USE datindiv, ONLY: nmescur
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0
  INTEGER, INTENT(IN) :: ntre0, npm0
  INTEGER, INTENT(IN) :: nquad_eta0, nquad_re0
  INTEGER, INTENT(IN) :: re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0
  INTEGER, INTENT(IN) :: update_gh0, fit_id0
  INTEGER, INTENT(IN) :: Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0
  INTEGER, INTENT(IN) :: delta_x_col0, delta_z_col0

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0),       INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0),    INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0

  DOUBLE PRECISION, DIMENSION(npm0),       INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*),          INTENT(IN) :: nodes_eta0
  DOUBLE PRECISION, DIMENSION(nquad_eta0), INTENT(IN) :: weights_eta0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: weights_re0

  DOUBLE PRECISION, DIMENSION(nsub0), INTENT(OUT) :: loglik_i

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.0D9

  INTEGER :: k, t, r, row, col, pos, jj, pos_tau, info
  INTEGER :: qidx, dimEta, dimU, ntheta_est, expected_npm, ntre_calc, nX, nZ
  INTEGER :: n_lambda_free, n_gamma_free, n_Lu_est
  INTEGER :: nf, h
  DOUBLE PRECISION :: tmp_tau, tmp_i
  INTEGER :: i, sm
  LOGICAL :: ok_i
  LOGICAL :: ok_grid

  INTEGER, ALLOCATABLE :: cm3(:)
  DOUBLE PRECISION, ALLOCATABLE :: delta_g(:)
  DOUBLE PRECISION, ALLOCATABLE :: Lu(:,:), Sigma_u(:,:), zeta_var(:,:)

  loglik_i(:) = -BAD_OBJ

  IF (n0 < 2 .OR. Q0 < 1 .OR. ABS(p_fix0) < 1 .OR. p_ran0 < 1) RETURN

  nvar = nvar0
  nsub = nsub0
  nobs = nobs0
  n    = n0
  Q    = Q0
  npm  = npm0

  dimEta = n * Q
  dimU   = Q * p_ran0
  dim_eta = dimEta

  nquad     = nquad_eta0
  nquad_eta = nquad_eta0
  nquad_re  = nquad_re0

  re_item_flag    = re_item0
  theta_time_flag = theta_time0
  phi_struct_flag = 0
  growth_subject_specific_design = MERGE(1, 0, p_fix0 > 0)
  growth_design_nsub = MERGE(nsub0, 1, growth_subject_specific_design == 1)
  p_fix = ABS(p_fix0)
  p_ran = p_ran0
  growth_delta_x_col = delta_x_col0
  growth_delta_z_col = delta_z_col0

  IF (fit_id0 /= agh_tt_fit_id_cached) THEN
     agh_tt_grid_ready = .FALSE.
     agh_tt_fit_id_cached = fit_id0

     IF (ALLOCATED(agh_tt_z)) DEALLOCATE(agh_tt_z)
     IF (ALLOCATED(agh_tt_modes)) DEALLOCATE(agh_tt_modes)
     IF (ALLOCATED(agh_tt_logconst)) DEALLOCATE(agh_tt_logconst)

     agh_tt_dim_cached = -1
     agh_tt_nquad_cached = -1
     agh_tt_nsub_cached = -1
  END IF

  CALL cleanup_growth_tt_allocations()

  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs))
  ALLOCATE(c(nvar), cm1(nvar))
  ALLOCATE(nodes(nquad*dimEta), weights(nquad))
  ALLOCATE(nodes_re(nquad_re), weights_re(nquad_re), logweights_re(nquad_re))
  ALLOCATE(loadmask(nvar,Q))
  ALLOCATE(lambda(nvar,Q), gamma(nvar))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(theta_long(nvar,n), invsd_long(nvar,n))
  ALLOCATE(sigmak(nvar), sqrt_sigmak(nvar))
  ALLOCATE(kappa(dimEta), L_phi(dimEta,dimEta))
  ALLOCATE(kk_base(nvar))
  ALLOCATE(tau(SUM(cm10)))
  ALLOCATE(cm3(nvar))
  ALLOCATE(delta_g(Q))
  ALLOCATE(beta(Q,p_fix), X_eta(growth_design_nsub,n,Q,p_fix), &
           Z_eta(growth_design_nsub,n,Q,p_ran))
  ALLOCATE(growth_delta_mask(growth_design_nsub,n,Q))
  ALLOCATE(Lu(dimU,dimU), Sigma_u(dimU,dimU), zeta_var(n,Q))

  nmes     = nmes0
  Y        = Y0
  nitem    = nitem0
  c        = c0
  cm1      = cm10
  loadmask = loadmask0

  nX = growth_design_nsub * n * Q * p_fix
  nZ = growth_design_nsub * n * Q * p_ran
  X_eta = RESHAPE(X0(1:nX), (/ growth_design_nsub, n, Q, p_fix /))
  Z_eta = RESHAPE(Z0(1:nZ), (/ growth_design_nsub, n, Q, p_ran /))
  growth_delta_mask = RESHAPE(delta_mask0(1:growth_design_nsub*n*Q), &
                              (/ growth_design_nsub, n, Q /))

  nodes(1:nquad*dimEta) = nodes_eta0(1:nquad*dimEta)
  weights(1:nquad) = weights_eta0(1:nquad)

  nodes_re(1:nquad_re) = nodes_re0(1:nquad_re)
  weights_re(1:nquad_re) = weights_re0(1:nquad_re)

  DO r = 1, nquad_re
     IF (weights_re(r) > 0.d0) THEN
        logweights_re(r) = LOG(weights_re(r))
     ELSE
        logweights_re(r) = -HUGE(1.d0)
     END IF
  END DO

  CALL setup_load_indices_from_loadmask_growth(info)
  IF (info /= 0) GOTO 999

  DO k = 1, nvar
     cm3(k) = cm1(k) - 2
     IF (cm3(k) < 0) cm3(k) = 0
  END DO

  ntre_calc = SUM(cm3)
  IF (ntre_calc /= ntre0) GOTO 999
  CALL decode_growth_thresholds(b0, npm0, ntre0, monotone_tau0, info)
  IF (info /= 0) GOTO 999

  pos = ntre0

  CALL decode_lambda_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999

  CALL refresh_lambda_act_growth()

  CALL decode_gamma_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999

  theta_k(:) = 1.d0
  invsd_k(:) = 1.d0
  theta_long(:,:) = 1.d0
  invsd_long(:,:) = 1.d0

  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) THEN
           pos = pos + 1
           theta_k(k) = EXP(b0(pos))
           invsd_k(k) = 1.d0 / SQRT(theta_k(k))

           DO t = 1, n
              theta_long(k,t) = theta_k(k)
              invsd_long(k,t) = invsd_k(k)
           END DO
        END IF
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) THEN
           DO t = 1, n
              pos = pos + 1
              theta_long(k,t) = EXP(b0(pos))
              invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
           END DO

           theta_k(k) = theta_long(k,1)
           invsd_k(k) = invsd_long(k,1)
        END IF
     END DO
  END IF
  CALL decode_growth_beta(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999
   CALL decode_growth_delta(b0, npm0, pos, delta_free0, delta_g, info)
   IF (info /= 0) GOTO 999
  CALL decode_growth_random_structure(b0, npm0, pos, dimU, Lu, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 999

  IF (re_item_flag == 1) THEN
     DO k = 1, nvar
        pos = pos + 1
        sigmak(k) = EXP(b0(pos))
        sqrt_sigmak(k) = SQRT(sigmak(k))
     END DO
  ELSE
     sigmak(:) = 0.d0
     sqrt_sigmak(:) = 0.d0
  END IF

  ntheta_est = 0

  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) ntheta_est = ntheta_est + 1
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) ntheta_est = ntheta_est + n
     END DO
  END IF
  n_lambda_free = growth_n_lambda_par()
  n_gamma_free = growth_n_gamma_par()
  n_Lu_est = dimU * (dimU + 1) / 2

  expected_npm = ntre0                                      &
               + n_lambda_free                              &
               + n_gamma_free                               &
               + ntheta_est                                 &
               + growth_n_beta_par()                        &
               + MERGE(Q, 0, delta_free0 == 1)              &
               + n_Lu_est                                   &
               + growth_n_zeta_par()                        &
               + MERGE(nvar, 0, re_item_flag == 1)

  IF (npm0 /= expected_npm) GOTO 999
  IF (pos /= npm0) GOTO 999

  CALL build_growth_eta_moments_design(delta_g, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 999

  IF (update_gh0 == 1 .OR. .NOT. agh_tt_grid_ready) THEN
     CALL build_agh_grid_eta_nested(ok_grid)

     IF (.NOT. ok_grid .OR. .NOT. agh_tt_grid_ready) THEN
        GOTO 999
        GOTO 999
     END IF
  END IF

  nmescur = 0

  DO i = 1, nsub
     sm = SUM(nmes(i,:))

     CALL vrais_long_i_tt_agh_fixed(i, sm, tmp_i, ok_i)

     IF (ok_i .AND. tmp_i == tmp_i .AND. ABS(tmp_i) <= 1.d300) THEN
        loglik_i(i) = tmp_i
     ELSE
        loglik_i(i) = -BAD_OBJ
     END IF

     nmescur = nmescur + sm
  END DO

999 CONTINUE

  CALL cleanup_growth_tt_allocations()

  IF (ALLOCATED(cm3)) DEALLOCATE(cm3)
  IF (ALLOCATED(delta_g)) DEALLOCATE(delta_g)
  IF (ALLOCATED(Lu)) DEALLOCATE(Lu)
  IF (ALLOCATED(Sigma_u)) DEALLOCATE(Sigma_u)
  IF (ALLOCATED(zeta_var)) DEALLOCATE(zeta_var)

END SUBROUTINE loglik_growth_tt_agh_i_vec




SUBROUTINE get_agh_tt_grid_info(ready0, dim0, nquad0, nsub0, fit_id0) &
     BIND(C, name="get_agh_tt_grid_info")
  USE, INTRINSIC :: ISO_C_BINDING
  USE shared
  IMPLICIT NONE

  INTEGER(C_INT), INTENT(OUT) :: ready0, dim0, nquad0, nsub0, fit_id0

  IF (agh_tt_grid_ready .AND. ALLOCATED(agh_tt_z) .AND. &
      ALLOCATED(agh_tt_modes) .AND. ALLOCATED(agh_tt_logconst)) THEN
     ready0 = 1_C_INT
     dim0 = INT(agh_tt_dim_cached, C_INT)
     nquad0 = INT(agh_tt_nquad_cached, C_INT)
     nsub0 = INT(agh_tt_nsub_cached, C_INT)
     fit_id0 = INT(agh_tt_fit_id_cached, C_INT)
  ELSE
     ready0 = 0_C_INT
     dim0 = 0_C_INT
     nquad0 = 0_C_INT
     nsub0 = 0_C_INT
     fit_id0 = INT(agh_tt_fit_id_cached, C_INT)
  END IF
END SUBROUTINE get_agh_tt_grid_info


SUBROUTINE get_agh_tt_grid_copy(z_out, modes_out, logconst_out, &
                                ready0, dim0, nquad0, nsub0, fit_id0) &
     BIND(C, name="get_agh_tt_grid_copy")
  USE, INTRINSIC :: ISO_C_BINDING
  USE shared
  IMPLICIT NONE

  REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT) :: z_out
  REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT) :: modes_out
  REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT) :: logconst_out
  INTEGER(C_INT), INTENT(OUT) :: ready0, dim0, nquad0, nsub0, fit_id0

  INTEGER :: d, l, i, idx

  IF (.NOT. agh_tt_grid_ready .OR. .NOT. ALLOCATED(agh_tt_z) .OR. &
      .NOT. ALLOCATED(agh_tt_modes) .OR. .NOT. ALLOCATED(agh_tt_logconst)) THEN
     ready0 = 0_C_INT
     dim0 = 0_C_INT
     nquad0 = 0_C_INT
     nsub0 = 0_C_INT
     fit_id0 = INT(agh_tt_fit_id_cached, C_INT)
     RETURN
  END IF

  ready0 = 1_C_INT
  dim0 = INT(agh_tt_dim_cached, C_INT)
  nquad0 = INT(agh_tt_nquad_cached, C_INT)
  nsub0 = INT(agh_tt_nsub_cached, C_INT)
  fit_id0 = INT(agh_tt_fit_id_cached, C_INT)

  idx = 0
  DO i = 1, agh_tt_nsub_cached
     DO l = 1, agh_tt_nquad_cached
        DO d = 1, agh_tt_dim_cached
           idx = idx + 1
           z_out(idx) = agh_tt_z(d,l,i)
        END DO
     END DO
  END DO

  idx = 0
  DO i = 1, agh_tt_nsub_cached
     DO d = 1, agh_tt_dim_cached
        idx = idx + 1
        modes_out(idx) = agh_tt_modes(d,i)
     END DO
  END DO

  DO i = 1, agh_tt_nsub_cached
     logconst_out(i) = agh_tt_logconst(i)
  END DO
END SUBROUTINE get_agh_tt_grid_copy


SUBROUTINE set_agh_tt_grid_copy(z_in, modes_in, logconst_in, &
                                dim0, nquad0, nsub0, fit_id0, ready0) &
     BIND(C, name="set_agh_tt_grid_copy")
  USE, INTRINSIC :: ISO_C_BINDING
  USE shared
  IMPLICIT NONE

  REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: z_in
  REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: modes_in
  REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: logconst_in
  INTEGER(C_INT), INTENT(IN) :: dim0, nquad0, nsub0, fit_id0
  INTEGER(C_INT), INTENT(OUT) :: ready0

  INTEGER :: d, l, i, idx
  INTEGER :: dimF, nquadF, nsubF

  dimF = INT(dim0)
  nquadF = INT(nquad0)
  nsubF = INT(nsub0)

  ready0 = 0_C_INT

  IF (dimF <= 0 .OR. nquadF <= 0 .OR. nsubF <= 0) RETURN

  IF (ALLOCATED(agh_tt_z)) DEALLOCATE(agh_tt_z)
  IF (ALLOCATED(agh_tt_modes)) DEALLOCATE(agh_tt_modes)
  IF (ALLOCATED(agh_tt_logconst)) DEALLOCATE(agh_tt_logconst)

  ALLOCATE(agh_tt_z(dimF, nquadF, nsubF))
  ALLOCATE(agh_tt_modes(dimF, nsubF))
  ALLOCATE(agh_tt_logconst(nsubF))

  idx = 0
  DO i = 1, nsubF
     DO l = 1, nquadF
        DO d = 1, dimF
           idx = idx + 1
           agh_tt_z(d,l,i) = z_in(idx)
        END DO
     END DO
  END DO

  idx = 0
  DO i = 1, nsubF
     DO d = 1, dimF
        idx = idx + 1
        agh_tt_modes(d,i) = modes_in(idx)
     END DO
  END DO

  DO i = 1, nsubF
     agh_tt_logconst(i) = logconst_in(i)
  END DO

  agh_tt_dim_cached = dimF
  agh_tt_nquad_cached = nquadF
  agh_tt_nsub_cached = nsubF
  agh_tt_fit_id_cached = INT(fit_id0)
  agh_tt_grid_ready = .TRUE.

  ready0 = 1_C_INT
END SUBROUTINE set_agh_tt_grid_copy




!===========================================================
! Wu & Estabrook (2016) longitudinal ordered-CFA extension
! exact stage-specific identification for the AGH full-equivalent likelihood
!
! wu_stage0:
!   0 = configural ordinal baseline: Condition 7 separately at each time (delta parametrization)
!   1 = threshold invariance, Condition 15
!   2 = threshold + loading invariance, Condition 19 / Proposition 7
!   3 = threshold + loading + intercept invariance, Condition 27 / Proposition 11
!   4 = full measurement invariance: stage 3 + unique/residual variance invariance
!
! exact Wu-Estabrook mode assumes:
!   - ordered polytomous items with at least 3 categories
!   - clustered/simple loading structure: exactly one factor per item
!   - no item random effect
!   - unstructured joint factor covariance across times
!
! baseline stage 0 identification at every time t:
!   kappa_t = 0, nu_t = 0, diag(Phi_tt) = I, diag(Sigma_y*_tt) = I
!   all finite thresholds are free and time-specific
!   with simple structure, theta_kt = 1 - lambda_kt^2
!===========================================================

DOUBLE PRECISION FUNCTION wu_invlogit(x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: x

  IF (x >= 0.d0) THEN
     wu_invlogit = 1.d0 / (1.d0 + EXP(-x))
  ELSE
     wu_invlogit = EXP(x) / (1.d0 + EXP(x))
  END IF
END FUNCTION wu_invlogit


SUBROUTINE build_Lphi_wu2016(b0, npm0, pos, n0, Q0, wu_stage0, Lout, info)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: npm0, n0, Q0, wu_stage0
  INTEGER, INTENT(INOUT) :: pos
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(n0*Q0,n0*Q0), INTENT(OUT) :: Lout
  INTEGER, INTENT(OUT) :: info

  INTEGER :: d, row, col, idx
  DOUBLE PRECISION :: rho, prod, rem
  DOUBLE PRECISION, ALLOCATABLE :: Lcorr(:,:), sdev(:)

  info = 0
  d = n0 * Q0
  Lout = 0.d0

  IF (d <= 0) THEN
     info = -1
     RETURN
  END IF

  ALLOCATE(Lcorr(d,d), sdev(d))
  Lcorr = 0.d0
  sdev = 1.d0

  Lcorr(1,1) = 1.d0

  ! unconstrained partial correlations -> valid correlation Cholesky
  DO row = 2, d
     prod = 1.d0
     DO col = 1, row-1
        pos = pos + 1
        IF (pos > npm0) THEN
           info = -2
           DEALLOCATE(Lcorr, sdev)
           RETURN
        END IF

        rho = TANH(b0(pos))
        rem = MAX(1.d-14, 1.d0 - rho*rho)
        Lcorr(row,col) = rho * prod
        prod = prod * SQRT(rem)
     END DO
     Lcorr(row,row) = prod
  END DO

  ! stages 0/1: all factor variances fixed to one
  ! stages 2/3/4: only time-1 factor variances fixed to one
  IF (wu_stage0 >= 2) THEN
     DO idx = Q0 + 1, d
        pos = pos + 1
        IF (pos > npm0) THEN
           info = -3
           DEALLOCATE(Lcorr, sdev)
           RETURN
        END IF
        sdev(idx) = EXP(b0(pos))
        IF (sdev(idx) <= 0.d0 .OR. sdev(idx) /= sdev(idx)) THEN
           info = -4
           DEALLOCATE(Lcorr, sdev)
           RETURN
        END IF
     END DO
  END IF

  DO row = 1, d
     DO col = 1, row
        Lout(row,col) = sdev(row) * Lcorr(row,col)
     END DO
  END DO

  DEALLOCATE(Lcorr, sdev)
END SUBROUTINE build_Lphi_wu2016


SUBROUTINE loglik_long_tt_agh_full_equiv_wu(b0, Y0, nvar0, nsub0, nobs0, n0, &
                              nmes0, c0, cm10, ntre0,                         &
                              nquad_eta0, nitem0,                             &
                              nodes_eta0, weights_eta0,                       &
                              nquad_re0, nodes_re0, weights_re0,              &
                              npm0, re_item0,                                 &
                              Q0, loadmask0, phi_struct0,                     &
                              monotone_tau0, wu_stage0,                       &
                              update_gh0, fit_id0, vrais)                     &
        BIND(C, name="loglik_long_tt_agh_full_equiv_wu")

  USE shared
  USE utilities
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0
  INTEGER, INTENT(IN) :: ntre0, nquad_eta0, nquad_re0, npm0
  INTEGER, INTENT(IN) :: re_item0, Q0, phi_struct0
  INTEGER, INTENT(IN) :: monotone_tau0, wu_stage0, update_gh0, fit_id0

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0), INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0), INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0

  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_eta0
  DOUBLE PRECISION, DIMENSION(nquad_eta0), INTENT(IN) :: weights_eta0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: weights_re0
  DOUBLE PRECISION, INTENT(OUT) :: vrais

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.d9
  DOUBLE PRECISION, PARAMETER :: LAM_EPS = 1.d-8

  INTEGER :: k, t, qq, pos, j, h, nf, refq
  INTEGER :: dimEta, dimZ, ntre_calc, ntre_base, n_corr, n_sd
  INTEGER :: n_lambda_est, n_gamma_est, n_theta_est, n_kappa_est
  INTEGER :: expected_npm, info, nload
  DOUBLE PRECISION :: prev, lam, ll_agh
  LOGICAL :: ok_grid
  INTEGER, ALLOCATABLE :: cm3(:)
  DOUBLE PRECISION, ALLOCATABLE :: ta(:)

  vrais = BAD_OBJ

  IF (wu_stage0 < 0 .OR. wu_stage0 > 4) RETURN
  IF (re_item0 /= 0) THEN
     RETURN
  END IF
  IF (phi_struct0 /= 0) THEN
     RETURN
  END IF
  IF (n0 < 2 .OR. Q0 < 1 .OR. nvar0 < 1) RETURN
  IF (ANY(c0 < 3)) THEN
     RETURN
  END IF

  ! exact propositions assume clustered/simple loadings
  DO k = 1, nvar0
     nload = COUNT(loadmask0(k,1:Q0) == 1)
     IF (nload /= 1) THEN
        RETURN
     END IF
  END DO
  DO qq = 1, Q0
     IF (COUNT(loadmask0(:,qq) == 1) < 1) THEN
        RETURN
     END IF
  END DO

  nvar = nvar0
  nsub = nsub0
  nobs = nobs0
  n = n0
  Q = Q0
  nquad = nquad_eta0
  nquad_re = nquad_re0
  npm = npm0
  dimEta = n * Q
  dim_eta = dimEta
  phi_struct_flag = 0
  re_item_flag = 0
  theta_time_flag = 1
  p_fix = 0
  p_ran = 0

  SELECT CASE (wu_stage0)
  CASE (0)
     tau_invar_flag = 0
     gamma_invar_flag = 0
     lambda_invar_flag = 0
     theta_invar_flag = 0
  CASE (1)
     tau_invar_flag = 1
     gamma_invar_flag = 0
     lambda_invar_flag = 0
     theta_invar_flag = 0
  CASE (2)
     tau_invar_flag = 1
     gamma_invar_flag = 0
     lambda_invar_flag = 1
     theta_invar_flag = 0
  CASE (3)
     tau_invar_flag = 1
     gamma_invar_flag = 1
     lambda_invar_flag = 1
     theta_invar_flag = 0
  CASE (4)
     tau_invar_flag = 1
     gamma_invar_flag = 1
     lambda_invar_flag = 1
     theta_invar_flag = 1
  END SELECT

  IF (agh_tt_fe_fit_id_cached /= fit_id0) THEN
     agh_tt_fe_grid_ready = .FALSE.
     agh_tt_fe_fit_id_cached = fit_id0
  END IF

  ! data
  IF (ALLOCATED(nmes)) DEALLOCATE(nmes)
  IF (ALLOCATED(Y)) DEALLOCATE(Y)
  IF (ALLOCATED(nitem)) DEALLOCATE(nitem)
  IF (ALLOCATED(c)) DEALLOCATE(c)
  IF (ALLOCATED(cm1)) DEALLOCATE(cm1)
  IF (ALLOCATED(loadmask)) DEALLOCATE(loadmask)

  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs), c(nvar), cm1(nvar))
  ALLOCATE(loadmask(nvar,Q))
  nmes = nmes0
  Y = Y0
  nitem = nitem0
  c = c0
  cm1 = cm10
  loadmask = loadmask0

  IF (nobs <= 0 .OR. SUM(nmes) <= 0) RETURN
  IF (MINVAL(nitem) < 1 .OR. MAXVAL(nitem) > n) RETURN

  ! quadrature
  dimZ = dimEta
  IF (ALLOCATED(nodes)) DEALLOCATE(nodes)
  IF (ALLOCATED(weights)) DEALLOCATE(weights)
  IF (ALLOCATED(nodes_re)) DEALLOCATE(nodes_re)
  IF (ALLOCATED(weights_re)) DEALLOCATE(weights_re)

  ALLOCATE(nodes(nquad*dimEta), weights(nquad))
  nodes(1:nquad*dimEta) = nodes_eta0(1:nquad*dimEta)
  weights(1:nquad) = weights_eta0(1:nquad)

  ALLOCATE(nodes_re(nquad_re), weights_re(nquad_re))
  nodes_re(1:nquad_re) = nodes_re0(1:nquad_re)
  weights_re(1:nquad_re) = weights_re0(1:nquad_re)

  ! parameters
  IF (ALLOCATED(lambda)) DEALLOCATE(lambda)
  IF (ALLOCATED(lambda_long)) DEALLOCATE(lambda_long)
  IF (ALLOCATED(gamma)) DEALLOCATE(gamma)
  IF (ALLOCATED(gamma_long)) DEALLOCATE(gamma_long)
  IF (ALLOCATED(theta_k)) DEALLOCATE(theta_k)
  IF (ALLOCATED(invsd_k)) DEALLOCATE(invsd_k)
  IF (ALLOCATED(theta_long)) DEALLOCATE(theta_long)
  IF (ALLOCATED(invsd_long)) DEALLOCATE(invsd_long)
  IF (ALLOCATED(sigmak)) DEALLOCATE(sigmak)
  IF (ALLOCATED(sqrt_sigmak)) DEALLOCATE(sqrt_sigmak)
  IF (ALLOCATED(kappa)) DEALLOCATE(kappa)
  IF (ALLOCATED(L_phi)) DEALLOCATE(L_phi)
  IF (ALLOCATED(tau)) DEALLOCATE(tau)
  IF (ALLOCATED(tau_long)) DEALLOCATE(tau_long)

  ALLOCATE(lambda(nvar,Q), lambda_long(nvar,Q,n))
  ALLOCATE(gamma(nvar), gamma_long(nvar,n))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(theta_long(nvar,n), invsd_long(nvar,n))
  ALLOCATE(sigmak(nvar), sqrt_sigmak(nvar))
  ALLOCATE(kappa(dimEta), L_phi(dimEta,dimEta))
  ALLOCATE(tau(SUM(cm1)), tau_long(SUM(cm1),n))
  ALLOCATE(cm3(nvar))

  DO k = 1, nvar
     cm3(k) = MAX(cm1(k)-2, 0)
  END DO
  ntre_base = SUM(cm3)

  IF (wu_stage0 == 0) THEN
     ntre_calc = n * SUM(cm1)
  ELSE
     ntre_calc = SUM(cm1)
  END IF

  IF (ntre0 /= ntre_calc) THEN
     RETURN
  END IF

  ALLOCATE(ta(MAX(1,ntre0)))
  IF (ntre0 > 0) ta(1:ntre0) = b0(1:ntre0)

  tau = 0.d0
  tau_long = 0.d0
  j = 1

  IF (wu_stage0 == 0) THEN
     ! configural baseline, Condition 7 separately at each time:
     ! all finite thresholds are free and time-specific
     DO t = 1, n
        pos = 0
        DO k = 1, nvar
           nf = cm1(k)
           IF (monotone_tau0 == 1) THEN
              tau_long(pos+1,t) = ta(j)
              prev = ta(j)
              j = j + 1
              DO h = 2, nf
                 prev = prev + EXP(ta(j))
                 tau_long(pos+h,t) = prev
                 j = j + 1
              END DO
           ELSE
              tau_long(pos+1:pos+nf,t) = ta(j:j+nf-1)
              j = j + nf
           END IF
           pos = pos + nf
        END DO
     END DO
  ELSE
     ! Conditions 15/19/27 and full invariance: all finite thresholds are common and free
     pos = 0
     DO k = 1, nvar
        nf = cm1(k)
        IF (monotone_tau0 == 1) THEN
           tau_long(pos+1,:) = ta(j)
           prev = ta(j)
           j = j + 1
           DO h = 2, nf
              prev = prev + EXP(ta(j))
              tau_long(pos+h,:) = prev
              j = j + 1
           END DO
        ELSE
           DO t = 1, n
              tau_long(pos+1:pos+nf,t) = ta(j:j+nf-1)
           END DO
           j = j + nf
        END IF
        pos = pos + nf
     END DO
  END IF
  tau = tau_long(:,1)

  pos = ntre0

  ! loadings: all nonzero loadings are estimated; no marker loading is used in Wu-Estabrook stages
  lambda = 0.d0
  lambda_long = 0.d0

  IF (wu_stage0 == 0) THEN
     DO t = 1, n
        DO k = 1, nvar
           refq = 0
           DO qq = 1, Q
              IF (loadmask(k,qq) == 1) THEN
                 refq = qq
                 EXIT
              END IF
           END DO
           pos = pos + 1
           IF (pos > npm0) GOTO 900
           ! Condition 7 has diag(Phi_tt)=I and diag(Sigma_y*_tt)=I.
           ! Under simple structure theta_kt=1-lambda_kt^2, hence 0<lambda_kt<1.
           lam = wu_invlogit(b0(pos))
           lam = MIN(1.d0-LAM_EPS, MAX(LAM_EPS, lam))
           lambda_long(k,refq,t) = lam
        END DO
     END DO

  ELSE IF (wu_stage0 == 1) THEN
     DO t = 1, n
        DO k = 1, nvar
           refq = 0
           DO qq = 1, Q
              IF (loadmask(k,qq) == 1) THEN
                 refq = qq
                 EXIT
              END IF
           END DO
           pos = pos + 1
           IF (pos > npm0) GOTO 900
           IF (t == 1) THEN
              lam = wu_invlogit(b0(pos))
              lam = MIN(1.d0-LAM_EPS, MAX(LAM_EPS, lam))
           ELSE
              lam = EXP(b0(pos))
           END IF
           lambda_long(k,refq,t) = lam
        END DO
     END DO

  ELSE
     DO k = 1, nvar
        refq = 0
        DO qq = 1, Q
           IF (loadmask(k,qq) == 1) THEN
              refq = qq
              EXIT
           END IF
        END DO
        pos = pos + 1
        IF (pos > npm0) GOTO 900
        lam = wu_invlogit(b0(pos))
        lam = MIN(1.d0-LAM_EPS, MAX(LAM_EPS, lam))
        DO t = 1, n
           lambda_long(k,refq,t) = lam
        END DO
     END DO
  END IF
  lambda = lambda_long(:,:,1)

  ! item intercepts nu
  gamma = 0.d0
  gamma_long = 0.d0
  IF (wu_stage0 == 0) THEN
     ! Condition 7 at each time fixes all item intercepts nu_t to zero.
     gamma_long = 0.d0
  ELSE IF (wu_stage0 <= 2) THEN
     gamma_long(:,1) = 0.d0
     DO t = 2, n
        DO k = 1, nvar
           pos = pos + 1
           IF (pos > npm0) GOTO 900
           gamma_long(k,t) = b0(pos)
        END DO
     END DO
  ELSE
     gamma_long = 0.d0
  END IF
  gamma = gamma_long(:,1)

  ! unique/residual variances Theta
  theta_long = 1.d0
  invsd_long = 1.d0

  IF (wu_stage0 == 0) THEN
     ! Condition 7 at each time fixes diag(Sigma_y*_tt)=I and diag(Phi_tt)=I.
     ! With one loading per item this determines theta_kt = 1 - lambda_kt^2.
     DO t = 1, n
        DO k = 1, nvar
           refq = 0
           DO qq = 1, Q
              IF (loadmask(k,qq) == 1) THEN
                 refq = qq
                 EXIT
              END IF
           END DO
           theta_long(k,t) = 1.d0 - lambda_long(k,refq,t)**2
           IF (theta_long(k,t) <= 1.d-10 .OR. theta_long(k,t) /= theta_long(k,t)) GOTO 900
        END DO
     END DO

  ELSE
     ! diag(Sigma_1)=I and diag(Phi_11)=I imply theta_k1 = 1-lambda_k1^2 for simple structure
     DO k = 1, nvar
        refq = 0
        DO qq = 1, Q
           IF (loadmask(k,qq) == 1) THEN
              refq = qq
              EXIT
           END IF
        END DO
        theta_long(k,1) = 1.d0 - lambda_long(k,refq,1)**2
        IF (theta_long(k,1) <= 1.d-10 .OR. theta_long(k,1) /= theta_long(k,1)) GOTO 900
     END DO

     IF (wu_stage0 <= 3) THEN
        DO t = 2, n
           DO k = 1, nvar
              pos = pos + 1
              IF (pos > npm0) GOTO 900
              theta_long(k,t) = EXP(b0(pos))
              IF (theta_long(k,t) <= 0.d0 .OR. theta_long(k,t) /= theta_long(k,t)) GOTO 900
           END DO
        END DO
     ELSE
        DO t = 2, n
           theta_long(:,t) = theta_long(:,1)
        END DO
     END IF
  END IF

  DO k = 1, nvar
     theta_k(k) = theta_long(k,1)
     invsd_k(k) = 1.d0 / SQRT(theta_k(k))
     DO t = 1, n
        invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
     END DO
  END DO

  sigmak = 0.d0
  sqrt_sigmak = 0.d0

  ! factor means kappa
  kappa = 0.d0
  IF (wu_stage0 >= 3) THEN
     DO t = 2, n
        DO qq = 1, Q
           pos = pos + 1
           IF (pos > npm0) GOTO 900
           kappa((t-1)*Q + qq) = b0(pos)
        END DO
     END DO
  END IF

  ! unstructured joint factor covariance with stage-specific variance constraints
  CALL build_Lphi_wu2016(b0, npm0, pos, n, Q, wu_stage0, L_phi, info)
  IF (info /= 0) GOTO 900

  ! parameter count
  n_corr = dimEta * (dimEta - 1) / 2
  IF (wu_stage0 >= 2) THEN
     n_sd = dimEta - Q
  ELSE
     n_sd = 0
  END IF

  SELECT CASE (wu_stage0)
  CASE (0)
     n_lambda_est = n*nvar
     n_gamma_est = 0
     n_theta_est = 0
     n_kappa_est = 0
  CASE (1)
     n_lambda_est = n*nvar
     n_gamma_est = (n-1)*nvar
     n_theta_est = (n-1)*nvar
     n_kappa_est = 0
  CASE (2)
     n_lambda_est = nvar
     n_gamma_est = (n-1)*nvar
     n_theta_est = (n-1)*nvar
     n_kappa_est = 0
  CASE (3)
     n_lambda_est = nvar
     n_gamma_est = 0
     n_theta_est = (n-1)*nvar
     n_kappa_est = (n-1)*Q
  CASE (4)
     n_lambda_est = nvar
     n_gamma_est = 0
     n_theta_est = 0
     n_kappa_est = (n-1)*Q
  END SELECT

  expected_npm = ntre_calc + n_lambda_est + n_gamma_est + n_theta_est + &
                 n_kappa_est + n_corr + n_sd

  IF (npm0 /= expected_npm .OR. pos /= npm0) THEN
     GOTO 900
  END IF

  IF (update_gh0 == 1 .OR. .NOT. agh_tt_fe_grid_ready) THEN
     CALL build_agh_grid_full_equiv_tt(ok_grid)
     IF (.NOT. ok_grid .OR. .NOT. agh_tt_fe_grid_ready) GOTO 900
  END IF

  CALL vrais_long_tt_agh_full_equiv_fixed(ll_agh)
  IF (ll_agh /= ll_agh .OR. ABS(ll_agh) > 1.d300) GOTO 900

  vrais = -ll_agh
  DEALLOCATE(cm3, ta)
  RETURN

900 CONTINUE
  vrais = BAD_OBJ
  IF (ALLOCATED(cm3)) DEALLOCATE(cm3)
  IF (ALLOCATED(ta)) DEALLOCATE(ta)
  RETURN

END SUBROUTINE loglik_long_tt_agh_full_equiv_wu

END MODULE core














! PROGRAM debug_one_iter
!   USE core
!   IMPLICIT NONE

!   INTEGER, PARAMETER :: nvar0       = 4      ! 4 item
!   INTEGER, PARAMETER :: nsub0       = 1      ! 1 soggetto
!   INTEGER, PARAMETER :: Q0          = 2      ! 2 fattori
!   INTEGER, PARAMETER :: n0          = 4      ! 4 occasioni (t=1..4)
!   INTEGER, PARAMETER :: nquad0      = 1      ! 1 draw QMC
!   INTEGER, PARAMETER :: p_fix0      = 1      ! solo slope (t-1)
!   INTEGER, PARAMETER :: p_ran0      = 2      ! random intercept + random slope
!   INTEGER, PARAMETER :: ident_type0 = 1      ! 0 = anchor λ=1; 1 = unit-variance (σ_eta^2=1)
!   INTEGER, PARAMETER :: re_item0 = 1
!   INTEGER :: ntre0, npm0, nobs0
!   INTEGER :: i, k, j, dimZ, idx
!   DOUBLE PRECISION :: vrais

!   INTEGER :: nmes0(nsub0,nvar0)
!   INTEGER :: c0(nvar0), cm10(nvar0), fac_item0(nvar0)
!   INTEGER, ALLOCATABLE :: Y0(:), nitem0(:)
!   DOUBLE PRECISION, ALLOCATABLE :: nodes0(:), weights0(:), b0(:)
!   DOUBLE PRECISION :: X0(nsub0,n0,Q0,p_fix0)
!   DOUBLE PRECISION :: Z0(nsub0,n0,Q0,p_ran0)

!   !-------------------------
!   ! Struttura categorie
!   !-------------------------
!   c0   = 4                  ! 4 categorie per ogni item
!   cm10 = c0 - 1             ! C-1 soglie
!   ntre0 = SUM(cm10 - 2)     ! ntre0 = 1 soglia libera per item -> 4

!   ! item 1,2 -> fattore 1 ; item 3,4 -> fattore 2
!   fac_item0 = (/1,1,2,2/)

!   !-------------------------
!   ! nmes, Y, nitem (ordine i,k,j)
!   !-------------------------
!   nmes0(1,:) = (/2,2,2,4/)   ! numero misure per item (soggetto 1)

!   nobs0 = SUM(nmes0)
!   ALLOCATE(Y0(nobs0), nitem0(nobs0))

!   idx = 0
!   DO i = 1, nsub0
!      DO k = 1, nvar0
!         DO j = 1, nmes0(i,k)
!            idx = idx + 1
!            nitem0(idx) = j          ! tempo = j

!            SELECT CASE (k)
!            CASE (1)
!               Y0(idx) = j           ! 1,2
!            CASE (2)
!               Y0(idx) = 3 - j       ! 2,1
!            CASE (3)
!               Y0(idx) = j + 1       ! 2,3
!            CASE (4)
!               Y0(idx) = j           ! 1,2
!            CASE DEFAULT
!               WRITE(*,*) 'ERRORE: k fuori range in costruzione Y0, k=', k
!               STOP
!            END SELECT
!         END DO
!      END DO
!   END DO

!   !-------------------------
!   ! Nodi QMC (dummy, giusto per test)
!   !-------------------------
!   dimZ = Q0*p_ran0 + nvar0         ! 2*2 + 4 = 8
!   ALLOCATE(nodes0(nquad0*dimZ), weights0(nquad0))

!   nodes0 = (/ 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.7d0 /)
!   weights0(1) = 1.d0

!   !-------------------------
!   ! Design X e Z: X = (t-1), Z = (1, t-1)
!   !-------------------------
!   DO i = 1, nsub0
!      DO j = 1, n0
!         DO k = 1, Q0
!            X0(i,j,k,1) = DBLE(j-1)

!            Z0(i,j,k,1) = 1.d0
!            Z0(i,j,k,2) = DBLE(j-1)
!         END DO
!      END DO
!   END DO

!   !-------------------------
!   ! Parametri b0 nel formato di loglik
!   ! ident_type0 = 1 (unit-variance):
!   !   - tutti i lambda(k,q) liberi (per item che caricano su q)
!   !   - sigma_eta2(q) fissate a 1 dentro loglik (nessun parametro in b0)
!   !
!   ! NPM = ntre
!   !       + n_lambda (qui = nvar0 per simple structure)
!   !       + 3*nvar0  (gamma, log theta_k, log sigmak)
!   !       + Q0*p_fix0
!   !       + Q0*p_ran0*(p_ran0+1)/2
!   !
!   !     = ntre0 + 4*nvar0 + Q0*p_fix0 + Q0*p_ran0*(p_ran0+1)/2
!   !-------------------------
!   npm0 = ntre0 + 4*nvar0 + Q0*p_fix0 + Q0*p_ran0*(p_ran0+1)/2
!   ALLOCATE(b0(npm0))

!   idx = 0

!   ! (1) soglie libere: per ciascun item 1 soglia libera (la terza)
!   idx = idx + 1 ; b0(idx) = 1.1d0    ! item 1
!   idx = idx + 1 ; b0(idx) = 1.2d0    ! item 2
!   idx = idx + 1 ; b0(idx) = 1.3d0    ! item 3
!   idx = idx + 1 ; b0(idx) = 1.4d0    ! item 4

!   ! (2) loadings λ(k,q) liberi per tutti gli item (simple structure)
!   ! ordine in loglik (ident_type0=1):
!   !   q=1: k=1,2 con fac_item0(k)=1
!   !   q=2: k=3,4 con fac_item0(k)=2
!   idx = idx + 1 ; b0(idx) = 1.3d0    ! lambda(1,1)
!   idx = idx + 1 ; b0(idx) = 0.7d0    ! lambda(2,1)
!   idx = idx + 1 ; b0(idx) = 1.4d0    ! lambda(3,2)
!   idx = idx + 1 ; b0(idx) = 0.7d0    ! lambda(4,2)

!   ! (3) gamma(k) intercette item
!   idx = idx + 1 ; b0(idx) = 0.10d0
!   idx = idx + 1 ; b0(idx) = 0.20d0
!   idx = idx + 1 ; b0(idx) = 0.30d0
!   idx = idx + 1 ; b0(idx) = 0.40d0

!   ! (4) log theta_k(k) varianza errore di misura
!   idx = idx + 1 ; b0(idx) = LOG(0.5d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.6d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.7d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.8d0)

!   ! (5) beta(q,1) (p_fix0 = 1)
!   idx = idx + 1 ; b0(idx) = 0.20d0       ! beta(1,1)
!   idx = idx + 1 ; b0(idx) = 0.10d0       ! beta(2,1)

!   ! (6) L_chol per ogni fattore q: diag log, poi sotto-diagonale (p_ran0 = 2)
!   ! fattore 1
!   idx = idx + 1 ; b0(idx) = LOG(0.8d0)   ! L11(1)
!   idx = idx + 1 ; b0(idx) = LOG(0.9d0)   ! L22(1)
!   idx = idx + 1 ; b0(idx) = 0.10d0       ! L21(1)

!   ! fattore 2
!   idx = idx + 1 ; b0(idx) = LOG(0.7d0)   ! L11(2)
!   idx = idx + 1 ; b0(idx) = LOG(0.6d0)   ! L22(2)
!   idx = idx + 1 ; b0(idx) = -0.05d0      ! L21(2)

!   ! (7) log sigmak(k) (varianze random sugli item)
!   idx = idx + 1 ; b0(idx) = LOG(0.3d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.4d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.5d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.6d0)

!   IF (idx /= npm0) THEN
!      PRINT *, 'ERRORE: idx != npm0', idx, npm0
!   END IF

!   !-------------------------
!   ! Chiamata singola a loglik
!   !-------------------------
!   CALL loglik(b0        = b0,    &
!               Y0        = Y0,    &
!               nvar0     = nvar0, &
!               nsub0     = nsub0, &
!               nobs0     = nobs0, &
!               Q0        = Q0,    &
!               n0        = n0,    &
!               nmes0     = nmes0, &
!               c0        = c0,    &
!               cm10      = cm10,  &
!               ntre0     = ntre0, &
!               nquad0    = nquad0,&
!               nitem0    = nitem0,&
!               nodes0    = nodes0,&
!               weights0  = weights0,&
!               npm0      = npm0,  &
!               p_fix0    = p_fix0,&
!               p_ran0    = p_ran0,&
!               X0        = X0,    &
!               Z0        = Z0,    &
!               fac_item0 = fac_item0, &
!               ident_type0 = ident_type0, &
!               re_item0 = re_item0, &
!               vrais     = vrais)

!   PRINT *, '*** -loglikelihood (unit-variance) =', vrais

!   DEALLOCATE(b0, Y0, nitem0, nodes0, weights0)

! END PROGRAM debug_one_iter


! PROGRAM debug_one_iter
!   USE core
!   IMPLICIT NONE

!   INTEGER, PARAMETER :: nvar0       = 4      ! 4 item
!   INTEGER, PARAMETER :: nsub0       = 1      ! 1 soggetto
!   INTEGER, PARAMETER :: Q0          = 2      ! 2 fattori
!   INTEGER, PARAMETER :: n0          = 4      ! 4 occasioni (t=1..4)
!   INTEGER, PARAMETER :: p_fix0      = 1      ! solo slope (t-1)
!   INTEGER, PARAMETER :: p_ran0      = 2      ! random intercept + random slope
!   INTEGER, PARAMETER :: ident_type0 = 1      ! 0 = anchor λ=1; 1 = unit-variance (σ_eta^2=1)
!   INTEGER, PARAMETER :: re_item0    = 1      ! 1 = random sugli item

!   ! Gauss–Hermite 1D per normale standard: 3 punti
!   INTEGER, PARAMETER :: nGH1d       = 3

!   INTEGER :: ntre0, npm0, nobs0
!   INTEGER :: i, k, j, dimZ, idx
!   INTEGER :: nquad0
!   INTEGER :: l, d, tmp, idx1d

!   DOUBLE PRECISION :: vrais, w_tot

!   INTEGER :: nmes0(nsub0,nvar0)
!   INTEGER :: c0(nvar0), cm10(nvar0), fac_item0(nvar0)
!   INTEGER, ALLOCATABLE :: Y0(:), nitem0(:)
!   DOUBLE PRECISION, ALLOCATABLE :: nodes0(:), weights0(:), b0(:)
!   DOUBLE PRECISION :: X0(nsub0,n0,Q0,p_fix0)
!   DOUBLE PRECISION :: Z0(nsub0,n0,Q0,p_ran0)

!   DOUBLE PRECISION, DIMENSION(nGH1d) :: x1, w1

!   !-------------------------
!   ! Struttura categorie
!   !-------------------------
!   c0   = 4                  ! 4 categorie per ogni item
!   cm10 = c0 - 1             ! C-1 soglie
!   ntre0 = SUM(cm10 - 2)     ! ntre0 = 1 soglia libera per item -> 4

!   ! item 1,2 -> fattore 1 ; item 3,4 -> fattore 2
!   fac_item0 = (/1,1,2,2/)

!   !-------------------------
!   ! nmes, Y, nitem (ordine i,k,j)
!   !-------------------------
!   nmes0(1,:) = (/2,2,2,4/)   ! numero misure per item (soggetto 1)

!   nobs0 = SUM(nmes0)
!   ALLOCATE(Y0(nobs0), nitem0(nobs0))

!   idx = 0
!   DO i = 1, nsub0
!      DO k = 1, nvar0
!         DO j = 1, nmes0(i,k)
!            idx = idx + 1
!            nitem0(idx) = j          ! tempo = j

!            SELECT CASE (k)
!            CASE (1)
!               Y0(idx) = j           ! 1,2
!            CASE (2)
!               Y0(idx) = 3 - j       ! 2,1
!            CASE (3)
!               Y0(idx) = j + 1       ! 2,3
!            CASE (4)
!               Y0(idx) = j           ! 1,2
!            CASE DEFAULT
!               WRITE(*,*) 'ERRORE: k fuori range in costruzione Y0, k=', k
!               STOP
!            END SELECT
!         END DO
!      END DO
!   END DO

!   !-------------------------
!   ! Dimensione totale random
!   ! dimZ = Q*p_ran + nvar (perché re_item0=1)
!   !-------------------------
!   dimZ = Q0*p_ran0 + nvar0

!   !-------------------------
!   ! Gauss–Hermite 1D per NORMALE STANDARD
!   ! nodi = (-sqrt(3), 0, +sqrt(3))
!   ! pesi = (1/6, 4/6, 1/6)
!   !-------------------------
!   x1(1) = -1.7320508075688772d0
!   x1(2) =  0.d0
!   x1(3) =  1.7320508075688772d0

!   w1(1) = 0.1666666666666667d0
!   w1(2) = 0.6666666666666667d0
!   w1(3) = 0.1666666666666667d0

!   !-------------------------
!   ! Numero totale di nodi GH multidimensionale:
!   ! nquad0 = nGH1d^dimZ (qui = 3^8 = 6561)
!   !-------------------------
!   nquad0 = 1
!   DO d = 1, dimZ
!      nquad0 = nquad0 * nGH1d
!   END DO

!   ALLOCATE(nodes0(nquad0*dimZ), weights0(nquad0))

!   !-------------------------
!   ! Costruzione griglia tensoriale:
!   ! per ogni l = 0..nquad0-1, si rappresenta l in base 3
!   ! e si prende il prodotto dei pesi 1D
!   !-------------------------
!   DO l = 0, nquad0-1
!      tmp   = l
!      w_tot = 1.d0

!      DO d = 1, dimZ
!         idx1d = MOD(tmp, nGH1d)      ! 0,1,2
!         tmp   = tmp / nGH1d          ! divisione intera

!         ! indice Fortran 1-based
!         nodes0(l*dimZ + d) = x1(idx1d+1)
!         w_tot = w_tot * w1(idx1d+1)
!      END DO

!      weights0(l+1) = w_tot
!   END DO

!   !-------------------------
!   ! Design X e Z: X = (t-1), Z = (1, t-1)
!   !-------------------------
!   DO i = 1, nsub0
!      DO j = 1, n0
!         DO k = 1, Q0
!            X0(i,j,k,1) = DBLE(j-1)

!            Z0(i,j,k,1) = 1.d0
!            Z0(i,j,k,2) = DBLE(j-1)
!         END DO
!      END DO
!   END DO

!   !-------------------------
!   ! Parametri b0 nel formato di loglik
!   ! ident_type0 = 1 (unit-variance):
!   !   - tutti i lambda(k,q) liberi (per item che caricano su q)
!   !   - sigma_eta2(q) fissate a 1 dentro loglik (nessun parametro in b0)
!   !
!   ! NPM = ntre
!   !       + n_lambda (qui = nvar0 per simple structure)
!   !       + 3*nvar0  (gamma, log theta_k, log sigmak)
!   !       + Q0*p_fix0
!   !       + Q0*p_ran0*(p_ran0+1)/2
!   !
!   !     = ntre0 + 4*nvar0 + Q0*p_fix0 + Q0*p_ran0*(p_ran0+1)/2
!   !-------------------------
!   npm0 = ntre0 + 4*nvar0 + Q0*p_fix0 + Q0*p_ran0*(p_ran0+1)/2
!   ALLOCATE(b0(npm0))

!   idx = 0

!   ! (1) soglie libere: per ciascun item 1 soglia libera (la terza)
!   idx = idx + 1 ; b0(idx) = 1.1d0    ! item 1
!   idx = idx + 1 ; b0(idx) = 1.2d0    ! item 2
!   idx = idx + 1 ; b0(idx) = 1.3d0    ! item 3
!   idx = idx + 1 ; b0(idx) = 1.4d0    ! item 4

!   ! (2) loadings λ(k,q) liberi per tutti gli item (simple structure)
!   ! ordine in loglik (ident_type0=1):
!   !   q=1: k=1,2 con fac_item0(k)=1
!   !   q=2: k=3,4 con fac_item0(k)=2
!   idx = idx + 1 ; b0(idx) = 1.3d0    ! lambda(1,1)
!   idx = idx + 1 ; b0(idx) = 0.7d0    ! lambda(2,1)
!   idx = idx + 1 ; b0(idx) = 1.4d0    ! lambda(3,2)
!   idx = idx + 1 ; b0(idx) = 0.7d0    ! lambda(4,2)

!   ! (3) gamma(k) intercette item
!   idx = idx + 1 ; b0(idx) = 0.10d0
!   idx = idx + 1 ; b0(idx) = 0.20d0
!   idx = idx + 1 ; b0(idx) = 0.30d0
!   idx = idx + 1 ; b0(idx) = 0.40d0

!   ! (4) log theta_k(k) varianza errore di misura
!   idx = idx + 1 ; b0(idx) = LOG(0.5d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.6d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.7d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.8d0)

!   ! (5) beta(q,1) (p_fix0 = 1)
!   idx = idx + 1 ; b0(idx) = 0.20d0       ! beta(1,1)
!   idx = idx + 1 ; b0(idx) = 0.10d0       ! beta(2,1)

!   ! (6) L_chol per ogni fattore q: diag log, poi sotto-diagonale (p_ran0 = 2)
!   ! fattore 1
!   idx = idx + 1 ; b0(idx) = LOG(0.8d0)   ! L11(1)
!   idx = idx + 1 ; b0(idx) = LOG(0.9d0)   ! L22(1)
!   idx = idx + 1 ; b0(idx) = 0.10d0       ! L21(1)

!   ! fattore 2
!   idx = idx + 1 ; b0(idx) = LOG(0.7d0)   ! L11(2)
!   idx = idx + 1 ; b0(idx) = LOG(0.6d0)   ! L22(2)
!   idx = idx + 1 ; b0(idx) = -0.05d0      ! L21(2)

!   ! (7) log sigmak(k) (varianze random sugli item)
!   idx = idx + 1 ; b0(idx) = LOG(0.3d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.4d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.5d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.6d0)

!   IF (idx /= npm0) THEN
!      PRINT *, 'ERRORE: idx != npm0', idx, npm0
!   END IF

!   !-------------------------
!   ! Chiamata singola a loglik
!   !-------------------------
!   CALL loglik(b0        = b0,    &
!               Y0        = Y0,    &
!               nvar0     = nvar0, &
!               nsub0     = nsub0, &
!               nobs0     = nobs0, &
!               Q0        = Q0,    &
!               n0        = n0,    &
!               nmes0     = nmes0, &
!               c0        = c0,    &
!               cm10      = cm10,  &
!               ntre0     = ntre0, &
!               nquad0    = nquad0,&
!               nitem0    = nitem0,&
!               nodes0    = nodes0,&
!               weights0  = weights0,&
!               npm0      = npm0,  &
!               p_fix0    = p_fix0,&
!               p_ran0    = p_ran0,&
!               X0        = X0,    &
!               Z0        = Z0,    &
!               fac_item0 = fac_item0, &
!               ident_type0 = ident_type0, &
!               re_item0 = re_item0, &
!               vrais     = vrais)

!   PRINT *, '*** -loglikelihood (unit-variance, GH 3 punti) =', vrais

!   DEALLOCATE(b0, Y0, nitem0, nodes0, weights0)

! END PROGRAM debug_one_iter

! PROGRAM debug_one_iter_agh
!   USE core
!   IMPLICIT NONE

!   INTEGER, PARAMETER :: nvar0       = 4      ! 4 item
!   INTEGER, PARAMETER :: nsub0       = 1      ! 1 soggetto
!   INTEGER, PARAMETER :: Q0          = 2      ! 2 fattori
!   INTEGER, PARAMETER :: n0          = 4      ! 4 occasioni (t=1..4)
!   INTEGER, PARAMETER :: p_fix0      = 1      ! solo slope (t-1)
!   INTEGER, PARAMETER :: p_ran0      = 2      ! random intercept + random slope
!   INTEGER, PARAMETER :: ident_type0 = 1      ! 0 = anchor λ=1; 1 = unit-variance
!   INTEGER, PARAMETER :: re_item0    = 1      ! 1 = random sugli item

!   ! Gauss–Hermite 1D per normale "standardizzata" (qui 3 punti)
!   INTEGER, PARAMETER :: nGH1d       = 3

!   INTEGER :: ntre0, npm0, nobs0
!   INTEGER :: i, k, j, dimZ, idx
!   INTEGER :: nquad0
!   INTEGER :: l, d, tmp, idx1d

!   DOUBLE PRECISION :: vrais, vrais_agh, w_tot

!   INTEGER :: nmes0(nsub0,nvar0)
!   INTEGER :: c0(nvar0), cm10(nvar0), fac_item0(nvar0)
!   INTEGER, ALLOCATABLE :: Y0(:), nitem0(:)
!   DOUBLE PRECISION, ALLOCATABLE :: nodes0(:), weights0(:), b0(:)
!   DOUBLE PRECISION :: X0(nsub0,n0,Q0,p_fix0)
!   DOUBLE PRECISION :: Z0(nsub0,n0,Q0,p_ran0)

!   DOUBLE PRECISION, DIMENSION(nGH1d) :: x1, w1

!   !-------------------------
!   ! Struttura categorie
!   !-------------------------
!   c0   = 4                  ! 4 categorie per ogni item
!   cm10 = c0 - 1             ! C-1 soglie
!   ntre0 = SUM(cm10 - 2)     ! ntre0 = 1 soglia libera per item -> 4

!   ! item 1,2 -> fattore 1 ; item 3,4 -> fattore 2
!   fac_item0 = (/1,1,2,2/)

!   !-------------------------
!   ! nmes, Y, nitem (ordine i,k,j)
!   !-------------------------
!   nmes0(1,:) = (/2,2,2,4/)   ! numero misure per item (soggetto 1)

!   nobs0 = SUM(nmes0)
!   ALLOCATE(Y0(nobs0), nitem0(nobs0))

!   idx = 0
!   DO i = 1, nsub0
!      DO k = 1, nvar0
!         DO j = 1, nmes0(i,k)
!            idx = idx + 1
!            nitem0(idx) = j          ! tempo = j

!            SELECT CASE (k)
!            CASE (1)
!               Y0(idx) = j           ! 1,2
!            CASE (2)
!               Y0(idx) = 3 - j       ! 2,1
!            CASE (3)
!               Y0(idx) = j + 1       ! 2,3
!            CASE (4)
!               Y0(idx) = j           ! 1,2
!            CASE DEFAULT
!               WRITE(*,*) 'ERRORE: k fuori range in costruzione Y0, k=', k
!               STOP
!            END SELECT
!         END DO
!      END DO
!   END DO

!   !-------------------------
!   ! Dimensione totale random
!   ! dimZ = Q*p_ran + nvar (perché re_item0=1)
!   !-------------------------
!   dimZ = Q0*p_ran0 + nvar0

!   !-------------------------
!   ! Gauss–Hermite 1D (3 punti)
!   ! NB: qui li uso solo come test; per AGH serio
!   !     in R userai fastGHQuad con la scalatura corretta.
!   !-------------------------
!   x1(1) = -1.7320508075688772d0
!   x1(2) =  0.d0
!   x1(3) =  1.7320508075688772d0

!   w1(1) = 0.1666666666666667d0
!   w1(2) = 0.6666666666666667d0
!   w1(3) = 0.1666666666666667d0

!   !-------------------------
!   ! Numero totale di nodi GH multidimensionale:
!   ! nquad0 = nGH1d^dimZ (qui = 3^8 = 6561)
!   !-------------------------
!   nquad0 = 1
!   DO d = 1, dimZ
!      nquad0 = nquad0 * nGH1d
!   END DO

!   ALLOCATE(nodes0(nquad0*dimZ), weights0(nquad0))

!   !-------------------------
!   ! Costruzione griglia tensoriale:
!   ! per ogni l = 0..nquad0-1, rappresento l in base 3
!   ! e prendo il prodotto dei pesi 1D
!   !-------------------------
!   DO l = 0, nquad0-1
!      tmp   = l
!      w_tot = 1.d0

!      DO d = 1, dimZ
!         idx1d = MOD(tmp, nGH1d)      ! 0,1,2
!         tmp   = tmp / nGH1d          ! divisione intera

!         ! indice Fortran 1-based
!         nodes0(l*dimZ + d) = x1(idx1d+1)
!         w_tot = w_tot * w1(idx1d+1)
!      END DO

!      weights0(l+1) = w_tot
!   END DO

!   !-------------------------
!   ! Design X e Z: X = (t-1), Z = (1, t-1)
!   !-------------------------
!   DO i = 1, nsub0
!      DO j = 1, n0
!         DO k = 1, Q0
!            X0(i,j,k,1) = DBLE(j-1)

!            Z0(i,j,k,1) = 1.d0
!            Z0(i,j,k,2) = DBLE(j-1)
!         END DO
!      END DO
!   END DO

!   !-------------------------
!   ! Parametri b0 nel formato di loglik / loglik_agh
!   ! ident_type0 = 1 (unit-variance):
!   !   - tutti i lambda(k,q) liberi (per item che caricano su q)
!   !   - sigma_eta2(q) fissate a 1 dentro loglik
!   !
!   ! NPM = ntre
!   !       + n_lambda (qui = nvar0)
!   !       + 3*nvar0  (gamma, log theta_k, log sigmak)
!   !       + Q0*p_fix0
!   !       + Q0*p_ran0*(p_ran0+1)/2
!   !-------------------------
!   npm0 = ntre0 + 4*nvar0 + Q0*p_fix0 + Q0*p_ran0*(p_ran0+1)/2
!   ALLOCATE(b0(npm0))

!   idx = 0

!   ! (1) soglie libere: per ciascun item 1 soglia libera (la terza)
!   idx = idx + 1 ; b0(idx) = 1.1d0    ! item 1
!   idx = idx + 1 ; b0(idx) = 1.2d0    ! item 2
!   idx = idx + 1 ; b0(idx) = 1.3d0    ! item 3
!   idx = idx + 1 ; b0(idx) = 1.4d0    ! item 4

!   ! (2) loadings λ(k,q) liberi (simple structure)
!   idx = idx + 1 ; b0(idx) = 1.3d0    ! lambda(1,1)
!   idx = idx + 1 ; b0(idx) = 0.7d0    ! lambda(2,1)
!   idx = idx + 1 ; b0(idx) = 1.4d0    ! lambda(3,2)
!   idx = idx + 1 ; b0(idx) = 0.7d0    ! lambda(4,2)

!   ! (3) gamma(k) intercette item
!   idx = idx + 1 ; b0(idx) = 0.10d0
!   idx = idx + 1 ; b0(idx) = 0.20d0
!   idx = idx + 1 ; b0(idx) = 0.30d0
!   idx = idx + 1 ; b0(idx) = 0.40d0

!   ! (4) log theta_k(k) varianza errore di misura
!   idx = idx + 1 ; b0(idx) = LOG(0.5d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.6d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.7d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.8d0)

!   ! (5) beta(q,1) (p_fix0 = 1)
!   idx = idx + 1 ; b0(idx) = 0.20d0       ! beta(1,1)
!   idx = idx + 1 ; b0(idx) = 0.10d0       ! beta(2,1)

!   ! (6) L_chol per ogni fattore q: diag log, poi sotto-diagonale (p_ran0 = 2)
!   ! fattore 1
!   idx = idx + 1 ; b0(idx) = LOG(0.8d0)   ! L11(1)
!   idx = idx + 1 ; b0(idx) = LOG(0.9d0)   ! L22(1)
!   idx = idx + 1 ; b0(idx) = 0.10d0       ! L21(1)

!   ! fattore 2
!   idx = idx + 1 ; b0(idx) = LOG(0.7d0)   ! L11(2)
!   idx = idx + 1 ; b0(idx) = LOG(0.6d0)   ! L22(2)
!   idx = idx + 1 ; b0(idx) = -0.05d0      ! L21(2)

!   ! (7) log sigmak(k) (varianze random sugli item)
!   idx = idx + 1 ; b0(idx) = LOG(0.3d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.4d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.5d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.6d0)

!   IF (idx /= npm0) THEN
!      PRINT *, 'ERRORE: idx != npm0', idx, npm0
!   END IF

!   !-------------------------
!   ! Chiamata a loglik (GH "globale")
!   !-------------------------
!   CALL loglik(b0        = b0,    &
!               Y0        = Y0,    &
!               nvar0     = nvar0, &
!               nsub0     = nsub0, &
!               nobs0     = nobs0, &
!               Q0        = Q0,    &
!               n0        = n0,    &
!               nmes0     = nmes0, &
!               c0        = c0,    &
!               cm10      = cm10,  &
!               ntre0     = ntre0, &
!               nquad0    = nquad0,&
!               nitem0    = nitem0,&
!               nodes0    = nodes0,&
!               weights0  = weights0,&
!               npm0      = npm0,  &
!               p_fix0    = p_fix0,&
!               p_ran0    = p_ran0,&
!               X0        = X0,    &
!               Z0        = Z0,    &
!               fac_item0 = fac_item0, &
!               ident_type0 = ident_type0, &
!               re_item0  = re_item0, &
!               vrais     = vrais)

!   PRINT *, '*** -loglik (GH globale) =', vrais

!   !-------------------------
!   ! Chiamata a loglik_agh (stessi nodi/pesi,
!   ! ma interpretati come soggetto-specifici)
!   ! Con nsub0=1 i due valori DEVONO coincidere.
!   !-------------------------
!   CALL loglik_agh(b0        = b0,    &
!                   Y0        = Y0,    &
!                   nvar0     = nvar0, &
!                   nsub0     = nsub0, &
!                   nobs0     = nobs0, &
!                   Q0        = Q0,    &
!                   n0        = n0,    &
!                   nmes0     = nmes0, &
!                   c0        = c0,    &
!                   cm10      = cm10,  &
!                   ntre0     = ntre0, &
!                   nquad0    = nquad0,&
!                   nitem0    = nitem0,&
!                   nodes0    = nodes0,&   ! lunghezza = nsub0*nquad0*dimZ
!                   weights0  = weights0,& ! lunghezza = nsub0*nquad0
!                   npm0      = npm0,  &
!                   p_fix0    = p_fix0,&
!                   p_ran0    = p_ran0,&
!                   X0        = X0,    &
!                   Z0        = Z0,    &
!                   fac_item0 = fac_item0, &
!                   ident_type0 = ident_type0, &
!                   re_item0  = re_item0, &
!                   vrais     = vrais_agh)

!   PRINT *, '*** -loglik (loglik_agh, stessi nodi) =', vrais_agh

!   DEALLOCATE(b0, Y0, nitem0, nodes0, weights0)

! END PROGRAM debug_one_iter_agh




! PROGRAM debug_long_one_iter
!   USE core
!   IMPLICIT NONE

!   INTEGER, PARAMETER :: nvar0   = 2      ! 2 item
!   INTEGER, PARAMETER :: nsub0   = 1      ! 1 soggetto
!   INTEGER, PARAMETER :: n0      = 3      ! 3 occasioni (t=1..3)
!   INTEGER, PARAMETER :: nquad0  = 1      ! 1 nodo di integrazione
!   INTEGER, PARAMETER :: re_item0 = 1     ! random sugli item
!   INTEGER, PARAMETER :: theta_time0 = 0     ! invarinza theta

!   INTEGER :: ntre0, npm0, nobs0, dimZ
!   INTEGER :: nmes0(nsub0,nvar0)
!   INTEGER :: c0(nvar0), cm10(nvar0)
!   INTEGER, ALLOCATABLE :: Y0(:), nitem0(:)
!   DOUBLE PRECISION, ALLOCATABLE :: nodes0(:), weights0(:), b0(:)
!   DOUBLE PRECISION :: vrais
!   INTEGER :: i, k, j, idx

!   !-------------------------
!   ! Struttura categorie
!   !-------------------------
!   c0   = 4                  ! 4 categorie per ogni item
!   cm10 = c0 - 1             ! C-1 soglie = 3
!   ! cm3 = cm1 - 2 = 1 soglia libera per item
!   ntre0 = SUM(cm10 - 2)     ! ntre0 = 2

!   !-------------------------
!   ! nmes, Y, nitem (ordine i,k,j)
!   !-------------------------
!   nmes0(1,1) = 3    ! item 1 misurato 3 volte
!   nmes0(1,2) = 3    ! item 2 misurato 3 volte

!   nobs0 = SUM(nmes0)
!   ALLOCATE(Y0(nobs0), nitem0(nobs0))

!   idx = 0
!   DO k = 1, nvar0
!      DO j = 1, nmes0(1,k)
!         idx = idx + 1
!         nitem0(idx) = j          ! tempo = j (1..3)

!         SELECT CASE (k)
!         CASE (1)
!            ! item 1: risposte crescenti 1,2,3
!            Y0(idx) = j
!         CASE (2)
!            ! item 2: risposte 3,2,1
!            Y0(idx) = 4 - j
!         END SELECT
!      END DO
!   END DO

!   !-------------------------
!   ! Nodi (eta(1..n) + a_ik)
!   ! dimZ = n + nvar se re_item0=1
!   !-------------------------
!   dimZ = n0 + nvar0          ! 3 + 2 = 5
!   ALLOCATE(nodes0(nquad0*dimZ), weights0(nquad0))

!   ! un solo nodo, valori fissi
!   nodes0 = (/ 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0 /)
!   weights0(1) = 1.d0

!   !-------------------------
!   ! Parametri b0 nel formato di loglik_long:
!   !
!   !  1) soglie libere (ntre0 = 2)
!   !  2) lambda_k        [nvar0 = 2]
!   !  3) gamma_k         [nvar0 = 2]
!   !  4) log theta_k     [nvar0 = 2]
!   !  5) log sigmak_k    [nvar0 = 2] (perché re_item0=1)
!   !  6) kappa_2..n0     [n0-1 = 2]  (kappa(1)=0 fissata in Fortran)
!   !  7) L_phi (n0 x n0):
!   !        - log diag L(t,t), t=1..n0 (3)
!   !        - sotto-diagonale: n0*(n0-1)/2 = 3
!   !
!   ! => npm0 = 2 + 2 + 2 + 2 + 2 + 2 + 6 = 18
!   !-------------------------
!   npm0 = ntre0 + nvar0 + nvar0 + nvar0 + nvar0 + (n0-1) + n0*(n0+1)/2
!   ALLOCATE(b0(npm0))

!   idx = 0

!   ! (1) soglie libere (terza soglia di ciascun item, ad esempio)
!   idx = idx + 1 ; b0(idx) = 1.1d0   ! item 1
!   idx = idx + 1 ; b0(idx) = 1.2d0   ! item 2

!   ! (2) lambda_k
!   idx = idx + 1 ; b0(idx) = 0.8d0   ! lambda(1)
!   idx = idx + 1 ; b0(idx) = 1.0d0   ! lambda(2)

!   ! (3) gamma_k
!   idx = idx + 1 ; b0(idx) = 0.0d0   ! gamma(1)
!   idx = idx + 1 ; b0(idx) = 0.2d0   ! gamma(2)

!   ! (4) log theta_k
!   idx = idx + 1 ; b0(idx) = LOG(0.5d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.7d0)

!   ! (5) log sigmak_k (random sugli item)
!   idx = idx + 1 ; b0(idx) = LOG(0.3d0)
!   idx = idx + 1 ; b0(idx) = LOG(0.4d0)

!   ! (6) kappa_2..n0 (kappa_1 = 0)
!   idx = idx + 1 ; b0(idx) = 0.1d0   ! kappa_2
!   idx = idx + 1 ; b0(idx) = 0.2d0   ! kappa_3

!   ! (7) L_phi: diag log, poi sotto-diagonale (n0=3)
!   ! diagonali
!   idx = idx + 1 ; b0(idx) = LOG(1.0d0)   ! L11
!   idx = idx + 1 ; b0(idx) = LOG(0.9d0)   ! L22
!   idx = idx + 1 ; b0(idx) = LOG(0.8d0)   ! L33

!   ! sotto-diagonale: col=1 -> row=2,3 ; col=2 -> row=3
!   idx = idx + 1 ; b0(idx) = 0.10d0       ! L21
!   idx = idx + 1 ; b0(idx) = 0.00d0       ! L31
!   idx = idx + 1 ; b0(idx) = 0.05d0       ! L32

!   IF (idx /= npm0) THEN
!      PRINT *, 'ERRORE: idx != npm0', idx, npm0
!   END IF

!   !-------------------------
!   ! Chiamata singola a loglik_long
!   !-------------------------
!   CALL loglik_long(b0        = b0,      &
!                    Y0        = Y0,      &
!                    nvar0     = nvar0,   &
!                    nsub0     = nsub0,   &
!                    nobs0     = nobs0,   &
!                    n0        = n0,      &
!                    nmes0     = nmes0,   &
!                    c0        = c0,      &
!                    cm10      = cm10,    &
!                    ntre0     = ntre0,   &
!                    nquad0    = nquad0,  &
!                    nitem0    = nitem0,  &
!                    nodes0    = nodes0,  &
!                    weights0  = weights0,&
!                    npm0      = npm0,    &
!                    re_item0  = re_item0,&
!                    theta_time0 = theta_time0, &
!                    vrais     = vrais)

!   PRINT *, '*** -loglikelihood (CFA longitudinale) =', vrais

!   DEALLOCATE(b0, Y0, nitem0, nodes0, weights0)

! END PROGRAM debug_long_one_iter


!=======================================================================
! ADD-ONLY EXTENSION: growth AGH two-tier equivalent to full joint AGH
! The original source above is preserved byte-for-byte and unchanged.
!=======================================================================
MODULE growth_agh_two_tier_equivalent_extension
  USE shared
  IMPLICIT NONE
CONTAINS

SUBROUTINE loglik_growth_tt_agh_full_equiv(b0, Y0, nvar0, nsub0, nobs0, n0,        &
                                nmes0, c0, cm10, ntre0,               &
                                nquad_eta0, nodes_eta0, weights_eta0, &
                                nquad_re0, nodes_re0, weights_re0,    &
                                nitem0, npm0, re_item0, theta_time0,  &
                                delta_free0, monotone_tau0,           &
                                update_gh0, fit_id0,                  &
                                Q0, loadmask0,                        &
                            p_fix0, p_ran0, X0, Z0, delta_mask0,    &
                            delta_x_col0, delta_z_col0,              &
                                vrais)                                &
        BIND(C, name="loglik_growth_tt_agh_full_equiv")

  USE shared
  USE core, ONLY: cleanup_growth_tt_allocations, setup_load_indices_from_loadmask_growth, &
                  decode_lambda_sum_constraints_growth, refresh_lambda_act_growth,       &
                  decode_gamma_sum_constraints_growth, decode_growth_thresholds,          &
                  decode_growth_beta, decode_growth_delta, decode_growth_random_structure, &
                  growth_n_lambda_par, growth_n_gamma_par, growth_n_beta_par,              &
                  growth_n_zeta_par, build_growth_eta_moments_design,                      &
                  build_agh_grid_full_equiv_tt, vrais_long_tt_agh_full_equiv_fixed
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0
  INTEGER, INTENT(IN) :: ntre0, npm0
  INTEGER, INTENT(IN) :: nquad_eta0, nquad_re0
  INTEGER, INTENT(IN) :: re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0
  INTEGER, INTENT(IN) :: update_gh0, fit_id0
  INTEGER, INTENT(IN) :: Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0
  INTEGER, INTENT(IN) :: delta_x_col0, delta_z_col0

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0),       INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0),    INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0

  DOUBLE PRECISION, DIMENSION(npm0),       INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*),          INTENT(IN) :: nodes_eta0
  DOUBLE PRECISION, DIMENSION(nquad_eta0), INTENT(IN) :: weights_eta0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: weights_re0

  DOUBLE PRECISION, INTENT(OUT) :: vrais

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.0D9

  INTEGER :: k, t, r, row, col, pos, jj, pos_tau, info
  INTEGER :: cache_key
  INTEGER :: qidx, dimEta, dimU, ntheta_est, expected_npm, ntre_calc, nX, nZ
  INTEGER :: n_lambda_free, n_gamma_free, n_Lu_est
  INTEGER :: nf, h
  DOUBLE PRECISION :: tmp_tau, ll_agh
  LOGICAL :: ok_grid

  INTEGER, ALLOCATABLE :: cm3(:)
  DOUBLE PRECISION, ALLOCATABLE :: delta_g(:)
  DOUBLE PRECISION, ALLOCATABLE :: Lu(:,:), Sigma_u(:,:), zeta_var(:,:)

  vrais = BAD_OBJ

  IF (n0 < 2 .OR. Q0 < 1 .OR. ABS(p_fix0) < 1 .OR. p_ran0 < 1) RETURN

  nvar = nvar0
  nsub = nsub0
  nobs = nobs0
  n    = n0
  Q    = Q0
  npm  = npm0

  dimEta = n * Q
  dimU   = Q * p_ran0
  dim_eta = dimEta

  nquad     = nquad_eta0
  nquad_eta = nquad_eta0
  nquad_re  = nquad_re0

  re_item_flag    = re_item0
  theta_time_flag = theta_time0
  phi_struct_flag = 0
  growth_subject_specific_design = MERGE(1, 0, p_fix0 > 0)
  growth_design_nsub = MERGE(nsub0, 1, growth_subject_specific_design == 1)
  p_fix = ABS(p_fix0)
  p_ran = p_ran0
  growth_delta_x_col = delta_x_col0
  growth_delta_z_col = delta_z_col0


  cache_key = IEOR(fit_id0, INT(Z'40000000'))

  IF (cache_key /= agh_tt_fe_fit_id_cached) THEN
     agh_tt_fe_grid_ready = .FALSE.
     agh_tt_fe_fit_id_cached = cache_key

     IF (ALLOCATED(agh_tt_fe_zeta_eta)) DEALLOCATE(agh_tt_fe_zeta_eta)
     IF (ALLOCATED(agh_tt_fe_alpha)) DEALLOCATE(agh_tt_fe_alpha)
     IF (ALLOCATED(agh_tt_fe_sdev)) DEALLOCATE(agh_tt_fe_sdev)
     IF (ALLOCATED(agh_tt_fe_logconst)) DEALLOCATE(agh_tt_fe_logconst)
     IF (ALLOCATED(agh_tt_fe_modes)) DEALLOCATE(agh_tt_fe_modes)

     agh_tt_fe_dimEta_cached = -1
     agh_tt_fe_nvar_cached = -1
     agh_tt_fe_nquad_eta_cached = -1
     agh_tt_fe_nquad_re_cached = -1
     agh_tt_fe_nsub_cached = -1
  END IF

  IF (ALLOCATED(lambda_long)) DEALLOCATE(lambda_long)
  IF (ALLOCATED(gamma_long))  DEALLOCATE(gamma_long)
  IF (ALLOCATED(tau_long))    DEALLOCATE(tau_long)

  CALL cleanup_growth_tt_allocations()

  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs))
  ALLOCATE(c(nvar), cm1(nvar))
  ALLOCATE(nodes(nquad*dimEta), weights(nquad))
  ALLOCATE(nodes_re(nquad_re), weights_re(nquad_re), logweights_re(nquad_re))
  ALLOCATE(loadmask(nvar,Q))
  ALLOCATE(lambda(nvar,Q), gamma(nvar))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(theta_long(nvar,n), invsd_long(nvar,n))
  ALLOCATE(sigmak(nvar), sqrt_sigmak(nvar))
  ALLOCATE(kappa(dimEta), L_phi(dimEta,dimEta))
  ALLOCATE(kk_base(nvar))
  ALLOCATE(tau(SUM(cm10)))
  ALLOCATE(cm3(nvar))
  ALLOCATE(delta_g(Q))
  ALLOCATE(beta(Q,p_fix), X_eta(growth_design_nsub,n,Q,p_fix), &
           Z_eta(growth_design_nsub,n,Q,p_ran))
  ALLOCATE(growth_delta_mask(growth_design_nsub,n,Q))
  ALLOCATE(Lu(dimU,dimU), Sigma_u(dimU,dimU), zeta_var(n,Q))

  nmes     = nmes0
  Y        = Y0
  nitem    = nitem0
  c        = c0
  cm1      = cm10
  loadmask = loadmask0

  nX = growth_design_nsub * n * Q * p_fix
  nZ = growth_design_nsub * n * Q * p_ran
  X_eta = RESHAPE(X0(1:nX), (/ growth_design_nsub, n, Q, p_fix /))
  Z_eta = RESHAPE(Z0(1:nZ), (/ growth_design_nsub, n, Q, p_ran /))
  growth_delta_mask = RESHAPE(delta_mask0(1:growth_design_nsub*n*Q), &
                              (/ growth_design_nsub, n, Q /))

  nodes(1:nquad*dimEta) = nodes_eta0(1:nquad*dimEta)
  weights(1:nquad) = weights_eta0(1:nquad)

  nodes_re(1:nquad_re) = nodes_re0(1:nquad_re)
  weights_re(1:nquad_re) = weights_re0(1:nquad_re)

  DO r = 1, nquad_re
     IF (weights_re(r) > 0.d0) THEN
        logweights_re(r) = LOG(weights_re(r))
     ELSE
        logweights_re(r) = -HUGE(1.d0)
     END IF
  END DO

  CALL setup_load_indices_from_loadmask_growth(info)
  IF (info /= 0) GOTO 999

  DO k = 1, nvar
     cm3(k) = cm1(k) - 2
     IF (cm3(k) < 0) cm3(k) = 0
  END DO

  ntre_calc = SUM(cm3)
  IF (ntre_calc /= ntre0) GOTO 999
  CALL decode_growth_thresholds(b0, npm0, ntre0, monotone_tau0, info)
  IF (info /= 0) GOTO 999

  pos = ntre0

  CALL decode_lambda_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999

  CALL refresh_lambda_act_growth()

  CALL decode_gamma_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999

  theta_k(:) = 1.d0
  invsd_k(:) = 1.d0
  theta_long(:,:) = 1.d0
  invsd_long(:,:) = 1.d0

  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) THEN
           pos = pos + 1
           theta_k(k) = EXP(b0(pos))
           invsd_k(k) = 1.d0 / SQRT(theta_k(k))

           DO t = 1, n
              theta_long(k,t) = theta_k(k)
              invsd_long(k,t) = invsd_k(k)
           END DO
        END IF
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) THEN
           DO t = 1, n
              pos = pos + 1
              theta_long(k,t) = EXP(b0(pos))
              invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
           END DO

           theta_k(k) = theta_long(k,1)
           invsd_k(k) = invsd_long(k,1)
        END IF
     END DO
  END IF
  CALL decode_growth_beta(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999
   CALL decode_growth_delta(b0, npm0, pos, delta_free0, delta_g, info)
   IF (info /= 0) GOTO 999
  CALL decode_growth_random_structure(b0, npm0, pos, dimU, Lu, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 999

  IF (re_item_flag == 1) THEN
     DO k = 1, nvar
        pos = pos + 1
        sigmak(k) = EXP(b0(pos))
        sqrt_sigmak(k) = SQRT(sigmak(k))
     END DO
  ELSE
     sigmak(:) = 0.d0
     sqrt_sigmak(:) = 0.d0
  END IF

  ntheta_est = 0

  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) ntheta_est = ntheta_est + 1
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) ntheta_est = ntheta_est + n
     END DO
  END IF
  n_lambda_free = growth_n_lambda_par()
  n_gamma_free = growth_n_gamma_par()
  n_Lu_est = dimU * (dimU + 1) / 2

  expected_npm = ntre0                                      &
               + n_lambda_free                              &
               + n_gamma_free                               &
               + ntheta_est                                 &
               + growth_n_beta_par()                        &
               + MERGE(Q, 0, delta_free0 == 1)              &
               + n_Lu_est                                   &
               + growth_n_zeta_par()                        &
               + MERGE(nvar, 0, re_item_flag == 1)

  IF (npm0 /= expected_npm) GOTO 999
  IF (pos /= npm0) GOTO 999

  CALL build_growth_eta_moments_design(delta_g, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 999

  IF (update_gh0 == 1 .OR. .NOT. agh_tt_fe_grid_ready) THEN
     CALL build_agh_grid_full_equiv_tt(ok_grid)

     IF (.NOT. ok_grid .OR. .NOT. agh_tt_fe_grid_ready) THEN
        vrais = BAD_OBJ
        GOTO 999
     END IF
  END IF

  CALL vrais_long_tt_agh_full_equiv_fixed(ll_agh)

  IF (ll_agh /= ll_agh .OR. ABS(ll_agh) > 1.d300) THEN
     vrais = BAD_OBJ
  ELSE
     vrais = -ll_agh
  END IF

999 CONTINUE

  IF (ALLOCATED(lambda_long)) DEALLOCATE(lambda_long)
  IF (ALLOCATED(gamma_long))  DEALLOCATE(gamma_long)
  IF (ALLOCATED(tau_long))    DEALLOCATE(tau_long)

  CALL cleanup_growth_tt_allocations()

  IF (ALLOCATED(cm3)) DEALLOCATE(cm3)
  IF (ALLOCATED(delta_g)) DEALLOCATE(delta_g)
  IF (ALLOCATED(Lu)) DEALLOCATE(Lu)
  IF (ALLOCATED(Sigma_u)) DEALLOCATE(Sigma_u)
  IF (ALLOCATED(zeta_var)) DEALLOCATE(zeta_var)

END SUBROUTINE loglik_growth_tt_agh_full_equiv


SUBROUTINE loglik_growth_tt_agh_full_equiv_i_vec(b0, Y0, nvar0, nsub0, nobs0, n0,        &
                                nmes0, c0, cm10, ntre0,               &
                                nquad_eta0, nodes_eta0, weights_eta0, &
                                nquad_re0, nodes_re0, weights_re0,    &
                                nitem0, npm0, re_item0, theta_time0,  &
                                delta_free0, monotone_tau0,           &
                                update_gh0, fit_id0,                  &
                                Q0, loadmask0,                        &
                            p_fix0, p_ran0, X0, Z0, delta_mask0,    &
                            delta_x_col0, delta_z_col0,              &
                                loglik_i)                             &
        BIND(C, name="loglik_growth_tt_agh_full_equiv_i_vec")

  USE shared
  USE core, ONLY: cleanup_growth_tt_allocations, setup_load_indices_from_loadmask_growth, &
                  decode_lambda_sum_constraints_growth, refresh_lambda_act_growth,       &
                  decode_gamma_sum_constraints_growth, decode_growth_thresholds,          &
                  decode_growth_beta, decode_growth_delta, decode_growth_random_structure, &
                  growth_n_lambda_par, growth_n_gamma_par, growth_n_beta_par,              &
                  growth_n_zeta_par, build_growth_eta_moments_design,                      &
                  build_agh_grid_full_equiv_tt,                                          &
                  vrais_long_i_tt_agh_full_equiv_fixed
  USE datindiv, ONLY: nmescur
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0
  INTEGER, INTENT(IN) :: ntre0, npm0
  INTEGER, INTENT(IN) :: nquad_eta0, nquad_re0
  INTEGER, INTENT(IN) :: re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0
  INTEGER, INTENT(IN) :: update_gh0, fit_id0
  INTEGER, INTENT(IN) :: Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0
  INTEGER, INTENT(IN) :: delta_x_col0, delta_z_col0

  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0),       INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0),       INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0),    INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0

  DOUBLE PRECISION, DIMENSION(npm0),       INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*),          INTENT(IN) :: nodes_eta0
  DOUBLE PRECISION, DIMENSION(nquad_eta0), INTENT(IN) :: weights_eta0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0),  INTENT(IN) :: weights_re0

  DOUBLE PRECISION, DIMENSION(nsub0), INTENT(OUT) :: loglik_i

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.0D9

  INTEGER :: k, t, r, row, col, pos, jj, pos_tau, info
  INTEGER :: cache_key
  INTEGER :: qidx, dimEta, dimU, ntheta_est, expected_npm, ntre_calc, nX, nZ
  INTEGER :: n_lambda_free, n_gamma_free, n_Lu_est
  INTEGER :: nf, h
  DOUBLE PRECISION :: tmp_tau, tmp_i
  INTEGER :: i, sm
  LOGICAL :: ok_i
  LOGICAL :: ok_grid

  INTEGER, ALLOCATABLE :: cm3(:)
  DOUBLE PRECISION, ALLOCATABLE :: delta_g(:)
  DOUBLE PRECISION, ALLOCATABLE :: Lu(:,:), Sigma_u(:,:), zeta_var(:,:)

  loglik_i(:) = -BAD_OBJ

  IF (n0 < 2 .OR. Q0 < 1 .OR. ABS(p_fix0) < 1 .OR. p_ran0 < 1) RETURN

  nvar = nvar0
  nsub = nsub0
  nobs = nobs0
  n    = n0
  Q    = Q0
  npm  = npm0

  dimEta = n * Q
  dimU   = Q * p_ran0
  dim_eta = dimEta

  nquad     = nquad_eta0
  nquad_eta = nquad_eta0
  nquad_re  = nquad_re0

  re_item_flag    = re_item0
  theta_time_flag = theta_time0
  phi_struct_flag = 0
  growth_subject_specific_design = MERGE(1, 0, p_fix0 > 0)
  growth_design_nsub = MERGE(nsub0, 1, growth_subject_specific_design == 1)
  p_fix = ABS(p_fix0)
  p_ran = p_ran0
  growth_delta_x_col = delta_x_col0
  growth_delta_z_col = delta_z_col0


  cache_key = IEOR(fit_id0, INT(Z'40000000'))

  IF (cache_key /= agh_tt_fe_fit_id_cached) THEN
     agh_tt_fe_grid_ready = .FALSE.
     agh_tt_fe_fit_id_cached = cache_key

     IF (ALLOCATED(agh_tt_fe_zeta_eta)) DEALLOCATE(agh_tt_fe_zeta_eta)
     IF (ALLOCATED(agh_tt_fe_alpha)) DEALLOCATE(agh_tt_fe_alpha)
     IF (ALLOCATED(agh_tt_fe_sdev)) DEALLOCATE(agh_tt_fe_sdev)
     IF (ALLOCATED(agh_tt_fe_logconst)) DEALLOCATE(agh_tt_fe_logconst)
     IF (ALLOCATED(agh_tt_fe_modes)) DEALLOCATE(agh_tt_fe_modes)

     agh_tt_fe_dimEta_cached = -1
     agh_tt_fe_nvar_cached = -1
     agh_tt_fe_nquad_eta_cached = -1
     agh_tt_fe_nquad_re_cached = -1
     agh_tt_fe_nsub_cached = -1
  END IF

  IF (ALLOCATED(lambda_long)) DEALLOCATE(lambda_long)
  IF (ALLOCATED(gamma_long))  DEALLOCATE(gamma_long)
  IF (ALLOCATED(tau_long))    DEALLOCATE(tau_long)

  CALL cleanup_growth_tt_allocations()

  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs))
  ALLOCATE(c(nvar), cm1(nvar))
  ALLOCATE(nodes(nquad*dimEta), weights(nquad))
  ALLOCATE(nodes_re(nquad_re), weights_re(nquad_re), logweights_re(nquad_re))
  ALLOCATE(loadmask(nvar,Q))
  ALLOCATE(lambda(nvar,Q), gamma(nvar))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(theta_long(nvar,n), invsd_long(nvar,n))
  ALLOCATE(sigmak(nvar), sqrt_sigmak(nvar))
  ALLOCATE(kappa(dimEta), L_phi(dimEta,dimEta))
  ALLOCATE(kk_base(nvar))
  ALLOCATE(tau(SUM(cm10)))
  ALLOCATE(cm3(nvar))
  ALLOCATE(delta_g(Q))
  ALLOCATE(beta(Q,p_fix), X_eta(growth_design_nsub,n,Q,p_fix), &
           Z_eta(growth_design_nsub,n,Q,p_ran))
  ALLOCATE(growth_delta_mask(growth_design_nsub,n,Q))
  ALLOCATE(Lu(dimU,dimU), Sigma_u(dimU,dimU), zeta_var(n,Q))

  nmes     = nmes0
  Y        = Y0
  nitem    = nitem0
  c        = c0
  cm1      = cm10
  loadmask = loadmask0

  nX = growth_design_nsub * n * Q * p_fix
  nZ = growth_design_nsub * n * Q * p_ran
  X_eta = RESHAPE(X0(1:nX), (/ growth_design_nsub, n, Q, p_fix /))
  Z_eta = RESHAPE(Z0(1:nZ), (/ growth_design_nsub, n, Q, p_ran /))
  growth_delta_mask = RESHAPE(delta_mask0(1:growth_design_nsub*n*Q), &
                              (/ growth_design_nsub, n, Q /))

  nodes(1:nquad*dimEta) = nodes_eta0(1:nquad*dimEta)
  weights(1:nquad) = weights_eta0(1:nquad)

  nodes_re(1:nquad_re) = nodes_re0(1:nquad_re)
  weights_re(1:nquad_re) = weights_re0(1:nquad_re)

  DO r = 1, nquad_re
     IF (weights_re(r) > 0.d0) THEN
        logweights_re(r) = LOG(weights_re(r))
     ELSE
        logweights_re(r) = -HUGE(1.d0)
     END IF
  END DO

  CALL setup_load_indices_from_loadmask_growth(info)
  IF (info /= 0) GOTO 999

  DO k = 1, nvar
     cm3(k) = cm1(k) - 2
     IF (cm3(k) < 0) cm3(k) = 0
  END DO

  ntre_calc = SUM(cm3)
  IF (ntre_calc /= ntre0) GOTO 999
  CALL decode_growth_thresholds(b0, npm0, ntre0, monotone_tau0, info)
  IF (info /= 0) GOTO 999

  pos = ntre0

  CALL decode_lambda_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999

  CALL refresh_lambda_act_growth()

  CALL decode_gamma_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999

  theta_k(:) = 1.d0
  invsd_k(:) = 1.d0
  theta_long(:,:) = 1.d0
  invsd_long(:,:) = 1.d0

  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) THEN
           pos = pos + 1
           theta_k(k) = EXP(b0(pos))
           invsd_k(k) = 1.d0 / SQRT(theta_k(k))

           DO t = 1, n
              theta_long(k,t) = theta_k(k)
              invsd_long(k,t) = invsd_k(k)
           END DO
        END IF
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) THEN
           DO t = 1, n
              pos = pos + 1
              theta_long(k,t) = EXP(b0(pos))
              invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
           END DO

           theta_k(k) = theta_long(k,1)
           invsd_k(k) = invsd_long(k,1)
        END IF
     END DO
  END IF
  CALL decode_growth_beta(b0, npm0, pos, info)
  IF (info /= 0) GOTO 999
   CALL decode_growth_delta(b0, npm0, pos, delta_free0, delta_g, info)
   IF (info /= 0) GOTO 999
  CALL decode_growth_random_structure(b0, npm0, pos, dimU, Lu, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 999

  IF (re_item_flag == 1) THEN
     DO k = 1, nvar
        pos = pos + 1
        sigmak(k) = EXP(b0(pos))
        sqrt_sigmak(k) = SQRT(sigmak(k))
     END DO
  ELSE
     sigmak(:) = 0.d0
     sqrt_sigmak(:) = 0.d0
  END IF

  ntheta_est = 0

  IF (theta_time_flag == 0) THEN
     DO k = 1, nvar
        IF (cm1(k) > 1) ntheta_est = ntheta_est + 1
     END DO
  ELSE
     DO k = 1, nvar
        IF (cm1(k) > 1) ntheta_est = ntheta_est + n
     END DO
  END IF
  n_lambda_free = growth_n_lambda_par()
  n_gamma_free = growth_n_gamma_par()
  n_Lu_est = dimU * (dimU + 1) / 2

  expected_npm = ntre0                                      &
               + n_lambda_free                              &
               + n_gamma_free                               &
               + ntheta_est                                 &
               + growth_n_beta_par()                        &
               + MERGE(Q, 0, delta_free0 == 1)              &
               + n_Lu_est                                   &
               + growth_n_zeta_par()                        &
               + MERGE(nvar, 0, re_item_flag == 1)

  IF (npm0 /= expected_npm) GOTO 999
  IF (pos /= npm0) GOTO 999

  CALL build_growth_eta_moments_design(delta_g, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 999

  IF (update_gh0 == 1 .OR. .NOT. agh_tt_fe_grid_ready) THEN
     CALL build_agh_grid_full_equiv_tt(ok_grid)

     IF (.NOT. ok_grid .OR. .NOT. agh_tt_fe_grid_ready) THEN
        GOTO 999
     END IF
  END IF

  nmescur = 0

  DO i = 1, nsub
     sm = SUM(nmes(i,:))

     CALL vrais_long_i_tt_agh_full_equiv_fixed(i, sm, tmp_i, ok_i)

     IF (ok_i .AND. tmp_i == tmp_i .AND. ABS(tmp_i) <= 1.d300) THEN
        loglik_i(i) = tmp_i
     ELSE
        loglik_i(i) = -BAD_OBJ
     END IF

     nmescur = nmescur + sm
  END DO

999 CONTINUE

  IF (ALLOCATED(lambda_long)) DEALLOCATE(lambda_long)
  IF (ALLOCATED(gamma_long))  DEALLOCATE(gamma_long)
  IF (ALLOCATED(tau_long))    DEALLOCATE(tau_long)

  CALL cleanup_growth_tt_allocations()

  IF (ALLOCATED(cm3)) DEALLOCATE(cm3)
  IF (ALLOCATED(delta_g)) DEALLOCATE(delta_g)
  IF (ALLOCATED(Lu)) DEALLOCATE(Lu)
  IF (ALLOCATED(Sigma_u)) DEALLOCATE(Sigma_u)
  IF (ALLOCATED(zeta_var)) DEALLOCATE(zeta_var)

END SUBROUTINE loglik_growth_tt_agh_full_equiv_i_vec


SUBROUTINE get_growth_agh_tt_fe_grid_info(ready0, dim_eta0, dim_z0, nvar0, &
                                          nquad_eta0, nquad_re0, nsub0, fit_id0) &
     BIND(C, name="get_growth_agh_tt_fe_grid_info")
  USE, INTRINSIC :: ISO_C_BINDING
  USE shared
  IMPLICIT NONE
  INTEGER(C_INT), INTENT(OUT) :: ready0, dim_eta0, dim_z0, nvar0
  INTEGER(C_INT), INTENT(OUT) :: nquad_eta0, nquad_re0, nsub0, fit_id0

  ready0 = MERGE(1_C_INT, 0_C_INT, agh_tt_fe_grid_ready)
  dim_eta0 = INT(agh_tt_fe_dimEta_cached, C_INT)
  nvar0 = INT(agh_tt_fe_nvar_cached, C_INT)
  IF (agh_tt_fe_dimEta_cached >= 0 .AND. agh_tt_fe_nvar_cached >= 0) THEN
     dim_z0 = INT(agh_tt_fe_dimEta_cached + agh_tt_fe_nvar_cached, C_INT)
  ELSE
     dim_z0 = -1_C_INT
  END IF
  nquad_eta0 = INT(agh_tt_fe_nquad_eta_cached, C_INT)
  nquad_re0 = INT(agh_tt_fe_nquad_re_cached, C_INT)
  nsub0 = INT(agh_tt_fe_nsub_cached, C_INT)
  fit_id0 = INT(IEOR(agh_tt_fe_fit_id_cached, INT(Z'40000000')), C_INT)
END SUBROUTINE get_growth_agh_tt_fe_grid_info




SUBROUTINE get_growth_agh_tt_fe_grid_copy(zeta_eta_out, alpha_out, sdev_out, &
                                           logconst_out, modes_out,          &
                                           ready0, dim_eta0, dim_z0, nvar0,  &
                                           nquad_eta0, nquad_re0, nsub0, fit_id0) &
     BIND(C, name="get_growth_agh_tt_fe_grid_copy")
  USE, INTRINSIC :: ISO_C_BINDING
  USE shared
  IMPLICIT NONE

  REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT) :: zeta_eta_out
  REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT) :: alpha_out
  REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT) :: sdev_out
  REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT) :: logconst_out
  REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT) :: modes_out

  INTEGER(C_INT), INTENT(OUT) :: ready0, dim_eta0, dim_z0, nvar0
  INTEGER(C_INT), INTENT(OUT) :: nquad_eta0, nquad_re0, nsub0, fit_id0

  INTEGER :: d, k, l, i, idx, dimZ

  ready0 = 0_C_INT
  dim_eta0 = 0_C_INT
  dim_z0 = 0_C_INT
  nvar0 = 0_C_INT
  nquad_eta0 = 0_C_INT
  nquad_re0 = 0_C_INT
  nsub0 = 0_C_INT
  fit_id0 = INT(IEOR(agh_tt_fe_fit_id_cached, INT(Z'40000000')), C_INT)

  IF (.NOT. agh_tt_fe_grid_ready) RETURN
  IF (.NOT. ALLOCATED(agh_tt_fe_zeta_eta)) RETURN
  IF (.NOT. ALLOCATED(agh_tt_fe_alpha)) RETURN
  IF (.NOT. ALLOCATED(agh_tt_fe_sdev)) RETURN
  IF (.NOT. ALLOCATED(agh_tt_fe_logconst)) RETURN
  IF (.NOT. ALLOCATED(agh_tt_fe_modes)) RETURN

  IF (agh_tt_fe_dimEta_cached <= 0 .OR. agh_tt_fe_nvar_cached <= 0 .OR. &
      agh_tt_fe_nquad_eta_cached <= 0 .OR. agh_tt_fe_nquad_re_cached <= 0 .OR. &
      agh_tt_fe_nsub_cached <= 0) RETURN

  dimZ = agh_tt_fe_dimEta_cached + agh_tt_fe_nvar_cached

  IF (SIZE(agh_tt_fe_zeta_eta,1) /= agh_tt_fe_dimEta_cached .OR. &
      SIZE(agh_tt_fe_zeta_eta,2) /= agh_tt_fe_nquad_eta_cached .OR. &
      SIZE(agh_tt_fe_zeta_eta,3) /= agh_tt_fe_nsub_cached) RETURN

  IF (SIZE(agh_tt_fe_alpha,1) /= agh_tt_fe_nvar_cached .OR. &
      SIZE(agh_tt_fe_alpha,2) /= agh_tt_fe_nquad_eta_cached .OR. &
      SIZE(agh_tt_fe_alpha,3) /= agh_tt_fe_nsub_cached) RETURN

  IF (SIZE(agh_tt_fe_sdev,1) /= agh_tt_fe_nvar_cached .OR. &
      SIZE(agh_tt_fe_sdev,2) /= agh_tt_fe_nsub_cached) RETURN

  IF (SIZE(agh_tt_fe_logconst) /= agh_tt_fe_nsub_cached) RETURN

  IF (SIZE(agh_tt_fe_modes,1) /= dimZ .OR. &
      SIZE(agh_tt_fe_modes,2) /= agh_tt_fe_nsub_cached) RETURN

  ready0 = 1_C_INT
  dim_eta0 = INT(agh_tt_fe_dimEta_cached, C_INT)
  dim_z0 = INT(dimZ, C_INT)
  nvar0 = INT(agh_tt_fe_nvar_cached, C_INT)
  nquad_eta0 = INT(agh_tt_fe_nquad_eta_cached, C_INT)
  nquad_re0 = INT(agh_tt_fe_nquad_re_cached, C_INT)
  nsub0 = INT(agh_tt_fe_nsub_cached, C_INT)
  fit_id0 = INT(IEOR(agh_tt_fe_fit_id_cached, INT(Z'40000000')), C_INT)

  idx = 0
  DO i = 1, agh_tt_fe_nsub_cached
     DO l = 1, agh_tt_fe_nquad_eta_cached
        DO d = 1, agh_tt_fe_dimEta_cached
           idx = idx + 1
           zeta_eta_out(idx) = agh_tt_fe_zeta_eta(d,l,i)
        END DO
     END DO
  END DO

  idx = 0
  DO i = 1, agh_tt_fe_nsub_cached
     DO l = 1, agh_tt_fe_nquad_eta_cached
        DO k = 1, agh_tt_fe_nvar_cached
           idx = idx + 1
           alpha_out(idx) = agh_tt_fe_alpha(k,l,i)
        END DO
     END DO
  END DO

  idx = 0
  DO i = 1, agh_tt_fe_nsub_cached
     DO k = 1, agh_tt_fe_nvar_cached
        idx = idx + 1
        sdev_out(idx) = agh_tt_fe_sdev(k,i)
     END DO
  END DO

  DO i = 1, agh_tt_fe_nsub_cached
     logconst_out(i) = agh_tt_fe_logconst(i)
  END DO

  idx = 0
  DO i = 1, agh_tt_fe_nsub_cached
     DO d = 1, dimZ
        idx = idx + 1
        modes_out(idx) = agh_tt_fe_modes(d,i)
     END DO
  END DO
END SUBROUTINE get_growth_agh_tt_fe_grid_copy


SUBROUTINE set_growth_agh_tt_fe_grid_copy(zeta_eta_in, alpha_in, sdev_in, &
                                           logconst_in, modes_in,          &
                                           dim_eta0, dim_z0, nvar0,        &
                                           nquad_eta0, nquad_re0, nsub0,   &
                                           fit_id0, ready0)                &
     BIND(C, name="set_growth_agh_tt_fe_grid_copy")
  USE, INTRINSIC :: ISO_C_BINDING
  USE shared
  IMPLICIT NONE

  REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: zeta_eta_in
  REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: alpha_in
  REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: sdev_in
  REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: logconst_in
  REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: modes_in

  INTEGER(C_INT), INTENT(IN) :: dim_eta0, dim_z0, nvar0
  INTEGER(C_INT), INTENT(IN) :: nquad_eta0, nquad_re0, nsub0, fit_id0
  INTEGER(C_INT), INTENT(OUT) :: ready0

  INTEGER :: dimEtaF, dimZF, nvarF, nquadEtaF, nquadReF, nsubF
  INTEGER :: d, k, l, i, idx

  ready0 = 0_C_INT

  dimEtaF = INT(dim_eta0)
  dimZF = INT(dim_z0)
  nvarF = INT(nvar0)
  nquadEtaF = INT(nquad_eta0)
  nquadReF = INT(nquad_re0)
  nsubF = INT(nsub0)

  IF (dimEtaF <= 0 .OR. nvarF <= 0 .OR. nquadEtaF <= 0 .OR. &
      nquadReF <= 0 .OR. nsubF <= 0) RETURN
  IF (dimZF /= dimEtaF + nvarF) RETURN

  agh_tt_fe_grid_ready = .FALSE.

  IF (ALLOCATED(agh_tt_fe_zeta_eta)) DEALLOCATE(agh_tt_fe_zeta_eta)
  IF (ALLOCATED(agh_tt_fe_alpha)) DEALLOCATE(agh_tt_fe_alpha)
  IF (ALLOCATED(agh_tt_fe_sdev)) DEALLOCATE(agh_tt_fe_sdev)
  IF (ALLOCATED(agh_tt_fe_logconst)) DEALLOCATE(agh_tt_fe_logconst)
  IF (ALLOCATED(agh_tt_fe_modes)) DEALLOCATE(agh_tt_fe_modes)

  ALLOCATE(agh_tt_fe_zeta_eta(dimEtaF, nquadEtaF, nsubF))
  ALLOCATE(agh_tt_fe_alpha(nvarF, nquadEtaF, nsubF))
  ALLOCATE(agh_tt_fe_sdev(nvarF, nsubF))
  ALLOCATE(agh_tt_fe_logconst(nsubF))
  ALLOCATE(agh_tt_fe_modes(dimZF, nsubF))

  idx = 0
  DO i = 1, nsubF
     DO l = 1, nquadEtaF
        DO d = 1, dimEtaF
           idx = idx + 1
           agh_tt_fe_zeta_eta(d,l,i) = zeta_eta_in(idx)
        END DO
     END DO
  END DO

  idx = 0
  DO i = 1, nsubF
     DO l = 1, nquadEtaF
        DO k = 1, nvarF
           idx = idx + 1
           agh_tt_fe_alpha(k,l,i) = alpha_in(idx)
        END DO
     END DO
  END DO

  idx = 0
  DO i = 1, nsubF
     DO k = 1, nvarF
        idx = idx + 1
        agh_tt_fe_sdev(k,i) = sdev_in(idx)
     END DO
  END DO

  DO i = 1, nsubF
     agh_tt_fe_logconst(i) = logconst_in(i)
  END DO

  idx = 0
  DO i = 1, nsubF
     DO d = 1, dimZF
        idx = idx + 1
        agh_tt_fe_modes(d,i) = modes_in(idx)
     END DO
  END DO

  agh_tt_fe_dimEta_cached = dimEtaF
  agh_tt_fe_nvar_cached = nvarF
  agh_tt_fe_nquad_eta_cached = nquadEtaF
  agh_tt_fe_nquad_re_cached = nquadReF
  agh_tt_fe_nsub_cached = nsubF
  agh_tt_fe_fit_id_cached = IEOR(INT(fit_id0), INT(Z'40000000'))
  agh_tt_fe_grid_ready = .TRUE.

  ready0 = 1_C_INT
END SUBROUTINE set_growth_agh_tt_fe_grid_copy


SUBROUTINE reset_growth_agh_tt_fe_grid() BIND(C, name="reset_growth_agh_tt_fe_grid")
  USE shared
  IMPLICIT NONE

  agh_tt_fe_grid_ready = .FALSE.
  agh_tt_fe_fit_id_cached = -999999999
  agh_tt_fe_dimEta_cached = -1
  agh_tt_fe_nvar_cached = -1
  agh_tt_fe_nquad_eta_cached = -1
  agh_tt_fe_nquad_re_cached = -1
  agh_tt_fe_nsub_cached = -1

  IF (ALLOCATED(agh_tt_fe_zeta_eta)) DEALLOCATE(agh_tt_fe_zeta_eta)
  IF (ALLOCATED(agh_tt_fe_alpha)) DEALLOCATE(agh_tt_fe_alpha)
  IF (ALLOCATED(agh_tt_fe_sdev)) DEALLOCATE(agh_tt_fe_sdev)
  IF (ALLOCATED(agh_tt_fe_logconst)) DEALLOCATE(agh_tt_fe_logconst)
  IF (ALLOCATED(agh_tt_fe_modes)) DEALLOCATE(agh_tt_fe_modes)
END SUBROUTINE reset_growth_agh_tt_fe_grid

END MODULE growth_agh_two_tier_equivalent_extension

!=======================================================================
! ADD-ONLY EXTENSION
! Growth likelihood parameterized by standardized growth random effects
! z_u and occasion-specific innovations z_zeta instead of eta directly.
!
! Nested AGH:
!   adaptive outer quadrature in z_u;
!   fixed GH marginalization in z_zeta;
!   adaptive univariate GH for each item random effect.
!
! Full-equivalent AGH:
!   joint mode/Hessian in (z_u,z_zeta,z_a), written in two-tier form;
!   outer block (z_u,z_zeta), conditional univariate item blocks z_a.
!=======================================================================
MODULE growth_re_marginal_agh_extension
  USE shared
  USE utilities, ONLY: alnorm_fast
  USE core, ONLY: cleanup_growth_tt_allocations,                         &
                  setup_load_indices_from_loadmask_growth,               &
                  decode_lambda_sum_constraints_growth,                  &
                  refresh_lambda_act_growth,                             &
                  decode_gamma_sum_constraints_growth,                   &
                  decode_growth_thresholds, decode_growth_beta,           &
                  decode_growth_delta, decode_growth_random_structure,    &
                  growth_n_lambda_par, growth_n_gamma_par,                &
                  growth_n_beta_par, growth_n_zeta_par,                   &
                  build_growth_eta_moments_design,                       &
                  inneragh_item_logint_cond_eta,                          &
                  inneragh_item_logint_full_equiv,                        &
                  blockroot_tt_from_precision_long,                      &
                  tau_at, gamma_at, lambda_at, theta_at
  USE datindiv, ONLY: nmescur
  IMPLICIT NONE

  PRIVATE
   PUBLIC :: loglik_growth_re_tt_agh_nested
   PUBLIC :: loglik_growth_re_tt_agh_full_equiv
   PUBLIC :: get_growth_re_agh_full_equiv_grid_info
   PUBLIC :: get_growth_re_agh_full_equiv_grid_copy

  INTEGER, SAVE :: gr_dim_u = 0
  INTEGER, SAVE :: gr_dim_eta = 0
  INTEGER, SAVE :: gr_dim_w = 0
  INTEGER, SAVE :: gr_nquad_u = 0
  INTEGER, SAVE :: gr_nquad_zeta = 0

  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_map_u(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_zeta_sd(:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_nodes_u(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_weights_u(:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_nodes_zeta(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_weights_zeta(:)

  LOGICAL, SAVE :: gr_nested_ready = .FALSE.
  INTEGER, SAVE :: gr_nested_fit_id = -999999999
  INTEGER, SAVE :: gr_nested_dim_u = -1
  INTEGER, SAVE :: gr_nested_nsub = -1
  INTEGER, SAVE :: gr_nested_nquad_u = -1
  INTEGER, SAVE :: gr_nested_nquad_zeta = -1
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_nested_mode(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_nested_root(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_nested_logdet(:)

  LOGICAL, SAVE :: gr_full_ready = .FALSE.
  INTEGER, SAVE :: gr_full_fit_id = -999999999
  INTEGER, SAVE :: gr_full_dim_w = -1
  INTEGER, SAVE :: gr_full_dim_z = -1
  INTEGER, SAVE :: gr_full_nsub = -1
  INTEGER, SAVE :: gr_full_nquad_u = -1
  INTEGER, SAVE :: gr_full_nquad_zeta = -1
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_full_mode(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_full_B(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_full_logdet_w(:)

  INTERFACE
    SUBROUTINE DPOTRF(UPLO, N, A, LDA, INFO)
      CHARACTER(LEN=*), INTENT(IN) :: UPLO
      INTEGER, INTENT(IN) :: N, LDA
      DOUBLE PRECISION, DIMENSION(LDA,*), INTENT(INOUT) :: A
      INTEGER, INTENT(OUT) :: INFO
    END SUBROUTINE DPOTRF

    SUBROUTINE DPOTRS(UPLO, N, NRHS, A, LDA, B, LDB, INFO)
      CHARACTER(LEN=*), INTENT(IN) :: UPLO
      INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
      DOUBLE PRECISION, DIMENSION(LDA,*), INTENT(IN) :: A
      DOUBLE PRECISION, DIMENSION(LDB,*), INTENT(INOUT) :: B
      INTEGER, INTENT(OUT) :: INFO
    END SUBROUTINE DPOTRS
  END INTERFACE

CONTAINS

SUBROUTINE clear_growth_re_work()
  IMPLICIT NONE

  IF (ALLOCATED(gr_map_u)) DEALLOCATE(gr_map_u)
  IF (ALLOCATED(gr_zeta_sd)) DEALLOCATE(gr_zeta_sd)
  IF (ALLOCATED(gr_nodes_u)) DEALLOCATE(gr_nodes_u)
  IF (ALLOCATED(gr_weights_u)) DEALLOCATE(gr_weights_u)
  IF (ALLOCATED(gr_nodes_zeta)) DEALLOCATE(gr_nodes_zeta)
  IF (ALLOCATED(gr_weights_zeta)) DEALLOCATE(gr_weights_zeta)

  gr_dim_u = 0
  gr_dim_eta = 0
  gr_dim_w = 0
  gr_nquad_u = 0
  gr_nquad_zeta = 0
END SUBROUTINE clear_growth_re_work


SUBROUTINE clear_growth_re_model()
  IMPLICIT NONE

  IF (ALLOCATED(lambda_long)) DEALLOCATE(lambda_long)
  IF (ALLOCATED(gamma_long)) DEALLOCATE(gamma_long)
  IF (ALLOCATED(tau_long)) DEALLOCATE(tau_long)

  CALL cleanup_growth_tt_allocations()
  CALL clear_growth_re_work()
END SUBROUTINE clear_growth_re_model


SUBROUTINE invalidate_nested_cache(fit_id0)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: fit_id0

  IF (fit_id0 /= gr_nested_fit_id) THEN
    gr_nested_ready = .FALSE.
    gr_nested_fit_id = fit_id0
    IF (ALLOCATED(gr_nested_mode)) DEALLOCATE(gr_nested_mode)
    IF (ALLOCATED(gr_nested_root)) DEALLOCATE(gr_nested_root)
    IF (ALLOCATED(gr_nested_logdet)) DEALLOCATE(gr_nested_logdet)
    gr_nested_dim_u = -1
    gr_nested_nsub = -1
    gr_nested_nquad_u = -1
    gr_nested_nquad_zeta = -1
  END IF
END SUBROUTINE invalidate_nested_cache


SUBROUTINE invalidate_full_cache(fit_id0)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: fit_id0

  IF (fit_id0 /= gr_full_fit_id) THEN
    gr_full_ready = .FALSE.
    gr_full_fit_id = fit_id0
    IF (ALLOCATED(gr_full_mode)) DEALLOCATE(gr_full_mode)
    IF (ALLOCATED(gr_full_B)) DEALLOCATE(gr_full_B)
    IF (ALLOCATED(gr_full_logdet_w)) DEALLOCATE(gr_full_logdet_w)
    gr_full_dim_w = -1
    gr_full_dim_z = -1
    gr_full_nsub = -1
    gr_full_nquad_u = -1
    gr_full_nquad_zeta = -1
  END IF
END SUBROUTINE invalidate_full_cache


SUBROUTINE prepare_growth_re_model(                                      &
    b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,           &
    nquad_u0, nodes_u0, weights_u0,                                     &
    nquad_zeta0, nodes_zeta0, weights_zeta0,                            &
    nquad_re0, nodes_re0, weights_re0,                                  &
    nitem0, npm0, re_item0, theta_time0, delta_free0,                   &
    monotone_tau0, Q0, loadmask0,                                     &
    p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0, ok)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0, ntre0
  INTEGER, INTENT(IN) :: nquad_u0, nquad_zeta0, nquad_re0
  INTEGER, INTENT(IN) :: npm0, re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0, Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0, delta_x_col0, delta_z_col0
  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0), INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0), INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_u0
  DOUBLE PRECISION, DIMENSION(nquad_u0), INTENT(IN) :: weights_u0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_zeta0), INTENT(IN) :: weights_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: weights_re0
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: k, t, r, row, col, pos, jj, pos_tau, info
  INTEGER :: i, qidx, s, dimEta, dimU, nX, nZ, rr, zcol
  INTEGER :: ntheta_est, expected_npm, ntre_calc
  INTEGER :: n_lambda_free, n_gamma_free, n_Lu_est
  INTEGER :: nf, h, base
  DOUBLE PRECISION :: tmp_tau, zval
  INTEGER, ALLOCATABLE :: cm3(:)
  DOUBLE PRECISION, ALLOCATABLE :: delta_g(:)
  DOUBLE PRECISION, ALLOCATABLE :: Lu(:,:), Sigma_u(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: zeta_var(:,:)

  ok = .FALSE.

  IF (n0 < 2 .OR. Q0 < 1 .OR. ABS(p_fix0) < 1 .OR. p_ran0 < 1) RETURN
  IF (nquad_u0 < 1 .OR. nquad_zeta0 < 1 .OR. nquad_re0 < 1) RETURN

  CALL clear_growth_re_model()

  nvar = nvar0
  nsub = nsub0
  nobs = nobs0
  n = n0
  Q = Q0
  npm = npm0

  dimEta = n * Q
  dimU = Q * p_ran0
  dim_eta = dimEta

  gr_dim_u = dimU
  gr_dim_eta = dimEta
  gr_dim_w = dimU + dimEta
  gr_nquad_u = nquad_u0
  gr_nquad_zeta = nquad_zeta0

  nquad_re = nquad_re0
  re_item_flag = re_item0
  theta_time_flag = theta_time0
  phi_struct_flag = 0
  growth_subject_specific_design = MERGE(1, 0, p_fix0 > 0)
  growth_design_nsub = MERGE(nsub0, 1, growth_subject_specific_design == 1)
  p_fix = ABS(p_fix0)
  p_ran = p_ran0
  growth_delta_x_col = delta_x_col0
  growth_delta_z_col = delta_z_col0

  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs))
  ALLOCATE(c(nvar), cm1(nvar))
  ALLOCATE(nodes_re(nquad_re), weights_re(nquad_re))
  ALLOCATE(logweights_re(nquad_re))
  ALLOCATE(loadmask(nvar,Q))
  ALLOCATE(lambda(nvar,Q), gamma(nvar))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(theta_long(nvar,n), invsd_long(nvar,n))
  ALLOCATE(sigmak(nvar), sqrt_sigmak(nvar))
  ALLOCATE(kappa(dimEta), L_phi(dimEta,dimEta))
  ALLOCATE(kk_base(nvar))
  ALLOCATE(tau(SUM(cm10)))
  ALLOCATE(cm3(nvar))
  ALLOCATE(delta_g(Q))
  ALLOCATE(beta(Q,p_fix), X_eta(growth_design_nsub,n,Q,p_fix), &
           Z_eta(growth_design_nsub,n,Q,p_ran))
  ALLOCATE(growth_delta_mask(growth_design_nsub,n,Q))
  ALLOCATE(Lu(dimU,dimU), Sigma_u(dimU,dimU))
  ALLOCATE(zeta_var(n,Q))

  ALLOCATE(gr_map_u(dimEta,dimU,growth_design_nsub), gr_zeta_sd(dimEta))
  ALLOCATE(gr_nodes_u(dimU,nquad_u0), gr_weights_u(nquad_u0))
  ALLOCATE(gr_nodes_zeta(dimEta,nquad_zeta0))
  ALLOCATE(gr_weights_zeta(nquad_zeta0))

  nmes = nmes0
  Y = Y0
  nitem = nitem0
  c = c0
  cm1 = cm10
  loadmask = loadmask0

  nX = growth_design_nsub * n * Q * p_fix
  nZ = growth_design_nsub * n * Q * p_ran
  X_eta = RESHAPE(X0(1:nX), (/ growth_design_nsub, n, Q, p_fix /))
  Z_eta = RESHAPE(Z0(1:nZ), (/ growth_design_nsub, n, Q, p_ran /))
  growth_delta_mask = RESHAPE(delta_mask0(1:growth_design_nsub*n*Q), &
                              (/ growth_design_nsub, n, Q /))

  DO r = 1, nquad_u0
    base = (r - 1) * dimU
    gr_nodes_u(:,r) = nodes_u0(base+1:base+dimU)
  END DO
  gr_weights_u = weights_u0

  DO r = 1, nquad_zeta0
    base = (r - 1) * dimEta
    gr_nodes_zeta(:,r) = nodes_zeta0(base+1:base+dimEta)
  END DO
  gr_weights_zeta = weights_zeta0

  nodes_re = nodes_re0
  weights_re = weights_re0
  DO r = 1, nquad_re
    IF (weights_re(r) > 0.d0) THEN
      logweights_re(r) = LOG(weights_re(r))
    ELSE
      logweights_re(r) = -HUGE(1.d0)
    END IF
  END DO

  CALL setup_load_indices_from_loadmask_growth(info)
  IF (info /= 0) GOTO 900

  DO k = 1, nvar
    cm3(k) = MAX(cm1(k) - 2, 0)
  END DO
  ntre_calc = SUM(cm3)
  IF (ntre_calc /= ntre0) GOTO 900
  CALL decode_growth_thresholds(b0, npm0, ntre0, monotone_tau0, info)
  IF (info /= 0) GOTO 900

  pos = ntre0
  CALL decode_lambda_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 900
  CALL refresh_lambda_act_growth()

  CALL decode_gamma_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 900

  theta_k = 1.d0
  invsd_k = 1.d0
  theta_long = 1.d0
  invsd_long = 1.d0

  IF (theta_time_flag == 0) THEN
    DO k = 1, nvar
      IF (cm1(k) > 1) THEN
        pos = pos + 1
        theta_k(k) = EXP(b0(pos))
        invsd_k(k) = 1.d0 / SQRT(theta_k(k))
        DO t = 1, n
          theta_long(k,t) = theta_k(k)
          invsd_long(k,t) = invsd_k(k)
        END DO
      END IF
    END DO
  ELSE
    DO k = 1, nvar
      IF (cm1(k) > 1) THEN
        DO t = 1, n
          pos = pos + 1
          theta_long(k,t) = EXP(b0(pos))
          invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
        END DO
        theta_k(k) = theta_long(k,1)
        invsd_k(k) = invsd_long(k,1)
      END IF
    END DO
  END IF
  CALL decode_growth_beta(b0, npm0, pos, info)
  IF (info /= 0) GOTO 900
  CALL decode_growth_delta(b0, npm0, pos, delta_free0, delta_g, info)
  IF (info /= 0) GOTO 900
  CALL decode_growth_random_structure(b0, npm0, pos, dimU, Lu, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 900

  IF (re_item_flag == 1) THEN
    DO k = 1, nvar
      pos = pos + 1
      sigmak(k) = EXP(b0(pos))
      sqrt_sigmak(k) = SQRT(sigmak(k))
    END DO
  ELSE
    sigmak = 0.d0
    sqrt_sigmak = 0.d0
  END IF

  ntheta_est = 0
  IF (theta_time_flag == 0) THEN
    DO k = 1, nvar
      IF (cm1(k) > 1) ntheta_est = ntheta_est + 1
    END DO
  ELSE
    DO k = 1, nvar
      IF (cm1(k) > 1) ntheta_est = ntheta_est + n
    END DO
  END IF
  n_lambda_free = growth_n_lambda_par()
  n_gamma_free = growth_n_gamma_par()
  n_Lu_est = dimU * (dimU + 1) / 2

  expected_npm = ntre0                                      &
               + n_lambda_free                              &
               + n_gamma_free                               &
               + ntheta_est                                 &
               + growth_n_beta_par()                        &
               + MERGE(Q, 0, delta_free0 == 1)              &
               + n_Lu_est                                   &
               + growth_n_zeta_par()                        &
               + MERGE(nvar, 0, re_item_flag == 1)

  IF (npm0 /= expected_npm) GOTO 900
  IF (pos /= npm0) GOTO 900

  CALL build_growth_eta_moments_design(delta_g, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 900

  gr_map_u = 0.d0
  DO i = 1, growth_design_nsub
    DO t = 1, n
      DO qidx = 1, Q
        row = (t - 1) * Q + qidx
        DO rr = 1, p_ran
          zval = Z_eta(i,t,qidx,rr)
          IF (rr == growth_delta_z_col .AND. growth_delta_mask(i,t,qidx) /= 0.d0) THEN
            IF (growth_delta_mode == 1) THEN
              zval = delta_g(qidx)
            ELSE
              zval = zval + delta_g(qidx) * growth_delta_mask(i,t,qidx)
            END IF
          END IF
          zcol = (rr - 1) * Q + qidx
          DO s = 1, dimU
            gr_map_u(row,s,i) = gr_map_u(row,s,i) + zval * Lu(zcol,s)
          END DO
        END DO
        gr_zeta_sd(row) = SQRT(zeta_var(t,qidx))
      END DO
    END DO
  END DO

  ok = .TRUE.

900 CONTINUE
  IF (ALLOCATED(cm3)) DEALLOCATE(cm3)
  IF (ALLOCATED(delta_g)) DEALLOCATE(delta_g)
  IF (ALLOCATED(Lu)) DEALLOCATE(Lu)
  IF (ALLOCATED(Sigma_u)) DEALLOCATE(Sigma_u)
  IF (ALLOCATED(zeta_var)) DEALLOCATE(zeta_var)

  IF (.NOT. ok) CALL clear_growth_re_model()
END SUBROUTINE prepare_growth_re_model


SUBROUTINE eta_from_u_zeta(i, z_u, z_zeta, eta_vec)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  INTEGER :: design_i
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(IN) :: z_u
  DOUBLE PRECISION, DIMENSION(gr_dim_eta), INTENT(IN) :: z_zeta
  DOUBLE PRECISION, DIMENSION(gr_dim_eta), INTENT(OUT) :: eta_vec

  design_i = MERGE(i, 1, growth_subject_specific_design == 1)
  eta_vec = kappa_sub(:,design_i) + MATMUL(gr_map_u(:,:,design_i), z_u) + &
            gr_zeta_sd * z_zeta
END SUBROUTINE eta_from_u_zeta


SUBROUTINE eta_from_w(i, z_w, eta_vec)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  INTEGER :: design_i
  DOUBLE PRECISION, DIMENSION(gr_dim_w), INTENT(IN) :: z_w
  DOUBLE PRECISION, DIMENSION(gr_dim_eta), INTENT(OUT) :: eta_vec

  design_i = MERGE(i, 1, growth_subject_specific_design == 1)
  eta_vec = kappa_sub(:,design_i) + &
            MATMUL(gr_map_u(:,:,design_i), z_w(1:gr_dim_u)) + &
            gr_zeta_sd * z_w(gr_dim_u+1:gr_dim_w)
END SUBROUTINE eta_from_w


SUBROUTINE log_items_integrated_a(i, eta_vec, logg, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(gr_dim_eta), INTENT(IN) :: eta_vec
  DOUBLE PRECISION, INTENT(OUT) :: logg
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: k_item
  DOUBLE PRECISION :: logint
  LOGICAL :: ok_item

  ok = .TRUE.
  logg = 0.d0

  DO k_item = 1, nvar
    CALL inneragh_item_logint_cond_eta(i, k_item, eta_vec, logint, ok_item)
    IF (.NOT. ok_item .OR. logint /= logint) THEN
      ok = .FALSE.
      logg = -1.d9
      RETURN
    END IF
    logg = logg + logint
  END DO
END SUBROUTINE log_items_integrated_a


SUBROUTINE log_marginal_zeta_given_u(i, z_u, logg, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(IN) :: z_u
  DOUBLE PRECISION, INTENT(OUT) :: logg
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: lz, nvalid
  DOUBLE PRECISION :: logitem, term, maxterm, sumexp
  DOUBLE PRECISION, ALLOCATABLE :: eta_vec(:)
  LOGICAL :: ok_item

  ok = .TRUE.
  logg = -1.d9
  nvalid = 0
  maxterm = 0.d0
  sumexp = 0.d0

  ALLOCATE(eta_vec(gr_dim_eta))

  DO lz = 1, gr_nquad_zeta
    IF (gr_weights_zeta(lz) <= 0.d0) CYCLE

    CALL eta_from_u_zeta(i, z_u, gr_nodes_zeta(:,lz), eta_vec)
    CALL log_items_integrated_a(i, eta_vec, logitem, ok_item)
    IF (.NOT. ok_item) CYCLE

    term = LOG(gr_weights_zeta(lz)) + logitem
    IF (term /= term) CYCLE

    IF (nvalid == 0) THEN
      maxterm = term
      sumexp = 1.d0
      nvalid = 1
    ELSE IF (term > maxterm) THEN
      sumexp = sumexp * EXP(maxterm - term) + 1.d0
      maxterm = term
      nvalid = nvalid + 1
    ELSE
      sumexp = sumexp + EXP(term - maxterm)
      nvalid = nvalid + 1
    END IF
  END DO

  IF (nvalid == 0 .OR. sumexp <= 0.d0 .OR. sumexp /= sumexp) THEN
    ok = .FALSE.
  ELSE
    logg = maxterm + LOG(sumexp)
  END IF

  DEALLOCATE(eta_vec)
END SUBROUTINE log_marginal_zeta_given_u


SUBROUTINE nested_outer_obj(i, z_u, f, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(IN) :: z_u
  DOUBLE PRECISION, INTENT(OUT) :: f
  LOGICAL, INTENT(OUT) :: ok
  DOUBLE PRECISION :: logg

  CALL log_marginal_zeta_given_u(i, z_u, logg, ok)
  IF (.NOT. ok) THEN
    f = 1.d9
  ELSE
    f = -logg + 0.5d0 * DOT_PRODUCT(z_u,z_u)
  END IF
END SUBROUTINE nested_outer_obj


SUBROUTINE nested_fd_grad_hess(i, z, f0, grad, H, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(IN) :: z
  DOUBLE PRECISION, INTENT(OUT) :: f0
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(OUT) :: grad
  DOUBLE PRECISION, DIMENSION(gr_dim_u,gr_dim_u), INTENT(OUT) :: H
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: r, s
  DOUBLE PRECISION :: hr, hs, fp, fm, fpp, fpm, fmp, fmm
  DOUBLE PRECISION, PARAMETER :: epsg = 2.d-4
  DOUBLE PRECISION, PARAMETER :: epsh = 8.d-4
  DOUBLE PRECISION, ALLOCATABLE :: zp(:), zm(:), zpp(:), zpm(:)
  DOUBLE PRECISION, ALLOCATABLE :: zmp(:), zmm(:)
  LOGICAL :: ok0, okp, okm, okpp, okpm, okmp, okmm

  ALLOCATE(zp(gr_dim_u), zm(gr_dim_u))
  ALLOCATE(zpp(gr_dim_u), zpm(gr_dim_u))
  ALLOCATE(zmp(gr_dim_u), zmm(gr_dim_u))

  CALL nested_outer_obj(i, z, f0, ok0)
  IF (.NOT. ok0 .OR. f0 /= f0) THEN
    ok = .FALSE.
    GOTO 800
  END IF

  grad = 0.d0
  H = 0.d0

  DO r = 1, gr_dim_u
    hr = epsg * MAX(ABS(z(r)), 1.d0)
    zp = z
    zm = z
    zp(r) = zp(r) + hr
    zm(r) = zm(r) - hr
    CALL nested_outer_obj(i, zp, fp, okp)
    CALL nested_outer_obj(i, zm, fm, okm)
    IF (.NOT. okp .OR. .NOT. okm) THEN
      ok = .FALSE.
      GOTO 800
    END IF
    grad(r) = (fp - fm) / (2.d0 * hr)

    hr = epsh * MAX(ABS(z(r)), 1.d0)
    zp = z
    zm = z
    zp(r) = zp(r) + hr
    zm(r) = zm(r) - hr
    CALL nested_outer_obj(i, zp, fp, okp)
    CALL nested_outer_obj(i, zm, fm, okm)
    IF (.NOT. okp .OR. .NOT. okm) THEN
      ok = .FALSE.
      GOTO 800
    END IF
    H(r,r) = (fp - 2.d0*f0 + fm) / (hr*hr)
  END DO

  DO r = 1, gr_dim_u - 1
    hr = epsh * MAX(ABS(z(r)), 1.d0)
    DO s = r + 1, gr_dim_u
      hs = epsh * MAX(ABS(z(s)), 1.d0)
      zpp = z
      zpm = z
      zmp = z
      zmm = z
      zpp(r) = zpp(r) + hr
      zpp(s) = zpp(s) + hs
      zpm(r) = zpm(r) + hr
      zpm(s) = zpm(s) - hs
      zmp(r) = zmp(r) - hr
      zmp(s) = zmp(s) + hs
      zmm(r) = zmm(r) - hr
      zmm(s) = zmm(s) - hs
      CALL nested_outer_obj(i, zpp, fpp, okpp)
      CALL nested_outer_obj(i, zpm, fpm, okpm)
      CALL nested_outer_obj(i, zmp, fmp, okmp)
      CALL nested_outer_obj(i, zmm, fmm, okmm)
      IF (.NOT. okpp .OR. .NOT. okpm .OR. .NOT. okmp .OR. .NOT. okmm) THEN
        ok = .FALSE.
        GOTO 800
      END IF
      H(r,s) = (fpp - fpm - fmp + fmm) / (4.d0*hr*hs)
      H(s,r) = H(r,s)
    END DO
  END DO

  ok = .TRUE.

800 CONTINUE
  DEALLOCATE(zp, zm, zpp, zpm, zmp, zmm)
END SUBROUTINE nested_fd_grad_hess


SUBROUTINE solve_spd_step(Hin, grad, step, ok)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Hin
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: grad
  DOUBLE PRECISION, DIMENSION(SIZE(grad)), INTENT(OUT) :: step
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: d, info
  DOUBLE PRECISION, ALLOCATABLE :: A(:,:), rhs(:,:)

  d = SIZE(grad)
  ALLOCATE(A(d,d), rhs(d,1))
  ok = .FALSE.

  A = 0.5d0 * (Hin + TRANSPOSE(Hin))
  CALL DPOTRF('U', d, A, d, info)
  IF (info == 0) THEN
    rhs(:,1) = -grad
    CALL DPOTRS('U', d, 1, A, d, rhs, d, info)
    IF (info == 0 .AND. ALL(rhs(:,1) == rhs(:,1))) THEN
      step = rhs(:,1)
      ok = .TRUE.
    END IF
  END IF

  DEALLOCATE(A, rhs)
END SUBROUTINE solve_spd_step


SUBROUTINE invert_upper(U, Uinv, ok)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: U
  DOUBLE PRECISION, DIMENSION(SIZE(U,1),SIZE(U,2)), INTENT(OUT) :: Uinv
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: d, i, j, k
  DOUBLE PRECISION :: rhs

  d = SIZE(U,1)
  Uinv = 0.d0
  ok = .FALSE.
  IF (SIZE(U,2) /= d) RETURN

  DO j = 1, d
    DO i = d, 1, -1
      IF (ABS(U(i,i)) <= TINY(1.d0) .OR. U(i,i) /= U(i,i)) RETURN
      rhs = 0.d0
      IF (i == j) rhs = 1.d0
      DO k = i + 1, d
        rhs = rhs - U(i,k) * Uinv(k,j)
      END DO
      Uinv(i,j) = rhs / U(i,i)
    END DO
  END DO

  ok = .TRUE.
END SUBROUTINE invert_upper


SUBROUTINE nested_find_mode_root(i, start, mode, R, logdetR, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(IN) :: start
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(OUT) :: mode
  DOUBLE PRECISION, DIMENSION(gr_dim_u,gr_dim_u), INTENT(OUT) :: R
  DOUBLE PRECISION, INTENT(OUT) :: logdetR
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: iter, ls, info, j
  DOUBLE PRECISION :: f, fnew, alpha, gnorm
  DOUBLE PRECISION, ALLOCATABLE :: grad(:), step(:), znew(:)
  DOUBLE PRECISION, ALLOCATABLE :: H(:,:), A(:,:)
  LOGICAL :: okgh, okstep, oknew

  ALLOCATE(grad(gr_dim_u), step(gr_dim_u), znew(gr_dim_u))
  ALLOCATE(H(gr_dim_u,gr_dim_u), A(gr_dim_u,gr_dim_u))

  mode = start
  ok = .FALSE.

  DO iter = 1, 35
    CALL nested_fd_grad_hess(i, mode, f, grad, H, okgh)
    IF (.NOT. okgh) GOTO 900

    gnorm = SQRT(DOT_PRODUCT(grad,grad))
    IF (gnorm < 2.d-5) EXIT

    CALL solve_spd_step(H, grad, step, okstep)
    IF (.NOT. okstep) step = -grad / MAX(gnorm,1.d0)

    alpha = 1.d0
    oknew = .FALSE.
    DO ls = 1, 14
      znew = mode + alpha * step
      CALL nested_outer_obj(i, znew, fnew, oknew)
      IF (oknew .AND. fnew < f) EXIT
      alpha = alpha * 0.5d0
    END DO

    IF (.NOT. oknew .OR. fnew >= f) THEN
      step = -grad / MAX(gnorm,1.d0)
      alpha = 0.25d0
      DO ls = 1, 14
        znew = mode + alpha * step
        CALL nested_outer_obj(i, znew, fnew, oknew)
        IF (oknew .AND. fnew < f) EXIT
        alpha = alpha * 0.5d0
      END DO
    END IF

    IF (.NOT. oknew .OR. fnew >= f) EXIT
    mode = znew
    IF (SQRT(DOT_PRODUCT(alpha*step,alpha*step)) < 1.d-6) EXIT
  END DO

  CALL nested_fd_grad_hess(i, mode, f, grad, H, okgh)
  IF (.NOT. okgh) GOTO 900

  A = 0.5d0 * (H + TRANSPOSE(H))
  CALL DPOTRF('U', gr_dim_u, A, gr_dim_u, info)
  IF (info /= 0) GOTO 900

  logdetR = 0.d0
  DO j = 1, gr_dim_u
    IF (A(j,j) <= 0.d0) GOTO 900
    logdetR = logdetR - LOG(A(j,j))
  END DO

  IF (gr_dim_u > 1) THEN
    DO j = 1, gr_dim_u-1
      A(j+1:gr_dim_u,j) = 0.d0
    END DO
  END IF
  CALL invert_upper(A, R, okgh)
  IF (.NOT. okgh) GOTO 900

  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(grad, step, znew, H, A)
END SUBROUTINE nested_find_mode_root


SUBROUTINE build_nested_cache(ok)
  IMPLICIT NONE
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: i, sm
  DOUBLE PRECISION, ALLOCATABLE :: start(:), mode(:), R(:,:)
  DOUBLE PRECISION :: logdetR
  LOGICAL :: ok_i

  ok = .FALSE.

  IF (ALLOCATED(gr_nested_mode)) THEN
    IF (gr_nested_dim_u /= gr_dim_u .OR. gr_nested_nsub /= nsub .OR. &
        gr_nested_nquad_u /= gr_nquad_u .OR.                      &
        gr_nested_nquad_zeta /= gr_nquad_zeta) THEN
      DEALLOCATE(gr_nested_mode, gr_nested_root, gr_nested_logdet)
      gr_nested_ready = .FALSE.
    END IF
  END IF

  IF (.NOT. ALLOCATED(gr_nested_mode)) THEN
    ALLOCATE(gr_nested_mode(gr_dim_u,nsub))
    ALLOCATE(gr_nested_root(gr_dim_u,gr_dim_u,nsub))
    ALLOCATE(gr_nested_logdet(nsub))
    gr_nested_mode = 0.d0
    gr_nested_root = 0.d0
    gr_nested_logdet = 0.d0
    gr_nested_dim_u = gr_dim_u
    gr_nested_nsub = nsub
    gr_nested_nquad_u = gr_nquad_u
    gr_nested_nquad_zeta = gr_nquad_zeta
  END IF

  ALLOCATE(start(gr_dim_u), mode(gr_dim_u), R(gr_dim_u,gr_dim_u))
  nmescur = 0

  DO i = 1, nsub
    sm = SUM(nmes(i,:))
    IF (gr_nested_ready) THEN
      start = gr_nested_mode(:,i)
    ELSE
      start = 0.d0
    END IF

    CALL nested_find_mode_root(i, start, mode, R, logdetR, ok_i)
    IF (.NOT. ok_i) THEN
      gr_nested_ready = .FALSE.
      GOTO 900
    END IF

    gr_nested_mode(:,i) = mode
    gr_nested_root(:,:,i) = R
    gr_nested_logdet(i) = logdetR
    nmescur = nmescur + sm
  END DO

  gr_nested_ready = .TRUE.
  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(start, mode, R)
END SUBROUTINE build_nested_cache


SUBROUTINE eval_nested_loglik(out, ok)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(OUT) :: out
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: i, lu, sm, nvalid
  DOUBLE PRECISION :: logg, term, maxterm, sumexp, qi, qv
  DOUBLE PRECISION, ALLOCATABLE :: z_u(:), v_u(:)
  LOGICAL :: okg

  out = 0.d0
  ok = .FALSE.
  IF (.NOT. gr_nested_ready) RETURN

  ALLOCATE(z_u(gr_dim_u), v_u(gr_dim_u))
  nmescur = 0

  DO i = 1, nsub
    sm = SUM(nmes(i,:))
    nvalid = 0
    maxterm = 0.d0
    sumexp = 0.d0

    DO lu = 1, gr_nquad_u
      IF (gr_weights_u(lu) <= 0.d0) CYCLE
      v_u = gr_nodes_u(:,lu)
      z_u = gr_nested_mode(:,i) + MATMUL(gr_nested_root(:,:,i),v_u)

      CALL log_marginal_zeta_given_u(i, z_u, logg, okg)
      IF (.NOT. okg) CYCLE

      qi = DOT_PRODUCT(z_u,z_u)
      qv = DOT_PRODUCT(v_u,v_u)
      term = LOG(gr_weights_u(lu)) + logg - 0.5d0*(qi-qv)
      IF (term /= term) CYCLE

      IF (nvalid == 0) THEN
        maxterm = term
        sumexp = 1.d0
        nvalid = 1
      ELSE IF (term > maxterm) THEN
        sumexp = sumexp * EXP(maxterm-term) + 1.d0
        maxterm = term
        nvalid = nvalid + 1
      ELSE
        sumexp = sumexp + EXP(term-maxterm)
        nvalid = nvalid + 1
      END IF
    END DO

    IF (nvalid == 0 .OR. sumexp <= 0.d0) GOTO 900
    out = out + gr_nested_logdet(i) + maxterm + LOG(sumexp)
    nmescur = nmescur + sm
  END DO

  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(z_u, v_u)
END SUBROUTINE eval_nested_loglik


SUBROUTINE full_grad_hess(i, z, f, grad, H, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: z
  DOUBLE PRECISION, INTENT(OUT) :: f
  DOUBLE PRECISION, DIMENSION(SIZE(z)), INTENT(OUT) :: grad
  DOUBLE PRECISION, DIMENSION(SIZE(z),SIZE(z)), INTENT(OUT) :: H
  LOGICAL, INTENT(OUT) :: ok

  DOUBLE PRECISION, PARAMETER :: INV_SQRT_2PI = &
       0.39894228040143267794d0
  DOUBLE PRECISION, PARAMETER :: P_FLOOR = 1.d-300

  INTEGER :: dimZ, k_item, j_mes, idx_obs, yijk, tuk
  INTEGER :: t_idx, kk_thr, qq, row, d1, d2, design_i
  DOUBLE PRECISION :: mu, linpred, var_e, invsd
  DOUBLE PRECISION :: aa, bb, p, phia, phib
  DOUBLE PRECISION :: dp, d2p, smu, hmu, loglik
  DOUBLE PRECISION, ALLOCATABLE :: eta_vec(:), a_item(:), dmu(:)
  DOUBLE PRECISION, ALLOCATABLE :: z_w(:)

  dimZ = SIZE(z)
  design_i = MERGE(i, 1, growth_subject_specific_design == 1)
  ok = .FALSE.
  f = 1.d9
  grad = 0.d0
  H = 0.d0

  IF (re_item_flag == 1) THEN
    IF (dimZ /= gr_dim_w + nvar) RETURN
  ELSE
    IF (dimZ /= gr_dim_w) RETURN
  END IF

  ALLOCATE(eta_vec(gr_dim_eta), a_item(nvar), dmu(dimZ))
  ALLOCATE(z_w(gr_dim_w))

  z_w = z(1:gr_dim_w)
  CALL eta_from_w(i, z_w, eta_vec)

  IF (re_item_flag == 1) THEN
    DO k_item = 1, nvar
      a_item(k_item) = sqrt_sigmak(k_item) * z(gr_dim_w+k_item)
    END DO
  ELSE
    a_item = 0.d0
  END IF

  loglik = 0.d0
  kk_thr = 0
  idx_obs = nmescur

  DO k_item = 1, nvar
    DO j_mes = 1, nmes(i,k_item)
      idx_obs = idx_obs + 1
      t_idx = nitem(idx_obs)
      var_e = theta_at(k_item,t_idx)
      IF (var_e <= 0.d0 .OR. var_e /= var_e) GOTO 900
      invsd = 1.d0 / SQRT(var_e)

      linpred = 0.d0
      DO qq = 1, Q
        IF (loadmask(k_item,qq) == 0) CYCLE
        row = (t_idx-1)*Q + qq
        linpred = linpred + lambda_at(k_item,qq,t_idx) * eta_vec(row)
      END DO
      mu = gamma_at(k_item,t_idx) + a_item(k_item) + linpred

      yijk = Y(idx_obs)
      tuk = yijk + kk_thr

      IF (c(k_item) == 2) THEN
        aa = (tau_at(kk_thr+1,t_idx)-mu) * invsd
        phia = INV_SQRT_2PI * EXP(-0.5d0*aa*aa)
        IF (yijk == 1) THEN
          p = alnorm_fast(aa,.FALSE.)
          dp = -phia*invsd
          d2p = -aa*phia*invsd*invsd
        ELSE
          p = 1.d0-alnorm_fast(aa,.FALSE.)
          dp = phia*invsd
          d2p = aa*phia*invsd*invsd
        END IF
      ELSE IF (yijk == 1) THEN
        aa = (tau_at(tuk,t_idx)-mu) * invsd
        p = alnorm_fast(aa,.FALSE.)
        phia = INV_SQRT_2PI * EXP(-0.5d0*aa*aa)
        dp = -phia*invsd
        d2p = -aa*phia*invsd*invsd
      ELSE IF (yijk == c(k_item)) THEN
        bb = (tau_at(tuk-1,t_idx)-mu) * invsd
        p = 1.d0-alnorm_fast(bb,.FALSE.)
        phib = INV_SQRT_2PI * EXP(-0.5d0*bb*bb)
        dp = phib*invsd
        d2p = bb*phib*invsd*invsd
      ELSE
        aa = (tau_at(tuk,t_idx)-mu) * invsd
        bb = (tau_at(tuk-1,t_idx)-mu) * invsd
        p = alnorm_fast(aa,.FALSE.)-alnorm_fast(bb,.FALSE.)
        phia = INV_SQRT_2PI * EXP(-0.5d0*aa*aa)
        phib = INV_SQRT_2PI * EXP(-0.5d0*bb*bb)
        dp = (phib-phia)*invsd
        d2p = (bb*phib-aa*phia)*invsd*invsd
      END IF

      IF (p <= 0.d0 .OR. p /= p) GOTO 900
      p = MAX(p,P_FLOOR)
      loglik = loglik + LOG(p)
      smu = dp/p
      hmu = d2p/p-smu*smu
      IF (smu /= smu .OR. hmu /= hmu) GOTO 900

      dmu = 0.d0
      DO qq = 1, Q
        IF (loadmask(k_item,qq) == 0) CYCLE
        row = (t_idx-1)*Q + qq
        DO d1 = 1, gr_dim_u
          dmu(d1) = dmu(d1) + &
               lambda_at(k_item,qq,t_idx)*gr_map_u(row,d1,design_i)
        END DO
        d1 = gr_dim_u + row
        dmu(d1) = dmu(d1) + &
             lambda_at(k_item,qq,t_idx)*gr_zeta_sd(row)
      END DO
      IF (re_item_flag == 1) THEN
        dmu(gr_dim_w+k_item) = sqrt_sigmak(k_item)
      END IF

      grad = grad + smu*dmu
      DO d1 = 1, dimZ
        IF (dmu(d1) == 0.d0) CYCLE
        DO d2 = 1, dimZ
          IF (dmu(d2) == 0.d0) CYCLE
          H(d1,d2) = H(d1,d2) + hmu*dmu(d1)*dmu(d2)
        END DO
      END DO
    END DO
    kk_thr = kk_thr + c(k_item)-1
  END DO

  f = -loglik + 0.5d0*DOT_PRODUCT(z,z)
  grad = -grad + z
  H = -H
  DO d1 = 1, dimZ
    H(d1,d1) = H(d1,d1) + 1.d0
  END DO
  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(eta_vec, a_item, dmu, z_w)
END SUBROUTINE full_grad_hess


SUBROUTINE full_find_mode_B(i, start, mode, B, logdetW, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: start
  DOUBLE PRECISION, DIMENSION(SIZE(start)), INTENT(OUT) :: mode
  DOUBLE PRECISION, DIMENSION(SIZE(start),SIZE(start)), INTENT(OUT) :: B
  DOUBLE PRECISION, INTENT(OUT) :: logdetW
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: dimZ, iter, ls
  DOUBLE PRECISION :: f, fnew, alpha, gnorm, logdetB
  DOUBLE PRECISION, ALLOCATABLE :: grad(:), step(:), znew(:)
  DOUBLE PRECISION, ALLOCATABLE :: H(:,:), Hreg(:,:)
  LOGICAL :: okgh, okstep, oknew, okroot

  dimZ = SIZE(start)
  ALLOCATE(grad(dimZ), step(dimZ), znew(dimZ))
  ALLOCATE(H(dimZ,dimZ), Hreg(dimZ,dimZ))

  mode = start
  ok = .FALSE.

  DO iter = 1, 50
    CALL full_grad_hess(i, mode, f, grad, H, okgh)
    IF (.NOT. okgh) GOTO 900
    gnorm = SQRT(DOT_PRODUCT(grad,grad))
    IF (gnorm < 1.d-6) EXIT

    CALL solve_spd_step(H, grad, step, okstep)
    IF (.NOT. okstep) step = -grad/MAX(gnorm,1.d0)

    alpha = 1.d0
    oknew = .FALSE.
    DO ls = 1, 16
      znew = mode + alpha*step
      CALL full_grad_hess(i, znew, fnew, grad, Hreg, oknew)
      IF (oknew .AND. fnew < f) EXIT
      alpha = alpha*0.5d0
    END DO

    IF (.NOT. oknew .OR. fnew >= f) THEN
      CALL full_grad_hess(i, mode, f, grad, H, okgh)
      step = -grad/MAX(SQRT(DOT_PRODUCT(grad,grad)),1.d0)
      alpha = 0.25d0
      DO ls = 1, 16
        znew = mode + alpha*step
        CALL full_grad_hess(i, znew, fnew, grad, Hreg, oknew)
        IF (oknew .AND. fnew < f) EXIT
        alpha = alpha*0.5d0
      END DO
    END IF

    IF (.NOT. oknew .OR. fnew >= f) EXIT
    mode = znew
    IF (SQRT(DOT_PRODUCT(alpha*step,alpha*step)) < 1.d-7) EXIT
  END DO

  CALL full_grad_hess(i, mode, f, grad, H, okgh)
  IF (.NOT. okgh) GOTO 900

  Hreg = 0.5d0*(H+TRANSPOSE(H))
  CALL blockroot_tt_from_precision_long(gr_dim_w, dimZ, Hreg, &
       B, logdetB, logdetW, okroot)
  IF (.NOT. okroot) GOTO 900

  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(grad, step, znew, H, Hreg)
END SUBROUTINE full_find_mode_B


SUBROUTINE build_full_cache(ok)
  IMPLICIT NONE
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: i, sm, dimZ
  DOUBLE PRECISION, ALLOCATABLE :: start(:), mode(:), B(:,:)
  DOUBLE PRECISION :: logdetW
  LOGICAL :: ok_i

  ok = .FALSE.
  dimZ = gr_dim_w + MERGE(nvar,0,re_item_flag == 1)

  IF (ALLOCATED(gr_full_mode)) THEN
    IF (gr_full_dim_w /= gr_dim_w .OR. gr_full_dim_z /= dimZ .OR. &
        gr_full_nsub /= nsub .OR.                                &
        gr_full_nquad_u /= gr_nquad_u .OR.                       &
        gr_full_nquad_zeta /= gr_nquad_zeta) THEN
      DEALLOCATE(gr_full_mode, gr_full_B, gr_full_logdet_w)
      gr_full_ready = .FALSE.
    END IF
  END IF

  IF (.NOT. ALLOCATED(gr_full_mode)) THEN
    ALLOCATE(gr_full_mode(dimZ,nsub))
    ALLOCATE(gr_full_B(dimZ,dimZ,nsub))
    ALLOCATE(gr_full_logdet_w(nsub))
    gr_full_mode = 0.d0
    gr_full_B = 0.d0
    gr_full_logdet_w = 0.d0
    gr_full_dim_w = gr_dim_w
    gr_full_dim_z = dimZ
    gr_full_nsub = nsub
    gr_full_nquad_u = gr_nquad_u
    gr_full_nquad_zeta = gr_nquad_zeta
  END IF

  ALLOCATE(start(dimZ), mode(dimZ), B(dimZ,dimZ))
  nmescur = 0

  DO i = 1, nsub
    sm = SUM(nmes(i,:))
    IF (gr_full_ready) THEN
      start = gr_full_mode(:,i)
    ELSE
      start = 0.d0
    END IF

    CALL full_find_mode_B(i, start, mode, B, logdetW, ok_i)
    IF (.NOT. ok_i) THEN
      gr_full_ready = .FALSE.
      GOTO 900
    END IF

    gr_full_mode(:,i) = mode
    gr_full_B(:,:,i) = B
    gr_full_logdet_w(i) = logdetW
    nmescur = nmescur + sm
  END DO

  gr_full_ready = .TRUE.
  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(start, mode, B)
END SUBROUTINE build_full_cache


SUBROUTINE eval_full_loglik(out, ok)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(OUT) :: out
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: i, lu, lz, k_item, sm, nvalid, dimZ
  DOUBLE PRECISION :: logg, logint, term, maxterm, sumexp
  DOUBLE PRECISION :: qw, qv, alpha_k, sdev_k
  DOUBLE PRECISION, ALLOCATABLE :: v_w(:), z_w(:), eta_vec(:)
  LOGICAL :: ok_item

  out = 0.d0
  ok = .FALSE.
  IF (.NOT. gr_full_ready) RETURN

  dimZ = gr_full_dim_z
  ALLOCATE(v_w(gr_dim_w), z_w(gr_dim_w), eta_vec(gr_dim_eta))
  nmescur = 0

  DO i = 1, nsub
    sm = SUM(nmes(i,:))
    nvalid = 0
    maxterm = 0.d0
    sumexp = 0.d0

    DO lu = 1, gr_nquad_u
      IF (gr_weights_u(lu) <= 0.d0) CYCLE
      DO lz = 1, gr_nquad_zeta
        IF (gr_weights_zeta(lz) <= 0.d0) CYCLE

        v_w(1:gr_dim_u) = gr_nodes_u(:,lu)
        v_w(gr_dim_u+1:gr_dim_w) = gr_nodes_zeta(:,lz)
        z_w = gr_full_mode(1:gr_dim_w,i) + &
              MATMUL(gr_full_B(1:gr_dim_w,1:gr_dim_w,i),v_w)
        CALL eta_from_w(i, z_w, eta_vec)

        logg = 0.d0
        DO k_item = 1, nvar
          IF (re_item_flag == 1) THEN
            alpha_k = gr_full_mode(gr_dim_w+k_item,i) + &
                 DOT_PRODUCT(gr_full_B(gr_dim_w+k_item,1:gr_dim_w,i),v_w)
            sdev_k = gr_full_B(gr_dim_w+k_item,gr_dim_w+k_item,i)
          ELSE
            alpha_k = 0.d0
            sdev_k = 1.d0
          END IF

          CALL inneragh_item_logint_full_equiv(i, k_item, eta_vec, &
               alpha_k, sdev_k, logint, ok_item)
          IF (.NOT. ok_item) THEN
            logg = -1.d9
            EXIT
          END IF
          logg = logg + logint
        END DO
        IF (logg <= -1.d9) CYCLE

        qw = DOT_PRODUCT(z_w,z_w)
        qv = DOT_PRODUCT(v_w,v_w)
        term = LOG(gr_weights_u(lu)) + LOG(gr_weights_zeta(lz)) + &
               logg - 0.5d0*(qw-qv)
        IF (term /= term) CYCLE

        IF (nvalid == 0) THEN
          maxterm = term
          sumexp = 1.d0
          nvalid = 1
        ELSE IF (term > maxterm) THEN
          sumexp = sumexp*EXP(maxterm-term)+1.d0
          maxterm = term
          nvalid = nvalid+1
        ELSE
          sumexp = sumexp+EXP(term-maxterm)
          nvalid = nvalid+1
        END IF
      END DO
    END DO

    IF (nvalid == 0 .OR. sumexp <= 0.d0) GOTO 900
    out = out + gr_full_logdet_w(i) + maxterm + LOG(sumexp)
    nmescur = nmescur + sm
  END DO

  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(v_w, z_w, eta_vec)
END SUBROUTINE eval_full_loglik


SUBROUTINE loglik_growth_re_tt_agh_nested(                                &
    b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,            &
    nquad_u0, nodes_u0, weights_u0,                                      &
    nquad_zeta0, nodes_zeta0, weights_zeta0,                             &
    nquad_re0, nodes_re0, weights_re0,                                   &
    nitem0, npm0, re_item0, theta_time0, delta_free0,                    &
    monotone_tau0, update_gh0, fit_id0, Q0, loadmask0,             &
    p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0, vrais)            &
    BIND(C, name="loglik_growth_re_tt_agh_nested")

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0, ntre0
  INTEGER, INTENT(IN) :: nquad_u0, nquad_zeta0, nquad_re0
  INTEGER, INTENT(IN) :: npm0, re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0
  INTEGER, INTENT(IN) :: update_gh0, fit_id0, Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0, delta_x_col0, delta_z_col0
  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0), INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0), INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_u0
  DOUBLE PRECISION, DIMENSION(nquad_u0), INTENT(IN) :: weights_u0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_zeta0), INTENT(IN) :: weights_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: weights_re0
  DOUBLE PRECISION, INTENT(OUT) :: vrais

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.d9
  DOUBLE PRECISION :: ll
  LOGICAL :: ok_model, ok_grid, ok_ll

  vrais = BAD_OBJ
  CALL invalidate_nested_cache(fit_id0)

  CALL prepare_growth_re_model(                                       &
       b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,      &
       nquad_u0, nodes_u0, weights_u0,                                &
       nquad_zeta0, nodes_zeta0, weights_zeta0,                       &
       nquad_re0, nodes_re0, weights_re0,                             &
       nitem0, npm0, re_item0, theta_time0, delta_free0,              &
       monotone_tau0, Q0, loadmask0,                              &
       p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0, ok_model)
  IF (.NOT. ok_model) GOTO 900

  IF (update_gh0 == 1 .OR. .NOT. gr_nested_ready) THEN
    CALL build_nested_cache(ok_grid)
    IF (.NOT. ok_grid) GOTO 900
  END IF

  CALL eval_nested_loglik(ll, ok_ll)
  IF (ok_ll .AND. ll == ll .AND. ABS(ll) < 1.d300) vrais = -ll

900 CONTINUE
  CALL clear_growth_re_model()
END SUBROUTINE loglik_growth_re_tt_agh_nested


SUBROUTINE loglik_growth_re_tt_agh_full_equiv(                             &
    b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,             &
    nquad_u0, nodes_u0, weights_u0,                                       &
    nquad_zeta0, nodes_zeta0, weights_zeta0,                              &
    nquad_re0, nodes_re0, weights_re0,                                    &
    nitem0, npm0, re_item0, theta_time0, delta_free0,                     &
    monotone_tau0, update_gh0, fit_id0, Q0, loadmask0,             &
    p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0, vrais)             &
    BIND(C, name="loglik_growth_re_tt_agh_full_equiv")

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0, ntre0
  INTEGER, INTENT(IN) :: nquad_u0, nquad_zeta0, nquad_re0
  INTEGER, INTENT(IN) :: npm0, re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0
  INTEGER, INTENT(IN) :: update_gh0, fit_id0, Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0, delta_x_col0, delta_z_col0
  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0), INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0), INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_u0
  DOUBLE PRECISION, DIMENSION(nquad_u0), INTENT(IN) :: weights_u0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_zeta0), INTENT(IN) :: weights_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: weights_re0
  DOUBLE PRECISION, INTENT(OUT) :: vrais

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.d9
  DOUBLE PRECISION :: ll
  LOGICAL :: ok_model, ok_grid, ok_ll

  vrais = BAD_OBJ
  CALL invalidate_full_cache(fit_id0)

  CALL prepare_growth_re_model(                                       &
       b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,      &
       nquad_u0, nodes_u0, weights_u0,                                &
       nquad_zeta0, nodes_zeta0, weights_zeta0,                       &
       nquad_re0, nodes_re0, weights_re0,                             &
       nitem0, npm0, re_item0, theta_time0, delta_free0,              &
       monotone_tau0, Q0, loadmask0,                              &
       p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0, ok_model)
  IF (.NOT. ok_model) GOTO 900

  IF (update_gh0 == 1 .OR. .NOT. gr_full_ready) THEN
    CALL build_full_cache(ok_grid)
    IF (.NOT. ok_grid) GOTO 900
  END IF

  CALL eval_full_loglik(ll, ok_ll)
  IF (ok_ll .AND. ll == ll .AND. ABS(ll) < 1.d300) vrais = -ll

900 CONTINUE
  CALL clear_growth_re_model()
END SUBROUTINE loglik_growth_re_tt_agh_full_equiv


SUBROUTINE get_growth_re_agh_full_equiv_grid_info(                       &
    ready0, fit_id0, dim_w0, dim_z0, nsub0, nquad_u0, nquad_zeta0)       &
    BIND(C, name="get_growth_re_agh_full_equiv_grid_info")

  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE

  INTEGER(C_INT), INTENT(OUT) :: ready0
  INTEGER(C_INT), INTENT(OUT) :: fit_id0
  INTEGER(C_INT), INTENT(OUT) :: dim_w0
  INTEGER(C_INT), INTENT(OUT) :: dim_z0
  INTEGER(C_INT), INTENT(OUT) :: nsub0
  INTEGER(C_INT), INTENT(OUT) :: nquad_u0
  INTEGER(C_INT), INTENT(OUT) :: nquad_zeta0

  ready0 = MERGE(1_C_INT, 0_C_INT, gr_full_ready)
  fit_id0 = INT(gr_full_fit_id, C_INT)
  dim_w0 = INT(gr_full_dim_w, C_INT)
  dim_z0 = INT(gr_full_dim_z, C_INT)
  nsub0 = INT(gr_full_nsub, C_INT)
  nquad_u0 = INT(gr_full_nquad_u, C_INT)
  nquad_zeta0 = INT(gr_full_nquad_zeta, C_INT)

END SUBROUTINE get_growth_re_agh_full_equiv_grid_info


SUBROUTINE get_growth_re_agh_full_equiv_grid_copy(                       &
    mode_out, B_out, logdet_w_out,                                       &
    ready0, fit_id0, dim_w0, dim_z0, nsub0, nquad_u0, nquad_zeta0)       &
    BIND(C, name="get_growth_re_agh_full_equiv_grid_copy")

  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE

  REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT) :: mode_out
  REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT) :: B_out
  REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT) :: logdet_w_out

  INTEGER(C_INT), INTENT(OUT) :: ready0
  INTEGER(C_INT), INTENT(OUT) :: fit_id0
  INTEGER(C_INT), INTENT(OUT) :: dim_w0
  INTEGER(C_INT), INTENT(OUT) :: dim_z0
  INTEGER(C_INT), INTENT(OUT) :: nsub0
  INTEGER(C_INT), INTENT(OUT) :: nquad_u0
  INTEGER(C_INT), INTENT(OUT) :: nquad_zeta0

  INTEGER :: i, j, k, pos

  ready0 = 0_C_INT
  fit_id0 = INT(gr_full_fit_id, C_INT)
  dim_w0 = INT(gr_full_dim_w, C_INT)
  dim_z0 = INT(gr_full_dim_z, C_INT)
  nsub0 = INT(gr_full_nsub, C_INT)
  nquad_u0 = INT(gr_full_nquad_u, C_INT)
  nquad_zeta0 = INT(gr_full_nquad_zeta, C_INT)

  IF (.NOT. gr_full_ready) RETURN
  IF (.NOT. ALLOCATED(gr_full_mode)) RETURN
  IF (.NOT. ALLOCATED(gr_full_B)) RETURN
  IF (.NOT. ALLOCATED(gr_full_logdet_w)) RETURN

  pos = 0
  DO i = 1, gr_full_nsub
    DO j = 1, gr_full_dim_z
      pos = pos + 1
      mode_out(pos) = gr_full_mode(j,i)
    END DO
  END DO

  pos = 0
  DO i = 1, gr_full_nsub
    DO k = 1, gr_full_dim_z
      DO j = 1, gr_full_dim_z
        pos = pos + 1
        B_out(pos) = gr_full_B(j,k,i)
      END DO
    END DO
  END DO

  DO i = 1, gr_full_nsub
    logdet_w_out(i) = gr_full_logdet_w(i)
  END DO

  ready0 = 1_C_INT

END SUBROUTINE get_growth_re_agh_full_equiv_grid_copy


END MODULE growth_re_marginal_agh_extension













MODULE growth_re_marginal_agh_score_extension
  USE shared
  USE utilities, ONLY: alnorm_fast
  USE core, ONLY: cleanup_growth_tt_allocations,                         &
                  setup_load_indices_from_loadmask_growth,               &
                  decode_lambda_sum_constraints_growth,                  &
                  refresh_lambda_act_growth,                             &
                  decode_gamma_sum_constraints_growth,                   &
                  decode_growth_thresholds, decode_growth_beta,           &
                  decode_growth_delta, decode_growth_random_structure,    &
                  growth_n_lambda_par, growth_n_gamma_par,                &
                  growth_n_beta_par, growth_n_zeta_par,                   &
                  build_growth_eta_moments_design,                       &
                  inneragh_item_logint_cond_eta,                          &
                  inneragh_item_logint_full_equiv,                        &
                  blockroot_tt_from_precision_long,                      &
                  activate_growth_subject, cfa_eta_mean,                  &
                  ensure_cfa_beta_storage,                                  &
                  tau_at, gamma_at, lambda_at, theta_at
  USE datindiv, ONLY: nmescur
  IMPLICIT NONE

   PRIVATE
   PUBLIC :: loglik_growth_re_tt_agh_full_equiv_score_q1
   PUBLIC :: set_growth_re_agh_full_equiv_score_grid
   PUBLIC :: loglik_growth_tt_agh_full_equiv_score_q1
   PUBLIC :: loglik_long_tt_agh_full_equiv_score_q1

  INTEGER, SAVE :: sc_npm = 0
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: sc_dtau(:,:), sc_dlambda(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: sc_dgamma(:,:,:), sc_dtheta(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: sc_dkappa(:,:,:), sc_dmap(:,:,:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: sc_dzeta_sd(:,:), sc_dsqrt_sigmak(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: sc_dL(:,:,:,:)

  INTEGER, SAVE :: gr_dim_u = 0
  INTEGER, SAVE :: gr_dim_eta = 0
  INTEGER, SAVE :: gr_dim_w = 0
  INTEGER, SAVE :: gr_nquad_u = 0
  INTEGER, SAVE :: gr_nquad_zeta = 0

  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_map_u(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_zeta_sd(:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_nodes_u(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_weights_u(:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_nodes_zeta(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_weights_zeta(:)

  LOGICAL, SAVE :: gr_nested_ready = .FALSE.
  INTEGER, SAVE :: gr_nested_fit_id = -999999999
  INTEGER, SAVE :: gr_nested_dim_u = -1
  INTEGER, SAVE :: gr_nested_nsub = -1
  INTEGER, SAVE :: gr_nested_nquad_u = -1
  INTEGER, SAVE :: gr_nested_nquad_zeta = -1
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_nested_mode(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_nested_root(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_nested_logdet(:)

  LOGICAL, SAVE :: gr_full_ready = .FALSE.
  INTEGER, SAVE :: gr_full_fit_id = -999999999
  INTEGER, SAVE :: gr_full_dim_w = -1
  INTEGER, SAVE :: gr_full_dim_z = -1
  INTEGER, SAVE :: gr_full_nsub = -1
  INTEGER, SAVE :: gr_full_nquad_u = -1
  INTEGER, SAVE :: gr_full_nquad_zeta = -1
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_full_mode(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_full_B(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: gr_full_logdet_w(:)

  INTERFACE
    SUBROUTINE DPOTRF(UPLO, N, A, LDA, INFO)
      CHARACTER(LEN=*), INTENT(IN) :: UPLO
      INTEGER, INTENT(IN) :: N, LDA
      DOUBLE PRECISION, DIMENSION(LDA,*), INTENT(INOUT) :: A
      INTEGER, INTENT(OUT) :: INFO
    END SUBROUTINE DPOTRF

    SUBROUTINE DPOTRS(UPLO, N, NRHS, A, LDA, B, LDB, INFO)
      CHARACTER(LEN=*), INTENT(IN) :: UPLO
      INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
      DOUBLE PRECISION, DIMENSION(LDA,*), INTENT(IN) :: A
      DOUBLE PRECISION, DIMENSION(LDB,*), INTENT(INOUT) :: B
      INTEGER, INTENT(OUT) :: INFO
    END SUBROUTINE DPOTRS
  END INTERFACE

CONTAINS

SUBROUTINE clear_growth_re_work()
  IMPLICIT NONE

  IF (ALLOCATED(gr_map_u)) DEALLOCATE(gr_map_u)
  IF (ALLOCATED(gr_zeta_sd)) DEALLOCATE(gr_zeta_sd)
  IF (ALLOCATED(gr_nodes_u)) DEALLOCATE(gr_nodes_u)
  IF (ALLOCATED(gr_weights_u)) DEALLOCATE(gr_weights_u)
  IF (ALLOCATED(gr_nodes_zeta)) DEALLOCATE(gr_nodes_zeta)
  IF (ALLOCATED(gr_weights_zeta)) DEALLOCATE(gr_weights_zeta)

  gr_dim_u = 0
  gr_dim_eta = 0
  gr_dim_w = 0
  gr_nquad_u = 0
  gr_nquad_zeta = 0
END SUBROUTINE clear_growth_re_work


SUBROUTINE clear_growth_re_model()
  IMPLICIT NONE

  IF (ALLOCATED(lambda_long)) DEALLOCATE(lambda_long)
  IF (ALLOCATED(gamma_long)) DEALLOCATE(gamma_long)
  IF (ALLOCATED(tau_long)) DEALLOCATE(tau_long)

  CALL cleanup_growth_tt_allocations()
  CALL clear_growth_re_work()
END SUBROUTINE clear_growth_re_model


SUBROUTINE invalidate_nested_cache(fit_id0)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: fit_id0

  IF (fit_id0 /= gr_nested_fit_id) THEN
    gr_nested_ready = .FALSE.
    gr_nested_fit_id = fit_id0
    IF (ALLOCATED(gr_nested_mode)) DEALLOCATE(gr_nested_mode)
    IF (ALLOCATED(gr_nested_root)) DEALLOCATE(gr_nested_root)
    IF (ALLOCATED(gr_nested_logdet)) DEALLOCATE(gr_nested_logdet)
    gr_nested_dim_u = -1
    gr_nested_nsub = -1
    gr_nested_nquad_u = -1
    gr_nested_nquad_zeta = -1
  END IF
END SUBROUTINE invalidate_nested_cache


SUBROUTINE invalidate_full_cache(fit_id0)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: fit_id0

  IF (fit_id0 /= gr_full_fit_id) THEN
    gr_full_ready = .FALSE.
    gr_full_fit_id = fit_id0
    IF (ALLOCATED(gr_full_mode)) DEALLOCATE(gr_full_mode)
    IF (ALLOCATED(gr_full_B)) DEALLOCATE(gr_full_B)
    IF (ALLOCATED(gr_full_logdet_w)) DEALLOCATE(gr_full_logdet_w)
    gr_full_dim_w = -1
    gr_full_dim_z = -1
    gr_full_nsub = -1
    gr_full_nquad_u = -1
    gr_full_nquad_zeta = -1
  END IF
END SUBROUTINE invalidate_full_cache


SUBROUTINE prepare_growth_re_model(                                      &
    b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,           &
    nquad_u0, nodes_u0, weights_u0,                                     &
    nquad_zeta0, nodes_zeta0, weights_zeta0,                            &
    nquad_re0, nodes_re0, weights_re0,                                  &
    nitem0, npm0, re_item0, theta_time0, delta_free0,                   &
    monotone_tau0, Q0, loadmask0,                                     &
    p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0, ok)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0, ntre0
  INTEGER, INTENT(IN) :: nquad_u0, nquad_zeta0, nquad_re0
  INTEGER, INTENT(IN) :: npm0, re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0, Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0, delta_x_col0, delta_z_col0
  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0), INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0), INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_u0
  DOUBLE PRECISION, DIMENSION(nquad_u0), INTENT(IN) :: weights_u0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_zeta0), INTENT(IN) :: weights_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: weights_re0
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: k, t, r, row, col, pos, jj, pos_tau, info
  INTEGER :: i, qidx, s, dimEta, dimU, nX, nZ, rr, zcol
  INTEGER :: ntheta_est, expected_npm, ntre_calc
  INTEGER :: n_lambda_free, n_gamma_free, n_Lu_est
  INTEGER :: nf, h, base
  DOUBLE PRECISION :: tmp_tau, zval
  INTEGER, ALLOCATABLE :: cm3(:)
  DOUBLE PRECISION, ALLOCATABLE :: delta_g(:)
  DOUBLE PRECISION, ALLOCATABLE :: Lu(:,:), Sigma_u(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: zeta_var(:,:)

  ok = .FALSE.

  IF (n0 < 2 .OR. Q0 < 1 .OR. ABS(p_fix0) < 1 .OR. p_ran0 < 1) RETURN
  IF (nquad_u0 < 1 .OR. nquad_zeta0 < 1 .OR. nquad_re0 < 1) RETURN

  CALL clear_growth_re_model()

  nvar = nvar0
  nsub = nsub0
  nobs = nobs0
  n = n0
  Q = Q0
  npm = npm0

  dimEta = n * Q
  dimU = Q * p_ran0
  dim_eta = dimEta

  gr_dim_u = dimU
  gr_dim_eta = dimEta
  gr_dim_w = dimU + dimEta
  gr_nquad_u = nquad_u0
  gr_nquad_zeta = nquad_zeta0

  nquad_re = nquad_re0
  re_item_flag = re_item0
  theta_time_flag = theta_time0
  phi_struct_flag = 0
  growth_subject_specific_design = MERGE(1, 0, p_fix0 > 0)
  growth_design_nsub = MERGE(nsub0, 1, growth_subject_specific_design == 1)
  p_fix = ABS(p_fix0)
  p_ran = p_ran0
  growth_delta_x_col = delta_x_col0
  growth_delta_z_col = delta_z_col0

  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs))
  ALLOCATE(c(nvar), cm1(nvar))
  ALLOCATE(nodes_re(nquad_re), weights_re(nquad_re))
  ALLOCATE(logweights_re(nquad_re))
  ALLOCATE(loadmask(nvar,Q))
  ALLOCATE(lambda(nvar,Q), gamma(nvar))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(theta_long(nvar,n), invsd_long(nvar,n))
  ALLOCATE(sigmak(nvar), sqrt_sigmak(nvar))
  ALLOCATE(kappa(dimEta), L_phi(dimEta,dimEta))
  ALLOCATE(kk_base(nvar))
  ALLOCATE(tau(SUM(cm10)))
  ALLOCATE(cm3(nvar))
  ALLOCATE(delta_g(Q))
  ALLOCATE(beta(Q,p_fix), X_eta(growth_design_nsub,n,Q,p_fix), &
           Z_eta(growth_design_nsub,n,Q,p_ran))
  ALLOCATE(growth_delta_mask(growth_design_nsub,n,Q))
  ALLOCATE(Lu(dimU,dimU), Sigma_u(dimU,dimU))
  ALLOCATE(zeta_var(n,Q))

  ALLOCATE(gr_map_u(dimEta,dimU,growth_design_nsub), gr_zeta_sd(dimEta))
  ALLOCATE(gr_nodes_u(dimU,nquad_u0), gr_weights_u(nquad_u0))
  ALLOCATE(gr_nodes_zeta(dimEta,nquad_zeta0))
  ALLOCATE(gr_weights_zeta(nquad_zeta0))

  nmes = nmes0
  Y = Y0
  nitem = nitem0
  c = c0
  cm1 = cm10
  loadmask = loadmask0

  nX = growth_design_nsub * n * Q * p_fix
  nZ = growth_design_nsub * n * Q * p_ran
  X_eta = RESHAPE(X0(1:nX), (/ growth_design_nsub, n, Q, p_fix /))
  Z_eta = RESHAPE(Z0(1:nZ), (/ growth_design_nsub, n, Q, p_ran /))
  growth_delta_mask = RESHAPE(delta_mask0(1:growth_design_nsub*n*Q), &
                              (/ growth_design_nsub, n, Q /))

  DO r = 1, nquad_u0
    base = (r - 1) * dimU
    gr_nodes_u(:,r) = nodes_u0(base+1:base+dimU)
  END DO
  gr_weights_u = weights_u0

  DO r = 1, nquad_zeta0
    base = (r - 1) * dimEta
    gr_nodes_zeta(:,r) = nodes_zeta0(base+1:base+dimEta)
  END DO
  gr_weights_zeta = weights_zeta0

  nodes_re = nodes_re0
  weights_re = weights_re0
  DO r = 1, nquad_re
    IF (weights_re(r) > 0.d0) THEN
      logweights_re(r) = LOG(weights_re(r))
    ELSE
      logweights_re(r) = -HUGE(1.d0)
    END IF
  END DO

  CALL setup_load_indices_from_loadmask_growth(info)
  IF (info /= 0) GOTO 900

  DO k = 1, nvar
    cm3(k) = MAX(cm1(k) - 2, 0)
  END DO
  ntre_calc = SUM(cm3)
  IF (ntre_calc /= ntre0) GOTO 900
  CALL decode_growth_thresholds(b0, npm0, ntre0, monotone_tau0, info)
  IF (info /= 0) GOTO 900

  pos = ntre0
  CALL decode_lambda_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 900
  CALL refresh_lambda_act_growth()

  CALL decode_gamma_sum_constraints_growth(b0, npm0, pos, info)
  IF (info /= 0) GOTO 900

  theta_k = 1.d0
  invsd_k = 1.d0
  theta_long = 1.d0
  invsd_long = 1.d0

  IF (theta_time_flag == 0) THEN
    DO k = 1, nvar
      IF (cm1(k) > 1) THEN
        pos = pos + 1
        theta_k(k) = EXP(b0(pos))
        invsd_k(k) = 1.d0 / SQRT(theta_k(k))
        DO t = 1, n
          theta_long(k,t) = theta_k(k)
          invsd_long(k,t) = invsd_k(k)
        END DO
      END IF
    END DO
  ELSE
    DO k = 1, nvar
      IF (cm1(k) > 1) THEN
        DO t = 1, n
          pos = pos + 1
          theta_long(k,t) = EXP(b0(pos))
          invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
        END DO
        theta_k(k) = theta_long(k,1)
        invsd_k(k) = invsd_long(k,1)
      END IF
    END DO
  END IF
  CALL decode_growth_beta(b0, npm0, pos, info)
  IF (info /= 0) GOTO 900
  CALL decode_growth_delta(b0, npm0, pos, delta_free0, delta_g, info)
  IF (info /= 0) GOTO 900
  CALL decode_growth_random_structure(b0, npm0, pos, dimU, Lu, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 900

  IF (re_item_flag == 1) THEN
    DO k = 1, nvar
      pos = pos + 1
      sigmak(k) = EXP(b0(pos))
      sqrt_sigmak(k) = SQRT(sigmak(k))
    END DO
  ELSE
    sigmak = 0.d0
    sqrt_sigmak = 0.d0
  END IF

  ntheta_est = 0
  IF (theta_time_flag == 0) THEN
    DO k = 1, nvar
      IF (cm1(k) > 1) ntheta_est = ntheta_est + 1
    END DO
  ELSE
    DO k = 1, nvar
      IF (cm1(k) > 1) ntheta_est = ntheta_est + n
    END DO
  END IF
  n_lambda_free = growth_n_lambda_par()
  n_gamma_free = growth_n_gamma_par()
  n_Lu_est = dimU * (dimU + 1) / 2

  expected_npm = ntre0                                      &
               + n_lambda_free                              &
               + n_gamma_free                               &
               + ntheta_est                                 &
               + growth_n_beta_par()                        &
               + MERGE(Q, 0, delta_free0 == 1)              &
               + n_Lu_est                                   &
               + growth_n_zeta_par()                        &
               + MERGE(nvar, 0, re_item_flag == 1)

  IF (npm0 /= expected_npm) GOTO 900
  IF (pos /= npm0) GOTO 900

  CALL build_growth_eta_moments_design(delta_g, Sigma_u, zeta_var, info)
  IF (info /= 0) GOTO 900

  gr_map_u = 0.d0
  DO i = 1, growth_design_nsub
    DO t = 1, n
      DO qidx = 1, Q
        row = (t - 1) * Q + qidx
        DO rr = 1, p_ran
          zval = Z_eta(i,t,qidx,rr)
          IF (rr == growth_delta_z_col .AND. growth_delta_mask(i,t,qidx) /= 0.d0) THEN
            IF (growth_delta_mode == 1) THEN
              zval = delta_g(qidx)
            ELSE
              zval = zval + delta_g(qidx) * growth_delta_mask(i,t,qidx)
            END IF
          END IF
          zcol = (rr - 1) * Q + qidx
          DO s = 1, dimU
            gr_map_u(row,s,i) = gr_map_u(row,s,i) + zval * Lu(zcol,s)
          END DO
        END DO
        gr_zeta_sd(row) = SQRT(zeta_var(t,qidx))
      END DO
    END DO
  END DO

  ok = .TRUE.

900 CONTINUE
  IF (ALLOCATED(cm3)) DEALLOCATE(cm3)
  IF (ALLOCATED(delta_g)) DEALLOCATE(delta_g)
  IF (ALLOCATED(Lu)) DEALLOCATE(Lu)
  IF (ALLOCATED(Sigma_u)) DEALLOCATE(Sigma_u)
  IF (ALLOCATED(zeta_var)) DEALLOCATE(zeta_var)

  IF (.NOT. ok) CALL clear_growth_re_model()
END SUBROUTINE prepare_growth_re_model


SUBROUTINE eta_from_u_zeta(i, z_u, z_zeta, eta_vec)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  INTEGER :: design_i
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(IN) :: z_u
  DOUBLE PRECISION, DIMENSION(gr_dim_eta), INTENT(IN) :: z_zeta
  DOUBLE PRECISION, DIMENSION(gr_dim_eta), INTENT(OUT) :: eta_vec

  design_i = MERGE(i, 1, growth_subject_specific_design == 1)
  eta_vec = kappa_sub(:,design_i) + MATMUL(gr_map_u(:,:,design_i), z_u) + &
            gr_zeta_sd * z_zeta
END SUBROUTINE eta_from_u_zeta


SUBROUTINE eta_from_w(i, z_w, eta_vec)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  INTEGER :: design_i
  DOUBLE PRECISION, DIMENSION(gr_dim_w), INTENT(IN) :: z_w
  DOUBLE PRECISION, DIMENSION(gr_dim_eta), INTENT(OUT) :: eta_vec

  design_i = MERGE(i, 1, growth_subject_specific_design == 1)
  eta_vec = kappa_sub(:,design_i) + &
            MATMUL(gr_map_u(:,:,design_i), z_w(1:gr_dim_u)) + &
            gr_zeta_sd * z_w(gr_dim_u+1:gr_dim_w)
END SUBROUTINE eta_from_w


SUBROUTINE log_items_integrated_a(i, eta_vec, logg, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(gr_dim_eta), INTENT(IN) :: eta_vec
  DOUBLE PRECISION, INTENT(OUT) :: logg
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: k_item
  DOUBLE PRECISION :: logint
  LOGICAL :: ok_item

  ok = .TRUE.
  logg = 0.d0

  DO k_item = 1, nvar
    CALL inneragh_item_logint_cond_eta(i, k_item, eta_vec, logint, ok_item)
    IF (.NOT. ok_item .OR. logint /= logint) THEN
      ok = .FALSE.
      logg = -1.d9
      RETURN
    END IF
    logg = logg + logint
  END DO
END SUBROUTINE log_items_integrated_a


SUBROUTINE log_marginal_zeta_given_u(i, z_u, logg, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(IN) :: z_u
  DOUBLE PRECISION, INTENT(OUT) :: logg
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: lz, nvalid
  DOUBLE PRECISION :: logitem, term, maxterm, sumexp
  DOUBLE PRECISION, ALLOCATABLE :: eta_vec(:)
  LOGICAL :: ok_item

  ok = .TRUE.
  logg = -1.d9
  nvalid = 0
  maxterm = 0.d0
  sumexp = 0.d0

  ALLOCATE(eta_vec(gr_dim_eta))

  DO lz = 1, gr_nquad_zeta
    IF (gr_weights_zeta(lz) <= 0.d0) CYCLE

    CALL eta_from_u_zeta(i, z_u, gr_nodes_zeta(:,lz), eta_vec)
    CALL log_items_integrated_a(i, eta_vec, logitem, ok_item)
    IF (.NOT. ok_item) CYCLE

    term = LOG(gr_weights_zeta(lz)) + logitem
    IF (term /= term) CYCLE

    IF (nvalid == 0) THEN
      maxterm = term
      sumexp = 1.d0
      nvalid = 1
    ELSE IF (term > maxterm) THEN
      sumexp = sumexp * EXP(maxterm - term) + 1.d0
      maxterm = term
      nvalid = nvalid + 1
    ELSE
      sumexp = sumexp + EXP(term - maxterm)
      nvalid = nvalid + 1
    END IF
  END DO

  IF (nvalid == 0 .OR. sumexp <= 0.d0 .OR. sumexp /= sumexp) THEN
    ok = .FALSE.
  ELSE
    logg = maxterm + LOG(sumexp)
  END IF

  DEALLOCATE(eta_vec)
END SUBROUTINE log_marginal_zeta_given_u


SUBROUTINE nested_outer_obj(i, z_u, f, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(IN) :: z_u
  DOUBLE PRECISION, INTENT(OUT) :: f
  LOGICAL, INTENT(OUT) :: ok
  DOUBLE PRECISION :: logg

  CALL log_marginal_zeta_given_u(i, z_u, logg, ok)
  IF (.NOT. ok) THEN
    f = 1.d9
  ELSE
    f = -logg + 0.5d0 * DOT_PRODUCT(z_u,z_u)
  END IF
END SUBROUTINE nested_outer_obj


SUBROUTINE nested_fd_grad_hess(i, z, f0, grad, H, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(IN) :: z
  DOUBLE PRECISION, INTENT(OUT) :: f0
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(OUT) :: grad
  DOUBLE PRECISION, DIMENSION(gr_dim_u,gr_dim_u), INTENT(OUT) :: H
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: r, s
  DOUBLE PRECISION :: hr, hs, fp, fm, fpp, fpm, fmp, fmm
  DOUBLE PRECISION, PARAMETER :: epsg = 2.d-4
  DOUBLE PRECISION, PARAMETER :: epsh = 8.d-4
  DOUBLE PRECISION, ALLOCATABLE :: zp(:), zm(:), zpp(:), zpm(:)
  DOUBLE PRECISION, ALLOCATABLE :: zmp(:), zmm(:)
  LOGICAL :: ok0, okp, okm, okpp, okpm, okmp, okmm

  ALLOCATE(zp(gr_dim_u), zm(gr_dim_u))
  ALLOCATE(zpp(gr_dim_u), zpm(gr_dim_u))
  ALLOCATE(zmp(gr_dim_u), zmm(gr_dim_u))

  CALL nested_outer_obj(i, z, f0, ok0)
  IF (.NOT. ok0 .OR. f0 /= f0) THEN
    ok = .FALSE.
    GOTO 800
  END IF

  grad = 0.d0
  H = 0.d0

  DO r = 1, gr_dim_u
    hr = epsg * MAX(ABS(z(r)), 1.d0)
    zp = z
    zm = z
    zp(r) = zp(r) + hr
    zm(r) = zm(r) - hr
    CALL nested_outer_obj(i, zp, fp, okp)
    CALL nested_outer_obj(i, zm, fm, okm)
    IF (.NOT. okp .OR. .NOT. okm) THEN
      ok = .FALSE.
      GOTO 800
    END IF
    grad(r) = (fp - fm) / (2.d0 * hr)

    hr = epsh * MAX(ABS(z(r)), 1.d0)
    zp = z
    zm = z
    zp(r) = zp(r) + hr
    zm(r) = zm(r) - hr
    CALL nested_outer_obj(i, zp, fp, okp)
    CALL nested_outer_obj(i, zm, fm, okm)
    IF (.NOT. okp .OR. .NOT. okm) THEN
      ok = .FALSE.
      GOTO 800
    END IF
    H(r,r) = (fp - 2.d0*f0 + fm) / (hr*hr)
  END DO

  DO r = 1, gr_dim_u - 1
    hr = epsh * MAX(ABS(z(r)), 1.d0)
    DO s = r + 1, gr_dim_u
      hs = epsh * MAX(ABS(z(s)), 1.d0)
      zpp = z
      zpm = z
      zmp = z
      zmm = z
      zpp(r) = zpp(r) + hr
      zpp(s) = zpp(s) + hs
      zpm(r) = zpm(r) + hr
      zpm(s) = zpm(s) - hs
      zmp(r) = zmp(r) - hr
      zmp(s) = zmp(s) + hs
      zmm(r) = zmm(r) - hr
      zmm(s) = zmm(s) - hs
      CALL nested_outer_obj(i, zpp, fpp, okpp)
      CALL nested_outer_obj(i, zpm, fpm, okpm)
      CALL nested_outer_obj(i, zmp, fmp, okmp)
      CALL nested_outer_obj(i, zmm, fmm, okmm)
      IF (.NOT. okpp .OR. .NOT. okpm .OR. .NOT. okmp .OR. .NOT. okmm) THEN
        ok = .FALSE.
        GOTO 800
      END IF
      H(r,s) = (fpp - fpm - fmp + fmm) / (4.d0*hr*hs)
      H(s,r) = H(r,s)
    END DO
  END DO

  ok = .TRUE.

800 CONTINUE
  DEALLOCATE(zp, zm, zpp, zpm, zmp, zmm)
END SUBROUTINE nested_fd_grad_hess


SUBROUTINE solve_spd_step(Hin, grad, step, ok)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Hin
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: grad
  DOUBLE PRECISION, DIMENSION(SIZE(grad)), INTENT(OUT) :: step
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: d, info
  DOUBLE PRECISION, ALLOCATABLE :: A(:,:), rhs(:,:)

  d = SIZE(grad)
  ALLOCATE(A(d,d), rhs(d,1))
  ok = .FALSE.

  A = 0.5d0 * (Hin + TRANSPOSE(Hin))
  CALL DPOTRF('U', d, A, d, info)
  IF (info == 0) THEN
    rhs(:,1) = -grad
    CALL DPOTRS('U', d, 1, A, d, rhs, d, info)
    IF (info == 0 .AND. ALL(rhs(:,1) == rhs(:,1))) THEN
      step = rhs(:,1)
      ok = .TRUE.
    END IF
  END IF

  DEALLOCATE(A, rhs)
END SUBROUTINE solve_spd_step


SUBROUTINE invert_upper(U, Uinv, ok)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: U
  DOUBLE PRECISION, DIMENSION(SIZE(U,1),SIZE(U,2)), INTENT(OUT) :: Uinv
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: d, i, j, k
  DOUBLE PRECISION :: rhs

  d = SIZE(U,1)
  Uinv = 0.d0
  ok = .FALSE.
  IF (SIZE(U,2) /= d) RETURN

  DO j = 1, d
    DO i = d, 1, -1
      IF (ABS(U(i,i)) <= TINY(1.d0) .OR. U(i,i) /= U(i,i)) RETURN
      rhs = 0.d0
      IF (i == j) rhs = 1.d0
      DO k = i + 1, d
        rhs = rhs - U(i,k) * Uinv(k,j)
      END DO
      Uinv(i,j) = rhs / U(i,i)
    END DO
  END DO

  ok = .TRUE.
END SUBROUTINE invert_upper


SUBROUTINE nested_find_mode_root(i, start, mode, R, logdetR, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(IN) :: start
  DOUBLE PRECISION, DIMENSION(gr_dim_u), INTENT(OUT) :: mode
  DOUBLE PRECISION, DIMENSION(gr_dim_u,gr_dim_u), INTENT(OUT) :: R
  DOUBLE PRECISION, INTENT(OUT) :: logdetR
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: iter, ls, info, j
  DOUBLE PRECISION :: f, fnew, alpha, gnorm
  DOUBLE PRECISION, ALLOCATABLE :: grad(:), step(:), znew(:)
  DOUBLE PRECISION, ALLOCATABLE :: H(:,:), A(:,:)
  LOGICAL :: okgh, okstep, oknew

  ALLOCATE(grad(gr_dim_u), step(gr_dim_u), znew(gr_dim_u))
  ALLOCATE(H(gr_dim_u,gr_dim_u), A(gr_dim_u,gr_dim_u))

  mode = start
  ok = .FALSE.

  DO iter = 1, 35
    CALL nested_fd_grad_hess(i, mode, f, grad, H, okgh)
    IF (.NOT. okgh) GOTO 900

    gnorm = SQRT(DOT_PRODUCT(grad,grad))
    IF (gnorm < 2.d-5) EXIT

    CALL solve_spd_step(H, grad, step, okstep)
    IF (.NOT. okstep) step = -grad / MAX(gnorm,1.d0)

    alpha = 1.d0
    oknew = .FALSE.
    DO ls = 1, 14
      znew = mode + alpha * step
      CALL nested_outer_obj(i, znew, fnew, oknew)
      IF (oknew .AND. fnew < f) EXIT
      alpha = alpha * 0.5d0
    END DO

    IF (.NOT. oknew .OR. fnew >= f) THEN
      step = -grad / MAX(gnorm,1.d0)
      alpha = 0.25d0
      DO ls = 1, 14
        znew = mode + alpha * step
        CALL nested_outer_obj(i, znew, fnew, oknew)
        IF (oknew .AND. fnew < f) EXIT
        alpha = alpha * 0.5d0
      END DO
    END IF

    IF (.NOT. oknew .OR. fnew >= f) EXIT
    mode = znew
    IF (SQRT(DOT_PRODUCT(alpha*step,alpha*step)) < 1.d-6) EXIT
  END DO

  CALL nested_fd_grad_hess(i, mode, f, grad, H, okgh)
  IF (.NOT. okgh) GOTO 900

  A = 0.5d0 * (H + TRANSPOSE(H))
  CALL DPOTRF('U', gr_dim_u, A, gr_dim_u, info)
  IF (info /= 0) GOTO 900

  logdetR = 0.d0
  DO j = 1, gr_dim_u
    IF (A(j,j) <= 0.d0) GOTO 900
    logdetR = logdetR - LOG(A(j,j))
  END DO

  IF (gr_dim_u > 1) THEN
    DO j = 1, gr_dim_u-1
      A(j+1:gr_dim_u,j) = 0.d0
    END DO
  END IF
  CALL invert_upper(A, R, okgh)
  IF (.NOT. okgh) GOTO 900

  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(grad, step, znew, H, A)
END SUBROUTINE nested_find_mode_root


SUBROUTINE build_nested_cache(ok)
  IMPLICIT NONE
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: i, sm
  DOUBLE PRECISION, ALLOCATABLE :: start(:), mode(:), R(:,:)
  DOUBLE PRECISION :: logdetR
  LOGICAL :: ok_i

  ok = .FALSE.

  IF (ALLOCATED(gr_nested_mode)) THEN
    IF (gr_nested_dim_u /= gr_dim_u .OR. gr_nested_nsub /= nsub .OR. &
        gr_nested_nquad_u /= gr_nquad_u .OR.                      &
        gr_nested_nquad_zeta /= gr_nquad_zeta) THEN
      DEALLOCATE(gr_nested_mode, gr_nested_root, gr_nested_logdet)
      gr_nested_ready = .FALSE.
    END IF
  END IF

  IF (.NOT. ALLOCATED(gr_nested_mode)) THEN
    ALLOCATE(gr_nested_mode(gr_dim_u,nsub))
    ALLOCATE(gr_nested_root(gr_dim_u,gr_dim_u,nsub))
    ALLOCATE(gr_nested_logdet(nsub))
    gr_nested_mode = 0.d0
    gr_nested_root = 0.d0
    gr_nested_logdet = 0.d0
    gr_nested_dim_u = gr_dim_u
    gr_nested_nsub = nsub
    gr_nested_nquad_u = gr_nquad_u
    gr_nested_nquad_zeta = gr_nquad_zeta
  END IF

  ALLOCATE(start(gr_dim_u), mode(gr_dim_u), R(gr_dim_u,gr_dim_u))
  nmescur = 0

  DO i = 1, nsub
    sm = SUM(nmes(i,:))
    IF (gr_nested_ready) THEN
      start = gr_nested_mode(:,i)
    ELSE
      start = 0.d0
    END IF

    CALL nested_find_mode_root(i, start, mode, R, logdetR, ok_i)
    IF (.NOT. ok_i) THEN
      gr_nested_ready = .FALSE.
      GOTO 900
    END IF

    gr_nested_mode(:,i) = mode
    gr_nested_root(:,:,i) = R
    gr_nested_logdet(i) = logdetR
    nmescur = nmescur + sm
  END DO

  gr_nested_ready = .TRUE.
  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(start, mode, R)
END SUBROUTINE build_nested_cache


SUBROUTINE eval_nested_loglik(out, ok)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(OUT) :: out
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: i, lu, sm, nvalid
  DOUBLE PRECISION :: logg, term, maxterm, sumexp, qi, qv
  DOUBLE PRECISION, ALLOCATABLE :: z_u(:), v_u(:)
  LOGICAL :: okg

  out = 0.d0
  ok = .FALSE.
  IF (.NOT. gr_nested_ready) RETURN

  ALLOCATE(z_u(gr_dim_u), v_u(gr_dim_u))
  nmescur = 0

  DO i = 1, nsub
    sm = SUM(nmes(i,:))
    nvalid = 0
    maxterm = 0.d0
    sumexp = 0.d0

    DO lu = 1, gr_nquad_u
      IF (gr_weights_u(lu) <= 0.d0) CYCLE
      v_u = gr_nodes_u(:,lu)
      z_u = gr_nested_mode(:,i) + MATMUL(gr_nested_root(:,:,i),v_u)

      CALL log_marginal_zeta_given_u(i, z_u, logg, okg)
      IF (.NOT. okg) CYCLE

      qi = DOT_PRODUCT(z_u,z_u)
      qv = DOT_PRODUCT(v_u,v_u)
      term = LOG(gr_weights_u(lu)) + logg - 0.5d0*(qi-qv)
      IF (term /= term) CYCLE

      IF (nvalid == 0) THEN
        maxterm = term
        sumexp = 1.d0
        nvalid = 1
      ELSE IF (term > maxterm) THEN
        sumexp = sumexp * EXP(maxterm-term) + 1.d0
        maxterm = term
        nvalid = nvalid + 1
      ELSE
        sumexp = sumexp + EXP(term-maxterm)
        nvalid = nvalid + 1
      END IF
    END DO

    IF (nvalid == 0 .OR. sumexp <= 0.d0) GOTO 900
    out = out + gr_nested_logdet(i) + maxterm + LOG(sumexp)
    nmescur = nmescur + sm
  END DO

  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(z_u, v_u)
END SUBROUTINE eval_nested_loglik


SUBROUTINE full_grad_hess(i, z, f, grad, H, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: z
  DOUBLE PRECISION, INTENT(OUT) :: f
  DOUBLE PRECISION, DIMENSION(SIZE(z)), INTENT(OUT) :: grad
  DOUBLE PRECISION, DIMENSION(SIZE(z),SIZE(z)), INTENT(OUT) :: H
  LOGICAL, INTENT(OUT) :: ok

  DOUBLE PRECISION, PARAMETER :: INV_SQRT_2PI = &
       0.39894228040143267794d0
  DOUBLE PRECISION, PARAMETER :: P_FLOOR = 1.d-300

  INTEGER :: dimZ, k_item, j_mes, idx_obs, yijk, tuk
  INTEGER :: t_idx, kk_thr, qq, row, d1, d2, design_i
  DOUBLE PRECISION :: mu, linpred, var_e, invsd
  DOUBLE PRECISION :: aa, bb, p, phia, phib
  DOUBLE PRECISION :: dp, d2p, smu, hmu, loglik
  DOUBLE PRECISION, ALLOCATABLE :: eta_vec(:), a_item(:), dmu(:)
  DOUBLE PRECISION, ALLOCATABLE :: z_w(:)

  dimZ = SIZE(z)
  design_i = MERGE(i, 1, growth_subject_specific_design == 1)
  ok = .FALSE.
  f = 1.d9
  grad = 0.d0
  H = 0.d0

  IF (re_item_flag == 1) THEN
    IF (dimZ /= gr_dim_w + nvar) RETURN
  ELSE
    IF (dimZ /= gr_dim_w) RETURN
  END IF

  ALLOCATE(eta_vec(gr_dim_eta), a_item(nvar), dmu(dimZ))
  ALLOCATE(z_w(gr_dim_w))

  z_w = z(1:gr_dim_w)
  CALL eta_from_w(i, z_w, eta_vec)

  IF (re_item_flag == 1) THEN
    DO k_item = 1, nvar
      a_item(k_item) = sqrt_sigmak(k_item) * z(gr_dim_w+k_item)
    END DO
  ELSE
    a_item = 0.d0
  END IF

  loglik = 0.d0
  kk_thr = 0
  idx_obs = nmescur

  DO k_item = 1, nvar
    DO j_mes = 1, nmes(i,k_item)
      idx_obs = idx_obs + 1
      t_idx = nitem(idx_obs)
      var_e = theta_at(k_item,t_idx)
      IF (var_e <= 0.d0 .OR. var_e /= var_e) GOTO 900
      invsd = 1.d0 / SQRT(var_e)

      linpred = 0.d0
      DO qq = 1, Q
        IF (loadmask(k_item,qq) == 0) CYCLE
        row = (t_idx-1)*Q + qq
        linpred = linpred + lambda_at(k_item,qq,t_idx) * eta_vec(row)
      END DO
      mu = gamma_at(k_item,t_idx) + a_item(k_item) + linpred

      yijk = Y(idx_obs)
      tuk = yijk + kk_thr

      IF (c(k_item) == 2) THEN
        aa = (tau_at(kk_thr+1,t_idx)-mu) * invsd
        phia = INV_SQRT_2PI * EXP(-0.5d0*aa*aa)
        IF (yijk == 1) THEN
          p = alnorm_fast(aa,.FALSE.)
          dp = -phia*invsd
          d2p = -aa*phia*invsd*invsd
        ELSE
          p = 1.d0-alnorm_fast(aa,.FALSE.)
          dp = phia*invsd
          d2p = aa*phia*invsd*invsd
        END IF
      ELSE IF (yijk == 1) THEN
        aa = (tau_at(tuk,t_idx)-mu) * invsd
        p = alnorm_fast(aa,.FALSE.)
        phia = INV_SQRT_2PI * EXP(-0.5d0*aa*aa)
        dp = -phia*invsd
        d2p = -aa*phia*invsd*invsd
      ELSE IF (yijk == c(k_item)) THEN
        bb = (tau_at(tuk-1,t_idx)-mu) * invsd
        p = 1.d0-alnorm_fast(bb,.FALSE.)
        phib = INV_SQRT_2PI * EXP(-0.5d0*bb*bb)
        dp = phib*invsd
        d2p = bb*phib*invsd*invsd
      ELSE
        aa = (tau_at(tuk,t_idx)-mu) * invsd
        bb = (tau_at(tuk-1,t_idx)-mu) * invsd
        p = alnorm_fast(aa,.FALSE.)-alnorm_fast(bb,.FALSE.)
        phia = INV_SQRT_2PI * EXP(-0.5d0*aa*aa)
        phib = INV_SQRT_2PI * EXP(-0.5d0*bb*bb)
        dp = (phib-phia)*invsd
        d2p = (bb*phib-aa*phia)*invsd*invsd
      END IF

      IF (p <= 0.d0 .OR. p /= p) GOTO 900
      p = MAX(p,P_FLOOR)
      loglik = loglik + LOG(p)
      smu = dp/p
      hmu = d2p/p-smu*smu
      IF (smu /= smu .OR. hmu /= hmu) GOTO 900

      dmu = 0.d0
      DO qq = 1, Q
        IF (loadmask(k_item,qq) == 0) CYCLE
        row = (t_idx-1)*Q + qq
        DO d1 = 1, gr_dim_u
          dmu(d1) = dmu(d1) + &
               lambda_at(k_item,qq,t_idx)*gr_map_u(row,d1,design_i)
        END DO
        d1 = gr_dim_u + row
        dmu(d1) = dmu(d1) + &
             lambda_at(k_item,qq,t_idx)*gr_zeta_sd(row)
      END DO
      IF (re_item_flag == 1) THEN
        dmu(gr_dim_w+k_item) = sqrt_sigmak(k_item)
      END IF

      grad = grad + smu*dmu
      DO d1 = 1, dimZ
        IF (dmu(d1) == 0.d0) CYCLE
        DO d2 = 1, dimZ
          IF (dmu(d2) == 0.d0) CYCLE
          H(d1,d2) = H(d1,d2) + hmu*dmu(d1)*dmu(d2)
        END DO
      END DO
    END DO
    kk_thr = kk_thr + c(k_item)-1
  END DO

  f = -loglik + 0.5d0*DOT_PRODUCT(z,z)
  grad = -grad + z
  H = -H
  DO d1 = 1, dimZ
    H(d1,d1) = H(d1,d1) + 1.d0
  END DO
  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(eta_vec, a_item, dmu, z_w)
END SUBROUTINE full_grad_hess


SUBROUTINE full_find_mode_B(i, start, mode, B, logdetW, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: start
  DOUBLE PRECISION, DIMENSION(SIZE(start)), INTENT(OUT) :: mode
  DOUBLE PRECISION, DIMENSION(SIZE(start),SIZE(start)), INTENT(OUT) :: B
  DOUBLE PRECISION, INTENT(OUT) :: logdetW
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: dimZ, iter, ls
  DOUBLE PRECISION :: f, fnew, alpha, gnorm, logdetB
  DOUBLE PRECISION, ALLOCATABLE :: grad(:), step(:), znew(:)
  DOUBLE PRECISION, ALLOCATABLE :: H(:,:), Hreg(:,:)
  LOGICAL :: okgh, okstep, oknew, okroot

  dimZ = SIZE(start)
  ALLOCATE(grad(dimZ), step(dimZ), znew(dimZ))
  ALLOCATE(H(dimZ,dimZ), Hreg(dimZ,dimZ))

  mode = start
  ok = .FALSE.

  DO iter = 1, 50
    CALL full_grad_hess(i, mode, f, grad, H, okgh)
    IF (.NOT. okgh) GOTO 900
    gnorm = SQRT(DOT_PRODUCT(grad,grad))
    IF (gnorm < 1.d-6) EXIT

    CALL solve_spd_step(H, grad, step, okstep)
    IF (.NOT. okstep) step = -grad/MAX(gnorm,1.d0)

    alpha = 1.d0
    oknew = .FALSE.
    DO ls = 1, 16
      znew = mode + alpha*step
      CALL full_grad_hess(i, znew, fnew, grad, Hreg, oknew)
      IF (oknew .AND. fnew < f) EXIT
      alpha = alpha*0.5d0
    END DO

    IF (.NOT. oknew .OR. fnew >= f) THEN
      CALL full_grad_hess(i, mode, f, grad, H, okgh)
      step = -grad/MAX(SQRT(DOT_PRODUCT(grad,grad)),1.d0)
      alpha = 0.25d0
      DO ls = 1, 16
        znew = mode + alpha*step
        CALL full_grad_hess(i, znew, fnew, grad, Hreg, oknew)
        IF (oknew .AND. fnew < f) EXIT
        alpha = alpha*0.5d0
      END DO
    END IF

    IF (.NOT. oknew .OR. fnew >= f) EXIT
    mode = znew
    IF (SQRT(DOT_PRODUCT(alpha*step,alpha*step)) < 1.d-7) EXIT
  END DO

  CALL full_grad_hess(i, mode, f, grad, H, okgh)
  IF (.NOT. okgh) GOTO 900

  Hreg = 0.5d0*(H+TRANSPOSE(H))
  CALL blockroot_tt_from_precision_long(gr_dim_w, dimZ, Hreg, &
       B, logdetB, logdetW, okroot)
  IF (.NOT. okroot) GOTO 900

  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(grad, step, znew, H, Hreg)
END SUBROUTINE full_find_mode_B


SUBROUTINE build_full_cache(ok)
  IMPLICIT NONE
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: i, sm, dimZ
  DOUBLE PRECISION, ALLOCATABLE :: start(:), mode(:), B(:,:)
  DOUBLE PRECISION :: logdetW
  LOGICAL :: ok_i

  ok = .FALSE.
  dimZ = gr_dim_w + MERGE(nvar,0,re_item_flag == 1)

  IF (ALLOCATED(gr_full_mode)) THEN
    IF (gr_full_dim_w /= gr_dim_w .OR. gr_full_dim_z /= dimZ .OR. &
        gr_full_nsub /= nsub .OR.                                &
        gr_full_nquad_u /= gr_nquad_u .OR.                       &
        gr_full_nquad_zeta /= gr_nquad_zeta) THEN
      DEALLOCATE(gr_full_mode, gr_full_B, gr_full_logdet_w)
      gr_full_ready = .FALSE.
    END IF
  END IF

  IF (.NOT. ALLOCATED(gr_full_mode)) THEN
    ALLOCATE(gr_full_mode(dimZ,nsub))
    ALLOCATE(gr_full_B(dimZ,dimZ,nsub))
    ALLOCATE(gr_full_logdet_w(nsub))
    gr_full_mode = 0.d0
    gr_full_B = 0.d0
    gr_full_logdet_w = 0.d0
    gr_full_dim_w = gr_dim_w
    gr_full_dim_z = dimZ
    gr_full_nsub = nsub
    gr_full_nquad_u = gr_nquad_u
    gr_full_nquad_zeta = gr_nquad_zeta
  END IF

  ALLOCATE(start(dimZ), mode(dimZ), B(dimZ,dimZ))
  nmescur = 0

  DO i = 1, nsub
    sm = SUM(nmes(i,:))
    IF (gr_full_ready) THEN
      start = gr_full_mode(:,i)
    ELSE
      start = 0.d0
    END IF

    CALL full_find_mode_B(i, start, mode, B, logdetW, ok_i)
    IF (.NOT. ok_i) THEN
      gr_full_ready = .FALSE.
      GOTO 900
    END IF

    gr_full_mode(:,i) = mode
    gr_full_B(:,:,i) = B
    gr_full_logdet_w(i) = logdetW
    nmescur = nmescur + sm
  END DO

  gr_full_ready = .TRUE.
  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(start, mode, B)
END SUBROUTINE build_full_cache


SUBROUTINE eval_full_loglik(out, ok)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(OUT) :: out
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: i, lu, lz, k_item, sm, nvalid, dimZ
  DOUBLE PRECISION :: logg, logint, term, maxterm, sumexp
  DOUBLE PRECISION :: qw, qv, alpha_k, sdev_k
  DOUBLE PRECISION, ALLOCATABLE :: v_w(:), z_w(:), eta_vec(:)
  LOGICAL :: ok_item

  out = 0.d0
  ok = .FALSE.
  IF (.NOT. gr_full_ready) RETURN

  dimZ = gr_full_dim_z
  ALLOCATE(v_w(gr_dim_w), z_w(gr_dim_w), eta_vec(gr_dim_eta))
  nmescur = 0

  DO i = 1, nsub
    sm = SUM(nmes(i,:))
    nvalid = 0
    maxterm = 0.d0
    sumexp = 0.d0

    DO lu = 1, gr_nquad_u
      IF (gr_weights_u(lu) <= 0.d0) CYCLE
      DO lz = 1, gr_nquad_zeta
        IF (gr_weights_zeta(lz) <= 0.d0) CYCLE

        v_w(1:gr_dim_u) = gr_nodes_u(:,lu)
        v_w(gr_dim_u+1:gr_dim_w) = gr_nodes_zeta(:,lz)
        z_w = gr_full_mode(1:gr_dim_w,i) + &
              MATMUL(gr_full_B(1:gr_dim_w,1:gr_dim_w,i),v_w)
        CALL eta_from_w(i, z_w, eta_vec)

        logg = 0.d0
        DO k_item = 1, nvar
          IF (re_item_flag == 1) THEN
            alpha_k = gr_full_mode(gr_dim_w+k_item,i) + &
                 DOT_PRODUCT(gr_full_B(gr_dim_w+k_item,1:gr_dim_w,i),v_w)
            sdev_k = gr_full_B(gr_dim_w+k_item,gr_dim_w+k_item,i)
          ELSE
            alpha_k = 0.d0
            sdev_k = 1.d0
          END IF

          CALL inneragh_item_logint_full_equiv(i, k_item, eta_vec, &
               alpha_k, sdev_k, logint, ok_item)
          IF (.NOT. ok_item) THEN
            logg = -1.d9
            EXIT
          END IF
          logg = logg + logint
        END DO
        IF (logg <= -1.d9) CYCLE

        qw = DOT_PRODUCT(z_w,z_w)
        qv = DOT_PRODUCT(v_w,v_w)
        term = LOG(gr_weights_u(lu)) + LOG(gr_weights_zeta(lz)) + &
               logg - 0.5d0*(qw-qv)
        IF (term /= term) CYCLE

        IF (nvalid == 0) THEN
          maxterm = term
          sumexp = 1.d0
          nvalid = 1
        ELSE IF (term > maxterm) THEN
          sumexp = sumexp*EXP(maxterm-term)+1.d0
          maxterm = term
          nvalid = nvalid+1
        ELSE
          sumexp = sumexp+EXP(term-maxterm)
          nvalid = nvalid+1
        END IF
      END DO
    END DO

    IF (nvalid == 0 .OR. sumexp <= 0.d0) GOTO 900
    out = out + gr_full_logdet_w(i) + maxterm + LOG(sumexp)
    nmescur = nmescur + sm
  END DO

  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(v_w, z_w, eta_vec)
END SUBROUTINE eval_full_loglik


SUBROUTINE loglik_growth_re_tt_agh_nested(                                &
    b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,            &
    nquad_u0, nodes_u0, weights_u0,                                      &
    nquad_zeta0, nodes_zeta0, weights_zeta0,                             &
    nquad_re0, nodes_re0, weights_re0,                                   &
    nitem0, npm0, re_item0, theta_time0, delta_free0,                    &
    monotone_tau0, update_gh0, fit_id0, Q0, loadmask0,             &
    p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0, vrais)            &
    BIND(C, name="score_aux_loglik_growth_re_tt_agh_nested")

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0, ntre0
  INTEGER, INTENT(IN) :: nquad_u0, nquad_zeta0, nquad_re0
  INTEGER, INTENT(IN) :: npm0, re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0
  INTEGER, INTENT(IN) :: update_gh0, fit_id0, Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0, delta_x_col0, delta_z_col0
  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0), INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0), INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_u0
  DOUBLE PRECISION, DIMENSION(nquad_u0), INTENT(IN) :: weights_u0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_zeta0), INTENT(IN) :: weights_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: weights_re0
  DOUBLE PRECISION, INTENT(OUT) :: vrais

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.d9
  DOUBLE PRECISION :: ll
  LOGICAL :: ok_model, ok_grid, ok_ll

  vrais = BAD_OBJ
  CALL invalidate_nested_cache(fit_id0)

  CALL prepare_growth_re_model(                                       &
       b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,      &
       nquad_u0, nodes_u0, weights_u0,                                &
       nquad_zeta0, nodes_zeta0, weights_zeta0,                       &
       nquad_re0, nodes_re0, weights_re0,                             &
       nitem0, npm0, re_item0, theta_time0, delta_free0,              &
       monotone_tau0, Q0, loadmask0,                              &
       p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0, ok_model)
  IF (.NOT. ok_model) GOTO 900

  IF (update_gh0 == 1 .OR. .NOT. gr_nested_ready) THEN
    CALL build_nested_cache(ok_grid)
    IF (.NOT. ok_grid) GOTO 900
  END IF

  CALL eval_nested_loglik(ll, ok_ll)
  IF (ok_ll .AND. ll == ll .AND. ABS(ll) < 1.d300) vrais = -ll

900 CONTINUE
  CALL clear_growth_re_model()
END SUBROUTINE loglik_growth_re_tt_agh_nested


SUBROUTINE loglik_growth_re_tt_agh_full_equiv(                             &
    b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,             &
    nquad_u0, nodes_u0, weights_u0,                                       &
    nquad_zeta0, nodes_zeta0, weights_zeta0,                              &
    nquad_re0, nodes_re0, weights_re0,                                    &
    nitem0, npm0, re_item0, theta_time0, delta_free0,                     &
    monotone_tau0, update_gh0, fit_id0, Q0, loadmask0,             &
    p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0, vrais)             &
    BIND(C, name="score_aux_loglik_growth_re_tt_agh_full_equiv")

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0, ntre0
  INTEGER, INTENT(IN) :: nquad_u0, nquad_zeta0, nquad_re0
  INTEGER, INTENT(IN) :: npm0, re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0
  INTEGER, INTENT(IN) :: update_gh0, fit_id0, Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0, delta_x_col0, delta_z_col0
  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0), INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0), INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_u0
  DOUBLE PRECISION, DIMENSION(nquad_u0), INTENT(IN) :: weights_u0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_zeta0), INTENT(IN) :: weights_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: weights_re0
  DOUBLE PRECISION, INTENT(OUT) :: vrais

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.d9
  DOUBLE PRECISION :: ll
  LOGICAL :: ok_model, ok_grid, ok_ll

  vrais = BAD_OBJ
  CALL invalidate_full_cache(fit_id0)

  CALL prepare_growth_re_model(                                       &
       b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,      &
       nquad_u0, nodes_u0, weights_u0,                                &
       nquad_zeta0, nodes_zeta0, weights_zeta0,                       &
       nquad_re0, nodes_re0, weights_re0,                             &
       nitem0, npm0, re_item0, theta_time0, delta_free0,              &
       monotone_tau0, Q0, loadmask0,                              &
       p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0, ok_model)
  IF (.NOT. ok_model) GOTO 900

  IF (update_gh0 == 1 .OR. .NOT. gr_full_ready) THEN
    CALL build_full_cache(ok_grid)
    IF (.NOT. ok_grid) GOTO 900
  END IF

  CALL eval_full_loglik(ll, ok_ll)
  IF (ok_ll .AND. ll == ll .AND. ABS(ll) < 1.d300) vrais = -ll

900 CONTINUE
  CALL clear_growth_re_model()
END SUBROUTINE loglik_growth_re_tt_agh_full_equiv

SUBROUTINE clear_score_maps_q1()
  IMPLICIT NONE
  IF (ALLOCATED(sc_dtau)) DEALLOCATE(sc_dtau)
  IF (ALLOCATED(sc_dlambda)) DEALLOCATE(sc_dlambda)
  IF (ALLOCATED(sc_dgamma)) DEALLOCATE(sc_dgamma)
  IF (ALLOCATED(sc_dtheta)) DEALLOCATE(sc_dtheta)
  IF (ALLOCATED(sc_dkappa)) DEALLOCATE(sc_dkappa)
  IF (ALLOCATED(sc_dmap)) DEALLOCATE(sc_dmap)
  IF (ALLOCATED(sc_dzeta_sd)) DEALLOCATE(sc_dzeta_sd)
  IF (ALLOCATED(sc_dsqrt_sigmak)) DEALLOCATE(sc_dsqrt_sigmak)
  IF (ALLOCATED(sc_dL)) DEALLOCATE(sc_dL)
  sc_npm = 0
END SUBROUTINE clear_score_maps_q1


SUBROUTINE build_score_maps_q1(b0, npm0, ntre0, monotone_tau0, delta_free0, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npm0, ntre0, monotone_tau0, delta_free0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: ntau, dimU, pos, idx, pos_tau, k, m, h, j0, j1
  INTEGER :: qidx, t, r, rr, row, s, i, ref_item
  INTEGER :: ivar, isd, ipc
  DOUBLE PRECISION :: x, v1, sd1, sd2, pc, cs, zsd
  DOUBLE PRECISION :: delta_val, mask, xbase, xval, zbase, zval
  DOUBLE PRECISION :: dxcoef, dzcoef, score_base, score_adj, dscore
  DOUBLE PRECISION, ALLOCATABLE :: dbeta(:,:,:), ddelta(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: dLu(:,:,:), Lu0(:,:), acc(:)

  ok = .FALSE.

  IF (Q /= 1) RETURN
  IF (p_ran /= 2) RETURN
  IF (growth_ident_flag /= 1 .AND. growth_ident_flag /= 2 .AND. growth_ident_flag /= 3) RETURN

  dimU = Q * p_ran
  IF (dimU /= 2) RETURN

  CALL clear_score_maps_q1()
  sc_npm = npm0
  ntau = SUM(cm1)

  ALLOCATE(sc_dtau(ntau,npm0))
  ALLOCATE(sc_dlambda(nvar,Q,npm0))
  ALLOCATE(sc_dgamma(nvar,n,npm0))
  ALLOCATE(sc_dtheta(nvar,n,npm0))
  ALLOCATE(sc_dkappa(gr_dim_eta,growth_design_nsub,npm0))
  ALLOCATE(sc_dmap(gr_dim_eta,gr_dim_u,growth_design_nsub,npm0))
  ALLOCATE(sc_dzeta_sd(gr_dim_eta,npm0))
  ALLOCATE(sc_dsqrt_sigmak(nvar,npm0))

  ALLOCATE(dbeta(Q,p_fix,npm0), ddelta(Q,npm0))
  ALLOCATE(dLu(dimU,dimU,npm0), Lu0(dimU,dimU), acc(npm0))

  sc_dtau = 0.d0
  sc_dlambda = 0.d0
  sc_dgamma = 0.d0
  sc_dtheta = 0.d0
  sc_dkappa = 0.d0
  sc_dmap = 0.d0
  sc_dzeta_sd = 0.d0
  sc_dsqrt_sigmak = 0.d0
  dbeta = 0.d0
  ddelta = 0.d0
  dLu = 0.d0
  Lu0 = 0.d0

  idx = 0
  pos_tau = 0
  DO k = 1, nvar
    m = cm1(k)
    IF (m > 1) THEN
      IF (growth_threshold_id_flag == 0) THEN
        IF (monotone_tau0 == 1) THEN
          acc = 0.d0
          DO h = 3, m
            idx = idx + 1
            IF (idx > ntre0) GOTO 900
            acc(idx) = acc(idx) + EXP(b0(idx))
            sc_dtau(pos_tau+h,:) = acc
          END DO
        ELSE
          DO h = 3, m
            idx = idx + 1
            IF (idx > ntre0) GOTO 900
            sc_dtau(pos_tau+h,idx) = 1.d0
          END DO
        END IF
      ELSE
        j0 = (m + 1) / 2
        j1 = j0 + 1
        IF (monotone_tau0 == 1) THEN
          acc = 0.d0
          DO h = j0-1, 1, -1
            idx = idx + 1
            IF (idx > ntre0) GOTO 900
            acc(idx) = acc(idx) - EXP(b0(idx))
            sc_dtau(pos_tau+h,:) = acc
          END DO
          acc = 0.d0
          DO h = j1+1, m
            idx = idx + 1
            IF (idx > ntre0) GOTO 900
            acc(idx) = acc(idx) + EXP(b0(idx))
            sc_dtau(pos_tau+h,:) = acc
          END DO
        ELSE
          DO h = j0-1, 1, -1
            idx = idx + 1
            IF (idx > ntre0) GOTO 900
            sc_dtau(pos_tau+h,idx) = 1.d0
          END DO
          DO h = j1+1, m
            idx = idx + 1
            IF (idx > ntre0) GOTO 900
            sc_dtau(pos_tau+h,idx) = 1.d0
          END DO
        END IF
      END IF
    END IF
    pos_tau = pos_tau + m
  END DO
  IF (idx /= ntre0) GOTO 900

  pos = ntre0

  IF (growth_ident_flag == 2 .OR. growth_ident_flag == 3) THEN
    IF (.NOT. ALLOCATED(growth_marker_item)) GOTO 900
    IF (SIZE(growth_marker_item) /= Q) GOTO 900
    DO qidx = 1, Q
      ref_item = growth_marker_item(qidx)
      IF (ref_item < 1 .OR. ref_item > nvar) GOTO 900
      IF (loadmask(ref_item,qidx) /= 1) GOTO 900
      DO k = 1, nvar
        IF (loadmask(k,qidx) == 0 .OR. k == ref_item) CYCLE
        pos = pos + 1
        IF (pos > npm0) GOTO 900
        sc_dlambda(k,qidx,pos) = 1.d0
      END DO
    END DO
  ELSE
    DO qidx = 1, Q
      DO k = 1, nvar
        IF (loadmask(k,qidx) == 0) CYCLE
        pos = pos + 1
        IF (pos > npm0) GOTO 900
        sc_dlambda(k,qidx,pos) = 1.d0
      END DO
    END DO
  END IF

  IF (growth_ident_flag == 3) THEN
    IF (.NOT. ALLOCATED(growth_marker_item)) GOTO 900
    ref_item = growth_marker_item(1)
    IF (ref_item < 1 .OR. ref_item > nvar) GOTO 900
    DO k = 1, nvar
      IF (k == ref_item) CYCLE
      pos = pos + 1
      IF (pos > npm0) GOTO 900
      DO t = 1, n
        sc_dgamma(k,t,pos) = 1.d0
      END DO
    END DO
  ELSE
    DO k = 1, nvar
      pos = pos + 1
      IF (pos > npm0) GOTO 900
      DO t = 1, n
        sc_dgamma(k,t,pos) = 1.d0
      END DO
    END DO
  END IF

  DO t = 2, n
    DO k = 1, nvar
      IF (growth_gamma_time_free(k,t) == 1) THEN
        pos = pos + 1
        IF (pos > npm0) GOTO 900
        sc_dgamma(k,t,:) = 0.d0
        sc_dgamma(k,t,pos) = 1.d0
      END IF
    END DO
  END DO

  IF (theta_time_flag == 0) THEN
    DO k = 1, nvar
      IF (cm1(k) > 1) THEN
        pos = pos + 1
        IF (pos > npm0) GOTO 900
        DO t = 1, n
          sc_dtheta(k,t,pos) = theta_long(k,t)
        END DO
      END IF
    END DO
  ELSE
    DO k = 1, nvar
      IF (cm1(k) > 1) THEN
        DO t = 1, n
          pos = pos + 1
          IF (pos > npm0) GOTO 900
          sc_dtheta(k,t,pos) = theta_long(k,t)
        END DO
      END IF
    END DO
  END IF

  DO r = 1, p_fix
    DO qidx = 1, Q
      IF ((growth_ident_flag == 1 .OR. growth_ident_flag == 2) .AND. r == 1) CYCLE
      pos = pos + 1
      IF (pos > npm0) GOTO 900
      dbeta(qidx,r,pos) = 1.d0
    END DO
  END DO

  delta_val = 0.d0
  IF (delta_free0 == 1) THEN
    pos = pos + 1
    IF (pos > npm0) GOTO 900
    delta_val = b0(pos)
    ddelta(1,pos) = 1.d0
  ELSE IF (growth_delta_mode == 1 .AND. growth_rho_fixed_flag == 1) THEN
    IF (.NOT. ALLOCATED(growth_rho_fixed)) GOTO 900
    IF (SIZE(growth_rho_fixed) /= Q) GOTO 900
    delta_val = growth_rho_fixed(1)
  END IF

  IF (growth_ident_flag == 2 .OR. growth_ident_flag == 3) THEN
    pos = pos + 1
    IF (pos > npm0) GOTO 900
    ivar = pos
    Lu0(1,1) = EXP(b0(ivar))
    dLu(1,1,ivar) = Lu0(1,1)

    pos = pos + 1
    IF (pos > npm0) GOTO 900
    isd = pos
    Lu0(2,1) = b0(isd)
    dLu(2,1,isd) = 1.d0

    pos = pos + 1
    IF (pos > npm0) GOTO 900
    ipc = pos
    Lu0(2,2) = EXP(b0(ipc))
    dLu(2,2,ipc) = Lu0(2,2)

    DO t = 1, n
      pos = pos + 1
      IF (pos > npm0) GOTO 900
      row = t
      zsd = EXP(0.5d0*b0(pos))
      sc_dzeta_sd(row,pos) = 0.5d0 * zsd
    END DO
  ELSE
    pos = pos + 1
    IF (pos > npm0) GOTO 900
    ivar = pos
    x = b0(ivar)
    IF (x >= 0.d0) THEN
      v1 = 1.d0 / (1.d0 + EXP(-x))
    ELSE
      v1 = EXP(x) / (1.d0 + EXP(x))
    END IF
    sd1 = SQRT(v1)

    pos = pos + 1
    IF (pos > npm0) GOTO 900
    isd = pos
    sd2 = EXP(b0(isd))

    pos = pos + 1
    IF (pos > npm0) GOTO 900
    ipc = pos
    pc = TANH(b0(ipc))
    cs = SQRT(MAX(1.d-12, 1.d0-pc*pc))

    Lu0(1,1) = sd1
    Lu0(2,1) = sd2 * pc
    Lu0(2,2) = sd2 * cs

    dLu(1,1,ivar) = 0.5d0 * sd1 * (1.d0-v1)
    dLu(2,1,isd) = Lu0(2,1)
    dLu(2,1,ipc) = sd2 * (1.d0-pc*pc)
    dLu(2,2,isd) = Lu0(2,2)
    dLu(2,2,ipc) = -pc * Lu0(2,2)

    zsd = SQRT(MAX(1.d-12, 1.d0-v1))
    sc_dzeta_sd(1,ivar) = -0.5d0 * v1 * zsd

    DO t = 2, n
      pos = pos + 1
      IF (pos > npm0) GOTO 900
      row = t
      zsd = EXP(0.5d0*b0(pos))
      sc_dzeta_sd(row,pos) = 0.5d0 * zsd
    END DO
  END IF

  IF (re_item_flag == 1) THEN
    DO k = 1, nvar
      pos = pos + 1
      IF (pos > npm0) GOTO 900
      sc_dsqrt_sigmak(k,pos) = 0.5d0 * sqrt_sigmak(k)
    END DO
  END IF

  IF (pos /= npm0) GOTO 900

  DO i = 1, growth_design_nsub
    DO t = 1, n
      row = t
      mask = growth_delta_mask(i,t,1)

      DO r = 1, p_fix
        xbase = X_eta(i,t,1,r)
        xval = xbase
        dxcoef = 0.d0

        IF (growth_delta_x_col > 0 .AND. r == growth_delta_x_col) THEN
          IF (mask /= 0.d0) THEN
            IF (growth_delta_mode == 1) THEN
              xval = delta_val
              dxcoef = 1.d0
            ELSE
              xval = xbase + delta_val * mask
              dxcoef = mask
            END IF
          END IF
        ELSE IF (ALLOCATED(growth_x_slope_col)) THEN
          IF (growth_delta_x_col > 0 .AND. growth_x_slope_col(r) == 1) THEN
            score_base = X_eta(i,t,1,growth_delta_x_col)
            score_adj = score_base
            dscore = 0.d0
            IF (mask /= 0.d0) THEN
              IF (growth_delta_mode == 1) THEN
                score_adj = delta_val
                dscore = 1.d0
              ELSE
                score_adj = score_base + delta_val * mask
                dscore = mask
              END IF
            END IF
            xval = xbase * score_adj
            dxcoef = xbase * dscore
          END IF
        END IF

        sc_dkappa(row,i,:) = sc_dkappa(row,i,:) + &
             xval * dbeta(1,r,:) + beta(1,r) * dxcoef * ddelta(1,:)
      END DO

      DO rr = 1, p_ran
        zbase = Z_eta(i,t,1,rr)
        zval = zbase
        dzcoef = 0.d0
        IF (growth_delta_z_col > 0 .AND. rr == growth_delta_z_col) THEN
          IF (mask /= 0.d0) THEN
            IF (growth_delta_mode == 1) THEN
              zval = delta_val
              dzcoef = 1.d0
            ELSE
              zval = zbase + delta_val * mask
              dzcoef = mask
            END IF
          END IF
        END IF

        DO s = 1, dimU
          sc_dmap(row,s,i,:) = sc_dmap(row,s,i,:) + &
               zval * dLu(rr,s,:) + Lu0(rr,s) * dzcoef * ddelta(1,:)
        END DO
      END DO
    END DO
  END DO

  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(dbeta, ddelta, dLu, Lu0, acc)
  IF (.NOT. ok) CALL clear_score_maps_q1()
END SUBROUTINE build_score_maps_q1


SUBROUTINE item_logh_score_q1(i, k_item, eta_vec, deta, z_a, npm0, logh, score, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i, k_item, npm0
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: eta_vec
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: deta
  DOUBLE PRECISION, INTENT(IN) :: z_a
  DOUBLE PRECISION, INTENT(OUT) :: logh
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(OUT) :: score
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: kk_thr, idx_obs, kk, j_mes, t_idx, yijk, tuk, qq, idx_tq
  DOUBLE PRECISION :: a_item, mu, var_e, invsd, au, bl, p, phiu, phil
  DOUBLE PRECISION, PARAMETER :: INV_SQRT_2PI = 0.39894228040143267794d0
  DOUBLE PRECISION, DIMENSION(npm0) :: da_item, dmu, dinvsd, dau, dbl, dp

  ok = .FALSE.
  logh = 0.d0
  score = 0.d0

  kk_thr = 0
  idx_obs = nmescur
  IF (k_item > 1) THEN
    DO kk = 1, k_item-1
      kk_thr = kk_thr + (c(kk)-1)
      idx_obs = idx_obs + nmes(i,kk)
    END DO
  END IF

  IF (nmes(i,k_item) == 0) THEN
    ok = .TRUE.
    RETURN
  END IF

  a_item = 0.d0
  da_item = 0.d0
  IF (re_item_flag == 1) THEN
    a_item = sqrt_sigmak(k_item) * z_a
    da_item = z_a * sc_dsqrt_sigmak(k_item,:)
  END IF

  DO j_mes = 1, nmes(i,k_item)
    idx_obs = idx_obs + 1
    t_idx = nitem(idx_obs)

    var_e = theta_at(k_item,t_idx)
    IF (var_e <= 0.d0 .OR. var_e /= var_e) RETURN
    invsd = 1.d0 / SQRT(var_e)
    dinvsd = -0.5d0 * invsd / var_e * sc_dtheta(k_item,t_idx,:)

    mu = gamma_at(k_item,t_idx) + a_item
    dmu = sc_dgamma(k_item,t_idx,:) + da_item

    DO qq = 1, Q
      IF (loadmask(k_item,qq) == 0) CYCLE
      idx_tq = (t_idx-1)*Q + qq
      mu = mu + lambda_at(k_item,qq,t_idx) * eta_vec(idx_tq)
      dmu = dmu + sc_dlambda(k_item,qq,:) * eta_vec(idx_tq) + &
             lambda_at(k_item,qq,t_idx) * deta(idx_tq,:)
    END DO

    yijk = Y(idx_obs)
    tuk = yijk + kk_thr
    dp = 0.d0

    IF (yijk == 1) THEN
      au = (tau_at(kk_thr+1,t_idx)-mu) * invsd
      dau = (sc_dtau(kk_thr+1,:) - dmu) * invsd + &
            (tau_at(kk_thr+1,t_idx)-mu) * dinvsd
      p = alnorm_fast(au,.FALSE.)
      phiu = INV_SQRT_2PI * EXP(-0.5d0*au*au)
      dp = phiu * dau
    ELSE IF (yijk == c(k_item)) THEN
      bl = (tau_at(tuk-1,t_idx)-mu) * invsd
      dbl = (sc_dtau(tuk-1,:) - dmu) * invsd + &
            (tau_at(tuk-1,t_idx)-mu) * dinvsd
      p = alnorm_fast(bl,.TRUE.)
      phil = INV_SQRT_2PI * EXP(-0.5d0*bl*bl)
      dp = -phil * dbl
    ELSE
      au = (tau_at(tuk,t_idx)-mu) * invsd
      bl = (tau_at(tuk-1,t_idx)-mu) * invsd
      dau = (sc_dtau(tuk,:) - dmu) * invsd + &
            (tau_at(tuk,t_idx)-mu) * dinvsd
      dbl = (sc_dtau(tuk-1,:) - dmu) * invsd + &
            (tau_at(tuk-1,t_idx)-mu) * dinvsd
      IF (bl > 0.d0) THEN
        p = alnorm_fast(bl,.TRUE.) - alnorm_fast(au,.TRUE.)
      ELSE
        p = alnorm_fast(au,.FALSE.) - alnorm_fast(bl,.FALSE.)
      END IF
      phiu = INV_SQRT_2PI * EXP(-0.5d0*au*au)
      phil = INV_SQRT_2PI * EXP(-0.5d0*bl*bl)
      dp = phiu*dau - phil*dbl
    END IF

    IF (p <= 0.d0 .OR. p /= p) RETURN
    logh = logh + LOG(p)
    score = score + dp / p
  END DO

  ok = .TRUE.
END SUBROUTINE item_logh_score_q1


SUBROUTINE inneragh_item_logint_full_equiv_score_q1( &
    i, k_item, eta_vec, deta, alpha_z, sdev_z, npm0, logint, score, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i, k_item, npm0
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: eta_vec
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: deta
  DOUBLE PRECISION, INTENT(IN) :: alpha_z, sdev_z
  DOUBLE PRECISION, INTENT(OUT) :: logint
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(OUT) :: score
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: r, nvalid
  DOUBLE PRECISION :: v, z_a, logh, term, maxterm, sumexp, w
  DOUBLE PRECISION, DIMENSION(npm0) :: sr, ssum
  LOGICAL :: okh
  DOUBLE PRECISION, PARAMETER :: LOG_ZERO = -1.d9

  ok = .FALSE.
  logint = LOG_ZERO
  score = 0.d0

  IF (nmes(i,k_item) == 0) THEN
    logint = 0.d0
    ok = .TRUE.
    RETURN
  END IF

  IF (re_item_flag == 0) THEN
    CALL item_logh_score_q1(i,k_item,eta_vec,deta,0.d0,npm0,logint,score,ok)
    RETURN
  END IF

  IF (sdev_z <= 0.d0 .OR. sdev_z /= sdev_z) RETURN

  nvalid = 0
  maxterm = 0.d0
  sumexp = 0.d0
  ssum = 0.d0

  DO r = 1, nquad_re
    IF (weights_re(r) <= 0.d0) CYCLE
    v = nodes_re(r)
    z_a = alpha_z + sdev_z*v

    CALL item_logh_score_q1(i,k_item,eta_vec,deta,z_a,npm0,logh,sr,okh)
    IF (.NOT. okh .OR. logh /= logh) CYCLE

    term = LOG(weights_re(r)) + LOG(sdev_z) + logh - &
           0.5d0*(z_a*z_a-v*v)
    IF (term /= term) CYCLE

    IF (nvalid == 0) THEN
      maxterm = term
      sumexp = 1.d0
      ssum = sr
      nvalid = 1
    ELSE IF (term > maxterm) THEN
      w = EXP(maxterm-term)
      sumexp = sumexp*w + 1.d0
      ssum = ssum*w + sr
      maxterm = term
      nvalid = nvalid + 1
    ELSE
      w = EXP(term-maxterm)
      sumexp = sumexp + w
      ssum = ssum + w*sr
      nvalid = nvalid + 1
    END IF
  END DO

  IF (nvalid == 0 .OR. sumexp <= 0.d0 .OR. sumexp /= sumexp) RETURN

  logint = maxterm + LOG(sumexp)
  score = ssum / sumexp
  ok = .TRUE.
END SUBROUTINE inneragh_item_logint_full_equiv_score_q1


SUBROUTINE eval_full_loglik_score_q1(npm0, out, score, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npm0
  DOUBLE PRECISION, INTENT(OUT) :: out
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(OUT) :: score
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: i, lu, lz, k_item, sm, nvalid, dimZ, row, s, design_i
  DOUBLE PRECISION :: logg, logint, term, maxterm, sumexp
  DOUBLE PRECISION :: qw, qv, alpha_k, sdev_k, w
  DOUBLE PRECISION, ALLOCATABLE :: v_w(:), z_w(:), eta_vec(:), deta(:,:)
  DOUBLE PRECISION, DIMENSION(npm0) :: snode, sitem, ssum
  LOGICAL :: ok_item

  out = 0.d0
  score = 0.d0
  ok = .FALSE.
  IF (.NOT. gr_full_ready) RETURN
  IF (sc_npm /= npm0) RETURN

  dimZ = gr_full_dim_z
  ALLOCATE(v_w(gr_dim_w), z_w(gr_dim_w), eta_vec(gr_dim_eta))
  ALLOCATE(deta(gr_dim_eta,npm0))
  nmescur = 0

  DO i = 1, nsub
    sm = SUM(nmes(i,:))
    design_i = MERGE(i,1,growth_subject_specific_design == 1)
    nvalid = 0
    maxterm = 0.d0
    sumexp = 0.d0
    ssum = 0.d0

    DO lu = 1, gr_nquad_u
      IF (gr_weights_u(lu) <= 0.d0) CYCLE
      DO lz = 1, gr_nquad_zeta
        IF (gr_weights_zeta(lz) <= 0.d0) CYCLE

        v_w(1:gr_dim_u) = gr_nodes_u(:,lu)
        v_w(gr_dim_u+1:gr_dim_w) = gr_nodes_zeta(:,lz)
        z_w = gr_full_mode(1:gr_dim_w,i) + &
              MATMUL(gr_full_B(1:gr_dim_w,1:gr_dim_w,i),v_w)
        CALL eta_from_w(i,z_w,eta_vec)

        deta = 0.d0
        DO row = 1, gr_dim_eta
          deta(row,:) = sc_dkappa(row,design_i,:) + &
                        sc_dzeta_sd(row,:) * z_w(gr_dim_u+row)
          DO s = 1, gr_dim_u
            deta(row,:) = deta(row,:) + sc_dmap(row,s,design_i,:) * z_w(s)
          END DO
        END DO

        logg = 0.d0
        snode = 0.d0

        DO k_item = 1, nvar
          IF (re_item_flag == 1) THEN
            alpha_k = gr_full_mode(gr_dim_w+k_item,i) + &
                 DOT_PRODUCT(gr_full_B(gr_dim_w+k_item,1:gr_dim_w,i),v_w)
            sdev_k = gr_full_B(gr_dim_w+k_item,gr_dim_w+k_item,i)
          ELSE
            alpha_k = 0.d0
            sdev_k = 1.d0
          END IF

          CALL inneragh_item_logint_full_equiv_score_q1( &
               i,k_item,eta_vec,deta,alpha_k,sdev_k,npm0,logint,sitem,ok_item)
          IF (.NOT. ok_item) THEN
            logg = -1.d9
            EXIT
          END IF
          logg = logg + logint
          snode = snode + sitem
        END DO
        IF (logg <= -1.d9) CYCLE

        qw = DOT_PRODUCT(z_w,z_w)
        qv = DOT_PRODUCT(v_w,v_w)
        term = LOG(gr_weights_u(lu)) + LOG(gr_weights_zeta(lz)) + &
               logg - 0.5d0*(qw-qv)
        IF (term /= term) CYCLE

        IF (nvalid == 0) THEN
          maxterm = term
          sumexp = 1.d0
          ssum = snode
          nvalid = 1
        ELSE IF (term > maxterm) THEN
          w = EXP(maxterm-term)
          sumexp = sumexp*w + 1.d0
          ssum = ssum*w + snode
          maxterm = term
          nvalid = nvalid + 1
        ELSE
          w = EXP(term-maxterm)
          sumexp = sumexp + w
          ssum = ssum + w*snode
          nvalid = nvalid + 1
        END IF
      END DO
    END DO

    IF (nvalid == 0 .OR. sumexp <= 0.d0) GOTO 900
    out = out + gr_full_logdet_w(i) + maxterm + LOG(sumexp)
    score = score + ssum/sumexp
    nmescur = nmescur + sm
  END DO

  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(v_w,z_w,eta_vec,deta)
END SUBROUTINE eval_full_loglik_score_q1


SUBROUTINE loglik_growth_re_tt_agh_full_equiv_score_q1(                    &
    b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,              &
    nquad_u0, nodes_u0, weights_u0,                                        &
    nquad_zeta0, nodes_zeta0, weights_zeta0,                               &
    nquad_re0, nodes_re0, weights_re0,                                     &
    nitem0, npm0, re_item0, theta_time0, delta_free0,                      &
    monotone_tau0, update_gh0, fit_id0, Q0, loadmask0,                     &
    p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0,       &
    vrais, score0, ok0)                                                     &
    BIND(C, name="loglik_growth_re_tt_agh_full_equiv_score_q1")

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0, ntre0
  INTEGER, INTENT(IN) :: nquad_u0, nquad_zeta0, nquad_re0
  INTEGER, INTENT(IN) :: npm0, re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0, update_gh0, fit_id0, Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0, delta_x_col0, delta_z_col0
  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0), INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0), INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_u0
  DOUBLE PRECISION, DIMENSION(nquad_u0), INTENT(IN) :: weights_u0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_zeta0), INTENT(IN) :: weights_zeta0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: weights_re0
  DOUBLE PRECISION, INTENT(OUT) :: vrais
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(OUT) :: score0
  INTEGER, INTENT(OUT) :: ok0

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.d9
  DOUBLE PRECISION :: ll
  DOUBLE PRECISION, DIMENSION(npm0) :: score_ll
  LOGICAL :: ok_model, ok_grid, ok_maps, ok_ll

  vrais = BAD_OBJ
  score0 = 0.d0
  ok0 = 0

  CALL invalidate_full_cache(fit_id0)

  CALL prepare_growth_re_model(                                           &
       b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,          &
       nquad_u0, nodes_u0, weights_u0,                                    &
       nquad_zeta0, nodes_zeta0, weights_zeta0,                           &
       nquad_re0, nodes_re0, weights_re0,                                 &
       nitem0, npm0, re_item0, theta_time0, delta_free0,                  &
       monotone_tau0, Q0, loadmask0,                                      &
       p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0,   &
       ok_model)
  IF (.NOT. ok_model) GOTO 900

  IF (Q /= 1 .OR. p_ran /= 2) GOTO 900
  IF (growth_ident_flag /= 1 .AND. growth_ident_flag /= 2 .AND. growth_ident_flag /= 3) GOTO 900

  IF (update_gh0 == 1 .OR. .NOT. gr_full_ready) THEN
    CALL build_full_cache(ok_grid)
    IF (.NOT. ok_grid) GOTO 900
  END IF

  CALL build_score_maps_q1(b0,npm0,ntre0,monotone_tau0,delta_free0,ok_maps)
  IF (.NOT. ok_maps) GOTO 900

  CALL eval_full_loglik_score_q1(npm0,ll,score_ll,ok_ll)
  IF (.NOT. ok_ll) GOTO 900
  IF (ll /= ll .OR. ABS(ll) > 1.d300) GOTO 900

  vrais = -ll
  score0 = -score_ll
  ok0 = 1

900 CONTINUE
  CALL clear_score_maps_q1()
  CALL clear_growth_re_model()
END SUBROUTINE loglik_growth_re_tt_agh_full_equiv_score_q1


SUBROUTINE set_growth_re_agh_full_equiv_score_grid(                       &
    mode_in, B_in, logdet_w_in,                                           &
    dim_w0, dim_z0, nsub0, nquad_u0, nquad_zeta0, fit_id0, ready0)        &
    BIND(C, name="set_growth_re_agh_full_equiv_score_grid")

  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE

  REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: mode_in
  REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: B_in
  REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: logdet_w_in

  INTEGER(C_INT), INTENT(IN) :: dim_w0
  INTEGER(C_INT), INTENT(IN) :: dim_z0
  INTEGER(C_INT), INTENT(IN) :: nsub0
  INTEGER(C_INT), INTENT(IN) :: nquad_u0
  INTEGER(C_INT), INTENT(IN) :: nquad_zeta0
  INTEGER(C_INT), INTENT(IN) :: fit_id0
  INTEGER(C_INT), INTENT(OUT) :: ready0

  INTEGER :: i, j, k, pos

  ready0 = 0_C_INT

  IF (dim_w0 < 1 .OR. dim_z0 < 1 .OR. nsub0 < 1) RETURN
  IF (dim_w0 > dim_z0) RETURN

  IF (ALLOCATED(gr_full_mode)) DEALLOCATE(gr_full_mode)
  IF (ALLOCATED(gr_full_B)) DEALLOCATE(gr_full_B)
  IF (ALLOCATED(gr_full_logdet_w)) DEALLOCATE(gr_full_logdet_w)

  ALLOCATE(gr_full_mode(dim_z0,nsub0))
  ALLOCATE(gr_full_B(dim_z0,dim_z0,nsub0))
  ALLOCATE(gr_full_logdet_w(nsub0))

  pos = 0
  DO i = 1, nsub0
    DO j = 1, dim_z0
      pos = pos + 1
      gr_full_mode(j,i) = mode_in(pos)
    END DO
  END DO

  pos = 0
  DO i = 1, nsub0
    DO k = 1, dim_z0
      DO j = 1, dim_z0
        pos = pos + 1
        gr_full_B(j,k,i) = B_in(pos)
      END DO
    END DO
  END DO

  DO i = 1, nsub0
    gr_full_logdet_w(i) = logdet_w_in(i)
  END DO

  gr_full_dim_w = dim_w0
  gr_full_dim_z = dim_z0
  gr_full_nsub = nsub0
  gr_full_nquad_u = nquad_u0
  gr_full_nquad_zeta = nquad_zeta0
  gr_full_fit_id = fit_id0
  gr_full_ready = .TRUE.

  ready0 = 1_C_INT

END SUBROUTINE set_growth_re_agh_full_equiv_score_grid



SUBROUTINE prepare_growth_full_equiv_score_model_q1(                       &
    b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,             &
    nquad_eta0, nodes_eta0, weights_eta0,                                  &
    nquad_re0, nodes_re0, weights_re0,                                     &
    nitem0, npm0, re_item0, theta_time0, delta_free0,                      &
    monotone_tau0, Q0, loadmask0,                                          &
    p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0, ok)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0, ntre0
  INTEGER, INTENT(IN) :: nquad_eta0, nquad_re0
  INTEGER, INTENT(IN) :: npm0, re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0, Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0, delta_x_col0, delta_z_col0
  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0), INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0), INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_eta0
  DOUBLE PRECISION, DIMENSION(nquad_eta0), INTENT(IN) :: weights_eta0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: nodes_re0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: weights_re0
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: k, t, r, row, pos, info, qidx, dimEta, dimU, nX, nZ, rr, s
  INTEGER :: ntheta_est, expected_npm, ntre_calc, n_lambda_free, n_gamma_free
  INTEGER :: n_Lu_est, base, zcol
  INTEGER, ALLOCATABLE :: cm3(:)
  DOUBLE PRECISION :: zval
  DOUBLE PRECISION, ALLOCATABLE :: delta_g(:), Lu(:,:), Sigma_u(:,:), zeta_var(:,:)

  ok = .FALSE.
  IF (n0 < 2 .OR. Q0 < 1 .OR. ABS(p_fix0) < 1 .OR. p_ran0 < 1) RETURN
  IF (nquad_eta0 < 1 .OR. nquad_re0 < 1) RETURN

  CALL clear_growth_re_model()

  nvar = nvar0
  nsub = nsub0
  nobs = nobs0
  n = n0
  Q = Q0
  npm = npm0

  dimEta = n * Q
  dimU = Q * p_ran0
  dim_eta = dimEta
  gr_dim_eta = dimEta
  gr_dim_u = dimU
  gr_dim_w = dimEta

  nquad = nquad_eta0
  nquad_eta = nquad_eta0
  nquad_re = nquad_re0

  re_item_flag = re_item0
  theta_time_flag = theta_time0
  phi_struct_flag = 0
  growth_subject_specific_design = MERGE(1, 0, p_fix0 > 0)
  growth_design_nsub = MERGE(nsub0, 1, growth_subject_specific_design == 1)
  p_fix = ABS(p_fix0)
  p_ran = p_ran0
  growth_delta_x_col = delta_x_col0
  growth_delta_z_col = delta_z_col0

  ALLOCATE(nmes(nsub,nvar), Y(nobs), nitem(nobs))
  ALLOCATE(c(nvar), cm1(nvar))
  ALLOCATE(nodes(nquad*dimEta), weights(nquad))
  ALLOCATE(nodes_re(nquad_re), weights_re(nquad_re), logweights_re(nquad_re))
  ALLOCATE(loadmask(nvar,Q))
  ALLOCATE(lambda(nvar,Q), gamma(nvar))
  ALLOCATE(theta_k(nvar), invsd_k(nvar))
  ALLOCATE(theta_long(nvar,n), invsd_long(nvar,n))
  ALLOCATE(sigmak(nvar), sqrt_sigmak(nvar))
  ALLOCATE(kappa(dimEta), L_phi(dimEta,dimEta))
  ALLOCATE(kk_base(nvar))
  ALLOCATE(tau(SUM(cm10)))
  ALLOCATE(cm3(nvar))
  ALLOCATE(delta_g(Q))
  ALLOCATE(beta(Q,p_fix), X_eta(growth_design_nsub,n,Q,p_fix), &
           Z_eta(growth_design_nsub,n,Q,p_ran))
  ALLOCATE(growth_delta_mask(growth_design_nsub,n,Q))
  ALLOCATE(Lu(dimU,dimU), Sigma_u(dimU,dimU), zeta_var(n,Q))
  ALLOCATE(gr_map_u(dimEta,dimU,growth_design_nsub), gr_zeta_sd(dimEta))

  nmes = nmes0
  Y = Y0
  nitem = nitem0
  c = c0
  cm1 = cm10
  loadmask = loadmask0

  nX = growth_design_nsub * n * Q * p_fix
  nZ = growth_design_nsub * n * Q * p_ran
  X_eta = RESHAPE(X0(1:nX), (/ growth_design_nsub, n, Q, p_fix /))
  Z_eta = RESHAPE(Z0(1:nZ), (/ growth_design_nsub, n, Q, p_ran /))
  growth_delta_mask = RESHAPE(delta_mask0(1:growth_design_nsub*n*Q), &
                              (/ growth_design_nsub, n, Q /))

  nodes(1:nquad*dimEta) = nodes_eta0(1:nquad*dimEta)
  weights = weights_eta0
  nodes_re = nodes_re0
  weights_re = weights_re0
  DO r = 1, nquad_re
    IF (weights_re(r) > 0.d0) THEN
      logweights_re(r) = LOG(weights_re(r))
    ELSE
      logweights_re(r) = -HUGE(1.d0)
    END IF
  END DO

  CALL setup_load_indices_from_loadmask_growth(info)
  IF (info /= 0) GOTO 900

  DO k = 1, nvar
    cm3(k) = MAX(cm1(k)-2,0)
  END DO
  ntre_calc = SUM(cm3)
  IF (ntre_calc /= ntre0) GOTO 900

  CALL decode_growth_thresholds(b0,npm0,ntre0,monotone_tau0,info)
  IF (info /= 0) GOTO 900

  pos = ntre0
  CALL decode_lambda_sum_constraints_growth(b0,npm0,pos,info)
  IF (info /= 0) GOTO 900
  CALL refresh_lambda_act_growth()

  CALL decode_gamma_sum_constraints_growth(b0,npm0,pos,info)
  IF (info /= 0) GOTO 900

  theta_k = 1.d0
  invsd_k = 1.d0
  theta_long = 1.d0
  invsd_long = 1.d0

  IF (theta_time_flag == 0) THEN
    DO k = 1, nvar
      IF (cm1(k) > 1) THEN
        pos = pos + 1
        theta_k(k) = EXP(b0(pos))
        invsd_k(k) = 1.d0 / SQRT(theta_k(k))
        DO t = 1, n
          theta_long(k,t) = theta_k(k)
          invsd_long(k,t) = invsd_k(k)
        END DO
      END IF
    END DO
  ELSE
    DO k = 1, nvar
      IF (cm1(k) > 1) THEN
        DO t = 1, n
          pos = pos + 1
          theta_long(k,t) = EXP(b0(pos))
          invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
        END DO
        theta_k(k) = theta_long(k,1)
        invsd_k(k) = invsd_long(k,1)
      END IF
    END DO
  END IF

  CALL decode_growth_beta(b0,npm0,pos,info)
  IF (info /= 0) GOTO 900
  CALL decode_growth_delta(b0,npm0,pos,delta_free0,delta_g,info)
  IF (info /= 0) GOTO 900
  CALL decode_growth_random_structure(b0,npm0,pos,dimU,Lu,Sigma_u,zeta_var,info)
  IF (info /= 0) GOTO 900

  IF (re_item_flag == 1) THEN
    DO k = 1, nvar
      pos = pos + 1
      IF (pos > npm0) GOTO 900
      sigmak(k) = EXP(b0(pos))
      sqrt_sigmak(k) = SQRT(sigmak(k))
    END DO
  ELSE
    sigmak = 0.d0
    sqrt_sigmak = 0.d0
  END IF

  ntheta_est = 0
  IF (theta_time_flag == 0) THEN
    DO k = 1, nvar
      IF (cm1(k) > 1) ntheta_est = ntheta_est + 1
    END DO
  ELSE
    DO k = 1, nvar
      IF (cm1(k) > 1) ntheta_est = ntheta_est + n
    END DO
  END IF

  n_lambda_free = growth_n_lambda_par()
  n_gamma_free = growth_n_gamma_par()
  n_Lu_est = dimU * (dimU + 1) / 2
  expected_npm = ntre0 + n_lambda_free + n_gamma_free + ntheta_est &
       + growth_n_beta_par() + MERGE(Q,0,delta_free0 == 1) + n_Lu_est &
       + growth_n_zeta_par() + MERGE(nvar,0,re_item_flag == 1)
  IF (npm0 /= expected_npm .OR. pos /= npm0) GOTO 900

  CALL build_growth_eta_moments_design(delta_g,Sigma_u,zeta_var,info)
  IF (info /= 0) GOTO 900

  gr_map_u = 0.d0
  DO k = 1, growth_design_nsub
    DO t = 1, n
      DO qidx = 1, Q
        row = (t-1)*Q + qidx
        DO rr = 1, p_ran
          zval = Z_eta(k,t,qidx,rr)
          IF (rr == growth_delta_z_col .AND. growth_delta_mask(k,t,qidx) /= 0.d0) THEN
            IF (growth_delta_mode == 1) THEN
              zval = delta_g(qidx)
            ELSE
              zval = zval + delta_g(qidx) * growth_delta_mask(k,t,qidx)
            END IF
          END IF
          zcol = (rr-1)*Q + qidx
          DO s = 1, dimU
            gr_map_u(row,s,k) = gr_map_u(row,s,k) + zval * Lu(zcol,s)
          END DO
        END DO
        gr_zeta_sd(row) = SQRT(zeta_var(t,qidx))
      END DO
    END DO
  END DO

  ok = .TRUE.

900 CONTINUE
  IF (ALLOCATED(cm3)) DEALLOCATE(cm3)
  IF (ALLOCATED(delta_g)) DEALLOCATE(delta_g)
  IF (ALLOCATED(Lu)) DEALLOCATE(Lu)
  IF (ALLOCATED(Sigma_u)) DEALLOCATE(Sigma_u)
  IF (ALLOCATED(zeta_var)) DEALLOCATE(zeta_var)
  IF (.NOT. ok) CALL clear_growth_re_model()
END SUBROUTINE prepare_growth_full_equiv_score_model_q1


SUBROUTINE chol_deriv_lower_q1(L, dV, dL, ok)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: L, dV
  DOUBLE PRECISION, DIMENSION(SIZE(L,1),SIZE(L,2)), INTENT(OUT) :: dL
  LOGICAL, INTENT(OUT) :: ok
  INTEGER :: d, i, j, k
  DOUBLE PRECISION :: tmp

  d = SIZE(L,1)
  ok = .FALSE.
  dL = 0.d0
  IF (SIZE(L,2) /= d .OR. SIZE(dV,1) /= d .OR. SIZE(dV,2) /= d) RETURN

  DO j = 1, d
    IF (L(j,j) <= 0.d0 .OR. L(j,j) /= L(j,j)) RETURN
    tmp = dV(j,j)
    DO k = 1, j-1
      tmp = tmp - 2.d0 * L(j,k) * dL(j,k)
    END DO
    dL(j,j) = tmp / (2.d0 * L(j,j))

    DO i = j+1, d
      tmp = dV(i,j) - L(i,j) * dL(j,j)
      DO k = 1, j-1
        tmp = tmp - dL(i,k)*L(j,k) - L(i,k)*dL(j,k)
      END DO
      dL(i,j) = tmp / L(j,j)
    END DO
  END DO
  ok = .TRUE.
END SUBROUTINE chol_deriv_lower_q1


SUBROUTINE build_chol_score_maps_q1(npm0, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npm0
  LOGICAL, INTENT(OUT) :: ok
  INTEGER :: i, j, r, c0, s, d
  DOUBLE PRECISION, ALLOCATABLE :: dV(:,:), dLc(:,:)
  LOGICAL :: okc

  ok = .FALSE.
  d = gr_dim_eta
  IF (d < 1 .OR. sc_npm /= npm0) RETURN
  IF (.NOT. ALLOCATED(L_phi_sub)) RETURN
  IF (.NOT. ALLOCATED(gr_map_u) .OR. .NOT. ALLOCATED(gr_zeta_sd)) RETURN
  IF (SIZE(L_phi_sub,1) /= d .OR. SIZE(L_phi_sub,2) /= d) RETURN

  IF (ALLOCATED(sc_dL)) DEALLOCATE(sc_dL)
  ALLOCATE(sc_dL(d,d,growth_design_nsub,npm0))
  sc_dL = 0.d0
  ALLOCATE(dV(d,d), dLc(d,d))

  DO i = 1, growth_design_nsub
    DO j = 1, npm0
      IF (.NOT. ANY(sc_dmap(:,:,i,j) /= 0.d0) .AND. &
          .NOT. ANY(sc_dzeta_sd(:,j) /= 0.d0)) CYCLE

      dV = 0.d0
      DO r = 1, d
        DO c0 = 1, r
          DO s = 1, gr_dim_u
            dV(r,c0) = dV(r,c0) + sc_dmap(r,s,i,j)*gr_map_u(c0,s,i) &
                 + gr_map_u(r,s,i)*sc_dmap(c0,s,i,j)
          END DO
          IF (r == c0) dV(r,r) = dV(r,r) + &
               2.d0*gr_zeta_sd(r)*sc_dzeta_sd(r,j)
          dV(c0,r) = dV(r,c0)
        END DO
      END DO

      CALL chol_deriv_lower_q1(L_phi_sub(:,:,i),dV,dLc,okc)
      IF (.NOT. okc) GOTO 900
      sc_dL(:,:,i,j) = dLc
    END DO
  END DO

  ok = .TRUE.
900 CONTINUE
  DEALLOCATE(dV,dLc)
  IF (.NOT. ok) THEN
    IF (ALLOCATED(sc_dL)) DEALLOCATE(sc_dL)
  END IF
END SUBROUTINE build_chol_score_maps_q1


SUBROUTINE eval_eta_full_loglik_score_q1(npm0, out, score, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npm0
  DOUBLE PRECISION, INTENT(OUT) :: out
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(OUT) :: score
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: i, l, d, j, k_item, sm, nvalid, dimEta, baseZ, design_i
  DOUBLE PRECISION :: q_eta, q_u, logg, logint, term, maxterm, sumexp, w
  DOUBLE PRECISION :: alpha_k, sdev_k
  DOUBLE PRECISION, ALLOCATABLE :: z_eta(:), u_eta(:), eta_vec(:), deta(:,:)
  DOUBLE PRECISION, DIMENSION(npm0) :: snode, sitem, ssum
  LOGICAL :: ok_item

  out = 0.d0
  score = 0.d0
  ok = .FALSE.
  IF (.NOT. agh_tt_fe_grid_ready) RETURN
  IF (sc_npm /= npm0 .OR. .NOT. ALLOCATED(sc_dL)) RETURN

  dimEta = n * Q
  IF (agh_tt_fe_dimEta_cached /= dimEta) RETURN
  IF (agh_tt_fe_nquad_eta_cached /= nquad) RETURN
  IF (agh_tt_fe_nquad_re_cached /= nquad_re) RETURN
  IF (agh_tt_fe_nsub_cached /= nsub) RETURN

  ALLOCATE(z_eta(dimEta),u_eta(dimEta),eta_vec(dimEta),deta(dimEta,npm0))
  nmescur = 0

  DO i = 1, nsub
    sm = SUM(nmes(i,:))
    design_i = MERGE(i,1,growth_subject_specific_design == 1)
    CALL activate_growth_subject(i)
    nvalid = 0
    maxterm = 0.d0
    sumexp = 0.d0
    ssum = 0.d0

    DO l = 1, nquad
      IF (weights(l) <= 0.d0) CYCLE
      baseZ = (l-1)*dimEta
      DO d = 1, dimEta
        u_eta(d) = nodes(baseZ+d)
      END DO
      z_eta = agh_tt_fe_zeta_eta(:,l,i)

      CALL cfa_eta_mean(i,eta_vec)
      DO d = 1, dimEta
        DO j = 1, d
          eta_vec(d) = eta_vec(d) + L_phi(d,j)*z_eta(j)
        END DO
      END DO

      deta = sc_dkappa(:,design_i,:)
      DO d = 1, dimEta
        DO j = 1, dimEta
          deta(d,:) = deta(d,:) + sc_dL(d,j,design_i,:)*z_eta(j)
        END DO
      END DO

      logg = 0.d0
      snode = 0.d0
      DO k_item = 1, nvar
        IF (re_item_flag == 1) THEN
          alpha_k = agh_tt_fe_alpha(k_item,l,i)
          sdev_k = agh_tt_fe_sdev(k_item,i)
        ELSE
          alpha_k = 0.d0
          sdev_k = 1.d0
        END IF

        CALL inneragh_item_logint_full_equiv_score_q1( &
             i,k_item,eta_vec,deta,alpha_k,sdev_k,npm0,logint,sitem,ok_item)
        IF (.NOT. ok_item) THEN
          logg = -1.d9
          EXIT
        END IF
        logg = logg + logint
        snode = snode + sitem
      END DO
      IF (logg <= -1.d9) CYCLE

      q_eta = DOT_PRODUCT(z_eta,z_eta)
      q_u = DOT_PRODUCT(u_eta,u_eta)
      term = LOG(weights(l)) + logg - 0.5d0*(q_eta-q_u)
      IF (term /= term) CYCLE

      IF (nvalid == 0) THEN
        maxterm = term
        sumexp = 1.d0
        ssum = snode
        nvalid = 1
      ELSE IF (term > maxterm) THEN
        w = EXP(maxterm-term)
        sumexp = sumexp*w + 1.d0
        ssum = ssum*w + snode
        maxterm = term
        nvalid = nvalid + 1
      ELSE
        w = EXP(term-maxterm)
        sumexp = sumexp + w
        ssum = ssum + w*snode
        nvalid = nvalid + 1
      END IF
    END DO

    IF (nvalid == 0 .OR. sumexp <= 0.d0 .OR. sumexp /= sumexp) GOTO 900
    out = out + agh_tt_fe_logconst(i) + maxterm + LOG(sumexp)
    score = score + ssum/sumexp
    nmescur = nmescur + sm
  END DO

  ok = .TRUE.
900 CONTINUE
  DEALLOCATE(z_eta,u_eta,eta_vec,deta)
END SUBROUTINE eval_eta_full_loglik_score_q1


SUBROUTINE loglik_growth_tt_agh_full_equiv_score_q1(                       &
    b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,              &
    nquad_eta0, nodes_eta0, weights_eta0,                                  &
    nquad_re0, nodes_re0, weights_re0,                                     &
    nitem0, npm0, re_item0, theta_time0, delta_free0,                      &
    monotone_tau0, update_gh0, fit_id0, Q0, loadmask0,                     &
    p_fix0, p_ran0, X0, Z0, delta_mask0, delta_x_col0, delta_z_col0,       &
    vrais, score0, ok0)                                                     &
    BIND(C, name="loglik_growth_tt_agh_full_equiv_score_q1")
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0, ntre0
  INTEGER, INTENT(IN) :: nquad_eta0, nquad_re0, npm0, re_item0, theta_time0
  INTEGER, INTENT(IN) :: delta_free0, monotone_tau0, update_gh0, fit_id0, Q0
  INTEGER, INTENT(IN) :: p_fix0, p_ran0, delta_x_col0, delta_z_col0
  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0), INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0), INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: X0, Z0, delta_mask0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_eta0
  DOUBLE PRECISION, DIMENSION(nquad_eta0), INTENT(IN) :: weights_eta0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: nodes_re0, weights_re0
  DOUBLE PRECISION, INTENT(OUT) :: vrais
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(OUT) :: score0
  INTEGER, INTENT(OUT) :: ok0

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.d9
  DOUBLE PRECISION :: ll
  DOUBLE PRECISION, DIMENSION(npm0) :: score_ll
  INTEGER :: cache_key
  LOGICAL :: ok_model, ok_maps, ok_chol, ok_ll

  vrais = BAD_OBJ
  score0 = 0.d0
  ok0 = 0

  IF (update_gh0 /= 0) GOTO 900

  CALL prepare_growth_full_equiv_score_model_q1(                           &
       b0,Y0,nvar0,nsub0,nobs0,n0,nmes0,c0,cm10,ntre0,                    &
       nquad_eta0,nodes_eta0,weights_eta0,nquad_re0,nodes_re0,weights_re0,&
       nitem0,npm0,re_item0,theta_time0,delta_free0,monotone_tau0,         &
       Q0,loadmask0,p_fix0,p_ran0,X0,Z0,delta_mask0,                       &
       delta_x_col0,delta_z_col0,ok_model)
  IF (.NOT. ok_model) GOTO 900

  IF (Q /= 1 .OR. p_ran /= 2) GOTO 900
  IF (growth_ident_flag /= 1 .AND. growth_ident_flag /= 2 .AND. growth_ident_flag /= 3) GOTO 900

  cache_key = IEOR(fit_id0,INT(Z'40000000'))
  IF (.NOT. agh_tt_fe_grid_ready) GOTO 900
  IF (agh_tt_fe_fit_id_cached /= cache_key) GOTO 900
  IF (agh_tt_fe_dimEta_cached /= n*Q) GOTO 900
  IF (agh_tt_fe_nvar_cached /= nvar) GOTO 900
  IF (agh_tt_fe_nquad_eta_cached /= nquad) GOTO 900
  IF (agh_tt_fe_nquad_re_cached /= nquad_re) GOTO 900
  IF (agh_tt_fe_nsub_cached /= nsub) GOTO 900

  CALL build_score_maps_q1(b0,npm0,ntre0,monotone_tau0,delta_free0,ok_maps)
  IF (.NOT. ok_maps) GOTO 900
  CALL build_chol_score_maps_q1(npm0,ok_chol)
  IF (.NOT. ok_chol) GOTO 900

  CALL eval_eta_full_loglik_score_q1(npm0,ll,score_ll,ok_ll)
  IF (.NOT. ok_ll) GOTO 900
  IF (ll /= ll .OR. ABS(ll) > 1.d300) GOTO 900

  vrais = -ll
  score0 = -score_ll
  ok0 = 1

900 CONTINUE
  CALL clear_score_maps_q1()
  CALL clear_growth_re_model()
END SUBROUTINE loglik_growth_tt_agh_full_equiv_score_q1


SUBROUTINE build_cfa_score_maps_q1(b0, npm0, ntre0, monotone_tau0, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npm0, ntre0, monotone_tau0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: ntau, dimEta, pos, idx, pos_tau, k, m, nf, h, j
  INTEGER :: t, r, row, col, ref_item, i
  DOUBLE PRECISION :: prev

  ok = .FALSE.

  IF (Q /= 1) RETURN
  IF (phi_struct_flag /= 0) RETURN
  IF (tau_invar_flag /= 1) RETURN
  IF (lambda_invar_flag /= 1) RETURN
  IF (.NOT. ALLOCATED(cfa_marker_item)) RETURN
  IF (SIZE(cfa_marker_item) /= Q) RETURN
  IF (.NOT. ALLOCATED(tau) .OR. .NOT. ALLOCATED(tau_long)) RETURN
  IF (.NOT. ALLOCATED(lambda) .OR. .NOT. ALLOCATED(lambda_long)) RETURN
  IF (.NOT. ALLOCATED(gamma) .OR. .NOT. ALLOCATED(gamma_long)) RETURN
  IF (.NOT. ALLOCATED(theta_k) .OR. .NOT. ALLOCATED(theta_long)) RETURN
  IF (.NOT. ALLOCATED(invsd_k) .OR. .NOT. ALLOCATED(invsd_long)) RETURN
  IF (.NOT. ALLOCATED(sigmak) .OR. .NOT. ALLOCATED(sqrt_sigmak)) RETURN
  IF (.NOT. ALLOCATED(kappa) .OR. .NOT. ALLOCATED(L_phi)) RETURN

  IF (ALLOCATED(kappa_sub)) DEALLOCATE(kappa_sub)
  IF (ALLOCATED(L_phi_sub)) DEALLOCATE(L_phi_sub)

  dimEta = n * Q
  ntau = SUM(cm1)

  CALL clear_score_maps_q1()
  sc_npm = npm0

  ALLOCATE(sc_dtau(ntau,npm0))
  ALLOCATE(sc_dlambda(nvar,Q,npm0))
  ALLOCATE(sc_dgamma(nvar,n,npm0))
  ALLOCATE(sc_dtheta(nvar,n,npm0))
  ALLOCATE(sc_dkappa(dimEta,nsub,npm0))
  ALLOCATE(sc_dsqrt_sigmak(nvar,npm0))
  ALLOCATE(sc_dL(dimEta,dimEta,1,npm0))

  sc_dtau = 0.d0
  sc_dlambda = 0.d0
  sc_dgamma = 0.d0
  sc_dtheta = 0.d0
  sc_dkappa = 0.d0
  sc_dsqrt_sigmak = 0.d0
  sc_dL = 0.d0

  tau = 0.d0
  tau_long = 0.d0
  idx = 0
  pos_tau = 0
  DO k = 1, nvar
    m = cm1(k)
    nf = m - 2
    IF (nf < 1) GOTO 900

    SELECT CASE (cfa_threshold_id_flag)

    CASE (0)
      tau_long(pos_tau+2,:) = 0.d0
      tau_long(pos_tau+3,:) = 1.d0

      idx = idx + 1
      IF (idx > ntre0) GOTO 900
      IF (monotone_tau0 == 0) THEN
        tau_long(pos_tau+1,:) = b0(idx)
        sc_dtau(pos_tau+1,idx) = 1.d0
      ELSE
        tau_long(pos_tau+1,:) = -EXP(b0(idx))
        sc_dtau(pos_tau+1,idx) = -EXP(b0(idx))
      END IF

      IF (nf > 1) THEN
        IF (monotone_tau0 == 0) THEN
          DO h = 2, nf
            idx = idx + 1
            IF (idx > ntre0) GOTO 900
            tau_long(pos_tau+h+2,:) = b0(idx)
            sc_dtau(pos_tau+h+2,idx) = 1.d0
          END DO
        ELSE
          prev = 1.d0
          DO h = 2, nf
            idx = idx + 1
            IF (idx > ntre0) GOTO 900
            prev = prev + EXP(b0(idx))
            tau_long(pos_tau+h+2,:) = prev
            DO j = h, nf
              sc_dtau(pos_tau+j+2,idx) = EXP(b0(idx))
            END DO
          END DO
        END IF
      END IF

    CASE (1)
      tau_long(pos_tau+1,:) = 0.d0
      tau_long(pos_tau+2,:) = 1.d0

      IF (monotone_tau0 == 0) THEN
        DO h = 1, nf
          idx = idx + 1
          IF (idx > ntre0) GOTO 900
          tau_long(pos_tau+h+2,:) = b0(idx)
          sc_dtau(pos_tau+h+2,idx) = 1.d0
        END DO
      ELSE
        prev = 1.d0
        DO h = 1, nf
          idx = idx + 1
          IF (idx > ntre0) GOTO 900
          prev = prev + EXP(b0(idx))
          tau_long(pos_tau+h+2,:) = prev
          DO j = h, nf
            sc_dtau(pos_tau+j+2,idx) = EXP(b0(idx))
          END DO
        END DO
      END IF

    CASE DEFAULT
      GOTO 900
    END SELECT

    pos_tau = pos_tau + m
  END DO
  IF (idx /= ntre0) GOTO 900
  tau = tau_long(:,1)

  pos = ntre0

  lambda = 0.d0
  lambda_long = 0.d0
  ref_item = cfa_marker_item(1)
  IF (ref_item < 1 .OR. ref_item > nvar) GOTO 900
  IF (loadmask(ref_item,1) /= 1) GOTO 900
  lambda(ref_item,1) = 1.d0
  DO k = 1, nvar
    IF (loadmask(k,1) == 0 .OR. k == ref_item) CYCLE
    pos = pos + 1
    IF (pos > npm0) GOTO 900
    lambda(k,1) = b0(pos)
    sc_dlambda(k,1,pos) = 1.d0
  END DO
  DO t = 1, n
    lambda_long(:,:,t) = lambda
  END DO

  gamma = 0.d0
  gamma_long = 0.d0
  IF (gamma_invar_flag == 1) THEN
    DO k = 1, nvar
      IF (cfa_mean_ident_flag == 1 .AND. k == ref_item) CYCLE
      pos = pos + 1
      IF (pos > npm0) GOTO 900
      gamma(k) = b0(pos)
      gamma_long(k,:) = gamma(k)
      DO t = 1, n
        sc_dgamma(k,t,pos) = 1.d0
      END DO
    END DO
    IF (ALLOCATED(cfa_gamma_time_free)) THEN
      IF (SIZE(cfa_gamma_time_free,1) == nvar .AND. SIZE(cfa_gamma_time_free,2) == n) THEN
        DO t = 2, n
          DO k = 1, nvar
            IF (cfa_gamma_time_free(k,t) /= 1) CYCLE
            IF (cfa_mean_ident_flag == 1 .AND. k == ref_item) CYCLE
            pos = pos + 1
            IF (pos > npm0) GOTO 900
            gamma_long(k,t) = b0(pos)
            sc_dgamma(k,t,:) = 0.d0
            sc_dgamma(k,t,pos) = 1.d0
          END DO
        END DO
      END IF
    END IF
  ELSE
    DO t = 1, n
      DO k = 1, nvar
        IF (cfa_mean_ident_flag == 1 .AND. k == ref_item) CYCLE
        pos = pos + 1
        IF (pos > npm0) GOTO 900
        gamma_long(k,t) = b0(pos)
        sc_dgamma(k,t,pos) = 1.d0
      END DO
    END DO
    gamma = gamma_long(:,1)
  END IF

  theta_long = 1.d0
  invsd_long = 1.d0
  IF (theta_invar_flag == 1) THEN
    DO k = 1, nvar
      IF (cm1(k) <= 1) THEN
        theta_k(k) = 1.d0
        invsd_k(k) = 1.d0
      ELSE
        pos = pos + 1
        IF (pos > npm0) GOTO 900
        theta_k(k) = EXP(b0(pos))
        invsd_k(k) = 1.d0 / SQRT(theta_k(k))
        theta_long(k,:) = theta_k(k)
        invsd_long(k,:) = invsd_k(k)
        DO t = 1, n
          sc_dtheta(k,t,pos) = theta_long(k,t)
        END DO
      END IF
    END DO
  ELSE
    DO k = 1, nvar
      IF (cm1(k) <= 1) THEN
        theta_k(k) = 1.d0
        invsd_k(k) = 1.d0
        theta_long(k,:) = 1.d0
        invsd_long(k,:) = 1.d0
      ELSE
        DO t = 1, n
          pos = pos + 1
          IF (pos > npm0) GOTO 900
          theta_long(k,t) = EXP(b0(pos))
          invsd_long(k,t) = 1.d0 / SQRT(theta_long(k,t))
          sc_dtheta(k,t,pos) = theta_long(k,t)
        END DO
        theta_k(k) = theta_long(k,1)
        invsd_k(k) = invsd_long(k,1)
      END IF
    END DO
  END IF

  IF (re_item_flag == 1) THEN
    DO k = 1, nvar
      pos = pos + 1
      IF (pos > npm0) GOTO 900
      sigmak(k) = EXP(b0(pos))
      sqrt_sigmak(k) = SQRT(sigmak(k))
      sc_dsqrt_sigmak(k,pos) = 0.5d0 * sqrt_sigmak(k)
    END DO
  ELSE
    sigmak = 0.d0
    sqrt_sigmak = 0.d0
  END IF

  kappa = 0.d0
  IF (cfa_mean_ident_flag == 1) THEN
    DO t = 1, n
      row = t
      pos = pos + 1
      IF (pos > npm0) GOTO 900
      kappa(row) = b0(pos)
      DO i = 1, nsub
        sc_dkappa(row,i,pos) = 1.d0
      END DO
    END DO
  ELSE
    IF (n > 1) THEN
      DO t = 2, n
        row = t
        pos = pos + 1
        IF (pos > npm0) GOTO 900
        kappa(row) = b0(pos)
        DO i = 1, nsub
          sc_dkappa(row,i,pos) = 1.d0
        END DO
      END DO
    END IF
  END IF

  CALL ensure_cfa_beta_storage()
  IF (p_fix > 0) THEN
    IF (.NOT. ALLOCATED(cfa_X)) GOTO 900
    DO r = 1, p_fix
      pos = pos + 1
      IF (pos > npm0) GOTO 900
      cfa_beta(1,r) = b0(pos)
      DO i = 1, nsub
        DO t = 1, n
          row = t
          sc_dkappa(row,i,pos) = cfa_X(i,t,1,r)
        END DO
      END DO
    END DO
  END IF

  L_phi = 0.d0
  DO row = 1, dimEta
    pos = pos + 1
    IF (pos > npm0) GOTO 900
    L_phi(row,row) = EXP(b0(pos))
    sc_dL(row,row,1,pos) = L_phi(row,row)
  END DO

  DO col = 1, dimEta
    DO row = col+1, dimEta
      pos = pos + 1
      IF (pos > npm0) GOTO 900
      L_phi(row,col) = b0(pos)
      sc_dL(row,col,1,pos) = 1.d0
    END DO
  END DO

  IF (pos /= npm0) GOTO 900

  ok = .TRUE.

900 CONTINUE
  IF (.NOT. ok) CALL clear_score_maps_q1()
END SUBROUTINE build_cfa_score_maps_q1


SUBROUTINE eval_cfa_eta_full_loglik_score_q1(npm0, out, score, ok)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npm0
  DOUBLE PRECISION, INTENT(OUT) :: out
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(OUT) :: score
  LOGICAL, INTENT(OUT) :: ok

  INTEGER :: i, l, d, j, k_item, sm, nvalid, dimEta, baseZ
  DOUBLE PRECISION :: q_eta, q_u, logg, logint, term, maxterm, sumexp, w
  DOUBLE PRECISION :: alpha_k, sdev_k
  DOUBLE PRECISION, ALLOCATABLE :: z_eta(:), u_eta(:), eta_vec(:), deta(:,:)
  DOUBLE PRECISION, DIMENSION(npm0) :: snode, sitem, ssum
  LOGICAL :: ok_item

  out = 0.d0
  score = 0.d0
  ok = .FALSE.

  IF (.NOT. agh_tt_fe_grid_ready) RETURN
  IF (sc_npm /= npm0) RETURN
  IF (.NOT. ALLOCATED(sc_dL)) RETURN
  IF (.NOT. ALLOCATED(sc_dkappa)) RETURN

  dimEta = n * Q
  IF (Q /= 1) RETURN
  IF (agh_tt_fe_dimEta_cached /= dimEta) RETURN
  IF (agh_tt_fe_nvar_cached /= nvar) RETURN
  IF (agh_tt_fe_nquad_eta_cached /= nquad) RETURN
  IF (agh_tt_fe_nquad_re_cached /= nquad_re) RETURN
  IF (agh_tt_fe_nsub_cached /= nsub) RETURN

  ALLOCATE(z_eta(dimEta),u_eta(dimEta),eta_vec(dimEta),deta(dimEta,npm0))
  nmescur = 0

  DO i = 1, nsub
    sm = SUM(nmes(i,:))
    nvalid = 0
    maxterm = 0.d0
    sumexp = 0.d0
    ssum = 0.d0

    DO l = 1, nquad
      IF (weights(l) <= 0.d0) CYCLE
      baseZ = (l-1)*dimEta
      DO d = 1, dimEta
        u_eta(d) = nodes(baseZ+d)
      END DO
      z_eta = agh_tt_fe_zeta_eta(:,l,i)

      CALL cfa_eta_mean(i,eta_vec)
      DO d = 1, dimEta
        DO j = 1, d
          eta_vec(d) = eta_vec(d) + L_phi(d,j)*z_eta(j)
        END DO
      END DO

      deta = sc_dkappa(:,i,:)
      DO d = 1, dimEta
        DO j = 1, d
          deta(d,:) = deta(d,:) + sc_dL(d,j,1,:)*z_eta(j)
        END DO
      END DO

      logg = 0.d0
      snode = 0.d0
      DO k_item = 1, nvar
        IF (re_item_flag == 1) THEN
          alpha_k = agh_tt_fe_alpha(k_item,l,i)
          sdev_k = agh_tt_fe_sdev(k_item,i)
        ELSE
          alpha_k = 0.d0
          sdev_k = 1.d0
        END IF

        CALL inneragh_item_logint_full_equiv_score_q1( &
             i,k_item,eta_vec,deta,alpha_k,sdev_k,npm0,logint,sitem,ok_item)
        IF (.NOT. ok_item) THEN
          logg = -1.d9
          EXIT
        END IF
        logg = logg + logint
        snode = snode + sitem
      END DO
      IF (logg <= -1.d9) CYCLE

      q_eta = DOT_PRODUCT(z_eta,z_eta)
      q_u = DOT_PRODUCT(u_eta,u_eta)
      term = LOG(weights(l)) + logg - 0.5d0*(q_eta-q_u)
      IF (term /= term) CYCLE

      IF (nvalid == 0) THEN
        maxterm = term
        sumexp = 1.d0
        ssum = snode
        nvalid = 1
      ELSE IF (term > maxterm) THEN
        w = EXP(maxterm-term)
        sumexp = sumexp*w + 1.d0
        ssum = ssum*w + snode
        maxterm = term
        nvalid = nvalid + 1
      ELSE
        w = EXP(term-maxterm)
        sumexp = sumexp + w
        ssum = ssum + w*snode
        nvalid = nvalid + 1
      END IF
    END DO

    IF (nvalid == 0 .OR. sumexp <= 0.d0 .OR. sumexp /= sumexp) GOTO 900
    out = out + agh_tt_fe_logconst(i) + maxterm + LOG(sumexp)
    score = score + ssum/sumexp
    nmescur = nmescur + sm
  END DO

  ok = .TRUE.

900 CONTINUE
  DEALLOCATE(z_eta,u_eta,eta_vec,deta)
END SUBROUTINE eval_cfa_eta_full_loglik_score_q1


SUBROUTINE loglik_long_tt_agh_full_equiv_score_q1(                        &
    b0, Y0, nvar0, nsub0, nobs0, n0, nmes0, c0, cm10, ntre0,              &
    nquad_eta0, nitem0, nodes_eta0, weights_eta0,                          &
    nquad_re0, nodes_re0, weights_re0, npm0, re_item0, theta_time0,       &
    meas_invar0, Q0, loadmask0, phi_struct0, monotone_tau0,               &
    update_gh0, fit_id0, vrais, score0, ok0)                               &
    BIND(C, name="loglik_long_tt_agh_full_equiv_score_q1")
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nvar0, nsub0, nobs0, n0, ntre0
  INTEGER, INTENT(IN) :: nquad_eta0, nquad_re0, npm0, re_item0, theta_time0
  INTEGER, DIMENSION(4), INTENT(IN) :: meas_invar0
  INTEGER, INTENT(IN) :: Q0, phi_struct0, monotone_tau0, update_gh0, fit_id0
  INTEGER, DIMENSION(nsub0,nvar0), INTENT(IN) :: nmes0
  INTEGER, DIMENSION(nvar0), INTENT(IN) :: c0, cm10
  INTEGER, DIMENSION(nobs0), INTENT(IN) :: Y0, nitem0
  INTEGER, DIMENSION(nvar0,Q0), INTENT(IN) :: loadmask0
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(IN) :: b0
  DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: nodes_eta0
  DOUBLE PRECISION, DIMENSION(nquad_eta0), INTENT(IN) :: weights_eta0
  DOUBLE PRECISION, DIMENSION(nquad_re0), INTENT(IN) :: nodes_re0, weights_re0
  DOUBLE PRECISION, INTENT(OUT) :: vrais
  DOUBLE PRECISION, DIMENSION(npm0), INTENT(OUT) :: score0
  INTEGER, INTENT(OUT) :: ok0

  DOUBLE PRECISION, PARAMETER :: BAD_OBJ = 1.d9
  DOUBLE PRECISION :: ll
  DOUBLE PRECISION, DIMENSION(npm0) :: score_ll
  LOGICAL :: ok_maps, ok_ll

  vrais = BAD_OBJ
  score0 = 0.d0
  ok0 = 0

  IF (update_gh0 /= 0) RETURN
  IF (Q0 /= 1) RETURN
  IF (phi_struct0 /= 0) RETURN
  IF (meas_invar0(1) /= 1 .OR. meas_invar0(3) /= 1) RETURN

  IF (.NOT. agh_tt_fe_grid_ready) RETURN
  IF (agh_tt_fe_fit_id_cached /= fit_id0) RETURN
  IF (agh_tt_fe_dimEta_cached /= n0*Q0) RETURN
  IF (agh_tt_fe_nvar_cached /= nvar0) RETURN
  IF (agh_tt_fe_nquad_eta_cached /= nquad_eta0) RETURN
  IF (agh_tt_fe_nquad_re_cached /= nquad_re0) RETURN
  IF (agh_tt_fe_nsub_cached /= nsub0) RETURN

  CALL build_cfa_score_maps_q1(b0,npm0,ntre0,monotone_tau0,ok_maps)
  IF (.NOT. ok_maps) GOTO 900

  CALL eval_cfa_eta_full_loglik_score_q1(npm0,ll,score_ll,ok_ll)
  IF (.NOT. ok_ll) GOTO 900
  IF (ll /= ll .OR. ABS(ll) > 1.d300) GOTO 900

  vrais = -ll
  score0 = -score_ll
  ok0 = 1

900 CONTINUE
  CALL clear_score_maps_q1()
END SUBROUTINE loglik_long_tt_agh_full_equiv_score_q1



END MODULE growth_re_marginal_agh_score_extension
