*> \brief \b DPOTRS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DPOTRS + dependencies
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpotrs.f">
*> [TGZ]</a>
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpotrs.f">
*> [ZIP]</a>
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpotrs.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DPOTRS solves a system of linear equations A*X = B with a symmetric
*> positive definite matrix A using the Cholesky factorization
*> A = U**T*U or A = L*L**T computed by DPOTRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The triangular factor U or L from the Cholesky
*factorization
*>          A = U**T*U or A = L*L**T, as computed by DPOTRF.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*>          On entry, the right hand side matrix B.
*>          On exit, the solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doublePOcomputational
*
*  =====================================================================
      SUBROUTINE dpotrs( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,
*  --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG
*  Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           dtrsm, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      upper = lsame( uplo, 'U' )
      IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( nrhs.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( ldb.LT.max( 1, n ) ) THEN
         info = -7
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DPOTRS', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 .OR. nrhs.EQ.0 )
     $   RETURN
*
      IF( upper ) THEN
*
*        Solve A*X = B where A = U**T *U.
*
*        Solve U**T *X = B, overwriting B with X.
*
         CALL dtrsm( 'Left', 'Upper', 'Transpose', 'Non-unit', n, nrhs,
     $               one, a, lda, b, ldb )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL dtrsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n,
     $               nrhs, one, a, lda, b, ldb )
      ELSE
*
*        Solve A*X = B where A = L*L**T.
*
*        Solve L*X = B, overwriting B with X.
*
         CALL dtrsm( 'Left', 'Lower', 'No transpose', 'Non-unit', n,
     $               nrhs, one, a, lda, b, ldb )
*
*        Solve L**T *X = B, overwriting B with X.
*
         CALL dtrsm( 'Left', 'Lower', 'Transpose', 'Non-unit', n, nrhs,
     $               one, a, lda, b, ldb )
      END IF
*
      RETURN
*
*     End of DPOTRS
*
      END
