*> \brief \b DGECON
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DGECON + dependencies
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgecon.f">
*> [TGZ]</a>
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgecon.f">
*> [ZIP]</a>
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgecon.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK,
*                          INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          NORM
*       INTEGER            INFO, LDA, N
*       DOUBLE PRECISION   ANORM, RCOND
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   A( LDA, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGECON estimates the reciprocal of the condition number of a general
*> real matrix A, in either the 1-norm or the infinity-norm, using
*> the LU factorization computed by DGETRF.
*>
*> An estimate is obtained for norm(inv(A)), and the reciprocal of the
*> condition number is computed as
*>    RCOND = 1 / ( norm(A) * norm(inv(A)) ).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NORM
*> \verbatim
*>          NORM is CHARACTER*1
*>          Specifies whether the 1-norm condition number or the
*>          infinity-norm condition number is required:
*>          = '1' or 'O':  1-norm;
*>          = 'I':         Infinity-norm.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The factors L and U from the factorization A = P*L*U
*>          as computed by DGETRF.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] ANORM
*> \verbatim
*>          ANORM is DOUBLE PRECISION
*>          If NORM = '1' or 'O', the 1-norm of the original matrix A.
*>          If NORM = 'I', the infinity-norm of the original matrix A.
*> \endverbatim
*>
*> \param[out] RCOND
*> \verbatim
*>          RCOND is DOUBLE PRECISION
*>          The reciprocal of the condition number of the matrix A,
*>          computed as RCOND = 1/(norm(A) * norm(inv(A))).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (4*N)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (N)
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
*> \ingroup doubleGEcomputational
*
*  =====================================================================
      SUBROUTINE dgecon( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK,
     $                   INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,
*  --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG
*  Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   ANORM, RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ONENRM
      CHARACTER          NORMIN
      INTEGER            IX, KASE, KASE1
      DOUBLE PRECISION   AINVNM, SCALE, SL, SMLNUM, SU
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           lsame, idamax, dlamch
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlacn2, dlatrs, drscl, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      onenrm = norm.EQ.'1' .OR. lsame( norm, 'O' )
      IF( .NOT.onenrm .AND. .NOT.lsame( norm, 'I' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( anorm.LT.zero ) THEN
         info = -5
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGECON', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      rcond = zero
      IF( n.EQ.0 ) THEN
         rcond = one
         RETURN
      ELSE IF( anorm.EQ.zero ) THEN
         RETURN
      END IF
*
      smlnum = dlamch( 'Safe minimum' )
*
*     Estimate the norm of inv(A).
*
      ainvnm = zero
      normin = 'N'
      IF( onenrm ) THEN
         kase1 = 1
      ELSE
         kase1 = 2
      END IF
      kase = 0
   10 CONTINUE
      CALL dlacn2( n, work( n+1 ), work, iwork, ainvnm, kase, isave )
      IF( kase.NE.0 ) THEN
         IF( kase.EQ.kase1 ) THEN
*
*           Multiply by inv(L).
*
            CALL dlatrs( 'Lower', 'No transpose', 'Unit', normin, n, a,
     $                   lda, work, sl, work( 2*n+1 ), info )
*
*           Multiply by inv(U).
*
            CALL dlatrs( 'Upper', 'No transpose', 'Non-unit', normin, n,
     $                   a, lda, work, su, work( 3*n+1 ), info )
         ELSE
*
*           Multiply by inv(U**T).
*
            CALL dlatrs( 'Upper', 'Transpose', 'Non-unit', normin, n, a,
     $                   lda, work, su, work( 3*n+1 ), info )
*
*           Multiply by inv(L**T).
*
            CALL dlatrs( 'Lower', 'Transpose', 'Unit', normin, n, a,
     $                   lda, work, sl, work( 2*n+1 ), info )
         END IF
*
*        Divide X by 1/(SL*SU) if doing so will not cause overflow.
*
         scale = sl*su
         normin = 'Y'
         IF( scale.NE.one ) THEN
            ix = idamax( n, work, 1 )
            IF( scale.LT.abs( work( ix ) )*smlnum .OR. scale.EQ.zero )
     $         GO TO 20
            CALL drscl( n, scale, work, 1 )
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( ainvnm.NE.zero )
     $   rcond = ( one / ainvnm ) / anorm
*
   20 CONTINUE
      RETURN
*
*     End of DGECON
*
      END
