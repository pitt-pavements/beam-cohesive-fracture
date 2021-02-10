	  subroutine syminv ( a, n, c, w, nullty, ifault )

!*****************************************************************************80
!
!! SYMINV computes the inverse of a symmetric matrix.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    FORTRAN77 version by Michael Healy
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Michael Healy,
!    Algorithm AS 7:
!    Inversion of a Positive Semi-Definite Symmetric Matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 198-199.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A((N*(N+1))/2), a positive definite matrix stored
!    by rows in lower triangular form as a one dimensional array, in the 
!    sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Output, real ( kind = 8 ) C((N*(N+1))/2), the inverse of A, or generalized
!    inverse if A is singular, stored using the same storage scheme employed
!    for A.  The program is written in such a way that A and U can share 
!    storage.
!
!    Workspace, real ( kind = 8 ) W(N).
!
!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.  If NULLTY 
!    is zero, the matrix is judged to have full rank.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no error detected.
!    1, N < 1.
!    2, A is not positive semi-definite.
!
	  implicit none

	  integer ( kind = 4 ) n

	  real ( kind = 8 ) a((n*(n+1))/2)
	  real ( kind = 8 ) c((n*(n+1))/2)
	  integer ( kind = 4 ) i
	  integer ( kind = 4 ) icol
	  integer ( kind = 4 ) ifault
	  integer ( kind = 4 ) irow
	  integer ( kind = 4 ) j
	  integer ( kind = 4 ) jcol
	  integer ( kind = 4 ) k
	  integer ( kind = 4 ) l
	  integer ( kind = 4 ) mdiag
	  integer ( kind = 4 ) ndiag
	  integer ( kind = 4 ) nn
	  integer ( kind = 4 ) nrow
	  integer ( kind = 4 ) nullty
	  real ( kind = 8 ) w(n)
	  real ( kind = 8 ) x

	  ifault = 0

	  if ( n <= 0 ) then
	    ifault = 1
	    return
	  end if

	  nrow = n
	!
	!  Compute the Cholesky factorization of A.
	!  The result is stored in C.
	!
	  nn = ( n* ( n + 1 ) ) / 2

	  call cholesky ( a, n, nn, c, nullty, ifault )

	  if ( ifault /= 0 ) then
	    return
	  end if
	!
	!  Invert C and form the product (Cinv)'* Cinv, where Cinv is the inverse
	!  of C, row by row starting with the last row.
	!  IROW = the row number,
	!  NDIAG = location of last element in the row.
	!
	  irow = nrow
	  ndiag = nn

	  do
	!
	!  Special case, zero diagonal element.
	!
	    if ( c(ndiag) == 0.0D+00 ) then

		l = ndiag
		do j = irow, nrow
		  c(l) = 0.0D+00
		  l = l + j
		end do

	    else

		l = ndiag
		do i = irow, nrow
		  w(i) = c(l)
		  l = l + i
		end do

		icol = nrow
		jcol = nn
		mdiag = nn

		do

		  l = jcol

		  if ( icol == irow ) then
			x = 1.0D+00 / w(irow)
		  else
			x = 0.0D+00
		  end if

		  k = nrow

		  do while ( irow < k )

			x = x - w(k)* c(l)
			k = k - 1
			l = l - 1

			if ( mdiag < l ) then
			  l = l - k + 1
			end if

		  end do

		  c(l) = x / w(irow)

		  if ( icol <= irow ) then
			exit
		  end if

		  mdiag = mdiag - icol
		  icol = icol - 1
		  jcol = jcol - 1

		end do

	    end if

	    ndiag = ndiag - irow
	    irow = irow - 1

	    if ( irow <= 0 ) then
		exit
	    end if

	  end do

	  return
		end
		subroutine timestamp ( )

	!*****************************************************************************80
	!
	!! TIMESTAMP prints the current YMDHMS date as a time stamp.
	!
	!  Example:
	!
	!    31 May 2001   9:45:54.872 AM
	!
	!  Licensing:
	!
	!    This code is distributed under the GNU LGPL license.
	!
	!  Modified:
	!
	!    18 May 2013
	!
	!  Author:
	!
	!    John Burkardt
	!
	!  Parameters:
	!
	!    None
	!
	  implicit none

	  character ( len = 8 ) ampm
	  integer ( kind = 4 ) d
	  integer ( kind = 4 ) h
	  integer ( kind = 4 ) m
	  integer ( kind = 4 ) mm
	  character ( len = 9 ), parameter, dimension(12) :: month = 
	1(/ 
	1    'January  ', 'February ', 'March    ', 'April    ', 
	1    'May      ', 'June     ', 'July     ', 'August   ', 
	1    'September', 'October  ', 'November ', 'December ' /)
	  integer ( kind = 4 ) n
	  integer ( kind = 4 ) s
	  integer ( kind = 4 ) values(8)
	  integer ( kind = 4 ) y

	  call date_and_time ( values = values )

	  y = values(1)
	  m = values(2)
	  d = values(3)
	  h = values(5)
	  n = values(6)
	  s = values(7)
	  mm = values(8)

	  if ( h < 12 ) then
	    ampm = 'AM'
	  else if ( h == 12 ) then
	    if ( n == 0 .and. s == 0 ) then
		ampm = 'Noon'
	    else
		ampm = 'PM'
	    end if
	  else
	    h = h - 12
	    if ( h < 12 ) then
		ampm = 'PM'
	    else if ( h == 12 ) then
		if ( n == 0 .and. s == 0 ) then
		  ampm = 'Midnight'
		else
		  ampm = 'AM'
		end if
	    end if
	  end if

C		write (*,'(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') 
C	 1 	d, trim ( month(m) ), y, h, ':', n, 
C     1 ':', s, '.', mm, trim ( ampm )

	  return
		end
		subroutine cholesky ( a, n, nn, u, nullty, ifault )

	!*****************************************************************************80
	!
	!! CHOLESKY computes the Cholesky factorization of a PDS matrix.
	!
	!  Discussion:
	!
	!    For a positive definite symmetric matrix A, the Cholesky factor U
	!    is an upper triangular matrix such that A = U'* U.
	!
	!    This routine was originally named "CHOL", but that conflicted with
	!    a built in MATLAB routine name.
	!
	!    The missing initialization "II = 0" has been added to the code.
	!
	!  Modified:
	!
	!    11 February 2008
	!
	!  Author:
	!
	!    Original FORTRAN77 version by Michael Healy.
	!    Modifications by AJ Miller.
	!    FORTRAN90 version by John Burkardt.
	!
	!  Reference:
	!
	!    Michael Healy,
	!    Algorithm AS 6:
	!    Triangular decomposition of a symmetric matrix,
	!    Applied Statistics,
	!    Volume 17, Number 2, 1968, pages 195-197.
	!
	!  Parameters:
	!
	!    Input, real ( kind = 8 ) A((N*(N+1))/2), a positive definite matrix
	!    stored by rows in lower triangular form as a one dimensional array,
	!    in the sequence
	!    A(1,1),
	!    A(2,1), A(2,2),
	!    A(3,1), A(3,2), A(3,3), and so on.
	!
	!    Input, integer ( kind = 4 ) N, the order of A.
	!
	!    Input, integer ( kind = 4 ) NN, the dimension of A, (N*(N+1))/2.
	!
	!    Output, real ( kind = 8 ) U((N*(N+1))/2), an upper triangular matrix,
	!    stored by columns, which is the Cholesky factor of A.  The program is
	!    written in such a way that A and U can share storage.
	!
	!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.
	!    If NULLTY is zero, the matrix is judged to have full rank.
	!
	!    Output, integer ( kind = 4 ) IFAULT, an error indicator.
	!    0, no error was detected;
	!    1, if N < 1;
	!    2, if A is not positive semi-definite.
	!    3, if NN < (N*(N+1))/2.
	!
	!  Local Parameters:
	!
	!    Local, real ( kind = 8 ) ETA, should be set equal to the smallest positive
	!    value such that 1.0 + ETA is calculated as being greater than 1.0 in the
	!    accuracy being used.
	!
	  implicit none

	  integer ( kind = 4 ) nn

	  real ( kind = 8 ) a(nn)
	  real ( kind = 8 ), parameter :: eta = 1.0D-09
	  integer ( kind = 4 ) i
	  integer ( kind = 4 ) icol
	  integer ( kind = 4 ) ifault
	  integer ( kind = 4 ) ii
	  integer ( kind = 4 ) irow
	  integer ( kind = 4 ) j
	  integer ( kind = 4 ) k
	  integer ( kind = 4 ) kk
	  integer ( kind = 4 ) l
	  integer ( kind = 4 ) m
	  integer ( kind = 4 ) n
	  integer ( kind = 4 ) nullty
	  real ( kind = 8 ) u(nn)
	  real ( kind = 8 ) w
	  real ( kind = 8 ) x

	  ifault = 0
	  nullty = 0

	  if ( n <= 0 ) then
	    ifault = 1
	    return
	  end if

	  if ( nn < ( n* ( n + 1 ) ) / 2 ) then
	    ifault = 3
	    return
	  end if

	  j = 1
	  k = 0
	  ii = 0
	!
	!  Factorize column by column, ICOL = column number.
	!
	  do icol = 1, n

	    ii = ii + icol
	    x = eta* eta* a(ii)
	    l = 0
	    kk = 0
	!
	!  IROW = row number within column ICOL.
	!
	    do irow = 1, icol

		kk = kk + irow
		k = k + 1
		w = a(k)
		m = j

		do i = 1, irow - 1
		  l = l + 1
		  w = w - u(l)* u(m)
		  m = m + 1
		end do

		l = l + 1

		if ( irow == icol ) then
		  exit
		end if

		if ( u(l) /= 0.0D+00 ) then

		  u(k) = w / u(l)

		else

		  u(k) = 0.0D+00

		  if ( abs ( x* a(k) ) < w* w ) then
			ifault = 2
			return
		  end if

		end if

	    end do
	!
	!  End of row, estimate relative accuracy of diagonal element.
	!
	    if ( abs ( w ) <= abs ( eta* a(k) ) ) then

		u(k) = 0.0D+00
		nullty = nullty + 1

	    else

		if ( w < 0.0D+00 ) then
		  ifault = 2
		  return
		end if

		u(k) = sqrt ( w )

	    end if

	    j = j + icol

	  end do

	  return
		end
c
c
c
      SUBROUTINE dsymv(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
c
c  -- Reference BLAS level2 routine (version 3.7.0) --
c  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     December 2016
c
c     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA,BETA
       INTEGER INCX,INCY,LDA,N
       CHARACTER UPLO
c     ..
c     .. Array Arguments ..
       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
c     ..

c  =====================================================================
c
c     .. Parameters ..
       DOUBLE PRECISION ONE,ZERO
       parameter(one=1.0d+0,zero=0.0d+0)
C     ..
C     .. Local Scalars ..
       DOUBLE PRECISION TEMP1,TEMP2
       INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
C     ..
C     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
C     ..
C     .. External Subroutines ..
       EXTERNAL xerbla
C     ..
C     .. Intrinsic Functions ..
       INTRINSIC max
C     ..
C
C     Test the input parameters.
C
       info = 0
       IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
           info = 1
       ELSE IF (n.LT.0) THEN
           info = 2
       ELSE IF (lda.LT.max(1,n)) THEN
           info = 5
       ELSE IF (incx.EQ.0) THEN
           info = 7
       ELSE IF (incy.EQ.0) THEN
           info = 10
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('DSYMV ',info)
           RETURN
       END IF
C
C     Quick return if possible.
C
       IF ((n.EQ.0) .OR. ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
C
C     Set up the start points in  X  and  Y.
C
       IF (incx.GT.0) THEN
           kx = 1
       ELSE
           kx = 1 - (n-1)*incx
       END IF
       IF (incy.GT.0) THEN
           ky = 1
       ELSE
           ky = 1 - (n-1)*incy
       END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
C     First form  y := beta*y.
C
       IF (beta.NE.one) THEN
           IF (incy.EQ.1) THEN
               IF (beta.EQ.zero) THEN
                   DO 10 i = 1,n
                       y(i) = zero
   10             CONTINUE
               ELSE
                   DO 20 i = 1,n
                       y(i) = beta*y(i)
   20             CONTINUE
               END IF
           ELSE
               iy = ky
               IF (beta.EQ.zero) THEN
                   DO 30 i = 1,n
                       y(iy) = zero
                       iy = iy + incy
   30             CONTINUE
               ELSE
                   DO 40 i = 1,n
                       y(iy) = beta*y(iy)
                       iy = iy + incy
   40             CONTINUE
               END IF
          END IF
       END IF
       IF (alpha.EQ.zero) RETURN
       IF (lsame(uplo,'U')) THEN
C
C        Form  y  when A is stored in upper triangle.
C
           IF ((incx.EQ.1) .AND. (incy.EQ.1)) THEN
               DO 60 j = 1,n
                   temp1 = alpha*x(j)
                   temp2 = zero
                   DO 50 i = 1,j - 1
                       y(i) = y(i) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(i)
   50             CONTINUE
                   y(j) = y(j) + temp1*a(j,j) + alpha*temp2
   60         CONTINUE
           ELSE
               jx = kx
               jy = ky
               DO 80 j = 1,n
                   temp1 = alpha*x(jx)
                   temp2 = zero
                   ix = kx
                   iy = ky
                   DO 70 i = 1,j - 1
                       y(iy) = y(iy) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(ix)
                       ix = ix + incx
                       iy = iy + incy
   70             CONTINUE
                   y(jy) = y(jy) + temp1*a(j,j) + alpha*temp2
                   jx = jx + incx
                   jy = jy + incy
   80         CONTINUE
           END IF
       ELSE
C
C        Form  y  when A is stored in lower triangle.
C
           IF ((incx.EQ.1) .AND. (incy.EQ.1)) THEN
               DO 100 j = 1,n
                   temp1 = alpha*x(j)
                   temp2 = zero
                   y(j) = y(j) + temp1*a(j,j)
                   DO 90 i = j + 1,n
                       y(i) = y(i) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(i)
   90             CONTINUE
                   y(j) = y(j) + alpha*temp2
  100         CONTINUE
           ELSE
               jx = kx
               jy = ky
               DO 120 j = 1,n
                   temp1 = alpha*x(jx)
                   temp2 = zero
                   y(jy) = y(jy) + temp1*a(j,j)
                   ix = jx
                   iy = jy
                   DO 110 i = j + 1,n
                       ix = ix + incx
                       iy = iy + incy
                       y(iy) = y(iy) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(ix)
  110             CONTINUE
                   y(jy) = y(jy) + alpha*temp2
                   jx = jx + incx
                   jy = jy + incy
  120         CONTINUE
           END IF
       END IF
C
       RETURN
C
C     End of DSYMV .
C
      END
