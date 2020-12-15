c
c---- SUBROUTINE SEGPLANEINT -----------------------------------------
c
      SUBROUTINE SEGPLANEINT(a,b,c,q,r,p,m,iflag)
C
      IMPLICIT NONE
c
c.....Input and Output arrays.
      REAL    D,DOTPROD
      REAL    A(3),B(3),C(3),N(3),Q(3),R(3),P(3),RQ(3)
      INTEGER I,M,IFLAG
c
c.....Local arrays.
      REAL    NUM, DENOM, T      
c
      CALL PLANECOEFF(A,B,C,N,D,M)
      NUM = D - DOTPROD(Q,N)
      CALL SUBVEC(R,Q,RQ)
      DENOM = DOTPROD(RQ,N)

c      WRITE(6,*) NUM,DENOM
      IF(DENOM==0.0) THEN ! Segment is parallel to plane
        IF(NUM==0.0) THEN ! The segment lies wholly withing the plane.
          IFLAG = 2
          RETURN
        ELSE
          IFLAG = 0
          RETURN  ! The segment lies strictly to one side or the other of the plane.
        ENDIF
      ELSE
        t=num/denom
      ENDIF

      DO I=1,3
        P(I) = Q(I)+T*(R(I)-Q(I))
      ENDDO

      IF ( (t>0.0).AND.(t<1.0) ) THEN
        IFLAG = 1
        RETURN
      ELSEIF (NUM==0.0) THEN !The first q endpoint is on the plane.
        IFLAG = 3
        RETURN
      ELSEIF (NUM==DENOM) THEN !The second r endpoint is on the plane.
        IFLAG = 4
c        WRITE(6,*) 'IFLAG=',IFLAG
        RETURN
      ELSE 
        IFLAG = 0
        RETURN !The segment lies strictly to one side ot the other of the plane.
      ENDIF

      RETURN

      END
c---------------------------------------------------------------------


c
c---- SUBROUTINE PLANECOEFF ------------------------------------------
c
      SUBROUTINE PLANECOEFF(a,b,c,n,d,m)
C
      IMPLICIT NONE
c
      REAL    D
      REAL    A(3),B(3),C(3),N(3)
      INTEGER I,M
c
c.....Local arrays.
      REAL    MAX,DOTPROD
c
      CALL NORMALVEC(A,B,C,N)
      D = DOTPROD(A,N)
c
c.....Find the largest component of the normal
      MAX = 0.0
      DO I=1,3
        IF(ABS(N(I))>MAX) M=I
      ENDDO
      IF(M<1 .OR. M>3) WRITE(0,*) 'WARNING M=',M

      RETURN

      END

c---------------------------------------------------------------------


c
c---- SUBROUTINE NORMALVEC -------------------------------------------

      SUBROUTINE NORMALVEC(a,b,c,n)
c
      IMPLICIT NONE
c
      REAL    A(3),B(3),C(3),N(3)
c
c.....Local arrays
c
      N(1) = (B(2)-A(2))*(C(3)-A(3))
     &      -(B(3)-A(3))*(C(2)-A(2))

      N(2) = (C(1)-A(1))*(B(3)-A(3))
     &      -(C(3)-A(3))*(B(1)-A(1))

      N(3) = (B(1)-A(1))*(C(2)-A(2))
     &      -(B(2)-A(2))*(C(1)-A(1))

      RETURN

      END

c---------------------------------------------------------------------


c---- FUNCTION DOTPROD------------------------------------------------
c
      FUNCTION DOTPROD(a,b)
c
      IMPLICIT NONE
c
      REAL    DOTPROD
      REAL    A(3),B(3)
c
c.....Local arrays
      INTEGER I
c
      DOTPROD = 0.
      DO I=1,3
        DOTPROD = DOTPROD + A(I)*B(I)
      ENDDO

      RETURN

      END

c---------------------------------------------------------------------

c---- SUBROUTINE SUBVEC ----------------------------------------------
c
      SUBROUTINE SUBVEC(a,b,c)
c
      IMPLICIT NONE
c
      REAL    A(3),B(3),C(3)
c
      INTEGER I
c
      DO I=1,3
        C(I) = A(I)-B(I)
      ENDDO

      RETURN

      END

c---------------------------------------------------------------------
c
c
c
c---- SUBROUTINE INTRI3D ----------------------------------------------
c
      SUBROUTINE INTRI3D(a,b,c,p,m,iflag)
C
      IMPLICIT NONE
c
c.....Input and Output arrays.
      REAL    A(3),B(3),C(3),P(3)
      INTEGER M,IFLAG
c
c.....Local arrays.
      INTEGER I,J,K
      REAL    AP(2),BP(2),CP(2),PP(2)
c
      J=1
      DO I=1,3
        IF(I/=M) THEN
          PP(J)=P(I) !Project point.
          AP(J)=A(I)
          BP(J)=B(I)
          CP(J)=C(I)
          J=J+1
        ENDIF
      ENDDO

c      WRITE(6,'(I4,A,8(1X,F14.8))') M,'PROJECTION'
c     &       ,AP(1),AP(2),BP(1),BP(2)
c     &       ,CP(1),CP(2),PP(1),PP(2)
      CALL INTRI2D(AP,BP,CP,PP,IFLAG)
c      WRITE(6,*) 'IFLAG=',IFLAG
     
      RETURN

      END
c
c---------------------------------------------------------------------
c
c
c
c---- SUBROUTINE INTRI2D ----------------------------------------------
c
      SUBROUTINE INTRI2D(a,b,c,p,iflag)
C
      IMPLICIT NONE
c
c.....Input and Output arrays.
      REAL    A(2),B(2),C(2),P(2)
      INTEGER IFLAG
c
c.....Local arrays.
      INTEGER I,J,K
c      INTEGER AREA0,AREA1,AREA2
c      INTEGER AREASIGN
      REAL    AREA0,AREA1,AREA2
      REAL    AREASIGN

      AREA0 = AREASIGN(p,a,b)
      AREA1 = AREASIGN(p,b,c)
      AREA2 = AREASIGN(p,c,a)

c      WRITE(6,*) 'AREA0=',AREA0
c      WRITE(6,*) 'AREA1=',AREA1
c      WRITE(6,*) 'AREA2=',AREA2


      IF ( ((AREA0 == 0) .AND. (AREA1>0) .AND. (AREA2>0)) .OR.
     &     ((AREA1 == 0) .AND. (AREA0>0) .AND. (AREA2>0)) .OR.
     &     ((AREA2 == 0) .AND. (AREA0>0) .AND. (AREA1>0)) ) THEN
        IFLAG = 2
        RETURN
      ENDIF

      IF ( ((AREA0 == 0) .AND. (AREA1<0) .AND. (AREA2>0)) .OR.
     &     ((AREA1 == 0) .AND. (AREA0<0) .AND. (AREA2>0)) .OR.
     &     ((AREA2 == 0) .AND. (AREA0<0) .AND. (AREA1>0)) ) THEN
        IFLAG = 2
        RETURN
      ENDIF

      IF ( ((AREA0 > 0) .AND. (AREA1>0) .AND. (AREA2>0)) .OR.
     &     ((AREA1 < 0) .AND. (AREA0<0) .AND. (AREA2<0)) ) THEN
        IFLAG = 1
        RETURN
      ENDIF

      IF ( (AREA0==0) .AND. (AREA1==0) .AND. (AREA2==0) ) THEN
        IFLAG = -1
        RETURN
      ENDIF

      IF ( ((AREA0 == 0) .AND. (AREA1==0)) .OR.
     &     ((AREA0 == 0) .AND. (AREA2==0)) .OR.
     &     ((AREA1 == 0) .AND. (AREA2==0))  ) THEN
        IFLAG = 3
        RETURN
      ELSE
        IFLAG = 0
        RETURN
      ENDIF

      END
c
c---------------------------------------------------------------------
c
c
c
c---- FUNCTION AREASIGN ----------------------------------------------
c
      FUNCTION AREASIGN(a,b,c)
C
      IMPLICIT NONE
c
c.....Input and Output arrays.
      REAL    AREA
      REAL    AREASIGN
      REAL    A(2),B(2),C(2)

      AREASIGN = 0.5*( (B(1)-A(1))*(C(2)-A(2))
     &                -(C(1)-A(1))*(B(2)-A(2)) )

c      WRITE(6,'(7(1X,F14.8))') A(1),A(2),B(1),B(2),C(1),C(2),AREASIGN

c      IF(0) THEN
      !The area should be an integer.
c      IF(AREA>0.5) THEN
c        AREASIGN=AREA
c        RETURN
c      ELSEIF(AREA<-0.5) THEN
c        AREASIGN=-1
c        RETURN
c      ELSE
c        AREASIGN=0
c        RETURN
c      ENDIF
c      ELSE
c        AREASIGN = SIGN(1.0,AREA)
c        RETURN
c      ENDIF
c       AREA = -0.001
c       AREASIGN = CEILING(AREA)

      RETURN
      
      END
c
c---------------------------------------------------------------------
c
c
c
c---- SUBROUTINE SEGTRIINT1 ------------------------------------------
c
      SUBROUTINE SEGTRIINT1(a,b,c,q,r,p,iflag)
C
      IMPLICIT NONE
c
c.....Input and Output Arrays.
      INTEGER IFLAG
      REAL    A(3),B(3),C(3),Q(3),R(3),P(3)
c
c.....Local Arrays.
      INTEGER M,CODE
c
      CALL SEGPLANEINT(A,B,C,Q,R,P,M,CODE)
c      WRITE(6,'(A,I4,1X,A,3(1X,F14.8))') 'CODE=',CODE,'P=',P(1),P(2),P(3)

      IF(CODE.EQ.3) THEN
        CALL INTRI3D(A,B,C,Q,M,IFLAG)
c        SEGTRIINT1 = IFLAG
        RETURN
      ELSEIF(CODE.EQ.4) THEN
        CALL INTRI3D(A,B,C,R,M,IFLAG)
c        SEGTRIINT1 = IFLAG
        RETURN
      ELSEIF(CODE.EQ.1) THEN
        CALL INTRI3D(A,B,C,P,M,IFLAG)
c        SEGTRIINT1 = IFLAG
c        WRITE(6,*) IFLAG,SEGTRIINT1
        RETURN
      ELSEIF(CODE.EQ.2) THEN
        CALL INTRI3D(A,B,C,Q,M,IFLAG)
        IF(IFLAG>0) THEN
c          SEGTRIINT1 = IFLAG
          RETURN
        ELSE
          CALL INTRI3D(A,B,C,R,M,IFLAG)
          IF(IFLAG>0) THEN
c            SEGTRIINT1 = IFLAG
            RETURN
          ELSE
c            SEGTRIINT1 = 0
          ENDIF
        ENDIF
      ELSEIF(CODE.EQ.0) THEN
        IFLAG = 0         
c        SEGTRIINT1 = 0
      ENDIF

      END
   
