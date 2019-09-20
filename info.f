C
C 
C-----SUBROUTINE-SCRINFO------------------------------------------------
C
      SUBROUTINE SCRINFO(u,v,w,p,divv,tv,dens,dtm,time,icycle,id,jd,kd,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - CHECK FOR MAX. AND MIN. VELOCITY IN THE FLOW FIELD
C                                                       
C-----------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
C 
C-----------------------------------------------------------------------
C
      integer nx,ny,nz
      real    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), p(nx,ny,nz)
      real    divv(nx,ny,nz), tv(nx,ny,nz), dens(nx,ny,nz)
      real    dtm, time
      INTEGER icycle, id, jd, kd
C
      REAL,    DIMENSION(:,:), ALLOCATABLE :: MINL,MING
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: INDL,INDG
      REAL    valin(14),MINV(14),MINRANK(14)
      INTEGER i,j,k,ind(3,14),jndi,INDC(3,14)
c
c.....calculate divergence      
c
      DIVV=0.0
      CALL DIVERGENCE(U,V,W,DIVV,NX,NY,NZ)
C 
C.....check for maximum and minimum values
C      
      IND(:, 1)=MINLOC(U(IX1-1:IX2,JY1:JY2,KZ1:KZ2))
      IND(:, 2)=MAXLOC(U(IX1-1:IX2,JY1:JY2,KZ1:KZ2))
      IND(:, 3)=MINLOC(V(IX1:IX2,JY1-1:JY2,KZ1:KZ2))
      IND(:, 4)=MAXLOC(V(IX1:IX2,JY1-1:JY2,KZ1:KZ2))
      IND(:, 5)=MINLOC(W(IX1:IX2,JY1:JY2,KZ1-1:KZ2))
      IND(:, 6)=MAXLOC(W(IX1:IX2,JY1:JY2,KZ1-1:KZ2))
      IND(:, 7)=MINLOC(P   (IX1:IX2,JY1:JY2,KZ1:KZ2))
      IND(:, 8)=MAXLOC(P   (IX1:IX2,JY1:JY2,KZ1:KZ2))
      IND(:, 9)=MINLOC(DIVV(IX1:IX2,JY1:JY2,KZ1:KZ2))
      IND(:,10)=MAXLOC(DIVV(IX1:IX2,JY1:JY2,KZ1:KZ2))
      IND(:,13)=MINLOC(DENS(IX1:IX2,JY1:JY2,KZ1:KZ2))
      IND(:,14)=MAXLOC(DENS(IX1:IX2,JY1:JY2,KZ1:KZ2))
C
      IF(ISGS/=0) THEN
        IND(:,11)=MINLOC(TV(IX1:IX2,JY1:JY2,KZ1:KZ2))
        IND(:,12)=MAXLOC(TV(IX1:IX2,JY1:JY2,KZ1:KZ2))
      ENDIF
C
C.....INDICES SHIFT
C
      IND(:,:)=IND(:,:)+1
      IND(1,1:2)=IND(1,1:2)-1
      IND(2,3:4)=IND(2,3:4)-1
      IND(3,5:6)=IND(3,5:6)-1
C
      VALIN( 1)= U   (IND(1, 1),IND(2, 1),IND(3, 1))
      VALIN( 2)=-U   (IND(1, 2),IND(2, 2),IND(3, 2))
      VALIN( 3)= V   (IND(1, 3),IND(2, 3),IND(3, 3))
      VALIN( 4)=-V   (IND(1, 4),IND(2, 4),IND(3, 4))
      VALIN( 5)= W   (IND(1, 5),IND(2, 5),IND(3, 5))
      VALIN( 6)=-W   (IND(1, 6),IND(2, 6),IND(3, 6))
      VALIN( 7)= P   (IND(1, 7),IND(2, 7),IND(3, 7))
      VALIN( 8)=-P   (IND(1, 8),IND(2, 8),IND(3, 8))
      VALIN( 9)= DIVV(IND(1, 9),IND(2, 9),IND(3, 9))
      VALIN(10)=-DIVV(IND(1,10),IND(2,10),IND(3,10))
      VALIN(13)= DENS(IND(1,13),IND(2,13),IND(3,13))
      VALIN(14)=-DENS(IND(1,14),IND(2,14),IND(3,14))
C
      IF(ISGS/=0) THEN
        VALIN(11)= TV  (IND(1,11),IND(2,11),IND(3,11))
        VALIN(12)=-TV  (IND(1,12),IND(2,12),IND(3,12))
      ENDIF
C
      ALLOCATE(MINL(14,MYSIZE),INDL(14*3,MYSIZE))
      ALLOCATE(MING(14,MYSIZE),INDG(14*3,MYSIZE))
      MINL = 0.
      INDL = 0.
      DO I=1,10
        MINL(I,MYRANK+1)=VALIN(I)
        K = (I-1)*3
        INDL(1+K,MYRANK+1)=IND(1,I)
        INDL(2+K,MYRANK+1)=IND(2,I)
        INDL(3+K,MYRANK+1)=IND(3,I)+MYRANK*(NZ-2)
      ENDDO
      IF(ISGS/=0) THEN
        DO I=11,12
          MINL(I,MYRANK+1)=VALIN(I)
          K = (I-1)*3
          INDL(1+K,MYRANK+1)=IND(1,I)
          INDL(2+K,MYRANK+1)=IND(2,I)
          INDL(3+K,MYRANK+1)=IND(3,I)+MYRANK*(NZ-2)
        ENDDO
      ENDIF
      DO I=13,14
        MINL(I,MYRANK+1)=VALIN(I)
        K = (I-1)*3
        INDL(1+K,MYRANK+1)=IND(1,I)
        INDL(2+K,MYRANK+1)=IND(2,I)
        INDL(3+K,MYRANK+1)=IND(3,I)+MYRANK*(NZ-2)
      ENDDO

      CALL MPI_REDUCE(MINL,MING,MYSIZE*14,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(INDL,INDG,MYSIZE*3*14,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
c
c.....only root writes to screen
c      
      IF(MYRANK==0) THEN
c the 0 processor establishes for each variable the rank corresponding
c to the minimum value and then generates the vectors MINV and INDC
c
        MINRANK = MINLOC(MING,DIM=2)
        DO I=1,14
          J = MINRANK(I)
          MINV(I) = MING(I,J)
          K = (I-1)*3
          INDC(1,I) = INDG(1+K,J)
          INDC(2,I) = INDG(2+K,J)
          INDC(3,I) = INDG(3+K,J)
        ENDDO
c
c.....write screen information
c.....the minimum and maximum values are written on screen
c.....for each variable
c
        write(6,'(79("*"))')
        write(6,'(A,I10)'   ) '*..Iteration  = ',icycle
        write(6,'(A,G24.16)') '*..Time Step  = ',dtm
        write(6,'(A,G24.16)') '*..Total Time = ',time           

        WRITE(6,100) '*..u   =', 
     &        MINV(1),' @(',INDC(1,1),INDC(2,1),INDC(3,1),');',
     &       -MINV(2),' @(',INDC(1,2),INDC(2,2),INDC(3,2),')'
        WRITE(6,100) '*..v   =',
     &        MINV(3),' @(',INDC(1,3),INDC(2,3),INDC(3,3),');',
     &       -MINV(4),' @(',INDC(1,4),INDC(2,4),INDC(3,4),')'
        WRITE(6,100) '*..w   =',
     &        MINV(5),' @(',INDC(1,5),INDC(2,5),INDC(3,5),');',
     &       -MINV(6),' @(',INDC(1,6),INDC(2,6),INDC(3,6),')'
        WRITE(6,100) '*..p   =',
     &        MINV(7),' @(',INDC(1,7),INDC(2,7),INDC(3,7),');',
     &       -MINV(8),' @(',INDC(1,8),INDC(2,8),INDC(3,8),')'
        IF(ISGS/=0) THEN
          WRITE(6,100) '*..tv  =',
     &        MINV(11),' @(',INDC(1,11),INDC(2,11),INDC(3,11),');',
     &       -MINV(12),' @(',INDC(1,12),INDC(2,12),INDC(3,12),')'
        ENDIF
        WRITE(6,100) '*..divv=',
     &        MINV( 9),' @(',INDC(1, 9),INDC(2, 9),INDC(3, 9),');',
     &       -MINV(10),' @(',INDC(1,10),INDC(2,10),INDC(3,10),')'
        WRITE(6,100) '*..dens=',
     &        MINV(13),' @(',INDC(1,13),INDC(2,13),INDC(3,13),');',
     &       -MINV(14),' @(',INDC(1,14),INDC(2,14),INDC(3,14),')'
C
      ENDIF
      DEALLOCATE(MINL,MING,INDL,INDG)
C
 100  FORMAT(A,G19.12,A,3I5,A,G19.12,A,3I5,A)
      RETURN 
      END 
C
C 
C-----SUBROUTINE-ENDJOB-------------------------------------------------
C
      SUBROUTINE ENDJOB(ROUTINE,ERRMSG,NUM)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - OUTPUT ERROR MESSAGE AND STOP THE RUN
C                                                       
C-----------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
C 
C-----------------------------------------------------------------------
C
      INTEGER NUM
      CHARACTER ROUTINE*20,ERRMSG*70
c
      IF(MYRANK==0) THEN
        WRITE(6,*) ' IN SUBROUTINE ',TRIM(ROUTINE)
        WRITE(6,*) TRIM(ERRMSG),NUM
        WRITE(6,*) ' ABORTING NOW......'
      ENDIF
c
      CALL MPI_FINALIZE(IERR)
      STOP

      RETURN
      END
C
C 
C-----SUBROUTINE-SCRINFO1-----------------------------------------------
C
      SUBROUTINE SCRINFO1(u,v,w,p,tv,dtm,time,icycle,id,jd,kd,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - CHECK FOR MAX. AND MIN. VELOCITY IN THE FLOW FIELD
C                                                       
C-----------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
C 
C-----------------------------------------------------------------------
C
      integer nx,ny,nz
      real    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), p(nx,ny,nz)
      real    divv(nx,ny,nz), tv(nx,ny,nz)
      real    dtm, time
      INTEGER icycle,id,jd,kd
C
      REAL,    DIMENSION(:,:), ALLOCATABLE :: MINL,MING
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: INDL,INDG
      REAL    valin(12),MINV(12),MINRANK(12)
      INTEGER i,j,k,ind(3,12),jndi,INDC(3,12)
c
c.....calculate divergence      
c
      DIVV = 0.0
      CALL DIVERGENCE(U,V,W,DIVV,NX,NY,NZ)
C 
C.....check for maximum and minimum values
C      
      IND(:, 1)=MINLOC(U(1:NX,1:NY,1:NZ))
      IND(:, 2)=MAXLOC(U(1:NX,1:NY,1:NZ))
      IND(:, 3)=MINLOC(V(1:NX,1:NY,1:NZ))
      IND(:, 4)=MAXLOC(V(1:NX,1:NY,1:NZ))
      IND(:, 5)=MINLOC(W(1:NX,1:NY,1:NZ))
      IND(:, 6)=MAXLOC(W(1:NX,1:NY,1:NZ))
      IND(:, 7)=MINLOC(P   (1:NX,1:NY,1:NZ))
      IND(:, 8)=MAXLOC(P   (1:NX,1:NY,1:NZ))
      IND(:, 9)=MINLOC(DIVV(1:NX,1:NY,1:NZ))
      IND(:,10)=MAXLOC(DIVV(1:NX,1:NY,1:NZ))
C
      IF(ISGS/=0) THEN
        IND(:,11)=MINLOC(TV(1:NX,1:NY,1:NZ))
        IND(:,12)=MAXLOC(TV(1:NX,1:NY,1:NZ))
      ENDIF
C
C.....INDICES SHIFT
C
      IND(:,:)=IND(:,:)
c      IND(1,1:2)=IND(1,1:2)-1
c      IND(2,3:4)=IND(2,3:4)-1
c      IND(3,5:6)=IND(3,5:6)-1
C
      VALIN( 1)= U   (IND(1, 1),IND(2, 1),IND(3, 1))
      VALIN( 2)=-U   (IND(1, 2),IND(2, 2),IND(3, 2))
      VALIN( 3)= V   (IND(1, 3),IND(2, 3),IND(3, 3))
      VALIN( 4)=-V   (IND(1, 4),IND(2, 4),IND(3, 4))
      VALIN( 5)= W   (IND(1, 5),IND(2, 5),IND(3, 5))
      VALIN( 6)=-W   (IND(1, 6),IND(2, 6),IND(3, 6))
      VALIN( 7)= P   (IND(1, 7),IND(2, 7),IND(3, 7))
      VALIN( 8)=-P   (IND(1, 8),IND(2, 8),IND(3, 8))
      VALIN( 9)= DIVV(IND(1, 9),IND(2, 9),IND(3, 9))
      VALIN(10)=-DIVV(IND(1,10),IND(2,10),IND(3,10))
C
      IF(ISGS/=0) THEN
        VALIN(11)= TV  (IND(1,11),IND(2,11),IND(3,11))
        VALIN(12)=-TV  (IND(1,12),IND(2,12),IND(3,12))
      ENDIF
C
      ALLOCATE(MINL(12,MYSIZE),INDL(12*3,MYSIZE))
      ALLOCATE(MING(12,MYSIZE),INDG(12*3,MYSIZE))
      MINL = 0.
      INDL = 0.
      DO I=1,10
        MINL(I,MYRANK+1)=VALIN(I)
        K = (I-1)*3
        INDL(1+K,MYRANK+1)=IND(1,I)
        INDL(2+K,MYRANK+1)=IND(2,I)
        INDL(3+K,MYRANK+1)=IND(3,I)+MYRANK*(NZ-2)
      ENDDO
      IF(ISGS/=0) THEN
        DO I=11,12
          MINL(I,MYRANK+1)=VALIN(I)
          K = (I-1)*3
          INDL(1+K,MYRANK+1)=IND(1,I)
          INDL(2+K,MYRANK+1)=IND(2,I)
          INDL(3+K,MYRANK+1)=IND(3,I)+MYRANK*(NZ-2)
        ENDDO
      ENDIF

      CALL MPI_REDUCE(MINL,MING,MYSIZE*12,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(INDL,INDG,MYSIZE*3*12,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
c
c.....only root write to screen
c      
      IF(MYRANK==0) THEN

        MINRANK = MINLOC(MING,DIM=2)
        DO I=1,12
          J = MINRANK(I)
          MINV(I) = MING(I,J)
          K = (I-1)*3
          INDC(1,I) = INDG(1+K,J)
          INDC(2,I) = INDG(2+K,J)
          INDC(3,I) = INDG(3+K,J)
        ENDDO
c
c.....write screen infomation
c
        write(6,'(79("*"))')
        write(6,'(A,I10)'   ) '*..Iteration  = ',icycle
        write(6,'(A,G24.16)') '*..Time Step  = ',dtm
        write(6,'(A,G24.16)') '*..Total Time = ',time           

        WRITE(6,100) '*..u   =', 
     &        MINV(1),' @(',INDC(1,1),INDC(2,1),INDC(3,1),');',
     &       -MINV(2),' @(',INDC(1,2),INDC(2,2),INDC(3,2),')'
        WRITE(6,100) '*..v   =',
     &        MINV(3),' @(',INDC(1,3),INDC(2,3),INDC(3,3),');',
     &       -MINV(4),' @(',INDC(1,4),INDC(2,4),INDC(3,4),')'
        WRITE(6,100) '*..w   =',
     &        MINV(5),' @(',INDC(1,5),INDC(2,5),INDC(3,5),');',
     &       -MINV(6),' @(',INDC(1,6),INDC(2,6),INDC(3,6),')'
        WRITE(6,100) '*..p   =',
     &        MINV(7),' @(',INDC(1,7),INDC(2,7),INDC(3,7),');',
     &       -MINV(8),' @(',INDC(1,8),INDC(2,8),INDC(3,8),')'
        WRITE(6,100) '*..divv=',
     &        MINV( 9),' @(',INDC(1, 9),INDC(2, 9),INDC(3, 9),');',
     &       -MINV(10),' @(',INDC(1,10),INDC(2,10),INDC(3,10),')'
C
        IF(ISGS/=0) THEN
          WRITE(6,100) '*..tv  =',
     &        MINV(11),' @(',INDC(1,11),INDC(2,11),INDC(3,11),');',
     &       -MINV(12),' @(',INDC(1,12),INDC(2,12),INDC(3,12),')'
        ENDIF

      ENDIF
      DEALLOCATE(MINL,MING,INDL,INDG)
C
 100  FORMAT(A,G19.12,A,3I5,A,G19.12,A,3I5,A)
      RETURN 
      END 


*---  real function tclock ----------------------------------------
*
*     PURPOSE: Returns the cpu time 
*
*------------------------------------------------------------------
      REAL FUNCTION TCLOCK()

      include 'mpif.h'

      tclock = mpi_wtime()
c      INTEGER COUNTS,RATE

c      CALL SYSTEM_CLOCK(COUNTS,RATE)
c      TCLOCK = real(counts)/real(rate)

      RETURN

      END
*------------------------------------------------------------------
