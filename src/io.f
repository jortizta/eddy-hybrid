c
c     ***************************************************************
c
C     subroutine for I/O - read or write a vector quantity
c
c     ***************************************************************
c
      SUBROUTINE IOVECTOR(NAME,UO,VO,WO,UT,VT,WT,NX,NY,NZ,DIR,TIME)
      
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'

      CHARACTER NAME*(*)
      INTEGER DIR,NX,NY,NZ
      REAL    UO(NX,NY,NZ),VO(NX,NY,NZ),WO(NX,NY,NZ)
      REAL    UT(NX,NY,NZ),VT(NX,NY,NZ),WT(NX,NY,NZ)
      REAL    TIME
C
      INTEGER I,J,K,JP,STATUS(MPI_STATUS_SIZE)
      
c     if dir = -1 read else write

      IF(DIR==-1) THEN
        OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
c
        READ(19) I,J,K,JP,TIME
c
        DO K=1,MYRANK*(NZ-2)
          READ(19)
        ENDDO
c
        DO K=1,NZ
          READ(19) ((UO(I,J,K),I=1,NX),J=1,NY)
     &         ,   ((VO(I,J,K),I=1,NX),J=1,NY)
     &         ,   ((WO(I,J,K),I=1,NX),J=1,NY)
        ENDDO
c
        DO K=1,(MYSIZE-(MYRANK+1))*(NZ-2)
          READ(19)
        ENDDO

        CLOSE(19)

      ELSEIF(DIR== 1) THEN
c     note that all files read the input, but only root writes to output
        IF(MYRANK==0) THEN
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
c
          WRITE(19) NX,NY,MYSIZE*(NZ-2)+2,3
          IF(MYSIZE==1) THEN
            DO K=1,NZ
              WRITE(19) ((UO(I,J,K),I=1,NX),J=1,NY)
     &             ,    ((VO(I,J,K),I=1,NX),J=1,NY)
     &             ,    ((WO(I,J,K),I=1,NX),J=1,NY)
            ENDDO
          ELSE
c the processor 0 writes the local values on a file
            DO K=1,KZ2
              WRITE(19) ((UO(I,J,K),I=1,NX),J=1,NY)
     &             ,    ((VO(I,J,K),I=1,NX),J=1,NY)
     &             ,    ((WO(I,J,K),I=1,NX),J=1,NY)
            ENDDO
c
c the processor 0 receives the values from the other processors
c and writes them
c                               
            DO JP=1,MYSIZE-1
              CALL MPI_RECV(UT(1,1,1),NX*NY*NZ,MTYPE,JP,0,
     &             MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(VT(1,1,1),NX*NY*NZ,MTYPE,JP,0,
     &             MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(WT(1,1,1),NX*NY*NZ,MTYPE,JP,0,
     &             MPI_COMM_EDDY,STATUS,IERR)
              DO K=KZ1,KZ2+(JP+1)/MYSIZE
                WRITE(19) ((UT(I,J,K),I=1,NX),J=1,NY)
     &               ,    ((VT(I,J,K),I=1,NX),J=1,NY)
     &               ,    ((WT(I,J,K),I=1,NX),J=1,NY)
              ENDDO
            ENDDO
c
          ENDIF
c
          CLOSE(19)
        ELSE
c the local values are sent to the processor 0
          CALL MPI_SEND(UO(1,1,1),NX*NY*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
          CALL MPI_SEND(VO(1,1,1),NX*NY*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
          CALL MPI_SEND(WO(1,1,1),NX*NY*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        ENDIF
      ELSEIF(DIR==2) THEN
c
c     note that all files write to output
c
        OPEN(19,FILE=NAME//PROC(MYRANK),
     $       STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(19) NX,NY,NZ,3,TIME
        DO K=1,NZ
          WRITE(19) ((UO(I,J,K),I=1,NX),J=1,NY)
     &         ,    ((VO(I,J,K),I=1,NX),J=1,NY)
     &         ,    ((WO(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        CLOSE(19)
      ENDIF
c
      return
      end
c
c     ***************************************************************
c
C     subroutine for I/O - read or write a scalar quantity
c
c     ***************************************************************
c
      SUBROUTINE IOSCALAR(NAME,P,DP,NX,NY,NZ,DIR,TIME,DTM1,nstep)
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
      CHARACTER NAME*(*)
      INTEGER DIR,NX,NY,NZ,nstep
      REAL P(NX,NY,NZ),DP(NX,NY,NZ)
c      REAL*4 P1(NX,NY,NZ)
      REAL TIME,DTM1
C
      INTEGER I,J,K,JP,STATUS(MPI_STATUS_SIZE)

C     IF DIR = -1 READ ELSE WRITE

      IF(DIR==-1) THEN
        OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
        READ(19) I,J,K,JP 
        IF(MYSIZE/=1) THEN
          DO K=1,MYRANK*(NZ-2)
            READ(19)
          ENDDO
        ENDIF
        DO K=1,NZ
          READ(19) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        DO K=1,(MYSIZE-(MYRANK+1))*(NZ-2)
          READ(19)
        ENDDO
        READ(19) nstep
        READ(19) TIME
        if(myrank.eq.0) write(6,*) 'time=',time
        READ(19) DTM1,grav
        CLOSE(19)

      ELSEIF(DIR==1) THEN
c     note that all files read the input, but only root writes to output
        IF(MYRANK==0) THEN
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
          WRITE(19) NX,NY,MYSIZE*(NZ-2)+2,1
          IF(MYSIZE==1) THEN
c case of a serial run
            DO K=1,NZ
              WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
            ENDDO
          ELSE            
C
c the processor 0 writes the local values on a file
            DO K=1,KZ2
              WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
            ENDDO
c
c the processor 0 receives the values from the other processors
c and writes them
c
            DO JP=1,MYSIZE-1
              CALL MPI_RECV(DP(1,1,1),NX*NY*NZ,MTYPE,JP,0,
     &             MPI_COMM_EDDY,STATUS,IERR)
              DO K=KZ1,KZ2+(JP+1)/MYSIZE
                WRITE(19) ((DP(I,J,K),I=1,NX),J=1,NY)
              ENDDO
            ENDDO
C            
          ENDIF
          WRITE(19) nstep
          WRITE(19) TIME
          WRITE(19) DTM1,grav
          CLOSE(19)

        ELSE
c the local values are sent to the processor 0
          CALL MPI_SEND(P(1,1,1),NX*NY*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        ENDIF

      ELSEIF(DIR==2) THEN
c     note that all files write to output
        OPEN(19,FILE=NAME//PROC(MYRANK),
     $       STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(19) NX,NY,NZ,1,TIME,DTM1,grav
        DO K=1,NZ
          WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        CLOSE(19)
      ENDIF
c
      return
      end

!THIS SUBROUTINE WRITES DATA OF 3 SET OF PLANES FOR POSTPROCESSING!
      SUBROUTINE IOSCALAR_POST(NAME,P,DP,NX,NY,NZ,DIR,TIME,DTM1,nstep)
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
      CHARACTER NAME*(*)
      INTEGER DIR,NX,NY,NZ,nstep
      REAL P(NX,NY,NZ),DP(NX,NY,NZ)
c      REAL*4 P1(NX,NY,NZ)
      REAL TIME,DTM1
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: TMPF,TEMP
C
      INTEGER I,J,K,KS,JP,STATUS(MPI_STATUS_SIZE)
      INTEGER KX,KXX0,KXX1,KXX2,KXX3,MK,NKR,COUNTER,KMIN,KMAX,NK,NKRR

C     IF DIR = -1 READ ELSE WRITE

      IF(DIR==-1) THEN
        OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
        READ(19) I,J,K,JP 
        IF(MYSIZE/=1) THEN
          DO K=1,MYRANK*(NZ-2)
            READ(19)
          ENDDO
        ENDIF
        DO K=1,NZ
          READ(19) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        DO K=1,(MYSIZE-(MYRANK+1))*(NZ-2)
          READ(19)
        ENDDO
        READ(19) nstep
        READ(19) TIME
        write(6,*) 'time=',time
        READ(19) DTM1,grav
        CLOSE(19)

      ELSEIF(DIR==1) THEN
c     note that all files read the input, but only root writes to output
        IF(MYRANK==0) THEN
        KMAX=1114;KMIN=151;
        NKRR = KMAX-KMIN+1+4*(FLOOR((6144-KMAX)/20.0))
        ALLOCATE (TMPF(NX,NY,6144),TEMP(NX-10,NY,NKRR))
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
          WRITE(19) NX,NY,MYSIZE*(NZ-2)+2,1
          IF(MYSIZE==1) THEN
c case of a serial run
            DO K=1,NZ
              WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
            ENDDO
          ELSE            
C
c the processor 0 writes the local values on a file
            DO K=1,KZ2
             DO I=1,NX
              DO J=1,NY
              TMPF(I,J,K) = P(I,J,K)
              ENDDO
             ENDDO
            ENDDO
             !WRITE(*,'(2(4x,e15.8))')TMPF(2,33,6),P(2,33,6)

c
c the processor 0 receives the values from the other processors
c and writes them
c
            DO JP=1,MYSIZE-1
              CALL MPI_RECV(DP(1,1,1),NX*NY*NZ,MTYPE,JP,0,
     &             MPI_COMM_EDDY,STATUS,IERR)
              DO K=KZ1,KZ2+(JP+1)/MYSIZE
              KS = (KZ2-1)*JP + K
!              write(*,*)"ISSUE HERE1"
               DO I=1,NX
                DO J=1,NY
                TMPF(I,J,KS) = DP(I,J,K)
                ENDDO
               ENDDO
              ENDDO
!               IF(JP==1) WRITE(*,'(i6,2(4x,e15.8))')JP,TMPF(2,33,15),DP(2,33,6)
!              IF(JP==256) WRITE(*,'(i6,2(4x,e15.8))')JP,TMPF(2,33,2311),DP(2,33,7)

            ENDDO
C            
          ENDIF
            
            DO K=KMIN,KMAX
            KX = K-(KMIN-1)
             DO I=1,NX-10
              DO J=1,NY
!               if(K.gt.KMIN.or.K.lt.KMAX) then
              TEMP(I,J,KX) = TMPF(I,J,K)
!               endif
              ENDDO
             ENDDO
            ENDDO

!           write(*,*) 'ALREADY HERE1'
           MK = KX;COUNTER = 1;NK = 1;
           DO NKR=1,FLOOR((6144-KMAX)/20.0)
            MK = MK+1
            KXX0 = KMAX+NK*20-2
            KXX1 = KMAX+NK*20-1
            KXX2 = KMAX+NK*20
            KXX3 = KMAX+NK*20+1
!            write(*,*) 'MK',MK,'KXX1',KXX1,'KXX2',KXX2,'KXX3',KXX3
             DO I=1,NX-10
              DO J=1,NY
              TEMP(I,J,MK  ) = TMPF(I,J,KXX0)
              TEMP(I,J,MK+1) = TMPF(I,J,KXX1)
              TEMP(I,J,MK+2) = TMPF(I,J,KXX2)
              TEMP(I,J,MK+3) = TMPF(I,J,KXX3)
              ENDDO
             ENDDO
            MK = MK+3
            NK = NK+1
            ENDDO

!          WRITE(*,'(6(4x,e15.8))') TMPF(2,33,992),TMPF(2,33,993),TMPF(2,33,994),TMPF(2,33,1002),TMPF(2,33,1003),TMPF(2,33,1004)  
c          WRITE(*,'(6(4x,e15.8))') TEMP(2,33,767),TEMP(2,33,768),TEMP(2,33,769),TEMP(2,33,770),TEMP(2,33,771),TEMP(2,33,772)
!          WRITE(*,'(6(4x,e15.8))') TEMP(2,33,10),TEMP(2,33,1000),TEMP(2,33,1200),TEMP(2,33,1300),TEMP(2,33,1400),TEMP(2,33,1500)
          DO K=1,NKRR
          WRITE(19) ((TEMP(I,J,K),I=1,NX-10),J=1,NY)
          ENDDO
          WRITE(19) nstep
          WRITE(19) TIME
          WRITE(19) DTM1,grav
          CLOSE(19)

        ELSE
c the local values are sent to the processor 0
          CALL MPI_SEND(P(1,1,1),NX*NY*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        ENDIF

      ELSEIF(DIR==2) THEN
c     note that all files write to output
        OPEN(19,FILE=NAME//PROC(MYRANK),
     $       STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(19) NX,NY,NZ,1,TIME,DTM1,grav
        DO K=1,NZ
          WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        CLOSE(19)
      ENDIF
c
      return
      end


C   THIS  IS FOR GRID  
C$$$        KMAX=983;KMIN=218;
C$$$        NKRR = KMAX-KMIN+1+4*(FLOOR((4608-KMAX)/20.0))

C$$$!THIS SUBROUTINE WRITES DATA OF 3 SET OF PLANES FOR POSTPROCESSING!
C$$$      SUBROUTINE IOSCALAR_POST(NAME,P,DP,NX,NY,NZ,DIR,TIME,DTM1,nstep)
C$$$c
C$$$      INCLUDE 'common.h'
C$$$      INCLUDE 'mpif.h'
C$$$c
C$$$      CHARACTER NAME*(*)
C$$$      INTEGER DIR,NX,NY,NZ,nstep
C$$$      REAL P(NX,NY,NZ),DP(NX,NY,NZ)
C$$$c      REAL*4 P1(NX,NY,NZ)
C$$$      REAL TIME,DTM1
C$$$      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: TMPF,TEMP
C$$$C
C$$$      INTEGER I,J,K,KS,JP,STATUS(MPI_STATUS_SIZE)
C$$$      INTEGER KX,KXX0,KXX1,KXX2,KXX3,MK,NKR,COUNTER,KMIN,KMAX,NK,NKRR

C$$$C     IF DIR = -1 READ ELSE WRITE

C$$$      IF(DIR==-1) THEN
C$$$        OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
C$$$        READ(19) I,J,K,JP 
C$$$        IF(MYSIZE/=1) THEN
C$$$          DO K=1,MYRANK*(NZ-2)
C$$$            READ(19)
C$$$          ENDDO
C$$$        ENDIF
C$$$        DO K=1,NZ
C$$$          READ(19) ((P(I,J,K),I=1,NX),J=1,NY)
C$$$        ENDDO
C$$$        DO K=1,(MYSIZE-(MYRANK+1))*(NZ-2)
C$$$          READ(19)
C$$$        ENDDO
C$$$        READ(19) nstep
C$$$        READ(19) TIME
C$$$        write(6,*) 'time=',time
C$$$        READ(19) DTM1,grav
C$$$        CLOSE(19)

C$$$      ELSEIF(DIR==1) THEN
C$$$c     note that all files read the input, but only root writes to output
C$$$        IF(MYRANK==0) THEN
C$$$        KMAX=983;KMIN=218;
C$$$        NKRR = KMAX-KMIN+1+4*(FLOOR((4608-KMAX)/20.0))
C$$$        ALLOCATE (TMPF(NX,NY,4608),TEMP(NX-10,NY,NKRR))
C$$$          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
C$$$          WRITE(19) NX,NY,MYSIZE*(NZ-2)+2,1
C$$$          IF(MYSIZE==1) THEN
C$$$c case of a serial run
C$$$            DO K=1,NZ
C$$$              WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
C$$$            ENDDO
C$$$          ELSE            
C$$$C
C$$$c the processor 0 writes the local values on a file
C$$$            DO K=1,KZ2
C$$$             DO I=1,NX
C$$$              DO J=1,NY
C$$$              TMPF(I,J,K) = P(I,J,K)
C$$$              ENDDO
C$$$             ENDDO
C$$$            ENDDO
C$$$             !WRITE(*,'(2(4x,e15.8))')TMPF(2,33,6),P(2,33,6)

C$$$c
C$$$c the processor 0 receives the values from the other processors
C$$$c and writes them
C$$$c
C$$$            DO JP=1,MYSIZE-1
C$$$              CALL MPI_RECV(DP(1,1,1),NX*NY*NZ,MTYPE,JP,0,
C$$$     &             MPI_COMM_EDDY,STATUS,IERR)
C$$$              DO K=KZ1,KZ2+(JP+1)/MYSIZE
C$$$              KS = (KZ2-1)*JP + K
C$$$!              write(*,*)"ISSUE HERE1"
C$$$               DO I=1,NX
C$$$                DO J=1,NY
C$$$                TMPF(I,J,KS) = DP(I,J,K)
C$$$                ENDDO
C$$$               ENDDO
C$$$              ENDDO
C$$$!               IF(JP==1) WRITE(*,'(i6,2(4x,e15.8))')JP,TMPF(2,33,15),DP(2,33,6)
C$$$!              IF(JP==256) WRITE(*,'(i6,2(4x,e15.8))')JP,TMPF(2,33,2311),DP(2,33,7)

C$$$            ENDDO
C$$$C            
C$$$          ENDIF
            
C$$$            DO K=KMIN,KMAX
C$$$            KX = K-(KMIN-1)
C$$$             DO I=1,NX-10
C$$$              DO J=1,NY
C$$$!               if(K.gt.KMIN.or.K.lt.KMAX) then
C$$$              TEMP(I,J,KX) = TMPF(I,J,K)
C$$$!               endif
C$$$              ENDDO
C$$$             ENDDO
C$$$            ENDDO

C$$$!           write(*,*) 'ALREADY HERE1'
C$$$           MK = KX;COUNTER = 1;NK = 1;
C$$$           DO NKR=1,FLOOR((4608-KMAX)/20.0)
C$$$            MK = MK+1
C$$$            KXX0 = KMAX+NK*20-2
C$$$            KXX1 = KMAX+NK*20-1
C$$$            KXX2 = KMAX+NK*20
C$$$            KXX3 = KMAX+NK*20+1
C$$$!            write(*,*) 'MK',MK,'KXX1',KXX1,'KXX2',KXX2,'KXX3',KXX3
C$$$             DO I=1,NX-10
C$$$              DO J=1,NY
C$$$              TEMP(I,J,MK  ) = TMPF(I,J,KXX0)
C$$$              TEMP(I,J,MK+1) = TMPF(I,J,KXX1)
C$$$              TEMP(I,J,MK+2) = TMPF(I,J,KXX2)
C$$$              TEMP(I,J,MK+3) = TMPF(I,J,KXX3)
C$$$              ENDDO
C$$$             ENDDO
C$$$            MK = MK+3
C$$$            NK = NK+1
C$$$            ENDDO

C$$$!          WRITE(*,'(6(4x,e15.8))') TMPF(2,33,992),TMPF(2,33,993),TMPF(2,33,994),TMPF(2,33,1002),TMPF(2,33,1003),TMPF(2,33,1004)  
C$$$c          WRITE(*,'(6(4x,e15.8))') TEMP(2,33,767),TEMP(2,33,768),TEMP(2,33,769),TEMP(2,33,770),TEMP(2,33,771),TEMP(2,33,772)
C$$$!          WRITE(*,'(6(4x,e15.8))') TEMP(2,33,10),TEMP(2,33,1000),TEMP(2,33,1200),TEMP(2,33,1300),TEMP(2,33,1400),TEMP(2,33,1500)
C$$$          DO K=1,NKRR
C$$$          WRITE(19) ((TEMP(I,J,K),I=1,NX-10),J=1,NY)
C$$$          ENDDO
C$$$          WRITE(19) nstep
C$$$          WRITE(19) TIME
C$$$          WRITE(19) DTM1,grav
C$$$          CLOSE(19)

C$$$        ELSE
C$$$c the local values are sent to the processor 0
C$$$          CALL MPI_SEND(P(1,1,1),NX*NY*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
C$$$        ENDIF

C$$$      ELSEIF(DIR==2) THEN
C$$$c     note that all files write to output
C$$$        OPEN(19,FILE=NAME//PROC(MYRANK),
C$$$     $       STATUS='UNKNOWN',FORM='UNFORMATTED')
C$$$        WRITE(19) NX,NY,NZ,1,TIME,DTM1,grav
C$$$        DO K=1,NZ
C$$$          WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
C$$$        ENDDO
C$$$        CLOSE(19)
C$$$      ENDIF
C$$$c
C$$$      return
C$$$      end

      SUBROUTINE IOSCALAR_POST_5P(NAME,P,DP,NX,NY,NZ,DIR,TIME,DTM1,nstep)
!    THIS SUBROUTINE WRITES DATA OF 5 SET OF PLANES FOR POSTPROCESSING

      INCLUDE 'common.h'
      INCLUDE 'mpif.h'

      CHARACTER NAME*(*)
      INTEGER DIR,NX,NY,NZ,NZG,nstep
      REAL P(NX,NY,NZ),DP(NX,NY,NZ)
      REAL TIME,DTM1
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: TMPF,TEMP

      INTEGER I,J,K,KS,JP,STATUS(MPI_STATUS_SIZE)
      INTEGER KX,KXX0,KXX1,KXX2,KXX3,KXX4,MK,NKR,COUNTER,NK,NKRR
!     IF DIR = -1 READ ELSE WRITE

      NZG = MYSIZE*(NZ-2)+2
      NKRR = KMAX5P-KMIN5P+1+5*(FLOOR((NZG-2-KMAX5P)/20.0))

      IF(DIR==-1) THEN
        OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
        READ(19) I,J,NKRR,JP 
        IF(MYSIZE/=1) THEN
          DO K=1,MYRANK*(NZ-2)
            READ(19)
          ENDDO
        ENDIF
        DO K=1,NZ
          READ(19) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        DO K=1,(MYSIZE-(MYRANK+1))*(NZ-2)
          READ(19)
        ENDDO
        READ(19) nstep
        READ(19) TIME
        write(6,*) 'time=',time
        READ(19) DTM1,grav
        READ(19) KMIN5P,KMAX5P,NZG
        CLOSE(19)

      ELSEIF(DIR==1) THEN
!     note that all files read the input, but only root writes to output

        IF(MYRANK==0) THEN

        ALLOCATE (TMPF(NX,NY,NZG),TEMP(NX,NY,NKRR))
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
          WRITE(19) NX,NY,NKRR,1
          IF(MYSIZE==1) THEN
! case of a serial run
            DO K=1,NZ
              WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
            ENDDO
          ELSE            

! the processor 0 writes the local values on a file
            DO K=1,KZ2
             DO I=1,NX
              DO J=1,NY
              TMPF(I,J,K) = P(I,J,K)
              ENDDO
             ENDDO
            ENDDO

! the processor 0 receives the values from the other processors
! and writes them

            DO JP=1,MYSIZE-1
              CALL MPI_RECV(DP(1,1,1),NX*NY*NZ,MTYPE,JP,0,
     &             MPI_COMM_EDDY,STATUS,IERR)
              DO K=KZ1,KZ2+(JP+1)/MYSIZE
              KS = (KZ2-1)*JP + K
               DO I=1,NX
                DO J=1,NY
                TMPF(I,J,KS) = DP(I,J,K)
                ENDDO
               ENDDO
              ENDDO

            ENDDO
            
          ENDIF
            
            DO K=KMIN5P,KMAX5P

            KX = K-(KMIN5P-1)

             DO I=1,NX
              DO J=1,NY
              TEMP(I,J,KX) = TMPF(I,J,K)
              ENDDO
             ENDDO
            ENDDO

           MK = KX;COUNTER = 1;NK = 1;
           DO NKR=1,FLOOR((NZG-2-KMAX5P)/20.0)
            MK = MK+1
            KXX0 = KMAX5P+NK*20-2
            KXX1 = KMAX5P+NK*20-1
            KXX2 = KMAX5P+NK*20
            KXX3 = KMAX5P+NK*20+1
            KXX4 = KMAX5P+NK*20+2
             
             DO I=1,NX
              DO J=1,NY
              TEMP(I,J,MK  ) = TMPF(I,J,KXX0)
              TEMP(I,J,MK+1) = TMPF(I,J,KXX1)
              TEMP(I,J,MK+2) = TMPF(I,J,KXX2)
              TEMP(I,J,MK+3) = TMPF(I,J,KXX3)
              TEMP(I,J,MK+4) = TMPF(I,J,KXX4)
              ENDDO
             ENDDO
            MK = MK+4
            NK = NK+1
            ENDDO


          DO K=1,NKRR
          WRITE(19) ((TEMP(I,J,K),I=1,NX),J=1,NY)
          ENDDO
          WRITE(19) nstep
          WRITE(19) TIME
          WRITE(19) DTM1,grav
          WRITE(19) KMIN5P,KMAX5P,NZG
          CLOSE(19)

        ELSE
! the local values are sent to the processor 0
          CALL MPI_SEND(P(1,1,1),NX*NY*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        ENDIF

      ELSEIF(DIR==2) THEN
!     note that all files write to output
        OPEN(19,FILE=NAME//PROC(MYRANK),
     $       STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(19) NX,NY,NZ,1,TIME,DTM1,grav
        DO K=1,NZ
          WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        CLOSE(19)
      ENDIF

      return
      end


c     ***************************************************************
c
C     subroutine for I/O - read or write an immersed boundary arrays
c
c     ***************************************************************
c
      SUBROUTINE IOMFORC1(FNAME,IU,JU,KU,IMRKU,IUMTRX,JUMTRX,KUMTRX
     &      ,UINDX,UMTRX,UIM,XNU,YNU,ZNU,NXU,NYU,NZU,LIMU,MIMU,NZ,IDIR)

      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'mpif.h'


      CHARACTER FNAME*(*)
      INTEGER LIMU,MIMU,NZ,IDIR
      INTEGER IU(NFCMAX),JU(NFCMAX),KU(NFCMAX),IMRKU(NFCMAX)
      INTEGER IUMTRX(NSM,NFCMAX),JUMTRX(NSM,NFCMAX),KUMTRX(NSM,NFCMAX)
      INTEGER UINDX(NSP,NFCMAX)
      REAL    UMTRX(NSP,NSP,NFCMAX)
      REAL    UIM(NFCMAX)
      REAL    XNU(NFCMAX),YNU(NFCMAX),ZNU(NFCMAX)
      REAL    NXU(NFCMAX),NYU(NFCMAX),NZU(NFCMAX)

      REAL    UMTRX1(NSP,NSP,NFCMAX)
      INTEGER IU1(NFCMAX),JU1(NFCMAX),KU1(NFCMAX),IMRKU1(NFCMAX)
      INTEGER IUMTRX1(NSM,NFCMAX),JUMTRX1(NSM,NFCMAX),KUMTRX1(NSM,NFCMAX)
      INTEGER UINDX1(NSP,NFCMAX)
      REAL    UIM1(NFCMAX)
      REAL    XNU1(NFCMAX),YNU1(NFCMAX),ZNU1(NFCMAX)
      REAL    NXU1(NFCMAX),NYU1(NFCMAX),NZU1(NFCMAX)


      INTEGER MIMUG,ICOUNT,IPROC,MU,I,J,K,JR,II,STATUS(MPI_STATUS_SIZE)
      REAL    UMTRX2(NSP,NSP)
      INTEGER IU2,JU2,KU2,IMRKU2
      INTEGER IUMTRX2(NSM),JUMTRX2(NSM),KUMTRX2(NSM)
      INTEGER UINDX2(NSP)
      REAL    UIM2
      REAL    XNU2,YNU2,ZNU2
      REAL    NXU2,NYU2,NZU2


c      WRITE(6,'(4(1X,A,I9))') 'MYRANK=',MYRANK,'NFCMAX=',NFCMAX,'NSM=',NSM,'NSP=',NSP
      IF(IDIR==1) THEN

c        WRITE(6,'(2(1X,A,I6))')
c     &              'MYRANK=',MYRANK,'MIMU=',MIMU
        ICOUNT = MIMU
        CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
c        WRITE(6,'(3(1X,A,I6))')
c     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG

c        DO I=LIMU+1,LIMU+MIMU
c          IF(IU(i)==8 .AND. JU(I)==8 .AND. KU(I)+MYRANK*(NZ-2)==137) THEN
c            write(6,'(A,8(1x,I4),3(1X,F14.8))') 'TEST MFORC:',MYRANK,I
c     &               ,IU(I),JU(I),KU(I)+MYRANK*(NZ-2)
c     &               ,IUMTRX1(1,I),JUMTRX1(1,I),KUMTRX1(1,I)+MYRANK*(NZ-2)
c     &               ,XNU(I),YNU(I),ZNU(I)
c          ENDIF
c        ENDDO

        IF(MYRANK.EQ.0) THEN
          OPEN(20,FILE=FNAME,FORM='UNFORMATTED',STATUS='UNKNOWN')
          WRITE(20) MIMUG !'(I6)'
          DO I=LIMU+1,LIMU+MIMU
            WRITE(20) IU(I),JU(I),KU(I)
     &              , IMRKU(I)
     &              , (IUMTRX(J,I),J=1,NSM)
     &              , (JUMTRX(J,I),J=1,NSM)
     &              , (KUMTRX(J,I),J=1,NSM)
     &              , (UINDX(J,I),J=1,NSP)
     &              , ((UMTRX(K,J,I),J=1,NSP),K=1,NSP)
     &              , UIM(I),XNU(I),YNU(I),ZNU(I),NXU(I),NYU(I),NZU(I)
          ENDDO
          DO IPROC=1,MYSIZE-1
            CALL MPI_RECV(MU,1,MPI_INTEGER,IPROC,1,MPI_COMM_EDDY,STATUS,IERR)
            IF(MU/=0) THEN
              CALL MPI_RECV(IU1(1),MU,MPI_INTEGER,IPROC,2,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(JU1(1),MU,MPI_INTEGER,IPROC,3,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(KU1(1),MU,MPI_INTEGER,IPROC,4,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(IUMTRX1(1,1),NSM*MU,MPI_INTEGER,IPROC,5,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(JUMTRX1(1,1),NSM*MU,MPI_INTEGER,IPROC,6,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(KUMTRX1(1,1),NSM*MU,MPI_INTEGER,IPROC,7,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(UINDX1(1,1),NSP*MU,MPI_INTEGER,IPROC,8,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(UMTRX1(1,1,1),NSP*NSP*MU,MTYPE,IPROC,9,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(IMRKU1(1),MU,MPI_INTEGER,IPROC,10,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(UIM1(1),MU,MTYPE,IPROC,11,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(XNU1(1),MU,MTYPE,IPROC,12,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(YNU1(1),MU,MTYPE,IPROC,13,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZNU1(1),MU,MTYPE,IPROC,14,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NXU1(1),MU,MTYPE,IPROC,15,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NYU1(1),MU,MTYPE,IPROC,16,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NZU1(1),MU,MTYPE,IPROC,17,MPI_COMM_EDDY,STATUS,IERR)

              DO I=1,MU
              WRITE(20) IU1(I),JU1(I),KU1(I)      
     &                  , IMRKU1(I)
     &                  , (IUMTRX1(J,I),J=1,NSM)
     &                  , (JUMTRX1(J,I),J=1,NSM)
     &                  , (KUMTRX1(J,I),J=1,NSM)
     &                  , (UINDX1(J,I),J=1,NSP)
     &                  , ((UMTRX1(K,J,I),J=1,NSP),K=1,NSP)
     &                  , UIM1(I),XNU1(I),YNU1(I),ZNU1(I),NXU1(I),NYU1(I),NZU1(I)
              ENDDO
            ENDIF
          ENDDO
          CLOSE(20)

          ELSE
            CALL MPI_SEND(MIMU,1,MPI_INTEGER,0,1,MPI_COMM_EDDY,IERR)
            IF(MIMU/=0) THEN
              IU1(1:MIMU) = IU(LIMU+1:LIMU+MIMU)
              JU1(1:MIMU) = JU(LIMU+1:LIMU+MIMU)
              KU1(1:MIMU) = KU(LIMU+1:LIMU+MIMU) + (NZ-2)*MYRANK
              IMRKU1(1:MIMU) = IMRKU(LIMU+1:LIMU+MIMU)
              IUMTRX1(1:NSM,1:MIMU) = IUMTRX(1:NSM,LIMU+1:LIMU+MIMU)
              JUMTRX1(1:NSM,1:MIMU) = JUMTRX(1:NSM,LIMU+1:LIMU+MIMU)
              KUMTRX1(1:NSM,1:MIMU) = KUMTRX(1:NSM,LIMU+1:LIMU+MIMU) + (NZ-2)*MYRANK
              UINDX1(1:NSP,1:MIMU) = UINDX(1:NSP,LIMU+1:LIMU+MIMU)
              UMTRX1(1:NSP,1:NSP,1:MIMU) = UMTRX(1:NSP,1:NSP,LIMU+1:LIMU+MIMU)
              UIM1(1:MIMU) = UIM(LIMU+1:LIMU+MIMU)
              XNU1(1:MIMU) = XNU(LIMU+1:LIMU+MIMU)
              YNU1(1:MIMU) = YNU(LIMU+1:LIMU+MIMU)
              ZNU1(1:MIMU) = ZNU(LIMU+1:LIMU+MIMU)
              NXU1(1:MIMU) = NXU(LIMU+1:LIMU+MIMU)
              NYU1(1:MIMU) = NYU(LIMU+1:LIMU+MIMU)
              NZU1(1:MIMU) = NZU(LIMU+1:LIMU+MIMU)

c           IF(IU1(I)==3.AND.JU1(I)==6.AND.KU1(I)==12) THEN
c             WRITE(6,'(7(1x,I4))') IPROC,IU1(I),JU1(I),KU1(I),IMRKU1(I)
c           ENDIF

              CALL MPI_SEND(IU1(1),MIMU,MPI_INTEGER,0,2,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(JU1(1),MIMU,MPI_INTEGER,0,3,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(KU1(1),MIMU,MPI_INTEGER,0,4,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(IUMTRX1(1,1),NSM*MIMU,MPI_INTEGER,0,5,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(JUMTRX1(1,1),NSM*MIMU,MPI_INTEGER,0,6,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(KUMTRX1(1,1),NSM*MIMU,MPI_INTEGER,0,7,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(UINDX1(1,1),NSP*MIMU,MPI_INTEGER,0,8,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(UMTRX1(1,1,1),NSP*NSP*MIMU,MTYPE,0,9,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(IMRKU1(1),MIMU,MPI_INTEGER,0,10,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(UIM1(1),MIMU,MTYPE,0,11,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(XNU1(1),MIMU,MTYPE,0,12,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(YNU1(1),MIMU,MTYPE,0,13,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZNU1(1),MIMU,MTYPE,0,14,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NXU1(1),MIMU,MTYPE,0,15,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NYU1(1),MIMU,MTYPE,0,16,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NZU1(1),MIMU,MTYPE,0,17,MPI_COMM_EDDY,IERR)
            ENDIF

          ENDIF

        ELSEIF(IDIR==-1) THEN

          OPEN(20,FILE=FNAME,FORM='UNFORMATTED',STATUS='OLD')
          READ(20) MIMUG !'(I6)'
          MIMU = 0
          II = 0
          DO I=1,MIMUG
            READ(20) IU2,JU2,KU2
     &             , IMRKU2
     &             , (IUMTRX2(J),J=1,NSM)
     &             , (JUMTRX2(J),J=1,NSM)
     &             , (KUMTRX2(J),J=1,NSM)
     &             , (UINDX2(J),J=1,NSP)
     &             , ((UMTRX2(K,J),J=1,NSP),K=1,NSP)
     &             , UIM2,XNU2,YNU2,ZNU2,NXU2,NYU2,NZU2

            IF( KU2.GE.KZ1+(NZ-2)*MYRANK .AND. KU2.LE.KZ2+(NZ-2)*MYRANK ) THEN
              MIMU=MIMU+1
              II = LIMU+MIMU

              IU(II) = IU2
              JU(II) = JU2
              KU(II) = KU2 - (NZ-2)*MYRANK
              IMRKU(II) = IMRKU2
              IUMTRX(1:NSM,II) = IUMTRX2(1:NSM)
              JUMTRX(1:NSM,II) = JUMTRX2(1:NSM)
              KUMTRX(1:NSM,II) = KUMTRX2(1:NSM) - (NZ-2)*MYRANK
              UINDX(1:NSP,II) = UINDX2(1:NSP)
              UMTRX(1:NSP,1:NSP,II) = UMTRX2(1:NSP,1:NSP)
              UIM(II) = UIM2
              XNU(II) = XNU2
              YNU(II) = YNU2
              ZNU(II) = ZNU2
              NXU(II) = NXU2
              NYU(II) = NYU2
              NZU(II) = NZU2
c              WRITE(6,*) MIMU,IU(MIMU),IU2

            ENDIF
          ENDDO
          CLOSE(20)
          ICOUNT = MIMU
          CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
c          WRITE(6,'(3(1X,A,I6))')
c     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG
        ENDIF

        RETURN
 120    FORMAT(17(1X,I4),23(1X,F14.8))
c 120  FORMAT(9(1X,I4),11(1X,F12.5))
 122    FORMAT(13(1X,I8),17(1X,F14.8))

        END


c     ***************************************************************
c
C     subroutine for I/O - read or write an integer 3D array
c
c     ***************************************************************
c
      SUBROUTINE IOSCALAR2(NAME,U1,NX,NY,NZ,DIR)
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
      CHARACTER NAME*(*)
      INTEGER DIR,NX,NY,NZ
      INTEGER U1(NX,NY,NZ)
c      REAL*4 P1(NX,NY,NZ)
      REAL TIME
C
      INTEGER U2(NX,NY,NZ)
      INTEGER I,J,K,JP,STATUS(MPI_STATUS_SIZE)

C     IF DIR = -1 READ ELSE WRITE
c      WRITE(6,*) 'INSIDE IOSCALAR2'
      IF(DIR==-1) THEN
        OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
        READ(19) I,J,K,JP
        IF(MYSIZE/=1) THEN
          DO K=1,MYRANK*(NZ-2)
          DO J=1,NY
            READ(19)
          ENDDO
          ENDDO
        ENDIF
        DO K=1,NZ
         DO J=1,NY
          READ(19) (U1(I,J,K),I=1,NX)
         ENDDO
        ENDDO

        CLOSE(19)

      ELSEIF(DIR==1) THEN
c     note that all files read the input, but only root writes to output
        U2 = U1
        IF(MYRANK==0) THEN
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
          WRITE(19) NX,NY,MYSIZE*(NZ-2)+2,1 !'(4(1X,I4))'
          IF(MYSIZE==1) THEN
            DO K=1,NZ
             WRITE(19) ((REAL(U1(I,J,K),4),I=1,NX),J=1,NY)
c             DO J=1,NY
c              WRITE(19) (U1(I,J,K),I=1,NX)
c             ENDDO
            ENDDO
          ELSE
C
            DO K=1,KZ2
             WRITE(19) ((REAL(U1(I,J,K),4),I=1,NX),J=1,NY)
c             DO J=1,NY
c              WRITE(19) (U1(I,J,K),I=1,NX)
c             ENDDO
c            DO I=1,NX
c              WRITE(19,104) I,J,K,U1(I,J,K)
c            ENDDO
            ENDDO

            DO JP=1,MYSIZE-1
              CALL MPI_RECV(U2(1,1,1),NX*NY*NZ,MPI_INTEGER,JP,0,
     &             MPI_COMM_EDDY,STATUS,IERR)
              DO K=KZ1,KZ2+(JP+1)/MYSIZE
                WRITE(19) ((REAL(U2(I,J,K),4),I=1,NX),J=1,NY)
c               DO J=1,NY
c                WRITE(19) (U2(I,J,K),I=1,NX)
c               ENDDO
c              DO I=1,NX
c                WRITE(19,104) I,J,K+JP*(NZ-2),U2(I,J,K)
c              ENDDO
              ENDDO
            ENDDO
C
          ENDIF
          CLOSE(19)
c          WRITE(6,*) 'Closing file'
        ELSE
          CALL MPI_SEND(U1(1,1,1),NX*NY*NZ,MPI_INTEGER,0,0,MPI_COMM_EDDY,IERR)
        ENDIF
      ENDIF
c
      return
104   FORMAT(4(1X,I5))
      end






c     ***************************************************************
c
C     subroutine for I/O - read or write immersed boundary arrays.
C     Format of file is read my TECPLOT
c
c     ***************************************************************
c
      SUBROUTINE IOMFORC2(FNAME,XU,YU,ZU,NX,NY,NZ,IU,JU,KU
     &      ,MRKU,XNU,YNU,ZNU,NXU,NYU,NZU,LIMU,MIMU,IDIR)

      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'mpif.h'


      CHARACTER FNAME*(*)
      INTEGER LIMU,MIMU,NX,NY,NZ,IDIR
      INTEGER IU(NFCMAX),JU(NFCMAX),KU(NFCMAX),IMRKU(NFCMAX)
c      INTEGER IUMTRX(NSM,NFCMAX),JUMTRX(NSM,NFCMAX),KUMTRX(NSM,NFCMAX)
      INTEGER MRKU(NFCMAX)
c      INTEGER UINDX(NSP,NFCMAX)
c      REAL    UMTRX(NSP,NSP,NFCMAX)
c      REAL    UIM(NFCMAX)
      REAL    XNU(NFCMAX),YNU(NFCMAX),ZNU(NFCMAX)
      REAL    NXU(NFCMAX),NYU(NFCMAX),NZU(NFCMAX)
      REAL    XU(NX,NY),YU(NX,NY),ZU(NZ)

c      REAL    UMTRX1(NSP,NSP,NFCMAX)
      INTEGER LIMU1,MIMU1
      INTEGER IU1(NFCMAX),JU1(NFCMAX),KU1(NFCMAX),IMRKU1(NFCMAX)
c      INTEGER IUMTRX1(NSM,NFCMAX),JUMTRX1(NSM,NFCMAX),KUMTRX1(NSM,NFCMAX)
      INTEGER MRKU1(NFCMAX)
c      INTEGER UINDX1(NSP,NFCMAX)
c      REAL    UIM1(NFCMAX)
      REAL    XNU1(NFCMAX),YNU1(NFCMAX),ZNU1(NFCMAX)
      REAL    NXU1(NFCMAX),NYU1(NFCMAX),NZU1(NFCMAX)
      REAL    ZU1(NZ)

      INTEGER MIMUG,ICOUNT,IPROC,MU,I,J,K,JR,II,IM,STATUS(MPI_STATUS_SIZE)
      real    anglerad
c      REAL    UMTRX2(NSP,NSP)
c      INTEGER IU2,JU2,KU2,IMRKU2
c      INTEGER IUMTRX2(NSM),JUMTRX2(NSM),KUMTRX2(NSM)
c      INTEGER UINDX2(NSP)
c      REAL    UIM2
c      REAL    XNU2,YNU2,ZNU2
c      REAL    NXU2,NYU2,NZU2


c      WRITE(6,'(4(1X,A,I9))') 'MYRANK=',MYRANK,'NFCMAX=',NFCMAX,'NSM=',NSM,'NSP=',NSP
      IF(IDIR==1) THEN

c        WRITE(6,'(2(1X,A,I6))')
c     &              'MYRANK=',MYRANK,'MIMU=',MIMU
        ICOUNT = MIMU
        CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
        WRITE(6,'(3(1X,A,I6))')
     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG


        IF(MYRANK.EQ.0) THEN
          OPEN(20,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
          WRITE(20, '(A)') 'VARIABLES = "X" "Y" "Z"'
          DO IM=LIMU+1,LIMU+MIMU
c            WRITE(20, '(A)') 'ZONE I=2 DATAPACKING=POINT'
            if(mrku(im)==1) then
              WRITE(20, '(A)') 'ZONE I=2 DATAPACKING=POINT'
              if(icyl==0) then
                write(20,'(3(1x,F10.6))') nxu(im),nyu(im),nzu(im)
              else
                write(20,'(3(1x,F10.6))') sqrt(nxu(im)**2.+nyu(im)**2.),anglerad(nxu(im),nyu(im)),nzu(im)
              endif

            i = iu(im)
            j = ju(im)
            k = ku(im)
c            write(20, '(3(1X,F10.6))') xu(i,j),yu(i,j),zu(k)
            if(icyl==0) then
              write(20, '(3(1X,F10.6))') xnu(im),ynu(im),znu(im)
            else
              write(20, '(3(1X,F10.6))') sqrt(xnu(im)**2.+ynu(im)**2.),anglerad(xnu(im),ynu(im)),znu(im)
            endif

            elseif(.false.) then
c              i = iumtrx(1,im)
c              j = jumtrx(1,im)
c              k = kumtrx(1,im)
              i = iu(im)+int(nxu(im))
              j = ju(im)+int(nyu(im))
              k = ku(im)+int(nzu(im))
              if(icyl==0) then
                write(20, '(3(1X,F10.6))') xu(i,j),yu(i,j),zu(k)
              else
                write(20, '(3(1X,F10.6))') sqrt(xu(i,j)**2.+yu(i,j)**2.),anglerad(xu(i,j),yu(i,j)),zu(k)
              endif
c              if(iu(im)==11) then
c                write(6,*) 'mforc2:',iu(im),ju(im),ku(im),i,j,k,sqrt(xnu(im)**2.+ynu(im)**2.),anglerad(xnu(im),ynu(im)),znu(im)
c              endif
c            endif
            i = iu(im)
            j = ju(im)
            k = ku(im)
c            write(20, '(3(1X,F10.6))') xu(i,j),yu(i,j),zu(k)
            if(icyl==0) then
              write(20, '(3(1X,F10.6))') xnu(im),ynu(im),znu(im)
            else
              write(20, '(3(1X,F10.6))') sqrt(xnu(im)**2.+ynu(im)**2.),anglerad(xnu(im),ynu(im)),znu(im)
            endif
          endif
c            write(20, *) i,j,k,mrku(im)
          ENDDO

          DO IPROC=1,MYSIZE-1
            CALL MPI_RECV(MU,1,MPI_INTEGER,IPROC,1,MPI_COMM_EDDY,STATUS,IERR)
            IF(MU/=0) THEN
              CALL MPI_RECV(IU1(1),MU,MPI_INTEGER,IPROC,2,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(JU1(1),MU,MPI_INTEGER,IPROC,3,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(KU1(1),MU,MPI_INTEGER,IPROC,4,MPI_COMM_EDDY,STATUS,IERR)
c              CALL MPI_RECV(IUMTRX1(1,1),NSM*MU,MPI_INTEGER,IPROC,5,MPI_COMM_EDDY,STATUS,IERR)
c              CALL MPI_RECV(JUMTRX1(1,1),NSM*MU,MPI_INTEGER,IPROC,6,MPI_COMM_EDDY,STATUS,IERR)
c              CALL MPI_RECV(KUMTRX1(1,1),NSM*MU,MPI_INTEGER,IPROC,7,MPI_COMM_EDDY,STATUS,IERR)
c              CALL MPI_RECV(UINDX1(1,1),NSP*MU,MPI_INTEGER,IPROC,8,MPI_COMM_EDDY,STATUS,IERR)
c              CALL MPI_RECV(UMTRX1(1,1,1),NSP*NSP*MU,MTYPE,IPROC,9,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(MRKU1(1),MU,MPI_INTEGER,IPROC,10,MPI_COMM_EDDY,STATUS,IERR)
c              CALL MPI_RECV(UIM1(1),MU,MTYPE,IPROC,11,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(XNU1(1),MU,MTYPE,IPROC,12,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(YNU1(1),MU,MTYPE,IPROC,13,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZNU1(1),MU,MTYPE,IPROC,14,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NXU1(1),MU,MTYPE,IPROC,15,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NYU1(1),MU,MTYPE,IPROC,16,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NZU1(1),MU,MTYPE,IPROC,17,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZU1(1),NZ,MTYPE,IPROC,18,MPI_COMM_EDDY,STATUS,IERR)

              DO IM=1,MU
                WRITE(20, '(A)') 'ZONE I=2 DATAPACKING=POINT'
                if(mrku1(im)==1) then
                  write(20,'(3(1x,F10.6))') nxu1(im),nyu1(im),nzu1(im)
                else
c                  I=IUMTRX1(1,IM)
c                  J=JUMTRX1(1,IM)
c                  K=KUMTRX1(1,IM)
                  i=iu1(im)+int(nxu1(im))
                  j=ju1(im)+int(nyu1(im))
                  k=ku1(im)+int(nzu1(im))
                  WRITE(20, '(3(1X,F10.6))') XU(I,J),YU(I,J),ZU1(K)
                endif
                i = iu1(im)
                j = ju1(im)
                k = ku1(im)
c                write(20, '(3(1X,F10.6))') xu(i,j),yu(i,j),zu1(k)
                WRITE(20, '(3(1X,F10.6))') XNU1(IM),YNU1(IM),ZNU1(IM)
                write(20, *) i,j,k+myrank*(nz-2)
              ENDDO

            ENDIF
          ENDDO
          CLOSE(20)

          ELSE
            CALL MPI_SEND(MIMU,1,MPI_INTEGER,0,1,MPI_COMM_EDDY,IERR)
            IF(MIMU/=0) THEN
              IU1(1:MIMU) = IU(LIMU+1:LIMU+MIMU)
              JU1(1:MIMU) = JU(LIMU+1:LIMU+MIMU)
              KU1(1:MIMU) = KU(LIMU+1:LIMU+MIMU)! + (NZ-2)*MYRANK
              MRKU1(1:MIMU) = MRKU(LIMU+1:LIMU+MIMU)
c              IUMTRX1(1:NSM,1:MIMU) = IUMTRX(1:NSM,LIMU+1:LIMU+MIMU)
c              JUMTRX1(1:NSM,1:MIMU) = JUMTRX(1:NSM,LIMU+1:LIMU+MIMU)
c              KUMTRX1(1:NSM,1:MIMU) = KUMTRX(1:NSM,LIMU+1:LIMU+MIMU)! + (NZ-2)*MYRANK
c              UINDX1(1:NSP,1:MIMU) = UINDX(1:NSP,LIMU+1:LIMU+MIMU)
c              UMTRX1(1:NSP,1:NSP,1:MIMU) = UMTRX(1:NSP,1:NSP,LIMU+1:LIMU+MIMU)
c              UIM1(1:MIMU) = UIM(LIMU+1:LIMU+MIMU)
              XNU1(1:MIMU) = XNU(LIMU+1:LIMU+MIMU)
              YNU1(1:MIMU) = YNU(LIMU+1:LIMU+MIMU)
              ZNU1(1:MIMU) = ZNU(LIMU+1:LIMU+MIMU)
              NXU1(1:MIMU) = NXU(LIMU+1:LIMU+MIMU)
              NYU1(1:MIMU) = NYU(LIMU+1:LIMU+MIMU)
              NZU1(1:MIMU) = NZU(LIMU+1:LIMU+MIMU)

              CALL MPI_SEND(IU1(1),MIMU,MPI_INTEGER,0,2,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(JU1(1),MIMU,MPI_INTEGER,0,3,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(KU1(1),MIMU,MPI_INTEGER,0,4,MPI_COMM_EDDY,IERR)
c              CALL MPI_SEND(IUMTRX1(1,1),NSM*MIMU,MPI_INTEGER,0,5,MPI_COMM_EDDY,IERR)
c              CALL MPI_SEND(JUMTRX1(1,1),NSM*MIMU,MPI_INTEGER,0,6,MPI_COMM_EDDY,IERR)
c              CALL MPI_SEND(KUMTRX1(1,1),NSM*MIMU,MPI_INTEGER,0,7,MPI_COMM_EDDY,IERR)
c              CALL MPI_SEND(UINDX1(1,1),NSP*MIMU,MPI_INTEGER,0,8,MPI_COMM_EDDY,IERR)
c              CALL MPI_SEND(UMTRX1(1,1,1),NSP*NSP*MIMU,MTYPE,0,9,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(MRKU1(1),MIMU,MPI_INTEGER,0,10,MPI_COMM_EDDY,IERR)
c              CALL MPI_SEND(UIM1(1),MIMU,MTYPE,0,11,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(XNU1(1),MIMU,MTYPE,0,12,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(YNU1(1),MIMU,MTYPE,0,13,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZNU1(1),MIMU,MTYPE,0,14,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NXU1(1),MIMU,MTYPE,0,15,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NYU1(1),MIMU,MTYPE,0,16,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NZU1(1),MIMU,MTYPE,0,17,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZU(1),NZ,MTYPE,0,18,MPI_COMM_EDDY,IERR)
            ENDIF

          ENDIF

          ICOUNT = MIMU
          CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
c          WRITE(6,'(3(1X,A,I6))')
c     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG
        ENDIF

        RETURN
 120    FORMAT(17(1X,I4),23(1X,F14.8))
c 120  FORMAT(9(1X,I4),11(1X,F12.5))
 122    FORMAT(13(1X,I8),17(1X,F14.8))

        END

c
c     ***************************************************************
c
C     subroutine for I/O - read or write a scalar quantity in single precision
c
c     ***************************************************************
c
      SUBROUTINE IOSCALAR3(NAME,P,NX,NY,NZ,DIR,TIME)
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
      CHARACTER NAME*(*)
      INTEGER DIR,NX,NY,NZ
      REAL P(NX,NY,NZ)
      REAL TIME
C
      INTEGER I,J,K,JP,STATUS(MPI_STATUS_SIZE)
      REAL    DP(NX,NY,NZ)
      REAL*4  P1(NX,NY,NZ)

C     IF DIR = -1 READ ELSE WRITE

      IF(DIR==-1) THEN
        OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
        READ(19) I,J,K,JP
        IF(MYSIZE/=1) THEN
          DO K=1,MYRANK*(NZ-2)
            READ(19)
          ENDDO
        ENDIF
        DO K=1,NZ
          READ(19) ((P1(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        DO K=1,(MYSIZE-(MYRANK+1))*(NZ-2)
          READ(19)
        ENDDO
        READ(19) TIME
        CLOSE(19)
        P = P1

      ELSEIF(DIR==1) THEN
c     note that all files read the input, but only root writes to output
        IF(MYRANK==0) THEN
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
          WRITE(19) NX,NY,MYSIZE*(NZ-2)+2,1
          IF(MYSIZE==1) THEN
            DO K=1,NZ
              WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
            ENDDO
          ELSE            
C
            DO K=1,KZ2
              WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
            ENDDO

            DO JP=1,MYSIZE-1
              CALL MPI_RECV(DP(1,1,1),NX*NY*NZ,MTYPE,JP,0,
     &             MPI_COMM_EDDY,STATUS,IERR)
              DO K=KZ1,KZ2+(JP+1)/MYSIZE
                WRITE(19) ((DP(I,J,K),I=1,NX),J=1,NY)
              ENDDO
            ENDDO
C            
          ENDIF
          WRITE(19) TIME
          CLOSE(19)

        ELSE
          CALL MPI_SEND(P(1,1,1),NX*NY*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        ENDIF
      ELSEIF(DIR==2) THEN
c     note that all files write to output
        OPEN(19,FILE=NAME//PROC(MYRANK),
     $       STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(19) NX,NY,NZ,1
        DO K=1,NZ
          WRITE(19) ((REAL(P(I,J,K),4),I=1,NX),J=1,NY)
        ENDDO
        CLOSE(19)
      ENDIF
c
      return
      end




c     ***************************************************************
c
C     subroutine for I/O - read or write immersed boundary arrays.
C     Format of file is read my TECPLOT
c
c     ***************************************************************
c
      SUBROUTINE IOMFORC3(FNAME,XU,YU,ZU,NX,NY,NZ,IU,JU,KU
     &      ,IUMTRX,JUMTRX,KUMTRX,MRKU
     &      ,XNU,YNU,ZNU,NXU,NYU,NZU,LIMU,MIMU,IDIR)

      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'mpif.h'


      CHARACTER FNAME*(*)
      INTEGER LIMU,MIMU,NX,NY,NZ,IDIR
      INTEGER IU(NFCMAX),JU(NFCMAX),KU(NFCMAX),IMRKU(NFCMAX)
      INTEGER IUMTRX(NSM,NFCMAX),JUMTRX(NSM,NFCMAX),KUMTRX(NSM,NFCMAX)
      INTEGER MRKU(NFCMAX)
      REAL    XNU(NFCMAX),YNU(NFCMAX),ZNU(NFCMAX)
      REAL    NXU(NFCMAX),NYU(NFCMAX),NZU(NFCMAX)
      REAL    XU(NX,NY),YU(NX,NY),ZU(NZ)

      INTEGER LIMU1,MIMU1
      INTEGER IU1(NFCMAX),JU1(NFCMAX),KU1(NFCMAX),IMRKU1(NFCMAX)
      INTEGER IUMTRX1(NSM,NFCMAX),JUMTRX1(NSM,NFCMAX),KUMTRX1(NSM,NFCMAX)
      INTEGER MRKU1(NFCMAX)
      REAL    XNU1(NFCMAX),YNU1(NFCMAX),ZNU1(NFCMAX)
      REAL    NXU1(NFCMAX),NYU1(NFCMAX),NZU1(NFCMAX)
      REAL    ZU1(NZ)

      INTEGER MIMUG,ICOUNT,IPROC,MU,I,J,K,JR,II,IM,STATUS(MPI_STATUS_SIZE)

      IF(IDIR==1) THEN

        ICOUNT = MIMU
        CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
        WRITE(6,'(3(1X,A,I6))')
     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG


        IF(MYRANK.EQ.0) THEN
          OPEN(20,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
          WRITE(20, '(A)') 'VARIABLES = "X" "Y" "Z"'
          DO IM=LIMU+1,LIMU+MIMU
            WRITE(20, '(A)') 'ZONE I=2 DATAPACKING=POINT'
            i = iumtrx(1,im)
            j = jumtrx(1,im)
            k = kumtrx(1,im)
            WRITE(20, '(3(1X,F10.6))') xu(i,j),yu(i,j),zu(k)
            write(20, '(3(1X,F10.6))') xnu(im),ynu(im),znu(im)
          ENDDO

          DO IPROC=1,MYSIZE-1
            CALL MPI_RECV(MU,1,MPI_INTEGER,IPROC,1,MPI_COMM_EDDY,STATUS,IERR)
            IF(MU/=0) THEN
              CALL MPI_RECV(IU1(1),MU,MPI_INTEGER,IPROC,2,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(JU1(1),MU,MPI_INTEGER,IPROC,3,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(KU1(1),MU,MPI_INTEGER,IPROC,4,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(IUMTRX1(1,1),NSM*MU,MPI_INTEGER,IPROC,5,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(JUMTRX1(1,1),NSM*MU,MPI_INTEGER,IPROC,6,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(KUMTRX1(1,1),NSM*MU,MPI_INTEGER,IPROC,7,MPI_COMM_EDDY,STATUS,IERR)
c              CALL MPI_RECV(UINDX1(1,1),NSP*MU,MPI_INTEGER,IPROC,8,MPI_COMM_EDDY,STATUS,IERR)
c              CALL MPI_RECV(UMTRX1(1,1,1),NSP*NSP*MU,MTYPE,IPROC,9,MPI_COMM_EDDY,STATUS,IERR)
c              CALL MPI_RECV(MRKU1(1),MU,MPI_INTEGER,IPROC,10,MPI_COMM_EDDY,STATUS,IERR)
c              CALL MPI_RECV(UIM1(1),MU,MTYPE,IPROC,11,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(XNU1(1),MU,MTYPE,IPROC,12,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(YNU1(1),MU,MTYPE,IPROC,13,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZNU1(1),MU,MTYPE,IPROC,14,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NXU1(1),MU,MTYPE,IPROC,15,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NYU1(1),MU,MTYPE,IPROC,16,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NZU1(1),MU,MTYPE,IPROC,17,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZU1(1),NZ,MTYPE,IPROC,18,MPI_COMM_EDDY,STATUS,IERR)

              DO IM=1,MU
                WRITE(20, '(A)') 'ZONE I=2 DATAPACKING=POINT'
                I=IUMTRX1(1,IM)
                J=JUMTRX1(1,IM)
                K=KUMTRX1(1,IM)
                WRITE(20, '(3(1X,F10.6))') XU(I,J),YU(I,J),ZU1(K)
                WRITE(20, '(3(1X,F10.6))') XNU1(IM),YNU1(IM),ZNU1(IM)
              ENDDO

            ENDIF
          ENDDO
          CLOSE(20)

          ELSE
            CALL MPI_SEND(MIMU,1,MPI_INTEGER,0,1,MPI_COMM_EDDY,IERR)
            IF(MIMU/=0) THEN
              IU1(1:MIMU) = IU(LIMU+1:LIMU+MIMU)
              JU1(1:MIMU) = JU(LIMU+1:LIMU+MIMU)
              KU1(1:MIMU) = KU(LIMU+1:LIMU+MIMU)! + (NZ-2)*MYRANK
c              MRKU1(1:MIMU) = MRKU(LIMU+1:LIMU+MIMU)
              IUMTRX1(1:NSM,1:MIMU) = IUMTRX(1:NSM,LIMU+1:LIMU+MIMU)
              JUMTRX1(1:NSM,1:MIMU) = JUMTRX(1:NSM,LIMU+1:LIMU+MIMU)
              KUMTRX1(1:NSM,1:MIMU) = KUMTRX(1:NSM,LIMU+1:LIMU+MIMU)! + (NZ-2)*MYRANK
c              UINDX1(1:NSP,1:MIMU) = UINDX(1:NSP,LIMU+1:LIMU+MIMU)
c              UMTRX1(1:NSP,1:NSP,1:MIMU) = UMTRX(1:NSP,1:NSP,LIMU+1:LIMU+MIMU)
c              UIM1(1:MIMU) = UIM(LIMU+1:LIMU+MIMU)
              XNU1(1:MIMU) = XNU(LIMU+1:LIMU+MIMU)
              YNU1(1:MIMU) = YNU(LIMU+1:LIMU+MIMU)
              ZNU1(1:MIMU) = ZNU(LIMU+1:LIMU+MIMU)
              NXU1(1:MIMU) = NXU(LIMU+1:LIMU+MIMU)
              NYU1(1:MIMU) = NYU(LIMU+1:LIMU+MIMU)
              NZU1(1:MIMU) = NZU(LIMU+1:LIMU+MIMU)

              CALL MPI_SEND(IU1(1),MIMU,MPI_INTEGER,0,2,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(JU1(1),MIMU,MPI_INTEGER,0,3,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(KU1(1),MIMU,MPI_INTEGER,0,4,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(IUMTRX1(1,1),NSM*MIMU,MPI_INTEGER,0,5,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(JUMTRX1(1,1),NSM*MIMU,MPI_INTEGER,0,6,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(KUMTRX1(1,1),NSM*MIMU,MPI_INTEGER,0,7,MPI_COMM_EDDY,IERR)
c              CALL MPI_SEND(UINDX1(1,1),NSP*MIMU,MPI_INTEGER,0,8,MPI_COMM_EDDY,IERR)
c              CALL MPI_SEND(UMTRX1(1,1,1),NSP*NSP*MIMU,MTYPE,0,9,MPI_COMM_EDDY,IERR)
c              CALL MPI_SEND(MRKU1(1),MIMU,MPI_INTEGER,0,10,MPI_COMM_EDDY,IERR)
c              CALL MPI_SEND(UIM1(1),MIMU,MTYPE,0,11,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(XNU1(1),MIMU,MTYPE,0,12,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(YNU1(1),MIMU,MTYPE,0,13,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZNU1(1),MIMU,MTYPE,0,14,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NXU1(1),MIMU,MTYPE,0,15,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NYU1(1),MIMU,MTYPE,0,16,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NZU1(1),MIMU,MTYPE,0,17,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZU(1),NZ,MTYPE,0,18,MPI_COMM_EDDY,IERR)
            ENDIF

          ENDIF

          ICOUNT = MIMU
          CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
c          WRITE(6,'(3(1X,A,I6))')
c     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG
        ENDIF

        RETURN
 120    FORMAT(17(1X,I4),23(1X,F14.8))
c 120  FORMAT(9(1X,I4),11(1X,F12.5))
 122    FORMAT(13(1X,I8),17(1X,F14.8))

        END



c
c     ***************************************************************
c
C     subroutine for I/O - read or write a scalar quantity
c
c     ***************************************************************
c
      SUBROUTINE IOSCALAR2D(NAME,P,NX,NZ,DIR,TIME)
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
      CHARACTER NAME*(*)
      INTEGER DIR,NX,NZ
      REAL P(NX,NZ)
      REAL TIME
C
      INTEGER I,J,K,JP,STATUS(MPI_STATUS_SIZE)
      REAL    DP(NX,NZ)

C     IF DIR = -1 READ ELSE WRITE

      IF(DIR==-1) THEN

      ELSEIF(DIR==1) THEN
c     note that all files read the input, but only root writes to output
        IF(MYRANK==0) THEN
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
          WRITE(19) NX,MYSIZE*(NZ-2)+2,1
          IF(MYSIZE==1) THEN
            DO K=1,NZ
              WRITE(19) (P(I,K),I=1,NX)
            ENDDO
c            WRITE(19) ((P(I,K),I=1,NX),K=1,NZ)
          ELSE            
C
            WRITE(19) ((P(I,K),I=1,NX),K=1,KZ2)

            DO JP=1,MYSIZE-1
              CALL MPI_RECV(DP(1,1),NX*NZ,MTYPE,JP,0,
     &             MPI_COMM_EDDY,STATUS,IERR)
             WRITE(19) ((P(I,K),I=1,NX),K=KZ1,KZ2+(JP+1)/MYSIZE)
            ENDDO
C            
          ENDIF
          WRITE(19) TIME
          CLOSE(19)

        ELSE
          CALL MPI_SEND(P(1,1),NX*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        ENDIF
      ELSEIF(DIR==2) THEN
c     note that all files write to output
        OPEN(19,FILE=NAME//PROC(MYRANK),
     $       STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(19) NX,NZ,1
        WRITE(19) ((P(I,K),I=1,NX),K=1,NZ)
        CLOSE(19)
      ENDIF
c
      return
      end



c
c     ***************************************************************
c
C     subroutine for I/O - read or write a scalar quantity
c
c     ***************************************************************
c
      SUBROUTINE IOSCALAR4(NAME,P,NX,NY,NZ,DIR,TIME)
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
      CHARACTER NAME*(*)
      INTEGER DIR,NX,NY,NZ
      REAL P(NX,NY,NZ),DP(NX,NY,NZ)
c      REAL*4 P1(NX,NY,NZ)
      REAL TIME
C
      INTEGER I,J,K,JP,STATUS(MPI_STATUS_SIZE)

C     IF DIR = -1 READ ELSE WRITE

      IF(DIR==-1) THEN
        OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
        READ(19) I,J,K,JP
        IF(MYSIZE/=1) THEN
          DO K=1,MYRANK*(NZ-2)
            READ(19)
          ENDDO
        ENDIF
        DO K=1,NZ
          READ(19) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        DO K=1,(MYSIZE-(MYRANK+1))*(NZ-2)
          READ(19)
        ENDDO
        READ(19) TIME
        CLOSE(19)

      ELSEIF(DIR==1) THEN
c     note that all files read the input, but only root writes to output
        IF(MYRANK==0) THEN
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
          WRITE(19) NX,NY,MYSIZE*(NZ-2)+2,1
          IF(MYSIZE==1) THEN
            DO K=1,NZ
              WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
            ENDDO
          ELSE            
C
            DO K=1,KZ2
              WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
            ENDDO

            DO JP=1,MYSIZE-1
              CALL MPI_RECV(DP(1,1,1),NX*NY*NZ,MTYPE,JP,0,
     &             MPI_COMM_EDDY,STATUS,IERR)
              DO K=KZ1,KZ2+(JP+1)/MYSIZE
                WRITE(19) ((DP(I,J,K),I=1,NX),J=1,NY)
              ENDDO
            ENDDO
C            
          ENDIF
          WRITE(19) TIME
          CLOSE(19)

        ELSE
          CALL MPI_SEND(P(1,1,1),NX*NY*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        ENDIF
      ELSEIF(DIR==2) THEN
c     note that all files write to output
        OPEN(19,FILE=NAME//PROC(MYRANK),
     $       STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(19) NX,NY,NZ,1
        DO K=1,NZ
          WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        CLOSE(19)
      ENDIF
c
      return
      end





c
c     ***************************************************************
c
C     subroutine for I/O - read or write a real 3D array in ASCI format
c
c     ***************************************************************
c
      SUBROUTINE IOSCALAR5(NAME,P,NX,NY,NZ,DIR,TIME)
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
      CHARACTER NAME*(*)
      INTEGER DIR,NX,NY,NZ
      REAL P(NX,NY,NZ),DP(NX,NY,NZ)
c      REAL*4 P1(NX,NY,NZ)
      REAL TIME
C
      INTEGER I,J,K,JP,STATUS(MPI_STATUS_SIZE)

C     IF DIR = -1 READ ELSE WRITE

      IF(DIR==-1) THEN
        OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='FORMATTED')
        READ(19,*) I,J,K,JP
        IF(MYSIZE/=1) THEN
          DO K=1,MYRANK*(NZ-2)
            READ(19,*)
          ENDDO
        ENDIF
        DO K=1,NZ
          READ(19,*) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        DO K=1,(MYSIZE-(MYRANK+1))*(NZ-2)
          READ(19,*)
        ENDDO
        READ(19,*) TIME
        CLOSE(19)

      ELSEIF(DIR==1) THEN
c     note that all files read the input, but only root writes to output
        IF(MYRANK==0) THEN
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='FORMATTED')
          WRITE(19,*) NX,NY,MYSIZE*(NZ-2)+2,1
          IF(MYSIZE==1) THEN
            DO K=1,NZ
              DO J=1,NY
                DO I=1,NX
                  WRITE(19,*) I,J,K,P(I,J,K)
                ENDDO
              ENDDO
            ENDDO
          ELSE            
C
            DO K=1,KZ2
              DO J=1,NY
                DO I=1,NX
                  WRITE(19,*) I,J,K,P(I,J,K)
                ENDDO
              ENDDO
            ENDDO

            DO JP=1,MYSIZE-1
              CALL MPI_RECV(DP(1,1,1),NX*NY*NZ,MTYPE,JP,0,
     &             MPI_COMM_EDDY,STATUS,IERR)
              DO K=KZ1,KZ2+(JP+1)/MYSIZE
                DO J=1,NY
                  DO I=1,NX
                    WRITE(19,*) I,J,K+(NZ-2)*JP,DP(I,J,K)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
C            
          ENDIF
          WRITE(19,*) TIME
          CLOSE(19)

        ELSE
          CALL MPI_SEND(P(1,1,1),NX*NY*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        ENDIF
      ELSEIF(DIR==2) THEN
c     note that all files write to output
        OPEN(19,FILE=NAME//PROC(MYRANK),
     $       STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(19) NX,NY,NZ,1
        DO K=1,NZ
          WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        CLOSE(19)
      ENDIF
c
      return
      end

c
c     ***************************************************************
c
C     subroutine for I/O - read or write an integer 3D array in ASCI format.
c
c     ***************************************************************
c
      SUBROUTINE IOSCALAR6(NAME,P,NX,NY,NZ,DIR)
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
      CHARACTER NAME*(*)
      INTEGER DIR,NX,NY,NZ
      INTEGER P(NX,NY,NZ),DP(NX,NY,NZ)
c      REAL*4 P1(NX,NY,NZ)
      REAL TIME
C
      INTEGER I,J,K,JP,STATUS(MPI_STATUS_SIZE)

C     IF DIR = -1 READ ELSE WRITE

      IF(DIR==-1) THEN
        OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='FORMATTED')
        READ(19,*) I,J,K,JP
        IF(MYSIZE/=1) THEN
          DO K=1,MYRANK*(NZ-2)
            READ(19,*)
          ENDDO
        ENDIF
        DO K=1,NZ
          READ(19,*) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        DO K=1,(MYSIZE-(MYRANK+1))*(NZ-2)
          READ(19,*)
        ENDDO
        READ(19,*) TIME
        CLOSE(19)

      ELSEIF(DIR==1) THEN
c     note that all files read the input, but only root writes to output
        IF(MYRANK==0) THEN
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='FORMATTED')
          WRITE(19,*) NX,NY,MYSIZE*(NZ-2)+2,1
          IF(MYSIZE==1) THEN
            DO K=1,NZ
              DO J=1,NY
                DO I=1,NX
c                  WRITE(19,*) I,J,K,REAL(P(I,J,K),4)
                  WRITE(19,*) I,J,K,P(I,J,K)
                ENDDO
              ENDDO
            ENDDO
          ELSE            
C
            DO K=1,KZ2
              DO J=1,NY
                DO I=1,NX
c                  WRITE(19,*) I,J,K,REAL(P(I,J,K),4)
                  WRITE(19,*) I,J,K,P(I,J,K)
                ENDDO
              ENDDO
            ENDDO

            DO JP=1,MYSIZE-1
              CALL MPI_RECV(DP(1,1,1),NX*NY*NZ,MPI_INTEGER,JP,0,
     &             MPI_COMM_EDDY,STATUS,IERR)
              DO K=KZ1,KZ2+(JP+1)/MYSIZE
                DO J=1,NY
                  DO I=1,NX
c                    WRITE(19,*) I,J,K+(NZ-2)*JP,REAL(DP(I,J,K),4)
                    WRITE(19,*) I,J,K+(NZ-2)*JP,DP(I,J,K)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
C            
          ENDIF
c          WRITE(19,*) TIME
          CLOSE(19)

        ELSE
          CALL MPI_SEND(P(1,1,1),NX*NY*NZ,MPI_INTEGER,0,0,MPI_COMM_EDDY,IERR)
        ENDIF
      ELSEIF(DIR==2) THEN
c     note that all files write to output
        OPEN(19,FILE=NAME//PROC(MYRANK),
     $       STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(19) NX,NY,NZ,1
        DO K=1,NZ
          WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        CLOSE(19)
      ENDIF
c
      return
      end

C---- subroutine ioshearvertexc -----------------N. Beratlis-11 May 2009-
C
C     PURPOSE: Write shear stresses on center of triangles.
C
C------------------------------------------------------------------------
      subroutine ioshearvertexc(filename,nbd,trino,pbd,dudxb,dudyb,dudzb
     &     ,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb,cf_aux,mrkp,mrks,nfacet,nfacetot,ilb,ile)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer nbd,nfacet,nfacetot,ilb,ile
      integer mrkp(nfacet),mrks(9,nfacet)
      integer trino(nfacet)
c      real    vertexc(3,nfacet)
      real    pbd(nfacet)
      real    dudxb(nfacet),dudyb(nfacet),dudzb(nfacet)
     &       ,dvdxb(nfacet),dvdyb(nfacet),dvdzb(nfacet)
     &       ,dwdxb(nfacet),dwdyb(nfacet),dwdzb(nfacet)
     &       ,cf_aux(nfacet,6),ut11(nfacet),ut12(nfacet)
     &       ,ut13(nfacet),dF1(nfacet),dF2(nfacet),uF1dF1(nfacet)

      character*(*) filename
c
c.... Local arrays
      integer nfacetmax

      CALL MPI_ALLREDUCE(nfacet,nfacetmax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_EDDY,IERR)

      call ioscaltrinobd('pres.'//trim(filename),trino,pbd,mrkp,nfacet,nfacetmax,nfacetot,ilb,ile)
!       call ioscaltrinobd('pres.'//trim(filename),trino,pbd,mrkp,nfacet,nfacetmax,nfacetot,nbd)
      call ioscaltrinobd('dudx.'//trim(filename),trino,dudxb,mrks(1,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)
      call ioscaltrinobd('dudy.'//trim(filename),trino,dudyb,mrks(2,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)
      call ioscaltrinobd('dudz.'//trim(filename),trino,dudzb,mrks(3,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)
      call ioscaltrinobd('dvdx.'//trim(filename),trino,dvdxb,mrks(4,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)
      call ioscaltrinobd('dvdy.'//trim(filename),trino,dvdyb,mrks(5,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)
      call ioscaltrinobd('dvdz.'//trim(filename),trino,dvdzb,mrks(6,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)
      call ioscaltrinobd('dwdx.'//trim(filename),trino,dwdxb,mrks(7,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)
      call ioscaltrinobd('dwdy.'//trim(filename),trino,dwdyb,mrks(8,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)
      call ioscaltrinobd('dwdz.'//trim(filename),trino,dwdzb,mrks(9,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)


  
      ut11   = cf_aux(:,1)
      ut12   = cf_aux(:,2)
      ut13   = cf_aux(:,3)
      dF1    = cf_aux(:,4)
      dF2    = cf_aux(:,5) 
      uF1dF1 = cf_aux(:,6)



      call ioscaltrinobd('ut11.'//trim(filename),trino,ut11,mrks(4,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)
      call ioscaltrinobd('ut12.'//trim(filename),trino,ut12,mrks(5,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)
      call ioscaltrinobd('ut13.'//trim(filename),trino,ut13,mrks(6,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)
      call ioscaltrinobd('dF1.'//trim(filename),trino,dF1,mrks(7,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)
      call ioscaltrinobd('dF2.'//trim(filename),trino,dF2,mrks(8,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)
      call ioscaltrinobd('uF1dF1.'//trim(filename),trino,uF1dF1,mrks(9,1:nfacet),nfacet,nfacetmax,nfacetot,ilb,ile)

      
      return

      end

C------------------------------------------------------------------------

C---- subroutine ioshearvertex -----------------N. Beratlis-11 May 2009-
C
C     PURPOSE: Write shear stresses on vertices of triangles.
C
C------------------------------------------------------------------------
      subroutine ioshearvertex(filename,nbd,trino,pbd,dudxb,dudyb,dudzb
     &     ,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb,mrkp,mrks,nfacet,nfacetot,ilb,ile)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer nbd,nfacet,nfacetot,ilb,ile
      integer mrkp(3,nfacet),mrks(9,3,nfacet)
      integer trino(nfacet)
      real    pbd(3,nfacet)
      real    dudxb(3,nfacet),dudyb(3,nfacet),dudzb(3,nfacet)
     &       ,dvdxb(3,nfacet),dvdyb(3,nfacet),dvdzb(3,nfacet)
     &       ,dwdxb(3,nfacet),dwdyb(3,nfacet),dwdzb(3,nfacet)
      character*(*) filename
c
c.... Local arrays
      integer nfacetmax

      CALL MPI_ALLREDUCE(nfacet,nfacetmax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_EDDY,IERR)

c      call ioscaltrinovtx('pres.'//trim(filename),trino,pbd,mrkp,nfacet,nfacetmax,nfacetot,ilb,ile)
c      call ioscaltrinobd('dudx.'//trim(filename),trino,dudxb,mrks(1,1),nfacet,nfacetmax,nfacetot,ilb,ile)
c      call ioscaltrinobd('dudy.'//trim(filename),trino,dudyb,mrks(1,2),nfacet,nfacetmax,nfacetot,ilb,ile)
c      call ioscaltrinobd('dudz.'//trim(filename),trino,dudzb,mrks(1,3),nfacet,nfacetmax,nfacetot,ilb,ile)
c      call ioscaltrinobd('dvdx.'//trim(filename),trino,dvdxb,mrks(1,4),nfacet,nfacetmax,nfacetot,ilb,ile)
c      call ioscaltrinobd('dvdy.'//trim(filename),trino,dvdyb,mrks(1,5),nfacet,nfacetmax,nfacetot,ilb,ile)
c      call ioscaltrinobd('dvdz.'//trim(filename),trino,dvdzb,mrks(1,6),nfacet,nfacetmax,nfacetot,ilb,ile)
c      call ioscaltrinobd('dwdx.'//trim(filename),trino,dwdxb,mrks(1,7),nfacet,nfacetmax,nfacetot,ilb,ile)
c      call ioscaltrinobd('dwdy.'//trim(filename),trino,dwdyb,mrks(1,8),nfacet,nfacetmax,nfacetot,ilb,ile)
c      call ioscaltrinobd('dwdz.'//trim(filename),trino,dwdzb,mrks(1,9),nfacet,nfacetmax,nfacetot,ilb,ile)

      return

      end

C------------------------------------------------------------------------



C---- call iotaubd--------------------------N. Beratlis-31 March 2010---
C
C     PURPOSE: Write components of tau located on center of triangel
C
C------------------------------------------------------------------------

!      subroutine iotaubd(filename,vertexc,pbd,mrkp,nfacet,nfacetmax,nfacetot,nbd)
!
!      include 'common.h'
!      include 'immersed.h'
!      include 'mpif.h'
c
c...  Input/Output arrays
!      integer nbd,nfacet,nfacetmax,nfacetot
!      integer mrkp(nfacet)
!      real    vertexc(3,nfacet)
!      real    pbd(nfacet)
!      character*(*) filename
c
c.... Local arrays
!      integer ibd,n,i,j,jp
!      integer STATUS(MPI_STATUS_SIZE)
!      real    tmp1(nfacetmax)
!      real    vertexc1(3,nfacetmax)
!      real    psi
c
c...... Functions
!      real    anglerad
        
!      IF(MYRANK==0) THEN
!        open(unit=20,file=trim(filename),form='formatted')
!        write(20,'(A,A)') 'VARIABLES = "X", "Y", "Z", "VAR"'
!        do ibd=1,nbd
!          write(20,'(A,I8,A)') 'ZONE I=',nfacetot,', DATAPACKING=POINT'
!          do i=lb(ibd)+1,lb(ibd)+mb(ibd)
!            if(mrkp(i)==1) then
!              psi = anglerad(-vertexc(3,i),sqrt(vertexc(1,i)**2.+vertexc(2,i)**2.))
!              write(20,'(12(1X,E18.11))') vertexc(1,i),vertexc(2,i),vertexc(3,i),pbd(i)
!            endif
!          enddo
!
!          IF(MYSIZE>1) THEN
!            
!            DO JP=1,MYSIZE-1
!              CALL MPI_RECV(n,1,MPI_INTEGER,JP,JP,MPI_COMM_EDDY,STATUS,IERR)
!              IF(n>0) THEN
!                CALL MPI_RECV(vertexc1,3*n,MTYPE,JP,1,MPI_COMM_EDDY,STATUS,IERR)
!                CALL MPI_RECV(tmp1,n,MTYPE,JP,2,MPI_COMM_EDDY,STATUS,IERR)
!                do i=1,n
!                  psi = anglerad(-vertexc1(3,i),sqrt(vertexc1(1,i)**2.+vertexc1(2,i)**2.))
!                  write(20,'(12(1X,E18.11))') vertexc1(1,i),vertexc1(2,i),vertexc1(3,i),tmp1(i)
!                enddo             
!              ENDIF
!            ENDDO
!          ENDIF
       ! ENDDO
!        close(20)
!
!      ELSE
!
!        do ibd=1,nbd
!          j=0
!          do i=lb(ibd)+1,lb(ibd)+mb(ibd)
!            if(mrkp(i)==1) then
!              j=j+1
!              vertexc1(:,j)=vertexc(:,i)
!              tmp1(j)=pbd(i)
!            endif
!          enddo
!          n=j
!
!          CALL MPI_SEND(n,1,MPI_INTEGER,0,MYRANK,MPI_COMM_EDDY,IERR)
!          if(n>0) then
!            CALL MPI_SEND(vertexc1,3*n,MTYPE,0,1,MPI_COMM_EDDY,IERR)
!            CALL MPI_SEND(tmp1,n,MTYPE,0,2,MPI_COMM_EDDY,IERR)
!          endif
!        enddo
!      ENDIF
!        
!      return
!
!      end
C------------------------------------------------------------------------

      subroutine iotaubd(filename,vertexc,pbd,mrkp,nfacet,nfacetmax,nfacetot,nbd)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c...  Input/Output arrays
      integer nbd,nfacet,nfacetmax,nfacetot
      integer mrkp(nfacet)
      real    vertexc(3,nfacet)
      real    pbd(nfacet)
      character*(*) filename
c
c.... Local arrays
      integer ibd,n,i,j,jp
      integer STATUS(MPI_STATUS_SIZE)
      real    tmp1(nfacetmax)
      real    vertexc1(3,nfacetmax)
      real    psi
      integer counter(180)
c
c...... Functions
      real    anglerad
        
      IF(MYRANK==0) THEN
        open(unit=20,file=trim(filename),form='formatted')
         write(20,'(A,A)') 'VARIABLES = "X", "Y", "Z", "VAR"'
!          write(20,'(A,A)') 'VARIABLES = "THETA", "VAR"'
        do ibd=1,nbd
          write(20,'(A,I8,A)') 'ZONE I=',nfacetot,', DATAPACKING=POINT'
          do i=lb(ibd)+1,lb(ibd)+mb(ibd)
            if(mrkp(i)==1) then
              psi = anglerad(-vertexc(3,i),sqrt(vertexc(1,i)**2.+vertexc(2,i)**2.))
                write(20,'(12(1X,E18.11))') vertexc(1,i),vertexc(2,i),vertexc(3,i),pbd(i)
!                 write(20,'(12(1X,E18.11))') psi*180.0d0/3.140,pbd(i)
            endif
          enddo

          IF(MYSIZE>1) THEN

            DO JP=1,MYSIZE-1
              CALL MPI_RECV(n,1,MPI_INTEGER,JP,JP,MPI_COMM_EDDY,STATUS,IERR)
              IF(n>0) THEN
                CALL MPI_RECV(vertexc1,3*n,MTYPE,JP,1,MPI_COMM_EDDY,STATUS,IERR)
                CALL MPI_RECV(tmp1,n,MTYPE,JP,2,MPI_COMM_EDDY,STATUS,IERR)
                do i=1,n
                  psi = anglerad(-vertexc1(3,i),sqrt(vertexc1(1,i)**2.+vertexc1(2,i)**2.))
                  write(20,'(12(1X,E18.11))') vertexc1(1,i),vertexc1(2,i),vertexc1(3,i),tmp1(i)
!                   write(20,'(12(1X,E18.11))') psi*180.0d0/3.140,tmp1(i)
                enddo             
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        close(20)

      ELSE

        do ibd=1,nbd
          j=0
          do i=lb(ibd)+1,lb(ibd)+mb(ibd)
            if(mrkp(i)==1) then
              j=j+1
              vertexc1(:,j)=vertexc(:,i)
              tmp1(j)=pbd(i)
            endif
          enddo
          n=j

          CALL MPI_SEND(n,1,MPI_INTEGER,0,MYRANK,MPI_COMM_EDDY,IERR)
          if(n>0) then
            CALL MPI_SEND(vertexc1,3*n,MTYPE,0,1,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(tmp1,n,MTYPE,0,2,MPI_COMM_EDDY,IERR)
          endif
        enddo
      ENDIF
        
      return

      end




!
C---- subroutine iopressnodes-------------------N. Beratlis-11 May 2009-
C
C     PURPOSE: Write shear stresses on center of triangles.
C
C------------------------------------------------------------------------
      subroutine iopressnodes(filename,node,pbd,nnodes,ilb,ile)

      implicit none
      include 'immersed.h'
c
c.... Input/Output Arrays
      integer nnodes,ilb,ile
      real    node(3,nnodes),pbd(nnodes)
      character*(*) filename
c
c.... Local arrays
      integer i,ibd,n
      real    psi,r
c
c.... Functions
      real    anglerad

      n = ile-ilb+1

      open(unit=20,file=trim(filename),form='formatted')
c      write(20,*) n
      write(20,'(A,A)') 'VARIABLES = "ANGLE", "P"'
      write(20,'(A,I8,A)') 'ZONE I=',n,', DATAPACKING=POINT'
      do i=ilb,ile
        r = sqrt(node(1,i)**2.+node(2,i)**2.)
        psi = anglerad(-node(3,i),r)
c        write(20,'(12(1X,F14.8))') node(1,i),node(2,i),node(3,i),psi,pbd(i)
        write(20,'(12(1X,F14.8))') psi,pbd(i)
      enddo
      close(20)

      return

      end

C------------------------------------------------------------------------

c     ***************************************************************
c
C     subroutine for I/O - read or write immersed boundary arrays.
C     Format of file is read my TECPLOT
c
c     ***************************************************************
c
      SUBROUTINE IOMFORC_CAR(FNAME,XU,YU,ZU,NX,NY,NZ,IU,JU,KU,MRKU
     &      ,XNU,YNU,ZNU,NXU,NYU,NZU,LIMU,MIMU,IDIR)

      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'mpif.h'

      CHARACTER FNAME*(*)
      INTEGER LIMU,MIMU,NX,NY,NZ,IDIR
      INTEGER IU(NFCMAX),JU(NFCMAX),KU(NFCMAX),IMRKU(NFCMAX)
      INTEGER MRKU(NFCMAX)
      REAL    XNU(NFCMAX),YNU(NFCMAX),ZNU(NFCMAX)
      REAL    NXU(NFCMAX),NYU(NFCMAX),NZU(NFCMAX)
      REAL    XU(NX,NY),YU(NX,NY),ZU(NZ)

      INTEGER LIMU1,MIMU1
      INTEGER IU1(NFCMAX),JU1(NFCMAX),KU1(NFCMAX),IMRKU1(NFCMAX)
      INTEGER MRKU1(NFCMAX)
      REAL    XNU1(NFCMAX),YNU1(NFCMAX),ZNU1(NFCMAX)
      REAL    NXU1(NFCMAX),NYU1(NFCMAX),NZU1(NFCMAX)
      REAL    ZU1(NZ)

      INTEGER MIMUG,ICOUNT,IPROC,MU,I,J,K,JR,II,IM,STATUS(MPI_STATUS_SIZE)

      IF(IDIR==1) THEN

        ICOUNT = MIMU
        CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
        WRITE(6,'(3(1X,A,I6))')
     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG

        IF(MYRANK.EQ.0) THEN
          OPEN(20,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
          WRITE(20, '(A)') 'VARIABLES = "X" "Y" "Z"'
          DO IM=LIMU+1,LIMU+MIMU
            WRITE(20, '(A)') 'ZONE I=2 DATAPACKING=POINT'
            WRITE(20, '(3(1X,F10.6))') sqrt(xnu(im)**2.0 + ynu(im)**2.0),real(ju(im),8),znu(im)
            IF(mrku(im)==1) THEN
              write(20, '(3(1X,F10.6))') sqrt(nxu(im)**2.0 + nyu(im)**2.0),real(ju(im),8),nzu(im)
            ELSE
              i = iu(im)+int(nxu(im))
              j = ju(im)+int(nyu(im))
              k = ku(im)+int(nzu(im))
              if(ju(im)==351) then
                 write(6,*) iu(im),ju(im),ku(im),i,j,k,int(nxu(im))
     $                ,int(nyu(im)),int(nzu(im)),sqrt(xnu(im)**2.0 +
     $                ynu(im)**2.0),sqrt(xu(i,j)**2.0 + yu(i,j)**2.0)
              endif
              write(20, '(3(1X,F10.6))') sqrt(xu(i,j)**2.0 + yu(i,j)**2.0),real(ju(im),8),zu(k)
            ENDIF
          ENDDO

          DO IPROC=1,MYSIZE-1
            CALL MPI_RECV(MU,1,MPI_INTEGER,IPROC,1,MPI_COMM_EDDY,STATUS,IERR)
            IF(MU/=0) THEN
              CALL MPI_RECV(IU1(1),MU,MPI_INTEGER,IPROC,2,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(JU1(1),MU,MPI_INTEGER,IPROC,3,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(KU1(1),MU,MPI_INTEGER,IPROC,4,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(MRKU1(1),MU,MPI_INTEGER,IPROC,5,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(XNU1(1),MU,MTYPE,IPROC,12,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(YNU1(1),MU,MTYPE,IPROC,13,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZNU1(1),MU,MTYPE,IPROC,14,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NXU1(1),MU,MTYPE,IPROC,15,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NYU1(1),MU,MTYPE,IPROC,16,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NZU1(1),MU,MTYPE,IPROC,17,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZU1(1),NZ,MTYPE,IPROC,18,MPI_COMM_EDDY,STATUS,IERR)

              DO IM=1,MU
                WRITE(20, '(A)') 'ZONE I=2 DATAPACKING=POINT'
                WRITE(20, '(3(1X,F10.6))') sqrt(XNU1(IM)**2.0 + YNU1(IM)**2.0),JU1(IM),ZNU1(IM)
                IF(mrku1(im)==1) THEN
                  write(20, '(3(1X,F10.6))') sqrt(nxu1(im)**2.0 + nyu1(im)**2.0),JU1(IM),nzu1(im)
                ELSE
                  i = iu1(im)+int(nxu1(im))
                  j = ju1(im)+int(nyu1(im))
                  k = ku1(im)+int(nzu1(im))
                  write(20, '(3(1X,F10.6))') sqrt(xu(i,j)**2.0 + yu(i,j)**2.0),JU1(IM),zu1(k)
                ENDIF
              ENDDO

            ENDIF
          ENDDO
          CLOSE(20)

          ELSE
            CALL MPI_SEND(MIMU,1,MPI_INTEGER,0,1,MPI_COMM_EDDY,IERR)
            IF(MIMU/=0) THEN
              IU1(1:MIMU) = IU(LIMU+1:LIMU+MIMU)
              JU1(1:MIMU) = JU(LIMU+1:LIMU+MIMU)
              KU1(1:MIMU) = KU(LIMU+1:LIMU+MIMU)! + (NZ-2)*MYRANK
              MRKU1(1:MIMU) = MRKU(LIMU+1:LIMU+MIMU)
              XNU1(1:MIMU) = XNU(LIMU+1:LIMU+MIMU)
              YNU1(1:MIMU) = YNU(LIMU+1:LIMU+MIMU)
              ZNU1(1:MIMU) = ZNU(LIMU+1:LIMU+MIMU)
              NXU1(1:MIMU) = NXU(LIMU+1:LIMU+MIMU)
              NYU1(1:MIMU) = NYU(LIMU+1:LIMU+MIMU)
              NZU1(1:MIMU) = NZU(LIMU+1:LIMU+MIMU)

              CALL MPI_SEND(IU1(1),MIMU,MPI_INTEGER,0,2,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(JU1(1),MIMU,MPI_INTEGER,0,3,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(KU1(1),MIMU,MPI_INTEGER,0,4,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(MRKU1(1),MIMU,MPI_INTEGER,0,5,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(XNU1(1),MIMU,MTYPE,0,12,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(YNU1(1),MIMU,MTYPE,0,13,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZNU1(1),MIMU,MTYPE,0,14,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NXU1(1),MIMU,MTYPE,0,15,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NYU1(1),MIMU,MTYPE,0,16,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NZU1(1),MIMU,MTYPE,0,17,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZU(1),NZ,MTYPE,0,18,MPI_COMM_EDDY,IERR)
            ENDIF

          ENDIF

          ICOUNT = MIMU
          CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
c          WRITE(6,'(3(1X,A,I6))')
c     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG
        ENDIF

        RETURN
 120    FORMAT(17(1X,I4),23(1X,F14.8))
c 120  FORMAT(9(1X,I4),11(1X,F12.5))
 122    FORMAT(13(1X,I8),17(1X,F14.8))

        END



c     ***************************************************************
c
C     subroutine for I/O - read or write immersed boundary arrays.
C     Format of file is read my TECPLOT
c
c     ***************************************************************
c
      SUBROUTINE IOMFORC_CYL(FNAME,XU,YU,ZU,NX,NY,NZ,IU,JU,KU,MRKU
     &      ,XNU,YNU,ZNU,NXU,NYU,NZU,LIMU,MIMU,IDIR)

      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'mpif.h'


      CHARACTER FNAME*(*)
      INTEGER LIMU,MIMU,NX,NY,NZ,IDIR
      INTEGER IU(NFCMAX),JU(NFCMAX),KU(NFCMAX),IMRKU(NFCMAX)
      INTEGER MRKU(NFCMAX)
      REAL    XNU(NFCMAX),YNU(NFCMAX),ZNU(NFCMAX)
      REAL    NXU(NFCMAX),NYU(NFCMAX),NZU(NFCMAX)
      REAL    XU(NX,NY),YU(NX,NY),ZU(NZ)

      INTEGER LIMU1,MIMU1
      INTEGER IU1(NFCMAX),JU1(NFCMAX),KU1(NFCMAX),IMRKU1(NFCMAX)
      INTEGER MRKU1(NFCMAX)
      REAL    XNU1(NFCMAX),YNU1(NFCMAX),ZNU1(NFCMAX)
      REAL    NXU1(NFCMAX),NYU1(NFCMAX),NZU1(NFCMAX)
      REAL    ZU1(NZ)

      INTEGER MIMUG,ICOUNT,IPROC,MU,I,J,K,JR,II,IM,STATUS(MPI_STATUS_SIZE)

      IF(IDIR==1) THEN

        ICOUNT = MIMU
        CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
        WRITE(6,'(3(1X,A,I6))')
     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG


        IF(MYRANK.EQ.0) THEN
          OPEN(20,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
          WRITE(20, '(A)') 'VARIABLES = "X" "Y" "Z"'
          DO IM=LIMU+1,LIMU+MIMU
            WRITE(20, '(A)') 'ZONE I=2 DATAPACKING=POINT'
            WRITE(20, '(3(1X,F10.6))') xnu(im),ynu(im),znu(im)
            IF(mrku(im)==1) THEN
              write(20, '(3(1X,F10.6))') nxu(im),nyu(im),nzu(im)
            ELSE
              i = iu(im)+int(nxu(im))
              j = ju(im)+int(nyu(im))
              k = ku(im)+int(nzu(im))
              write(20, '(3(1X,F10.6))') xu(i,j),yu(i,j),zu(k)
            ENDIF
          ENDDO

          DO IPROC=1,MYSIZE-1
            CALL MPI_RECV(MU,1,MPI_INTEGER,IPROC,1,MPI_COMM_EDDY,STATUS,IERR)
            IF(MU/=0) THEN
              CALL MPI_RECV(IU1(1),MU,MPI_INTEGER,IPROC,2,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(JU1(1),MU,MPI_INTEGER,IPROC,3,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(KU1(1),MU,MPI_INTEGER,IPROC,4,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(MRKU1(1),MU,MPI_INTEGER,IPROC,5,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(XNU1(1),MU,MTYPE,IPROC,12,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(YNU1(1),MU,MTYPE,IPROC,13,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZNU1(1),MU,MTYPE,IPROC,14,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NXU1(1),MU,MTYPE,IPROC,15,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NYU1(1),MU,MTYPE,IPROC,16,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NZU1(1),MU,MTYPE,IPROC,17,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZU1(1),NZ,MTYPE,IPROC,18,MPI_COMM_EDDY,STATUS,IERR)

              DO IM=1,MU
                WRITE(20, '(A)') 'ZONE I=2 DATAPACKING=POINT'
                WRITE(20, '(3(1X,F10.6))') XNU1(IM),YNU1(IM),ZNU1(IM)
                IF(mrku1(im)==1) THEN
                  write(20, '(3(1X,F10.6))') nxu1(im),nyu1(im),nzu1(im)
                ELSE
                  i = iu1(im)+int(nxu1(im))
                  j = ju1(im)+int(nyu1(im))
                  k = ku1(im)+int(nzu1(im))
                  write(20, '(3(1X,F10.6))') xu(i,j),yu(i,j),zu1(k)
                ENDIF
              ENDDO

            ENDIF
          ENDDO
          CLOSE(20)

          ELSE
            CALL MPI_SEND(MIMU,1,MPI_INTEGER,0,1,MPI_COMM_EDDY,IERR)
            IF(MIMU/=0) THEN
              IU1(1:MIMU) = IU(LIMU+1:LIMU+MIMU)
              JU1(1:MIMU) = JU(LIMU+1:LIMU+MIMU)
              KU1(1:MIMU) = KU(LIMU+1:LIMU+MIMU)! + (NZ-2)*MYRANK
              MRKU1(1:MIMU) = MRKU(LIMU+1:LIMU+MIMU)
              XNU1(1:MIMU) = XNU(LIMU+1:LIMU+MIMU)
              YNU1(1:MIMU) = YNU(LIMU+1:LIMU+MIMU)
              ZNU1(1:MIMU) = ZNU(LIMU+1:LIMU+MIMU)
              NXU1(1:MIMU) = NXU(LIMU+1:LIMU+MIMU)
              NYU1(1:MIMU) = NYU(LIMU+1:LIMU+MIMU)
              NZU1(1:MIMU) = NZU(LIMU+1:LIMU+MIMU)

              CALL MPI_SEND(IU1(1),MIMU,MPI_INTEGER,0,2,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(JU1(1),MIMU,MPI_INTEGER,0,3,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(KU1(1),MIMU,MPI_INTEGER,0,4,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(MRKU1(1),MIMU,MPI_INTEGER,0,5,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(XNU1(1),MIMU,MTYPE,0,12,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(YNU1(1),MIMU,MTYPE,0,13,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZNU1(1),MIMU,MTYPE,0,14,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NXU1(1),MIMU,MTYPE,0,15,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NYU1(1),MIMU,MTYPE,0,16,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NZU1(1),MIMU,MTYPE,0,17,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZU(1),NZ,MTYPE,0,18,MPI_COMM_EDDY,IERR)
            ENDIF

          ENDIF

          ICOUNT = MIMU
          CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
c          WRITE(6,'(3(1X,A,I6))')
c     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG
        ENDIF

        RETURN
 120    FORMAT(17(1X,I4),23(1X,F14.8))
c 120  FORMAT(9(1X,I4),11(1X,F12.5))
 122    FORMAT(13(1X,I8),17(1X,F14.8))

        END

c     ***************************************************************
c
C     subroutine for I/O - read or write immersed boundary arrays.
C     Format of file is read my TECPLOT
c
c     ***************************************************************
c
      SUBROUTINE IOMFORC_CT_CYL(FNAME,XU,YU,ZU,NX,NY,NZ,IU,JU,KU,MRKU
     &      ,XNU,YNU,ZNU,NXU,NYU,NZU,LIMU,MIMU,IDIR)

      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'mpif.h'


      CHARACTER FNAME*(*)
      INTEGER LIMU,MIMU,NX,NY,NZ,IDIR
      INTEGER IU(NFCMAX),JU(NFCMAX),KU(NFCMAX),IMRKU(NFCMAX)
      INTEGER MRKU(NFCMAX)
      REAL    XNU(NFCMAX),YNU(NFCMAX),ZNU(NFCMAX)
      REAL    NXU(NFCMAX),NYU(NFCMAX),NZU(NFCMAX)
      REAL    XU(NX,NY),YU(NX,NY),ZU(NZ)

      INTEGER LIMU1,MIMU1
      INTEGER IU1(NFCMAX),JU1(NFCMAX),KU1(NFCMAX),IMRKU1(NFCMAX)
      INTEGER MRKU1(NFCMAX)
      REAL    XNU1(NFCMAX),YNU1(NFCMAX),ZNU1(NFCMAX)
      REAL    NXU1(NFCMAX),NYU1(NFCMAX),NZU1(NFCMAX)
      REAL    ZU1(NZ)

      INTEGER MIMUG,ICOUNT,IPROC,MU,I,J,K,JR,II,IM,STATUS(MPI_STATUS_SIZE)

      IF(IDIR==1) THEN

        ICOUNT = COUNT(MRKU(LIMU+1:LIMU+MIMU)==3)
        CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
        WRITE(6,'(3(1X,A,I6))')
     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG


        IF(MYRANK.EQ.0) THEN
          OPEN(20,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
          WRITE(20,'(A)') 'VARIABLES = "X" "Y" "Z"'
          DO IM=LIMU+1,LIMU+MIMU
            IF(mrku(im)==3) THEN
              WRITE(20,'(A)') 'ZONE I=2 DATAPACKING=POINT'
              WRITE(20,'(3(1X,F10.6))') xnu(im),ynu(im),znu(im)
              write(20,'(3(1X,F10.6))') nxu(im),nyu(im),nzu(im)
            ENDIF
          ENDDO

          DO IPROC=1,MYSIZE-1
            CALL MPI_RECV(MU,1,MPI_INTEGER,IPROC,1,MPI_COMM_EDDY,STATUS,IERR)
            IF(MU/=0) THEN
              CALL MPI_RECV(IU1(1),MU,MPI_INTEGER,IPROC,2,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(JU1(1),MU,MPI_INTEGER,IPROC,3,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(KU1(1),MU,MPI_INTEGER,IPROC,4,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(MRKU1(1),MU,MPI_INTEGER,IPROC,5,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(XNU1(1),MU,MTYPE,IPROC,12,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(YNU1(1),MU,MTYPE,IPROC,13,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZNU1(1),MU,MTYPE,IPROC,14,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NXU1(1),MU,MTYPE,IPROC,15,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NYU1(1),MU,MTYPE,IPROC,16,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NZU1(1),MU,MTYPE,IPROC,17,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZU1(1),NZ,MTYPE,IPROC,18,MPI_COMM_EDDY,STATUS,IERR)

              DO IM=1,MU
                IF(mrku1(im)==3) THEN
                  WRITE(20,'(A)') 'ZONE I=2 DATAPACKING=POINT'
                  WRITE(20,'(3(1X,F10.6))') XNU1(IM),YNU1(IM),ZNU1(IM)
                  write(20,'(3(1X,F10.6))') nxu1(im),nyu1(im),nzu1(im)
                ENDIF
              ENDDO

            ENDIF
          ENDDO
          CLOSE(20)

          ELSE
            CALL MPI_SEND(MIMU,1,MPI_INTEGER,0,1,MPI_COMM_EDDY,IERR)
            IF(MIMU/=0) THEN
              IU1(1:MIMU) = IU(LIMU+1:LIMU+MIMU)
              JU1(1:MIMU) = JU(LIMU+1:LIMU+MIMU)
              KU1(1:MIMU) = KU(LIMU+1:LIMU+MIMU)! + (NZ-2)*MYRANK
              MRKU1(1:MIMU) = MRKU(LIMU+1:LIMU+MIMU)
              XNU1(1:MIMU) = XNU(LIMU+1:LIMU+MIMU)
              YNU1(1:MIMU) = YNU(LIMU+1:LIMU+MIMU)
              ZNU1(1:MIMU) = ZNU(LIMU+1:LIMU+MIMU)
              NXU1(1:MIMU) = NXU(LIMU+1:LIMU+MIMU)
              NYU1(1:MIMU) = NYU(LIMU+1:LIMU+MIMU)
              NZU1(1:MIMU) = NZU(LIMU+1:LIMU+MIMU)

              CALL MPI_SEND(IU1(1),MIMU,MPI_INTEGER,0,2,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(JU1(1),MIMU,MPI_INTEGER,0,3,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(KU1(1),MIMU,MPI_INTEGER,0,4,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(MRKU1(1),MIMU,MPI_INTEGER,0,5,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(XNU1(1),MIMU,MTYPE,0,12,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(YNU1(1),MIMU,MTYPE,0,13,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZNU1(1),MIMU,MTYPE,0,14,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NXU1(1),MIMU,MTYPE,0,15,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NYU1(1),MIMU,MTYPE,0,16,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NZU1(1),MIMU,MTYPE,0,17,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZU(1),NZ,MTYPE,0,18,MPI_COMM_EDDY,IERR)
            ENDIF

          ENDIF

          ICOUNT = MIMU
          CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
c          WRITE(6,'(3(1X,A,I6))')
c     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG
        ENDIF

        RETURN
 120    FORMAT(17(1X,I4),23(1X,F14.8))
c 120  FORMAT(9(1X,I4),11(1X,F12.5))
 122    FORMAT(13(1X,I8),17(1X,F14.8))

        END

c     ***************************************************************
c
C     subroutine for I/O - read or write immersed boundary arrays.
C     Format of file is read my TECPLOT
c
c     ***************************************************************
c
      SUBROUTINE IOMFORC_WM_CAR(FNAME,XU,YU,ZU,NX,NY,NZ,IU,JU,KU,MRKU,XNU,YNU
     $     ,ZNU,NXU,NYU,NZU,UNVECT,NFACET,LIMU,MIMU,IDIR)

      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'mpif.h'
c
c.... Input/Output arrays
      CHARACTER FNAME*(*)
      INTEGER LIMU,MIMU,NX,NY,NZ,IDIR,NFACET
      INTEGER IU(NFCMAX),JU(NFCMAX),KU(NFCMAX),IMRKU(NFCMAX)
      INTEGER MRKU(NFCMAX)
      REAL    UNVECT(3,NFACET)
c      INTEGER UINDX(NSP,NFCMAX)
      REAL    XNU(NFCMAX),YNU(NFCMAX),ZNU(NFCMAX)
      REAL    NXU(NFCMAX),NYU(NFCMAX),NZU(NFCMAX)
      REAL    XU(NX,NY),YU(NX,NY),ZU(NZ)
c
c.... Local arrays
      INTEGER LIMU1,MIMU1,ITR
      REAL    DCELL
      INTEGER IU1(NFCMAX),JU1(NFCMAX),KU1(NFCMAX),IMRKU1(NFCMAX)
      INTEGER MRKU1(NFCMAX)
      REAL    XNU1(NFCMAX),YNU1(NFCMAX),ZNU1(NFCMAX)
      REAL    NXU1(NFCMAX),NYU1(NFCMAX),NZU1(NFCMAX)
      REAL    ZU1(NZ)

      INTEGER MIMUG,ICOUNT,IPROC,MU,I,J,K,JR,II,IM,STATUS(MPI_STATUS_SIZE)


      IF(IDIR==1) THEN

        ICOUNT = MIMU
        CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
        WRITE(6,'(3(1X,A,I6))')
     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG

        IF(MYRANK.EQ.0) THEN
          OPEN(20,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
          WRITE(20, '(A)') 'VARIABLES = "X" "Y" "Z"'
          DO IM=LIMU+1,LIMU+MIMU
            WRITE(20, '(A)') 'ZONE I=2 DATAPACKING=POINT'
            write(20,'(3(1x,F10.6))') sqrt(xnu(im)**2.+ynu(im)**2.),real(ju(im),8),znu(im)
            itr = int(nyu(im))
            dcell = nxu(im)
            write(20, '(3(1X,F10.6))') sqrt((xnu(im)+uext1*dcell*unvect(1,itr))**2.
     &       + (ynu(im)+uext1*dcell*unvect(2,itr))**2.),real(ju(im),8),znu(im)+uext1*dcell*unvect(3,itr)
          ENDDO

          DO IPROC=1,MYSIZE-1
            CALL MPI_RECV(MU,1,MPI_INTEGER,IPROC,1,MPI_COMM_EDDY,STATUS,IERR)
            IF(MU/=0) THEN
              CALL MPI_RECV(IU1(1),MU,MPI_INTEGER,IPROC,2,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(JU1(1),MU,MPI_INTEGER,IPROC,3,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(KU1(1),MU,MPI_INTEGER,IPROC,4,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(MRKU1(1),MU,MPI_INTEGER,IPROC,10,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(XNU1(1),MU,MTYPE,IPROC,12,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(YNU1(1),MU,MTYPE,IPROC,13,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZNU1(1),MU,MTYPE,IPROC,14,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NXU1(1),MU,MTYPE,IPROC,15,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NYU1(1),MU,MTYPE,IPROC,16,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NZU1(1),MU,MTYPE,IPROC,17,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZU1(1),NZ,MTYPE,IPROC,18,MPI_COMM_EDDY,STATUS,IERR)

              DO IM=1,MU
                WRITE(20, '(A)') 'ZONE I=2 DATAPACKING=POINT'
                write(20,'(3(1x,F10.6))') sqrt(xnu1(im)**2.+ynu1(im)**2.),real(ju(im),8),znu1(im)
                itr = int(nyu1(im))
                dcell = nxu1(im)
                write(20, '(3(1X,F10.6))') sqrt((xnu1(im)+uext1*dcell*unvect(1,itr))**2.
     &                + (ynu1(im)+uext1*dcell*unvect(2,itr))**2.),real(ju(im),8),znu1(im)+uext1*dcell*unvect(3,itr)
              ENDDO

            ENDIF
          ENDDO
          CLOSE(20)

          ELSE
            CALL MPI_SEND(MIMU,1,MPI_INTEGER,0,1,MPI_COMM_EDDY,IERR)
            IF(MIMU/=0) THEN
              IU1(1:MIMU) = IU(LIMU+1:LIMU+MIMU)
              JU1(1:MIMU) = JU(LIMU+1:LIMU+MIMU)
              KU1(1:MIMU) = KU(LIMU+1:LIMU+MIMU)! + (NZ-2)*MYRANK
              MRKU1(1:MIMU) = MRKU(LIMU+1:LIMU+MIMU)
              XNU1(1:MIMU) = XNU(LIMU+1:LIMU+MIMU)
              YNU1(1:MIMU) = YNU(LIMU+1:LIMU+MIMU)
              ZNU1(1:MIMU) = ZNU(LIMU+1:LIMU+MIMU)
              NXU1(1:MIMU) = NXU(LIMU+1:LIMU+MIMU)
              NYU1(1:MIMU) = NYU(LIMU+1:LIMU+MIMU)
              NZU1(1:MIMU) = NZU(LIMU+1:LIMU+MIMU)

              CALL MPI_SEND(IU1(1),MIMU,MPI_INTEGER,0,2,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(JU1(1),MIMU,MPI_INTEGER,0,3,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(KU1(1),MIMU,MPI_INTEGER,0,4,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(MRKU1(1),MIMU,MPI_INTEGER,0,10,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(XNU1(1),MIMU,MTYPE,0,12,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(YNU1(1),MIMU,MTYPE,0,13,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZNU1(1),MIMU,MTYPE,0,14,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NXU1(1),MIMU,MTYPE,0,15,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NYU1(1),MIMU,MTYPE,0,16,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NZU1(1),MIMU,MTYPE,0,17,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZU(1),NZ,MTYPE,0,18,MPI_COMM_EDDY,IERR)
            ENDIF

          ENDIF

          ICOUNT = MIMU
          CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
c          WRITE(6,'(3(1X,A,I6))')
c     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG
        ENDIF

        RETURN
 120    FORMAT(17(1X,I4),23(1X,F14.8))
c 120  FORMAT(9(1X,I4),11(1X,F12.5))
 122    FORMAT(13(1X,I8),17(1X,F14.8))

        END
C------------------------------------------------------------------------


c     ***************************************************************
c
C     subroutine for I/O - read or write immersed boundary arrays.
C     Format of file is read my TECPLOT
c
c     ***************************************************************
c
      SUBROUTINE IOMFORC_WM_CYL(FNAME,XU,YU,ZU,NX,NY,NZ,IU,JU,KU,MRKU,XNU,YNU
     $     ,ZNU,NXU,NYU,NZU,UNVECT,NFACET,LIMU,MIMU,IDIR)

      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'mpif.h'
c
c.... Input/Output arrays
      CHARACTER FNAME*(*)
      INTEGER LIMU,MIMU,NX,NY,NZ,IDIR,NFACET
      INTEGER IU(NFCMAX),JU(NFCMAX),KU(NFCMAX),IMRKU(NFCMAX)
      INTEGER MRKU(NFCMAX)
      REAL    UNVECT(3,NFACET)
c      INTEGER UINDX(NSP,NFCMAX)
      REAL    XNU(NFCMAX),YNU(NFCMAX),ZNU(NFCMAX)
      REAL    NXU(NFCMAX),NYU(NFCMAX),NZU(NFCMAX)
      REAL    XU(NX,NY),YU(NX,NY),ZU(NZ)
c
c.... Local arrays
      INTEGER LIMU1,MIMU1,ITR
      REAL    DCELL
      INTEGER IU1(NFCMAX),JU1(NFCMAX),KU1(NFCMAX),IMRKU1(NFCMAX)
      INTEGER MRKU1(NFCMAX)
      REAL    XNU1(NFCMAX),YNU1(NFCMAX),ZNU1(NFCMAX)
      REAL    NXU1(NFCMAX),NYU1(NFCMAX),NZU1(NFCMAX)
      REAL    ZU1(NZ)

      INTEGER MIMUG,ICOUNT,IPROC,MU,I,J,K,JR,II,IM,STATUS(MPI_STATUS_SIZE)


      IF(IDIR==1) THEN

        ICOUNT = MIMU
        CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
        WRITE(6,'(3(1X,A,I6))')
     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG

        IF(MYRANK.EQ.0) THEN
          OPEN(20,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
          WRITE(20, '(A)') 'VARIABLES = "X" "Y" "Z"'
          DO IM=LIMU+1,LIMU+MIMU
            WRITE(20, '(A)') 'ZONE I=2 DATAPACKING=POINT'
            write(20,'(3(1x,F10.6))') xnu(im),ynu(im),znu(im)
            itr = int(nyu(im))
            dcell = nxu(im)
            write(20, '(3(1X,F10.6))') xnu(im)+uext1*dcell*unvect(1,itr)
     &           ,ynu(im)+uext1*dcell*unvect(2,itr)
     &           ,znu(im)+uext1*dcell*unvect(3,itr)
          ENDDO

          DO IPROC=1,MYSIZE-1
            CALL MPI_RECV(MU,1,MPI_INTEGER,IPROC,1,MPI_COMM_EDDY,STATUS,IERR)
            IF(MU/=0) THEN
              CALL MPI_RECV(IU1(1),MU,MPI_INTEGER,IPROC,2,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(JU1(1),MU,MPI_INTEGER,IPROC,3,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(KU1(1),MU,MPI_INTEGER,IPROC,4,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(MRKU1(1),MU,MPI_INTEGER,IPROC,10,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(XNU1(1),MU,MTYPE,IPROC,12,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(YNU1(1),MU,MTYPE,IPROC,13,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZNU1(1),MU,MTYPE,IPROC,14,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NXU1(1),MU,MTYPE,IPROC,15,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NYU1(1),MU,MTYPE,IPROC,16,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(NZU1(1),MU,MTYPE,IPROC,17,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(ZU1(1),NZ,MTYPE,IPROC,18,MPI_COMM_EDDY,STATUS,IERR)

              DO IM=1,MU
                WRITE(20, '(A)') 'ZONE I=2 DATAPACKING=POINT'
                write(20,'(3(1x,F10.6))') sqrt(xnu1(im)**2.+ynu1(im)**2.),ynu1(im),znu1(im)
                itr = int(nyu1(im))
                dcell = nxu1(im)
                write(20, '(3(1X,F10.6))') xnu1(im)+uext1*dcell*unvect(1,itr)
     &               ,ynu1(im)+uext1*dcell*unvect(2,itr)
     &               ,znu1(im)+uext1*dcell*unvect(3,itr)
              ENDDO

            ENDIF
          ENDDO
          CLOSE(20)

          ELSE
            CALL MPI_SEND(MIMU,1,MPI_INTEGER,0,1,MPI_COMM_EDDY,IERR)
            IF(MIMU/=0) THEN
              IU1(1:MIMU) = IU(LIMU+1:LIMU+MIMU)
              JU1(1:MIMU) = JU(LIMU+1:LIMU+MIMU)
              KU1(1:MIMU) = KU(LIMU+1:LIMU+MIMU)! + (NZ-2)*MYRANK
              MRKU1(1:MIMU) = MRKU(LIMU+1:LIMU+MIMU)
              XNU1(1:MIMU) = XNU(LIMU+1:LIMU+MIMU)
              YNU1(1:MIMU) = YNU(LIMU+1:LIMU+MIMU)
              ZNU1(1:MIMU) = ZNU(LIMU+1:LIMU+MIMU)
              NXU1(1:MIMU) = NXU(LIMU+1:LIMU+MIMU)
              NYU1(1:MIMU) = NYU(LIMU+1:LIMU+MIMU)
              NZU1(1:MIMU) = NZU(LIMU+1:LIMU+MIMU)

              CALL MPI_SEND(IU1(1),MIMU,MPI_INTEGER,0,2,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(JU1(1),MIMU,MPI_INTEGER,0,3,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(KU1(1),MIMU,MPI_INTEGER,0,4,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(MRKU1(1),MIMU,MPI_INTEGER,0,10,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(XNU1(1),MIMU,MTYPE,0,12,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(YNU1(1),MIMU,MTYPE,0,13,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZNU1(1),MIMU,MTYPE,0,14,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NXU1(1),MIMU,MTYPE,0,15,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NYU1(1),MIMU,MTYPE,0,16,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(NZU1(1),MIMU,MTYPE,0,17,MPI_COMM_EDDY,IERR)
              CALL MPI_SEND(ZU(1),NZ,MTYPE,0,18,MPI_COMM_EDDY,IERR)
            ENDIF

          ENDIF

          ICOUNT = MIMU
          CALL MPI_ALLREDUCE(ICOUNT,MIMUG,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
c          WRITE(6,'(3(1X,A,I6))')
c     &              'MYRANK=',MYRANK,'MIMU=',MIMU,'/',MIMUG
        ENDIF

        RETURN
 120    FORMAT(17(1X,I4),23(1X,F14.8))
c 120  FORMAT(9(1X,I4),11(1X,F12.5))
 122    FORMAT(13(1X,I8),17(1X,F14.8))

        END
C------------------------------------------------------------------------



C------------------------------------------------------------------------
      SUBROUTINE IODIR_SCALAR(NAME,P,DP,NX,NY,NZ,DIR,TIME)
c
      INCLUDE 'common.h'
      INCLUDE 'io.h'
      INCLUDE 'mpif.h'
c
      CHARACTER NAME*(*)
      INTEGER DIR,NX,NY,NZ
      REAL P(NX,NY,NZ),DP(NX,NY,NZ)
c      REAL*4 P1(NX,NY,NZ)
      REAL TIME
C
      INTEGER I,J,K,JP,STATUS(MPI_STATUS_SIZE)

C     IF DIR = -1 READ ELSE WRITE

      IF(DIR==-1) THEN
        OPEN(UNIT=19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &        access='direct',recl=nx*ny*recunit*2)
        DO K=1,NZ
          READ(19,rec=k+myrank*(nz-2)) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        CLOSE(19)

      ELSEIF(DIR==1) THEN
c     note that all files read the input, but only root writes to output
        IF(MYRANK==0) THEN
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED'
     &        ,access='direct',recl=nx*ny*recunit*2)
          IF(MYSIZE==1) THEN
            DO K=1,NZ
              WRITE(19,rec=k) ((P(I,J,K),I=1,NX),J=1,NY)
            ENDDO
          ELSE            
C
            DO K=1,KZ2
              WRITE(19,rec=k) ((P(I,J,K),I=1,NX),J=1,NY)
            ENDDO

            DO JP=1,MYSIZE-1
              CALL MPI_RECV(DP(1,1,1),NX*NY*NZ,MTYPE,JP,0,
     &             MPI_COMM_EDDY,STATUS,IERR)
              DO K=KZ1,KZ2+(JP+1)/MYSIZE
                WRITE(19,rec=k+jp*(nz-2)) ((DP(I,J,K),I=1,NX),J=1,NY)
              ENDDO
            ENDDO
C            
          ENDIF
          CLOSE(19)

        ELSE
          CALL MPI_SEND(P(1,1,1),NX*NY*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        ENDIF

      ELSEIF(DIR==2) THEN
c     note that all files write to output
        OPEN(19,FILE=NAME//PROC(MYRANK),
     $       STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(19) NX,NY,NZ,1
        DO K=1,NZ
          WRITE(19) ((P(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        CLOSE(19)
      ENDIF
c
      return
      end
C------------------------------------------------------------------------



C---- subroutine iompi_3dscalar -----------------------------------------
C      
      subroutine iompi_3dscalar(filename,U,nx,ny,nz,io)
C
C     PURPOSE: Read from or write to a file a 3D real array using
C     MPI subroutines.
C
C------------------------------------------------------------------------
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,io
      real    U(nx,ny,nz)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,k1,k2
      integer fh,filemode
      integer newtype,newtype2
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      if(io==1) then

        if(mysize==1) then
          call mpi_type_contiguous(nx*ny*nz,mtype,newtype,ierr)
        else
          if(myrank==0 .OR. myrank==mysize-1) then
c the number of points for the first and the last processors is equal to
c nx*ny*(nz-1): also the inflow and the outflow layers are considered
c respectively
            call mpi_type_contiguous(nx*ny*(nz-1),mtype,newtype,ierr)
          else
c the number of points for each intermediate processor is equal to
c nx*ny*(nz-2)
            call mpi_type_contiguous(nx*ny*(nz-2),mtype,newtype,ierr)
          endif
        endif

c a new type of variable is generated
        call mpi_type_commit(newtype,ierr)

        if(mysize==1) then
          k1=1
        else
          k1=kz1
          if(myrank==0) k1=1
        endif

        filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

        call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

        disp = 0
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)
        offset = nx*ny*(k1+myrank*(nz-2)-1)
c        write(6,*) 'myrank=',myrank,'k1=',k1,', k2=',k2,k1+myrank*(nz-2)
c
c the variable U is written on a file using a different OFFSET for each
c processor
        call mpi_file_write_at_all(fh,offset,U(:,:,k1),1,newtype,status,ierr)

        call mpi_file_close(fh,ierr)

      else

        if(mysize==1) then
          offset = 0
        else
          offset = nx*ny*(myrank*(nz-2))
          if(myrank==0) offset = 0
        endif

        filemode = MPI_MODE_RDONLY

        call mpi_type_contiguous(nx*ny*nz,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

        disp = 0
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)
        call mpi_file_read_at_all(fh,offset,U(:,:,1),1,newtype,status,ierr)

        call mpi_file_close(fh,ierr)

      endif

c the variable NEWTYPE is deallocated
      call mpi_type_free(newtype,ierr)

      return

      end
*------------------------------------------------------------------------



C---- subroutine iofield2d_plt-------------------N. Beratlis-13 Jan 10---
C
C     PURPOSE: Write 2D field to tecplot formatted file
C
C------------------------------------------------------------------------
      SUBROUTINE IOFIELD2D_PLT(U,V,W,P,XU,XC,YV,YC,ZW,ZC,ZWG,ZCG,NX,NY,NZ,NZG,TIME)

      include 'common.h'
c
c.... Input/Output variables      
      integer nx,ny,nz,nzg
      real    time
      real    xu(nx),xc(nx),yv(ny),yc(ny),zw(nz),zc(nz),zwg(nzg),zcg(nzg)
      real    u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),p(nx,ny,nz)
c
c.... Local variables
      integer   i,j,k,i1,i2,iskip,j1,j2,jskip,k1g,k2g,k1,k2,kskip,nxl,nyl,nzl
      CHARACTER dirwrite*120,filename*120
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: UC,WC,OY
      integer, save :: isave=0

      isave = isave+1
      nwrite2d = nwrite2d+1

c the input file defines the limit of the area and the steps
      OPEN(UNIT=9,FILE='io2d.input',FORM='FORMATTED')
      READ(9,'(A)') DIRWRITE
      READ(9,*) I1,I2,ISKIP
      READ(9,*) J1,J2,JSKIP
      READ(9,*) K1G,K2G,KSKIP
      CLOSE(9)

      NXL = I2-I1+1
      NYL = J2-J1+1

c the local Z limits are established
      K1 = K1G-MYRANK*(NZ-2)
      K2 = K2G-MYRANK*(NZ-2)

c the local Z dimension is found
      IF(K1>1 .AND. K2<NZ) THEN
        NZL = K2-K1+1
      ELSEIF(K1<1 .AND. K2>1 .AND. K2<NZ) THEN
        K1 = 2
        NZL = K2-K1+1
      ELSEIF(K1<1 .AND. K2>NZ-1) THEN
        K1 = 2
        K2 = NZ-1
        NZL = K2-K1+1
      ELSEIF(K1>1 .AND. K1<NZ .AND. K2>NZ-1) THEN
        K2 = NZ-1
        NZL = K2-K1+1
      ELSE
        K1 = NZ
        K2 = 0
        NZL = 0
      ENDIF

      IF(MYRANK.EQ.0) THEN
        IF(ISAVE.EQ.1) THEN
          OPEN(UNIT=9,FILE='time.2d',FORM='FORMATTED')
          WRITE(9,*) 'VARIABLES = "INDEX", "TIME"'
          WRITE(9,*) NWRITE2D,TIME
          CLOSE(9)
        ELSE
          OPEN(UNIT=9,FILE='time.2d',FORM='FORMATTED',POSITION='APPEND')
          WRITE(9,*) NWRITE2D,TIME
          CLOSE(9)
        ENDIF
          
      ENDIF

c
c.... Write gridfile
      IF(MYRANK.EQ.0 .AND. ISAVE.EQ.1) THEN
        DO J=J1,J2,JSKIP
          filename = 'grid3dc.i'//INDEX(I1)//'-'//INDEX(I2)//'.j'//INDEX(J)
     &          //'.k'//INDEX(K1G)//'-'//INDEX(K2G)
          CALL GRID3DASCI(XC(I1:I2),YC(J:J),ZCG(K1G:K2G)
     &        ,TRIM(DIRWRITE)//trim(filename),NXL,1,K2G-K1G+1,ISKIP,KSKIP,icyl)
        ENDDO
        
      ENDIF

      ALLOCATE(UC(NX,1,NZ),WC(NX,1,NZ),OY(NX,1,NZ))

c The variables are evaluated at the centers of the grid cells
      DO J=J1,J2,JSKIP

        UC(2:NX,1,:) = 0.5*(U(2:NX,J,:)+U(1:NX-1,J,:))
        UC(1,1,:) = UC(2,1,:)

        WC(:,1,2:NZ) = 0.5*(W(:,J,2:NZ)+W(:,J,1:NZ-1))
        WC(:,1,1) = WC(:,1,2)

c Y Vorticity Component
        DO K=2,NZ-1
        DO I=2,NX-1
          OY(I,1,K) = (WC(I+1,1, K )-WC(I-1,1, K ))/(XC(I+1)-XC(I-1))
     &             - (UC( I ,1,K+1)-UC( I ,1,K-1))/(ZC(K+1)-ZC(K-1))
        ENDDO
        ENDDO

        OY(1 ,1,: ) = OY( 2  ,1, :  )
        OY(NX,1,: ) = OY(NX-1,1, :  )
        OY(: ,1,1 ) = OY( :  ,1, 2  )
        OY(: ,1,NZ) = OY( :  ,1,NZ-1)

c        filename = 'uc.'//index(nwrite2d)
c     &        //'.i'//INDEX(I1)//'-'//INDEX(I2)//'.j'//INDEX(J)
c     &        //'.k'//INDEX(K1G)//'-'//INDEX(K2G)
c        IF(MYRANK.EQ.0) WRITE(6,'(A)') 'Writing '//TRIM(DIRWRITE)//TRIM(FILENAME)
c        CALL writeplot3dasci(TRIM(DIRWRITE)//TRIM(filename)
c     &      ,UC(I1:I2,J,1:NZ),NXL,1,NZ,K2G-K1G+1,K1,K2,K1G,ISKIP,KSKIP,TIME)

c        filename = 'wc.'//index(nwrite2d)
c     &        //'.i'//INDEX(I1)//'-'//INDEX(I2)//'.j'//INDEX(J)
c     &        //'.k'//INDEX(K1G)//'-'//INDEX(K2G)
c        IF(MYRANK.EQ.0) WRITE(6,'(A)') 'Writing '//TRIM(DIRWRITE)//TRIM(FILENAME)
c        CALL writeplot3dasci(TRIM(DIRWRITE)//TRIM(FILENAME)
c     &      ,WC(I1:I2,J,1:NZ),NXL,1,NZ,K2G-K1G+1,K1,K2,K1G,ISKIP,KSKIP,TIME)

c The Y vorticity component is written on a file
        filename = 'oy.'//index(nwrite2d)
     &        //'.i'//INDEX(I1)//'-'//INDEX(I2)//'.j'//INDEX(J)
     &        //'.k'//INDEX(K1G)//'-'//INDEX(K2G)
c        IF(MYRANK.EQ.0) WRITE(6,'(A)') 'Writing '//TRIM(DIRWRITE)//TRIM(FILENAME)
        CALL writeplot3dasci(TRIM(DIRWRITE)//trim(filename)
     &      ,OY(I1:I2,1,1:NZ),NXL,1,NZ,K2G-K1G+1,K1,K2,K1G,ISKIP,KSKIP,TIME)

c        filename = 'pres.'//index(nwrite2d)
c     &        //'.i'//INDEX(I1)//'-'//INDEX(I2)//'.j'//INDEX(J)
c     &        //'.k'//INDEX(K1G)//'-'//INDEX(K2G)
c        IF(MYRANK.EQ.0) WRITE(6,'(A)') 'Writing '//TRIM(DIRWRITE)//TRIM(FILENAME)
c        call writeplot3dasci(TRIM(DIRWRITE)//TRIM(FILENAME)
c     &      ,P(I1:I2,J,1:NZ),NXL,1,NZ,K2G-K1G+1,K1,K2,K1G,ISKIP,KSKIP,TIME)

      ENDDO

      DEALLOCATE(UC,WC,OY)

      RETURN

      END
C------------------------------------------------------------------------



c     ***************************************************************
c
C     subroutine for I/O - read or write a scalar quantity
c
c     ***************************************************************
c
      SUBROUTINE writeplot3d(NAME,P,NX,NY,NZ,NZG,K1,K2,K1G,ISKIP,KSKIP,TIME)
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
      CHARACTER NAME*(*)
      INTEGER DIR,NX,NY,NZ,NZG,K1,K2,ISKIP,KSKIP,KG,K1G
      REAL P(NX,NY,NZ),DP(NX,NY,NZ)
c      REAL*4 P1(NX,NY,NZ)
      REAL TIME
C
      INTEGER I,J,K,M,N,NZL,KL1,KL2,JP,STATUS(MPI_STATUS_SIZE)

      N=0
      DO I=1,NX,ISKIP
        N=N+1
      ENDDO
      M=0
      DO K=1,NZG,KSKIP
        M=M+1
      ENDDO


c     note that all files read the input, but only root writes to output
        IF(MYRANK==0) THEN
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
          WRITE(19) N,NY,M,1
c     WRITE(6,*) 'Writing to file header:',N,NY,M,1

          IF(MYSIZE==1) THEN
            DO K=K1,K2,KSKIP
              WRITE(19) ((REAL(P(I,J,K),4),I=1,NX,ISKIP),J=1,NY)
            ENDDO

          ELSE
C
c            DO K=K1,K2
            M=0
            K=K1
            DO WHILE(K<=K2)
              IF(MOD(K-K1,KSKIP)==0) THEN
                WRITE(19) ((REAL(P(I,J,K),4),I=1,NX,ISKIP),J=1,NY)
                N=0
                M=M+1
c                write(6,*) k,k-k1,mod(k-k1,kskip),m
              ELSE
                N=N+1
              ENDIF
              K=K+1
            ENDDO

            DO JP=1,MYSIZE-1
              CALL MPI_RECV(NZL,1,MPI_INTEGER,JP,0,MPI_COMM_EDDY,STATUS,IERR)
              IF(NZL>0) THEN
                CALL MPI_RECV(KL1,1,MPI_INTEGER,JP,1,MPI_COMM_EDDY,STATUS,IERR)
                CALL MPI_RECV(KL2,1,MPI_INTEGER,JP,2,MPI_COMM_EDDY,STATUS,IERR)
                CALL MPI_RECV(DP(1,1,1),NX*NY*NZ,MTYPE,JP,3,MPI_COMM_EDDY,STATUS,IERR)
c                WRITE(6,*) 'Receiving from',JP,KL1,KL2 
c                K=KL1+(KSKIP-N-1)
c                write(6,*) kl1,kl2,m,n,k
c                DO K=KL1,KL2
                K=KL1
                DO WHILE(K<=KL2)
                  KG = K+JP*(NZ-2)
                  IF(MOD(KG-K1G,KSKIP)==0) THEN
                    WRITE(19) ((REAL(DP(I,J,K),4),I=1,NX,ISKIP),J=1,NY)
                    N=0
                    M=M+1
c                    write(6,*) kg,kg-k1g,mod(kg-k1g,kskip),m
                  ELSE
                    N=N+1
                  ENDIF
                  K=K+1
                ENDDO
              ENDIF
            ENDDO
C
c            WRITE(6,*) 'M=',M
          ENDIF
          CLOSE(19)

        ELSE
          NZL = K2-K1+1
c          write(6,*) 'myrank=',myrank,'sending',k1,k2,nzl
          CALL MPI_SEND(NZL,1,MPI_INTEGER,0,0,MPI_COMM_EDDY,IERR)
          IF(NZL>0) THEN
             CALL MPI_SEND(K1,1,MPI_INTEGER,0,1,MPI_COMM_EDDY,IERR)
             CALL MPI_SEND(K2,1,MPI_INTEGER,0,2,MPI_COMM_EDDY,IERR)
             CALL MPI_SEND(P(1,1,1),NX*NY*NZ,MTYPE,0,3,MPI_COMM_EDDY,IERR)
          ENDIF
        ENDIF
c
      return
      end
*------------------------------------------------------------------------

c     ***************************************************************
c
C     subroutine for I/O - read or write a scalar quantity
c
c     ***************************************************************
c
      SUBROUTINE writeplot3dasci(NAME,P,NX,NY,NZ,NZG,K1,K2,K1G,ISKIP,KSKIP,TIME)
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
      CHARACTER NAME*(*)
      INTEGER DIR,NX,NY,NZ,NZG,K1,K2,ISKIP,KSKIP,KG,K1G
      REAL P(NX,NY,NZ),DP(NX,NY,NZ)
c      REAL*4 P1(NX,NY,NZ)
      REAL TIME
C
      INTEGER I,J,K,M,N,NZL,KL1,KL2,JP,STATUS(MPI_STATUS_SIZE)

      N=0
      DO I=1,NX,ISKIP
        N=N+1
      ENDDO
      M=0
      DO K=1,NZG,KSKIP
        M=M+1
      ENDDO


c     note that all files read the input, but only root writes to output
        IF(MYRANK==0) THEN
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='FORMATTED')
          WRITE(19,*) N,NY,M,1

          IF(MYSIZE==1) THEN
            DO K=K1,K2,KSKIP
              WRITE(19,*) ((P(I,J,K),I=1,NX,ISKIP),J=1,NY)
            ENDDO

          ELSE
C
c The processor of rank 0 writes on a file the local variables
            M=0
            K=K1
            DO WHILE(K<=K2)
              IF(MOD(K-K1,KSKIP)==0) THEN
                WRITE(19,*) ((P(I,J,K),I=1,NX,ISKIP),J=1,NY)
                N=0
                M=M+1
              ELSE
                N=N+1
              ENDIF
              K=K+1
            ENDDO

            DO JP=1,MYSIZE-1
c The processor of rank 0 receives NZL from the other processors
              CALL MPI_RECV(NZL,1,MPI_INTEGER,JP,0,MPI_COMM_EDDY,STATUS,IERR)
              IF(NZL>0) THEN
c The processor of rank 0 receives the Z limits and the variable from the
c other processors
                CALL MPI_RECV(KL1,1,MPI_INTEGER,JP,1,MPI_COMM_EDDY,STATUS,IERR)
                CALL MPI_RECV(KL2,1,MPI_INTEGER,JP,2,MPI_COMM_EDDY,STATUS,IERR)
                CALL MPI_RECV(DP(1,1,1),NX*NY*NZ,MTYPE,JP,3,MPI_COMM_EDDY,STATUS,IERR)
                K=KL1
                DO WHILE(K<=KL2)
                  KG = K+JP*(NZ-2)
                  IF(MOD(KG-K1G,KSKIP)==0) THEN
                    WRITE(19,*) ((DP(I,J,K),I=1,NX,ISKIP),J=1,NY)
                    N=0
                    M=M+1
                  ELSE
                    N=N+1
                  ENDIF
                  K=K+1
                ENDDO
              ENDIF
            ENDDO
C
          ENDIF
          CLOSE(19)

        ELSE
c The other processors send NZL,K1,K2,P to the processor 0
          NZL = K2-K1+1
          CALL MPI_SEND(NZL,1,MPI_INTEGER,0,0,MPI_COMM_EDDY,IERR)
          IF(NZL>0) THEN
             CALL MPI_SEND(K1,1,MPI_INTEGER,0,1,MPI_COMM_EDDY,IERR)
             CALL MPI_SEND(K2,1,MPI_INTEGER,0,2,MPI_COMM_EDDY,IERR)
             CALL MPI_SEND(P(1,1,1),NX*NY*NZ,MTYPE,0,3,MPI_COMM_EDDY,IERR)
          ENDIF
        ENDIF
c
      return
      end
*------------------------------------------------------------------------


C-----------------------------------------------------------------------
      SUBROUTINE GRID3D(X,Y,Z,NAME,NX,NY,NZ,ISKIP,KSKIP,ICYL)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    Write 3D grid file.
C
C-----------------------------------------------------------------------
C
      implicit none
      INCLUDE 'mpif.h'
C
C-----------------------------------------------------------------------
C
c...Input variables and arrays.
      INTEGER NX,NY,NZ,ICYL,ISKIP,KSKIP
      REAL    X(NX),Y(NY),Z(NZ)
      CHARACTER NAME*(*)
c
c...Local variables and arrays.
      INTEGER I,J,K,JP,M,N
C
C
      N=0
      DO I=1,NX,ISKIP
        N=N+1
      ENDDO
      M=0
      DO K=1,NZ,KSKIP
        M=M+1
      ENDDO

c      IF (MYRANK==0) THEN
        OPEN(166, FILE=NAME, FORM='UNFORMATTED')
        WRITE(166) N,NY,M
c        WRITE(6,*) 'Writing file',NAME,N,NY,M
        IF(ICYL==0) THEN
          WRITE(166) (((REAL(X(I),4),I=1,NX,ISKIP),J=1,NY),K=1,NZ,KSKIP),
     &               (((REAL(Y(J),4),I=1,NX,ISKIP),J=1,NY),K=1,NZ,KSKIP),
     &               (((REAL(Z(K),4),I=1,NX,ISKIP),J=1,NY),K=1,NZ,KSKIP)
        ELSE
          WRITE(166) (((REAL(X(I)*COS(Y(J)),4),I=1,NX,ISKIP),J=1,NY),K=1,NZ,KSKIP),
     &               (((REAL(X(I)*SIN(Y(J)),4),I=1,NX,ISKIP),J=1,NY),K=1,NZ,KSKIP),
     &               (((REAL(Z(K),4),I=1,NX,ISKIP),J=1,NY),K=1,NZ,KSKIP)
        ENDIF
        CLOSE(166)
c       ENDIF

       RETURN

       END

C-----------------------------------------------------------------------
 


C-----------------------------------------------------------------------
      SUBROUTINE GRID3DASCI(X,Y,Z,NAME,NX,NY,NZ,ISKIP,KSKIP,ICYL)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    Write 3D grid file.
C
C-----------------------------------------------------------------------
C
      implicit none
      INCLUDE 'mpif.h'
C
C-----------------------------------------------------------------------
C
c...Input variables and arrays.
      INTEGER NX,NY,NZ,ICYL,ISKIP,KSKIP
      REAL    X(NX),Y(NY),Z(NZ)
      CHARACTER NAME*(*)
c
c...Local variables and arrays.
      INTEGER I,J,K,JP,M,N
C
C
      N=0
      DO I=1,NX,ISKIP
        N=N+1
      ENDDO
      M=0
      DO K=1,NZ,KSKIP
        M=M+1
      ENDDO

        OPEN(166, FILE=NAME, FORM='FORMATTED')
        WRITE(166,*) N,NY,M
        IF(ICYL==0) THEN
          WRITE(166,*) (((REAL(X(I),4),I=1,NX,ISKIP),J=1,NY),K=1,NZ,KSKIP),
     &               (((REAL(Y(J),4),I=1,NX,ISKIP),J=1,NY),K=1,NZ,KSKIP),
     &               (((REAL(Z(K),4),I=1,NX,ISKIP),J=1,NY),K=1,NZ,KSKIP)
        ELSE
          WRITE(166,*) (((REAL(X(I)*COS(Y(J)),4),I=1,NX,ISKIP),J=1,NY),K=1,NZ,KSKIP),
     &               (((REAL(X(I)*SIN(Y(J)),4),I=1,NX,ISKIP),J=1,NY),K=1,NZ,KSKIP),
     &               (((REAL(Z(K),4),I=1,NX,ISKIP),J=1,NY),K=1,NZ,KSKIP)
        ENDIF
        CLOSE(166)
c       ENDIF

       RETURN

       END

C-----------------------------------------------------------------------
 



C---- call ioscaltrinobd ------------------------N. Beratlis-31 March 2010---
C
C     PURPOSE: Write components of tau located on center of triangel
C
C------------------------------------------------------------------------
      subroutine ioscaltrinobd(filename,trino,pbd,mrkp,nfacet,nfacetmax,nfacetot,ilb,ile)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c...  Input/Output arrays
      integer nbd,nfacet,nfacetmax,nfacetot,ilb,ile
      integer mrkp(nfacet),trino(nfacet)
c      real    vertexc(3,nfacet)
      real    pbd(nfacet)
      character*(*) filename
c
c.... Local arrays
      integer ibd,n,i,j,jp,ng
      integer STATUS(MPI_STATUS_SIZE)
      real    tmp1(nfacetmax)
      integer trino1(nfacetmax)
c      real    vertexc1(3,nfacetmax)
      real    psi
c
c...... Functions
      real    anglerad

c Local number of valid probes
      n = count(mrkp(ilb:ile)==1)    !!!!!! instead of mrkp==1
c Global number of valid probes
      CALL MPI_ALLREDUCE(n,ng,1,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
c      CALL MPI_ALLREDUCE(nfacet,nfacetmax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_EDDY,IERR)

      IF(MYRANK==0) THEN
        open(unit=20,file=trim(filename),form='formatted')
c        write(20,'(A,A)') 'VARIABLES = "TRINO", "VAR"'
c        do ibd=1,nbd
c          write(20,'(A,I8,A)') 'ZONE I=',ng,', DATAPACKING=POINT'
          write(20,*) ng
c          do i=lb(ibd)+1,lb(ibd)+mb(ibd)
c The processor 0 writes on a file the local values of the triangle
c number and of the relevant stress
          do i=ilb,ile
            if(mrkp(i)==1) then
              write(20,'((I7,1X,E18.11))') trino(i),pbd(i)
            endif
          enddo

          IF(MYSIZE>1) THEN
            
            DO JP=1,MYSIZE-1
c Processor 0 receives the value of N from the processor JP
              CALL MPI_RECV(n,1,MPI_INTEGER,JP,JP,MPI_COMM_EDDY,STATUS,IERR)
              IF(n>0) THEN
c Processor 0 receives the values of TRINO and PBD from the processor JP
                CALL MPI_RECV(trino1,n,MPI_INTEGER,JP,1,MPI_COMM_EDDY,STATUS,IERR)
                CALL MPI_RECV(tmp1,n,MTYPE,JP,2,MPI_COMM_EDDY,STATUS,IERR)
                do i=1,n
                  write(20,'((I7,1X,E18.11))') trino1(i),tmp1(i)
                enddo             
              ENDIF
            ENDDO
          ENDIF
c        ENDDO
        close(20)

      ELSE

c        do ibd=1,nbd
          j=0
c          do i=lb(ibd)+1,lb(ibd)+mb(ibd)
c The values of TRINO and PBD are copied on temporary arrays
          do i=ilb,ile
            if(mrkp(i)==1) then
              j=j+1
c              vertexc1(:,j)=vertexc(:,i)
              trino1(j) = trino(i)
              tmp1(j) = pbd(i)
            endif
          enddo
          n=j
c The value of N is sent to the processor 0
          CALL MPI_SEND(n,1,MPI_INTEGER,0,MYRANK,MPI_COMM_EDDY,IERR)
          if(n>0) then
c The values of TRINO and PBD are sent to the processor 0
            CALL MPI_SEND(trino1,n,MPI_INTEGER,0,1,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(tmp1,n,MTYPE,0,2,MPI_COMM_EDDY,IERR)
          endif
c        enddo
      ENDIF
        
      return

      end


      SUBROUTINE IOSCALAR_MEAN(NAME,TEMP,NX,NY,NZ,NKRR,DIR,TIME,DTM1,nstep)

      INCLUDE 'common.h'
      INCLUDE 'mpif.h'

      CHARACTER NAME*(*)
      INTEGER DIR,NX,NY,NZ,NKRR,nstep
      REAL TEMP(NX,NY,NKRR)
      REAL TIME,DTM1

      INTEGER I,J,K,JP,STATUS(MPI_STATUS_SIZE)
      INTEGER KX,KXX1,KXX2,KXX3,MK,NKR,COUNTER,KMIN,KMAX,NK

      IF(DIR==-1) THEN
        OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
        READ(19) I,J,K,JP
        DO k=1,NKRR
        READ(19) ((TEMP(I,J,K),I=1,NX),J=1,NY)
        ENDDO
        READ(19) nstep
        READ(19) TIME
        READ(19) DTM1,grav
        WRITE(*,*)"READING DONE"
        CLOSE(19)

      ELSEIF(DIR==1) THEN
          OPEN(19,FILE=NAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
          WRITE(19) NX,NY,MYSIZE*(NZ-2)+2,1
          DO k=1,NKRR
          WRITE(19) ((TEMP(I,J,K),I=1,NX),J=1,NY)
          ENDDO
          WRITE(19) nstep
          WRITE(19) TIME
          WRITE(19) DTM1,grav
          CLOSE(19)
      ENDIF

      return
      end






C------------------------------------------------------------------------


