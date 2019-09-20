C


C-----SUBROUTINE-INPALL-------------------------P. FLOHR--30/10/1993----
C
      SUBROUTINE INPALL(nx,ny,nz,nzg)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - READ INPUT-DATA FROM UNIT 1
C                                                                       
C     INPUT:      - NAMELIST DIRECTED DATA FROM UNIT 1
C                                                                       
C-----------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'io.h'
      INCLUDE 'mpif.h'
C                                                                       
C-----------------------------------------------------------------------
C
      CHARACTER   TEXT(72)*1,ROUTINE*20,ERRMSG*70
      INTEGER     i,j,nx,ny,nz,nzg     
C
C-----------------------------------------------------------------------
C                                                    NAMELIST STATEMENTS
C-----------------------------------------------------------------------
C
      NAMELIST /STRPRM/ STR1,STR2,STR3,STR4,STR5,STR6
      NAMELIST /INTGRD/ IGRID,JGRID,KGRID,IOGRID
      NAMELIST /CMPDMN/ XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
      NAMELIST /BNDPRM/ ITYPE
      NAMELIST /SOLPRM/ ISCHM,IBM,ISGS,IPRES,IFIELD,IDENS,INOISE
!      NAMELIST /ITRPRM/ ITMAX,ITINI,ITSCR,ISTAT,ITSER,ITRES,ITPOST,IT5P,ITPLN,IVRTX,ITBDY,
!    &                  ITCALF,IT2D,
!    &                  ITVPFIELD,NWRITE,NWRITE2D,NWRITEBDY,IOLVL,IOCLOCK
!      NAMELIST /SOLPRM/ ISCHM,IBM,ISGS,IPRES,IFIELD,IDENS,INOISE
      NAMELIST /ITRPRM/ ITMAX,ITINI,ITSCR,ISTAT,ITSER,ITRES,ITPOST,IT5P,
     &   KMIN5P,KMAX5P,ITPLN,IVRTX,ITBDY,
     &   ITCALF,IT2D,ITVPFIELD,NWRITE,NWRITE2D,NWRITEBDY,IOLVL,IOCLOCK
      NAMELIST /CONTRL/ ICYL,IMPLX,IMPLY,IMPLZ,IMPLCX,IMPLCY,IMPLCZ
     &     ,IMPLXDCMP,IMPLYDCMP,IMPLZDCMP,ICALF,ICFL
      NAMELIST /TIMEST/ CFLC,TSTEP,DTSAVE,TINI,RESSTEP
      NAMELIST /FLOPRM/ RU1,CSS,EPS,DPDX,DPDY,DPDZ,FRN,PRN,DENP1,RHO_0,
     & gt1, gt2, g_orig
      NAMELIST /MOVWAL/ UMOV1,UMOV2,VMOV1,VMOV2,WMOV1,WMOV2
      NAMELIST /DUMPRM/ FRE,AMP,CPHS,KWAV,DUMMYFP1,DUMMYFP2,DUMMYFP3,
     &                   DUMMYFP4,DUMMYFP5
      NAMELIST /GRIDDM/ NX,NY,NZ,NXPROCS,NYPROCS,NZPROCS
      NAMELIST /IMBPRM/ IPRSDV,IVELRCN,UEXT1,UEXT2,UEXT3,NFLU,NFLP,IDOMY
     &     ,IMB3D,IMBDCMP,IMBOVLP,NFCMAX,NINTMAX,NNTR2D,NNTR3D,ITAGY,
     &      IGRDINT,IMLS,IZERO
      NAMELIST /IOPRMS/ RECUNIT, BACKGROUND_TIME
      NAMELIST /STATS/  ITMPRB,ITKE,YSPCTR,ITMAVG,IREGTMAVG
      NAMELIST /HYBPRM/  start_feed,end_feed,stride_feed,index_feed,path_feed


      ROUTINE='inpall'
c
c-----------------------------------------------------------------------
c                                              define MPI REAL precision
c-----------------------------------------------------------------------
      IF(MYRANK==0) THEN
        IF(KIND(RU1)==4) THEN
          MTYPE  = MPI_REAL
          MTYPE2 = MPI_2REAL
        ELSEIF(KIND(RU1)==8) THEN
          MTYPE  = MPI_DOUBLE_PRECISION
          MTYPE2 = MPI_2DOUBLE_PRECISION
        ENDIF
      ENDIF
C
      CALL MPI_BCAST(MTYPE ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(MTYPE2,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
C
C-----------------------------------------------------------------------
C                                                            READ UNIT 1
C-----------------------------------------------------------------------
C
      IF(MYRANK==0) THEN
         OPEN (1,STATUS='OLD',FILE='eddy.input')
         OPEN (2,STATUS='SCRATCH')
         OPEN (3,STATUS='UNKNOWN'
     &        ,FILE='eddy.input.backup')
C 
C-----------------------------------------------------------------------
C                            COPY NAMELIST-DIRECTED INPUT DATA TO UNIT 2
C-----------------------------------------------------------------------
C
         REWIND 1 
         REWIND 2 
         REWIND 3
 10      READ (1,'(72(A1))',END=20) (TEXT(I),I=1,72)
         WRITE(3,'(72(A1))') (TEXT(I),I=1,72)
         IF (TEXT(1) .EQ. '*') GOTO 10
         WRITE(2,'(72(A1))') (TEXT(I),I=1,72)
         GOTO 10
 20      REWIND 2
      ENDIF

C 
C-----------------------------------------------------------------------
C                                    SET DEFAULT-VALUES & READ NAMELISTS
C-----------------------------------------------------------------------
C
C....."STRPRM" 
C
c
c names of the grid and geometry files
c 
      IF(MYRANK==0) THEN
         READ (2,STRPRM) 
         WRITE(6,STRPRM)
      ENDIF

      CALL MPI_BCAST(STR1,60,MPI_CHARACTER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(STR2,60,MPI_CHARACTER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(STR3,60,MPI_CHARACTER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(STR4,60,MPI_CHARACTER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(STR5,60,MPI_CHARACTER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(STR6,60,MPI_CHARACTER,0,MPI_COMM_EDDY,IERR)

C
C....."INTGRD"
C
c
c grid type in each direction
c     
      IF(MYRANK==0) THEN
         READ (2,INTGRD)
         WRITE(6,INTGRD)
      ENDIF
C
      CALL MPI_BCAST(IGRID,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(JGRID,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(KGRID,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IOGRID,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
C
C....."CMPDMN"
C
c
c domain size in each direction
c
      IF(MYRANK==0) THEN
        READ (2,CMPDMN) 
        WRITE(6,CMPDMN)
      ENDIF

      CALL MPI_BCAST(XMIN,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(XMAX,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(YMIN,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(YMAX,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ZMIN,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ZMAX,1,MTYPE,0,MPI_COMM_EDDY,IERR)
C
C....."BNDPRM" 
C
c
c boundary conditions
c 
      IF (MYRANK==0) THEN
         READ (2,BNDPRM)
         WRITE(6,BNDPRM)
C
         MP = 1
         LP = 1
         NP = 1
C
         IF(ITYPE(2)==500) MP = 0
         IF(ITYPE(4)==500) LP = 0
         IF(ITYPE(6)==500) NP = 0
C
       ENDIF
C
      CALL MPI_BCAST(ITYPE,6,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(MP   ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(LP   ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NP   ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
C
      IF(ITYPE(2)==500.AND.ITYPE(1)/=500) THEN
        ERRMSG=' ITYPE(1)/=ITYPE(2)== '
        CALL ENDJOB(ROUTINE,ERRMSG,500)
      ENDIF
      IF(ITYPE(4)==500.AND.ITYPE(3)/=500) THEN
        ERRMSG=' ITYPE(3)/=ITYPE(4)== '
        CALL ENDJOB(ROUTINE,ERRMSG,500)
      ENDIF
      IF(ITYPE(6)==500.AND.ITYPE(5)/=500) THEN
        ERRMSG=' ITYPE(5)/=ITYPE(6)== '
        CALL ENDJOB(ROUTINE,ERRMSG,500)
      ENDIF
C
      IF(MYSIZE>1) THEN
c
c periodic boundary conditions are set on Z boundaries between
c processors
c         
C     THE UPPER Z B.C. IS SET TO NULL FOR PROC 0
         IF(MYRANK==0) THEN
            ITYPE(6)=0
C     THE LOWER Z B.C. IS SET TO NULL FOR PROC MYSIZE-1
         ELSEIF(MYRANK==MYSIZE-1) THEN
            ITYPE(5)=0
         ELSE
            ITYPE(5)=0
            ITYPE(6)=0
         ENDIF
      ENDIF
C
C....."SOLPRM" 
C
      IF(MYRANK==0) THEN
        READ (2,SOLPRM)
        WRITE(6,SOLPRM)
C
c.....second order Adams-Bashfort Scheme      
        IF(ISCHM==1) THEN
          GAM(1) =  1.5
          GAM(2) =  0.0
          GAM(3) =  0.0
          RHO(1) = -0.5
          RHO(2) =  0.0
          RHO(3) =  0.0
c.....third order Runge-Kutta scheme
        ELSEIF(ISCHM==3) THEN
          GAM(1) =  8./15.
          GAM(2) =  5./12.
          GAM(3) =  3./4.
          RHO(1) =  0.
          RHO(2) = -17./60. 
          RHO(3) = -5./12.             
        ENDIF
        
        ALF(:) = GAM(:)+RHO(:)

      ENDIF

      CALL MPI_BCAST(ISCHM ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IBM   ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ISGS  ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IPRES ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IFIELD,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR) 
      CALL MPI_BCAST(IDENS,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR) 
      CALL MPI_BCAST(INOISE,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR) 

      IF(ISCHM/=1.AND.ISCHM/=3) THEN
        ERRMSG=' UNKNOWN TIME-ADVANCEMENT SCHEME '
        CALL ENDJOB(ROUTINE,ERRMSG,ISCHM)
      ENDIF
C
!!!!!!      IF(ISGS<0.AND.ISGS>4) THEN      
      IF(ISGS<0.OR.ISGS>6) THEN   
        ERRMSG=' UNKNOWN SGS MODEL '
        CALL ENDJOB(ROUTINE,ERRMSG,ISGS)
      ENDIF
C
      IF(IPRES<1.OR.IPRES>6) THEN
        ERRMSG=' UNKNOWN POISSON SOLVER NUMBER '
        CALL ENDJOB(ROUTINE,ERRMSG,IPRES)
      ENDIF
C
      CALL MPI_BCAST(GAM,3,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(RHO,3,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ALF,3,MTYPE,0,MPI_COMM_EDDY,IERR)
C 
C....."ITRPRM" 
C
c
c steps of the output
c
      IF(MYRANK==0) THEN
        READ (2,ITRPRM) 
        WRITE(6,ITRPRM)
      ENDIF
C
      CALL MPI_BCAST(ITMAX ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ITINI ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ITSCR ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ISTAT ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ITSER ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ITRES ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ITPOST ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IT5P ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(KMIN5P ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(KMAX5P ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR) 
      CALL MPI_BCAST(ITPLN ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IVRTX ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ITBDY ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ITCALF,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IT2D  ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ITVPFIELD,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NWRITE,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NWRITE2D,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NWRITEBDY,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IOLVL  ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IOCLOCK,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
C      
      CALL MPI_BARRIER(MPI_COMM_EDDY,IERR)
C
C....."FLOPRM"
C
c
c flow parameters
c
      IF(MYRANK==0) THEN
        READ (2,FLOPRM) 
        WRITE(6,FLOPRM)
      ENDIF

      CALL MPI_BCAST(RU1 ,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(CSS ,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(EPS ,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(DPDX,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(DPDY,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(DPDZ,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(FRN, 1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(PRN, 1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(DENP1, 1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(RHO_0, 1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(gt1, 1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(gt2, 1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(g_orig, 1,MTYPE,0,MPI_COMM_EDDY,IERR)

C
C....."CONTRL"
C
c
c kind of grid and criterion of discretization for viscous and
c convective terms
c
      IF(MYRANK==0) THEN
        READ (2,CONTRL)
        WRITE(6,CONTRL)
      ENDIF

      CALL MPI_BCAST(ICYL ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IMPLX,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IMPLY,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IMPLZ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IMPLCX,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IMPLCY,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IMPLCZ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IMPLXDCMP,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IMPLYDCMP,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IMPLZDCMP,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ICALF,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ICFL ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)

      IF(IMPLXDCMP>0 .AND. IMPLX==0 .AND. IMPLCX==0) THEN
        ERRMSG='PARTIALLY IMPLICIT TREATMENT IN X REQUIRES '
     &        //'IMPLX OR IMPLCX>0'
        CALL ENDJOB(ROUTINE,ERRMSG,MYSIZE)
      ENDIF

      IF(IMPLYDCMP>0 .AND. IMPLY==0 .AND. IMPLCY==0) THEN
        ERRMSG='PARTIALLY IMPLICIT TREATMENT IN Y REQUIRES '
     &        //'IMPLY OR IMPLCY>0'
        CALL ENDJOB(ROUTINE,ERRMSG,MYSIZE)
      ENDIF

      IF(IMPLZDCMP>0 .AND. IMPLZ==0 .AND. IMPLCZ==0) THEN
        ERRMSG='PARTIALLY IMPLICIT TREATMENT IN Z REQUIRES '
     &        //'IMPLZ OR IMPLCZ>0'
        CALL ENDJOB(ROUTINE,ERRMSG,MYSIZE)
      ENDIF

c      IF(IIMPLX>0 .AND. IIMPLY==0) THEN
c        ERRMSG='IIMPLX>0 REQUIRES IIMPLY>0 TOO. '
c        CALL ENDJOB(ROUTINE,ERRMSG,MYSIZE)
c      ENDIF

      IF(ICYL+IGRID+JGRID+KGRID.EQ.0) THEN
         STRETCH = .False.
      ELSE
         STRETCH = .True.
      ENDIF

C
C....."IMBPRM"
C
c
c parameters related to the immersed-boundary
c
      IF(MYRANK==0) THEN
        READ (2,IMBPRM)
        WRITE(6,IMBPRM)
      ENDIF

      CALL MPI_BCAST(IPRSDV,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IVELRCN,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(UEXT1,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(UEXT2,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(UEXT3,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NFLU,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NFLP,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IDOMY,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IMB3D,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IMBDCMP,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IMBOVLP,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NFCMAX,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NINTMAX,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NNTR2D,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NNTR3D,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ITAGY,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IGRDINT,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IMLS,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IZERO,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
C
C....."TIMEST"
C
c
c CFL number and TIME STEP
c
      IF(MYRANK==0) THEN
        READ (2,TIMEST)
        WRITE(6,TIMEST)
      ENDIF

      CALL MPI_BCAST(CFLC  ,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(TSTEP ,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(DTSAVE,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(TINI  ,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(RESSTEP  ,1,MTYPE,0,MPI_COMM_EDDY,IERR)
C
C....."STATS"
C
c
c probes
c
      IF(MYRANK==0) THEN
        READ (2,STATS)
        WRITE(6,STATS)
      ENDIF

      CALL MPI_BCAST(ITMPRB,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ITKE  ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(YSPCTR,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ITMAVG,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IREGTMAVG,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
C
C....."MOVWAL" 
C
c
c moving walls
c
      IF(MYRANK==0) THEN
        READ (2,MOVWAL)
        WRITE(6,MOVWAL)
      ENDIF

      CALL MPI_BCAST(UMOV1,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(UMOV2,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(VMOV1,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(VMOV2,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(WMOV1,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(WMOV2,1,MTYPE,0,MPI_COMM_EDDY,IERR)
C
C....."DUMPRM" 
C
c
c dummy parameters for various cases
c
      IF(MYRANK==0) THEN
        READ (2,DUMPRM)
        WRITE(6,DUMPRM)
      ENDIF

      CALL MPI_BCAST(FRE ,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(AMP ,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(CPHS,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(KWAV,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(DUMMYFP1,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(DUMMYFP2,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(DUMMYFP3,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(DUMMYFP4,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(DUMMYFP5,1,MTYPE,0,MPI_COMM_EDDY,IERR)
C
C....."IOPRMS" 
C
      IF(MYRANK==0) THEN
        READ (2,IOPRMS)
        WRITE(6,IOPRMS)
      ENDIF

      CALL MPI_BCAST(RECUNIT,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(BACKGROUND_TIME,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)

! ....."HYBPRM" 
! hybrid parameters

 
      IF(MYRANK==0) THEN
         READ (2,HYBPRM) 
         WRITE(6,HYBPRM)
      ENDIF

      CALL MPI_BCAST(start_feed,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(end_feed,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(stride_feed,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(index_feed,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(path_feed,600,MPI_CHARACTER,0,MPI_COMM_EDDY,IERR)

C....."GRIDDM"
C 
c
c points of the computational grid
c    
      IF(MYRANK==0) THEN
         READ (2,GRIDDM)
         WRITE(6,GRIDDM)
C
         NZG=NZ
         NZ=(NZG-2)/MYSIZE+2
c
c NZG: total number of grid points along the Z direction
c (including the ghost layers)
c NZ: number of grid points along the Z direction for each processor
c (including the ghost layers)
C
      ENDIF
C
      CALL MPI_BCAST(NX ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NY ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NZ ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NZG,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NXPROCS ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NYPROCS ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(NZPROCS ,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
C
      IF(MOD(NZG-2,MYSIZE)>0) THEN
        ERRMSG='NZG IS NOT DIVISIBLE BY # PROC'
        CALL ENDJOB(ROUTINE,ERRMSG,MYSIZE)
      ENDIF

      IF(IPRES<4 .AND. MOD(NY-2,MYSIZE)>0) THEN
        ERRMSG='NYG IS NOT DIVISIBLE BY # PROC'
c        CALL ENDJOB(ROUTINE,ERRMSG,MYSIZE)
      ENDIF
C
      IF(MYSIZE>NY-2 .AND. IPRES/=5) THEN
        if(myrank.eq.0) write(6,*)
     &        'Switching to parallel poisson solver, IPRES=5'
        IPRES = 5
      ENDIF

      IF(IPRES==5 .AND. MYSIZE<=NY-2) THEN
        IPRES = 4
      ENDIF

      IF(MYRANK==0) THEN
         CLOSE (1)
         CLOSE (2)
         CLOSE (3)
      ENDIF
C
C..... Partially implicit treatment
C
c
c Variables for the local implicit treatment of the U and V
c derivatives along the X and Y directions
c
      RUIMPLX = 1.0
      RVIMPLX = 1.0
      RUIMPLY = 1.0
      RVIMPLY = 1.0
      IMPLXULIM(1,1) = IBU
      IMPLXULIM(1,2) = IEU
      IMPLXVLIM(1,1) = 2
      IMPLXVLIM(1,2) = IX2
      IMPLYULIM(1,1) = 2
      IMPLYULIM(1,2) = IX2
      IMPLYVLIM(1,1) = 2
      IMPLYVLIM(1,2) = IX2
      NIMPLX = 0
      NIMPLY = 0

      IF(IMPLXDCMP>0 .OR. IMPLYDCMP>0) THEN
c
c partial implicit treatment for X or Y derivatives
c
        IF(MYRANK.EQ.0) THEN
          open(unit=10,file='impldecomp.input',form='formatted')
          read(10,'(A)')
c
c number of implicit subdomains for the X derivatives
c
          read(10,*) NIMPLX
          IF(IMPLXDCMP>0) THEN
            read(10,'(A)')
c
c boundaries of the implicit subdomains for the X derivatives
c
            DO I=1,NIMPLX
              read(10,*) IMPLXVLIM(I,:)
              IF(IMPLXVLIM(I,1)<2) IMPLXVLIM(I,1)=2             
              IF(IMPLXVLIM(I,2)>NX-1) IMPLXVLIM(I,2)=NX-1
            ENDDO
          ENDIF
          read(10,'(A)')
c
c number of implicit subdomains for the Y derivatives
c
          read(10,*) NIMPLY
          IF(IMPLYDCMP>0) THEN
            read(10,'(A)')
c
c boundaries of the implicit subdomains for the Y derivatives
c
            DO I=1,NIMPLY
              read(10,*) IMPLYVLIM(I,:)
              IF(IMPLYVLIM(I,1)<2) IMPLYVLIM(I,1)=2
              IF(IMPLYVLIM(I,2)>NX-1) IMPLYVLIM(I,2)=NX-1
            ENDDO
          ENDIF
          close(10)
        ENDIF
        

        CALL MPI_BCAST(NIMPLX,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(NIMPLY,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)

        IF(IMPLXDCMP==0 .AND. NIMPLX>0) THEN
          ERRMSG='NIMPLX SHOULD BE ZERO '
          CALL ENDJOB(ROUTINE,ERRMSG,MYSIZE)          
         ENDIF

        IF(IMPLYDCMP==0 .AND. NIMPLY>0) THEN
          ERRMSG='NIMPLY SHOULD BE ZERO '
          CALL ENDJOB(ROUTINE,ERRMSG,MYSIZE)          
        ENDIF

        IMPLXULIM = IMPLXVLIM
        IMPLYULIM = IMPLYVLIM
        DO I=1,NIMPLX
          IF(IMPLXULIM(I,2)>NX-2) IMPLXULIM(I,2)=NX-2
        ENDDO
        DO I=1,NIMPLY
          IF(IMPLYULIM(I,2)>NX-2) IMPLYULIM(I,2)=NX-2
        ENDDO

        CALL MPI_BCAST(IMPLXULIM,2*NIMPLXMAX,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(IMPLXVLIM,2*NIMPLXMAX,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(IMPLYULIM,2*NIMPLYMAX,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(IMPLYVLIM,2*NIMPLYMAX,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)

        !Check for overlapping regions
        RVIMPLX = 0.0
        RVIMPLY = 0.0
        DO I=1,NIMPLX
          RVIMPLX(IMPLXVLIM(I,1):IMPLXVLIM(I,2)) = 1.0 
        ENDDO
        DO I=1,NIMPLY
          RVIMPLY(IMPLYVLIM(I,1):IMPLYVLIM(I,2)) = 1.0 
        ENDDO
        

        I = MAX(MAXVAL(IMPLXVLIM(:,2)),MAXVAL(IMPLYVLIM(:,2)))
        IF(ANY(RVIMPLX(1:I)+RVIMPLY(1:I)>1)) THEN
          ERRMSG='PARTIALLY IMPLICIT DOMAINS OVERLAP '
          CALL ENDJOB(ROUTINE,ERRMSG,MYSIZE)          
        ENDIF

        !Check for adjacent region and adjust index for U
        DO J=1,NIMPLY
          DO I=1,NIMPLX
            IF(IMPLXVLIM(I,1)-IMPLYVLIM(J,2)==1) THEN
              IMPLXULIM(I,1) = IMPLXULIM(I,1)-1
              IMPLYULIM(J,2) = IMPLYULIM(J,2)-1
            ENDIF
          ENDDO
c            IMPLYULIM(J,2) = IMPLYULIM(J,2)-1
        ENDDO        

        if(myrank.eq.0) then
          do i=1,nimplx
            write(6,'(A,10(1x,I4))') 'IMPLXULIM,IMPLXVLIM:',i,implxulim(i,:),implxvlim(i,:)
          enddo
        
          do i=1,nimply
            write(6,'(A,10(1x,I4))') 'IMPLYULIM,IMPLYVLIM:',i,implyulim(i,:),implyvlim(i,:)
          enddo
        endif

        RUIMPLX = 0.0
        RVIMPLX = 0.0
        DO I=1,NIMPLX
          RUIMPLX(IMPLXULIM(I,1):IMPLXULIM(I,2)) = 1.0 
          RVIMPLX(IMPLXVLIM(I,1):IMPLXVLIM(I,2)) = 1.0 
        ENDDO
        
        RUIMPLY = 0.0
        RVIMPLY = 0.0
        DO I=1,NIMPLY
          RUIMPLY(IMPLYULIM(I,1):IMPLYULIM(I,2)) = 1.0 
          RVIMPLY(IMPLYVLIM(I,1):IMPLYVLIM(I,2)) = 1.0 
        ENDDO

      ENDIF

c
c Fully implicit treatment of the X derivatives
c
      IF(IMPLXDCMP==0) THEN
        NIMPLX = 1
        IMPLXULIM(1,1) = 2
        IMPLXULIM(1,2) = NX-2
        IMPLXVLIM(1,1) = 2
        IMPLXVLIM(1,2) = NX-1
      ENDIF

c
c Fully implicit treatment of the Y derivatives
c
      IF(IMPLYDCMP==0) THEN
        NIMPLY = 1
        IMPLYULIM(1,1) = 2
        IMPLYULIM(1,2) = NX-2
        IMPLYVLIM(1,1) = 2
        IMPLYVLIM(1,2) = NX-1
      ENDIF

c      CALL MPI_FINALIZE(IERR)
c      STOP
C
      RETURN 
      END

      SUBROUTINE INPALL_POST
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - READ INPUT-DATA FROM UNIT 1
C                                                                       
C     INPUT:      - NAMELIST DIRECTED DATA FROM UNIT 1
C                                                                       
C-----------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'io.h'
      INCLUDE 'mpif.h'
      CHARACTER   TEXT(72)*1
      INTEGER     I    
     
       NAMELIST /POST/ RES1,RES2,PSTRIDE
      

C-----------------------------------------------------------------------
C                                                            READ UNIT 1
C-----------------------------------------------------------------------
C 

      IF(MYRANK==0) THEN
         OPEN (1,STATUS='OLD',FILE='post.input')
         OPEN (2,STATUS='SCRATCH')
         OPEN (3,STATUS='UNKNOWN'
     &        ,FILE='post.input.backup')
C 
C-----------------------------------------------------------------------
C                            COPY NAMELIST-DIRECTED INPUT DATA TO UNIT 2
C-----------------------------------------------------------------------

         REWIND 1 
         REWIND 2 
         REWIND 3
 10      READ (1,'(72(A1))',END=20) (TEXT(I),I=1,72)
         WRITE(3,'(72(A1))') (TEXT(I),I=1,72)
         IF (TEXT(1) .EQ. '*') GOTO 10
         WRITE(2,'(72(A1))') (TEXT(I),I=1,72)
         GOTO 10
 20      REWIND 2
      ENDIF

C 
C-----------------------------------------------------------------------
C                                    SET DEFAULT-VALUES & READ NAMELISTS
C-----------------------------------------------------------------------
C
C....."STRPRM" 
C
c
c names of the grid and geometry files
c 
      IF(MYRANK==0) THEN
         READ (2,POST) 
         WRITE(*,POST)
      ENDIF

!      CALL MPI_BCAST(STR1,60,MPI_CHARACTER,0,MPI_COMM_EDDY,IERR)
!      CALL MPI_BCAST(STR2,60,MPI_CHARACTER,0,MPI_COMM_EDDY,IERR)

	RETURN
	END 











 
