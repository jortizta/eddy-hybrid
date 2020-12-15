C-----------------------------------------------------------------------
C                 ***************************
C                 *         E D D Y         *
C                 ***************************
C-----------------------------------------------------------------------
C
C     THREEDIMENSIONAL INCOMPRESSIBLE NAVIER-STOKES SOLVER
C     (LARGE EDDY SIMULATION)
C     ----------------------------
C
C     PURPOSE:            - TIME ACCURATE NAVIER-STOKES SOLVER FOR
C                           TURBULENT FLOW FIELDS
C                         - FINITE DIFFERENCE SCHEMES
C
C     SUBGRID MODELS:     - SMAGORINSKI
C                         - DYNAMIC MODEL
C
C
C     GRID:               - CARTESIAN RECTANGULAR GRID (SINGLE BLOCK)
C
C     SCHEMES:            - ADAMS-BASHFORTH (FRACTIONAL STEP)
C                         - DIRECT POISSON SOLVER
C                         - MULTIGRID POISSON SOLVER
C
C     NON-DIMENSIONAL:    - BOUNDARY LAYER:
C                           U_INFTY, DISPLACEMENT THICKNESS
C                           DENSITY = 1.
C                         - CHANNEL:
C                           U_TAU, HALF DIAMETER
C                         - DUCT:
C                           U_TAU, HALF DIAMETER
C
C     INCLUDES:           - common.h
C
C     INPUT:              - <dir>:    ../calc/
C                         - <dir>:    ../init/
C                         - NAMELIST
C
C     OUTPUT:             - <dir>:    ../calc/
C                         - RESTART-FILE
C                         - STATISTICS
C                         - INST. FLOW FIELD
C
C
C     COMPUTING TIME:     - in a typical run the single steps take
C                           the following comp. time (Adams Bashforth):
C
C                         - RHS                  1.707024       22.2 %
C                         - velocity update      0.7407840       9.8 %
C                         - b.c.                 4.2943999E-02   0.7 %
C                         - divergence           0.2957280       3.9 %
C                         - Poisson              0.5641280       7.5 %
C                         - velocity correct.    0.6744159       8.8 %
C                         - space average        0.1864160       2.6 %
C                         - twall                2.9279999E-03   0.1 %
C                         - b.c.                 4.1967999E-02   0.6 %
C                         - strain               0.7076000       3.7 %
C                         - turvis (Smagor.)     1.554768       20.2 %
C                         - space time average   0.4792160       6.4 %
C                         - screen info          0.7446880       9.8 %
C                         - total (sec)          7.7426078
C
C-----------------------------------------------------------------------
C
C     (C) I. BALARAS
C         JIANMING YANG
C
C         DATE OF CREATION :  20/10/1993
C         LAST UPDATE      :  05/10/2000
C
C-----------------------------------------------------------------------
C
      PROGRAM EDDY
C
C-----------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
!!!!!!
      INCLUDE 'averages.h'
      INCLUDE 'fields.h'
      integer iplant2,iplant2c,iplant2d,iplant2u
      integer iplant9c,iplant9d,iplant9p,ifields
!!!!!!

c.....parameters for the immersed boundary ..............................

      include 'immersed.h'
C
C-----------------------------------------------------------------------
C                                                           declarations
C-----------------------------------------------------------------------
C
      INTEGER STATUS(MPI_STATUS_SIZE)
      INTEGER nx,ny,nz,nyg,nzg,nyl,nzl
      INTEGER id,jd,kd
      INTEGER idir
C
C.....arrays for the solver
C
C.....xu, yv, wz:		cell faces
C.....xc, yc, zc:               cell centers
C.....uo, vo, wo:		inst. velocity field
C.....us, vs, ws:        	predicted velocities
C.....p, dp:	                pressure, correction (or divergence)
C.....ua, va, wa:		RHS predictor step
C.....ub, vb, wb:		RHS corrector step
C.....tv:                       turbulent viscosity
C.....sxx,sxy,sxz,syy,syz,szz:  strain rates
C
C     allocatable arrays
      REAL, DIMENSION(:),ALLOCATABLE :: xu,yv,zw,zwg,xc,yc,zc,zcg

      REAL, DIMENSION(:,:),ALLOCATABLE :: xu_car,xv_car,yu_car,yv_car,xc_car,yc_car

      REAL, DIMENSION(:,:,:), ALLOCATABLE ::
     &     UO, VO, WO,
     &     US, VS, WS,
     &     TV, P , DP,
     &     UA, VA, WA,
     &     UB, VB, WB
     &     ,VP,KSGS,utmp1,vtmp1,wtmp1,ptmp1,denstmp1
c
      REAL, DIMENSION(:,:,:), ALLOCATABLE ::
     &     UC , VC , WC ,
     &     SXX, SYY, SZZ,
     &     SXY, SXZ, SYZ,
     &     LM, MM, ILM, IMM,
     &     G , DXDYDZ, CLES
     &    ,CLESP,CLESN
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE ::
     &     NILMP,NILMN
c
      REAL, DIMENSION(:,:,:), ALLOCATABLE ::
     &     DENS, RHA, RHB, DENSO,RS
C
      INTEGER i, j, k,kg,n,nstep
C
      INTEGER icycle
      INTEGER IS
      REAL    TIME,DTM1,TLEVEL
      REAL    CFLM,ttotal
      REAL    ALFXDT,GAMXDT,RHOXDT
C
      integer iclock,nclocks,counts,rate,nclocks2
      real clocktemp,clocktemp1,tclock
      real, dimension(:), allocatable :: CLOCK,CLOCKG,CLOCKGMIN,CLOCKGMAX
c
c.....immersed boundary
      INTEGER nbd,nfacet,nfacetmax,nfacetot,ibd,mbd,m,l,im,jm,km,jys,kzs,mbdn    !!!!!!!
      INTEGER nvtx,nvtxtot,nvtxmax
      REAL    dsbmax
      integer, dimension(:,:), allocatable  :: limu,limv,limw,limp
     &     ,                                   mimu,mimv,mimw,mimp
     &                                        ,limtv,mimtv

      integer, dimension(:), allocatable :: limu_wm,mimu_wm,limv_wm,mimv_wm,limw_wm,mimw_wm

      integer, dimension(:), allocatable  :: limuint,limvint,limwint,limpint
     &     ,                                 mimuint,mimvint,mimwint,mimpint
     &                                      ,limtvint,mimtvint

      integer, dimension(:,:), allocatable  :: limpbd,mimpbd
      integer, dimension(:,:), allocatable  :: limpnd,mimpnd
      integer, dimension(:,:), allocatable  :: limuvtx,limvvtx,limwvtx
     &                                        ,mimuvtx,mimvvtx,mimwvtx

      real, dimension(:), allocatable ::
     &     uim,xnu,ynu,znu,nxu,nyu,nzu,
     &     vim,xnv,ynv,znv,nxv,nyv,nzv,
     &     wim,xnw,ynw,znw,nxw,nyw,nzw,
     &     pim,xnp,ynp,znp,nxp,nyp,nzp,rhoim
     &    ,tvim,xntv,yntv,zntv,nxtv,nytv,nztv

      real, dimension(:), allocatable :: uwm,vwm,wwm

      Integer, dimension(:), allocatable :: iu,ju,ku,iv,jv,kv
     &     ,                                iw,jw,kw,ip,jp,kp
     &     ,                                itv,jtv,ktv

      Integer, dimension(:), allocatable :: iuint,juint,kuint
      Integer, dimension(:), allocatable :: ivint,jvint,kvint
      Integer, dimension(:), allocatable :: iwint,jwint,kwint
      Integer, dimension(:), allocatable :: ipint,jpint,kpint
      Integer, dimension(:), allocatable :: itvint,jtvint,ktvint

      INTEGER, DIMENSION(:), ALLOCATABLE :: fp,fpu,fpv,fpw,fptv,trino,mbg    !!!!!!

      Integer, dimension(:), allocatable :: mrku,mrkv,mrkw,mrkp,mrktv
      Integer, dimension(:), allocatable :: diru,dirv,dirw,dirp,dirtv
      Integer, dimension(:), allocatable :: mrkuwm,mrkvwm,mrkwwm

      real, dimension(:,:,:), allocatable :: vertex
      real, dimension(:,:), allocatable :: unvect,vertexc,node

      integer, dimension(:), allocatable :: mrkpb,mrkcf
      integer, dimension(:,:), allocatable :: mrksb
      real, dimension(:), allocatable :: pbd,dudxb,dudyb,dudzb
     &                ,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb,areaf,cf
      real, dimension(:,:), allocatable :: dsb,cf_aux

      real, dimension(:,:), allocatable :: trvtx,vtxnrm
      integer, dimension(:), allocatable :: mrkpvtx,vtxno
      integer, dimension(:,:), allocatable :: mrksvtx
      real, dimension(:), allocatable :: pvtx,dudxvtx,dudyvtx,dudzvtx
     &                ,dvdxvtx,dvdyvtx,dvdzvtx,dwdxvtx,dwdyvtx,dwdzvtx

      real, dimension(:,:,:), allocatable :: umtrx,vmtrx,wmtrx,pmtrx,tvmtrx
      Integer, dimension(:,:), allocatable :: uindx,vindx,windx,pindx,tvindx
C
      Integer, dimension(:,:,:,:), allocatable :: flaguo,flagvo,flagwo,flagpo,flagpbd,flagtv

      Integer, dimension(:,:,:), allocatable :: flag
c
      REAL, DIMENSION(:), ALLOCATABLE :: dpdnn,drhodnn

c.... Cf probes
      integer njcf,nkcf,nkcfprev,nkcfmax,icfprb,itcfprb,ncfprb,ncfsmpl
      character*160 filecfprb
      integer, dimension(:), allocatable :: jcf,kcf
      real, dimension(:,:,:), allocatable :: cfprb
      real, dimension(:), allocatable :: tcf
c
c.... Average Velocity at Top Boundary
      integer flag_wavg_top,freq_wavg_top,wavg_top_nsmpl
      real, dimension(:,:), allocatable :: wavg_top
      logical file_exists
      character*160 file_wavg_top

      real cvlim1(4),cvlim2(4),cvlim3(4)
      real ubd,vbd,wbd,ddn,dpds,dpdxx2,dpdyy2,dpdzz2,dpdxx3,dpdyy3,dpdzz3   !vbdm,vbdm2
      real interp_cellface
      integer nimu,nimv,nimw,nimp,nimtv,ilb,ile,icom
!!!!!!     integer nflumax,nflvmax,nflwmax,nflpmax,nflpbdmax,nfltvmax
      integer, dimension(:), allocatable :: nflumax,nflvmax,nflwmax,nflpmax,nflpbdmax,
     %nfltvmax

      real a(nsp,nsp),b(nsp),c(nsp)
      integer indx(nsp)
      integer ism

      integer iflag

      integer imn,imx,jmn,jmx,kmn,kmx

      real    angle,radius,anglerad
      real    fxim,fyim,fzim,fx,fy,fz,rbmax
      REAL    FB(3),FB2(3),FB3(3)
      real tloct,tmove,pref,temp(2)

      REAL RAN1,VPER
      INTEGER IDUM
      integer ii,iid,iiu,jj,kk,i0,j0,k0,i1,j1,k1,k1g,i2,j2,k2
      INTEGER iflag_psb,iflag_ptr,iflag_pnd
      integer NSTAR,imarker,IND(3,6)
      real distance,cellsize,sext,area,areaij,qin(2),uwall,dudyw
      REAL TIME_ITERATION,TIME_ITERATIONS,TIME_LAST_ITERATION
      real XTEMP,ZTEMP
      real uvel_taylor_green,wvel_taylor_green,pres_taylor_green
c
c.... Time probe
      integer ntprb,ntprbmax,tmprbindx,nprbser,iprbser,iapnd,ntprbsmpl
c      real     xprb,yprb,zprb
      integer, dimension(:), allocatable :: itprb,jtprb,ktprb
      real   , dimension(:,:), allocatable :: utprb,vtprb,wtprb,ptprb
      real   , dimension(:), allocatable :: xprb,yprb,zprb
c
c.... Azimuthal spectra and autocorrelation
      integer nyprb,nyprbmax,yprbindx,nycorvar,nyesdvar,itysamp,yesdsmpl,ycorsmpl
      real    ns,invs
      integer, dimension(:), allocatable :: iyprb,jyprb,kyprb
      real, dimension(:,:,:), allocatable :: esdprb,esd
      real, dimension(:,:,:), allocatable :: ycorprb,ycor
c
c.... Time average sampling of prime variables
      integer limtmavg(8),ftmavg,nxt,nyt,nzt
      integer nvarreg,nreg,fregtmavg,kmintmavg,kmaxtmavg,nregtot,nvartmavg
      integer nzreg,nzgreg,nregprev,zindreg,nsmplregtmavg,j1reg,j2reg,njreg
      integer, dimension(:), allocatable :: nsmpltmavg
      integer, dimension(:,:), allocatable :: iindtmavg
      real, dimension(:,:), allocatable :: regtmavg
      real, dimension(:,:,:), allocatable :: utmavg,vtmavg,wtmavg,ptmvavg
      real*4, dimension(:,:,:,:), allocatable :: vartmavg
c
c.... Variables for writing VPfield
      integer limVPfield(8),fVPfield,nxVPf,nyVPf,nzVPf
      integer kminVPf,kmaxVPf
c
c.... Blasius solution
      integer neta
      real  xo,zo
      real  deta,eta1,eta2,ueta1,ueta2,veta1,veta2
      real, dimension(:), allocatable :: eta,f1eta,f2eta
c
c.... Variable for region VPfield
      logical VPreg_input_exists
      character*150 VPreg_ind_file,fileVPreg
      integer kminVPreg,kmaxVPreg,nVPreg,nzVPreg,nzgVPreg,nprevVPreg,zindVPreg
     &     ,j1VPreg,j2VPreg,njVPreg,VPreg_freq,iVPreg,nVPregtot
      integer, dimension(:,:), allocatable :: iindVPreg
c
c.... Variables for average CFL
      integer nsmplcfl,icflave
      logical cfl_input_exists
      character*200 filecflave
      real, dimension(:,:), allocatable :: CFLAVE
c
c.... Variables for stat1D
      integer istat1D,itstat1D,nstat1D,nstat1Dave,nvarstat1D
      logical stat1D_input_exists
      real, dimension(:,:), allocatable :: stat1Dave,stat1D
c
c.... Pipe flow
      real    ubold
      integer iutau,itutau,ntutau,itu,utaudim
      logical utau_input_exists
      real, dimension(:), allocatable :: timeutau,ubulk,dpdzutau
      real, dimension(:,:), allocatable :: utau
      real, dimension(:), allocatable :: cfg
      real    vecmag
      real    unorm,a1,a2,a3
      real    sxxb,syyb,szzb,sxyb,sxzb,syzb

      logical flag_top_pos
      real    osc_freq
      character*160 hspdir,filename,filename1,fileregind,ICfile
      character*160 filetprb,fileycor,fileyspc,fileprbtmp,filetmavg,filereg

      character*160 icycle_string

      real,allocatable,dimension(:,:) :: planeU,planeV,planeW,planeD
      real,allocatable,dimension(:,:) :: planeURef,planeVRef,planeWRef,planeDRef
      character*600 filenameRef 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                 directories


c-----------------------------------------------------------------------
c                                                 input data and setup
c-----------------------------------------------------------------------

      flag_top_pos=.false.
      CALL MPI_INIT(IERR)

      MPI_COMM_EDDY=MPI_COMM_WORLD
c...start mpi, get rank and number of processors
      CALL MPI_COMM_SIZE(MPI_COMM_EDDY,MYSIZE,IERR)
      CALL MPI_COMM_RANK(MPI_COMM_EDDY,MYRANK,IERR)
c
c...input parameters (namelist-statements)
c

      begin_time=tclock()
      save_res = .FALSE.
      ttotal=1.0

      CALL INPALL(nx,ny,nz,nzg)

      hspdir = './'
      IF(MYRANK.EQ.0) write(6,'(A)') 'hspdir='//trim(hspdir)

      CALL MPI_SETUP(nx,ny,nz,nzg,ierr)
      MYLEFT  = MYRANK-1
      MYRIGHT = MYRANK+1
      IF(MYRIGHT==MYSIZE) MYRIGHT = -1
c      IF(MYRANK==0       .AND.ITYPE(5)==500) MYLEFT  = MYSIZE-1
c      IF(MYRANK==MYSIZE-1.AND.ITYPE(6)==500) MYRIGHT = 0

!
!     valve case: refresh flag at inlet and outlet through periodicity
!
      IF(MYRANK==0       ) MYLEFT  = MYSIZE-1
      IF(MYRANK==MYSIZE-1) MYRIGHT = 0
C
      IF(MYRANK==0) THEN
        WRITE(6,'(A,I4)') '*.. MYSIZE  = ',mysize
        WRITE(6,'(A,I4)') '*.. MYRANK  = ',myrank
        WRITE(6,'(A,I4)') '*.. MYLEFT  = ',myleft
        WRITE(6,'(A,I4)') '*.. MYRIGHT = ',myright
        WRITE(6,'(A)') '*..input data read'
      ENDIF
c
c...  Allocate clock array
c
      NCLOCKS=100
      ALLOCATE(CLOCK(NCLOCKS),CLOCKG(NCLOCKS),CLOCKGMIN(NCLOCKS),CLOCKGMAX(NCLOCKS))
      CLOCK=0.0
c
c...  fill index array to create a series of restart files
c
      m=0
      do i=0,9
      do j=0,9
      do k=0,9
      do l=0,9
        index(m)=CHAR(i+48)//CHAR(j+48)//CHAR(k+48)//CHAR(l+48)
        m=m+1
      enddo
      enddo
      enddo
      enddo
c
      do m=0,999
        proc(m)=index(m)(2:4)
      enddo
c
c...constants for info output
c
      kd=(nz-2)*MYRANK

      do j=1,9
        jd=10**j
        if(mod(nzg,jd)==nzg) goto 200
      enddo
 200  continue
      do i=1,9
        id=10**i
        if(mod(ny ,id)==ny ) goto 100
      enddo
 100  continue
c
      id=id*jd
c
C     ALLOCATE ARRAYS
c
      clock(1) = tclock()

      ALLOCATE(XU(NX),YV(NY),ZW(NZ),ZWG(NZG))
      ALLOCATE(XC(NX),YC(NY),ZC(NZ),ZCG(NZG))

      ALLOCATE(XU_CAR(NX,NY),YU_CAR(NX,NY),XV_CAR(NX,NY),YV_CAR(NX,NY),XC_CAR(NX,NY),YC_CAR(NX,NY))
C
      ALLOCATE(UO(NX,NY,NZ),VO(NX,NY,NZ),WO(NX,NY,NZ))
      ALLOCATE(US(NX,NY,NZ),VS(NX,NY,NZ),WS(NX,NY,NZ))
      ALLOCATE(P (NX,NY,NZ),DP(NX,NY,NZ))
      ALLOCATE(UA(NX,NY,NZ),VA(NX,NY,NZ),WA(NX,NY,NZ))
      ALLOCATE(UB(NX,NY,NZ),VB(NX,NY,NZ),WB(NX,NY,NZ))
      ALLOCATE(utmp1(NX,NY,NZ),vtmp1(NX,NY,NZ),wtmp1(NX,NY,NZ),ptmp1(NX,NY,NZ),denstmp1(NX,NY,NZ))
      ALLOCATE(TV(NX,NY,NZ),KSGS(NX,NY,NZ))
      IF(ISGS>0 .AND. ISGS<5) THEN
        ALLOCATE(UC (NX,NY,NZ),VC (NX,NY,NZ),WC (NX,NY,NZ))
        ALLOCATE(SXX(NX,NY,NZ),SYY(NX,NY,NZ),SZZ(NX,NY,NZ))
        ALLOCATE(SXY(NX,NY,NZ),SXZ(NX,NY,NZ),SYZ(NX,NY,NZ))
        ALLOCATE(DXDYDZ(NX,NY,NZ))
        ALLOCATE(G(NX,NY,NZ),LM(NX,NY,NZ),MM(NX,NY,NZ))
        ALLOCATE(ILM(NX,NY,NZ),IMM(NX,NY,NZ))
        ALLOCATE(CLES(NX,NY,NZ))
        ALLOCATE(CLESP(NX,NY,NZ),CLESN(NX,NY,NZ))
        ALLOCATE(NILMP(NX,NY,NZ),NILMN(NX,NY,NZ))
      ENDIF
        ALLOCATE(RHA(NX,NY,NZ),RHB(NX,NY,NZ),RS(NX,NY,NZ))
        ALLOCATE(DENS(NX,NY,NZ))
        allocate(planeU(nx,ny),planeV(nx,ny),planeW(nx,ny),planeD(nx,ny))
        allocate(planeURef(nx,ny),planeVRef(nx,ny),planeWRef(nx,ny),planeDRef(nx,ny))

      clock(1) = tclock() - clock(1)
c
c...grid and matrix coefficients
c
      clock(2) = tclock()

      IF(MYRANK==0) WRITE(6,*) 'Reading grid'
      CALL Grid(xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nzg)
      IF(MYRANK==0) WRITE(6,*) 'Grid read'
c
c.....Cartesian grid coordinates
c
      if(icyl==1) then
        do j=1,ny
        do i=1,nx
          xc_car(i,j) = rp(i)*cos(yc(j))
          yc_car(i,j) = rp(i)*sin(yc(j))
          xu_car(i,j) = ru(i)*cos(yc(j))
          yu_car(i,j) = ru(i)*sin(yc(j))
          xv_car(i,j) = rp(i)*cos(yv(j))
          yv_car(i,j) = rp(i)*sin(yv(j))
        enddo
        enddo
      else
        do j=1,ny
        do i=1,nx
c           rp(i) = 1.0

          xc_car(i,j) = xc(i)
          yc_car(i,j) = yc(j)
          xu_car(i,j) = xu(i)
          yu_car(i,j) = yc(j)
          xv_car(i,j) = xc(i)
          yv_car(i,j) = yv(j)
        enddo
        enddo
      endif

      clock(2) = tclock() - clock(2)

      IF(MYRANK==0) THEN
        WRITE(6,'(A)') '*..grid read'
        write(6,'(A,I5,A,I5)') '*..ix1 =',ix1,' ix2 =',ix2
        write(6,'(A,I5,A,I5)') '*..jy1 =',jy1,' jy2 =',jy2
        write(6,'(A,I5,A,I5)') '*..kz1 =',kz1,' kz2 =',kz2
      ENDIF
c
!!!!!!
      clock(39) = tclock()
c
c.....Setup for output
      ifields=0
      call setup_indices(ny)
      call read_sections(nz)
      call read_sections_plot2D(nz)
      call read_probes_all(nz)
!      call allocate_output(nx,ny,nz)
      call allocate_output_xyz(nx,ny,nz,icyl)
!      if(ifield.eq.0) then
!         call setup_output(nx,ny,nz)
!      else
         call inirea(nx,ny,nz)
!      endif
      iplant2=0
c
      clock(39) = tclock() - clock(39)
!!!!!!
c
c.....solver setup
c
      clock(3) = tclock()
c
c setup the starting and the ending indices for each boundary and
c for u, v and w momentum
c
      CALL SETUP
      if (myrank.eq.0) then
        WRITE(6,'(A)') '*..solver setup'
      endif

      clock(3) = tclock() - clock(3)

c
c.....initialize field
c
      clock(4) = tclock()

      time = 0.

      UO = 0.
      VO = 0.
      WO = 0.
      US = 0.
      VS = 0.
      WS = 0.
      P  = 0.
      TV = 0.
c
      UA = 0.
      VA = 0.
      WA = 0.
      UB = 0.
      VB = 0.
      WB = 0.

      planeU = 0.
      planeV = 0.
      planeW = 0.
      planeD = 0. 

c
c.....begin from scratch
c
      IF(ifield==0) THEN


      call hybrid_init(uo,vo,wo,dens,xc,yc,nx,ny,nz,
     & planeU,planeV,planeW,planeD)

        IF (MYRANK==256) THEN

        WRITE(*,*)'WO(100,65,5)',WO(100,65,5)

        ENDIF


	CALL BOUNDARY(UO,VO,WO,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TIME)  
	CALL BOUNDARY_DENS(DENS,XC,YC,NX,NY,NZ)

	IF(MP==0) THEN
          P(1 ,:,:) = P(IX2,:,:) !Neumann
          P(NX,:,:) = P(IX1,:,:)
        ELSE
          P(1 ,:,:) = -P(IX1,:,:) !Dirichlet
          P(NX,:,:) = -P(IX2,:,:)
        ENDIF

	IF(ITYPE(1)==300) THEN
          DO J=1,NY
            P(1,J,:) =  P(IX1,JSYM(J),:)
          ENDDO
        ENDIF
C-----------------------------------------------------------------------
c     Cyclic b.c. in y direction
C-----------------------------------------------------------------------
        IF(LP==0) THEN
          P(:,1 ,:) = P(:,JY2,:)
          P(:,NY,:) = P(:,JY1,:)
        ELSE
          P(:,1 ,:) = P(:,JY1,:)
          P(:,NY,:) = P(:,JY2,:)
        ENDIF
C-----------------------------------------------------------------------
c     Neumann b.c. in z direction
c N.B. ITYPE(5)/=0 only for the process 0
c N.B. ITYPE(6)/=0 only for the process MYSIZE-1
C-----------------------------------------------------------------------
        IF(ITYPE(5)==500) P(:,:,1 ) = P(:,:,KZ2)
        IF(ITYPE(6)==500) P(:,:,NZ) = P(:,:,KZ1)

        IF(ITYPE(5)/=0.AND.ITYPE(5)/=500) P(:,:,1 ) = P(:,:,KZ1)
        IF(ITYPE(6)/=0.AND.ITYPE(6)/=500) P(:,:,NZ) = P(:,:,KZ2)
!	CALL SPONGE_SETUP(UO,VO,WO,DENS,P,NX,NY,NZ,XC,IERR)

        CALL SPONGE_SETUP(NX,NY,NZ,NZG,XC,XU,ZWG,ZCG,IERR)

!         if(idens.eq.1) then
!           do i=1,nx
!             dens(i,:,:) = drdx*xc(i)
!           enddo
!         endif

c.... Taylor Green vortex
c        if(.false.) then
c          do k=1,nz
c          do j=1,ny
c          do i=1,nx
c            uo(i,j,k) = uvel_taylor_green((xu(i)-xmin),(yc(j)-ymin),0.0,ru1)
c            vo(i,j,k) = wvel_taylor_green((xc(i)-xmin),(yv(j)-ymin),0.0,ru1)
cc            vo(i,j,k) = uvel_taylor_green((yv(j)-ymin),(xc(i)-ymin),0.0,ru1)
cc            uo(i,j,k) = wvel_taylor_green((yc(j)-xmin),(xu(i)-ymin),0.0,ru1)
c            p(i,j,k)  = pres_taylor_green((xc(i)-xmin),(yc(j)-zmin),0.0,ru1)
c          enddo
c          enddo
c          enddo
c
c        endif
c
c.... Parabolic
c        if(.false.) then
c
c          wo = 0.0

c$$$             call add_random_noise(uo(:,:,:),nx,ny,nz,10.0)
c$$$             call add_random_noise(vo(:,:,:),nx,ny,nz,10.0)
c$$$             call add_random_noise(wo(:,:,:),nx,ny,nz,10.0)

c          do i=1,nx
c            uo(i,:,:) = uo(i,:,:)*sin( (xu(i)-xu(1))/xlen*2.0*pi)
c            vo(i,:,:) = vo(i,:,:)*sin( (xc(i)-xu(1))/xlen*2.0*pi)
c            wo(i,:,:) = wo(i,:,:)*sin( (xc(i)-xu(1))/xlen*2.0*pi)
c          enddo
c          call locate(xc,nx,-0.2,i1)
c          call locate(xc,nx, 0.2,i2)
c          uo(i1:i2,:,:) = 0.0
c          vo(i1:i2,:,:) = 0.0
c          wo(i1:i2,:,:) = 0.0
c
c... Channel laminar parabolic profile
c          do i=1,nx
c            wo(i,:,:) = wo(i,:,:) + 1.5*(1-(2.0*xc(i)/xlen)**2.0)
c          enddo
c        endif
c
c        if(.true.) then
cc          do i=1,nx
cc            wo(i,:,:) = 2.0*(1-(xc(i)/xu(nx-1))**2.0)
cc          enddo
c
c          call turb_pipe_vel(uo,vo,wo,xu,xc,nx,ny,nz)
cc          call turb_chan_vel(uo,vo,wo,xu,xc,nx,ny,nz,0)
c
c        endif
c
c...... Blasius
c        if(.false.) then
c          deta = 0.05
c          xo = 0.0
c          zo = dummyfp3
c          if(ibm>0) xo = dummyfp4
c          eta1 = (xu(nx)-xo)*sqrt(0.5/ru1/(zo+zc(1)-zmin))
c          eta2 = (xu(nx)-xo)*sqrt(0.5/ru1/(zo+zw(nz)-zmin))
c          neta = ceiling(max(eta1,eta2)/deta)+2
c          allocate(eta(neta),f1eta(neta),f2eta(neta))
c          eta = 0.0
c          do i=2,neta
c            eta(i) = deta*real(i-1)
c          enddo
c          call falkner_skan(eta,f1eta,f2eta,neta,dummyfp1,dummyfp2)
c
c          i1 = 1
c          if(ibm>0) call locate(xc,nx,xo,i1)
c          i1 = i1+1
c          do k=1,nz
c          do i=i1,nx
c            eta2 = (xc(i)-xo)*sqrt(0.5/ru1/(zo+zw(k)-zmin))
c            call locate(eta,neta,eta2,m)
c            ueta1 = 1.0*f2eta(m)
c            ueta2 = 1.0*f2eta(m+1)
c            wo(i,:, k ) = ueta1 + (ueta2-ueta1)/deta*(eta2-eta(m))
c
c            eta2 = (xu(i)-xo)*sqrt(0.5/ru1/(zo+zc(k)-zmin))
c            call locate(eta,neta,eta2,m)
c            veta1 = sqrt(0.5*ru1/(zo+zc(k)-zmin))*(eta(m)*f2eta(m)-f1eta(m))
c            veta2 = sqrt(0.5*ru1/(zo+zc(k)-zmin))*(eta(m+1)*f2eta(m+1)-f1eta(m+1))
c            uo(i,:, k ) = veta1 + (veta2-veta1)/deta*(eta2-eta(m))
c          enddo
c          enddo
c
c        endif
c
C-----------------------------------------------------------------------
C.....refresh block interfaces (and periodic boundary in z direction)
C-----------------------------------------------------------------------
c.....boundary conditions

c
c.....information transfer to ghost layers
c
        CALL REFRESHBC(UO,NX*NY,NZ)
        CALL REFRESHBC(VO,NX*NY,NZ)
        CALL REFRESHBC(WO,NX*NY,NZ)

        if(icfl==1) then
c.....courant number calculation for running at cfl=constant
c.....the implicit discretization of the convective and viscous terms is
c.....taked into account
          CALL CALCFL(UO,VO,WO,TV,DP,NX,NY,NZ,CFLM,0)
          DTM1=CFLC/CFLM
          TSTEP=DTM1
          if(tstep.le.1.e-06) then
             call mpi_finalize(ierr)
             stop
          endif
        ELSE
          DTM1=TSTEP
        ENDIF
        ALFXDT=ALF(1)*DTM1
c
c.....definition of the inflow and outflow conditions
c
        CALL BOUNDINOUT(UO,VO,WO,UO,VO,WO,XU,XC,YC,ZW,ZC,ALFXDT,TIME,
     & NX,NY,NZ,1,planeU,planeV,planeW)
	CALL BOUNDARY(UO,VO,WO,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TIME)  	    
        CALL BOUNDINOUTD(DENS,WO,ZCG,ALFXDT,NX,NY,NZ,NZG,planeD)

        IF(ISCHM==1) THEN
          CALL RHS(UO,VO,WO,DENS,US,VS,WS,TV,UB,VB,WB,DP,NX,NY,NZ,XC,XU,YC,YV,TIME)
c
c US,VS,WS: terms treated by C.N.
c UB,VB,WB: all other terms of the RHS
c
          UA=UB
          VA=VB
          WA=WB
c
!           IF(IDENS.EQ.1) THEN
            CALL RHS_DENSITY(UO,VO,WO,DENS,TV,RHB,RS,NX,NY,NZ,YC,YV)
            RHA=RHB
!           ENDIF
c
        ENDIF
        nstep=0
c.....start from saved field
c.....it is necessary to modify the name of the restart file from
c.....(variable).res to res.ini.(variable)
c
      ELSE
!**************************************************************************************************
! Original
	call add_density(xc,yc,nx,ny,nz,dens,ierr)
        CALL SPONGE_SETUP(NX,NY,NZ,NZG,XC,XU,ZWG,ZCG,IERR)
        IDIR  = -1 ! IOMPI_3DSCALAR reads restart files
        IF(MYRANK.EQ.0) write(6,*) 'Reading u restart file ...'
        write(ICfile,'(a,i8.8,a)')trim(hspdir)//"u_",resstep,".res"
        CALL IOSCALAR(ICfile,UO,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
c        CALL IOMPI_3DSCALAR(trim(hspdir)//'res.ini.u',UO,NX,NY,NZ,IDIR)
!         CALL HDF5_MPI_3DREAL(trim(hspdir)//'res.ini.u',UO,NX,NY,NZ,IDIR)
        IF(MYRANK.EQ.0) write(6,*) 'Reading v restart file ...'
        write(ICfile,'(a,i8.8,a)')trim(hspdir)//"v_",resstep,".res"
        CALL IOSCALAR(ICfile,VO,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
c        CALL IOMPI_3DSCALAR(trim(hspdir)//'res.ini.v',VO,NX,NY,NZ,IDIR)
!         CALL HDF5_MPI_3DREAL(trim(hspdir)//'res.ini.v',VO,NX,NY,NZ,IDIR)
        IF(MYRANK.EQ.0) write(6,*) 'Reading w restart file ...'
        write(ICfile,'(a,i8.8,a)')trim(hspdir)//"w_",resstep,".res"
        CALL IOSCALAR(ICfile,WO,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
c        CALL IOMPI_3DSCALAR(trim(hspdir)//'res.ini.w',WO,NX,NY,NZ,IDIR)
!         CALL HDF5_MPI_3DREAL(trim(hspdir)//'res.ini.w',WO,NX,NY,NZ,IDIR)
        IF(MYRANK.EQ.0) write(6,*) 'Reading p restart file ...'
        write(ICfile,'(a,i8.8,a)')trim(hspdir)//"p_",resstep,".res"
        CALL IOSCALAR(ICfile,P ,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
c        CALL IOMPI_3DSCALAR(trim(hspdir)//'res.ini.p',P,NX,NY,NZ,IDIR)
!         CALL HDF5_MPI_3DREAL(trim(hspdir)//'res.ini.p',P,NX,NY,NZ,IDIR)

        IF(IDENS.EQ.1) then
        write(ICfile,'(a,i8.8,a)')trim(hspdir)//"dens_",resstep,".res"
        CALL IOSCALAR(ICfile,DENS ,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
!         CALL HDF5_MPI_3DREAL(trim(hspdir)//'res.ini.dens',DENS,NX,NY,NZ,IDIR)
        ENDIF

        write(filenameRef,'(2a,i4.4,a,i8.8,a)')
     &  trim(path_feed),"/hybrid_k",
     &  index_feed,"_n", end_feed,".interp"

        call read_hybrid(filenameRef,planeURef,planeVRef,planeWRef,planeDRef,nx,ny)

        call sequence_hybrid(planeU,planeV,planeW,planeD,
     &  planeURef,planeVRef,planeWRef,planeDRef,
     &  nstep-1,nx,ny)         

        


! If restarting after interpolation add
!        CALL BOUNDARY_DENS(DENS,XC,YC,NX,NY,NZ)
!        CALL REFRESHBC(DENS,NX*NY,NZ)
!
!        do k=kz1,kz2
!        do j=jy1,jy2
!        do i=ix1,ix2
!
!        dens(i,j,k)   =  dens(i,j,k)+denP1*rp(i)*sin(yc(j))
!
!        enddo
!        enddo
!        enddo
!
!        CALL BOUNDARY_DENS(DENS,XC,YC,NX,NY,NZ)
!        CALL REFRESHBC(DENS,NX*NY,NZ)
!
!        CALL REFRESHBC(UO,NX*NY,NZ)
!        CALL REFRESHBC(VO,NX*NY,NZ)
!        CALL REFRESHBC(WO,NX*NY,NZ)
!
!        CALL BOUNDARY(UO,VO,WO,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TLEVEL)


! Karu adds ------------------------------------
        CALL BOUNDARY_DENS(DENS,XC,YC,NX,NY,NZ)
        CALL REFRESHBC(DENS,NX*NY,NZ)
!------------------------------------------------
        TSTEP=DTM1

        IF(MYRANK.EQ.0) write(6,*) 'Restart files read. time=',time

        IF(ISCHM==1) THEN
          IF(IFIELD>1) THEN
c
c The RHS for the predictor step is read from file in the case ifield > 1
c
c            CALL IOMPI_3DSCALAR(trim(hspdir)//'res.ini.ua',UA,NX,NY,NZ,IDIR)
c            CALL IOMPI_3DSCALAR(trim(hspdir)//'res.ini.va',VA,NX,NY,NZ,IDIR)
c            CALL IOMPI_3DSCALAR(trim(hspdir)//'res.ini.wa',WA,NX,NY,NZ,IDIR)
!            CALL HDF5_MPI_3DREAL(trim(hspdir)//'res.ini.ua',UA,NX,NY,NZ,IDIR)
!            CALL HDF5_MPI_3DREAL(trim(hspdir)//'res.ini.va',VA,NX,NY,NZ,IDIR)
!            CALL HDF5_MPI_3DREAL(trim(hspdir)//'res.ini.wa',WA,NX,NY,NZ,IDIR)
c            CALL IOSCALAR('res.ini.ua',UA,DP,NX,NY,NZ,IDIR,TIME)
c            CALL IOSCALAR('res.ini.va',VA,DP,NX,NY,NZ,IDIR,TIME)
c            CALL IOSCALAR('res.ini.wa',WA,DP,NX,NY,NZ,IDIR,TIME)
!            IF(IDENS.EQ.1) THEN
!            CALL HDF5_MPI_3DREAL(trim(hspdir)//'res.ini.rha',RHA,NX,NY,NZ,IDIR)
            ELSE
            CALL RHS(UO,VO,WO,DENS,US,VS,WS,TV,UB,VB,WB,DP,NX,NY,NZ,XC,XU,YC,YV,TIME)
            UA=UB
            VA=VB
            WA=WB
!             IF(IDENS.EQ.1) THEN
              CALL RHS_DENSITY(UO,VO,WO,DENS,TV,RHB,RS,NX,NY,NZ,YC,YV)
              RHA=RHB
! Karu adds ------------------------------------
              
              CALL BOUNDARY_DENS(RHA,XC,YC,NX,NY,NZ)
              CALL BOUNDARY_DENS(RHB,XC,YC,NX,NY,NZ)
              CALL BOUNDARY_DENS(RS,XC,YC,NX,NY,NZ)

              CALL REFRESHBC(RHA,NX*NY,NZ)
              CALL REFRESHBC(RHB,NX*NY,NZ)
              CALL REFRESHBC(RS,NX*NY,NZ)
! ----------------------------------------------

!             ENDIF
          ENDIF
        ENDIF

!        IF(ISGS>0) THEN
!c          CALL IOMPI_3DSCALAR(trim(hspdir)//'res.ini.tv',TV,NX,NY,NZ,IDIR)
!          CALL HDF5_MPI_3DREAL(trim(hspdir)//'res.ini.tv',TV,NX,NY,NZ,IDIR)
!        ENDIF

        IF(ISGS==4) THEN
c
c in the case the simulation is restarted from a field given by a
c LES with the Lagrangian SGS model
c
c          CALL IOMPI_3DSCALAR(trim(hspdir)//'res.ini.ilm',ILM,NX,NY,NZ,IDIR)
c          CALL IOMPI_3DSCALAR(trim(hspdir)//'res.ini.imm',IMM,NX,NY,NZ,IDIR)
!          CALL HDF5_MPI_3DREAL(trim(hspdir)//'res.ini.ilm',ILM,NX,NY,NZ,IDIR)
!          CALL HDF5_MPI_3DREAL(trim(hspdir)//'res.ini.imm',IMM,NX,NY,NZ,IDIR)
c          CALL IOSCALAR('res.ini.ilm',ILM,DP,NX,NY,NZ,IDIR,TIME)
c          CALL IOSCALAR('res.ini.imm',IMM,DP,NX,NY,NZ,IDIR,TIME)
        ENDIF

!         IF(IFIELD>1) THEN
!           open(unit=10,file='res.time',form='formatted')
!           read(10,*) time
!           close(10)
!         ELSE
!           TIME = TINI
!         ENDIF

      END IF

      clock(4) = tclock() - clock(4)

      TLEVEL=TIME
      if(imedie.eq.1)timep=time    !!!!!!

c        if(.false.) then
c        call pipe_turb_prof(dp(:,1,1),xc,nx)
c
c-----------------------------------------------------------------------
c     setup external forcing
c-----------------------------------------------------------------------
      if(ibm/=0) then
c
c.....open file to write immersed boundary information
c
        clock(5) = tclock()
        IF(IOLVL>0 .AND. MYRANK.EQ.0) THEN
          OPEN(UNIT=16,FILE='stats_imb.dat',FORM='FORMATTED'
!!!!!!     &      ,STATUS='REPLACE')
     &,POSITION='APPEND')
          CLOSE(16)
        ENDIF
c
c.....read stl file
c
        if(myrank==0) then
          write(6,'(A)') '*..immersed solids in stl format'
          open(unit=19,file=str5,form='formatted',status='old')
          read(19,*) nbd
          write(6,'(A,I5)') '*..number of solids: ',nbd
          do ibd=1,nbd
            read(19,*) solid(ibd)
            write(6,'(A,I5,3X,A)') '*..solid no. ',ibd,solid(ibd)
            read(19,*) mb(ibd)
            write(6,'(A,I8)') '*..triangles in this solid: ',mb(ibd)
c            read(19,*) mnd(ibd)
c            write(6,'(A,I8)') '*..nodes in this solid: ',mnd(ibd)
            IF(ITBDY>0 .AND. IVRTX==1) THEN
              read(19,*) mv(ibd)
              write(6,'(A,I8)') '*..vertices in this solid: ',mv(ibd)
            ENDIF
            read(19,*) mrb(ibd)
            write(6,'(A,I3)') '*..this solid does rigid body rotation?'
     &          //' (yes - 1 or 2; no - 0)',mrb(ibd)
            read(19,*) crb(1,ibd),crb(2,ibd),crb(3,ibd)
            write(6,'(A,3G18.6)') '*..the axis of the rotation is '
     &           ,crb(1,ibd),crb(2,ibd),crb(3,ibd)
          enddo
          close(19)
        endif

        CALL MPI_BCAST(NBD,1  ,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(MB ,NBD,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(MRB,NBD,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(CRB,NBD*3,MTYPE,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(SOLID,NBDMAX*80,MPI_CHARACTER,0,MPI_COMM_EDDY,IERR)

        lb(1) = 0
        do ibd=2,nbd
          lb(ibd) = lb(ibd-1)+mb(ibd-1)
        enddo
c
c definition of the total number of faces of the immersed boundaries
c
        nfacet = sum(mb(1:nbd))
        nfacetot = nfacet

        IF(ITBDY>0 .AND. IVRTX==1) THEN
          CALL MPI_BCAST(MV ,NBD,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
          lv(1) = 0
          do ibd=2,nbd
            lv(ibd) = lv(ibd-1)+mv(ibd-1)
          enddo
          nvtx = sum(mv(1:nbd))
          nvtxtot = nvtx
          nvtxmax = nvtx
        ENDIF

        IF(IMBDCMP==1 .AND. MYSIZE>1) THEN
          nfacet = sum(mb(1:nbd))
          nvtx = sum(mv(1:nbd))
        ENDIF

!c
!c N.B. The moving boundaries have to be provided after the stationary ones
!c From 1 to mbd-1 there are the stationary boundaries; from mbd to nbd the
!c moving ones
!c
!!!!!!        mbd = nbd+1
!!!!!!        do ibd=1,nbd
!!!!!!          if(mrb(ibd)>0) then
!!!!!!            mbd=ibd
!!!!!!            if(any(mrb(mbd:nbd)/=1)) then
!!!!!!c              write(8,'(A)') '*..wrong parameter for stl.input!'
!!!!!!            endif
!!!!!!c            write(8,*) 'those are moving parts: ',(i,i=mbd,nbd)
!!!!!!            exit
!!!!!!          endif
!!!!!!        enddo
c
c From 1 to mbdn-1 there are the stationary boundaries; from mbdn to mbd-1 the
c moving constant ones and from mbd to nbd the moving variable ones
c
        mbdn = nbd+1
        mbd = nbd+1
        do ibd=1,nbd
          if(mrb(ibd).eq.1) then
            mbdn=ibd
            if(any(mrb(mbdn:nbd).eq.0)) then
               write(6,*)'wrong boundary input'
               call mpi_finalize(ierr)
               stop
            endif
            exit
          endif
        enddo
        do ibd=1,nbd
          if(mrb(ibd).eq.2) then
            mbd=ibd
            if(any(mrb(mbd:nbd).ne.2)) then
               write(6,*)'wrong boundary input'
               call mpi_finalize(ierr)
               stop
            endif
            if(mbdn.gt.mbd)mbdn=mbd
            exit
          endif
        enddo

        clock(5) = tclock()-clock(5)

c
c.....allocate arrays
c
        clock(6) = tclock()
c
        if(idomy==0) then
          nyl = ny
          nzl = nz
          jys = 1
          kzs = (nz-2)*myrank+1
        else
          nyl = (ny-2)/mysize+2
          nzl = nzg
          jys = (ny-2)*myrank/mysize + 1
          kzs = 1
        endif
c
        nfacetmax = nfacet
c
        ALLOCATE(unvect(3,nfacet),vertex(3,3,nfacet),vertexc(3,nfacet))

        IF(ITBDY>0 .OR.ICALF/=0) THEN
           ALLOCATE(limpbd(nbd,nflp),mimpbd(nbd,nflp))
          IF(IVRTX==0) THEN
            ALLOCATE(pbd(nfacet),mrkpb(nfacet),mrksb(9,nfacet)    !!!!!! instead of mrksb(nfacet,9)
     &          ,dudxb(nfacet),dudyb(nfacet),dudzb(nfacet)
     &          ,dvdxb(nfacet),dvdyb(nfacet),dvdzb(nfacet)
     &          ,dwdxb(nfacet),dwdyb(nfacet),dwdzb(nfacet),cf_aux(nfacet,6))
          ELSEIF(IVRTX==1) THEN
            ALLOCATE(limpnd(nbd,nflp),mimpnd(nbd,nflp))
            IF(ICALF>0) THEN
              n = max(nfacet,nvtx)
            ELSE
              n = nvtx
            ENDIF
            ALLOCATE(trvtx(3,nvtx),vtxnrm(3,nvtx),vtxno(nvtx)) !vtx2tr(nvtx))
            ALLOCATE(pbd(n),mrkpb(n),mrksb(9,n))
            ALLOCATE(dudxb(n),dudyb(n),dudzb(n)
     &            ,dvdxb(n),dvdyb(n),dvdzb(n)
     &            ,dwdxb(n),dwdyb(n),dwdzb(n),cf_aux(n,6))
          ENDIF
        ENDIF
        ALLOCATE(areaf(nfacet),trino(nfacet))

        ALLOCATE(dsb(NBD,NX),mbg(nbd))

        ALLOCATE(limu(nbd,nflu),limv(nbd,nflu),limw(nbd,nflu),limp(nbd,nflp)
     &       ,   mimu(nbd,nflu),mimv(nbd,nflu),mimw(nbd,nflu),mimp(nbd,nflp))
        IF(ISGS>0) THEN
          ALLOCATE(limtv(nbd,nflu),mimtv(nbd,nflu))
        ENDIF

c
        ALLOCATE(flaguo(nx,nyl,nzl,nbd),flagvo(nx,nyl,nzl,nbd))
        ALLOCATE(flagwo(nx,nyl,nzl,nbd),flagpo(nx,nyl,nzl,nbd))
        ALLOCATE(flag(nx,nyl,nzl))
        IF(isgs>0) ALLOCATE(flagtv(nx,nyl,nzl,nbd))
        IF(itcalf>0 .OR. itbdy>0) ALLOCATE(flagpbd(nx,nyl,nzl,nbd))
c
        ALLOCATE(iu(nfcmax),ju(nfcmax),ku(nfcmax))
        ALLOCATE(iv(nfcmax),jv(nfcmax),kv(nfcmax))
        ALLOCATE(iw(nfcmax),jw(nfcmax),kw(nfcmax))
        ALLOCATE(ip(nfcmax),jp(nfcmax),kp(nfcmax))
        ALLOCATE(fp(nfcmax),fpu(nfcmax),fpv(nfcmax),fpw(nfcmax),fptv(nfcmax))

        IF(ISGS>0) THEN
          ALLOCATE(itv(nfcmax),jtv(nfcmax),ktv(nfcmax))
        ENDIF

        IF(idomy==1 .AND. mysize>1) THEN
          ALLOCATE(iuint(nintmax),juint(nintmax),kuint(nintmax))
          ALLOCATE(ivint(nintmax),jvint(nintmax),kvint(nintmax))
          ALLOCATE(iwint(nintmax),jwint(nintmax),kwint(nintmax))
          ALLOCATE(ipint(nintmax),jpint(nintmax),kpint(nintmax))

          ALLOCATE(limuint(nbd),limvint(nbd),limwint(nbd),limpint(nbd)
     &            ,mimuint(nbd),mimvint(nbd),mimwint(nbd),mimpint(nbd))
          ALLOCATE(itvint(nintmax),jtvint(nintmax),ktvint(nintmax))
          ALLOCATE(limtvint(nbd),mimtvint(nbd))
        ENDIF

        ALLOCATE(mrku(nfcmax),mrkv(nfcmax),mrkw(nfcmax),mrkp(nfcmax))
        ALLOCATE(diru(nfcmax),dirv(nfcmax),dirw(nfcmax),dirp(nfcmax))

        ALLOCATE(xnu(nfcmax),ynu(nfcmax),znu(nfcmax)
     &       ,   nxu(nfcmax),nyu(nfcmax),nzu(nfcmax),uim(nfcmax))
        ALLOCATE(xnv(nfcmax),ynv(nfcmax),znv(nfcmax)
     &       ,   nxv(nfcmax),nyv(nfcmax),nzv(nfcmax),vim(nfcmax))
        ALLOCATE(xnw(nfcmax),ynw(nfcmax),znw(nfcmax)
     &       ,   nxw(nfcmax),nyw(nfcmax),nzw(nfcmax),wim(nfcmax))
        ALLOCATE(xnp(nfcmax),ynp(nfcmax),znp(nfcmax),pim(nfcmax),rhoim(nfcmax))
        ALLOCATE(nxp(nfcmax),nyp(nfcmax),nzp(nfcmax))
c
        ALLOCATE(umtrx(nsp,nsp,nfcmax),vmtrx(nsp,nsp,nfcmax)
     &       ,   wmtrx(nsp,nsp,nfcmax),pmtrx(nsp,nsp,nfcmax))
        ALLOCATE(uindx(nsp,nfcmax),vindx(nsp,nfcmax)
     &       ,   windx(nsp,nfcmax),pindx(nsp,nfcmax))
c
        ALLOCATE(dpdnn(nfcmax),drhodnn(nfcmax))
        ALLOCATE(nflumax(nbd),nflvmax(nbd),nflwmax(nbd),nflpmax(nbd),nflpbdmax(nbd)
     &,nfltvmax(nbd))
        IF(ISGS>0) THEN
          ALLOCATE(mrktv(nfcmax),dirtv(nfcmax))
          ALLOCATE(xntv(nfcmax),yntv(nfcmax),zntv(nfcmax)
     &            ,nxtv(nfcmax),nytv(nfcmax),nztv(nfcmax),tvim(nfcmax))
          ALLOCATE(tvmtrx(nsp,nsp,nfcmax),tvindx(nsp,nfcmax))
        ENDIF
c
        clock(6) = tclock() - clock(6)
c.... parametric description of immersed object
c
        clock(7) = tclock()

        mbg(1:nbd) = mb(1:nbd)
c
c Read the STL files of the immersed-boundaries
c
!!!!!!        CALL READSTL(UNVECT,VERTEX,VERTEXC,AREAF,TRINO,NFACET,ZC,NZ,NBD,MBD,0)
        CALL READSTL_MOD(UNVECT,VERTEX,VERTEXC,AREAF,TRINO,NFACET,ZC,NZ,NBD,MBD,1)
        IF(ITBDY>0 .AND. IVRTX==1) CALL READTRVTX(TRVTX,VTXNRM,VTXNO,NVTX,NBD,MBD)   !!!!!! NBD,MBD

c        ibd = 1
c        CALL WRITESTL('mbd.'//index(ibd)//'.res',unvect,vertex,vertexc,zc,nz,nfacet,nfacetot,ibd,0)
c        CALL IOMPI_IMB(UNVECT,VERTEX,VERTEXC,AREAF,TRINO,NFACET,ZC,NZ,NBD,MBD,MBG,0)
c        CALL IOMPI_IMB(UNVECT,VERTEX,VERTEXC,AREAF,TRINO,NFACET,ZC,NZ,NBD,MBD,MBG,1)
        clock(7) = tclock() - clock(7)

        erb0 = 0.

        erb = 0.
        rrb = 0.
        arb = 0.
        iflag_psb = 0
        iflag_ptr = 0
        iflag_pnd = 0

c
c.... obstacle location on the grid
c
        NIMU = 0
        NIMV = 0
        NIMW = 0
        NIMP = 0
        NIMTV = 0

        MIMU = 0
        MIMV = 0
        MIMW = 0
        MIMP = 0
        IF(ISGS>0) MIMTV = 0

c        do icycle=1,itmax
c        tlevel = tlevel+tstep
c        write(6,*) 'icycle=',icycle,', tlevel=',tlevel,',tstep=',tstep
c        call rbm(unvect,vertex,vertexc,nfacet,mbd,nbd,tstep,tlevel)

        icycle = 0
        do ibd=1,nbd

          ilb = lb(ibd)+1
          ile = lb(ibd)+mb(ibd)

          clock(8) = tclock()
c
c evaluation of the max length of the triangles edges (DSB)
c
          call calc_ds_triangles(vertex,vertexc,nfacet,ilb,ile,dsb(ibd,:),xu,nx,icyl)
          clock(8) = tclock() - clock(8)

          clock(9) = tclock()
c
c definition of the limits of the subdomains enclosing the immersed
c boundaries along the X and the Z directions: IBMIN,IBMAX,JBMIN,JBMAX,
c KBMIN,KBMAX
c
!!!!!!          call limimb(xu,yv(jys),zwg(kzs),zcg(kzs),nx,nyl,nzl,vertexc,nfacet,ibd)
          call limimb
     %(xu,yv(jys),zwg(kzs),zcg(kzs),nx,nyl,nzl,vertexc(:,ilb:ile),mb(ibd),ibd)
          clock(9) = tclock() - clock(9)
c
c... Find forcing points and build stencil for U
c....FLAG is 1 for the exterior points, -1 for the interior ones
          clock(10) = tclock()
!!!!!!          call tag3d(xu_car(:,jys),yu_car(:,jys),zcg(kzs),nx,nyl,nzl
!!!!!!     &         ,vertex,vertexc,unvect,nfacet,ilb,ile,dsb(ibd,1),flag,ibd,0)
           call tag3d_mod
     %(xu_car(:,jys),yu_car(:,jys),zcg(kzs),nx,nyl,nzl,
     %vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile),
     %mb(ibd),dsb(ibd,:),flag,ibd,0)
          clock(10) = tclock() - clock(10)

          clock(11) = tclock()
c....Definition of the interior (flaguo=ibd), the exterior (flaguo=0) and
c....the interface points (flaguo=-ibd)
          CALL taguo(flag,flaguo,nx,nyl,nzl,nbd,ibd)
          clock(11) = tclock() - clock(11)
c
c... Find forcing points and build stencil for V
          clock(12) = tclock()
!!!!!!          call tag3d(xv_car(:,jys),yv_car(:,jys),zcg(kzs),nx,nyl,nzl
!!!!!!     &         ,vertex,vertexc,unvect,nfacet,ilb,ile,dsb(ibd,1),flag,ibd,1)
          call tag3d_mod
     %(xv_car(:,jys),yv_car(:,jys),zcg(kzs),nx,nyl,nzl,
     %vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile),
     %mb(ibd),dsb(ibd,:),flag,ibd,1)
          clock(12) = tclock() - clock(12)

          clock(13) = tclock()
          CALL taguo(flag,flagvo,nx,nyl,nzl,nbd,ibd)
          clock(13) = tclock() - clock(13)
c
c... Find forcing points and build stencil for W
          clock(14) = tclock()
          icom = 2
          IF(ISGS>0) icom=1
!!!!!!          call tag3d(xc_car(:,jys),yc_car(:,jys),zwg(kzs),nx,nyl,nzl
!!!!!!     &         ,vertex,vertexc,unvect,nfacet,ilb,ile,dsb(ibd,1),flag,ibd,icom)
          call tag3d_mod
     %(xc_car(:,jys),yc_car(:,jys),zwg(kzs),nx,nyl,nzl,
     %vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile),
     %mb(ibd),dsb(ibd,:),flag,ibd,icom)
          clock(14) = tclock() - clock(14)

          clock(15) = tclock()
          CALL taguo(flag,flagwo,nx,nyl,nzl,nbd,ibd)
          clock(15) = tclock() - clock(15)
c
c... Find forcing points and build stencil for TV
          IF(ISGS>0) THEN
            clock(16) = tclock()
!!!!!!            call tag3d(xc_car(:,jys),yc_car(:,jys),zcg(kzs),nx,nyl,nzl
!!!!!!     &         ,vertex,vertexc,unvect,nfacet,ilb,ile,dsb(ibd,1),flag,ibd,2)
            call tag3d_mod
     %(xc_car(:,jys),yc_car(:,jys),zcg(kzs),nx,nyl,nzl,
     %vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile),
     %mb(ibd),dsb(ibd,:),flag,ibd,2)
            clock(16) = tclock() - clock(16)

            clock(15) = tclock()
            CALL taguo(flag,flagtv,nx,nyl,nzl,nbd,ibd)
            clock(15) = tclock() - clock(15)
          ENDIF

!!!!!!          nimu = 0
!!!!!!          nimv = 0
!!!!!!          nimw = 0
!!!!!!          nimtv = 0
          icom = 0

          clock(20) = tclock()
c....Interpolation stencils for each interface point
!!!!!!          call geom
!!!!!!     &(xu,yc(jys),xu_car(:,jys),yu_car(:,jys),zcg(kzs),vertex,vertexc,unvect
!!!!!!     &,limu(ibd,:),mimu(ibd,:),iu,ju,ku,nimu,xnu,ynu,znu,nxu,nyu,nzu
!!!!!!     &,mrku,diru,fp,flag,flaguo,nx,nyl,nzl,nbd,nfacet,ibd,icom,nflu,nflumax,1,icycle,tlevel)
ccc          call geom
          call geom_mod
     %(xu,yc(jys),xu_car(:,jys),yu_car(:,jys),zcg(kzs),
     %vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile),
     %limu(ibd,:),mimu(ibd,:),iu,ju,ku,nimu,xnu,ynu,znu,
     %nxu,nyu,nzu,mrku,diru,fpu,flag,flaguo,nx,nyl,nzl,nbd,
     %mb(ibd),ibd,icom,nflu,nflumax(ibd),1,icycle,tlevel)
          clock(20) = tclock() - clock(20)

          clock(26) = tclock()
c....Definition of the matrix of coefficients for the outer velocity points
          call mtrx(xu,yc(jys),zcg(kzs),limu(ibd,1),sum(mimu(ibd,:)),mrku
     &         ,iu,ju,ku,xnu,ynu,znu,nxu,nyu,nzu,umtrx,uindx,nx,nyl,nzl,nbd,ibd,1)
          clock(26) = tclock() - clock(26)

          icom = 1
          clock(22) = tclock()
!!!!!!          call geom
!!!!!!     &(xc,yv(jys),xv_car(:,jys),yv_car(:,jys),zcg(kzs),vertex,vertexc,unvect
!!!!!!     &,limv(ibd,1),mimv(ibd,1),iv,jv,kv,nimv,xnv,ynv,znv,nxv,nyv,nzv
!!!!!!     &,mrkv,dirv,fp,flag,flagvo,nx,nyl,nzl,nbd,nfacet,ibd,icom,nflu,nflvmax,2,icycle,tlevel)
ccc          call geom
          call geom_mod
     %(xc,yv(jys),xv_car(:,jys),yv_car(:,jys),zcg(kzs),
     %vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile),
     %limv(ibd,:),mimv(ibd,:),iv,jv,kv,nimv,xnv,ynv,znv,
     %nxv,nyv,nzv,mrkv,dirv,fpv,flag,flagvo,nx,nyl,nzl,nbd,
     %mb(ibd),ibd,icom,nflu,nflvmax(ibd),2,icycle,tlevel)
          clock(22) = tclock() - clock(22)

          clock(28) = tclock()
          call mtrx(xc,yv(jys),zcg(kzs),limv(ibd,1),sum(mimv(ibd,:)),mrkv
     &         ,iv,jv,kv,xnv,ynv,znv,nxv,nyv,nzv,vmtrx,vindx,nx,nyl,nzl,nbd,ibd,1)
          clock(28) = tclock() - clock(28)

          clock(24) = tclock()
!!!!!!          call geom
!!!!!!     &(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zwg(kzs),vertex,vertexc,unvect
!!!!!!     &,limw(ibd,1),mimw(ibd,1),iw,jw,kw,nimw,xnw,ynw,znw,nxw,nyw,nzw
!!!!!!     &,mrkw,dirw,fp,flag,flagwo,nx,nyl,nzl,nbd,nfacet,ibd,icom,nflu,nflwmax,3,icycle,tlevel)
ccc          call geom
          call geom_mod
     %(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zwg(kzs),
     %vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile),
     %limw(ibd,:),mimw(ibd,:),iw,jw,kw,nimw,xnw,ynw,znw,
     %nxw,nyw,nzw,mrkw,dirw,fpw,flag,flagwo,nx,nyl,nzl,nbd,
     %mb(ibd),ibd,icom,nflu,nflwmax(ibd),3,icycle,tlevel)
          clock(24) = tclock() - clock(24)
c
          clock(30) = tclock()
          call mtrx(xc,yc(jys),zwg(kzs),limw(ibd,1),sum(mimw(ibd,:)),mrkw
     &         ,iw,jw,kw,xnw,ynw,znw,nxw,nyw,nzw,wmtrx,windx,nx,nyl,nzl,nbd,ibd,1)
          clock(30) = tclock() - clock(30)

          IF(ISGS>0) THEN
            clocktemp = tclock()
!!!!!!            call
!!!!!!     &geom(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs),vertex,vertexc,unvect
!!!!!!     &,limtv(ibd,:),mimtv(ibd,:),itv,jtv,ktv,nimtv,xntv,yntv,zntv,nxtv,nytv,nztv
!!!!!!     &,mrktv,dirtv,fp,flag,flagtv,nx,nyl,nzl,nbd,nfacet,ibd,icom,nflu,nfltvmax,5,icycle,tlevel)
            call
ccc     &geom(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs),
     &geom_mod(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs),
     &vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile),
     &limtv(ibd,:),mimtv(ibd,:),itv,jtv,ktv,nimtv,xntv,yntv,zntv,nxtv,nytv,nztv,
     &mrktv,dirtv,fptv,flag,flagtv,nx,nyl,nzl,nbd,mb(ibd),ibd,icom,nflu,nfltvmax(ibd),5,icycle,tlevel)
            clock(24) = clock(24) + tclock() - clocktemp

            clocktemp = tclock()
            call mtrx(xc,yc(jys),zcg(kzs),limtv(ibd,1),sum(mimtv(ibd,:)),mrktv
     &,itv,jtv,ktv,xntv,yntv,zntv,nxtv,nytv,nztv,tvmtrx,tvindx,nx,nyl,nzl,nbd,ibd,1)
            clock(30) = clock(30) + tclock() - clocktemp
          ENDIF

          clock(17) = tclock()
c....Definition of the pressure exterior (flagpo=0),
c....interface (flagpo=-ibd) and interior (flagpo=ibd) points
          call flagp(flagpo,flaguo,flagvo,flagwo,nx,nyl,nzl,nbd,ibd)
c          call flagp1(flagpo,flaguo,flagvo,flagwo,nx,nyl,nzl,nbd,iu,ju,ku,iv,jv,kv,iw,jw,kw
c     &         ,limu(ibd,1),sum(mimu(ibd,:)),limv(ibd,1),sum(mimv(ibd,:)),limw(ibd,1),sum(mimw(ibd,:)),ibd)
          clock(17) = tclock() - clock(17)

          clock(32) = tclock()
!!!!!!          nimp = 0
          icom = 2
!!!!!!          IF(IBM==1 .AND. (ITBDY>0 .OR. ICALF>0)) icom=1 ! a flag for COMMITTEE3
!!!!!!          call geom
!!!!!!     &(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs),vertex,vertexc,unvect
!!!!!!     &,limp(ibd,1),mimp(ibd,1),ip,jp,kp,nimp,xnp,ynp,znp,nxp,nyp,nzp
!!!!!!     &,mrkp,dirp,fp,flag,flagpo,nx,nyl,nzl,nbd,nfacet,ibd,icom,nflp,nflpmax,4,icycle,tlevel)
ccc          call geom
          call geom_mod
     %(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs),
     %vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile),
     %limp(ibd,:),mimp(ibd,:),ip,jp,kp,nimp,xnp,ynp,znp,
     %nxp,nyp,nzp,mrkp,dirp,fp,flag,flagpo,nx,nyl,nzl,nbd,
     %mb(ibd),ibd,icom,nflp,nflpmax(ibd),4,icycle,tlevel)
          clock(32) = tclock() - clock(32)

          clock(34) = tclock()
          call mtrx(xc,yc(jys),zcg(kzs),limp(ibd,1),sum(mimp(ibd,:)),mrkp
     &         ,ip,jp,kp,xnp,ynp,znp,nxp,nyp,nzp,pmtrx,pindx,nx,nyl,nzl,nbd,ibd,0)
          clock(34) = tclock() - clock(34)

        enddo

c        stop
c in the case the matrices of the tagging are not large enough
        i=sum(mimu)
        if(i>nfcmax) then
          write(6,'(A,1X,I10)') 'ERROR:mimu increase nfcmax to at least',i,mimu
          call mpi_finalize(ierr)
          stop
        endif

        i=sum(mimv)
        if(i>nfcmax) then
          write(6,'(A,1X,I10)') 'ERROR:mimv increase nfcmax to at least',i
          call mpi_finalize(ierr)
          stop
        endif

        i=sum(mimw)
        if(i>nfcmax) then
          write(6,'(A,1X,I10)') 'ERROR:mimw increase nfcmax to at least',i
          call mpi_finalize(ierr)
          stop
        endif

        i=sum(mimp)
        if(i>nfcmax) then
          write(6,'(A,1X,I10)') 'ERROR:mimp increase nfcmax to at least',i
          call mpi_finalize(ierr)
          stop
        endif

        clock(35) = tclock()
c
c...Set surface boundary velocity for non-moving bodies
c...A modification allows to take into account the moving bodies
c...whose geometry is constant
c...MBDN is the index of the first moving body
c...MBD is the index of the first variable moving body
!!!!!!        do ibd=1,mbd-1
        do ibd=1,mbdn-1
          do im=limu(ibd,1)+1,limu(ibd,1)+sum(mimu(ibd,:))
            uim(im) = 0.0
          enddo
          do jm=limv(ibd,1)+1,limv(ibd,1)+sum(mimv(ibd,:))
            vim(jm) = 0.0
          enddo
          do km=limw(ibd,1)+1,limw(ibd,1)+sum(mimw(ibd,:))
            wim(km) = 0.0
          enddo
          do i=limp(ibd,1)+1,limp(ibd,1)+sum(mimp(ibd,:))
            dpdnn(i) = 0.0
          enddo
          do i=limp(ibd,1)+1,limp(ibd,1)+sum(mimp(ibd,:))
            drhodnn(i) = 0.0
          enddo
          IF(ISGS>0) THEN
            do im=limtv(ibd,1)+1,limtv(ibd,1)+sum(mimtv(ibd,:))
              tvim(im) = 0.0
            enddo
          ENDIF
        enddo
c
c...Set surface boundary velocity for moving bodies
c...Input: Cartesian coordinates of the immersed-boundary points
c...Output: velocities in Cartesian or cylindrical coordinates
c...The angular velocity is AMP in EDDY.INPUT
c...A modification allows to take into account the moving bodies
c...whose geometry is constant
!!!!!!        do ibd=mbd,nbd
        do ibd=mbdn,nbd
          do im=limu(ibd,1)+1,limu(ibd,1)+sum(mimu(ibd,:))
            uim(im) = ubd(xnu(im),ynu(im),znu(im),tlevel,ibd)
          enddo
          do jm=limv(ibd,1)+1,limv(ibd,1)+sum(mimv(ibd,:))
            vim(jm) = vbd(xnv(jm),ynv(jm),znv(jm),tlevel,ibd)
          enddo
          do km=limw(ibd,1)+1,limw(ibd,1)+sum(mimw(ibd,:))
            wim(km) = wbd(xnw(km),ynw(km),znw(km),tlevel,ibd)
          enddo
          do i=limp(ibd,1)+1,limp(ibd,1)+sum(mimp(ibd,:))
            dpdnn(i) = 0.0
          enddo
          do i=limp(ibd,1)+1,limp(ibd,1)+sum(mimp(ibd,:))
            drhodnn(i) = 0.0
          enddo
          IF(ISGS>0) THEN
            do im=limtv(ibd,1)+1,limtv(ibd,1)+sum(mimtv(ibd,:))
              tvim(im) = 0.0
            enddo
          ENDIF
        enddo

        clock(35) = tclock() - clock(35)

        IF(idomy==1 .AND. mysize>1) THEN
c
c in the case there is a decomposition along y for tagging purposes: the
c indices of the interior points are defined for each variable
c
          mimuint = 0
          mimvint = 0
          mimwint = 0
          mimpint = 0
          DO ibd=1,nbd
            limuint(ibd) = sum(mimuint(1:ibd-1))
            call flag2int(flaguo(:,:,:,ibd),nx,nyl,nzl,ibd,iuint,juint,kuint,limuint(ibd),mimuint(ibd))
            limvint(ibd) = sum(mimvint(1:ibd-1))
            call flag2int(flagvo(:,:,:,ibd),nx,nyl,nzl,ibd,ivint,jvint,kvint,limvint(ibd),mimvint(ibd))
            limwint(ibd) = sum(mimwint(1:ibd-1))
            call flag2int(flagwo(:,:,:,ibd),nx,nyl,nzl,ibd,iwint,jwint,kwint,limwint(ibd),mimwint(ibd))
            limpint(ibd) = sum(mimpint(1:ibd-1))
            call flag2int(flagpo(:,:,:,ibd),nx,nyl,nzl,ibd,ipint,jpint,kpint,limpint(ibd),mimpint(ibd))
          ENDDO
        ENDIF

        if(idomy==1 .AND. mysize>1) then
c
c Conversion of the tagging variables from the decomposition along Y to the
c one along Z
c
!!!!!!          call imb_domy2z(iu,ju,ku,xnu,ynu,znu,nxu,nyu,nzu,mrku,diru,fp
!!!!!!     &         ,uim,umtrx,uindx,zcg,nzg,zc,nz,nyl,limu(1,1),mimu(1,1),nbd,nflu,0)
          call imb_domy2z(iu,ju,ku,xnu,ynu,znu,nxu,nyu,nzu,mrku,diru,fpu
     &         ,uim,umtrx,uindx,zcg,nzg,zc,nz,nyl,limu(1,1),mimu(1,1),nbd,nflu,0)
!!!!!!          call imb_domy2z(iv,jv,kv,xnv,ynv,znv,nxv,nyv,nzv,mrkv,dirv,fp
!!!!!!     &         ,vim,vmtrx,vindx,zcg,nzg,zc,nz,nyl,limv(1,1),mimv(1,1),nbd,nflu,0)
          call imb_domy2z(iv,jv,kv,xnv,ynv,znv,nxv,nyv,nzv,mrkv,dirv,fpv
     &         ,vim,vmtrx,vindx,zcg,nzg,zc,nz,nyl,limv(1,1),mimv(1,1),nbd,nflu,0)
!!!!!!          call imb_domy2z(iw,jw,kw,xnw,ynw,znw,nxw,nyw,nzw,mrkw,dirw,fp
!!!!!!     &         ,wim,wmtrx,windx,zcg,nzg,zc,nz,nyl,limw(1,1),mimw(1,1),nbd,nflu,0)
          call imb_domy2z(iw,jw,kw,xnw,ynw,znw,nxw,nyw,nzw,mrkw,dirw,fpw
     &         ,wim,wmtrx,windx,zcg,nzg,zc,nz,nyl,limw(1,1),mimw(1,1),nbd,nflu,0)
          call imb_domy2z(ip,jp,kp,xnp,ynp,znp,nxp,nyp,nzp,mrkp,dirp,fp
     &         ,pim,pmtrx,pindx,zcg,nzg,zc,nz,nyl,limp(1,1),mimp(1,1),nbd,nflp,1)

          call intind_domy2z(iuint,juint,kuint,limuint,mimuint,nbd,nyl,nz)
          call intind_domy2z(ivint,jvint,kvint,limvint,mimvint,nbd,nyl,nz)
          call intind_domy2z(iwint,jwint,kwint,limwint,mimwint,nbd,nyl,nz)
          call intind_domy2z(ipint,jpint,kpint,limpint,mimpint,nbd,nyl,nz)
        endif

        clocktemp = tclock()
c        do ibd=1,nbd
c....A PLT TECPLOT file is written for each immersed-boundary
c....It could be useful to modify this subroutine to write an STL file !!!!!!
c          call imb2plt('mbd.'//index(ibd)//'.plt',vertex,nfacet,ibd)
c
c          WRITE(6,*) 'Calling IOSCALAR2'
c          CALL IOSCALAR6('flaguo.'//INDEX(IBD),FLAGUO(1,1,1,IBD),NX,NY,NZ,1)
c          CALL IOSCALAR2(trim(hspdir)//'flaguo.'//INDEX(IBD),FLAGUO(1,1,1,IBD),NX,NY,NZ,1)
c          CALL IOSCALAR2(trim(hspdir)//'flagvo.'//INDEX(IBD),FLAGVO(1,1,1,IBD),NX,NY,NZ,1)
c          CALL IOSCALAR2(trim(hspdir)//'flagwo.'//INDEX(IBD),FLAGWO(1,1,1,IBD),NX,NY,NZ,1)
c          CALL IOSCALAR2(trim(hspdir)//'flagp.'//INDEX(IBD),FLAGPO(1,1,1,IBD),NX,NY,NZ,1)
c          flag = 0
c          where(flaguo(:,:,:,ibd)==-ibd) flag=1
c          CALL IOSCALAR2('flag_ufrc.'//INDEX(IBD),FLAG,NX,NY,NZ,1)
c          flag = 0
c          where(flaguo(:,:,:,ibd)==0) flag=1
c          CALL IOSCALAR2('flag_ufld.'//INDEX(IBD),FLAG,NX,NY,NZ,1)
c          flag = 0
c          where(flaguo(:,:,:,ibd)==ibd) flag=1
c          CALL IOSCALAR2('flag_ubdy.'//INDEX(IBD),FLAG,NX,NY,NZ,1)
c
c          flag = 0
c          where(flagpo(:,:,:,ibd)==-ibd) flag=1
c          CALL IOSCALAR2('flag_pfrc.'//INDEX(IBD),FLAG,NX,NY,NZ,1)
c          flag = 0
c          where(flagpo(:,:,:,ibd)==0) flag=1
c          CALL IOSCALAR2('flag_pfld.'//INDEX(IBD),FLAG,NX,NY,NZ,1)
c          flag = 0
c          where(flagpo(:,:,:,ibd)==ibd) flag=1
c          CALL IOSCALAR2('flag_pbdy.'//INDEX(IBD),FLAG,NX,NY,NZ,1)
c
c          CALL IOSCALAR2('flagui.'//INDEX(IBD),FLAGUI(1,1,1,IBD),NX,NY,NZ,1)
c
c          CALL IOMFORC1('mforc.uo.'//INDEX(IBD),IU,JU,KU,MRKU
c     &      ,IUMTRX,JUMTRX,KUMTRX,UINDX,UMTRX,UIM,XNU,YNU,ZNU,NXU,NYU,NZU,LIMU(IBD,1),MIMU(IBD,1),NZ,1)
c          CALL IOMFORC1('mforc.ui.'//INDEX(IBD),IU,JU,KU,MRKU
c     &      ,IUMTRX,JUMTRX,KUMTRX,UINDX,UMTRX,UIM,XNU,YNU,ZNU,NXU,NYU,NZU,LIMU(IBD,2),MIMU(IBD,2),NZ,1)
c          CALL IOMFORC2('mforc2.uo.'//INDEX(IBD)//'.plt',XU_CAR,YU_CAR
c     &         ,ZC,NX,NY,NZ,IU,JU,KU,MRKU,XNU,YNU,ZNU,NXU,NYU,NZU,LIMU(IBD,1),SUM(MIMU(IBD,:)),1)
c
c          CALL IOSCALAR6('flagvo.'//INDEX(IBD),FLAGVO(1,1,1,IBD),NX,NY,NZ,1)
c          CALL IOSCALAR2('flagvo.'//INDEX(IBD),FLAGVO(1,1,1,IBD),NX,NY,NZ,1)
c          CALL IOMFORC1('mforc.vo.'//INDEX(IBD),IV,JV,KV,MRKV
c     &      ,IVMTRX,JVMTRX,KVMTRX,VINDX,VMTRX,VIM,XNV,YNV,ZNV,NXV,NYV,NZV,LIMV(IBD,1),MIMV(IBD,1),NZ,1)
c          CALL IOMFORC1('mforc.vi.'//INDEX(IBD),IV,JV,KV,MRKV
c     &      ,IVMTRX,JVMTRX,KVMTRX,VINDX,VMTRX,VIM,XNV,YNV,ZNV,NXV,NYV,NZV,LIMV(IBD,2),MIMV(IBD,2),NZ,1)
c          CALL IOMFORC2('mforc2.vo.'//INDEX(IBD)//'.plt',XV_CAR,YV_CAR
c     &         ,ZC,NX,NY,NZ,IV,JV,KV,MRKV,XNV,YNV,ZNV,NXV,NYV,NZV,LIMV(IBD,1),SUM(MIMV(IBD,:)),1)
c
c          CALL IOSCALAR6('flagwo.'//INDEX(IBD),FLAGWO(1,1,1,IBD),NX,NY,NZ,1,0.0)
c          CALL IOSCALAR2('flagwo.'//INDEX(IBD),FLAGWO(1,1,1,IBD),NX,NY,NZ,1)
c          CALL IOMFORC1('mforc.wo.'//INDEX(IBD),IW,JW,KW,MRKW
c     &      ,IWMTRX,JWMTRX,KWMTRX,WINDX,WMTRX,WIM,XNW,YNW,ZNW,NXW,NYW,NZW,LIMW(IBD,1),MIMW(IBD,1),NZ,1)
c          CALL IOMFORC1('mforc.wi.'//INDEX(IBD),IW,JW,KW,MRKW
c     &      ,IWMTRX,JWMTRX,KWMTRX,WINDX,WMTRX,WIM,XNW,YNW,ZNW,NXW,NYW,NZW,LIMW(IBD,2),MIMW(IBD,2),NZ,1)
c          CALL IOMFORC2('mforc2.wo.'//INDEX(IBD)//'.plt',XC_CAR,YC_CAR
c     &         ,ZW,NX,NY,NZ,IW,JW,KW,MRKW,XNW,YNW,ZNW,NXW,NYW,NZW,LIMW(IBD,1),SUM(MIMW(IBD,:)),1)
c
c          CALL IOSCALAR2('flagpo.'//INDEX(IBD),FLAGPO(1,1,1,IBD),NX,NY,NZ,1)
c          CALL IOSCALAR6('flagpo.'//INDEX(IBD),FLAGPO(1,1,1,IBD),NX,NY,NZ,1,0.0)
c          call mpi_finalize(ierr)
c          stop
c          CALL IOMFORC1('mforc.po.'//INDEX(IBD),IP,JP,KP,MRKP
c     &      ,IPMTRX,JPMTRX,KPMTRX,PINDX,PMTRX,DPDNN,XNP,YNP,ZNP,NXP,NYP,NZP,LIMP(IBD,1),MIMP(IBD,1),NZ,1)
c          CALL IOMFORC1('mforc.pi.'//INDEX(IBD),IP,JP,KP,MRKP
c     &      ,IPMTRX,JPMTRX,KPMTRX,PINDX,PMTRX,DPDNN,XNP,YNP,ZNP,NXP,NYP,NZP,LIMP(IBD,2),MIMP(IBD,2),NZ,1)
c          CALL IOMFORC2('mforc2.po.'//INDEX(IBD)//'.plt',XC_CAR,YC_CAR
c     &         ,ZC,NX,NY,NZ,IP,JP,KP,MRKP,XNP,YNP,ZNP,NXP,NYP,NZP,LIMP(IBD,1),SUM(MIMP(IBD,:)),1)
c          CALL IOMFORC2('mforc2.pi.'//INDEX(IBD)//'.plt',XC_CAR,YC_CAR
c     &         ,ZC,NX,NY,NZ,IP,JP,KP,IPMTRX,JPMTRX,KPMTRX,MRKP
c     &         ,XNP,YNP,ZNP,NXP(:,2),NYP(:,2),NZP(:,2),LIMP(IBD,2),MIMP(IBD,2),1)
c          CALL IOSCALAR2('flagpbdo.'//INDEX(IBD),FLAGPBDO(1,1,1,IBD),NX,NY,NZ,1)
c          CALL IOSCALAR2('flagpbdi.'//INDEX(IBD),FLAGPBDI(1,1,1,IBD),NX,NY,NZ,1)
c
c          CALL IOMFORC2('mforc2.pbdo.'//INDEX(IBD)//'.plt',XC_CAR,YC_CAR
c     &         ,ZC,NX,NY,NZ,IP,JP,KP,IPMTRX,JPMTRX,KPMTRX,MRKP
c     &         ,XNP,YNP,ZNP,NXP(:,2),NYP(:,2),NZP(:,2),LIMPBD(IBD,1),MIMPBD(IBD,1),1)
c          CALL IOMFORC2('mforc2.pbdi.'//INDEX(IBD)//'.plt',XC_CAR,YC_CAR
c     &         ,ZC,NX,NY,NZ,IP,JP,KP,IPMTRX,JPMTRX,KPMTRX,MRKP
c     &         ,XNP,YNP,ZNP,NXP(:,2),NYP(:,2),NZP(:,2),LIMPBD(IBD,2),MIMPBD(IBD,2),1)
c
c        enddo

c        call mpi_finalize(ierr)
c        stop

c
        clock(80) = clock(80) + tclock() - clocktemp

      endif
c
c.....open file to write force information
c
      IF(IBM/=0.AND.ICALF>0) THEN
        ibd=1
        cvlim1(1)=-1.0
        cvlim1(2)= 1.0
        cvlim1(3)=-1.0
        cvlim1(4)= 1.0

c        call calcfrc(uo,vo,wo,p,flaguo(1,1,1,ibd),flagvo(1,1,1,ibd),flagwo(1,1,1,ibd)
c     &       ,nx,ny,nz,xu,yv,zw,xc,yc,zc,fb,tstep,0,time,'force',cvlim1)
      ENDIF
c
      IF(MYRANK==0) THEN
        write(6,'(A,E18.6)') '*..starting time=',TIME
        write(6,'(A)') '*..Initial field: '
      ENDIF

      open(49,file='time.out',position='append',form='formatted')

!      open(121,file='time_averaged_values_vmod_vr_vt_vz_pr_prtot_tv.txt',
!     %position='append')
!      open(122,file='time_averaged_values_vorr_voraz_vorz.txt',
!     %position='append')
!      open(123,file='time_averaged_values_2_vmod_vr_vt_vz_pr_prtot_tv.txt',
!     %position='append')
!      open(124,file='time_averaged_values_2_vorr_voraz_vorz.txt',
!     %position='append')
!      open(125,file='time_averaged_values_3_vmod_vr_vt_vz_pr_prtot_tv.txt',
!     %position='append')
!      open(126,file='time_averaged_values_3_vorr_voraz_vorz.txt',
!     %position='append')
!      if(irms.eq.1) then
!         open(127,file='RMS_vr_vt_vz_pr.txt',position='append')
!         open(128,file='time_averaged_values_uv_uw_vw.txt',
!     %position='append')
!         open(129,file='RMS_2_vr_vt_vz_pr.txt',position='append')
!         open(130,file='time_averaged_values_2_uv_uw_vw.txt',
!     %position='append')
!         open(131,file='RMS_3_vr_vt_vz_pr.txt',position='append')
!         open(132,file='time_averaged_values_3_uv_uw_vw.txt',
!     %position='append')
!      endif
! 334  format(a145)
!      do i=1,nprbmax2
!         ii=prbindx2+(i-1)
!         iid=ii/10
!         iiu=mod(ii,10)
!         open(700+ii,
!     %file='probes_'//char(48+iid)//char(48+iiu)//'.out',
!     %position='append')
!         open(800+ii,
!     %file='probes_'//char(48+iid)//char(48+iiu)//'_stat.out',
!     %position='append')
!      enddo

c      STOP
c
c-----------------------------------------------------------------------
c                                                   Check Starting Field
c-----------------------------------------------------------------------
      ICYCLE = nstep
c-----------------------------------------------------------------------
c                                            compute turbulent viscosity
c-----------------------------------------------------------------------
c

      IF(ISGS/=0) THEN

        clock(37) = tclock()
        IF(ISGS==5) THEN
!!!!!!          IF(IFIELD==0) CALL STRUCTUREFUNCTION(UO,VO,WO,TV,DP,NX,NY,NZ)
          IF(IFIELD>=0) CALL STRUCTUREFUNCTION(UO,VO,WO,TV,DP,NX,NY,NZ)
        ELSEIF(ISGS==6) THEN
!!!!!!          IF(IFIELD==0) CALL WALE(UO,VO,WO,TV,NX,NY,NZ)
          IF(IFIELD>=0) CALL WALE(UO,VO,WO,TV,NX,NY,NZ)
        ELSE
c....Evaluation of the filter size and the coefficients of the test
c....filter for each direction (WHX,WHY,WHZ)
          CALL filtersize(dxdydz,nx,ny,nz)
          CALL weights
c
          IF(IFIELD>=0) THEN
c....Evaluation of the eddy viscosity
            CALL TURVIS(UO,VO,WO,TV,G,DP,SXX,SYY,SZZ,SXY,SYZ,SXZ,
     &         US,VS,WS,UB,VB,WB,ILM,IMM,LM,MM,UC,VC,WC,
     &         DXDYDZ,CLES,CLESP,CLESN,NILMP,NILMN,XC,YC,ZC,TSTEP,ICYCLE,NX,NY,NZ)
          ENDIF
        ENDIF
c
        if(ibm/=0) then

ccc          do ibd=1,nbd
!!!!!!            call momforc(tv,tv,xc,yc,zc,nx,ny,nz,itv,jtv,ktv,tvmtrx,tvindx
!!!!!!     &,tvim,mrktv,dirtv,nxtv,nytv,nztv,nfltvmax,limtv(ibd,1),mimtv(ibd,1),0,3)
ccc            call momforc(tv,tv,xc,yc,zc,nx,ny,nz,itv,jtv,ktv,tvmtrx,tvindx
ccc     &,tvim,mrktv,dirtv,nxtv,nytv,nztv,nfltvmax(ibd),limtv(ibd,:),mimtv(ibd,:),0,3)
          call momforc_mod2(tv,tv,xc,yc,zc,nx,ny,nz,itv,jtv,ktv,tvmtrx,tvindx
     &,tvim,mrktv,dirtv,nxtv,nytv,nztv,nfltvmax,limtv,mimtv,0,3,nbd,sum(mimtv))

            tv(:,1,:) = tv(:,ny-1,:)
            tv(:,ny,:) = tv(:,2,:)
c
c......Set turbulent viscosity inside bodies
          do ibd=1,nbd
            IF(idomy==1 .AND. mysize>1) THEN
              call tvinterior1(tv,itv,jtv,ktv,xc,yc,zc,nx,ny,nz,limtvint(ibd),mimtvint(ibd))
            ELSE
              call tvinterior(tv,flagtv,nx,ny,nz,nbd,ibd)
            ENDIF
          enddo

        endif

        clock(37) = tclock() - clock(37)

      ENDIF
c....The minimum and maximum values of the velocities, pressure,
c....divergence and eddy viscosity are found and written on screen
      CALL SCRINFO(UO,VO,WO,P,DP,TV,DENS,TSTEP,TIME,ICYCLE,ID,JD,KD,NX,NY,NZ)
c
c.... Time-average top velocity (use to calc. Kappa acceleration parameter)
      flag_wavg_top = 0
      INQUIRE(FILE='wavg_top.input', EXIST=file_exists)
      IF(file_exists .eqv. .true.) THEN
        call read_wavg_top_prms(flag_wavg_top,freq_wavg_top)
        if(flag_wavg_top>0) then
          allocate(wavg_top(ny,nz))
          file_wavg_top = trim(hspdir)//'wavg_top.bin'
          wavg_top_nsmpl = 0
          call write_header_wavg_top(file_wavg_top,flag_wavg_top,wavg_top,ny,nz)
        endif
      ENDIF
c
c.... Read cf location for shear stress
      icfprb = 0
      INQUIRE(FILE='cf.input', EXIST=file_exists)
      IF(file_exists .eqv. .true.) THEN
        call read_no_cf_probes(njcf,nkcf,icfprb,itcfprb)
        ncfprb = itres/itcfprb
        ncfsmpl = 0
        if(icfprb>0) then
          allocate(jcf(njcf),kcf(nkcf))
          allocate(cfprb(njcf,nkcf,ncfprb),tcf(ncfprb))
          filecfprb = trim(hspdir)//'cfprb.bin'
          call setup_cf_probes(filecfprb,jcf,kcf,njcf,nkcf,nkcfmax,nkcfprev,xc,yc,zw,nx,ny,nz,icfprb)
        endif
      endif
c
c...  Read time probe location
      if(itmprb>=1) then

c.....Number of time probes
        call read_no_tprobes(ntprb)
        allocate(itprb(ntprb),jtprb(ntprb),ktprb(ntprb))
        allocate(utprb(ntprb,itres),vtprb(ntprb,itres),wtprb(ntprb,itres),ptprb(ntprb,itres))
        allocate(xprb(ntprb),yprb(ntprb),zprb(ntprb))

        iprbser = 1
        filetprb = trim(hspdir)//'probes.bin'
c.....Coordinates of the time probes (Cartesian or cylindrical)
c.....The probes are ordered according to increasing axial positions
c.....The restart files are eventually read
        call setup_tprobes(filetprb,itprb,jtprb,ktprb,ntprb,xc,yc,zc,nx,ny,nz,ntprbmax,tmprbindx,4,ntprbsmpl,time)
      endif
c
c...  Read probe locations for the azimuthal spectra
      if(yspctr>=1) then
        nycorvar = 4+4
        nyesdvar = 4
c.....Number of probes
        call read_no_yprobes(nyprb)
        allocate(iyprb(nyprb),jyprb(nyprb),kyprb(nyprb))
        allocate(esdprb(ny/2,nyprb,nyesdvar),ycorprb(ny-2,nyprb,nycorvar))
        allocate(esd(ny/2,nyprb,nyesdvar),ycor(ny-2,nyprb,nycorvar))
        esd = 0.0
        ycor = 0.0
        esdprb = 0.0
        ycorprb = 0.0
        fileycor = trim(hspdir)//'ycor.bin'
        fileyspc = trim(hspdir)//'yesd.bin'

        iyprb = 0
        kyprb = 0
        yprbindx = 0
        itysamp = 0
        yesdsmpl = 0
        ycorsmpl = 0
        nyprbmax = 0
c.....Locations of the probes for the evaluation of the correlations and
c.....the spectra along the Y direction
        call setup_yprobes(fileycor,iyprb,kyprb,nyprb,xc,zc,nx,nz
     &     ,nyprbmax,yprbindx,itysamp,nycorvar,ycorsmpl,ny-2)

        call setup_yprobes(fileyspc,iyprb,kyprb,nyprb,xc,zc,nx,nz
     &     ,nyprbmax,yprbindx,itysamp,nyesdvar,yesdsmpl,ny/2)
c.....The restart files of the correlations and the spectra are read
        if(yspctr>1) then
          call read_ycorrelations(fileycor,ycorprb,ny-2,nyprb,nycorvar,nyprbmax,yprbindx)
          call read_yspectra(fileyspc,esdprb,ny/2,nyprb,nyesdvar,nyprbmax,yprbindx)
        endif
      endif

      if(iregtmavg>0) then
        allocate(iindtmavg(nz,2))
c.....The indices limiting a particular portion of the computational domain
c.....are established
        call region2indices(kmintmavg,kmaxtmavg,iindtmavg,xc,zc,nx,ny,nz,nreg,nregprev,zindreg,j1reg,j2reg)
        nzreg = 0
        if(kmaxtmavg>=kmintmavg) nzreg = kmaxtmavg-kmintmavg+1
        njreg = j2reg-j1reg+1
        nreg = nreg*njreg
        nregprev = nregprev*njreg
        nvarreg = 4
c.....Information on the number of points of the region of interest
        CALL MPI_ALLREDUCE(nreg,nregtot,1,mpi_integer,MPI_SUM,MPI_COMM_EDDY,IERR)
        CALL MPI_ALLREDUCE(nzreg,nzgreg,1,mpi_integer,MPI_SUM,MPI_COMM_EDDY,IERR)
        allocate(regtmavg(nreg,nvarreg))
        filereg = trim(hspdir)//'reg_tmavg.bin'
c.....Iteration step of the averaging on the defined region
        call read_regtmavgprm(fregtmavg)
c        write(6,*) 'fregtmavg=',fregtmavg
        nsmplregtmavg = 0
        if(iregtmavg==1) then
          fileregind = trim(hspdir)//'reg_tmavg.ind'
c.....The indices of the region boundaries are written on a file
          call write_regindices(fileregind,kmintmavg,kmaxtmavg,iindtmavg,nz,zindreg,nzgreg,j1reg,j2reg)
c.....The number of points and of variables are written on a file
          if(myrank.eq.0) call write_regtmavg_header(filereg,nregtot,nvarreg)
        endif

        if(iregtmavg>1) then
c.....The stored time-averaged values are read from the restart file
          call read_regvar_tmavg(filereg,regtmavg,nreg,nvarreg,iindtmavg,nz
     &          ,kmintmavg,kmaxtmavg,nsmplregtmavg,nregprev,nregtot)
        endif
      endif

      if(itmavg>0) then
c.....The step of and the area of the time-averaging are read from
c.....file
        call read_tmavgprm(limtmavg,nxt,nyt,nzt,ftmavg,xc,yc,zc,zcg,nx,ny,nz,nzg)
        nvartmavg = 11
        if(isgs.gt.0) nvartmavg = 12
        if(nzt<1) then
          allocate(vartmavg(nxt,nyt,1,nvartmavg))
        else
          allocate(vartmavg(nxt,nyt,nzt,nvartmavg))
        endif
        vartmavg = 0.0
        filetmavg = trim(hspdir)//'vars_tmavg.bin'
c.....Some information on the area and the number of time-averaged variables
c.....is written on a file
c        if(itmavg==1 .AND. myrank==0) call write_tmavg_headers(hspdir,nxt,nyt,nzt,nvartmavg,limtmavg)
c        if(itmavg==1 .AND. myrank==0) call write_tmavg_header(filetmavg,nxt,nyt,nzt,nvartmavg,limtmavg)
        allocate(nsmpltmavg(nvartmavg))
        nsmpltmavg = 0
c.....The grid associated with the area of the time-averaging is written
        if(itmavg.eq.1)
c     %call write_grid_subdom(trim(hspdir)//'grid_tmavg.sp',xu,yv,zw,xc,yc,zc,nx,ny,nz,limtmavg)
     %call write_grid_subdom_sp(trim(hspdir)//'grid_tmavg.sp',xu,yv,zw,xc,yc,zc,nx,ny,nz,limtmavg)
c.....The restart file of the time-averaging is read
c        if(itmavg>1) call read_primevar_tmavg(filetmavg,vartmavg,nxt,nyt,nzt,nvartmavg,nsmpltmavg,limtmavg,nz)
c        if(itmavg>1) call read_primevars_tmavg(hspdir,vartmavg,nxt,nyt,nzt,8,nsmpltmavg,limtmavg,nz)
        if(itmavg>1) call read_primevars_tmavg(hspdir,vartmavg,nxt,nyt,nzt,nvartmavg,nsmpltmavg,limtmavg,nz)

      endif

      if(itVPfield>0) then
c.....The limits of the area for the instantaneous fields are read
c.....also the step is defined
        call read_VPfieldprm(limVPfield,nxVPf,nyVPf,nzVPf,xc,yc,zc,zcg,nx,ny,nz,nzg)
c.....The grid file is generated
c        call write_grid_subdom(trim(hspdir)//'gridVPfield.g',xu,yv,zw,xc,yc,zc,nx,ny,nz,limVPfield)
        call write_grid_subdom_sp(trim(hspdir)//'grid_VPfield.sp',xu,yv,zw,xc,yc,zc,nx,ny,nz,limVPfield)
      endif

c
c.... Instantaneous region files
      VPreg_freq = 0
      INQUIRE(FILE='VPreg.input', EXIST=VPreg_input_exists)
      IF(VPreg_input_exists .eqv. .true.) THEN
        call read_VPreg_prms(VPreg_freq,iVPreg)
      ENDIF

      IF(VPreg_freq>0) THEN
        allocate(iindVPreg(nz,2))
        call VPregion2indices(kminVPreg,kmaxVPreg,iindVPreg,xc,zc,nx,ny,nz,nVPreg,nprevVPreg,zindVPreg,j1VPreg,j2VPreg)
        nzVPreg = 0
        if(kmaxVPreg>=kminVPreg) nzVPreg = kmaxVPreg-kminVPreg+1
        njVPreg = j2VPreg-j1VPreg+1
        nVPreg = nVPreg*njVPreg
        nprevVPreg = nprevVPreg*njVPreg
        VPreg_ind_file = trim(hspdir)//'VPregind.bin'
        CALL MPI_ALLREDUCE(nVPreg,nVPregtot,1,mpi_integer,MPI_SUM,MPI_COMM_EDDY,IERR)
        if(myrank.eq.0) write(6,*) 'nVPregtot=',nVPregtot
        call write_regindices(VPreg_ind_file,kminVPreg,kmaxVPreg,iindVPreg,nz,zindVPreg,nzgVPreg,j1VPreg,j2VPreg)
        call write_grid_VPreg(trim(hspdir)//'gridVPreg.bin',xu,yv,zw,xc
     $       ,yc,zc,nx,ny,nz,iindVPreg,j1VPreg,j2VPreg,kminVPreg,kmaxVPreg,zindVPreg)
      ENDIF
c
c.... CFL input file
      icflave = 0
      INQUIRE(FILE='cflave.input', EXIST=cfl_input_exists)
      IF(cfl_input_exists .eqv. .true.) THEN
        open(unit=10,file='cflave.input',form='formatted')
        read(10,'(A)')
        read(10,*) icflave
        close(10)
      ENDIF

      IF(icflave>0) then
        ALLOCATE(CFLAVE(NX,6))
        CFLAVE = 0.0
        filecflave = trim(hspdir)//'cflave.bin'
        nsmplcfl = 0
        IF(icflave==2) call read_cflave(filecflave,cflave,nx,nsmplcfl)
      ENDIF
c
c.... Stat 1D average file (for flows homogeneous in x and z directions, ie. pipe, channel)
      istat1D = 0
      INQUIRE(FILE='stat1D.input', EXIST=stat1D_input_exists)
      IF(stat1D_input_exists .eqv. .true.) THEN
        call read_stat1D_input(istat1D,itstat1D)
      ENDIF

      IF(istat1D>0) THEN
c if istat1D>0 the 1D statistics are evaluated
        nvarstat1D = 23
        if(isgs>0 .AND. isgs<5) nvarstat1D = 37
        ALLOCATE(stat1Dave(nx,nvarstat1D),stat1D(nx,nvarstat1D))
        stat1Dave = 0.0
        stat1D = 0.0
        nstat1D = 0
        nstat1Dave = 0
c if istat1D=2 there is a file of statistics based on a previous run
        IF(istat1D==2) call read_stat1D_ave(stat1Dave,nx,nvarstat1D,nstat1Dave)
      ENDIF
c
c.... Bulk velocity for calculating dpdz
C       IF(MYRANK==0 .AND. ITYPE(5)==500) THEN
C         ubold = 1.0
C         IF(IFIELD>0) THEN
C           open(unit=10,file='res.ubulk_old',form='formatted')
C           read(10,*) ubold
C           read(10,*) dpdz
C           close(10)
C         ENDIF
C       ENDIF
c
c.... Utau file
      iutau = 0
      INQUIRE(FILE='utau.input', EXIST=utau_input_exists)
      IF(utau_input_exists .eqv. .true.) THEN
c the values of IUTAU and ITUTAU are read from an input file
c if IUTAU>0 the skin friction velocity is evaluated
c ITUTAU defines the step of the evaluation of the skin-friction velocity
        call read_utau_input(iutau,itutau)
      ENDIF

      IF(iutau>0) THEN
        ntutau = ITRES/ITUTAU
        ALLOCATE(timeutau(ntutau),dpdzutau(ntutau))
        utaudim = 1
        IF(icyl==0) utaudim = 2
        ALLOCATE(utau(ntutau,utaudim))
        itu = 0
      ENDIF

! Some set up for the subroutine NOISE
!      call noise(nx,ny,nz,jy1,jy2,kz1,kz2,uo,vo,wo,myrank,0)
c
c-----------------------------------------------------------------------
c
c.....Periodic boundary conditions along the Z direction
      IF(MYRANK==0       .AND.ITYPE(5)==500) MYLEFT  = MYSIZE-1
      IF(MYRANK==MYSIZE-1.AND.ITYPE(6)==500) MYRIGHT = 0

      NSTAR=0
c
C---- Write time stats to file
      IF(IOCLOCK>0) THEN
        clock(99) = sum(clock)
        clock(50) = clock(1) + clock(6) !Allocate arrays
        clock(51) = clock(10) + clock(12) + clock(14) + clock(16) !Tag 3D
        clock(52) = clock(11) + clock(13) + clock(15) !Flag-U
        clock(53) = clock(20) + clock(22) + clock(24) !Geom
        clock(54) = clock(26) + clock(28) + clock(30) !Mtrx
        CALL MPI_REDUCE(CLOCK,CLOCKG,NCLOCKS,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(CLOCK,CLOCKGMIN,NCLOCKS,MTYPE,MPI_MIN,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(CLOCK,CLOCKGMAX,NCLOCKS,MTYPE,MPI_MAX,0,MPI_COMM_EDDY,IERR)
        CLOCKG=CLOCKG/REAL(MYSIZE)

        IF(MYRANK==0) THEN
C        OPEN(UNIT=16,FILE='clock.csv',FORM='FORMATTED')
C        WRITE(16,*) 'Task/Time'
C        CLOSE(16)
          OPEN(UNIT=16,FILE='clock.dat',FORM='FORMATTED',
     &POSITION='APPEND')
          WRITE(16,'(2A)') ' Task/Time  ',
     &                  '      Ave.     Max.     Min.'
          WRITE(16,'(A,3(1x,F8.4))') 'Total                  :',CLOCKG(99),CLOCKGMAX(99),CLOCKGMIN(99)
          WRITE(16,'(A,3(1x,F8.4))') 'Allocate arrays        :',CLOCKG(50),CLOCKGMAX(50),CLOCKGMIN(50)
          WRITE(16,'(A,3(1x,F8.4))') 'Grid and metrics       :',CLOCKG(2),CLOCKGMAX(2),CLOCKGMIN(2)
          WRITE(16,'(A,3(1x,F8.4))') 'Setup                  :',CLOCKG(3),CLOCKGMAX(3),CLOCKGMIN(3)
          WRITE(16,'(A,3(1x,F8.4))') 'Initialization         :',CLOCKG(4),CLOCKGMAX(4),CLOCKGMIN(4)
          IF(ibm/=0) THEN
            WRITE(16,'(A,3(1x,F8.4))') 'Stl Input              :',CLOCKG(5),CLOCKGMAX(5),CLOCKGMIN(5)
            WRITE(16,'(A,3(1x,F8.4))') 'Read STL               :',CLOCKG(7),CLOCKGMAX(7),CLOCKGMIN(7)
            WRITE(16,'(A,3(1x,F8.4))') 'Calc DS triangles      :',CLOCKG(8),CLOCKGMAX(8),CLOCKGMIN(8)
            WRITE(16,'(A,3(1x,F8.4))') 'Lim Im. Bound.         :',CLOCKG(9),CLOCKGMAX(9),CLOCKGMIN(9)
            WRITE(16,'(A,3(1x,F8.4))') 'Tag3D (U,V,W,P)        :',CLOCKG(51),CLOCKGMAX(51),CLOCKGMIN(51)
            WRITE(16,'(A,3(1x,F8.4))') 'FlagU(U,V,W)           :',CLOCKG(52),CLOCKGMAX(52),CLOCKGMIN(52)
            WRITE(16,'(A,3(1x,F8.4))') 'FlagP                  :',CLOCKG(17),CLOCKGMAX(17),CLOCKGMIN(17)
            WRITE(16,'(A,3(1x,F8.4))') 'RefreshFlag            :',CLOCKG(70),CLOCKGMAX(70),CLOCKGMIN(70)
            WRITE(16,'(A,3(1x,F8.4))') 'I/O                    :',CLOCKG(80),CLOCKGMAX(80),CLOCKGMIN(80)
            WRITE(16,'(A,3(1x,F8.4))') 'Geom(U,V,W)            :',CLOCKG(53),CLOCKGMAX(53),CLOCKGMIN(53)
            WRITE(16,'(A,3(1x,F8.4))') 'Mtrx(U,V,W)            :',CLOCKG(54),CLOCKGMAX(54),CLOCKGMIN(54)
            WRITE(16,'(A,3(1x,F8.4))') 'GeomP                  :',CLOCKG(32),CLOCKGMAX(32),CLOCKGMIN(32)
            WRITE(16,'(A,3(1x,F8.4))') 'MtrxP                  :',CLOCKG(34),CLOCKGMAX(34),CLOCKGMIN(34)
            WRITE(16,'(A,3(1x,F8.4))') 'Set U forcing          :',CLOCKG(35),CLOCKGMAX(35),CLOCKGMIN(35)
            WRITE(16,'(A,3(1x,F8.4))') 'Geom_WM(U,V,W)         :',CLOCKG(18),CLOCKGMAX(18),CLOCKGMIN(18)
          ENDIF
          IF(isgs/=0) THEN
            WRITE(16,'(A,3(1x,F8.4))') 'Init. turb. visc.      :',CLOCKG(37),CLOCKGMAX(37),CLOCKGMIN(37)
          ENDIF
          WRITE(16,'(A,3(1x,F8.4))') 'Setup for output       :',CLOCKG(39),CLOCKGMAX(39),CLOCKGMIN(39)
          CLOSE(16)
        ENDIF
      ENDIF

      clock = 0.0
      TIME_ITERATIONS = tclock()

      if(ifield==1 .and. inoise==1) then
c          call add_random_noise_2(uo(:,:,:),nx,ny,nz,1,xu,26.0)
c          call add_random_noise_2(vo(:,:,:),nx,ny,nz,1,xc,26.0)
c          call add_random_noise_2(wo(:,:,:),nx,ny,nz,1,xc,26.0)

          call add_random_noise_2(uo(:,:,:),nx,ny,nz,1,xu,30.0)
          call add_random_noise_2(vo(:,:,:),nx,ny,nz,1,xc,30.0)
          call add_random_noise_2(wo(:,:,:),nx,ny,nz,1,xc,30.0)
      endif

c$$$      call PERIODIC_2D_AVERAGE_RMS_1CPU(NX,NY,NZ,NZG,NBD,ICYCLE,time,DTM1,XC,XU,YC,YV,ZCG,ZWG,xc_car,
c$$$     &                              yc_car,xu_car,yu_car,xv_car,yv_car,P,VO,UO,WO,DENS)
c$$$      if(myrank==0) write(*,*) "DONE PERIODIC_2D_AVERAGE_RS"

c
C-----------------------------------------------------------------------
C     begin main loop
C-----------------------------------------------------------------------
C
      do icycle=1+nstep,itmax+nstep
c
      if(mod(icycle-1,itscr)==0) clock=0.0

      TIME_ITERATION = tclock()
c
c.....time step
c

      clocktemp = tclock()

      if((tclock()-begin_time) .ge. ttotal*background_time) then
      IF(MYRANK==0) WRITE(*,*) 'CHANGING SAVE RES'
      SAVE_RES=.TRUE.
      ttotal=ttotal+1.0
      endif


      if(icfl==1) then
c
c.....courant number calculation for running at cfl=constant
c.....some statistics on the CFL number are written on screen
c
        CALL CALCFL(UO,VO,WO,TV,DP,NX,NY,NZ,CFLM,ICYCLE)
        if(icflave>0) then
          call COMPUTE_CFLAVE(cflave,uo,vo,wo,tv,dp,nx,ny,nz,nsmplcfl)
        endif
        DTM1=CFLC/CFLM
        TSTEP=DTM1
        if(tstep.le.1.e-06) then
           call mpi_finalize(ierr)
           stop
        endif

      ELSE

        DTM1=TSTEP

      ENDIF

      clock(1) = clock(1) + tclock() - clocktemp

!      CALL FORTERM(WO(:,:,1),1.0,UBOLD,DTM1,NX,NY)
c
c-----------------------------------------------------------------------
c     setup external forcing
c-----------------------------------------------------------------------
c
c.....time advancement: ischm=1, Adams-Bashforth; ischm=3, Runge-Kutta
c
      if(idens.eq.1)CALL GRAVITY(time)

      do is=1,ischm

        ALFXDT=ALF(IS)*DTM1
        GAMXDT=GAM(IS)*DTM1
        RHOXDT=RHO(IS)*DTM1

        TLEVEL=TIME+SUM(ALF(1:IS))*DTM1

!        CALL FORTERM(WO(:,:,1),0.0,UBOLD,ALFXDT,NX,NY)
C-----------------------------------------------------------------------
C                                               predictor-corrector step
C-----------------------------------------------------------------------
c
        clocktemp1 = tclock()


        IF(IBM>1) THEN
c
c...Rotate body
c...The angular velocity is defined by the function rotspeed
c...Currently a rotation around Z is considered !!!!!!
          clocktemp = tclock()
          call rbm(unvect,vertex,vertexc,nfacet,mbd,nbd,alfxdt,tlevel)
          IF(ITBDY>0 .AND. IVRTX==1) call rbm_nodes(trvtx,vtxnrm,nvtx,mbd,nbd,alfxdt,tlevel)
!!!!!!          clock(21) = clock(21) + tclock()-clocktemp
          clock(51) = clock(51) + tclock()-clocktemp

          IF(IMBDCMP>0 .AND. MYSIZE>1) THEN
c...In the case the immersed-boundary is decomposed
            call MOVTRI(UNVECT,VERTEX,VERTEXC,AREAF,TRINO,NFACET,NFACETMAX,ZC,NZ,MBD,NBD)
          ENDIF

c...Calculate new forcing points stencil considering the moving boundaries
          nimu=limu(mbd,1)
          nimv=limv(mbd,1)
          nimw=limw(mbd,1)
          nimp=limp(mbd,1)
          if(isgs.ne.0)nimtv=limtv(mbd,1)
          DO ibd=mbd,nbd

            ilb = lb(ibd)+1
            ile = lb(ibd)+mb(ibd)
c
c...        Calculate distance of triangles
            clocktemp = tclock()
!!!!!!!            call calc_ds_triangles(vertex,vertexc,nfacet,ilb,ile,dsb(ibd,1),xu,nx,icyl)
            call calc_ds_triangles(vertex,vertexc,nfacet,ilb,ile,dsb(ibd,:),xu,nx,icyl)
            clock(22) = clock(22) + tclock()-clocktemp
c
c...        Calculate limits of subdomain enclosing immersed body
c...They are stored in IBMIN,IBMAX,JBMIN,JBMAX,KBMIN,KBMAX
            clocktemp = tclock()
!!!!!!            call limimb(xu,yv(jys),zwg(kzs),zcg(kzs),nx,nyl,nzl,vertexc,nfacet,ibd)
            call limimb
     %(xu,yv(jys),zwg(kzs),zcg(kzs),nx,nyl,nzl,vertexc(:,ilb:ile),mb(ibd),ibd)

            clock(23) = clock(23) + tclock()-clocktemp
c
c... Find forcing points and build stencil for U
            clocktemp = tclock()
!!!!!!            call tag3d(xu_car(:,jys),yu_car(:,jys),zcg(kzs),nx,nyl,nzl
!!!!!!     &           ,vertex,vertexc,unvect,nfacet,ilb,ile,dsb(ibd,1),flag,ibd,0)
            call tag3d_mod
     %(xu_car(:,jys),yu_car(:,jys),zcg(kzs),nx,nyl,nzl
     %,vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile)
     %,mb(ibd),dsb(ibd,:),flag,ibd,0)
            clock(24) = clock(24) + tclock()-clocktemp
c
            clocktemp = tclock()
            CALL taguo(flag,flaguo,nx,nyl,nzl,nbd,ibd)
            clock(25) = clock(25) + tclock() - clocktemp
c
c... Find forcing points and build stencil for V
            clocktemp = tclock()
!!!!!!            call tag3d(xv_car(:,jys),yv_car(:,jys),zcg(kzs),nx,nyl,nzl
!!!!!!     &           ,vertex,vertexc,unvect,nfacet,ilb,ile,dsb(ibd,1),flag,ibd,1)
            call tag3d_mod
     %(xv_car(:,jys),yv_car(:,jys),zcg(kzs),nx,nyl,nzl
     %,vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile)
     %,mb(ibd),dsb(ibd,:),flag,ibd,1)
            clock(28) = clock(28) + tclock()-clocktemp

            clocktemp = tclock()
            call taguo(flag,flagvo,nx,nyl,nzl,nbd,ibd)
            clock(29) = clock(29) + tclock()-clocktemp
c
c... Find forcing points and build stencil for W
            clocktemp = tclock()
            icom = 2
            IF(ISGS>0) icom=1
!!!!!!            call tag3d(xc_car(:,jys),yc_car(:,jys),zwg(kzs),nx,nyl,nzl
!!!!!!     &           ,vertex,vertexc,unvect,nfacet,ilb,ile,dsb(ibd,1),flag,ibd,icom)
            call tag3d_mod
     %(xc_car(:,jys),yc_car(:,jys),zwg(kzs),nx,nyl,nzl
     %,vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile)
     %,mb(ibd),dsb(ibd,:),flag,ibd,icom)
            clock(32) = clock(32) + tclock()-clocktemp

            clocktemp = tclock()
            call taguo(flag,flagwo,nx,nyl,nzl,nbd,ibd)
            clock(33) = clock(33) + tclock()-clocktemp
c
c... Find forcing points and build stencil for TV
            IF(ISGS>0) THEN
              clocktemp = tclock()
!!!!!!              call tag3d(xc_car(:,jys),yc_car(:,jys),zcg(kzs),nx,nyl,nzl
!!!!!!     &         ,vertex,vertexc,unvect,nfacet,ilb,ile,dsb(ibd,1),flag,ibd,2)
              call tag3d_mod
     %(xc_car(:,jys),yc_car(:,jys),zcg(kzs),nx,nyl,nzl
     %,vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile)
     %,mb(ibd),dsb(ibd,:),flag,ibd,2)
              clock(32) = clock(32) + tclock()-clocktemp

              clocktemp = tclock()
              CALL taguo(flag,flagtv,nx,nyl,nzl,nbd,ibd)
              clock(33) = clock(33) + tclock()-clocktemp
            ENDIF

            clocktemp = tclock()
!!!!!!            nimu = 0
!!!!!!            nimv = 0
!!!!!!            nimw = 0
!!!!!!            nimtv = 0
            icom = 0

            clocktemp = tclock()
!!!!!!            call geom
!!!!!!     &(xu,yc(jys),xu_car(:,jys),yu_car(:,jys),zcg(kzs),vertex,vertexc,unvect
!!!!!!     &,limu(ibd,1),mimu(ibd,1),iu,ju,ku,nimu,xnu,ynu,znu,nxu,nyu,nzu
!!!!!!     &,mrku,diru,fp,flag,flaguo,nx,nyl,nzl,nbd,nfacet,ibd,icom,nflu,nflumax,1,icycle,tlevel)
ccc            call geom
            call geom_mod
     &(xu,yc(jys),xu_car(:,jys),yu_car(:,jys),zcg(kzs)
     &,vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile)
     &,limu(ibd,:),mimu(ibd,:),iu,ju,ku,nimu,xnu,ynu,znu,nxu,nyu,nzu
     &,mrku,diru,fpu,flag,flaguo,nx,nyl,nzl,nbd,mb(ibd),ibd,icom,nflu,nflumax(ibd),1,icycle,tlevel)
            clock(26) = clock(26) + tclock() - clocktemp

            clocktemp = tclock()
            call mtrx(xu,yc(jys),zcg(kzs),limu(ibd,1),sum(mimu(ibd,:)),mrku
     &           ,iu,ju,ku,xnu,ynu,znu,nxu,nyu,nzu,umtrx,uindx,nx,nyl,nzl,nbd,ibd,1)
            clock(27) = clock(27) + tclock() - clocktemp

            icom = 1
            clocktemp = tclock()
!!!!!!            call geom
!!!!!!     &(xc,yv(jys),xv_car(:,jys),yv_car(:,jys),zcg(kzs),vertex,vertexc,unvect
!!!!!!     &,limv(ibd,1),mimv(ibd,1),iv,jv,kv,nimv,xnv,ynv,znv,nxv,nyv,nzv
!!!!!!     &,mrkv,dirv,fp,flag,flagvo,nx,nyl,nzl,nbd,nfacet,ibd,icom,nflu,nflvmax,2,icycle,tlevel)
ccc            call geom
            call geom_mod
     &(xc,yv(jys),xv_car(:,jys),yv_car(:,jys),zcg(kzs)
     &,vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile)
     &,limv(ibd,:),mimv(ibd,:),iv,jv,kv,nimv,xnv,ynv,znv,nxv,nyv,nzv
     &,mrkv,dirv,fpv,flag,flagvo,nx,nyl,nzl,nbd,mb(ibd),ibd,icom,nflu,nflvmax(ibd),2,icycle,tlevel)
            clock(30) = clock(30) + tclock()-clocktemp
c
            clocktemp = tclock()
            call mtrx(xc,yv(jys),zcg(kzs),limv(ibd,1),sum(mimv(ibd,:)),mrkv
     &           ,iv,jv,kv,xnv,ynv,znv,nxv,nyv,nzv,vmtrx,vindx,nx,nyl,nzl,nbd,ibd,1)
            clock(31) = clock(31) + tclock()-clocktemp

            clocktemp = tclock()
!!!!!!            call geom
!!!!!!     &(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zwg(kzs),vertex,vertexc,unvect
!!!!!!     &,limw(ibd,1),mimw(ibd,1),iw,jw,kw,nimw,xnw,ynw,znw,nxw,nyw,nzw
!!!!!!     &,mrkw,dirw,fp,flag,flagwo,nx,nyl,nzl,nbd,nfacet,ibd,icom,nflu,nflwmax,3,icycle,tlevel)
ccc            call geom
            call geom_mod
     &(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zwg(kzs)
     &,vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile)
     &,limw(ibd,:),mimw(ibd,:),iw,jw,kw,nimw,xnw,ynw,znw,nxw,nyw,nzw
     &,mrkw,dirw,fpw,flag,flagwo,nx,nyl,nzl,nbd,mb(ibd),ibd,icom,nflu,nflwmax(ibd),3,icycle,tlevel)
            clock(34) = clock(34) + tclock()-clocktemp
c
            clocktemp = tclock()
            call mtrx(xc,yc(jys),zwg(kzs),limw(ibd,1),sum(mimw(ibd,:)),mrkw
     &           ,iw,jw,kw,xnw,ynw,znw,nxw,nyw,nzw,wmtrx,windx,nx,nyl,nzl,nbd,ibd,1)
            clock(35) = clock(35) + tclock()-clocktemp

            IF(ISGS>0) THEN
              clocktemp = tclock()
!!!!!!              call geom
!!!!!!     &(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs),vertex,vertexc,unvect
!!!!!!     &,limtv(ibd,:),mimtv(ibd,:),itv,jtv,ktv,nimtv,xntv,yntv,zntv,nxtv,nytv,nztv
!!!!!!     &,mrktv,dirtv,fp,flag,flagtv,nx,nyl,nzl,nbd,nfacet,ibd,icom,nflu,nfltvmax,5,icycle,tlevel)
ccc              call geom
              call geom_mod
     &(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs)
     &,vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile)
     &,limtv(ibd,:),mimtv(ibd,:),itv,jtv,ktv,nimtv,xntv,yntv,zntv,nxtv,nytv,nztv
     &,mrktv,dirtv,fptv,flag,flagtv,nx,nyl,nzl,nbd,mb(ibd),ibd,icom,nflu,nfltvmax(ibd),5,icycle,tlevel)
              clock(34) = clock(34) + tclock()-clocktemp

              clocktemp = tclock()
              call mtrx(xc,yc(jys),zcg(kzs),limtv(ibd,1),sum(mimtv(ibd,:)),mrktv
     &,itv,jtv,ktv,xntv,yntv,zntv,nxtv,nytv,nztv,tvmtrx,tvindx,nx,nyl,nzl,nbd,ibd,1)
              clock(35) = clock(35) + tclock()-clocktemp
            ENDIF

c... Find forcing points and buld stencil for p
            clocktemp = tclock()
c... The flagp subroutine tags the pressure boundary points
            call flagp(flagpo,flaguo,flagvo,flagwo,nx,nyl,nzl,nbd,ibd)
c            call flagp1(flagpo,flaguo,flagvo,flagwo,nx,nyl,nzl,nbd,iu,ju,ku,iv,jv,kv,iw,jw,kw
c     &         ,limu(ibd,1),sum(mimu(ibd,:)),limv(ibd,1),sum(mimv(ibd,:)),limw(ibd,1),sum(mimw(ibd,:)),ibd)
            clock(36) = clock(36) + tclock()-clocktemp

            clocktemp = tclock()
!!!!!!            nimp = 0
            icom = 2
!!!!!!            IF(IS==ISCHM .AND. ((ITBDY>0 .AND. MOD(ICYCLE,ITBDY)==0) .OR. (ICALF>0 .AND. MOD(ICYCLE,ITCALF)==0)) ) icom=1
!!!!!!            call geom
!!!!!!     &(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs),vertex,vertexc,unvect
!!!!!!     &,limp(ibd,1),mimp(ibd,1),ip,jp,kp,nimp,xnp,ynp,znp,nxp,nyp,nzp
!!!!!!     &,mrkp,dirp,fp,flag,flagpo,nx,nyl,nzl,nbd,nfacet,ibd,icom,nflp,nflpmax,4,icycle,tlevel)
ccc            call geom
            call geom_mod
     &(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs)
     &,vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile)
     &,limp(ibd,:),mimp(ibd,:),ip,jp,kp,nimp,xnp,ynp,znp,nxp,nyp,nzp
     &,mrkp,dirp,fp,flag,flagpo,nx,nyl,nzl,nbd,mb(ibd),ibd,icom,nflp,nflpmax(ibd),4,icycle,tlevel)
            clock(37) = clock(37) + tclock()-clocktemp

            clocktemp = tclock()
            call mtrx(xc,yc(jys),zcg(kzs),limp(ibd,1),sum(mimp(ibd,:)),mrkp
     &           ,ip,jp,kp,xnp,ynp,znp,nxp,nyp,nzp,pmtrx,pindx,nx,nyl,nzl,nbd,ibd,0)
            clock(38) = clock(38) + tclock()-clocktemp

          ENDDO
c
c...Set surface boundary velocity and pressure gradients for moving bodies
c...The coordinates are always Cartesian, but the velocity components
c...UIM,VIM,WIM can be Cartesian or cylindrical
          clocktemp = tclock()
          DO ibd=mbd,nbd

            DO im=limu(ibd,1)+1,limu(ibd,1)+sum(mimu(ibd,:))
              uim(im) = ubd(xnu(im),ynu(im),znu(im),tlevel,ibd)
            ENDDO
C
            DO jm=limv(ibd,1)+1,limv(ibd,1)+sum(mimv(ibd,:))
              vim(jm) = vbd(xnv(jm),ynv(jm),znv(jm),tlevel,ibd)
            ENDDO

            DO km=limw(ibd,1)+1,limw(ibd,1)+sum(mimw(ibd,:))
              wim(km) = wbd(xnw(km),ynw(km),znw(km),tlevel,ibd)
            ENDDO

            IF(ISGS>0) THEN
              DO im=limtv(ibd,1)+1,limtv(ibd,1)+sum(mimtv(ibd,:))
                tvim(im) = 0.0
              ENDDO
            ENDIF
          ENDDO

          clock(56) = clock(56) + tclock() - clocktemp

c....The subroutines FLAG2INT, IMB_DOMY2Z and INTIND_DOMY2z are used to
c....convert the tagging variables from the decomposition along Y to the
c....one along Z
          IF(idomy==1 .AND. mysize>1) THEN
            DO ibd=mbd,nbd
              limuint(ibd) = sum(mimuint(1:ibd-1))
              call flag2int(flaguo(:,:,:,ibd),nx,nyl,nzl,ibd,iuint,juint,kuint,limuint(ibd),mimuint(ibd))
              limvint(ibd) = sum(mimvint(1:ibd-1))
              call flag2int(flagvo(:,:,:,ibd),nx,nyl,nzl,ibd,ivint,jvint,kvint,limvint(ibd),mimvint(ibd))
              limwint(ibd) = sum(mimwint(1:ibd-1))
              call flag2int(flagwo(:,:,:,ibd),nx,nyl,nzl,ibd,iwint,jwint,kwint,limwint(ibd),mimwint(ibd))
              limpint(ibd) = sum(mimpint(1:ibd-1))
              call flag2int(flagpo(:,:,:,ibd),nx,nyl,nzl,ibd,ipint,jpint,kpint,limpint(ibd),mimpint(ibd))
            ENDDO
          ENDIF

          if(idomy==1 .AND. mysize>1) then
!!!!!!            call imb_domy2z(iu,ju,ku,xnu,ynu,znu,nxu,nyu,nzu,mrku,diru,fp
!!!!!!     &           ,uim,umtrx,uindx,zcg,nzg,zc,nz,nyl,limu(mbd,1),mimu(mbd,1),nbd-mbd+1,nflu,0)
            call imb_domy2z(iu,ju,ku,xnu,ynu,znu,nxu,nyu,nzu,mrku,diru,fpu
     &           ,uim,umtrx,uindx,zcg,nzg,zc,nz,nyl,limu(mbd:nbd,:),mimu(mbd:nbd,:),nbd-mbd+1,nflu,0)
!!!!!!            call imb_domy2z(iv,jv,kv,xnv,ynv,znv,nxv,nyv,nzv,mrkv,dirv,fp
!!!!!!     &           ,vim,vmtrx,vindx,zcg,nzg,zc,nz,nyl,limv(mbd,1),mimv(mbd,1),nbd-mbd+1,nflu,0)
            call imb_domy2z(iv,jv,kv,xnv,ynv,znv,nxv,nyv,nzv,mrkv,dirv,fpv
     &           ,vim,vmtrx,vindx,zcg,nzg,zc,nz,nyl,limv(mbd:nbd,:),mimv(mbd:nbd,:),nbd-mbd+1,nflu,0)
!!!!!!            call imb_domy2z(iw,jw,kw,xnw,ynw,znw,nxw,nyw,nzw,mrkw,dirw,fp
!!!!!!     &           ,wim,wmtrx,windx,zcg,nzg,zc,nz,nyl,limw(mbd,1),mimw(mbd,1),nbd-mbd+1,nflu,0)
            call imb_domy2z(iw,jw,kw,xnw,ynw,znw,nxw,nyw,nzw,mrkw,dirw,fpw
     &           ,wim,wmtrx,windx,zcg,nzg,zc,nz,nyl,limw(mbd:nbd,:),mimw(mbd:nbd,:),nbd-mbd+1,nflu,0)
!!!!!!            call imb_domy2z(ip,jp,kp,xnp,ynp,znp,nxp,nyp,nzp,mrkp,dirp,fp
!!!!!!     &           ,pim,pmtrx,pindx,zcg,nzg,zc,nz,nyl,limp(mbd,1),mimp(mbd,1),nbd-mbd+1,nflp,1)
            call imb_domy2z(ip,jp,kp,xnp,ynp,znp,nxp,nyp,nzp,mrkp,dirp,fp
     &           ,pim,pmtrx,pindx,zcg,nzg,zc,nz,nyl,limp(mbd:nbd,:),mimp(mbd:nbd,:),nbd-mbd+1,nflp,1)

            call intind_domy2z(iuint,juint,kuint,limuint,mimuint,nbd,nyl,nz)
            call intind_domy2z(ivint,jvint,kvint,limvint,mimvint,nbd,nyl,nz)
            call intind_domy2z(iwint,jwint,kwint,limwint,mimwint,nbd,nyl,nz)
            call intind_domy2z(ipint,jpint,kpint,limpint,mimpint,nbd,nyl,nz)

          endif

        ENDIF
        clock(2) = clock(2) + tclock() - clocktemp1

! A noise is added in a particular area of the computational domain
!        call noise(nx,ny,nz,jy1,jy2,kz1,kz2,uo,vo,wo,myrank,1)

        clocktemp = tclock()
        CALL RHS(UO,VO,WO,DENS,US,VS,WS,TV,UB,VB,WB,DP,NX,NY,NZ,XC,XU,YC,YV,TIME)
        clock(3) = clock(3) + tclock() - clocktemp

        clocktemp = tclock()
        CALL REFRESHBC(US,NX*NY,NZ)
        CALL REFRESHBC(VS,NX*NY,NZ)
        CALL REFRESHBC(WS,NX*NY,NZ)
        CALL REFRESHBC(UA,NX*NY,NZ)
        CALL REFRESHBC(VA,NX*NY,NZ)
        CALL REFRESHBC(WA,NX*NY,NZ)
        CALL REFRESHBC(UB,NX*NY,NZ)
        CALL REFRESHBC(VB,NX*NY,NZ)
        CALL REFRESHBC(WB,NX*NY,NZ)
c N.B. UA,VA,WA are the values of RHS at the time level L-2,
c whereas UB,VB,WB are the values at the time level L-1
        clock(8) = clock(8) + tclock() - clocktemp
C
C.....Predictor step + boundary conditions for us
        clocktemp = tclock()
!         IF(IDENS.EQ.0) THEN
!           CALL PREDICTOR(UO,VO,WO,P,UA,VA,WA,UB,VB,WB,US,VS,WS,
!      &       ALFXDT,GAMXDT,RHOXDT,NX,NY,NZ,clock(40),9)
!         ELSE
!          DENSO=DENS
!          CALL RHS_DENSITY(UO,VO,WO,DENS,RHB,NX,NY,NZ)
!          CALL REFRESHBC(RHB,NX*NY,NZ)
!          CALL DENSITY(DENS,RHA,RHB,XC,ALFXDT,GAMXDT,RHOXDT,NX,NY,NZ)
!          CALL BOUNDARY_DENS(DENS,XC,NX,NY,NZ)
!          CALL REFRESHBC(DENS,NX*NY,NZ)
!          CALL PREDICTORD(UO,VO,WO,P,UA,VA,WA,UB,VB,WB,US,VS,WS,DENSO,XU,

             CALL PREDICTORD(UO,VO,WO,P,UA,VA,WA,UB,VB,WB,US,VS,WS,DENS,XU,XC,YC,YV,
     &       ZWG,ZCG,ALFXDT,GAMXDT,RHOXDT,NX,NY,NZ,NZG,TIME,clock(40),9)

!         ENDIF
c RHS at the time level L-1 in UA,VA,WA
c implicit terms in UB,VB,WB
c provisional solution in US,VS,WS
        clock(4) = clock(4) + tclock() - clocktemp

c.....compute forcing in momentum
        clocktemp1 = tclock()
        IF(ibm/=0.and.(implx/=0.or.imply/=0)) THEN
c
c.....predicted velocity from explicit step
c
          UB = US+0.5*ALFXDT*UB
          VB = VS+0.5*ALFXDT*VB
          WB = WS+0.5*ALFXDT*WB
c in UB,VB,WB is evaluated the explicit velocity
c
c.....other boundary conditions
c
          CALL BOUNDARY(UB,VB,WB,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TLEVEL)

C-----------------------------------------------------------------------
C.....refresh block interfaces (and periodic boundary in z direction)
C-----------------------------------------------------------------------
          CALL REFRESHBC(UB,NX*NY,NZ)
          CALL REFRESHBC(VB,NX*NY,NZ)
          CALL REFRESHBC(WB,NX*NY,NZ)
c
c.....inlet and outlet boundary conditions
c
          CALL BOUNDINOUT(UB,VB,WB,UO,VO,WO,XU,XC,YC,ZW,ZC,ALFXDT,TLEVEL,
     & NX,NY,NZ,1,planeU,planeV,planeW)
	  CALL BOUNDARY(UB,VB,WB,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TIME)  


C-----momentum forcing from an explicit step
C-----------------------------------------------------------------------

c.....compute forcing in x momentum
          call momforc_mod2(us,ub,xu,yc,zc,nx,ny,nz,iu,ju,ku,umtrx,uindx
     &           ,uim,mrku,diru,nxu,nyu,nzu,nflumax,limu,mimu,1,1,nbd,sum(mimu))
c.....compute forcing in y momentum
          call momforc_mod2(vs,vb,xc,yv,zc,nx,ny,nz,iv,jv,kv,vmtrx,vindx
     &           ,vim,mrkv,dirv,nxv,nyv,nzv,nflvmax,limv,mimv,1,2,nbd,sum(mimv))
c.....compute forcing in z momentum
          call momforc_mod2(ws,wb,xc,yc,zw,nx,ny,nz,iw,jw,kw,wmtrx,windx
     &           ,wim,mrkw,dirw,nxw,nyw,nzw,nflwmax,limw,mimw,1,3,nbd,sum(mimw))

c
c......Set velocity inside stationary bodies
          IF(idomy==1 .AND. mysize>1) THEN
!!!!!!            DO ibd=1,mbd-1
            DO ibd=1,mbdn-1
              call velinterior1
     %(us,vs,ws,ub,vb,wb,iuint,juint,kuint,ivint,jvint,kvint,
     %iwint,jwint,kwint,xu,yv,zw,xc,yc,zc,nx,ny,nz,
     %limuint(ibd),mimuint(ibd),limvint(ibd),mimvint(ibd),
     %limwint(ibd),mimwint(ibd),tlevel,ibd,0,1)
            ENDDO
          ELSE
!!!!!!            DO ibd=1,mbd-1
            DO ibd=1,mbdn-1
              call velinterior
     %(us,vs,ws,ub,vb,wb,flaguo,flagvo,flagwo,xu,yv,zw,
     %xc,yc,zc,nx,ny,nz,nbd,tlevel,ibd,0,1)
            ENDDO
c            call obstacle(us,vs,ws,p,nx,ny,nz,0)
          ENDIF
c
c.....Set velocity inside moving bodies
          IF(idomy==1 .AND. mysize>1) THEN
!!!!!!            DO ibd=mbd,nbd
            DO ibd=mbdn,nbd
              call velinterior1
     &(us,vs,ws,ub,vb,wb,iuint,juint,kuint,ivint,jvint,kvint
     &,iwint,jwint,kwint,xu,yv,zw,xc,yc,zc,nx,ny,nz
     &,limuint(ibd),mimuint(ibd),limvint(ibd),mimvint(ibd)
     &,limwint(ibd),mimwint(ibd),tlevel,ibd,1,1)
            ENDDO
          ELSE
!!!!!!            DO ibd=mbd,nbd
            DO ibd=mbdn,nbd
              call velinterior
     &(us,vs,ws,ub,vb,wb,flaguo,flagvo,flagwo,xu,yv,zw
     &,xc,yc,zc,nx,ny,nz,nbd,tlevel,ibd,1,1)
            ENDDO
         ENDIF

      ENDIF
      clock(5) = clock(5) + tclock() - clocktemp1

c
c.....solve tridiagonal symstems in y and/or x directions
c
        clocktemp = tclock()
        IF(IBM==0 .AND. (IMPLX==1 .OR. IMPLCX==1)) THEN
          UB(IX1,:,:) = US(IX1,:,:)+0.5*ALFXDT*UB(IX1,:,:)
          UB(IX2,:,:) = US(IX2,:,:)+0.5*ALFXDT*UB(IX2,:,:)
          VB(IX1,:,:) = VS(IX1,:,:)+0.5*ALFXDT*VB(IX1,:,:)
          VB(IX2,:,:) = VS(IX2,:,:)+0.5*ALFXDT*VB(IX2,:,:)
          WB(IX1,:,:) = WS(IX1,:,:)+0.5*ALFXDT*WB(IX1,:,:)
          WB(IX2,:,:) = WS(IX2,:,:)+0.5*ALFXDT*WB(IX2,:,:)
          CALL BOUNDARY(UB,VB,WB,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TLEVEL)
         CALL BOUNDARY(UB,VB,WB,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TLEVEL)
c        write(*,*), "Eddy6.f u(nx-1,1,nz-1)
c       &  u(nx-1,ny-1,nz-1) 2574",WO(nx-1,1,nz-1), WO(nx-1,ny-1,nz-1)  
        ENDIF
        clock(8) = clock(8) + tclock() - clocktemp

        clocktemp = tclock()
        IF(imply==1 .OR. implcy==1) THEN
            CALL INVERSEY(US,VS,WS,TV,VO,ALFXDT,GAMXDT,RHOXDT,NX,NY,NZ)
        ENDIF

        IF(implx==1 .OR. implcx==1) THEN
          IF(ITYPE(1)==500 .AND. ITYPE(2)==500 .AND. implxdcmp==0 .AND. implydcmp==0) THEN
            CALL INVERSEX_PER(US,VS,WS,TV,UO,ALFXDT,GAMXDT,RHOXDT,NX,NY,NZ)
          ELSE
            CALL INVERSEX(US,VS,WS,TV,UO,UB,VB,WB,ALFXDT,GAMXDT,RHOXDT,NX,NY,NZ)
          ENDIF
        ENDIF

        clock(6) = clock(6) + tclock() - clocktemp
c
c....compute momentum forcing for explicit time advancement
c....in the case of explicit treatment
        clocktemp1 = tclock()
        IF(ibm/=0.and.implx==0.and.imply==0) THEN
c
c.....other boundary conditions
c
          CALL BOUNDARY(US,VS,WS,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TLEVEL)
c
C-----------------------------------------------------------------------
C.....refresh block interfaces (and periodic boundary in z direction)
C-----------------------------------------------------------------------
          CALL REFRESHBC(US,NX*NY,NZ)
          CALL REFRESHBC(VS,NX*NY,NZ)
          CALL REFRESHBC(WS,NX*NY,NZ)

c.....inlet and outlet boundary conditions
          CALL BOUNDINOUT(US,VS,WS,UO,VO,WO,XU,XC,YC,ZW,ZC,ALFXDT,TLEVEL,
     &    NX,NY,NZ,1,planeU,planeV,planeW)
          CALL BOUNDARY(US,VS,WS,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TIME)  

c.....moving bodies

c.....compute forcing in x momentum
          call momforc_mod2(us,ub,xu,yc,zc,nx,ny,nz,iu,ju,ku,umtrx,uindx
     &           ,uim,mrku,diru,nxu,nyu,nzu,nflumax,limu,mimu,0,1,nbd,sum(mimu))
c.....compute forcing in y momentum
          call momforc_mod2(vs,vb,xc,yv,zc,nx,ny,nz,iv,jv,kv,vmtrx,vindx
     &           ,vim,mrkv,dirv,nxv,nyv,nzv,nflvmax,limv,mimv,0,2,nbd,sum(mimv))
c.....compute forcing in z momentum
          call momforc_mod2(ws,wb,xc,yc,zw,nx,ny,nz,iw,jw,kw,wmtrx,windx
     &           ,wim,mrkw,dirw,nxw,nyw,nzw,nflwmax,limw,mimw,0,3,nbd,sum(mimw))

c
c...Set velocity inside moving bodies
          IF(idomy==1 .AND. mysize>1) THEN
!!!!!!            DO ibd=mbd,nbd
            DO ibd=mbdn,nbd
              call velinterior1
     &(us,vs,ws,ub,vb,wb,iuint,juint,kuint,ivint,jvint,kvint
     &,iwint,jwint,kwint,xu,yv,zw,xc,yc,zc,nx,ny,nz
     &,limuint(ibd),mimuint(ibd),limvint(ibd),mimvint(ibd),limwint(ibd),mimwint(ibd)
     &,tlevel,ibd,1,0)
            ENDDO
          ELSE
!!!!!!            DO ibd=mbd,nbd
            DO ibd=mbdn,nbd
              call velinterior
     &(us,vs,ws,ub,vb,wb,flaguo,flagvo,flagwo,xu,yv,zw
     &,xc,yc,zc,nx,ny,nz,nbd,tlevel,ibd,1,0)
            ENDDO
          ENDIF
c
c.....Set velocity inside stationary bodies
          IF(idomy==1 .AND. mysize>1) THEN
!!!!!!            DO ibd=1,mbd-1
            DO ibd=1,mbdn-1
              call velinterior1
     &(us,vs,ws,ub,vb,wb,iuint,juint,kuint,ivint,jvint,kvint
     &,iwint,jwint,kwint,xu,yv,zw,xc,yc,zc,nx,ny,nz
     &,limuint(ibd),mimuint(ibd),limvint(ibd),mimvint(ibd),limwint(ibd),mimwint(ibd)
     &,tlevel,ibd,0,0)
            ENDDO
          ELSE
!!!!!!            DO ibd=1,mbd-1
            DO ibd=1,mbdn-1
              call velinterior
     &(us,vs,ws,ub,vb,wb,flaguo,flagvo,flagwo,xu,yv,zw
     &,xc,yc,zc,nx,ny,nz,nbd,tlevel,ibd,0,0)
            ENDDO
c            call obstacle(us,vs,ws,p,nx,ny,nz,0)
          ENDIF

        ENDIF
        clock(7) = clock(7) + tclock() - clocktemp1
c
c.....other boundary conditions
        clocktemp = tclock()
        CALL BOUNDARY(US,VS,WS,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TLEVEL)
c
C-----------------------------------------------------------------------
C.....refresh block interfaces (and periodic boundary in z direction)
C-----------------------------------------------------------------------
        CALL REFRESHBC(US,NX*NY,NZ)
        CALL REFRESHBC(VS,NX*NY,NZ)
        CALL REFRESHBC(WS,NX*NY,NZ)
        CALL REFRESHBC(UB,NX*NY,NZ)
        CALL REFRESHBC(VB,NX*NY,NZ)
        CALL REFRESHBC(WB,NX*NY,NZ)
c
c.....inlet and outlet boundary conditions
        CALL BOUNDINOUT(US,VS,WS,UO,VO,WO,XU,XC,YC,ZW,ZC,ALFXDT,TLEVEL,
     &  NX,NY,NZ,1,planeU,planeV,planeW)
	CALL BOUNDARY(US,VS,WS,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TIME)  

        clock(8) = clock(8) + tclock() - clocktemp
c
C-----------------------------------------------------------------------
c .....divergence
C-----------------------------------------------------------------------
        clocktemp = tclock()
        CALL DIVERGENCE(US,VS,WS,DP,NX,NY,NZ)
        clock(9) = clock(9) + tclock() - clocktemp

        DP = DP/ALFXDT
C
C-----------------------------------------------------------------------
C.....Poisson equation
C-----------------------------------------------------------------------

        clocktemp = tclock()
        CALL DIRECT(DP,UB,VB,DELYSQ,DELZSQ,AP,AU,CPG,CWG,RP,RU,NX,NY,NZ,
     &       NZG,MYSIZE,MYRANK,MPI_COMM_EDDY,MTYPE,IPRES,MP,LP,NP)
        clock(10) = clock(10) + tclock() - clocktemp

c
C-----------------------------------------------------------------------
C.....refresh block interfaces (and periodic boundary in z direction)
C-----------------------------------------------------------------------
        CALL REFRESHBC(DP,NX*NY,NZ)
c
C-----------------------------------------------------------------------
C.....Update velocity values , compute residuals and store RHS
c.....The corrected velocities are stored in UO,VO,WO
C-----------------------------------------------------------------------

        clocktemp = tclock()
        CALL CORRECT(US,VS,WS,UO,VO,WO,DP,NX,NY,NZ,ALFXDT)
        clock(11) = clock(11) + tclock() - clocktemp

C-----------------------------------------------------------------------
C                              boundary conditions to the updated values
C-----------------------------------------------------------------------
c.....other boundary conditions
        clocktemp = tclock()
        CALL BOUNDARY(UO,VO,WO,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TLEVEL)
c        !write(*,*), "Eddy6.f WO(nx-1,1,nz-1)
c     !&  WO(nx-1,ny-1,nz-1) 2725",WO(nx-1,1,nz-1), WO(nx-1,ny-1,nz-1)  
cc
C-----------------------------------------------------------------------
C.....refresh block interfaces (and periodic boundary in z direction)
C-----------------------------------------------------------------------
        CALL REFRESHBC(UO,NX*NY,NZ)
        CALL REFRESHBC(VO,NX*NY,NZ)
        CALL REFRESHBC(WO,NX*NY,NZ)
        

c.....inlet and outlet boundary conditions
        CALL BOUNDINOUT(UO,VO,WO,US,VS,WS,XU,XC,YC,ZW,ZC,ALFXDT,TLEVEL,
     &  NX,NY,NZ,0,planeU,planeV,planeW)
	CALL BOUNDARY(UO,VO,WO,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TIME)  

        clock(8) = clock(8) + tclock() - clocktemp
        

C-----------------------------------------------------------------------
C     Update pressure field
C-----------------------------------------------------------------------

        clocktemp = tclock()

        P = P+DP

        !Assume u*=un + dt*dp/dx, v*=vn+dt*dp/dy, w*=wn+dt*dp/dz
        IF(IMPLY==1) THEN

c If there is an implicit treatment another correction on the pressure
c field is necessary
          DO I=2,IX2
            P(I,2:JY2,2:KZ2) = P(I,2:JY2,2:KZ2)
     &           -0.5*ALFXDT*RU1/(DELY*RP(I))**2*
     &           (DP(I,1:NY-2,2:KZ2)-2.*DP(I,2:JY2,2:KZ2)+DP(I,3:NY,2:KZ2))*RVIMPLY(I)
          END DO

        ENDIF

        IF(IMPLX==1) THEN

          DO I=2,IX2
            P(I,2:JY2,2:KZ2) = P(I,2:JY2,2:KZ2)
     &           -0.5*ALFXDT*RU1*AP(I)/RP(I)*
     &           ((DP(I+1,2:JY2,2:KZ2)-DP(I  ,2:JY2,2:KZ2))*AU(I  )*RU(I  )
     &           -(DP(I  ,2:JY2,2:KZ2)-DP(I-1,2:JY2,2:KZ2))*AU(I-1)*RU(I-1))*RVIMPLX(I)
          END DO
c
        ENDIF

c A reference value for the pressure is established: this is sent
c from the processor 0 to all processors !!!!!!
c All processors scale the pressure value using PREF
!!!!!!        pref = p(nx/2,ny/2,2)
        pref = sum(p(ix1,jy1:jy2,2))/real(ny-2)

        CALL MPI_BCAST(PREF,1,MTYPE,0,MPI_COMM_EDDY,IERR)

        P = P-PREF
c        write(*,*), "Eddy6.f u(nx-1,1,nz-1)
c     &  u(nx-1,ny-1,nz-1) 2782",WO(nx-1,1,nz-1), WO(nx-1,ny-1,nz-1)  
c
C
C-----------------------------------------------------------------------
c     Cyclic b.c. in x direction
C-----------------------------------------------------------------------
        IF(MP==0) THEN
          P(1 ,:,:) = P(IX2,:,:)
          P(NX,:,:) = P(IX1,:,:)
        ELSE
          P(1 ,:,:) = -P(IX1,:,:)
          P(NX,:,:) = -P(IX2,:,:)
        ENDIF
	IF(ITYPE(1)==300) THEN
          DO J=1,NY
            P(1,J,:) =  P(IX1,JSYM(J),:)
          ENDDO
        ENDIF

C-----------------------------------------------------------------------
c     Cyclic b.c. in y direction
C-----------------------------------------------------------------------
        IF(LP==0) THEN
          P(:,1 ,:) = P(:,JY2,:)
          P(:,NY,:) = P(:,JY1,:)
        ELSE
          P(:,1 ,:) = P(:,JY1,:)
          P(:,NY,:) = P(:,JY2,:)
        ENDIF
c        write(*,*), "Eddy6.f u(nx-1,1,nz-1)
c     &  u(nx-1,ny-1,nz-1) 2808",WO(nx-1,1,nz-1), WO(nx-1,ny-1,nz-1)  
cC-----------------------------------------------------------------------
c     Neumann b.c. in z direction
c N.B. ITYPE(5)/=0 only for the process 0
c N.B. ITYPE(6)/=0 only for the process MYSIZE-1
C-----------------------------------------------------------------------
	IF(ITYPE(5)==500) P(:,:,1 ) = P(:,:,KZ2)
        IF(ITYPE(6)==500) P(:,:,NZ) = P(:,:,KZ1)

        IF(ITYPE(5)/=0.AND.ITYPE(5)/=500) P(:,:,1 ) = P(:,:,KZ1)
        IF(ITYPE(6)/=0.AND.ITYPE(6)/=500) P(:,:,NZ) = P(:,:,KZ2)
! 	CALL BOUNDARY_P(P,XC,NX,NY,NZ)
C-----------------------------------------------------------------------
C.....refresh block interfaces (and periodic boundary in z direction)
C-----------------------------------------------------------------------
        CALL REFRESHBC(P,NX*NY,NZ)

C
C-----------------------------------------------------------------------
c     boundary pressure and interior points
C-----------------------------------------------------------------------
        clock(12) = clock(12) + tclock()-clocktemp

        clocktemp = tclock()

        IF(ibm>=1) then
c
c.....field extension --- pressure
c
c
ccc          do ibd=1,nbd
c the pressure at the interface points is corrected using the normal
c pressure gradient
c for the circumferential moving bodies the normal pressure gradient is 0 !!!!!!
!!!!!!            call correctpres(p,nx,ny,nz,ip,jp,kp,pmtrx,pim,nxp,nyp,nzp
!!!!!!     &            ,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect,fp,nfacet
!!!!!!     &            ,xc_car,yc_car,xc,yc,zc,limp(ibd,1),mimp(ibd,1),nflpmax,tlevel,alfxdt,ibd,mbd)
ccc            ilb = lb(ibd)+1
ccc            ile = lb(ibd)+mb(ibd)
ccc            call correctpres(p,nx,ny,nz,ip,jp,kp,pmtrx,pim,nxp,nyp,nzp
ccc     &,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect(:,ilb:ile),fp,mb(ibd)
ccc     &,xc_car,yc_car,xc,yc,zc,limp(ibd,:),mimp(ibd,:),nflpmax(ibd)
ccc     &,tlevel,alfxdt,ibd,mbd)
          call correctpres_mod2(p,nx,ny,nz,ip,jp,kp,pmtrx,pim,nxp,nyp,nzp
     &,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect,fp,nfacet
     &,xc_car,yc_car,xc,yc,zc,limp,mimp,nflpmax
     &,tlevel,alfxdt,mbd,nbd,sum(mimp))

!           call correctpres_mod2(dens,nx,ny,nz,ip,jp,kp,pmtrx,pim,nxp,nyp,nzp
!      &,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect,fp,nfacet
!      &,xc_car,yc_car,xc,yc,zc,limp,mimp,nflpmax
!      &,tlevel,alfxdt,mbd,nbd,sum(mimp))
c the pressure at the body points is set equal to 0
          do ibd=1,nbd
            if(idomy==1 .AND. mysize>1) then
              call presinterior1(p,nx,ny,nz,ipint,jpint,kpint,limpint(ibd),mimpint(ibd))
            else
              call presinterior(p,nx,ny,nz,flagpo,ibd,nbd)
!               call presinterior(dens,nx,ny,nz,flagpo,ibd,nbd)
c              call obstacle(us,vs,ws,p,nx,ny,nz,1)
            endif

          enddo

        ENDIF
        clock(13) = clock(13) + tclock() - clocktemp
c        write(*,*), "Eddy6.f u(nx-1,1,nz-1)
c     &  u(nx-1,ny-1,nz-1) 2873",WO(nx-1,1,nz-1), WO(nx-1,ny-1,nz-1)  
cC
C-----------------------------------------------------------------------
C.....refresh block interfaces (and periodic boundary in z direction)
C-----------------------------------------------------------------------
        clocktemp = tclock()
        CALL REFRESHBC(P,NX*NY,NZ)
C
C-----------------------------------------------------------------------
c     Neumann b.c. in z direction
C-----------------------------------------------------------------------
C

        IF(ITYPE(5)/=0.AND.ITYPE(5)/=500) P(:,:,1 ) = P(:,:,KZ1)
        IF(ITYPE(6)/=0.AND.ITYPE(6)/=500) P(:,:,NZ) = P(:,:,KZ2)
C
c
c.....axis
C
        IF(ITYPE(1)==300) THEN
          DO J=1,NY
            P(1,J,:) =  P(IX1,JSYM(J),:)
          ENDDO
        ENDIF
C-----------------------------------------------------------------------
c     Cyclic b.c. in x direction
C-----------------------------------------------------------------------
        IF(MP==0) THEN
          P(1 ,:,:) = P(IX2,:,:)
          P(NX,:,:) = P(IX1,:,:)
        ELSE
          P(1 ,:,:) = P(IX1,:,:)
          P(NX,:,:) = P(IX2,:,:)
        ENDIF

C-----------------------------------------------------------------------
c     Cyclic b.c. in y direction
C-----------------------------------------------------------------------

        IF(LP==0) THEN
          P(:,1 ,:) = P(:,JY2,:)
          P(:,NY,:) = P(:,JY1,:)
        ELSE
          P(:,1 ,:) = P(:,JY1,:)
          P(:,NY,:) = P(:,JY2,:)
        ENDIF
        clock(12) = clock(12) + tclock() - clocktemp
C
C Evaluates the density field based on the Boussinesq approximation
C
!         IF(IDENS.EQ.1) THEN
         
          CALL BOUNDINOUTD(DENS,WO,ZCG,ALFXDT,NX,NY,NZ,NZG,planeD)
          CALL BOUNDARY_DENS(DENS,XC,YC,NX,NY,NZ)
          CALL RHS_DENSITY(UO,VO,WO,DENS,TV,RHB,RS,NX,NY,NZ,YC,YV)
          CALL BOUNDARY_DENS(RHB,XC,YC,NX,NY,NZ)
          CALL BOUNDARY_DENS(RS,XC,YC,NX,NY,NZ)
	  CALL REFRESHBC(RHB,NX*NY,NZ)
          CALL REFRESHBC(RS,NX*NY,NZ)
         
          CALL DENSITY(DENS,UO,VO,TV,RHA,RHB,RS,XC,YC,ALFXDT,GAMXDT,RHOXDT,NX,NY,NZ)
	 
          CALL BOUNDINOUTD(DENS,WO,ZCG,ALFXDT,NX,NY,NZ,NZG,planeD)
          CALL BOUNDARY_DENS(DENS,XC,YC,NX,NY,NZ)
          CALL REFRESHBC(DENS,NX*NY,NZ)

          IF(ibm>=1) then
          call correctdens_mod2(dens,nx,ny,nz,ip,jp,kp,pmtrx,rhoim,nxp,nyp,nzp
     &,xnp,ynp,znp,dirp,mrkp,pindx,drhodnn,unvect,fp,nfacet
     &,xc_car,yc_car,xc,yc,zc,limp,mimp,nflpmax
     &,tlevel,alfxdt,mbd,nbd,sum(mimp))

            do ibd=1,nbd
              if(idomy==1 .AND. mysize>1) then
              call densinterior1(dens,nx,ny,nz,ipint,jpint,kpint,limpint(ibd),mimpint(ibd))
              else
              call densinterior(dens,nx,ny,nz,flagpo,ibd,nbd)
              endif
            enddo
            ENDIF

	  IF(MYRANK.EQ.0) then 

          WRITE(*,*)"**********END OF RK substep**************",is
          call sequence_hybrid(planeU,planeV,planeW,planeD,
     &    planeURef,planeVRef,planeWRef,planeDRef,
     &    icycle,nx,ny)
          endif
C
C-----------------------------------------------------------------------
C                                                 time marching loop end
C-----------------------------------------------------------------------
      enddo
c
      TIME = TLEVEL
      write(49,'(2E25.12)')time,dtm1
c        write(*,*), "Eddy6.f u(nx-1,1,nz-1)
c     &  u(nx-1,ny-1,nz-1) 2958",WO(nx-1,1,nz-1), WO(nx-1,ny-1,nz-1)  

cc
c-----------------------------------------------------------------------
c                                            compute turbulent viscosity
c-----------------------------------------------------------------------
c
      clocktemp = tclock()
      IF(ISGS/=0) THEN
c
        IF(ISGS==5) THEN
          CALL STRUCTUREFUNCTION(UO,VO,WO,TV,DP,NX,NY,NZ)
        ELSEIF(ISGS==6) THEN
          CALL WALE(UO,VO,WO,TV,NX,NY,NZ)
        ELSE
          CALL TURVIS(UO,VO,WO,TV,G,DP,SXX,SYY,SZZ,SXY,SYZ,SXZ,
     &         US,VS,WS,UB,VB,WB,ILM,IMM,LM,MM,UC,VC,WC,
     &         DXDYDZ,CLES,CLESP,CLESN,NILMP,NILMN,XC,YC,ZC,DTM1,ICYCLE,NX,NY,NZ)


        ENDIF

        if(ibm/=0) then
c the eddy viscosity is evaluated at the interface points

ccc          do ibd=1,nbd
!!!!!!            call momforc(tv,tv,xc,yc,zc,nx,ny,nz,itv,jtv,ktv,tvmtrx,tvindx
!!!!!!     &            ,tvim,mrktv,dirtv,nxtv,nytv,nztv,nfltvmax,limtv(ibd,1),mimtv(ibd,1),0,3)
ccc            call momforc(tv,tv,xc,yc,zc,nx,ny,nz,itv,jtv,ktv,tvmtrx,tvindx
ccc     &,tvim,mrktv,dirtv,nxtv,nytv,nztv,nfltvmax(ibd),limtv(ibd,:),mimtv(ibd,:),0,3)
            call momforc_mod2(tv,tv,xc,yc,zc,nx,ny,nz,itv,jtv,ktv,tvmtrx,tvindx
     &,tvim,mrktv,dirtv,nxtv,nytv,nztv,nfltvmax,limtv,mimtv,0,3,nbd,sum(mimtv))
            tv(:,1,:) = tv(:,ny-1,:)
            tv(:,ny,:) = tv(:,2,:)
c
c......Set turbulent viscosity inside bodies
          do ibd=1,nbd
            IF(idomy==1 .AND. mysize>1) THEN
              call tvinterior1(tv,itv,jtv,ktv,xc,yc,zc,nx,ny,nz,limtvint(ibd),mimtvint(ibd))
            ELSE
              call tvinterior(tv,flagtv,nx,ny,nz,nbd,ibd)
            ENDIF
          enddo

        endif

      ENDIF
      clock(14) = clock(14) + tclock() - clocktemp
C
C-----------------------------------------------------------------------
c                                                pressure on the bodies
C-----------------------------------------------------------------------
      clocktemp = tclock()

      IF(IBM>0) iflag_psb=0
      IF(IBM>1) THEN
        iflag_ptr = 0
        iflag_pnd = 0
      ENDIF

      IF(IBM>0 .AND. ITBDY/=0.AND. MOD(ICYCLE,ITBDY)==0) THEN

c        CALL SHEARIMB(flagpo,flag,nx,ny,nz,nbd,xu,xc,yv,yc,zw,zc,zwg,zcg,nzg
c     &        ,xc_car,yc_car,vertexc,unvect,nfacet,icycle)

        IF(iflag_psb==0) THEN
          iflag_psb=1

        nimp = sum(mimp(:,:))   !!!!!!
        if(ivrtx.eq.0) then
          mimpbd = 0
        else
          mimpnd = 0
        endif
        DO IBD=1,NBD

          IF(IVRTX==0) THEN
            ilb = lb(ibd)+1
            ile = lb(ibd)+mb(ibd)
          ELSE
            ilb = lv(ibd)+1
            ile = lv(ibd)+mv(ibd)
          ENDIF

!!!!!!          IF(ICALF==0) THEN
!!!!!!            icom = 2
!!!!!!          ELSE
!!!!!!            icom = 2
!!!!!!            IF(MOD(ICYCLE,ITCALF)==0 .AND. IVRTX==1) icom=1
!!!!!!          ENDIF
          icom=-1

          IF(iflag_ptr==0 .AND. ivrtx==0) THEN
!!!!!!            call TAGPBD2(FLAGPO(:,:,:,ibd),FLAG,XC,YC,XC_CAR,YC_CAR
!!!!!!     &            ,ZC,ZCG,NX,NY,NZ,NZG,IBD,VERTEXC,UNVECT,NFACET)
            call TAGPBD2(FLAGPO(:,:,:,ibd),FLAG,XC,YC,XC_CAR,YC_CAR
     &,ZC,ZCG,NX,NY,NZ,NZG,IBD,VERTEXC(:,ilb:ile),UNVECT(:,ilb:ile),mb(ibd))
c FLAGPBD finds the points of the computational grid where it is
c necessary to evaluate the physical pressure
            flagpbd(:,:,:,ibd) = flag(:,:,:)
!!!!!!            nimp = sum(mimp)
!!!!!!            call geom(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs),vertex,vertexc,unvect
!!!!!!     &,limpbd(ibd,1),mimpbd(ibd,1),ip,jp,kp,nimp,xnp,ynp,znp,nxp,nyp,nzp
!!!!!!     &,mrkp,dirp,fp,flag,flagpbd,nx,nyl,nzl,nbd,nfacet,ibd,icom,nflp,nflpbdmax,4,icycle,tlevel)
ccc            call geom(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs)
            call geom_mod(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs)
     &,vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile)
     &,limpbd(ibd,:),mimpbd(ibd,:),ip,jp,kp,nimp,xnp,ynp,znp,nxp,nyp,nzp
     &,mrkp,dirp,fp,flag,flagpbd,nx,nyl,nzl,nbd,mb(ibd),ibd,icom,nflp
     &,nflpbdmax(ibd),4,icycle,tlevel)
            call mtrx(xc,yc(jys),zcg(kzs),limpbd(ibd,1),sum(mimpbd(ibd,:)),mrkp
     &,ip,jp,kp,xnp,ynp,znp,nxp,nyp,nzp,pmtrx,pindx,nx,nyl,nzl,nbd,ibd,0)
            iflag_ptr = 1
          ENDIF

          IF(iflag_pnd==0 .AND. ivrtx==1) THEN
!!!!!!            call TAGPBD2(FLAGPO(:,:,:,ibd),FLAG,XC,YC,XC_CAR,YC_CAR
!!!!!!     &            ,ZC,ZCG,NX,NY,NZ,NZG,IBD,TRVTX,VTXNRM,NVTX)
            call TAGPBD2(FLAGPO(:,:,:,ibd),FLAG,XC,YC,XC_CAR,YC_CAR
     &,ZC,ZCG,NX,NY,NZ,NZG,IBD,TRVTX(:,ilb:ile),VTXNRM(:,ilb:ile),mv(ibd))
            flagpbd(:,:,:,ibd) = flag(:,:,:)
!!!!!!            nimp = sum(mimp)
!!!!!!            call geom(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs),vertex,vertexc,unvect
!!!!!!     &,limpnd(ibd,1),mimpnd(ibd,1),ip,jp,kp,nimp,xnp,ynp,znp,nxp,nyp,nzp
!!!!!!     &,mrkp,dirp,fp,flag,flagpbd,nx,nyl,nzl,nbd,nfacet,ibd,icom,nflp,nflpbdmax,4,icycle,tlevel)
            ii = lb(ibd)+1
            jj = lb(ibd)+mb(ibd)
ccc            call geom(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs)
            call geom_mod(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs)
     &,vertex(:,:,ii:jj),vertexc(:,ii:jj),unvect(:,ii:jj)
     &,limpnd(ibd,:),mimpnd(ibd,:),ip,jp,kp,nimp,xnp,ynp,znp,nxp,nyp,nzp
     &,mrkp,dirp,fp,flag,flagpbd,nx,nyl,nzl,nbd,mb(ibd),ibd,icom,nflp
     &,nflpbdmax(ibd),4,icycle,tlevel)
            call mtrx(xc,yc(jys),zcg(kzs),limpnd(ibd,1),sum(mimpnd(ibd,:)),mrkp
     &,ip,jp,kp,xnp,ynp,znp,nxp,nyp,nzp,pmtrx,pindx,nx,nyl,nzl,nbd,ibd,0)
            iflag_pnd=1
          ENDIF
        ENDDO

        dp = p

        call correctpres_mod2(dp,nx,ny,nz,ip,jp,kp,pmtrx,pim,nxp,nyp,nzp
     &,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect,fp,nfacet
     &,xc_car,yc_car,xc,yc,zc,limpbd,mimpbd,nflpbdmax
     &,tlevel,alfxdt,mbd,nbd,sum(mimpbd))

c the pressure in the points found by the TAGPBD2 is corrected using the
c normal gradient
        IF(IVRTX==0) THEN
!!!!!!            call correctpres(dp,nx,ny,nz,ip,jp,kp,pmtrx,pim,nxp,nyp,nzp
!!!!!!     &,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect,fp,nfacet
!!!!!!     &,xc_car,yc_car,xc,yc,zc,limpbd(ibd,1),mimpbd(ibd,1),nflpbdmax
!!!!!!     &,tlevel,alfxdt,ibd,mbd)
ccc            call correctpres(dp,nx,ny,nz,ip,jp,kp,pmtrx,pim,nxp,nyp,nzp
ccc     &,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect(:,ilb:ile),fp,mb(ibd)
ccc     &,xc_car,yc_car,xc,yc,zc,limpbd(ibd,:),mimpbd(ibd,:),nflpbdmax(ibd)
ccc     &,tlevel,alfxdt,ibd,mbd)
c PRESS_BODY3 interpolates the pressure at the center of the triangles
c using the surrounding grid nodes
c output:PBD,MRKPB
          DO IBD=1,NBD
            ilb = lb(ibd)+1
            ile = lb(ibd)+mb(ibd)
            call press_body3(dp,flagpbd,xc,yc,zc,xc_car,yc_car,nx,ny,nz,nbd,pbd,mrkpb
     &          ,vertexc,unvect,nfacet,ibd,mbd,ilb,ile,tlevel,alfxdt)
          ENDDO
        ELSE
!!!!!!            call correctpres(dp,nx,ny,nz,ip,jp,kp,pmtrx,pim,nxp,nyp,nzp
!!!!!!     &,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect,fp,nfacet
!!!!!!     &,xc_car,yc_car,xc,yc,zc,limpnd(ibd,1),mimpnd(ibd,1),nflpbdmax
!!!!!!     &,tlevel,alfxdt,ibd,mbd)
ccc            call correctpres(dp,nx,ny,nz,ip,jp,kp,pmtrx,pim,nxp,nyp,nzp
ccc     &,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect(:,ii:jj),fp,mb(ibd)
ccc     &,xc_car,yc_car,xc,yc,zc,limpnd(ibd,:),mimpnd(ibd,:),nflpbdmax(ibd)
ccc     &,tlevel,alfxdt,ibd,mbd)
          DO IBD=1,NBD
            ilb = lv(ibd)+1
            ile = lv(ibd)+mv(ibd)
            call press_body3(dp,flagpbd,xc,yc,zc,xc_car,yc_car,nx,ny,nz
     &            ,nbd,pbd,mrkpb,trvtx,vtxnrm,nvtx,ibd,mbd,ilb,ile,tlevel,alfxdt)
          ENDDO
        ENDIF

        ENDIF

        ub = uo
        vb = vo
        wb = wo

c.....compute forcing in x momentum
        call momforc_mod2(ub,ub,xu,yc,zc,nx,ny,nz,iu,ju,ku,umtrx,uindx,uim
     &         ,mrku,diru,nxu,nyu,nzu,nflumax,limu,mimu,0,1,nbd,sum(mimu))
c.....compute forcing in y momentum
        call momforc_mod2(vb,vb,xc,yv,zc,nx,ny,nz,iv,jv,kv,vmtrx,vindx,vim
     &         ,mrkv,dirv,nxv,nyv,nzv,nflvmax,limv,mimv,0,2,nbd,sum(mimv))
c.....compute forcing in z momentum
        call momforc_mod2(wb,wb,xc,yc,zw,nx,ny,nz,iw,jw,kw,wmtrx,windx,wim
     &         ,mrkw,dirw,nxw,nyw,nzw,nflwmax,limw,mimw,0,3,nbd,sum(mimw))

        ub(:,1,:) = ub(:,ny-1,:)
        ub(:,ny,:) = ub(:,2,:)
        vb(:,1,:) = vb(:,ny-1,:)
        vb(:,ny,:) = vb(:,2,:)
        wb(:,1,:) = wb(:,ny-1,:)
        wb(:,ny,:) = wb(:,2,:)

c the SHEAR is evaluated by SHEAR_BODY using the velocities UB,VB,WB
c established by MOMFORC
c output: components along X,Y,Z of the normal derivatives of U,V,W
c (DUDXB,DUDYB,DUDZB,DVDXB,DVDYB,DVDZB,DWDXB,DWDYB,DWDZB)





        IF(IVRTX==0) THEN
          DO ibd=1,mbdn-1
            ilb = lb(ibd)+1
            ile = lb(ibd)+mb(ibd)
            CALL SHEAR_BODY(ub,vb,wb,nx,ny,nz,xu,yv,zw,xc,yc,zc
     &            ,xu_car,yu_car,xv_car,yv_car,xc_car,yc_car,vertexc,unvect,nfacet
     &            ,dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb,cf_aux,mrksb,ibd,nbd,ilb,ile,tlevel,0)
          ENDDO
          DO ibd=mbdn,nbd
            ilb = lb(ibd)+1
            ile = lb(ibd)+mb(ibd)
            CALL SHEAR_BODY(ub,vb,wb,nx,ny,nz,xu,yv,zw,xc,yc,zc
     &            ,xu_car,yu_car,xv_car,yv_car,xc_car,yc_car,vertexc,unvect,nfacet
     &            ,dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb,cf_aux,mrksb,ibd,nbd,ilb,ile,tlevel,1)
          ENDDO
        ELSE
          DO ibd=1,mbdn-1
            ilb = lv(ibd)+1
            ile = lv(ibd)+mv(ibd)
            CALL SHEAR_BODY(ub,vb,wb,nx,ny,nz,xu,yv,zw,xc,yc,zc
     &            ,xu_car,yu_car,xv_car,yv_car,xc_car,yc_car,trvtx,vtxnrm,nvtx
     &            ,dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb,cf_aux,mrksb,ibd,nbd,ilb,ile,tlevel,0)
          ENDDO
          DO ibd=mbdn,nbd
            ilb = lv(ibd)+1
            ile = lv(ibd)+mv(ibd)
            CALL SHEAR_BODY(ub,vb,wb,nx,ny,nz,xu,yv,zw,xc,yc,zc
     &            ,xu_car,yu_car,xv_car,yv_car,xc_car,yc_car,trvtx,vtxnrm,nvtx
     &            ,dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb,cf_aux,mrksb,ibd,nbd,ilb,ile,tlevel,1)
          ENDDO
        ENDIF

        nwritebdy = nwritebdy+1

c The values of pressure and shear stress are written on a file




        DO ibd=1,nbd
          IF(IVRTX==0) THEN
            ilb = lb(ibd)+1
            ile = lb(ibd)+mb(ibd)


c$$$            CALL IOSHEARVERTEXC('imb'//index(ibd)//'.'//index(nwritebdy)
c$$$     &         //'.tno',nbd,trino,pbd,dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb
c$$$     &         ,dwdxb,dwdyb,dwdzb,mrkpb,mrksb,nfacet,nfacetot,ilb,ile)

            write(icycle_string,'(i8.8)') icycle

            CALL IOSHEARVERTEXC('imb'//index(ibd)//'.'//trim(icycle_string)
     &         //'.tno',nbd,trino,pbd,dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb
     &         ,dwdxb,dwdyb,dwdzb,cf_aux,mrkpb,mrksb,nfacet,nfacetot,ilb,ile)



          ELSE
!!!!!            ilb = lb(ibd)+1
!!!!!            ile = lb(ibd)+mv(ibd)
            ilb = lv(ibd)+1
            ile = lv(ibd)+mv(ibd)
            CALL IOSHEARVERTEXC('imb.'//index(ibd)//'.'//index(nwritebdy)
     &         //'.vno',nbd,vtxno,pbd,dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb
     &         ,dwdxb,dwdyb,dwdzb,cf_aux,mrkpb,mrksb,nvtx,nvtxtot,ilb,ile)
          ENDIF
        ENDDO
!        if(myrank.eq.0) then
!           open(unit=10,file='time.imb',form='formatted',position
!     $          ='append')
!           write(10,*) nwritebdy,tlevel
!           close(10)
!        endif

c        call calc_cf(dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb
c     &     ,mrksb,cf,mrkcf,unvect,nfacet)
!        if(myrank==0) write(*,*)'******CALCULATE Cf******'
!        allocate(cfg(nfacet),cf(nfacet))
!        cf = 0.0
!        do i=1,nfacet

!        unorm = vecmag(unvect(:,i),3)
!        a1 = unvect(1,i)/unorm
!        a2 = unvect(2,i)/unorm
!        a3 = unvect(3,i)/unorm

!        sxxb = 2.*dudxb(i)
!        syyb = 2.*dvdyb(i)
!        szzb = 2.*dwdzb(i)
!        sxyb = dudyb(i)+dvdxb(i)
!        sxzb = dudzb(i)+dwdxb(i)
!        syzb = dvdzb(i)+dwdyb(i)

c        if(sum(mrks(i,1:6))==6) then
c          mrkcf(i) = 1
!          cf(i) = ru1*(sxzb*a1 + syzb*a2 + szzb*a3)
!           cf(i) = (ru1*(dwdxb(i)*a1+dwdyb(i)*a2+dwdzb(i)*a3))
c        endif

!         enddo

!         CALL MPI_REDUCE(cf,cfg,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)

!         call calc_cf(dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb
!      &     ,mrksb,cf,mrkcf,unvect,nfacet)
!        write(ICfile,'(a,i8.8,a)')trim(hspdir)//"cf_",icycle,".plt"
!       call iotaubd(ICfile,vertexc,dwdzb,mrkpb,nfacet,nfacetmax,nfacetot,nbd)
!       write(ICfile,'(a,i8.8,a)')trim(hspdir)//"pbd_",icycle,".plt"
!       call iotaubd(ICfile,vertexc,pbd,mrkpb,nfacet,nfacetmax,nfacetot,nbd)
!
!        deallocate(cfg,cf)
      ENDIF

      clock(15) = clock(15) + tclock() - clocktemp
C
C-----------------------------------------------------------------------
c.... write restart file
C-----------------------------------------------------------------------
C
!      utmp1=utmp1+uo
!      vtmp1=vtmp1+vo
!      wtmp1=wtmp1+wo
!      ptmp1=ptmp1+p
!      denstmp1=denstmp1+dens
!      ttotal=ttotal+dtm1
      clocktemp = tclock()

      IF(MOD(ICYCLE,ITPOST)==0.OR.MOD(ICYCLE,ITPOST)==0) THEN
        if(it5p == 1) then
          nstep=icycle
          IFIELD = 2
          IDIR  = 1  ! IOMPI_3DSCALAR works in writing mode

          write(ICfile,'(a,i8.8,a)')trim(hspdir)//"up_",icycle,".5p"
          CALL IOSCALAR_POST_5P(ICfile,UO,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
          if(myrank==0) write(6,*) "up.5p done"

          write(ICfile,'(a,i8.8,a)')trim(hspdir)//"vp_",icycle,".5p"
          CALL IOSCALAR_POST_5P(ICfile,VO,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
          if(myrank==0) write(6,*) "vp.5p done"

          write(ICfile,'(a,i8.8,a)')trim(hspdir)//"wp_",icycle,".5p"
          CALL IOSCALAR_POST_5P(ICfile,WO,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
          if(myrank==0) write(6,*) "wp.5p done"

          write(ICfile,'(a,i8.8,a)')trim(hspdir)//"pp_",icycle,".5p"
          CALL IOSCALAR_POST_5P(ICfile,P ,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
          if(myrank==0) write(6,*) "pp.5p done"

          IF(IDENS.EQ.1) then
          write(ICfile,'(a,i8.8,a)')trim(hspdir)//"densp_",icycle,".5p"
          CALL IOSCALAR_POST_5P(ICfile,DENS,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
          if(myrank==0) write(6,*) "densp.5p done"
          endif
        endif
      ENDIF

      IF(ITRES/=0.AND.MOD(ICYCLE,ITRES)==0 .OR. SAVE_RES) THEN

        IF(MYRANK==0) THEN
          write(6,*) ' Begin write restarting file ...',icycle,time
          open(unit=10,file='time.res',form='formatted')
          write(10,*) time
          close(10)
        ENDIF
        nstep=icycle
        IFIELD = 2
        IDIR  = 1  ! IOMPI_3DSCALAR works in writing mode
c        CALL IOMPI_3DSCALAR(trim(hspdir)//'u.res',UO,NX,NY,NZ,IDIR)
c        CALL IOMPI_3DSCALAR(trim(hspdir)//'v.res',VO,NX,NY,NZ,IDIR)
c        CALL IOMPI_3DSCALAR(trim(hspdir)//'w.res',WO,NX,NY,NZ,IDIR)
c        CALL IOMPI_3DSCALAR(trim(hspdir)//'p.res', P,NX,NY,NZ,IDIR)

        write(ICfile,'(a,i8.8,a)')trim(hspdir)//"u_",icycle,".res"
        CALL IOSCALAR(ICfile,UO,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
!         CALL HDF5_MPI_3DREAL(ICfile,UO,NX,NY,NZ,IDIR)
        write(ICfile,'(a,i8.8,a)')trim(hspdir)//"v_",icycle,".res"
        CALL IOSCALAR(ICfile,VO,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
!         CALL HDF5_MPI_3DREAL(ICfile,VO,NX,NY,NZ,IDIR)
        write(ICfile,'(a,i8.8,a)')trim(hspdir)//"w_",icycle,".res"
        CALL IOSCALAR(ICfile,WO,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
!         CALL HDF5_MPI_3DREAL(ICfile,WO,NX,NY,NZ,IDIR)
        write(ICfile,'(a,i8.8,a)')trim(hspdir)//"p_",icycle,".res"
        CALL IOSCALAR(ICfile,P ,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
!         CALL HDF5_MPI_3DREAL(ICfile, P,NX,NY,NZ,IDIR)
        write(ICfile,'(a,i8.8,a)')trim(hspdir)//"tv_",icycle,".res"
        CALL IOSCALAR(ICfile,TV ,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
!         CALL HDF5_MPI_3DREAL(ICfile, TV,NX,NY,NZ,IDIR)

        IF(IDENS.EQ.1) then
        write(ICfile,'(a,i8.8,a)')trim(hspdir)//"dens_",icycle,".res"
        CALL IOSCALAR(ICfile,DENS,DP,NX,NY,NZ,IDIR,TIME,DTM1,nstep)
!         CALL HDF5_MPI_3DREAL(ICfile,DENS,NX,NY,NZ,IDIR)
        endif

        SAVE_RES=.FALSE.

c        call write_var_arc(p,xc,yc,zc,nx,ny,nz)
c        DP(:,1,:) = SUM(WO(:,2:NY-1,:),2)/REAL(NY-2)
c        CALL IOSCALAR('w2d.res',DP(:,1,:),DP(:,2,:),NX,1,NZ,IDIR,TIME)


        IF(ISCHM==1) THEN !in the case AB is used
c UA,VA,WA are the  components at the time level L-2
c          CALL IOMPI_3DSCALAR(trim(hspdir)//'ua.res',UA,NX,NY,NZ,IDIR)
c          CALL IOMPI_3DSCALAR(trim(hspdir)//'va.res',VA,NX,NY,NZ,IDIR)
c          CALL IOMPI_3DSCALAR(trim(hspdir)//'wa.res',WA,NX,NY,NZ,IDIR)
!          CALL HDF5_MPI_3DREAL(trim(hspdir)//'ua.res',UA,NX,NY,NZ,IDIR)
!          CALL HDF5_MPI_3DREAL(trim(hspdir)//'va.res',VA,NX,NY,NZ,IDIR)
!          CALL HDF5_MPI_3DREAL(trim(hspdir)//'wa.res',WA,NX,NY,NZ,IDIR)
c          CALL IOSCALAR('ua.res',UA,DP,NX,NY,NZ,IDIR,TIME)
c          CALL IOSCALAR('va.res',VA,DP,NX,NY,NZ,IDIR,TIME)
c          CALL IOSCALAR('wa.res',WA,DP,NX,NY,NZ,IDIR,TIME)
!          IF(IDENS.EQ.1)
!          CALL HDF5_MPI_3DREAL(trim(hspdir)//'rha.res',RHA,NX,NY,NZ,IDIR)
        ENDIF
        IF(ISGS/=0) THEN !if a SGS model is used

!c          CALL IOMPI_3DSCALAR(trim(hspdir)//'tv.res',TV,NX,NY,NZ,IDIR)
!          CALL HDF5_MPI_3DREAL(trim(hspdir)//'tv.res',TV,NX,NY,NZ,IDIR)
!c          CALL IOSCALAR('tv.res',TV,DP,NX,NY,NZ,IDIR,TIME)

          IF(ISGS==4) THEN !if a Lagrangian model is used
c            CALL IOMPI_3DSCALAR(trim(hspdir)//'ilm.res',ILM,NX,NY,NZ,IDIR)
c            CALL IOMPI_3DSCALAR(trim(hspdir)//'imm.res',IMM,NX,NY,NZ,IDIR)
            !CALL HDF5_MPI_3DREAL(trim(hspdir)//'ilm.res',ILM,NX,NY,NZ,IDIR)
            !CALL HDF5_MPI_3DREAL(trim(hspdir)//'imm.res',IMM,NX,NY,NZ,IDIR)
c            CALL IOSCALAR('ilm.res',ILM,DP,NX,NY,NZ,IDIR,TIME)
c            CALL IOSCALAR('imm.res',IMM,DP,NX,NY,NZ,IDIR,TIME)
          ENDIF
        ENDIF

c in the case of moving boundaries the STL files are written
        if(ibm>1) then
             do ibd=mbd,nbd
c          CALL IOMPI_IMB(UNVECT,VERTEX,VERTEXC,AREAF,TRINO,NFACET,ZC,NZ,NBD,MBD,MBG,1)
                CALL WRITESTL('mbd.'//index(ibd)//'.res',unvect,vertex,vertexc,zc,nz,nfacet,nfacetot,ibd,0)
                IF(ITBDY>0 .AND. IVRTX==1) CALL WRITETRVTX
     %('mbd_nodes_conn.'//index(ibd)//'.res',TRVTX,VTXNRM,NVTX,IBD)   !!!!!! IBD
c            call imb2plt('mbd.'//index(ibd)//'.plt.res',vertex,nfacet,ibd)
             enddo
        endif

        if(myrank==0) write(6,*) ' End write restarting file'
c
      ENDIF

      !Flat surface dimple, shear stress
      IF(icfprb>0) THEN
c         write(6,*) 'icycle=',icycle,', icfprb=',icfprb,', itcfprb=',itcfprb
         IF(MOD(ICYCLE,ITCFPRB)==0) THEN
            CALL RECORD_CF_FLATSURF(wo,xc,yc,zw,nx,ny,nz,jcf,kcf,njcf
     $           ,nkcf,nkcfmax,cfprb,tcf,ncfprb,ncfsmpl,tlevel)
         ENDIF
         IF(MOD(ICYCLE,ITRES)==0) THEN
           CALL WRITE_CF_FLATSURF(filecfprb,cfprb,tcf,njcf,nkcf,nkcfmax,nkcfprev,ncfprb,ncfsmpl)
         ENDIF
      ENDIF

      if(flag_wavg_top>0 .AND. mod(icycle,freq_wavg_top)==0) then
        call calc_wavg_top(wavg_top,wo,nx,ny,nz,wavg_top_nsmpl)
        if(mod(icycle,itres)==0) then
          call write_wavg_top(file_wavg_top,wavg_top,ny,nz,wavg_top_nsmpl)
        endif
      endif


      IF(IT2D/=0 .AND. MOD(ICYCLE,IT2D)==0) THEN
c Plot of the Y vorticity fields on Y=const surfaces in Tecplot format
c This subroutine is ready to print also UC, WC and P  !!!!!!
        CALL IOFIELD2D_PLT(UO,VO,WO,P,XU,XC,YV,YC,ZW,ZC,ZWG,ZCG,NX,NY,NZ,NZG,TLEVEL)
      ENDIF
      clock(16) = clock(16) + tclock() - clocktemp

C
C-----------------------------------------------------------------------
C                                                    drag and lift coef.
C-----------------------------------------------------------------------
      clocktemp1 = tclock()

      IF(ibm/=0.and.icalf>0 .AND. itcalf>0 .AND. mod(icycle,itcalf)==0) THEN
c        ibd=1
c        clock(60) = tclock()
c        call calcfrc(ub,vb,wb,p,flaguo(1,1,1,ibd),flagvo(1,1,1,ibd)
c     &        ,flagwo(1,1,1,ibd),nx,ny,nz,xu,yv,zw,xc,yc,zc,fb
c     &        ,dtm1,icycle,time,'force',cvlim1)
c        clock(60) = tclock()-clock(60)

c IFLAG_PSB=0 if the values of pressure and shear stress have not been
c estimated above
        if(iflag_psb==0 .OR. (iflag_psb==1 .AND. ivrtx==1)) then

          nimp = sum(mimp(:,:))   !!!!!!
          mimpbd = 0
          clocktemp = tclock()
          DO IBD=1,NBD
ccc            clocktemp = tclock()

            ilb = lb(ibd)+1
            ile = lb(ibd)+mb(ibd)

c              iflag_ptr=1
              if(iflag_ptr==0) then
!!!!!!              call TAGPBD2(FLAGPO(:,:,:,ibd),FLAG,XC,YC,XC_CAR,YC_CAR
!!!!!!     &            ,ZC,ZCG,NX,NY,NZ,NZG,IBD,VERTEXC,UNVECT,NFACET)
              call TAGPBD2(FLAGPO(:,:,:,ibd),FLAG,XC,YC,XC_CAR,YC_CAR
     &,ZC,ZCG,NX,NY,NZ,NZG,IBD,VERTEXC(:,ilb:ile),UNVECT(:,ilb:ile),mb(ibd))

              flagpbd(:,:,:,ibd) = flag(:,:,:)

c              if(ibd.eq.nbd)icom=2   !!!!!!
              icom=-1    !!!!!!
c              iflag_ptr = 1
              IF(IBM==1) THEN
                iflag_ptr = 1
c                icom=1
!!!!!!                IF(ITBDY>0 .AND. ITBDY>ITCALF .AND. IVRTX==1) icom=1
              ENDIF

!!!!!!              nimp = sum(mimp)
!!!!!!              call geom
!!!!!!     &(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs),vertex,vertexc,unvect
!!!!!!     &,limpbd(ibd,1),mimpbd(ibd,1),ip,jp,kp,nimp,xnp,ynp,znp,nxp,nyp,nzp
!!!!!!     &,mrkp,dirp,fp,flag,flagpbd,nx,nyl,nzl,nbd,nfacet,ibd,icom,nflp,nflpbdmax,4,icycle,tlevel)
ccc              call geom
              call geom_mod
     &(xc,yc(jys),xc_car(:,jys),yc_car(:,jys),zcg(kzs)
     &,vertex(:,:,ilb:ile),vertexc(:,ilb:ile),unvect(:,ilb:ile)
     &,limpbd(ibd,:),mimpbd(ibd,:),ip,jp,kp,nimp,xnp,ynp,znp,nxp,nyp,nzp
     &,mrkp,dirp,fp,flag,flagpbd,nx,nyl,nzl,nbd,mb(ibd),ibd,icom,nflp,nflpbdmax(ibd),4,icycle,tlevel)
              call mtrx(xc,yc(jys),zcg(kzs),limpbd(ibd,1),sum(mimpbd(ibd,:)),mrkp
     &,ip,jp,kp,xnp,ynp,znp,nxp,nyp,nzp,pmtrx,pindx,nx,nyl,nzl,nbd,ibd,0)
              endif
          ENDDO

          dp = p
!!!!!!            call correctpres(dp,nx,ny,nz,ip,jp,kp,pmtrx,pim,nxp,nyp,nzp
!!!!!!     &,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect,fp,nfacet
!!!!!!     &,xc_car,yc_car,xc,yc,zc,limpbd(ibd,1),mimpbd(ibd,1),nflpbdmax,tlevel,alfxdt,ibd,mbd)
ccc            call correctpres(dp,nx,ny,nz,ip,jp,kp,pmtrx,pim,nxp,nyp,nzp
ccc     &,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect(:,ilb:ile),fp,mb(ibd)
ccc     &,xc_car,yc_car,xc,yc,zc,limpbd(ibd,:),mimpbd(ibd,:),nflpbdmax(ibd),tlevel,alfxdt,ibd,mbd)
c          call correctpres_mod2(dp,nx,ny,nz,ip,jp,kp,pmtrx,pim,nxp,nyp,nzp
          call correctpres_mod3(dp,nx,ny,nz,ip,jp,kp,pmtrx,pim,nxp,nyp,nzp
     &,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect,fp,nfacet
     &,xc_car,yc_car,xc,yc,zc,limpbd,mimpbd,nflpbdmax,tlevel,alfxdt,mbd,nbd,sum(mimpbd(1:nbd,:)))
          DO IBD=1,NBD
            ilb = lb(ibd)+1
            ile = lb(ibd)+mb(ibd)
            call press_body3(dp,flagpbd,xc,yc,zc,xc_car,yc_car,nx,ny,nz,nbd
     &,pbd,mrkpb,vertexc,unvect,nfacet,ibd,mbd,ilb,ile,tlevel,alfxdt)
ccc            clock(61) = clock(61) + tclock()-clocktemp
          ENDDO
          clock(61) = clock(61) + tclock()-clocktemp

          ub = uo
          vb = vo
          wb = wo

          clocktemp = tclock()

c.....compute forcing in x momentum
!          call momforc_mod2(ub,ub,xu,yc,zc,nx,ny,nz,iu,ju,ku,umtrx,uindx
          call momforc_mod3(ub,ub,xu,yc,zc,nx,ny,nz,iu,ju,ku,umtrx,uindx
     &           ,uim,mrku,diru,nxu,nyu,nzu,nflumax,limu,mimu,0,1,nbd,sum(mimu))
c.....compute forcing in y momentum
!          call momforc_mod2(vb,vb,xc,yv,zc,nx,ny,nz,iv,jv,kv,vmtrx,vindx
          call momforc_mod3(vb,vb,xc,yv,zc,nx,ny,nz,iv,jv,kv,vmtrx,vindx
     &           ,vim,mrkv,dirv,nxv,nyv,nzv,nflvmax,limv,mimv,0,2,nbd,sum(mimv))
c.....compute forcing in z momentum
!          call momforc_mod2(wb,wb,xc,yc,zw,nx,ny,nz,iw,jw,kw,wmtrx,windx
          call momforc_mod3(wb,wb,xc,yc,zw,nx,ny,nz,iw,jw,kw,wmtrx,windx
     &           ,wim,mrkw,dirw,nxw,nyw,nzw,nflwmax,limw,mimw,0,3,nbd,sum(mimw))

          ub(:,1,:) = ub(:,ny-1,:)
          ub(:,ny,:) = ub(:,2,:)
          vb(:,1,:) = vb(:,ny-1,:)
          vb(:,ny,:) = vb(:,2,:)
          wb(:,1,:) = wb(:,ny-1,:)
          wb(:,ny,:) = wb(:,2,:)

          clock(62) = clock(62) + tclock()-clocktemp

          clocktemp = tclock()




          DO ibd=1,mbdn-1
            ilb = lb(ibd)+1
            ile = lb(ibd)+mb(ibd)
            CALL SHEAR_BODY(ub,vb,wb,nx,ny,nz,xu,yv,zw,xc,yc,zc
     &,xu_car,yu_car,xv_car,yv_car,xc_car,yc_car,vertexc,unvect,nfacet
     &,dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb,cf_aux,mrksb,ibd,nbd,ilb,ile,tlevel,0)
          ENDDO
          DO ibd=mbdn,nbd
            ilb = lb(ibd)+1
            ile = lb(ibd)+mb(ibd)
            CALL SHEAR_BODY(ub,vb,wb,nx,ny,nz,xu,yv,zw,xc,yc,zc
     &,xu_car,yu_car,xv_car,yv_car,xc_car,yc_car,vertexc,unvect,nfacet
     &,dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb,cf_aux,mrksb,ibd,nbd,ilb,ile,tlevel,1)
          ENDDO
          clock(63) = clock(63) + tclock()-clocktemp

        endif

c The forces on the immersed boundaries are evaluated and written on a file
c Their components are evaluated in Cartesian coordinates
c The file format is Tecplot !!!!!!
        do ibd=1,nbd
          clocktemp = tclock()
          CALL calcforce(unvect,vertex,vertexc,areaf,pbd,mrkpb,dudxb,dudyb,dudzb
     &         ,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb,mrksb,nfacet,ibd,tlevel,icycle)
          clock(64) = clock(64) + tclock()-clocktemp
        enddo

      ENDIF
      clock(17) = clock(17) + tclock() - clocktemp1
c
C-----------------------------------------------------------------------
C     Record velocity and pressure time signal  at probe locations
C-----------------------------------------------------------------------
C
      call rec_timeprobes_all(uo,vo,wo,p,nx,ny,nz,TIME)
      clocktemp = tclock()
!      if(itmprb>0) then
c the centered values of velocity and pressure are stored in
c UTPRB,VTPRB,WTPRB,PTPRB
c IPRBSER: is the counter associated with the specific time step
!        call rec_timeprobes(uo,vo,wo,p,nx,ny,nz,itprb,jtprb,ktprb
!     &          ,utprb(1,iprbser),vtprb(1,iprbser),wtprb(1,iprbser),ptprb(1,iprbser),ntprb,ntprbmax)

!        if(mod(icycle,itres)==0) then
c the time evolutions of velocity and pressure at the probes are written
c N.B. IO_TIMEPROBES uses MPI_MODE_CREATE instead of MPI_MODE_APPEND
!          call io_timeprobes(filetprb,utprb,vtprb,wtprb,ptprb,ntprb,itres,tmprbindx,ntprbsmpl,ntprbmax)
c the counter IPRBSER is set at 0 after the print
!          iprbser=0
!        endif

c the counter IPRBSER is increased for the next record
!        iprbser = iprbser+1
!      endif
      clock(18) = clock(18) + tclock() - clocktemp

C
C-----------------------------------------------------------------------
C     Compute correlation and eds of azimuthal signal
C-----------------------------------------------------------------------
      clocktemp = tclock()
      if(yspctr>0 .AND. mod(icycle,itysamp)==0) then

        if(nyprbmax>0) then
c correlation of the velocity components and pressure along the direction Y
          call compute_correlations(uo,vo,wo,p,nx,ny,nz,iyprb,jyprb,kyprb,nyprb,nyprbmax,ycor,nycorvar)
        endif

c the azimuthal spectra are evaluated
        if(nyprbmax>0) then
           if(icycle==itysamp) then
             call compute_esd(uo,vo,wo,p,nx,ny,nz,iyprb,kyprb,nyprb,nyprbmax,esd,nyesdvar,0)
           else
             call compute_esd(uo,vo,wo,p,nx,ny,nz,iyprb,kyprb,nyprb,nyprbmax,esd,nyesdvar,1)
           endif
        endif

        if(mod(icycle,itres)==0) then
          ns = real(ycorsmpl)
!number of samples for correlations
          ycorsmpl = ycorsmpl+(itres/itysamp)
          invs = 1./real(ycorsmpl)
!total value of correlations
          ycorprb = (ycorprb*ns+ycor)*invs
!correlations are written
          call write_correlations(fileycor,ycorprb,ny-2,nyprb,nycorvar,yprbindx,nyprbmax,ycorsmpl)
          ycor = 0.0

          ns = real(yesdsmpl)
!number of samples for spectra
          yesdsmpl = yesdsmpl+(itres/itysamp)
          invs = 1./real(yesdsmpl)
!total value of spectra
          esdprb = (esdprb*ns+esd)*invs
!spectra are written
          call write_correlations(fileyspc,esdprb,ny/2,nyprb,nyesdvar,yprbindx,nyprbmax,yesdsmpl)
          esd = 0.0

        endif


      endif
      clock(19) = clock(19) + tclock() - clocktemp

c
c-----------------------------------------------------------------------
C     Write VPfield files
c The instantaneous fields of velocity and pressure are written
c-----------------------------------------------------------------------
c
      clocktemp = tclock()
      if(itVPfield>0 .AND. mod(icycle,itVPfield)==0) then
        nVPfield = nVPfield+1
!        call write_VPfield(trim(hspdir),'field.'//index(nVPfield),uo,vo,wo,p,nx,ny,nz,limVPfield)
!        call write_VPfield_hdf5_sp(trim(hspdir),'field.'//index(nVPfield),uo,vo,wo,p,nx,ny,nz,limVPfield)
        if(myrank.eq.0) then
           open(unit=10,file='time.VPfield',form='formatted',position
     $          ='append')
           write(10,*) nVPfield,tlevel
           close(10)
        endif
      endif
      clock(20) = clock(20) + tclock() - clocktemp
c
C-----------------------------------------------------------------------
C     Compute time averaged primary variables
C-----------------------------------------------------------------------
      clocktemp = tclock()
      if(itmavg>0 .AND. mod(icycle,ftmavg)==0) then
c the averaged velocity and pressure are evaluated
        call compute_primevar_tmavg(vartmavg,nxt,nyt,nzt,4,uo,vo,wo,p,nx,ny,nz,limtmavg,nsmpltmavg)
c the averaged squared velocity and pressure are evaluated
        call compute_Restresses_tmavg(vartmavg(1,1,1,5),nxt,nyt,nzt,4,uo,vo,wo,p,nx,ny,nz,limtmavg,nsmpltmavg(5))
        call compute_crossvel_tmavg(vartmavg(1,1,1,9),nxt,nyt,nzt,nvartmavg-8,uo,vo,wo,tv,nx,ny,nz,limtmavg,nsmpltmavg(9))
        nsmpltmavg = nsmpltmavg+1
        if(myrank.eq.0)write(6,*)'output3',nsmpltmavg

        if(mod(icycle,itres)==0) then
c the averaged values are written on a file
c          call write_primevar_tmavg(filetmavg,vartmavg,nxt,nyt,nzt,nvartmavg,nsmpltmavg,limtmavg,nz)
          call write_primevars_tmavg(hspdir,vartmavg,nxt,nyt,nzt,nvartmavg,nsmpltmavg,limtmavg,nz)
        endif

      endif

      if(iregtmavg>0 .AND. mod(icycle,fregtmavg)==0) then
c the averaged velocity and pressure are evaluated
        call compute_regvar_tmavg(regtmavg,nreg,nvarreg,iindtmavg
     &        ,kmintmavg,kmaxtmavg,j1reg,j2reg,uo,vo,wo,p,nx,ny,nz,nsmplregtmavg)

        if(mod(icycle,itres)==0) then
c the averaged values are written on a file
          call write_regvar_tmavg(filereg,regtmavg,nreg,nvarreg,iindtmavg,nz
     &          ,kmintmavg,kmaxtmavg,nsmplregtmavg,nregprev,nregtot)
        endif
      endif
      clock(21) = clock(21) + tclock() - clocktemp

c
C-----------------------------------------------------------------------
C     Write files for VP region
C-----------------------------------------------------------------------
      if(VPreg_freq>0 .AND. mod(icycle,VPreg_freq)==0) then
        fileVPreg = trim(hspdir)//'VPreg.'//index(iVPreg)//'.bin'
        call write_var_VPreg(fileVPreg,uo,vo,wo,p,nx,ny,nz,iindVPreg
     $     ,j1VPreg,j2VPreg,kminVPreg,kmaxVPreg,nprevVPreg,nVPregtot,tlevel)
        iVPreg = iVPreg+1
      endif

c
C-----------------------------------------------------------------------
C     Write CFL ave
C-----------------------------------------------------------------------
      if(icflave>0 .AND. mod(icycle,ITRES)==0) then
        call write_cflave(filecflave,cflave,nx,nsmplcfl)
      endif
c
C-----------------------------------------------------------------------
C     Write stat 1D
C-----------------------------------------------------------------------
      if(istat1D>0) then
        if(mod(icycle,itstat1D)==0) then
c the 1 dimensional statistics (along the X direction) are evaluated
c at the centered positions
c the values are averaged along the Y and Z directions
          call calc_stat1d_ave(uo,vo,wo,p,tv,stat1d,nvarstat1D,nx,ny,nz,nzg,nstat1D)
          if(isgs>0 .AND. isgs<5) call calc_stat1d_ave_LES(cles,clesp
     &         ,clesn,nilmp,nilmn,ilm,imm,g,sxx,syy,szz,sxy,sxz,syz,stat1d,nvarstat1D,nx,ny,nz,nzg,nstat1D)
          nstat1D = nstat1D + 1
        endif

        if(mod(icycle,itres)==0) then
c the file of the statistics is written taking into account the current and the previous values
          call write_stat1d_ave(stat1dave,stat1d,xc,nx,nvarstat1D,nstat1Dave,nstat1D)
          stat1d = 0.0
          nstat1D = 0
        endif
      endif
c
c.... Bulk velocity for calculating dpdz
      if(myrank==0 .AND. itype(5)==500 .AND. mod(icycle,ITRES)==0) then
        open(unit=10,file='ubulk_old.res',form='formatted')
        write(10,*) ubold
        write(10,*) dpdz
        close(10)
      endif

C-----------------------------------------------------------------------
C     Write utau
C-----------------------------------------------------------------------
      if(iutau>0) then
        if(mod(icycle,itutau)==0) then
          itu = itu+1
c the skin-friction velocity is evaluated as mean value along the Y and Z directions
c the cylindrical (utaudim=1) and the Cartesian case (utaudim=2) are distinguished
          call calc_utau(wo,xu,xc,nx,ny,nz,utau(itu,:),utaudim,dpdzutau(itu))
c the bulk velocity is evaluated
c          call calc_ubulk(wo(:,:,1),nx,ny,1,utau(itu,1))
          timeutau(itu) = tlevel
        endif

        if(mod(icycle,itutau)==0) then
c the values of the skin-friction velocity are written on a file
          call write_utau(utau,utaudim,dpdzutau,timeutau,ntutau,iutau)
          itu = 0
        endif
      endif
c
C-----------------------------------------------------------------------
C     create time series about mean flow
C-----------------------------------------------------------------------
C

      IF(ICYCLE>=ITINI) THEN

C-----------------------------------------------------------------------
c     one or two dimensional statistics, three dimensional fields
C-----------------------------------------------------------------------
c
        IF(ICFL==1.AND.ISTAT/=0.AND.MOD(ICYCLE,ISTAT)==0) GOTO 9999

c        IF(ICFL/=0.AND.MOD(TIME-DTM1,DTSAVE)<DTSAVE
c     &       .AND.MOD(TIME-DTM1,DTSAVE)+DTM1>=DTSAVE) GOTO 9999
c        IF(ICFL==0.AND.ISTAT/=0.AND.MOD(ICYCLE,ISTAT)==0) GOTO 9999

        GOTO 8888

 9999   CONTINUE

          NWRITE=NWRITE+1

c NWRITE can vary between 0 and 100    !!!!!!
          IF(NWRITE>100) NWRITE = NWRITE-100
c
          IF(MYRANK==0) write(6,*) ' Begin write field ...',icycle,time,nwrite

c Tecplot files of the moving immersed-boundaries are written
!          do ibd=mbd,nbd
!            call imb2plt('mbd.'//index(ibd)//'.'//index(nwrite)//'.plt',vertex,nfacet,ibd)
!          enddo

c STAT1D evaluates the one-dimensional statistics along X
c        CALL STAT1D(UO,VO,WO,P,TV,'stat1d.'//INDEX(NWRITE),NX,NY,NZ,TIME)
c STAT2D evaluates the two-dimensional statistics on XZ surfaces
c        CALL STAT2D(UO,VO,WO,P,TV,'stat2d.'//INDEX(NWRITE),NX,NY,NZ,TIME)
          IF(MYRANK==0) write(6,*) ' End write field'

 8888     CONTINUE

C        ENDIF

      ENDIF

!      if((time-timep).gt.((periodo)*(iplant8+1))) then
!         clocktemp = tclock()
c
c        Evaluation of the vorticity values for the subroutine
c        mediestaggavg
c
!         call vort1staggavg(uo,vo,wo,nx,ny,nz,xc,xu,zc)
!         call vort1staggnoavg(uo,vo,wo,nx,ny,nz,xc,xu,zc)
c
c        Evaluation of the time averaged values along each direction
c
!         call mediestaggavg(nx,ny,nz,uo,vo,wo,p,tv,dtm1)
!         call mediestaggnoavg(nx,ny,nz,uo,vo,wo,p,tv,dtm1,nbd,flagpo,time)
!         call mediestaggnoavgcart(nx,ny,nz,yc,yv,uo,vo,wo,p,tv,dtm1,nbd,flagpo,time)
c
!         if(myrank.eq.0)write(6,*)'output1',tempo
!         clock(50) = clock(50) + tclock() - clocktemp
!      endif
c
!      if(mod(icycle,10000000).eq.0) then
!         clocktemp = tclock()
c
c     Print of the instantaneous fields (velocities, pressure,
c     vorticities and stagnation pressure)
c
!         call vort2(nx,ny,nz,xc,xu,zc,uo,vo,wo)
!         call plot2D(iplant2,iplant2c,iplant2d,iplant2u,
!     %nx,ny,nz,nzg,icycle,nbd,mbd,nfacet,unvect,vertex,xc,yc,zcg,
!     %xc_car,yc_car,flagpo,tv,vo,uo,wo,p)
!         call plot2Dvort(iplant2c,iplant2d,iplant2u,
!     %nx,ny,nz,nzg,nbd,icycle,flagpo,vo,uo,wo,p)
!         iplant2=iplant2+1
!         clock(49) = clock(49) + tclock() - clocktemp
!      endif
c
!      if((time-timep).gt.((periodo5)*(iplant9+1))) then
         if(mod(icycle,itpln).eq.0) then
         clocktemp = tclock()
         iplant9=iplant9+1
c
c     Print of the instantaneous fields in the staggered
c     positions
c
         if(idens.eq.1) then
!         call plotinststagg(iplant9c,iplant9d,iplant9p,
!         call plotinststagg_xyz(iplant9c,iplant9d,iplant9p,
!     %nx,ny,nz,nzg,nbd,icycle,xc,xu,yc,yv,zcg,zwg,p,vo,uo,wo)
!     %,flaguo,flagvo,flagwo)
        call WRITE_PLANES(NX,NY,NZ,NZG,NBD,ICYCLE,time,DTM1,XC,XU,YC,YV,ZC,ZCG,ZWG,xc_car,
     &                  yc_car,xu_car,yu_car,xv_car,yv_car,P,VO,UO,WO,DENS)

!         call plotinststaggcart(iplant9c,iplant9d,iplant9p,
!     %nx,ny,nz,nzg,nbd,icycle,xc,xu,yc,yv,zcg,zwg,p,vo,uo,wo,
!     %flaguo,flagvo,flagwo)
!         call vort2stagg(nx,ny,nz,xc,xu,zc,uo,vo,wo)
! !         call vort2stagg_xyz(nx,ny,nz,xc,xu,zc,uo,vo,wo)
! !         call vort2staggcart(nx,ny,nz,xc,xu,yc,yv,zc,uo,vo,wo)
! !         call plotvortstagg(iplant9c,iplant9d,iplant9p,
!             call plotvortstagg_xyz(iplant9c,iplant9d,iplant9p,
!       %nx,ny,nz,nzg,nbd,icycle,p,uo,vo,wo,tv)
!     %,flaguo,flagvo,flagwo)
         else
!           call plotinststagg_d(iplant9c,iplant9d,iplant9p,
!     %nx,ny,nz,nzg,nbd,icycle,xc,xu,yc,yv,zcg,zwg,p,vo,uo,wo,dens)
!           call plotinststagg_xyz_d(iplant9c,iplant9d,iplant9p,
!     %nx,ny,nz,nzg,nbd,icycle,xc,xu,yc,yv,zcg,zwg,p,vo,uo,wo)
!     %,flaguo,flagvo,flagwo)
c        write(*,*), "Eddy6.f u(nx-1,1,nz-1)
c     &  u(nx-1,ny-1,nz-1) 3898",WO(nx-1,1,nz-1), WO(nx-1,ny-1,nz-1)  
c
 	  call WRITE_PLANES(NX,NY,NZ,NZG,NBD,ICYCLE,time,DTM1,XC,XU,YC,YV,ZC,ZCG,ZWG,xc_car,
     &                  yc_car,xu_car,yu_car,xv_car,yv_car,P,VO,UO,WO,DENS)
!           call vort2stagg_xyz(nx,ny,nz,xc,xu,zc,uo,vo,wo)
!           call plotvortstagg_xyz_d(iplant9c,iplant9d,iplant9p,
!     %nx,ny,nz,nzg,nbd,icycle,dens,uo,vo,wo,tv,xc,xu)
!      %,flaguo,flagvo,flagwo)


         endif
c
         clock(49) = clock(49) + tclock() - clocktemp
         endif

c$$$          if(mod(icycle,10000).eq.0) then
c$$$          call PERIODIC_1D_AVERAGE_RMS(NX,NY,NZ,NZG,NBD,ICYCLE,time,DTM1,XC,XU,YC,YV,ZCG,ZWG,xc_car,
c$$$     &                              yc_car,xu_car,yu_car,xv_car,yv_car,P,VO,UO,WO,DENS)
c$$$          endif


C           if(mod(icycle,100).eq.0) then
C           call PERIODIC_2D_AVERAGE_RMS_NCPU(NX,NY,NZ,NZG,NBD,ICYCLE,time,DTM1,XC,XU,YC,YV,ZCG,ZWG,xc_car,
C      &                              yc_car,xu_car,yu_car,xv_car,yv_car,P,VO,UO,WO,DENS,TV)
C           if(myrank==0) write(*,*) "DONE PERIODIC_2D_AVERAGE_RS"
C           endif




C           if(mod(icycle,100).eq.0) then
C           call PERIODIC_2D_DISSIPATION_NCPU(NX,NY,NZ,NZG,NBD,ICYCLE,time,DTM1,XC,XU,YC,YV,ZCG,ZWG,xc_car,
C      &                              yc_car,xu_car,yu_car,xv_car,yv_car,P,VO,UO,WO,DENS,TV)
C           if(myrank==0) write(*,*) "DONE PERIODIC_2D_DISSIPATION_NCPU"
C           endif
c




!      if((time-timep).gt.((periodo)*(iplant4+1))) then
!         clocktemp = tclock()
!         call mediecalcstaggavg
!     %(nx,ny,nz,nzg,nbd,time,dtm1,vo,uo,wo,p,tv,xc,xu,yc,yv,zc,
!     %zcg,zwg,flaguo,flagvo,flagwo,icycle)
!!    %mbd,nfacet,unvect,vertex)
!         call mediecalcstagg
!     %(nx,ny,nz,nzg,nbd,time,dtm1,vo,uo,wo,p,tv,xc,xu,yc,yv,zc,
!     %zcg,zwg,flagpo,flaguo,flagvo,flagwo,icycle)
!!    %mbd,nfacet,unvect,vertex)
!         call mediecalcstaggcart
!     %(nx,ny,nz,nzg,nbd,time,dtm1,vo,uo,wo,p,tv,xc,xu,yc,yv,zc,
!     %zcg,zwg,flagpo,flaguo,flagvo,flagwo,icycle)
!!    %mbd,nfacet,unvect,vertex)
!         if(myrank.eq.0)write(6,*)'output2',tempo3
!         clock(55) = clock(55) + tclock() - clocktemp
!      endif
c
      clocktemp = tclock()
      call sondestagg(nx,ny,nz,nbd,time,dtm1,vo,uo,wo,p)    !,flagpo)
!      call sondestaggcart(nx,ny,nz,nbd,time,dtm1,yv,vo,uo,wo,p,flagpo)
      clock(57) = clock(57) + tclock() - clocktemp
c
      if((itres.ne.0).and.(mod(icycle,itres).eq.0)) then
         clocktemp = tclock()
         call continua(nx,ny,nz)
         clock(58) = clock(58) + tclock() - clocktemp
      endif
c
!!!!!!
C
C-----------------------------------------------------------------------
C     end main loop
C-----------------------------------------------------------------------
C
      TIME_ITERATION = tclock() - TIME_ITERATION
      TIME_LAST_ITERATION = TIME_ITERATION
      clock(99) = TIME_LAST_ITERATION
      clock(100) = clock(100) + time_iteration
C
C-----------------------------------------------------------------------
C     screen information
C-----------------------------------------------------------------------

      IF(MOD(ICYCLE,ITSCR)==0) THEN

        CALL SCRINFO(UO,VO,WO,P,DP,TV,DENS,DTM1,TIME,ICYCLE,ID,JD,KD,NX,NY,NZ)

c        TIME_ITERATIONS = (tclock() - TIME_ITERATIONS)/REAL(ITSCR)
c        CLOCK(100) = TIME_ITERATIONS

        CLOCK(91) = CLOCK(24)+CLOCK(28)+CLOCK(32) !TAG3D
        CLOCK(92) = CLOCK(25)+CLOCK(29)+CLOCK(33) !TAGU
        CLOCK(93) = CLOCK(26)+CLOCK(30)+CLOCK(34) !GEOMU
        CLOCK(94) = CLOCK(27)+CLOCK(31)+CLOCK(35) !MTRXU

        clock(1:98) = clock(1:98)/real(itscr)
        clock(100) = clock(100)/real(itscr)
        CALL MPI_REDUCE(CLOCK,CLOCKG,NCLOCKS,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(CLOCK,CLOCKGMIN,NCLOCKS,MTYPE,MPI_MIN,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(CLOCK,CLOCKGMAX,NCLOCKS,MTYPE,MPI_MAX,0,MPI_COMM_EDDY,IERR)
        CLOCKG=CLOCKG/REAL(MYSIZE)

        IF(MYRANK==0 .AND. IOCLOCK>0) THEN
c the statistics on the computational times are written
          OPEN(UNIT=16,FILE='clock.dat',FORM='FORMATTED'
     &        ,POSITION='APPEND')
c global time statistics
          WRITE(16,'(2(A,1x,I6))') 'Cycles=',icycle-itscr,'-',icycle
          WRITE(16,'(2A)') ' Task/Time          '
     &         ,'        Ave.(sec/%)     Max.(sec/%)     Min.(sec/%)'
          WRITE(16,'(A,3(10x,F8.4))') 'Total                  :',CLOCKG(100),CLOCKGMAX(100),CLOCKGMIN(100)
c time statistics for single sections of the code
          i=1
          WRITE(16,102) 'Calc. CFL              :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          IF(IBM.GT.0) THEN
            i=2
            WRITE(16,102) 'Moving Boundary Tasks  :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
!!!!!!          i=21
            i=51
            WRITE(16,102) '  - RBM                :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(2),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(2),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(2),'%)'
            i=22
            WRITE(16,102) '  - CALC_DS_TRIANGLES  :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(2),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(2),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(2),'%)'
            i=23
            WRITE(16,102) '  - LIMIMB             :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(2),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(2),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(2),'%)'
            i=91
            WRITE(16,102) '  - TAG3D (U,V,W)      :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(2),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(2),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(2),'%)'
            i=92
            WRITE(16,102) '  - TAGU  (U,V,W)      :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(2),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(2),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(2),'%)'
            i=93
            WRITE(16,102) '  - GEOM  (U,V,W)      :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(2),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(2),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(2),'%)'
            i=94
            WRITE(16,102) '  - MTRX (U,V,W)       :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(2),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(2),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(2),'%)'
            i=36
            WRITE(16,102) '  - FLAGP (P)          :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(2),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(2),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(2),'%)'
            i=37
            WRITE(16,102) '  - GEOMP (P)          :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(2),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(2),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(2),'%)'
            i=38
            WRITE(16,102) '  - MTRX (P)           :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(2),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(2),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(2),'%)'
            i=56
            WRITE(16,102) '  - Set surface vel.   :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(2),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(2),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(2),'%)'
          ENDIF
          i=3
          WRITE(16,102) 'Calc. RHS              :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=4
          WRITE(16,102) 'Predictor              :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=40
          WRITE(16,102) '  - U loop             :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(4),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(4),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(4),'%)'
          i=41
          WRITE(16,102) '    - 1st part         :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(40),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(40),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(40),'%)'
          i=42
          WRITE(16,102) '    - 2nd part         :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(40),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(40),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(40),'%)'
          i=43
          WRITE(16,102) '  - V loop             :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(4),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(4),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(4),'%)'
          i=44
          WRITE(16,102) '    - 1st part         :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(43),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(43),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(43),'%)'
          i=45
          WRITE(16,102) '    - 2nd part         :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(43),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(43),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(43),'%)'
          i=46
          WRITE(16,102) '  - W loop             :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(4),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(4),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(4),'%)'
          i=47
          WRITE(16,102) '    - 1st part         :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(46),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(46),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(46),'%)'
          i=48
          WRITE(16,102) '    - 2nd part         :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(46),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(46),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(46),'%)'
          IF (implx/=0.OR.imply/=0) THEN
            i=5
            WRITE(16,102) 'Mom. Forcing (Implicit):'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
            i=6
            WRITE(16,102) 'Inverse LHS (Implicit) :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          ELSE
            i=7
            WRITE(16,102) 'MoM. Forcing (Explicit):'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          ENDIF
          i=8
          WRITE(16,102) 'Boundary               :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=9
          WRITE(16,102) 'Divergence             :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=10
          WRITE(16,102) 'Pressure Direct        :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=11
          WRITE(16,102) 'Corrector              :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=12
          WRITE(16,102) 'Update/Refresh Pressure:'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=13
          WRITE(16,102) 'Press. Field Extension :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          IF (ISGS/=0) THEN
            i=14
            WRITE(16,102) 'Turbulent viscosity    :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          ENDIF
          i=15
          WRITE(16,102) 'Shear stress on body   :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=16
          WRITE(16,102) 'Input/Output (restart) :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=17
          WRITE(16,102) 'Drag and Lift          :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=61
          WRITE(16,102) '  - Correct+Press_body :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=62
          WRITE(16,102) '  - Momforc            :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=63
          WRITE(16,102) '  - Shear_body         :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=64
          WRITE(16,102) '  - Calcforce          :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=18
          WRITE(16,102) 'Time probes            :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=19
          WRITE(16,102) 'Y correlations         :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=20
          WRITE(16,102) 'VPfield files          :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=21
          WRITE(16,102) 'Time averaged files    :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=50
          WRITE(16,102) 'Medie                  :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=49
          WRITE(16,102) 'Instantaneous fields   :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=55
          WRITE(16,102) 'Time-averaged fields   :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=57
          WRITE(16,102) 'Sonde                  :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=58
          WRITE(16,102) 'Continua               :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          i=59
          WRITE(16,102) '3D fields              :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(100),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(100),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(100),'%)'
          WRITE(16,'(A,I4,A,3(1x,F8.4))') 'Time per ',ITSCR,' iters:',CLOCKG(100),CLOCKGMAX(100),CLOCKGMIN(100)
          CLOSE(16)
        ENDIF

c time statistics on the last iteration
        IF(MYRANK.EQ.0) THEN
          WRITE(6,'(2A,3(1x,F8.4))') 'Time for last iteration'
     &              ,'(Ave-Max-Min):',CLOCKG(99),CLOCKGMAX(99),CLOCKGMIN(99)
          WRITE(6,'(A,I4,A,3(1x,F8.4))') 'Time per ',ITSCR,' iters:',CLOCKG(100),CLOCKGMAX(100),CLOCKGMIN(100)
        ENDIF

c time statistics on another file
C        IF(MYRANK.EQ.0 .AND. IOCLOCK>0) THEN
C          OPEN(UNIT=16,FILE='clock.csv',FORM='FORMATTED'
C     &         ,POSITION='APPEND')
C          WRITE(16,'(2(A,1x,I6))') 'Cycles=',icycle-itscr,'-',icycle
C          WRITE(16,'(2A)') ' Task/Time          '
C     &         ,'        Ave.(sec/%)     Max.(sec/%)     Min.(sec/%)'
C          i=100   !!!!!! 100 instead of 99
C          WRITE(16,103) clockg(i),',',clockgmax(i),',',clockgmin(i)
C          i=1
C          WRITE(16,103) clockg(i),',',clockgmax(i),',',clockgmin(i)
C          i=3
C          WRITE(16,103) clockg(i),',',clockgmax(i),',',clockgmin(i)
C          i=4
C          WRITE(16,103) clockg(i),',',clockgmax(i),',',clockgmin(i)
C          IF (implx/=0.OR.imply/=0) THEN
C            i=6
C            WRITE(16,103) clockg(i),',',clockgmax(i),',',clockgmin(i)
C          ENDIF
C          i=8
C          WRITE(16,103) clockg(i),',',clockgmax(i),',',clockgmin(i)
C          i=9
C          WRITE(16,103) clockg(i),',',clockgmax(i),',',clockgmin(i)
C          i=10
C          WRITE(16,103) clockg(i),',',clockgmax(i),',',clockgmin(i)
C          i=11
C          WRITE(16,103) clockg(i),',',clockgmax(i),',',clockgmin(i)
C          i=12
C          WRITE(16,103) clockg(i),',',clockgmax(i),',',clockgmin(i)
C          i=16
C          WRITE(16,103) clockg(i),',',clockgmax(i),',',clockgmin(i)

C          write(16,104) 'Ave:',clockg(100),clockg(1),clockg(3),clockg(4),clockg(6)
C     &,clockg(8),clockg(9),clockg(10),clockg(11),clockg(12),clockg(16)  !!!!!!
C          write(16,104) 'Max:',clockgmax(100),clockgmax(1),clockgmax(3),clockgmax(4),clockgmax(6)
C     &,clockgmax(8),clockgmax(9),clockgmax(10),clockgmax(11),clockgmax(12),clockgmax(16)  !!!!!!
C          write(16,104) 'Min:',clockgmin(100),clockgmin(1),clockgmin(3),clockgmin(4),clockgmin(6)
C     &,clockgmin(8),clockgmin(9),clockgmin(10),clockgmin(11),clockgmin(12),clockgmin(16)  !!!!!!
C          CLOSE(16)
C
C        ENDIF
        TIME_ITERATIONS = tclock()
        clock = 0.0
c
      ENDIF
C
      enddo
c
!!!!!!
      close(49)
!      do i=121,132
!        close(i)
!      enddo
!      do i=1,nprbmax2
!        ii=prbindx2+(i-1)
!        close(700+ii)
!        close(800+ii)
!      enddo
!!!!!!
      IF(MYRANK==0) THEN
         WRITE(6,'(A)') '*.........................'
         WRITE(6,'(A)') '*..EDDY finished.'
      ENDIF

      CALL MPI_FINALIZE(IERR)
C
      STOP
101   FORMAT(A,3(I4,1x),E15.8)
102   FORMAT(A,3(1x,F8.4,A,F4.1,A))
103   FORMAT(F8.4,A,F8.4,A,F8.4)
104   FORMAT(A,20(1x,F8.4))

      END
