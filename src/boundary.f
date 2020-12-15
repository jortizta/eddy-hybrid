* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*                 boundary conditions
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      SUBROUTINE BOUNDARY(US,VS,WS,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TLEVEL)
c
      INCLUDE 'common.h'
c
      INTEGER NX,NY,NZ,IDIV
      REAL    TLEVEL
      REAL    XU(NX),XC(NX),YV(NY),YC(NY),ZW(NZ),ZC(NZ)
      REAL US(NX,NY,NZ),VS(NX,NY,NZ),WS(NX,NY,NZ)
c
      INTEGER IBND
c
c.....wall boundary conditions
c
      DO IBND=1,6

        SELECT CASE (ITYPE(IBND))
c
c.....non-slip wall           
        CASE ( 50)
          CALL BOUND050(US,VS,WS,NX,NY,NZ,IBND)
c
c.....slip wall           
        CASE ( 60)
          CALL BOUND060(US,VS,WS,NX,NY,NZ,IBND)
c
c.....moving wall           
        CASE ( 70)
          CALL BOUND070(US,VS,WS,NX,NY,NZ,IBND)
c
c.....zero vorticity wall           
        CASE ( 80)
          CALL BOUND080(US,VS,WS,XU,XC,ZW,ZC,NX,NY,NZ,IBND)
c
c.....taylor-green velocity
        CASE ( 81)
          CALL BOUND081(US,VS,WS,XU,XC,YV,YC,ZW,ZC,NX,NY,NZ,TLEVEL,IBND)
c
c.....propagation of waves in a linearly stratified fluid
        CASE ( 82)
          CALL BOUND082(US,VS,WS,NX,NY,NZ,IBND)
c     
c.....axis           
        CASE (300)
          CALL BOUND300(US,VS,WS,NX,NY,NZ,IBND)

        END SELECT

      ENDDO
c
c.....periodic boundary conditions
c
      DO IBND=1,3,2

        IF(ITYPE(IBND)==500) CALL BOUND500(US,VS,WS,NX,NY,NZ,IBND)

      ENDDO
c
      RETURN
      END


c----------------------------------------------------------------------
      SUBROUTINE BOUND050(ubc,vbc,wbc,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c     boundary conditions for non-slip wall
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    ubc(nx,ny,nz),vbc(nx,ny,nz),wbc(nx,ny,nz)
c
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (1) 
        ubc(ib1,:,:) =  0.
        vbc(ib1,:,:) = -vbc(ix1,:,:)
        wbc(ib1,:,:) = -wbc(ix1,:,:)      

      CASE (2)
        ubc(ix2,:,:) =  0.
        vbc(ib2,:,:) = -vbc(ix2,:,:)
        wbc(ib2,:,:) = -wbc(ix2,:,:)              

      CASE (3)
        ubc(:,jb1,:) = -ubc(:,jy1,:)
        vbc(:,jb1,:) =  0.
        wbc(:,jb1,:) = -wbc(:,jy1,:)      
         
      CASE (4)
        ubc(:,jb2,:) = -ubc(:,jy2,:)
        vbc(:,jy2,:) =  0.
        wbc(:,jb2,:) = -wbc(:,jy2,:)              
        
      CASE (5)
        ubc(:,:,kb1) = -ubc(:,:,kz1)
        vbc(:,:,kb1) = -vbc(:,:,kz1)
        wbc(:,:,kb1) =  0.
        
      CASE (6)
        ubc(:,:,kb2) = -ubc(:,:,kz2)
        vbc(:,:,kb2) = -vbc(:,:,kz2)
        wbc(:,:,kz2) =  0.              
        
      END SELECT
c
      RETURN
      END

c----------------------------------------------------------------------
      SUBROUTINE BOUND060(ubc,vbc,wbc,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c     boundary conditions for slip wall
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    ubc(nx,ny,nz),vbc(nx,ny,nz),wbc(nx,ny,nz)
c
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (1) 
        ubc(ib1,:,:) = 0.
        vbc(ib1,:,:) = vbc(ix1,:,:)
        wbc(ib1,:,:) = wbc(ix1,:,:)      

      CASE (2)
        ubc(ix2,:,:) = 0.
        vbc(ib2,:,:) = vbc(ix2,:,:)
        wbc(ib2,:,:) = wbc(ix2,:,:)              

      CASE (3)
        ubc(:,jb1,:) = ubc(:,jy1,:)
        vbc(:,jb1,:) = 0.
        wbc(:,jb1,:) = wbc(:,jy1,:)      
         
      CASE (4)
        ubc(:,jb2,:) = ubc(:,jy2,:)
        vbc(:,jy2,:) = 0.
        wbc(:,jb2,:) = wbc(:,jy2,:)              
        
      CASE (5)
        ubc(:,:,kb1) = ubc(:,:,kz1)
        vbc(:,:,kb1) = vbc(:,:,kz1)
        wbc(:,:,kb1) = 0.
        
      CASE (6)
        ubc(:,:,kb2) = ubc(:,:,kz2)
        vbc(:,:,kb2) = vbc(:,:,kz2)
        wbc(:,:,kz2) = 0.              
        
      END SELECT
c
      RETURN
      END
c




c----------------------------------------------------------------------
      SUBROUTINE BOUND070(ubc,vbc,wbc,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c     boundary conditions for moving wall
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    ubc(nx,ny,nz),vbc(nx,ny,nz),wbc(nx,ny,nz)
c
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (1) 
        ubc(ib1,:,:) =    UMOV1
        vbc(ib1,:,:) = 2.*VMOV1*RU(1  )-vbc(ix1,:,:)
        wbc(ib1,:,:) = 2.*WMOV1        -wbc(ix1,:,:)      

      CASE (2)
        ubc(ix2,:,:) =    UMOV2
        vbc(ib2,:,:) = 2.*VMOV2*RU(IX2)-vbc(ix2,:,:)
        wbc(ib2,:,:) = 2.*WMOV2        -wbc(ix2,:,:)              

      CASE (3)
        ubc(:,jb1,:) = 2.*UMOV1-ubc(:,jy1,:)
        vbc(:,jb1,:) =    VMOV1
        wbc(:,jb1,:) = 2.*WMOV1-wbc(:,jy1,:)      
         
      CASE (4)
        ubc(:,jb2,:) = 2.*UMOV2-ubc(:,jy2,:)
        vbc(:,jy2,:) =    VMOV2
        wbc(:,jb2,:) = 2.*WMOV2-wbc(:,jy2,:)              
        
      CASE (5)
        ubc(:,:,kb1) = 2.*UMOV1-ubc(:,:,kz1)
        vbc(:,:,kb1) = 2.*VMOV1-vbc(:,:,kz1)
        wbc(:,:,kb1) =    WMOV1
        
      CASE (6)
        ubc(:,:,kb2) = 2.*UMOV2-ubc(:,:,kz2)
        vbc(:,:,kb2) = 2.*VMOV2-vbc(:,:,kz2)
        wbc(:,:,kz2) =    WMOV2              
        
      END SELECT
c
      RETURN
      END
c


c----------------------------------------------------------------------
      SUBROUTINE BOUND080(ubc,vbc,wbc,xu,xc,zw,zc,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c     Blasius boundary conditions for top wall
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
c.... Input/Output arrays
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    xu(nx),xc(nx),zw(nz),zc(nz)
      REAL    ubc(nx,ny,nz),vbc(nx,ny,nz),wbc(nx,ny,nz)
c
c.... Local arrays
      INTEGER i,j,k
c
c.... Blasius
      INTEGER, SAVE :: NETA
      REAL XO,ZO,DETA,ETA1,ETA2,ALPHA,BETA,UTEMP,VETA1,VETA2
      REAL wtemp(nz)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: ETA,F1,F2
      INTEGER, SAVE :: ISAVE=0
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (2)

c...... Blasius
        xo = dummyfp4
        if(ibm==0) xo=0.0
        zo = dummyfp3 !Re=100, eta=5, x=1
        deta = 0.05

c...... Initialize parameters for Blasius
        if(isave==0) then
          isave = 1
          eta1 = (xu(nx)-xo)*sqrt(0.5/ru1/(zo+zc(1)-zmin))
          eta2 = (xu(nx)-xo)*sqrt(0.5/ru1/(zo+zw(nz)-zw(1)))
          neta = ceiling(max(eta1,eta2)/deta)+2
          allocate(eta(neta),f1(neta),f2(neta))
          eta = 0.0
          do i=2,neta
            eta(i) = deta*real(i-1)
          enddo
          alpha = dummyfp1
          beta = dummyfp2
          call falkner_skan(eta,f1,f2,neta,alpha,beta)
        endif
        
        do k=1,nz
          eta2 = (xu(ix2)-xo)*sqrt(0.5/ru1/(zo+zc(k)-zmin))
          call locate(eta,neta,eta2,i)
          veta1 = sqrt(0.5*ru1/(zo+zc(k)-zmin))*(eta(i)*f2(i)-f1(i))
          veta2 = sqrt(0.5*ru1/(zo+zc(k)-zmin))*(eta(i+1)*f2(i+1)-f1(i+1))
          ubc(ix2,:,k) = veta1 + (veta2-veta1)/deta*(eta2-eta(i))          
        enddo

        VBC(ib2,:,:) = VBC(ix2,:,:)
        WBC(ib2,:,:) = WBC(ix2,:,:)

      END SELECT
c
      RETURN
      END
c
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE BOUND081(ubc,vbc,wbc,xu,xc,yv,yc,zw,zc,nx,ny,nz,time,ibnd)
c----------------------------------------------------------------------
c
c     Taylor Green vortex boundary condition for u-velocity
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
c.... Input/Output arrays
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    time
      REAL    xu(nx),xc(nx),yv(ny),yc(ny),zw(nz),zc(nz)
      REAL    ubc(nx,ny,nz),vbc(nx,ny,nz),wbc(nx,ny,nz)
c
c.... Local arrays
      INTEGER J,K
c
c.... Functions
      REAL   uvel_taylor_green,wvel_taylor_green

      SELECT CASE(IBND)

      CASE (1)
        DO K=1,NZ
        DO J=1,NY
          ubc(1,j,k) = uvel_taylor_green((xu(1)-xmin),(yc(j)-ymin),time,ru1)
c          vbc(ib1,:,:) = -vbc(ix1,:,:)+2.0*wvel_taylor_green((xc(1)-xmin),(yv(j)-ymin),time,ru1)
          vbc(ib1,:,:) = vbc(ix1,:,:)
          wbc(ib1,:,:) =-wbc(ix1,:,:)
        ENDDO
        ENDDO
      CASE (2)
        DO K=1,NZ
        DO J=1,NY
          ubc(ix2,j,k) = uvel_taylor_green((xu(ix2)-xmin),(yc(j)-ymin),time,ru1)
c          vbc(ib2,:,:) = -vbc(ix2,:,:)+2.0*wvel_taylor_green((xc(ix2)-xmin),(yv(j)-ymin),time,ru1)
          vbc(ib2,:,:) = vbc(ix2,:,:)
          wbc(ib2,:,:) =-wbc(ix2,:,:)
        ENDDO
        ENDDO

      END SELECT

      RETURN

      END
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE BOUND082(ubc,vbc,wbc,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c  c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    ubc(nx,ny,nz),vbc(nx,ny,nz),wbc(nx,ny,nz)
c
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (1)
        ubc(ib1,:,:) =  ubc(ix1,:,:)
        vbc(ib1,:,:) = -vbc(ix1,:,:)
        wbc(ib1,:,:) = -wbc(ix1,:,:)

      CASE (2)
        ubc(ib2,:,:) = ubc(ix2,:,:)
        vbc(ib2,:,:) = vbc(ix2,:,:)
        wbc(ib2,:,:) = wbc(ix2,:,:)              

      END SELECT
c
      RETURN
      END
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE BOUND500(ubc,vbc,wbc,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c     boundary conditions for periodic boundaries
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    ubc(nx,ny,nz),vbc(nx,ny,nz),wbc(nx,ny,nz)
c
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (1) 
        ubc(ib1,:,:) = ubc(ix2,:,:)
        vbc(ib1,:,:) = vbc(ix2,:,:)
        wbc(ib1,:,:) = wbc(ix2,:,:)      

c      CASE (2)
        ubc(ib2,:,:) = ubc(ix1,:,:)
        vbc(ib2,:,:) = vbc(ix1,:,:)
        wbc(ib2,:,:) = wbc(ix1,:,:)              

      CASE (3)
        ubc(:,jb1,:) = ubc(:,jy2,:)
        vbc(:,jb1,:) = vbc(:,jy2,:)   !! jb1 = 1, jy2 = ny-1
        wbc(:,jb1,:) = wbc(:,jy2,:)      
         
c      CASE (4)
        ubc(:,jb2,:) = ubc(:,jy1,:)
        vbc(:,jb2,:) = vbc(:,jy1,:)  !! jb2 = ny, jy1 = 2
        wbc(:,jb2,:) = wbc(:,jy1,:)              
        
      CASE (5)
        ubc(:,:,kb1) = ubc(:,:,kz2)
        vbc(:,:,kb1) = vbc(:,:,kz2)
        wbc(:,:,kb1) = wbc(:,:,kz2)
        
c      CASE (6)
        ubc(:,:,kb2) = ubc(:,:,kz1)
        vbc(:,:,kb2) = vbc(:,:,kz1)
        wbc(:,:,kb2) = wbc(:,:,kz1)    
        
      END SELECT
c
      RETURN
      END
c


c----------------------------------------------------------------------
      SUBROUTINE BOUND300(ubc,vbc,wbc,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c     boundary conditions for axis boundary in cylindrical coordinates
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    ubc(nx,ny,nz),vbc(nx,ny,nz),wbc(nx,ny,nz)
      INTEGER J
c
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (1) 
        DO J=1,NY
          UBC(1,J,:) =(-UBC(IX1,JSYM(J),:)+UBC(IX1,J,:))*.5
          VBC(1,J,:) = -VBC(IX1,JSYM(J),:)
          WBC(1,J,:) =  WBC(IX1,JSYM(J),:)
        ENDDO

      CASE DEFAULT

        WRITE(6,'(A)') '*..Not a defined case!'
        
      END SELECT
c
      RETURN
      END
c

      SUBROUTINE BOUNDINOUT(UBC,VBC,WBC,UOF,VOF,WOF,XU,XC,YC,ZW,ZC,
     & ALFXDT,TIME,NX,NY,NZ,IDIV,PLANEU,PLANEV,PLANEW)

      INCLUDE 'common.h'
      INCLUDE 'mpif.h'

      INTEGER NX,NY,NZ
      INTEGER IDIV
      REAL XU(NX),XC(NX),YC(NY),ZW(NZ),ZC(NZ)
      REAL UBC(NX,NY,NZ),VBC(NX,NY,NZ),WBC(NX,NY,NZ)
      REAL UOF(NX,NY,NZ),VOF(NX,NY,NZ),WOF(NX,NY,NZ)
      REAL TIME
      REAL PLANEU(NX,NY),PLANEV(NX,NY),PLANEW(NX,NY)

      INTEGER STATUS(MPI_STATUS_SIZE)

      INTEGER I,J,K,I1,I2,N
      INTEGER IINQ,IOUQ,IINQX
      REAL AREA,AREAIJ,AREAJK,QIN(2),QOUT,COEF,UCONV,RR,QINX,QINXG
      REAL ALFXDT,WAMP,TLOCT
c
c.... Blasius
      INTEGER, SAVE :: NETA
      REAL XO,ZO,DETA,ETA1,ETA2,ALPHA,BETA,UTEMP,UETA1,UETA2,VETA1,VETA2
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: F1ETA,F2ETA,ETA
      REAL, SAVE :: TEMP
      REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: NOISE
c      REAL, DIMENSION(NX,NY,3), SAVE :: NOISE
      INTEGER, SAVE :: ISAVE=0,INS=1
c
c.....constants
c
      IINQ = IDIV 
      IOUQ = IDIV
c
c IDIV establishes the difference from the predictor and the corrector
c steps (1: predictor step; 0: corrector step)
c
C-----------------------------------------------------------------------
C-----Lower/Upper X boundary condition
C-----------------------------------------------------------------------
      IINQX = IDIV
      QINX = 0.0
      QINXG = 0.0
      IF(ITYPE(1)/=80 .AND. ITYPE(2)/=80 .AND. ITYPE(2)/=81 .AND. 
     %ITYPE(1)/=82 .AND. ITYPE(2)/=82) IINQX=0

      IF((ITYPE(1)==80 .OR. ITYPE(1)==82) .AND. IINQX==1) THEN
c Area and flow rate through the lower X boundary
        AREA = 0.0
        DO K=KZ1,KZ2
        DO J=JY1,JY2
!!!!!!          AREAJK = 1.0/(CP(K)*BP(J))
          AREAJK = RU(IB1)/(CP(K)*BP(J))
          QINX = QINX+UBC(IB1,J,K)*AREAJK
          AREA = AREA+AREAJK
        ENDDO
        ENDDO
      ENDIF

      IF((ITYPE(2)==80 .OR. ITYPE(2)==81 .OR. ITYPE(2)==82).AND. IINQX==1) THEN

        IF(MYRANK>0) THEN
          IF(ITYPE(2).NE.82)CALL MPI_RECV(QINX,1,MTYPE,MYLEFT,0,MPI_COMM_EDDY,STATUS,IERR)
        ENDIF

c Area and flow rate through the upper X boundary
        AREA = 0.0
        DO K=KZ1,KZ2
        DO J=JY1,JY2
!!!!!!          AREAJK = 1.0/(CP(K)*BP(J))
          AREAJK = RU(IX2)/(CP(K)*BP(J))    !!!!!!
          QINX = QINX-UBC(IX2,J,K)*AREAJK
          AREA = AREA+AREAJK
        ENDDO
        ENDDO

c The flow rates through the X boundaries are sent to the right processor
        IF(MYRANK<MYSIZE-1) THEN
          IF(ITYPE(2).NE.82)CALL MPI_SEND(QINX,1,MTYPE,MYRIGHT,0,MPI_COMM_EDDY,IERR)
        ENDIF

        IF(ITYPE(2).EQ.82) THEN
          COEF=1.-QINX/AREA
          UBC(IX2,:,:)=UBC(IX2,:,:)*COEF
        ENDIF

        QINXG = QINX
c        write(6,*) 'myrank=',myrank,', Qinx=',Qinx

      ENDIF

c      IF(IINQX==1) THEN
c        CALL MPI_REDUCE(QINX,QINXG,1,MTYPE,MPI_SUM,MYSIZE-1,MPI_COMM_EDDY,IERR)         
c      ENDIF
C-----------------------------------------------------------------------
C-----south boundary/inlet boundary condition
C-----------------------------------------------------------------------
      SELECT CASE(ITYPE(5)) ! ITYPE(5) is different from 0 only for MYRANK=0

c.....flat velocity profile for flatsurf dimple
      CASE (303)
         
        xo = dummyfp4
        call locate(xc,nx,xo,i)

        UBC(:,:,KB1) =  -UBC(:,:,KZ1)    !!!!!!
        VBC(:,:,KB1) =  -VBC(:,:,KZ1)    !!!!!!
        WBC(:,:,KB1) = 1.0               !!!!!!
        WBC(1:i,:,KB1) = 0.0             !!!!!!

c.....Blasius solution
      CASE (304)
         
        xo = dummyfp4
        if(ibm==0) xo = 0.0
        zo = dummyfp3 !Re=100, eta=5, x=1
        deta = 0.05
         
c...... Initialize parameters for Blasius
        if(isave==0) then
          allocate(noise(nx,ny,3))
          noise = 0.0
          isave = 1
          eta1 = (xu(nx)-xo)*sqrt(0.5/ru1/(zo+zc(1)-zmin))
          eta2 = (xu(nx)-xo)*sqrt(0.5/ru1/(zo+zw(nz)-zmin))
          neta = ceiling(max(eta1,eta2)/deta)+2
          allocate(eta(neta),f1eta(neta),f2eta(neta))
          eta = 0.0
          do i=2,neta
            eta(i) = deta*real(i-1)
          enddo
          alpha = dummyfp1
          beta = dummyfp2
          call falkner_skan(eta,f1eta,f2eta,neta,alpha,beta)
        endif
c
c...... Streamwise velocity
        WBC(:,:,kb1) = 0.0
        i1 = 1
        if(ibm>0) call locate(xc,nx,xo,i1)
        i1 = i1+1
        do i=i1,nx
          eta2 = (xc(i)-xo)*sqrt(0.5/ru1/zo)
          call locate(eta,neta,eta2,j)
          ueta1 = f2eta(j)
          ueta2 = f2eta(j+1)
          WBC(i,:,kb1) = ueta1 + (ueta2-ueta1)/deta*(eta2-eta(j))

          eta2 = (xu(i)-xo)*sqrt(0.5/ru1/(zo+zc(kb1)-zmin))
          call locate(eta,neta,eta2,j)
          veta1 = sqrt(0.5*ru1/(zo+zc(kb1)-zmin))*(eta(j)*f2eta(j)-f1eta(j))
          veta2 = sqrt(0.5*ru1/(zo+zc(kb1)-zmin))*(eta(j+1)*f2eta(j+1)-f1eta(j+1))
          UBC(i,:,kb1) = veta1 + (veta2-veta1)/deta*(eta2-eta(j))
        enddo
          
        VBC(:,:,KB1) =  -VBC(:,:,KZ1)
c
c... Add noise
        if(.false.) then
          i1 = 1
          call locate(xc,nx,2.0,i2)
          if(ibm>0) call locate(xu,nx,1.0,i1)
          i1=i1+1

          if(ins==1) then
            n = i2-i1+1
            noise = 0.0
            call add_random_noise(noise(i1:i2,:,1),n,ny,1,0.01)
            call add_random_noise(noise(i1:i2,:,2),n,ny,1,0.01)
            call add_random_noise(noise(i1:i2,:,3),n,ny,1,0.02)
          endif

          do j=jy1,jy2
          do i=i1,i2
            ubc(i,j,kb1) = ubc(i,j,kb1) + noise(i,j,1)*sin((xu(i)-xu(i1))/(xu(i2)-xu(i1))*pi)
            vbc(i,j,kb1) = vbc(i,j,kb1) + noise(i,j,2)*sin((xc(i)-xc(i1))/(xc(i2)-xc(i1))*pi)
            wbc(i,j,kb1) = wbc(i,j,kb1) + noise(i,j,3)*sin((xc(i)-xc(i1))/(xc(i2)-xc(i1))*pi)
          enddo
          enddo

          if(idiv==0) then
            ins = 1
          else
            ins = 0
          endif
        endif

c.....hybrid inlet
      CASE (305)

!	do j=1,ny
!	do i=1,nx
!        ubc(i,j,kb1) = planeU(i,j)
!        vbc(i,j,kb1) = planeV(i,j)
!       	wbc(i,j,kb1) = planeW(i,j)
!	enddo
!	enddo

        ubc(:,:,kb1) = planeU(:,:)
        vbc(:,:,kb1) = planeV(:,:)
       	wbc(:,:,kb1) = planeW(:,:)






c.....uniform inlet
      CASE (310)
        UBC(:,:,KB1) =  -UBC(:,:,KZ1)    !!!!!!
        VBC(:,:,KB1) =  -VBC(:,:,KZ1)    !!!!!!

c        WAMP=1.

        if(ibm<2) then
          tloct = 0.
        else
          tloct = mod(time,2.)
        endif

        if(    tloct>=0.0.and.tloct<=0.4) then
          WAMP = 1.0
        elseif(tloct> 0.4.and.tloct< 0.9) then
          WAMP = 0.5*(1.0+COS(2.*PI*(tloct-0.4)))
        elseif(tloct>=0.9.and.tloct<=1.1) then
          WAMP = 0.0
        elseif(tloct> 1.1.and.tloct< 1.6) then
          WAMP = 0.5*(1.0-COS(2.*PI*(tloct-1.1)))
        elseif(tloct>=1.6) then
          WAMP = 1.0
        endif

        do j=jy1,jy2
        do i=ix1,ix2
          rr = sqrt(xc(i)**2+yc(j)**2)*(1-icyl)+xc(i)*icyl
          if(rr<=0.5) then
            wbc(i,j,kb1) = wamp    !!!!!!
          else
            wbc(i,j,kb1) = 0.      !!!!!!
          endif
        enddo
        enddo

c.....parabolic profile
      CASE (320)
         
        UBC(:,:,KB1) =  -UBC(:,:,KZ1)
        VBC(:,:,KB1) =  -VBC(:,:,KZ1)

        do i=ix1,ix2
        do j=jy1,jy2
          wbc(i,j,kb1) = 1.5*(1-(2.0*xc(i)/xlen)**2.)
        enddo
        enddo

        if(isave==0) then
          isave = 1
c          allocate(noise(nx,ny,3))
        endif

        if(idiv==1 .AND. .false.) then
          noise = 0.0
          call add_random_noise(noise(:,:,1),nx,ny,1,0.01)
          call add_random_noise(noise(:,:,2),nx,ny,1,0.01)
          call add_random_noise(noise(:,:,3),nx,ny,1,0.01)
        endif

        do j=jy1,jy2
        do i=ix1,ix2
          ubc(i,j,kb1) = ubc(i,j,kb1) + noise(i,j,1)
          vbc(i,j,kb1) = vbc(i,j,kb1) + noise(i,j,2)
          wbc(i,j,kb1) = wbc(i,j,kb1) + noise(i,j,3)
        enddo
        enddo

      CASE DEFAULT

        IINQ = 0
        
      END SELECT
c
c.....mass flux at the inlet
c
      IF(IINQ==1) THEN
        AREA = 0.
        QIN = 0.0

        DO J=JY1,JY2
        DO I=IX1,IX2
          AREAIJ = RP(I)/(AP(I)*BP(J))
          QIN(1) = QIN(1)+WBC(I,J,KB1)*AREAIJ
          AREA = AREA+AREAIJ
        ENDDO 
        ENDDO

c        write(6,*) '1. qin(1)=',qin(1),areaij
        
        QIN(2) = QIN(1)/AREA

c.....send it to the last block (containing the outlet boundary) 
        IF(MYSIZE/=1) 
     &       CALL MPI_SEND(QIN,2,MTYPE,MYSIZE-1,0,MPI_COMM_EDDY,IERR)

      ENDIF
C-----------------------------------------------------------------------
C-----north boundary/outlet boundary
C-----------------------------------------------------------------------
c
c.....predictor step
c
      IF(IDIV==1) THEN

        SELECT CASE(ITYPE(6))

c.....homogeneous Neumann boundary condition
        CASE (710)

          WBC(:,:,KZ2) = WBC(:,:,KZ2-1)

c.....extrapolation from the interior
        CASE (720)

          COEF = CP(KZ2-1)/CP(KZ2)
          WBC(:,:,KZ2) = (1+COEF)*WBC(:,:,KZ2-1)-COEF*WBC(:,:,KZ2-2)

c.....convective boundary condition
        CASE (730)

c.....convective boundary condition
        CASE (731)

c     delayed

        CASE DEFAULT

          IOUQ = 0

        END SELECT

        IF(IOUQ==1) THEN

          IF(MYSIZE/=1) 
     &         CALL MPI_RECV(QIN,2,MTYPE,0,0,MPI_COMM_EDDY,STATUS,IERR)
c          write(6,*) '2. QIN(1)=',QIN(1),QINXG
          QIN(1) = QIN(1)+QINXG

c.....convective boundary condition
          IF(ITYPE(6)==730) THEN

            COEF = QIN(2)*ALFXDT*CP(KZ2)

            DO J=JY1,JY2
            DO I=IX1,IX2
                UBC(I,J,NZ ) = UOF(I,J,NZ )-COEF*(UOF(I,J,NZ )-UOF(I,J,KZ2  ))
                VBC(I,J,NZ ) = VOF(I,J,NZ )-COEF*(VOF(I,J,NZ )-VOF(I,J,KZ2  ))
                WBC(I,J,KZ2) = WOF(I,J,KZ2)-COEF*(WOF(I,J,KZ2)-WOF(I,J,KZ2-1))
            ENDDO
            ENDDO

            UBC(:,1 ,:) = UBC(:,ny-1,:) ! update gc 
            UBC(:,ny,:) = UBC(:, 2  ,:)
            VBC(:,1 ,:) = VBC(:,ny-1,:)
            VBC(:,ny,:) = VBC(:, 2  ,:)
            WBC(:,1 ,:) = WBC(:,ny-1,:)
            WBC(:,ny,:) = WBC(:, 2  ,:)


          ENDIF

c.....convective boundary condition
          IF(ITYPE(6)==731) THEN

            COEF = QIN(2)*ALFXDT*CP(KZ2)

            DO J=JY1,JY2
            DO I=IX1,IX2
                UBC(I,J,NZ ) = UOF(I,J,NZ )-COEF*(UOF(I,J,NZ )-UOF(I,J,KZ2  ))
                VBC(I,J,NZ ) = VOF(I,J,NZ )-COEF*(VOF(I,J,NZ )-VOF(I,J,KZ2  ))
                WBC(I,J,KZ2) = WOF(I,J,KZ2)-COEF*(WOF(I,J,KZ2)-WOF(I,J,KZ2-1))
            ENDDO
            ENDDO            
	WRITE(*,*)"****************************=",KZ2
            xo = dummyfp4

            call locate(xc,nx,xo,i1)
            i1 = i1+1
            ubc(1:i1,:,nz) = 0.0

            call locate(xc,nx,xo,i1)
            i1 = i1+1
            vbc(1:i1,:,nz) = 0.0
            wbc(1:i1,:,kz2:nz) = 0.0

          ENDIF
c
c.....mass balance: the outflow is evaluated
c
          AREA = 0.
          QOUT = 0.
          DO J=JY1,JY2
          DO I=IX1,IX2
            AREAIJ = RP(I)/(AP(I)*BP(J))
            QOUT = QOUT+WBC(I,J,KZ2)*AREAIJ
            AREA = AREA+AREAIJ
          ENDDO 
          ENDDO

          IF(QOUT.NE.0.) THEN !!!!!!

!!!!!!          COEF = QIN(1)/(QOUT+TINY(0.))
             COEF = QIN(1)/QOUT

             DO J=JY1,JY2
             DO I=IX1,IX2
               WBC(I,J,KZ2) = WBC(I,J,KZ2)*COEF
             ENDDO 
             ENDDO

          ELSE

             DO J=JY1,JY2
             DO I=IX1,IX2
                WBC(I,J,KZ2) = QIN(1)/AREA
             ENDDO 
             ENDDO
            
          ENDIF

        ENDIF
c
c.....corrector step
c
      ELSE

        SELECT CASE(ITYPE(6))

c.....homogeneous Neumann boundary condition
        CASE (710)

          UBC(:,:,NZ ) = UOF(:,:,KZ2)
          VBC(:,:,NZ ) = VOF(:,:,KZ2)
          WBC(:,:,KZ2) = WOF(:,:,KZ2)

c.....extrapolation from the interior
        CASE (720)

          COEF = CW(KZ2-1)/CW(KZ2)

          UBC(:,:,NZ ) = (1+COEF)*UOF(:,:,KZ2)-COEF*UOF(:,:,KZ2-1)
          VBC(:,:,NZ ) = (1+COEF)*VOF(:,:,KZ2)-COEF*VOF(:,:,KZ2-1)
          WBC(:,:,KZ2) = WOF(:,:,KZ2)

c.....convective boundary condition
        CASE (730)
 
          UBC(:,:,NZ ) = UOF(:,:,NZ )
          VBC(:,:,NZ ) = VOF(:,:,NZ )
          WBC(:,:,KZ2) = WOF(:,:,KZ2)

c.....convective boundary condition
        CASE (731)

          UBC(:,:,NZ ) = UOF(:,:,NZ )
          VBC(:,:,NZ ) = VOF(:,:,NZ )
          WBC(:,:,KZ2) = WOF(:,:,KZ2)

        CASE DEFAULT

          IOUQ = 0

        END SELECT

        SELECT CASE(ITYPE(2))

        CASE(82)

          UBC(IB1,:,:) = UOF(IX1,:,:)
          UBC(IX2,:,:) = UOF(IX2-1,:,:)

        END SELECT    

      ENDIF

c
      RETURN
      END


C---- function u_taylor_green--------------N. Beratlis-9 Dec. 2009 ---

      real function uvel_taylor_green(x,y,t,nu)
C
C     PURPOSE: Analytical u velocity for taylor green vortex. Assumes
C     0<x,y<pi
c---------------------------------------------------------------------

      implicit none
c      include 'common.h'

      real    x,y,t,nu
      
      uvel_taylor_green = exp(-2.0*nu*t)*sin(x)*cos(y)

      return

      end
c---------------------------------------------------------------------



C---- function w_taylor_green--------------N. Beratlis-9 Dec. 2009 ---

      real function wvel_taylor_green(x,y,t,nu)
C
C     PURPOSE: Analytical u velocity for taylor green vortex. Assumes
C     0<x,y<pi
c---------------------------------------------------------------------

      implicit none
c      include 'common.h'

      real    x,y,t,nu
      
      wvel_taylor_green =-exp(-2.0*nu*t)*cos(x)*sin(y)

      return

      end
c---------------------------------------------------------------------

C---- function pres_taylor_green ----------N. Beratlis-9 Dec. 2009 ---

      real function pres_taylor_green(x,y,t,nu)
C
C     PURPOSE: Analytical u velocity for taylor green vortex. Assumes
C     0<x,y<pi
c---------------------------------------------------------------------

      implicit none
c      include 'common.h'

      real    x,y,t,nu
      
      pres_taylor_green =-0.25*exp(-4.0*nu*t)*(cos(2.0*x) + cos(2.0*y))

      return

      end
c---------------------------------------------------------------------



c---- subroutine falkner_skan ---- N. Beratlis, J. Mode -14 Jul 2010 ---
C
C     PURPOSE: Compute the streamwise velocity for the Falkner Skan boundary
C     layer equation.
C
C-----------------------------------------------------------------------
      subroutine falkner_skan(eta,ya,yb,n,alpha,beta)
c
      implicit none
c
c.... Input/Output arrays
      integer n
      real    beta,alpha
      real    eta(n),ya(n),yb(n)
c
c.... Local arrays
      integer i,j
      real    deta,ybprev,error,der,ycprev
      real    k1a,k1b,k1c,k2a,k2b,k2c,k3a,k3b,k3c,k4a,k4b,k4c
      real    yc(n)
c
c.... Functions
      real    boundval_fp_falkner_skan
      ya = 0
      yb = 0
      yc = 0
c      yc(1) = boundval_fp_falkner_skan(beta)
c      write(6,*) '1.',alpha,beta
      yc(1) = 0.332
c      write(6,*) 'beta=',beta,', yc(1)=',yc(1)

      ycprev = 0
      ybprev = 0

      do i=1,5

      do j=1,n-1

        deta = eta(j+1)-eta(j)

       !k1
        k1c=-deta*( alpha*ya(j)*yc(j) + Beta*(1.0-yb(j)*yb(j)) ) 
        k1b=deta*yc(j) 
        k1a=deta*yb(j) 

       !k2
        k2c=-deta*( alpha*(ya(j)+0.5*k1a)*(yc(j)+0.5*k1c) + Beta*(1.0 - (yb(j)+0.5*k1b)*(yb(j)+0.5*k1b)) ) 
        k2b=deta*(yc(j)+0.5*k1c) 
        k2a=deta*(yb(j)+0.5*k1b) 

       !k3
        k3c=-deta*( alpha*(ya(j)+0.5*k2a)*(yc(j)+0.5*k2c) + Beta*(1.0 - (yb(j)+0.5*k2b)*(yb(j)+0.5*k2b))) 
        k3b=deta*(yc(j)+0.5*k2c) 
        k3a=deta*(yb(j)+0.5*k2b) 

       !k4
        k4c=-deta*( alpha*(ya(j)+0.5*k3a)*(yc(j)+0.5*k3c) + Beta*(1.0 - (yb(j)+0.5*k3b)*(yb(j)+0.5*k3b)))
        k4b=deta*(yc(j)+k3c)
        k4a=deta*(yb(j)+k3b)

        yc(j+1)=yc(j)+(1.0/6.0)*(k1c+2.0*k2c+2.0*k3c+k4c)
        yb(j+1)=yb(j)+(1.0/6.0)*(k1b+2.0*k2b+2.0*k3b+k4b)
        ya(j+1)=ya(j)+(1.0/6.0)*(k1a+2.0*k2a+2.0*k3a+k4a)
      enddo
        error = yb(n)-1.0
        der = (yb(n)-ybprev)/(yc(1)-ycprev)
        ycprev = yc(1)
        ybprev = yb(n)
c        write(6,*) error,der,yc(1),yc(1) - error/der,yb(n)
        yc(1) = yc(1) - error/der
        if(abs(error)<1.e-12) exit
      enddo

c      do j=1,n
c       u(j)=yb(j) !u_FK=u/U_inf
c       write(6,*) j,yb(j)
c      enddo  

      return

      end
C------------------------------------------------------------------------------


C---- function boundval_fp_falkner_skan --------- N. Beratlis - 16 Jul 2010 ---
C
C     PURPOSE: Given a value of parameter beta choose F'(eta=0) at the boundary
C
C------------------------------------------------------------------------------
      real function boundval_fp_falkner_skan(beta)
c
      implicit none
c
c.... Input/Output arrays
      real beta
c
c.... Local arrays
      real fp,f1,f2,b1,b2
c
      if(beta<0.0) then
        b1 = -0.5
        b2 = 0.0
        f1 = 0.1
        f2 = 0.4696
        fp = f1+(f2-f1)*(beta-b1)/(b2-b1)
      elseif(beta<=0.3) then
        b1 = 0.0
        b2 = 0.3
        f1 = 0.4696
        f2 = 0.77476
        fp = f1+(f2-f1)*(beta-b1)/(b2-b1)
      elseif(beta<=0.4) then
        b1 = 0.3
        b2 = 0.4
        f1 = 0.77476
        f2 = 0.8
        fp = f1+(f2-f1)*(beta-b1)/(b2-b1)
      elseif(beta<=0.5) then
        b1 = 0.4
        b2 = 0.5
        f1 = 0.8
        f2 = 0.9
        fp = f1+(f2-f1)*(beta-b1)/(b2-b1)
      elseif(beta<=0.6) then
        b1 = 0.5
        b2 = 0.6
        f1 = 0.9
        f2 = 1.0
        fp = f1+(f2-f1)*(beta-b1)/(b2-b1)
      elseif(beta<=0.8) then
        b1 = 0.6
        b2 = 0.8
        f1 = 1.0
        f2 = 1.15
        fp = f1+(f2-f1)*(beta-b1)/(b2-b1)
      elseif(beta<=1.0) then
        b1 = 0.8
        b2 = 1.0
        f1 = 1.15
        f2 = 1.23459
        fp = f1+(f2-f1)*(beta-b1)/(b2-b1)
      else
        fp = 1.23459
      endif

      boundval_fp_falkner_skan = fp

      return

      end
C------------------------------------------------------------------------------



* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
      SUBROUTINE FORTERM(UO,UBREF,UBOLD,DTM,NX,NY)
*                                                                 *
*       Computes the forcing term DPDX which represents           *
*       the overall pressure gradient in the streamwise           *
*       direction                                                 *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
c.... Input/Output arrays
      INTEGER NX,NY
      REAL    DTM,UBREF,UBOLD
      REAL    UO(NX,NY)
c
c.... Local arrays
      INTEGER I,J,K
      REAL    UB,AREAIJ,AREA
*
****  Bulk velocity
*
      IF(MYRANK.EQ.0) THEN

        IF(ITYPE(5)==500) THEN

          UB = 0.0
          AREA = 0.0
          DO J=JY1,JY2
          DO I=IX1,IX2
            AREAIJ = RP(I)/(AP(I)*BP(J))
            UB = UB + UO(I,J)*AREAIJ
            AREA = AREA+AREAIJ
          ENDDO
          ENDDO
          UB = UB/AREA
C
C.... Average streamwise presure gradient
C
          DPDZ = DPDZ + ( 2.*(UB-UBREF)/DTM - 1.*(UBOLD-UBREF)/DTM )
          write(6,*) 'dpdz=',dpdz,', Ub=',Ub,', Ubold=',Ubold
          open(81,file='forterm.out',position='append')
          write(81,'(3f20.10)')dpdz,Ub,Ubold
          close(81)
          UBOLD = UB
        ELSE
          DPDZ = 0.0
        ENDIF
      ENDIF

      CALL MPI_BCAST(DPDZ,1,MTYPE,0,MPI_COMM_EDDY,IERR)

      RETURN
      END

C---- subroutine pipe_turb_prof-----------------------N. Beratlis-20 Sep. 2011---
C
C     PURPOSE: Velocity profile normalized by bulk velocity for turbulent 
C     flow (Re=5300=Ub*D/nu).
C
C--------------------------------------------------------------------------------
      subroutine pipe_turb_prof(u,x,n)
c
      IMPLICIT NONE
c
c.... Input/Output arrays
      integer n
      real    x(n),u(n)
c
c.... Local arrays
      integer i,i1
      real    a
      integer, parameter :: n1=43
      real    x1(n1),u1(n1)
      data x1 /0.0, 0.00895963, 0.027478, 0.0436061, 0.0669008, 0.0895981, 0.120657
     &  , 0.143351, 0.162463, 0.18217, 0.202476, 0.227556, 0.245471, 0.265772, 0.280701
     &  , 0.292643, 0.307571, 0.32429, 0.339216, 0.359515, 0.375036, 0.38996, 0.401895
     &  , 0.413235, 0.42278, 0.431724, 0.438282, 0.444836, 0.453172, 0.459122, 0.463281
     &  , 0.466847, 0.470999, 0.474556, 0.478119, 0.481083, 0.483445, 0.486402, 0.489362
     &  , 0.491725, 0.494087, 0.497048, 0.5/
      data u1 /1.0, 0.99842, 0.998254, 0.996609, 0.991902, 0.9872, 0.979424, 0.970225
     & , 0.962557, 0.953384, 0.944206, 0.928988, 0.918332, 0.903156, 0.894027, 0.884924
     & , 0.874295, 0.86215, 0.847024, 0.82885, 0.810719, 0.792594, 0.772996, 0.754903
     & , 0.733828, 0.706761, 0.684213, 0.655669, 0.610616, 0.570083, 0.532564, 0.501048
     & , 0.451535, 0.406526, 0.372011, 0.333004, 0.286506, 0.237004, 0.192, 0.147001
     & , 0.102002, 0.0584973, 0.0/
c
      u1 = u1*1.31

      do i=1,n
        if(x(i)<x1(1)) then
          u(i) = 0.0
        elseif(x(i)>x1(n1)) then
          u(i) = 2.0
        else
          CALL LOCATE(x1,n1,x(i),i1)
          a = (x1(i1+1)-x(i))/(x1(i1+1)-x1(i1))
          u(i) = a*u1(i1) + (1.-a)*u1(i1+1)
        endif
      enddo

      return

      end
C--------------------------------------------------------------------------------



C---- subroutine pipe_turb_prof_wrms ------------------N. Beratlis-20 Sep. 2011---
C
C     PURPOSE: Velocity profile normalized by bulk velocity for turbulent 
C     flow (Re=5300=Ub*D/nu).
C
C--------------------------------------------------------------------------------
      subroutine pipe_turb_prof_wrms(u,x,n)
c
      IMPLICIT NONE
c
c.... Input/Output arrays
      integer n
      real    x(n),u(n)
c
c.... Local arrays
      integer i,i1
      real    a
      integer, parameter :: n1=83
      real    x1(n1),u1(n1)

      data x1 /0.0, 0.00112512, 0.00201875, 0.00314226 , 0.0196783 , 0.0261469 , 0.0270368
     &     , 0.0393571 , 0.0569811 , 0.061453, 0.0706795, 0.0864861, 0.086948, 0.100226
     &     , 0.118229 ,0.11827, 0.130677, 0.150028, 0.15132, 0.159792, 0.176428, 0.190709
     &     , 0.192018, 0.209091, 0.218483, 0.234499, 0.241754, 0.246708, 0.272264, 0.274308
     &     , 0.278446, 0.29379, 0.306047, 0.306198, 0.314423, 0.326781, 0.337763, 0.338672
     &     , 0.347354, 0.353034, 0.362346, 0.366071, 0.367025, 0.377328, 0.377962, 0.387691
     &     , 0.390259, 0.39802, 0.399421, 0.408452, 0.409458, 0.416604, 0.417597, 0.422404
     &     , 0.4234, 0.432454, 0.434031, 0.435782, 0.442866, 0.44867, 0.450035, 0.450545
     &     , 0.457728, 0.461145, 0.464847, 0.465044, 0.469991, 0.470107, 0.473103, 0.476523
     &     , 0.477932, 0.479909, 0.482765, 0.484499, 0.485376, 0.485637, 0.488493, 0.490588
     &     , 0.492801, 0.494672, 0.494898, 0.496408, 0.5/

      data u1 /0.059093143, 0.05910516, 0.059117108, 0.059132111, 0.059353157, 0.059746911
     &     , 0.059801086, 0.06055112, 0.062181942, 0.062595791, 0.063449559, 0.065395384
     &     , 0.065452206, 0.067086762, 0.069671419, 0.069677529, 0.071458927, 0.074528853
     &     , 0.074733876, 0.076078072, 0.079090292, 0.081676171, 0.081894094, 0.084735234
     &     , 0.086298031, 0.089059063, 0.090309572, 0.091163612, 0.096769857, 0.09725594
     &     , 0.098241005, 0.101891378, 0.104934148, 0.104971487, 0.107013578, 0.110633401
     &     , 0.113849966, 0.114155465, 0.117072641, 0.118980312, 0.122647658, 0.124114732
     &     , 0.12459131, 0.129741344, 0.130086219, 0.135369993, 0.136978276, 0.141836388
     &     , 0.142714189, 0.148835709, 0.149712152, 0.155939579, 0.156834352, 0.161166327
     &     , 0.162063815, 0.169657841, 0.170993211, 0.172476578, 0.1784759, 0.182547862
     &     , 0.183505092, 0.183862865, 0.1860611, 0.185113374, 0.184087576, 0.183838425
     &     , 0.177600136, 0.177454175, 0.168861507, 0.158795655, 0.153687033, 0.146522064
     &     , 0.128849966, 0.11877461, 0.113677529, 0.112160217, 0.094488798, 0.075272912
     &     , 0.0549759, 0.040973591, 0.039281263, 0.027978683, 0.0/

      do i=1,n
        if(x(i)<x1(1)) then
          u(i) = u1(1)
        elseif(x(i)>x1(n1)) then
          u(i) = u1(n1)
        else
          CALL LOCATE(x1,n1,x(i),i1)
          a = (x1(i1+1)-x(i))/(x1(i1+1)-x1(i1))
          u(i) = a*u1(i1) + (1.-a)*u1(i1+1)
        endif
      enddo

      return

      end
C--------------------------------------------------------------------------------


C---- subroutine pipe_turb_prof_vrms ------------------N. Beratlis-20 Sep. 2011---
C
C     PURPOSE: Velocity profile normalized by bulk velocity for turbulent 
C     flow (Re=5300=Ub*D/nu).
C
C--------------------------------------------------------------------------------
      subroutine pipe_turb_prof_vrms(u,x,n)
c
      IMPLICIT NONE
c
c.... Input/Output arrays
      integer n
      real    x(n),u(n)
c
c.... Local arrays
      integer i,i1
      real    a
      integer, parameter :: n1=83
      real    x1(n1),u1(n1)

      data x1 /0.0, 0.00112512, 0.00201875, 0.00314226 , 0.0196783 , 0.0261469 , 0.0270368
     &     , 0.0393571 , 0.0569811 , 0.061453, 0.0706795, 0.0864861, 0.086948, 0.100226
     &     , 0.118229 ,0.11827, 0.130677, 0.150028, 0.15132, 0.159792, 0.176428, 0.190709
     &     , 0.192018, 0.209091, 0.218483, 0.234499, 0.241754, 0.246708, 0.272264, 0.274308
     &     , 0.278446, 0.29379, 0.306047, 0.306198, 0.314423, 0.326781, 0.337763, 0.338672
     &     , 0.347354, 0.353034, 0.362346, 0.366071, 0.367025, 0.377328, 0.377962, 0.387691
     &     , 0.390259, 0.39802, 0.399421, 0.408452, 0.409458, 0.416604, 0.417597, 0.422404
     &     , 0.4234, 0.432454, 0.434031, 0.435782, 0.442866, 0.44867, 0.450035, 0.450545
     &     , 0.457728, 0.461145, 0.464847, 0.465044, 0.469991, 0.470107, 0.473103, 0.476523
     &     , 0.477932, 0.479909, 0.482765, 0.484499, 0.485376, 0.485637, 0.488493, 0.490588
     &     , 0.492801, 0.494672, 0.494898, 0.496408, 0.5/

      data u1 /0.043922811, 0.04416687, 0.04416558, 0.044163883, 0.044139579, 0.044130007
     &     , 0.044141073, 0.044294229, 0.044513238, 0.044568839, 0.04517685, 0.046218534
     &     , 0.046248948, 0.047477665, 0.049143585, 0.049147386, 0.049991921, 0.051309097
     &     , 0.051427223, 0.05220224, 0.053724033, 0.054990428, 0.055106449, 0.056620502
     &     , 0.05745336, 0.058873659, 0.059516972, 0.059973523, 0.062328581, 0.062516904
     &     , 0.062898303, 0.064639443, 0.066030278, 0.066047386, 0.066819688, 0.067980312
     &     , 0.068619145, 0.068671419, 0.069176511, 0.069380855, 0.069715547, 0.069849287
     &     , 0.06988391, 0.069868975, 0.069867617, 0.06985336, 0.069849966, 0.069441276
     &     , 0.069367957, 0.068892736, 0.068839783, 0.068422946, 0.068365241, 0.068084861
     &     , 0.067955193, 0.066779566, 0.066574745, 0.066347318, 0.064873388, 0.063665716
     &     , 0.063381806, 0.063200883, 0.060652003, 0.059439308, 0.056969586, 0.056837746
     &     , 0.053537203, 0.053440462, 0.050936592, 0.04807814, 0.046900204, 0.044238086
     &     , 0.04039131, 0.038056823, 0.036500747, 0.036037475, 0.030968907, 0.02725112
     &     , 0.021955804, 0.017479769, 0.016938764, 0.012720502, 0.0/

      do i=1,n
        if(x(i)<x1(1)) then
          u(i) = u1(1)
        elseif(x(i)>x1(n1)) then
          u(i) = u1(n1)
        else
          CALL LOCATE(x1,n1,x(i),i1)
          a = (x1(i1+1)-x(i))/(x1(i1+1)-x1(i1))
          u(i) = a*u1(i1) + (1.-a)*u1(i1+1)
        endif
      enddo

      return

      end
C--------------------------------------------------------------------------------



C---- subroutine pipe_turb_prof_urms ------------------N. Beratlis-20 Sep. 2011---
C
C     PURPOSE: Velocity profile normalized by bulk velocity for turbulent 
C     flow (Re=5300=Ub*D/nu).
C
C--------------------------------------------------------------------------------
      subroutine pipe_turb_prof_urms(u,x,n)
c
      IMPLICIT NONE
c
c.... Input/Output arrays
      integer n
      real    x(n),u(n)
c
c.... Local arrays
      integer i,i1
      real    a
      integer, parameter :: n1=83
      real    x1(n1),u1(n1)

      data x1 /0.0, 0.00112512, 0.00201875, 0.00314226 , 0.0196783 , 0.0261469 , 0.0270368
     &     , 0.0393571 , 0.0569811 , 0.061453, 0.0706795, 0.0864861, 0.086948, 0.100226
     &     , 0.118229 ,0.11827, 0.130677, 0.150028, 0.15132, 0.159792, 0.176428, 0.190709
     &     , 0.192018, 0.209091, 0.218483, 0.234499, 0.241754, 0.246708, 0.272264, 0.274308
     &     , 0.278446, 0.29379, 0.306047, 0.306198, 0.314423, 0.326781, 0.337763, 0.338672
     &     , 0.347354, 0.353034, 0.362346, 0.366071, 0.367025, 0.377328, 0.377962, 0.387691
     &     , 0.390259, 0.39802, 0.399421, 0.408452, 0.409458, 0.416604, 0.417597, 0.422404
     &     , 0.4234, 0.432454, 0.434031, 0.435782, 0.442866, 0.44867, 0.450035, 0.450545
     &     , 0.457728, 0.461145, 0.464847, 0.465044, 0.469991, 0.470107, 0.473103, 0.476523
     &     , 0.477932, 0.479909, 0.482765, 0.484499, 0.485376, 0.485637, 0.488493, 0.490588
     &     , 0.492801, 0.494672, 0.494898, 0.496408, 0.5/

      data u1 /0.044185743, 0.044175628, 0.04416558, 0.044152885, 0.043966327, 0.043893347
     &     , 0.043883367, 0.044067142, 0.044330007, 0.044472166, 0.044765513, 0.045268092
     &     , 0.045285268, 0.045778955, 0.046448269, 0.046450373, 0.047076103, 0.048052071
     &     , 0.048117244, 0.048615547, 0.049594094, 0.050434148, 0.050511134, 0.051373523
     &     , 0.051847997, 0.052656959, 0.053138221, 0.053466802, 0.055161982, 0.055297556
     &     , 0.055419416, 0.055871351, 0.056232383, 0.05623442, 0.056346029, 0.056513714
     &     , 0.056662729, 0.056675085, 0.056482349, 0.056356212, 0.056149491, 0.055792804
     &     , 0.055701426, 0.054714936, 0.054654175, 0.053092532, 0.052680312, 0.051434691
     &     , 0.051046096, 0.048542634, 0.048263815, 0.046282688, 0.046007536, 0.044134216
     &     , 0.043745893, 0.040218058, 0.039603394, 0.038778955, 0.035443856, 0.032711202
     &     , 0.031850238, 0.031528649, 0.026997013, 0.024841005, 0.022505771, 0.022381127
     &     , 0.018851324, 0.018768771, 0.016631093, 0.014190699, 0.013185064, 0.011774542
     &     , 0.009736388, 0.008499457, 0.007873727, 0.0077148, 0.005976544, 0.004701494
     &     , 0.003354705, 0.002216273, 0.002078676, 0.001159735, 0.0/

      do i=1,n
        if(x(i)<x1(1)) then
          u(i) = u1(1)
        elseif(x(i)>x1(n1)) then
          u(i) = u1(n1)
        else
          CALL LOCATE(x1,n1,x(i),i1)
          a = (x1(i1+1)-x(i))/(x1(i1+1)-x1(i1))
          u(i) = a*u1(i1) + (1.-a)*u1(i1+1)
        endif
      enddo

      return

      end
C--------------------------------------------------------------------------------



C---- subroutine turb_pipe_vel -----------------------N. Beratlis-21 Sep. 2011---
C
C     PURPOSE: Turbulent pipe velocity
C
C--------------------------------------------------------------------------------
      subroutine turb_pipe_vel(uo,vo,wo,xu,xc,nx,ny,nz)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz
      real    xu(nx),xc(nx)
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k,i1,ni
      real    wrk1d(nx),uprms(nx),vprms(nx),wprms(nx)

      uo = 0.0
      vo = 0.0
      wo = 0.0

      call locate(xc,nx,0.2,i1)
      ni = nx-i1+1
c      call pipe_turb_prof(wrk1D,xc,nx)
c      do i=1,nx
c        wo(i,:,:) = wo(i,:,:) + wrk1D(i)
c      enddo
      call add_random_noise(uo(i1:nx,:,:),ni,ny,nz,1.0)
      call add_random_noise(vo(i1:nx,:,:),ni,ny,nz,1.0)
      call add_random_noise(wo(i1:nx,:,:),ni,ny,nz,1.0)
c
c.... Compute rms of noise
      uprms = 0.0
      vprms = 0.0
      wprms = 0.0
      do i=i1,nx
        uprms(i) = uprms(i) + sum(uo(i,2:ny-1,2:nz-1)**2.0)/real(ny-2)/real(nz-2)
        vprms(i) = vprms(i) + sum(vo(i,2:ny-1,2:nz-1)**2.0)/real(ny-2)/real(nz-2)
        wprms(i) = wprms(i) + sum((wo(i,2:ny-1,2:nz-1))**2.0)/real(ny-2)/real(nz-2)
      enddo

      wrk1D = uprms/real(mysize)
      CALL MPI_ALLREDUCE(wrk1D,uprms,nx,mtype,mpi_sum,mpi_comm_eddy,ierr)
      uprms = sqrt(uprms)
      wrk1D = vprms/real(mysize)
      CALL MPI_ALLREDUCE(wrk1D,vprms,nx,mtype,mpi_sum,mpi_comm_eddy,ierr)
      vprms = sqrt(vprms)
      wrk1D = wprms/real(mysize)
      CALL MPI_ALLREDUCE(wrk1D,wprms,nx,mtype,mpi_sum,mpi_comm_eddy,ierr)
      wprms = sqrt(wprms)
c
c.... Scale noise
      call pipe_turb_prof_urms(wrk1D,xu,nx)
      do i=i1,nx-1
        uo(i,:,:) = uo(i,:,:)*wrk1D(i)/uprms(i)
      enddo
          
      call pipe_turb_prof_vrms(wrk1D,xc,nx)
      do i=i1,nx-1
        vo(i,:,:) = vo(i,:,:)*wrk1D(i)/vprms(i)
      enddo

      call pipe_turb_prof_wrms(wrk1D,xc,nx)
      do i=i1,nx-1
        wo(i,:,:) = wo(i,:,:)*wrk1D(i)/wprms(i)
      enddo

      call pipe_turb_prof(wrk1D,xc,nx)
      do i=i1,nx
        wo(i,:,:) = wo(i,:,:) + wrk1D(i)
      enddo
c          
      uprms = 0.0
      vprms = 0.0
      wprms = 0.0
c      call pipe_turb_prof(wrk1D,xc,nx)
      do i=i1,nx
        uprms(i) = uprms(i) + sum(uo(i,2:ny-1,2:nz-1)**2.0)/real(ny-2)/real(nz-2)
        vprms(i) = vprms(i) + sum(vo(i,2:ny-1,2:nz-1)**2.0)/real(ny-2)/real(nz-2)
        wprms(i) = wprms(i) + sum((wo(i,2:ny-1,2:nz-1)-wrk1d(i))**2.0)/real(ny-2)/real(nz-2)
      enddo

      wrk1D = uprms/real(mysize)
      CALL MPI_ALLREDUCE(wrk1D,uprms,nx,mtype,mpi_sum,mpi_comm_eddy,ierr)
      uprms = sqrt(uprms)
      wrk1D = vprms/real(mysize)
      CALL MPI_ALLREDUCE(wrk1D,vprms,nx,mtype,mpi_sum,mpi_comm_eddy,ierr)
      vprms = sqrt(vprms)
      wrk1D = wprms/real(mysize)
      CALL MPI_ALLREDUCE(wrk1D,wprms,nx,mtype,mpi_sum,mpi_comm_eddy,ierr)
      wprms = sqrt(wprms)

      uo(nx,:,:) = -vo(nx-1,:,:)
      vo(nx,:,:) = -vo(nx-1,:,:)
      wo(nx,:,:) = -wo(nx-1,:,:)

c      do i=1,nx
c        write(6,'(I4,10(1x,F10.5))') i,xc(i),wprms(i),vprms(i),uprms(i),wo(i,2,2),vo(i,2,2),uo(i,2,2)
c      enddo

      return

      end
C--------------------------------------------------------------------------------



C---- subroutine turb_chan_vel -----------------------N. Beratlis-21 Sep. 2011---
C
C     PURPOSE: Turbulent pipe velocity
C
C--------------------------------------------------------------------------------
      subroutine turb_chan_vel(uo,vo,wo,xu,xc,nx,ny,nz,flag)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,flag
      real    xu(nx),xc(nx)
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k
      real    ub,tmp
      real    wrk1d(nx),uprms(nx),vprms(nx),wprms(nx)

      uo = 0.0
      vo = 0.0
      wo = 0.0
      call add_random_noise(uo(:,:,:),nx,ny,nz,1.0)
      call add_random_noise(vo(:,:,:),nx,ny,nz,1.0)
      call add_random_noise(wo(:,:,:),nx,ny,nz,1.0)
c
c.... Compute rms of noise
      uprms = 0.0
      vprms = 0.0
      wprms = 0.0
      do i=1,nx
        uprms(i) = uprms(i) + sum(uo(i,2:ny-1,2:nz-1)**2.0)/real(ny-2)/real(nz-2)
        vprms(i) = vprms(i) + sum(vo(i,2:ny-1,2:nz-1)**2.0)/real(ny-2)/real(nz-2)
        wprms(i) = wprms(i) + sum((wo(i,2:ny-1,2:nz-1))**2.0)/real(ny-2)/real(nz-2)
      enddo

      wrk1D = uprms/real(mysize)
      CALL MPI_ALLREDUCE(wrk1D,uprms,nx,mtype,mpi_sum,mpi_comm_eddy,ierr)
      uprms = sqrt(uprms)
      wrk1D = vprms/real(mysize)
      CALL MPI_ALLREDUCE(wrk1D,vprms,nx,mtype,mpi_sum,mpi_comm_eddy,ierr)
      vprms = sqrt(vprms)
      wrk1D = wprms/real(mysize)
      CALL MPI_ALLREDUCE(wrk1D,wprms,nx,mtype,mpi_sum,mpi_comm_eddy,ierr)
      wprms = sqrt(wprms)
c
c.... Scale noise
      call chan_turb_prof_urms(wrk1D,xu,nx)
      do i=1,nx-1
        uo(i,:,:) = uo(i,:,:)*wrk1D(i)/uprms(i)
      enddo
          
      call chan_turb_prof_vrms(wrk1D,xc,nx)
      do i=1,nx-1
        vo(i,:,:) = vo(i,:,:)*wrk1D(i)/vprms(i)
      enddo

      call chan_turb_prof_wrms(wrk1D,xc,nx)
      do i=1,nx-1
        wo(i,:,:) = wo(i,:,:)*wrk1D(i)/wprms(i)
      enddo

      call chan_turb_prof(wrk1D,xc,nx)
      do i=1,nx
        wo(i,:,:) = wo(i,:,:) + wrk1D(i)
      enddo
c          
      uprms = 0.0
      vprms = 0.0
      wprms = 0.0
      do i=1,nx
        uprms(i) = uprms(i) + sum(uo(i,2:ny-1,2:nz-1)**2.0)/real(ny-2)/real(nz-2)
        vprms(i) = vprms(i) + sum(vo(i,2:ny-1,2:nz-1)**2.0)/real(ny-2)/real(nz-2)
        wprms(i) = wprms(i) + sum((wo(i,2:ny-1,2:nz-1)-wrk1d(i))**2.0)/real(ny-2)/real(nz-2)
      enddo

      wrk1D = uprms/real(mysize)
      CALL MPI_ALLREDUCE(wrk1D,uprms,nx,mtype,mpi_sum,mpi_comm_eddy,ierr)
      uprms = sqrt(uprms)
      wrk1D = vprms/real(mysize)
      CALL MPI_ALLREDUCE(wrk1D,vprms,nx,mtype,mpi_sum,mpi_comm_eddy,ierr)
      vprms = sqrt(vprms)
      wrk1D = wprms/real(mysize)
      CALL MPI_ALLREDUCE(wrk1D,wprms,nx,mtype,mpi_sum,mpi_comm_eddy,ierr)
      wprms = sqrt(wprms)

c      uo(nx,:,:) =-uo(nx-1,:,:)
      vo(nx,:,:) =-vo(nx-1,:,:)
      wo(nx,:,:) =-wo(nx-1,:,:)
      vo(1,:,:) = -vo(2,:,:)
      wo(1,:,:) = -wo(2,:,:)

      if(flag==0) then
        if(myrank.eq.0) call calc_ubulk(wo(:,:,1),nx,ny,1,ub)
        CALL MPI_BCAST(ub,1,mtype,0,mpi_comm_eddy,ierr)
        wo = wo/ub
        vo = vo/ub
        uo = uo/ub
      endif

      return

      end
C--------------------------------------------------------------------------------


C---- subroutine chan_turb_prof-----------------------N. Beratlis-20 Sep. 2011---
C
C     PURPOSE: Velocity profile normalized by bulk velocity for turbulent 
C     flow (Re=360=Utau*H/nu).
C
C--------------------------------------------------------------------------------
      subroutine chan_turb_prof(u,x,n)
c
      IMPLICIT NONE
c
c.... Input/Output arrays
      integer n
      real    x(n),u(n)
c
c.... Local arrays
      integer i,i1
      real    a
      integer, parameter :: n1=66
      real    x1(n1),u1(n1)
      data x1 /-1.0, -0.994266333, -0.993023222, -0.991510556, -0.989615111, -0.987888889, -0.985496278
     &     , -0.983085667, -0.979420611, -0.975616444, -0.972309389, -0.967365722, -0.961743222
     &     , -0.956094222, -0.947706722, -0.939985556, -0.923425, -0.910708333, -0.902808333
     &     , -0.893083333, -0.87986, -0.865002222, -0.844232222, -0.821212778, -0.783615556
     &     , -0.713408889, -0.632307222, -0.550246111, -0.44694, -0.330633333, -0.202627778
     &     , -0.118105556, -0.045072222, 0.045072222, 0.118105556, 0.202627778, 0.330633333
     &     , 0.44694 ,0.550246111, 0.632307222, 0.713408889, 0.783615556, 0.821212778, 0.844232222
     &     , 0.865002222, 0.87986, 0.893083333, 0.902808333, 0.910708333, 0.923425, 0.939985556
     &     , 0.947706722, 0.956094222, 0.961743222, 0.967365722, 0.972309389, 0.975616444
     &     , 0.979420611, 0.983085667, 0.985496278, 0.987888889, 0.989615111, 0.991510556
     &     , 0.993023222, 0.994266333, 1.0/

      data u1 /0.0, 1.04301, 1.30777, 1.57252, 1.86995, 2.19895, 2.59365, 2.98779
     &     , 3.64342, 4.23334, 4.78979, 5.54246, 6.32759, 7.14484, 8.12575
     &     , 8.97558, 10.3488, 11.199, 11.6242, 12.1149, 12.5408, 12.9993
     &     , 13.5235, 13.9499, 14.4751, 15.1977, 15.8219, 16.4125, 16.9706
     &     , 17.4632, 17.8253, 18.0228, 18.0897, 18.0897, 18.0228, 17.8253
     &     , 17.4632, 16.9706, 16.4125, 15.8219, 15.1977, 14.4751, 13.9499
     &     , 13.5235, 12.9993, 12.5408, 12.1149, 11.6242, 11.199, 10.3488
     &     , 8.97558, 8.12575, 7.14484, 6.32759, 5.54246, 4.78979, 4.23334
     &     , 3.64342, 2.98779, 2.59365, 2.19895, 1.86995, 1.57252, 1.30777
     &     , 1.04301, 0.0/
c
      do i=1,n
        if(x(i)<x1(1)) then
          u(i) = 0.0
        elseif(x(i)>x1(n1)) then
          u(i) = 2.0
        else
          CALL LOCATE(x1,n1,x(i),i1)
          a = (x1(i1+1)-x(i))/(x1(i1+1)-x1(i1))
          u(i) = a*u1(i1) + (1.-a)*u1(i1+1)
        endif
      enddo

      return

      end
C--------------------------------------------------------------------------------



C---- subroutine chan_turn_prof_wrms ------------------N. Beratlis-20 Sep. 2011---
C
C     PURPOSE: Streamwise r.m.s velocity fluctuations normalized by utau for turbulent 
C     flow (Re=360=ut*H/nu).
C
C--------------------------------------------------------------------------------
      subroutine chan_turb_prof_wrms(u,x,n)
c
      IMPLICIT NONE
c
c.... Input/Output arrays
      integer n
      real    x(n),u(n)
c
c.... Local arrays
      integer i,i1
      real    a
      integer, parameter :: n1=92
      real    x1(n1),u1(n1)

      data x1 /-1.0, -0.986017, -0.979919, -0.97383, -0.973709, -0.965932, -0.961959
     &     , -0.95463, -0.947273, -0.941804, -0.937193, -0.931201, -0.918981
     &     , -0.917303, -0.906845, -0.896867, -0.886209, -0.868142, -0.84159
     &     , -0.835832, -0.831407, -0.78951, -0.778027, -0.762725, -0.745922
     &     , -0.72633, -0.700687, -0.671034, -0.653303, -0.627664, -0.598309
     &     , -0.544263, -0.540439, -0.454423, -0.424397, -0.403313, -0.308806
     &     , -0.302773, -0.207105, -0.158541, -0.157953, -0.10116, -0.0729155
     &     , -0.0646677, -0.0316345, -0.0138334, 0.0138334, 0.0316345, 0.0646677
     &     , 0.0729155, 0.10116, 0.157953, 0.158541, 0.207105, 0.302773, 0.308806
     &     , 0.403313, 0.424397, 0.454423, 0.540439, 0.544263, 0.598309, 0.627664
     &     , 0.653303, 0.671034, 0.700687, 0.72633, 0.745922, 0.762725, 0.778027
     &     , 0.78951, 0.831407, 0.835832, 0.84159, 0.868142, 0.886209, 0.896867
     &     , 0.906845, 0.917303, 0.918981, 0.931201, 0.937193, 0.941804, 0.947273
     &     , 0.95463, 0.961959, 0.965932, 0.973709, 0.97383, 0.979919, 0.986017, 1.0/

      data u1 /0.0, 0.779992, 1.15822, 1.53588, 1.54211, 1.94195, 2.07335, 2.31581
     &     , 2.48864, 2.59715, 2.62858, 2.66941, 2.68937, 2.68713, 2.67314
     &     , 2.62441, 2.57235, 2.46792, 2.31444, 2.28377, 2.2602, 2.03704
     &     , 1.97587, 1.91498, 1.84812, 1.77016, 1.69898, 1.61667, 1.58312
     &     , 1.53461, 1.47907, 1.4044, 1.39912, 1.28028, 1.23846, 1.20909
     &     , 1.07745, 1.06966, 0.946089, 0.883361, 0.882603, 0.830888, 0.805169
     &     , 0.802263, 0.790622, 0.784349, 0.784349, 0.790622, 0.802263, 0.805169
     &     , 0.830888, 0.882603, 0.883361, 0.946089, 1.06966, 1.07745, 1.20909
     &     , 1.23846, 1.28028, 1.39912, 1.4044, 1.47907, 1.53461, 1.58312
     &     , 1.61667, 1.69898, 1.77016, 1.84812, 1.91498, 1.97587, 2.03704
     &     , 2.2602, 2.28377, 2.31444, 2.46792, 2.57235, 2.62441, 2.67314
     &     , 2.68713, 2.68937, 2.66941, 2.62858, 2.59715, 2.48864, 2.31581
     &     , 2.07335, 1.94195, 1.54211, 1.53588, 1.15822, 0.779992, 0.0/

      do i=1,n
        if(x(i)<x1(1)) then
          u(i) = u1(1)
        elseif(x(i)>x1(n1)) then
          u(i) = u1(n1)
        else
          CALL LOCATE(x1,n1,x(i),i1)
          a = (x1(i1+1)-x(i))/(x1(i1+1)-x1(i1))
          u(i) = a*u1(i1) + (1.-a)*u1(i1+1)
        endif
      enddo

      return

      end
C--------------------------------------------------------------------------------


C---- subroutine chan_turb_prof_vrms ------------------N. Beratlis-20 Sep. 2011---
C
C     PURPOSE: Spanwise r.m.s velocity fluctuations normalized by utau for turbulent 
C     flow (Re=360=ut*H/nu).
C
C--------------------------------------------------------------------------------
      subroutine chan_turb_prof_vrms(u,x,n)
c
      IMPLICIT NONE
c
c.... Input/Output arrays
      integer n
      real    x(n),u(n)
c
c.... Local arrays
      integer i,i1
      real    a
      integer, parameter :: n1=92
      real    x1(n1),u1(n1)

      data x1 /-1.0, -0.986017, -0.979919, -0.97383, -0.973709, -0.965932, -0.961959
     &     , -0.95463, -0.947273, -0.941804, -0.937193, -0.931201, -0.918981
     &     , -0.917303, -0.906845, -0.896867, -0.886209, -0.868142, -0.84159
     &     , -0.835832, -0.831407, -0.78951, -0.778027, -0.762725, -0.745922
     &     , -0.72633, -0.700687, -0.671034, -0.653303, -0.627664, -0.598309
     &     , -0.544263, -0.540439, -0.454423, -0.424397, -0.403313, -0.308806
     &     , -0.302773, -0.207105, -0.158541, -0.157953, -0.10116, -0.0729155
     &     , -0.0646677, -0.0316345, -0.0138334, 0.0138334, 0.0316345, 0.0646677
     &     , 0.0729155, 0.10116, 0.157953, 0.158541, 0.207105, 0.302773, 0.308806
     &     , 0.403313, 0.424397, 0.454423, 0.540439, 0.544263, 0.598309, 0.627664
     &     , 0.653303, 0.671034, 0.700687, 0.72633, 0.745922, 0.762725, 0.778027
     &     , 0.78951, 0.831407, 0.835832, 0.84159, 0.868142, 0.886209, 0.896867
     &     , 0.906845, 0.917303, 0.918981, 0.931201, 0.937193, 0.941804, 0.947273
     &     , 0.95463, 0.961959, 0.965932, 0.973709, 0.97383, 0.979919, 0.986017, 1.0/

      data u1 /0.0, 0.326721, 0.409945, 0.493042, 0.494697, 0.600827, 0.655037, 0.708502
     &     , 0.762172, 0.802067, 0.835703, 0.856543, 0.899049, 0.904887, 0.941262
     &     , 0.975968, 0.988479, 1.00969, 1.04086, 1.04762, 1.04927, 1.06489
     &     , 1.06918, 1.07489, 1.06542, 1.05438, 1.03994, 1.02323, 1.01324
     &     , 0.995577, 0.97535, 0.938109, 0.935474, 0.863804, 0.838787, 0.82122
     &     , 0.756158, 0.752005, 0.686145, 0.657908, 0.657566, 0.624545, 0.617669
     &     , 0.615661, 0.60762, 0.603287, 0.603287, 0.60762, 0.615661, 0.617669
     &     , 0.624545, 0.657566, 0.657908, 0.686145, 0.752005, 0.756158, 0.82122
     &     , 0.838787, 0.863804, 0.935474, 0.938109, 0.97535, 0.995577, 1.01324
     &     , 1.02323, 1.03994, 1.05438, 1.06542, 1.07489, 1.06918, 1.06489
     &     , 1.04927, 1.04762, 1.04086, 1.00969, 0.988479, 0.975968, 0.941262
     &     , 0.904887, 0.899049, 0.856543, 0.835703, 0.802067, 0.762172, 0.708502
     &     , 0.655037, 0.600827, 0.494697, 0.493042, 0.409945, 0.326721, 0.0/

      do i=1,n
        if(x(i)<x1(1)) then
          u(i) = u1(1)
        elseif(x(i)>x1(n1)) then
          u(i) = u1(n1)
        else
          CALL LOCATE(x1,n1,x(i),i1)
          a = (x1(i1+1)-x(i))/(x1(i1+1)-x1(i1))
          u(i) = a*u1(i1) + (1.-a)*u1(i1+1)
        endif
      enddo

      return

      end
C--------------------------------------------------------------------------------


C---- subroutine chan_turb_prof_urms ------------------N. Beratlis-20 Sep. 2011---
C
C     PURPOSE: Wall normal r.m.s velocity fluctuations normalized by utau for turbulent 
C     flow (Re=360=ut*H/nu).
C
C--------------------------------------------------------------------------------
      subroutine chan_turb_prof_urms(u,x,n)
c
      IMPLICIT NONE
c
c.... Input/Output arrays
      integer n
      real    x(n),u(n)
c
c.... Local arrays
      integer i,i1
      real    a
      integer, parameter :: n1=92
      real    x1(n1),u1(n1)

      data x1 /-1.0, -0.986017, -0.979919, -0.97383, -0.973709, -0.965932, -0.961959
     &     , -0.95463, -0.947273, -0.941804, -0.937193, -0.931201, -0.918981
     &     , -0.917303, -0.906845, -0.896867, -0.886209, -0.868142, -0.84159
     &     , -0.835832, -0.831407, -0.78951, -0.778027, -0.762725, -0.745922
     &     , -0.72633, -0.700687, -0.671034, -0.653303, -0.627664, -0.598309
     &     , -0.544263, -0.540439, -0.454423, -0.424397, -0.403313, -0.308806
     &     , -0.302773, -0.207105, -0.158541, -0.157953, -0.10116, -0.0729155
     &     , -0.0646677, -0.0316345, -0.0138334, 0.0138334, 0.0316345, 0.0646677
     &     , 0.0729155, 0.10116, 0.157953, 0.158541, 0.207105, 0.302773, 0.308806
     &     , 0.403313, 0.424397, 0.454423, 0.540439, 0.544263, 0.598309, 0.627664
     &     , 0.653303, 0.671034, 0.700687, 0.72633, 0.745922, 0.762725, 0.778027
     &     , 0.78951, 0.831407, 0.835832, 0.84159, 0.868142, 0.886209, 0.896867
     &     , 0.906845, 0.917303, 0.918981, 0.931201, 0.937193, 0.941804, 0.947273
     &     , 0.95463, 0.961959, 0.965932, 0.973709, 0.97383, 0.979919, 0.986017, 1.0/

      data u1 /0.0, 0.0373709, 0.0573709, 0.0773709, 0.0973709, 0.132972, 0.155864, 0.198101
     &     , 0.2405, 0.272017, 0.29859, 0.333116, 0.403539, 0.413211, 0.455859, 0.496548
     &     , 0.540013, 0.613691, 0.680225, 0.694654, 0.705741, 0.7736, 0.786175, 0.802931
     &     , 0.821331, 0.824579, 0.828829, 0.825205, 0.823038, 0.819905, 0.808227, 0.786727
     &     , 0.784885, 0.743446, 0.728981, 0.720361, 0.681723, 0.679256, 0.646093, 0.629259
     &     , 0.629151, 0.618734, 0.613553, 0.61204, 0.605981, 0.602716, 0.602716, 0.605981
     &     , 0.61204, 0.613553, 0.618734, 0.629151, 0.629259, 0.646093, 0.679256, 0.681723
     &     , 0.720361, 0.728981, 0.743446, 0.784885, 0.786727, 0.808227, 0.819905, 0.823038
     &     , 0.825205, 0.828829, 0.824579, 0.821331, 0.802931, 0.786175, 0.7736, 0.705741
     &     , 0.694654, 0.680225, 0.613691, 0.540013, 0.496548, 0.455859, 0.413211, 0.403539
     &     , 0.333116, 0.29859, 0.272017, 0.2405, 0.198101, 0.155864, 0.132972, 0.0973709
     &     , 0.0773709, 0.0573709, 0.0373709, 0.0/

      do i=1,n
        if(x(i)<x1(1)) then
          u(i) = u1(1)
        elseif(x(i)>x1(n1)) then
          u(i) = u1(n1)
        else
          CALL LOCATE(x1,n1,x(i),i1)
          a = (x1(i1+1)-x(i))/(x1(i1+1)-x1(i1))
          u(i) = a*u1(i1) + (1.-a)*u1(i1+1)
        endif
      enddo

      return

      end
C--------------------------------------------------------------------------------

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*                 boundary conditions
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      SUBROUTINE BOUNDARY_DENS(DENS,XC,YC,NX,NY,NZ)
c
      INCLUDE 'common.h'
c
      INTEGER NX,NY,NZ,I,J
      REAL XC(NX),DENS(NX,NY,NZ),YC(NY)
c
      INTEGER IBND
c
c.....wall boundary conditions
c
      
      DO IBND=1,6

        SELECT CASE (ITYPE(IBND))

        CASE(50)
          CALL BOUND050D(dens,NX,NY,NZ,IBND)           
c
c.....propagation of waves in a linearly stratified fluid
c Adds case(80), for some reason, setting boundary=82 for velocity
c give high div.
        CASE ( 80)
          CALL BOUND082D(DENS,XC,YC,NX,NY,NZ,IBND)

        CASE ( 82)
          CALL BOUND082D(DENS,XC,YC,NX,NY,NZ,IBND)
	
	CASE (300)
	  CALL BOUND300D(DENS,NX,NY,NZ,IBND)

        END SELECT

      ENDDO
c
c.....periodic boundary conditions
c
      DO IBND=1,5,2

        IF(ITYPE(IBND)==500) CALL BOUND500D(DENS,NX,NY,NZ,IBND)

      ENDDO
c
      RETURN
      END

c----------------------------------------------------------------------
      SUBROUTINE BOUND050D(dens,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c     boundary conditions for density
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    dens(nx,ny,nz)
c
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (1) 
        dens(ib1,:,:) = 2.0*( 0.5) - dens(ib1+1,:,:)

      CASE (2)
        dens(ib2,:,:) = 2.0*(-0.5) - dens(ib2-1,:,:)
 
C       CASE (3)
C         ubc(:,jb1,:) = -ubc(:,jy1,:)
C         vbc(:,jb1,:) =  0.
C         wbc(:,jb1,:) = -wbc(:,jy1,:)      
         
C       CASE (4)
C         ubc(:,jb2,:) = -ubc(:,jy2,:)
C         vbc(:,jy2,:) =  0.
C         wbc(:,jb2,:) = -wbc(:,jy2,:)              
        
C       CASE (5)
C         ubc(:,:,kb1) = -ubc(:,:,kz1)
C         vbc(:,:,kb1) = -vbc(:,:,kz1)
C         wbc(:,:,kb1) =  0.
        
C       CASE (6)
C         ubc(:,:,kb2) = -ubc(:,:,kz2)
C         vbc(:,:,kb2) = -vbc(:,:,kz2)
c        wbc(:,:,kz2) =  0.              
        
      END SELECT
c
      RETURN
      END

c

c----------------------------------------------------------------------
      SUBROUTINE BOUND082D(dens,xc,yc,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c     boundary conditions for non-slip wall
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER ibnd
      INTEGER nx,ny,nz,j
      REAL    xc(nx),dens(nx,ny,nz),yc(ny)
c
*
**** Main loop
*
	
      SELECT CASE(IBND)

      CASE (1)
!         dens(ix1,:,:) = dens(ix1+1,:,:) + denP1 * (xc(ix1+1)-xc(ix1))
! 	dens(ib1,:,:) = dens(ix1,:,:) + denP1 * (xc(ib1)-xc(ix1))

! 	dens(ix1,:,:) = dens(ix1+1,:,:) 
	dens(ib1,:,:) = dens(ix1,:,:) 

      CASE (2)
!         dens(ix2,:,:) = dens(ix2-1,:,:) + denP1 * (xc(ix2)-xc(ix2-1))
! 	dens(ib2,:,:) = dens(ix2,:,:) + denP1 * (xc(ib2)-xc(ix2))
! 	dens(ix2,:,:) = dens(ix2-1,:,:) 
! 	dens(ib2,:,:) = -dens(ix2,:,:)+2.0d0 
        if(idens.eq.1) then
	do j=1,ny
	dens(ix2,j,:) = dens(ix2-1,j,:)+ denP1*(rp(ix2)-rp(ix2-1))*sin(yc(j))
 	dens(ib2,j,:) = dens(ix2  ,j,:)+ denP1*(rp(ib2)-rp(ix2  ))*sin(yc(j))
! 	dens(ib2,j,:) = dens(ix2,j,:)
	enddo
        else
        do j=1,ny
	dens(ix2,j,:) = dens(ix2-1,j,:)
 	dens(ib2,j,:) = dens(ix2-1,j,:)
	enddo
        endif
      END SELECT
	
      RETURN
      END
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE BOUND500D(dens,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c     boundary conditions for periodic boundaries
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    dens(nx,ny,nz)
c
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (1)
        dens(ib1,:,:) = dens(ix2,:,:)

c      CASE (2)
        dens(ib2,:,:) = dens(ix1,:,:)

      CASE (3)
        dens(:,jb1,:) = dens(:,jy2,:)

c      CASE (4)
        dens(:,jb2,:) = dens(:,jy1,:)

      CASE (5)
        dens(:,:,kb1) = dens(:,:,kz2)

c      CASE (4)
        dens(:,:,kb2) = dens(:,:,kz1)

      END SELECT
c
      RETURN
      END
c
	SUBROUTINE BOUND300D(dens,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c     boundary conditions for axis boundary in cylindrical coordinates
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    dens(nx,ny,nz)
      INTEGER J
c
*
**** Main loop
*
	
      SELECT CASE(IBND)

      CASE (1) 

!  	dens(ix1,:,:) = dens(ix1+1,:,:) 
!   	dens(1 ,:,:) = dens(ix1,:,:) 
          DO J=1,NY
!            DENS(1,J,:) =(DENS(IX1,JSYM(J),:)+DENS(IX1,J,:))*.5
!            DENS(1,J,:) = -DENS(IX1,JSYM(J),:)
! 	   DENS(IX1,JSYM(J),:) =  DENS(IX1+1,JSYM(J),:)
            DENS(1,J,:) =  DENS(IX1,JSYM(J),:)
          ENDDO

      CASE DEFAULT

        WRITE(6,'(A)') '*..Not a defined case!'
        
      END SELECT
c
      RETURN
      END


      SUBROUTINE BOUNDINOUTD(DENS,WO,ZCG,ALFXDT,NX,NY,NZ,NZG,PLANED)
      use density_bg
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'

      INTEGER NX,NY,NZ,I,J,NZG
      REAL DENS(NX,NY,NZ),WO(NX,NY,NZ),ZCG(NZG),ALFXDT,UCONVT,AREA,AREAIJ,COEFF,QIN(2)
      INTEGER IBND
      INTEGER STATUS(MPI_STATUS_SIZE)
      REAL PLANED(NX,NY)    
  
      IF(ITYPE(5)==305) THEN

!      do i=1,nx
!      do j=1,ny

!      dens(i,j,KB1)=planeD(i,j)
      dens(:,:,KB1)=planeD(:,:)

!      enddo
!      enddo

      dens(:,:,2)=dens(:,:,1)

      ENDIF

      IF(ITYPE(6)==730) THEN
      dens(:,:,KZ2) = dens(:,:,KZ2-1)
      dens(:,:,NZ) = dens(:,:,KZ2) 

!       WRITE(*,*)"****************=",KZ2,NZ
!         AREA = 0.0
!         UCONVT = 0.0
! !       IF(ez+1.eq.NZG) THEN
!         DO J=JY1,JY2
!         DO I=IX1,IX2
!           AREAIJ = RP(I)/(AP(I)*BP(J))
!           QIN(1) = QIN(1)+WO(I,J,KB1)*AREAIJ
!           AREA = AREA+AREAIJ
!         ENDDO 
!         ENDDO
!        QIN(2)=QIN(1)/AREA
! !       IF(MYSIZE/=1) 
! !      &       CALL MPI_SEND(QIN,2,MTYPE,MYSIZE-1,0,MPI_COMM_EDDY,IERR)
! ! 
! !       IF(MYSIZE/=1) 
! !      &         CALL MPI_RECV(QIN,2,MTYPE,0,0,MPI_COMM_EDDY,STATUS,IERR)
! 
!        COEFF = QIN(2)*ALFXDT*CP(KZ2)


! extrapolating bc
!         DO J=JY1,JY2
!         DO I=IX1,IX2
!         dens(I,J,KZ2)=dens(I,J,KZ2)-0.50*(WO(I,J,KZ2)+WO(I-1,J,KZ2))*ALFXDT*(dens(I,J,KZ2)-dens(I,J,KZ2-1))/(ZCG(KZ2)-ZCG(KZ2-1))
!         dens(I,J,NZ)=dens(I,J,NZ)-0.50*(WO(I,J,NZ)+WO(I-1,J,NZ))*ALFXDT*(dens(I,J,NZ)-dens(I,J,KZ2))/(ZCG(NZ)-ZCG(KZ2))
!         ENDDO
!         ENDDO
!       ENDIF
      ENDIF
      
      IF(ITYPE(6)==710) THEN
      dens(:,:,KZ2) = dens(:,:,KZ2-1)
      ENDIF
      


      RETURN
      END

      SUBROUTINE BOUNDINOUTD_2(DENS,ALFXDT,NX,NY,NZ)
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'

      INTEGER NX,NY,NZ,I,J
      REAL DENS(NX,NY,NZ),WO(NX,NY,NZ),ALFXDT,UCONVT,AREA,AREAIJ,COEFF,QIN(2)
      INTEGER IBND
      INTEGER STATUS(MPI_STATUS_SIZE)
      
      IF(ITYPE(5)==305) THEN
      dens(:,:,KB1)=dens(:,:,KZ1)
      ENDIF

      IF(ITYPE(6)==730) THEN
      dens(:,:,KZ2) = dens(:,:,KZ2-1)
      dens(:,:,NZ) = dens(:,:,KZ2) 
      ENDIF
      
      IF(ITYPE(6)==710) THEN
      dens(:,:,KZ2) = dens(:,:,KZ2-1)
      ENDIF
      


      RETURN
      END

c----------------------------------------------------------------------
      SUBROUTINE BOUNDARY_P(P,XC,NX,NY,NZ)
c
      INCLUDE 'common.h'
c
      INTEGER NX,NY,NZ
      REAL XC(NX),P(NX,NY,NZ)
c
      INTEGER IBND
c
c.....wall boundary conditions
c
      DO IBND=1,6
	SELECT CASE (ITYPE(IBND))
	CASE(82)
	  CALL BOUND082P(P,NX,NY,NZ,IBND)

	CASE (300)
	  CALL BOUND300P(P,NX,NY,NZ,IBND)

        END SELECT

      ENDDO
c
c.....periodic boundary conditions
c
      DO IBND=1,5,2

        IF(ITYPE(IBND)==500) CALL BOUND500P(P,NX,NY,NZ,IBND)

      ENDDO
c
      RETURN
      END

      SUBROUTINE BOUND082P(P,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c     boundary conditions for non-slip wall
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    p(nx,ny,nz)
c
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (1)
	p(ib1,:,:) = -p(ix1,:,:) 

      CASE (2)
	p(ib2,:,:) = -p(ix2,:,:) 

      END SELECT
c
      RETURN
      END

      SUBROUTINE BOUND500P(p,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c     boundary conditions for periodic boundaries
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    p(nx,ny,nz)
c
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (1)
        p(ib1,:,:) = p(ix2,:,:)

c      CASE (2)
        p(ib2,:,:) = p(ix1,:,:)

      CASE (3)
        p(:,jb1,:) = p(:,jy2,:)

c      CASE (4)
        p(:,jb2,:) = p(:,jy1,:)

      CASE (5)
        p(:,:,kb1) = p(:,:,kz2)

c      CASE (4)
        p(:,:,kb2) = p(:,:,kz1)

      END SELECT
c
      RETURN
      END
c
	SUBROUTINE BOUND300P(P,nx,ny,nz,ibnd)
c----------------------------------------------------------------------
c
c     boundary conditions for axis boundary in cylindrical coordinates
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER ibnd
      INTEGER nx,ny,nz
      REAL    P(nx,ny,nz)
      INTEGER J
c
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (1) 
	DO J=1,NY
	P(1,J,:) =  P(IX1,JSYM(J),:) 
	ENDDO

      CASE DEFAULT

        WRITE(6,'(A)') '*..Not a defined case!'
        
      END SELECT

      RETURN
      END

