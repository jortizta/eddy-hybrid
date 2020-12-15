C---- SUBROUTINE FLAGU -----------------------------------------------
C
C     PURPOSE: Flag points as interior, exterior and forcing
c
      SUBROUTINE FLAGU(FLAG,FLAGO,FLAGI,NX,NY,NZ,NBD,IBD)
C
c---------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'mpif.h'
c
      integer nx,ny,nz,nbd
      integer ibd
      integer flag(nx,ny,nz)
      integer flago(nx,ny,nz,nbd),flagi(nx,ny,nz,nbd)

c....local variables
      integer nxm,nym,nzm,i,j,k,ic,i1,i2,j1,j2
      integer nbor(6)
c....Timing variables
      INTEGER COUNTS,RATE
      REAL    CLOCK
c
c      WRITE(6,*) 'INSIDE FLAGU'
      nxm = nx-1
      nym = ny-1
      nzm = nz-1
c
      flago(:,:,:,ibd) = 0
      flagi(:,:,:,ibd) = 0
c     
c.....define interior
      DO K=KBMIN(IBD),KBMAX(IBD)
      DO J=JBMIN(IBD),JBMAX(IBD)
      DO I=IBMIN(IBD),IBMAX(IBD)
        IF(FLAG(I,J,K)==-1) FLAGO(I,J,K,IBD)=IBD
      ENDDO
      ENDDO
      ENDDO
      FLAGO(:,1,:,IBD) = FLAGO(:,NY-1,:,IBD)
      FLAGO(:,NY,:,IBD) = FLAGO(:,2,:,IBD)
      FLAGI(:,1,:,IBD) = FLAGI(:,NY-1,:,IBD)
      FLAGI(:,NY,:,IBD) = FLAGI(:,2,:,IBD)
      FLAGI = FLAGO

      i1=ibmin(ibd)
      i2=ibmax(ibd)
      j1=jbmin(ibd)
      j2=jbmax(ibd)
      IF(icyl==1) THEN
        j1=2
        j2=ny-1
      ENDIF

c.....define outer and inner boundary
      DO I=MAX(2,i1),i2         !2,nxm
      DO J=j1,j2
      DO K=KBMIN(IBD),KBMAX(IBD) !2,nzm

c        WRITE(6,*) i,j,k
        nbor=1
        nbor(1)=flag(i,j,k)+flag(i+1,j,k)
        nbor(2)=flag(i,j,k)+flag(i-1,j,k)
        nbor(3)=flag(i,j,k)+flag(i,j,k+1)
        nbor(4)=flag(i,j,k)+flag(i,j,k-1)
        nbor(5)=flag(i,j,k)+flag(i,j-1,k)
        nbor(6)=flag(i,j,k)+flag(i,j+1,k)

        if(product(nbor)==0) then
c          WRITE(6,*) 'INSIDE product(nbor)'
c.....outer boundary
          if(flag(i,j,k)==1) then
            flago(i,j,k,ibd)=-ibd
c.....inner boundary             
          elseif(flag(i,j,k)==-1) then
            flagi(i,j,k,ibd)=-ibd
          endif    
        endif
c
      ENDDO
      ENDDO
      ENDDO

      IF(icyl==1 .AND. ibmin(ibd)==1) THEN

        i=1
        DO j=j1,j2              !2,nym
        DO k=KBMIN(IBD),KBMAX(IBD) !2,nzm
c
c            nbor=1
          nbor(1)=flag(i,j,k)+flag(i+1,j,k)
          nbor(2)=flag(i,j,k)+flag(i,j,k+1)
          nbor(3)=flag(i,j,k)+flag(i,j,k-1)
          nbor(4)=flag(i,j,k)+flag(i,j+1,k)
          nbor(5)=flag(i,j,k)+flag(i,j-1,k)
c
          if(product(nbor(1:5))==0) then
c.....outer boundary
            if(flag(i,j,k)==1) then
              flago(i,j,k,ibd)=-ibd
c.....inner boundary             
            elseif(flag(i,j,k)==-1) then
              flagi(i,j,k,ibd)=-ibd
            endif    
          endif
c
        ENDDO
        ENDDO
      ENDIF

c
c...    Periodic boudnary conditions in Y-direction
      FLAGO(:,1 ,:,ibd) = FLAGO(:,JY2,:,ibd)
      FLAGO(:,NY,:,ibd) = FLAGO(:,JY1,:,ibd)
      FLAGI(:,1 ,:,ibd) = FLAGI(:,JY2,:,ibd)
      FLAGI(:,NY,:,ibd) = FLAGI(:,JY1,:,ibd)

      call refreshflag(flago(1,1,1,ibd),nx,ny,nz)
      call refreshflag(flagi(1,1,1,ibd),nx,ny,nz)

c     
      return
      end
c---------------------------------------------------------------------


C---- SUBROUTINE FLAGP -----------------------------------------------
c
c flagp=0: outer pressure points
c flagp=-ibd: interface pressure points
c flagp=ibd: inner pressure points
c The interface points can be inside or outside the immersed-boundary
C
      subroutine flagp(flagpo,flaguo,flagvo,flagwo,nx,ny,nz,nbd,ibd)
C
C     locates i,k points on the cartecian grid near
C     the immersed boundary to be modified  
C     General routine that work for any boundary shape   
C
c      implicit none
      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'mpif.h'
c
      integer nx,ny,nz,nbd,ibd
      integer flaguo(nx,ny,nz,nbd)
      integer flagvo(nx,ny,nz,nbd)
      integer flagwo(nx,ny,nz,nbd)
      integer flagpo(nx,ny,nz,nbd)
c
c....local variables
      integer i,j,k,i1,i2,j1,j2,k1,k2
      integer itst(6),ifld,ibnd,ibdy
      integer nxm,nym,nzm
      integer status(mpi_status_size)
c
      nxm=nx-1
      nym=ny-1
      nzm=nz-1
c
c.....set up the flag for true pressure cells
      
      FLAGPO(:,:,:,ibd) = 0    !!!!!! instead of FLAGPO = 0

      i1 = ibmin(ibd)
      i2 = ibmax(ibd)
      j1 = jbmin(ibd)
      j2 = jbmax(ibd)
      k1 = max(2,kbmin(ibd))
      k2 = kbmax(ibd)
c      k2 = min(nz-1,kbmax(ibd))
      IF(icyl==1) THEN
        j1=2
        j2=ny-1
      ENDIF
c
c.....define the pressure boundary points
      DO K=KBMIN(IBD),1
      DO J=j1,j2
      DO I=max(2,i1),i2
c        
        itst(1)=flaguo(i  ,j,k,ibd)
        itst(2)=flaguo(i-1,j,k,ibd)
        itst(3)=flagwo(i,j,k  ,ibd)
        itst(4)=flagvo(i,j  ,k,ibd)
        itst(5)=flagvo(i,j-1,k,ibd)
c
        ifld=count(itst(1:5)==0)
        ibdy=count(itst(1:5)> 0)
        ibnd=count(itst(1:5)< 0)
c
c.....outer pressure boundary points
        IF((ibdy+ibnd)/=0.AND.ifld/=0) THEN
          flagpo(i,j,k,ibd)=-ibd
        ENDIF
c
        IF(ifld==0.AND.ibnd/=0) THEN
c.....outer pressure body points
          flagpo(i,j,k,ibd)=-ibd
        ENDIF
c
        IF(ifld==0.AND.ibnd==0.AND.ibdy/=0) THEN
c.....outer pressure body points
          flagpo(i,j,k,ibd)=ibd
        ENDIF

c        if(k+myrank*(nz-2)==1 .AND. j==3) then
c          write(6,*) 'flagp:',i,j,k,ibd,ifld,ibdy,ibnd,itst(1:5),flagpo(i,j,k,ibd)
c        endif

      ENDDO
      ENDDO
      ENDDO
c
c.....define the pressure boundary points
      DO K=k1,k2
      DO J=j1,j2
      DO I=max(2,i1),i2
c        
        itst(1)=flaguo(i  ,j,k,ibd)
        itst(2)=flaguo(i-1,j,k,ibd)
        itst(3)=flagwo(i,j,k  ,ibd)
        itst(4)=flagwo(i,j,k-1,ibd)
        itst(5)=flagvo(i,j  ,k,ibd)
        itst(6)=flagvo(i,j-1,k,ibd)
c
        ifld=count(itst==0)
        ibdy=count(itst> 0)
        ibnd=count(itst< 0)

c        if(k+myrank*(nz-2)==2 .AND. j==3) then
c          write(6,*) 'flagp:',i,j,k,ibd,ifld,ibdy,ibnd,itst
c        endif
c
c.....outer pressure boundary points
        IF((ibdy+ibnd)/=0.AND.ifld/=0) THEN
          flagpo(i,j,k,ibd)=-ibd
        ENDIF
c
        IF(ifld==0.AND.ibnd/=0) THEN
c.....outer pressure body points
          flagpo(i,j,k,ibd)=-ibd
        ENDIF
c
        IF(ifld==0.AND.ibnd==0.AND.ibdy/=0) THEN
c.....outer pressure body points
          flagpo(i,j,k,ibd)=ibd
        ENDIF

      ENDDO
      ENDDO
      ENDDO

      IF(icyl==1 .AND. ibmin(ibd)==1) THEN
c.....define the pressure boundary points
        i = 1

        DO k=KBMIN(IBD),1
        DO j=j1,j2
c        
          itst(1)=flaguo(i  ,j,k,ibd)
          itst(2)=flagwo(i,j,k  ,ibd)
          itst(3)=flagvo(i,j  ,k,ibd)
          itst(4)=flagvo(i,j-1,k,ibd)

          ifld=count(itst(1:4)==0)
          ibdy=count(itst(1:4)> 0)
          ibnd=count(itst(1:4)< 0)
c
c.....outer pressure boundary points
          IF(ibdy==0.AND.ibnd/=0.AND.ifld/=0) THEN
            flagpo(i,j,k,ibd)=-ibd
          ENDIF
c
c          IF(ifld==0.AND.ibnd/=0.AND.ibdy/=0) THEN
          IF(ifld==0.AND.ibnd/=0) THEN
c.....outer pressure body points
            flagpo(i,j,k,ibd)=ibd
          ENDIF
c
          IF(ifld==0.AND.ibnd==0.AND.ibdy/=0) THEN
c.....outer pressure body points
            flagpo(i,j,k,ibd)=ibd
          ENDIF

        ENDDO
        ENDDO

        DO k=k1,k2
        DO j=j1,j2
c        
          itst(1)=flaguo(i  ,j,k,ibd)
          itst(2)=flagwo(i,j,k  ,ibd)
          itst(3)=flagwo(i,j,k-1,ibd)
          itst(4)=flagvo(i,j  ,k,ibd)
          itst(5)=flagvo(i,j-1,k,ibd)

          ifld=count(itst(1:5)==0)
          ibdy=count(itst(1:5)> 0)
          ibnd=count(itst(1:5)< 0)
c
c.....outer pressure boundary points
          IF(ibdy==0.AND.ibnd/=0.AND.ifld/=0) THEN
            flagpo(i,j,k,ibd)=-ibd
          ENDIF
c
c          IF(ifld==0.AND.ibnd/=0.AND.ibdy/=0) THEN
          IF(ifld==0.AND.ibnd/=0) THEN
c.....outer pressure body points
            flagpo(i,j,k,ibd)=ibd
          ENDIF
c
          IF(ifld==0.AND.ibnd==0.AND.ibdy/=0) THEN
c.....outer pressure body points
            flagpo(i,j,k,ibd)=ibd
          ENDIF

        ENDDO
        ENDDO
      ENDIF
c
c.....Periodic boundary conditions in Y-direction
      if(idomy==0) then
        FLAGPO(:,1,:,IBD) = FLAGPO(:,NY-1,:,IBD)
        FLAGPO(:,NY,:,IBD) = FLAGPO(:,2,:,IBD)
c
c.....Axis
        IF(ITYPE(1)==300) THEN
          DO J=1,NY
            FLAGPO(1,J,:,IBD) = FLAGPO(IX1,JSYM(J),:,IBD)
          ENDDO
        ENDIF
      else
        if(itype(1)==300) then
          call mpi_sendrecv(flagpo(2,:,:,ibd),ny*nz,mpi_integer,ranksym,1
     &          ,           flagpo(1,:,:,ibd),ny*nz,mpi_integer,ranksym,1
     &          ,           mpi_comm_eddy,status,ierr)
        endif         
      endif

      call refreshflag(flagpo(1,1,1,ibd),nx,ny,nz)
     
      RETURN

      END
c
c---------------------------------------------------------------------



C---- SUBROUTINE TAGPBD --------------------N. Beratlis-16 Apr. 2009-
C
C     PURPOSE: Flag points as interior, exterior and forcing
c
      SUBROUTINE TAGPBD(FLAGP,FLAGPBD,XC,YC,ZC,NX,NY,NZ,IBD
     &                 ,VERTEXC,NFACET)
C
c---------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      include 'mpif.h'
c
      integer nx,ny,nz,ibd,nfacet
      integer flagp(nx,ny,nz),flagpbd(nx,ny,nz)
      real    vertexc(3,nfacet)
      real    xc(nx),yc(ny),zc(nz)
c
c....local variables
      integer i,j,k,ii,j1,j2,k1
      real    xp,yp,zp
      logical cond

      flagpbd = 0

      WHERE(flagp(:,:,:)==ibd) flagpbd=ibd
      
      DO ii=1,nfacet

        xp = vertexc(1,ii)
        yp = vertexc(2,ii)
        zp = vertexc(3,ii)

        if(icyl==0) then
          cond = zp>=zc(1).AND.zp<zc(nz-1).AND.yp>=yc(1).AND.yp<yc(ny)
        else
          cond = zp>=zc(1).AND.zp<zc(nz-1)
        endif

        IF(cond) THEN
          call ijk_xyz(xp,yp,zp,xc,yc,zc,nx,ny,nz,i,j,k,icyl)

          flagpbd(i  ,j  ,k  )=-ibd
          flagpbd(i  ,j+1,k  )=-ibd
          flagpbd(i  ,j  ,k+1)=-ibd
          flagpbd(i  ,j+1,k+1)=-ibd
          flagpbd(i+1,j  ,k  )=-ibd
          flagpbd(i+1,j+1,k  )=-ibd
          flagpbd(i+1,j  ,k+1)=-ibd
          flagpbd(i+1,j+1,k+1)=-ibd

        ENDIF
      ENDDO

      WHERE(flagp(:,:,:)<=0) flagpbd=0      
c
c.... Periodic boundary conditions in Y-direction
      FLAGPBD(:,1 ,:) = FLAGPBD(:,JY2,:)
      FLAGPBD(:,NY,:) = FLAGPBD(:,JY1,:)
c
c.... Axis
      IF(ICYL==1 .AND. ITYPE(1)==300) THEN
        DO J=1,NY
          FLAGPBD(1,J,:) = FLAGPBD(2,JSYM(J),:)
        ENDDO
      ENDIF
      
c      write(6,*) 'tagpbd, n=',count(flagpbd(2:nx,2:ny-1,:)==-ibd)
c     
      return
      end
c---------------------------------------------------------------------

C---- SUBROUTINE TAGPBD1 -------------------N. Beratlis-16 Apr. 2009-
C
C     PURPOSE: Flag points as interior, exterior and forcing
c
      SUBROUTINE TAGPBD1(FLAGP,FLAGPBD,XC,YC,XC_CAR,YC_CAR,ZC,NX,NY,NZ,IBD
     &                 ,VERTEXC,UNVECT,NFACET)
C
c---------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      include 'mpif.h'
c
      integer nx,ny,nz,ibd,nfacet
      integer flagp(nx,ny,nz),flagpbd(nx,ny,nz)
      real    vertexc(3,nfacet),unvect(3,nfacet)
      real    xc(nx),yc(ny),zc(nz)
      real    xc_car(nx,ny),yc_car(nx,ny)
c
c....local variables
      integer i,j,k,ii,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4,dir,ksign
      real    xp,yp,zp,xint,yint,zint
      logical cond
c
c.... Functions
      integer body2grid

      flagpbd = 0

      WHERE(flagp(:,:,:)==ibd) flagpbd=ibd    ! interior points
      
      DO ii=1,nfacet

        xp = vertexc(1,ii)
        yp = vertexc(2,ii)
        zp = vertexc(3,ii)

        ksign = 0
        if(unvect(3,ii)>0) ksign=1

        if(icyl==0) then
c          cond = zp>=zc(1).AND.zp<zc(nz-1).AND.yp>=yc(1).AND.yp<yc(ny)
          cond = zp>=zc(2-ksign).AND.zp<zc(nz-ksign).AND.yp>=yc(1).AND.yp<yc(ny)
        else
c          cond = zp>=zc(1).AND.zp<zc(nz-1)
          cond = zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
        endif

        if(cond) then
c from the immersed-boundary an intersection with a grid face is found
c along the outward normal direction
          if(body2grid(vertexc(:,ii),unvect(:,ii),xc_car,yc_car,xc,yc,zc
     &          ,nx,ny,nz,xint,yint,zint,i,j,k,dir,icyl)==1) then

c the grid points corresponding to the intersection are tagged as -IBD
            if(dir==1) then
              call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
              call index_bnd(i  ,j+1,k  ,i2,j2,k2,nx,ny,nz)
              call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
              call index_bnd(i  ,j+1,k+1,i4,j4,k4,nx,ny,nz)
              flagpbd(i1,j1,k1)=-ibd
              flagpbd(i2,j2,k2)=-ibd
              flagpbd(i3,j3,k3)=-ibd
              flagpbd(i4,j4,k4)=-ibd              
            elseif(dir==2) then
              call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
              call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
              call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
              call index_bnd(i+1,j  ,k+1,i4,j4,k4,nx,ny,nz)
              flagpbd(i1,j1,k1)=-ibd
              flagpbd(i2,j2,k2)=-ibd
              flagpbd(i3,j3,k3)=-ibd
              flagpbd(i4,j4,k4)=-ibd
            elseif(dir==3) then
              call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
              call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
              call index_bnd(i  ,j+1,k  ,i3,j3,k3,nx,ny,nz)
              call index_bnd(i+1,j+1,k  ,i4,j4,k4,nx,ny,nz)
              flagpbd(i1,j1,k1)=-ibd
              flagpbd(i2,j2,k2)=-ibd
              flagpbd(i3,j3,k3)=-ibd
              flagpbd(i4,j4,k4)=-ibd              
            endif

c            if(ii==17971) then
c              write(6,*) '1.tagpbd1',ii,xint,yint,zint,i,j,k,dir
c     &              ,flagpbd(i1,j1,k1),flagpbd(i2,j2,k2),flagpbd(i3,j3,k3),flagpbd(i4,j4,k4)
c            endif

c            write(6,*) ii,flagpbd(100,5,276),flagpbd(101,5,276),flagpbd(100,6,276),flagpbd(101,6,276)

          endif
        endif
      enddo

c      i=100
c      j=5
c      k=276
c      write(6,*) '1. tagbd1',flagpbd(i,j,k),flagpbd(i+1,j,k),flagpbd(i,j+1,k),flagpbd(i+1,j+1,k)
c      write(6,*) '2. tagbd1',flagp(i,j,k),flagp(i+1,j,k),flagp(i,j+1,k),flagp(i+1,j+1,k)

c FLAGP=IBD for the interior points
c FLAGP=-IBD for the interface points
c FLAGP=0 for the fluid points
      where(flagp(:,:,:)<=0) flagpbd=0
c FLAGPBD is set 0 at the interface and fluid points
c finally FLAGPBD is IBD at the body points, 0 at the interface and fluid
c points and -IBD at the points found by BODY2GRID which are body points      

c      write(6,*) '3. tagbd1',flagpbd(i,j,k),flagpbd(i+1,j,k),flagpbd(i,j+1,k),flagpbd(i+1,j+1,k)
c
c.... Periodic boundary conditions in Y-direction
      FLAGPBD(:,1 ,:) = FLAGPBD(:,JY2,:)
      FLAGPBD(:,NY,:) = FLAGPBD(:,JY1,:)

c      write(6,*) '4. tagbd1',flagpbd(i,j,k),flagpbd(i+1,j,k),flagpbd(i,j+1,k),flagpbd(i+1,j+1,k)

c
c.... Axis
      IF(ICYL==1 .AND. ITYPE(1)==300) THEN
        DO J=1,NY
          FLAGPBD(1,J,:) = FLAGPBD(2,JSYM(J),:)
        ENDDO
      ENDIF      
c     
      call refreshflag(flagpbd(1,1,1),nx,ny,nz)

      return
      end
c---------------------------------------------------------------------



C---- SUBROUTINE TAGPBD2 -------------------N. Beratlis-16 Apr. 2009-
C
C     PURPOSE: Flag points as interior, exterior and forcing
c
      SUBROUTINE TAGPBD2(FLAGP,FLAGPBD,XC,YC,XC_CAR,YC_CAR,ZC,ZCG
     &     ,NX,NY,NZ,NZG,IBD,VERTEXC,UNVECT,NFACET)
C
c---------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      include 'mpif.h'
c
      integer nx,ny,nz,nzg,ibd,nfacet
      integer flagp(nx,ny,nz),flagpbd(nx,ny,nz)
      real    vertexc(3,nfacet),unvect(3,nfacet)
      real    xc(nx),yc(ny),zc(nz),zcg(nzg)
      real    xc_car(nx,ny),yc_car(nx,ny)
c
c....local variables
      integer i,j,k,ii,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4,dir,ksign,kg
      real    xp,yp,zp,xint,yint,zint
      logical cond
c
c.... Functions
      integer body2grid

      flagpbd = 0

      WHERE(flagp(:,:,:)==ibd) flagpbd=ibd   ! body points
      
      DO ii=1,nfacet

        xp = vertexc(1,ii)
        yp = vertexc(2,ii)
        zp = vertexc(3,ii)

        if(icyl==0) then
          cond = zp>=zc(1).AND.zp<=zc(nz).AND.yp>=yc(1).AND.yp<yc(ny)
        else
          cond = zp>=zc(1).AND.zp<=zc(nz)
        endif

        if(cond) then
c from the immersed-boundary an intersection with a grid face is found
c along the outward normal direction
          if(body2grid(vertexc(:,ii),unvect(:,ii),xc_car,yc_car,xc,yc,zcg
     &          ,nx,ny,nzg,xint,yint,zint,i,j,kg,dir,icyl)==1) then
             
            k=kg-myrank*(nz-2)

            if( (dir/=3 .AND. k>=1 .AND. k<nz) .OR. (dir==3 .AND. k>=1 .AND. k<=nz)) then

c the grid points corresponding to the intersection are tagged as -IBD
              if(dir==1) then
                call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
                call index_bnd(i  ,j+1,k  ,i2,j2,k2,nx,ny,nz)
                call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
                call index_bnd(i  ,j+1,k+1,i4,j4,k4,nx,ny,nz)
                flagpbd(i1,j1,k1)=-ibd
                flagpbd(i2,j2,k2)=-ibd
                flagpbd(i3,j3,k3)=-ibd
                flagpbd(i4,j4,k4)=-ibd              
              elseif(dir==2) then
                call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
                call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
                call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
                call index_bnd(i+1,j  ,k+1,i4,j4,k4,nx,ny,nz)
                flagpbd(i1,j1,k1)=-ibd
                flagpbd(i2,j2,k2)=-ibd
                flagpbd(i3,j3,k3)=-ibd
                flagpbd(i4,j4,k4)=-ibd
              elseif(dir==3) then
                call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
                call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
                call index_bnd(i  ,j+1,k  ,i3,j3,k3,nx,ny,nz)
                call index_bnd(i+1,j+1,k  ,i4,j4,k4,nx,ny,nz)
                flagpbd(i1,j1,k1)=-ibd
                flagpbd(i2,j2,k2)=-ibd
                flagpbd(i3,j3,k3)=-ibd
                flagpbd(i4,j4,k4)=-ibd              
              endif
            endif

          endif
        endif
      enddo

c at the interface and fluid points FLAGPBD is set equal to 0
      where(flagp(:,:,:)<=0) flagpbd=0      
c FLAGPBD=0: interface and fluid points
c FLAGPBD=IBD: interior points
c FLAGPBD=-IBD: interior points of the interpolation stencils

c
c.... Periodic boundary conditions in Y-direction
      FLAGPBD(:,1 ,:) = FLAGPBD(:,JY2,:)
      FLAGPBD(:,NY,:) = FLAGPBD(:,JY1,:)

c
c.... Axis
      IF(ICYL==1 .AND. ITYPE(1)==300) THEN
        DO J=1,NY
          FLAGPBD(1,J,:) = FLAGPBD(2,JSYM(J),:)
        ENDDO
      ENDIF      
c     
      call refreshflag(flagpbd(1,1,1),nx,ny,nz)

      return
      end
c---------------------------------------------------------------------



C---- SUBROUTINE TAGPBD3 -------------------N. Beratlis-31 Oct. 2010-
C
C     PURPOSE: Flag points as interior, exterior and forcing
c
      SUBROUTINE TAGPBD3(FLAGP,FLAGPBD,XC,YC,XC_CAR,YC_CAR,ZC,ZCG
     &     ,NX,NY,NZ,NZG,IBD,VERTEX,UNVECT,NFACET)
C
c---------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      include 'mpif.h'
c
      integer nx,ny,nz,nzg,ibd,nfacet
      integer flagp(nx,ny,nz),flagpbd(nx,ny,nz)
      real    vertex(3,3,nfacet),unvect(3,nfacet)
      real    xc(nx),yc(ny),zc(nz),zcg(nzg)
      real    xc_car(nx,ny),yc_car(nx,ny)
c
c....local variables
      integer i,j,k,ii,jj,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4,dir,ksign,kg
      real    xp,yp,zp,xint,yint,zint
      logical cond
c
c.... Functions
      integer body2grid

      flagpbd = 0

      WHERE(flagp(:,:,:)==ibd) flagpbd=ibd   ! body points
      
      DO ii=1,nfacet
      DO jj=1,3
         
        xp = vertex(1,jj,ii)
        yp = vertex(2,jj,ii)
        zp = vertex(3,jj,ii)

        if(icyl==0) then
          cond = zp>=zc(1).AND.zp<=zc(nz).AND.yp>=yc(1).AND.yp<yc(ny)
        else
          cond = zp>=zc(1).AND.zp<=zc(nz)
        endif

        if(cond) then
c from the immersed-boundary an intersection with a grid face is found
c along the outward normal direction
          if(body2grid(vertex(:,jj,ii),unvect(:,ii),xc_car,yc_car,xc,yc,zcg
     &          ,nx,ny,nzg,xint,yint,zint,i,j,kg,dir,icyl)==1) then
             
            k=kg-myrank*(nz-2)

            if( (dir/=3 .AND. k>=1 .AND. k<nz) .OR. (dir==3 .AND. k>=1 .AND. k<=nz)) then

c the grid points corresponding to the intersection are tagged as -IBD
              if(dir==1) then
                call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
                call index_bnd(i  ,j+1,k  ,i2,j2,k2,nx,ny,nz)
                call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
                call index_bnd(i  ,j+1,k+1,i4,j4,k4,nx,ny,nz)
                flagpbd(i1,j1,k1)=-ibd
                flagpbd(i2,j2,k2)=-ibd
                flagpbd(i3,j3,k3)=-ibd
                flagpbd(i4,j4,k4)=-ibd              
              elseif(dir==2) then
                call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
                call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
                call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
                call index_bnd(i+1,j  ,k+1,i4,j4,k4,nx,ny,nz)
                flagpbd(i1,j1,k1)=-ibd
                flagpbd(i2,j2,k2)=-ibd
                flagpbd(i3,j3,k3)=-ibd
                flagpbd(i4,j4,k4)=-ibd
              elseif(dir==3) then
                call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
                call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
                call index_bnd(i  ,j+1,k  ,i3,j3,k3,nx,ny,nz)
                call index_bnd(i+1,j+1,k  ,i4,j4,k4,nx,ny,nz)
                flagpbd(i1,j1,k1)=-ibd
                flagpbd(i2,j2,k2)=-ibd
                flagpbd(i3,j3,k3)=-ibd
                flagpbd(i4,j4,k4)=-ibd              
              endif
            endif

          endif
        endif

      ENDDO
      ENDDO

c at the interface and fluid points FLAGPBD is set equal to 0
      where(flagp(:,:,:)<=0) flagpbd=0      
c FLAGPBD=0: interface and fluid points
c FLAGPBD=IBD: interior points
c FLAGPBD=-IBD: interior points of the interpolation stencils

c
c.... Periodic boundary conditions in Y-direction
      FLAGPBD(:,1 ,:) = FLAGPBD(:,JY2,:)
      FLAGPBD(:,NY,:) = FLAGPBD(:,JY1,:)

c
c.... Axis
      IF(ICYL==1 .AND. ITYPE(1)==300) THEN
        DO J=1,NY
          FLAGPBD(1,J,:) = FLAGPBD(2,JSYM(J),:)
        ENDDO
      ENDIF      
c     
      call refreshflag(flagpbd(1,1,1),nx,ny,nz)

      return
      end
c---------------------------------------------------------------------


C---- subroutine grid_intr------------------------N. Beratlis-04 Jan. 2009--
C
C     PURPOSE: Given a point im,jm,km find intesection with immersed body
C     along grid points (gridlines plus diagonals).
c     RETURNS: 1 if triangle intrs is found and stencil is physical - SUCCESS
C              0 if no triangle intrs is found
C             -1 if triangle intrs is found, but stencil is forcing point
C             -2 if triangle intrs is found, but stencil is body point
C
C---------------------------------------------------------------------------
      integer function grid_intr(flag,nx,ny,nz,ibd,vertex,unvect,nfacet
     &        ,xu,yu,zu,inn,nn,xnim,ynim,znim,nxim,nyim,nzim,im,jm,km,isrt,fp,io,itri)
C
c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'

      INTEGER nx,ny,nz,nfacet,im,jm,km,nn,ibd,isrt,fp,io,itri
      INTEGER inn(nn)
      INTEGER flag(nx,ny,nz)
      REAL    xu(nx,ny),yu(nx,ny),zu(nz)
      REAL    vertex(3,3,nfacet),unvect(3,nfacet)
      REAL    xnim,ynim,znim,nxim,nyim,nzim

      INTEGER iflagt,intrs,ibint,in,imin,imax,i,ii,isten,n,nitrs,ext,nphy
      REAL    xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,a
      REAL    sdiag,nmag,nnx,nny,nnz
      REAL    stencil_normal_unitdot
      REAL    rvec(3),rvec2(3)
      LOGICAL cond
      REAL    XINT(26),YINT(26),ZINT(26),DS(26)
      INTEGER NXINT(26),NYINT(26),NZINT(26),IORDER(26),ITRIANGLE(26)
      INTEGER ind(26,3),icand(26),ord(26)
      REAL    crit(26)
      INTEGER id,jd,kd,nind,ncand
c
c.... Functions
      REAL    vecmag,dotprod

      INTRS = 0
      ISTEN = 0
      IBINT = 0
      ncand = 0

      ext = 10.0
      itri = 0

      nind = 26
      ind = 0

      ind(1,1) = 1
      ind(1,2) = 0
      ind(1,3) = 0

      ind(2,1) =-1
      ind(2,2) = 0
      ind(2,3) = 0

      ind(3,1) = 0
      ind(3,2) = 0
      ind(3,3) = 1

      ind(4,1) = 0
      ind(4,2) = 0
      ind(4,3) =-1

      ind(5,1) = 1
      ind(5,2) = 0
      ind(5,3) = 1

      ind(6,1) = 1
      ind(6,2) = 0
      ind(6,3) =-1

      ind(7,1) =-1
      ind(7,2) = 0
      ind(7,3) = 1

      ind(8,1) =-1
      ind(8,2) = 0
      ind(8,3) =-1

      ind(9,1) = 0
      ind(9,2) = 1
      ind(9,3) = 0

      ind(10,1) = 0
      ind(10,2) =-1
      ind(10,3) = 0

      ind(11,1) = 1
      ind(11,2) = 1
      ind(11,3) = 0

      ind(12,1) =-1
      ind(12,2) = 1
      ind(12,3) = 0

      ind(13,1) = 0
      ind(13,2) = 1
      ind(13,3) = 1

      ind(14,1) = 0
      ind(14,2) = 1
      ind(14,3) =-1

      ind(15,1) = 1
      ind(15,2) = 1
      ind(15,3) = 1

      ind(16,1) = 1
      ind(16,2) = 1
      ind(16,3) =-1

      ind(17,1) =-1
      ind(17,2) = 1
      ind(17,3) = 1

      ind(18,1) =-1
      ind(18,2) = 1
      ind(18,3) =-1

      ind(19,1) = 1
      ind(19,2) =-1
      ind(19,3) = 0

      ind(20,1) =-1
      ind(20,2) =-1
      ind(20,3) = 0

      ind(21,1) = 0
      ind(21,2) =-1
      ind(21,3) = 1

      ind(22,1) = 0
      ind(22,2) =-1
      ind(22,3) =-1

      ind(23,1) = 1
      ind(23,2) =-1
      ind(23,3) = 1

      ind(24,1) = 1
      ind(24,2) =-1
      ind(24,3) =-1

      ind(25,1) =-1
      ind(25,2) =-1
      ind(25,3) = 1

      ind(26,1) =-1
      ind(26,2) =-1
      ind(26,3) =-1

      !Determine order of grid ray tracing
      if(isrt==0) then
        crit = -2.0
        do i=1,nind
          id = ind(i,1)
          jd = ind(i,2)
          kd = ind(i,3)
          !Give priority to stencil with fluid points
          if(flag(im+id,jm+jd,km+kd)==0) then
c fluid points
             crit(i) = 1.0
          elseif(flag(im+id,jm+jd,km+kd)<0) then
c boundary points
            crit(i) = 0.5
          endif
        enddo
      elseif(isrt==1) then
        crit = -2.0
        do i=1,nind
          id = ind(i,1)
          jd = ind(i,2)
          kd = ind(i,3)
          rvec(1) = xu(im+id,jm+jd)-xu(im,jm)
          rvec(2) = yu(im+id,jm+jd)-yu(im,jm)
          rvec(3) = zu(km+kd)-zu(km)
          rvec = rvec/vecmag(rvec,3)
          !Give priority to stencil with fluid points
          if(flag(im+id,jm+jd,km+kd)==0) then
            crit(i) = dotprod(rvec,unvect(:,fp))    !!!!!!
          endif
        enddo
      else
        crit = -2.0
        do i=1,nind
          id = ind(i,1)
          jd = ind(i,2)
          kd = ind(i,3)
          rvec(1) = xu(im+id,jm+jd)-xu(im,jm)
          rvec(2) = yu(im+id,jm+jd)-yu(im,jm)
          rvec(3) = zu(km+kd)-zu(km)
          rvec = rvec/vecmag(rvec,3)
          !Give priority to stencil with fluid points
          if(flag(im+id,jm+jd,km+kd)<=0) then
            crit(i) = dotprod(rvec,unvect(:,fp))    !!!!!!
          endif
        enddo
      endif

      nphy = 0
      do i=1,nind
        ord(i) = maxloc(crit(1:nind),1)
        crit(ord(i)) = -10.0
c only fluid and boundary points are taken into account
        if(flag(im+ind(i,1),jm+ind(i,2),km+ind(i,3))<=0) nphy = nphy+1
      enddo

      do i=1,nphy

        ii = ord(i)
        id = ind(ii,1)
        jd = ind(ii,2)
        kd = ind(ii,3)
c        if(io==1) then
c          cond = flag(im-id,jm-jd,km-kd)==ibd !Opposite is body TRUE WHEN Dx~Dy~Dz, not good when Dy is large
          cond = flag(im+id,jm+jd,km+kd)<=0 !Extension point is physical MORE ROBUST
c        else
c          cond = flag(im+id,jm+jd,km+kd)<=0
c        endif

        if(cond) then
c          isten = isten+1
          do in=1,nn
            n = inn(in)
            xa = VERTEX(1,1,n)
            ya = VERTEX(2,1,n)
            za = VERTEX(3,1,n)
            xb = VERTEX(1,2,n)
            yb = VERTEX(2,2,n)
            zb = VERTEX(3,2,n)
            xc = VERTEX(1,3,n)
            yc = VERTEX(2,3,n)
            zc = VERTEX(3,3,n)

            rvec(1) = xu(im+id,jm+jd)-xu(im,jm)
            rvec(2) = yu(im+id,jm+jd)-yu(im,jm)
            rvec(3) = zu(km+kd)-zu(km)
            rvec = rvec*ext
            call segtriint(xu(im,jm)-rvec(1),yu(im,jm)-rvec(2),zu(km)-rvec(3)
     &           ,xu(im  ,jm)+rvec(1),yu(im  ,jm)+rvec(2),zu(km)+rvec(3)
     &           ,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)

c            call segtriint(xu(im-id,jm-jd),yu(im-id,jm-jd),zu(km-kd)
c     &           ,xu(im  ,jm),yu(im  ,jm),zu(km)
c     &           ,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)
c            if(id==0 .AND. jd==1 .AND. kd==0 .AND.

            IF(IFLAGT==1) THEN

              if(ibint==0) itri = in
              ibint = ibint+1
                            
c              rvec(1) = xu(im+id,jm+jd)-xu(im,jm)
c              rvec(2) = yu(im+id,jm+jd)-yu(im,jm)
c              rvec(3) = zu(km+kd)-zu(km)
c              rvec = rvec/vecmag(rvec,3)
c              rvec2(1) = xmt-xu(im,jm)
c              rvec2(2) = ymt-yu(im,jm)
c              rvec2(3) = zmt-zu(km)
c              rvec2 = rvec2/vecmag(rvec2,3)
c              a = dotprod(rvec,rvec2,3)
c
c evaluation of a dot product between unit vectors to determine if the
c extension point is along the outward direction
c 
              a = stencil_normal_unitdot(xu(im+id,jm+jd),yu(im+id,jm+jd)
     $             ,zu(km+kd),xu(im,jm),yu(im,jm),zu(km),unvect(1,n)
     $             ,unvect(2,n),unvect(3,n))
c              if(io==1 .OR. (io==0 .AND. a<0.0)) then
              if(a>0.0) then
c the extension point is on the side of the outward direction
                if(flag(im+id,jm+jd,km+kd)==0) then
c the extension point is fluid
                  isten = isten+1
                  INTRS = INTRS + 1
                  ITRIANGLE(INTRS) = N
                  XINT (INTRS) = XMT
                  YINT (INTRS) = YMT
                  ZINT (INTRS) = ZMT
                  NXINT(INTRS) = id
                  NYINT(INTRS) = jd
                  NZINT(INTRS) = kd
                  icand(intrs) = 0
                  GO TO 300
                elseif(flag(im+id,jm+jd,km+kd)<0) then
c the extension point is boundary
                  ncand = ncand+1
                  isten = isten+1
                  INTRS = INTRS + 1
                  ITRIANGLE(INTRS) = N
                  XINT (INTRS) = XMT
                  YINT (INTRS) = YMT
                  ZINT (INTRS) = ZMT
                  NXINT(INTRS) = id
                  NYINT(INTRS) = jd
                  NZINT(INTRS) = kd
                  icand(intrs) = 1
                  GO TO 300
c                grid_intr = -2
                else
                  grid_intr =-2
                endif
              ENDIF
            ENDIF  
          enddo

 300      continue
          
        endif
      enddo

c 300  continue

      IF(INTRS==0) THEN
        if(ibint==0) then
          grid_intr=0
        else
          grid_intr=-2
        endif
      ELSE
        
        grid_intr = 1
        IF(INTRS>nind) THEN
          WRITE(6,*) 'PROBLEM: MORE THAN 26 INTERSECTIONS'
          INTRS=nind
        ENDIF
        DO nitrs=1,intrs
          ii = itriangle(nitrs)
          id = int(nxint(nitrs))
          jd = int(nyint(nitrs))
          kd = int(nzint(nitrs))
c          a = stencil_normal_unitdot(xu(im,jm),yu(im,jm),zu(km)
c     &          ,xint(nitrs),yint(nitrs),zint(nitrs)
c     &          ,unvect(1,ii),unvect(2,ii),unvect(3,ii))
c          ds(nitrs)= abs(stencil_normal_unitdot(xu(im,jm),yu(im,jm),zu(km)
c     &          ,xint(nitrs),yint(nitrs),zint(nitrs)
c     &          ,unvect(1,ii),unvect(2,ii),unvect(3,ii)))
c
c dot product between the unit extension vector and the unit vector normal to the boundary
c
          ds(nitrs)= stencil_normal_unitdot(xu(im+id,jm+jd),yu(im+id,jm+jd),zu(km+kd)
     &          ,xu(im,jm),yu(im,jm),zu(km)
     &          ,unvect(1,ii),unvect(2,ii),unvect(3,ii))
          if(igrdint==1) then
c fluid points before boundary points
            ds(nitrs) = ds(nitrs) + flag(im+id,jm+jd,km+kd)*nitrs
          endif

c          DS(NITRS) = SQRT( (XU(IM,JM)-XINT(NITRS))**2. 
c     &          +  (YU(IM,JM)-YINT(NITRS))**2.
c     &          +  (ZU(KM)-ZINT(NITRS))**2. )
        ENDDO
        ii = maxloc(ds(1:intrs),1)

        XNIM = XINT(II)
        YNIM = YINT(II)
        ZNIM = ZINT(II)
        NXIM = REAL(NXINT(II))
        NYIM = REAL(NYINT(II))
        NZIM = REAL(NZINT(II))
        FP = ITRIANGLE(II)

        if(icand(ii)==1) grid_intr=-1

      ENDIF


      return
      end
C---------------------------------------------------------------------------



C---- subroutine grid_intr_p----------------------N. Beratlis-08 Jan. 2009--
C
C     PURPOSE: Given a point im,jm,km find 2 physical (fluid) points
C     along grid points (gridlines plus diagonals). It returns:
c     RETURNS: 1 if triangle intrs is found and stencil is physical - SUCCESS
C              0 if triangle intrs is found, but stencil is body point
C             -1 if no triangle intrs is found
C             -2 if triangle intrs is found, but stencil is forcing point
C
C---------------------------------------------------------------------------
C
      integer function grid_intr_p(flag,flag2,nx,ny,nz,ibd,nbd
     &        ,vertex,unvect,nfacet,xu,yu,zu,inn,nn
     &        ,xnim,ynim,znim,nxim,nyim,nzim,im,jm,km,itr)
C
c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'

      INTEGER nx,ny,nz,nbd,nfacet,im,jm,km,nn,ibd,itr
      INTEGER inn(nn)
      INTEGER flag(nx,ny,nz),flag2(nx,ny,2)
      REAL    xu(nx,ny),yu(nx,ny),zu(nz)
      REAL    vertex(3,3,nfacet),unvect(3,nfacet)
      REAL    xnim,ynim,znim,nxim,nyim,nzim

      INTEGER iflagt,intrs,in,imin,imax,i,ii,isten,n,nitrs,ibint
      REAL    xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,ext
      REAL    rvec(3)
      REAL    stencil_normal_unitdot
      INTEGER ind(6,3)
      INTEGER im1,jm1,km1,im2,jm2,km2,id,jd,kd,nd,ncand
      LOGICAL cond1,cond2
      REAL    XINT(6),YINT(6),ZINT(6),DS(6)
      INTEGER NXINT(6),NYINT(6),NZINT(6),IORDER(6),ITRIANGLE(6),ICAND(6)
      INTEGER crit(6),ord(6)

      INTRS = 0
      ISTEN = 0
      ibint = 0
      ncand = 0

      ext = 10.0

      nd=4

      ind(1,1) = 1
      ind(1,2) = 0
      ind(1,3) = 0
      
      ind(2,1) =-1
      ind(2,2) = 0
      ind(2,3) = 0

      ind(3,1) = 0
      ind(3,2) = 0
      ind(3,3) = 1

      ind(4,1) = 0
      ind(4,2) = 0
      ind(4,3) =-1

      ind(5,1) = 0
      ind(5,2) = 1
      ind(5,3) = 0

      ind(6,1) = 0
      ind(6,2) =-1
      ind(6,3) = 0

      !Determine order of grid ray tracing
      crit = 0
      do i=1,nd
        id = ind(i,1)
        jd = ind(i,2)
        kd = ind(i,3)
        !Give priority to stencil with fluid points
        if(flag(im+id,jm+jd,km+kd)==0) then
          crit(i) = 1
        endif
      enddo

      do i=1,nd
        ord(i) = maxloc(crit(1:nd),1)
        crit(ord(i)) = -1
      enddo

      do i=1,nd

        ii=ord(i) 

        id = ind(ii,1)
        jd = ind(ii,2)
        kd = ind(ii,3)

        im1 = im + id
        jm1 = jm + jd
        km1 = km + kd

        im2 = im1 + id
        jm2 = jm1 + jd
        km2 = km1 + kd

        call index_bnd(im2,jm2,km2,im2,jm2,km2,nx,ny,nz)

        do in=1,nn
          n=inn(in)
          xa = VERTEX(1,1,n)
          ya = VERTEX(2,1,n)
          za = VERTEX(3,1,n)
          xb = VERTEX(1,2,n)
          yb = VERTEX(2,2,n)
          zb = VERTEX(3,2,n)
          xc = VERTEX(1,3,n)
          yc = VERTEX(2,3,n)
          zc = VERTEX(3,3,n)

          call unitvector(xu(im,jm),xu(im1,jm1),yu(im,jm),yu(im1,jm1),zu(km),zu(km1),rvec)
          call segtriint(xu(im,jm)-ext*rvec(1),yu(im,jm)-ext*rvec(2),zu(km)-ext*rvec(3)
     &         ,xu(im,jm)+ext*rvec(1),yu(im,jm)+ext*rvec(2),zu(km)+ext*rvec(3)
     &         ,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)

          if(IFLAGT==1) then
            ibint = ibint+1

            if(km==nz-1 .AND. kd==1) then
              cond1 = flag(im1,jm1,km1)==0 .AND. flag2(im2,jm2,2)==0
              cond2 = flag(im1,jm1,km1)<=0 .AND. flag2(im2,jm2,2)<=0
            elseif(km==2 .AND. kd==-1) then
              cond1 = flag(im1,jm1,km1)==0 .AND. flag2(im2,jm2,1)==0 
              cond2 = flag(im1,jm1,km1)<=0 .AND. flag2(im2,jm2,1)<=0
            else
              cond1 = flag(im1,jm1,km1)==0 .AND. flag(im2,jm2,km2)==0
              cond2 = flag(im1,jm1,km1)<=0 .AND. flag(im2,jm2,km2)<=0
            endif

c            if(im==101 .AND. jm==5 .AND. km==272) then
c              write(6,*) in,n,'grid_intr_p:',im,jm,km,im1,jm1,km1,im2,jm2,km2,cond1,cond2
c     &              ,flag(im1,jm1,km1),flag(im2,jm2,km2)
c            endif

            if(cond1) then
              INTRS = INTRS + 1
              ITRIANGLE(INTRS) = N
              XINT (INTRS) = XMT
              YINT (INTRS) = YMT
              ZINT (INTRS) = ZMT
              NXINT(INTRS) = id
              NYINT(INTRS) = jd
              NZINT(INTRS) = kd
              icand(intrs) = 0
              grid_intr_p = 1
              go to 300
            elseif(cond2) then
              ncand = ncand+1
              INTRS = INTRS + 1
              ITRIANGLE(INTRS) = N
              XINT (INTRS) = XMT
              YINT (INTRS) = YMT
              ZINT (INTRS) = ZMT
              NXINT(INTRS) = id
              NYINT(INTRS) = jd
              NZINT(INTRS) = kd
              icand(intrs) = 1
c              grid_intr_p = -2
              go to 300
            else
              grid_intr_p = 0
            endif
          endif

        enddo

c 300    continue
        
      enddo

 300  continue
c
c.....Shoot rays in x+,z+ direction
c      if(im==88 .AND. jm==3 .AND. km==224) then
c        write(6,*) im,jm,km,intrs
c        do i=1,intrs
c          write(6,*) i,xint(i),yint(i),zint(i),nxint(i),nyint(i),nzint(i),icand(i)
c        enddo
c      endif

      IF(INTRS==0) THEN
        if(ibint==0) then
          grid_intr_p = -1
c        elseif(ncand>0) then
c          grid_intr_p = -2
        endif
      ELSE
        grid_intr_p = 1
        IF(INTRS>6) THEN
          WRITE(6,*) 'PROBLEM: MORE THAN 6 INTERSECTIONS'
          INTRS=2
        ENDIF

        DO nitrs=1,intrs
          ii = itriangle(nitrs)
          id = int(nxint(nitrs))
          jd = int(nyint(nitrs))
          kd = int(nzint(nitrs))
          ds(nitrs)= stencil_normal_unitdot(
     &         xu(im+id,jm+jd),yu(im+id,jm+jd),zu(km+kd)
     &         ,xu(im,jm),yu(im,jm),zu(km)
     &         ,unvect(1,ii),unvect(2,ii),unvect(3,ii))
c          DS(NITRS) = SQRT( (XU(IM,JM)-XINT(NITRS))**2. 
c     &          +  (YU(IM,JM)-YINT(NITRS))**2.
c     &          +  (ZU(KM)-ZINT(NITRS))**2. )
        ENDDO
        ii = maxloc(ds(1:intrs),1)       
c        if(im==88 .AND. jm==3 .AND. km==224) then
c          write(6,*) ii,ds(1:intrs)
c        endif
        XNIM = XINT(II)
        YNIM = YINT(II)
        ZNIM = ZINT(II)
        NXIM = REAL(NXINT(II))
        NYIM = REAL(NYINT(II))
        NZIM = REAL(NZINT(II))
        itr = itriangle(ii)

        if(icand(ii)==1) grid_intr_p=-2
      ENDIF


      return
      end
C---------------------------------------------------------------------------




C---- function norm_intr------------------------N. Beratlis-05 Jan. 2009--
c
      integer function norm_intr(flag,nx,ny,nz,ibd,nbd,vertex,unvect,nfacet
     &          ,xu,yu,zu,inn,nn,xnim,ynim,znim,nxim,nyim,nzim,fp,im,jm,km,dir,io,itri)
C
C     PURPOSE: Find normal intersection of a forcing point with immersed body.
C     RETURNS:
C      1 if inters. with physical stencil is found (SUCCESS)
C      0 if no inters. with immersed body is found.
C     -1 if intrs. with face is found and face contains forcing pts but no body pts.
C     -2 if intrs. with face is found and face contains body pts.
C     -3 physical stencil but no intersection with face
C
C-------------------------------------------------------------------------
C
c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'
c
c.... Input/Output Arrays
      INTEGER nx,ny,nz,nbd,nfacet,im,jm,km,nn,ibd,io,dir,fp,itri
      REAL    xnim,ynim,znim,nxim,nyim,nzim
      INTEGER flag(nx,ny,nz)
      INTEGER inn(nn)
      REAL    xu(nx,ny),yu(nx,ny),zu(nz)
      REAL    vertex(3,3,nfacet),unvect(3,nfacet)
c
c.... Local arrays
      INTEGER in,n,iflagt,ibint,intrs,nitrs,i,j,ii
      REAL    xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,xinp,yinp,zinp,ext
      REAL    a,theta,dz,jy,rm
      INTEGER iext,jext,kext,iext1,jext1,kext1
      integer ifint1,iphf1,ifint2,iphf2,ifint3,iphf3,i_no_fint
      INTEGER itriangle(nn),idir(nn)
      REAL    rvec(3),q(3)
      REAL    xint(nn),yint(nn),zint(nn),nxint(nn),nyint(nn),nzint(nn)
     &       ,ds(nn)
      INTEGER icand(nn)
c      INTEGER ii(nn),jj(nn),kk(nn)
c
c.... Functions
      INTEGER fluidface,ray_face_int
      real    extmag,anglerad
c
      norm_intr = 0

      ext = 10.0
      ibint = 0
      intrs = 0
      i_no_fint = 0
      itri = 0

      q(1) = xu(im,jm)
      q(2) = yu(im,jm)
      q(3) = zu(km)

      jy = 0.0
      if(icyl==1) then
        jy = -1.5
        if(yu(2,1)==0.0) jy=-1.0
      endif
      theta = dely*(real(jm)+jy)

      DO in=1,nn
        n=inn(in)
        xa = vertex(1,1,n)
        ya = vertex(2,1,n)
        za = vertex(3,1,n)
        xb = vertex(1,2,n)
        yb = vertex(2,2,n)
        zb = vertex(3,2,n)
        xc = vertex(1,3,n)
        yc = vertex(2,3,n)
        zc = vertex(3,3,n)
        rvec(1:3)=-unvect(1:3,n)
           
        ifint1 = 0
        ifint2 = 0
        ifint3 = 0
        iphf1 = 0
        iphf2 = 0
        iphf3 = 0
c in xmt, ymt and zmt the coordinates of the intersection points with
c the immersed-boundary are stored
        call segtriint(xu(im,jm)-ext*rvec(1),yu(im,jm)-ext*rvec(2),zu(km)-ext*rvec(3)
     &       ,xu(im,jm)+ext*rvec(1),yu(im,jm)+ext*rvec(2),zu(km)+ext*rvec(3)
     &       ,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)

        if(iflagt==1) then

          if(ibint==0) itri = in
          ibint = ibint+1
 
          rvec(1:3) = abs(ext)*unvect(1:3,n)
          fp = n

          rm = sqrt(xu(im-1,jm)**2. + yu(im-1,jm)**2.)
          if(icyl==1 .AND. im==2) rm=0.0
c the outward extensions iext, jext and kext from the interface points are
c found (in terms of grid indices)
          call vec_ijkext_gridpoint(q,rvec,iext,jext,kext,icyl,1,rm,dely) !1,-1
c
c the coordinates xinp, yinp and zinp of the intersection of the
c outward vector with a grid face are found
          ifint1 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &         ,im+iext,jm+jext,km+kext,1,xinp,yinp,zinp,icyl)          

          if(ifint1==1) then
c
c fluidface establishes if the grid face is physical
c 1: all points are fluid (OK)
c 0: there is at least one boundary point
c -1: there is at least one interior point
c
            iphf1 = fluidface(flag,nx,ny,nz,im+iext,jm+jext,km+kext,1)

            if(iphf1==1) then
              INTRS = INTRS+1
              IDIR(INTRS)=1
              ITRIANGLE(INTRS)=N
              XINT (INTRS) = XMT
              YINT (INTRS) = YMT
              ZINT (INTRS) = ZMT
              NXINT(INTRS) = xinp
              NYINT(INTRS) = yinp
              NZINT(INTRS) = zinp
              icand(intrs) = 0
              go to 401
            elseif(iphf1==0) then
              norm_intr = -1
              INTRS = INTRS+1
              IDIR(INTRS)=1
              ITRIANGLE(INTRS)=N
              XINT (INTRS) = XMT
              YINT (INTRS) = YMT
              ZINT (INTRS) = ZMT
              NXINT(INTRS) = xinp
              NYINT(INTRS) = yinp
              NZINT(INTRS) = zinp
              icand(intrs) = 1
              go to 401
            elseif(iphf1==-1) then
              norm_intr = -2
            endif

          else
c
c if ray_face_int is not able to find an intersection with the grid faces
c the direction along which the code looks for an intersection is modified
c from 1 (x) to 2 (y)
c
            call vec_ijkext_gridpoint(q,unvect(1:3,n),iext,jext,kext,icyl,2,theta,dely) !1,-1

            ifint2 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &            ,im+iext,jm+jext,km+kext,2,xinp,yinp,zinp,icyl)

            if(ifint2==1) then

              iphf2 = fluidface(flag,nx,ny,nz,im+iext,jm+jext,km+kext,2)

              if(iphf2==1) then
                INTRS = INTRS + 1
                IDIR(INTRS)=2
                ITRIANGLE(INTRS)=N
                XINT (INTRS) = XMT
                YINT (INTRS) = YMT
                ZINT (INTRS) = ZMT
                NXINT(INTRS) = xinp
                NYINT(INTRS) = yinp
                NZINT(INTRS) = zinp
                icand(intrs) = 0
                go to 401
              elseif(iphf2==0) then
                INTRS = INTRS + 1
                IDIR(INTRS)=2
                ITRIANGLE(INTRS)=N
                XINT (INTRS) = XMT
                YINT (INTRS) = YMT
                ZINT (INTRS) = ZMT
                NXINT(INTRS) = xinp
                NYINT(INTRS) = yinp
                NZINT(INTRS) = zinp
                norm_intr = -1
                icand(intrs) = 1
                go to 401
              elseif(iphf2==-1) then
                norm_intr = -2
              endif
            else

              dz = zu(km + int(sign(1.0,unvect(3,n))))-zu(km)
              call vec_ijkext_gridpoint(q,unvect(1:3,n),iext,jext,kext,icyl,3,dz,dely) !1,-1

              ifint3 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &              ,im+iext,jm+jext,km+kext,3,xinp,yinp,zinp,icyl)

              if(ifint3==1) then

                iphf3 = fluidface(flag,nx,ny,nz,im+iext,jm+jext,km+kext,3)

                if(iphf3==1) then
                  INTRS = INTRS + 1
                  IDIR(INTRS)=3
                  ITRIANGLE(INTRS)=N
                  XINT (INTRS) = XMT
                  YINT (INTRS) = YMT
                  ZINT (INTRS) = ZMT
                  NXINT(INTRS) = xinp
                  NYINT(INTRS) = yinp
                  NZINT(INTRS) = zinp
                  icand(intrs) = 0
                  go to 401
                elseif(iphf3==0) then
                  INTRS = INTRS + 1
                  IDIR(INTRS)=3
                  ITRIANGLE(INTRS)=N
                  XINT (INTRS) = XMT
                  YINT (INTRS) = YMT
                  ZINT (INTRS) = ZMT
                  NXINT(INTRS) = xinp
                  NYINT(INTRS) = yinp
                  NZINT(INTRS) = zinp
                  norm_intr = -1
                  icand(intrs) = 1
                  go to 401
                elseif(iphf3==-1) then
                  norm_intr = -2
                endif

              endif
            endif
          endif

        endif

      ENDDO

 401  continue    

      if(intrs==0) then
        if(ibint==0) then
          norm_intr = 0
        else
          if(i_no_fint>0) then
            write(6,*) 'i_no_fint=',i_no_fint
            norm_intr =-3
          endif
        endif
      else
        norm_intr = 1

c
c.... Find fluid with vector most aligned with normal
        DO NITRS=1,INTRS
c          IF(ICAND(NITRS)==0) THEN
c            II = NITRS
c            exit
c          ELSE
c            II = 1
c          ENDIF
c distance of the boundary point from the intersections with the immersed
c boundary
          DS(NITRS) = SQRT( (XU(IM,JM)-XINT(NITRS))**2. 
     &        +  (YU(IM,JM)-YINT(NITRS))**2.
     &        +  (ZU(KM)-ZINT(NITRS))**2. )
        ENDDO
c among the intersection points with the immersed boundary the one closest
c to the interface node is chosen
        II = MINLOC(DS(1:intrs),1)

        xnim = xint(ii)
        ynim = yint(ii)
        znim = zint(ii)
        nxim = REAL(nxint(ii))
        nyim = REAL(nyint(ii))
        nzim = REAL(nzint(ii))
        dir = idir(ii)

        if(icand(ii)==1) then
          norm_intr=-1
        endif
      endif

      return

      end
C---------------------------------------------------------------------------


C---- function closest_intr------------------------N. Beratlis-05 Jan. 2009--
c
      integer function closest_intr(flag,nx,ny,nz,ibd,nbd,vertex,unvect,nfacet
     &          ,xu,yu,zu,inn,nn,xnim,ynim,znim,nxim,nyim,nzim,fp,im,jm,km,dir,io)
C
C     PURPOSE: Find intersection of a forcing point with closest point with immersed body.
C     RETURNS:
C      1 if inters. with physical stencil is found (SUCCESS)
C      0 if no inters. with immersed body is found.
C     -1 if intrs. with face is found and face contains forcing pts but no body pts.
C     -2 if intrs. with face is found and face contains body pts.
C     -3 physical stencil but no intersection with face
C
C-------------------------------------------------------------------------
C
c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'
c
c.... Input/Output Arrays
      INTEGER nx,ny,nz,nbd,nfacet,im,jm,km,nn,ibd,io,dir,fp
      REAL    xnim,ynim,znim,nxim,nyim,nzim
      INTEGER flag(nx,ny,nz)
      INTEGER inn(nn)
      REAL    xu(nx,ny),yu(nx,ny),zu(nz)
      REAL    vertex(3,3,nfacet),unvect(3,nfacet)
c
c.... Local arrays
      INTEGER in,n,iflagt,ibint,intrs,nitrs,i,j,ii
      REAL    xmt,ymt,zmt,xinp,yinp,zinp,ext
      REAL    a,theta,dz,jy,rm
      INTEGER iext,jext,kext,iext1,jext1,kext1
      integer ifint1,iphf1,ifint2,iphf2,ifint3,iphf3,i_no_fint
      INTEGER itriangle(nn),idir(nn),iphf(nn)
      REAL    rvec(3),q(3),p(3)
      REAL    xint(nn),yint(nn),zint(nn),nxint(nn),nyint(nn),nzint(nn)
     &       ,ds(nn)
      INTEGER icand(nn)
c      INTEGER ii(nn),jj(nn),kk(nn)
c
c.... Functions
      INTEGER fluidface,ray_face_int
      real    extmag,anglerad,dotproduct,stencil_normal_unitdot,vecmag
c
      closest_intr = 0

      ext = 100.0
      ibint = 0
      intrs = 0
      i_no_fint = 0
      iphf = -1

      q(1) = xu(im,jm)
      q(2) = yu(im,jm)
      q(3) = zu(km)

      jy = 0.0
      if(icyl==1) then
        jy = -1.5
        if(yu(2,1)==0.0) jy=-1.0
      endif
      theta = dely*(real(jm)+jy)

      DO in=1,nn

        n=inn(in)
           
        ifint1 = 0
        ifint2 = 0
        ifint3 = 0
        iphf1 = 0
        iphf2 = 0
        iphf3 = 0

c the closest point to the triangle N is found and its coordinates are stored in P
        call tri_closest_pt(vertex(:,:,n),q,p)
        iflagt = 1

        if(iflagt==1) then

          ibint = ibint+1          

          rvec = q-p
          rvec = abs(ext)*rvec/vecmag(rvec,3)
          if(dotproduct(q-p,unvect(:,n),3)<0) rvec = -rvec
          fp = n

c the point P is taken as intersection point with the immersed-boundary surface
          xmt = p(1)
          ymt = p(2)
          zmt = p(3)

          rm = sqrt(xu(im-1,jm)**2. + yu(im-1,jm)**2.)
          if(icyl==1 .AND. im==2) rm=0.0
c the outward extensions iext, jext and kext from the interface point Q are
c found (in terms of grid indices)
          call vec_ijkext_gridpoint(q,rvec,iext,jext,kext,icyl,1,rm,dely) !1,-1
c
c the coordinates xinp, yinp and zinp of the intersection of the
c outward vector with a grid face are found
          ifint1 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &         ,im+iext,jm+jext,km+kext,1,xinp,yinp,zinp,icyl)          

          if(ifint1==1) then
c
c fluidface establishes if the grid face is physical
c 1: all points are fluid (OK)
c 0: there is at least one boundary point
c -1: there is at least one interior point
c
            iphf1 = fluidface(flag,nx,ny,nz,im+iext,jm+jext,km+kext,1)

            if(iphf1==1) then
              INTRS = INTRS+1
              IPHF(INTRS) = 1
              IDIR(INTRS)= 1
              ITRIANGLE(INTRS)=N
              XINT (INTRS) = XMT
              YINT (INTRS) = YMT
              ZINT (INTRS) = ZMT
              NXINT(INTRS) = xinp
              NYINT(INTRS) = yinp
              NZINT(INTRS) = zinp
              icand(intrs) = 0
              go to 401
            elseif(iphf1==0) then
              closest_intr = -1
              INTRS = INTRS+1
              IPHF(INTRS) = 0
              IDIR(INTRS)=1
              ITRIANGLE(INTRS)=N
              XINT (INTRS) = XMT
              YINT (INTRS) = YMT
              ZINT (INTRS) = ZMT
              NXINT(INTRS) = xinp
              NYINT(INTRS) = yinp
              NZINT(INTRS) = zinp
              icand(intrs) = 1
              go to 401
            elseif(iphf1==-1) then
              closest_intr = -2
            endif

          else
c
c if ray_face_int is not able to find an intersection with the grid faces
c the direction along which the code looks for an intersection is modified
c from 1 (x) to 2 (y)
c
            call vec_ijkext_gridpoint(q,rvec,iext,jext,kext,icyl,2,theta,dely) !1,-1

            ifint2 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &            ,im+iext,jm+jext,km+kext,2,xinp,yinp,zinp,icyl)

            if(ifint2==1) then

              iphf2 = fluidface(flag,nx,ny,nz,im+iext,jm+jext,km+kext,2)

              if(iphf2==1) then
                INTRS = INTRS + 1
                IPHF(INTRS) = 1     
                IDIR(INTRS)=2
                ITRIANGLE(INTRS)=N
                XINT (INTRS) = XMT
                YINT (INTRS) = YMT
                ZINT (INTRS) = ZMT
                NXINT(INTRS) = xinp
                NYINT(INTRS) = yinp
                NZINT(INTRS) = zinp
                icand(intrs) = 0
                go to 401
              elseif(iphf2==0) then
                INTRS = INTRS + 1
                IPHF(INTRS) = 0
                IDIR(INTRS)=2
                ITRIANGLE(INTRS)=N
                XINT (INTRS) = XMT
                YINT (INTRS) = YMT
                ZINT (INTRS) = ZMT
                NXINT(INTRS) = xinp
                NYINT(INTRS) = yinp
                NZINT(INTRS) = zinp
                closest_intr = -1
                icand(intrs) = 1
                go to 401
              elseif(iphf2==-1) then
                closest_intr = -2
              endif
            else

              dz = zu(km + int(sign(1.0,rvec(3))))-zu(km)
              call vec_ijkext_gridpoint(q,rvec,iext,jext,kext,icyl,3,dz,dely) !1,-1

              ifint3 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &              ,im+iext,jm+jext,km+kext,3,xinp,yinp,zinp,icyl)

              if(ifint3==1) then

                iphf3 = fluidface(flag,nx,ny,nz,im+iext,jm+jext,km+kext,3)

                if(iphf3==1) then
                  INTRS = INTRS + 1
                  IPHF(INTRS) = 1 
                  IDIR(INTRS)=3
                  ITRIANGLE(INTRS)=N
                  XINT (INTRS) = XMT
                  YINT (INTRS) = YMT
                  ZINT (INTRS) = ZMT
                  NXINT(INTRS) = xinp
                  NYINT(INTRS) = yinp
                  NZINT(INTRS) = zinp
                  icand(intrs) = 0
                  go to 401
                elseif(iphf3==0) then
                  INTRS = INTRS + 1
                  IPHF(INTRS) = 0
                  IDIR(INTRS)=3
                  ITRIANGLE(INTRS)=N
                  XINT (INTRS) = XMT
                  YINT (INTRS) = YMT
                  ZINT (INTRS) = ZMT
                  NXINT(INTRS) = xinp
                  NYINT(INTRS) = yinp
                  NZINT(INTRS) = zinp
                  closest_intr = -1
                  icand(intrs) = 1
                  go to 401
                elseif(iphf3==-1) then
                  closest_intr = -2
                endif

              endif
            endif
          endif

        endif

 401  continue

      ENDDO

c 401  continue
      if(intrs==0) then
        if(ibint==0) then
          closest_intr = 0
        else
          if(i_no_fint>0) then
            write(6,*) 'i_no_fint=',i_no_fint
            closest_intr =-3
          else
            closest_intr =-2
          endif
        endif
      else
        closest_intr = 1
        DO NITRS=1,INTRS
          IF(ICAND(NITRS)==0) THEN
c the first case with a valid fluid face is selected
            II = NITRS
            exit
          ELSE
            II = 1
          ENDIF
c           a = stencil_normal_unitdot(nxint(nitrs),nyint(nitrs)
c     $          ,nzint(nitrs),xu(im,jm),yu(im,jm),zu(km),unvect(1,n)
c     $          ,unvect(2,n),unvect(3,n))
c          DS(NITRS) = a+iphf(nitrs)*nitrs
cc          DS(NITRS) = SQRT( (XU(IM,JM)-XINT(NITRS))**2. 
cc     &        +  (YU(IM,JM)-YINT(NITRS))**2.
cc     &        +  (ZU(KM)-ZINT(NITRS))**2. )
        ENDDO
cc        II = MINLOC(DS(1:intrs),1)
c        II = MAXLOC(DS(1:intrs),1)
c        write(6,*) '3. closest intr:',im,jm,km,ibint,closest_intr,intrs,icand(1:intrs)
        xnim = xint(ii)
        ynim = yint(ii)
        znim = zint(ii)
        nxim = REAL(nxint(ii))
        nyim = REAL(nyint(ii))
        nzim = REAL(nzint(ii))
        dir = idir(ii)

        if(icand(ii)==1) then
          closest_intr=-1
        endif
      endif

      return

      end
C---------------------------------------------------------------------------





C---- function norm_intr_p----------------------N. Beratlis-06 Jan. 2009--
C
C     PURPOSE: Find normal intersection of a forcing point with immersed 
C     body and 2 physical points in the fluid along the normal vector.
C     RETURNS:
C     1 - Intrs. with body is found and intrs. face contain only fluid pts (SUCCESS)   
C     0 - No normal inters. with immers. body is found
C    -1 - If intrs. with body is found and intrs. faces contain forcing pts but no body pts.
C    -2 - If intrs. with body is found but intrs. faces contain body pts.
C
C-------------------------------------------------------------------------
      integer function norm_intr_p(flag,flag2,nx,ny,nz,ibd,nbd
     &        ,vertex,unvect,nfacet,xu,yu,zu,zug,nzg,inn,nn
     &        ,xnim,ynim,znim,nxim,nyim,nzim,im,jm,km,dir,itr,icntl,nt,clock,nclocks)
C
      IMPLICIT NONE
      include 'immersed.h'

      INTEGER nx,ny,nz,nzg,nbd,nfacet,im,jm,km,nn,ibd,itr,icntl,nt,nclocks
      REAL    xnim,ynim,znim
      REAL    clock(nclocks)
      REAL    nxim(2),nyim(2),nzim(2)
      INTEGER flag(nx,ny,nz),flag2(nx,ny,2)
      INTEGER inn(nn),dir(2)
      REAL    xu(nx,ny),yu(nx,ny),zu(nz),zug(nzg)
      REAL    vertex(3,3,nfacet),unvect(3,nfacet)
c
      INTEGER i,in,n,iflagt,iflag,ibint,intrs,ii,nitrs
      REAL    xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,ext
      REAL    xinp,yinp,zinp,xinp2,yinp2,zinp2
      INTEGER iext,jext,kext,iext1,jext1,kext1,ncand
      INTEGER interp_points_pres
      INTEGER itriangle(nn),idir(nn,2),iplane(nn,2)
      REAL    rvec(3),q(3),q2(3)
      INTEGER icand(nn)
      REAL    xint(nn),yint(nn),zint(nn),nxint(nn,2),nyint(nn,2),nzint(nn,2)
     &       ,ds(nn)
c      REAL    clock(2)
      REAL    clocktemp
c
c.... Functions
      REAL    tclock
c
      ext = 10.0
      ibint = 0
      intrs = 0
      norm_intr_p=0
      ncand = 0
      icand = 0

      DO in=1,nn
        nt = nt+1
        n=inn(in)
        xa = vertex(1,1,n)
        ya = vertex(2,1,n)
        za = vertex(3,1,n)
        xb = vertex(1,2,n)
        yb = vertex(2,2,n)
        zb = vertex(3,2,n)
        xc = vertex(1,3,n)
        yc = vertex(2,3,n)
        zc = vertex(3,3,n)
        rvec(1:3)=unvect(1:3,n)
  
c        clocktemp = tclock()
        call segtriint(xu(im,jm)-ext*rvec(1),yu(im,jm)-ext*rvec(2),zu(km)-ext*rvec(3)
     &       ,xu(im,jm)+ext*rvec(1),yu(im,jm)+ext*rvec(2),zu(km)+ext*rvec(3)
     &       ,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)
c        clock(1) = clock(1) + tclock() - clocktemp

        if(iflagt==1) then

          ibint = ibint+1
          q(1) = xu(im,jm)
          q(2) = yu(im,jm)
          q(3) = zu(km)
c          clocktemp = tclock()
          iflag = interp_points_pres(flag,flag2,xu,yu,zu,zug,rvec,nx,ny,nz,nzg
     &      ,xinp,yinp,zinp,xinp2,yinp2,zinp2,dir,im,jm,km,icntl,clock(3),3)
c          iflag=0
c          clock(2) = clock(2) + tclock() - clocktemp
          if(iflag>0) then
            norm_intr_p = 1
            intrs = intrs+1
            xint(intrs) = xmt
            yint(intrs) = ymt
            zint(intrs) = zmt
            nxint(intrs,1) = xinp
            nyint(intrs,1) = yinp
            nzint(intrs,1) = zinp
            nxint(intrs,2) = xinp2
            nyint(intrs,2) = yinp2
            nzint(intrs,2) = zinp2
            idir(intrs,:) = dir(:)
            itriangle(intrs) = n
            icand(intrs) = 0
            go to 402
          elseif(iflag==0) then
            ncand = ncand+1
            norm_intr_p = -1
            intrs = intrs+1
            xint(intrs) = xmt
            yint(intrs) = ymt
            zint(intrs) = zmt
            nxint(intrs,1) = xinp
            nyint(intrs,1) = yinp
            nzint(intrs,1) = zinp
            nxint(intrs,2) = xinp2
            nyint(intrs,2) = yinp2
            nzint(intrs,2) = zinp2
            idir(intrs,:) = dir(:)
            itriangle(intrs) = n
            icand(intrs) = 1
            go to 402
          else
            norm_intr_p = -2
          endif
        endif

      enddo

 402  continue

      if(intrs==0) then
        if(ibint==0) then
          norm_intr_p = 0
        endif
      elseif(intrs>0) then
        do i=1,intrs
          ds(i) = sqrt( (xu(im,jm)-xint(i))**2. 
     &                + (yu(im,jm)-yint(i))**2.
     &                + (zu(km)-zint(i))**2. )
        enddo

        ii = minloc(ds(1:intrs),1)

        xnim = xint(ii)
        ynim = yint(ii)
        znim = zint(ii)
        nxim(:) = nxint(ii,:)
        nyim(:) = nyint(ii,:)
        nzim(:) = nzint(ii,:)
        dir(:) = idir(ii,:)
        itr = itriangle(ii)

        if(icand(ii)==1) norm_intr_p = -1
c        write(6,*) ii,itr,itriangle(1:intrs)
      endif

      return
      end
C-----------------------------------------------------------------------


C---- function diag_intr_p---------------------N. Beratlis-06 Jan. 2009--
C
C     PURPOSE: Find normal intersection of a forcing point with immersed 
C     body and 2 physical points in the fluid along the normal vector.
C     OUTPUT:
C     1 - Intersection with immersed body and stencil contains physical pts (SUCCESS)
C     0 - No normal inters. with immers. body is found
C    -1 - If intrs. with body is found and intrs. faces contain forcing pts but no body pts.
C    -2 - If intrs. with body is found but intrs. faces contain body pts.
C
C------------------------------------------------------------------------
      integer function diag_intr_p(flag,flag2,nx,ny,nz,ibd,nbd
     &        ,vertex,unvect,nfacet,xu,yu,zu,zug,nzg,inn,nn
     &        ,xnim,ynim,znim,nxim,nyim,nzim,im,jm,km,dir,itr)
C
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'

      INTEGER nx,ny,nz,nzg,nbd,nfacet,im,jm,km,nn,ibd,itr
      REAL    xnim,ynim,znim
      REAL    nxim(2),nyim(2),nzim(2)
      INTEGER flag(nx,ny,nz),flag2(nx,ny,2)
      INTEGER inn(nn),dir(2)
      REAL    xu(nx,ny),yu(nx,ny),zu(nz),zug(nzg)
      REAL    vertex(3,3,nfacet),unvect(3,nfacet)
c
c.... Local arrays
      INTEGER i,in,n,iflagt,iflag,ibint,intrs,ii,nitrs,dirn
      INTEGER imd,jmd,kmd,kmdg,imd2,jmd2,kmd2,iph1
      INTEGER ifint(3),ifph(3)
      INTEGER flage(nx,ny,2)
      REAL    xa,ya,za,xb,yb,zb,xc,yc,zc,ext,rm,theta,dz
      REAL    xmt,ymt,zmt,xinp,yinp,zinp
      INTEGER iext,jext,kext,nd,ncand,nphy,nfor,nbdy,nphy1,nfor1,nbdy1
      INTEGER itriangle(nn),idir(nn,1)
      INTEGER id(20,3),icand(20),crit(20),ord(20)
      REAL    rvec(3),q(3),q2(3),diag(3)
      REAL    xint(nn),yint(nn),zint(nn),nxint(nn,2),nyint(nn,2),nzint(nn,2)
     &       ,ds(nn)
c
c.... Functions
      INTEGER ray_face_int,fluidface

      ext = 10.0
      ifint = 0
      iph1 = 0
      nphy = 0
      nfor = 0
      nbdy = 0
      iflag = 0
      ibint = 0
      intrs = 0
      ncand = 0
      diag_intr_p=0

      nd = 4

      !Define diagonal extensions
      id(1,1) = 1
      id(1,2) = 0
      id(1,3) = 1

      id(2,1) =-1
      id(2,2) = 0
      id(2,3) = 1

      id(3,1) =-1
      id(3,2) = 0
      id(3,3) =-1

      id(4,1) = 1
      id(4,2) = 0
      id(4,3) =-1
      
      id(5,1) = 1
      id(5,2) = 1
      id(5,3) = 1

      id(6,1) =-1
      id(6,2) = 1
      id(6,3) = 1

      id(7,1) =-1
      id(7,2) = 1
      id(7,3) =-1

      id(8,1) = 1
      id(8,2) = 1
      id(8,3) =-1
      
      id(9,1) = 1
      id(9,2) =-1
      id(9,3) = 1

      id(10,1) =-1
      id(10,2) =-1
      id(10,3) = 1

      id(11,1) =-1
      id(11,2) =-1
      id(11,3) =-1

      id(12,1) = 1
      id(12,2) =-1
      id(12,3) =-1
      
      id(13,1) = 1
      id(13,2) = 1
      id(13,3) = 0

      id(14,1) =-1
      id(14,2) = 1
      id(14,3) = 0

      id(15,1) = 0
      id(15,2) = 1
      id(15,3) = 1

      id(16,1) = 0
      id(16,2) = 1
      id(16,3) =-1

      id(17,1) = 1
      id(17,2) =-1
      id(17,3) = 0

      id(18,1) =-1
      id(18,2) =-1
      id(18,3) = 0

      id(19,1) = 0
      id(19,2) =-1
      id(19,3) = 1

      id(20,1) = 0
      id(20,2) =-1
      id(20,3) =-1
     
      !Determine order of grid ray tracing
      crit = 0
      do i=1,nd
        imd = im+id(i,1)
        jmd = jm+id(i,2)
        kmd = km+id(i,3)
        call index_bnd(imd,jmd,kmd,imd,jmd,kmd,nx,ny,nz)
        !Give priority to stencil with fluid points
        if(flag(imd,jmd,kmd)==0) then
          crit(i) = 1
        endif
      enddo

      do i=1,nd
        ord(i) = maxloc(crit(1:nd),1)
        crit(ord(i)) = -1
      enddo


      DO i=1,nd

        ii = ord(i)
        nfor = 0
        nbdy = 0
        nphy = 0

        imd = im+id(ii,1)
        jmd = jm+id(ii,2)
        kmd = km+id(ii,3)
        kmdg = kmd+myrank*(nz-2)

c        call index_bnd(imd,jmd,kmd,imd,jmd,kmd,nx,ny,nz)

        diag(1) = xu(imd,jmd)-xu(im,jm)
        diag(2) = yu(imd,jmd)-yu(im,jm)
        diag(3) = zu(kmd)-zu(km)
        
        if(flag(imd,jmd,kmd)>0) then
          nbdy = nbdy+1
        else

          if(flag(imd,jmd,kmd)==0) then
            nphy = nphy+1
          elseif(flag(imd,jmd,kmd)<0) then
            nfor = nfor+1
          endif
          iph1 = iph1+1

c          write(6,*) 'diag:',im,jm,km,imd,jmd,kmd,nn

          nfor1 = nfor
          nbdy1 = nbdy
          nphy1 = nphy

          DO in=1,nn

            nbdy = nbdy1
            nfor = nfor1
            nphy = nphy1
             
            n=inn(in)
            xa = vertex(1,1,n)
            ya = vertex(2,1,n)
            za = vertex(3,1,n)
            xb = vertex(1,2,n)
            yb = vertex(2,2,n)
            zb = vertex(3,2,n)
            xc = vertex(1,3,n)
            yc = vertex(2,3,n)
            zc = vertex(3,3,n)

            rvec=diag
  
            call segtriint(xu(im,jm)-ext*rvec(1),yu(im,jm)-ext*rvec(2),zu(km)-ext*rvec(3)
     &           ,xu(im,jm)+ext*rvec(1),yu(im,jm)+ext*rvec(2),zu(km)+ext*rvec(3)
     &           ,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)

            if(iflagt==1) then

              ibint = ibint+1
              q(1) = xu(imd,jmd)
              q(2) = yu(imd,jmd)
              q(3) = zu(kmd)

              do dirn=1,3

                if(dirn==1) then
c                  if(imd<=1) then
c                    write(6,*) imd,jmd,kmd
c                  endif
                  if(icyl==1 .AND. imd<=1) then 
                    rm=0.0
                  else
                    rm = sqrt(xu(imd-1,jmd)**2. + yu(imd-1,jmd)**2.)
                  endif
                  call vec_ijkext_gridpoint(q,rvec,iext,jext,kext,icyl,dirn,rm,dely)
                elseif(dirn==2) then
                  theta = dely*(real(jmd)-1.5)
                  call vec_ijkext_gridpoint(q,rvec,iext,jext,kext,icyl,dirn,theta,dely)
                else
                  dz = zug(kmdg + int(sign(1.0,rvec(3))))-zug(kmdg)
                  call vec_ijkext_gridpoint(q,rvec,iext,jext,kext,icyl,dirn,dz,dely)
                endif

                call index_bnd(imd+iext,jmd+jext,kmd,imd2,jmd2,kmd,nx,ny,nz)

                if(kmd==1 .AND. kext<0) then
                  ifint(dirn) = ray_face_int(q,rvec,xu,yu,zug,nx,ny,nzg
     &                  ,imd2,jmd2,kmdg+kext,dirn,xinp,yinp,zinp,icyl)
                  flage(:,:,1) = flag2(:,:,1)
                  flage(:,:,2) = flag(:,:,1)
                  ifph(dirn) = fluidface(flage,nx,ny,2,imd2,jmd2,1,dirn)
                elseif(kmd==nz) then
                  ifint(dirn) = ray_face_int(q,rvec,xu,yu,zug,nx,ny,nzg
     &                  ,imd2,jmd2,kmdg+kext,dirn,xinp,yinp,zinp,icyl)
                  flage(:,:,1) = flag(:,:,nz)
                  flage(:,:,2) = flag2(:,:,2)
                  ifph(dirn) = fluidface(flage,nx,ny,2,imd2,jmd2,1,dirn)
                else
                  ifint(dirn) = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &                  ,imd2,jmd2,kmd+kext,dirn,xinp,yinp,zinp,icyl)
                  ifph(dirn) = fluidface(flag,nx,ny,nz,imd2,jmd2,kmd+kext,dirn)
                endif

                if(ifint(dirn)==1) then

                  if(ifph(dirn)==1) then
                    if(nphy==1) then
                      nphy = nphy+1
                      intrs = intrs+1
                      xint(intrs) = xmt
                      yint(intrs) = ymt
                      zint(intrs) = zmt
                      nxint(intrs,1) = real(id(ii,1))
                      nyint(intrs,1) = real(id(ii,2))
                      nzint(intrs,1) = real(id(ii,3))
                      nxint(intrs,2) = xinp
                      nyint(intrs,2) = yinp
                      nzint(intrs,2) = zinp
                      idir(intrs,1) = dirn
                      itriangle(intrs) = n
                      icand(intrs) = 0
                      goto 350
                    else
                      ncand = ncand+1
                      intrs = intrs+1
                      xint(intrs) = xmt
                      yint(intrs) = ymt
                      zint(intrs) = zmt
                      nxint(intrs,1) = real(id(i,1))
                      nyint(intrs,1) = real(id(i,2))
                      nzint(intrs,1) = real(id(i,3))
                      nxint(intrs,2) = xinp
                      nyint(intrs,2) = yinp
                      nzint(intrs,2) = zinp
                      idir(intrs,1) = dirn
                      itriangle(intrs) = n
                      icand(intrs) = 1
                      goto 350
                    endif
                  elseif(ifph(dirn)==0) then
                    ncand = ncand+1
                    intrs = intrs+1
                    xint(intrs) = xmt
                    yint(intrs) = ymt
                    zint(intrs) = zmt
                    nxint(intrs,1) = real(id(ii,1))
                    nyint(intrs,1) = real(id(ii,2))
                    nzint(intrs,1) = real(id(ii,3))
                    nxint(intrs,2) = xinp
                    nyint(intrs,2) = yinp
                    nzint(intrs,2) = zinp
                    idir(intrs,1) = dirn
                    itriangle(intrs) = n
                    icand(intrs) = 1
                    go to 350
                  endif
                  goto 349
                endif
              enddo
 349          continue

            endif
          enddo

        endif

c 350    continue

      enddo

 350  continue
c      if(im==86 .AND. jm==2 .AND. km==323) then
c        write(6,*) '2. diag_intr_p',im,jm,km,nphy,nfor,nbdy
c      endif

      if(intrs==0) then
        if(ibint==0) then
          diag_intr_p = 0
c        elseif(ncand>0) then
c          diag_intr_p = -1
        else
          diag_intr_p = -2
        endif
c        write(6,*) 'diagonal:',im,jm,km,intrs,ibint,iph1
      elseif(intrs>0) then
        diag_intr_p = 1
        do i=1,intrs
          ds(i) = sqrt( (xu(imd,jmd)-xint(i))**2. 
     &                + (yu(imd,jmd)-yint(i))**2.
     &                + (zu(kmd)-zint(i))**2. )
        enddo

        ii = minloc(ds(1:intrs),1)

        xnim = xint(ii)
        ynim = yint(ii)
        znim = zint(ii)
        nxim(:) = nxint(ii,:)
        nyim(:) = nyint(ii,:)
        nzim(:) = nzint(ii,:)
        dir(1) = idir(ii,1)
        itr = itriangle(ii)
        
        if(icand(ii)==1) diag_intr_p=-1
c        write(6,*) ii,itr,itriangle(1:intrs)
      endif

      return
      end
C-----------------------------------------------------------------------


C---- function interp_points_pres-------------N. Beratlis-06 Jan. 2009--
C
C     PURPOSE: Find 2 interpolations points for pressure in the fluid along
C     the normal intersection with immersred
C     OUTPUT: 
C      1 - Intrs. faces contains only fluid pts (SUCCESS)
C      0 - Intrs. faces contain forcing pts but no body pts
C     -1 - Intrs. faces contain body pts.
C
C-----------------------------------------------------------------------
      integer function interp_points_pres(flag,flag2,xu,yu,zu,zug,rvec
     &        ,nx,ny,nz,nzg,xinp,yinp,zinp,xinp2,yinp2,zinp2,idir
     &        ,im,jm,km,icntl,clock,nclocks)

      include 'common.h'

      INTEGER im,jm,km,nx,ny,nz,nzg,icntl,nclocks
      REAL    xinp,yinp,zinp,xinp2,yinp2,zinp2
      REAL    clock(nclocks)
      INTEGER idir(2)
      REAL    rvec(3)
      INTEGER flag(nx,ny,nz),flag2(nx,ny,2)
      REAL    xu(nx,ny),yu(nx,ny),zu(nz),zug(nzg)
c
c.... Input/Ouput Arrays
      INTEGER im1,jm1,km1,im2,jm2,km2,km1g,km2g
      INTEGER iext,jext,kext,iext1,jext1,kext1,iext2,jext2,kext2
      INTEGER km31,km31g,jm11
      integer ifph,ifph1,ifph2,ifph3,nphy,nfor,ncand
      integer  ifint1,ifint11,ifint12,ifint13
      integer  ifint2,ifint21,ifint22,ifint23
      integer  ifint3,ifint31,ifint32,ifint33
      real    a,theta,dz,rm,clocktemp
      INTEGER flagt(nx,ny,2)
      REAL    q(3)
c      REAL    clock(3)
c
c.... Functions
      INTEGER ray_face_int,fluidface
      real    extmag,anglerad,tclock

c      clock = 0.0

      interp_points_pres=1
      ifint1 = 0
      ifint2 = 0
      ifint3 = 0
      ifint11 = 0
      ifint12 = 0
      ifint13 = 0
      ifint21 = 0
      ifint22 = 0
      ifint23 = 0
      ifint31 = 0
      ifint32 = 0
      ifint33 = 0
      nphy = 0
      nfor = 0
      ncand = 0

      q(1) = xu(im,jm)
      q(2) = yu(im,jm)
      q(3) = zu(km)

      theta = dely*(real(jm)-1.5)

c      clocktemp = tclock()
      rm = sqrt(xu(im-1,jm)**2. + yu(im-1,jm)**2.)
      if(icyl==1 .AND. im==2) rm=0.0
      call vec_ijkext_gridpoint(q,rvec,iext,jext,kext,icyl,1,rm,dely) !1,-1
c      clock(1) = clock(1) + tclock() - clocktemp

c      clocktemp = tclock()
      ifint1 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &  ,im+iext,jm+jext,km+kext,1,xinp,yinp,zinp,icyl)
c      clock(2) = clock(2) + tclock() - clocktemp

c      if(im==2 .AND. jm==2 .AND. km==173) then
c         write(6,*) '1. intep_points_pres:',im,jm,km,ifint1
c      endif

      if(ifint1==1) then

c        clocktemp = tclock()
        ifph1 = fluidface(flag,nx,ny,nz,im+iext,jm+jext,km+kext,1)
c        clock(3) = clock(3) + tclock() - clocktemp

        if(ifph1==-1) then
          interp_points_pres=-1
          return
        else
          if(ifph1==1) then
            interp_points_pres=1
            nphy = 1
          else
            interp_points_pres=0
          endif
          idir(1) = 1

          q(1) = xinp
          q(2) = yinp
          q(3) = zinp

          !Set new indices im1,jm1,km1 for new q
c          clocktemp = tclock()
          im1 = im+iext
          jm1 = jm+jext
          km1 = km+kext
          call per_index(jm1,ny,jm1)

          rm = sqrt(xu(im1-1,jm)**2. + yu(im1-1,jm)**2.)
          if(icyl==1 .AND. im1==2) rm=0.0
          call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,1,rm,dely,1) !0,-1
c          clock(1) = clock(1) + tclock() - clocktemp

          im2 = im1+iext

c          clocktemp = tclock()
          ifint11 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &         ,im2,jm1,km1,1,xinp2,yinp2,zinp2,icyl)
c          clock(2) = clock(2) + tclock() - clocktemp

          if(ifint11==1) then

c            clocktemp = tclock()
            ifph = fluidface(flag,nx,ny,nz,im2,jm1,km1,1)
c            clock(3) = clock(3) + tclock() - clocktemp

            idir(2)=1
            if(ifph==1) then
              if(nphy==1) then
                interp_points_pres=1
                return                
              else
                interp_points_pres=0
                return
              endif
            elseif(ifph==0) then
              interp_points_pres=0
              return
            else
              interp_points_pres=-1
              return
            endif

          else

c            clocktemp = tclock()
            theta = dely*(real(jm1)-1.5)
            call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,2,theta,dely,1) !0,-1

            jm2 = jm1+jext
            call per_index(jm2,ny,jm2)
c            clock(1) = clock(1) + tclock() - clocktemp

c            clocktemp = tclock()
            ifint12 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &           ,im1+iext,jm2,km1,2,xinp2,yinp2,zinp2,icyl)
c            clock(2) = clock(2) + tclock() - clocktemp

            if(ifint12==1) then

c              clocktemp = tclock()
              ifph = fluidface(flag,nx,ny,nz,im1+iext,jm2,km1,2)
c              clock(3) = clock(3) + tclock() - clocktemp

              idir(2)=2
              if(ifph==1) then
                if(nphy==1) then
                  interp_points_pres=1
                  return
                else
                  interp_points_pres=0
                  return
                endif
              elseif(ifph==0) then
                interp_points_pres=0 
                return
              else
                interp_points_pres=-1
                return
              endif             
            else

c              clocktemp = tclock()
              dz = zu(km1)-q(3)
              if(rvec(3)>=0.0) then
                dz = zu(km1+1)-q(3)
              endif
              call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,3,dz,dely,1)
              km2 = km1+kext
c              clock(1) = clock(1) + tclock() - clocktemp

c              clocktemp = tclock()
              ifint13 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &             ,im1+iext,jm1,km2,3,xinp2,yinp2,zinp2,icyl)
c              clock(2) = clock(2) + tclock() - clocktemp

              if(ifint13==1) then

c                clocktemp = tclock()
                ifph = fluidface(flag,nx,ny,nz,im1+iext,jm1,km2,3)
c                clock(3) = clock(3) + tclock() - clocktemp

                idir(2)=3
                if(ifph==1) then
                  if(nphy==1) then
                    interp_points_pres=1
                    return
                  else
                    interp_points_pres=0
                    return
                  endif
                elseif(ifph==0) then
                  interp_points_pres=0
                  return
                else
                  interp_points_pres=-1
                  return
                endif
              endif

            endif

          endif

        endif

      else

c        clocktemp = tclock()
        call vec_ijkext_gridpoint(q,rvec,iext,jext,kext,icyl,2,theta,dely)
c        clock(1) = clock(1) + tclock() - clocktemp

c        clocktemp = tclock()
        ifint2 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &  ,im+iext,jm+jext,km+kext,2,xinp,yinp,zinp,icyl)
c        clock(2) = clock(2) + tclock() - clocktemp

        if(ifint2==1) then

c          clocktemp = tclock()
          ifph2 = fluidface(flag,nx,ny,nz,im+iext,jm+jext,km+kext,2)
c          clock(3) = clock(3) + tclock() - clocktemp

          if(ifph2==-1) then
            interp_points_pres = -1
            return
          else
            idir(1) = 2
            if(ifph2==1) then
              nphy=1
              interp_points_pres=1
            else
              interp_points_pres=0
            endif

            q(1) = xinp
            q(2) = yinp
            q(3) = zinp
          
            !Set new indices im1,jm1,km1 for new q
c            clocktemp = tclock()
            im1 = im+iext
            jm1 = jm+jext
            km1 = km+kext
            call per_index(jm1,ny,jm1)

            rm = sqrt(xu(im1,jm)**2. + yu(im1,jm)**2.)
            if(icyl==1 .AND. im1==1) rm=0.0
            call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,1,rm,dely,2) !0,-1
            
            im2 = im1+iext
            jm11 = jm1+jext
            call per_index(jm11,ny,jm11)
c            clock(1) = clock(1) + tclock() - clocktemp

c            clocktemp = tclock()
            ifint21 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &           ,im2,jm11,km1,1,xinp2,yinp2,zinp2,icyl)
c            clock(2) = clock(2) + tclock() - clocktemp

            if(ifint21==1) then

c              clocktemp = tclock()
              ifph = fluidface(flag,nx,ny,nz,im2,jm11,km1,1)
c              clock(3) = clock(3) + tclock() - clocktemp

              if(ifph==1) then
                if(nphy==1) then
                  interp_points_pres=1
                  return
                else
                  ncand = ncand+1
                  interp_points_pres=0
                  return
                endif
              elseif(ifph==0) then
                ncand = ncand+1
                interp_points_pres=0
                return
              else
                interp_points_pres=-1
                return
              endif
              idir(2)=1
            else

c              clocktemp = tclock()
              theta = dely*(real(jm11)-1.5)
              call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,2,rm,dely,2) !0,-1
              jm2 = jm1+jext
              call per_index(jm2,ny,jm2)
c              clock(1) = clock(1) + tclock() - clocktemp

c              clocktemp = tclock()
              ifint22 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &             ,im1,jm2,km1,2,xinp2,yinp2,zinp2,icyl)
c              clock(2) = clock(2) + tclock() - clocktemp

              if(ifint22==1) then

c                clocktemp = tclock()
                ifph = fluidface(flag,nx,ny,nz,im1,jm2,km1,2)
c                clock(3) = clock(3) + tclock() - clocktemp

                idir(2)=2
                if(ifph==1) then
                  if(nphy==1) then
                    interp_points_pres=1
                    return
                  else
                    ncand = ncand+1
                    interp_points_pres=0
                    return
                  endif
                elseif(ifph==0) then
                  ncand = ncand+1
                  interp_points_pres=0
                  return
                else
                  interp_points_pres=-1
                  return
                endif
              else

c              dz = zu(km1 + int(sign(1.0,rvec(3))))-zu(km1)
c                clocktemp = tclock()
                dz = zu(km1)-q(3)
                if(rvec(3)>=0.0) then
                  dz = zu(km1+1)-q(3)
                endif
                call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,3,theta,dely,2) !0,-1
       
                km2=km1+kext              
                jm11 = jm1+jext
                call per_index(jm11,ny,jm11)
c                clock(1) = clock(1) + tclock() - clocktemp

c                clocktemp = tclock()
                ifint23 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &               ,im1,jm11,km2,3,xinp2,yinp2,zinp2,icyl)
c                clock(2) = clock(2) + tclock() - clocktemp

                if(ifint23==1) then

c                  clocktemp = tclock()
                  ifph = fluidface(flag,nx,ny,nz,im1,jm11,km2,3)
c                  clock(3) = clock(3) + tclock() - clocktemp

                  idir(2)=3
                  if(ifph==1) then
                    if(nphy==1) then
                      interp_points_pres=1
                      return
                    else
                      ncand = ncand+1
                      interp_points_pres=0
                      return
                    endif
                  elseif(ifph==0) then
                    ncand = ncand+1
                    interp_points_pres=0
                    return
                  else
                    interp_points_pres=-1
                    return
                  endif
                endif

              endif

            endif

          endif

        else

c          clocktemp = tclock()
          dz = zu(km + int(sign(1.0,rvec(3))))-zu(km)
          call vec_ijkext_gridpoint(q,rvec,iext,jext,kext,icyl,3,dz,dely) !1,-1
c          clock(1) = clock(1) + tclock() - clocktemp

c          clocktemp = tclock()
          ifint3 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &          ,im+iext,jm+jext,km+kext,3,xinp,yinp,zinp,icyl)
c          clock(2) = clock(2) + tclock() - clocktemp

c          if(im==2 .AND. jm==2 .AND. km==173) then
c            write(6,*) '3. intep_points_pres:',im,jm,km,ifint3
c          endif


          if(ifint3==1) then

            idir(1) = 3

c            clocktemp = tclock()
            ifph3 = fluidface(flag,nx,ny,nz,im+iext,jm+jext,km+kext,3)
c            clock(3) = clock(3) + tclock() - clocktemp

            if(ifph3==-1) then
              interp_points_pres=-1
              return
            else
              if(ifph3==1) then
                nphy = 1
                interp_points_pres=1
              else
                interp_points_pres=0
              endif

              q(1) = xinp
              q(2) = yinp
              q(3) = zinp

c              clocktemp = tclock()
              im1 = im+iext
              jm1 = jm+jext
              km1 = km+kext
              call per_index(jm1,ny,jm1)

c            theta = dely*(real(jm1)-1.5)
              rm = sqrt(xu(im1,jm1)**2. + yu(im1,jm1)**2.)
              if(icyl==1 .AND. im1==1) rm=0.0
              call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,1,rm,dely,3) !0,-1
c              clock(1) = clock(1) + tclock() - clocktemp

              im2 = im1+iext
              km31 = km1+kext
        
              if(km31>=1.AND.km31<nz) then

c                clocktemp = tclock()
                ifint31 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &               ,im1+iext,jm1,km1+kext,1,xinp2,yinp2,zinp2,icyl)
c                clock(2) = clock(2) + tclock() - clocktemp

c                if(im==2 .AND. jm==2 .AND. km==173) then
c                  write(6,*) '1. intep_points_pres:',im,jm,km,ifint31
c                endif
              
                if(ifint31==1) then

c                  clocktemp = tclock()
                  ifph = fluidface(flag,nx,ny,nz,im1+iext,jm1,km1+kext,1)
c                  clock(3) = clock(3) + tclock() - clocktemp

                  idir(2)=1
                  if(ifph==1) then
                    if(nphy==1) then
                      interp_points_pres=1
                      return
                    else
                      ncand = ncand+1
                      interp_points_pres=0
                      return
                    endif
                  elseif(ifph==0) then
                    ncand = ncand+1
                    interp_points_pres=0
                    return
                  else
                    interp_points_pres=-1
                    return
                  endif
c                  goto 300
                else

c                  clocktemp = tclock()
                  theta = dely*(real(jm1)-1.5)
                  call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,2,theta,dely,3) !0,-1

                  jm2 =jm1+jext
                  call per_index(jm2,ny,jm2)
c                  clock(1) = clock(1) + tclock() - clocktemp

c                  clocktemp = tclock()
                  ifint32 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &                 ,im1,jm2,km1+kext,2,xinp2,yinp2,zinp2,icyl)
c                  clock(2) = clock(2) + tclock() - clocktemp

c                  if(im==2 .AND. jm==2 .AND. km==173) then
c                    write(6,*) '2. intep_points_pres:',im,jm,km,ifint32
c                  endif

                  if(ifint32==1) then

c                    clocktemp = tclock()
                    ifph = fluidface(flag,nx,ny,nz,im1,jm2,km1+kext,2)
c                    clock(3) = clock(3) + tclock() - clocktemp

                    idir(2)=2
                    if(ifph==1) then
                      if(nphy==1) then
                        interp_points_pres=1
                        return
                      else
                        ncand = ncand+1
                        interp_points_pres=0
                        return
                      endif
                    elseif(ifph==0) then
                      ncand = ncand+1
                      interp_points_pres=0
                      return
                    else
                      interp_points_pres=-1
                      return
                    endif
c                    goto 300
                  endif
                endif

              elseif(km31==0) then

                if(myrank>0) then
                  km31g = km31+(nz-2)*myrank
                  flagt(:,:,1) = flag2(:,:,1)
                  flagt(:,:,2) = flag(:,:,1)

c                  clocktemp = tclock()
                  call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,1,rm,dely,3) !0,-1
                  im2 = im1+iext
c                  clock(1) = clock(1) + tclock() - clocktemp

c                  clocktemp = tclock()
                  ifint31 = ray_face_int(q,rvec,xu,yu,zug,nx,ny,nzg
     &                 ,im2,jm1,km31g,1,xinp2,yinp2,zinp2,icyl)
c                  clock(2) = clock(2) + tclock() - clocktemp

                  if(ifint31==1) then

c                    clocktemp = tclock()
                    ifph = fluidface(flagt,nx,ny,2,im2,jm1,1,1)
c                    clock(3) = clock(3) + tclock() - clocktemp

                    idir(2)=1
                    if(ifph==1) then
                      if(nphy==1) then
                        interp_points_pres=1
                        return
                      else
                        ncand = ncand+1
                        interp_points_pres=0
                        return
                      endif
                    elseif(ifph==0) then
                      ncand = ncand+1
                      interp_points_pres=0
                      return
                    else
                      interp_points_pres=-1
                      return
                    endif
c                    goto 300
                  else

c                    clocktemp = tclock()
                    theta = dely*(real(jm1)-1.5)
                    call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,2,theta,dely,3) !0,-1
c                    clock(1) = clock(1) + tclock() - clocktemp

                    jm2=jm1+jext
c                    clocktemp = tclock()
                    ifint32 = ray_face_int(q,rvec,xu,yu,zug,nx,ny,nzg
     &                   ,im1,jm2,km31g,2,xinp2,yinp2,zinp2,icyl)
c                    clock(2) = clock(2) + tclock() - clocktemp

                    if(ifint32==1) then

c                      clocktemp = tclock()
                      ifph = fluidface(flagt,nx,ny,2,im1,jm2,1,2)
c                      clock(3) = clock(3) + tclock() - clocktemp

                      idir(2)=2
                      if(ifph==1) then
                        if(nphy==1) then
                          interp_points_pres=1
                          return
                        else
                          ncand = ncand+1
                          interp_points_pres=0
                          return
                        endif
                      elseif(ifph==0) then
                        ncand = ncand+1
                        interp_points_pres=0
                        return
                      else
                        interp_points_pres=-1
                        return
                      endif
c                      goto 300
                    endif
                  endif
                endif

              elseif(km31==nz) then

                if(myrank<mysize-1) then
                  km31g = km31+(nz-2)*myrank
                  flagt(:,:,1) = flag(:,:,nz)
                  flagt(:,:,2) = flag2(:,:,2)

c                  clocktemp = tclock()
                  call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,1,rm,dely,3) !0,-1
                  im2 = im1+iext
c                  clock(1) = clock(1) + tclock() - clocktemp

c                  clocktemp = tclock()
                  ifint31 = ray_face_int(q,rvec,xu,yu,zug,nx,ny,nzg
     &                 ,im2,jm1,km31g,1,xinp2,yinp2,zinp2,icyl)
c                  clock(2) = clock(2) + tclock() - clocktemp
                
                  if(ifint31==1) then

c                    clocktemp = tclock()
                    ifph = fluidface(flagt,nx,ny,2,im2,jm1,1,1)
c                    clock(3) = clock(3) + tclock() - clocktemp

                    idir(2)=1
                    if(ifph==1) then
                      if(nphy==1) then
                        interp_points_pres=1
                        return
                      else
                        ncand = ncand+1
                        interp_points_pres=0
                        return
                      endif
                    elseif(ifph==0) then
                      ncand = ncand+1
                      interp_points_pres=0
                      return
                    else
                      interp_points_pres=-1
                      return
                    endif
c                    goto 300
                  else

c                    clocktemp = tclock()
                    theta = dely*(real(jm1)-1.5)
                    call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,2,theta,dely,3) !0,-1
                    jm2 = jm1+jext
c                    clock(1) = clock(1) + tclock() - clocktemp

c                    clocktemp = tclock()
                    ifint32 = ray_face_int(q,rvec,xu,yu,zug,nx,ny,nzg
     &                   ,im1,jm2,km31g,2,xinp2,yinp2,zinp2,icyl)
c                    clock(2) = clock(2) + tclock() - clocktemp

                    if(ifint32==1) then

c                      clocktemp = tclock()
                      ifph =  fluidface(flagt,nx,ny,2,im1,jm2,1,2)
c                      clock(3) = clock(3) + tclock() - clocktemp

                      idir(2)=2
                      if(ifph==1) then
                        if(nphy==1) then
                          interp_points_pres=1
                          return
                        else
                          ncand = ncand+1
                          interp_points_pres=0
                          return
                        endif
                      elseif(ifph==0) then
                        ncand = ncand+1
                        interp_points_pres=0
                        return
                      else
                        interp_points_pres=-1
                        return
                      endif
c                      goto 300
                    endif
                  endif
                endif
 
              endif

c              clocktemp = tclock()
              km1g = km1+myrank*(nz-2)
              km2g = km1g+kext
              dz = zug(km2g)-zug(km1g)
              call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,3,dz,dely,3) !0,-1
              km2 = km1+kext
c              clock(1) = clock(1) + tclock() - clocktemp

              if(km2>=1.AND.km2<=nz) then

c                clocktemp = tclock()
                ifint33 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &              ,im1,jm1,km2,3,xinp2,yinp2,zinp2,icyl)
c                clock(2) = clock(2) + tclock() - clocktemp

c                if(im==2 .AND. jm==2 .AND. km==173) then
c                  write(6,*) '3. intep_points_pres:',im,jm,km,ifint33
c                endif

                if(ifint33==1) then
c                  clocktemp = tclock()
                  ifph = fluidface(flag,nx,ny,nz,im1,jm1,km2,3)
c                  clock(3) = clock(3) + tclock() - clocktemp

                  idir(2)=3
                  if(ifph==1) then
                    if(nphy==1) then
                      interp_points_pres=1
                      return
                    else
                      ncand = ncand+1
                      interp_points_pres=0
                      return
                    endif
                  elseif(ifph==0) then
                    ncand = ncand+1
                    interp_points_pres=0
                    return
                  else
                    interp_points_pres=-1
                    return
                  endif
                endif

              elseif(km2==0.AND.myrank>0) then
                km2g=km2+myrank*(nz-2)
c                clocktemp = tclock()
                ifint33 = ray_face_int(q,rvec,xu,yu,zug(km2g),nx,ny,1
     &               ,im1,jm1,1,3,xinp2,yinp2,zinp2,icyl)
c                clock(2) = clock(2) + tclock() - clocktemp

                if(ifint33==1) then

c                  clocktemp = tclock()
                  ifph = fluidface(flag2(1,1,1),nx,ny,1,im1,jm1,1,3)
c                  clock(3) = clock(3) + tclock() - clocktemp

                  idir(2)=3

                  if(ifph==1) then
                    if(nphy==1) then
                      interp_points_pres=1
                      return
                    else
                      ncand = ncand+1
                      interp_points_pres=0
                      return
                    endif
                  elseif(ifph==0) then
                    ncand = ncand+1
                    interp_points_pres=0
                    return
                  else
                    interp_points_pres=-1
                    return
                  endif
                endif
              elseif(km2==nz+1.AND.myrank<mysize-1) then
                km2g=km2+myrank*(nz-2)

c                clocktemp = tclock(2)
                ifint33 = ray_face_int(q,rvec,xu,yu,zug(km2g),nx,ny,1
     &               ,im1,jm1,1,3,xinp2,yinp2,zinp2,icyl)
c                clock(2) = clock(2) + tclock() - clocktemp

                if(ifint33==1) then

c                  clocktemp = tclock()
                  ifph = fluidface(flag2(1,1,2),nx,ny,1,im1,jm1,1,3)
c                  clock(3) = clock(3) + tclock() - clocktemp
                  idir(2)=3

                  if(ifph==1) then
                    if(nphy==1) then
                      interp_points_pres=1
                      return
                    else
                      ncand = ncand+1
                      interp_points_pres=0
                      return
                    endif
                  elseif(ifph==0) then
                    ncand = ncand+1
                    interp_points_pres=0
                    return
                  else
                    interp_points_pres=-1
                    return
                  endif
                endif
              endif
            endif
          endif
        endif
      endif

 300  continue


      return      
      end
C-----------------------------------------------------------------------



C---- function bndpts-----------------------------N. Beratlis-31 Dec. 2008--
C
      INTEGER function bndpts(flag,nx,ny,nz,ibd,nbd,iim,jim,kim,nim)
C
C     PURPOSE: Count and identify indices for boudnary points. Return number 
C     of pts.
C
C---------------------------------------------------------------------------
C
      IMPLICIT NONE
      include 'immersed.h'

      INTEGER nx,ny,nz,ibd,nbd,nim
      INTEGER flag(nx,ny,nz,nbd)
      INTEGER iim(nfcmax),jim(nfcmax),kim(nfcmax)

      INTEGER i,j,k,n

      n=0
c      write(6,*) jbmin,jbmax
c.....outer boundary points
      DO I=max(2,IBMIN(IBD)),IBMAX(IBD)
      DO J=JBMIN(IBD),JBMAX(IBD)
      DO K=KBMIN(IBD),KBMAX(IBD)
        if(flag(i,j,k,ibd)==-ibd) then
          n=n+1
          nim=nim+1
          iim(nim)=i
          jim(nim)=j
          kim(nim)=k
        endif
      ENDDO
      ENDDO
      ENDDO
      
      bndpts = n
c      
      RETURN
      END
C---------------------------------------------------------------------------




C---- function forcpts ---------------------------N. Beratlis-15 Aug. 2009--
C
      INTEGER function forcpts(flag,nx,ny,nz,ibd,nbd,iim,jim,kim,nim)
C
C     PURPOSE: Count and identify indices for boudnary points. Return number 
C     of pts.
C
C---------------------------------------------------------------------------
C
c flag: input; ibd interior points; 0 exterior points; -ibd interface points
c iim, jim, kim: output; coordinates of the interface points
c nim: number of the interface points
c
c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'
c
c.... Input/Output arrays
      INTEGER nx,ny,nz,ibd,nbd,nim
      INTEGER flag(nx,ny,nz,nbd)
      INTEGER iim(nfcmax),jim(nfcmax),kim(nfcmax)
c
c.... Local arrays
      INTEGER i,j,k,n,k1,k2

      n=0
      k1 = max(2,kbmin(ibd))
      k2 = min(nz-1,kbmax(ibd))

c.....outer boundary points
      DO K=k1,k2
      DO J=JBMIN(IBD),JBMAX(IBD)
      DO I=max(2,IBMIN(IBD)),IBMAX(IBD)
        if(flag(i,j,k,ibd)==-ibd) then
          n=n+1
          nim=nim+1
          iim(nim)=i
          jim(nim)=j
          kim(nim)=k
        endif
      ENDDO
      ENDDO
      ENDDO
      
      forcpts = n
c      
      RETURN
      END
C---------------------------------------------------------------------------


C---- function rdfn_forcpts------------------------N. Beratlis-31 Dec. 2008--
C
      INTEGER function rdfn_forcpts(flag,nx,ny,nz,ibd,iim,jim,kim
     &     ,xnim,ynim,znim,nxim,nyim,nzim,mrk,dir,ord,nq,lim,mim)
C
C     PURPOSE: Count and re-identify indices for boundary points. Return number 
C     of pts.
C
C---------------------------------------------------------------------------
C
      IMPLICIT NONE
      include 'immersed.h'

      INTEGER nx,ny,nz,ibd,nbd,nq,nn
      INTEGER lim,mim
      INTEGER flag(nx,ny,nz)
      INTEGER iim(nfcmax),jim(nfcmax),kim(nfcmax),mrk(nfcmax)
     &     ,dir(nfcmax)
      INTEGER ord(nq)
      REAL    nxim(nfcmax),nyim(nfcmax),nzim(nfcmax),
     &        xnim(nfcmax),ynim(nfcmax),znim(nfcmax)

      INTEGER i,j,k,n,im,jm,km,nim,ii,nvld,n1,n2,n3

      nvld = 0
      nim=0
      DO i=lim+1,lim+mim
        im = iim(i)
        jm = jim(i)
        km = kim(i)
        IF(flag(im,jm,km)==-ibd) THEN
          nim = nim+1
          ii = lim+nim
          iim(ii) = im
          jim(ii) = jm
          kim(ii) = km
          xnim(ii) = xnim(i)
          ynim(ii) = ynim(i)
          znim(ii) = znim(i)
          nxim(ii) = nxim(i)
          nyim(ii) = nyim(i)
          nzim(ii) = nzim(i)
          mrk (ii) = mrk (i)
          dir (ii) = dir (i)
          ord (ii) = ord (i)
c          innf(ii,:) = innf(i,:)
        ELSE
          nvld = nvld+1
        ENDIF

      ENDDO

      n1 = lim+nim+1
      n2 = nfcmax-nvld
      n3 = lim+mim+1
      iim(n1:n2)  = iim(n3:nfcmax)
      jim(n1:n2)  = jim(n3:nfcmax)
      kim(n1:n2)  = kim(n3:nfcmax)
      xnim(n1:n2) = xnim(n3:nfcmax)
      ynim(n1:n2) = ynim(n3:nfcmax)
      znim(n1:n2) = znim(n3:nfcmax)
      nxim(n1:n2) = nxim(n3:nfcmax)
      nyim(n1:n2) = nyim(n3:nfcmax)
      nzim(n1:n2) = nzim(n3:nfcmax)
      mrk(n1:n2)  = mrk(n3:nfcmax)
      dir(n1:n2)  = dir(n3:nfcmax)
      ord(lim+nim+1:nq-nvld) = ord(lim+mim+1:nq)
      mim = nim
      
      rdfn_forcpts = nim
c      
      RETURN
      END
C---------------------------------------------------------------------------



C---- function rdfn_forcpts_p ---------------------N. Beratlis-31 Dec. 2008--
C
      INTEGER function rdfn_forcpts_p(flag,nx,ny,nz,ibd,iim,jim,kim
     &     ,xnim,ynim,znim,nxim,nyim,nzim,mrk,dir,itr,ord,nq,lim,mim)
C
C     PURPOSE: Count and re-identify indices for boundary points. Return number 
C     of pts.
C
C---------------------------------------------------------------------------
C
      IMPLICIT NONE
      include 'immersed.h'

      INTEGER nx,ny,nz,ibd,nbd,nq
      INTEGER lim,mim
      INTEGER flag(nx,ny,nz)
      INTEGER iim(nfcmax),jim(nfcmax),kim(nfcmax),mrk(nfcmax)
     &     ,dir(nfcmax,2),itr(nfcmax)
      INTEGER ord(nq)
      REAL    nxim(nfcmax,2),nyim(nfcmax,2),nzim(nfcmax,2),
     &        xnim(nfcmax),ynim(nfcmax),znim(nfcmax)

      INTEGER i,j,k,n,im,jm,km,nim,ii,nvld,n1,n2,n3,n4,i2

      nvld = 0
      nim=0
      DO i=lim+1,lim+mim
        i2 = ord(i)
        im = iim(i)
        jm = jim(i)
        km = kim(i)
        IF(flag(im,jm,km)==-ibd) THEN
          nim = nim+1
          ii = lim+nim
          iim(ii) = im
          jim(ii) = jm
          kim(ii) = km
          xnim(ii) = xnim(i)
          ynim(ii) = ynim(i)
          znim(ii) = znim(i)
          nxim(ii,:) = nxim(i,:)
          nyim(ii,:) = nyim(i,:)
          nzim(ii,:) = nzim(i,:)
          mrk (ii)   = mrk (i)
          dir (ii,:) = dir (i,:)
          itr (ii) = itr (i)
          ord(ii) = ord(i)
c          where(ord>=i) ord=ord-1
c          where(ord>=i) ord=ord-1
c            ord(ii) = ord(i)-1
c          else
c            ord(ii) = ord(i)
c          endif
c          write(6,*) i,ii,ord(ii)
        ELSE
          nvld = nvld+1
        ENDIF

      ENDDO

      n1 = lim+nim+1
      n2 = nfcmax-nvld
      n3 = lim+mim+1
      n4 = nfcmax
      if(n2>n1) then
        iim(n1:n2)  = iim(n3:n4)
        jim(n1:n2)  = jim(n3:n4)
        kim(n1:n2)  = kim(n3:n4)
        xnim(n1:n2) = xnim(n3:n4)
        ynim(n1:n2) = ynim(n3:n4)
        znim(n1:n2) = znim(n3:n4)
        nxim(n1:n2,:) = nxim(n3:n4,:)
        nyim(n1:n2,:) = nyim(n3:n4,:)
        nzim(n1:n2,:) = nzim(n3:n4,:)
        mrk(n1:n2)  = mrk(n3:n4)
        dir(n1:n2,:)  = dir(n3:n4,:)
        itr(n1:n2)  = itr(n3:n4)
        ord(n1:nq-nvld) = ord(n3:nq)
      endif

      mim = nim
      
      rdfn_forcpts_p = nim
c      
      RETURN
      END
C---------------------------------------------------------------------------


C---- function rdfn_forcpts_p1 --------------------N. Beratlis-31 Dec. 2008--
C
      INTEGER function rdfn_forcpts_p1(flag,nx,ny,nz,ibd,iim,jim,kim
     &     ,xnim,ynim,znim,nxim,nyim,nzim,mrk,dir,itr,ord,nq,lim,mim)
C
C     PURPOSE: Count and re-identify indices for boundary points. Return number 
C     of pts.
C
C---------------------------------------------------------------------------
C
      IMPLICIT NONE
      include 'immersed.h'

      INTEGER nx,ny,nz,ibd,nbd,nq
      INTEGER lim,mim
      INTEGER flag(nx,ny,nz)
      INTEGER iim(nfcmax),jim(nfcmax),kim(nfcmax),mrk(nfcmax)
     &     ,dir(nfcmax,2),itr(nfcmax)
      INTEGER ord(nq)
      REAL    nxim(nfcmax,2),nyim(nfcmax,2),nzim(nfcmax,2),
     &        xnim(nfcmax),ynim(nfcmax),znim(nfcmax)
c
c.... Local arrays
      INTEGER i,j,k,n,im,jm,km,nim,ii,nvld,n1,n2,n3,n4,i2
      INTEGER iim2(nq),jim2(nq),kim2(nq),mrk2(nq),dir2(nq,2),itr2(nq)
      INTEGER ord2(nq)
      REAL    nxim2(nq,2),nyim2(nq,2),nzim2(nq,2),xnim2(nq),ynim2(nq),znim2(nq)

      iim2(1:nq) = iim(1:nq)
      jim2(1:nq) = jim(1:nq)
      kim2(1:nq) = kim(1:nq)
      mrk2(1:nq) = mrk(1:nq)
      dir2(1:nq,:) = dir(1:nq,:)
      nxim2(1:nq,:) = nxim(1:nq,:)
      nyim2(1:nq,:) = nyim(1:nq,:)
      nzim2(1:nq,:) = nzim(1:nq,:)
      xnim2(1:nq) = xnim(1:nq)
      ynim2(1:nq) = ynim(1:nq)
      znim2(1:nq) = znim(1:nq)
      itr2(1:nq) = itr(1:nq)
      ord2(1:nq) = ord(1:nq)
      
      nvld = 0
      nim=0
      DO i=lim+1,lim+mim
        i2 = ord2(i)
        im = iim2(i2)
        jm = jim2(i2)
        km = kim2(i2)
        IF(flag(im,jm,km)==-ibd) THEN
          nim = nim+1
          ii = lim+nim
          iim(i2) = im
          jim(i2) = jm
          kim(i2) = km
          xnim(i2) = xnim2(i2)
          ynim(i2) = ynim2(i2)
          znim(i2) = znim2(i2)
          nxim(i2,:) = nxim2(i2,:)
          nyim(i2,:) = nyim2(i2,:)
          nzim(i2,:) = nzim2(i2,:)
          mrk (i2)   = mrk2 (i2)
          dir (i2,:) = dir2 (i2,:)
          itr (i2) = itr2 (i2)
          ord (i2) = i2
c          ord(ii) = ord(i)
c          where(ord>=i) ord=ord-1
c          where(ord>=i) ord=ord-1
c            ord(ii) = ord(i)-1
c          else
c            ord(ii) = ord(i)
c          endif
c          write(6,*) i,ii,ord(ii)
        ELSE
          nvld = nvld+1
        ENDIF

      ENDDO

      n1 = lim+nim+1
      n2 = nfcmax-nvld
      n3 = lim+mim+1
      n4 = nfcmax
      if(n2>n1) then
        iim(n1:n2)  = iim(n3:n4)
        jim(n1:n2)  = jim(n3:n4)
        kim(n1:n2)  = kim(n3:n4)
        xnim(n1:n2) = xnim(n3:n4)
        ynim(n1:n2) = ynim(n3:n4)
        znim(n1:n2) = znim(n3:n4)
        nxim(n1:n2,:) = nxim(n3:n4,:)
        nyim(n1:n2,:) = nyim(n3:n4,:)
        nzim(n1:n2,:) = nzim(n3:n4,:)
        mrk(n1:n2)  = mrk(n3:n4)
        dir(n1:n2,:)  = dir(n3:n4,:)
        itr(n1:n2)  = itr(n3:n4)
        ord(n1:nq-nvld) = ord(n3:nq)
      endif

      mim = nim
      
      rdfn_forcpts_p1 = nim
c      
      RETURN
      END
C---------------------------------------------------------------------------



C---- function rdfn_bndpts------------------------N. Beratlis-31 Dec. 2008--
C
      INTEGER function rdfn_bndpts(flag,nx,ny,nz,ibd,nbd,iim,jim,kim
     &     ,xnim,ynim,znim,nxim,nyim,nzim,mrk,dir,lim,mim)
c         mim(ibd,2) = rdfn_bndpts(flagi,nx,ny,nz,ibd,nbd,iim,jim,kim
c     &     ,xnim,ynim,znim,nxim,nyim,nzim,mrk,dir,lim(ibd,io),mim(ibd,io))
c        nim = mim(ibd,1)+mim(ibd,2)
C
C     PURPOSE: Count and re-identify indices for boundary points. Return number 
C     of pts.
C
C---------------------------------------------------------------------------
C
      IMPLICIT NONE
      include 'immersed.h'

      INTEGER nx,ny,nz,ibd,nbd
      INTEGER lim,mim
      INTEGER flag(nx,ny,nz,nbd)
      INTEGER iim(nfcmax),jim(nfcmax),kim(nfcmax),mrk(nfcmax)
     &     ,dir(nfcmax)
      REAL    nxim(nfcmax),nyim(nfcmax),znim(nfcmax),
     &        xnim(nfcmax),ynim(nfcmax),nzim(nfcmax)

      INTEGER i,j,k,n,im,jm,km,nim,ii

      nim=0
      DO i=lim+1,lim+mim
        im = iim(i)
        jm = jim(i)
        km = kim(i)
        IF(flag(im,jm,km,ibd)==-ibd) THEN
          nim = nim+1
          ii = lim+nim
          iim(ii) = im
          jim(ii) = jm
          kim(ii) = km
          xnim(ii) = xnim(i)
          ynim(ii) = ynim(i)
          znim(ii) = znim(i)
          nxim(ii) = nxim(i)
          nyim(ii) = nyim(i)
          nzim(ii) = nzim(i)
          mrk (ii) = mrk (i)
          dir (ii) = dir (i)
        ENDIF

      ENDDO
c      mim = nim-lim
c      mim = lim+nim
      mim = nim
      
      rdfn_bndpts = nim
c      
      RETURN
      END
C---------------------------------------------------------------------------



C---- function rdfn_bndpts_p----------------------N. Beratlis-31 Dec. 2008--
C
      INTEGER function rdfn_bndpts_p(flag,nx,ny,nz,ibd,nbd,iim,jim,kim
     &        ,xnim,ynim,znim,nxim,nyim,nzim,mrk,dir,itr,lim,mim)
C
C     PURPOSE: Count and re-identify indices for boundary points. Return number 
C     of pts.
C
C---------------------------------------------------------------------------
C
      IMPLICIT NONE
      include 'immersed.h'

      INTEGER nx,ny,nz,ibd,nbd
      INTEGER lim,mim
      INTEGER flag(nx,ny,nz,nbd)
      INTEGER iim(nfcmax),jim(nfcmax),kim(nfcmax),mrk(nfcmax)
      INTEGER dir(nfcmax,2),itr(nfcmax)
      REAL    xnim(nfcmax),ynim(nfcmax),znim(nfcmax),
     &        nxim(nfcmax,2),nyim(nfcmax,2),nzim(nfcmax,2)

      INTEGER i,j,k,n,im,jm,km,nim,ii

      nim=0
      DO i=lim+1,lim+mim
        im = iim(i)
        jm = jim(i)
        km = kim(i)
        IF(flag(im,jm,km,ibd)==-ibd) THEN
          nim = nim+1
          ii = lim+nim
          iim(ii) = im
          jim(ii) = jm
          kim(ii) = km
          xnim(ii) = xnim(i)
          ynim(ii) = ynim(i)
          znim(ii) = znim(i)
          nxim(ii,:) = nxim(i,:)
          nyim(ii,:) = nyim(i,:)
          nzim(ii,:) = nzim(i,:)
          mrk (ii) = mrk (i)
          dir(ii,:)= dir(i,:)
          itr(ii) = itr(i)
        ENDIF

      ENDDO
      mim = nim

      
      rdfn_bndpts_p = nim
c      
      RETURN
      END
C---------------------------------------------------------------------------



C---- subroutine univector----------------N. Beratlis-Jun 11 2009---
C
C     PURPOSE: Calculate unit vector from (x1,y1,z1) to (x2,y2,z2)
C
C-------------------------------------------------------------------
      subroutine unitvector(x1,x2,y1,y2,z1,z2,vec)

      REAL x1,x2,y1,y2,z1,z2
      REAL vec(3)

      REAL mag

      vec(1) = x2-x1
      vec(2) = y2-y1
      vec(3) = z2-z1
      mag = sqrt(vec(1)**2. + vec(2)**2. + vec(3)**2.)
      if(mag>0.0) vec = vec/mag

      return
      end
C-------------------------------------------------------------------



C-----SUBROUTINE PRESGRAD_FLAGP-------------------N. Beratlis-18 Dec. 2008--
C
      SUBROUTINE PRESGRAD_FLAGP(dpdl,mrkp,xnp,ynp,znp,nxp,nyp,nzp,ip,jp,kp
     &     ,xc_car,yc_car,zc,nx,ny,nz,unvect,fp,nfacet
     &     ,limp,mimp,tlevel,dt,ibd)
C
C     PURPOSE: Calculate pressure gradient on surface of immersed body
C
C---------------------------------------------------------------------------
      
      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'mpif.h'

      integer limp,mimp,nfacet,ibd,nx,ny,nz
      real    tlevel,dt
      integer fp(nfcmax),mrkp(nfcmax)
      real    unvect(3,nfacet)
      real    dpdl(nfcmax),xnp(nfcmax),ynp(nfcmax),znp(nfcmax)
      real    nxp(nfcmax),nyp(nfcmax),nzp(nfcmax)
      integer ip(nfcmax),jp(nfcmax),kp(nfcmax)
      real    xc_car(nx,ny),yc_car(nx,ny),zc(nz)
c
c.... Local arrays
      integer im,i,j,k,i1,j1,k1,itr
      real    xn,yn,zn
      real    a,b,c,d
      real    r(3)
c
c.... Functions
      real    dpdn,dpdx1,dpdx2,dpdx3,vecmag

      dpdl(limp+1:limp+mimp) = 0.0     !!!!!! instead of dpdl = 0.0

      do im=limp+1,limp+mimp
        i = ip(im)
        j = jp(im)
        k = kp(im)

        if(mrkp(im)==1) then
c extension point from the NORM_INTR
          itr = fp(im)
c the function DPDN evaluates the pressure gradient along the normal to
c the body
          dpdl(im) = dpdn(xnp(im),ynp(im),znp(im),unvect(1,itr),unvect(2,itr),unvect(3,itr),tlevel,dt,ibd)
        elseif(mrkp(im)==2) then
c extension point from the GRID_INTR
          i1 = i+int(nxp(im))
          j1 = j+int(nyp(im))
          k1 = k+int(nzp(im))
          call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
          r(1) = xc_car(i1,j1)-xc_car(i,j)
          r(2) = yc_car(i1,j1)-yc_car(i,j)
          r(3) = zc(k1)-zc(k)
          d = vecmag(r,3)
          r = r/d
c R is the unit vector between the interface point and the fluid point
c DPDX1, DPDX2 and DPDX3 are the pressure gradients along the directions
c X, Y and Z
c DPDL is the total pressure gradient along the direction defined by the
c unit vector R
          dpdl(im) = dpdx1(xnp(im),ynp(im),znp(im),tlevel,dt,ibd)*r(1)
     &         + dpdx2(xnp(im),ynp(im),znp(im),tlevel,dt,ibd)*r(2)
     &         + dpdx3(xnp(im),ynp(im),znp(im),tlevel,dt,ibd)*r(3)
        else
c extension point from the CLOSEST_INTR
          r(1) = nxp(im)-xnp(im)
          r(2) = nyp(im)-ynp(im)
          r(3) = nzp(im)-znp(im)
          r = r/vecmag(r,3)
c R is the unit vector between the intersection (closest) point and the extension point
c this criterion is similar to the one used above for MRKP=2 (extension point from GRID_INTR)
          dpdl(im) = dpdx1(xnp(im),ynp(im),znp(im),tlevel,dt,ibd)*r(1)
     &         + dpdx2(xnp(im),ynp(im),znp(im),tlevel,dt,ibd)*r(2)
     &         + dpdx3(xnp(im),ynp(im),znp(im),tlevel,dt,ibd)*r(3)
        endif

      enddo

      return

      end
C---------------------------------------------------------------------------


C---- subroutine determine_flagpopi---------------N. Beratlis-31 Jan. 2008--
C
      SUBROUTINE determine_flagpopi(flagpo,flagpi,nx,ny,nz,ibd,nbd,i,j,k,ifld,ibnd,ibdy)
C
C     PURPOSE: Determine if pressure point i,j,k should be fluid, forcing
C     or body point.
C---------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER nx,ny,nz,ibd,nbd,i,j,k,ifld,ibnd,ibdy
      INTEGER flagpo(nx,ny,nz,nbd),flagpi(nx,ny,nz,nbd)
c
c.....outer pressure boundary points
      IF(ibdy==0.AND.ibnd/=0.AND.ifld/=0) THEN
        flagpo(i,j,k,ibd)=-ibd
      ENDIF
c
      IF(ifld==0.AND.ibnd/=0.AND.ibdy/=0) THEN
c.....inner pressure boundary points
        flagpi(i,j,k,ibd)=-ibd
c.....outer pressure body points
        flagpo(i,j,k,ibd)=ibd
      ENDIF
c
      IF(ifld==0.AND.ibnd==0.AND.ibdy/=0) THEN
c.....inner pressure body points
        flagpi(i,j,k,ibd)=ibd
c.....outer pressure body points
        flagpo(i,j,k,ibd)=ibd
      ENDIF

      RETURN
      END
C---------------------------------------------------------------------------



C---- subtourine sort_cand_stencil----------------N. Beratlis-31 Dec. 2008--
C
      SUBROUTINE sort_cand_stencil(DS,INTRS,IORDER,ORD)
C
C     PURPOSE: Sort cand stencil point by variable DS. If ORD=1, sorting
C     is in ascending order, if ORD=2 sorting is in descending order.
C
C---------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER INTRS,ORD
      REAL DS(INTRS)

      INTEGER NITRS,IMIN,IMAX
      INTEGER IORDER(INTRS)
      REAL    DS1(INTRS)

      DS1=DS
      IF(ORD==1) THEN
        DO NITRS=1,INTRS
          IMIN = MINLOC(DS1(1:INTRS),1)
          IORDER(NITRS) = IMIN
          DS1(IMIN) = 10000.
        ENDDO
      ELSE
        DO NITRS=1,INTRS
          IMAX = MAXLOC(DS1(1:INTRS),1)
          IORDER(NITRS) = IMAX
          DS1(IMAX) = -10000.
        ENDDO
      ENDIF

      RETURN
      
      END
C---------------------------------------------------------------------------



C---- function stencil_normal_unitdot-------------N. Beratlis-31 Dec. 2008--
C
      REAL function stencil_normal_unitdot(x,y,z,xn,yn,zn,nx,ny,nz)
C
C     PURPOSE: Calculate dot product between unit vector (x-xn,y-yn,z-zn) 
C     and (nx,ny,nz).
C
C---------------------------------------------------------------------------
      IMPLICIT NONE

      REAL    x,y,z,xn,yn,zn,nx,ny,nz

      REAL    a(3),b(3),vecmag,dotprod

      a(1) = x-xn
      a(2) = y-yn
      a(3) = z-zn
      a = a/vecmag(a,3)

      b(1) = nx
      b(2) = ny
      b(3) = nz
      b = b/vecmag(b,3)

      stencil_normal_unitdot = dotprod(a,b)

      RETURN
      END
C---------------------------------------------------------------------------



C---- function physicalface--------------------N. Beratlis-26 Dec. 2008-
C
      INTEGER function physicalface(flago,nx,ny,nz,i,j,k,nbd,dir)
C
C     PURPOSE: Determine if face of cell is physical.
C     Returns: 0 - if it conatins body point.
C
C-----------------------------------------------------------------------
      INTEGER i,j,k,nx,ny,nz,dir
      INTEGER flago(nx,ny,nz,nbd)

      INTEGER i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4,ibd
      INTEGER itst(4,nbd)

c      write(6,*) 'physicalface:',nx,ny,nz,i,j,k,ibd,nbd,dir


      i1=i
      j1=j
      k1=k      
      
      IF(dir==1) THEN
        i2 = i
        j2 = j+1
        k2 = k
        i3 = i
        j3 = j
        k3 = k+1
        i4 = i
        j4 = j+1
        k4 = k+1
      ELSEIF(dir==2) THEN
        i2 = i+1
        j2 = j
        k2 = k
        i3 = i
        j3 = j
        k3 = k+1
        i4 = i+1
        j4 = j
        k4 = k+1
      ELSEIF(dir==3) THEN
        i2 = i+1
        j2 = j
        k2 = k
        i3 = i
        j3 = j+1
        k3 = k
        i4 = i+1
        j4 = j+1
        k4 = k
      ENDIF

      call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
      call index_bnd(i2,j2,k2,i2,j2,k2,nx,ny,nz)
      call index_bnd(i3,j3,k3,i3,j3,k3,nx,ny,nz)
      call index_bnd(i4,j4,k4,i4,j4,k4,nx,ny,nz)

      itst = 0
      do ibd=1,nbd
        itst(1,ibd) = (flago(i1,j1,k1,ibd)-ibd)
        itst(2,ibd) = (flago(i2,j2,k2,ibd)-ibd)
        itst(3,ibd) = (flago(i3,j3,k3,ibd)-ibd)
        itst(4,ibd) = (flago(i4,j4,k4,ibd)-ibd)
      enddo

      physicalface = product(itst)

c      write(6,*) 'physicalface=',physicalface

      RETURN

      END
C-----------------------------------------------------------------------



C---- function fluidface-----------------------N. Beratlis-26 Dec. 2008-
C
      INTEGER function fluidface(flag,nx,ny,nz,i,j,k,dir)
C
C     PURPOSE: Determine if face of cell is physical.
C     RETURNS: 1 if all points are fluid points (0) - SUCCESS
C              0 if face contains any forcing points (negative flag)
C             -1 if face contains any body points (postive flag)
C
C-----------------------------------------------------------------------
      INTEGER i,j,k,nx,ny,nz,dir
      INTEGER flag(nx,ny,nz)

      INTEGER ii,itst(4)

      if(i<1.OR.i>nx.OR.j<1.OR.j>ny) then
        fluidface=0
        return
      endif

      IF(dir==1) THEN
        itst(1) = flag(i,j  ,k  )
        itst(2) = flag(i,j+1,k  )
        itst(3) = flag(i,j  ,k+1)
        itst(4) = flag(i,j+1,k+1)
      ELSEIF(dir==2) THEN
        itst(1) = flag(i  ,j,k  )
        itst(2) = flag(i+1,j,k  )
        itst(3) = flag(i  ,j,k+1)
        itst(4) = flag(i+1,j,k+1)
      ELSEIF(dir==3) THEN
        itst(1) = flag(i  ,j  ,k)
        itst(2) = flag(i+1,j  ,k)
        itst(3) = flag(i  ,j+1,k)
        itst(4) = flag(i+1,j+1,k)
      ENDIF


      fluidface=1

      if(any(itst/=0)) then
        fluidface=0
        if(any(itst>0)) then
          fluidface=-1
        endif
      endif

c      DO ii=1,4        
c        IF(itst(ii)/=0) fluidface=0      
c      ENDDO

      RETURN

      END
C-----------------------------------------------------------------------



C-----function physicalcell-----------------------N. Beratlis-21 Dec. 2008--
C
      INTEGER function physicalcell(flagpo,flagpi,nx,ny,nz,ibd,nbd,i,j,k)
C
C     PURPOSE: Determine if cell i,j,k is surrounded by non-body nodes
C
C---------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER i,j,k,ibd,nbd,nx,ny,nz
      INTEGER flagpo(nx,ny,nz,nbd),flagpi(nx,ny,nz,nbd)

      INTEGER itst(8)

      itst(1) = (flagpo(i  ,j  ,k  ,ibd)-ibd)*(flagpi(i  ,j  ,k  ,ibd)-ibd)
      itst(2) = (flagpo(i+1,j  ,k  ,ibd)-ibd)*(flagpi(i+1,j  ,k  ,ibd)-ibd)
      itst(3) = (flagpo(i  ,j+1,k  ,ibd)-ibd)*(flagpi(i  ,j+1,k  ,ibd)-ibd)
      itst(4) = (flagpo(i  ,j  ,k+1,ibd)-ibd)*(flagpi(i  ,j  ,k+1,ibd)-ibd)
      itst(5) = (flagpo(i  ,j+1,k+1,ibd)-ibd)*(flagpi(i  ,j+1,k+1,ibd)-ibd)
      itst(6) = (flagpo(i+1,j+1,k  ,ibd)-ibd)*(flagpi(i+1,j+1,k  ,ibd)-ibd)
      itst(7) = (flagpo(i+1,j  ,k+1,ibd)-ibd)*(flagpi(i+1,j  ,k+1,ibd)-ibd)
      itst(8) = (flagpo(i+1,j+1,k+1,ibd)-ibd)*(flagpi(i+1,j+1,k+1,ibd)-ibd)
      
      physicalcell = product(itst)

c      IF(I==NX.OR.J==NY.OR.K==NZ) WRITE(6,*) i,j,k,physicalcell

      RETURN
      END
C---------------------------------------------------------------------------

C-----function physicalnodes----------------------N. Beratlis-21 Dec. 2008--
C
      INTEGER function physicalnodes(flag,nx,ny,nz,i,j,k)
C
C     PURPOSE: Determine if cell i,j,k is surrounded by non-body nodes
c     RETURN: 0 if all nodes are fluids
c            -1 if no body fluids but at least one forcing node
C             1 if at least one body node
C
C---------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER i,j,k,ibd,nbd,nx,ny,nz
      INTEGER flag(nx,ny,nz)

      INTEGER itst(8)

      itst(1) = flag( i , j , k )
      itst(2) = flag(i+1, j , k )
      itst(3) = flag( i ,j+1, k )
      itst(4) = flag( i ,j+1,k+1)
      itst(5) = flag(i+1,j+1, k )
      itst(6) = flag( i ,j+1,k+1)
      itst(7) = flag(i+1, j ,k+1)
      itst(8) = flag(i+1,j+1,k+1)

      physicalnodes = 0
      if(any(itst>0)) then
        physicalnodes = 1
      elseif(any(itst<0)) then
        physicalnodes = -1
      endif

      RETURN
      END
C---------------------------------------------------------------------------



C-----function point_fld -------------------------N. Beratlis-21 Dec. 2008--
C
      INTEGER function point_fld(xp,yp,zp,flag,xu,yu,zu,nx,ny,nz,icyl)
C
C     PURPOSE: Determine if cell i,j,k is surrounded by non-body nodes
c     RETURN: 0 if all nodes are fluids
c            -1 if no body fluids but at least one forcing node
C             1 if at least one body node
C
C---------------------------------------------------------------------------
      IMPLICIT NONE
c
c.... Input/Output arrays
      INTEGER nx,ny,nz,icyl
      REAL    xp,yp,zp
      REAL    xu(nx),yu(ny),zu(nz)
      INTEGER flag(nx,ny,nz)
c
c.... Local arrays
      INTEGER i,j,k
      INTEGER itst(8)
c
c.... Functions
      REAL    anglerad

      call ijk_xyz(xp,yp,zp,xu,yu,zu,nx,ny,nz,i,j,k,icyl)

      itst(1) = flag( i , j , k )
      itst(2) = flag( i ,j+1 ,k )
      itst(3) = flag( i , j ,k+1)
      itst(4) = flag( i ,j+1,k+1)
      itst(5) = flag(i+1, j , k )
      itst(6) = flag(i+1,j+1 ,k )
      itst(7) = flag(i+1, j ,k+1)
      itst(8) = flag(i+1,j+1,k+1)

      point_fld = 0
      if(any(itst>0)) then
        point_fld = 1
      elseif(any(itst<0)) then
        point_fld = -1
      endif

      RETURN
      END
C---------------------------------------------------------------------------


C---- function physicalface1-----------------------N. Beratlis-26 Dec. 2008-
C
      INTEGER function physicalface1(flago,flagi,nx,ny,nz,i,j,k,nbd,dir)
C
C     PURPOSE: Determine if face of cell is physical. Assumes only one
C     immersed body in the domain
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER i,j,k,nx,ny,nz,nbd,dir
      INTEGER flagi(nx,ny,nz,nbd),flago(nx,ny,nz,nbd)

      INTEGER i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4,ibd
      INTEGER itst(4)

      i1=i
      j1=j
      k1=k
      
      IF(dir==1) THEN
        i2 = i
        j2 = j+1
        k2 = k
        i3 = i
        j3 = j
        k3 = k+1
        i4 = i
        j4 = j+1
        k4 = k+1
      ELSEIF(dir==2) THEN
        i2 = i+1
        j2 = j
        k2 = k
        i3 = i
        j3 = j
        k3 = k+1
        i4 = i+1
        j4 = j
        k4 = k+1
      ELSEIF(dir==3) THEN
        i2 = i+1
        j2 = j
        k2 = k
        i3 = i
        j3 = j+1
        k3 = k
        i4 = i+1
        j4 = j+1
        k4 = k
      ENDIF

      call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
      call index_bnd(i2,j2,k2,i2,j2,k2,nx,ny,nz)
      call index_bnd(i3,j3,k3,i3,j3,k3,nx,ny,nz)
      call index_bnd(i4,j4,k4,i4,j4,k4,nx,ny,nz)

      itst = 0
      do ibd=1,nbd
        itst(1) = itst(1) + (flago(i1,j1,k1,ibd)-ibd)
        itst(2) = itst(2) + (flago(i2,j2,k2,ibd)-ibd)
        itst(3) = itst(3) + (flago(i3,j3,k3,ibd)-ibd)
        itst(4) = itst(4) + (flago(i4,j4,k4,ibd)-ibd)
      enddo

      physicalface1 = product(itst)

      RETURN

      END
C-----------------------------------------------------------------------



C---- subroutine index_bnd --------------------N. Beratlis-17 Jan. 2009-
C
C     PURPOSE: Rearrange indices that are not in the interior of the domain
C
C-----------------------------------------------------------------------
      subroutine index_bnd(i1,j1,k1,i2,j2,k2,nx,ny,nz)

      include 'common.h'
      
      integer i1,j1,k1,i2,j2,k2,nx,ny,nz

c      write(6,*) i1,j1,k1,i2,j2,k2,nx,ny,nz
      i2=i1
      j2=j1
      k2=k1

      if(j1<=1) then
        j2 = ny-2+j1
      elseif(j1>=ny) then
        j2 = j1-ny+2
      endif

c      write(6,*) i2,j2,k2

      if(i1<1 .AND. icyl==1) then
        i2 = 3-i1
        j2 = jsym(j1)
      endif
c      write(6,*) i2,j2,k2

      return
      end
C-----------------------------------------------------------------------




C---- subroutine fillextra_cell_z--------------N. Beratlis-7 Jan. 2009---
C
C     PURPOSE: Communicate extra k=2,nz+1 cell information among procs.
C
C-------------------------------------------------------------------------
      subroutine fillextra_cell_z(flag,flag2,nx,ny,nz,myleft,myright
     &           ,mpi_comm_eddy)

      IMPLICIT NONE
      include 'mpif.h'

      INTEGER nx,ny,nz,myleft,myright,mpi_comm_eddy
      INTEGER flag(nx,ny,nz),flag2(nx,ny,2)
      
      INTEGER ierr,status(mpi_status_size)

      CALL MPI_SENDRECV( flag(1,1,3),nx*ny,MPI_INTEGER,MYRIGHT,10
     &     ,            flag2(1,1,1),nx*ny,MPI_INTEGER,MYLEFT ,10
     &     ,            MPI_COMM_EDDY,STATUS,IERR)
      
      CALL MPI_SENDRECV(flag(1,1,nz-2),nx*ny,MPI_INTEGER,MYLEFT ,11
     &     ,            flag2(1,1,2  ),nx*ny,MPI_INTEGER,MYRIGHT,11
     &     ,            MPI_COMM_EDDY,STATUS,IERR)


      RETURN

      END
C----------------------------------------------------------------------



C-----SUBROUTINE PRESGRAD_FLAGP1--------------N. Beratlis-11 Jan. 2009-
C
C     PURPOSE: Calculate pressure gradient on surface of immersed body
C     by extrapolating using 2 fluid pressure points.
C
C----------------------------------------------------------------------
      SUBROUTINE PRESGRAD_FLAGP1(P,NX,NY,NZ,IP,JP,KP,PIM,NXP,NYP,NZP
     &     ,DIR,MRKP,DPDSS,XC,YC,ZC,LIMP,MIMP)
      
      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'mpif.h'
      
      INTEGER nx,ny,nz
      INTEGER ip(nfcmax),jp(nfcmax),kp(nfcmax)
      INTEGER dir(nfcmax,2),mrkp(nfcmax)
      INTEGER limp,mimp
      REAL    XC(NX),YC(NY),ZC(NZ)
      REAL    P(NX,NY,NZ)
      REAL    NXP(nfcmax,2),NYP(nfcmax,2),NZP(nfcmax,2)
      REAL    DPDSS(nfcmax),PIM(nfcmax)
c
c.... Local arrays
      INTEGER I,J,K,IM,II,I1,J1,K1,I2,J2,K2,IL,IR,IL2,IR2,NVS,NIM
      INTEGER imn,img,idg,nnorm,ngrid,ndiag
      INTEGER STATUS(MPI_STATUS_SIZE)
      REAL    DS,xcar,ycar,pint1,pint2,x1car,y1car
      INTEGER IML(nfcmax),IMR(nfcmax),LI(3,nfcmax),RI(3,nfcmax)
      INTEGER dirl(nfcmax),dirr(nfcmax)
      REAL    XL(nfcmax),YL(nfcmax),ZL(nfcmax)
     &       ,XR(nfcmax),YR(nfcmax),ZR(nfcmax)
c      REAL    p1l(nfcmax),p1r(nfcmax)
      INTEGER NL_REQ,NR_REQ,NL_SUP,NR_SUP
      REAL, DIMENSION(:), ALLOCATABLE :: PL_REQ,PL_SUP,PR_REQ,PR_SUP
     &                                  ,ZL_REQ,ZL_SUP,ZR_REQ,ZR_SUP
     &                                  ,XL_REQ,YL_REQ,XR_REQ,YR_REQ
     &                                  ,XL_SUP,YL_SUP,XR_SUP,YR_SUP
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IL_REQ,IL_SUP,IR_REQ,IR_SUP
      INTEGER, DIMENSION(:), ALLOCATABLE :: DL_REQ,DL_SUP,DR_REQ,DR_SUP
      INTEGER, DIMENSION(:), ALLOCATABLE :: INRM,IGRD,IDIA
c      REAL, DIMENSION(:), ALLOCATABLE :: p1l,p1r

c      INTEGER IML1(3,nfcmax),IMR1(3,nfcmax),LI1(3,nfcmax),RI1(3,nfcmax)
c
c.... Functions
      REAL    calc_dpdn,interp_cellface,calc_dpdn1,interp_cellface1,anglerad

      nim = mimp
      ALLOCATE(INRM(nim),IGRD(nim),IDIA(nim))

      IL=0
      IR=0

      INRM=0
      IGRD=0
      i=0
      j=0
      k=0

      DO IM=LIMP+1,LIMP+MIMP
        if(mrkp(im)==1) then
          i=i+1
          inrm(i)=im
        elseif(mrkp(im)==2) then
          j=j+1
          igrd(j)=im
        elseif(mrkp(im)==3) then
          k=k+1
          idia(k)=im
        endif
      ENDDO
      ngrid = j
      nnorm = i
      ndiag = k



      DO imn=1,nnorm

        im = inrm(imn)
        i = ip(im)
        j = jp(im)
        k = kp(im)

c        if(i==70 .AND. j==4 .AND. k+myrank*(nz-2)==202) then
c          write(6,*) 'correctpres',i,j,k,im
c        endif
 
        pim(im) = interp_cellface(nxp(im,1),nyp(im,1),nzp(im,1)
     &       ,xc,yc,zc,p,nx,ny,nz,icyl,dir(im,1))

        if(abs(dir(im,1))==3) then
          if(nzp(im,2)<zc(1)) then
            il = il+1
            iml(il) = im
            xl(il) = nxp(im,2)
            yl(il) = nyp(im,2)
            zl(il) = nzp(im,2)
            dirl(il) = dir(im,2)

          elseif(nzp(im,2)>zc(nz)) then
            ir = ir+1
            imr(ir) = im
            xr(ir) = nxp(im,2)
            yr(ir) = nyp(im,2)
            zr(ir) = nzp(im,2)
            dirr(ir) = dir(im,2)
          else
            pint2 = interp_cellface(nxp(im,2),nyp(im,2),nzp(im,2)
     &            ,xc,yc,zc,p,nx,ny,nz,icyl,dir(im,2))
            dpdss(im) = calc_dpdn(nxp(im,1),nyp(im,1),nzp(im,1)
     &           ,nxp(im,2),nyp(im,2),nzp(im,2),pim(im),pint2,0)
          endif
        else
          pint2 = interp_cellface(nxp(im,2),nyp(im,2),nzp(im,2)
     &          ,xc,yc,zc,p,nx,ny,nz,icyl,dir(im,2))
          dpdss(im) = calc_dpdn(nxp(im,1),nyp(im,1),nzp(im,1)
     &         ,nxp(im,2),nyp(im,2),nzp(im,2),pim(im),pint2,0)
        endif

      ENDDO


      NL_REQ = IL !Number of points required from left processor
      NR_REQ = IR !Number of points required from right processor


      NVS=20

c....Check if info from adjacent domain is needed
      IF(MYSIZE>1) THEN

        !First send data to left and right processor
        IF(MYRANK>0) THEN
          CALL MPI_SENDRECV(NL_REQ,1,MPI_INTEGER,MYLEFT,MYRANK*NVS,
     &                      NL_SUP,1,MPI_INTEGER,MYLEFT,MYLEFT*NVS+10,MPI_COMM_EDDY,STATUS,IERR)
          IF(NL_REQ>0) THEN
            ALLOCATE(XL_REQ(NL_REQ),YL_REQ(NL_REQ),ZL_REQ(NL_REQ)
     &              ,PL_REQ(NL_REQ),DL_REQ(NL_REQ))
            DO I=1,NL_REQ
              XL_REQ(I)=XL(I)
              YL_REQ(I)=YL(I)
              ZL_REQ(I)=ZL(I)
              DL_REQ(I)=DIRL(I)
            ENDDO
            CALL MPI_SEND(XL_REQ,NL_REQ,MTYPE,MYLEFT,MYRANK*NVS+1,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(YL_REQ,NL_REQ,MTYPE,MYLEFT,MYRANK*NVS+2,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(ZL_REQ,NL_REQ,MTYPE,MYLEFT,MYRANK*NVS+3,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(DL_REQ,NL_REQ,MPI_INTEGER,MYLEFT,MYRANK*NVS+4,MPI_COMM_EDDY,IERR)
          ENDIF
        ENDIF

        IF(MYRANK<MYSIZE-1) THEN
          CALL MPI_SENDRECV(NR_REQ,1,MPI_INTEGER,MYRIGHT,MYRANK*NVS+10,
     &                      NR_SUP,1,MPI_INTEGER,MYRIGHT,MYRIGHT*NVS,MPI_COMM_EDDY,STATUS,IERR)
          IF(NR_REQ>0) THEN
            ALLOCATE(XR_REQ(NR_REQ),YR_REQ(NR_REQ),ZR_REQ(NR_REQ)
     &              ,PR_REQ(NR_REQ),DR_REQ(NR_REQ))
            DO I=1,NR_REQ
              XR_REQ(I)=XR(I)
              YR_REQ(I)=YR(I)
              ZR_REQ(I)=ZR(I)
              DR_REQ(I)=DIRR(I)
            ENDDO
            CALL MPI_SEND(XR_REQ,NR_REQ,MTYPE,MYRIGHT,MYRANK*NVS+11,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(YR_REQ,NR_REQ,MTYPE,MYRIGHT,MYRANK*NVS+12,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(ZR_REQ,NR_REQ,MTYPE,MYRIGHT,MYRANK*NVS+13,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(DR_REQ,NR_REQ,MPI_INTEGER,MYRIGHT,MYRANK*NVS+14,MPI_COMM_EDDY,IERR)
          ENDIF
        ENDIF

        !Supply info if it has been requested
        IF(MYRANK>0) THEN
          IF(NL_SUP>0) THEN
            ALLOCATE(XL_SUP(NL_SUP),YL_SUP(NL_SUP),ZL_SUP(NL_SUP)
     &              ,PL_SUP(NL_SUP),DL_SUP(NL_SUP))
            CALL MPI_RECV(XL_SUP,NL_SUP,MTYPE,MYLEFT,MYLEFT*NVS+11,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(YL_SUP,NL_SUP,MTYPE,MYLEFT,MYLEFT*NVS+12,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(ZL_SUP,NL_SUP,MTYPE,MYLEFT,MYLEFT*NVS+13,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(DL_SUP,NL_SUP,MPI_INTEGER,MYLEFT,MYLEFT*NVS+14,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NL_SUP
              PL_SUP(II) = interp_cellface(xl_sup(ii),yl_sup(ii),zl_sup(ii)
     &            ,xc,yc,zc,p,nx,ny,nz,icyl,dl_sup(ii))
            ENDDO
            CALL MPI_SEND(PL_SUP,NL_SUP,MTYPE,MYLEFT,MYRANK*NVS+15,MPI_COMM_EDDY,IERR)
          ENDIF
        ENDIF


        IF(MYRANK<MYSIZE-1) THEN
          IF(NR_SUP>0) THEN
            ALLOCATE(XR_SUP(NR_SUP),YR_SUP(NR_SUP),ZR_SUP(NR_SUP)
     &              ,PR_SUP(NR_SUP),DR_SUP(NR_SUP))
            CALL MPI_RECV(XR_SUP,NR_SUP,MTYPE,MYRIGHT,MYRIGHT*NVS+1,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(YR_SUP,NR_SUP,MTYPE,MYRIGHT,MYRIGHT*NVS+2,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(ZR_SUP,NR_SUP,MTYPE,MYRIGHT,MYRIGHT*NVS+3,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(DR_SUP,NR_SUP,MPI_INTEGER,MYRIGHT,MYRIGHT*NVS+4,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NR_SUP
              PR_SUP(II) = interp_cellface(xr_sup(ii),yr_sup(ii),zr_sup(ii)
     &            ,xc,yc,zc,p,nx,ny,nz,icyl,dr_sup(ii))
            ENDDO
            CALL MPI_SEND(PR_SUP,NR_SUP,MTYPE,MYRIGHT,MYRANK*NVS+5,MPI_COMM_EDDY,IERR)
          ENDIF

        ENDIF

       !Receive data from left and right processor
        IF(MYRANK>0) THEN
          IF(NL_REQ>0) THEN
            CALL MPI_RECV(PL_REQ,NL_REQ,MTYPE,MYLEFT,MYLEFT*NVS+5,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NL_REQ
              im = iml(ii)
              dpdss(im) = calc_dpdn(nxp(im,1),nyp(im,1),nzp(im,1)
     &             ,nxp(im,2),nyp(im,2),nzp(im,2),pim(im),pl_req(ii),0)
            ENDDO
            DEALLOCATE(XL_REQ,YL_REQ,ZL_REQ,PL_REQ,DL_REQ)
          ENDIF
          IF(NL_SUP>0) DEALLOCATE(XL_SUP,YL_SUP,ZL_SUP,PL_SUP,DL_SUP)
        ENDIF

        IF(MYRANK<MYSIZE-1) THEN
          IF(NR_REQ>0) THEN
            CALL MPI_RECV(PR_REQ,NR_REQ,MTYPE,MYRIGHT,MYRIGHT*NVS+15,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NR_REQ
              im = imr(ii)
              dpdss(im) = calc_dpdn(nxp(im,1),nyp(im,1),nzp(im,1)
     &             ,nxp(im,2),nyp(im,2),nzp(im,2),pim(im),pr_req(ii),0)

           ENDDO
           DEALLOCATE(XR_REQ,YR_REQ,ZR_REQ,PR_REQ,DR_REQ)
          ENDIF
          IF(NR_SUP>0) DEALLOCATE(XR_SUP,YR_SUP,ZR_SUP,PR_SUP,DR_SUP)
        ENDIF

      ENDIF

c      NL_REQ = 0
c      NR_REQ = 0
      IL = 0
      IR = 0
      
      DO img=1,ngrid

        im=igrd(img)
        I=IP(IM)
        J=JP(IM)
        K=KP(IM)

        I1=IP(IM)+NXP(IM,1)
        J1=JP(IM)+NYP(IM,1)
        K1=KP(IM)+NZP(IM,1)
        I2=I1+NXP(IM,1)
        J2=J1+NYP(IM,1)
        K2=K1+NZP(IM,1)

        IF(I2<1) THEN
          open(unit=10,file='stats_imb.dat',form='formatted'
     &         ,position='append')
          write(10,'(2(A,3(1x,I4)))')
     &         'WARNING: 2nd pres. stencil point past centerline',i2,j2,k2
     &        ,' changing to',abs(i2)+1,jsym(j2),k2 
          close(10)

          call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
          call index_bnd(i2,j2,k2,i2,j2,k2,nx,ny,nz)
c          I2=ABS(I2)+1
c          J2=JSYM(J2)
        ENDIF


        IF(K2<1) THEN
          IL=IL+1
          LI(1,IL)=I2
          LI(2,IL)=J2
          LI(3,IL)=NZ-2-K2
          IML(IL)=IM
        ELSEIF(K2>NZ) THEN
          IR=IR+1
          RI(1,IR)=I2
          RI(2,IR)=J2
          RI(3,IR)=K2+2-NZ
          IMR(IR)=IM
        ELSE
          dpdss(im) = calc_dpdn(xc(i1),yc(j1),zc(k1),xc(i2),yc(j2),zc(k2),p(i1,j1,k1),p(i2,j2,k2),icyl)
        ENDIF

      ENDDO

      NL_REQ = IL !Number of points required from left processor
      NR_REQ = IR !Number of points required from right processor


      NVS=10

c....Check if info from adjacent domain is needed
      IF(MYSIZE>1) THEN

        !First send data to left and right processor
        IF(MYRANK>0) THEN
          CALL MPI_SENDRECV(NL_REQ,1,MPI_INTEGER,MYLEFT,MYRANK*NVS,
     &                      NL_SUP,1,MPI_INTEGER,MYLEFT,MYLEFT*NVS+5,MPI_COMM_EDDY,STATUS,IERR)
          IF(NL_REQ>0) THEN
            ALLOCATE(IL_REQ(3,NL_REQ),PL_REQ(NL_REQ),ZL_REQ(NL_REQ))
            DO I=1,NL_REQ
              IL_REQ(:,I)=LI(:,I)
            ENDDO
            CALL MPI_SEND(IL_REQ,3*NL_REQ,MPI_INTEGER,MYLEFT,MYRANK*NVS+1,MPI_COMM_EDDY,IERR)
          ENDIF
        ENDIF

        IF(MYRANK<MYSIZE-1) THEN
          CALL MPI_SENDRECV(NR_REQ,1,MPI_INTEGER,MYRIGHT,MYRANK*NVS+5,
     &                      NR_SUP,1,MPI_INTEGER,MYRIGHT,MYRIGHT*NVS,MPI_COMM_EDDY,STATUS,IERR)
          IF(NR_REQ>0) THEN
            ALLOCATE(IR_REQ(3,NR_REQ),PR_REQ(NR_REQ),ZR_REQ(NR_REQ))
            DO I=1,NR_REQ
              IR_REQ(:,I)=RI(:,I)
            ENDDO
            CALL MPI_SEND(IR_REQ,3*NR_REQ,MPI_INTEGER,MYRIGHT,MYRANK*NVS+6,MPI_COMM_EDDY,IERR)
          ENDIF
        ENDIF

        !Supply info if it has been requested
        IF(MYRANK>0) THEN
          IF(NL_SUP>0) THEN
            ALLOCATE(IL_SUP(3,NL_SUP),PL_SUP(NL_SUP),ZL_SUP(NL_SUP))
            CALL MPI_RECV(IL_SUP,3*NL_SUP,MPI_INTEGER,MYLEFT,MYLEFT*NVS+6,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NL_SUP
              I=IL_SUP(1,II)
              J=IL_SUP(2,II)
              K=IL_SUP(3,II)
              PL_SUP(II)=P(I,J,K)
              ZL_SUP(II)=ZC(K)
            ENDDO
            CALL MPI_SEND(PL_SUP,NL_SUP,MTYPE,MYLEFT,MYRANK*NVS+7,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(ZL_SUP,NL_SUP,MTYPE,MYLEFT,MYRANK*NVS+8,MPI_COMM_EDDY,IERR)
          ENDIF
        ENDIF


        IF(MYRANK<MYSIZE-1) THEN
          IF(NR_SUP>0) THEN
            ALLOCATE(IR_SUP(3,NR_SUP),PR_SUP(NR_SUP),ZR_SUP(NR_SUP))
            CALL MPI_RECV(IR_SUP,3*NR_SUP,MPI_INTEGER,MYRIGHT,MYRIGHT*NVS+1,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NR_SUP
              I=IR_SUP(1,II)
              J=IR_SUP(2,II)
              K=IR_SUP(3,II)
              PR_SUP(II)=P(I,J,K)
              ZR_SUP(II)=ZC(K)
            ENDDO
            CALL MPI_SEND(PR_SUP,NR_SUP,MTYPE,MYRIGHT,MYRANK*NVS+2,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(ZR_SUP,NR_SUP,MTYPE,MYRIGHT,MYRANK*NVS+3,MPI_COMM_EDDY,IERR)
          ENDIF

        ENDIF

        !Receive data from left and right processor
        IF(MYRANK>0) THEN
          IF(NL_REQ>0) THEN
            CALL MPI_RECV(PL_REQ,NL_REQ,MTYPE,MYLEFT,MYLEFT*NVS+2,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(ZL_REQ,NL_REQ,MTYPE,MYLEFT,MYLEFT*NVS+3,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NL_REQ
              im = iml(ii)
              i  = ip(im)
              j  = jp(im)
              k  = kp(im)
              i1 = i+nxp(im,1)
              j1 = j+nyp(im,1)
              k1 = k+nzp(im,1)
              i2 = i1+nxp(im,1)
              j2 = j1+nyp(im,1)
              dpdss(im) = calc_dpdn(xc(i1),yc(j1),zc(k1),xc(i2),yc(j2),zl_req(ii),p(i1,j1,k1),pl_req(ii),icyl)

           ENDDO
           DEALLOCATE(IL_REQ,PL_REQ,ZL_REQ)
          ENDIF
          IF(NL_SUP>0) DEALLOCATE(IL_SUP,PL_SUP,ZL_SUP)
        ENDIF


        IF(MYRANK<MYSIZE-1) THEN
          IF(NR_REQ>0) THEN
            CALL MPI_RECV(PR_REQ,NR_REQ,MTYPE,MYRIGHT,MYRIGHT*NVS+7,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(ZR_REQ,NR_REQ,MTYPE,MYRIGHT,MYRIGHT*NVS+8,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NR_REQ
              im = imr(ii)
              i  = ip(im)
              j  = jp(im)
              k  = kp(im)
              i1 = i+nxp(im,1)
              j1 = j+nyp(im,1)
              k1 = k+nzp(im,1)
              i2 = i1+nxp(im,1)
              j2 = j1+nyp(im,1)
              dpdss(im) = calc_dpdn(xc(i1),yc(j1),zc(k1),xc(i2),yc(j2),zr_req(ii),p(i1,j1,k1),pr_req(ii),icyl)

           ENDDO
           DEALLOCATE(IR_REQ,PR_REQ,ZR_REQ)
          ENDIF
          IF(NR_SUP>0) DEALLOCATE(IR_SUP,PR_SUP,ZR_SUP)
        ENDIF

      ENDIF


      il=0
      ir=0
      DO idg=1,ndiag

        im = idia(idg)
        i = ip(im)
        j = jp(im)
        k = kp(im)

        i1=ip(im)+int(nxp(im,1))
        j1=jp(im)+int(nyp(im,1))
        k1=kp(im)+int(nzp(im,1))

        x1car = xc(i1)
        y1car = yc(j1)
        if(icyl==1) then
          x1car = xc(i1)*cos(yc(j1))
          y1car = xc(i1)*sin(yc(j1))
        endif

        IF(i1<1) THEN
          call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
          open(unit=10,file='stats_imb.dat',form='formatted'
     &         ,position='append')
          write(10,'(A,3(1x,I4))')
     &         'WARNING: 2nd pres. stencil point past centerline
     &         , changing to',i1,j1,k1
          close(10)

        ENDIF
 
        pim(im) = p(i1,j1,k1)

        if(nzp(im,2)<zc(1)) then
          il = il+1
          iml(il) = im
          xl(il) = nxp(im,2)
          yl(il) = nyp(im,2)
          zl(il) = nzp(im,2)
          dirl(il) = dir(im,1)
        elseif(nzp(im,2)>zc(nz)) then
          ir = ir+1
          imr(ir) = im
          xr(ir) = nxp(im,2)
          yr(ir) = nyp(im,2)
          zr(ir) = nzp(im,2)
          dirr(ir) = dir(im,1)
        else
          pint2 = interp_cellface(nxp(im,2),nyp(im,2),nzp(im,2)
     &          ,xc,yc,zc,p,nx,ny,nz,icyl,dir(im,1))
          dpdss(im) = calc_dpdn(x1car,y1car,zc(k1)
     &         ,nxp(im,2),nyp(im,2),nzp(im,2),pim(im),pint2,0)
        endif


      ENDDO


      NL_REQ = IL !Number of points required from left processor
      NR_REQ = IR !Number of points required from right processor


c      write(6,*) 'diag dpdss:',nl_req,nr_req

      NVS=20

c....Check if info from adjacent domain is needed
      IF(MYSIZE>1) THEN

        !First send data to left and right processor
        IF(MYRANK>0) THEN
          CALL MPI_SENDRECV(NL_REQ,1,MPI_INTEGER,MYLEFT,MYRANK*NVS,
     &                      NL_SUP,1,MPI_INTEGER,MYLEFT,MYLEFT*NVS+10,MPI_COMM_EDDY,STATUS,IERR)
          IF(NL_REQ>0) THEN
            ALLOCATE(XL_REQ(NL_REQ),YL_REQ(NL_REQ),ZL_REQ(NL_REQ)
     &              ,PL_REQ(NL_REQ),DL_REQ(NL_REQ))
            DO I=1,NL_REQ
              XL_REQ(I)=XL(I)
              YL_REQ(I)=YL(I)
              ZL_REQ(I)=ZL(I)
              DL_REQ(I)=DIRL(I)
            ENDDO
            CALL MPI_SEND(XL_REQ,NL_REQ,MTYPE,MYLEFT,MYRANK*NVS+1,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(YL_REQ,NL_REQ,MTYPE,MYLEFT,MYRANK*NVS+2,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(ZL_REQ,NL_REQ,MTYPE,MYLEFT,MYRANK*NVS+3,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(DL_REQ,NL_REQ,MPI_INTEGER,MYLEFT,MYRANK*NVS+4,MPI_COMM_EDDY,IERR)
          ENDIF
        ENDIF

        IF(MYRANK<MYSIZE-1) THEN
          CALL MPI_SENDRECV(NR_REQ,1,MPI_INTEGER,MYRIGHT,MYRANK*NVS+10,
     &                      NR_SUP,1,MPI_INTEGER,MYRIGHT,MYRIGHT*NVS,MPI_COMM_EDDY,STATUS,IERR)
          IF(NR_REQ>0) THEN
            ALLOCATE(XR_REQ(NR_REQ),YR_REQ(NR_REQ),ZR_REQ(NR_REQ)
     &              ,PR_REQ(NR_REQ),DR_REQ(NR_REQ))
            DO I=1,NR_REQ
              XR_REQ(I)=XR(I)
              YR_REQ(I)=YR(I)
              ZR_REQ(I)=ZR(I)
              DR_REQ(I)=DIRR(I)
            ENDDO
            CALL MPI_SEND(XR_REQ,NR_REQ,MTYPE,MYRIGHT,MYRANK*NVS+11,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(YR_REQ,NR_REQ,MTYPE,MYRIGHT,MYRANK*NVS+12,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(ZR_REQ,NR_REQ,MTYPE,MYRIGHT,MYRANK*NVS+13,MPI_COMM_EDDY,IERR)
            CALL MPI_SEND(DR_REQ,NR_REQ,MPI_INTEGER,MYRIGHT,MYRANK*NVS+14,MPI_COMM_EDDY,IERR)
          ENDIF
        ENDIF

        !Supply info if it has been requested
        IF(MYRANK>0) THEN
          IF(NL_SUP>0) THEN
            ALLOCATE(XL_SUP(NL_SUP),YL_SUP(NL_SUP),ZL_SUP(NL_SUP)
     &              ,PL_SUP(NL_SUP),DL_SUP(NL_SUP))
            CALL MPI_RECV(XL_SUP,NL_SUP,MTYPE,MYLEFT,MYLEFT*NVS+11,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(YL_SUP,NL_SUP,MTYPE,MYLEFT,MYLEFT*NVS+12,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(ZL_SUP,NL_SUP,MTYPE,MYLEFT,MYLEFT*NVS+13,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(DL_SUP,NL_SUP,MPI_INTEGER,MYLEFT,MYLEFT*NVS+14,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NL_SUP
              PL_SUP(II) = interp_cellface(xl_sup(ii),yl_sup(ii),zl_sup(ii)
     &            ,xc,yc,zc,p,nx,ny,nz,icyl,dl_sup(ii))
            ENDDO
            CALL MPI_SEND(PL_SUP,NL_SUP,MTYPE,MYLEFT,MYRANK*NVS+15,MPI_COMM_EDDY,IERR)
          ENDIF
        ENDIF


        IF(MYRANK<MYSIZE-1) THEN
          IF(NR_SUP>0) THEN
            ALLOCATE(XR_SUP(NR_SUP),YR_SUP(NR_SUP),ZR_SUP(NR_SUP)
     &              ,PR_SUP(NR_SUP),DR_SUP(NR_SUP))
            CALL MPI_RECV(XR_SUP,NR_SUP,MTYPE,MYRIGHT,MYRIGHT*NVS+1,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(YR_SUP,NR_SUP,MTYPE,MYRIGHT,MYRIGHT*NVS+2,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(ZR_SUP,NR_SUP,MTYPE,MYRIGHT,MYRIGHT*NVS+3,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(DR_SUP,NR_SUP,MPI_INTEGER,MYRIGHT,MYRIGHT*NVS+4,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NR_SUP
              PR_SUP(II) = interp_cellface(xr_sup(ii),yr_sup(ii),zr_sup(ii)
     &            ,xc,yc,zc,p,nx,ny,nz,icyl,dr_sup(ii))
            ENDDO
            CALL MPI_SEND(PR_SUP,NR_SUP,MTYPE,MYRIGHT,MYRANK*NVS+5,MPI_COMM_EDDY,IERR)
          ENDIF

        ENDIF

       !Receive data from left and right processor
        IF(MYRANK>0) THEN
          IF(NL_REQ>0) THEN
            CALL MPI_RECV(PL_REQ,NL_REQ,MTYPE,MYLEFT,MYLEFT*NVS+5,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NL_REQ
              im = iml(ii)
              i1=ip(im)+nxp(im,1)
              j1=jp(im)+nyp(im,1)
              k1=kp(im)+nzp(im,1)

              call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)

              x1car = xc(i1)
              y1car = yc(j1)
              if(icyl==1) then
                x1car = xc(i1)*cos(yc(j1))
                y1car = xc(i1)*sin(yc(j1))
              endif
              dpdss(im) = calc_dpdn(x1car,y1car,zc(k1)
     &             ,nxp(im,2),nyp(im,2),nzp(im,2),pim(im),pl_req(ii),0)
            ENDDO
            DEALLOCATE(XL_REQ,YL_REQ,ZL_REQ,PL_REQ,DL_REQ)
          ENDIF
          IF(NL_SUP>0) DEALLOCATE(XL_SUP,YL_SUP,ZL_SUP,PL_SUP,DL_SUP)
        ENDIF

        IF(MYRANK<MYSIZE-1) THEN
          IF(NR_REQ>0) THEN
            CALL MPI_RECV(PR_REQ,NR_REQ,MTYPE,MYRIGHT,MYRIGHT*NVS+15,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NR_REQ
              im = imr(ii)
              i1=ip(im)+nxp(im,1)
              j1=jp(im)+nyp(im,1)
              k1=kp(im)+nzp(im,1)

              call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)

              x1car = xc(i1)
              y1car = yc(j1)
              if(icyl==1) then
                x1car = xc(i1)*cos(yc(j1))
                y1car = xc(i1)*sin(yc(j1))
              endif
              dpdss(im) = calc_dpdn(x1car,y1car,zc(k1)
     &             ,nxp(im,2),nyp(im,2),nzp(im,2),pim(im),pr_req(ii),0)
           ENDDO
           DEALLOCATE(XR_REQ,YR_REQ,ZR_REQ,PR_REQ,DR_REQ)
          ENDIF
          IF(NR_SUP>0) DEALLOCATE(XR_SUP,YR_SUP,ZR_SUP,PR_SUP,DR_SUP)
        ENDIF

      ENDIF

      DEALLOCATE(inrm,igrd,idia)

      RETURN
      END
C------------------------------------------------------------------------


C---- function calc_dpdn-----------------------N. Beratlis-11 Jan 2009---
C
C     PURPOSE: Calculate pressure gradient of variable p given values
C     at points (x1,y1,z1) and (x2,y2,z2).
C
C------------------------------------------------------------------------
      REAL function calc_dpdn(x1,y1,z1,x2,y2,z2,p1,p2,icyl)
C
      IMPLICIT NONE
      INTEGER icyl
      REAL    x1,y1,z1,x2,y2,z2,p1,p2

      REAL    ds

      if(icyl==0) then
        ds = sqrt((x2-x1)**2. + (y2-y1)**2. + (z2-z1)**2.)
      else
        ds = sqrt((x2*cos(y2)-x1*cos(y1))**2.
     &          + (x2*sin(y2)-x1*sin(y1))**2.
     &          + (z2-z1)**2. )
      endif

      calc_dpdn = (p2-p1)/ds

      RETURN

      END
C------------------------------------------------------------------------


C---- function calc_dpdn1 ---------------------N. Beratlis-11 Jan 2009---
C
C     PURPOSE: Calculate pressure gradient of variable p given values
C     at points (x1,y1,z1) and (x2,y2,z2).
C
C------------------------------------------------------------------------
      REAL function calc_dpdn1(x1,y1,z1,x2,y2,z2,p1,p2,icyl,ie)
C
      IMPLICIT NONE
      INTEGER icyl,ie
      REAL    x1,y1,z1,x2,y2,z2,p1,p2

      REAL    ds

      if(icyl==0) then
        ds = sqrt((x2-x1)**2. + (y2-y1)**2. + (z2-z1)**2.)
      else
        ds = sqrt((x2*cos(y2)-x1*cos(y1))**2.
     &          + (x2*sin(y2)-x1*sin(y1))**2.
     &          + (z2-z1)**2. )
      endif

      calc_dpdn1 = (p2-p1)/ds

c      if(ie==385) then
c        write(6,*) 'dpdn:',ie,x1,y1,z1,x2,y2,z2,p1,p2,ds,calc_dpdn1
c      endif

      RETURN

      END
C------------------------------------------------------------------------


C---- function body2grid-----------------N. Beratlis-14 Jan. 2009--
C
C     PURPOSE: From a point q on the immersed boundary surface find 1
C     intersection along surface normal r and interpolate value of uo
C     at the intersection
C
C-----------------------------------------------------------------------
      integer function body2grid(q,rvec,xu_car,yu_car,xu,yu,zu
     &     ,nx,ny,nz,xint,yint,zint,i1,j1,k1,dir,icyl)

      IMPLICIT NONE

      REAL, PARAMETER :: pi=3.141592653589793
c
c.... Input/Output variables
      integer nx,ny,nz,dir,icyl
      real    xint,yint,zint
      real    q(3),rvec(3)
      real    xu(nx),yu(ny),zu(nz)
      real    xu_car(nx,ny),yu_car(nx,ny)
c
c.... Local variables
      integer i,j,k,i1,j1,k1,iext,jext,kext
      integer iflag1,iflag2,iflag3
      real    xp,yp,zp
      real    dy,rm,dz
c
c.... Functions
      integer ray_face_int

      xp = q(1)
      yp = q(2)
      zp = q(3)
c IJK_XYZ finds the indices of the cell containing the input coordinates (xp,yp,zp)
      call ijk_xyz(xp,yp,zp,xu,yu,zu,nx,ny,nz,i,j,k,icyl)

c      theta = yu(j)
      dy = 2.0*pi/real(ny-2)
      rm = xu(i)
      if(icyl==1 .AND. i==1) rm=0.0
c VEC_IJKEXT finds the extensions of the grid indices along the vector RVEC
      call vec_ijkext(q,rvec,iext,jext,kext,icyl,1,rm,dy)

      i1 = i+iext
      j1 = j
      k1 = k
      dir = 1
c RAY_FACE_INT looks for an intersection of the vector Q+RVEC with
c a grid face normal to the X direction
c the coordinates of the eventual intersection are stored in XINT,YINT,ZINT
      iflag1 = ray_face_int(q,rvec,xu_car,yu_car,zu,nx,ny,nz
     &     ,i1,j1,k1,1,xint,yint,zint,icyl)

      if(iflag1==0) then
c the intersection has not been found with DIR=1 (faces normal to X)
        i1 = i
        j1 = j+jext
        k1 = k
        dir = 2
        iflag2 = ray_face_int(q,rvec,xu_car,yu_car,zu,nx,ny,nz
     &     ,i,j+jext,k,2,xint,yint,zint,icyl)
        if(iflag2==0) then
c the intersection has not been found also with DIR=2 (faces normal to Y)
          dir = 3
          k1 = k+kext
          dz = zu(k1)-zu(k)
          call vec_ijkext(q,rvec,iext,jext,kext,icyl,3,dz,dy)
          iflag3 = ray_face_int(q,rvec,xu_car,yu_car,zu,nx,ny,nz
     &         ,i,j,k+kext,3,xint,yint,zint,icyl)
          i1 = i
          j1 = j
          k1 = k+kext
        endif
      endif

      if(iflag1==0.AND.iflag2==0.AND.iflag3==0) then
        body2grid = 0	 ! no intersection found
      else
        body2grid = 1	 ! an intersection has been found
      endif

      return

      end
C-----------------------------------------------------------------------




C---- function vel_body2fluid-----------------N. Beratlis-14 Jan. 2009--
C
C     PURPOSE: From a point q on the immersed boundary surface find 1
C     intersection along surface normal r and interpolate value of uo
C     at the intersection
C
C-----------------------------------------------------------------------
      integer function vel_body2fluid(q,rvec,xu_car,yu_car,xu,yu,zu,zug,flaguo
     &     ,uo,nx,ny,nz,nzg,xint,yint,zint,uint,n1,nbd,myrank,icyl,ifacet,amin)
c
      IMPLICIT NONE

      REAL, PARAMETER :: pi=3.141592653589793
c
c.... Input/Output Arrays      
      INTEGER nx,ny,nz,nzg,dir,nbd,n1,icyl,myrank,ifacet
      REAL    xint,yint,zint
      INTEGER flaguo(nx,ny,nz,nbd)
      REAL    q(3),rvec(3)
      REAL    uo(nx,ny,nz)
      REAL    xu_car(nx,ny),yu_car(nx,ny),xu(nx),yu(ny),zu(nz),zug(nzg)
c
c.... Local Arrays
      INTEGER i,j,k,kg,i1,j1,k1,k1g,i2,j2,k2,k2g,iext,jext,kext,next
      INTEGER intrs,ii,iphy,iflagt,iflag,iflag1,iflag2,iflag3
      REAL    xp,yp,zp,xint1,yint1,zint1,xint2,yint2,zint2,uint
      REAL    a,theta,dy,rm,dz,s,ds,amin

      INTEGER nclock
      real    tclock,clocktemp
      REAL    clock(10)
c
c.... Functions      
      INTEGER ray_face_int,physicalface,grid2grid_intr
      REAL    interp_cellface,extmag,anglerad,mindxdydz

      clock = 0.0
      nclock = 10
c      amin = 0.01
      xint1 = 0.0
      yint1 = 0.0
      zint1 = 0.0

      xp = q(1)
      yp = q(2)
      zp = q(3)
      
      iphy = 0
      intrs = 1


      clocktemp = tclock()
      call ijk_xyz(xp,yp,zp,xu,yu,zug,nx,ny,nzg,i,j,k,icyl)

      ds = mindxdydz(xu,yu,zug,nx,ny,nzg,i,j,k,icyl)

      theta = yu(j)
      dy = 2.0*pi/real(ny-2)
      rm = xu(i)
      if(icyl==1 .AND. i==1) rm=0.0
      call vec_ijkext(q,rvec,iext,jext,kext,icyl,1,rm,dy)

      i1 = i+iext
      j1 = j
      k1 = k

      clock(1) = clock(1) + tclock() - clocktemp

      clocktemp = tclock()
      iflag1 = ray_face_int(q,rvec,xu_car,yu_car,zug,nx,ny,nzg
     &     ,i1,j1,k1,1,xint1,yint1,zint1,icyl)
      clock(2) = clock(2) + tclock() - clocktemp

      !Compute distance from body to intersection
      s = sqrt( (q(1)-xint1)**2. + (q(2)-yint1)**2. + (q(3)-zint1)**2. )

      if(iflag1==1) then
        i1 = i+iext
        j1 = j
        k1g = k
        dir = 1
        if(s>amin*ds) then

          if(zint1>=zu(1).AND.zint1<zu(nz-1)) then
            k1 = k1g-myrank*(nz-2)
            iflag = physicalface(flaguo,nx,ny,nz,i1,j1,k1,nbd,1)
            if(iflag/=0) then
              iphy = iphy+1
              uint = interp_cellface(xint1,yint1,zint1,xu,yu,zu,uo
     &                ,nx,ny,nz,icyl,dir)
              xint = xint1
              yint = yint1
              zint = zint1
            endif
          endif

        endif

        if(iphy==0) then
          iflagt=0
          do while(iflagt==0.AND.intrs<n1)
            intrs = intrs+1

c            if(ifacet==16378) then
c              write(6,*) '1. vel_body2fluid',ifacet,rvec
c            endif

            if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &         ,rvec,xu_car,yu_car,zug,nx,ny,nzg
     &         ,i1,j1,k1g,i2,j2,k2g,dir,dir,icyl,ifacet)==1) then
              if(zint2>=zu(1).AND.zint2<zu(nz-1)) then
                k2 = k2g-myrank*(nz-2)                
                if(physicalface(flaguo,nx,ny,nz,i2,j2,k2,nbd,dir)/=0) then
                  iphy = iphy+1
                  iflagt=1
                  uint = interp_cellface(xint2,yint2,zint2,xu,yu,zu,uo
     &                 ,nx,ny,nz,icyl,dir)
                  xint = xint2
                  yint = yint2
                  zint = zint2
                endif
              endif
            endif
            xint1 = xint2
            yint1 = yint2
            zint1 = zint2  
            i1 = i2
            j1 = j2
            k1g= k2g
            if(intrs>=n1) iflagt=1
          enddo
        endif

      else

        i1 = i
        j1 = j+jext
        k1 = k

        iflag2 = ray_face_int(q,rvec,xu_car,yu_car,zug,nx,ny,nzg
     &     ,i,j+jext,k,2,xint1,yint1,zint1,icyl)
        s = sqrt( (q(1)-xint1)**2. + (q(2)-yint1)**2. + (q(3)-zint1)**2. )

c        write(6,*) 's=',s
c        write(6,*) 'iflag2=',iflag2

        if(iflag2==1) then
          i1 = i
          j1 = j+jext
          k1g = k
          dir = 2

          if(s>amin*ds) then
          if(zint1>=zu(1).AND.zint1<zu(nz-1)) then
            k1 = k1g-myrank*(nz-2)

            if(physicalface(flaguo,nx,ny,nz,i1,j1,k1,nbd,2)/=0) then
              iphy = iphy+1
              uint = interp_cellface(xint1,yint1,zint1,xu,yu,zu,uo
     &             ,nx,ny,nz,icyl,dir)
              xint = xint1
              yint = yint1
              zint = zint1
            endif
          endif
          endif

          if(iphy==0) then
            iflagt=0
            do while(iflagt==0.AND.intrs<n1)
              intrs = intrs+1
              
c              if(ifacet==16378) then
c                write(6,*) '2. vel_body2fluid',ifacet,rvec
c              endif

              if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &             ,rvec,xu_car,yu_car,zug,nx,ny,nzg
     &             ,i1,j1,k1g,i2,j2,k2g,dir,dir,icyl,ifacet)==1) then
                if(zint2>=zu(1).AND.zint2<zu(nz-1)) then
                  k2 = k2g-myrank*(nz-2)

                  if(physicalface(flaguo,nx,ny,nz,i2,j2,k2,nbd,dir)/=0) then
                    iphy = iphy+1
                    iflagt=1
                    uint = interp_cellface(xint2,yint2,zint2,xu,yu,zu,uo
     &                   ,nx,ny,nz,icyl,dir)
                    xint = xint2
                    yint = yint2
                    zint = zint2
                  endif
                endif
              endif
              xint1 = xint2
              yint1 = yint2
              zint1 = zint2  
              i1 = i2
              j1 = j2
              k1g= k2g
              if(intrs>=n1) iflagt=1
            enddo
          endif

        else

          k1g = k+kext 
          dz = zug(k1g)-zug(k)
          call vec_ijkext(q,rvec,iext,jext,kext,icyl,3,dz,dy)

          iflag3 = ray_face_int(q,rvec,xu_car,yu_car,zug,nx,ny,nzg
     &         ,i,j,k+kext,3,xint1,yint1,zint1,icyl)
          s = sqrt( (q(1)-xint1)**2. + (q(2)-yint1)**2. + (q(3)-zint1)**2. )

          if(iflag3==1) then
            i1 = i
            j1 = j
            k1g = k+kext
            dir = 3

            if(s>amin*ds) then
            if(zint1>=zu(1).AND.zint1<zu(nz-1)) then
              k1 = k1g-myrank*(nz-2)
              if(physicalface(flaguo,nx,ny,nz,i1,j1,k1,nbd,3)/=0) then
                iphy = iphy+1
                uint = interp_cellface(xint1,yint1,zint1,xu,yu,zu,uo
     &               ,nx,ny,nz,icyl,dir)
                xint = xint1
                yint = yint1
                zint = zint1
              endif
            endif
            endif

            if(iphy==0) then
              iflagt=0
              do while(iflagt==0.AND.intrs<n1)
                intrs = intrs+1

c                if(ifacet==16378) then
c                  write(6,*) '3. vel_body2fluid',ifacet,rvec
c                endif

                iflag = grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &               ,rvec,xu_car,yu_car,zug,nx,ny,nzg
     &               ,i1,j1,k1g,i2,j2,k2g,dir,dir,icyl,ifacet)

                if(iflag==1) then

                  if(zint2>=zu(1).AND.zint2<zu(nz-1)) then
                    k2 = k2g-myrank*(nz-2)

                    if(physicalface(flaguo,nx,ny,nz,i2,j2,k2,nbd,dir)/=0) then
                      iphy = iphy+1
                      iflagt=1
                      uint = interp_cellface(xint2,yint2,zint2,xu,yu,zu,uo
     &                     ,nx,ny,nz,icyl,dir)
                      xint = xint2
                      yint = yint2
                      zint = zint2
                    endif
                  endif
                endif
                xint1 = xint2
                yint1 = yint2
                zint1 = zint2  
                i1 = i2
                j1 = j2
                k1g= k2g
                if(intrs>=n1) iflagt=1
              enddo
            endif
          endif
        endif
      endif


      if(iflag1==0.AND.iflag2==0.AND.iflag3==0) then
        write(6,*) ifacet,'vel_body2fluid no face intersection',iphy,q,rvec,amin
c        stop
      endif

c      if(ifacet==766) then
c        write(6,*) 'velbody2fluid, myrank=',myrank,ifacet,iflag1,iflag2,iflag3
c      endif

      if(iphy==1) then
        vel_body2fluid = intrs
      else
        vel_body2fluid = n1+1
      endif


      RETURN

      END

C------------------------------------------------------------------------



C---- function press_body2fluid---------------N. Beratlis-16 Jan. 2009--
C
C     PURPOSE: From a point q on the immersed boundary surface find 2
C     intersections along surface normal r and interpolate value of p
C     at the intersections
C
C-----------------------------------------------------------------------
      integer function press_body2fluid(q,rvec,xc_car,yc_car,xc,yc,zc,zcg
     &     ,flagpo,flagpi,p,nx,ny,nz,nzg,nbd,xint,yint,zint,pint,n1,n2,icyl,myrank,ifacet)

      IMPLICIT NONE
      
      integer nx,ny,nz,nzg,nbd,icyl,dir,n1,n2,myrank,ifacet
      real    xint(2),yint(2),zint(2),pint(2)
      real    q(3),rvec(3)
      integer flagpo(nx,ny,nz,nbd),flagpi(nx,ny,nz,nbd)
      real    p(nx,ny,nz)
      real    xc_car(nx,ny),yc_car(nx,ny),xc(nx),yc(ny),zc(nz),zcg(nzg)
      real    unvect(3)

      integer i,j,k,i1,j1,k1,i2,j2,k2,k1g,k2g
      integer iext,jext,kext,iflagt,iphy,intr,it,intrs
      real    xp,yp,zp,theta,dy
      real    xint1,yint1,zint1,xint2,yint2,zint2
     &       ,pint1,pint2

      integer ray_face_int,physicalface1,grid2grid_intr
      real    interp_cellface

      xp = q(1)
      yp = q(2)
      zp = q(3)

      iphy=0
      intr=0

      xint = 0.0
      yint = 0.0
      zint = 0.0
      pint = 0.0
      xint1 = 0.0
      yint1 = 0.0
      zint1 = 0.0
      xint2 = 0.0
      yint2 = 0.0
      zint2 = 0.0
      pint1 = 0.0
      pint2 = 0.0

      call ijk_xyz(xp,yp,zp,xc,yc,zcg,nx,ny,nzg,i,j,k,icyl)

      theta = yc(j)
      dy = yc(2)-yc(1)
      call vec_iext(q,rvec,iext,jext,kext,icyl,theta,dy)

      if(ray_face_int(q,rvec,xc_car,yc_car,zcg,nx,ny,nzg
     &     ,i+iext,j,k,1,xint1,yint1,zint1,icyl)==1) then
        intr = intr+1
        i1 = i+iext
        j1 = j
        k1g = k
        dir = 1
        if(zint1>=zc(1).AND.zint1<zc(nz-1)) then
          k1 = k1g-myrank*(nz-2)
          if(physicalface1(flagpo,flagpi,nx,ny,nz,i1,j1,k1,nbd,1)/=0) then
            iphy = iphy+1
            pint1=interp_cellface(xint1,yint1,zint1,xc,yc,zc,p
     &                ,nx,ny,nz,icyl,dir)
            xint(1) = xint1
            yint(1) = yint1
            zint(1) = zint1
            pint(1) = pint1
          endif
        endif

        if(iphy==0) then
          intrs = 0
          iflagt=0
          do while(iflagt==0)
            intrs = intrs+1
            if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &         ,rvec,xc_car,yc_car,zcg,nx,ny,nzg
     &         ,i1,j1,k1g,i2,j2,k2g,dir,dir,icyl,ifacet)==1) then
              if(zint2>=zc(1).AND.zint2<zc(nz-1)) then
                k2 = k2g-myrank*(nz-2)
                if(physicalface1(flagpo,flagpi,nx,ny,nz,i2,j2,k2,nbd,dir)/=0) then
                  iphy = iphy+1
                  iflagt=1
                  pint1=interp_cellface(xint2,yint2,zint2,xc,yc,zc,p
     &                 ,nx,ny,nz,icyl,dir)
                  xint(1) = xint2
                  yint(1) = yint2
                  zint(1) = zint2
                  pint(1) = pint1
                endif
              endif
            endif
            xint1 = xint2
            yint1 = yint2
            zint1 = zint2  
            i1 = i2
            j1 = j2
            k1g= k2g
            if(intrs>=n1) iflagt=1
          enddo
        endif

        if(iphy==1) then
          intrs = 0
          iflagt=0
c          call ijk_xyz(xint1,yint1,zint1,xc,yc,zcg,nx,ny,nzg,i1,j1,k1g,icyl)
          do while(iflagt==0)
            intrs=intrs+1
            if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &        ,rvec,xc_car,yc_car,zcg,nx,ny,nzg
     &        ,i1,j1,k1g,i2,j2,k2g,dir,dir,icyl,ifacet)==1) then
              if(zint2>=zc(1).AND.zint2<zc(nz-1)) then
                k2 = k2g-myrank*(nz-2)
                if(physicalface1(flagpo,flagpi,nx,ny,nz,i2,j2,k2,nbd,dir)/=0) then
                  iphy = iphy+1
                  iflagt=1
                  pint2=interp_cellface(xint2,yint2,zint2,xc,yc,zc,p
     &                 ,nx,ny,nz,icyl,dir)
                  xint(2) = xint2
                  yint(2) = yint2
                  zint(2) = zint2
                  pint(2) = pint2
                endif
              endif
            endif
            xint1=xint2
            yint1=yint2
            zint1=zint2
            i1 = i2
            j1 = j2
            k1g = k2g
            if(intrs>=n2) iflagt=1
          enddo
        endif

      elseif(ray_face_int(q,rvec,xc_car,yc_car,zcg,nx,ny,nzg
     &     ,i,j+jext,k,2,xint1,yint1,zint1,icyl)==1) then
        intr = intr+1
        i1 = i
        j1 = j+jext
        k1g = k
        dir = 2
        if(zint1>=zc(1).AND.zint1<zc(nz-1)) then
          k1 = k1g-myrank*(nz-2)
          if(physicalface1(flagpo,flagpi,nx,ny,nz,i1,j1,k1,nbd,2)/=0) then
            iphy = iphy+1
            pint1=interp_cellface(xint1,yint1,zint1,xc,yc,zc,p
     &                ,nx,ny,nz,icyl,dir)
            xint(1) = xint1
            yint(1) = yint1
            zint(1) = zint1
            pint(1) = pint1
          endif
        endif

        if(iphy==0) then
          intrs = 0
          iflagt=0
          do while(iflagt==0)
            intrs = intrs+1
            if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &         ,rvec,xc_car,yc_car,zcg,nx,ny,nzg
     &         ,i1,j1,k1g,i2,j2,k2g,dir,dir,icyl,ifacet)==1) then
              if(zint2>=zc(1).AND.zint2<zc(nz-1)) then
                k2 = k2g-myrank*(nz-2)
                if(physicalface1(flagpo,flagpi,nx,ny,nz,i2,j2,k2,nbd,dir)/=0) then
                  iphy = iphy+1
                  iflagt=1
                  pint1=interp_cellface(xint2,yint2,zint2,xc,yc,zc,p
     &                 ,nx,ny,nz,icyl,dir)
                  xint(1) = xint2
                  yint(1) = yint2
                  zint(1) = zint2
                  pint(1) = pint1
                endif
              endif
            endif
            xint1 = xint2
            yint1 = yint2
            zint1 = zint2 
            i1 = i2
            j1 = j2
            k1g = k2g
            if(intrs>=n1) iflagt=1
          enddo
        endif

        if(iphy==1) then
          intrs = 0
          iflagt=0
c          call ijk_xyz(xint1,yint1,zint1,xc,yc,zcg,nx,ny,nzg,i1,j1,k1g,icyl)
          do while(iflagt==0)
            intrs=intrs+1
            if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &        ,rvec,xc_car,yc_car,zcg,nx,ny,nzg
     &        ,i1,j1,k1g,i2,j2,k2g,dir,dir,icyl,ifacet)==1) then
              if(zint2>=zc(1).AND.zint2<zc(nz-1)) then
                k2 = k2g-myrank*(nz-2)
                if(physicalface1(flagpo,flagpi,nx,ny,nz,i2,j2,k2,nbd,dir)/=0) then
                  iphy = iphy+1
                  iflagt=1
                  pint2=interp_cellface(xint2,yint2,zint2,xc,yc,zc,p
     &                 ,nx,ny,nz,icyl,dir)
                  xint(2) = xint2
                  yint(2) = yint2
                  zint(2) = zint2
                  pint(2) = pint2
                endif
              endif
            endif
            xint1=xint2
            yint1=yint2
            zint1=zint2
            i1 = i2
            j1 = j2
            k1g = k2g
            if(intrs>=n2) iflagt=1
          enddo
        endif

      elseif(ray_face_int(q,rvec,xc_car,yc_car,zcg,nx,ny,nzg
     &     ,i,j,k+kext,3,xint1,yint1,zint1,icyl)==1) then
        intr = intr+1
        i1 = i
        j1 = j
        k1g = k+kext
        dir = 3
        if(zint1>=zc(1).AND.zint1<zc(nz-1)) then
          k1 = k1g-myrank*(nz-2)
          if(physicalface1(flagpo,flagpi,nx,ny,nz,i1,j1,k1,nbd,3)/=0) then
            iphy = iphy+1
            pint1=interp_cellface(xint1,yint1,zint1,xc,yc,zc,p
     &                ,nx,ny,nz,icyl,dir)
            xint(1) = xint1
            yint(1) = yint1
            zint(1) = zint1
            pint(1) = pint1
          endif
        endif

        if(iphy==0) then
          intrs = 0
          iflagt=0
          do while(iflagt==0)
            intrs = intrs+1
            if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &         ,rvec,xc_car,yc_car,zcg,nx,ny,nzg
     &         ,i1,j1,k1g,i2,j2,k2g,dir,dir,icyl,ifacet)==1) then
c               if(ifacet==498) then
c                 write(6,*) intrs,iflagt,xint1,yint1,zint1,xint2,yint2,zint2,i1,j1,k1,i2,j2,k2g
c               endif

              if(zint2>=zc(1).AND.zint2<zc(nz-1)) then
                k2 = k2g-myrank*(nz-2)
                if(physicalface1(flagpo,flagpi,nx,ny,nz,i2,j2,k2,nbd,dir)/=0) then
                  iphy = iphy+1
                  iflagt=1
                  pint1=interp_cellface(xint2,yint2,zint2,xc,yc,zc,p
     &                 ,nx,ny,nz,icyl,dir)
                  xint(1) = xint2
                  yint(1) = yint2
                  zint(1) = zint2
                  pint(1) = pint1
                endif
              endif
            endif
            xint1 = xint2
            yint1 = yint2
            zint1 = zint2            
            i1 = i2
            j1 = j2
            k1g = k2g
            if(intrs>=n1) iflagt=1
          enddo
        endif

        if(iphy==1) then
          intrs = 0
          iflagt=0
c          call ijk_xyz(xint1,yint1,zint1,xc,yc,zcg,nx,ny,nzg,i1,j1,k1g,icyl)
          do while(iflagt==0)
            intrs=intrs+1
            if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &        ,rvec,xc_car,yc_car,zcg,nx,ny,nzg
     &        ,i1,j1,k1g,i2,j2,k2g,dir,dir,icyl,ifacet)==1) then
              if(zint2>=zc(1).AND.zint2<zc(nz-1)) then
                k2 = k2g-myrank*(nz-2)
                if(physicalface1(flagpo,flagpi,nx,ny,nz,i2,j2,k2,nbd,dir)/=0) then
                  iphy = iphy+1
                  iflagt=1
                  pint2=interp_cellface(xint2,yint2,zint2,xc,yc,zc,p
     &                 ,nx,ny,nz,icyl,dir)
                  xint(2) = xint2
                  yint(2) = yint2
                  zint(2) = zint2
                  pint(2) = pint2
                endif
              endif
            endif
            xint1=xint2
            yint1=yint2
            zint1=zint2
            i1 = i2
            j1 = j2
            k1g = k2g
            if(intrs>=n2) iflagt=1
          enddo
        endif

      endif

      press_body2fluid = iphy

      return

      end
C-----------------------------------------------------------------------

