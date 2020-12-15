C
C----------------------------------------------------------------------- 


 
C-----SUBROUTINE-Grid---------------------------------------------------
C
      SUBROUTINE Grid(xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nzg)
C
C-----------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
C
C-----------------------------------------------------------------------
C
      INTEGER nx,ny,nz,nzg
      REAL xu(nx),yv(ny),zw(nz),zwg(nzg),xc(nx),yc(ny),zc(nz),zcg(nzg)
C
      INTEGER i,j,k

      real rdelx,rdely,rdelz
      real, allocatable, dimension(:) :: cug,cvg

      ALLOCATE(cug(nzg),cvg(nzg))
c
c.....Open file
c
      IF(MYRANK==0) THEN
c
c.....input non-uniform x grid
c.....the grid files contains nx-1 coordinates: the nx coordinate
c.....comes from the ix2 and ix2-1 coordinates
c
        if(igrid==1) then
c
          ix1=2
          OPEN(UNIT=1,FILE=STR1,STATUS='OLD',FORM='FORMATTED')
          read(1,*) ix2
          read(1,*) (j,xu(i),i=1,ix2)
          close(1)

          xlen=xu(ix2)-xu(1)

          xu(nx) = 2.*xu(ix2)-xu(ix2-1)
c     
          xc(ix1:ix2) = .5*(xu(ix1-1:ix2-1)+xu(ix1:ix2))
          xc(1 ) = 2.*xu(1  )-xc(2  )
          xc(nx) = 2.*xu(ix2)-xc(ix2)
c     
          au(ix1:ix2) = 2./(xu(ix1+1:ix2+1)-xu(ix1-1:ix2-1))
          au(1  ) = 2./(-3.*xu(1  )+4.*xu(2    )-xu(3    ))
          au(ix2) = 2./( 3.*xu(ix2)-4.*xu(ix2-1)+xu(ix2-2))
          au(nx ) = au(ix2)
c
          ap(ix1:ix2) = 1./(xu(ix1:ix2)-xu(ix1-1:ix2-1))
          ap(1 ) = ap(ix1)
          ap(nx) = ap(ix2)
c
          av(ix1:ix2) = 2./(xc(ix1+1:ix2+1)-xc(ix1-1:ix2-1))
          av(1 ) = 2./(-3.*xc(1 )+4.*xc(2   )-xc(3   ))
          av(nx) = 2./( 3.*xc(nx)-4.*xc(nx-1)+xc(nx-2))
c
          aw = av
c
c.....generate uniform x grid

        else
c
          ix1=2
          ix2=nx-1

          xlen = xmax-xmin
          delx = xlen/real(nx-2)
          rdelx = 1./delx

          open(81,file='xgrid.dat')
          do i=1,nx
            xu(i)=real(i-1)*delx+xmin
            write(81,*)i,xu(i)
            ap(i)=rdelx
            au(i)=rdelx
            av(i)=rdelx
            aw(i)=rdelx
          enddo
          close(81)
c
        endif
c
c.....generate uniform y grid
c
        jy1 = 2
        jy2 = ny-1

        if(icyl==1) then
          ymin = 0.0
          ymax = 2.0*pi
        endif
        ylen = ymax-ymin        
        dely = ylen/real(ny-2)
        rdely = 1./dely

        open(81,file='ygrid.dat')
        do j=1,ny
          yv(j)=real(j-1)*dely+ymin
          write(81,*)j,yv(j)
          bp(j)=rdely
          bu(j)=rdely
          bv(j)=rdely
          bw(j)=rdely
        enddo
        close(81)
c
c.....generate non-uniform z grid
c.....the grid files contains nzg-1 coordinates: the nzg coordinate
c.....comes from the kz2 and kz2-1 coordinates
c
        if(kgrid==1) then
c
          kz1=2
          OPEN(UNIT=3,FILE=STR3,STATUS='OLD',FORM='FORMATTED')
          read(3,*) kz2
          read(3,*) (j,zwg(k),k=1,kz2)
          close(3)

          zmin = zwg(1)
          zmax = zwg(nzg)
          zlen=zwg(kz2)-zwg(1)

          zwg(nzg) = 2.*zwg(kz2)-zwg(kz2-1)
c     
          zcg(kz1:kz2) = .5*(zwg(kz1-1:kz2-1)+zwg(kz1:kz2))
          zcg(1  ) = 2.*zwg(1  )-zcg(2  )
          zcg(nzg) = 2.*zwg(kz2)-zcg(kz2)
c     
          cwg(kz1:kz2) = 2./(zwg(kz1+1:kz2+1)-zwg(kz1-1:kz2-1))
          cwg(1  ) = 2./(-3.*zwg(1  )+4.*zwg(2    )-zwg(3    ))
          cwg(kz2) = 2./( 3.*zwg(kz2)-4.*zwg(kz2-1)+zwg(kz2-2))
          cwg(nzg) = cw(kz2)
c
          cpg(kz1:kz2) = 1./(zwg(kz1:kz2)-zwg(kz1-1:kz2-1))
          cpg(1  ) = cpg(kz1)
          cpg(nzg) = cpg(kz2)
c
          cvg(kz1:kz2) = 2./(zcg(kz1+1:kz2+1)-zcg(kz1-1:kz2-1))
          cvg(1  ) = 2./(-3.*zcg(1  )+4.*zcg(2    )-zcg(3    ))
          cvg(nzg) = 2./( 3.*zcg(nzg)-4.*zcg(nzg-1)+zcg(nzg-2))
c
          cug = cvg
c
c.....generate uniform z grid
c
        else

          kz1=2
          kz2=nzg-1
c
          zlen = zmax-zmin
          delz = zlen/real(nzg-2)
          rdelz = 1./delz
c
          open(81,file='zgrid.dat')
          do k=1,nzg
            zwg(k)=real(k-1)*delz+zmin
            write(81,*)k,zwg(k)
            cpg(k)=rdelz
            cug(k)=rdelz
            cvg(k)=rdelz
            cwg(k)=rdelz
          enddo
          close(81)

        endif
c
c.....grid center
c
        do i=ix1,ix2
          xc(i)=0.5*(xu(i)+xu(i-1))
        enddo
        xc(ix1-1)=2.*xu(ix1-1)-xc(ix1)
        xc(ix2+1)=2.*xu(ix2  )-xc(ix2)
c     
        do j=jy1,jy2
          yc(j)=0.5*(yv(j)+yv(j-1))
        enddo
        yc(jy1-1)=2.*yv(jy1-1)-yc(jy1)
        yc(jy2+1)=2.*yv(jy2  )-yc(jy2)
c
        do k=kz1,kz2
          zcg(k)=0.5*(zwg(k)+zwg(k-1))
        enddo
        zcg(kz1-1)=2.*zwg(kz1-1)-zcg(kz1)
        zcg(kz2+1)=2.*zwg(kz2  )-zcg(kz2)
c
        if(IOGRID==1) then
          
        OPEN(UNIT=4,FILE=STR4,FORM='FORMATTED')
C
C.....indices in three directions
C
        write(4,*) ix1,ix2
        write(4,*) jy1,jy2
        write(4,*) kz1,kz2
        write(4,*)
c
        write(4,*) (xu(i),i=1,nx)
        write(4,*) (xc(i),i=1,nx)
        write(4,*) (ap(i),i=1,nx)
        write(4,*) (au(i),i=1,nx)
        write(4,*) (av(i),i=1,nx)
        write(4,*) (aw(i),i=1,nx)
        write(4,*) 
c
        write(4,*) (yv(j),j=1,ny)
        write(4,*) (yc(j),j=1,ny)
        write(4,*) (bp(j),j=1,ny)
        write(4,*) (bu(j),j=1,ny)
        write(4,*) (bv(j),j=1,ny)
        write(4,*) (bw(j),j=1,ny)
        write(4,*) 
c
        write(4,*) (zwg(k),k=1,nzg)
        write(4,*) (zcg(k),k=1,nzg)
        write(4,*) (cpg(k),k=1,nzg)
        write(4,*) (cug(k),k=1,nzg)
        write(4,*) (cvg(k),k=1,nzg)
        write(4,*) (cwg(k),k=1,nzg)
c
        CLOSE(4)

        endif
c
        kz2=(kz2-1)/mysize+1    !! Parallelization in the streamwise
c                                  direction
c
        if(IOGRID==1) then

        open(120,file='grid3dc.g',form='unformatted')
        write(120) nx,ny,nzg
        do k=1,nzg
          write(120) ((REAL(xc (i),4),i=1,nx),j=1,ny)
     &         ,     ((REAL(yc (j),4),i=1,nx),j=1,ny)
     &         ,     ((REAL(zcg(k),4),i=1,nx),j=1,ny)
        enddo
        close(120)

!        call hdf5_3Dgrid_sp('grid3du.h5sp',xu,yc,zcg,nx,ny,nzg)
!        call hdf5_3Dgrid_sp('grid3dv.h5sp',xc,yv,zcg,nx,ny,nzg)
!        call hdf5_3Dgrid_sp('grid3dw.h5sp',xc,yc,zwg,nx,ny,nzg)
!        call hdf5_3Dgrid_sp('grid3dc.h5sp',xc,yc,zcg,nx,ny,nzg)
        
c
c        ENDIF

c        IF(.TRUE.) then

        open(120,file='grid3du.g',form='unformatted')
        write(120) nx,ny,nzg
        do k=1,nzg
          write(120) ((REAL(xu (i),4),i=1,nx),j=1,ny)
     &         ,     ((REAL(yc (j),4),i=1,nx),j=1,ny)
     &         ,     ((REAL(zcg(k),4),i=1,nx),j=1,ny)
        enddo
        close(120)
c        ENDIF

c        IF(.TRUE.) THEN
c
        open(120,file='grid3dv.g',form='unformatted')
        write(120) nx,ny,nzg
        do k=1,nzg
          write(120) ((REAL(xc (i),4),i=1,nx),j=1,ny)
     &         ,     ((REAL(yv (j),4),i=1,nx),j=1,ny)
     &         ,     ((REAL(zcg(k),4),i=1,nx),j=1,ny)
        enddo
        close(120)
c
        open(120,file='grid3dw.g',form='unformatted')
        write(120) nx,ny,nzg
        do k=1,nzg
          write(120) ((REAL(xc (i),4),i=1,nx),j=1,ny)
     &         ,     ((REAL(yc (j),4),i=1,nx),j=1,ny)
     &         ,     ((REAL(zwg(k),4),i=1,nx),j=1,ny)
        enddo
        close(120)

        open(120,file='grid2dc.g',form='unformatted')
        write(120) nx,1,nzg
        do k=1,nzg
           write(120) ((REAL(xc (i),4),i=1,nx),j=1,1)
     &          ,     ((REAL(yc (j),4),i=1,nx),j=1,1)
     &          ,     ((REAL(zcg(k),4),i=1,nx),j=1,1)
        enddo
        close(120)

c        open(120,file='grid3d.g',form='unformatted')
c        write(120) nx,ny,nzg
c        do k=1,nzg
c          write(120) ((REAL(xu (i),4),i=1,nx),j=1,ny)
c     &         ,     ((REAL(yv (j),4),i=1,nx),j=1,ny)
c     &         ,     ((REAL(zwg(k),4),i=1,nx),j=1,ny)
c        enddo
c        close(120)


        if(icyl==1 .AND. .true.) then

          open(120,file='grid3dc_cyl.g',form='unformatted')
          write(120) nx,ny,nzg
          do k=1,nzg
            write(120) ((REAL(xc (i)*cos(yc(j)),4),i=1,nx),j=1,ny)
     &            ,    ((REAL(xc (i)*sin(yc(j)),4),i=1,nx),j=1,ny)
     &            ,    ((REAL(zcg(k),4),i=1,nx),j=1,ny)
          enddo
          close(120)
c
          open(120,file='grid3du_cyl.g',form='unformatted')
          write(120) nx,ny,nzg
          do k=1,nzg
            write(120) ((REAL(xu (i)*cos(yc(j)),4),i=1,nx),j=1,ny)
     &            ,    ((REAL(xu (i)*sin(yc(j)),4),i=1,nx),j=1,ny)
     &            ,    ((REAL(zcg(k),4),i=1,nx),j=1,ny)
          enddo
          close(120)
c
          open(120,file='grid3dv_cyl.g',form='unformatted')
          write(120) nx,ny,nzg
          do k=1,nzg
            write(120) ((REAL(xc (i)*cos(yv(j)),4),i=1,nx),j=1,ny)
     &            ,    ((REAL(xc (i)*sin(yv(j)),4),i=1,nx),j=1,ny)
     &            ,    ((REAL(zcg(k),4),i=1,nx),j=1,ny)
          enddo
          close(120)
c
          open(120,file='grid3dw_cyl.g',form='unformatted')
          write(120) nx,ny,nzg
          do k=1,nzg
            write(120) ((REAL(xc (i)*cos(yc(j)),4),i=1,nx),j=1,ny)
     &            ,    ((REAL(xc (i)*sin(yc(j)),4),i=1,nx),j=1,ny)
     &            ,    ((REAL(zwg(k),4),i=1,nx),j=1,ny)
          enddo
          close(120)

        endif

        IF(IT2D/=0) THEN
          J=NY/2
          open(120,file='grid3dc_'//INDEX(J)//'.g',form='unformatted')
          write(120) nx,1,nzg
          do k=1,nzg
            write(120) (REAL(xc (i),4),i=1,nx)
     &           ,     (REAL(yc (j),4),i=1,nx)
     &           ,     (REAL(zcg(k),4),i=1,nx)
          enddo
          close(120)
c
          open(120,file='grid3du_'//INDEX(J)//'.g',form='unformatted')
          write(120) nx,1,nzg
          do k=1,nzg
            write(120) (REAL(xu (i),4),i=1,nx)
     &           ,     (REAL(yc (j),4),i=1,nx)
     &           ,     (REAL(zcg(k),4),i=1,nx)
          enddo
          close(120)
c
          open(120,file='grid3dv_'//INDEX(J)//'.g',form='unformatted')
          write(120) nx,1,nzg
          do k=1,nzg
            write(120) (REAL(xc (i),4),i=1,nx)
     &           ,     (REAL(yv (j),4),i=1,nx)
     &           ,     (REAL(zcg(k),4),i=1,nx)
          enddo
          close(120)
c
          open(120,file='grid3dw_'//INDEX(J)//'.g',form='unformatted')
          write(120) nx,1,nzg
          do k=1,nzg
            write(120) (REAL(xc (i),4),i=1,nx)
     &           ,     (REAL(yc (j),4),i=1,nx)
     &           ,     (REAL(zwg(k),4),i=1,nx)
          enddo
          close(120)
        ENDIF

        endif

      ENDIF

      CALL MPI_BCAST(IX1,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(IX2,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(JY1,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(JY2,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(KZ1,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(KZ2,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)

      CALL MPI_BCAST(XLEN,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(YLEN,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ZLEN,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ZMIN,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ZMAX,1,MTYPE,0,MPI_COMM_EDDY,IERR)

      CALL MPI_BCAST(XU,NX,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(XC,NX,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(AP,NX,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(AU,NX,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(AV,NX,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(AW,NX,MTYPE,0,MPI_COMM_EDDY,IERR)

      CALL MPI_BCAST(YV,NY,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(YC,NY,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(BP,NY,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(BU,NY,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(BV,NY,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(BW,NY,MTYPE,0,MPI_COMM_EDDY,IERR)

      CALL MPI_BCAST(ZWG,NZG,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(ZCG,NZG,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(CPG,NZG,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(CUG,NZG,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(CVG,NZG,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(CWG,NZG,MTYPE,0,MPI_COMM_EDDY,IERR)

      CALL MPI_BCAST(DELX,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(DELY,1,MTYPE,0,MPI_COMM_EDDY,IERR)
      CALL MPI_BCAST(DELZ,1,MTYPE,0,MPI_COMM_EDDY,IERR)

      delxsq=delx*delx
      delysq=dely*dely
      delzsq=delz*delz
      
C     SET THE LOCAL METRIC
c
c the local values of the coordinates and of the coefficient
c are set
c
      DO K=KZ1-1,KZ2+1
         IERR = MYRANK*(NZ-2)+K
         ZW(K) = ZWG(IERR)
         ZC(K) = ZCG(IERR)
         CP(K) = CPG(IERR)
         CU(K) = CUG(IERR)
         CV(K) = CVG(IERR)
         CW(K) = CWG(IERR)
      ENDDO
         
      IF(MYRANK==0) THEN
         WRITE(6,'(A,3I5)')'*..local domain size   :',nx,ny,nz
         WRITE(6,'(A,3I5)')'*..global domain size  :',nx,ny,nzg
      endif
c
c the values of the staggered and centered radii are set
c also jsym and ranksym are set
c
      if(icyl==1) then
        ru(1:nx)=xu(1:nx)
c
        rp(1:nx)=xc(1:nx)
c
        do j=1,ny
          jsym(j)=j+(ny-2)/2
          if(jsym(j).gt.jy2) jsym(j)=jsym(j)-(ny-2)
        end do

c        do j=0,mysize-1
          ranksym = myrank+(mysize)/2
          if(ranksym.gt.mysize-1) ranksym = myrank-(mysize/2)
c        enddo
c        write(6,*) 'myrank=',myrank,',ranksym=',ranksym
c
      else
          ru(:)=1.
          rp(:)=1.
c
          do j=1,ny
            jsym(j)=j
          enddo

      endif
C 
      DEALLOCATE(cug,cvg)

      RETURN 
      END 



C-----function mindxdydz--------------------------N. Beratlis-21 Dec. 2008--
C
      REAL function mindxdydz(xc,yc,zc,nx,ny,nz,i,j,k,icyl)
C
C     PURPOSE: Find min dimension of cell i,j,k
C
C---------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER i,j,k,nx,ny,nz,icyl
      REAL    xc(nx),yc(ny),zc(nz)

      REAL    dx,dy,dz

      dy=abs(yc(j+1)-yc(j))
      IF(icyl==1) THEN
        dy = abs(xc(i)*dy)
      ENDIF
      dx = abs(xc(i+1)-xc(i))
      dz = abs(zc(k+1)-zc(k))

      mindxdydz = 100000.
      IF(dx<mindxdydz) mindxdydz=dx
      IF(dy<mindxdydz) mindxdydz=dy
      IF(dz<mindxdydz) mindxdydz=dz

      RETURN
      END
C---------------------------------------------------------------------------


C-----subroutine ijk_xyz--------------------------N. Beratlis-21 Dec. 2008--
C
      SUBROUTINE ijk_xyz(xp,yp,zp,xc,yc,zc,nx,ny,nz,ic,jc,kc,icyl)
C
C     PURPOSE: Find the low indices of the cell containing xp,yp,zp
C---------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER ic,jc,kc,nx,ny,nz,icyl
      REAL    xp,yp,zp
      REAL    xc(nx),yc(ny),zc(nz)

      REAL    rp,thetap

      REAL    anglerad

      IF(icyl==0) THEN
        CALL LOCATE(xc,nx,xp,ic)
        CALL LOCATE(yc,ny,yp,jc)
      ELSE
        rp=sqrt(xp**2.+yp**2.)
        thetap = anglerad(xp,yp)
        CALL LOCATE(xc,nx,rp,ic)
        CALL LOCATE(yc,ny,thetap,jc)
      ENDIF
      CALL LOCATE(zc,nz,zp,kc)

      RETURN

      END
C------------------------------------------------------------------------




C---- subroutine calc_zg-------------------------N. Beratlis-7 Jan. 2009-
C
C     PURPOSE: Calculate global z coordinate
C
C------------------------------------------------------------------------
      subroutine calc_zg(z,zg,nz,nzg)

      include 'common.h'
      include 'mpif.h'
      INTEGER nz,nzg
      REAL    z(nz),zg(nzg)

      INTEGER k
      REAL    zg1(nzg)

      if(mysize==1) then
        zg=z
      else
        zg1=0.0
        do k=2,nz-1
          zg1(k+(nz-2)*myrank)=z(k)
        enddo

        if(myrank==0) zg1(1)=z(1)
        if(myrank==mysize-1) zg1(nzg)=z(nz)
        CALL MPI_ALLREDUCE(zg1,zg,nzg,MTYPE,MPI_SUM,MPI_COMM_EDDY,IERR)
      endif

      return

      end
C------------------------------------------------------------------------


C---- SUBROUTINE CENTERVEL----------------------N. Beratlis--8 Apr. 2009-
C
C     PURPOSE: Center velocities U,V,W
C
C------------------------------------------------------------------------
C
      SUBROUTINE CENTERVEL(U,UC,NX,NY,NZ,DIR)

      include 'common.h'

      INTEGER NX,NY,NZ,DIR
      REAL    U(NX,NY,NZ),UC(NX,NY,NZ)

      INTEGER I,J,K

      IF(DIR==1) THEN

        DO I=2,NX
        DO J=1,NY
        DO K=1,NZ
          UC(i,j,k) = 0.5*(U(i,j,k)+U(i-1,j,k))
        ENDDO
        ENDDO
        ENDDO

        UC(1,:,:) = UC(2,:,:)

        if(icyl==1) then
          do j=1,ny
            uc(1,j,:) = -uc(2,jsym(j),:)     ! Treatment of the pizza
c                                            cell to ensure zero axis
c                                            velocity
          enddo
        endif

      ELSEIF(DIR==2) THEN

        DO I=1,NX
        DO J=2,NY
        DO K=1,NZ
          UC(i,j,k) = 0.5*(U(i,j,k)+U(i,j-1,k))
        ENDDO
        ENDDO
        ENDDO

        UC(:,1,:) = UC(:,NY-1,:)             ! Periodicity

        if(icyl==1) then
          do j=1,ny
            uc(1,j,:) = -uc(2,jsym(j),:)     ! Same as u
          enddo
        endif

      ELSEIF(DIR==3) THEN

        DO I=1,NX
        DO J=1,NY
        DO K=2,NZ
          UC(i,j,k) = 0.5*(U(i,j,k)+U(i,j,k-1))
        ENDDO
        ENDDO
        ENDDO

        UC(:,:,1) = UC(:,:,2)

        if(icyl==1) then
          do j=1,ny
            uc(1,j,:) = uc(2,jsym(j),:)   ! Ensure proper w at axis
          enddo
        endif

      ENDIF

      RETURN

      END
C-------------------------------------------------------------------------------


c---- function point_grid_size_ray------------------N. Beratlis-07 Dec. 2010 ---
C
C     PURPOSE: Estimate length of line formed by intersections of ray
C     with faces of cell. Ray passes through center of cell.
C
C-------------------------------------------------------------------------------
      real function point_grid_size_ray(xp,yp,zp,rnrm,xc,yc,zc,xc_car,yc_car
     $     ,nx,ny,nz,icyl)
c
      IMPLICIT NONE
c
c.... Input/Output arrays
      integer nx,ny,nz,icyl
      real    xp,yp,zp
      real    rnrm(3)
      real    xc(nx),yc(ny),zc(nz)
      real    xc_car(nx,ny),yc_car(nx,ny)
c
c.... Local arrays
      integer i,j,k
      integer iflag1,iflag2,dir1,dir2
      integer i1,j1,k1,i2,j2,k2
      real    xint1,yint1,zint1,xint2,yint2,zint2,dcell
      real    q(3),rvec(3)
c
c.... Functions
      integer body2grid

      call ijk_xyz(xp,yp,zp,xc,yc,zc,nx,ny,nz,i,j,k,icyl)
      if(icyl==0) then
        q(1) = 0.5*(xc(i)+xc(i+1))
        q(2) = 0.5*(yc(j)+yc(j+1))
      else
         q(1) = 0.5*(xc(i)+xc(i+1))*cos(0.5*(yc(j)+yc(j+1)))+1.0e-5
         q(2) = 0.5*(xc(i)+xc(i+1))*sin(0.5*(yc(j)+yc(j+1)))
       endif
       q(3) = 0.5*(zc(k)+zc(k+1))

       rvec = 1000.0*rnrm

       iflag1 = body2grid(q,rvec,xc_car,yc_car,xc,yc,zc
     &      ,nx,ny,nz,xint1,yint1,zint1,i1,j1,k1,dir1,icyl)

       iflag2 = body2grid(q,-rvec,xc_car,yc_car,xc,yc,zc
     &      ,nx,ny,nz,xint2,yint2,zint2,i2,j2,k2,dir2,icyl)

       if(iflag1==1 .AND. iflag2==1 .AND. dir1/=2 .AND. dir2/=2) then
         dcell = sqrt( (zint2-zint1)**2. + (yint2-yint1)**2. + (xint2-xint1)**2. )
       else
         call locate(zc,nz,zp,k)
         dcell = zc(k+1)-zc(k)
       endif

       point_grid_size_ray = dcell

       return

       end
C-------------------------------------------------------------------------------


c---- function ijk_grid_size_ray--------------------N. Beratlis-07 Dec. 2010 ---
C
C     PURPOSE: Estimate length of line formed by intersections of ray
C     with faces of cell. Ray passes through center of cell.
C
C-------------------------------------------------------------------------------
      real function ijk_grid_size_ray(i,j,k,rnrm,xc,yc,zc,xc_car,yc_car
     $     ,nx,ny,nz,icyl)
c
      IMPLICIT NONE
c
c.... Input/Output arrays
      integer nx,ny,nz,i,j,k,icyl
      real    rnrm(3)
      real    xc(nx),yc(ny),zc(nz)
      real    xc_car(nx,ny),yc_car(nx,ny)
c
c.... Local arrays
      integer iflag1,iflag2,dir1,dir2
      integer i1,j1,k1,i2,j2,k2
      real    xint1,yint1,zint1,xint2,yint2,zint2,dcell,a,b,dy
      real    q(3),rvec(3)
c
c.... Functions
      integer body2grid

c in Q the centered position is considered
      if(icyl==0) then
        q(1) = 0.5*(xc(i)+xc(i+1))
        q(2) = 0.5*(yc(j)+yc(j+1))
      else
         q(1) = 0.5*(xc(i)+xc(i+1))*cos(0.5*(yc(j)+yc(j+1)))+1.0e-5
         q(2) = 0.5*(xc(i)+xc(i+1))*sin(0.5*(yc(j)+yc(j+1)))
       endif
       q(3) = 0.5*(zc(k)+zc(k+1))

       rvec = 1000.0*rnrm
       !Rotate vector by half dtheta for icyl=1
       if(icyl==1) then
         dy = 0.5*(yc(j+1)-yc(j))
         a = rvec(1)
         b = rvec(2)
         rvec(1) = -b*sin(dy) + a*cos(dy)
         rvec(2) =  b*cos(dy) + a*sin(dy)
       endif

c from the point Q the intersections with the cell faces are found on
c both sides

       iflag1 = body2grid(q,rvec,xc_car,yc_car,xc,yc,zc
     &      ,nx,ny,nz,xint1,yint1,zint1,i1,j1,k1,dir1,icyl)

       iflag2 = body2grid(q,-rvec,xc_car,yc_car,xc,yc,zc
     &      ,nx,ny,nz,xint2,yint2,zint2,i2,j2,k2,dir2,icyl)

c if the intersections are found along the Y direction they are not
c used: in this case the axial dimension of the computational cell is
c taken 
c on the contary the distance between the intersection points is
c considered

!!!!!!       if(iflag1==1 .AND. iflag2==1 .AND. dir1/=2 .AND. dir2/=2) then
       if(iflag1==1 .AND. iflag2==1 .AND. icyl*dir1/=2 .AND. icyl*dir2/=2) then
         dcell = sqrt( (zint2-zint1)**2. + (yint2-yint1)**2. + (xint2-xint1)**2. )
       else
         dcell = zc(k+1)-zc(k)
       endif

       ijk_grid_size_ray = dcell

       return

       end
C-------------------------------------------------------------------------------
