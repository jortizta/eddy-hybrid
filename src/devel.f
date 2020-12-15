C---- subroutine taguo-------------------- N. Beratlis-22 May 2009 ---
C
C     PURPOSE: Tag outer velocity points
C
c---------------------------------------------------------------------
c flag: input; -1 interior points; 1 exterior points.
c flago: output; ibd interior points; 0 fluid points; -ibd interface points
      subroutine taguo(flag,flago,nx,ny,nz,nbd,ibd)
C
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer nx,ny,nz,ibd,nbd
      integer flag(nx,ny,nz),flago(nx,ny,nz,nbd)
c
c...  Local arrays
      integer i,j,k,i1,i2,j1,j2,i3,j3,k3,k1,k2,nn
      integer nbor(6)

      nn=0
c if itagy=0 the neighbour points along the Y direction are considered
c for tagging purposes; on the contrary only the ones along the X and Z
c directions are taken into account
      if(itagy==0) nn=2

      flago(:,:,:,ibd) = 0
c     
c.....define interior points
      do k=kbmin(ibd),kbmax(ibd)
      do j=jbmin(ibd),jbmax(ibd)
      do i=ibmin(ibd),ibmax(ibd)
        if(flag(i,j,k)==-1) flago(i,j,k,ibd)=ibd
      enddo
      enddo
      enddo
      flago(:,1,:,ibd) = flago(:,ny-1,:,ibd)
      flago(:,ny,:,ibd) = flago(:,2,:,ibd)

      i1 = ibmin(ibd)
      i2 = ibmax(ibd)
      j1 = jbmin(ibd)
      j2 = jbmax(ibd)
      k1 = max(2,kbmin(ibd))
      k2 = min(nz-1,kbmax(ibd))
      IF(icyl==1) THEN
        j1=2
        j2=ny-1
      ENDIF
      
c.....define outer boundary points
      do k=k1,k2
      do j=j1,j2
      do i=max(2,i1),i2
        nbor=1
        nbor(1)=flag(i,j,k)+flag(i+1,j,k)
        nbor(2)=flag(i,j,k)+flag(i-1,j,k)
        nbor(3)=flag(i,j,k)+flag(i,j,k+1)
        nbor(4)=flag(i,j,k)+flag(i,j,k-1)
        nbor(5)=flag(i,j,k)+flag(i,j-1,k)
        nbor(6)=flag(i,j,k)+flag(i,j+1,k)
          
        if(product(nbor(1:6-nn))==0) then
c.....outer boundary
          if(flag(i,j,k)==1) then
            flago(i,j,k,ibd)=-ibd
          endif    
        endif
c
      enddo
      enddo
      enddo

c      i=2
c      j=2
c      k = 434 - myrank*(nz-2)
c      if(k>0 .AND. k<=nz) then
c        write(6,*) 'taguo, myrank=',myrank,i,j,k,k+myrank*(nz-2),k1,k2,kbmin(ibd),kbmax(ibd),flag(i,j,k),flago(i,j,k,1)
c      endif
      
      if(icyl==1 .and. ibmin(ibd)==1) then

        i=1
        do k=k1,k2
        do j=j1,j2              !2,nym
c
          nbor(1)=flag(i,j,k)+flag(i+1,j,k)
          nbor(2)=flag(i,j,k)+flag(i,j,k+1)
          nbor(3)=flag(i,j,k)+flag(i,j,k-1)
          nbor(4)=flag(i,j,k)+flag(i,j+1,k)
          nbor(5)=flag(i,j,k)+flag(i,j-1,k)
c
          if(product(nbor(1:5-nn))==0) then
c.....outer boundary
            if(flag(i,j,k)==1) then
              flago(i,j,k,ibd)=-ibd
            endif    
          endif
c
        enddo
        enddo
      endif
c
c...  Periodic boundary conditions in Y-direction
      if(idomy==0) then
        flago(:,1,:,ibd) = flago(:,ny-1,:,ibd)
        flago(:,ny,:,ibd) = flago(:,2,:,ibd)
      endif
c
      call refreshflag(flago(1,1,1,ibd),nx,ny,nz)

      return
      end

c---------------------------------------------------------------------


C---- subroutine tagui---------------------N. Beratlis-22 May 2009---
C
C     PURPOSE: Tag inner velocity points for calculating derivatives.
C
c---------------------------------------------------------------------
      subroutine tagui(flag,flagui,flagvi,flagwi,flaguo,flagvo,flagwo,nx,ny,nz,nbd,ibd)
C
      include 'common.h'
      include 'immersed.h'

      integer nx,ny,nz,ibd,nbd
      integer flag(nx,ny,nz)
      integer flagui(nx,ny,nz,nbd),flagvi(nx,ny,nz,nbd)
     &       ,flagwi(nx,ny,nz,nbd),flaguo(nx,ny,nz,nbd)
     &       ,flagvo(nx,ny,nz,nbd),flagwo(nx,ny,nz,nbd)
c
c...  Local arrays
      integer i,j,k,k1,k2,i1,i2,j1,j2,nfe

      nfe=0
      flag = 0

      k1 = kbmin(ibd)
      k2 = kbmax(ibd)
      i1 = ibmin(ibd)
      i2 = ibmax(ibd)
      j1 = jbmin(ibd)
      j2 = jbmax(ibd)
      if(icyl==1) then
        j1 = 2
        j2 = ny-1
      endif
c     
c.....define true field extension points      
      do k=k1,k2
      do j=j1,j2
      do i=i1,i2
        if(flagui(i,j,k,ibd)/=0.AND.flaguo(i,j,k,ibd)==0) then

          if(flagui(i+1,j  ,k  ,ibd)==ibd) flag(i+1,j  ,k  )=-ibd
          if(flagui(i-1,j  ,k  ,ibd)==ibd) flag(i-1,j  ,k  )=-ibd
          if(flagui(i  ,j+1,k  ,ibd)==ibd) flag(i  ,j+1,k  )=-ibd
          if(flagui(i  ,j-1,k  ,ibd)==ibd) flag(i  ,j-1,k  )=-ibd
          if(flagui(i  ,j  ,k+1,ibd)==ibd) flag(i  ,j  ,k+1)=-ibd
          if(flagui(i  ,j  ,k-1,ibd)==ibd) flag(i  ,j  ,k-1)=-ibd

        endif

        if(flagvi(i,j,k,ibd)/=0.AND.flagvo(i,j,k,ibd)==0) then

          if(flagui(i  ,j  ,k  ,ibd)==ibd) flag(i  ,j  ,k  )=-ibd
          if(flagui(i  ,j+1,k  ,ibd)==ibd) flag(i  ,j+1,k  )=-ibd
          if(flagui(i-1,j+1,k  ,ibd)==ibd) flag(i-1,j+1,k  )=-ibd
          if(flagui(i-1,j  ,k  ,ibd)==ibd) flag(i-1,j  ,k  )=-ibd

        endif

        if(flagwi(i,j,k,ibd)/=0.AND.flagwo(i,j,k,ibd)==0) then

          if(flagui(i  ,j  ,k  ,ibd)==ibd) flag(i  ,j  ,k  )=-ibd
          if(flagui(i  ,j  ,k+1,ibd)==ibd) flag(i  ,j  ,k+1)=-ibd
          if(flagui(i-1,j  ,k+1,ibd)==ibd) flag(i-1,j  ,k+1)=-ibd
          if(flagui(i-1,j  ,k  ,ibd)==ibd) flag(i-1,j  ,k  )=-ibd

        endif

      enddo
      enddo
      enddo

      flag(:,1,:)=flag(:,ny-1,:)
      flag(:,ny,:) = flag(:,2,:)

      nfe = count(flag==-ibd)
c     
      return
      end

c---------------------------------------------------------------------


C---- subroutine tagvi---------------------N. Beratlis-22 May 2009---
C
C     PURPOSE: Tag inner velocity points for calculating derivatives.
C
c---------------------------------------------------------------------
      subroutine tagvi(flag,flagui,flagvi,flagwi,flaguo,flagvo,flagwo,nx,ny,nz,nbd,ibd)
C
      include 'common.h'
      include 'immersed.h'

      integer nx,ny,nz,ibd,nbd
      integer flag(nx,ny,nz)
      integer flagui(nx,ny,nz,nbd),flagvi(nx,ny,nz,nbd)
     &       ,flagwi(nx,ny,nz,nbd),flaguo(nx,ny,nz,nbd)
     &       ,flagvo(nx,ny,nz,nbd),flagwo(nx,ny,nz,nbd)
c
c...  Local arrays
      integer i,j,k,k1,k2,i1,i2,j1,j2,nfe

      nfe=0
      flag = 0

      k1 = kbmin(ibd)
      k2 = kbmax(ibd)
      i1 = ibmin(ibd)
      i2 = ibmax(ibd)
      j1 = jbmin(ibd)
      j2 = jbmax(ibd)
      if(icyl==1) then
        j1 = 2
        j2 = ny-1
      endif
c     
c.....define true field extension points      
      do k=k1,k2
      do j=j1,j2
      do i=i1,i2

        if(flagui(i,j,k,ibd)/=0.AND.flaguo(i,j,k,ibd)==0) then

          if(flagvi(i  ,j  ,k  ,ibd)==ibd) flag(i  ,j  ,k  )=-ibd
          if(flagvi(i+1,j  ,k  ,ibd)==ibd) flag(i+1,j  ,k  )=-ibd
          if(flagvi(i+1,j-1,k  ,ibd)==ibd) flag(i+1,j-1,k  )=-ibd
          if(flagvi(i  ,j-1,k  ,ibd)==ibd) flag(i  ,j-1,k  )=-ibd

        endif

        if(flagvi(i,j,k,ibd)/=0.AND.flagvo(i,j,k,ibd)==0) then

          if(flagvi(i+1,j  ,k  ,ibd)==ibd) flag(i+1,j  ,k  )=-ibd
          if(flagvi(i-1,j  ,k  ,ibd)==ibd) flag(i-1,j  ,k  )=-ibd
          if(flagvi(i  ,j+1,k  ,ibd)==ibd) flag(i  ,j+1,k  )=-ibd
          if(flagvi(i  ,j-1,k  ,ibd)==ibd) flag(i  ,j-1,k  )=-ibd
          if(flagvi(i  ,j  ,k+1,ibd)==ibd) flag(i  ,j  ,k+1)=-ibd
          if(flagvi(i  ,j  ,k-1,ibd)==ibd) flag(i  ,j  ,k-1)=-ibd

        endif

        if(flagwi(i,j,k,ibd)/=0.AND.flagwo(i,j,k,ibd)==0) then

          if(flagvi(i  ,j  ,k  ,ibd)==ibd) flag(i  ,j  ,k  )=-ibd
          if(flagvi(i  ,j  ,k+1,ibd)==ibd) flag(i  ,j  ,k+1)=-ibd
          if(flagvi(i  ,j-1,k+1,ibd)==ibd) flag(i  ,j-1,k+1)=-ibd
          if(flagvi(i  ,j-1,k  ,ibd)==ibd) flag(i  ,j-1,k  )=-ibd

        endif

      enddo
      enddo
      enddo

      flag(:,1,:)=flag(:,ny-1,:)
      flag(:,ny,:) = flag(:,2,:)

      nfe = count(flag==-ibd)
c     
      return
      end

c---------------------------------------------------------------------


C---- subroutine tagwi---------------------N. Beratlis-22 May 2009---
C
C     PURPOSE: Tag inner velocity points for calculating derivatives.
C
c---------------------------------------------------------------------
      subroutine tagwi(flag,flagui,flagvi,flagwi,flaguo,flagvo,flagwo,nx,ny,nz,nbd,ibd)
C
      include 'common.h'
      include 'immersed.h'

      integer nx,ny,nz,ibd,nbd
      integer flag(nx,ny,nz)
      integer flagui(nx,ny,nz,nbd),flagvi(nx,ny,nz,nbd)
     &       ,flagwi(nx,ny,nz,nbd),flaguo(nx,ny,nz,nbd)
     &       ,flagvo(nx,ny,nz,nbd),flagwo(nx,ny,nz,nbd)
c
c...  Local arrays
      integer i,j,k,k1,k2,i1,i2,j1,j2,nfe

      nfe=0
      flag = 0

      k1 = kbmin(ibd)
      k2 = kbmax(ibd)
      i1 = ibmin(ibd)
      i2 = ibmax(ibd)
      j1 = jbmin(ibd)
      j2 = jbmax(ibd)
      if(icyl==1) then
        j1 = 2
        j2 = ny-1
      endif
c     
c.....define true field extension points      
      do k=k1,k2
      do j=j1,j2
      do i=i1,i2

        if(flagui(i,j,k,ibd)/=0.AND.flaguo(i,j,k,ibd)==0) then

          if(flagwi(i  ,j  ,k  ,ibd)==ibd) flag(i  ,j  ,k  )=-ibd
          if(flagwi(i+1,j  ,k  ,ibd)==ibd) flag(i+1,j  ,k  )=-ibd
          if(flagwi(i+1,j  ,k-1,ibd)==ibd) flag(i+1,j  ,k-1)=-ibd
          if(flagwi(i  ,j  ,k-1,ibd)==ibd) flag(i  ,j  ,k-1)=-ibd

        endif

        if(flagvi(i,j,k,ibd)/=0.AND.flagvo(i,j,k,ibd)==0) then

          if(flagwi(i  ,j  ,k  ,ibd)==ibd) flag(i  ,j  ,k  )=-ibd
          if(flagwi(i  ,j+1,k  ,ibd)==ibd) flag(i  ,j+1,k  )=-ibd
          if(flagwi(i  ,j+1,k-1,ibd)==ibd) flag(i  ,j+1,k-1)=-ibd
          if(flagwi(i  ,j  ,k-1,ibd)==ibd) flag(i  ,j  ,k-1)=-ibd

        endif

        if(flagwi(i,j,k,ibd)/=0.AND.flagwo(i,j,k,ibd)==0) then

          if(flagwi(i+1,j  ,k  ,ibd)==ibd) flag(i+1,j  ,k  )=-ibd
          if(flagwi(i-1,j  ,k  ,ibd)==ibd) flag(i-1,j  ,k  )=-ibd
          if(flagwi(i  ,j+1,k  ,ibd)==ibd) flag(i  ,j+1,k  )=-ibd
          if(flagwi(i  ,j-1,k  ,ibd)==ibd) flag(i  ,j-1,k  )=-ibd
          if(flagwi(i  ,j  ,k+1,ibd)==ibd) flag(i  ,j  ,k+1)=-ibd
          if(flagwi(i  ,j  ,k-1,ibd)==ibd) flag(i  ,j  ,k-1)=-ibd

        endif

      enddo
      enddo
      enddo

      flag(:,1,:)=flag(:,ny-1,:)
      flag(:,ny,:) = flag(:,2,:)

      nfe = count(flag==-ibd)
c     
      return
      end

c---------------------------------------------------------------------


c---- subroutine tagpfe---------------------N. Beratlis-May 24 2009---
C
C     PURPOSE: Find true probelmatic pressure points. Array flapo contains
C     potential problematic points from previous time step.
C
c---------------------------------------------------------------------
      subroutine tagpfe(flagpo,flaguo,flagvo,flagwo,flag,nx,ny,nz,nbd,ibd)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c...input and output arrays
      integer nx,ny,nz,nbd,ibd
      integer flagpo(nx,ny,nz,nbd),flaguo(nx,ny,nz,nbd)
     &       ,flagvo(nx,ny,nz,nbd),flagwo(nx,ny,nz,nbd)
      integer flag(nx,ny,nz)

      integer i,j,k

      flag = 0

      do i=2,nx-1
      do j=2,ny-1
      do k=2,nx-1

        if(flaguo(i,j,k,ibd)==0) then
          if(flagpo(i+1,j,k,ibd)/=0) flag(i+1,j,k)=-ibd
          if(flagpo(i  ,j,k,ibd)/=0) flag(i  ,j,k)=-ibd
        endif

        if(flagvo(i,j,k,ibd)==0) then
          if(flagpo(i,j+1,k,ibd)/=0) flag(i,j+1,k)=-ibd
          if(flagpo(i,j  ,k,ibd)/=0) flag(i,j  ,k)=-ibd
        endif

        if(flagwo(i,j,k,ibd)==0) then
          if(flagpo(i,j,k+1,ibd)/=0) flag(i,j,k+1)=-ibd
          if(flagpo(i,j,k  ,ibd)/=0) flag(i,j,k  )=-ibd
        endif

      enddo
      enddo
      enddo

      flag(:,1 ,:) = flag(:,ny-1,:)
      flag(:,ny,:) = flag(:,2,:)

      flagpo(:,:,:,ibd) = flag(:,:,:)

      return

      end
c---------------------------------------------------------------------



C---- subroutine geom ----------------------N. Beratlis-May 22 2009---
C
C     PURPOSE: Find intersections and interpolation stencils for
C     the forcing nodes 
C     io = 0, velocity forcing points
C     io = 1, pressure and eddy viscosity forcing points
C
C---------------------------------------------------------------------
      subroutine geom(xu,yu,xu_car,yu_car,zu,vertex,vertexc,unvect,lim,mim
     &    ,iim,jim,kim,nim,xnim,ynim,znim,nxim,nyim,nzim,mrk,dir,fp
     &    ,flag,flagu,nx,ny,nz,nbd,nfacet,ibd,icom,nfu,nfumax,ivar,icycle,tlevel)     !!!!!! nfumaxg
C
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
      integer nx,ny,nz,nbd,nfacet,ibd,nfu,nfumax,nfumaxg,icycle,ivar
      integer lim(nfu),mim(nfu)
      integer iim(nfcmax),jim(nfcmax),kim(nfcmax)
      integer mrk(nfcmax),dir(nfcmax),fp(nfcmax)
      real    tlevel
      real    xu_car(nx,ny),yu_car(nx,ny)
      real    xu(nx),yu(ny),zu(nz)
      real    xnim(nfcmax),ynim(nfcmax),znim(nfcmax)
      real    nxim(nfcmax),nyim(nfcmax),nzim(nfcmax)
      integer flagu(nx,ny,nz,nbd),flag(nx,ny,nz)
      real    vertex(3,3,NFACET),vertexc(3,nfacet),unvect(3,nfacet)
      integer nim,icom,srtcrt
c
c....local variables and arrays
      integer iflagt,imin,intrs,gintrs,nitrs,n,nn,ibint,ilp,ip,intrs1,cintrs,io
      integer im,jm,km,i,j,k,ilb,ile,ic,ii,in,nimu,mimg,mimmax,mimmin
     &     ,mimu,mimug,nord,nordg,iord,mimprev,mimprevg,mimdif,ndif,itri
      integer i1,j1,k1,i2,icm,ncm
      real    xint,yint,zint,anglerad,q(3)
      integer icountnorm,icountgrid,icountnostenc,icountnointrs
      integer icountnormg,icountgridg,icountnostencg,icountnointrsg
      integer icountnrmnostenc,icountnrmnostencg,icountgrdnostenc,icountgrdnostencg,icountclonostenc,icountclonostencg
      integer icountnrmnointrs,icountnrmnointrsg,icountgrdnointrs,icountgrdnointrsg
      integer icountnovld(nfu),icountnovldg
      integer kmg,nq,inonvld
      integer nnrm,nnrmrchk,ngrd,ngrdrchk,nvld,nmrk
      integer nnrmg,nnrmrchkg,ngrdg,ngrdrchkg
      character*80 strio
 
c.... functions
      INTEGER forcpts,rdfn_forcpts,fluidface,ray_face_int,grid_intr
     $     ,norm_intr,recheckface,grid_intr_vel,closest_intr

      REAL    CLOCKTEMP,tclock
      INTEGER NCLOCKS
      REAL, DIMENSION(:), ALLOCATABLE :: CLOCK,CLOCKG,CLOCKGMAX,CLOCKGMIN

      INTEGER, DIMENSION(:), ALLOCATABLE :: ORD,ORDC,ORDN,NTRI,INN,NNTRI
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: INNF
      REAL, DIMENSION(:,:), ALLOCATABLE :: QUERYF

      nclocks=21
      ALLOCATE(CLOCK(NCLOCKS),CLOCKG(NCLOCKS),CLOCKGMAX(NCLOCKS),CLOCKGMIN(NCLOCKS))
      clock = 0.0
      clockg = 0.0
      clockgmin = 0.0
      clockgmax = 0.0

      clock(1) = tclock()

      ncm=3
      ALLOCATE(nntri(ncm))
      nntri(1) = 3
      nntri(2) = 5
      nntri(3) = nntr3D

      ALLOCATE(ntri(nntr3d))
      ntri = 0

      nimu = nim
      nfumax = 0
      nfumaxg = 0

      icountnovld = 0
      icountnovldg = 0
      icountgrid = 0
      icountgridg = 0
      icountnorm = 0
      icountnormg = 0
      icountnostenc = 0
      icountgrdnointrs = 0
      icountgrdnointrsg = 0
      icountnrmnostenc = 0
      icountnrmnostencg = 0
      icountnrmnointrs = 0
      icountnrmnointrsg = 0
      icountgrdnostenc = 0
      icountgrdnostencg = 0
      icountclonostenc = 0

      mimg = 0
      mimmax = 0
      mimmin = 0

      nnrm = 0
      nnrmg = 0
      ngrd = 0
      ngrdg = 0
      nnrmrchk = 0
      nnrmrchkg = 0
      ngrdrchk = 0
      ngrdrchkg = 0
      nord = 0
      nvld = 0
      nmrk = 0

      io = 0
      if(ivar>=4) io=1

      clock(2) = tclock()
      lim(1) = nim
      if(ivelrcn==0 .OR. io>0) then
c case for velocity forcing points when the GEOM_WM has been previously called
c and for non-velocity forcing points
        lim = 0
        mim = 0
        lim(1) = nim
c the number of the interface points is established in MIM(1)
c using the variable FLAGU (which is -IBD for the boundary points)
c in IIM,JIM,KIM the indices of the interface points are stored
        mim(1) = forcpts(flagu,nx,ny,nz,ibd,nbd,iim,jim,kim,nim)
        flag(:,:,:) = flagu(:,:,:,ibd)
      else
        lim(2:nfu) = 0
        mim(2:nfu) = 0
        call refreshflag(flagu(:,:,:,ibd),nx,ny,nz)    !!!!!!
        flag(:,:,:) = flagu(:,:,:,ibd)
        do i=1,lim(1)
          im = iim(i)
          jm = jim(i)
          km = kim(i)
          flag(im,jm,km) = 0
        enddo
        call refreshflag(flag(1,1,1),nx,ny,nz)
      endif
      nim = nimu + mim(1)

      inonvld = 0
      if(io==1) inonvld=ibd 

      if(nim>nfcmax) then
        write(6,'(A,1X,I9)') 'ERROR: Increase nfcmax to at least',nim
        call mpi_finalize(ierr)
        stop
      endif
      mrk(lim(1)+1:lim(1)+mim(1))=0
      clock(2) = tclock() - clock(2)

      nq = mim(1)
      allocate(innf(nq,nnTR3D+2),queryf(3,nq),ord(nq),ordc(nq),ordn(nq),inn(nntr3D))

      do i=1,nq
        ii = lim(1)+i
        queryf(1,i) = xu_car(iim(ii),jim(ii))
        queryf(2,i) = yu_car(iim(ii),jim(ii))
        queryf(3,i) = zu(kim(ii)) 
        ordc(i) = i+lim(1)
      enddo
     
      clock(3) = tclock()
!!!!!!      if(icom==0 .AND. nfacet>0) then
      if(icom<1 .AND. nfacet>0) then
c the subroutine COMMITTEE3 is initialized in the case ICOM is less than 1
        call committee3(QUERYF, VERTEXC, NFACET, INNF, NN, 0, MYRANK, imb3D)
        if(icom==-1) icom=2   !!!!!!
c at the end of the subroutine the variables of COMMITTEE3 are deallocated
      endif
      clock(3) = tclock() - clock(3)

      clock(5) = tclock()
      
      mimprev = 0
      ndif=0

c NFU: available number of loops for the search of valid interface points
      do ilp=1,nfu

c MIM: number of interface points for which to find a valid stencil
c MIM takes into account the points for which the stencil at the previous loop
c was composed of at least one interface point (the points for which no
c intersection was found or for which the stencil involved one body point
c were discarded at the previous loops)
        mimu = mim(ilp)
        call mpi_allreduce(mimu,mimmax,1,mpi_integer,mpi_max,mpi_comm_eddy,ierr)
        call mpi_allreduce(mimu,mimug,1,mpi_integer,mpi_sum,mpi_comm_eddy,ierr)
        call mpi_allreduce(mimprev,mimprevg,1,mpi_integer,mpi_sum,mpi_comm_eddy,ierr)
        
        mimdif = mimug-mimprevg
        mimprev = mimu
        if(mimdif==0) ndif=ndif+1

        if(ilp>1) then
          lim(ilp)=nimu + sum(mim(1:ilp-1))
        endif
        nmrk = 0

        if(mimmax==0) exit  ! no left interface points

        nfumax = nfumax+1  ! max number of loops to find valid interface points

        ii=lim(ilp)+1
        do while(ii<=lim(ilp)+mimu)

          i = ordc(ii-lim(ilp))

          im = iim(i)
          jm = jim(i)
          km = kim(i)

          if(mrk(i)==-1) then
c at the previous loop the external stencil along the normal (norm_intr) had
c at least 1 interface point
            !Check if face is intersecting fluid
c  1: all face points are fluid
c -1: at least one face point is an interface point
c -2: at least one face point is an interior point
            clocktemp = tclock()
            intrs = recheckface(flag,xu,yu,zu,nxim(i),nyim(i),nzim(i),nx,ny,nz,dir(i),icyl,idomy) 
            nnrmrchk = nnrmrchk+1
            clock(16) = clock(16) + tclock() - clocktemp
          elseif(mrk(i)==0) then

            icm = 1
            do while(icm<=ncm)
              nn = nntri(icm)
              innf(i-lim(1),nntr3d+1) = nn
              innf(i-lim(1),nntr3d+2) = icm
              q(:) = queryf(:,i-lim(1))
c committee3 finds the points of the array VERTEXC closest to Q considering the 3D distance
              call committee3(q, vertexc, nfacet, inn, nn, 1, myrank, imb3D)
              innf(i-lim(1),1:nn) = inn(1:nn)
              clocktemp = tclock()
c intrs=1 if a physical stencil is found; the output of committee3 is used
c in xnim, ynim and znim the coordinates of the intersection with
c the immersed-boundary are written
c in nxim, nyim and nzim the coordinates of the intersection with a grid face
c are provided by the function NORM_INTR
c fp is the index of the triangle where the intersection is found
c dir is the direction along which the intersection is found
c intrs=0: no intersection with the immersed boundary
c intrs=-1: the stencil contains interface points
c intrs=-2: the stencil contains body points
              intrs = norm_intr(flag,nx,ny,nz,ibd,nbd,vertex,unvect,nfacet
     &             ,xu_car,yu_car,zu,innf(i-lim(1),:),nn,xnim(i),ynim(i),znim(i),nxim(i),nyim(i)
     &             ,nzim(i),fp(i),im,jm,km,dir(i),1,itri)
              clock(12) = clock(12) + tclock() - clocktemp
              if(intrs/=0) exit   ! an intersection with the immersed-boundary has been found
              icm = icm+1
            enddo
            if(itri>0) ntri(itri) = ntri(itri)+1
            if(intrs==0) then
c no intersections found along the normal direction: fp is set equal to -1 for that interface point
              icountnrmnointrs = icountnrmnointrs+1
              fp(i)=-1
            endif
            nnrm = nnrm+1
          else
            intrs = 0
          endif
         
          if(im<=2) intrs=0

          if(intrs==1) then
c the stencil along the direction normal to the immersed-boundary is OK
            mrk(i)=1
            icountnorm = icountnorm+1
            nvld = nvld+1
            ord(nvld)=i
          elseif(intrs==-1 .AND. ilp<nfu .AND. mimdif/=0) then
c the stencil has at least one interface point: it is not valid
c another attempt is planned for the following loop
c the number of interface points at the current loop is decreased
c the number of interface points at the next loop is increased
            mrk(i) =-1
            nmrk = nmrk+1
            ordn(nmrk) = i
            mim(ilp) = mim(ilp)-1
            mim(ilp+1)=mim(ilp+1)+1
          else

            if(mrk(i)==-2 .AND. ndif<3) then
c if at the previous loop the exterior point found by GRID_INTR was an
c interface point (MRK(I)=-2)
              clocktemp = tclock()
              if(flag(im+int(nxim(i)),jm+int(nyim(i)),km+int(nzim(i)))==0) then
c the point along the outward grid direction is fluid: it can be used
                gintrs = 1
              else
c the point along the outward grid direction is a boundary point: it is not good
                gintrs = -1
              endif
              clock(17) = clock(17) + tclock() - clocktemp
              ngrdrchk = ngrdrchk+1
            else!if(mrk(i)/=-3) then
              clocktemp = tclock()
              srtcrt = 0
              if(fp(i)>0) then
                srtcrt=2
                if(ndif<2) srtcrt=1
              endif
                nn = innf(i-lim(1),nntr3d+1)
                icm = innf(i-lim(1),nntr3d+2)
                do while (icm<=ncm)
                  clocktemp = tclock()
c
c the intersection and the interpolation stencil are sought along the grid directions
c  1: the intersection is found and the stencil is physical;
c  0: the intersection is not found;
c -1: the intersection is found, but the projection point is boundary
c -2: the intersection is found, but the projection point is interior
c also in this case xnim, ynim and znim are the coordinates of the intersection
c point with the immersed-boundary, but nxim, nyim and nzim are the extensions
c of the indices of the exterior point (not the coordinates!)
c
                  gintrs = grid_intr(flag,nx,ny,nz,ibd,vertex,unvect,nfacet
     &                 ,xu_car,yu_car,zu,innf(i-lim(1),:),nn,xnim(i),ynim(i),znim(i)
     &                 ,nxim(i),nyim(i),nzim(i),im,jm,km,srtcrt,fp(i),io,itri)
                  clock(13) = clock(13) + tclock() - clocktemp

                  if(gintrs/=0 .AND. gintrs/=-2) exit
c
c if the intersection with the immersed-boundary has not been found or if the projection point
c is a body point: the search is performed on a larger number of triangles
c
                  icm = icm+1
                  if(icm<=ncm) then
                    nn = nntri(icm)
                    innf(i-lim(1),nntr3d+1) = nn
                    innf(i-lim(1),nntr3d+2) = icm
                    q(:) = queryf(:,i-lim(1))
                    call committee3(q, vertexc, nfacet, inn, nn, 1, myrank, imb3D)
                    innf(i-lim(1),1:nn) = inn(1:nn)
                  endif
                enddo
                if(itri>0) ntri(itri) = ntri(itri)+1
                if(mrk(i)==0 .AND. gintrs==0) icountgrdnointrs = icountgrdnointrs+1
c no intersections with the surface of the immersed-boundary along the grid directions
c
                ngrd = ngrd+1
c            else
c              gintrs = 0
            endif

            if(gintrs==1) then
c the exterior point along the grid direction is a fluid point: it is OK
              mrk(i) = 2
              nvld = nvld+1
              ord(nvld) = i
            elseif(gintrs==-1 .AND. ilp<nfu) then !stencil is forcing point
c the exterior point is a boundary point
c another attempt is planned for the following loop
c also in this case the number of interface points at the current loop
c is decreased and the one at the next loop is increased
              mrk(i) = -2
              nmrk = nmrk+1
              ordn(nmrk) = i
              mim(ilp) = mim(ilp)-1
              mim(ilp+1)=mim(ilp+1)+1
            else            !stencil is body point or no intrs. with body was found
              if(io==3) then
c the number of interface points is decreased: the current interface point is not valid
c this part seems not working in the current version
                write(6,*) 'point not used',ilp,i,mrk(i),im,jm,km+myrank*(nz-2),intrs,gintrs,io
                mrk(i) = 0
                flagu(im,jm,km,ibd) = inonvld
                icountnovld(ilp) = icountnovld(ilp)+1
                mim(ilp) = mim(ilp)-1

              else
                if(mrk(i)==-3) then
               !Check if face is intersecting fluid
c this check is associated with a previous call of the function CLOSEST_INTR
c  1: all face points are fluid
c -1: at least one face point is an interface point
c -2: at least one face point is an interior point
                   clocktemp = tclock()
                   cintrs = recheckface(flag,xu,yu,zu,nxim(i),nyim(i),nzim(i),nx,ny,nz,dir(i),icyl,idomy) 
                   nnrmrchk = nnrmrchk+1
                   clock(16) = clock(16) + tclock() - clocktemp
                else
                  clocktemp = tclock()
c
c the intersection of a forcing point with the immersed-boundary is found using the node closest to
c the surface of the body
c -2: the fluid face contains at least 1 body point
c -1: the fluid face contains at least 1 interface point
c  1: the fluid face contains only fluid points and therefore it is valid 
c
                  cintrs = closest_intr(flag,nx,ny,nz,ibd,nbd,vertex,unvect,nfacet
     &                 ,xu_car,yu_car,zu,innf(i-lim(1),:),innf(i-lim(1),nntr3d+1),xnim(i),ynim(i),znim(i),nxim(i),nyim(i)
     &                 ,nzim(i),fp(i),im,jm,km,dir(i),1)
c                  write(6,*) 'geom, closest_intr:',ilp,i,mrk(i),im,jm,km,intrs,gintrs,cintrs,io
                  clock(12) = clock(12) + tclock() - clocktemp
                  nnrm = nnrm+1
                endif

                if(cintrs==1) then
c the interface point is OK
                  mrk(i)=3
                  icountnorm = icountnorm+1
                  nvld = nvld+1
                  ord(nvld)=i
                elseif(cintrs==-1 .AND. ilp<nfu .AND. ndif<2) then
c the stencil is not OK: it involves at least 1 interface point
c another attempt is planned for the following loop
c also in this case the number of interface points at the current loop
c is decreased and the one at the next loop is increased
                  mrk(i) =-3
                  nmrk = nmrk+1
                  ordn(nmrk) = i
                  mim(ilp) = mim(ilp)-1
                  mim(ilp+1)=mim(ilp+1)+1
                else
c the stencil involves at least 1 body point: the interface point is discarded
                  mrk(i) = 0
                  flagu(im,jm,km,ibd) = inonvld
                  icountnovld(ilp) = icountnovld(ilp)+1
                  mim(ilp) = mim(ilp)-1
                  icountclonostenc = icountclonostenc+1
                  write(6,'(A,20(1x,I7))') 'geom, not valid:',ibd,ivar,ilp,i,mrk(i),im,jm,km,km+myrank*(nz-2),intrs,gintrs,cintrs
     &                 ,io,mimdif,ndif,innf(i-lim(1),nntr3d+1)
                endif
              endif
            endif
          endif

          ii = ii+1

        enddo

        clocktemp = tclock()
c
c flag is set to 0 in the points where mrk>0
c
        call physical_mrk2flag1(flag,nx,ny,nz,mrk,iim,jim,kim,ord,nq,lim(ilp)-nimu,mim(ilp))
        clock(18) = clock(18) + tclock() - clocktemp

        clocktemp = tclock()
        if(mimu>0) then
          if(idomy==0) then
c
c update of flag at the ghost layers along the Y direction
c
            call yrefresh3darray(flag,nx,ny,nz)
          endif
        endif

        if(mimmax>0) then
c.....Axis
c
c update of flag at the ghost layers at the axis
c
          if(itype(1)==300) call refreshcntln(flag,nx,ny,nz)
c
c update of flag at the ghost layers between processors
c
          call refreshflag(flag(1,1,1),nx,ny,nz)
        endif
        clock(14) = clock(14) + tclock() - clocktemp

        ordc = ordn
   
C.....redefine boundary points
c.....the number of interface points is updated in the case icountnovld>0
c     (some forcing points were discarded)
        clocktemp = tclock()
        IF(icountnovld(ilp)>0) THEN
          nim = nimu + sum(mim)
        ENDIF
        clock(15) = clock(15) + tclock() - clocktemp
c        if(ndif==2) exit

      enddo

      nim = nimu + sum(mim)
      CALL MPI_ALLREDUCE(nfumax,nfumaxg,1,MPI_INTEGER,MPI_MAX,MPI_COMM_EDDY,IERR)

 514  continue
      
      clocktemp = tclock()
c
c the vectors of the valid interface points are ordered
c ORD is a vector associated with the points for which INTRS, GINTRS or CINTRS
c are equal to 1 (valid boundary points)
c
      call rord_forcpts1(iim,jim,kim,xnim,ynim,znim,nxim,nyim,nzim
     &     ,mrk,dir,fp,ord,nq,lim(1),sum(mim))

      clock(21) = clock(21) + tclock() - clocktemp

      clock(5) = tclock() - clock(5)

      clock(8) = tclock()
      if(icom==2 .AND. nfacet>0) then
c the variables of COMMITTEE3 are deallocated
        CALL COMMITTEE3(QUERYF, VERTEXC, NFACET, INNF, NN, 2, MYRANK, imb3D)
      endif
      clock(8) = tclock() - clock(8)

      deallocate(queryf,innf,ord,ordc,ordn,inn)
c
c number of forcing points for which a valid stencil is found along the
c normal direction
      icountnorm = count(mrk(lim(1)+1:lim(1)+sum(mim))==1)
c
c number of forcing points for which a valid stencil is found along a
c grid direction
      icountgrid = count(mrk(lim(1)+1:lim(1)+sum(mim))==2)
c
c number of forcing points for which the stencil along the normal direction
c is not valid
      icountnrmnostenc = nq-icountnorm-icountnrmnointrs
c
c number of forcing points for which the stencil along a grid direction
c is not valid
      icountgrdnostenc = nq-icountnorm-icountgrid-icountgrdnointrs

      clocktemp = tclock()
c
c reduce to 0 process the statistics of the tagging procedure
c
      CALL MPI_REDUCE(SUM(ICOUNTNOVLD),ICOUNTNOVLDG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTNORM,ICOUNTNORMG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTGRID,ICOUNTGRIDG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(MIM),MIMG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(MIM),MIMMAX,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(MIM),MIMMIN,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTGRDNOSTENC,ICOUNTGRDNOSTENCG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTNRMNOINTRS,ICOUNTNRMNOINTRSG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTNRMNOSTENC,ICOUNTNRMNOSTENCG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTGRDNOINTRS,ICOUNTGRDNOINTRSG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTCLONOSTENC,ICOUNTCLONOSTENCG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(nnrm,nnrmg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ngrd,ngrdg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(nnrmrchk,nnrmrchkg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ngrdrchk,ngrdrchkg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(nord,nordg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      clock(9) = tclock() - clocktemp

      clocktemp = tclock()
      strio = 'forcing points'
      if(ivar==1) then 
        strio = 'U velocity '//strio
      elseif(ivar==2) then 
        strio = 'V velocity '//strio
      elseif(ivar==3) then 
        strio = 'W velocity '//strio
      elseif(ivar==4) then 
        strio = 'Pressure '//strio
      elseif(ivar==5) then 
        strio = 'viscosity '//strio
      endif
c
c the 0 process writes the statistics on a file
c
      IF(MYRANK.EQ.0 .AND. iolvl>0) THEN
        OPEN(UNIT=16, FILE='stats_imb.dat', FORM='FORMATTED'
     &        ,POSITION='APPEND')
        write(16,'(A,I6,A,E16.8,1x,2(A,1x,I9),A,F6.2,A,A,I8,A,F6.2,A)')
     &       'Cycle no=',icycle,', time=',tlevel
     &       ,trim(strio),MIMG
     &       ,', forced along normal:',ICOUNTNORMG
     &       ,' (',100.*REAL(ICOUNTNORMG,8)/REAL(MIMG,8),'%)'
     &       ,', forced along grid points:',ICOUNTGRIDG  
     &       ,' (',100.*REAL(ICOUNTGRIDG,8)/REAL(MIMG,8),'%)'
        write(16,'(A,F14.2,3(A,I9))') 'No. of forc. pts., ave:',real(mimg,8)/real(mysize,8)
     &       ,' ,max:',mimmax,', min:',mimmin,', Max no. of loops:',nfumaxg
        write(16,'(2(A,I9))') 'No. nrm. intrs calls=',nnrmg
     &       ,',no grd intrs calls=',ngrdg
        write(16,'(2(A,I9))') 'No. nrm. intrs rechecks=',nnrmrchkg
     &       ,',no grd intrs rechecks=',ngrdrchkg
c        write(16,'(A)') 'Triangle intersection stats:'
c        do i=1,nn
c          write(16,'(I4,1x,I8)') i,ntri(i)
c        enddo
        IF(ICOUNTNOVLDG>0) THEN
          WRITE(16,'(2(A,1X,I9,1X),A)')
     &          'WARNING:',ICOUNTNOVLDG,' of ',MIMG+ICOUNTNOVLDG
     &          ,trim(strio)//' forcing points can not be used:'
          write(16,'(A,2(I9,A))') 'along normal:'
     &          ,icountnrmnostencg,' due to no stencil,'
     &          ,icountnrmnointrsg,' due to no intrs.'
          write(16,'(A,2(I9,A))') 'along gridlines:'
     &          ,icountgrdnostencg,' due to no stencil,'
     &          ,icountgrdnointrsg,' due to no intrs.'
          write(16,'(A,1(I9,A))') 'along closest:'
     &          ,icountclonostencg,' due to no stencil.'
        ENDIF
        CLOSE(16)
      ENDIF
      clock(10) = tclock()-clocktemp

      clock(1) = tclock() - clock(1)

      IF(ioclock>0) THEN
c
c the 0 process writes on a file the computational times for the tagging procedure
c
        CALL MPI_REDUCE(CLOCK,CLOCKG,NCLOCKS,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(CLOCK,CLOCKGMIN,NCLOCKS,MTYPE,MPI_MIN,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(CLOCK,CLOCKGMAX,NCLOCKS,MTYPE,MPI_MAX,0,MPI_COMM_EDDY,IERR)
        clockg = clockg/real(mysize)

        where(clockg==0.0) clockg=1.e-8
        where(clockgmin==0.0) clockgmin=1.e-8
        where(clockgmax==0.0) clockgmax=1.e-8

        IF(MYRANK==0) THEN
          OPEN(UNIT=16,FILE='clock.dat',FORM='FORMATTED'
     &          ,POSITION='APPEND')
          write(16,'(A,I2,A)') '---- geom: icom=',icom,
     &       '---------------------------------------------------------'
          write(16,'(A)') '    Task/Time           Ave. (sec/%)        Max. (sec/%)   
     &     Min (sec/%)'
          write(16,'(A,3(6x,F8.4))') 'Total               :',clockg(1),clockgmax(1),clockgmin(1)
          i=2
          write(16,905) 'Bndpts points       :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=3
          write(16,905) 'Initialize ANN      :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=11
          write(16,905) 'Find near. neighbor.:'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=5
          write(16,905) 'Find intrs./stencil :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=12
          write(16,905) '   -Norm. intrs.    :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          write(16,'(A,I6,A)') '    (No. of calls=',nnrmg,')'
          i=16
          write(16,905) '   -Norm. intrs.rchk:'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          write(16,'(A,I6,A)') '    (No. of calls=',nnrmrchkg,')'
          i=13
          write(16,905) '   -Grid intrs.     :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          write(16,'(A,I6,A)') '    (No. of calls=',ngrdg,')'
          i=17
          write(16,905) '   -Grid intrs.rchk :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          write(16,'(A,I6,A)') '    (No. of calls=',ngrdrchkg,')'
          i=19
          write(16,905) '   -Reorder ord arr.:'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          write(16,'(A,I6,A)') '    (No. of calls=',nordg,')'
          i=20
          write(16,905) '   -Reorder innf    :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          i=14
          write(16,905) '   -Refresh flag    :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          i=18
          write(16,905) '   -Mrk2flag        :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          i=15
          write(16,905) '   -Redefine bndpt  :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          i=21
          write(16,905) 'Reorder forc. pts   :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=8
          write(16,905) 'Delete ANN structure:'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=9
          write(16,905) 'MPI reduce opertns. :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=10
          write(16,905) 'I/O opertns.        :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          write(16,'(2A)') '-----------------------------------------',
     &         '-----'
          close(16)
        ENDIF

      ENDIF

      deallocate(ntri,nntri)
      deallocate(clock,clockg,clockgmin,clockgmax)

      return

 905  format(A,3(1x,F8.4,A,F6.2,A,2x))

      end
c-----------------------------------------------------------------------



c--- subroutine mtrx -------------------------N. Beratlis-23 May 2009---
c
c    PURPOSE: Build interpolation stencil for forcing points
c    - imn = 0 : pressure forcing points
c    - imn = 1 : outer velocity points
c    - imn = 2 : inner velocity points
c
c-----------------------------------------------------------------------
      subroutine mtrx(xu,yu,zu,lim,mim,mrk,iim,jim,kim,xnim,ynim,znim,
     &     nxim,nyim,nzim,vmtrx,vindx,nx,ny,nz,nbd,ibd,imn)
C
      include 'common.h'
      include 'immersed.h'
      INCLUDE 'mpif.h'
c
      integer nx,ny,nz,nbd,ibd,imn
      integer lim,mim
      integer mrk(nfcmax),dir(nfcmax)
      integer iim(nfcmax),jim(nfcmax),kim(nfcmax)
      integer vindx(nsp,nfcmax)
      real    vmtrx(nsp,nsp,nfcmax)
      real    xnim(nfcmax),ynim(nfcmax),znim(nfcmax)
      real    nxim(nfcmax),nyim(nfcmax),nzim(nfcmax)
      real    xu(nx),yu(ny),zu(nz)
c
c....local variables and arrays
      integer i,ii,jj,kk,in,jn,kn,in1,jn1,kn1,nyg
      integer isign1,isign2
      integer indx(nsp)
      real    amtrx(nsp,nsp),d,so,s1
      real    xucar,yucar,zucar,xu1car,yu1car,zu1car,dtp

c
c...loop over outer boundary points  
      if(imn==0)then
        isign1 = 1
        isign2 = 1
      elseif(imn==1) then
        isign1 = -1
        isign2 = 1
      elseif(imn==2) then
        isign1 = 1
        isign2 = 1
      endif

      do i=lim+1,lim+mim
c
c Indices of the interface points
c
        in = iim(i)
        jn = jim(i)
        kn = kim(i)

        ii = int(nxim(i))
        jj = int(nyim(i))
        kk = int(nzim(i))
c
c This subroutine works always using Cartesian coordinates
c
        if(icyl==1) then
c
c Cartesian coordinates of the interface points
c
          xucar = xu(in)*cos(yu(jn))
          yucar = xu(in)*sin(yu(jn))
          zucar = zu(kn)
c
c If mrk=1 nxim, nyim and nzim are the Cartesian coordinates
c of the intersection points with a fluid face, while if mrk=2
c they are the relevant incremental indices
c
          if(mrk(i)==2) then
            xu1car = xu(in+ii)*cos(yu(jn+jj))
            yu1car = xu(in+ii)*sin(yu(jn+jj))
            zu1car = zu(kn+kk)
          else
            xu1car = nxim(i)
            yu1car = nyim(i)
            zu1car = nzim(i)
          endif

        else
          xucar = xu(in)
          yucar = yu(jn)
          zucar = zu(kn)
          
          if(mrk(i)==2) then
            xu1car = xu(in+ii)
            yu1car = yu(jn+jj)
            zu1car = zu(kn+kk)
          else
            xu1car = nxim(i)
            yu1car = nyim(i)
            zu1car = nzim(i)
          endif

        endif
c
c so: distance between the interface point and the immersed-boundary
c s1: distance between the interface point and the fluid point
c N.B. the origin of the local coordinates system is at the
c interface point
c
        so = sqrt( (xucar-xnim(i))**2.
     &       + (yucar-ynim(i))**2.
     &       + (zucar-znim(i))**2. )*real(isign1)
        s1 = sqrt( (xucar-xu1car)**2.
     &       + (yucar-yu1car)**2.
     &       + (zucar-zu1car)**2. )*real(isign2)
c.....matrix of coefficients (amtrx)
c.....imn/=0 is for the velocity interpolation, imn=0 is for the pressure one
        if(imn/=0) then
          amtrx(1,1)=1.
          amtrx(1,2:nsp)=so
        elseif(imn==0) then
          amtrx(1,1)=0.
          amtrx(1,2:nsp)=1.0
        endif
        amtrx(2:nsp,1)=1.
        amtrx(2:nsp,2)=s1

c        if(in==12 .AND. jn==5 .AND. kn==175) then
c          write(6,*) '1. mtrx',i,mrk(i),in,jn,jn,amtrx,indx,d,so,xnim(i),ynim(i),znim(i),xucar,yucar,zucar
c        endif
c
c LU decomposition of amtrx, whose output is stored in vmtrx and vindx
c
        call ludcmp(amtrx,nsp,nsp,indx,d)

c        if(in==12 .AND. jn==5 .AND. kn==175) then
c          write(6,*) '2. mtrx',i,mrk(i),in,jn,jn,amtrx,indx,d
c        endif

c.....store arrays
        vmtrx(:,:,i)=amtrx(:,:)
        vindx(:,i)=indx(:)

      enddo
    
      return
      end
C-------------------------------------------------------------------------





C---- subroutine calcdpdx-------------------N. Beratlis-15 May 2009---
C
C     PURPOSE: Calculate dpdx
C
C---------------------------------------------------------------------
      subroutine calcdpdx(P,DP,NX,NY,NZ)

      include 'common.h'

      integer nx,ny,nz
      real    p(nx,ny,nz),dp(nx,ny,nz)

      integer i,j,k
      
      dp = 0.0

      do i=1,nx-1
      do j=1,ny
      do k=1,nz
        dp(i,j,k) = au(i)*(p(i+1,j,k)-p(i,j,k))
      enddo
      enddo
      enddo

      return

      end
c---------------------------------------------------------------------


C---- subroutine calcdpdy-------------------N. Beratlis-15 May 2009---
C
C     PURPOSE: Calculate dpdy
C
C---------------------------------------------------------------------
      subroutine calcdpdy(P,DP,NX,NY,NZ)

      include 'common.h'

      integer nx,ny,nz
      real    p(nx,ny,nz),dp(nx,ny,nz)

      integer i,j,k
      
      dp = 0.0

      do i=1,nx
      do j=1,ny-1
      do k=1,nz
        dp(i,j,k) = bv(j)*(p(i,j+1,k)-p(i,j,k))
      enddo
      enddo
      enddo

      return

      end
c---------------------------------------------------------------------


C---- subroutine calcdpdz-------------------N. Beratlis-15 May 2009---
C
C     PURPOSE: Calculate dpdz
C
C---------------------------------------------------------------------
      subroutine calcdpdz(P,DP,NX,NY,NZ)

      include 'common.h'

      integer nx,ny,nz
      real    p(nx,ny,nz),dp(nx,ny,nz)

      integer i,j,k
      
      dp = 0.0

      do i=1,nx
      do j=1,ny
      do k=1,nz-1
        dp(i,j,k) = cw(k)*(p(i,j,k+1)-p(i,j,k))
      enddo
      enddo
      enddo

      return

      end
c---------------------------------------------------------------------



c---- subroutine correctpres2---------------N. Beratlis-18 May 2009---
C
C     PURPOSE: Correct pressure assuming value at body is known
C
C---------------------------------------------------------------------
      subroutine correctpres2(p,nx,ny,nz,ip,jp,kp,ipmtrx,jpmtrx,kpmtrx,pmtrx
     &     ,pim,nxp,nyp,nzp,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect,fp,elem
     &     ,nfacet,pb,node,nnodes,xc,yc,zc,limp,mimp,tlevel,dt,ibd,mbd)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer nx,ny,nz,limp,mimp,nfacet,ibd,mbd,nnodes
      integer ip(nfcmax),jp(nfcmax),kp(nfcmax)
      integer ipmtrx(nsm,nfcmax),jpmtrx(nsm,nfcmax),kpmtrx(nsm,nfcmax)
      integer mrkp(nfcmax),dirp(nfcmax,2)
      integer pindx(nsp,nfcmax)
      integer fp(nfcmax),elem(3,nfcmax)
      real    tlevel,dt
      real    xc(nx),yc(ny),zc(nz)
      real    p(nx,ny,nz)
      real    node(3,nnodes),pb(nnodes)
      real    pim(nfcmax),xnp(nfcmax),ynp(nfcmax),znp(nfcmax),dpdnn(nfcmax)
      real    unvect(3,nfcmax)
      real    nxp(nfcmax,2),nyp(nfcmax,2),nzp(nfcmax,2)
      real    pmtrx(nsp,nsp,nfcmax)
c
c.... Local arrays
      integer i,j,k,im,itr
      integer indx(nsp)
      real    pint,psi,d
      real    xpcar,ypcar,zpcar,xp2car,yp2car,zp2car,s1,s2
      real    a(nsp,nsp),b(nsp)
      real    dpdnn1(nfcmax)
      real    vt(3,3),pt(3),v1(3),v2(3)
c
c.... Functions
      real triangle_interp,anglerad,vecmag,dotproduct
c
      
c      open(unit=20,file='pintrs.plt',form='formatted')
c      write(20,'(A,A)') 'VARIABLES = "X", "Y", "Z", "ANGLE", "P"'
c      write(20,'(A,I8,A)') 'ZONE I=',mimp,', DATAPACKING=POINT'

      do im=limp+1,limp+mimp

        i = ip(im)
        j = jp(im)
        k = kp(im)

        if(icyl==0) then
          xpcar = xc(i)
          ypcar = yc(j)
          if(mrkp(im)==1) then
            xp2car = nxp(im,1)
            yp2car = nyp(im,1)
            zp2car = nzp(im,1)
          elseif(mrkp(im)==2) then
            xp2car = xc(ipmtrx(1,im))
            yp2car = yc(jpmtrx(1,im))
            zp2car = zc(kpmtrx(1,im))
          endif
        else
          xpcar = xc(i)*cos(yc(j))
          ypcar = xc(i)*sin(yc(j))
          if(mrkp(im)==1) then
            xp2car = nxp(im,1)
            yp2car = nyp(im,1)
            zp2car = nzp(im,1)
          elseif(mrkp(im)==2) then
            xp2car = xc(ipmtrx(1,im))*cos(yc(jpmtrx(1,im)))
            yp2car = xc(ipmtrx(1,im))*sin(yc(jpmtrx(1,im)))
            zp2car = zc(kpmtrx(1,im))
          endif
        endif
        zpcar = zc(k)
       
        !Find value of pressure at triangle intersection
        itr = fp(im)
        vt(1,1) = node(1,elem(1,itr))
        vt(1,2) = node(2,elem(1,itr))
        vt(1,3) = node(3,elem(1,itr))
        vt(2,1) = node(1,elem(2,itr))
        vt(2,2) = node(2,elem(2,itr))
        vt(2,3) = node(3,elem(2,itr))
        vt(3,1) = node(1,elem(3,itr))
        vt(3,2) = node(2,elem(3,itr))
        vt(3,3) = node(3,elem(3,itr))
        pt(1) = pb(elem(1,itr))
        pt(2) = pb(elem(2,itr))
        pt(3) = pb(elem(3,itr))
        pint = triangle_interp(xnp(im),ynp(im),znp(im),vt,pt,im)

c        psi = anglerad(-znp(i),sqrt(xnp(i)**2. + ynp(i)**2.))
c        write(20,'(12(1X,F14.8))') xnp(i),ynp(i),znp(i),psi,pint

        b(1) = pint
        b(2) = pim(im)

        v1(1) = xnp(im) - xpcar
        v1(2) = ynp(im) - ypcar
        v1(3) = znp(im) - zpcar

        v2(1) = xp2car - xpcar
        v2(2) = yp2car - ypcar
        v2(3) = zp2car - zpcar

        s1 = vecmag(v1,3)*sign(1.0,dotproduct(v1,unvect(:,itr),3))
        s2 = vecmag(v2,3)*sign(1.0,dotproduct(v2,unvect(:,itr),3))

        a(1:2,1) = 1
        a(1,2) = s1
        a(2,2) = s2

        call ludcmp(a,nsp,nsp,indx,d)
        call lubksb(a,2,2,indx,b)

        p(i,j,k) = b(1)

      enddo

c      close(20)

      return

      end
c---------------------------------------------------------------------


C---- subroutine presinterior----------------N. Beratlis-18 May 2009--
C
C     PURPOSE: Set the pressure in the interior to zero
C
c---------------------------------------------------------------------
      subroutine presinterior(p,nx,ny,nz,flagp,ibd,nbd)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'

      integer nx,ny,nz,ibd,nbd
      integer flagp(nx,ny,nz,nbd)
      real    p(nx,ny,nz)
      
      integer i,j,k

c      write(6,*) 'presinterior:',myrank,ibmin(ibd),ibmax(ibd),jbmin(ibd),jbmax(ibd),kbmin(ibd),kbmax(ibd)
c      write(6,*) 'flagpo:',myrank,flagp(4,2,126,ibd)
      do k=kbmin(ibd),kbmax(ibd)
      do j=jbmin(ibd),jbmax(ibd)
      do i=ibmin(ibd),ibmax(ibd)
        if((flagp(i,j,k,ibd)==ibd)) p(i,j,k) = 0.
c        if(i==4 .AND. k+myrank*(nz-2)==126) write(6,*) 'myrank=',myrank,flagp(i,j,k,ibd),p(i,j,k)
      enddo
      enddo
      enddo

      CALL REFRESHBC(P,NX*NY,NZ)
      
      return

      end
c---------------------------------------------------------------------

C---- SUBROUTINE FLAGFE ---------------------------------------------
C
C     PURPOSE: Flag points as interior, exterior and forcing
c
      SUBROUTINE FLAGFE(FLAG,FLAGUO,FLAGUFE,NX,NY,NZ,NBD,IBD)
C
c---------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
c
      integer nx,ny,nz,nbd
      integer ibd
      integer flag(nx,ny,nz)
      integer flaguo(nx,ny,nz,nbd),flagufe(nx,ny,nz)
c
c....local variables
      integer nxm,nym,nzm,i,j,k,nfe
c....Timing variables
      INTEGER COUNTS,RATE
      REAL    CLOCK
c
c      WRITE(6,*) 'INSIDE FLAGU'
      nxm = nx-1
      nym = ny-1
      nzm = nz-1
c
      nfe=0
      flag = 0
c     
c.....define true field extension points
      do k=kbmin(ibd),kbmax(ibd)
      do j=jbmin(ibd),jbmax(ibd)
      do i=ibmin(ibd),ibmax(ibd)
        if(flagufe(i,j,k)/=0.AND.flaguo(i,j,k,ibd)==0) then
          flag(i,j,k) = 1
          nfe = nfe+1
        endif
      enddo
      enddo
      enddo

      if(nfe>0) write(6,*) 'WARNING: FIELD EXTENSION',nfe

      flagufe = flag

      flagufe(:,1,:)  = flagufe(:,NY-1,:)
      flagufe(:,ny,:) = flagufe(:,2,:)
c
      return
      end
c---------------------------------------------------------------------



C---- subroutine geomp --------------------N. Beratlis-23 Jan. 2009---
C
C     PURPOSE: Find interpolation points and intersection with immersed 
C     object for all forcing points.
C
C---------------------------------------------------------------------
      subroutine geomp(xu,yu,xu_car,yu_car,zu,vertex,vertexc,unvect,lim,mim
     &    ,iim,jim,kim,itr,nim,xnim,ynim,znim,nxim,nyim,nzim,mrk,dir
     &    ,flag,flagp,nx,ny,nz,nbd,nfacet,ibd,icom,nfp)
C
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
      integer nx,ny,nz,nbd,nfacet,ibd,nim,icom,nfp
      integer lim(nfp),mim(nfp)
      integer iim(nfcmax),jim(nfcmax),kim(nfcmax),itr(nfcmax)
      integer mrk(nfcmax),dir(nfcmax,2)
      real    xu(nx),yu(ny),zu(nz)
      real    xu_car(nx,ny),yu_car(nx,ny)
      real    xnim(nfcmax),ynim(nfcmax),znim(nfcmax)
      real    nxim(nfcmax,2),nyim(nfcmax,2),nzim(nfcmax,2)
      integer flagp(nx,ny,nz,nbd),flag(nx,ny,nz)
      real    vertex(3,3,nfacet),vertexc(3,nfacet),unvect(3,nfacet)
c
c....local variables
      integer ic,nn,nimp,icntl,nt,ntg,ntmaxg,ntming,nfps,mimg,mimmax,mimmin,nq,iord,mimp
      integer im,jm,km,kmg,i,j,k,ilb,ile,im1,ii,in,nzg,intrs,gintrs,dintrs,ilp,intrs1
      integer iext,jext,kext,iext1,jext1,kext1
      integer icount1,icount2,icount3,icount4,isten
      integer nnrm,nnrmg,ndia,ndiag,ngrd,ngrdg,nvld,nmrk
      integer nnrmrchk,nnrmrchkg,ndiarchk,ndiarchkg,ngrdrchk,ngrdrchkg
      integer flag2(nx,ny,2)
      integer icount,icountnorm,icountgrid,icountdiag
      integer icountnormg,icountgridg,icountdiagg
      integer icountnovld(nfp),icountnovldg
      integer icountnrmnostenc(nfp),icountnrmnostencg
      integer icountnrmnointrs(nfp),icountnrmnointrsg
      integer icountgrdnostenc(nfp),icountgrdnostencg
      integer icountgrdnointrs(nfp),icountgrdnointrsg
      integer icountdianostenc(nfp),icountdianostencg
      integer icountdianointrs(nfp),icountdianointrsg

      integer nclocks
      real    clocktemp
      real, dimension(:), allocatable :: clock,clockg,clockgmin,clockgmax
      REAL, DIMENSION(:), ALLOCATABLE :: zug
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: INNF
      REAL, DIMENSION(:,:), ALLOCATABLE :: queryf
      INTEGER, DIMENSION(:), ALLOCATABLE :: ORD,ORDC,ORDN
c      real    queryf(3,)
c
c.... functions
      INTEGER forcpts,rdfn_forcpts_p,grid_intr_p,norm_intr_p,diag_intr_p
      INTEGER recheckfaces_p,recheckgrdpts_p,recheckdiagintr_p
      REAL    tclock
c
c      open(unit=10,file='temp.'//char(myrank+48)//'.txt')
c      write(10,*) myrank,'inside geomp, icom',icom
c      close(10)

      nclocks = 27
      allocate(clock(nclocks),clockg(nclocks),clockgmin(nclocks),clockgmax(nclocks))

      clock = 0.0
      clockg = 0.0
      clockgmin = 0.0
      clockgmax = 0.0

      clock(1) = tclock()

      nnrm = 0
      ngrd = 0
      ndia = 0
      nnrmg = 0
      ngrdg = 0
      ndiag = 0
      nnrmrchk = 0
      ngrdrchk = 0
      ndiarchk = 0
      nnrmrchkg = 0
      ngrdrchkg = 0
      ndiarchkg = 0
      nvld = 0
      nmrk = 0

      NN=50 !50 !50

      mimg = 0
      mimmax = 0
      mimmin = 0

      nimp = nim
      nt = 0
      icntl = 0
      if(icyl==1) then
        if(abs(xu_car(1,1))/=0.0) icntl=1
      endif
      
      !Initialize
      icount = 0
      icountnorm = 0
      icountgrid = 0
      icountdiag = 0
      icountnovld = 0
      icountnrmnostenc = 0
      icountnrmnointrs = 0
      icountgrdnostenc = 0
      icountgrdnointrs = 0
      icountdianostenc = 0
      icountdianointrs = 0

      clocktemp = tclock()
      nzg = (nz-2)*mysize+2
      ALLOCATE(zug(nzg))
      call calc_zg(zu,zug,nz,nzg)
      clock(2) = tclock() - clocktemp

      lim=0
      mim=0

      clocktemp = tclock()
      lim(1) = nim
      mim(1) = forcpts(flagp,nx,ny,nz,ibd,nbd,iim,jim,kim,nim)      
      nim = nimp+mim(1)
      nfps = nim

      IF(nim>nfcmax) THEN
        WRITE(6,'(A,1X,I9)') 'ERROR: Increase nfcmax to at least',nim
        CALL MPI_FINALIZE(IERR)
      ENDIF
      clock(3) = tclock() - clocktemp


      mrk(1:nim)=0
      nq = nim
      allocate(innf(nq,nn),queryf(3,nq),ord(nq),ordc(nq),ordn(nq))
      do i=1,nq
        queryf(1,i) = xu_car(iim(i),jim(i))
        queryf(2,i) = yu_car(iim(i),jim(i))
        queryf(3,i) = zu(kim(i)) 
        ordc(i) = i
      enddo

      clocktemp = tclock()
      if(icom==0 .OR. icom==3) then
        CALL COMMITTEE(QUERYF, NQ, VERTEXC, NFACET, INNF, NN, 0, MYRANK)
      endif
      clock(4) = tclock() - clocktemp

      clocktemp = tclock()
      call committee(queryf, nq, vertexc, nfacet, innf, nn, 1, myrank)
      clock(6) = clock(6) + tclock() - clocktemp

      flag(:,:,:) = flagp(:,:,:,ibd)
      
      clock(12) = tclock()
      do ilp=1,nfp

        clocktemp = tclock()
        call fillextra_cell_z(flag,flag2,nx,ny,nz,myleft,myright,mpi_comm_eddy)
        clock(5) = clock(5) + tclock() - clocktemp

        clocktemp = tclock()
        mimp = mim(ilp)
        call mpi_allreduce(mimp,mimmax,1,mpi_integer,mpi_max,mpi_comm_eddy,ierr)
        clock(25) = clock(25) + tclock() - clocktemp

        if(ilp>1) then
          lim(ilp)=sum(mim(1:ilp-1))
        endif

        nmrk = 0

        ii = lim(ilp)+1
        do while(ii<=lim(ilp)+mimp)

          i = ordc(ii-lim(ilp))
          im = iim(i)
          jm = jim(i)
          km = kim(i)     

          if(mrk(i)==-1) then
            clocktemp = tclock()
            intrs = recheckfaces_p(flag,flag2,xu,yu,zu,nx,ny,nz,nxim(i,:),nyim(i,:),nzim(i,:),dir(i,:),icyl)            
            nnrmrchk = nnrmrchk+1
            clock(21) = clock(21) + tclock() - clocktemp
          elseif(mrk(i)==0) then
            clocktemp = tclock()
            intrs=norm_intr_p(flag,flag2,nx,ny,nz,ibd,nbd,vertex,unvect,nfacet
     &      ,xu_car,yu_car,zu,zug,nzg,innf(i,:),nn,xnim(i),ynim(i),znim(i)
     &      ,nxim(i,:),nyim(i,:),nzim(i,:),im,jm,km,dir(i,:),itr(i),icntl,nt,clock(16),5)
            nnrm = nnrm+1
            clock(7) = clock(7) + tclock() - clocktemp
          else
            intrs = 0
          endif

          kmg = km+myrank*(nz-2)

          if(intrs==1) then
            mrk(i)=1
            icountnorm = icountnorm+1
            nvld = nvld+1
            ord(nvld) = i
          elseif(intrs==-1 .AND. ilp<nfp) then
            mrk(i)=-1

c            clocktemp = tclock()
c            iord = ord(ii)
c            ord(ii:nq-1) = ord(ii+1:nq)
c            ord(nq) = iord
c            clock(26) = clock(26) + tclock()-clocktemp
c            ii = ii-1
            mim(ilp) = mim(ilp)-1
            mim(ilp+1)=mim(ilp+1)+1
            nmrk = nmrk+1
            ordn(nmrk) = i
          else

            if(intrs==0) then
              icountnrmnointrs(ilp) = icountnrmnointrs(ilp)+1
            else
              icountnrmnostenc(ilp) = icountnrmnostenc(ilp)+1
            endif


            if(mrk(i)==-2) then
              clocktemp = tclock()
              gintrs = recheckgrdpts_p(flag,flag2,nx,ny,nz,im,jm,km,nxim(i,1),nyim(i,1),nzim(i,1))              
              clock(22) = clock(22) + tclock() - clocktemp
              ngrdrchk = ngrdrchk+1
            elseif(mrk(i)==0) then
              clocktemp = tclock()
              gintrs = grid_intr_p(flag,flag2,nx,ny,nz,ibd,nbd,vertex,unvect,nfacet
     &             ,xu_car,yu_car,zu,innf(i,:),nn,xnim(i),ynim(i),znim(i)
     &             ,nxim(i,1),nyim(i,1),nzim(i,1),im,jm,km,itr(i))
              ngrd = ngrd+1
              clock(8) = clock(8) + tclock() - clocktemp
            else
              gintrs = 0
            endif

            if(gintrs==1) then
              mrk(i)=2
              icountgrid = icountgrid+1
              nvld = nvld+1
              ord(nvld) = i
            elseif(gintrs==-2 .AND. ilp<nfp) then
              mrk(i) = -2
c              clocktemp = tclock()
c              iord = ord(ii)
c              ord(ii:nq-1) = ord(ii+1:nq)
c              ord(nq) = iord
c              clock(26) = clock(26) + tclock()-clocktemp
c              ii = ii-1
              mim(ilp) = mim(ilp)-1
              mim(ilp+1)=mim(ilp+1)+1
              nmrk = nmrk+1
              ordn(nmrk) = i
            else !stencil has body points or no intrs. with body was found

              if(gintrs==0) then
                icountgrdnointrs(ilp) = icountgrdnointrs(ilp)+1
              else
                icountgrdnostenc(ilp) = icountgrdnostenc(ilp)+1
              endif

              if(mrk(i)==-3) then
                clocktemp = tclock()
                dintrs=recheckdiagintr_p(flag,flag2,xu,yu,zu,nx,ny,nz,im,jm,km,nxim(i,:),nyim(i,:),nzim(i,:),dir(i,1),icyl)
                ndiarchk= ndiarchk+1
                clock(23) = clock(23) + tclock() - clocktemp
              else
                clocktemp = tclock()
                dintrs=diag_intr_p(flag,flag2,nx,ny,nz,ibd,nbd,vertex,unvect,nfacet
     &               ,xu_car,yu_car,zu,zug,nzg,innf(i,:),nn,xnim(i),ynim(i),znim(i)
     &               ,nxim(i,:),nyim(i,:),nzim(i,:),im,jm,km,dir(i,:),itr(i))
                ndia = ndia+1
                clock(9) = clock(9) + tclock() - clocktemp
              endif

              if(dintrs==1) then
                mrk(i)=3
                icountdiag = icountdiag+1
                nvld = nvld+1
                ord(nvld) = i
              elseif(dintrs==-1 .AND. ilp<nfp) then
                mrk(i) = -3
c                clocktemp = tclock()
c                iord = ord(ii)
c                ord(ii:nq-1) = ord(ii+1:nq)
c                ord(nq) = iord
c                clock(26) = clock(26) + tclock()-clocktemp
c                ii = ii-1
                mim(ilp) = mim(ilp)-1
                mim(ilp+1)=mim(ilp+1)+1
                nmrk = nmrk+1
                ordn(nmrk)= i
              else
c                write(6,*) 'point cannot be used',im,jm,km,mrk(i),intrs,gintrs,dintrs

                mrk(i) = 0
                flagp(im,jm,km,ibd) = ibd

c                clocktemp = tclock()
c                ord(ii:nq-1) = ord(ii+1:nq)
c                clock(26) = clock(26) + tclock()-clocktemp
c                ii = ii-1

                mim(ilp) = mim(ilp)-1

                if(dintrs==0) then
                  icountdianointrs(ilp) = icountdianointrs(ilp)+1
                else
                  icountdianostenc(ilp) = icountdianostenc(ilp)+1
                endif

              endif
            endif
          endif
          
          ii = ii+1

        enddo

        clocktemp = tclock()
        call physical_mrk2flag1(flag,nx,ny,nz,mrk,iim,jim,kim,ord,nq,lim(ilp),mim(ilp))
        clock(24) = clock(24) + tclock() - clocktemp

        clocktemp = tclock()
        if(mimp>0) then
          call yrefresh3darray(flag,nx,ny,nz)
        endif

        if(mimmax>0) then
          call refreshflag(flag(1,1,1),nx,ny,nz)
        endif
        clock(10) = clock(10) + tclock() - clocktemp

        icountnovld(ilp) = count(mrk(lim(ilp)+1:lim(ilp)+mim(ilp))==0)

        ordc = ordn

        clocktemp = tclock()
        if(icountnovld(ilp)>0) then  
c          mim(ilp) = rdfn_forcpts_p(flagp(:,:,:,ibd),nx,ny,nz,ibd,iim,jim,kim,xnim,ynim,znim
c     &          ,nxim,nyim,nzim,mrk,dir,itr,ord,nq,lim(ilp),mim(ilp))
          nim = sum(mim)        
        endif
        clock(11) = tclock() - clocktemp

      enddo

      clocktemp = tclock()
      call rord_forcpts_p(iim,jim,kim,xnim,ynim,znim,nxim,nyim,nzim
     &     ,mrk,dir,itr,ord,nq,lim(1),lim(nfp)+mim(nfp))
      clock(27) = clock(27) + tclock() - clocktemp
c      write(6,*) '3. iim=',ilp,iim(839),ord(839),iim(ord(839)),lim,mim,lim+mim

      if(nfps>0) nt = int(nt/real(nfps))

      clock(12) = tclock() - clock(12)

      clocktemp = tclock()
      if(icom>=2) then
        CALL COMMITTEE(QUERYF, NQ, VERTEXC, NFACET, INNF, NN, 2, MYRANK)
      endif
      clock(13) = tclock() - clocktemp

      DEALLOCATE(INNF,queryf,ord,ordn,ordc)
      DEALLOCATE(zug)

      clocktemp = tclock()
      CALL MPI_REDUCE(SUM(MIM),ICOUNT,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(MIM),MIMMAX,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(MIM),MIMMIN,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTNORM,ICOUNTNORMG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTGRID,ICOUNTGRIDG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTDIAG,ICOUNTDIAGG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(ICOUNTNOVLD),ICOUNTNOVLDG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(ICOUNTNRMNOSTENC),ICOUNTNRMNOSTENCG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(ICOUNTGRDNOSTENC),ICOUNTGRDNOSTENCG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(ICOUNTDIANOSTENC),ICOUNTDIANOSTENCG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(ICOUNTNRMNOINTRS),ICOUNTNRMNOINTRSG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(ICOUNTGRDNOINTRS),ICOUNTGRDNOINTRSG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(ICOUNTDIANOINTRS),ICOUNTDIANOINTRSG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(NT,NTG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(NT,NTMAXG,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(NT,NTMING,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(nnrm,nnrmg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ngrd,ngrdg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ndia,ndiag,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(nnrmrchk,nnrmrchkg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ngrdrchk,ngrdrchkg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ndiarchk,ndiarchkg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      clock(14) = tclock() - clocktemp

      clocktemp = tclock()
      IF(MYRANK.EQ.0 .AND. .true.) THEN
        OPEN(UNIT=16, FILE='stats_imb.dat', FORM='FORMATTED'
     &      ,POSITION='APPEND')
        write(16,*) 'Ave. no. of triangles searched for normal intrs.',nt
        write(16,'(2(A,1x,I7),A,F6.2,2A,I7,A,F6.2,2A,I7,A,F6.2,A)') 
     &       ' pressure forcing points:',ICOUNT
     &       ,', forced along normal:',ICOUNTNORMG
     &       ,' (',100.*REAL(ICOUNTNORMG)/REAL(ICOUNT),'/100)'
     &       ,', forced along grid points:',ICOUNTGRIDG
     &       ,' (',100.*REAL(ICOUNTGRIDG)/REAL(ICOUNT),'/100)'
     &       ,', forced along diagonals:',ICOUNTDIAGG
     &       ,' (',100.*REAL(ICOUNTDIAGG)/REAL(ICOUNT),'/100)'
        write(16,'(A,F8.2,2(A,I7))') 'No. of forc. pts., ave:',real(icount/real(mysize))
     &       ,' ,max:',mimmax,', min:',mimmin
        
        if(icountnovldg>0) then
          write(16,'(A,I7,A)') 'WARNING:',icountnovldg
     &          ,' pressure forcing points cannot be used:'
          write(16,'(A,2(I7,A))') '  along normal '
     &         ,icountnrmnointrsg,' due to no intrs., '
     &         ,icountnrmnostencg,' due to no physical stencil, ' 
          write(16,'(A,2(I7,A))') '  along gridlines '
     &         ,icountgrdnointrsg,' due to no intrs., '
     &         ,icountgrdnostencg,' due to no physical stencil' 
          write(16,'(A,2(I7,A))') '  along diagonal '
     &         ,icountdianointrsg,' due to no intrs.,'
     &         ,icountdianostencg,' due to no physical stencil' 
        endif
        CLOSE(16)
      ENDIF
      clock(15) = tclock() - clocktemp
        
      clock(1) = tclock() - clock(1)

      IF(.TRUE.) THEN
      CALL MPI_REDUCE(CLOCK,CLOCKG,NCLOCKS,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(CLOCK,CLOCKGMIN,NCLOCKS,MTYPE,MPI_MIN,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(CLOCK,CLOCKGMAX,NCLOCKS,MTYPE,MPI_MAX,0,MPI_COMM_EDDY,IERR)
      clockg = clockg/real(mysize)

      where(clockg==0.0) clockg=1.e-8
      where(clockgmin==0.0) clockgmin=1.e-8
      where(clockgmax==0.0) clockgmax=1.e-8
      

      IF(MYRANK==0) THEN
        OPEN(UNIT=16,FILE='clock.dat',FORM='FORMATTED',POSITION='APPEND')
        write(16,'(A)') '---- geomp: ----------------------------------'
        write(16,'(A)') '    Task/Time           Ave. (sec/%)        Max. (sec/%)   
     &     Min (sec/%)'
        write(16,'(A,3(6x,F8.4))') 'Total               :',clockg(1),clockgmax(1),clockgmin(1)
        i=2
        write(16,906) 'Calc z global       :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=3
        write(16,906) 'Forcpts points      :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=4
        write(16,906) 'Initialize ANN      :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=6
        write(16,906) 'Find nearest neigh. :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=12
        write(16,906) 'Find intrs./stencil :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=5
        write(16,906) '  -Fill extra flag2 :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(12),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(12),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(12),'%)'
        i=7
        write(16,906) '  -Norm. intrs.     :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(12),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(12),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(12),'%)'
        write(16,'(A,3(1x,I8))') '     No. of triangles searched:',ntg,ntmaxg,ntming
        i=16
        write(16,906) '     -segtriint     :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(7),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(7),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(7),'%)'
        i=17
        write(16,906) '     -interp_point_p:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(7),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(7),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(7),'%)'
        i=18
        write(16,906) '       -ray ext fnts:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(7),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(7),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(7),'%)'
        i=19
        write(16,906) '       -ray_face_int:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(7),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(7),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(7),'%)'
        i=20
        write(16,906) '       -fluidface   :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(7),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(7),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(7),'%)'
        i=8
        write(16,906) '  -Grid intrs.      :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(12),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(12),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(12),'%)'
        i=9
        write(16,906) '  -Diag. intrs.     :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(12),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(12),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(12),'%)'
        i=21
        write(16,906) '  -Norm intrs rchk  :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(12),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(12),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(12),'%)'
        write(16,'(3(A,1x,I8),A)') '    (No. nrm=',nnrmg,',no. grd=',ngrdg
     &       ,', no. dia=',ndiag,')'
        write(16,'(3(A,1x,I8),A)') '    (Rchk: no. nrm=',nnrmrchkg
     &       ,',no. grd=',ngrdrchkg,', no. dia=',ndiarchkg,')'
        i=22
        write(16,906) '  -Grid intrs rchk  :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(12),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(12),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(12),'%)'
        i=23
        write(16,906) '  -Diag intrs rchk  :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(12),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(12),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(12),'%)'
        i=24
        write(16,906) '  -Physical mrk2flag:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(12),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(12),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(12),'%)'
        i=25
        write(16,906) '  -MPI_allreduce MIM:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(12),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(12),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(12),'%)'
        i=26
        write(16,906) '  -Reorder ord array:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(12),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(12),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(12),'%)'
        i=10
        write(16,906) '  -Refresh flag arr :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(12),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(12),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(12),'%)'
        i=11
        write(16,906) '  -Redefine forcpts :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(12),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(12),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(12),'%)'
        i=27
        write(16,906) 'Reorder forpts      :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=13
        write(16,906) 'Delete ANN structure:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=14
        write(16,906) 'MPI_REDUCE count arr:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=15
        write(16,906) 'I/O operations      :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'

        write(16,'(A)') '----------------------------------------------'
        close(16)
      ENDIF

      ENDIF

      deallocate(clock,clockg,clockgmin,clockgmax)

      RETURN

 906  format(A,3(1x,F8.4,A,F6.2,A,2x))

      END
C-------------------------------------------------------------------------



C---- subroutine flagp_uoi ---------------------N. Beratlis-May 28 2009---
C
C     PURPOSE: locates i,j,k pressure points on the cartesian grid near
C     the immersed boundary to be modified. Takes into account both outer
C     inner velocity forcing points. 
C     General routine that work for any boundary shape   
C
C-------------------------------------------------------------------------
      subroutine flagp_uoi(flaguo,flagvo,flagwo,flagui,flagvi,flagwi,flagpo
     &     ,flagpi,nx,ny,nz,nbd,ibd)
C
c      implicit none
      INCLUDE 'common.h'
      INCLUDE 'immersed.h'
      INCLUDE 'mpif.h'
c
      integer nx,ny,nz,nbd,ibd
      integer flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),flagwo(nx,ny,nz,nbd)
      integer flagui(nx,ny,nz,nbd),flagvi(nx,ny,nz,nbd),flagwi(nx,ny,nz,nbd)
      integer flagpo(nx,ny,nz,nbd),flagpi(nx,ny,nz,nbd)
c
c....local variables
      integer i,j,k,i1,i2,j1,j2
      integer itst(6),ifld,ibnd,ibdy
      integer nxm,nym,nzm
c
      nxm=nx-1
      nym=ny-1
      nzm=nz-1
c
c.....set up the flag for true pressure cells
      
      FLAGPO = 0 !
      FLAGPI = 0
c      WHERE(FLAGPO==IBD) FLAGPI=FLAGPO
      i1 = ibmin(ibd)
      i2 = ibmax(ibd)
      j1 = jbmin(ibd)
      j2 = jbmax(ibd)
      IF(icyl==1) THEN
        j1=2
        j2=ny-1
      ENDIF
c
c.....define the pressure boundary points
      DO I=max(2,i1),i2
      DO J=j1,j2
      DO K=KBMIN(IBD),KBMAX(IBD)
c        
        itst(1) = min(flaguo(i  ,j,k,ibd),flagui(i  ,j,k,ibd))
        itst(2) = min(flaguo(i-1,j,k,ibd),flagui(i-1,j,k,ibd))
        itst(3) = min(flagwo(i,j,k  ,ibd),flagwi(i,j,k  ,ibd))
        itst(4) = min(flagwo(i,j,k-1,ibd),flagwi(i,j,k-1,ibd))
        itst(5) = min(flagvo(i,j  ,k,ibd),flagvi(i,j  ,k,ibd))
        itst(6) = min(flagvo(i,j-1,k,ibd),flagvi(i,j-1,k,ibd))
c
        ifld=count(itst==0)
        ibdy=count(itst> 0)
        ibnd=count(itst< 0)
c
c... If pressure forcing point is surrounded by fluid or forcing velocity points it is correct (make it fluid)
c        IF((IBDY+IBND)==0) THEN
c          IF(FLAGPO(I,J,K,IBD)==-IBD) FLAGPO(I,J,K,IBD) = 0 
c          IF(FLAGPI(I,J,K,IBD)==-IBD) FLAGPI(I,J,K,IBD) = 0
c        ENDIF
c
c.....outer pressure boundary points
        IF((ibdy+ibnd)/=0.AND.ifld/=0) THEN
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

c        if(i==215.AND.k==146) then
c          write(6,*) 'flagp:',i,j,k,ifld,ibdy,ibnd,flagpo(i,j,k,ibd),flagpi(i,j,k,ibd)
c        endif

      ENDDO
      ENDDO
      ENDDO

      IF(icyl==1 .AND. ibmin(ibd)==1) THEN
c.....define the pressure boundary points
        i = 1
        DO j=j1,j2
        DO k=KBMIN(IBD),KBMAX(IBD)
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
c... If pressure forcing point is surrounded by fluid or forcing velocity points it is correct (make it fluid)
c          IF((IBDY+IBND)==0) THEN
c            IF(FLAGPO(I,J,K,IBD)==-IBD) FLAGPO(I,J,K,IBD) = 0
c            IF(FLAGPI(I,J,K,IBD)==-IBD) FLAGPI(I,J,K,IBD) = 0
c          ENDIF
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

        ENDDO
        ENDDO
      ENDIF

c
c.....Periodic boundary conditions in Y-direction
      FLAGPI(:,1,:,IBD) = FLAGPI(:,NY-1,:,IBD)
      FLAGPO(:,1,:,IBD) = FLAGPO(:,NY-1,:,IBD)
      FLAGPI(:,NY,:,IBD) = FLAGPI(:,2,:,IBD)
      FLAGPO(:,NY,:,IBD) = FLAGPO(:,2,:,IBD)
c
c.....Axis
c      IF(ITYPE(1)==300) THEN
c        DO J=1,NY
c          FLAGPO(1,J,:,IBD) = FLAGPO(IX1,JSYM(J),:,IBD)
c          FLAGPI(1,J,:,IBD) = FLAGPI(IX1,JSYM(J),:,IBD)
c        ENDDO
c      ENDIF


     
      RETURN

      END
c
c---------------------------------------------------------------------


c---- function extmag----------------------N. Beratlis-01 Jun. 2009---
C
      real function extmag(xu,yu,zu,i,j,k,nx,ny,nz,icyl)
c
      implicit none
c
c.... Input/Output Arrays
      integer nx,ny,nz,i,j,k,icyl
      real    xu(nx,ny),yu(nx,ny),zu(nz)
c
c.... Local arrays
      real    dr,dr1,dr2,darc,dz,r1,r2

      if(icyl==1) then
c
        if(i==1) then
          dr  = sqrt(xu(i+1,j)**2.+yu(i+1,j)**2)
     &        - sqrt(xu(i  ,j)**2.+yu(i  ,j)**2)
        elseif(i==nx) then
          dr  = sqrt(xu(i  ,j)**2.+yu(i  ,j)**2)
     &        - sqrt(xu(i-1,j)**2.+yu(i-1,j)**2)
        else
          dr1 = sqrt(xu(i+1,j)**2.+yu(i+1,j)**2)
     &        - sqrt(xu(i  ,j)**2.+yu(i  ,j)**2)
          dr2 = sqrt(xu(i  ,j)**2.+yu(i  ,j)**2)
     &        - sqrt(xu(i-1,j)**2.+yu(i-1,j)**2)
          if(dr2>0) then
            dr = min(dr1,dr2)
          else
            dr = dr1
          endif

        endif

        if(j==1) then
          darc = sqrt( (xu(i,j)-xu(i,j+1))**2. 
     &               + (yu(i,j)-yu(i,j+1))**2.)
        elseif(j==ny) then
          darc = sqrt( (xu(i,j)-xu(i,j-1))**2. 
     &               + (yu(i,j)-yu(i,j-1))**2.)
        else
          darc = sqrt( (xu(i,j)-xu(i,j+1))**2. 
     &               + (yu(i,j)-yu(i,j+1))**2.)
        endif

        if(k==1) then
          dz = abs(zu(k+1)-zu(k))
        elseif(k==nz) then
          dz = abs(zu(k)-zu(k-1))
        else
          dz = min(abs(zu(k+1)-zu(k)),abs(zu(k)-zu(k-1)))
        endif

        if(darc>0.0) then
          if(dr>0.0) then
            extmag = min(dr,darc,dz)
          else
            extmag = min(darc,dz)
          endif
        else
          if(dr>0.0) then
            extmag = min(dr,dz)
          else
            extmag = dz
          endif
        endif

c        if(extmag==0.0) then
c          write(6,*) 'extmag:',i,j,k,dr,darc,dz
c          stop
c        endif

      endif

      return

      end

c---------------------------------------------------------------------



C---- subroutine shear_body2 ------------------N. Beratlis-16 Mar. 2010-
C
C     PURPOSE: Calculate shear strain on surface on immersed body
C
C-----------------------------------------------------------------------
      subroutine shear_body2(uo,vo,wo,nx,ny,nz,xu,yv,zw,xc,yc,zc
     &     ,xu_car,yu_car,xv_car,yv_car,xc_car,yc_car
     &     ,vertexc,unvect,nfacet,sxxb,syyb,szzb,sxyb,sxzb,syzb,oyb,mrks
     &     ,ibd,nbd,tlevel,io)

c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'      
      include 'mpif.h'
c
c.... Input/Output Arrays
      INTEGER nx,ny,nz,ibd,nbd,nfacet,io
      REAL    tlevel
      INTEGER mrks(nfacet,6)
      REAL    xu(nx),yv(ny),zw(nz),xc(nx),yc(ny),zc(nz)
      REAL    xu_car(nx,ny),yu_car(nx,ny),xv_car(nx,ny),yv_car(nx,ny)
     &       ,xc_car(nx,ny),yc_car(nx,ny)
      REAL    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
      REAL    vertexc(3,nfacet),unvect(3,nfacet)
      REAL    sxxb(nfacet),syyb(nfacet),szzb(nfacet)
     &       ,sxyb(nfacet),sxzb(nfacet),syzb(nfacet)
     &       ,oyb(nfacet)
c
c.... Local arrays
      INTEGER i,j,k,ii,ifacet,iext,jext,kext,nzg,iflag,n1,dir,ksign
      INTEGER i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4
      REAL    sxx,syy,szz,sxy,syx,szy,syz,sxz,szx,up
      REAL    dudn,dvdn,dwdn,r,theta
      REAL    dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dn,dx,dy,dz
      REAL    xp,yp,zp,xint,yint,zint,uint,vint,wint,ub,vb,wb,psi,a1,a2,a3,unorm,amin
      INTEGER icount(6),icountg(6)
      REAL    q(3),rvec(3)
      REAL    sxxb1(nfacet),syyb1(nfacet),szzb1(nfacet)
     &       ,sxyb1(nfacet),sxzb1(nfacet),syzb1(nfacet)
     &       ,oyb1(nfacet)
      INTEGER, SAVE :: icall=0
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: mrk,mrkg

      LOGICAL condx,condy,condz
c
c.... Functions
      INTEGER body2fluid,body2grid2
      REAL    vecmag,ubd,vbd,wbd,interp_cellface,anglerad,interp3d
      REAL    calc_sxx,calc_syy,calc_szz,calc_sxy,calc_syx,calc_syz,calc_szy,calc_sxz,calc_szx
      

      amin = 0.75
      icount=0

      n1=4

      ALLOCATE(mrk(mb(ibd),3),mrkg(mb(ibd),3))
      mrk = n1+1
      mrkg = n1+1

      mrks = 1

      sxxb1 = 0.0
      syyb1 = 0.0
      szzb1 = 0.0
      sxyb1 = 0.0
      sxzb1 = 0.0
      syzb1 = 0.0

      q = 0.0
      rvec = 0.0
      xint = 0.0
      yint = 0.0
      zint = 0.0

      DO ifacet=lb(ibd)+1,lb(ibd)+mb(ibd)


        ii = ifacet-lb(ibd)

        xp = vertexc(1,ifacet)
        yp = vertexc(2,ifacet)
        zp = vertexc(3,ifacet)

        q(1) = xp
        q(2) = yp
        q(3) = zp   

        rvec(1) = unvect(1,ifacet)
        rvec(2) = unvect(2,ifacet)
        rvec(3) = unvect(3,ifacet)
        ksign = 0
        if(rvec(3)>0) ksign=1

        if(icyl==0) then
          condx = yp>=ymin.AND.yp<=ymax .AND. zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
          condz = yp>=ymin.AND.yp<=ymax .AND. zp>=zw(2-ksign) .AND. zp<zw(nz-ksign)
        else
          condx = zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
          condz = zp>=zw(2-ksign) .AND. zp<zw(nz-ksign)
        endif

        if(icyl==0) then
          theta = yp
          r = xp
        else
          theta = anglerad(xp,yp)
          r = sqrt(xp**2. + yp**2.) 
        endif

        if(condx) then

          icount(1) = icount(1)+1
          !Calculate derivatives of u
          mrk(ii,1) = body2grid2(q,rvec,xu_car,yu_car,xu,yc,zc
     &         ,uo,nx,ny,nz,xint,yint,zint,uint,n1,nbd,icyl,ifacet,amin)           
          if(mrk(ii,1)<=n1) then
            if(io==1) then
              ub = ubd(xp,yp,zp,tlevel,ibd)
            else
              ub = 0.0
            endif

            dn = sqrt( (xp-xint)**2. + (yp-yint)**2. + (zp-zint)**2.)
            sxxb1(ifacet) = (uint-ub)/dn
          endif

          !Calculate derivatives of v
          mrk(ii,2) = body2grid2(q,rvec,xv_car,yv_car,xc,yv,zc
     &         ,vo,nx,ny,nz,xint,yint,zint,vint,n1,nbd,icyl,ifacet,amin)           
          if(mrk(ii,2)<=n1) then
            if(io==1) then
              vb = vbd(xp,yp,zp,tlevel,ibd)
            else
              vb = 0.0
            endif

            dn = sqrt( (xp-xint)**2. + (yp-yint)**2. + (zp-zint)**2.)
            syyb1(ifacet) = (vint-vb)/dn
          endif

        endif

        if(condz) then

          icount(2) = icount(2)+1
          !Calculate derivatives of w
          mrk(ii,3) = body2grid2(q,rvec,xc_car,yc_car,xc,yc,zw
     &         ,wo,nx,ny,nz,xint,yint,zint,wint,n1,nbd,icyl,ifacet,amin)           
          if(mrk(ii,3)<=n1) then
            if(io==1) then
              wb = wbd(xp,yp,zp,tlevel,ibd)
            else
              wb = 0.0
            endif

            dn = sqrt( (xp-xint)**2. + (yp-yint)**2. + (zp-zint)**2.)
            szzb1(ifacet) = (wint-wb)/dn
          endif

        endif

      ENDDO

      IF(mysize>1) THEN
        CALL MPI_ALLREDUCE(mrk,mrkg,3*mb(ibd),MPI_INTEGER,MPI_MIN,MPI_COMM_EDDY,IERR)
      ELSE
        mrkg = mrk
      ENDIF

      DO ifacet=lb(ibd)+1,lb(ibd)+mb(ibd)

        ii = ifacet-lb(ibd)

        dudx = 0.0
        dudy = 0.0
        dudz = 0.0
        dvdx = 0.0
        dvdy = 0.0
        dvdz = 0.0
        dwdx = 0.0
        dwdy = 0.0
        dwdz = 0.0

        unorm = vecmag(unvect(:,ifacet),3)
        a1 = unvect(1,ifacet)/unorm
        a2 = unvect(2,ifacet)/unorm
        a3 = unvect(3,ifacet)/unorm

        psi = anglerad(unvect(1,ifacet),unvect(2,ifacet))

        if(mrk(ii,1)<=mrkg(ii,1) .AND. mrk(ii,1)<=n1) then
          mrk(ii,1) = 1
          icount(2) = icount(2)+1
        else
          sxxb1(ifacet) = 0.0
          mrk(ii,1) = 0
        endif

        if(mrk(ii,2)<=mrkg(ii,2) .AND. mrk(ii,2)<=n1) then
          mrk(ii,2) = 1
          icount(4) = icount(4)+1
        else
          syyb1(ifacet) = 0.0
          mrk(ii,2) = 0
        endif

        if(mrk(ii,3)<=mrkg(ii,3) .AND. mrk(ii,3)<=n1) then
          mrk(ii,3) = 1
          icount(6) = icount(6)+1
        else
          szzb1(ifacet) = 0.0
          mrk(ii,3) = 0
        endif

        if(icyl==0) then

          dudn = sxxb1(ifacet)
          dvdn = syyb1(ifacet)

          dudx = dudn*a1
          dudy = dudn*a2
          dudz = dudn*a3

          dvdx = dvdn*a1
          dvdy = dvdn*a2
          dvdz = dvdn*a3

        else

          ub = sxxb1(ifacet)
          vb = syyb1(ifacet)

          dudn = ub*cos(psi)-vb*sin(psi)
          dudx = dudn*a1
          dudy = dudn*a2
          dudz = dudn*a3

          dvdn = ub*sin(psi)+vb*cos(psi)
          dvdx = dvdn*a1
          dvdy = dvdn*a2
          dvdz = dvdn*a3

        endif
        
        dwdn = szzb1(ifacet)

        dwdx = dwdn*a1
        dwdy = dwdn*a2
        dwdz = dwdn*a3


        sxxb1(ifacet) = 2.*dudx
        syyb1(ifacet) = 2.*dvdy
        szzb1(ifacet) = 2.*dwdz
        sxyb1(ifacet) = dudy+dvdx
        sxzb1(ifacet) = dudz+dwdx
        syzb1(ifacet) = dvdz+dwdy
        oyb1(ifacet) = dudz-dwdx

      ENDDO



      IF(mysize>1) THEN
        CALL MPI_REDUCE(sxxb1,sxxb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(syyb1,syyb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(szzb1,szzb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(sxyb1,sxyb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(sxzb1,sxzb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(syzb1,syzb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(oyb1,oyb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(mrk,mrkg,3*mb(ibd),MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(icount,icountg,6,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      ELSE
        sxxb = sxxb1
        syyb = syyb1
        szzb = szzb1
        sxyb = sxyb1
        sxzb = sxzb1
        syzb = syzb1
        oyb = oyb1
        mrkg = mrk
        icountg = icount
      ENDIF

      IF(myrank==0 .AND. iolvl>0) then

        if(ibm>1 .OR. icall==0) then
          if(icount(1)/=icountg(2) .OR. icount(3)/=icountg(4) .OR. icount(5)/=icountg(6)) then
            open(unit=16,file='stats_imb.dat',position='append'
     &            ,form='formatted')
            write(16,*) 'Number of triangles=',mb(ibd)
            write(16,*) 'U-velocity:'
            write(16,*) 'Number of tri. within domain',icount(1)
            write(16,*) 'Number of dudn points extrap.',icountg(2)
            write(16,*) 'V-velocity:'
            write(16,*) 'Number of tri. within domain',icount(3)
            write(16,*) 'Number of dvdn points extrap.',icountg(4)
            write(16,*) 'W-velocity:'
            write(16,*) 'Number of tri. within domain',icount(5)
            write(16,*) 'Number of dwdn points extrap.',icountg(6)
            close(16)
          endif
        endif
      ENDIF

      icall = icall+1

      DEALLOCATE(mrk,mrkg)

      RETURN

      END
C------------------------------------------------------------------------




C---- subroutine shear_body3 ------------------N. Beratlis-16 Mar. 2010-
C
C     PURPOSE: Calculate shear strain on surface on immersed body
C
C-----------------------------------------------------------------------
      subroutine shear_body3(uo,vo,wo,nx,ny,nz,xu,yv,zw,xc,yc,zc
     &     ,xu_car,yu_car,xv_car,yv_car,xc_car,yc_car,vertexc,unvect
     &     ,nfacet,dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb
     &     ,mrks,ibd,nbd,ilb,ile,tlevel,io)

c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      INTEGER nx,ny,nz,ibd,nbd,nfacet,io,ilb,ile
      REAL    tlevel
      INTEGER mrks(9,nfacet)
      REAL    xu(nx),yv(ny),zw(nz),xc(nx),yc(ny),zc(nz)
      REAL    xu_car(nx,ny),yu_car(nx,ny),xv_car(nx,ny),yv_car(nx,ny)
     &       ,xc_car(nx,ny),yc_car(nx,ny)
      REAL    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
      REAL    vertexc(3,nfacet),unvect(3,nfacet)
      REAL    dudxb(nfacet),dudyb(nfacet),dudzb(nfacet)
     &       ,dvdxb(nfacet),dvdyb(nfacet),dvdzb(nfacet)
     &       ,dwdxb(nfacet),dwdyb(nfacet),dwdzb(nfacet)
c
c.... Local arrays
      INTEGER i,j,k,ii,ifacet,iext,jext,kext,nzg,iflag,n1,dir,ksign
      REAL    sxx,syy,szz,sxy,syx,szy,syz,sxz,szx,up
      REAL    dudn,dvdn,dwdn,r,theta
      REAL    dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dn,dx,dy,dz
      REAL    xp,yp,zp,xint,yint,zint,uint,vint,wint,ub,vb,wb,psi,a1,a2,a3,unorm,amin
      INTEGER icount(6),icountg(6)
      REAL    q(3),rvec(3)
      INTEGER, SAVE :: icall=0

      LOGICAL condx,condy,condz
c
c.... Functions
      INTEGER body2fluid,body2grid2
      REAL    vecmag,ubd,vbd,wbd,interp_cellface,anglerad,interp3d      

      amin = 0.75
      icount=0

      n1=4

      mrks(:,ilb:ile) = 0      !!!!!!! instead of mrks = 0
      dudxb(ilb:ile) = 0.0     !!!!!!! instead of dudxb = 0.0
      dudyb(ilb:ile) = 0.0
      dudzb(ilb:ile) = 0.0
      dvdxb(ilb:ile) = 0.0
      dvdyb(ilb:ile) = 0.0
      dvdzb(ilb:ile) = 0.0
      dwdxb(ilb:ile) = 0.0
      dwdyb(ilb:ile) = 0.0
      dwdzb(ilb:ile) = 0.0

      q = 0.0
      rvec = 0.0
      xint = 0.0
      yint = 0.0
      zint = 0.0

      DO ifacet=ilb,ile

        dudn=0
        dvdn=0
        dwdn=0

        xp = vertexc(1,ifacet)
        yp = vertexc(2,ifacet)
        zp = vertexc(3,ifacet)

        q(1) = xp
        q(2) = yp
        q(3) = zp   

        rvec(1) = unvect(1,ifacet)
        rvec(2) = unvect(2,ifacet)
        rvec(3) = unvect(3,ifacet)
        ksign = 0
        if(rvec(3)>0) ksign=1

        if(icyl==0) then
          condx = yp>=ymin.AND.yp<=ymax .AND. zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
          condz = yp>=ymin.AND.yp<=ymax .AND. zp>=zw(2-ksign) .AND. zp<zw(nz-ksign)
        else
          condx = zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
          condz = zp>=zw(2-ksign) .AND. zp<zw(nz-ksign)
        endif

        if(icyl==0) then
          theta = yp
          r = xp
        else
          theta = anglerad(xp,yp)
          r = sqrt(xp**2. + yp**2.) 
        endif

        if(condx) then

          icount(1) = icount(1)+1
          !Calculate derivatives of u
c in BODY2GRID2 the interpolated value for UO is UINT
          mrks(1,ifacet) = body2grid2(q,rvec,xu_car,yu_car,xu,yc,zc
     &         ,uo,nx,ny,nz,xint,yint,zint,uint,n1,nbd,icyl,ifacet,amin)
          if(mrks(1,ifacet)<=n1) then
c a valid intersection has been found in BODY2GRID2
            icount(2) = icount(2)+1
            if(io==1) then
c moving boundary
              ub = ubd(xp,yp,zp,tlevel,ibd)
            else
c stationary boundary
              ub = 0.0
            endif

            dn = sqrt( (xp-xint)**2. + (yp-yint)**2. + (zp-zint)**2.)
c the normal derivative for U is evaluated
            dudn = (uint-ub)/dn
          endif

          icount(3) = icount(3)+1     !!!!!!!
          !Calculate derivatives of v
          mrks(2,ifacet) = body2grid2(q,rvec,xv_car,yv_car,xc,yv,zc
     &         ,vo,nx,ny,nz,xint,yint,zint,vint,n1,nbd,icyl,ifacet,amin)           
          if(mrks(2,ifacet)<=n1) then
            icount(4) = icount(4)+1
            if(io==1) then
              vb = vbd(xp,yp,zp,tlevel,ibd)
            else
              vb = 0.0
            endif

            dn = sqrt( (xp-xint)**2. + (yp-yint)**2. + (zp-zint)**2.)
            dvdn = (vint-vb)/dn
          endif

        endif

        if(condz) then

          icount(5) = icount(5)+1
          !Calculate derivatives of w
          mrks(3,ifacet) = body2grid2(q,rvec,xc_car,yc_car,xc,yc,zw
     &         ,wo,nx,ny,nz,xint,yint,zint,wint,n1,nbd,icyl,ifacet,amin)
          if(mrks(3,ifacet)<=n1) then
            mrks(7:9,ifacet)=1
            icount(6) = icount(6)+1
            if(io==1) then
              wb = wbd(xp,yp,zp,tlevel,ibd)
            else
              wb = 0.0
            endif

            dn = sqrt( (xp-xint)**2. + (yp-yint)**2. + (zp-zint)**2.)
            dwdn = (wint-wb)/dn
          endif

        endif

        unorm = vecmag(unvect(:,ifacet),3)
c components of the unit normal vector
        a1 = unvect(1,ifacet)/unorm
        a2 = unvect(2,ifacet)/unorm
        a3 = unvect(3,ifacet)/unorm

        psi = anglerad(unvect(1,ifacet),unvect(2,ifacet))

c the components of the normal derivatives along X,Y,Z are evaluated
        if(icyl==0) then

          if(mrks(1,ifacet)>0) mrks(1:3,ifacet)=1
          dudxb(ifacet) = dudn*a1
          dudyb(ifacet) = dudn*a2
          dudzb(ifacet) = dudn*a3

          if(mrks(2,ifacet)>0) mrks(4:6,ifacet)=1
          dvdxb(ifacet) = dvdn*a1
          dvdyb(ifacet) = dvdn*a2
          dvdzb(ifacet) = dvdn*a3

        else

          if(mrks(1,ifacet)>0 .AND. mrks(2,ifacet)>0) mrks(1:6,ifacet)=1

          ub = dudn
          vb = dvdn

          dudn = ub*cos(psi)-vb*sin(psi)
          dudxb(ifacet) = dudn*a1
          dudyb(ifacet) = dudn*a2
          dudzb(ifacet) = dudn*a3

          dvdn = ub*sin(psi)+vb*cos(psi)
          dvdxb(ifacet) = dvdn*a1
          dvdyb(ifacet) = dvdn*a2
          dvdzb(ifacet) = dvdn*a3

        endif
        
        dwdxb(ifacet) = dwdn*a1
        dwdyb(ifacet) = dwdn*a2
        dwdzb(ifacet) = dwdn*a3

      ENDDO

      IF(mysize>1) THEN
       CALL MPI_REDUCE(icount,icountg,6,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      ELSE
        icountg = icount
      ENDIF

      IF(myrank==0 .AND. iolvl>0) then

        if(ibm>1 .OR. icall==0) then
!          if(icount(1)/=icountg(2) .OR. icount(3)/=icountg(4) .OR. icount(5)/=icountg(6)) then    !!!!!!
            open(unit=16,file='stats_imb.dat',position='append'
     &            ,form='formatted')
            write(16,*) 'Number of triangles=',mb(ibd)
            write(16,*) 'U-velocity:'
            write(16,*) 'Number of tri. within domain',icountg(1)    !!!!!! instead of icount(1)
            write(16,*) 'Number of dudn points extrap.',icountg(2)
            write(16,*) 'V-velocity:'
            write(16,*) 'Number of tri. within domain',icountg(3)    !!!!!! instead of icount(3)
            write(16,*) 'Number of dvdn points extrap.',icountg(4)
            write(16,*) 'W-velocity:'
            write(16,*) 'Number of tri. within domain',icountg(5)    !!!!!! instead of icount(5)
            write(16,*) 'Number of dwdn points extrap.',icountg(6)
            close(16)
!          endif    !!!!!!
        endif
      ENDIF

      icall = icall+1

      RETURN

      END
C------------------------------------------------------------------------



C---- function body2grid2---------------------N. Beratlis-16 Mar. 2010---
C
C     PURPOSE: From a point q on the immersed boundary surface find 1
C     intersection along surface normal r and interpolate value of uo
C     at the intersection. Intersection cannot be outside local subdomain.
C
C-----------------------------------------------------------------------
      integer function body2grid2(q,rvec,xu_car,yu_car,xu,yu,zu
     &     ,uo,nx,ny,nz,xint,yint,zint,uint,n1,nbd,icyl,ifacet,amin)
c
      IMPLICIT NONE

      REAL, PARAMETER :: pi=3.141592653589793
c
c.... Input/Output Arrays      
      INTEGER nx,ny,nz,nzg,dir,nbd,n1,icyl,ifacet
      REAL    xint,yint,zint
      REAL    q(3),rvec(3)
      REAL    uo(nx,ny,nz)
      REAL    xu_car(nx,ny),yu_car(nx,ny),xu(nx),yu(ny),zu(nz)
c
c.... Local Arrays
      INTEGER i,j,k,kg,i1,j1,k1,i2,j2,k2,iext,jext,kext,next
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
!      amin = 0.
      xint1 = 0.0
      yint1 = 0.0
      zint1 = 0.0

      xp = q(1)
      yp = q(2)
      zp = q(3)
      
      iphy = 0
      intrs = 1

      iflag1=0
      iflag2=0
      iflag3=0

      clocktemp = tclock()
      call ijk_xyz(xp,yp,zp,xu,yu,zu,nx,ny,nz,i,j,k,icyl)

      ds = mindxdydz(xu,yu,zu,nx,ny,nz,i,j,k,icyl)

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
      iflag1 = ray_face_int(q,rvec,xu_car,yu_car,zu,nx,ny,nz
     &     ,i1,j1,k1,1,xint1,yint1,zint1,icyl)
      clock(2) = clock(2) + tclock() - clocktemp

      !Compute distance from body to intersection
      s = sqrt( (q(1)-xint1)**2. + (q(2)-yint1)**2. + (q(3)-zint1)**2. )

      if(iflag1==1) then
        i1 = i+iext
        j1 = j
        k1 = k
        dir = 1
        if(s>amin*ds) then
          iphy = iphy+1
          uint = interp_cellface(xint1,yint1,zint1,xu,yu,zu,uo
     &         ,nx,ny,nz,icyl,dir)
          xint = xint1
          yint = yint1
          zint = zint1
        endif

        if(iphy==0) then
          iflagt=0
          do while(iflagt==0)
            intrs = intrs+1

            if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &        ,rvec,xu_car,yu_car,zu,nx,ny,nz
     &        ,i1,j1,k1,i2,j2,k2,dir,dir,icyl,ifacet)==1) then
              s = sqrt( (q(1)-xint2)**2. + (q(2)-yint2)**2. + (q(3)-zint2)**2. )
              if(s>amin*ds) then
                iphy = iphy+1
                iflagt=1
                uint = interp_cellface(xint2,yint2,zint2,xu,yu,zu,uo
     &               ,nx,ny,nz,icyl,dir)
                xint = xint2
                yint = yint2
                zint = zint2
              endif
            endif
            xint1 = xint2
            yint1 = yint2
            zint1 = zint2  
            i1 = i2
            j1 = j2
            k1 = k2
            if(intrs>=n1 .OR. (dir==3 .AND. (k2-k==-1 .OR. k2-k==2)) ) iflagt=1
          enddo
        endif

      else

        i1 = i
        j1 = j+jext
        k1 = k

        iflag2 = ray_face_int(q,rvec,xu_car,yu_car,zu,nx,ny,nz
     &     ,i,j+jext,k,2,xint1,yint1,zint1,icyl)
        s = sqrt( (q(1)-xint1)**2. + (q(2)-yint1)**2. + (q(3)-zint1)**2. )

        if(iflag2==1) then
          i1 = i
          j1 = j+jext
          k1 = k
          dir = 2

          if(s>amin*ds) then
            iphy = iphy+1
            uint = interp_cellface(xint1,yint1,zint1,xu,yu,zu,uo
     &           ,nx,ny,nz,icyl,dir)
            xint = xint1
            yint = yint1
            zint = zint1
          endif

          if(iphy==0) then
            iflagt=0
            do while(iflagt==0)
              intrs = intrs+1
              
              if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &             ,rvec,xu_car,yu_car,zu,nx,ny,nz
     &             ,i1,j1,k1,i2,j2,k2,dir,dir,icyl,ifacet)==1) then
                s = sqrt( (q(1)-xint2)**2. + (q(2)-yint2)**2. + (q(3)-zint2)**2. )
                if(s>amin*ds) then
                  iphy = iphy+1
                  iflagt=1
                  uint = interp_cellface(xint2,yint2,zint2,xu,yu,zu,uo
     &                 ,nx,ny,nz,icyl,dir)
                  xint = xint2
                  yint = yint2
                  zint = zint2
                endif
              endif
              xint1 = xint2
              yint1 = yint2
              zint1 = zint2  
              i1 = i2
              j1 = j2
              k1= k2
              if(intrs>=n1 .OR. (dir==3 .AND. (k2-k==-1 .OR. k2-k==2)) ) iflagt=1
            enddo
          endif

        else

          k1 = k+kext 
          dz = zu(k1)-zu(k)
          call vec_ijkext(q,rvec,iext,jext,kext,icyl,3,dz,dy)

          iflag3 = ray_face_int(q,rvec,xu_car,yu_car,zu,nx,ny,nz
     &         ,i,j,k+kext,3,xint1,yint1,zint1,icyl)
          s = sqrt( (q(1)-xint1)**2. + (q(2)-yint1)**2. + (q(3)-zint1)**2. )

          if(iflag3==1) then
            i1 = i
            j1 = j
            k1 = k+kext
            dir = 3

            if(s>amin*ds) then
              iphy = iphy+1
              uint = interp_cellface(xint1,yint1,zint1,xu,yu,zu,uo
     &             ,nx,ny,nz,icyl,dir)
              xint = xint1
              yint = yint1
              zint = zint1
            endif

            if(iphy==0) then
              iflagt=0
              do while(iflagt==0)
                intrs = intrs+1

                if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &               ,rvec,xu_car,yu_car,zu,nx,ny,nz
     &               ,i1,j1,k1,i2,j2,k2,dir,dir,icyl,ifacet)==1) then
                  s = sqrt( (q(1)-xint2)**2. + (q(2)-yint2)**2. + (q(3)-zint2)**2. )
                  if(s>amin*ds) then
                    iphy = iphy+1
                    iflagt=1
                    uint = interp_cellface(xint2,yint2,zint2,xu,yu,zu,uo
     &                   ,nx,ny,nz,icyl,dir)
                    xint = xint2
                    yint = yint2
                    zint = zint2
                  endif
                endif
                xint1 = xint2
                yint1 = yint2
                zint1 = zint2  
                i1 = i2
                j1 = j2
                k1= k2
                if(intrs>=n1 .OR. (dir==3 .AND. (k2-k==-1 .OR. k2-k==2))) iflagt=1
              enddo
            endif
          endif
        endif
      endif


c      if(q(3)>0.0099 .AND. q(3)<0.01) then
c        write(6,*) 'body2grid2',iflag1,iflag2,iflag3,dir
c      endif

!      if(iflag1==0.AND.iflag2==0.AND.iflag3==0) then
!        write(6,*) ifacet,'vel_body2grid2 no face intersection',iphy,q,rvec,amin
!      endif


      if(iphy==1) then
        body2grid2 = intrs
      else
        body2grid2 = n1+1
      endif


      RETURN

      END

C------------------------------------------------------------------------


C---- function body2fluid2 -------------------N. Beratlis-16 Mar. 2010---
C
C     PURPOSE: From a point q on the immersed boundary surface find 1
C     intersection along surface normal r and interpolate value of uo
C     at the intersection. Intersection cannot be outside local subdomain.
C
C-----------------------------------------------------------------------
      integer function body2fluid2(q,rvec,xu_car,yu_car,xu,yu,zu
     &     ,uo,flaguo,nx,ny,nz,xint,yint,zint,uint,n1,nbd,icyl,ifacet,amin)
c
      IMPLICIT NONE

      REAL, PARAMETER :: pi=3.141592653589793
c
c.... Input/Output Arrays      
      INTEGER nx,ny,nz,nzg,dir,nbd,n1,icyl,ifacet
      REAL    xint,yint,zint
      INTEGER flaguo(nx,ny,nz,nbd)
      REAL    q(3),rvec(3)
      REAL    uo(nx,ny,nz)
      REAL    xu_car(nx,ny),yu_car(nx,ny),xu(nx),yu(ny),zu(nz)
c
c.... Local Arrays
      INTEGER i,j,k,kg,i1,j1,k1,i2,j2,k2,iext,jext,kext,next
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
      amin = 0.
      xint1 = 0.0
      yint1 = 0.0
      zint1 = 0.0

      xp = q(1)
      yp = q(2)
      zp = q(3)
      
      iphy = 0
      intrs = 1

      iflag1=0
      iflag2=0
      iflag3=0

      clocktemp = tclock()
c the computational cell containing Q is found
      call ijk_xyz(xp,yp,zp,xu,yu,zu,nx,ny,nz,i,j,k,icyl)

      ds = mindxdydz(xu,yu,zu,nx,ny,nz,i,j,k,icyl)

      theta = yu(j)
      dy = 2.0*pi/real(ny-2)
      rm = xu(i)
      if(icyl==1 .AND. i==1) rm=0.0
c the extensions of the indices are established
      call vec_ijkext(q,rvec,iext,jext,kext,icyl,1,rm,dy)

      i1 = i+iext
      j1 = j
      k1 = k

      clock(1) = clock(1) + tclock() - clocktemp

      clocktemp = tclock()
c an intersection with a grid face normal to X is eventually found
      iflag1 = ray_face_int(q,rvec,xu_car,yu_car,zu,nx,ny,nz
     &     ,i1,j1,k1,1,xint1,yint1,zint1,icyl)
      clock(2) = clock(2) + tclock() - clocktemp

      !Compute distance from body to intersection
      s = sqrt( (q(1)-xint1)**2. + (q(2)-yint1)**2. + (q(3)-zint1)**2. )

      if(iflag1==1) then
        i1 = i+iext
        j1 = j
        k1 = k
        dir = 1

c another check: physicalface is 0 if the grid face involves body points
        if(s>amin*ds .AND. physicalface(flaguo,nx,ny,nz,i1,j1,k1,nbd,dir)/=0) then
          iphy = iphy+1
c the variable is interpolated at the intersection point
          uint = interp_cellface(xint1,yint1,zint1,xu,yu,zu,uo
     &         ,nx,ny,nz,icyl,dir)
          xint = xint1
          yint = yint1
          zint = zint1
        endif

        if(iphy==0) then
c if the intersection is not acceptable
          iflagt=0
          do while(iflagt==0)
            intrs = intrs+1
c it looks for another intersection along the same direction
            if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &        ,rvec,xu_car,yu_car,zu,nx,ny,nz
     &        ,i1,j1,k1,i2,j2,k2,dir,dir,icyl,ifacet)==1) then
              s = sqrt( (q(1)-xint2)**2. + (q(2)-yint2)**2. + (q(3)-zint2)**2. )
c the same check of above is used
              if(s>amin*ds .AND. physicalface(flaguo,nx,ny,nz,i2,j2,k2,nbd,dir)/=0) then
                iphy = iphy+1
                iflagt=1
                uint = interp_cellface(xint2,yint2,zint2,xu,yu,zu,uo
     &               ,nx,ny,nz,icyl,dir)
                xint = xint2
                yint = yint2
                zint = zint2
              endif
            endif
            xint1 = xint2
            yint1 = yint2
            zint1 = zint2  
            i1 = i2
            j1 = j2
            k1 = k2
            if(intrs>=n1 .OR. (dir==3 .AND. (k2-k==-1 .OR. k2-k==2)) ) iflagt=1
          enddo
        endif

      else

c if an intersection with a grid face normal to X is not found
        i1 = i
        j1 = j+jext
        k1 = k

        iflag2 = ray_face_int(q,rvec,xu_car,yu_car,zu,nx,ny,nz
     &     ,i,j+jext,k,2,xint1,yint1,zint1,icyl)
        s = sqrt( (q(1)-xint1)**2. + (q(2)-yint1)**2. + (q(3)-zint1)**2. )

        if(iflag2==1) then
          i1 = i
          j1 = j+jext
          k1 = k
          dir = 2

          if(s>amin*ds .AND. physicalface(flaguo,nx,ny,nz,i1,j1,k1,nbd,dir)/=0) then
            iphy = iphy+1
            uint = interp_cellface(xint1,yint1,zint1,xu,yu,zu,uo
     &           ,nx,ny,nz,icyl,dir)
            xint = xint1
            yint = yint1
            zint = zint1
          endif

          if(iphy==0) then
            iflagt=0
            do while(iflagt==0)
              intrs = intrs+1
              
              if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &             ,rvec,xu_car,yu_car,zu,nx,ny,nz
     &             ,i1,j1,k1,i2,j2,k2,dir,dir,icyl,ifacet)==1) then
                s = sqrt( (q(1)-xint2)**2. + (q(2)-yint2)**2. + (q(3)-zint2)**2. )
                if(s>amin*ds .AND. physicalface(flaguo,nx,ny,nz,i2,j2,k2,nbd,dir)/=0) then
                  iphy = iphy+1
                  iflagt=1
                  uint = interp_cellface(xint2,yint2,zint2,xu,yu,zu,uo
     &                 ,nx,ny,nz,icyl,dir)
                  xint = xint2
                  yint = yint2
                  zint = zint2
                endif
              endif
              xint1 = xint2
              yint1 = yint2
              zint1 = zint2  
              i1 = i2
              j1 = j2
              k1= k2
              if(intrs>=n1 .OR. (dir==3 .AND. (k2-k==-1 .OR. k2-k==2)) ) iflagt=1
            enddo
          endif

        else

c if an intersection with a grid face normal to X or Y is not found
          k1 = k+kext 
          dz = zu(k1)-zu(k)
          call vec_ijkext(q,rvec,iext,jext,kext,icyl,3,dz,dy)

          iflag3 = ray_face_int(q,rvec,xu_car,yu_car,zu,nx,ny,nz
     &         ,i,j,k+kext,3,xint1,yint1,zint1,icyl)
          s = sqrt( (q(1)-xint1)**2. + (q(2)-yint1)**2. + (q(3)-zint1)**2. )

          if(iflag3==1) then
            i1 = i
            j1 = j
            k1 = k+kext
            dir = 3

            if(s>amin*ds .AND. physicalface(flaguo,nx,ny,nz,i1,j1,k1,nbd,dir)/=0) then
              iphy = iphy+1
              uint = interp_cellface(xint1,yint1,zint1,xu,yu,zu,uo
     &             ,nx,ny,nz,icyl,dir)
              xint = xint1
              yint = yint1
              zint = zint1
            endif

            if(iphy==0) then
              iflagt=0
              do while(iflagt==0)
                intrs = intrs+1

                if(grid2grid_intr(xint1,yint1,zint1,xint2,yint2,zint2
     &               ,rvec,xu_car,yu_car,zu,nx,ny,nz
     &               ,i1,j1,k1,i2,j2,k2,dir,dir,icyl,ifacet)==1) then
                  s = sqrt( (q(1)-xint2)**2. + (q(2)-yint2)**2. + (q(3)-zint2)**2. )
                  if(s>amin*ds .AND. physicalface(flaguo,nx,ny,nz,i2,j2,k2,nbd,dir)/=0) then
                    iphy = iphy+1
                    iflagt=1
                    uint = interp_cellface(xint2,yint2,zint2,xu,yu,zu,uo
     &                   ,nx,ny,nz,icyl,dir)
                    xint = xint2
                    yint = yint2
                    zint = zint2
                  endif
                endif
                xint1 = xint2
                yint1 = yint2
                zint1 = zint2  
                i1 = i2
                j1 = j2
                k1= k2
                if(intrs>=n1 .OR. (dir==3 .AND. (k2-k==-1 .OR. k2-k==2))) iflagt=1
              enddo
            endif
          endif
        endif
      endif


!      if(iflag1==0.AND.iflag2==0.AND.iflag3==0) then
!        write(6,*) ifacet,'vel_body2fluid no face intersection',iphy,q,rvec,amin
!      endif


      if(iphy==1) then
        body2fluid2 = intrs
      else
        body2fluid2 = n1+1
      endif


      RETURN

      END

C------------------------------------------------------------------------



C---- subroutine shear_body1 ------------------N. Beratlis-26 Dec. 2008-
C
C     PURPOSE: Calculate shear strain on surface on immersed body
C
C-----------------------------------------------------------------------
      subroutine shear_body1(uo,vo,wo,nx,ny,nz,xu,yv,zw,xc,yc,zc
     &     ,xu_car,yu_car,xv_car,yv_car,xc_car,yc_car
     &     ,vertexc,unvect,nfacet,sxxb,syyb,szzb,sxyb,sxzb,syzb,oyb,mrks
     &     ,ibd,nbd,tlevel,io)

c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'      
      include 'mpif.h'
c
c.... Input/Output Arrays
      INTEGER nx,ny,nz,ibd,nbd,nfacet,io
      REAL    tlevel
      INTEGER mrks(nfacet,6)
      REAL    xu(nx),yv(ny),zw(nz),xc(nx),yc(ny),zc(nz)
      REAL    xu_car(nx,ny),yu_car(nx,ny),xv_car(nx,ny),yv_car(nx,ny)
     &       ,xc_car(nx,ny),yc_car(nx,ny)
      REAL    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
      REAL    vertexc(3,nfacet),unvect(3,nfacet)
      REAL    sxxb(nfacet),syyb(nfacet),szzb(nfacet)
     &       ,sxyb(nfacet),sxzb(nfacet),syzb(nfacet)
     &       ,oyb(nfacet)
c
c.... Local arrays
      INTEGER i,j,k,ii,ifacet,iext,jext,kext,nzg,iflag,n1,dir
      INTEGER i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4
      REAL    sxx,syy,szz,sxy,syx,szy,syz,sxz,szx,up
      REAL    dudn,dvdn,dwdn,r,theta
      REAL    dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dn,dx,dy,dz
      REAL    xp,yp,zp,xint,yint,zint,uint,vint,wint,ub,vb,wb,psi,a1,a2,a3,unorm,amin
      INTEGER icount(6),icountg(6)
      REAL    q(3),rvec(3)
      REAL    sxxb1(nfacet),syyb1(nfacet),szzb1(nfacet)
     &       ,sxyb1(nfacet),sxzb1(nfacet),syzb1(nfacet)
     &       ,oyb1(nfacet)
      INTEGER, SAVE :: icall=0

      LOGICAL condx,condy,condz
c
c.... Functions
      INTEGER vel_body2fluid,body2grid
      REAL    vecmag,ubd,vbd,wbd,interp_cellface,anglerad,interp3d
      REAL    calc_sxx,calc_syy,calc_szz,calc_sxy,calc_syx,calc_syz,calc_szy,calc_sxz,calc_szx
      

      icount=0

      n1=1
      mrks = 1

      sxxb1 = 0.0
      syyb1 = 0.0
      szzb1 = 0.0
      sxyb1 = 0.0
      sxzb1 = 0.0
      syzb1 = 0.0

      q = 0.0
      rvec = 0.0
      xint = 0.0
      yint = 0.0
      zint = 0.0

      DO ifacet=lb(ibd)+1,lb(ibd)+mb(ibd)


        ii = ifacet-lb(ibd)

        xp = vertexc(1,ifacet)
        yp = vertexc(2,ifacet)
        zp = vertexc(3,ifacet)

        q(1) = xp
        q(2) = yp
        q(3) = zp   

        rvec(1) = unvect(1,ifacet)
        rvec(2) = unvect(2,ifacet)
        rvec(3) = unvect(3,ifacet)


        if(icyl==0) then
          condx = yp>=ymin.AND.yp<=ymax .AND. zp>=zc(2) .AND. zp<zc(nz)
          condz = yp>=ymin.AND.yp<=ymax .AND. zp>=zw(2) .AND. zp<zw(nz)
        else
          condx = zp>=zc(2) .AND. zp<zc(nz)
          condz = zp>=zw(2) .AND. zp<zw(nz)
        endif

        sxx = 0.0
        syy = 0.0
        szz = 0.0
        sxy = 0.0
        syx = 0.0
        sxz = 0.0
        szx = 0.0
        syz = 0.0
        szy = 0.0

        if(icyl==0) then
          theta = yp
          r = xp
        else
          theta = anglerad(xp,yp)
          r = sqrt(xp**2. + yp**2.) 
        endif

        if(condx) then
           
          sxx = calc_sxx(uo,xp,yp,zp,xu,yc,zc,nx,ny,nz,icyl)
          syy = calc_syy(vo,xp,yp,zp,xc,yv,zc,nx,ny,nz,icyl)
          if(icyl==1) then 
            up = interp3d(xp,yp,zp,xu,yc,zc,uo,nx,ny,nz,icyl)
            syy = syy+up/r
           endif
          syx = calc_syx(vo,xp,yp,zp,xc,yv,zc,nx,ny,nz,icyl)
          sxy = calc_sxy(uo,xp,yp,zp,xu,yc,zc,nx,ny,nz,icyl)
          syz = calc_syz(vo,xp,yp,zp,xc,yv,zc,nx,ny,nz,icyl)
          sxz = calc_sxz(uo,xp,yp,zp,xu,yc,zc,nx,ny,nz,icyl)
        endif

        if(condz) then
          szz = calc_szz(wo,xp,yp,zp,xc,yc,zw,nx,ny,nz,icyl)
          szy = calc_szy(wo,xp,yp,zp,xc,yc,zw,nx,ny,nz,icyl)
          szx = calc_szx(wo,xp,yp,zp,xc,yc,zw,nx,ny,nz,icyl)
        endif

        sxxb1(ifacet) = sxx
        syyb1(ifacet) = syy
        szzb1(ifacet) = szz
        sxyb1(ifacet) = sxy+syx
        sxzb1(ifacet) = sxz+szx
        syzb1(ifacet) = syz+szy
        oyb1(ifacet) = 2.0*(sxz-szx)

        if(abs(xp)>0.495) then
          write(6,*) 'shear_body1',ifacet,xp,yp,zp,sxz,szx,oyb1(ifacet)
        endif

      ENDDO

      IF(mysize>1) THEN
        CALL MPI_REDUCE(sxxb1,sxxb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(syyb1,syyb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(szzb1,szzb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(sxyb1,sxyb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(sxzb1,sxzb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(syzb1,syzb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(oyb1,oyb,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      ELSE
        sxxb = sxxb1
        syyb = syyb1
        szzb = szzb1
        sxyb = sxyb1
        sxzb = sxzb1
        syzb = syzb1
        oyb = oyb1
      ENDIF

      icall = icall+1

      RETURN

      END
C------------------------------------------------------------------------


c---- function calc_sxx ----------------------N. Beratlis-15 Mar. 2010---
C
C     PURPOSE: Given a point xp,yp,zp inside a grid cell calculate shear 
C     stress Sxx.
C------------------------------------------------------------------------
      real function calc_sxx(u,xp,yp,zp,xu,yc,zc,nx,ny,nz,icyl)
c
      implicit none
c
c.... Input/Output Arrays
      integer nx,ny,nz,icyl
      real    xp,yp,zp
      real    xu(nx),yc(ny),zc(nz)
      real    u(nx,ny,nz)
c
c.... Local Arrays
      integer i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4
      real    b1,c1,dx,dy,dz,theta,u1,u2
c
c.... Functions
      real    anglerad

      calc_sxx = 0.0

      call ijk_xyz(xp,yp,zp,xu,yc,zc,nx,ny,nz,i,j,k,icyl)

      if(icyl==0) then
        theta = yp        
      else
        theta = anglerad(xp,yp)
      endif

      dx = xu(i+1)-xu(i)
      dy = yc(j+1)-yc(j)
      dz = zc(k+1)-zc(k)

      b1 = (theta-yc(j))/dy
      c1 = (zp-zc(k))/dz

      call index_bnd(i+1,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j+1,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i+1,j  ,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k+1,i4,j4,k4,nx,ny,nz)
      u2 = (1-b1)*(1-c1)*u(i1,j1,k1) + b1*(1-c1)*u(i2,j2,k2)
     &   + (1-b1)*  c1  *u(i3,j3,k3) + b1*  c1  *u(i4,j4,k4)

      call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i  ,j+1,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i  ,j+1,k+1,i4,j4,k4,nx,ny,nz)
      u1 = (1-b1)*(1-c1)*u(i1,j1,k1) + b1*(1-c1)*u(i2,j2,k2)
     &   + (1-b1)*  c1  *u(i3,j3,k3) + b1*  c1  *u(i4,j4,k4)

      calc_sxx = (u2-u1)/dx

      return

      end
C------------------------------------------------------------------------


c---- function calc_syy ----------------------N. Beratlis-15 Mar. 2010---
C
C     PURPOSE: Given a point xp,yp,zp inside a grid cell calculate shear 
C     stress Sxx.
C------------------------------------------------------------------------
      real function calc_syy(v,xp,yp,zp,xc,yv,zc,nx,ny,nz,icyl)
c
      implicit none
c
c.... Input/Output Arrays
      integer nx,ny,nz,icyl
      real    xp,yp,zp
      real    xc(nx),yv(ny),zc(nz)
      real    v(nx,ny,nz)
c
c.... Local Arrays
      integer i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4
      real    a1,c1,dx,dy,dz,r,v1,v2
c
c.... Functions
      real    anglerad

      calc_syy = 0.0

      call ijk_xyz(xp,yp,zp,xc,yv,zc,nx,ny,nz,i,j,k,icyl)

      if(icyl==0) then
        r = xp
      else
        r = sqrt(xp**2. + yp**2.)
      endif

      dx = xc(i+1)-xc(i)
      dy = yv(j+1)-yv(j)
      dz = zc(k+1)-zc(k)

      a1 = (r-xc(i))/dx
      c1 = (zp-zc(k))/dz

      call index_bnd(i  ,j+1,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j+1,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j+1,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k+1,i4,j4,k4,nx,ny,nz)
      v2 = (1-a1)*(1-c1)*v(i1,j1,k1) + a1*(1-c1)*v(i2,j2,k2)
     &   + (1-a1)*  c1  *v(i3,j3,k3) + a1*  c1  *v(i4,j4,k4)

      call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j  ,k+1,i4,j4,k4,nx,ny,nz)
      v1 = (1-a1)*(1-c1)*v(i1,j1,k1) + a1*(1-c1)*v(i2,j2,k2)
     &   + (1-a1)*  c1  *v(i3,j3,k3) + a1*  c1  *v(i4,j4,k4)
      if(icyl==0) then
        calc_syy = (v2-v1)/dy
      else
        calc_syy = (v2-v1)/(r*dy)
      endif

      return

      end
C------------------------------------------------------------------------


c---- function calc_szz ----------------------N. Beratlis-15 Mar. 2010---
C
C     PURPOSE: Given a point xp,yp,zp inside a grid cell calculate shear 
C     stress Szz.
C------------------------------------------------------------------------
      real function calc_szz(w,xp,yp,zp,xc,yc,zw,nx,ny,nz,icyl)
c
      implicit none
c
c.... Input/Output Arrays
      integer nx,ny,nz,icyl
      real    xp,yp,zp
      real    xc(nx),yc(ny),zw(nz)
      real    w(nx,ny,nz)
c
c.... Local Arrays
      integer i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4
      real    a1,b1,dx,dy,dz,r,theta,w1,w2
c
c.... Functions
      real    anglerad

      calc_szz=0.0

      call ijk_xyz(xp,yp,zp,xc,yc,zw,nx,ny,nz,i,j,k,icyl)

      if(icyl==0) then
        r = xp
        theta = yp
      else
        r = sqrt(xp**2. + yp**2.)
        theta = anglerad(xp,yp)
      endif

      dx = xc(i+1)-xc(i)
      dy = yc(j+1)-yc(j)
      dz = zw(k+1)-zw(k)

      a1 = (r-xc(i))/dx
      b1 = (theta-yc(j))/dy

      call index_bnd(i  ,j  ,k+1,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j  ,k+1,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j+1,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k+1,i4,j4,k4,nx,ny,nz)
      w2 = (1-a1)*(1-b1)*w(i1,j1,k1) + a1*(1-b1)*w(i2,j2,k2)
     &   + (1-a1)*  b1  *w(i3,j3,k3) + a1*  b1  *w(i4,j4,k4)

      call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j+1,k  ,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k  ,i4,j4,k4,nx,ny,nz)
      w1 = (1-a1)*(1-b1)*w(i1,j1,k1) + a1*(1-b1)*w(i2,j2,k2)
     &   + (1-a1)*  b1  *w(i3,j3,k3) + a1*  b1  *w(i4,j4,k4)

      calc_szz = (w2-w1)/dz

      return

      end
C------------------------------------------------------------------------


c---- function calc_syx ----------------------N. Beratlis-15 Mar. 2010---
C
C     PURPOSE: Given a point xp,yp,zp inside a grid cell calculate shear 
C     stress Syx.
C------------------------------------------------------------------------
      real function calc_syx(v,xp,yp,zp,xc,yv,zc,nx,ny,nz,icyl)
c
      implicit none
c
c.... Input/Output Arrays
      integer nx,ny,nz,icyl
      real    xp,yp,zp
      real    xc(nx),yv(ny),zc(nz)
      real    v(nx,ny,nz)
c
c.... Local Arrays
      integer i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4
      real    b1,c1,dx,dy,dz,r,theta,v1,v2
c
c.... Functions
      real    anglerad

      calc_syx = 0.0

      call ijk_xyz(xp,yp,zp,xc,yv,zc,nx,ny,nz,i,j,k,icyl)

      if(icyl==0) then
        theta = yp        
      else
        r = sqrt(xp**2.+yp**2.)
        theta = anglerad(xp,yp)
      endif

      dx = xc(i+1)-xc(i)
      dy = yv(j+1)-yv(j)
      dz = zc(k+1)-zc(k)

      b1 = (theta-yv(j))/dy
      c1 = (zp-zc(k))/dz

      call index_bnd(i+1,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j+1,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i+1,j  ,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k+1,i4,j4,k4,nx,ny,nz)
      v2 = (1-b1)*(1-c1)*v(i1,j1,k1) + b1*(1-c1)*v(i2,j2,k2)
     &   + (1-b1)*  c1  *v(i3,j3,k3) + b1*  c1  *v(i4,j4,k4)

      call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i  ,j+1,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i  ,j+1,k+1,i4,j4,k4,nx,ny,nz)
      v1 = (1-b1)*(1-c1)*v(i1,j1,k1) + b1*(1-c1)*v(i2,j2,k2)
     &   + (1-b1)*  c1  *v(i3,j3,k3) + b1*  c1  *v(i4,j4,k4)

      if(icyl==0) then
        calc_syx = 0.5*(v2-v1)/dx
      else
        calc_syx = 0.5*r*((v2/xc(i+1)-v1/xc(i))/dx)
      endif

      return

      end
C------------------------------------------------------------------------


c---- function calc_syx ----------------------N. Beratlis-15 Mar. 2010---
C
C     PURPOSE: Given a point xp,yp,zp inside a grid cell calculate shear 
C     stress Syx.
C------------------------------------------------------------------------
      real function calc_sxy(u,xp,yp,zp,xu,yc,zc,nx,ny,nz,icyl)
c
      implicit none
c
c.... Input/Output Arrays
      integer nx,ny,nz,icyl
      real    xp,yp,zp
      real    xu(nx),yc(ny),zc(nz)
      real    u(nx,ny,nz)
c
c.... Local Arrays
      integer i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4
      real    a1,c1,dx,dy,dz,r,u1,u2
c
c.... Functions
      real    anglerad

      calc_sxy = 0.0

      call ijk_xyz(xp,yp,zp,xu,yc,zc,nx,ny,nz,i,j,k,icyl)

      if(icyl==0) then
        r = xp
      else
        r = sqrt(xp**2. + yp**2.)
      endif

      dx = xu(i+1)-xu(i)
      dy = yc(j+1)-yc(j)
      dz = zc(k+1)-zc(k)

      a1 = (r-xu(i))/dx
      c1 = (zp-zc(k))/dz

      call index_bnd(i  ,j+1,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j+1,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j+1,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k+1,i4,j4,k4,nx,ny,nz)
      u2 = (1-a1)*(1-c1)*u(i1,j1,k1) + a1*(1-c1)*u(i2,j2,k2)
     &   + (1-a1)*  c1  *u(i3,j3,k3) + a1*  c1  *u(i4,j4,k4)

      call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j  ,k+1,i4,j4,k4,nx,ny,nz)
      u1 = (1-a1)*(1-c1)*u(i1,j1,k1) + a1*(1-c1)*u(i2,j2,k2)
     &   + (1-a1)*  c1  *u(i3,j3,k3) + a1*  c1  *u(i4,j4,k4)
      if(icyl==0) then
        calc_sxy = 0.5*(u2-u1)/dy
      else
        calc_sxy = 0.5*(u2-u1)/(r*dy)
      endif

      return

      end
C------------------------------------------------------------------------

c---- function calc_szy ----------------------N. Beratlis-15 Mar. 2010---
C
C     PURPOSE: Given a point xp,yp,zp inside a grid cell calculate shear 
C     stress Sxx.
C------------------------------------------------------------------------
      real function calc_szy(w,xp,yp,zp,xc,yc,zw,nx,ny,nz,icyl)
c
      implicit none
c
c.... Input/Output Arrays
      integer nx,ny,nz,icyl
      real    xp,yp,zp
      real    xc(nx),yc(ny),zw(nz)
      real    w(nx,ny,nz)
c
c.... Local Arrays
      integer i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4
      real    a1,c1,dx,dy,dz,r,w1,w2
c
c.... Functions
      real    anglerad

      calc_szy = 0.0

      call ijk_xyz(xp,yp,zp,xc,yc,zw,nx,ny,nz,i,j,k,icyl)

      if(icyl==0) then
        r = xp
      else
        r = sqrt(xp**2. + yp**2.)
      endif

      dx = xc(i+1)-xc(i)
      dy = yc(j+1)-yc(j)
      dz = zw(k+1)-zw(k)

      a1 = (r-xc(i))/dx
      c1 = (zp-zw(k))/dz

      call index_bnd(i  ,j+1,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j+1,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j+1,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k+1,i4,j4,k4,nx,ny,nz)
      w2 = (1-a1)*(1-c1)*w(i1,j1,k1) + a1*(1-c1)*w(i2,j2,k2)
     &   + (1-a1)*  c1  *w(i3,j3,k3) + a1*  c1  *w(i4,j4,k4)

      call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j  ,k+1,i4,j4,k4,nx,ny,nz)
      w1 = (1-a1)*(1-c1)*w(i1,j1,k1) + a1*(1-c1)*w(i2,j2,k2)
     &   + (1-a1)*  c1  *w(i3,j3,k3) + a1*  c1  *w(i4,j4,k4)
      if(icyl==0) then
        calc_szy = 0.5*(w2-w1)/dy
      else
        calc_szy = 0.5*(w2-w1)/(r*dy)
      endif

      return

      end
C------------------------------------------------------------------------

c---- function calc_syz ----------------------N. Beratlis-15 Mar. 2010---
C
C     PURPOSE: Given a point xp,yp,zp inside a grid cell calculate shear 
C     stress Syz.
C------------------------------------------------------------------------
      real function calc_syz(v,xp,yp,zp,xc,yv,zc,nx,ny,nz,icyl)
c
      implicit none
c
c.... Input/Output Arrays
      integer nx,ny,nz,icyl
      real    xp,yp,zp
      real    xc(nx),yv(ny),zc(nz)
      real    v(nx,ny,nz)
c
c.... Local Arrays
      integer i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4
      real    a1,b1,dx,dy,dz,r,theta,v1,v2
c
c.... Functions
      real    anglerad

      calc_syz=0.0

      call ijk_xyz(xp,yp,zp,xc,yv,zc,nx,ny,nz,i,j,k,icyl)

      if(icyl==0) then
        r = xp
        theta = yp
      else
        r = sqrt(xp**2. + yp**2.)
        theta = anglerad(xp,yp)
      endif

      dx = xc(i+1)-xc(i)
      dy = yv(j+1)-yv(j)
      dz = zc(k+1)-zc(k)

      a1 = (r-xc(i))/dx
      b1 = (theta-yv(j))/dy

      call index_bnd(i  ,j  ,k+1,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j  ,k+1,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j+1,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k+1,i4,j4,k4,nx,ny,nz)
      v2 = (1-a1)*(1-b1)*v(i1,j1,k1) + a1*(1-b1)*v(i2,j2,k2)
     &   + (1-a1)*  b1  *v(i3,j3,k3) + a1*  b1  *v(i4,j4,k4)

      call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j+1,k  ,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k  ,i4,j4,k4,nx,ny,nz)
      v1 = (1-a1)*(1-b1)*v(i1,j1,k1) + a1*(1-b1)*v(i2,j2,k2)
     &   + (1-a1)*  b1  *v(i3,j3,k3) + a1*  b1  *v(i4,j4,k4)

      calc_syz = 0.5*(v2-v1)/dz

      return

      end
C------------------------------------------------------------------------


c---- function calc_sxz ----------------------N. Beratlis-15 Mar. 2010---
C
C     PURPOSE: Given a point xp,yp,zp inside a grid cell calculate shear 
C     stress Sxz.
C------------------------------------------------------------------------
      real function calc_sxz(u,xp,yp,zp,xu,yc,zc,nx,ny,nz,icyl)
c
      implicit none
c
c.... Input/Output Arrays
      integer nx,ny,nz,icyl
      real    xp,yp,zp
      real    xu(nx),yc(ny),zc(nz)
      real    u(nx,ny,nz)
c
c.... Local Arrays
      integer i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4
      real    a1,b1,dx,dy,dz,r,theta,u1,u2
c
c.... Functions
      real    anglerad

      calc_sxz=0.0

      call ijk_xyz(xp,yp,zp,xu,yc,zc,nx,ny,nz,i,j,k,icyl)

      if(icyl==0) then
        r = xp
        theta = yp
      else
        r = sqrt(xp**2. + yp**2.)
        theta = anglerad(xp,yp)
      endif

      dx = xu(i+1)-xu(i)
      dy = yc(j+1)-yc(j)
      dz = zc(k+1)-zc(k)

      a1 = (r-xu(i))/dx
      b1 = (theta-yc(j))/dy

      call index_bnd(i  ,j  ,k+1,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j  ,k+1,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j+1,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k+1,i4,j4,k4,nx,ny,nz)
      u2 = (1-a1)*(1-b1)*u(i1,j1,k1) + a1*(1-b1)*u(i2,j2,k2)
     &   + (1-a1)*  b1  *u(i3,j3,k3) + a1*  b1  *u(i4,j4,k4)

      call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j+1,k  ,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k  ,i4,j4,k4,nx,ny,nz)
      u1 = (1-a1)*(1-b1)*u(i1,j1,k1) + a1*(1-b1)*u(i2,j2,k2)
     &   + (1-a1)*  b1  *u(i3,j3,k3) + a1*  b1  *u(i4,j4,k4)

      calc_sxz = 0.5*(u2-u1)/dz

      return

      end
C------------------------------------------------------------------------


c---- function calc_szx ----------------------N. Beratlis-15 Mar. 2010---
C
C     PURPOSE: Given a point xp,yp,zp inside a grid cell calculate shear 
C     stress Szx.
C------------------------------------------------------------------------
      real function calc_szx(w,xp,yp,zp,xc,yc,zw,nx,ny,nz,icyl)
c
      implicit none
c
c.... Input/Output Arrays
      integer nx,ny,nz,icyl
      real    xp,yp,zp
      real    xc(nx),yc(ny),zw(nz)
      real    w(nx,ny,nz)
c
c.... Local Arrays
      integer i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4
      real    b1,c1,dx,dy,dz,r,theta,w1,w2
c
c.... Functions
      real    anglerad

      calc_szx = 0.0

      call ijk_xyz(xp,yp,zp,xc,yc,zw,nx,ny,nz,i,j,k,icyl)

      if(icyl==0) then
        theta = yp        
      else
        r = sqrt(xp**2.+yp**2.)
        theta = anglerad(xp,yp)
      endif

      dx = xc(i+1)-xc(i)
      dy = yc(j+1)-yc(j)
      dz = zw(k+1)-zw(k)

      b1 = (theta-yc(j))/dy
      c1 = (zp-zw(k))/dz

      call index_bnd(i+1,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j+1,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i+1,j  ,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k+1,i4,j4,k4,nx,ny,nz)
      w2 = (1-b1)*(1-c1)*w(i1,j1,k1) + b1*(1-c1)*w(i2,j2,k2)
     &   + (1-b1)*  c1  *w(i3,j3,k3) + b1*  c1  *w(i4,j4,k4)

c      if(abs(xp)>0.495) then
c        write(6,*) '2. szx',xp,yp,zp,b1,c1,w(i1,j1,k1),w(i2,j2,k2),w(i3,j3,k3),w(i4,j4,k4),w2,dx
c      endif

      call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i  ,j+1,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i  ,j+1,k+1,i4,j4,k4,nx,ny,nz)
      w1 = (1-b1)*(1-c1)*w(i1,j1,k1) + b1*(1-c1)*w(i2,j2,k2)
     &   + (1-b1)*  c1  *w(i3,j3,k3) + b1*  c1  *w(i4,j4,k4)

      calc_szx = 0.5*(w2-w1)/dx

      if(abs(xp)>0.495) then
        write(6,*) '1. szx',xp,yp,zp,b1,c1,w(i1,j1,k1),w(i2,j2,k2),w(i3,j3,k3),w(i4,j4,k4),w1,dx,calc_szx
      endif

      return

      end
C------------------------------------------------------------------------


C---- subroutine calc_derv_grid --------------N. Beratlis-15 Mar. 2010---
C
C     PURPOSE: Given a point xp,yp,zp inside grid cell calculate variable
C     derivatives d/dx,d/dy,d/dz using cell values. The derivatives can be
C     in cartesian or cylindrical coordinates.
C------------------------------------------------------------------------
      subroutine calc_derv_grid(u,xp,yp,zp,xu,yu,zu,nx,ny,nz,dudx,dudy,dudz,icyl)
c
      implicit none
c
c...  Input/Output Arrays
      integer nx,ny,nz,icyl
      real    xp,yp,zp,dudx,dudy,dudz,sxx,syy,szz,sxy,sxz,syz
      real    xu(nx),yu(ny),zu(nz)
      real    u(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k
      integer i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4
      real    u1,u2,dx,dy,dz,a1,b1,c1,r,theta

      call ijk_xyz(xp,yp,zp,xu,yu,zu,nx,ny,nz,i,j,k,icyl)

      dx = xu(i+1)-xu(i)
      dy = yu(j+1)-yu(j)
      dz = zu(k+1)-zu(k)

      a1 = (r-xu(i))/dx
      b1 = (theta-yu(j))/dy
      c1 = (zp-zu(k))/dz
c
c.... Compute du/dx
      call index_bnd(i+1,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j+1,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i+1,j  ,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k+1,i4,j4,k4,nx,ny,nz)
      u2 = (1-b1)*(1-c1)*u(i1,j1,k1) + b1*(1-c1)*u(i2,j2,k2)
     &   + (1-b1)*  c1  *u(i3,j3,k3) + b1*  c1  *u(i4,j4,k4)

      call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i  ,j+1,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i  ,j+1,k+1,i4,j4,k4,nx,ny,nz)
      u1 = (1-b1)*(1-c1)*u(i1,j1,k1) + b1*(1-c1)*u(i2,j2,k2)
     &   + (1-b1)*  c1  *u(i3,j3,k3) + b1*  c1  *u(i4,j4,k4)

      dudx = (u2-u1)/dx
c
c.... Compute du/dy
      call index_bnd(i  ,j+1,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j+1,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j+1,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k+1,i4,j4,k4,nx,ny,nz)
      u2 = (1-a1)*(1-c1)*u(i1,j1,k1) + a1*(1-c1)*u(i2,j2,k2)
     &   + (1-a1)*  c1  *u(i3,j3,k3) + a1*  c1  *u(i4,j4,k4)

      call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j  ,k+1,i4,j4,k4,nx,ny,nz)
      u1 = (1-a1)*(1-c1)*u(i1,j1,k1) + a1*(1-c1)*u(i2,j2,k2)
     &   + (1-a1)*  c1  *u(i3,j3,k3) + a1*  c1  *u(i4,j4,k4)
      if(icyl==0) then
        dudy = (u2-u1)/dy
      else
        dudy = (u2-u1)/(r*dy)
      endif
c
c.... Compute du/dz
      call index_bnd(i  ,j  ,k+1,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j  ,k+1,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j+1,k+1,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k+1,i4,j4,k4,nx,ny,nz)
      u2 = (1-a1)*(1-b1)*u(i1,j1,k1) + a1*(1-b1)*u(i2,j2,k2)
     &   + (1-a1)*  b1  *u(i3,j3,k3) + a1*  b1  *u(i4,j4,k4)

      call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
      call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
      call index_bnd(i  ,j+1,k  ,i3,j3,k3,nx,ny,nz)
      call index_bnd(i+1,j+1,k  ,i4,j4,k4,nx,ny,nz)
      u1 = (1-a1)*(1-b1)*u(i1,j1,k1) + a1*(1-b1)*u(i2,j2,k2)
     &   + (1-a1)*  b1  *u(i3,j3,k3) + a1*  b1  *u(i4,j4,k4)
      dudz = (u2-u1)/dz

      return

      end
C------------------------------------------------------------------------

C---- subroutine centerflag-------------------N. Beratlis-02 Jun. 2009---
C
C     PURPOSE: Center flag of velocity
C
C------------------------------------------------------------------------
      subroutine centerflag(flaguo,flaguc,nx,ny,nz,nbd,idir)
c
      implicit none
c
c.... Input/Output Arrays
      integer nx,ny,nz,nbd,idir
      integer flaguo(nx,ny,nz,nbd),flaguc(nx,ny,nz,nbd)
c
c.... Local arrays
      integer i,j,k

      if(idir==1) then

        do i=2,nx
        do j=1,ny
        do k=1,nz
          flaguc(i,j,k,:)=flaguo(i,j,k,:)
        enddo
        enddo
        enddo

      endif

      return

      end
C------------------------------------------------------------------------



C---- subroutine writeflagp-------------------N. Beratlis-04 Jun. 2009---

C------------------------------------------------------------------------
      subroutine writeflagp(stem,flagp,p,xc,yc,zc,nx,ny,nz,nbd,ibd)

      include 'common.h'
      include 'immersed.h'
c
c.... Input/Output Arrays
      integer nx,ny,nz,nbd,ibd
      integer flagp(nx,ny,nz,nbd)
      real    p(nx,ny,nz)
      real    xc(nx),yc(ny),zc(nz)
      character*(*) stem
c
c.... Local Arrays
      integer i,j,k
      real    radius,angle
c
c.... Function
      real    anglerad

      do j=1,ny
c        write(6,*) trim(stem)
        open(unit=10,file=trim(stem)//'.'//index(ibd)//'.'//index(j)
     &        //'.plt',form='formatted')
        write(10,*) 'VARIABLES = "r", "z", "angle", "p"'
        radius = count(flagp(:,j,:,ibd)==-ibd)
        write(10,*) 'ZONE I=',radius,', DATAPACKING=POINT'
        do i=ibmin(ibd),ibmax(ibd)
          do k=jbmin(ibd),kbmax(ibd)
            if(flagp(i,j,k,ibd)==-ibd) then
              angle = anglerad(-zc(k),xc(i))
c              write(6,*) i,j,k,xc(i),zc(k),angle,p(i,j,k)
              write(10,'(5(F14.8,1X))') xc(i),zc(k),angle,p(i,j,k)

c              if(angle>0.16.AND.angle<0.17) then
              if(i==2.AND.k==149) then
                write(6,*) i,j,k,angle,p(i,j,k)
              endif

            endif
          enddo
        enddo
        close(10)
      enddo

      return

      end
C------------------------------------------------------------------------

C---- subroutine writeflagp-------------------N. Beratlis-04 Jun. 2009---

C------------------------------------------------------------------------
      subroutine writeflagu(stem,flagu,u,xc,yc,zc,nx,ny,nz,nbd,ibd)

      include 'common.h'
      include 'immersed.h'
c
c.... Input/Output Arrays
      integer nx,ny,nz,nbd,ibd
      integer flagu(nx,ny,nz,nbd)
      real    u(nx,ny,nz)
      real    xc(nx),yc(ny),zc(nz)
      character*(*) stem
c
c.... Local Arrays
      integer i,j,k
      real    radius,angle
c
c.... Function
      real    anglerad

      do j=1,ny
c        write(6,*) trim(stem)
        open(unit=10,file=trim(stem)//'.'//index(ibd)//'.'//index(j)
     &        //'.plt',form='formatted')
        write(10,*) 'VARIABLES = "r", "z", "angle", "u"'
        radius = count(flagu(:,j,:,ibd)==-ibd)
        write(10,*) 'ZONE I=',radius,', DATAPACKING=POINT'
        do i=ibmin(ibd),ibmax(ibd)
          do k=jbmin(ibd),kbmax(ibd)
            if(flagu(i,j,k,ibd)==-ibd) then
              angle = anglerad(-zc(k),xc(i))
c              write(6,*) i,j,k,xc(i),zc(k),angle,p(i,j,k)
              write(10,'(5(F14.8,1X))') xc(i),zc(k),angle,u(i,j,k)
c              if(angle>0.16.AND.angle<0.17) then
c                write(6,*) i,j,k,angle,p(i,j,k)
c              endif
            endif
          enddo
        enddo
        close(10)
      enddo

      return

      end
C------------------------------------------------------------------------


C---- subroutine write_ubody---------------------N. Beratlis-16 Jun 09---
C
C     PURPOSE: Wirte velocity at surface of body
C
C------------------------------------------------------------------------
      subroutine write_ubody(stem,limu,mimu,u,xu,yu,zu,nx,ny,nz
     &     ,iu,ju,ku,iumtrx,jumtrx,kumtrx,mrku,diru,nxu,nyu,nzu,xnu,ynu,znu)

      include 'common.h'
      include 'immersed.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,limu,mimu
      real    u(nx,ny,nz)
      integer iu(nfcmax),ju(nfcmax),ku(nfcmax)
      integer iumtrx(nsm,nfcmax),jumtrx(nsm,nfcmax),kumtrx(nsm,nfcmax)
      integer mrku(nfcmax),diru(nfcmax)
      real    nxu(nfcmax),nyu(nfcmax),nzu(nfcmax)
      real    xnu(nfcmax),ynu(nfcmax),znu(nfcmax)
      real    xu(nx),yu(ny),zu(nz)
      character*(*) stem
c
c.... Local arrays
      integer i,j,k,im,ii,jj,kk
      real    s1,s2,u1,u2,psi,a,b
c
c.... Functions
      real    anglerad,interp_cellface

      open(unit=10,file=trim(stem)//'.forc.plt',form='formatted')
      write(10,*) 'VARIABLES = "x", "y", "z", "r", "angle", "u"'
      write(10,*) 'ZONE I=',mimu-limu,', DATAPACKING=POINT'

      if(icyl==1) then
        do im=limu+1,mimu

          i = iu(im)
          j = ju(im)
          k = ku(im)

          ii = int(nxu(i))
          jj = int(nyu(i))
          kk = int(nzu(i))

c          write(6,*) im,i,j,k

          s1 = sqrt( (xnu(im) - xu(i)*cos(yu(j)))**2.
     &          + (ynu(im) - xu(i)*sin(yu(j)))**2.
     &          + (znu(im) - zu(k))**2. )
          u1 = u(i,j,k)

          if(mrku(im)==1) then          
            s2 = sqrt( (xnu(im) - nxu(im))**2.
     &               + (ynu(im) - nyu(im))**2.
     &               + (znu(im) - nzu(im))**2.)
            u2 = interp_cellface(nxu(im),nyu(im),nzu(im),xu,yu,zu
     &          ,u,nx,ny,nz,icyl,diru(im))
          elseif(mrku(im)==2) then
            s2 = sqrt( (xnu(im) - xu(iumtrx(1,im))*cos(yu(jumtrx(1,im))))**2.
     &               + (ynu(im) - xu(iumtrx(1,im))*sin(yu(jumtrx(1,im))))**2.
     &               + (znu(im) - zu(kumtrx(1,im)))**2. ) 
            u2 = u(iumtrx(1,im),jumtrx(1,im),kumtrx(1,im))
          endif

          a = (u1-u2)/(s1-s2)
          b = u1-a*s1

c          write(6,*) im,s1,s2,u1,u2,a,b,mrku(im),i,j,k,iu(mtrx(1,im),jumtrx(1,im),kumtrx(1,im)

          psi = anglerad(-zu(k),xu(i))

          write(10,*) xu(i)*cos(yu(j)),xu(i)*sin(yu(j)),zu(k),xu(i),psi,b
 
        enddo
      endif

      close(10)

      return

      end
C------------------------------------------------------------------------



C---- subroutine vec_iext -----------------------N. Beratlis-14 Jan. 09--
C
C     PURPOSE: Given a point q incide a cell and a vector r determine 
C     indices of extension along grid in direction of r. Returns 0 or 1
C
C-----------------------------------------------------------------------
      subroutine vec_iext(q,r,iext,jext,kext,icyl,theta,dtheta)

      IMPLICIT NONE
      REAL, PARAMETER :: pi=3.141592653589793

      INTEGER iext,jext,kext,icyl
      REAL    theta,dtheta
      REAL    q(3),r(3)

      REAL    rp1,thetap1,rp2,thetap2,a,t
      real    x,y,xp,yp,xp2,yp2,psi,psi1,psi2,rp,az,nom,denom
      real    v1(2),v2(2),v3(2)
      integer jext2

      REAL    anglerad

      iext = 0
      jext = 0
      kext = 0

      if(icyl==0) then
        if(r(1)>0) iext=1
        if(r(2)>0) jext=1
      else
        psi = anglerad(q(1),q(2))
        x = q(1)+r(1)
        y = q(2)+r(2)
        call cord_2drot(x,y,psi,xp,yp)
        if(yp>0) jext=1
c
        if(icyl<3) then

          rp = sqrt(q(1)**2. + q(2)**2.)

          psi = theta+0.5*dtheta
          call cord_2drot(q(1),q(2),psi,xp,yp)
          call cord_2drot(q(1)+r(1),q(2)+r(2),psi,xp2,yp2)

          v1(1) = rp*cos(-0.5*dtheta)-xp
          v1(2) = rp*sin(-0.5*dtheta)-yp

          v2(1) = rp*cos(0.5*dtheta)-xp
          v2(2) = rp*sin(0.5*dtheta)-yp
          
          psi1 = anglerad(v1(2),v1(1))
          psi2 = anglerad(v2(2),v2(1))

          v3(1) = xp2-xp
          v3(2) = yp2-yp
            
          psi = anglerad(v3(2),v3(1))

          if(psi>psi1 .AND. psi<psi2) then
            iext=0
          else
            iext=1
          endif

        else

          if(r(3)/=0.0) then
            az = theta/r(3)
            xp = q(1) + az*r(1)
            yp = q(2) + az*r(2)
            rp1 = sqrt(q(1)**2. + q(2)**2.)
            rp2 = sqrt(xp**2. + yp**2.)
            if(rp2>rp1) iext=1
          endif
        endif

      endif

      if(r(3)>0) kext=1

      return
      end
C-----------------------------------------------------------------------



C---- subroutine vec_iext1-----------------------N. Beratlis-14 Jan. 09-
C
C     PURPOSE: Given a point q that lies on a face and a vector r
C     determine indices of extension in direction of r. Returns
C     0 or -1
C
C-----------------------------------------------------------------------
      subroutine vec_iext1(q,r,iext,jext,kext,icyl,theta,dtheta)

      IMPLICIT NONE

      INTEGER iext,jext,kext,icyl
      REAL    theta,dtheta
      REAL    q(3),r(3)

      integer jext2
      REAL    rp1,thetap1,rp2,thetap2,psi,x,y,xp,yp,xp2,yp2,az,psi1,psi2,rp
      REAL    v1(2),v2(2),v3(2)
      REAL    nom,denom,a,t

      REAL    anglerad

      iext = 0
      jext = 0
      kext = 0

      if(icyl==0) then
        if(r(1)<0) iext=-1
        if(r(2)<0) jext=-1
      else
        psi = anglerad(q(1),q(2))
        call cord_2drot(q(1)+r(1),q(2)+r(2),psi,xp,yp)
        if(yp<0) jext=-1

        if(icyl<3) then

          rp = sqrt(q(1)**2. + q(2)**2.)

          psi = theta+0.5*dtheta
          call cord_2drot(q(1),q(2),psi,xp,yp)
          call cord_2drot(q(1)+r(1),q(2)+r(2),psi,xp2,yp2)

          v1(1) = rp*cos(-1.0*dtheta)-xp
          v1(2) = rp*sin(-1.0*dtheta)-yp

          v2(1) = rp*cos(1.0*dtheta)-xp
          v2(2) = rp*sin(1.0*dtheta)-yp
          
          psi1 = anglerad(v1(2),v1(1))
          psi2 = anglerad(v2(2),v2(1))

          v3(1) = xp2-xp
          v3(2) = yp2-yp
            
          psi = anglerad(v3(2),v3(1))

          if(psi>psi1 .AND. psi<psi2) then
            iext=-1
          else
            iext=0
          endif

c            psi = theta+0.5*dtheta

c            call cord_2drot(q(1),q(2),psi,xp,yp)
c            call cord_2drot(q(1)+r(1),q(2)+r(2),psi,xp2,yp2)
c            if(xp2<xp) iext=-1
        else
          if(r(3)/=0.0) then
            az = theta/r(3)
            xp = q(1) + az*r(1)
            yp = q(2) + az*r(2)
            rp1 = sqrt(q(1)**2. + q(2)**2.)
            rp2 = sqrt(xp**2. + yp**2.)
            if(rp2<rp1) iext=-1
          endif
        endif

      endif

      if(r(3)<0) kext=-1

      return
      end
C-----------------------------------------------------------------------


C---- subroutine vec_iext2 ----------------------N. Beratlis-14 Jan. 09-
C
C     PURPOSE: Given a point q that coincides with an eulerian grid point
C     and a vector r determine indices of extension in direction of r.
C     Returns -1 or 1
C
C-----------------------------------------------------------------------
      subroutine vec_iext2(q,r,iext,jext,kext,icyl,theta,dtheta)

      IMPLICIT NONE

      INTEGER iext,jext,kext,icyl
      REAL    theta,dtheta
      REAL    q(3),r(3)

      REAL    rp1,thetap1,rp2,thetap2
      real    x,y,xp,yp,xp2,yp2,psi,az,a,t,nom,denom,psi1,psi2,rp
      real    v1(2),v2(2),v3(2)
      integer jext2

      REAL    anglerad


      iext = 1
      jext = 1
      kext = 1

      if(icyl==0) then
        if(r(1)<0) iext=-1
        if(r(2)<0) jext=-1
      else

        psi = anglerad(q(1),q(2))
        call cord_2drot(q(1)+r(1),q(2)+r(2),psi,xp,yp)
        if(yp<0) jext=-1

        if(icyl<3) then

          rp = sqrt(q(1)**2. + q(2)**2.)
          psi = theta+0.5*dtheta
          call cord_2drot(q(1),q(2),psi,xp,yp)
          call cord_2drot(q(1)+r(1),q(2)+r(2),psi,xp2,yp2)

          v1(1) = rp*cos(dtheta)-xp
          v1(2) = rp*sin(dtheta)-yp

          v2(1) = rp*cos(dtheta)-xp
          v2(2) = rp*sin(dtheta)-yp
          
          psi1 = anglerad(v1(2),v1(1))
          psi2 = anglerad(v2(2),v2(1))
          if(psi2<psi1) psi2 = 4.0*acos(0.0)+psi2

          v3(1) = xp2-xp
          v3(2) = yp2-yp
            
          psi = anglerad(v3(2),v3(1))

          if(psi>psi1 .AND. psi<psi2) then
            iext=0
          else
            iext=1
          endif

c          psi = theta+0.5*dtheta

c          call cord_2drot(q(1),q(2),psi,xp,yp)
c          call cord_2drot(q(1)+r(1),q(2)+r(2),psi,xp2,yp2)
c          if(xp2<xp) iext=-1
        else
          az = theta/r(3)
          xp = q(1) + az*r(1)
          yp = q(2) + az*r(2)
          rp1 = sqrt(q(1)**2. + q(2)**2.)
          rp2 = sqrt(xp**2. + yp**2.)
          if(rp2<rp1) iext=-1
        endif
      endif

      if(r(3)<0) kext=-1

      return
      end
C-----------------------------------------------------------------------


C---- subroutine per_index ------------------N. Beratlis-Jun. 11 2009---
C
C     PURPOSE: Corresponding periodic index when index is out of bounds
C
C-----------------------------------------------------------------------
C
      subroutine per_index(j1,ny,j2)

      implicit none
c
c.... Input/Output
      integer j1,ny,j2

      if(j1<1) then
        j2 = ny-2+j1
      elseif(j1>ny) then
        j2 = j1-ny+2
      endif

      return

      end
C-----------------------------------------------------------------------

      
C---- subroutine add_random_noise -----------N. Beratlis-15 Jun. 2009---
C
C     PURPOSE: Add random noise to field variable u(nx,ny,nz)
C
C-----------------------------------------------------------------------
      subroutine add_random_noise(u,nx,ny,nz,m)
c
c      USE IFPORT
      implicit none
c
c.... Input/Output arrays      
      integer nx,ny,nz
      real    m
      real    u(nx,ny,nz)
C
C.... Local arrays
      integer i,j,k
c
c.... Function
c      real    rand   
      integer iseed   !!!!!!
      real x   !!!!!!
      write(*,*) "add_random_noise"
      do i=1,nx
      do j=1,ny
      do k=1,nz
        CALL RANDOM_SEED(iseed)   !!!!!!
        CALL RANDOM_NUMBER(x)  !!!!!!
        u(i,j,k) = u(i,j,k) + m*(x-0.5)  !!!!!!
!!!!!!        u(i,j,k) = u(i,j,k) + m*(rand()-0.5)
c        u(i,j,k) = u(i,j,k) + m*(ranf()-0.5)
      enddo
      enddo
      enddo

      return

      end
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------

      
C---- subroutine add_random_noise_2 --------------
C
C     PURPOSE: Add random noise to field variable u(nx,ny,nz)
C
C-----------------------------------------------------------------------
      subroutine add_random_noise_2(u,nx,ny,nz,tag,loc,mag)
c
c      USE IFPORT
      implicit none
c
c.... Input/Output arrays      
      integer nx,ny,nz
      integer tag
      real    u(nx,ny,nz)
      real    m(nx,ny,nz),loc(nx)
      real    mag
C
C.... Local arrays
      integer i,j,k
c
c.... Function
c      real    rand   
      integer iseed   !!!!!!
      real x   !!!!!!

      if(tag == 1) then
c        parabolic
c        zero at centerline
         do i=1,nx
            m(i,:,:) = (loc(i)-1.)**2
         enddo

      do i=1,nx
      do j=1,ny
      do k=1,nz
        CALL RANDOM_SEED(iseed)   !!!!!!
        CALL RANDOM_NUMBER(x)  !!!!!!
        u(i,j,k) = u(i,j,k) + mag*m(i,j,k)*(x-0.5)  !!!!!!
      enddo
      enddo
      enddo

      endif

      return

      end


c---- subroutine iopressrad------------------N. Beratlis-30 Jun. 2009---
C
      subroutine iopressrad(filename,p,nx,ny,nz,xc,zc,jsl,r,n)
C
C     PURPOSE: Interpolate pressure along a circle of radius r and write
C     to file
C
C-----------------------------------------------------------------------
      implicit none
c
      REAL, PARAMETER :: pi=3.141592653589793
c
c.... Input/Output Array
      integer nx,ny,nz,n,jsl
      real    r
      real    xc(nx),zc(nz)
      real    p(nx,ny,nz)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,ip
      real    xp,zp,pint,dtheta,theta
      real    p2d(nx,nz)
c
c.... Functions
      real    interp2d,anglerad

      dtheta = pi/real(n-1)

      p2d(:,:) = p(:,jsl,:)

      open(unit=10,file=trim(filename),form='formatted')
c      write(10,'(A,A)') 'VARIABLES = "ANGLE", "P"'
c      write(10,'(A,I8,A)') 'ZONE I=',n,', DATAPACKING=POINT'
      do ip=1,n
        theta = dtheta*real(ip-1)
        xp = r*sin(theta)
        zp = r*cos(theta)

        pint = interp2d(xp,zp,xc,zc,p2d,nx,nz)
c`        write(6,*) xp,zp,theta,pint
        write(10,'(12(1X,F14.8))') pi-theta,pint
      enddo

      close(10)
      return

      end
C-----------------------------------------------------------------------





C---- subroutine vec_ijkext --------------------N. Beratlis-14 Jan. 09--
C
C     PURPOSE: Given a point q incide a cell and a vector r determine 
C     indices of extension along grid in direction of r.
C
C-----------------------------------------------------------------------
      subroutine vec_ijkext(q,r,iext,jext,kext,icyl,dir,theta,dtheta)

      IMPLICIT NONE
      REAL, PARAMETER :: pi=3.141592653589793

      INTEGER iext,jext,kext,icyl,dir
      REAL    theta,dtheta
      REAL    q(3),r(3)

      REAL    rp1,thetap1,rp2,thetap2,a,b,c,d,det,R2,t1,t2,t
      real    x,y,xp,yp,xp2,yp2,psi,psi1,psi2,rp,az,nom,denom
      real    v1(2),v2(2),v3(2)
      integer jext2

      REAL    anglerad

      iext = 0
      jext = 0
      kext = 0

      if(icyl==0) then
        if(r(1)>0) iext=1
        if(r(2)>0) jext=1
      else
        psi = anglerad(q(1),q(2))
        x = q(1)+r(1)
        y = q(2)+r(2)
        call cord_2drot(x,y,psi,xp,yp)
        if(yp>0) jext=1
c
        if(dir<3) then

          iext=1

          R2 = theta**2.
          a = r(1)**2.+r(2)**2.
          b = 2.*(q(1)*r(1)+q(2)*r(2))
          c = q(1)**2.0 + q(2)**2.- R2

          det = (b**2.0) - 4.0*a*c
          if(det>0.0) then
            t1 = (-b+sqrt(det))/(2*a)        
            t2 = (-b-sqrt(det))/(2*a)
            if(t1>0.0 .OR. t2>0.0) iext=0
          endif

        else

          if(r(3)/=0.0) then
            az = theta/r(3)
            xp = q(1) + az*r(1)
            yp = q(2) + az*r(2)
            rp1 = sqrt(q(1)**2. + q(2)**2.)
            rp2 = sqrt(xp**2. + yp**2.)
            if(rp2>rp1) iext=1
          endif
        endif

      endif

      if(r(3)>0) kext=1

c      if(dir==1) then
c        jext=0
c        kext=0
c      elseif(dir==2) then
c        iext=0
c        kext=0
c      elseif(dir==3) then
c        iext=0
c        jext=0
c      endif

      return
      end
C-----------------------------------------------------------------------





C---- subroutine vec_ijkext_face ----------------N. Beratlis-14 Jan. 09-
C
C     PURPOSE: Given a point q that lies on a face and a vector r
C     determine indices of extension in direction of r.
C
C-----------------------------------------------------------------------
      subroutine vec_ijkext_face(q,r,iext,jext,kext,icyl,dir,theta,dtheta,face)

      IMPLICIT NONE

      INTEGER iext,jext,kext,dir,icyl,face
      REAL    theta,dtheta
      REAL    q(3),r(3)

      integer jext2
      REAL    rp1,thetap1,rp2,thetap2,psi,x,y,xp,yp,xp2,yp2,az,psi1,psi2,rp
      real    a,b,c,d,det,t1,t2,R2
      REAL    v1(2),v2(2),v3(2)
      REAL    nom,denom,t

      REAL    anglerad

      iext = 0
      jext = 0
      kext = 0

      if(icyl==0) then
        if(r(1)>0) iext=1
        if(r(2)>0) jext=1
      else
        psi = anglerad(q(1),q(2))
        call cord_2drot(q(1)+r(1),q(2)+r(2),psi,xp,yp)
        if(yp>0) jext=1

        if(dir==1) then

          iext=1

          R2 = theta**2.
          a = r(1)**2.+r(2)**2.
          b = 2.*(q(1)*r(1)+q(2)*r(2))
          c = q(1)**2.0 + q(2)**2.- R2

          det = (b**2.0) - 4.0*a*c
          if(det>0.0) then
            t1 = (-b+sqrt(det))/(2*a)        
            t2 = (-b-sqrt(det))/(2*a)
c            write(6,*) 't1=',t1,',t2=',t2
            if(t1>0.0 .OR. t2>0.0) iext=0
          endif

        elseif(dir==2) then

          if(face==1) then
            iext=1
            a=tan(theta)
            if(jext==1) a = tan(theta+dtheta)

            d = r(2)-a*r(1)
            if(d/=0.0) then
              t = (a*q(1)-q(2))/d
              if(t>0.0) then
                b = q(1)+t*r(1)
                c = q(2)+t*r(2)
                rp = sqrt(b**2. + c**2.)
                if(rp<sqrt(q(1)**2.+q(2)**2.)) iext=0
              endif
            endif
          endif

c          rp = sqrt(q(1)**2. + q(2)**2.)
c          if(jext==1) then 
c            psi = theta+0.5*dtheta
c          else
c            psi = theta-0.5*dtheta
c          endif
c          call cord_2drot(q(1),q(2),psi,xp,yp)
c          call cord_2drot(q(1)+r(1),q(2)+r(2),psi,xp2,yp2)
c
c          if(xp2>xp) then
c            iext=1
c          endif

        else
          if(r(3)/=0.0) then
            az = theta/r(3)
            xp = q(1) + az*r(1)
            yp = q(2) + az*r(2)
            rp1 = sqrt(q(1)**2. + q(2)**2.)
            rp2 = sqrt(xp**2. + yp**2.)
            if(rp2>rp1) iext=1
          endif
        endif

      endif

      if(r(3)>0) kext=1

      if(face==1) then
        if(dir==1) then
          jext=0
          kext=0
          if(iext==0) iext=-1
        elseif(dir==2) then
          if(iext==0) then
            iext=-1
          else
            iext=0
          endif
          kext=0
        else
          if(iext==0) then
            iext=-1
          else
            iext=0
          endif
          jext=0
        endif

      elseif(face==2) then
        if(dir==1) then
          if(jext==0) then
            jext=-1
          else
            jext=0
          endif
          kext=0
        elseif(dir==2) then
          iext=0
          kext=0
          if(jext==0) jext=-1
        else
          iext=0
          if(jext==0) then
            jext=-1
          else
            jext=0
          endif
        endif

      else !face=3
        if(dir==1) then
          jext=0
          if(kext==0) then
            kext=-1
          else
            kext=0
          endif
        elseif(dir==2) then
          iext=0
          if(kext==0) then
            kext=-1
          else
            kext=0
          endif
        else
          iext=0
          jext=0
          if(kext==0) kext=-1
        endif
      endif


      return
      end
C-----------------------------------------------------------------------



C---- subroutine vec_ijkext_gridpoint -----------N. Beratlis-14 Jan. 09-
C
C     PURPOSE: Given a point q that coincides with an eulerian grid point
C     and a vector r determine indices of extension in direction of r.
C     Returns -1 or 1
C
C-----------------------------------------------------------------------
      subroutine vec_ijkext_gridpoint(q,r,iext,jext,kext,icyl,dir,theta,dtheta)

      IMPLICIT NONE

      INTEGER iext,jext,kext,icyl,dir
      REAL    theta,dtheta
      REAL    q(3),r(3)

      REAL    rp1,thetap1,rp2,thetap2,a,b,c,d,det,t1,t2,R2
      real    x,y,xp,yp,xp2,yp2,psi,az,t,nom,denom,psi1,psi2,rp
      real    v1(2),v2(2),v3(2)
      integer jext2

      REAL    anglerad


      iext = 1
      jext = 1
      kext = 1

      if(icyl==0) then
        if(r(1)<0) iext=0
        if(r(2)<0) jext=0
      else

        psi = anglerad(q(1),q(2))
        call cord_2drot(q(1)+r(1),q(2)+r(2),psi,xp,yp)
        if(yp<0.0) jext=0

        if(dir==1) then
          R2 = theta**2.
          a = r(1)**2.+r(2)**2.
          b = 2.*(q(1)*r(1)+q(2)*r(2))
          c = q(1)**2.0 + q(2)**2.- R2

          det = (b**2.0) - 4.0*a*c
          if(det>0.0) then
            t1 = (-b+sqrt(det))/(2*a)        
            t2 = (-b-sqrt(det))/(2*a)
            if(t1>0.0 .OR. t2>0.0) iext=0
          endif
        elseif(dir==2) then

          a=tan(theta-dtheta)
          if(jext==1) a = tan(theta+dtheta)

          d = r(2)-a*r(1)
          if(d/=0.0) then
            t = (a*q(1)-q(2))/d
            if(t>0.0) then
              b = q(1)+t*r(1)
              c = q(2)+t*r(2)
              rp = sqrt(b**2. + c**2.)
              if(rp<sqrt(q(1)**2.+q(2)**2.)) iext=0
            endif
          endif

        else
          az = theta/r(3)
          xp = q(1) + az*r(1)
          yp = q(2) + az*r(2)
          rp1 = sqrt(q(1)**2. + q(2)**2.)
          rp2 = sqrt(xp**2. + yp**2.)
          if(rp2<rp1) iext=0
        endif
      endif

      if(r(3)<0) kext=0

      if(dir==1) then
        if(iext==0) iext=-1

        if(jext==1) then
          jext=0
        else
          jext=-1
        endif

        if(kext==1) then
          kext=0
        else
          kext=-1
        endif

      elseif(dir==2) then
        if(iext==1) then
          iext=0
        else
          iext=-1
        endif

        if(jext==0) jext=-1

        if(kext==1) then
          kext=0
        else
          kext=-1
        endif
      else !dir=3
        if(iext==1) then
          iext=0
        else
          iext=-1
        endif

        if(jext==1) then
          jext=0
        else
          jext=-1
        endif

        if(kext==0) kext=-1
      endif

      return
      end
C-----------------------------------------------------------------------


C---- function vel_body2grid-----------------N. Beratlis-14 Jan. 2009--
C
C     PURPOSE: From a point q on the immersed boundary surface find 1
C     intersection along grilines and interpolate value of uo
C     at the intersection
C
C-----------------------------------------------------------------------
      integer function vel_body2grid(q,xu_car,yu_car,xu,yu,zu,zug,flaguo
     &     ,uo,nx,ny,nz,nzg,xint,yint,zint,uint,n1,nbd,myrank,icyl,dir1,ifacet)
c
      IMPLICIT NONE
      REAL, PARAMETER :: pi=3.141592653589793
c
c.... Input/Output Arrays      
      INTEGER nx,ny,nz,nzg,dir,nbd,n1,icyl,myrank,ifacet,dir1
      REAL    xint,yint,zint,uint
      INTEGER flaguo(nx,ny,nz,nbd)
      REAL    q(3)
      REAL    uo(nx,ny,nz)
      REAL    xu_car(nx,ny),yu_car(nx,ny),xu(nx),yu(ny),zu(nz),zug(nzg)
c
c.... Local Arrays
      INTEGER i,j,k,kg,iext,jext,kext,iext1,jext1,kext1,next,sgn
      INTEGER intrs,ii,iphy,iflagt,iflag
      REAL    xp,yp,zp,r
      INTEGER nclock
      real    tclock,clocktemp
      REAL    clock(10)
c
c.... Functions      
      INTEGER ray_face_int,physicalface,grid2grid_intr
      REAL    interp_cellface,extmag,anglerad

      
      nclock = 10

      xp = q(1)
      yp = q(2)
      zp = q(3)
      
      iphy = 0
      intrs = 1

      sgn = isign(1,dir1)
      dir = abs(dir1)

      call ijk_xyz(xp,yp,zp,xu,yu,zug,nx,ny,nzg,i,j,kg,icyl)

      if(dir==1) then

        r = sqrt(xp**2.0 + yp**2.0)

        if(r==0.0) then
          vel_body2grid=0
          return
        endif

        k = kg - myrank*(nz-2)

        if(k>=1 .AND. k<nz-1) then

          iext = 0
          iext1 = -1
          if(sgn==1) then
            iext=1
            iext1=1
          endif
          
          zint = zp

          i = i+iext
          if(icyl==1) call index_bnd(i,j,k,i,j,k,nx,ny,nz)
          iflag = physicalface(flaguo,nx,ny,nz,i,j,k,nbd,dir)
          if(iflag/=0) then
            iphy = iphy+1
            if(icyl==0) then
              xint = xu(i)
              yint = yp
            else
              r = sqrt(xp**2.0 + yp**2.0)
              xint = abs(xu(i))*xp/r
              yint = abs(xu(i))*yp/r
            endif
            uint = interp_cellface(xint,yint,zint,xu,yu,zu,uo
     &                ,nx,ny,nz,icyl,dir)
          endif

          if(iphy==0) then
            iflagt=0
            do while(iflagt==0.AND.intrs<n1)
              intrs = intrs+1
              i = i+iext1
              if(icyl==1) call index_bnd(i,j,k,i,j,k,nx,ny,nz)
              if(physicalface(flaguo,nx,ny,nz,i,j,k,nbd,dir)/=0) then
                iphy = iphy+1
                iflagt=1
                if(icyl==0) then
                  xint = xu(i)
                  yint = yp
                else
                  r = sqrt(xp**2.0 + yp**2.0)
                  xint = abs(xu(i))*xp/r
                  yint = abs(xu(i))*yp/r
                endif
                uint = interp_cellface(xint,yint,zint,xu,yu,zu,uo
     &               ,nx,ny,nz,icyl,dir)

              endif
              if(intrs>=n1) iflagt=1
            enddo
          endif
        endif

      elseif(dir==2) then

        r = sqrt(xp**2.0 + yp**2.)
        if(r==0.0) then
          vel_body2grid=0
          return
        endif

        k = kg - myrank*(nz-2)

        if(k>=1 .AND. k<nz-1) then

          jext = 0
          jext1 = -1
          if(sgn==1) then
            jext=1
            jext1=1
          endif
          
          zint = zp

          j = j+jext
          if(icyl==1) call index_bnd(i,j,k,i,j,k,nx,ny,nz)
          iflag = physicalface(flaguo,nx,ny,nz,i,j,k,nbd,dir)

          if(iflag/=0) then
            iphy = iphy+1            
            if(icyl==0) then
              xint = xp
              yint = yu(j)
            else
              r = sqrt(xp**2.0 + yp**2.0)
              xint = r*cos(yu(j))
              yint = r*sin(yu(j))
            endif
            uint = interp_cellface(xint,yint,zint,xu,yu,zu,uo
     &                ,nx,ny,nz,icyl,dir)
          endif

          if(iphy==0) then
            iflagt=0
            do while(iflagt==0.AND.intrs<n1)
              intrs = intrs+1
              j = j+jext1
              if(icyl==1) call index_bnd(i,j,k,i,j,k,nx,ny,nz)
              if(physicalface(flaguo,nx,ny,nz,i,j,k,nbd,dir)/=0) then
                iphy = iphy+1
                iflagt=1
                if(icyl==0) then
                  xint = xp
                  yint = yu(j)
                else
                  r = sqrt(xp**2.0 + yp**2.0)
                  xint = r*cos(yu(j))
                  yint = r*sin(yu(j))
                endif
                uint = interp_cellface(xint,yint,zint,xu,yu,zu,uo
     &               ,nx,ny,nz,icyl,dir)
              endif
              if(intrs>=n1) iflagt=1
            enddo
          endif
        endif

      elseif(dir==3) then

        kext = 0
        kext1 = -1
        if(sgn==1) then
          kext=1
          kext1=1
        endif
          
        xint = xp
        yint = yp

        kg = kg + kext
        k  = kg - myrank*(nz-2)
        if(k>=1 .AND. k<nz-1) then

          iflag = physicalface(flaguo,nx,ny,nz,i,j,k,nbd,dir)
          if(iflag/=0) then
            iphy = iphy+1            
            zint = zu(k)
            uint = interp_cellface(xint,yint,zint,xu,yu,zu,uo
     &                ,nx,ny,nz,icyl,dir)
          endif
        endif


        if(iphy==0) then
          iflagt=0
          do while(iflagt==0.AND.intrs<n1)
            intrs = intrs+1
            kg = kg+kext1
            k = kg - myrank*(nz-2)
            if(k>=1 .AND. k<nz-1) then
              if(physicalface(flaguo,nx,ny,nz,i,j,k,nbd,dir)/=0) then
                iphy = iphy+1
                iflagt=1
                zint = zu(k)
                uint = interp_cellface(xint,yint,zint,xu,yu,zu,uo
     &               ,nx,ny,nz,icyl,dir)
              endif
            endif
            if(intrs>=n1) iflagt=1
          enddo
        endif

      endif

      if(iphy==1) then
        vel_body2grid = intrs
      else
        vel_body2grid = n1+1
      endif

      return

      end
C------------------------------------------------------------------------





C---- subroutine press_body1 ------------------N. Beratlis-26 Dec. 2008-
C
C     PURPOSE: Calculate pressure on surface on immersed body by inter-
C     polating the value of pressure on the center of triangles using
C     the nodes of the cell that contains the vertex center
C     WARNING: Works assuming flagpbdo and flagpbdi have been corrected.
C
C-----------------------------------------------------------------------
      subroutine press_body1(p,flagp,xc,yc,zc,xc_car,yc_car,nx,ny,nz,nbd
     &     ,pbd,mrkpb,vertexc,unvect,nfacet,ilb,ile,ibd,t,dt)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c      
c.... Input/Output arrays
      integer nx,ny,nz,nfacet,ilb,ile,ibd,nbd
      real    t,dt
      integer mrkpb(nfacet)
      integer flagp(nx,ny,nz,nbd)
      real    xc(nx),yc(ny),zc(nz)
      real    xc_car(nx,ny),yc_car(nx,ny)
      real    p(nx,ny,nz)
      real    vertexc(3,nfacet),unvect(3,nfacet),pbd(nfacet)
c
c.... Local arrays
      integer i,j,k,ii,nb,ib,nzg,n1,ifacet,nvld,nnvld,nvldg,nnvldg,n1max
      real    xp,yp,zp,xint,yint,zint,pint,s,amin,psi
      real    q(3),rvec(3)
      logical cond

      integer, dimension(:), allocatable :: mrkpb1,mrkpbg1
      real, dimension(:), allocatable :: pbd1,zcg
c
c.... Functions
      integer vel_body2fluid
      real    interp3d,dpdn,anglerad

c      write(6,*) 'inside press_body1'
      
      n1max = -1
      n1 = 8
      amin = 0.0
      nnvld = 0
      nvld = 0

      ALLOCATE(mrkpb1(mb(ibd)),mrkpbg1(mb(ibd)))
      ALLOCATE(pbd1(mb(ibd)))
      mrkpb1 = n1+1
      mrkpbg1 = n1+1

      nzg = mysize*(nz-2)+2
      ALLOCATE(zcg(nzg))
      call calc_zg(zc,zcg,nz,nzg)

      do ii=ilb,ile

        xp = vertexc(1,ii)
        yp = vertexc(2,ii)
        zp = vertexc(3,ii)

        q(1) = xp
        q(2) = yp
        q(3) = zp

        rvec(1) = unvect(1,ii)
        rvec(2) = unvect(2,ii)
        rvec(3) = unvect(3,ii)

        if(icyl==0) then
          cond = yp>ymin .AND. yp<ymax
        else
          cond = .true.
        endif

        if(cond) then
          mrkpb1(ii) = vel_body2fluid(q,rvec,xc_car,yc_car,xc,yc,zc,zcg,flagp,p
     &     ,nx,ny,nz,nzg,xint,yint,zint,pint,n1,nbd,myrank,icyl,ifacet,amin)

          if(mrkpb1(ii)>n1max) n1max=mrkpb1(ii)

          if(mrkpb1(ii)<=n1) then
            s = sqrt( (q(1)-xint)**2. + (q(2)-yint)**2. + (q(3)-zint)**2.)
            pbd1(ii) = pint+s*dpdn(xp,yp,zp,unvect(1,ii),unvect(2,ii),unvect(3,ii),t,dt,ibd)
          endif
        endif
      enddo

c      write(6,*) 'n1max=',n1max

      IF(mysize>1) THEN
        CALL MPI_ALLREDUCE(mrkpb1,mrkpbg1,mb(ibd),MPI_INTEGER,MPI_MIN,MPI_COMM_EDDY,IERR)
      ELSE
        mrkpbg1 = mrkpb1
      ENDIF

      do ifacet=ilb,ile
        ii = ifacet-ilb+1

        if(mrkpb1(ii)<=mrkpbg1(ii) .AND. mrkpb1(ii)<=n1) then           
          mrkpb1(ii) = 1
          nvld = nvld+1
        else
          mrkpb1(ii) = 0
          pbd1(ii) = 0.0
        endif

      enddo

      if(mysize>1) then
        call mpi_reduce(pbd1,pbd(ilb),mb(ibd),MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        call mpi_reduce(mrkpb1,mrkpb(ilb),mb(ibd),MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        call mpi_reduce(nvld,nvldg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      else
        pbd(ilb:ile) = pbd1(1:mb(ibd))
        mrkpb(ilb:ile) = mrkpb1(1:mb(ibd))
        nvldg = nvld
      endif

      nnvldg = ile-ilb+1-nvldg

      DEALLOCATE(pbd1,mrkpb1,mrkpbg1)

      if(myrank==0) then
        if(nnvldg>0) then
          open(unit=16,file='stats_imb.dat',position='append'
     &        ,form='formatted')
          write(16,*) 'Number of valid pressure triangles',nvldg
          write(16,*) 'Number of non-valid pressure triangles',nnvldg
          close(16)
        endif
      endif

      return

      end
C-----------------------------------------------------------------------


C---- subroutine press_body2 ------------------N. Beratlis-26 Dec. 2008-
C
C     PURPOSE: Calculate pressure on surface on immersed body by inter-
C     polating the value of pressure on the center of triangles using
C     the nodes of the cell that contains the vertex center
C     WARNING: Works assuming flagpbdo and flagpbdi have been corrected.
C
C-----------------------------------------------------------------------
      subroutine press_body2(p,flagp,xc,yc,zc,xc_car,yc_car,nx,ny,nz,nbd
     &     ,pbd,mrkpb,vertexc,unvect,nfacet,ilb,ile,ibd,t,dt)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c      
c.... Input/Output arrays
      integer nx,ny,nz,nfacet,ilb,ile,ibd,nbd
      real    t,dt
      integer mrkpb(nfacet)
      integer flagp(nx,ny,nz,nbd)
      real    xc(nx),yc(ny),zc(nz)
      real    xc_car(nx,ny),yc_car(nx,ny)
      real    p(nx,ny,nz)
      real    vertexc(3,nfacet),unvect(3,nfacet),pbd(nfacet)
c
c.... Local arrays
      integer i,j,k,ii,nb,ib,nzg,n1,ifacet,nvld,nnvld,nvldg,nnvldg,n1max,ksign
      real    xp,yp,zp,xint,yint,zint,pint,s,amin,psi,ds
      real    q(3),rvec(3)
      logical cond

      integer, dimension(:), allocatable :: mrkpb1,mrkpbg1
      real, dimension(:), allocatable :: pbd1,zcg
c
c.... Functions
      integer vel_body2fluid,body2grid2,body2fluid2
      real    interp3d,dpdn,anglerad,mindxdydz

c      write(6,*) 'inside press_body1'
      
      n1max = -1
      n1 = 8
      amin = 0.0
      nnvld = 0
      nvld = 0

      ALLOCATE(mrkpb1(mb(ibd)),mrkpbg1(mb(ibd)))
      ALLOCATE(pbd1(mb(ibd)))
      mrkpb1 = n1+1
      mrkpbg1 = n1+1

      nzg = mysize*(nz-2)+2
      ALLOCATE(zcg(nzg))
      call calc_zg(zc,zcg,nz,nzg)

      do ii=ilb,ile

        xp = vertexc(1,ii)
        yp = vertexc(2,ii)
        zp = vertexc(3,ii)

        q(1) = xp
        q(2) = yp
        q(3) = zp

        rvec(1) = unvect(1,ii)
        rvec(2) = unvect(2,ii)
        rvec(3) = unvect(3,ii)

        ksign = 0
        if(rvec(3)>0) ksign=1

        if(icyl==0) then
          cond = yp>=ymin.AND.yp<=ymax .AND. zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
        else
          cond = zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
        endif

        if(cond) then
          mrkpb1(ii) = body2fluid2(q,rvec,xc_car,yc_car,xc,yc,zc
     &         ,p,flagp,nx,ny,nz,xint,yint,zint,pint,n1,nbd,icyl,ifacet,amin)

          if(mrkpb1(ii)>n1max) n1max=mrkpb1(ii)

          if(mrkpb1(ii)<=n1) then
            s = sqrt( (q(1)-xint)**2. + (q(2)-yint)**2. + (q(3)-zint)**2.)
            pbd1(ii) = pint+s*dpdn(xp,yp,zp,unvect(1,ii),unvect(2,ii),unvect(3,ii),t,dt,ibd)
          endif

        endif
      enddo



      IF(mysize>1) THEN
        CALL MPI_ALLREDUCE(mrkpb1,mrkpbg1,mb(ibd),MPI_INTEGER,MPI_MIN,MPI_COMM_EDDY,IERR)
      ELSE
        mrkpbg1 = mrkpb1
      ENDIF

      do ifacet=ilb,ile
        ii = ifacet-ilb+1

        if(mrkpb1(ii)<=mrkpbg1(ii) .AND. mrkpb1(ii)<=n1) then           
          mrkpb1(ii) = 1
          nvld = nvld+1
        else
          mrkpb1(ii) = 0
          pbd1(ii) = 0.0
        endif

      enddo

      if(mysize>1) then
        call mpi_reduce(pbd1,pbd(ilb),mb(ibd),MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        call mpi_reduce(mrkpb1,mrkpb(ilb),mb(ibd),MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        call mpi_reduce(nvld,nvldg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      else
        pbd(ilb:ile) = pbd1(1:mb(ibd))
        mrkpb(ilb:ile) = mrkpb1(1:mb(ibd))
        nvldg = nvld
      endif

      nnvldg = ile-ilb+1-nvldg

      DEALLOCATE(pbd1,mrkpb1,mrkpbg1)

      if(myrank==0) then
        if(nnvldg>0) then
          open(unit=16,file='stats_imb.dat',position='append'
     &        ,form='formatted')
          write(16,*) 'Number of valid pressure triangles',nvldg
          write(16,*) 'Number of non-valid pressure triangles',nnvldg
          close(16)
        endif
      endif

      return

      end
C-----------------------------------------------------------------------


C---- subroutine press_body3 ------------------N. Beratlis-26 Dec. 2008-
C
C     PURPOSE: Calculate pressure on surface on immersed body by inter-
C     polating the value of pressure on the center of triangles using
C     the nodes of the cell that contains the vertex center
C     WARNING: Works assuming flagpbdo and flagpbdi have been corrected.
C
C-----------------------------------------------------------------------
      subroutine press_body3(p,flagp,xc,yc,zc,xc_car,yc_car,nx,ny,nz,nbd
     &     ,pbd,mrkpb,vertexc,unvect,nfacet,ibd,mbd,ilb,ile,t,dt)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c      
c.... Input/Output arrays
      integer nx,ny,nz,nfacet,ibd,nbd,mbd,ilb,ile
      real    t,dt
      integer mrkpb(nfacet)
      integer flagp(nx,ny,nz,nbd)
      real    xc(nx),yc(ny),zc(nz)
      real    xc_car(nx,ny),yc_car(nx,ny)
      real    p(nx,ny,nz)
      real    vertexc(3,nfacet),unvect(3,nfacet),pbd(nfacet)
c
c.... Local arrays
      integer imov,n
      integer i,j,k,ii,nb,ib,nzg,n1,ifacet,nvld,nnvld,nvldg,nnvldg,n1max,ksign
      real    xp,yp,zp,xint,yint,zint,pint,s,amin,psi,ds
      real    q(3),rvec(3)
      logical cond

      integer, dimension(:), allocatable :: mrkpb1,mrkpbg1
      real, dimension(:), allocatable :: pbd1
c
c.... Functions
      integer vel_body2fluid,body2grid2,body2fluid2
      real    interp3d,dpdn,anglerad,mindxdydz

      
      n1max = -1
      n1 = 8
      amin = 0.0
      nnvld = 0
      nnvldg = 0
      nvld = 0

      mrkpb(ilb:ile)=0    !!!!!! instead of mrkpb=0

c      IF(IVRTX==0) THEN
c        ilb = lb(ibd)+1
c        ile = lb(ibd)+mb(ibd)
c      ELSE
c        ilb = lv(ibd)+1
c        ile = lv(ibd)+mv(ibd)
c      ENDIF

      imov=0
      if(ibd>=mbd) imov=1

      do ii=ilb,ile
  
        xp = vertexc(1,ii)
        yp = vertexc(2,ii)
        zp = vertexc(3,ii)

        q(1) = xp
        q(2) = yp
        q(3) = zp

        rvec(1) = unvect(1,ii)
        rvec(2) = unvect(2,ii)
        rvec(3) = unvect(3,ii)

        ksign = 0
        if(rvec(3)>0) ksign=1

        if(icyl==0) then
          cond = yp>=ymin.AND.yp<=ymax .AND. zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
        else
          cond = zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
        endif

        if(cond) then
c if the center of the triangle is inside the local subdomain the function
c BODY2FLUID2 is used
c this function finds an intersection with the grid along RVEC and
c interpolates the pressure
          mrkpb(ii) = body2fluid2(q,rvec,xc_car,yc_car,xc,yc,zc
     &         ,p,flagp,nx,ny,nz,xint,yint,zint,pint,n1,nbd,icyl,ii,amin)

          if(mrkpb(ii)<=n1) then
c a valid intersection has been found by BODY2FLUID2
            mrkpb(ii)=1
            s = sqrt( (q(1)-xint)**2. + (q(2)-yint)**2. + (q(3)-zint)**2.)
c in the case of stationary boundaries the normal pressure gradient is
c equal to 0
            if(imov==0) then
              pbd(ii) = pint
            else
              pbd(ii) = pint-s*dpdn(xp,yp,zp,unvect(1,ii),unvect(2,ii),unvect(3,ii),t,dt,ibd)
            endif
          endif

        endif
      enddo
      
      if(imbdcmp==0 .AND. mysize>1 .AND. .false.) then
         
        IF(IVRTX==0) THEN
          n = mb(ibd)
        ELSE
          n = mv(ibd)
        ENDIF

        ALLOCATE(mrkpb1(n),mrkpbg1(n))
        ALLOCATE(pbd1(n))

        mrkpb1(1:n)=mrkpb(ilb:ile)
        pbd1(1:n)=pbd(ilb:ile)

        mrkpbg1=n1+1

        CALL MPI_ALLREDUCE(mrkpb1,mrkpbg1,n,MPI_INTEGER,MPI_MIN,MPI_COMM_EDDY,IERR)
c the minimum value of MRKPB among all processors is stored in MRKPBG1

        do ifacet=ilb,ile
          ii = ifacet-ilb+1

          if(mrkpb1(ii)<=mrkpbg1(ii) .AND. mrkpb1(ii)<=n1) then 
c valid intersections and pressure values
            mrkpb1(ii) = 1
            nvld = nvld+1
          else
c discarded intersections and pressure values
            mrkpb1(ii) = 0
            pbd1(ii) = 0.0
          endif

        enddo

c the final values of PBD,MRKPB,NVLDG are established by addiction among
c all processors
        call mpi_reduce(pbd1,pbd(ilb),n,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        call mpi_reduce(mrkpb1,mrkpb(ilb),n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        call mpi_reduce(nvld,nvldg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)

        nnvldg = ile-ilb+1-nvldg

        DEALLOCATE(pbd1,mrkpb1,mrkpbg1)

      endif


      if(myrank==0 .AND. iolvl>0) then
        if(nnvldg>0) then
          open(unit=16,file='stats_imb.dat',position='append'
     &        ,form='formatted')
          write(16,*) 'Number of valid pressure triangles',nvldg
          write(16,*) 'Number of non-valid pressure triangles',nnvldg
          close(16)
        endif
      endif

      return

      end
C-----------------------------------------------------------------------




C---- subroutine press_bdvtx ------------------N. Beratlis-26 Dec. 2008-
C
C     PURPOSE: Calculate pressure on surface on immersed body by inter-
C     polating the value of pressure on the center of triangles using
C     the nodes of the cell that contains the vertex center
C     WARNING: Works assuming flagpbdo and flagpbdi have been corrected.
C
C-----------------------------------------------------------------------
      subroutine press_bdvtx(p,flagp,xc,yc,zc,xc_car,yc_car,nx,ny,nz,nbd
     &     ,pvtx,mrkpvtx,vertex,unvect,nfacet,ibd,mbd,t,dt)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c      
c.... Input/Output arrays
      integer nx,ny,nz,nfacet,ibd,nbd,mbd
      real    t,dt
      integer mrkpvtx(3,nfacet)
      integer flagp(nx,ny,nz,nbd)
      real    xc(nx),yc(ny),zc(nz)
      real    xc_car(nx,ny),yc_car(nx,ny)
      real    p(nx,ny,nz)
      real    vertex(3,3,nfacet),unvect(3,nfacet),pvtx(3,nfacet)
c
c.... Local arrays
      integer ilb,ile,imov
      integer i,j,k,ii,jj,nb,ib,nzg,n1,ifacet,nvld,nnvld,nvldg,nnvldg,n1max,ksign
      real    xp,yp,zp,xint,yint,zint,pint,s,amin,psi,ds
      real    q(3),rvec(3)
      logical cond

      integer, dimension(:,:), allocatable :: mrkpvtx1,mrkpvtxg1
      real, dimension(:,:), allocatable :: pvtx1
c
c.... Functions
      integer vel_body2fluid,body2grid2,body2fluid2
      real    interp3d,dpdn,anglerad,mindxdydz

      
      n1max = -1
      n1 = 8
      amin = 0.0
      nnvld = 0
      nnvldg = 0
      nvld = 0


      mrkpvtx=0

      ilb = lb(ibd)+1
      ile = lb(ibd)+mb(ibd)

      imov=0
      if(ibd>=mbd) imov=1

      do ii=ilb,ile
      do jj=1,3
  
        xp = vertex(1,jj,ii)
        yp = vertex(2,jj,ii)
        zp = vertex(3,jj,ii)

        q(1) = xp
        q(2) = yp
        q(3) = zp

        rvec(1) = unvect(1,ii)
        rvec(2) = unvect(2,ii)
        rvec(3) = unvect(3,ii)

        ksign = 0
        if(rvec(3)>0) ksign=1

        if(icyl==0) then
          cond = yp>=ymin.AND.yp<=ymax .AND. zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
        else
          cond = zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
        endif

        if(cond) then
          mrkpvtx(jj,ii) = body2fluid2(q,rvec,xc_car,yc_car,xc,yc,zc
     &         ,p,flagp,nx,ny,nz,xint,yint,zint,pint,n1,nbd,icyl,ii,amin)           

          if(mrkpvtx(jj,ii)<=n1) then
            mrkpvtx(jj,ii)=1
            s = sqrt( (q(1)-xint)**2. + (q(2)-yint)**2. + (q(3)-zint)**2.)
            if(imov==0) then
              pvtx(jj,ii) = pint
            else
              pvtx(jj,ii) = pint+s*dpdn(xp,yp,zp,unvect(1,ii),unvect(2,ii),unvect(3,ii),t,dt,ibd)
            endif

          endif

        endif
      enddo
      enddo
      
      if(imbdcmp==0 .AND. mysize>1 .AND. .false.) then
         
        ALLOCATE(mrkpvtx1(3,mb(ibd)),mrkpvtxg1(3,mb(ibd)))
        ALLOCATE(pvtx1(3,mb(ibd)))

        mrkpvtx1(:,1:mb(ibd))=mrkpvtx(:,ilb:ile)
        pvtx1(:,1:mb(ibd))=pvtx(:,ilb:ile)

        mrkpvtxg1=n1+1

        CALL MPI_ALLREDUCE(mrkpvtx1,mrkpvtxg1,3*mb(ibd),MPI_INTEGER,MPI_MIN,MPI_COMM_EDDY,IERR)

        do ifacet=ilb,ile        
          ii = ifacet-ilb+1
          do jj=1,3

          if(mrkpvtx1(jj,ii)<=mrkpvtxg1(jj,ii) .AND. mrkpvtx1(jj,ii)<=n1) then           
            mrkpvtx1(jj,ii) = 1
            nvld = nvld+1
          else
            mrkpvtx1(jj,ii) = 0
            pvtx1(jj,ii) = 0.0
          endif

          enddo
        enddo

        call mpi_reduce(pvtx1,pvtx(1,ilb),3*mb(ibd),MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        call mpi_reduce(mrkpvtx1,mrkpvtx(1,ilb),3*mb(ibd),MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        call mpi_reduce(nvld,nvldg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)

        nnvldg = ile-ilb+1-nvldg

        DEALLOCATE(pvtx1,mrkpvtx1,mrkpvtxg1)

      endif


      if(myrank==0 .AND. iolvl>0) then
        if(nnvldg>0) then
          open(unit=16,file='stats_imb.dat',position='append'
     &        ,form='formatted')
          write(16,*) 'Number of valid pressure triangles',nvldg
          write(16,*) 'Number of non-valid pressure triangles',nnvldg
          close(16)
        endif
      endif

      return

      end
C-----------------------------------------------------------------------




C---- subroutine rord_forcpts ----------------N. Beratlis-16 Aug. 2008--
C
      subroutine rord_forcpts(iim,jim,kim,xnim,ynim,znim,nxim,nyim,nzim
     &     ,mrk,dir,ord,nq,lim,mim)
C
C     PURPOSE: Shift forcing points from lim-mim down by 1, and set forcing
C     point lim to end of array
C
C---------------------------------------------------------------------------
C
      IMPLICIT NONE
      include 'immersed.h'

      INTEGER lim,mim,nq,nn
      INTEGER iim(nfcmax),jim(nfcmax),kim(nfcmax),mrk(nfcmax)
     &     ,dir(nfcmax)
      INTEGER ord(nq)
      REAL    nxim(nfcmax),nyim(nfcmax),nzim(nfcmax),
     &        xnim(nfcmax),ynim(nfcmax),znim(nfcmax)

      INTEGER i,j,k,idir,imrk,iord
      REAL    xn,yn,zn
      REAL    nx,ny,nz

      i = iim(lim)
      j = jim(lim)
      k = kim(lim)
      xn = xnim(lim)
      yn = ynim(lim)
      zn = znim(lim)
      nx = nxim(lim)
      ny = nyim(lim)
      nz = nzim(lim)
      imrk = mrk(lim)
      idir = dir(lim)
      iord = ord(lim)

      iim(lim:mim-1) = iim(lim+1:mim)
      jim(lim:mim-1) = jim(lim+1:mim)
      kim(lim:mim-1) = kim(lim+1:mim)
      xnim(lim:mim-1) = xnim(lim+1:mim)
      ynim(lim:mim-1) = ynim(lim+1:mim)
      znim(lim:mim-1) = znim(lim+1:mim)
      nxim(lim:mim-1) = nxim(lim+1:mim)
      nyim(lim:mim-1) = nyim(lim+1:mim)
      nzim(lim:mim-1) = nzim(lim+1:mim)
      mrk(lim:mim-1) = mrk(lim+1:mim)
      dir(lim:mim-1) = dir(lim+1:mim)
      ord(lim:mim-1) = ord(lim+1:mim)

      iim(mim) = i
      jim(mim) = j
      kim(mim) = k
      xnim(mim) = xn
      ynim(mim) = yn
      znim(mim) = zn
      nxim(mim) = nx
      nyim(mim) = ny
      nzim(mim) = nz
      mrk(mim) = imrk
      dir(mim) = idir
      ord(mim) = iord
c      
      RETURN
      END
C---------------------------------------------------------------------------



C---- subroutine rord_forcpts ----------------N. Beratlis-16 Aug. 2008--
C
      subroutine rord_forcpts1(iim,jim,kim,xnim,ynim,znim,nxim,nyim,nzim
     &     ,mrk,dir,fp,ord,nq,lim,mim)
C
C     PURPOSE: Shift forcing points from lim-mim down by 1, and set forcing
C     point lim to end of array
C
C---------------------------------------------------------------------------
C
      IMPLICIT NONE
      include 'immersed.h'
c
c.... Input/Output Arrays
      INTEGER lim,mim,nq
      INTEGER ord(nq)
      INTEGER iim(nfcmax),jim(nfcmax),kim(nfcmax),mrk(nfcmax)
     &     ,dir(nfcmax),fp(nfcmax)
      REAL    nxim(nfcmax),nyim(nfcmax),nzim(nfcmax),
     &        xnim(nfcmax),ynim(nfcmax),znim(nfcmax)
c
c.... Local Arrays
      INTEGER i,ii,i1,i2
      INTEGER ord2(nq),fp2(nq)
      INTEGER iim2(nq),jim2(nq),kim2(nq),mrk2(nq)
     &     ,dir2(nq)
      REAL    nxim2(nq),nyim2(nq),nzim2(nq),
     &        xnim2(nq),ynim2(nq),znim2(nq)

      i1 = lim+1
      i2 = lim+nq
      iim2(1:nq) = iim(i1:i2)
      jim2(1:nq) = jim(i1:i2)
      kim2(1:nq) = kim(i1:i2)
      mrk2(1:nq) = mrk(i1:i2)
      dir2(1:nq) = dir(i1:i2)
      fp2 (1:nq) = fp (i1:i2)
      nxim2(1:nq) = nxim(i1:i2)
      nyim2(1:nq) = nyim(i1:i2)
      nzim2(1:nq) = nzim(i1:i2)
      xnim2(1:nq) = xnim(i1:i2)
      ynim2(1:nq) = ynim(i1:i2)
      znim2(1:nq) = znim(i1:i2)

      do ii=lim+1,lim+mim
        i=ord(ii-lim)-lim

        iim(ii) = iim2(i)
        jim(ii) = jim2(i)
        kim(ii) = kim2(i)
        mrk(ii) = mrk2(i)
        dir(ii) = dir2(i)
        fp (ii) = fp2 (i)
        nxim(ii) = nxim2(i)
        nyim(ii) = nyim2(i)
        nzim(ii) = nzim2(i)
        xnim(ii) = xnim2(i)
        ynim(ii) = ynim2(i)
        znim(ii) = znim2(i)

      enddo
c      
      RETURN
      END
C---------------------------------------------------------------------------




C---- subroutine rord_forcpts_p ------------------N. Beratlis-16 Aug. 2008--
C
      subroutine rord_forcpts_p(iim,jim,kim,xnim,ynim,znim,nxim,nyim,nzim
     &     ,mrk,dir,itr,ord,nq,lim,mim)
C
C     PURPOSE: Shift forcing points from lim-mim down by 1, and set forcing
C     point lim to end of array
C
C---------------------------------------------------------------------------
C
      IMPLICIT NONE
      include 'immersed.h'
c
c.... Input/Output arrays
      INTEGER lim,mim,nq,nn
      INTEGER iim(nfcmax),jim(nfcmax),kim(nfcmax),mrk(nfcmax),itr(nfcmax)
      INTEGER dir(nfcmax,2)
      INTEGER ord(nq)
      REAL    nxim(nfcmax,2),nyim(nfcmax,2),nzim(nfcmax,2),
     &        xnim(nfcmax),ynim(nfcmax),znim(nfcmax)
c
c.... Local arrays
      INTEGER i,ii

      INTEGER iim2(nq),jim2(nq),kim2(nq),mrk2(nq),itr2(nq)
      INTEGER dir2(nq,2)
      INTEGER ord2(nq)
      REAL    nxim2(nq,2),nyim2(nq,2),nzim2(nq,2),
     &        xnim2(nq),ynim2(nq),znim2(nq)

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

      do ii=lim+1,lim+mim
        i=ord(ii)

        iim(ii) = iim2(i)
        jim(ii) = jim2(i)
        kim(ii) = kim2(i)

        mrk(ii) = mrk2(i)
        nxim(ii,:) = nxim2(i,:)
        nyim(ii,:) = nyim2(i,:)
        nzim(ii,:) = nzim2(i,:)
        xnim(ii) = xnim2(i)
        ynim(ii) = ynim2(i)
        znim(ii) = znim2(i)
        dir(ii,:) = dir2(i,:)
        itr(ii) = itr2(i)

c        if((dir(ii,1)<1 .OR. dir(ii,2)<1 .OR. dir(ii,1)>3 .OR. dir(ii,2)>3) .AND. mrk(ii)/=2) then
c          write(6,*) 'WARNING rdfn_forcpts_p:',ii,i,dir(ii,:),mrk(ii),iim(ii),jim(ii),kim(ii)
c          stop
c        endif

c        write(6,*) ii,i,iim2(i),jim2(i),kim2(i)
      enddo
c      
      RETURN
      END
C---------------------------------------------------------------------------




C---- subroutine geomp1 -------------------N. Beratlis-23 Jan. 2009---
C
C     PURPOSE: Find interpolation points and intersection with immersed 
C     object for all forcing points.
C
C---------------------------------------------------------------------
      subroutine geomp1(xu,yu,zu,vertex,vertexc,unvect,lim,mim
     &    ,iim,jim,kim,itr,nim,xnim,ynim,znim,nxim,nyim,nzim,mrk,dir
     &    ,flag,flagp,flagp2,nx,ny,nz,nbd,nfacet,ibd,icom,ino)
C
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
      integer nx,ny,nz,nbd,nfacet,ibd,nim,icom,ino
      integer lim,mim
      integer iim(nfcmax),jim(nfcmax),kim(nfcmax),itr(nfcmax)
      integer mrk(nfcmax),dir(nfcmax,2)
      real    xu(nx,ny),yu(nx,ny),zu(nz)
      real    xnim(nfcmax),ynim(nfcmax),znim(nfcmax)
      real    nxim(nfcmax,2),nyim(nfcmax,2),nzim(nfcmax,2)
      integer flagp(nx,ny,nz,nbd),flagp2(nx,ny,nz,nbd),flag(nx,ny,nz)
      real    vertex(3,3,nfacet),vertexc(3,nfacet),unvect(3,nfacet)
c
c....local variables
      integer ic,nn,nimp,icntl,nt,ntg,ntmaxg,ntming,nfps,mimg,mimmax,mimmin
      integer im,jm,km,i,j,k,ilb,ile,im1,ii,in,nzg,intrs,gintrs,dintrs
      integer iext,jext,kext,iext1,jext1,kext1
      integer icount1,icount2,icount3,icount4,isten
      integer flag2(nx,ny,2)
      integer icount
      INTEGER ICOUNTSTEN,ICOUNTINTR,ICOUNTSTENG,ICOUNTINTRG
      INTEGER ICOUNTSTEN2,ICOUNTSTEN2G
      INTEGER ICOUNTSTEN_GRD,ICOUNTSTEN_GRDG,icountsten_dia,icountintr_dia
      INTEGER ICOUNTINTR_GRD,ICOUNTINTR_GRDG,icountsten_diag,icountintr_diag
      INTEGER ICOUNTGRID,ICOUNTGRIDG,ICOUNTNORM,ICOUNTNORMG
      INTEGER ICOUNTDIAG,ICOUNTDIAGG
      INTEGER ICOUNTNVLD,ICOUNTNVLDG

      integer nclocks
      REAL    clocktemp,clock(22),clockg(22),clockgmin(22),clockgmax(22)
      REAL, DIMENSION(:), ALLOCATABLE :: zug
      INTEGER, DIMENSION(:), ALLOCATABLE :: INN
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: INNF
      REAL, DIMENSION(:,:), ALLOCATABLE :: QUERYF
c
c.... functions
      INTEGER bndpts,rdfn_bndpts_p,grid_intr_p1,norm_intr_p1,diag_intr_p1
      REAL    tclock
c
c      open(unit=10,file='temp.'//char(myrank+48)//'.txt')
c      write(10,*) myrank,'inside geomp, icom',icom
c      close(10)

      nclocks = 22
      clock = 0.0
      clockg = 0.0
      clockgmin = 0.0
      clockgmax = 0.0

      clock(1) = tclock()

      NN=50 !50

      mimg = 0
      mimmax = 0
      mimmin = 0

      nimp = nim
      nt = 0
      icntl = 0
      if(icyl==1) then
        if(abs(xu(1,1))/=0.0) icntl=1
      endif

      icount = 0
      icountsten = 0
      icountsten2= 0
      icountintr = 0
      icountintr_grd = 0
      icountsten_grd = 0
      icountintr_dia = 0
      icountsten_dia = 0
      icountdiag = 0
      icountgrid = 0
      icountnorm = 0
      icountnvld = 0

      clocktemp = tclock()
      nzg = (nz-2)*mysize+2
      ALLOCATE(zug(nzg))
      call calc_zg(zu,zug,nz,nzg)
      clock(2) = tclock() - clocktemp

      clocktemp = tclock()
      lim = nim
      mim = bndpts(flagp,nx,ny,nz,ibd,nbd,iim,jim,kim,nim)
      
      nim = nimp+mim
      nfps = mim

      IF(nim>nfcmax) THEN
        WRITE(6,'(A,1X,I9)') 'ERROR: Increase nfcmax to at least',nim
        CALL MPI_FINALIZE(IERR)
      ENDIF
      clock(3) = tclock() - clocktemp


      clocktemp = tclock()
      mrk(lim+1:lim+mim)=0

      ALLOCATE(QUERYF(3,mim),INNF(mim,NN))

      DO ic=1,mim
        i = iim(ic+lim)
        j = jim(ic+lim)
        k = kim(ic+lim)
        queryf(1,ic) = xu(i,j)
        queryf(2,ic) = yu(i,j)
        queryf(3,ic) = zu(k)
      ENDDO

      clock(4) = tclock() - clocktemp

      clocktemp = tclock()
      if(icom==0 .OR. icom==3) then
        CALL COMMITTEE(QUERYF, MIM, VERTEXC, NFACET, INNF, NN, 0, MYRANK)
      endif
      clock(5) = tclock() - clocktemp

      clocktemp = tclock()
      CALL COMMITTEE(QUERYF, MIM, VERTEXC, NFACET, INNF, NN, 1, MYRANK)
      clock(6) = tclock() - clocktemp

      clocktemp = tclock()
      call fillextra_cell_z(flagp(:,:,:,ibd),flag2,nx,ny,nz,myleft,myright,mpi_comm_eddy)
      clock(7) = tclock() - clocktemp
      
      clock(8) = tclock()
      IF(mim/=0) THEN

        flag(:,:,:) = flagp(:,:,:,ibd)

        do i=lim+1,lim+mim
          ii = i-lim
          im = iim(i)
          jm = jim(i)
          km = kim(i)     

          clocktemp = tclock()
          intrs=norm_intr_p1(flag,flag2,nx,ny,nz,ibd,nbd,vertex,unvect,nfacet
     &    ,xu,yu,zu,zug,nzg,innf(ii,:),nn,xnim(i),ynim(i),znim(i)
     &    ,nxim(i,:),nyim(i,:),nzim(i,:),im,jm,km,dir(i,:),itr(i),icntl,nt,clock(16),5)
          clock(9) = clock(9) + tclock() - clocktemp

          if(intrs==3) then
            mrk(i)=1
          else
            if(intrs==0) then
              icountintr=icountintr+1
            elseif(intrs==1) then
              icountsten=icountsten+1
            elseif(intrs==2) then
              icountsten2=icountsten2+1
            endif

            clocktemp = tclock()
            gintrs = grid_intr_p1(flag,flag2,nx,ny,nz,ibd,nbd,vertex,unvect,nfacet
     &           ,xu,yu,zu,innf(ii,:),nn,xnim(i),ynim(i),znim(i)
     &           ,nxim(i,1),nyim(i,1),nzim(i,1),im,jm,km,itr(i))
            clock(10) = clock(10) + tclock() - clocktemp

            if(gintrs==1) then
              mrk(i)=2
            else

              if(gintrs==0) then
                icountintr_grd = icountintr_grd+1
              elseif(gintrs==-1) then
                icountsten_grd = icountsten_grd+1
              endif

              clocktemp = tclock()
              dintrs=diag_intr_p1(flag,flag2,nx,ny,nz,ibd,nbd,vertex,unvect,nfacet
     &             ,xu,yu,zu,zug,nzg,innf(ii,:),nn,xnim(i),ynim(i),znim(i)
     &             ,nxim(i,:),nyim(i,:),nzim(i,:),im,jm,km,dir(i,:),itr(i))
              clock(11) = clock(11) + tclock() - clocktemp
c              dintrs = 0 
              if(dintrs==3) then
                mrk(i)=3

c                if(nyim(i,1)/=0.0) then
c                if(im==100 .AND. jm==4 .AND. km+myrank*(nz-2)==273) then
c                  write(6,*) 'diag:',i,im,jm,km,int(nxim(i,1)),int(nyim(i,1)),int(nzim(i,1)),flag(im,jm,km)
c     &                  ,flag(im+int(nxim(i,1)),jm+int(nyim(i,1)),km+int(nzim(i,1)))
c                endif
              else

                if(dintrs==0) then
                  icountintr_dia = icountintr_dia+1
                else
                  icountsten_dia = icountsten_dia+1
                endif

                icountnvld = icountnvld+1
                if(ino<2) then
                  flagp(im,jm,km,ibd)=0
                elseif(ino==2) then
                  flagp(im,jm,km,ibd)=ibd
                endif

                if(ino==1) then
                  flagp2(im,jm,km,ibd)=-ibd
                endif

              endif
            endif
          endif

        enddo

        icountnorm = count(mrk(lim+1:lim+mim)==1)
        icountgrid = count(mrk(lim+1:lim+mim)==2)
        icountdiag = count(mrk(lim+1:lim+mim)==3)
        icountnvld = mim-icountnorm-icountgrid-icountdiag

        clocktemp = tclock()
        if(icountnvld>0) then          
          mim = rdfn_bndpts_p(flagp,nx,ny,nz,ibd,nbd,iim,jim,kim,xnim,ynim,znim
     &         ,nxim,nyim,nzim,mrk,dir,itr,lim,mim)
          nim = nimp+mim
        endif
        clock(12) = tclock() - clocktemp

      ENDIF

      if(nfps>0) nt = int(nt/real(nfps))
      clock(8) = tclock() - clock(8)

      clocktemp = tclock()
      if(icom>=2) then
        CALL COMMITTEE(QUERYF, MIM, VERTEXC, NFACET, INNF, NN, 2, MYRANK)
      endif
      clock(13) = tclock() - clocktemp

      DEALLOCATE(QUERYF,INNF)
      DEALLOCATE(zug)


      clocktemp = tclock()
      CALL MPI_REDUCE(MIM,ICOUNT,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(MIM,MIMMAX,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(MIM,MIMMIN,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTNORM,ICOUNTNORMG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTGRID,ICOUNTGRIDG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTDIAG,ICOUNTDIAGG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTINTR,ICOUNTINTRG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTSTEN,ICOUNTSTENG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTSTEN2,ICOUNTSTEN2G,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTSTEN_GRD,ICOUNTSTEN_GRDG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTINTR_GRD,ICOUNTINTR_GRDG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTSTEN_DIA,ICOUNTSTEN_DIAG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTINTR_DIA,ICOUNTINTR_DIAG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTNVLD,ICOUNTNVLDG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(NT,NTG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(NT,NTMAXG,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(NT,NTMING,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_EDDY,IERR)
      clock(14) = tclock() - clocktemp

      clocktemp = tclock()
      IF(MYRANK.EQ.0) THEN
        OPEN(UNIT=16, FILE='stats_imb.dat', FORM='FORMATTED'
     &      ,POSITION='APPEND')
        write(16,*) 'Ave. no. of triangles searched for normal intrs.',nt
        write(16,'(2(A,1x,I5),A,F6.2,2A,I5,A,F6.2,2A,I5,A,F6.2,A)') 
     &       ' pressure forcing points:',ICOUNT
     &       ,', forced along normal:',ICOUNTNORMG
     &       ,' (',100.*REAL(ICOUNTNORMG)/REAL(ICOUNT),'/100)'
     &       ,', forced along grid points:',ICOUNTGRIDG
     &       ,' (',100.*REAL(ICOUNTGRIDG)/REAL(ICOUNT),'/100)'
     &       ,', forced along diagonals:',ICOUNTDIAGG
     &       ,' (',100.*REAL(ICOUNTDIAGG)/REAL(ICOUNT),'/100)'
        write(16,'(A,F8.2,2(A,I6))') 'No. of forc. pts., ave:',real(icount/real(mysize))
     &       ,' ,max:',mimmax,', min:',mimmin
        
        if(icountnvldg>0) then
          write(16,'(A,I5,A)') 'WARNING:',icountnvldg
     &          ,' pressure forcing points cannot be used:'
          write(16,'(A,3(I6,A))') '  along normal '
     &         ,icountintrg,' due to no intrs., '
     &         ,icountsteng,' due to no physical 1st stencil, ' 
     &         ,icountsten2g,' due to no physical 2nd stencil'
          write(16,'(A,2(I6,A))') '  along gridlines '
     &         ,icountintr_grdg,' due to no intrs., '
     &         ,icountsten_grdg,' due to no physical stencil' 
          write(16,'(A,2(I6,A))') '  along diagonal '
     &         ,icountintr_diag,' due to no intrs.,'
     &         ,icountsten_diag,' due to no physical stencil' 
        endif
        CLOSE(16)
      ENDIF
      clock(15) = tclock() - clocktemp
        
      clock(1) = tclock() - clock(1)

      clocktemp = tclock()
      CALL MPI_REDUCE(CLOCK,CLOCKG,NCLOCKS,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(CLOCK,CLOCKGMIN,NCLOCKS,MTYPE,MPI_MIN,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(CLOCK,CLOCKGMAX,NCLOCKS,MTYPE,MPI_MAX,0,MPI_COMM_EDDY,IERR)
      clock(14) = clock(14) + tclock() - clocktemp
      clockg = clockg/real(mysize)

      where(clockg==0.0) clockg=1.e-8
      where(clockgmin==0.0) clockgmin=1.e-8
      where(clockgmax==0.0) clockgmax=1.e-8
      

      IF(MYRANK==0 .AND. .true.) THEN
        OPEN(UNIT=16,FILE='clock.dat',FORM='FORMATTED',POSITION='APPEND')
        write(16,'(A)') '---- geomp: ----------------------------------'
        write(16,'(A)') '    Task/Time           Ave. (sec/%)        Max. (sec/%)   
     &     Min (sec/%)'
        write(16,'(A,3(6x,F8.4))') 'Total               :',clockg(1),clockgmax(1),clockgmin(1)
        i=2
        write(16,905) 'Calc z global       :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=3
        write(16,905) 'Bndpts points       :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=4
        write(16,905) 'Fill query array    :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=5
        write(16,905) 'Initialize ANN      :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=6
        write(16,905) 'Find nearest neigh. :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=7
        write(16,905) 'Fill extra flag cell:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=8
        write(16,905) 'Find intrs./stencil :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=9
        write(16,905) '   -Norm. intrs.    :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        write(16,'(A,3(1x,I8))') '     No. of triangles searched:',ntg,ntmaxg,ntming
        i=16
        write(16,905) '     -segtriint     :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(9),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(9),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(9),'%)'
        i=17
        write(16,905) '     -interp_point_p:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(9),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(9),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(9),'%)'
        i=18
        write(16,905) '       -ray ext fnts:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(17),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(17),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(17),'%)'
        i=19
        write(16,905) '       -ray_face_int:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(17),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(17),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(17),'%)'
        i=20
        write(16,905) '       -fluidface   :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(17),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(17),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(17),'%)'

        i=10
        write(16,905) '   -Grid intrs.     :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=11
        write(16,905) '   -Diag. intrs.    :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=12
        write(16,905) '   -Redefine bndpt  :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=13
        write(16,905) 'Delete ANN structure:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=14
        write(16,905) 'MPI_REDUCE count arr:'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
        i=15
        write(16,905) 'I/O operations      :'
     &       ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &       ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &       ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'

        write(16,'(A)') '----------------------------------------------'
        close(16)
      ENDIF

c      open(unit=10,file='temp.'//char(myrank+48)//'.txt'
c     &     ,position='append')
c      write(10,*) myrank,'end of geomp'
c      close(10)


      RETURN

 905  format(A,3(1x,F8.4,A,F5.2,A,2x))

      END
C-------------------------------------------------------------------------


C---- subroutine grid_intr_p1 --------------------N. Beratlis-08 Jan. 2009--
C
C     PURPOSE: Given a point im,jm,km find 2 physical (fluid) points
C     along grid points (gridlines plus diagonals). It returns:
C     -1, no physical first stencil
C      0, no intersection with immersed body
C
C---------------------------------------------------------------------------
C
      integer function grid_intr_p1(flag,flag2,nx,ny,nz,ibd,nbd
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

      INTEGER iflagt,intrs,in,imin,imax,i,ii,isten,n,nitrs
      REAL    xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,ext
      REAL    rvec(3)
      REAL    stencil_normal_unitdot
      INTEGER ind(6,3)
      INTEGER im1,jm1,km1,im2,jm2,km2,id,jd,kd,nd
      LOGICAL cond
      REAL    XINT(6),YINT(6),ZINT(6),DS(6)
      INTEGER NXINT(6),NYINT(6),NZINT(6),IORDER(6),ITRIANGLE(6)


      INTRS = 0
      ISTEN = 0

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

      do i=1,nd
        id = ind(i,1)
        jd = ind(i,2)
        kd = ind(i,3)

        im1 = im + id
        jm1 = jm + jd
        km1 = km + kd

        im2 = im1 + id
        jm2 = jm1 + jd
        km2 = km1 + kd

        call index_bnd(im2,jm2,km2,im2,jm2,km2,nx,ny,nz)

        if(km==nz-1 .AND. kd==1) then
          cond = flag(im1,jm1,km1)==0 .AND. flag2(im2,jm2,2)==0
        elseif(km==2 .AND. kd==-1) then
          cond = flag(im1,jm1,km1)==0 .AND. flag2(im2,jm2,1)==0 
        else
          cond = flag(im1,jm1,km1)==0 .AND. flag(im2,jm2,km2)==0
        endif

        if(cond) then
          isten = isten+1
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
     &           ,xu(im,jm)+ext*rvec(1),yu(im,jm)+ext*rvec(2),zu(km)+ext*rvec(3)
     &           ,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)
            if(IFLAGT==1) then
              INTRS = INTRS + 1
              ITRIANGLE(INTRS) = N
              XINT (INTRS) = XMT
              YINT (INTRS) = YMT
              ZINT (INTRS) = ZMT
              NXINT(INTRS) = id
              NYINT(INTRS) = jd
              NZINT(INTRS) = kd
              go to 300
            endif

          enddo
        endif

 300    continue
        
      enddo

c 300  continue
c
c.....Shoot rays in x+,z+ direction
      IF(INTRS==0) THEN
        if(isten==0) then
          grid_intr_p1 =-1
        else
          grid_intr_p1 = 0
        endif
      ELSE
        grid_intr_p1 = 1
        IF(INTRS>6) THEN
          WRITE(6,*) 'PROBLEM: MORE THAN 6 INTERSECTIONS'
          INTRS=2
        ENDIF
        DO nitrs=1,intrs
          ii = itriangle(nitrs)
          ds(nitrs)= abs(stencil_normal_unitdot(xu(im,jm),yu(im,jm),zu(km)
     &          ,xint(nitrs),yint(nitrs),zint(nitrs)
     &          ,unvect(1,ii),unvect(2,ii),unvect(3,ii)))
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
        itr = itriangle(ii)
      ENDIF


      return
      end
C---------------------------------------------------------------------------


C---- function norm_intr_p1 --------------------N. Beratlis-06 Jan. 2009--
C
C     PURPOSE: Find normal intersection of a forcing point with immersed 
C     body and 2 physical points in the fluid along the normal vector.
C     OUTPUT:
C     0 - No normal inters. with immers. body is found
C     1 - Normal intres. with immers. body but no physical stencil points     
C     3 - Normal inters. with immers. body and 2 physical stencil points
C
C-------------------------------------------------------------------------
      integer function norm_intr_p1(flag,flag2,nx,ny,nz,ibd,nbd
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
      INTEGER iext,jext,kext,iext1,jext1,kext1
      INTEGER interp_points_pres1
      INTEGER itriangle(nn),idir(nn,2),iplane(nn,2)
      REAL    rvec(3),q(3),q2(3)
      REAL    xint(nn),yint(nn),zint(nn),nxint(nn,2),nyint(nn,2),nzint(nn,2)
     &       ,ds(nn)
c      REAL    clock(2)
      REAL    clocktemp
c
c.... Functions
      REAL    tclock
c
c      write(6,*) 'inside norm_intr_p'
c      clock = 0.0
      ext = 10.0
      ibint = 0
      intrs = 0
      norm_intr_p1=0

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
  
        clocktemp = tclock()
        call segtriint(xu(im,jm)-ext*rvec(1),yu(im,jm)-ext*rvec(2),zu(km)-ext*rvec(3)
     &       ,xu(im,jm)+ext*rvec(1),yu(im,jm)+ext*rvec(2),zu(km)+ext*rvec(3)
     &       ,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)
        clock(1) = clock(1) + tclock() - clocktemp
        if(iflagt==1) then

          ibint = ibint+1
          q(1) = xu(im,jm)
          q(2) = yu(im,jm)
          q(3) = zu(km)
          clocktemp = tclock()
          iflag = interp_points_pres1(flag,flag2,xu,yu,zu,zug,rvec,nx,ny,nz,nzg
     &      ,xinp,yinp,zinp,xinp2,yinp2,zinp2,dir,im,jm,km,icntl,clock(3),3)
c          iflag=0
          clock(2) = clock(2) + tclock() - clocktemp
c          write(6,*) in,tclock()-clocktemp,clock(2)
          if(iflag==3) then
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

            go to 402

          elseif(iflag<0) then
            write(6,*) 'No face intesrection for pressure',im,jm,km,q,rvec
c            stop
          endif
        endif

      enddo

 402  continue

      if(intrs==0) then
        if(ibint==0) then
          norm_intr_p1 = 0
        else
          norm_intr_p1 = 1
        endif
      elseif(intrs>0) then
        norm_intr_p1 = 3
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

c        write(6,*) ii,itr,itriangle(1:intrs)
      endif

      return
      end
C-----------------------------------------------------------------------


C---- function diag_intr_p1 -------------------N. Beratlis-06 Jan. 2009--
C
C     PURPOSE: Find normal intersection of a forcing point with immersed 
C     body and 2 physical points in the fluid along the normal vector.
C     OUTPUT:
C     0 - No diag inters. with immers. body is found
C     1 - Diag intres. with immers. body but no physical stencil points
C     2 - Diag intres. with immers. body but no physical stencil points     
C     3 - Diag inters. with immers. body and 2 physical stencil points
C
C------------------------------------------------------------------------
      integer function diag_intr_p1(flag,flag2,nx,ny,nz,ibd,nbd
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
      REAL    xa,ya,za,xb,yb,zb,xc,yc,zc,ext,rm,theta,dz
      REAL    xmt,ymt,zmt,xinp,yinp,zinp
      INTEGER iext,jext,kext,nd
      INTEGER itriangle(nn),idir(nn,1),id(20,3)
      REAL    rvec(3),q(3),q2(3),diag(3)
      REAL    xint(nn),yint(nn),zint(nn),nxint(nn,2),nyint(nn,2),nzint(nn,2)
     &       ,ds(nn)
c
c.... Functions
      INTEGER ray_face_int,fluidface

      ext = 10.0
      ifint = 0
      iph1 = 0
      iflag = 0
      ibint = 0
      intrs = 0
      diag_intr_p1=0

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
     

      DO i=1,nd
        imd = im+id(i,1)
        jmd = jm+id(i,2)
        kmd = km+id(i,3)
        kmdg = kmd+myrank*(nz-2)

c        call index_bnd(imd,jmd,kmd,imd,jmd,kmd,nx,ny,nz)

        diag(1) = xu(imd,jmd)-xu(im,jm)
        diag(2) = yu(imd,jmd)-yu(im,jm)
        diag(3) = zu(kmd)-zu(km)

        if(flag(imd,jmd,kmd)==0) then

          iph1 = iph1+1

c          write(6,*) 'diag:',im,jm,km,imd,jmd,kmd,nn

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

            rvec=diag
  
            call segtriint(xu(im,jm)-ext*rvec(1),yu(im,jm)-ext*rvec(2),zu(km)-ext*rvec(3)
     &           ,xu(im,jm)+ext*rvec(1),yu(im,jm)+ext*rvec(2),zu(km)+ext*rvec(3)
     &           ,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)

c            write(6,*) 'diag',in,n,im,jm,km,iflagt,rvec

            if(iflagt==1) then

              ibint = ibint+1
              q(1) = xu(imd,jmd)
              q(2) = yu(imd,jmd)
              q(3) = zu(kmd)

              do dirn=1,3

                if(dirn==1) then
                  rm = sqrt(xu(imd-1,jmd)**2. + yu(imd-1,jmd)**2.)
                  if(icyl==1 .AND. imd==2) rm=0.0
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
                  ifph(dirn) = fluidface(flag2(1,1,1),nx,ny,1,imd2,jmd2,1,dirn)
                elseif(kmd==nz .AND. kext>0) then
                  ifint(dirn) = ray_face_int(q,rvec,xu,yu,zug,nx,ny,nzg
     &                  ,imd2,jmd2,kmdg+kext,dirn,xinp,yinp,zinp,icyl)
                  ifph(dirn) = fluidface(flag2(1,1,2),nx,ny,2,imd+iext,jmd+jext,1,dirn)
                else
                  ifint(dirn) = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &                  ,imd2,jmd2,kmd+kext,dirn,xinp,yinp,zinp,icyl)
                  ifph(dirn) = fluidface(flag,nx,ny,nz,imd+iext,jmd+jext,kmd+kext,dirn)
                endif

                if(ifint(dirn)==1 .AND. ifph(dirn)==1) then
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
                  goto 350
                endif
              enddo
c 350          continue

            endif
          enddo

 350      continue

        endif
      enddo

      if(intrs==0) then
        if(ibint==0) then
          diag_intr_p1 = 0
        else
          diag_intr_p1 = 1
        endif
c        write(6,*) 'diagonal:',im,jm,km,intrs,ibint,iph1
      elseif(intrs>0) then
        diag_intr_p1 = 3
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
c        write(6,*) ii,itr,itriangle(1:intrs)
      endif

      return
      end
C-----------------------------------------------------------------------


C---- function interp_points_pres1 -----------N. Beratlis-06 Jan. 2009--
C
C     PURPOSE: Find 2 interpolations points for pressure in the fluid along
C     the normal intersection with immersred
C     OUTPUT: 
C     1 - No possible inters. with grid or physical interp. points
C     2 - Intersection with body and 1 physical interp. point
C     3 - Intersection with body and 2 physical interp. points
C
C-----------------------------------------------------------------------
      integer function interp_points_pres1(flag,flag2,xu,yu,zu,zug,rvec
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
      integer ifph,ifph1,ifph2,ifph3
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

      interp_points_pres1=1
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

c      write(6,*) im,jm,km

      q(1) = xu(im,jm)
      q(2) = yu(im,jm)
      q(3) = zu(km)

      theta = dely*(real(jm)-1.5)

      clocktemp = tclock()
      rm = sqrt(xu(im-1,jm)**2. + yu(im-1,jm)**2.)
      if(icyl==1 .AND. im==2) rm=0.0
      call vec_ijkext_gridpoint(q,rvec,iext,jext,kext,icyl,1,rm,dely) !1,-1
      clock(1) = clock(1) + tclock() - clocktemp

      clocktemp = tclock()
      ifint1 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &  ,im+iext,jm+jext,km+kext,1,xinp,yinp,zinp,icyl)
      clock(2) = clock(2) + tclock() - clocktemp

      clocktemp = tclock()
      ifph1 = fluidface(flag,nx,ny,nz,im+iext,jm+jext,km+kext,1)
      clock(3) = clock(3) + tclock() - clocktemp

      if(ifint1==1 .AND. ifph1==1) then

        interp_points_pres1=interp_points_pres1+1
        idir(1) = 1

        q(1) = xinp
        q(2) = yinp
        q(3) = zinp

        !Set new indices im1,jm1,km1 for new q
        clocktemp = tclock()
        im1 = im+iext
        jm1 = jm+jext
        km1 = km+kext
        call per_index(jm1,ny,jm1)

        rm = sqrt(xu(im1-1,jm)**2. + yu(im1-1,jm)**2.)
        if(icyl==1 .AND. im1==2) rm=0.0
        call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,1,rm,dely,1) !0,-1
        clock(1) = clock(1) + tclock() - clocktemp

        im2 = im1+iext

        clocktemp = tclock()
        ifint11 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &    ,im2,jm1,km1,1,xinp2,yinp2,zinp2,icyl)
        clock(2) = clock(2) + tclock() - clocktemp
        clocktemp = tclock()
        ifph = fluidface(flag,nx,ny,nz,im2,jm1,km1,1)
        clock(3) = clock(3) + tclock() - clocktemp

        if(ifint11==1 .AND. ifph==1) then
          idir(2)=1
          interp_points_pres1=interp_points_pres1+1
        else

          clocktemp = tclock()
          theta = dely*(real(jm1)-1.5)
          call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,2,theta,dely,1) !0,-1

          jm2 = jm1+jext
          call per_index(jm2,ny,jm2)
          clock(1) = clock(1) + tclock() - clocktemp

          clocktemp = tclock()
          ifint12 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &    ,im1+iext,jm2,km1,2,xinp2,yinp2,zinp2,icyl)
          clock(2) = clock(2) + tclock() - clocktemp

          clocktemp = tclock()
          ifph = fluidface(flag,nx,ny,nz,im1+iext,jm2,km1,2)
          clock(3) = clock(3) + tclock() - clocktemp

          if(ifint12==1 .AND. ifph==1) then
            idir(2)=2
            interp_points_pres1=interp_points_pres1+1
          else

            clocktemp = tclock()
            dz = zu(km1)-q(3)
            if(rvec(3)>=0.0) then
              dz = zu(km1+1)-q(3)
            endif
            call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,3,dz,dely,1)
            km2 = km1+kext
            clock(1) = clock(1) + tclock() - clocktemp

            clocktemp = tclock()
            ifint13 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &            ,im1+iext,jm1,km2,3,xinp2,yinp2,zinp2,icyl)
            clock(2) = clock(2) + tclock() - clocktemp
            clocktemp = tclock()
            ifph = fluidface(flag,nx,ny,nz,im1+iext,jm1,km2,3)
            clock(3) = clock(3) + tclock() - clocktemp
            if(ifint13==1 .AND. ifph==1) then
              idir(2)=3
              interp_points_pres1=interp_points_pres1+1
            endif
          endif
        endif

      else

        clocktemp = tclock()
        call vec_ijkext_gridpoint(q,rvec,iext,jext,kext,icyl,2,theta,dely)
        clock(1) = clock(1) + tclock() - clocktemp

        clocktemp = tclock()
        ifint2 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &  ,im+iext,jm+jext,km+kext,2,xinp,yinp,zinp,icyl)
        clock(2) = clock(2) + tclock() - clocktemp

        clocktemp = tclock()
        ifph2 = fluidface(flag,nx,ny,nz,im+iext,jm+jext,km+kext,2)
        clock(3) = clock(3) + tclock() - clocktemp

        if(ifint2==1 .AND. ifph2==1) then

          interp_points_pres1=interp_points_pres1+1
          idir(1) = 2

          q(1) = xinp
          q(2) = yinp
          q(3) = zinp
          
          !Set new indices im1,jm1,km1 for new q
          clocktemp = tclock()
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
          clock(1) = clock(1) + tclock() - clocktemp

          clocktemp = tclock()
          ifint21 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &         ,im2,jm11,km1,1,xinp2,yinp2,zinp2,icyl)
          clock(2) = clock(2) + tclock() - clocktemp

          clocktemp = tclock()
          ifph = fluidface(flag,nx,ny,nz,im2,jm11,km1,1)
          clock(3) = clock(3) + tclock() - clocktemp
          if(ifint21==1 .AND. ifph==1) then
            interp_points_pres1=interp_points_pres1+1
            idir(2)=1
          else

            clocktemp = tclock()
            theta = dely*(real(jm11)-1.5)
            call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,2,rm,dely,2) !0,-1
            jm2 = jm1+jext
            call per_index(jm2,ny,jm2)
            clock(1) = clock(1) + tclock() - clocktemp

            clocktemp = tclock()
            ifint22 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &            ,im1,jm2,km1,2,xinp2,yinp2,zinp2,icyl)
            clock(2) = clock(2) + tclock() - clocktemp
            clocktemp = tclock()
            ifph = fluidface(flag,nx,ny,nz,im1,jm2,km1,2)
            clock(3) = clock(3) + tclock() - clocktemp
            if(ifint22==1 .AND. ifph==1) then
              interp_points_pres1=interp_points_pres1+1
              idir(2)=2
            else

c              dz = zu(km1 + int(sign(1.0,rvec(3))))-zu(km1)
              clocktemp = tclock()
              dz = zu(km1)-q(3)
              if(rvec(3)>=0.0) then
                dz = zu(km1+1)-q(3)
              endif
              call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,3,theta,dely,2) !0,-1              
              km2=km1+kext              
              jm11 = jm1+jext
              call per_index(jm11,ny,jm11)
              clock(1) = clock(1) + tclock() - clocktemp

              clocktemp = tclock()
              ifint23 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &              ,im1,jm11,km2,3,xinp2,yinp2,zinp2,icyl)
              clock(2) = clock(2) + tclock() - clocktemp
              clocktemp = tclock()
              ifph = fluidface(flag,nx,ny,nz,im1,jm11,km2,3)
              clock(3) = clock(3) + tclock() - clocktemp
              if(ifint23==1 .AND. ifph==1) then
                interp_points_pres1=interp_points_pres1+1
                idir(2)=3
              endif
            endif
          endif

        else

          clocktemp = tclock()
          dz = zu(km + int(sign(1.0,rvec(3))))-zu(km)
          call vec_ijkext_gridpoint(q,rvec,iext,jext,kext,icyl,3,dz,dely) !1,-1
          clock(1) = clock(1) + tclock() - clocktemp

          clocktemp = tclock()
          ifint3 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &          ,im+iext,jm+jext,km+kext,3,xinp,yinp,zinp,icyl)
          clock(2) = clock(2) + tclock() - clocktemp
          clocktemp = tclock()
          ifph3 = fluidface(flag,nx,ny,nz,im+iext,jm+jext,km+kext,3)
          clock(3) = clock(3) + tclock() - clocktemp

          if(ifint3==1 .AND. ifph3==1) then
            interp_points_pres1=interp_points_pres1+1        

            idir(1) = 3

            q(1) = xinp
            q(2) = yinp
            q(3) = zinp

            clocktemp = tclock()
            im1 = im+iext
            jm1 = jm+jext
            km1 = km+kext
            call per_index(jm1,ny,jm1)

c            theta = dely*(real(jm1)-1.5)
            rm = sqrt(xu(im1,jm1)**2. + yu(im1,jm1)**2.)
            if(icyl==1 .AND. im1==1) rm=0.0
            call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,1,rm,dely,3) !0,-1
            clock(1) = clock(1) + tclock() - clocktemp

c            if(im==11 .AND. jm==4 .AND. km==189) then
c              write(6,*) '     q(1)=',q(1)
c              write(6,*) '     q(2)=',q(2)
c              write(6,*) '     q(3)=',q(3)
c              write(6,*) '     r(1)=',rvec(1)
c              write(6,*) '     r(2)=',rvec(2)
c              write(6,*) '     r(3)=',rvec(3)
c              write(6,*) '1.',im,jm,km,im1,jm1,km1,iext,jext,kext,rm
c            endif

            im2 = im1+iext
            km31 = km1+kext
        
            if(km31>=1.AND.km31<nz) then

              clocktemp = tclock()
              ifint31 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &              ,im1+iext,jm1,km1+kext,1,xinp2,yinp2,zinp2,icyl)
              clock(2) = clock(2) + tclock() - clocktemp
              clocktemp = tclock()
              ifph = fluidface(flag,nx,ny,nz,im1+iext,jm1,km1+kext,1)
              clock(3) = clock(3) + tclock() - clocktemp
              
              if(ifint31==1 .AND. ifph==1) then
                idir(2)=1
                interp_points_pres1=interp_points_pres1+1
                goto 300
              else

                clocktemp = tclock()
                theta = dely*(real(jm1)-1.5)
                call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,2,theta,dely,3) !0,-1

                jm2 =jm1+jext
                call per_index(jm2,ny,jm2)
                clock(1) = clock(1) + tclock() - clocktemp

                clocktemp = tclock()
                ifint32 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &                ,im1,jm2,km1+kext,2,xinp2,yinp2,zinp2,icyl)
                clock(2) = clock(2) + tclock() - clocktemp
                clocktemp = tclock()
                ifph = fluidface(flag,nx,ny,nz,im1,jm2,km1+kext,2)
                clock(3) = clock(3) + tclock() - clocktemp
                if(ifint32==1 .AND. ifph==1) then
                  idir(2)=2
                  interp_points_pres1=interp_points_pres1+1
                  goto 300
                endif
              endif

            elseif(km31==0) then

              if(myrank>0) then
                km31g = km31+(nz-2)*myrank
                flagt(:,:,1) = flag2(:,:,1)
                flagt(:,:,2) = flag(:,:,1)

                clocktemp = tclock()
                call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,1,rm,dely,3) !0,-1
                im2 = im1+iext
                clock(1) = clock(1) + tclock() - clocktemp

                clocktemp = tclock()
                ifint31 = ray_face_int(q,rvec,xu,yu,zug,nx,ny,nzg
     &               ,im2,jm1,km31g,1,xinp2,yinp2,zinp2,icyl)
                clock(2) = clock(2) + tclock() - clocktemp

                clocktemp = tclock()
                ifph = fluidface(flagt,nx,ny,2,im2,jm1,1,1)
                clock(3) = clock(3) + tclock() - clocktemp

                if(ifint31==1 .AND. ifph==1) then
                  idir(2)=1
                  interp_points_pres1=interp_points_pres1+1
                  goto 300
                else

                  clocktemp = tclock()
                  theta = dely*(real(jm1)-1.5)
                  call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,2,theta,dely,3) !0,-1
                  clock(1) = clock(1) + tclock() - clocktemp

                  jm2=jm1+jext
                  clocktemp = tclock()
                  ifint32 = ray_face_int(q,rvec,xu,yu,zug,nx,ny,nzg
     &                  ,im1,jm2,km31g,2,xinp2,yinp2,zinp2,icyl)
                  clock(2) = clock(2) + tclock() - clocktemp

                  clocktemp = tclock()
                  ifph = fluidface(flagt,nx,ny,2,im1,jm2,1,2)
                  clock(3) = clock(3) + tclock() - clocktemp
                  if(ifint32==1 .AND. ifph==1) then
                    idir(2)=2
                    interp_points_pres1=interp_points_pres1+1
                    goto 300
                  endif
                endif
              endif

            elseif(km31==nz) then

              if(myrank<mysize-1) then
                km31g = km31+(nz-2)*myrank
                flagt(:,:,1) = flag(:,:,nz)
                flagt(:,:,2) = flag2(:,:,2)

                clocktemp = tclock()
                call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,1,rm,dely,3) !0,-1
                im2 = im1+iext
                clock(1) = clock(1) + tclock() - clocktemp

                clocktemp = tclock()
                ifint31 = ray_face_int(q,rvec,xu,yu,zug,nx,ny,nzg
     &               ,im2,jm1,km31g,1,xinp2,yinp2,zinp2,icyl)
                clock(2) = clock(2) + tclock() - clocktemp

                clocktemp = tclock()
                ifph = fluidface(flagt,nx,ny,2,im2,jm1,1,1)
                clock(3) = clock(3) + tclock() - clocktemp
                
                if(ifint31==1 .AND. ifph==1) then
                  idir(2)=1
                  interp_points_pres1=interp_points_pres1+1
                  goto 300
                else

                  clocktemp = tclock()
                  theta = dely*(real(jm1)-1.5)
                  call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,2,theta,dely,3) !0,-1
                  jm2 = jm1+jext
                  clock(1) = clock(1) + tclock() - clocktemp

                  clocktemp = tclock()
                  ifint32 = ray_face_int(q,rvec,xu,yu,zug,nx,ny,nzg
     &                  ,im1,jm2,km31g,2,xinp2,yinp2,zinp2,icyl)
                  clock(2) = clock(2) + tclock() - clocktemp

                  clocktemp = tclock()
                  ifph =  fluidface(flagt,nx,ny,2,im1,jm2,1,2)
                  clock(3) = clock(3) + tclock() - clocktemp

                  if(ifint32==1 .AND. ifph==1) then
                    idir(2)=2
                    interp_points_pres1=interp_points_pres1+1
                    goto 300
                  endif
                endif
              endif
 
            endif

            clocktemp = tclock()
            km1g = km1+myrank*(nz-2)
            km2g = km1g+kext
            dz = zug(km2g)-zug(km1g)
            call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,3,dz,dely,3) !0,-1
            km2 = km1+kext
            clock(1) = clock(1) + tclock() - clocktemp

c            if(im==11 .AND. jm==4 .AND. km==189) then
c              write(6,*) '3.',im,jm,km,im1,jm1,km1,iext,jext,kext,dz
c            endif

            if(km2>=1.AND.km2<=nz) then

c              dz = zu(km1 + int(sign(1.0,rvec(3))))-zu(km1)
c              call vec_ijkext_face(q,rvec,iext,jext,kext,icyl,3,dz,dely,3) !0,-1
c              km2 = km1+kext

              clocktemp = tclock()
              ifint33 = ray_face_int(q,rvec,xu,yu,zu,nx,ny,nz
     &              ,im1,jm1,km2,3,xinp2,yinp2,zinp2,icyl)
              clock(2) = clock(2) + tclock() - clocktemp
              clocktemp = tclock()
              ifph = fluidface(flag,nx,ny,nz,im1,jm1,km2,3)
              clock(3) = clock(3) + tclock() - clocktemp
              if(ifint33==1 .AND. ifph==1) then
                idir(2)=3
                interp_points_pres1=interp_points_pres1+1
              endif
            elseif(km2==0.AND.myrank>0) then
              km2g=km2+myrank*(nz-2)
              clocktemp = tclock()
              ifint33 = ray_face_int(q,rvec,xu,yu,zug(km2g),nx,ny,1
     &             ,im1,jm1,1,3,xinp2,yinp2,zinp2,icyl)
              clock(2) = clock(2) + tclock() - clocktemp
              clocktemp = tclock()
              ifph = fluidface(flag2(1,1,1),nx,ny,1,im1,jm1,1,3)
              clock(3) = clock(3) + tclock() - clocktemp
              if(ifint33==1 .AND. ifph==1) then
                idir(2)=3
                interp_points_pres1=interp_points_pres1+1
              endif
            elseif(km2==nz+1.AND.myrank<mysize-1) then
              km2g=km2+myrank*(nz-2)
              clocktemp = tclock()
              ifint33 = ray_face_int(q,rvec,xu,yu,zug(km2g),nx,ny,1
     &             ,im1,jm1,1,3,xinp2,yinp2,zinp2,icyl)
              clock(2) = clock(2) + tclock() - clocktemp

              clocktemp = tclock()
              ifph = fluidface(flag2(1,1,2),nx,ny,1,im1,jm1,1,3)
              clock(3) = clock(3) + tclock() - clocktemp
              if(ifint33==1 .AND. ifph==1) then
                idir(2)=3
                interp_points_pres1=interp_points_pres1+1
              endif
            endif
          endif
        endif

      endif

 300  continue

c      if(im==23.AND.jm==40.AND.km==180) then
c        write(6,*) 'inside inter_points_pres',im,jm,km,ifint1,ifint2,ifint3,ifph1,ifph2,ifph3
c     &        ,ifint11,ifint12,ifint13,ifint21,ifint22,ifint23,ifint31,ifint32,ifint33
c      endif

c      open(unit=9,file='stats_imb.dat',form='formatted'
c     &     ,position='append')
c      if(ifint1==0.AND.ifint2==0.AND.ifint3==0) then
c        write(9,*) 'WARNING: pressure no 1st face instersection',im,jm,km
c        interp_points_pres = -1
c      endif

c      if(ifint1==1 .AND. ifph1==1) then
c        if(ifint11==0 .AND. ifint12==0 .AND. ifint13==0) then
c          write(9,*) 'WARNING: pressure no 2nd face instersection, 1:',im,jm,km
c          interp_points_pres = -2
c        endif
c      endif

c      if(ifint2==1 .AND. ifph2==1) then
c        if(ifint21==0 .AND. ifint22==0 .AND. ifint23==0) then
c          write(9,*) 'WARNING: pressure no 2nd face instersection, 2:',im,jm,km
c          interp_points_pres = -2
c        endif
c      endif

c      if(ifint3==1 .AND. ifph3==1) then
c        if(ifint31==0 .AND. ifint32==0 .AND. ifint33==0) then
c          write(9,*) 'WARNING: pressure no 2nd face instersection, 3:',im,jm,km
c          interp_points_pres = -2
c        endif
c      endif
c      close(9)

      return      
      end
C-----------------------------------------------------------------------

C---- subroutine physical_mrk2flag----------N . Beratlis-22 Aug. 2009---
C
C     PURPOSE: Set flag to physical value (0) based on value of array mrk
C
C-----------------------------------------------------------------------      
      subroutine physical_mrk2flag(flag,nx,ny,nz,mrk,iim,jim,kim,lim,mim)

      implicit none
      include 'immersed.h'
      
      integer nx,ny,nz,lim,mim
      integer mrk(nfcmax),iim(nfcmax),jim(nfcmax),kim(nfcmax)
      integer flag(nx,ny,nz)

      integer im,i,j,k

      do i=lim+1,lim+mim
        if(mrk(i)>0) then
          flag(iim(i),jim(i),kim(i))=0
        endif
      enddo
      
      return

      end
C-----------------------------------------------------------------------



C---- subroutine physical_mrk2flag1---------N . Beratlis-22 Aug. 2009---
C
C     PURPOSE: Set flag to physical value (0) based on value of array mrk
C
C-----------------------------------------------------------------------
      subroutine physical_mrk2flag1(flag,nx,ny,nz,mrk,iim,jim,kim,ord,nq,lim,mim)

      implicit none
      include 'immersed.h'
      
      integer nx,ny,nz,lim,mim,nq
      integer ord(nq)
      integer mrk(nfcmax),iim(nfcmax),jim(nfcmax),kim(nfcmax)
      integer flag(nx,ny,nz)

      integer i,ii

      do ii=lim+1,lim+mim
        i=ord(ii)
c        if(iim(i)==102 .AND. jim(i)==4 .AND. kim(i)==264) then
c          write(6,*) 'inside mrk2flag1',lim,mim,ii,i,iim(i),jim(i),kim(i),mrk(i)
c        endif
        if(mrk(i)>0) then
          flag(iim(i),jim(i),kim(i))=0
        endif
      enddo
      
      return

      end
C-----------------------------------------------------------------------


C---- function recheckface -----------------N . Beratlis-22 Aug. 2009---
C
C     PURPOSE: Check if face is physical
C     RETURNS: 1 if all points are fluid points - SUCCESS
C             -1 if face contains any forcing points (negative flag)
C             -2 if face contains any body points (postive flag)
C
C-----------------------------------------------------------------------      
      integer function recheckface(flag,xu,yu,zu,xp,yp,zp,nx,ny,nz,dir,icyl,idomy)

      implicit none
c      include 'common.h'
c
      REAL, PARAMETER :: pi=3.141592653589793

c.... Input/Output Arrays
      integer nx,ny,nz,dir,icyl,idomy
      integer flag(nx,ny,nz)
      real    xp,yp,zp
      real    xu(nx),yu(ny),zu(nz)
c
c.... Local Arrays
      integer i,j,k
      real    xint,yint,zint
c
c.... Function
      integer fluidface
      real    anglerad

      IF(icyl==1) THEN
        xint = sqrt(xp**2. + yp**2.)
        yint = anglerad(xp,yp)
        if(idomy==1) then
          if(yint<yu(1)) then
            yint = yint+2.*pi
          elseif(yint>yu(ny)) then
            yint = yint-2.*pi
          endif
        endif
      ELSE
        xint = xp
        yint = yp
      ENDIF
      zint = zp

      if(dir==1) then
        call closest(xu,nx,xint,i)
        call locate(yu,ny,yint,j)
        call locate(zu,nz,zint,k)
c        call per_index(j,ny,j)
        recheckface = fluidface(flag,nx,ny,nz,i,j,k,dir)
      elseif(dir==2) then
        call closest(yu,ny,yint,j)
        call locate(xu,nx,xint,i)
        call locate(zu,nz,zint,k)
c        call per_index(j,ny,j)
        recheckface = fluidface(flag,nx,ny,nz,i,j,k,dir)
      else
        call closest(zu,nz,zint,k)
        call locate(xu,nx,xint,i)
        call locate(yu,ny,yint,j)
c        call per_index(j,ny,j)
        recheckface = fluidface(flag,nx,ny,nz,i,j,k,dir)
      endif


      if(recheckface==0) then
        recheckface=-1
      elseif(recheckface==-1) then
        recheckface=-2
      endif

      return

      end
C-----------------------------------------------------------------------      


C---- function recheckfaces_p --------------N . Beratlis-22 Aug. 2009---
C
C     PURPOSE: Check if faces are physical
C
C-----------------------------------------------------------------------      
      integer function recheckfaces_p(flag,flag2,xu,yu,zu,nx,ny,nz,xp,yp,zp,dir,icyl)

      implicit none
c
c.... Input/Output Arrays
      integer nx,ny,nz,icyl
      integer dir(2)
      integer flag(nx,ny,nz),flag2(nx,ny,2)
      real    xp(2),yp(2),zp(2)
      real    xu(nx),yu(ny),zu(nz)
c
c.... Local Arrays
      integer i,j,k,iflag1,iflag2
      integer flagt(nx,ny,2)
      real    xint,yint,zint
c
c.... Function
      integer fluidface
      real    anglerad

      IF(icyl==1) THEN
        xint = sqrt(xp(1)**2. + yp(1)**2.)
        yint = anglerad(xp(1),yp(1))
      ELSE
        xint = xp(1)
        yint = yp(1)
      ENDIF
      zint = zp(1)

      if(dir(1)==1) then
        call closest(xu,nx,xint,i)
        call locate(yu,ny,yint,j)
        call locate(zu,nz,zint,k)
        call per_index(j,ny,j)
        iflag1 = fluidface(flag,nx,ny,nz,i,j,k,dir(1))
c        write(6,*) 'recheckface_p:',i,j,k,iflag1
c     &       ,flag(i,j,k),flag(i,j+1,k),flag(i,j,k+1),flag(i,j+1,k+1)
      elseif(dir(1)==2) then
        call closest(yu,ny,yint,j)
        call locate(xu,nx,xint,i)
        call locate(zu,nz,zint,k)
        call per_index(j,ny,j)
        iflag1 = fluidface(flag,nx,ny,nz,i,j,k,dir(1))
      else
        call closest(zu,nz,zint,k)
        call locate(xu,nx,xint,i)
        call locate(yu,ny,yint,j)
        call per_index(j,ny,j)
        iflag1 = fluidface(flag,nx,ny,nz,i,j,k,dir(1))
      endif


      !Check second point
      if(icyl==1) then
        xint = sqrt(xp(2)**2. + yp(2)**2.)
        yint = anglerad(xp(2),yp(2))
      else
        xint = xp(2)
        yint = yp(2)
      endif
      zint = zp(2)

      if(zp(2)<zu(1)) then

        flagt(:,:,1) = flag2(:,:,1)
        flagt(:,:,2) = flag(:,:,1)

        if(dir(2)==1) then
          call closest(xu,nx,xint,i)
          call locate(yu,ny,yint,j)
          call per_index(j,ny,j)
          iflag2 = fluidface(flagt,nx,ny,2,i,j,1,dir(2))
        elseif(dir(2)==2) then
          call closest(yu,ny,yint,j)
          call locate(xu,nx,xint,i)
          call per_index(j,ny,j)
          iflag2 = fluidface(flagt,nx,ny,2,i,j,1,dir(2))
        else
          call locate(xu,nx,xint,i)
          call locate(yu,ny,yint,j)
          call per_index(j,ny,j)
          iflag2 = fluidface(flagt,nx,ny,2,i,j,1,dir(2))
        endif

      elseif(zp(2)>zu(nz)) then

        flagt(:,:,1) = flag(:,:,nz)
        flagt(:,:,2) = flag2(:,:,2)

        if(dir(2)==1) then
          call closest(xu,nx,xint,i)
          call locate(yu,ny,yint,j)
          call per_index(j,ny,j)
          iflag2 = fluidface(flagt,nx,ny,2,i,j,1,dir(2))
        elseif(dir(2)==2) then
          call closest(yu,ny,yint,j)
          call locate(xu,nx,xint,i)
          call per_index(j,ny,j)
          iflag2 = fluidface(flagt,nx,ny,2,i,j,1,dir(2))
        else
          call locate(xu,nx,xint,i)
          call locate(yu,ny,yint,j)
          call per_index(j,ny,j)
          iflag2 = fluidface(flagt,nx,ny,2,i,j,2,dir(2))
        endif

      else

        if(dir(2)==1) then
          call closest(xu,nx,xint,i)
          call locate(yu,ny,yint,j)
          call locate(zu,nz,zint,k)
          call per_index(j,ny,j)
         iflag2 = fluidface(flag,nx,ny,nz,i,j,k,dir(2))
        elseif(dir(2)==2) then
          call closest(yu,ny,yint,j)
          call locate(xu,nx,xint,i)
          call locate(zu,nz,zint,k)
          call per_index(j,ny,j)
          iflag2 = fluidface(flag,nx,ny,nz,i,j,k,dir(2))
        else
          call closest(zu,nz,zint,k)
          call locate(xu,nx,xint,i)
          call locate(yu,ny,yint,j)
          call per_index(j,ny,j)
          iflag2 = fluidface(flag,nx,ny,nz,i,j,k,dir(2))
        endif

      endif

c      write(6,*) 'inside recheckfaces_p:',iflag1,iflag2
       
      if(iflag1==1 .AND. iflag2==1) then
        recheckfaces_p= 1
      else
        recheckfaces_p=-1
      endif

      return

      end
C-----------------------------------------------------------------------      




C---- function recheckgrdpts_p --------------N . Beratlis-22 Aug. 2009---
C
C     PURPOSE: Check if faces are physical
C
C-----------------------------------------------------------------------      
      integer function recheckgrdpts_p(flag,flag2,nx,ny,nz,i,j,k,xext,yext,zext)

      implicit none
c
c.... Input/Output Arrays
      integer i,j,k,nx,ny,nz,icyl
      real    xext,yext,zext
      integer flag(nx,ny,nz),flag2(nx,ny,2)

c.... Local Arrays
      integer i1,j1,k1,i2,j2,k2
      integer iflag1,iflag2
      real    xint,yint,zint
c
c.... Function
      integer fluidface
      real    anglerad

      i1 = i+int(xext)
      j1 = j+int(yext)
      k1 = k+int(zext)

      i2 = i1+int(xext)
      j2 = j1+int(yext)
      k2 = k1+int(zext)

      call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)

      !Check 1st grid point
      if(flag(i1,j1,k1)==0) then
        iflag1=1
      else
        iflag1=0
      endif

      !Check 2nd grid point
      call index_bnd(i2,j2,k2,i2,j2,k2,nx,ny,nz)

      if(k2<1) then
        if(flag2(i2,j2,1)==0) then
          iflag2=1
        else
          iflag2=0
        endif
      elseif(k2>nz) then
        if(flag2(i2,j2,2)==0) then
          iflag2=1
        else
          iflag2=0
        endif
      else
        if(flag(i2,j2,k2)==0) then
          iflag2=1
        else
          iflag2=0
        endif
      endif

c      write(6,*) iflag1,iflag2,i,j,k,i1,j1,k1,i2,j2,k2,flag(i1,j1,k1),flag(i2,j2,k2)

      if(iflag1==1 .AND. iflag2==1) then
        recheckgrdpts_p=1
      else
        recheckgrdpts_p=-2
      endif


      return

      end
C-----------------------------------------------------------------------      



C---- function recheckdiagintr_p --------------N . Beratlis-22 Aug. 2009---
C
C     PURPOSE: Check if pressure stencil along diagonal is physical
C
C-----------------------------------------------------------------------      
      integer function recheckdiagintr_p(flag,flag2,xu,yu,zu,nx,ny,nz,i,j,k,xp,yp,zp,dir,icyl)

      implicit none
c
c.... Input/Output Arrays
      integer nx,ny,nz,i,j,k,icyl
      integer dir
      integer flag(nx,ny,nz),flag2(nx,ny,2)
      real    xp(2),yp(2),zp(2)
      real    xu(nx),yu(ny),zu(nz)
c
c.... Local Arrays
      integer i1,j1,k1,iflag1,iflag2
      integer flagt(nx,ny,2)
      real    xint,yint,zint
c
c.... Function
      integer fluidface
      real    anglerad

c      if(i1>357) write(6,*) 'WARNING',i1,j1,k1,xp,yp,zp,nx,ny,nz

      i1 = i+int(xp(1))
      j1 = j+int(yp(1))
      k1 = k+int(zp(1))

      call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)


      !Check 1st grid point
      if(flag(i1,j1,k1)==0) then
        iflag1=1
      else
        iflag1=0
      endif

      !Check second point
      if(icyl==1) then
        xint = sqrt(xp(2)**2. + yp(2)**2.)
        yint = anglerad(xp(2),yp(2))
      else
        xint = xp(2)
        yint = yp(2)
      endif
      zint = zp(2)

      if(zp(2)<zu(1)) then

        flagt(:,:,1) = flag2(:,:,1)
        flagt(:,:,2) = flag(:,:,1)

        if(dir==1) then
          call closest(xu,nx,xint,i)
          call locate(yu,ny,yint,j)
          call per_index(j,ny,j)
          iflag2 = fluidface(flagt,nx,ny,2,i,j,1,dir)
        elseif(dir==2) then
          call closest(yu,ny,yint,j)
          call locate(xu,nx,xint,i)
          call per_index(j,ny,j)
          iflag2 = fluidface(flagt,nx,ny,2,i,j,1,dir)
        else
          call locate(xu,nx,xint,i)
          call locate(yu,ny,yint,j)
          call per_index(j,ny,j)
          iflag2 = fluidface(flag2,nx,ny,2,i,j,1,dir)
        endif

      elseif(zp(2)>zu(nz)) then

        flagt(:,:,1) = flag(:,:,nz)
        flagt(:,:,2) = flag2(:,:,2)

        if(dir==1) then
          call closest(xu,nx,xint,i)
          call locate(yu,ny,yint,j)
          call per_index(j,ny,j)
          iflag2 = fluidface(flagt,nx,ny,2,i,j,1,dir)
        elseif(dir==2) then
          call closest(yu,ny,yint,j)
          call locate(xu,nx,xint,i)
          call per_index(j,ny,j)
          iflag2 = fluidface(flagt,nx,ny,2,i,j,1,dir)
        else
          call locate(xu,nx,xint,i)
          call locate(yu,ny,yint,j)
          call per_index(j,ny,j)
          iflag2 = fluidface(flag2,nx,ny,2,i,j,2,dir)
        endif

      else

        if(dir==1) then
          call closest(xu,nx,xint,i)
          call locate(yu,ny,yint,j)
          call locate(zu,nz,zint,k)
          call per_index(j,ny,j)
         iflag2 = fluidface(flag,nx,ny,nz,i,j,k,dir)
        elseif(dir==2) then
          call closest(yu,ny,yint,j)
          call locate(xu,nx,xint,i)
          call locate(zu,nz,zint,k)
          call per_index(j,ny,j)
          iflag2 = fluidface(flag,nx,ny,nz,i,j,k,dir)
        else
          call closest(zu,nz,zint,k)
          call locate(xu,nx,xint,i)
          call locate(yu,ny,yint,j)
          call per_index(j,ny,j)
          iflag2 = fluidface(flag,nx,ny,nz,i,j,k,dir)
        endif

      endif

c      write(6,*) 'inside recheckfaces_p:',iflag1,iflag2
       
      if(iflag1==1 .AND. iflag2==1) then
        recheckdiagintr_p= 1
      else
        recheckdiagintr_p=-1
      endif

      return

      end
C-----------------------------------------------------------------------      


C---- subroutine flagp1(flagp,nx,ny,nz,nbd,iu,ju,ku,iv,jv,kv,iw,jw,kw
C     ,limu,mimu,limv,mimv,limw,mimw,ibd)
C
C     PURPOSE:
C
C-----------------------------------------------------------------------      
      subroutine flagp1(flagp,flagu,flagv,flagw,nx,ny,nz,nbd
     &     ,iu,ju,ku,iv,jv,kv,iw,jw,kw,limu,mimu,limv,mimv,limw,mimw,ibd)

c      implicit none
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer nx,ny,nz,nbd,ibd,limu,limv,limw,mimu,mimv,mimw
      integer flagp(nx,ny,nz,nbd)
      integer flagu(nx,ny,nz,nbd),flagv(nx,ny,nz,nbd),flagw(nx,ny,nz,nbd)
      integer iu(nfcmax),ju(nfcmax),ku(nfcmax)
      integer iv(nfcmax),jv(nfcmax),kv(nfcmax)
      integer iw(nfcmax),jw(nfcmax),kw(nfcmax)
c
c.... Local Arrays
      integer i,j,k,im,i1,j1,k1,i2,j2,k2,kg
      integer status(mpi_status_size)

      flagp = 0
      where(flagu==ibd) flagp=ibd
      where(flagv==ibd) flagp=ibd
      where(flagw==ibd) flagp=ibd

      do im=limu+1,limu+mimu
        i = iu(im)
        j = ju(im)
        k = ku(im)
        flagp(i  ,j,k,ibd) = -ibd
        flagp(i+1,j,k,ibd) = -ibd
c        call index_bnd(i,j,k,i1,j1,k1,nx,ny,nz)
c        flagp(i1,j1,k1,ibd) = -ibd
c        call index_bnd(i+1,j,k,i2,j2,k2,nx,ny,nz)
c        flagp(i2,j2,k2,ibd) = -ibd
c        if(i==18 .AND. j==2) then
c          write(6,*) '1. myrank=',myrank,i,j,k,flagw(i,j,k-1,ibd)
c        endif
      enddo

      do im=limv+1,limv+mimv
        i = iv(im)
        j = jv(im)
        k = kv(im)
        flagp(i,j  ,k,ibd) = -ibd
        flagp(i,j+1,k,ibd) = -ibd
c        call index_bnd(i,j,k,i1,j1,k1,nx,ny,nz)
c        flagp(i1,j1,k1,ibd) = -ibd
c        call index_bnd(i,j+1,k,i2,j2,k2,nx,ny,nz)
c        flagp(i2,j2,k2,ibd) = -ibd
        if(flagv(i,j-1,k,ibd)==-ibd) flagp(i,j,k,ibd) = -ibd
c        if(i==18 .AND. j==2) then
c          write(6,*) '2. myrank=',myrank,i,j,k,flagw(i,j,k-1,ibd)
c        endif
      enddo
      
      do im=limw+1,limw+mimw
        i = iw(im)
        j = jw(im)
        k = kw(im)
        flagp(i,j,k  ,ibd) = -ibd
        flagp(i,j,k+1,ibd) = -ibd
        if(flagw(i,j,k-1,ibd)==-ibd) flagp(i,j,k,ibd) = -ibd
      enddo

      if(idomy==0) then
        FLAGP(:, 1,:,IBD) = FLAGP(:,NY-1,:,IBD)
        FLAGP(:,NY,:,IBD) = FLAGP(:, 2  ,:,IBD)
c
c.....Axis
        IF(ITYPE(1)==300) THEN
          DO J=1,NY
            FLAGP(1,J,:,IBD) = FLAGP(IX1,JSYM(J),:,IBD)
          ENDDO
        ENDIF
      else        
        if(itype(1)==300) then
          call mpi_sendrecv(flagp(2,:,:,ibd),ny*nz,mpi_integer,ranksym,1
     &          ,           flagp(1,:,:,ibd),ny*nz,mpi_integer,ranksym,1
     &          ,           mpi_comm_eddy,status,ierr)
        endif
      endif

      call refreshflag(flagp(1,1,1,ibd),nx,ny,nz)

      return

      end
C-----------------------------------------------------------------------      
c
c Conversion of the tagging variables from the decomposition along Y
c to the one along Z
c
C---- subroutine imb_domy2z------------------N. Beratlis-22 Nov. 2009---
C
      subroutine imb_domy2z(iu,ju,ku,xnu,ynu,znu,nxu,nyu,nzu,mrku,diru
     &     ,fp,uim,umtrx,uindx,zug,nzg,zu,nz,ny,limu,mimu,nbd,nfl,ipr)
C
C-----------------------------------------------------------------------      
c      implicit none
      include 'common.h'
      include 'mpif.h'
      include 'immersed.h'
c
c.... Input/Outpu arrays
      integer nzg,nz,nbd,nfl,ipr,ny
      integer limu(nbd,nfl),mimu(nbd,nfl)
      integer iu(nfcmax),ju(nfcmax),ku(nfcmax),fp(nfcmax)
      integer mrku(nfcmax),diru(nfcmax)
      real    xnu(nfcmax),ynu(nfcmax),znu(nfcmax)
      real    nxu(nfcmax),nyu(nfcmax),nzu(nfcmax)
      real    uim(nfcmax)
      real    zug(nzg),zu(nz)
      real    umtrx(nsp,nsp,nfcmax)
      integer uindx(nsp,nfcmax)
c
c.... Local arrays
      integer im,i,j,k,nimu,il,ml,ip,ibd,ilp,rank,i1,i2,ii,nsend,is,ie,nim
      integer status(mpi_status_size)
c      integer mim(nbd,nfl)

      integer ind
      integer, dimension(:), allocatable :: indrcv
      integer, dimension(:), allocatable :: iu2,ju2,ku2,iu3,ju3,ku3,iu4,ju4,ku4
      integer, dimension(:), allocatable :: mrku2,diru2,mrku3,diru3,mrku4,diru4
      real, dimension(:), allocatable :: xnu2,ynu2,znu2,xnu3,ynu3,znu3,xnu4,ynu4,znu4
      real, dimension(:), allocatable :: nxu2,nyu2,nzu2,nxu3,nyu3,nzu3,nxu4,nyu4,nzu4
      real, dimension(:,:,:), allocatable :: umtrx2,umtrx3,umtrx4
      real, dimension(:), allocatable :: uim2,uim3,uim4
      integer, dimension(:,:), allocatable :: mimu2,mimu3,mimu4
      integer, dimension(:,:), allocatable :: uindx2,uindx3,uindx4
      integer, dimension(:), allocatable :: fp2,fp3,fp4

      nimu = sum(mimu)

      allocate(indrcv(mysize))
      allocate(iu2(nimu),ju2(nimu),ku2(nimu))
      allocate(mrku2(nimu),diru2(nimu))
      allocate(xnu2(nimu),ynu2(nimu),znu2(nimu))
      allocate(nxu2(nimu),nyu2(nimu),nzu2(nimu))
      allocate(umtrx2(nsp,nsp,nimu))
      allocate(uim2(nimu))
      allocate(mimu2(nbd,nfl))
      allocate(uindx2(nsp,nimu))
      if(ipr==1) allocate(fp2(nimu))

      allocate(iu3(nimu),ju3(nimu),ku3(nimu))
      allocate(mrku3(nimu),diru3(nimu))
      allocate(xnu3(nimu),ynu3(nimu),znu3(nimu))
      allocate(nxu3(nimu),nyu3(nimu),nzu3(nimu))
      allocate(umtrx3(nsp,nsp,nimu))
      allocate(uim3(nimu))
      allocate(mimu3(nbd,nfl))
      allocate(uindx3(nsp,nimu))
      if(ipr==1) allocate(fp3(nimu))

      allocate(iu4(nimu),ju4(nimu),ku4(nimu))
      allocate(mrku4(nimu),diru4(nimu))
      allocate(xnu4(nimu),ynu4(nimu),znu4(nimu))
      allocate(nxu4(nimu),nyu4(nimu),nzu4(nimu))
      allocate(umtrx4(nsp,nsp,nimu))
      allocate(uim4(nimu))
      allocate(mimu4(nbd,nfl))
      allocate(uindx4(nsp,nimu))
      if(ipr==1) allocate(fp4(nimu))

      iu4(1:nimu) = iu(1:nimu)
      ju4(1:nimu) = ju(1:nimu)
      ku4(1:nimu) = ku(1:nimu)
      mrku4(1:nimu) = mrku(1:nimu)
      diru4(1:nimu) = diru(1:nimu)
      xnu4(1:nimu) = xnu(1:nimu)
      ynu4(1:nimu) = ynu(1:nimu)
      znu4(1:nimu) = znu(1:nimu)
      nxu4(1:nimu) = nxu(1:nimu)
      nyu4(1:nimu) = nyu(1:nimu)
      nzu4(1:nimu) = nzu(1:nimu)
      umtrx4(:,:,1:nimu) = umtrx(:,:,1:nimu)
      uim4(1:nimu) = uim(1:nimu)
      mimu4 = mimu
      uindx4(:,1:nimu) = uindx(:,1:nimu)
      if(ipr==1) fp4(1:nimu) = fp(1:nimu)

      il = limu(1,1)

      ind = 0
      mimu2 = 0

c      write(6,*) myrank,'domy2z',limu,mimu

c      do ibd=1,nbd
c      do ilp=1,nfl
c        do im=limu(ibd,ilp)+1,limu(ibd,ilp)+mimu(ibd,ilp)
c          if(iu(im)==32 .AND. ku(im)==16) then
c             write(6,*) '0. domy2z',ilp,ind,im,myrank,iu(im),ju(im),ku(im),rank,iu2(ind),ju2(ind),ku2(ind),ny-2,nz-2             
c          endif
c        enddo
c      enddo
c      enddo

      do ibd=1,nbd
        do ilp=1,nfl
        
          do ip=0,mysize-1
            mimu2(ibd,ilp) = 0
            ind = 0
            do im=limu(ibd,ilp)+1,limu(ibd,ilp)+mimu(ibd,ilp)

              rank= int((ku4(im)-2)/(nz-2))

              if(rank==ip) then
                if(ip/=myrank) then
                  mimu2(ibd,ilp) = mimu2(ibd,ilp)+1
                  ind = ind+1
                  iu2(ind) = iu4(im)
                  ju2(ind) = ju4(im)+myrank*(ny-2)
                  ku2(ind) = ku4(im)-rank*(nz-2)
                  mrku2(ind) = mrku4(im)
                  diru2(ind) = diru4(im)
                  xnu2(ind) = xnu4(im)
                  ynu2(ind) = ynu4(im)
                  znu2(ind) = znu4(im)
                  nxu2(ind) = nxu4(im)
                  nyu2(ind) = nyu4(im)
                  nzu2(ind) = nzu4(im)
                  umtrx2(:,:,ind) = umtrx4(:,:,im)
                  uindx2(:,ind) = uindx4(:,im)
                  uim2(ind) = uim4(im)
                  if(ipr==1) fp2(ind) = fp4(im)

c                  if(iu2(ind)==32 .AND. ju2(ind)==4) then
c                  if(ku2(ind)>nz) then
c                    write(6,*) '1. domy2z',ilp,ind,im,myrank,rank,iu4(im),ju4(im),ku4(im),iu2(ind),ju2(ind),ku2(ind),ny-2,nz-2
c                  endif

                else
                  ind = ind+1
                  mimu2(ibd,ilp) = mimu2(ibd,ilp)+1
                  iu3(ind) = iu4(im)
                  ju3(ind) = ju4(im)+myrank*(ny-2)
                  ku3(ind) = ku4(im)-rank*(nz-2)
                  mrku3(ind) = mrku4(im)
                  diru3(ind) = diru4(im)
                  xnu3(ind) = xnu4(im)
                  ynu3(ind) = ynu4(im)
                  znu3(ind) = znu4(im)
                  nxu3(ind) = nxu4(im)
                  nyu3(ind) = nyu4(im)
                  nzu3(ind) = nzu4(im)
                  umtrx3(:,:,ind) = umtrx4(:,:,im)
                  uindx3(:,ind) = uindx4(:,im)
                  uim3(ind) = uim4(im)
                  if(ipr==1) fp3(ind) = fp4(im)

c                  if(ku3(ind)>nz) then
c                    write(6,*) '2. domy2z',ilp,ind,im,myrank,iu(im),ju(im),ku(im),rank,iu2(ind),ju2(ind),ku2(ind),ny-2,nz-2
c                  endif

c                  if(iu3(ind)==32 .AND. ju3(ind)==4) then
c                    write(6,*) '2. domy2z',ilp,ind,myrank,iu(im),ju(im),ku(im)
c                  endif

                endif
              endif
            enddo

            rank = ip

            call mpi_sendrecv(mimu2(ibd,ilp),1,mpi_integer,rank,rank,
     &           indrcv(ip+1),1,mpi_integer,rank,myrank,mpi_comm_eddy,status,ierr)
c            open(unit=10,file=char(myrank+48)//'.txt',form='formatted',position='append')
c            write(10,*) ilp,'myrank=',myrank,'sending to ',rank,'mimu=',mimu2(ibd,ilp)
c     &           ,'receiving',indrcv(ip+1)
c            close(10)


            if(mimu2(ibd,ilp)>0) then
              if(ip/=myrank) then
c                open(unit=10,file=char(myrank+48)//'.txt'
c     &                ,form='formatted',position='append')
c                write(10,*) 'mpi_send from',myrank,'to ',rank
c                close(10)
                nsend = mimu2(ibd,ilp)
                is = 1
                ie = mimu2(ibd,ilp)
                call mpi_send(iu2(is:ie),nsend,mpi_integer,rank,rank+mysize+1,mpi_comm_eddy,ierr)
                call mpi_send(ju2(is:ie),nsend,mpi_integer,rank,rank+mysize+2,mpi_comm_eddy,ierr)
                call mpi_send(ku2(is:ie),nsend,mpi_integer,rank,rank+mysize+3,mpi_comm_eddy,ierr)
                call mpi_send(mrku2(is:ie),nsend,mpi_integer,rank,rank+mysize+7,mpi_comm_eddy,ierr)
                call mpi_send(diru2(is:ie),nsend,mpi_integer,rank,rank+mysize+8,mpi_comm_eddy,ierr)
                call mpi_send(xnu2(is:ie),nsend,mtype,rank,rank+mysize+9,mpi_comm_eddy,ierr)
                call mpi_send(ynu2(is:ie),nsend,mtype,rank,rank+mysize+10,mpi_comm_eddy,ierr)
                call mpi_send(znu2(is:ie),nsend,mtype,rank,rank+mysize+11,mpi_comm_eddy,ierr)
                call mpi_send(nxu2(is:ie),nsend,mtype,rank,rank+mysize+12,mpi_comm_eddy,ierr)
                call mpi_send(nyu2(is:ie),nsend,mtype,rank,rank+mysize+13,mpi_comm_eddy,ierr)
                call mpi_send(nzu2(is:ie),nsend,mtype,rank,rank+mysize+14,mpi_comm_eddy,ierr)
                call mpi_send(umtrx2(:,:,is:ie),nsp*nsp*nsend,mtype,rank,rank+mysize+15,mpi_comm_eddy,ierr)
                call mpi_send(uindx2(:,is:ie),nsp*nsend,mpi_integer,rank,rank+mysize+16,mpi_comm_eddy,ierr)
                call mpi_send(uim2(is:ie),nsend,mtype,rank,rank+mysize+17,mpi_comm_eddy,ierr)
                if(ipr==1) call mpi_send(fp2(is:ie),nsend,mpi_integer,rank,rank+mysize+18,mpi_comm_eddy,ierr)
              endif
            endif
          enddo

          mimu(ibd,ilp) = 0

          do ip=0,mysize-1
            nim=0
            rank=ip
            if(ip/=myrank) then
              if(indrcv(ip+1)>0) then
c                open(unit=10,file=char(myrank+48)//'.txt'
c     &                ,form='formatted',position='append')
c                write(10,*) myrank,ibd,ilp,'mpi_receive from',rank,indrcv(ip+1)
c                close(10)
                nim = indrcv(ip+1)
                i1 = il+1
                i2 = il+indrcv(ip+1)
                call mpi_recv(iu(i1),indrcv(ip+1),mpi_integer,rank,myrank+mysize+1,mpi_comm_eddy,status,ierr)
                call mpi_recv(ju(i1),indrcv(ip+1),mpi_integer,rank,myrank+mysize+2,mpi_comm_eddy,status,ierr)
                call mpi_recv(ku(i1),indrcv(ip+1),mpi_integer,rank,myrank+mysize+3,mpi_comm_eddy,status,ierr)
                call mpi_recv(mrku(i1),indrcv(ip+1),mpi_integer,rank,myrank+mysize+7,mpi_comm_eddy,status,ierr)
                call mpi_recv(diru(i1),indrcv(ip+1),mpi_integer,rank,myrank+mysize+8,mpi_comm_eddy,status,ierr)
                call mpi_recv(xnu(i1),indrcv(ip+1),mtype,rank,myrank+mysize+9,mpi_comm_eddy,status,ierr)
                call mpi_recv(ynu(i1),indrcv(ip+1),mtype,rank,myrank+mysize+10,mpi_comm_eddy,status,ierr)
                call mpi_recv(znu(i1),indrcv(ip+1),mtype,rank,myrank+mysize+11,mpi_comm_eddy,status,ierr)
                call mpi_recv(nxu(i1),indrcv(ip+1),mtype,rank,myrank+mysize+12,mpi_comm_eddy,status,ierr)
                call mpi_recv(nyu(i1),indrcv(ip+1),mtype,rank,myrank+mysize+13,mpi_comm_eddy,status,ierr)
                call mpi_recv(nzu(i1),indrcv(ip+1),mtype,rank,myrank+mysize+14,mpi_comm_eddy,status,ierr)
                call mpi_recv(umtrx(:,:,i1:i2),nsp*nsp*indrcv(ip+1),mtype,rank,myrank+mysize+15,mpi_comm_eddy,status,ierr)
                call mpi_recv(uindx(:,i1:i2),nsp*indrcv(ip+1),mpi_integer,rank,myrank+mysize+16,mpi_comm_eddy,status,ierr)
                call mpi_recv(uim(i1),indrcv(ip+1),mtype,rank,myrank+mysize+17,mpi_comm_eddy,status,ierr)
                if(ipr==1) call mpi_recv(fp(i1),indrcv(ip+1),mpi_integer,rank,myrank+mysize+18,mpi_comm_eddy,status,ierr)
                il = il+nim
              endif
            else
              if(indrcv(ip+1)>0) then
c                open(unit=10,file=char(myrank+48)//'.txt'
c     &                ,form='formatted',position='append')
c                write(10,*) myrank,ibd,ilp,'mpi_receive from',rank,indrcv(ip+1)
c                close(10)
                nim = indrcv(ip+1)
                i1 = il+1
                i2 = il+nim
                is = 1
                ie = nim
                iu(i1:i2) = iu3(is:ie)
                ju(i1:i2) = ju3(is:ie)
                ku(i1:i2) = ku3(is:ie)
                mrku(i1:i2) = mrku3(is:ie)
                diru(i1:i2) = diru3(is:ie)
                xnu(i1:i2) = xnu3(is:ie)
                ynu(i1:i2) = ynu3(is:ie)
                znu(i1:i2) = znu3(is:ie)
                nxu(i1:i2) = nxu3(is:ie)
                nyu(i1:i2) = nyu3(is:ie)
                nzu(i1:i2) = nzu3(is:ie)
                umtrx(:,:,i1:i2) = umtrx3(:,:,is:ie)
                uindx(:,i1:i2) = uindx3(:,is:ie)
                uim(i1:i2) = uim3(is:ie)
                if(ipr==1) fp(i1:i2) = fp3(is:ie)
                il = il+nim
              endif
            endif
            mimu(ibd,ilp) = mimu(ibd,ilp)+nim
          enddo
c          write(6,*) 'domy2z:',ibd,ilp,myrank,mimu(ibd,ilp)
        enddo
      enddo

c      mimu = 0
      nim = 0

c      write(6,*) 'domy2z:',myrank,mimu


      do ibd=1,nbd
      do ilp=1,nfl
        if(ibd>1 .AND. ilp==1) then
          limu(ibd,ilp)=limu(ibd-1,1)+sum(mimu(ibd-1,:))
        elseif(ilp>1) then
          limu(ibd,ilp)=limu(ibd,1)+sum(mimu(ibd,1:ilp-1))
        endif
      enddo
      enddo


      deallocate(indrcv)
      deallocate(iu2,ju2,ku2,iu3,ju3,ku3,iu4,ju4,ku4)
      deallocate(mrku2,diru2,mrku3,diru3,mrku4,diru4)
      deallocate(xnu2,ynu2,znu2,xnu3,ynu3,znu3,xnu4,ynu4,znu4)
      deallocate(nxu2,nyu2,nzu2,nxu3,nyu3,nzu3,nxu4,nyu4,nzu4)
      deallocate(umtrx2,umtrx3,umtrx4)
      deallocate(uim2,uim3,uim4)
      deallocate(mimu2,mimu3,mimu4)
      deallocate(uindx2,uindx3,uindx4)
      if(ipr==1) deallocate(fp2,fp3,fp4)

      return

      end
C-----------------------------------------------------------------------      


C---- subroutine flag2int--------------------N. Beratlis-24 Nov. 2009---
C
      subroutine flag2int(flag,nx,ny,nz,ibd,iuint,juint,kuint,limu,mimu)
C
C     PURPOSE: Rearrange 3D flag array among processors.
C
C-----------------------------------------------------------------------      
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer nx,ny,nz,ibd
      integer limu,mimu
      integer iuint(nintmax),juint(nintmax),kuint(nintmax)
      integer flag(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k,int,nint,n1

      mimu=0
      int = limu     !!!!!! instead of int = 0
      
      do k=kbmin(ibd),kbmax(ibd)
      do j=jbmin(ibd),jbmax(ibd)
      do i=ibmin(ibd),ibmax(ibd)
        if(flag(i,j,k)==ibd) then
          int = int+1     !!!!!! instead of int = int+limu+1
          mimu = mimu+1
          iuint(int) = i
          juint(int) = j
          kuint(int) = k
        endif
      enddo
      enddo
      enddo

      if(int>nintmax) then
        write(6,*) 'WARNING, int>',nintmax
      endif

      return

      end
C-----------------------------------------------------------------------      


C---- subroutine intind_domy2z---------------N. Beratlis-24 Nov. 2009---
C
      subroutine intind_domy2z(iu,ju,ku,limu,mimu,nbd,ny,nz)
C
C     PURPOSE: Rearrange indices arrays for interior points from y to z
C     domain decomposition
C
C----------------------------------------------------------------------- 
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
C.... Input/Output arrays
      integer nbd,ny,nz
      integer iu(nintmax),ju(nintmax),ku(nintmax)
      integer limu(nbd),mimu(nbd)
c
c.... Local arrays
      integer rank,ip,i1,i2,il,ibd,im,nint,is,ie,nsend,nim
      integer ind,indrcv(mysize)
      integer mim(nbd)
      integer status(mpi_status_size)

      integer, dimension(:), allocatable :: iu2,ju2,ku2,mimu2
      integer, dimension(:), allocatable :: iu3,ju3,ku3,mimu3

      nint = sum(mimu)
      allocate(iu2(nint),ju2(nint),ku2(nint),mimu2(nbd))
      allocate(iu3(nint),ju3(nint),ku3(nint),mimu3(nbd))

      il = limu(1)
      ind=0
      indrcv=0
      mimu2=0

      do ibd=1,nbd
        do ip=0,mysize-1
          mimu2(ibd) = 0    !!!!!! instead of mimu2 = 0
          ind   = 0
          do im=limu(ibd)+1,limu(ibd)+mimu(ibd)

            rank= int((ku(im)-2)/(nz-2))
            if(rank==ip) then

              if(ip/=myrank) then

                ind = ind+1
                mimu2(ibd) = mimu2(ibd)+1

                iu2(ind) = iu(im)
                ju2(ind) = ju(im)+myrank*(ny-2)
                ku2(ind) = ku(im)-rank*(nz-2)

              else

                ind = ind+1
                mimu2(ibd) = mimu2(ibd)+1

                iu3(ind) = iu(im)
                ju3(ind) = ju(im)+myrank*(ny-2)
                ku3(ind) = ku(im)-rank*(nz-2)

              endif
            endif
c            if(ku(im)-rank*(nz-2)>nz) write(6,*) 'WARNING',myrank,im,iu(im),ju(im),ku(im),rank,ku(im)-rank*(nz-2)
          enddo

          rank = ip
          call mpi_sendrecv(ind,1,mpi_integer,rank,rank,
     &         indrcv(ip+1),1,mpi_integer,rank,myrank,mpi_comm_eddy,status,ierr)

          if(mimu2(ibd)>0) then
            if(ip/=myrank) then
              nsend = mimu2(ibd)
              is = 1
              ie = mimu2(ibd)
              call mpi_send(iu2(is:ie),nsend,mpi_integer,rank,rank+mysize+1,mpi_comm_eddy,ierr)
              call mpi_send(ju2(is:ie),nsend,mpi_integer,rank,rank+mysize+2,mpi_comm_eddy,ierr)
              call mpi_send(ku2(is:ie),nsend,mpi_integer,rank,rank+mysize+3,mpi_comm_eddy,ierr)
            endif
          endif
        enddo

        mimu(ibd) = 0
        do ip=0,mysize-1
          nim=0
          rank = ip
          if(ip/=myrank) then
            if(indrcv(ip+1)>0) then
              nim = indrcv(ip+1)
              i1 = il+1
              i2 = il+indrcv(ip+1)
              if(i2>nintmax) then
                write(6,*) 'WARNING, nintmax>',nintmax
c                call mpi_finalize(ierr)
                stop
              endif
              call mpi_recv(iu(i1),indrcv(ip+1),mpi_integer,rank,myrank+mysize+1,mpi_comm_eddy,status,ierr)
              call mpi_recv(ju(i1),indrcv(ip+1),mpi_integer,rank,myrank+mysize+2,mpi_comm_eddy,status,ierr)
              call mpi_recv(ku(i1),indrcv(ip+1),mpi_integer,rank,myrank+mysize+3,mpi_comm_eddy,status,ierr)
              il = il+nim
            endif
          else
            if(indrcv(ip+1)>0) then
              nim = indrcv(ip+1)
              i1 = il+1
              i2 = il+nim
              is = 1
              ie = nim
              if(i2>nintmax) then
                write(6,*) 'WARNING, nintmax>',nintmax
c                call mpi_finalize(ierr)
                stop
              endif
              iu(i1:i2) = iu3(is:ie)
              ju(i1:i2) = ju3(is:ie)
              ku(i1:i2) = ku3(is:ie)
              il = il+nim
            endif
          endif
          mimu(ibd) = mimu(ibd)+nim
        enddo
      enddo

      do ibd=1,nbd
        if(ibd>1) then
          limu(ibd)=sum(mimu(1:ibd-1))
        endif
      enddo


      deallocate(iu2,ju2,ku2,mimu2)
      deallocate(iu3,ju3,ku3,mimu3)
      
      return

      end
C-----------------------------------------------------------------------      



C---- subroutine calc_cf ------------------ N. Beratlis -27 Aug 2010 ---
C
C     PURPOSE: Calculate skin friction coefficient
C
C-----------------------------------------------------------------------      
      subroutine calc_cf(dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb
     &     ,mrks,cf,mrkcf,unvect,nfacet)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      INTEGER nfacet
      INTEGER mrks(nfacet,6),mrkcf(nfacet)
      REAL    unvect(3,nfacet)
      REAL    dudxb(nfacet),dudyb(nfacet),dudzb(nfacet)
     &       ,dvdxb(nfacet),dvdyb(nfacet),dvdzb(nfacet)
     &       ,dwdxb(nfacet),dwdyb(nfacet),dwdzb(nfacet),cf(nfacet)
c
c.... Local arrays
      integer i
      real    unorm,a1,a2,a3
      real    sxxb,syyb,szzb,sxyb,sxzb,syzb
      real    cfg(nfacet)
c
c.... Functions
      real    vecmag
c
      cf = 0.0

      do i=1,nfacet

        unorm = vecmag(unvect(:,i),3)
        a1 = unvect(1,i)/unorm
        a2 = unvect(2,i)/unorm
        a3 = unvect(3,i)/unorm

        sxxb = 2.*dudxb(i)
        syyb = 2.*dvdyb(i)
        szzb = 2.*dwdzb(i)
        sxyb = dudyb(i)+dvdxb(i)
        sxzb = dudzb(i)+dwdxb(i)
        syzb = dvdzb(i)+dwdyb(i)

c        if(sum(mrks(i,1:6))==6) then
c          mrkcf(i) = 1
          cf(i) = 2.0*ru1*(sxzb*a1 + syzb*a2 + szzb*a3)
c        endif

      enddo

      CALL MPI_REDUCE(cf,cfg,nfacet,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)


      return
c
      end
C-----------------------------------------------------------------------      


C---- function node_norm_intrs -------------------N. Beratlis-19 Dec. 2010 ---
C
C     PURPOSE: Find normal intersection with immersed body.
C     RETURN: 1 If intersection is found
C             0 If no intersection is found
C
C-----------------------------------------------------------------------------
      integer function node_norm_intrs(xp,yp,zp,vertex,unvect
     $     ,nfacet,inn,nn,xint,yint,zint,itr,ntr)
c
      IMPLICIT NONE
c
c.... Input/Output arrays
      integer nn,nfacet,itr,ntr
      real    xp,yp,zp,xint,yint,zint
      integer inn(nn)
      real    unvect(3,nfacet),vertex(3,3,nfacet)
c
c.... Local arrays
      integer in,n,iflag
      real    xa,ya,za,xb,yb,zb,xc,yc,zc,ext
      real    rvec(3)

      node_norm_intrs = 0
      ext = 10.0
      ntr = 0

      do in=1,nn

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

        call segtriint(xp,yp,zp,xp+ext*rvec(1),yp+ext*rvec(2),zp+ext
     $       *rvec(3),xa,ya,za,xb,yb,zb,xc,yc,zc,xint,yint,zint,iflag)

        if(iflag==1) then
          itr = n
          ntr = in
          node_norm_intrs = 1
          return
        endif

      enddo

      RETURN

      END
C-----------------------------------------------------------------------------

C---- subroutine geom_mod ------------------N. Beratlis-May 22 2009---
C---------------------------------------------A. Posa - 5 Nov 2012----
C
C     PURPOSE: Find intersections and interpolation stencils for
C     the forcing nodes 
C     io = 0, velocity forcing points
C     io = 1, pressure and eddy viscosity forcing points
C
C---------------------------------------------------------------------
      subroutine geom_mod(xu,yu,xu_car,yu_car,zu,vertex,vertexc,unvect,lim,mim
     &    ,iim,jim,kim,nim,xnim,ynim,znim,nxim,nyim,nzim,mrk,dir,fp
     &    ,flag,flagu,nx,ny,nz,nbd,nfacet,ibd,icom,nfu,nfumax,ivar,icycle,tlevel)     !!!!!! nfumaxg
C
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
      integer nx,ny,nz,nbd,nfacet,ibd,nfu,nfumax,nfumaxg,icycle,ivar
      integer lim(nfu),mim(nfu)
      integer iim(nfcmax),jim(nfcmax),kim(nfcmax)
      integer mrk(nfcmax),dir(nfcmax),fp(nfcmax)
      real    tlevel
      real    xu_car(nx,ny),yu_car(nx,ny)
      real    xu(nx),yu(ny),zu(nz)
      real    xnim(nfcmax),ynim(nfcmax),znim(nfcmax)
      real    nxim(nfcmax),nyim(nfcmax),nzim(nfcmax)
      integer flagu(nx,ny,nz,nbd),flag(nx,ny,nz)
      real    vertex(3,3,NFACET),vertexc(3,nfacet),unvect(3,nfacet)
      integer nim,icom,srtcrt
c
c....local variables and arrays
      integer iflagt,imin,intrs,gintrs,nitrs,n,nn,ibint,ilp,ip,intrs1,cintrs,io
      integer im,jm,km,i,j,k,ilb,ile,ic,ii,in,nimu,mimg,mimmax,mimmin
     &     ,mimu,mimug,nord,nordg,iord,mimprev,mimprevg,mimdif,ndif,itri
      integer i1,j1,k1,i2,icm,ncm
      real    xint,yint,zint,anglerad,q(3)
      integer icountnorm,icountgrid,icountnostenc,icountnointrs
      integer icountnormg,icountgridg,icountnostencg,icountnointrsg
      integer icountnrmnostenc,icountnrmnostencg,icountgrdnostenc,icountgrdnostencg,icountclonostenc,icountclonostencg
      integer icountnrmnointrs,icountnrmnointrsg,icountgrdnointrs,icountgrdnointrsg
      integer icountnovld(nfu),icountnovldg
      integer kmg,nq,inonvld
      integer nnrm,nnrmrchk,ngrd,ngrdrchk,nvld,nmrk
      integer nnrmg,nnrmrchkg,ngrdg,ngrdrchkg
      character*80 strio
      integer ilu,iru
      integer, dimension(:,:), allocatable :: indexl,indexr
 
c.... functions
      INTEGER forcpts,rdfn_forcpts,fluidface,ray_face_int,grid_intr  
     $     ,norm_intr,recheckface,grid_intr_vel,closest_intr

      REAL    CLOCKTEMP,tclock
      INTEGER NCLOCKS
      REAL, DIMENSION(:), ALLOCATABLE :: CLOCK,CLOCKG,CLOCKGMAX,CLOCKGMIN

      INTEGER, DIMENSION(:), ALLOCATABLE :: ORD,ORDC,ORDN,NTRI,INN,NNTRI
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: INNF
      REAL, DIMENSION(:,:), ALLOCATABLE :: QUERYF

      nclocks=21
      ALLOCATE(CLOCK(NCLOCKS),CLOCKG(NCLOCKS),CLOCKGMAX(NCLOCKS),CLOCKGMIN(NCLOCKS))
      clock = 0.0
      clockg = 0.0
      clockgmin = 0.0
      clockgmax = 0.0

      clock(1) = tclock()

      ncm=3
      ALLOCATE(nntri(ncm))
      nntri(1) = 3
      nntri(2) = 5
      nntri(3) = nntr3D

      ALLOCATE(ntri(nntr3d))
      ntri = 0

      nimu = nim
      nfumax = 0
      nfumaxg = 0

      icountnovld = 0
      icountnovldg = 0
      icountgrid = 0
      icountgridg = 0
      icountnorm = 0
      icountnormg = 0
      icountnostenc = 0
      icountgrdnointrs = 0
      icountgrdnointrsg = 0
      icountnrmnostenc = 0
      icountnrmnostencg = 0
      icountnrmnointrs = 0
      icountnrmnointrsg = 0
      icountgrdnostenc = 0
      icountgrdnostencg = 0
      icountclonostenc = 0

      mimg = 0
      mimmax = 0
      mimmin = 0

      nnrm = 0
      nnrmg = 0
      ngrd = 0
      ngrdg = 0
      nnrmrchk = 0
      nnrmrchkg = 0
      ngrdrchk = 0
      ngrdrchkg = 0
      nord = 0
      nvld = 0
      nmrk = 0

      io = 0
      if(ivar>=4) io=1

      clock(2) = tclock()
      lim(1) = nim
      if(ivelrcn==0 .OR. io>0) then
c case for velocity forcing points when the GEOM_WM has been previously called
c and for non-velocity forcing points
        lim = 0
        mim = 0
        lim(1) = nim
c the number of the interface points is established in MIM(1)
c using the variable FLAGU (which is -IBD for the boundary points)
c in IIM,JIM,KIM the indices of the interface points are stored
        mim(1) = forcpts(flagu,nx,ny,nz,ibd,nbd,iim,jim,kim,nim)  
        flag(:,:,:) = flagu(:,:,:,ibd)
      else
        lim(2:nfu) = 0
        mim(2:nfu) = 0
        call refreshflag(flagu(:,:,:,ibd),nx,ny,nz)    !!!!!!
        flag(:,:,:) = flagu(:,:,:,ibd)
        do i=1,lim(1)
          im = iim(i)
          jm = jim(i)
          km = kim(i)
          flag(im,jm,km) = 0
        enddo
        call refreshflag(flag(1,1,1),nx,ny,nz)
      endif
      nim = nimu + mim(1)

      inonvld = 0
      if(io==1) inonvld=ibd 

      if(nim>nfcmax) then
        write(6,'(A,1X,I9)') 'ERROR: Increase nfcmax to at least',nim
        call mpi_finalize(ierr)
        stop
      endif
      mrk(lim(1)+1:lim(1)+mim(1))=0
      clock(2) = tclock() - clock(2)

      nq = mim(1)
      allocate(innf(nq,nnTR3D+2),queryf(3,nq),ord(nq),ordc(nq),ordn(nq),inn(nntr3D))
      
      allocate(indexl(nq,2),indexr(nq,2))

      do i=1,nq
        ii = lim(1)+i
        queryf(1,i) = xu_car(iim(ii),jim(ii))
        queryf(2,i) = yu_car(iim(ii),jim(ii))
        queryf(3,i) = zu(kim(ii)) 
        ordc(i) = i+lim(1)
      enddo
     
      clock(3) = tclock()
!!!!!!      if(icom==0 .AND. nfacet>0) then
      if(icom<1 .AND. nfacet>0) then
c the subroutine COMMITTEE3 is initialized in the case ICOM is less than 1
        call committee3(QUERYF, VERTEXC, NFACET, INNF, NN, 0, MYRANK, imb3D)
        if(icom==-1) icom=2   !!!!!!
c at the end of the subroutine the variables of COMMITTEE3 are deallocated
      endif
      clock(3) = tclock() - clock(3)

      clock(5) = tclock()
      
      mimprev = 0
      ndif=0

c NFU: available number of loops for the search of valid interface points
      do ilp=1,nfu

c MIM: number of interface points for which to find a valid stencil
c MIM takes into account the points for which the stencil at the previous loop
c was composed of at least one interface point (the points for which no
c intersection was found or for which the stencil involved one body point
c were discarded at the previous loops)
        mimu = mim(ilp)
        call mpi_allreduce(mimu,mimmax,1,mpi_integer,mpi_max,mpi_comm_eddy,ierr)
        call mpi_allreduce(mimu,mimug,1,mpi_integer,mpi_sum,mpi_comm_eddy,ierr)
        call mpi_allreduce(mimprev,mimprevg,1,mpi_integer,mpi_sum,mpi_comm_eddy,ierr)
        
        mimdif = mimug-mimprevg
        mimprev = mimu
        if(mimdif==0) ndif=ndif+1

        if(ilp>1) then
          lim(ilp)=nimu + sum(mim(1:ilp-1))
        endif
        nmrk = 0

        if(mimmax==0) exit  ! no left interface points

        nfumax = nfumax+1  ! max number of loops to find valid interface points

        ii=lim(ilp)+1
        do while(ii<=lim(ilp)+mimu)

          i = ordc(ii-lim(ilp))

          im = iim(i)
          jm = jim(i)
          km = kim(i)

          if(mrk(i)==-1) then
c at the previous loop the external stencil along the normal (norm_intr) had
c at least 1 interface point
            !Check if face is intersecting fluid
c  1: all face points are fluid
c -1: at least one face point is an interface point
c -2: at least one face point is an interior point
            clocktemp = tclock()
            intrs = recheckface(flag,xu,yu,zu,nxim(i),nyim(i),nzim(i),nx,ny,nz,dir(i),icyl,idomy) 
            nnrmrchk = nnrmrchk+1
            clock(16) = clock(16) + tclock() - clocktemp
          elseif(mrk(i)==0) then

            icm = 1
            do while(icm<=ncm)
              nn = nntri(icm)
              innf(i-lim(1),nntr3d+1) = nn
              innf(i-lim(1),nntr3d+2) = icm
              q(:) = queryf(:,i-lim(1))
c committee3 finds the points of the array VERTEXC closest to Q considering the 3D distance
              call committee3(q, vertexc, nfacet, inn, nn, 1, myrank, imb3D)
              innf(i-lim(1),1:nn) = inn(1:nn)
              clocktemp = tclock()
c intrs=1 if a physical stencil is found; the output of committee3 is used
c in xnim, ynim and znim the coordinates of the intersection with
c the immersed-boundary are written
c in nxim, nyim and nzim the coordinates of the intersection with a grid face
c are provided by the function NORM_INTR
c fp is the index of the triangle where the intersection is found
c dir is the direction along which the intersection is found
c intrs=0: no intersection with the immersed boundary
c intrs=-1: the stencil contains interface points
c intrs=-2: the stencil contains body points
              intrs = norm_intr(flag,nx,ny,nz,ibd,nbd,vertex,unvect,nfacet
     &             ,xu_car,yu_car,zu,innf(i-lim(1),:),nn,xnim(i),ynim(i),znim(i),nxim(i),nyim(i)
     &             ,nzim(i),fp(i),im,jm,km,dir(i),1,itri)
              clock(12) = clock(12) + tclock() - clocktemp
              if(intrs/=0) exit   ! an intersection with the immersed-boundary has been found
              icm = icm+1
            enddo
            if(itri>0) ntri(itri) = ntri(itri)+1
            if(intrs==0) then
c no intersections found along the normal direction: fp is set equal to -1 for that interface point
              icountnrmnointrs = icountnrmnointrs+1
              fp(i)=-1
            endif
            nnrm = nnrm+1
          else
            intrs = 0
          endif
         
          if(im<=2) intrs=0
!          if(((ibd.eq.3).or.(ibd.eq.4)).and.(ivar.ne.4)) intrs=0  !!!!!! 

          if(intrs==1) then
c the stencil along the direction normal to the immersed-boundary is OK
            mrk(i)=1
            icountnorm = icountnorm+1
            nvld = nvld+1
            ord(nvld)=i
          elseif(intrs==-1 .AND. ilp<nfu .AND. mimdif/=0) then
c the stencil has at least one interface point: it is not valid
c another attempt is planned for the following loop
c the number of interface points at the current loop is decreased
c the number of interface points at the next loop is increased
            mrk(i) =-1
            nmrk = nmrk+1
            ordn(nmrk) = i
            mim(ilp) = mim(ilp)-1
            mim(ilp+1)=mim(ilp+1)+1
          else

            if(mrk(i)==-2 .AND. ndif<3) then
c if at the previous loop the exterior point found by GRID_INTR was an
c interface point (MRK(I)=-2)
              clocktemp = tclock()
              if(flag(im+int(nxim(i)),jm+int(nyim(i)),km+int(nzim(i)))==0) then
c the point along the outward grid direction is fluid: it can be used
                gintrs = 1
              else
c the point along the outward grid direction is a boundary point: it is not good
                gintrs = -1
              endif
              clock(17) = clock(17) + tclock() - clocktemp
              ngrdrchk = ngrdrchk+1
            else!if(mrk(i)/=-3) then
              clocktemp = tclock()
              srtcrt = 0
              if(fp(i)>0) then
                srtcrt=2
                if(ndif<2) srtcrt=1
              endif
                nn = innf(i-lim(1),nntr3d+1)
                icm = innf(i-lim(1),nntr3d+2)
                do while (icm<=ncm)
                  clocktemp = tclock()
c
c the intersection and the interpolation stencil are sought along the grid directions
c  1: the intersection is found and the stencil is physical;
c  0: the intersection is not found;
c -1: the intersection is found, but the projection point is boundary
c -2: the intersection is found, but the projection point is interior
c also in this case xnim, ynim and znim are the coordinates of the intersection
c point with the immersed-boundary, but nxim, nyim and nzim are the extensions
c of the indices of the exterior point (not the coordinates!)
c
                  gintrs = grid_intr(flag,nx,ny,nz,ibd,vertex,unvect,nfacet
     &                 ,xu_car,yu_car,zu,innf(i-lim(1),:),nn,xnim(i),ynim(i),znim(i)
     &                 ,nxim(i),nyim(i),nzim(i),im,jm,km,srtcrt,fp(i),io,itri)

                  clock(13) = clock(13) + tclock() - clocktemp

                  if(gintrs/=0 .AND. gintrs/=-2) exit
c
c if the intersection with the immersed-boundary has not been found or if the projection point
c is a body point: the search is performed on a larger number of triangles
c
                  icm = icm+1
                  if(icm<=ncm) then
                    nn = nntri(icm)
                    innf(i-lim(1),nntr3d+1) = nn
                    innf(i-lim(1),nntr3d+2) = icm
                    q(:) = queryf(:,i-lim(1))
                    call committee3(q, vertexc, nfacet, inn, nn, 1, myrank, imb3D)
                    innf(i-lim(1),1:nn) = inn(1:nn)
                  endif
                enddo
                if(itri>0) ntri(itri) = ntri(itri)+1
                if(mrk(i)==0 .AND. gintrs==0) icountgrdnointrs = icountgrdnointrs+1
c no intersections with the surface of the immersed-boundary along the grid directions
c
                ngrd = ngrd+1
c            else
c              gintrs = 0
            endif

            if(gintrs==1) then
c the exterior point along the grid direction is a fluid point: it is OK
              mrk(i) = 2
              nvld = nvld+1
              ord(nvld) = i
            elseif(gintrs==-1 .AND. ilp<nfu) then !stencil is forcing point
c the exterior point is a boundary point
c another attempt is planned for the following loop
c also in this case the number of interface points at the current loop
c is decreased and the one at the next loop is increased
              mrk(i) = -2
              nmrk = nmrk+1
              ordn(nmrk) = i
              mim(ilp) = mim(ilp)-1
              mim(ilp+1)=mim(ilp+1)+1
            else            !stencil is body point or no intrs. with body was found
              if(io==3) then
c the number of interface points is decreased: the current interface point is not valid
c this part seems not working in the current version
                write(6,*) 'point not used',ilp,i,mrk(i),im,jm,km+myrank*(nz-2),intrs,gintrs,io
                mrk(i) = 0
                flagu(im,jm,km,ibd) = inonvld
                icountnovld(ilp) = icountnovld(ilp)+1
                mim(ilp) = mim(ilp)-1

              else
                if(mrk(i)==-3) then
               !Check if face is intersecting fluid
c this check is associated with a previous call of the function CLOSEST_INTR
c  1: all face points are fluid
c -1: at least one face point is an interface point
c -2: at least one face point is an interior point
                   clocktemp = tclock()
                   cintrs = recheckface(flag,xu,yu,zu,nxim(i),nyim(i),nzim(i),nx,ny,nz,dir(i),icyl,idomy) 
                   nnrmrchk = nnrmrchk+1
                   clock(16) = clock(16) + tclock() - clocktemp
                else
                  clocktemp = tclock()
c
c the intersection of a forcing point with the immersed-boundary is found using the node closest to
c the surface of the body
c -2: the fluid face contains at least 1 body point
c -1: the fluid face contains at least 1 interface point
c  1: the fluid face contains only fluid points and therefore it is valid 
c
                  cintrs = closest_intr(flag,nx,ny,nz,ibd,nbd,vertex,unvect,nfacet
     &                 ,xu_car,yu_car,zu,innf(i-lim(1),:),innf(i-lim(1),nntr3d+1),xnim(i),ynim(i),znim(i),nxim(i),nyim(i)
     &                 ,nzim(i),fp(i),im,jm,km,dir(i),1)
c                  write(6,*) 'geom, closest_intr:',ilp,i,mrk(i),im,jm,km,intrs,gintrs,cintrs,io
                  clock(12) = clock(12) + tclock() - clocktemp
                  nnrm = nnrm+1
                endif

                if(cintrs==1) then
c the interface point is OK
                  mrk(i)=3
                  icountnorm = icountnorm+1
                  nvld = nvld+1
                  ord(nvld)=i
                elseif(cintrs==-1 .AND. ilp<nfu .AND. ndif<2) then
c the stencil is not OK: it involves at least 1 interface point
c another attempt is planned for the following loop
c also in this case the number of interface points at the current loop
c is decreased and the one at the next loop is increased
                  mrk(i) =-3
                  nmrk = nmrk+1
                  ordn(nmrk) = i
                  mim(ilp) = mim(ilp)-1
                  mim(ilp+1)=mim(ilp+1)+1
                else
c the stencil involves at least 1 body point: the interface point is discarded
                  mrk(i) = 0
                  flagu(im,jm,km,ibd) = inonvld
                  icountnovld(ilp) = icountnovld(ilp)+1
                  mim(ilp) = mim(ilp)-1
                  icountclonostenc = icountclonostenc+1
                  write(6,'(A,20(1x,I7))') 'geom, not valid:',ibd,ivar,ilp,i,mrk(i),im,jm,km,km+myrank*(nz-2),intrs,gintrs,cintrs
     &                 ,io,mimdif,ndif,innf(i-lim(1),nntr3d+1)
                endif
              endif
            endif
          endif

          ii = ii+1

        enddo

        clocktemp = tclock()
c
c flag is set to 0 in the points where mrk>0
c
!        call physical_mrk2flag1(flag,nx,ny,nz,mrk,iim,jim,kim,ord,nq,lim(ilp)-nimu,mim(ilp))
        call physical_mrk2flag1_mod(flag,nx,ny,nz,mrk,iim,jim,kim,ord,nq,lim(ilp)-nimu,mim(ilp),ilu,iru,indexl,indexr)
        clock(18) = clock(18) + tclock() - clocktemp

        clocktemp = tclock()
!        if(mimu>0) then
!          if(idomy==0) then
!c
!c update of flag at the ghost layers along the Y direction
!c
!            call yrefresh3darray(flag,nx,ny,nz)
!          endif
!        endif
!
!        if(mimmax>0) then
!c.....Axis
!c
!c update of flag at the ghost layers at the axis
!c
!          if(itype(1)==300) call refreshcntln(flag,nx,ny,nz)
!c
!c update of flag at the ghost layers between processors
!c
!          call refreshflag(flag(1,1,1),nx,ny,nz)
!        endif

        if(mimmax>0) then
c
c update of flag at the ghost layers between processors
c
!          call refreshflag(flag(1,1,1),nx,ny,nz)
          call local_refreshflag_smp(indexl(1:ilu,:),ilu,indexr(1:iru,:),iru,flag,nx,ny,nz)
c.....Axis
c
c update of flag at the ghost layers at the axis
c
          if(itype(1)==300) call refreshcntln(flag,nx,ny,nz)
        endif

        if(mimu>0) then
          if(idomy==0) then
c
c update of flag at the ghost layers along the Y direction
c
            call yrefresh3darray(flag,nx,ny,nz)
          endif
        endif

        clock(14) = clock(14) + tclock() - clocktemp

        ordc = ordn
   
C.....redefine boundary points
c.....the number of interface points is updated in the case icountnovld>0
c     (some forcing points were discarded)
        clocktemp = tclock()
        IF(icountnovld(ilp)>0) THEN
          nim = nimu + sum(mim)
        ENDIF
        clock(15) = clock(15) + tclock() - clocktemp
c        if(ndif==2) exit

      enddo

      nim = nimu + sum(mim)
      CALL MPI_ALLREDUCE(nfumax,nfumaxg,1,MPI_INTEGER,MPI_MAX,MPI_COMM_EDDY,IERR)

 514  continue
      
      clocktemp = tclock()
c
c the vectors of the valid interface points are ordered
c ORD is a vector associated with the points for which INTRS, GINTRS or CINTRS
c are equal to 1 (valid boundary points)
c
      call rord_forcpts1(iim,jim,kim,xnim,ynim,znim,nxim,nyim,nzim
     &     ,mrk,dir,fp,ord,nq,lim(1),sum(mim))

      clock(21) = clock(21) + tclock() - clocktemp

      clock(5) = tclock() - clock(5)

      clock(8) = tclock()
      if(icom==2 .AND. nfacet>0) then
c the variables of COMMITTEE3 are deallocated
        CALL COMMITTEE3(QUERYF, VERTEXC, NFACET, INNF, NN, 2, MYRANK, imb3D)
      endif
      clock(8) = tclock() - clock(8)

      deallocate(queryf,innf,ord,ordc,ordn,inn)
c
c number of forcing points for which a valid stencil is found along the
c normal direction
      icountnorm = count(mrk(lim(1)+1:lim(1)+sum(mim))==1)
c
c number of forcing points for which a valid stencil is found along a
c grid direction
      icountgrid = count(mrk(lim(1)+1:lim(1)+sum(mim))==2)
c
c number of forcing points for which the stencil along the normal direction
c is not valid
      icountnrmnostenc = nq-icountnorm-icountnrmnointrs
c
c number of forcing points for which the stencil along a grid direction
c is not valid
      icountgrdnostenc = nq-icountnorm-icountgrid-icountgrdnointrs

      clocktemp = tclock()
c
c reduce to 0 process the statistics of the tagging procedure
c
      CALL MPI_REDUCE(SUM(ICOUNTNOVLD),ICOUNTNOVLDG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTNORM,ICOUNTNORMG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTGRID,ICOUNTGRIDG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(MIM),MIMG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(MIM),MIMMAX,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(SUM(MIM),MIMMIN,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTGRDNOSTENC,ICOUNTGRDNOSTENCG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTNRMNOINTRS,ICOUNTNRMNOINTRSG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTNRMNOSTENC,ICOUNTNRMNOSTENCG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTGRDNOINTRS,ICOUNTGRDNOINTRSG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ICOUNTCLONOSTENC,ICOUNTCLONOSTENCG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(nnrm,nnrmg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ngrd,ngrdg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(nnrmrchk,nnrmrchkg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(ngrdrchk,ngrdrchkg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(nord,nordg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      clock(9) = tclock() - clocktemp

      clocktemp = tclock()
      strio = 'forcing points'
      if(ivar==1) then 
        strio = 'U velocity '//strio
      elseif(ivar==2) then 
        strio = 'V velocity '//strio
      elseif(ivar==3) then 
        strio = 'W velocity '//strio
      elseif(ivar==4) then 
        strio = 'Pressure '//strio
      elseif(ivar==5) then 
        strio = 'viscosity '//strio
      endif
c
c the 0 process writes the statistics on a file
c
      IF(MYRANK.EQ.0 .AND. iolvl>0 .AND. MIMG.NE.0) THEN
        OPEN(UNIT=16, FILE='stats_imb.dat', FORM='FORMATTED'
     &        ,POSITION='APPEND')
        write(16,'(A,I6,A,E16.8,1x,2(A,1x,I9),A,F6.2,A,A,I8,A,F6.2,A)')
     &       'Cycle no=',icycle,', time=',tlevel
     &       ,trim(strio),MIMG
     &       ,', forced along normal:',ICOUNTNORMG
     &       ,' (',100.*REAL(ICOUNTNORMG,8)/REAL(MIMG,8),'%)'
     &       ,', forced along grid points:',ICOUNTGRIDG  
     &       ,' (',100.*REAL(ICOUNTGRIDG,8)/REAL(MIMG,8),'%)'
        write(16,'(A,F14.2,3(A,I9))') 'No. of forc. pts., ave:',real(mimg,8)/real(mysize,8)
     &       ,' ,max:',mimmax,', min:',mimmin,', Max no. of loops:',nfumaxg
        write(16,'(2(A,I9))') 'No. nrm. intrs calls=',nnrmg
     &       ,',no grd intrs calls=',ngrdg
        write(16,'(2(A,I9))') 'No. nrm. intrs rechecks=',nnrmrchkg
     &       ,',no grd intrs rechecks=',ngrdrchkg
c        write(16,'(A)') 'Triangle intersection stats:'
c        do i=1,nn
c          write(16,'(I4,1x,I8)') i,ntri(i)
c        enddo
        IF(ICOUNTNOVLDG>0) THEN
          WRITE(16,'(2(A,1X,I9,1X),A)')
     &          'WARNING:',ICOUNTNOVLDG,' of ',MIMG+ICOUNTNOVLDG
     &          ,trim(strio)//' forcing points can not be used:'
          write(16,'(A,2(I9,A))') 'along normal:'
     &          ,icountnrmnostencg,' due to no stencil,'
     &          ,icountnrmnointrsg,' due to no intrs.'
          write(16,'(A,2(I9,A))') 'along gridlines:'
     &          ,icountgrdnostencg,' due to no stencil,'
     &          ,icountgrdnointrsg,' due to no intrs.'
          write(16,'(A,1(I9,A))') 'along closest:'
     &          ,icountclonostencg,' due to no stencil.'
        ENDIF
        CLOSE(16)
      ENDIF
      clock(10) = tclock()-clocktemp

      clock(1) = tclock() - clock(1)

      IF(ioclock>0) THEN
c
c the 0 process writes on a file the computational times for the tagging procedure
c
        CALL MPI_REDUCE(CLOCK,CLOCKG,NCLOCKS,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(CLOCK,CLOCKGMIN,NCLOCKS,MTYPE,MPI_MIN,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(CLOCK,CLOCKGMAX,NCLOCKS,MTYPE,MPI_MAX,0,MPI_COMM_EDDY,IERR)
        clockg = clockg/real(mysize)

        where(clockg==0.0) clockg=1.e-8
        where(clockgmin==0.0) clockgmin=1.e-8
        where(clockgmax==0.0) clockgmax=1.e-8

        IF(MYRANK==0) THEN
          OPEN(UNIT=16,FILE='clock.dat',FORM='FORMATTED'
     &          ,POSITION='APPEND')
          write(16,'(A,I2,A)') '---- geom: icom=',icom,
     &       '---------------------------------------------------------'
          write(16,'(A)') '    Task/Time           Ave. (sec/%)        Max. (sec/%)   
     &     Min (sec/%)'
          write(16,'(A,3(6x,F8.4))') 'Total               :',clockg(1),clockgmax(1),clockgmin(1)
          i=2
          write(16,905) 'Bndpts points       :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=3
          write(16,905) 'Initialize ANN      :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=11
          write(16,905) 'Find near. neighbor.:'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=5
          write(16,905) 'Find intrs./stencil :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=12
          write(16,905) '   -Norm. intrs.    :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          write(16,'(A,I6,A)') '    (No. of calls=',nnrmg,')'
          i=16
          write(16,905) '   -Norm. intrs.rchk:'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          write(16,'(A,I6,A)') '    (No. of calls=',nnrmrchkg,')'
          i=13
          write(16,905) '   -Grid intrs.     :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          write(16,'(A,I6,A)') '    (No. of calls=',ngrdg,')'
          i=17
          write(16,905) '   -Grid intrs.rchk :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          write(16,'(A,I6,A)') '    (No. of calls=',ngrdrchkg,')'
          i=19
          write(16,905) '   -Reorder ord arr.:'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          write(16,'(A,I6,A)') '    (No. of calls=',nordg,')'
          i=20
          write(16,905) '   -Reorder innf    :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          i=14
          write(16,905) '   -Refresh flag    :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          i=18
          write(16,905) '   -Mrk2flag        :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          i=15
          write(16,905) '   -Redefine bndpt  :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(5),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(5),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(5),'%)'
          i=21
          write(16,905) 'Reorder forc. pts   :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=8
          write(16,905) 'Delete ANN structure:'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=9
          write(16,905) 'MPI reduce opertns. :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          i=10
          write(16,905) 'I/O opertns.        :'
     &         ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &         ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &         ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
          write(16,'(2A)') '-----------------------------------------',
     &         '-----'
          close(16)
        ENDIF

      ENDIF

      deallocate(ntri,nntri)
      deallocate(clock,clockg,clockgmin,clockgmax)

      return

 905  format(A,3(1x,F8.4,A,F6.2,A,2x))

      end
c-----------------------------------------------------------------------


C---- subroutine physical_mrk2flag1_mod------N . Beratlis-22 Aug. 2009--
C----------------------------------------------A. Posa - 5 Nov 2012-----
C
C     PURPOSE: Set flag to physical value (0) based on value of array mrk
C
C-----------------------------------------------------------------------
      subroutine physical_mrk2flag1_mod(flag,nx,ny,nz,mrk,iim,jim,kim,ord,nq,lim,mim,ilu,iru,indexl,indexr)

c      implicit none
      include 'common.h'
      include 'immersed.h'
      
      integer nx,ny,nz,lim,mim,nq
      integer ord(nq)
      integer mrk(nfcmax),iim(nfcmax),jim(nfcmax),kim(nfcmax)
      integer flag(nx,ny,nz)
      integer ilu,iru,indexl(nq,2),indexr(nq,2)

      integer i,ii

      ilu = 0
      iru = 0

      do ii=lim+1,lim+mim
        i=ord(ii)
c        if(iim(i)==102 .AND. jim(i)==4 .AND. kim(i)==264) then
c          write(6,*) 'inside mrk2flag1',lim,mim,ii,i,iim(i),jim(i),kim(i),mrk(i)
c        endif
        if(mrk(i)>0) then
          flag(iim(i),jim(i),kim(i))=0
          if(kim(i).eq.kz1) then
            ilu = ilu + 1
            indexl(ilu,1) = iim(i)
            indexl(ilu,2) = jim(i)
          endif
          if(kim(i).eq.kz2) then
            iru = iru + 1
            indexr(iru,1) = iim(i)
            indexr(iru,2) = jim(i)
          endif
        endif
      enddo
      
      return

      end
C-----------------------------------------------------------------------

