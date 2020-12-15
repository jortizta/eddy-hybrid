	subroutine write_hybrid(nx,ny,nz,icycle,time,dtm1,xc,xu,yc,yv,zcg,zwg
     &  ,UO,VO,WO,P,DENS,index1)

! This routine saves the variables without centering them and in double precision
        
	include 'common.h'
        include 'mpif.h'
	integer :: nx,ny,nz,icycle,iu,iv,iw,dir,index1,index2,iProc,dummy
	real	:: dtm1,time
	real    :: xu(nx),yv(ny)
	real    :: xc(nx),yc(ny)
	real    :: zcg(nz),zwg(nz)
	real    :: UO(nx,ny,nz),VO(nx,ny,nz),WO(nx,ny,nz),DENS(nx,ny,nz)
        real    :: P(nx,ny,nz) 
        character(len=600) :: filename        
        real :: planeU(nx,ny),planeV(nx,ny),planeW(nx,ny),planeD(nx,ny), planeP(nx,ny)
	
        dir=3

        if(index1.gt.sz.and.index1.lt.ez.or.index1.eq.sz.or.index1.eq.ez)

     &  then 

        iProc=1

        endif
 
        if(iProc.eq.1) then

        index2=index1-(sz-2)
        planeU(:,:)=UO(:,:,index2)     
        planeV(:,:)=VO(:,:,index2)     
        planeW(:,:)=WO(:,:,index2)     
        planeP(:,:)=P(:,:,index2)     
        planeD(:,:)=DENS(:,:,index2)     

        write(filename,'(2a,i4.4,a,i8.8,a)')
     &  trim(path_feed),"/hybrid_k",
     &  index1,"_n", icycle,".interp"

        open(210,file=filename,form='unformatted',status='unknown',readonly)
        write(210) icycle,TIME,dtm1,grav,1.0,1/ru1,1.0
        write(210) dir,index1
        write(210) zcg(index1), zwg(index1)
        write(210) nx, ny
        write(210) xc, xu
        write(210) yc, yv
        write(210)  planeU
        write(210)  planeV
        write(210)  planeW
        write(210)  planeP
        write(210)  planeD

        close(210)

        write(*,*)"*************** HYB plane ********************"
        endif

        return
        end

        subroutine hybrid_init(uo,vo,wo,dens,xc,yc,nx,ny,nz,
     &  planeU,planeV,planeW,planeD)
        
        use density_bg
        include 'common.h'
        include 'mpif.h'

        integer nx,ny,nz,i,j,k,icycle,s1,dummy
        character(len=600) :: filename        
        real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),dens(nx,ny,nz)
        real xc(nx),yc(ny)

        integer :: nxRead,nyRead,nzRead,icycleRead,dirRead,index1,indexRead
	real	:: dtm1Read,timeRead
	real,allocatable,dimension(:) :: xuRead,yvRead
	real,allocatable,dimension(:) :: xcRead,ycRead
	real :: zcgRead,zwgRead
        real(8) planeU(nx,ny), planeV(nx,ny), planeW(nx,ny),
     &  planeP(nx,ny),planeD(nx,ny)
        real(8) :: gravRead,rhoRead,ReRead,PrRead	

        allocate(dens_bg(nx,ny))

        write(filename,'(2a,i4.4,a,i8.8,a)')
     &  trim(path_feed),"/hybrid_k",
     &  index_feed,"_n", start_feed,".interp"

        IF (myrank==0) THEN

        allocate(xuRead(nx),yvRead(ny),xcRead(nx),ycRead(ny))

        open(211,file=trim(filename),form='unformatted',
     &  status='old',action='read',iostat=s1)
        read(211) icycleRead
        read(211) timeRead,dtm1Read,gravRead,rhoRead,ReRead,PrRead
        read(211) dirRead,dummy,indexRead
        read(211) zcgRead, zwgRead
        read(211) nxRead,dummy, nyRead
        read(211) xcRead, xuRead
        read(211) ycRead, yvRead
        read(211) planeU
        read(211) planeV
        read(211) planeW
        read(211) planeP
        read(211) planeD

        close(211)

        deallocate(xuRead,yvRead,xcRead,ycRead)        


        !planeU=0.0
        !planeV=0.0
        !planeW=1.0

        ENDIF
   
        CALL MPI_BCAST(planeU(1,1),nx*ny,mpi_double_precision,0,mpi_comm_eddy,ierr)
        CALL MPI_BCAST(planeV(1,1),nx*ny,mpi_double_precision,0,mpi_comm_eddy,ierr)
        CALL MPI_BCAST(planeW(1,1),nx*ny,mpi_double_precision,0,mpi_comm_eddy,ierr)
        CALL MPI_BCAST(planeD(1,1),nx*ny,mpi_double_precision,0,mpi_comm_eddy,ierr)

        uo=0.0d0
        vo=0.0d0       
        wo=0.0d0
        dens=0.0d0        

        do j=1,ny
        do i=1,nx
        uo(i,j,:) = planeU(i,j)
        vo(i,j,:) = planeV(i,j)
        wo(i,j,:) = planeW(i,j)
        enddo
        enddo


        if(idens.eq.1) then

        do j=1,ny
        do i=1,nx
        dens(i,j,:)   =  planeD(i,j)
        dens_bg(i,j)  =  denP1*(rp(i))*sin(yc(j))
        enddo
        enddo

        else
        dens=0.0d0
        dens_bg=0.0d0
        endif

        if (myrank==0) then
            write(*,*) "    hybrid field initialized"
        endif

         CALL BOUNDARY_DENS(DENS,XC,YC,NX,NY,NZ)
         CALL REFRESHBC(DENS,NX*NY,NZ)        

        return
        end 

	subroutine read_hybrid(filename,planeU,planeV,planeW,planeD,nx,ny)
        
	include 'common.h'
        include 'mpif.h'
        integer :: nx,ny,icycle,s1,status,dummy
	integer :: nxRead,nyRead,nzRead,icycleRead,dirRead,index1,indexRead
	real	:: dtm1Read,timeRead
	real,allocatable,dimension(:) :: xuRead,yvRead
	real,allocatable,dimension(:) :: xcRead,ycRead
	real :: zcgRead,zwgRead
        character(len=600) :: filename        
        real(8) planeU(nx,ny),planeV(nx,ny),planeW(nx,ny),planeP(nx,ny),planeD(nx,ny)
        real(8) :: gravRead,rhoRead,ReRead,PrRead	


        IF (myrank==0) THEN

        allocate(xuRead(nx),yvRead(ny),xcRead(nx),ycRead(ny))

        open(211,file=trim(filename),form='unformatted',
     &  status='old',action='read',iostat=s1)
        read(211) icycleRead
        read(211) timeRead,dtm1Read,gravRead,rhoRead,ReRead,PrRead
        read(211) dirRead,dummy,indexRead
        read(211) zcgRead, zwgRead
        read(211) nxRead,dummy,nyRead
        read(211) xcRead, xuRead
        read(211) ycRead, yvRead
        read(211) planeU
        read(211) planeV
        read(211) planeW
        read(211) planeP
        read(211) planeD

        close(211)

        deallocate(xuRead,yvRead,xcRead,ycRead)        


        write(*,*)"*************** read HYB plane ********************"

        ENDIF
   
!        CALL MPI_BCAST(planeU(1,1),nx*ny,mpi_double_precision,0,mpi_comm_eddy,ierr)
!        CALL MPI_BCAST(planeV(1,1),nx*ny,mpi_double_precision,0,mpi_comm_eddy,ierr)
!        CALL MPI_BCAST(planeW(1,1),nx*ny,mpi_double_precision,0,mpi_comm_eddy,ierr)
!        CALL MPI_BCAST(planeD(1,1),nx*ny,mpi_double_precision,0,mpi_comm_eddy,ierr)

        return
        end

	subroutine sequence_hybrid(planeU,planeV,planeW,planeD,
     &  planeURef,planeVRef,planeWRef,planeDRef,
     &  icycle,nx,ny)

        
	include 'common.h'
        include 'mpif.h'
        integer :: nx,ny,icycle,index1
        character(len=600) :: filename,filenameRef        
        real :: planeU(nx,ny),planeV(nx,ny),planeW(nx,ny),planeD(nx,ny)
        real :: planeURef(nx,ny),planeVRef(nx,ny),planeWRef(nx,ny),planeDRef(nx,ny)
        integer :: num_feed, n_feed_cycles,icycle_feed
        integer :: i,j
 
        num_feed=(end_feed-start_feed)/stride_feed+1
        
        icycle_feed=start_feed+icycle
         
        n_feed_cycles=FLOOR(dble(icycle)/dble(num_feed))

        icycle_feed=INT(start_feed + stride_feed*icycle
     &  -stride_feed*n_feed_cycles*num_feed)
     
     
        IF(MYRANK==0) write(*,*) 'HYB Plane #:',icycle_feed
        
        write(filename,'(2a,i4.4,a,i8.8,a)')
     &  trim(path_feed),"/hybrid_k",
     &  index_feed,"_n", icycle_feed,".interp"


        call read_hybrid(filename,planeU,planeV,planeW,planeD,nx,ny)

        do j=1,ny
        do i=116,nx  !from r=1 impose the steady state IGW field

        planeU(i,j) = planeURef(i,j)
        planeV(i,j) = planeVRef(i,j)
        planeW(i,j) = planeWRef(i,j)
        planeD(i,j) = planeDRef(i,j)

        enddo
        enddo

      
        return
        end

