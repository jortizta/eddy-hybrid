	subroutine mpi_setup(nx,ny,nz,nzg,stat)
	INCLUDE 'common.h'
        INCLUDE 'mpif.h'
!Passed Variables
	integer,intent(in)        :: nx,ny,nz,nzg 
 	integer,intent(out)       :: stat  !0=sucess 
        
!Local Variables
 	logical                       :: belongs(3)
	integer                       :: nxdim, nydim, nzdim,nxp,nyp,nzp
	integer                       :: IOUTL
 	logical                       :: debug = .false.

 	stat=0


!Create Cartesian topology with periodic boundary conditions
 	dims(1)=nxprocs
	dims(2)=nyprocs
 	dims(3)=nzprocs 
	nxp=nx-2
	nyp=ny-2
	nzp=nzg-2
 	periods(1)=.true.
 	periods(2)=.true.
 	periods(3)=.true.
	call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,periods,.true.,comm3d,stat)
 	call MPI_COMM_DUP(comm3d,comm3dp,stat)
 	call MPI_COMM_DUP(comm3d,comm3dl,stat)
 	call MPI_COMM_DUP(comm3d,comm3dc,stat)
 	call MPI_COMM_RANK(comm3d,myid,stat)
 	call MPI_CART_GET(comm3d,3,dims,periods,coords,stat)
	
! find the six neighboring blocks, m1 is shorthand for - 1, p1 for + 1
 	call MPI_CART_SHIFT(comm3d,0,1,nbrx1m1, nbrx1p1, stat)
 	call MPI_CART_SHIFT(comm3d,1,1,nbrx2m1, nbrx2p1, stat)
 	call MPI_CART_SHIFT(comm3d,2,1,nbrx3m1, nbrx3p1, stat)
	
!Create 1D communicators
 !3.2 create X1 sub-communicator
 	belongs(1) = .true.
 	belongs(2) = .false.
 	belongs(3) = .false.
 	call MPI_CART_SUB(comm3d,belongs,commx1,stat)
 	call MPI_COMM_RANK(commx1,rankx1,stat)
 	call MPI_COMM_SIZE(commx1,sizex1,stat)

 !3.3 create X2 sub-communicator
 	belongs(1) = .false.
 	belongs(2) = .true.
 	belongs(3) = .false.
 	call MPI_CART_SUB(comm3D,belongs,commx2,stat)
 	call MPI_COMM_RANK(commx2,rankx2,stat)
 	call MPI_COMM_SIZE(commx2,sizex2,stat)

 !3.4 create X3 sub-communicator
 	belongs(1) = .false.
 	belongs(2) = .false.
 	belongs(3) = .true.
 	call MPI_CART_SUB(comm3D,belongs,commx3,stat)
 	call MPI_COMM_RANK(commx3,rankx3,stat)
 	call MPI_COMM_SIZE(commx3,sizex3,stat)

!Create 2D communicators
 !3.5 create X1-X2 sub-communicator
 	belongs(1) = .true.
 	belongs(2) = .true.
 	belongs(3) = .false.
 	call MPI_CART_SUB(comm3D,belongs,commx1x2,stat)
 	call MPI_COMM_RANK(commx1x2,rankx1x2,stat)
 	call MPI_COMM_SIZE(commx1x2,sizex1x2,stat)
 !3.6 create X2-X3 sub-communicator
 	belongs(1) = .false.
 	belongs(2) = .true.
 	belongs(3) = .true.
 	call MPI_CART_SUB(comm3D,belongs,commx2x3,stat)
 	call MPI_COMM_RANK(commx2x3,rankx2x3,stat)
 	call MPI_COMM_SIZE(commx2x3,sizex2x3,stat)

 !3.7 create X1-X3 sub-communicator
 	belongs(1) = .true.
 	belongs(2) = .false.
 	belongs(3) = .true.
 	call MPI_CART_SUB(comm3D,belongs,commx1x3,stat)
 	call MPI_COMM_RANK(commx1x3,rankx1x3,stat)
 	call MPI_COMM_SIZE(commx1x3,sizex1x3,stat)

 	call MPI_BARRIER(comm3d,stat) !wait for all processes to finish 


!Find indices of subdomain 
 	call MPE_DECOMP1D(nxp,dims(1),coords(1),sx,ex)
  	sx=sx+1
  	ex=ex+1
  	cex=sx+(ex-sx)/2

 	call MPE_DECOMP1D(nyp,dims(2),coords(2),sy,ey)
 	 sy=sy+1
 	 ey=ey+1
 	call MPE_DECOMP1D(nzp,dims(3),coords(3),sz,ez)
 	 sz=sz+1
 	 ez=ez+1

        IOUTL=6

 	write(IOUTL,'(a)') "COORDINATES"
 	write(IOUTL,130) coords(1),' ',coords(2),' ',coords(3)
 	write(IOUTL,'(a)') "LOCAL INDICIES"
 	write(IOUTL,140) 'sx=',sx,' ex=',ex,' sy=',sy,' ey=',ey,' sz=',sz,' ez=',ez
 	write(IOUTL,'(a)') ""
	write(IOUTL,'(a)') "MPI DOMAIN DECOMPOSITION SETUP COMPLETED"
                                                                                                                             
 	stat = 0
 	return

130     FORMAT( 10(1x,i4,a1) )
140     FORMAT( 10(1x,a4,i5) )

	end subroutine mpi_setup



 	subroutine MPE_DECOMP1D(n,numprocs,myid,s,e)
	INCLUDE 'common.h'
        INCLUDE 'mpif.h'

	
!Passed Variables
 	integer :: nlocal
 	integer :: deficit
 	integer :: n, numprocs, s, e

!------------------------------------------------------------------------
!  From the MPE library
!  This file contains a routine for producing a decomposition of a 1-d 
!  array when given a number of processors.  It may be used in "direct" 
!  product decomposition.  The values returned assume a "global" domain 
!  in [1:n]
!
! Code      : tmgd3
! Called in : mpi_setup
! Calls     : --
!------------------------------------------------------------------------

 	nlocal  = n / numprocs
 	s       = myid * nlocal + 1
 	deficit = mod(n,numprocs) 
 	s       = s + min(myid,deficit)

 	if (myid .lt. deficit) then
  	nlocal = nlocal + 1
	endif

 	e       = s + nlocal - 1

 	if (e .gt. n .or. myid .eq. numprocs-1) e = n

 	return
	end
