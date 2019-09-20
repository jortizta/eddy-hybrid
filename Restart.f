        subroutine read_flow(in_file,in_iter,stat)
!@t
! \textbf{subroutine read\_flow(in\_file,in\_iter,stat)}
!@h
!   Description:
!     Reads in a restart file.
!@q
!   Modification History
!     Version   Date     Comment 
!     -------   ----     ------- 
!     1.0       07/2008  Original code. [Kyle A. Brucker] 

        use IO, only: IObig

 	implicit none

 !Passed Variables
 	character(len=*),intent(in) :: in_file
 	integer,intent(in)          :: in_iter
 	integer,intent(out)         :: stat

 !Local Variables
 	integer                     :: ok 

 	ok=0
	#ifdef PARALLEL
 	if (IObig) then
  	call start_big(in_file,in_iter,ok)
 	else
  	call start_small(in_file,in_iter,ok)
 	endif
	#else
 	call start_small(in_file,in_iter,ok)
	#endif

 	stat=ok

	return
	end subroutine read_flow

	subroutine write_flow(out_file,stat)
!@t
! \textbf{subroutine write\_flow(out\_file,stat)}
!@h
!   Description:
!     Writes a restart file.
!@q
!   Modification History
!     Version   Date     Comment 
!     -------   ----     ------- 
!     1.0       07/2008  Original code. [Kyle A. Brucker] 

 	use IO, only: IObig
 	use Parameters, only: Spatial
 	implicit none

 !Passed Variables
 	character(len=*),intent(in) :: out_file
 	integer,intent(out)         :: stat

 !Local Variables
 	integer                     :: ok 
 	ok=0

	#ifdef PARALLEL
    	call write_flow_big(out_file,ok)
    	if(Spatial) then
    	call write_flow_small(out_file,ok)
    	endif
!  if (IObig) then
!     call write_flow_big(out_file,ok)
!  elseif(IObig.and.Spatial) then
!     call write_flow_big(out_file,ok)
!     call write_flow_small(out_file,ok)
!  else
!     call write_flow_small(out_file,ok)
!  endif
	#else
 	call write_flow_small(out_file,ok)
	#endif

 	stat=ok

	return
	end subroutine write_flow

	subroutine Greadfile(ICfile,u,v,w,p,r,nx,ny,nz,  &
                         n_time,time,delt,g,rho_0,Re,Pr,stat)
!@t
! \textbf{subroutine Greadfile(ICfile,u,v,w,p,r,nx,ny,nz,n\_time,time,delt,g,rho\_0,Re,Pr,stat)}
!@h
!   Description:
!     Reads in a restart file line by line.
!@q
!   Modification History
!     Version   Date     Comment 
!     -------   ----     ------- 
!     1.0       07/2008  Original code. [Kyle A. Brucker] 

 	use ntypes, only: r8
 	use IO,     only: IOUT
 	implicit none

!Passed Variables
 	character(len=*),intent(in)     :: ICfile
 	integer,intent(in)              :: nx, ny, nz
 	integer,intent(out)             :: n_time
 	real(r8),intent(out)            :: time,delt
 	real(r8),intent(out)            :: g,rho_0,Re,Pr
 	real(r8),intent(out)            :: u(1:nx,1:ny,1:nz)
 	real(r8),intent(out)            :: v(1:nx,1:ny,1:nz)
 	real(r8),intent(out)            :: w(1:nx,1:ny,1:nz)
 	real(r8),intent(out)            :: p(1:nx,1:ny,1:nz)
 	real(r8),intent(out)            :: r(1:nx,1:ny,1:nz)
 	integer,intent(out)             :: stat

!Local Variables
 	integer                     :: n_r(1:3), s1

!Direct Access read
 	open(310,file=ICfile,form='unformatted',status='old',iostat=s1)
  	if (s1.NE.0) then
   	write(IOUT,'(a,i4)') "ERROR OPENING FILE: "//trim(Icfile)//" IOSTAT=",s1
   	stat=1
   	goto 2000
  	endif
  	read(310) n_time,time,delt,g,rho_0,Re,Pr
  	read(310) n_r 
   	if (n_r(1).NE.nx)   goto 1000 !Check to make sure dump is correct size
   	if (n_r(2).NE.ny)   goto 1000
   	if (n_r(3).NE.nz)   goto 1000
 	read(310) u 
  	read(310) v
  	read(310) w
  	read(310) p
  	read(310) r 
 	close(310)

 	write(IOUT,'(a)') "READ OF INITIAL FIELD: "//trim(Icfile)//" COMPLETED"
 	return

 	1000 continue
 	close(310)
 	write(IOUT,'(a36,(3x,i4,a1))') "ERROR: Data file size is [nx,ny,nz]", &
                                      n_r(1),' ',n_r(2),' ',n_r(3)
 	write(IOUT,'(a36,(3x,i4,a1))') "ERROR: Run  size is [nx,ny,nz]", nx, &
                                   ' ',ny,' ',nz
 	stat = 1
 	stop

 	2000 continue
 	close(310)
 	stat = 2
 	stop
	end subroutine Greadfile

	#ifdef PARALLEL
	subroutine start_big(basename,ntime,err1)
!@t
! \textbf{subroutine start\_big(basename,ntime,err1)}
!@h
!   Description:
!     Starts a simulation with data from multiple distributed restart files.
!@q
!   Modification History
!     Version   Date     Comment 
!     -------   ----     ------- 
!     1.0       07/2008  Original code. [Kyle A. Brucker] 

 	use ntypes,     only: r8
 	use Flow,       only: u,v,w,p,rho
 	use Domain,     only: sx,ex,sy,ey,sz,ez
 	use Parameters, only: g, rho_0, Re, Pr, time, delt, nstep, rRe, rPr, icparam
 	use IO,         only: IOUT
 	use dd,         only: coords
 	implicit none

!Passed Variables
 	character(len=*)       :: basename 
 	integer,intent(in)     :: ntime
 	integer,intent(out)    :: err1

!Local Variables
 	integer  :: s1,stat, n_r(1:6),i
 	real(r8) :: g_IN, rho_0_IN, Re_IN, Pr_IN, time_IN, delt_IN
 	integer  :: nstep_IN 
 	character (len = 300 )  :: ICfile 

 	ICfile=basename
 	call concati(ICfile,coords(1))
 	call concat(ICfile,'_')
 	call concati(ICfile,coords(2))
 	call concat(ICfile,'_')
 	call concati(ICfile,coords(3))
 	call concat(ICfile,'.')
 	call concati(ICfile,ntime)

!Direct Access read
 	open(310,file=ICfile,form='unformatted',status='old',iostat=s1,action='read')
  	if (s1.NE.0) then
   	write(IOUT,'(a,i4)') "ERROR OPENING FILE: "//trim(Icfile)//" IOSTAT=",s1
   	stat=1
   	stop
  	endif

  	read(310) nstep_IN,time_IN,delt_IN,g_IN,rho_0_IN,Re_IN,Pr_IN
  	read(310) n_r
   	if (n_r(1).NE.sx)   goto 1000 !Check to make sure dump is correct size
   	if (n_r(2).NE.ex)   goto 1000
   	if (n_r(3).NE.sy)   goto 1000
   	if (n_r(4).NE.ey)   goto 1000 
   	if (n_r(5).NE.sz)   goto 1000
   	if (n_r(6).NE.ez)   goto 1000
  	read(310) u
  	read(310) v
  	read(310) w
  	read(310) p
  	read(310) rho
 	close(310)
 
 	if (icparam) then
  	Re=Re_IN
  	rRe = 1.d0/Re
  	Pr=Pr_IN
 	 rPr=1.d0/Pr
  	nstep=nstep_IN
  	time=time_IN
  	delt=delt_IN
 	 g=g_IN
	rho_0=rho_0_IN
	endif
 	err1=0
 	write(IOUT,'(a)') "READ OF INITIAL FIELD: "//trim(ICfile)//" COMPLETED"
 	return
                                                                                                                             
 	1000 continue
 	close(310)
 	write(IOUT,'(a36,3(1x,i4))') "ERROR: Data file size is [sx:ex,sy:ey,ez:ez]", &
                                      (n_r(i),i=1,6)
 	write(IOUT,'(a36,3(1x,i4))') "ERROR: Run  size is [sx:ex,sy:ey,sz:ez]",sx,ex,&
                                  sy,ey,sz,ez
 
 	close(unit=310)
 	stop
	end subroutine start_big
	#endif

	subroutine start_small(basename,ntime,err1)
!@t
! \textbf{subroutine start\_small(basename,ntime,err1)}
!@h
!   Description:
!     Starts a simulation from one large restart file.
!@q
!   Modification History
!     Version   Date     Comment 
!     -------   ----     ------- 
!     1.0       07/2008  Original code. [Kyle A. Brucker] 

 	use ntypes,     only: r8
 	use Flow,       only: u,v,w,p,rho
 	use Domain,     only: nxp2,nyp2,nzp2,sx,ex,sy,ey,sz,ez
 	use Parameters, only: g, rho_0, Re, Pr, time, delt, nstep, rRe, rPr, icparam
 	use IO,         only: IOUT
#ifdef PARALLEL
 	use dd,         only: myid,comm3d,realtype, inttype
#endif
 	implicit none

 !Passed Variables
 	character(len=*)       :: basename
 	integer,intent(in)     :: ntime 
 	integer,intent(out)    :: err1

 !Local Variables
 	integer  :: s1,stat
 	real(r8) :: g_IN, rho_0_IN, Re_IN, Pr_IN, time_IN, delt_IN
 	integer  :: nstep_IN 
 	character (len = 300 )  :: ICfile

#ifdef PARALLEL
	real(r8),allocatable,dimension(:,:,:) :: uF,vF,wF,pF,rF,tmp1
 	if (myid.eq.0) then
  	call deallocate_temps(stat)
  	if (err1.NE.0) then
   	write(IOUT,'(a,i4)') "ERROR: Deallocating Temps, stat:",s1
   	stat=1
   	goto 1000
  	endif
  	allocate( uF(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  	allocate( vF(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  	allocate( wF(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  	allocate( pF(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  	allocate( rF(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  	allocate( tmp1(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1), stat=s1 )
  	if (s1.NE.0) then
   	write(IOUT,'(a,i4)') "ERROR: Allocating full fields, stat:",s1
   	stat=2
  	goto 1000
  	endif

  	ICfile=basename
  	call concati(ICfile,ntime)

  	call Greadfile(ICfile,uF,vF,wF,pF,rF,nxp2,nyp2,nzp2,&
                  nstep_IN,time_IN,delt_IN,g_IN,rho_0_IN,Re_IN,Pr_IN,stat)

 	endif 

 	if (myid.eq.0) then
   	call Distribute3dM(uF,u,tmp1,0,stat)
 	else
   	call Distribute3dS(u,0,stat)
 	endif
 	write(IOUT,*) "U1 DISTRIBUTED"
 	call MPI_BARRIER(comm3d,stat)

 	if (myid.eq.0) then
  	call Distribute3dM(vF,v,tmp1,0,stat)
 	else
  	call Distribute3dS(v,0,stat)
 	endif
 	write(IOUT,*) "U2 DISTRIBUTED"
 	call MPI_BARRIER(comm3d,stat)

 	if (myid.eq.0) then
  	call Distribute3dM(wF,w,tmp1,0,stat)
 	else
  	call Distribute3dS(w,0,stat)
 	endif
 	write(IOUT,*) "U3 DISTRIBUTED"
 	call MPI_BARRIER(comm3d,stat)

 	if (myid.eq.0) then
  	call Distribute3dM(pF,p,tmp1,0,stat)
 	else
  	call Distribute3dS(p,0,stat)
 	endif
 	write(IOUT,*) "P DISTRIBUTED"
 	call MPI_BARRIER(comm3d,stat)

 	if (myid.eq.0) then
  	call Distribute3dM(rF,rho,tmp1,0,stat)
 	else
  	call Distribute3dS(rho,0,stat)
 	endif
 	write(IOUT,*) "RHO DISTRIBUTED"
 	call MPI_BARRIER(comm3d,stat)


 	if (myid.EQ.0) then
   	deallocate( uF, stat=s1 )
   	deallocate( vF, stat=s1 )
   	deallocate( wF, stat=s1 )
   	deallocate( pF, stat=s1 )
   	deallocate( rF, stat=s1 )
   	deallocate( tmp1, stat=s1 )
   	if (s1.NE.0) then
   	 write(IOUT,'(a,i4)') "ERROR: Deallocating full fields in  &
                  & Distribute3d, stat:",s1

    	stat=2
    	goto 1000
   	endif
   	call allocate_temps(stat)
  	if (err1.NE.0) then
   	write(IOUT,'(a,i4)') "ERROR ALLOCATING TEMPS IN `startup` FATAL!"
   	goto 2000
  	endif
 	endif

 
 	if (icparam) then
  	if (myid.eq.0) then
   	Re=Re_IN
   	rRe = 1.d0/Re
   	Pr=Pr_IN
   	rPr = 1.d0/Pr
   	nstep=nstep_IN
   	time=time_IN
   	delt=delt_IN
   	g=g_IN
   	rho_0=rho_0_IN
  	endif

  	call MPI_BCAST(Re,1,realtype,0,comm3d,stat)
  	rRe=1.d0/Re
  	call MPI_BCAST(Pr,1,realtype,0,comm3d,stat)
  	rPr=1.d0/Pr
  	call MPI_BCAST(nstep,1,inttype,0,comm3d,stat)
  	call MPI_BCAST(time,1,realtype,0,comm3d,stat)
  	call MPI_BCAST(delt,1,realtype,0,comm3d,stat)
  	call MPI_BCAST(g,1,realtype,0,comm3d,stat)
  	call MPI_BCAST(rho_0,1,realtype,0,comm3d,stat)
 	endif

	#else

 	ICfile=basename
 	call concati(ICfile,ntime)

 	call Greadfile(ICfile,u,v,w,p,rho,nxp2,nyp2,nzp2,&
                  nstep_IN,time_IN,delt_IN,g_IN,rho_0_IN,Re_IN,Pr_IN,stat)
 	if (icparam) then
  	Re=Re_IN
  	rRe=1.d0/Re
  	Pr=Pr_IN
  	rPr=1.d0/Pr
  	nstep=nstep_IN
  	time=time_IN
  	delt=delt_IN
  	g=g_IN
  	rho_0=rho_0_IN
  	rRe=1.d0/Re
  	rPr=1.d0/Pr
 	endif

#endif

 	call ghost(u,'u',stat     )
 	call ghost(v,'v',stat     )
 	call ghost(w,'w',stat     )
 	call ghost(p,'p',stat     )
 	call ghost(rho,'rho',stat )

 	write(IOUT,'(a)') "START FROM SERIAL DUMP COMPLETED" 
 	err1=stat
 	return
 	1000 continue
 	write(IOUT,'(a)') "START FROM SERIAL DUMP FAILED" 
	 err1=stat
 	return
 
 	2000 stop
	end subroutine start_small

#ifdef PARALLEL
	subroutine write_flow_big(basename,err1)
!@t
! \textbf{subroutine write\_flow\_big(basename,err1)}
!@h
!   Description:
!     Write a restart file for each processor for a parallel simulation.
!@q
!   Modification History
!     Version   Date     Comment 
!     -------   ----     ------- 
!     1.0       07/2008  Original code. [Kyle A. Brucker] 

 	use ntypes,     only: r8
 	use Flow,       only: u,v,w,p,rho
 	use Domain,     only: sx,ex,sy,ey,sz,ez
 	use Parameters, only: g, rho_0, Re, Pr, time, delt, nstep
 	use IO,         only: IOUT,bigDIR
 	use dd,         only: coords
 	implicit none

!Passed Variables
 	character(len=*),intent(in)       :: basename 
 	integer,intent(out)               :: err1

!Local Variables
 	integer  :: s1
 	character (len = 250 )            :: ICfile 
 	logical,parameter                 :: debug=.false.

 	if (debug) call check_point('write_flow_big#1',.false.)


 	err1=0

 	ICfile=bigDIR
 	call concat(ICfile,basename)
 	call concati(ICfile,coords(1))
 	call concat(ICfile,'_')
 	call concati(ICfile,coords(2))
 	call concat(ICfile,'_')
 	call concati(ICfile,coords(3))
 	call concat(ICfile,'.')
 	call concati(ICfile,nstep)
                                                                                                                             
!Direct Access Fortran Binary
 	open(310,file=ICfile,form='unformatted',status='new',iostat=s1,action='write')
  	if (s1.NE.0) then
  	 write(IOUT,'(a,i4)') "ERROR OPENING FILE: "//trim(Icfile)//" IOSTAT=",s1
   	err1=1
   	goto 2000
  	endif
   	if (debug) call check_point('write_flow_big#2',.false.)

  	write(310) nstep,time,delt,g,rho_0,Re,Pr
  	write(310) sx,ex,sy,ey,sz,ez 
  	write(310) u
  	write(310) v
  	write(310) w
  	write(310) p
  	write(310) rho
 	close(310,iostat=s1)
  	if (s1.NE.0) then
   	write(IOUT,'(a,i4)') "ERROR CLOSING FILE: "//trim(Icfile)//" IOSTAT=",s1
   	err1=1
   	goto 2000
  	endif

    	if (debug) call check_point('write_flow_big#3',.false.)

 	err1=max(err1,s1)
 	write(IOUT,'(a)') "WRITE OF BIG RESTART :"//trim(ICfile)//" COMPLETED"
 	return

 	2000 continue
 	write(IOUT,'(a)') "WRITE OF BIG RESTART :"//trim(ICfile)//" FAILED"
 	return
	end subroutine write_flow_big
#endif

	subroutine write_flow_small(basename,stat)
!@t
! \textbf{subroutine write\_flow\_small(basename,stat)}
!@h
!   Description:
!     Writes one large restart file for multiple processors.
!@q
!   Modification History
!     Version   Date     Comment 
!     -------   ----     ------- 
!     1.0       07/2008  Original code. [Kyle A. Brucker] 

 	use ntypes, only: r8
 	use Flow
 	use Domain
 	use IO, only: IOUT, smallDIR,flowDIR
 	use Parameters, only: Re, Pr, g, rho_0, time, nstep, delt,displn
#ifdef PARALLEL
 	use dd, only: comm3d, myid, MPI_STATUS_SIZE
#endif
 	implicit none

 !Passed Variables
 	integer,intent(out)                   :: stat
 	character(len=*),intent(in)           :: basename

 !Local Variables
 	real(4),allocatable,dimension(:,:,:) :: uF,vF,wF,pF,rhoF,tmpu,tmpv,tmpw,tmpp,tmprho,tmp1
 	character(len=250)                     :: ICfile
 	integer                               :: s1, err,ix,i,j,nxr,n
 	logical,parameter                     :: debug=.false.

#ifdef PARALLEL
 	stat=0
 	err=0
 	s1=0
   	if (debug) call check_point('write_flow_small#mpi_1',.false.)

 	if (myid.eq.0) then
  	call deallocate_temps(err)
  	if (err.NE.0) then
   	write(IOUT,'(a)') "ERROR: DEALLOCATION OF TEMPS IN `write_flow` FAILED"
   	stat=1
   	goto 1000
  	endif


  	allocate( tmpu(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  	allocate( tmpv(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  	allocate( tmpw(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  	allocate( tmpp(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  	allocate( tmprho(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  
  	nxr=3.0d0*floor((dble(nxp2)/dble(displn))+1)

  	allocate( uF(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  	allocate( vF(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  	allocate( wF(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  	allocate( pF(1:nxp2,1:nyp2,1:nzp2), stat=s1 )
  	allocate( rhoF(1:nxp2,1:nyp2,1:nzp2), stat=s1 )

  
!   allocate( uF(1:nxr,1:nyp2,1:nzp2), stat=s1 )
!   allocate( vF(1:nxr,1:nyp2,1:nzp2), stat=s1 )
!   allocate( wF(1:nxr,1:nyp2,1:nzp2), stat=s1 )
!   allocate( pF(1:nxr,1:nyp2,1:nzp2), stat=s1 )
!   allocate( rhoF(1:nxr,1:nyp2,1:nzp2), stat=s1 )


  	if (s1.NE.0) then
   	write(IOUT,'(a,i4)') "ERROR: Allocating full field, stat:",s1
   	stat=2
   	goto 1000
  	endif

  	allocate( tmp1(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1), stat=s1 )
  	if (s1.NE.0) then
  	 write(IOUT,'(a,i4)') "ERROR: Allocating temp field, stat:",s1
   	stat=2
  	 goto 1000
  	endif

  	write(ICfile,'(a,i5.5,a)') trim(flowDIR)//"myout",nstep
  	open(310,file=ICfile,form='unformatted',status='unknown',iostat=s1)
  	if (s1.NE.0) then
   	write(IOUT,'(a,i4)') "ERROR OPENING FILE: "//trim(ICfile)//" IOSTAT: ",s1
   	stat=3
   	goto 1000
  	endif
  	write(310) nstep,time,delt,g,rho_0,Re,Pr
  	write(310) nxp2, nyp2, nzp2, nxr 
 	endif
   	if (debug) call check_point('write_flow_small#mpi_2',.false.)
 


 !Gather U
  	if (myid.eq.0) then
   	call gather3d_StatM(u,tmpu,0,err)
   	call gather3d_StatM(v,tmpv,0,err)
   	call gather3d_StatM(w,tmpw,0,err)
   	call gather3d_StatM(p,tmpp,0,err)
   	call gather3d_StatM(rho,tmprho,0,err)
  	else
	call gather3d_statS(u,0,err)
   	call gather3d_statS(v,0,err)
   	call gather3d_statS(w,0,err)
   	call gather3d_statS(p,0,err)
   	call gather3d_statS(rho,0,err)
  	endif
   	call MPI_BARRIER(comm3d,err)
   	if (debug) call check_point('write_flow_small#mpi_3',.false.)

!!!!!!!!!!!!!!!!!!!!!!Testing low storage!!!!!!!!!!!!!!!!!!!!!!!!

	uF(:,:,:)=tmpu(:,:,:)
	vF(:,:,:)=tmpv(:,:,:)
	wF(:,:,:)=tmpw(:,:,:)
	pF(:,:,:)=tmpp(:,:,:)
	rhoF(:,:,:)=tmprho(:,:,:)
! 
! n=1
! do i=1,nxr,3
! ix=1+(i-n)*displn
! do j=1,nyp2
! uF(i,j,:)=tmpu(ix,j,:)
! vF(i,j,:)=tmpv(ix,j,:)
! wF(i,j,:)=tmpw(ix,j,:)
! pF(i,j,:)=tmpp(ix,j,:)
! rhoF(i,j,:)=tmprho(ix,j,:)
! enddo
! n=n+2
! enddo
! 
! n=2
! do i=2,nxr,3
! ix=2+(i-n)*displn
! do j=1,nyp2
! uF(i,j,:)=tmpu(ix,j,:)
! vF(i,j,:)=tmpv(ix,j,:)
! wF(i,j,:)=tmpw(ix,j,:)
! pF(i,j,:)=tmpp(ix,j,:)
! rhoF(i,j,:)=tmprho(ix,j,:)
! enddo
! n=n+2
! enddo
! 
! n=3
! do i=3,nxr,3
! ix=3+(i-n)*displn
! do j=1,nyp2
! uF(i,j,:)=tmpu(ix,j,:)
! vF(i,j,:)=tmpv(ix,j,:)
! wF(i,j,:)=tmpw(ix,j,:)
! pF(i,j,:)=tmpp(ix,j,:)
! rhoF(i,j,:)=tmprho(ix,j,:)
! enddo
! n=n+2
! enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !Write U
  	if(myid.eq.0) then 
!    if (displn.eq.1) then
    	write(310) uF 
    	write(310) vF
    	write(310) wF
    	write(310) pF
    	write(310) rhoF
!    else
!     do i=1,nxr
!      do j=1,nyp2
!      write(310) uF(i,j,:) 
!      write(310) vF(i,j,:) 
!      write(310) wF(i,j,:) 
!      write(310) pF(i,j,:) 
!      write(310) rhoF(i,j,:) 
!      enddo
!     enddo
!    endif
 	 endif
  	 call MPI_BARRIER(comm3d,err)
  	 if (debug) call check_point('write_flow_small#mpi_4',.false.)

  	write(IOUT,*) "Planes Written"
  	if(myid.eq.0) then
   	deallocate( uF, stat=s1 )
   	deallocate( vF, stat=s1 )
   	deallocate( wF, stat=s1 )
   	deallocate( pF, stat=s1 )
   	deallocate( rhoF, stat=s1 )
   	if (s1.NE.0) then
    	write(IOUT,'(a,i4)') "ERROR: Deallocating full field, stat:",s1
    	goto 1000
   	endif

   	deallocate( tmp1, stat=s1 )
   	if (s1.NE.0) then
    	write(IOUT,'(a,i4)') "ERROR: Deallocating tmp1, stat:",s1
    	goto 1000
   	endif
   	call allocate_temps(err)
   	if (err.NE.0) then
    	write(IOUT,'(a)') "ERROR ALLOCATING TEMPS IN `write_flow` FATAL!"
    	stat=4
    	goto 2000
   	endif
    	close(310)
    	write(IOUT,'(a)') "WRITE OF SMALL RESTART :"//trim(ICfile)//" COMPLETED"
   	endif

   	call MPI_BARRIER(comm3d,err)
  
#else
 	err=0
 	ICfile=flowDIR
 	call concat(ICfile,basename)
 	call concati(ICfile,nstep)
  	open(310,file=ICfile,form='unformatted',status='unknown',iostat=s1)
  	if (s1.NE.0) then
   	write(IOUT,'(a,i4)') "ERROR OPENING FILE: "//trim(ICfile)//" IOSTAT: ",s1
   	stat=3
   	goto 1000
  	endif
  	write(310) nstep,time,delt,g,rho_0,Re,Pr
  	write(310) nxp2, nyp2, nzp2 
  	write(310) u
  	write(310) v
  	write(310) w
  	write(310) p
  	write(310) rho
 	close(310,iostat=s1)
  	if (s1.NE.0) then
   	write(IOUT,'(a,i4)') "ERROR CLOSING FILE: "//trim(ICfile)//" IOSTAT: ",s1
   	stat=4
   	goto 1000
  	endif

#endif

 	stat=max(err,s1)
 	1000 continue
 	return
 	2000 continue 
 	write(IOUT,'(a,i4)') "WRITE OF SMALL RESTART :"//trim(ICfile)//" FAILED STAT: ",s1
 	return
	end subroutine write_flow_small



	subroutine gather3d_StatM(varLt,OutPlane,myidM,ok)
 	use ntypes, only: r8
 	use dd,     only: myid, commX1X2X3,sizeX1X2X3,coords,nxprocs,nyprocs,nzprocs,&
                   MPI_STATUS_SIZE, inttype, realtype,sizex1,sizex2,sizex3,comm3d
 	use Domain, only: sx,ex,sy,ey,sz,ez,nxp2,nyp2,nzp2
 	use IO,     only: IOUT
 	implicit none
                                                                                                                             
!Passed Variables
 	integer,intent(in)                  :: myidM
 	real(r8),intent(in)                 :: varLt(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
 	real(4),intent(out)                 :: OutPlane(1:nxp2,1:nyp2,1:nzp2)
 	integer,intent(out)                 :: ok
                                                                                                                             
!Local Variables
 	integer                             :: Tsize, Rcoords(3), s1, status1(MPI_STATUS_SIZE),ierr
 	integer                             :: i,j,k,n, istart,jstart, kstart, i2,j2, k2
 	real(r8),allocatable,dimension(:,:,:) :: Temp_Recv
 	integer                             :: is,ie,js,ks,je,ke
 
 	logical,parameter                   :: debug=.false.
 	s1=0
 	ierr=0

                                                                     
 	allocate( Temp_Recv(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1), STAT=s1 )
 	Tsize=size(Temp_recv)
                                                                    
 	do n=0,sizeX1X2X3-1 !1
                                                                   
 	 if (n.Eq.myidM) then
   	Rcoords=coords
   	Temp_Recv(:,:,:)=varLt(:,:,:)
   	else
   	call MPI_RECV(Rcoords,3,inttype,n,2,comm3d,status1,ierr)
  	call MPI_RECV(Temp_Recv,Tsize,realtype,n,1,comm3d,status1,ierr)
  	endif
                                                                  
  !Determine Block of Data to recieve
   	istart = Rcoords(1)*(nxp2-2)/nxprocs
   	jstart = Rcoords(2)*(nyp2-2)/nyprocs
  	 kstart = Rcoords(3)*(nzp2-2)/nzprocs

  !Determine if there is boundary data 
   	is=0
   	ie=0
   	js=0
   	je=0
   	ks=0
   	ke=0
   	if ( Rcoords(1).EQ.0      )   is=1
   	if ( Rcoords(2).EQ.0      )   js=1
   	if ( Rcoords(3).EQ.0      )   ks=1
   	if ( Rcoords(1).EQ.sizex1-1 ) ie=1
   	if ( Rcoords(2).EQ.sizex2-1 ) je=1
   	if ( Rcoords(3).EQ.sizex3-1 ) ke=1

   !UnPack Data
   	do k=sz-ks,ez+ke
    	do j=sy-js,ey+je
     	do i=sx-is,ex+ie
      	i2=istart+i
      	j2=jstart+j
      	k2=kstart+k
      	Outplane(i2,j2,k2)=Temp_Recv(i,j,k)
     	enddo
    	enddo
  	 enddo
  	if (debug) write(6,*) "master:",n,maxval(Temp_Recv)
  	enddo

                                                                                                                            
 	deallocate(Temp_Recv,STAT=s1)
                                                                                                                             
 	ok=max(ierr,s1)
 	return
	end subroutine gather3d_statM


	subroutine gather3d_statS(varLt,myidM,ok)
 	use ntypes, only: r8
 	use dd,     only: comm3d, coords,MPI_STATUS_SIZE, inttype, realtype,myid
	 use domain, only: sx,ex,sz,ez,sy,ey
 	implicit none
                                                                                                                             
!Passed Variables
 	integer,intent(in)        :: myidM
 	real(r8),intent(in)       :: varLt(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
 	integer,intent(out)       :: ok
                                                                                                                             
!Local Variables
 	integer                            :: Tsize,ierr

 	logical,parameter                   :: debug=.false.

 	ierr=0
 	Tsize=size(varLt)
 	call MPI_SEND(coords,3,inttype,myidM,2,comm3d,ierr)
 	call MPI_SEND(varLt,Tsize,realtype,myidM,1,comm3d,ierr)
 	if (debug) write(6,*) "slave:",myid,maxval(varLt)  
	ok=ierr
	return
	end subroutine gather3d_statS
