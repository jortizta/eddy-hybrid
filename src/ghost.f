	subroutine ghost(var,vtype,stat)
!@t
! \textbf{subroutine ghost(var,vtype,stat)}
!@h
!   Description:
!     Updates ghost values for parallel version or does the equivalent for 
!     the serial version.
!@h  
!   Comments:
!     The original boundary conditions are overwritten by ghosting and are
!     reapplied after ghosting is completed.
!@q

 	use ntypes, only: r8
 	use Domain, only: sx,ex,sy,ey,sz,ez
	 use Parameters, only: flow_solver
	#ifdef PARALLEL
	 use dd,     only: neighbor, bd, comm3d,MPI_STATUS_SIZE,gtype 
	#endif
	use IO,     only: IOUT
 	use boundC
 	implicit none
 
!Passed Variables
 real(r8),intent(inout)     :: var(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
 character(len=*)           :: vtype
 integer,intent(out)        :: stat

!Local Variables
 integer                    :: l, m, n, rt 
 integer                    :: err1,err2

 err1=0
 err2=0

#ifdef PARALLEL
 call ghostALL(var,vtype,err1)
#else
!Make all BCs periodic just like ghosting would, they will be corrected 
 !in the rest of the subroutine
 !X1
 var(sx-1,:,:)=var(ex,:,:)
 var(ex+1,:,:)=var(sx,:,:)
 !X2
 var(:,sy-1,:)=var(:,ey,:)
 var(:,ey+1,:)=var(:,sy,:)
 !X3
 var(:,:,sz-1)=var(:,:,ez)
 var(:,:,ez+1)=var(:,:,sz)
#endif

!UPDATE DIRICHLET/NEUMANN BOUNDARY CONDITIONS SINCE GHOSTING MAKES THEM ALL PERIODIC
rt=0
select case(vtype)
 case('u','uc','Ucont')
  l=1
 case('v','vc','Vcont')
  l=2
 case('w','wc','Wcont')
  l=3
 case('p')
  l=4
 case('r','rho','rf')
  l=5
  rt=1
 case('cfluc','fluc')
  l=6
 case('rp')
  l=6
  rt=1
 case('psource')
  l=7
 case('utemp','vtemp','wtemp','rtemp','Colltemp')
  l=7
  !call calc_bound(var,vtype,VB(:,:,7),err1)
 case DEFAULT
  write(IOUT,'(a)') "BC TYPE NOT AVAILABLE (ghost.f90): "//trim(vtype)
end select

 select case(flow_solver)
  case('explicit_RK3')  ! BCs for staggered grid 
    do n=1,3            !1 is X1, 2 is X2, 3 is X3
      do m=1,2          !1 is min, 2 is max
        if ( TB(n,1,l).NE.3) then
          call bound_update_stag(var,TB(n,m,l),m,n,VB(n,m,l),vtype,rt,err2 )
        elseif ( TB(n,1,l).EQ.3.AND.m.EQ.1 ) then ! Periodic boundary, no need to call twice
          call bound_update_stag(var,TB(n,m,l),m,n,VB(n,m,l),vtype,rt,err2 )
        endif
#ifdef PARALLEL
        call MPI_BARRIER(comm3d,err1)
        if (err1.NE.0.or.err2.NE.0) goto 9999
#endif
      enddo
    enddo

  case('mixed_RK3_ADI_PC_coll')
    do n=1,3            !1 is X1, 2 is X2, 3 is X3
      do m=1,2          !1 is min, 2 is max
        if ( TB(n,1,l).NE.3) then
          call bound_update_coll(var,TB(n,m,l),m,n,VB(n,m,l),vtype,rt,err2 )
        elseif ( TB(n,1,l).EQ.3.AND.m.EQ.1 ) then ! Periodic boundary, no need to call twice
          call bound_update_coll(var,TB(n,m,l),m,n,VB(n,m,l),vtype,rt,err2 )
        endif
#ifdef PARALLEL
        call MPI_BARRIER(comm3d,err1)
        if (err1.NE.0.or.err2.NE.0) goto 9999
#endif
      enddo
    enddo

  case DEFAULT
    write(IOUT,'(a)') "ABORTING! BCS FOR FLOW SOLVER "//trim(flow_solver)//" NOT IMPLEMENTED"
    err1=1
    goto 9999
 end select

 !normal exit
 stat=max(err1,err2)
 return

 !Unclean
 9999 continue
 write(IOUT,'(a,2(1x,i4))') "ERROR (ghost.f90) for var: "//trim(vtype)," err1/err2: ",err1,err2
 stat=max(err1,err2)
 return

end subroutine ghost
