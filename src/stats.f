c---- subroutine sortprobes -------------------N. Beratlis-06 May 2010---
C
C     PURPOSE: Sort probes in increasing z coordinate.
C
C------------------------------------------------------------------------
      subroutine sortprobes(xp,yp,zp,n)
c
      implicit none
c
c...  Input/Output Arrays
      integer n
      real    xp(n),yp(n),zp(n)
c
c.... Local arrays
      integer i,imin
      real    temp
      real    x2(n),y2(n),z2(n)

      x2 = xp
      y2 = yp
      z2 = zp

      do i=1,n
        imin = i-1 + minloc(z2(i:n),1)
        temp = z2(i)
        z2(i) = z2(imin)
        z2(imin) = temp
        
        temp = x2(i)
        x2(i) = x2(imin)
        x2(imin) = temp

        temp = y2(i)
        y2(i) = y2(imin)
        y2(imin) = temp

      enddo

      zp = z2
      yp = y2
      xp = x2

      return

      end


      subroutine read_probes_all(nz)
c
      use points
      include 'common.h'
      include 'averages.h'
c
c.... Input/Output Arrays
      integer nz,nprobes2
c
c.... Local arrays
      integer i,ii,jj,kk,k,kk1,kk2,
     %iprb2g(nprobes),jprb2g(nprobes),kprb2g(nprobes)
c
      open(unit=10,file='probes.input',form='formatted')
      do i=1,nprobes
	read(10,*) iprb2g(i),jprb2g(i),kprb2g(i)
!      if(MYRANK==0)write(*,*) iprb2g(i),jprb2g(i),kprb2g(i)
      enddo
      close(10)
c
      call sortprobes2(iprb2g,jprb2g,kprb2g,nprobes)
c
      if(myrank.eq.0) then

         open(50,file='probes.dat',form='formatted')

         do i=1,nprobes

            write(50,'(4i5)')i,iprb2g(i),jprb2g(i),kprb2g(i)

         enddo

         close(50)

      endif
c
      if(mysize.eq.1) then
	 kk1=1
	 kk2=nz
      else
	 if(myrank.eq.0) then
	    kk1=1
	    kk2=kz2
	 elseif(myrank.eq.mysize-1) then
	    kk1=kz1
	    kk2=nz
	 else
	    kk1=kz1
	    kk2=kz2
	 endif
      endif
c
      jj=0
      prbindx2=0
      do ii=1,nprobes
	kk=kprb2g(ii)
	k=kk-myrank*(nz-2)
	if(k>=kk1 .AND. k<=kk2) then
	  jj = jj+1
	  if(jj==1) prbindx2=ii
!	   iprb2(jj)=iprb2g(ii)
!	   jprb2(jj)=jprb2g(ii)
!	   kprb2(jj)=k
	endif
      enddo
      nprbmax2 = jj
      nprobes2 = nprbmax2
c
      allocate(iprb2(nprbmax2),jprb2(nprbmax2),kprb2(nprbmax2))
c
      jj=0
      do ii=prbindx2,prbindx2+nprbmax2-1
	jj=jj+1
	iprb2(jj)=iprb2g(ii)
	jprb2(jj)=jprb2g(ii)
	kk=kprb2g(ii)
	k=kk-myrank*(nz-2)
	kprb2(jj)=k
!        write(*,*),"iprb2=",iprb2(jj),"jprb2=",jprb2(jj),"kprb2=",kprb2(jj)
      enddo
c
      return
      end

      subroutine rec_timeprobes_all(uo,vo,wo,p,nx,ny,nz,time)
c
      use points
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
      include 'averages.h'
c
 
c.... Input/Output array
      integer nx,ny,nz,flag,ii,iid,iiu
!       integer iprb(nprobe2),jprb(nprobe2),kprb(nprobe2)
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
      real    uprb,vprb,wprb,pprb,time
c
c.... Local arrays
      integer i
      character*120 filename
c N.B. NPRBMAX is the number of probes inside the domain
c the centered values of velocity and pressure are evaluated
!       write(*,*),"iprb2=",iprb2,"jprb2=",jprb2,"kprb2=",kprb2
      do i=1,nprbmax2
         ii=prbindx2+(i-1)
         iid=ii/10
         iiu=mod(ii,10)
         open(750+ii,file='probes_'//char(48+iid)//char(48+iiu)//'.out',position='append')
!         write(*,*) i,iprb(i),jprb(i),kprb(i)
        uprb = uo(iprb2(i)  ,jprb2(i)  ,kprb2(i))

        vprb = vo(iprb2(i)  ,jprb2(i)  ,kprb2(i))

        wprb = wo(iprb2(i)  ,jprb2(i)  ,kprb2(i))
!         uprb = 0.5*(uo(iprb2(i)  ,jprb2(i)  ,kprb2(i)  )
!      &                     + uo(iprb2(i)-1,jprb2(i)  ,kprb2(i)  ) )
! 
!         vprb = 0.5*(vo(iprb2(i)  ,jprb2(i)  ,kprb2(i)  )
!      &                     + vo(iprb2(i)  ,jprb2(i)-1,kprb2(i)  ) )
! 
!         wprb = 0.5*(wo(iprb2(i)  ,jprb2(i)  ,kprb2(i)  )
!      &                     + wo(iprb2(i)  ,jprb2(i)  ,kprb2(i)-1) )

        pprb = p(iprb2(i)  ,jprb2(i)  ,kprb2(i))
         write(750+ii,'(5(4x,e15.8),i6,i6,i6)') time,uprb,vprb,wprb,pprb,iprb2(i),jprb2(i),kprb2(i)+myrank*(nz-2)
!         write(*,*) "WRITING PROBE DATA"
         close(750+ii)
!         write(*,*) time,uprb,vprb,wprb,pprb,iprb2(i),jprb2(i),kprb2(i)
      enddo

      return

      end

C------------------------------------------------------------------------


c---- subroutine rec_timeprobes ---------------N. Beratlis-03 May 2010---
C
C     PURPOSE: Record velocity and pressure at time probes
C
C------------------------------------------------------------------------
      subroutine rec_timeprobes(uo,vo,wo,p,nx,ny,nz,iprb,jprb,kprb,uprb,vprb,wprb,pprb,nprb,nprbmax)
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output array
      integer nprb,nx,ny,nz,nprbmax
      integer iprb(nprb),jprb(nprb),kprb(nprb)
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
      real    uprb(nprb),vprb(nprb),wprb(nprb),pprb(nprb)
c
c.... Local arrays
      integer i
      character*120 filename
c N.B. NPRBMAX is the number of probes inside the domain
c the centered values of velocity and pressure are evaluated
      do i=1,nprbmax

c        write(6,*) i,iprb(i),jprb(i),kprb(i)
        uprb(i) = 0.5*(uo(iprb(i)  ,jprb(i)  ,kprb(i)  )
     &                     + uo(iprb(i)-1,jprb(i)  ,kprb(i)  ) )

        vprb(i) = 0.5*(vo(iprb(i)  ,jprb(i)  ,kprb(i)  )
     &                     + vo(iprb(i)  ,jprb(i)-1,kprb(i)  ) )

        wprb(i) = 0.5*(wo(iprb(i)  ,jprb(i)  ,kprb(i)  )
     &                     + wo(iprb(i)  ,jprb(i)  ,kprb(i)-1) )

        pprb(i) = p(iprb(i)  ,jprb(i)  ,kprb(i))

      enddo

      return

      end
C------------------------------------------------------------------------


c---- subroutine write_probesheader -----------N. Beratlis-04 May 2010---
C
C     PURPOSE: Write header information to file containing probes' signal
C
C------------------------------------------------------------------------
      subroutine write_probesheader(filename,xc,yc,zc,nx,ny,nz,iprb,jprb,kprb,nprb,nprbmax,iloc,nvar,ti)
c      
      include 'common.h'
      include 'mpif.h'
c
      integer nx,ny,nz,nprb,nprbmax,iloc,nvar
      real    ti
      integer iprb(nprb),jprb(nprb),kprb(nprb)
      real    xc(nx),yc(ny),zc(nz)
      character*(*) filename
c
c.... Local array
      integer i,fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)


c      if(nprbmax>0 .AND. .false.) then
        call mpi_type_contiguous(1,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)      

        filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
        call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 8
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

        offset = 0
        call mpi_file_write_at(fh,offset,real(nprb),1,newtype,status,ierr)
        offset = 1
        call mpi_file_write_at(fh,offset,real(nvar),1,newtype,status,ierr)
        offset = 2
        call mpi_file_write_at(fh,offset,tstep,1,newtype,status,ierr)
        offset = 3
        call mpi_file_write_at(fh,offset,ti,1,newtype,status,ierr)

        offset = 4+(iloc-1)*6
        do i=1,nprbmax
          call mpi_file_write_at(fh,offset,real(iprb(i)),1,newtype,status,ierr)
          offset = offset+1
          call mpi_file_write_at(fh,offset,real(jprb(i)),1,newtype,status,ierr)
          offset = offset+1
          call mpi_file_write_at(fh,offset,real(kprb(i)+myrank*(nz-2)),1,newtype,status,ierr)
          offset = offset+1
          call mpi_file_write_at(fh,offset,xc(iprb(i)),1,newtype,status,ierr)
          offset = offset+1
          call mpi_file_write_at(fh,offset,yc(jprb(i)),1,newtype,status,ierr)
          offset = offset+1
          call mpi_file_write_at(fh,offset,zc(kprb(i)),1,newtype,status,ierr)
          offset = offset+1
        enddo
        call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------
c


c---- subroutine write_yprobesheader ----------N. Beratlis-04 May 2010---
C
C     PURPOSE: Write header information to file containing probes' signal
C
C------------------------------------------------------------------------
      subroutine write_yprobesheader(filename,iprb,kprb,nprb,nprbmax,iloc,nvar,n,nz)
c
      include 'common.h'
      include 'mpif.h'
c
      integer n,nz,nprb,nprbmax,iloc,nvar
      integer iprb(nprb),kprb(nprb)
      character*(*) filename
c
c.... Local array
      integer i,fh,nsample
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      nsample = 0.0

      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)      

      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 8
      call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

      offset = 0
      call mpi_file_write_at(fh,offset,real(nprb),1,newtype,status,ierr)
      offset = 1
      call mpi_file_write_at(fh,offset,real(nvar),1,newtype,status,ierr)
      offset = 2
      call mpi_file_write_at(fh,offset,real(n),1,newtype,status,ierr)

      offset = 3+(iloc-1)*2
      do i=1,nprbmax
        call mpi_file_write_at(fh,offset,real(iprb(i)),1,newtype,status,ierr)
        offset = offset+1
        call mpi_file_write_at(fh,offset,real(kprb(i)+myrank*(nz-2)),1,newtype,status,ierr)
        offset = offset+1
      enddo

      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------
c

c---- subroutine io_timeprobes ----------------N. Beratlis-03 May 2010---
C
C     PURPOSE: Record velocity and pressure at time probes
C
C------------------------------------------------------------------------
      subroutine io_timeprobes(filename,uprb,vprb,wprb,pprb,nprb,nser,iloc,it,nprbmax)
c
      include 'common.h'
c      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer nprb,nser,iloc,it,flag,nprbmax,nsample
      real    uprb(nprb,nser),vprb(nprb,nser),wprb(nprb,nser),pprb(nprb,nser) 
      character*(*) filename
c
c.... Local arrays
      integer i,j,fh,filemode,nvar
      real    var(4)
      real, dimension(:,:), allocatable :: iovar
      integer newtype,newtype2
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      nvar = 4

      if(nprbmax>0) then

c A new type is generated
        call mpi_type_contiguous(nvar*nprbmax,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
c disp is dependent on the number of probes and on the number of samples
        disp = 5*8+nprb*6*8+(it)*nprb*4*8
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)
      
        allocate(iovar(nvar,nprbmax))

        do i=1,nser
          do j=1,nprbmax
            iovar(1,j) = uprb(j,i)
            iovar(2,j) = vprb(j,i)
            iovar(3,j) = wprb(j,i)
            iovar(4,j) = pprb(j,i)
          enddo

c the offset is dependent on the number of probes NPRB, the number of
c variables NVAR and the index ILOC of the first probe for the specific
c processor
          offset = (i-1)*nprb*nvar+(iloc-1)*nvar
          call mpi_file_write_at(fh,offset,iovar,1,newtype,status,ierr)
          offset = offset + nprbmax*nvar
        enddo
        
        deallocate(iovar)

c the variable NEWTYPE is deallocated
        call mpi_type_free(newtype,ierr)
        call mpi_file_close(fh,ierr)

      endif

c the number of samples is increased considering the step ITRES
c (number of time instants for each plot)
      it = it+nser

c the processor 0 writes IT, the number of samples
      if(myrank.eq.0) call write_tprbsample(filename,it)

      return

      end
C------------------------------------------------------------------------

c---- subroutine write_tprbsample -------------N. Beratlis-03 May 2010---
C
C     PURPOSE: Record velocity and pressure at time probes
C
C------------------------------------------------------------------------
      subroutine write_tprbsample(filename,it)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer it
      character*(*) filename
c
c.... Local arrays
      integer i,j,fh,filemode,nvar
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

      call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)
      offset = 0

      call mpi_file_write_at(fh,offset,real(it),1,newtype,status,ierr)

      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------


c---- subroutine cp_timeprobes ----------------N. Beratlis-03 May 2010---
C
C     PURPOSE: Record velocity and pressure at time probes
C
C------------------------------------------------------------------------
      subroutine cp_timeprobes(filename1,filename2,nprb,append)
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer append
      character*(*) filename1,filename2
c
c.... Local arrays
      integer i,j,k,sample,nprb
      real    var(4),tmp
      integer fh,fh1,fh2
      integer filemode
      integer newtype,newtype2
      INTEGER(KIND=MPI_OFFSET_KIND) disp,disp1,disp2,offset,offset1,offset2
      integer status(mpi_status_size)

      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)      


c      write(6,*) 'inside, iapnd=',append

      sample = 0
      if(append==1) then
        filemode = MPI_MODE_RDONLY
        call mpi_file_open(mpi_comm_self,trim(filename2),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     &       ,'native',mpi_info_null,ierr)
        offset = 0
        call mpi_file_read_at(fh,offset,tmp,1,newtype,status,ierr)
        sample = int(tmp)
        call mpi_file_close(fh,ierr)
      endif
c
c.... Open filename1 to read
      filemode = MPI_MODE_RDONLY
      call mpi_file_open(mpi_comm_self,trim(filename1),filemode,mpi_info_null,fh1,ierr)
      disp1 = 0
      call mpi_file_set_view(fh1,disp1,mpi_double_precision,newtype
     &       ,'native',mpi_info_null,ierr)

c
c.... Open filename2 to write
      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      call mpi_file_open(mpi_comm_self,trim(filename2),filemode,mpi_info_null,fh2,ierr)
      disp2 = 3*8+nprb*6*8+sample*nprb*4*8
      call mpi_file_set_view(fh2,disp2,mpi_double_precision,newtype
     &       ,'native',mpi_info_null,ierr)

      offset = 0
      do i=1,itres
      do j=1,nprb
      do k=1,4
        call mpi_file_read_at(fh1,offset,tmp,1,newtype,status,ierr)
        call mpi_file_write_at(fh2,offset,tmp,1,newtype,status,ierr)
        offset = offset+1
      enddo
      enddo
      enddo
      call mpi_file_close(fh1,ierr)
      call mpi_file_close(fh2,ierr)

      sample = sample+itres

      call mpi_file_open(mpi_comm_self,trim(filename2),filemode,mpi_info_null,fh2,ierr)
      disp2 = 0
      call mpi_file_set_view(fh2,disp2,mpi_double_precision,newtype
     &     ,'native',mpi_info_null,ierr)
      offset = 0
      call mpi_file_write_at(fh2,offset,real(sample),1,newtype,status,ierr)
      call mpi_file_close(fh2,ierr)


      return

      end
C------------------------------------------------------------------------



C---- subroutine read_no_tprobes---------------N. Beratlis-10 May 2010---
C
C     PURPOSE: Read number of probes for azimuthal spectra.
C
C------------------------------------------------------------------------
      subroutine read_no_tprobes(n)
c
      implicit none
c
c.... Input/Output Arrays
      integer n
c
      open(unit=10,file='timeprobes.input',form='formatted')
      read(10,*)
      read(10,*) n
      close(10)

      return

      end
C------------------------------------------------------------------------


C---- subroutine read_no_yprobes---------------N. Beratlis-10 May 2010---
C
C     PURPOSE: Read number of probes for azimuthal spectra.
C
C------------------------------------------------------------------------
      subroutine read_no_yprobes(n)
c
      implicit none
c
c.... Input/Output Arrays
      integer n
c
      open(unit=10,file='yspctrprobes.input',form='formatted')
      read(10,*)
      read(10,*) n
      close(10)

      return

      end
C------------------------------------------------------------------------


C---- subroutine set_yprobes ------------------N. Beratlis-10 May 2010---
C
C     PURPOSE: Setup probes for azimuthal spectra.
C
C------------------------------------------------------------------------
      subroutine setup_tprobes(filename,iprb,jprb,kprb,nprb,xc,yc,zc,nx,ny,nz
     &     ,nprbmax,prbindx,nvar,nsample,ti)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer n,nprb,nx,ny,nz,nprbmax,prbindx,nvar,nsample
      real    ti
      integer iprb(nprb),jprb(nprb),kprb(nprb)
      real    xc(nx),yc(nz),zc(nz)
      character*(*) filename
c
c.... Local arrays
      integer i,ii,jj
      real    tmp
      integer fh,filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)
      real    xprb(nprb),yprb(nprb),zprb(nprb)

      open(unit=10,file='timeprobes.input',form='formatted')
      read(10,*)
      read(10,*) nprb
      do i=1,nprb
        read(10,*) xprb(i),yprb(i),zprb(i)
      enddo
      close(10)

      call sortprobes(xprb,yprb,zprb,nprb)

      jj=0
      prbindx=0
      do ii=1,nprb
!        if(zprb(ii)>zc(1) .AND. zprb(ii)<=zc(nz-1)) then           
        if(zprb(ii)>zc(2) .AND. zprb(ii)<=zc(nz)) then           
          jj = jj+1
          call closest(zc,nz,zprb(ii),kprb(jj))
          call closest(yc,ny,yprb(ii),jprb(jj))
          call closest(xc,nx,xprb(ii),iprb(jj))
          if(jprb(jj)<2) jprb(jj)=ny-jprb(jj)
          if(iprb(jj)==1 .AND. icyl==1) then
            iprb(jj)=2
            jprb(jj)=jsym(jprb(jj))
          endif
          if(jj==1) prbindx=ii
        endif
      enddo
      nprbmax = jj


      if(itmprb==1) call write_probesheader(trim(filename),xc,yc,zc,nx,ny,nz
     &     ,iprb,jprb,kprb,nprb,nprbmax,prbindx,nvar,ti)

      nsample = 0
      if(itmprb>1) then
        call mpi_type_contiguous(1,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)
        filemode = MPI_MODE_RDONLY
        call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     &       ,'native',mpi_info_null,ierr)
        offset = 0
        call mpi_file_read_at(fh,offset,tmp,1,newtype,status,ierr)
        call mpi_file_close(fh,ierr)
        nsample = int(tmp)
      endif

      return

      end
C------------------------------------------------------------------------


C---- subroutine set_yprobes ------------------N. Beratlis-10 May 2010---
C
C     PURPOSE: Setup probes for azimuthal spectra.
C
C------------------------------------------------------------------------
      subroutine setup_yprobes(filename,iprb,kprb,nprb,xc,zc,nx,nz
     &     ,nprbmax,prbindx,itsamp,nvar,nsample,n)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer n,nprb,itsamp,nx,nz,nprbmax,prbindx,nvar,nsample
      integer iprb(nprb),kprb(nprb)
      real    xc(nx),zc(nz)
      character*(*) filename
c
c.... Local arrays
      integer i,ii,jj
      real    tmp
      integer fh,filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)
      real    xprb(nprb),zprb(nprb)

      open(unit=10,file='yspctrprobes.input',form='formatted')
      read(10,*)
      read(10,*) nprb
      read(10,*)
      read(10,*) itsamp
      do i=1,nprb
        read(10,*) xprb(i),zprb(i)
      enddo
      close(10)

      call sortprobes(xprb,xprb,zprb,nprb)

      jj=0
      prbindx=0
      do ii=1,nprb
        if(zprb(ii)>zc(1) .AND. zprb(ii)<=zc(nz-1)) then
          jj = jj+1
          call closest(zc,nz,zprb(ii),kprb(jj))
          call closest(xc,nx,xprb(ii),iprb(jj))
          if(jj==1) prbindx=ii
        endif
      enddo
      nprbmax = jj

      if(yspctr==1) call write_yprobesheader(trim(filename),iprb,kprb,nprb,nprbmax,prbindx,nvar,n,nz)

      nsample = 0
      if(yspctr>1) then
        call mpi_type_contiguous(1,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)
        filemode = MPI_MODE_RDONLY
        call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     &       ,'native',mpi_info_null,ierr)
        offset = 0
        call mpi_file_read_at_all(fh,offset,tmp,1,newtype,status,ierr)
        call mpi_file_close(fh,ierr)
        nsample = int(tmp)
      endif

      return

      end
C------------------------------------------------------------------------


C---- subroutine correlation ------------------N. Beratlis-11 May 2010---
C
C     PURPOSE: Compute the correlation of 2 signals
C
C------------------------------------------------------------------------
      subroutine correlation(u1,u2,corr,n)
C
      implicit none
c
c.... Input/Output arrays
      integer n
      real    u1(n),u2(n),corr(n)
c
c.... Local arrays
      integer i,j,j1,j2

      corr = 0.0
      
      do i=1,n
      do j=1,n
        j1 = j
        j2 = j+i-1
        if(j2>n) j2=j2-n
        corr(i) = corr(i) + u1(j1)*u2(j2)
      enddo
      enddo

      return

      end
C------------------------------------------------------------------------

C---- subroutine compute_correlations ------- N. Beratlis-11 May 2010 ---
C
C     PURPOSE: Compute correlations of various terms at selected probe 
C     locations.
C
C------------------------------------------------------------------------
      subroutine compute_correlations(uo,vo,wo,p,nx,ny,nz
     &     ,iprb,jprb,kprb,nprb,nprbmax,ycorprb,ncorrvar)
C
      implicit none
C
c.... Input/Output Arrays
      integer nx,ny,nz,nprb,nprbmax,ncorrvar,nsample
      integer iprb(nprb),jprb(nprb),kprb(nprb)
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
      real    ycorprb(ny-2,nprb,4+ncorrvar)
c
c.... Local array
      integer i,j,k,ip,j1,j2
      real    ns,invs
      real    corr(ny-2)

      j1 = 2
      j2 = ny-1
      do ip=1,nprbmax
        i = iprb(ip)
        k = kprb(ip)

        !Store prime variables
        ycorprb(:,ip,1) = (ycorprb(:,ip,1) + Uo(i,j1:j2,k))
        ycorprb(:,ip,2) = (ycorprb(:,ip,2) + Vo(i,j1:j2,k))
        ycorprb(:,ip,3) = (ycorprb(:,ip,3) + Wo(i,j1:j2,k))
        ycorprb(:,ip,4) = (ycorprb(:,ip,4) + P (i,j1:j2,k))

        call correlation(uo(i,j1:j2,k),uo(i,j1:j2,k),corr,ny-2)
        ycorprb(:,ip,5) = (ycorprb(:,ip,5) + corr)
        call correlation(vo(i,j1:j2,k),vo(i,j1:j2,k),corr,ny-2)
        ycorprb(:,ip,6) = (ycorprb(:,ip,6) + corr)
        call correlation(wo(i,j1:j2,k),wo(i,j1:j2,k),corr,ny-2)
        ycorprb(:,ip,7) = (ycorprb(:,ip,7) + corr)
        call correlation(p (i,j1:j2,k),p (i,j1:j2,k),corr,ny-2)
        ycorprb(:,ip,8) = (ycorprb(:,ip,8) + corr)

      enddo

      return

      end
C------------------------------------------------------------------------



C---- subroutine write_correlations -----------N. Beratlis-03 May 2010---
C
C     PURPOSE: Record velocity and pressure correlations at selected probe
C     locations.
C
C------------------------------------------------------------------------
      subroutine write_correlations(filename,ycor,n,nprb,nvar,iloc,nprbmax,nsample)
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer n,nprb,nvar,nprbmax,iloc,nsample
      real    ycor(n,nprb,nvar)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,fh,filemode
      integer newtype,newtype2
      real, dimension(:,:,:), allocatable :: iovar
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)


      if(nprbmax>0) then

        call mpi_type_contiguous(n*nvar*nprbmax,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)


        allocate(iovar(n,nvar,nprbmax))

c the offset is defined using
c NPRB: number of probes
c ILOC: index of the first probe inside the specific processor
c N: number of points along the direction Y
c NVAR: number of variables for each probe
        offset = 4+2*nprb+(iloc-1)*n*nvar
        do i=1,nprbmax
        do j=1,nvar
        do k=1,n
           iovar(k,j,i) = ycor(k,i,j)
        enddo
        enddo
        enddo

        call mpi_file_write_at(fh,offset,iovar,1,newtype,status,ierr)

        deallocate(iovar)
        
        call mpi_type_free(newtype,ierr)
        call mpi_file_close(fh,ierr)

      endif

      if(myrank.eq.0) then
c the processor of rank 0 writes also the number of samples

        call mpi_type_contiguous(1,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)
        filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

        offset = 0
        call mpi_file_write_at(fh,offset,real(nsample),1,newtype,status,ierr)

        call mpi_type_free(newtype,ierr)
        call mpi_file_close(fh,ierr)

      endif
      
      return

      end
C------------------------------------------------------------------------


C---- subroutine read_ycorrelations -----------N. Beratlis-03 May 2010---
C
C     PURPOSE: Record velocity and pressure correlations at selected probe
C     locations.
C
C------------------------------------------------------------------------
      subroutine read_ycorrelations(filename,ycor,n,nprb,nvar,nprbmax,iloc)
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer n,nprb,nvar,iloc,nprbmax
      real    ycor(n,nprb,nvar)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,sample
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = MPI_MODE_RDONLY

      call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     &     ,'native',mpi_info_null,ierr)

      if(nprbmax>0) then
        offset = 4+2*nprb+(iloc-1)*n*nvar
        do i=1,nprbmax
        do j=1,nvar
        do k=1,n
          call mpi_file_read_at(fh,offset,ycor(k,i,j),1,newtype,status,ierr)
          offset = offset+1
        enddo
        enddo
        enddo
      endif
        
      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------




C---- subroutine compute_esd ---------------- N. Beratlis-11 May 2010 ---
C
C     PURPOSE: Compute correlations of various terms at selected probe 
C     locations.
C
C------------------------------------------------------------------------
      subroutine compute_esd(uo,vo,wo,p,nx,ny,nz
     &     ,iprb,kprb,nprb,nprbmax,esd,nvar,ini)
C
      implicit none
C
c.... Input/Output Arrays
      integer nx,ny,nz,nprb,nprbmax,nvar,ini
      integer iprb(nprb),kprb(nprb)
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
      real    esd(ny/2,nprb,nvar)
c
c.... Local array
      integer i,j,k,ip,j1,j2
      real    fft(ny-2)
      real    power(ny/2)
c      real    wsave(ny)
      real, dimension(:), allocatable, save :: wsave

      if(ini==0) then
c FFT is initialized
        allocate(wsave(2*ny+11))
        call rffti(ny-2,wsave)
      endif

      j1 = 2
      j2 = ny-1
      do ip=1,nprbmax
        i = iprb(ip)
        k = kprb(ip)

        !Compute and store esd of prime variables
        fft(:) = Uo(i,j1:j2,k)
        call rfftf(ny-2,fft,wsave)
        call fftpowcoef(fft,power,ny-2)
        esd(:,ip,1) = esd(:,ip,1) + power(:)

        fft(:) = Vo(i,j1:j2,k)
        call rfftf(ny-2,fft,wsave)
        call fftpowcoef(fft,power,ny-2)
        esd(:,ip,2) = esd(:,ip,2) + power(:)

        fft(:) = Wo(i,j1:j2,k)
        call rfftf(ny-2,fft,wsave)
        call fftpowcoef(fft,power,ny-2)
        esd(:,ip,3) = esd(:,ip,3) + power(:)

        fft(:) = P(i,j1:j2,k)
        call rfftf(ny-2,fft,wsave)
        call fftpowcoef(fft,power,ny-2)
        esd(:,ip,4) = esd(:,ip,4) + power(:)


      enddo


      return

      end
C------------------------------------------------------------------------




C---- subroutine fftpowcoef ----------------- N. Beratlis-14 May 2010 ---
C
C     PURPOSE: Compute the power coefficients of a FFT.
C
C------------------------------------------------------------------------
      subroutine fftpowcoef(fft,pow,n)

      implicit none
c
c.... Input/Output Arrays
      integer n
      real    fft(n),pow(n/2+1)
c
c.... Local arrays
      integer i

      pow = 0.0

      pow(1) = fft(1)**2.
      do i=2,n/2
        pow(i) = (fft(2*i-2)**2. + fft(2*i-1)**2.)
      enddo
      pow(n/2+1) = fft(n)**2.

      return

      end
C------------------------------------------------------------------------



C---- subroutine read_yspectra ----------------N. Beratlis-03 May 2010---
C
C     PURPOSE: Record velocity and pressure correlations at selected probe
C     locations.
C
C------------------------------------------------------------------------
      subroutine read_yspectra(filename,ycor,n,nprb,nvar,nprbmax,iloc)
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer n,nprb,nvar,iloc,nprbmax
      real    ycor(n,nprb,nvar)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,sample
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = MPI_MODE_RDONLY

      call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     &     ,'native',mpi_info_null,ierr)

      if(nprbmax>0) then
        offset = 4+2*nprb+(iloc-1)*n*nvar
        do i=1,nprbmax
        do j=1,nvar
        do k=1,n
          call mpi_file_read_at(fh,offset,ycor(k,i,j),1,newtype,status,ierr)
          offset = offset+1
        enddo
        enddo
        enddo
      endif
        
      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------



C---- subroutine read_tmavgprm --------------- N. Beratlis-29 May 2010---
C
C     PURPOSE: Read input parameters for time average sampling
C
C------------------------------------------------------------------------
      subroutine read_tmavgprm(lim,nxt,nyt,nzt,ftmavg,xc,yc,zc,zcg,nx,ny,nz,nzg)
C
c iflag=0: area defined by coordinates
c iflag=/0: area defined by indices
c also the frequency ftmavg of the time-averaging is defined
c      implicit none
      include 'common.h'
C
C.... Input/Output arrays
      integer lim(8),nxt,nyt,nzt,nzgt,nx,ny,nz,nzg,ftmavg
      real    xc(nx),yc(ny),zc(nz),zcg(nzg)
C
C.... Local arrays
      integer i1,i2,j1,j2,k1,k2,iflag,k1g,k2g,kl,ku
      real    x1,x2,y1,y2,z1,z2
      
      open(unit=10,file='tmavg.input',form='formatted')
      read(10,*)
      read(10,*) ftmavg
      read(10,*)
      read(10,*) iflag
      if(iflag==0) then
        read(10,*)
        read(10,*) x1,x2
        read(10,*)
        read(10,*) y1,y2
        read(10,*)
        read(10,*) z1,z2
      else
        read(10,*) 
        read(10,*) i1,i2
        read(10,*) 
        read(10,*) j1,j2
        read(10,*)
        read(10,*) k1g,k2g
      endif
      close(10)

      if(iflag==0) then
        call closest(xc,nx,x1,i1)
        call closest(xc,nx,x2,i2)
        call closest(yc,ny,y1,j1)
        call closest(yc,ny,y2,j2)

        k2 = 0
        k1 = nz

        if(z1>=zc(2) .AND. z2<=zc(nz-1)) then
          call closest(zc,nz,z1,k1)
          call closest(zc,nz,z2,k2)
        elseif(z1>=zc(2) .AND. z2>zc(nz-1)) then
          call closest(zc,nz,z1,k1)
          k2 = nz-1
        elseif(z1<zc(2) .AND. z2<=zc(nz-1)) then
          k1 = 2
          call closest(zc,nz,z2,k2)
        elseif(z1<zc(2) .AND. z2>zc(nz-1)) then
          k1 = 2
          k2 = nz-1
        endif
      else
        k2 = 0
        k1 = nz
        kl = k1g-myrank*(nz-2)
        ku = k2g-myrank*(nz-2)
        if(kl>=2 .AND. ku<=nz-1) then
          k1 = kl
          k2 = ku
        elseif(kl>=2 .AND. kl<=nz-1 .AND. ku>nz-1) then
          k1 = kl
          k2 = nz-1
        elseif(kl<2 .AND. ku>=2 .AND. ku<=nz-1) then
          k1 = 2
          k2 = ku
        elseif(kl<2 .AND. ku>nz-1) then
          k1 = 2
          k2 = nz-1
        endif
      endif
      
      if(i1<1) i1=1
      if(i2>nx) i2=nx
      if(j1<1) j1=1
      if(j2>ny) j2=ny
      if(k1<2) k1=2
      if(k2>nz-1) k2=nz-1

      if(k1g==1 .AND. myrank==0) k1=1
      if(k2g==nzg .AND. myrank==mysize-1) k2=nz

      lim(1) = i1
      lim(2) = i2
      lim(3) = j1
      lim(4) = j2
      lim(5) = k1
      lim(6) = k2

      if(iflag==0) then
        call closest(zcg,nzg,z1,lim(7))
        call closest(zcg,nzg,z2,lim(8))
      else
        lim(7) = k1g
        lim(8) = k2g
      endif

      nxt = i2-i1+1
      nyt = j2-j1+1
      nzt = k2-k1+1

      if(nzt<0) nzt=0

      return

      end
C------------------------------------------------------------------------


C---- compute_primevar_tmavg ---------------- N. Beratlis-May 31 2010 ---
C
C     PURPOSE: Compute time average prime variables
C
C------------------------------------------------------------------------
      subroutine compute_primevar_tmavg(var,nxt,nyt,nzt,nvar,uo,vo,wo,p,nx,ny,nz,lim,n)
C
      implicit none
c
c.... Input/Output arrays
      integer nx,ny,nz,nxt,nyt,nzt,nvar,n(nvar),lim(8)
      real*4  var(nxt,nyt,nzt,nvar)
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k,ii,jj,kk,is,js,ks,ie,je,ke
      real    invs(nvar),s(nvar)

      do i=1,nvar
        s(i) = real(n(i))
        invs(i) = 1/real(n(i)+1)
      enddo

      is = 1
      ie = nxt
      js = 1
      je = nyt
      ks = 1
      ke = nzt
      if(lim(1)==1) is=2
      if(lim(2)==nx) ie=nx-1
      if(lim(3)==1) js=2
      if(lim(4)==ny) je=ny-1
      if(lim(5)==1) ks=2
      if(lim(6)==nz) ke=nz-1

      do kk=ks,ke
      do jj=js,je
      do ii=is,ie
        i = lim(1)+ii-1
        j = lim(3)+jj-1
        k = lim(5)+kk-1
c.... Center collocation
c        var(ii,jj,kk,1) = (var(ii,jj,kk,1)*s(1) + 0.5*(uo(i,j,k)+uo(i-1,j,k)))*invs(1)
c        var(ii,jj,kk,2) = (var(ii,jj,kk,2)*s(2) + 0.5*(vo(i,j,k)+vo(i,j-1,k)))*invs(2)
c        var(ii,jj,kk,3) = (var(ii,jj,kk,3)*s(3) + 0.5*(wo(i,j,k)+wo(i,j,k-1)))*invs(3)
c        var(ii,jj,kk,4) = (var(ii,jj,kk,4)*s(4) + p(i,j,k))*invs(4)
c
c.... Face collocation
        var(ii,jj,kk,1) = (var(ii,jj,kk,1)*s(1) + uo(i,j,k))*invs(1)
        var(ii,jj,kk,2) = (var(ii,jj,kk,2)*s(2) + vo(i,j,k))*invs(2)
        var(ii,jj,kk,3) = (var(ii,jj,kk,3)*s(3) + wo(i,j,k))*invs(3)
        var(ii,jj,kk,4) = (var(ii,jj,kk,4)*s(4) + p(i,j,k))*invs(4)
      enddo
      enddo
      enddo

c      n = n+1

      return

      end
C------------------------------------------------------------------------


C---- compute_Restresses_tmavg ---------------- N. Beratlis-May 31 2010 ---
C
C     PURPOSE: Compute the normal stresses
C
C------------------------------------------------------------------------
      subroutine compute_Restresses_tmavg(var,nxt,nyt,nzt,nvar,uo,vo,wo,p,nx,ny,nz,lim,n)
C
      implicit none
c
c.... Input/Output arrays
      integer nx,ny,nz,nxt,nyt,nzt,nvar,n(nvar),lim(8)
      real*4  var(nxt,nyt,nzt,nvar)
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k,i1,j1,k1
      real    invs(nvar),s(nvar)

      i1 = lim(1)
      j1 = lim(3)
      k1 = lim(5)

      do i=1,nvar
        s(i) = real(n(i))
        invs(i) = 1/real(n(i)+1)
      enddo

      do k=1,nzt
      do j=1,nyt
      do i=1,nxt
        i1 = lim(1)+i-1
        j1 = lim(3)+j-1
        k1 = lim(5)+k-1
c.... Center collocation
c
c.... Face collocation
        var(i,j,k,1) =(var(i,j,k,1)*s(1) + uo(i1,j1,k1)*uo(i1,j1,k1))*invs(1)
        var(i,j,k,2) =(var(i,j,k,2)*s(2) + vo(i1,j1,k1)*vo(i1,j1,k1))*invs(2)
        var(i,j,k,3) =(var(i,j,k,3)*s(3) + wo(i1,j1,k1)*wo(i1,j1,k1))*invs(3)
        var(i,j,k,4) =(var(i,j,k,4)*s(4) + p(i1,j1,k1)*p(i1,j1,k1))*invs(4)
      enddo
      enddo
      enddo

      return

      end
C------------------------------------------------------------------------


C---- compute_crossvel_tmavg ---------------- N. Beratlis-May 31 2010 ---
C----------------------------------------------A. Posa - Feb 2014 -------
C
C     PURPOSE: Compute the shear stresses
C
C------------------------------------------------------------------------
      subroutine compute_crossvel_tmavg(var,nxt,nyt,nzt,nvar,uo,vo,wo,tv,nx,ny,nz,lim,n)
C
c      implicit none
      include 'common.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,nxt,nyt,nzt,nvar,n(nvar),lim(8)
      real*4  var(nxt,nyt,nzt,nvar)
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
      real    tv(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k,i1,j1,k1,k2
      real    uc,vc,wc
      real    invs(nvar),s(nvar)

      i1 = lim(1)
      j1 = lim(3)
      k1 = lim(5)

      do i=1,nvar
        s(i) = real(n(i))
        invs(i) = 1/real(n(i)+1)
      enddo

      do k=1,nzt
      k1 = lim(5)+k-1
      k2 = max(1,k1)
      do j=1,nyt
      j1 = lim(3)+j-1
      do i=1,nxt
      i1 = lim(1)+i-1
c.... Center collocation
        if(i1==1) then
          if(icyl==1) then
            uc = 0.5*(uo(1,j1,k1)-uo(2,jsym(j1),k1))
          else
            uc =-0.5*(uo(2,j1,k1)+uo(1,j1,k1))
          endif
        else
          uc = 0.5*(uo(i1,j1,k1)+uo(i1-1,j1,k1))
        endif

        if(j1==1) then
          vc = 0.5*(vo(i1,1,k1) + vo(i1,ny-2,k1))
        else
          vc = 0.5*(vo(i1,j1,k1)+vo(i1,j1-1,k1))
        endif

        wc= 0.5*(wo(i1,j1,k1)+wo(i1,j1,k2))
c
        var(i,j,k,1) =(var(i,j,k,1)*s(1) + uc*wc)*invs(1)
        var(i,j,k,2) =(var(i,j,k,2)*s(2) + uc*vc)*invs(2)
        var(i,j,k,3) =(var(i,j,k,3)*s(3) + vc*wc)*invs(3)
c        var(i,j,k,4) =(var(i,j,k,4)*s(4) + uo(i1,j1,k1)*wo(i1,j1,k1))*invs(1)
      enddo
      enddo
      enddo

      if(isgs.gt.0) then

        do k=1,nzt
        k1 = lim(5)+k-1
        do j=1,nyt
        j1 = lim(3)+j-1
        do i=1,nxt
        i1 = lim(1)+i-1
        var(i,j,k,4) =(var(i,j,k,4)*s(4) + tv(i1,j1,k1)/ru1)*invs(4)
        enddo
        enddo
        enddo

      endif

      return

      end
C------------------------------------------------------------------------


C---- write_tmavg_header -------------------- N. Beratlis-31 May 2010 ---
c
c     PURPOSE: Write time average variables to file
c
C------------------------------------------------------------------------
      subroutine write_tmavg_header(filename,nx,ny,nz,nvar,lim)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer nx,ny,nz,nvar,lim(8)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k
      real    n
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

      call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     &     ,'native',mpi_info_null,ierr)
      offset = 0 
      n = 0.0
      call mpi_file_write_at(fh,offset,n,1,newtype,status,ierr)
      offset = 1
      n = nvar
      call mpi_file_write_at(fh,offset,n,1,newtype,status,ierr)
      offset = 2
      n = lim(2)-lim(1)+1
      call mpi_file_write_at(fh,offset,n,1,newtype,status,ierr)
      offset = 3
      n = lim(4)-lim(3)+1
      call mpi_file_write_at(fh,offset,n,1,newtype,status,ierr)
      offset = 4
      n = lim(8)-lim(7)+1
      call mpi_file_write_at(fh,offset,n,1,newtype,status,ierr)

      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------


C---- write_tmavg_headers ------------------- N. Beratlis-31 May 2010 ---
c
c     PURPOSE: Write time average variables to file
c
C------------------------------------------------------------------------
      subroutine write_tmavg_headers(hspdir,nx,ny,nz,nvar,lim)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer nx,ny,nz,nvar,lim(8)
      character*(*) hspdir
      character*150 filename
c
c.... Local arrays
      integer i,j,k,iv
      real    n
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

      do iv=1,nvar

        if(iv==1) then
          filename = trim(hspdir)//'u_tmavg.bin.sp'
        elseif(iv==2) then
          filename = trim(hspdir)//'v_tmavg.bin.sp'
        elseif(iv==3) then
          filename = trim(hspdir)//'w_tmavg.bin.sp'
        elseif(iv==4) then
          filename = trim(hspdir)//'p_tmavg.bin.sp'
        elseif(iv==5) then
          filename = trim(hspdir)//'uu_tmavg.bin.sp'
        elseif(iv==6) then
          filename = trim(hspdir)//'vv_tmavg.bin.sp'
        elseif(iv==7) then
          filename = trim(hspdir)//'ww_tmavg.bin.sp'
        elseif(iv==8) then
          filename = trim(hspdir)//'pp_tmavg.bin.sp'
        elseif(iv==9) then
          filename = trim(hspdir)//'uw_tmavg.bin.sp'
        elseif(iv==10) then
          filename = trim(hspdir)//'uv_tmavg.bin.sp'
        elseif(iv==11) then
          filename = trim(hspdir)//'vw_tmavg.bin.sp'
        elseif(iv==12) then
          filename = trim(hspdir)//'uw1_tmavg.bin.sp'
        endif
        
        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     $       ,'native',mpi_info_null,ierr)
        offset = 0 
        n = 0.0
        call mpi_file_write_at(fh,offset,n,1,newtype,status,ierr)
        offset = 1
        n = nvar
        call mpi_file_write_at(fh,offset,n,1,newtype,status,ierr)
        offset = 2
        n = lim(2)-lim(1)+1
        call mpi_file_write_at(fh,offset,n,1,newtype,status,ierr)
        offset = 3
        n = lim(4)-lim(3)+1
        call mpi_file_write_at(fh,offset,n,1,newtype,status,ierr)
        offset = 4
        n = lim(8)-lim(7)+1
        call mpi_file_write_at(fh,offset,n,1,newtype,status,ierr)

        call mpi_file_close(fh,ierr)

      enddo

      call mpi_type_free(newtype,ierr)

      return

      end
C------------------------------------------------------------------------




C---- write_primevars_tmavg ------------------ N. Beratlis-31 May 2010 ---
c
c     PURPOSE: Write time average variables to file
c
C------------------------------------------------------------------------
      subroutine write_primevars_tmavg(hspdir,var,nxt,nyt,nzt,nvar,nsmpl,lim,nz)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer nxt,nyt,nzt,nz,nvar,nsmpl(nvar),lim(8)
      real*4  var(nxt,nyt,nzt,nvar)
      character*(*) hspdir
      character*150 filename
c
c.... Local arrays
      integer i,j,k,iv
      integer fh
      integer filemode
      integer newtype
c      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      INTEGER(KIND=8) disp,offset
      integer status(mpi_status_size)

      do iv=1,nvar

        call mpi_type_contiguous(1,MPI_REAL,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

        if(iv==1) then
          filename = trim(hspdir)//'u_tmavg.bin.sp'
        elseif(iv==2) then
          filename = trim(hspdir)//'v_tmavg.bin.sp'
        elseif(iv==3) then
          filename = trim(hspdir)//'w_tmavg.bin.sp'
        elseif(iv==4) then
          filename = trim(hspdir)//'p_tmavg.bin.sp'
        elseif(iv==5) then
          filename = trim(hspdir)//'uu_tmavg.bin.sp'
        elseif(iv==6) then
          filename = trim(hspdir)//'vv_tmavg.bin.sp'
        elseif(iv==7) then
          filename = trim(hspdir)//'ww_tmavg.bin.sp'
        elseif(iv==8) then
          filename = trim(hspdir)//'pp_tmavg.bin.sp'
        elseif(iv==9) then
          filename = trim(hspdir)//'uw_tmavg.bin.sp'
        elseif(iv==10) then
          filename = trim(hspdir)//'uv_tmavg.bin.sp'
        elseif(iv==11) then
          filename = trim(hspdir)//'vw_tmavg.bin.sp'
        elseif(iv==12) then
          filename = trim(hspdir)//'tv_tmavg.bin.sp'
        endif

        call mpi_file_open(mpi_comm_world,trim(filename),filemode,mpi_info_null,fh,ierr)

        disp = 0
        call mpi_file_set_view(fh,disp,MPI_REAL,newtype
     &     ,'native',mpi_info_null,ierr)

        offset = 0
        if(myrank==0) call mpi_file_write_at(fh,offset,real(nsmpl(iv),4),1,newtype,status,ierr)
        call mpi_type_free(newtype,ierr)
        call mpi_file_close(fh,ierr)

        if(nzt>0) then

          call mpi_type_contiguous(nxt*nyt*nzt,MPI_REAL,newtype,ierr)
          call mpi_type_commit(newtype,ierr)
          filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

          call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)

          disp = 0
          call mpi_file_set_view(fh,disp,MPI_REAL,newtype
     &         ,'native',mpi_info_null,ierr)

          offset = 5+nxt*nyt*((lim(5)+myrank*(nz-2))-lim(7))
          call mpi_file_write_at(fh,offset,var(:,:,:,iv),1,newtype,status,ierr)
        
          call mpi_file_close(fh,ierr)
          call mpi_type_free(newtype,ierr)

        endif

      enddo

      return

      end
C------------------------------------------------------------------------




C---- read_primevar_tmavg ------------------ N. Beratlis-31 May 2010 ---
c
c     PURPOSE: Write time average variables to file
c
C------------------------------------------------------------------------
      subroutine read_primevar_tmavg(filename,var,nxt,nyt,nzt,nvar,nsmpl,lim,nz)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer nxt,nyt,nzt,nvar,nsmpl,lim(8),nz
!!!!!!      real    var(nxt,nyt,nzt,nvar)
      real*4  var(nxt,nyt,nzt,nvar)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,iv
      real*4  temp
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      call mpi_type_contiguous(1,MPI_REAL,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = MPI_MODE_RDONLY

      call mpi_file_open(mpi_comm_world,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,MPI_REAL,newtype
     &     ,'native',mpi_info_null,ierr)

      offset = 0
      call mpi_file_read_at(fh,offset,temp,1,newtype,status,ierr)
      nsmpl = int(temp)

      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      if(nzt>0) then

        call mpi_type_contiguous(nxt*nyt*nzt,MPI_REAL,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        filemode = MPI_MODE_RDONLY

        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,MPI_REAL,newtype
     &       ,'native',mpi_info_null,ierr)

        do iv=1,nvar
          offset = 5+nxt*nyt*(lim(8)-lim(7)+1)*(iv-1)+nxt*nyt*((lim(5)+myrank*(nz-2))-lim(7))
          call mpi_file_read_at(fh,offset,var(:,:,:,iv),1,newtype,status,ierr)
          offset = offset + nxt*nyt*nzt
        enddo

        call mpi_type_free(newtype,ierr)
        call mpi_file_close(fh,ierr)

      endif

      return

      end
C------------------------------------------------------------------------


C---- write_primevar_tmavg ------------------ N. Beratlis-31 May 2010 ---
c
c     PURPOSE: Write time average variables to file
c
C------------------------------------------------------------------------
      subroutine write_primevar_tmavg(filename,var,nxt,nyt,nzt,nvar,nsmpl,lim,nz)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer nxt,nyt,nzt,nz,nvar,nsmpl,lim(8)
      real*4  var(nxt,nyt,nzt,nvar)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,iv
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

!!!!!!      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_contiguous(1,MPI_REAL,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

      call mpi_file_open(mpi_comm_world,trim(filename),filemode,mpi_info_null,fh,ierr)

      disp = 0
!!!!!!      call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
      call mpi_file_set_view(fh,disp,MPI_REAL,newtype
     &     ,'native',mpi_info_null,ierr)

      offset = 0
c the processor of rank 0 writes on file the number of samples
!!!!!!      if(myrank==0) call mpi_file_write_at(fh,offset,real(nsmpl),1,newtype,status,ierr)
      if(myrank==0) call mpi_file_write_at(fh,offset,real(nsmpl,4),1,newtype,status,ierr)
      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      if(nzt>0) then
!!!!!!        call mpi_type_contiguous(nxt*nyt*nzt,mtype,newtype,ierr)
        call mpi_type_contiguous(nxt*nyt*nzt,MPI_REAL,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)

        disp = 0
!!!!!!        call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
        call mpi_file_set_view(fh,disp,MPI_REAL,newtype
     &       ,'native',mpi_info_null,ierr)

c every processor writes on file
         do iv=1,nvar
           offset = 5+nxt*nyt*(lim(8)-lim(7)+1)*(iv-1)+nxt*nyt*((lim(5)+myrank*(nz-2))-lim(7))
           call mpi_file_write_at(fh,offset,var(:,:,:,iv),1,newtype,status,ierr)
           offset = offset + nxt*nyt*nzt
        enddo

        call mpi_type_free(newtype,ierr)
        call mpi_file_close(fh,ierr)

      endif

      return

      end
C------------------------------------------------------------------------


C---- read_primevars_tmavg ----------------- N. Beratlis-31 May 2010 ---
c
c     PURPOSE: Write time average variables to file
c
C------------------------------------------------------------------------
      subroutine read_primevars_tmavg(hspdir,var,nxt,nyt,nzt,nvar,nsmpl,lim,nz)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer nxt,nyt,nzt,nvar,nsmpl(nvar),lim(8),nz
      real*4  var(nxt,nyt,nzt,nvar)
      character*(*) hspdir
      character*150 filename
c
c.... Local arrays
      integer i,j,k,iv
      real*4  temp
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      do iv=1,nvar

        call mpi_type_contiguous(1,MPI_REAL,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        filemode = MPI_MODE_RDONLY

        if(iv==1) then
          filename = trim(hspdir)//'res.u_tmavg.bin.sp'
        elseif(iv==2) then
          filename = trim(hspdir)//'res.v_tmavg.bin.sp'
        elseif(iv==3) then
          filename = trim(hspdir)//'res.w_tmavg.bin.sp'
        elseif(iv==4) then
          filename = trim(hspdir)//'res.p_tmavg.bin.sp'
        elseif(iv==5) then
          filename = trim(hspdir)//'res.uu_tmavg.bin.sp'
        elseif(iv==6) then
          filename = trim(hspdir)//'res.vv_tmavg.bin.sp'
        elseif(iv==7) then
          filename = trim(hspdir)//'res.ww_tmavg.bin.sp'
        elseif(iv==8) then
          filename = trim(hspdir)//'res.pp_tmavg.bin.sp'
        elseif(iv==9) then
          filename = trim(hspdir)//'res.uw_tmavg.bin.sp'
        elseif(iv==10) then
          filename = trim(hspdir)//'res.uv_tmavg.bin.sp'
        elseif(iv==11) then
          filename = trim(hspdir)//'res.vw_tmavg.bin.sp'
        elseif(iv==12) then
          filename = trim(hspdir)//'res.tv_tmavg.bin.sp'
        endif

        call mpi_file_open(mpi_comm_world,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,MPI_REAL,newtype
     &       ,'native',mpi_info_null,ierr)

        offset = 0
        call mpi_file_read_at(fh,offset,temp,1,newtype,status,ierr)
        nsmpl(iv) = int(temp)

        call mpi_type_free(newtype,ierr)
        call mpi_file_close(fh,ierr)

        if(nzt>0) then

          call mpi_type_contiguous(nxt*nyt*nzt,MPI_REAL,newtype,ierr)
          call mpi_type_commit(newtype,ierr)

          filemode = MPI_MODE_RDONLY

          call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
          disp = 0
          call mpi_file_set_view(fh,disp,MPI_REAL,newtype
     &       ,'native',mpi_info_null,ierr)

          offset = 5+nxt*nyt*((lim(5)+myrank*(nz-2))-lim(7))
          call mpi_file_read_at(fh,offset,var(:,:,:,iv),1,newtype,status,ierr)

          call mpi_type_free(newtype,ierr)
          call mpi_file_close(fh,ierr)

        endif

      enddo

      return

      end
C------------------------------------------------------------------------


C---- subroutine region2index---------------- N. Beratlis-23 Jun 2010 ---
C
C     PURPOSE: Compute indices bounding region of interest
C
C------------------------------------------------------------------------
      subroutine region2indices(kmin,kmax,iind,xc,zc,nx,ny,nz,n,nprev,nzprev,j1,j2)
c
c      implicit none
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,kmin,kmax,n,nprev,nzprev,j1,j2
      integer iind(nz,2)      
      real    xc(nx),zc(nz)
c
c.... Local arrays
      integer i,j,k,jp,ip,imin,imax,n1,n2
      real    r,rmin,rmax
      integer mpi_comm_prev,mpi_group_prev
      integer status(mpi_status_size)


      rmin = 0.49
      rmax = 0.53

      j1 = 2
      j2 = ny-1

      n = 0
      n2 = 0
      nzprev = 0

      kmin = 100000
      kmax =-100000

      do k=2,nz-1
        imin = 100000
        imax =-100000

        do i=2,nx
          if(zc(k)>=-0.51 .AND. zc(k)<0.51) then
            r = sqrt(zc(k)**2. + xc(i)**2.)
            if(r<=rmax .AND. r>=rmin) then
              if(i<=imin) imin = i
              if(i>=imax) imax = i
              if(k<=kmin) kmin = k
              if(k>=kmax) kmax = k
              n2 = n2+1
            endif
c          elseif(zc(k)>=-0.1 .AND. zc(k)<rmin) then
c            r = sqrt(zc(k)**2. + xc(i)**2.)
c            if(r>=rmin .AND. xc(i)<0.8) then
c             if(i<=imin) imin = i
c             if(i>=imax) imax = i
c             if(k<=kmin) kmin = k
c             if(k>=kmax) kmax = k 
cc             n = n+1
c            endif
c          elseif(zc(k)>=rmin .AND. zc(k)<1.5) then 
c            imin = 2
c            call locate(xc,nx,0.8,imax) 
c            if(k<=kmin) kmin = k
c            if(k>=kmax) kmax = k 
cc            n = n+1
          endif

        enddo

        iind(k,1) = imin
        iind(k,2) = imax

c        if(imax>=imin) n2 = n2+imax-imin+1
c        write(6,*) k,zc(k),imin,imax,xc(imin),xc(imax),n,n2
      enddo

      write(6,*) 'n2=',n2
      
      kmin = kmin
c      write(6,*) 'myrank=',myrank,',kmin=',kmin,',kmax=',kmax

      n = 0
      do k=kmin,kmax
        n = n + iind(k,2)-iind(k,1)+1
      enddo

      nzprev = 0
      nprev = 0

      if(mysize>0) then

        if(myrank>0) then
c          do ip=1,mysize
c            if(ip==myleft) then
c            CALL MPI_RECV(jp,1,MPI_INTEGER,myleft,0,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(i,1,MPI_INTEGER,myleft,1,MPI_COMM_EDDY,STATUS,IERR)
              nprev = nprev+i
              CALL MPI_RECV(j,1,MPI_INTEGER,myleft,2,MPI_COMM_EDDY,STATUS,IERR)
              nzprev = nzprev+j
c              write(6,*) 'myrank=',myrank,'receiving',i,j,nprev,nzprev
c            endif
c          enddo
        endif

        if(myrank<mysize-1) then
c          CALL MPI_SEND(myrank,1,MPI_INTEGER,myright,0,MPI_COMM_EDDY,IERR)
          i = n
          CALL MPI_SEND(nprev+i,1,MPI_INTEGER,myright,1,MPI_COMM_EDDY,IERR)
          j = 0
          if(kmax>=kmin) j = kmax-kmin+1
          CALL MPI_SEND(nzprev+j,1,MPI_INTEGER,myright,2,MPI_COMM_EDDY,IERR)
c          write(6,*) 'myrank=',myrank,'sending',i,j,nprev+n,nzprev+j
        endif

      endif

      return

      end
C------------------------------------------------------------------------



C---- subroutine read_regtmavgprm --------------- N. Beratlis-29 May 2010---
C
C     PURPOSE: Read input parameters for time average sampling
C
C------------------------------------------------------------------------
      subroutine read_regtmavgprm(ftmavg)
C
      implicit none
C
C.... Input/Output arrays
      integer ftmavg
      
      open(unit=10,file='regtmavg.input',form='formatted')
      read(10,*)
      read(10,*) ftmavg
      close(10)

      return

      end
C------------------------------------------------------------------------




C---- write_tmavg_header -------------------- N. Beratlis-31 May 2010 ---
c
c     PURPOSE: Write time average variables to file
c
C------------------------------------------------------------------------
      subroutine write_regtmavg_header(filename,nreg,nvar)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer nreg,nvar
      character*(*) filename
c
c.... Local arrays
      integer i,j,k
      real    n
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

      call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     &     ,'native',mpi_info_null,ierr)
      offset = 0 
      n = 0.0
      call mpi_file_write_at(fh,offset,n,1,newtype,status,ierr)
      offset = 1
      n = real(nvar)
      call mpi_file_write_at(fh,offset,n,1,newtype,status,ierr)
      offset = 2
      n = real(nreg)
      call mpi_file_write_at(fh,offset,n,1,newtype,status,ierr)
        
      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------



C---- write_regindices ---------------------- N. Beratlis-31 May 2010 ---
c
c     PURPOSE: Write time average variables to file
c
C------------------------------------------------------------------------
      subroutine write_regindices(filename,kmin,kmax,ind,nz,nzprev,nzgreg,j1,j2)
c
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer nz,nzprev,nzgreg,kmin,kmax,j1,j2
      integer ind(nz,2)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,k1,k2,nk
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      call mpi_type_contiguous(1,mpi_integer,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

      call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)      
      disp = 0
      call mpi_file_set_view(fh,disp,mpi_integer,mpi_integer
     &     ,'native',mpi_info_null,ierr)

      nk = 0
      if(kmax>=kmin) nk=kmax-kmin+1
      CALL MPI_ALLREDUCE(nk,nzgreg,1,mpi_integer,MPI_SUM,MPI_COMM_EDDY,IERR)

      if(myrank.eq.0) then
        offset = 0 
        call mpi_file_write_at(fh,offset,j1,1,newtype,status,ierr)
        offset = offset+1
        call mpi_file_write_at(fh,offset,j2,1,newtype,status,ierr)
        offset = offset+1
        call mpi_file_write_at(fh,offset,nzgreg,1,newtype,status,ierr)
      endif
c      write(6,*) 'nzgreg=',nzgreg
c      call mpi_finalize(ierr)
c      stop

      offset = nzprev*3+3
c      write(6,*) 'myrank=',myrank,'offset=',offset,',nzprev=',nzprev,kmin,kmax
      do k=kmin,kmax
        call mpi_file_write_at(fh,offset,k+myrank*(nz-2),1,newtype,status,ierr)
        offset = offset+1
        call mpi_file_write_at(fh,offset,ind(k,1),1,newtype,status,ierr)
        offset = offset+1
        call mpi_file_write_at(fh,offset,ind(k,2),1,newtype,status,ierr)
        offset = offset+1
      enddo
        
      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------


C---- compute_primevar_tmavg ---------------- N. Beratlis-May 31 2010 ---
C
C     PURPOSE: Compute time average prime variables
C
C------------------------------------------------------------------------
      subroutine compute_regvar_tmavg(var,npts,nvar,ind,kmin,kmax,j1,j2,uo,vo,wo,p,nx,ny,nz,nsmpl)
C
      implicit none
c
c.... Input/Output arrays
      integer nx,ny,nz,nvar,npts,nsmpl,kmin,kmax,j1,j2
      integer ind(nz,2)
      real    var(npts,nvar)
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k,ipt
      real    invs,s

      s = real(nsmpl)
      invs = 1/real(nsmpl+1)

      ipt = 1
      do k=kmin,kmax
      do j=j1,j2
      do i=ind(k,1),ind(k,2)
c        var(ipt,1) =(var(ipt,1)*s + 0.5*(uo(i,j,k)+uo(i-1,j,k)))*invs
c        var(ipt,2) =(var(ipt,2)*s + 0.5*(vo(i,j,k)+vo(i,j-1,k)))*invs
c        var(ipt,3) =(var(ipt,3)*s + 0.5*(wo(i,j,k)+wo(i,j,k-1)))*invs
        var(ipt,1) =(var(ipt,1)*s + uo(i,j,k))*invs
        var(ipt,2) =(var(ipt,2)*s + vo(i,j,k))*invs
        var(ipt,3) =(var(ipt,3)*s + wo(i,j,k))*invs
        var(ipt,4) =(var(ipt,4)*s + p(i,j,k))*invs
        ipt = ipt+1
      enddo
      enddo
      enddo

c      write(6,*) var(:,3)

      nsmpl = nsmpl+1

      return

      end
C------------------------------------------------------------------------




C---- write_regvar_tmavg -------------------- N. Beratlis-31 May 2010 ---
c
c     PURPOSE: Write time average variables to file
c
C------------------------------------------------------------------------
      subroutine write_regvar_tmavg(filename,var,npts,nvar,ind,nz,kmin,kmax
     &     ,nsmpl,nprev,nptsg)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer npts,nvar,nsmpl,nprev,nz,kmin,kmax,nptsg
      integer ind(nz,2)
      real    var(npts,nvar)
      character*(*) filename
c
c.... Local arrays
      integer iv,ip
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

      call mpi_file_open(mpi_comm_world,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     &     ,'native',mpi_info_null,ierr)

      offset = 0
c the 0 processor writes the number of samples
      if(myrank==0) call mpi_file_write_at(fh,offset,real(nsmpl),1,newtype,status,ierr)

c      write(6,*) 'myrank=',myrank,kmin,kmax
      if(kmax-kmin+1>0) then
c        write(6,*) 'offset=',offset,myrank,nvar,npts,nprev
c        offset = 3+nprev*iv
c
c every processor writes the averaged values
c the offset is evaluated considering the total number of points NPTSG,
c the number of points at the left processors NPREV and the number of
c the current variable IV 
        do iv=1,nvar
          offset = 3+(iv-1)*nptsg+nprev
        do ip=1,npts
          call mpi_file_write_at(fh,offset,var(ip,iv),1,newtype,status,ierr)
          offset = offset+1
        enddo
        enddo

c        do ip=1,10
c          write(6,*) var(ip,1),var(ip,2),var(ip,3),var(ip,4)
c        enddo
      endif

      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------



C---- read_regvar_tmavg -------------------- N. Beratlis-31 May 2010 ---
c
c     PURPOSE: Read time average variables from file
c
C------------------------------------------------------------------------
      subroutine read_regvar_tmavg(filename,var,npts,nvar,ind,nz,kmin,kmax
     &     ,nsmpl,nprev,nptsg)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer npts,nvar,nsmpl,nprev,nz,kmin,kmax,nptsg
      integer ind(nz,2)
      real    var(npts,nvar)
      character*(*) filename
c
c.... Local arrays
      integer iv,ip
      real    tmp
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = MPI_MODE_RDONLY

      call mpi_file_open(mpi_comm_world,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     &     ,'native',mpi_info_null,ierr)

      offset = 0
      call mpi_file_read_at(fh,offset,tmp,1,newtype,status,ierr)
      nsmpl = int(tmp)

      if(kmax-kmin+1>0) then
        do iv=1,nvar
          offset = 3+(iv-1)*nptsg+nprev
        do ip=1,npts
          call mpi_file_read_at(fh,offset,var(ip,iv),1,newtype,status,ierr)
          offset = offset+1
        enddo
        enddo
      endif

      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------





C---- subroutine read_VPfieldprm --------------- N. Beratlis-29 May 2010---
C
C     PURPOSE: Read input parameters for time average sampling
C
C------------------------------------------------------------------------
      subroutine read_VPfieldprm(lim,nxt,nyt,nzt,xc,yc,zc,zcg,nx,ny,nz,nzg)
C
c      implicit none
      include 'common.h'
C
C.... Input/Output arrays
      integer lim(8),nxt,nyt,nzt,nzgt,nx,ny,nz,nzg,ftmavg
      real    xc(nx),yc(ny),zc(nz),zcg(nzg)
C
C.... Local arrays
      integer i1,i2,j1,j2,k1,k2,iflag,k1g,k2g,kl,ku
      real    x1,x2,y1,y2,z1,z2
      
      open(unit=10,file='VPfield.input',form='formatted')
      read(10,*)
      read(10,*) nVPfield   ! the index of the next print
      read(10,*)
      read(10,*) iflag
      if(iflag==0) then
        read(10,*)
        read(10,*) x1,x2
        read(10,*)
        read(10,*) y1,y2
        read(10,*)
        read(10,*) z1,z2
      else
        read(10,*) 
        read(10,*) i1,i2
        read(10,*) 
        read(10,*) j1,j2
        read(10,*)
        read(10,*) k1g,k2g
      endif
      close(10)

      if(iflag==0) then
        call closest(xc,nx,x1,i1)
        call closest(xc,nx,x2,i2)
        call closest(yc,ny,y1,j1)
        call closest(yc,ny,y2,j2)

        k2 = 0
        k1 = nz

        if(z1>=zc(2) .AND. z2<=zc(nz-1)) then
          call closest(zc,nz,z1,k1)
          call closest(zc,nz,z2,k2)
        elseif(z1>=zc(2) .AND. z2>zc(nz-1)) then
          call closest(zc,nz,z1,k1)
          k2 = nz-1
        elseif(z1<zc(2) .AND. z2<=zc(nz-1)) then
          k1 = 2
          call closest(zc,nz,z2,k2)
        elseif(z1<zc(2) .AND. z2>zc(nz-1)) then
          k1 = 2
          k2 = nz-1
        endif
      else
        k2 = 0
        k1 = nz
        kl = k1g-myrank*(nz-2)
        ku = k2g-myrank*(nz-2)
        if(kl>=2 .AND. ku<=nz-1) then
          k1 = kl
          k2 = ku
        elseif(kl>=2 .AND. kl<=nz-1 .AND. ku>nz-1) then
          k1 = kl
          k2 = nz-1
        elseif(kl<2 .AND. ku>=2 .AND. ku<=nz-1) then
          k1 = 2
          k2 = ku
        elseif(kl<2 .AND. ku>nz-1) then
          k1 = 2
          k2 = nz-1
        endif
      endif
      
      if(i1<1) i1=1
      if(i2>nx) i2=nx
      if(j1<1) j1=1
      if(j2>ny) j2=ny
      if(k1<2) k1=2
      if(k2>nz-1) k2=nz-1

      if(k1g==1 .AND. myrank==0) k1=1
      if(k2g==nzg .AND. myrank==mysize-1) k2=nz
      
      lim(1) = i1
      lim(2) = i2
      lim(3) = j1
      lim(4) = j2
      lim(5) = k1
      lim(6) = k2

      if(iflag==0) then
        call closest(zcg,nzg,z1,lim(7))
        call closest(zcg,nzg,z2,lim(8))
      else
        lim(7) = k1g
        lim(8) = k2g
      endif

      nxt = i2-i1+1
      nyt = j2-j1+1
      nzt = k2-k1+1

      return

      end
C------------------------------------------------------------------------


C---- subroutine write_grid_subdom ------- N. Beratlis - 27 Sep. 2010 ---
C
C     PURPOSE: Write grid for VPfield
C
C------------------------------------------------------------------------
      subroutine write_grid_subdom(filename,xu,yv,zw,xc,yc,zc,nx,ny,nz,lim)
c
c      implicit none
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,lim(8)
      real    xu(nx),yv(ny),zw(nz),xc(nx),yc(ny),zc(nz)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,nxl,nyl,nzl,i1,i2,j1,j2,k1,k2
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      nxl = lim(2)-lim(1)+1
      nyl = lim(4)-lim(3)+1
      nzl = lim(8)-lim(7)+1
      i1 = lim(1)
      i2 = lim(2)
      j1 = lim(3)
      j2 = lim(4)
      k1 = lim(5)
      k2 = lim(6)

      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)      
      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

      offset = 0
      if(myrank.eq.0) then
        offset = 0
        call mpi_file_write_at(fh,offset,real(nxl),1,newtype,status,ierr)
        offset = 1
        call mpi_file_write_at(fh,offset,real(nyl),1,newtype,status,ierr)
        offset = 2
        call mpi_file_write_at(fh,offset,real(nzl),1,newtype,status,ierr)
        offset = 3

        do i=i1,i2
          call mpi_file_write_at(fh,offset,xu(i),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do i=i1,i2
          call mpi_file_write_at(fh,offset,xc(i),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do i=i1,i2
          call mpi_file_write_at(fh,offset,au(i),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do i=i1,i2
          call mpi_file_write_at(fh,offset,av(i),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do i=i1,i2
          call mpi_file_write_at(fh,offset,ap(i),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,yv(j),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,yc(j),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,bu(j),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,bv(j),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,bp(j),1,newtype,status,ierr)
          offset = offset+1
        enddo

      endif

      if( (lim(6)-lim(5)+1)>0 ) then

        offset = 3 + nxl*5 + nyl*5 + (lim(5)+myrank*(nz-2)-lim(7))
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,zw(k),1,newtype,status,ierr)
          offset = offset+1
        enddo

        offset = 3 + nxl*5 + nyl*5 + (lim(5)+myrank*(nz-2)-lim(7)) + nzl
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,zc(k),1,newtype,status,ierr)
          offset = offset+1
        enddo

        offset = 3 + nxl*5 + nyl*5 + (lim(5)+myrank*(nz-2)-lim(7)) + 2*nzl
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,cu(k),1,newtype,status,ierr)
          offset = offset+1
        enddo

        offset = 3 + nxl*5 + nyl*5 + (lim(5)+myrank*(nz-2)-lim(7)) + 3*nzl
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,cw(k),1,newtype,status,ierr)
          offset = offset+1
        enddo

        offset = 3 + nxl*5 + nyl*5 + (lim(5)+myrank*(nz-2)-lim(7)) + 4*nzl
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,cp(k),1,newtype,status,ierr)
          offset = offset+1
        enddo

      endif

      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------


C---- subroutine write_grid_subdom_sp ---- N. Beratlis - 27 Sep. 2010 ---
C
C     PURPOSE: Write grid for VPfield
C
C------------------------------------------------------------------------
      subroutine write_grid_subdom_sp(filename,xu,yv,zw,xc,yc,zc,nx,ny,nz,lim)
c
c      implicit none
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,lim(8)
      real    xu(nx),yv(ny),zw(nz),xc(nx),yc(ny),zc(nz)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,nxl,nyl,nzl,i1,i2,j1,j2,k1,k2
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      nxl = lim(2)-lim(1)+1
      nyl = lim(4)-lim(3)+1
      nzl = lim(8)-lim(7)+1
      i1 = lim(1)
      i2 = lim(2)
      j1 = lim(3)
      j2 = lim(4)
      k1 = lim(5)
      k2 = lim(6)

      call mpi_type_contiguous(1,MPI_REAL,newtype,ierr)
      call mpi_type_commit(newtype,ierr)      
      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,MPI_REAL,newtype,'native',mpi_info_null,ierr)

      offset = 0
      if(myrank.eq.0) then
        offset = 0
        call mpi_file_write_at(fh,offset,real(nxl,4),1,newtype,status,ierr)
        offset = 1
        call mpi_file_write_at(fh,offset,real(nyl,4),1,newtype,status,ierr)
        offset = 2
        call mpi_file_write_at(fh,offset,real(nzl,4),1,newtype,status,ierr)
        offset = 3

        do i=i1,i2
          call mpi_file_write_at(fh,offset,real(xu(i),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do i=i1,i2
          call mpi_file_write_at(fh,offset,real(xc(i),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do i=i1,i2
          call mpi_file_write_at(fh,offset,real(au(i),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do i=i1,i2
          call mpi_file_write_at(fh,offset,real(av(i),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do i=i1,i2
          call mpi_file_write_at(fh,offset,real(ap(i),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,real(yv(j),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,real(yc(j),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,real(bu(j),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,real(bv(j),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,real(bp(j),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

      endif

      if( (lim(6)-lim(5)+1)>0 ) then

        offset = 3 + nxl*5 + nyl*5 + (lim(5)+myrank*(nz-2)-lim(7))
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,real(zw(k),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

        offset = 3 + nxl*5 + nyl*5 + (lim(5)+myrank*(nz-2)-lim(7)) + nzl
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,real(zc(k),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

        offset = 3 + nxl*5 + nyl*5 + (lim(5)+myrank*(nz-2)-lim(7)) + 2*nzl
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,real(cu(k),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

        offset = 3 + nxl*5 + nyl*5 + (lim(5)+myrank*(nz-2)-lim(7)) + 3*nzl
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,real(cw(k),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

        offset = 3 + nxl*5 + nyl*5 + (lim(5)+myrank*(nz-2)-lim(7)) + 4*nzl
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,real(cp(k),4),1,newtype,status,ierr)
          offset = offset+1
        enddo

      endif

      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------


C---- subroutine write_VPfield ----------- N. Beratlis - 28 Sep. 2010 ---
C
C     PURPOSE: Write VPfield file
C
C------------------------------------------------------------------------
      subroutine write_VPfield(hspdir,fileend,uo,vo,wo,p,nx,ny,nz,lim)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,lim(8)
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
      character*(*) fileend,hspdir
      character*300 filename
c      
c.... Local arrays
      integer k,i1,i2,j1,j2,nxl,nyl,nzl
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)
c
      nxl = lim(2)-lim(1)+1
      nyl = lim(4)-lim(3)+1
      nzl = lim(6)-lim(5)+1

      i1 = lim(1)
      i2 = lim(2)
      j1 = lim(3)
      j2 = lim(4)

      
      call mpi_type_contiguous(nxl*nyl,mpi_real,newtype,ierr)
      call mpi_type_commit(newtype,ierr)
      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

      filename = trim(hspdir)//'U'//trim(fileend)
      call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mpi_real,newtype,'native',mpi_info_null,ierr)

      if(nzl>0) then

        offset = (lim(5)+myrank*(nz-2)-lim(7))*nxl*nyl
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,real(uo(i1:i2,j1:j2,k),4),1,newtype,status,ierr)
          offset = offset + nxl*nyl
        enddo

      endif
      call mpi_file_close(fh,ierr)

      filename = trim(hspdir)//'V'//trim(fileend)
      call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mpi_real,newtype,'native',mpi_info_null,ierr)

      if(nzl>0) then

        offset = (lim(5)+myrank*(nz-2)-lim(7))*nxl*nyl
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,real(vo(i1:i2,j1:j2,k),4),1,newtype,status,ierr)
          offset = offset + nxl*nyl
        enddo

      endif
      call mpi_file_close(fh,ierr)

      filename = trim(hspdir)//'W'//trim(fileend)
      call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mpi_real,newtype,'native',mpi_info_null,ierr)

      if(nzl>0) then

        offset = (lim(5)+myrank*(nz-2)-lim(7))*nxl*nyl
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,real(wo(i1:i2,j1:j2,k),4),1,newtype,status,ierr)
          offset = offset + nxl*nyl
        enddo

      endif
      call mpi_file_close(fh,ierr)


      filename = trim(hspdir)//'P'//trim(fileend)
      call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mpi_real,newtype,'native',mpi_info_null,ierr)

      if(nzl>0) then

        offset = (lim(5)+myrank*(nz-2)-lim(7))*nxl*nyl
        do k=lim(5),lim(6)
          call mpi_file_write_at(fh,offset,real(p(i1:i2,j1:j2,k),4),1,newtype,status,ierr)
          offset = offset + nxl*nyl
        enddo

      endif
      call mpi_file_close(fh,ierr)


      call mpi_type_free(newtype,ierr)

      return

      end
C------------------------------------------------------------------------


C---- subroutine read_VPreg_prms---------------N. Beratlis-15 May 2011---
C
C     PURPOSE: Read parameters for VP region files
C
C------------------------------------------------------------------------
      subroutine read_VPreg_prms(freq,ifile)
c
      implicit none
c
      integer freq,ifile
c
      open(unit=10,file='VPreg.input',form='formatted')
      read(10,'(A)')
      read(10,*) freq
      read(10,'(A)')
      read(10,*) ifile
      close(10)

      return
c
      end
C------------------------------------------------------------------------


C---- subroutine VPregion2indices---------------- N. Beratlis-23 Jun 2010 ---
C
C     PURPOSE: Compute indices bounding region of interest
C
C------------------------------------------------------------------------
      subroutine VPregion2indices(kmin,kmax,iind,xc,zc,nx,ny,nz,n,nprev,nzprev,j1,j2)
c
c      implicit none
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,kmin,kmax,n,nprev,nzprev,j1,j2
      integer iind(nz,2)      
      real    xc(nx),zc(nz)
c
c.... Local arrays
      integer i,j,k,jp,ip,imin,imax,n1,n2
      real    r,rmin,rmax
      integer mpi_comm_prev,mpi_group_prev
      integer status(mpi_status_size)

      rmin = 0.49
      rmax = 0.53

      j1 = 2
      j2 = ny-1

      n = 0
      n2 = 0
      nzprev = 0

      kmin = 100000
      kmax =-100000

      do k=2,nz-1
        imin = 100000
        imax =-100000

        do i=2,nx
          if(zc(k)>=-0.53 .AND. zc(k)<0.53) then
            r = sqrt(zc(k)**2. + xc(i)**2.)
            if(r<=rmax .AND. r>=rmin) then
              if(i<=imin) imin = i
              if(i>=imax) imax = i
              if(k<=kmin) kmin = k
              if(k>=kmax) kmax = k
              n2 = n2+1
            endif
          endif

        enddo

        iind(k,1) = imin
        iind(k,2) = imax

c        if(imax>=imin) n2 = n2+imax-imin+1
c        write(6,*) k,zc(k),imin,imax,xc(imin),xc(imax),n,n2
      enddo

      write(6,*) 'n2=',n2
      
c      kmin = kmin
c      write(6,*) 'myrank=',myrank,',kmin=',kmin,',kmax=',kmax

      n = 0
      do k=kmin,kmax
        n = n + iind(k,2)-iind(k,1)+1
      enddo

      nzprev = 0
      nprev = 0

      if(mysize>0) then

        if(myrank>0) then
c          do ip=1,mysize
c            if(ip==myleft) then
c            CALL MPI_RECV(jp,1,MPI_INTEGER,myleft,0,MPI_COMM_EDDY,STATUS,IERR)
              CALL MPI_RECV(i,1,MPI_INTEGER,myleft,1,MPI_COMM_EDDY,STATUS,IERR)
              nprev = nprev+i
              CALL MPI_RECV(j,1,MPI_INTEGER,myleft,2,MPI_COMM_EDDY,STATUS,IERR)
              nzprev = nzprev+j
c              write(6,*) 'myrank=',myrank,'receiving',i,j,nprev,nzprev
c            endif
c          enddo
        endif

        if(myrank<mysize-1) then
c          CALL MPI_SEND(myrank,1,MPI_INTEGER,myright,0,MPI_COMM_EDDY,IERR)
          i = n
          CALL MPI_SEND(nprev+i,1,MPI_INTEGER,myright,1,MPI_COMM_EDDY,IERR)
          j = 0
          if(kmax>=kmin) j = kmax-kmin+1
          CALL MPI_SEND(nzprev+j,1,MPI_INTEGER,myright,2,MPI_COMM_EDDY,IERR)
c          write(6,*) 'myrank=',myrank,'sending',i,j,nprev+n,nzprev+j
        endif

      endif



      return

      end
C------------------------------------------------------------------------



C---- write_var_VPreg ----------------------- N. Beratlis-15 May 2011 ---
c
c     PURPOSE: Write time average variables to file
c
C------------------------------------------------------------------------
      subroutine write_var_VPreg(filename,uo,vo,wo,p,nx,ny,nz,ind
     $     ,j1,j2,kmin,kmax,nprev,nptsg,tlevel)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer npts,nsmpl,nprev,nx,ny,nz,kmin,kmax,nptsg,j1,j2
      real    tlevel
      integer ind(nz,2)
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
      character*(*) filename
c
c.... Local arrays
      integer i,k,ip,nj
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)
c
c.... Write time in header
      if(myrank.eq.0) then
        call mpi_type_contiguous(1,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,mtype,newtype
     &       ,'native',mpi_info_null,ierr)
        offset = 0
        call mpi_file_write_at(fh,offset,tlevel,1,newtype,status,ierr)

        call mpi_type_free(newtype,ierr)
        call mpi_file_close(fh,ierr)
      endif
c
c.... Write velocities and pressure
      nj = j2-j1+1
      call mpi_type_contiguous(nj,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

      call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mtype,newtype
     &     ,'native',mpi_info_null,ierr)

      if(kmax-kmin+1>0) then

        offset = 1+nprev
        do k=kmin,kmax
        do i=ind(k,1),ind(k,2)
          call mpi_file_write_at(fh,offset,uo(i,j1:j2,k),1,newtype,status,ierr)
          offset = offset+nj
        enddo
        enddo

        offset = 1+nptsg+nprev
        do k=kmin,kmax
        do i=ind(k,1),ind(k,2)
          call mpi_file_write_at(fh,offset,vo(i,j1:j2,k),1,newtype,status,ierr)
          offset = offset+nj
        enddo
        enddo

        offset = 1+2*nptsg+nprev
        do k=kmin,kmax
        do i=ind(k,1),ind(k,2)
          call mpi_file_write_at(fh,offset,wo(i,j1:j2,k),1,newtype,status,ierr)
          offset = offset+nj
        enddo
        enddo

        offset = 1+3*nptsg+nprev
        do k=kmin,kmax
        do i=ind(k,1),ind(k,2)
          call mpi_file_write_at(fh,offset,p(i,j1:j2,k),1,newtype,status,ierr)
          offset = offset+nj
        enddo
        enddo

      endif

      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------



C---- subroutine write_grid_VPreg -------- N. Beratlis - 27 Sep. 2010 ---
C
C     PURPOSE: Write grid for VP region
C
C------------------------------------------------------------------------
      subroutine write_grid_VPreg(filename,xu,yv,zw,xc,yc,zc,nx,ny,nz,ind,j1,j2,kmin,kmax,nzprev)
c
c      implicit none
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,j1,j2,kmin,kmax,nzprev
      integer ind(nz,2)
      real    xu(nx),yv(ny),zw(nz),xc(nx),yc(ny),zc(nz)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,nxl,nyl,nzl,imin,iming,imax,imaxg,nk,nkg,i1,i2,k1,k2
      integer fh
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)
c
      imin = nx
      imax = 0
      if(kmax>kmin) then
        imin = minval(ind(kmin:kmax,1))
        imax = maxval(ind(kmin:kmax,2))
      endif
c
      CALL MPI_REDUCE(imin,iming,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(imax,imaxg,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_EDDY,IERR)
c
      nk = 0
      if(kmax>kmin) nk = kmax-kmin+1
      CALL MPI_REDUCE(nk,nkg,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)      

      nxl = imax-imin+1
      nyl = j2-j1+1
      nzl = nkg
      i1 = iming
      i2 = imaxg

      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)      
      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

      offset = 0
      if(myrank.eq.0) then
        offset = 0
        call mpi_file_write_at(fh,offset,real(nxl),1,newtype,status,ierr)
        offset = 1
        call mpi_file_write_at(fh,offset,real(nyl),1,newtype,status,ierr)
        offset = 2
        call mpi_file_write_at(fh,offset,real(nzl),1,newtype,status,ierr)
        offset = 3

        do i=i1,i2
          call mpi_file_write_at(fh,offset,xu(i),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do i=i1,i2
          call mpi_file_write_at(fh,offset,xc(i),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do i=i1,i2
          call mpi_file_write_at(fh,offset,au(i),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do i=i1,i2
          call mpi_file_write_at(fh,offset,av(i),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do i=i1,i2
          call mpi_file_write_at(fh,offset,ap(i),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,yv(j),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,yc(j),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,bu(j),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,bv(j),1,newtype,status,ierr)
          offset = offset+1
        enddo

        do j=j1,j2
          call mpi_file_write_at(fh,offset,bp(j),1,newtype,status,ierr)
          offset = offset+1
        enddo

      endif

      if(nk>0) then

        offset = 3 + nxl*5 + nyl*5 + nzprev
        do k=kmin,kmax
c           write(6,*) k,offset,zw(k),nzl
          call mpi_file_write_at(fh,offset,zw(k),1,newtype,status,ierr)
          offset = offset+1
        enddo

        offset = 3 + nxl*5 + nyl*5 + nzprev + nzl
        do k=kmin,kmax
          call mpi_file_write_at(fh,offset,zc(k),1,newtype,status,ierr)
          offset = offset+1
        enddo

        offset = 3 + nxl*5 + nyl*5 + nzprev + 2*nzl
        do k=kmin,kmax
          call mpi_file_write_at(fh,offset,cu(k),1,newtype,status,ierr)
          offset = offset+1
        enddo

        offset = 3 + nxl*5 + nyl*5 + nzprev + 3*nzl
        do k=kmin,kmax
          call mpi_file_write_at(fh,offset,cw(k),1,newtype,status,ierr)
          offset = offset+1
        enddo

        offset = 3 + nxl*5 + nyl*5 + nzprev + 4*nzl
        do k=kmin,kmax
          call mpi_file_write_at(fh,offset,cp(k),1,newtype,status,ierr)
          offset = offset+1
        enddo

      endif

      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------

C---- subroutine read_wavg_top_prms-----------N. Beratlis-17 Feb. 2011---
C
C     PURPOSE: Read parameters for subroutines that average w on top boundary
C
C------------------------------------------------------------------------
      subroutine read_wavg_top_prms(flag,freq)
c
      implicit none
c
c.... Passed arrays
      integer flag,freq

      open(unit=10,file='wavg_top.input',form='formatted')
      read(10,'(A)') 
      read(10,*) flag
      read(10,'(A)')
      read(10,*) freq
      close(10)

      return

      end
C------------------------------------------------------------------------


C---- subroutine write_header_wavg_top--------N. Beratlis-17 Feb. 2011---
C
C     PURPOSE: Write header for wavg file
C
C------------------------------------------------------------------------
      subroutine write_header_wavg_top(filename,flag,w,ny,nz)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Passed arrays
      integer ny,nz,flag
      real    w(ny,nz)
      character*(*) filename
c
c.... Local arrays
      real    tmp
      integer filemode,fh
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      if(myrank.eq.0) then
        call mpi_type_contiguous(1,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)
        offset = 1
        call mpi_file_write_at(fh,offset,real(ny),1,newtype,status,ierr)
        offset = 2
        call mpi_file_write_at(fh,offset,real(myrank*(nz-2)+2),1,newtype,status,ierr)
        if(myrank.eq.0) then
          tmp = 0.0
          offset = 0
          call mpi_file_write_at(fh,offset,tmp,1,newtype,status,ierr)
        endif
        call mpi_file_close(fh,ierr)
        call mpi_type_free(newtype,ierr)
      endif
c
c.... Initialize or read w velocity
      if(flag==1) then
        w = 0.0
      elseif(flag>1) then
        offset = 3+ny*myrank*(nz-2)

        filemode = MPI_MODE_RDONLY
        call mpi_type_contiguous(ny*nz,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)
        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)
        call mpi_file_read_at_all(fh,offset,w,1,newtype,status,ierr)
        call mpi_file_close(fh,ierr)
        call mpi_type_free(newtype,ierr)
      endif
      
      return

      end
C------------------------------------------------------------------------



C---- subroutine calc_wavg_top----------------N. Beratlis-17 Dec. 2011---
C
C     PURPOSE: Calculate wavg at top boundary
C
C------------------------------------------------------------------------
      subroutine calc_wavg_top(wavg,w,nx,ny,nz,nsmpl)
c
c      implicit none
      include 'common.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,nsmpl
      real    w(nx,ny,nz),wavg(ny,nz)
c
c.... Local arrays

      wavg(:,:) = (wavg(:,:)*nsmpl + 0.5*(w(nx-1,:,:)+w(nx,:,:)) )/real(nsmpl+1)
      nsmpl = nsmpl+1
c      write(6,*) 'calc_wavg_top:, nsmpl=',nsmpl,myrank,w(nx-1,2,3),w(nx,2,3),wavg(2,3)      

      return

      end
C------------------------------------------------------------------------



C---- subroutine write_wavg_top--------------N. Beratlis-17 Dec. 2011----
C
C     PURPOSE: Write averaged velocity at top boundary to file
C
C------------------------------------------------------------------------
      subroutine write_wavg_top(filename,w,ny,nz,nsmpl)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer ny,nz,nsmpl
      real    w(ny,nz)
      character*(*) filename
c
c.... Local arrays
      integer n,j,k
      real    tmp
      integer filemode,fh
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      if(myrank.eq.0) then
        call mpi_type_contiguous(1,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        filemode = MPI_MODE_RDONLY
        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)
        offset = 0
        call mpi_file_read_at_all(fh,offset,tmp,1,newtype,status,ierr)
        n = int(tmp)
        call mpi_file_close(fh,ierr)

        filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)
        offset = 0
        call mpi_file_write_at(fh,offset,real(n+nsmpl),1,newtype,status,ierr)
        call mpi_file_close(fh,ierr)
        call mpi_type_free(newtype,ierr)
      endif

      call mpi_type_contiguous(ny*nz,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)
      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)
      offset = 3+myrank*ny*(nz-2)
c      write(6,*) 'myrank=',myrank,ny,nz,offset
      call mpi_file_write_at(fh,offset,w,1,newtype,status,ierr)
      call mpi_file_close(fh,ierr)
      call mpi_type_free(newtype,ierr)

c      open(unit=10,file='wavg_top.'//index(myrank+1),form='formatted')
c      do k=1,nz
c      do j=1,ny
c        write(10,*) k+myrank*(nz-2),j,w(j,k)
c      enddo
c      enddo
c      close(10)

      return

      end
C------------------------------------------------------------------------


C---- subroutine read_no_cf_probes-------------N. Beratlis-3 Feb. 2011---
C
C     PURPOSE: Read no of probes for shear stress for dimpled flat surface.
C
C------------------------------------------------------------------------
      subroutine read_no_cf_probes(nj,nk,icf,itcf)
c
      implicit none
c
c.... Input/Output arrays
      integer nj,nk,icf,itcf
c
      open(unit=10,file='cf.input',form='formatted')
      read(10,'(A)')
      read(10,*) icf
      read(10,'(A)')
      read(10,*) itcf
      read(10,'(A)')
      read(10,*) nj,nk
c      write(6,*) nj,nk
      close(10)

      return

      end
C------------------------------------------------------------------------


C---- subroutine setup_cf_probes---------------N. Beratlis-3 Feb. 2011---
C
C     PURPOSE: Setp ptobes for shear stress for dimpled flatsurf
C
C------------------------------------------------------------------------
      subroutine setup_cf_probes(filename,jcf,kcf,nj,nk,nkmax,nkprev,xc,yc,zw,nx,ny,nz,icf)
C
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nj,nk,nkmax,nkprev,nx,ny,nz,icf
      integer jcf(nj),kcf(nk)
      real    xc(nx),yc(ny),zw(nz)
      character*(*) filename
c
c.... Local arrays
      integer j,k,kk
      real    ycf(nj),zcf(nk)

      open(unit=10,file='cf.input',form='formatted')
      read(10,'(A)')
      read(10,*) icf
      read(10,'(A)')
      read(10,*)
      read(10,'(A)')
      read(10,*) 
      read(10,*) ycf
      read(10,*) zcf      
      close(10)

      do j=1,nj
        call closest(yc,ny,ycf(j),jcf(j))
      enddo

      kk=0
      do k=1,nk
        if(zcf(k)>=zw(1) .AND. zcf(k)<zw(nz-1)) then
          kk = kk+1
          call closest(zw,nz,zcf(k),kcf(kk))
          if(kk==1) nkprev=k
        endif
      enddo
      nkmax = kk

      call write_cf_probesheader(trim(filename),jcf,kcf,nj,nk,nkmax,nkprev,nz,icf)

      return

      end
C------------------------------------------------------------------------


C---- subroutine write_cf_probesheader---------N. Beratlis-3 Feb. 2011---
C
C     PURPOSE: Write header information for cf probe file
C
C------------------------------------------------------------------------
      subroutine write_cf_probesheader(filename,jcf,kcf,nj,nk,nkmax,nkprev,nz,icf)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nj,nk,nkmax,nkprev,nz,icf
      integer jcf(nj),kcf(nk)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,fh
      real    tmp
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)
      
      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)      

      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

      if(myrank.eq.0) then
        offset = 1
        call mpi_file_write_at(fh,offset,real(nj),1,newtype,status,ierr)
        offset = 2
        call mpi_file_write_at(fh,offset,real(nk),1,newtype,status,ierr)
        offset = 3

        do j=1,nj
          call mpi_file_write_at(fh,offset,real(jcf(j)),1,newtype,status,ierr)
          offset = offset+1
        enddo
        if(icf==1) then
          offset = 0
          tmp = 0.0
          call mpi_file_write_at(fh,offset,tmp,1,newtype,status,ierr)
        endif
      endif

      offset = 3+nj+(nkprev-1)
      do k=1,nkmax
        call mpi_file_write_at(fh,offset,real(kcf(k)+myrank*(nz-2)),1,newtype,status,ierr)
        offset = offset+1
      enddo
      call mpi_file_close(fh,ierr)

      return

      end
C------------------------------------------------------------------------


C---- subroutine record_cf_flatsurf------------N. Beratlis-3 Feb. 2011---
C
C     PURPOSE: Record cf at probes
C
C------------------------------------------------------------------------
      subroutine record_cf_flatsurf(wo,xc,yc,zw,nx,ny,nz,jcf,kcf,nj,nk
     $     ,nkmax,cfprb,tcf,ncf,it,tlevel)
c
c      implicit none
      include 'common.h'
      include 'immersed.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,nj,nk,nkmax,ncf,it
      real    tlevel
      integer jcf(nj),kcf(nk)
      real    xc(nx),yc(ny),zw(nz)
      real    wo(nx,ny,nz)
      real    cfprb(nj,nk,ncf),tcf(ncf)
c
c.... Local arrays
      integer i,j,k,i1,i2
      real    a1,a2,xplate,dx,x1,x2,w1,w2
      integer, dimension(:), allocatable :: indx
      real, dimension(:,:), allocatable :: amtrx
      real, dimension(:), allocatable :: b,rmtrx

      if(ivelrcn>1) allocate(indx(ivelrcn),amtrx(ivelrcn,ivelrcn),b(ivelrcn),rmtrx(ivelrcn))

      xplate = dummyfp4
      call locate(xc,nx,xplate,i)
      dx = xc(i+1)-xc(i)
      x1 = xplate+uext1*dx
      call locate(xc,nx,x1,i1)
      a1 = (xc(i1+1)-x1)/(xc(i1+1)-xc(i1))
      if(ivelrcn>1) then
        x2 = xplate+uext2*dx
        call locate(xc,nx,x2,i2)
        a2 = (xc(i2+1)-x2)/(xc(i2+1)-xc(i2))
        amtrx(1,1) = (x1-xplate)
        amtrx(1,2) = (x1-xplate)**2.0
        amtrx(2,1) = (x2-xplate)
        amtrx(2,2) = (x2-xplate)**2.0
c        amtrx(2,1) = 1.0
c        amtrx(2,2) = 2.0*(x1-xplate)
        call ludcmp(amtrx,ivelrcn,ivelrcn,indx,b)
      endif

c      if(nkmax>0) then
c        open(unit=10,file='record_cf.'//index(myrank)//'.txt',form
c     $     ='formatted',position='append')
c      endif

      it = it+1
      tcf(it) = tlevel
      do k=1,nkmax
      do j=1,nj
        w1 = a1*wo(i1,jcf(j),kcf(k)) + (1.0-a1)*wo(i1+1,jcf(j),kcf(k))
        if(ivelrcn>1) then
          w2 = a2*wo(i2,jcf(j),kcf(k)) + (1.0-a2)*wo(i2+1,jcf(j),kcf(k))
          rmtrx(1) = w1
          rmtrx(2) = w2
          amtrx(1,1) = (x1-xplate)
          amtrx(1,2) = (x1-xplate)**2.0
          amtrx(2,1) = (x2-xplate)
          amtrx(2,2) = (x2-xplate)**2.0
          call ludcmp(amtrx,ivelrcn,ivelrcn,indx,b)
          call lubksb(amtrx,ivelrcn,ivelrcn,indx,rmtrx)
          cfprb(j,k,it) = ru1*rmtrx(1)
c          write(10,'(3(1x,I6),8(1x,F16.10))') it,k,j,w1,w2,x1,x2,rmtrx,cfprb(j,k,it),ru1*w1/(x1-xplate)
        else
          cfprb(j,k,it) = ru1*w1/(x1-xplate)
        endif
      enddo
      enddo
      
c      if(nkmax>0) close(10)

      if(ivelrcn>1) deallocate(indx,amtrx,rmtrx,b)

      return

      end
C------------------------------------------------------------------------


C---- subroutine write_cf_flatsurf-----------N. Beratlis-03 Feb. 2011 ---
C
C     PURPOSE: Write Cf  to file
C
C------------------------------------------------------------------------
      subroutine write_cf_flatsurf(filename,cfprb,tcf,nj,nk,nkmax,nkprev,ncfprb,ncfsmpl)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nj,nk,nkmax,nkprev,ncfprb,ncfsmpl
      real    cfprb(nj,nk,ncfprb),tcf(ncfprb)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,nt,fh,it
      real    tmp
      integer filemode
      integer newtype
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)
c
c.... Read sample size
      if(myrank.eq.0) then
        call mpi_type_contiguous(1,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        filemode = MPI_MODE_RDONLY
        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)
        offset = 0
        call mpi_file_read_at(fh,offset,tmp,1,newtype,status,ierr)
        it = int(tmp)
        call mpi_file_close(fh,ierr)
        call mpi_type_free(newtype,ierr)
      endif

      CALL MPI_BCAST(IT,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)

      if(nkmax>0) then
c        open(unit=10,file='write_cf.'//index(myrank)//'.txt',form
c     $        ='formatted',position='append')
       

        call mpi_type_contiguous(nj,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 3*8 + (nj+nk)*8
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

        do i=1,ncfprb
          offset = it*(nj*nk+1) + (i-1)*(nj*nk+1) + (nkprev-1)*nj + 1
        do k=1,nkmax
c          write(10,'(7(I6),10(1x,E14.6))') i,k,nj,nk,nkprev,it,offset,cfprb(:,k,i)
          call mpi_file_write_at(fh,offset,cfprb(:,k,i),1,newtype,status,ierr)
          offset = offset+nj
        enddo
        enddo

        call mpi_type_free(newtype,ierr)
        call mpi_file_close(fh,ierr)
c        close(10)
      endif
c
c.... Write sample size to header      
      if(myrank.eq.0) then
        call mpi_type_contiguous(1,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
        call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        offset = 0
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)
        call mpi_file_write_at(fh,offset,real(it+ncfprb),1,newtype,status,ierr)
        offset = 3+nj+nk + it*(nj*nk+1)
        do i=1,ncfprb
c          write(6,*) 'timecf',i,nj,nk,it,offset,tcf(i)
          call mpi_file_write_at(fh,offset,tcf(i),1,newtype,status,ierr)
          offset = offset + nj*nk + 1
        enddo
        call mpi_type_free(newtype,ierr)
        call mpi_file_close(fh,ierr)
      endif

      ncfsmpl = 0

      return

      end
C------------------------------------------------------------------------
