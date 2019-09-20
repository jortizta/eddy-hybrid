ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
      subroutine plotinststagg_xyz(iplant9c,iplant9d,iplant9p,
     %nx,ny,nz,nzg,nbd,ntime,xc,xu,yc,yv,zcg,zwg,p,vo,uo,wo)
!     %,flaguo,flagvo,flagwo)
c
      use sections
      include'common.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer iplant9c,iplant9d,iplant9p,nx,ny,nz,nzg,nbd,ntime
!      integer flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),
!     %flagwo(nx,ny,nz,nbd)
      real xc(nx),xu(nx),yc(ny),yv(ny),zcg(nzg),zwg(nzg)
      real p(nx,ny,nz),vo(nx,ny,nz),uo(nx,ny,nz),wo(nx,ny,nz)
c
c Local variables
      integer iic,iicd,iicu,is,im,ip,kkc,kkcd,kkcu,ks,km,kp,
     %jjc,jjcd,jjcu,js,jm,jp,i,j,k,l,ir,jr,kr,lr,
     %i1,jjp,kk1,kk2,kk,kkd,kku
      integer STATUS(MPI_STATUS_SIZE)
!      integer flaguor(ny,nzg,nbd),
!     %flagvor2(nx,nzg,nbd)
      real pstagg,vstagg,ustagg,wstagg
      real pr(ny,nzg),vor(ny,nzg),uor(ny,nzg),wor(ny,nzg),
     %pr2(nx,nzg),vor2(nx,nzg),uor2(nx,nzg),wor2(nx,nzg)
      integer kg
c
c Parameters
      integer lamb2,nxp,nyp,nzp,nxpc,nypc,nzpc
      real deltasect
      parameter (lamb2=1)
      parameter (nyp=1,nxp=1,nzp=1)
      parameter (deltasect=2.*pi)
c
   8  format(10f25.8)
   9  format(10f25.8)
 101  format(3i8)
 102  format(i8,3i2)
 103  format(i2)
c
      if(mod(nx-2,nxp).eq.0) then
        nxpc = (nx-2)/nxp
      else
        nxpc = (nx-2)/nxp+1
      endif
      if(mod(ny-2,nyp).eq.0) then
        nypc = (ny-2)/nyp
      else
        nypc = (ny-2)/nyp+1
      endif
      if(mod(nz-2,nzp).eq.0) then
        nzpc = (nzg-2)/nzp
      else
        nzpc = ((nz-2)/nzp+1)*mysize
      endif
c
      if(iplant9.eq.1) then
         if(myrank.eq.0) then
           do iic=1,nsez1
           iicd=iic/10
           iicu=mod(iic,10)
           open
     %(81,file='grid1_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'_plotinst.xyz')
           is=isez1(iic)
!           write(81,101)lamb2*(ny-2)/nyp,1,(nzg-2)/nzp
           write(81,101)lamb2*nypc,1,nzpc
!           do 10 k=2,nzg-1,nzp
           do 10 k=1,nzpc
           do 10 l=1,lamb2
           do 10 j=jy1,jy2,nyp
   10      write(81,9)xu(is)
!           do 20 k=2,nzg-1,nzp
           do 20 k=1,nzpc
           do 20 l=1,lamb2
           do 20 j=jy1,jy2,nyp
   20      write(81,9)yc(j)
!           do 30 k=2,nzg-1,nzp
           do 30 i=1,mysize
           do 30 k=2,nz-1,nzp
           kg=k+(i-1)*(nz-2)
           do 30 l=1,lamb2
           do 30 j=jy1,jy2,nyp
!   30      write(81,9)zcg(k)
   30      write(81,9)zcg(kg)
           close(81)
           enddo
           do kkc=1,nsez1
           kkcd=kkc/10
           kkcu=mod(kkc,10)
           open(81,file='grid2_2D_slice'//char(48+kkcd)
     %//char(48+kkcu)//'_plotinst.xyz')
           ks=ksez2g(kkc)
!           write(81,101)lamb2*(ny-2)/nyp,(nx-2)/nxp,1
           write(81,101)lamb2*nypc,nxpc,1
           do 110 i=ix1,ix2,nxp
           do 110 l=1,lamb2
           do 110 j=jy1,jy2,nyp
  110      write(81,9)xc(i)
           do 120 i=ix1,ix2,nxp
           do 120 l=1,lamb2
           do 120 j=jy1,jy2,nyp
  120      write(81,9)yc(j)
           do 130 i=ix1,ix2,nxp
           do 130 l=1,lamb2
           do 130 j=jy1,jy2,nyp 
  130      write(81,9)zwg(ks)
           close(81)
           enddo
           do jjc=1,nsez14
           jjcd=jjc/10
           jjcu=mod(jjc,10)
           open(81,file='grid3_2D_slice'//char(48+jjcd)
     %//char(48+jjcu)//'_plotinst.xyz')
           js=jsez14(jjc)
!           write(81,101)1,(nx-2)/nxp,(nzg-2)/nzp
           write(81,101)1,nxpc,nzpc
!           do 210 k=2,nzg-1,nzp
           do 210 k=1,nzpc
           do 212 i=ix1,ix2,nxp
  212      write(81,9)xc(i)
  210      continue
!           do 220 k=2,nzg-1,nzp
           do 220 k=1,nzpc 
           do 222 i=ix1,ix2,nxp
  222      write(81,9)yv(js)
  220      continue
!           do 230 k=2,nzg-1,nzp
           do 230 l=1,mysize
           do 230 k=2,nz-1,nzp
           kg=k+(l-1)*(nz-2)
           do 230 i=1,nx-2,nxp
!  230	       write(81,9)zcg(k)
  230      write(81,9)zcg(kg)
           close(81)
           enddo
        endif
      endif
c
      iplant9p=iplant9
      iplant9c=iplant9p/100
      iplant9p=mod(iplant9p,100)
      iplant9d=iplant9p/10
      iplant9p=mod(iplant9p,10)
c
      do iic=1,nsez1
	 is=isez1(iic)
	 ip=is+1
	 if(myrank.eq.0) then
	    iicd=iic/10
	    iicu=mod(iic,10)
	    open(82,file='results1_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'_plotinst.q')
!        write(82,101)lamb2*(ny-2)/nyp,1,(nzg-2)/nzp
        write(82,101)lamb2*nypc,1,nzpc
	    i1=1
	    write(82,102)ntime,i1,i1,i1
!
! Stagnation pressure
!
	    do 40 k=kz1,kz2,nzp
	    km=k-1
	    do 41 l=1,lamb2
	    do 41 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       jm=jmv(j)
	       pstagg=
     %0.5*(p(is,j,k)+p(ip,j,k))
	       vstagg=
     %0.25*(vo(is,j,k)+vo(ip,j,k)+vo(is,jm,k)+vo(ip,jm,k))
	       ustagg=
     %uo(is,j,k)
	       wstagg=
     %0.25*(wo(is,j,k)+wo(ip,j,k)+wo(is,j,km)+wo(ip,j,km))
	       write(82,9)
     %(pstagg+0.5*(vstagg**2.+ustagg**2.+
     %wstagg**2.))	!*ros
c	    else
c	       write(82,9)1.e+03
c	    endif
 41	    continue
 40	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
c		  CALL MPI_RECV
c    %(flaguor(jy1:jy2,kk1:kk2,:),
c    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,jjp,1,
c    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(pr(jy1:jy2,kk1:kk2),1*(ny-2)*(nz-2),MTYPE,jjp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(vor(jy1:jy2,kk1:kk2),1*(ny-2)*(nz-2),MTYPE,jjp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(uor(jy1:jy2,kk1:kk2),1*(ny-2)*(nz-2),MTYPE,jjp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(wor(jy1:jy2,kk1:kk2),1*(ny-2)*(nz-2),MTYPE,jjp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 42 kr=kk1,kk2,nzp
		  do 43 lr=1,lamb2
		  do 43 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)
     %(pr(jr,kr)+0.5*(vor(jr,kr)**2.+uor(jr,kr)**2.+
     %wor(jr,kr)**2.))	    !*ros
c		  else
c		     write(82,9)1.e+03
c		  endif
 43		  continue
 42		  continue
	       enddo
	    endif
!	    write(82,*)	  !!!!!!
!	    write(82,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp),i=1,1),k=1,(nzg-2)/nzp)
!
! X velocity
!
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
c           if(all(flaguo(is,j,k,:)).le.0) then
               write(82,9)uo(is,j,k)
c           else
c              write(82,9)1.e+03
c           endif
   51	    continue
   50	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 52 kr=kk1,kk2,nzp
		  do 53 lr=1,lamb2
		  do 53 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)uor(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   53		  continue
   52		  continue
	       enddo
	    endif
!
! Y velocity
!
	    do 60 k=kz1,kz2,nzp
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2,nyp
c           if(all(flaguo(is,j,k,:)).le.0) then
               jm=jmv(j)
               write(82,9)
     %0.25*(vo(is,j,k)+vo(ip,j,k)+vo(is,jm,k)+vo(ip,jm,k))
c           else
c              write(82,9)1.e+03
c           endif
   61	    continue
   60	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 62 kr=kk1,kk2,nzp
		  do 63 lr=1,lamb2
		  do 63 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)vor(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   63		  continue
   62		  continue
	       enddo
	    endif
!
!Axial velocity
!
	    do 70 k=kz1,kz2,nzp
	    km=k-1
	    do 71 l=1,lamb2
	    do 71 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)
     %0.25*(wo(is,j,k)+wo(ip,j,k)+wo(is,j,km)+wo(ip,j,km))
c	    else
c	       write(82,9)1.e+03
c	    endif
   71	    continue
   70	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 72 kr=kk1,kk2,nzp
		  do 73 lr=1,lamb2
		  do 73 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)wor(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   73		  continue
   72		  continue
	       enddo
	    endif
!
!Static pressure
!
	    do 80 k=kz1,kz2,nzp
	    do 81 l=1,lamb2
	    do 81 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,8)0.5*(p(is,j,k)+p(ip,j,k))   !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
   81	    continue
   80	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 82 kr=kk1,kk2,nzp
		  do 83 lr=1,lamb2
		  do 83 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,8)pr(jr,kr)   !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
   83		  continue
   82		  continue
	       enddo
	    endif
	    close(82)
	 else
c	    CALL MPI_SEND
c    %(flaguo(is,jy1:jy2,kz1:kz2,:),
c    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,0,1,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(0.5*(p(is,jy1:jy2,kz1:kz2)+p(ip,jy1:jy2,kz1:kz2)),
     %1*(ny-2)*(nz-2),MTYPE,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(0.25*(vo(is,1:jy2-1,kz1:kz2)+vo(ip,1:jy2-1,kz1:kz2)+
     %vo(is,jy1:jy2,kz1:kz2)+vo(ip,jy1:jy2,kz1:kz2)),
     %1*(ny-2)*(nz-2),MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uo(is,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(0.25*(wo(is,jy1:jy2,1:kz2-1)+wo(ip,jy1:jy2,1:kz2-1)+
     %wo(is,jy1:jy2,kz1:kz2)+wo(ip,jy1:jy2,kz1:kz2)),
     %1*(ny-2)*(nz-2),MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
      do kkc=1,nsez1
	 ks=ksez2(kkc)
	 kp=ks+1
	 kk=isezindx2+(kkc-1)
	 kkd=kk/10
	 kku=mod(kk,10)
	 open(82,file='results2_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_plotinst.q')
!         write(82,101)lamb2*(ny-2)/nyp,(nx-2)/nxp,1
         write(82,101)lamb2*nypc,nxpc,1
	 i1=1
	 write(82,102)ntime,i1,i1,i1
!
!Stagnation pressure
!
	 do 140 i=ix1,ix2,nxp
	 im=i-1
	 do 141 l=1,lamb2
	 do 141 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    jm=jmv(j)
	    pstagg=
     %0.5*(p(i,j,ks)+p(i,j,kp))
	    vstagg=
     %0.25*(vo(i,j,ks)+vo(i,jm,ks)+vo(i,j,kp)+vo(i,jm,kp))
	    ustagg=
     %0.25*(uo(i,j,ks)+uo(im,j,ks)+uo(i,j,kp)+uo(im,j,kp))
	    wstagg=
     %wo(i,j,ks)
	    write(82,9)
     %(pstagg+0.5*(vstagg**2.+ustagg**2.+wstagg**2.))	!*ros
c	 else
c	    write(82,9)1.e+03
c	 endif
 141	 continue
 140	 continue
!	 write(82,*)	!!!!!!
!	 write(82,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp),i=1,(nx-2)/nxp),k=1,1)
!
! X velocity
!
	 do 150 i=ix1,ix2,nxp
	 im=i-1
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
c        if(all(flagwo(i,j,ks,:)).le.0) then
            write(82,9)
     %0.25*(uo(i,j,ks)+uo(im,j,ks)+uo(i,j,kp)+uo(im,j,kp))
c        else
c           write(82,9)1.e+03
c        endif
 151	 continue
 150	 continue
!
! Y velocity
!
	 do 160 i=ix1,ix2,nxp
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2,nyp
c        if(all(flagwo(i,j,ks,:)).le.0) then
            jm=jmv(j)
            write(82,9)
     %0.25*(vo(i,j,ks)+vo(i,jm,ks)+vo(i,j,kp)+vo(i,jm,kp))
c        else
c           write(82,9)1.e+03
c        endif
 161	 continue
 160	 continue
!
!Axial velocity
!
	 do 170 i=ix1,ix2,nxp
	 do 171 l=1,lamb2
	 do 171 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)wo(i,j,ks)
c	 else
c	    write(82,9)1.e+03
c	 endif
 171	 continue
 170	 continue
!
!Static pressure
!
	 do 180 i=ix1,ix2,nxp
	 do 181 l=1,lamb2
	 do 181 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,8)0.5*(p(i,j,ks)+p(i,j,kp))    !*ros
c	 else
c	    write(82,8)1.e+03
c	 endif
 181	 continue
 180	 continue
	 close(82)
      enddo
c

      do jjc=1,nsez14
	 js=jsez14(jjc)
	 jp=jpv(js)
	 if(myrank.eq.0) then
	    jjcd=jjc/10
	    jjcu=mod(jjc,10)
	    open(82,file='results3_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_plotinst.q')
!        write(82,101)1,(nx-2)/nxp,(nzg-2)/nzp
        write(82,101)1,nxpc,nzpc
	    i1=1
	    write(82,102)ntime,i1,i1,i1
!
!Stagnation pressure
!
	    do 240 k=kz1,kz2,nzp
	    km=k-1
	    do 242 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       im=i-1
	       pstagg=
     %0.5*(p(i,js,k)+p(i,jp,k))
	       vstagg=
     %vo(i,js,k)
	       ustagg=
     %0.25*(uo(i,js,k)+uo(im,js,k)+uo(i,jp,k)+uo(im,jp,k))
	       wstagg=
     %0.25*(wo(i,js,k)+wo(i,jp,k)+wo(i,js,km)+wo(i,jp,km))
	       write(82,9)
     %(pstagg+0.5*(vstagg**2.+ustagg**2.+wstagg**2.))	 !*ros
c	    else
c	       write(82,9)1.e+03
c	    endif
 242	    continue
 240	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
c		  CALL MPI_RECV
c    %(flagvor2(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,6,
c    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(pr2(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,7,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(vor2(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(uor2(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(wor2(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 243 kr=kk1,kk2,nzp
		  do 245 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)
     %(pr2(ir,kr)+0.5*(vor2(ir,kr)**2.+uor2(ir,kr)**2.
     %+wor2(ir,kr)**2.))    !*ros
c		  else
c		     write(82,9)1.e+03
c		  endif
 245		  continue
 243		  continue
	       enddo
	    endif


! X velocity
!
	    do 250 k=kz1,kz2,nzp
	    do 252 i=ix1,ix2,nxp
c           if(all(flagvo(i,js,k,:)).le.0) then
               im=i-1
               write(82,9)
     %0.25*(uo(i,js,k)+uo(im,js,k)+uo(i,jp,k)+uo(im,jp,k))
c           else
c              write(82,9)1.e+03
c           endif
 252	    continue
 250	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 253 kr=kk1,kk2,nzp
		  do 255 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)uor2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 255		  continue
 253		  continue
	       enddo
	    endif
!
! Y velocity
!
	    do 260 k=kz1,kz2,nzp
	    do 262 i=ix1,ix2,nxp
c           if(all(flagvo(i,js,k,:)).le.0) then
               write(82,9)vo(i,js,k)
c           else
c              write(82,9)1.e+03
c           endif
 262	    continue
 260	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 263 kr=kk1,kk2,nzp
		  do 265 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)vor2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 265		  continue
 263		  continue
	       enddo
	    endif
!
!Axial velocity
!
	    do 270 k=kz1,kz2,nzp
	    km=k-1
	    do 272 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,9)
     %0.25*(wo(i,js,k)+wo(i,jp,k)+wo(i,js,km)+wo(i,jp,km))
c	    else
c	       write(82,9)1.e+03
c	    endif
 272	    continue
 270	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 273 kr=kk1,kk2,nzp
		  do 275 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)wor2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 275		  continue
 273		  continue
	       enddo
	    endif
!
!Static pressure
!
	    do 280 k=kz1,kz2,nzp
	    do 282 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,8)0.5*(p(i,js,k)+p(i,jp,k))    !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
 282	    continue
 280	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 283 kr=kk1,kk2,nzp
		  do 285 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,8)pr2(ir,kr)    !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
 285		  continue
 283		  continue
	       enddo
	    endif
	    close(82)
	 else
c	    CALL MPI_SEND
c    %(flagvo(ix1:ix2,js,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,6,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((p(ix1:ix2,js,kz1:kz2)+p(ix1:ix2,jp,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vo(ix1:ix2,js,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((uo(1:ix2-1,js,kz1:kz2)+uo(ix1:ix2,js,kz1:kz2)+
     %uo(1:ix2-1,jp,kz1:kz2)+uo(ix1:ix2,jp,kz1:kz2))*0.25,
     %(nx-2)*1*(nz-2),MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((wo(ix1:ix2,js,1:kz2-1)+wo(ix1:ix2,jp,1:kz2-1)+
     %wo(ix1:ix2,js,kz1:kz2)+wo(ix1:ix2,jp,kz1:kz2))*0.25,
     %(nx-2)*1*(nz-2),MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vort2stagg_xyz(nx,ny,nz,xc,xu,zc,uo,vo,wo)
c
      use sections
      use vorticity
      include'common.h'
      include'averages.h'
c
c Global variables
      integer nx,ny,nz
      real xc(nx),xu(nx),zc(nz),
     %uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
c
c Local variables
      integer iic,ic,ip,im,jjc,jc,jp,jm,jsy,jpsy,kkc,kc,kp,km
      real dtheta,dr,dz,dq1x3,dq3x1,dq2x3,dq3x2,dq2x1,dq1x2
c
      ! Y vorticity - Circumferential sections !
      do 101 kc=kz1,kz2
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	do 101 iic=1,nsez1
	  ic=isez1(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 101 jc=jy1,jy2
	    dq1x3=(uo(ic,jc,kp)-uo(ic,jc,km))/dz
	    dq3x1=(0.5*(wo(ip,jc,kc)+wo(ip,jc,km))
     %-0.5*(wo(ic,jc,kc)+wo(ic,jc,km)))/dr
	    voraz1(iic,jc,kc)=dq1x3-dq3x1
  101 continue
c
      ! Y vorticity - Cross sections !
      do 111 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
        if(icyl.eq.0) goto 114
	ic=ix1
	ip=ic+1
	im=ic-1
	dr=xc(ip)-xc(im)
	do 112 jc=jy1,jy2
	  jsy=jsym(jc)
	  dq1x3=(0.5*(uo(ic,jc,kp)+uo(im,jc,kp))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dz
	  dq3x1=(wo(ip,jc,kc)
     %-wo(ic,jsy,kc))/dr
	  voraz2(ic,jc,kkc)=dq1x3-dq3x1
  112 continue
  114 continue 
	do 113 ic=ix1+icyl,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 113 jc=jy1,jy2
	    dq1x3=(0.5*(uo(ic,jc,kp)+uo(im,jc,kp))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dz
	    dq3x1=(wo(ip,jc,kc)
     %-wo(im,jc,kc))/dr
	    voraz2(ic,jc,kkc)=dq1x3-dq3x1
  113 continue
  111 continue
c
      ! Y vorticity - Meridian sections !
      do 121 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jsy=jsym(jc)
	jpsy=jsym(jp)
	do 121 kc=kz1,kz2
	  kp=kc+1
	  km=kc-1
	  dz=zc(kp)-zc(km)
          if(icyl.eq.0) goto 122
	  ic=ix1
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  dq1x3=(0.25*(uo(ic,jc,kp)+uo(im,jc,kp)+
     %uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,km)+uo(im,jc,km)+
     %uo(ic,jp,km)+uo(im,jp,km)))/dz
	  dq3x1=(0.25*(wo(ip,jc,kc)+wo(ip,jp,kc)+
     %wo(ip,jc,km)+wo(ip,jp,km))
     %-0.25*(wo(ic,jsy,kc)+wo(ic,jpsy,kc)+
     %wo(ic,jsy,km)+wo(ic,jpsy,km)))/dr
	  voraz14(ic,jjc,kc)=dq1x3-dq3x1
 122      continue
	  do 121 ic=ix1+icyl,ix2
	    ip=ic+1
	    im=ic-1
	    dr=xc(ip)-xc(im)
	    dq1x3=(0.25*(uo(ic,jc,kp)+uo(im,jc,kp)+
     %uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,km)+uo(im,jc,km)+
     %uo(ic,jp,km)+uo(im,jp,km)))/dz
	    dq3x1=(0.25*(wo(ip,jc,kc)+wo(ip,jp,kc)+
     %wo(ip,jc,km)+wo(ip,jp,km))
     %-0.25*(wo(im,jc,kc)+wo(im,jp,kc)+
     %wo(im,jc,km)+wo(im,jp,km)))/dr
	    voraz14(ic,jjc,kc)=dq1x3-dq3x1
  121 continue
c
      ! X vorticity - Circumferential sections !
      do 131 kc=kz1,kz2
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	do 131 iic=1,nsez1
	  ic=isez1(iic)
	  ip=ic+1
	  do 131 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq3x2=(0.25*(wo(ic,jp,kc)+wo(ip,jp,kc)+
     %wo(ic,jp,km)+wo(ip,jp,km))
     %-0.25*(wo(ic,jm,kc)+wo(ip,jm,kc)+
     %wo(ic,jm,km)+wo(ip,jm,km)))/dtheta     !
	    dq2x3=(0.25*(vo(ic,jc,kp)+vo(ip,jc,kp)+
     %vo(ic,jm,kp)+vo(ip,jm,kp))
     %-0.25*(vo(ic,jc,km)+vo(ip,jc,km)+
     %vo(ic,jm,km)+vo(ip,jm,km)))/dz
	    vorr1(iic,jc,kc)=dq3x2/ru(ic)-dq2x3
  131 continue
c
      ! X vorticity - Cross sections !
      do 141 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
	do 141 ic=ix1,ix2
	  do 141 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq3x2=(wo(ic,jp,kc)
     %-wo(ic,jm,kc))/dtheta	!
	    dq2x3=(0.5*(vo(ic,jc,kp)+vo(ic,jm,kp))
     %-0.5*(vo(ic,jc,kc)+vo(ic,jm,kc)))/dz
	    vorr2(ic,jc,kkc)=dq3x2/rp(ic)-dq2x3
  141 continue
c
      ! X vorticity - Meridian sections !
      do 151 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	dtheta=dely
	do 151 kc=kz1,kz2
	  kp=kc+1
	  km=kc-1
	  dz=zc(kp)-zc(km)
	  do 151 ic=ix1,ix2
	    dq3x2=(0.5*(wo(ic,jp,kc)+wo(ic,jp,km))
     %-0.5*(wo(ic,jc,kc)+wo(ic,jc,km)))/dtheta	  !
	    dq2x3=(vo(ic,jc,kp)
     %-vo(ic,jc,km))/dz
	    vorr14(ic,jjc,kc)=dq3x2/rp(ic)-dq2x3
  151 continue
c
      ! Z vorticity - Circumferential sections !
      do 161 kc=kz1,kz2
	do 161 iic=1,nsez1
	  ic=isez1(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 161 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.5*(vo(ip,jc,kc)+vo(ip,jm,kc))*rp(ip)
     %-0.5*(vo(ic,jc,kc)+vo(ic,jm,kc))*rp(ic))/dr
	    dq1x2=(uo(ic,jp,kc)
     %-uo(ic,jm,kc))/dtheta    !
	    vorz1(iic,jc,kc)=(dq2x1-dq1x2)/ru(ic)
  161 continue
c
      ! Z vorticity - Cross sections !
      do 171 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
        if(icyl.eq.0) goto 174
	ic=ix1
	im=ic-1
	ip=ic+1
	dr=xu(ic)-xu(im)
	do 172 jc=jy1,jy2
	  jp=jpv(jc)
	  jm=jmv(jc)
	  dtheta=2.*dely
	  dq2x1=(0.125*(vo(ic,jc,kc)+vo(ip,jc,kc)
     %+vo(ic,jm,kc)+vo(ip,jm,kc)
     %+vo(ic,jc,kp)+vo(ip,jc,kp)
     %+vo(ic,jm,kp)+vo(ip,jm,kp))*ru(ic))/dr
	  dq1x2=(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)+
     %uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jm,kc)+uo(im,jm,kc)+
     %uo(ic,jm,kp)+uo(im,jm,kp)))/dtheta    !
	  vorz2(ic,jc,kkc)=(dq2x1-dq1x2)/rp(ic)
  172	continue
  174   continue
c
	do 173 ic=ix1+icyl,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 173 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.25*(vo(ip,jc,kc)+vo(ip,jm,kc)+
     %vo(ip,jc,kp)+vo(ip,jm,kp))*rp(ip)
     %-0.25*(vo(im,jc,kc)+vo(im,jm,kc)+
     %vo(im,jc,kp)+vo(im,jm,kp))*rp(im))/dr
	    dq1x2=
     %(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)+
     %uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jm,kc)+uo(im,jm,kc)+
     %uo(ic,jm,kp)+uo(im,jm,kp)))/dtheta   !
	    vorz2(ic,jc,kkc)=(dq2x1-dq1x2)/rp(ic)
  173 continue
  171 continue
c
      ! Z vorticity - Meridian sections !
      do 181 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	dtheta=dely
	do 181 kc=kz1,kz2
          if(icyl.eq.0) goto 182
	  ic=ix1
	  im=ic-1
	  ip=ic+1
	  dr=xu(ic)-xu(im)
	  dq2x1=(0.5*(vo(ip,jc,kc)+vo(ic,jc,kc))*ru(ic))/dr
	  dq1x2=(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dtheta	  !
	  vorz14(ic,jjc,kc)=(dq2x1-dq1x2)/rp(ic)
 182      continue
	  do 181 ic=ix1+icyl,ix2
	    ip=ic+1
	    im=ic-1
	    dr=xc(ip)-xc(im)
	    dq2x1=(vo(ip,jc,kc)*rp(ip)
     %-vo(im,jc,kc)*rp(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dtheta	  !
	    vorz14(ic,jjc,kc)=(dq2x1-dq1x2)/rp(ic)
  181 continue
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine plotvortstagg_xyz(iplant9c,iplant9d,iplant9p,
     %nx,ny,nz,nzg,nbd,ntime,p,uo,vo,wo,tv)
!     %,flaguo,flagvo,flagwo)
c
      use sections
      use vorticity
      include'common.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer iplant9c,iplant9d,iplant9p,nx,ny,nz,nzg,nbd,ntime
!      integer flaguo(nx,ny,nz,nbd),flagwo(nx,ny,nz,nbd),
!     %flagvo(nx,ny,nz,nbd)
      real p(nx,ny,nz),uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),
     %tv(nx,ny,nz)
c
c Local variables
      integer iic,iicd,iicu,is,ip,jjc,jjcd,jjcu,js,jp,
     %kkc,ks,kp,i,j,k,l,ir,jr,kr,lr,i1,jjp,kk1,kk2,
     %kk,kkd,kku,krg
      integer STATUS(MPI_STATUS_SIZE)
!      integer flaguor(ny,nzg,nbd),
!     %flagvor2(nx,nzg,nbd)
      real recvar(ny,nz),recvar2(nx,nz)
c
c Parameters
      integer lamb2,nxp,nyp,nzp,nxpc,nypc,nzpc
      parameter (lamb2=1)
      parameter (nyp=1,nxp=1,nzp=1)
c
   8  format(10f25.8)
   9  format(10f25.8)
 101  format(3i8)
 102  format(i8,3i2)
 103  format(i2)
c
      if(mod(nx-2,nxp).eq.0) then
        nxpc = (nx-2)/nxp
      else
        nxpc = (nx-2)/nxp+1
      endif
      if(mod(ny-2,nyp).eq.0) then
        nypc = (ny-2)/nyp
      else
        nypc = (ny-2)/nyp+1
      endif
      if(mod(nz-2,nzp).eq.0) then
        nzpc = (nzg-2)/nzp
      else
        nzpc = ((nz-2)/nzp+1)*mysize
      endif
c
      i1=1
      do iic=1,nsez1
	 is=isez1(iic)
	 ip=is+1
	 if(myrank.eq.0) then
	    iicd=iic/10
	    iicu=mod(iic,10)
	    open(83,file='results1_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'_vort.q')
!        write(83,101)lamb2*(ny-2)/nyp,1,(nzg-2)/nzp
        write(83,101)lamb2*nypc,1,nzpc
	    write(83,102)ntime,i1,i1,i1
!!!!!!		  write(83,*)
	    write(83,103)
!     %(((i1,j=1,lamb2*(ny-2)/nyp),i=1,1),k=1,(nzg-2)/nzp)
     %(((i1,j=1,lamb2*nypc),i=1,1),k=1,nzpc)
	    do 40 k=kz1,kz2,nzp
	    do 41 l=1,lamb2
	    do 41 j=jy1,jy2,nyp
c           if(all(flaguo(is,j,k,:)).le.0) then
               write(83,9)vorr1(iic,j,k)
c           else
c              write(83,9)1.e+03
c           endif
   41	    continue
   40	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
c		  kk1=kz1+jjp*(nz-2)
c		  kk2=kz2+jjp*(nz-2)
c		  CALL MPI_RECV
c    %(flaguor(jy1:jy2,kk1:kk2,:),
c    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,jjp,1,
c    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 42 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 43 lr=1,lamb2
		  do 43 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,9)recvar(jr,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
   43		  continue
   42		  continue
	       enddo
	    endif
c
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
c           if(all(flaguo(is,j,k,:)).le.0) then
               write(83,9)voraz1(iic,j,k)
c           else
c              write(83,9)1.e+03
c           endif
   51	    continue
   50	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 52 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 53 lr=1,lamb2
		  do 53 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,9)recvar(jr,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
   53		  continue
   52		  continue
	       enddo
	    endif
c
	    do 60 k=kz1,kz2,nzp
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)vorz1(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
   61	    continue
   60	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 62 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 63 lr=1,lamb2
		  do 63 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,9)recvar(jr,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
   63		  continue
   62		  continue
	       enddo
	    endif
c
	    do 70 k=kz1,kz2,nzp
	    do 71 l=1,lamb2
	    do 71 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)0.5*(tv(is,j,k)+tv(ip,j,k))/ru1
c	    else
c	       write(83,9)1.e+03
c	    endif
   71	    continue
   70	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 72 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 73 lr=1,lamb2
		  do 73 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,9)recvar(jr,kr)/ru1
c		  else
c		     write(83,9)1.e+03
c		  endif
   73		  continue
   72		  continue
	       enddo
	    endif
	    close(83)
	 else
c	    CALL MPI_SEND
c    %(flaguo(is,jy1:jy2,kz1:kz2,:),
c    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,0,1,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorr1(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(voraz1(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorz1(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((tv(is,jy1:jy2,kz1:kz2)+tv(ip,jy1:jy2,kz1:kz2))*0.5,
     %1*(ny-2)*(nz-2),MTYPE,0,5,MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
      do kkc=1,nsez2max
	 ks=ksez2(kkc)
	 kp=ks+1
	 kk=isezindx2+(kkc-1)
	 kkd=kk/10
	 kku=mod(kk,10)
	 open(83,file='results2_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_vort.q')
!         write(83,101)lamb2*(ny-2)/nyp,(nx-2)/nxp,1
         write(83,101)lamb2*nypc,nxpc,1
	 write(83,102)ntime,i1,i1,i1
!!!!!!	       write(83,*)
	 write(83,103)
!     %(((i1,j=1,lamb2*(ny-2)/nyp),i=1,(nx-2)/nxp),k=1,1)
     %(((i1,j=1,lamb2*nypc),i=1,nxpc),k=1,1)
	 do 140 i=ix1,ix2,nxp
	 do 141 l=1,lamb2
	 do 141 j=jy1,jy2,nyp
c        if(all(flagwo(i,j,ks,:)).le.0) then
            write(83,9)vorr2(i,j,kkc)
c        else
c           write(83,9)1.e+03
c        endif
  141	 continue
  140	 continue
	 do 150 i=ix1,ix2,nxp
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
c        if(all(flagwo(i,j,ks,:)).le.0) then
            write(83,9)voraz2(i,j,kkc)
c        else
c           write(83,9)1.e+03
c        endif
  151	 continue
  150	 continue
	 do 160 i=ix1,ix2,nxp
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)vorz2(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
  161	 continue
  160	 continue
	 do 170 i=ix1,ix2,nxp
	 do 171 l=1,lamb2
	 do 171 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)0.5*(tv(i,j,ks)+tv(i,j,kp))/ru1
c	 else
c	    write(83,9)1.e+03
c	 endif
  171	 continue
  170	 continue
	 close(83)
      enddo
c
      do jjc=1,nsez14
	 js=jsez14(jjc)
	 jp=jpv(js)
	 if(myrank.eq.0) then
	    jjcd=jjc/10
	    jjcu=mod(jjc,10)
	    open(83,file='results3_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_vort.q')
!        write(83,101)1,(nx-2)/nxp,(nzg-2)/nzp
        write(83,101)1,nxpc,nzpc
	    write(83,102)ntime,i1,i1,i1
!!!!!!		  write(83,*)
	    write(83,103)
!     %(((i1,j=1,1),i=1,(nx-2)/nxp),k=1,(nzg-2)/nzp)
     %(((i1,j=1,1),i=1,nxpc),k=1,nzpc)
	    do 240 k=kz1,kz2,nzp
	    do 242 i=ix1,ix2,nxp
c           if(all(flagvo(i,js,k,:)).le.0) then
               write(83,9)vorr14(i,jjc,k)
c           else
c              write(83,9)1.e+03
c           endif
  242	    continue
  240	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
c		  kk1=kz1+jjp*(nz-2)
c		  kk2=kz2+jjp*(nz-2)
c		  CALL MPI_RECV
c    %(flagvor2(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,6,
c    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 243 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 245 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg,:)).le.0) then
		     write(83,9)recvar2(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
  245		  continue
  243		  continue
	       enddo
	    endif
c
	    do 250 k=kz1,kz2,nzp
	    do 252 i=ix1,ix2,nxp
c           if(all(flagvo(i,js,k,:)).le.0) then
               write(83,9)voraz14(i,jjc,k)
c           else
c              write(83,9)1.e+03
c           endif
  252	    continue
  250	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 253 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 255 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg,:)).le.0) then
		     write(83,9)recvar2(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
  255		  continue
  253		  continue
	       enddo
	    endif
c
	    do 260 k=kz1,kz2,nzp
	    do 262 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(83,9)vorz14(i,jjc,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
  262	    continue
  260	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 263 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 265 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg,:)).le.0) then
		     write(83,9)recvar2(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
  265		  continue
  263		  continue
	       enddo
	    endif
c
	    do 270 k=kz1,kz2,nzp
	    do 272 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(83,9)0.5*(tv(i,js,k)+tv(i,jp,k))/ru1
c	    else
c	       write(83,9)1.e+03
c	    endif
  272	    continue
  270	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 273 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 275 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg,:)).le.0) then
		     write(83,9)recvar2(ir,kr)/ru1
c		  else
c		     write(83,9)1.e+03
c		  endif
  275		  continue
  273		  continue
	       enddo
	    endif
	    close(83)
	 else
c	    CALL MPI_SEND
c    %(flagvo(ix1:ix2,js,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,6,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorr14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(voraz14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorz14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((tv(ix1:ix2,js,kz1:kz2)+tv(ix1:ix2,jp,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,14,MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
      return
      end
c
C---- subroutine allocate_output_xyz--------------------------------
C
      subroutine allocate_output_xyz(nx,ny,nz,icyl)
C
C-------------------------------------------------------------------
      use mediey
      use mediex
      use mediez
      use rmsy
      use rmsx
      use rmsz
      use vorticity
      use mediecalc
      use mediecalc1
      use mediecalc2
      use mediecalc3
      use vortavg
      use points
      implicit none
      include 'averages.h'

      integer nx,ny,nz,icyl

      ALLOCATE(voraz1(nsezmax19,ny,nz),vorr1(nsezmax19,ny,nz),
     %vorz1(nsezmax19,ny,nz))
      ALLOCATE(voraz2(nx,ny,nsez2max),vorr2(nx,ny,nsez2max),
     %vorz2(nx,ny,nsez2max))
      ALLOCATE(tempo1(nsez3,ny,nsez4max),uavg(nsez3,ny,nsez4max),
     %vavg(nsez3,ny,nsez4max),wavg(nsez3,ny,nsez4max),
     %vmodavg(nsez3,ny,nsez4max),vorazavg(nsez3,ny,nsez4max),
     %voraz(nsez3,ny,nsez4max),vorravg(nsez3,ny,nsez4max),
     %vorr(nsez3,ny,nsez4max),vorzavg(nsez3,ny,nsez4max),
     %vorz(nsez3,ny,nsez4max),pavg(nsez3,ny,nsez4max),
     %ptotavg(nsez3,ny,nsez4max),vrms(nsez3,ny,nsez4max),
     %urms(nsez3,ny,nsez4max),wrms(nsez3,ny,nsez4max),
     %prms(nsez3,ny,nsez4max),uvs(nsez3,ny,nsez4max),
     %uws(nsez3,ny,nsez4max),vws(nsez3,ny,nsez4max),
     %tvavg(nsez3,ny,nsez4max))
      ALLOCATE(tempo2(nx,nsez5,nsez6max),ruavg(nx,nsez5,nsez6max),
     %rvavg(nx,nsez5,nsez6max),rwavg(nx,nsez5,nsez6max),
     %rvmodavg(nx,nsez5,nsez6max),rvorazavg(nx,nsez5,nsez6max),
     %rvoraz(nx,nsez5,nsez6max),rvorravg(nx,nsez5,nsez6max),
     %rvorr(nx,nsez5,nsez6max),rvorzavg(nx,nsez5,nsez6max),
     %rvorz(nx,nsez5,nsez6max),rpavg(nx,nsez5,nsez6max),
     %rptotavg(nx,nsez5,nsez6max),rvrms(nx,nsez5,nsez6max),
     %rurms(nx,nsez5,nsez6max),rwrms(nx,nsez5,nsez6max),
     %rprms(nx,nsez5,nsez6max),ruvs(nx,nsez5,nsez6max),
     %ruws(nx,nsez5,nsez6max),rvws(nx,nsez5,nsez6max),
     %rtvavg(nx,nsez5,nsez6max))
      ALLOCATE(ztempo(nsez9,nsez91,nz),zuavg(nsez9,nsez91,nz),
     %zvavg(nsez9,nsez91,nz),zwavg(nsez9,nsez91,nz),
     %zvmodavg(nsez9,nsez91,nz),zvorazavg(nsez9,nsez91,nz),
     %zvorravg(nsez9,nsez91,nz),zvorzavg(nsez9,nsez91,nz),
     %zpavg(nsez9,nsez91,nz),zptotavg(nsez9,nsez91,nz),
     %zvrms(nsez9,nsez91,nz),zurms(nsez9,nsez91,nz),
     %zwrms(nsez9,nsez91,nz),zprms(nsez9,nsez91,nz),
     %zuvs(nsez9,nsez91,nz),zuws(nsez9,nsez91,nz),
     %zvws(nsez9,nsez91,nz),ztvavg(nsez9,nsez91,nz))
      ALLOCATE(tempo4(nsez1,ny,nz),tempo5(nx,ny,nsez2max))
      ALLOCATE(vplot(nsez1,ny,nz),uplot(nsez1,ny,nz),
     %wplot(nsez1,ny,nz),pplot(nsez1,ny,nz),vrmsplot(nsez1,ny,nz),
     %urmsplot(nsez1,ny,nz),wrmsplot(nsez1,ny,nz),
     %prmsplot(nsez1,ny,nz),uvplot(nsez1,ny,nz),vwplot(nsez1,ny,nz),
     %uwplot(nsez1,ny,nz),vmodplot(nsez1,ny,nz),tvplot(nsez1,ny,nz))
      ALLOCATE(vplot2(nx,ny,nsez2max),uplot2(nx,ny,nsez2max),
     %wplot2(nx,ny,nsez2max),pplot2(nx,ny,nsez2max),
     %vrmsplot2(nx,ny,nsez2max),urmsplot2(nx,ny,nsez2max),
     %wrmsplot2(nx,ny,nsez2max),prmsplot2(nx,ny,nsez2max),
     %uvplot2(nx,ny,nsez2max),vwplot2(nx,ny,nsez2max),
     %uwplot2(nx,ny,nsez2max),vmodplot2(nx,ny,nsez2max),
     %tvplot2(nx,ny,nsez2max))
      ALLOCATE(voraz1med(nsez1,ny,nz),vorr1med(nsez1,ny,nz),
     %vorz1med(nsez1,ny,nz),voraz2med(nx,ny,nsez2max),
     %vorr2med(nx,ny,nsez2max),vorz2med(nx,ny,nsez2max))
      ALLOCATE(timepoints1(nprbmax2))
      ALLOCATE(avrpoints1(nprbmax2,4),rmspoints1(nprbmax2,4))

      if(icyl.eq.0) then

        ALLOCATE(voraz14(nx,nsez14,nz),vorr14(nx,nsez14,nz),
     %vorz14(nx,nsez14,nz))
        ALLOCATE(tempo14(nx,nsez14,nz))
        ALLOCATE(vplot14(nx,nsez14,nz),uplot14(nx,nsez14,nz),
     %wplot14(nx,nsez14,nz),pplot14(nx,nsez14,nz),
     %vrmsplot14(nx,nsez14,nz),urmsplot14(nx,nsez14,nz),
     %wrmsplot14(nx,nsez14,nz),prmsplot14(nx,nsez14,nz),
     %uvplot14(nx,nsez14,nz),vwplot14(nx,nsez14,nz),
     %uwplot14(nx,nsez14,nz),vmodplot14(nx,nsez14,nz),
     %tvplot14(nx,nsez14,nz))
        ALLOCATE(voraz14med(nx,nsez14,nz),vorr14med(nx,nsez14,nz),
     %vorz14med(nx,nsez14,nz))

      else

        ALLOCATE(voraz14(nx,2*nsez14,nz),vorr14(nx,2*nsez14,nz),
     %vorz14(nx,2*nsez14,nz))
        ALLOCATE(tempo14(nx,2*nsez14,nz))
        ALLOCATE(vplot14(nx,2*nsez14,nz),uplot14(nx,2*nsez14,nz),
     %wplot14(nx,2*nsez14,nz),pplot14(nx,2*nsez14,nz),
     %vrmsplot14(nx,2*nsez14,nz),urmsplot14(nx,2*nsez14,nz),
     %wrmsplot14(nx,2*nsez14,nz),prmsplot14(nx,2*nsez14,nz),
     %uvplot14(nx,2*nsez14,nz),vwplot14(nx,2*nsez14,nz),
     %uwplot14(nx,2*nsez14,nz),vmodplot14(nx,2*nsez14,nz),
     %tvplot14(nx,2*nsez14,nz))
        ALLOCATE(voraz14med(nx,2*nsez14,nz),vorr14med(nx,2*nsez14,nz),
     %vorz14med(nx,2*nsez14,nz))

      endif

      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C---- subroutine setup_output_xyz ----------A. Posa - Jan 2013------
c
      subroutine setup_output_xyz(nx,ny,nz)
c
      use mediey
      use mediex
      use mediez
      use rmsy
      use rmsx
      use rmsz
      use mediecalc
      use mediecalc1
      use mediecalc2
      use mediecalc3
      use vortavg
      use points
      include 'common.h'
      include 'averages.h'
c
      integer iic,i,jjc,j,kkc,k,nx,ny,nz
c
      timep=0.
      iplant4=0
      iplant5=0
      iplant8=0
      iplant9=0
      tempo=0.
      do 10 iic=1,nsez3
      do 10 kkc=1,nsez4max
      do 10 j=1,ny
      tempo1(iic,j,kkc)=0.
      vmodavg(iic,j,kkc)=0.
      vavg(iic,j,kkc)=0.
      vorazavg(iic,j,kkc)=0.
      uavg(iic,j,kkc)=0.
      vorravg(iic,j,kkc)=0.
      wavg(iic,j,kkc)=0.
      vorzavg(iic,j,kkc)=0.
      pavg(iic,j,kkc)=0.
      ptotavg(iic,j,kkc)=0.
      vrms(iic,j,kkc)=0.
      urms(iic,j,kkc)=0.
      wrms(iic,j,kkc)=0.
      prms(iic,j,kkc)=0.
      uvs(iic,j,kkc)=0.
      uws(iic,j,kkc)=0.
      vws(iic,j,kkc)=0.
      tvavg(iic,j,kkc)=0.
  10  continue
      do 11 jjc=1,nsez5
      do 11 kkc=1,nsez6max
      do 11 i=1,nx
      tempo2(i,jjc,kkc)=0.
      rvmodavg(i,jjc,kkc)=0.
      rvavg(i,jjc,kkc)=0.
      rvorazavg(i,jjc,kkc)=0.
      ruavg(i,jjc,kkc)=0.
      rvorravg(i,jjc,kkc)=0.
      rwavg(i,jjc,kkc)=0.
      rvorzavg(i,jjc,kkc)=0.
      rpavg(i,jjc,kkc)=0.
      rptotavg(i,jjc,kkc)=0.
      rvrms(i,jjc,kkc)=0.
      rurms(i,jjc,kkc)=0.
      rwrms(i,jjc,kkc)=0.
      rprms(i,jjc,kkc)=0.
      ruvs(i,jjc,kkc)=0.
      ruws(i,jjc,kkc)=0.
      rvws(i,jjc,kkc)=0.
      rtvavg(i,jjc,kkc)=0.
  11  continue
      do 12 jjc=1,nsez91
      do 12 iic=1,nsez9
      do 12 k=1,nz
      zvmodavg(iic,jjc,k)=0.
      zvavg(iic,jjc,k)=0.
      zvorazavg(iic,jjc,k)=0.
      zuavg(iic,jjc,k)=0.
      zvorravg(iic,jjc,k)=0.
      zwavg(iic,jjc,k)=0.
      zvorzavg(iic,jjc,k)=0.
      zpavg(iic,jjc,k)=0.
      zptotavg(iic,jjc,k)=0.
      zvrms(iic,jjc,k)=0.
      zurms(iic,jjc,k)=0.
      zwrms(iic,jjc,k)=0.
      zprms(iic,jjc,k)=0.
      zuvs(iic,jjc,k)=0.
      zuws(iic,jjc,k)=0.
      zvws(iic,jjc,k)=0.
      ztvavg(iic,jjc,k)=0.
      ztempo(iic,jjc,k)=0.
  12  continue
      tempo3=0.
      do 13 k=1,nz
      do 13 j=1,ny
      do 13 iic=1,nsez1
      tempo4(iic,j,k)=0.
      vplot(iic,j,k)=0.
      uplot(iic,j,k)=0.
      wplot(iic,j,k)=0.
      pplot(iic,j,k)=0.
      vrmsplot(iic,j,k)=0.
      urmsplot(iic,j,k)=0.
      wrmsplot(iic,j,k)=0.
      prmsplot(iic,j,k)=0.
      uvplot(iic,j,k)=0.
      uwplot(iic,j,k)=0.
      vwplot(iic,j,k)=0.
      voraz1med(iic,j,k)=0.
      vorr1med(iic,j,k)=0.
      vorz1med(iic,j,k)=0.
      vmodplot(iic,j,k)=0.
      tvplot(iic,j,k)=0.
  13  continue
      do 14 kkc=1,nsez2max
      do 14 j=1,ny
      do 14 i=1,nx
      tempo5(i,j,kkc)=0.
      vplot2(i,j,kkc)=0.
      uplot2(i,j,kkc)=0.
      wplot2(i,j,kkc)=0.
      pplot2(i,j,kkc)=0.
      vrmsplot2(i,j,kkc)=0.
      urmsplot2(i,j,kkc)=0.
      wrmsplot2(i,j,kkc)=0.
      prmsplot2(i,j,kkc)=0.
      uvplot2(i,j,kkc)=0.
      uwplot2(i,j,kkc)=0.
      vwplot2(i,j,kkc)=0.
      voraz2med(i,j,kkc)=0.
      vorr2med(i,j,kkc)=0.
      vorz2med(i,j,kkc)=0.
      vmodplot2(i,j,kkc)=0.
      tvplot2(i,j,kkc)=0.
  14  continue
      do 15 k=1,nz
      do 15 jjc=1,(icyl+1)*nsez14
      do 15 i=1,nx
      tempo14(i,jjc,k)=0.
      vplot14(i,jjc,k)=0.
      uplot14(i,jjc,k)=0.
      wplot14(i,jjc,k)=0.
      pplot14(i,jjc,k)=0.
      vrmsplot14(i,jjc,k)=0.
      urmsplot14(i,jjc,k)=0.
      wrmsplot14(i,jjc,k)=0.
      prmsplot14(i,jjc,k)=0.
      uvplot14(i,jjc,k)=0.
      uwplot14(i,jjc,k)=0.
      vwplot14(i,jjc,k)=0.
      voraz14med(i,jjc,k)=0.
      vorr14med(i,jjc,k)=0.
      vorz14med(i,jjc,k)=0.
      vmodplot14(i,jjc,k)=0.
      tvplot14(i,jjc,k)=0.
  15  continue
      timepoints=0.
      do 16 i=1,nprbmax2
      timepoints1(i)=0.
      do 16 j=1,4
      avrpoints1(i,j)=0.
      rmspoints1(i,j)=0.
  16  continue
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine plotvortstagg_xyz_d(iplant9c,iplant9d,iplant9p,
     %nx,ny,nz,nzg,nbd,ntime,dens,uo,vo,wo,tv,xc,xu)
!     %,flaguo,flagvo,flagwo)
c
      use sections
      use vorticity
      include'common.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer iplant9c,iplant9d,iplant9p,nx,ny,nz,nzg,nbd,ntime
!      integer flaguo(nx,ny,nz,nbd),flagwo(nx,ny,nz,nbd),
!     %flagvo(nx,ny,nz,nbd)
      real dens(nx,ny,nz),uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),
     %tv(nx,ny,nz)
      real xc(nx),xu(nx)
c
c Local variables
      integer iic,iicd,iicu,is,ip,jjc,jjcd,jjcu,js,jp,
     %kkc,ks,kp,i,j,k,l,ir,jr,kr,lr,i1,jjp,kk1,kk2,
     %kk,kkd,kku,krg
      integer STATUS(MPI_STATUS_SIZE)
!      integer flaguor(ny,nzg,nbd),
!     %flagvor2(nx,nzg,nbd)
      real recvar(ny,nz),recvar2(nx,nz)
c
c Parameters
      integer lamb2,nxp,nyp,nzp,nxpc,nypc,nzpc
      parameter (lamb2=1)
      parameter (nyp=1,nxp=1,nzp=1)
c
   8  format(10f25.8)
   9  format(10f25.8)
 101  format(3i8)
 102  format(i8,3i2)
 103  format(i2)
c
      if(mod(nx-2,nxp).eq.0) then
        nxpc = (nx-2)/nxp
      else
        nxpc = (nx-2)/nxp+1
      endif
      if(mod(ny-2,nyp).eq.0) then
        nypc = (ny-2)/nyp
      else
        nypc = (ny-2)/nyp+1
      endif
      if(mod(nz-2,nzp).eq.0) then
        nzpc = (nzg-2)/nzp
      else
        nzpc = ((nz-2)/nzp+1)*mysize
      endif
c
      i1=1
      do iic=1,nsez1
         is=isez1(iic)
         ip=is+1
         if(myrank.eq.0) then
            iicd=iic/10
            iicu=mod(iic,10)
            open(83,file='results1_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'_vort.q')
!        write(83,101)lamb2*(ny-2)/nyp,1,(nzg-2)/nzp
        write(83,101)lamb2*nypc,1,nzpc
            write(83,102)ntime,i1,i1,i1
!!!!!!            write(83,*)
!            write(83,103)
!!     %(((i1,j=1,lamb2*(ny-2)/nyp),i=1,1),k=1,(nzg-2)/nzp)
!     %(((i1,j=1,lamb2*nypc),i=1,1),k=1,nzpc)
            do 30 k=kz1,kz2,nzp
            do 31 l=1,lamb2
            do 31 j=jy1,jy2,nyp
c           if(all(flaguo(is,j,k,:)).le.0) then
               write(83,9)0.5*(dens(is,j,k)+dens(ip,j,k))-denP1*xu(is)
c           else
c              write(83,9)1.e+03
c           endif
   31       continue
   30       continue
            if(mysize.ne.1) then
               do jjp=1,mysize-1
c                 kk1=kz1+jjp*(nz-2)
c                 kk2=kz2+jjp*(nz-2)
c                 CALL MPI_RECV
c    %(flaguor(jy1:jy2,kk1:kk2,:),
c    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,jjp,1,
c    %MPI_COMM_EDDY,STATUS,IERR)
                  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,21,
     %MPI_COMM_EDDY,STATUS,IERR)
                  do 32 kr=kz1,kz2,nzp
c                 krg=kr+jjp*(nz-2)
                  do 33 lr=1,lamb2
                  do 33 jr=jy1,jy2,nyp
c                 if(all(flaguor(jr,krg,:)).le.0) then
                     write(83,9)recvar(jr,kr)-denP1*xu(is)
c                 else
c                    write(83,9)1.e+03
c                 endif
   33             continue
   32             continue
               enddo
            endif
c
            do 40 k=kz1,kz2,nzp
            do 41 l=1,lamb2
            do 41 j=jy1,jy2,nyp
c           if(all(flaguo(is,j,k,:)).le.0) then
               write(83,9)vorr1(iic,j,k)
c           else
c              write(83,9)1.e+03
c           endif
   41       continue
   40       continue
            if(mysize.ne.1) then
               do jjp=1,mysize-1
c                 kk1=kz1+jjp*(nz-2)
c                 kk2=kz2+jjp*(nz-2)
                  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
                  do 42 kr=kz1,kz2,nzp
c                 krg=kr+jjp*(nz-2)
                  do 43 lr=1,lamb2
                  do 43 jr=jy1,jy2,nyp
c                 if(all(flaguor(jr,krg,:)).le.0) then
                     write(83,9)recvar(jr,kr)
c                 else
c                    write(83,9)1.e+03
c                 endif
   43             continue
   42             continue
               enddo
            endif
c
            do 50 k=kz1,kz2,nzp
            do 51 l=1,lamb2
            do 51 j=jy1,jy2,nyp
c           if(all(flaguo(is,j,k,:)).le.0) then
               write(83,9)voraz1(iic,j,k)
c           else
c              write(83,9)1.e+03
c           endif
   51       continue
   50       continue
            if(mysize.ne.1) then
               do jjp=1,mysize-1
                  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
                  do 52 kr=kz1,kz2,nzp
c                 krg=kr+jjp*(nz-2)
                  do 53 lr=1,lamb2
                  do 53 jr=jy1,jy2,nyp
c                 if(all(flaguor(jr,krg,:)).le.0) then
                     write(83,9)recvar(jr,kr)
c                 else
c                    write(83,9)1.e+03
c                 endif
   53             continue
   52             continue
               enddo
            endif
c
            do 60 k=kz1,kz2,nzp
            do 61 l=1,lamb2
            do 61 j=jy1,jy2,nyp
c           if(all(flaguo(is,j,k,:)).le.0) then
               write(83,9)vorz1(iic,j,k)
c           else
c              write(83,9)1.e+03
c           endif
   61       continue
   60       continue
            if(mysize.ne.1) then
               do jjp=1,mysize-1
                  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
                  do 62 kr=kz1,kz2,nzp
c                 krg=kr+jjp*(nz-2)
                  do 63 lr=1,lamb2
                  do 63 jr=jy1,jy2,nyp
c                 if(all(flaguor(jr,krg,:)).le.0) then
                     write(83,9)recvar(jr,kr)
c                 else
c                    write(83,9)1.e+03
c                 endif
   63             continue
   62             continue
               enddo
            endif
c
            do 70 k=kz1,kz2,nzp
            do 71 l=1,lamb2
            do 71 j=jy1,jy2,nyp
c           if(all(flaguo(is,j,k,:)).le.0) then
               write(83,9)0.5*(tv(is,j,k)+tv(ip,j,k))/ru1
c           else
c              write(83,9)1.e+03
c           endif
   71       continue
   70       continue
            if(mysize.ne.1) then
               do jjp=1,mysize-1
                  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
                  do 72 kr=kz1,kz2,nzp
c                 krg=kr+jjp*(nz-2)
                  do 73 lr=1,lamb2
                  do 73 jr=jy1,jy2,nyp
c                 if(all(flaguor(jr,krg,:)).le.0) then
                     write(83,9)recvar(jr,kr)/ru1
c                 else
c                    write(83,9)1.e+03
c                 endif
   73             continue
   72             continue
               enddo
            endif
            close(83)
         else
c           CALL MPI_SEND
c    %(flaguo(is,jy1:jy2,kz1:kz2,:),
c    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,0,1,
c    %MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_SEND
     %((dens(is,jy1:jy2,kz1:kz2)+dens(ip,jy1:jy2,kz1:kz2))*0.5,
     %1*(ny-2)*(nz-2),MTYPE,0,21,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_SEND
     %(vorr1(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_SEND
     %(voraz1(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_SEND
     %(vorz1(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_SEND
     %((tv(is,jy1:jy2,kz1:kz2)+tv(ip,jy1:jy2,kz1:kz2))*0.5,
     %1*(ny-2)*(nz-2),MTYPE,0,5,MPI_COMM_EDDY,STATUS,IERR)
         endif
      enddo
c
      do kkc=1,nsez2max
         ks=ksez2(kkc)
         kp=ks+1
         kk=isezindx2+(kkc-1)
         kkd=kk/10
         kku=mod(kk,10)
         open(83,file='results2_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_vort.q')
!         write(83,101)lamb2*(ny-2)/nyp,(nx-2)/nxp,1
         write(83,101)lamb2*nypc,nxpc,1
         write(83,102)ntime,i1,i1,i1
!!!!!!         write(83,*)
!         write(83,103)
!!     %(((i1,j=1,lamb2*(ny-2)/nyp),i=1,(nx-2)/nxp),k=1,1)
!     %(((i1,j=1,lamb2*nypc),i=1,nxpc),k=1,1)
         do 130 i=ix1,ix2,nxp
         do 131 l=1,lamb2
         do 131 j=jy1,jy2,nyp
c        if(all(flagwo(i,j,ks,:)).le.0) then
            write(83,9)0.5*(dens(i,j,ks)+dens(i,j,kp))-denP1*xc(i)
c        else
c           write(83,9)1.e+03
c        endif
  131    continue
  130    continue
         do 140 i=ix1,ix2,nxp
         do 141 l=1,lamb2
         do 141 j=jy1,jy2,nyp
c        if(all(flagwo(i,j,ks,:)).le.0) then
            write(83,9)vorr2(i,j,kkc)
c        else
c           write(83,9)1.e+03
c        endif
  141    continue
  140    continue
         do 150 i=ix1,ix2,nxp
         do 151 l=1,lamb2
         do 151 j=jy1,jy2,nyp
c        if(all(flagwo(i,j,ks,:)).le.0) then
            write(83,9)voraz2(i,j,kkc)
c        else
c           write(83,9)1.e+03
c        endif
  151    continue
  150    continue
         do 160 i=ix1,ix2,nxp
         do 161 l=1,lamb2
         do 161 j=jy1,jy2,nyp
c        if(all(flagwo(i,j,ks,:)).le.0) then
            write(83,9)vorz2(i,j,kkc)
c        else
c           write(83,9)1.e+03
c        endif
  161    continue
  160    continue
         do 170 i=ix1,ix2,nxp
         do 171 l=1,lamb2
         do 171 j=jy1,jy2,nyp
c        if(all(flagwo(i,j,ks,:)).le.0) then
            write(83,9)0.5*(tv(i,j,ks)+tv(i,j,kp))/ru1
c        else
c           write(83,9)1.e+03
c        endif
  171    continue
  170    continue
         close(83)
      enddo
c
      do jjc=1,nsez14
         js=jsez14(jjc)
         jp=jpv(js)
         if(myrank.eq.0) then
            jjcd=jjc/10
            jjcu=mod(jjc,10)
            open(83,file='results3_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_vort.q')
!        write(83,101)1,(nx-2)/nxp,(nzg-2)/nzp
        write(83,101)1,nxpc,nzpc
            write(83,102)ntime,i1,i1,i1
!!!!!!            write(83,*)
!            write(83,103)
!!     %(((i1,j=1,1),i=1,(nx-2)/nxp),k=1,(nzg-2)/nzp)
!     %(((i1,j=1,1),i=1,nxpc),k=1,nzpc)
            do 230 k=kz1,kz2,nzp
            do 232 i=ix1,ix2,nxp
c           if(all(flagvo(i,js,k,:)).le.0) then
               write(83,9)0.5*(dens(i,js,k)+dens(i,jp,k))-denP1*xc(i)
c           else
c              write(83,9)1.e+03
c           endif
  232       continue
  230       continue
            if(mysize.ne.1) then
               do jjp=1,mysize-1
c                 kk1=kz1+jjp*(nz-2)
c                 kk2=kz2+jjp*(nz-2)
c                 CALL MPI_RECV
c    %(flagvor2(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,6,
c    %MPI_COMM_EDDY,STATUS,IERR)
                  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,22,
     %MPI_COMM_EDDY,STATUS,IERR)
                  do 233 kr=kz1,kz2,nzp
c                 krg=kr+jjp*(nz-2)
                  do 235 ir=ix1,ix2,nxp
c                 if(all(flagvor2(ir,krg,:)).le.0) then
                     write(83,9)recvar2(ir,kr)-denP1*xc(ir)
c                 else
c                    write(83,9)1.e+03
c                 endif
  235             continue
  233             continue
               enddo
            endif
c
            do 240 k=kz1,kz2,nzp
            do 242 i=ix1,ix2,nxp
c           if(all(flagvo(i,js,k,:)).le.0) then
               write(83,9)vorr14(i,jjc,k)
c           else
c              write(83,9)1.e+03
c           endif
  242       continue
  240       continue
            if(mysize.ne.1) then
               do jjp=1,mysize-1
c                 kk1=kz1+jjp*(nz-2)
c                 kk2=kz2+jjp*(nz-2)
                  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
                  do 243 kr=kz1,kz2,nzp
c                 krg=kr+jjp*(nz-2)
                  do 245 ir=ix1,ix2,nxp
c                 if(all(flagvor2(ir,krg,:)).le.0) then
                     write(83,9)recvar2(ir,kr)
c                 else
c                    write(83,9)1.e+03
c                 endif
  245             continue
  243             continue
               enddo
            endif
c
            do 250 k=kz1,kz2,nzp
            do 252 i=ix1,ix2,nxp
c           if(all(flagvo(i,js,k,:)).le.0) then
               write(83,9)voraz14(i,jjc,k)
c           else
c              write(83,9)1.e+03
c           endif
  252       continue
  250       continue
            if(mysize.ne.1) then
               do jjp=1,mysize-1
                  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
                  do 253 kr=kz1,kz2,nzp
c                 krg=kr+jjp*(nz-2)
                  do 255 ir=ix1,ix2,nxp
c                 if(all(flagvor2(ir,krg,:)).le.0) then
                     write(83,9)recvar2(ir,kr)
c                 else
c                    write(83,9)1.e+03
c                 endif
  255             continue
  253             continue
               enddo
            endif
c
            do 260 k=kz1,kz2,nzp
            do 262 i=ix1,ix2,nxp
c           if(all(flagvo(i,js,k,:)).le.0) then
               write(83,9)vorz14(i,jjc,k)
c           else
c              write(83,9)1.e+03
c           endif
  262       continue
  260       continue
            if(mysize.ne.1) then
               do jjp=1,mysize-1
                  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
                  do 263 kr=kz1,kz2,nzp
c                 krg=kr+jjp*(nz-2)
                  do 265 ir=ix1,ix2,nxp
c                 if(all(flagvor2(ir,krg,:)).le.0) then
                     write(83,9)recvar2(ir,kr)
c                 else
c                    write(83,9)1.e+03
c                 endif
  265             continue
  263             continue
               enddo
            endif
c
            do 270 k=kz1,kz2,nzp
            do 272 i=ix1,ix2,nxp
c           if(all(flagvo(i,js,k,:)).le.0) then
               write(83,9)0.5*(tv(i,js,k)+tv(i,jp,k))/ru1
c           else
c              write(83,9)1.e+03
c           endif
  272       continue
  270       continue
            if(mysize.ne.1) then
               do jjp=1,mysize-1
                  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
                  do 273 kr=kz1,kz2,nzp
c                 krg=kr+jjp*(nz-2)
                  do 275 ir=ix1,ix2,nxp
c                 if(all(flagvor2(ir,krg,:)).le.0) then
                     write(83,9)recvar2(ir,kr)/ru1
c                 else
c                    write(83,9)1.e+03
c                 endif
  275             continue
  273             continue
               enddo
            endif
            close(83)
         else
c           CALL MPI_SEND
c    %(flagvo(ix1:ix2,js,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,6,
c    %MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_SEND
     %((dens(ix1:ix2,js,kz1:kz2)+dens(ix1:ix2,jp,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,22,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_SEND
     %(vorr14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_SEND
     %(voraz14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_SEND
     %(vorz14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_SEND
     %((tv(ix1:ix2,js,kz1:kz2)+tv(ix1:ix2,jp,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,14,MPI_COMM_EDDY,STATUS,IERR)
         endif
      enddo
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine plotinststagg_xyz_d(iplant9c,iplant9d,iplant9p,
     %nx,ny,nz,nzg,nbd,ntime,xc,xu,yc,yv,zcg,zwg,p,vo,uo,wo)
!     %,flaguo,flagvo,flagwo)
c
      use sections
      include'common.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer iplant9c,iplant9d,iplant9p,nx,ny,nz,nzg,nbd,ntime
!      integer flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),
!     %flagwo(nx,ny,nz,nbd)
      real xc(nx),xu(nx),yc(ny),yv(ny),zcg(nzg),zwg(nzg)
      real p(nx,ny,nz),vo(nx,ny,nz),uo(nx,ny,nz),wo(nx,ny,nz)
c
c Local variables
      integer iic,iicd,iicu,is,im,ip,kkc,kkcd,kkcu,ks,km,kp,
     %jjc,jjcd,jjcu,js,jm,jp,i,j,k,l,ir,jr,kr,lr,
     %i1,jjp,kk1,kk2,kk,kkd,kku
      integer STATUS(MPI_STATUS_SIZE)
!      integer flaguor(ny,nzg,nbd),
!     %flagvor2(nx,nzg,nbd)
      real pr(ny,nzg),vor(ny,nzg),uor(ny,nzg),wor(ny,nzg),der(ny,nzg),
     %pr2(nx,nzg),vor2(nx,nzg),uor2(nx,nzg),wor2(nx,nzg),der2(nx,nzg)
      integer kg
c
c Parameters
      integer lamb2,nxp,nyp,nzp,nxpc,nypc,nzpc
      real deltasect
      parameter (lamb2=1)
      parameter (nyp=1,nxp=1,nzp=1)
      parameter (deltasect=2.*pi)
c
   8  format(10f25.8)
   9  format(10f25.8)
 101  format(3i8)
 102  format(i8,3i2)
 103  format(i2)
c
      if(mod(nx-2,nxp).eq.0) then
        nxpc = (nx-2)/nxp
      else
        nxpc = (nx-2)/nxp+1
      endif
      if(mod(ny-2,nyp).eq.0) then
        nypc = (ny-2)/nyp
      else
        nypc = (ny-2)/nyp+1
      endif
      if(mod(nz-2,nzp).eq.0) then
        nzpc = (nzg-2)/nzp
      else
        nzpc = ((nz-2)/nzp+1)*mysize
      endif
c
      if(iplant9.eq.1) then
         if(myrank.eq.0) then
           do iic=1,nsez1
           iicd=iic/10
           iicu=mod(iic,10)
           open
     %(81,file='grid1_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'_plotinst.xyz')
           is=isez1(iic)
!           write(81,101)lamb2*(ny-2)/nyp,1,(nzg-2)/nzp
           write(81,101)lamb2*nypc,1,nzpc
!           do 10 k=2,nzg-1,nzp
           do 10 k=1,nzpc
           do 10 l=1,lamb2
           do 10 j=jy1,jy2,nyp
   10      write(81,9)xu(is)
!           do 20 k=2,nzg-1,nzp
           do 20 k=1,nzpc
           do 20 l=1,lamb2
           do 20 j=jy1,jy2,nyp
   20      write(81,9)yc(j)
!           do 30 k=2,nzg-1,nzp
           do 30 i=1,mysize
           do 30 k=2,nz-1,nzp
           kg=k+(i-1)*(nz-2)
           do 30 l=1,lamb2
           do 30 j=jy1,jy2,nyp
!   30      write(81,9)zcg(k)
   30      write(81,9)zcg(kg)
           close(81)
           enddo
           do kkc=1,nsez2
           kkcd=kkc/10
           kkcu=mod(kkc,10)
           open(81,file='grid2_2D_slice'//char(48+kkcd)
     %//char(48+kkcu)//'_plotinst.xyz')
           ks=ksez2g(kkc)
!           write(81,101)lamb2*(ny-2)/nyp,(nx-2)/nxp,1
           write(81,101)lamb2*nypc,nxpc,1
           do 110 i=ix1,ix2,nxp
           do 110 l=1,lamb2
           do 110 j=jy1,jy2,nyp
  110      write(81,9)xc(i)
           do 120 i=ix1,ix2,nxp
           do 120 l=1,lamb2
           do 120 j=jy1,jy2,nyp
  120      write(81,9)yc(j)
           do 130 i=ix1,ix2,nxp
           do 130 l=1,lamb2
           do 130 j=jy1,jy2,nyp
  130      write(81,9)zwg(ks)
           close(81)
           enddo
           do jjc=1,nsez14
           jjcd=jjc/10
           jjcu=mod(jjc,10)
           open(81,file='grid3_2D_slice'//char(48+jjcd)
     %//char(48+jjcu)//'_plotinst.xyz')
           js=jsez14(jjc)
!           write(81,101)1,(nx-2)/nxp,(nzg-2)/nzp
           write(81,101)1,nxpc,nzpc
!           do 210 k=2,nzg-1,nzp
           do 210 k=1,nzpc
           do 212 i=ix1,ix2,nxp
  212      write(81,9)xc(i)
  210      continue
!           do 220 k=2,nzg-1,nzp
           do 220 k=1,nzpc
           do 222 i=ix1,ix2,nxp
  222      write(81,9)yv(js)
  220      continue
!           do 230 k=2,nzg-1,nzp
           do 230 l=1,mysize
           do 230 k=2,nz-1,nzp
           kg=k+(l-1)*(nz-2)
           do 230 i=1,nx-2,nxp
!  230	       write(81,9)zcg(k)
  230      write(81,9)zcg(kg)
           close(81)
           enddo
        endif
      endif
c
      iplant9p=iplant9
      iplant9c=iplant9p/100
      iplant9p=mod(iplant9p,100)
      iplant9d=iplant9p/10
      iplant9p=mod(iplant9p,10)
c
      do iic=1,nsez1
	 is=isez1(iic)
	 ip=is+1
         im=is-1
	 if(myrank.eq.0) then
	    iicd=iic/10
	    iicu=mod(iic,10)
	    open(82,file='results1_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'_plotinst.q')
!        write(82,101)lamb2*(ny-2)/nyp,1,(nzg-2)/nzp
        write(82,101)lamb2*nypc,1,nzpc
	    i1=1
	    write(82,102)ntime,i1,i1,i1
!
! Cross-stream derivative of the cross-stream velocity field
!
	    do 40 k=kz1,kz2,nzp
	    do 41 l=1,lamb2
	    do 41 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)(uo(ip,j,k)-uo(im,j,k))/(xu(ip)-xu(im))
c	    else
c	       write(82,9)1.e+03
c	    endif
 41	    continue
 40	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
c		  CALL MPI_RECV
c    %(flaguor(jy1:jy2,kk1:kk2,:),
c    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,jjp,1,
c    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(pr(jy1:jy2,kk1:kk2),1*(ny-2)*(nz-2),MTYPE,jjp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(vor(jy1:jy2,kk1:kk2),1*(ny-2)*(nz-2),MTYPE,jjp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(uor(jy1:jy2,kk1:kk2),1*(ny-2)*(nz-2),MTYPE,jjp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(wor(jy1:jy2,kk1:kk2),1*(ny-2)*(nz-2),MTYPE,jjp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(der(jy1:jy2,kk1:kk2),1*(ny-2)*(nz-2),MTYPE,jjp,21,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 42 kr=kk1,kk2,nzp
		  do 43 lr=1,lamb2
		  do 43 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)der(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 43		  continue
 42		  continue
	       enddo
	    endif
!	    write(82,*)	  !!!!!!
!	    write(82,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp),i=1,1),k=1,(nzg-2)/nzp)
!
! X velocity
!
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
c           if(all(flaguo(is,j,k,:)).le.0) then
               write(82,9)uo(is,j,k)
c           else
c              write(82,9)1.e+03
c           endif
   51	    continue
   50	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 52 kr=kk1,kk2,nzp
		  do 53 lr=1,lamb2
		  do 53 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)uor(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   53		  continue
   52		  continue
	       enddo
	    endif
!
! Y velocity
!
	    do 60 k=kz1,kz2,nzp
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2,nyp
c           if(all(flaguo(is,j,k,:)).le.0) then
               jm=jmv(j)
               write(82,9)
     %0.25*(vo(is,j,k)+vo(ip,j,k)+vo(is,jm,k)+vo(ip,jm,k))
c           else
c              write(82,9)1.e+03
c           endif
   61	    continue
   60	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 62 kr=kk1,kk2,nzp
		  do 63 lr=1,lamb2
		  do 63 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)vor(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   63		  continue
   62		  continue
	       enddo
	    endif
!
!Axial velocity
!
	    do 70 k=kz1,kz2,nzp
	    km=k-1
	    do 71 l=1,lamb2
	    do 71 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)
     %0.25*(wo(is,j,k)+wo(ip,j,k)+wo(is,j,km)+wo(ip,j,km))
c	    else
c	       write(82,9)1.e+03
c	    endif
   71	    continue
   70	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 72 kr=kk1,kk2,nzp
		  do 73 lr=1,lamb2
		  do 73 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)wor(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   73		  continue
   72		  continue
	       enddo
	    endif
!
!Static pressure
!
	    do 80 k=kz1,kz2,nzp
	    do 81 l=1,lamb2
	    do 81 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,8)0.5*(p(is,j,k)+p(ip,j,k))   !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
   81	    continue
   80	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 82 kr=kk1,kk2,nzp
		  do 83 lr=1,lamb2
		  do 83 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,8)pr(jr,kr)   !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
   83		  continue
   82		  continue
	       enddo
	    endif
	    close(82)
	 else
c	    CALL MPI_SEND
c    %(flaguo(is,jy1:jy2,kz1:kz2,:),
c    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,0,1,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(0.5*(p(is,jy1:jy2,kz1:kz2)+p(ip,jy1:jy2,kz1:kz2)),
     %1*(ny-2)*(nz-2),MTYPE,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(0.25*(vo(is,1:jy2-1,kz1:kz2)+vo(ip,1:jy2-1,kz1:kz2)+
     %vo(is,jy1:jy2,kz1:kz2)+vo(ip,jy1:jy2,kz1:kz2)),
     %1*(ny-2)*(nz-2),MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uo(is,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(0.25*(wo(is,jy1:jy2,1:kz2-1)+wo(ip,jy1:jy2,1:kz2-1)+
     %wo(is,jy1:jy2,kz1:kz2)+wo(ip,jy1:jy2,kz1:kz2)),
     %1*(ny-2)*(nz-2),MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_SEND
     %((uo(ip,jy1:jy2,kz1:kz2)-uo(im,jy1:jy2,kz1:kz2))/(xu(ip)-xu(im)),
     %1*(ny-2)*(nz-2),MTYPE,0,21,MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo 
c
      do kkc=1,nsez2max
	 ks=ksez2(kkc) 
	 kp=ks+1
	 kk=isezindx2+(kkc-1)
	 kkd=kk/10
	 kku=mod(kk,10) 
	 open(82,file='results2_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_plotinst.q')
!         write(82,101)lamb2*(ny-2)/nyp,(nx-2)/nxp,1
         write(82,101)lamb2*nypc,nxpc,1
	 i1=1
	 write(82,102)ntime,i1,i1,i1
!
! Cross-stream derivative of the cross-stream velocity
!
	 do 140 i=ix1,ix2,nxp
	 im=i-1
	 do 141 l=1,lamb2
	 do 141 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)
     %0.5*(uo(i,j,ks)+uo(i,j,kp)-uo(im,j,ks)-uo(im,j,kp))/(xu(i)-xu(im))
c	 else
c	    write(82,9)1.e+03
c	 endif
 141	 continue
 140	 continue
!	 write(82,*)	!!!!!!
!	 write(82,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp),i=1,(nx-2)/nxp),k=1,1)
!
! X velocity
!
	 do 150 i=ix1,ix2,nxp
	 im=i-1
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
c        if(all(flagwo(i,j,ks,:)).le.0) then
            write(82,9)
     %0.25*(uo(i,j,ks)+uo(im,j,ks)+uo(i,j,kp)+uo(im,j,kp))
c        else
c           write(82,9)1.e+03
c        endif
 151	 continue
 150	 continue
!
! Y velocity
!
	 do 160 i=ix1,ix2,nxp
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2,nyp
c        if(all(flagwo(i,j,ks,:)).le.0) then
            jm=jmv(j)
            write(82,9)
     %0.25*(vo(i,j,ks)+vo(i,jm,ks)+vo(i,j,kp)+vo(i,jm,kp))
c        else
c           write(82,9)1.e+03
c        endif
 161	 continue
 160	 continue
!
!Axial velocity
!
	 do 170 i=ix1,ix2,nxp
	 do 171 l=1,lamb2
	 do 171 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)wo(i,j,ks)
c	 else
c	    write(82,9)1.e+03
c	 endif
 171	 continue
 170	 continue
!
!Static pressure
!
	 do 180 i=ix1,ix2,nxp
	 do 181 l=1,lamb2
	 do 181 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,8)0.5*(p(i,j,ks)+p(i,j,kp))    !*ros
c	 else
c	    write(82,8)1.e+03
c	 endif
 181	 continue
 180	 continue
	 close(82)
      enddo
c
      do jjc=1,nsez14
	 js=jsez14(jjc)
	 jp=jpv(js)
	 if(myrank.eq.0) then
	    jjcd=jjc/10
	    jjcu=mod(jjc,10)
	    open(82,file='results3_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_plotinst.q')
!        write(82,101)1,(nx-2)/nxp,(nzg-2)/nzp
        write(82,101)1,nxpc,nzpc
	    i1=1
	    write(82,102)ntime,i1,i1,i1
!
! Cross-stream derivative of the cross-stream velocity
!
	    do 240 k=kz1,kz2,nzp
	    do 242 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       im=i-1
	       write(82,9)
     %0.5*(uo(i,js,k)+uo(i,jp,k)-uo(im,js,k)-uo(im,jp,k))/(xu(i)-xu(im))
c	    else
c	       write(82,9)1.e+03
c	    endif
 242	    continue
 240	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
c		  CALL MPI_RECV
c    %(flagvor2(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,6,
c    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(pr2(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,7,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(vor2(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(uor2(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(wor2(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(der2(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,22,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 243 kr=kk1,kk2,nzp
		  do 245 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)der2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 245		  continue
 243		  continue
	       enddo
	    endif
!
! X velocity
!
	    do 250 k=kz1,kz2,nzp
	    do 252 i=ix1,ix2,nxp
c           if(all(flagvo(i,js,k,:)).le.0) then
               im=i-1
               write(82,9)
     %0.25*(uo(i,js,k)+uo(im,js,k)+uo(i,jp,k)+uo(im,jp,k))
c           else
c              write(82,9)1.e+03
c           endif
 252	    continue
 250	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 253 kr=kk1,kk2,nzp
		  do 255 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)uor2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 255		  continue
 253		  continue
	       enddo
	    endif
!
! Y velocity
!
	    do 260 k=kz1,kz2,nzp
	    do 262 i=ix1,ix2,nxp
c           if(all(flagvo(i,js,k,:)).le.0) then
               write(82,9)vo(i,js,k)
c           else
c              write(82,9)1.e+03
c           endif
 262	    continue
 260	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 263 kr=kk1,kk2,nzp
		  do 265 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)vor2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 265		  continue
 263		  continue
	       enddo
	    endif
!
!Axial velocity
!
	    do 270 k=kz1,kz2,nzp
	    km=k-1
	    do 272 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,9)
     %0.25*(wo(i,js,k)+wo(i,jp,k)+wo(i,js,km)+wo(i,jp,km))
c	    else
c	       write(82,9)1.e+03
c	    endif
 272	    continue
 270	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 273 kr=kk1,kk2,nzp
		  do 275 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)wor2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 275		  continue
 273		  continue
	       enddo
	    endif
!
!Static pressure
!
	    do 280 k=kz1,kz2,nzp
	    do 282 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,8)0.5*(p(i,js,k)+p(i,jp,k))    !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
 282	    continue
 280	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 283 kr=kk1,kk2,nzp
		  do 285 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,8)pr2(ir,kr)    !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
 285		  continue
 283		  continue
	       enddo
	    endif
	    close(82)
	 else
c	    CALL MPI_SEND
c    %(flagvo(ix1:ix2,js,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,6,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((p(ix1:ix2,js,kz1:kz2)+p(ix1:ix2,jp,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vo(ix1:ix2,js,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((uo(1:ix2-1,js,kz1:kz2)+uo(ix1:ix2,js,kz1:kz2)+
     %uo(1:ix2-1,jp,kz1:kz2)+uo(ix1:ix2,jp,kz1:kz2))*0.25,
     %(nx-2)*1*(nz-2),MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((wo(ix1:ix2,js,1:kz2-1)+wo(ix1:ix2,jp,1:kz2-1)+
     %wo(ix1:ix2,js,kz1:kz2)+wo(ix1:ix2,jp,kz1:kz2))*0.25,
     %(nx-2)*1*(nz-2),MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
           do k=kz1,kz2
             do i=ix1,ix2
               im=i-1
               der2(i,k)=
     %0.5*(-uo(im,js,k)+uo(i,js,k)-uo(im,jp,k)+uo(i,jp,k))/(xu(i)-xu(im))
             enddo
           enddo
           CALL MPI_SEND
     %(der2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,22,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
      return
      end
c
