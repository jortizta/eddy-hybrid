cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C---- subroutine setup_output --------------A. Posa - Oct 2011------
c
      subroutine setup_output(nx,ny,nz)
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
      do 15 jjc=1,2*nsez14
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C---- subroutine read_sections_plot2D -------A. Posa - Dec 2010------
C
C     PURPOSE: Read sections for plot2D and plot2Dvort.
C
C--------------------------------------------------------------------
      subroutine read_sections_plot2D(nz)
c
      use sections
      include 'common.h'
      include 'averages.h'
c
c.... Input/Output Arrays
      integer nz
c
c.... Local arrays
      integer i,ii,jj,k,kg,kk1,kk2

      open(unit=10,file='plot2D.input',form='formatted')
      read(10,*) (isez1(i),i=1,nsez1)
      if(MYRANK==0) write(*,*) (isez1(i),i=1,nsez1)
      read(10,*) (ksez2g(i),i=1,nsez2)
      if(MYRANK==0) write(*,*) (ksez2g(i),i=1,nsez2)
      read(10,*) (jsez14(i),i=1,nsez14)
      if(MYRANK==0) write(*,*) (jsez14(i),i=1,nsez14)
      close(10)

      call sortsections(ksez2g,nsez2)

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

      jj=0
      isezindx2=0
      do ii=1,nsez2
	kg=ksez2g(ii)
	k=kg-myrank*(nz-2)
	if(k>=kk1 .AND. k<=kk2) then
	  jj = jj+1
!	   ksez2(jj)=k
	  if(jj.eq.1)isezindx2=ii
	endif
      enddo
      nsez2max = jj
c
      allocate(ksez2(nsez2max))
c
      jj=0
      do ii=isezindx2,isezindx2+nsez2max-1
	jj=jj+1
	kg=ksez2g(ii)
	k=kg-myrank*(nz-2)
	ksez2(jj)=k
      enddo

      return

      end
C--------------------------------------------------------------------
c
C---- subroutine read_sections --------------A. Posa - Dec 2010------
C
C     PURPOSE: Read sections for mediestaggavg.
C
C--------------------------------------------------------------------
      subroutine read_sections(nz)
c
      use sections
      include 'common.h'
      include 'averages.h'
c
c.... Input/Output Arrays
      integer nz
c
c.... Local arrays
      integer i,ii,jj,k,kg,kk1,kk2
      integer ksez4g(nsez4),ksez6g(nsez6)
      logical file_exists

      inquire(file='medie.input',exist=file_exists)
      if(file_exists.eqv..true.) then
        open(unit=10,file='medie.input',form='formatted')
        read(10,*) (isez3(i),i=1,nsez3)
        read(10,*) (ksez4g(i),i=1,nsez4)
        read(10,*) (jsez5(i),i=1,nsez5)
        read(10,*) (ksez6g(i),i=1,nsez6)
        read(10,*) (isez9(i),i=1,nsez9)
        read(10,*) (jsez91(i),i=1,nsez91)
        close(10)
      endif

      call sortsections(ksez4g,nsez4)
      call sortsections(ksez6g,nsez6)

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

      jj=0
      isezindx4=0
      do ii=1,nsez4
	kg=ksez4g(ii)
	k=kg-myrank*(nz-2)
	if(k>=kk1 .AND. k<=kk2) then
	  jj = jj+1
!	   ksez4(jj)=k
	  if(jj.eq.1)isezindx4=ii
	endif
      enddo
      nsez4max = jj
c
      allocate(ksez4(nsez4max))
c
      jj=0
      do ii=isezindx4,isezindx4+nsez4max-1
	jj=jj+1
	kg=ksez4g(ii)
	k=kg-myrank*(nz-2)
	ksez4(jj)=k
      enddo

      jj=0
      isezindx6=0
      do ii=1,nsez6
	kg=ksez6g(ii)
	k=kg-myrank*(nz-2)
	if(k>=kk1 .AND. k<=kk2) then
	  jj = jj+1
!	   ksez6(jj)=k
	  if(jj.eq.1)isezindx6=ii
	endif
      enddo
      nsez6max = jj
c
      allocate(ksez6(nsez6max))
c
      jj=0
      do ii=isezindx6,isezindx6+nsez6max-1
	jj=jj+1
	kg=ksez6g(ii)
	k=kg-myrank*(nz-2)
	ksez6(jj)=k
      enddo

      return

      end
c
C--------------------------------------------------------------------
c
c---- subroutine sortsections ------------N. Beratlis-06 May 2010----
c--------------------------------------------A. Posa - Dec 2010------
C
C     PURPOSE: Sort sections in increasing z coordinate.
C
C--------------------------------------------------------------------
      subroutine sortsections(zp,n)
c
      implicit none
c
c...  Input/Output Arrays
      integer n
      integer zp(n)
c
c.... Local arrays
      integer i,imin
      integer temp
      integer z2(n)

      z2 = zp

      do i=1,n
	imin = i-1 + minloc(z2(i:n),1)
	temp = z2(i)
	z2(i) = z2(imin)
	z2(imin) = temp
      enddo

      zp = z2

      return

      end
C----------------------------------------------------------------------
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vort1staggavg(uo,vo,wo,nx,ny,nz,xc,xu,zc)
c
      use sections
      use mediey
      use mediex
      use vorticity
      include 'common.h'
      include 'averages.h'
c
c.... Input/Output array
      integer nx,ny,nz
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
      real xc(nx),xu(nx),zc(nz)
c
c.... Local arrays
      integer kkc,kc,kp,km,iic,ic,ip,im,jjc,jc,jp,jm,jsy,jpsy
      real dz,dtheta,dr,dq1x3,dq3x1,dq2x3,dq3x2,dq2x1,dq1x2
c
!     if(icyl.eq.1) then
c
c    ! Azimuthal vorticity - Distribution along y !
c
      do 101 kkc=1,nsez4max
	kc=ksez4(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
	do 101 iic=1,nsez3
	  ic=isez3(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 101 jc=jy1,jy2
	    dq1x3=(uo(ic,jc,kp)-uo(ic,jc,kc))/dz
	    dq3x1=(wo(ip,jc,kc)-wo(ic,jc,kc))/dr
	    voraz(iic,jc,kkc)=dq1x3-dq3x1
  101 continue
c
c    ! Azimuthal vorticity - Distribution along x (azimuthal average) !
c
      do 111 kkc=1,nsez6max
	kc=ksez6(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
	ic=ix1
	ip=ic+1
	im=ic-1
	dr=xc(ip)-xc(im)
	jjc=1
	rvoraz(ic,jjc,kkc)=0.
	do 112 jc=jy1,jy2
	   jp=jpv(jc)
	   jsy=jsym(jc)
	   jpsy=jsym(jp)
	   dq1x3=(0.25*(uo(ic,jc,kp)+uo(im,jc,kp)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,kc)+uo(im,jc,kc)
     %+uo(ic,jp,kc)+uo(im,jp,kc)))/dz
	   dq3x1=(0.5*(wo(ip,jc,kc)+wo(ip,jp,kc))
     %-0.5*(wo(ic,jsy,kc)+wo(ic,jpsy,kc)))/dr
	   rvoraz(ic,jjc,kkc)=
     %rvoraz(ic,jjc,kkc)+(dq1x3-dq3x1)
  112 continue
	do 113 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  jjc=1
	  rvoraz(ic,jjc,kkc)=0.
	  do 113 jc=jy1,jy2
	     jp=jpv(jc)
	     dq1x3=(0.25*(uo(ic,jc,kp)+uo(im,jc,kp)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,kc)+uo(im,jc,kc)
     %+uo(ic,jp,kc)+uo(im,jp,kc)))/dz
	     dq3x1=(0.5*(wo(ip,jc,kc)+wo(ip,jp,kc))
     %-0.5*(wo(im,jc,kc)+wo(im,jp,kc)))/dr
	     rvoraz(ic,jjc,kkc)=
     %rvoraz(ic,jjc,kkc)+(dq1x3-dq3x1)
  113 continue
  111 continue
c
c    ! Radial vorticity - Distribution along y !
c
      do 121 kkc=1,nsez4max
	kc=ksez4(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
	do 121 iic=1,nsez3
	  ic=isez3(iic)
	  ip=ic+1
	  do 121 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq3x2=(0.5*(wo(ip,jp,kc)+wo(ic,jp,kc))
     %-0.5*(wo(ip,jm,kc)+wo(ic,jm,kc)))/dtheta	    !
	    dq2x3=(0.25*(vo(ip,jc,kp)+vo(ic,jc,kp)
     %+vo(ip,jm,kp)+vo(ic,jm,kp))
     %-0.25*(vo(ip,jc,kc)+vo(ic,jc,kc)
     %+vo(ip,jm,kc)+vo(ic,jm,kc)))/dz
	    vorr(iic,jc,kkc)=dq3x2/xu(ic)-dq2x3
  121 continue
c
c    ! Radial vorticity - Distribution along x (azimuthal average) !
c
      do 131 kkc=1,nsez6max
	kc=ksez6(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
	do 131 ic=ix1,ix2
	  jjc=1
	  rvorr(ic,jjc,kkc)=0.
	  do 131 jc=jy1,jy2
	    jp=jpv(jc)
	    dtheta=dely
	    dq3x2=(wo(ic,jp,kc)-wo(ic,jc,kc))/dtheta	!
	    dq2x3=(vo(ic,jc,kp)-vo(ic,jc,kc))/dz
	    rvorr(ic,jjc,kkc)=
     %rvorr(ic,jjc,kkc)+(dq3x2/xc(ic)-dq2x3)
  131 continue
c
c    ! Axial vorticity - Distribution along y !
c
      do 141 kkc=1,nsez4max
	kc=ksez4(kkc)
	kp=kc+1
	do 141 iic=1,nsez3
	  ic=isez3(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 141 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.25*(vo(ip,jc,kp)+vo(ip,jm,kp)
     %+vo(ip,jc,kc)+vo(ip,jm,kc))*xc(ip)
     %-0.25*(vo(ic,jc,kp)+vo(ic,jm,kp)
     %+vo(ic,jc,kc)+vo(ic,jm,kc))*xc(ic))/dr
	    dq1x2=(0.5*(uo(ic,jp,kp)+uo(ic,jp,kc))
     %-0.5*(uo(ic,jm,kp)+uo(ic,jm,kc)))/dtheta	  !
	    vorz(iic,jc,kkc)=(dq2x1-dq1x2)/xu(ic)
  141 continue
c
c    ! Axial vorticity - Distribution along x (azimuthal average) !
c
      do 151 kkc=1,nsez6max
	kc=ksez6(kkc)
	kp=kc+1
	ic=ix1
	ip=ic+1
	im=ic-1
	dr=xu(ic)-xu(im)
	jjc=1
	rvorz(ic,jjc,kkc)=0.
	do 152 jc=jy1,jy2
	   jp=jpv(jc)
	   dtheta=dely
	   dq2x1=(0.25*(vo(ip,jc,kp)+vo(ip,jc,kc)
     %+vo(ic,jc,kp)+vo(ic,jc,kc))*xu(ic))/dr
	   dq1x2=(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,kc)+uo(im,jc,kc)
     %+uo(ic,jc,kp)+uo(im,jc,kp)))/dtheta    !
	   rvorz(ic,jjc,kkc)=
     %rvorz(ic,jjc,kkc)+(dq2x1-dq1x2)/xc(ic)
  152	continue
c
	do 153 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  jjc=1
	  rvorz(ic,jjc,kkc)=0.
	  do 153 jc=jy1,jy2
	    jp=jpv(jc)
	    dtheta=dely
	    dq2x1=(0.5*(vo(ip,jc,kc)+vo(ip,jc,kp))*xc(ip)
     %-0.5*(vo(im,jc,kc)+vo(im,jc,kp))*xc(im))/dr
	    dq1x2=
     %(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,kc)+uo(im,jc,kc)
     %+uo(ic,jc,kp)+uo(im,jc,kp)))/dtheta    !
	    rvorz(ic,jjc,kkc)=
     %rvorz(ic,jjc,kkc)+(dq2x1-dq1x2)/xc(ic)
  153 continue
  151 continue
c
c    ! Azimuthal vorticity - Distribution along z (azimuthal average) !
c
      do 161 kc=kz1,kz2
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	do 161 iic=1,nsez9
	  ic=isez9(iic)
	  ip=ic+1
	  jjc=1
	  voraz1(iic,jjc,kc)=0.
	  dr=xc(ip)-xc(ic)
	  do 161 jc=jy1,jy2
	    jp=jpv(jc)
	    dq1x3=(0.5*(uo(ic,jc,kp)+uo(ic,jp,kp))
     %-0.5*(uo(ic,jc,km)+uo(ic,jp,km)))/dz
	    dq3x1=(0.25*(wo(ip,jc,kc)+wo(ip,jp,kc)
     %+wo(ip,jc,km)+wo(ip,jp,km))
     %-0.25*(wo(ic,jc,kc)+wo(ic,jp,kc)
     %+wo(ic,jc,km)+wo(ic,jp,km)))/dr
	    voraz1(iic,jjc,kc)=
     %voraz1(iic,jjc,kc)+(dq1x3-dq3x1)
  161 continue
c
c    ! Radial vorticity - Distribution along z (azimuthal average) !
c
      do 171 kc=kz1,kz2
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	do 171 iic=1,nsez9
	  ic=isez9(iic)
	  ip=ic+1
	  jjc=1
	  vorr1(iic,jjc,kc)=0.
	  do 171 jc=jy1,jy2
	    jp=jpv(jc)
	    dtheta=dely
	    dq3x2=(0.25*(wo(ic,jp,kc)+wo(ip,jp,kc)
     %+wo(ic,jp,km)+wo(ip,jp,km))
     %-0.25*(wo(ic,jc,kc)+wo(ip,jc,kc)
     %+wo(ic,jc,km)+wo(ip,jc,km)))/dtheta   !
	    dq2x3=(0.5*(vo(ic,jc,kp)+vo(ip,jc,kp))
     %-0.5*(vo(ic,jc,km)+vo(ip,jc,km)))/dz
	    vorr1(iic,jjc,kc)=
     %vorr1(iic,jjc,kc)+(dq3x2/xu(ic)-dq2x3)
  171 continue
c
c    ! Axial vorticity - Distribution along z (azimuthal average) !
c
      do 181 kc=kz1,kz2
	do 181 iic=1,nsez9
	  ic=isez9(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  jjc=1
	  vorz1(iic,jjc,kc)=0.
	  do 181 jc=jy1,jy2
	    jp=jpv(jc)
	    dtheta=dely
	    dq2x1=(vo(ip,jc,kc)*xc(ip)
     %-vo(ic,jc,kc)*xc(ic))/dr
	    dq1x2=(uo(ic,jp,kc)
     %-uo(ic,jc,kc))/dtheta   !
	    vorz1(iic,jjc,kc)=
     %vorz1(iic,jjc,kc)+(dq2x1-dq1x2)/xu(ic)
  181 continue
c
!     endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine mediestaggavg(nx,ny,nz,uo,vo,wo,p,tv,dt)
c
      use sections
      use mediey
      use mediex
      use mediez
      use vorticity
      include'common.h'
      include'averages.h'
      include'mpif.h'
c
c Global variables
      integer nx,ny,nz
      real dt
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),
     %p(nx,ny,nz),tv(nx,ny,nz)
c
c Local variables
      integer iic,i,is,ip,im,jjc,j,jp,jm,jsy,jpsy,jmsy,
     %kkc,k,ks,kp,km,mm1,mm2,kg,nsez4maxr,nsez6maxr
      integer, dimension(:), allocatable :: ksez4r,ksez6r
      integer STATUS(MPI_STATUS_SIZE)
      real ustagg,vstagg,wstagg,pstagg,vmod,tvstagg
      real, dimension(:,:,:), allocatable :: vmodavgr,uavgr,vavgr,
     %wavgr,pavgr,ptotavgr,vorravgr,vorazavgr,vorzavgr,tvavgr,
     %rvmodavgr,ruavgr,rvavgr,rwavgr,rpavgr,rptotavgr,rvorravgr,
     %rvorazavgr,rvorzavgr,rtvavgr
      real zuavgr(nz),zvavgr(nz),zwavgr(nz),zpavgr(nz)
c
  10  format(4i6,10f25.8)
  11  format(3i6,10f25.8)
c
      tempo=tempo+dt
!
!     Time averages at some (x,z) positions along y
!
      do 101 iic=1,nsez3
      is=isez3(iic)
      ip=is+1
      if(is.eq.1) then
	do 1010 kkc=1,nsez4max
	ks=ksez4(kkc)
	kp=ks+1
	do 1011 j=jy1,jy2
	jm=jmv(j)
	jsy=jsym(j)
	jmsy=jsym(jm)
	tempo1(iic,j,kkc)=tempo1(iic,j,kkc)+dt
	vstagg=(vo(ip,j,ks)-vo(ip,jsy,ks)+
     %vo(ip,jm,ks)-vo(ip,jmsy,ks)+
     %vo(ip,j,kp)-vo(ip,jsy,kp)+
     %vo(ip,jm,kp)-vo(ip,jmsy,kp))*0.125
	ustagg=(uo(is,j,ks)+uo(is,j,kp))*0.5
	wstagg=(wo(ip,j,ks)+wo(ip,jsy,ks))*0.5
	pstagg=(p(ip,j,ks)+p(ip,j,kp)+
     %p(ip,jsy,ks)+p(ip,jsy,kp))*0.25
	tvstagg=(tv(ip,j,ks)+tv(ip,j,kp)+
     %tv(ip,jsy,ks)+tv(ip,jsy,kp))*0.25
	vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
	vavg(iic,j,kkc)=vavg(iic,j,kkc)
     %+vstagg*dt
	uavg(iic,j,kkc)=uavg(iic,j,kkc)
     %+ustagg*dt
	wavg(iic,j,kkc)=wavg(iic,j,kkc)
     %+wstagg*dt
	vmodavg(iic,j,kkc)=vmodavg(iic,j,kkc)
     %+vmod*dt
	vorazavg(iic,j,kkc)=vorazavg(iic,j,kkc)
     %+voraz(iic,j,kkc)*dt
	vorravg(iic,j,kkc)=vorravg(iic,j,kkc)
     %+vorr(iic,j,kkc)*dt
	vorzavg(iic,j,kkc)=vorzavg(iic,j,kkc)
     %+vorz(iic,j,kkc)*dt
	pavg(iic,j,kkc)=pavg(iic,j,kkc)
     %+pstagg*dt    !*ros
	ptotavg(iic,j,kkc)=ptotavg(iic,j,kkc)
     %+(pstagg+0.5*(vmod**2.))*dt   !*ros
	tvavg(iic,j,kkc)=tvavg(iic,j,kkc)
     %+tvstagg*dt
 1011	continue
 1010	continue
      else
	do 1012 kkc=1,nsez4max
	ks=ksez4(kkc)
	kp=ks+1
	do 1013 j=jy1,jy2
	jm=jmv(j)
	tempo1(iic,j,kkc)=tempo1(iic,j,kkc)+dt
	vstagg=(vo(is,j,ks)+vo(ip,j,ks)+
     %vo(is,jm,ks)+vo(ip,jm,ks)+
     %vo(is,j,kp)+vo(ip,j,kp)+
     %vo(is,jm,kp)+vo(ip,jm,kp))*0.125
	ustagg=(uo(is,j,ks)+uo(is,j,kp))*0.5
	wstagg=(wo(is,j,ks)+wo(ip,j,ks))*0.5
	pstagg=(p(is,j,ks)+p(ip,j,ks)+
     %p(is,j,kp)+p(ip,j,kp))*0.25
	tvstagg=(tv(is,j,ks)+tv(ip,j,ks)+
     %tv(is,j,kp)+tv(ip,j,kp))*0.25
	vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
	vavg(iic,j,kkc)=vavg(iic,j,kkc)
     %+vstagg*dt
	uavg(iic,j,kkc)=uavg(iic,j,kkc)
     %+ustagg*dt
	wavg(iic,j,kkc)=wavg(iic,j,kkc)
     %+wstagg*dt
	vmodavg(iic,j,kkc)=vmodavg(iic,j,kkc)
     %+vmod*dt
	vorazavg(iic,j,kkc)=vorazavg(iic,j,kkc)
     %+voraz(iic,j,kkc)*dt
	vorravg(iic,j,kkc)=vorravg(iic,j,kkc)
     %+vorr(iic,j,kkc)*dt
	vorzavg(iic,j,kkc)=vorzavg(iic,j,kkc)
     %+vorz(iic,j,kkc)*dt
	pavg(iic,j,kkc)=pavg(iic,j,kkc)
     %+pstagg*dt    !*ros
	ptotavg(iic,j,kkc)=ptotavg(iic,j,kkc)
     %+(pstagg+0.5*(vmod**2.))*dt   !*ros
	tvavg(iic,j,kkc)=tvavg(iic,j,kkc)
     %+tvstagg*dt
 1013	continue
 1012	continue
      endif
  101 continue
!
!     Time averages at some z positions along x (averaged along y)
!
      jjc=1
      do 104 kkc=1,nsez6max
      ks=ksez6(kkc)
      kp=ks+1
      do 103 i=ix1,ix2
      im=i-1
      tempo2(i,jjc,kkc)=tempo2(i,jjc,kkc)+dt
      do j=jy1,jy2
	 jp=jpv(j)
	 vstagg=(vo(i,j,ks)+vo(i,j,kp))*0.5
	 ustagg=(uo(i,j,ks)+uo(im,j,ks)+
     %uo(i,jp,ks)+uo(im,jp,ks)+
     %uo(i,j,kp)+uo(im,j,kp)+
     %uo(i,jp,kp)+uo(im,jp,kp))*0.125
	 wstagg=(wo(i,j,ks)+wo(i,jp,ks))*0.5
	 pstagg=(p(i,j,ks)+p(i,jp,ks)+
     %p(i,j,kp)+p(i,jp,kp))*0.25
	 tvstagg=(tv(i,j,ks)+tv(i,jp,ks)+
     %tv(i,j,kp)+tv(i,jp,kp))*0.25
	 rvavg(i,jjc,kkc)=rvavg(i,jjc,kkc)
     %+vstagg*dt
	 ruavg(i,jjc,kkc)=ruavg(i,jjc,kkc)
     %+ustagg*dt
	 rwavg(i,jjc,kkc)=rwavg(i,jjc,kkc)
     %+wstagg*dt
	 vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
	 rvmodavg(i,jjc,kkc)=rvmodavg(i,jjc,kkc)
     %+vmod*dt
	 rpavg(i,jjc,kkc)=rpavg(i,jjc,kkc)+
     %pstagg*dt	  !*ros
	 rptotavg(i,jjc,kkc)=rptotavg(i,jjc,kkc)
     %+(pstagg+0.5*(vmod**2.))*dt   !*ros
	 rtvavg(i,jjc,kkc)=rtvavg(i,jjc,kkc)+
     %tvstagg*dt
      enddo
      rvorazavg(i,jjc,kkc)=rvorazavg(i,jjc,kkc)
     %+rvoraz(i,jjc,kkc)*dt
      rvorravg(i,jjc,kkc)=rvorravg(i,jjc,kkc)
     %+rvorr(i,jjc,kkc)*dt
      rvorzavg(i,jjc,kkc)=rvorzavg(i,jjc,kkc)
     %+rvorz(i,jjc,kkc)*dt
  103 continue
  104 continue
!
!     Time averages at some x positions along z (averaged along y)
!
      jjc=1
      do 106 iic=1,nsez9
      is=isez9(iic)
      ip=is+1
      if(is.eq.1) then
	do 1060 k=kz1,kz2
	km=k-1
	do j=jy1,jy2
	  jp=jpv(j)
	  jsy=jsym(j)
	  jpsy=jsym(jp)
	  vstagg=(vo(ip,j,k)-vo(ip,jsy,k))*0.5
	  ustagg=(uo(is,j,k)+uo(is,jp,k))*0.5
	  wstagg=(wo(ip,j,k)+wo(ip,jsy,k)+
     %wo(ip,jp,k)+wo(ip,jpsy,k)+wo(ip,j,km)+wo(ip,jsy,km)+
     %wo(ip,jp,km)+wo(ip,jpsy,km))*0.125
	  pstagg=(p(ip,j,k)+p(ip,jsy,k)+
     %p(ip,jp,k)+p(ip,jpsy,k))*0.25
	  tvstagg=(tv(ip,j,k)+tv(ip,jsy,k)+
     %tv(ip,jp,k)+tv(ip,jpsy,k))*0.25
	  zvavg(iic,jjc,k)=zvavg(iic,jjc,k)
     %+vstagg*dt
	  zuavg(iic,jjc,k)=zuavg(iic,jjc,k)
     %+ustagg*dt
	  zwavg(iic,jjc,k)=zwavg(iic,jjc,k)
     %+wstagg*dt
	  vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
	  zvmodavg(iic,jjc,k)=
     %zvmodavg(iic,jjc,k)+vmod*dt
	  zpavg(iic,jjc,k)=zpavg(iic,jjc,k)+
     %pstagg*dt	  !*ros
	  zptotavg(iic,jjc,k)=zptotavg(iic,jjc,k)+
     %(pstagg+0.5*(vmod**2.))*dt   !*ros
	  ztvavg(iic,jjc,k)=ztvavg(iic,jjc,k)+
     %tvstagg*dt
	enddo
	zvorazavg(iic,jjc,k)=
     %zvorazavg(iic,jjc,k)+voraz1(iic,jjc,k)*dt
	zvorravg(iic,jjc,k)=
     %zvorravg(iic,jjc,k)+vorr1(iic,jjc,k)*dt
	zvorzavg(iic,jjc,k)=
     %zvorzavg(iic,jjc,k)+vorz1(iic,jjc,k)*dt
 1060	continue
      else
	do 1061 k=kz1,kz2
	km=k-1
	do j=jy1,jy2
	  jp=jpv(j)
	  vstagg=(vo(is,j,k)+vo(ip,j,k))*0.5
	  ustagg=(uo(is,j,k)+uo(is,jp,k))*0.5
	  wstagg=(wo(is,j,k)+wo(ip,j,k)+
     %wo(is,jp,k)+wo(ip,jp,k)+wo(is,j,km)+wo(ip,j,km)+
     %wo(is,jp,km)+wo(ip,jp,km))*0.125
	  pstagg=(p(is,j,k)+p(ip,j,k)+
     %p(is,jp,k)+p(ip,jp,k))*0.25
	  tvstagg=(tv(is,j,k)+tv(ip,j,k)+
     %tv(is,jp,k)+tv(ip,jp,k))*0.25
	  zvavg(iic,jjc,k)=zvavg(iic,jjc,k)
     %+vstagg*dt
	  zuavg(iic,jjc,k)=zuavg(iic,jjc,k)
     %+ustagg*dt
	  zwavg(iic,jjc,k)=zwavg(iic,jjc,k)
     %+wstagg*dt
	  vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
	  zvmodavg(iic,jjc,k)=
     %zvmodavg(iic,jjc,k)+vmod*dt
	  zpavg(iic,jjc,k)=zpavg(iic,jjc,k)+
     %pstagg*dt	  !*ros
	  zptotavg(iic,jjc,k)=zptotavg(iic,jjc,k)+
     %(pstagg+0.5*(vmod**2.))*dt   !*ros
	  ztvavg(iic,jjc,k)=ztvavg(iic,jjc,k)+
     %tvstagg*dt
	enddo
	zvorazavg(iic,jjc,k)=
     %zvorazavg(iic,jjc,k)+voraz1(iic,jjc,k)*dt
	zvorravg(iic,jjc,k)=
     %zvorravg(iic,jjc,k)+vorr1(iic,jjc,k)*dt
	zvorzavg(iic,jjc,k)=
     %zvorzavg(iic,jjc,k)+vorz1(iic,jjc,k)*dt
 1061	continue
      endif
  106 continue
      if(tempo.ge.periodo4) then
	 iplant8=iplant8+1
	 do 200 iic=1,nsez3
	 do 200 kkc=1,nsez4max
	 do 200 j=jy1,jy2
	 vmodavg(iic,j,kkc)=vmodavg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	 vavg(iic,j,kkc)=vavg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	 vorazavg(iic,j,kkc)=vorazavg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	 uavg(iic,j,kkc)=uavg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	 vorravg(iic,j,kkc)=vorravg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	 wavg(iic,j,kkc)=wavg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	 vorzavg(iic,j,kkc)=vorzavg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	 pavg(iic,j,kkc)=pavg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	 ptotavg(iic,j,kkc)=ptotavg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	 tvavg(iic,j,kkc)=tvavg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	 tvavg(iic,j,kkc)=tvavg(iic,j,kkc)/ru1
  200	 continue
	 if(myrank.eq.0) then
	    mm1=121
	    mm2=122
	    do 201 kkc=1,nsez4max
	    do 201 iic=1,nsez3
	    write(mm1,333)'i=',isez3(iic),'____k=',ksez4(kkc)
	    write(mm2,333)'i=',isez3(iic),'____k=',ksez4(kkc)
	    do 201 j=jy1,jy2
	    write(mm1,10)
     %iplant8,isez3(iic),j,ksez4(kkc),vmodavg(iic,j,kkc),
     %uavg(iic,j,kkc),vavg(iic,j,kkc),wavg(iic,j,kkc),
     %pavg(iic,j,kkc),ptotavg(iic,j,kkc),tvavg(iic,j,kkc)
	    write(mm2,10)
     %iplant8,isez3(iic),j,ksez4(kkc),vorravg(iic,j,kkc),
     %vorazavg(iic,j,kkc),vorzavg(iic,j,kkc)
  201	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(nsez4maxr,1,MPI_INTEGER,jp,1,
     %MPI_COMM_EDDY,STATUS,IERR)
		  ALLOCATE(ksez4r(nsez4maxr))
		  ALLOCATE(vmodavgr(nsez3,ny,nsez4maxr),
     %uavgr(nsez3,ny,nsez4maxr),vavgr(nsez3,ny,nsez4maxr),
     %wavgr(nsez3,ny,nsez4maxr),pavgr(nsez3,ny,nsez4maxr),
     %ptotavgr(nsez3,ny,nsez4maxr),vorravgr(nsez3,ny,nsez4maxr),
     %vorazavgr(nsez3,ny,nsez4maxr),vorzavgr(nsez3,ny,nsez4maxr),
     %tvavgr(nsez3,ny,nsez4maxr))
		  CALL MPI_RECV(ksez4r,nsez4maxr,MPI_INTEGER,jp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vmodavgr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(uavgr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vavgr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(wavgr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,6,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(pavgr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,7,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(ptotavgr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vorravgr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vorazavgr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vorzavgr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(tvavgr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 203 kkc=1,nsez4maxr
		  k=ksez4r(kkc)
		  kg=jp*(nz-2)+k
		  do 203 iic=1,nsez3
		  write(mm1,333)'i=',isez3(iic),'____k=',kg
		  write(mm2,333)'i=',isez3(iic),'____k=',kg
		  do 203 j=jy1,jy2
		  write(mm1,10)
     %iplant8,isez3(iic),j,kg,vmodavgr(iic,j,kkc),
     %uavgr(iic,j,kkc),vavgr(iic,j,kkc),wavgr(iic,j,kkc),
     %pavgr(iic,j,kkc),ptotavgr(iic,j,kkc),tvavgr(iic,j,kkc)
		  write(mm2,10)
     %iplant8,isez3(iic),j,kg,vorravgr(iic,j,kkc),
     %vorazavgr(iic,j,kkc),vorzavgr(iic,j,kkc)
  203		  continue
		  DEALLOCATE(ksez4r)
		  DEALLOCATE(vmodavgr,uavgr,vavgr,wavgr,pavgr,
     %ptotavgr,vorravgr,vorazavgr,vorzavgr,tvavgr)
	       enddo
	    endif
	 else
	    CALL MPI_SEND(nsez4max,1,MPI_INTEGER,0,1,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ksez4,nsez4max,MPI_INTEGER,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vmodavg(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(uavg(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vavg(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(wavg(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(pavg(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ptotavg(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vorravg(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vorazavg(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vorzavg(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(tvavg(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
c
	 jjc=1
	 do 204 kkc=1,nsez6max
	 do 204 i=ix1,ix2
	 rvmodavg(i,jjc,kkc)=rvmodavg(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2))
	 rvavg(i,jjc,kkc)=rvavg(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2))
	 rvorazavg(i,jjc,kkc)=rvorazavg(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2))
	 ruavg(i,jjc,kkc)=ruavg(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2))
	 rvorravg(i,jjc,kkc)=rvorravg(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2))
	 rwavg(i,jjc,kkc)=rwavg(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2))
	 rvorzavg(i,jjc,kkc)=rvorzavg(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2))
	 rpavg(i,jjc,kkc)=rpavg(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2))
	 rptotavg(i,jjc,kkc)=rptotavg(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2))
	 rtvavg(i,jjc,kkc)=rtvavg(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2))
	 rtvavg(i,jjc,kkc)=rtvavg(i,jjc,kkc)/ru1
 204	 continue
	 if(myrank.eq.0) then
	    mm1=123
	    mm2=124
	    do 205 kkc=1,nsez6max
	    write(mm1,334)'averages_along_y____k=',ksez6(kkc)
	    write(mm2,334)'averages_along_y____k=',ksez6(kkc)
	    do 205 i=ix1,ix2
	    write(mm1,11)
     %iplant8,i,ksez6(kkc),rvmodavg(i,jjc,kkc),
     %ruavg(i,jjc,kkc),rvavg(i,jjc,kkc),rwavg(i,jjc,kkc),
     %rpavg(i,jjc,kkc),rptotavg(i,jjc,kkc),rtvavg(i,jjc,kkc)
	    write(mm2,11)
     %iplant8,i,ksez6(kkc),rvorravg(i,jjc,kkc),
     %rvorazavg(i,jjc,kkc),rvorzavg(i,jjc,kkc)
 205	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(nsez6maxr,1,MPI_INTEGER,jp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
		  ALLOCATE(ksez6r(nsez6maxr))
		  ALLOCATE(rvmodavgr(nx,nsez5,nsez6maxr),
     %ruavgr(nx,nsez5,nsez6maxr),rvavgr(nx,nsez5,nsez6maxr),
     %rwavgr(nx,nsez5,nsez6maxr),rpavgr(nx,nsez5,nsez6maxr),
     %rptotavgr(nx,nsez5,nsez6maxr),rvorravgr(nx,nsez5,nsez6maxr),
     %rvorazavgr(nx,nsez5,nsez6maxr),rvorzavgr(nx,nsez5,nsez6maxr),
     %rtvavgr(nx,nsez5,nsez6maxr))
		  CALL MPI_RECV(ksez6r,nsez6maxr,MPI_INTEGER,jp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvmodavgr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(ruavgr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,16,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvavgr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,17,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rwavgr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,18,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rpavgr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,19,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rptotavgr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,20,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvorravgr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,21,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvorazavgr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,22,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvorzavgr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,23,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rtvavgr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,24,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 207 kkc=1,nsez6maxr
		  k=ksez6r(kkc)
		  kg=jp*(nz-2)+k
		  write(mm1,334)'averages_along_y____k=',kg
		  write(mm2,334)'averages_along_y____k=',kg
		  do 207 i=ix1,ix2
		  write(mm1,11)
     %iplant8,i,kg,rvmodavgr(i,jjc,kkc),
     %ruavgr(i,jjc,kkc),rvavgr(i,jjc,kkc),rwavgr(i,jjc,kkc),
     %rpavgr(i,jjc,kkc),rptotavgr(i,jjc,kkc),rtvavgr(i,jjc,kkc)
		  write(mm2,11)
     %iplant8,i,kg,rvorravgr(i,jjc,kkc),
     %rvorazavgr(i,jjc,kkc),rvorzavgr(i,jjc,kkc)
 207		  continue
		  DEALLOCATE(ksez6r)
		  DEALLOCATE(rvmodavgr,ruavgr,rvavgr,rwavgr,rpavgr,
     %rptotavgr,rvorravgr,rvorazavgr,rvorzavgr,rtvavgr)
	       enddo
	    endif
	 else
	    CALL MPI_SEND(nsez6max,1,MPI_INTEGER,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ksez6,nsez6max,MPI_INTEGER,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvmodavg(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ruavg(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,16,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvavg(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,17,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rwavg(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,18,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rpavg(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,19,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rptotavg(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,20,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvorravg(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,21,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvorazavg(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,22,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvorzavg(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,23,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rtvavg(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,24,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
c
	 jjc=1
	 do 208 iic=1,nsez9
	 do 208 k=kz1,kz2
	 zvmodavg(iic,jjc,k)=
     %zvmodavg(iic,jjc,k)/(tempo*real(ny-2))
	 zvavg(iic,jjc,k)=
     %zvavg(iic,jjc,k)/(tempo*real(ny-2))
	 zvorazavg(iic,jjc,k)=
     %zvorazavg(iic,jjc,k)/(tempo*real(ny-2))
	 zuavg(iic,jjc,k)=
     %zuavg(iic,jjc,k)/(tempo*real(ny-2))
	 zvorravg(iic,jjc,k)=
     %zvorravg(iic,jjc,k)/(tempo*real(ny-2))
	 zwavg(iic,jjc,k)=
     %zwavg(iic,jjc,k)/(tempo*real(ny-2))
	 zvorzavg(iic,jjc,k)=
     %zvorzavg(iic,jjc,k)/(tempo*real(ny-2))
	 zpavg(iic,jjc,k)=
     %zpavg(iic,jjc,k)/(tempo*real(ny-2))
	 zptotavg(iic,jjc,k)=
     %zptotavg(iic,jjc,k)/(tempo*real(ny-2))
	 ztvavg(iic,jjc,k)=
     %ztvavg(iic,jjc,k)/(tempo*real(ny-2))
	 ztvavg(iic,jjc,k)=
     %ztvavg(iic,jjc,k)/ru1
 208	 continue
	 mm1=125
	 mm2=126
	 do 209 iic=1,nsez9
	 if(myrank.eq.0) then
	    write(mm1,334)'averages_along_y____i=',isez9(iic)
	    write(mm2,334)'averages_along_y____i=',isez9(iic)
	    do 210 k=kz1,kz2
	    write(mm1,11)iplant8,isez9(iic),k,
     %zvmodavg(iic,jjc,k),zuavg(iic,jjc,k),zvavg(iic,jjc,k),
     %zwavg(iic,jjc,k),zpavg(iic,jjc,k),zptotavg(iic,jjc,k),
     %ztvavg(iic,jjc,k)
	    write(mm2,11)iplant8,isez9(iic),k,
     %zvorravg(iic,jjc,k),zvorazavg(iic,jjc,k),zvorzavg(iic,jjc,k)
 210	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(zvmodavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,25,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zuavgr(kz1:kz2),(nz-2),MTYPE,jp,26,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvavgr(kz1:kz2),(nz-2),MTYPE,jp,27,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zwavgr(kz1:kz2),(nz-2),MTYPE,jp,28,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zpavgr(kz1:kz2),(nz-2),MTYPE,jp,29,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zptotavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,30,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvorravg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,31,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvorazavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,32,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvorzavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,33,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(ztvavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,34,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 211 k=kz1,kz2
		  kg=jp*(nz-2)+k
		  write(mm1,11)iplant8,isez9(iic),kg,
     %zvmodavg(iic,jjc,k),zuavgr(k),zvavgr(k),zwavgr(k),
     %zpavgr(k),zptotavg(iic,jjc,k),ztvavg(iic,jjc,k)
		  write(mm2,11)iplant8,isez9(iic),kg,
     %zvorravg(iic,jjc,k),zvorazavg(iic,jjc,k),zvorzavg(iic,jjc,k)
 211		  continue
	       enddo
	    endif
	 else
	    CALL MPI_SEND(zvmodavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,25,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zuavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,26,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,27,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zwavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,28,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zpavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,29,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zptotavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,30,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvorravg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,31,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvorazavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,32,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvorzavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,33,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ztvavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,34,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
 209	 continue
      endif
      if(irms.eq.1)call rmsmodstaggavg(nx,ny,nz,uo,vo,wo,p,dt)
      if(tempo.ge.periodo4) then
	 do 110 iic=1,nsez3
	 do 110 kkc=1,nsez4max
	 do 110 j=jy1,jy2
	 vmodavg(iic,j,kkc)=0.
	 vavg(iic,j,kkc)=0.
	 vorazavg(iic,j,kkc)=0.
	 uavg(iic,j,kkc)=0.
	 vorravg(iic,j,kkc)=0.
	 wavg(iic,j,kkc)=0.
	 vorzavg(iic,j,kkc)=0.
	 pavg(iic,j,kkc)=0.
	 ptotavg(iic,j,kkc)=0.
	 tvavg(iic,j,kkc)=0.
	 tempo1(iic,j,kkc)=0.
  110	 continue
	 jjc=1
	 do 112 kkc=1,nsez6max
	 do 112 i=ix1,ix2
	 rvmodavg(i,jjc,kkc)=0.
	 rvavg(i,jjc,kkc)=0.
	 rvorazavg(i,jjc,kkc)=0.
	 ruavg(i,jjc,kkc)=0.
	 rvorravg(i,jjc,kkc)=0.
	 rwavg(i,jjc,kkc)=0.
	 rvorzavg(i,jjc,kkc)=0.
	 rpavg(i,jjc,kkc)=0.
	 rptotavg(i,jjc,kkc)=0.
	 rtvavg(i,jjc,kkc)=0.
	 tempo2(i,jjc,kkc)=0.
  112	 continue
	 jjc=1
	 do 114 iic=1,nsez9
	 do 114 k=kz1,kz2
	 zvmodavg(iic,jjc,k)=0.
	 zvavg(iic,jjc,k)=0.
	 zvorazavg(iic,jjc,k)=0.
	 zuavg(iic,jjc,k)=0.
	 zvorravg(iic,jjc,k)=0.
	 zwavg(iic,jjc,k)=0.
	 zvorzavg(iic,jjc,k)=0.
	 zpavg(iic,jjc,k)=0.
	 zptotavg(iic,jjc,k)=0.
	 ztvavg(iic,jjc,k)=0.
  114	 continue
 333	 format(a2,i6,a6,i6)
 334	 format(a25,i6)
	 tempo=0.
      endif
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine rmsmodstaggavg(nx,ny,nz,uo,vo,wo,p,dt)
c
      use sections
      use mediey
      use mediex
      use mediez
      use rmsy
      use rmsx
      use rmsz
      include'common.h'
      include'averages.h'
      include'mpif.h'
c
c Global variables
      integer nx,ny,nz
      real dt
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
c
c Local variables
      integer iic,i,is,ip,im,jjc,j,jp,jm,jsy,jpsy,jmsy,
     %kkc,k,ks,kp,km,mm1,mm2,kg,nsez4maxr,nsez6maxr
      integer, dimension(:), allocatable :: ksez4r,ksez6r
      integer STATUS(MPI_STATUS_SIZE)
      real ustagg,vstagg,wstagg,pstagg
      real, dimension(:,:,:), allocatable :: urmsr,vrmsr,wrmsr,
     %prmsr,uvsr,uwsr,vwsr,rurmsr,rvrmsr,rwrmsr,rprmsr,ruvsr,
     %ruwsr,rvwsr
c
  10  format(4i6,10f25.8)
  11  format(3i6,10f25.8)
c
      do 101 iic=1,nsez3
      is=isez3(iic)
      ip=is+1
      if(is.eq.1) then
	do 1010 kkc=1,nsez4max
	ks=ksez4(kkc)
	kp=ks+1
	do 1011 j=jy1,jy2
	jm=jmv(j)
	jsy=jsym(j)
	jmsy=jsym(jm)
	vstagg=0.125*(vo(ip,j,ks)-vo(ip,jsy,ks)+
     %vo(ip,jm,ks)-vo(ip,jmsy,ks)+
     %vo(ip,j,kp)-vo(ip,jsy,kp)+
     %vo(ip,jm,kp)-vo(ip,jmsy,kp))
	ustagg=0.5*(uo(is,j,ks)+uo(is,j,kp))
	wstagg=0.5*(wo(ip,j,ks)+wo(ip,jsy,ks))
	pstagg=0.25*(p(ip,j,ks)+p(ip,jsy,ks)+
     %p(ip,j,kp)+p(ip,jsy,kp))
	vrms(iic,j,kkc)=
     %vrms(iic,j,kkc)+dt*(vstagg**2.)
	urms(iic,j,kkc)=
     %urms(iic,j,kkc)+dt*(ustagg**2.)
	wrms(iic,j,kkc)=
     %wrms(iic,j,kkc)+dt*(wstagg**2.)
	prms(iic,j,kkc)=
     %prms(iic,j,kkc)+dt*(pstagg**2.)	!*ros**2.
	uvs(iic,j,kkc)=
     %uvs(iic,j,kkc)+dt*ustagg*vstagg
	uws(iic,j,kkc)=
     %uws(iic,j,kkc)+dt*ustagg*wstagg
	vws(iic,j,kkc)=
     %vws(iic,j,kkc)+dt*vstagg*wstagg
 1011	continue
 1010	continue
      else
	do 1012 kkc=1,nsez4max
	ks=ksez4(kkc)
	kp=ks+1
	do 1013 j=jy1,jy2
	jm=jmv(j)
	vstagg=0.125*(vo(is,j,ks)+vo(ip,j,ks)+
     %vo(is,jm,ks)+vo(ip,jm,ks)+
     %vo(is,j,kp)+vo(ip,j,kp)+
     %vo(is,jm,kp)+vo(ip,jm,kp))
	ustagg=0.5*(uo(is,j,ks)+uo(is,j,kp))
	wstagg=0.5*(wo(is,j,ks)+wo(ip,j,ks))
	pstagg=0.25*(p(is,j,ks)+p(ip,j,ks)+
     %p(is,j,kp)+p(ip,j,kp))
	vrms(iic,j,kkc)=
     %vrms(iic,j,kkc)+dt*(vstagg**2.)
	urms(iic,j,kkc)=
     %urms(iic,j,kkc)+dt*(ustagg**2.)
	wrms(iic,j,kkc)=
     %wrms(iic,j,kkc)+dt*(wstagg**2.)
	prms(iic,j,kkc)=
     %prms(iic,j,kkc)+dt*(pstagg**2.)	!*ros**2.
	uvs(iic,j,kkc)=
     %uvs(iic,j,kkc)+dt*ustagg*vstagg
	uws(iic,j,kkc)=
     %uws(iic,j,kkc)+dt*ustagg*wstagg
	vws(iic,j,kkc)=
     %vws(iic,j,kkc)+dt*vstagg*wstagg
 1013	continue
 1012	continue
      endif
  101 continue
c
      jjc=1
      do 104 kkc=1,nsez6max
      ks=ksez6(kkc)
      kp=ks+1
      do 103 i=ix1,ix2
      im=i-1
      do 103 j=jy1,jy2
	 jp=jpv(j)
	 vstagg=0.5*(vo(i,j,ks)+vo(i,j,kp))
	 ustagg=0.125*(uo(i,j,ks)+uo(im,j,ks)+
     %uo(i,jp,ks)+uo(im,jp,ks)+
     %uo(i,j,kp)+uo(im,j,kp)+
     %uo(i,jp,kp)+uo(im,jp,kp))
	 wstagg=0.5*(wo(i,j,ks)+wo(i,jp,ks))
	 pstagg=0.25*(p(i,j,ks)+p(i,jp,ks)+
     %p(i,j,kp)+p(i,jp,kp))
	 rvrms(i,jjc,kkc)=
     %rvrms(i,jjc,kkc)+dt*(vstagg**2.)
	 rurms(i,jjc,kkc)=
     %rurms(i,jjc,kkc)+dt*(ustagg**2.)
	 rwrms(i,jjc,kkc)=
     %rwrms(i,jjc,kkc)+dt*(wstagg**2.)
	 rprms(i,jjc,kkc)=
     %rprms(i,jjc,kkc)+dt*(pstagg**2.)	 !*ros**2.
	 ruvs(i,jjc,kkc)=
     %ruvs(i,jjc,kkc)+dt*ustagg*vstagg
	 ruws(i,jjc,kkc)=
     %ruws(i,jjc,kkc)+dt*ustagg*wstagg
	 rvws(i,jjc,kkc)=
     %rvws(i,jjc,kkc)+dt*vstagg*wstagg
  103 continue
  104 continue
c
      jjc=1
      do 106 iic=1,nsez9
      is=isez9(iic)
      ip=is+1
      if(is.eq.1) then
	do 1060 k=kz1,kz2
	km=k-1
	do 1061 j=jy1,jy2
	jp=jpv(j)
	jsy=jsym(j)
	jpsy=jsym(jp)
	vstagg=0.5*(vo(ip,j,k)-vo(ip,jsy,k))
	ustagg=0.5*(uo(is,j,k)+uo(is,jp,k))
	wstagg=0.125*(wo(ip,j,k)+wo(ip,jsy,k)+
     %wo(ip,jp,k)+wo(ip,jpsy,k)+
     %wo(ip,j,km)+wo(ip,jsy,km)+
     %wo(ip,jp,km)+wo(ip,jpsy,km))
	pstagg=0.25*(p(ip,j,k)+p(ip,jsy,k)+
     %p(ip,jp,k)+p(ip,jpsy,k))
	zvrms(iic,jjc,k)=zvrms(iic,jjc,k)
     %+dt*(vstagg**2.)
	zurms(iic,jjc,k)=zurms(iic,jjc,k)
     %+dt*(ustagg**2.)
	zwrms(iic,jjc,k)=zwrms(iic,jjc,k)
     %+dt*(wstagg**2.)
	zprms(iic,jjc,k)=zprms(iic,jjc,k)
     %+dt*(pstagg**2.)	 !*ros**2.
	zuvs(iic,jjc,k)=zuvs(iic,jjc,k)
     %+dt*ustagg*vstagg
	zuws(iic,jjc,k)=zuws(iic,jjc,k)
     %+dt*ustagg*wstagg
	zvws(iic,jjc,k)=zvws(iic,jjc,k)
     %+dt*vstagg*wstagg
 1061	continue
 1060	continue
      else
	do 1062 k=kz1,kz2
	km=k-1
	do 1063 j=jy1,jy2
	jp=jpv(j)
	vstagg=0.5*(vo(is,j,k)+vo(ip,j,k))
	ustagg=0.5*(uo(is,j,k)+uo(is,jp,k))
	wstagg=0.125*(wo(is,j,k)+wo(ip,j,k)+
     %wo(is,jp,k)+wo(ip,jp,k)+
     %wo(is,j,km)+wo(ip,j,km)+
     %wo(is,jp,km)+wo(ip,jp,km))
	pstagg=0.25*(p(is,j,k)+p(ip,j,k)+
     %p(is,jp,k)+p(ip,jp,k))
	zvrms(iic,jjc,k)=zvrms(iic,jjc,k)
     %+dt*(vstagg**2.)
	zurms(iic,jjc,k)=zurms(iic,jjc,k)
     %+dt*(ustagg**2.)
	zwrms(iic,jjc,k)=zwrms(iic,jjc,k)
     %+dt*(wstagg**2.)
	zprms(iic,jjc,k)=zprms(iic,jjc,k)
     %+dt*(pstagg**2.)	  !*ros**2.
	zuvs(iic,jjc,k)=zuvs(iic,jjc,k)
     %+dt*ustagg*vstagg
	zuws(iic,jjc,k)=zuws(iic,jjc,k)
     %+dt*ustagg*wstagg
	zvws(iic,jjc,k)=zvws(iic,jjc,k)
     %+dt*vstagg*wstagg
 1063	continue
 1062	continue
      endif
  106 continue
c
      if(tempo.ge.periodo4) then
	 do 121 iic=1,nsez3
	 do 121 kkc=1,nsez4max
	 do 120 j=jy1,jy2
	 vrms(iic,j,kkc)=
     %sqrt(abs((vrms(iic,j,kkc)/tempo1(iic,j,kkc))
     %-vavg(iic,j,kkc)**2.))
	 urms(iic,j,kkc)=
     %sqrt(abs((urms(iic,j,kkc)/tempo1(iic,j,kkc))
     %-uavg(iic,j,kkc)**2.))
	 wrms(iic,j,kkc)=
     %sqrt(abs((wrms(iic,j,kkc)/tempo1(iic,j,kkc))
     %-wavg(iic,j,kkc)**2.))
	 prms(iic,j,kkc)=
     %sqrt(abs((prms(iic,j,kkc)/tempo1(iic,j,kkc))
     %-pavg(iic,j,kkc)**2.))
	 uvs(iic,j,kkc)=
     %uvs(iic,j,kkc)/tempo1(iic,j,kkc)
     %-uavg(iic,j,kkc)*vavg(iic,j,kkc)
	 uws(iic,j,kkc)=
     %uws(iic,j,kkc)/tempo1(iic,j,kkc)
     %-uavg(iic,j,kkc)*wavg(iic,j,kkc)
	 vws(iic,j,kkc)=
     %vws(iic,j,kkc)/tempo1(iic,j,kkc)
     %-vavg(iic,j,kkc)*wavg(iic,j,kkc)
  120	 continue
  121	 continue
	 if(myrank.eq.0) then
	    mm1=127
	    mm2=128
	    do 122 kkc=1,nsez4max
	    do 122 iic=1,nsez3
	    write(mm1,333)'i=',isez3(iic),'____k=',ksez4(kkc)
	    write(mm2,333)'i=',isez3(iic),'____k=',ksez4(kkc)
	    do 122 j=jy1,jy2
	    write(mm1,10)iplant8,isez3(iic),j,ksez4(kkc),
     %urms(iic,j,kkc),vrms(iic,j,kkc),
     %wrms(iic,j,kkc),prms(iic,j,kkc)
	    write(mm2,10)iplant8,isez3(iic),j,ksez4(kkc),
     %uvs(iic,j,kkc),uws(iic,j,kkc),vws(iic,j,kkc)
  122	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(nsez4maxr,1,MPI_INTEGER,jp,1,
     %MPI_COMM_EDDY,STATUS,IERR)
		  ALLOCATE(ksez4r(nsez4maxr))
		  ALLOCATE(urmsr(nsez3,ny,nsez4maxr),
     %vrmsr(nsez3,ny,nsez4maxr),wrmsr(nsez3,ny,nsez4maxr),
     %prmsr(nsez3,ny,nsez4maxr),uvsr(nsez3,ny,nsez4maxr),
     %uwsr(nsez3,ny,nsez4maxr),vwsr(nsez3,ny,nsez4maxr))
		  CALL MPI_RECV(ksez4r,nsez4maxr,MPI_INTEGER,jp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(urmsr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vrmsr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(wrmsr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(prmsr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,6,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(uvsr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,7,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(uwsr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vwsr(:,jy1:jy2,:),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 123 kkc=1,nsez4maxr
		  k=ksez4r(kkc)
		  kg=jp*(nz-2)+k
		  do 123 iic=1,nsez3
		  write(mm1,333)'i=',isez3(iic),'____k=',kg
		  write(mm2,333)'i=',isez3(iic),'____k=',kg
		  do 123 j=jy1,jy2
		  write(mm1,10)
     %iplant8,isez3(iic),j,kg,
     %urmsr(iic,j,kkc),vrmsr(iic,j,kkc),
     %wrmsr(iic,j,kkc),prmsr(iic,j,kkc)
		  write(mm2,10)
     %iplant8,isez3(iic),j,kg,
     %uvsr(iic,j,kkc),uwsr(iic,j,kkc),vwsr(iic,j,kkc)
  123		  continue
		  DEALLOCATE(ksez4r)
		  DEALLOCATE(urmsr,vrmsr,wrmsr,prmsr,uvsr,uwsr,vwsr)
	       enddo
	    endif
	 else
	    CALL MPI_SEND(nsez4max,1,MPI_INTEGER,0,1,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ksez4,nsez4max,MPI_INTEGER,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(urms(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vrms(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(wrms(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(prms(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(uvs(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(uws(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vws(:,jy1:jy2,:),nsez3*(ny-2)*nsez4max,MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
c
	 do kkc=1,nsez4max
	    do iic=1,nsez3
	       do j=jy1,jy2
		  vrms(iic,j,kkc)=0.
		  urms(iic,j,kkc)=0.
		  wrms(iic,j,kkc)=0.
		  prms(iic,j,kkc)=0.
		  uvs(iic,j,kkc)=0.
		  uws(iic,j,kkc)=0.
		  vws(iic,j,kkc)=0.
	       enddo
	    enddo
	 enddo
c
	 jjc=1
	 do 124 kkc=1,nsez6max
	 do 124 i=ix1,ix2
	 rvrms(i,jjc,kkc)=
     %sqrt(abs((rvrms(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2)))
     %-rvavg(i,jjc,kkc)**2.))
	 rurms(i,jjc,kkc)=
     %sqrt(abs((rurms(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2)))
     %-ruavg(i,jjc,kkc)**2.))
	 rwrms(i,jjc,kkc)=
     %sqrt(abs((rwrms(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2)))
     %-rwavg(i,jjc,kkc)**2.))
	 rprms(i,jjc,kkc)=
     %sqrt(abs((rprms(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2)))
     %-rpavg(i,jjc,kkc)**2.))
	 ruvs(i,jjc,kkc)=
     %ruvs(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2))
     %-ruavg(i,jjc,kkc)*rvavg(i,jjc,kkc)
	 ruws(i,jjc,kkc)=
     %ruws(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2))
     %-ruavg(i,jjc,kkc)*rwavg(i,jjc,kkc)
	 rvws(i,jjc,kkc)=
     %rvws(i,jjc,kkc)/
     %(tempo2(i,jjc,kkc)*real(ny-2))
     %-rvavg(i,jjc,kkc)*rwavg(i,jjc,kkc)
  124	 continue
c
	 jjc=1
	 if(myrank.eq.0) then
	    mm1=129
	    mm2=130
	    do 125 kkc=1,nsez6max
	    write(mm1,334)'averages_along_y____k=',ksez6(kkc)
	    write(mm2,334)'averages_along_y____k=',ksez6(kkc)
	    do 125 i=ix1,ix2
	    write(mm1,11)iplant8,i,ksez6(kkc),
     %rurms(i,jjc,kkc),rvrms(i,jjc,kkc),
     %rwrms(i,jjc,kkc),rprms(i,jjc,kkc)
	    write(mm2,11)iplant8,i,ksez6(kkc),
     %ruvs(i,jjc,kkc),ruws(i,jjc,kkc),rvws(i,jjc,kkc)
  125	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(nsez6maxr,1,MPI_INTEGER,jp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  ALLOCATE(ksez6r(nsez6maxr))
		  ALLOCATE(rurmsr(nx,nsez5,nsez6maxr),
     %rvrmsr(nx,nsez5,nsez6maxr),rwrmsr(nx,nsez5,nsez6maxr),
     %rprmsr(nx,nsez5,nsez6maxr),ruvsr(nx,nsez5,nsez6maxr),
     %ruwsr(nx,nsez5,nsez6maxr),rvwsr(nx,nsez5,nsez6maxr))
		  CALL MPI_RECV(ksez6r,nsez6maxr,MPI_INTEGER,jp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rurmsr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvrmsr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rwrmsr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rprmsr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(ruvsr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,16,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(ruwsr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,17,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvwsr(ix1:ix2,:,:),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,18,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 126 kkc=1,nsez6maxr
		  k=ksez6r(kkc)
		  kg=jp*(nz-2)+k
		  write(mm1,334)'averages_along_y____k=',kg
		  write(mm2,334)'averages_along_y____k=',kg
		  do 126 i=ix1,ix2
		  write(mm1,11)
     %iplant8,i,kg,
     %rurmsr(i,jjc,kkc),rvrmsr(i,jjc,kkc),
     %rwrmsr(i,jjc,kkc),rprmsr(i,jjc,kkc)
		  write(mm2,11)
     %iplant8,i,kg,
     %ruvsr(i,jjc,kkc),ruwsr(i,jjc,kkc),rvwsr(i,jjc,kkc)
  126		  continue
		  DEALLOCATE(ksez6r)
		  DEALLOCATE(rurmsr,rvrmsr,rwrmsr,rprmsr,
     %ruvsr,ruwsr,rvwsr)
	       enddo
	    endif
	 else
	    CALL MPI_SEND(nsez6max,1,MPI_INTEGER,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ksez6,nsez6max,MPI_INTEGER,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rurms(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvrms(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rwrms(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rprms(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ruvs(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,16,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ruws(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,17,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvws(ix1:ix2,:,:),(nx-2)*nsez5*nsez6max,MTYPE,0,18,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
c
	 jjc=1
	 do kkc=1,nsez6max
	    do i=ix1,ix2
	       rvrms(i,jjc,kkc)=0.
	       rurms(i,jjc,kkc)=0.
	       rwrms(i,jjc,kkc)=0.
	       rprms(i,jjc,kkc)=0.
	       ruvs(i,jjc,kkc)=0.
	       ruws(i,jjc,kkc)=0.
	       rvws(i,jjc,kkc)=0.
	    enddo
	 enddo
c
	 jjc=1
	 do 127 iic=1,nsez9
	 do 127 k=kz1,kz2
	 zvrms(iic,jjc,k)=
     %sqrt(abs(zvrms(iic,jjc,k)/(tempo*real(ny-2))
     %-zvavg(iic,jjc,k)**2.))
	 zurms(iic,jjc,k)=
     %sqrt(abs(zurms(iic,jjc,k)/(tempo*real(ny-2))
     %-zuavg(iic,jjc,k)**2.))
	 zwrms(iic,jjc,k)=
     %sqrt(abs(zwrms(iic,jjc,k)/(tempo*real(ny-2))
     %-zwavg(iic,jjc,k)**2.))
	 zprms(iic,jjc,k)=
     %sqrt(abs(zprms(iic,jjc,k)/(tempo*real(ny-2))
     %-zpavg(iic,jjc,k)**2.))
	 zuvs(iic,jjc,k)=
     %zuvs(iic,jjc,k)/(tempo*real(ny-2))
     %-zuavg(iic,jjc,k)*zvavg(iic,jjc,k)
	 zuws(iic,jjc,k)=
     %zuws(iic,jjc,k)/(tempo*real(ny-2))
     %-zuavg(iic,jjc,k)*zwavg(iic,jjc,k)
	 zvws(iic,jjc,k)=
     %zvws(iic,jjc,k)/(tempo*real(ny-2))
     %-zvavg(iic,jjc,k)*zwavg(iic,jjc,k)
  127	 continue
c
	 mm1=131
	 mm2=132
	 do 128 iic=1,nsez9
	 if(myrank.eq.0) then
	    write(mm1,334)'averages_along_y____i=',isez9(iic)
	    write(mm2,334)'averages_along_y____i=',isez9(iic)
	    do 129 k=kz1,kz2
	    write(mm1,11)iplant8,isez9(iic),k,
     %zurms(iic,jjc,k),zvrms(iic,jjc,k),zwrms(iic,jjc,k),
     %zprms(iic,jjc,k)
	    write(mm2,11)iplant8,isez9(iic),k,
     %zuvs(iic,jjc,k),zuws(iic,jjc,k),zvws(iic,jjc,k)
  129	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(zurms(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,19,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvrms(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,20,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zwrms(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,21,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zprms(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,22,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zuvs(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,23,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zuws(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,24,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvws(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,25,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 130 k=kz1,kz2
		  kg=jp*(nz-2)+k
		  write(mm1,11)iplant8,isez9(iic),kg,
     %zurms(iic,jjc,k),zvrms(iic,jjc,k),zwrms(iic,jjc,k),
     %zprms(iic,jjc,k)
		  write(mm2,11)iplant8,isez9(iic),kg,
     %zuvs(iic,jjc,k),zuws(iic,jjc,k),zvws(iic,jjc,k)
 130		  continue
	       enddo
	    endif
	 else
	    CALL MPI_SEND(zurms(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,19,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvrms(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,20,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zwrms(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,21,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zprms(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,22,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zuvs(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,23,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zuws(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,24,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvws(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,25,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
 128	 continue
c
	 jjc=1
	 do iic=1,nsez9
	    do k=kz1,kz2
	       zvrms(iic,jjc,k)=0.
	       zurms(iic,jjc,k)=0.
	       zwrms(iic,jjc,k)=0.
	       zprms(iic,jjc,k)=0.
	       zuvs(iic,jjc,k)=0.
	       zuws(iic,jjc,k)=0.
	       zvws(iic,jjc,k)=0.
	    enddo
	 enddo
      endif
 333  format(a2,i6,a6,i6)
 334  format(a25,i6)
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vort2(nx,ny,nz,xc,xu,zc,uo,vo,wo)
c
      use sections
      use vorticity
      include'common.h'
      include'averages.h'
c
c Global variables
      integer nx,ny,nz
      real xc(nx),xu(nx),zc(nz)
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
c
c Local variables
      integer iic,ic,ip,im,jjc,jc,jp,jm,jsy,jpsy,jmsy,kkc,kc,kp,km
      real dtheta,dr,dz,dq1x3,dq3x1,dq2x3,dq3x2,dq1x2,dq2x1
c
      ! Azimuthal vorticity - Circumferential sections !
      do 101 kc=kz1,kz2
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	do 101 iic=1,nsez1
	  ic=isez1(iic)
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 101 jc=jy1,jy2
	    dq1x3=(0.5*(uo(ic,jc,kp)+uo(im,jc,kp))
     %-0.5*(uo(ic,jc,km)+uo(im,jc,km)))/dz
	    dq3x1=(0.5*(wo(ip,jc,kc)+wo(ip,jc,km))
     %-0.5*(wo(im,jc,kc)+wo(im,jc,km)))/dr
	    voraz1(iic,jc,kc)=dq1x3-dq3x1
  101 continue
c
      ! Azimuthal vorticity - Cross sections !
      do 111 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	ic=ix1
	ip=ic+1
	im=ic-1
	dr=xc(ip)+xc(ic)
	do 112 jc=jy1,jy2
	   jsy=jsym(jc)
	   dq1x3=(0.5*(uo(ic,jc,kp)+uo(im,jc,kp))
     %-0.5*(uo(ic,jc,km)+uo(im,jc,km)))/dz
	   dq3x1=(0.5*(wo(ip,jc,kc)+wo(ip,jc,km))
     %-0.5*(wo(ic,jsy,kc)+wo(ic,jsy,km)))/dr
	   voraz2(ic,jc,kkc)=dq1x3-dq3x1
  112	continue
c
	do 113 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 113 jc=jy1,jy2
	    dq1x3=(0.5*(uo(ic,jc,kp)+uo(im,jc,kp))
     %-0.5*(uo(ic,jc,km)+uo(im,jc,km)))/dz
	    dq3x1=(0.5*(wo(ip,jc,kc)+wo(ip,jc,km))
     %-0.5*(wo(im,jc,kc)+wo(im,jc,km)))/dr
	    voraz2(ic,jc,kkc)=dq1x3-dq3x1
  113 continue
  111 continue
c
      ! Azimuthal vorticity - Meridian sections !
      do 121 jjc=1,nsez14
	jc=jsez14(jjc)
	jsy=jsym(jc)
	do 121 kc=kz1,kz2
	  kp=kc+1
	  km=kc-1
	  dz=zc(kp)-zc(km)
	  ic=ix1
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)+xc(ic)
	  dq1x3=(0.5*(uo(ic,jc,kp)+uo(im,jc,kp))
     %-0.5*(uo(ic,jc,km)+uo(im,jc,km)))/dz
	  dq3x1=(0.5*(wo(ip,jc,kc)+wo(ip,jc,km))
     %-0.5*(wo(ic,jsy,kc)+wo(ic,jsy,km)))/dr
	  voraz14(ic,jjc,kc)=dq1x3-dq3x1
	  dq1x3=(0.5*(uo(ic,jsy,kp)+uo(im,jsy,kp))
     %-0.5*(uo(ic,jsy,km)+uo(im,jsy,km)))/dz
	  dq3x1=(0.5*(wo(ip,jsy,kc)+wo(ip,jsy,km))
     %-0.5*(wo(ic,jc,kc)+wo(ic,jc,km)))/dr
	  voraz14(ic,nsez14+jjc,kc)=dq1x3-dq3x1
c
	  do 121 ic=ix1+1,ix2
	    ip=ic+1
	    im=ic-1
	    dr=xc(ip)-xc(im)
	    dq1x3=(0.5*(uo(ic,jc,kp)+uo(im,jc,kp))
     %-0.5*(uo(ic,jc,km)+uo(im,jc,km)))/dz
	    dq3x1=(0.5*(wo(ip,jc,kc)+wo(ip,jc,km))
     %-0.5*(wo(im,jc,kc)+wo(im,jc,km)))/dr
	    voraz14(ic,jjc,kc)=dq1x3-dq3x1
	    dq1x3=(0.5*(uo(ic,jsy,kp)+uo(im,jsy,kp))
     %-0.5*(uo(ic,jsy,km)+uo(im,jsy,km)))/dz
	    dq3x1=(0.5*(wo(ip,jsy,kc)+wo(ip,jsy,km))
     %-0.5*(wo(im,jsy,kc)+wo(im,jsy,km)))/dr
	    voraz14(ic,nsez14+jjc,kc)=dq1x3-dq3x1
  121 continue
c
      ! Radial vorticity - Circumferential sections !
      do 131 kc=kz1,kz2
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	do 131 iic=1,nsez1
	  ic=isez1(iic)
	  do 131 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq3x2=(0.5*(wo(ic,jp,kc)+wo(ic,jp,km))
     %-0.5*(wo(ic,jm,kc)+wo(ic,jm,km)))/dtheta	  !
	    dq2x3=(0.5*(vo(ic,jc,kp)+vo(ic,jm,kp))
     %-0.5*(vo(ic,jc,km)+vo(ic,jm,km)))/dz
	    vorr1(iic,jc,kc)=dq3x2/xc(ic)-dq2x3
  131 continue
c
      ! Radial vorticity - Cross sections !
      do 141 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	do 141 ic=ix1,ix2
	  do 141 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq3x2=(0.5*(wo(ic,jp,kc)+wo(ic,jp,km))
     %-0.5*(wo(ic,jm,kc)+wo(ic,jm,km)))/dtheta	 !
	    dq2x3=(0.5*(vo(ic,jc,kp)+vo(ic,jm,kp))
     %-0.5*(vo(ic,jc,km)+vo(ic,jm,km)))/dz
	    vorr2(ic,jc,kkc)=dq3x2/xc(ic)-dq2x3
  141 continue
c
      ! Radial vorticity - Meridian sections !
      do 151 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jm=jmv(jc)
	jsy=jsym(jc)
	jpsy=jsym(jp)
	jmsy=jsym(jm)
	dtheta=2.*dely
	do 151 kc=kz1,kz2
	  kp=kc+1
	  km=kc-1
	  dz=zc(kp)-zc(km)
	  do 151 ic=ix1,ix2
	    dq3x2=(0.5*(wo(ic,jp,kc)+wo(ic,jp,km))
     %-0.5*(wo(ic,jm,kc)+wo(ic,jm,km)))/dtheta	 !
	    dq2x3=(0.5*(vo(ic,jc,kp)+vo(ic,jm,kp))
     %-0.5*(vo(ic,jc,km)+vo(ic,jm,km)))/dz
	    vorr14(ic,jjc,kc)=dq3x2/xc(ic)-dq2x3
	    dq3x2=(0.5*(wo(ic,jpsy,kc)+wo(ic,jpsy,km))
     %-0.5*(wo(ic,jmsy,kc)+wo(ic,jmsy,km)))/dtheta    !
	    dq2x3=(0.5*(vo(ic,jsy,kp)+vo(ic,jmsy,kp))
     %-0.5*(vo(ic,jsy,km)+vo(ic,jmsy,km)))/dz
	    vorr14(ic,nsez14+jjc,kc)=dq3x2/xc(ic)-dq2x3
  151 continue
c
      ! Axial vorticity - Circumferential sections !
      do 161 kc=kz1,kz2
	do 161 iic=1,nsez1
	  ic=isez1(iic)
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 161 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=
     %(0.5*(vo(ip,jc,kc)+vo(ip,jm,kc))*xc(ip)
     %-0.5*(vo(im,jc,kc)+vo(im,jm,kc))*xc(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jm,kc)+uo(im,jm,kc)))/dtheta	 !
	    vorz1(iic,jc,kc)=(dq2x1-dq1x2)/xc(ic)
  161 continue
c
      ! Axial vorticity - Cross sections !
      do 171 kkc=1,nsez2max
	kc=ksez2(kkc)
	ic=ix1
	ip=ic+1
	im=ic-1
	dr=xu(ic)-xu(im)
	do 172 jc=jy1,jy2
	  jp=jpv(jc)
	  jm=jmv(jc)
	  dtheta=2.*dely
	  dq2x1=(0.25*(vo(ic,jc,kc)+vo(ip,jc,kc)
     %+vo(ic,jm,kc)+vo(ip,jm,kc))*xu(ic))/dr
	  dq1x2=
     %(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jm,kc)+uo(im,jm,kc)))/dtheta	  !
	  vorz2(ic,jc,kkc)=(dq2x1-dq1x2)/xc(ic)
  172 continue
c
	do 173 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 173 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=
     %(0.5*(vo(ip,jc,kc)+vo(ip,jm,kc))*xc(ip)
     %-0.5*(vo(im,jc,kc)+vo(im,jm,kc))*xc(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jm,kc)+uo(im,jm,kc)))/dtheta	 !
	    vorz2(ic,jc,kkc)=(dq2x1-dq1x2)/xc(ic)
  173 continue
  171 continue
c
      ! Axial vorticity - Meridian sections !
      do 181 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jm=jmv(jc)
	jsy=jsym(jc)
	jpsy=jsym(jp)
	jmsy=jsym(jm)
	dtheta=2.*dely
	do 181 kc=kz1,kz2
	  ic=ix1
	  ip=ic+1
	  im=ic-1
	  dr=xu(ic)-xu(im)
	  dq2x1=(0.25*(vo(ip,jc,kc)+vo(ic,jc,kc)
     %+vo(ip,jm,kc)+vo(ic,jm,kc))*xu(ic))/dr
	  dq1x2=(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jm,kc)+uo(im,jm,kc)))/dtheta	  !
	  vorz14(ic,jjc,kc)=(dq2x1-dq1x2)/xc(ic)
	  dq2x1=(0.25*(vo(ip,jsy,kc)+vo(ic,jsy,kc)
     %+vo(ip,jmsy,kc)+vo(ic,jmsy,kc))*xu(ic))/dr
	  dq1x2=(0.5*(uo(ic,jpsy,kc)+uo(im,jpsy,kc))
     %-0.5*(uo(ic,jmsy,kc)+uo(im,jmsy,kc)))/dtheta   !
	  vorz14(ic,nsez14+jjc,kc)=(dq2x1-dq1x2)/xc(ic)
	  do 181 ic=ix1+1,ix2
	    ip=ic+1
	    im=ic-1
	    dr=xc(ip)-xc(im)
	    dq2x1=
     %(0.5*(vo(ip,jc,kc)+vo(ip,jm,kc))*xc(ip)
     %-0.5*(vo(im,jc,kc)+vo(im,jm,kc))*xc(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jm,kc)+uo(im,jm,kc)))/dtheta	 !
	    vorz14(ic,jjc,kc)=(dq2x1-dq1x2)/xc(ic)
	    dq2x1=
     %(0.5*(vo(ip,jsy,kc)+vo(ip,jmsy,kc))*xc(ip)
     %-0.5*(vo(im,jsy,kc)+vo(im,jmsy,kc))*xc(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jpsy,kc)+uo(im,jpsy,kc))
     %-0.5*(uo(ic,jmsy,kc)+uo(im,jmsy,kc)))/dtheta   !
	    vorz14(ic,nsez14+jjc,kc)=(dq2x1-dq1x2)/xc(ic)
  181 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine plot2D(iplant2,iplant2c,iplant2d,iplant2p,
     %nx,ny,nz,nzg,ntime,nbd,mbd,nfacet,unvect,vertex,xc,yc,zcg,
     %xc_car,yc_car,flagpo,tv,vo,uo,wo,p)
c
      use sections
      include'common.h'
      include'immersed.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer iplant2,iplant2c,iplant2d,iplant2p
      integer nx,ny,nz,nzg,ntime,nbd,mbd,nfacet
      integer flagpo(nx,ny,nz,nbd)
      real unvect(3,nfacet),vertex(3,3,nfacet)
      real xc(nx),yc(ny),zcg(nzg),xc_car(nx,ny),yc_car(nx,ny),
     %tv(nx,ny,nz),vo(nx,ny,nz),uo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
c
c Local variables
      integer iic,iicd,iicu,is,im,i,jjc,jjcd,jjcu,js,jm,j,jsy,jmsy,
     %kkc,kkcd,kkcu,ks,km,k,l,i1,jjp,ir,jr,kr,lr,krg,kk,kkd,kku,
     %kk1,kk2,ibd
      integer STATUS(MPI_STATUS_SIZE)
!      integer flagpor(ny,nzg,nbd),
!     %flagpor2(nx,nzg,nbd),flagpor2sym(nx,nzg,nbd)
      real recvar(ny,nz),recvar2(nx,nz),recvar2sym(nx,nz)
c
c Parameters
      integer lamb2
      real deltasect
      parameter (lamb2=1,deltasect=2.*pi)
c
   8  format(10f25.8)
   9  format(10f25.8)
 100  format(3i8)
 101  format(i8,3i2)
 102  format(i2)
      if(iplant2.eq.0) then
	 if(myrank.eq.0) then
!	     write(60,*)'iplant2',iplant2
	    do iic=1,nsez1
	       iicd=iic/10
	       iicu=mod(iic,10)
	       open
     %(81,file='grid1_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'.xyz')
	       open
     %(82,file='grid1_new_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'.xyz')
	       is=isez1(iic)
	       write(81,100)(ny-2)*lamb2+1,1,nzg-2
	       write(82,100)(ny-2)*lamb2+1,1,nzg-2
	       do 1 k=2,nzg-1
	       do 10 l=1,lamb2
	       do 10 j=jy1,jy2
	       write(81,9)xc_car(is,jy1)
   10	       write(82,9)xc(is)*cos(yc(j)+(l-1)*2.*pi/lamb2)
	       write(81,9)xc_car(is,jy1)
   1	       write(82,9)xc(is)*cos(yc(jy1)+lamb2*deltasect)
	       do 2 k=2,nzg-1
	       do 20 l=1,lamb2
	       do 20 j=jy1,jy2
	       write(81,9)
     %yc_car(is,jy1)+xc(is)*(yc(j)+(l-1)*2.*pi/lamb2)
   20	       write(82,9)xc(is)*sin(yc(j)+(l-1)*2.*pi/lamb2)
	       write(81,9)
     %yc_car(is,jy1)+xc(is)*(yc(jy1)+deltasect*lamb2)
   2	       write(82,9)xc(is)*sin(yc(jy1)+lamb2*deltasect)
	       do 3 k=2,nzg-1
	       do 30 l=1,lamb2
	       do 30 j=jy1,jy2
	       write(81,9)zcg(k)
   30	       write(82,9)zcg(k)
	       write(81,9)zcg(k)
   3	       write(82,9)zcg(k)
	       close(81)
	       close(82)
	    enddo
	    do kkc=1,nsez2
	       kkcd=kkc/10
	       kkcu=mod(kkc,10)
	       open(81,file='grid2_2D_slice'//char(48+kkcd)
     %//char(48+kkcu)//'.xyz')
	       ks=ksez2g(kkc)
	       write(81,100)(ny-2)*lamb2+1,nx-2,1
	       do 11 i=ix1,ix2
	       do 110 l=1,lamb2
	       do 110 j=jy1,jy2
  110	       write(81,9)xc(i)*cos(yc(j)+(l-1)*2.*pi/lamb2)
  11	       write(81,9)xc(i)*cos(yc(jy1)+lamb2*deltasect)
	       do 12 i=ix1,ix2
	       do 120 l=1,lamb2
	       do 120 j=jy1,jy2
  120	       write(81,9)xc(i)*sin(yc(j)+(l-1)*2.*pi/lamb2)
  12	       write(81,9)xc(i)*sin(yc(jy1)+lamb2*deltasect)
	       do 13 i=ix1,ix2
	       do 130 l=1,lamb2
	       do 130 j=jy1,jy2
  130	       write(81,9)zcg(ks)
  13	       write(81,9)zcg(ks)
	       close(81)
	    enddo
	    do jjc=1,nsez14
	       jjcd=jjc/10
	       jjcu=mod(jjc,10)
	       open(81,file='grid3_2D_slice'//char(48+jjcd)
     %//char(48+jjcu)//'.xyz')
	       js=jsez14(jjc)
	       jsy=jsym(js)
	       write(81,100)1,2*(nx-2),nzg-2
	       do 210 k=2,nzg-1
	       do 211 i=ix2,ix1,-1
  211	       write(81,9)xc(i)*cos(yc(jsy))
	       do 212 i=ix1,ix2,1
  212	       write(81,9)xc(i)*cos(yc(js))
  210	       continue
	       do 220 k=2,nzg-1
	       do 221 i=ix2,ix1,-1
  221	       write(81,9)xc(i)*sin(yc(jsy))
	       do 222 i=ix1,ix2,1
  222	       write(81,9)xc(i)*sin(yc(js))
  220	       continue
	       do 230 k=2,nzg-1
	       do 230 i=1,2*(nx-2)
  230	       write(81,9)zcg(k)
	       close(81)
	    enddo
	 endif
	 iplant2=1
!	  return
      endif
c
      iplant2p=iplant2
      iplant2c=iplant2p/100
      iplant2p=mod(iplant2p,100)
      iplant2d=iplant2p/10
      iplant2p=mod(iplant2p,10)
c
      do iic=1,nsez1
	 is=isez1(iic)
	 im=is-1
	 if(myrank.eq.0) then
	    iicd=iic/10
	    iicu=mod(iic,10)
	    open(82,file='results1_2D_'//char(48+iplant2c)
     %//char(48+iplant2d)//char(48+iplant2p)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'.q')
	    write(82,100)(ny-2)*lamb2+1,1,nzg-2
	    i1=1
	    write(82,101)ntime,i1,i1,i1
	    do 40 k=kz1,kz2
	    do 41 l=1,lamb2
	    do 41 j=jy1,jy2
!	    if(all(flagpo(is,j,k,:)).le.0) then
!	       jm=jmv(j)
	       write(82,9)
!     %(vo(is,j,k)+vo(is,jm,k))*0.5-amp*xc(is)
     %tv(is,j,k)/ru1
!	    else
!	       write(82,9)1.e+03
!	    endif
 41	    continue
!	    if(all(flagpo(is,jy1,k,:)).le.0) then
	       write(82,9)
!     %(vo(is,jy1,k)+vo(is,jy2,k))*0.5-amp*xc(is)
     %tv(is,jy1,k)/ru1
!	    else
!	       write(82,9)1.e+03
!	    endif
 40	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
!		  kk1=kz1+jjp*(nz-2)
!		  kk2=kz2+jjp*(nz-2)
!		  CALL MPI_RECV
!    %(flagpor(jy1:jy2,kk1:kk2,:),
!    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,jjp,1,
!    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),(ny-2)*(nz-2),MTYPE,jjp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 42 kr=kz1,kz2
!		  krg=kr+jjp*(nz-2)
		  do 43 lr=1,lamb2
		  do 43 jr=jy1,jy2
!		  if(all(flagpor(jr,krg,:)).le.0) then
		     write(82,9)
!     %recvar(jr,kr)-amp*xc(is)
     %recvar(jr,kr)/ru1
!		  else
!		     write(82,9)1.e+03
!		  endif
  43		  continue
!		  if(all(flagpor(jy1,krg,:)).le.0) then
		     write(82,9)
!     %recvar(jy1,kr)-amp*xc(is)
     %recvar(jy1,kr)/ru1
!		  else
!		     write(82,9)1.e+03
!		  endif
  42		  continue
	       enddo
	    endif
!	    write(82,*)(((i1,j=1,ny-1),i=1,1),k=2,nzg-1)   !!!!!!
!	    write(82,102)(((i1,j=1,ny-1),i=1,1),k=2,nzg-1)
c
	    do 50 k=kz1,kz2
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2
!	    if(all(flagpo(is,j,k,:)).le.0) then
	       jm=jmv(j)
	       write(82,9)(vo(is,j,k)+vo(is,jm,k))*0.5
!	    else
!	       write(82,9)1.e+03
!	    endif
   51	    continue
!	    if(all(flagpo(is,jy1,k,:)).le.0) then
	       write(82,9)(vo(is,jy1,k)+vo(is,jy2,k))*0.5
!	    else
!	       write(82,9)1.e+03
!	    endif
   50	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),(ny-2)*(nz-2),MTYPE,jjp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 52 kr=kz1,kz2
!		  krg=kr+jjp*(nz-2)
		  do 53 lr=1,lamb2
		  do 53 jr=jy1,jy2
!		  if(all(flagpor(jr,krg,:)).le.0) then
		     write(82,9)
     %recvar(jr,kr)
!		  else
!		     write(82,9)1.e+03
!		  endif
  53		  continue
!		  if(all(flagpor(jy1,krg,:)).le.0) then
		     write(82,9)
     %recvar(jy1,kr)
!		  else
!		     write(82,9)1.e+03
!		  endif
  52		  continue
	       enddo
	    endif
c
	    do 60 k=kz1,kz2
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2
!	    if(all(flagpo(is,j,k,:)).le.0) then
	       write(82,9)(uo(is,j,k)+uo(im,j,k))*0.5
!	    else
!	       write(82,9)1.e+03
!	    endif
  61	    continue
!	    if(all(flagpo(is,jy1,k,:)).le.0) then
	       write(82,9)(uo(is,jy1,k)+uo(im,jy1,k))*0.5
!	    else
!	       write(82,9)1.e+03
!	    endif
  60	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),(ny-2)*(nz-2),MTYPE,jjp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 62 kr=kz1,kz2
!		  krg=kr+jjp*(nz-2)
		  do 63 lr=1,lamb2
		  do 63 jr=jy1,jy2
!		  if(all(flagpor(jr,krg,:)).le.0) then
		     write(82,9)recvar(jr,kr)
!		  else
!		     write(82,9)1.e+03
!		  endif
  63		  continue
!		  if(all(flagpor(jy1,krg,:)).le.0) then
		     write(82,9)recvar(jy1,kr)
!		  else
!		     write(82,9)1.e+03
!		  endif
  62		  continue
	       enddo
	    endif
c
	    do 70 k=kz1,kz2
	    km=k-1
	    do 71 l=1,lamb2
	    do 71 j=jy1,jy2
!	    if(all(flagpo(is,j,k,:)).le.0) then
	       write(82,9)(wo(is,j,k)+wo(is,j,km))*0.5
!	    else
!	       write(82,9)1.e+03
!	    endif
  71	    continue
!	    if(all(flagpo(is,jy1,k,:)).le.0) then
	      write(82,9)(wo(is,jy1,k)+wo(is,jy1,km))*0.5
!	    else
!	      write(82,9)1.e+03
!	    endif
  70	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),(ny-2)*(nz-2),MTYPE,jjp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 72 kr=kz1,kz2
!		  krg=kr+jjp*(nz-2)
		  do 73 lr=1,lamb2
		  do 73 jr=jy1,jy2
!		  if(all(flagpor(jr,krg,:)).le.0) then
		     write(82,9)
     %recvar(jr,kr)
!		  else
!		     write(82,9)1.e+03
!		  endif
  73		  continue
!		  if(all(flagpor(jy1,krg,:)).le.0) then
		     write(82,9)
     %recvar(jy1,kr)
!		  else
!		     write(82,9)1.e+03
!		  endif
  72		  continue
	       enddo
	    endif
c
	    do 80 k=kz1,kz2
	    do 81 l=1,lamb2
	    do 81 j=jy1,jy2
!	    if(all(flagpo(is,j,k,:)).le.0) then
	       write(82,9)p(is,j,k)    !*ros
!	    else
!	       write(82,9)1.e+03
!	    endif
  81	    continue
!	    if(all(flagpo(is,jy1,k,:)).le.0) then
	      write(82,9)p(is,jy1,k)   !*ros
!	    else
!	      write(82,9)1.e+03
!	    endif
  80	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),(ny-2)*(nz-2),MTYPE,jjp,6,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 82 kr=kz1,kz2
!		  krg=kr+jjp*(nz-2)
		  do 83 lr=1,lamb2
		  do 83 jr=jy1,jy2
!		  if(all(flagpor(jr,krg,:)).le.0) then
		     write(82,9)recvar(jr,kr)	!*ros
!		  else
!		     write(82,9)1.e+03
!		  endif
  83		  continue
!		  if(all(flagpor(jy1,krg,:)).le.0) then
		     write(82,9)recvar(jy1,kr)	!*ros
!		  else
!		     write(82,9)1.e+03
!		  endif
  82		  continue
	       enddo
	    endif
	    close(82)
	 else
!	    CALL MPI_SEND
!    %(flagpo(is,jy1:jy2,kz1:kz2,:),
!    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,0,1,
!    %MPI_COMM_EDDY,STATUS,IERR)
!	    CALL MPI_SEND
!    %((vo(is,1:jy2-1,kz1:kz2)+vo(is,jy1:jy2,kz1:kz2))*0.5,
!    %1*(ny-2)*(nz-2),MTYPE,0,2,
!    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(tv(is,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((vo(is,1:jy2-1,kz1:kz2)+vo(is,jy1:jy2,kz1:kz2))*0.5,
     %1*(ny-2)*(nz-2),MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((uo(is,jy1:jy2,kz1:kz2)+uo(im,jy1:jy2,kz1:kz2))*0.5,
     %1*(ny-2)*(nz-2),MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((wo(is,jy1:jy2,1:kz2-1)+wo(is,jy1:jy2,kz1:kz2))*0.5,
     %1*(ny-2)*(nz-2),MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(p(is,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
      do kkc=1,nsez2max
	 ks=ksez2(kkc)
	 km=ks-1
	 kk=isezindx2+(kkc-1)
	 kkd=kk/10
	 kku=mod(kk,10)
	 open(82,file='results2_2D_'//char(48+iplant2c)
     %//char(48+iplant2d)//char(48+iplant2p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'.q')
	 write(82,100)(ny-2)*lamb2+1,nx-2,1
	 i1=1
	 write(82,101)ntime,i1,i1,i1
	 do 140 i=ix1,ix2
	 do 141 l=1,lamb2
	 do 141 j=jy1,jy2
!	 if(all(flagpo(i,j,ks,:)).le.0) then
!	    jm=jmv(j)
	    write(82,9)tv(i,j,ks)/ru1
!    %(vo(i,j,ks)+vo(i,jm,ks))*0.5-amp*xc(i)
!	 else
!	    write(82,9)1.e+03
!	 endif
 141	 continue
!	 if(all(flagpo(i,jy1,ks,:)).le.0) then
	    write(82,9)tv(i,jy1,ks)/ru1
!    %(vo(i,jy1,ks)+vo(i,jy2,ks))*0.5-amp*xc(i)
!	 else
!	    write(82,9)1.e+03
!	 endif
 140	 continue
!	 write(82,*)(((i1,j=1,ny-1),i=ix1,ix2),k=1,1)	!!!!!!
!	 write(82,102)(((i1,j=1,ny-1),i=ix1,ix2),k=1,1)
	 do 150 i=ix1,ix2
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2
!	 if(all(flagpo(i,j,ks,:)).le.0) then
	    jm=jmv(j)
	    write(82,9)(vo(i,j,ks)+vo(i,jm,ks))*0.5
!	 else
!	    write(82,9)1.e+03
!	 endif
 151	 continue
!	 if(all(flagpo(i,jy1,ks,:)).le.0) then
	    write(82,9)(vo(i,jy1,ks)+vo(i,jy2,ks))*0.5
!	 else
!	    write(82,9)1.e+03
!	 endif
 150	 continue
	 do 160 i=ix1,ix2
	 im=i-1
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2
!	 if(all(flagpo(i,j,ks,:)).le.0) then
	    write(82,9)(uo(i,j,ks)+uo(im,j,ks))*0.5
!	 else
!	    write(82,9)1.e+03
!	 endif
 161	 continue
!	 if(all(flagpo(i,jy1,ks,:)).le.0) then
	    write(82,9)(uo(i,jy1,ks)+uo(im,jy1,ks))*0.5
!	 else
!	    write(82,9)1.e+03
!	 endif
 160	 continue
	 do 170 i=ix1,ix2
	 do 171 l=1,lamb2
	 do 171 j=jy1,jy2
!	 if(all(flagpo(i,j,ks,:)).le.0) then
	    write(82,9)(wo(i,j,ks)+wo(i,j,km))*0.5
!	 else
!	    write(82,9)1.e+03
!	 endif
 171	 continue
!	 if(all(flagpo(i,jy1,ks,:)).le.0) then
	    write(82,9)(wo(i,jy1,ks)+wo(i,jy1,km))*0.5
!	 else
!	    write(82,9)1.e+03
!	 endif
 170	 continue
	 do 180 i=ix1,ix2
	 do 181 l=1,lamb2
	 do 181 j=jy1,jy2
!	 if(all(flagpo(i,j,ks,:)).le.0) then
	    write(82,8)p(i,j,ks)   !*ros
!	 else
!	    write(82,8)1.e+03
!	 endif
 181	 continue
!	 if(all(flagpo(i,jy1,ks,:)).le.0) then
	    write(82,8)p(i,jy1,ks)    !*ros
!	 else
!	    write(82,8)1.e+03
!	 endif
 180	 continue
	 close(82)
      enddo
c
      do jjc=1,nsez14
	 js=jsez14(jjc)
	 jsy=jsym(js)
	 jm=jmv(js)
	 jmsy=jmv(jsy)
	 if(myrank.eq.0) then
	    jjcd=jjc/10
	    jjcu=mod(jjc,10)
	    open(82,file='results3_2D_'//char(48+iplant2c)
     %//char(48+iplant2d)//char(48+iplant2p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'.q')
	    write(82,100)1,2*(nx-2),nzg-2
	    i1=1
	    write(82,101)ntime,i1,i1,i1
	    do 240 k=kz1,kz2
	    do 241 i=ix2,ix1,-1
!	    if(all(flagpo(i,jsy,k,:)).le.0) then
	       write(82,9)
!    %(vo(i,jsy,k)+vo(i,jmsy,k))*0.5-amp*xc(i)
     %tv(i,jsy,k)/ru1
!	    else
!	       write(82,9)1.e+03
!	    endif
  241	    continue
	    do 242 i=ix1,ix2,1
!	    if(all(flagpo(i,js,k,:)).le.0) then
	       write(82,9)
!    %(vo(i,js,k)+vo(i,jm,k))*0.5-amp*xc(i)
     %tv(i,js,k)/ru1
!	    else
!	       write(82,9)1.e+03
!	    endif
  242	    continue
  240	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
!		  kk1=kz1+jjp*(nz-2)
!		  kk2=kz2+jjp*(nz-2)
!		  CALL MPI_RECV
!    %(flagpor2(ix1:ix2,kk1:kk2,:),
!    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,7,
!    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
!		  CALL MPI_RECV
!    %(flagpor2sym(ix1:ix2,kk1:kk2,:),
!    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,9,
!    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 243 kr=kz1,kz2
!		  krg=kr+jjp*(nz-2)
		  do 244 ir=ix2,ix1,-1
!		  if(all(flagpor2sym(ir,krg,:)).le.0) then
		     write(82,9)
!    %recvar2sym(ir,kr)-amp*xc(ir)
     %recvar2sym(ir,kr)/ru1
!		  else
!		     write(82,9)1.e+03
!		  endif
  244		  continue
		  do 245 ir=ix1,ix2,1
!		  if(all(flagpor2(ir,krg,:)).le.0) then
		     write(82,9)
!    %recvar2(ir,kr)*0.5-amp*xc(ir)
     %recvar2(ir,kr)/ru1
!		  else
!		    write(82,9)1.e+03
!		  endif
  245		  continue
  243		  continue
	       enddo
	    endif
!	    write(82,*)(((i1,j=1,1),i=1,2*(nx-2)),k=2,nzg-1)   !!!!!!
!	    write(82,102)(((i1,j=1,1),i=1,2*(nx-2)),k=2,nzg-1)
c
	    do 250 k=kz1,kz2
	    do 251 i=ix2,ix1,-1
!	    if(all(flagpo(i,jsy,k,:)).le.0) then
	       write(82,9)(vo(i,jsy,k)+vo(i,jmsy,k))*0.5
!	    else
!	       write(82,9)1.e+03
!	    endif
  251	    continue
	    do 252 i=ix1,ix2,1
!	    if(all(flagpo(i,js,k,:)).le.0) then
	       write(82,9)(vo(i,js,k)+vo(i,jm,k))*0.5
!	    else
!	       write(82,9)1.e+03
!	    endif
  252	    continue
  250	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 253 kr=kz1,kz2
!		  krg=kr+jjp*(nz-2)
		  do 254 ir=ix2,ix1,-1
!		  if(all(flagpor2sym(ir,krg,:)).le.0) then
		     write(82,9)recvar2sym(ir,kr)
!		  else
!		     write(82,9)1.e+03
!		  endif
  254		  continue
		  do 255 ir=ix1,ix2,1
!		  if(all(flagpor2(ir,krg,:)).le.0) then
		     write(82,9)recvar2(ir,kr)
!		  else
!		     write(82,9)1.e+03
!		  endif
  255		  continue
  253		  continue
	       enddo
	    endif
c
	    do 260 k=kz1,kz2
	    do 261 i=ix2,ix1,-1
!	    if(all(flagpo(i,jsy,k,:)).le.0) then
	       im=i-1
	       write(82,9)(uo(i,jsy,k)+uo(im,jsy,k))*0.5
!	    else
!	       write(82,9)1.e+03
!	    endif
  261	    continue
	    do 262 i=ix1,ix2,1
!	    if(all(flagpo(i,js,k,:)).le.0) then
	       im=i-1
	       write(82,9)(uo(i,js,k)+uo(im,js,k))*0.5
!	    else
!	       write(82,9)1.e+03
!	    endif
  262	    continue
  260	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 263 kr=kz1,kz2
!		  krg=kr+jjp*(nz-2)
		  do 264 ir=ix2,ix1,-1
!		  if(all(flagpor2sym(ir,krg,:)).le.0) then
		     write(82,9)
     %recvar2sym(ir,kr)
!		  else
!		     write(82,9)1.e+03
!		  endif
  264		  continue
		  do 265 ir=ix1,ix2,1
!		  if(all(flagpor2(ir,krg,:)).le.0) then
		     write(82,9)
     %recvar2(ir,kr)
!		  else
!		     write(82,9)1.e+03
!		  endif
  265		  continue
  263		  continue
	       enddo
	    endif
c
	    do 270 k=kz1,kz2
	    km=k-1
	    do 271 i=ix2,ix1,-1
!	    if(all(flagpo(i,jsy,k,:)).le.0) then
	       write(82,9)(wo(i,jsy,k)+wo(i,jsy,km))*0.5
!	    else
!	       write(82,9)1.e+03
!	    endif
  271	    continue
	    do 272 i=ix1,ix2,1
!	    if(all(flagpo(i,js,k,:)).le.0) then
	       write(82,9)(wo(i,js,k)+wo(i,js,km))*0.5
!	    else
!	       write(82,9)1.e+03
!	    endif
  272	    continue
  270	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,16,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 273 kr=kz1,kz2
!		  krg=kr+jjp*(nz-2)
		  do 274 ir=ix2,ix1,-1
!		  if(all(flagpor2sym(ir,krg,:)).le.0) then
		     write(82,9)
     %recvar2sym(ir,kr)
!		  else
!		     write(82,9)1.e+03
!		  endif
  274		  continue
		  do 275 ir=ix1,ix2,1
!		  if(all(flagpor2(ir,krg,:)).le.0) then
		     write(82,9)
     %recvar2(ir,kr)
!		  else
!		     write(82,9)1.e+03
!		  endif
  275		  continue
  273		  continue
	       enddo
	    endif
c
	    do 280 k=kz1,kz2
	    do 281 i=ix2,ix1,-1
!	    if(all(flagpo(i,jsy,k,:)).le.0) then
	       write(82,8)p(i,jsy,k)	!*ros
!	    else
!	       write(82,8)1.e+03
!	    endif
  281	    continue
	    do 282 i=ix1,ix2,1
!	    if(all(flagpo(i,js,k,:)).le.0) then
	       write(82,8)p(i,js,k)    !*ros
!	    else
!	       write(82,8)1.e+03
!	    endif
  282	    continue
  280	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,17,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,18,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 283 kr=kz1,kz2
!		  krg=kr+jjp*(nz-2)
		  do 284 ir=ix2,ix1,-1
!		  if(all(flagpor2sym(ir,krg,:)).le.0) then
		     write(82,8)recvar2sym(ir,kr)    !*ros
!		  else
!		     write(82,8)1.e+03
!		  endif
  284		  continue
		  do 285 ir=ix1,ix2,1
!		  if(all(flagpor2(ir,krg,:)).le.0) then
		     write(82,8)recvar2(ir,kr)	  !*ros
!		  else
!		     write(82,8)1.e+03
!		  endif
  285		  continue
  283		  continue
	       enddo
	    endif
	    close(82)
	 else
c	    CALL MPI_SEND
c    %(flagpo(ix1:ix2,js,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,7,
c    %MPI_COMM_EDDY,STATUS,IERR)
c	    CALL MPI_SEND
c    %((vo(ix1:ix2,js,kz1:kz2)+vo(ix1:ix2,jm,kz1:kz2))*0.5,
c    %(nx-2)*1*(nz-2),MTYPE,0,8,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(tv(ix1:ix2,js,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
c	    CALL MPI_SEND
c    %(flagpo(ix1:ix2,jsy,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,9,
c    %MPI_COMM_EDDY,STATUS,IERR)
c	    CALL MPI_SEND
c    %((vo(ix1:ix2,jsy,kz1:kz2)+vo(ix1:ix2,jmsy,kz1:kz2))*0.5,
c    %(nx-2)*1*(nz-2),MTYPE,0,10,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(tv(ix1:ix2,jsy,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((vo(ix1:ix2,js,kz1:kz2)+vo(ix1:ix2,jm,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((vo(ix1:ix2,jsy,kz1:kz2)+vo(ix1:ix2,jmsy,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((uo(1:ix2-1,js,kz1:kz2)+uo(ix1:ix2,js,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((uo(1:ix2-1,jsy,kz1:kz2)+uo(ix1:ix2,jsy,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((wo(ix1:ix2,js,1:kz2-1)+wo(ix1:ix2,js,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((wo(ix1:ix2,jsy,1:kz2-1)+wo(ix1:ix2,jsy,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,16,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(p(ix1:ix2,js,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,17,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(p(ix1:ix2,jsy,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,18,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
 103  format(2x,a12,1x,e14.7,1x,e14.7,1x,e14.7)
 104  format(6x,a6,3x,e14.7,1x,e14.7,1x,e14.7)
      if((ibm.eq.2).and.(myrank.eq.0)) then
	 do ibd=mbd,nbd
	    open(83,file='ib_n.'//index(ibd)//'_'//char(48+iplant2c)
     %//char(48+iplant2d)//char(48+iplant2p)//'.stl')
	    write(83,*)'solid  OBJECT'
	    do j=lb(ibd)+1,lb(ibd)+mb(ibd)
	       write(83,103)'facet normal',unvect(1,j),unvect(2,j),
     %unvect(3,j)
	       write(83,*)'outer loop'
	       write(83,104)
     %'vertex',vertex(1,1,j),vertex(2,1,j),vertex(3,1,j)
	       write(83,104)
     %'vertex',vertex(1,2,j),vertex(2,2,j),vertex(3,2,j)
	       write(83,104)
     %'vertex',vertex(1,3,j),vertex(2,3,j),vertex(3,3,j)
	       write(83,*)'endloop'
	       write(83,*)'endfacet'
	    enddo
	    write(83,*)'endsolid  OBJECT'
	    close(83)
	 enddo
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine plot2Dvort(iplant2c,iplant2d,iplant2u,
     %nx,ny,nz,nzg,nbd,ntime,flagpo,vo,uo,wo,p)
c
      use sections
      use vorticity
      include'common.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer iplant2c,iplant2d,iplant2u
      integer nx,ny,nz,nzg,nbd,ntime
      integer flagpo(nx,ny,nz,nbd)
      real p(nx,ny,nz),vo(nx,ny,nz),uo(nx,ny,nz),wo(nx,ny,nz)
c
c Local variables
      integer iic,iicd,iicu,i,is,im,jjc,jjcd,jjcu,j,js,jm,jsy,jmsy,
     %kkc,k,ks,km,i1,l,ir,jr,kr,lr,jjp,kk,kkd,kku,kk1,kk2,krg
      integer STATUS(MPI_STATUS_SIZE)
!      integer flagpor(ny,nzg,nbd),
!     %flagpor2(nx,nzg,nbd),flagpor2sym(nx,nzg,nbd)
      real recvar(ny,nz),recvar2(nx,nz),recvar2sym(nx,nz)
c
c Parameters
      integer lamb2,nxp,nyp,nzp
      parameter (lamb2=1)
      parameter (nyp=1,nxp=1,nzp=1)
c
   8  format(10f25.8)
   9  format(10f25.8)
 100  format(3i8)
 101  format(i8,3i2)
 102  format(i2)
c
      do iic=1,nsez1
	 is=isez1(iic)
	 im=is-1
	 if(myrank.eq.0) then
	    iicd=iic/10
	    iicu=mod(iic,10)
	    open(83,file='results1_2D_'//char(48+iplant2c)
     %//char(48+iplant2d)//char(48+iplant2u)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'_vort2D.q')
	    write(83,100)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	    i1=1
	    write(83,101)ntime,i1,i1,i1
!!!!!!		  write(83,*)
	    write(83,102)
     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
	    do 40 k=kz1,kz2,nzp
	    do 41 l=1,lamb2
	    do 41 j=jy1,jy2,nyp
!	    if(all(flagpo(is,j,k,:)).le.0) then
	       write(83,9)voraz1(iic,j,k)
!	    else
!	       write(83,9)1.e+03
!	    endif
   41	    continue
!	    if(all(flagpo(is,jy1,k,:)).le.0) then
	       write(83,9)voraz1(iic,jy1,k)
!	    else
!	       write(83,9)1.e+03
!	    endif
   40	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
!		  kk1=kz1+jjp*(nz-2)
!		  kk2=kz2+jjp*(nz-2)
!		  CALL MPI_RECV
!    %(flagpor(jy1:jy2,kk1:kk2,:),
!    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,jjp,1,
!    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),(ny-2)*(nz-2),MTYPE,jjp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 42 kr=kz1,kz2,nzp
!		  krg=kr+jjp*(nz-2)
		  do 43 lr=1,lamb2
		  do 43 jr=jy1,jy2,nyp
!		  if(all(flagpor(jr,krg,:)).le.0) then
		     write(83,9)recvar(jr,kr)
!		  else
!		     write(83,9)1.e+03
!		  endif
   43		  continue
!		  if(all(flagpor(jy1,krg,:)).le.0) then
		     write(83,9)recvar(jy1,kr)
!		  else
!		     write(83,9)1.e+03
!		  endif
   42		  continue
	       enddo
	    endif
c
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
!	    if(all(flagpo(is,j,k,:)).le.0) then
	       write(83,9)vorr1(iic,j,k)
!	    else
!	       write(83,9)1.e+03
!	    endif
   51	    continue
!	    if(all(flagpo(is,jy1,k,:)).le.0) then
	       write(83,9)vorr1(iic,jy1,k)
!	    else
!	       write(83,9)1.e+03
!	    endif
   50	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),(ny-2)*(nz-2),MTYPE,jjp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 52 kr=kz1,kz2,nzp
!		  krg=kr+jjp*(nz-2)
		  do 53 lr=1,lamb2
		  do 53 jr=jy1,jy2,nyp
!		  if(all(flagpor(jr,krg,:)).le.0) then
		     write(83,9)recvar(jr,kr)
!		  else
!		     write(83,9)1.e+03
!		  endif
   53		  continue
!		  if(all(flagpor(jy1,krg,:)).le.0) then
		     write(83,9)recvar(jy1,kr)
!		  else
!		     write(83,9)1.e+03
!		  endif
   52		  continue
	       enddo
	    endif
c
	    do 60 k=kz1,kz2,nzp
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2,nyp
!	    if(all(flagpo(is,j,k,:)).le.0) then
	       write(83,9)vorz1(iic,j,k)
!	    else
!	       write(83,9)1.e+03
!	    endif
   61	    continue
!	    if(all(flagpo(is,jy1,k,:)).le.0) then
	       write(83,9)vorz1(iic,jy1,k)
!	    else
!	       write(83,9)1.e+03
!	    endif
   60	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),(ny-2)*(nz-2),MTYPE,jjp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 62 kr=kz1,kz2,nzp
!		  krg=kr+jjp*(nz-2)
		  do 63 lr=1,lamb2
		  do 63 jr=jy1,jy2,nyp
!		  if(all(flagpor(jr,krg,:)).le.0) then
		     write(83,9)recvar(jr,kr)
!		  else
!		     write(83,9)1.e+03
!		  endif
   63		  continue
!		  if(all(flagpor(jy1,krg,:)).le.0) then
		     write(83,9)recvar(jy1,kr)
!		  else
!		     write(83,9)1.e+03
!		  endif
   62		  continue
	       enddo
	    endif
c
	    do 70 k=kz1,kz2,nzp
	    km=k-1
	    do 71 l=1,lamb2
	    do 71 j=jy1,jy2,nyp
!	    if(all(flagpo(is,j,k,:)).le.0) then
	       jm=jmv(j)
	       write(83,8)
     %(p(is,j,k)+0.5*0.25*((vo(is,j,k)+vo(is,jm,k))**2.+
     %(uo(is,j,k)+uo(im,j,k))**2.+
     %(wo(is,j,k)+wo(is,j,km))**2.))	 !*ros
!	    else
!	       write(83,8)1.e+03
!	    endif
   71	    continue
!	    if(all(flagpo(is,jy1,k,:)).le.0) then
	       write(83,8)
     %(p(is,jy1,k)+0.5*0.25*((vo(is,jy1,k)+vo(is,jy2,k))**2.+
     %(uo(is,jy1,k)+uo(im,jy1,k))**2.+
     %(wo(is,jy1,k)+wo(is,jy1,km))**2.))    !*ros
!	    else
!	       write(83,8)1.e+03
!	    endif
   70	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),(ny-2)*(nz-2),MTYPE,jjp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 72 kr=kz1,kz2,nzp
!		  krg=kr+jjp*(nz-2)
		  do 73 lr=1,lamb2
		  do 73 jr=jy1,jy2,nyp
!		  if(all(flagpor(jr,krg,:)).le.0) then
		     write(83,8)recvar(jr,kr)
!		  else
!		     write(83,8)1.e+03
!		  endif
   73		  continue
!		  if(all(flagpor(jy1,krg,:)).le.0) then
		     write(83,8)recvar(jy1,kr)
!		  else
!		     write(83,8)1.e+03
!		  endif
   72		  continue
	       enddo
	    endif
	    close(83)
	 else
!	    CALL MPI_SEND
!    %(flagpo(is,jy1:jy2,kz1:kz2,:),
!    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,0,1,
!    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(voraz1(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorr1(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorz1(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    do 700 k=kz1,kz2,nzp
	    km=k-1
	    do 710 l=1,lamb2
	    do 710 j=jy1,jy2,nyp
!	    if(all(flagpo(is,j,k,:)).le.0) then
	       jm=jmv(j)
	       recvar(j,k)=
     %(p(is,j,k)+0.5*0.25*((vo(is,j,k)+vo(is,jm,k))**2.+
     %(uo(is,j,k)+uo(im,j,k))**2.+
     %(wo(is,j,k)+wo(is,j,km))**2.))	 !*ros
!	    endif
  710	    continue
!	    if(all(flagpo(is,jy1,k,:)).le.0) then
	       recvar(jy1,k)=
     %(p(is,jy1,k)+0.5*0.25*((vo(is,jy1,k)+vo(is,jy2,k))**2.+
     %(uo(is,jy1,k)+uo(im,jy1,k))**2.+
     %(wo(is,jy1,k)+wo(is,jy1,km))**2.))    !*ros
!	    endif
  700	    continue
	    CALL MPI_SEND
     %(recvar(jy1:jy2,kz1:kz2),(ny-2)*(nz-2),MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
      do kkc=1,nsez2max
	 ks=ksez2(kkc)
	 km=ks-1
	 kk=isezindx2+(kkc-1)
	 kkd=kk/10
	 kku=mod(kk,10)
	 open(83,file='results2_2D_'//char(48+iplant2c)
     %//char(48+iplant2d)//char(48+iplant2u)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_vort2D.q')
	 write(83,100)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	 i1=1
	 write(83,101)ntime,i1,i1,i1
!!!!!!	       write(83,*)
	 write(83,102)
     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
c
	 do 140 i=ix1,ix2,nxp
	 do 141 l=1,lamb2
	 do 141 j=jy1,jy2,nyp
!	 if(all(flagpo(i,j,ks,:)).le.0) then
	    write(83,9)voraz2(i,j,kkc)
!	 else
!	    write(83,9)1.e+03
!	 endif
  141	 continue
!	 if(all(flagpo(i,jy1,ks,:)).le.0) then
	    write(83,9)voraz2(i,jy1,kkc)
!	 else
!	    write(83,9)1.e+03
!	 endif
  140	 continue
c
	 do 150 i=ix1,ix2,nxp
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
!	 if(all(flagpo(i,j,ks,:)).le.0) then
	    write(83,9)vorr2(i,j,kkc)
!	 else
!	    write(83,9)1.e+03
!	 endif
  151	 continue
!	 if(all(flagpo(i,jy1,ks,:)).le.0) then
	    write(83,9)vorr2(i,jy1,kkc)
!	 else
!	    write(83,9)1.e+03
!	 endif
  150	 continue
c
	 do 160 i=ix1,ix2,nxp
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2,nyp
!	 if(all(flagpo(i,j,ks,:)).le.0) then
	    write(83,9)vorz2(i,j,kkc)
!	 else
!	    write(83,9)1.e+03
!	 endif
  161	 continue
!	 if(all(flagpo(i,jy1,ks,:)).le.0) then
	    write(83,9)vorz2(i,jy1,kkc)
!	 else
!	    write(83,9)1.e+03
!	 endif
  160	 continue
c
	 do 170 i=ix1,ix2,nxp
	 im=i-1
	 do 171 l=1,lamb2
	 do 171 j=jy1,jy2,nyp
!	 if(all(flagpo(i,j,ks,:)).le.0) then
	    jm=jmv(j)
	    write(83,8)
     %(p(i,j,ks)+0.5*0.25*((vo(i,j,ks)+vo(i,jm,ks))**2.+
     %(uo(i,j,ks)+uo(im,j,ks))**2.+
     %(wo(i,j,ks)+wo(i,j,km))**2.))    !*ros
!	 else
!	    write(83,8)1.e+03
!	 endif
  171	 continue
!	 if(all(flagpo(i,jy1,ks,:)).le.0) then
	    write(83,8)
     %(p(i,jy1,ks)+0.5*0.25*((vo(i,jy1,ks)+vo(i,jy2,ks))**2.+
     %(uo(i,jy1,ks)+uo(im,jy1,ks))**2.+
     %(wo(i,jy1,ks)+wo(i,jy1,km))**2.))	   !*ros
!	 else
!	    write(83,8)1.e+03
!	 endif
  170	 continue
c
	 close(83)
      enddo
c
      do jjc=1,nsez14
	 js=jsez14(jjc)
	 jsy=jsym(js)
	 jm=jmv(js)
	 jmsy=jmv(jsy)
	 if(myrank.eq.0) then
	    jjcd=jjc/10
	    jjcu=mod(jjc,10)
	    open(83,file='results3_2D_'//char(48+iplant2c)
     %//char(48+iplant2d)//char(48+iplant2u)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_vort2D.q')
	    write(83,100)1,2*(nx-2)/nxp,(nzg-2)/nzp
	    i1=1
	    write(83,101)ntime,i1,i1,i1
!!!!!!		  write(83,*)
	    write(83,102)
     %(((i1,j=1,1),i=1,2*(nx-2)/nxp),k=1,(nzg-2)/nzp)
c
	    do 240 k=kz1,kz2,nzp
	    do 241 i=ix2,ix1,-nxp
!	    if(all(flagpo(i,jsy,k,:)).le.0) then
	       write(83,9)voraz14(i,nsez14+jjc,k)
!	    else
!	       write(83,9)1.e+03
!	    endif
  241	    continue
	    do 242 i=ix1,ix2,nxp
!	    if(all(flagpo(i,js,k,:)).le.0) then
	       write(83,9)voraz14(i,jjc,k)
!	    else
!	       write(83,9)1.e+03
!	    endif
  242	    continue
  240	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
!		  kk1=kz1+jjp*(nz-2)
!		  kk2=kz2+jjp*(nz-2)
!		  CALL MPI_RECV
!    %(flagpor2(ix1:ix2,kk1:kk2,:),
!    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,6,
!    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,7,
     %MPI_COMM_EDDY,STATUS,IERR)
!		  CALL MPI_RECV
!    %(flagpor2sym(ix1:ix2,kk1:kk2,:),
!    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,8,
!    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 243 kr=kz1,kz2,nzp
!		  krg=kr+jjp*(nz-2)
		  do 244 ir=ix2,ix1,-nxp
!		  if(all(flagpor2sym(ir,krg,:)).le.0) then
		     write(83,9)recvar2sym(ir,kr)
!		  else
!		     write(83,9)1.e+03
!		  endif
  244		  continue
		  do 245 ir=ix1,ix2,nxp
!		  if(all(flagpor2(ir,krg,:)).le.0) then
		     write(83,9)recvar2(ir,kr)
!		  else
!		     write(83,9)1.e+03
!		  endif
  245		  continue
  243		  continue
	       enddo
	    endif
c
	    do 250 k=kz1,kz2,nzp
	    do 251 i=ix2,ix1,-nxp
!	    if(all(flagpo(i,jsy,k,:)).le.0) then
	       write(83,9)vorr14(i,nsez14+jjc,k)
!	    else
!	       write(83,9)1.e+03
!	    endif
  251	    continue
	    do 252 i=ix1,ix2,nxp
!	    if(all(flagpo(i,js,k,:)).le.0) then
	       write(83,9)vorr14(i,jjc,k)
!	    else
!	       write(83,9)1.e+03
!	    endif
  252	    continue
  250	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 253 kr=kz1,kz2,nzp
!		  krg=kr+jjp*(nz-2)
		  do 254 ir=ix2,ix1,-nxp
!		  if(all(flagpor2sym(ir,krg,:)).le.0) then
		     write(83,9)recvar2sym(ir,kr)
!		  else
!		     write(83,9)1.e+03
!		  endif
  254		  continue
		  do 255 ir=ix1,ix2,nxp
!		  if(all(flagpor2(ir,krg,:)).le.0) then
		     write(83,9)recvar2(ir,kr)
!		  else
!		     write(83,9)1.e+03
!		  endif
  255		  continue
  253		  continue
	       enddo
	    endif
c
	    do 260 k=kz1,kz2,nzp
	    do 261 i=ix2,ix1,-nxp
!	    if(all(flagpo(i,jsy,k,:)).le.0) then
	       write(83,9)vorz14(i,nsez14+jjc,k)
!	    else
!	       write(83,9)1.e+03
!	    endif
  261	    continue
	    do 262 i=ix1,ix2,nxp
!	    if(all(flagpo(i,js,k,:)).le.0) then
	       write(83,9)vorz14(i,jjc,k)
!	    else
!	       write(83,9)1.e+03
!	    endif
  262	    continue
  260	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 263 kr=kz1,kz2,nzp
!		  krg=kr+jjp*(nz-2)
		  do 264 ir=ix2,ix1,-nxp
!		  if(all(flagpor2sym(ir,krg,:)).le.0) then
		     write(83,9)recvar2sym(ir,kr)
!		  else
!		     write(83,9)1.e+03
!		  endif
  264		  continue
		  do 265 ir=ix1,ix2,nxp
!		  if(all(flagpor2(ir,krg,:)).le.0) then
		     write(83,9)recvar2(ir,kr)
!		  else
!		     write(83,9)1.e+03
!		  endif
  265		  continue
  263		  continue
	       enddo
	    endif
c
	    do 270 k=kz1,kz2,nzp
	    km=k-1
	    do 271 i=ix2,ix1,-nxp
!	    if(all(flagpo(i,jsy,k,:)).le.0) then
	       im=i-1
	       write(83,8)
     %(p(i,jsy,k)+0.5*0.25*((vo(i,jsy,k)+vo(i,jmsy,k))**2.+
     %(uo(i,jsy,k)+uo(im,jsy,k))**2.+
     %(wo(i,jsy,k)+wo(i,jsy,km))**2.))	   !*ros
!	    else
!	       write(83,8)1.e+03
!	    endif
  271	    continue
	    do 272 i=ix1,ix2,nxp
!	    if(all(flagpo(i,js,k,:)).le.0) then
	       im=i-1
	       write(83,8)
     %(p(i,js,k)+0.5*0.25*((vo(i,js,k)+vo(i,jm,k))**2.+
     %(uo(i,js,k)+uo(im,js,k))**2.+
     %(wo(i,js,k)+wo(i,js,km))**2.))	 !*ros
!	    else
!	       write(83,8)1.e+03
!	    endif
  272	    continue
  270	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,jjp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 273 kr=kz1,kz2,nzp
!		  krg=kr+jjp*(nz-2)
		  do 274 ir=ix2,ix1,-nxp
!		  if(all(flagpor2sym(ir,krg,:)).le.0) then
		     write(83,8)recvar2sym(ir,kr)
!		  else
!		     write(83,8)1.e+03
!		  endif
  274		  continue
		  do 275 ir=ix1,ix2,nxp
!		  if(all(flagpor2(ir,krg,:)).le.0) then
		     write(83,8)recvar2(ir,kr)
!		  else
!		     write(83,8)1.e+03
!		  endif
  275		  continue
  273		  continue
	       enddo
	    endif
	    close(83)
	 else
c	    CALL MPI_SEND
c    %(flagpo(ix1:ix2,js,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,6,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(voraz14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,7,
     %MPI_COMM_EDDY,STATUS,IERR)
c	    CALL MPI_SEND
c    %(flagpo(ix1:ix2,jsy,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,8,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(voraz14(ix1:ix2,nsez14+jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorr14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorr14(ix1:ix2,nsez14+jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorz14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorz14(ix1:ix2,nsez14+jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	    do 701 k=kz1,kz2,nzp
	    km=k-1
	    do 711 i=ix2,ix1,-nxp
!	    if(all(flagpo(i,jsy,k,:)).le.0) then
	       im=i-1
	       recvar2sym(i,k)=
     %(p(i,jsy,k)+0.5*0.25*((vo(i,jsy,k)+vo(i,jmsy,k))**2.+
     %(uo(i,jsy,k)+uo(im,jsy,k))**2.+
     %(wo(i,jsy,k)+wo(i,jsy,km))**2.))	   !*ros
!	    endif
  711	    continue
	    do 721 i=ix1,ix2,nxp
!	    if(all(flagpo(i,js,k,:)).le.0) then
	       im=i-1
	       recvar2(i,k)=
     %(p(i,js,k)+0.5*0.25*((vo(i,js,k)+vo(i,jm,k))**2.+
     %(uo(i,js,k)+uo(im,js,k))**2.+
     %(wo(i,js,k)+wo(i,js,km))**2.))	 !*ros
!	    endif
  721	    continue
  701	    continue
	    CALL MPI_SEND
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*(nz-2),MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine plotinststagg(iplant9c,iplant9d,iplant9p,
     %nx,ny,nz,nzg,nbd,ntime,xc,xu,yc,yv,zcg,zwg,p,vo,uo,wo)
!     %flaguo,flagvo,flagwo)
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
     %jjc,jjcd,jjcu,js,jm,jp,jsy,jpsy,i,j,k,l,ir,jr,kr,lr,
     %i1,jjp,kk1,kk2,kk,kkd,kku
      integer STATUS(MPI_STATUS_SIZE)
!      integer flaguor(ny,nzg,nbd),
!     %flagvor2(nx,nzg,nbd),flagvor2sym(nx,nzg,nbd)
      real pstagg,vstagg,ustagg,wstagg
      real pr(ny,nzg),vor(ny,nzg),uor(ny,nzg),wor(ny,nzg),
     %pr2(nx,nzg),pr2sym(nx,nzg),vor2(nx,nzg),vor2sym(nx,nzg),
     %uor2(nx,nzg),uor2sym(nx,nzg),wor2(nx,nzg),wor2sym(nx,nzg)
c
c Parameters
      integer lamb2,nxp,nyp,nzp
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
      if(iplant9.eq.1) then
	 if(myrank.eq.0) then
	    do iic=1,nsez1
	       iicd=iic/10
	       iicu=mod(iic,10)
	       open
     %(81,file='grid1_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'_plotinst.xyz')
	       open
     %(82,file='grid1_new_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'_plotinst.xyz')
	       is=isez1(iic)
	       write(81,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	       write(82,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	       do 1 k=2,nzg-1,nzp
	       do 10 l=1,lamb2
	       do 10 j=jy1,jy2,nyp
	       write(81,9)xu(is)*cos(yc(jy1))
   10	       write(82,9)xu(is)*cos(yc(j)+(l-1)*2.*pi/lamb2)
	       write(81,9)xu(is)*cos(yc(jy1))
   1	       write(82,9)xu(is)*cos(yc(jy1)+lamb2*deltasect)
	       do 2 k=2,nzg-1,nzp
	       do 20 l=1,lamb2
	       do 20 j=jy1,jy2,nyp
	       write(81,9)
     %xu(is)*sin(yc(jy1))+xu(is)*(yc(j)+(l-1)*2.*pi/lamb2)
   20	       write(82,9)xu(is)*sin(yc(j)+(l-1)*2.*pi/lamb2)
	       write(81,9)
     %xu(is)*sin(yc(jy1))+xu(is)*deltasect*lamb2
   2	       write(82,9)xu(is)*sin(yc(jy1)+lamb2*deltasect)
	       do 3 k=2,nzg-1,nzp
	       do 30 l=1,lamb2
	       do 30 j=jy1,jy2,nyp
	       write(81,9)zcg(k)
   30	       write(82,9)zcg(k)
	       write(81,9)zcg(k)
   3	       write(82,9)zcg(k)
	       close(81)
	       close(82)
	    enddo
	    do kkc=1,nsez2
	       kkcd=kkc/10
	       kkcu=mod(kkc,10)
	       open(81,file='grid2_2D_slice'//char(48+kkcd)
     %//char(48+kkcu)//'_plotinst.xyz')
	       ks=ksez2g(kkc)
	       write(81,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	       do 11 i=ix1,ix2,nxp
	       do 110 l=1,lamb2
	       do 110 j=jy1,jy2,nyp
  110	       write(81,9)xc(i)*cos(yc(j)+(l-1)*2.*pi/lamb2)
  11	       write(81,9)xc(i)*cos(yc(jy1)+lamb2*deltasect)
	       do 12 i=ix1,ix2,nxp
	       do 120 l=1,lamb2
	       do 120 j=jy1,jy2,nyp
  120	       write(81,9)xc(i)*sin(yc(j)+(l-1)*2.*pi/lamb2)
  12	       write(81,9)xc(i)*sin(yc(jy1)+lamb2*deltasect)
	       do 13 i=ix1,ix2,nxp
	       do 130 l=1,lamb2
	       do 130 j=jy1,jy2,nyp
  130	       write(81,9)zwg(ks)
  13	       write(81,9)zwg(ks)
	       close(81)
	    enddo
	    do jjc=1,nsez14
	       jjcd=jjc/10
	       jjcu=mod(jjc,10)
	       open(81,file='grid3_2D_slice'//char(48+jjcd)
     %//char(48+jjcu)//'_plotinst.xyz')
	       js=jsez14(jjc)
	       jsy=jsym(js)
	       write(81,101)1,2*(nx-2)/nxp,(nzg-2)/nzp
	       do 210 k=2,nzg-1,nzp
	       do 211 i=ix2,ix1,-nxp
  211	       write(81,9)xc(i)*cos(yv(jsy))
	       do 212 i=ix1,ix2,nxp
  212	       write(81,9)xc(i)*cos(yv(js))
  210	       continue
	       do 220 k=2,nzg-1,nzp
	       do 221 i=ix2,ix1,-nxp
  221	       write(81,9)xc(i)*sin(yv(jsy))
	       do 222 i=ix1,ix2,nxp
  222	       write(81,9)xc(i)*sin(yv(js))
  220	       continue
	       do 230 k=2,nzg-1,nzp
	       do 230 i=1,2*(nx-2),nxp
  230	       write(81,9)zcg(k)
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
	    write(82,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	    i1=1
	    write(82,102)ntime,i1,i1,i1
!
!Stagnation pressure
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
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       pstagg=
     %0.5*(p(is,jy1,k)+p(ip,jy1,k))
	       vstagg=
     %0.25*(vo(is,jy1,k)+vo(ip,jy1,k)+vo(is,jy2,k)+vo(ip,jy2,k))
	       ustagg=
     %uo(is,jy1,k)
	       wstagg=
     %0.25*(wo(is,jy1,k)+wo(ip,jy1,k)+wo(is,jy1,km)+wo(ip,jy1,km))
	       write(82,9)
     %(pstagg+0.5*(vstagg**2.+ustagg**2.+
     %wstagg**2.))     !*ros
c	    else
c	       write(82,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,9)
     %(pr(jy1,kr)+0.5*(vor(jy1,kr)**2.+uor(jy1,kr)**2.+
     %wor(jy1,kr)**2.))	    !*ros
c		  else
c		     write(82,9)1.e+03
c		  endif
 42		  continue
	       enddo
	    endif
!	    write(82,*)	  !!!!!!
!	    write(82,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
!
!Azimuthal velocity
!
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       jm=jmv(j)
	       write(82,9)
     %0.25*(vo(is,j,k)+vo(ip,j,k)+vo(is,jm,k)+vo(ip,jm,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
   51	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)
     %0.25*(vo(is,jy1,k)+vo(ip,jy1,k)+vo(is,jy2,k)+vo(ip,jy2,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
   50	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 52 kr=kk1,kk2,nzp
		  do 53 lr=1,lamb2
		  do 53 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)vor(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   53		  continue
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,9)vor(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   52		  continue
	       enddo
	    endif
!
!Radial velocity
!
	    do 60 k=kz1,kz2,nzp
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)uo(is,j,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
   61	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)uo(is,jy1,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
   60	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 62 kr=kk1,kk2,nzp
		  do 63 lr=1,lamb2
		  do 63 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)uor(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   63		  continue
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,9)uor(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
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
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)
     %0.25*(wo(is,jy1,k)+wo(ip,jy1,k)+wo(is,jy1,km)+wo(ip,jy1,km))
c	    else
c	       write(82,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,9)wor(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
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
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,8)0.5*(p(is,jy1,k)+p(ip,jy1,k))	  !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,8)pr(jy1,kr)   !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
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
      do kkc=1,nsez2max
	 ks=ksez2(kkc)
	 kp=ks+1
	 kk=isezindx2+(kkc-1)
	 kkd=kk/10
	 kku=mod(kk,10)
	 open(82,file='results2_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_plotinst.q')
	 write(82,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
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
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    pstagg=
     %0.5*(p(i,jy1,ks)+p(i,jy1,kp))
	    vstagg=
     %0.25*(vo(i,jy1,ks)+vo(i,jy2,ks)+vo(i,jy1,kp)+vo(i,jy2,kp))
	    ustagg=
     %0.25*(uo(i,jy1,ks)+uo(im,jy1,ks)+uo(i,jy1,kp)+uo(im,jy1,kp))
	    wstagg=
     %wo(i,jy1,ks)
	    write(82,9)
     %(pstagg+0.5*(vstagg**2.+ustagg**2.+wstagg**2.))	 !*ros
c	 else
c	    write(82,9)1.e+03
c	 endif
 140	 continue
!	 write(82,*)	!!!!!!
!	 write(82,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
!
!Azimuthal velocity
!
	 do 150 i=ix1,ix2,nxp
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    jm=jmv(j)
	    write(82,9)
     %0.25*(vo(i,j,ks)+vo(i,jm,ks)+vo(i,j,kp)+vo(i,jm,kp))
c	 else
c	    write(82,9)1.e+03
c	 endif
 151	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)
     %0.25*(vo(i,jy1,ks)+vo(i,jy2,ks)+vo(i,jy1,kp)+vo(i,jy2,kp))
c	 else
c	    write(82,9)1.e+03
c	 endif
 150	 continue
!
!Radial velocity
!
	 do 160 i=ix1,ix2,nxp
	 im=i-1
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)
     %0.25*(uo(i,j,ks)+uo(im,j,ks)+uo(i,j,kp)+uo(im,j,kp))
c	 else
c	    write(82,9)1.e+03
c	 endif
 161	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)
     %0.25*(uo(i,jy1,ks)+uo(im,jy1,ks)+uo(i,jy1,kp)+uo(im,jy1,kp))
c	 else
c	    write(82,9)1.e+03
c	 endif
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
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)wo(i,jy1,ks)
c	 else
c	    write(82,9)1.e+03
c	 endif
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
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,8)0.5*(p(i,jy1,ks)+p(i,jy1,kp))   !*ros
c	 else
c	    write(82,8)1.e+03
c	 endif
 180	 continue
	 close(82)
      enddo
c
      do jjc=1,nsez14
	 js=jsez14(jjc)
	 jsy=jsym(js)
	 jp=jpv(js)
	 jpsy=jpv(jsy)
	 if(myrank.eq.0) then
	    jjcd=jjc/10
	    jjcu=mod(jjc,10)
	    open(82,file='results3_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_plotinst.q')
	    write(82,101)1,2*(nx-2)/nxp,(nzg-2)/nzp
	    i1=1
	    write(82,102)ntime,i1,i1,i1
!
!Stagnation pressure
!
	    do 240 k=kz1,kz2,nzp
	    km=k-1
	    do 241 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       im=i-1
	       pstagg=
     %0.5*(p(i,jsy,k)+p(i,jpsy,k))
	       vstagg=
     %vo(i,jsy,k)
	       ustagg=
     %0.25*(uo(i,jsy,k)+uo(im,jsy,k)+uo(i,jpsy,k)+uo(im,jpsy,k))
	       wstagg=
     %0.25*(wo(i,jsy,k)+wo(i,jpsy,k)+wo(i,jsy,km)+wo(i,jpsy,km))
	       write(82,9)
     %(pstagg+0.5*(vstagg**2.+ustagg**2.+wstagg**2.))	 !*ros
c	    else
c	       write(82,9)1.e+03
c	    endif
 241	    continue
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
c		  CALL MPI_RECV
c    %(flagvor2sym(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,11,
c    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(pr2sym(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(vor2sym(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(uor2sym(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(wor2sym(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 243 kr=kk1,kk2,nzp
		  do 244 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,9)
     %(pr2sym(ir,kr)+0.5*(vor2sym(ir,kr)**2.+uor2sym(ir,kr)**2.
     %+wor2sym(ir,kr)**2.))    !*ros
c		  else
c		     write(82,9)1.e+03
c		  endif
 244		  continue
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
!
!Azimuthal velocity
!
	    do 250 k=kz1,kz2,nzp
	    do 251 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,9)vo(i,jsy,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 251	    continue
	    do 252 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,9)vo(i,js,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 252	    continue
 250	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 253 kr=kk1,kk2,nzp
		  do 254 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,9)vor2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 254		  continue
		  do 255 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)vor2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 255		  continue
 253		  continue
	       enddo
	    endif
!
!Radial velocity
!
	    do 260 k=kz1,kz2,nzp
	    do 261 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       im=i-1
	       write(82,9)
     %0.25*(uo(i,jsy,k)+uo(im,jsy,k)+uo(i,jpsy,k)+uo(im,jpsy,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
 261	    continue
	    do 262 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       im=i-1
	       write(82,9)
     %0.25*(uo(i,js,k)+uo(im,js,k)+uo(i,jp,k)+uo(im,jp,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
 262	    continue
 260	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 263 kr=kk1,kk2,nzp
		  do 264 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,9)uor2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 264		  continue
		  do 265 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)uor2(ir,kr)
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
	    do 271 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,9)
     %0.25*(wo(i,jsy,k)+wo(i,jpsy,k)+wo(i,jsy,km)+wo(i,jpsy,km))
c	    else
c	       write(82,9)1.e+03
c	    endif
 271	    continue
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
		  do 274 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,9)wor2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 274		  continue
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
	    do 281 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,8)0.5*(p(i,jsy,k)+p(i,jpsy,k))	  !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
 281	    continue
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
		  do 284 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,8)pr2sym(ir,kr)	 !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
 284		  continue
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
c	    CALL MPI_SEND
c    %(flagvo(ix1:ix2,jsy,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,11,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((p(ix1:ix2,jsy,kz1:kz2)+p(ix1:ix2,jpsy,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vo(ix1:ix2,jsy,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((uo(1:ix2-1,jsy,kz1:kz2)+uo(ix1:ix2,jsy,kz1:kz2)+
     %uo(1:ix2-1,jpsy,kz1:kz2)+uo(ix1:ix2,jpsy,kz1:kz2))*0.25,
     %(nx-2)*1*(nz-2),MTYPE,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((wo(ix1:ix2,jsy,1:kz2-1)+wo(ix1:ix2,jpsy,1:kz2-1)+
     %wo(ix1:ix2,jsy,kz1:kz2)+wo(ix1:ix2,jpsy,kz1:kz2))*0.25,
     %(nx-2)*1*(nz-2),MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vort2stagg(nx,ny,nz,xc,xu,zc,uo,vo,wo)
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
      ! Azimuthal vorticity - Circumferential sections !
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
      ! Azimuthal vorticity - Cross sections !
      do 111 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
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
	do 113 ic=ix1+1,ix2
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
      ! Azimuthal vorticity - Meridian sections !
      do 121 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jsy=jsym(jc)
	jpsy=jsym(jp)
	do 121 kc=kz1,kz2
	  kp=kc+1
	  km=kc-1
	  dz=zc(kp)-zc(km)
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
	  dq1x3=(0.25*(uo(ic,jsy,kp)+uo(im,jsy,kp)+
     %uo(ic,jpsy,kp)+uo(im,jpsy,kp))
     %-0.25*(uo(ic,jsy,km)+uo(im,jsy,km)+
     %uo(ic,jpsy,km)+uo(im,jpsy,km)))/dz
	  dq3x1=(0.25*(wo(ip,jsy,kc)+wo(ip,jpsy,kc)+
     %wo(ip,jsy,km)+wo(ip,jpsy,km))
     %-0.25*(wo(ic,jc,kc)+wo(ic,jp,kc)+
     %wo(ic,jc,km)+wo(ic,jp,km)))/dr
	  voraz14(ic,nsez14+jjc,kc)=dq1x3-dq3x1
	  do 121 ic=ix1+1,ix2
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
	    dq1x3=(0.25*(uo(ic,jsy,kp)+uo(im,jsy,kp)+
     %uo(ic,jpsy,kp)+uo(im,jpsy,kp))
     %-0.25*(uo(ic,jsy,km)+uo(im,jsy,km)+
     %uo(ic,jpsy,km)+uo(im,jpsy,km)))/dz
	    dq3x1=(0.25*(wo(ip,jsy,kc)+wo(ip,jpsy,kc)+
     %wo(ip,jsy,km)+wo(ip,jpsy,km))
     %-0.25*(wo(im,jsy,kc)+wo(im,jpsy,kc)+
     %wo(im,jsy,km)+wo(im,jpsy,km)))/dr
	    voraz14(ic,nsez14+jjc,kc)=dq1x3-dq3x1
  121 continue
c
      ! Radial vorticity - Circumferential sections !
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
	    vorr1(iic,jc,kc)=dq3x2/xu(ic)-dq2x3
  131 continue
c
      ! Radial vorticity - Cross sections !
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
	    vorr2(ic,jc,kkc)=dq3x2/xc(ic)-dq2x3
  141 continue
c
      ! Radial vorticity - Meridian sections !
      do 151 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jsy=jsym(jc)
	jpsy=jpv(jsy)
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
	    vorr14(ic,jjc,kc)=dq3x2/xc(ic)-dq2x3
	    dq3x2=(0.5*(wo(ic,jpsy,kc)+wo(ic,jpsy,km))
     %-0.5*(wo(ic,jsy,kc)+wo(ic,jsy,km)))/dtheta   !
	    dq2x3=(vo(ic,jsy,kp)
     %-vo(ic,jsy,km))/dz
	    vorr14(ic,nsez14+jjc,kc)=dq3x2/xc(ic)-dq2x3
  151 continue
c
      ! Axial vorticity - Circumferential sections !
      do 161 kc=kz1,kz2
	do 161 iic=1,nsez1
	  ic=isez1(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 161 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.5*(vo(ip,jc,kc)+vo(ip,jm,kc))*xc(ip)
     %-0.5*(vo(ic,jc,kc)+vo(ic,jm,kc))*xc(ic))/dr
	    dq1x2=(uo(ic,jp,kc)
     %-uo(ic,jm,kc))/dtheta    !
	    vorz1(iic,jc,kc)=(dq2x1-dq1x2)/xu(ic)
  161 continue
c
      ! Axial vorticity - Cross sections !
      do 171 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
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
     %+vo(ic,jm,kp)+vo(ip,jm,kp))*xu(ic))/dr
	  dq1x2=(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)+
     %uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jm,kc)+uo(im,jm,kc)+
     %uo(ic,jm,kp)+uo(im,jm,kp)))/dtheta    !
	  vorz2(ic,jc,kkc)=(dq2x1-dq1x2)/xc(ic)
  172	continue
c
	do 173 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 173 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.25*(vo(ip,jc,kc)+vo(ip,jm,kc)+
     %vo(ip,jc,kp)+vo(ip,jm,kp))*xc(ip)
     %-0.25*(vo(im,jc,kc)+vo(im,jm,kc)+
     %vo(im,jc,kp)+vo(im,jm,kp))*xc(im))/dr
	    dq1x2=
     %(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)+
     %uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jm,kc)+uo(im,jm,kc)+
     %uo(ic,jm,kp)+uo(im,jm,kp)))/dtheta   !
	    vorz2(ic,jc,kkc)=(dq2x1-dq1x2)/xc(ic)
  173 continue
  171 continue
c
      ! Axial vorticity - Meridian sections !
      do 181 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jsy=jsym(jc)
	jpsy=jsym(jp)
	dtheta=dely
	do 181 kc=kz1,kz2
	  ic=ix1
	  im=ic-1
	  ip=ic+1
	  dr=xu(ic)-xu(im)
	  dq2x1=(0.5*(vo(ip,jc,kc)+vo(ic,jc,kc))*xu(ic))/dr
	  dq1x2=(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dtheta	  !
	  vorz14(ic,jjc,kc)=(dq2x1-dq1x2)/xc(ic)
	  dq2x1=(0.5*(vo(ip,jsy,kc)+vo(ic,jsy,kc))*xu(ic))/dr
	  dq1x2=(0.5*(uo(ic,jpsy,kc)+uo(im,jpsy,kc))
     %-0.5*(uo(ic,jsy,kc)+uo(im,jsy,kc)))/dtheta   !
	  vorz14(ic,nsez14+jjc,kc)=(dq2x1-dq1x2)/xc(ic)
	  do 181 ic=ix1+1,ix2
	    ip=ic+1
	    im=ic-1
	    dr=xc(ip)-xc(im)
	    dq2x1=(vo(ip,jc,kc)*xc(ip)
     %-vo(im,jc,kc)*xc(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dtheta	  !
	    vorz14(ic,jjc,kc)=(dq2x1-dq1x2)/xc(ic)
	    dq2x1=(vo(ip,jsy,kc)*xc(ip)
     %-vo(im,jsy,kc)*xc(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jpsy,kc)+uo(im,jpsy,kc))
     %-0.5*(uo(ic,jsy,kc)+uo(im,jsy,kc)))/dtheta   !
	    vorz14(ic,nsez14+jjc,kc)=(dq2x1-dq1x2)/xc(ic)
  181 continue
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine plotvortstagg(iplant9c,iplant9d,iplant9p,
     %nx,ny,nz,nzg,nbd,ntime,flaguo,flagwo,flagvo,p,uo,vo,wo,tv)
c
      use sections
      use vorticity
      include'common.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer iplant9c,iplant9d,iplant9p,nx,ny,nz,nzg,nbd,ntime
      integer flaguo(nx,ny,nz,nbd),flagwo(nx,ny,nz,nbd),
     %flagvo(nx,ny,nz,nbd)
      real p(nx,ny,nz),uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),
     %tv(nx,ny,nz)
c
c Local variables
      integer iic,iicd,iicu,is,ip,jjc,jjcd,jjcu,js,jp,
     %jsy,jpsy,kkc,ks,kp,i,j,k,l,ir,jr,kr,lr,i1,jjp,kk1,kk2,
     %kk,kkd,kku,krg
      integer STATUS(MPI_STATUS_SIZE)
!      integer flaguor(ny,nzg,nbd),
!     %flagvor2(nx,nzg,nbd),flagvor2sym(nx,nzg,nbd)
      real recvar(ny,nz),recvar2(nx,nz),recvar2sym(nx,nz)
c
c Parameters
      integer lamb2,nxp,nyp,nzp
      parameter (lamb2=1)
      parameter (nyp=1,nxp=1,nzp=1)
c
   8  format(10f25.8)
   9  format(10f25.8)
 101  format(3i8)
 102  format(i8,3i2)
 103  format(i2)
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
	    write(83,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	    write(83,102)ntime,i1,i1,i1
!!!!!!		  write(83,*)
	    write(83,103)
     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
	    do 40 k=kz1,kz2,nzp
	    do 41 l=1,lamb2
	    do 41 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)voraz1(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
   41	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)voraz1(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)recvar(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
   42		  continue
	       enddo
	    endif
c
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)vorr1(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
   51	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)vorr1(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)recvar(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
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
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)vorz1(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)recvar(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
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
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)0.5*(tv(is,jy1,k)+tv(ip,jy1,k))/ru1
c	    else
c	       write(83,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)recvar(jy1,kr)/ru1
c		  else
c		     write(83,9)1.e+03
c		  endif
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
     %(voraz1(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorr1(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,3,
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
	 write(83,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	 write(83,102)ntime,i1,i1,i1
!!!!!!	       write(83,*)
	 write(83,103)
     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
	 do 140 i=ix1,ix2,nxp
	 do 141 l=1,lamb2
	 do 141 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)voraz2(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
  141	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)voraz2(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
  140	 continue
	 do 150 i=ix1,ix2,nxp
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)vorr2(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
  151	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)vorr2(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
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
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)vorz2(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
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
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)0.5*(tv(i,jy1,ks)+tv(i,jy1,kp))/ru1
c	 else
c	    write(83,9)1.e+03
c	 endif
  170	 continue
	 close(83)
      enddo
c
      do jjc=1,nsez14
	 js=jsez14(jjc)
	 jsy=jsym(js)
	 jp=jpv(js)
	 jpsy=jpv(jsy)
	 if(myrank.eq.0) then
	    jjcd=jjc/10
	    jjcu=mod(jjc,10)
	    open(83,file='results3_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_vort.q')
	    write(83,101)1,2*(nx-2)/nxp,(nzg-2)/nzp
	    write(83,102)ntime,i1,i1,i1
!!!!!!		  write(83,*)
	    write(83,103)
     %(((i1,j=1,1),i=1,2*(nx-2)/nxp),k=1,(nzg-2)/nzp)
	    do 240 k=kz1,kz2,nzp
	    do 241 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(83,9)voraz14(i,nsez14+jjc,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
  241	    continue
	    do 242 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(83,9)voraz14(i,jjc,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
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
c		  CALL MPI_RECV
c    %(flagvor2sym(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,7,
c    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 243 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 244 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(83,9)recvar2sym(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
  244		  continue
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
	    do 251 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(83,9)vorr14(i,nsez14+jjc,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
  251	    continue
	    do 252 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(83,9)vorr14(i,jjc,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
  252	    continue
  250	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 253 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 254 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(83,9)recvar2sym(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
  254		  continue
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
	    do 261 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(83,9)vorz14(i,nsez14+jjc,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
  261	    continue
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
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 263 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 264 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(83,9)recvar2sym(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
  264		  continue
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
	    do 271 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(83,9)0.5*(tv(i,jsy,k)+tv(i,jpsy,k))/ru1
c	    else
c	       write(83,9)1.e+03
c	    endif
  271	    continue
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
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 273 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 274 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(83,9)recvar2sym(ir,kr)/ru1
c		  else
c		     write(83,9)1.e+03
c		  endif
  274		  continue
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
c	    CALL MPI_SEND
c    %(flagvo(ix1:ix2,jsy,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,7,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(voraz14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(voraz14(ix1:ix2,nsez14+jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorr14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorr14(ix1:ix2,nsez14+jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorz14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorz14(ix1:ix2,nsez14+jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((tv(ix1:ix2,js,kz1:kz2)+tv(ix1:ix2,jp,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,14,MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((tv(ix1:ix2,jsy,kz1:kz2)+tv(ix1:ix2,jpsy,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,15,MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine mediecalcstaggavg
     %(nx,ny,nz,nzg,nbd,time,dt,vo,uo,wo,p,tv,xc,xu,yc,yv,zc,
     %zcg,zwg,flaguo,flagvo,flagwo,ntime)
!    %mbd,nfacet,unvect,vertex)
c
      use sections
      use mediecalc
      use mediecalc1
      use mediecalc2
      use mediecalc3
      use vortavg
      include'common.h'
      include'averages.h'
!     include'immersed.h'
c
c Global variables
      integer nx,ny,nz,nzg,nbd,ntime
      integer flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),
     %flagwo(nx,ny,nz,nbd)
      real time,dt
      real vo(nx,ny,nz),uo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz),
     %xc(nx),xu(nx),yc(ny),yv(ny),zc(nz),zcg(nzg),zwg(nzg),
     %tv(nx,ny,nz)
!     integer mbd,nfacet
!     real unvect(3,nfacet),vertex(3,3,nfacet)
c
c Local variables
      integer iic,is,ip,im,i,jjc,jp,jm,j,kkc,ks,kp,km,k,
     %iplant4c,iplant4d,iplant4p,iplant5c,iplant5d,iplant5p
      real vstagg,ustagg,wstagg,pstagg,vmod,tvstagg
!     integer ibd
c
      tempo3=tempo3+dt
!
!     Time averaged values and root mean squares
!     at some circumferential sections
!
      do 101 iic=1,nsez1
      is=isez1(iic)
      ip=is+1
      do 100 k=kz1,kz2
      km=k-1
      do 100 j=jy1,jy2
      jm=jmv(j)
      tempo4(iic,j,k)=tempo4(iic,j,k)+dt
      vstagg=(vo(is,j,k)+vo(ip,j,k)+
     %vo(is,jm,k)+vo(ip,jm,k))*0.25
      vplot(iic,j,k)=vplot(iic,j,k)+
     %vstagg*dt
      ustagg=uo(is,j,k)
      uplot(iic,j,k)=uplot(iic,j,k)+
     %ustagg*dt
      wstagg=(wo(is,j,k)+wo(ip,j,k)+
     %wo(is,j,km)+wo(ip,j,km))*0.25
      wplot(iic,j,k)=wplot(iic,j,k)+
     %wstagg*dt
      pstagg=(p(is,j,k)+p(ip,j,k))*0.5
      tvstagg=(tv(is,j,k)+tv(ip,j,k))*0.5
      pplot(iic,j,k)=pplot(iic,j,k)+
     %pstagg*dt
      vrmsplot(iic,j,k)=vrmsplot(iic,j,k)+
     %dt*(vstagg**2.)
      urmsplot(iic,j,k)=urmsplot(iic,j,k)+
     %dt*(ustagg**2.)
      wrmsplot(iic,j,k)=wrmsplot(iic,j,k)+
     %dt*(wstagg**2.)
      prmsplot(iic,j,k)=prmsplot(iic,j,k)+
     %dt*(pstagg**2.)
      vmod=sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
      vmodplot(iic,j,k)=vmodplot(iic,j,k)+
     %vmod*dt
      tvplot(iic,j,k)=tvplot(iic,j,k)+
     %tvstagg*dt
  100 continue
  101 continue
!
!     Time averaged values and root mean squares
!     at some cross sections
!
      do 103 kkc=1,nsez2max
      ks=ksez2(kkc)
      kp=ks+1
      do 102 i=ix1,ix2
      im=i-1
      do 102 j=jy1,jy2
      jm=jmv(j)
      tempo5(i,j,kkc)=tempo5(i,j,kkc)+dt
      vstagg=(vo(i,j,ks)+vo(i,jm,ks)+
     %vo(i,j,kp)+vo(i,jm,kp))*0.25
      vplot2(i,j,kkc)=vplot2(i,j,kkc)+
     %vstagg*dt
      ustagg=(uo(i,j,ks)+uo(im,j,ks)+
     %uo(i,j,kp)+uo(im,j,kp))*0.25
      uplot2(i,j,kkc)=uplot2(i,j,kkc)+
     %ustagg*dt
      wstagg=wo(i,j,ks)
      wplot2(i,j,kkc)=wplot2(i,j,kkc)+
     %wstagg*dt
      pstagg=(p(i,j,ks)+p(i,j,kp))*0.5
      tvstagg=(tv(i,j,ks)+tv(i,j,kp))*0.5
      pplot2(i,j,kkc)=pplot2(i,j,kkc)+
     %pstagg*dt
      vrmsplot2(i,j,kkc)=vrmsplot2(i,j,kkc)+
     %dt*(vstagg**2.)
      urmsplot2(i,j,kkc)=urmsplot2(i,j,kkc)+
     %dt*(ustagg**2.)
      wrmsplot2(i,j,kkc)=wrmsplot2(i,j,kkc)+
     %dt*(wstagg**2.)
      prmsplot2(i,j,kkc)=prmsplot2(i,j,kkc)+
     %dt*(pstagg**2.)
      vmod=sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
      vmodplot2(i,j,kkc)=vmodplot2(i,j,kkc)+
     %vmod*dt
      tvplot2(i,j,kkc)=tvplot2(i,j,kkc)+
     %tvstagg*dt
  102 continue
  103 continue
!
!     Time averaged values and root mean squares
!     at some meridian sections
!     (averaged along the azimuthal direction)
!
      jjc=1
      do 111 k=kz1,kz2
      km=k-1
      do 110 i=ix1,ix2
      im=i-1
      tempo14(i,jjc,k)=tempo14(i,jjc,k)+dt
      do 110 j=jy1,jy2
      jp=jpv(j)
      vstagg=vo(i,j,k)
      vplot14(i,jjc,k)=vplot14(i,jjc,k)+
     %vstagg*dt
      ustagg=(uo(i,j,k)+uo(im,j,k)+
     %uo(i,jp,k)+uo(im,jp,k))*0.25
      uplot14(i,jjc,k)=uplot14(i,jjc,k)+
     %ustagg*dt
      wstagg=(wo(i,j,k)+wo(i,jp,k)+
     %wo(i,j,km)+wo(i,jp,km))*0.25
      wplot14(i,jjc,k)=wplot14(i,jjc,k)+
     %wstagg*dt
      pstagg=(p(i,j,k)+p(i,jp,k))*0.5
      tvstagg=(tv(i,j,k)+tv(i,jp,k))*0.5
      pplot14(i,jjc,k)=pplot14(i,jjc,k)+
     %pstagg*dt
      vrmsplot14(i,jjc,k)=vrmsplot14(i,jjc,k)+
     %dt*(vstagg**2.)
      urmsplot14(i,jjc,k)=urmsplot14(i,jjc,k)+
     %dt*(ustagg**2.)
      wrmsplot14(i,jjc,k)=wrmsplot14(i,jjc,k)+
     %dt*(wstagg**2.)
      prmsplot14(i,jjc,k)=prmsplot14(i,jjc,k)+
     %dt*(pstagg**2.)
      vmod=sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
      vmodplot14(i,jjc,k)=vmodplot14(i,jjc,k)+
     %vmod*dt
      tvplot14(i,jjc,k)=tvplot14(i,jjc,k)+
     %tvstagg*dt
  110 continue
  111 continue
!
!!!	 if((time-timep).gt.((periodo)*(iplant5+1))) then
c      if((abs(romega)*time_plot2*180./pig).gt.(rplot1*(iplant5+1)))
	 call reynoldsstaggavg(1,nx,ny,nz,uo,vo,wo,dt)
	 call avgvort2staggavg(1,nx,ny,nz,dt,xc,xu,zc,uo,vo,wo)
!!!	 endif
!
      if(tempo3.ge.periodo2) then
c
	 do 105 iic=1,nsez1
	 do 104 k=kz1,kz2
	 do 104 j=jy1,jy2
	    vplot(iic,j,k)=
     %vplot(iic,j,k)/tempo4(iic,j,k)
	    uplot(iic,j,k)=
     %uplot(iic,j,k)/tempo4(iic,j,k)
	    wplot(iic,j,k)=
     %wplot(iic,j,k)/tempo4(iic,j,k)
	    pplot(iic,j,k)=
     %pplot(iic,j,k)/tempo4(iic,j,k)
	    vrmsplot(iic,j,k)=
     %vrmsplot(iic,j,k)/tempo4(iic,j,k)
     %-vplot(iic,j,k)**2.
	    urmsplot(iic,j,k)=
     %urmsplot(iic,j,k)/tempo4(iic,j,k)
     %-uplot(iic,j,k)**2.
	    wrmsplot(iic,j,k)=
     %wrmsplot(iic,j,k)/tempo4(iic,j,k)
     %-wplot(iic,j,k)**2.
	    prmsplot(iic,j,k)=
     %prmsplot(iic,j,k)/tempo4(iic,j,k)
     %-pplot(iic,j,k)**2.
	    vmodplot(iic,j,k)=
     %vmodplot(iic,j,k)/tempo4(iic,j,k)
	    tvplot(iic,j,k)=
     %tvplot(iic,j,k)/tempo4(iic,j,k)
	    tvplot(iic,j,k)=
     %tvplot(iic,j,k)/ru1
  104	 continue
  105	 continue
c
	 do 107 kkc=1,nsez2max
	 do 106 i=ix1,ix2
	 do 106 j=jy1,jy2
	    vplot2(i,j,kkc)=
     %vplot2(i,j,kkc)/tempo5(i,j,kkc)
	    uplot2(i,j,kkc)=
     %uplot2(i,j,kkc)/tempo5(i,j,kkc)
	    wplot2(i,j,kkc)=
     %wplot2(i,j,kkc)/tempo5(i,j,kkc)
	    pplot2(i,j,kkc)=
     %pplot2(i,j,kkc)/tempo5(i,j,kkc)
	    vrmsplot2(i,j,kkc)=
     %vrmsplot2(i,j,kkc)/tempo5(i,j,kkc)
     %-vplot2(i,j,kkc)**2.
	    urmsplot2(i,j,kkc)=
     %urmsplot2(i,j,kkc)/tempo5(i,j,kkc)
     %-uplot2(i,j,kkc)**2.
	    wrmsplot2(i,j,kkc)=
     %wrmsplot2(i,j,kkc)/tempo5(i,j,kkc)
     %-wplot2(i,j,kkc)**2.
	    prmsplot2(i,j,kkc)=
     %prmsplot2(i,j,kkc)/tempo5(i,j,kkc)
     %-pplot2(i,j,kkc)**2.
	    vmodplot2(i,j,kkc)=
     %vmodplot2(i,j,kkc)/tempo5(i,j,kkc)
	    tvplot2(i,j,kkc)=
     %tvplot2(i,j,kkc)/tempo5(i,j,kkc)
	    tvplot2(i,j,kkc)=
     %tvplot2(i,j,kkc)/ru1
  106	 continue
  107	 continue
c
	 jjc=1
	 do 115 k=kz1,kz2
	 do 114 i=ix1,ix2
	    vplot14(i,jjc,k)=
     %vplot14(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))
	    uplot14(i,jjc,k)=
     %uplot14(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))
	    wplot14(i,jjc,k)=
     %wplot14(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))
	    pplot14(i,jjc,k)=
     %pplot14(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))
	    vrmsplot14(i,jjc,k)=
     %vrmsplot14(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))
     %-vplot14(i,jjc,k)**2.
	    urmsplot14(i,jjc,k)=
     %urmsplot14(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))
     %-uplot14(i,jjc,k)**2.
	    wrmsplot14(i,jjc,k)=
     %wrmsplot14(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))
     %-wplot14(i,jjc,k)**2.
	    prmsplot14(i,jjc,k)=
     %prmsplot14(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))
     %-pplot14(i,jjc,k)**2.
	    vmodplot14(i,jjc,k)=
     %vmodplot14(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))
	    tvplot14(i,jjc,k)=
     %tvplot14(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))
	    tvplot14(i,jjc,k)=
     %tvplot14(i,jjc,k)/ru1
  114	 continue
  115	 continue
!
	 iplant4=iplant4+1
	 call kplotstaggavg(iplant4c,iplant4d,iplant4p,
     %nx,ny,nz,nzg,nbd,ntime,xu,xc,yv,yc,zwg,zcg,
     %flaguo,flagvo,flagwo)
!
ccc 13	    format(2x,a12,1x,e14.7,1x,e14.7,1x,e14.7)
ccc 14	    format(6x,a6,3x,e14.7,1x,e14.7,1x,e14.7)
ccc	    if((ibm.gt.1).and.(myrank.eq.0)) then
ccc	       do ibd=mbd,nbd
ccc		  open(83,file='ib_n.'//index(ibd)//'_'
ccc	%//char(48+iplant4c)//char(48+iplant4d)//char(48+iplant4p)//
ccc	%'_kplot1.stl')
ccc		  write(83,*)'solid  OBJECT'
ccc		  do j=lb(ibd)+1,lb(ibd)+mb(ibd)
ccc		     write(83,13)'facet normal',unvect(1,j),unvect(2,j),
ccc	%unvect(3,j)
ccc		     write(83,*)'outer loop'
ccc		     write(83,14)
ccc	%'vertex',vertex(1,1,j),vertex(2,1,j),vertex(3,1,j)
ccc		     write(83,14)
ccc	%'vertex',vertex(1,2,j),vertex(2,2,j),vertex(3,2,j)
ccc		     write(83,14)
ccc	%'vertex',vertex(1,3,j),vertex(2,3,j),vertex(3,3,j)
ccc		     write(83,*)'endloop'
ccc		     write(83,*)'endfacet'
ccc		  enddo
ccc		  write(83,*)'endsolid	OBJECT'
ccc		  close(83)
ccc	       enddo
ccc	    endif
!
!!!	    if((time-timep).gt.((periodo)*(iplant5+1))) then
c	  if((abs(romega)*time_plot2*180./pig).gt.
c     %(rplot1*(iplant5+1))) then
	    iplant5=iplant5+1
	    call reynoldsstaggavg(2,nx,ny,nz,uo,vo,wo,dt)
	    call kplot2avg(iplant5c,iplant5d,iplant5p,
     %nx,ny,nz,nzg,nbd,ntime,flaguo,flagvo,flagwo)
	    call avgvort2staggavg(2,nx,ny,nz,dt,xc,xu,zc,uo,vo,wo)
	    call plotavgvortavg(iplant5c,iplant5d,iplant5p,
     %nx,ny,nz,nzg,nbd,ntime,flaguo,flagvo,flagwo)
c	     call vort2
c	     call plotvort(iplant5c,iplant5d,iplant5p)
!
ccc	    if((ibm.gt.1).and.(myrank.eq.0)) then
ccc	       do ibd=mbd,nbd
ccc		  open(83,file='ib_n.'//index(ibd)//'_'
ccc	%//char(48+iplant4c)//char(48+iplant4d)//char(48+iplant4p)//
ccc	%'_kplot2.stl')
ccc		  write(83,*)'solid  OBJECT'
ccc		  do j=lb(ibd)+1,lb(ibd)+mb(ibd)
ccc		     write(83,13)'facet normal',unvect(1,j),unvect(2,j),
ccc	%unvect(3,j)
ccc		     write(83,*)'outer loop'
ccc		     write(83,14)
ccc	%'vertex',vertex(1,1,j),vertex(2,1,j),vertex(3,1,j)
ccc		     write(83,14)
ccc	%'vertex',vertex(1,2,j),vertex(2,2,j),vertex(3,2,j)
ccc		     write(83,14)
ccc	%'vertex',vertex(1,3,j),vertex(2,3,j),vertex(3,3,j)
ccc		     write(83,*)'endloop'
ccc		     write(83,*)'endfacet'
ccc		  enddo
ccc		  write(83,*)'endsolid	OBJECT'
ccc		  close(83)
ccc	       enddo
ccc	    endif
!
!!!	    endif
!
	 do iic=1,nsez1
	    do k=kz1,kz2
	       do j=jy1,jy2
		  vrmsplot(iic,j,k)=0.
		  urmsplot(iic,j,k)=0.
		  wrmsplot(iic,j,k)=0.
		  prmsplot(iic,j,k)=0.
		  vplot(iic,j,k)=0.
		  uplot(iic,j,k)=0.
		  wplot(iic,j,k)=0.
		  pplot(iic,j,k)=0.
		  tempo4(iic,j,k)=0.
		  uvplot(iic,j,k)=0.
		  vwplot(iic,j,k)=0.
		  uwplot(iic,j,k)=0.
		  voraz1med(iic,j,k)=0.
		  vorr1med(iic,j,k)=0.
		  vorz1med(iic,j,k)=0.
		  vmodplot(iic,j,k)=0.
		  tvplot(iic,j,k)=0.
	       enddo
	    enddo
	 enddo
c
	 do kkc=1,nsez2max
	    do i=ix1,ix2
	       do j=jy1,jy2
		  vrmsplot2(i,j,kkc)=0.
		  urmsplot2(i,j,kkc)=0.
		  wrmsplot2(i,j,kkc)=0.
		  prmsplot2(i,j,kkc)=0.
		  vplot2(i,j,kkc)=0.
		  uplot2(i,j,kkc)=0.
		  wplot2(i,j,kkc)=0.
		  pplot2(i,j,kkc)=0.
		  tempo5(i,j,kkc)=0.
		  uvplot2(i,j,kkc)=0.
		  vwplot2(i,j,kkc)=0.
		  uwplot2(i,j,kkc)=0.
		  voraz2med(i,j,kkc)=0.
		  vorr2med(i,j,kkc)=0.
		  vorz2med(i,j,kkc)=0.
		  vmodplot2(i,j,kkc)=0.
		  tvplot2(i,j,kkc)=0.
	       enddo
	    enddo
	 enddo
c
	 jjc=1
	    do k=kz1,kz2
	       do i=ix1,ix2
		  vrmsplot14(i,jjc,k)=0.
		  urmsplot14(i,jjc,k)=0.
		  wrmsplot14(i,jjc,k)=0.
		  prmsplot14(i,jjc,k)=0.
		  vplot14(i,jjc,k)=0.
		  uplot14(i,jjc,k)=0.
		  wplot14(i,jjc,k)=0.
		  pplot14(i,jjc,k)=0.
		  tempo14(i,jjc,k)=0.
		  uvplot14(i,jjc,k)=0.
		  vwplot14(i,jjc,k)=0.
		  uwplot14(i,jjc,k)=0.
		  voraz14med(i,jjc,k)=0.
		  vorr14med(i,jjc,k)=0.
		  vorz14med(i,jjc,k)=0.
		  vmodplot14(i,jjc,k)=0.
		  tvplot14(i,jjc,k)=0.
	       enddo
	    enddo
	 tempo3=0.
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine reynoldsstaggavg(ind,nx,ny,nz,uo,vo,wo,dt)
c
      use sections
      use mediecalc
      use mediecalc1
      use mediecalc2
      use mediecalc3
      include'common.h'
      include'averages.h'
c
c Global variables
      integer ind,nx,ny,nz
      real dt,uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
c
c Local variables
      integer iic,i,is,ip,im,jjc,j,jp,jm,kkc,k,ks,kp,km
      real ustagg,vstagg,wstagg
c
      if(ind.eq.2) goto 200
!
!     Circumferential sections
!
      do 101 iic=1,nsez1
      is=isez1(iic)
      ip=is+1
      do 100 k=kz1,kz2
      km=k-1
      do 100 j=jy1,jy2
      jm=jmv(j)
      ustagg=
     %uo(is,j,k)
      vstagg=
     %(vo(is,j,k)+vo(ip,j,k)+vo(is,jm,k)+vo(ip,jm,k))*0.25
      wstagg=
     %(wo(is,j,k)+wo(ip,j,k)+wo(is,j,km)+wo(ip,j,km))*0.25
      uvplot(iic,j,k)=uvplot(iic,j,k)+
     %ustagg*vstagg*dt
      uwplot(iic,j,k)=uwplot(iic,j,k)+
     %ustagg*wstagg*dt
      vwplot(iic,j,k)=vwplot(iic,j,k)+
     %vstagg*wstagg*dt
  100 continue
  101 continue
!
!     Cross sections
!
      do 103 kkc=1,nsez2max
      ks=ksez2(kkc)
      kp=ks+1
      do 102 i=ix1,ix2
      im=i-1
      do 102 j=jy1,jy2
      jm=jmv(j)
      ustagg=
     %(uo(i,j,ks)+uo(im,j,ks)+uo(i,j,kp)+uo(im,j,kp))*0.25
      vstagg=
     %(vo(i,j,ks)+vo(i,jm,ks)+vo(i,j,kp)+vo(i,jm,kp))*0.25
      wstagg=
     %wo(i,j,ks)
      uvplot2(i,j,kkc)=uvplot2(i,j,kkc)+
     %ustagg*vstagg*dt
      uwplot2(i,j,kkc)=uwplot2(i,j,kkc)+
     %ustagg*wstagg*dt
      vwplot2(i,j,kkc)=vwplot2(i,j,kkc)+
     %vstagg*wstagg*dt
  102 continue
  103 continue
!
!     Meridian sections (azimuthal averages)
!
      jjc=1
      do 111 k=kz1,kz2
      km=k-1
      do 110 i=ix1,ix2
      im=i-1
      do 110 j=jy1,jy2
      jp=jpv(j)
      ustagg=
     %(uo(i,j,k)+uo(im,j,k)+uo(i,jp,k)+uo(im,jp,k))*0.25
      vstagg=
     %vo(i,j,k)
      wstagg=
     %(wo(i,j,k)+wo(i,jp,k)+wo(i,j,km)+wo(i,jp,km))*0.25
      uvplot14(i,jjc,k)=uvplot14(i,jjc,k)+
     %ustagg*vstagg*dt
      uwplot14(i,jjc,k)=uwplot14(i,jjc,k)+
     %ustagg*wstagg*dt
      vwplot14(i,jjc,k)=vwplot14(i,jjc,k)+
     %vstagg*wstagg*dt
  110 continue
  111 continue
      return
  200 continue
!
!     Circumferential sections
!
      do 105 iic=1,nsez1
      do 104 k=kz1,kz2
      do 104 j=jy1,jy2
	 uvplot(iic,j,k)=
     %uvplot(iic,j,k)/tempo4(iic,j,k)-
     %uplot(iic,j,k)*vplot(iic,j,k)
	 uwplot(iic,j,k)=
     %uwplot(iic,j,k)/tempo4(iic,j,k)-
     %uplot(iic,j,k)*wplot(iic,j,k)
	 vwplot(iic,j,k)=
     %vwplot(iic,j,k)/tempo4(iic,j,k)-
     %vplot(iic,j,k)*wplot(iic,j,k)
  104 continue
  105 continue
!
!     Cross sections
!
      do 107 kkc=1,nsez2max
      do 106 i=ix1,ix2
      do 106 j=jy1,jy2
	 uvplot2(i,j,kkc)=
     %uvplot2(i,j,kkc)/tempo5(i,j,kkc)-
     %uplot2(i,j,kkc)*vplot2(i,j,kkc)
	 uwplot2(i,j,kkc)=
     %uwplot2(i,j,kkc)/tempo5(i,j,kkc)-
     %uplot2(i,j,kkc)*wplot2(i,j,kkc)
	 vwplot2(i,j,kkc)=
     %vwplot2(i,j,kkc)/tempo5(i,j,kkc)-
     %vplot2(i,j,kkc)*wplot2(i,j,kkc)
  106 continue
  107 continue
!
!     Meridian sections (azimuthal averages)
!
      jjc=1
      do 115 k=kz1,kz2
      do 114 i=ix1,ix2
	 uvplot14(i,jjc,k)=
     %uvplot14(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))-
     %uplot14(i,jjc,k)*vplot14(i,jjc,k)
	 uwplot14(i,jjc,k)=
     %uwplot14(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))-
     %uplot14(i,jjc,k)*wplot14(i,jjc,k)
	 vwplot14(i,jjc,k)=
     %vwplot14(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))-
     %vplot14(i,jjc,k)*wplot14(i,jjc,k)
  114 continue
  115 continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine avgvort2staggavg
     %(ind,nx,ny,nz,dt,xc,xu,zc,uo,vo,wo)
c
      use sections
      use mediecalc
      use vortavg
      include'common.h'
      include'averages.h'
c
c Global variables
      integer ind,nx,ny,nz
      real dt,xc(nx),xu(nx),zc(nz),
     %uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
c
c Local variables
      integer iic,ic,ip,im,jjc,jc,jp,jm,jsy,jpsy,kkc,kc,kp,km,
     %i,j,k
      real dtheta,dr,dz,dq1x3,dq3x1,dq3x2,dq2x3,dq1x2,dq2x1
c
      if(ind.eq.2) goto 200
c
c    ! Time-averaged azimuthal vorticity - Circumferential sections !
      do 101 kc=kz1,kz2
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	do 101 iic=1,nsez1
	  ic=isez1(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 101 jc=jy1,jy2
	    dq1x3=(uo(ic,jc,kp)
     %-uo(ic,jc,km))/dz
	    dq3x1=(0.5*(wo(ip,jc,kc)+wo(ip,jc,km))
     %-0.5*(wo(ic,jc,kc)+wo(ic,jc,km)))/dr
	    voraz1med(iic,jc,kc)=
     %voraz1med(iic,jc,kc)+(dq1x3-dq3x1)*dt
  101 continue
c
c    ! Time-averaged azimuthal vorticity - Cross sections !
      do 111 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
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
	  voraz2med(ic,jc,kkc)=
     %voraz2med(ic,jc,kkc)+(dq1x3-dq3x1)*dt
  112 continue
	do 113 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 113 jc=jy1,jy2
	    dq1x3=(0.5*(uo(ic,jc,kp)+uo(im,jc,kp))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dz
	    dq3x1=(wo(ip,jc,kc)
     %-wo(im,jc,kc))/dr
	    voraz2med(ic,jc,kkc)=
     %voraz2med(ic,jc,kkc)+(dq1x3-dq3x1)*dt
  113 continue
  111 continue
c
c    ! Time-averaged azimuthal vorticity - Meridian sections !
c    ! Azimuthal averages !
      jjc=1
      do 121 kc=kz1,kz2
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	ic=ix1
	ip=ic+1
	im=ic-1
	dr=xc(ip)-xc(im)
	do 122 jc=jy1,jy2
	  jp=jpv(jc)
	  jsy=jsym(jc)
	  jpsy=jsym(jp)
	  dq1x3=(0.25*(uo(ic,jc,kp)+uo(im,jc,kp)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,km)+uo(im,jc,km)
     %+uo(ic,jp,km)+uo(im,jp,km)))/dz
	  dq3x1=(0.25*(wo(ip,jc,kc)+wo(ip,jp,kc)
     %+wo(ip,jc,km)+wo(ip,jp,km))
     %-0.25*(wo(ic,jsy,kc)+wo(ic,jpsy,kc)
     %+wo(ic,jsy,km)+wo(ic,jpsy,km)))/dr
	  voraz14med(ic,jjc,kc)=
     %voraz14med(ic,jjc,kc)+(dq1x3-dq3x1)*dt
  122 continue
	do 123 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 123 jc=jy1,jy2
	    jp=jpv(jc)
	    dq1x3=(0.25*(uo(ic,jc,kp)+uo(im,jc,kp)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,km)+uo(im,jc,km)
     %+uo(ic,jp,km)+uo(im,jp,km)))/dz
	    dq3x1=(0.25*(wo(ip,jc,kc)+wo(ip,jp,kc)
     %+wo(ip,jc,km)+wo(ip,jp,km))
     %-0.25*(wo(im,jc,kc)+wo(im,jp,kc)
     %+wo(im,jc,km)+wo(im,jp,km)))/dr
	    voraz14med(ic,jjc,kc)=
     %voraz14med(ic,jjc,kc)+(dq1x3-dq3x1)*dt
  123 continue
  121 continue
c
c    ! Time-averaged radial vorticity - Circumferential sections !
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
	    dq3x2=(0.25*(wo(ic,jp,kc)+wo(ip,jp,kc)
     %+wo(ic,jp,km)+wo(ip,jp,km))
     %-0.25*(wo(ic,jm,kc)+wo(ip,jm,kc)
     %+wo(ic,jm,km)+wo(ip,jm,km)))/dtheta    !
	    dq2x3=(0.25*(vo(ic,jc,kp)+vo(ip,jc,kp)
     %+vo(ic,jm,kp)+vo(ip,jm,kp))
     %-0.25*(vo(ic,jc,km)+vo(ip,jc,km)
     %+vo(ic,jm,km)+vo(ip,jm,km)))/dz
	    vorr1med(iic,jc,kc)=
     %vorr1med(iic,jc,kc)+(dq3x2/xu(ic)-dq2x3)*dt
  131 continue
c
c    ! Time-averaged radial vorticity - Cross sections !
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
     %-wo(ic,jm,kc))/dtheta    !
	    dq2x3=(0.5*(vo(ic,jc,kp)+vo(ic,jm,kp))
     %-0.5*(vo(ic,jc,kc)+vo(ic,jm,kc)))/dz
	    vorr2med(ic,jc,kkc)=
     %vorr2med(ic,jc,kkc)+(dq3x2/xc(ic)-dq2x3)*dt
  141 continue
c
c    ! Time-averaged radial vorticity - Meridian sections !
c    ! Azimuthal averages !
      jjc=1
      dtheta=dely
      do 151 kc=kz1,kz2
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	do 151 ic=ix1,ix2
	  do 151 jc=jy1,jy2
	    jp=jpv(jc)
	    dq3x2=(0.5*(wo(ic,jp,kc)+wo(ic,jp,km))
     %-0.5*(wo(ic,jc,kc)+wo(ic,jc,km)))/dtheta	  !
	    dq2x3=(vo(ic,jc,kp)
     %-vo(ic,jc,km))/dz
	    vorr14med(ic,jjc,kc)=
     %vorr14med(ic,jjc,kc)+(dq3x2/xc(ic)-dq2x3)*dt
  151 continue
c
c    ! Time-averaged axial vorticity - Circumferential sections !
      do 161 kc=kz1,kz2
	do 161 iic=1,nsez1
	  ic=isez1(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 161 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.5*(vo(ip,jc,kc)+vo(ip,jm,kc))*xc(ip)
     %-0.5*(vo(ic,jc,kc)+vo(ic,jm,kc))*xc(ic))/dr
	    dq1x2=(uo(ic,jp,kc)
     %-uo(ic,jm,kc))/dtheta    !
	    vorz1med(iic,jc,kc)=
     %vorz1med(iic,jc,kc)+((dq2x1-dq1x2)/xu(ic))*dt
  161 continue
c
c    ! Time-averaged axial vorticity - Cross sections !
      do 171 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
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
     %+vo(ic,jm,kp)+vo(ip,jm,kp))*xu(ic))/dr
	  dq1x2=(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jm,kc)+uo(im,jm,kc)
     %+uo(ic,jm,kp)+uo(im,jm,kp)))/dtheta   !
	  vorz2med(ic,jc,kkc)=
     %vorz2med(ic,jc,kkc)+((dq2x1-dq1x2)/xc(ic))*dt
  172 continue
c
	do 173 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 173 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.25*(vo(ip,jc,kc)+vo(ip,jm,kc)
     %+vo(ip,jc,kp)+vo(ip,jm,kp))*xc(ip)
     %-0.25*(vo(im,jc,kc)+vo(im,jm,kc)
     %+vo(im,jc,kp)+vo(im,jm,kp))*xc(im))/dr
	    dq1x2=
     %(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jm,kc)+uo(im,jm,kc)
     %+uo(ic,jm,kp)+uo(im,jm,kp)))/dtheta    !
	    vorz2med(ic,jc,kkc)=
     %vorz2med(ic,jc,kkc)+((dq2x1-dq1x2)/xc(ic))*dt
  173 continue
  171 continue
c
c    ! Time-averaged axial vorticity - Meridian sections !
c    ! Azimuthal averages !
      jjc=1
      dtheta=dely
      do 181 kc=kz1,kz2
	ic=ix1
	im=ic-1
	ip=ic+1
	dr=xu(ic)-xu(im)
	do 182 jc=jy1,jy2
	  jp=jpv(jc)
	  dq2x1=(0.5*(vo(ic,jc,kc)+vo(ip,jc,kc))*xu(ic))/dr
	  dq1x2=(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dtheta	  !
	  vorz14med(ic,jjc,kc)=
     %vorz14med(ic,jjc,kc)+((dq2x1-dq1x2)/xc(ic))*dt
  182	continue
c
	do 183 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 183 jc=jy1,jy2
	    jp=jpv(jc)
	    dq2x1=(vo(ip,jc,kc)*xc(ip)
     %-vo(im,jc,kc)*xc(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dtheta	 !
	    vorz14med(ic,jjc,kc)=
     %vorz14med(ic,jjc,kc)+((dq2x1-dq1x2)/xc(ic))*dt
  183 continue
  181 continue
      return
c
  200 continue
c
!
!     Circumferential sections
!
      do 205 iic=1,nsez1
      do 204 k=kz1,kz2
      do 204 j=jy1,jy2
	 voraz1med(iic,j,k)=
     %voraz1med(iic,j,k)/tempo4(iic,j,k)
	 vorr1med(iic,j,k)=
     %vorr1med(iic,j,k)/tempo4(iic,j,k)
	 vorz1med(iic,j,k)=
     %vorz1med(iic,j,k)/tempo4(iic,j,k)
  204 continue
  205 continue
c
!
!     Cross sections
!
      do 207 kkc=1,nsez2max
      do 206 i=ix1,ix2
      do 206 j=jy1,jy2
	 voraz2med(i,j,kkc)=
     %voraz2med(i,j,kkc)/tempo5(i,j,kkc)
	 vorr2med(i,j,kkc)=
     %vorr2med(i,j,kkc)/tempo5(i,j,kkc)
	 vorz2med(i,j,kkc)=
     %vorz2med(i,j,kkc)/tempo5(i,j,kkc)
  206 continue
  207 continue
c
!
!     Meridian sections
!
      jjc=1
      do 215 k=kz1,kz2
      do 214 i=ix1,ix2
	 voraz14med(i,jjc,k)=
     %voraz14med(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))
	 vorr14med(i,jjc,k)=
     %vorr14med(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))
	 vorz14med(i,jjc,k)=
     %vorz14med(i,jjc,k)/(tempo14(i,jjc,k)*real(ny-2))
  214 continue
  215 continue
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine kplotstaggavg(iplant4c,iplant4d,iplant4p,
     %nx,ny,nz,nzg,nbd,ntime,xu,xc,yv,yc,zwg,zcg,flaguo,flagvo,
     %flagwo)
c
      use sections
      use mediecalc1
      use mediecalc2
      use mediecalc3
      include'common.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer nx,ny,nz,nzg,nbd,iplant4p,iplant4d,iplant4c,ntime,
     %flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),flagwo(nx,ny,nz,nbd)
      real xu(nx),xc(nx),yv(ny),yc(ny),zwg(nzg),zcg(nzg)
c
c Local variables
      integer iic,iicd,iicu,is,jjc,jjcd,jjcu,js,kkc,kkcd,kkcu,ks,
     %i,j,k,l,i1,jjp,ir,jr,kr,lr,kk,kkd,kku
!     integer kk1,kk2,krg,flaguor(ny,nzg,nbd),flagvor2(nx,nzg,nbd)
      integer STATUS(MPI_STATUS_SIZE)
      real recvar(ny,nz),recvar2(nx,nz)
c
c Parameters
      integer lamb2,nxp,nyp,nzp
      real deltasect
      parameter (lamb2=1,nyp=1,nxp=1,nzp=1)
      parameter (deltasect=2.*pi)
c
   8  format(10f25.8)
   9  format(10f25.8)
 101  format(3i8)
 102  format(i8,3i2)
 103  format(i2)
c
      i1=1
      if(iplant4.eq.1) then
	 if(myrank.eq.0) then
	    do iic=1,nsez1
	       iicd=iic/10
	       iicu=mod(iic,10)
	       open
     %(81,file='grid1_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'_kplot.xyz')
	       open
     %(82,file='grid1_new_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'_kplot.xyz')
	       is=isez1(iic)
	       write(81,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	       write(82,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	       do 1 k=2,nzg-1,nzp
	       do 10 l=1,lamb2
	       do 10 j=jy1,jy2,nyp
	       write(81,9)xu(is)*cos(yc(jy1))
   10	       write(82,9)xu(is)*cos(yc(j)+(l-1)*2.*pi/lamb2)
	       write(81,9)xu(is)*cos(yc(jy1))
    1	       write(82,9)xu(is)*cos(yc(jy1)+lamb2*deltasect)
	       do 2 k=2,nzg-1,nzp
	       do 20 l=1,lamb2
	       do 20 j=jy1,jy2,nyp
	       write(81,9)
     %xu(is)*sin(yc(jy1))+xu(is)*(yc(j)+(l-1)*2.*pi/lamb2)
   20	       write(82,9)xu(is)*sin(yc(j)+(l-1)*2.*pi/lamb2)
	       write(81,9)
     %xu(is)*sin(yc(jy1))+xu(is)*deltasect*lamb2
   2	       write(82,9)xu(is)*sin(yc(jy1)+lamb2*deltasect)
	       do 3 k=2,nzg-1,nzp
	       do 30 l=1,lamb2
	       do 30 j=jy1,jy2,nyp
	       write(81,9)zcg(k)
   30	       write(82,9)zcg(k)
	       write(81,9)zcg(k)
   3	       write(82,9)zcg(k)
	       close(81)
	       close(82)
	    enddo
	    do kkc=1,nsez2
	       kkcd=kkc/10
	       kkcu=mod(kkc,10)
	       open(81,file='grid2_2D_slice'//char(48+kkcd)
     %//char(48+kkcu)//'_kplot.xyz')
	       ks=ksez2g(kkc)
	       write(81,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	       do 11 i=ix1,ix2,nxp
	       do 110 l=1,lamb2
	       do 110 j=jy1,jy2,nyp
 110	       write(81,9)xc(i)*cos(yc(j)+(l-1)*2.*pi/lamb2)
  11	       write(81,9)xc(i)*cos(yc(jy1)+lamb2*deltasect)
	       do 12 i=ix1,ix2,nxp
	       do 120 l=1,lamb2
	       do 120 j=jy1,jy2,nyp
 120	       write(81,9)xc(i)*sin(yc(j)+(l-1)*2.*pi/lamb2)
  12	       write(81,9)xc(i)*sin(yc(jy1)+lamb2*deltasect)
	       do 13 i=ix1,ix2,nxp
	       do 130 l=1,lamb2
	       do 130 j=jy1,jy2,nyp
 130	       write(81,9)zwg(ks)
  13	       write(81,9)zwg(ks)
	       close(81)
	    enddo
	    jjc=1
	    jjcd=jjc/10
	    jjcu=mod(jjc,10)
	    open(81,file='grid3_2D_slice'//char(48+jjcd)
     %//char(48+jjcu)//'_kplot.xyz')
	    js=jsez14(jjc)
	    write(81,101)1,(nx-2)/nxp,(nzg-2)/nzp
	    do 210 k=2,nzg-1,nzp
	    do 210 i=ix1,ix2,nxp
  210	    write(81,9)xc(i)*cos(yv(js))
	    do 220 k=2,nzg-1,nzp
	    do 220 i=ix1,ix2,nxp
  220	    write(81,9)xc(i)*sin(yv(js))
	    do 230 k=2,nzg-1,nzp
	    do 230 i=ix1,ix2,nxp
  230	    write(81,9)zcg(k)
	    close(81)
	 endif
      endif
c
      iplant4p=iplant4
      iplant4c=iplant4p/100
      iplant4p=mod(iplant4p,100)
      iplant4d=iplant4p/10
      iplant4p=mod(iplant4p,10)
c
      do iic=1,nsez1
	 is=isez1(iic)
	 if(myrank.eq.0) then
	    iicd=iic/10
	    iicu=mod(iic,10)
	    open(82,file='results1_2D_'//char(48+iplant4c)
     %//char(48+iplant4d)//char(48+iplant4p)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'_kplot1.q')
	    write(82,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	    write(82,102)ntime,i1,i1,i1
!
! Turbulent kinetic energy
!
	    do 40 k=kz1,kz2,nzp
	    do 41 l=1,lamb2
	    do 41 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)
     %0.5*(vrmsplot(iic,j,k)+urmsplot(iic,j,k)+
     %wrmsplot(iic,j,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
 41	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)
     %0.5*(vrmsplot(iic,jy1,k)+urmsplot(iic,jy1,k)+
     %wrmsplot(iic,jy1,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
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
		     write(82,9)recvar(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 43		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(82,9)recvar(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 42		  continue
	       enddo
	    endif
c	    write(82,*)	  !!!!!!
c	    write(82,103)
c    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
!
! Time-averaged azimuthal velocity
!
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)vplot(iic,j,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 51	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)vplot(iic,jy1,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
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
		     write(82,9)recvar(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 53		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(82,9)recvar(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 52		  continue
	       enddo
	    endif
!
! Time-averaged radial velocity
!
	    do 60 k=kz1,kz2,nzp
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)uplot(iic,j,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 61	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)uplot(iic,jy1,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
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
		     write(82,9)recvar(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 63		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(82,9)recvar(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 62		  continue
	       enddo
	    endif
!
! Time-averaged axial velocity
!
	    do 70 k=kz1,kz2,nzp
	    do 71 l=1,lamb2
	    do 71 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)wplot(iic,j,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 71	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)wplot(iic,jy1,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
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
		     write(82,9)recvar(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 73		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(82,9)recvar(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 72		  continue
	       enddo
	    endif
!
! Time-averaged static pressure
!
	    do 80 k=kz1,kz2,nzp
	    do 81 l=1,lamb2
	    do 81 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,8)pplot(iic,j,k)    !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
 81	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,8)pplot(iic,jy1,k)   !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
 80	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,6,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 82 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 83 lr=1,lamb2
		  do 83 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(82,8)recvar(jr,kr)	 !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
 83		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(82,8)recvar(jy1,kr)	 !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
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
     %(0.5*(vrmsplot(iic,jy1:jy2,kz1:kz2)+
     %urmsplot(iic,jy1:jy2,kz1:kz2)+
     %wrmsplot(iic,jy1:jy2,kz1:kz2)),1*(ny-2)*(nz-2),MTYPE,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(wplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(pplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
      do kkc=1,nsez2max
	 ks=ksez2(kkc)
	 kk=isezindx2+(kkc-1)
	 kkd=kk/10
	 kku=mod(kk,10)
	 open(82,file='results2_2D_'//char(48+iplant4c)
     %//char(48+iplant4d)//char(48+iplant4p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_kplot1.q')
	 write(82,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	 write(82,102)ntime,i1,i1,i1
!
! Turbulent kinetic energy
!
	 do 140 i=ix1,ix2,nxp
	 do 141 l=1,lamb2
	 do 141 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)
     %0.5*(vrmsplot2(i,j,kkc)+urmsplot2(i,j,kkc)+
     %wrmsplot2(i,j,kkc))
c	 else
c	    write(82,9)1.e+03
c	 endif
 141	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)
     %0.5*(vrmsplot2(i,jy1,kkc)+urmsplot2(i,jy1,kkc)+
     %wrmsplot2(i,jy1,kkc))
c	 else
c	    write(82,9)1.e+03
c	 endif
 140	 continue
c	 write(82,*)   !!!!!!
c	 write(82,103)
c    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
!
! Time-averaged azimuthal velocity
!
	 do 150 i=ix1,ix2,nxp
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)vplot2(i,j,kkc)
c	 else
c	    write(82,9)1.e+03
c	 endif
 151	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)vplot2(i,jy1,kkc)
c	 else
c	    write(82,9)1.e+03
c	 endif
 150	 continue
!
! Time-averaged radial velocity
!
	 do 160 i=ix1,ix2,nxp
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)uplot2(i,j,kkc)
c	 else
c	    write(82,9)1.e+03
c	 endif
 161	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)uplot2(i,jy1,kkc)
c	 else
c	    write(82,9)1.e+03
c	 endif
 160	 continue
!
! Time-averaged axial velocity
!
	 do 170 i=ix1,ix2,nxp
	 do 171 l=1,lamb2
	 do 171 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)wplot2(i,j,kkc)
c	 else
c	    write(82,9)1.e+03
c	 endif
 171	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)wplot2(i,jy1,kkc)
c	 else
c	    write(82,9)1.e+03
c	 endif
 170	 continue
!
! Time-averaged static pressure
!
	 do 180 i=ix1,ix2,nxp
	 do 181 l=1,lamb2
	 do 181 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,8)pplot2(i,j,kkc)	  !*ros
c	 else
c	    write(82,8)1.e+03
c	 endif
 181	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,8)pplot2(i,jy1,kkc)   !*ros
c	 else
c	    write(82,8)1.e+03
c	 endif
 180	 continue
	 close(82)
      enddo
c
      jjc=1
      js=jsez14(jjc)
      if(myrank.eq.0) then
	 jjcd=jjc/10
	 jjcu=mod(jjc,10)
	 open(82,file='results3_2D_'//char(48+iplant4c)
     %//char(48+iplant4d)//char(48+iplant4p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_kplot1.q')
	 write(82,101)1,(nx-2)/nxp,(nzg-2)/nzp
	 write(82,102)ntime,i1,i1,i1
!
! Turbulent kinetic energy (azimuthal average)
!
	 do 240 k=kz1,kz2,nzp
	 do 240 i=ix1,ix2,nxp
c	 if(all(flagvo(i,js,k).le.0) then
	    write(82,9)
     %0.5*(vrmsplot14(i,jjc,k)+urmsplot14(i,jjc,k)+
     %wrmsplot14(i,jjc,k))
c	 else
c	    write(82,9)1.e+03
c	 endif
 240	 continue
	 if(mysize.ne.1) then
	    do jjp=1,mysize-1
c	       kk1=kz1+jjp*(nz-2)
c	       kk2=kz2+jjp*(nz-2)
c	       CALL MPI_RECV
c    %(flagvor2(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,7,
c    %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 241 kr=kz1,kz2,nzp
c	       krg=kr+jjp*(nz-2)
	       do 241 ir=ix1,ix2,nxp
c	       if(all(flagvor2(ir,krg).le.0) then
		  write(82,9)recvar2(ir,kr)
c	       else
c		 write(82,9)1.e+03
c	       endif
 241	       continue
	    enddo
	 endif
!
! Time-averaged azimuthal velocity (azimuthal average)
!
	 do 250 k=kz1,kz2,nzp
	 do 250 i=ix1,ix2,nxp
c	 if(all(flagvo(i,js,k,:)).le.0) then
	    write(82,9)vplot14(i,jjc,k)
c	 else
c	    write(82,9)1.e+03
c	 endif
 250	 continue
	 if(mysize.ne.1) then
	    do jjp=1,mysize-1
	       CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 251 kr=kz1,kz2,nzp
c	       krg=kr+jjp*(nz-2)
	       do 251 ir=ix1,ix2,nxp
c	       if(all(flagvor2(ir,krg,:)).le.0) then
		  write(82,9)recvar2(ir,kr)
c	       else
c		  write(82,9)1.e+03
c	       endif
 251	       continue
	    enddo
	 endif
!
! Time-averaged radial velocity (azimuthal average)
!
	 do 260 k=kz1,kz2,nzp
	 do 260 i=ix1,ix2,nxp
c	 if(all(flagvo(i,js,k,:)).le.0) then
	    write(82,9)uplot14(i,jjc,k)
c	 else
c	    write(82,9)1.e+03
c	 endif
 260	 continue
	 if(mysize.ne.1) then
	    do jjp=1,mysize-1
	       CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 261 kr=kz1,kz2,nzp
c	       krg=kr+jjp*(nz-2)
	       do 261 ir=ix1,ix2,nxp
c	       if(all(flagvor2(ir,krg,:)).le.0) then
		  write(82,9)recvar2(ir,kr)
c	       else
c		  write(82,9)1.e+03
c	       endif
 261	       continue
	    enddo
	 endif
!
! Time-averaged axial velocity (azimuthal average)
!
	 do 270 k=kz1,kz2,nzp
	 do 270 i=ix1,ix2,nxp
c	 if(all(flagvo(i,js,k,:)).le.0) then
	    write(82,9)wplot14(i,jjc,k)
c	 else
c	    write(82,9)1.e+03
c	 endif
 270	 continue
	 if(mysize.ne.1) then
	    do jjp=1,mysize-1
	       CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 271 kr=kz1,kz2,nzp
c	       krg=kr+jjp*(nz-2)
	       do 271 ir=ix1,ix2,nxp
c	       if(all(flagvor2(ir,krg,:)).le.0) then
		  write(82,9)recvar2(ir,kr)
c	       else
c		  write(82,9)1.e+03
c	       endif
 271	       continue
	    enddo
	 endif
!
! Time-averaged static pressure (azimuthal average)
!
	 do 280 k=kz1,kz2,nzp
	 do 280 i=ix1,ix2,nxp
c	 if(all(flagvo(i,js,k,:)).le.0) then
	    write(82,8)pplot14(i,jjc,k)	   !*ros
c	 else
c	    write(82,8)1.e+03
c	 endif
 280	 continue
	 if(mysize.ne.1) then
	    do jjp=1,mysize-1
	       CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 281 kr=kz1,kz2,nzp
c	       krg=kr+jjp*(nz-2)
	       do 281 ir=ix1,ix2,nxp
c	       if(all(flagvor2(ir,krg,:)).le.0) then
		  write(82,9)recvar2(ir,kr)
c	       else
c		  write(82,9)1.e+03
c	       endif
 281	       continue
	    enddo
	 endif
	 close(82)
      else
c	 CALL MPI_SEND
c    %(flagvo(ix1:ix2,js,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,7,
c    %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(0.5*(vrmsplot14(ix1:ix2,jjc,kz1:kz2)+
     %urmsplot14(ix1:ix2,jjc,kz1:kz2)+
     %wrmsplot14(ix1:ix2,jjc,kz1:kz2)),(nx-2)*1*(nz-2),MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(vplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(uplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(wplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(pplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine kplot2avg(iplant5c,iplant5d,iplant5p,
     %nx,ny,nz,nzg,nbd,ntime,flaguo,flagvo,flagwo)
c
      use sections
      use mediecalc1
      use mediecalc2
      use mediecalc3
      include'common.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer iplant5p,iplant5d,iplant5c,nx,ny,nz,nzg,nbd,ntime
      integer flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),
     %flagwo(nx,ny,nz,nbd)
c
c Local variables
      integer iic,iicd,iicu,is,jjc,jjcd,jjcu,js,kkc,ks,kk,kkd,kku,
     %i,j,k,l,ir,jr,kr,lr,jjp,i1
      integer STATUS(MPI_STATUS_SIZE)
!     integer kk1,kk2,krg,flaguor(ny,nzg,nbd),flagvor2(nx,nzg,nbd)
      real recvar1(ny,nz),recvar2(ny,nz)
      real recvar3(nx,nz),recvar4(nx,nz)
c
c Parameters
      integer lamb2,nxp,nyp,nzp
      parameter (lamb2=1)
      parameter (nyp=1,nxp=1,nzp=1)
c
   8  format(10f25.8)
   9  format(10f25.8)
 101  format(3i8)
 102  format(i8,3i2)
 103  format(i2)
c
      iplant5p=iplant5
      iplant5c=iplant5p/100
      iplant5p=mod(iplant5p,100)
      iplant5d=iplant5p/10
      iplant5p=mod(iplant5p,10)
c
      i1=1
      do iic=1,nsez1
	 is=isez1(iic)
	 if(myrank.eq.0) then
	    iicd=iic/10
	    iicu=mod(iic,10)
	    open(83,file='results1_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'_kplot2.q')
	    open(84,file='results1_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'_kplot3.q')
	    write(83,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	    write(84,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	    write(83,102)ntime,i1,i1,i1
	    write(84,102)ntime,i1,i1,i1
!
! Time-averaged velocity magnitude and stagnation pressure
!
	    do 40 k=kz1,kz2,nzp
	    do 41 l=1,lamb2
	    do 41 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)vmodplot(iic,j,k)
	       write(84,8)
     %(pplot(iic,j,k)+0.5*(vmodplot(iic,j,k)**2.))     !*ros
c	    else
c	       write(83,9)1.e+03
c	       write(84,8)1.e+03
c	    endif
 41	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)vmodplot(iic,jy1,k)
	       write(84,8)
     %(pplot(iic,jy1,k)+0.5*(vmodplot(iic,jy1,k)**2.))	   !*ros
c	    else
c	       write(83,9)1.e+03
c	       write(84,8)1.e+03
c	    endif
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
     %(recvar1(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 42 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 43 lr=1,lamb2
		  do 43 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,9)recvar1(jr,kr)
		     write(84,8)
     %(recvar2(jr,kr)+0.5*(recvar1(jr,kr)**2.))	     !*ros
c		  else
c		     write(83,9)1.e+03
c		     write(84,8)1.e+03
c		  endif
 43		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)recvar1(jy1,kr)
		     write(84,8)
     %(recvar2(jy1,kr)+0.5*(recvar1(jy1,kr)**2.))     !*ros
c		  else
c		     write(83,9)1.e+03
c		     write(84,8)1.e+03
c		  endif
 42		  continue
	       enddo
	    endif
c
!	    write(83,*)	  !!!!!!
!	    write(83,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
!	    write(84,*)	  !!!!!!
!	    write(84,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
!
! VRMS and UVPLOT
!
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)sqrt(abs(vrmsplot(iic,j,k)))
	       write(84,9)uvplot(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 51	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)sqrt(abs(vrmsplot(iic,jy1,k)))
	       write(84,9)uvplot(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 50	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar1(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 52 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 53 lr=1,lamb2
		  do 53 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar1(jr,kr)))
		     write(84,9)recvar2(jr,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 53		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar1(jy1,kr)))
		     write(84,9)recvar2(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 52		  continue
	       enddo
	    endif
!
! URMS and UWPLOT
!
	    do 60 k=kz1,kz2,nzp
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)sqrt(abs(urmsplot(iic,j,k)))
	       write(84,9)uwplot(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 61	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)sqrt(abs(urmsplot(iic,jy1,k)))
	       write(84,9)uwplot(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 60	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar1(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,6,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,7,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 62 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 63 lr=1,lamb2
		  do 63 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar1(jr,kr)))
		     write(84,9)recvar2(jr,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 63		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar1(jy1,kr)))
		     write(84,9)recvar2(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 62		  continue
	       enddo
	    endif
!
! WRMS and VWPLOT
!
	    do 70 k=kz1,kz2,nzp
	    do 71 l=1,lamb2
	    do 71 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)sqrt(abs(wrmsplot(iic,j,k)))
	       write(84,9)vwplot(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 71	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)sqrt(abs(wrmsplot(iic,jy1,k)))
	       write(84,9)vwplot(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 70	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar1(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 72 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 73 lr=1,lamb2
		  do 73 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar1(jr,kr)))
		     write(84,9)recvar2(jr,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 73		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar1(jy1,kr)))
		     write(84,9)recvar2(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 72		  continue
	       enddo
	    endif
!
! Static pressure root mean squares and time-averaged ratio
! eddy viscosity / fluid viscosity
!
	    do 80 k=kz1,kz2,nzp
	    do 81 l=1,lamb2
	    do 81 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,8)sqrt(abs(prmsplot(iic,j,k)))	  !*ros
	       write(84,9)tvplot(iic,j,k)
c	    else
c	       write(83,8)1.e+03
c	       write(84,9)1.e+03
c	    endif
 81	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,8)sqrt(abs(prmsplot(iic,jy1,k)))    !*ros
	       write(84,9)tvplot(iic,jy1,k)
c	    else
c	       write(83,8)1.e+03
c	       write(84,9)1.e+03
c	    endif
 80	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar1(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 82 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 83 lr=1,lamb2
		  do 83 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,8)sqrt(abs(recvar1(jr,kr)))    !*ros
		     write(84,9)recvar2(jr,kr)
c		  else
c		     write(83,8)1.e+03
c		     write(84,9)1.e+03
c		  endif
 83		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,8)sqrt(abs(recvar1(jy1,kr)))    !*ros
		     write(84,9)recvar2(jy1,kr)
c		  else
c		     write(83,8)1.e+03
c		     write(84,9)1.e+03
c		  endif
 82		  continue
	       enddo
	    endif
c	     write(84,*)   !!!!!!
c	     write(84,103)
c     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
	    close(83)
	    close(84)
	 else
c	    CALL MPI_SEND
c    %(flaguo(is,jy1:jy2,kz1:kz2,:),
c    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,0,1,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vmodplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(pplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vrmsplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uvplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(urmsplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uwplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(wrmsplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vwplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(prmsplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(tvplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
      do kkc=1,nsez2max
	 ks=ksez2(kkc)
	 kk=isezindx2+(kkc-1)
	 kkd=kk/10
	 kku=mod(kk,10)
	 open(83,file='results2_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_kplot2.q')
	 open(84,file='results2_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_kplot3.q')
	 write(83,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	 write(84,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	 write(83,102)ntime,i1,i1,i1
	 write(84,102)ntime,i1,i1,i1
!
! Time-averaged velocity magnitude and stagnation pressure
!
	 do 140 i=ix1,ix2,nxp
	 do 141 l=1,lamb2
	 do 141 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)vmodplot2(i,j,kkc)
	    write(84,8)
     %(pplot2(i,j,kkc)+0.5*(vmodplot2(i,j,kkc)**2.))	 !*ros
c	 else
c	    write(83,9)1.e+03
c	    write(84,8)1.e+03
c	 endif
 141	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)vmodplot2(i,jy1,kkc)
	    write(84,8)
     %(pplot2(i,jy1,kkc)+0.5*(vmodplot2(i,jy1,kkc)**2.))     !*ros
c	 else
c	    write(83,9)1.e+03
c	    write(84,8)1.e+03
c	 endif
 140	 continue
!	 write(83,*)   !!!!!!
!	 write(83,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
!	 write(84,*)   !!!!!!
!	 write(84,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
!
! VRMS and UVPLOT
!
	 do 150 i=ix1,ix2,nxp
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)sqrt(abs(vrmsplot2(i,j,kkc)))
	    write(84,9)uvplot2(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 151	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)sqrt(abs(vrmsplot2(i,jy1,kkc)))
	    write(84,9)uvplot2(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 150	 continue
!
! URMS and UWPLOT
!
	 do 160 i=ix1,ix2,nxp
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)sqrt(abs(urmsplot2(i,j,kkc)))
	    write(84,9)uwplot2(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 161	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)sqrt(abs(urmsplot2(i,jy1,kkc)))
	    write(84,9)uwplot2(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 160	 continue
!
! WRMS and VWPLOT
!
	 do 170 i=ix1,ix2,nxp
	 do 171 l=1,lamb2
	 do 171 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)sqrt(abs(wrmsplot2(i,j,kkc)))
	    write(84,9)vwplot2(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 171	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)sqrt(abs(wrmsplot2(i,jy1,kkc)))
	    write(84,9)vwplot2(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 170	 continue
!
! Static pressure root mean squares and time-averaged ratio
! eddy viscosity / fluid viscosity
!
	 do 180 i=ix1,ix2,nxp
	 do 181 l=1,lamb2
	 do 181 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,8)sqrt(abs(prmsplot2(i,j,kkc)))	!*ros
	    write(84,9)tvplot2(i,j,kkc)
c	 else
c	    write(83,8)1.e+03
c	    write(84,9)1.e+03
c	 endif
 181	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,8)sqrt(abs(prmsplot2(i,jy1,kkc)))	  !*ros
	    write(84,9)tvplot2(i,jy1,kkc)
c	 else
c	    write(83,8)1.e+03
c	    write(84,9)1.e+03
c	 endif
 180	 continue
c	  write(84,*)	!!!!!!
c	  write(84,103)
c     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
	 close(83)
	 close(84)
      enddo
c
      jjc=1
      js=jsez14(jjc)
      if(myrank.eq.0) then
	 jjcd=jjc/10
	 jjcu=mod(jjc,10)
	 open(83,file='results3_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_kplot2.q')
	 open(84,file='results3_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_kplot3.q')
	 write(83,101)1,(nx-2)/nxp,(nzg-2)/nzp
	 write(84,101)1,(nx-2)/nxp,(nzg-2)/nzp
	 write(83,102)ntime,i1,i1,i1
	 write(84,102)ntime,i1,i1,i1
!
! Time-averaged velocity magnitude and stagnation pressure
!
	 do 240 k=kz1,kz2,nzp
	 do 240 i=ix1,ix2,nxp
c	 if(all(flagvo(i,js,k,:)).le.0) then
	    write(83,9)vmodplot14(i,jjc,k)
	    write(84,8)
     %(pplot14(i,jjc,k)+0.5*(vmodplot14(i,jjc,k)**2.))	  !*ros
c	 else
c	    write(83,9)1.e+03
c	    write(84,8)1.e+03
c	 endif
 240	 continue
	 if(mysize.ne.1) then
	    do jjp=1,mysize-1
c	       kk1=kz1+jjp*(nz-2)
c	       kk2=kz2+jjp*(nz-2)
c	       CALL MPI_RECV
c    %(flagvor2(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,12,
c    %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV
     %(recvar3(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV
     %(recvar4(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 241 kr=kz1,kz2,nzp
c	       krg=kr+jjp*(nz-2)
	       do 241 ir=ix1,ix2,nxp
c	       if(all(flagvor2(ir,krg,:)).le.0) then
		  write(83,9)recvar3(ir,kr)
		  write(84,8)
     %(recvar4(ir,kr)+0.5*(recvar3(ir,kr)**2.))	   !*ros
c	       else
c		  write(83,9)1.e+03
c		  write(84,8)1.e+03
c	       endif
 241	       continue
	    enddo
	 endif
!	 write(84,*)(((i1,j=1,1),i=1,(nx-2)/nxp),k=1,(nzg-2)/nzp)   !!!!!!
!	 write(84,103)(((i1,j=1,1),i=1,(nx-2)/nxp),k=1,(nzg-2)/nzp)
!
! VRMS and UVPLOT
!
	 do 250 k=kz1,kz2,nzp
	 do 250 i=ix1,ix2,nxp
c	 if(all(flagvo(i,js,k,:)).le.0) then
	    write(83,9)sqrt(abs(vrmsplot14(i,jjc,k)))
	    write(84,9)uvplot14(i,jjc,k)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 250	 continue
	 if(mysize.ne.1) then
	    do jjp=1,mysize-1
	       CALL MPI_RECV
     %(recvar3(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV
     %(recvar4(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,16,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 251 kr=kz1,kz2,nzp
c	       krg=kr+jjp*(nz-2)
	       do 251 ir=ix1,ix2,nxp
c	       if(all(flagvor2(ir,krg,:)).le.0) then
		  write(83,9)sqrt(abs(recvar3(ir,kr)))
		  write(84,9)recvar4(ir,kr)
c	       else
c		  write(83,9)1.e+03
c		  write(84,9)1.e+03
c	       endif
 251	       continue
	    enddo
	 endif
!
! URMS and UWPLOT
!
	 do 260 k=kz1,kz2,nzp
	 do 260 i=ix1,ix2,nxp
c	 if(all(flagvo(i,js,k,:)).le.0) then
	    write(83,9)sqrt(abs(urmsplot14(i,jjc,k)))
	    write(84,9)uwplot14(i,jjc,k)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 260	 continue
	 if(mysize.ne.1) then
	    do jjp=1,mysize-1
	       CALL MPI_RECV
     %(recvar3(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,17,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV
     %(recvar4(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,18,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 261 kr=kz1,kz2,nzp
c	       krg=kr+jjp*(nz-2)
	       do 261 ir=ix1,ix2,nxp
c	       if(all(flagvor2(ir,krg,:)).le.0) then
		  write(83,9)sqrt(abs(recvar3(ir,kr)))
		  write(84,9)recvar4(ir,kr)
c	       else
c		  write(83,9)1.e+03
c		  write(84,9)1.e+03
c	       endif
 261	       continue
	    enddo
	 endif
!
! WRMS and VWPLOT
!
	 do 270 k=kz1,kz2,nzp
	 do 270 i=ix1,ix2,nxp
c	 if(all(flagvo(i,js,k,:)).le.0) then
	    write(83,9)sqrt(abs(wrmsplot14(i,jjc,k)))
	    write(84,9)vwplot14(i,jjc,k)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 270	 continue
	 if(mysize.ne.1) then
	    do jjp=1,mysize-1
	       CALL MPI_RECV
     %(recvar3(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,19,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV
     %(recvar4(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,20,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 271 kr=kz1,kz2,nzp
c		krg=kr+jjp*(nz-2)
	       do 271 ir=ix1,ix2,nxp
c	       if(all(flagvor2(ir,krg,:)).le.0) then
		  write(83,9)sqrt(abs(recvar3(ir,kr)))
		  write(84,9)recvar4(ir,kr)
c	       else
c		  write(83,9)1.e+03
c		  write(84,9)1.e+03
c	       endif
 271	       continue
	    enddo
	 endif
!
! Static pressure root mean squares and time-averaged ratio
! eddy viscosity / fluid viscosity
!
	 do 280 k=kz1,kz2,nzp
	 do 280 i=ix1,ix2,nxp
c	 if(all(flagvo(i,js,k,:)).le.0) then
	    write(83,8)sqrt(abs(prmsplot14(i,jjc,k)))	 !*ros
	    write(84,9)tvplot14(i,jjc,k)
c	 else
c	    write(83,8)1.e+03
c	    write(84,9)1.e+03
c	 endif
 280	 continue
	 if(mysize.ne.1) then
	    do jjp=1,mysize-1
	       CALL MPI_RECV
     %(recvar3(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,21,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV
     %(recvar4(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,22,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 281 kr=kz1,kz2,nzp
c	       krg=kr+jjp*(nz-2)
	       do 281 ir=ix1,ix2,nxp
c	       if(all(flagvor2(ir,krg,:)).le.0) then
		  write(83,8)sqrt(abs(recvar3(ir,kr)))	  !*ros
		  write(84,9)recvar4(ir,kr)
c	       else
c		  write(83,8)1.e+03
c		  write(84,9)1.e+03
c	       endif
 281	       continue
	    enddo
	 endif
c	  write(84,*)(((i1,j=1,1),i=1,(nx-2)/nxp),k=1,(nzg-2)/nzp)  !!!!!!
c	  write(84,103)(((i1,j=1,1),i=1,(nx-2)/nxp),k=1,(nzg-2)/nzp)
	 close(83)
	 close(84)
      else
c	 CALL MPI_SEND
c    %(flagvo(ix1:ix2,js,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,12,
c    %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(vmodplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(pplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(vrmsplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(uvplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,16,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(urmsplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,17,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(uwplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,18,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(wrmsplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,19,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(vwplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,20,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(prmsplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,21,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(tvplot14(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,22,
     %MPI_COMM_EDDY,STATUS,IERR)
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine plotavgvortavg(iplant5c,iplant5d,iplant5p,
     %nx,ny,nz,nzg,nbd,ntime,flaguo,flagvo,flagwo)
c
!     use sections
      use vortavg
      include'common.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer iplant5p,iplant5d,iplant5c,nx,ny,nz,nzg,nbd,ntime
      integer flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),
     %flagwo(nx,ny,nz,nbd)
c
c Local variables
      integer iic,iicd,iicu,jjc,jjcd,jjcu,kkc,kk,
     %kkd,kku,i1,i,j,k,l,ir,jr,kr,lr,jjp
      integer STATUS(MPI_STATUS_SIZE)
!     integer is,js,ks,kk1,kk2,krg,flaguor(ny,nzg,nbd),
!    %flagvor2(nx,nzg,nbd)
      real recvar(ny,nz),recvar2(nx,nz)
c
c Parameters
      integer lamb2,nxp,nyp,nzp
      parameter (lamb2=1)
      parameter (nyp=1,nxp=1,nzp=1)
c
   8  format(10f25.8)
   9  format(10f25.8)
 101  format(3i8)
 102  format(i8,3i2)
 103  format(i2)
c
      i1=1
      do iic=1,nsez1
!	 is=isez1(iic)
	 if(myrank.eq.0) then
	    iicd=iic/10
	    iicu=mod(iic,10)
	    open(83,file='results1_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'_avgvort.q')
	    write(83,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	    write(83,102)ntime,i1,i1,i1
!!!!!!		  write(83,*)
	    write(83,103)
     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
	    do 40 k=kz1,kz2,nzp
	    do 41 l=1,lamb2
	    do 41 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)voraz1med(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
 41	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)voraz1med(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)recvar(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
 42		  continue
	       enddo
	    endif
c
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)vorr1med(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
 51	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)vorr1med(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)recvar(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
 52		  continue
	       enddo
	    endif
c
	    do 60 k=kz1,kz2,nzp
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)vorz1med(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
 61	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)vorz1med(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)recvar(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
 62		  continue
	       enddo
	    endif
!!!!!!		  write(83,*)
	    write(83,103)
     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
	    close(83)
	 else
c	    CALL MPI_SEND
c    %(flaguo(is,jy1:jy2,kz1:kz2,:),
c    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,0,1,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(voraz1med(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorr1med(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorz1med(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
      do kkc=1,nsez2max
!	 ks=ksez2(kkc)
	 kk=isezindx2+(kkc-1)
	 kkd=kk/10
	 kku=mod(kk,10)
	 open(83,file='results2_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_avgvort.q')
	 write(83,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	 write(83,102)ntime,i1,i1,i1
!!!!!!	       write(83,*)
	 write(83,103)
     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
	 do 140 i=ix1,ix2,nxp
	 do 141 l=1,lamb2
	 do 141 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)voraz2med(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
 141	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)voraz2med(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
 140	 continue
	 do 150 i=ix1,ix2,nxp
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)vorr2med(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
 151	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)vorr2med(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
 150	 continue
	 do 160 i=ix1,ix2,nxp
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)vorz2med(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
 161	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)vorz2med(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
 160	 continue
!!!!!!	       write(83,*)
	 write(83,103)
     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
	 close(83)
      enddo
c
      jjc=1
!     js=jsez14(jjc)
      if(myrank.eq.0) then
	 jjcd=jjc/10
	 jjcu=mod(jjc,10)
	 open(83,file='results3_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_avgvort.q')
	 write(83,101)1,(nx-2)/nxp,(nzg-2)/nzp
	 write(83,102)ntime,i1,i1,i1
!!!!!!	       write(83,*)(((i1,j=1,1),i=1,(nx-2)/nxp),k=1,(nzg-2)/nzp)
	 write(83,103)(((i1,j=1,1),i=1,(nx-2)/nxp),k=1,(nzg-2)/nzp)
	 do 240 k=kz1,kz2,nzp
	 do 240 i=ix1,ix2,nxp
c	 if(all(flagvo(i,js,k,:)).le.0) then
	    write(83,9)voraz14med(i,jjc,k)
c	 else
c	    write(83,9)1.e+03
c	 endif
 240	 continue
	 if(mysize.ne.1) then
	    do jjp=1,mysize-1
c	       kk1=kz1+jjp*(nz-2)
c	       kk2=kz2+jjp*(nz-2)
c	       CALL MPI_RECV
c    %(flagvor2(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,5,
c    %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 241 kr=kz1,kz2,nzp
c	       krg=kr+jjp*(nz-2)
	       do 241 ir=ix1,ix2,nxp
c	       if(all(flagvor2(ir,krg,:)).le.0) then
		  write(83,9)recvar2(ir,kr)
c	       else
c		  write(83,9)1.e+03
c	       endif
 241	       continue
	    enddo
	 endif
c
	 do 250 k=kz1,kz2,nzp
	 do 250 i=ix1,ix2,nxp
c	 if(all(flagvo(i,js,k,:)).le.0) then
	    write(83,9)vorr14med(i,jjc,k)
c	 else
c	    write(83,9)1.e+03
c	 endif
 250	 continue
	 if(mysize.ne.1) then
	    do jjp=1,mysize-1
	       CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 251 kr=kz1,kz2,nzp
c	       krg=kr+jjp*(nz-2)
	       do 251 ir=ix1,ix2,nxp
c	       if(all(flagvor2(ir,krg,:)).le.0) then
		  write(83,9)recvar2(ir,kr)
c	       else
c		  write(83,9)1.e+03
c	       endif
 251	       continue
	    enddo
	 endif
c
	 do 260 k=kz1,kz2,nzp
	 do 260 i=ix1,ix2,nxp
c	 if(all(flagvo(i,js,k,:)).le.0) then
	    write(83,9)vorz14med(i,jjc,k)
c	 else
c	    write(83,9)1.e+03
c	 endif
 260	 continue
	 if(mysize.ne.1) then
	    do jjp=1,mysize-1
	       CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 261 kr=kz1,kz2,nzp
c	       krg=kr+jjp*(nz-2)
	       do 261 ir=ix1,ix2,nxp
c	       if(all(flagvor2(ir,krg,:)).le.0) then
		  write(83,9)recvar2(ir,kr)
c	       else
c		  write(83,9)1.e+03
c	       endif
 261	       continue
	    enddo
	 endif
!!!!!!	       write(83,*)(((i1,j=1,1),i=1,(nx-2)/nxp),k=1,(nzg-2)/nzp)
	 write(83,103)(((i1,j=1,1),i=1,(nx-2)/nxp),k=1,(nzg-2)/nzp)
	 close(83)
      else
c	 CALL MPI_SEND
c    %(flagvo(ix1:ix2,js,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,5,
c    %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND 
     %(voraz14med(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(vorr14med(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND
     %(vorz14med(ix1:ix2,jjc,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
      endif
c
      return
      end
c
C---- subroutine read_probes ----------------N. Beratlis-10 May 2010---
c------------------------------------------------A. Posa - Dec 2010----
C
C     PURPOSE: Setup probes.
C
C------------------------------------------------------------------------
      subroutine read_probes(nz)
c
      use points
      include 'common.h'
      include 'averages.h'
c
c.... Input/Output Arrays
      integer nz
c
c.... Local arrays
      integer i,ii,jj,kk,k,kk1,kk2,
     %iprb2g(nprobes),jprb2g(nprobes),kprb2g(nprobes)
c
      open(unit=10,file='probes.input',form='formatted')
      do i=1,nprobes
	read(10,*) iprb2g(i),jprb2g(i),kprb2g(i)
      if(MYRANK==0)write(*,*) iprb2g(i),jprb2g(i),kprb2g(i)
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
      enddo
c
      return
      end
c
c---- subroutine sortprobes2 --------------N. Beratlis-06 May 2010---
c---------------------------------------------A. Posa - Dec 2010-----
C
C     PURPOSE: Sort probes in increasing k index.
C
C------------------------------------------------------------------------
      subroutine sortprobes2(ip,jp,kp,n)
c
      implicit none
c
c...  Input/Output Arrays
      integer n
      integer ip(n),jp(n),kp(n)
c
c.... Local arrays
      integer i,imin,temp
      integer i2(n),j2(n),k2(n)

      i2 = ip
      j2 = jp
      k2 = kp

      do i=1,n
	imin = i-1 + minloc(k2(i:n),1)
	temp = k2(i)
	k2(i) = k2(imin)
	k2(imin) = temp

	temp = i2(i)
	i2(i) = i2(imin)
	i2(imin) = temp

	temp = j2(i)
	j2(i) = j2(imin)
	j2(imin) = temp

      enddo

      kp = k2
      jp = j2
      ip = i2

      return

      end
C--------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sondestagg(nx,ny,nz,nbd,time,dt,vo,uo,wo,p)    !,flagpo)
c
      use points
      implicit none
      include 'averages.h'
c
c Global variables
      integer nx,ny,nz,nbd
!      integer flagpo(nx,ny,nz,nbd)
      real time,dt
      real vo(nx,ny,nz),uo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
c
c Local variables
      integer ii,i,ic,ip,jc,jp,kc,kp
      real vstagg,ustagg,wstagg,pstagg
c
 10   format(5f25.8)
 11   format(9f25.8)
!
!     Instantaneous values, time averages and root mean squares
!
      timepoints=timepoints+dt
      do ii=1,nprbmax2
	 i=prbindx2+(ii-1)
	 ic=iprb2(ii)
	 ip=ic+1
	 jc=jprb2(ii)
	 jp=jpv(jc)
	 kc=kprb2(ii)
	 kp=kc+1
	 vstagg=0.25*(vo(ic,jc,kc)+vo(ip,jc,kc)+
     %vo(ic,jc,kp)+vo(ip,jc,kp))
	 ustagg=0.25*(uo(ic,jc,kc)+uo(ic,jp,kc)+
     %uo(ic,jc,kp)+uo(ic,jp,kp))
	 wstagg=0.25*(wo(ic,jc,kc)+wo(ip,jc,kc)+
     %wo(ic,jp,kc)+wo(ip,jp,kc))
	 pstagg=0.125*(p(ic,jc,kc)+p(ip,jc,kc)+
     %p(ic,jp,kc)+p(ip,jp,kc)+
     %p(ic,jc,kp)+p(ip,jc,kp)+
     %p(ic,jp,kp)+p(ip,jp,kp))
         if((time-timep).le.periodo) goto 100
!         if(any(flagpo(ic:ip,jc:jp,kc:kp,:).ne.0)) goto 100
	 avrpoints1(ii,1)=avrpoints1(ii,1)+dt*vstagg
	 avrpoints1(ii,2)=avrpoints1(ii,2)+dt*ustagg
	 avrpoints1(ii,3)=avrpoints1(ii,3)+dt*wstagg
	 avrpoints1(ii,4)=avrpoints1(ii,4)+dt*pstagg
	 rmspoints1(ii,1)=rmspoints1(ii,1)+dt*(vstagg**2.)
	 rmspoints1(ii,2)=rmspoints1(ii,2)+dt*(ustagg**2.)
	 rmspoints1(ii,3)=rmspoints1(ii,3)+dt*(wstagg**2.)
	 rmspoints1(ii,4)=rmspoints1(ii,4)+dt*(pstagg**2.)
	 timepoints1(ii)=timepoints1(ii)+dt
 100	 continue
	 write(700+i,10)time,vstagg,ustagg,wstagg,pstagg   !*ros
      enddo
c
      if(timepoints.ge.periodo3) then
!
!     Time averages and root mean squares on file
!
	 do ii=1,nprbmax2
	    i=prbindx2+(ii-1)
	    avrpoints1(ii,1)=avrpoints1(ii,1)/timepoints1(ii)
	    avrpoints1(ii,2)=avrpoints1(ii,2)/timepoints1(ii)
	    avrpoints1(ii,3)=avrpoints1(ii,3)/timepoints1(ii)
	    avrpoints1(ii,4)=avrpoints1(ii,4)/timepoints1(ii)
	    rmspoints1(ii,1)=
     %sqrt(rmspoints1(ii,1)/timepoints1(ii)-avrpoints1(ii,1)**2.)
	    rmspoints1(ii,2)=
     %sqrt(rmspoints1(ii,2)/timepoints1(ii)-avrpoints1(ii,2)**2.)
	    rmspoints1(ii,3)=
     %sqrt(rmspoints1(ii,3)/timepoints1(ii)-avrpoints1(ii,3)**2.)
	    rmspoints1(ii,4)=
     %sqrt(rmspoints1(ii,4)/timepoints1(ii)-avrpoints1(ii,4)**2.)
	    write(800+i,11)time-periodo3*0.5,
     %avrpoints1(ii,1),rmspoints1(ii,1),
     %avrpoints1(ii,2),rmspoints1(ii,2),
     %avrpoints1(ii,3),rmspoints1(ii,3),
     %avrpoints1(ii,4),rmspoints1(ii,4)	   !*ros     !*ros
	    avrpoints1(ii,1)=0.
	    avrpoints1(ii,2)=0.
	    avrpoints1(ii,3)=0.
	    avrpoints1(ii,4)=0.
	    rmspoints1(ii,1)=0.
	    rmspoints1(ii,2)=0.
	    rmspoints1(ii,3)=0.
	    rmspoints1(ii,4)=0.
	    timepoints1(ii)=0.
	 enddo
	 timepoints=0.
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine continua(nx,ny,nz)
c
      include 'common.h'
      include 'averages.h'
c
      integer nx,ny,nz
c
! 999  format(4f20.10)
      if(myrank.eq.0) then
         open(81,file='tempi.res',form='unformatted')
         write(81)tempo,tempo3,timepoints,timep
         write(81)iplant4,iplant5,iplant8,iplant9
         close(81)
         return
      endif
      return
      CALL IOMPI_3DSCALAR_mediey
     %('mediey.res',ny,1)
      CALL IOMPI_3DSCALAR_mediex
     %('mediex.res',nx,1)
      CALL IOMPI_3DSCALAR_mediez
     %('mediez.res',nz,1)
      CALL IOMPI_3DSCALAR_mediecalc1
     %('mediecalc1.res',ny,nz,1)
      CALL IOMPI_3DSCALAR_mediecalc2
     %('mediecalc2.res',nx,ny,1)
      CALL IOMPI_3DSCALAR_mediecalc3
     %('mediecalc3.res',nx,nz,1)
      CALL IOMPI_3DSCALAR_sonde
     %('sonde.res',1)
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine inirea(nx,ny,nz)
c
      include 'common.h'
      include 'averages.h'
      include 'mpif.h'
c
      integer nx,ny,nz
c
! 999  format(4f20.10)
!      if(myrank.eq.0) then
         open(81,file='res.tempi',form='unformatted')
         read(81,end=100)tempo,tempo3,timepoints,timep
         read(81)iplant4,iplant5,iplant8,iplant9
         close(81)
         return
!      endif
!      CALL MPI_BCAST(tempo,1,MTYPE,0,MPI_COMM_EDDY,IERR)
!      CALL MPI_BCAST(tempo3,1,MTYPE,0,MPI_COMM_EDDY,IERR)
!      CALL MPI_BCAST(timepoints,1,MTYPE,0,MPI_COMM_EDDY,IERR)
!      CALL MPI_BCAST(timep,1,MTYPE,0,MPI_COMM_EDDY,IERR)
!      CALL MPI_BCAST(iplant4,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
!      CALL MPI_BCAST(iplant5,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
!      CALL MPI_BCAST(iplant8,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
!      CALL MPI_BCAST(iplant9,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)
      CALL IOMPI_3DSCALAR_mediey
     %('res.mediey',ny,0)
      CALL IOMPI_3DSCALAR_mediex
     %('res.mediex',nx,0)
      CALL IOMPI_3DSCALAR_mediez
     %('res.mediez',nz,0)
      CALL IOMPI_3DSCALAR_mediecalc1
     %('res.mediecalc1',ny,nz,0)
      CALL IOMPI_3DSCALAR_mediecalc2
     %('res.mediecalc2',nx,ny,0)
      CALL IOMPI_3DSCALAR_mediecalc3
     %('res.mediecalc3',nx,nz,0)
      CALL IOMPI_3DSCALAR_sonde
     %('res.sonde',0)
c
      return
 100  continue
      close(81)
!      call setup_output(nx,ny,nz)
      return
      end
c
C---- subroutine iompi_3dscalar_medie-------------------------------
C
      subroutine iompi_3dscalar_medie(filename,nx,ny,nz,io)
C
C     PURPOSE: Read from or write to a file a 3D real array using
C     MPI subroutines.
C
C-------------------------------------------------------------------
      use mediey
      use mediex
      use mediez
      use rmsy
      use rmsx
      use rmsz
      include 'common.h'
      include 'averages.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,io
      character*(*) filename
c
c.... Local arrays
      integer fh,filemode
      integer newtype,ii,cost0,cost1,cost2
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      if(io==1) then

	filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

	call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %(nsez3*(ny-2)*nsez4max,mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=nsez3*(ny-2)*nsez4
	cost2=nsez3*(ny-2)*max(0,isezindx4-1)

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,tempo1(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,uavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,wavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vmodavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vorazavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vorravg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vorzavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,pavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,ptotavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vrms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,urms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,wrms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,prms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,uvs(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,uws(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vws(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=18
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,tvavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_type_contiguous
     %((nx-2)*nsez5*nsez6max,mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	cost0=ii*nsez3*(ny-2)*nsez4
	cost1=(nx-2)*nsez5*nsez6
	cost2=(nx-2)*nsez5*max(0,isezindx6-1)

	ii=1
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,tempo2(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=2
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=3
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,ruavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=4
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rwavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=5
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvmodavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=6
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rpavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=7
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rptotavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=8
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvorazavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=9
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvorravg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=10
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvorzavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=11
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvrms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=12
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rurms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=13
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rwrms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=14
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rprms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=15
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,ruvs(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=16
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,ruws(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=17
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvws(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=18
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rtvavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_type_contiguous
     %(nsez9*nsez91*(nz-2),mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	cost0=ii*nsez3*(ny-2)*nsez4+ii*(nx-2)*nsez5*nsez6
	cost1=nsez9*nsez91*(nz-2)*mysize
	cost2=nsez9*nsez91*(nz-2)*myrank

	ii=1
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=2
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zuavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=3
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zwavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=4
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvmodavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=5
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zpavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=6
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zptotavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=7
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvorazavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=8
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvorravg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=9
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvorzavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=10
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvrms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=11
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zurms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=12
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zwrms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=13
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zprms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=14
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zuvs(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=15
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zuws(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=16
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvws(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=17
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,ztempo(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=18
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,ztvavg(:,:,kz1:kz2),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      else

	filemode = MPI_MODE_RDONLY

	call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %(nsez3*(ny-2)*nsez4max,mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=nsez3*(ny-2)*nsez4
	cost2=nsez3*(ny-2)*max(0,isezindx4-1)

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,tempo1(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,uavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,wavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vmodavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vorazavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vorravg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vorzavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,pavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,ptotavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vrms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,urms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,wrms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,prms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,uvs(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,uws(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vws(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=18
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,tvavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_type_contiguous
     %((nx-2)*nsez5*nsez6max,mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	cost0=ii*nsez3*(ny-2)*nsez4
	cost1=(nx-2)*nsez5*nsez6
	cost2=(nx-2)*nsez5*max(0,isezindx6-1)

	ii=1
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,tempo2(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=2
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=3
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,ruavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=4
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rwavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=5
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvmodavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=6
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rpavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=7
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rptotavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=8
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvorazavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=9
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvorravg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=10
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvorzavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=11
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvrms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=12
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rurms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=13
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rwrms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=14
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rprms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=15
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,ruvs(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=16
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,ruws(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=17
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvws(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=18
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rtvavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_type_contiguous
     %(nsez9*nsez91*(nz-2),mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	cost0=ii*nsez3*(ny-2)*nsez4+ii*(nx-2)*nsez5*nsez6
	cost1=nsez9*nsez91*(nz-2)*mysize
	cost2=nsez9*nsez91*(nz-2)*myrank

	ii=1
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=2
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zuavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=3
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zwavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=4
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvmodavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=5
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zpavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=6
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zptotavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=7
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvorazavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=8
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvorravg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=9
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvorzavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=10
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvrms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=11
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zurms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=12
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zwrms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=13
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zprms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=14
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zuvs(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=15
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zuws(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=16
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvws(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=17
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,ztempo(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=18
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,ztvavg(:,:,kz1:kz2),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      endif

      return

      end
c
C---- subroutine iompi_3dscalar_mediecalc1---------------------------
C
      subroutine iompi_3dscalar_mediecalc1(filename,ny,nz,io)
C
C     PURPOSE: Read from or write to a file a 3D real array using
C     MPI subroutines.
C
C-------------------------------------------------------------------
      use mediecalc
      use mediecalc1
      use vortavg
      include 'common.h'
      include 'averages.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer ny,nz,io
      character*(*) filename
c
c.... Local arrays
      integer fh,filemode
      integer newtype,ii,cost1,cost2
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      if(io==1) then

	filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

	call mpi_file_open
     %(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %(nsez1*(ny-2)*(nz-2),mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view
     %(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=nsez1*(ny-2)*(nz-2)*mysize
	cost2=nsez1*(ny-2)*(nz-2)*myrank

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,tempo4(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,uplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,wplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,pplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vrmsplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,urmsplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,wrmsplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,prmsplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,uvplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,uwplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vwplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,voraz1med(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vorr1med(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vorz1med(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vmodplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,tvplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      else

	filemode = MPI_MODE_RDONLY

	call mpi_file_open
     %(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %(nsez1*(ny-2)*(nz-2),mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view
     %(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=nsez1*(ny-2)*(nz-2)*mysize
	cost2=nsez1*(ny-2)*(nz-2)*myrank

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,tempo4(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,uplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,wplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,pplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vrmsplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,urmsplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,wrmsplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,prmsplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,uvplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,uwplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vwplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,voraz1med(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vorr1med(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vorz1med(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vmodplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,tvplot(:,jy1:jy2,kz1:kz2),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      endif

      return

      end
*-------------------------------------------------------------------
C---- subroutine iompi_3dscalar_mediecalc2---------------------------
C
      subroutine iompi_3dscalar_mediecalc2(filename,nx,ny,io)
C
C     PURPOSE: Read from or write to a file a 3D real array using
C     MPI subroutines.
C
C-------------------------------------------------------------------
      use mediecalc
      use mediecalc2
      use vortavg
      include 'common.h'
      include 'averages.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,io
      character*(*) filename
c
c.... Local arrays
      integer fh,filemode
      integer newtype,ii,cost1,cost2
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      if(io==1) then

	filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

	call mpi_file_open
     %(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %((nx-2)*(ny-2)*nsez2max,mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view
     %(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=(nx-2)*(ny-2)*nsez2
	cost2=(nx-2)*(ny-2)*max(0,isezindx2-1)

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,tempo5(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,uplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,wplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,pplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vrmsplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,urmsplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,wrmsplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,prmsplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,uvplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,uwplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vwplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,voraz2med(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vorr2med(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vorz2med(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vmodplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,tvplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      else

	filemode = MPI_MODE_RDONLY

	call mpi_file_open
     %(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %((nx-2)*(ny-2)*nsez2max,mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view
     %(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=(nx-2)*(ny-2)*nsez2
	cost2=(nx-2)*(ny-2)*max(0,isezindx2-1)

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,tempo5(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,uplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,wplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,pplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vrmsplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,urmsplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,wrmsplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,prmsplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,uvplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,uwplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vwplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,voraz2med(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vorr2med(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vorz2med(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vmodplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,tvplot2(ix1:ix2,jy1:jy2,1:nsez2max),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      endif

      return

      end
*-------------------------------------------------------------------
C---- subroutine iompi_3dscalar_mediecalc3---------------------------
C
      subroutine iompi_3dscalar_mediecalc3(filename,nx,nz,io)
C
C     PURPOSE: Read from or write to a file a 3D real array using
C     MPI subroutines.
C
C-------------------------------------------------------------------
      use mediecalc
      use mediecalc3
      use vortavg
      include 'common.h'
      include 'averages.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,nz,io
      character*(*) filename
c
c.... Local arrays
      integer fh,filemode
      integer newtype,ii,cost1,cost2
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      if(io==1) then

	filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

	call mpi_file_open
     %(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %((nx-2)*2*nsez14*(nz-2),mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view
     %(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=(nx-2)*2*nsez14*(nz-2)*mysize
	cost2=(nx-2)*2*nsez14*(nz-2)*myrank

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,tempo14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,uplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,wplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,pplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vrmsplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,urmsplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,wrmsplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,prmsplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,uvplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,uwplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vwplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,voraz14med(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vorr14med(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vorz14med(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,vmodplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all
     %(fh,offset,tvplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      else

	filemode = MPI_MODE_RDONLY

	call mpi_file_open
     %(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %((nx-2)*2*nsez14*(nz-2),mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view
     %(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=(nx-2)*2*nsez14*(nz-2)*mysize
	cost2=(nx-2)*2*nsez14*(nz-2)*myrank

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,tempo14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,uplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,wplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,pplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vrmsplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,urmsplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,wrmsplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,prmsplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,uvplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,uwplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vwplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,voraz14med(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vorr14med(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vorz14med(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,vmodplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all
     %(fh,offset,tvplot14(ix1:ix2,:,kz1:kz2),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      endif

      return

      end
*-------------------------------------------------------------------
C---- subroutine iompi_3dscalar_sonde-------------------------------
C
      subroutine iompi_3dscalar_sonde(filename,io)
C
C     PURPOSE: Read from or write to a file a 3D real array using
C     MPI subroutines.
C
C-------------------------------------------------------------------
      use points
      include 'common.h'
      include 'averages.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer io
      character*(*) filename
c
c.... Local arrays
      integer fh,filemode
      integer newtype,i
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      if(io==1) then

	filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

	call mpi_file_open
     %(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %(nprbmax2,mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view
     %(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	do 100 i=1,4

	offset = max(0,prbindx2-1) + (i-1)*nprobes

	call mpi_file_write_at_all
     %(fh,offset,avrpoints1(1:nprbmax2,i),1,newtype,status,ierr)

	offset = offset + 4*nprobes

	call mpi_file_write_at_all
     %(fh,offset,rmspoints1(1:nprbmax2,i),1,newtype,status,ierr)

  100	continue

	offset = max(0,prbindx2-1) + 8*nprobes

	call mpi_file_write_at_all
     %(fh,offset,timepoints1(1:nprbmax2),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      else

	filemode = MPI_MODE_RDONLY

	call mpi_file_open
     %(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %(nprbmax2,mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view
     %(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	do 200 i=1,4

	offset = max(0,prbindx2-1) + (i-1)*nprobes

	call mpi_file_read_at_all
     %(fh,offset,avrpoints1(1:nprbmax2,i),1,newtype,status,ierr)

	offset = offset + 4*nprobes

	call mpi_file_read_at_all
     %(fh,offset,rmspoints1(1:nprbmax2,i),1,newtype,status,ierr)

  200	continue

	offset = max(0,prbindx2-1) + 8*nprobes

	call mpi_file_read_at_all
     %(fh,offset,timepoints1(1:nprbmax2),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      endif

      return

      end
*-------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vort1staggnoavg(uo,vo,wo,nx,ny,nz,xc,xu,zc)
c
      use sections
      use mediey
      use mediex
      use vorticity
      include 'common.h'
      include 'averages.h'
c
c.... Input/Output array
      integer nx,ny,nz
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
      real xc(nx),xu(nx),zc(nz)
c
c.... Local arrays
      integer kkc,kc,kp,km,iic,ic,ip,im,jjc,jc,jp,jm,jsy,jpsy
      real dz,dtheta,dr,dq1x3,dq3x1,dq2x3,dq3x2,dq2x1,dq1x2
c
!     if(icyl.eq.1) then
c
c    ! Azimuthal vorticity - Distribution along y !
c
      do 101 kkc=1,nsez4max
	kc=ksez4(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
	do 101 iic=1,nsez3
	  ic=isez3(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 101 jc=jy1,jy2
	    dq1x3=(uo(ic,jc,kp)-uo(ic,jc,kc))/dz
	    dq3x1=(wo(ip,jc,kc)-wo(ic,jc,kc))/dr
	    voraz(iic,jc,kkc)=dq1x3-dq3x1
  101 continue
c
c    ! Azimuthal vorticity - Distribution along x !
c
      do 111 kkc=1,nsez6max
	kc=ksez6(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
	ic=ix1
	ip=ic+1
	im=ic-1
	dr=xc(ip)-xc(im)
	do 112 jjc=1,nsez5
	   jc=jsez5(jjc)
	   jp=jpv(jc)
	   jsy=jsym(jc)
	   jpsy=jsym(jp)
	   dq1x3=(0.25*(uo(ic,jc,kp)+uo(im,jc,kp)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,kc)+uo(im,jc,kc)
     %+uo(ic,jp,kc)+uo(im,jp,kc)))/dz
	   dq3x1=(0.5*(wo(ip,jc,kc)+wo(ip,jp,kc))
     %-0.5*(wo(ic,jsy,kc)+wo(ic,jpsy,kc)))/dr
	   rvoraz(ic,jjc,kkc)=(dq1x3-dq3x1)
  112 continue
	do 113 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 113 jjc=1,nsez5
	     jc=jsez5(jjc)
	     jp=jpv(jc)
	     dq1x3=(0.25*(uo(ic,jc,kp)+uo(im,jc,kp)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,kc)+uo(im,jc,kc)
     %+uo(ic,jp,kc)+uo(im,jp,kc)))/dz
	     dq3x1=(0.5*(wo(ip,jc,kc)+wo(ip,jp,kc))
     %-0.5*(wo(im,jc,kc)+wo(im,jp,kc)))/dr
	     rvoraz(ic,jjc,kkc)=(dq1x3-dq3x1)
  113 continue
  111 continue
c
c    ! Radial vorticity - Distribution along y !
c
      do 121 kkc=1,nsez4max
	kc=ksez4(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
	do 121 iic=1,nsez3
	  ic=isez3(iic)
	  ip=ic+1
	  do 121 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq3x2=(0.5*(wo(ip,jp,kc)+wo(ic,jp,kc))
     %-0.5*(wo(ip,jm,kc)+wo(ic,jm,kc)))/dtheta
	    dq2x3=(0.25*(vo(ip,jc,kp)+vo(ic,jc,kp)
     %+vo(ip,jm,kp)+vo(ic,jm,kp))
     %-0.25*(vo(ip,jc,kc)+vo(ic,jc,kc)
     %+vo(ip,jm,kc)+vo(ic,jm,kc)))/dz
	    vorr(iic,jc,kkc)=dq3x2/xu(ic)-dq2x3
  121 continue
c
c    ! Radial vorticity - Distribution along x !
c
      do 131 kkc=1,nsez6max
	kc=ksez6(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
	do 131 ic=ix1,ix2
	  do 131 jjc=1,nsez5
	    jc=jsez5(jjc)
	    jp=jpv(jc)
	    dtheta=dely
	    dq3x2=(wo(ic,jp,kc)-wo(ic,jc,kc))/dtheta
	    dq2x3=(vo(ic,jc,kp)-vo(ic,jc,kc))/dz
	    rvorr(ic,jjc,kkc)=(dq3x2/xc(ic)-dq2x3)
  131 continue
c
c    ! Axial vorticity - Distribution along y !
c
      do 141 kkc=1,nsez4max
	kc=ksez4(kkc)
	kp=kc+1
	do 141 iic=1,nsez3
	  ic=isez3(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 141 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.25*(vo(ip,jc,kp)+vo(ip,jm,kp)
     %+vo(ip,jc,kc)+vo(ip,jm,kc))*xc(ip)
     %-0.25*(vo(ic,jc,kp)+vo(ic,jm,kp)
     %+vo(ic,jc,kc)+vo(ic,jm,kc))*xc(ic))/dr
	    dq1x2=(0.5*(uo(ic,jp,kp)+uo(ic,jp,kc))
     %-0.5*(uo(ic,jm,kp)+uo(ic,jm,kc)))/dtheta
	    vorz(iic,jc,kkc)=(dq2x1-dq1x2)/xu(ic)
  141 continue
c
c    ! Axial vorticity - Distribution along x !
c
      do 151 kkc=1,nsez6max
	kc=ksez6(kkc)
	kp=kc+1
	ic=ix1
	ip=ic+1
	im=ic-1
	dr=xu(ic)-xu(im)
	do 152 jjc=1,nsez5
	   jc=jsez5(jjc)
	   jp=jpv(jc)
	   dtheta=dely
	   dq2x1=(0.25*(vo(ip,jc,kp)+vo(ip,jc,kc)
     %+vo(ic,jc,kp)+vo(ic,jc,kc))*xu(ic))/dr
	   dq1x2=(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,kc)+uo(im,jc,kc)
     %+uo(ic,jc,kp)+uo(im,jc,kp)))/dtheta
	    rvorz(ic,jjc,kkc)=(dq2x1-dq1x2)/xc(ic)
  152	continue
c
	do 153 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 153 jjc=1,nsez5
	    jc=jsez5(jjc)
	    jp=jpv(jc)
	    dtheta=dely
	    dq2x1=(0.5*(vo(ip,jc,kc)+vo(ip,jc,kp))*xc(ip)
     %-0.5*(vo(im,jc,kc)+vo(im,jc,kp))*xc(im))/dr
	    dq1x2=
     %(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,kc)+uo(im,jc,kc)
     %+uo(ic,jc,kp)+uo(im,jc,kp)))/dtheta
	    rvorz(ic,jjc,kkc)=(dq2x1-dq1x2)/xc(ic)
  153 continue
  151 continue
c
c    ! Azimuthal vorticity - Distribution along z !
c
      do 161 kc=kz1,kz2
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	do 161 iic=1,nsez9
	  ic=isez9(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 161 jjc=1,nsez91
	    jc=jsez91(jjc)
	    jp=jpv(jc)
	    dq1x3=(0.5*(uo(ic,jc,kp)+uo(ic,jp,kp))
     %-0.5*(uo(ic,jc,km)+uo(ic,jp,km)))/dz
	    dq3x1=(0.25*(wo(ip,jc,kc)+wo(ip,jp,kc)
     %+wo(ip,jc,km)+wo(ip,jp,km))
     %-0.25*(wo(ic,jc,kc)+wo(ic,jp,kc)
     %+wo(ic,jc,km)+wo(ic,jp,km)))/dr
	    voraz1(iic,jjc,kc)=(dq1x3-dq3x1)
  161 continue
c
c    ! Radial vorticity - Distribution along z !
c
      do 171 kc=kz1,kz2
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	do 171 iic=1,nsez9
	  ic=isez9(iic)
	  ip=ic+1
	  do 171 jjc=1,nsez91
	    jc=jsez91(jjc)
	    jp=jpv(jc)
	    dtheta=dely
	    dq3x2=(0.25*(wo(ic,jp,kc)+wo(ip,jp,kc)
     %+wo(ic,jp,km)+wo(ip,jp,km))
     %-0.25*(wo(ic,jc,kc)+wo(ip,jc,kc)
     %+wo(ic,jc,km)+wo(ip,jc,km)))/dtheta
	    dq2x3=(0.5*(vo(ic,jc,kp)+vo(ip,jc,kp))
     %-0.5*(vo(ic,jc,km)+vo(ip,jc,km)))/dz
	    vorr1(iic,jjc,kc)=(dq3x2/xu(ic)-dq2x3)
  171 continue
c
c    ! Axial vorticity - Distribution along z !
c
      do 181 kc=kz1,kz2
	do 181 iic=1,nsez9
	  ic=isez9(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 181 jjc=1,nsez91
	    jc=jsez91(jjc)
	    jp=jpv(jc)
	    dtheta=dely
	    dq2x1=(vo(ip,jc,kc)*xc(ip)
     %-vo(ic,jc,kc)*xc(ic))/dr
	    dq1x2=(uo(ic,jp,kc)
     %-uo(ic,jc,kc))/dtheta
	    vorz1(iic,jjc,kc)=(dq2x1-dq1x2)/xu(ic)
  181 continue
c
!     endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine mediestaggnoavg(nx,ny,nz,uo,vo,wo,p,tv,dt,nbd,flagpo,time)
c
      use sections
      use mediey
      use mediex
      use mediez
      use vorticity
      include'common.h'
      include'averages.h'
      include'mpif.h'
c
c Global variables
      integer nx,ny,nz,nbd,flagpo(nx,ny,nz,nbd)
      real dt,time
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),
     %p(nx,ny,nz),tv(nx,ny,nz)
c
c Local variables
      integer iic,i,is,ip,im,jjc,j,js,jp,jm,jsy,jpsy,jmsy,
     %kkc,k,ks,kp,km,mm1,mm2,kg,nsez4maxr,nsez6maxr
      integer, dimension(:), allocatable :: ksez4r,ksez6r
      integer STATUS(MPI_STATUS_SIZE)
      real ustagg,vstagg,wstagg,pstagg,vmod,tvstagg
      real, dimension(:,:,:), allocatable :: vmodavgr,uavgr,vavgr,
     %wavgr,pavgr,ptotavgr,vorravgr,vorazavgr,vorzavgr,tvavgr,
     %rvmodavgr,ruavgr,rvavgr,rwavgr,rpavgr,rptotavgr,rvorravgr,
     %rvorazavgr,rvorzavgr,rtvavgr
      real zuavgr(nz),zvavgr(nz),zwavgr(nz),zpavgr(nz)
c
  10  format(4i6,10f25.8)
c
      tempo=tempo+dt
!
!     Time averages at some (x,z) positions along y
!
      do 101 iic=1,nsez3
      is=isez3(iic)
      ip=is+1
      if(is.eq.1) then
	do 1010 kkc=1,nsez4max
	ks=ksez4(kkc)
	kp=ks+1
	do 1011 j=jy1,jy2
	jm=jmv(j)
	jsy=jsym(j)
	jmsy=jsym(jm)
	vstagg=(vo(ip,j,ks)-vo(ip,jsy,ks)+
     %vo(ip,jm,ks)-vo(ip,jmsy,ks)+
     %vo(ip,j,kp)-vo(ip,jsy,kp)+
     %vo(ip,jm,kp)-vo(ip,jmsy,kp))*0.125
	ustagg=(uo(is,j,ks)+uo(is,j,kp))*0.5
	wstagg=(wo(ip,j,ks)+wo(ip,jsy,ks))*0.5
	vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
	vavg(iic,j,kkc)=vavg(iic,j,kkc)
     %+vstagg*dt
	uavg(iic,j,kkc)=uavg(iic,j,kkc)
     %+ustagg*dt
	wavg(iic,j,kkc)=wavg(iic,j,kkc)
     %+wstagg*dt
	vmodavg(iic,j,kkc)=vmodavg(iic,j,kkc)
     %+vmod*dt
	vorazavg(iic,j,kkc)=vorazavg(iic,j,kkc)
     %+voraz(iic,j,kkc)*dt
	vorravg(iic,j,kkc)=vorravg(iic,j,kkc)
     %+vorr(iic,j,kkc)*dt
	vorzavg(iic,j,kkc)=vorzavg(iic,j,kkc)
     %+vorz(iic,j,kkc)*dt
	if(all(flagpo(is:ip,j,ks:kp,:).le.0)) then
	   tempo1(iic,j,kkc)=tempo1(iic,j,kkc)+dt
	   pstagg=(p(ip,j,ks)+p(ip,j,kp)+
     %p(ip,jsy,ks)+p(ip,jsy,kp))*0.25
	   pavg(iic,j,kkc)=pavg(iic,j,kkc)
     %+pstagg*dt    !*ros
	   ptotavg(iic,j,kkc)=ptotavg(iic,j,kkc)
     %+(pstagg+0.5*(vmod**2.))*dt   !*ros
	endif
	tvstagg=(tv(ip,j,ks)+tv(ip,j,kp)+
     %tv(ip,jsy,ks)+tv(ip,jsy,kp))*0.25
	tvavg(iic,j,kkc)=tvavg(iic,j,kkc)
     %+tvstagg*dt
 1011	continue
 1010	continue
      else
	do 1012 kkc=1,nsez4max
	ks=ksez4(kkc)
	kp=ks+1
	do 1013 j=jy1,jy2
	jm=jmv(j)
	vstagg=(vo(is,j,ks)+vo(ip,j,ks)+
     %vo(is,jm,ks)+vo(ip,jm,ks)+
     %vo(is,j,kp)+vo(ip,j,kp)+
     %vo(is,jm,kp)+vo(ip,jm,kp))*0.125
	ustagg=(uo(is,j,ks)+uo(is,j,kp))*0.5
	wstagg=(wo(is,j,ks)+wo(ip,j,ks))*0.5
	vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
	vavg(iic,j,kkc)=vavg(iic,j,kkc)
     %+vstagg*dt
	uavg(iic,j,kkc)=uavg(iic,j,kkc)
     %+ustagg*dt
	wavg(iic,j,kkc)=wavg(iic,j,kkc)
     %+wstagg*dt
	vmodavg(iic,j,kkc)=vmodavg(iic,j,kkc)
     %+vmod*dt
	vorazavg(iic,j,kkc)=vorazavg(iic,j,kkc)
     %+voraz(iic,j,kkc)*dt
	vorravg(iic,j,kkc)=vorravg(iic,j,kkc)
     %+vorr(iic,j,kkc)*dt
	vorzavg(iic,j,kkc)=vorzavg(iic,j,kkc)
     %+vorz(iic,j,kkc)*dt
	if(all(flagpo(is:ip,j,ks:kp,:).le.0)) then
	   tempo1(iic,j,kkc)=tempo1(iic,j,kkc)+dt
	   pstagg=(p(is,j,ks)+p(ip,j,ks)+
     %p(is,j,kp)+p(ip,j,kp))*0.25
	   pavg(iic,j,kkc)=pavg(iic,j,kkc)
     %+pstagg*dt    !*ros
	   ptotavg(iic,j,kkc)=ptotavg(iic,j,kkc)
     %+(pstagg+0.5*(vmod**2.))*dt   !*ros
	endif
	tvstagg=(tv(is,j,ks)+tv(ip,j,ks)+
     %tv(is,j,kp)+tv(ip,j,kp))*0.25
	tvavg(iic,j,kkc)=tvavg(iic,j,kkc)
     %+tvstagg*dt
 1013	continue
 1012	continue
      endif
  101 continue
!
!     Time averages at some (y,z) positions along x
!
      do 104 jjc=1,nsez5
      js=jsez5(jjc)
      jp=jpv(js)
      do 104 kkc=1,nsez6max
      ks=ksez6(kkc)
      kp=ks+1
      do 103 i=ix1,ix2
      im=i-1
      vstagg=(vo(i,js,ks)+vo(i,js,kp))*0.5
      ustagg=(uo(i,js,ks)+uo(im,js,ks)+
     %uo(i,jp,ks)+uo(im,jp,ks)+
     %uo(i,js,kp)+uo(im,js,kp)+
     %uo(i,jp,kp)+uo(im,jp,kp))*0.125
      wstagg=(wo(i,js,ks)+wo(i,jp,ks))*0.5
      rvavg(i,jjc,kkc)=rvavg(i,jjc,kkc)
     %+vstagg*dt
      ruavg(i,jjc,kkc)=ruavg(i,jjc,kkc)
     %+ustagg*dt
      rwavg(i,jjc,kkc)=rwavg(i,jjc,kkc)
     %+wstagg*dt
      vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
      rvmodavg(i,jjc,kkc)=rvmodavg(i,jjc,kkc)
     %+vmod*dt
      rvorazavg(i,jjc,kkc)=rvorazavg(i,jjc,kkc)
     %+rvoraz(i,jjc,kkc)*dt
      rvorravg(i,jjc,kkc)=rvorravg(i,jjc,kkc)
     %+rvorr(i,jjc,kkc)*dt
      rvorzavg(i,jjc,kkc)=rvorzavg(i,jjc,kkc)
     %+rvorz(i,jjc,kkc)*dt
      if(all(flagpo(i,js:jp,ks:kp,:).le.0)) then
	 tempo2(i,jjc,kkc)=tempo2(i,jjc,kkc)+dt
	 pstagg=(p(i,js,ks)+p(i,jp,ks)+
     %p(i,js,kp)+p(i,jp,kp))*0.25
	 rpavg(i,jjc,kkc)=rpavg(i,jjc,kkc)+
     %pstagg*dt	  !*ros
	 rptotavg(i,jjc,kkc)=rptotavg(i,jjc,kkc)
     %+(pstagg+0.5*(vmod**2.))*dt   !*ros
      endif
      tvstagg=(tv(i,js,ks)+tv(i,jp,ks)+
     %tv(i,js,kp)+tv(i,jp,kp))*0.25
      rtvavg(i,jjc,kkc)=rtvavg(i,jjc,kkc)+
     %tvstagg*dt
  103 continue
  104 continue
!
!     Time averages at some (x,y) positions along z
!
      do 106 jjc=1,nsez91
      js=jsez91(jjc)
      jp=jpv(js)
      jsy=jsym(js)
      jpsy=jsym(jp)
      do 106 iic=1,nsez9
      is=isez9(iic)
      ip=is+1
      if(is.eq.1) then
	do 1060 k=kz1,kz2
	km=k-1
	vstagg=(vo(ip,js,k)-vo(ip,jsy,k))*0.5
	ustagg=(uo(is,js,k)+uo(is,jp,k))*0.5
	wstagg=(wo(ip,js,k)+wo(ip,jsy,k)+
     %wo(ip,jp,k)+wo(ip,jpsy,k)+wo(ip,js,km)+wo(ip,jsy,km)+
     %wo(ip,jp,km)+wo(ip,jpsy,km))*0.125
	zvavg(iic,jjc,k)=zvavg(iic,jjc,k)
     %+vstagg*dt
	zuavg(iic,jjc,k)=zuavg(iic,jjc,k)
     %+ustagg*dt
	zwavg(iic,jjc,k)=zwavg(iic,jjc,k)
     %+wstagg*dt
	vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
	zvmodavg(iic,jjc,k)=
     %zvmodavg(iic,jjc,k)+vmod*dt
	zvorazavg(iic,jjc,k)=
     %zvorazavg(iic,jjc,k)+voraz1(iic,jjc,k)*dt
	zvorravg(iic,jjc,k)=
     %zvorravg(iic,jjc,k)+vorr1(iic,jjc,k)*dt
	zvorzavg(iic,jjc,k)=
     %zvorzavg(iic,jjc,k)+vorz1(iic,jjc,k)*dt
	if(all(flagpo(is:ip,js:jp,k,:).le.0)) then
	  ztempo(iic,jjc,k)=ztempo(iic,jjc,k)+dt
	  pstagg=(p(ip,js,k)+p(ip,jsy,k)+
     %p(ip,jp,k)+p(ip,jpsy,k))*0.25
	  zpavg(iic,jjc,k)=zpavg(iic,jjc,k)+
     %pstagg*dt	  !*ros
	  zptotavg(iic,jjc,k)=zptotavg(iic,jjc,k)+
     %(pstagg+0.5*(vmod**2.))*dt   !*ros
	endif
	tvstagg=(tv(ip,js,k)+tv(ip,jsy,k)+
     %tv(ip,jp,k)+tv(ip,jpsy,k))*0.25
	ztvavg(iic,jjc,k)=ztvavg(iic,jjc,k)+
     %tvstagg*dt
 1060	continue
      else
	do 1061 k=kz1,kz2
	km=k-1
	vstagg=(vo(is,js,k)+vo(ip,js,k))*0.5
	ustagg=(uo(is,js,k)+uo(is,jp,k))*0.5
	wstagg=(wo(is,js,k)+wo(ip,js,k)+
     %wo(is,jp,k)+wo(ip,jp,k)+wo(is,js,km)+wo(ip,js,km)+
     %wo(is,jp,km)+wo(ip,jp,km))*0.125
	zvavg(iic,jjc,k)=zvavg(iic,jjc,k)
     %+vstagg*dt
	zuavg(iic,jjc,k)=zuavg(iic,jjc,k)
     %+ustagg*dt
	zwavg(iic,jjc,k)=zwavg(iic,jjc,k)
     %+wstagg*dt
	vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
	zvmodavg(iic,jjc,k)=
     %zvmodavg(iic,jjc,k)+vmod*dt
	zvorazavg(iic,jjc,k)=
     %zvorazavg(iic,jjc,k)+voraz1(iic,jjc,k)*dt
	zvorravg(iic,jjc,k)=
     %zvorravg(iic,jjc,k)+vorr1(iic,jjc,k)*dt
	zvorzavg(iic,jjc,k)=
     %zvorzavg(iic,jjc,k)+vorz1(iic,jjc,k)*dt
	if(all(flagpo(is:ip,js:jp,k,:).le.0)) then
	  ztempo(iic,jjc,k)=ztempo(iic,jjc,k)+dt
	  pstagg=(p(is,js,k)+p(ip,js,k)+
     %p(is,jp,k)+p(ip,jp,k))*0.25
	  zpavg(iic,jjc,k)=zpavg(iic,jjc,k)+
     %pstagg*dt	  !*ros
	  zptotavg(iic,jjc,k)=zptotavg(iic,jjc,k)+
     %(pstagg+0.5*(vmod**2.))*dt   !*ros
	endif
	tvstagg=(tv(is,js,k)+tv(ip,js,k)+
     %tv(is,jp,k)+tv(ip,jp,k))*0.25
	ztvavg(iic,jjc,k)=ztvavg(iic,jjc,k)+
     %tvstagg*dt
 1061	continue
      endif
  106 continue
      if(tempo.ge.periodo4) then
	 iplant8=iplant8+1
	 do 200 kkc=1,nsez4max
	 do 200 iic=1,nsez3
	 do 200 j=jy1,jy2
	 vmodavg(iic,j,kkc)=vmodavg(iic,j,kkc)/
     %tempo
	 vavg(iic,j,kkc)=vavg(iic,j,kkc)/
     %tempo
	 vorazavg(iic,j,kkc)=vorazavg(iic,j,kkc)/
     %tempo
	 uavg(iic,j,kkc)=uavg(iic,j,kkc)/
     %tempo
	 vorravg(iic,j,kkc)=vorravg(iic,j,kkc)/
     %tempo
	 wavg(iic,j,kkc)=wavg(iic,j,kkc)/
     %tempo
	 vorzavg(iic,j,kkc)=vorzavg(iic,j,kkc)/
     %tempo
	 if(tempo1(iic,j,kkc).ne.0.) then
	    pavg(iic,j,kkc)=pavg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	    ptotavg(iic,j,kkc)=ptotavg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	 endif
	 tvavg(iic,j,kkc)=tvavg(iic,j,kkc)/
     %tempo
	 tvavg(iic,j,kkc)=tvavg(iic,j,kkc)/ru1
  200	 continue
	 if(myrank.eq.0) then
	    mm1=121
	    mm2=122
	    do 201 kkc=1,nsez4max
	    do 201 iic=1,nsez3
	    write(mm1,333)'i=',isez3(iic),'____k=',ksez4(kkc),time-tempo*0.5
	    write(mm2,333)'i=',isez3(iic),'____k=',ksez4(kkc),time-tempo*0.5
	    do 201 j=jy1,jy2
	    write(mm1,10)
     %iplant8,isez3(iic),j,ksez4(kkc),vmodavg(iic,j,kkc),
     %uavg(iic,j,kkc),vavg(iic,j,kkc),wavg(iic,j,kkc),
     %pavg(iic,j,kkc),ptotavg(iic,j,kkc),tvavg(iic,j,kkc)
	    write(mm2,10)
     %iplant8,isez3(iic),j,ksez4(kkc),vorravg(iic,j,kkc),
     %vorazavg(iic,j,kkc),vorzavg(iic,j,kkc)
  201	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(nsez4maxr,1,MPI_INTEGER,jp,1,
     %MPI_COMM_EDDY,STATUS,IERR)
		  ALLOCATE(ksez4r(nsez4maxr))
		  ALLOCATE(vmodavgr(nsez3,ny,nsez4maxr),
     %uavgr(nsez3,ny,nsez4maxr),vavgr(nsez3,ny,nsez4maxr),
     %wavgr(nsez3,ny,nsez4maxr),pavgr(nsez3,ny,nsez4maxr),
     %ptotavgr(nsez3,ny,nsez4maxr),vorravgr(nsez3,ny,nsez4maxr),
     %vorazavgr(nsez3,ny,nsez4maxr),vorzavgr(nsez3,ny,nsez4maxr),
     %tvavgr(nsez3,ny,nsez4maxr))
		  CALL MPI_RECV(ksez4r(1:nsez4maxr),nsez4maxr,MPI_INTEGER,jp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vmodavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(uavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(wavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,6,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(pavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,7,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(ptotavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vorravgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vorazavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vorzavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(tvavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 203 kkc=1,nsez4maxr
		  k=ksez4r(kkc)
		  kg=jp*(nz-2)+k
		  do 203 iic=1,nsez3
		  write(mm1,333)'i=',isez3(iic),'____k=',kg,time-tempo*0.5
		  write(mm2,333)'i=',isez3(iic),'____k=',kg,time-tempo*0.5
		  do 203 j=jy1,jy2
		  write(mm1,10)
     %iplant8,isez3(iic),j,kg,vmodavgr(iic,j,kkc),
     %uavgr(iic,j,kkc),vavgr(iic,j,kkc),wavgr(iic,j,kkc),
     %pavgr(iic,j,kkc),ptotavgr(iic,j,kkc),tvavgr(iic,j,kkc)
		  write(mm2,10)
     %iplant8,isez3(iic),j,kg,vorravgr(iic,j,kkc),
     %vorazavgr(iic,j,kkc),vorzavgr(iic,j,kkc)
  203		  continue
		  DEALLOCATE(ksez4r)
		  DEALLOCATE(vmodavgr,uavgr,vavgr,wavgr,pavgr,
     %ptotavgr,vorravgr,vorazavgr,vorzavgr,tvavgr)
	       enddo
	    endif
	 else
	    CALL MPI_SEND(nsez4max,1,MPI_INTEGER,0,1,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ksez4(1:nsez4max),nsez4max,MPI_INTEGER,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vmodavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(uavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(wavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(pavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ptotavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vorravg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vorazavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vorzavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(tvavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
c
	 do 204 kkc=1,nsez6max
	 do 204 jjc=1,nsez5
	 do 204 i=ix1,ix2
	 rvmodavg(i,jjc,kkc)=rvmodavg(i,jjc,kkc)/
     %tempo
	 rvavg(i,jjc,kkc)=rvavg(i,jjc,kkc)/
     %tempo
	 rvorazavg(i,jjc,kkc)=rvorazavg(i,jjc,kkc)/
     %tempo
	 ruavg(i,jjc,kkc)=ruavg(i,jjc,kkc)/
     %tempo
	 rvorravg(i,jjc,kkc)=rvorravg(i,jjc,kkc)/
     %tempo
	 rwavg(i,jjc,kkc)=rwavg(i,jjc,kkc)/
     %tempo
	 rvorzavg(i,jjc,kkc)=rvorzavg(i,jjc,kkc)/
     %tempo
	 if(tempo2(i,jjc,kkc).ne.0.) then
	    rpavg(i,jjc,kkc)=rpavg(i,jjc,kkc)/
     %tempo2(i,jjc,kkc)
	    rptotavg(i,jjc,kkc)=rptotavg(i,jjc,kkc)/
     %tempo2(i,jjc,kkc)
	 endif
	 rtvavg(i,jjc,kkc)=rtvavg(i,jjc,kkc)/
     %tempo
	 rtvavg(i,jjc,kkc)=rtvavg(i,jjc,kkc)/ru1
 204	 continue
	 if(myrank.eq.0) then
	    mm1=123
	    mm2=124
	    do 205 kkc=1,nsez6max
	    do 205 jjc=1,nsez5
	    write(mm1,333)'j=',jsez5(jjc),'____k=',ksez6(kkc),time-tempo*0.5
	    write(mm2,333)'j=',jsez5(jjc),'____k=',ksez6(kkc),time-tempo*0.5
	    do 205 i=ix1,ix2
	    write(mm1,10)
     %iplant8,i,jsez5(jjc),ksez6(kkc),rvmodavg(i,jjc,kkc),
     %ruavg(i,jjc,kkc),rvavg(i,jjc,kkc),rwavg(i,jjc,kkc),
     %rpavg(i,jjc,kkc),rptotavg(i,jjc,kkc),rtvavg(i,jjc,kkc)
	    write(mm2,10)
     %iplant8,i,jsez5(jjc),ksez6(kkc),rvorravg(i,jjc,kkc),
     %rvorazavg(i,jjc,kkc),rvorzavg(i,jjc,kkc)
 205	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(nsez6maxr,1,MPI_INTEGER,jp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
		  ALLOCATE(ksez6r(nsez6maxr))
		  ALLOCATE(rvmodavgr(nx,nsez5,nsez6maxr),
     %ruavgr(nx,nsez5,nsez6maxr),rvavgr(nx,nsez5,nsez6maxr),
     %rwavgr(nx,nsez5,nsez6maxr),rpavgr(nx,nsez5,nsez6maxr),
     %rptotavgr(nx,nsez5,nsez6maxr),rvorravgr(nx,nsez5,nsez6maxr),
     %rvorazavgr(nx,nsez5,nsez6maxr),rvorzavgr(nx,nsez5,nsez6maxr),
     %rtvavgr(nx,nsez5,nsez6maxr))
		  CALL MPI_RECV(ksez6r(1:nsez6maxr),nsez6maxr,MPI_INTEGER,jp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvmodavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(ruavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,16,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,17,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rwavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,18,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rpavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,19,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rptotavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,20,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvorravgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,21,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvorazavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,22,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvorzavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,23,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rtvavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,24,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 207 kkc=1,nsez6maxr
		  k=ksez6r(kkc)
		  kg=jp*(nz-2)+k
		  do 207 jjc=1,nsez5
		  write(mm1,333)'j=',jsez5(jjc),'____k=',kg,time-tempo*0.5
		  write(mm2,333)'j=',jsez5(jjc),'____k=',kg,time-tempo*0.5
		  do 207 i=ix1,ix2
		  write(mm1,10)
     %iplant8,i,jsez5(jjc),kg,rvmodavgr(i,jjc,kkc),
     %ruavgr(i,jjc,kkc),rvavgr(i,jjc,kkc),rwavgr(i,jjc,kkc),
     %rpavgr(i,jjc,kkc),rptotavgr(i,jjc,kkc),rtvavgr(i,jjc,kkc)
		  write(mm2,10)
     %iplant8,i,jsez5(jjc),kg,rvorravgr(i,jjc,kkc),
     %rvorazavgr(i,jjc,kkc),rvorzavgr(i,jjc,kkc)
 207		  continue
		  DEALLOCATE(ksez6r)
		  DEALLOCATE(rvmodavgr,ruavgr,rvavgr,rwavgr,rpavgr,
     %rptotavgr,rvorravgr,rvorazavgr,rvorzavgr,rtvavgr)
	       enddo
	    endif
	 else
	    CALL MPI_SEND(nsez6max,1,MPI_INTEGER,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ksez6(1:nsez6max),nsez6max,MPI_INTEGER,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvmodavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ruavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,16,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,17,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rwavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,18,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rpavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,19,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rptotavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,20,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvorravg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,21,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvorazavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,22,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvorzavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,23,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rtvavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,24,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
c
	 do 208 jjc=1,nsez91
	 do 208 iic=1,nsez9
	 do 208 k=kz1,kz2
	 zvmodavg(iic,jjc,k)=zvmodavg(iic,jjc,k)/tempo
	 zvavg(iic,jjc,k)=zvavg(iic,jjc,k)/tempo
	 zvorazavg(iic,jjc,k)=zvorazavg(iic,jjc,k)/tempo
	 zuavg(iic,jjc,k)=zuavg(iic,jjc,k)/tempo
	 zvorravg(iic,jjc,k)=zvorravg(iic,jjc,k)/tempo
	 zwavg(iic,jjc,k)=zwavg(iic,jjc,k)/tempo
	 zvorzavg(iic,jjc,k)=zvorzavg(iic,jjc,k)/tempo
	 if(ztempo(iic,jjc,k).ne.0.) then
	    zpavg(iic,jjc,k)=
     %zpavg(iic,jjc,k)/ztempo(iic,jjc,k)
	    zptotavg(iic,jjc,k)=
     %zptotavg(iic,jjc,k)/ztempo(iic,jjc,k)
	 endif
	 ztvavg(iic,jjc,k)=ztvavg(iic,jjc,k)/tempo
	 ztvavg(iic,jjc,k)=ztvavg(iic,jjc,k)/ru1
 208	 continue
	 mm1=125
	 mm2=126
	 do 209 jjc=1,nsez91
	 do 209 iic=1,nsez9
	 if(myrank.eq.0) then
	    write(mm1,333)'i=',isez9(iic),'____j=',jsez91(jjc),time-tempo*0.5
	    write(mm2,333)'i=',isez9(iic),'____j=',jsez91(jjc),time-tempo*0.5
	    do 210 k=kz1,kz2
	    write(mm1,10)iplant8,isez9(iic),jsez91(jjc),k,
     %zvmodavg(iic,jjc,k),zuavg(iic,jjc,k),zvavg(iic,jjc,k),
     %zwavg(iic,jjc,k),zpavg(iic,jjc,k),zptotavg(iic,jjc,k),
     %ztvavg(iic,jjc,k)
	    write(mm2,10)iplant8,isez9(iic),jsez91(jjc),k,
     %zvorravg(iic,jjc,k),zvorazavg(iic,jjc,k),zvorzavg(iic,jjc,k)
 210	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(zvmodavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,25,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zuavgr(kz1:kz2),(nz-2),MTYPE,jp,26,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvavgr(kz1:kz2),(nz-2),MTYPE,jp,27,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zwavgr(kz1:kz2),(nz-2),MTYPE,jp,28,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zpavgr(kz1:kz2),(nz-2),MTYPE,jp,29,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zptotavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,30,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvorravg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,31,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvorazavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,32,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvorzavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,33,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(ztvavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,34,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 211 k=kz1,kz2
		  kg=jp*(nz-2)+k
		  write(mm1,10)iplant8,isez9(iic),jsez91(jjc),kg,
     %zvmodavg(iic,jjc,k),zuavgr(k),zvavgr(k),zwavgr(k),
     %zpavgr(k),zptotavg(iic,jjc,k),ztvavg(iic,jjc,k)
		  write(mm2,10)iplant8,isez9(iic),jsez91(jjc),kg,
     %zvorravg(iic,jjc,k),zvorazavg(iic,jjc,k),zvorzavg(iic,jjc,k)
 211		  continue
	       enddo
	    endif
	 else
	    CALL MPI_SEND(zvmodavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,25,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zuavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,26,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,27,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zwavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,28,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zpavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,29,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zptotavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,30,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvorravg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,31,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvorazavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,32,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvorzavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,33,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ztvavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,34,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
 209	 continue
      endif
      if(irms.eq.1)call rmsmodstaggnoavg(nx,ny,nz,uo,vo,wo,p,dt,nbd,flagpo,time)
      if(tempo.ge.periodo4) then
	 do 110 kkc=1,nsez4max
	 do 110 iic=1,nsez3
	 do 110 j=jy1,jy2
	 vmodavg(iic,j,kkc)=0.
	 vavg(iic,j,kkc)=0.
	 vorazavg(iic,j,kkc)=0.
	 uavg(iic,j,kkc)=0.
	 vorravg(iic,j,kkc)=0.
	 wavg(iic,j,kkc)=0.
	 vorzavg(iic,j,kkc)=0.
	 pavg(iic,j,kkc)=0.
	 ptotavg(iic,j,kkc)=0.
	 tvavg(iic,j,kkc)=0.
	 tempo1(iic,j,kkc)=0.
  110	 continue
	 do 112 kkc=1,nsez6max
	 do 112 jjc=1,nsez5
	 do 112 i=ix1,ix2
	 rvmodavg(i,jjc,kkc)=0.
	 rvavg(i,jjc,kkc)=0.
	 rvorazavg(i,jjc,kkc)=0.
	 ruavg(i,jjc,kkc)=0.
	 rvorravg(i,jjc,kkc)=0.
	 rwavg(i,jjc,kkc)=0.
	 rvorzavg(i,jjc,kkc)=0.
	 rpavg(i,jjc,kkc)=0.
	 rptotavg(i,jjc,kkc)=0.
	 rtvavg(i,jjc,kkc)=0.
	 tempo2(i,jjc,kkc)=0.
  112	 continue
	 do 114 jjc=1,nsez91
	 do 114 iic=1,nsez9
	 do 114 k=kz1,kz2
	 zvmodavg(iic,jjc,k)=0.
	 zvavg(iic,jjc,k)=0.
	 zvorazavg(iic,jjc,k)=0.
	 zuavg(iic,jjc,k)=0.
	 zvorravg(iic,jjc,k)=0.
	 zwavg(iic,jjc,k)=0.
	 zvorzavg(iic,jjc,k)=0.
	 zpavg(iic,jjc,k)=0.
	 zptotavg(iic,jjc,k)=0.
	 ztvavg(iic,jjc,k)=0.
	 ztempo(iic,jjc,k)=0.
  114	 continue
 333	 format(a2,i6,a6,i6,f20.10)
	 tempo=0.
      endif
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine rmsmodstaggnoavg(nx,ny,nz,uo,vo,wo,p,dt,nbd,flagpo,time)
c
      use sections
      use mediey
      use mediex
      use mediez
      use rmsy
      use rmsx
      use rmsz
      include'common.h'
      include'averages.h'
      include'mpif.h'
c
c Global variables
      integer nx,ny,nz,nbd,flagpo(nx,ny,nz,nbd)
      real dt,time
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
c
c Local variables
      integer iic,i,is,ip,im,jjc,j,js,jp,jm,jsy,jpsy,jmsy,
     %kkc,k,ks,kp,km,mm1,mm2,kg,nsez4maxr,nsez6maxr
      integer, dimension(:), allocatable :: ksez4r,ksez6r
      integer STATUS(MPI_STATUS_SIZE)
      real ustagg,vstagg,wstagg,pstagg
      real, dimension(:,:,:), allocatable :: urmsr,vrmsr,wrmsr,
     %prmsr,uvsr,uwsr,vwsr,rurmsr,rvrmsr,rwrmsr,rprmsr,ruvsr,
     %ruwsr,rvwsr
c
  10  format(4i6,10f25.8)
c
      do 101 iic=1,nsez3
      is=isez3(iic)
      ip=is+1
      if(is.eq.1) then
	do 1010 kkc=1,nsez4max
	ks=ksez4(kkc)
	kp=ks+1
	do 1011 j=jy1,jy2
	jm=jmv(j)
	jsy=jsym(j)
	jmsy=jsym(jm)
	vstagg=0.125*(vo(ip,j,ks)-vo(ip,jsy,ks)+
     %vo(ip,jm,ks)-vo(ip,jmsy,ks)+
     %vo(ip,j,kp)-vo(ip,jsy,kp)+
     %vo(ip,jm,kp)-vo(ip,jmsy,kp))
	ustagg=0.5*(uo(is,j,ks)+uo(is,j,kp))
	wstagg=0.5*(wo(ip,j,ks)+wo(ip,jsy,ks))
	vrms(iic,j,kkc)=
     %vrms(iic,j,kkc)+dt*(vstagg**2.)
	urms(iic,j,kkc)=
     %urms(iic,j,kkc)+dt*(ustagg**2.)
	wrms(iic,j,kkc)=
     %wrms(iic,j,kkc)+dt*(wstagg**2.)
	uvs(iic,j,kkc)=
     %uvs(iic,j,kkc)+dt*ustagg*vstagg
	uws(iic,j,kkc)=
     %uws(iic,j,kkc)+dt*ustagg*wstagg
	vws(iic,j,kkc)=
     %vws(iic,j,kkc)+dt*vstagg*wstagg
	if(all(flagpo(is:ip,j,ks:kp,:).le.0)) then
	   pstagg=0.25*(p(ip,j,ks)+p(ip,jsy,ks)+
     %p(ip,j,kp)+p(ip,jsy,kp))
	   prms(iic,j,kkc)=
     %prms(iic,j,kkc)+dt*(pstagg**2.)	!*ros**2.
	endif
 1011	continue
 1010	continue
      else
	do 1012 kkc=1,nsez4max
	ks=ksez4(kkc)
	kp=ks+1
	do 1013 j=jy1,jy2
	jm=jmv(j)
	vstagg=0.125*(vo(is,j,ks)+vo(ip,j,ks)+
     %vo(is,jm,ks)+vo(ip,jm,ks)+
     %vo(is,j,kp)+vo(ip,j,kp)+
     %vo(is,jm,kp)+vo(ip,jm,kp))
	ustagg=0.5*(uo(is,j,ks)+uo(is,j,kp))
	wstagg=0.5*(wo(is,j,ks)+wo(ip,j,ks))
	vrms(iic,j,kkc)=
     %vrms(iic,j,kkc)+dt*(vstagg**2.)
	urms(iic,j,kkc)=
     %urms(iic,j,kkc)+dt*(ustagg**2.)
	wrms(iic,j,kkc)=
     %wrms(iic,j,kkc)+dt*(wstagg**2.)
	uvs(iic,j,kkc)=
     %uvs(iic,j,kkc)+dt*ustagg*vstagg
	uws(iic,j,kkc)=
     %uws(iic,j,kkc)+dt*ustagg*wstagg
	vws(iic,j,kkc)=
     %vws(iic,j,kkc)+dt*vstagg*wstagg
	if(all(flagpo(is:ip,j,ks:kp,:).le.0)) then
	   pstagg=0.25*(p(is,j,ks)+p(ip,j,ks)+
     %p(is,j,kp)+p(ip,j,kp))
	   prms(iic,j,kkc)=
     %prms(iic,j,kkc)+dt*(pstagg**2.)	!*ros**2.
	endif
 1013	continue
 1012	continue
      endif
  101 continue
c
      do 104 jjc=1,nsez5
      js=jsez5(jjc)
      jp=jpv(js)
      do 104 kkc=1,nsez6max
      ks=ksez6(kkc)
      kp=ks+1
      do 103 i=ix1,ix2
      im=i-1
      vstagg=0.5*(vo(i,js,ks)+vo(i,js,kp))
      ustagg=0.125*(uo(i,js,ks)+uo(im,js,ks)+
     %uo(i,jp,ks)+uo(im,jp,ks)+
     %uo(i,js,kp)+uo(im,js,kp)+
     %uo(i,jp,kp)+uo(im,jp,kp))
      wstagg=0.5*(wo(i,js,ks)+wo(i,jp,ks))
      rvrms(i,jjc,kkc)=
     %rvrms(i,jjc,kkc)+dt*(vstagg**2.)
      rurms(i,jjc,kkc)=
     %rurms(i,jjc,kkc)+dt*(ustagg**2.)
      rwrms(i,jjc,kkc)=
     %rwrms(i,jjc,kkc)+dt*(wstagg**2.)
      ruvs(i,jjc,kkc)=
     %ruvs(i,jjc,kkc)+dt*ustagg*vstagg
      ruws(i,jjc,kkc)=
     %ruws(i,jjc,kkc)+dt*ustagg*wstagg
      rvws(i,jjc,kkc)=
     %rvws(i,jjc,kkc)+dt*vstagg*wstagg
      if(all(flagpo(i,js:jp,ks:kp,:).le.0)) then
	 pstagg=0.25*(p(i,js,ks)+p(i,jp,ks)+
     %p(i,js,kp)+p(i,jp,kp))
	 rprms(i,jjc,kkc)=
     %rprms(i,jjc,kkc)+dt*(pstagg**2.)	 !*ros**2.
      endif
  103 continue
  104 continue
c
      do 106 jjc=1,nsez91
      js=jsez91(jjc)
      jp=jpv(js)
      jsy=jsym(js)
      jpsy=jsym(jp)
      do 106 iic=1,nsez9
      is=isez9(iic)
      ip=is+1
      if(is.eq.1) then
	do 1060 k=kz1,kz2
	km=k-1
	vstagg=0.5*(vo(ip,js,k)-vo(ip,jsy,k))
	ustagg=0.5*(uo(is,js,k)+uo(is,jp,k))
	wstagg=0.125*(wo(ip,js,k)+wo(ip,jsy,k)+
     %wo(ip,jp,k)+wo(ip,jpsy,k)+
     %wo(ip,js,km)+wo(ip,jsy,km)+
     %wo(ip,jp,km)+wo(ip,jpsy,km))
	zvrms(iic,jjc,k)=zvrms(iic,jjc,k)
     %+dt*(vstagg**2.)
	zurms(iic,jjc,k)=zurms(iic,jjc,k)
     %+dt*(ustagg**2.)
	zwrms(iic,jjc,k)=zwrms(iic,jjc,k)
     %+dt*(wstagg**2.)
	zuvs(iic,jjc,k)=zuvs(iic,jjc,k)
     %+dt*ustagg*vstagg
	zuws(iic,jjc,k)=zuws(iic,jjc,k)
     %+dt*ustagg*wstagg
	zvws(iic,jjc,k)=zvws(iic,jjc,k)
     %+dt*vstagg*wstagg
	if(all(flagpo(is:ip,js:jp,k,:).le.0)) then
	   pstagg=0.25*(p(ip,js,k)+p(ip,jsy,k)+
     %p(ip,jp,k)+p(ip,jpsy,k))
	   zprms(iic,jjc,k)=zprms(iic,jjc,k)
     %+dt*(pstagg**2.)	 !*ros**2.
	endif
 1061	continue
 1060	continue
      else
	do 1062 k=kz1,kz2
	km=k-1
	vstagg=0.5*(vo(is,js,k)+vo(ip,js,k))
	ustagg=0.5*(uo(is,js,k)+uo(is,jp,k))
	wstagg=0.125*(wo(is,js,k)+wo(ip,js,k)+
     %wo(is,jp,k)+wo(ip,jp,k)+
     %wo(is,js,km)+wo(ip,js,km)+
     %wo(is,jp,km)+wo(ip,jp,km))
	zvrms(iic,jjc,k)=zvrms(iic,jjc,k)
     %+dt*(vstagg**2.)
	zurms(iic,jjc,k)=zurms(iic,jjc,k)
     %+dt*(ustagg**2.)
	zwrms(iic,jjc,k)=zwrms(iic,jjc,k)
     %+dt*(wstagg**2.)
	zuvs(iic,jjc,k)=zuvs(iic,jjc,k)
     %+dt*ustagg*vstagg
	zuws(iic,jjc,k)=zuws(iic,jjc,k)
     %+dt*ustagg*wstagg
	zvws(iic,jjc,k)=zvws(iic,jjc,k)
     %+dt*vstagg*wstagg
	if(all(flagpo(is:ip,js:jp,k,:).le.0)) then
	   pstagg=0.25*(p(is,js,k)+p(ip,js,k)+
     %p(is,jp,k)+p(ip,jp,k))
	   zprms(iic,jjc,k)=zprms(iic,jjc,k)
     %+dt*(pstagg**2.)	  !*ros**2.
	endif
 1063	continue
 1062	continue
      endif
  106 continue
c
      if(tempo.ge.periodo4) then
	 do 121 kkc=1,nsez4max
	 do 121 iic=1,nsez3
	 do 120 j=jy1,jy2
	 vrms(iic,j,kkc)=
     %sqrt(abs((vrms(iic,j,kkc)/tempo)
     %-vavg(iic,j,kkc)**2.))
	 urms(iic,j,kkc)=
     %sqrt(abs((urms(iic,j,kkc)/tempo)
     %-uavg(iic,j,kkc)**2.))
	 wrms(iic,j,kkc)=
     %sqrt(abs((wrms(iic,j,kkc)/tempo)
     %-wavg(iic,j,kkc)**2.))
	 uvs(iic,j,kkc)=
     %uvs(iic,j,kkc)/tempo
     %-uavg(iic,j,kkc)*vavg(iic,j,kkc)
	 uws(iic,j,kkc)=
     %uws(iic,j,kkc)/tempo
     %-uavg(iic,j,kkc)*wavg(iic,j,kkc)
	 vws(iic,j,kkc)=
     %vws(iic,j,kkc)/tempo
     %-vavg(iic,j,kkc)*wavg(iic,j,kkc)
	 if(tempo1(iic,j,kkc).ne.0.)
     %prms(iic,j,kkc)=
     %sqrt(abs((prms(iic,j,kkc)/tempo1(iic,j,kkc))
     %-pavg(iic,j,kkc)**2.))
  120	 continue
  121	 continue
	 if(myrank.eq.0) then
	    mm1=127
	    mm2=128
	    do 122 kkc=1,nsez4max
	    do 122 iic=1,nsez3
	    write(mm1,333)'i=',isez3(iic),'____k=',ksez4(kkc),time-0.5*tempo
	    write(mm2,333)'i=',isez3(iic),'____k=',ksez4(kkc),time-0.5*tempo
	    do 122 j=jy1,jy2
	    write(mm1,10)iplant8,isez3(iic),j,ksez4(kkc),
     %urms(iic,j,kkc),vrms(iic,j,kkc),wrms(iic,j,kkc),prms(iic,j,kkc)
	    write(mm2,10)iplant8,isez3(iic),j,ksez4(kkc),
     %uvs(iic,j,kkc),uws(iic,j,kkc),vws(iic,j,kkc)
  122	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(nsez4maxr,1,MPI_INTEGER,jp,1,
     %MPI_COMM_EDDY,STATUS,IERR)
		  ALLOCATE(ksez4r(nsez4maxr))
		  ALLOCATE(urmsr(nsez3,ny,nsez4maxr),
     %vrmsr(nsez3,ny,nsez4maxr),wrmsr(nsez3,ny,nsez4maxr),
     %prmsr(nsez3,ny,nsez4maxr),uvsr(nsez3,ny,nsez4maxr),
     %uwsr(nsez3,ny,nsez4maxr),vwsr(nsez3,ny,nsez4maxr))
		  CALL MPI_RECV(ksez4r(1:nsez4maxr),nsez4maxr,MPI_INTEGER,jp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(urmsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vrmsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(wrmsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(prmsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,6,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(uvsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,7,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(uwsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vwsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 123 kkc=1,nsez4maxr
		  k=ksez4r(kkc)
		  kg=jp*(nz-2)+k
		  do 123 iic=1,nsez3
		  write(mm1,333)'i=',isez3(iic),'____k=',kg,time-0.5*tempo
		  write(mm2,333)'i=',isez3(iic),'____k=',kg,time-0.5*tempo
		  do 123 j=jy1,jy2
		  write(mm1,10)
     %iplant8,isez3(iic),j,kg,urmsr(iic,j,kkc),vrmsr(iic,j,kkc),
     %wrmsr(iic,j,kkc),prmsr(iic,j,kkc)
		  write(mm2,10)
     %iplant8,isez3(iic),j,kg,
     %uvsr(iic,j,kkc),uwsr(iic,j,kkc),vwsr(iic,j,kkc)
  123		  continue
		  DEALLOCATE(ksez4r)
		  DEALLOCATE(urmsr,vrmsr,wrmsr,prmsr,uvsr,uwsr,vwsr)
	       enddo
	    endif
	 else
	    CALL MPI_SEND(nsez4max,1,MPI_INTEGER,0,1,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ksez4(1:nsez4max),nsez4max,MPI_INTEGER,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(urms(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vrms(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(wrms(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(prms(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(uvs(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(uws(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vws(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
c
	 do kkc=1,nsez4max
	    do iic=1,nsez3
	       do j=jy1,jy2
		  vrms(iic,j,kkc)=0.
		  urms(iic,j,kkc)=0.
		  wrms(iic,j,kkc)=0.
		  prms(iic,j,kkc)=0.
		  uvs(iic,j,kkc)=0.
		  uws(iic,j,kkc)=0.
		  vws(iic,j,kkc)=0.
	       enddo
	    enddo
	 enddo
c
	 do 124 kkc=1,nsez6max
	 do 124 jjc=1,nsez5
	 do 124 i=ix1,ix2
	 rvrms(i,jjc,kkc)=
     %sqrt(abs((rvrms(i,jjc,kkc)/tempo)
     %-rvavg(i,jjc,kkc)**2.))
	 rurms(i,jjc,kkc)=
     %sqrt(abs((rurms(i,jjc,kkc)/tempo)
     %-ruavg(i,jjc,kkc)**2.))
	 rwrms(i,jjc,kkc)=
     %sqrt(abs((rwrms(i,jjc,kkc)/tempo)
     %-rwavg(i,jjc,kkc)**2.))
	 ruvs(i,jjc,kkc)=
     %ruvs(i,jjc,kkc)/tempo
     %-ruavg(i,jjc,kkc)*rvavg(i,jjc,kkc)
	 ruws(i,jjc,kkc)=
     %ruws(i,jjc,kkc)/tempo
     %-ruavg(i,jjc,kkc)*rwavg(i,jjc,kkc)
	 rvws(i,jjc,kkc)=
     %rvws(i,jjc,kkc)/tempo
     %-rvavg(i,jjc,kkc)*rwavg(i,jjc,kkc)
	 if(tempo2(i,jjc,kkc).ne.0.)
     %rprms(i,jjc,kkc)=
     %sqrt(abs((rprms(i,jjc,kkc)/tempo2(i,jjc,kkc))
     %-rpavg(i,jjc,kkc)**2.))
  124	 continue
c
	 if(myrank.eq.0) then
	    mm1=129
	    mm2=130
	    do 125 kkc=1,nsez6max
	    do 125 jjc=1,nsez5
	    write(mm1,333)'j=',jsez5(jjc),'____k=',ksez6(kkc),time-0.5*tempo
	    write(mm2,333)'j=',jsez5(jjc),'____k=',ksez6(kkc),time-0.5*tempo
	    do 125 i=ix1,ix2
	    write(mm1,10)iplant8,i,jsez5(jjc),ksez6(kkc),
     %rurms(i,jjc,kkc),rvrms(i,jjc,kkc),
     %rwrms(i,jjc,kkc),rprms(i,jjc,kkc)
	    write(mm2,10)iplant8,i,jsez5(jjc),ksez6(kkc),
     %ruvs(i,jjc,kkc),ruws(i,jjc,kkc),rvws(i,jjc,kkc)
  125	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(nsez6maxr,1,MPI_INTEGER,jp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  ALLOCATE(ksez6r(nsez6maxr))
		  ALLOCATE(rurmsr(nx,nsez5,nsez6maxr),
     %rvrmsr(nx,nsez5,nsez6maxr),rwrmsr(nx,nsez5,nsez6maxr),
     %rprmsr(nx,nsez5,nsez6maxr),ruvsr(nx,nsez5,nsez6maxr),
     %ruwsr(nx,nsez5,nsez6maxr),rvwsr(nx,nsez5,nsez6maxr))
		  CALL MPI_RECV(ksez6r(1:nsez6maxr),nsez6maxr,MPI_INTEGER,jp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rurmsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvrmsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rwrmsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rprmsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(ruvsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,16,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(ruwsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,17,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvwsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,18,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 126 kkc=1,nsez6maxr
		  k=ksez6r(kkc)
		  kg=jp*(nz-2)+k
		  do 126 jjc=1,nsez5
		  write(mm1,333)'j=',jsez5(jjc),'____k=',kg,time-0.5*tempo
		  write(mm2,333)'j=',jsez5(jjc),'____k=',kg,time-0.5*tempo
		  do 126 i=ix1,ix2
		  write(mm1,10)
     %iplant8,i,jsez5(jjc),kg,
     %rurmsr(i,jjc,kkc),rvrmsr(i,jjc,kkc),
     %rwrmsr(i,jjc,kkc),rprmsr(i,jjc,kkc)
		  write(mm2,10)
     %iplant8,i,jsez5(jjc),kg,
     %ruvsr(i,jjc,kkc),ruwsr(i,jjc,kkc),rvwsr(i,jjc,kkc)
  126		  continue
		  DEALLOCATE(ksez6r)
		  DEALLOCATE(rurmsr,rvrmsr,rwrmsr,rprmsr,
     %ruvsr,ruwsr,rvwsr)
	       enddo
	    endif
	 else
	    CALL MPI_SEND(nsez6max,1,MPI_INTEGER,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ksez6(1:nsez6max),nsez6max,MPI_INTEGER,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rurms(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvrms(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rwrms(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rprms(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ruvs(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,16,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ruws(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,17,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvws(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,18,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
c
	 do kkc=1,nsez6max
	    do jjc=1,nsez5
	       do i=ix1,ix2
		  rvrms(i,jjc,kkc)=0.
		  rurms(i,jjc,kkc)=0.
		  rwrms(i,jjc,kkc)=0.
		  rprms(i,jjc,kkc)=0.
		  ruvs(i,jjc,kkc)=0.
		  ruws(i,jjc,kkc)=0.
		  rvws(i,jjc,kkc)=0.
	       enddo
	    enddo
	 enddo
c
	 do 127 jjc=1,nsez91
	 do 127 iic=1,nsez9
	 do 127 k=kz1,kz2
	 zvrms(iic,jjc,k)=
     %sqrt(abs(zvrms(iic,jjc,k)/tempo
     %-zvavg(iic,jjc,k)**2.))
	 zurms(iic,jjc,k)=
     %sqrt(abs(zurms(iic,jjc,k)/tempo
     %-zuavg(iic,jjc,k)**2.))
	 zwrms(iic,jjc,k)=
     %sqrt(abs(zwrms(iic,jjc,k)/tempo
     %-zwavg(iic,jjc,k)**2.))
	 zuvs(iic,jjc,k)=
     %zuvs(iic,jjc,k)/tempo
     %-zuavg(iic,jjc,k)*zvavg(iic,jjc,k)
	 zuws(iic,jjc,k)=
     %zuws(iic,jjc,k)/tempo
     %-zuavg(iic,jjc,k)*zwavg(iic,jjc,k)
	 zvws(iic,jjc,k)=
     %zvws(iic,jjc,k)/tempo
     %-zvavg(iic,jjc,k)*zwavg(iic,jjc,k)
	 if(ztempo(iic,jjc,k).ne.0.)
     %zprms(iic,jjc,k)=
     %sqrt(abs(zprms(iic,jjc,k)/ztempo(iic,jjc,k)
     %-zpavg(iic,jjc,k)**2.))
  127	 continue
c
	 mm1=131
	 mm2=132
	 do 128 jjc=1,nsez91
	 do 128 iic=1,nsez9
	 if(myrank.eq.0) then
	    write(mm1,333)'i=',isez9(iic),'____j=',jsez91(jjc),time-0.5*tempo
	    write(mm2,333)'i=',isez9(iic),'____j=',jsez91(jjc),time-0.5*tempo
	    do 129 k=kz1,kz2
	    write(mm1,10)iplant8,isez9(iic),jsez91(jjc),k,
     %zurms(iic,jjc,k),zvrms(iic,jjc,k),zwrms(iic,jjc,k),
     %zprms(iic,jjc,k)
	    write(mm2,10)iplant8,isez9(iic),jsez91(jjc),k,
     %zuvs(iic,jjc,k),zuws(iic,jjc,k),zvws(iic,jjc,k)
  129	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(zurms(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,19,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvrms(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,20,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zwrms(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,21,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zprms(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,22,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zuvs(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,23,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zuws(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,24,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvws(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,25,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 130 k=kz1,kz2
		  kg=jp*(nz-2)+k
		  write(mm1,10)iplant8,isez9(iic),jsez91(jjc),kg,
     %zurms(iic,jjc,k),zvrms(iic,jjc,k),zwrms(iic,jjc,k),
     %zprms(iic,jjc,k)
		  write(mm2,10)iplant8,isez9(iic),jsez91(jjc),kg,
     %zuvs(iic,jjc,k),zuws(iic,jjc,k),zvws(iic,jjc,k)
 130		  continue
	       enddo
	    endif
	 else
	    CALL MPI_SEND(zurms(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,19,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvrms(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,20,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zwrms(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,21,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zprms(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,22,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zuvs(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,23,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zuws(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,24,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvws(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,25,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
 128	 continue
c
	 do jjc=1,nsez91
	    do iic=1,nsez9
	       do k=kz1,kz2
		  zvrms(iic,jjc,k)=0.
		  zurms(iic,jjc,k)=0.
		  zwrms(iic,jjc,k)=0.
		  zprms(iic,jjc,k)=0.
		  zuvs(iic,jjc,k)=0.
		  zuws(iic,jjc,k)=0.
		  zvws(iic,jjc,k)=0.
	       enddo
	    enddo
	 enddo
      endif
 333  format(a2,i6,a6,i6,f20.10)
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine mediecalcstagg
     %(nx,ny,nz,nzg,nbd,time,dt,vo,uo,wo,p,tv,xc,xu,yc,yv,zc,
     %zcg,zwg,flagpo,flaguo,flagvo,flagwo,ntime)
!    %mbd,nfacet,unvect,vertex)
c
      use sections
      use mediecalc
      use mediecalc1
      use mediecalc2
      use mediecalc3
      use vortavg
      include'common.h'
      include'averages.h'
!     include'immersed.h'
c
c Global variables
      integer nx,ny,nz,nzg,nbd,ntime
      integer flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),
     %flagwo(nx,ny,nz,nbd),flagpo(nx,ny,nz,nbd)
      real time,dt
      real vo(nx,ny,nz),uo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz),
     %xc(nx),xu(nx),yc(ny),yv(ny),zc(nz),zcg(nzg),zwg(nzg),
     %tv(nx,ny,nz)
!     integer mbd,nfacet
!     real unvect(3,nfacet),vertex(3,3,nfacet)
c
c Local variables
      integer iic,is,ip,im,i,jjc,jp,jm,j,js,jsy,jpsy,kkc,ks,kp,km,k,
     %iplant4c,iplant4d,iplant4p,iplant5c,iplant5d,iplant5p
      real vstagg,ustagg,wstagg,pstagg,vmod,tvstagg
!     integer ibd
c
      tempo3=tempo3+dt
!
!     Time averaged values and root mean squares
!     at some circumferential sections
!
      do 101 iic=1,nsez1
      is=isez1(iic)
      ip=is+1
      do 100 k=kz1,kz2
      km=k-1
      do 100 j=jy1,jy2
      jm=jmv(j)
      vstagg=(vo(is,j,k)+vo(ip,j,k)+
     %vo(is,jm,k)+vo(ip,jm,k))*0.25
      vplot(iic,j,k)=vplot(iic,j,k)+
     %vstagg*dt
      ustagg=uo(is,j,k)
      uplot(iic,j,k)=uplot(iic,j,k)+
     %ustagg*dt
      wstagg=(wo(is,j,k)+wo(ip,j,k)+
     %wo(is,j,km)+wo(ip,j,km))*0.25
      wplot(iic,j,k)=wplot(iic,j,k)+
     %wstagg*dt
      vrmsplot(iic,j,k)=vrmsplot(iic,j,k)+
     %dt*(vstagg**2.)
      urmsplot(iic,j,k)=urmsplot(iic,j,k)+
     %dt*(ustagg**2.)
      wrmsplot(iic,j,k)=wrmsplot(iic,j,k)+
     %dt*(wstagg**2.)
      if(all(flagpo(is:ip,j,k,:).le.0)) then
	 tempo4(iic,j,k)=tempo4(iic,j,k)+dt
	 pstagg=(p(is,j,k)+p(ip,j,k))*0.5
	 pplot(iic,j,k)=pplot(iic,j,k)+
     %pstagg*dt
	 prmsplot(iic,j,k)=prmsplot(iic,j,k)+
     %dt*(pstagg**2.)
      endif
      vmod=sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
      vmodplot(iic,j,k)=vmodplot(iic,j,k)+
     %vmod*dt
      tvstagg=(tv(is,j,k)+tv(ip,j,k))*0.5
      tvplot(iic,j,k)=tvplot(iic,j,k)+
     %tvstagg*dt
  100 continue
  101 continue
!
!     Time averaged values and root mean squares
!     at some cross sections
!
      do 103 kkc=1,nsez2max
      ks=ksez2(kkc)
      kp=ks+1
      do 102 i=ix1,ix2
      im=i-1
      do 102 j=jy1,jy2
      jm=jmv(j)
      vstagg=(vo(i,j,ks)+vo(i,jm,ks)+
     %vo(i,j,kp)+vo(i,jm,kp))*0.25
      vplot2(i,j,kkc)=vplot2(i,j,kkc)+
     %vstagg*dt
      ustagg=(uo(i,j,ks)+uo(im,j,ks)+
     %uo(i,j,kp)+uo(im,j,kp))*0.25
      uplot2(i,j,kkc)=uplot2(i,j,kkc)+
     %ustagg*dt
      wstagg=wo(i,j,ks)
      wplot2(i,j,kkc)=wplot2(i,j,kkc)+
     %wstagg*dt
      vrmsplot2(i,j,kkc)=vrmsplot2(i,j,kkc)+
     %dt*(vstagg**2.)
      urmsplot2(i,j,kkc)=urmsplot2(i,j,kkc)+
     %dt*(ustagg**2.)
      wrmsplot2(i,j,kkc)=wrmsplot2(i,j,kkc)+
     %dt*(wstagg**2.)
      if(all(flagpo(i,j,ks:kp,:).le.0)) then
	 tempo5(i,j,kkc)=tempo5(i,j,kkc)+dt
	 pstagg=(p(i,j,ks)+p(i,j,kp))*0.5
	 pplot2(i,j,kkc)=pplot2(i,j,kkc)+
     %pstagg*dt
	 prmsplot2(i,j,kkc)=prmsplot2(i,j,kkc)+
     %dt*(pstagg**2.)
      endif
      vmod=sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
      vmodplot2(i,j,kkc)=vmodplot2(i,j,kkc)+
     %vmod*dt
      tvstagg=(tv(i,j,ks)+tv(i,j,kp))*0.5
      tvplot2(i,j,kkc)=tvplot2(i,j,kkc)+
     %tvstagg*dt
  102 continue
  103 continue
!
!     Time averaged values and root mean squares
!     at some meridian sections
!
      do 111 jjc=1,nsez14
      js=jsez14(jjc)
      jp=jpv(js)
      jsy=jsym(js)
      jpsy=jsym(jp)
      do 110 k=kz1,kz2
      km=k-1
      do 110 i=ix1,ix2
      im=i-1
      vstagg=vo(i,js,k)
      vplot14(i,jjc,k)=vplot14(i,jjc,k)+
     %vstagg*dt
      ustagg=(uo(i,js,k)+uo(im,js,k)+
     %uo(i,jp,k)+uo(im,jp,k))*0.25
      uplot14(i,jjc,k)=uplot14(i,jjc,k)+
     %ustagg*dt
      wstagg=(wo(i,js,k)+wo(i,jp,k)+
     %wo(i,js,km)+wo(i,jp,km))*0.25
      wplot14(i,jjc,k)=wplot14(i,jjc,k)+
     %wstagg*dt
      vrmsplot14(i,jjc,k)=vrmsplot14(i,jjc,k)+
     %dt*(vstagg**2.)
      urmsplot14(i,jjc,k)=urmsplot14(i,jjc,k)+
     %dt*(ustagg**2.)
      wrmsplot14(i,jjc,k)=wrmsplot14(i,jjc,k)+
     %dt*(wstagg**2.)
      if(all(flagpo(i,js:jp,k,:).le.0)) then
	 tempo14(i,jjc,k)=tempo14(i,jjc,k)+dt
	 pstagg=(p(i,js,k)+p(i,jp,k))*0.5
	 pplot14(i,jjc,k)=pplot14(i,jjc,k)+
     %pstagg*dt
	 prmsplot14(i,jjc,k)=prmsplot14(i,jjc,k)+
     %dt*(pstagg**2.)
      endif
      vmod=sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
      vmodplot14(i,jjc,k)=vmodplot14(i,jjc,k)+
     %vmod*dt
      tvstagg=(tv(i,js,k)+tv(i,jp,k))*0.5
      tvplot14(i,jjc,k)=tvplot14(i,jjc,k)+
     %tvstagg*dt
c
      vstagg=vo(i,jsy,k)
      vplot14(i,nsez14+jjc,k)=vplot14(i,nsez14+jjc,k)+
     %vstagg*dt
      ustagg=(uo(i,jsy,k)+uo(im,jsy,k)+
     %uo(i,jpsy,k)+uo(im,jpsy,k))*0.25
      uplot14(i,nsez14+jjc,k)=uplot14(i,nsez14+jjc,k)+
     %ustagg*dt
      wstagg=(wo(i,jsy,k)+wo(i,jpsy,k)+
     %wo(i,jsy,km)+wo(i,jpsy,km))*0.25
      wplot14(i,nsez14+jjc,k)=wplot14(i,nsez14+jjc,k)+
     %wstagg*dt
      vrmsplot14(i,nsez14+jjc,k)=vrmsplot14(i,nsez14+jjc,k)+
     %dt*(vstagg**2.)
      urmsplot14(i,nsez14+jjc,k)=urmsplot14(i,nsez14+jjc,k)+
     %dt*(ustagg**2.)
      wrmsplot14(i,nsez14+jjc,k)=wrmsplot14(i,nsez14+jjc,k)+
     %dt*(wstagg**2.)
      if(all(flagpo(i,jsy:jpsy,k,:).le.0)) then
	 tempo14(i,nsez14+jjc,k)=tempo14(i,nsez14+jjc,k)+dt
	 pstagg=(p(i,jsy,k)+p(i,jpsy,k))*0.5
	 pplot14(i,nsez14+jjc,k)=pplot14(i,nsez14+jjc,k)+
     %pstagg*dt
	 prmsplot14(i,nsez14+jjc,k)=prmsplot14(i,nsez14+jjc,k)+
     %dt*(pstagg**2.)
      endif
      vmod=sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
      vmodplot14(i,nsez14+jjc,k)=vmodplot14(i,nsez14+jjc,k)+
     %vmod*dt
      tvstagg=(tv(i,jsy,k)+tv(i,jpsy,k))*0.5
      tvplot14(i,nsez14+jjc,k)=tvplot14(i,nsez14+jjc,k)+
     %tvstagg*dt
  110 continue
  111 continue
!
!!!	 if((time-timep).gt.((periodo)*(iplant5+1))) then
c      if((abs(romega)*time_plot2*180./pig).gt.(rplot1*(iplant5+1)))
	 call reynoldsstagg(1,nx,ny,nz,uo,vo,wo,dt)
	 call avgvort2stagg(1,nx,ny,nz,dt,xc,xu,zc,uo,vo,wo)
!!!	 endif
!
      if(tempo3.ge.periodo2) then
c
	 do 105 iic=1,nsez1
	 do 104 k=kz1,kz2
	 do 104 j=jy1,jy2
	    vplot(iic,j,k)=vplot(iic,j,k)/tempo3
	    uplot(iic,j,k)=uplot(iic,j,k)/tempo3
	    wplot(iic,j,k)=wplot(iic,j,k)/tempo3
	    vrmsplot(iic,j,k)=
     %vrmsplot(iic,j,k)/tempo3-vplot(iic,j,k)**2.
	    urmsplot(iic,j,k)=
     %urmsplot(iic,j,k)/tempo3-uplot(iic,j,k)**2.
	    wrmsplot(iic,j,k)=
     %wrmsplot(iic,j,k)/tempo3-wplot(iic,j,k)**2.
	    if(tempo4(iic,j,k).ne.0.) then
	       pplot(iic,j,k)=pplot(iic,j,k)/tempo4(iic,j,k)
	       prmsplot(iic,j,k)=
     %prmsplot(iic,j,k)/tempo4(iic,j,k)-pplot(iic,j,k)**2.
	    endif
	    vmodplot(iic,j,k)=vmodplot(iic,j,k)/tempo3
	    tvplot(iic,j,k)=tvplot(iic,j,k)/tempo3
	    tvplot(iic,j,k)=tvplot(iic,j,k)/ru1
  104	 continue
  105	 continue
c
	 do 107 kkc=1,nsez2max
	 do 106 i=ix1,ix2
	 do 106 j=jy1,jy2
	    vplot2(i,j,kkc)=vplot2(i,j,kkc)/tempo3
	    uplot2(i,j,kkc)=uplot2(i,j,kkc)/tempo3
	    wplot2(i,j,kkc)=wplot2(i,j,kkc)/tempo3
	    vrmsplot2(i,j,kkc)=
     %vrmsplot2(i,j,kkc)/tempo3-vplot2(i,j,kkc)**2.
	    urmsplot2(i,j,kkc)=
     %urmsplot2(i,j,kkc)/tempo3-uplot2(i,j,kkc)**2.
	    wrmsplot2(i,j,kkc)=
     %wrmsplot2(i,j,kkc)/tempo3-wplot2(i,j,kkc)**2.
	    if(tempo5(i,j,kkc).ne.0.) then
	       pplot2(i,j,kkc)=pplot2(i,j,kkc)/tempo5(i,j,kkc)
	       prmsplot2(i,j,kkc)=
     %prmsplot2(i,j,kkc)/tempo5(i,j,kkc)-pplot2(i,j,kkc)**2.
	    endif
	    vmodplot2(i,j,kkc)=vmodplot2(i,j,kkc)/tempo3
	    tvplot2(i,j,kkc)=tvplot2(i,j,kkc)/tempo3
	    tvplot2(i,j,kkc)=tvplot2(i,j,kkc)/ru1
  106	 continue
  107	 continue
c
	 do 115 jjc=1,2*nsez14
	 do 114 k=kz1,kz2
	 do 114 i=ix1,ix2
	    vplot14(i,jjc,k)=vplot14(i,jjc,k)/tempo3
	    uplot14(i,jjc,k)=uplot14(i,jjc,k)/tempo3
	    wplot14(i,jjc,k)=wplot14(i,jjc,k)/tempo3
	    vrmsplot14(i,jjc,k)=
     %vrmsplot14(i,jjc,k)/tempo3-vplot14(i,jjc,k)**2.
	    urmsplot14(i,jjc,k)=
     %urmsplot14(i,jjc,k)/tempo3-uplot14(i,jjc,k)**2.
	    wrmsplot14(i,jjc,k)=
     %wrmsplot14(i,jjc,k)/tempo3-wplot14(i,jjc,k)**2.
	    if(tempo14(i,jjc,k).ne.0.) then
	       pplot14(i,jjc,k)=
     %pplot14(i,jjc,k)/tempo14(i,jjc,k)
	       prmsplot14(i,jjc,k)=
     %prmsplot14(i,jjc,k)/tempo14(i,jjc,k)-pplot14(i,jjc,k)**2.
	    endif
	    vmodplot14(i,jjc,k)=vmodplot14(i,jjc,k)/tempo3
	    tvplot14(i,jjc,k)=tvplot14(i,jjc,k)/tempo3
	    tvplot14(i,jjc,k)=tvplot14(i,jjc,k)/ru1
  114	 continue
  115	 continue
!
	 iplant4=iplant4+1
	 call kplotstagg(iplant4c,iplant4d,iplant4p,
     %nx,ny,nz,nzg,nbd,ntime,xu,xc,yv,yc,zwg,zcg,
     %flaguo,flagvo,flagwo)
!
ccc 13	    format(2x,a12,1x,e14.7,1x,e14.7,1x,e14.7)
ccc 14	    format(6x,a6,3x,e14.7,1x,e14.7,1x,e14.7)
ccc	    if((ibm.gt.1).and.(myrank.eq.0)) then
ccc	       do ibd=mbd,nbd
ccc		  open(83,file='ib_n.'//index(ibd)//'_'
ccc	%//char(48+iplant4c)//char(48+iplant4d)//char(48+iplant4p)//
ccc	%'_kplot1.stl')
ccc		  write(83,*)'solid  OBJECT'
ccc		  do j=lb(ibd)+1,lb(ibd)+mb(ibd)
ccc		     write(83,13)'facet normal',unvect(1,j),unvect(2,j),
ccc	%unvect(3,j)
ccc		     write(83,*)'outer loop'
ccc		     write(83,14)
ccc	%'vertex',vertex(1,1,j),vertex(2,1,j),vertex(3,1,j)
ccc		     write(83,14)
ccc	%'vertex',vertex(1,2,j),vertex(2,2,j),vertex(3,2,j)
ccc		     write(83,14)
ccc	%'vertex',vertex(1,3,j),vertex(2,3,j),vertex(3,3,j)
ccc		     write(83,*)'endloop'
ccc		     write(83,*)'endfacet'
ccc		  enddo
ccc		  write(83,*)'endsolid	OBJECT'
ccc		  close(83)
ccc	       enddo
ccc	    endif
!
!!!	    if((time-timep).gt.((periodo)*(iplant5+1))) then
c	  if((abs(romega)*time_plot2*180./pig).gt.
c     %(rplot1*(iplant5+1))) then
	    iplant5=iplant5+1
	    call reynoldsstagg(2,nx,ny,nz,uo,vo,wo,dt)
	    call kplot2(iplant5c,iplant5d,iplant5p,
     %nx,ny,nz,nzg,nbd,ntime,flaguo,flagvo,flagwo)
	    call avgvort2stagg(2,nx,ny,nz,dt,xc,xu,zc,uo,vo,wo)
	    call plotavgvort(iplant5c,iplant5d,iplant5p,
     %nx,ny,nz,nzg,nbd,ntime,flaguo,flagvo,flagwo)
c	     call vort2
c	     call plotvort(iplant5c,iplant5d,iplant5p)
!
ccc	    if((ibm.gt.1).and.(myrank.eq.0)) then
ccc	       do ibd=mbd,nbd
ccc		  open(83,file='ib_n.'//index(ibd)//'_'
ccc	%//char(48+iplant4c)//char(48+iplant4d)//char(48+iplant4p)//
ccc	%'_kplot2.stl')
ccc		  write(83,*)'solid  OBJECT'
ccc		  do j=lb(ibd)+1,lb(ibd)+mb(ibd)
ccc		     write(83,13)'facet normal',unvect(1,j),unvect(2,j),
ccc	%unvect(3,j)
ccc		     write(83,*)'outer loop'
ccc		     write(83,14)
ccc	%'vertex',vertex(1,1,j),vertex(2,1,j),vertex(3,1,j)
ccc		     write(83,14)
ccc	%'vertex',vertex(1,2,j),vertex(2,2,j),vertex(3,2,j)
ccc		     write(83,14)
ccc	%'vertex',vertex(1,3,j),vertex(2,3,j),vertex(3,3,j)
ccc		     write(83,*)'endloop'
ccc		     write(83,*)'endfacet'
ccc		  enddo
ccc		  write(83,*)'endsolid	OBJECT'
ccc		  close(83)
ccc	       enddo
ccc	    endif
!
!!!	    endif
!
	 do iic=1,nsez1
	    do k=kz1,kz2
	       do j=jy1,jy2
		  vrmsplot(iic,j,k)=0.
		  urmsplot(iic,j,k)=0.
		  wrmsplot(iic,j,k)=0.
		  prmsplot(iic,j,k)=0.
		  vplot(iic,j,k)=0.
		  uplot(iic,j,k)=0.
		  wplot(iic,j,k)=0.
		  pplot(iic,j,k)=0.
		  tempo4(iic,j,k)=0.
		  uvplot(iic,j,k)=0.
		  vwplot(iic,j,k)=0.
		  uwplot(iic,j,k)=0.
		  voraz1med(iic,j,k)=0.
		  vorr1med(iic,j,k)=0.
		  vorz1med(iic,j,k)=0.
		  vmodplot(iic,j,k)=0.
		  tvplot(iic,j,k)=0.
	       enddo
	    enddo
	 enddo
c
	 do kkc=1,nsez2max
	    do i=ix1,ix2
	       do j=jy1,jy2
		  vrmsplot2(i,j,kkc)=0.
		  urmsplot2(i,j,kkc)=0.
		  wrmsplot2(i,j,kkc)=0.
		  prmsplot2(i,j,kkc)=0.
		  vplot2(i,j,kkc)=0.
		  uplot2(i,j,kkc)=0.
		  wplot2(i,j,kkc)=0.
		  pplot2(i,j,kkc)=0.
		  tempo5(i,j,kkc)=0.
		  uvplot2(i,j,kkc)=0.
		  vwplot2(i,j,kkc)=0.
		  uwplot2(i,j,kkc)=0.
		  voraz2med(i,j,kkc)=0.
		  vorr2med(i,j,kkc)=0.
		  vorz2med(i,j,kkc)=0.
		  vmodplot2(i,j,kkc)=0.
		  tvplot2(i,j,kkc)=0.
	       enddo
	    enddo
	 enddo
c
	 do jjc=1,2*nsez14
	    do k=kz1,kz2
	       do i=ix1,ix2
		  vrmsplot14(i,jjc,k)=0.
		  urmsplot14(i,jjc,k)=0.
		  wrmsplot14(i,jjc,k)=0.
		  prmsplot14(i,jjc,k)=0.
		  vplot14(i,jjc,k)=0.
		  uplot14(i,jjc,k)=0.
		  wplot14(i,jjc,k)=0.
		  pplot14(i,jjc,k)=0.
		  tempo14(i,jjc,k)=0.
		  uvplot14(i,jjc,k)=0.
		  vwplot14(i,jjc,k)=0.
		  uwplot14(i,jjc,k)=0.
		  voraz14med(i,jjc,k)=0.
		  vorr14med(i,jjc,k)=0.
		  vorz14med(i,jjc,k)=0.
		  vmodplot14(i,jjc,k)=0.
		  tvplot14(i,jjc,k)=0.
	       enddo
	    enddo
	 enddo
	 tempo3=0.
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine reynoldsstagg(ind,nx,ny,nz,uo,vo,wo,dt)
c
      use sections
      use mediecalc1
      use mediecalc2
      use mediecalc3
      include'common.h'
      include'averages.h'
c
c Global variables
      integer ind,nx,ny,nz
      real dt,uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
c
c Local variables
      integer iic,i,is,ip,im,jjc,j,js,jp,jsy,jpsy,jm,kkc,k,ks,kp,km
      real ustagg,vstagg,wstagg
c
      if(ind.eq.2) goto 200
!
!     Circumferential sections
!
      do 101 iic=1,nsez1
      is=isez1(iic)
      ip=is+1
      do 100 k=kz1,kz2
      km=k-1
      do 100 j=jy1,jy2
      jm=jmv(j)
      ustagg=
     %uo(is,j,k)
      vstagg=
     %(vo(is,j,k)+vo(ip,j,k)+vo(is,jm,k)+vo(ip,jm,k))*0.25
      wstagg=
     %(wo(is,j,k)+wo(ip,j,k)+wo(is,j,km)+wo(ip,j,km))*0.25
      uvplot(iic,j,k)=uvplot(iic,j,k)+
     %ustagg*vstagg*dt
      uwplot(iic,j,k)=uwplot(iic,j,k)+
     %ustagg*wstagg*dt
      vwplot(iic,j,k)=vwplot(iic,j,k)+
     %vstagg*wstagg*dt
  100 continue
  101 continue
!
!     Cross sections
!
      do 103 kkc=1,nsez2max
      ks=ksez2(kkc)
      kp=ks+1
      do 102 i=ix1,ix2
      im=i-1
      do 102 j=jy1,jy2
      jm=jmv(j)
      ustagg=
     %(uo(i,j,ks)+uo(im,j,ks)+uo(i,j,kp)+uo(im,j,kp))*0.25
      vstagg=
     %(vo(i,j,ks)+vo(i,jm,ks)+vo(i,j,kp)+vo(i,jm,kp))*0.25
      wstagg=
     %wo(i,j,ks)
      uvplot2(i,j,kkc)=uvplot2(i,j,kkc)+
     %ustagg*vstagg*dt
      uwplot2(i,j,kkc)=uwplot2(i,j,kkc)+
     %ustagg*wstagg*dt
      vwplot2(i,j,kkc)=vwplot2(i,j,kkc)+
     %vstagg*wstagg*dt
  102 continue
  103 continue
!
!     Meridian sections
!
      do 111 jjc=1,nsez14
      js=jsez14(jjc)
      jp=jpv(js)
      jsy=jsym(js)
      jpsy=jsym(jp)
      do 110 k=kz1,kz2
      km=k-1
      do 110 i=ix1,ix2
      im=i-1
      ustagg=
     %(uo(i,js,k)+uo(im,js,k)+uo(i,jp,k)+uo(im,jp,k))*0.25
      vstagg=
     %vo(i,js,k)
      wstagg=
     %(wo(i,js,k)+wo(i,jp,k)+wo(i,js,km)+wo(i,jp,km))*0.25
      uvplot14(i,jjc,k)=uvplot14(i,jjc,k)+
     %ustagg*vstagg*dt
      uwplot14(i,jjc,k)=uwplot14(i,jjc,k)+
     %ustagg*wstagg*dt
      vwplot14(i,jjc,k)=vwplot14(i,jjc,k)+
     %vstagg*wstagg*dt
c
      ustagg=
     %(uo(i,jsy,k)+uo(im,jsy,k)+uo(i,jpsy,k)+uo(im,jpsy,k))*0.25
      vstagg=
     %vo(i,jsy,k)
      wstagg=
     %(wo(i,jsy,k)+wo(i,jpsy,k)+wo(i,jsy,km)+wo(i,jpsy,km))*0.25
      uvplot14(i,nsez14+jjc,k)=uvplot14(i,nsez14+jjc,k)+
     %ustagg*vstagg*dt
      uwplot14(i,nsez14+jjc,k)=uwplot14(i,nsez14+jjc,k)+
     %ustagg*wstagg*dt
      vwplot14(i,nsez14+jjc,k)=vwplot14(i,nsez14+jjc,k)+
     %vstagg*wstagg*dt
  110 continue
  111 continue
      return
  200 continue
!
!     Circumferential sections
!
      do 105 iic=1,nsez1
      do 104 k=kz1,kz2
      do 104 j=jy1,jy2
	 uvplot(iic,j,k)=
     %uvplot(iic,j,k)/tempo3-uplot(iic,j,k)*vplot(iic,j,k)
	 uwplot(iic,j,k)=
     %uwplot(iic,j,k)/tempo3-uplot(iic,j,k)*wplot(iic,j,k)
	 vwplot(iic,j,k)=
     %vwplot(iic,j,k)/tempo3-vplot(iic,j,k)*wplot(iic,j,k)
  104 continue
  105 continue
!
!     Cross sections
!
      do 107 kkc=1,nsez2max
      do 106 i=ix1,ix2
      do 106 j=jy1,jy2
	 uvplot2(i,j,kkc)=
     %uvplot2(i,j,kkc)/tempo3-uplot2(i,j,kkc)*vplot2(i,j,kkc)
	 uwplot2(i,j,kkc)=
     %uwplot2(i,j,kkc)/tempo3-uplot2(i,j,kkc)*wplot2(i,j,kkc)
	 vwplot2(i,j,kkc)=
     %vwplot2(i,j,kkc)/tempo3-vplot2(i,j,kkc)*wplot2(i,j,kkc)
  106 continue
  107 continue
!
!     Meridian sections
!
      do 115 jjc=1,2*nsez14
      do 114 k=kz1,kz2
      do 114 i=ix1,ix2
	 uvplot14(i,jjc,k)=
     %uvplot14(i,jjc,k)/tempo3-uplot14(i,jjc,k)*vplot14(i,jjc,k)
	 uwplot14(i,jjc,k)=
     %uwplot14(i,jjc,k)/tempo3-uplot14(i,jjc,k)*wplot14(i,jjc,k)
	 vwplot14(i,jjc,k)=
     %vwplot14(i,jjc,k)/tempo3-vplot14(i,jjc,k)*wplot14(i,jjc,k)
  114 continue
  115 continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine avgvort2stagg
     %(ind,nx,ny,nz,dt,xc,xu,zc,uo,vo,wo)
c
      use sections
      use vortavg
      include'common.h'
      include'averages.h'
c
c Global variables
      integer ind,nx,ny,nz
      real dt,xc(nx),xu(nx),zc(nz),
     %uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
c
c Local variables
      integer iic,ic,ip,im,jjc,jc,jp,jm,jsy,jpsy,kkc,kc,kp,km,
     %i,j,k
      real dtheta,dr,dz,dq1x3,dq3x1,dq3x2,dq2x3,dq1x2,dq2x1
c
      if(ind.eq.2) goto 200
c
c    ! Time-averaged azimuthal vorticity - Circumferential sections !
      do 101 kc=kz1,kz2
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	do 101 iic=1,nsez1
	  ic=isez1(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 101 jc=jy1,jy2
	    dq1x3=(uo(ic,jc,kp)
     %-uo(ic,jc,km))/dz
	    dq3x1=(0.5*(wo(ip,jc,kc)+wo(ip,jc,km))
     %-0.5*(wo(ic,jc,kc)+wo(ic,jc,km)))/dr
	    voraz1med(iic,jc,kc)=
     %voraz1med(iic,jc,kc)+(dq1x3-dq3x1)*dt
  101 continue
c
c    ! Time-averaged azimuthal vorticity - Cross sections !
      do 111 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
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
	  voraz2med(ic,jc,kkc)=
     %voraz2med(ic,jc,kkc)+(dq1x3-dq3x1)*dt
  112 continue
	do 113 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 113 jc=jy1,jy2
	    dq1x3=(0.5*(uo(ic,jc,kp)+uo(im,jc,kp))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dz
	    dq3x1=(wo(ip,jc,kc)
     %-wo(im,jc,kc))/dr
	    voraz2med(ic,jc,kkc)=
     %voraz2med(ic,jc,kkc)+(dq1x3-dq3x1)*dt
  113 continue
  111 continue
c
c    ! Time-averaged azimuthal vorticity - Meridian sections !
      do 121 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jsy=jsym(jc)
	jpsy=jsym(jp)
	do 121 kc=kz1,kz2
	  kp=kc+1
	  km=kc-1
	  dz=zc(kp)-zc(km)
	  ic=ix1
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  dq1x3=(0.25*(uo(ic,jc,kp)+uo(im,jc,kp)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,km)+uo(im,jc,km)
     %+uo(ic,jp,km)+uo(im,jp,km)))/dz
	  dq3x1=(0.25*(wo(ip,jc,kc)+wo(ip,jp,kc)
     %+wo(ip,jc,km)+wo(ip,jp,km))
     %-0.25*(wo(ic,jsy,kc)+wo(ic,jpsy,kc)
     %+wo(ic,jsy,km)+wo(ic,jpsy,km)))/dr
	  voraz14med(ic,jjc,kc)=
     %voraz14med(ic,jjc,kc)+(dq1x3-dq3x1)*dt
	  dq1x3=(0.25*(uo(ic,jsy,kp)+uo(im,jsy,kp)
     %+uo(ic,jpsy,kp)+uo(im,jpsy,kp))
     %-0.25*(uo(ic,jsy,km)+uo(im,jsy,km)
     %+uo(ic,jpsy,km)+uo(im,jpsy,km)))/dz
	  dq3x1=(0.25*(wo(ip,jsy,kc)+wo(ip,jpsy,kc)
     %+wo(ip,jsy,km)+wo(ip,jpsy,km))
     %-0.25*(wo(ic,jc,kc)+wo(ic,jp,kc)
     %+wo(ic,jc,km)+wo(ic,jp,km)))/dr
	  voraz14med(ic,nsez14+jjc,kc)=
     %voraz14med(ic,nsez14+jjc,kc)+(dq1x3-dq3x1)*dt
	  do 121 ic=ix1+1,ix2
	    ip=ic+1
	    im=ic-1
	    dr=xc(ip)-xc(im)
	    dq1x3=(0.25*(uo(ic,jc,kp)+uo(im,jc,kp)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,km)+uo(im,jc,km)
     %+uo(ic,jp,km)+uo(im,jp,km)))/dz
	    dq3x1=(0.25*(wo(ip,jc,kc)+wo(ip,jp,kc)
     %+wo(ip,jc,km)+wo(ip,jp,km))
     %-0.25*(wo(im,jc,kc)+wo(im,jp,kc)
     %+wo(im,jc,km)+wo(im,jp,km)))/dr
	    voraz14med(ic,jjc,kc)=
     %voraz14med(ic,jjc,kc)+(dq1x3-dq3x1)*dt
	    dq1x3=(0.25*(uo(ic,jsy,kp)+uo(im,jsy,kp)
     %+uo(ic,jpsy,kp)+uo(im,jpsy,kp))
     %-0.25*(uo(ic,jsy,km)+uo(im,jsy,km)
     %+uo(ic,jpsy,km)+uo(im,jpsy,km)))/dz
	    dq3x1=(0.25*(wo(ip,jsy,kc)+wo(ip,jpsy,kc)
     %+wo(ip,jsy,km)+wo(ip,jpsy,km))
     %-0.25*(wo(im,jsy,kc)+wo(im,jpsy,kc)
     %+wo(im,jsy,km)+wo(im,jpsy,km)))/dr
	    voraz14med(ic,nsez14+jjc,kc)=
     %voraz14med(ic,nsez14+jjc,kc)+(dq1x3-dq3x1)*dt
  121 continue
c
c    ! Time-averaged radial vorticity - Circumferential sections !
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
	    dq3x2=(0.25*(wo(ic,jp,kc)+wo(ip,jp,kc)
     %+wo(ic,jp,km)+wo(ip,jp,km))
     %-0.25*(wo(ic,jm,kc)+wo(ip,jm,kc)
     %+wo(ic,jm,km)+wo(ip,jm,km)))/dtheta
	    dq2x3=(0.25*(vo(ic,jc,kp)+vo(ip,jc,kp)
     %+vo(ic,jm,kp)+vo(ip,jm,kp))
     %-0.25*(vo(ic,jc,km)+vo(ip,jc,km)
     %+vo(ic,jm,km)+vo(ip,jm,km)))/dz
	    vorr1med(iic,jc,kc)=
     %vorr1med(iic,jc,kc)+(dq3x2/xu(ic)-dq2x3)*dt
  131 continue
c
c    ! Time-averaged radial vorticity - Cross sections !
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
     %-wo(ic,jm,kc))/dtheta
	    dq2x3=(0.5*(vo(ic,jc,kp)+vo(ic,jm,kp))
     %-0.5*(vo(ic,jc,kc)+vo(ic,jm,kc)))/dz
	    vorr2med(ic,jc,kkc)=
     %vorr2med(ic,jc,kkc)+(dq3x2/xc(ic)-dq2x3)*dt
  141 continue
c
c    ! Time-averaged radial vorticity - Meridian sections !
      do 151 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jsy=jsym(jc)
	jpsy=jsym(jp)
	dtheta=dely
	do 151 kc=kz1,kz2
	  kp=kc+1
	  km=kc-1
	  dz=zc(kp)-zc(km)
	  do 151 ic=ix1,ix2
	    dq3x2=(0.5*(wo(ic,jp,kc)+wo(ic,jp,km))
     %-0.5*(wo(ic,jc,kc)+wo(ic,jc,km)))/dtheta
	    dq2x3=(vo(ic,jc,kp)
     %-vo(ic,jc,km))/dz
	    vorr14med(ic,jjc,kc)=
     %vorr14med(ic,jjc,kc)+(dq3x2/xc(ic)-dq2x3)*dt
	    dq3x2=(0.5*(wo(ic,jpsy,kc)+wo(ic,jpsy,km))
     %-0.5*(wo(ic,jsy,kc)+wo(ic,jsy,km)))/dtheta
	    dq2x3=(vo(ic,jsy,kp)
     %-vo(ic,jsy,km))/dz
	    vorr14med(ic,nsez14+jjc,kc)=
     %vorr14med(ic,nsez14+jjc,kc)+(dq3x2/xc(ic)-dq2x3)*dt
  151 continue
c
c    ! Time-averaged axial vorticity - Circumferential sections !
      do 161 kc=kz1,kz2
	do 161 iic=1,nsez1
	  ic=isez1(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 161 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.5*(vo(ip,jc,kc)+vo(ip,jm,kc))*xc(ip)
     %-0.5*(vo(ic,jc,kc)+vo(ic,jm,kc))*xc(ic))/dr
	    dq1x2=(uo(ic,jp,kc)
     %-uo(ic,jm,kc))/dtheta
	    vorz1med(iic,jc,kc)=
     %vorz1med(iic,jc,kc)+((dq2x1-dq1x2)/xu(ic))*dt
  161 continue
c
c    ! Time-averaged axial vorticity - Cross sections !
      do 171 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
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
     %+vo(ic,jm,kp)+vo(ip,jm,kp))*xu(ic))/dr
	  dq1x2=(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jm,kc)+uo(im,jm,kc)
     %+uo(ic,jm,kp)+uo(im,jm,kp)))/dtheta
	  vorz2med(ic,jc,kkc)=
     %vorz2med(ic,jc,kkc)+((dq2x1-dq1x2)/xc(ic))*dt
  172 continue
c
	do 173 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 173 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.25*(vo(ip,jc,kc)+vo(ip,jm,kc)
     %+vo(ip,jc,kp)+vo(ip,jm,kp))*xc(ip)
     %-0.25*(vo(im,jc,kc)+vo(im,jm,kc)
     %+vo(im,jc,kp)+vo(im,jm,kp))*xc(im))/dr
	    dq1x2=
     %(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jm,kc)+uo(im,jm,kc)
     %+uo(ic,jm,kp)+uo(im,jm,kp)))/dtheta
	    vorz2med(ic,jc,kkc)=
     %vorz2med(ic,jc,kkc)+((dq2x1-dq1x2)/xc(ic))*dt
  173 continue
  171 continue
c
c    ! Time-averaged axial vorticity - Meridian sections !
      do 181 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jsy=jsym(jc)
	jpsy=jsym(jp)
	dtheta=dely
	do 181 kc=kz1,kz2
	  ic=ix1
	  im=ic-1
	  ip=ic+1
	  dr=xu(ic)-xu(im)
	  dq2x1=(0.5*(vo(ic,jc,kc)+vo(ip,jc,kc))*xu(ic))/dr
	  dq1x2=(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dtheta
	  vorz14med(ic,jjc,kc)=
     %vorz14med(ic,jjc,kc)+((dq2x1-dq1x2)/xc(ic))*dt
	  dq2x1=(0.5*(vo(ic,jsy,kc)+vo(ip,jsy,kc))*xu(ic))/dr
	  dq1x2=(0.5*(uo(ic,jpsy,kc)+uo(im,jpsy,kc))
     %-0.5*(uo(ic,jsy,kc)+uo(im,jsy,kc)))/dtheta
	  vorz14med(ic,nsez14+jjc,kc)=
     %vorz14med(ic,nsez14+jjc,kc)+((dq2x1-dq1x2)/xc(ic))*dt
c
	  do 181 ic=ix1+1,ix2
	    ip=ic+1
	    im=ic-1
	    dr=xc(ip)-xc(im)
	    dq2x1=(vo(ip,jc,kc)*xc(ip)
     %-vo(im,jc,kc)*xc(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dtheta
	    vorz14med(ic,jjc,kc)=
     %vorz14med(ic,jjc,kc)+((dq2x1-dq1x2)/xc(ic))*dt
	    dq2x1=(vo(ip,jsy,kc)*xc(ip)
     %-vo(im,jsy,kc)*xc(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jpsy,kc)+uo(im,jpsy,kc))
     %-0.5*(uo(ic,jsy,kc)+uo(im,jsy,kc)))/dtheta
	    vorz14med(ic,nsez14+jjc,kc)=
     %vorz14med(ic,nsez14+jjc,kc)+((dq2x1-dq1x2)/xc(ic))*dt
  181 continue
      return
c
  200 continue
c
!
!     Circumferential sections
!
      do 205 iic=1,nsez1
      do 204 k=kz1,kz2
      do 204 j=jy1,jy2
	 voraz1med(iic,j,k)=voraz1med(iic,j,k)/tempo3
	 vorr1med(iic,j,k)=vorr1med(iic,j,k)/tempo3
	 vorz1med(iic,j,k)=vorz1med(iic,j,k)/tempo3
  204 continue
  205 continue
c
!
!     Cross sections
!
      do 207 kkc=1,nsez2max
      do 206 i=ix1,ix2
      do 206 j=jy1,jy2
	 voraz2med(i,j,kkc)=voraz2med(i,j,kkc)/tempo3
	 vorr2med(i,j,kkc)=vorr2med(i,j,kkc)/tempo3
	 vorz2med(i,j,kkc)=vorz2med(i,j,kkc)/tempo3
  206 continue
  207 continue
c
!
!     Meridian sections
!
      do 215 jjc=1,2*nsez14
      do 214 k=kz1,kz2
      do 214 i=ix1,ix2
	 voraz14med(i,jjc,k)=voraz14med(i,jjc,k)/tempo3
	 vorr14med(i,jjc,k)=vorr14med(i,jjc,k)/tempo3
	 vorz14med(i,jjc,k)=vorz14med(i,jjc,k)/tempo3
  214 continue
  215 continue
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine kplotstagg(iplant4c,iplant4d,iplant4p,
     %nx,ny,nz,nzg,nbd,ntime,xu,xc,yv,yc,zwg,zcg,flaguo,flagvo,
     %flagwo)
c
      use sections
      use mediecalc1
      use mediecalc2
      use mediecalc3
      include'common.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer nx,ny,nz,nzg,nbd,iplant4p,iplant4d,iplant4c,ntime,
     %flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),flagwo(nx,ny,nz,nbd)
      real xu(nx),xc(nx),yv(ny),yc(ny),zwg(nzg),zcg(nzg)
c
c Local variables
      integer iic,iicd,iicu,is,jjc,jjcd,jjcu,js,jsy,kkc,kkcd,kkcu,ks,
     %i,j,k,l,i1,jjp,ir,jr,kr,lr,kk,kkd,kku
!     integer kk1,kk2,krg,flaguor(ny,nzg,nbd),flagvor2(nx,nzg,nbd),
!    %flagvor2sym(nx,nzg,nbd)
      integer STATUS(MPI_STATUS_SIZE)
      real recvar(ny,nz),recvar2(nx,nz),recvar2sym(nx,nz)
c
c Parameters
      integer lamb2,nxp,nyp,nzp
      real deltasect
      parameter (lamb2=1,nyp=1,nxp=1,nzp=1)
      parameter (deltasect=2.*pi)
c
   8  format(10f25.8)
   9  format(10f25.8)
 101  format(3i8)
 102  format(i8,3i2)
 103  format(i2)
c
      i1=1
      if(iplant4.eq.1) then
	 if(myrank.eq.0) then
	    do iic=1,nsez1
	       iicd=iic/10
	       iicu=mod(iic,10)
	       open
     %(81,file='grid1_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'_kplot.xyz')
	       open
     %(82,file='grid1_new_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'_kplot.xyz')
	       is=isez1(iic)
	       write(81,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	       write(82,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	       do 1 k=2,nzg-1,nzp
	       do 10 l=1,lamb2
	       do 10 j=jy1,jy2,nyp
	       write(81,9)xu(is)*cos(yc(jy1))
   10	       write(82,9)xu(is)*cos(yc(j)+(l-1)*2.*pi/lamb2)
	       write(81,9)xu(is)*cos(yc(jy1))
    1	       write(82,9)xu(is)*cos(yc(jy1)+lamb2*deltasect)
	       do 2 k=2,nzg-1,nzp
	       do 20 l=1,lamb2
	       do 20 j=jy1,jy2,nyp
	       write(81,9)
     %xu(is)*sin(yc(jy1))+xu(is)*(yc(j)+(l-1)*2.*pi/lamb2)
   20	       write(82,9)xu(is)*sin(yc(j)+(l-1)*2.*pi/lamb2)
	       write(81,9)
     %xu(is)*sin(yc(jy1))+xu(is)*deltasect*lamb2
   2	       write(82,9)xu(is)*sin(yc(jy1)+lamb2*deltasect)
	       do 3 k=2,nzg-1,nzp
	       do 30 l=1,lamb2
	       do 30 j=jy1,jy2,nyp
	       write(81,9)zcg(k)
   30	       write(82,9)zcg(k)
	       write(81,9)zcg(k)
   3	       write(82,9)zcg(k)
	       close(81)
	       close(82)
	    enddo
	    do kkc=1,nsez2
	       kkcd=kkc/10
	       kkcu=mod(kkc,10)
	       open(81,file='grid2_2D_slice'//char(48+kkcd)
     %//char(48+kkcu)//'_kplot.xyz')
	       ks=ksez2g(kkc)
	       write(81,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	       do 11 i=ix1,ix2,nxp
	       do 110 l=1,lamb2
	       do 110 j=jy1,jy2,nyp
 110	       write(81,9)xc(i)*cos(yc(j)+(l-1)*2.*pi/lamb2)
  11	       write(81,9)xc(i)*cos(yc(jy1)+lamb2*deltasect)
	       do 12 i=ix1,ix2,nxp
	       do 120 l=1,lamb2
	       do 120 j=jy1,jy2,nyp
 120	       write(81,9)xc(i)*sin(yc(j)+(l-1)*2.*pi/lamb2)
  12	       write(81,9)xc(i)*sin(yc(jy1)+lamb2*deltasect)
	       do 13 i=ix1,ix2,nxp
	       do 130 l=1,lamb2
	       do 130 j=jy1,jy2,nyp
 130	       write(81,9)zwg(ks)
  13	       write(81,9)zwg(ks)
	       close(81)
	    enddo
	    do jjc=1,nsez14
	       jjcd=jjc/10
	       jjcu=mod(jjc,10)
	       open(81,file='grid3_2D_slice'//char(48+jjcd)
     %//char(48+jjcu)//'_kplot.xyz')
	       js=jsez14(jjc)
	       jsy=jsym(js)
	       write(81,101)1,2*(nx-2)/nxp,(nzg-2)/nzp
	       do 210 k=2,nzg-1,nzp
	       do 211 i=ix2,ix1,-nxp
  211	       write(81,9)xc(i)*cos(yv(jsy))
	       do 212 i=ix1,ix2,nxp
  212	       write(81,9)xc(i)*cos(yv(js))
  210	       continue
	       do 220 k=2,nzg-1,nzp
	       do 221 i=ix2,ix1,-nxp
  221	       write(81,9)xc(i)*sin(yv(jsy))
	       do 222 i=ix1,ix2,nxp
  222	       write(81,9)xc(i)*sin(yv(js))
  220	       continue
	       do 230 k=2,nzg-1,nzp
	       do 231 i=ix2,ix1,-nxp
  231	       write(81,9)zcg(k)
	       do 232 i=ix1,ix2,nxp
  232	       write(81,9)zcg(k)
  230	       continue
	       close(81)
	    enddo
	 endif
      endif
c
      iplant4p=iplant4
      iplant4c=iplant4p/100
      iplant4p=mod(iplant4p,100)
      iplant4d=iplant4p/10
      iplant4p=mod(iplant4p,10)
c
      do iic=1,nsez1
	 is=isez1(iic)
	 if(myrank.eq.0) then
	    iicd=iic/10
	    iicu=mod(iic,10)
	    open(82,file='results1_2D_'//char(48+iplant4c)
     %//char(48+iplant4d)//char(48+iplant4p)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'_kplot1.q')
	    write(82,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	    write(82,102)ntime,i1,i1,i1
!
! Turbulent kinetic energy
!
	    do 40 k=kz1,kz2,nzp
	    do 41 l=1,lamb2
	    do 41 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)
     %0.5*(vrmsplot(iic,j,k)+urmsplot(iic,j,k)+
     %wrmsplot(iic,j,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
 41	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)
     %0.5*(vrmsplot(iic,jy1,k)+urmsplot(iic,jy1,k)+
     %wrmsplot(iic,jy1,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
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
		     write(82,9)recvar(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 43		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(82,9)recvar(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 42		  continue
	       enddo
	    endif
c	    write(82,*)	  !!!!!!
c	    write(82,103)
c    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
!
! Time-averaged azimuthal velocity
!
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)vplot(iic,j,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 51	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)vplot(iic,jy1,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
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
		     write(82,9)recvar(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 53		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(82,9)recvar(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 52		  continue
	       enddo
	    endif
!
! Time-averaged radial velocity
!
	    do 60 k=kz1,kz2,nzp
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)uplot(iic,j,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 61	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)uplot(iic,jy1,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
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
		     write(82,9)recvar(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 63		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(82,9)recvar(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 62		  continue
	       enddo
	    endif
!
! Time-averaged axial velocity
!
	    do 70 k=kz1,kz2,nzp
	    do 71 l=1,lamb2
	    do 71 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)wplot(iic,j,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 71	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)wplot(iic,jy1,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
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
		     write(82,9)recvar(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 73		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(82,9)recvar(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 72		  continue
	       enddo
	    endif
!
! Time-averaged static pressure
!
	    do 80 k=kz1,kz2,nzp
	    do 81 l=1,lamb2
	    do 81 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,8)pplot(iic,j,k)    !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
 81	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,8)pplot(iic,jy1,k)   !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
 80	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,6,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 82 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 83 lr=1,lamb2
		  do 83 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(82,8)recvar(jr,kr)	 !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
 83		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(82,8)recvar(jy1,kr)	 !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
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
     %(0.5*(vrmsplot(iic,jy1:jy2,kz1:kz2)+
     %urmsplot(iic,jy1:jy2,kz1:kz2)+
     %wrmsplot(iic,jy1:jy2,kz1:kz2)),1*(ny-2)*(nz-2),MTYPE,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(wplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(pplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
      do kkc=1,nsez2max
	 ks=ksez2(kkc)
	 kk=isezindx2+(kkc-1)
	 kkd=kk/10
	 kku=mod(kk,10)
	 open(82,file='results2_2D_'//char(48+iplant4c)
     %//char(48+iplant4d)//char(48+iplant4p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_kplot1.q')
	 write(82,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	 write(82,102)ntime,i1,i1,i1
!
! Turbulent kinetic energy
!
	 do 140 i=ix1,ix2,nxp
	 do 141 l=1,lamb2
	 do 141 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)
     %0.5*(vrmsplot2(i,j,kkc)+urmsplot2(i,j,kkc)+
     %wrmsplot2(i,j,kkc))
c	 else
c	    write(82,9)1.e+03
c	 endif
 141	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)
     %0.5*(vrmsplot2(i,jy1,kkc)+urmsplot2(i,jy1,kkc)+
     %wrmsplot2(i,jy1,kkc))
c	 else
c	    write(82,9)1.e+03
c	 endif
 140	 continue
c	 write(82,*)   !!!!!!
c	 write(82,103)
c    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
!
! Time-averaged azimuthal velocity
!
	 do 150 i=ix1,ix2,nxp
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)vplot2(i,j,kkc)
c	 else
c	    write(82,9)1.e+03
c	 endif
 151	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)vplot2(i,jy1,kkc)
c	 else
c	    write(82,9)1.e+03
c	 endif
 150	 continue
!
! Time-averaged radial velocity
!
	 do 160 i=ix1,ix2,nxp
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)uplot2(i,j,kkc)
c	 else
c	    write(82,9)1.e+03
c	 endif
 161	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)uplot2(i,jy1,kkc)
c	 else
c	    write(82,9)1.e+03
c	 endif
 160	 continue
!
! Time-averaged axial velocity
!
	 do 170 i=ix1,ix2,nxp
	 do 171 l=1,lamb2
	 do 171 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)wplot2(i,j,kkc)
c	 else
c	    write(82,9)1.e+03
c	 endif
 171	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)wplot2(i,jy1,kkc)
c	 else
c	    write(82,9)1.e+03
c	 endif
 170	 continue
!
! Time-averaged static pressure
!
	 do 180 i=ix1,ix2,nxp
	 do 181 l=1,lamb2
	 do 181 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,8)pplot2(i,j,kkc)	  !*ros
c	 else
c	    write(82,8)1.e+03
c	 endif
 181	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,8)pplot2(i,jy1,kkc)   !*ros
c	 else
c	    write(82,8)1.e+03
c	 endif
 180	 continue
	 close(82)
      enddo
c
      do jjc=1,nsez14
	 js=jsez14(jjc)
	 jsy=jsym(js)
	 if(myrank.eq.0) then
	    jjcd=jjc/10
	    jjcu=mod(jjc,10)
	    open(82,file='results3_2D_'//char(48+iplant4c)
     %//char(48+iplant4d)//char(48+iplant4p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_kplot1.q')
	    write(82,101)1,2*(nx-2)/nxp,(nzg-2)/nzp
	    write(82,102)ntime,i1,i1,i1
!
! Turbulent kinetic energy
!
	    do 240 k=kz1,kz2,nzp
	    do 241 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k).le.0) then
	       write(82,9)
     %0.5*(vrmsplot14(i,nsez14+jjc,k)+urmsplot14(i,nsez14+jjc,k)+
     %wrmsplot14(i,nsez14+jjc,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
 241	    continue
	    do 242 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k).le.0) then
	       write(82,9)
     %0.5*(vrmsplot14(i,jjc,k)+urmsplot14(i,jjc,k)+
     %wrmsplot14(i,jjc,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
 242	    continue
 240	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
c		  kk1=kz1+jjp*(nz-2)
c		  kk2=kz2+jjp*(nz-2)
c		  CALL MPI_RECV
c    %(flagvor2(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,7,
c    %MPI_COMM_EDDY,STATUS,IERR)
c		  CALL MPI_RECV
c    %(flagvor2sym(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,8,
c    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 243 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 244 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg).le.0) then
		     write(82,9)recvar2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 244		  continue
		  do 245 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg).le.0) then
		     write(82,9)recvar2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 245		  continue
 243		  continue
	       enddo
	    endif
!
! Time-averaged azimuthal velocity
!
	    do 250 k=kz1,kz2,nzp
	    do 251 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,9)vplot14(i,nsez14+jjc,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 251	    continue
	    do 252 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,9)vplot14(i,jjc,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 252	    continue
 250	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 253 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 254 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(82,9)recvar2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 254		  continue
		  do 255 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg,:)).le.0) then
		     write(82,9)recvar2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 255		  continue
 253		  continue
	       enddo
	    endif
!
! Time-averaged radial velocity
!
	    do 260 k=kz1,kz2,nzp
	    do 261 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,9)uplot14(i,nsez14+jjc,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 261	    continue
	    do 262 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,9)uplot14(i,jjc,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 262	    continue
 260	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 263 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 264 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(82,9)recvar2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 264		  continue
		  do 265 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg,:)).le.0) then
		     write(82,9)recvar2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 265		  continue
 263		  continue
	       enddo
	    endif
!
! Time-averaged axial velocity
!
	    do 270 k=kz1,kz2,nzp
	    do 271 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,9)wplot14(i,nsez14+jjc,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 271	    continue
	    do 272 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,9)wplot14(i,jjc,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 272	    continue
 270	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,16,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 273 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 274 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(82,9)recvar2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 274		  continue
		  do 275 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg,:)).le.0) then
		     write(82,9)recvar2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 275		  continue
 273		  continue
	       enddo
	    endif
!
! Time-averaged static pressure
!
	    do 280 k=kz1,kz2,nzp
	    do 281 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,8)pplot14(i,nsez14+jjc,k)    !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
 281	    continue
	    do 282 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,8)pplot14(i,jjc,k)    !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
 282	    continue
 280	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,17,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,18,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 283 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 284 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(82,9)recvar2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 284		  continue
		  do 285 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg,:)).le.0) then
		     write(82,9)recvar2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 285		  continue
 283		  continue
	       enddo
	    endif
	    close(82)
	 else
c	    CALL MPI_SEND
c    %(flagvo(ix1:ix2,js,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,7,
c    %MPI_COMM_EDDY,STATUS,IERR)
c	    CALL MPI_SEND
c    %(flagvo(ix1:ix2,jsy,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,8,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(0.5*(vrmsplot14(ix1:ix2,jjc,kz1:kz2)+
     %urmsplot14(ix1:ix2,jjc,kz1:kz2)+
     %wrmsplot14(ix1:ix2,jjc,kz1:kz2)),
     %(nx-2)*1*(nz-2),MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(0.5*(vrmsplot14(ix1:ix2,nsez14+jjc,kz1:kz2)+
     %urmsplot14(ix1:ix2,nsez14+jjc,kz1:kz2)+
     %wrmsplot14(ix1:ix2,nsez14+jjc,kz1:kz2)),
     %(nx-2)*1*(nz-2),MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(wplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(wplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,16,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(pplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,17,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(pplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,18,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine kplot2(iplant5c,iplant5d,iplant5p,
     %nx,ny,nz,nzg,nbd,ntime,flaguo,flagvo,flagwo)
c
      use sections
      use mediecalc1
      use mediecalc2
      use mediecalc3
      include'common.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer iplant5p,iplant5d,iplant5c,nx,ny,nz,nzg,nbd,ntime
      integer flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),
     %flagwo(nx,ny,nz,nbd)
c
c Local variables
      integer iic,iicd,iicu,is,jjc,jjcd,jjcu,js,kkc,ks,kk,kkd,kku,
     %i,j,k,l,ir,jr,kr,lr,jjp,i1,jsy
      integer STATUS(MPI_STATUS_SIZE)
!     integer kk1,kk2,krg,flaguor(ny,nzg,nbd),flagvor2(nx,nzg,nbd),
!    %flagvor2sym(nx,nzg,nbd)
      real recvar1(ny,nz),recvar2(ny,nz)
      real recvar3(nx,nz),recvar4(nx,nz)
      real recvar3sym(nx,nz),recvar4sym(nx,nz)
c
c Parameters
      integer lamb2,nxp,nyp,nzp
      parameter (lamb2=1)
      parameter (nyp=1,nxp=1,nzp=1)
c
   8  format(10f25.8)
   9  format(10f25.8)
 101  format(3i8)
 102  format(i8,3i2)
 103  format(i2)
c
      iplant5p=iplant5
      iplant5c=iplant5p/100
      iplant5p=mod(iplant5p,100)
      iplant5d=iplant5p/10
      iplant5p=mod(iplant5p,10)
c
      i1=1
      do iic=1,nsez1
	 is=isez1(iic)
	 if(myrank.eq.0) then
	    iicd=iic/10
	    iicu=mod(iic,10)
	    open(83,file='results1_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'_kplot2.q')
	    open(84,file='results1_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'_kplot3.q')
	    write(83,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	    write(84,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	    write(83,102)ntime,i1,i1,i1
	    write(84,102)ntime,i1,i1,i1
!
! Time-averaged velocity magnitude and stagnation pressure
!
	    do 40 k=kz1,kz2,nzp
	    do 41 l=1,lamb2
	    do 41 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)vmodplot(iic,j,k)
	       write(84,8)
     %(pplot(iic,j,k)+0.5*(vmodplot(iic,j,k)**2.))     !*ros
c	    else
c	       write(83,9)1.e+03
c	       write(84,8)1.e+03
c	    endif
 41	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)vmodplot(iic,jy1,k)
	       write(84,8)
     %(pplot(iic,jy1,k)+0.5*(vmodplot(iic,jy1,k)**2.))	   !*ros
c	    else
c	       write(83,9)1.e+03
c	       write(84,8)1.e+03
c	    endif
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
     %(recvar1(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 42 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 43 lr=1,lamb2
		  do 43 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,9)recvar1(jr,kr)
		     write(84,8)
     %(recvar2(jr,kr)+0.5*(recvar1(jr,kr)**2.))	     !*ros
c		  else
c		     write(83,9)1.e+03
c		     write(84,8)1.e+03
c		  endif
 43		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)recvar1(jy1,kr)
		     write(84,8)
     %(recvar2(jy1,kr)+0.5*(recvar1(jy1,kr)**2.))     !*ros
c		  else
c		     write(83,9)1.e+03
c		     write(84,8)1.e+03
c		  endif
 42		  continue
	       enddo
	    endif
c
!	    write(83,*)	  !!!!!!
!	    write(83,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
!	    write(84,*)	  !!!!!!
!	    write(84,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
!
! VRMS and UVPLOT
!
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)sqrt(abs(vrmsplot(iic,j,k)))
	       write(84,9)uvplot(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 51	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)sqrt(abs(vrmsplot(iic,jy1,k)))
	       write(84,9)uvplot(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 50	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar1(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 52 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 53 lr=1,lamb2
		  do 53 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar1(jr,kr)))
		     write(84,9)recvar2(jr,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 53		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar1(jy1,kr)))
		     write(84,9)recvar2(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 52		  continue
	       enddo
	    endif
!
! URMS and UWPLOT
!
	    do 60 k=kz1,kz2,nzp
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)sqrt(abs(urmsplot(iic,j,k)))
	       write(84,9)uwplot(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 61	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)sqrt(abs(urmsplot(iic,jy1,k)))
	       write(84,9)uwplot(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 60	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar1(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,6,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,7,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 62 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 63 lr=1,lamb2
		  do 63 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar1(jr,kr)))
		     write(84,9)recvar2(jr,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 63		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar1(jy1,kr)))
		     write(84,9)recvar2(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 62		  continue
	       enddo
	    endif
!
! WRMS and VWPLOT
!
	    do 70 k=kz1,kz2,nzp
	    do 71 l=1,lamb2
	    do 71 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)sqrt(abs(wrmsplot(iic,j,k)))
	       write(84,9)vwplot(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 71	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)sqrt(abs(wrmsplot(iic,jy1,k)))
	       write(84,9)vwplot(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 70	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar1(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 72 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 73 lr=1,lamb2
		  do 73 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar1(jr,kr)))
		     write(84,9)recvar2(jr,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 73		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar1(jy1,kr)))
		     write(84,9)recvar2(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 72		  continue
	       enddo
	    endif
!
! Static pressure root mean squares and time-averaged ratio
! eddy viscosity / fluid viscosity
!
	    do 80 k=kz1,kz2,nzp
	    do 81 l=1,lamb2
	    do 81 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,8)sqrt(abs(prmsplot(iic,j,k)))	  !*ros
	       write(84,9)tvplot(iic,j,k)
c	    else
c	       write(83,8)1.e+03
c	       write(84,9)1.e+03
c	    endif
 81	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,8)sqrt(abs(prmsplot(iic,jy1,k)))    !*ros
	       write(84,9)tvplot(iic,jy1,k)
c	    else
c	       write(83,8)1.e+03
c	       write(84,9)1.e+03
c	    endif
 80	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar1(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,jjp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 82 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 83 lr=1,lamb2
		  do 83 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,krg,:)).le.0) then
		     write(83,8)sqrt(abs(recvar1(jr,kr)))    !*ros
		     write(84,9)recvar2(jr,kr)
c		  else
c		     write(83,8)1.e+03
c		     write(84,9)1.e+03
c		  endif
 83		  continue
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,8)sqrt(abs(recvar1(jy1,kr)))    !*ros
		     write(84,9)recvar2(jy1,kr)
c		  else
c		     write(83,8)1.e+03
c		     write(84,9)1.e+03
c		  endif
 82		  continue
	       enddo
	    endif
c	     write(84,*)   !!!!!!
c	     write(84,103)
c     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
	    close(83)
	    close(84)
	 else
c	    CALL MPI_SEND
c    %(flaguo(is,jy1:jy2,kz1:kz2,:),
c    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,0,1,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vmodplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(pplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vrmsplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uvplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(urmsplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uwplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(wrmsplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vwplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(prmsplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(tvplot(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
      do kkc=1,nsez2max
	 ks=ksez2(kkc)
	 kk=isezindx2+(kkc-1)
	 kkd=kk/10
	 kku=mod(kk,10)
	 open(83,file='results2_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_kplot2.q')
	 open(84,file='results2_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_kplot3.q')
	 write(83,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	 write(84,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	 write(83,102)ntime,i1,i1,i1
	 write(84,102)ntime,i1,i1,i1
!
! Time-averaged velocity magnitude and stagnation pressure
!
	 do 140 i=ix1,ix2,nxp
	 do 141 l=1,lamb2
	 do 141 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)vmodplot2(i,j,kkc)
	    write(84,8)
     %(pplot2(i,j,kkc)+0.5*(vmodplot2(i,j,kkc)**2.))	 !*ros
c	 else
c	    write(83,9)1.e+03
c	    write(84,8)1.e+03
c	 endif
 141	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)vmodplot2(i,jy1,kkc)
	    write(84,8)
     %(pplot2(i,jy1,kkc)+0.5*(vmodplot2(i,jy1,kkc)**2.))     !*ros
c	 else
c	    write(83,9)1.e+03
c	    write(84,8)1.e+03
c	 endif
 140	 continue
!	 write(83,*)   !!!!!!
!	 write(83,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
!	 write(84,*)   !!!!!!
!	 write(84,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
!
! VRMS and UVPLOT
!
	 do 150 i=ix1,ix2,nxp
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)sqrt(abs(vrmsplot2(i,j,kkc)))
	    write(84,9)uvplot2(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 151	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)sqrt(abs(vrmsplot2(i,jy1,kkc)))
	    write(84,9)uvplot2(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 150	 continue
!
! URMS and UWPLOT
!
	 do 160 i=ix1,ix2,nxp
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)sqrt(abs(urmsplot2(i,j,kkc)))
	    write(84,9)uwplot2(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 161	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)sqrt(abs(urmsplot2(i,jy1,kkc)))
	    write(84,9)uwplot2(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 160	 continue
!
! WRMS and VWPLOT
!
	 do 170 i=ix1,ix2,nxp
	 do 171 l=1,lamb2
	 do 171 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)sqrt(abs(wrmsplot2(i,j,kkc)))
	    write(84,9)vwplot2(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 171	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)sqrt(abs(wrmsplot2(i,jy1,kkc)))
	    write(84,9)vwplot2(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	    write(84,9)1.e+03
c	 endif
 170	 continue
!
! Static pressure root mean squares and time-averaged ratio
! eddy viscosity / fluid viscosity
!
	 do 180 i=ix1,ix2,nxp
	 do 181 l=1,lamb2
	 do 181 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,8)sqrt(abs(prmsplot2(i,j,kkc)))	!*ros
	    write(84,9)tvplot2(i,j,kkc)
c	 else
c	    write(83,8)1.e+03
c	    write(84,9)1.e+03
c	 endif
 181	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,8)sqrt(abs(prmsplot2(i,jy1,kkc)))	  !*ros
	    write(84,9)tvplot2(i,jy1,kkc)
c	 else
c	    write(83,8)1.e+03
c	    write(84,9)1.e+03
c	 endif
 180	 continue
c	  write(84,*)	!!!!!!
c	  write(84,103)
c     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
	 close(83)
	 close(84)
      enddo
c
      do jjc=1,nsez14
	 js=jsez14(jjc)
	 jsy=jsym(js)
	 if(myrank.eq.0) then
	    jjcd=jjc/10
	    jjcu=mod(jjc,10)
	    open(83,file='results3_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_kplot2.q')
	    open(84,file='results3_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_kplot3.q')
	    write(83,101)1,2*(nx-2)/nxp,(nzg-2)/nzp
	    write(84,101)1,2*(nx-2)/nxp,(nzg-2)/nzp
	    write(83,102)ntime,i1,i1,i1
	    write(84,102)ntime,i1,i1,i1
!
! Time-averaged velocity magnitude and stagnation pressure
!
	    do 240 k=kz1,kz2,nzp
	    do 241 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(83,9)vmodplot14(i,nsez14+jjc,k)
	       write(84,8)
     %(pplot14(i,nsez14+jjc,k)+0.5*(vmodplot14(i,nsez14+jjc,k)**2.))	!*ros
c	    else
c	       write(83,9)1.e+03
c	       write(84,8)1.e+03
c	    endif
 241	    continue
	    do 242 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(83,9)vmodplot14(i,jjc,k)
	       write(84,8)
     %(pplot14(i,jjc,k)+0.5*(vmodplot14(i,jjc,k)**2.))	  !*ros
c	    else
c	       write(83,9)1.e+03
c	       write(84,8)1.e+03
c	    endif
 242	    continue
 240	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
c		  kk1=kz1+jjp*(nz-2)
c		  kk2=kz2+jjp*(nz-2)
c		  CALL MPI_RECV
c    %(flagvor2(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,12,
c    %MPI_COMM_EDDY,STATUS,IERR)
c		  CALL MPI_RECV
c    %(flagvor2sym(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,13,
c    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar3(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar3sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar4(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,16,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar4sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,17,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 243 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 244 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(83,9)recvar3sym(ir,kr)
		     write(84,8)
     %(recvar4sym(ir,kr)+0.5*(recvar3sym(ir,kr)**2.))	 !*ros
c		  else
c		     write(83,9)1.e+03
c		     write(84,8)1.e+03
c		  endif
 244		  continue
		  do 245 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg,:)).le.0) then
		     write(83,9)recvar3(ir,kr)
		     write(84,8)
     %(recvar4(ir,kr)+0.5*(recvar3(ir,kr)**2.))	   !*ros
c		  else
c		     write(83,9)1.e+03
c		     write(84,8)1.e+03
c		  endif
 245		  continue
 243		  continue
	       enddo
	    endif
!	    write(84,*)(((i1,j=1,1),i=1,2*(nx-2)/nxp),k=1,(nzg-2)/nzp)	 !!!!!!
!	    write(84,103)(((i1,j=1,1),i=1,2*(nx-2)/nxp),k=1,(nzg-2)/nzp)
!
! VRMS and UVPLOT
!
	    do 250 k=kz1,kz2,nzp
	    do 251 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(83,9)sqrt(abs(vrmsplot14(i,nsez14+jjc,k)))
	       write(84,9)uvplot14(i,nsez14+jjc,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 251	    continue
	    do 252 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(83,9)sqrt(abs(vrmsplot14(i,jjc,k)))
	       write(84,9)uvplot14(i,jjc,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 252	    continue
 250	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar3(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,18,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar3sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,19,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar4(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,20,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar4sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,21,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 253 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 254 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar3sym(ir,kr)))
		     write(84,9)recvar4sym(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 254		  continue
		  do 255 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar3(ir,kr)))
		     write(84,9)recvar4(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 255		  continue
 253		  continue
	       enddo
	    endif
!
! URMS and UWPLOT
!
	    do 260 k=kz1,kz2,nzp
	    do 261 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(83,9)sqrt(abs(urmsplot14(i,nsez14+jjc,k)))
	       write(84,9)uwplot14(i,nsez14+jjc,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 261	    continue
	    do 262 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(83,9)sqrt(abs(urmsplot14(i,jjc,k)))
	       write(84,9)uwplot14(i,jjc,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 262	    continue
 260	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar3(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,22,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar3sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,23,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar4(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,24,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar4sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,25,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 263 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 264 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar3sym(ir,kr)))
		     write(84,9)recvar4sym(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 264		  continue
		  do 265 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar3(ir,kr)))
		     write(84,9)recvar4(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 265		  continue
 263		  continue
	       enddo
	    endif
!
! WRMS and VWPLOT
!
	    do 270 k=kz1,kz2,nzp
	    do 271 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(83,9)sqrt(abs(wrmsplot14(i,nsez14+jjc,k)))
	       write(84,9)vwplot14(i,nsez14+jjc,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 271	    continue
	    do 272 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(83,9)sqrt(abs(wrmsplot14(i,jjc,k)))
	       write(84,9)vwplot14(i,jjc,k)
c	    else
c	       write(83,9)1.e+03
c	       write(84,9)1.e+03
c	    endif
 272	    continue
 270	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar3(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,26,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar3sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,27,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar4(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,28,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar4sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,29,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 273 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 274 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar3sym(ir,kr)))
		     write(84,9)recvar4sym(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 274		  continue
		  do 275 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg,:)).le.0) then
		     write(83,9)sqrt(abs(recvar3(ir,kr)))
		     write(84,9)recvar4(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		     write(84,9)1.e+03
c		  endif
 275		  continue
 273		  continue
	       enddo
	    endif
!
! Static pressure root mean squares anf time-averaged ratio
! eddy viscosity / fluid viscosity
!
	    do 280 k=kz1,kz2,nzp
	    do 281 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(83,8)sqrt(abs(prmsplot14(i,nsez14+jjc,k))) !*ros
	       write(84,9)tvplot14(i,nsez14+jjc,k)
c	    else
c	       write(83,8)1.e+03
c	       write(84,9)1.e+03
c	    endif
 281	    continue
	    do 282 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(83,8)sqrt(abs(prmsplot14(i,jjc,k))) !*ros
	       write(84,9)tvplot14(i,jjc,k)
c	    else
c	       write(83,8)1.e+03
c	       write(84,9)1.e+03
c	    endif
 282	    continue
 280	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar3(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,30,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar3sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,31,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar4(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,32,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar4sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,33,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 283 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 284 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(83,8)sqrt(abs(recvar3sym(ir,kr)))	!*ros
		     write(84,9)recvar4sym(ir,kr)
c		  else
c		     write(83,8)1.e+03
c		     write(84,9)1.e+03
c		  endif
 284		  continue
		  do 285 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,krg,:)).le.0) then
		     write(83,8)sqrt(abs(recvar3(ir,kr)))    !*ros
		     write(84,9)recvar4(ir,kr)
c		  else
c		     write(83,8)1.e+03
c		     write(84,9)1.e+03
c		  endif
 285		  continue
 283		  continue
	       enddo
	    endif
c	     write(84,*)(((i1,j=1,1),i=1,2*(nx-2)/nxp),k=1,(nzg-2)/nzp)	  !!!!!!
c	     write(84,103)(((i1,j=1,1),i=1,2*(nx-2)/nxp),k=1,(nzg-2)/nzp)
	    close(83)
	    close(84)
	 else
c	    CALL MPI_SEND
c    %(flagvo(ix1:ix2,js,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,12,
c    %MPI_COMM_EDDY,STATUS,IERR)
c	    CALL MPI_SEND
c    %(flagvo(ix1:ix2,jsy,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,13,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vmodplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vmodplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(pplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,16,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(pplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,17,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vrmsplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,18,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vrmsplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,19,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uvplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,20,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uvplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,21,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(urmsplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,22,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(urmsplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,23,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uwplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,24,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(uwplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,25,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(wrmsplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,26,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(wrmsplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,27,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vwplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,28,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vwplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,29,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(prmsplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,30,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(prmsplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,31,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(tvplot14(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,32,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(tvplot14(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,33,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine plotavgvort(iplant5c,iplant5d,iplant5p,
     %nx,ny,nz,nzg,nbd,ntime,flaguo,flagvo,flagwo)
c
!     use sections
      use vortavg
      include'common.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer iplant5p,iplant5d,iplant5c,nx,ny,nz,nzg,nbd,ntime
      integer flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),
     %flagwo(nx,ny,nz,nbd)
c
c Local variables
      integer iic,iicd,iicu,jjc,jjcd,jjcu,kkc,kk,
     %kkd,kku,i1,i,j,k,l,ir,jr,kr,lr,jjp
      integer STATUS(MPI_STATUS_SIZE)
!     integer is,js,jsy,ks,kk1,kk2,krg,flaguor(ny,nzg,nbd),
!    %flagvor2(nx,nzg,nbd),flagvor2sym(nx,nzg,nbd)
      real recvar(ny,nz),recvar2(nx,nz),recvar2sym(nx,nz)
c
c Parameters
      integer lamb2,nxp,nyp,nzp
      parameter (lamb2=1)
      parameter (nyp=1,nxp=1,nzp=1)
c
   8  format(10f25.8)
   9  format(10f25.8)
 101  format(3i8)
 102  format(i8,3i2)
 103  format(i2)
c
      i1=1
      do iic=1,nsez1
!	 is=isez1(iic)
	 if(myrank.eq.0) then
	    iicd=iic/10
	    iicu=mod(iic,10)
	    open(83,file='results1_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+iicd)//char(48+iicu)//'_avgvort.q')
	    write(83,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	    write(83,102)ntime,i1,i1,i1
!!!!!!		  write(83,*)
	    write(83,103)
     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
	    do 40 k=kz1,kz2,nzp
	    do 41 l=1,lamb2
	    do 41 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)voraz1med(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
 41	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)voraz1med(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)recvar(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
 42		  continue
	       enddo
	    endif
c
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)vorr1med(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
 51	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)vorr1med(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)recvar(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
 52		  continue
	       enddo
	    endif
c
	    do 60 k=kz1,kz2,nzp
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(83,9)vorz1med(iic,j,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
 61	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(83,9)vorz1med(iic,jy1,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,krg,:)).le.0) then
		     write(83,9)recvar(jy1,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
 62		  continue
	       enddo
	    endif
!!!!!!		  write(83,*)
	    write(83,103)
     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
	    close(83)
	 else
c	    CALL MPI_SEND
c    %(flaguo(is,jy1:jy2,kz1:kz2,:),
c    %1*(ny-2)*(nz-2)*nbd,MPI_INTEGER,0,1,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(voraz1med(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorr1med(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorz1med(iic,jy1:jy2,kz1:kz2),1*(ny-2)*(nz-2),MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
      do kkc=1,nsez2max
!	 ks=ksez2(kkc)
	 kk=isezindx2+(kkc-1)
	 kkd=kk/10
	 kku=mod(kk,10)
	 open(83,file='results2_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_avgvort.q')
	 write(83,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	 write(83,102)ntime,i1,i1,i1
!!!!!!	       write(83,*)
	 write(83,103)
     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
	 do 140 i=ix1,ix2,nxp
	 do 141 l=1,lamb2
	 do 141 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)voraz2med(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
 141	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)voraz2med(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
 140	 continue
	 do 150 i=ix1,ix2,nxp
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)vorr2med(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
 151	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)vorr2med(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
 150	 continue
	 do 160 i=ix1,ix2,nxp
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(83,9)vorz2med(i,j,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
 161	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(83,9)vorz2med(i,jy1,kkc)
c	 else
c	    write(83,9)1.e+03
c	 endif
 160	 continue
!!!!!!	       write(83,*)
	 write(83,103)
     %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
	 close(83)
      enddo
c
      do jjc=1,nsez14
!	 js=jsez14(jjc)
!	 jsy=jsym(js)
	 if(myrank.eq.0) then
	    jjcd=jjc/10
	    jjcu=mod(jjc,10)
	    open(83,file='results3_2D_'//char(48+iplant5c)
     %//char(48+iplant5d)//char(48+iplant5p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_avgvort.q')
	    write(83,101)1,2*(nx-2)/nxp,(nzg-2)/nzp
	    write(83,102)ntime,i1,i1,i1
!!!!!!		  write(83,*)(((i1,j=1,1),i=1,2*(nx-2)/nxp),k=1,(nzg-2)/nzp)
	    write(83,103)(((i1,j=1,1),i=1,2*(nx-2)/nxp),k=1,(nzg-2)/nzp)
	    do 240 k=kz1,kz2,nzp
	    do 241 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(83,9)voraz14med(i,nsez14+jjc,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
 241	    continue
	    do 242 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(83,9)voraz14med(i,jjc,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
 242	    continue
 240	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
c		  kk1=kz1+jjp*(nz-2)
c		  kk2=kz2+jjp*(nz-2)
c		  CALL MPI_RECV
c    %(flagvor2(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,5,
c    %MPI_COMM_EDDY,STATUS,IERR)
c		  CALL MPI_RECV
c    %(flagvor2sym(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,6,
c    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,7,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 243 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 244 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(83,9)recvar2sym(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
 244		  continue
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
	    do 251 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(83,9)vorr14med(i,nsez14+jjc,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
 251	    continue
	    do 252 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(83,9)vorr14med(i,jjc,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
 252	    continue
 250	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 253 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 254 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(83,9)recvar2sym(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
 254		  continue
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
	    do 261 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(83,9)vorz14med(i,nsez14+jjc,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
 261	    continue
	    do 262 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(83,9)vorz14med(i,jjc,k)
c	    else
c	       write(83,9)1.e+03
c	    endif
 262	    continue
 260	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  CALL MPI_RECV
     %(recvar2(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(recvar2sym(ix1:ix2,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,jjp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 263 kr=kz1,kz2,nzp
c		  krg=kr+jjp*(nz-2)
		  do 264 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,krg,:)).le.0) then
		     write(83,9)recvar2sym(ir,kr)
c		  else
c		     write(83,9)1.e+03
c		  endif
 264		  continue
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
!!!!!!		  write(83,*)(((i1,j=1,1),i=1,2*(nx-2)/nxp),k=1,(nzg-2)/nzp)
	    write(83,103)(((i1,j=1,1),i=1,2*(nx-2)/nxp),k=1,(nzg-2)/nzp)
	    close(83)
	 else
c	    CALL MPI_SEND
c    %(flagvo(ix1:ix2,js,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,5,
c    %MPI_COMM_EDDY,STATUS,IERR)
c	    CALL MPI_SEND
c    %(flagvo(ix1:ix2,jsy,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,6,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(voraz14med(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(voraz14med(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorr14med(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorr14med(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorz14med(ix1:ix2,jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vorz14med(ix1:ix2,nsez14+jjc,kz1:kz2),
     %(nx-2)*1*(nz-2),MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
c
      return
      end
c
C---- subroutine iompi_3dfields-------------------------------------
C
      subroutine iompi_3dfields(nx,ny,nz,ifields,uo,vo,wo,p,tv,
     %mbd,nbd,flaguo,flagvo,flagwo,flagpo,xc,yc,zc)
C
C     PURPOSE: Write to a file 3D fields using MPI subroutines
C
C-------------------------------------------------------------------
      include 'common.h'
      include 'averages.h'
      include 'fields.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,ifields,mbd,nbd
      integer flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),
     %flagwo(nx,ny,nz,nbd),flagpo(nx,ny,nz,nbd)
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),
     %p(nx,ny,nz),tv(nx,ny,nz)
      real xc(nx),yc(ny),zc(nz)
c
c.... Local arrays
      character*90 filename
      integer kz1g,kz2g,ibd
      integer fh,filemode
      integer newtype,ii,cost0,cost1,cost2
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)
      integer, dimension(:), allocatable :: nnzvect,nnzvectg
      ALLOCATE(nnzvect(mysize),nnzvectg(mysize))
c
c
      ifields=ifields+1
c
      if(ifields.eq.1) then

	 nnzvect = 0

	 nnx=limx2-limx1+1
	 if(mod(nnx,nnxp).eq.0) then
	    nnx=nnx/nnxp
	 else
	    nnx=nnx/nnxp+1
	 endif

	 nny=limy2-limy1+1
	 if(mod(nny,nnyp).eq.0) then
	    nny=nny/nnyp
	 else
	    nny=nny/nnyp+1
	 endif

	 kz1g=kz1+myrank*(nz-2)
	 kz2g=kz2+myrank*(nz-2)
	 if((limz1.le.kz2g).and.(limz2.ge.kz1g)) then
	    limz1l=max(kz1,limz1-myrank*(nz-2))
	    limz2l=min(kz2,limz2-myrank*(nz-2))
	    nnz=limz2l-limz1l+1
	    if(mod(nnz,nnzp).eq.0) then
	       nnz=nnz/nnzp
	    else
	       nnz=nnz/nnzp+1
	    endif
	 else
	    limz1l=nz
	    limz2l=1
	    nnz=0
	 endif

	 nnzvect(myrank+1)=nnz
	 CALL MPI_ALLREDUCE(nnzvect,nnzvectg,mysize,
     %MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
	 nnzg=sum(nnzvectg)
	 nnzprev=sum(nnzvectg(1:myrank))

	 if(isgs.eq.0) then
	    ivar=4
	 else
	    ivar=5
	 endif

	 if(myrank.eq.0) then

	    open(10,file='ngrid.dat',form='formatted')
	    write(10,*)nnx,nny,nnzg
	    close(10)

	    open(10,file='xgrid.bin',form='unformatted')
	    write(10)xc(limx1:limx2:nnxp)
	    close(10)

	    open(10,file='ygrid.bin',form='unformatted')
	    write(10)yc(limy1:limy2:nnyp)
	    close(10)

	 endif

	 filename='zgrid.bin'

	 filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

	 call mpi_file_open
     %(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	 call mpi_type_contiguous
     %(nnz,MTYPE,newtype,ierr)

	 call mpi_type_commit(newtype,ierr)

	 disp = 0
	 call mpi_file_set_view
     %(fh,disp,MTYPE,newtype,'native',mpi_info_null,ierr)

	 offset = nnzprev
	 call mpi_file_write_at_all(fh,offset,zc(limz1l:limz2l:nnzp),
     %1,newtype,status,ierr)

	 call mpi_type_free(newtype,ierr)

	 call mpi_file_close(fh,ierr)

	 if(itagvar.eq.1) then

	    filename='itagvar.bin'

	    filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

	    call mpi_file_open
     %(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	    if(mbd.ne.1) then

	       call mpi_type_contiguous
     %(nnx*nny*nnz,MPI_INTEGER,newtype,ierr)

	       call mpi_type_commit(newtype,ierr)

	       disp = 0
	       call mpi_file_set_view
     %(fh,disp,MPI_INTEGER,newtype,'native',mpi_info_null,ierr)

	       do ibd=1,mbd-1

!		   cost0=nnx*nny*nnzg*4*(ibd-1)
!		   cost1=nnx*nny*nnzg
!		   cost2=nnx*nny*nnzprev
!
!		   ii=1
!		   offset = cost0 + (ii-1)*cost1 + cost2
!		   call mpi_file_write_at_all(fh,offset,
!     %flaguo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp,
!     %ibd),1,newtype,status,ierr)
!		   ii=2
!		   offset = cost0 + (ii-1)*cost1 + cost2
!		   call mpi_file_write_at_all(fh,offset,
!     %flagvo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp,
!     %ibd),1,newtype,status,ierr)
!		   ii=3
!		   offset = cost0 + (ii-1)*cost1 + cost2
!		   call mpi_file_write_at_all(fh,offset,
!     %flagwo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp,
!     %ibd),1,newtype,status,ierr)
!		   ii=4
!		   offset = cost0 + (ii-1)*cost1 + cost2
!		   call mpi_file_write_at_all(fh,offset,
!     %flagpo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp,
!     %ibd),1,newtype,status,ierr)

		  cost0=nnx*nny*nnzg*(ibd-1)
		  cost2=nnx*nny*nnzprev

		  offset = cost0 + cost2
		  call mpi_file_write_at_all(fh,offset,
     %flagpo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp,
     %ibd),1,newtype,status,ierr)

	       enddo

	       call mpi_type_free(newtype,ierr)

	    endif

	    call mpi_file_close(fh,ierr)

	 endif

	 filename='fields.bin'

	 filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

	 call mpi_file_open
     %(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	 call mpi_file_close(fh,ierr)

      endif

      filename='fields.bin'

      filemode = IOR(MPI_MODE_APPEND,MPI_MODE_WRONLY)

      call mpi_file_open
     %(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

      call mpi_type_contiguous
     %(nnx*nny*nnz,mtype,newtype,ierr)

      call mpi_type_commit(newtype,ierr)

      disp = 0
      call mpi_file_set_view
     %(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

      cost0=nnx*nny*nnzg*ivar*(ifields-1)
      cost1=nnx*nny*nnzg
      cost2=nnx*nny*nnzprev

      ii=1
      offset = cost0 + (ii-1)*cost1 + cost2
      call mpi_file_write_at_all(fh,offset,
     %(uo(limx1-1:limx2-1:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp)+
     %uo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp))*0.5,
     %1,newtype,status,ierr)
      ii=2
      offset = cost0 + (ii-1)*cost1 + cost2
      call mpi_file_write_at_all(fh,offset,
     %(vo(limx1:limx2:nnxp,limy1-1:limy2-1:nnyp,limz1l:limz2l:nnzp)+
     %vo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp))*0.5,
     %1,newtype,status,ierr)
      ii=3
      offset = cost0 + (ii-1)*cost1 + cost2
      call mpi_file_write_at_all(fh,offset,
     %(wo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l-1:limz2l-1:nnzp)+
     %wo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp))*0.5,
     %1,newtype,status,ierr)
      ii=4
      offset = cost0 + (ii-1)*cost1 + cost2
      call mpi_file_write_at_all(fh,offset,
     %p(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp),
     %1,newtype,status,ierr)
      if(isgs.gt.0) then
	ii=5
	offset = cost0 + (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,
     %tv(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp),
     %1,newtype,status,ierr)
      endif

      call mpi_type_free(newtype,ierr)

      call mpi_file_close(fh,ierr)

      if((itagvar.eq.1).and.(mbd.le.nbd)) then

	 filename='itagvar.bin'

	 filemode = IOR(MPI_MODE_APPEND,MPI_MODE_WRONLY)

	 call mpi_file_open
     %(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	 call mpi_type_contiguous
     %(nnx*nny*nnz,MPI_INTEGER,newtype,ierr)

	 call mpi_type_commit(newtype,ierr)

	 disp = 0
	 call mpi_file_set_view
     %(fh,disp,MPI_INTEGER,newtype,'native',mpi_info_null,ierr)

	 do ibd=mbd,nbd

!	     cost0=nnx*nny*nnzg*4*
!     %((mbd-1)+(nbd-mbd+1)*(ifields-1)+(ibd-mbd))
!	     cost1=nnx*nny*nnzg
!	     cost2=nnx*nny*nnzprev
!
!	     ii=1
!	     offset = cost0 + (ii-1)*cost1 + cost2
!	     call mpi_file_write_at_all(fh,offset,
!     %flaguo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp,
!     %ibd),1,newtype,status,ierr)
!	     ii=2
!	     offset = cost0 + (ii-1)*cost1 + cost2
!	     call mpi_file_write_at_all(fh,offset,
!     %flagvo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp,
!     %ibd),1,newtype,status,ierr)
!	     ii=3
!	     offset = cost0 + (ii-1)*cost1 + cost2
!	     call mpi_file_write_at_all(fh,offset,
!     %flagwo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp,
!     %ibd),1,newtype,status,ierr)
!	     ii=4
!	     offset = cost0 + (ii-1)*cost1 + cost2
!	     call mpi_file_write_at_all(fh,offset,
!     %flagpo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp,
!     %ibd),1,newtype,status,ierr)

	    cost0=nnx*nny*nnzg*
     %((mbd-1)+(nbd-mbd+1)*(ifields-1)+(ibd-mbd))
	    cost2=nnx*nny*nnzprev

	    offset = cost0 + cost2
	    call mpi_file_write_at_all(fh,offset,
     %flagpo(limx1:limx2:nnxp,limy1:limy2:nnyp,limz1l:limz2l:nnzp,
     %ibd),1,newtype,status,ierr)

	 enddo

	 call mpi_type_free(newtype,ierr)

	 call mpi_file_close(fh,ierr)

      endif

      return

      end

C---- subroutine allocate_output------------------------------------
C
      subroutine allocate_output(nx,ny,nz)
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

      integer nx,ny,nz

      ALLOCATE(voraz1(nsezmax19,ny,nz),vorr1(nsezmax19,ny,nz),
     %vorz1(nsezmax19,ny,nz))
      ALLOCATE(voraz2(nx,ny,nsez2max),vorr2(nx,ny,nsez2max),
     %vorz2(nx,ny,nsez2max))
      ALLOCATE(voraz14(nx,2*nsez14,nz),vorr14(nx,2*nsez14,nz),
     %vorz14(nx,2*nsez14,nz))
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
      ALLOCATE(tempo4(nsez1,ny,nz),tempo5(nx,ny,nsez2max),
     %tempo14(nx,2*nsez14,nz))
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
      ALLOCATE(vplot14(nx,2*nsez14,nz),uplot14(nx,2*nsez14,nz),
     %wplot14(nx,2*nsez14,nz),pplot14(nx,2*nsez14,nz),
     %vrmsplot14(nx,2*nsez14,nz),urmsplot14(nx,2*nsez14,nz),
     %wrmsplot14(nx,2*nsez14,nz),prmsplot14(nx,2*nsez14,nz),
     %uvplot14(nx,2*nsez14,nz),vwplot14(nx,2*nsez14,nz),
     %uwplot14(nx,2*nsez14,nz),vmodplot14(nx,2*nsez14,nz),
     %tvplot14(nx,2*nsez14,nz))
      ALLOCATE(voraz1med(nsez1,ny,nz),vorr1med(nsez1,ny,nz),
     %vorz1med(nsez1,ny,nz),voraz2med(nx,ny,nsez2max),
     %vorr2med(nx,ny,nsez2max),vorz2med(nx,ny,nsez2max),
     %voraz14med(nx,2*nsez14,nz),vorr14med(nx,2*nsez14,nz),
     %vorz14med(nx,2*nsez14,nz))
      ALLOCATE(timepoints1(nprbmax2))
      ALLOCATE(avrpoints1(nprbmax2,4),rmspoints1(nprbmax2,4))

      return
      end

C---- subroutine iompi_3dscalar_mediey------------------------------
C
      subroutine iompi_3dscalar_mediey(filename,ny,io)
C
C     PURPOSE: Read from or write to a file a 3D real array using
C     MPI subroutines.
C
C-------------------------------------------------------------------
      use mediey
      use rmsy
      include 'common.h'
      include 'averages.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer ny,io
      character*(*) filename
c
c.... Local arrays
      integer fh,filemode
      integer newtype,ii,cost1,cost2
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      if(io==1) then

	filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

	call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %(nsez3*(ny-2)*nsez4max,mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=nsez3*(ny-2)*nsez4
	cost2=nsez3*(ny-2)*max(0,isezindx4-1)

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,tempo1(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,uavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,wavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vmodavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vorazavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vorravg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vorzavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,pavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,ptotavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vrms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,urms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,wrms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,prms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,uvs(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,uws(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,vws(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=18
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,tvavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      else

	filemode = MPI_MODE_RDONLY

	call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %(nsez3*(ny-2)*nsez4max,mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=nsez3*(ny-2)*nsez4
	cost2=nsez3*(ny-2)*max(0,isezindx4-1)

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,tempo1(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,uavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,wavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vmodavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vorazavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vorravg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vorzavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,pavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,ptotavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vrms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,urms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,wrms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,prms(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,uvs(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,uws(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,vws(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)
	ii=18
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,tvavg(:,jy1:jy2,1:nsez4max),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      endif

      return

      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C---- subroutine iompi_3dscalar_mediex------------------------------
C
      subroutine iompi_3dscalar_mediex(filename,nx,io)
C
C     PURPOSE: Read from or write to a file a 3D real array using
C     MPI subroutines.
C
C-------------------------------------------------------------------
      use mediex
      use rmsx
      include 'common.h'
      include 'averages.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,io
      character*(*) filename
c
c.... Local arrays
      integer fh,filemode
      integer newtype,ii,cost1,cost2
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      if(io==1) then

	filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

	call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %((nx-2)*nsez5*nsez6max,mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=(nx-2)*nsez5*nsez6
	cost2=(nx-2)*nsez5*max(0,isezindx6-1)

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,tempo2(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,ruavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rwavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvmodavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rpavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rptotavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvorazavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvorravg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvorzavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvrms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rurms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rwrms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rprms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,ruvs(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,ruws(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rvws(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=18
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,rtvavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      else

	filemode = MPI_MODE_RDONLY

	call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %((nx-2)*nsez5*nsez6max,mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=(nx-2)*nsez5*nsez6
	cost2=(nx-2)*nsez5*max(0,isezindx6-1)

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,tempo2(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,ruavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rwavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvmodavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rpavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rptotavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvorazavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvorravg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvorzavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvrms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rurms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rwrms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rprms(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,ruvs(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,ruws(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rvws(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)
	ii=18
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,rtvavg(ix1:ix2,:,1:nsez6max),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      endif

      return

      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C---- subroutine iompi_3dscalar_mediez------------------------------
C
      subroutine iompi_3dscalar_mediez(filename,nz,io)
C
C     PURPOSE: Read from or write to a file a 3D real array using
C     MPI subroutines.
C
C-------------------------------------------------------------------
      use mediez
      use rmsz
      include 'common.h'
      include 'averages.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nz,io
      character*(*) filename
c
c.... Local arrays
      integer fh,filemode
      integer newtype,ii,cost1,cost2
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      if(io==1) then

	filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

	call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %(nsez9*nsez91*(nz-2),mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=nsez9*nsez91*(nz-2)*mysize
	cost2=nsez9*nsez91*(nz-2)*myrank

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zuavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zwavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvmodavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zpavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zptotavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvorazavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvorravg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvorzavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvrms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zurms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zwrms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zprms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zuvs(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zuws(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,zvws(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,ztempo(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=18
	offset = (ii-1)*cost1 + cost2
	call mpi_file_write_at_all(fh,offset,ztvavg(:,:,kz1:kz2),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      else

	filemode = MPI_MODE_RDONLY

	call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

	call mpi_type_contiguous
     %(nsez9*nsez91*(nz-2),mtype,newtype,ierr)

	call mpi_type_commit(newtype,ierr)

	disp = 0
	call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

	cost1=nsez9*nsez91*(nz-2)*mysize
	cost2=nsez9*nsez91*(nz-2)*myrank

	ii=1
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=2
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zuavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=3
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zwavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=4
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvmodavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=5
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zpavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=6
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zptotavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=7
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvorazavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=8
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvorravg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=9
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvorzavg(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=10
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvrms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=11
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zurms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=12
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zwrms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=13
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zprms(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=14
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zuvs(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=15
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zuws(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=16
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,zvws(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=17
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,ztempo(:,:,kz1:kz2),1,newtype,status,ierr)
	ii=18
	offset = (ii-1)*cost1 + cost2
	call mpi_file_read_at_all(fh,offset,ztvavg(:,:,kz1:kz2),1,newtype,status,ierr)

	call mpi_type_free(newtype,ierr)

	call mpi_file_close(fh,ierr)

      endif

      return

      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C---- subroutine setup_indices -------------------A. Posa - Apr 2012------
c
      subroutine setup_indices(ny)
c
      include 'common.h'
      include 'averages.h'
c
      integer j,ny
c
      do 10 j=jy1+1,jy2-1
      jpv(j)=j+1
      jmv(j)=j-1
  10  continue
c
      jpv(jy1)=jy1+1
      jmv(jy1)=jy2
      jpv(jy2)=jy1
      jmv(jy2)=jy2-1
c
      jpv(1)=jy1
      jmv(1)=jy2-1
      jpv(ny)=jy1+1
      jmv(ny)=jy2
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Subroutine for the evaluation of the statistics on the cylindrical
c grid in Cartesian coordinates
c
      subroutine mediestaggnoavgcart(nx,ny,nz,yc,yv,uo,vo,wo,p,tv,dt,nbd,flagpo,time)
c
      use sections
      use mediey
      use mediex
      use mediez
      use vorticity
      use rmsy
      use rmsx
      use rmsz
      include'common.h'
      include'averages.h'
      include'mpif.h'
c
c Global variables
      integer nx,ny,nz,nbd,flagpo(nx,ny,nz,nbd)
      real dt,time
      real yc(ny),yv(ny)
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),
     %p(nx,ny,nz),tv(nx,ny,nz)
c
c Local variables
      integer iic,i,is,ip,im,jjc,j,js,jp,jm,jsy,jpsy,jmsy,
     %kkc,k,ks,kp,km,mm1,mm2,kg,nsez4maxr,nsez6maxr
      integer, dimension(:), allocatable :: ksez4r,ksez6r
      integer STATUS(MPI_STATUS_SIZE)
      real ustagg,vstagg,wstagg,pstagg,vmod,tvstagg,xstagg,ystagg
      real, dimension(:,:,:), allocatable :: vmodavgr,uavgr,vavgr,
     %wavgr,pavgr,ptotavgr,vorravgr,vorazavgr,vorzavgr,tvavgr,
     %rvmodavgr,ruavgr,rvavgr,rwavgr,rpavgr,rptotavgr,rvorravgr,
     %rvorazavgr,rvorzavgr,rtvavgr
      real zuavgr(nz),zvavgr(nz),zwavgr(nz),zpavgr(nz)
c
  10  format(4i6,10f25.8)
c
      tempo=tempo+dt
!
!     Time averages at some (x,z) positions along y
!
      do 101 iic=1,nsez3
      is=isez3(iic)
      ip=is+1
      if(is.eq.1) then
	do 1010 kkc=1,nsez4max
	ks=ksez4(kkc)
	kp=ks+1
	do 1011 j=jy1,jy2
	jm=jmv(j)
	jsy=jsym(j)
	jmsy=jsym(jm)
	vstagg=(vo(ip,j,ks)-vo(ip,jsy,ks)+
     %vo(ip,jm,ks)-vo(ip,jmsy,ks)+
     %vo(ip,j,kp)-vo(ip,jsy,kp)+
     %vo(ip,jm,kp)-vo(ip,jmsy,kp))*0.125
	ustagg=(uo(is,j,ks)+uo(is,j,kp))*0.5
	wstagg=(wo(ip,j,ks)+wo(ip,jsy,ks))*0.5
	xstagg=ustagg*cos(yc(j))-vstagg*sin(yc(j))
	ystagg=ustagg*sin(yc(j))+vstagg*cos(yc(j))
	vstagg=xstagg
	ustagg=ystagg
	vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
	vavg(iic,j,kkc)=vavg(iic,j,kkc)
     %+vstagg*dt
	uavg(iic,j,kkc)=uavg(iic,j,kkc)
     %+ustagg*dt
	wavg(iic,j,kkc)=wavg(iic,j,kkc)
     %+wstagg*dt
	vmodavg(iic,j,kkc)=vmodavg(iic,j,kkc)
     %+vmod*dt
	vorazavg(iic,j,kkc)=vorazavg(iic,j,kkc)
     %+(vorr(iic,j,kkc)*cos(yc(j))
     %-voraz(iic,j,kkc)*sin(yc(j)))*dt
	vorravg(iic,j,kkc)=vorravg(iic,j,kkc)
     %+(vorr(iic,j,kkc)*sin(yc(j))
     %+voraz(iic,j,kkc)*cos(yc(j)))*dt
	vorzavg(iic,j,kkc)=vorzavg(iic,j,kkc)
     %+vorz(iic,j,kkc)*dt
	vrms(iic,j,kkc)=
     %vrms(iic,j,kkc)+dt*(vstagg**2.)
	urms(iic,j,kkc)=
     %urms(iic,j,kkc)+dt*(ustagg**2.)
	wrms(iic,j,kkc)=
     %wrms(iic,j,kkc)+dt*(wstagg**2.)
	uvs(iic,j,kkc)=
     %uvs(iic,j,kkc)+dt*ustagg*vstagg
	uws(iic,j,kkc)=
     %uws(iic,j,kkc)+dt*ustagg*wstagg
	vws(iic,j,kkc)=
     %vws(iic,j,kkc)+dt*vstagg*wstagg
	if(all(flagpo(is:ip,j,ks:kp,:).le.0)) then
	   tempo1(iic,j,kkc)=tempo1(iic,j,kkc)+dt
	   pstagg=(p(ip,j,ks)+p(ip,j,kp)+
     %p(ip,jsy,ks)+p(ip,jsy,kp))*0.25
	   pavg(iic,j,kkc)=pavg(iic,j,kkc)
     %+pstagg*dt    !*ros
	   ptotavg(iic,j,kkc)=ptotavg(iic,j,kkc)
     %+(pstagg+0.5*(vmod**2.))*dt   !*ros
	   prms(iic,j,kkc)=
     %prms(iic,j,kkc)+dt*(pstagg**2.)	!*ros**2.
	endif
	tvstagg=(tv(ip,j,ks)+tv(ip,j,kp)+
     %tv(ip,jsy,ks)+tv(ip,jsy,kp))*0.25
	tvavg(iic,j,kkc)=tvavg(iic,j,kkc)
     %+tvstagg*dt
 1011	continue
 1010	continue
      else
	do 1012 kkc=1,nsez4max
	ks=ksez4(kkc)
	kp=ks+1
	do 1013 j=jy1,jy2
	jm=jmv(j)
	vstagg=(vo(is,j,ks)+vo(ip,j,ks)+
     %vo(is,jm,ks)+vo(ip,jm,ks)+
     %vo(is,j,kp)+vo(ip,j,kp)+
     %vo(is,jm,kp)+vo(ip,jm,kp))*0.125
	ustagg=(uo(is,j,ks)+uo(is,j,kp))*0.5
	wstagg=(wo(is,j,ks)+wo(ip,j,ks))*0.5
	xstagg=ustagg*cos(yc(j))-vstagg*sin(yc(j))
	ystagg=ustagg*sin(yc(j))+vstagg*cos(yc(j))
	vstagg=xstagg
	ustagg=ystagg
	vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
	vavg(iic,j,kkc)=vavg(iic,j,kkc)
     %+vstagg*dt
	uavg(iic,j,kkc)=uavg(iic,j,kkc)
     %+ustagg*dt
	wavg(iic,j,kkc)=wavg(iic,j,kkc)
     %+wstagg*dt
	vmodavg(iic,j,kkc)=vmodavg(iic,j,kkc)
     %+vmod*dt
	vorazavg(iic,j,kkc)=vorazavg(iic,j,kkc)
     %+(vorr(iic,j,kkc)*cos(yc(j))
     %-voraz(iic,j,kkc)*sin(yc(j)))*dt
	vorravg(iic,j,kkc)=vorravg(iic,j,kkc)
     %+(vorr(iic,j,kkc)*sin(yc(j))
     %+voraz(iic,j,kkc)*cos(yc(j)))*dt
	vorzavg(iic,j,kkc)=vorzavg(iic,j,kkc)
     %+vorz(iic,j,kkc)*dt
	vrms(iic,j,kkc)=
     %vrms(iic,j,kkc)+dt*(vstagg**2.)
	urms(iic,j,kkc)=
     %urms(iic,j,kkc)+dt*(ustagg**2.)
	wrms(iic,j,kkc)=
     %wrms(iic,j,kkc)+dt*(wstagg**2.)
	uvs(iic,j,kkc)=
     %uvs(iic,j,kkc)+dt*ustagg*vstagg
	uws(iic,j,kkc)=
     %uws(iic,j,kkc)+dt*ustagg*wstagg
	vws(iic,j,kkc)=
     %vws(iic,j,kkc)+dt*vstagg*wstagg
	if(all(flagpo(is:ip,j,ks:kp,:).le.0)) then
	   tempo1(iic,j,kkc)=tempo1(iic,j,kkc)+dt
	   pstagg=(p(is,j,ks)+p(ip,j,ks)+
     %p(is,j,kp)+p(ip,j,kp))*0.25
	   pavg(iic,j,kkc)=pavg(iic,j,kkc)
     %+pstagg*dt    !*ros
	   ptotavg(iic,j,kkc)=ptotavg(iic,j,kkc)
     %+(pstagg+0.5*(vmod**2.))*dt   !*ros
	   prms(iic,j,kkc)=
     %prms(iic,j,kkc)+dt*(pstagg**2.)	!*ros**2.
	endif
	tvstagg=(tv(is,j,ks)+tv(ip,j,ks)+
     %tv(is,j,kp)+tv(ip,j,kp))*0.25
	tvavg(iic,j,kkc)=tvavg(iic,j,kkc)
     %+tvstagg*dt
 1013	continue
 1012	continue
      endif
  101 continue
!
!     Time averages at some (y,z) positions along x
!
      do 104 jjc=1,nsez5
      js=jsez5(jjc)
      jp=jpv(js)
      do 104 kkc=1,nsez6max
      ks=ksez6(kkc)
      kp=ks+1
      do 103 i=ix1,ix2
      im=i-1
      vstagg=(vo(i,js,ks)+vo(i,js,kp))*0.5
      ustagg=(uo(i,js,ks)+uo(im,js,ks)+
     %uo(i,jp,ks)+uo(im,jp,ks)+
     %uo(i,js,kp)+uo(im,js,kp)+
     %uo(i,jp,kp)+uo(im,jp,kp))*0.125
      wstagg=(wo(i,js,ks)+wo(i,jp,ks))*0.5
      xstagg=ustagg*cos(yv(js))-vstagg*sin(yv(js))
      ystagg=ustagg*sin(yv(js))+vstagg*cos(yv(js))
      vstagg=xstagg
      ustagg=ystagg
      rvavg(i,jjc,kkc)=rvavg(i,jjc,kkc)
     %+vstagg*dt
      ruavg(i,jjc,kkc)=ruavg(i,jjc,kkc)
     %+ustagg*dt
      rwavg(i,jjc,kkc)=rwavg(i,jjc,kkc)
     %+wstagg*dt
      vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
      rvmodavg(i,jjc,kkc)=rvmodavg(i,jjc,kkc)
     %+vmod*dt
      rvorazavg(i,jjc,kkc)=rvorazavg(i,jjc,kkc)
     %+(rvorr(i,jjc,kkc)*cos(yv(js))
     %-rvoraz(i,jjc,kkc)*sin(yv(js)))*dt
      rvorravg(i,jjc,kkc)=rvorravg(i,jjc,kkc)
     %+(rvorr(i,jjc,kkc)*sin(yv(js))
     %+rvoraz(i,jjc,kkc)*cos(yv(js)))*dt
      rvorzavg(i,jjc,kkc)=rvorzavg(i,jjc,kkc)
     %+rvorz(i,jjc,kkc)*dt
      rvrms(i,jjc,kkc)=
     %rvrms(i,jjc,kkc)+dt*(vstagg**2.)
      rurms(i,jjc,kkc)=
     %rurms(i,jjc,kkc)+dt*(ustagg**2.)
      rwrms(i,jjc,kkc)=
     %rwrms(i,jjc,kkc)+dt*(wstagg**2.)
      ruvs(i,jjc,kkc)=
     %ruvs(i,jjc,kkc)+dt*ustagg*vstagg
      ruws(i,jjc,kkc)=
     %ruws(i,jjc,kkc)+dt*ustagg*wstagg
      rvws(i,jjc,kkc)=
     %rvws(i,jjc,kkc)+dt*vstagg*wstagg
      if(all(flagpo(i,js:jp,ks:kp,:).le.0)) then
	 tempo2(i,jjc,kkc)=tempo2(i,jjc,kkc)+dt
	 pstagg=(p(i,js,ks)+p(i,jp,ks)+
     %p(i,js,kp)+p(i,jp,kp))*0.25
	 rpavg(i,jjc,kkc)=rpavg(i,jjc,kkc)+
     %pstagg*dt	  !*ros
	 rptotavg(i,jjc,kkc)=rptotavg(i,jjc,kkc)
     %+(pstagg+0.5*(vmod**2.))*dt   !*ros
	 rprms(i,jjc,kkc)=
     %rprms(i,jjc,kkc)+dt*(pstagg**2.)	 !*ros**2.
      endif
      tvstagg=(tv(i,js,ks)+tv(i,jp,ks)+
     %tv(i,js,kp)+tv(i,jp,kp))*0.25
      rtvavg(i,jjc,kkc)=rtvavg(i,jjc,kkc)+
     %tvstagg*dt
  103 continue
  104 continue
!
!     Time averages at some (x,y) positions along z
!
      do 106 jjc=1,nsez91
      js=jsez91(jjc)
      jp=jpv(js)
      jsy=jsym(js)
      jpsy=jsym(jp)
      do 106 iic=1,nsez9
      is=isez9(iic)
      ip=is+1
      if(is.eq.1) then
	do 1060 k=kz1,kz2
	km=k-1
	vstagg=(vo(ip,js,k)-vo(ip,jsy,k))*0.5
	ustagg=(uo(is,js,k)+uo(is,jp,k))*0.5
	wstagg=(wo(ip,js,k)+wo(ip,jsy,k)+
     %wo(ip,jp,k)+wo(ip,jpsy,k)+wo(ip,js,km)+wo(ip,jsy,km)+
     %wo(ip,jp,km)+wo(ip,jpsy,km))*0.125
	xstagg=ustagg*cos(yv(js))-vstagg*sin(yv(js))
	ystagg=ustagg*sin(yv(js))+vstagg*cos(yv(js))
	vstagg=xstagg
	ustagg=ystagg
	zvavg(iic,jjc,k)=zvavg(iic,jjc,k)
     %+vstagg*dt
	zuavg(iic,jjc,k)=zuavg(iic,jjc,k)
     %+ustagg*dt
	zwavg(iic,jjc,k)=zwavg(iic,jjc,k)
     %+wstagg*dt
	vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
	zvmodavg(iic,jjc,k)=
     %zvmodavg(iic,jjc,k)+vmod*dt
	zvorazavg(iic,jjc,k)=zvorazavg(iic,jjc,k)
     %+(vorr1(iic,jjc,k)*cos(yv(js))
     %-voraz1(iic,jjc,k)*sin(yv(js)))*dt
	zvorravg(iic,jjc,k)=zvorravg(iic,jjc,k)
     %+(vorr1(iic,jjc,k)*sin(yv(js))
     %+voraz1(iic,jjc,k)*cos(yv(js)))*dt
	zvorzavg(iic,jjc,k)=zvorzavg(iic,jjc,k)
     %+vorz1(iic,jjc,k)*dt
	zvrms(iic,jjc,k)=zvrms(iic,jjc,k)
     %+dt*(vstagg**2.)
	zurms(iic,jjc,k)=zurms(iic,jjc,k)
     %+dt*(ustagg**2.)
	zwrms(iic,jjc,k)=zwrms(iic,jjc,k)
     %+dt*(wstagg**2.)
	zuvs(iic,jjc,k)=zuvs(iic,jjc,k)
     %+dt*ustagg*vstagg
	zuws(iic,jjc,k)=zuws(iic,jjc,k)
     %+dt*ustagg*wstagg
	zvws(iic,jjc,k)=zvws(iic,jjc,k)
     %+dt*vstagg*wstagg
	if(all(flagpo(is:ip,js:jp,k,:).le.0)) then
	   ztempo(iic,jjc,k)=ztempo(iic,jjc,k)+dt
	   pstagg=(p(ip,js,k)+p(ip,jsy,k)+
     %p(ip,jp,k)+p(ip,jpsy,k))*0.25
	   zpavg(iic,jjc,k)=zpavg(iic,jjc,k)+
     %pstagg*dt	  !*ros
	   zptotavg(iic,jjc,k)=zptotavg(iic,jjc,k)+
     %(pstagg+0.5*(vmod**2.))*dt   !*ros
	   zprms(iic,jjc,k)=zprms(iic,jjc,k)
     %+dt*(pstagg**2.)	 !*ros**2.
	endif
	tvstagg=(tv(ip,js,k)+tv(ip,jsy,k)+
     %tv(ip,jp,k)+tv(ip,jpsy,k))*0.25
	ztvavg(iic,jjc,k)=ztvavg(iic,jjc,k)+
     %tvstagg*dt
 1060	continue
      else
	do 1061 k=kz1,kz2
	km=k-1
	vstagg=(vo(is,js,k)+vo(ip,js,k))*0.5
	ustagg=(uo(is,js,k)+uo(is,jp,k))*0.5
	wstagg=(wo(is,js,k)+wo(ip,js,k)+
     %wo(is,jp,k)+wo(ip,jp,k)+wo(is,js,km)+wo(ip,js,km)+
     %wo(is,jp,km)+wo(ip,jp,km))*0.125
	xstagg=ustagg*cos(yv(js))-vstagg*sin(yv(js))
	ystagg=ustagg*sin(yv(js))+vstagg*cos(yv(js))
	vstagg=xstagg
	ustagg=ystagg
	zvavg(iic,jjc,k)=zvavg(iic,jjc,k)
     %+vstagg*dt
	zuavg(iic,jjc,k)=zuavg(iic,jjc,k)
     %+ustagg*dt
	zwavg(iic,jjc,k)=zwavg(iic,jjc,k)
     %+wstagg*dt
	vmod=
     %sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
	zvmodavg(iic,jjc,k)=
     %zvmodavg(iic,jjc,k)+vmod*dt
	zvorazavg(iic,jjc,k)=zvorazavg(iic,jjc,k)
     %+(vorr1(iic,jjc,k)*cos(yv(js))
     %-voraz1(iic,jjc,k)*sin(yv(js)))*dt
	zvorravg(iic,jjc,k)=zvorravg(iic,jjc,k)
     %+(vorr1(iic,jjc,k)*sin(yv(js))
     %+voraz1(iic,jjc,k)*cos(yv(js)))*dt
	zvorzavg(iic,jjc,k)=zvorzavg(iic,jjc,k)
     %+vorz1(iic,jjc,k)*dt
	zvrms(iic,jjc,k)=zvrms(iic,jjc,k)
     %+dt*(vstagg**2.)
	zurms(iic,jjc,k)=zurms(iic,jjc,k)
     %+dt*(ustagg**2.)
	zwrms(iic,jjc,k)=zwrms(iic,jjc,k)
     %+dt*(wstagg**2.)
	zuvs(iic,jjc,k)=zuvs(iic,jjc,k)
     %+dt*ustagg*vstagg
	zuws(iic,jjc,k)=zuws(iic,jjc,k)
     %+dt*ustagg*wstagg
	zvws(iic,jjc,k)=zvws(iic,jjc,k)
     %+dt*vstagg*wstagg
	if(all(flagpo(is:ip,js:jp,k,:).le.0)) then
	   ztempo(iic,jjc,k)=ztempo(iic,jjc,k)+dt
	   pstagg=(p(is,js,k)+p(ip,js,k)+
     %p(is,jp,k)+p(ip,jp,k))*0.25
	   zpavg(iic,jjc,k)=zpavg(iic,jjc,k)+
     %pstagg*dt	  !*ros
	   zptotavg(iic,jjc,k)=zptotavg(iic,jjc,k)+
     %(pstagg+0.5*(vmod**2.))*dt   !*ros
	   zprms(iic,jjc,k)=zprms(iic,jjc,k)
     %+dt*(pstagg**2.)	  !*ros**2.
	endif
	tvstagg=(tv(is,js,k)+tv(ip,js,k)+
     %tv(is,jp,k)+tv(ip,jp,k))*0.25
	ztvavg(iic,jjc,k)=ztvavg(iic,jjc,k)+
     %tvstagg*dt
 1061	continue
      endif
  106 continue
      if(tempo.ge.periodo4) then
	 iplant8=iplant8+1
	 do 200 kkc=1,nsez4max
	 do 200 iic=1,nsez3
	 do 200 j=jy1,jy2
	 vmodavg(iic,j,kkc)=vmodavg(iic,j,kkc)/
     %tempo
	 vavg(iic,j,kkc)=vavg(iic,j,kkc)/
     %tempo
	 vorazavg(iic,j,kkc)=vorazavg(iic,j,kkc)/
     %tempo
	 uavg(iic,j,kkc)=uavg(iic,j,kkc)/
     %tempo
	 vorravg(iic,j,kkc)=vorravg(iic,j,kkc)/
     %tempo
	 wavg(iic,j,kkc)=wavg(iic,j,kkc)/
     %tempo
	 vorzavg(iic,j,kkc)=vorzavg(iic,j,kkc)/
     %tempo
	 if(tempo1(iic,j,kkc).ne.0.) then
	    pavg(iic,j,kkc)=pavg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	    ptotavg(iic,j,kkc)=ptotavg(iic,j,kkc)/
     %tempo1(iic,j,kkc)
	 endif
	 tvavg(iic,j,kkc)=tvavg(iic,j,kkc)/
     %tempo
	 tvavg(iic,j,kkc)=tvavg(iic,j,kkc)/ru1
  200	 continue
	 if(myrank.eq.0) then
	    mm1=121
	    mm2=122
	    do 201 kkc=1,nsez4max
	    do 201 iic=1,nsez3
	    write(mm1,333)'i=',isez3(iic),'____k=',ksez4(kkc),time-tempo*0.5
	    write(mm2,333)'i=',isez3(iic),'____k=',ksez4(kkc),time-tempo*0.5
	    do 201 j=jy1,jy2
	    write(mm1,10)
     %iplant8,isez3(iic),j,ksez4(kkc),vmodavg(iic,j,kkc),
     %vavg(iic,j,kkc),uavg(iic,j,kkc),wavg(iic,j,kkc),
     %pavg(iic,j,kkc),ptotavg(iic,j,kkc),tvavg(iic,j,kkc)
	    write(mm2,10)
     %iplant8,isez3(iic),j,ksez4(kkc),vorazavg(iic,j,kkc),
     %vorravg(iic,j,kkc),vorzavg(iic,j,kkc)
  201	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(nsez4maxr,1,MPI_INTEGER,jp,1,
     %MPI_COMM_EDDY,STATUS,IERR)
		  ALLOCATE(ksez4r(nsez4maxr))
		  ALLOCATE(vmodavgr(nsez3,ny,nsez4maxr),
     %uavgr(nsez3,ny,nsez4maxr),vavgr(nsez3,ny,nsez4maxr),
     %wavgr(nsez3,ny,nsez4maxr),pavgr(nsez3,ny,nsez4maxr),
     %ptotavgr(nsez3,ny,nsez4maxr),vorravgr(nsez3,ny,nsez4maxr),
     %vorazavgr(nsez3,ny,nsez4maxr),vorzavgr(nsez3,ny,nsez4maxr),
     %tvavgr(nsez3,ny,nsez4maxr))
		  CALL MPI_RECV(ksez4r(1:nsez4maxr),nsez4maxr,MPI_INTEGER,jp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vmodavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(uavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(wavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,6,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(pavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,7,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(ptotavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vorravgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vorazavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(vorzavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(tvavgr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 203 kkc=1,nsez4maxr
		  k=ksez4r(kkc)
		  kg=jp*(nz-2)+k
		  do 203 iic=1,nsez3
		  write(mm1,333)'i=',isez3(iic),'____k=',kg,time-tempo*0.5
		  write(mm2,333)'i=',isez3(iic),'____k=',kg,time-tempo*0.5
		  do 203 j=jy1,jy2
		  write(mm1,10)
     %iplant8,isez3(iic),j,kg,vmodavgr(iic,j,kkc),
     %vavgr(iic,j,kkc),uavgr(iic,j,kkc),wavgr(iic,j,kkc),
     %pavgr(iic,j,kkc),ptotavgr(iic,j,kkc),tvavgr(iic,j,kkc)
		  write(mm2,10)
     %iplant8,isez3(iic),j,kg,vorazavgr(iic,j,kkc),
     %vorravgr(iic,j,kkc),vorzavgr(iic,j,kkc)
  203		  continue
		  DEALLOCATE(ksez4r)
		  DEALLOCATE(vmodavgr,uavgr,vavgr,wavgr,pavgr,
     %ptotavgr,vorravgr,vorazavgr,vorzavgr,tvavgr)
	       enddo
	    endif
	 else
	    CALL MPI_SEND(nsez4max,1,MPI_INTEGER,0,1,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ksez4(1:nsez4max),nsez4max,MPI_INTEGER,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vmodavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(uavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(wavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(pavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ptotavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vorravg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vorazavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(vorzavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(tvavg(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
c
	 do 204 kkc=1,nsez6max
	 do 204 jjc=1,nsez5
	 do 204 i=ix1,ix2
	 rvmodavg(i,jjc,kkc)=rvmodavg(i,jjc,kkc)/
     %tempo
	 rvavg(i,jjc,kkc)=rvavg(i,jjc,kkc)/
     %tempo
	 rvorazavg(i,jjc,kkc)=rvorazavg(i,jjc,kkc)/
     %tempo
	 ruavg(i,jjc,kkc)=ruavg(i,jjc,kkc)/
     %tempo
	 rvorravg(i,jjc,kkc)=rvorravg(i,jjc,kkc)/
     %tempo
	 rwavg(i,jjc,kkc)=rwavg(i,jjc,kkc)/
     %tempo
	 rvorzavg(i,jjc,kkc)=rvorzavg(i,jjc,kkc)/
     %tempo
	 if(tempo2(i,jjc,kkc).ne.0.) then
	    rpavg(i,jjc,kkc)=rpavg(i,jjc,kkc)/
     %tempo2(i,jjc,kkc)
	    rptotavg(i,jjc,kkc)=rptotavg(i,jjc,kkc)/
     %tempo2(i,jjc,kkc)
	 endif
	 rtvavg(i,jjc,kkc)=rtvavg(i,jjc,kkc)/
     %tempo
	 rtvavg(i,jjc,kkc)=rtvavg(i,jjc,kkc)/ru1
 204	 continue
	 if(myrank.eq.0) then
	    mm1=123
	    mm2=124
	    do 205 kkc=1,nsez6max
	    do 205 jjc=1,nsez5
	    write(mm1,333)'j=',jsez5(jjc),'____k=',ksez6(kkc),time-tempo*0.5
	    write(mm2,333)'j=',jsez5(jjc),'____k=',ksez6(kkc),time-tempo*0.5
	    do 205 i=ix1,ix2
	    write(mm1,10)
     %iplant8,i,jsez5(jjc),ksez6(kkc),rvmodavg(i,jjc,kkc),
     %rvavg(i,jjc,kkc),ruavg(i,jjc,kkc),rwavg(i,jjc,kkc),
     %rpavg(i,jjc,kkc),rptotavg(i,jjc,kkc),rtvavg(i,jjc,kkc)
	    write(mm2,10)
     %iplant8,i,jsez5(jjc),ksez6(kkc),rvorazavg(i,jjc,kkc),
     %rvorravg(i,jjc,kkc),rvorzavg(i,jjc,kkc)
 205	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(nsez6maxr,1,MPI_INTEGER,jp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
		  ALLOCATE(ksez6r(nsez6maxr))
		  ALLOCATE(rvmodavgr(nx,nsez5,nsez6maxr),
     %ruavgr(nx,nsez5,nsez6maxr),rvavgr(nx,nsez5,nsez6maxr),
     %rwavgr(nx,nsez5,nsez6maxr),rpavgr(nx,nsez5,nsez6maxr),
     %rptotavgr(nx,nsez5,nsez6maxr),rvorravgr(nx,nsez5,nsez6maxr),
     %rvorazavgr(nx,nsez5,nsez6maxr),rvorzavgr(nx,nsez5,nsez6maxr),
     %rtvavgr(nx,nsez5,nsez6maxr))
		  CALL MPI_RECV(ksez6r(1:nsez6maxr),nsez6maxr,MPI_INTEGER,jp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvmodavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(ruavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,16,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,17,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rwavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,18,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rpavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,19,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rptotavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,20,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvorravgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,21,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvorazavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,22,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rvorzavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,23,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(rtvavgr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,24,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 207 kkc=1,nsez6maxr
		  k=ksez6r(kkc)
		  kg=jp*(nz-2)+k
		  do 207 jjc=1,nsez5
		  write(mm1,333)'j=',jsez5(jjc),'____k=',kg,time-tempo*0.5
		  write(mm2,333)'j=',jsez5(jjc),'____k=',kg,time-tempo*0.5
		  do 207 i=ix1,ix2
		  write(mm1,10)
     %iplant8,i,jsez5(jjc),kg,rvmodavgr(i,jjc,kkc),
     %rvavgr(i,jjc,kkc),ruavgr(i,jjc,kkc),rwavgr(i,jjc,kkc),
     %rpavgr(i,jjc,kkc),rptotavgr(i,jjc,kkc),rtvavgr(i,jjc,kkc)
		  write(mm2,10)
     %iplant8,i,jsez5(jjc),kg,rvorazavgr(i,jjc,kkc),
     %rvorravgr(i,jjc,kkc),rvorzavgr(i,jjc,kkc)
 207		  continue
		  DEALLOCATE(ksez6r)
		  DEALLOCATE(rvmodavgr,ruavgr,rvavgr,rwavgr,rpavgr,
     %rptotavgr,rvorravgr,rvorazavgr,rvorzavgr,rtvavgr)
	       enddo
	    endif
	 else
	    CALL MPI_SEND(nsez6max,1,MPI_INTEGER,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ksez6(1:nsez6max),nsez6max,MPI_INTEGER,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvmodavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ruavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,16,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,17,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rwavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,18,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rpavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,19,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rptotavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,20,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvorravg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,21,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvorazavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,22,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rvorzavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,23,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(rtvavg(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,24,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
c
	 do 208 jjc=1,nsez91
	 do 208 iic=1,nsez9
	 do 208 k=kz1,kz2
	 zvmodavg(iic,jjc,k)=zvmodavg(iic,jjc,k)/tempo
	 zvavg(iic,jjc,k)=zvavg(iic,jjc,k)/tempo
	 zvorazavg(iic,jjc,k)=zvorazavg(iic,jjc,k)/tempo
	 zuavg(iic,jjc,k)=zuavg(iic,jjc,k)/tempo
	 zvorravg(iic,jjc,k)=zvorravg(iic,jjc,k)/tempo
	 zwavg(iic,jjc,k)=zwavg(iic,jjc,k)/tempo
	 zvorzavg(iic,jjc,k)=zvorzavg(iic,jjc,k)/tempo
	 if(ztempo(iic,jjc,k).ne.0.) then
	    zpavg(iic,jjc,k)=
     %zpavg(iic,jjc,k)/ztempo(iic,jjc,k)
	    zptotavg(iic,jjc,k)=
     %zptotavg(iic,jjc,k)/ztempo(iic,jjc,k)
	 endif
	 ztvavg(iic,jjc,k)=ztvavg(iic,jjc,k)/tempo
	 ztvavg(iic,jjc,k)=ztvavg(iic,jjc,k)/ru1
 208	 continue
	 mm1=125
	 mm2=126
	 do 209 jjc=1,nsez91
	 do 209 iic=1,nsez9
	 if(myrank.eq.0) then
	    write(mm1,333)'i=',isez9(iic),'____j=',jsez91(jjc),time-tempo*0.5
	    write(mm2,333)'i=',isez9(iic),'____j=',jsez91(jjc),time-tempo*0.5
	    do 210 k=kz1,kz2
	    write(mm1,10)iplant8,isez9(iic),jsez91(jjc),k,
     %zvmodavg(iic,jjc,k),zvavg(iic,jjc,k),zuavg(iic,jjc,k),
     %zwavg(iic,jjc,k),zpavg(iic,jjc,k),zptotavg(iic,jjc,k),
     %ztvavg(iic,jjc,k)
	    write(mm2,10)iplant8,isez9(iic),jsez91(jjc),k,
     %zvorazavg(iic,jjc,k),zvorravg(iic,jjc,k),zvorzavg(iic,jjc,k)
 210	    continue
	    if(mysize.ne.1) then
	       do jp=1,mysize-1
		  CALL MPI_RECV(zvmodavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,25,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zuavgr(kz1:kz2),(nz-2),MTYPE,jp,26,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvavgr(kz1:kz2),(nz-2),MTYPE,jp,27,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zwavgr(kz1:kz2),(nz-2),MTYPE,jp,28,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zpavgr(kz1:kz2),(nz-2),MTYPE,jp,29,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zptotavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,30,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvorravg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,31,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvorazavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,32,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(zvorzavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,33,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV(ztvavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,34,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 211 k=kz1,kz2
		  kg=jp*(nz-2)+k
		  write(mm1,10)iplant8,isez9(iic),jsez91(jjc),kg,
     %zvmodavg(iic,jjc,k),zvavgr(k),zuavgr(k),zwavgr(k),
     %zpavgr(k),zptotavg(iic,jjc,k),ztvavg(iic,jjc,k)
		  write(mm2,10)iplant8,isez9(iic),jsez91(jjc),kg,
     %zvorazavg(iic,jjc,k),zvorravg(iic,jjc,k),zvorzavg(iic,jjc,k)
 211		  continue
	       enddo
	    endif
	 else
	    CALL MPI_SEND(zvmodavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,25,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zuavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,26,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,27,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zwavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,28,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zpavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,29,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zptotavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,30,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvorravg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,31,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvorazavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,32,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(zvorzavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,33,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND(ztvavg(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,34,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
 209	 continue
c
	 call rmsmodstaggnoavgcart(nx,ny,nz,time)
!      endif
!      if(irms.eq.1)call rmsmodstaggnoavg(nx,ny,nz,uo,vo,wo,p,dt,nbd,flagpo,time)
!      if(tempo.ge.periodo4) then
	 do 110 kkc=1,nsez4max
	 do 110 iic=1,nsez3
	 do 110 j=jy1,jy2
	 vmodavg(iic,j,kkc)=0.
	 vavg(iic,j,kkc)=0.
	 vorazavg(iic,j,kkc)=0.
	 uavg(iic,j,kkc)=0.
	 vorravg(iic,j,kkc)=0.
	 wavg(iic,j,kkc)=0.
	 vorzavg(iic,j,kkc)=0.
	 pavg(iic,j,kkc)=0.
	 ptotavg(iic,j,kkc)=0.
	 tvavg(iic,j,kkc)=0.
	 tempo1(iic,j,kkc)=0.
  110	 continue
	 do 112 kkc=1,nsez6max
	 do 112 jjc=1,nsez5
	 do 112 i=ix1,ix2
	 rvmodavg(i,jjc,kkc)=0.
	 rvavg(i,jjc,kkc)=0.
	 rvorazavg(i,jjc,kkc)=0.
	 ruavg(i,jjc,kkc)=0.
	 rvorravg(i,jjc,kkc)=0.
	 rwavg(i,jjc,kkc)=0.
	 rvorzavg(i,jjc,kkc)=0.
	 rpavg(i,jjc,kkc)=0.
	 rptotavg(i,jjc,kkc)=0.
	 rtvavg(i,jjc,kkc)=0.
	 tempo2(i,jjc,kkc)=0.
  112	 continue
	 do 114 jjc=1,nsez91
	 do 114 iic=1,nsez9
	 do 114 k=kz1,kz2
	 zvmodavg(iic,jjc,k)=0.
	 zvavg(iic,jjc,k)=0.
	 zvorazavg(iic,jjc,k)=0.
	 zuavg(iic,jjc,k)=0.
	 zvorravg(iic,jjc,k)=0.
	 zwavg(iic,jjc,k)=0.
	 zvorzavg(iic,jjc,k)=0.
	 zpavg(iic,jjc,k)=0.
	 zptotavg(iic,jjc,k)=0.
	 ztvavg(iic,jjc,k)=0.
	 ztempo(iic,jjc,k)=0.
  114	 continue
 333	 format(a2,i6,a6,i6,f20.10)
	 tempo=0.
      endif
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Subroutine that prints the RMS and the Reynolds stresses in
c Cartesian coordinates on a cylindrical grid
c
      subroutine rmsmodstaggnoavgcart(nx,ny,nz,time)
c
      use sections
      use mediey
      use mediex
      use mediez
      use rmsy
      use rmsx
      use rmsz
      include'common.h'
      include'averages.h'
      include'mpif.h'
c
c Global variables
      integer nx,ny,nz
      real time
c
c Local variables
      integer iic,i,jjc,j,jp,kkc,k,mm1,mm2,kg,nsez4maxr,nsez6maxr
      integer, dimension(:), allocatable :: ksez4r,ksez6r
      integer STATUS(MPI_STATUS_SIZE)
      real, dimension(:,:,:), allocatable :: urmsr,vrmsr,wrmsr,
     %prmsr,uvsr,uwsr,vwsr,rurmsr,rvrmsr,rwrmsr,rprmsr,ruvsr,
     %ruwsr,rvwsr
c
  10  format(4i6,10f25.8)
c
      do 121 kkc=1,nsez4max
      do 121 iic=1,nsez3
      do 120 j=jy1,jy2
      vrms(iic,j,kkc)=
     %sqrt(abs((vrms(iic,j,kkc)/tempo)
     %-vavg(iic,j,kkc)**2.))
      urms(iic,j,kkc)=
     %sqrt(abs((urms(iic,j,kkc)/tempo)
     %-uavg(iic,j,kkc)**2.))
      wrms(iic,j,kkc)=
     %sqrt(abs((wrms(iic,j,kkc)/tempo)
     %-wavg(iic,j,kkc)**2.))
      uvs(iic,j,kkc)=
     %uvs(iic,j,kkc)/tempo
     %-uavg(iic,j,kkc)*vavg(iic,j,kkc)
      uws(iic,j,kkc)=
     %uws(iic,j,kkc)/tempo
     %-uavg(iic,j,kkc)*wavg(iic,j,kkc)
      vws(iic,j,kkc)=
     %vws(iic,j,kkc)/tempo
     %-vavg(iic,j,kkc)*wavg(iic,j,kkc)
      if(tempo1(iic,j,kkc).ne.0.)
     %prms(iic,j,kkc)=
     %sqrt(abs((prms(iic,j,kkc)/tempo1(iic,j,kkc))
     %-pavg(iic,j,kkc)**2.))
  120 continue
  121 continue
      if(myrank.eq.0) then
	 mm1=127
	 mm2=128
	 do 122 kkc=1,nsez4max
	 do 122 iic=1,nsez3
	 write(mm1,333)'i=',isez3(iic),'____k=',ksez4(kkc),time-0.5*tempo
	 write(mm2,333)'i=',isez3(iic),'____k=',ksez4(kkc),time-0.5*tempo
	 do 122 j=jy1,jy2
	 write(mm1,10)iplant8,isez3(iic),j,ksez4(kkc),
     %vrms(iic,j,kkc),urms(iic,j,kkc),wrms(iic,j,kkc),prms(iic,j,kkc)
	 write(mm2,10)iplant8,isez3(iic),j,ksez4(kkc),
     %uvs(iic,j,kkc),uws(iic,j,kkc),vws(iic,j,kkc)
  122	 continue
	 if(mysize.ne.1) then
	    do jp=1,mysize-1
	       CALL MPI_RECV(nsez4maxr,1,MPI_INTEGER,jp,1,
     %MPI_COMM_EDDY,STATUS,IERR)
	       ALLOCATE(ksez4r(nsez4maxr))
	       ALLOCATE(urmsr(nsez3,ny,nsez4maxr),
     %vrmsr(nsez3,ny,nsez4maxr),wrmsr(nsez3,ny,nsez4maxr),
     %prmsr(nsez3,ny,nsez4maxr),uvsr(nsez3,ny,nsez4maxr),
     %uwsr(nsez3,ny,nsez4maxr),vwsr(nsez3,ny,nsez4maxr))
	       CALL MPI_RECV(ksez4r(1:nsez4maxr),nsez4maxr,MPI_INTEGER,jp,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(urmsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(vrmsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(wrmsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(prmsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(uvsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(uwsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(vwsr(:,jy1:jy2,1:nsez4maxr),nsez3*(ny-2)*nsez4maxr,MTYPE,jp,9,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 123 kkc=1,nsez4maxr
	       k=ksez4r(kkc)
	       kg=jp*(nz-2)+k
	       do 123 iic=1,nsez3
	       write(mm1,333)'i=',isez3(iic),'____k=',kg,time-0.5*tempo
	       write(mm2,333)'i=',isez3(iic),'____k=',kg,time-0.5*tempo
	       do 123 j=jy1,jy2
	       write(mm1,10)
     %iplant8,isez3(iic),j,kg,vrmsr(iic,j,kkc),urmsr(iic,j,kkc),
     %wrmsr(iic,j,kkc),prmsr(iic,j,kkc)
	       write(mm2,10)
     %iplant8,isez3(iic),j,kg,
     %uvsr(iic,j,kkc),uwsr(iic,j,kkc),vwsr(iic,j,kkc)
  123	       continue
	       DEALLOCATE(ksez4r)
	       DEALLOCATE(urmsr,vrmsr,wrmsr,prmsr,uvsr,uwsr,vwsr)
	    enddo
	 endif
      else
	 CALL MPI_SEND(nsez4max,1,MPI_INTEGER,0,1,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(ksez4(1:nsez4max),nsez4max,MPI_INTEGER,0,2,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(urms(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,3,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(vrms(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,4,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(wrms(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,5,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(prms(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,6,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(uvs(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,7,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(uws(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,8,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(vws(:,jy1:jy2,1:nsez4max),nsez3*(ny-2)*nsez4max,MTYPE,0,9,
     %MPI_COMM_EDDY,STATUS,IERR)
      endif
c
      do kkc=1,nsez4max
	 do iic=1,nsez3
	    do j=jy1,jy2
	       vrms(iic,j,kkc)=0.
	       urms(iic,j,kkc)=0.
	       wrms(iic,j,kkc)=0.
	       prms(iic,j,kkc)=0.
	       uvs(iic,j,kkc)=0.
	       uws(iic,j,kkc)=0.
	       vws(iic,j,kkc)=0.
	    enddo
	 enddo
      enddo
c
      do 124 kkc=1,nsez6max
      do 124 jjc=1,nsez5
      do 124 i=ix1,ix2
      rvrms(i,jjc,kkc)=
     %sqrt(abs((rvrms(i,jjc,kkc)/tempo)
     %-rvavg(i,jjc,kkc)**2.))
      rurms(i,jjc,kkc)=
     %sqrt(abs((rurms(i,jjc,kkc)/tempo)
     %-ruavg(i,jjc,kkc)**2.))
      rwrms(i,jjc,kkc)=
     %sqrt(abs((rwrms(i,jjc,kkc)/tempo)
     %-rwavg(i,jjc,kkc)**2.))
      ruvs(i,jjc,kkc)=
     %ruvs(i,jjc,kkc)/tempo
     %-ruavg(i,jjc,kkc)*rvavg(i,jjc,kkc)
      ruws(i,jjc,kkc)=
     %ruws(i,jjc,kkc)/tempo
     %-ruavg(i,jjc,kkc)*rwavg(i,jjc,kkc)
      rvws(i,jjc,kkc)=
     %rvws(i,jjc,kkc)/tempo
     %-rvavg(i,jjc,kkc)*rwavg(i,jjc,kkc)
      if(tempo2(i,jjc,kkc).ne.0.)
     %rprms(i,jjc,kkc)=
     %sqrt(abs((rprms(i,jjc,kkc)/tempo2(i,jjc,kkc))
     %-rpavg(i,jjc,kkc)**2.))
  124 continue
c
      if(myrank.eq.0) then
	 mm1=129
	 mm2=130
	 do 125 kkc=1,nsez6max
	 do 125 jjc=1,nsez5
	 write(mm1,333)'j=',jsez5(jjc),'____k=',ksez6(kkc),time-0.5*tempo
	 write(mm2,333)'j=',jsez5(jjc),'____k=',ksez6(kkc),time-0.5*tempo
	 do 125 i=ix1,ix2
	 write(mm1,10)iplant8,i,jsez5(jjc),ksez6(kkc),
     %rvrms(i,jjc,kkc),rurms(i,jjc,kkc),
     %rwrms(i,jjc,kkc),rprms(i,jjc,kkc)
	 write(mm2,10)iplant8,i,jsez5(jjc),ksez6(kkc),
     %ruvs(i,jjc,kkc),ruws(i,jjc,kkc),rvws(i,jjc,kkc)
  125	 continue
	 if(mysize.ne.1) then
	    do jp=1,mysize-1
	       CALL MPI_RECV(nsez6maxr,1,MPI_INTEGER,jp,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	       ALLOCATE(ksez6r(nsez6maxr))
	       ALLOCATE(rurmsr(nx,nsez5,nsez6maxr),
     %rvrmsr(nx,nsez5,nsez6maxr),rwrmsr(nx,nsez5,nsez6maxr),
     %rprmsr(nx,nsez5,nsez6maxr),ruvsr(nx,nsez5,nsez6maxr),
     %ruwsr(nx,nsez5,nsez6maxr),rvwsr(nx,nsez5,nsez6maxr))
	       CALL MPI_RECV(ksez6r(1:nsez6maxr),nsez6maxr,MPI_INTEGER,jp,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(rurmsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(rvrmsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(rwrmsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(rprmsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(ruvsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,16,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(ruwsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,17,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(rvwsr(ix1:ix2,:,1:nsez6maxr),(nx-2)*nsez5*nsez6maxr,MTYPE,jp,18,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 126 kkc=1,nsez6maxr
	       k=ksez6r(kkc)
	       kg=jp*(nz-2)+k
	       do 126 jjc=1,nsez5
	       write(mm1,333)'j=',jsez5(jjc),'____k=',kg,time-0.5*tempo
	       write(mm2,333)'j=',jsez5(jjc),'____k=',kg,time-0.5*tempo
	       do 126 i=ix1,ix2
	       write(mm1,10)
     %iplant8,i,jsez5(jjc),kg,
     %rvrmsr(i,jjc,kkc),rurmsr(i,jjc,kkc),
     %rwrmsr(i,jjc,kkc),rprmsr(i,jjc,kkc)
	       write(mm2,10)
     %iplant8,i,jsez5(jjc),kg,
     %ruvsr(i,jjc,kkc),ruwsr(i,jjc,kkc),rvwsr(i,jjc,kkc)
  126	       continue
	       DEALLOCATE(ksez6r)
	       DEALLOCATE(rurmsr,rvrmsr,rwrmsr,rprmsr,
     %ruvsr,ruwsr,rvwsr)
	    enddo
	 endif
      else
	 CALL MPI_SEND(nsez6max,1,MPI_INTEGER,0,10,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(ksez6(1:nsez6max),nsez6max,MPI_INTEGER,0,11,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(rurms(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(rvrms(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(rwrms(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(rprms(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(ruvs(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,16,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(ruws(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,17,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(rvws(ix1:ix2,:,1:nsez6max),(nx-2)*nsez5*nsez6max,MTYPE,0,18,
     %MPI_COMM_EDDY,STATUS,IERR)
      endif
c
      do kkc=1,nsez6max
	 do jjc=1,nsez5
	    do i=ix1,ix2
	       rvrms(i,jjc,kkc)=0.
	       rurms(i,jjc,kkc)=0.
	       rwrms(i,jjc,kkc)=0.
	       rprms(i,jjc,kkc)=0.
	       ruvs(i,jjc,kkc)=0.
	       ruws(i,jjc,kkc)=0.
	       rvws(i,jjc,kkc)=0.
	    enddo
	 enddo
      enddo
c
      do 127 jjc=1,nsez91
      do 127 iic=1,nsez9
      do 127 k=kz1,kz2
      zvrms(iic,jjc,k)=
     %sqrt(abs(zvrms(iic,jjc,k)/tempo
     %-zvavg(iic,jjc,k)**2.))
      zurms(iic,jjc,k)=
     %sqrt(abs(zurms(iic,jjc,k)/tempo
     %-zuavg(iic,jjc,k)**2.))
      zwrms(iic,jjc,k)=
     %sqrt(abs(zwrms(iic,jjc,k)/tempo
     %-zwavg(iic,jjc,k)**2.))
      zuvs(iic,jjc,k)=
     %zuvs(iic,jjc,k)/tempo
     %-zuavg(iic,jjc,k)*zvavg(iic,jjc,k)
      zuws(iic,jjc,k)=
     %zuws(iic,jjc,k)/tempo
     %-zuavg(iic,jjc,k)*zwavg(iic,jjc,k)
      zvws(iic,jjc,k)=
     %zvws(iic,jjc,k)/tempo
     %-zvavg(iic,jjc,k)*zwavg(iic,jjc,k)
      if(ztempo(iic,jjc,k).ne.0.)
     %zprms(iic,jjc,k)=
     %sqrt(abs(zprms(iic,jjc,k)/ztempo(iic,jjc,k)
     %-zpavg(iic,jjc,k)**2.))
  127 continue
c
      mm1=131
      mm2=132
      do 128 jjc=1,nsez91
      do 128 iic=1,nsez9
      if(myrank.eq.0) then
	 write(mm1,333)'i=',isez9(iic),'____j=',jsez91(jjc),time-0.5*tempo
	 write(mm2,333)'i=',isez9(iic),'____j=',jsez91(jjc),time-0.5*tempo
	 do 129 k=kz1,kz2
	 write(mm1,10)iplant8,isez9(iic),jsez91(jjc),k,
     %zvrms(iic,jjc,k),zurms(iic,jjc,k),zwrms(iic,jjc,k),
     %zprms(iic,jjc,k)
	 write(mm2,10)iplant8,isez9(iic),jsez91(jjc),k,
     %zuvs(iic,jjc,k),zuws(iic,jjc,k),zvws(iic,jjc,k)
  129	 continue
	 if(mysize.ne.1) then
	    do jp=1,mysize-1
	       CALL MPI_RECV(zurms(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,19,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(zvrms(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,20,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(zwrms(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,21,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(zprms(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,22,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(zuvs(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,23,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(zuws(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,24,
     %MPI_COMM_EDDY,STATUS,IERR)
	       CALL MPI_RECV(zvws(iic,jjc,kz1:kz2),(nz-2),MTYPE,jp,25,
     %MPI_COMM_EDDY,STATUS,IERR)
	       do 130 k=kz1,kz2
	       kg=jp*(nz-2)+k
	       write(mm1,10)iplant8,isez9(iic),jsez91(jjc),kg,
     %zvrms(iic,jjc,k),zurms(iic,jjc,k),zwrms(iic,jjc,k),
     %zprms(iic,jjc,k)
	       write(mm2,10)iplant8,isez9(iic),jsez91(jjc),kg,
     %zuvs(iic,jjc,k),zuws(iic,jjc,k),zvws(iic,jjc,k)
 130	       continue
	    enddo
	 endif
      else
	 CALL MPI_SEND(zurms(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,19,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(zvrms(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,20,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(zwrms(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,21,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(zprms(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,22,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(zuvs(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,23,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(zuws(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,24,
     %MPI_COMM_EDDY,STATUS,IERR)
	 CALL MPI_SEND(zvws(iic,jjc,kz1:kz2),(nz-2),MTYPE,0,25,
     %MPI_COMM_EDDY,STATUS,IERR)
      endif
 128  continue
c
      do jjc=1,nsez91
	 do iic=1,nsez9
	    do k=kz1,kz2
	       zvrms(iic,jjc,k)=0.
	       zurms(iic,jjc,k)=0.
	       zwrms(iic,jjc,k)=0.
	       zprms(iic,jjc,k)=0.
	       zuvs(iic,jjc,k)=0.
	       zuws(iic,jjc,k)=0.
	       zvws(iic,jjc,k)=0.
	    enddo
	 enddo
      enddo
 333  format(a2,i6,a6,i6,f20.10)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Subroutine that prints the velocity components in Cartesian
c coordinates on a cylindrical grid
c
      subroutine plotinststaggcart(iplant9c,iplant9d,iplant9p,
     %nx,ny,nz,nzg,nbd,ntime,xc,xu,yc,yv,zcg,zwg,p,vo,uo,wo,
     %flaguo,flagvo,flagwo)
c
      use sections
      include'common.h'
      include'mpif.h'
      include'averages.h'
c
c Global variables
      integer iplant9c,iplant9d,iplant9p,nx,ny,nz,nzg,nbd,ntime
      integer flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),
     %flagwo(nx,ny,nz,nbd)
      real xc(nx),xu(nx),yc(ny),yv(ny),zcg(nzg),zwg(nzg)
      real p(nx,ny,nz),vo(nx,ny,nz),uo(nx,ny,nz),wo(nx,ny,nz)
c
c Local variables
      integer iic,iicd,iicu,is,im,ip,kkc,kkcd,kkcu,ks,km,kp,
     %jjc,jjcd,jjcu,js,jm,jp,jsy,jpsy,i,j,k,l,ir,jr,kr,lr,
     %i1,jjp,kk1,kk2,kk,kkd,kku
      integer STATUS(MPI_STATUS_SIZE)
!      integer flaguor(ny,nzg,nbd),
!     %flagvor2(nx,nzg,nbd),flagvor2sym(nx,nzg,nbd)
      real pstagg,vstagg,ustagg,wstagg
      real pr(ny,nzg),vor(ny,nzg),uor(ny,nzg),wor(ny,nzg),
     %pr2(nx,nzg),pr2sym(nx,nzg),vor2(nx,nzg),vor2sym(nx,nzg),
     %uor2(nx,nzg),uor2sym(nx,nzg),wor2(nx,nzg),wor2sym(nx,nzg)
     %,pbuff(nx,ny),xbuff(nx,ny),ybuff(nx,ny)
c
c Parameters
      integer lamb2,nxp,nyp,nzp
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
      if(iplant9.eq.1) then
	 if(myrank.eq.0) then
	    do iic=1,nsez1
	       iicd=iic/10
	       iicu=mod(iic,10)
	       open
     %(81,file='grid1_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'_plotinst.xyz')
	       open
     %(82,file='grid1_new_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'_plotinst.xyz')
	       is=isez1(iic)
	       write(81,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	       write(82,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	       do 1 k=2,nzg-1,nzp
	       do 10 l=1,lamb2
	       do 10 j=jy1,jy2,nyp
	       write(81,9)xu(is)*cos(yc(jy1))
   10	       write(82,9)xu(is)*cos(yc(j)+(l-1)*2.*pi/lamb2)
	       write(81,9)xu(is)*cos(yc(jy1))
   1	       write(82,9)xu(is)*cos(yc(jy1)+lamb2*deltasect)
	       do 2 k=2,nzg-1,nzp
	       do 20 l=1,lamb2
	       do 20 j=jy1,jy2,nyp
	       write(81,9)
     %xu(is)*sin(yc(jy1))+xu(is)*(yc(j)+(l-1)*2.*pi/lamb2)
   20	       write(82,9)xu(is)*sin(yc(j)+(l-1)*2.*pi/lamb2)
	       write(81,9)
     %xu(is)*sin(yc(jy1))+xu(is)*deltasect*lamb2
   2	       write(82,9)xu(is)*sin(yc(jy1)+lamb2*deltasect)
	       do 3 k=2,nzg-1,nzp
	       do 30 l=1,lamb2
	       do 30 j=jy1,jy2,nyp
	       write(81,9)zcg(k)
   30	       write(82,9)zcg(k)
	       write(81,9)zcg(k)
   3	       write(82,9)zcg(k)
	       close(81)
	       close(82)
	    enddo
	    do kkc=1,nsez2
	       kkcd=kkc/10
	       kkcu=mod(kkc,10)
	       open(81,file='grid2_2D_slice'//char(48+kkcd)
     %//char(48+kkcu)//'_plotinst.xyz')
	       ks=ksez2g(kkc)
	       write(81,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	       do 11 i=ix1,ix2,nxp
	       do 110 l=1,lamb2
	       do 110 j=jy1,jy2,nyp
  110	       write(81,9)xc(i)*cos(yc(j)+(l-1)*2.*pi/lamb2)
  11	       write(81,9)xc(i)*cos(yc(jy1)+lamb2*deltasect)
	       do 12 i=ix1,ix2,nxp
	       do 120 l=1,lamb2
	       do 120 j=jy1,jy2,nyp
  120	       write(81,9)xc(i)*sin(yc(j)+(l-1)*2.*pi/lamb2)
  12	       write(81,9)xc(i)*sin(yc(jy1)+lamb2*deltasect)
	       do 13 i=ix1,ix2,nxp
	       do 130 l=1,lamb2
	       do 130 j=jy1,jy2,nyp
  130	       write(81,9)zwg(ks)
  13	       write(81,9)zwg(ks)
	       close(81)
	    enddo
	    do jjc=1,nsez14
	       jjcd=jjc/10
	       jjcu=mod(jjc,10)
	       open(81,file='grid3_2D_slice'//char(48+jjcd)
     %//char(48+jjcu)//'_plotinst.xyz')
	       js=jsez14(jjc)
	       jsy=jsym(js)
	       write(81,101)1,2*(nx-2)/nxp,(nzg-2)/nzp
	       do 210 k=2,nzg-1,nzp
	       do 211 i=ix2,ix1,-nxp
  211	       write(81,9)xc(i)*cos(yv(jsy))
	       do 212 i=ix1,ix2,nxp
  212	       write(81,9)xc(i)*cos(yv(js))
  210	       continue
	       do 220 k=2,nzg-1,nzp
	       do 221 i=ix2,ix1,-nxp
  221	       write(81,9)xc(i)*sin(yv(jsy))
	       do 222 i=ix1,ix2,nxp
  222	       write(81,9)xc(i)*sin(yv(js))
  220	       continue
	       do 230 k=2,nzg-1,nzp
	       do 230 i=1,2*(nx-2),nxp
  230	       write(81,9)zcg(k)
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
	    write(82,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	    i1=1
	    write(82,102)ntime,i1,i1,i1
!
!Stagnation pressure
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
	       vor(j,k)=ustagg*cos(yc(j))-vstagg*sin(yc(j))
	       uor(j,k)=ustagg*sin(yc(j))+vstagg*cos(yc(j))
	       wor(j,k)=wstagg
	       pr(j,k)=pstagg
	       write(82,9)
     %(pstagg+0.5*(vstagg**2.+ustagg**2.+
     %wstagg**2.))	!*ros
c	    else
c	       write(82,9)1.e+03
c	    endif
 41	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       pstagg=
     %0.5*(p(is,jy1,k)+p(ip,jy1,k))
	       vstagg=
     %0.25*(vo(is,jy1,k)+vo(ip,jy1,k)+vo(is,jy2,k)+vo(ip,jy2,k))
	       ustagg=
     %uo(is,jy1,k)
	       wstagg=
     %0.25*(wo(is,jy1,k)+wo(ip,jy1,k)+wo(is,jy1,km)+wo(ip,jy1,km))
	       vor(jy1,k)=ustagg*cos(yc(jy1))-vstagg*sin(yc(jy1))
	       uor(jy1,k)=ustagg*sin(yc(jy1))+vstagg*cos(yc(jy1))
	       wor(jy1,k)=wstagg
	       pr(jy1,k)=pstagg
	       write(82,9)
     %(pstagg+0.5*(vstagg**2.+ustagg**2.+
     %wstagg**2.))     !*ros
c	    else
c	       write(82,9)1.e+03
c	    endif
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
		     vstagg=uor(jr,kr)*cos(yc(jr))-vor(jr,kr)*sin(yc(jr))
		     ustagg=uor(jr,kr)*sin(yc(jr))+vor(jr,kr)*cos(yc(jr))
		     vor(jr,kr)=vstagg
		     uor(jr,kr)=ustagg
		     write(82,9)
     %(pr(jr,kr)+0.5*(vor(jr,kr)**2.+uor(jr,kr)**2.+
     %wor(jr,kr)**2.))	    !*ros
c		  else
c		     write(82,9)1.e+03
c		  endif
 43		  continue
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     vstagg=uor(jy1,kr)*cos(yc(jy1))-vor(jy1,kr)*sin(yc(jy1))
		     ustagg=uor(jy1,kr)*sin(yc(jy1))+vor(jy1,kr)*cos(yc(jy1))
		     vor(jy1,kr)=vstagg
		     uor(jy1,kr)=ustagg
		     write(82,9)
     %(pr(jy1,kr)+0.5*(vor(jy1,kr)**2.+uor(jy1,kr)**2.+
     %wor(jy1,kr)**2.))	    !*ros
c		  else
c		     write(82,9)1.e+03
c		  endif
 42		  continue
	       enddo
	    endif
!	    write(82,*)	  !!!!!!
!	    write(82,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
!
! X velocity
!
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)vor(j,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
   51	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)vor(jy1,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
   50	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 52 kr=kk1,kk2,nzp
		  do 53 lr=1,lamb2
		  do 53 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)vor(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   53		  continue
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,9)vor(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   52		  continue
	       enddo
	    endif
!
! Y velocity
!
	    do 60 k=kz1,kz2,nzp
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)uor(j,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
   61	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)uor(jy1,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
   60	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 62 kr=kk1,kk2,nzp
		  do 63 lr=1,lamb2
		  do 63 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)uor(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   63		  continue
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,9)uor(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   62		  continue
	       enddo
	    endif
!
! Z velocity
!
	    do 70 k=kz1,kz2,nzp
	    do 71 l=1,lamb2
	    do 71 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)wor(j,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
   71	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)wor(jy1,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,9)wor(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
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
	       write(82,8)pr(j,k)   !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
   81	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,8)pr(jy1,k)   !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,8)pr(jy1,kr)   !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
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
      do kkc=1,nsez2max
	 ks=ksez2(kkc)
	 kp=ks+1
	 kk=isezindx2+(kkc-1)
	 kkd=kk/10
	 kku=mod(kk,10)
	 open(82,file='results2_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+kkd)//char(48+kku)//'_plotinst.q')
	 write(82,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
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
	    xbuff(i,j)=ustagg*cos(yc(j))-vstagg*sin(yc(j))
	    ybuff(i,j)=ustagg*sin(yc(j))+vstagg*cos(yc(j))
	    pbuff(i,j)=pstagg
	    write(82,9)
     %(pstagg+0.5*(vstagg**2.+ustagg**2.+wstagg**2.))	!*ros
c	 else
c	    write(82,9)1.e+03
c	 endif
 141	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    pstagg=
     %0.5*(p(i,jy1,ks)+p(i,jy1,kp))
	    vstagg=
     %0.25*(vo(i,jy1,ks)+vo(i,jy2,ks)+vo(i,jy1,kp)+vo(i,jy2,kp))
	    ustagg=
     %0.25*(uo(i,jy1,ks)+uo(im,jy1,ks)+uo(i,jy1,kp)+uo(im,jy1,kp))
	    wstagg=
     %wo(i,jy1,ks)
	    xbuff(i,jy1)=ustagg*cos(yc(jy1))-vstagg*sin(yc(jy1))
	    ybuff(i,jy1)=ustagg*sin(yc(jy1))+vstagg*cos(yc(jy1))
	    pbuff(i,jy1)=pstagg
	    write(82,9)
     %(pstagg+0.5*(vstagg**2.+ustagg**2.+wstagg**2.))	 !*ros
c	 else
c	    write(82,9)1.e+03
c	 endif
 140	 continue
!	 write(82,*)	!!!!!!
!	 write(82,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
!
! X velocity
!
	 do 150 i=ix1,ix2,nxp
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)xbuff(i,j)
c	 else
c	    write(82,9)1.e+03
c	 endif
 151	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)xbuff(i,jy1)
c	 else
c	    write(82,9)1.e+03
c	 endif
 150	 continue
!
! Y velocity
!
	 do 160 i=ix1,ix2,nxp
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)ybuff(i,j)
c	 else
c	    write(82,9)1.e+03
c	 endif
 161	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)ybuff(i,jy1)
c	 else
c	    write(82,9)1.e+03
c	 endif
 160	 continue
!
! Z velocity
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
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)wo(i,jy1,ks)
c	 else
c	    write(82,9)1.e+03
c	 endif
 170	 continue
!
!Static pressure
!
	 do 180 i=ix1,ix2,nxp
	 do 181 l=1,lamb2
	 do 181 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,8)pbuff(i,j)     !*ros
c	 else
c	    write(82,8)1.e+03
c	 endif
 181	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,8)pbuff(i,jy1)    !*ros
c	 else
c	    write(82,8)1.e+03
c	 endif
 180	 continue
	 close(82)
      enddo
c
      do jjc=1,nsez14
	 js=jsez14(jjc)
	 jsy=jsym(js)
	 jp=jpv(js)
	 jpsy=jpv(jsy)
	 if(myrank.eq.0) then
	    jjcd=jjc/10
	    jjcu=mod(jjc,10)
	    open(82,file='results3_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_plotinst.q')
	    write(82,101)1,2*(nx-2)/nxp,(nzg-2)/nzp
	    i1=1
	    write(82,102)ntime,i1,i1,i1
!
!Stagnation pressure
!
	    do 240 k=kz1,kz2,nzp
	    km=k-1
	    do 241 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       im=i-1
	       pstagg=
     %0.5*(p(i,jsy,k)+p(i,jpsy,k))
	       vstagg=
     %vo(i,jsy,k)
	       ustagg=
     %0.25*(uo(i,jsy,k)+uo(im,jsy,k)+uo(i,jpsy,k)+uo(im,jpsy,k))
	       wstagg=
     %0.25*(wo(i,jsy,k)+wo(i,jpsy,k)+wo(i,jsy,km)+wo(i,jpsy,km))
	       vor2sym(i,k)=ustagg*cos(yv(jsy))-vstagg*sin(yv(jsy))
	       uor2sym(i,k)=ustagg*sin(yv(jsy))+vstagg*cos(yv(jsy))
	       wor2sym(i,k)=wstagg
	       pr2sym(i,k)=pstagg
	       write(82,9)
     %(pstagg+0.5*(vstagg**2.+ustagg**2.+wstagg**2.))	 !*ros
c	    else
c	       write(82,9)1.e+03
c	    endif
 241	    continue
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
	       vor2(i,k)=ustagg*cos(yv(js))-vstagg*sin(yv(js))
	       uor2(i,k)=ustagg*sin(yv(js))+vstagg*cos(yv(js))
	       wor2(i,k)=wstagg
	       pr2(i,k)=pstagg
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
c		  CALL MPI_RECV
c    %(flagvor2sym(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,11,
c    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(pr2sym(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(vor2sym(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(uor2sym(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(wor2sym(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 243 kr=kk1,kk2,nzp
		  do 244 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     vstagg=uor2sym(ir,kr)*cos(yv(jsy))-vor2sym(ir,kr)*sin(yv(jsy))
		     ustagg=uor2sym(ir,kr)*sin(yv(jsy))+vor2sym(ir,kr)*cos(yv(jsy))
		     vor2sym(ir,kr)=vstagg
		     uor2sym(ir,kr)=ustagg
		     write(82,9)
     %(pr2sym(ir,kr)+0.5*(vor2sym(ir,kr)**2.+uor2sym(ir,kr)**2.
     %+wor2sym(ir,kr)**2.))    !*ros
c		  else
c		     write(82,9)1.e+03
c		  endif
 244		  continue
		  do 245 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     vstagg=uor2(ir,kr)*cos(yv(js))-vor2(ir,kr)*sin(yv(js))
		     ustagg=uor2(ir,kr)*sin(yv(js))+vor2(ir,kr)*cos(yv(js))
		     vor2(ir,kr)=vstagg
		     uor2(ir,kr)=ustagg
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
!
! X velocity
!
	    do 250 k=kz1,kz2,nzp
	    do 251 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,9)vor2sym(i,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 251	    continue
	    do 252 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,9)vor2(i,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 252	    continue
 250	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 253 kr=kk1,kk2,nzp
		  do 254 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,9)vor2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 254		  continue
		  do 255 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)vor2(ir,kr)
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
	    do 261 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,9)uor2sym(i,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 261	    continue
	    do 262 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,9)uor2(i,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 262	    continue
 260	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 263 kr=kk1,kk2,nzp
		  do 264 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,9)uor2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 264		  continue
		  do 265 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)uor2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 265		  continue
 263		  continue
	       enddo
	    endif
!
! Z velocity
!
	    do 270 k=kz1,kz2,nzp
	    do 271 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,9)wor2sym(i,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 271	    continue
	    do 272 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,9)wor2(i,k)
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
		  do 274 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,9)wor2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 274		  continue
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
	    do 281 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,8)pr2sym(i,k)	 !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
 281	    continue
	    do 282 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,8)pr2(i,k)    !*ros
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
		  do 284 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,8)pr2sym(ir,kr)	 !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
 284		  continue
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
c	    CALL MPI_SEND
c    %(flagvo(ix1:ix2,jsy,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,11,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((p(ix1:ix2,jsy,kz1:kz2)+p(ix1:ix2,jpsy,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vo(ix1:ix2,jsy,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((uo(1:ix2-1,jsy,kz1:kz2)+uo(ix1:ix2,jsy,kz1:kz2)+
     %uo(1:ix2-1,jpsy,kz1:kz2)+uo(ix1:ix2,jpsy,kz1:kz2))*0.25,
     %(nx-2)*1*(nz-2),MTYPE,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((wo(ix1:ix2,jsy,1:kz2-1)+wo(ix1:ix2,jpsy,1:kz2-1)+
     %wo(ix1:ix2,jsy,kz1:kz2)+wo(ix1:ix2,jpsy,kz1:kz2))*0.25,
     %(nx-2)*1*(nz-2),MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Subroutines that evaluates the instantaneous vorticity components
c in Cartesian coordinates on a cylindrical grid
c
      subroutine vort2staggcart(nx,ny,nz,xc,xu,yc,yv,zc,uo,vo,wo)
c
      use sections
      use vorticity
      include'common.h'
      include'averages.h'
c
c Global variables
      integer nx,ny,nz
      real xc(nx),xu(nx),yc(ny),yv(ny),zc(nz),
     %uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
c
c Local variables
      integer iic,ic,ip,im,jjc,jc,jp,jm,jsy,jpsy,kkc,kc,kp,km
      real dtheta,dr,dz,dq1x3,dq3x1,dq2x3,dq3x2,dq2x1,dq1x2
      real xstagg,ystagg
c
      ! Azimuthal vorticity - Circumferential sections !
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
      ! Azimuthal vorticity - Cross sections !
      do 111 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
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
	do 113 ic=ix1+1,ix2
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
      ! Azimuthal vorticity - Meridian sections !
      do 121 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jsy=jsym(jc)
	jpsy=jsym(jp)
	do 121 kc=kz1,kz2
	  kp=kc+1
	  km=kc-1
	  dz=zc(kp)-zc(km)
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
	  dq1x3=(0.25*(uo(ic,jsy,kp)+uo(im,jsy,kp)+
     %uo(ic,jpsy,kp)+uo(im,jpsy,kp))
     %-0.25*(uo(ic,jsy,km)+uo(im,jsy,km)+
     %uo(ic,jpsy,km)+uo(im,jpsy,km)))/dz
	  dq3x1=(0.25*(wo(ip,jsy,kc)+wo(ip,jpsy,kc)+
     %wo(ip,jsy,km)+wo(ip,jpsy,km))
     %-0.25*(wo(ic,jc,kc)+wo(ic,jp,kc)+
     %wo(ic,jc,km)+wo(ic,jp,km)))/dr
	  voraz14(ic,nsez14+jjc,kc)=dq1x3-dq3x1
	  do 121 ic=ix1+1,ix2
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
	    dq1x3=(0.25*(uo(ic,jsy,kp)+uo(im,jsy,kp)+
     %uo(ic,jpsy,kp)+uo(im,jpsy,kp))
     %-0.25*(uo(ic,jsy,km)+uo(im,jsy,km)+
     %uo(ic,jpsy,km)+uo(im,jpsy,km)))/dz
	    dq3x1=(0.25*(wo(ip,jsy,kc)+wo(ip,jpsy,kc)+
     %wo(ip,jsy,km)+wo(ip,jpsy,km))
     %-0.25*(wo(im,jsy,kc)+wo(im,jpsy,kc)+
     %wo(im,jsy,km)+wo(im,jpsy,km)))/dr
	    voraz14(ic,nsez14+jjc,kc)=dq1x3-dq3x1
  121 continue
c
      ! Radial vorticity - Circumferential sections !
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
	    vorr1(iic,jc,kc)=dq3x2/xu(ic)-dq2x3
  131 continue
c
      ! Radial vorticity - Cross sections !
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
	    vorr2(ic,jc,kkc)=dq3x2/xc(ic)-dq2x3
  141 continue
c
      ! Radial vorticity - Meridian sections !
      do 151 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jsy=jsym(jc)
	jpsy=jpv(jsy)
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
	    vorr14(ic,jjc,kc)=dq3x2/xc(ic)-dq2x3
	    dq3x2=(0.5*(wo(ic,jpsy,kc)+wo(ic,jpsy,km))
     %-0.5*(wo(ic,jsy,kc)+wo(ic,jsy,km)))/dtheta   !
	    dq2x3=(vo(ic,jsy,kp)
     %-vo(ic,jsy,km))/dz
	    vorr14(ic,nsez14+jjc,kc)=dq3x2/xc(ic)-dq2x3
  151 continue
c
      ! Axial vorticity - Circumferential sections !
      do 161 kc=kz1,kz2
	do 161 iic=1,nsez1
	  ic=isez1(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 161 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.5*(vo(ip,jc,kc)+vo(ip,jm,kc))*xc(ip)
     %-0.5*(vo(ic,jc,kc)+vo(ic,jm,kc))*xc(ic))/dr
	    dq1x2=(uo(ic,jp,kc)
     %-uo(ic,jm,kc))/dtheta    !
	    vorz1(iic,jc,kc)=(dq2x1-dq1x2)/xu(ic)
  161 continue
c
      ! Axial vorticity - Cross sections !
      do 171 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
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
     %+vo(ic,jm,kp)+vo(ip,jm,kp))*xu(ic))/dr
	  dq1x2=(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)+
     %uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jm,kc)+uo(im,jm,kc)+
     %uo(ic,jm,kp)+uo(im,jm,kp)))/dtheta    !
	  vorz2(ic,jc,kkc)=(dq2x1-dq1x2)/xc(ic)
  172	continue
c
	do 173 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 173 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.25*(vo(ip,jc,kc)+vo(ip,jm,kc)+
     %vo(ip,jc,kp)+vo(ip,jm,kp))*xc(ip)
     %-0.25*(vo(im,jc,kc)+vo(im,jm,kc)+
     %vo(im,jc,kp)+vo(im,jm,kp))*xc(im))/dr
	    dq1x2=
     %(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)+
     %uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jm,kc)+uo(im,jm,kc)+
     %uo(ic,jm,kp)+uo(im,jm,kp)))/dtheta   !
	    vorz2(ic,jc,kkc)=(dq2x1-dq1x2)/xc(ic)
  173 continue
  171 continue
c
      ! Axial vorticity - Meridian sections !
      do 181 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jsy=jsym(jc)
	jpsy=jsym(jp)
	dtheta=dely
	do 181 kc=kz1,kz2
	  ic=ix1
	  im=ic-1
	  ip=ic+1
	  dr=xu(ic)-xu(im)
	  dq2x1=(0.5*(vo(ip,jc,kc)+vo(ic,jc,kc))*xu(ic))/dr
	  dq1x2=(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dtheta	  !
	  vorz14(ic,jjc,kc)=(dq2x1-dq1x2)/xc(ic)
	  dq2x1=(0.5*(vo(ip,jsy,kc)+vo(ic,jsy,kc))*xu(ic))/dr
	  dq1x2=(0.5*(uo(ic,jpsy,kc)+uo(im,jpsy,kc))
     %-0.5*(uo(ic,jsy,kc)+uo(im,jsy,kc)))/dtheta   !
	  vorz14(ic,nsez14+jjc,kc)=(dq2x1-dq1x2)/xc(ic)
	  do 181 ic=ix1+1,ix2
	    ip=ic+1
	    im=ic-1
	    dr=xc(ip)-xc(im)
	    dq2x1=(vo(ip,jc,kc)*xc(ip)
     %-vo(im,jc,kc)*xc(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dtheta	  !
	    vorz14(ic,jjc,kc)=(dq2x1-dq1x2)/xc(ic)
	    dq2x1=(vo(ip,jsy,kc)*xc(ip)
     %-vo(im,jsy,kc)*xc(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jpsy,kc)+uo(im,jpsy,kc))
     %-0.5*(uo(ic,jsy,kc)+uo(im,jsy,kc)))/dtheta   !
	    vorz14(ic,nsez14+jjc,kc)=(dq2x1-dq1x2)/xc(ic)
  181 continue
c
c
      ! Circumferential sections - conversion from cylindrical to Cartesian !
      do 201 kc=kz1,kz2
	do 201 iic=1,nsez1
	  do 201 jc=jy1,jy2
	    xstagg=
     %vorr1(iic,jc,kc)*cos(yc(jc))-voraz1(iic,jc,kc)*sin(yc(jc))
	    ystagg=
     %vorr1(iic,jc,kc)*sin(yc(jc))+voraz1(iic,jc,kc)*cos(yc(jc))
	    voraz1(iic,jc,kc)=xstagg
	    vorr1(iic,jc,kc)=ystagg
  201 continue
c
      ! Cross sections - conversion from cylindrical to Cartesian !
      do 202 kkc=1,nsez2max
	do 202 ic=ix1,ix2
	  do 202 jc=jy1,jy2
	    xstagg=
     %vorr2(ic,jc,kkc)*cos(yc(jc))-voraz2(ic,jc,kkc)*sin(yc(jc))
	    ystagg=
     %vorr2(ic,jc,kkc)*sin(yc(jc))+voraz2(ic,jc,kkc)*cos(yc(jc))
	    voraz2(ic,jc,kkc)=xstagg
	    vorr2(ic,jc,kkc)=ystagg
  202 continue
c
      ! Meridian sections - conversion from cylindrical to Cartesian !
      do 203 jjc=1,nsez14
	jc=jsez14(jjc)
	jsy=jsym(jc)
	do 203 kc=kz1,kz2
	  do 203 ic=ix1,ix2
	    xstagg=
     %vorr14(ic,jjc,kc)*cos(yv(jc))-voraz14(ic,jjc,kc)*sin(yv(jc))
	    ystagg=
     %vorr14(ic,jjc,kc)*sin(yv(jc))+voraz14(ic,jjc,kc)*cos(yv(jc))
	    voraz14(ic,jjc,kc)=xstagg
	    vorr14(ic,jjc,kc)=ystagg
	    xstagg=
     %vorr14(ic,nsez14+jjc,kc)*cos(yv(jsy))-voraz14(ic,nsez14+jjc,kc)*sin(yv(jsy))
	    ystagg=
     %vorr14(ic,nsez14+jjc,kc)*sin(yv(jsy))+voraz14(ic,nsez14+jjc,kc)*cos(yv(jsy))
	    voraz14(ic,nsez14+jjc,kc)=xstagg
	    vorr14(ic,nsez14+jjc,kc)=ystagg
  203 continue
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Subroutine for the evaluation of the statistics on the cylindrical
c grid in Cartesian coordinates
c
      subroutine mediecalcstaggcart
     %(nx,ny,nz,nzg,nbd,time,dt,vo,uo,wo,p,tv,xc,xu,yc,yv,zc,
     %zcg,zwg,flagpo,flaguo,flagvo,flagwo,ntime)
!    %mbd,nfacet,unvect,vertex)
c
      use sections
      use mediecalc
      use mediecalc1
      use mediecalc2
      use mediecalc3
      use vortavg
      include'common.h'
      include'averages.h'
!     include'immersed.h'
c
c Global variables
      integer nx,ny,nz,nzg,nbd,ntime
      integer flaguo(nx,ny,nz,nbd),flagvo(nx,ny,nz,nbd),
     %flagwo(nx,ny,nz,nbd),flagpo(nx,ny,nz,nbd)
      real time,dt
      real vo(nx,ny,nz),uo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz),
     %xc(nx),xu(nx),yc(ny),yv(ny),zc(nz),zcg(nzg),zwg(nzg),
     %tv(nx,ny,nz)
!     integer mbd,nfacet
!     real unvect(3,nfacet),vertex(3,3,nfacet)
c
c Local variables
      integer iic,is,ip,im,i,jjc,jp,jm,j,js,jsy,jpsy,kkc,ks,kp,km,k,
     %iplant4c,iplant4d,iplant4p,iplant5c,iplant5d,iplant5p
      real vstagg,ustagg,wstagg,pstagg,vmod,tvstagg,xstagg,ystagg
!     integer ibd
c
      tempo3=tempo3+dt
!
!     Time averaged values and root mean squares
!     at some circumferential sections
!
      do 101 iic=1,nsez1
      is=isez1(iic)
      ip=is+1
      do 100 k=kz1,kz2
      km=k-1
      do 100 j=jy1,jy2
      jm=jmv(j)
      vstagg=(vo(is,j,k)+vo(ip,j,k)+
     %vo(is,jm,k)+vo(ip,jm,k))*0.25
      ustagg=uo(is,j,k)
      wstagg=(wo(is,j,k)+wo(ip,j,k)+
     %wo(is,j,km)+wo(ip,j,km))*0.25
      xstagg=ustagg*cos(yc(j))-vstagg*sin(yc(j))
      ystagg=ustagg*sin(yc(j))+vstagg*cos(yc(j))
      vstagg=xstagg
      ustagg=ystagg
      vplot(iic,j,k)=vplot(iic,j,k)+
     %vstagg*dt
      uplot(iic,j,k)=uplot(iic,j,k)+
     %ustagg*dt
      wplot(iic,j,k)=wplot(iic,j,k)+
     %wstagg*dt
      vrmsplot(iic,j,k)=vrmsplot(iic,j,k)+
     %dt*(vstagg**2.)
      urmsplot(iic,j,k)=urmsplot(iic,j,k)+
     %dt*(ustagg**2.)
      wrmsplot(iic,j,k)=wrmsplot(iic,j,k)+
     %dt*(wstagg**2.)
      if(all(flagpo(is:ip,j,k,:).le.0)) then
	 tempo4(iic,j,k)=tempo4(iic,j,k)+dt
	 pstagg=(p(is,j,k)+p(ip,j,k))*0.5
	 pplot(iic,j,k)=pplot(iic,j,k)+
     %pstagg*dt
	 prmsplot(iic,j,k)=prmsplot(iic,j,k)+
     %dt*(pstagg**2.)
      endif
      vmod=sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
      vmodplot(iic,j,k)=vmodplot(iic,j,k)+
     %vmod*dt
      tvstagg=(tv(is,j,k)+tv(ip,j,k))*0.5
      tvplot(iic,j,k)=tvplot(iic,j,k)+
     %tvstagg*dt
      uvplot(iic,j,k)=uvplot(iic,j,k)+
     %ustagg*vstagg*dt
      uwplot(iic,j,k)=uwplot(iic,j,k)+
     %ustagg*wstagg*dt
      vwplot(iic,j,k)=vwplot(iic,j,k)+
     %vstagg*wstagg*dt
  100 continue
  101 continue
!
!     Time averaged values and root mean squares
!     at some cross sections
!
      do 103 kkc=1,nsez2max
      ks=ksez2(kkc)
      kp=ks+1
      do 102 i=ix1,ix2
      im=i-1
      do 102 j=jy1,jy2
      jm=jmv(j)
      vstagg=(vo(i,j,ks)+vo(i,jm,ks)+
     %vo(i,j,kp)+vo(i,jm,kp))*0.25
      ustagg=(uo(i,j,ks)+uo(im,j,ks)+
     %uo(i,j,kp)+uo(im,j,kp))*0.25
      wstagg=wo(i,j,ks)
      xstagg=ustagg*cos(yc(j))-vstagg*sin(yc(j))
      ystagg=ustagg*sin(yc(j))+vstagg*cos(yc(j))
      vstagg=xstagg
      ustagg=ystagg
      vplot2(i,j,kkc)=vplot2(i,j,kkc)+
     %vstagg*dt
      uplot2(i,j,kkc)=uplot2(i,j,kkc)+
     %ustagg*dt
      wplot2(i,j,kkc)=wplot2(i,j,kkc)+
     %wstagg*dt
      vrmsplot2(i,j,kkc)=vrmsplot2(i,j,kkc)+
     %dt*(vstagg**2.)
      urmsplot2(i,j,kkc)=urmsplot2(i,j,kkc)+
     %dt*(ustagg**2.)
      wrmsplot2(i,j,kkc)=wrmsplot2(i,j,kkc)+
     %dt*(wstagg**2.)
      if(all(flagpo(i,j,ks:kp,:).le.0)) then
	 tempo5(i,j,kkc)=tempo5(i,j,kkc)+dt
	 pstagg=(p(i,j,ks)+p(i,j,kp))*0.5
	 pplot2(i,j,kkc)=pplot2(i,j,kkc)+
     %pstagg*dt
	 prmsplot2(i,j,kkc)=prmsplot2(i,j,kkc)+
     %dt*(pstagg**2.)
      endif
      vmod=sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
      vmodplot2(i,j,kkc)=vmodplot2(i,j,kkc)+
     %vmod*dt
      tvstagg=(tv(i,j,ks)+tv(i,j,kp))*0.5
      tvplot2(i,j,kkc)=tvplot2(i,j,kkc)+
     %tvstagg*dt
      uvplot2(i,j,kkc)=uvplot2(i,j,kkc)+
     %ustagg*vstagg*dt
      uwplot2(i,j,kkc)=uwplot2(i,j,kkc)+
     %ustagg*wstagg*dt
      vwplot2(i,j,kkc)=vwplot2(i,j,kkc)+
     %vstagg*wstagg*dt
  102 continue
  103 continue
!
!     Time averaged values and root mean squares
!     at some meridian sections
!
      do 111 jjc=1,nsez14
      js=jsez14(jjc)
      jp=jpv(js)
      jsy=jsym(js)
      jpsy=jsym(jp)
      do 110 k=kz1,kz2
      km=k-1
      do 110 i=ix1,ix2
      im=i-1
      vstagg=vo(i,js,k)
      ustagg=(uo(i,js,k)+uo(im,js,k)+
     %uo(i,jp,k)+uo(im,jp,k))*0.25
      wstagg=(wo(i,js,k)+wo(i,jp,k)+
     %wo(i,js,km)+wo(i,jp,km))*0.25
      xstagg=ustagg*cos(yv(js))-vstagg*sin(yv(js))
      ystagg=ustagg*sin(yv(js))+vstagg*cos(yv(js))
      vstagg=xstagg
      ustagg=ystagg
      vplot14(i,jjc,k)=vplot14(i,jjc,k)+
     %vstagg*dt
      uplot14(i,jjc,k)=uplot14(i,jjc,k)+
     %ustagg*dt
      wplot14(i,jjc,k)=wplot14(i,jjc,k)+
     %wstagg*dt
      vrmsplot14(i,jjc,k)=vrmsplot14(i,jjc,k)+
     %dt*(vstagg**2.)
      urmsplot14(i,jjc,k)=urmsplot14(i,jjc,k)+
     %dt*(ustagg**2.)
      wrmsplot14(i,jjc,k)=wrmsplot14(i,jjc,k)+
     %dt*(wstagg**2.)
      if(all(flagpo(i,js:jp,k,:).le.0)) then
	 tempo14(i,jjc,k)=tempo14(i,jjc,k)+dt
	 pstagg=(p(i,js,k)+p(i,jp,k))*0.5
	 pplot14(i,jjc,k)=pplot14(i,jjc,k)+
     %pstagg*dt
	 prmsplot14(i,jjc,k)=prmsplot14(i,jjc,k)+
     %dt*(pstagg**2.)
      endif
      vmod=sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
      vmodplot14(i,jjc,k)=vmodplot14(i,jjc,k)+
     %vmod*dt
      tvstagg=(tv(i,js,k)+tv(i,jp,k))*0.5
      tvplot14(i,jjc,k)=tvplot14(i,jjc,k)+
     %tvstagg*dt
      uvplot14(i,jjc,k)=uvplot14(i,jjc,k)+
     %ustagg*vstagg*dt
      uwplot14(i,jjc,k)=uwplot14(i,jjc,k)+
     %ustagg*wstagg*dt
      vwplot14(i,jjc,k)=vwplot14(i,jjc,k)+
     %vstagg*wstagg*dt
c
      vstagg=vo(i,jsy,k)
      ustagg=(uo(i,jsy,k)+uo(im,jsy,k)+
     %uo(i,jpsy,k)+uo(im,jpsy,k))*0.25
      wstagg=(wo(i,jsy,k)+wo(i,jpsy,k)+
     %wo(i,jsy,km)+wo(i,jpsy,km))*0.25
      xstagg=ustagg*cos(yv(jsy))-vstagg*sin(yv(jsy))
      ystagg=ustagg*sin(yv(jsy))+vstagg*cos(yv(jsy))
      vstagg=xstagg
      ustagg=ystagg
      vplot14(i,nsez14+jjc,k)=vplot14(i,nsez14+jjc,k)+
     %vstagg*dt
      uplot14(i,nsez14+jjc,k)=uplot14(i,nsez14+jjc,k)+
     %ustagg*dt
      wplot14(i,nsez14+jjc,k)=wplot14(i,nsez14+jjc,k)+
     %wstagg*dt
      vrmsplot14(i,nsez14+jjc,k)=vrmsplot14(i,nsez14+jjc,k)+
     %dt*(vstagg**2.)
      urmsplot14(i,nsez14+jjc,k)=urmsplot14(i,nsez14+jjc,k)+
     %dt*(ustagg**2.)
      wrmsplot14(i,nsez14+jjc,k)=wrmsplot14(i,nsez14+jjc,k)+
     %dt*(wstagg**2.)
      if(all(flagpo(i,jsy:jpsy,k,:).le.0)) then
	 tempo14(i,nsez14+jjc,k)=tempo14(i,nsez14+jjc,k)+dt
	 pstagg=(p(i,jsy,k)+p(i,jpsy,k))*0.5
	 pplot14(i,nsez14+jjc,k)=pplot14(i,nsez14+jjc,k)+
     %pstagg*dt
	 prmsplot14(i,nsez14+jjc,k)=prmsplot14(i,nsez14+jjc,k)+
     %dt*(pstagg**2.)
      endif
      vmod=sqrt(vstagg**2.+ustagg**2.+wstagg**2.)
      vmodplot14(i,nsez14+jjc,k)=vmodplot14(i,nsez14+jjc,k)+
     %vmod*dt
      tvstagg=(tv(i,jsy,k)+tv(i,jpsy,k))*0.5
      tvplot14(i,nsez14+jjc,k)=tvplot14(i,nsez14+jjc,k)+
     %tvstagg*dt
      uvplot14(i,nsez14+jjc,k)=uvplot14(i,nsez14+jjc,k)+
     %ustagg*vstagg*dt
      uwplot14(i,nsez14+jjc,k)=uwplot14(i,nsez14+jjc,k)+
     %ustagg*wstagg*dt
      vwplot14(i,nsez14+jjc,k)=vwplot14(i,nsez14+jjc,k)+
     %vstagg*wstagg*dt
  110 continue
  111 continue
!
!!!	 if((time-timep).gt.((periodo)*(iplant5+1))) then
c      if((abs(romega)*time_plot2*180./pig).gt.(rplot1*(iplant5+1)))
!	  call reynoldsstagg(1,nx,ny,nz,uo,vo,wo,dt)
	 call avgvort2staggcart(1,nx,ny,nz,dt,xc,xu,yc,yv,zc,uo,vo,wo)
!!!	 endif
!
      if(tempo3.ge.periodo2) then
c
	 do 105 iic=1,nsez1
	 do 104 k=kz1,kz2
	 do 104 j=jy1,jy2
	    vplot(iic,j,k)=vplot(iic,j,k)/tempo3
	    uplot(iic,j,k)=uplot(iic,j,k)/tempo3
	    wplot(iic,j,k)=wplot(iic,j,k)/tempo3
	    vrmsplot(iic,j,k)=
     %vrmsplot(iic,j,k)/tempo3-vplot(iic,j,k)**2.
	    urmsplot(iic,j,k)=
     %urmsplot(iic,j,k)/tempo3-uplot(iic,j,k)**2.
	    wrmsplot(iic,j,k)=
     %wrmsplot(iic,j,k)/tempo3-wplot(iic,j,k)**2.
	    if(tempo4(iic,j,k).ne.0.) then
	       pplot(iic,j,k)=pplot(iic,j,k)/tempo4(iic,j,k)
	       prmsplot(iic,j,k)=
     %prmsplot(iic,j,k)/tempo4(iic,j,k)-pplot(iic,j,k)**2.
	    endif
	    vmodplot(iic,j,k)=vmodplot(iic,j,k)/tempo3
	    tvplot(iic,j,k)=tvplot(iic,j,k)/tempo3
	    tvplot(iic,j,k)=tvplot(iic,j,k)/ru1
	    uvplot(iic,j,k)=
     %uvplot(iic,j,k)/tempo3-uplot(iic,j,k)*vplot(iic,j,k)
	    uwplot(iic,j,k)=
     %uwplot(iic,j,k)/tempo3-uplot(iic,j,k)*wplot(iic,j,k)
	    vwplot(iic,j,k)=
     %vwplot(iic,j,k)/tempo3-vplot(iic,j,k)*wplot(iic,j,k)
  104	 continue
  105	 continue
c
	 do 107 kkc=1,nsez2max
	 do 106 i=ix1,ix2
	 do 106 j=jy1,jy2
	    vplot2(i,j,kkc)=vplot2(i,j,kkc)/tempo3
	    uplot2(i,j,kkc)=uplot2(i,j,kkc)/tempo3
	    wplot2(i,j,kkc)=wplot2(i,j,kkc)/tempo3
	    vrmsplot2(i,j,kkc)=
     %vrmsplot2(i,j,kkc)/tempo3-vplot2(i,j,kkc)**2.
	    urmsplot2(i,j,kkc)=
     %urmsplot2(i,j,kkc)/tempo3-uplot2(i,j,kkc)**2.
	    wrmsplot2(i,j,kkc)=
     %wrmsplot2(i,j,kkc)/tempo3-wplot2(i,j,kkc)**2.
	    if(tempo5(i,j,kkc).ne.0.) then
	       pplot2(i,j,kkc)=pplot2(i,j,kkc)/tempo5(i,j,kkc)
	       prmsplot2(i,j,kkc)=
     %prmsplot2(i,j,kkc)/tempo5(i,j,kkc)-pplot2(i,j,kkc)**2.
	    endif
	    vmodplot2(i,j,kkc)=vmodplot2(i,j,kkc)/tempo3
	    tvplot2(i,j,kkc)=tvplot2(i,j,kkc)/tempo3
	    tvplot2(i,j,kkc)=tvplot2(i,j,kkc)/ru1
	    uvplot2(i,j,kkc)=
     %uvplot2(i,j,kkc)/tempo3-uplot2(i,j,kkc)*vplot2(i,j,kkc)
	    uwplot2(i,j,kkc)=
     %uwplot2(i,j,kkc)/tempo3-uplot2(i,j,kkc)*wplot2(i,j,kkc)
	    vwplot2(i,j,kkc)=
     %vwplot2(i,j,kkc)/tempo3-vplot2(i,j,kkc)*wplot2(i,j,kkc)
  106	 continue
  107	 continue
c
	 do 115 jjc=1,2*nsez14
	 do 114 k=kz1,kz2
	 do 114 i=ix1,ix2
	    vplot14(i,jjc,k)=vplot14(i,jjc,k)/tempo3
	    uplot14(i,jjc,k)=uplot14(i,jjc,k)/tempo3
	    wplot14(i,jjc,k)=wplot14(i,jjc,k)/tempo3
	    vrmsplot14(i,jjc,k)=
     %vrmsplot14(i,jjc,k)/tempo3-vplot14(i,jjc,k)**2.
	    urmsplot14(i,jjc,k)=
     %urmsplot14(i,jjc,k)/tempo3-uplot14(i,jjc,k)**2.
	    wrmsplot14(i,jjc,k)=
     %wrmsplot14(i,jjc,k)/tempo3-wplot14(i,jjc,k)**2.
	    if(tempo14(i,jjc,k).ne.0.) then
	       pplot14(i,jjc,k)=
     %pplot14(i,jjc,k)/tempo14(i,jjc,k)
	       prmsplot14(i,jjc,k)=
     %prmsplot14(i,jjc,k)/tempo14(i,jjc,k)-pplot14(i,jjc,k)**2.
	    endif
	    vmodplot14(i,jjc,k)=vmodplot14(i,jjc,k)/tempo3
	    tvplot14(i,jjc,k)=tvplot14(i,jjc,k)/tempo3
	    tvplot14(i,jjc,k)=tvplot14(i,jjc,k)/ru1
	    uvplot14(i,jjc,k)=
     %uvplot14(i,jjc,k)/tempo3-uplot14(i,jjc,k)*vplot14(i,jjc,k)
	    uwplot14(i,jjc,k)=
     %uwplot14(i,jjc,k)/tempo3-uplot14(i,jjc,k)*wplot14(i,jjc,k)
	    vwplot14(i,jjc,k)=
     %vwplot14(i,jjc,k)/tempo3-vplot14(i,jjc,k)*wplot14(i,jjc,k)
  114	 continue
  115	 continue
!
	 iplant4=iplant4+1
	 call kplotstagg(iplant4c,iplant4d,iplant4p,
     %nx,ny,nz,nzg,nbd,ntime,xu,xc,yv,yc,zwg,zcg,
     %flaguo,flagvo,flagwo)
!
ccc 13	    format(2x,a12,1x,e14.7,1x,e14.7,1x,e14.7)
ccc 14	    format(6x,a6,3x,e14.7,1x,e14.7,1x,e14.7)
ccc	    if((ibm.gt.1).and.(myrank.eq.0)) then
ccc	       do ibd=mbd,nbd
ccc		  open(83,file='ib_n.'//index(ibd)//'_'
ccc	%//char(48+iplant4c)//char(48+iplant4d)//char(48+iplant4p)//
ccc	%'_kplot1.stl')
ccc		  write(83,*)'solid  OBJECT'
ccc		  do j=lb(ibd)+1,lb(ibd)+mb(ibd)
ccc		     write(83,13)'facet normal',unvect(1,j),unvect(2,j),
ccc	%unvect(3,j)
ccc		     write(83,*)'outer loop'
ccc		     write(83,14)
ccc	%'vertex',vertex(1,1,j),vertex(2,1,j),vertex(3,1,j)
ccc		     write(83,14)
ccc	%'vertex',vertex(1,2,j),vertex(2,2,j),vertex(3,2,j)
ccc		     write(83,14)
ccc	%'vertex',vertex(1,3,j),vertex(2,3,j),vertex(3,3,j)
ccc		     write(83,*)'endloop'
ccc		     write(83,*)'endfacet'
ccc		  enddo
ccc		  write(83,*)'endsolid	OBJECT'
ccc		  close(83)
ccc	       enddo
ccc	    endif
!
!!!	    if((time-timep).gt.((periodo)*(iplant5+1))) then
c	  if((abs(romega)*time_plot2*180./pig).gt.
c     %(rplot1*(iplant5+1))) then
	    iplant5=iplant5+1
!	     call reynoldsstagg(2,nx,ny,nz,uo,vo,wo,dt)
	    call kplot2(iplant5c,iplant5d,iplant5p,
     %nx,ny,nz,nzg,nbd,ntime,flaguo,flagvo,flagwo)
	    call avgvort2staggcart(2,nx,ny,nz,dt,xc,xu,yc,yv,zc,uo,vo,wo)
	    call plotavgvort(iplant5c,iplant5d,iplant5p,
     %nx,ny,nz,nzg,nbd,ntime,flaguo,flagvo,flagwo)
c	     call vort2
c	     call plotvort(iplant5c,iplant5d,iplant5p)
!
ccc	    if((ibm.gt.1).and.(myrank.eq.0)) then
ccc	       do ibd=mbd,nbd
ccc		  open(83,file='ib_n.'//index(ibd)//'_'
ccc	%//char(48+iplant4c)//char(48+iplant4d)//char(48+iplant4p)//
ccc	%'_kplot2.stl')
ccc		  write(83,*)'solid  OBJECT'
ccc		  do j=lb(ibd)+1,lb(ibd)+mb(ibd)
ccc		     write(83,13)'facet normal',unvect(1,j),unvect(2,j),
ccc	%unvect(3,j)
ccc		     write(83,*)'outer loop'
ccc		     write(83,14)
ccc	%'vertex',vertex(1,1,j),vertex(2,1,j),vertex(3,1,j)
ccc		     write(83,14)
ccc	%'vertex',vertex(1,2,j),vertex(2,2,j),vertex(3,2,j)
ccc		     write(83,14)
ccc	%'vertex',vertex(1,3,j),vertex(2,3,j),vertex(3,3,j)
ccc		     write(83,*)'endloop'
ccc		     write(83,*)'endfacet'
ccc		  enddo
ccc		  write(83,*)'endsolid	OBJECT'
ccc		  close(83)
ccc	       enddo
ccc	    endif
!
!!!	    endif
!
	 do iic=1,nsez1
	    do k=kz1,kz2
	       do j=jy1,jy2
		  vrmsplot(iic,j,k)=0.
		  urmsplot(iic,j,k)=0.
		  wrmsplot(iic,j,k)=0.
		  prmsplot(iic,j,k)=0.
		  vplot(iic,j,k)=0.
		  uplot(iic,j,k)=0.
		  wplot(iic,j,k)=0.
		  pplot(iic,j,k)=0.
		  tempo4(iic,j,k)=0.
		  uvplot(iic,j,k)=0.
		  vwplot(iic,j,k)=0.
		  uwplot(iic,j,k)=0.
		  voraz1med(iic,j,k)=0.
		  vorr1med(iic,j,k)=0.
		  vorz1med(iic,j,k)=0.
		  vmodplot(iic,j,k)=0.
		  tvplot(iic,j,k)=0.
	       enddo
	    enddo
	 enddo
c
	 do kkc=1,nsez2max
	    do i=ix1,ix2
	       do j=jy1,jy2
		  vrmsplot2(i,j,kkc)=0.
		  urmsplot2(i,j,kkc)=0.
		  wrmsplot2(i,j,kkc)=0.
		  prmsplot2(i,j,kkc)=0.
		  vplot2(i,j,kkc)=0.
		  uplot2(i,j,kkc)=0.
		  wplot2(i,j,kkc)=0.
		  pplot2(i,j,kkc)=0.
		  tempo5(i,j,kkc)=0.
		  uvplot2(i,j,kkc)=0.
		  vwplot2(i,j,kkc)=0.
		  uwplot2(i,j,kkc)=0.
		  voraz2med(i,j,kkc)=0.
		  vorr2med(i,j,kkc)=0.
		  vorz2med(i,j,kkc)=0.
		  vmodplot2(i,j,kkc)=0.
		  tvplot2(i,j,kkc)=0.
	       enddo
	    enddo
	 enddo
c
	 do jjc=1,2*nsez14
	    do k=kz1,kz2
	       do i=ix1,ix2
		  vrmsplot14(i,jjc,k)=0.
		  urmsplot14(i,jjc,k)=0.
		  wrmsplot14(i,jjc,k)=0.
		  prmsplot14(i,jjc,k)=0.
		  vplot14(i,jjc,k)=0.
		  uplot14(i,jjc,k)=0.
		  wplot14(i,jjc,k)=0.
		  pplot14(i,jjc,k)=0.
		  tempo14(i,jjc,k)=0.
		  uvplot14(i,jjc,k)=0.
		  vwplot14(i,jjc,k)=0.
		  uwplot14(i,jjc,k)=0.
		  voraz14med(i,jjc,k)=0.
		  vorr14med(i,jjc,k)=0.
		  vorz14med(i,jjc,k)=0.
		  vmodplot14(i,jjc,k)=0.
		  tvplot14(i,jjc,k)=0.
	       enddo
	    enddo
	 enddo
	 tempo3=0.
      endif
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Subroutine that evaluates the time-averaged Cartesian components
c of the vorticity on a cylindrical grid
c
      subroutine avgvort2staggcart
     %(ind,nx,ny,nz,dt,xc,xu,yc,yv,zc,uo,vo,wo)
c
      use sections
      use vortavg
      include'common.h'
      include'averages.h'
c
c Global variables
      integer ind,nx,ny,nz
      real yc(ny),yv(ny)
      real dt,xc(nx),xu(nx),zc(nz),
     %uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
c
c Local variables
      integer iic,ic,ip,im,jjc,jc,jp,jm,jsy,jpsy,kkc,kc,kp,km,
     %i,j,k
      real dtheta,dr,dz,dq1x3,dq3x1,dq3x2,dq2x3,dq1x2,dq2x1
      real xstagg,ystagg
c
      if(ind.eq.2) goto 200
c
c    ! Time-averaged azimuthal vorticity - Circumferential sections !
      do 101 kc=kz1,kz2
	kp=kc+1
	km=kc-1
	dz=zc(kp)-zc(km)
	do 101 iic=1,nsez1
	  ic=isez1(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 101 jc=jy1,jy2
	    dq1x3=(uo(ic,jc,kp)
     %-uo(ic,jc,km))/dz
	    dq3x1=(0.5*(wo(ip,jc,kc)+wo(ip,jc,km))
     %-0.5*(wo(ic,jc,kc)+wo(ic,jc,km)))/dr
	    voraz1med(iic,jc,kc)=
     %voraz1med(iic,jc,kc)+(dq1x3-dq3x1)*dt
  101 continue
c
c    ! Time-averaged azimuthal vorticity - Cross sections !
      do 111 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
	dz=zc(kp)-zc(kc)
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
	  voraz2med(ic,jc,kkc)=
     %voraz2med(ic,jc,kkc)+(dq1x3-dq3x1)*dt
  112 continue
	do 113 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 113 jc=jy1,jy2
	    dq1x3=(0.5*(uo(ic,jc,kp)+uo(im,jc,kp))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dz
	    dq3x1=(wo(ip,jc,kc)
     %-wo(im,jc,kc))/dr
	    voraz2med(ic,jc,kkc)=
     %voraz2med(ic,jc,kkc)+(dq1x3-dq3x1)*dt
  113 continue
  111 continue
c
c    ! Time-averaged azimuthal vorticity - Meridian sections !
      do 121 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jsy=jsym(jc)
	jpsy=jsym(jp)
	do 121 kc=kz1,kz2
	  kp=kc+1
	  km=kc-1
	  dz=zc(kp)-zc(km)
	  ic=ix1
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  dq1x3=(0.25*(uo(ic,jc,kp)+uo(im,jc,kp)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,km)+uo(im,jc,km)
     %+uo(ic,jp,km)+uo(im,jp,km)))/dz
	  dq3x1=(0.25*(wo(ip,jc,kc)+wo(ip,jp,kc)
     %+wo(ip,jc,km)+wo(ip,jp,km))
     %-0.25*(wo(ic,jsy,kc)+wo(ic,jpsy,kc)
     %+wo(ic,jsy,km)+wo(ic,jpsy,km)))/dr
	  voraz14med(ic,jjc,kc)=
     %voraz14med(ic,jjc,kc)+(dq1x3-dq3x1)*dt
	  dq1x3=(0.25*(uo(ic,jsy,kp)+uo(im,jsy,kp)
     %+uo(ic,jpsy,kp)+uo(im,jpsy,kp))
     %-0.25*(uo(ic,jsy,km)+uo(im,jsy,km)
     %+uo(ic,jpsy,km)+uo(im,jpsy,km)))/dz
	  dq3x1=(0.25*(wo(ip,jsy,kc)+wo(ip,jpsy,kc)
     %+wo(ip,jsy,km)+wo(ip,jpsy,km))
     %-0.25*(wo(ic,jc,kc)+wo(ic,jp,kc)
     %+wo(ic,jc,km)+wo(ic,jp,km)))/dr
	  voraz14med(ic,nsez14+jjc,kc)=
     %voraz14med(ic,nsez14+jjc,kc)+(dq1x3-dq3x1)*dt
	  do 121 ic=ix1+1,ix2
	    ip=ic+1
	    im=ic-1
	    dr=xc(ip)-xc(im)
	    dq1x3=(0.25*(uo(ic,jc,kp)+uo(im,jc,kp)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jc,km)+uo(im,jc,km)
     %+uo(ic,jp,km)+uo(im,jp,km)))/dz
	    dq3x1=(0.25*(wo(ip,jc,kc)+wo(ip,jp,kc)
     %+wo(ip,jc,km)+wo(ip,jp,km))
     %-0.25*(wo(im,jc,kc)+wo(im,jp,kc)
     %+wo(im,jc,km)+wo(im,jp,km)))/dr
	    voraz14med(ic,jjc,kc)=
     %voraz14med(ic,jjc,kc)+(dq1x3-dq3x1)*dt
	    dq1x3=(0.25*(uo(ic,jsy,kp)+uo(im,jsy,kp)
     %+uo(ic,jpsy,kp)+uo(im,jpsy,kp))
     %-0.25*(uo(ic,jsy,km)+uo(im,jsy,km)
     %+uo(ic,jpsy,km)+uo(im,jpsy,km)))/dz
	    dq3x1=(0.25*(wo(ip,jsy,kc)+wo(ip,jpsy,kc)
     %+wo(ip,jsy,km)+wo(ip,jpsy,km))
     %-0.25*(wo(im,jsy,kc)+wo(im,jpsy,kc)
     %+wo(im,jsy,km)+wo(im,jpsy,km)))/dr
	    voraz14med(ic,nsez14+jjc,kc)=
     %voraz14med(ic,nsez14+jjc,kc)+(dq1x3-dq3x1)*dt
  121 continue
c
c    ! Time-averaged radial vorticity - Circumferential sections !
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
	    dq3x2=(0.25*(wo(ic,jp,kc)+wo(ip,jp,kc)
     %+wo(ic,jp,km)+wo(ip,jp,km))
     %-0.25*(wo(ic,jm,kc)+wo(ip,jm,kc)
     %+wo(ic,jm,km)+wo(ip,jm,km)))/dtheta
	    dq2x3=(0.25*(vo(ic,jc,kp)+vo(ip,jc,kp)
     %+vo(ic,jm,kp)+vo(ip,jm,kp))
     %-0.25*(vo(ic,jc,km)+vo(ip,jc,km)
     %+vo(ic,jm,km)+vo(ip,jm,km)))/dz
	    vorr1med(iic,jc,kc)=
     %vorr1med(iic,jc,kc)+(dq3x2/xu(ic)-dq2x3)*dt
  131 continue
c
c    ! Time-averaged radial vorticity - Cross sections !
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
     %-wo(ic,jm,kc))/dtheta
	    dq2x3=(0.5*(vo(ic,jc,kp)+vo(ic,jm,kp))
     %-0.5*(vo(ic,jc,kc)+vo(ic,jm,kc)))/dz
	    vorr2med(ic,jc,kkc)=
     %vorr2med(ic,jc,kkc)+(dq3x2/xc(ic)-dq2x3)*dt
  141 continue
c
c    ! Time-averaged radial vorticity - Meridian sections !
      do 151 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jsy=jsym(jc)
	jpsy=jsym(jp)
	dtheta=dely
	do 151 kc=kz1,kz2
	  kp=kc+1
	  km=kc-1
	  dz=zc(kp)-zc(km)
	  do 151 ic=ix1,ix2
	    dq3x2=(0.5*(wo(ic,jp,kc)+wo(ic,jp,km))
     %-0.5*(wo(ic,jc,kc)+wo(ic,jc,km)))/dtheta
	    dq2x3=(vo(ic,jc,kp)
     %-vo(ic,jc,km))/dz
	    vorr14med(ic,jjc,kc)=
     %vorr14med(ic,jjc,kc)+(dq3x2/xc(ic)-dq2x3)*dt
	    dq3x2=(0.5*(wo(ic,jpsy,kc)+wo(ic,jpsy,km))
     %-0.5*(wo(ic,jsy,kc)+wo(ic,jsy,km)))/dtheta
	    dq2x3=(vo(ic,jsy,kp)
     %-vo(ic,jsy,km))/dz
	    vorr14med(ic,nsez14+jjc,kc)=
     %vorr14med(ic,nsez14+jjc,kc)+(dq3x2/xc(ic)-dq2x3)*dt
  151 continue
c
c    ! Time-averaged axial vorticity - Circumferential sections !
      do 161 kc=kz1,kz2
	do 161 iic=1,nsez1
	  ic=isez1(iic)
	  ip=ic+1
	  dr=xc(ip)-xc(ic)
	  do 161 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.5*(vo(ip,jc,kc)+vo(ip,jm,kc))*xc(ip)
     %-0.5*(vo(ic,jc,kc)+vo(ic,jm,kc))*xc(ic))/dr
	    dq1x2=(uo(ic,jp,kc)
     %-uo(ic,jm,kc))/dtheta
	    vorz1med(iic,jc,kc)=
     %vorz1med(iic,jc,kc)+((dq2x1-dq1x2)/xu(ic))*dt
  161 continue
c
c    ! Time-averaged axial vorticity - Cross sections !
      do 171 kkc=1,nsez2max
	kc=ksez2(kkc)
	kp=kc+1
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
     %+vo(ic,jm,kp)+vo(ip,jm,kp))*xu(ic))/dr
	  dq1x2=(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jm,kc)+uo(im,jm,kc)
     %+uo(ic,jm,kp)+uo(im,jm,kp)))/dtheta
	  vorz2med(ic,jc,kkc)=
     %vorz2med(ic,jc,kkc)+((dq2x1-dq1x2)/xc(ic))*dt
  172 continue
c
	do 173 ic=ix1+1,ix2
	  ip=ic+1
	  im=ic-1
	  dr=xc(ip)-xc(im)
	  do 173 jc=jy1,jy2
	    jp=jpv(jc)
	    jm=jmv(jc)
	    dtheta=2.*dely
	    dq2x1=(0.25*(vo(ip,jc,kc)+vo(ip,jm,kc)
     %+vo(ip,jc,kp)+vo(ip,jm,kp))*xc(ip)
     %-0.25*(vo(im,jc,kc)+vo(im,jm,kc)
     %+vo(im,jc,kp)+vo(im,jm,kp))*xc(im))/dr
	    dq1x2=
     %(0.25*(uo(ic,jp,kc)+uo(im,jp,kc)
     %+uo(ic,jp,kp)+uo(im,jp,kp))
     %-0.25*(uo(ic,jm,kc)+uo(im,jm,kc)
     %+uo(ic,jm,kp)+uo(im,jm,kp)))/dtheta
	    vorz2med(ic,jc,kkc)=
     %vorz2med(ic,jc,kkc)+((dq2x1-dq1x2)/xc(ic))*dt
  173 continue
  171 continue
c
c    ! Time-averaged axial vorticity - Meridian sections !
      do 181 jjc=1,nsez14
	jc=jsez14(jjc)
	jp=jpv(jc)
	jsy=jsym(jc)
	jpsy=jsym(jp)
	dtheta=dely
	do 181 kc=kz1,kz2
	  ic=ix1
	  im=ic-1
	  ip=ic+1
	  dr=xu(ic)-xu(im)
	  dq2x1=(0.5*(vo(ic,jc,kc)+vo(ip,jc,kc))*xu(ic))/dr
	  dq1x2=(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dtheta
	  vorz14med(ic,jjc,kc)=
     %vorz14med(ic,jjc,kc)+((dq2x1-dq1x2)/xc(ic))*dt
	  dq2x1=(0.5*(vo(ic,jsy,kc)+vo(ip,jsy,kc))*xu(ic))/dr
	  dq1x2=(0.5*(uo(ic,jpsy,kc)+uo(im,jpsy,kc))
     %-0.5*(uo(ic,jsy,kc)+uo(im,jsy,kc)))/dtheta
	  vorz14med(ic,nsez14+jjc,kc)=
     %vorz14med(ic,nsez14+jjc,kc)+((dq2x1-dq1x2)/xc(ic))*dt
c
	  do 181 ic=ix1+1,ix2
	    ip=ic+1
	    im=ic-1
	    dr=xc(ip)-xc(im)
	    dq2x1=(vo(ip,jc,kc)*xc(ip)
     %-vo(im,jc,kc)*xc(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jp,kc)+uo(im,jp,kc))
     %-0.5*(uo(ic,jc,kc)+uo(im,jc,kc)))/dtheta
	    vorz14med(ic,jjc,kc)=
     %vorz14med(ic,jjc,kc)+((dq2x1-dq1x2)/xc(ic))*dt
	    dq2x1=(vo(ip,jsy,kc)*xc(ip)
     %-vo(im,jsy,kc)*xc(im))/dr
	    dq1x2=
     %(0.5*(uo(ic,jpsy,kc)+uo(im,jpsy,kc))
     %-0.5*(uo(ic,jsy,kc)+uo(im,jsy,kc)))/dtheta
	    vorz14med(ic,nsez14+jjc,kc)=
     %vorz14med(ic,nsez14+jjc,kc)+((dq2x1-dq1x2)/xc(ic))*dt
  181 continue
      return
c
  200 continue
c
!
!     Circumferential sections
!
      do 205 iic=1,nsez1
      do 204 k=kz1,kz2
      do 204 j=jy1,jy2
	 voraz1med(iic,j,k)=voraz1med(iic,j,k)/tempo3
	 vorr1med(iic,j,k)=vorr1med(iic,j,k)/tempo3
	 vorz1med(iic,j,k)=vorz1med(iic,j,k)/tempo3
	 xstagg=
     %vorr1med(iic,j,k)*cos(yc(j))-voraz1med(iic,j,k)*sin(yc(j))
	 ystagg=
     %vorr1med(iic,j,k)*sin(yc(j))+voraz1med(iic,j,k)*cos(yc(j))
	 voraz1med(iic,j,k)=xstagg
	 vorr1med(iic,j,k)=ystagg
  204 continue
  205 continue
c
!
!     Cross sections
!
      do 207 kkc=1,nsez2max
      do 206 i=ix1,ix2
      do 206 j=jy1,jy2
	 voraz2med(i,j,kkc)=voraz2med(i,j,kkc)/tempo3
	 vorr2med(i,j,kkc)=vorr2med(i,j,kkc)/tempo3
	 vorz2med(i,j,kkc)=vorz2med(i,j,kkc)/tempo3
	 xstagg=
     %vorr2med(i,j,kkc)*cos(yc(j))-voraz2med(i,j,kkc)*sin(yc(j))
	 ystagg=
     %vorr2med(i,j,kkc)*sin(yc(j))+voraz2med(i,j,kkc)*cos(yc(j))
	 voraz2med(i,j,kkc)=xstagg
	 vorr2med(i,j,kkc)=ystagg
  206 continue
  207 continue
c
!
!     Meridian sections
!
      do 215 jjc=1,nsez14
      jc=jsez14(jjc)
      jsy=jsym(jc)
      do 214 k=kz1,kz2
      do 214 i=ix1,ix2
	 voraz14med(i,jjc,k)=voraz14med(i,jjc,k)/tempo3
	 vorr14med(i,jjc,k)=vorr14med(i,jjc,k)/tempo3
	 vorz14med(i,jjc,k)=vorz14med(i,jjc,k)/tempo3
	 xstagg=
     %vorr14med(i,jjc,k)*cos(yv(jc))
     %-voraz14med(i,jjc,k)*sin(yv(jc))
	 ystagg=
     %vorr14med(i,jjc,k)*sin(yv(jc))
     %+voraz14med(i,jjc,k)*cos(yv(jc))
	 voraz14med(i,jjc,k)=xstagg
	 vorr14med(i,jjc,k)=ystagg
	 voraz14med(i,nsez14+jjc,k)=voraz14med(i,nsez14+jjc,k)/tempo3
	 vorr14med(i,nsez14+jjc,k)=vorr14med(i,nsez14+jjc,k)/tempo3
	 vorz14med(i,nsez14+jjc,k)=vorz14med(i,nsez14+jjc,k)/tempo3
	 xstagg=
     %vorr14med(i,nsez14+jjc,k)*cos(yv(jsy))
     %-voraz14med(i,nsez14+jjc,k)*sin(yv(jsy))
	 ystagg=
     %vorr14med(i,nsez14+jjc,k)*sin(yv(jsy))
     %+voraz14med(i,nsez14+jjc,k)*cos(yv(jsy))
	 voraz14med(i,nsez14+jjc,k)=xstagg
	 vorr14med(i,nsez14+jjc,k)=ystagg
  214 continue
  215 continue
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sondestaggcart(nx,ny,nz,nbd,time,dt,yv,vo,uo,wo,p,flagpo)
c
      use points
      implicit none
      include 'averages.h'
c
c Global variables
      integer nx,ny,nz,nbd
      integer flagpo(nx,ny,nz,nbd)
      real time,dt
      real yv(ny)
      real vo(nx,ny,nz),uo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
c
c Local variables
      integer ii,i,ic,ip,jc,jp,kc,kp
      real vstagg,ustagg,wstagg,pstagg,xstagg,ystagg
c
 10   format(5f25.8)
 11   format(9f25.8)
!
!     Instantaneous values, time averages and root mean squares
!
      timepoints=timepoints+dt
      do ii=1,nprbmax2
	 i=prbindx2+(ii-1)
	 ic=iprb2(ii)
	 ip=ic+1
	 jc=jprb2(ii)
	 jp=jpv(jc)
	 kc=kprb2(ii)
	 kp=kc+1
	 vstagg=0.25*(vo(ic,jc,kc)+vo(ip,jc,kc)+
     %vo(ic,jc,kp)+vo(ip,jc,kp))
	 ustagg=0.25*(uo(ic,jc,kc)+uo(ic,jp,kc)+
     %uo(ic,jc,kp)+uo(ic,jp,kp))
	 wstagg=0.25*(wo(ic,jc,kc)+wo(ip,jc,kc)+
     %wo(ic,jp,kc)+wo(ip,jp,kc))
	 pstagg=0.125*(p(ic,jc,kc)+p(ip,jc,kc)+
     %p(ic,jp,kc)+p(ip,jp,kc)+
     %p(ic,jc,kp)+p(ip,jc,kp)+
     %p(ic,jp,kp)+p(ip,jp,kp))
	 xstagg=ustagg*cos(yv(jc))-vstagg*sin(yv(jc))
	 ystagg=ustagg*sin(yv(jc))+vstagg*cos(yv(jc))
	 vstagg=xstagg
	 ustagg=ystagg
         if((time-timep).le.periodo) goto 100
	 if(any(flagpo(ic:ip,jc:jp,kc:kp,:).ne.0)) goto 100
	 avrpoints1(ii,1)=avrpoints1(ii,1)+dt*vstagg
	 avrpoints1(ii,2)=avrpoints1(ii,2)+dt*ustagg
	 avrpoints1(ii,3)=avrpoints1(ii,3)+dt*wstagg
	 avrpoints1(ii,4)=avrpoints1(ii,4)+dt*pstagg
	 rmspoints1(ii,1)=rmspoints1(ii,1)+dt*(vstagg**2.)
	 rmspoints1(ii,2)=rmspoints1(ii,2)+dt*(ustagg**2.)
	 rmspoints1(ii,3)=rmspoints1(ii,3)+dt*(wstagg**2.)
	 rmspoints1(ii,4)=rmspoints1(ii,4)+dt*(pstagg**2.)
	 timepoints1(ii)=timepoints1(ii)+dt
 100	 continue
	 write(700+i,10)time,vstagg,ustagg,wstagg,pstagg   !*ros
      enddo
c
      if(timepoints.ge.periodo3) then
!
!     Time averages and root mean squares on file
!
	 do ii=1,nprbmax2
	    i=prbindx2+(ii-1)
	    avrpoints1(ii,1)=avrpoints1(ii,1)/timepoints1(ii)
	    avrpoints1(ii,2)=avrpoints1(ii,2)/timepoints1(ii)
	    avrpoints1(ii,3)=avrpoints1(ii,3)/timepoints1(ii)
	    avrpoints1(ii,4)=avrpoints1(ii,4)/timepoints1(ii)
	    rmspoints1(ii,1)=
     %sqrt(rmspoints1(ii,1)/timepoints1(ii)-avrpoints1(ii,1)**2.)
	    rmspoints1(ii,2)=
     %sqrt(rmspoints1(ii,2)/timepoints1(ii)-avrpoints1(ii,2)**2.)
	    rmspoints1(ii,3)=
     %sqrt(rmspoints1(ii,3)/timepoints1(ii)-avrpoints1(ii,3)**2.)
	    rmspoints1(ii,4)=
     %sqrt(rmspoints1(ii,4)/timepoints1(ii)-avrpoints1(ii,4)**2.)
	    write(800+i,11)time-periodo3*0.5,
     %avrpoints1(ii,1),rmspoints1(ii,1),
     %avrpoints1(ii,2),rmspoints1(ii,2),
     %avrpoints1(ii,3),rmspoints1(ii,3),
     %avrpoints1(ii,4),rmspoints1(ii,4)	   !*ros     !*ros
	    avrpoints1(ii,1)=0.
	    avrpoints1(ii,2)=0.
	    avrpoints1(ii,3)=0.
	    avrpoints1(ii,4)=0.
	    rmspoints1(ii,1)=0.
	    rmspoints1(ii,2)=0.
	    rmspoints1(ii,3)=0.
	    rmspoints1(ii,4)=0.
	    timepoints1(ii)=0.
	 enddo
	 timepoints=0.
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine plotinststagg_d(iplant9c,iplant9d,iplant9p,
     %nx,ny,nz,nzg,nbd,ntime,xc,xu,yc,yv,zcg,zwg,p,vo,uo,wo,dens)
!     %flaguo,flagvo,flagwo)
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
      real p(nx,ny,nz),vo(nx,ny,nz),uo(nx,ny,nz),wo(nx,ny,nz),dens(nx,ny,nz)
c
c Local variables
      integer iic,iicd,iicu,is,im,ip,kkc,kkcd,kkcu,ks,km,kp,
     %jjc,jjcd,jjcu,js,jm,jp,jsy,jpsy,i,j,k,l,ir,jr,kr,lr,
     %i1,jjp,kk1,kk2,kk,kkd,kku
      integer STATUS(MPI_STATUS_SIZE)
!      integer flaguor(ny,nzg,nbd),
!     %flagvor2(nx,nzg,nbd),flagvor2sym(nx,nzg,nbd)
      real pr(ny,nzg),vor(ny,nzg),uor(ny,nzg),wor(ny,nzg),densr(ny,nzg),
     %pr2(nx,nzg),pr2sym(nx,nzg),vor2(nx,nzg),vor2sym(nx,nzg),
     %uor2(nx,nzg),uor2sym(nx,nzg),wor2(nx,nzg),wor2sym(nx,nzg),
     %densr2(nx,nzg),densr2sym(nx,nzg)
c
c Parameters
      integer lamb2,nxp,nyp,nzp
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
      if(iplant9.eq.1) then
	 if(myrank.eq.0) then
	    do iic=1,nsez1
	       iicd=iic/10
	       iicu=mod(iic,10)
	       open
     %(81,file='grid1_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'_plotinst.xyz')
	       open
     %(82,file='grid1_new_2D_slice'//char(48+iicd)
     %//char(48+iicu)//'_plotinst.xyz')
	       is=isez1(iic)
	       write(81,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	       write(82,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	       do 1 k=2,nzg-1,nzp
	       do 10 l=1,lamb2
	       do 10 j=jy1,jy2,nyp
	       write(81,9)xu(is)*cos(yc(jy1))
   10	       write(82,9)xu(is)*cos(yc(j)+(l-1)*2.*pi/lamb2)
	       write(81,9)xu(is)*cos(yc(jy1))
   1	       write(82,9)xu(is)*cos(yc(jy1)+lamb2*deltasect)
	       do 2 k=2,nzg-1,nzp
	       do 20 l=1,lamb2
	       do 20 j=jy1,jy2,nyp
	       write(81,9)
     %xu(is)*sin(yc(jy1))+xu(is)*(yc(j)+(l-1)*2.*pi/lamb2)
   20	       write(82,9)xu(is)*sin(yc(j)+(l-1)*2.*pi/lamb2)
	       write(81,9)
     %xu(is)*sin(yc(jy1))+xu(is)*deltasect*lamb2
   2	       write(82,9)xu(is)*sin(yc(jy1)+lamb2*deltasect)
	       do 3 k=2,nzg-1,nzp
	       do 30 l=1,lamb2
	       do 30 j=jy1,jy2,nyp
	       write(81,9)zcg(k)
   30	       write(82,9)zcg(k)
	       write(81,9)zcg(k)
   3	       write(82,9)zcg(k)
	       close(81)
	       close(82)
	    enddo
	    do kkc=1,nsez2
	       kkcd=kkc/10
	       kkcu=mod(kkc,10)
	       open(81,file='grid2_2D_slice'//char(48+kkcd)
     %//char(48+kkcu)//'_plotinst.xyz')
	       ks=ksez2g(kkc)
	       write(81,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	       do 11 i=ix1,ix2,nxp
	       do 110 l=1,lamb2
	       do 110 j=jy1,jy2,nyp
  110	       write(81,9)xc(i)*cos(yc(j)+(l-1)*2.*pi/lamb2)
  11	       write(81,9)xc(i)*cos(yc(jy1)+lamb2*deltasect)
	       do 12 i=ix1,ix2,nxp
	       do 120 l=1,lamb2
	       do 120 j=jy1,jy2,nyp
  120	       write(81,9)xc(i)*sin(yc(j)+(l-1)*2.*pi/lamb2)
  12	       write(81,9)xc(i)*sin(yc(jy1)+lamb2*deltasect)
	       do 13 i=ix1,ix2,nxp
	       do 130 l=1,lamb2
	       do 130 j=jy1,jy2,nyp
  130	       write(81,9)zwg(ks)
  13	       write(81,9)zwg(ks)
	       close(81)
	    enddo
	    do jjc=1,nsez14
	       jjcd=jjc/10
	       jjcu=mod(jjc,10)
	       open(81,file='grid3_2D_slice'//char(48+jjcd)
     %//char(48+jjcu)//'_plotinst.xyz')
	       js=jsez14(jjc)
	       jsy=jsym(js)
	       write(81,101)1,2*(nx-2)/nxp,(nzg-2)/nzp
	       do 210 k=2,nzg-1,nzp
	       do 211 i=ix2,ix1,-nxp
  211	       write(81,9)xc(i)*cos(yv(jsy))
	       do 212 i=ix1,ix2,nxp
  212	       write(81,9)xc(i)*cos(yv(js))
  210	       continue
	       do 220 k=2,nzg-1,nzp
	       do 221 i=ix2,ix1,-nxp
  221	       write(81,9)xc(i)*sin(yv(jsy))
	       do 222 i=ix1,ix2,nxp
  222	       write(81,9)xc(i)*sin(yv(js))
  220	       continue
	       do 230 k=2,nzg-1,nzp
	       do 230 i=1,2*(nx-2),nxp
  230	       write(81,9)zcg(k)
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
	    write(82,101)lamb2*(ny-2)/nyp+1,1,(nzg-2)/nzp
	    i1=1
	    write(82,102)ntime,i1,i1,i1
!
!Density
!
	    do 40 k=kz1,kz2,nzp
	    do 41 l=1,lamb2
	    do 41 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)0.5*(dens(is,j,k)+dens(ip,j,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
 41	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)0.5*(dens(is,jy1,k)+dens(ip,jy1,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
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
     %(densr(jy1:jy2,kk1:kk2),1*(ny-2)*(nz-2),MTYPE,jjp,21,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 42 kr=kk1,kk2,nzp
		  do 43 lr=1,lamb2
		  do 43 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)densr(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 43		  continue
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,9)densr(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 42		  continue
	       enddo
	    endif
!	    write(82,*)	  !!!!!!
!	    write(82,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,1),k=1,(nzg-2)/nzp)
!
!Azimuthal velocity
!
	    do 50 k=kz1,kz2,nzp
	    do 51 l=1,lamb2
	    do 51 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       jm=jmv(j)
	       write(82,9)
     %0.25*(vo(is,j,k)+vo(ip,j,k)+vo(is,jm,k)+vo(ip,jm,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
   51	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)
     %0.25*(vo(is,jy1,k)+vo(ip,jy1,k)+vo(is,jy2,k)+vo(ip,jy2,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
   50	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 52 kr=kk1,kk2,nzp
		  do 53 lr=1,lamb2
		  do 53 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)vor(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   53		  continue
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,9)vor(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   52		  continue
	       enddo
	    endif
!
!Radial velocity
!
	    do 60 k=kz1,kz2,nzp
	    do 61 l=1,lamb2
	    do 61 j=jy1,jy2,nyp
c	    if(all(flaguo(is,j,k,:)).le.0) then
	       write(82,9)uo(is,j,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
   61	    continue
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)uo(is,jy1,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
   60	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 62 kr=kk1,kk2,nzp
		  do 63 lr=1,lamb2
		  do 63 jr=jy1,jy2,nyp
c		  if(all(flaguor(jr,kr,:)).le.0) then
		     write(82,9)uor(jr,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
   63		  continue
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,9)uor(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
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
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,9)
     %0.25*(wo(is,jy1,k)+wo(ip,jy1,k)+wo(is,jy1,km)+wo(ip,jy1,km))
c	    else
c	       write(82,9)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,9)wor(jy1,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
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
c	    if(all(flaguo(is,jy1,k,:)).le.0) then
	       write(82,8)0.5*(p(is,jy1,k)+p(ip,jy1,k))	  !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
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
c		  if(all(flaguor(jy1,kr,:)).le.0) then
		     write(82,8)pr(jy1,kr)   !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
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
     %(0.5*(dens(is,jy1:jy2,kz1:kz2)+dens(ip,jy1:jy2,kz1:kz2)),
     %1*(ny-2)*(nz-2),MTYPE,0,21,
     %MPI_COMM_EDDY,STATUS,IERR)
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
	 write(82,101)lamb2*(ny-2)/nyp+1,(nx-2)/nxp,1
	 i1=1
	 write(82,102)ntime,i1,i1,i1
!
!Density
!
	 do 140 i=ix1,ix2,nxp
	 do 141 l=1,lamb2
	 do 141 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)0.5*(dens(i,j,ks)+dens(i,j,kp))
c	 else
c	    write(82,9)1.e+03
c	 endif
 141	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)0.5*(dens(i,jy1,ks)+dens(i,jy1,kp))
c	 else
c	    write(82,9)1.e+03
c	 endif
 140	 continue
!	 write(82,*)	!!!!!!
!	 write(82,103)
!    %(((i1,j=1,lamb2*(ny-2)/nyp+1),i=1,(nx-2)/nxp),k=1,1)
!
!Azimuthal velocity
!
	 do 150 i=ix1,ix2,nxp
	 do 151 l=1,lamb2
	 do 151 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    jm=jmv(j)
	    write(82,9)
     %0.25*(vo(i,j,ks)+vo(i,jm,ks)+vo(i,j,kp)+vo(i,jm,kp))
c	 else
c	    write(82,9)1.e+03
c	 endif
 151	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)
     %0.25*(vo(i,jy1,ks)+vo(i,jy2,ks)+vo(i,jy1,kp)+vo(i,jy2,kp))
c	 else
c	    write(82,9)1.e+03
c	 endif
 150	 continue
!
!Radial velocity
!
	 do 160 i=ix1,ix2,nxp
	 im=i-1
	 do 161 l=1,lamb2
	 do 161 j=jy1,jy2,nyp
c	 if(all(flagwo(i,j,ks,:)).le.0) then
	    write(82,9)
     %0.25*(uo(i,j,ks)+uo(im,j,ks)+uo(i,j,kp)+uo(im,j,kp))
c	 else
c	    write(82,9)1.e+03
c	 endif
 161	 continue
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)
     %0.25*(uo(i,jy1,ks)+uo(im,jy1,ks)+uo(i,jy1,kp)+uo(im,jy1,kp))
c	 else
c	    write(82,9)1.e+03
c	 endif
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
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,9)wo(i,jy1,ks)
c	 else
c	    write(82,9)1.e+03
c	 endif
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
c	 if(all(flagwo(i,jy1,ks,:)).le.0) then
	    write(82,8)0.5*(p(i,jy1,ks)+p(i,jy1,kp))   !*ros
c	 else
c	    write(82,8)1.e+03
c	 endif
 180	 continue
	 close(82)
      enddo
c
      do jjc=1,nsez14
	 js=jsez14(jjc)
	 jsy=jsym(js)
	 jp=jpv(js)
	 jpsy=jpv(jsy)
	 if(myrank.eq.0) then
	    jjcd=jjc/10
	    jjcu=mod(jjc,10)
	    open(82,file='results3_2D_'//char(48+iplant9c)
     %//char(48+iplant9d)//char(48+iplant9p)//
     %'_slice'//char(48+jjcd)//char(48+jjcu)//'_plotinst.q')
	    write(82,101)1,2*(nx-2)/nxp,(nzg-2)/nzp
	    i1=1
	    write(82,102)ntime,i1,i1,i1
!
!Density
!
	    do 240 k=kz1,kz2,nzp
	    do 241 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,9)0.5*(dens(i,jsy,k)+dens(i,jpsy,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
 241	    continue
	    do 242 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,9)0.5*(dens(i,js,k)+dens(i,jp,k))
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
     %(densr2(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,22,
     %MPI_COMM_EDDY,STATUS,IERR)
c		  CALL MPI_RECV
c    %(flagvor2sym(ix1:ix2,kk1:kk2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,jjp,11,
c    %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(pr2sym(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,12,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(vor2sym(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,13,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(uor2sym(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,14,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(wor2sym(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,15,
     %MPI_COMM_EDDY,STATUS,IERR)
		  CALL MPI_RECV
     %(densr2sym(ix1:ix2,kk1:kk2),(nx-2)*1*(nz-2),MTYPE,jjp,23,
     %MPI_COMM_EDDY,STATUS,IERR)
		  do 243 kr=kk1,kk2,nzp
		  do 244 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,9)densr2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 244		  continue
		  do 245 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)densr2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 245		  continue
 243		  continue
	       enddo
	    endif
!
!Azimuthal velocity
!
	    do 250 k=kz1,kz2,nzp
	    do 251 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,9)vo(i,jsy,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 251	    continue
	    do 252 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       write(82,9)vo(i,js,k)
c	    else
c	       write(82,9)1.e+03
c	    endif
 252	    continue
 250	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 253 kr=kk1,kk2,nzp
		  do 254 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,9)vor2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 254		  continue
		  do 255 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)vor2(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 255		  continue
 253		  continue
	       enddo
	    endif
!
!Radial velocity
!
	    do 260 k=kz1,kz2,nzp
	    do 261 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       im=i-1
	       write(82,9)
     %0.25*(uo(i,jsy,k)+uo(im,jsy,k)+uo(i,jpsy,k)+uo(im,jpsy,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
 261	    continue
	    do 262 i=ix1,ix2,nxp
c	    if(all(flagvo(i,js,k,:)).le.0) then
	       im=i-1
	       write(82,9)
     %0.25*(uo(i,js,k)+uo(im,js,k)+uo(i,jp,k)+uo(im,jp,k))
c	    else
c	       write(82,9)1.e+03
c	    endif
 262	    continue
 260	    continue
	    if(mysize.ne.1) then
	       do jjp=1,mysize-1
		  kk1=kz1+jjp*(nz-2)
		  kk2=kz2+jjp*(nz-2)
		  do 263 kr=kk1,kk2,nzp
		  do 264 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,9)uor2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 264		  continue
		  do 265 ir=ix1,ix2,nxp
c		  if(all(flagvor2(ir,kr,:)).le.0) then
		     write(82,9)uor2(ir,kr)
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
	    do 271 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,9)
     %0.25*(wo(i,jsy,k)+wo(i,jpsy,k)+wo(i,jsy,km)+wo(i,jpsy,km))
c	    else
c	       write(82,9)1.e+03
c	    endif
 271	    continue
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
		  do 274 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,9)wor2sym(ir,kr)
c		  else
c		     write(82,9)1.e+03
c		  endif
 274		  continue
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
	    do 281 i=ix2,ix1,-nxp
c	    if(all(flagvo(i,jsy,k,:)).le.0) then
	       write(82,8)0.5*(p(i,jsy,k)+p(i,jpsy,k))	  !*ros
c	    else
c	       write(82,8)1.e+03
c	    endif
 281	    continue
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
		  do 284 ir=ix2,ix1,-nxp
c		  if(all(flagvor2sym(ir,kr,:)).le.0) then
		     write(82,8)pr2sym(ir,kr)	 !*ros
c		  else
c		     write(82,8)1.e+03
c		  endif
 284		  continue
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
            CALL MPI_SEND
     %((dens(ix1:ix2,js,kz1:kz2)+dens(ix1:ix2,jp,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,22,
     %MPI_COMM_EDDY,STATUS,IERR)
c	    CALL MPI_SEND
c    %(flagvo(ix1:ix2,jsy,kz1:kz2,:),
c    %(nx-2)*1*(nz-2)*nbd,MPI_INTEGER,0,11,
c    %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((p(ix1:ix2,jsy,kz1:kz2)+p(ix1:ix2,jpsy,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,12,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %(vo(ix1:ix2,jsy,kz1:kz2),(nx-2)*1*(nz-2),MTYPE,0,13,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((uo(1:ix2-1,jsy,kz1:kz2)+uo(ix1:ix2,jsy,kz1:kz2)+
     %uo(1:ix2-1,jpsy,kz1:kz2)+uo(ix1:ix2,jpsy,kz1:kz2))*0.25,
     %(nx-2)*1*(nz-2),MTYPE,0,14,
     %MPI_COMM_EDDY,STATUS,IERR)
	    CALL MPI_SEND
     %((wo(ix1:ix2,jsy,1:kz2-1)+wo(ix1:ix2,jpsy,1:kz2-1)+
     %wo(ix1:ix2,jsy,kz1:kz2)+wo(ix1:ix2,jpsy,kz1:kz2))*0.25,
     %(nx-2)*1*(nz-2),MTYPE,0,15,
     %MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_SEND
     %((dens(ix1:ix2,jsy,kz1:kz2)+dens(ix1:ix2,jpsy,kz1:kz2))*0.5,
     %(nx-2)*1*(nz-2),MTYPE,0,23,
     %MPI_COMM_EDDY,STATUS,IERR)
	 endif
      enddo
      return
      end
c
