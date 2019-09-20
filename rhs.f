C 
C-----SUBROUTINE-Rhs---------------------------E. Balaras  4/12/98------
C 
      subroutine rhs(uo,vo,wo,dens,ut,vt,wt,tv,ur,vr,wr,qr,nx,ny,nz,xc,xu,yc,yv,time)
c
c      creates the right side of the mom. eq.                         
c      (advection and viscous terms)                                  
c
c uo,vo,wo,tv: input (velocities and eddy viscosity)
c ut,vt,wt: output (associated with the implicit discretization of the
c convective and viscous terms)
c ur,vr,wr: RHS (without the terms treated by C.N.)
c the viscous and convective terms of derivative along the X and Y directions
c are eventually treated implicitly
c
c-----------------------------------------------------------------------
c
      use density_bg
      include 'common.h'
      
c-----------------------------------------------------------------------
c
      integer nx,ny,nz
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),tv(nx,ny,nz),dens(nx,ny,nz)
      real    ur(nx,ny,nz),vr(nx,ny,nz),wr(nx,ny,nz)
      real    ut(nx,ny,nz),vt(nx,ny,nz),wt(nx,ny,nz)
      real    qr(nx,ny,nz),time,xc(nx),xu(nx),yc(ny),yv(ny)

      integer i,j,k,ok,face,err,var,index1,index2
      real    rflagx,rflagy,rflagcx,rflagcy
      real    uxp,uxm,vxp,vxm,wxp,wxm
      real    uyp,uym,vyp,vym,wyp,wym
      real    uzp,uzm,vzp,vzm,wzp,wzm
      real    dvdxp,dvdxm,dwdxp,dwdxm
      real    dudyp,dudym,dvdyp,dvdym,dwdyp,dwdym
      real    dudzp,dudzm,dvdzp,dvdzm,dwdzp,dwdzm
      real    tvip,tvim,tvjp,tvjm,tvkp,tvkm
      real    txxp,txxm,tyyp,tyym,tzzp,tzzm
      real    txyp,txym,txzp,txzm,tyzp,tyzm
      real    qxp,qxm,qyp,qym,qzp,qzm
      real    dqdxp,dqdxm
      real    dtvdx,dtvdy
      real    dudy,dvdy,vxc,uyc
      real    h11,h22,h33,t11,t22,t33,s11,s22,s33
      real    cyl

      integer nzg,kg
c
      rflagx=real(implx)
      rflagcx=real(implcx)
      rflagy=real(imply)
      rflagcy=real(implcy)
c
c      write(6,*) 'rflagx=',rflagx,', rflagy=',rflagy,', rflagcx=',rflagcx,', rflagcy=',rflagcy
c
     
      IF(ICYL==0) THEN
        QR = UO
      ELSE
        DO I=1,NX
          QR(I,:,:) = RU(I)*UO(I,:,:)
        ENDDO
      ENDIF        

      ut = 0.0
      ur = 0.0
      vt = 0.0
      vr = 0.0
      wt = 0.0
      wr = 0.0
	index1=65
c
cc
ccccccccccccccccccc     for the u component     cccccccccccccccccc
cc
c
      do k=kz1,kz2
      do j=jy1,jy2
      do i=ibu,ieu
c...get velocities at 1/2 locations
c
        qxp = (qr(i+1,j,k)+qr(i  ,j,k))*0.5
        qxm = (qr(i  ,j,k)+qr(i-1,j,k))*0.5
c
        uxp = (uo(i+1,j,k)+uo(i  ,j,k))*0.5
        uxm = (uo(i  ,j,k)+uo(i-1,j,k))*0.5
c
        uyp = (uo(i,j+1,k)+uo(i,j  ,k))*0.5
        uym = (uo(i,j  ,k)+uo(i,j-1,k))*0.5
c
        uzp = (uo(i,j,k+1)+uo(i,j,k  ))*0.5
        uzm = (uo(i,j,k  )+uo(i,j,k-1))*0.5
c
        vxp = (vo(i+1,j  ,k)+vo(i,j  ,k))*0.5
        vxm = (vo(i+1,j-1,k)+vo(i,j-1,k))*0.5
c
        wxp = (wo(i+1,j,k  )+wo(i,j,k  ))*0.5
        wxm = (wo(i+1,j,k-1)+wo(i,j,k-1))*0.5
c
c...get derivatives at 1/2 locations
        dqdxp = ap(i+1)*(qr(i+1,j,k)-qr(i  ,j,k))
        dqdxm = ap(i  )*(qr(i  ,j,k)-qr(i-1,j,k))
c
        dudyp = bv(j  )*(uo(i,j+1,k)-uo(i,j  ,k))
        dudym = bv(j-1)*(uo(i,j  ,k)-uo(i,j-1,k))
c
        dudzp = cw(k  )*(uo(i,j,k+1)-uo(i,j,k  ))
        dudzm = cw(k-1)*(uo(i,j,k  )-uo(i,j,k-1))
c
        dvdxp = au(i)*(vo(i+1,j  ,k)-vo(i,j  ,k))
        dvdxm = au(i)*(vo(i+1,j-1,k)-vo(i,j-1,k))
c
        dwdxp = au(i)*(wo(i+1,j,k  )-wo(i,j,k  ))
        dwdxm = au(i)*(wo(i+1,j,k-1)-wo(i,j,k-1))
c
c...get nu_t where needed
c
        tvjp = 0.25*(tv(i+1,j,k)+tv(i,j,k)+tv(i+1,j+1,k)+tv(i,j+1,k))
        tvjm = 0.25*(tv(i+1,j,k)+tv(i,j,k)+tv(i+1,j-1,k)+tv(i,j-1,k))
c
        tvkp = 0.25*(tv(i+1,j,k)+tv(i,j,k)+tv(i+1,j,k+1)+tv(i,j,k+1))
        tvkm = 0.25*(tv(i+1,j,k)+tv(i,j,k)+tv(i+1,j,k-1)+tv(i,j,k-1))
c
c...flux of normal total stresses
c
        txxp=(ru1+2.*tv(i+1,j,k))*dqdxp/rp(i+1)
        txxm=(ru1+2.*tv(i  ,j,k))*dqdxm/rp(i  )
c
        tyyp=(ru1+tvjp)*dudyp
        tyym=(ru1+tvjm)*dudym
c
        tzzp=(ru1+tvkp)*dudzp
        tzzm=(ru1+tvkm)*dudzm
c...flux of cross sgs stresses
        txyp=tvjp*dvdxp
        txym=tvjm*dvdxm
c
        txzp=tvkp*dwdxp
        txzm=tvkm*dwdxm
c..advective term in conservative formulation
        h11 = -au(i)*(qxp*uxp-qxm*uxm)/ru(i)
        h22 = -bu(j)*(vxp*uyp-vxm*uym)/ru(i)
        h33 = -cu(k)*(wxp*uzp-wxm*uzm)
c..viscous+part of sgs diffusion
        t11 = au(i)*(txxp-txxm)
        t22 = bu(j)*(tyyp-tyym)/ru(i)**2
        t33 = cu(k)*(tzzp-tzzm)
c..rest of sgs diffusion
        s22 = bu(j)*(txyp-txym)/ru(i)
        s33 = cu(k)*(txzp-txzm)
c
        if(icyl==1) then
c
          tvip = 0.5*(tv(i+1,j,k)+tv(i,j,k))
c
          vxc = 0.5*(vxp+vxm)
c
          dvdy = bu(j)*(vxp-vxm)
c
          dtvdx = au(i)*(tv(i+1,j,k)-tv(i,j,k))
          dtvdy = bu(j)*(tvjp       -tvjm     )
c
          cyl = (vxc**2
c     
     &         -(3.*tvip+2.*ru1)*dvdy/ru(i)
     &         -vxc*dtvdy/ru(i)
     &         -2*uo(i,j,k)*dtvdx
     &         )/ru(i)
        else
          cyl = 0.
        endif
c
c......calculate rhs for u-momentum
c
        ur(i,j,k) = (1.-rflagcx*ruimplx(i))*h11+(1.-rflagcy*ruimply(i))*h22+h33
     &       +(1.-rflagx*ruimplx(i))*t11+(1.-rflagy*ruimply(i))*t22+t33
c     &       +t11+t22+t33
     &       +s22+s33
     &       +cyl
! 	IF(idens.EQ.1) ur(i,j,k)=ur(i,j,k)-0.50d0*grav*(dens(i+1,j,k)+dens(i,j,k))
! 	IF(icyl.eq.1.and.idens==1)  ur(i,j,k)=ur(i,j,k)-0.50d0*grav*sin(yc(j))*((dens(i+1,j,k)-1.0d0-denP1*rp(i+1)*sin(yc(j)))+
!      &		(dens(i,j,k)-1.0d0-denP1*rp(i)*sin(yc(j))))
! 	if(index1.gt.sz.and.index1.lt.ez.or.index1.eq.sz.or.index1.eq.ez)index2=index1-(sz-2)
! 	if(i.eq.2.and.j.eq.2.and.k.eq.index2) ur(i,j,index2)=ur(i,j,index2)+0.1d0*cos(0.8d0*time)*sin(yc(j))
c
c..store explicit part for C.N.
! 	if(i==2.and.k==5) write(*,*), ur(I,J,K)
! 	ur(i,j,k)=ur(i,j,k)-0.50d0*grav*sin(yc(j))*((dens(i+1,j,k)-dens_bg(i+1,j))+
!      &		(dens(i,j,k)-dens_bg(i,j)))
! 	if(i==2.and.k==5) write(*,*)"*************",ur(I,J,K)
        ut(i,j,k) = rflagx*ruimplx(i)*t11 + rflagy*ruimply(i)*t22 + rflagcx*ruimplx(i)*h11 + rflagcy*ruimply(i)*h22

	
c
!        nzg=130
!        kg=k+myrank*(nz-2)
!        if((i.eq.nx/2+1).and.(j.eq.ny/2).and.(kg.eq.nzg/2))
!     %write(60,*)i,j,k,myrank,ur(i,j,k)    
!        if((i.eq.nx/2).and.(j.eq.ny/2).and.(kg.eq.nzg/2))
!     %write(60,*)i,j,k,myrank,ur(i,j,k)
!        if((i.eq.nx/2-1).and.(j.eq.ny/2).and.(kg.eq.nzg/2))
!     %write(60,*)i,j,k,myrank,ur(i,j,k)
      enddo
      enddo
      enddo

	var=1 !u=Flow(:,:,:,5)
!  	do dir=1,3
!  	do face=1,2
!    	 call Sponge(ur,uo,var,nx,ny,nz,err)
! 	 call Sponge(ut,uo,var,nx,ny,nz,err)
!  	enddo
!  	enddo
c
cc
ccccccccccccccccccc     for the v component     cccccccccccccccccc
cc
c
      do k=kz1,kz2
      do j=jbv,jev
      do i=ix1,ix2
c...get velocities at 1/2 locations
        vxp=(vo(i+1,j,k)+vo(i  ,j,k))*0.5
        vxm=(vo(i  ,j,k)+vo(i-1,j,k))*0.5
c     
        vyp=(vo(i,j+1,k)+vo(i,j  ,k))*0.5
        vym=(vo(i,j  ,k)+vo(i,j-1,k))*0.5
c     
        vzp=(vo(i,j,k+1)+vo(i,j,k  ))*0.5
        vzm=(vo(i,j,k  )+vo(i,j,k-1))*0.5
c     
        qyp=(qr(i  ,j+1,k)+qr(i  ,j,k))*0.5
        qym=(qr(i-1,j+1,k)+qr(i-1,j,k))*0.5
c     
        wyp=(wo(i,j+1,k  )+wo(i,j,k  ))*0.5   
        wym=(wo(i,j+1,k-1)+wo(i,j,k-1))*0.5

c...get derivatives at 1/2 locations
        dvdxp= au(i  )*(vo(i+1,j,k)-vo(i  ,j,k))
        dvdxm= au(i-1)*(vo(i  ,j,k)-vo(i-1,j,k))
c
        dvdyp= bp(j+1)*(vo(i,j+1,k)-vo(i,j  ,k))
        dvdym= bp(j  )*(vo(i,j  ,k)-vo(i,j-1,k))
c
        dvdzp= cw(k  )*(vo(i,j,k+1)-vo(i,j,  k))
        dvdzm= cw(k-1)*(vo(i,j,k  )-vo(i,j,k-1))
c
        dudyp= bv(j)*(uo(i  ,j+1,k)-uo(i  ,j,k))
        dudym= bv(j)*(uo(i-1,j+1,k)-uo(i-1,j,k))
c
        dwdyp= bv(j)*(wo(i,j+1,k  )-wo(i,j,k  ))
        dwdym= bv(j)*(wo(i,j+1,k-1)-wo(i,j,k-1))
c
c...get nu_t where needed
c
        tvip = 0.25*(tv(i,j+1,k)+tv(i,j,k)+tv(i+1,j+1,k)+tv(i+1,j,k))
        tvim = 0.25*(tv(i,j+1,k)+tv(i,j,k)+tv(i-1,j+1,k)+tv(i-1,j,k))
c
        tvkp = 0.25*(tv(i,j+1,k)+tv(i,j,k)+tv(i,j+1,k+1)+tv(i,j,k+1))
        tvkm = 0.25*(tv(i,j+1,k)+tv(i,j,k)+tv(i,j+1,k-1)+tv(i,j,k-1))
c
c...flux of normal total stresses
        txxp=(ru1+tvip)*dvdxp*ru(i  )
        txxm=(ru1+tvim)*dvdxm*ru(i-1)
c
        tyyp=(ru1+2.*tv(i,j+1,k))*dvdyp
        tyym=(ru1+2.*tv(i,j  ,k))*dvdym
c
        tzzp=(ru1+tvkp)*dvdzp
        tzzm=(ru1+tvkm)*dvdzm
c...flux of cross sgs stresses
        txyp=tvip*dudyp
        txym=tvim*dudym
c
        tyzp=tvkp*dwdyp
        tyzm=tvkm*dwdym
c..advective term in conservative formulation
        h11 = -av(i)*(qyp*vxp-qym*vxm)/rp(i)
        h22 = -bv(j)*(vyp*vyp-vym*vym)/rp(i)
        h33 = -cv(k)*(wyp*vzp-wym*vzm)
c
c..viscous+part of sgs diffusion
        t11 = av(i)*(txxp-txxm)/rp(i)
        t22 = bv(j)*(tyyp-tyym)/rp(i)**2
        t33 = cv(k)*(tzzp-tzzm)
c..rest of sgs diffusion
        s11 = av(i)*(txyp-txym)/rp(i)
        s33 = cv(k)*(tyzp-tyzm)/rp(i)
c
        if(icyl==1) then
c
          tvjp = 0.5*(tv(i,j+1,k)+tv(i,j,k)) 
c
          uyp = (uo(i,j+1,k)+uo(i-1,j+1,k))*0.5
          uym = (uo(i,j  ,k)+uo(i-1,j  ,k))*0.5
c
          uyc = 0.5*(uyp+uym)
c
          dudy = bv(j)*(uyp-uym)
c
          dtvdx = av(i)*(tvip       -tvim     )
          dtvdy = bv(j)*(tv(i,j+1,k)-tv(i,j,k))

          cyl = (-vo(i,j,k)*uyc
c
     &         -dtvdx*vo(i,j,k)
     &         -(tvjp+ru1)*vo(i,j,k)/rp(i)
     &         +(3.*tvjp+2.*ru1)*dudy/rp(i)
     &         +2.*uyc*dtvdy/rp(i)
     &         )/rp(i)
        else
          cyl = 0.
        endif
c
c......calculate rhs for v-momentum
c
        vr(i,j,k) = (1.-rflagcx*rvimplx(i))*h11+(1.-rflagcy*rvimply(i))*h22+h33
     &       +(1.-rflagx*rvimplx(i))*t11+(1.-rflagy*rvimply(i))*t22+t33
c     &       +t11+t22+t33
     &       +s11+s33
     &       +cyl

! 	if(i.eq.2.and.j.eq.2.and.k.eq.index2) vr(i,j,index2)=vr(i,j,index2)+0.1d0*cos(0.8d0*time)*cos(yc(j))

! 	IF(icyl.eq.1.and.idens==1)  vr(i,j,k)=vr(i,j,k)-0.50d0*grav*cos(yv(j))*((dens(i,j+1,k)-1.0d0-denP1*rp(i)*sin(yc(j+1)))
!      &                              +(dens(i,j,k)-1.0d0-denP1*rp(i)*sin(yc(j))))

! 	vr(i,j,k)=vr(i,j,k)-0.50d0*grav*cos(yv(j))*((dens(i,j+1,k)-dens_bg(i,j+1))
!      &                              +(dens(i,j,k)-dens_bg(i,j)))
c..store explicit part for C.N.
        vt(i,j,k) = rflagx*rvimplx(i)*t11 + rflagy*rvimply(i)*t22 + rflagcx*rvimplx(i)*h11 + rflagcy*rvimply(i)*h22
c
      enddo
      enddo
      enddo
     	var=2 !v=Flow(:,:,:,2)
!  	do dir=1,3
!  	do face=1,2
!    	 call Sponge(vr,vo,var,nx,ny,nz,err)
! 	 call Sponge(vt,vo,var,nx,ny,nz,err)
!  	enddo
!  	enddo
cc
ccccccccccccccccccc     for the w component     cccccccccccccccccc
cc
c
      do k=kbw,kew
      do j=jy1,jy2
      do i=ix1,ix2
c...get velcities at 1/2 locations
        wxp=(wo(i+1,j,k)+wo(i  ,j,k))*0.5
        wxm=(wo(i  ,j,k)+wo(i-1,j,k))*0.5
c     
        wyp=(wo(i,j+1,k)+wo(i,j  ,k))*0.5   
        wym=(wo(i,j  ,k)+wo(i,j-1,k))*0.5
c     
        wzp=(wo(i,j,k+1)+wo(i,j,k  ))*0.5
        wzm=(wo(i,j,k  )+wo(i,j,k-1))*0.5
c     
        qzp=(qr(i  ,j,k+1)+qr(i  ,j,k))*0.5
        qzm=(qr(i-1,j,k+1)+qr(i-1,j,k))*0.5
c     
        vzp=(vo(i,j  ,k+1)+vo(i,j  ,k))*0.5
        vzm=(vo(i,j-1,k+1)+vo(i,j-1,k))*0.5
c
c...get derivatives at 1/2 locations
        dwdxp= au(i  )*(wo(i+1,j,k)-wo(i  ,j,k))
        dwdxm= au(i-1)*(wo(i  ,j,k)-wo(i-1,j,k))
c
        dwdyp= bv(j  )*(wo(i,j+1,k)-wo(i,j  ,k))
        dwdym= bv(j-1)*(wo(i,j  ,k)-wo(i,j-1,k))
c
        dwdzp= cp(k+1)*(wo(i,j,k+1)-wo(i,j,k  ))
        dwdzm= cp(k  )*(wo(i,j,k  )-wo(i,j,k-1))
c
        dudzp= cw(k)*(uo(i  ,j,k+1)-uo(i  ,j,k))
        dudzm= cw(k)*(uo(i-1,j,k+1)-uo(i-1,j,k))
c
        dvdzp= cw(k)*(vo(i,j  ,k+1)-vo(i,j  ,k))
        dvdzm= cw(k)*(vo(i,j-1,k+1)-vo(i,j-1,k))
c
c...get nu_t where needed
c
        tvip = 0.25*(tv(i,j,k+1)+tv(i,j,k)+tv(i+1,j,k+1)+tv(i+1,j,k))
        tvim = 0.25*(tv(i,j,k+1)+tv(i,j,k)+tv(i-1,j,k+1)+tv(i-1,j,k))
c
        tvjp = 0.25*(tv(i,j,k+1)+tv(i,j,k)+tv(i,j+1,k+1)+tv(i,j+1,k))
        tvjm = 0.25*(tv(i,j,k+1)+tv(i,j,k)+tv(i,j-1,k+1)+tv(i,j-1,k))
c
c...flux of normal total stresses
        txxp=(ru1+tvip)*dwdxp*ru(i  )
        txxm=(ru1+tvim)*dwdxm*ru(i-1)
c
        tyyp=(ru1+tvjp)*dwdyp
        tyym=(ru1+tvjm)*dwdym
c
        tzzp=(ru1+2.*tv(i,j,k+1))*dwdzp
        tzzm=(ru1+2.*tv(i,j,k  ))*dwdzm
c...flux of cross sgs stresses
        txzp=tvip*dudzp*ru(i  )
        txzm=tvim*dudzm*ru(i-1)
c
        tyzp=tvjp*dvdzp
        tyzm=tvjm*dvdzm
c..advective term in conservative formulation
        h11 = -aw(i)*(qzp*wxp-qzm*wxm)/rp(i)
        h22 = -bw(j)*(vzp*wyp-vzm*wym)/rp(i)
        h33 = -cw(k)*(wzp*wzp-wzm*wzm)
c..viscous+part of sgs diffusion
        t11 = aw(i)*(txxp-txxm)/rp(i)
        t22 = bw(j)*(tyyp-tyym)/rp(i)**2
        t33 = cw(k)*(tzzp-tzzm)
c..rest of sgs diffusion
        s11 = aw(i)*(txzp-txzm)/rp(i)
        s22 = bw(j)*(tyzp-tyzm)/rp(i)
c
c......calculate rhs for w-momentum
c
        wr(i,j,k) = (1.-rflagcx*rvimplx(i))*h11+(1.-rflagcy*rvimply(i))*h22+h33
     &       +(1.-rflagx*rvimplx(i))*t11+(1.-rflagy*rvimply(i))*t22+t33
c     &       +t11+t22+t33
     &       +s11+s22
c..store explicit part of wall normal diffusio for c.n.
        wt(i,j,k) = rflagx*rvimplx(i)*t11 + rflagy*rvimply(i)*t22 + rflagcx*rvimplx(i)*h11 + rflagcy*rvimply(i)*h22
c
      enddo
      enddo
      enddo
c
	var=3 !w=Flow(:,:,:,3)
!  	do dir=1,3
!  	do face=1,2
!    	 call Sponge(wr,wo,var,nx,ny,nz,err)
! 	 call Sponge(wt,wo,var,nx,ny,nz,err)
!  	enddo
!  	enddo
      return
      end


	subroutine gravity(time1)
!@h
!   Comments:
!     Ramp function should bring g smoothly from 0 to g_z_nd_orig in time
!     gt2-gt1. Using tanh this should be
!     g(t) = g*tanh( 1.5 * (time-gt1)/(gt2-gt1) )
!@q
	include 'common.h'
 
 !Passed Variables
 	real,intent(in)  :: time1
! 	real,intent(out)    :: g1
!	real		    :: gt1,gt2,g_orig
 	

!	 gt1=0.0d0
!	 gt2=1.0d0
!	 g_orig=9.81d0

!	 g_orig=18.0d0

c Re=10000 & Fr = 3. 
c use unstratified restart file starting at around 126
c restart file *_000314000.res
	 !gt1=0.0d0
	 !gt2=50.0d0

	 !g_orig = 0.0d0


	 if (time1.LT.gt1) then
   	 grav = 0.d0
 	 elseif (time1.GT.gt1.AND.time1.LT.gt2) then
! tanh ramping
   	 grav = g_orig * dtanh( 1.5d0*(time1-gt1)/(gt2-gt1) )
! linear ramping
!   g = g_orig * (time-gt1)/(gt2-gt1) 
 	 elseif (time1.GE.gt2) then
   	 grav = g_orig
 	 endif
         if(MYRANK==0) then
c            if (dabs(grav).LT.dabs(g_orig)) then 
               write(*,'(a,2(2x,f12.6))') "   GRAVITY: ",grav
c            endif
         endif
 	 return
	 end subroutine gravity
