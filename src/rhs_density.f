C 
C-----SUBROUTINE-Rhs_density------------------A. Posa  03/12/2014------
C 
      subroutine rhs_density(uo,vo,wo,dens,kv,rr,rt,nx,ny,nz,yc,yv)
c
c      creates the right side of the density equation                        
c      (advection and viscous terms)                                  
c
c-----------------------------------------------------------------------
c
      include 'common.h'
c 
c-----------------------------------------------------------------------
c
      integer nx,ny,nz,err
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),dens(nx,ny,nz)
      real    rr(nx,ny,nz),rt(nx,ny,nz),kv(nx,ny,nz)
      real    qr(nx,ny,nz),yc(ny),yv(ny)

      integer i,j,k,var
      real    rflagx,rflagy,rflagcx,rflagcy
      real    uxp,uxm
      real    vyp,vym
      real    wzp,wzm
      real    rxp,rxm,ryp,rym,rzp,rzm
      real    drdxp,drdxm
      real    drdyp,drdym
      real    drdzp,drdzm
      real    h11,h22,h33,t11,t22,t33
      real    repr1,qxp,qxm
      real    txxp,txxm,tyyp,tyym,tzzp,tzzm
      real    kvip,kvim,kvjp,kvjm,kvkp,kvkm
c
		
! 	uo=0.0d0
! 	vo=0.0d0
! 	wo=0.0d0
      repr1 = ru1/prn
      rflagx=real(implx)
      rflagcx=real(implcx)
      rflagy=real(imply)
      rflagcy=real(implcy)

	
      IF(ICYL==0) THEN
        QR = UO
      ELSE
        DO I=1,NX
          QR(I,:,:) = RU(I)*UO(I,:,:)
        ENDDO
      ENDIF
	
      rr = 0.0
c
      do k=kz1,kz2
      do j=jy1,jy2
      do i=ix1,ix2
c
c...get density at staggered locations
c
        qxp = qr(i,j,k)
        qxm = qr(i-1,j,k)
        rxp = (dens(i+1,j,k)+dens(i,j,k))*0.5
        rxm = (dens(i,j,k)+dens(i-1,j,k))*0.5
c
        vyp = vo(i,j,k)
        vym = vo(i,j-1,k)
        ryp = (dens(i,j+1,k)+dens(i,j,k))*0.5
        rym = (dens(i,j,k)+dens(i,j-1,k))*0.5
c
        wzp = wo(i,j,k)
        wzm = wo(i,j,k-1)
        rzp = (dens(i,j,k+1)+dens(i,j,k))*0.5
        rzm = (dens(i,j,k)+dens(i,j,k-1))*0.5
c
c...get derivatives at 1/2 locations
c
        drdxp = au(i)*(dens(i+1,j,k)-dens(i,j,k))
        drdxm = au(i-1)*(dens(i,j,k)-dens(i-1,j,k))
c
        drdyp = bv(j)*(dens(i,j+1,k)-dens(i,j,k))
        drdym = bv(j-1)*(dens(i,j,k)-dens(i,j-1,k))
c
        drdzp = cw(k)*(dens(i,j,k+1)-dens(i,j,k))
        drdzm = cw(k-1)*(dens(i,j,k)-dens(i,j,k-1))

c...get ksgs at staggered locations

        kvip = 0.5*(kv(i+1,j,k)+kv(i,j,k))
        kvim = 0.5*(kv(i-1,j,k)+kv(i,j,k))

        kvjp = 0.5*(kv(i,j+1,k)+kv(i,j,k))
        kvjm = 0.5*(kv(i,j-1,k)+kv(i,j,k))

        kvkp = 0.5*(kv(i,j,k+1)+kv(i,j,k))
        kvkm = 0.5*(kv(i,j,k-1)+kv(i,j,k))

         if(myid==0) then
C            write(6,*) "kvip = ", kvip
C            write(6,*) "kvim = ", kvim
C            write(6,*) "kvjp = ", kvjp
C            write(6,*) "kvjm = ", kvjm
C            write(6,*) "kvkp = ", kvkp
C            write(6,*) "kvkm = ", kvkm

c            write(6,*) "rp(i) = ", rp(i)
c             write(6,*) "ru(i) = ", ru(i)
c            write(6,*) "repr1 = ", repr1



c            write(6,*) "rflagcx = ", rflagcx
c            write(6,*) "ruimplx(i) = ", ruimplx(i)

c            write(6,*) "rflagcy = ", rflagcy
c            write(6,*) "ruimply(i) = ", ruimply(i)


         endif
c
c
c
	txxp=(repr1+kvip)*drdxp*ru(i)
        txxm=(repr1+kvim)*drdxm*ru(i-1)
c
        tyyp=(repr1+kvjp)*drdyp
        tyym=(repr1+kvjm)*drdym
c
        tzzp=(repr1+kvkp)*drdzp
        tzzm=(repr1+kvkm)*drdzm

c...advective terms in conservative formulation
        h11 = -ap(i)*(rxp*qxp-rxm*qxm)/rp(i)
        h22 = -bp(j)*(ryp*vyp-rym*vym)/rp(i)
        h33 = -cp(k)*(rzp*wzp-rzm*wzm)
c
c...viscous terms
c 
        t11 = ap(i)*(txxp-txxm)/rp(i)
        t22 = bp(j)*(tyyp-tyym)/rp(i)**2
        t33 = cp(k)*(tzzp-tzzm)
c
c......calculate rhs for density equation
c
!         rr(i,j,k) = h11+h22+h33+t11+t22+t33

	rr(i,j,k) =(1.-rflagcx*ruimplx(i))*h11+(1.-rflagcy*ruimply(i))*h22+h33
     &            +(1.-rflagx *ruimplx(i))*t11+(1.-rflagy *ruimply(i))*t22+t33
! 	IF(icyl==1.and.idens==1) then
! ! 	rr(i,j,k)= rr(i,j,k)-denP1*(0.50d0*(uo(i,j,k)+uo(i-1,j,k))*sin(yc(j))+0.50d0*(vo(i,j,k)*cos(yv(j))+vo(i,j-1,k)*cos(yv(j-1))))
! 
! 	rr(i,j,k)=rr(i,j,k)-denP1*(0.50d0*(uo(i,j,k)+uo(i-1,j,k))*sin(yc(j))+0.50d0*(vo(i,j,k)+vo(i,j-1,k))*cos(yc(j)))
! 	ENDIF
! 	if(i==2.and.k==5) write(*,*), rr(i,j,k)
                
	rt(i,j,k) = rflagx*ruimplx(i)*t11 + rflagy*ruimply(i)*t22 
     &           + rflagcx*ruimplx(i)*h11 + rflagcy*ruimply(i)*h22

c$$$        if(myrank==0) then 
c$$$           write(*,*) "rflagx*ruimplx(",i,") = ", rflagx*ruimplx(i)
c$$$           write(*,*) "rflagy*ruimply(",i,") = ", rflagy*ruimply(i)
c$$$        endif

      enddo
      enddo
      enddo
!  	WRITE(*,*),maxval(rr(2,:,:)), minval(rr(2,:,:))
! 	WRITE(*,*),dens(1,jsym(1),17),dens(2,1,17)
	var=5 !u=Flow(:,:,:,5)
!  	do dir=1,3
!  	do face=1,2
!    	call Sponge(rr,dens,var,nx,ny,nz,err)

      return
      end



