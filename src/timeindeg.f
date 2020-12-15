* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*                                                                 *
*           Predictor step  :   USTARAB                           *
*           Divergence      :   DIVAB                             *
*           Velocity corect.:   CORAB                             *
*                                                                 *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


      SUBROUTINE PREDICTOR(UO,VO,WO,P,UA,VA,WA,UB,VB,WB,
     &     US,VS,WS,ALFXDT,GAMXDT,RHOXDT,IM,JM,KM,CLOCK,NCLOCK)
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
      INTEGER IM,JM,KM
      REAL UO(IM,JM,KM),VO(IM,JM,KM),WO(IM,JM,KM),P(IM,JM,KM),
     &     UA(IM,JM,KM),VA(IM,JM,KM),WA(IM,JM,KM),
     &     UB(IM,JM,KM),VB(IM,JM,KM),WB(IM,JM,KM),
     &     US(IM,JM,KM),VS(IM,JM,KM),WS(IM,JM,KM)
c
      INTEGER I,J,K
      REAL ALFXDT,GAMXDT,RHOXDT
      REAL USTAR,VSTAR,WSTAR
      REAL COEF
c
      integer nclock
      real clock(nclock)
      real clocktemp,clocktemp1
      real tclock
*
**** Compute provisional values for each momentum eq.
*
      COEF=0.5*ALFXDT
c------------------------------------------------------------------
c                                                    compute ustar
c------------------------------------------------------------------
      clocktemp1 = tclock()
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IBU,IEU
        USTAR = UO(I,J,K)
     &       +GAMXDT*UB(I,J,K)
     &       +RHOXDT*UA(I,J,K)
     &       -ALFXDT*AU(I)*(P(I+1,J,K)-P(I,J,K))
     &       +COEF*US(I,J,K)
     &       -ALFXDT*DPDX
c UO is the velocity at the time level L-1
c UB is the RHS (explicit part) at the time level L-1
c UA is the RHS (explicit part) at the time level L-2
c P is the pressure (for the pressure gradient) at the time level L-1
c US is the semi-explicit part for the terms treated by C.N.
C
        UA(I,J,K) = UB(I,J,K)
        UB(I,J,K) = US(I,J,K)
        US(I,J,K) = USTAR
c the RHS at the time level L-1 is stored in UA
c the C.N. terms are stored in UB
c the provisional solution is stored in US
      ENDDO
      ENDDO
      ENDDO
      clocktemp1 = tclock() - clocktemp1
      clock(1) = clock(1) + clocktemp1

c------------------------------------------------------------------
c                                                    compute vstar
c------------------------------------------------------------------
      clocktemp1 = tclock()
      DO K=KZ1,KZ2
      DO J=JBV,JEV
      DO I=IX1,IX2
        VSTAR = VO(I,J,K)
     &       +GAMXDT*VB(I,J,K)
     &       +RHOXDT*VA(I,J,K)
     &       -ALFXDT*BV(J)*(P(I,J+1,K)-P(I,J,K))/RP(I)
     &       +COEF*VS(I,J,K)
     &       -ALFXDT*DPDY
C
        VA(I,J,K) = VB(I,J,K)
        VB(I,J,K) = VS(I,J,K)
        VS(I,J,K) = VSTAR
      ENDDO
      ENDDO
      ENDDO
      clock(4) = clock(4) + tclock()-clocktemp1
c      write(6,*) myrank,'Vloop:',kz1,kz2,jbv,jev,ix1,ix2,
c     &     (kz2-kz1+1)*(jev-jbv+1)*(ix2-ix1+1)
c------------------------------------------------------------------
c                                                    compute wstar
c------------------------------------------------------------------
      clocktemp1 = tclock()
      DO K=KBW,KEW
      DO J=JY1,JY2
      DO I=IX1,IX2
        WSTAR = WO(I,J,K)
     &       +GAMXDT*WB(I,J,K)
     &       +RHOXDT*WA(I,J,K)
     &       -ALFXDT*CW(K)*(P(I,J,K+1)-P(I,J,K))
     &       +COEF*WS(I,J,K) 
     &       -ALFXDT*DPDZ
        
C
        WA(I,J,K) = WB(I,J,K)
        WB(I,J,K) = WS(I,J,K)
        WS(I,J,K) = WSTAR
      ENDDO
      ENDDO
      ENDDO
      clock(7) = clock(7) + tclock()-clocktemp1
c      write(6,*) myrank,'Wloop:',kbw,kew,jy1,jy2,ix1,ix2,
c     &     (kew-kbw+1)*(jy2-jy1+1)*(ix2-ix1+1)
c
      RETURN
      END    


      SUBROUTINE INVERSEY(US,VS,WS,TV,VO,ALFXDT,GAMXDT,RHOXDT,IM,JM,KM)
c
      INCLUDE 'common.h'
c
      INTEGER IM,JM,KM
      REAL US(IM,JM,KM),VS(IM,JM,KM),WS(IM,JM,KM),TV(IM,JM,KM)
      REAL VO(IM,JM,KM)
c
      INTEGER I,J,K,JJ,I1,I2,NI,ID
      REAL AY(MMY),BY(MMY),CY(MMY),RY(MMY),V1D(MMY)
C
      REAL ALFXDT,GAMXDT,RHOXDT
      REAL TVJP,TVJM,RFLAGY,RFLAGCY
      REAL COEF,COEFY,AYV,BYV,CYV,AYC,BYC,CYC
      REAL ALPHA,BETA

      REAL vxp,vxm,vzp,vzm,vyp,vym
*
**** Compute provisional values for each momentum eq.
*
      COEF=0.5*ALFXDT
c
      rflagy  = real(imply)  ! 1 if the viscous terms are treated implicitly
      rflagcy = real(implcy)   ! 1 if the convective terms are treated implicitly
c
c------------------------------------------------------------------
c                             If fully implicit compute vstar first!
c------------------------------------------------------------------
c
      do k=kz1,kz2
      do id=1,nimply
        i1 = implyvlim(id,1)  ! left border of the implicit region
        i2 = implyvlim(id,2)  ! right border of the implicit region
      do i=i1,i2
      do j=jbv,jev
c..fill arrays for 3-diagonal system
        jj=j-1
c...calculate total viscosity where needed
        tvjp=2.*tv(i,j+1,k)+ru1
        tvjm=2.*tv(i,j  ,k)+ru1
c...rhs
        ry(jj)=vs(i,j,k)
c
c...Viscous diffusion
c...coefy
        coefy  = rflagy*coef*bv(j)/rp(i)**2
c...v(j-1)
        ayv =-coefy*bp(j  )*tvjm
c...v(j+1)
        cyv =-coefy*bp(j+1)*tvjp
c...v(j)
c        by(jj) = by(jj)+coefy*(bp(j)*tvjm + bp(j+1)*tvjp)
        byv = -ayv - cyv
c
c...Convective terms
c...coefy
        coefy = rflagcy*coef*bv(j)/rp(i)
c...v(j-1)
        ayc =-0.25*coefy*(vo(i,j-1,k)+vo(i,j,k))*2.0
c...v(j+1)
        cyc = 0.25*coefy*(vo(i,j+1,k)+vo(i,j,k))*2.0
c...v(j)
c        by(jj) =by(jj)+0.25*coefy*(vo(i,j+1,k)-vo(i,j-1,k))*2.0
        byc = ayc+cyc

        ay(jj) = ayv+ayc
        cy(jj) = cyv+cyc
        by(jj) = 1.0 + byv + byc 

        vyp=(vo(i,j+1,k)+vo(i,j  ,k))*0.5
        vym=(vo(i,j  ,k)+vo(i,j-1,k))*0.5

c...the RHS of the system is modified in the case of the implicit
c...treatment of the convective terms along Y
        ry(jj) = ry(jj) + coefy*(vyp*vyp-vym*vym)

      enddo
c...periodic boundary condition
      beta  =ay(1)
      alpha =cy(jev-1)
c...solve the system
      call cyclic(ay,by,cy,alpha,beta,ry,v1d,jev-1)
c
c...the provisional velocity is stored in VS
      do j=jbv,jev
        jj=j-1
        vs(i,j,k)=v1d(jj)
      enddo
c
      enddo
      enddo
      enddo

c...REFRESHBC is called to allow the solution of the equations
c...for U and W using the provisional values of V evaluated above
      call refreshbc(vs,im*jm,km)
c------------------------------------------------------------------
c                                                    compute wstar
c------------------------------------------------------------------
      do k=kbw,kew
      do id=1,nimply
        i1 = implyvlim(id,1)
        i2 = implyvlim(id,2)
      do i=i1,i2
      do j=jy1,jy2
c..fill arrays for 3-diagonl system
        jj=j-1
c...calculate total viscosity where needed
        tvjp = 0.25*(tv(i,j,k+1)+tv(i,j,k)+tv(i,j+1,k+1)+tv(i,j+1,k))+ru1
        tvjm = 0.25*(tv(i,j,k+1)+tv(i,j,k)+tv(i,j-1,k+1)+tv(i,j-1,k))+ru1
c...rhs
        ry(jj)=ws(i,j,k)
c
c...Viscous diffusion
c...coefy
        coefy  = rflagy*coef*bw(j)/rp(i)**2
c...w(j-1)
        ayv =-coefy*bv(j-1)*tvjm
c...w(j+1)
        cyv =-coefy*bv(j  )*tvjp
c...w(j)
        byv =-ayv - cyv
c     
c...Convective term (the component V is known)
        vzp=(vs(i,j  ,k+1)+vs(i,j  ,k))*0.5
        vzm=(vs(i,j-1,k+1)+vs(i,j-1,k))*0.5
c...coefy
        coefy = rflagcy*0.5*coef*bw(j)/rp(i)
c...w(j-1)
        ayc =-coefy*vzm
c...w(j+1)
        cyc = coefy*vzp
c...w(j)
        byc = ayc + cyc

        ay(jj) = ayv + ayc
        cy(jj) = cyv + cyc
        by(jj) = 1.0  + byv + byc
      enddo
c...periodic boundary conditions
      beta =ay(1)
      alpha=cy(jy2-1)
c...solve the system
      call cyclic(ay,by,cy,alpha,beta,ry,v1d,jy2-1)
c
c...the provisional velocity is stored in WS
      do j=jy1,jy2
        jj=j-1
        ws(i,j,k)=v1d(jj)
      enddo

      enddo
      enddo
      enddo
c
c------------------------------------------------------------------
c                                                    compute ustar
c------------------------------------------------------------------
c
      do k=kz1,kz2
      do id=1,nimply
        i1 = implyulim(id,1)
        i2 = implyulim(id,2)
      do i=i1,i2
      do j=jy1,jy2
c..fill arrays for 3-diagonal system
        jj=j-1
c...calculate total viscosity where needed
        tvjp = 0.25*(tv(i+1,j,k)+tv(i,j,k)+tv(i+1,j+1,k)+tv(i,j+1,k))+ru1
        tvjm = 0.25*(tv(i+1,j,k)+tv(i,j,k)+tv(i+1,j-1,k)+tv(i,j-1,k))+ru1
c...rhs
        ry(jj)=us(i,j,k)
c
c...Viscous diffusion
c...coefy
        coefy = rflagy*coef*bu(j)/ru(i)**2
c...u(j-1)
        ayv =-coefy*bv(j-1)*tvjm
c...u(j+1)
        cyv =-coefy*bv(j  )*tvjp
c...u(j)
        byv =-ayv-cyv
c
c...Convective terms (the component V is known)
        vxp = (vs(i+1,j  ,k)+vs(i,j  ,k))*0.5
        vxm = (vs(i+1,j-1,k)+vs(i,j-1,k))*0.5
c...coefy
        coefy = rflagcy*0.5*coef*bu(j)/ru(i)
c...u(j-1)
        ayc =-coefy*vxm
c...u(j+1)
        cyc = coefy*vxp
c...u(j)
        byc = ayc+cyc

        ay(jj) = ayv + ayc
        cy(jj) = cyv + cyc
        by(jj) = 1.0 + byv + byc

      enddo
c...periodic boundary condition
      beta =ay(1)
      alpha=cy(jy2-1)
c...solve the system
      call cyclic(ay,by,cy,alpha,beta,ry,v1d,jy2-1)
c
c...the provisional velocity is stored in US
      do j=jy1,jy2
        jj=j-1
        us(i,j,k)=v1d(jj)
      enddo
c
      enddo
      enddo
      enddo

      call refreshbc(us,im*jm,km)
      call refreshbc(ws,im*jm,km)

      RETURN
      END    



      SUBROUTINE INVERSEX(US,VS,WS,TV,UO,UB,VB,WB,ALFXDT,GAMXDT,RHOXDT,IM,JM,KM)
c
      INCLUDE 'common.h'
c
      INTEGER IM,JM,KM
      REAL US(IM,JM,KM),VS(IM,JM,KM),WS(IM,JM,KM),TV(IM,JM,KM),UO(IM,JM,KM)
      REAL UB(IM,JM,KM),VB(IM,JM,KM),WB(IM,JM,KM)
c
      INTEGER I,J,K,II,I1,I2,NI,ID
      REAL AX(MMX),BX(MMX),CX(MMX),RX(MMX),U1D(MMX)
      REAL RPERB(NIMPLX,3,2),RDIRB(NIMPLX,3,2),RDIR(NIMPLX,3,2)
     &    ,BDIRB(NIMPLX,3,2),BNEUB(NIMPLX,3,2)
C
      REAL ALFXDT,GAMXDT,RHOXDT
      REAL TVIP,TVIM,UXP,UXM,QXP,QXM,QYP,QYM,QZP,QZM,ALPHA,BETA
      REAL AXV,BXV,CXV,AXC,BXC,CXC
      REAL VBC,WBC,RVBC(2),RWBC(2)
      REAL COEF,COEFX,RFLAGX,RFLAGCX
c
c------------------------------------------------------------------
c                                                    compute ustar
c------------------------------------------------------------------
      COEF=0.5*ALFXDT

      rflagx = real(implx)  ! 1 if the viscous terms are treated implicitly
      rflagcx = real(implcx)  ! 1 if the convective terms are treated implicitly
c
      rdir  = 1.0
      rperb = 0.0
      rdirb = 0.0
      bneub = 0.0
      bdirb = 0.0

C       if(myid == 0) then
C       write(6,*) "nimplx = ", nimplx
C       write(6,*) "implxvlim(1,1) = ", implxvlim(1,1)
C       write(6,*) "implxvlim(1,2) = ", implxvlim(1,2)
C       write(6,*) "ibu = ", ibu
C       endif


      do id=1,nimplx
        i1 = implxvlim(id,1)   ! left border of the implicit region
        i2 = implxvlim(id,2)   ! right border of the implicit region
        if(i1==ibu) then
          rdir(id,:,1) = 0.0
          if(itype(1)==500) then
            rperb(id,:,1) = 1.0 !Dirichlet for u
          else
            rperb(id,1,1) = 1.0
            do j=2,3
              if(ibc(j,1)==1) then !Dirichlet
                rdirb(id,j,1) = 2.0
                bdirb(id,j,1) =-1.0
              else !Neumann
                rdirb(id,j,1) = 0.0
                bdirb(id,j,1) = 1.0
              endif

c$$$         if(myrank==0) then
c$$$            write(*,*) "rperb(id=",id,",1,1)=",rperb(id,1,1)
c$$$            write(*,*) "rdirb(id=",id,",j=",j,",1)=",rdirb(id,j,1)
c$$$            write(*,*) "bdirb(id=",id,",j=",j,",1)=",bdirb(id,j,1)
c$$$         endif


            enddo
          endif
        endif
        if(i2==ix2) then
          rdir(id,:,2) = 0.0
          if(itype(2)==500) then
            rperb(id,:,2) = 1.0
            ub(im,:,:) = us(2,:,:)
            vb(im,:,:) = vs(2,:,:)
            wb(im,:,:) = ws(2,:,:)
          else
            rperb(id,1,2) = 1.0 !Dirichlet for u
            do j=2,3
              if(ibc(j,2)==1) then !Dirichlet
                rdirb(id,j,2) = 2.0
                bdirb(id,j,2) =-1.0
              else !Neumann
                rdirb(id,j,2) = 0.0
                bdirb(id,j,2) = 1.0
              endif
            enddo
          endif
        endif
      enddo

C c------------------------------------------------------------------
C c                                                   
C c------------------------------------------------------------------

      do k=kz1,kz2
      do j=jy1,jy2
      do id=1,nimplx
        i1 = implxulim(id,1)   ! left border of the implicit region
        i2 = implxulim(id,2)   ! right border of the implicit region
        ni = i2-i1+1
      do i=i1,i2
c..fill arrays for 3-diagonal system
        ii=i-i1+1
c...calculate total viscosity where needed
        tvip=2.*tv(i+1,j,k)+ru1
        tvim=2.*tv(i  ,j,k)+ru1
c...rhs
        rx(ii) = us(i,j,k)
c...coefx
        coefx = rflagx*coef*au(i)
c...u(i-1)
        axv =-coefx*ap(i  )*tvim*ru(i-1)/rp(i)
c...u(i+1)
        cxv =-coefx*ap(i+1)*tvip*ru(i+1)/rp(i+1)
c...u(i)
        bxv =-axv - cxv
c
c...Convective terms
c...coefx
        coefx = rflagcx*coef*au(i)/ru(i)
c
        qxm = (ru(i)*uo(i,j,k)+ru(i-1)*uo(i-1,j,k))*0.5
        qxp = (ru(i)*uo(i,j,k)+ru(i+1)*uo(i+1,j,k))*0.5
c
c...u(i-1,j,k)
        axc =-coefx*qxm
c        axc =-coefx*0.25*rp(i)*(uo(i,j,k)+uo(i-1,j,k))*2.0
c...u(i+1,j,k)
        cxc = coefx*qxp
c        cxc = coefx*0.25*rp(i+1)*(uo(i,j,k)+uo(i+1,j,k))*2.0
c...u(i,j,k)
        bxc = axc + cxc
c        bx(ii) = bx(ii) + coefx*0.25*(uo(i+1,j,k)-uo(i-1,j,k))*2.0

        ax(ii) = axv + axc
        cx(ii) = cxv + cxc
        bx(ii) = 1.0 + bxv + bxc

        uxp = (uo(i+1,j,k)+uo(i,j,k))*0.5
        uxm = (uo(i-1,j,k)+uo(i,j,k))*0.5

c...the RHS of the system is modified in the case of the implicit
c...treatment of the convective terms along X
        rx(ii) = rx(ii) + coefx*(qxp*uxp - qxm*uxm)

      enddo

c...solve the system
c      bx(1)  = bx(1)  + (bdirb(id,1,1)+bneub(id,1,1))*ax(1)
c      bx(ni) = bx(ni) + (bdirb(id,1,2)+bneub(id,1,2))*ax(ni)

C      if(myid == 0) then
C          write(6,*) "ni = ", ni
C          write(6,*) "i1 = ", i1
C          write(6,*) "i2 = ", i2
C      endif

      rx(1)  = rx(1)  - ax(1) *(us(i1-1,j,k)*rdir (id,1,1)
     &                        + ub(i1-1,j,k)*rperb(id,1,1))
      rx(ni) = rx(ni) - cx(ni)*(us(i2+1,j,k)*rdir (id,1,2)
     &                        + ub(i2+1,j,k)*rperb(id,1,2))
      call tridag(ax,bx,cx,rx,u1d,ni)
c
c...the provisional velocity is stored in US
      do i=i1,i2
        ii=i-i1+1
        us(i,j,k)=u1d(ii) !/ru(i)
      enddo
c
      enddo
      enddo
      enddo

c...REFRESHBC is called to allow the solution of the equations
c...for V and W using the provisional values of U evaluated above
      call refreshbc(us,im*jm,km)

c------------------------------------------------------------------
c                                                    compute vstar
c------------------------------------------------------------------

      do k=kz1,kz2
      do j=jbv,jev
      do id=1,nimplx
        i1 = implxvlim(id,1)   ! left border of the implicit region
        i2 = implxvlim(id,2)   ! right border of the implicit region
        ni = i2-i1+1
      do i=i1,i2
c..fill arrays for 3-diagonal system
        ii=i-i1+1
c...Viscous diffusion terms
c...calculate total viscosity where needed
        tvip = 0.25*(tv(i,j+1,k)+tv(i,j,k)+tv(i+1,j+1,k)+tv(i+1,j,k))+ru1
        tvim = 0.25*(tv(i,j+1,k)+tv(i,j,k)+tv(i-1,j+1,k)+tv(i-1,j,k))+ru1
c...rhs
        rx(ii)=vs(i,j,k)
c...coefx
        coefx =rflagx*coef*av(i)/rp(i)
c...v(i-1)
        axv =-coefx*au(i-1)*tvim*ru(i-1)
c...v(i+1)
        cxv =-coefx*au(i  )*tvip*ru(i  )
c...v(i)
        bxv =-axv-cxv
c
c...Convective terms
        coefx = rflagcx*0.5*coef*av(i)/rp(i)
c
        qyp=ru(i  )*(us(i  ,j+1,k)+us(i  ,j,k))*0.5
        qym=ru(i-1)*(us(i-1,j+1,k)+us(i-1,j,k))*0.5
c
c...v(i-1)
        axc =-coefx*qym
c...v(i+1)
        cxc = coefx*qyp
c...v(i)
        bxc = axc + cxc
c
        ax(ii) = axv + axc
        cx(ii) = cxv + cxc
        bx(ii) = 1.0 + bxv + bxc
      enddo
c...solve the system
      bx(1)  = bx(1)  + (bdirb(id,2,1)+bneub(id,2,1))*ax(1)
      bx(ni) = bx(ni) + (bdirb(id,2,2)+bneub(id,2,2))*cx(ni)

      rx(1)  = rx(1)  - ax(1) *(vs(i1-1,j,k)*rdir (id,2,1)
     &                        + vb(i1-1,j,k)*rperb(id,2,1)
     &       + 0.5*(vb(i1,j,k)+vb(i1-1,j,k))*rdirb(id,2,1))
              
      rx(ni) = rx(ni) - cx(ni)*(vs(i2+1,j,k)*rdir (id,2,2)
     &                        + vb(i2+1,j,k)*rperb(id,2,2)
     &       + 0.5*(vb(i2,j,k)+vb(i2+1,j,k))*rdirb(id,2,2))
      call tridag(ax,bx,cx,rx,u1d,ni)
c
c...the provisional velocity is stored in VS
      do i=i1,i2
        ii=i-i1+1
        vs(i,j,k)=u1d(ii)
      enddo

      enddo
      enddo
      enddo
c------------------------------------------------------------------
c                                                    compute wstar
c------------------------------------------------------------------
      do k=kbw,kew
      do j=jy1,jy2
      do id=1,nimplx
        i1 = implxvlim(id,1)   ! left border of the implicit region
        i2 = implxvlim(id,2)   ! right border of the implicit region
        ni = i2-i1+1
      do i=i1,i2
c..fill arrays for 3-diagonal system
        ii=i-i1+1
c...Viscous diffusion terms
c...calculate total viscosity where needed
        tvip = 0.25*(tv(i,j,k+1)+tv(i,j,k)+tv(i+1,j,k+1)+tv(i+1,j,k))+ru1
        tvim = 0.25*(tv(i,j,k+1)+tv(i,j,k)+tv(i-1,j,k+1)+tv(i-1,j,k))+ru1
c...rhs
        rx(ii)=ws(i,j,k)
c...coefx
        coefx = rflagx*coef*aw(i)/rp(i)
c...w(i-1)
        axv =-coefx*au(i-1)*tvim*ru(i-1)
c...w(i+1)
        cxv =-coefx*au(i  )*tvip*ru(i  )
c...w(i)
        bxv =-axv - cxv
c...Convective terms
        coefx = rflagcx*0.5*coef*aw(i)/rp(i)
c
        qzp=ru(i  )*(us(i  ,j,k+1)+us(i  ,j,k))*0.5
        qzm=ru(i-1)*(us(i-1,j,k+1)+us(i-1,j,k))*0.5
c...w(i-1)
        axc =-coefx*qzm
c...w(i+1)
        cxc = coefx*qzp
c...w(i)
        bxc = axc + cxc

        ax(ii) = axv + axc
        cx(ii) = cxv + cxc
        bx(ii) = 1.0 + bxv + bxc
      enddo
c...solve the system
c
      
       if(myid==0) then
C          write(6,*) "id = ", id
C          write(6,*) "rdir(id,3,1) = ", rdir(id,3,1)
C          write(6,*) "rdir(id,3,2) = ", rdir(id,3,2)

C          write(6,*) "rperb(id,3,1) = ", rperb(id,3,1)
C          write(6,*) "rperb(id,3,2) = ", rperb(id,3,2)

C          write(6,*) "rdirb(id,3,1) = ", rdirb(id,3,1)
C          write(6,*) "rdirb(id,3,2) = ", rdirb(id,3,2)

C          write(6,*) "bdirb(id,3,1) = ", bdirb(id,3,1)
C          write(6,*) "bdirb(id,3,2) = ", bdirb(id,3,2)

C          write(6,*) "bneub(id,3,1) = ", bneub(id,3,1)
C          write(6,*) "bneub(id,3,2) = ", bneub(id,3,2)

c          write(6,*) "wb(i1-1) + wb(i1) = ", wb(i1-1,j,k) + wb(i1  ,j,k)

       endif



      bx(1)  = bx(1)  + (bdirb(id,3,1)+bneub(id,3,1))*ax(1)
      bx(ni) = bx(ni) + (bdirb(id,3,2)+bneub(id,3,2))*cx(ni)

      rx(1)  = rx(1)  - ax(1) *(ws(i1-1,j,k)*rdir (id,3,1)
     &                        + wb(i1-1,j,k)*rperb(id,3,1)
     &       + 0.5*(wb(i1,j,k)+wb(i1-1,j,k))*rdirb(id,3,1))

c$$$      if(myrank==0) then
c$$$          write(*,*) "rdir (id,3,1) = ", rdir (id,3,1)  ! << 0
c$$$          write(*,*) "rperb(id,3,1) = ", rperb(id,3,1)  ! << 0
c$$$          write(*,*) "rdirb(id,3,1) = ", rdirb(id,3,1)  ! << 2
c$$$       endif
              
      rx(ni) = rx(ni) - cx(ni)*(ws(i2+1,j,k)*rdir (id,3,2)
     &                        + wb(i2+1,j,k)*rperb(id,3,2)
     &       + 0.5*(wb(i2,j,k)+wb(i2+1,j,k))*rdirb(id,3,2))
      call tridag(ax,bx,cx,rx,u1d,ni)
c
c...the provisional velocity is stored in WS
      do i=i1,i2
        ii=i-i1+1
        ws(i,j,k)=u1d(ii)
      enddo

      enddo
      enddo
      enddo

      call refreshbc(vs,im*jm,km)
      call refreshbc(ws,im*jm,km)

c
      RETURN
      END    

      SUBROUTINE INVERSEX_PER(US,VS,WS,TV,UO,ALFXDT,GAMXDT,RHOXDT,IM,JM,KM)
c
      INCLUDE 'common.h'
c
      INTEGER IM,JM,KM
      REAL US(IM,JM,KM),VS(IM,JM,KM),WS(IM,JM,KM),TV(IM,JM,KM),UO(IM,JM,KM)
c
      INTEGER I,J,K,II
      REAL AX(MMX),BX(MMX),CX(MMX),RX(MMX),U1D(MMX)
C
      REAL ALFXDT,GAMXDT,RHOXDT
      REAL TVIP,TVIM,UXP,UXM,QXP,QXM,QYP,QYM,QZP,QZM,ALPHA,BETA,AXV,BXV,CXV,AXC,BXC,CXC
      REAL COEF,COEFX,RFLAGX,RFLAGCX
c
c------------------------------------------------------------------
c                                                    compute ustar
c------------------------------------------------------------------
      COEF=0.5*ALFXDT

      rflagx = real(implx)   ! if 1 the viscous terms are treated implicitly
      rflagcx = real(implcx)   ! if 1 the convective terms are treated implicitly

      do k=kz1,kz2
      do j=jy1,jy2
      do i=ibu,ieu   ! no decomposition of the computational domain for the implicit treatment
c..fill arrays for 3-diagonal system
        ii=i-1
c...calculate total viscosity where needed
        tvip=2.*tv(i+1,j,k)+ru1
        tvim=2.*tv(i  ,j,k)+ru1
c...rhs
        rx(ii) = us(i,j,k)
c...coefx
        coefx = rflagx*coef*au(i)
c...u(i-1)
!!!!!!        axv =-coefx*ap(i  )*tvim/rp(i  )
        axv =-coefx*ap(i  )*tvim*ru(i-1)/rp(i)
c...u(i+1)
!!!!!!        cxv =-coefx*ap(i+1)*tvip/rp(i+1)
        cxv =-coefx*ap(i+1)*tvip*ru(i+1)/rp(i+1)
c...u(i)
        bxv =-axv - cxv
c
c...Convective terms
c...coefx
        coefx = rflagcx*coef*au(i)/ru(i)
c
        qxm = (ru(i)*uo(i,j,k)+ru(i-1)*uo(i-1,j,k))*0.5    !!!!!!
        qxp = (ru(i)*uo(i,j,k)+ru(i+1)*uo(i+1,j,k))*0.5    !!!!!!
c
c...u(i-1,j,k)
!!!!!!        axc =-coefx*0.25*rp(i)*(uo(i,j,k)+uo(i-1,j,k))*2.0
        axc =-coefx*qxm
c...u(i+1,j,k)
!!!!!!        cxc = coefx*0.25*rp(i+1)*(uo(i,j,k)+uo(i+1,j,k))*2.0
        cxc = coefx*qxp
c...u(i,j,k)
c        bx(ii) = bx(ii) + coefx*0.25*(uo(i+1,j,k)-uo(i-1,j,k))*2.0
        bxc = axc + cxc

        ax(ii) = axv + axc
        cx(ii) = cxv + cxc
!!!!!!        bx(ii) = 1.0/ru(i) + bxv + bxc
        bx(ii) = 1.0 + bxv + bxc

        uxp = (uo(i+1,j,k)+uo(i,j,k))*0.5
        uxm = (uo(i-1,j,k)+uo(i,j,k))*0.5

c...the RHS of the system is modified in the case of the implicit
c...treatment of the convective terms along X
!!!!!!        rx(ii) = rx(ii) + coefx*(rp(i+1)*uxp*uxp - rp(i)*uxm*uxm)
        rx(ii) = rx(ii) + coefx*(qxp*uxp - qxm*uxm)

      enddo

c...solve the system
c      For periodic boundary conditions
      beta  =ax(1)
      alpha =cx(ieu-1)
      call cyclic(ax,bx,cx,alpha,beta,rx,u1d,ieu-1)
c
c...the provisional velocity is stored in US
      do i=ibu,ieu
        ii=i-1
!!!!!!        us(i,j,k)=u1d(ii)/ru(i)
        us(i,j,k)=u1d(ii)
      enddo
c
      enddo
      enddo

c...REFRESHBC is called to allow the solution of the equations
c...for V and W using the provisional values of U evaluated above
      call refreshbc(us,im*jm,km)

c------------------------------------------------------------------
c                                                    compute vstar
c------------------------------------------------------------------
      do k=kz1,kz2
      do j=jbv,jev
      do i=ix1,ix2
c..fill arrays for 3-diagonal system
        ii=i-1
c...Viscous diffusion terms
c...calculate total viscosity where needed
        tvip=0.25*(tv(i,j  ,k)+tv(i+1,j  ,k)
     %       +     tv(i,j+1,k)+tv(i+1,j+1,k))+ru1
        tvim=0.25*(tv(i,j  ,k)+tv(i-1,j  ,k)
     %       +     tv(i,j+1,k)+tv(i-1,j+1,k))+ru1
c...rhs
        rx(ii)=vs(i,j,k)
c...coefx
        coefx =rflagx*coef*av(i)/rp(i)
c...v(i-1)
        axv =-coefx*au(i-1)*tvim*ru(i-1)
c...v(i+1)
        cxv =-coefx*au(i  )*tvip*ru(i  )
c...v(i)
        bxv =-axv-cxv
c
c...Convective terms
        coefx = rflagcx*0.5*coef*av(i)/rp(i)
c
        qyp=ru(i  )*(us(i  ,j+1,k)+us(i  ,j,k))*0.5
        qym=ru(i-1)*(us(i-1,j+1,k)+us(i-1,j,k))*0.5
c
c...v(i-1)
        axc =-coefx*qym
c...v(i+1)
        cxc = coefx*qyp
c...v(i)
        bxc = axc + cxc
c
        ax(ii) = axv + axc
        cx(ii) = cxv + cxc
        bx(ii) = 1.0 + bxv + bxc
      enddo
c...solve the system
c      Periodic boundary conditions
      beta  = ax(1)
      alpha = cx(ieu-1)
      call cyclic(ax,bx,cx,alpha,beta,rx,u1d,ieu-1)
c
c...the provisional velocity is stored in VS
      do i=ix1,ix2
        ii=i-1
        vs(i,j,k)=u1d(ii)
      enddo
c
      enddo
      enddo
c------------------------------------------------------------------
c                                                    compute wstar
c------------------------------------------------------------------
      do k=kbw,kew
      do j=jy1,jy2
      do i=ix1,ix2
c..fill arrays for 3-diagonal system
        ii=i-1
c...Viscous diffusion terms
c...calculate total viscosity where needed
        tvip=0.25*(tv(i  ,j,k)+tv(i  ,j,k+1)
     %       +     tv(i+1,j,k)+tv(i+1,j,k+1))+ru1
        tvim=0.25*(tv(i  ,j,k)+tv(i  ,j,k+1)
     %       +     tv(i-1,j,k)+tv(i-1,j,k+1))+ru1
c...rhs
        rx(ii)=ws(i,j,k)
c...coefx
        coefx = rflagx*coef*aw(i)/rp(i)
c...w(i-1)
        axv =-coefx*au(i-1)*tvim*ru(i-1)
c...w(i+1)
        cxv =-coefx*au(i  )*tvip*ru(i  )
c...w(i)
        bxv =-axv - cxv
c...Convective terms
        coefx = rflagcx*0.5*coef*aw(i)/rp(i)
c
        qzp=ru(i  )*(us(i  ,j,k+1)+us(i  ,j,k))*0.5
        qzm=ru(i-1)*(us(i-1,j,k+1)+us(i-1,j,k))*0.5
c...w(i-1)
        axc =-coefx*qzm
c...w(i+1)
        cxc = coefx*qzp
c...w(i)
        bxc = axc + cxc

        ax(ii) = axv + axc
        cx(ii) = cxv + cxc
        bx(ii) = 1.0 + bxv + bxc
      enddo
c...solve the system
c      Periodic boundary conditions
      beta  = ax(1)
      alpha = cx(ieu-1)
      call cyclic(ax,bx,cx,alpha,beta,rx,u1d,ieu-1)
c
c...the provisional velocity is stored in WS
      do i=ix1,ix2
        ii=i-1
        ws(i,j,k)=u1d(ii)
      enddo
c
      enddo
      enddo

      call refreshbc(vs,im*jm,km)
      call refreshbc(ws,im*jm,km)

c
      RETURN
      END    


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
      SUBROUTINE DIVERGENCE(US,VS,WS,DIV,IM,JM,KM)
*                                                                 *
*          Computes the divergence for an AB step                 *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      INCLUDE 'common.h'

      INTEGER IM,JM,KM
      REAL US(IM,JM,KM),VS(IM,JM,KM),WS(IM,JM,KM),DIV(IM,JM,KM)

      INTEGER I,J,K
*
**** Compute divergence
*
c      DO K=KZ1,KZ2
c      DO J=JY1,JY2
c      DO I=IX1,IX2
c        DIV(I,J,K) = AP(I)*(RU(I)*US(I,J,K)-RU(I-1)*US(I-1,J,K))/RP(I)
c      ENDDO
c      ENDDO
c      ENDDO
c      CALL IOSCALAR4('dudx.res',DIV,IM,JM,KM,1,0.0)

c      DO K=KZ1,KZ2
c      DO J=JY1,JY2
c      DO I=IX1,IX2
c        DIV(I,J,K) = CP(K)*(      WS(I,J,K)-        WS(I,J,K-1))
c      ENDDO
c      ENDDO
c      ENDDO
c      CALL IOSCALAR4('dwdz.res',DIV,IM,JM,KM,1,0.0)

      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
        DIV(I,J,K) = (AP(I)*(RU(I)*US(I,J,K)-RU(I-1)*US(I-1,J,K))/RP(I)
     &       +        BP(J)*(      VS(I,J,K)-        VS(I,J-1,K))/RP(I)
     &       +        CP(K)*(      WS(I,J,K)-        WS(I,J,K-1)))        
      ENDDO
      ENDDO
      ENDDO

c      I=2
c      J=2
c      K=6
c      IF(MYRANK.EQ.0) THEN
c          WRITE(6,'(15(1X,E16.8))') DIV(I,J,K),US(I,J,K),US(I-1,J,K)
c     &     ,VS(I,J,K),VS(I,J-1,K),WS(I,J,K),WS(I,J,K-1)
c     &     ,AP(I)*(RU(I)*US(I,J,K)-RU(I-1)*US(I-1,J,K))/RP(I)
c     &     ,BP(J)*(      VS(I,J,K)-        VS(I,J-1,K))/RP(I)
c     &     ,CP(K)*(      WS(I,J,K)-        WS(I,J,K-1))
c     &     ,AP(I),BP(J),CP(K),RU(I-1),RP(I)
c
c      ENDIF
c      I=11
c      J=2
c      K=6
c      IF(MYRANK.EQ.0) THEN
c          WRITE(6,'(15(1X,E16.8))') DIV(I,J,K),US(I,J,K),US(I-1,J,K)
c     &     ,VS(I,J,K),VS(I,J-1,K),WS(I,J,K),WS(I,J,K-1)
c     &     ,AP(I)*(RU(I)*US(I,J,K)-RU(I-1)*US(I-1,J,K))/RP(I)
c     &     ,BP(J)*(      VS(I,J,K)-        VS(I,J-1,K))/RP(I)
c     &     ,CP(K)*(      WS(I,J,K)-        WS(I,J,K-1))
c     &     ,AP(I),BP(J),CP(K),RU(I-1),RP(I)
c
c      ENDIF
C
      RETURN
      END

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * *
*                                                                    *
      SUBROUTINE CORRECT(US,VS,WS,UO,VO,WO,DP,IM,JM,KM,ALFXDT)
*                                                                    *
*                Performs an AB corector step                        *
*                                                                    *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * *

      INCLUDE 'common.h'

      INTEGER IM,JM,KM
      REAL UO(IM,JM,KM),VO(IM,JM,KM),WO(IM,JM,KM),
     &     US(IM,JM,KM),VS(IM,JM,KM),WS(IM,JM,KM),
     &     DP(IM,JM,KM)
      REAL ALFXDT

      INTEGER I,J,K
*
**** Update u values
*
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IBU,IEU
        UO(I,J,K) = US(I,J,K)-AU(I)*(DP(I+1,J,K)-DP(I,J,K))*ALFXDT
      ENDDO
      ENDDO
      ENDDO
*
**** Update v values
*
      DO K=KZ1,KZ2
      DO J=JBV,JEV
      DO I=IX1,IX2
        VO(I,J,K) = VS(I,J,K)-BV(J)*(DP(I,J+1,K)-DP(I,J,K))*ALFXDT/RP(I)
      ENDDO
      ENDDO
      ENDDO
*
**** Update w values
*
      DO K=KBW,KEW
      DO J=JY1,JY2
      DO I=IX1,IX2
        WO(I,J,K) = WS(I,J,K)-CW(K)*(DP(I,J,K+1)-DP(I,J,K))*ALFXDT        
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*                                                                 *
*           Predictor step  :   USTARAB                           *
*           Divergence      :   DIVAB                             *
*           Velocity corect.:   CORAB                             *
*                                                                 *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


      SUBROUTINE PREDICTORD(UO,VO,WO,P,UA,VA,WA,UB,VB,WB,
     &     US,VS,WS,DENS,XU,XC,YC,YV,ZWG,ZCG,ALFXDT,GAMXDT,RHOXDT,IM,JM,KM,KMG,TLEVEL,
     &     CLOCK,NCLOCK)

      use density_bg
      use spnge
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'

      
      INTEGER IM,JM,KM,dir,face,var,err
      REAL UO(IM,JM,KM),VO(IM,JM,KM),WO(IM,JM,KM),P(IM,JM,KM),
     &     UA(IM,JM,KM),VA(IM,JM,KM),WA(IM,JM,KM),
     &     UB(IM,JM,KM),VB(IM,JM,KM),WB(IM,JM,KM),
     &     US(IM,JM,KM),VS(IM,JM,KM),WS(IM,JM,KM)
      REAL DENS(IM,JM,KM),XU(IM),XC(IM),YC(JM),YV(JM),ZWG(KM),ZCG(KM)
!       REAL dfxug(KMG),dfxwg(KMG),dfxul(KM),dfxwl(KM)

      INTEGER I,J,K
      REAL ALFXDT,GAMXDT,RHOXDT
      REAL USTAR,VSTAR,WSTAR
      REAL COEF,COEFD,TLEVEL
      REAL CF1,CF2,CF3,DF,DFC,DFD
      PARAMETER (DFC=2.,DFD=20.)
      INTEGER K1G,K2G,KMG,KL

      integer nclock
      real clock(nclock)
      real clocktemp,clocktemp1
      real tclock

!*** Compute provisional values for each momentum eq.
      
      
      COEF=0.5*ALFXDT
      COEFD=ALFXDT
      CF1=DFC/(DFD**2.)
      CF2=-2.*DFD*CF1
      CF3=DFC
!------------------------------------------------------------------
!                                                    compute ustar
!------------------------------------------------------------------
      clocktemp1 = tclock()
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IBU,IEU
        USTAR = UO(I,J,K)
     &       +GAMXDT*UB(I,J,K)
     &       +RHOXDT*UA(I,J,K)
     &       -ALFXDT*AU(I)*(P(I+1,J,K)-P(I,J,K))
     &       +COEF*US(I,J,K)
     &       -ALFXDT*DPDX



       IF(IDENS.EQ.1) then
           USTAR=USTAR-COEFD*0.50d0*grav*sin(yc(j))*
     &           ((DENS(I+1,J,K)-DENS_BG(I+1,J))+(DENS(I,J,K)-DENS_BG(I,J)))
        endif


! (dens(i+1,j,k)+dens(i,j,k))

        !if(myid==0) write(6,*) "grav = ", grav
!	IF(IDENS.EQ.1) USTAR=USTAR-COEFD*0.50d0*grav*(DENS(I+1,J,K)+DENS(I,J,K))

! UO is the velocity at the time level L-1
! UB is the RHS (explicit part) at the time level L-1
! UA is the RHS (explicit part) at the time level L-2
! P is the pressure (for the pressure gradient) at the time level L-1
! US is the semi-explicit part for the terms treated by C.N.

	
        UA(I,J,K) = UB(I,J,K)
        UB(I,J,K) = US(I,J,K)
        US(I,J,K) = USTAR
! the RHS at the time level L-1 is stored in UA
! the C.N. terms are stored in UB
! the provisional solution is stored in US
      ENDDO
      ENDDO
      ENDDO


       var=1 
	
	      

       DO K=KZ1,KZ2
       DO J=JY1,JY2
       DO I=IEU-vspngx1,IEU

         US(I,J,K) = US(I,J,K) + dfxu(i)*(0.-US(I,J,K))*ALFXDT*1.0

       ENDDO
       ENDDO
       ENDDO


       IF (ivspngl==1) THEN

       DO K=KZ1,KZ2
       DO J=JY1,JY2
       DO I=IBU,IEU

         US(I,J,K) = US(I,J,K) + dfxul(K)*(0.-US(I,J,K))*ALFXDT*1.0

       ENDDO
       ENDDO
       ENDDO

       ENDIF

      clocktemp1 = tclock() - clocktemp1
      clock(1) = clock(1) + clocktemp1

!------------------------------------------------------------------
!                                                    compute vstar
!------------------------------------------------------------------
      clocktemp1 = tclock()
      DO K=KZ1,KZ2
      DO J=JBV,JEV
      DO I=IX1,IX2
        VSTAR = VO(I,J,K)
     &       +GAMXDT*VB(I,J,K)
     &       +RHOXDT*VA(I,J,K)
     &       -ALFXDT*BV(J)*(P(I,J+1,K)-P(I,J,K))/RP(I)
     &       +COEF*VS(I,J,K)
     &       -ALFXDT*DPDY

        IF(IDENS.eq.1) then
          VSTAR=VSTAR-COEFD*0.50d0*grav*cos(yv(j))*
     & ((DENS(I,J+1,K)-DENS_BG(I,J+1))+(DENS(I,J,K)-DENS_BG(I,J)))
       endif

        VA(I,J,K) = VB(I,J,K)
        VB(I,J,K) = VS(I,J,K)
        VS(I,J,K) = VSTAR

      ENDDO
      ENDDO
      ENDDO

       var=2 

       DO K=KZ1,KZ2
       DO J=JBV,JEV
       DO I=IX2-vspngx1,IX2

         VS(I,J,K) = VS(I,J,K) + dfxw(i)*(0.-VS(I,J,K))*ALFXDT*1.0

       ENDDO
       ENDDO
       ENDDO

       
       IF (ivspngl==1) THEN

       DO K=KZ1,KZ2
       DO J=JY1,JY2
       DO I=IBU,IEU

         VS(I,J,K) = VS(I,J,K) + dfxwl(K)*(0.-VS(I,J,K))*ALFXDT*1.0

       ENDDO
       ENDDO
       ENDDO
 
       ENDIF

      clock(4) = clock(4) + tclock()-clocktemp1

!------------------------------------------------------------------
!                                                    compute wstar
!------------------------------------------------------------------
      clocktemp1 = tclock()
      DO K=KBW,KEW
      DO J=JY1,JY2
      DO I=IX1,IX2

        WSTAR = WO(I,J,K)
     &       +GAMXDT*WB(I,J,K)
     &       +RHOXDT*WA(I,J,K)
     &       -ALFXDT*CW(K)*(P(I,J,K+1)-P(I,J,K))
     &       +COEF*WS(I,J,K)
     &       -ALFXDT*DPDZ

        WA(I,J,K) = WB(I,J,K)
        WB(I,J,K) = WS(I,J,K)
        WS(I,J,K) = WSTAR
      ENDDO
      ENDDO
      ENDDO

        var=3 

       DO K=KBW,KEW
       DO J=JY1,JY2
       DO I=IX2-vspngx1,IX2

         WS(I,J,K) = WS(I,J,K) + dfxw(i)*(1.-WS(I,J,K))*ALFXDT*1.0

       ENDDO
       ENDDO
       ENDDO


       IF (ivspngl==1) THEN

       DO K=KBW,KEW 
       DO J=JY1,JY2
       DO I=IX1,IX2

        WS(I,J,K) = WS(I,J,K) + dfxwl(k)*(1.-WS(I,J,K))*ALFXDT*1.0

       ENDDO
       ENDDO
       ENDDO

       ENDIF

      clock(7) = clock(7) + tclock()-clocktemp1

      RETURN
      END

 
