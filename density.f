C
C -------------------------------------- A. Posa - 03/11/2014 ----
C
C     -------------------------------------- Karu
      
      SUBROUTINE DENSITY(DENS,UO,VO,TV,RHA,RHB,RS,XC,YC,ALFXDT,GAMXDT,RHOXDT,IM,JM,KM)
c
      USE density_bg
      USE spnge
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
      INTEGER IM,JM,KM,var,dir,face,err
      REAL XC(IM),YC(JM)
      REAL DENS(IM,JM,KM),RHA(IM,JM,KM),RHB(IM,JM,KM),UO(IM,JM,KM),VO(IM,JM,KM)
      REAL RS(IM,JM,KM),DUM1(IM,JM,KM),TV(IM,JM,KM)
c
      INTEGER I,J,K
      REAL ALFXDT,GAMXDT,RHOXDT,COEF
      REAL CF1,CF2,CF3,DF,DFC,DFD
      PARAMETER (DFC=2.,DFD=20.)
c
      DUM1=DENS
      COEF=0.5*ALFXDT
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
        DENS(I,J,K) = DENS(I,J,K)
     &       +GAMXDT*RHB(I,J,K)
     &       +RHOXDT*RHA(I,J,K)
     &       +COEF*RS(I,J,K)
!      &       -denP1*(0.50d0*(uo(i,j,k)+uo(i-1,j,k))*sin(yc(j))+0.50d0*(vo(i,j,k)+vo(i,j-1,k))*cos(yc(j)))
! 	if(idens.eq.1) then
! 	DENS(I,J,K)=DENS(I,J,K)-(0.50d0*sin(YC(J))*(UO(I,J,K)+UO(I-1,J,K))+0.50d0*cos(YC(J))*(VO(I,J+1,K)+VO(I,J-1,K)))*denP1
! 	endif
! 	if(i==2.and.k==5) write(*,*), RHA(I,J,K),RHB(I,J,K),i,j,k
c
c DENS is the density at the time level L-1 (updated at the time level L)
c RHB is the RHS (explicit part) at the time level L-1
c RHA is the RHS (explicit part) at the time level L-2
c
        RHA(I,J,K) = RHB(I,J,K)
c Karu adds ----------------------
        RHB(I,J,K) = RS(I,J,K)
c -------------------------------        
	RS(I,J,K)  = DENS(I,J,K)
	
c	DENS(I,J,K) = RS(I,J,K)
c the RHS at the time level L-1 is stored in RHA
      ENDDO
      ENDDO
      ENDDO
c Karu adds----------------------------------------------      
c.....predicted density from explicit step
      RHB = RS+0.5*ALFXDT*RHB
      
      CALL BOUNDARY_DENS(RHB,XC,YC,IM,JM,KM)
      CALL REFRESHBC(RHB,IM*JM,KM)

      ! Karu adds ... (Karu: I should be more organized)
      IF(imply==1 .AND. implcy==1) then          
         CALL INVERSEYD_NEW(RS,VO,TV,ALFXDT,GAMXDT,RHOXDT,IM,JM,KM)
         DENS = RS
      ELSEIF(imply==1) then         
         CALL INVERSEYD(RS,TV,ALFXDT,GAMXDT,RHOXDT,IM,JM,KM)
         DENS = RS         
      ENDIF
      
      IF(implx==1 .OR. implcx==1) then
         CALL INVERSEXD(RS,RHB,ALFXDT,GAMXDT,RHOXDT,IM,JM,KM)
         DENS = RS
      ENDIF


! Sponge density: inlet, outlet, radial

      IF (idspngl==1) then
     
       DO I=IX1,IX2
       DO J=JY1,JY2
       DO K=KZ1,KZ2
                  
         DENS(I,J,K) = DENS(I,J,K) + dfxdl(K)*(DENS_BG(I,J) - DENS(I,J,K))*ALFXDT*1.0

       ENDDO
       ENDDO
       ENDDO
      
      ENDIF  
     

       DO I=IX2-dspngx1,IX2
       DO J=JY1,JY2
       DO K=KZ1,KZ2
                  
         DENS(I,J,K) = DENS(I,J,K) + dfxd(I)*(DENS_BG(I,J) - DENS(I,J,K))*ALFXDT*1.0

       ENDDO
       ENDDO
       ENDDO

c ----------------------------------------------      
c ----------------------------------------------      
 	var=5 !u=Flow(:,:,:,5)
! !  	do dir=1,3
! !  	do face=1,2
!      	 call Sponge(DENS,DUM1,var,IM,JM,KM,err)
!  	enddo
!  	enddo
!       cf1=DFC/(DFD**2.)
!       cf2=-2.*DFD*cf1
!       cf3=DFC
!       DO K=KZ1,KZ2
!       DO J=JY1,JY2
!       DO I=IX1,IX1+20
!         df=cf1*real(i-ix1)**2.+cf2*real(i-ix1)+cf3
!         DENS(I,J,K) = DENS(I,J,K) + df*(denP1*XC(I)-DENS(I,J,K))*ALFXDT
!       ENDDO
!       DO I=IX2-20,IX2
!         df=cf1*real(ix2-i)**2.+cf2*real(ix2-i)+cf3
!         DENS(I,J,K) = DENS(I,J,K) + df*(denP1*XC(I)-DENS(I,J,K))*ALFXDT
!       ENDDO
!       ENDDO
!       ENDDO

      RETURN
      END

       SUBROUTINE INVERSEXD(RS,RHB,ALFXDT,GAMXDT,RHOXDT,IM,JM,KM)
c
       INCLUDE 'common.h'
c
       INTEGER IM,JM,KM
       REAL RS(IM,JM,KM),RHB(IM,JM,KM)
c
       INTEGER I,J,K,II,I1,I2,NI,ID
       REAL AX(MMX),BX(MMX),CX(MMX),RX(MMX),U1D(MMX)
       REAL RPERB(NIMPLX,3,2),RDIRB(NIMPLX,3,2),RDIR(NIMPLX,3,2)
     &     ,BDIRB(NIMPLX,3,2),BNEUB(NIMPLX,3,2)
C
       REAL ALFXDT,GAMXDT,RHOXDT
       REAL TVIP,TVIM
       REAL AXV,BXV,CXV,AXC,BXC,CXC
       REAL ALPHA,BETA
       REAL COEF,COEFX,RFLAGX,RFLAGCX

C c------------------------------------------------------------------
       COEF=0.5*ALFXDT

       rflagx = real(implx)  ! 1 if the viscous terms are treated implicitly
c
        rdir  = 1.0
        rperb = 0.0
        rdirb = 0.0
        bneub = 0.0
        bdirb = 0.0

       do id=1,nimplx
         i1 = implxvlim(id,1)   ! left border of the implicit region
         i2 = implxvlim(id,2)   ! right border of the implicit region

c         write(*,*) "myrank,id,nimplx,i1,i2,ibu=",myrank,id,nimplx,i1,i2,ibu

         if(i1==ibu) then
           rdir(id,:,1) = 0.0
           if(itype(1)==500) then
             rperb(id,:,1) = 1.0 
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
             enddo
           endif
         endif
         
         


        if(i2==ix2) then
           rdir(id,:,2) = 0.0
           if(itype(2)==500) then
c              rperb(id,:,2) = 1.0
c              rhb(im,:,:) = rs(2,:,:)
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
C c                                                    compute wstar
C c------------------------------------------------------------------

       do k=kbw,kew
       do j=jy1,jy2
       do id=1,nimplx
         i1 = implxvlim(id,1)   ! left border of the implicit region
         i2 = implxvlim(id,2)   ! right border of the implicit region
         ni = i2-i1+1

       do i=i1,i2
C c..fill arrays for 3-diagonal system
          ii=i-i1+1
c...Viscous diffusion terms
c...calculate total viscosity where needed
          tvip = ru1/prn
          tvim = ru1/prn
c... rhs
          rx(ii)=rs(i,j,k)
c...coefx
          coefx = rflagx*coef*aw(i)/rp(i)
c...w(i-1)
          axv =-coefx*au(i-1)*tvim*ru(i-1)
c...w(i+1)
          cxv =-coefx*au(i  )*tvip*ru(i  )
c...w(i)
          bxv =-axv - cxv

! march convective terms explicitly        
          ax(ii) = axv !+ axc
          cx(ii) = cxv !+ cxc
          bx(ii) = 1.0 + bxv !+ bxc
      enddo
      
c...solve the system
c Channel flow: Karu (Septemper 2015)
c bdirb(id,3,1) = -1, bneub(id,3,1) = 0
c bdirb(id,3,2) = -1, bneub(id,3,2) = 0
c  rdir(id,3,1) =  0
c  rdir(id,3,2) =  0
c rperb(id,3,1) =  0
c rperb(id,3,2) =  0

c rdirb(id,3,1) =  2
c rdirb(id,3,1) =  2

       bx(1)  = bx(1)  + (bdirb(id,3,1)+bneub(id,3,1))*ax(1)
       bx(ni) = bx(ni) + (bdirb(id,3,2)+bneub(id,3,2))*cx(ni)
       
      rx(1)  = rx(1)  - ax(1) * (rs(i1-1,j,k)*rdir (id,3,1)
     &                        + rhb(i1-1,j,k)*rperb(id,3,1)
     &       + 0.5*(rhb(i1,j,k)+rhb(i1-1,j,k))*rdirb(id,3,1))

c$$$      if(myrank==0) then
c$$$          write(*,*) "rhb(2,2,40) = ", rhb(2,2,40)
c$$$          write(*,*) "rhb(2,18,40) = ", rhb(2,18,40)
c$$$       endif

              
      rx(ni) = rx(ni) - cx(ni)* (rs(i2+1,j,k)*rdir (id,3,2)
     &                        + rhb(i2+1,j,k)*rperb(id,3,2)
     &       + 0.5*(rhb(i2,j,k)+rhb(i2+1,j,k))*rdirb(id,3,2))

      


      call tridag(ax,bx,cx,rx,u1d,ni)
c
c...the provisional velocity is stored in WS
       do i=i1,i2
         ii=i-i1+1
         rs(i,j,k)=u1d(ii)
       enddo

       enddo
       enddo
       enddo

      call refreshbc(rs,im*jm,km)
C c
      RETURN
      END    

      SUBROUTINE INVERSEYD_NEW(RS,VO,TV,ALFXDT,GAMXDT,RHOXDT,IM,JM,KM)
c Karu adds
      INCLUDE 'common.h'
c
      INTEGER IM,JM,KM
      REAL RS(IM,JM,KM),VO(IM,JM,KM),TV(IM,JM,KM)
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
      rflagy  = real(imply)     ! 1 if the viscous terms are treated implicitly
      rflagcy = real(implcy)   ! 1 if the convective terms are treated implicitly
c
c------------------------------------------------------------------
c                             If fully implicit compute vstar first!
c------------------------------------------------------------------
c

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
	tvjp = 0.5*(tv(i,j,k)+tv(i,j+1,k)) + ru1/prn
	tvjm = 0.5*(tv(i,j,k)+tv(i,j-1,k)) + ru1/prn
c...rhs
        ry(jj)=rs(i,j,k)
c
c...Viscous diffusion
c...coefy
        coefy  = rflagy*coef*bp(j)/rp(i)**2
c...w(j-1)
        ayv =-coefy*bv(j-1)*tvjm
c...w(j+1)
        cyv =-coefy*bv(j  )*tvjp
c...w(j)
        byv =-ayv - cyv
c
c...Convective term (the component V is known)
        vzp=vo(i,j  ,k)
        vzm=vo(i,j-1,k)
c...coefy
        coefy = rflagcy*0.5*coef*bp(j)/rp(i)
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
        rs(i,j,k)=v1d(jj)
      enddo

      enddo
      enddo
      enddo
      call refreshbc(rs,im*jm,km)

      RETURN
      END    
      
      SUBROUTINE INVERSEYD(RS,TV,ALFXDT,GAMXDT,RHOXDT,IM,JM,KM)
c
      INCLUDE 'common.h'
c
      INTEGER IM,JM,KM
      REAL RS(IM,JM,KM),TV(IM,JM,KM)
c
      INTEGER I,J,K,JJ,I1,I2,NI,ID
      REAL AY(MMY),BY(MMY),CY(MMY),RY(MMY),V1D(MMY)
C
      REAL ALFXDT,GAMXDT,RHOXDT
      REAL TVJP,TVJM,RFLAGY
      REAL COEF,COEFY,AYV,BYV,CYV,AYC,BYC,CYC
      REAL ALPHA,BETA

      REAL vxp,vxm,vzp,vzm,vyp,vym
*
**** Compute provisional values for each momentum eq.
*
      COEF=0.5*ALFXDT
c
      rflagy  = real(imply)  ! 1 if the viscous terms are treated implicitly
c
c------------------------------------------------------------------
c                             If fully implicit compute vstar first!
c------------------------------------------------------------------
c

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
	tvjp = 0.5*(tv(i,j,k)+tv(i,j+1,k)) + ru1/prn
	tvjm = 0.5*(tv(i,j,k)+tv(i,j+1,k)) + ru1/prn
c...rhs
        ry(jj)=rs(i,j,k)
c
c...Viscous diffusion
c...coefy
        coefy  = rflagy*coef*bp(j)/rp(i)**2
c...w(j-1)
        ayv =-coefy*bp(j-1)*tvjm
c...w(j+1)
        cyv =-coefy*bp(j  )*tvjp
c...w(j)
        byv =-ayv - cyv
c     
        ay(jj) = ayv 
        cy(jj) = cyv 
        by(jj) = 1.0  + byv 
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
        rs(i,j,k)=v1d(jj)
      enddo

      enddo
      enddo
      enddo
      call refreshbc(rs,im*jm,km)

      RETURN
      END    

