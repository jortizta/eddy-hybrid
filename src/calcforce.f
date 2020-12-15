C-----SUBROUTINE-calcfrc-----------------N. Beratlis  01/31/07------
C
      subroutine calcfrc(u,v,w,p,flagu,flagv,flagw,nx,ny,nz,x,y,z
     &     ,xc,yc,zc,fb,dt,icycle,time,filestem,cvlim)
c
c     PURPOSE: Calculates the forces on the immersed boundary using 
c              the control volume approach
c     
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
c.... Input/Output Arrays
      INTEGER NX,NY,NZ,IBD,ICYCLE
      REAL    DT,TIME
      REAL    FB(3)
      REAL    CVLIM(4)
      INTEGER FLAGW(NX,NY,NZ),FLAGV(NX,NY,NZ),FLAGU(NX,NY,NZ)
      REAL    U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ)
      REAL    X(NX),Y(NY),Z(NZ)
      REAL    XC(NX),YC(NY),ZC(NZ)
      CHARACTER filestem*(*)
c
      INTEGER I,J,K,ICOUNT,ICOUNTG
      INTEGER ICB,ICE,KCB,KCE,KBFLAG,KEFLAG
      REAL    TEMP(20),TEMPG
      REAL    ZLOC(2),XLOC(2)
      REAL    UC,USQ,VC,VSQ,WC,WSQ
      REAL    U1,V1,W1,UV,UW,VW
      REAL    DUDX,DUDY,DUDZ,DWDX,DWDY,DWDZ,DVDX,DVDY,DVDZ
      REAL    SXX,SYY,SZZ,SXY,SXZ,SYZ,PRES,TXX,TYY,TZZ,TXY,TXZ,TYZ,DA,DV
      REAL    DUDY1,DUDY2,DUDZ1,DUDZ2
      REAL    DVDX1,DVDX2,DVDZ1,DVDZ2
      REAL    DWDX1,DWDX2,DWDY1,DWDY2
      REAL    FBX,FBXC,FBXCG,FBY,FBYC,FBYCG,FBZ,FBZC,FBZCG
      REAL    FX(6),FXG(6),FZ(6),FZG(6),FY(6),FYG(6)
      REAL    Finlet(8),Foutlet(8),FGinlet(8),FGoutlet(8)
      REAL    Ftop(8),FGtop(8),Fbot(8),FGbot(8)

      xloc(1) = cvlim(1)
      xloc(2) = cvlim(2)
      zloc(1) = cvlim(3)
      zloc(2) = cvlim(4)

      
      if(icyl==1) then
        icb = 1
c        call closest(X,NX,0.0,ICB)
      else
        call closest(X,NX,xloc(1),ICB)
      endif
c      IF(ICB<=1) ICB=2
      CALL closest(X,NX,xloc(2),ICE)
c      write(6,'(A,3(1x,I4))') 'MYRANK=',MYRANK,ICB,ICE

      KCB=1
      KCE=0
      KBFLAG=0
      KEFLAG=0
c.....Locate inlet(KCB) and outlet(KCE).
      CALL LOCATE(Z,NZ-1,zloc(1),KCB)
      CALL LOCATE(Z,NZ-1,zloc(2),KCE)



      IF(KCB==0 .AND. KCE==0 ) THEN !CV Inlet and Outlet are in front of the domain.
        KBFLAG=0
        KEFLAG=0
        KCB=NZ !Do not integrate over Y-Z plane
        KCE=0 !Do not integrate over Y-Z plane
      ELSEIF(KCB==0 .AND. ( KCE>0.AND.KCE<NZ-1) ) THEN !CV inlet is ahead and CV outlet is inside the domain.
        KBFLAG=0 
        KEFLAG=1
        KCB=1 !For Y-Z plane
        IF(KCE==1) KCE=2
      ELSEIF( (KCB>0 .AND. KCB<NZ-1) .AND. (KCE>0 .AND. KCE<NZ-1) ) THEN !CV is inside the domain.
        KBFLAG=1
        KEFLAG=1
        KCB=KCB+1 
      ELSEIF(KCB==0 .AND. KCE==NZ-1) THEN !CV is larger than the domain.
        KBFLAG=0
        KEFLAG=0
        KCB=1 !Integrate over the entire Y-Z plane.
      ELSEIF( (KCB>0 .AND. KCB<NZ-1) .AND. KCE==NZ-1) THEN !CV inlet is inside but CV outlet is outside the domain.
        KBFLAG=1 
        KEFLAG=0
        KCB=KCB+1 !Integrate over the entire Y-Z plane ( KCB+1 <=K<= NZ-1 )
      ELSEIF( KCB==NZ-1 .AND. KCE==NZ-1 ) THEN !CV is behind the domain.
        KBFLAG=0
        KEFLAG=0
        KCB=NZ
        KCE=0
      ENDIF

      FBXC = 0.
      FBYC = 0.
      FBZC = 0.
      
      DO K=KCB+1,KCE
        DO I=ICB+1,ICE
          DO J=JY1,JY2
            dV = RP(I)/(BP(J)*AP(I)*CP(K))

            if(flagu(i,j,k)<=0 .OR. flagu(i-1,j,k)<=0) then
              uc = 0.5*(u(i,j,k)+u(i-1,j,k))
              FBXC = FBXC + uc*dV
            endif

            if(flagv(i,j,k)<=0 .OR. flagv(i,j-1,k)<=0) then
              vc = 0.5*(v(i,j,k)+v(i,j-1,k))
              FBYC = FBYC + vc*dV
            endif

            if(flagw(i,j,k)<=0 .OR. flagw(i,j,k-1)<=0) then
              wc = 0.5*(w(i,j,k)+w(i,j,k-1))
              FBZC = FBZC + wc*dV
            endif

          ENDDO
        ENDDO
      ENDDO

      CALL MPI_REDUCE(FBXC,FBXCG,1,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(FBYC,FBYCG,1,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(FBZC,FBZCG,1,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)

      IF(MYRANK.EQ.0) THEN
        IF(ICYCLE.EQ.0) THEN 
          FBX = 0.0
          FBY = 0.0
          FBZ = 0.0
        ELSE
          FBX = (FBXCG-FB(1))/DT
          FBY = (FBYCG-FB(2))/DT
          FBZ = (FBZCG-FB(3))/DT
        ENDIF
        FB(1) = FBXCG
        FB(2) = FBYCG
        FB(3) = FBZCG
      ENDIF
c
c---  X-Y planes
c
      FX = 0.
      FY = 0.
      FZ = 0.
      Finlet = 0.
      Foutlet = 0.
      Ftop = 0.
      Fbot = 0.
      temp = 0.

c      open(unit=10,file='p_inlet'//index(icycle)//'.dat'
c     &     ,form='formatted')
c      write(10,*) 'VARIABLES = "x", "y", "p"'
c      write(10,*) 'ZONE I=',(ice-icb)*(jy2-jy1+1)
c     &     ,', DATAPACKING=POINT'

c      open(unit=20,file='p_outlet'//index(icycle)//'.dat'
c     &     ,form='formatted')
c      write(20,*) 'VARIABLES = "x", "y", "p"'
c      write(20,*) 'ZONE I=',(ice-icb)*(jy2-jy1+1)
c     &     ,', DATAPACKING=POINT'

c 
c---  X-Y planes

      if(icyl==1) then


        DO I=ICB+1,ICE
          DO J=JY1,JY2
            dA = RP(I)/BP(J)/AP(I)

c            write(10,*) xc(i),yc(j),p(i,j,kcb)

            IF(KBFLAG>0) THEN
c...Inlet boundary.
              k=kcb

              !r-component
              DWDX = (w(i+1,j,k)-w(i-1,j,k))*aw(i)/2.
              DUDZ1 = (u(i  ,j,k+1)-u(i  ,j,k))*cw(k)
              DUDZ2 = (u(i-1,j,k+1)-u(i-1,j,k))*cw(k)
              DUDZ = 0.5*(DUDZ1+DUDZ2)
              SXZ = DUDZ+DWDX
              TXZ = RU1*SXZ
              U1 = 0.25*(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1))
              UW = U1*W(i,j,k)

              !Decompose r into x-component
              FX(1) = FX(1) + (UW - TXZ)*(-dA)*cos(yc(j))
              Finlet(1) = Finlet(1) + UW*(-dA)*cos(yc(j))
              Finlet(2) = Finlet(2) - TXZ*(-dA)*cos(yc(j))

              !Decompose r into y-component
              FY(1) = FY(1) + (UW - TXZ)*(-dA)*sin(yc(j))
              Finlet(3) = Finlet(3) + UW*(-dA)*sin(yc(j))
              Finlet(4) = Finlet(4) - TXZ*(-dA)*sin(yc(j))

              !theta-component
              DWDY = (w(i,j+1,k)-w(i,j-1,k))*bw(j)*0.5/rp(i)
              DVDZ1 = (v(i,j  ,k+1)-v(i,j  ,k  ))*cw(k)
              DVDZ2 = (v(i,j-1,k+1)-v(i,j-1,k  ))*cw(k)
              DVDZ = 0.5*(DVDZ1+DVDZ2)
              SYZ = DWDY+DVDZ
              TYZ = RU1*SYZ
              V1 = 0.25*(v(i,j,k)+v(i,j,k+1)+v(i,j-1,k+1)+v(i,j-1,k))
              VW = V1*W(i,j,k)

              !Decompose theta into x-component
              FX(1) = FX(1) + (VW - TYZ)*(-dA)*(-sin(yc(j)))
              Finlet(1) = Finlet(1) + VW*(-dA)*(-sin(yc(j)))
              Finlet(2) = Finlet(2) - TYZ*(-dA)*(-sin(yc(j)))

              !Decompose theta into y-component
              FY(1) = FY(1) + (VW - TYZ)*(-dA)*cos(yc(j))
              Finlet(3) = Finlet(3) + VW*(-dA)*cos(yc(j))
              Finlet(4) = Finlet(4) - TYZ*(-dA)*cos(yc(j))

              !z-component
              DWDZ    = (W(I,J,KCB+1)-W(I,J,KCB-1))*CW(KCB)/2.
              SZZ     = 2.*RU1*DWDZ
              PRES    = 0.5*(P(I,J,KCB+1)+P(I,J,KCB))
              TZZ     = PRES-SZZ
              WSQ     = W(I,J,KCB)**2.
              FZ(1)  = FZ(1) + (WSQ + TZZ)*(-dA)
              Finlet(5) = Finlet(5) + PRES*(-dA)
              Finlet(6) = Finlet(6) - SZZ*(-dA)
              Finlet(7) = Finlet(7) + WSQ*(-dA)

            ENDIF

            IF(KEFLAG>0) THEN
c...Outlet boundary.
              k=kce

c              write(20,*) xc(i),yc(j),p(i,j,k)

              !x-component
              DWDX = (w(i+1,j,k)-w(i-1,j,k))*aw(i)/2.
              DUDZ1 = (u(i  ,j,k+1)-u(i  ,j,k))*cw(k)
              DUDZ2 = (u(i-1,j,k+1)-u(i-1,j,k))*cw(k)
              DUDZ = 0.5*(DUDZ1+DUDZ2)
              SXZ = DUDZ+DWDX
              TXZ = RU1*SXZ
              U1 = 0.25*(U(i,j,k)+U(i,j,k+1)+U(i-1,j,k)+U(i-1,j,k+1))
              UW = U1*W(i,j,k)

              !Decompose r into x-component
              FX(2) = FX(2) + (UW - TXZ)*dA*cos(yc(j))
              Foutlet(1) = Foutlet(1) + UW*dA*cos(yc(j))
              Foutlet(2) = Foutlet(2) - TXZ*dA*cos(yc(j))

              !Decompose r into y-component
              FY(2) = FY(2) + (UW - TXZ)*dA*sin(yc(j))
              Foutlet(3) = Foutlet(3) + UW*dA*sin(yc(j))
              Foutlet(4) = Foutlet(4) - TXZ*dA*sin(yc(j))

              !theta-component
              DWDY = (w(i,j+1,k)-w(i,j-1,k))*bw(j)*0.5/rp(i)
              DVDZ1 = (v(i,j  ,k+1)-v(i,j  ,k  ))*cw(k)
              DVDZ2 = (v(i,j-1,k+1)-v(i,j-1,k  ))*cw(k)
              DVDZ = 0.5*(DVDZ1+DVDZ2)
              SYZ = DWDY+DVDZ
              TYZ = RU1*SYZ
              V1 = 0.25*(v(i,j,k)+v(i,j,k+1)+v(i,j-1,k+1)+v(i,j-1,k))
              VW = V1*W(i,j,k)

              !Decompose theta into x-component
              FX(2) = FX(2) + (VW - TYZ)*(dA)*(-sin(yc(j)))
              Foutlet(1) = Foutlet(1) + VW*(dA)*(-sin(yc(j)))
              Foutlet(2) = Foutlet(2) - TYZ*(dA)*(-sin(yc(j)))

              !Decompose theta into y-component
              FY(2) = FY(2) + (VW - TYZ)*(dA)*cos(yc(j))
              Foutlet(3) = Foutlet(3) + VW*(dA)*cos(yc(j))
              Foutlet(4) = Foutlet(4) - TYZ*(dA)*cos(yc(j))

              DWDZ    = (W(I,J,KCE+1)-W(I,J,KCE-1))*CW(KCE)/2.
              SZZ     = 2.*RU1*DWDZ
              PRES    = 0.5*(P(I,J,KCE+1)+P(I,J,KCE))
              TZZ     = PRES-SZZ
              WSQ     = W(I,J,KCE)**2.

              !z-component.
              FZ(2)  = FZ(2) + (WSQ + TZZ)*dA
              Foutlet(5) = Foutlet(5) + PRES*dA
              Foutlet(6) = Foutlet(6) -SZZ*dA
              Foutlet(7) = Foutlet(7) + WSQ*dA

            ENDIF

          ENDDO
        ENDDO

        i = ice
        DO K=KCB+1,KCE
          DO J=JY1,JY2
c...Top plane
            dA = RU(I)/BP(J)/CP(K)

            !r-component
            DUDX = (u(i+1,j,k)-u(i-1,j,k))*au(i)/2.0
            SXX = 2.0*ru1*dudx
            PRES = 0.5*(p(i+1,j,k)+p(i,j,k))
            TXX = PRES - SXX
            USQ = u(i,j,k)**2.

            !Decompose r into x-component
            FX(3) = FX(3) + (USQ + TXX)*dA*cos(yc(j))
            Ftop(1) = Ftop(1) + USQ*dA*cos(yc(j))
            Ftop(2) = Ftop(2) + PRES*dA*cos(yc(j))
            Ftop(3) = Ftop(3) - SXX*dA*cos(yc(j))

            !Decompose r into y-component
            FY(3) = FY(3) + (USQ + TXX)*dA*sin(yc(j))
            Ftop(4) = Ftop(4) + USQ*dA*sin(yc(j))
            Ftop(5) = Ftop(5) + PRES*dA*sin(yc(j))
            Ftop(6) = Ftop(6) - SXX*dA*sin(yc(j))

            !theta-component
            V1 = 0.25*(v(i,j,k)+v(i,j-1,k)+v(i+1,j-1,k)+v(i+1,j,k))
c            DUDY = (u(i,j+1,k)-u(i,j-1,k))*bv(j)/(2.0) - V1/ru(i)
            DUDY = (u(i,j+1,k)-u(i,j-1,k))*bv(j)/(2.0*ru(i)) - V1/ru(i)
            DVDX1 = au(i)*(v(i+1,j  ,k)-v(i,j  ,k))
            DVDX2 = au(i)*(v(i+1,j-1,k)-v(i,j-1,k))
            DVDX = 0.5*(DVDX1+DVDX2)
            SXY = RU1*(DUDY+DVDX)
            UV = u(i,j,k)*v1

            !Decompose theta into x-component
            FX(3) = FX(3) + (UV - SXY)*(dA)*(-sin(yc(j)))
            Ftop(1) = Ftop(1) + UV*(dA)*(-sin(yc(j)))
            Ftop(2) = Ftop(2) - SXY*(dA)*(-sin(yc(j)))

            !Decompose theta into y-component
            FY(3) = FY(3) + (UV - SXY)*dA*cos(yc(j))
            Ftop(4) = Ftop(4) + UV*(dA)*cos(yc(j))
            Ftop(6) = Ftop(6) - SXY*(dA)*cos(yc(j))            

            !z-component
            DUDZ  = (U(I  ,J,K+1) - U(I,J,K-1))*CU(K)/2.
            DWDX1 = (W(I+1,J,K  ) - W(I,J,K  ))*AU(I)
            DWDX2 = (W(I+1,J,K-1) - W(I,J,K-1))*AU(I)
            DWDX  = 0.5*(DWDX1+DWDX2)
            W1    = 0.25*(W(I,J,K) + W(I,J,K-1)
     &           + W(I+1,J,K) + W(I+1,J,K-1) )
            UW    = U(I,J,K)*W1

            SXZ   = RU1*(DUDZ+DWDX)

            FZ(3)    = FZ(3) + (UW - SXZ)*dA
            Ftop(7) = Ftop(7) + UW*dA
            Ftop(8) = Ftop(8) - SXZ*dA

          ENDDO
        ENDDO

      
      else !Cartesian Coordinates
c
        DO I=ICB+1,ICE
          DO J=JY1,JY2
            dA = RP(I)/BP(J)/AP(I)

            IF(KBFLAG>0) THEN
c...Inlet boundary.
              k=kcb

c              write(10,*) xc(i),yc(j),p(i,j,k)

              !x-component
              DWDX = (w(i+1,j,k)-w(i-1,j,k))*aw(i)/2.
              DUDZ1 = (u(i  ,j,k+1)-u(i  ,j,k))*cw(k)
              DUDZ2 = (u(i-1,j,k+1)-u(i-1,j,k))*cw(k)
              DUDZ = 0.5*(DUDZ1+DUDZ2)
              SXZ = DUDZ+DWDX
              TXZ = RU1*SXZ
              U1 = 0.25*(U(i,j,k)+U(i,j,k+1)+U(i-1,j,k)+U(i-1,j,k+1))
              UW = U1*W(i,j,k)
              FX(1) = FX(1) + (UW - TXZ)*(-dA)
              Finlet(1) = Finlet(1) + UW*(-dA)
              Finlet(2) = Finlet(2) - TXZ*(-dA)

              !z-component
              DWDZ    = (W(I,J,KCB+1)-W(I,J,KCB-1))*CW(KCB)/2.
              SZZ     = 2.*RU1*DWDZ
              PRES    = 0.5*(P(I,J,KCB+1)+P(I,J,KCB))
              TZZ     = PRES-SZZ
              WSQ     = W(I,J,KCB)**2.
              FZ(1)  = FZ(1) + (WSQ + TZZ)*(-dA)
              Finlet(3) = Finlet(3) + PRES*(-dA)
              Finlet(4) = Finlet(4) - SZZ*(-dA)
              Finlet(5) = Finlet(5) + WSQ*(-dA)
            ENDIF

            IF(KEFLAG>0) THEN
c...Outlet boundary.
              k=kce

c              write(20,*) xc(i),yc(j),p(i,j,k)

              !x-component
              DWDX = (w(i+1,j,k)-w(i-1,j,k))*aw(i)/2.
              DUDZ1 = (u(i  ,j,k+1)-u(i  ,j,k))*cw(k)
              DUDZ2 = (u(i-1,j,k+1)-u(i-1,j,k))*cw(k)
              DUDZ = 0.5*(DUDZ1+DUDZ2)
              SXZ = DUDZ+DWDX
              TXZ = RU1*SXZ
              U1 = 0.25*(U(i,j,k)+U(i,j,k+1)+U(i-1,j,k)+U(i-1,j,k+1))
              UW = U1*W(i,j,k)
              FX(2) = FX(2) + (UW - TXZ)*dA
              Foutlet(1) = Foutlet(1) + UW*dA
              Foutlet(2) = Foutlet(2) - TXZ*dA

              DWDZ    = (W(I,J,KCE+1)-W(I,J,KCE-1))*CW(KCE)/2.
              SZZ     = 2.*RU1*DWDZ
              PRES    = 0.5*(P(I,J,KCE+1)+P(I,J,KCE))
              TZZ     = PRES-SZZ
              WSQ     = W(I,J,KCE)**2.

              !z-component.
              FZ(2)  = FZ(2) + (WSQ + TZZ)*dA
              Foutlet(3) = Foutlet(3) + PRES*dA
              Foutlet(4) = Foutlet(4) -SZZ*dA
              Foutlet(5) = Foutlet(5) + WSQ*dA
            ENDIF

          ENDDO
        ENDDO
c
        i = ice
        DO K=KCB+1,KCE
          DO J=JY1,JY2
c...Top plane
            dA = RU(I)/BP(J)/CP(K)

            !x-component
            DUDX = (u(i+1,j,k)-u(i-1,j,k))*au(i)/2
            SXX = 2.*ru1*dudx
            PRES = 0.5*(p(i+1,j,k)+p(i,j,k))
            TXX = PRES - SXX
c          TXX = 0.0
            USQ = u(i,j,k)**2.
c          USQ = 0.0
            FX(3) = FX(3) + (USQ + TXX)*dA
            Ftop(1) = Ftop(1) + USQ*dA
            Ftop(2) = Ftop(2) + PRES*dA
            Ftop(3) = Ftop(3) - SXX*dA

            !z-component
            DUDZ  = (U(I  ,J,K+1) - U(I,J,K-1))*CU(K)/2.
            DWDX1 = (W(I+1,J,K  ) - W(I,J,K  ))*AU(I)
            DWDX2 = (W(I+1,J,K-1) - W(I,J,K-1))*AU(I)
            DWDX  = 0.5*(DWDX1+DWDX2)
            W1    = 0.25*(W(I,J,K) + W(I,J,K-1)
     &           + W(I+1,J,K) + W(I+1,J,K-1) )
            UW    = U(I,J,K)*W1
            
            SXZ   = RU1*(DUDZ+DWDX)

            FZ(3)    = FZ(3) + (UW - SXZ)*dA
            Ftop(4) = Ftop(4) + UW*dA
            Ftop(5) = Ftop(5) - SXZ*dA

          ENDDO
        ENDDO

        i=icb
        DO K=KCB+1,KCE
          DO J=JY1,JY2
            dA = 1./BP(J)/CP(K)

            !x-component
            DUDX = (u(i+1,j,k)-u(i-1,j,k))*au(i)/2
            SXX = 2.*ru1*dudx
            PRES = 0.5*(p(i+1,j,k)+p(i,j,k))
            TXX = PRES - SXX
c            TXX = 0.0
            USQ = u(i,j,k)**2.
c            USQ = 0.0
            FX(4) = FX(4) + (USQ + TXX)*(-dA)
            Fbot(1) = Fbot(1) + USQ*(-dA)
            Fbot(2) = Fbot(2) + PRES*(-dA)
            Fbot(3) = Fbot(3) - SXX*(-dA)

            DUDZ  = (U(I,J,K+1) - U(I,J,K-1))*CU(K)/2.
            DWDX1 = (W(I+1,J,K  ) - W(I,J,K  ))*AU(ICB)
            DWDX2 = (W(I+1,J,K-1) - W(I,J,K-1))*AU(ICB)
            DWDX  = 0.5*(DWDX1+DWDX2)
            W1    = 0.25*(W(I,J,K) + W(I,J,K-1)
     &           + W(I+1,J,K) + W(I+1,J,K-1) )
            UW    = U(I,J,K)*W1
            SXZ   = RU1*(DUDZ+DWDX)
            FZ(4)    = FZ(4) + (UW - SXZ)*(-dA)
            Fbot(4) = Fbot(4) + UW*(-dA)
            Fbot(5) = Fbot(5) - SXZ*(-dA)

          ENDDO
        ENDDO

        DO K=KCB+1,KCE
          DO I=ICB+1,ICE
            dA = 1./CP(K)/AP(I)
            
            J = 1
            if(flagu(i,j,k)<=0.OR.flagv(i,j,k)<=0) then
              DUDY1 = (u(i,j+1,k)-u(i,j,k))*bu(j)
              DUDY2 = (u(i-1,j+1,k)-u(i-1,j,k))*bu(j)
              DUDY = 0.5*(DUDY1+DUDY2)
              DVDX = (v(i+1,j,k)-v(i-1,j,k))*av(i)/2.
              SXY = RU1*(DUDY+DVDX)
              U1 = 0.25*(u(i,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i-1,j+1,k))
              UV = v(i,j,k)*U1
              FX(5) = FX(5) + (UV - SXY)*(-dA)
            endif


            if(flagv(i,j,k)<=0.OR.flagw(i,j,k)<=0) then
              DVDZ = (V(I,J,K+1)-V(I,J,K-1))*CU(K)/2.
              DWDY1 = (W(I,J+1,K  )-W(I,J,K  ))*BW(J)
              DWDY2 = (W(I,J+1,K-1)-W(I,J,K-1))*BW(J)
              DWDY  = 0.5*(DWDY1+DWDY2)
              SYZ = RU1*(DVDZ + DWDY)
              W1 = 0.25*( W(I,J  ,K) + W(I,J  ,K-1)
     &             +   W(I,J+1,K) + W(I,J+1,K-1))
              VW = V(I,J,K)*W1

              FZ(5) = FZ(5) + (VW - SYZ)*(-dA)
            endif

            J = JY2
            if(flagu(i,j,k)<=0.OR.flagv(i,j,k)<=0) then
              DUDY1 = (u(i,j+1,k)-u(i,j,k))*bu(j)
              DUDY2 = (u(i-1,j+1,k)-u(i-1,j,k))*bu(j)
              DUDY = 0.5*(DUDY1+DUDY2)
              DVDX = (v(i+1,j,k)-v(i-1,j,k))*av(i)/2.
              SXY = RU1*(DUDY+DVDX)
              U1 = 0.25*(u(i,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i-1,j+1,k))
              UV = v(i,j,k)*U1
              FX(6) = FX(6) + (UV - SXY)*dA
            endif

            if(flagv(i,j,k)<=0.OR.flagw(i,j,k)<=0) then
              DWDY1 = (W(I,J+1,K  )-W(I,J,K  ))*BW(J)
              DWDY1 = (W(I,J+1,K-1)-W(I,J,K-1))*BW(J)
              DWDY  = 0.5*(DWDY1+DWDY2)
              DVDZ = (V(I,J,K+1)-V(I,J,K-1))*CU(K)/2.
              SYZ = RU1*(DVDZ + DWDY)
              W1 = 0.25*( W(I,J  ,K) + W(I,J  ,K-1)
     &             +   W(I,J+1,K) + W(I,J+1,K-1))
              VW = V(I,J,K)*W1

              FZ(6) = FZ(6) + (VW - SYZ)*dA
            endif

          ENDDO
        ENDDO
      endif


c      close(10)
c      close(20)

      Finlet(6) = SUM(Finlet(1:2))
      Finlet(7) = SUM(Finlet(3:5))
      Foutlet(6) = SUM(Foutlet(1:2))
      Foutlet(7) = SUM(Foutlet(3:5))
      Ftop(6) = SUM(Ftop(1:3))
      Ftop(7) = SUM(Ftop(4:5))
      Fbot(6) = SUM(Fbot(1:3))
      Fbot(7) = SUM(Fbot(4:5))


      CALL MPI_REDUCE(FX,FXG,6,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(FY,FYG,6,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(FZ,FZG,6,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(Finlet,FGinlet,8,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(Foutlet,FGoutlet,8,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(Ftop,FGtop,8,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(Fbot,FGbot,8,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)


      IF(MYRANK.EQ.0) THEN

        IF(icycle==0 .AND. icalf<2) THEN
          OPEN(UNIT=91,FILE=trim(filestem)//'-z.plt',STATUS='REPLACE')
          WRITE(91,'(3A)')  'VARIABLES = "T", "Fz-Inlet", "Fz-Outlet"'
     &         ,', "Fz-Top", "Fz-Bottom", "Fz-Left", "Fz-Right"'
     &         ,', "Fz-Unst.", "Fz", ZONE DATAPACKING=POINT'
          CLOSE(91)

          OPEN(UNIT=91,FILE=trim(filestem)//'-x.plt',STATUS='REPLACE')
          WRITE(91,'(3A)')  'VARIABLES = "T", "Fx-Inlet", "Fx-Outlet"'
     &         ,', "Fx-Top", "Fx-Bottom", "Fx-Left", "Fx-Right"'
     &         ,', "Fx-Unst.", "Fx", ZONE DATAPACKING=POINT'
          CLOSE(91)

          OPEN(UNIT=91,FILE=trim(filestem)//'-y.plt',STATUS='REPLACE')
          WRITE(91,'(3A)')  'VARIABLES = "T", "Fy-Inlet", "Fy-Outlet"'
     &         ,', "Fy-Top", "Fy-Bottom", "Fy-Left", "Fy-Right"'
     &         ,', "Fy-Unst.", "Fy", ZONE DATAPACKING=POINT'
          CLOSE(91)

          OPEN(UNIT=91,FILE=trim(filestem)//'-inlet.plt'
     &         ,STATUS='REPLACE')
          WRITE(91,'(3A)')  'VARIABLES = "T", "Fx-UW,VW", "Fx-SXZ,SYZ"'
     &         ,', "Fy-UW,VW", "Fy-SXZ,SYZ", "Fz-Pres", "Fz-Szz",'
     &         ,' "Fz-WW", ZONE DATAPACKING=POINT'
          CLOSE(91)

          OPEN(UNIT=91,FILE=trim(filestem)//'-outet.plt'
     &         ,STATUS='REPLACE')
          WRITE(91,'(3A)')  'VARIABLES = "T", "Fx-UW", "Fx-SXZ"'
     &         ,', "Fz-Pres", "Fz-SZZ", "Fz-WW", "Fx-Outlet",' 
     &         ,'"Fz-Outlet", ZONE DATAPACKING=POINT'
          CLOSE(91)

          OPEN(UNIT=91,FILE=trim(filestem)//'-top.plt'
     &         ,STATUS='REPLACE')
          WRITE(91,'(3A)')  'VARIABLES = "T", "Fx-UU", "Fx-Pres"'
     &         ,', "Fx-SXX", "Fz-UW", "Fz-SXZ", "Fx-Top",'
     &         ,' "Fz-Top", ZONE DATAPACKING=POINT'
          CLOSE(91)

          OPEN(UNIT=91,FILE=trim(filestem)//'-bottom.plt'
     &         ,STATUS='REPLACE')
          WRITE(91,'(3A)')  'VARIABLES = "T", "Fx-UU", "Fx-Pres"'
     &         ,', "Fx-SXX", "Fz-UW", "Fz-SXZ", "Fx-Bottom",'
     &         ,' "Fz-Bottom", ZONE DATAPACKING=POINT'
          CLOSE(91)
        ELSE


          open(unit=91,file=trim(filestem)//'-x.plt',position='append'
     &          ,form='formatted')
          write(91,'(10(E15.8,1x))') TIME,FXG(1),FXG(2),FXG(3)
     &         ,FXG(4),FXG(5),FXG(6),FBX,-(SUM(FXG)+FBX)
          close(91)

          open(unit=91,file=trim(filestem)//'-y.plt',position='append'
     &          ,form='formatted')
          write(91,'(10(E15.8,1x))') TIME,FYG(1),FYG(2),FYG(3)
     &         ,FYG(4),FYG(5),FYG(6),FBY,-(SUM(FYG)+FBY)
          close(91)

          open(unit=91,file=trim(filestem)//'-z.plt',position='append'
     &          ,form='formatted')
          write(91,'(10(E15.8,1x))') TIME,FZG(1),FZG(2),FZG(3)
     &         ,FZG(4),FZG(5),FZG(6),FBZ,-(SUM(FZG)+FBZ)
          close(91)

          open(unit=91,file=trim(filestem)//'-inlet.plt'
     &         ,position='append',form='formatted')
          write(91,'(10(E15.8,1x))') TIME,FGinlet(1),FGinlet(2)
     &         ,FGinlet(3),FGinlet(4),FGinlet(5),FGinlet(6),FGinlet(7)
          close(91)

          open(unit=91,file=trim(filestem)//'-outet.plt'
     &         ,position='append',form='formatted')
          write(91,'(10(E15.8,1x))') TIME,FGoutlet(1),FGoutlet(2)
     &         ,FGoutlet(3),FGoutlet(4),FGoutlet(5),FGoutlet(6),FGoutlet(7)
          close(91)

          open(unit=91,file=trim(filestem)//'-top.plt'
     &         ,position='append',form='formatted')
          write(91,'(10(E15.8,1x))') TIME,FGtop(1),FGtop(2)
     &         ,FGtop(3),FGtop(4),FGtop(5),FGtop(6),FGtop(7)
          close(91)

          open(unit=91,file=trim(filestem)//'-bottom.plt'
     &         ,position='append',form='formatted')
          write(91,'(10(E15.8,1x))') TIME,FGbot(1),FGbot(2)
     &         ,FGbot(3),FGbot(4),FGbot(5),FGbot(6),FGbot(7)
          close(91)


        ENDIF

c        WRITE(9,'(A,E15.8)') 'Force in Z at the inlet:',FZG(1)
c        WRITE(9,'(A,E15.8)') 'Force in Z at the outlet:',FZG(2)
c        WRITE(9,'(A,E15.8)') 'Force in Z at the top:',FZG(3)
c        if(icyl==0) then
c          WRITE(9,'(A,E15.8)') 'Force in Z at the bottom:',FZG(4)
c          WRITE(9,'(A,E15.8)') 'Force in Z at the left:',FZG(5)
c          WRITE(9,'(A,E15.8)') 'Force in Z at the right:',FZG(6)
c        endif
c        WRITE(9,'(A,E15.8)') 'Total Drag:',SUM(FZG)
      ENDIF
c
      RETURN

      END
c------------------------------------------------------------------------


C---- subroutine press_body--------------------N. Beratlis-26 Dec. 2008-
C
C     PURPOSE: Calculate pressure on surface on immersed body by inter-
C     polating the value of pressure on the center of triangles using
C     the nodes of the cell that contains the vertex center
C     WARNING: Works assuming additional pressure nodes have been corrected.
C
C-----------------------------------------------------------------------
      subroutine press_body(p,xc,yc,zc,nx,ny,nz,pbd,mrkpb,vertexc,nfacet,ilb,ile)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
      
      integer nx,ny,nz,nfacet,ilb,ile

      integer mrkpb(nfacet)
      real    xc(nx),yc(ny),zc(nz)
      real    p(nx,ny,nz)
      real    vertexc(3,nfacet),pbd(nfacet)

      integer i,j,k,ii,nb,ib
      real    xp,yp,zp

      integer, dimension(:), allocatable :: mrkpb1
      real, dimension(:), allocatable :: pbd1

      real    interp3d

      nb = ile-ilb+1

      ALLOCATE(pbd1(nb),mrkpb1(nb))
      pbd1 = 0.0
      mrkpb1 = 0

      ib = 0
      
      if(icyl==0) then

        do ii=ilb,ile
          ib=ib+1
          xp = vertexc(1,ii)
          yp = vertexc(2,ii)
          zp = vertexc(3,ii)
          if(zp>=zc(1).AND.zp<zc(nz-1).AND.yp>=yc(1).AND.yp<yc(ny)) THEN
            pbd1(ib) = interp3d(xp,yp,zp,xc,yc,zc,p,nx,ny,nz,icyl)
            mrkpb1(ib) = 1
          endif
        enddo

      else

        do ii=ilb,ile
          ib=ib+1
          xp = vertexc(1,ii)
          yp = vertexc(2,ii)
          zp = vertexc(3,ii)
          if(zp>=zc(1).AND.zp<zc(nz-1)) then
            pbd1(ib) = interp3d(xp,yp,zp,xc,yc,zc,p,nx,ny,nz,icyl)
            mrkpb1(ib) = 1
c            write(6,*) ib,'pbd1=',pbd1(ib)
          endif
        enddo

      endif

      CALL MPI_ALLREDUCE(pbd1,pbd(ilb),nb,MTYPE,MPI_SUM,MPI_COMM_EDDY,IERR)
      CALL MPI_ALLREDUCE(mrkpb1,mrkpb(ilb),nb,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)

      DEALLOCATE(pbd1,mrkpb1)

      return

      end
C-----------------------------------------------------------------------



C---- subroutine shear_body  ------------------N. Beratlis-16 Mar. 2010-
C
C     PURPOSE: Calculate shear strain on surface on immersed body
C
C-----------------------------------------------------------------------
      subroutine shear_body(uo,vo,wo,nx,ny,nz,xu,yv,zw,xc,yc,zc
     &     ,xu_car,yu_car,xv_car,yv_car,xc_car,yc_car,vertexc,unvect
     &     ,nfacet,dudxb,dudyb,dudzb,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb
     &     ,cf_aux,mrks,ibd,nbd,ilb,ile,tlevel,io)

c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      INTEGER nx,ny,nz,ibd,nbd,nfacet,io,ilb,ile
      REAL    tlevel
      INTEGER mrks(9,nfacet)
      REAL    xu(nx),yv(ny),zw(nz),xc(nx),yc(ny),zc(nz)
      REAL    xu_car(nx,ny),yu_car(nx,ny),xv_car(nx,ny),yv_car(nx,ny)
     &       ,xc_car(nx,ny),yc_car(nx,ny)
      REAL    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
      REAL    vertexc(3,nfacet),unvect(3,nfacet)
      REAL    dudxb(nfacet),dudyb(nfacet),dudzb(nfacet)
     &       ,dvdxb(nfacet),dvdyb(nfacet),dvdzb(nfacet)
     &       ,dwdxb(nfacet),dwdyb(nfacet),dwdzb(nfacet)
      REAL    cf_aux(nfacet,6)

c
c.... Local arrays
      INTEGER i,j,k,ii,ifacet,iext,jext,kext,nzg,iflag,n1,dir,ksign,n
      REAL    sxx,syy,szz,sxy,syx,szy,syz,sxz,szx,up
      REAL    dudn,dvdn,dwdn,r,theta
      REAL    dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dn,dx,dy,dz
      REAL    xp,yp,zp,xint,yint,zint,uint,vint,wint,ub,vb,wb,psi,a1,a2,a3,unorm,amin
      REAL    uF1,uF2,dF1,dF2,c1,c2,c3,ub1,a,b,c,unF,utbd,unbd,unF1,dutdn,dundn
      REAL    un1(3),ut1(3),ut1nrm(3),ubdt(3),uext(2),ubdn(3),ubdy(3)
      REAL    amtrx(ivelrcn,ivelrcn),rmtrx(ivelrcn),utF(ivelrcn)
      INTEGER icount(6),icountg(6),indx(ivelrcn)
      REAL    q(3),rvec(3)
      INTEGER nl_req,nr_req,ilu,iru,ilw,irw,nmax
      INTEGER, DIMENSION(:), ALLOCATABLE :: uordl,uordr,wordl,wordr
      REAL, DIMENSION(:,:), ALLOCATABLE :: ucordl,ucordr,wcordl,wcordr
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: uvec,xyzu
      REAL, DIMENSION(:), ALLOCATABLE :: ul,ur,wl,wr,dcell
      INTEGER, SAVE :: icall=0

      LOGICAL condx,condy,condz
c
c.... Functions
      INTEGER body2fluid,body2grid2
      REAL    vecmag,ubd,vbd,wbd,interp_cellface,anglerad,interp3d
     $     ,dotproduct,point_grid_size_ray

      amin = 0.75
      icount=0

      n1=4

      mrks(:,ilb:ile) = 0      !!!!!! instead of mrks = 0
      dudxb(ilb:ile) = 0.0     !!!!!! instead of dudxb = 0.0
      dudyb(ilb:ile) = 0.0     !!!!!!
      dudzb(ilb:ile) = 0.0     !!!!!!
      dvdxb(ilb:ile) = 0.0     !!!!!!
      dvdyb(ilb:ile) = 0.0     !!!!!!
      dvdzb(ilb:ile) = 0.0     !!!!!!
      dwdxb(ilb:ile) = 0.0     !!!!!!
      dwdyb(ilb:ile) = 0.0     !!!!!!
      dwdzb(ilb:ile) = 0.0     !!!!!!

      cf_aux = 0.0


      q = 0.0
      rvec = 0.0
      xint = 0.0
      yint = 0.0
      zint = 0.0

      IF(ivelrcn==0) THEN

        DO ifacet=ilb,ile

          dudn=0
          dvdn=0
          dwdn=0

          xp = vertexc(1,ifacet)
          yp = vertexc(2,ifacet)
          zp = vertexc(3,ifacet)

          q(1) = xp
          q(2) = yp
          q(3) = zp   

          rvec(1) = unvect(1,ifacet)
          rvec(2) = unvect(2,ifacet)
          rvec(3) = unvect(3,ifacet)
          ksign = 0
          if(rvec(3)>0) ksign=1

          if(icyl==0) then
            condx = yp>=ymin.AND.yp<=ymax .AND. zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
            condz = yp>=ymin.AND.yp<=ymax .AND. zp>=zw(2-ksign) .AND. zp<zw(nz-ksign)
          else
            condx = zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
            condz = zp>=zw(2-ksign) .AND. zp<zw(nz-ksign)
          endif

          if(icyl==0) then
            theta = yp
            r = xp
          else
            theta = anglerad(xp,yp)
            r = sqrt(xp**2. + yp**2.) 
          endif

          if(condx) then

            icount(1) = icount(1)+1
            !Calculate derivatives of u
c in BODY2GRID2 the interpolated value for UO is UINT
            mrks(1,ifacet) = body2grid2(q,rvec,xu_car,yu_car,xu,yc,zc
     &           ,uo,nx,ny,nz,xint,yint,zint,uint,n1,nbd,icyl,ifacet,amin)
            if(mrks(1,ifacet)<=n1) then
c a valid intersection has been found in BODY2GRID2
              icount(2) = icount(2)+1
              if(io==1) then
c moving boundary
                ub = ubd(xp,yp,zp,tlevel,ibd)
              else
c stationary boundary
                ub = 0.0
              endif

              dn = sqrt( (xp-xint)**2. + (yp-yint)**2. + (zp-zint)**2.)
c the normal derivative for U is evaluated
              dudn = (uint-ub)/dn
            endif

          icount(3) = icount(3)+1     !!!!!!!
            !Calculate derivatives of v
            mrks(2,ifacet) = body2grid2(q,rvec,xv_car,yv_car,xc,yv,zc
     &           ,vo,nx,ny,nz,xint,yint,zint,vint,n1,nbd,icyl,ifacet,amin)           
            if(mrks(2,ifacet)<=n1) then
              icount(4) = icount(4)+1
              if(io==1) then
                vb = vbd(xp,yp,zp,tlevel,ibd)
              else
                vb = 0.0
              endif

              dn = sqrt( (xp-xint)**2. + (yp-yint)**2. + (zp-zint)**2.)
              dvdn = (vint-vb)/dn
            endif

          endif

          if(condz) then

            icount(5) = icount(5)+1
            !Calculate derivatives of w
            mrks(3,ifacet) = body2grid2(q,rvec,xc_car,yc_car,xc,yc,zw
     &           ,wo,nx,ny,nz,xint,yint,zint,wint,n1,nbd,icyl,ifacet,amin)
            if(mrks(3,ifacet)<=n1) then
              mrks(7:9,ifacet)=1
              icount(6) = icount(6)+1
              if(io==1) then
                wb = wbd(xp,yp,zp,tlevel,ibd)
              else
                wb = 0.0
              endif

              dn = sqrt( (xp-xint)**2. + (yp-yint)**2. + (zp-zint)**2.)
              dwdn = (wint-wb)/dn
            endif

          endif

          unorm = vecmag(unvect(:,ifacet),3)
c components of the unit normal vector
          a1 = unvect(1,ifacet)/unorm
          a2 = unvect(2,ifacet)/unorm
          a3 = unvect(3,ifacet)/unorm

          psi = anglerad(unvect(1,ifacet),unvect(2,ifacet))

c the components of the normal derivatives along X,Y,Z are evaluated
          if(icyl==0) then

            if(mrks(1,ifacet)>0) mrks(1:3,ifacet)=1
            dudxb(ifacet) = dudn*a1
            dudyb(ifacet) = dudn*a2
            dudzb(ifacet) = dudn*a3

            if(mrks(2,ifacet)>0) mrks(4:6,ifacet)=1
            dvdxb(ifacet) = dvdn*a1
            dvdyb(ifacet) = dvdn*a2
            dvdzb(ifacet) = dvdn*a3

          else
  
            if(mrks(1,ifacet)>0 .AND. mrks(2,ifacet)>0) mrks(1:6,ifacet)=1

            ub = dudn
            vb = dvdn

            dudn = ub*cos(psi)-vb*sin(psi)
            dudxb(ifacet) = dudn*a1
            dudyb(ifacet) = dudn*a2
            dudzb(ifacet) = dudn*a3

            dvdn = ub*sin(psi)+vb*cos(psi)
            dvdxb(ifacet) = dvdn*a1
            dvdyb(ifacet) = dvdn*a2
            dvdzb(ifacet) = dvdn*a3

          endif
          
          dwdxb(ifacet) = dwdn*a1
          dwdyb(ifacet) = dwdn*a2
          dwdzb(ifacet) = dwdn*a3

        ENDDO
      ELSE  ! if ivelrcn==1 i.e the extension of the ibm is 1 then

        nmax = ile-ilb+1
        allocate(ucordl(3,nmax),ucordr(3,nmax))
        allocate(uordl(nmax),uordr(nmax))
        allocate(wcordl(3,nmax),wcordr(3,nmax))
        allocate(wordl(nmax),wordr(nmax))
        ucordl = 0.0      
        ucordr = 0.0
        wcordl = 0.0
        wcordr = 0.0
        uordl = 0
        uordr = 0
        wordl = 0
        wordr = 0

        n = ile-ilb+1
        allocate(uvec(3,n,2),xyzu(3,n,2))
        allocate(dcell(n))

        uvec = 0.0
        xyzu = 0.0
        dcell = 0.0

        uext(1) = uext1 
        uext(2) = uext2

        DO iext=1,2 
! iext=1 the velocities are interpolated only in one point over the body. 
! iext=2 the velocities are also interpolated further from the surface. 
! This is though to be more accurate although most of the time ext2 data is not used.    

          ilu = 0
          iru = 0
          ilw = 0
          irw = 0

          DO ifacet=ilb,ile

            ii = ifacet-ilb+1

            dudn=0
            dvdn=0
            dwdn=0

            xp = vertexc(1,ifacet)  !coordianates of the triangle defining the body.
            yp = vertexc(2,ifacet)
            zp = vertexc(3,ifacet)

            ksign = 0
            if(unvect(3,ifacet)>0) ksign=1 
 
            if(icyl==0) then    !check that the body is well defined and inside the domain.   
              condx = yp>=ymin.AND.yp<=ymax .AND. zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
            else
              condx = zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
            endif


            if(condx) then

              if(iext==1) then
               dcell(ii) = point_grid_size_ray(xp,yp,zp,unvect(:,ifacet),xc
     $             ,yc,zc,xc_car,yc_car,nx,ny,nz,icyl)
! Finds the distance between the two intersections of a normal ray (positive and negative senses) emerging from a cell center. 
!This cell center is the closest grid cell center to the center of the triangle defining the surface. 
!The direction is given by the normal vector to that triangle. 
!Any other distance could have been used,this is used only for convenience. 
 
              endif
              dn = dcell(ii)

              q(1) = xp
              q(2) = yp
              q(3) = zp
              rvec = unvect(:,ifacet)
              dn = uext(iext)*dcell(ii)
! A point (xyzu) along the normal to the surface of the boundary a distance dn.
              xyzu(:,ii,iext) = q + rvec*dn

! Check that the point is in the domain included in this processor. 

              if(xyzu(3,ii,iext)<zc(1)) then
! The point is on the left of the local domain.
                ilu = ilu + 1
                ucordl(:,ilu) = xyzu(:,ii,iext)
                uordl(ilu) = ii
              elseif(xyzu(3,ii,iext)>=zc(nz)) then
! The point is on the right of the local domain.
                iru = iru + 1
                ucordr(:,iru) = xyzu(:,ii,iext)
                uordr(iru) = ii
              else
! If the point is inside the local domain (the domain in this processor): the velocity components are evaluated by interpolations. Uvec is interpolated velocity over the surface of the body at a distance dcell.

                uvec(1,ii,iext) = interp3d(xyzu(1,ii,iext),xyzu(2,ii,iext),xyzu(3,ii,iext),xu,yc,zc,uo,nx,ny,nz,icyl)
                uvec(2,ii,iext) = interp3d(xyzu(1,ii,iext),xyzu(2,ii,iext),xyzu(3,ii,iext),xc,yv,zc,vo,nx,ny,nz,icyl)
              endif

              if(xyzu(3,ii,iext)<zw(1)) then
                ilw = ilw + 1
                wcordl(:,ilw) = xyzu(:,ii,iext)
                wordl(ilw) = ii
              elseif(xyzu(3,ii,iext)>=zw(nz)) then
                irw = irw + 1
                wcordr(:,irw) = xyzu(:,ii,iext)
                wordr(irw) = ii
              else
                uvec(3,ii,iext) = interp3d(xyzu(1,ii,iext),xyzu(2,ii,iext),xyzu(3,ii,iext),xc,yc,zw,wo,nx,ny,nz,icyl)
              endif

            endif
            
          ENDDO

          IF(mysize>1) THEN

! Interpolation and communication of variables for the left and the right processors.

            nl_req = ilu
            nr_req = iru

            allocate(ul(nl_req),ur(nr_req))
            ul = 0.0
            ur = 0.0
            !Interp. and comm. u comp.
            call  mpi_var_interp(ucordl,ul,nl_req,ucordr,ur,nr_req,uo,xu,yc,zc,nx,ny,nz)

            do i=1,nl_req
              uvec(1,uordl(i),iext) = ul(i)
            enddo

            do i=1,nr_req
              uvec(1,uordr(i),iext) = ur(i)
            enddo

            ul = 0.0
            ur = 0.0
            !Interp. and comm. v comp.
            call  mpi_var_interp(ucordl,ul,nl_req,ucordr,ur,nr_req,vo,xc,yv,zc,nx,ny,nz)

            do i=1,nl_req
              uvec(2,uordl(i),iext) = ul(i)
            enddo

            do i=1,nr_req
              uvec(2,uordr(i),iext) = ur(i)
            enddo
            deallocate(ul,ur)

            nl_req = ilw
            nr_req = irw

            allocate(wl(nl_req),wr(nr_req))
            wl = 0.0
            wr = 0.0

            !Interp. and comm. w comp.
            call  mpi_var_interp(wcordl,wl,nl_req,wcordr,wr,nr_req,wo,xc,yc,zw,nx,ny,nz)

            do i=1,nl_req
              uvec(3,wordl(i),iext) = wl(i)
            enddo

            do i=1,nr_req
              uvec(3,wordr(i),iext) = wr(i)
            enddo

            deallocate(wl,wr)

          ENDIF

        ENDDO

        deallocate(ucordl,ucordr,wcordl,wcordr)
        deallocate(uordl,uordr,wordl,wordr)

        DO ifacet=ilb,ile

          ii = ifacet-ilb+1

          ksign = 0
          if(unvect(3,ifacet)>0) ksign=1

          xp = vertexc(1,ifacet) !coordinates of the center of the triangle.
          yp = vertexc(2,ifacet)
          zp = vertexc(3,ifacet)

          if(icyl==0) then  !check that the body is inside the domain.
            condx = yp>=ymin.AND.yp<=ymax .AND. zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
          else
            condx = zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
          endif
          
          if(condx) then !convert velocity over the surface uvec to cartesian. 
            call convert_velcyl2car(uvec(1,ii,1),uvec(2,ii,1),uvec(1,ii,1),uvec(2,ii,1),xyzu(1,ii,1),xyzu(2,ii,1))
            call convert_velcyl2car(uvec(1,ii,2),uvec(2,ii,2),uvec(1,ii,2),uvec(2,ii,2),xyzu(1,ii,2),xyzu(2,ii,2))
          endif

        ENDDO


        DO ifacet=ilb,ile

          ii = ifacet-ilb+1

          ksign = 0
          if(unvect(3,ifacet)>0) ksign=1

          xp = vertexc(1,ifacet)
          yp = vertexc(2,ifacet)
          zp = vertexc(3,ifacet)

          if(icyl==0) then
            condx = yp>=ymin.AND.yp<=ymax .AND. zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
          else
            condx = zp>=zc(2-ksign) .AND. zp<zc(nz-ksign)
          endif
          
          if(condx) then

! (xp,yp,zp) is inside the local domain

            icount = icount + 1

            unF1 = dotproduct(uvec(:,ii,1),unvect(:,ifacet),3)
! normal vector associated with uvec(:,ii,1)
            un1 = unF1*unvect(:,ifacet)
         
! tangential vector associated with uvec(:,ii,1)
            ut1 = uvec(:,ii,1) - un1
! corresponding unitary tangential vector
            ut1nrm = ut1/vecmag(ut1,3)
          
            uF1 = vecmag(ut1,3)
            uF2 = dotproduct(uvec(:,ii,2),ut1nrm,3)

            ubdy = 0.0
            ubdt = 0.0
            ubdn = 0.0
            ub1 = 0.0
            if(io==1) then
! components of the boundary velocity
              ubdy(1) = ubd(vertexc(1,ifacet),vertexc(2,ifacet),vertexc(3,ifacet),tlevel,ibd)
              ubdy(2) = vbd(vertexc(1,ifacet),vertexc(2,ifacet),vertexc(3,ifacet),tlevel,ibd)
              ubdy(3) = wbd(vertexc(1,ifacet),vertexc(2,ifacet),vertexc(3,ifacet),tlevel,ibd)
! conversion from cylindrical to Cartesian velocity components
              if(icyl==1) call convert_velcyl2car(ubdy(1),ubdy(2),ubdy(1),ubdy(2),vertexc(1,ifacet),vertexc(2,ifacet))
            endif
            utbd = dotproduct(ubdy,ut1nrm,3)
! tangential component of the boundary velocity
            ubdt = utbd*ut1nrm
            unbd = dotproduct(ubdy,unvect(:,ifacet),3)
! normal component of the boundary velocity
            ubdn = unbd*unvect(:,ifacet)            

! distance over the body for both extensions
            dF1 = uext(1)*dcell(ii)
            dF2 = uext(2)*dcell(ii)

            !Tangential velocity component

            !One velocity and velocity gradient
c            dn = (uext(2)-uext(1))/(dF2-dF1)
c            dudn = (uF2-uF1)/dn
c            c1 = ( 1./dF*dudn - uF1/(dF**2.) )
c            c2 = (2.*uF1/dF1 - dudn)
c            c3 = ub1

            amtrx(1,1) = dF1
            amtrx(1,2) = dF1**2.0
            amtrx(2,1) = 1.0
            amtrx(2,2) = 2.0*dF1


            call ludcmp(amtrx,ivelrcn,ivelrcn,indx,b)

! Difference between the fluid and the boundary tangential velocity components
            rmtrx(1) = uF1-utbd

! Normal derivative of the tangential velocity between the fluid points 1 and 2
            rmtrx(2) = (uF2-uF1)/(dF2-dF1)
            call lubksb(amtrx,ivelrcn,ivelrcn,indx,rmtrx)
            c2 = rmtrx(1)

! The LU decomposition and solver are not being used since ivelrcn ==1, dutdn=c2=uF1-utbd)/dF1. 

        
! Tangential velocity derivative along the normal direction
            dutdn = c2

! Normal velocity derivative along the normal direction     

            dundn = 0.0 !Quadratic           
! dundn = (unF1-unbd)/dF1 Continuity shows this has to be zero. 
! The value is just numerical error that comes from non having infinite resolution. 
 
            a1 = unvect(1,ifacet)
            a2 = unvect(2,ifacet)
            a3 = unvect(3,ifacet)
            
! components of the normal derivative
            dudn = dutdn*ut1nrm(1) + dundn*unvect(1,ifacet)
            dvdn = dutdn*ut1nrm(2) + dundn*unvect(2,ifacet)
            dwdn = dutdn*ut1nrm(3) + dundn*unvect(3,ifacet)

! components of each velocity derivative along the x,y,z directions
            dudxb(ifacet) = dudn*a1
            dudyb(ifacet) = dudn*a2
            dudzb(ifacet) = dudn*a3

            dvdxb(ifacet) = dvdn*a1
            dvdyb(ifacet) = dvdn*a2
            dvdzb(ifacet) = dvdn*a3

            dwdxb(ifacet) = dwdn*a1
            dwdyb(ifacet) = dwdn*a2
            dwdzb(ifacet) = dwdn*a3
     

            cf_aux(ifacet,1)  = ut1(1)
            cf_aux(ifacet,2)  = ut1(2)
            cf_aux(ifacet,3)  = ut1(3)
            cf_aux(ifacet,4)  = dF1
            cf_aux(ifacet,5)  = dF2
            cf_aux(ifacet,6)  = uF1/dF1


! Turn on to compute C_tau using Aniskesh Matlab code;

!            dvdxb(ifacet) = ut1(1)
!            dvdyb(ifacet) = ut1(2) 
! This sign is needed to compute the friction coefficient and cannot be recovered from the tau_ij components. It has to be exported. 
!            dvdzb(ifacet) = ut1(3)
!            dwdxb(ifacet) = dF1 
!            dwdyb(ifacet) = dF2
!            dwdzb(ifacet) = uF1/dF1
            
              mrks(:,ifacet) = 1


          endif

        ENDDO

        deallocate(uvec,dcell,xyzu)

      ENDIF

      IF(mysize>1) THEN
       CALL MPI_REDUCE(icount,icountg,6,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      ELSE
        icountg = icount
      ENDIF

      IF(myrank==0 .AND. iolvl>0) then

        if(ibm>1 .OR. icall==0) then
!          if(icount(1)/=icountg(2) .OR. icount(3)/=icountg(4) .OR. icount(5)/=icountg(6)) then    !!!!!!
            open(unit=16,file='stats_imb.dat',position='append'
     &            ,form='formatted')
            write(16,*) 'Number of triangles=',mb(ibd)
            write(16,*) 'U-velocity:'
            write(16,*) 'Number of tri. within domain',icountg(1)    !!!!!! instead of icount(1)
            write(16,*) 'Number of dudn points extrap.',icountg(2)
            write(16,*) 'V-velocity:'
            write(16,*) 'Number of tri. within domain',icountg(3)    !!!!!! instead of icount(3)
            write(16,*) 'Number of dvdn points extrap.',icountg(4)
            write(16,*) 'W-velocity:'
            write(16,*) 'Number of tri. within domain',icountg(5)    !!!!!! instead of icount(5)
            write(16,*) 'Number of dwdn points extrap.',icountg(6)
            close(16)
!          endif    !!!!!!
        endif
      ENDIF

      icall = icall+1

      RETURN

      END
C------------------------------------------------------------------------



C---- subroutine press_node2center-----------N. Beratlis-29 Jun. 2009---
C
      subroutine press_node2center(pnd,node,elem,nnodes,pc,mrkp,vertexc,nfacet,ilb,ile)
C
C     PURPOSE: Interpolate pressure at center of triangle from values at 
C     nodes of triangle
C
C-----------------------------------------------------------------------
      implicit none
c
c.... Input/Output Arrays
      integer nnodes,nfacet,ilb,ile
      integer elem(3,nfacet),mrkp(nfacet)
      real    pnd(nnodes),pc(nfacet)
      real    node(3,nnodes),vertexc(3,nfacet)
c
c.... Local Arrays      
      integer i
      real    xc,yc,zc
      real    vt(3,3),pt(3)
c
c.... Functions
      real    triangle_interp

      mrkp = 1

      do i=ilb,ile

        xc = vertexc(1,i)
        yc = vertexc(2,i)
        zc = vertexc(3,i)
        
        vt(1,1) = node(1,elem(1,i))
        vt(1,2) = node(2,elem(1,i))
        vt(1,3) = node(3,elem(1,i))
        vt(2,1) = node(1,elem(2,i))
        vt(2,2) = node(2,elem(2,i))
        vt(2,3) = node(3,elem(2,i))
        vt(3,1) = node(1,elem(3,i))
        vt(3,2) = node(2,elem(3,i))
        vt(3,3) = node(3,elem(3,i))

        pt(1) = pnd(elem(1,i))
        pt(2) = pnd(elem(2,i))
        pt(3) = pnd(elem(3,i))

        pc(i) = triangle_interp(xc,yc,zc,vt,pt,i)

      enddo

      return

      end
C-----------------------------------------------------------------------

c-----subroutine calcforce ----------------------N. Beratlis-19 Jan. 2009-
c                                                S. Schroeder-5 Jan. 2012
c                                                A. Posa ------ Sep. 2013
c     PURPOSE: Calculates the forces + moments on the immersed boundary using
c     direct surface integration of the pressure and shear stress on
c     the immersed object.     
c
c------------------------------------------------------------------------
      subroutine calcforce(unvect,vertex,vertexc,areaf,pb,mrkp,dudxb,dudyb,dudzb
     &     ,dvdxb,dvdyb,dvdzb,dwdxb,dwdyb,dwdzb,mrks,nfacet,ibd,tlevel,icycle)
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
      INCLUDE 'immersed.h'
c
c.... Input/Output Arrays
      integer ibd,nfacet,icycle
      real    tlevel
      integer mrkp(nfacet),mrks(9,nfacet)    !!!!!! instead of mrks(nfacet,9)
      real    unvect(3,nfacet),areaf(nfacet),vertex(3,3,nfacet)
      real    vertexc(3,nfacet)
      real    dudxb(nfacet),dudyb(nfacet),dudzb(nfacet)
     &       ,dvdxb(nfacet),dvdyb(nfacet),dvdzb(nfacet)
     &       ,dwdxb(nfacet),dwdyb(nfacet),dwdzb(nfacet),pb(nfacet)
c
c.... Local arrays
      integer i,ifacet
      real    vxmin,vxmax,vymin,vymax,vyc,fr,ft,r,theta,areax,areay,areaz
      real    sxxb,syyb,szzb,sxyb,sxzb,syzb
      integer nf(2),nfg(2),icount(2),icountg(2)    !!!!!!
      real    miny(2),maxy(2),minx(2),maxx(2),minz(2),maxz(2)
      real    minyg(2),maxyg(2),minxg(2),maxxg(2),minzg(2),maxzg(2)
      real    area(2),areag(2),rvec(3)    !!!!!!
      real    fx(2),fy(2),fz(2),fxg(2),fyg(2),fzg(2)
      real    mx(2),my(2),mz(2),mxg(2),myg(2),mzg(2)
      logical condbound

      real    clocktemp,tclock
      real    clock(10)
      real    ymean
c
c.... Functions
      real    vecmag,anglerad,dotproduct
c
      clock = 0.

      clock(1) = tclock()

      clock(2) = tclock()
      area = 0.0
      fx = 0.0
      fy = 0.0
      fz = 0.0
      mx = 0.0
      my = 0.0
      mz = 0.0

      minx=10.
      miny=10.
      minz=10.
      maxx=-10.
      maxy=-10.
      maxz=-10.

      nf = 0
      icount = 0

      ymean = (ymin + ymax) * 0.5

      clock(2) = tclock()-clock(2)

      clock(8) = tclock()
      do i=lb(ibd)+1,lb(ibd)+mb(ibd)

        clocktemp = tclock()
        vxmin = minval(vertex(1,:,i))
        vxmax = maxval(vertex(1,:,i))
        vymin = minval(vertex(2,:,i))
        vymax = maxval(vertex(2,:,i))
        vyc = 1./3*sum(vertex(2,1:3,i))
        clock(3) = clock(3) + tclock() - clocktemp

c check for center points inside the domain
        if(icyl==0) then
          condbound = vyc>=ymin.AND.vyc<=ymax
        else
          condbound = .true.
        endif

        if(condbound) then
           
          icount(1) = icount(1)+1
          clocktemp = tclock()
          !Pressure contribution
c N.B. MRKP=1 when a valid probe has been found
          if(mrkp(i)==1) then

            area(1) = area(1)+areaf(i)
c pressure contributions
            fx(1) = fx(1) - pb(i)*areaf(i)*unvect(1,i)
            fy(1) = fy(1) - pb(i)*areaf(i)*unvect(2,i)
            fz(1) = fz(1) - pb(i)*areaf(i)*unvect(3,i)
            nf(1) = nf(1)+1
            
            !Moments
            mx(1) = mx(1) + pb(i)*areaf(i)*unvect(2,i)*vertexc(3,i)
            mx(1) = mx(1) - pb(i)*areaf(i)*unvect(3,i)*(vertexc(2,i)-ymean)
            my(1) = my(1) - pb(i)*areaf(i)*unvect(1,i)*vertexc(3,i)
            my(1) = my(1) + pb(i)*areaf(i)*unvect(3,i)*vertexc(1,i)
            mz(1) = mz(1) + pb(i)*areaf(i)*unvect(1,i)*(vertexc(2,i)-ymean)
            mz(1) = mz(1) - pb(i)*areaf(i)*unvect(2,i)*vertexc(1,i)

            !Find min and max coordinates of body
            if(vxmin<minx(1)) then
              minx(1)=vxmin
            endif
            if(vymin<miny(1)) then
              miny(1)=vymin
            endif

            if(vxmax>maxx(1)) then
              maxx(1)=vxmax
            endif
            if(vymax>maxy(1)) then
              maxy(1)=vymax
            endif

          endif
          clock(4) = clock(4) + tclock()-clocktemp

          clocktemp = tclock()
          !Shear contribution

          area(2) = area(2)+areaf(i)

c projection of the area on an XY surface
          areaz = areaf(i)*unvect(3,i)

!!!!!!          if(icyl<1) then    !!!!!!  1 instead of 2
            areax = areaf(i)*unvect(1,i)
            areay = areaf(i)*unvect(2,i)

c the estimate of the shear stress is based on
c tau_ij=nu*(du_i/dx_j+du_j/dx_i)
            sxxb = 2.*dudxb(i)
            syyb = 2.*dvdyb(i)
            szzb = 2.*dwdzb(i)
            sxyb = dudyb(i)+dvdxb(i)
            sxzb = dudzb(i)+dwdxb(i)
            syzb = dvdzb(i)+dwdyb(i)

c viscous forces
            fx(2) = fx(2) + ru1*sxxb*areax
            fx(2) = fx(2) + ru1*sxyb*areay
            fx(2) = fx(2) + ru1*sxzb*areaz

            fy(2) = fy(2) + ru1*syyb*areay
            fy(2) = fy(2) + ru1*sxyb*areax
            fy(2) = fy(2) + ru1*syzb*areaz
            
            !Moments
            mx(2) = mx(2) - vertexc(3,i)*ru1*(syyb*areay + sxyb*areax + syzb*areaz)
            my(2) = my(2) + vertexc(3,i)*ru1*(sxxb*areax + sxyb*areay + sxzb*areaz)
            mz(2) = mz(2) - (vertexc(2,i)-ymean)*ru1*(sxxb*areax + sxyb*areay + sxzb*areaz)
            mz(2) = mz(2) + vertexc(1,i)*ru1*(syyb*areay + sxyb*areax + syzb*areaz)
!!!!!!          else

!!!!!!            !compute area in cylindrical coordinates areax=a_r, areay=a_theta
!!!!!!            rvec(1) = vertexc(1,i)
!!!!!!            rvec(2) = vertexc(2,i)
!!!!!!            rvec(3) = 0.0
!!!!!!            rvec = rvec/vecmag(rvec,3)
!!!!!!            areax = areaf(i)*dotproduct(unvect(:,i),rvec,3)

!!!!!!            rvec(1) =-vertexc(2,i)
!!!!!!            rvec(2) = vertexc(1,i)
!!!!!!            rvec(3) = 0.0
!!!!!!            rvec = rvec/vecmag(rvec,3)
!!!!!!            areay = areaf(i)*dotproduct(unvect(:,i),rvec,3)
!!!!!!            if(sum(mrks(:,i))-mrks(3,i)==5) then
!!!!!!              fr = ru1*(sxxb*areax + sxyb*areay + sxzb*areaz )        
!!!!!!              ft = ru1*(syyb*areay + sxyb*areax + syzb*areaz )

!!!!!!              r = sqrt(vertexc(1,i)**2. + vertexc(2,i)**2.)
!!!!!!              if(r>0.0) then
!!!!!!                theta = anglerad(vertexc(1,i),vertexc(2,i))
!!!!!!                fx(2) = fx(2) + fr*cos(theta) - ft*sin(theta)
!!!!!!                fy(2) = fy(2) + fr*sin(theta) + ft*cos(theta)
!!!!!!              endif
!!!!!!            endif
!!!!!!          endif

          fz(2) = fz(2) + ru1*szzb*areaz
          fz(2) = fz(2) + ru1*sxzb*areax
          fz(2) = fz(2) + ru1*syzb*areay
          
          !Moments
          mx(2) = mx(2) + (vertexc(2,i)-ymean)*ru1*(szzb*areaz + sxzb*areax + syzb*areay)
          my(2) = my(2) - vertexc(1,i)*ru1*(szzb*areaz + sxzb*areax + syzb*areay)

          nf(2) = nf(2)+1

          !Find min and max coordinates of body
          if(vxmin<minx(2)) then
            minx(2)=vxmin
          endif
          if(vymin<miny(2)) then
            miny(2)=vymin
          endif

          if(vxmax>maxx(2)) then
            maxx(2)=vxmax
          endif
          if(vymax>maxy(2)) then
            maxy(2)=vymax
          endif

          clock(5) = clock(5) + tclock() - clocktemp

        endif

      enddo
      clock(8) = tclock() - clock(8)

      nfg=nf   !!!!!!
      areag=area    !!!!!!
      if(mysize>1) then
c Global forces and boundaries of the body are evaluated at the
c processor of rank 0
        CALL MPI_REDUCE(fx,fxg,2,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(fy,fyg,2,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(fz,fzg,2,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(minx,minxg,2,MTYPE,MPI_MIN,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(miny,minyg,2,MTYPE,MPI_MIN,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(minz,minzg,2,MTYPE,MPI_MIN,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(maxx,maxxg,2,MTYPE,MPI_MAX,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(maxy,maxyg,2,MTYPE,MPI_MAX,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(maxz,maxzg,2,MTYPE,MPI_MAX,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(mx,mxg,2,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(my,myg,2,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(mz,mzg,2,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
!!!!!!
        CALL MPI_REDUCE(nf(1),nfg(1),1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(area(1),areag(1),1,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
!!!!!!
      else
        fxg=fx
        fyg=fy
        fzg=fz
        mxg=mx
        myg=my
        mzg=mz
        minxg=minx
        minyg=miny
        minzg=minz
        maxxg=maxx
        maxyg=maxy
        maxzg=maxz
      endif

c Processor 0 writes the statistics
      if(myrank.eq.0 .AND. iolvl>0) then
        clock(6) = tclock()
        if(icycle==itcalf) then
c Only at the first time
          open(unit=16,file='stats_imb.dat',position='append'
     &          ,form='formatted')
          write(16,*) 'immersed-boundary no.=',ibd  !!!!!!
          write(16,*) 'no. of triangles within ymin-ymax',icount(1)
          write(16,*) 'area=',areag    !!!!!! areag instead of area
          write(16,*) 'nf=',nfg     !!!!!! nfg instead of nf
          write(16,*) 'X: min=',minxg,',max=',maxxg,',Lx=',maxxg-minxg
          write(16,*) 'Y: min=',minyg,',max=',maxyg,',Ly=',maxyg-minyg
          close(16)
        endif
      endif

      if(myrank.eq.0) then

        if(icycle==itcalf .AND. icalf<2) then
c Only at the first time
          open(unit=91, file='forces_xyz_imb'//index(ibd)//'.plt',status='replace')   !!!!!!
          write(91,'(3A)')  'VARIABLES = "T", "Px", "Sx", "Fx"' 
     &          ,', "Py", "Sy", "Fy", "Pz", "Sz", "Fz"'
     &          ,' ZONE DATAPACKING=POINT'
          close(91)

          open(unit=91, file='moments_xyz_imb'//index(ibd)//'.plt',status='replace')   !!!!!!
          write(91,'(3A)')  'VARIABLES = "T", "Mpx", "Msx", "Mx"'
     &          ,', "Mpy", "Msy", "My", "Mpz", "Msz", "Mz"'
     &          ,', ZONE DATAPACKING=POINT'
          close(91)
        endif

c Every time this subroutine is called
        open(unit=91,file='forces_xyz_imb'//index(ibd)//'.plt',position='append'
     &       ,form='formatted')    !!!!!!
        write(91,'(10(E15.8,1x))') tlevel,fxg(1),fxg(2),sum(fxg)
     &       ,fyg(1),fyg(2),sum(fyg),fzg(1),fzg(2),sum(fzg)
        close(91)

        open(unit=91,file='moments_xyz_imb'//index(ibd)//'.plt',position='append'
     &       ,form='formatted')    !!!!!!
        write(91,'(10(E15.8,1x))') tlevel,mxg(1),mxg(2),sum(mxg)
     &       ,myg(1),myg(2),sum(myg),mzg(1),mzg(2),sum(mzg)
        close(91)
        clock(6) = tclock() - clock(6)
      endif

      clock(1) = tclock() - clock(1)

c      open(unit=16, file='clock.dat', position='append'
c     &   , form='formatted')
c      write(16,'(A)') '----- CALCFORCE: -------------------------'
c      write(16,'(A,1x,F10.6)') 'Initialization        :',clock(2)
c      write(16,'(A,1x,F10.6)') 'Do loop               :',clock(8)      
c      write(16,'(A,1x,F10.6)') '  Min,max vertexc     :',clock(3)
c      write(16,'(A,1x,F10.6)') '  Pressure calculation:',clock(4)
c      write(16,'(A,1x,F10.6)') '  Shear calculation   :',clock(5)
c      write(16,'(A,1x,F10.6)') 'Input/Output          :',clock(5)
c      write(16,'(A,1x,F10.6)') 'Total time            :',clock(1)
c      write(16,'(A)') '------------------------------------------'
c      close(16)

      RETURN

      END
C------------------------------------------------------------------------


