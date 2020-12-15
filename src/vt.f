c
c-----------------------------------------------------------------------
c                 ***************************                         
c                 *         turvis.f        *                        
c                 ***************************                       
c----------------------------------------------------------------------- 
c
c	- Turvis:       calls the different eddy viscosity routines
c
c----------------------------------------------------------------------- 
c
c
c-----SUBROUTINE-Turvis------------------------E. Balaras 7/12/98-------
c
c
      SUBROUTINE TURVIS(UO,VO,WO,TV,G,GB,SXX,SYY,SZZ,SXY,SYZ,SXZ,
     &     MXX,MYY,MZZ,MXY,MYZ,MXZ,ILM,IMM,LM,MM,UC,VC,WC,
     &     DXDYDZ,CLES,CLESP,CLESN,NILMP,NILMN,XC,YC,ZC,DT,ICYCLE,IM,JM,KM)

c          Computes turbulent viscosity at the cell center 
c          using Lagrangian Dynamic model         
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
c----------------------------------------------------------------------
c
      INTEGER IM,JM,KM
      INTEGER ICYCLE
      REAL    XC(IM),YC(JM),ZC(KM)
      REAL    UO(IM,JM,KM),VO(IM,JM,KM),WO(IM,JM,KM)
      REAL    TV(IM,JM,KM),G(IM,JM,KM),GB(IM,JM,KM)
      REAL    SXX(IM,JM,KM),SYY(IM,JM,KM),SZZ(IM,JM,KM)
      REAL    SXY(IM,JM,KM),SXZ(IM,JM,KM),SYZ(IM,JM,KM)
      REAL    MXX(IM,JM,KM),MYY(IM,JM,KM),MZZ(IM,JM,KM)
      REAL    MXY(IM,JM,KM),MYZ(IM,JM,KM),MXZ(IM,JM,KM)
      REAL    ILM(IM,JM,KM),IMM(IM,JM,KM)
      REAL    LM(IM,JM,KM),MM(IM,JM,KM)
      REAL    DXDYDZ(IM,JM,KM),CLES(IM,JM,KM)
      REAL    CLESP(IM,JM,KM),CLESN(IM,JM,KM)
      INTEGER NILMP(IM,JM,KM),NILMN(IM,JM,KM)
c
      REAL    UC(IM,JM,KM),VC(IM,JM,KM),WC(IM,JM,KM)

      REAL    DT
c
      INTEGER I,J,K
      REAL    MMI(IM),LMI(IM),MMO(IM),LMO(IM)
      REAL    ALPHA,RPLUS,UTAU,C
c
c----------------------------------------------------------------------
c                      square of the ratio between test and grid filter
c----------------------------------------------------------------------
      REAL,PARAMETER :: R = 6.
c
c----------------------------------------------------------------------
c                       calculate Strain Rates and |S| at cell centers
c----------------------------------------------------------------------
c
c....Centered velocity components
c
      UC(IX1:IX2,:,:) = 0.5*(UO(1:IX2-1,:,:)+UO(IX1:IX2,:,:))
      UC(1 ,:,:) = 2.*UO(1  ,:,:)-UC(IX1,:,:)
      UC(IM,:,:) = 2.*UO(IX2,:,:)-UC(IX2,:,:)
c
      VC(:,JY1:JY2,:) = 0.5*(VO(:,1:JY2-1,:)+VO(:,JY1:JY2,:))
      VC(:,1 ,:) = 2.*VO(:,1  ,:)-VC(:,JY1,:)
      VC(:,JM,:) = 2.*VO(:,JY2,:)-VC(:,JY2,:)
c
      WC(:,:,KZ1:KZ2) = 0.5*(WO(:,:,1:KZ2-1)+WO(:,:,KZ1:KZ2))
      WC(:,:,1 ) = 2.*WO(:,:,1  )-WC(:,:,KZ1)
      WC(:,:,KM) = 2.*WO(:,:,KZ2)-WC(:,:,KZ2)
c
c in the case of cylindrical coordinates the values at I=1 come from the
c symmetric positions
c
      IF(ICYL==1) THEN
        DO J=1,JM
          UC(1,J,:) = -UC(2,JSYM(J),:)
          VC(1,J,:) = -VC(2,JSYM(J),:)
          WC(1,J,:) = WC(2,JSYM(J),:)
        ENDDO
      ENDIF

      CALL REFRESHBC(UC,IM*JM,KM)
      CALL REFRESHBC(VC,IM*JM,KM)
      CALL REFRESHBC(WC,IM*JM,KM)
c
c....Elements of the strain tensor at the cell centers
c....In G the value of |S| is stored
c
      IF(ICYL==0) THEN
        CALL STRAIN(UC,VC,WC,SXX,SYY,SZZ,SXY,SXZ,SYZ,G,IM,JM,KM) 
      ELSE
        CALL STRAINCYL(UC,VC,WC,SXX,SYY,SZZ,SXY,SXZ,SYZ,G,IM,JM,KM) 
      ENDIF

c.....Smagorinsky model

      IF(ISGS==1) THEN

!        DO I=IX1,IX2
!c
!c.....the van Driest correction is evaluated
!c
!        ALPHA = 1.0
!        IF(LP==0 .AND. NP==0) THEN
!c
!c.....case with periodic conditions along the azimuthal and the axial directions
!c UTAU: friction velocity
!c RPLUS: distance from the wall in wall units
!c ALPHA: van Driest damping function
!c 
!          IF(ICYL==1) THEN
!            UTAU = 0.5*SQRT(-DPDZ)
!            RPLUS = UTAU*(RU(IX2)-RP(I))/RU1
!            ALPHA = (1.-EXP(-RPLUS/25.))**2
!          ELSEIF(ICYL==0) THEN
!            UTAU = SQRT(-DPDZ)
!            RPLUS = UTAU*(1-ABS(XC(I)))/RU1
!            ALPHA = (1.-EXP(-RPLUS/25.))**2
!          ENDIF
!        ENDIF
!
!        TV(I,:,:) = CSS*CSS*ALPHA*DXDYDZ(I,:,:)*G(I,:,:)
!        CLES(I,:,:) = CSS*CSS*ALPHA
!        ENDDO

        ALPHA = 1.0

        IF(LP==0 .AND. NP==0) THEN

          UTAU = 0.5*SQRT(-DPDZ)

          IF(ICYL==1) THEN

            DO I=IX1,IX2

              RPLUS = UTAU*ABS(0.5-RP(I))/RU1
              ALPHA = (1.-EXP(-RPLUS/25.))**2

              TV(I,:,:) = CSS*CSS*ALPHA*DXDYDZ(I,:,:)*G(I,:,:)
              CLES(I,:,:) = CSS*CSS*ALPHA
            
            ENDDO

          ELSEIF(ICYL==0) THEN

            DO I=IX1,IX2
              DO J=JY1,JY2
c 
                RPLUS = UTAU*ABS(0.5-SQRT(XC(I)**2.+YC(J)**2.))/RU1
                ALPHA = (1.-EXP(-RPLUS/25.))**2

                TV(I,J,:) = CSS*CSS*ALPHA*DXDYDZ(I,J,:)*G(I,J,:)
                CLES(I,J,:) = CSS*CSS*ALPHA

              ENDDO
            ENDDO
 
          ENDIF

        ELSE

          TV(:,:,:) = CSS*CSS*ALPHA*DXDYDZ(:,:,:)*G(:,:,:)
          CLES(:,:,:) = CSS*CSS*ALPHA

        ENDIF

      ELSE
c
c A dynamic model is used
c
        IF(LP==0) THEN
          SXX(:,1 ,:) = SXX(:,JY2,:)
          SXX(:,JM,:) = SXX(:,2  ,:)
          SYY(:,1 ,:) = SYY(:,JY2,:)
          SYY(:,JM,:) = SYY(:,2  ,:)
          SZZ(:,1 ,:) = SZZ(:,JY2,:)
          SZZ(:,JM,:) = SZZ(:,2  ,:)
          SXY(:,1 ,:) = SXY(:,JY2,:)
          SXY(:,JM,:) = SXY(:,2  ,:)
          SXZ(:,1 ,:) = SXZ(:,JY2,:)
          SXZ(:,JM,:) = SXZ(:,2  ,:)
          SYZ(:,1 ,:) = SYZ(:,JY2,:)
          SYZ(:,JM,:) = SYZ(:,2  ,:)
          G  (:,1 ,:) = G  (:,JY2,:)
          G  (:,JM,:) = G  (:,2  ,:)
        ENDIF
      
        CALL REFRESHBC(SXX,IM*JM,KM)
        CALL REFRESHBC(SYY,IM*JM,KM)
        CALL REFRESHBC(SZZ,IM*JM,KM)
        CALL REFRESHBC(SXY,IM*JM,KM)
        CALL REFRESHBC(SYZ,IM*JM,KM)
        CALL REFRESHBC(SXZ,IM*JM,KM)
        CALL REFRESHBC(G  ,IM*JM,KM)

c----------------------------------------------------------------------
c                                   calculate numerator and denominator
c----------------------------------------------------------------------

c.....|S|Sxx --> mxx; |S|Syy --> myy; |S|Szz --> mzz; 
c.....|S|Sxy --> mxy; |S|Syz --> myz; |S|Sxz --> mxz
        MXX = G*SXX
        MYY = G*SYY
        MZZ = G*SZZ
        MXY = G*SXY
        MYZ = G*SYZ
        MXZ = G*SXZ

c.....\hat{|S|Sxx} --> mxx; \hat{|S|Syy} --> myy; \hat{|S|Szz} --> mzz;
c.....\hat{|S|Sxy} --> mxy; \hat{|S|Syz} --> myz; \hat{|S|Sxz} --> mxz
        CALL FILTER3D(MXX,MXX,LM,TV,IM,JM,KM) 
        CALL FILTER3D(MYY,MYY,LM,TV,IM,JM,KM) 
        CALL FILTER3D(MZZ,MZZ,LM,TV,IM,JM,KM) 
        CALL FILTER3D(MXY,MXY,LM,TV,IM,JM,KM) 
        CALL FILTER3D(MYZ,MYZ,LM,TV,IM,JM,KM) 
        CALL FILTER3D(MXZ,MXZ,LM,TV,IM,JM,KM) 

c.....\hat{Sxx} --> sxx; \hat{Syy} --> syy; \hat{Szz} --> szz; 
c.....\hat{Sxy} --> sxy; \hat{Syz} --> syz; \hat{Sxz} --> sxz
        CALL FILTER3D(SXX,SXX,LM,TV,IM,JM,KM)
        CALL FILTER3D(SYY,SYY,LM,TV,IM,JM,KM)
        CALL FILTER3D(SZZ,SZZ,LM,TV,IM,JM,KM)
        CALL FILTER3D(SXY,SXY,LM,TV,IM,JM,KM)
        CALL FILTER3D(SYZ,SYZ,LM,TV,IM,JM,KM)
        CALL FILTER3D(SXZ,SXZ,LM,TV,IM,JM,KM)

c.....\hat{|S|} --> gb
        CALL FILTER3D(G,GB,LM,TV,IM,JM,KM) 

c.....\Delta**2*(6*\hat{|S|}*\hat{Sij}-\hat{|S|Sij})
        C = 1.0
c        C= -2.0
        MXX = C*DXDYDZ*(R*GB*SXX-MXX)
        MYY = C*DXDYDZ*(R*GB*SYY-MYY)
        MZZ = C*DXDYDZ*(R*GB*SZZ-MZZ)
        MXY = C*DXDYDZ*(R*GB*SXY-MXY)
        MYZ = C*DXDYDZ*(R*GB*SYZ-MYZ)
        MXZ = C*DXDYDZ*(R*GB*SXZ-MXZ)

c.....MijMij --> mm
        MM = MXX**2+MYY**2+MZZ**2+2.*(MXY**2+MYZ**2+MXZ**2)

c.....\hat(u) --> sxx; \hat(v) --> syy; \hat(w) --> szz
        CALL FILTER3D(UC,SXX,LM,TV,IM,JM,KM)
        CALL FILTER3D(VC,SYY,LM,TV,IM,JM,KM)
        CALL FILTER3D(WC,SZZ,LM,TV,IM,JM,KM)

c.....uu --> gb; \hat{uu} --> gb; 
c.....(\hat{uu}-\hat(u)\hat(u))*mxx --> mxx
        GB = UC*UC
        CALL FILTER3D(GB,GB,LM,TV,IM,JM,KM)
        MXX = (GB-SXX*SXX)*MXX

c.....vv --> gb; \hat{vv} --> gb; 
c.....(\hat{vv}-\hat(v)\hat(v))*myy --> myy
        GB = VC*VC
        CALL FILTER3D(GB,GB,LM,TV,IM,JM,KM)
        MYY = (GB-SYY*SYY)*MYY

c.....ww --> gb; \hat{ww} --> gb; 
c.....(\hat{ww}-\hat(w)\hat(w))*mzz --> mzz
        GB = WC*WC
        CALL FILTER3D(GB,GB,LM,TV,IM,JM,KM)
        MZZ = (GB-SZZ*SZZ)*MZZ

c.....uv --> gb; \hat{uv} --> gb; 
c.....(\hat{uv}-\hat(u)\hat(v))*mxy --> mxy
        GB = UC*VC
        CALL FILTER3D(GB,GB,LM,TV,IM,JM,KM)
        MXY = (GB-SXX*SYY)*MXY

c.....vw --> gb; \hat{vw} --> gb; 
c.....(\hat{vw}-\hat(v)\hat(w))*myz --> myz
        GB = VC*WC
        CALL FILTER3D(GB,GB,LM,TV,IM,JM,KM)
        MYZ=(GB-SYY*SZZ)*MYZ

c.....uw --> gb; \hat{uw} --> gb; 
c.....(\hat{uw}-\hat(u)\hat(w))*mxz --> mxz
        GB = UC*WC
        CALL FILTER3D(GB,GB,LM,TV,IM,JM,KM)
        MXZ = (GB-SXX*SZZ)*MXZ

c.....LijMij --> lm
        LM = MXX+MYY+MZZ+2.*(MXY+MYZ+MXZ)

c----------------------------------------------------------------------
c                                         Homogeneous Direction Average
c----------------------------------------------------------------------

        IF(ISGS>1 .AND. ISGS<4) THEN
          
          IF(ISGS==3) ISGS=4  ! Lagrangian model, but restart with homogeneous direction average

c.....periodic boundary in y direction
c.....average along the Y direction
          IF(LP==0.AND.NP/=0) THEN
c         IF(.TRUE.) THEN

            WHERE (LM<0.) LM = 0.
c            WHERE (LM>0.) LM = 0.
c
            DO K=KZ1,KZ2
c            DO J=JY1,JY2
            DO I=IX1,IX2
              IMM(I,:,K) = SUM(MM(I,JY1:JY2,K))/REAL(JM-2)
              ILM(I,:,K) = SUM(LM(I,JY1:JY2,K))/REAL(JM-2)
c              IMM(I,J,:) = SUM(MM(I,J,KZ1:KZ2))/REAL(KM-2)
c              ILM(I,J,:) = SUM(LM(I,J,KZ1:KZ2))/REAL(KM-2)
            ENDDO
            ENDDO           
C
          ENDIF
C
c.....periodic boundary in y and z directions
c.....average along the Y and Z directions
          IF(LP==0.AND.NP==0) THEN
c          IF(.FALSE.) THEN
C
            WHERE (LM<0.) LM = 0.
c            WHERE (LM>0.) LM = 0.
c
            DO I=IX1,IX2
              MMI(I) = SUM(MM(I,JY1:JY2,KZ1:KZ2))/REAL((JM-2)*(KM-2))
              LMI(I) = SUM(LM(I,JY1:JY2,KZ1:KZ2))/REAL((JM-2)*(KM-2))
            ENDDO
C
            CALL MPI_ALLREDUCE(MMI,MMO,IM,MTYPE,MPI_SUM,MPI_COMM_EDDY,IERR)
            CALL MPI_ALLREDUCE(LMI,LMO,IM,MTYPE,MPI_SUM,MPI_COMM_EDDY,IERR)
C
            MMO = MMO/REAL(MYSIZE)!/REAL((JM-2)*(KM-2))
            LMO = LMO/REAL(MYSIZE)!/REAL((JM-2)*(KM-2))
C
            DO I=IX1,IX2
              IMM(I,:,:) = MMO(I)
              ILM(I,:,:) = LMO(I)
            ENDDO
c
          ENDIF

c----------------------------------------------------------------------
c                                                    Lagrangian Average
c----------------------------------------------------------------------
        ELSEIF(ISGS==4) THEN
c
c...The values of ILM and IMM are averaged along streamlines and
c...stored in MXX and MYY
          IF(.TRUE.) THEN
c          IF(MYSIZE==1) THEN
c            CALL LAGRANGE(ILM,IMM,MXX,MYY,UC,VC,WC,DT,IM,JM,KM)
c          ELSE
            CALL LAGRANGE_PARAL(ILM,IMM,MXX,MYY,UC,VC,WC,XC,YC,ZC,DT,IM,JM,KM)
c          ENDIF
c...Calculate T^n, eps

          TV = ILM*IMM
c          TV = -6.0*ILM*IMM
          GB = 0.
          WHERE(TV>EPS) GB = 1./(1.+1.5*SQRT(DXDYDZ)*TV**(-0.125)/DT)

c the new values of ILM and IMM are estimated
          ILM = GB*LM+(1.-GB)*MXX
          IMM = GB*MM+(1.-GB)*MYY

c          WHERE(ILM<0) ILM = 0.

          ELSE
c...The convective derivatives of ILM and IMM are evaluated and stored in MXX and MYY
            CALL CONV_DERV_SCALAR(MXX,ILM,UC,VC,WC,IM,JM,KM)
            CALL CONV_DERV_SCALAR(MYY,IMM,UC,VC,WC,IM,JM,KM)
            TV = ILM*IMM

            GB = 0.
            WHERE(TV>EPS) GB = 1./(1.+1.5*SQRT(DXDYDZ)*TV**(-0.125)/DT)

            UC = ILM
            VC = IMM
c the new values of ILM and IMM are estimated
            ILM = GB*LM + (1.0-GB)*(ILM-DT*MXX)
            IMM = GB*MM + (1.0-GB)*(IMM-DT*MYY)

c            WHERE(ILM<0) ILM = 0.

          ENDIF
        ENDIF
c 
c...calculate turbulent viscosity

        GB = 0.0
        CLESP = 0.0
        CLESN = 0.0
        WHERE(IMM>SQRT(EPS)) GB = ABS(ILM/IMM)
        WHERE(IMM>SQRT(EPS)) CLESP = 0.5*((ILM/IMM)+GB)
        WHERE(IMM>SQRT(EPS)) CLESN = 0.5*((ILM/IMM)-GB)        

        nilmn = 0
        nilmp = 0
        WHERE(ILM<0.0) NILMN = 1
        WHERE(ILM>=0.0) NILMP = 1
c
c...ILM and IMM at the ghost layers are updated
        CALL REFRESHBC(ILM,IM*JM,KM)
        CALL REFRESHBC(IMM,IM*JM,KM)
c
c...ILM and IMM at the boundaries of the computational domain are updated
        CALL BOUNDTURB(ILM,IM,JM,KM,0)
        CALL BOUNDTURB(IMM,IM,JM,KM,1)

c        write(6,*) 'min, max',minval(ilm),maxval(ilm),minval(imm),maxval(imm)

        WHERE(ILM<0) ILM = 0.

        TV = 0.
        CLES = 0.
c the value of the eddy viscosity is estimated
c if TV is negative it is not physical and is enforced equal to 0
        WHERE(IMM>SQRT(EPS)) TV = (ILM/IMM)*DXDYDZ*G
c        WHERE(IMM>SQRT(EPS)) TV = -0.5*(ILM/IMM)*DXDYDZ*G
        WHERE(IMM>SQRT(EPS)) CLES = (ILM/IMM)
c        WHERE(IMM>SQRT(EPS)) CLES = -0.5*(ILM/IMM)
c        WHERE(TV<-1.0*RU1) TV = -1.0*RU1
        WHERE(TV<0.0) TV = 0.0

!!!!!!        TV = 0.0
C
      ENDIF
C
c...eddy viscosity and LES coefficient at the ghost layers are updated
      CALL REFRESHBC(TV,IM*JM,KM)
      CALL REFRESHBC(CLES,IM*JM,KM)
c
c...eddy viscosity and LES coefficient at the boundaries of the computational domain are updated
      CALL BOUNDTURB(TV,IM,JM,KM,0)
      CALL BOUNDTURB(CLES,IM,JM,KM,0)
c
      RETURN
      END
c
c
c
c---------------------------------------------------------------------
      SUBROUTINE LAGRANGE(VIN1,VIN2,VOUT1,VOUT2,UC,VC,WC,DT,IM,JM,KM)
c
c     given intergal at i,j,k returns value at (n-1) along 
c     fluid pathline
c
      include 'common.h'
c
      INTEGER IM,JM,KM
      REAL    VIN1(IM,JM,KM),VOUT1(IM,JM,KM)
      REAL    VIN2(IM,JM,KM),VOUT2(IM,JM,KM)
      REAL    UC(IM,JM,KM),VC(IM,JM,KM),WC(IM,JM,KM)
      REAL    DT
C
      INTEGER I,J,K,II,JJ,KK,IC,JC,KC,IP,JP,KP
      REAL    FXC,FYC,FZC,FXP,FYP,FZP
      REAL    FN1,FN2,FN3,FN4,FN5,FN6,FN7,FN8
c
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
c.....indices shift
        II = INT(SIGN(1.,-UC(I,J,K)))
        JJ = INT(SIGN(1.,-VC(I,J,K)))
        KK = INT(SIGN(1.,-WC(I,J,K)))
c.....particle location
        IC = MIN(I,I+II)
        JC = MIN(J,J+JJ)
        KC = MIN(K,K+KK)
        IP = IC+1
        JP = JC+1
        KP = KC+1
c.....particle is located between (ic:ip,jc:jp,kc:kp)
c.....interpolation factors
        FXP = 0.5*REAL(1-II)-UC(I,J,K)*DT*AU(IC)
        FYP = 0.5*REAL(1-JJ)-VC(I,J,K)*DT*BV(JC)/RP(I)
        FZP = 0.5*REAL(1-KK)-WC(I,J,K)*DT*CW(KC)
        FXC = 1.-FXP
        FYC = 1.-FYP
        FZC = 1.-FZP
c.....weights:    point 1: ic,jc,kc;     point 2: ic,jp,kc
c                 point 3: ic,jp,kp;     point 4: ic,jc,kp
c                 point 5: ip,jc,kc;     point 6: ip,jp,kc
c                 point 7: ip,jp,kp;     point 8: ip,jc,kp
        FN1 = FXC*FYC*FZC
        FN2 = FXC*FYP*FZC
        FN3 = FXC*FYP*FZP
        FN4 = FXC*FYC*FZP
        FN5 = FXP*FYC*FZC
        FN6 = FXP*FYP*FZC
        FN7 = FXP*FYP*FZP
        FN8 = FXP*FYC*FZP
c     values
        VOUT1(I,J,K) = VIN1(IC,JC,KC)*FN1+VIN1(IC,JP,KC)*FN2
     &       +         VIN1(IC,JP,KP)*FN3+VIN1(IC,JC,KP)*FN4
     &       +         VIN1(IP,JC,KC)*FN5+VIN1(IP,JP,KC)*FN6
     &       +         VIN1(IP,JP,KP)*FN7+VIN1(IP,JC,KP)*FN8
C
        VOUT2(I,J,K) = VIN2(IC,JC,KC)*FN1+VIN2(IC,JP,KC)*FN2
     &       +         VIN2(IC,JP,KP)*FN3+VIN2(IC,JC,KP)*FN4
     &       +         VIN2(IP,JC,KC)*FN5+VIN2(IP,JP,KC)*FN6
     &       +         VIN2(IP,JP,KP)*FN7+VIN2(IP,JC,KP)*FN8

c        if(i==2 .AND. j==2 .AND. k==2) then
c          write(6,*) i,j,k,ic,jc,kc,fxp,fyp,fzp,-UC(I,J,K)*DT,AU(IC)
c     &          ,vout1(i,j,k),vout2(i,j,k)
c     &          ,vin1(ic:ip,jc:jp,kc:kp),vin2(ic:ip,jc:jp,kc:kp)
c          stop
c        endif

      ENDDO
      ENDDO
      ENDDO
c
      RETURN
      END
c
c
c---------------------------------------------------------------------
      SUBROUTINE LAGRANGE_PARAL(VIN1,VIN2,VOUT1,VOUT2,UC,VC,WC,XC,YC,ZC,DT,NX,NY,NZ)
c
c     given intergal at i,j,k returns value at (n-1) along 
c     fluid pathline
c
      include 'common.h'
c      include 'mpif.h'
c
c.... Input/Output arrays
      INTEGER NX,NY,NZ
      REAL    XC(NX),YC(NY),ZC(NZ)
      REAL    VIN1(NX,NY,NZ),VOUT1(NX,NY,NZ)
      REAL    VIN2(NX,NY,NZ),VOUT2(NX,NY,NZ)
      REAL    UC(NX,NY,NZ),VC(NX,NY,NZ),WC(NX,NY,NZ)
      REAL    DT
C
c.... Local arrays
      INTEGER I,J,K,II,JJ,KK,IC,JC,KC,IP,JP,KP,NMAX,NL,NR,ILU,IRU,NL_REQ,NR_REQ
      REAL    XP,YP,ZP,RCYL,RCAR,A,B
      real, dimension(:), allocatable :: ul,ur
      real, dimension(:,:), allocatable :: ucordl,ucordr
      integer, dimension(:), allocatable :: iordl,jordl,kordl,iordr,jordr,kordr
c
c.... Functions
      real    interp3d
c
      nmax = nx*ny
      nl = nmax
      nr = nmax
      allocate(ucordl(3,nl),ucordr(3,nr))
      allocate(iordl(nl),jordl(nl),kordl(nl),iordr(nr),jordr(nr),kordr(nr))

      rcyl = icyl
      rcar = 1.0-icyl

      ilu = 0
      iru = 0

      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
        XP = XC(I) - UC(I,J,K)*DT
        YP = YC(J) - VC(I,J,K)*DT/RP(I)
        ZP = ZC(K) - WC(I,J,K)*DT

        IF(YP<YC(1)) THEN
          YP = YC(NY)-(YC(2)-YP)
        ELSEIF(YP>=YC(NY)) THEN
          YP = YC(2)+(YP-YC(NY))
        ENDIF

        IF(XP>XC(NX)) THEN
          IF(ITYPE(2)==500) THEN
            XP = XC(2) + (XP-XC(NX))
          ELSE
            XP = XC(NX)-1.0e-8
          ENDIF
        ENDIF

        IF(XP<XC(1)) THEN
          IF(ITYPE(1)==500) THEN
            XP = XC(NX)-(XC(2)-XP)
          ELSE
            XP = XC(1)*RCAR + XP*RCYL
          ENDIF
        ENDIF

        !Convert cylindrical to cartesian coordinates
        a = xp
        b = yp
        xp = a*rcar + a*cos(b)*rcyl
        yp = b*rcar + a*sin(b)*rcyl

        IF(ZP<ZC(1)) THEN
          IF(MYSIZE==1 .AND. ITYPE(5)==500) THEN
            zp = zlen+zp
          ELSEIF(MYRANK==0 .AND. ITYPE(5)/=500) THEN
            ZP = ZC(1)
          ENDIF          
        ELSEIF(ZP>=ZC(NZ)) THEN
          IF(MYSIZE==1 .AND. ITYPE(6)==500) THEN
            zp = zp-zlen
          ELSEIF(MYRANK==MYSIZE-1 .AND. ITYPE(5)/=500) THEN
            ZP = ZC(NZ)-1.0e-8
          ENDIF
        ENDIF
        
        IF(ZP<ZC(1)) THEN
          if(myrank==0 .AND. itype(5)==500) zp = zlen+zp
          ilu = ilu+1
          ucordl(1,ilu) = xp
          ucordl(2,ilu) = yp
          ucordl(3,ilu) = zp
          iordl(ilu) = i
          jordl(ilu) = j
          kordl(ilu) = k
        ELSEIF(ZP>=ZC(NZ)) THEN
          if(myrank==mysize-1 .AND. itype(6)==500) zp = zp-zlen 
          iru = iru+1
          ucordr(1,iru) = xp
          ucordr(2,iru) = yp
          ucordr(3,iru) = zp
          iordr(iru) = i
          jordr(iru) = j
          kordr(iru) = k
        ELSE
          vout1(i,j,k) = interp3d(xp,yp,zp,xc,yc,zc,vin1,nx,ny,nz,icyl)
          vout2(i,j,k) = interp3d(xp,yp,zp,xc,yc,zc,vin2,nx,ny,nz,icyl)
        ENDIF

      ENDDO
      ENDDO
      ENDDO

c      write(6,*) 'myrank=',myrank,', ilu=',ilu,', iru=',iru
c      return

      nl_req = ilu
      nr_req = iru

      if(nl_req==0) then
        allocate(ul(1))
      else
        allocate(ul(nl_req))
      endif
        
      if(nr_req==0) then
        allocate(ur(1))
      else
        allocate(ur(nr_req))
      endif

      call  mpi_var_interp(ucordl,ul,nl_req,ucordr,ur,nr_req,vin1,xc,yc,zc,nx,ny,nz)
     
      do ii=1,nl_req
        i = iordl(ii)
        j = jordl(ii)
        k = kordl(ii)        
        vout1(i,j,k) = ul(ii)
      enddo

      do ii=1,nr_req
        i = iordr(ii)
        j = jordr(ii)
        k = kordr(ii)
        vout1(i,j,k) = ur(ii)
      enddo

      call  mpi_var_interp(ucordl,ul,nl_req,ucordr,ur,nr_req,vin2,xc,yc,zc,nx,ny,nz)
     
      do ii=1,nl_req
        i = iordl(ii)
        j = jordl(ii)
        k = kordl(ii)        
        vout2(i,j,k) = ul(ii)
      enddo

      do ii=1,nr_req
        i = iordr(ii)
        j = jordr(ii)
        k = kordr(ii)
        vout2(i,j,k) = ur(ii)
      enddo

      deallocate(ur,ul)
      deallocate(ucordl,ucordr)
      deallocate(iordl,jordl,kordl,iordr,jordr,kordr)
c
      RETURN
      END
c
c-----------------------------------------------------------------------



      SUBROUTINE BOUNDTURB(VT,IM,JM,KM,IFLG)
c
c-----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c-----------------------------------------------------------------------
c
      INTEGER IM,JM,KM
      INTEGER IFLG
      REAL    VT(IM,JM,KM)
c
      INTEGER IBND
c
c----------------------------------------------------------------------
c                        set boundary values for turbulent viscosity
c----------------------------------------------------------------------     
      
      DO IBND=1,6
c
        IF(IFLG==0) THEN
c.....axis
          IF(ITYPE(IBND)==300) THEN
            CALL BOUNDAXIS(VT,IM,JM,KM,IBND)
c.....free-stream boundaries
          ELSEIF(ITYPE(IBND)>700) THEN
            CALL BOUNDGRAD(VT,IM,JM,KM,IBND)
c.....wall
          ELSEIF(ITYPE(IBND)/=0.AND.ITYPE(IBND)/=500) THEN
            CALL BOUNDWALL(VT,IM,JM,KM,IBND)
          ENDIF
        ENDIF
c.....others except parallel or periodic boundaries
        IF(IFLG==1) THEN
c.....axis
          IF(ITYPE(IBND)==300) THEN
            CALL BOUNDAXIS(VT,IM,JM,KM,IBND)
c.....free-stream boundaries
          ELSEIF(ITYPE(IBND)/=0.AND.ITYPE(IBND)/=500) THEN
            CALL BOUNDGRAD(VT,IM,JM,KM,IBND)
          ENDIF
        ENDIF

      ENDDO
c
c.....periodic boundaries
      DO IBND=1,4

        IF(ITYPE(IBND)==500)
     &       CALL BOUNDPERIOD(VT,IM,JM,KM,IBND)

      ENDDO
c      
      RETURN
      END
c
c
c---------------------------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE BOUNDWALL(vt,im,jm,km,ibnd)

c       Sets the boundary values of nu_t for the case 
c       of a solid wall 
c
c           i.e.  nu_t(i,j,1)=-nu_t(i,j,2) 
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER im,jm,km
      REAL    vt(im,jm,km)
c
      INTEGER ibnd
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (1) 
        VT(1 ,:,:) = -VT(2  ,:,:)
         
      CASE(2)
        VT(IM,:,:) = -VT(IX2,:,:)
         
      CASE(3)
        VT(:,1 ,:) = -VT(:,2  ,:)
         
      CASE(4)
        VT(:,JM,:) = -VT(:,JY2,:)
        
      CASE(5)
        VT(:,:,1 ) = -VT(:,:,2  )
        
      CASE(6)
        VT(:,:,KM) = -VT(:,:,KZ2)
        
      END SELECT
c
      RETURN
      END
c
c
c---------------------------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE BOUNDGRAD(vt,im,jm,km,ibnd)

c       Sets the boundary values of nu_t for the case 
c       of a freestream boundaries 
c
c           i.e.  nu_t(i,j,1)= nu_t(i,j,2) 
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER im,jm,km
      REAL    vt(im,jm,km)
c
      INTEGER ibnd
*
**** Main loop
*
      SELECT CASE(IBND)

      CASE (1) 
        VT(1 ,:,:) = VT(2  ,:,:)
         
      CASE(2)
        VT(IM,:,:) = VT(IX2,:,:)
         
      CASE(3)
        VT(:,1 ,:) = VT(:,2  ,:)
         
      CASE(4)
        VT(:,JM,:) = VT(:,JY2,:)
        
      CASE(5)
        VT(:,:,1 ) = VT(:,:,2  )
        
      CASE(6)
        VT(:,:,KM) = VT(:,:,KZ2)
        
      END SELECT
c
      RETURN
      END
c
c
c---------------------------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE BOUNDPERIOD(vt,im,jm,km,ibnd)

c      Sets boundary values for nu_t for the case of
c      periodic boundaries
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER im,jm,km
      REAL    vt(im,jm,km)
c
      INTEGER ibnd
*
****  Main loop
*
      SELECT CASE(IBND)

      CASE (1) 
        VT(1 ,:,:) = VT(IX2,:,:)
         
      CASE(2)
        VT(IM,:,:) = VT(2  ,:,:)
         
      CASE(3)
        VT(:,1 ,:) = VT(:,JY2,:)
         
      CASE(4)
        VT(:,JM,:) = VT(:,2  ,:)
        
      CASE(5)
        VT(:,:,1 ) = VT(:,:,KZ2)
        
      CASE(6)
        VT(:,:,KM) = VT(:,:,2  )
        
      END SELECT
c
      RETURN
      END
c
c
c---------------------------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE BOUNDAXIS(vt,im,jm,km,ibnd)

c      Sets boundary values for nu_t for the case of
c      periodic boundaries
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER im,jm,km
      REAL    vt(im,jm,km)
c
      INTEGER ibnd,J
*
****  Main loop
*
      SELECT CASE(IBND)

      CASE(1) 
        DO J=1,JM
          VT(1 ,J,:) = VT(2  ,JSYM(J),:)
        ENDDO
         
      CASE(2)
        DO J=1,JM
          VT(IM,J,:) = VT(IX2,JSYM(J),:)
        ENDDO
         
      CASE(3)
        WRITE(6,*) 'WRONG BOUNDARY TYPE!'
        
      CASE(4)
        WRITE(6,*) 'WRONG BOUNDARY TYPE!'
        
      CASE(5)
        DO J=1,JM
          VT(:,J,1 ) = VT(:,JSYM(J),2  )
        ENDDO
        
      CASE(6)
        DO J=1,JM
          VT(:,J,KM) = VT(:,JSYM(J),KZ2)
        ENDDO
        
      END SELECT
c
      RETURN
      END
c
c
c-----SUBROUTINE-Strain-------------------------------------------------
c
      SUBROUTINE Strain(uc,vc,wc,sxx,syy,szz,sxy,sxz,syz,g,im,jm,km)
c
c       Calculate strain rates at cell centers
c
c
c-----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c-----------------------------------------------------------------------
c
      INTEGER IM,JM,KM
      REAL    UC(IM,JM,KM),VC(IM,JM,KM),WC(IM,JM,KM)
      REAL    SXX(IM,JM,KM),SYY(IM,JM,KM),SZZ(IM,JM,KM),
     &        SXY(IM,JM,KM),SXZ(IM,JM,KM),SYZ(IM,JM,KM)
      REAL    G(IM,JM,KM)
c
      INTEGER I,J,K
      REAL    DUDXP,DUDXM,DUDYP,DUDYM,DUDZP,DUDZM
      REAL    DVDXP,DVDXM,DVDYP,DVDYM,DVDZP,DVDZM
      REAL    DWDXP,DWDXM,DWDYP,DWDYM,DWDZP,DWDZM
c
c-----------------------------------------------------------------------
c
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
        DUDXP = AU(I  )*(UC(I+1,J,K)-UC(I  ,J,K))
        DUDXM = AU(I-1)*(UC(I  ,J,K)-UC(I-1,J,K))
        DUDYP = BV(J  )*(UC(I,J+1,K)-UC(I,J  ,K))
        DUDYM = BV(J-1)*(UC(I  ,J,K)-UC(I,J-1,K))
        DUDZP = CW(K  )*(UC(I,J,K+1)-UC(I,J,K  ))
        DUDZM = CW(K-1)*(UC(I,J,K  )-UC(I,J,K-1))

        DVDXP = AU(I  )*(VC(I+1,J,K)-VC(I  ,J,K))
        DVDXM = AU(I-1)*(VC(I  ,J,K)-VC(I-1,J,K))
        DVDYP = BV(J  )*(VC(I,J+1,K)-VC(I,J  ,K))
        DVDYM = BV(J-1)*(VC(I  ,J,K)-VC(I,J-1,K))
        DVDZP = CW(K  )*(VC(I,J,K+1)-VC(I,J,K  ))
        DVDZM = CW(K-1)*(VC(I,J,K  )-VC(I,J,K-1))

        DWDXP = AU(I  )*(WC(I+1,J,K)-WC(I  ,J,K))
        DWDXM = AU(I-1)*(WC(I  ,J,K)-WC(I-1,J,K))
        DWDYP = BV(J  )*(WC(I,J+1,K)-WC(I,J  ,K))
        DWDYM = BV(J-1)*(WC(I  ,J,K)-WC(I,J-1,K))
        DWDZP = CW(K  )*(WC(I,J,K+1)-WC(I,J,K  ))
        DWDZM = CW(K-1)*(WC(I,J,K  )-WC(I,J,K-1))

        SXX(I,J,K) = 0.5*(DUDXP+DUDXM)
        SYY(I,J,K) = 0.5*(DVDYP+DVDYM)
        SZZ(I,J,K) = 0.5*(DWDZP+DWDZM)
        SXY(I,J,K) = 0.25*(DUDYP+DUDYM+DVDXP+DVDXM)
        SXZ(I,J,K) = 0.25*(DUDZP+DUDZM+DWDXP+DWDXM)
        SYZ(I,J,K) = 0.25*(DWDYP+DWDYM+DVDZP+DVDZM)

        G(I,J,K) = SQRT(2.*(SXX(I,J,K)**2+SYY(I,J,K)**2+SZZ(I,J,K)**2)
     &       +          4.*(SXY(I,J,K)**2+SYZ(I,J,K)**2+SXZ(I,J,K)**2))

      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
c
c
c-----SUBROUTINE-Strain-------------------------------------------------
c
      SUBROUTINE StrainCyl(uc,vc,wc,sxx,syy,szz,sxy,sxz,syz,g,im,jm,km)
c
c       Calculate strain rates at cell centers
c
c
c-----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c-----------------------------------------------------------------------
c
      INTEGER IM,JM,KM
      REAL    UC(IM,JM,KM),VC(IM,JM,KM),WC(IM,JM,KM)
      REAL    SXX(IM,JM,KM),SYY(IM,JM,KM),SZZ(IM,JM,KM),
     &        SXY(IM,JM,KM),SXZ(IM,JM,KM),SYZ(IM,JM,KM)
      REAL    G(IM,JM,KM)
c
      INTEGER I,J,K
      REAL    DUDXP,DUDXM,DUDYP,DUDYM,DUDZP,DUDZM
      REAL    DVDXP,DVDXM,DVDYP,DVDYM,DVDZP,DVDZM
      REAL    DWDXP,DWDXM,DWDYP,DWDYM,DWDZP,DWDZM
c
c-----------------------------------------------------------------------
c
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
        DUDXP = AU(I  )*(UC(I+1,J,K)-UC(I  ,J,K))
        DUDXM = AU(I-1)*(UC(I  ,J,K)-UC(I-1,J,K))
        DUDYP = BV(J  )*(UC(I,J+1,K)-UC(I,J  ,K))
        DUDYM = BV(J-1)*(UC(I  ,J,K)-UC(I,J-1,K))
        DUDZP = CW(K  )*(UC(I,J,K+1)-UC(I,J,K  ))
        DUDZM = CW(K-1)*(UC(I,J,K  )-UC(I,J,K-1))

        DVDXP = AU(I  )*(VC(I+1,J,K)-VC(I  ,J,K))
        DVDXM = AU(I-1)*(VC(I  ,J,K)-VC(I-1,J,K))
        DVDYP = BV(J  )*(VC(I,J+1,K)-VC(I,J  ,K))
        DVDYM = BV(J-1)*(VC(I  ,J,K)-VC(I,J-1,K))
        DVDZP = CW(K  )*(VC(I,J,K+1)-VC(I,J,K  ))
        DVDZM = CW(K-1)*(VC(I,J,K  )-VC(I,J,K-1))

        DWDXP = AU(I  )*(WC(I+1,J,K)-WC(I  ,J,K))
        DWDXM = AU(I-1)*(WC(I  ,J,K)-WC(I-1,J,K))
        DWDYP = BV(J  )*(WC(I,J+1,K)-WC(I,J  ,K))
        DWDYM = BV(J-1)*(WC(I  ,J,K)-WC(I,J-1,K))
        DWDZP = CW(K  )*(WC(I,J,K+1)-WC(I,J,K  ))
        DWDZM = CW(K-1)*(WC(I,J,K  )-WC(I,J,K-1))

        SXX(I,J,K) = 0.5*(DUDXP+DUDXM)
        SYY(I,J,K) = 0.5*(DVDYP+DVDYM)/RP(I) + UC(I,J,K)/RP(I)
        SZZ(I,J,K) = 0.5*(DWDZP+DWDZM)
        SXY(I,J,K) = 0.25*((DUDYP+DUDYM)/RP(I) + (DVDXP+DVDXM)) - 0.5*VC(I,J,K)/RP(I)
        SXZ(I,J,K) = 0.25*(DUDZP+DUDZM+DWDXP+DWDXM)
        SYZ(I,J,K) = 0.25*((DWDYP+DWDYM)/RP(I)+DVDZP+DVDZM)

        G(I,J,K) = SQRT(2.*(SXX(I,J,K)**2+SYY(I,J,K)**2+SZZ(I,J,K)**2)
     &       +          4.*(SXY(I,J,K)**2+SYZ(I,J,K)**2+SXZ(I,J,K)**2))

      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
c
c
c---------------------------------------------E. Balaras 7/12/98-------

c
c
c---------------------------------------------E. Balaras 7/12/98-------

c
      SUBROUTINE FILTERSIZE(DXDYDZ,NX,NY,NZ)

c       Calculates filter size on the non-uniform grid 
c
c           delta^2 = (Dx_i * Dy_j * Dz_k)^(2/3)
c
c----------------------------------------------------------------------
c
      INCLUDE 'common.h'
c
c----------------------------------------------------------------------
c
      INTEGER NX,NY,NZ
      REAL    DXDYDZ(NX,NY,NZ)
      REAL,PARAMETER :: PPP = 2./3.
      INTEGER I,J,K
c
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
        DXDYDZ(I,J,K) = (ABS(RP(I))/(AP(I)*BP(J)*CP(K)))**PPP
      ENDDO
      ENDDO
      ENDDO

      IF(LP==0) THEN
        DXDYDZ(:,1 ,:) = DXDYDZ(:,NY-1,:)
        DXDYDZ(:,NY,:) = DXDYDZ(:,2   ,:)
      ENDIF
C
      RETURN
      END



c--------------------------------------------------------------
      subroutine filter3d(f,fb,wk1,wk2,im,jm,km)
c--------------------------------------------------------------
cc box filter routine (trapezodial rule for NON-UNIFORM grid)
cc On the boundaries backward and forward formulas are used
cc fb is not defined on the boundaries on output
cc 
c
      include 'common.h'
c
      integer im,jm,km
      real f(im,jm,km),fb(im,jm,km)
      real wk1(im,jm,km),wk2(im,jm,km)
      integer i,j,k

      if(icyl==1) then
        call filter3d_cyl(f,fb,im,jm,km)
        return 
      endif

c--- initialize work arrays
      wk1=0.0
      wk2=0.0

c--- filter in x ---
      DO K=1,KM
      DO J=1,JM

      do i=ix1+1,ix2-1
        wk1(i,j,k) = whx(i,-1)*f(i-1,j,k)
     &       +       whx(i, 0)*f(i  ,j,k)
     &       +       whx(i, 1)*f(i+1,j,k)
      enddo
c.....boundary points
      wk1(ix1,j,k) = whx(ix1,-1)*f(ix1  ,j,k)
     &     +         whx(ix1, 0)*f(ix1+1,j,k)
     &     +         whx(ix1, 1)*f(ix1+2,j,k)
      wk1(ix2,j,k) = whx(ix2,-1)*f(ix2  ,j,k)
     &     +         whx(ix2, 0)*f(ix2-1,j,k)
     &     +         whx(ix2, 1)*f(ix2-2,j,k)

      ENDDO
      ENDDO

c      wk1 = f

c--- filter in y ---
      IF(LP==0) THEN

        wk2(:,JY1:JY2,:) = 0.25*wk1(:,JY1-1:JY2-1,:)
     &       +             0.50*wk1(:,JY1  :JY2  ,:)
     &       +             0.25*wk1(:,JY1+1:JY2+1,:)

      ELSEIF(LP==1) THEN
                                                                                   
      DO K=1,KM
      DO I=1,IM

      do j=jy1+1,jy2-1      
        wk2(i,j,k) = why(j,-1)*wk1(i,j-1,k)
     &       +       why(j, 0)*wk1(i,j  ,k)
     &       +       why(j, 1)*wk1(i,j+1,k)
      enddo
c.....boundary points
      wk2(i,jy1,k) = why(jy1,-1)*wk1(i,jy1  ,k)
     &     +         why(jy1, 0)*wk1(i,jy1+1,k)
     &     +         why(jy1, 1)*wk1(i,jy1+2,k)
      wk2(i,jy2,k) = why(jy2,-1)*wk1(i,jy2  ,k)
     &     +         why(jy2, 0)*wk1(i,jy2-1,k)
     &     +         why(jy2, 1)*wk1(i,jy2-2,k)

      ENDDO
      ENDDO

      ENDIF

c      wk2 = wk1
c--- filter in z ---
      IF(NP==0) THEN

        fb(:,:,KZ1:KZ2) = 0.25*wk2(:,:,kZ1-1:KZ2-1)
     &       +            0.50*wk2(:,:,KZ1  :KZ2  )
     &       +            0.25*wk2(:,:,KZ1+1:KZ2+1)

      ELSEIF(NP==1) THEN      

      DO J=1,JM
      DO I=1,IM

      do k=kz1+1,kz2-1      
        fb(i,j,k) = whz(k,-1)*wk2(i,j,k-1)
     &       +      whz(k, 0)*wk2(i,j,k  )
     &       +      whz(k, 1)*wk2(i,j,k+1)
      enddo

      ENDDO
      ENDDO
c.....boundary points
      IF(itype(5)>0) THEN       !BOTTOM BOUNDARY
        DO J=1,JM
        DO I=1,IM
          fb(i,j,kz1) = whz(kz1,-1)*wk2(i,j,kz1  )
     &         +        whz(kz1, 0)*wk2(i,j,kz1+1)
     &         +        whz(kz1, 1)*wk2(i,j,kz1+2)
        ENDDO
        ENDDO
      ELSE                      !GHOST LAYER
        DO J=1,JM
        DO I=1,IM
          fb(i,j,kz1) = whz(kz1,-1)*wk2(i,j,kz1-1)
     &         +        whz(kz1, 0)*wk2(i,j,kz1  )
     &         +        whz(kz1, 1)*wk2(i,j,kz1+1)
        ENDDO
        ENDDO
      ENDIF
c
      IF(itype(6)>0) THEN       !TOP BOUNDARY
        DO J=1,JM
        DO I=1,IM
          fb(i,j,kz2) = whz(kz2,-1)*wk2(i,j,kz2  )
     &         +        whz(kz2, 0)*wk2(i,j,kz2-1)
     &         +        whz(kz2, 1)*wk2(i,j,kz2-2)
        ENDDO
        ENDDO
      ELSE                      !GHOST LAYER
        DO J=1,JM
        DO I=1,IM
          fb(i,j,kz2) = whz(kz2,-1)*wk2(i,j,kz2-1)
     &         +        whz(kz2, 0)*wk2(i,j,kz2  )
     &         +        whz(kz2, 1)*wk2(i,j,kz2+1)
        ENDDO
        ENDDO
      ENDIF

      ENDIF
c
      return
      end
c-------------------------------------------------------------

c--------------------------------------------------------------
      subroutine filter3d_cyl(f,fb,nx,ny,nz)
c--------------------------------------------------------------
cc box filter routine (trapezodial rule for NON-UNIFORM grid)
cc On the boundaries backward and forward formulas are used
cc fb is not defined on the boundaries on output
cc 
c
      include 'common.h'
c
      integer nx,ny,nz
      real f(nx,ny,nz),fb(nx,ny,nz)
      real wk1(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k
      real    a,f1,f2,f3,f4,f5,f6,f7,f8,dV
c
      a = 1./8.

      wk1 = f

      do k=kz1,kz2
      do j=jy1,jy2
      do i=ix1,ix2
! 1 / volume of the computational cell (i,j,k) / number of the positions (8)
        dV = av(i)*bv(j)*cv(k)/(8.0*rp(i))  

        f1 = a*(f(i-1,j-1,k-1)+f(i-1,j-1,k  )+f(i-1,j  ,k-1)+f(i-1,j  ,k  )
     &        + f(i  ,j-1,k-1)+f(i  ,j-1,k  )+f(i  ,j  ,k-1)+f(i  ,j  ,k  ))

        f2 = a*(f(i-1,j  ,k-1)+f(i-1,j  ,k  )+f(i-1,j+1,k-1)+f(i-1,j+1,k  )
     &        + f(i  ,j  ,k-1)+f(i  ,j  ,k  )+f(i  ,j+1,k-1)+f(i  ,j+1,k  ))

        f3 = a*(f(i-1,j-1,k  )+f(i-1,j-1,k+1)+f(i-1,j  ,k  )+f(i-1,j  ,k+1)
     &        + f(i  ,j-1,k  )+f(i  ,j-1,k+1)+f(i  ,j  ,k  )+f(i  ,j  ,k+1))

        f4 = a*(f(i-1,j  ,k  )+f(i-1,j  ,k+1)+f(i-1,j+1,k  )+f(i-1,j+1,k+1)
     &        + f(i  ,j  ,k  )+f(i  ,j  ,k+1)+f(i  ,j+1,k  )+f(i  ,j+1,k+1))

        f5 = a*(f(i  ,j-1,k-1)+f(i  ,j-1,k  )+f(i  ,j  ,k-1)+f(i  ,j  ,k  )
     &        + f(i+1,j-1,k-1)+f(i+1,j-1,k  )+f(i+1,j  ,k-1)+f(i+1,j  ,k  ))

        f6 = a*(f(i  ,j  ,k-1)+f(i  ,j  ,k  )+f(i  ,j+1,k-1)+f(i  ,j+1,k  )
     &        + f(i+1,j  ,k-1)+f(i+1,j  ,k  )+f(i+1,j+1,k-1)+f(i+1,j+1,k  ))

        f7 = a*(f(i  ,j-1,k  )+f(i  ,j-1,k+1)+f(i  ,j  ,k  )+f(i  ,j  ,k+1)
     &        + f(i+1,j-1,k  )+f(i+1,j-1,k+1)+f(i+1,j  ,k  )+f(i+1,j  ,k+1))

        f8 = a*(f(i  ,j  ,k  )+f(i  ,j  ,k+1)+f(i  ,j+1,k  )+f(i  ,j+1,k+1)
     &        + f(i+1,j  ,k  )+f(i+1,j  ,k+1)+f(i+1,j+1,k  )+f(i+1,j+1,k+1))
 
c at each position the weight is based on the local dimension of the computational cell       
        wk1(i,j,k) =(f1/(au(i-1)*bv(j-1)*cw(k-1))*0.5*(rp(i)+rp(i-1))
     &            + f2/(au(i-1)*bv(j  )*cw(k-1))*0.5*(rp(i)+rp(i-1))
     &            + f3/(au(i-1)*bv(j-1)*cw(k  ))*0.5*(rp(i)+rp(i-1))
     &            + f4/(au(i-1)*bv(j  )*cw(k  ))*0.5*(rp(i)+rp(i-1))
     &            + f5/(au(i  )*bv(j-1)*cw(k-1))*0.5*(rp(i)+rp(i+1))
     &            + f6/(au(i  )*bv(j  )*cw(k-1))*0.5*(rp(i)+rp(i+1))
     &            + f7/(au(i  )*bv(j-1)*cw(k  ))*0.5*(rp(i)+rp(i+1))
     &            + f8/(au(i  )*bv(j  )*cw(k  ))*0.5*(rp(i)+rp(i+1)))*dV
      enddo
      enddo
      enddo

      fb = wk1
      
      return
      end
c-------------------------------------------------------------

c--------------------------------------------------------------
      subroutine filter2d(f,fb,wk1,wk2,im,jm,km)
c--------------------------------------------------------------
cc box filter routine (trapezodial rule for NON-UNIFORM grid)
cc On the boundaries backward and forward formulas are used
cc fb is not defined on the boundaries on output
cc 
c
      include 'common.h'
c
      integer im,jm,km
      real f(im,jm,km),fb(im,jm,km)
      real wk1(im,jm,km),wk2(im,jm,km)
      integer i,j,k

c--- initialize work arrays
      wk1=f

c--- filter in y and z ---
      fb(:,JY1:JY2,KZ1:KZ2) = 0.125*wk1(:,JY1-1:JY2-1,KZ1:KZ2)
     &       +           0.50*wk1(:,JY1  :JY2  ,KZ1:KZ2)
     &       +          0.125*wk1(:,JY1+1:JY2+1,KZ1:KZ2)
     &       +          0.125*wk1(:,JY1:JY2,KZ1-1:KZ2-1)
     &       +          0.125*wk1(:,JY1:JY2,KZ1+1:KZ2+1)
c
      return
      end
c-------------------------------------------------------------



c--------------------------------------------------------------
      subroutine weights
c
cc calculates the weights on each direction 
cc for a box filter on non-uniform grids
cc trapezodial rule is used
c
      include 'common.h'
c
c--------------------------------------------------------------------
c                   Weights for hat filter (2*DELTA)
c--------------------------------------------------------------------
c        
c--------------------------  x-direction  ---------------------------

      whx(ix1  :ix2  , 0) = 0.50
c
      whx(ix1+1:ix2-1,-1) = 0.25*av(ix1+1:ix2-1)/au(ix1+1:ix2-1)
      whx(ix1+1:ix2-1, 1) = 0.25*av(ix1+1:ix2-1)/au(ix1  :ix2-2)
c
c...backward and forward formulas for bc points
      whx(ix1,-1) =  0.50*av(ix1+1)/au(ix1  )+0.25*av(ix1+1)/au(ix1+1)
      whx(ix1, 1) = -0.25*av(ix1+1)/au(ix1  )
c     
      whx(ix2,-1) =  0.50*av(ix2-1)/au(ix2-1)+0.25*av(ix2-1)/au(ix2-2)
      whx(ix2, 1) = -0.25*av(ix2-1)/au(ix2-1)

c--------------------------  y-direction  ---------------------------

      why(jy1  :jy2  , 0) = 0.50
c
      why(jy1+1:jy2-1,-1) = 0.25*bw(jy1+1:jy2-1)/bv(jy1+1:jy2-1)
      why(jy1+1:jy2-1, 1) = 0.25*bw(jy1+1:jy2-1)/bv(jy1  :jy2-2)
c
c...backward and forward formulas for bc points

      why(jy1,-1) =  0.50*bw(jy1+1)/bv(jy1  )+0.25*bw(jy1+1)/bv(jy1+1)
      why(jy1, 1) = -0.25*bw(jy1+1)/bv(jy1  )
c
      why(jy2,-1) =  0.50*bw(jy2-1)/bv(jy2-1)+0.25*bw(jy2-1)/bv(jy2-2)
      why(jy2, 1) = -0.25*bw(jy2-1)/bv(jy2-1)

c--------------------------  z-direction  ---------------------------

      whz(kz1  :kz2  , 0) = 0.50
c
      whz(kz1  :kz2  ,-1) = 0.25*cu(kz1  :kz2  )/cw(kz1  :kz2  )
      whz(kz1  :kz2  , 1) = 0.25*cu(kz1  :kz2  )/cw(kz1-1:kz2-1)
c
c...backward and forward formulas for bc points
      if(itype(5)/=0) then
        whz(kz1,-1) =  0.50*cu(kz1+1)/cw(kz1  )+0.25*cu(kz1+1)/cw(kz1+1)
        whz(kz1, 1) = -0.25*cu(kz1+1)/cw(kz1  )
      endif
c
      if(itype(6)/=0) then
        whz(kz2,-1) =  0.50*cu(kz2-1)/cw(kz2-1)+0.25*cu(kz2-1)/cw(kz2-2)
        whz(kz2, 1) = -0.25*cu(kz2-1)/cw(kz2-1)
      endif
c 
      return
      end


c---- subroutine filterLaplace------------------------------------------
C
C     PURPOSE: Laplacian filter
C
c-----------------------------------------------------------------------
      subroutine filterLaplace(uc,nx,ny,nz,flag,nf)
c
c flag=0: uniform case
c flag=1: non uniform case, but the curvature of the cylindrical grid is not taken into account
c otherwise: non uniform case with weights based also on the curvature of the cylindrical grid
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,flag,nf
      real    uc(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k,ii
      real    rcyl
      real    uf(nx,ny,nz)
c
      rcyl = real(icyl)
      if(flag==1) rcyl=0.0

      !Second order accurate on non-uniform grids
      if(flag==0) then
c uniform case
        do ii=1,nf
          do k=2,nz-1
          do j=2,ny-1
          do i=2,nx-1
            uf(i,j,k) = uc(i-1,j,k)+uc(i+1,j,k)+uc(i,j-1,k)+uc(i,j+1,k)+uc(i,j,k-1)+uc(i,j,k+1)-6.0*uc(i,j,k)
          enddo
          enddo
          enddo

          call refreshbc(uf,nx*ny,nz)
          uc = uf
        enddo
      else
c non uniform case: the weight is dependent on the distance of each component from the (i,j,k) position
c smaller distances are associated with larger weights
c also the cylindrical grid case is taken into account
        do ii=1,nf
          do k=2,nz-1
          do j=2,ny-1
          do i=2,nx-1
            uf(i,j,k) =
     &                  av(i)/au( i )*(1.0-rcyl*0.5/(rp(i)*au(i-1)))*uc(i-1,j,k)
     &                + av(i)/au(i-1)*(1.0+rcyl*0.5/(rp(i)*au( i )))*uc(i+1,j,k)
     &                + uc(i,j-1,k) + uc(i,j+1,k)
     &                + cv(k)/cw( k )*uc(i,j,k-1) + cv(k)/cw(k-1)*uc(i,j,k+1)
     &                - 6.0*uc(i,j,k)
c     &                  uc(i,j-1,k) + uc(i,j+1,k)
c     &                + cv(k)/cw( k )*uc(i,j,k-1) + cv(k)/cw(k-1)*uc(i,j,k+1)
c     &                - 4.0*uc(i,j,k)
          enddo
          enddo
          enddo          

          call refreshbc(uf,nx*ny,nz)
          uc = uf
        enddo

      endif

      return

      end
c-----------------------------------------------------------------------

C---- subroutine structurefunction-----------N. Beratlis-Sep. 10 2011---
C
C     PURPOSE: Calculate second order velocity structure function.
C
c-----------------------------------------------------------------------
      subroutine structurefunction(uo,vo,wo,tv,uc,nx,ny,nz)
c
      include 'common.h'
c
c.... Input/Output arrays
      integer nx,ny,nz
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),tv(nx,ny,nz),uc(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k,Lf,nf
      real    pw,pw2,inv,c,Ck

      tv = 0.0
      pw = 2.0/3.0    ! exponent for the grid steps in each direction
c      pw2 = 3.0/2.0
      pw2 = 4.0/9.0   ! exponent for the product of the local dimensions
      inv = 1.0/6.0   ! coefficient of the FSF
      Ck = 1.4    ! Kolmogorov constant
c      Ck = 0.5
      c = 0.0014*(Ck**(-3.0/2.0))  ! constant of the FSF(3)
      Lf = 2 !Laplace flag   
! 0 uniform grid; 1 non uniform grid without curvature; otherwise non uniform grid with curvature
      nf = 3 !Filter iterations

c evaluation of the X velocity component at the centered positions
      call centervel(uo,uc,nx,ny,nz,1)
c values to the ghost layers
      call refreshbc(uc,nx*ny,nz)
c the X velocity field is filtered by a Laplacian filter taking into account the non-uniformity of
c the computational grid
      call filterLaplace(uc,nx,ny,nz,Lf,nf)

c component of the filtered structure function from the X velocity
      do k=2,nz-1
      do j=2,ny-1
      do i=2,nx-1
        tv(i,j,k) =
     &              ((uc(i,j,k)-uc(i-1,j,k))**2.0 + (uc(i,j,k)-uc(i+1,j,k))**2.0)*(ap(i)**pw)
     &            + ((uc(i,j,k)-uc(i,j-1,k))**2.0 + (uc(i,j,k)-uc(i,j+1,k))**2.0)*((bp(j)/rp(i))**pw)
     &            + ((uc(i,j,k)-uc(i,j,k-1))**2.0 + (uc(i,j,k)-uc(i,j,k+1))**2.0)*(cp(k)**pw)
      enddo
      enddo
      enddo

c evaluation of the Y velocity component at the centered positions
      call centervel(vo,uc,nx,ny,nz,2)
c values to the ghost layers
      call refreshbc(uc,nx*ny,nz)
c the Y velocity field is filtered
      call filterLaplace(uc,nx,ny,nz,Lf,nf)

c component of the filtered structure function from the Y velocity
      do k=2,nz-1
      do j=2,ny-1
      do i=2,nx-1
        tv(i,j,k) = tv(i,j,k)
     &            + ((uc(i,j,k)-uc(i-1,j,k))**2.0 + (uc(i,j,k)-uc(i+1,j,k))**2.0)*(ap(i)**pw)
     &            + ((uc(i,j,k)-uc(i,j-1,k))**2.0 + (uc(i,j,k)-uc(i,j+1,k))**2.0)*((bp(j)/rp(i))**pw)
     &            + ((uc(i,j,k)-uc(i,j,k-1))**2.0 + (uc(i,j,k)-uc(i,j,k+1))**2.0)*(cp(k)**pw)
      enddo
      enddo
      enddo

c evaluation of the Z velocity component at the centered positions
      call centervel(wo,uc,nx,ny,nz,3)
c values to the ghost layers
      call refreshbc(uc,nx*ny,nz)
c the Z velocity field is filtered
      call filterLaplace(uc,nx,ny,nz,Lf,nf)

c component of the filtered structure function from the Z velocity and final
c value of the same function
      do k=2,nz-1
      do j=2,ny-1
      do i=2,nx-1
        tv(i,j,k) = c*((rp(i)/(ap(i)*bp(j)*cp(k)))**pw2)*((tv(i,j,k)
     &            + ((uc(i,j,k)-uc(i-1,j,k))**2.0 + (uc(i,j,k)-uc(i+1,j,k))**2.0)*(ap(i)**pw)
     &            + ((uc(i,j,k)-uc(i,j-1,k))**2.0 + (uc(i,j,k)-uc(i,j+1,k))**2.0)*((bp(j)/rp(i))**pw)
     &            + ((uc(i,j,k)-uc(i,j,k-1))**2.0 + (uc(i,j,k)-uc(i,j,k+1))**2.0)*(cp(k)**pw) )*inv)**0.5
      enddo
      enddo
      enddo
c
c eddy viscosity to the ghost layers
      CALL REFRESHBC(TV,NX*NY,NZ)
c
c eddy viscosity to the boundaries of the computational domain
      CALL BOUNDTURB(TV,NX,NY,NZ,0)

      return

      end
c-----------------------------------------------------------------------


C---- SUBROUTINE WALE------------------------N. Beratlis-09 Nov. 2011---
C
C     PURPOSE: Compute the WALE eddy viscosity model
C
c-----------------------------------------------------------------------
      subroutine wale(uo,vo,wo,tv,nx,ny,nz)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),tv(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k
      real    dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
      real    s2,o2,ss,sdsd,IVso,Cwale,c1,c2,p1,p2,p3,Delta2,ppp
      real    s(3,3),o(3,3),a(3,3),b(3,3),error
c
c.... Functions
      real    tensor_double_inner_product,tensor_double_product

      error = 1.0e-10
      c1 = 1.0/6.0
      c2 = 2.0/3.0
      p1 = 3.0/2.0
      p2 = 5.0/2.0
      p3 = 5.0/4.0
      ppp = 2.0/3.0
      Cwale = 0.5**2.0

      do k=2,nz-1
      do j=2,ny-1
      do i=2,nx-1

        dudx = ap(i)*(uo(i,j,k)-uo(i-1,j,k))

        dudy = 0.5*( 0.5*bv(j)*(uo(i  ,j+1, k )-uo(i  ,j-1, k )) 
     &             + 0.5*bv(j)*(uo(i-1,j+1, k )-uo(i-1,j-1, k )) )/rp(i)

        dudz = 0.5*( 0.5*cu(k)*(uo(i  , j ,k+1)-uo(i  , j ,k-1))
     &             + 0.5*cu(k)*(uo(i-1, j ,k+1)-uo(i-1, j ,k-1)) )

        dvdx = 0.5*( 0.5*av(i)*(vo(i+1,j  , k )-vo(i-1,j  , k ))
     &             + 0.5*av(i)*(vo(i+1,j-1, k )-vo(i-1,j-1, k )) )


        dvdy = bp(j)/rp(i)*(vo(i,j+1,k)-vo(i,j,k))

        dvdz = 0.5*( 0.5*cv(k)*(vo( i ,j  ,k+1)-vo( i ,j  ,k-1))
     &             + 0.5*cv(k)*(vo( i ,j-1,k+1)-vo( i ,j-1,k-1)) )

        dwdx = 0.5*( 0.5*aw(i)*(wo(i+1, j ,k  )-wo(i-1, j ,k  ))
     &             + 0.5*aw(i)*(wo(i+1, j ,k-1)-wo(i-1, j ,k-1)) )

        dwdy = 0.5*( 0.5*bw(j)*(wo( i ,j+1,k  )-wo( i ,j-1,k  ))
     &             + 0.5*bw(j)*(wo( i ,j+1,k-1)-wo( i ,j-1,k-1)) )/rp(i)

        dwdz = cp(k)*(wo(i,j,k)-wo(i,j,k-1))

        s(1,1) = dudx
        s(1,2) = 0.5*(dudy + dvdx - icyl*0.5*(vo(i,j,k)+vo(i,j-1,k))/rp(i))
        s(1,3) = 0.5*(dudz + dwdx)
        s(2,1) = s(1,2)
        s(2,2) = dvdy+icyl*0.5*(uo(i,j,k)+uo(i-1,j,k))/rp(i)
        s(2,3) = 0.5*(dvdz + dwdy)
        s(3,1) = s(1,3)
        s(3,2) = s(2,3)
        s(3,3)  = dwdz

        o(1,1) = 0.0
        o(1,2) = 0.5*(dudy - dvdx - icyl*0.5*(vo(i,j,k)+vo(i,j-1,k))/rp(i))
        o(1,3) = 0.5*(dudz - dwdx)
        o(2,1) =-o(1,2)
        o(2,2) = 0.0
        o(2,3) = 0.5*(dvdz - dwdy)
        o(3,1) =-o(1,3)
        o(3,2) =-o(2,3)
        o(3,3) = 0.0

        s2 = tensor_double_product(s,s,3)
        o2 = tensor_double_product(o,o,3)
      
        call tensor_inner_product(s,s,a,3)
        call tensor_inner_product(o,o,b,3)

        IVso = tensor_double_inner_product(a,b,3)


        sdsd = c1*(s2*s2 + o2*o2) + c2*s2*o2  +2.0*IVso 


        if(sdsd.le.0) sdsd=0

        Delta2 = (RP(I)/(AP(I)*BP(J)*CP(K)))**PPP        
        tv(i,j,k) = Cwale*Delta2*(sdsd**p1)/(s2**p2 + sdsd**p3 + error)

      enddo
      enddo
      enddo

      CALL REFRESHBC(TV,NX*NY,NZ)
c
      CALL BOUNDTURB(TV,NX,NY,NZ,0)

      return

      end
c-----------------------------------------------------------------------


C---- subroutine tensor_inner_product-------N. Beratlis-09 Nov. 2011----
C
C     PUUPOSE: Compute inner product P of two tensors T and S.
C
c-----------------------------------------------------------------------
      subroutine tensor_inner_product(T,S,P,n)
c
c.... Input/Output arrays
      integer n
      real    T(n,n),S(n,n),P(n,n)
c
c.... Local arrays
      integer i,j,k

      P = 0.0

      do i=1,n
      do j=1,n
      do k=1,n
        P(i,j) = P(i,j) + T(i,k)*S(k,j)
      enddo
      enddo
      enddo

      return

      end
c-----------------------------------------------------------------------


C---- function tensor_double_inner_product--N. Beratlis-09 Nov. 2011----
C
C     PUUPOSE: Compute inner product of two 3x3 tensors
C
c-----------------------------------------------------------------------
      real function tensor_double_inner_product(T,S,n)
c
c.... Input/Output arrays
      integer n
      real    T(n,n),S(n,n)
c
c.... Local arrays
      integer i,j
      real    P

      P = 0.0

      do i=1,n
      do j=1,n
        P = P + T(i,j)*S(j,i)
      enddo
      enddo

      tensor_double_inner_product = P

      return 

      end
c-----------------------------------------------------------------------


C---- function tensor_double_product--N. Beratlis-09 Nov. 2011----
C
C     PUUPOSE: Compute inner product of two 3x3 tensors
C
c-----------------------------------------------------------------------
      real function tensor_double_product(T,S,n)
c
c.... Input/Output arrays
      integer n
      real    T(n,n),S(n,n)
c
c.... Local arrays
      integer i,j
      real    P

      P = 0.0

      do i=1,n
      do j=1,n
        P = P + T(i,j)*S(i,j)
      enddo
      enddo

      tensor_double_product = P

      return 

      end
c-----------------------------------------------------------------------


C---- subroutine conv_derv_scalar------------N. Beratlis-01 Dec. 2011--
C
C     PURPOSE: Calculate convective derivative of a scalar quantity.
C
c-----------------------------------------------------------------------
      subroutine conv_derv_scalar(der,scal,uc,vc,wc,nx,ny,nz)
c
      include 'common.h'
c
c.... Input/Output arrays
      integer nx,ny,nz
      real    uc(nx,ny,nz),vc(nx,ny,nz),wc(nx,ny,nz),der(nx,ny,nz),scal(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k
      real    dsdxp,dsdxm,dsdx,dsdyp,dsdym,dsdy,dsdzp,dsdzm,dsdz

      DER = 0.0

      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
        DSDXP = AU(I  )*(SCAL(I+1,J,K)-SCAL(I  ,J,K))
        DSDXM = AU(I-1)*(SCAL(I  ,J,K)-SCAL(I-1,J,K))
        DSDX = 0.5*(DSDXP+DSDXM)
        DSDYP = BV(J  )*(SCAL(I,J+1,K)-SCAL(I,J  ,K))/RP(I)
        DSDYM = BV(J-1)*(SCAL(I,J  ,K)-SCAL(I,J-1,K))/RP(I)
        DSDY = 0.5*(DSDYP+DSDYM)
        DSDZP = CW(K  )*(SCAL(I,J,K+1)-SCAL(I,J,K  ))
        DSDZM = CW(K-1)*(SCAL(I,J,K  )-SCAL(I,J,K-1))
        DSDZ = 0.5*(DSDZP+DSDZM)
        DER(I,J,K) = UC(I,J,K)*DSDX + VC(I,J,K)*DSDY + WC(I,J,K)*DSDZ
      ENDDO
      ENDDO
      ENDDO
      
      RETURN

      END
c-----------------------------------------------------------------------
