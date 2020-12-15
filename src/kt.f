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
c-----SUBROUTINE-Turvis------------------------A.Pal 7/22/15-------
c
c
      SUBROUTINE KURVIS(UO,VO,WO,DENS,KV,G,GB,SXX,SYY,SZZ,SXY,SYZ,SXZ,
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
      REAL    UO(IM,JM,KM),VO(IM,JM,KM),WO(IM,JM,KM),DENS(IM,JM,KM)
      REAL    KV(IM,JM,KM),G(IM,JM,KM),GB(IM,JM,KM)
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
!        KV(I,:,:) = CSS*CSS*ALPHA*DXDYDZ(I,:,:)*G(I,:,:)
!        CLES(I,:,:) = CSS*CSS*ALPHA
!        ENDDO

        ALPHA = 1.0

        IF(LP==0 .AND. NP==0) THEN

          UTAU = 0.5*SQRT(-DPDZ)

          IF(ICYL==1) THEN

            DO I=IX1,IX2

              RPLUS = UTAU*ABS(0.5-RP(I))/RU1
              ALPHA = (1.-EXP(-RPLUS/25.))**2

              KV(I,:,:) = CSS*CSS*ALPHA*DXDYDZ(I,:,:)*G(I,:,:)
              CLES(I,:,:) = CSS*CSS*ALPHA
            
            ENDDO

          ELSEIF(ICYL==0) THEN

            DO I=IX1,IX2
              DO J=JY1,JY2
c 
                RPLUS = UTAU*ABS(0.5-SQRT(XC(I)**2.+YC(J)**2.))/RU1
                ALPHA = (1.-EXP(-RPLUS/25.))**2

                KV(I,J,:) = CSS*CSS*ALPHA*DXDYDZ(I,J,:)*G(I,J,:)
                CLES(I,J,:) = CSS*CSS*ALPHA

              ENDDO
            ENDDO
 
          ENDIF

        ELSE

          KV(:,:,:) = CSS*CSS*ALPHA*DXDYDZ(:,:,:)*G(:,:,:)
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
        CALL FILTER3D(MXX,MXX,LM,KV,IM,JM,KM) 
        CALL FILTER3D(MYY,MYY,LM,KV,IM,JM,KM) 
        CALL FILTER3D(MZZ,MZZ,LM,KV,IM,JM,KM) 
        CALL FILTER3D(MXY,MXY,LM,KV,IM,JM,KM) 
        CALL FILTER3D(MYZ,MYZ,LM,KV,IM,JM,KM) 
        CALL FILTER3D(MXZ,MXZ,LM,KV,IM,JM,KM) 

c.....\hat{Sxx} --> sxx; \hat{Syy} --> syy; \hat{Szz} --> szz; 
c.....\hat{Sxy} --> sxy; \hat{Syz} --> syz; \hat{Sxz} --> sxz
        CALL FILTER3D(SXX,SXX,LM,KV,IM,JM,KM)
        CALL FILTER3D(SYY,SYY,LM,KV,IM,JM,KM)
        CALL FILTER3D(SZZ,SZZ,LM,KV,IM,JM,KM)
        CALL FILTER3D(SXY,SXY,LM,KV,IM,JM,KM)
        CALL FILTER3D(SYZ,SYZ,LM,KV,IM,JM,KM)
        CALL FILTER3D(SXZ,SXZ,LM,KV,IM,JM,KM)

c.....\hat{|S|} --> gb
        CALL FILTER3D(G,GB,LM,KV,IM,JM,KM) 

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
        CALL FILTER3D(UC,SXX,LM,KV,IM,JM,KM)
        CALL FILTER3D(VC,SYY,LM,KV,IM,JM,KM)
        CALL FILTER3D(WC,SZZ,LM,KV,IM,JM,KM)

c.....uu --> gb; \hat{uu} --> gb; 
c.....(\hat{uu}-\hat(u)\hat(u))*mxx --> mxx
        GB = UC*UC
        CALL FILTER3D(GB,GB,LM,KV,IM,JM,KM)
        MXX = (GB-SXX*SXX)*MXX

c.....vv --> gb; \hat{vv} --> gb; 
c.....(\hat{vv}-\hat(v)\hat(v))*myy --> myy
        GB = VC*VC
        CALL FILTER3D(GB,GB,LM,KV,IM,JM,KM)
        MYY = (GB-SYY*SYY)*MYY

c.....ww --> gb; \hat{ww} --> gb; 
c.....(\hat{ww}-\hat(w)\hat(w))*mzz --> mzz
        GB = WC*WC
        CALL FILTER3D(GB,GB,LM,KV,IM,JM,KM)
        MZZ = (GB-SZZ*SZZ)*MZZ

c.....uv --> gb; \hat{uv} --> gb; 
c.....(\hat{uv}-\hat(u)\hat(v))*mxy --> mxy
        GB = UC*VC
        CALL FILTER3D(GB,GB,LM,KV,IM,JM,KM)
        MXY = (GB-SXX*SYY)*MXY

c.....vw --> gb; \hat{vw} --> gb; 
c.....(\hat{vw}-\hat(v)\hat(w))*myz --> myz
        GB = VC*WC
        CALL FILTER3D(GB,GB,LM,KV,IM,JM,KM)
        MYZ=(GB-SYY*SZZ)*MYZ

c.....uw --> gb; \hat{uw} --> gb; 
c.....(\hat{uw}-\hat(u)\hat(w))*mxz --> mxz
        GB = UC*WC
        CALL FILTER3D(GB,GB,LM,KV,IM,JM,KM)
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

          KV = ILM*IMM
c          KV = -6.0*ILM*IMM
          GB = 0.
          WHERE(KV>EPS) GB = 1./(1.+1.5*SQRT(DXDYDZ)*KV**(-0.125)/DT)

c the new values of ILM and IMM are estimated
          ILM = GB*LM+(1.-GB)*MXX
          IMM = GB*MM+(1.-GB)*MYY

c          WHERE(ILM<0) ILM = 0.

          ELSE
c...The convective derivatives of ILM and IMM are evaluated and stored in MXX and MYY
            CALL CONV_DERV_SCALAR(MXX,ILM,UC,VC,WC,IM,JM,KM)
            CALL CONV_DERV_SCALAR(MYY,IMM,UC,VC,WC,IM,JM,KM)
            KV = ILM*IMM

            GB = 0.
            WHERE(KV>EPS) GB = 1./(1.+1.5*SQRT(DXDYDZ)*KV**(-0.125)/DT)

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

        KV = 0.
        CLES = 0.
c the value of the eddy viscosity is estimated
c if KV is negative it is not physical and is enforced equal to 0
        WHERE(IMM>SQRT(EPS)) KV = (ILM/IMM)*DXDYDZ*G
c        WHERE(IMM>SQRT(EPS)) KV = -0.5*(ILM/IMM)*DXDYDZ*G
        WHERE(IMM>SQRT(EPS)) CLES = (ILM/IMM)
c        WHERE(IMM>SQRT(EPS)) CLES = -0.5*(ILM/IMM)
c        WHERE(KV<-1.0*RU1) KV = -1.0*RU1
        WHERE(KV<0.0) KV = 0.0

!!!!!!        KV = 0.0
C
      ENDIF
C
c...eddy viscosity and LES coefficient at the ghost layers are updated
      CALL REFRESHBC(KV,IM*JM,KM)
      CALL REFRESHBC(CLES,IM*JM,KM)
c
c...eddy viscosity and LES coefficient at the boundaries of the computational domain are updated
      CALL BOUNDTURB(KV,IM,JM,KM,0)
      CALL BOUNDTURB(CLES,IM,JM,KM,0)
c
      RETURN
      END