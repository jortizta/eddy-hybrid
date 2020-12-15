	SUBROUTINE WRITE_PLANES(NX,NY,NZ,NZG,NBD,ICYCLE,TIME,DTM1,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,P,VO,UO,WO,DENS)
        
        use vorticity
        use density_bg
	INCLUDE 'common.h'
        INCLUDE 'mpif.h'
	INTEGER :: NX,NY,NZ,NZG,NBD,ICYCLE
	REAL	:: DTM1,TIME
	REAL    :: XU(NX),YV(NY),ZW(NZ),ZWG(NZG)
	REAL    :: XC(NX),YC(NY),ZC(NZ),ZCG(NZG)
	REAL    :: UO(NX,NY,NZ),VO(NX,NY,NZ),WO(NX,NY,NZ),P(NX,NY,NZ),DENS(NX,NY,NZ)
	REAL	:: XU_CAR(NX,NY),YU_CAR(NX,NY),XV_CAR(NX,NY),YV_CAR(NX,NY),XC_CAR(NX,NY),YC_CAR(NX,NY)
	INTEGER :: I,J,K,STAT
        integer iic,ic,ip,im,jjc,jc,jp,jm,jsy,jpsy,kkc,kc,kp,km
        real dtheta,dr,dz,dq1x3,dq3x1,dq2x3,dq3x2,dq2x1,dq1x2
	REAL,ALLOCATABLE,DIMENSION(:,:,:)::U1_TMP2,U2_TMP2,U3_TMP2,rtmp1
	
!	if(mod(ICYCLE,10).eq.0) then
	ALLOCATE(U1_TMP2(NX,NY,NZ),U2_TMP2(NX,NY,NZ),U3_TMP2(NX,NY,NZ),rtmp1(nx,ny,nz))
 	CALL CENTER_VELOCITY(NX,NY,NZ,UO,U1_TMP2,1)
 	CALL CENTER_VELOCITY(NX,NY,NZ,VO,U2_TMP2,2)
 	CALL CENTER_VELOCITY(NX,NY,NZ,WO,U3_TMP2,3)
! 	WRITE(*,*), UO(:,:,:)
! 	IF(MYRANK.EQ.0) CALL GATHER_2D_SLAVES
! 	IF(MYRANK.EQ.0) THEN
!  	WRITE(*,*),P(:,5,1),P(:,5,KZ2),P(5,1 ,:),P(5,JY2,:)
! 	WRITE(*,*),WO(NX/2,NY/2,NZ-2),WO(NX/2,NY/2,NZ-1),0.50*(WO(NX/2,NY/2,NZ-1)+WO(NX/2,NY/2,NZ-2)),U3_TMP2(NX/2,NY/2,NZ-1)
! 	ELSE
! 	WRITE(*,*),NZ,KZ1,KZ2,MYRANK
! 	WRITE(*,*),WO(NX/2,NY/2,1),WO(NX/2,NY/2,2),0.50*(WO(NX/2,NY/2,NZ-1)+WO(NX/2,NY/2,NZ)),U3_TMP2(NX/2,NY/2,NZ-1)
! 	ENDIF
! 	CALL WRITE_PLANE(UO,2,NY/2,0,0,'UO_',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,ICYCLE,TIME,DTM1,STAT)
! 	CALL WRITE_PLANE(UO,2,2,0,0,'UO_',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,ICYCLE,TIME,DTM1,STAT)


!----------------------------------------------------------------------------------------------------------------



        CALL WRITE_PLANE(U1_TMP2,2,NY/4+1,0,0,'U0',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)
        CALL WRITE_PLANE(U2_TMP2,2,NY/4+1,0,0,'V0',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)
        CALL WRITE_PLANE(U3_TMP2,2,NY/4+1,0,0,'W0',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)





 	
        CALL WRITE_PLANE(U1_TMP2,2,NY/2+NY/4+1,0,0,'U0',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)
        CALL WRITE_PLANE(U2_TMP2,2,NY/2+NY/4+1,0,0,'V0',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)
        CALL WRITE_PLANE(U3_TMP2,2,NY/2+NY/4+1,0,0,'W0',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)






        CALL WRITE_PLANE(U1_TMP2,2,2,0,0,'U0',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)
        CALL WRITE_PLANE(U2_TMP2,2,2,0,0,'V0',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)
        CALL WRITE_PLANE(U3_TMP2,2,2,0,0,'W0',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)




        CALL WRITE_PLANE(U1_TMP2,2,NY/2+1,0,0,'U0',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)
        CALL WRITE_PLANE(U2_TMP2,2,NY/2+1,0,0,'V0',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)
        CALL WRITE_PLANE(U3_TMP2,2,NY/2+1,0,0,'W0',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)
!----------------------------------------------------------------------------------------------------------------



        do k=kz1,kz2
        do j=jy1,jy2
        do i=ix1,ix2
            dz=zcg(k)-zcg(k-1)
            dr=xc(i)-xc(i-1)
            dq1x3=(0.25*(uo(i,j,k)+uo(i,j,k+1)+uo(i-1,j,k)+uo(i-1,j,k+1))-0.25*(uo(i,j,k)+uo(i,j,k-1)+uo(i-1,j,k)+uo(i-1,j,k-1)))/dz
            dq3x1=(0.25*(wo(i,j,k)+wo(i,j,k-1)+wo(i+1,j,k)+wo(i+1,j,k-1))-0.25*(wo(i,j,k)+wo(i,j,k-1)+wo(i-1,j,k)+wo(i-1,j,k-1)))/dr
            rtmp1(i,j,k)=dq1x3-dq3x1
        enddo
        enddo
        enddo

        CALL WRITE_PLANE(rtmp1,2,NY/4+1,0,0,'omgt',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)


        CALL WRITE_PLANE(rtmp1,2,NY/2+NY/4+1,0,0,'omgt',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)


        CALL WRITE_PLANE(rtmp1,2,2,0,0,'omgt',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)


        CALL WRITE_PLANE(rtmp1,2,NY/2+1,0,0,'omgt',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)


        do k=kz1,kz2
        do j=jy1,jy2
        do i=ix1,ix2
            dr=xc(i)-xc(i-1)
            dq1x3=(uo(i,j,k)-uo(i-1,j,k))/dr
            rtmp1(i,j,k)=dq1x3
        enddo
        enddo
        enddo

        CALL WRITE_PLANE(rtmp1,2,NY/4+1,0,0,'DURDR',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &                          YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)
 
        CALL WRITE_PLANE(rtmp1,2,NY/2+NY/4+1,0,0,'DURDR',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &                          YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

        CALL WRITE_PLANE(rtmp1,2,2,0,0,'DURDR',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &                          YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

        CALL WRITE_PLANE(rtmp1,2,NY/2+1,0,0,'DURDR',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &                          YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)


!         call DERIVATIVE(uo,rtmp1,NX,NY,NZ,NZG,XC,ZCG,1,1)
!         CALL WRITE_PLANE(rtmp1,2,33,0,0,'DURDRD',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
!     &                          YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)


        do k=kz1,kz2
        do j=jy1,jy2
        do i=ix1,ix2
        dz=zcg(k)-zcg(k-1)
        dtheta=dely
        dq2x3=(0.25*(vo(i,j,k)+vo(i,j-1,k)+vo(i,j,k+1)+vo(i,j-1,k+1))-0.25*(vo(i,j,k)+vo(i,j-1,k)+vo(i,j,k-1)+vo(i,j-1,k-1)))/dz
        dq3x2=(0.25*(wo(i,j,k)+wo(i,j+1,k)+wo(i,j,k-1)+wo(i,j+1,k-1))-0.25*(wo(i,j,k)+wo(i,j-1,k)+wo(i,j,k-1)+wo(i,j-1,k-1)))/dtheta
        rtmp1(i,j,k)=dq3x2/rp(i)-dq2x3
        enddo
        enddo
        enddo

        CALL WRITE_PLANE(rtmp1,2,NY/4+1,0,0,'omgr',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

 
	CALL WRITE_PLANE(rtmp1,2,NY/2+NY/4+1,0,0,'omgr',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

        CALL WRITE_PLANE(rtmp1,2,2,0,0,'omgr',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

        CALL WRITE_PLANE(rtmp1,2,NY/2+1,0,0,'omgr',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

!
        do k=kz1,kz2
        do j=jy1,jy2
        do i=ix1,ix2
        dr=xc(i)-xc(i-1)
        dtheta=dely
        dq2x1=(0.25*((vo(i,j,k)+vo(i,j-1,k))*rp(i)+(vo(i+1,j,k)+vo(i+1,j-1,k))*rp(i+1))
     &  -0.25*((vo(i,j,k)+vo(i,j-1,k))*rp(i)+(vo(i-1,j,k)+vo(i-1,j-1,k))*rp(i-1)))/dz
        dq1x2=(0.25*(uo(i,j,k)+uo(i,j+1,k)+uo(i-1,j,k)+uo(i-1,j+1,k))-0.25*(uo(i,j,k)+uo(i,j-1,k)+uo(i-1,j,k)+uo(i-1,j-1,k)))/dtheta
        rtmp1(i,j,k)=(dq2x1-dq1x2)/rp(i)
        enddo
        enddo
        enddo

        CALL WRITE_PLANE(rtmp1,2,NY/4+1,0,0,'omgz',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

 
	CALL WRITE_PLANE(rtmp1,2,NY/2+NY/4+1,0,0,'omgz',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

	CALL WRITE_PLANE(rtmp1,2,2,0,0,'omgz',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

	CALL WRITE_PLANE(rtmp1,2,NY/2+1,0,0,'omgz',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)


!	CALL WRITE_PLANE(U3_TMP2,2,192,0,0,'omgaz1',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
!     &			YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)
	
!        if(MYRANK==0) write(*,*)"********WHAT IS WRONG    1***********"
!	CALL WRITE_PLANE(U3_TMP2,2,96,0,0,'omgaz1',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
!     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)



c$$$!        if(MYRANK==0) write(*,*)"********WHAT IS WRONG    2***********"

        do k=kz1,kz2
        do j=jy1,jy2
        do i=ix1,ix2
	U3_TMP2(I,J,K)=DENS(I,J,K)-DENS_BG(I,J)
	ENDDO
	ENDDO
	ENDDO
!------------------------------------------------------------------------------------------------------	
C$$$        CALL WRITE_PLANE(DENS,2,2,0,0,'RHO',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
C$$$     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

C$$$        CALL WRITE_PLANE(DENS,2,NY/2,0,0,'RHO',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
C$$$     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

C$$$        CALL WRITE_PLANE(DENS,2,NY/4,0,0,'RHO',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
C$$$     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

C$$$        CALL WRITE_PLANE(DENS,2,NY/2+NY/4,0,0,'RHO',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
C$$$     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

C$$$	CALL WRITE_PLANE(U3_TMP2,2,NY/4,0,0,'RHOP',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
C$$$     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

C$$$	CALL WRITE_PLANE(U3_TMP2,2,NY/2+NY/4,0,0,'RHOP',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
C$$$     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)
!------------------------------------------------------------------------------------------------------	



!        CALL WRITE_PLANE(U3_TMP2,3,768,0,0,'RHOP',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,XC_CAR,
!     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,ICYCLE,TIME,DTM1,STAT)

! 	CALL WRITE_PLANE(U3_TMP2,1,125,0,0,'WO',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,ICYCLE,TIME,DTM1,STAT)
! 
! 	CALL WRITE_PLANE(U3_TMP2,1,122,0,0,'WO',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,ICYCLE,TIME,DTM1,STAT)
! 
! 	CALL WRITE_PLANE(U3_TMP2,1,115,0,0,'WO',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,ICYCLE,TIME,DTM1,STAT)
! 
! 	CALL WRITE_PLANE(P,2,NY/2,0,0,'P',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,ICYCLE,TIME,DTM1,STAT)
! 	write(*,*),U3_TMP2(2,2,nz),U3_TMP2(2,2,nz-1),U3_TMP2(2,2,1),U3_TMP2(2,2,2)
! 	CALL WRITE_PLANE(WO,2,2,0,0,'WO_',.FALSE.,NX,NY,NZ,NZG,XC,XU,YC,YV,ZCG,ZWG,ICYCLE,TIME,DTM1,STAT)
!	endif
!------------------------------------------------------------------------------------------------------	

	END


       subroutine center_velocity(nx,ny,nz,Uin,Ucen,dir)
	INCLUDE 'common.h'
        INCLUDE 'mpif.h'
!Passed Variables
 	integer,intent(in)      :: dir,nx,ny,nz
 	real,intent(in)         :: Uin(nx,ny,nz) 
 	real,intent(out)        :: Ucen(nx,ny,nz) 
 	

!Local Variables
 	integer              	:: i,j,k,err
!Zero Output Array
 	Ucen=0.d0
!*************************************************
!********************X1***************************
!*************************************************
	if(dir.EQ.1) then
  !U
    	do k=2,nz-1
    	do j=2,ny-1
      	do i=2,nx-1
	Ucen(i,j,k)=0.50d0*(uin(i,j,k)+uin(i-1,j,k))
      	enddo
     	enddo
    	enddo
	
!*************************************************
!********************X2***************************
!*************************************************
	elseif (dir.EQ.2) then 
 !V
  !U
    	do k=2,nz-1
    	do j=2,ny-1
      	do i=2,nx-1
    	Ucen(i,j,k)=0.50d0*(uin(i,j,k)+uin(i,j-1,k))
    	enddo
   	enddo
 	enddo
!*************************************************
!********************X3***************************
!*************************************************
	elseif (dir.EQ.3) then 
 !W
  !U
    	do k=2,nz-1
    	do j=2,ny-1
      	do i=2,nx-1
    	Ucen(i,j,k)=0.50d0*(uin(i,j,k)+uin(i,j,k-1))
    	enddo
   	enddo
  	enddo
	
 	else
 !Invalid direction
  	write(*,'(a60,i2)') "INVALID DIRECTION IN center_velocities 
     &                      dir must be 1,2,3.  dir= ", dir
 	endif

	return
	end subroutine center_velocity


        subroutine write_plane(var,dir,index1,myidM,prec,varname,verbose,nx,ny,nz,nzg,xc,xu,yc,yv,zcg,zwg,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,icycle,TIME,dtm,stat)
	INCLUDE 'common.h'
        INCLUDE 'mpif.h'


!Passed Variables
	integer,intent(in)    :: dir,index1,myidM,prec,nx,ny,nz,nzg,icycle
 	real,intent(in)       :: var(nx,ny,nz),TIME,dtm
 	real,intent(in)       :: xc(nx),xu(nx),yc(ny),yv(ny),zcg(nzg),zwg(nzg)
 	integer,intent(out)   :: stat
 	character(len=*)      :: varname
 	logical,intent(in)    :: verbose

!Local Variables
 	integer :: Msub, MsubG
        integer,parameter      :: realtype  =MPI_DOUBLE_PRECISION
        integer,parameter      :: inttype   =MPI_INTEGER
	integer :: status1(MPI_STATUS_SIZE),index2
 	integer :: err1, s1,Rcoords(3),myidM2
 	character(len=300) :: filename
 	real(4),allocatable,dimension(:,:) :: SP_plane,TEMP
 	logical             :: cfile
 	integer             :: iu, iv, iw,N,Tsize,i,j,k,ks,kstart,sp1,sp2,is,js,istart,jstart
   	real,allocatable,dimension(:,:,:) :: PlnX1
   	real,allocatable,dimension(:,:,:) :: PlnX2
   	real,allocatable,dimension(:,:,:) :: PlnX3
	real,allocatable,dimension(:,:) :: Outplane
	real,allocatable,dimension(:,:) :: Temp_Recv, Temp_Send
	integer                            :: Xmaster,Xnode
	REAL	:: XU_CAR(NX,NY),YU_CAR(NX,NY),XV_CAR(NX,NY),YV_CAR(NX,NY),XC_CAR(NX,NY),YC_CAR(NX,NY)


	allocate(PlnX1(1:ny,1:nzg,2),PlnX2(1:nx,1:nzg,1:2),PlnX3(1:nx,1:ny,2))
	iu=0
 	iv=0
 	iw=0
	select case(varname)
 	case('uo')
 	iu=1 
 	case('vo')
 	iv=1
 	case('wo')
 	iw=1
	end select
	
 	if (MYRANK.eq.myidM) then
	if (dir.EQ.1) write(filename,'(a,i4.4,a,i8.8,a)') trim(varname)//"_i",index1,"_n",icycle,".pln"
  	if (dir.EQ.2) write(filename,'(a,i4.4,a,i8.8,a)') trim(varname)//"_j",index1,"_n",icycle,".pln"
!   	if (dir.EQ.3) write(filename,'(a,i4.4,a,i5.5,a)') trim(varname)//"_k",index1,"_n",icycle,".pln"

  	open(210,file=filename,form='unformatted',status='unknown',iostat=s1)
  	endif
	

	IF(DIR.EQ.1) THEN
	
	allocate(Temp_Recv(ny,nz))
	Temp_Recv(:,:)=var(index1,:,:)
	allocate(Outplane(1:ny,1:nzg))
	if(myrank.eq.0) then

	DO n=0,sizex2x3-1
	  if(n.eq.0) then
	    do j=2,ny-1
	     do k=2,nz-1
	      Outplane(j,k)=Temp_Recv(j,k)
	     enddo
	    enddo
	  else
 	      call MPI_RECV(Temp_Recv,1*ny*nz,MPI_DOUBLE_PRECISION,n,1,commx2x3,status1,stat)
	  endif
	 if(n.ne.0) then
	      kstart=n*(nzg-2)/nzprocs
! 	      istart=(nx-2)/nxprocs
!  	      write(*,*), kstart
	 endif
	 kstart=0
	      do j=2,ny-1
	       do k=2,nz-1
	        ks=kstart+k
! 		is=istart+i
 	        Outplane(j,ks)=Temp_Recv(j,k)
	       enddo
	      enddo
	ENDDO
	PlnX1(:,:,1)=Outplane(:,:)
	else
	call MPI_SEND(Temp_Recv,1*ny*nz,MPI_DOUBLE_PRECISION,0,1,commx2x3,status1,stat)

	endif
	
    	if (MYRANK.EQ.myidM) then
	 write(210) icycle,TIME,dtm,9.810d0,1.0d0,180.0,1.0
    	 write(210) dir,index1, iu, iv, iw
    	write(210) xc(index1), xu(index1)
    	write(210) ny, nzg
    	write(210) yc, yv
        write(210) zcg, zwg
    	if (prec.EQ.0) then
     	allocate(SP_plane(1:ny,1:nzg))
     	SP_plane=PlnX1(:,:,1) !store plane in single precsion
     	write(210)  SP_plane
     	deallocate(SP_plane)
    	else 
    	write(210) PlnX1(:,:,1)
    	endif
   	endif

	ELSEIF(DIR.EQ.2) THEN
	
	allocate(Temp_Recv(nx,nz))
	Temp_Recv(:,:)=var(:,index1,:)
	sp1=size(PlnX2,1)
	sp2=size(PlnX2,2)
	allocate(Outplane(sp1,sp2))
	if(myrank.eq.0) then

	DO n=0,sizex1x3-1
	  if(n.eq.0) then
	    do i=2,nx-1
	     do k=2,nz-1
	      Outplane(i,k)=Temp_Recv(i,k)
	     enddo
	    enddo

	  else
	      call MPI_RECV(Temp_Recv,nx*1*nz,MPI_DOUBLE_PRECISION,n,1,commx1x3,status1,stat)
	  endif
!	 if(n.ne.0) then
! 	      write(*,*), coords(1),coords(2),coords(3)
	      kstart=n*(nzg-2)/nzprocs
! 	      istart=(nx-2)/nxprocs
! 	      write(*,*), kstart
!	 endif
	 do i=2,nx-1
	    do k=2,nz-1
	       ks=kstart+k
!              is=istart+i
	       Outplane(i,ks)=Temp_Recv(i,k)
	    enddo
	 enddo
	ENDDO
	PlnX2(:,:,1)=Outplane(:,:)
	else
	    call MPI_SEND(Temp_Recv,nx*1*nz,MPI_DOUBLE_PRECISION,0,1,commx1x3,status1,stat)

	endif

	
 	if (MYRANK.EQ.myidM) then 
    	write(210) icycle,TIME,dtm,9.810d0,1.0d0,180.0,1.0
    	write(210) dir,index1, iu, iv, iw
    	write(210) yc(index1), yv(index1)
    	write(210) nx, nzg
    	write(210) xc, xu
    	write(210) zcg, zwg
    	if (prec.EQ.0) then
     	allocate( SP_plane(1:nx,1:nzg),STAT=s1 )
     	SP_plane=PlnX2(:,:,1) !store plane in single precsion
     	write(210)  SP_plane
     	deallocate(SP_plane)
    	else 
     	write(210) PlnX2(:,:,1)
    	endif
   	endif

	ELSEIF(DIR.EQ.3) THEN
!LET US FIND OUT WHICH PROCESSOR HAS THE REQUIRED INDEX1
! 	WRITE(*,*)MYRANK,NZ,NZG
! 	WRITE(*,*)"KLMIN=",sz,"KGMIN=",1,"KLMAX=",ez,"KGMAX=",NZG
	if(index1.gt.sz.and.index1.lt.ez.or.index1.eq.sz.or.index1.eq.ez) then
!	WRITE(*,*)"**FOUND THE PROCESSOR***", MYRANK
	Xmaster=1
	endif
	

 	if(Xmaster.eq.1) then 
	index2=index1-(sz-2)
!	write(*,*) index2
	myidM2=MYRANK
	allocate(Temp_Recv(nx,ny))
	Temp_Recv(:,:)=var(:,:,index2)
! 	call MPI_SEND(Temp_Recv,nx*ny*1,MPI_DOUBLE_PRECISION,0,1,commx1x2,status1,stat)
! 	endif
	

! 	if(MYRANK.eq.0) then
! 	allocate(Temp_Recv(nx,ny))
! 	Temp_Recv=0.0
! 	WRITE(*,*)"*********RECEIVE ERROR********"
! 	call MPI_RECV(Temp_Recv,nx*ny*1,MPI_DOUBLE_PRECISION,0,1,commx1x2,status1,stat)
! 	WRITE(*,*)"*********RECEIVE ERROR2********"
! 
! 	call MPI_BARRIER(MPI_COMM_EDDY,stat)
! 	
! 	elseif(MYRANK.eq.myidM2) then
! 	allocate(Temp_Recv(nx,ny))
! 	Temp_Recv(:,:)=var(:,:,nz-1)
! 	call MPI_SEND(Temp_Recv,nx*ny*1,MPI_DOUBLE_PRECISION,0,1,commx1x2,status1,stat)
! 	endif




! 	Xmaster = -1
!   	Xnode   = -1
!   	if (index1.GT.sz-1.AND.index1.LT.ez+1) then !1 
! !This node contains this plane
!    	if (rankx1x2.EQ.0) then 
!     	Xmaster = 1
! !     	Master  = myid
!    	else 
!     	Xnode   = 1
!    	endif 
! 
!   	elseif ( (index1.EQ.sz-1.AND.sz-1.EQ.1).OR.(index1.EQ.ez+1.AND.ez+1.EQ.NZG) )  then !1 
! !This node contains this plane and its a boundary plane
! 
!    	if (rankx1x2.EQ.0) then 
!     	Xmaster = 1
! !     	Master  = myid
!    	else 
!     	Xnode   = 1
!    	endif 
! 
! 	endif
! 	write(*,*)Xmaster, Xnode
! 	if(Xmaster.eq.1) then
! 	allocate(Temp_Recv(nx,ny))
! 	allocate(Temp_Send(nx,ny))
! 	do j=1,ny
!    	do i=1,nx
!     	Temp_Send(i,j) = var(i,j,index1)
!    	enddo
!   	enddo
! 	write(*,*)"****INDIR31*****"
! 	sp1=size(PlnX3,1)
! 	sp2=size(PlnX3,2)
! 	allocate(Outplane(1:sp1,1:sp2))
! 	
! ! 	if(myrank.eq.0) then
! 
! 	DO n=0,sizex1x2-1
! 	write(*,*)"****INDIR32*****"
! 	  if(n.eq.0) then
! 	    do i=2,nx-1
! 	     do j=2,ny-1
! 	      Outplane(i,j)=Temp_Recv(i,j)
! 	     enddo
! 	    enddo
! 	    Temp_Recv=Temp_Send
! 	  else
!  	      call MPI_RECV(Temp_Recv,nx*ny*1,MPI_DOUBLE_PRECISION,n,1,commx1x2,status1,stat)
! 	  endif
! 	 
! 	 if(n.ne.0) then
! 	      jstart=n*(ny-2)/nyprocs
!  	      istart=n*(nx-2)/nxprocs
! 	 endif
! 	      do i=2,nx-1
! 	       do j=2,ny-1
! 	        js=jstart+j
! 		is=istart+i
!  	        Outplane(is,js)=Temp_Recv(i,j)
! 	       enddo
! 	      enddo
! 	ENDDO
! 	PlnX3(:,:,1)=Outplane(:,:)
! 	elseif(Xnode.eq.1) then
! 	allocate( Temp_Send(nx,ny), STAT=s1)
!     	do j=1,ny
!     	do i=1,nx
!      	Temp_Send(i,j) = var(i,j,index1)
!     	enddo
!    	enddo
! 	call MPI_SEND(Temp_Recv,nx*ny*1,MPI_DOUBLE_PRECISION,0,1,commx1x2,status1,stat)
! 
! 	endif
! 	if (MYRANK.EQ.myidM2) then
	write(filename,'(a,i4.4,a,i8.8,a)') trim(varname)//"_k",index1,"_n",icycle,".pln"

  	open(210,file=filename,form='unformatted',status='unknown',iostat=s1) 
    	write(210) icycle,TIME,dtm,9.810d0,1.0d0,180.0,1.0
    	write(210) dir,index1, iu, iv, iw
    	write(210) zcg(index1), zwg(index1)
    	write(210) nx, ny
    	write(210) xc, xu
    	write(210) yc, yv
    	if (prec.EQ.0) then
     	allocate( SP_plane(1:nx,1:ny),STAT=s1 )
     	SP_plane=Temp_Recv !store plane in single precsion
     	write(210)  SP_plane
     	deallocate(SP_plane)
    	else 
     	write(210) Temp_Recv(:,:)
    	endif
	WRITE(*,*)"***************WROTE PLANE IN Z DIR********************"
   	endif 
	ENDIF
	
	close(210)
	err1=0
 	stat=err1
 	return
 
1000    continue
  
 	if (MYRANK.EQ.0) then
  	inquire(unit=210,opened=cfile)
  	if (cfile) then
   	write(210) "ERROR ERROR"
   	close(210)
  	endif
	endif
 	stat=-1
 	return
	end subroutine write_plane
! 
	subroutine Reduce_Plane_to_Master(var,dir,plane,Master,PX,sp1,sp2,nx,ny,nz,nzg,stat)
	INCLUDE 'common.h'
        INCLUDE 'mpif.h'

!Passed Variables
	integer,intent(in)             :: dir, plane, sp1, sp2,nx,ny,nz,nzg
 	real,intent(in)                :: var(nx,ny,nz)
	integer,intent(out)            :: Master
 	real,intent(out)               :: PX(1:sp1,1:sp2)
 	integer,intent(out)            :: stat
!Local Variables
 	integer                        :: Xmaster,Xnode
 	
 	Master = -1
 	if (dir.EQ.1) then     !0 
  !sub-communicator is x2x3
 	 Xmaster = -1
 	 Xnode   = -1
 	 if (plane.GT.sx-1.AND.plane.LT.ex+1) then !1 
!This node contains this plane

 	  if (rankx2x3.EQ.0) then 
 	   Xmaster = 1
 	   Master = MYRANK
  	 else 
 	   Xnode   = 1
 	  endif 

 	 elseif ( (plane.EQ.sx-1.AND.sx-1.EQ.1).OR.(plane.EQ.ex+1.AND.ex+1.EQ.nx) )  then !1 
!This node contains this plane and its a boundary plane

 	  if (rankx2x3.EQ.0) then 
 	   Xmaster = 1
  	  Master  = MYRANK
  	 else 
   	 Xnode   = 1
  	 endif 

  	elseif (plane.LT.1.or.plane.GT.nx) then  !1 
    !This is an invalid plane 
       	write(*,'(a,i1,a,i4,a)') "ILLEGAL PLANE REQUEST, DIRECTION:",dir,"PLANE",plane,"is not between 1 and nx"
  	 goto 1000

 	 else !1 
!This node does not contain this plane
  	 continue 
 	 endif !1

 	elseif (dir.EQ.2) then !0 
!sub-communicator is x2x3
  	Xmaster = -1
  	Xnode   = -1
  	if (plane.GT.sy-1.AND.plane.LT.ey+1) then !1 
!This node contains this plane

  	 if (rankx1x3.EQ.0) then 
   	 Xmaster = 1
   	 Master  = MYRANK
   	else 
   	 Xnode   = 1
   	endif 

  	elseif ( (plane.EQ.sy-1.AND.sy-1.EQ.1).OR.(plane.EQ.ey+1.AND.ey+1.EQ.ny) )  then !1 
!This node contains this plane and its a boundary plane

  	 if (rankx1x3.EQ.0) then 
   	 Xmaster = 1
  	  Master  = MYRANK
  	 else 
   	 Xnode   = 1
   	endif 

  	elseif (plane.LT.1.or.plane.GT.ny) then !1 
!This is an invalid plane 
  	 write(*,'(a,i1,a,i4,a)') "ILLEGAL PLANE REQUEST, DIRECTION:",dir,"PLANE",plane,"is not between 1 and nyp2"
   	goto 1000

  	else !1 
!This node does not contain this plane
   	continue 
  	endif !1

 	elseif (dir.EQ.3) then !0 !sub-communicator is x1x2
  	Xmaster = -1
  	Xnode   = -1
  	if (plane.GT.sz-1.AND.plane.LT.ez+1) then !1 
!This node contains this plane
   	if (rankx1x2.EQ.0) then 
   	 Xmaster = 1
    	Master  = MYRANK
   	else 
    	Xnode   = 1
   	endif 

  	elseif ( (plane.EQ.sz-1.AND.sz-1.EQ.1).OR.(plane.EQ.ez+1.AND.ez+1.EQ.nzg) )  then !1 
!This node contains this plane and its a boundary plane

   	if (rankx1x2.EQ.0) then 
    	Xmaster = 1
   	 Master  = MYRANK
  	 else 
   	 Xnode   = 1
   	endif 

  	elseif (plane.LT.1.or.plane.GT.nzg) then !1 
!This is an invalid plane 
   	write(*,'(a,i1,a,i4,a)') "ILLEGAL PLANE REQUEST, DIRECTION:",dir,"PLANE",plane,"is not between 1 and nzp2"
  	 goto 1000
  	else !1 
!This node does not contain this plane
   	continue 
 	 endif !1

        
 	 else
 
  	endif
	
  	if (Xmaster.EQ.1)call gather2dM(var,PX,sp1,sp2,dir,plane,nx,ny,nz,nzg,0,stat)
  	if (Xnode.EQ.1) call gather2dS(var,dir,plane,nx,ny,nz,nzg,0,stat)


!  stat=0 stat should contain the stat from the gather2d calls


  	return

1000    continue
  	stat=1
  	return
	end subroutine Reduce_Plane_to_Master 
! 
! 
	subroutine gather2dM(varL,OutPlane,sb1,sb2,dir,plane,nx,ny,nz,nzg,myidM,stat)
	include'common.h'
        include'mpif.h'
 	 	
!Passed Variables 
 	integer,intent(in)        :: myidM,dir,plane,sb1,sb2,nx,ny,nz,nzg
 	real,intent(in)           :: varL(nx,ny,nz)
 	real,intent(out)          :: OutPlane(1:sb1,1:sb2)
 	integer,intent(out)       :: stat

!Local Variables
 	integer                   :: Tsize, Rcoords(3), s1, status1(MPI_STATUS_SIZE)
 	integer                   :: i,j,k,n,istart, jstart, kstart, i2, j2, k2
 	real,allocatable,dimension(:,:) :: Temp_Recv, Temp_Send
 	integer                   :: is,js,ks,ie,je,ke
	integer,parameter      :: realtype  =MPI_DOUBLE_PRECISION
        integer,parameter      :: inttype   =MPI_INTEGER

 !Gather Data

	if (dir.eq.1) then !0
  	allocate( Temp_Recv(ny,nz), STAT=s1 )
  	Tsize=size(Temp_Recv)
  	allocate( Temp_Send(ny,nz), STAT=s1)
  	do k=1,nz
  	 do j=1,ny
   	 Temp_Send(j,k) = varL(plane,j,k)
  	 enddo
  	enddo
!  call MPI_SEND(coords,3,inttype,myidM,2,commx2x3,stat)
!  call MPI_SEND(Temp_Send,Tsize,realtype,myidM,1,commx2x3,stat)
 
  	do n=0,sizex2x3-1
   	if (n.EQ.myidM) then
    	Rcoords=coords
    	Temp_Recv=Temp_Send
   	else
    	call MPI_RECV(Rcoords,3,inttype,n,2,commx2x3,status1,stat)
    	call MPI_RECV(Temp_Recv,Tsize,realtype,n,1,commx2x3,status1,stat)
   	endif
  !Determine Block of Data to recieve
   	jstart = Rcoords(2)*(ny-2)/nyprocs
   	kstart = Rcoords(3)*(nzg-2)/nzprocs

  !Determine if there is boundary data 
   	js=0
   	je=0
   	ks=0
   	ke=0
   	if ( Rcoords(2).EQ.0      )   js=1
   	if ( Rcoords(3).EQ.0      )   ks=1
   	if ( Rcoords(2).EQ.sizex2-1 ) je=1
   	if ( Rcoords(3).EQ.sizex3-1 ) ke=1

   !UnPack Data
   	do k=2-ks,nz-1+ke  
   	 do j=2-js,ny-1+je 
      	j2=jstart+j
      	k2=kstart+k
      	Outplane(j2,k2)=Temp_Recv(j,k)
     	enddo
    	enddo
  	enddo

 	elseif (dir.eq.2) then !0
	allocate( Temp_Recv(nx,nz), STAT=s1 )
  	Tsize=size(Temp_Recv)
  	allocate( Temp_Send(nx,nz), STAT=s1 )
  	do k=1,nz
   	do i=1,nx
    	Temp_Send(i,k) = varL(i,plane,k)
	
  	enddo
  	enddo
	
  	do n=0,sizex1x3-1
   	if (n.EQ.myidM) then
    	Rcoords=coords
    	Temp_Recv=Temp_Send
	istart=0
	kstart=0
	do k=2,nz  
    	do i=2,nx 
     	 i2=istart+i
     	 k2=kstart+k
      	Outplane(i2,k2)=Temp_Recv(i,k)
     	enddo
    	enddo
   	else
    	call MPI_RECV(Rcoords,3,inttype,n,2,commx1x3,status1,stat)
    	call MPI_RECV(Temp_Recv,Tsize,realtype,n,1,commx1x3,status1,stat)
	istart = (nx-2)/nxprocs	 	!change if the number of processor changes in x dir
	kstart = n*(nzg-2)/nzprocs 
	do k=2,nz  
    	do i=2,nx 
     	 i2=istart*0+i
     	 k2=kstart+k
      	Outplane(i2,k2)=Temp_Recv(i,k)
     	enddo
    	enddo
	write(*,*) istart, kstart
   	endif

!   !Determine Block of Data to recieve
!    	istart = Rcoords(1)*(nx-2)/nxprocs
!    	kstart = Rcoords(3)*(nzg-2)/nzprocs
! ! 	WRITE(*,*)"********SHIT***********",ISTART,KSTART
!   !Determine if there is boundary data 
!    	is=0
!    	ie=0
!    	ks=0
!    	ke=0
!    	if ( Rcoords(1).EQ.0      )   is=1
!    	if ( Rcoords(3).EQ.0      )   ks=1
!    	if ( Rcoords(1).EQ.sizex1-1 ) ie=1
!    	if ( Rcoords(3).EQ.sizex3-1 ) ke=1
! 
!    !UnPack Data
!    	do k=2-ks,nz-1+ke  
!     	do i=2-is,nx-1+ie 
!      	 i2=istart+i
!      	 k2=kstart+k
!       	Outplane(i2,k2)=Temp_Recv(i,k)
!      	enddo
!     	enddo
   	enddo

	elseif (dir.eq.3) then !0
	allocate( Temp_Recv(sx-1:ex+1,sy-1:ey+1), STAT=s1 )
	Tsize=size(Temp_Recv)
  	allocate( Temp_Send(sx-1:ex+1,sy-1:ey+1), STAT=s1 )
	
  	do j=1,ny
   	do i=1,nx
    	Temp_Send(i,j) = varL(i,j,plane)
   	enddo
  	enddo
	do n=0,sizex1x2-1
   	if (n.EQ.myidM) then
    	Rcoords=coords
    	Temp_Recv=Temp_Send
   	else
    	call MPI_RECV(Rcoords,3,inttype,n,2,commx1x2,status1,stat)
    	call MPI_RECV(Temp_Recv,Tsize,realtype,n,1,commx1x2,status1,stat)
   	endif
	
  !Determine Block of Data to recieve
   	istart = Rcoords(1)*(nx-2)/nxprocs
   	jstart = Rcoords(2)*(ny-2)/nyprocs

  !Determine if there is boundary data 
   	is=0
   	ie=0
   	js=0
   	je=0
   	if ( Rcoords(1).EQ.0      )   is=1
   	if ( Rcoords(2).EQ.0      )   js=1
   	if ( Rcoords(1).EQ.sizex1-1 ) ie=1
   	if ( Rcoords(2).EQ.sizex2-1 ) je=1

   !UnPack Data
   	do j=2-js,ny-1+je  
    	do i=2-is,nx-1+ie 
      	i2=istart+i
      	j2=jstart+j
      	Outplane(i2,j2)=Temp_Recv(i,j)
     	enddo
    	enddo
   	enddo

 	else !0
  	write(*,'(a,i1,a)') "INVALID DIRECTION:", dir, "not (1,2,3)"
  	goto 2000

 	endif !0

 	deallocate(Temp_Send,STAT=s1)
 	deallocate(Temp_Recv,STAT=s1)


 	return

 2000 	continue
 	stat=1
 	return
	end subroutine gather2dM

	subroutine gather2dS(varL,dir,plane,nx,ny,nz,nzg,myidM,stat)
	include'common.h'
        include'mpif.h'
 	

!Passed Variables 
 	integer,intent(in)        :: myidM, dir, plane,nx,ny,nz,nzg
 	real,intent(in)       :: varL(nx,ny,nz)
 	integer,intent(out)       :: stat

!Local Variables
        integer,parameter      :: realtype  =MPI_DOUBLE_PRECISION
        integer,parameter      :: inttype   =MPI_INTEGER
 	integer                            :: Tsize
 	integer                            :: i,j,k 
 	integer                            :: s1
 	real,allocatable,dimension(:,:) :: Temp_Send
	if (dir.eq.1) then !0
   	allocate( Temp_Send(ny,nz), STAT=s1)
   	Tsize=size(Temp_Send)
   	do k=1,nz
    	do j=1,ny
     	Temp_Send(j,k) = varL(plane,j,k)
    	enddo
   	enddo
   	call MPI_SEND(coords,3,inttype,myidM,2,commx2x3,stat)
   	call MPI_SEND(Temp_Send,Tsize,realtype,myidM,1,commx2x3,stat)

 	elseif (dir.eq.2) then !0
	allocate( Temp_Send(nx,nz), STAT=s1)
   	Tsize=size(Temp_Send)
	do k=1,nz
    	do i=1,nx
     	Temp_Send(i,k) = varL(i,plane,k)
	enddo
   	enddo
	call MPI_SEND(coords,3,inttype,myidM,2,commx1x3,stat)
	call MPI_SEND(Temp_Send,Tsize,realtype,myidM,1,commx1x3,stat)
	elseif (dir.eq.3) then !0
   	allocate( Temp_Send(nx,ny), STAT=s1)
   	Tsize=size(Temp_Send)
   	do j=1,ny
    	do i=1,nx
     	Temp_Send(i,j) = varL(i,j,plane)
   	 enddo
   	enddo
   	call MPI_SEND(coords,3,inttype,myidM,2,commx1x2,stat)
   	call MPI_SEND(Temp_Send,Tsize,realtype,myidM,1,commx1x2,stat)
 	else !0
  	write(*,'(a,i1,a)') "INVALID DIRECTION:", dir, "not (1,2,3)"
  	goto 2000
 	endif !0

 	deallocate( Temp_Send, STAT=s1 )

 2000 	continue
 	stat=1 !Unsucessful Return
	return
	end subroutine gather2dS
! 


        SUBROUTINE DERIVATIVE(VAR,DERI,NX,NY,NZ,NZG,XC,ZCG,DIR,IVAR)
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'
        INTEGER :: NX,NY,NZ,NZG
        REAL    :: XU(NX),ZWG(NZG)
        REAL    :: XC(NX),YC(NY),ZC(NZ),ZCG(NZG)
        INTEGER :: I,J,K,DIR,IVAR
        REAL    :: VAR(NX,NY,NZ),DERI(NX,NY,NZ)
        REAL    :: dtheta,dr,dz

        IF(DIR.EQ.1) THEN
        do k=kz1,kz2
        do j=jy1,jy2
        do i=ix1,ix2
          dr=xc(i)-xc(i-1)
          IF(IVAR.EQ.1) THEN
          DERI(I,J,K) = (VAR(I,J,K)-VAR(I-1,J,K))/dr
          ELSEIF(IVAR.EQ.2) THEN
          DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j-1,k)+VAR(i+1,j,k)+VAR(i+1,j-1,k))-0.25*(VAR(i,j,k)+
     &    VAR(i,j-1,k)+VAR(i-1,j,k)+VAR(i-1,j-1,k)))/dr
          ELSEIF(IVAR.EQ.3) THEN
          DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j,k-1)+VAR(i+1,j,k)+VAR(i+1,j,k-1))-0.25*(VAR(i,j,k)+
     &    VAR(i,j,k-1)+VAR(i-1,j,k)+VAR(i-1,j,k-1)))/dr
          ENDIF
        enddo
        enddo
        enddo
        ELSEIF(DIR.EQ.2) THEN
        do k=kz1,kz2
        do j=jy1,jy2
        do i=ix1,ix2
          dtheta=dely
          IF(IVAR.EQ.1) THEN
          DERI(I,J,K) = (0.25*(VAR(i,j,k)+VAR(i,j+1,k)+VAR(i-1,j,k)+VAR(i-1,j+1,k))-0.25*(VAR(i,j,k)+
     &    VAR(i,j-1,k)+VAR(i-1,j,k)+VAR(i-1,j-1,k)))/dtheta
          ELSEIF(IVAR.EQ.2) THEN
          DERI(I,J,K) = (VAR(I,J,K)-VAR(I,J-1,K))/dtheta
          ELSEIF(IVAR.EQ.3) THEN
          DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j+1,k)+VAR(i,j,k-1)+VAR(i,j+1,k-1))-0.25*(VAR(i,j,k)+
     &    VAR(i,j-1,k)+VAR(i,j,k-1)+VAR(i,j-1,k-1)))/dtheta
          ENDIF
        enddo
        enddo
        enddo
        ELSEIF(DIR.EQ.3) THEN
        do k=kz1,kz2
        do j=jy1,jy2
        do i=ix1,ix2
          dz=zcg(k)-zcg(k-1)
          IF(IVAR.EQ.1) THEN
          DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j,k+1)+VAR(i-1,j,k)+VAR(i-1,j,k+1))-0.25*(VAR(i,j,k)+
     &    VAR(i,j,k-1)+VAR(i-1,j,k)+VAR(i-1,j,k-1)))/dz
          ELSEIF(IVAR.EQ.2) THEN
          DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j-1,k)+VAR(i,j,k+1)+VAR(i,j-1,k+1))-0.25*(VAR(i,j,k)+
     &    VAR(i,j-1,k)+VAR(i,j,k-1)+VAR(i,j-1,k-1)))/dz
          ELSEIF(IVAR.EQ.3) THEN
          DERI(I,J,K) = (VAR(I,J,K)-VAR(I,J,K-1))/dz
          ENDIF
        enddo
        enddo
        enddo
        ENDIF

        RETURN
        END


        SUBROUTINE DERIVATIVE_1CPU(VAR,DERI,NX,NY,NZ,NZG,XC,ZCG,DIR,IVAR)
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'
        INTEGER :: NX,NY,NZ,NZG
        REAL    :: XU(NX),ZWG(NZG)
        REAL    :: XC(NX),YC(NY),ZC(NZ),ZCG(NZG)
        INTEGER :: I,J,K,DIR,IVAR
        REAL    :: VAR(NX,NY,NZ),DERI(NX,NY,NZ)
        REAL    :: dtheta,dr,dz

        IF(DIR.EQ.1) THEN
        do k=1,NZG
        do j=1,NY
        do i=1,NX
          dr=xc(i)-xc(i-1)
          IF(IVAR.EQ.1) THEN
          DERI(I,J,K) = (VAR(I,J,K)-VAR(I-1,J,K))/dr
          ELSEIF(IVAR.EQ.2) THEN
          DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j-1,k)+VAR(i+1,j,k)+VAR(i+1,j-1,k))-0.25*(VAR(i,j,k)+
     &    VAR(i,j-1,k)+VAR(i-1,j,k)+VAR(i-1,j-1,k)))/dr
          ELSEIF(IVAR.EQ.3) THEN
          DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j,k-1)+VAR(i+1,j,k)+VAR(i+1,j,k-1))-0.25*(VAR(i,j,k)+
     &    VAR(i,j,k-1)+VAR(i-1,j,k)+VAR(i-1,j,k-1)))/dr
          ENDIF
        enddo
        enddo
        enddo
        ELSEIF(DIR.EQ.2) THEN
        do k=1,NZG
        do j=1,NY
        do i=1,NX
          dtheta=dely
          IF(IVAR.EQ.1) THEN
          DERI(I,J,K) = (0.25*(VAR(i,j,k)+VAR(i,j+1,k)+VAR(i-1,j,k)+VAR(i-1,j+1,k))-0.25*(VAR(i,j,k)+
     &    VAR(i,j-1,k)+VAR(i-1,j,k)+VAR(i-1,j-1,k)))/dtheta
          ELSEIF(IVAR.EQ.2) THEN
          DERI(I,J,K) = (VAR(I,J,K)-VAR(I,J-1,K))/dtheta
          ELSEIF(IVAR.EQ.3) THEN
          DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j+1,k)+VAR(i,j,k-1)+VAR(i,j+1,k-1))-0.25*(VAR(i,j,k)+
     &    VAR(i,j-1,k)+VAR(i,j,k-1)+VAR(i,j-1,k-1)))/dtheta
          ENDIF
        enddo
        enddo
        enddo
        ELSEIF(DIR.EQ.3) THEN
        do k=1,NZG
        do j=1,NY
        do i=1,NX
          dz=zcg(k)-zcg(k-1)
          IF(IVAR.EQ.1) THEN
          DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j,k+1)+VAR(i-1,j,k)+VAR(i-1,j,k+1))-0.25*(VAR(i,j,k)+
     &    VAR(i,j,k-1)+VAR(i-1,j,k)+VAR(i-1,j,k-1)))/dz
          ELSEIF(IVAR.EQ.2) THEN
          DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j-1,k)+VAR(i,j,k+1)+VAR(i,j-1,k+1))-0.25*(VAR(i,j,k)+
     &    VAR(i,j-1,k)+VAR(i,j,k-1)+VAR(i,j-1,k-1)))/dz
          ELSEIF(IVAR.EQ.3) THEN
          DERI(I,J,K) = (VAR(I,J,K)-VAR(I,J,K-1))/dz
          ENDIF
        enddo
        enddo
        enddo
        ENDIF

        RETURN
        END

! WRONG!!!1 DON'T USE
c$$$        SUBROUTINE DERIVATIVE_CENTER_1CPU(VAR,DERI,NX,NY,NZ,NZG,XU,YV,ZWG,DIR)
c$$$	! VAR is at cell-center
c$$$	! 1 CPU
c$$$        INCLUDE 'common.h'
c$$$        INCLUDE 'mpif.h'
c$$$        INTEGER :: NX,NY,NZ,NZG
c$$$        REAL    :: XU(NX),YV(NY),ZWG(NZG)
c$$$        INTEGER :: I,J,K,DIR,IVAR
c$$$        REAL    :: VAR(NX,NY,NZG),DERI(NX,NY,NZG)
c$$$        REAL    :: dtheta,dr,dz
c$$$
c$$$        IF(DIR.EQ.1) THEN
c$$$c$$$        do k=kz1,kz2
c$$$c$$$        do j=jy1,jy2
c$$$c$$$        do i=ix1,ix2
c$$$ 	do i = 2,NX-1
c$$$	do j = 2,NY-1
c$$$	do k = 2,NZG-1
c$$$
c$$$c          dr=xc(i)-xc(i-1)
c$$$          dr=xu(i)-xu(i-1)
c$$$	  
c$$$	  DERI(I,J,K) = (0.5*(VAR(I,J,K)+VAR(I+1,J,K)) - 0.5*(VAR(I,J,K)+VAR(I-1,J,K)))/dr
c$$$	  
c$$$        enddo
c$$$        enddo
c$$$        enddo
c$$$        ELSEIF(DIR.EQ.2) THEN
c$$$c$$$        do k=kz1,kz2
c$$$c$$$        do j=jy1,jy2
c$$$c$$$        do i=ix1,ix2
c$$$ 	do i = 2,NX-1
c$$$	do j = 2,NY-1
c$$$	do k = 2,NZG-1
c$$$
c$$$c          dtheta=dely
c$$$          dtheta=yv(j)-yv(j-1)
c$$$
c$$$	  DERI(I,J,K) = (0.5*(VAR(I,J,K)+VAR(I,J+1,K)) - 0.5*(VAR(I,J,K)+VAR(I,J-1,K)))/dtheta
c$$$
c$$$        enddo
c$$$        enddo
c$$$        enddo
c$$$        ELSEIF(DIR.EQ.3) THEN
c$$$c$$$        do k=kz1,kz2
c$$$c$$$        do j=jy1,jy2
c$$$c$$$        do i=ix1,ix2
c$$$ 	do i = 2,NX-1
c$$$	do j = 2,NY-1
c$$$	do k = 2,NZG-1
c$$$
c$$$          dz=zwg(k)-zwg(k-1)
c$$$
c$$$	  DERI(I,J,K) = (0.5*(VAR(I,J,K)+VAR(I,J,K+1)) - 0.5*(VAR(I,J,K)+VAR(I,J,K-1)))/dz
c$$$
c$$$        enddo
c$$$        enddo
c$$$        enddo
c$$$        ENDIF
c$$$
c$$$        RETURN
c$$$        END


	SUBROUTINE PERIODIC_1D_AVERAGE_RMS(NX,NY,NZ,NZG,NBD,ICYCLE,time,DTM1,XC,XU,YC,YV,ZCG,ZWG,xc_car,
     &                  yc_car,xu_car,yu_car,xv_car,yv_car,P,VO,UO,WO,DENS)

        use vorticity
        use density_bg
	INCLUDE 'common.h'
        INCLUDE 'mpif.h'
	INTEGER :: NX,NY,NZ,NZG,NBD,ICYCLE
	REAL	:: DTM1,TIME
	REAL    :: XU(NX),YV(NY),ZW(NZ),ZWG(NZG)
	REAL    :: XC(NX),YC(NY),ZC(NZ),ZCG(NZG)
	REAL    :: UO(NX,NY,NZ),VO(NX,NY,NZ),WO(NX,NY,NZ),P(NX,NY,NZ),DENS(NX,NY,NZ)
	REAL	:: XU_CAR(NX,NY),YU_CAR(NX,NY),XV_CAR(NX,NY),YV_CAR(NX,NY),XC_CAR(NX,NY),YC_CAR(NX,NY)
	INTEGER :: I,J,K,STAT
        integer iic,ic,ip,im,jjc,jc,jp,jm,jsy,jpsy,kkc,kc,kp,km
        real dtheta,dr,dz,dq1x3,dq3x1,dq2x3,dq3x2,dq2x1,dq1x2
	REAL,ALLOCATABLE,DIMENSION(:,:,:)::U1_TMP2,U2_TMP2,U3_TMP2,rtmp1

	REAL    :: UMEAN(NX,NZ), VMEAN(NX,NZ), WMEAN(NX,NZ), DENSMEAN(NX,NZ)
	REAL    ::  URMS(NX,NZ),  VRMS(NX,NZ),  WRMS(NX,NZ),  DENSRMS(NX,NZ)
	REAL    ::  TMP1(NX,NZ),  TMP2(NX,NZ),  TMP3(NX,NZ),     TMP4(NX,NZ)

c	if(myrank ==0 ) write(*,*) "NY = ", NY

	DO I=IX1,IX2
	DO K=KZ1,KZ2
	   UMEAN(I,K)=SUM(UO(I,JY1:JY2,K))/REAL(NY-2)
	   VMEAN(I,K)=SUM(VO(I,JY1:JY2,K))/REAL(NY-2)
	   WMEAN(I,K)=SUM(WO(I,JY1:JY2,K))/REAL(NY-2)
	DENSMEAN(I,K)=SUM(DENS(I,JY1:JY2,K))/REAL(NY-2)
	ENDDO
	ENDDO

	DO I=IX1,IX2
	DO K=KZ1,KZ2

	   TMP1(I,K) = 0.0
	   TMP2(I,K) = 0.0
	   TMP3(I,K) = 0.0
	   TMP4(I,K) = 0.0
	DO J=JY1,JY2
	   TMP1(I,K) = TMP1(I,K) + (UO(I,J,K) - UMEAN(I,K))**2.0
	   TMP2(I,K) = TMP2(I,K) + (VO(I,J,K) - VMEAN(I,K))**2.0
	   TMP3(I,K) = TMP3(I,K) + (WO(I,J,K) - WMEAN(I,K))**2.0

	   TMP4(I,K) = TMP4(I,K) + (DENS(I,J,K) - DENSMEAN(I,K))**2.0
	ENDDO
	
	URMS(I,K) = SQRT(TMP1(I,K)/REAL(NY-2))
	VRMS(I,K) = SQRT(TMP2(I,K)/REAL(NY-2))
	WRMS(I,K) = SQRT(TMP3(I,K)/REAL(NY-2))
      DENSRMS(I,K) = SQRT(TMP4(I,K)/REAL(NY-2))


	ENDDO
	ENDDO

 	call WRITE_2D_AVERAGE_PLANE(UMEAN,2,1,0,0,'UMEAN',.FALSE.,nx,ny,nz,nzg,xc,xu,yc,yv,zcg,zwg,XC_CAR,
     &	                    YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,icycle,TIME,dtm1,stat)
 	call WRITE_2D_AVERAGE_PLANE(VMEAN,2,1,0,0,'VMEAN',.FALSE.,nx,ny,nz,nzg,xc,xu,yc,yv,zcg,zwg,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,icycle,TIME,dtm1,stat)
 	call WRITE_2D_AVERAGE_PLANE(WMEAN,2,1,0,0,'WMEAN',.FALSE.,nx,ny,nz,nzg,xc,xu,yc,yv,zcg,zwg,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,icycle,TIME,dtm1,stat)
 	call WRITE_2D_AVERAGE_PLANE(DENSMEAN,2,1,0,0,'DENSMEAN',.FALSE.,nx,ny,nz,nzg,xc,xu,yc,yv,zcg,zwg,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,icycle,TIME,dtm1,stat)


 	call WRITE_2D_AVERAGE_PLANE(URMS,2,1,0,0,'URMS',.FALSE.,nx,ny,nz,nzg,xc,xu,yc,yv,zcg,zwg,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,icycle,TIME,dtm1,stat)
 	call WRITE_2D_AVERAGE_PLANE(VRMS,2,1,0,0,'VRMS',.FALSE.,nx,ny,nz,nzg,xc,xu,yc,yv,zcg,zwg,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,icycle,TIME,dtm1,stat)
 	call WRITE_2D_AVERAGE_PLANE(WRMS,2,1,0,0,'WRMS',.FALSE.,nx,ny,nz,nzg,xc,xu,yc,yv,zcg,zwg,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,icycle,TIME,dtm1,stat)
 	call WRITE_2D_AVERAGE_PLANE(DENSRMS,2,1,0,0,'DENSRMS',.FALSE.,nx,ny,nz,nzg,xc,xu,yc,yv,zcg,zwg,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,icycle,TIME,dtm1,stat)

	END

	SUBROUTINE WRITE_2D_AVERAGE_PLANE(var,dir,index1,myidM,prec,varname,verbose,nx,ny,nz,nzg,xc,xu,yc,yv,zcg,zwg,XC_CAR,
     &				YC_CAR,XU_CAR,YU_CAR,XV_CAR,YV_CAR,icycle,TIME,dtm,stat)

	INCLUDE 'common.h'
        INCLUDE 'mpif.h'


!Passed Variables
	integer,intent(in)    :: dir,index1,myidM,prec,nx,ny,nz,nzg,icycle
 	real,intent(in)       :: var(nx,nz),TIME,dtm
 	real,intent(in)       :: xc(nx),xu(nx),yc(ny),yv(ny),zcg(nzg),zwg(nzg)
 	integer,intent(out)   :: stat
 	character(len=*)      :: varname
 	logical,intent(in)    :: verbose

!Local Variables
 	integer :: Msub, MsubG
        integer,parameter      :: realtype  =MPI_DOUBLE_PRECISION
        integer,parameter      :: inttype   =MPI_INTEGER
	integer :: status1(MPI_STATUS_SIZE),index2
 	integer :: err1, s1,Rcoords(3),myidM2
 	character(len=300) :: filename
 	real(4),allocatable,dimension(:,:) :: SP_plane,TEMP
 	logical             :: cfile
 	integer             :: iu, iv, iw,N,Tsize,i,j,k,ks,kstart,sp1,sp2,is,js,istart,jstart
   	real,allocatable,dimension(:,:,:) :: PlnX1
   	real,allocatable,dimension(:,:,:) :: PlnX2
   	real,allocatable,dimension(:,:,:) :: PlnX3
	real,allocatable,dimension(:,:) :: Outplane
	real,allocatable,dimension(:,:) :: Temp_Recv, Temp_Send
	integer                            :: Xmaster,Xnode
	REAL	:: XU_CAR(NX,NY),YU_CAR(NX,NY),XV_CAR(NX,NY),YV_CAR(NX,NY),XC_CAR(NX,NY),YC_CAR(NX,NY)

	allocate(PlnX1(1:ny,1:nzg,2),PlnX2(1:nx,1:nzg,1:2),PlnX3(1:nx,1:ny,2))
	iu=0
 	iv=0
 	iw=0
	select case(varname)
 	case('uo')
 	iu=1 
 	case('vo')
 	iv=1
 	case('wo')
 	iw=1
	end select
	
 	if (MYRANK.eq.myidM) then
	if (dir.EQ.1) write(filename,'(a,i4.4,a,i8.8,a)') trim(varname)//"_i",index1,"_n",icycle,".pln"
  	if (dir.EQ.2) write(filename,'(a,i4.4,a,i8.8,a)') trim(varname)//"_j",index1,"_n",icycle,".pln"
!   	if (dir.EQ.3) write(filename,'(a,i4.4,a,i5.5,a)') trim(varname)//"_k",index1,"_n",icycle,".pln"

  	open(210,file=filename,form='unformatted',status='unknown',iostat=s1)
  	endif
	
	allocate(Temp_Recv(nx,nz))
	Temp_Recv(:,:)=var(:,:)
 	sp1=size(PlnX2,1)
 	sp2=size(PlnX2,2)
 	allocate(Outplane(sp1,sp2))
 	if(myrank.eq.0) then

 	DO n=0,sizex1x3-1
 	  if(n.eq.0) then
 	    do i=2,nx-1
 	     do k=2,nz-1
 	      Outplane(i,k)=Temp_Recv(i,k)
 	     enddo
 	    enddo
 	  else
 	      call MPI_RECV(Temp_Recv,nx*1*nz,MPI_DOUBLE_PRECISION,n,1,commx1x3,status1,stat)
 	  endif
 	 if(n.ne.0) then
c	    write(*,*), coords(1),coords(2),coords(3)
 	      kstart=n*(nzg-2)/nzprocs
 	      istart=(nx-2)/nxprocs
c 	      write(*,*), kstart
 	 endif
 	      do i=2,nx-1
 	       do k=2,nz-1
 	        ks=kstart+k
 		is=istart+i
  	        Outplane(i,ks)=Temp_Recv(i,k)
 	       enddo
 	      enddo
 	ENDDO
 	PlnX2(:,:,1)=Outplane(:,:)
 	else
 	    call MPI_SEND(Temp_Recv,nx*1*nz,MPI_DOUBLE_PRECISION,0,1,commx1x3,status1,stat)

 	endif

	
  	if (MYRANK.EQ.myidM) then 
      	write(210) icycle,TIME,dtm,9.810d0,1.0d0,180.0,1.0
      	write(210) dir,index1, iu, iv, iw
      	write(210) yc(index1), yv(index1)
      	write(210) nx, nzg
      	write(210) xc, xu
      	write(210) zcg, zwg
      	if (prec.EQ.0) then
	   allocate( SP_plane(1:nx,1:nzg),STAT=s1 )
	   SP_plane= OUTPLANE(:,:)  !PlnX2(:,:,1) !store plane in single precsion

	   
	   write(210)  SP_plane
	   deallocate(SP_plane)
      	else 
	   write(210) PlnX2(:,:,1)
      	endif
 	endif
	
	END

	SUBROUTINE PERIODIC_2D_AVERAGE_RMS_1CPU(NX,NY,NZ,NZG,NBD,ICYCLE,time,DTM1,XC,XU,YC,YV,ZCG,ZWG,xc_car,
     &                  yc_car,xu_car,yu_car,xv_car,yv_car,P,VO,UO,WO,DENS)

c       USE 1 CPU
c       perform average in periodic direction (y & z)
c       calculate reynolds stresses ww,uu,vv,wu
c        write 
c        x+(vertical), w+(streamwise vel), Rww, Ruu, Rvv, Rwu


        use vorticity
        use density_bg
	INCLUDE 'common.h'
        INCLUDE 'mpif.h'
	INTEGER :: NX,NY,NZ,NZG,NBD,ICYCLE
	REAL	:: DTM1,TIME
	REAL    :: XU(NX),YV(NY),ZW(NZ),ZWG(NZG)
	REAL    :: XC(NX),YC(NY),ZC(NZ),ZCG(NZG)
	REAL    :: UO(NX,NY,NZ),VO(NX,NY,NZ),WO(NX,NY,NZ),P(NX,NY,NZ),DENS(NX,NY,NZ)
	REAL	:: XU_CAR(NX,NY),YU_CAR(NX,NY),XV_CAR(NX,NY),YV_CAR(NX,NY),XC_CAR(NX,NY),YC_CAR(NX,NY)
	INTEGER :: I,J,K,STAT
        integer iic,ic,ip,im,jjc,jc,jp,jm,jsy,jpsy,kkc,kc,kp,km
        real dtheta,dr,dz,dq1x3,dq3x1,dq2x3,dq3x2,dq2x1,dq1x2
	REAL,ALLOCATABLE,DIMENSION(:,:,:)::U1_TMP2,U2_TMP2,U3_TMP2,rtmp1

	REAL    :: UMEAN(NX),VMEAN(NX),WMEAN(NX)
	REAL    :: RWW(NX),RUU(NX),RVV(NX),RWV(NX),RWU(NX)
	REAL    :: TMP1(NX),TMP2(NX),TMP3(NX),TMP4(NX),TMP5(NX)

	DO I=IX1,IX2
	   UMEAN(I)=SUM(UO(I,JY1:JY2,KZ1:KZ2))/(REAL(NY-2)*REAL(NZ-2))
	   VMEAN(I)=SUM(VO(I,JY1:JY2,KZ1:KZ2))/(REAL(NY-2)*REAL(NZ-2))
	   WMEAN(I)=SUM(WO(I,JY1:JY2,KZ1:KZ2))/(REAL(NY-2)*REAL(NZ-2))
	ENDDO

	
	DO I=IX1,IX2
	   TMP1(I) = 0.0
	   TMP2(I) = 0.0
	   TMP3(I) = 0.0
	   TMP4(I) = 0.0
	   TMP5(I) = 0.0
	   DO K=KZ1,KZ2
	   DO J=JY1,JY2	   
	      TMP1(I) = TMP1(I) + (WO(I,J,K) - WMEAN(I))*(WO(I,J,K) - WMEAN(I))
	      TMP2(I) = TMP2(I) + (UO(I,J,K) - UMEAN(I))*(UO(I,J,K) - UMEAN(I))
	      TMP3(I) = TMP3(I) + (VO(I,J,K) - VMEAN(I))*(VO(I,J,K) - VMEAN(I))
	      TMP4(I) = TMP4(I) + (WO(I,J,K) - WMEAN(I))*(VO(I,J,K) - VMEAN(I))	      
	      TMP5(I) = TMP5(I) + (WO(I,J,K) - WMEAN(I))*(UO(I,J,K) - UMEAN(I))	      
	   ENDDO
	   ENDDO
	   RWW(I) = TMP1(I)/(REAL(NY-2)*REAL(NZ-2))
	   RUU(I) = TMP2(I)/(REAL(NY-2)*REAL(NZ-2))
	   RVV(I) = TMP3(I)/(REAL(NY-2)*REAL(NZ-2))
	   RWV(I) = TMP4(I)/(REAL(NY-2)*REAL(NZ-2))	   
	   RWU(I) = TMP5(I)/(REAL(NY-2)*REAL(NZ-2))	   
	ENDDO
	
	open(unit=10,file='w_rs.dat',form='formatted')
	DO I=IX1,IX2
	   write(10,'(7F15.8)') xc(I)*(1.0/ru1),WMEAN(I),
     &                          RWW(I),RUU(I),RVV(I),RWV(I),RWU(I)
	ENDDO
	close(10)

	END


C$$$	SUBROUTINE PERIODIC_2D_AVERAGE_RMS_NCPU(NX,NY,NZ,NZG,NBD,ICYCLE,time,DTM1,XC,XU,YC,YV,ZCG,ZWG,xc_car,
C$$$     &                  yc_car,xu_car,yu_car,xv_car,yv_car,P,VO,UO,WO,DENS,TV)

C$$$C       CALCULATE SPATIAL AND TIME AVERAGE AND RMS OF VELOCITIES USING THE FOLLOWING STEPS
C$$$C       1) FOR EACH 100 STEPS, CALCULATE SPATIAL AVERAGE AND RMS
C$$$C       2) WRITE A FILE CONTAINING TIME, xc(I)*(1.0/ru1),WMEAN(I),
C$$$C                                  RWW(I),RUU(I),RVV(I),RWV(I),RWU(I)
C$$$C       3) REPEAT
C$$$C       4) NOW WE HAVE FILES (SNAPSHOT) CONTAINING ACUMULATE TIME AND STATS
C$$$C       5) USE MATLAB TO READ EACH FILE AND PERFORM TIME AVERAGE

C$$$        use vorticity
C$$$        use density_bg
C$$$	INCLUDE 'common.h'
C$$$        INCLUDE 'mpif.h'
C$$$	INTEGER :: n,NX,NY,NZ,NZG,NBD,ICYCLE
C$$$	REAL	:: DTM1,TIME
C$$$	REAL    :: XU(NX),YV(NY),ZW(NZ),ZWG(NZG)
C$$$	REAL    :: XC(NX),YC(NY),ZC(NZ),ZCG(NZG)
C$$$	REAL    :: UO(NX,NY,NZ),VO(NX,NY,NZ),WO(NX,NY,NZ),P(NX,NY,NZ)
C$$$	REAL    :: DENS(NX,NY,NZ),TV(NX,NY,NZ)
C$$$	REAL	:: XU_CAR(NX,NY),YU_CAR(NX,NY),XV_CAR(NX,NY),YV_CAR(NX,NY),XC_CAR(NX,NY),YC_CAR(NX,NY)
C$$$	INTEGER :: I,J,K,STAT,kstart
C$$$	integer :: status1(MPI_STATUS_SIZE),index2
C$$$        integer iic,ic,ip,im,jjc,jc,jp,jm,jsy,jpsy,kkc,kc,kp,km
C$$$        real dtheta,dr,dz,dq1x3,dq3x1,dq2x3,dq3x2,dq2x1,dq1x2

C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:)::U_TMP,V_TMP,W_TMP,DENS_TMP,DENSF_TMP,TV_TMP
C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:)::U_GLOBAL,V_GLOBAL,W_GLOBAL
C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:)::DENS_GLOBAL,DENSF_GLOBAL,TV_GLOBAL

C$$$	REAL    :: CF

C$$$	REAL    :: UMEAN(NX),VMEAN(NX),WMEAN(NX),DENSMEAN(NX)
C$$$	REAL    :: DENSF_UF_MEAN(NX),DENSRMS(NX),TVMEAN(NX)
C$$$	REAL    :: RWW(NX),RUU(NX),RVV(NX),RWV(NX),RWU(NX)
C$$$	REAL    :: TMP1(NX),TMP2(NX),TMP3(NX),TMP4(NX),TMP5(NX),TMP6(NX),TMP7(NX),TMP8(NX)

C$$$ 	character(len=300) :: filename

C$$$	real,allocatable,dimension(:,:,:) :: Temp_Recv_U,Temp_Recv_V,Temp_Recv_W
C$$$	real,allocatable,dimension(:,:,:) :: Temp_Recv_DENS,Temp_Recv_DENSF,Temp_Recv_TV

C$$$	ALLOCATE(U_TMP(NX,NY,NZ),V_TMP(NX,NY,NZ),W_TMP(NX,NY,NZ))
C$$$	ALLOCATE(DENS_TMP(NX,NY,NZ),DENSF_TMP(NX,NY,NZ))
C$$$ 	CALL CENTER_VELOCITY(NX,NY,NZ,UO,U_TMP,1)
C$$$ 	CALL CENTER_VELOCITY(NX,NY,NZ,VO,V_TMP,2)
C$$$ 	CALL CENTER_VELOCITY(NX,NY,NZ,WO,W_TMP,3)
	
C$$$	ALLOCATE(Temp_Recv_U(2:NX-1,2:NY-1,2:NZ-1),
C$$$     &           Temp_Recv_V(2:NX-1,2:NY-1,2:NZ-1),
C$$$     &           Temp_Recv_W(2:NX-1,2:NY-1,2:NZ-1),
C$$$     &           Temp_Recv_DENS(2:NX-1,2:NY-1,2:NZ-1),
C$$$     &           Temp_Recv_DENSF(2:NX-1,2:NY-1,2:NZ-1),
C$$$     &           Temp_Recv_TV(2:NX-1,2:NY-1,2:NZ-1))

C$$$	if(myrank.eq.0) then
C$$$	   ALLOCATE(U_GLOBAL(2:NX-1,2:NY-1,2:NZG-1),
C$$$     &              V_GLOBAL(2:NX-1,2:NY-1,2:NZG-1),
C$$$     &              W_GLOBAL(2:NX-1,2:NY-1,2:NZG-1),
C$$$     &              DENS_GLOBAL(2:NX-1,2:NY-1,2:NZG-1),
C$$$     &              DENSF_GLOBAL(2:NX-1,2:NY-1,2:NZG-1),
C$$$     &              TV_GLOBAL(2:NX-1,2:NY-1,2:NZG-1) )

C$$$	   U_GLOBAL(2:NX-1,2:NY-1,2:NZ-1) = U_TMP(2:NX-1,2:NY-1,2:NZ-1)
C$$$	   V_GLOBAL(2:NX-1,2:NY-1,2:NZ-1) = V_TMP(2:NX-1,2:NY-1,2:NZ-1)
C$$$	   W_GLOBAL(2:NX-1,2:NY-1,2:NZ-1) = W_TMP(2:NX-1,2:NY-1,2:NZ-1)
C$$$	   DENS_GLOBAL(2:NX-1,2:NY-1,2:NZ-1) = DENS(2:NX-1,2:NY-1,2:NZ-1)
C$$$	   TV_GLOBAL(2:NX-1,2:NY-1,2:NZ-1) = TV(2:NX-1,2:NY-1,2:NZ-1)

C$$$C 	   do i = 2,NX-1
C$$$C 	   do j = 2,NY-1
C$$$C 	   do k = 2,NZ-1
C$$$C 	   DENSF_GLOBAL(i,j,k) = DENS(i,j,k)-dens_bg(i,j)
C$$$C 	   enddo
C$$$C 	   enddo
C$$$C 	   enddo

C$$$	   do n=1,nzprocs-1

C$$$ 	      call MPI_RECV(Temp_Recv_U,(nx-2)*(ny-2)*(nz-2),MPI_DOUBLE_PRECISION,n,1,
C$$$     &                      commx3,status1,stat)
C$$$ 	      call MPI_RECV(Temp_Recv_V,(nx-2)*(ny-2)*(nz-2),MPI_DOUBLE_PRECISION,n,1,
C$$$     &                      commx3,status1,stat)
C$$$ 	      call MPI_RECV(Temp_Recv_W,(nx-2)*(ny-2)*(nz-2),MPI_DOUBLE_PRECISION,n,1,
C$$$     &                      commx3,status1,stat)

C$$$ 	      call MPI_RECV(Temp_Recv_DENS,(nx-2)*(ny-2)*(nz-2),MPI_DOUBLE_PRECISION,n,1,
C$$$     &                      commx3,status1,stat)
C$$$C  	      call MPI_RECV(Temp_Recv_DENSF,(nx-2)*(ny-2)*(nz-2),MPI_DOUBLE_PRECISION,n,1,
C$$$C      &                      commx3,status1,stat)
C$$$ 	      call MPI_RECV(Temp_Recv_TV,(nx-2)*(ny-2)*(nz-2),MPI_DOUBLE_PRECISION,n,1,
C$$$     &                      commx3,status1,stat)


C$$$ 	      kstart=n*(nzg-2)/nzprocs	      

C$$$ 	      U_GLOBAL(2:NX-1,2:NY-1,kstart+2:kstart+nz-1) = Temp_Recv_U(2:NX-1,2:NY-1,2:NZ-1)
C$$$ 	      V_GLOBAL(2:NX-1,2:NY-1,kstart+2:kstart+nz-1) = Temp_Recv_V(2:NX-1,2:NY-1,2:NZ-1)
C$$$ 	      W_GLOBAL(2:NX-1,2:NY-1,kstart+2:kstart+nz-1) = Temp_Recv_W(2:NX-1,2:NY-1,2:NZ-1)

C$$$ 	      DENS_GLOBAL(2:NX-1,2:NY-1,kstart+2:kstart+nz-1) = Temp_Recv_DENS(2:NX-1,2:NY-1,2:NZ-1)
C$$$c 	      DENSF_GLOBAL(2:NX-1,2:NY-1,kstart+2:kstart+nz-1) = Temp_Recv_DENSF(2:NX-1,2:NY-1,2:NZ-1)
C$$$ 	      TV_GLOBAL(2:NX-1,2:NY-1,kstart+2:kstart+nz-1) = Temp_Recv_TV(2:NX-1,2:NY-1,2:NZ-1)

C$$$	   enddo

C$$$	else
	   
C$$$ 	   Temp_Recv_U(2:NX-1,2:NY-1,2:NZ-1) = U_TMP(2:NX-1,2:NY-1,2:NZ-1)
C$$$	   Temp_Recv_V(2:NX-1,2:NY-1,2:NZ-1) = V_TMP(2:NX-1,2:NY-1,2:NZ-1)
C$$$ 	   Temp_Recv_W(2:NX-1,2:NY-1,2:NZ-1) = W_TMP(2:NX-1,2:NY-1,2:NZ-1)

C$$$ 	   Temp_Recv_DENS(2:NX-1,2:NY-1,2:NZ-1) = DENS(2:NX-1,2:NY-1,2:NZ-1)
C$$$ 	   Temp_Recv_TV(2:NX-1,2:NY-1,2:NZ-1) = TV(2:NX-1,2:NY-1,2:NZ-1)

C$$$C 	   do i = 2,NX-1
C$$$C 	   do j = 2,NY-1
C$$$C 	   do k = 2,NZ-1
C$$$C 	   Temp_Recv_DENSF(i,j,k) = DENS(i,j,k)-dens_bg(i,j)
C$$$C 	   enddo
C$$$C 	   enddo
C$$$C 	   enddo

C$$$ 	   call MPI_SEND(Temp_Recv_U,(nx-2)*(ny-2)*(nz-2),MPI_DOUBLE_PRECISION,0,1,
C$$$     &                   commx3,status1,stat)
C$$$ 	   call MPI_SEND(Temp_Recv_V,(nx-2)*(ny-2)*(nz-2),MPI_DOUBLE_PRECISION,0,1,
C$$$     &                   commx3,status1,stat)
C$$$ 	   call MPI_SEND(Temp_Recv_W,(nx-2)*(ny-2)*(nz-2),MPI_DOUBLE_PRECISION,0,1,
C$$$     &                   commx3,status1,stat)


C$$$ 	   call MPI_SEND(Temp_Recv_DENS,(nx-2)*(ny-2)*(nz-2),MPI_DOUBLE_PRECISION,0,1,
C$$$     &                   commx3,status1,stat)
C$$$C  	   call MPI_SEND(Temp_Recv_DENSF,(nx-2)*(ny-2)*(nz-2),MPI_DOUBLE_PRECISION,0,1,
C$$$C      &                   commx3,status1,stat)
C$$$ 	   call MPI_SEND(Temp_Recv_TV,(nx-2)*(ny-2)*(nz-2),MPI_DOUBLE_PRECISION,0,1,
C$$$     &                   commx3,status1,stat)


C$$$	endif
C$$$C       DONE SENDING W,U,V,DENS,DENSF TO CPU1
C$$$C       START CALCULATING WMEAN,REYNOLD STRESS,
C$$$c                         (DENS-DENS_cl)/delta(DENS),
C$$$ 	if(myrank.eq.0) then

C$$$	   DO I=2,NX-1
C$$$	      UMEAN(I)=SUM(U_GLOBAL(I,2:NY-1,2:NZG-1))/
C$$$     &                    (REAL(NY-2)*REAL(NZG-2))
C$$$	      VMEAN(I)=SUM(V_GLOBAL(I,2:NY-1,2:NZG-1))/
C$$$     &                    (REAL(NY-2)*REAL(NZG-2))
C$$$	      WMEAN(I)=SUM(W_GLOBAL(I,2:NY-1,2:NZG-1))/
C$$$     &                    (REAL(NY-2)*REAL(NZG-2))

C$$$	      DENSMEAN(I)=SUM(DENS_GLOBAL(I,2:NY-1,2:NZG-1))/
C$$$     &                    (REAL(NY-2)*REAL(NZG-2))
C$$$	      TVMEAN(I)=SUM(TV_GLOBAL(I,2:NY-1,2:NZG-1))/
C$$$     &                    (REAL(NY-2)*REAL(NZG-2))


C$$$	   ENDDO
	
C$$$ 	   DO I=2,NX-1
C$$$ 	      TMP1(I) = 0.0
C$$$ 	      TMP2(I) = 0.0
C$$$ 	      TMP3(I) = 0.0
C$$$ 	      TMP4(I) = 0.0
C$$$ 	      TMP5(I) = 0.0
C$$$	      TMP6(I) = 0.0
C$$$	      TMP7(I) = 0.0
C$$$ 	      DO K=2,NZG-1
C$$$ 		 DO J=2,NY-1
C$$$ 		    TMP1(I) = TMP1(I)+(W_GLOBAL(I,J,K)-WMEAN(I))*(W_GLOBAL(I,J,K)-WMEAN(I))
C$$$ 		    TMP2(I) = TMP2(I)+(U_GLOBAL(I,J,K)-UMEAN(I))*(U_GLOBAL(I,J,K)-UMEAN(I))
C$$$ 		    TMP3(I) = TMP3(I)+(V_GLOBAL(I,J,K)-VMEAN(I))*(V_GLOBAL(I,J,K)-VMEAN(I))
C$$$ 		    TMP4(I) = TMP4(I)+(W_GLOBAL(I,J,K)-WMEAN(I))*(V_GLOBAL(I,J,K)-VMEAN(I))
C$$$ 		    TMP5(I) = TMP5(I)+(W_GLOBAL(I,J,K)-WMEAN(I))*(U_GLOBAL(I,J,K)-UMEAN(I))

C$$$		    TMP6(I) = TMP6(I)+(DENS_GLOBAL(I,J,K)-DENSMEAN(I))*(U_GLOBAL(I,J,K)-UMEAN(I))
C$$$ 		    TMP7(I) = TMP7(I)+(DENS_GLOBAL(I,J,K)-DENSMEAN(I))*(DENS_GLOBAL(I,J,K)-DENSMEAN(I))

C$$$ 		 ENDDO
C$$$ 	      ENDDO
C$$$ 	      RWW(I) = TMP1(I)/(REAL(NY-2)*REAL(NZG-2))
C$$$ 	      RUU(I) = TMP2(I)/(REAL(NY-2)*REAL(NZG-2))
C$$$ 	      RVV(I) = TMP3(I)/(REAL(NY-2)*REAL(NZG-2))
C$$$	      RWV(I) = TMP4(I)/(REAL(NY-2)*REAL(NZG-2))
C$$$ 	      RWU(I) = TMP5(I)/(REAL(NY-2)*REAL(NZG-2))

C$$$              DENSF_UF_MEAN(I) = TMP6(I)/(REAL(NY-2)*REAL(NZG-2))
C$$$	      DENSRMS(I)=TMP7(I)/(REAL(NY-2)*REAL(NZG-2))
C$$$ 	   ENDDO

C$$$!	   write(filename,'(a,i8.8,a)') "/home/karu/code_Finfty_LES/run/S_LES_ALAMO/B2/w_rs/w_rs_",icycle,".dat"
C$$$c	   write(filename,'(a,i8.8,a)') "/home/karu/code_Finfty_LES/run/S_LES_ALAMO/B1_5/w_rs/w_rs_",icycle,".dat"
C$$$	   write(filename,'(a,i8.8,a)') "/home/karu/code_Finfty_LES/run/US_LES/w_rs/w_rs_",icycle,".dat"

C$$$ 	   open(unit=10,file=filename,form='formatted')
C$$$ 	   DO I=2,NX-1
C$$$ 	      write(10,'(11F15.8)')  xc(I)*(1.0/ru1),WMEAN(I),
C$$$     &	                            RWW(I),RUU(I),RVV(I),RWV(I),RWU(I),
C$$$     &                              DENSMEAN(I), DENSF_UF_MEAN(I), DENSRMS(I),TVMEAN(I)
C$$$ 	   ENDDO
C$$$ 	   close(10)

C$$$!	   write(filename,'(a)') "/home/karu/code_Finfty_LES/run/S_LES_ALAMO/B2/w_rs/time.dat"
C$$$!	   write(filename,'(a)') "/home/karu/code_Finfty_LES/run/S_LES_ALAMO/B1_5/w_rs/time.dat"
C$$$	   write(filename,'(a)') "/home/karu/code_Finfty_LES/run/US_LES/w_rs/time.dat"

C$$$ 	   open(unit=10,file=filename,form='formatted',access='append')
C$$$	   write(10,'(i8.8,F15.8)') icycle, time
C$$$ 	   close(10)

C$$$c$$$c Calculate friction coefficient Cf=2*(u_tau/u_bulk)^2
C$$$c$$$
C$$$c$$$	   CF = 0.0
C$$$c$$$	   DO I=2,NX-1
C$$$c$$$	      CF = CF + WMEAN(I)*(xu(I)-xu(I-1))
C$$$c$$$           ENDDO
C$$$c$$$	   CF = CF/2.0 ! Channel height = 2, NOW CF is u_bulk
C$$$c$$$	   CF = 2.0*(1.0/CF)*(1.0/CF)
C$$$c$$$
C$$$c$$$	   write(filename,'(a)') "/home/karu/code_Finfty_LES/run/S_LES_ALAMO/B2/w_rs/cf.dat"
C$$$c$$$
C$$$c$$$ 	   open(unit=10,file=filename,form='formatted',access='append')
C$$$c$$$	   write(10,'(2F15.8)') time, CF

C$$$ 	endif
	
C$$$	END

C$$$	SUBROUTINE PERIODIC_2D_DISSIPATION_NCPU(NX,NY,NZ,NZG,NBD,ICYCLE,time,DTM1,XC,XU,YC,YV,ZCG,ZWG,xc_car,
C$$$     &                  yc_car,xu_car,yu_car,xv_car,yv_car,P,VO,UO,WO,DENS,TV)

C$$$        use vorticity
C$$$        use density_bg
C$$$	INCLUDE 'common.h'
C$$$        INCLUDE 'mpif.h'
C$$$	INTEGER :: n,NX,NY,NZ,NZG,NBD,ICYCLE
C$$$	REAL	:: DTM1,TIME
C$$$	REAL    :: XU(NX),YV(NY),ZW(NZ),ZWG(NZG)
C$$$	REAL    :: XC(NX),YC(NY),ZC(NZ),ZCG(NZG)
C$$$	REAL    :: UO(NX,NY,NZ),VO(NX,NY,NZ),WO(NX,NY,NZ),P(NX,NY,NZ)
C$$$	REAL    :: DENS(NX,NY,NZ),TV(NX,NY,NZ)
C$$$	REAL	:: XU_CAR(NX,NY),YU_CAR(NX,NY),XV_CAR(NX,NY),YV_CAR(NX,NY),XC_CAR(NX,NY),YC_CAR(NX,NY)
C$$$	INTEGER :: I,J,K,STAT,kstart
C$$$	integer :: status1(MPI_STATUS_SIZE),index2
C$$$        integer iic,ic,ip,im,jjc,jc,jp,jm,jsy,jpsy,kkc,kc,kp,km
C$$$        real dtheta,dr,dz,dq1x3,dq3x1,dq2x3,dq3x2,dq2x1,dq1x2

C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:)::U_TMP,V_TMP,W_TMP,DENS_TMP,DENSF_TMP,TV_TMP
C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:)::U_GLOBAL,V_GLOBAL,W_GLOBAL
C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:)::DENS_GLOBAL,DENSF_GLOBAL,TV_GLOBAL

C$$$	REAL    :: CF

C$$$	REAL    :: UMEAN(NX),VMEAN(NX),WMEAN(NX)
C$$$	REAL    :: EPSILON_RES_MEAN(NX),EPSILON_SGS_MEAN(NX)

C$$$ 	character(len=300) :: filename

C$$$	real,allocatable,dimension(:,:,:) :: Temp_Recv_U,Temp_Recv_V,Temp_Recv_W
C$$$	real,allocatable,dimension(:,:,:) :: Temp_Recv_DENS,Temp_Recv_DENSF,Temp_Recv_TV

C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:):: EPSILON_RES ! resolved dissipation
C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:):: EPSILON_SGS ! subgrid dissipation

C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:):: UF, VF, WF ! fluctuations

C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:):: DU_DX, DU_DY, DU_DZ !derivative of 
C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:):: DV_DX, DV_DY, DV_DZ !velocities
C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:):: DW_DX, DW_DY, DW_DZ 

C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:):: DUF_DX, DUF_DY, DUF_DZ !derivative of 
C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:):: DVF_DX, DVF_DY, DVF_DZ !fluctuating velocities
C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:):: DWF_DX, DWF_DY, DWF_DZ 

C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:):: S11,S12,S13, SF11,SF12,SF13 
C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:):: S21,S22,S23, SF21,SF22,SF23
C$$$	REAL,ALLOCATABLE,DIMENSION(:,:,:):: S31,S32,S33, SF31,SF32,SF33

C$$$	ALLOCATE(Temp_Recv_U(NX,NY,NZ),
C$$$     &           Temp_Recv_V(NX,NY,NZ),
C$$$     &           Temp_Recv_W(NX,NY,NZ),
C$$$     &          Temp_Recv_TV(NX,NY,NZ))

C$$$	if(myrank.eq.0) then
C$$$	   ALLOCATE(U_GLOBAL(NX,NY,NZG),
C$$$     &              V_GLOBAL(NX,NY,NZG),
C$$$     &              W_GLOBAL(NX,NY,NZG),
C$$$     &             TV_GLOBAL(NX,NY,NZG))

C$$$	   U_GLOBAL(1:NX,1:NY,1:NZ) = UO(1:NX,1:NY,1:NZ)
C$$$	   V_GLOBAL(1:NX,1:NY,1:NZ) = VO(1:NX,1:NY,1:NZ)
C$$$	   W_GLOBAL(1:NX,1:NY,1:NZ) = WO(1:NX,1:NY,1:NZ)
C$$$	  TV_GLOBAL(1:NX,1:NY,1:NZ) = TV(1:NX,1:NY,1:NZ)

C$$$	   do n=1,nzprocs-1

C$$$ 	      call MPI_RECV(Temp_Recv_U,nx*ny*nz,MPI_DOUBLE_PRECISION,n,1,
C$$$     &                      commx3,status1,stat)
C$$$ 	      call MPI_RECV(Temp_Recv_V,nx*ny*nz,MPI_DOUBLE_PRECISION,n,1,
C$$$     &                      commx3,status1,stat)
C$$$ 	      call MPI_RECV(Temp_Recv_W,nx*ny*nz,MPI_DOUBLE_PRECISION,n,1,
C$$$     &                      commx3,status1,stat)
C$$$ 	      call MPI_RECV(Temp_Recv_TV,nx*ny*nz,MPI_DOUBLE_PRECISION,n,1,
C$$$     &                      commx3,status1,stat)

C$$$ 	      kstart=n*(nzg-2)/nzprocs	      

C$$$ 	    U_GLOBAL(1:NX,1:NY,kstart+1:kstart+nz)=Temp_Recv_U(1:NX,1:NY,1:NZ)
C$$$ 	    V_GLOBAL(1:NX,1:NY,kstart+1:kstart+nz)=Temp_Recv_V(1:NX,1:NY,1:NZ)
C$$$ 	    W_GLOBAL(1:NX,1:NY,kstart+1:kstart+nz)=Temp_Recv_W(1:NX,1:NY,1:NZ)
C$$$          TV_GLOBAL(1:NX,1:NY,kstart+1:kstart+nz)=Temp_Recv_TV(1:NX,1:NY,1:NZ)

C$$$	   enddo

C$$$	else
	   
C$$$ 	   Temp_Recv_U(1:NX,1:NY,1:NZ) = UO(1:NX,1:NY,1:NZ)
C$$$	   Temp_Recv_V(1:NX,1:NY,1:NZ) = VO(1:NX,1:NY,1:NZ)
C$$$ 	   Temp_Recv_W(1:NX,1:NY,1:NZ) = WO(1:NX,1:NY,1:NZ)
C$$$ 	   Temp_Recv_TV(1:NX,1:NY,1:NZ)= TV(1:NX,1:NY,1:NZ)

C$$$ 	   call MPI_SEND(Temp_Recv_U,nx*ny*nz,MPI_DOUBLE_PRECISION,0,1,
C$$$     &                   commx3,status1,stat)
C$$$ 	   call MPI_SEND(Temp_Recv_V,nx*ny*nz,MPI_DOUBLE_PRECISION,0,1,
C$$$     &                   commx3,status1,stat)
C$$$ 	   call MPI_SEND(Temp_Recv_W,nx*ny*nz,MPI_DOUBLE_PRECISION,0,1,
C$$$     &                   commx3,status1,stat)
C$$$ 	   call MPI_SEND(Temp_Recv_TV,nx*ny*nz,MPI_DOUBLE_PRECISION,0,1,
C$$$     &                   commx3,status1,stat)

C$$$	endif
C$$$C       DONE SENDING W,U,V,DENS,DENSF TO CPU1
C$$$C       START CALCULATING WMEAN,REYNOLD STRESS,
C$$$c                         (DENS-DENS_cl)/delta(DENS),
C$$$ 	if(myrank.eq.0) then

C$$$	   DO I=1,NX
C$$$	      UMEAN(I)=SUM(U_GLOBAL(I,1:NY,1:NZG))/
C$$$     &                    (REAL(NY)*REAL(NZG))
C$$$	      VMEAN(I)=SUM(V_GLOBAL(I,1:NY,1:NZG))/
C$$$     &                    (REAL(NY)*REAL(NZG))
C$$$	      WMEAN(I)=SUM(W_GLOBAL(I,1:NY,1:NZG))/
C$$$     &                    (REAL(NY)*REAL(NZG))

C$$$	   ENDDO

C$$$c---CALCULATE DISSIPATION	
C$$$	ALLOCATE(EPSILON_RES(NX,NY,NZG))
C$$$	ALLOCATE(EPSILON_SGS(NX,NY,NZG))

C$$$	ALLOCATE(UF(NX,NY,NZG),VF(NX,NY,NZG),WF(NX,NY,NZG))

C$$$	ALLOCATE(DU_DX(NX,NY,NZG),DU_DY(NX,NY,NZG),DU_DZ(NX,NY,NZG),
C$$$     &           DV_DX(NX,NY,NZG),DV_DY(NX,NY,NZG),DV_DZ(NX,NY,NZG),
C$$$     &           DW_DX(NX,NY,NZG),DW_DY(NX,NY,NZG),DW_DZ(NX,NY,NZG))
	
C$$$	ALLOCATE(DUF_DX(NX,NY,NZG),DUF_DY(NX,NY,NZG),DUF_DZ(NX,NY,NZG),
C$$$     &           DVF_DX(NX,NY,NZG),DVF_DY(NX,NY,NZG),DVF_DZ(NX,NY,NZG),
C$$$     &           DWF_DX(NX,NY,NZG),DWF_DY(NX,NY,NZG),DWF_DZ(NX,NY,NZG))

C$$$	ALLOCATE(S11(NX,NY,NZG),S12(NX,NY,NZG),S13(NX,NY,NZG),
C$$$     &           S21(NX,NY,NZG),S22(NX,NY,NZG),S23(NX,NY,NZG),
C$$$     &           S31(NX,NY,NZG),S32(NX,NY,NZG),S33(NX,NY,NZG))

C$$$	ALLOCATE(SF11(NX,NY,NZG),SF12(NX,NY,NZG),SF13(NX,NY,NZG),
C$$$     &           SF21(NX,NY,NZG),SF22(NX,NY,NZG),SF23(NX,NY,NZG),
C$$$     &           SF31(NX,NY,NZG),SF32(NX,NY,NZG),SF33(NX,NY,NZG))

C$$$	!------------------------------------------------------!
C$$$	DO I=1,NX
C$$$	DO J=1,NY
C$$$	DO K=1,NZG

C$$$	UF(I,J,K) = U_GLOBAL(I,J,K) - UMEAN(I)
C$$$	VF(I,J,K) = V_GLOBAL(I,J,K) - VMEAN(I)
C$$$	WF(I,J,K) = W_GLOBAL(I,J,K) - WMEAN(I)

C$$$	ENDDO
C$$$	ENDDO
C$$$	ENDDO
C$$$	!------------------------------------------------------!
C$$$        CALL DERIVATIVE_1CPU(U_GLOBAL,DU_DX,NX,NY,NZ,NZG,XC,ZCG,1,1)
C$$$        CALL DERIVATIVE_1CPU(V_GLOBAL,DV_DX,NX,NY,NZ,NZG,XC,ZCG,1,2)
C$$$        CALL DERIVATIVE_1CPU(W_GLOBAL,DW_DX,NX,NY,NZ,NZG,XC,ZCG,1,3)

C$$$        CALL DERIVATIVE_1CPU(U_GLOBAL,DU_DY,NX,NY,NZ,NZG,XC,ZCG,2,1)
C$$$        CALL DERIVATIVE_1CPU(V_GLOBAL,DV_DY,NX,NY,NZ,NZG,XC,ZCG,2,2)
C$$$        CALL DERIVATIVE_1CPU(W_GLOBAL,DW_DY,NX,NY,NZ,NZG,XC,ZCG,2,3)

C$$$        CALL DERIVATIVE_1CPU(U_GLOBAL,DU_DZ,NX,NY,NZ,NZG,XC,ZCG,3,1)
C$$$        CALL DERIVATIVE_1CPU(V_GLOBAL,DV_DZ,NX,NY,NZ,NZG,XC,ZCG,3,2)
C$$$        CALL DERIVATIVE_1CPU(W_GLOBAL,DW_DZ,NX,NY,NZ,NZG,XC,ZCG,3,3)
C$$$	!------------------------------------------------------!
C$$$        CALL DERIVATIVE_1CPU(UF,DUF_DX,NX,NY,NZ,NZG,XC,ZCG,1,1)
C$$$        CALL DERIVATIVE_1CPU(VF,DVF_DX,NX,NY,NZ,NZG,XC,ZCG,1,2)
C$$$        CALL DERIVATIVE_1CPU(WF,DWF_DX,NX,NY,NZ,NZG,XC,ZCG,1,3)

C$$$        CALL DERIVATIVE_1CPU(UF,DUF_DY,NX,NY,NZ,NZG,XC,ZCG,2,1)
C$$$        CALL DERIVATIVE_1CPU(VF,DVF_DY,NX,NY,NZ,NZG,XC,ZCG,2,2)
C$$$        CALL DERIVATIVE_1CPU(WF,DWF_DY,NX,NY,NZ,NZG,XC,ZCG,2,3)

C$$$        CALL DERIVATIVE_1CPU(UF,DUF_DZ,NX,NY,NZ,NZG,XC,ZCG,3,1)
C$$$        CALL DERIVATIVE_1CPU(VF,DVF_DZ,NX,NY,NZ,NZG,XC,ZCG,3,2)
C$$$        CALL DERIVATIVE_1CPU(WF,DWF_DZ,NX,NY,NZ,NZG,XC,ZCG,3,3)
C$$$	!------------------------------------------------------!
C$$$	!------------------------------------------------------!
C$$$	S11 = 0.5*(DU_DX + DU_DX)
C$$$	S12 = 0.5*(DU_DY + DV_DX)
C$$$	S13 = 0.5*(DU_DZ + DW_DX)
	
C$$$	S21 = 0.5*(DV_DX + DU_DY)
C$$$	S22 = 0.5*(DV_DY + DV_DY)
C$$$	S23 = 0.5*(DV_DZ + DW_DY)

C$$$	S31 = 0.5*(DW_DX + DU_DZ)
C$$$	S32 = 0.5*(DW_DY + DV_DZ)
C$$$	S33 = 0.5*(DW_DZ + DW_DZ)
C$$$	!------------------------------------------------------!
C$$$	SF11 = 0.5*(DUF_DX + DUF_DX)
C$$$	SF12 = 0.5*(DUF_DY + DVF_DX)
C$$$	SF13 = 0.5*(DUF_DZ + DWF_DX)
	
C$$$	SF21 = 0.5*(DVF_DX + DUF_DY)
C$$$	SF22 = 0.5*(DVF_DY + DVF_DY)
C$$$	SF23 = 0.5*(DVF_DZ + DWF_DY)

C$$$	SF31 = 0.5*(DWF_DX + DUF_DZ)
C$$$	SF32 = 0.5*(DWF_DY + DVF_DZ)
C$$$	SF33 = 0.5*(DWF_DZ + DWF_DZ)

C$$$	DEALLOCATE(DUF_DX,DUF_DY,DUF_DZ)
C$$$	DEALLOCATE(DVF_DX,DVF_DY,DVF_DZ)
C$$$	DEALLOCATE(DWF_DX,DWF_DY,DWF_DZ)
C$$$	!------------------------------------------------------!
C$$$	! COMPUTE RESOLVED DISSIPATION
C$$$	! EPSILON_RES = 2 \NU <SF_IJ*SF_IJ>
C$$$	EPSILON_RES = 2.0*ru1*( SF11*SF11 + SF12*SF12 + SF13*SF13 +
C$$$     &                          SF21*SF21 + SF22*SF22 + SF23*SF23 + 
C$$$     &                          SF31*SF31 + SF32*SF32 + SF33*SF33  )
C$$$	DO I=2,NX-1
C$$$	   EPSILON_RES_MEAN(I)=SUM(EPSILON_RES(I,2:NY-1,2:NZG-1))/
C$$$     &                                        (REAL(NY-2)*REAL(NZG-2))
C$$$	ENDDO

C$$$	! COMPUTE SGS DISSIPATION
C$$$	! EPSILON_SGS = -<-2 \NU_EDDY S_IJ>*SFIJ
C$$$	EPSILON_SGS = 2.0*TV_GLOBAL*( S11*SF11 + S12*SF12 + S13*SF13 +
C$$$     &                         S21*SF21 + S22*SF22 + S23*SF23 + 
C$$$     &                         S31*SF31 + S32*SF32 + S33*SF33  )

C$$$	DO I=2,NX-1
C$$$	   EPSILON_SGS_MEAN(I)=SUM(EPSILON_SGS(I,2:NY-1,2:NZG-1))/
C$$$     &                                        (REAL(NY-2)*REAL(NZG-2))
C$$$	ENDDO

C$$$!	write(filename,'(a,i8.8,a)') "/home/karu/code_Finfty_LES/run/US_LES/epsilon/eps_",icycle,".dat"
C$$$!	write(filename,'(a,i8.8,a)') "/home/karu/code_Finfty_LES/run/S_LES_ALAMO/B1_3/epsilon/eps_",icycle,".dat"
C$$$	write(filename,'(a,i8.8,a)') "/home/karu/code_Finfty_LES/run/S_LES_ALAMO/B2/epsilon/eps_",icycle,".dat"
	
C$$$	open(unit=10,file=filename,form='formatted')
C$$$	DO I=2,NX-1
C$$$C 	   write(10,'(4(4x,E15.8))')  xc(I),EPSILON_RES_MEAN(I),EPSILON_SGS_MEAN(I),
C$$$C      &                                 EPSILON_RES_MEAN(I)+EPSILON_SGS_MEAN(I)
C$$$	   write(10,'(4F15.8)')  xc(I),EPSILON_RES_MEAN(I),EPSILON_SGS_MEAN(I),
C$$$     &                                 EPSILON_RES_MEAN(I)+EPSILON_SGS_MEAN(I)

C$$$	ENDDO
C$$$	close(10)
	
C$$$c	write(filename,'(a)') "/home/karu/code_Finfty_LES/run/US_LES/epsilon/time.dat"
C$$$c	write(filename,'(a)') "/home/karu/code_Finfty_LES/run/S_LES_ALAMO/B1_3/epsilon/time.dat"
C$$$	write(filename,'(a)') "/home/karu/code_Finfty_LES/run/S_LES_ALAMO/B2/epsilon/time.dat"
	
C$$$	open(unit=10,file=filename,form='formatted',access='append')
C$$$	write(10,'(i8.8,F15.8)') icycle, time
C$$$	close(10)
	
C$$$ 	endif
	
C$$$	END
	
