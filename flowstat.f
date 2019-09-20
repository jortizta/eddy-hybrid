      SUBROUTINE STAT2D(UO,VO,WO,P,TV,NAME,NX,NY,NZ,TIME)

      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
      CHARACTER NAME*(*)
      INTEGER NX,NY,NZ,STATUS(MPI_STATUS_SIZE)

      REAL UO(NX,NY,NZ),VO(NX,NY,NZ),WO(NX,NY,NZ)
      REAL P(NX,NY,NZ),TV(NX,NY,NZ)
      REAL TIME
      REAL, DIMENSION(:,:), ALLOCATABLE :: UMEAN,VMEAN,WMEAN,PMEAN,TVMEAN
      REAL, DIMENSION(:,:), ALLOCATABLE :: UU,VV,WW,UV,UW,VW,PP
      INTEGER I,K,JP
      REAL SCLY

      ALLOCATE(UMEAN(NX,NZ),VMEAN(NX,NZ),WMEAN(NX,NZ),PMEAN(NX,NZ),TVMEAN(NX,NZ))
      ALLOCATE(UU(NX,NZ),VV(NX,NZ),WW(NX,NZ),UV(NX,NZ),UW(NX,NZ),VW(NX,NZ),PP(NX,NZ))

      UMEAN = 0.
      VMEAN = 0.
      WMEAN = 0.
      PMEAN = 0.
      TVMEAN = 0.
      UU(:,:) = 0.
      VV(:,:) = 0.
      WW(:,:) = 0.
      UV(:,:) = 0.
      UW(:,:) = 0.
      VW(:,:) = 0.
      PP(:,:) = 0.

      SCLY = 1./REAL(NY-2)

      IF(MYRANK==0) WRITE(6,*) ' writing stat at t= ',TIME
c
C     COMPUTE THE STATISTICS. 
c
      DO K=KZ1,KZ2
      DO I=IX1,IX2
        UMEAN(I,K) = SUM(UO(I,JY1:JY2,K)+UO(I-1,JY1  :JY2  ,K  ))
        VMEAN(I,K) = SUM(VO(I,JY1:JY2,K)+VO(I  ,JY1-1:JY2-1,K  ))
        WMEAN(I,K) = SUM(WO(I,JY1:JY2,K)+WO(I  ,JY1  :JY2  ,K-1))
        PMEAN(I,K) =  SUM(P (I,JY1:JY2,K))
        TVMEAN(I,K) = SUM(TV(I,JY1:JY2,K))
        UU(I,K) = SUM((UO(I,JY1:JY2,K)+UO(I-1,JY1  :JY2  ,K  ))**2)
        VV(I,K) = SUM((VO(I,JY1:JY2,K)+VO(I  ,JY1-1:JY2-1,K  ))**2)
        WW(I,K) = SUM((WO(I,JY1:JY2,K)+WO(I  ,JY1  :JY2  ,K-1))**2)
        UV(I,K) = SUM((UO(I,JY1:JY2,K)+UO(I-1,JY1  :JY2  ,K  ))
     &       *        (VO(I,JY1:JY2,K)+VO(I  ,JY1-1:JY2-1,K  )))
        UW(I,K) = SUM((UO(I,JY1:JY2,K)+UO(I-1,JY1  :JY2  ,K  ))
     &       *        (WO(I,JY1:JY2,K)+WO(I  ,JY1  :JY2  ,K-1)))
        VW(I,K) = SUM((VO(I,JY1:JY2,K)+VO(I  ,JY1-1:JY2-1,K  ))
     &       *        (WO(I,JY1:JY2,K)+WO(I  ,JY1  :JY2  ,K-1)))
        PP(I,K) = SUM(P (I,JY1:JY2,K)**2)
      ENDDO
      ENDDO
      UMEAN(:,:) = UMEAN(:,:)*SCLY*0.5
      VMEAN(:,:) = VMEAN(:,:)*SCLY*0.5
      WMEAN(:,:) = WMEAN(:,:)*SCLY*0.5
      PMEAN(:,:) = PMEAN(:,:)*SCLY
      TVMEAN(:,:) = TVMEAN(:,:)*SCLY
      UU(:,:) = UU(:,:)*SCLY*0.25
      VV(:,:) = VV(:,:)*SCLY*0.25
      WW(:,:) = WW(:,:)*SCLY*0.25
      UV(:,:) = UV(:,:)*SCLY*0.25
      UW(:,:) = UW(:,:)*SCLY*0.25
      VW(:,:) = VW(:,:)*SCLY*0.25
      PP(:,:) = PP(:,:)*SCLY

C     ROOT PROCESSOR WRITES THE STATISTICS
      IF(MYRANK==0) THEN
        OPEN(19,FILE=NAME,FORM='unformatted')
c        WRITE(19) TIME
        WRITE(19) NX,1,MYSIZE*(NZ-2)+2,12
        DO K=1,KZ2+1/MYSIZE     ! write my own data
          WRITE(19) (UMEAN(I,K),I=1,NX)
     &         ,    (VMEAN(I,K),I=1,NX)
     &         ,    (WMEAN(I,K),I=1,NX)
     &         ,    (PMEAN(I,K),I=1,NX)
     &         ,    (TVMEAN(I,K),I=1,NX)
     &         ,    (UU(I,K),I=1,NX)
     &         ,    (VV(I,K),I=1,NX)
     &         ,    (WW(I,K),I=1,NX)
     &         ,    (UW(I,K),I=1,NX)
     &         ,    (UV(I,K),I=1,NX)
     &         ,    (VW(I,K),I=1,NX)
     &         ,    (PP(I,K),I=1,NX)
        ENDDO
C.....receive data for others
        DO JP=1,MYSIZE-1
          CALL MPI_RECV(UMEAN(1,1),NX*NZ,MTYPE,JP,0,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(VMEAN(1,1),NX*NZ,MTYPE,JP,0,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(WMEAN(1,1),NX*NZ,MTYPE,JP,0,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(PMEAN(1,1),NX*NZ,MTYPE,JP,0,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(TVMEAN(1,1),NX*NZ,MTYPE,JP,0,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(UU(1,1),NX*NZ,MTYPE,JP,0,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(VV(1,1),NX*NZ,MTYPE,JP,0,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(WW(1,1),NX*NZ,MTYPE,JP,0,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(UV(1,1),NX*NZ,MTYPE,JP,0,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(UW(1,1),NX*NZ,MTYPE,JP,0,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(VW(1,1),NX*NZ,MTYPE,JP,0,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(PP(1,1),NX*NZ,MTYPE,JP,0,MPI_COMM_EDDY,STATUS,IERR)
C.....write others data
          DO K=KZ1,KZ2+(JP+1)/MYSIZE
            WRITE(19) (UMEAN(I,K),I=1,NX)
     &           ,    (VMEAN(I,K),I=1,NX)
     &           ,    (WMEAN(I,K),I=1,NX)
     &           ,    (PMEAN(I,K),I=1,NX)
     &           ,    (TVMEAN(I,K),I=1,NX)
     &           ,    (UU(I,K),I=1,NX)
     &           ,    (VV(I,K),I=1,NX)
     &           ,    (WW(I,K),I=1,NX)
     &           ,    (UW(I,K),I=1,NX)
     &           ,    (UV(I,K),I=1,NX)
     &           ,    (VW(I,K),I=1,NX)
     &           ,    (PP(I,K),I=1,NX)
          ENDDO
        ENDDO
        CLOSE(19)
C.....send data to root
      ELSE
        CALL MPI_SEND(UMEAN(1,1),NX*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND(VMEAN(1,1),NX*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND(WMEAN(1,1),NX*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND(PMEAN(1,1),NX*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND(TVMEAN(1,1),NX*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND(UU(1,1),NX*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND(VV(1,1),NX*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND(WW(1,1),NX*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND(UV(1,1),NX*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND(UW(1,1),NX*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND(VW(1,1),NX*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND(PP(1,1),NX*NZ,MTYPE,0,0,MPI_COMM_EDDY,IERR)
      ENDIF
C
      DEALLOCATE(UMEAN,VMEAN,WMEAN,PMEAN,TVMEAN)
      DEALLOCATE(UU,VV,WW,UV,UW,VW,PP)
C
      RETURN
      END


      SUBROUTINE STAT1D(UO,VO,WO,P,TV,NAME,NX,NY,NZ,NZG,TIME)

      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
      CHARACTER NAME*(*)
      INTEGER NX,NY,NZ,NZG,STATUS(MPI_STATUS_SIZE)

      REAL UO(NX,NY,NZ),VO(NX,NY,NZ),WO(NX,NY,NZ)
      REAL P(NX,NY,NZ),TV(NX,NY,NZ)
      REAL TIME
      INTEGER I,J

      REAL, DIMENSION(:,:), ALLOCATABLE :: STAT,STATG

      ALLOCATE(STAT(NX,12),STATG(NX,12))

c the centered values are added in STAT
      DO I=IX1,IX2
        STAT(I,1) = SUM(UO(I,JY1:JY2,KZ1:KZ2)+UO(I-1,JY1:JY2  ,KZ1:KZ2))
        STAT(I,2) = SUM(VO(I,JY1:JY2,KZ1:KZ2)+VO(I,JY1-1:JY2-1,KZ1:KZ2))
        STAT(I,3) = SUM(WO(I,JY1:JY2,KZ1:KZ2)+WO(I,JY1:JY2,KZ1-1:KZ2-1))
        STAT(I,4) = SUM(P (I,JY1:JY2,KZ1:KZ2))
        STAT(I,5) = SUM(TV(I,JY1:JY2,KZ1:KZ2))
        STAT(I,6) = SUM((UO(I,JY1:JY2,KZ1:KZ2)+UO(I-1,JY1:JY2  ,KZ1:KZ2))**2)
        STAT(I,7) = SUM((VO(I,JY1:JY2,KZ1:KZ2)+VO(I,JY1-1:JY2-1,KZ1:KZ2))**2)
        STAT(I,8) = SUM((WO(I,JY1:JY2,KZ1:KZ2)+WO(I,JY1:JY2,KZ1-1:KZ2-1))**2)
        STAT(I,9) = SUM((UO(I,JY1:JY2,KZ1:KZ2)+UO(I-1,JY1:JY2  ,KZ1:KZ2))
     &       *          (VO(I,JY1:JY2,KZ1:KZ2)+VO(I,JY1-1:JY2-1,KZ1:KZ2)))
        STAT(I,10) = SUM((VO(I,JY1:JY2,KZ1:KZ2)+VO(I,JY1-1:JY2-1,KZ1:KZ2))
     &       *           (WO(I,JY1:JY2,KZ1:KZ2)+WO(I,JY1:JY2,KZ1-1:KZ2-1)))
        STAT(I,11) = SUM((UO(I,JY1:JY2,KZ1:KZ2)+UO(I-1,JY1:JY2  ,KZ1:KZ2))
     &       *           (WO(I,JY1:JY2,KZ1:KZ2)+WO(I,JY1:JY2,KZ1-1:KZ2-1)))
        STAT(I,12) = SUM( P (I,JY1:JY2,KZ1:KZ2)**2)
!!!!!!        STAT(I,3) = SUM(WO(I,JY1:JY2,KZ1:KZ2))
      ENDDO

c The global sums are evaluated at the 0 processor
      CALL MPI_REDUCE(STAT,STATG,12*NX,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)

      IF(MYRANK==0) THEN

c The averages are established
        STATG(:,1:3 ) = STATG(:,1:3 )/REAL(2*(NY-2)*(NZG-2))
        STATG(:,6:11) = STATG(:,6:11)/REAL(4*(NY-2)*(NZG-2))

c The 0 processor writes the averaged values
        OPEN(19,FILE='stat1d.'//INDEX(NWRITE),FORM='unformatted')
        WRITE(19) NX-2,12
        WRITE(19) ((STATG(I,J),I=IX1,IX2),J=1,12)
        CLOSE(19)

      ENDIF

      DEALLOCATE(STAT,STATG)
C
      RETURN
      END

c
c--------------------------------------------------------------------------


C---- read_stat1d_input ------------------------N. Beratlis-19 Sep. 2011---
C
C     PURPOSE: Read parameters from stats1d.input file
C
c--------------------------------------------------------------------------
      subroutine read_stat1d_input(flag,it)
c
      IMPLICIT NONE
c
c.... Input/Output arrays
      INTEGER flag,it
      
      open(unit=10,file='stat1D.input',form='formatted')
      read(10,'(A)')
      read(10,*) flag
      read(10,'(A)')
      read(10,*) it
      close(10)

      return

      end
c--------------------------------------------------------------------------



C---- read_stat1d_ave --------------------------N. Beratlis-19 Sep. 2011---
C
C     PURPOSE: Read stats1d
C
c--------------------------------------------------------------------------
      subroutine read_stat1d_ave(stat,nx,nvar,nsmpl)
c
      INCLUDE 'common.h'
c
c.... Input/Output arrays
      integer nx,nvar,nsmpl
      real    stat(nx,nvar)
c
c.... Local arrays
      integer i,j
      real tmp
      
      if(myrank.eq.0) then
        open(unit=10,file='stat1dave.bin',form='unformatted')
        read(10) nsmpl
        DO I=1,NX
          read(10) TMP,(STAT(I,J),J=1,NVAR)
        ENDDO
        close(10)
      endif

      return
      end
c--------------------------------------------------------------------------


C---- subroutine calc_stat1d_ave ---------------N. Beratlis-19 Sep. 2011---
C
C     PURPOSE: Compute average 1D stat (flow homogenenous in y and z directions,
C      i.e channel, pipe flow)
C
c--------------------------------------------------------------------------
      SUBROUTINE calc_stat1d_ave(UO,VO,WO,P,TV,STAT,NVAR,NX,NY,NZ,NZG,NSMPL)
c
      include 'common.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,nzg,nsmpl,nvar
      real stat(nx,nvar)
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz),tv(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k
      real    inv,s,invs,tmp

      inv = 1./real(nz-2)/real(ny-2)
      s = real(nsmpl)
      invs = 1./(s+1.0)
      
      DO I=IX1,NX
        STAT(I,1) = (STAT(I,1)*s + SUM(UO(I,JY1:JY2,KZ1:KZ2)+UO(I-1,JY1:JY2  ,KZ1:KZ2))*inv*0.5)*invs
      ENDDO
      DO I=1,NX
        STAT(I,2) = (STAT(I,2)*s + SUM(VO(I,JY1:JY2,KZ1:KZ2)+VO(I,JY1-1:JY2-1,KZ1:KZ2))*inv*0.5)*invs
        STAT(I,3) = (STAT(I,3)*s + SUM(WO(I,JY1:JY2,KZ1:KZ2)+WO(I,JY1:JY2,KZ1-1:KZ2-1))*inv*0.5)*invs
        STAT(I,4) = (STAT(I,4)*s + SUM(P (I,JY1:JY2,KZ1:KZ2))*inv)*invs
        STAT(I,5) = (STAT(I,5)*s + SUM(TV(I,JY1:JY2,KZ1:KZ2))*inv)*invs
      ENDDO
      DO I=IX1,NX
        STAT(I,6) = (STAT(I,6)*s + SUM((UO(I,JY1:JY2,KZ1:KZ2)+UO(I-1,JY1:JY2  ,KZ1:KZ2))**2)*inv*0.25)*invs
      ENDDO
      DO I=1,NX
        STAT(I,7) = (STAT(I,7)*s + SUM((VO(I,JY1:JY2,KZ1:KZ2)+VO(I,JY1-1:JY2-1,KZ1:KZ2))**2)*inv*0.25)*invs
        STAT(I,8) = (STAT(I,8)*s + SUM((WO(I,JY1:JY2,KZ1:KZ2)+WO(I,JY1:JY2,KZ1-1:KZ2-1))**2)*inv*0.25)*invs
      ENDDO
      DO I=IX1,NX
        STAT(I,9) = (STAT(I,9)*s + SUM((UO(I,JY1:JY2,KZ1:KZ2)+UO(I-1,JY1:JY2  ,KZ1:KZ2))
     &       *          (VO(I,JY1:JY2,KZ1:KZ2)+VO(I,JY1-1:JY2-1,KZ1:KZ2)))*inv*0.25)*invs
        STAT(I,10) = (STAT(I,10)*s + SUM((UO(I,JY1:JY2,KZ1:KZ2)+UO(I-1,JY1:JY2  ,KZ1:KZ2))
     &       *           (WO(I,JY1:JY2,KZ1:KZ2)+WO(I,JY1:JY2,KZ1-1:KZ2-1)))*inv*0.25)*invs
      ENDDO
      DO I=1,NX
        STAT(I,11) = (STAT(I,11)*s + SUM((VO(I,JY1:JY2,KZ1:KZ2)+VO(I,JY1-1:JY2-1,KZ1:KZ2))
     &       *           (WO(I,JY1:JY2,KZ1:KZ2)+WO(I,JY1:JY2,KZ1-1:KZ2-1)))*inv*0.25)*invs
        STAT(I,12) = (STAT(I,12)*s + SUM( P (I,JY1:JY2,KZ1:KZ2)**2)*inv)*invs
        STAT(I,13) = (STAT(I,13)*s + dpdz)*invs
      ENDDO
      DO I=IX1,NX
        tmp = SUM( (AP(I)*(UO(I,JY1:JY2,KZ1:KZ2)-UO(I-1,JY1:JY2,KZ1:KZ2)))**2.0)*inv
        STAT(I,14) = (STAT(I,14)*s + tmp)*invs

        tmp = 0.0
        DO J=JY1,JY2
          tmp = tmp + SUM((( 1.0*(BU(J)*(UO( I ,J+1,KZ1:KZ2)-UO( I ,J-1,KZ1:KZ2)))
     &                   + 1.0*(BU(J)*(UO(I-1,J+1,KZ1:KZ2)-UO(I-1,J-1,KZ1:KZ2)))
     &                   - 0.5*(VO(I,J,KZ1:KZ2)+VO(I,J-1,KZ1:KZ2))*icyl )/rp(i))**2.0)
        ENDDO
        tmp = tmp*inv
        STAT(I,15) = (STAT(I,15)*s + tmp)*invs

        tmp = 0.0
        DO K=KZ1,KZ2
          tmp = tmp + SUM(( 1.0*(CU(K)*(UO( I ,JY1:JY2,K+1)-UO( I ,JY1:JY2,K-1)))
     &                    + 1.0*(CU(K)*(UO(I-1,JY1:JY2,K+1)-UO(I-1,JY1:JY2,K-1))) )**2.0)
        ENDDO
        tmp = tmp*inv
        STAT(I,16) = (STAT(I,16)*s + tmp)*invs
      ENDDO

      DO I=IX1,IX2
        tmp = SUM(( 1.0*(AU(I)*(VO(I+1,JY1:JY2,KZ1:KZ2)-VO(I-1,JY1:JY2,KZ1:KZ2)))
     &            + 1.0*(AU(I)*(VO(I+1,JY1-1:JY2-1,KZ1:KZ2)-VO(I-1,JY1-1:JY2-1,KZ1:KZ2))) )**2.0)*inv
        STAT(I,17) = (STAT(I,17)*s + tmp)*invs
      ENDDO

      DO I=IX1,NX
        tmp = 0.0
        DO J=JY1,JY2
          tmp = tmp + SUM( (BP(J)*(VO(I,J,KZ1:KZ2)-VO(I,J-1,KZ1:KZ2))
     &                   + icyl*0.5*(UO(I,J,KZ1:KZ2)+UO(I-1,J,KZ1:KZ2))/rp(i))**2.0)
        ENDDO
        tmp = tmp*inv
        STAT(I,18) = (STAT(I,18)*s + tmp)*invs
      ENDDO

      DO I=1,NX
        tmp = 0.0
        DO K=KZ1,KZ2
        tmp = tmp + SUM(( 1.0*(CV(K)*(VO(I,JY1:JY2,K+1)-VO(I,JY1:JY2,K-1)))
     &                  + 1.0*(CV(K)*(VO(I,JY1-1:JY2-1,K+1)-VO(I,JY1-1:JY2-1,K-1))) )**2.0)
        ENDDO
        tmp = tmp*inv
        STAT(I,19) = (STAT(I,19)*s + tmp)*invs
      ENDDO

      DO I=IX1,IX2
        tmp = SUM(( 1.0*(AW(I)*(WO(I+1,JY1:JY2, KZ1 : KZ2 )-WO(I-1,JY1:JY2, KZ1 : KZ2 )))
     &            + 1.0*(AW(I)*(WO(I+1,JY1:JY2,KZ1-1:KZ2-1)-WO(I-1,JY1:JY2,KZ1-1:KZ2-1))))**2.0)*inv
        STAT(I,20) = (STAT(I,20)*s + tmp)*invs
      ENDDO

      DO I=1,NX
        tmp = 0.0
        DO J=JY1,JY2
          tmp = tmp + SUM((( 1.0*(BW(J)*(WO(I,J+1, KZ1 : KZ2 )-WO(I,J-1, KZ1 : KZ2 )))
     &                     + 1.0*(BW(J)*(WO(I,J+1,KZ1-1:KZ2-1)-WO(I,J-1,KZ1-1:KZ2-1))))/rp(i))**2.0)
        ENDDO
        tmp = tmp*inv
        STAT(I,21) = (STAT(I,21)*s + tmp)*invs

        tmp = 0.0
        DO K=KZ1,KZ2
          tmp = tmp + SUM( (CP(K)*(WO(I,JY1:JY2,K)-WO(I,JY1:JY2,K-1)))**2.0)
        ENDDO
        tmp = tmp*inv
        STAT(I,22) = (STAT(I,22)*s + tmp)*invs
      ENDDO

      DO I=IX1,IX2
c        tmp = SUM(WO(I,JY1:JY2,KZ1:KZ2))*inv
        tmp = SUM(( 1.0*(AW(I)*(WO(I+1,JY1:JY2, KZ1 : KZ2 )-WO(I-1,JY1:JY2, KZ1 : KZ2 )))
     &            + 1.0*(AW(I)*(WO(I+1,JY1:JY2,KZ1-1:KZ2-1)-WO(I-1,JY1:JY2,KZ1-1:KZ2-1)))))*inv*0.5
        STAT(I,23) = (STAT(I,23)*s + tmp)*invs
      ENDDO

c      STAT(1,:) = STAT(2,:)
c      STAT(1,1) = -STAT(2,1)
      
      RETURN

      END
c--------------------------------------------------------------------------


C---- subroutine calc_stat1d_ave_les -----------N. Beratlis-19 Sep. 2011---
C
C     PURPOSE: Compute average 1D stat (flow homogenenous in y and z directions,
C      i.e channel, pipe flow)
C
c--------------------------------------------------------------------------
      SUBROUTINE calc_stat1d_ave_les(CLES,CLESP,CLESN,NILMP,NILMN,ILM,IMM
     &     ,G,SXX,SYY,SZZ,SXY,SXZ,SYZ,STAT,NVAR,NX,NY,NZ,NZG,NSMPL)
c
      include 'common.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,nzg,nsmpl,nvar
      real stat(nx,nvar)
      real cles(nx,ny,nz),ilm(nx,ny,nz),imm(nx,ny,nz),g(nx,ny,nz)
      real clesp(nx,ny,nz),clesn(nx,ny,nz)
      integer nilmp(nx,ny,nz),nilmn(nx,ny,nz)
      real sxx(nx,ny,nz),syy(nx,ny,nz),szz(nx,ny,nz)
     &    ,sxy(nx,ny,nz),sxz(nx,ny,nz),syz(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k
      real    inv,s,invs,tmp

      inv = 1./real(nz-2)/real(ny-2)
      s = real(nsmpl)
      invs = 1./(s+1.0)

      DO I=1,NX
        STAT(I,24) = (STAT(I,24)*s + SUM(CLES(I,JY1:JY2,KZ1:KZ2))*inv)*invs
        STAT(I,25) = (STAT(I,25)*s + SUM(ILM(I,JY1:JY2,KZ1:KZ2))*inv)*invs
        STAT(I,26) = (STAT(I,26)*s + SUM(IMM(I,JY1:JY2,KZ1:KZ2))*inv)*invs
        STAT(I,27) = (STAT(I,27)*s + SUM(G(I,JY1:JY2,KZ1:KZ2))*inv)*invs
        STAT(I,28) = (STAT(I,28)*s + SUM(CLESP(I,JY1:JY2,KZ1:KZ2))*inv)*invs
        STAT(I,29) = (STAT(I,29)*s + SUM(CLESN(I,JY1:JY2,KZ1:KZ2))*inv)*invs
        STAT(I,30) = (STAT(I,30)*s + SUM(NILMP(I,JY1:JY2,KZ1:KZ2))*inv)*invs
        STAT(I,31) = (STAT(I,31)*s + SUM(NILMN(I,JY1:JY2,KZ1:KZ2))*inv)*invs
        STAT(I,32) = (STAT(I,32)*s + SUM(ABS(SXX(I,JY1:JY2,KZ1:KZ2)))*inv)*invs
        STAT(I,33) = (STAT(I,33)*s + SUM(ABS(SYY(I,JY1:JY2,KZ1:KZ2)))*inv)*invs
        STAT(I,34) = (STAT(I,34)*s + SUM(ABS(SZZ(I,JY1:JY2,KZ1:KZ2)))*inv)*invs
        STAT(I,35) = (STAT(I,35)*s + SUM(ABS(SXY(I,JY1:JY2,KZ1:KZ2)))*inv)*invs
        STAT(I,36) = (STAT(I,36)*s + SUM(ABS(SXZ(I,JY1:JY2,KZ1:KZ2)))*inv)*invs
        STAT(I,37) = (STAT(I,37)*s + SUM(ABS(SYZ(I,JY1:JY2,KZ1:KZ2)))*inv)*invs
      ENDDO
      
      RETURN

      END
c--------------------------------------------------------------------------


C---- write_stat1d_ave -------------------------N. Beratlis-19 Sep. 2011---
C
C     PURPOSE: Write stats1d
C
c--------------------------------------------------------------------------
      subroutine write_stat1d_ave(stat1Dave,stat1D,xc,nx,nvar,nstat1Dave,nstat1D)
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
c
c.... Input/Output arrays
      integer nx,nstat1Dave,nstat1D,nvar
      real    stat1Dave(nx,nvar),stat1D(nx,nvar),xc(nx)
c
c.... Local arrays
      integer i,j
      real    n1,n2,n3,ppp
      real    stat1DG(nx,nvar),dwdx(nx)

      PPP = 2./3.

      IF(MYSIZE==1) THEN
        stat1dg = stat1d
      ELSE
        CALL MPI_REDUCE(stat1D,stat1DG,NVAR*NX,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      ENDIF

      if(myrank.eq.0) then
        n1 = real(nstat1Dave)
        n2 = real(nstat1D)
        n3 = real(nstat1Dave+nstat1D)
        stat1Dave = (stat1Dave*n1+stat1D*n2)/n3
        nstat1Dave = nstat1Dave + nstat1D

        do i=ix1,ix2
          stat1Dave(i,23) = aw(i)*(stat1dave(i+1,3)-stat1dave(i-1,3))
        enddo

        open(unit=10,file='stat1dave.bin',form='unformatted')
        write(10) nstat1Dave+1
        write(6,*) nstat1Dave,nvar
        DO I=1,NX
          write(10) XC(I),(STAT1DAVE(I,J),J=1,NVAR),(ABS(RP(I))/(AP(I)*BP(2)*CP(2)))**PPP
        ENDDO
        close(10)

        open(unit=10,file='stat1dave.plt',form='formatted')
        IF(ISGS==0 .OR. ISGS==5 .OR. ISGS==6) THEN
          write(10,'(A)') 'VARIABLES = "X", "U", "V", "W", "P", "TV", '
     &        //'"UU", "VV", "WW", "UV", "UW", "VW", "PP", "DPDZ", '
     &        //'"DUDXSQ", "DUDYSQ", "DUDZSQ", "DVDXSQ", "DVDYSQ", '
     &        //'"DVDZSQ", "DWDXSQ", "DWDYSQ", "DWDZSQ", "DWDX", '
     &        //'"DELTA2"'
        ELSE
          write(10,'(A)') 'VARIABLES = "X", "U", "V", "W", "P", "TV", '
     &        //'"UU", "VV", "WW", "UV", "UW", "VW", "PP", "DPDZ", '
     &        //'"DUDXSQ", "DUDYSQ", "DUDZSQ", "DVDXSQ", "DVDYSQ", '
     &        //'"DVDZSQ", "DWDXSQ", "DWDYSQ", "DWDZSQ", "DWDX" '
     &        //', "CLES", "ILM", "IMM", "G"'
     &        //', "CLESP", "CLESN", "NILMP", "NILMN"'
     &        //', "SXX", "SYY", "SZZ", "SXY", "SXZ", "SYZ"'
     &        //', "DELTA2"'
        ENDIF

        DO I=1,NX
          write(10,'(63(1x,E19.12))') XC(I),(STAT1DAVE(I,J),J=1,NVAR),(ABS(RP(I))/(AP(I)*BP(2)*CP(2)))**PPP
        ENDDO
        
        close(10)
      endif

      return
      end
c--------------------------------------------------------------------------



      SUBROUTINE TIMESERIES(UO,VO,WO,P,TV,NAME,NX,NY,NZ,TIME)

      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
      CHARACTER NAME*(*)
      INTEGER NX,NY,NZ,STATUS(MPI_STATUS_SIZE)

      REAL UO(NX,NY,NZ),VO(NX,NY,NZ),WO(NX,NY,NZ)
      REAL P(NX,NY,NZ),TV(NX,NY,NZ)
      REAL TIME
      INTEGER I,J,K
      REAL    AREANZ,AREAIJ
      REAL    QSTAR,UBUK,UBUKG,TAUW,TAUWT,TAUWB,TAUWG,UBAR,UBARG

        QSTAR = 0.0
        AREANZ = 0.
        DO I=IX1,IX2
        DO J=JY1,JY2
          AREAIJ = RP(I)/(AP(I)*BP(J))
          DO K=KZ1,KZ2
            QSTAR = QSTAR+WO(I,J,K)*AREAIJ
            AREANZ = AREANZ+AREAIJ
          ENDDO
        ENDDO
        ENDDO
  
        UBUK = QSTAR/AREANZ

        IF(ICYL==0.AND.ITYPE(1)==50.AND.ITYPE(2)==50) THEN
          TAUWT = RU1*SUM(WO(IX2,JY1:JY2,KZ1:KZ2))/(REAL(NY-2)*REAL(NZ-2))*AP(IX2)*2.
          TAUWB = RU1*SUM(WO(IX1,JY1:JY2,KZ1:KZ2))/(REAL(NY-2)*REAL(NZ-2))*AP(IX1)*2.
          TAUW  = .5*(TAUWT+TAUWB)
          UBAR = SUM(WO((NX+1)/2:NX/2+1,JY1:JY2,KZ1:KZ2))/(REAL(NY-2)*REAL(NZ-2))
          UBAR = UBAR/REAL(NX/2+1-(NX+1)/2+1)
        ELSEIF(ICYL==1) THEN
          TAUW = RU1*SUM(WO(IX2,JY1:JY2,KZ1:KZ2))/(REAL(NY-2)*REAL(NZ-2))*AP(IX2)*2.
          UBAR =     SUM(WO(IX1,JY1:JY2,KZ1:KZ2))/(REAL(NY-2)*REAL(NZ-2))
        ENDIF
     
        CALL MPI_REDUCE(UBUK,UBUKG,1,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(TAUW,TAUWG,1,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(UBAR,UBARG,1,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)

        IF(MYRANK==0) THEN
          UBUKG = UBUKG/REAL(MYSIZE)
          TAUWG = TAUWG/REAL(MYSIZE)
          UBARG = UBARG/REAL(MYSIZE)

          OPEN(UNIT=100,FILE='timser.dat',POSITION='APPEND')
          WRITE(100,'(5G14.6)') TIME,TAUWG,SQRT(TAUWG),UBARG,UBUKG
          WRITE(*,'(5G14.6)') TIME,TAUWG,SQRT(TAUWG),UBARG,UBUKG
          CLOSE(100)
        ENDIF
C
      RETURN
      END

      SUBROUTINE VPFIELD3D(uo,vo,wo,p,tlevel,NAME,im,jm,km,kmg)

C
C-----------------------------------------------------------------------
C
C     PURPOSE:    Write instantaneous 3D flow field at given time
C
C-----------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
C
C-----------------------------------------------------------------------
C
c...Input variables and arrays.
      INTEGER im,jm,km,kmg
      REAL    UO(IM,JM,KM),VO(IM,JM,KM),WO(IM,JM,KM),P(IM,JM,KM)
      REAL    tlevel
      CHARACTER NAME*(*)
c
c...Local variables and arrays.
      INTEGER i,j,k
      INTEGER jp,DSIZE
      INTEGER STATUS(MPI_STATUS_SIZE)
      REAL    U1(IM,JM,KM),V1(IM,JM,KM),W1(IM,JM,KM),P1(IM,JM,KM)
c
c...write on disk
      IF (MYRANK==0) THEN
        open(166,file=NAME,form='unformatted')
        write(166) IM,JM,KMG,4
        write(6,*) '******************* 3D tlevel=',tlevel
        DO K=1,KZ2
          WRITE(166) ((REAL(UO(I,J,K),4),I=1,IM),J=1,JM)
     &          ,    ((REAL(VO(I,J,K),4),I=1,IM),J=1,JM)
     &          ,    ((REAL(WO(I,J,K),4),I=1,IM),J=1,JM)
     &          ,    ((REAL(P(I,J,K),4),I=1,IM),J=1,JM)
        ENDDO

        DO JP=1,MYSIZE-1
          DSIZE = IM*JM*KM
          CALL MPI_RECV(U1(1,1,1),DSIZE,MTYPE,JP,JP*4+0,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(V1(1,1,1),DSIZE,MTYPE,JP,JP*4+1,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(W1(1,1,1),DSIZE,MTYPE,JP,JP*4+2,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(P1(1,1,1),DSIZE,MTYPE,JP,JP*4+3,MPI_COMM_EDDY,STATUS,IERR)
          DO K=KZ1,KZ2+(JP+1)/MYSIZE
            WRITE(166) ((REAL(U1(I,J,K),4),I=1,IM),J=1,JM)
     &            ,    ((REAL(V1(I,J,K),4),I=1,IM),J=1,JM)
     &            ,    ((REAL(W1(I,J,K),4),I=1,IM),J=1,JM)
     &            ,    ((REAL(P1(I,J,K),4),I=1,IM),J=1,JM)

          ENDDO
        ENDDO

        WRITE(166) REAL(tlevel,4)
        CLOSE(166)

      ELSE
        DSIZE = IM*JM*KM
        CALL MPI_SEND(UO(1,1,1),DSIZE,MTYPE,0,MYRANK*4+0,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND(VO(1,1,1),DSIZE,MTYPE,0,MYRANK*4+1,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND(WO(1,1,1),DSIZE,MTYPE,0,MYRANK*4+2,MPI_COMM_EDDY,IERR)
        CALL MPI_SEND( P(1,1,1),DSIZE,MTYPE,0,MYRANK*4+3,MPI_COMM_EDDY,IERR)
      ENDIF

      RETURN

      END
c
c---- subroutine calc_utau----------------------------N. Beratis-26 Sep. 2011---
C
C     PURPOSE: Calc utau for pipe flow
C
C-------------------------------------------------------------------------------
      subroutine calc_utau(WO,XU,XC,NX,NY,NZ,UTAU,DIM,DPDZUTAU)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,dim
      real    utau(dim),dpdzutau
      real    xu(nx),xc(nx)
      real    wo(nx,ny,nz)
c
c.... Local arrays
      real    inv,tmp,tmp1

      inv = 1./real(ny-2)/real(nz-2)

      if(icyl==1) then
        tmp = SUM(WO(NX-1,2:NY-1,2:NZ-1))*inv
        CALL MPI_REDUCE(tmp,utau,1,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)      
        utau(1) = utau(1)/(xu(nx-1)-xc(nx-1))/real(mysize)
      else
        tmp = SUM(WO(2,2:NY-1,2:NZ-1))*inv
        CALL MPI_REDUCE(tmp,utau(1),1,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)      
        utau(1) = utau(1)/(xc(2)-xu(1))/real(mysize)

        tmp = SUM(WO(NX-1,2:NY-1,2:NZ-1))*inv
        CALL MPI_REDUCE(tmp,utau(2),1,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)      
        utau(2) = utau(2)/(xu(nx-1)-xc(nx-1))/real(mysize)
      endif
      utau = sqrt(abs(utau)*ru1)
      dpdzutau = dpdz

      return

      end
C-------------------------------------------------------------------------------


c
c---- subroutine calc_ubulk --------------------------N. Beratis-26 Sep. 2011---
C
C     PURPOSE: Calc ubulk for pipe flow
C
C-------------------------------------------------------------------------------
      subroutine calc_ubulk(WO,NX,NY,NZ,UB)
c
      include 'common.h'
c
c.... Input/Output arrays
      integer nx,ny,nz
      real    ub
      real    wo(nx,ny)
c
c.... Local arrays
      integer i,j
      real    area,areaij

      UB = 0.0
      AREA = 0.0
      DO J=JY1,JY2
      DO I=IX1,IX2
        AREAIJ = RP(I)/(AP(I)*BP(J))
        UB = UB + WO(I,J)*AREAIJ
        AREA = AREA+AREAIJ
      ENDDO
      ENDDO
      UB = UB/AREA
C
      return

      end
C-------------------------------------------------------------------------------


C---- read_utau_input --------------------------N. Beratlis-19 Sep. 2011---
C
C     PURPOSE: Read parameters from utau.input file
C
c--------------------------------------------------------------------------
      subroutine read_utau_input(flag,it)
c
      IMPLICIT NONE
c
c.... Input/Output arrays
      INTEGER flag,it
      
      open(unit=10,file='utau.input',form='formatted')
!       read(10,'(A)')
      read(10,'(A)') flag
!       read(10,'(A)')
      read(10,'(A)') it
      close(10)

      return

      end
c--------------------------------------------------------------------------


C---- write_utau------------------------------ N. Beratlis-27 Sep. 2011 ---
C
C     PURPOSE: Write utau to file
C
C--------------------------------------------------------------------------
      subroutine write_utau(utau,dim,dpdzutau,time,nt,flag)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer icycle,flag,nt,dim
      real    utau(nt,dim),dpdzutau(nt),time(nt)
c
c.... Local arrays
      integer i
      INTEGER, SAVE :: isave=0

      IF(MYRANK.EQ.0) THEN
        IF(flag==1 .AND. isave==0) THEN
          open(unit=10,file='utau.plt',form='formatted')
          if(icyl==0) then
            write(10,'(A)') 'VARIABLES = "T", "UTAU_BOT", '
     &        //'"UTAU_TOP", "DPDZ"'
          else 
             write(10,'(A)') 'VARIABLES = "T", "UTAU", "DPDZ"'
          endif
          isave = 1
        ELSE
          open(unit=10,file='utau.plt',form='formatted'
     &          ,position='append')
        ENDIF
        do i=1,1
          write(10,'(4(E16.8))') time(i),utau(i,:),dpdzutau(i)
        enddo
        close(10)
      ENDIF

      return

      end
c--------------------------------------------------------------------------
