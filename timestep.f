C-----------------------------------------------------------------------
      SUBROUTINE CALCFL(UO,VO,WO,TV,CFL,IM,JM,KM,CFLM,ICYCLE)
C-----------------------------------------------------------------------
C
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
C 
C-----------------------------------------------------------------------
C
      INTEGER IM,JM,KM,ICYCLE
      REAL CFL(IM,JM,KM)
      REAL UO(IM,JM,KM),VO(IM,JM,KM),WO(IM,JM,KM),TV(IM,JM,KM)
      REAL CFLM
C
      REAL     CFLL,CFLX,CFLY,CFLZ,CFLVX,CFLVY,CFLVZ,CFLL2,CFLM2
      REAL     TOTVX
      REAL     RFLAGX,RFLAGY,RFLAGCX,RFLAGCY
      REAL     CONST
      REAL     CFL2(IM,JM,KM)!,CFL3(IM,JM,KM)
      REAL     CFL1XL,CFL1YL,CFL1ZL,CFL2XL,CFL2YL,CFL2ZL
      REAL     CFL1XG,CFL1YG,CFL1ZG,CFL2XG,CFL2YG,CFL2ZG
      REAL     MAXL(19,MYSIZE),MAXG(19,MYSIZE)
      INTEGER  INDL(21,MYSIZE),INDG(21,MYSIZE)
c      REAL,    DIMENSION(:,:), ALLOCATABLE :: MAXL,MAXG
c      INTEGER, DIMENSION(:,:), ALLOCATABLE :: INDL,INDG
      INTEGER  IND(3,7),INDC(3,7),MAXRANK(19)
      INTEGER  I,J,K,STATUS(MPI_STATUS_SIZE)
      REAL     MAXC(19),MAXV(19),VALIN(19)
      REAL     TEMP
C 
C-----------------------------------------------------------------------
C
      IF(ICYL.EQ.0) THEN
        CONST = 2.
      ELSE
        CONST = 4.
      ENDIF
!      IF((IMPLX+IMPLY)==0) CONST = 2.
      RFLAGX = REAL(IMPLX)
      RFLAGY = REAL(IMPLY)
      RFLAGCX = REAL(IMPLCX)
      RFLAGCY = REAL(IMPLCY)
C 
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
c
c     Convective constraint 
c
        CFLX = (1.-RFLAGCX*RVIMPLX(I))*0.5*ABS(UO(I-1,J,K)+UO(I,J,K))*AP(I)
        CFLY = (1.-RFLAGCY*RVIMPLY(I))*0.5*ABS(VO(I,J-1,K)+VO(I,J,K))*BP(J)/RP(I)
        CFLZ = 0.5*ABS(WO(I,J,K-1)+WO(I,J,K))*CP(K)
c
c     Viscous constraint 
c
        TOTVX = RU1+TV(I,J,K)
C
        CFLVX = (1.-RFLAGX*RVIMPLX(I))* AP(I)**2
        CFLVY = (1.-RFLAGY*RVIMPLY(I))*(BP(J)/RP(I))**2
        CFLVZ =              CP(K)**2
C
c        CFL(I,J,K)  = CFLX+CFLY+CFLZ+CONST*TOTVX*(CFLVX+CFLVY+CFLVZ)
        CFL(I,J,K) = MAX(CFLX,CFLY,CFLZ,CONST*TOTVX*CFLVX
     &       ,CONST*TOTVX*CFLVY,CONST*TOTVX*CFLVZ)
c        CFL3(I,J,K) = CFLX+CFLZ+CONST*TOTVX*(CFLVX+CFLVZ)
      ENDDO
      ENDDO
      ENDDO
C
      I=7
      IND(:, I)=MAXLOC(CFL(IX1:IX2,JY1:JY2,KZ1:KZ2))
      IND(1,I) = IND(1,I)+IX1-1
      IND(2,I) = IND(2,I)+JY1-1
      IND(3,I) = IND(3,I)+KZ1-1
      VALIN(I) = CFL( IND(1,I),IND(2,I),IND(3,I) )
C
      CFLL = MAXVAL(CFL(IX1:IX2,JY1:JY2,KZ1:KZ2))
      CALL MPI_ALLREDUCE(CFLL,CFLM,1,MTYPE,MPI_MAX,MPI_COMM_EDDY,IERR)
      CFL2 = 0.0
C
C...Individual components
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
c     Convective x constraint
        CFL2(I,J,K) = 0.5*ABS(UO(I-1,J,K)+UO(I,J,K))*AP(I)*(1.-RFLAGCX*RVIMPLX(I))
      ENDDO
      ENDDO
      ENDDO
C
      I=1
      IND(:, I)=MAXLOC(CFL2(IX1:IX2,JY1:JY2,KZ1:KZ2))
      IND(1,I) = IND(1,I)+IX1-1
      IND(2,I) = IND(2,I)+JY1-1
      IND(3,I) = IND(3,I)+KZ1-1
      VALIN(I)   = CFL2( IND(1,I),IND(2,I),IND(3,I) )
      VALIN(I+7) = CFL2( IND(1,7),IND(2,7),IND(3,7) )
      VALIN(I+13)= CFL ( IND(1,I),IND(2,I),IND(3,I) )
c
      CFL2 = 0.0
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
c     Convective y constraint 
        CFL2(I,J,K) = 0.5*ABS(VO(I,J-1,K)+VO(I,J,K))*BP(J)/RP(I)*(1.-RFLAGCY*RVIMPLY(I))
      ENDDO
      ENDDO
      ENDDO
C
      I=2
      IND(:, I)=MAXLOC(CFL2(IX1:IX2,JY1:JY2,KZ1:KZ2))
      IND(1,I) = IND(1,I)+IX1-1
      IND(2,I) = IND(2,I)+JY1-1
      IND(3,I) = IND(3,I)+KZ1-1
      VALIN(I)   = CFL2( IND(1,I),IND(2,I),IND(3,I) )
      VALIN(I+7) = CFL2( IND(1,7),IND(2,7),IND(3,7) )
      VALIN(I+13)= CFL ( IND(1,I),IND(2,I),IND(3,I) )

c
      CFL2 = 0.0
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
c     Convective z constraint 
        CFL2(I,J,K) = 0.5*ABS(WO(I,J,K-1)+WO(I,J,K))*CP(K)
      ENDDO
      ENDDO
      ENDDO
C
      I=3
      IND(:, I)=MAXLOC(CFL2(IX1:IX2,JY1:JY2,KZ1:KZ2))
      IND(1,I) = IND(1,I)+IX1-1
      IND(2,I) = IND(2,I)+JY1-1
      IND(3,I) = IND(3,I)+KZ1-1
      VALIN(I)   = CFL2( IND(1,I),IND(2,I),IND(3,I) )
      VALIN(I+7) = CFL2( IND(1,7),IND(2,7),IND(3,7) )
      VALIN(I+13)= CFL ( IND(1,I),IND(2,I),IND(3,I) )
c
      CFL2 = 0.0
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
c     Viscous x constraint 
        TOTVX = RU1+TV(I,J,K)
        CFL2(I,J,K) = (1.-RFLAGX*RVIMPLX(I))*CONST*TOTVX*(AP(I)**2)
      ENDDO
      ENDDO
      ENDDO
C
      I=4
      IND(:, I)=MAXLOC(CFL2(IX1:IX2,JY1:JY2,KZ1:KZ2))
      IND(1,I) = IND(1,I)+IX1-1
      IND(2,I) = IND(2,I)+JY1-1
      IND(3,I) = IND(3,I)+KZ1-1
      VALIN(I)   = CFL2( IND(1,I),IND(2,I),IND(3,I) )
      VALIN(I+7) = CFL2( IND(1,7),IND(2,7),IND(3,7) )
      VALIN(I+13)= CFL ( IND(1,I),IND(2,I),IND(3,I) )
c
      CFL2 = 0.0
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
c     Viscous y constraint 
        TOTVX = RU1+TV(I,J,K)
        CFL2(I,J,K) = (1.-RFLAGY*RVIMPLY(I))*CONST*TOTVX*(BP(J)/RP(I))**2
      ENDDO
      ENDDO
      ENDDO

      I=5
      IND(:, I)=MAXLOC(CFL2(IX1:IX2,JY1:JY2,KZ1:KZ2))
      IND(1,I) = IND(1,I)+IX1-1
      IND(2,I) = IND(2,I)+JY1-1
      IND(3,I) = IND(3,I)+KZ1-1
      VALIN(I)   = CFL2( IND(1,I),IND(2,I),IND(3,I) )
      VALIN(I+7) = CFL2( IND(1,7),IND(2,7),IND(3,7) )
      VALIN(I+13)= CFL ( IND(1,I),IND(2,I),IND(3,I) )

c
      DO K=KZ1,KZ2
      DO J=JY1,JY2
      DO I=IX1,IX2
c     Viscous z constraint 
        TOTVX = RU1+TV(I,J,K)
        CFL2(I,J,K) = CONST*TOTVX*CP(K)**2
      ENDDO
      ENDDO
      ENDDO
C
      I=6
      IND(:, I)=MAXLOC(CFL2(IX1:IX2,JY1:JY2,KZ1:KZ2))
      IND(1,I) = IND(1,I)+IX1-1
      IND(2,I) = IND(2,I)+JY1-1
      IND(3,I) = IND(3,I)+KZ1-1
      VALIN(I)   = CFL2( IND(1,I),IND(2,I),IND(3,I) )
      VALIN(I+7) = CFL2( IND(1,7),IND(2,7),IND(3,7) )
      VALIN(I+13)= CFL ( IND(1,I),IND(2,I),IND(3,I) )


c      ALLOCATE(MAXL(19,MYSIZE),INDL(7*3,MYSIZE))
c      ALLOCATE(MAXG(19,MYSIZE),INDG(7*3,MYSIZE))
c...Convert index for z into global index.
      MAXL = 0.
      INDL = 0.
      DO I=1,6
        MAXL(I,MYRANK+1)=VALIN(I)
        K = (I-1)*3
        INDL(1+K,MYRANK+1)=IND(1,I)
        INDL(2+K,MYRANK+1)=IND(2,I)
        INDL(3+K,MYRANK+1)=IND(3,I)+MYRANK*(KM-2)
      ENDDO
      DO I=7,19 !13
        MAXL(I,MYRANK+1)=VALIN(I)
      ENDDO

      INDL(19,MYRANK+1)=IND(1,7)
      INDL(20,MYRANK+1)=IND(2,7)
      INDL(21,MYRANK+1)=IND(3,7)+MYRANK*(KM-2)
 
      CALL MPI_REDUCE(MAXL,MAXG,MYSIZE*19,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
      CALL MPI_REDUCE(INDL,INDG,MYSIZE*3*7,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)

      IF(MYRANK.EQ.0) THEN
        MAXRANK = MAXLOC(MAXG,DIM=2)
        DO I=1,6
          J = MAXRANK(I)
          MAXV(I) = MAXG(I,J)
          K = (I-1)*3
          INDC(1,I) = INDG(1+K,J)
          INDC(2,I) = INDG(2+K,J)
          INDC(3,I) = INDG(3+K,J)
        ENDDO
        DO I=7,13
          J = MAXRANK(7)
          MAXV(I) = MAXG(I,J)
          INDC(1,7) = INDG(19,J)
          INDC(2,7) = INDG(20,J)
          INDC(3,7) = INDG(21,J)
        ENDDO

        DO I=14,19
          J = MAXRANK(7)
          MAXV(I) = MAXG(I,J)
        ENDDO

        IF(ITSCR>0 .AND. MOD(ICYCLE,ITSCR)==0) THEN
        WRITE(6,*) 'MAX of CFL Components :'

        IF(implcx==0 .OR. (implcx==1 .AND. implxdcmp==1)) THEN
          write(6,102) 'CFLX =',
     &          MAXV(1),', @(',INDC(1,1),INDC(2,1),INDC(3,1),')'
     &       ,', CFL=',MAXV(14)
        ENDIF

        IF(implcy==0 .OR. (implcy==1 .AND. implydcmp==1)) THEN
          write(6,102) 'CFLY =',
     &            MAXV(2),', @(',INDC(1,2),INDC(2,2),INDC(3,2),')'
     &         ,', CFL=',MAXV(15)
        ENDIF

        write(6,102) 'CFLZ =',
     &          MAXV(3),', @(',INDC(1,3),INDC(2,3),INDC(3,3),')'
     &       ,', CFL=',MAXV(16)

        IF(implx==0 .OR. (implx==1 .AND. implxdcmp==1)) THEN
          write(6,102) 'CFLVX=',
     &          MAXV(4),', @(',INDC(1,4),INDC(2,4),INDC(3,4),')'
     &       ,', CFL=',MAXV(17)
        ENDIF

        IF(imply==0 .OR. (imply==1 .AND. implydcmp==1)) THEN
          write(6,102) 'CFLVY=',
     &            MAXV(5),', @(',INDC(1,5),INDC(2,5),INDC(3,5),')'
     &         ,', CFL=',MAXV(18)
        ENDIF

        write(6,102) 'CFLVZ=',
     &          MAXV(6),', @(',INDC(1,6),INDC(2,6),INDC(3,6),')'
     &       ,', CFL=',MAXV(19)

        write(6,*) 'MAX CFL# :'
        write(6,101) 'CFL  =',
     &             CFLM,', @(',INDC(1,7),INDC(2,7),INDC(3,7),')'

        IF(implcx==0 .OR. (implcx==1 .AND. implxdcmp==1)) THEN
          write(6,101) 'CFLX =',
     &          MAXV(8),', @(',INDC(1,7),INDC(2,7),INDC(3,7),')'
        ENDIF

        IF(implcy==0 .OR. (implcy==1 .AND. implydcmp==1)) THEN
          write(6,101) 'CFLY =',
     &            MAXV(9),', @(',INDC(1,7),INDC(2,7),INDC(3,7),')'
        ENDIF

        write(6,101) 'CFLZ =',
     &          MAXV(10),', @(',INDC(1,7),INDC(2,7),INDC(3,7),')'

        IF(implx==0 .OR. (implx==1 .AND. implxdcmp==1)) THEN
          write(6,101) 'CFLVX=',
     &          MAXV(11),', @(',INDC(1,7),INDC(2,7),INDC(3,7),')'
        ENDIF

        IF(imply==0 .OR. (imply==1 .AND. implydcmp==1)) THEN
          write(6,101) 'CFLVY=',
     &            MAXV(12),', @(',INDC(1,7),INDC(2,7),INDC(3,7),')'
        ENDIF

        write(6,101) 'CFLVZ=',
     &          MAXV(13),', @(',INDC(1,7),INDC(2,7),INDC(3,7),')'
      ENDIF
      ENDIF

c      DEALLOCATE(MAXL,MAXG,INDL,INDG)

101   FORMAT(A,E15.8,A,3(1x,I6),A) 
102   FORMAT(A,E15.8,A,3(1x,I6),A,A,E15.8) 
c
      END 
c
c
C---- subroutine read_cflave------------------------------------------------
C
C     PURPOSE: Read average CFL from file
C
C---------------------------------------------------------------------------
      subroutine read_cflave(filename,cfl,nx,nsmpl)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output array
      integer nx,nsmpl
      real    cfl(nx,6)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,iv,fh,filemode,newtype,nzg
      real    temp
      INTEGER(KIND=8) disp,offset
      integer status(mpi_status_size)

      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = MPI_MODE_RDONLY

      call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mtype,newtype
     &     ,'native',mpi_info_null,ierr)

      offset = 0
      call mpi_file_read_at(fh,offset,temp,1,newtype,status,ierr)
      nsmpl = int(temp)

      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      call mpi_type_contiguous(nx,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = MPI_MODE_RDONLY

      call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     &       ,'native',mpi_info_null,ierr)

      do iv=1,6
        offset = 1+nx*(iv-1)
        call mpi_file_read_at(fh,offset,cfl(:,iv),1,newtype,status,ierr)
      enddo

      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      return

      end
C---------------------------------------------------------------------------


C---- compute_cflave-----------------------------N. Beratlis-Sep. 18 2011---
C
C     PURPOSE: Compute cfl ave
C
C---------------------------------------------------------------------------
      subroutine compute_cflave(cfl,uo,vo,wo,tv,dp,nx,ny,nz,nsmpl)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,nsmpl
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),tv(nx,ny,nz),dp(nx,ny,nz)
      real    cfl(nx,6)
c
c.... Local arrays
      integer i,j,k
      real    totvx,const,s,invs
      real    cflcx,cflcy,cflcz,cflvx,cflvy,cflvz
      real    cfll(6),cflm(6)
c      real    cflg(nx,6)

      s = real(nsmpl)
      invs = 1./(s+1.0)

      const = 4.0
      IF((IMPLX+IMPLY)==0) CONST = 2.

      do i=ix1,ix2
      do k=kz1,kz2      
      do j=jy1,jy2
c
c... Convective terms
        dp(1,j,k) = 0.5*ABS(UO(I-1,J,K)+UO(I,J,K))*AP(I)
        dp(2,j,k) = 0.5*ABS(VO(I,J-1,K)+VO(I,J,K))*BP(J)/RP(I)
        dp(3,j,k) = 0.5*ABS(WO(I,J,K-1)+WO(I,J,K))*CP(K)
c
c... Viscous terms
        TOTVX = RU1+TV(I,J,K)

        dp(4,j,k) = CONST*TOTVX*AP(I)**2
        dp(5,j,k) = CONST*TOTVX*(BP(J)/RP(I))**2
        dp(6,j,k) = CONST*TOTVX*CP(K)**2

      enddo
      enddo

      do j=1,6
        CFLL(j) = MAXVAL(DP(j,JY1:JY2,KZ1:KZ2))
      enddo
      CALL MPI_ALLREDUCE(CFLL,CFLM,6,MTYPE,MPI_MAX,MPI_COMM_EDDY,IERR)

      do j=1,6
        cfl(i,j) = (cfl(i,j)*s + cflm(j))*invs
      enddo

      enddo

      nsmpl = nsmpl+1

      return

      end
C---------------------------------------------------------------------------


C---- subroutine write_cflave -------------------N. Beratlis-18 Sep. 2011---
C
C     PURPOSE: Write average CFL to file
C
C---------------------------------------------------------------------------
      subroutine write_cflave(filename,cfl,nx,nsmpl)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,nsmpl
      real    cfl(nx,6)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,iv,fh,filemode,newtype,nzg,k1
      real    cflg(nx,6)
      INTEGER(KIND=8) disp,offset
      integer status(mpi_status_size)
c      
      call mpi_type_contiguous(1,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)

      call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)

      disp = 0
      call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     &     ,'native',mpi_info_null,ierr)

      offset = 0
      if(myrank==0) call mpi_file_write_at(fh,offset,real(nsmpl),1,newtype,status,ierr)
      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      call mpi_type_contiguous(nx,mtype,newtype,ierr)

      call mpi_type_commit(newtype,ierr)
      filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      call mpi_file_open(mpi_comm_self,trim(filename),filemode,mpi_info_null,fh,ierr)
      disp = 0
      call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     &       ,'native',mpi_info_null,ierr)

      if(mysize>0) then
        CALL MPI_ALLREDUCE(cfl,cflg,nx*6,mtype,mpi_sum,mpi_comm_eddy,ierr)
        cfl = cflg/real(mysize)
      endif

      if(myrank==0) then
        do iv=1,6
          offset = 1+nx*(iv-1)
          call mpi_file_write_at(fh,offset,cfl(:,iv),1,newtype,status,ierr)
        enddo

        open(unit=10,file='cflave.plt',form='formatted')
        write(10,'(A)') 'FILETYPE = "Solution"'
        write(10,'(A)') 'VARIABLES = "CFLCX", "CFLCY", "CFLCZ",'
     &       //' "CFLVX", "CFLVY", "CFLVZ"'
        do i=1,nx
          write(10,'(6(E14.7))') cfl(i,:)
        enddo
        close(10)
      endif

      call mpi_type_free(newtype,ierr)
      call mpi_file_close(fh,ierr)

      return

      end
C---------------------------------------------------------------------------
