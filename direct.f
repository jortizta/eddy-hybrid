C-----------------------------------------------------------------------
      SUBROUTINE DIRECT(DP,DUMMY,DUMMYT,DELYSQ,DELZSQ,AP,AU,CP,CW,RP,RU,
     &     NX,NY,NZ,NZG,MYSIZE,MYRANK,MPI_COMM_EDDY,MTYPE,IPRES,MP,LP,NP)
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
      include 'mpif.h'
C
C     ...ARRAY ARGUMENTS...
C
      INTEGER NX,NY,NZ,NZG
      INTEGER IPRES
      INTEGER MYRANK, MYSIZE, MPI_COMM_EDDY, MTYPE
      REAL AP(*), AU(*), CP(*), CW(*), RP(*), RU(*)
      REAL DP(NX,NY,NZ)
      REAL DUMMY(NX*NY*NZ),DUMMYT(NX*NY*NZ)
      REAL DELYSQ,DELZSQ
C
C     ...LOCAL SCALARS...
C
      INTEGER I, J, JL, K, ICR, IM, IN, NJSKIP, NJSET, JJ, PLSET
      INTEGER NWK, IERROR
      INTEGER MP, LP, NP
      INTEGER,SAVE :: L, M, N, NL, LL
      INTEGER, SAVE :: ICYC=0
C
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: SNDCOUNTS,SDISPL
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: RCVCOUNTS,RDISPL
C
      LOGICAL,SAVE :: INITP
      DATA INITP /.TRUE./
C
      REAL PI, DPM
C
C     LOCAL ARRAYS
C
      REAL, DIMENSION(:), ALLOCATABLE :: BMM,CMM,DPX,DPY,DPZ
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: RHS1D
      REAL, DIMENSION(:,:), ALLOCATABLE :: RHSP
C
c      REAL RHSP(NZG-2,NX-2)

      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: AM, BM, CM
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: AN, BN, CN
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: W
      REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: W2, BML
C
C     LOCAL ARRAYS FOR FFT
C
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: AK,AL
      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: WSAVE,WSAVF
c
c     LOCAL ARRAYS FOR SWAPPING
      INTEGER II,IP,ISWAP,IERR,STATUS(MPI_STATUS_SIZE)

      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: SENDRECVPROC
      REAL,   DIMENSION(:,:,:),ALLOCATABLE :: DUMMY1,DUMMY2
      REAL,   DIMENSION(:,:,:),ALLOCATABLE :: WRK1,WRK2
c
c.... PDC2D
      INTEGER, SAVE :: ldw, liw, ilf, iuf, nprocj, iproc, world_group, jpl, comm
      integer nlmx, nr, p, pe, k1proc, k2proc, k1, k2, k1g, k2g, kg,
     &     iset, nplanes, pow, itmp
      real    tmp,crit,critmax
      INTEGER, SAVE :: nsets,nsetsmax
      real, save ::    ch
      INTEGER, DIMENSION(:), ALLOCATABLE :: RANKSJ, GROUPJ
      INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: IW
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: DW
      integer, dimension(:), ALLOCATABLE, SAVE :: ILU,IMU,JORDER,JPLANE
      integer, save :: COMMJ1
      logical init(3), first
      parameter (nlmx = 6, nr = 2)
      integer, save :: pemn,pemx
      character index2(0:99)*2 
      data first /.true./      
c
c.... Timing arrays
      REAL clocktemp,clock(3)
      REAL tclock
      
C-----------------------------------------------------------------------
C     INITIALIZATION
C-----------------------------------------------------------------------
      dummy = 0.0
      dummyt = 0.0      

      k = 0
      do i=0,9
      do j=0,9        
        index2(k) = char(i+48)//char(j+48)
        k=k+1
      enddo
      enddo

      iswap = 0
      icr = 1
      if(ipres==3) icr=0
c      if(ipres==4) icr=0
      icyc = icyc+1

      clock = 0.0

      N = NZG-2
      NL= NZ-2
      M = NX-2
      L = NY-2
      IF(IPRES/=6) THEN
        LL= L/MYSIZE
        IF(IPRES==5) LL=1
      ENDIF

      IF(IPRES==2 .OR. IPRES==3) THEN
        if(icr==0) then
          ALLOCATE(RHSP(M,N))
        else
          ALLOCATE(RHSP(N,M))
        endif
      ENDIF

      IF(INITP) THEN

        INITP = .FALSE.        
C-----------------------------------------------------------------------
C     DIMENSIONS
C-----------------------------------------------------------------------
C
        IF(ISWAP==1) THEN
          ALLOCATE(SENDRECVPROC(0:MYSIZE-1))
          CALL CLOSESTRANK(SENDRECVPROC,MYRANK,MYSIZE)
        ENDIF

        ALLOCATE(SNDCOUNTS(0:MYSIZE-1),SDISPL(0:MYSIZE-1))
        ALLOCATE(RCVCOUNTS(0:MYSIZE-1),RDISPL(0:MYSIZE-1))
C
        K = M*LL*NL
        SNDCOUNTS = K
        RCVCOUNTS = K
        DO I=0,MYSIZE-1
          SDISPL(I) = I*K
          RDISPL(I) = I*K
        ENDDO

C-----------------------------------------------------------------------
C     DIMENSION OF ARRAY FOR WORK SPACE AND SAVED ARRAYS
C-----------------------------------------------------------------------
        IF(IPRES==1) THEN

          if(icr==0) then
            IF(NP==0) THEN
              K = INT(LOG10(REAL(N-1))/LOG10(REAL(2)))+1
              I = 2**(K+1)
              NWK = (K-2)*I+K+5+2*N+MAX(2*N,6*M)
            ELSEIF(NP==1) THEN
              K = INT(LOG10(REAL(N))/LOG10(REAL(2)))+1
              I = 2**(K+1)
              NWK = (K-2)*I+K+5+MAX(2*N,6*M)
            ENDIF
            ALLOCATE(W(NWK))
          else
            IF(NP==0) THEN
              K = INT(LOG10(REAL(M-1))/LOG10(REAL(2)))+1
              I = 2**(K+1)
              NWK = (K-2)*I+K+5+2*M+MAX(2*M,6*N)
            ELSEIF(NP==1) THEN
              K = INT(LOG10(REAL(M))/LOG10(REAL(2)))+1
              I = 2**(K+1)
              NWK = (K-2)*I+K+5+MAX(2*M,6*N)
            ENDIF
            ALLOCATE(W2(LL,NWK))
c            ALLOCATE(W(NWK))
          endif

        ELSEIF(IPRES==2) THEN

          NWK=4*N+(10+INT(LOG10(REAL(N))/LOG10(2.)))*M
          ALLOCATE(W(NWK))

        ELSEIF(IPRES==4 .OR. IPRES==5) THEN

          nsetsmax = 1
          nsets = 1
          pemn = 0
          pemx = 0

          IF(ICR==0) then
            I = 1+max(int(log10(real(M))/log10(4.0)),0)
            LDW = 6*I*M + max(9*M,11*N)
            LIW = 6*M + 2*I + (4**I - 1)/3 + 7 + int(log10(real(2*mysize))/log10(4.0))

            !Check for parallel solution
            IF(MYSIZE>NY-2) THEN
              tmp = REAL(MYSIZE)/REAL(NY-2)
              tmp = LOG10(tmp)/LOG10(2.0)

              IF(tmp==int(tmp)) THEN
                pemx = int(tmp)

                tmp = 2**pemx
                nprocj = tmp
                jpl = floor(real(myrank)/real(nprocj))+1

                ALLOCATE(GROUPJ(L),RANKSJ(0:NPROCJ-1))

                CALL MPI_COMM_GROUP(MPI_COMM_WORLD,WORLD_GROUP,ierr)
                do i=1,nprocj
                  ranksj(i-1) = (jpl-1)*nprocj+(i-1)
                enddo
                CALL MPI_GROUP_INCL(world_group,nprocj,ranksj,groupj,ierr)
                CALL MPI_COMM_CREATE(MPI_COMM_WORLD,groupj,comm,ierr)
                if(jpl==1) commj1 = comm
                DEALLOCATE(RANKSJ,GROUPJ)
              ELSE

                IPRES=6

                critmax = 100000.0
                DO I=1,3
                  tmp = REAL(MYSIZE)/real(2**I)
                  IF(tmp<L) THEN
                    nsetsmax = ceiling(real(L)/floor(tmp))
                    crit = real(nsetsmax)/real(2**I)
                    IF(crit<critmax) THEN
                      p = nsetsmax
                      critmax = crit
                      pow=i
                    ENDIF
                  ENDIF
                ENDDO
                nsetsmax = p
                
                nprocj = 2**pow
                njskip = CEILING(REAL(MYSIZE)/REAL(nprocj))
                njset = CEILING(REAL(MYSIZE)/REAL(nprocj))
                plset = (REAL(L)/REAL(NJSKIP)-FLOOR(REAL(L)/REAL(NJSKIP)))*MYSIZE
                nsets = nsetsmax
                if(myrank>plset) nsets = nsetsmax-1
                
                ALLOCATE(jplane(nsets))                
                do i=1,nsets
                  jplane(i) = floor(real(myrank)/real(nprocj))+1+njset*(i-1)
                  if(jplane(i)>L) then
                    jplane(i) = 0
                  endif
                enddo
                ll = nsets

                ALLOCATE(GROUPJ(L),RANKSJ(0:NPROCJ-1))
                CALL MPI_COMM_GROUP(MPI_COMM_WORLD,WORLD_GROUP,ierr)

                do i=1,nprocj
                  ranksj(i-1) = (jplane(1)-1)*nprocj+(i-1)
                enddo
                CALL MPI_GROUP_INCL(world_group,nprocj,ranksj,groupj,ierr)
                CALL MPI_COMM_CREATE(MPI_COMM_WORLD,groupj,comm,ierr)
                if(jplane(1)==1) commj1 = comm
                DEALLOCATE(RANKSJ,GROUPJ)
                K=0
                ALLOCATE(JORDER(L))
                DO I=1,NJSKIP
                DO J=1,NSETSMAX !CEILING(REAL(L)/REAL(NJSKIP))
                  P = I+(J-1)*NJSKIP
                  IF(P<=L) THEN
                    K=K+1
                    JORDER(K) = P
                  ENDIF
                ENDDO
                ENDDO
              ENDIF
            ELSE
              call MPI_COMM_SPLIT(MPI_COMM_WORLD,myrank,myrank,comm,ierr)
            ENDIF
            ALLOCATE(DW(LL,LDW),IW(LL,LIW))
            iw = 0
            dw = 0.0
            CH = 0.0
          ELSE
            I = 1+max(int(log10(real(N))/log10(4.0)),0)
            LDW = 6*I*N + max(9*N,11*M)
            LIW = 6*N + 2*I + (4**I - 1)/3 + 7 + int(log10(real(2*mysize))/log10(4.0))
            ALLOCATE(DW(1,LDW),IW(1,LIW))
            CH = 0.0

            !Check for parallel solution
            IF(MYSIZE>NY-2 .AND. MOD(MYSIZE,NY-2)==0) THEN
              tmp = (MYSIZE)/(NY-2)
              pemx = int(LOG10(tmp)/LOG10(2.0))

              tmp = 2**pemx
              nprocj = tmp
              jpl = floor(real(myrank)/real(nprocj))+1

              ALLOCATE(GROUPJ(L),RANKSJ(0:NPROCJ-1))

              CALL MPI_COMM_GROUP(MPI_COMM_WORLD,WORLD_GROUP,ierr)

              do i=1,nprocj
                ranksj(i-1) = (jpl-1)*nprocj+(i-1)
              enddo
              CALL MPI_GROUP_INCL(world_group,nprocj,ranksj,groupj,ierr)
              CALL MPI_COMM_CREATE(MPI_COMM_WORLD,groupj,comm,ierr)
              if(jpl==1) commj1 = comm

              DEALLOCATE(RANKSJ,GROUPJ)
            ELSE
              call MPI_COMM_SPLIT(MPI_COMM_WORLD,myrank,myrank,comm,ierr)
            ENDIF
            CH = 0.0
          ENDIF

        ENDIF
C
        ALLOCATE(WSAVE(2*L+15))
        ALLOCATE(AK(L))
        ALLOCATE(AM(M),BM(M),CM(M))
        ALLOCATE(BML(M,LL))
C     
C-----------------------------------------------------------------------
C     MODIFIED WAVE NUMBERS
C-----------------------------------------------------------------------
        PI = 4.0*ATAN(1.0)
       
        DO K=1, L/2
          AK(K) =  2.*PI*(K-1)
        ENDDO
C
        DO K=L/2+1, L
          AK(K) = -2.*PI*(L-K+1)
        ENDDO
C
        AK(1:L) = 2.*(1.-COS(AK(1:L)/REAL(L)))/DELYSQ
C
C-----------------------------------------------------------------------
C     CREATE MATRIX FOR BLKTRI 
C-----------------------------------------------------------------------
        IF(IPRES==4 .OR. IPRES==5. OR. IPRES==6) THEN
          AM(1:M) = -RU(1:M)*AU(1:M)
          BM(1:M) = (AU(1:M  )*RU(1:M  ) + AU(2:M+1)*RU(2:M+1))
          CM(1:M) = RP(2:M+1)/AP(2:M+1)
          !Set boundary coefficients later

        ELSE

          AM(1:M) = AP(2:M+1)/RP(2:M+1)*AU(1:M  )*RU(1:M  )
          CM(1:M) = AP(2:M+1)/RP(2:M+1)*AU(2:M+1)*RU(2:M+1)

          IF(MP==1) THEN
            AM(1)=0.
            CM(M)=0.
          ENDIF

          BM = - AM - CM
        ENDIF

        IF(IPRES==1 .OR. IPRES==4 .OR. IPRES==5 .OR. IPRES==6) THEN

          ALLOCATE(AN(N),BN(N),CN(N))
C
          IF(IPRES==4 .OR. IPRES==5 .OR. IPRES==6) THEN
            AN(1:N) = -CW(1:N)
            BN(1:N) = CW(1:N) + CW(2:N+1)
            CN(1:N) = 1.0/CP(2:N+1)

            IF(NP==1) THEN
              BN(1) = BN(1) + AN(1)
              BN(N) = BN(N) - CW(N+1)
            ENDIF
            AN(1) = 0.0
          ELSE
            AN(1:N) = CP(2:N+1)*CW(1:N  )
            CN(1:N) = CP(2:N+1)*CW(2:N+1)

            IF(NP==1) THEN
              AN(1)=0.
              CN(N)=0.
            ENDIF
         
            BN = - AN - CN
          ENDIF
C
C.... Initialize solver
C
          IF(IPRES==4 .OR. IPRES==5 .OR. IPRES==6) THEN

            init(1) = .not.first
            first = .false.
            init(2) = .true.
            init(3) = .true.

            DO JL = 1, LL       ! LOCAL INDEX
              IF(IPRES==5) THEN
                J = JPL 
              ELSEIF(IPRES==6) THEN
                J = JPLANE(JL)
              ELSE
                J = LL*MYRANK+JL ! GLOBAL INDEX
              ENDIF

              BML(1:M,JL) = BM(1:M) + CM(1:M)*AK(J/2+1)/RP(2:M+1)**2

              !Change coefficients on boundaries
              IF(MP==1) THEN
                BML(1,JL) = BML(1,JL) + AM(1)
                BML(M,JL) = BML(M,JL) - RU(M+1)*AU(M+1)
              ENDIF
            ENDDO
            AM(1) = 0.0
            
            IF(ICR==0) THEN
              ALLOCATE(RHSP(M,N))
              DO JL = 1, LL
                IF(IPRES==5) THEN
                  J = JPL 
                ELSEIF(IPRES==6) THEN
                  J = JPLANE(JL)
                ELSE
                  J=LL*MYRANK+JL ! GLOBAL INDEX
                ENDIF
                call pdc2dn(m,n,rhsp,n,ilf,iuf,am,bml(:,jl),cm,an,bn,cn,ch,
     &               dw(jl,:),ldw,iw(jl,:),liw,comm,init,ierr)

              ENDDO
              DEALLOCATE(RHSP)
              ALLOCATE(rhs1d((iuf-ilf+1)*n))

            ELSE
              ALLOCATE(RHSP(N,M))
              call pdc2dn(n,m,rhsp,m,ilf,iuf,an,bn,cn,am,bm,cm,ch,
     &              dw(1,:),ldw,iw(1,:),liw,comm,init,ierr)
              DEALLOCATE(RHSP)
              ALLOCATE(rhs1d((iuf-ilf+1)*m))
            ENDIF
c
c.... Setup arrays for domain swap
            IF(IPRES==6) THEN

              ALLOCATE(ilu(nprocj),imu(nprocj))
            
              ilu = 0
              imu = 0
              IF(myrank>=0 .AND. myrank<nprocj) THEN
                ilu(myrank+1) = ilf
                imu(myrank+1) = iuf
              ENDIF

              DO i=1,nprocj
                CALL MPI_ALLREDUCE(ilu(i),j,1,mpi_integer,MPI_SUM,MPI_COMM_EDDY,IERR)
                ilu(i) = j
                CALL MPI_ALLREDUCE(imu(i),j,1,mpi_integer,MPI_SUM,MPI_COMM_EDDY,IERR)
                imu(i) = j
              ENDDO
 
              IF(ICR==0) THEN
                DO J=1,L
                DO I=1,NPROCJ
                  IPROC = (J-1)*NPROCJ+(I-1)
                  IF(IPROC<MYSIZE) THEN
                    RCVCOUNTS(IPROC) = (IUF-ILF+1)*NL*NSETS
                    RDISPL(IPROC) = IPROC*(IUF-ILF+1)*NL*NSETS
                    IF(IPROC>PLSET) THEN
                      SNDCOUNTS(IPROC) = (IMU(I)-ILU(I)+1)*NL*(NSETSMAX-1)
                    ELSE
                      SNDCOUNTS(IPROC) = (IMU(I)-ILU(I)+1)*NL*NSETSMAX
                    ENDIF
                    SDISPL(IPROC) = FLOOR(REAL(IPROC)/REAL(NPROCJ))*(IMU(I)-ILU(I)+1)*NL + (ILU(I)-1)*NL*L
                  ENDIF
                ENDDO
                ENDDO
              ENDIF

            ELSEIF(IPRES==5) THEN

              ALLOCATE(ilu(nprocj),imu(nprocj))
            
              ilu = 0
              imu = 0
              IF(myrank>=0 .AND. myrank<nprocj) THEN
                ilu(myrank+1) = ilf
                imu(myrank+1) = iuf
              ENDIF

              DO i=1,nprocj
                CALL MPI_ALLREDUCE(ilu(i),j,1,mpi_integer,MPI_SUM,MPI_COMM_EDDY,IERR)
                ilu(i) = j
                CALL MPI_ALLREDUCE(imu(i),j,1,mpi_integer,MPI_SUM,MPI_COMM_EDDY,IERR)
                imu(i) = j
              ENDDO
            
              K = M*LL
              IF(ICR==0) THEN
                DO J=1,L
                DO I=1,NPROCJ
                  IPROC = (J-1)*NPROCJ+(I-1)
                  RCVCOUNTS(IPROC) = (IUF-ILF+1)*NL
                  SNDCOUNTS(IPROC) = (IMU(I)-ILU(I)+1)*NL
                  RDISPL(IPROC) = IPROC*(IUF-ILF+1)*NL
                  SDISPL(IPROC) = (J-1)*M*NL + (ILU(I)-1)*NL
                ENDDO
                ENDDO
              ELSE
                SDISPL = 0
                RDISPL = 0
                SNDCOUNTS = 0
                RCVCOUNTS = 0
                K1G = MYRANK*NL+1
                K2G = K1G+NL-1
                DO J=1,L
                DO I=1,NPROCJ
                  IPROC = (J-1)*NPROCJ+(I-1)
                  K1PROC = IPROC*NL+1
                  K2PROC = K1PROC+NL-1
                  IF(K1G>=ILU(I) .OR. K2G<=IMU(I)) THEN
                    K2 = MIN(K2G,IMU(I))
                    K1 = MAX(K1G,ILU(I))
                    IF(K2>K1) THEN
                      SNDCOUNTS(IPROC) = M*(K2-K1+1)
                      SDISPL(IPROC) = (J-1)*M*NL + (K1-K1G)*M
                    ENDIF
                  ENDIF
                  IF(K1PROC>=ILF .OR. K2PROC<=IUF) THEN
                    K2 = MIN(K2PROC,IUF)
                    K1 = MAX(K1PROC,ILF)
                    IF(K2>K1) THEN
                      RCVCOUNTS(IPROC) = M*(K2-K1+1)
                      RDISPL(IPROC) = (K1-ILF)*M
                    ENDIF
                  ENDIF
                ENDDO
                ENDDO
              ENDIF
            ENDIF            
          ELSE
            ALLOCATE(rhs1d(m*n))
            DO JL = 1, LL      ! LOCAL INDEX
              J=LL*MYRANK+JL   ! GLOBAL INDEX
              BML(1:M,JL) = BM(1:M) - AK(J/2+1)/RP(2:M+1)**2
            ENDDO
            if(icr==0) then
              CALL BLKTRI(0,NP,N,AN,BN,CN,MP,M,AM,BM,CM,M,RHS1D,IERROR,W)
            else
              DO JL = 1, LL      ! LOCAL INDEX
                CALL BLKTRI(0,MP,M,AM,BML(:,JL),CM,NP,N,AN,BN,CN,N,RHS1D,IERROR,W2(JL,:))
              ENDDO
            endif
            IF(IERROR/=0) WRITE(6,*) 'INIT. BLKTRI: IERROR = ',IERROR
          ENDIF

        ELSEIF(IPRES==2) THEN

          AM(:)=AM(:)*DELZSQ
          BM(:)=BM(:)*DELZSQ
          CM(:)=CM(:)*DELZSQ
          AK(:)=AK(:)*DELZSQ

        ELSEIF(IPRES==3) THEN

          DO JL = 1, LL        ! LOCAL INDEX
            J=LL*MYRANK+JL     ! GLOBAL INDEX
            BML(1:M,JL) = BM(1:M) - AK(J/2+1)/RP(2:M+1)**2
          ENDDO

          ALLOCATE(WSAVF(2*N+15))
          ALLOCATE(AL(N))

C-----------------------------------------------------------------------
C     MODIFIED WAVE NUMBERS
C-----------------------------------------------------------------------
          DO K=1,N/2
            AL(K)= 2.*PI*(K-1)
          ENDDO
          DO K=N/2+1,N
            AL(K)=-2.*PI*(N-K+1)
          ENDDO
          
          AL(1:N)=2.*(1.-COS(AL(1:N)/REAL(N)))/DELZSQ

          CALL RFFTI(N,WSAVF)

        ENDIF
C-----------------------------------------------------------------------
C     INITIALIZE FFT
C-----------------------------------------------------------------------
        CALL RFFTI(L,WSAVE)
C
      ENDIF
C
C-----------------------------------------------------------------------
C ALLOCATE VARIABLES UNSAVED LOCAL VARIABLES (DEALLOCATED AT THE END!)
C-----------------------------------------------------------------------      
      ALLOCATE(DPY(L))

      IF(IPRES==2) DP=DP*DELZSQ
C-----------------------------------------------------------------------
C     APPLY FFT IN THE Y DIRECTION
C-----------------------------------------------------------------------
      IF(ISWAP==0) THEN
        IF(IPRES==6) THEN
          DO K = 1, NL
          DO IPROC=1,NPROCJ
          DO I = ILU(IPROC),IMU(IPROC)
            DPY(1:L) = DP(I+1,2:L+1,K+1)
            CALL RFFTF(L,DPY,WSAVE)
            DO JJ=1,L
              J = JORDER(JJ)
              DUMMY(K+(I-ILU(IPROC))*NL+(J-1)*(IMU(IPROC)-ILU(IPROC)+1)*NL+NL*L*(ILU(IPROC)-1))=DPY(JJ)
            ENDDO
          ENDDO
          ENDDO
          ENDDO
        ELSEIF(IPRES==5 .AND. ICR==1) THEN
          DO K = 1, NL
          DO I = 1, M
            DPY(1:L) = DP(I+1,2:L+1,K+1)
            CALL RFFTF(L,DPY,WSAVE)
            DO J=1, L
              DUMMY(I+(K-1)*M+(J-1)*M*NL)=DPY(J)
            ENDDO
          ENDDO
          ENDDO
        ELSE
          DO K = 1, NL
          DO I = 1, M
            DPY(1:L) = DP(I+1,2:L+1,K+1)
            CALL RFFTF(L,DPY,WSAVE)
            DO J=1, L
              DUMMY(K+(I-1)*NL+(J-1)*M*NL)=DPY(J)
            ENDDO
          ENDDO
          ENDDO
        ENDIF
      ELSE
        ALLOCATE (DUMMY1(M,L,NL),DUMMY2(M,LL,N))
        DO K = 1, NL
        DO I = 1, M
          DPY(1:L) = DP(I+1,2:L+1,K+1)
          CALL RFFTF(L,DPY,WSAVE)
          DO J=1, L
            DUMMY(I+(J-1)*M+(K-1)*M*L)=DPY(J)
            DUMMY1(I,J,K) = DPY(J)
          ENDDO
        ENDDO
        ENDDO
      ENDIF
c
C-----------------------------------------------------------------------
c     SWAP
C-----------------------------------------------------------------------
      CLOCK(1) = tclock()

      IF(ISWAP==0) THEN
        CALL MPI_ALLTOALLV(DUMMY,SNDCOUNTS,SDISPL,MTYPE,
     &        DUMMYT,RCVCOUNTS,RDISPL,MTYPE,MPI_COMM_EDDY,IERROR)
      ELSEIF(ISWAP==1) THEN

        ALLOCATE(WRK1(M,LL,NL),WRK2(M,LL,NL))

        DO II=0,MYSIZE-1
          IP=SENDRECVPROC(II)
         DO JL=1,LL
          J = IP*LL+JL
          DO I=1,M
          DO K=1,NL
            WRK1(I,JL,K) = DUMMY1(I,J,K)
          ENDDO
          ENDDO
         ENDDO
          CALL MPI_SENDRECV(WRK1(1,1,1),M*LL*NL,MTYPE,IP,JL,
     &                      WRK2(1,1,1),M*LL*NL,MTYPE,IP,JL,MPI_COMM_EDDY,STATUS,IERROR)
         DO JL=1,LL
          DO I=1,M
          DO K=1,NL
            DUMMY2(I,JL,K+IP*NL) = WRK2(I,JL,K)
          ENDDO
          ENDDO
         ENDDO
        ENDDO

        DO I=1,M
        DO JL=1,LL
        DO K=1,N
          DUMMYT(I+(JL-1)*M+(K-1)*M*LL) = DUMMY2(I,JL,K)
        ENDDO
        ENDDO
        ENDDO
        DEALLOCATE(WRK1,WRK2)

      ENDIF

      CLOCK(1) = tclock() - CLOCK(1)

C-----------------------------------------------------------------------
C     FROM NOW ON VARIABLES ARE IN THE FOURIER SPACE
C-----------------------------------------------------------------------
C     INVERT THE SYSTEM OF EQUATIONS FOR EACH L AND DETERMINE THE 
C     UNKNOWN FOURIER COEFFICIENTS OF THE TRANSFORMED PRESSURE. 
C     THE COEFFICIENTS ARE NY-2, EACH DETERMINED BY THE COMPOSITION 
C     OF A COUPLE OF THE VALUES STORED IN DUMMY, SO THAT:
C     QR(0)          = DUMMY(1)
C     QR(1)          = DUMMY(2)+I*DUMMY(3)
C     ...
C     QR((NY-2)/2-1) = DUMMY((NY-2)-2)+I*DUMMY((NY-2)-1)
C     QR((NY-2)/2  ) = DUMMY( NY-2)
C     QR((NY-2)/2+1) = DUMMY((NY-2)-2)-I*DUMMY((NY-2)-1)
C     ...
C     QR((NY-2)-1)   = DUMMY(2)-I*DUMMY(3)
C     QR((NY-2)  )   = DUMMY(1)
C     
C-----------------------------------------------------------------------

      DO JL = 1, LL             ! LOCAL INDEX

        IF(IPRES==5) THEN
          J = JPL 
        ELSEIF(IPRES==6) THEN
          J = JPLANE(JL)
        ELSE
          J=LL*MYRANK+JL          ! GLOBAL INDEX
        ENDIF
C     
C-----------------------------------------------------------------------
C     CHANGE THE DIAGONAL COEFFICIENTS
C-----------------------------------------------------------------------
c
C     ASSIGN THE RHS
        IF(IPRES==6) THEN
          IF(ICR==0) THEN
            DO IPROC=0,MYSIZE-1
            DO K=1,NL
              KG = K+IPROC*NL
            DO I=ILF,IUF
              RHS1D(kg+n*(i - ilf)) = DUMMYT(K+(I-ILF)*NL+(JL-1)*(IUF-ILF+1)*NL+IPROC*(IUF-ILF+1)*NL*NSETS)
            ENDDO
            ENDDO
            ENDDO
          ENDIF
        ELSEIF(IPRES==5) THEN
          IF(ICR==0) THEN
            DO IPROC=0,MYSIZE-1
            DO K=1,NL
              KG = K+IPROC*NL
            DO I=ILF,IUF
              RHS1D(kg+n*(i - ilf)) = DUMMYT(K+(I-ILF)*NL+IPROC*(IUF-ILF+1)*NL)
            ENDDO
            ENDDO
            ENDDO
          ELSE
            DO K=ILF,IUF
            DO I=1,M
              RHS1D(i+m*(k - ilf)) = DUMMYT(I+(K-ILF)*M)
            ENDDO
            ENDDO
          ENDIF
        ELSEIF(IPRES/=3) THEN
          if((icr==0 .AND. ipres==4) .OR. (icr==1 .AND. ipres<4)) then
            DO IPROC=0,MYSIZE-1
            DO K=1,NL
              KG = K + IPROC*NL
            DO I=1,M
              RHS1D(kg+n*(i-1)) = DUMMYT(K + (I-1)*NL + (JL-1)*M*NL + IPROC*M*NL*LL)
            ENDDO
            ENDDO
            ENDDO
          else
            DO IPROC=0,MYSIZE-1
            DO K=1,NL
              KG = K + IPROC*NL
            DO I=1,M
              RHS1D(i+m*(kg-1)) = DUMMYT(K + (I-1)*NL + (JL-1)*M*NL + IPROC*M*NL*LL)
            ENDDO
            ENDDO
            ENDDO
          endif
        ELSE !IPRES=3
          DO IPROC=0,MYSIZE-1
          DO K=1,NL
            KG = K + IPROC*NL
            DO I=1,M
              RHSP(I,KG) = DUMMYT(K + (I-1)*NL + (JL-1)*M*NL + IPROC*M*NL*LL)
            ENDDO
          ENDDO
          ENDDO
        ENDIF

C-----------------------------------------------------------------------
C     CALL THE SOLVER
C-----------------------------------------------------------------------

        IF(IPRES==1) THEN
          clocktemp = tclock()

          if(icr==0) then
            CALL BLKTRI(1,NP,N,AN,BN,CN,MP,M,AM,BML(:,JL),CM,M,RHS1D,IERROR,W)
          else
            CALL BLKTRI(1,MP,M,AM,BML(:,JL),CM,NP,N,AN,BN,CN,N,RHS1D,IERROR,W2(JL,:))
          endif

          clock(2) = clock(2) + tclock() - clocktemp
        ELSEIF(IPRES==4 .OR. IPRES==5 .OR. IPRES==6) THEN
          init(1) = .false.
          init(2) = .false.
          init(3) = .false.

          IF(icr==0) THEN

            do im=ilf,iuf
            do in=1,n
              rhs1d(in+n*(im - ilf)) = -cm(im)*cn(in)*rhs1d(in+n*(im - ilf))
            end do
            end do

            clocktemp = tclock()
            CALL pdc2dn(m,n,rhs1d,n,ilf,iuf,am,bml(:,jl),cm,an,bn,cn,ch,
     &            dw(jl,:),ldw,iw(jl,:),liw,comm,init,ierr)
            clock(2) = clock(2) + tclock() - clocktemp

          ELSE

            do in=ilf,iuf
            do im=1,m
              rhs1d(im+m*(in - ilf)) = -cm(im)*cn(in)*rhs1d(im+m*(in - ilf))
            end do
            end do

            clocktemp = tclock()
            CALL pdc2dn(n,m,rhs1d,m,ilf,iuf,an,bn,cn,am,bml(:,jl),cm,ch,
     &            dw(1,:),ldw,iw(1,:),liw,comm,init,ierr)

            clock(2) = clock(2)+ tclock() - clocktemp

          ENDIF

        ELSEIF(IPRES==2) THEN
          CALL GENBUN(NP,N,MP,M,AM,BML(:,JL),CM,M,RHSP,IERROR,W)
        ELSEIF(IPRES==3) THEN

          ALLOCATE(BMM(M),CMM(M),DPX(M),DPZ(N))

          DO I=1,M
            DPZ(:)=RHSP(I,:)
            CALL RFFTF(N,DPZ,WSAVF)
            RHSP(I,:)=DPZ(:)
          ENDDO

          DO K=1,N
            BMM=BML(:,JL)-AL(K/2+1)
            DPX(:)=RHSP(:,K)
            IF(MP==1) THEN
              CMM=CM
              IF(J==1.AND.K==1) THEN
                CMM(1)=0.
                DPX(1)=0.
              ENDIF
              CALL TRIDAG(AM,BMM,CMM,DPX,DPX,M)
            ELSE
              CALL CYCLIC(AM,BMM,CM,CM(M),AM(1),DPX,DPX,M)
            ENDIF

            RHSP(:,K)=DPX(:)

          ENDDO

          DO I=1,M
            DPZ(:)=RHSP(I,:)/N
            CALL RFFTB(N,DPZ,WSAVF)
            RHSP(I,:)=DPZ(:)
          ENDDO

          DEALLOCATE(BMM,CMM,DPX,DPZ)
          
        ENDIF

        IF(IERROR/=0) WRITE(6,*) ' J,IERROR = ',J,IERROR
C-----------------------------------------------------------------------
C     FIX THE SUM OF PRESSURE CORRECTION OF WHOLE FIELD TO ZERO
C-----------------------------------------------------------------------

        IF(IPRES==6) THEN
          IF(J==1) THEN
            TMP = SUM(RHS1D(:))
            CALL MPI_ALLREDUCE(TMP,DPM,1,MTYPE,MPI_SUM,COMMJ1,IERR)
            RHS1D = RHS1D - DPM/REAL(M*N)
          ENDIF
        ELSEIF(IPRES==5) THEN
          IF(JPL==1) THEN
            TMP = SUM(RHS1D(:))
            CALL MPI_ALLREDUCE(TMP,DPM,1,MTYPE,MPI_SUM,COMMJ1,IERR)
            RHS1D = RHS1D - DPM/REAL(M*N)
          ENDIF
        ELSEIF(IPRES==3) THEN
          IF(J==1) THEN
            DPM = SUM(RHSP(:,:))/REAL(M*N)
            RHSP = RHSP-DPM
          ENDIF
        ELSE
          IF(J==1) THEN
            DPM = SUM(RHS1D(1:M*N))/REAL(M*N)
            RHS1D = RHS1D-DPM
          ENDIF
        ENDIF
C-----------------------------------------------------------------------
C     DUMP THE SOLUTION INTO DUMMY
C-----------------------------------------------------------------------
        IF(ISWAP==0) THEN
          IF(IPRES==6) THEN
            IF(icr==0) THEN
              DO IPROC=0,MYSIZE-1
              DO K=1,NL
                KG = K+IPROC*NL
                DO I=ILF,IUF
                  DUMMYT(K+(I-ILF)*NL+(JL-1)*(IUF-ILF+1)*NL+IPROC*(IUF-ILF+1)*NL*NSETS) = RHS1D(kg+n*(i - ilf))
                ENDDO
              ENDDO
              ENDDO
            ENDIF
          ELSEIF(IPRES==5) THEN
            if(icr==0) then
              DO IPROC=0,MYSIZE-1
              DO K=1,NL
                KG = K+IPROC*NL
              DO I=ILF,IUF
                DUMMYT(K+(I-ILF)*NL+IPROC*(IUF-ILF+1)*NL) = RHS1D(kg+n*(i - ilf))
              ENDDO
              ENDDO
              ENDDO
            else
              do k=ilf,iuf
              do I=1,m
                DUMMYT(I+(K-ILF)*M) = RHS1D(i+m*(k - ilf))
              end do
              end do
            endif
          ELSEIF(IPRES/=3) THEN
            if((icr==0 .AND. ipres==4) .OR. (icr==1 .AND. ipres<4)) then
              DO IPROC=0,MYSIZE-1
              DO K=1,NL
                KG = K + IPROC*NL
              DO I=1,M
                DUMMYT(K + (I-1)*NL + (JL-1)*M*NL + IPROC*M*NL*LL) = RHS1D(kg+n*(i-1))
              ENDDO
              ENDDO
              ENDDO
            else
              DO IPROC=0,MYSIZE-1
              DO K=1,NL
                KG = K + IPROC*NL
              DO I=1,M
                DUMMYT(K + (I-1)*NL + (JL-1)*M*NL + IPROC*M*NL*LL) = rhs1d(i+m*(kg-1))
              ENDDO
              ENDDO
              ENDDO
            endif
          ELSE !IPRES=3
            DO IPROC=0,MYSIZE-1
            DO K=1,NL
              KG = K + IPROC*NL
              DO I=1,M
                DUMMYT(K + (I-1)*NL + (JL-1)*M*NL + IPROC*M*NL*LL) = rhsp(i,kg)
              ENDDO
            ENDDO
            ENDDO
          ENDIF
        ELSE
          if(icr==0) then
            DO K=1,N
            DO I=1,M
              DUMMYT(I+(JL-1)*M+(K-1)*M*LL) = RHSP(I,K)
              DUMMY2(I,JL,K) = RHSP(I,K)
            ENDDO
            ENDDO
          else
            DO K=1,N
            DO I=1,M
              DUMMYT(I+(JL-1)*M+(K-1)*M*LL) = RHSP(K,I)
              DUMMY2(I,JL,K) = RHSP(K,I)
            ENDDO
            ENDDO
          endif
        ENDIF

C
      ENDDO
C
C-----------------------------------------------------------------------
C     SWAP BACK
C-----------------------------------------------------------------------
      CLOCK(3) = tclock()

      IF(ISWAP==0) THEN
        CALL MPI_ALLTOALLV(DUMMYT,RCVCOUNTS,RDISPL,MTYPE,
     &        DUMMY,SNDCOUNTS,SDISPL,MTYPE,MPI_COMM_EDDY,IERROR)
      ELSEIF(ISWAP==1) THEN

        ALLOCATE(WRK1(M,LL,NL),WRK2(M,LL,NL))

        DO II=0,MYSIZE-1
          IP=SENDRECVPROC(II)

         DO JL=1,LL
          DO I=1,M
          DO K=1,NL
            WRK2(I,JL,K) = DUMMY2(I,JL,K+IP*NL)
          ENDDO
          ENDDO
         ENDDO

          CALL MPI_SENDRECV(WRK2(1,1,1),M*LL*NL,MTYPE,IP,JL,
     &                      WRK1(1,1,1),M*LL*NL,MTYPE,IP,JL,MPI_COMM_EDDY,STATUS,IERROR)

         DO JL=1,LL
          J = IP*LL+JL
          DO I=1,M
          DO K=1,NL
            DUMMY1(I,J,K) = WRK1(I,JL,K)
          ENDDO
          ENDDO
         ENDDO

        ENDDO


        DO K = 1, NL
        DO I = 1, M
        DO J=1, L
          DUMMY(I+(J-1)*M+(K-1)*M*L)=DUMMY1(I,J,K)
        ENDDO
        ENDDO
        ENDDO

        DEALLOCATE(WRK1,WRK2)

      ENDIF
C
      CLOCK(3) = tclock()-CLOCK(3)
C
      DP=0.
      DUMMY=DUMMY/REAL(L)
C-----------------------------------------------------------------------
C     TRANSFORM BACK INTO PHYSICAL SPACE.
C-----------------------------------------------------------------------
      IF(IPRES==6) THEN
        IF(ICR==0) THEN
          DO K = 1, NL
          DO IPROC=1,NPROCJ
          DO I = ILU(IPROC),IMU(IPROC)
            DO JJ=1,L
              J = JORDER(JJ)
              DPY(JJ) = DUMMY(K+(I-ILU(IPROC))*NL+(J-1)*(IMU(IPROC)-ILU(IPROC)+1)*NL+NL*L*(ILU(IPROC)-1))
            ENDDO
            CALL RFFTB(L,DPY,WSAVE)
            DP(I+1,2:L+1,K+1)=DPY(1:L)
          ENDDO
          ENDDO
          ENDDO
        ENDIF
      ELSEIF(IPRES==5 .AND. ICR==1) THEN
        DO K = 1, NL
        DO I = 1, M
          DO J=1,L
            DPY(J) = DUMMY(I+(K-1)*M+(J-1)*M*NL)
          ENDDO
          CALL RFFTB(L,DPY,WSAVE)
          DP(I+1,2:L+1,K+1)=DPY(1:L)
        ENDDO
        ENDDO
      ELSEIF(IPRES==5 .AND. ICR==0) THEN
        DO K = 1, NL
        DO I = 1, M
          CALL RFFTF(L,DPY,WSAVE)
          DO J=1, L
            DPY(J) = DUMMY(K+(I-1)*NL+(J-1)*M*NL)
          ENDDO
          CALL RFFTB(L,DPY,WSAVE)
          DP(I+1,2:L+1,K+1) = DPY(1:L)
        ENDDO
        ENDDO
      ELSE
        DO K = 1, NL
        DO I = 1, M
          DO J=1,L
            DPY(J) = DUMMY(K+(I-1)*NL+(J-1)*M*NL)
          ENDDO
          CALL RFFTB(L,DPY,WSAVE)
          DP(I+1,2:L+1,K+1)=DPY(1:L)
        ENDDO
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
C     CYCLIC B.C. IN X DIRECTION
C-----------------------------------------------------------------------
      IF(MP==0) THEN
        DP(1 ,:,:) = DP(NX-1,:,:)
        DP(NX,:,:) = DP(2   ,:,:)
      ELSE
        DP(1 ,:,:) = DP(2   ,:,:)
        DP(NX,:,:) = DP(NX-1,:,:)
      ENDIF
C-----------------------------------------------------------------------
C     CYCLIC B.C. IN Y DIRECTION
C-----------------------------------------------------------------------
      IF(LP==0) THEN
        DP(:,1 ,:) = DP(:,NY-1,:)
        DP(:,NY,:) = DP(:,2   ,:)
      ELSE
        DP(:,1 ,:) = DP(:,2   ,:)
        DP(:,NY,:) = DP(:,NY-1,:)
      ENDIF

      IF(ISWAP==1) THEN
        DEALLOCATE(DUMMY1,DUMMY2)
      ENDIF
      DEALLOCATE(DPY)

      IF(IPRES==2 .AND. IPRES==3) THEN
        DEALLOCATE(RHSP)
      ENDIF

      RETURN
      END


      SUBROUTINE CLOSESTRANK(SENDRECVRANK,MYRANK,MYSIZE)


      INTEGER MYRANK,MYSIZE
      INTEGER SENDRECVRANK(0:MYSIZE-1)

      INTEGER LEFT,RIGHT,MINNEIGHBORS,I

      LEFT = MYRANK
      RIGHT = MYSIZE-1-MYRANK

      MINNEIGHBORS = MIN(LEFT,RIGHT)
      SENDRECVRANK(0) = MYRANK
      DO I=1,MINNEIGHBORS
        SENDRECVRANK(2*I-1) = MYRANK+I
        SENDRECVRANK(2*I  ) = MYRANK-I
      ENDDO

c      IF(MYRANK.EQ.1) THEN
c        WRITE(6,*) LEFT,RIGHT,MINNEIGHBORS
c        WRITE(6,*) SENDRECVRANK
c      ENDIF

      I=MINNEIGHBORS+1
      IF(LEFT>RIGHT) THEN
        DO WHILE(MYRANK-I>=0)
          SENDRECVRANK(MINNEIGHBORS+I) = MYRANK-I
          I=I+1
        ENDDO
      ELSE
        DO WHILE(MYRANK+I<=MYSIZE-1)
          SENDRECVRANK(MINNEIGHBORS+I)= MYRANK+I
c          IF(MYRANK.EQ.0) WRITE(6,*) MINNEIGHBORS+I,MYRANK+I
          I=I+1
        ENDDO
      ENDIF

c      IF(MYRANK.EQ.1) THEN
c        DO I=0,MYSIZE-1
c          WRITE(6,'(2(1X,I5))') I,SENDRECVRANK(I)
c        ENDDO
c      ENDIF

      RETURN

      END
