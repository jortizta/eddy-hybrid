C-----SUBROUTINE-------------------------------------------------------
C
C     PURPOSE: Reads an ASCII or BINARY stl file of the immersed object
C
      SUBROUTINE READSTL(UNVECT,VERTEX,VERTEXC,AREAF,TRINO,NFACET,ZC,NZ,NBD,MBD,IFLAG)
C
c
c outputs
c UNVECT: outward normal
c VERTEX: vertices of the triangles
c VERTEXC: centers of the triangles
c AREAF: areas of the triangles
c TRINO: face numbers
c
c inputs
c IFLAG: 0 (binary input file) or > 0 (formatted input file)
c
C---------------------------------------------------------------------
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
      INTEGER NFACET,NBD,IFLAG,MBD,NZ
      REAL    DSMIN,DSMAX,DSAVE
      REAL    ZC(NZ)
      INTEGER TRINO(NFACET)
      REAL*4, DIMENSION(:),ALLOCATABLE :: DS

      REAL    UNVECT(3,NFACET),VERTEX(3,3,NFACET),VERTEXC(3,NFACET)
      REAL    AREAF(NFACET)
      CHARACTER FILENAME*80

c..local arrays
      INTEGER i,ibd,ifacet,ilb,ile,nfctmax,nfctmin
      INTEGER iun,ivtx,ifct,nfct(nbd)
      INTEGER IVERTEX
      REAL    a(3),b(3),c(3)
      REAL    vecmag,dotprod
      REAL*4  UNVECT1(3,NFACET),VERTEX1(3,3,NFACET),VERTEX1C(3,NFACET)
      REAL*4  UN(3),VC(3)

      nfct =0

      IF(MYSIZE==1 .OR. IMBDCMP==0) THEN

c        write(6,*) myrank,'inside readstl'

        do ibd=1,nbd

          ilb = lb(ibd)+1
          ile = lb(ibd)+mb(ibd)

          if(myrank==0) then

            IF(iflag>0) THEN

              filename = solid(ibd)
              if(ibd>=mbd) then
                if(ifield>0) then
                  filename = 'res.mbd.'//index(ibd)
                endif
              endif

              open(unit=20,file=trim(filename),form='formatted'
     &         ,status='old')
              read(20,*) ifacet
              WRITE(6,*) IFACET,ILB,ILE,trim(filename)
              do ifacet=ilb,ile
                trino(ifacet) = ifacet
                READ(20,*) UNVECT1(1,IFACET),UNVECT1(2,IFACET),UNVECT1(3,IFACET)
                READ(20,*) VERTEX1C(1,IFACET),VERTEX1C(2,IFACET),VERTEX1C(3,IFACET)
                READ(20,*) VERTEX1(1,1,IFACET),VERTEX1(2,1,IFACET),VERTEX1(3,1,IFACET)
                READ(20,*) VERTEX1(1,2,IFACET),VERTEX1(2,2,IFACET),VERTEX1(3,2,IFACET)
                READ(20,*) VERTEX1(1,3,IFACET),VERTEX1(2,3,IFACET),VERTEX1(3,3,IFACET)
              enddo
              close(20)

              UNVECT = UNVECT1
              VERTEX = VERTEX1
              VERTEXC = VERTEX1C

            ELSE

              filename = solid(ibd)
              if(ibd>=mbd) then
                if(ifield>0) then
                  filename = 'res.mbd.'//index(ibd)
                endif
              endif

              open(unit=20,file=trim(filename),form='UNFORMATTED'
     &            ,status='OLD')
              READ(20) IFACET
              WRITE(6,*) IFACET,ILB,ILE
              DO IFACET=ILB,ILE
                trino(ifacet) = ifacet
                READ(20) UNVECT1(1,IFACET),UNVECT1(2,IFACET),UNVECT1(3,IFACET)
                READ(20) VERTEX1C(1,IFACET),VERTEX1C(2,IFACET),VERTEX1C(3,IFACET)
                READ(20) VERTEX1(1,1,IFACET),VERTEX1(2,1,IFACET),VERTEX1(3,1,IFACET)
                READ(20) VERTEX1(1,2,IFACET),VERTEX1(2,2,IFACET),VERTEX1(3,2,IFACET)
                READ(20) VERTEX1(1,3,IFACET),VERTEX1(2,3,IFACET),VERTEX1(3,3,IFACET)
              ENDDO
              CLOSE(20)
  
              UNVECT = UNVECT1
              VERTEX = VERTEX1
              VERTEXC = VERTEX1C

c              write(6,*) 'object:',unvect(:,1307)

            ENDIF

          ENDIF
        ENDDO

        CALL MPI_BCAST(UNVECT,NFACET*3,MTYPE,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(VERTEX,NFACET*9,MTYPE,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(VERTEXC,NFACET*3,MTYPE,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(TRINO,NFACET,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)    !!!!!!
c        CALL MPI_BCAST(AREAF,NFACET,MTYPE,0,MPI_COMM_EDDY,IERR)

      ELSE

        DO ibd=1,nbd

          ilb = lb(ibd)+1
          ile = lb(ibd)+mb(ibd)

          filename = solid(ibd)
          if(ibd>=mbd) then
            if(ifield>0) then
              filename = 'res.mbd.'//index(ibd)
            endif
          endif

          open(unit=20,file=trim(filename),form='UNFORMATTED'
     &         ,status='OLD')
          READ(20) IFACET

          DO IFACET=ILB,ILE
            READ(20) UN(1),UN(2),UN(3)
            READ(20) VC(1),VC(2),VC(3)
            IF(VC(3)>=ZC(1)-IMBOVLP .AND. VC(3)<ZC(NZ-1)+IMBOVLP) THEN
              nfct(ibd) = nfct(ibd)+1
              ifct = nfct(ibd)
              IF(ifacet>nfacet) THEN
                write(6,*) 'ERROR: Increase nfacet'
                call mpi_finalize(ierr)
                stop
              ENDIF
              TRINO(ifct)=ifacet
              UNVECT1(:,ifct) = UN(:)
              VERTEX1C(:,ifct) = VC(:)
              READ(20) VERTEX1(1,1,ifct),VERTEX1(2,1,ifct),VERTEX1(3,1,ifct)
              READ(20) VERTEX1(1,2,ifct),VERTEX1(2,2,ifct),VERTEX1(3,2,ifct)
              READ(20) VERTEX1(1,3,ifct),VERTEX1(2,3,ifct),VERTEX1(3,3,ifct)
            ELSE
              READ(20)
              READ(20)
              READ(20)
            ENDIF
          ENDDO

        ENDDO

        lb=0
        mb(1:nbd)=nfct(1:nbd)

        DO IBD=2,NBD
          lb(ibd)=lb(ibd-1)+mb(ibd-1)
        ENDDO
        nfacet = sum(mb(1:nbd))

        unvect(:,1:nfacet)    = unvect1(:,1:nfacet)
        vertex(:,:,1:nfacet)  = vertex1(:,:,1:nfacet)
        vertexc(:,1:nfacet)   = vertex1c(:,1:nfacet)

        CALL MPI_REDUCE(NFACET,NFCTMAX,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_EDDY,IERR)
        CALL MPI_REDUCE(NFACET,NFCTMIN,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_EDDY,IERR)

        IF(MYRANK.eq.0) THEN
          write(6,*) 'Max triangles in subdomain:',nfctmax,'min',nfctmin
        ENDIF

      ENDIF
  
   
      DO ibd=1,nbd

        ilb=lb(ibd)+1
        ile=lb(ibd)+mb(ibd)

        !Calculate area
        do ifacet=ilb,ile
          do i=1,3
            a(i) = vertex(i,1,ifacet)-vertex(i,2,ifacet)
            b(i) = vertex(i,3,ifacet)-vertex(i,2,ifacet)
          enddo

          call cross(c,a,b)
          areaf(ifacet) = 0.5*sqrt(c(1)**2.+c(2)**2.+c(3)**2.)
        enddo

        !Calculate normal
        IF(.false.) THEN
          do ifacet=ilb,ile

            do i=1,3
!!!!!!              a(i) = vertex(i,1,ifacet)-vertex(i,2,ifacet)
              a(i) = vertex(i,2,ifacet)-vertex(i,1,ifacet)
              b(i) = vertex(i,3,ifacet)-vertex(i,2,ifacet)
            enddo
           
            CALL cross(c,a,b)

!!!!!!            a(1) = vertexc(1,ifacet)
!!!!!!            a(2) = vertexc(2,ifacet)
!!!!!!            a(3) = vertexc(3,ifacet)
!!!!!!            if(dotprod(c,a)<0) c=-c
c             c=a
!!!!!!            c(2)=0.0
            c=c/vecmag(c,3)
            unvect(1,ifacet)=c(1)
            unvect(2,ifacet)=c(2)
            unvect(3,ifacet)=c(3)
          enddo
        ENDIF

      ENDDO

c      write(6,*) '2. object:',unvect(:,1307)


      RETURN

      END
C---------------------------------------------------------------------


C---- subroutine readtrvtx -------------- N. Beratlis-31 Oct. 2010 ---
C
C     PURPOSEL Read vertex coordinates of triangles and connectivity.
C
C---------------------------------------------------------------------
      subroutine readtrvtx(vertex,norm,vtxno,nv,nbd,mbd)   !!!!!! nbd,mbd
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nv,nbd,mbd   !!!!!!
      integer vtxno(nv)
      real    vertex(3,nv),norm(3,nv)
c
c.... Local arrays
      integer i,ibd,ilb,ile   !!!!!!
      CHARACTER FILENAME*80   !!!!!!

!!!!!!
      do ibd=1,nbd

        ilb = lv(ibd)+1
        ile = lv(ibd)+mv(ibd)

        if((ibd>=mbd).and.(ifield>0)) then
           filename = 'res.mbd_nodes_conn.'//index(ibd)
        else
           filename = 'vtx_conn.'//index(ibd)//'.bin'
        endif

        open(unit=10,file=trim(filename),form='unformatted')

        do i=ilb,ile
           read(10) vertex(:,i),norm(:,i)
           vtxno(i) = i
        enddo

        close(10)

      enddo
!!!!!!

      RETURN

      END
C---------------------------------------------------------------------


C---- subroutine writetrvtx ------------- N. Beratlis-31 Oct. 2010 ---
C
C     PURPOSEL Write vertex coordinates of triangles and connectivity.
C
C---------------------------------------------------------------------
      subroutine writetrvtx(filename,vertex,norm,nv,ibd)   !!!!!! ibd
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nv,ibd   !!!!!! ibd
      real    vertex(3,nv),norm(3,nv)
      character*(*) filename
c
c.... Local arrays
      integer i,ilb,ile   !!!!!! ilb,ile

      ilb = lv(ibd)+1   !!!!!!
      ile = lv(ibd)+mv(ibd)   !!!!!!

      open(unit=10,file=trim(filename),form='unformatted')
!!!!!!      write(10) mv(ibd)   !!!!!! instead of nv
      do i=ilb,ile   !!!!!! instead of 1,nv
        write(10) vertex(:,i),norm(:,i)
      enddo
      close(10)

      RETURN

      END
C---------------------------------------------------------------------


C---- subroutine iompi_imb------------------N. Beratlis-5 Apr. 2010---
C
      subroutine iompi_imb(unvect,vertex,vertexc,areaf,trino,nfacet
     &     ,zc,nz,nbd,mbd,mbg,iflag)
C
C     PURPOSE: Read or write file contating immersed body information.
C
C---------------------------------------------------------------------
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
      integer nfacet,iflag,nz,ibd,mbd,nbd,mbg(nbd)
      integer trino(nfacet)
      real    zc(nz)
      real    areaf(nfacet)
      real    unvect(3,nfacet),vertex(3,3,nfacet),vertexc(3,nfacet)
c
c.... Local arrays
      integer i,ilb,ile,ifacet,nb,ifct,nfctmax,nfctmin
      integer nfct(nbd)
      character*120 filename,filein
      real    un(3),vc(3),a(3),b(3),c(3),time
      integer, dimension(:), allocatable :: ntr,ntrg
      integer fh,filemode
      integer newtype,newtype2
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      nfct = 0

      call mpi_type_contiguous(3,mtype,newtype,ierr)
      call mpi_type_commit(newtype,ierr)

      if(iflag==0) then

        if(imbdcmp==0 .OR. mysize==1) then

          do ibd=1,nbd

            filename = solid(ibd)
            if(ibd>=mbd) then
              if(ifield>0) then
                filename = 'res.mbd.'//index(ibd)//'.mpi.bin'
              endif
            endif

            filename = trim(str6)//trim(filename)
c            write(6,*) 'str6=',trim(str6)
c            filename = 'mbd.'//index(ibd)//'.mpi.bin.res'

c            write(6,*) 'filename=',filename
c            write(6,*) 'filein=',filein
          
            filemode = MPI_MODE_RDONLY
            call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

            disp = 0
            call mpi_file_set_view(fh,disp,mpi_double_precision,newtype,'native',mpi_info_null,ierr)
            call mpi_file_read(fh,un,1,newtype,status,ierr)
            i = un(1)
c            write(6,*) 'i=',i,',mbg=',mbg(ibd)
            if(i/=mbg(ibd)) then
              write(6,*) 'WARNING:No. of triangles incorrect in header',filename
            endif

            ilb=lb(ibd)+1
            ile=lb(ibd)+mb(ibd)
            do i=ilb,ile
              call mpi_file_read(fh,unvect(:,i),1,newtype,status,ierr)
              call mpi_file_read(fh,vertexc(:,i),1,newtype,status,ierr)
              call mpi_file_read(fh,vertex(:,1,i),1,newtype,status,ierr)
              call mpi_file_read(fh,vertex(:,2,i),1,newtype,status,ierr)
              call mpi_file_read(fh,vertex(:,3,i),1,newtype,status,ierr)
c              write(6,*) i,unvect(:,i),vertexc(:,i)
            enddo
            call mpi_file_close(fh,ierr)

          enddo

        else

          do ibd=1,nbd

            filename = solid(ibd)
            if(ibd>=mbd) then
              if(ifield>0) then
                filename = 'res.mbd.'//index(ibd)//'.mpi.bin'
              endif
            endif

            filename = trim(str6)//trim(filename)

c            filename = 'mbd.'//index(ibd)//'.mpi.bin.res'

            filemode = MPI_MODE_RDONLY
            call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)

            disp = 0
            call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)
            call mpi_file_read(fh,un,1,newtype,status,ierr)

            ilb=lb(ibd)+1
            ile=lb(ibd)+mb(ibd)

c            write(6,*) 'ilb=',ilb,', ile=',ile

            do i=ilb,ile
              call mpi_file_read(fh,un,1,newtype,status,ierr)
              call mpi_file_read(fh,vc,1,newtype,status,ierr)
c              write(6,*) myrank,i,un,vc
              if(vc(3)>=zc(1)-imbovlp .AND. vc(3)<zc(nz-1)+imbovlp) then
                nfct(ibd) = nfct(ibd)+1
                ifct = nfct(ibd)
                if(i>nfacet) then
                  write(6,*) 'ERROR: Increase nfacet'
                  call mpi_finalize(ierr)
                  stop
                endif
                trino(ifct)=i
                unvect(:,ifct) = un(:)
                vertexc(:,ifct) = vc(:)
                call mpi_file_read(fh,vertex(:,1,ifct),1,newtype,status,ierr)
                call mpi_file_read(fh,vertex(:,2,ifct),1,newtype,status,ierr)
                call mpi_file_read(fh,vertex(:,3,ifct),1,newtype,status,ierr)
              else
                call mpi_file_read(fh,un,1,newtype,status,ierr)
                call mpi_file_read(fh,un,1,newtype,status,ierr)
                call mpi_file_read(fh,un,1,newtype,status,ierr)
              endif
            enddo
            call mpi_file_close(fh,ierr)

          enddo

          lb=0
          mb(1:nbd)=nfct(1:nbd)

          do ibd=2,nbd
            lb(ibd)=lb(ibd-1)+mb(ibd-1)
          enddo
          nfacet = sum(mb(1:nbd))

          CALL MPI_REDUCE(NFACET,NFCTMAX,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_EDDY,IERR)
          CALL MPI_REDUCE(NFACET,NFCTMIN,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_EDDY,IERR)

          if(myrank.eq.0) then
            write(6,*) 'Max triangles in subdomain:',nfctmax,'min',nfctmin
          endif
 
        endif
        
        do ibd=1,nbd

          ilb=lb(ibd)+1
          ile=lb(ibd)+mb(ibd)

          !Calculate area
          do ifacet=ilb,ile
            do i=1,3
              a(i) = vertex(i,1,ifacet)-vertex(i,2,ifacet)
              b(i) = vertex(i,3,ifacet)-vertex(i,2,ifacet)
            enddo

            call cross(c,a,b)
            areaf(ifacet) = 0.5*sqrt(c(1)**2.+c(2)**2.+c(3)**2.)
          enddo
        enddo

      elseif(iflag==1) then

        if(imbdcmp==1 .AND. mysize>1) then

          do ibd=1,nbd

            filename = trim(str6)//'mbd.'//index(ibd)//'.mpi.bin.res'
            filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
            call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
            disp = 0
            call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)
 
            if(myrank==0) then
              un(:) = real(mbg(ibd))
              call mpi_file_write(fh,un,1,newtype,status,ierr)
            endif

            CALL MPI_BARRIER(MPI_COMM_EDDY,IERR)

            ilb = lb(ibd)+1
            ile = lb(ibd)+mb(ibd)

            ALLOCATE(ntr(mysize),ntrg(mysize))
            ntr = 0
            do i=ilb,ile
              if(vertexc(3,i)>=zc(1) .AND. vertexc(3,i)<zc(nz-1)) then
                ntr(myrank+1) = ntr(myrank+1)+1
              endif
            enddo
            
            CALL MPI_ALLREDUCE(ntr,ntrg,mysize,MPI_INTEGER,MPI_MAX,MPI_COMM_EDDY,IERR)
          
            offset = sum(ntrg(1:myrank))*3*5+3
            do i=ilb,ile
              if(vertexc(3,i)>=zc(1) .AND. vertexc(3,i)<zc(nz-1)) then
                call mpi_file_write_at(fh,offset,unvect(:,i),1,newtype,status,ierr)                
                call mpi_file_write_at(fh,offset+3,vertexc(:,i),1,newtype,status,ierr)
                call mpi_file_write_at(fh,offset+6,vertex(:,1,i),1,newtype,status,ierr)
                call mpi_file_write_at(fh,offset+9,vertex(:,2,i),1,newtype,status,ierr)
                call mpi_file_write_at(fh,offset+12,vertex(:,3,i),1,newtype,status,ierr)
                offset = offset+15
              endif
            enddo

            DEALLOCATE(ntr,ntrg)

            call mpi_file_close(fh,ierr)
          enddo

        elseif(myrank==0) then
           
          do ibd=1,nbd

            filename = trim(str6)//'mbd.'//index(ibd)//'.mpi.bin.res'
            filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
            call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
            disp = 0
            call mpi_file_set_view(fh,disp,mpi_double_precision,newtype
     &           ,'native',mpi_info_null,ierr)
            un(:) = real(mbg(ibd))
            call mpi_file_write(fh,un,1,newtype,status,ierr)

            ilb = lb(ibd)+1
            ile = lb(ibd)+mb(ibd)
            do i=ilb,ile
              call mpi_file_write(fh,unvect(:,i),1,newtype,status,ierr)
              call mpi_file_write(fh,vertexc(:,i),1,newtype,status,ierr)
              call mpi_file_write(fh,vertex(:,1,i),1,newtype,status,ierr)
              call mpi_file_write(fh,vertex(:,2,i),1,newtype,status,ierr)
              call mpi_file_write(fh,vertex(:,3,i),1,newtype,status,ierr)
c              write(6,*) i,vertexc(:,i)
            enddo
            call mpi_file_close(fh,ierr)

          enddo

        endif
           
      endif

      call mpi_type_free(newtype,ierr)

      return

      end
C---------------------------------------------------------------------


C---- subroutine writestl -----------------N. Beratlis-30 Apr. 2009---
C
C     PURPOSE: Writes a BINARY stl file of the immersed object
C
      subroutine writestl(filename,unvect,vertex,vertexc,zc,nz,nfacet,nfacetot,ibd,iflag)
C
C---------------------------------------------------------------------
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
      integer nfacet,iflag,nz,nfacetot
      real    zc(nz)
      real    unvect(3,nfacet),vertex(3,3,nfacet),vertexc(3,nfacet)
      character*(*) filename

      integer i,ibd,ilb,ile,ifacet,nb
      integer fh,filemode
      integer newtype,newtype2
c      integer nft(mysize)
      INTEGER(KIND=MPI_OFFSET_KIND) disp,offset
      integer status(mpi_status_size)

      ilb = lb(ibd)+1
      ile = lb(ibd)+mb(ibd)

      if(iflag==0) then
        if(myrank.eq.0) then
          open(unit=20,file=trim(filename),form='unformatted')
          write(20) mb(ibd)     !!!!!! instead of nfacet
          do i=ilb,ile
            write(20) real(unvect(1,i),4),real(unvect(2,i),4),real(unvect(3,i),4)
            write(20) real(vertexc(1,i),4),real(vertexc(2,i),4),real(vertexc(3,i),4)
            write(20) real(vertex(1,1,i),4),real(vertex(2,1,i),4),real(vertex(3,1,i),4)
            write(20) real(vertex(1,2,i),4),real(vertex(2,2,i),4),real(vertex(3,2,i),4)
            write(20) real(vertex(1,3,i),4),real(vertex(2,3,i),4),real(vertex(3,3,i),4)
          enddo
          close(20)
        endif
      elseif(iflag==1) then

c a newtype is generated
        call mpi_type_contiguous(3,mtype,newtype,ierr)
        call mpi_type_commit(newtype,ierr)

        filemode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
        call mpi_file_open(mpi_comm_eddy,trim(filename),filemode,mpi_info_null,fh,ierr)
        disp = 0
        call mpi_file_set_view(fh,disp,mtype,newtype,'native',mpi_info_null,ierr)

c the header is written by the processor 0
        if(myrank==0) then
!!!!!!           call mpi_file_write(fh,nfacetot,1,mpi_integer,status,ierr)
           call mpi_file_write(fh,mb(ibd),1,mpi_integer,status,ierr)
        endif

c each processor writes the components of the normal vector (UNVECT),
c the coordinates of the centers of triangles (VERTEXC) and the ones
c of the vertices (VERTEX)
        do i=ilb,ile
          nb=0
          if(vertexc(3,i)>=zc(1) .AND. vertexc(3,i)<zc(nz-1)) then
            call mpi_file_write(fh,unvect(:,i),1,newtype,status,ierr)
            call mpi_file_write(fh,vertexc(:,i),1,newtype,status,ierr)
            call mpi_file_write(fh,vertex(:,1,i),1,newtype,status,ierr)
            call mpi_file_write(fh,vertex(:,2,i),1,newtype,status,ierr)
            call mpi_file_write(fh,vertex(:,3,i),1,newtype,status,ierr)
          endif
        enddo

        call mpi_file_close(fh,ierr)
        

      endif
      
      return

      end
C---------------------------------------------------------------------

C---- subroutine imb2plt ------------------N. Beratlis-30 Apr. 2009---
C
C     PURPOSE: Writes an plt TECPLOT file of the immersed object
c
c N.B. There is no need of communication among processors since all
c of them have information about the immersed-boundary
C
      subroutine imb2plt(filename,vertex,nfacet,ibd)
C
C---------------------------------------------------------------------
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
      integer nfacet
      real    vertex(3,3,nfacet)
      character*(*) filename

      integer i,j,ibd,ilb,ile,ifacet

      ilb = lb(ibd)+1
      ile = lb(ibd)+mb(ibd)

      if(myrank==0) then
        open(unit=20,file=trim(filename),form='formatted')
        write(20,'(A)') 'FILETYPE = Grid'
        write(20,'(A)') 'VARIABLES = "X", "Y", "Z"'
        write(20,'(A,A,I8,A,I8)') 'ZONE '
     &       ,'DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE,N='
     &       ,MB(IBD)*3,',E=',MB(IBD)
        write(20,'(E15.8)') ((VERTEX(1,I,J),I=1,3),J=LB(IBD)+1,LB(IBD)+MB(IBD))
        write(20,'(E15.8)') ((VERTEX(2,I,J),I=1,3),J=LB(IBD)+1,LB(IBD)+MB(IBD))
        write(20,'(E15.8)') ((VERTEX(3,I,J),I=1,3),J=LB(IBD)+1,LB(IBD)+MB(IBD))
        write(20,'(3I8)') ((I-1)*3+1,(I-1)*3+2,(I-1)*3+3,I=1,MB(IBD))
        close(20)
      endif
      
      return

      end
C---------------------------------------------------------------------


c
c---- SUBROUTINE LIMIMB ----------------------------------------------
c
c     PURPOSE: Find limits of subdomain that enclose immersed object
c
      SUBROUTINE LIMIMB(X,Y,Z,ZC,NX,NY,NZ,XYZB,NB,IBD)
c
c---------------------------------------------------------------------
c
      INCLUDE 'common.h'
      INCLUDE 'mpif.h'
      INCLUDE 'immersed.h'

c...Input Arrays
      INTEGER NX,NY,NZ,NB,IBD
      REAL    X(NX),Y(NY),Z(NZ),ZC(NZ),XYZB(3,NB)
c
c...Local arrays
      INTEGER K,IMIN
      REAL, DIMENSION(:), ALLOCATABLE :: RB
      REAL    MINZ,MAXZ

      KBMIN(IBD) = NZ
      KBMAX(IBD) = 1

      IF(ICYL==1) THEN

        ALLOCATE(RB(NB))
        RB = XYZB(1,:)**2. + XYZB(2,:)**2.
        CALL LOCATE(X,NX-1,SQRT(MINVAL(RB)),IBMIN(IBD))
        CALL LOCATE(X,NX-1,SQRT(MAXVAL(RB)),IBMAX(IBD))

        JBMIN(IBD) = 2
        JBMAX(IBD) = NY-1
        DEALLOCATE(RB)

      ELSE

        CALL LOCATE(X,NX,MINVAL(XYZB(1,:)),IBMIN(IBD))
        CALL LOCATE(X,NX,MAXVAL(XYZB(1,:)),IBMAX(IBD))

        JBMIN(IBD) = 2
        JBMAX(IBD) = NY-1

      ENDIF

      MINZ = MINVAL(XYZB(3,:))!-(Z(2)-Z(1))
      MAXZ = MAXVAL(XYZB(3,:))!+(Z(NZ)-Z(NZ-1))
      IF(MINZ<Z(NZ).AND.MAXZ>=ZC(1)) CALL LOCATE(Z,NZ,MINZ,KBMIN(IBD))
      IF(MINZ<Z(NZ).AND.MAXZ>=ZC(1)) CALL LOCATE(Z,NZ,MAXZ,KBMAX(IBD))
c      write(6,*) 'limimb:myrank=',myrank,minz,maxz,z(1),z(2),z(nz-1),z(nz),kbmin(ibd),kbmax(ibd)

      IBMIN(IBD) = IBMIN(IBD)-3
      IBMAX(IBD) = IBMAX(IBD)+3
      IF(IBMIN(IBD)<1 ) IBMIN(IBD)=1 !2
      IF(IBMAX(IBD)>NX-1) IBMAX(IBD)=NX-1

      IF(KBMAX(IBD)>=KBMIN(IBD)) THEN
        KBMIN(IBD) = KBMIN(IBD)-3
        KBMAX(IBD) = KBMAX(IBD)+3
      ENDIF

      IF(MYRANK.EQ.0) THEN
        IF(ITYPE(5)/=500) THEN
          IF(KBMIN(IBD)<1) KBMIN(IBD)=1
        ELSE
          IF(KBMIN(IBD)<2) KBMIN(IBD)=2
        ENDIF
      ELSE
        IF(KBMIN(IBD)<2) KBMIN(IBD)=2
      ENDIF

      IF(MYRANK.EQ.MYSIZE-1) THEN
        IF(ITYPE(6)/=500) THEN
          IF(KBMAX(IBD)>NZ) KBMAX(IBD)=NZ
        ELSE
          IF(KBMAX(IBD)>NZ-1) KBMAX(IBD)=NZ-1
        ENDIF
      ELSE
        IF(KBMAX(IBD)>NZ-1) KBMAX(IBD)=NZ-1
      ENDIF

      IF(IMBDCMP==1 .AND. MYSIZE>1) THEN
        CALL MPI_ALLREDUCE(IBMIN(IBD),IMIN,1,MPI_INTEGER,MPI_MIN,MPI_COMM_EDDY,IERR)
        IBMIN(IBD)=IMIN
      ENDIF

      RETURN

      END
c---------------------------------------------------------------------



c---- SUBROUTINE CALC_DS_TRIANGLES -----------------------------------
c
c     PURPOSE: Evaluation of the max length among the triangles edges
c
      SUBROUTINE CALC_DS_TRIANGLES(VERTEX,VERTEXC,NFACET,ILB,ILE,S,X,NX,ICYL)

      IMPLICIT NONE
C
C... Input/Output arrays
      INTEGER NFACET,ILB,ILE,NX,ICYL
      REAL    S(NX),X(NX)
      REAL    VERTEX(3,3,NFACET),VERTEXC(3,NFACET)
C
C... Local arrays
      INTEGER I,II,JJ,J,I1,I2,ILIM(3)
      REAL    DSMIN,DSMAX,DSAVE,SX,SY,SZ,a,b,c,eps,DS,R,RCYL,THETA
c
c... Functions
      REAL    anglerad

      rcyl = real(icyl)

      S = 0.0
      !Calculate triangle circumradius
      !R = a*b*c/sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c)) from http://mathworld.wolfram.com/Circumradius.html
      DO I=ILB,ILE
        a = sqrt ( (vertex(1,1,i)-vertex(1,2,i))**2.
     &           + (vertex(2,1,i)-vertex(2,2,i))**2. )
c     &           + (vertex(3,1,i)-vertex(3,2,i))**2. )
        b = sqrt ( (vertex(1,1,i)-vertex(1,3,i))**2.
     &           + (vertex(2,1,i)-vertex(2,3,i))**2. )
c     &           + (vertex(3,1,i)-vertex(3,3,i))**2. )
        c = sqrt ( (vertex(1,3,i)-vertex(1,2,i))**2.
     &           + (vertex(2,3,i)-vertex(2,2,i))**2. )
c     &           + (vertex(3,3,i)-vertex(3,2,i))**2. )
c        DS(i-ilb+1) = max(a,b,c)
        DS = max(a,b,c)

c the coordinates of the vertices are always in a Cartesian frame of reference
        do j=1,3
          r = sqrt(vertex(1,j,i)**2.0 + vertex(2,j,i)**2.0)*rcyl
     &        + vertex(1,j,i)*(1.0-rcyl)
          call locate(x,nx,r,ilim(j))
        enddo

c correction for bodies outside the X boundaries of the computational domain
        i1 = minval(ilim)
        i2 = maxval(ilim)
        if(i1<1) i1=1
        if(i2<1) i2=1
        if(i1>nx-1) i1=nx-1
        if(i2>nx-1) i2=nx-1

c for each coordinate along X the value of S is set equal to the maximum length of the
c edges of the triangles involving that coordinate 
        do ii=i1,i2
          if(ds>s(ii)) s(ii)=ds
          if(ds>s(ii+1)) s(ii+1)=ds
        enddo

      ENDDO

      if(icyl==1) s(1) = s(2)

      DSMAX = MAXVAL(S)
      DO I=1,NX
        IF(S(I)==0.0) S(I)=DSMAX
      ENDDO
c      stop
c      DSMIN = MINVAL(DS)
c      DSMAX = MAXVAL(DS)
c      DSAVE = SUM(DS)/REAL(ILE-ILB+1)      
c      S = DSMIN
c      S = DSMAX
c      write(6,*) 'DSMAX=',S,' @',MAXLOC(DS)
c      WRITE(6,*) 'DSMIN=',DSMIN,', DSMAX=',DSMAX,', DSAVE=',DSAVE
c      DEALLOCATE(DS)

      RETURN

      END
c
c---------------------------------------------------------------------

c---- SUBROUTINE MYFIRSTLAST -----------------N. Beratlis-25 Mar 2010-
c
c     PURPOSE: Find first and last processor containing body
c
      SUBROUTINE MYFIRSTLAST(VERTEX,NFACET,ILB,ILE)
      
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'

      INTEGER NFACET,ILB,ILE
      REAL    VERTEX(3,3,NFACET)

      RETURN

      END
c---------------------------------------------------------------------



C---- subroutine filterpbd------------------------N. Beratlis-25 Jun 2009---
C
C     PURPOSE: Filter variable on body using values at neighboring triangles.
C
C---------------------------------------------------------------------------

      subroutine filterpbd(pbd,nngh,inngh,nnodes,nnh,ilb,ile)

      implicit none
      include 'immersed.h'
c
c.... Input/Output arrays
      integer nnodes,nnh,nh
      integer nngh(nnodes,nnh),inngh(nnodes)
      real    pbd(nnodes)
c
c.... Local arrays
      integer i,ilb,ile,iloop,nloops,j,n
      real, dimension(:), allocatable :: pbd1

      allocate(pbd1(nnodes))
      nloops = 5

c      write(6,*) 'inside filterpbd'

      do iloop=1,nloops
        n=0
        do i=ilb,ile
          n=n+1
          nh = inngh(i)
          pbd1(n) = nh*pbd(i)
c          pbd1(n) = pbd(i)
          do j=1,nh
            pbd1(n) = pbd1(n) + pbd(nngh(i,j))
          enddo
          pbd1(n) = pbd1(n)/real(nh+nh)
c          write(6,*) i,pbd1(n),pbd(i)
        enddo
c        pbd1 = 0.0
        pbd(ilb:ile) = pbd1(1:n)
      enddo

      deallocate(pbd1)

      return

      end
C---------------------------------------------------------------------------



c---- SUBROUTINE MOVTRI ----------------------- N. Beratlis- 26 Mar. 2010---
c
c     PURPOSE: Triangles in domain
c
      SUBROUTINE MOVTRI(UNVECT,VERTEX,VERTEXC,AREA,TRINO,NFACET,NFACETMAX,ZC,NZ,MBD,NBD)
c
C---------------------------------------------------------------------------

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      INTEGER NFACET,NFACETMAX,MBD,NBD,NZ
      INTEGER TRINO(NFACETMAX)
      REAL    ZC(NZ)
      REAL    UNVECT(3,NFACETMAX),VERTEX(3,3,NFACETMAX),VERTEXC(3,NFACETMAX)
     &       ,AREA(NFACETMAX)
c
c.... Local arrays
      INTEGER I,J,K,IBD,NTR,NTRL,NTRR,ILB,ILE,RCVL,RCVR,I1,I2,REP,N,S
      REAL    a(3),b(3),c(3),d
      INTEGER STATUS(MPI_STATUS_SIZE)
      INTEGER TRINOL(NFACETMAX),TRINOR(NFACETMAX)
      REAL    VERTEXL(3,3,NFACETMAX),VERTEXR(3,3,NFACETMAX)
      REAL    UNVECTL(3,NFACETMAX),UNVECTR(3,NFACETMAX)
      INTEGER, DIMENSION(:), ALLOCATABLE :: MB1
c
c.... Functions
      REAL    calc_triarea

      ALLOCATE(MB1(NBD))

      mb1=0
c      lb1(1) = lb(mbd)
c      write(6,*) myrank,'inside movtri'

      DO IBD=MBD,NBD

        ntr=0
        ntrl=0
        ntrr=0

        ilb = lb(ibd)+1
        ile = lb(ibd)+mb(ibd)

        DO i=ilb,ile

          IF(VERTEXC(3,I)>=ZC(1)-IMBOVLP .AND. VERTEXC(3,I)<ZC(NZ-1)+IMBOVLP) THEN
            NTR=NTR+1
            MB1(IBD) = MB1(IBD)+1
            VERTEX(:,:,lb(ibd)+NTR) = VERTEX(:,:,I)
            UNVECT(:,lb(ibd)+NTR) = UNVECT(:,I)
c            VERTEXC(:,lb(ibd)+NTR) = VERTEXC(:,I)
            TRINO(lb(ibd)+NTR) = TRINO(I)
          ELSEIF(VERTEXC(3,I)<ZC(1)-IMBOVLP) THEN
            NTRL=NTRL+1
            VERTEXL(:,:,NTRL) = VERTEX(:,:,I)
            UNVECTL(:,NTRL) = UNVECT(:,I)
            TRINOL(NTRL) = TRINO(I)
          ELSEIF(VERTEXC(3,I)>=ZC(NZ-1)+IMBOVLP) THEN
            NTRR=NTRR+1
            VERTEXR(:,:,NTRR) = VERTEX(:,:,I)
            UNVECTR(:,NTRR) = UNVECT(:,I)
            TRINOR(NTRR) = TRINO(I)
          ENDIF

        ENDDO

        if(ntrl>0 .OR. ntrr>0) then
          vertexc(:,lb(ibd)+1:lb(ibd)+ntr)=sum(vertex(:,:,lb(ibd)+1:lb(ibd)+ntr),2)/3.0
          do i=lb(ibd)+1,lb(ibd)+ntr
            area(i) = calc_triarea(vertex(:,:,i))
          enddo
        endif

        CALL MPI_SENDRECV(NTRL,1,MPI_INTEGER,MYLEFT ,10
     &       ,            RCVR,1,MPI_INTEGER,MYRIGHT,10
     &       ,            MPI_COMM_EDDY,STATUS,IERR)

        CALL MPI_SENDRECV(NTRR,1,MPI_INTEGER,MYRIGHT,11
     &       ,            RCVL,1,MPI_INTEGER,MYLEFT ,11
     &       ,            MPI_COMM_EDDY,STATUS,IERR)

c        open(unit=10,file='movtri.'//index(myrank),form='formatted')
c        write(10,*) 'myrank=',myrank,nfacet,mb1(ibd)
c     &       ,'myleft=',myleft,',myright=',myright
c     &       ,'send to left:',ntrl,', to right:',ntrr
c     &       ,'recv from left:',rcvl,',recv from right:',rcvr
c        close(10)
        
        IF(NTRL>0) THEN
          CALL MPI_SEND(VERTEXL,9*NTRL,MTYPE,MYLEFT,MYRANK  ,MPI_COMM_EDDY,IERR)
          CALL MPI_SEND(UNVECTL,3*NTRL,MTYPE,MYLEFT,MYRANK+1,MPI_COMM_EDDY,IERR)
          CALL MPI_SEND(TRINOL ,NTRL  ,MTYPE,MYLEFT,MYRANK+2,MPI_COMM_EDDY,IERR)
        ENDIF

        IF(NTRR>0) THEN
          CALL MPI_SEND(VERTEXR,9*NTRR,MTYPE,MYRIGHT,MYRANK+10 ,MPI_COMM_EDDY,IERR)
          CALL MPI_SEND(UNVECTR,3*NTRR,MTYPE,MYRIGHT,MYRANK+11,MPI_COMM_EDDY,IERR)
          CALL MPI_SEND(TRINOR ,NTRR  ,MTYPE,MYRIGHT,MYRANK+12,MPI_COMM_EDDY,IERR)
        ENDIF

        IF(RCVL>0) THEN
          CALL MPI_RECV(VERTEXL,9*RCVL,MTYPE,MYLEFT,MYLEFT+10,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(UNVECTL,3*RCVL,MTYPE,MYLEFT,MYLEFT+11,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(TRINOL ,RCVL  ,MTYPE,MYLEFT,MYLEFT+12,MPI_COMM_EDDY,STATUS,IERR)
        ENDIF

        IF(RCVR>0) THEN
          CALL MPI_RECV(VERTEXR,9*RCVR,MTYPE,MYRIGHT,MYRIGHT  ,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(UNVECTR,3*RCVR,MTYPE,MYRIGHT,MYRIGHT+1,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(TRINOR ,RCVR  ,MTYPE,MYRIGHT,MYRIGHT+2,MPI_COMM_EDDY,STATUS,IERR)
        ENDIF

        
        !Eliminate common triangles
        K=0
        I=1
        DO WHILE(I<=RCVL)
          REP=0
          DO J=LB(IBD)+1,LB(IBD)+MB1(IBD)
            IF(TRINO(J)==TRINOL(I)) THEN
              REP=REP+1
              EXIT
            ENDIF
          ENDDO
          IF(REP==0) THEN
            K=K+1
            TRINOL(K)=TRINOL(I)
            VERTEXL(:,:,K)=VERTEXL(:,:,I)
            UNVECTL(:,K)=UNVECTL(:,I)
          ELSE
            TRINOL(I:RCVL-1)=TRINOL(I+1:RCVL)
            UNVECTL(:,I:RCVL-1)=UNVECTL(:,I+1:RCVL)
            VERTEXL(:,:,I:RCVL-1)=VERTEXL(:,:,I+1:RCVL)
            I=I-1
            RCVL=RCVL-1
          ENDIF
          I=I+1
        ENDDO
        RCVL=K

        K=0
        I=1
        DO WHILE(I<=RCVR)
          REP=0
          DO J=LB(IBD)+1,LB(IBD)+MB1(IBD)
            IF(TRINO(J)==TRINOR(I)) THEN
              REP=REP+1
              EXIT
            ENDIF
          ENDDO
          IF(REP==0) THEN
            K=K+1
            TRINOR(K)=TRINOR(I)
            UNVECTR(:,K)=UNVECTR(:,I)
            VERTEXR(:,:,K)=VERTEXR(:,:,I)
          ELSE
            TRINOR(I:RCVR-1)=TRINOR(I+1:RCVR)
            UNVECTR(:,I:RCVR-1)=UNVECTR(:,I+1:RCVR)
            VERTEXR(:,:,I:RCVR-1)=VERTEXR(:,:,I+1:RCVR)
            I=I-1
            RCVR=RCVR-1
          ENDIF
          I=I+1
        ENDDO
        RCVR=K

c        write(6,*) 'myrank=',myrank,'rcvl=',rcvl,',rcvr=',rcvr
c        call mpi_finalize(ierr)
c        stop

        IF(RCVL>0) THEN
          IF(mbd<nbd) THEN
            n=mb1(ibd)+rcvl
            IF(n>mb(ibd)) THEN
              s=n-mb(ibd)

              i1=lb(ibd+1)+1+s
              i2=lb(ibd+1)+mb(ibd+1)+s
              trino(i1:i2)=trino(lb(ibd+1)+1:lb(ibd+1)+mb(ibd+1))
              unvect(:,i1:i2)=unvect(:,lb(ibd+1)+1:lb(ibd+1)+mb(ibd+1))
              vertex(:,:,i1:i2)=vertex(:,:,lb(ibd+1)+1:lb(ibd+1)+mb(ibd+1))
              vertexc(:,i1:i2)=vertexc(:,lb(ibd+1)+1:lb(ibd+1)+mb(ibd+1))
              area(i1:i2)=area(lb(ibd+1)+1:lb(ibd+1)+mb(ibd+1))

              i1=lb(ibd)+mb1(ibd)+1
              i2=lb(ibd)+mb1(ibd)+rcvl
              trino(i1:i2)=trinol(1:rcvl)
              vertex(:,:,i1:i2)=vertexl(:,:,1:rcvl)
              unvect(:,i1:i2)=unvectl(:,1:rcvl)
              vertexc(:,i1:i2)=sum(vertex(:,:,i1:i2),2)/3.0
              do i=i1,i2
                area(i) = calc_triarea(vertex(:,:,i))
              enddo
            ELSEIF(n<mb(ibd)) THEN
              s=mb(ibd)-n

              i1=lb(ibd)+mb1(ibd)+1
              i2=lb(ibd)+mb1(ibd)+rcvl
              trino(i1:i2)=trinol(1:rcvl)
              unvect(:,i1:i2)=unvectl(:,1:rcvl)
              vertex(:,:,i1:i2)=vertexl(:,:,1:rcvl)
              vertexc(:,i1:i2)=sum(vertex(:,:,i1:i2),2)/3.0              
              do i=i1,i2
                area(i) = calc_triarea(vertex(:,:,i))
              enddo

              trino(i2+1:nfacet-s)=trino(i2+1+s:nfacet)
              unvect(:,i2+1:nfacet-s)=unvect(:,i2+1+s:nfacet)
              vertex(:,:,i2+1:nfacet-s)=vertex(:,:,i2+1+s:nfacet)
              vertexc(:,i2+1:nfacet-s)=vertexc(:,i2+1+s:nfacet)
              area(i2+1:nfacet-s)=area(i2+1+s:nfacet)
            ELSE
              i1=lb(ibd)+mb1(ibd)+1
              i2=lb(ibd)+mb1(ibd)+rcvl
              trino(i1:i2)=trinol(1:rcvl)
              unvect(:,i1:i2)=unvectl(:,1:rcvl)
              vertex(:,:,i1:i2)=vertexl(:,:,1:rcvl)
              vertexc(:,i1:i2)=sum(vertex(:,:,i1:i2),2)/3.0
              do i=i1,i2
                area(i) = calc_triarea(vertex(:,:,i))
              enddo
            ENDIF
          ELSE
            i1=lb(ibd)+mb1(ibd)+1
            i2=lb(ibd)+mb1(ibd)+rcvl
            trino(i1:i2)=trinol(1:rcvl)
            unvect(:,i1:i2)=unvectl(:,1:rcvl)
            vertex(:,:,i1:i2)=vertexl(:,:,1:rcvl)
            vertexc(:,i1:i2)=sum(vertex(:,:,i1:i2),2)/3.0
            do i=i1,i2
              area(i) = calc_triarea(vertex(:,:,i))
            enddo
          ENDIF

          mb1(ibd)=mb1(ibd)+rcvl
        ENDIF


        IF(RCVR>0) THEN
          IF(mbd<nbd) THEN
            n=mb1(ibd)+rcvr
            IF(n>mb(ibd)) THEN
              s=n-mb(ibd)

              i1=lb(ibd+1)+1+s
              i2=lb(ibd+1)+mb(ibd+1)+s
              trino(i1:i2)=trino(lb(ibd+1)+1:lb(ibd+1)+mb(ibd+1))
              unvect(:,i1:i2)=unvect(:,lb(ibd+1)+1:lb(ibd+1)+mb(ibd+1))
              vertex(:,:,i1:i2)=vertex(:,:,lb(ibd+1)+1:lb(ibd+1)+mb(ibd+1))
              vertexc(:,i1:i2)=vertexc(:,lb(ibd+1)+1:lb(ibd+1)+mb(ibd+1))
              area(i1:i2)=area(lb(ibd+1)+1:lb(ibd+1)+mb(ibd+1))

              i1=lb(ibd)+mb1(ibd)+1
              i2=lb(ibd)+mb1(ibd)+rcvr
              trino(i1:i2)=trinor(1:rcvr)
              unvect(:,i1:i2)=unvectr(:,1:rcvr)
              vertex(:,:,i1:i2)=vertexr(:,:,1:rcvr)
              vertexc(:,i1:i2)=sum(vertex(:,:,i1:i2),2)/3.0
              do i=i1,i2
                area(i) = calc_triarea(vertex(:,:,i))
              enddo
            ELSEIF(n<mb(ibd)) THEN
              s=mb(ibd)-n

              i1=lb(ibd)+mb1(ibd)+1
              i2=lb(ibd)+mb1(ibd)+rcvr
              trino(i1:i2)=trinor(1:rcvr)
              vertex(:,:,i1:i2)=vertexr(:,:,1:rcvr)
              unvect(:,i1:i2)=unvectr(:,1:rcvr)
              vertexc(:,i1:i2)=sum(vertex(:,:,i1:i2),2)/3.0
              do i=i1,i2
                area(i) = calc_triarea(vertex(:,:,i))
              enddo

              trino(i2+1:nfacet-s)=trino(i2+1+s:nfacet)
              unvect(:,i2+1:nfacet-s)=unvect(:,i2+1+s:nfacet)
              vertex(:,:,i2+1:nfacet-s)=vertex(:,:,i2+1+s:nfacet)
              vertexc(:,i2+1:nfacet-s)=vertexc(:,i2+1+s:nfacet)
              area(i2+1:nfacet-s)=area(i2+1+s:nfacet)

            ELSE
              i1=lb(ibd)+mb1(ibd)+1
              i2=lb(ibd)+mb1(ibd)+rcvr
              trino(i1:i2)=trinor(1:rcvr)
              vertex(:,:,i1:i2)=vertexr(:,:,1:rcvr)
              unvect(:,i1:i2)=unvectr(:,1:rcvr)
              vertexc(:,i1:i2)=sum(vertex(:,:,i1:i2),2)/3.0
              do i=i1,i2
                area(i) = calc_triarea(vertex(:,:,i))
              enddo
            ENDIF
          ELSE
            i1=lb(ibd)+mb1(ibd)+1
            i2=lb(ibd)+mb1(ibd)+rcvr
            trino(i1:i2)=trinor(1:rcvr)
            vertex(:,:,i1:i2)=vertexr(:,:,1:rcvr)
            unvect(:,i1:i2)=unvectr(:,1:rcvr)
            vertexc(:,i1:i2)=sum(vertex(:,:,i1:i2),2)/3.0
            do i=i1,i2
              area(i) = calc_triarea(vertex(:,:,i))
            enddo
          ENDIF

          mb1(ibd)=mb1(ibd)+rcvr
        ENDIF

        mb(ibd)=mb1(ibd)

        do i=ibd+1,nbd
          lb(i) = lb(i-1)+mb(i-1)
        enddo

      ENDDO
      
      do ibd=mbd+1,nbd
        lb(ibd) = lb(ibd-1)+mb(ibd-1)
      enddo
      nfacet = sum(mb(1:nbd))
 
      DEALLOCATE(mb1)

c      call mpi_finalize(ierr)
c      stop

      RETURN

      END
C---------------------------------------------------------------------------



C---- function calc_triarea----------------------N. Beratlis-01 Apr. 2010---
C
C     PURPOSE: Calculate area of triangle defined by 3 vertices
C
C---------------------------------------------------------------------------
      real function calc_triarea(vertex)

      implicit none
c
c.... Input/Output arrays
      real    vertex(3,3)
c
c.... Local arrays
      integer i
      real    a(3),b(3),c(3)

      calc_triarea=0.0

      do i=1,3
        a(i) = vertex(i,1)-vertex(i,2)
        b(i) = vertex(i,3)-vertex(i,2)
      enddo

      call cross(c,a,b)
      calc_triarea = 0.5*sqrt(c(1)**2.+c(2)**2.+c(3)**2.)      

      return
      end
C---------------------------------------------------------------------------


C---SUBROUTINE-------------------------------------------------------
C
C     PURPOSE: Reads an ASCII or BINARY stl file of the immersed object
C
      SUBROUTINE READSTL_MOD(UNVECT,VERTEX,VERTEXC,AREAF,TRINO,NFACET,ZC,NZ,NBD,MBD,IFLAG)
C
c
c outputs
c UNVECT: outward normal
c VERTEX: vertices of the triangles
c VERTEXC: centers of the triangles
c AREAF: areas of the triangles
c TRINO: face numbers
c
c inputs
c IFLAG: 0 (binary input file) or > 0 (formatted input file)
c
c The reading of the formatted files is modified
c To restart the unformatted version is used for the moving
c boundaries
c
C---------------------------------------------------------------------
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
      INTEGER NFACET,NBD,IFLAG,MBD,NZ
      REAL    DSMIN,DSMAX,DSAVE
      REAL    ZC(NZ)
      INTEGER TRINO(NFACET)
      REAL*4, DIMENSION(:),ALLOCATABLE :: DS

      REAL    UNVECT(3,NFACET),VERTEX(3,3,NFACET),VERTEXC(3,NFACET)
      REAL    AREAF(NFACET)
      CHARACTER FILENAME*80

c..local arrays
      INTEGER i,ibd,ifacet,ilb,ile,nfctmax,nfctmin
      INTEGER ifct,nfct(nbd),iflagl
      REAL    a(3),b(3),c(3)
      REAL    vecmag,dotprod
      REAL*4  UNVECT1(3,NFACET),VERTEX1(3,3,NFACET),VERTEX1C(3,NFACET)
      REAL*4  UN(3),VC(3)
      CHARACTER string*20

      nfct =0

      IF(MYSIZE==1 .OR. IMBDCMP==0) THEN

c	 write(6,*) myrank,'inside readstl'

	iflagl=iflag

	do ibd=1,nbd

	  ilb = lb(ibd)+1
	  ile = lb(ibd)+mb(ibd)

          if(myrank==0) then

	    filename = solid(ibd)
	    if(ibd>=mbd) then
	      if(ifield>0) then
		filename = 'res.mbd.'//index(ibd)
		iflagl=0
	      endif
	    endif

	    IF(iflagl>0) THEN

	      open(unit=20,file=trim(filename),form='formatted'
     &	       ,status='old')
	      READ(20,*)string
	      do ifacet=ilb,ile
		trino(ifacet) = ifacet
		READ(20,*) string,string,UNVECT1(1,IFACET),UNVECT1(2,IFACET),UNVECT1(3,IFACET)
		READ(20,*) string
		READ(20,*) string,VERTEX1(1,1,IFACET),VERTEX1(2,1,IFACET),VERTEX1(3,1,IFACET)
		READ(20,*) string,VERTEX1(1,2,IFACET),VERTEX1(2,2,IFACET),VERTEX1(3,2,IFACET)
		READ(20,*) string,VERTEX1(1,3,IFACET),VERTEX1(2,3,IFACET),VERTEX1(3,3,IFACET)
		DO I=1,3
		   VERTEX1C(I,IFACET)=
     %(VERTEX1(I,1,IFACET)+VERTEX1(I,2,IFACET)+VERTEX1(I,3,IFACET))
     %/3.
		ENDDO
		READ(20,*)string
		READ(20,*)string
	      enddo
	      close(20)

	      UNVECT = UNVECT1
	      VERTEX = VERTEX1
	      VERTEXC = VERTEX1C

	    ELSE

	      open(unit=20,file=trim(filename),form='UNFORMATTED'
     &		  ,status='OLD')
	      READ(20) IFACET
	      WRITE(6,*) IFACET,ILB,ILE
	      DO IFACET=ILB,ILE
		trino(ifacet) = ifacet
		READ(20) UNVECT1(1,IFACET),UNVECT1(2,IFACET),UNVECT1(3,IFACET)
		READ(20) VERTEX1C(1,IFACET),VERTEX1C(2,IFACET),VERTEX1C(3,IFACET)
		READ(20) VERTEX1(1,1,IFACET),VERTEX1(2,1,IFACET),VERTEX1(3,1,IFACET)
		READ(20) VERTEX1(1,2,IFACET),VERTEX1(2,2,IFACET),VERTEX1(3,2,IFACET)
		READ(20) VERTEX1(1,3,IFACET),VERTEX1(2,3,IFACET),VERTEX1(3,3,IFACET)
	      ENDDO
	      CLOSE(20)

	      UNVECT = UNVECT1
	      VERTEX = VERTEX1
	      VERTEXC = VERTEX1C

c	       write(6,*) 'object:',unvect(:,1307)

	    ENDIF

	  ENDIF
	ENDDO

        CALL MPI_BCAST(UNVECT,NFACET*3,MTYPE,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(VERTEX,NFACET*9,MTYPE,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(VERTEXC,NFACET*3,MTYPE,0,MPI_COMM_EDDY,IERR)
c	 CALL MPI_BCAST(AREAF,NFACET,MTYPE,0,MPI_COMM_EDDY,IERR)
        CALL MPI_BCAST(TRINO,NFACET,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)  !!!!!!

      ELSE

	DO ibd=1,nbd

	  ilb = lb(ibd)+1
	  ile = lb(ibd)+mb(ibd)

	  filename = solid(ibd)
	  if(ibd>=mbd) then
	    if(ifield>0) then
	      filename = 'res.mbd.'//index(ibd)
	    endif
	  endif

	  open(unit=20,file=trim(filename),form='UNFORMATTED'
     &	       ,status='OLD')
	  READ(20) IFACET

	  DO IFACET=ILB,ILE
	    READ(20) UN(1),UN(2),UN(3)
	    READ(20) VC(1),VC(2),VC(3)
	    IF(VC(3)>=ZC(1)-IMBOVLP .AND. VC(3)<ZC(NZ-1)+IMBOVLP) THEN
	      nfct(ibd) = nfct(ibd)+1
	      ifct = nfct(ibd)
	      IF(ifacet>nfacet) THEN
		write(6,*) 'ERROR: Increase nfacet'
		call mpi_finalize(ierr)
		stop
	      ENDIF
	      TRINO(ifct)=ifacet
	      UNVECT1(:,ifct) = UN(:)
	      VERTEX1C(:,ifct) = VC(:)
	      READ(20) VERTEX1(1,1,ifct),VERTEX1(2,1,ifct),VERTEX1(3,1,ifct)
	      READ(20) VERTEX1(1,2,ifct),VERTEX1(2,2,ifct),VERTEX1(3,2,ifct)
	      READ(20) VERTEX1(1,3,ifct),VERTEX1(2,3,ifct),VERTEX1(3,3,ifct)
	    ELSE
	      READ(20)
	      READ(20)
	      READ(20)
	    ENDIF
	  ENDDO

	ENDDO

	lb=0
	mb(1:nbd)=nfct(1:nbd)

	DO IBD=2,NBD
	  lb(ibd)=lb(ibd-1)+mb(ibd-1)
	ENDDO
	nfacet = sum(mb(1:nbd))

	unvect(:,1:nfacet)    = unvect1(:,1:nfacet)
	vertex(:,:,1:nfacet)  = vertex1(:,:,1:nfacet)
	vertexc(:,1:nfacet)   = vertex1c(:,1:nfacet)

	CALL MPI_REDUCE(NFACET,NFCTMAX,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_EDDY,IERR)
	CALL MPI_REDUCE(NFACET,NFCTMIN,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_EDDY,IERR)

	IF(MYRANK.eq.0) THEN
	  write(6,*) 'Max triangles in subdomain:',nfctmax,'min',nfctmin
	ENDIF

      ENDIF


      DO ibd=1,nbd

	ilb=lb(ibd)+1
	ile=lb(ibd)+mb(ibd)

	!Calculate area
	do ifacet=ilb,ile
	  do i=1,3
	    a(i) = vertex(i,1,ifacet)-vertex(i,2,ifacet)
	    b(i) = vertex(i,3,ifacet)-vertex(i,2,ifacet)
	  enddo

	  call cross(c,a,b)
	  areaf(ifacet) = 0.5*sqrt(c(1)**2.+c(2)**2.+c(3)**2.)
	enddo

	!Calculate normal
	IF(.false.) THEN
	  do ifacet=ilb,ile

	    do i=1,3
!!!!!!              a(i) = vertex(i,1,ifacet)-vertex(i,2,ifacet)
	      a(i) = vertex(i,2,ifacet)-vertex(i,1,ifacet)
	      b(i) = vertex(i,3,ifacet)-vertex(i,2,ifacet)
	    enddo

	    CALL cross(c,a,b)

!!!!!!            a(1) = vertexc(1,ifacet)
!!!!!!            a(2) = vertexc(2,ifacet)
!!!!!!            a(3) = vertexc(3,ifacet)
!!!!!!            if(dotprod(c,a)<0) c=-c
c	    c=a
!!!!!!            c(2)=0.0
	    c=c/vecmag(c,3)
	    unvect(1,ifacet)=c(1)
	    unvect(2,ifacet)=c(2)
	    unvect(3,ifacet)=c(3)
	  enddo
	ENDIF

      ENDDO

c      write(6,*) '2. object:',unvect(:,1307)


      RETURN

      END
