      SUBROUTINE REFRESHBC(X,IM,KM)

c      IMPLICIT NONE
      include 'common.h'
      include 'mpif.h'
      INTEGER IM,KM,IP
      INTEGER STATUS(MPI_STATUS_SIZE)

      REAL X(IM,*),Y(IM)

      IF(MYRANK==0 .AND. ITYPE(5)/=500) Y(1:IM) = X(1:IM,1)
      
      CALL MPI_SENDRECV(X(1,KM-1),IM,MTYPE,MYRIGHT,10
     &     ,            X(1,1   ),IM,MTYPE,MYLEFT ,10
     &     ,            MPI_COMM_EDDY,STATUS,IERR)

      IF(MYRANK==0 .AND. ITYPE(5)/=500) X(1:IM,1) = Y(1:IM)

      IF(MYRANK==MYSIZE-1 .AND. ITYPE(6)/=500) Y(1:IM) = X(1:IM,KM)

      CALL MPI_SENDRECV(X(1,2   ),IM,MTYPE,MYLEFT ,11
     &     ,            X(1,KM  ),IM,MTYPE,MYRIGHT,11
     &     ,            MPI_COMM_EDDY,STATUS,IERR)

      IF(MYRANK==MYSIZE-1 .AND. ITYPE(6)/=500) X(1:IM,KM) = Y(1:IM)
c
      RETURN
      END
C---------------------------------------------------------------------

C---- subroutine refreshflag----------------N. Beratlis-17 Aug 2009---
C
      SUBROUTINE REFRESHFLAG(X,IM,JM,KM)
C
C     PURPOSE: Refresh ghost cells of flag array between processors.
C
C---------------------------------------------------------------------

c      IMPLICIT NONE
      include 'common.h'
      include 'mpif.h'
      include 'immersed.h'
      INTEGER IM,JM,KM
      INTEGER STATUS(MPI_STATUS_SIZE)
      
      INTEGER X(IM,JM,KM)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: Y
      
      if(idomy==0) then
        allocate(y(im,jm))   !!!!!!
        IF(MYRANK==0 .AND. ITYPE(5)/=500) Y(1:IM,1:JM) = X(1:IM,1:JM,1)   !!!!!!
        
        CALL MPI_SENDRECV(X(1,1,KM-1),IM*JM,MPI_INTEGER,MYRIGHT,10
     &       ,            X(1,1,1   ),IM*JM,MPI_INTEGER,MYLEFT ,10
     &       ,            MPI_COMM_EDDY,STATUS,IERR)

        IF(MYRANK==0 .AND. ITYPE(5)/=500) X(1:IM,1:JM,1) = Y(1:IM,1:JM)   !!!!!!
      
        IF(MYRANK==MYSIZE-1 .AND. ITYPE(6)/=500) Y(1:IM,1:JM) = X(1:IM,1:JM,KM)   !!!!!!

        CALL MPI_SENDRECV(X(1,1,2   ),IM*JM,MPI_INTEGER,MYLEFT ,11
     &       ,            X(1,1,KM  ),IM*JM,MPI_INTEGER,MYRIGHT,11
     &       ,            MPI_COMM_EDDY,STATUS,IERR)

        IF(MYRANK==MYSIZE-1 .AND. ITYPE(6)/=500) X(1:IM,1:JM,KM) = Y(1:IM,1:JM)   !!!!!!
        deallocate(y)   !!!!!!

      else
c        write(6,*) 'refreshflag: sending',myright,jm-1,',receiving',myleft,1
        call mpi_sendrecv(x(:,jm-1,:),im*km,mpi_integer,myright,10
     &        ,           x(:,1,:)   ,im*km,mpi_integer,myleft ,10
     &        ,           mpi_comm_eddy,status,ierr)

c        write(6,*) 'refreshflag: sending',myright,2,',receiving',jm
        call mpi_sendrecv(x(:,2,:   ),im*km,mpi_integer,myleft,11
     &        ,           x(:,jm,:  ),im*km,mpi_integer,myright ,11
     &        ,           mpi_comm_eddy,status,ierr)
      endif
c

      RETURN
      END
C---------------------------------------------------------------------


c---- subroutine yrefresh3darray-----------N. Beratlis-17 Aug. 2009---
C
      SUBROUTINE YREFRESH3DARRAY(ARRAY,NX,NY,NZ)
C
C     PURPOSE: Refresh 3D array(nx,ny,nz) in the y direction. Assume array
C     is periodic in y direction.
C
C---------------------------------------------------------------------
      implicit none

      integer nx,ny,nz
      integer array(nx,ny,nz)

      array(:,1 ,:) = array(:,ny-1,:)
      array(:,ny,:) = array(:, 2  ,:)

      return

      end
C---------------------------------------------------------------------

c---- subroutine yrefresh3darray-----------N. Beratlis-17 Dec. 2009---
C
      SUBROUTINE REFRESHCNTLN(ARRAY,NX,NY,NZ)
C
C     PURPOSE: Refresh 3D array(nx,ny,nz) in the centerline. Assume array
C     is periodic in y direction.
C
C---------------------------------------------------------------------
c      implicit none
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
      integer nx,ny,nz
      integer array(nx,ny,nz)
c
      integer j
      integer status(mpi_status_size)
c
      if(idomy==0) then
        do j=1,ny
          array(1,j,:) = array(2,jsym(j),:)
        enddo
      else
         call mpi_sendrecv(array(2,:,:),ny*nz,mpi_integer,ranksym,1
     &          ,          array(1,:,:),ny*nz,mpi_integer,ranksym,1
     &          ,           mpi_comm_eddy,status,ierr)

      endif

      return

      end
C---------------------------------------------------------------------
