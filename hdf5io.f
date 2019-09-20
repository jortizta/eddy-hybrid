      SUBROUTINE HDF5_3DREAL(filename,uo,nx,ny,nz,idir)
c
      USE HDF5 ! This module contains all necessary modules
c
      IMPLICIT NONE
c
c... Input/Output arrays
      INTEGER NX,NY,NZ,IDIR
      REAL, DIMENSION(:,:,:)  :: UO(NX,NY,NZ)
      CHARACTER*(*) filename
c
c... Local arrays
      CHARACTER(LEN=5), PARAMETER :: dsetname = "var3d"     ! Dataset name
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

      INTEGER*4     ::   rank = 3                        ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(3) :: dims ! Dataset dimensions

      INTEGER*4     ::   error ! Error flag

      dims(1) = NX
      dims(2) = NY
      dims(3) = NZ

      ! Initialize FORTRAN interface.
      CALL h5open_f(error)

      ! Create a new file using default properties.
      if(idir==0) then
        CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
      else
        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      endif

      ! Create the dataspace.
      CALL h5screate_simple_f(rank, dims, dspace_id, error)

      ! Create the dataset with default properties.
      if(idir==0) then
        CALL h5dopen_f(file_id, dsetname, dset_id, error)
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, UO, dims, error)
      else
        CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE
     &     ,dspace_id, dset_id, error)
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, UO, dims, error)
      endif

      ! End access to the dataset and release resources used by it.
      CALL h5dclose_f(dset_id, error)

      ! Terminate access to the data space.
      CALL h5sclose_f(dspace_id, error)

      ! Close the file.
      CALL h5fclose_f(file_id, error)

      ! Close FORTRAN interface.
      CALL h5close_f(error)

      END SUBROUTINE




      SUBROUTINE HDF5_MPI_3DREAL(filename,uo,nx,ny,nz,idir)
c
      USE HDF5 ! This module contains all necessary modules

      include 'common.h'
      include 'mpif.h'
c
c... Input/Output arrays
      INTEGER :: NX,NY,NZ,IDIR
      REAL UO(NX,NY,NZ)
      CHARACTER*(*) filename
c
c... Local arrays
      CHARACTER(LEN=5), PARAMETER :: dsetname = "var3d"     ! Dataset name
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
      INTEGER(HID_T) :: filespace_id  ! Dataspace identifier
      INTEGER(HID_T) :: memspace_id   ! Dataspace identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier 

      INTEGER     ::   rank = 3                        ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(3) :: dimsf,dimsm,offset,count,block
     &     ,stride              ! Dataset dimensions

      INTEGER     ::   error,comm,info
      INTEGER     ::   k1 ! Error flag

      dimsf(1) = NX
      dimsf(2) = NY
      dimsf(3) = (NZ-2)*MYSIZE+2
      stride = 1

      info = MPI_INFO_NULL
      comm = mpi_comm_eddy
      
      ! Initialize FORTRAN interface
      CALL h5open_f(error) 
      IF(MYRANK==0) WRITE(*,*)filename
      ! Setup file access property list with parallel I/O access.
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)

      ! Create/Read the file collectively.
      IF(IDIR/=1) THEN
        CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error
     &     , access_prp = plist_id)
      ELSE
        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error
     &     , access_prp = plist_id)
      ENDIF
      CALL h5pclose_f(plist_id, error)

      !Create the data space for the file
      CALL h5screate_simple_f(rank, dimsf, filespace_id, error)
      offset(1) = 0
      offset(2) = 0 
      offset(3) = myrank*(nz-2)
      if(mysize>1 .AND. myrank>0 .AND. idir==1) then
        offset(3) = myrank*(nz-2)+1
      endif
      k1 = 2
      if(mysize==1 .OR. myrank==0) k1=1

      count = 1
      block(1) = NX
      block(2) = NY
      block(3) = NZ-2
      if(mysize==1) then
        block(3) = NZ
      else
        if(myrank==0 .OR. myrank==mysize-1) then
          block(3) = NZ-1
        endif
      endif
      if(idir/=1) block(3) = NZ

      CALL h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, offset
     &     , count, error, stride, block)

      !Create the data space for the memory
      CALL h5screate_simple_f(rank, block, memspace_id, error)

      ! Create the dataset with default properties.
      IF(idir/=1) THEN
        CALL h5dopen_f(file_id, dsetname, dset_id, error)
      ELSE
        CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE
     &     , filespace_id, dset_id, error)
      ENDIF

      ! Create property list for collective dataset write
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

      ! Write the dataset collectively. 
      IF(idir/=1) THEN
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, UO, block, error
     &        , memspace_id, filespace_id, xfer_prp = plist_id)
      ELSE
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, uo(:,:,k1), block
     &     , error, memspace_id, filespace_id, xfer_prp = plist_id)
      ENDIF

      ! Close resources.
      CALL h5sclose_f(memspace_id, error)
      CALL h5sclose_f(filespace_id, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5pclose_f(plist_id, error)
      CALL h5fclose_f(file_id, error)

      ! Close FORTRAN interface
      CALL h5close_f(error)

      END SUBROUTINE
C------------------------------------------------------------------------



C---- write_primevars_tmavg_hdf5 ------------ N. Beratlis-04 Jan 2011 ---
c
c     PURPOSE: Write time average variables to file
c
C------------------------------------------------------------------------
      subroutine write_primevars_tmavg_hdf5(hspdir,var,nxt,nyt,nzt,nvar,nsmpl,lim,nz)
c
      USE HDF5

      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer nxt,nyt,nzt,nz,nvar,nsmpl(nvar),lim(8)
      real*4  var(nxt,nyt,nzt,nvar)
      character*(*) hspdir
      character*150 filename
c
c.... Local arrays
      integer i,j,k,iv,error,info
      CHARACTER(LEN=5), PARAMETER :: dsetname = "var3d"     ! Dataset name
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
      INTEGER(HID_T) :: filespace_id  ! Dataspace identifier
      INTEGER(HID_T) :: memspace_id   ! Dataspace identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier 
      INTEGER(HID_T) :: group_id      ! Group identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute space identifier
      INTEGER(HID_T) :: attr_id     ! Attribute space identifier

      INTEGER     ::   rank = 3                        ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(3) :: dimsf,dimsm,offset,count,block
     &     ,stride              ! Dataset dimensions
      INTEGER     ::   arank = 1                        ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(1) :: adims ! Dataset dimensions
      

      do iv=1,nvar

        if(iv==1) then
          filename = trim(hspdir)//'u_tmavg.h5sp'
        elseif(iv==2) then
          filename = trim(hspdir)//'v_tmavg.h5sp'
        elseif(iv==3) then
          filename = trim(hspdir)//'w_tmavg.h5sp'
        elseif(iv==4) then
          filename = trim(hspdir)//'p_tmavg.h5sp'
        elseif(iv==5) then
          filename = trim(hspdir)//'uu_tmavg.h5sp'
        elseif(iv==6) then
          filename = trim(hspdir)//'vv_tmavg.h5sp'
        elseif(iv==7) then
          filename = trim(hspdir)//'ww_tmavg.h5sp'
        elseif(iv==8) then
          filename = trim(hspdir)//'pp_tmavg.h5sp'
        elseif(iv==9) then
          filename = trim(hspdir)//'uw_tmavg.h5sp'
        elseif(iv==10) then
          filename = trim(hspdir)//'uv_tmavg.h5sp'
        elseif(iv==11) then
          filename = trim(hspdir)//'vw_tmavg.h5sp'
        elseif(iv==12) then
          filename = trim(hspdir)//'uw1_tmavg.h5sp'
        endif
        
        dimsf(1) = nxt
        dimsf(2) = nyt
        dimsf(3) = lim(8)-lim(7)+1
        stride = 1

        info = MPI_INFO_NULL

        ! Initialize FORTRAN interface
        CALL h5open_f(error) 

        ! Setup file access property list with parallel I/O access.
        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        CALL h5pset_fapl_mpio_f(plist_id, mpi_comm_eddy, info, error)

        ! Create the file collectively.
        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error
     &      , H5P_DEFAULT_F, access_prp = plist_id)
        CALL h5pclose_f(plist_id, error)

        !Create the data space for the file        
        CALL h5screate_simple_f(rank, dimsf, filespace_id, error)

        ! Create the dataset with default properties.
        CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL
     &       , filespace_id, dset_id, error)
        CALL h5sclose_f(filespace_id, error)

        offset(1) = 0
        offset(2) = 0 
        offset(3) = (lim(5)+myrank*(nz-2))-lim(7)

        count = 1
        block(1) = nxt
        block(2) = nyt
        block(3) = nzt

        if(nzt==0) then
          offset = 0
          block = 1
        endif

        !Create the data space for the memory
        CALL h5screate_simple_f(rank, block, memspace_id, error)
        if(nzt==0) call h5Sselect_none_f(memspace_id, error)

        call h5dget_space_f(dset_id, filespace_id, error)

        CALL h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, offset
     &       , count, error, stride, block)
        if(nzt==0) call h5Sselect_none_f(filespace_id, error)

        ! Create property list for collective dataset write
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var(:,:,:,iv), block
     &        , error, memspace_id, filespace_id, xfer_prp = plist_id)

        ! Close resources.
        CALL h5sclose_f(memspace_id, error)
        CALL h5sclose_f(filespace_id, error)
        CALL h5pclose_f(plist_id, error)
        CALL h5dclose_f(dset_id, error)
          
        CALL h5fclose_f(file_id, error)

        ! Close FORTRAN interface
        CALL h5close_f(error) 

        !Root writes sample size as an attribute to "root (/)" group
        if(myrank.eq.0) then
          ! Initialize FORTRAN interface.
          CALL h5open_f(error)
          CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

         ! Create the dataspace for the attribute.
          adims(1) = 1
          CALL h5screate_simple_f(arank, adims, aspace_id, error)
          call h5gopen_f(file_id, '/', group_id, error)
          CALL h5acreate_f(group_id, 'sample', H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
          adims(1) = 1
          CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, nsmpl(iv), adims, error)
          CALL h5aclose_f(attr_id, error)
          call h5gclose_f(group_id, error)
          CALL h5sclose_f(aspace_id, error)
          ! Close the file.
          CALL h5fclose_f(file_id, error)

          ! Close FORTRAN interface.
          CALL h5close_f(error)

        endif

      enddo

      return

      end
C------------------------------------------------------------------------


C---- read_primevars_tmavg_hdf5 ------------ N. Beratlis-04 Jan 2011 ---
c
c     PURPOSE: Read time average variables to file
c
C------------------------------------------------------------------------
      subroutine read_primevars_tmavg_hdf5(hspdir,var,nxt,nyt,nzt,nvar,nsmpl,lim,nz)
c
      USE HDF5

      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer nxt,nyt,nzt,nz,nvar,nsmpl(nvar),lim(8)
      real*4  var(nxt,nyt,nzt,nvar)
      character*(*) hspdir
      character*150 filename
c
c.... Local arrays
      integer i,j,k,iv,error,info
      CHARACTER(LEN=5), PARAMETER :: dsetname = "var3d"     ! Dataset name
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
      INTEGER(HID_T) :: filespace_id  ! Dataspace identifier
      INTEGER(HID_T) :: memspace_id   ! Dataspace identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier 
      INTEGER(HID_T) :: group_id      ! Group identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute space identifier
      INTEGER(HID_T) :: attr_id     ! Attribute space identifier

      INTEGER     ::   rank = 3                        ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(3) :: dimsf,dimsm,offset,count,block
     &     ,stride              ! Dataset dimensions
      INTEGER     ::   arank = 1                        ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(1) :: adims ! Dataset dimensions

      do iv=1,nvar

        if(iv==1) then
          filename = trim(hspdir)//'u_tmavg.h5sp'
        elseif(iv==2) then
          filename = trim(hspdir)//'v_tmavg.h5sp'
        elseif(iv==3) then
          filename = trim(hspdir)//'w_tmavg.h5sp'
        elseif(iv==4) then
          filename = trim(hspdir)//'p_tmavg.h5sp'
        elseif(iv==5) then
          filename = trim(hspdir)//'uu_tmavg.h5sp'
        elseif(iv==6) then
          filename = trim(hspdir)//'vv_tmavg.h5sp'
        elseif(iv==7) then
          filename = trim(hspdir)//'ww_tmavg.h5sp'
        elseif(iv==8) then
          filename = trim(hspdir)//'pp_tmavg.h5sp'
        elseif(iv==9) then
          filename = trim(hspdir)//'uw_tmavg.h5sp'
        elseif(iv==10) then
          filename = trim(hspdir)//'uv_tmavg.h5sp'
        elseif(iv==11) then
          filename = trim(hspdir)//'vw_tmavg.h5sp'
        elseif(iv==12) then
          filename = trim(hspdir)//'uw1_tmavg.h5sp'
        endif
        
        dimsf(1) = nxt
        dimsf(2) = nyt
        dimsf(3) = lim(8)-lim(7)+1
        stride = 1

        info = MPI_INFO_NULL

        ! Initialize FORTRAN interface
        CALL h5open_f(error)

        ! Setup file access property list with parallel I/O access.
        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        CALL h5pset_fapl_mpio_f(plist_id, mpi_comm_self, info, error)

        ! Create the file collectively.
        CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error
     &     , access_prp = plist_id)
        CALL h5pclose_f(plist_id, error)

        !Create the data space for the file
        CALL h5screate_simple_f(rank, dimsf, filespace_id, error)

        ! Create the dataset with default properties.
c        CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL
c     &       , filespace_id, dset_id, error)
c        CALL h5sclose_f(filespace_id, error)

        offset(1) = 0
        offset(2) = 0 
        offset(3) = (lim(5)+myrank*(nz-2))-lim(7)

        count = 1
        block(1) = nxt
        block(2) = nyt
        block(3) = nzt

        if(nzt==0) then
          offset = 0
          block = 1
        endif

        !Create the data space for the memory
c        CALL h5screate_simple_f(rank, block, memspace_id, error)
c        if(nzt==0) call h5Sselect_none_f(memspace_id, error)
c        call h5dget_space_f(dset_id, filespace_id, error)

        CALL h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, offset
     &       , count, error, stride, block)      
c       call h5Sselect_none_f(filespace_id, error)

        !Create the data space for the memory
        CALL h5screate_simple_f(rank, block, memspace_id, error)
c       call h5Sselect_none_f(memspace_id, error)

        ! Create property list for collective dataset read
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! Create the dataset with default properties.
        CALL h5dopen_f(file_id, dsetname, dset_id, error)

        CALL h5dread_f(dset_id, H5T_NATIVE_REAL, var(:,:,:,iv), block, error
     &        , memspace_id, filespace_id, xfer_prp = plist_id)

        ! Close resources.
        CALL h5sclose_f(memspace_id, error)
        CALL h5sclose_f(filespace_id, error)
        CALL h5pclose_f(plist_id, error)
        CALL h5dclose_f(dset_id, error)

        CALL h5fclose_f(file_id, error)

        ! Close FORTRAN interface
        CALL h5close_f(error) 

        !Root writes sample size as an attribute to "root (/)" group
        if(myrank.eq.0) then
          ! Initialize FORTRAN interface.
          CALL h5open_f(error)
          CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

         ! Create the dataspace for the attribute.
          adims(1) = 1
          call h5gopen_f(file_id, '/', group_id, error)
          CALL h5aopen_f(group_id, "sample", attr_id, error)
          CALL h5aread_f(attr_id, H5T_NATIVE_INTEGER, nsmpl(iv), adims, error)
          CALL h5aclose_f(attr_id, error)
          call h5gclose_f(group_id, error)

          ! Close the file.
          CALL h5fclose_f(file_id, error)

          ! Close FORTRAN interface.
          CALL h5close_f(error)

        endif

      CALL MPI_BCAST(nsmpl(iv),1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)

      enddo

      return

      end
C------------------------------------------------------------------------


C---- read_primevar_tmavg_hdf5 -------------- N. Beratlis-04 Jan 2011 ---
c
c     PURPOSE: Read time average variables to file
c
C------------------------------------------------------------------------
      subroutine read_primevar_tmavg_hdf5(filename,var,nxt,nyt,nzt,nvar,nsmpl,lim,nz)
c
      USE HDF5

      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output variables
      integer nxt,nyt,nzt,nz,nvar,nsmpl,lim(8)
      real*4  var(nxt,nyt,nzt)
      character*(*) filename
c
c.... Local arrays
      integer i,j,k,iv,error,info
      CHARACTER(LEN=5), PARAMETER :: dsetname = "var3d"     ! Dataset name
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
      INTEGER(HID_T) :: filespace_id  ! Dataspace identifier
      INTEGER(HID_T) :: memspace_id   ! Dataspace identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier 
      INTEGER(HID_T) :: group_id      ! Group identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute space identifier
      INTEGER(HID_T) :: attr_id     ! Attribute space identifier

      INTEGER     ::   rank = 3                        ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(3) :: dimsf,dimsm,offset,count,block
     &     ,stride              ! Dataset dimensions
      INTEGER     ::   arank = 1                        ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(1) :: adims ! Dataset dimensions

        
      dimsf(1) = nxt
      dimsf(2) = nyt
      dimsf(3) = lim(8)-lim(7)+1
      stride = 1

      info = MPI_INFO_NULL

      ! Initialize FORTRAN interface
      CALL h5open_f(error)

      ! Setup file access property list with parallel I/O access.
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, mpi_comm_self, info, error)

      ! Create the file collectively.
      CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error
     &     , access_prp = plist_id)
      CALL h5pclose_f(plist_id, error)

      !Create the data space for the file
      CALL h5screate_simple_f(rank, dimsf, filespace_id, error)

        ! Create the dataset with default properties.
c        CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL
c     &       , filespace_id, dset_id, error)
c        CALL h5sclose_f(filespace_id, error)

      offset(1) = 0
      offset(2) = 0 
      offset(3) = (lim(5)+myrank*(nz-2))-lim(7)

      count = 1
      block(1) = nxt
      block(2) = nyt
      block(3) = nzt

      if(nzt==0) then
        offset = 0
        block = 1
      endif

        !Create the data space for the memory
c        CALL h5screate_simple_f(rank, block, memspace_id, error)
c        if(nzt==0) call h5Sselect_none_f(memspace_id, error)
c        call h5dget_space_f(dset_id, filespace_id, error)

      CALL h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, offset
     &       , count, error, stride, block)      
      if(nzt==0) call h5Sselect_none_f(filespace_id, error)

      !Create the data space for the memory
      CALL h5screate_simple_f(rank, block, memspace_id, error)
      if(nzt==0) call h5Sselect_none_f(memspace_id, error)

      ! Create property list for collective dataset read
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

      ! Create the dataset with default properties.
      CALL h5dopen_f(file_id, dsetname, dset_id, error)

      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, var, block, error
     &        , memspace_id, filespace_id, xfer_prp = plist_id)

      ! Close resources.
      CALL h5sclose_f(memspace_id, error)
      CALL h5sclose_f(filespace_id, error)
      CALL h5pclose_f(plist_id, error)
      CALL h5dclose_f(dset_id, error)

      CALL h5fclose_f(file_id, error)

      ! Close FORTRAN interface
      CALL h5close_f(error) 

      !Root writes sample size as an attribute to "root (/)" group
      if(myrank.eq.0) then
        ! Initialize FORTRAN interface.
        CALL h5open_f(error)
        CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

        ! Create the dataspace for the attribute.
        adims(1) = 1
        call h5gopen_f(file_id, '/', group_id, error)
        CALL h5aopen_f(group_id, "sample", attr_id, error)
        CALL h5aread_f(attr_id, H5T_NATIVE_INTEGER, nsmpl, adims, error)
        CALL h5aclose_f(attr_id, error)
        call h5gclose_f(group_id, error)

        ! Close the file.
        CALL h5fclose_f(file_id, error)

        ! Close FORTRAN interface.
        CALL h5close_f(error)
      endif

      CALL MPI_BCAST(nsmpl,1,MPI_INTEGER,0,MPI_COMM_EDDY,IERR)

      return

      end
C------------------------------------------------------------------------

C---- subroutine write_VPfield_hdf5_sp --- N. Beratlis - 28 Sep. 2010 ---
C
C     PURPOSE: Write VPfield file
C
C------------------------------------------------------------------------
      subroutine write_VPfield_hdf5_sp(hspdir,fileend,uo,vo,wo,p,nx,ny,nz,lim)
c
      USE HDF5

      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,lim(8)
      real    uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
      character*(*) fileend,hspdir
      character*300 filename
c
c.... Local arrays

      filename = trim(hspdir)//'U'//trim(fileend)//'.h5sp'
      call write_field_hdf5_sp(trim(filename),uo,nx,ny,nz,lim)

      filename = trim(hspdir)//'V'//trim(fileend)//'.h5sp'
      call write_field_hdf5_sp(trim(filename),vo,nx,ny,nz,lim)

      filename = trim(hspdir)//'W'//trim(fileend)//'.h5sp'
      call write_field_hdf5_sp(trim(filename),wo,nx,ny,nz,lim)

      filename = trim(hspdir)//'P'//trim(fileend)//'.h5sp'
      call write_field_hdf5_sp(trim(filename),p,nx,ny,nz,lim)

      return

      end
C------------------------------------------------------------------------


C---- subroutine write_field_hdf5_sp ----- N. Beratlis - 28 Sep. 2010 ---
C
C     PURPOSE: Write VPfield file
C
C------------------------------------------------------------------------
      subroutine write_field_hdf5_sp(filename,uo,nx,ny,nz,lim)
c
      USE HDF5

      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,lim(8)
      real    uo(nx,ny,nz)
      character*(*) filename
c
c.... Local arrays
      integer k,i1,i2,j1,j2,k1,k2,nxl,nyl,nzl,info,nzg
      CHARACTER*10,  dsetname     ! Dataset name
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
      INTEGER(HID_T) :: filespace_id  ! Dataspace identifier
      INTEGER(HID_T) :: memspace_id   ! Dataspace identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier 
      INTEGER(HID_T) :: group_id      ! Group identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute space identifier
      INTEGER(HID_T) :: attr_id     ! Attribute space identifier

      INTEGER     ::   rank = 3             ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(3) :: dimsf,offset,count,block
     &     ,stride              ! Dataset dimensions
      INTEGER     ::   arank = 1                        ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(1) :: adims ! Dataset dimensions
      INTEGER*4     ::   error ! Error flag
c
      nxl = lim(2)-lim(1)+1
      nyl = lim(4)-lim(3)+1
      nzl = lim(6)-lim(5)+1
      nzg = lim(8)-lim(7)+1

      dimsf(1) = nxl
      dimsf(2) = nyl
      dimsf(3) = nzg
      stride = 1

      info = MPI_INFO_NULL

      i1 = lim(1)
      i2 = lim(2)
      j1 = lim(3)
      j2 = lim(4)
      k1 = lim(5)
      k2 = lim(6)

      ! Initialize FORTRAN interface
      CALL h5open_f(error)
c
      ! Setup file access property list with parallel I/O access.
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, mpi_comm_eddy, info, error)
c
      ! Create the file collectively.
      CALL h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error
     &     , H5P_DEFAULT_F, access_prp = plist_id)
      CALL h5pclose_f(plist_id, error)

      !Create the data space for the file
      CALL h5screate_simple_f(rank, dimsf, filespace_id, error)

      dsetname = 'var3d'
      ! Create the dataset with default properties.
      CALL h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_REAL
     &     , filespace_id, dset_id, error)
      CALL h5sclose_f(filespace_id, error)

      offset(1) = 0
      offset(2) = 0
      offset(3) = (lim(5)+myrank*(nz-2))-lim(7)

      count = 1
      block(1) = nxl
      block(2) = nyl
      block(3) = nzl

      if(nzl<=0) then
        offset = 0
        block = 1
      endif

      !Create the data space for the memory
      CALL h5screate_simple_f(rank, block, memspace_id, error)
      if(nzl<=0) call h5Sselect_none_f(memspace_id, error)

      call h5dget_space_f(dset_id, filespace_id, error)

      CALL h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, offset
     &     , count, error, stride, block)
      if(nzl<=0) call h5Sselect_none_f(filespace_id, error)

      ! Create property list for collective dataset write
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, real(uo(i1:i2,j1:j2,k1:k2),4), block
     &     , error, memspace_id, filespace_id, xfer_prp = plist_id)

      ! Close resources.
      CALL h5sclose_f(memspace_id, error)
      CALL h5sclose_f(filespace_id, error)
      CALL h5pclose_f(plist_id, error)
      CALL h5dclose_f(dset_id, error)

      CALL h5fclose_f(file_id, error)

      ! Close FORTRAN interface
      CALL h5close_f(error)
c
      if(myrank.eq.0) then
        ! Initialize FORTRAN interface.
        CALL h5open_f(error)
        CALL h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error)
c
        ! Create the dataspace for the attribute.
        adims(1) = 1
        CALL h5screate_simple_f(arank, adims, dspace_id, error)
        CALL h5dcreate_f(file_id, 'NX', H5T_NATIVE_INTEGER
     &     ,dspace_id, dset_id, error)
        CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nxl, adims, error)
        CALL h5dclose_f(dset_id, error)
        CALL h5sclose_f(dspace_id, error)

        ! Create the dataspace for the attribute.
        adims(1) = 1
        CALL h5screate_simple_f(arank, adims, dspace_id, error)
        CALL h5dcreate_f(file_id, 'NY', H5T_NATIVE_INTEGER
     &     ,dspace_id, dset_id, error)
        CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nyl, adims, error)
        CALL h5dclose_f(dset_id, error)
        CALL h5sclose_f(dspace_id, error)

        ! Create the dataspace for the attribute.
        adims(1) = 1
        CALL h5screate_simple_f(arank, adims, dspace_id, error)
        CALL h5dcreate_f(file_id, 'NZ', H5T_NATIVE_INTEGER
     &     ,dspace_id, dset_id, error)
        CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nzg, adims, error)
        CALL h5dclose_f(dset_id, error)
        CALL h5sclose_f(dspace_id, error)

        CALL h5fclose_f(file_id, error)
        CALL h5close_f(error)
      endif

      return

      end
C------------------------------------------------------------------------

C---- subroutine hdf5_3dgrid_sp ------------ N. Beratlis-25 Feb. 2013
C---
C
      SUBROUTINE HDF5_3DGRID_SP(filename,X,Y,Z,nx,ny,nz)
C
C     PURPOSE: Write 3D grid to a HDF5 file as 3 1D arrays in single 
C     precision.
C------------------------------------------------------------------------
c
      USE HDF5 ! This module contains all necessary modules
c
      IMPLICIT NONE
c
c... Input/Output arrays
      INTEGER NX,NY,NZ
      REAL X(NX),Y(NY),Z(NZ)
      CHARACTER*(*) filename
c
c... Local arrays
      CHARACTER(LEN=1), PARAMETER :: xname = "X"     ! Dataset name
      CHARACTER(LEN=1), PARAMETER :: yname = "Y"     ! Dataset name
      CHARACTER(LEN=1), PARAMETER :: zname = "Z"     ! Dataset name
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

      INTEGER*4     ::   rank = 1                        ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(1) :: xdim,ydim,zdim ! Dataset dimensions

      INTEGER*4     ::   error ! Error flag
      REAL*4 X1(NX),Y1(NY),Z1(NZ)

      xdim(1) = NX
      ydim(1) = NY
      zdim(1) = NZ

      x1 = x
      y1 = y
      z1 = z

      ! Initialize FORTRAN interface.
      CALL h5open_f(error)

      ! Create a new file using default properties.
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      ! Create the dataspace for X grid.
      CALL h5screate_simple_f(rank, xdim, dspace_id, error)

      ! Create the dataset with default properties.
      CALL h5dcreate_f(file_id, xname, H5T_NATIVE_REAL
     &     ,dspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, x1, xdim, error)

      ! End access to the dataset and release resources used by it.
      CALL h5dclose_f(dset_id, error)

      ! Terminate access to the data space.
      CALL h5sclose_f(dspace_id, error)

      ! Create the dataspace for Y grid.
      CALL h5screate_simple_f(rank, ydim, dspace_id, error)

      ! Create the dataset with default properties.
      CALL h5dcreate_f(file_id, yname, H5T_NATIVE_REAL
     &     ,dspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, y1, ydim, error)

      ! End access to the dataset and release resources used by it.
      CALL h5dclose_f(dset_id, error)

      ! Terminate access to the data space.
      CALL h5sclose_f(dspace_id, error)

      ! Create the dataspace for Z grid.
      CALL h5screate_simple_f(rank, zdim, dspace_id, error)

      ! Create the dataset with default properties.
      CALL h5dcreate_f(file_id, zname, H5T_NATIVE_REAL
     &     ,dspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, z1, zdim, error)

      ! End access to the dataset and release resources used by it.
      CALL h5dclose_f(dset_id, error)

      ! Terminate access to the data space.
      CALL h5sclose_f(dspace_id, error)

      ! Close the file.
      CALL h5fclose_f(file_id, error)

      ! Close FORTRAN interface.
      CALL h5close_f(error)

      END SUBROUTINE

C------------------------------------------------------------------------

