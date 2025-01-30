module hdf5_lib
  use hdf5
  use mpi_lib
  implicit none

  integer(HID_T), private :: file_id, dspace_id, dset_id, group_id, plist_id  ! Identifiers
  integer(HID_T), private :: dataspace, memspace
  integer(HSIZE_T), dimension(3), private :: datasize, offset
  integer, private :: error

  interface write_hdf5
    module procedure :: write_hdf5_dataset_in_root, write_hdf5_dataset_in_group
  end interface write_hdf5

contains

  !============================================================================!
  subroutine create_hdf5_file(filename)
    implicit none
    character(len=*), intent(in) :: filename

    ! Initialize FORTRAN interface.
    call h5open_f(error)
    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    ! Create a new file using default properties.
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)
    ! Terminate access to the file.
    call h5fclose_f(file_id, error)
    ! Close FORTRAN interface.
    call h5close_f(error)

  end subroutine create_hdf5_file
  !============================================================================!
  subroutine create_hdf5_group(filename, groupname)
    implicit none
    character(len=*), intent(in) :: filename, groupname

    ! Initialize FORTRAN interface.
    call h5open_f(error)
    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    ! Create a new file using default properties.
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)
    ! Create all group (folder).
    call h5gcreate_f(file_id, groupname, group_id, error)
    ! Close group
    call h5gclose_f(group_id,error)
    ! Close file.
    call h5fclose_f(file_id, error)
    ! Close FORTRAN interface.
    call h5close_f(error)

  end subroutine create_hdf5_group
  !============================================================================!
  subroutine write_hdf5_dataset_in_group(filename, groupname, name, data)
    implicit none
    character(len=*), intent(in) :: filename, groupname, name
    real(kind=8), dimension(:,:,:) :: data
    integer, parameter :: ndim = 3
    integer(8), dimension(ndim) :: dim_local
    integer(8), dimension(ndim) :: dim_global

    dim_local = [size(data, 1), size(data, 2), size(data, 3)]
    dim_global = dim_local*[nproc,1,1]
    datasize = dim_local
    offset = dim_local *[rank,0,0]
    ! Initialize FORTRAN interface.
    call h5open_f(error)
    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    ! Create a new file using default properties.
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)
    ! Create all group (folder).
    call h5gopen_f(file_id, groupname, group_id, error)
    !==========================================================================!
    ! write in filename / groupname
    !==========================================================================!
    ! Create the data space for the dataset.
    call h5screate_simple_f(ndim, dim_global, dataspace, error)
    ! Create the dataset with default properties.
    call h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, dataspace, dset_id, error)
    call h5screate_simple_f(ndim, datasize, memspace, error)
    ! Select hyperslab in the file.
    call h5sselect_hyperslab_f (dataspace, H5S_SELECT_SET_F, offset, datasize, error)
    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, datasize, error, &
    file_space_id = dataspace, mem_space_id = memspace, xfer_prp = plist_id)
    ! Close dataspaces.
    CALL h5sclose_f(dataspace, error)
    CALL h5sclose_f(memspace, error)
    ! Close the dataset and property list.
    CALL h5dclose_f(dset_id, error)
    CALL h5pclose_f(plist_id, error)
    !==========================================================================!
    ! Close group
    call h5gclose_f(group_id,error)
    ! Close file.
    call h5fclose_f(file_id, error)
    ! Close FORTRAN interface.
    call h5close_f(error)

  end subroutine write_hdf5_dataset_in_group
  !============================================================================!
  subroutine write_hdf5_dataset_in_root(filename, name, data)
    implicit none
    character(len=*), intent(in) :: filename, name
    real(kind=8), dimension(:,:,:) :: data
    integer, parameter :: ndim = 3
    integer(8), dimension(ndim) :: dim_local
    integer(8), dimension(ndim) :: dim_global

    dim_local = [size(data, 1), size(data, 2), size(data, 3)]
    dim_global = dim_local*[nproc,1,1]
    datasize = dim_local
    offset = dim_local *[rank,0,0]
    ! Initialize FORTRAN interface.
    call h5open_f(error)
    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    ! Create a new file using default properties.
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)
    !==========================================================================!
    ! write in filename /
    !==========================================================================!
    ! Create the data space for the dataset.
    call h5screate_simple_f(ndim, dim_global, dataspace, error)
    ! Create the dataset with default properties.
    call h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace, dset_id, error)
    call h5screate_simple_f(ndim, datasize, memspace, error)
    ! Select hyperslab in the file.
    call h5sselect_hyperslab_f (dataspace, H5S_SELECT_SET_F, offset, datasize, error)
    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, datasize, error, &
    file_space_id = dataspace, mem_space_id = memspace, xfer_prp = plist_id)
    ! Close dataspaces.
    CALL h5sclose_f(dataspace, error)
    CALL h5sclose_f(memspace, error)
    ! Close the dataset and property list.
    CALL h5dclose_f(dset_id, error)
    CALL h5pclose_f(plist_id, error)
    !==========================================================================!
    ! Close file.
    call h5fclose_f(file_id, error)
    ! Close FORTRAN interface.
    call h5close_f(error)

  end subroutine write_hdf5_dataset_in_root
  !============================================================================!

end module hdf5_lib
