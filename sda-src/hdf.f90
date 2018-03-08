module hdf  
  use HDF5
  
  integer(hid_t)                    ::  ihdf5_in                !File identifier for HDF5 input file
  integer(hid_t)                    ::  ihdf5_out               !File identifier for HDF5 output file

  integer(hid_t)                    ::  sdab_id                 ! HDF handler for compound type
  integer(hid_t)                    ::  h5t_int                 ! HDF handler for interger type
  integer(hsize_t)                  ::  g_dims(2)               ! grid dimemsion
  integer(hsize_t)                  ::  c_dims(2)               ! cb dimemsion
contains
  subroutine saveHDFreal(file_id, arr, dsetname, dims)
    integer(hid_t)    :: file_id
    integer(hid_t)    :: dtype
    real              :: arr(*)
    integer(hsize_t)  :: dims(:)
    character(len=*)  :: dsetname
    integer           :: ierr
    integer(hid_t)    :: dspace_id
    integer(hid_t)    :: dset_id

    CALL h5tcopy_f(H5T_NATIVE_REAL, dtype, ierr)
    ! Create the dataspace.
    CALL h5screate_simple_f(size(dims), dims, dspace_id, ierr)
    ! Create the dataset with default properties.
    CALL h5dcreate_f(file_id, trim(dsetname), dtype, dspace_id, dset_id, ierr)
    ! Write the dataset.
    CALL h5dwrite_f(dset_id, dtype, arr, dims, ierr)
    CALL h5dclose_f(dset_id, ierr)
    CALL h5sclose_f(dspace_id, ierr)
    if (ierr /= 0) call USTOP('error when saving '//trim(dsetname))
  end subroutine

  subroutine saveHDFdbl(file_id, arr, dsetname, dims)
    integer(hid_t)    :: file_id
    integer(hid_t)    :: dtype
    doubleprecision   :: arr(*)
    integer(hsize_t)  :: dims(:)
    character(len=*)  :: dsetname
    integer           :: ierr
    integer(hid_t)    :: dspace_id
    integer(hid_t)    :: dset_id

    CALL h5tcopy_f(H5T_NATIVE_DOUBLE, dtype, ierr)
    ! Create the dataspace.
    CALL h5screate_simple_f(size(dims), dims, dspace_id, ierr)
    ! Create the dataset with default properties.
    CALL h5dcreate_f(file_id, trim(dsetname), dtype, dspace_id, dset_id, ierr)
    ! Write the dataset.
    CALL h5dwrite_f(dset_id, dtype, arr, dims, ierr)
    CALL h5dclose_f(dset_id, ierr)
    CALL h5sclose_f(dspace_id, ierr)
    if (ierr /= 0) call USTOP('error when saving '//trim(dsetname))
  end subroutine
  
  subroutine saveHDFstring(file_id, arr, dsetname, dims, str_size)
    integer(hid_t)    :: file_id
    integer(hsize_t)  :: str_size
    character(len=str_size)      :: arr(*)
    integer(hsize_t)  :: dims(:)
    character(len=*)  :: dsetname
    integer           :: ierr
    integer(hid_t)    :: dspace_id
    integer(hid_t)    :: dset_id
    integer(hid_t)    :: dtype

    ! Create string data type
    CALL h5tcopy_f(H5T_FORTRAN_S1, dtype, ierr)
    CALL h5tset_size_f(dtype, str_size, ierr)

    ! Create the dataspace.
    CALL h5screate_simple_f(size(dims), dims, dspace_id, ierr)
    ! Create the dataset with default properties.
    CALL h5dcreate_f(file_id, trim(dsetname), dtype, dspace_id, dset_id, ierr)
    ! Write the dataset.
    CALL h5dwrite_f(dset_id, dtype, arr, dims, ierr)
    CALL h5dclose_f(dset_id, ierr)
    CALL h5sclose_f(dspace_id, ierr)
    if (ierr /= 0) call USTOP('error when saving '//trim(dsetname))

  end subroutine

  subroutine saveHDFint(file_id, arr, dsetname, dims, dtype)
    integer(hid_t)    :: file_id
    integer(hid_t)    :: dtype
    integer           :: arr(*)
    integer(hsize_t)  :: dims(:)
    character(len=*)  :: dsetname
    integer           :: ierr
    integer(hid_t)    :: dspace_id
    integer(hid_t)    :: dset_id

    ! Create the dataspace.
    CALL h5screate_simple_f(size(dims), dims, dspace_id, ierr)
    ! Create the dataset with default properties.
    CALL h5dcreate_f(file_id, trim(dsetname), dtype, dspace_id, dset_id, ierr)
    ! Write the dataset.
    CALL h5dwrite_f(dset_id, dtype, arr, dims, ierr)
    CALL h5dclose_f(dset_id, ierr)
    CALL h5sclose_f(dspace_id, ierr)
    if (ierr /= 0) call USTOP('error when saving '//trim(dsetname))
  end subroutine

  subroutine readHDFint(file_id, arr, dsetname, dims)
    integer(hid_t)    :: file_id
    integer           :: arr(*)
    character(len=*)  :: dsetname
    integer(hsize_t)  :: dims(:)
    !local
    integer(hid_t)    :: dtype
    integer           :: ierr
    integer(hid_t)    :: dset_id


    CALL h5dopen_f(file_id, trim(dsetname), dset_id, ierr)
    CALL h5dget_type_f(dset_id, dtype, ierr)
    CALL h5dread_f(dset_id, dtype, arr, dims, ierr)
    CALL h5dclose_f(dset_id, ierr)
    if (ierr /= 0) call USTOP('error when reading '//trim(dsetname))

  end subroutine


  subroutine readHDFreal(file_id, arr, dsetname, dims)
    integer(hid_t)    :: file_id
    real              :: arr(*)
    character(len=*)  :: dsetname
    integer(hsize_t)  :: dims(:)
    ! local
    integer(hid_t)    :: dtype
    integer           :: ierr
    integer(hid_t)    :: dset_id


    CALL h5dopen_f(file_id, trim(dsetname), dset_id, ierr)
    CALL h5dget_type_f(dset_id, dtype, ierr)
    CALL h5dread_f(dset_id, dtype, arr, dims, ierr)
    CALL h5dclose_f(dset_id, ierr)
    if (ierr /= 0) call USTOP('error when reading '//trim(dsetname))

  end subroutine



  subroutine readHDFdbl(file_id, arr, dsetname, dims)
    integer(hid_t)    :: file_id
    doubleprecision   :: arr(*)
    character(len=*)  :: dsetname
    integer(hsize_t)  :: dims(:)
    ! local
    integer(hid_t)    :: dtype
    integer           :: ierr
    integer(hid_t)    :: dset_id


    CALL h5dopen_f(file_id, trim(dsetname), dset_id, ierr)
    CALL h5dget_type_f(dset_id, dtype, ierr)
    CALL h5dread_f(dset_id, dtype, arr, dims, ierr)
    CALL h5dclose_f(dset_id, ierr)
    if (ierr /= 0) call USTOP('error when reading '//trim(dsetname))

  end subroutine



  subroutine readHDFstring(file_id, arr, dsetname, dims, str_size)
    integer(hid_t)    :: file_id
    integer(hsize_t)  :: str_size
    character(*)      :: arr(*)
    integer(hsize_t)  :: dims(:)
    character(len=*)  :: dsetname
    integer           :: ierr
    integer(hid_t)    :: dset_id
    integer(hid_t)    :: dtype


    CALL h5dopen_f(file_id, trim(dsetname), dset_id, ierr)
    CALL h5dget_type_f(dset_id, dtype, ierr)
    CALL h5dread_f(dset_id, dtype, arr, dims, ierr)
    CALL h5dclose_f(dset_id, ierr)
    if (ierr /= 0) call USTOP('error when reading '//trim(dsetname))
  end subroutine

end module
