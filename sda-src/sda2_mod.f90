module SDA
  !use GLOBAL, only                  :   DEFAULT_KIND
  use HDF5
  use ISO_C_BINDING
  logical                           ::  debug = .FALSE.
  integer                           ::  idebug                  !output file number for scenario runs (2D head differebces)
  integer                           ::  IGRID = 1               !fake sub grid number

  doubleprecision                   ::  DTRead                  !Total time used for reading files
  doubleprecision                   ::  DTSolver                !Total time used to solve equation
  doubleprecision                   ::  DTStpIni                !Total time used to solve equation
  doubleprecision                   ::  DTBD                    !Total time used to finalize the scenario run
  doubleprecision                   ::  DTStpEnd                !Total time used for step calculations


  integer                           ::  NSEN = 0                !total number of scenario runs, > 0 no baseline run; < 0 include baseline run
  integer                           ::  QWEL                    !flag if the well package is activated in the scenario runs, 1 is activated, other is none
  integer                           ::  QRCH                    !flag if the recharge package is activated in the scenario runs, 1 is activated, other is none
  integer                           ::  QHED                    !flag for writting head difference, if QHED > 0, write as formatted, QHED < 0, write as unformated, QHED = 0, no write
  integer                           ::  QCCB                    !flag for writting cell by cell budgets


  integer                           ::  mCell                   !number of active cells


  integer                           ::  NCBMAX                  !maximum number of boundary
  integer                           ::  NCB                     !number of boundary each time step

  integer                           ::  iSTEP                   !current time step
  integer                           ::  NSTEP                   !total number of time steps


  integer                           ::  iTypeBC
  integer                           ::  nTypeBC



  real,allocatable    ::  SDACR(:,:,:)            !baseline CR
  real,allocatable    ::  SDACC(:,:,:)            !baseline CC
  real,allocatable    ::  SDACV(:,:,:)            !baseline CV
  real,allocatable       ::  SDAHCOF(:,:,:)          !baseline HCOF
  !real,allocatable,dimension(:,:,:) ::  SDARHS                 !baseline RHS
  real,allocatable,dimension(:,:,:) ::  SC                      !baseline Storage Capacity
  real,allocatable,dimension(:,:,:) ::  SDASC                   !Storage Capacity
  integer,allocatable,dimension(:,:,:) ::  SDAIBD               !baseline IBOUND



  doubleprecision,allocatable    ::  HCHG1(:,:,:)             !new head change
  doubleprecision,allocatable    ::  HCHG0(:,:,:)            !previous head change

  integer*1                         ::  CBINDX(100)             !the order of CB in outputs
  integer                           ::  CBNTOT(100)             !total number of CB during the simulation
  character(len=16)                 ::  CBTYPE(100)             !CB name

  doubleprecision                   ::  CBVBVL(4,100)           !In&Out, Rate and Cumulative


  integer                           ::  NBZONEMax               !maximum number of boundary ZONEs for all scenarios
  integer                           ::  NBZONE                  !number of boundary ZONEs of a scenario run
  integer,allocatable               ::  BZONE(:,:,:)            !boundary ZONE, defined at NewScen
  doubleprecision,allocatable       ::  RBZONE(:,:,:)           !flow rate of each zone, Calculated each NewStp

  ! prefix used to construct dataset name in HDF5, full name will be 'prefix + 0000 + 000'
  character*2, parameter               ::  pxIB = 'IB'              !hdf5 name prefix for iBOUND
  character*2, parameter               ::  pxCR = 'CR'                  !hdf5 name prefix for baseline CR
  character*2, parameter               ::  pxCC = 'CC'                 !hdf5 name prefix for baseline CC
  character*2, parameter               ::  pxCV = 'CV'                !hdf5 name prefix for baseline CV
  character*2, parameter               ::  pxSC = 'SC'               !hdf5 name prefix for baseline SC
  character*2, parameter               ::  pxCB = 'CB'               !hdf5 name prefix for baseline boundary condutance
  character*2, parameter               ::  pxHC = 'HC'               !hdf5 name prefix for baseline boundary condutance


  integer, parameter                   ::  nMaxStep = 9999         !maximum number of boundary
  integer,allocatable, dimension(:,:)  ::  newCR, newCC, newCV, newHC, newSC, newIB  ! ilayer, istep
  integer,allocatable, dimension(:)    ::  newCB ! istep

  integer,allocatable, dimension(:)    ::  actLay, actRow, actCol ! active layer row column
  integer                              ::  nActive

  integer                           ::  iHED =0                 !output file number for scenario runs (2D head differebces)
  integer(hid_t)                    ::  ihdf5_in                !File identifier for HDF5 input file
  integer(hid_t)                    ::  ihdf5_out               !File identifier for HDF5 output file

  integer(hid_t)                    ::  sdab_id                 ! HDF handler for compound type
  integer(hid_t)                    ::  h5t_int                 ! HDF handler for interger type
  integer(hsize_t)                  ::  g_dims(2)               ! grid dimemsion
  integer(hsize_t)                  ::  c_dims(2)               ! cb dimemsion


  character*80                      ::  ScenName                !Scenario Name
  integer                           ::  iSen                    !index of scenario run

  integer                           ::  iCBB =  0               !output file number for scenario runs

  character*10                      ::  FORM='BINARY'
  character*15                      ::  ACCESS='SEQUENTIAL'
  character*10                      ::  ACTION(2)=(/'READ      ','READWRITE '/)

  !cell of head-dependent flow boundary condition
  integer,      allocatable                           ::  cbCOL(:),cbCOL_o(:)                  !col
  integer,      allocatable                           ::  cbROW(:),cbROW_o(:)                  !row
  integer,      allocatable                           ::  cbLAY(:),cbLAY_o(:)                  !layer
  integer,      allocatable                           ::  cbTyp(:),cbTyp_o(:)                   !boundary text

  real,    allocatable                                ::  cbCond(:),cbCond_o(:)                    !conductance
  real,    allocatable                                ::  cbHead(:),cbHead_o(:)                    !boundary head

!  integer,                                            ::  nCBType
!  integer,      allocatable                           ::  cbiType(:)
!  character*16, allocatable                           ::  cbTypeList(:)

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
    character(len=str_size), target      :: arr(*)
    integer(hsize_t)  :: dims(:)
    character(len=*)  :: dsetname
    integer           :: ierr
    integer(hid_t)    :: dspace_id
    integer(hid_t)    :: dset_id
    integer(hid_t)    :: dtype
    TYPE(C_PTR)       :: f_ptr

    ! Create string data type
    CALL h5tcopy_f(H5T_FORTRAN_S1, dtype, ierr)
    CALL h5tset_size_f(dtype, str_size, ierr)

    ! Create the dataspace.
    CALL h5screate_simple_f(size(dims), dims, dspace_id, ierr)
    ! Create the dataset with default properties.
    CALL h5dcreate_f(file_id, trim(dsetname), dtype, dspace_id, dset_id, ierr)
    ! Write the dataset.
    f_ptr = C_LOC(arr(1)(1:1))
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
    if (ierr /= 0) call USTOP('error when saving '//trim(dsetname))
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
  end subroutine

end module
