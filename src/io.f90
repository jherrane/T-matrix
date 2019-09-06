module io
use hdf5
use common
implicit none

contains 

subroutine read_mesh(mesh,file)
type (mesh_struct) :: mesh 
CHARACTER(LEN=80) :: file
CHARACTER(LEN=5), PARAMETER :: coord_dataset = "coord"     ! Dataset name
CHARACTER(LEN=6), PARAMETER :: etopol_dataset = "etopol"     ! Dataset2 name
CHARACTER(LEN=7), PARAMETER :: param_r_dataset = "param_r"     ! Dataset3 name
CHARACTER(LEN=7), PARAMETER :: param_i_dataset = "param_i"     ! Dataset4 name

integer(HID_T) :: file_id 
integer(HID_T) :: coord_dataset_id, etopol_dataset_id, param_r_dataset_id, param_i_dataset_id
integer(HID_T) :: coord_dataspace_id, etopol_dataspace_id, param_r_dataspace_id, param_i_dataspace_id
integer(HSIZE_T), dimension(2) :: dims_out, coord_dims, etopol_dims, param_r_dims, param_i_dims

integer :: error

double precision, dimension(:,:), allocatable :: coord
integer, dimension(:,:), allocatable :: etopol
double precision, dimension(:), allocatable :: param_r, param_i
double precision :: vol, a_eff


call h5open_f(error) ! initialize interface
call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error) ! open file

!_________________ Read coord__________________________________________________!

call h5dopen_f(file_id, coord_dataset, coord_dataset_id, error) ! open dataset
call h5dget_space_f(coord_dataset_id, coord_dataspace_id, error) 
call H5sget_simple_extent_dims_f(coord_dataspace_id, dims_out, coord_dims, error)

allocate(coord(coord_dims(1), coord_dims(2)))
call h5dread_f(coord_dataset_id, H5T_NATIVE_DOUBLE, coord, coord_dims, error)
call h5dclose_f(coord_dataset_id, error) ! close dataset

!___________________ Read etopol_______________________________________________!

call h5dopen_f(file_id, etopol_dataset, etopol_dataset_id, error) ! open dataset
call h5dget_space_f(etopol_dataset_id, etopol_dataspace_id, error) 
call H5sget_simple_extent_dims_f(etopol_dataspace_id, dims_out, etopol_dims, error)

allocate(etopol(etopol_dims(1), etopol_dims(2)))
call h5dread_f(etopol_dataset_id, H5T_NATIVE_INTEGER, etopol, etopol_dims, error)
call h5dclose_f(etopol_dataset_id, error)


!___________________ Read param_______________________________________________!

call h5dopen_f(file_id, param_r_dataset, param_r_dataset_id, error) ! open dataset
call h5dget_space_f(param_r_dataset_id, param_r_dataspace_id, error) 
call H5sget_simple_extent_dims_f(param_r_dataspace_id, dims_out, param_r_dims, error)

allocate(param_r(param_r_dims(1)))
call h5dread_f(param_r_dataset_id, H5T_NATIVE_DOUBLE, param_r, param_r_dims, error)
call h5dclose_f(param_r_dataset_id, error)

! Imaginary part
call h5dopen_f(file_id, param_i_dataset, param_i_dataset_id, error) ! open dataset
call h5dget_space_f(param_i_dataset_id, param_i_dataspace_id, error) 
call H5sget_simple_extent_dims_f(param_i_dataspace_id, dims_out, param_i_dims, error)

allocate(param_i(param_i_dims(1)))
call h5dread_f(param_i_dataset_id, H5T_NATIVE_DOUBLE, param_i, param_i_dims, error)
call h5dclose_f(param_i_dataset_id, error)

!______________________________________________________________________________!

call h5fclose_f(file_id, error) ! close file
call h5close_f(error) ! close inteface

allocate(mesh%coord(3,size(coord,2)))
mesh%coord = coord
print *,'   Number of nodes      =',size(coord,2)
mesh%N_node = size(coord,2)

allocate(mesh%etopol(4,size(etopol,2)))
mesh%etopol = etopol
print *,'   Number of elements   =',size(etopol,2)
mesh%N_tet = size(etopol,2)

allocate(mesh%param(size(param_r)))
mesh%param = dcmplx(param_r,param_i)

vol = get_tetra_vol(mesh)
a_eff = (3d0*vol/4d0/pi)**(1d0/3d0)
mesh%coord = mesh%coord*mesh%a/a_eff

end subroutine read_mesh

!******************************************************************************

subroutine T_write2file(T_mat, Nmax, fname)
double complex, intent(in) :: T_mat(:,:)
integer, intent(in) :: Nmax

CHARACTER(LEN=80), intent(in) :: fname
character(len=80) :: filename

CHARACTER(LEN=5), PARAMETER :: dsetname1 = "Taa_r" ! Dataset name
CHARACTER(LEN=5), PARAMETER :: dsetname2 = "Taa_i" ! Dataset name
CHARACTER(LEN=5), PARAMETER :: dsetname3 = "Tab_r" ! Dataset name
CHARACTER(LEN=5), PARAMETER :: dsetname4 = "Tab_i" ! Dataset name
CHARACTER(LEN=5), PARAMETER :: dsetname5 = "Tba_r" ! Dataset name
CHARACTER(LEN=5), PARAMETER :: dsetname6 = "Tba_i" ! Dataset name
CHARACTER(LEN=5), PARAMETER :: dsetname7 = "Tbb_r" ! Dataset name
CHARACTER(LEN=5), PARAMETER :: dsetname8 = "Tbb_i" ! Dataset name

INTEGER(HID_T) :: file_id       ! File identifier
INTEGER(HID_T) :: dset_id1       ! Dataset identifier
INTEGER(HID_T) :: dset_id2       ! Dataset identifier
INTEGER(HID_T) :: dset_id3       ! Dataset identifier
INTEGER(HID_T) :: dset_id4       ! Dataset identifier
INTEGER(HID_T) :: dset_id5       ! Dataset identifier
INTEGER(HID_T) :: dset_id6       ! Dataset identifier
INTEGER(HID_T) :: dset_id7       ! Dataset identifier
INTEGER(HID_T) :: dset_id8       ! Dataset identifier

INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

INTEGER(HSIZE_T), DIMENSION(2) :: dims  ! Dataset dimensions
INTEGER     ::    rank = 2                       ! Dataset rank

INTEGER     ::   error, nm ! Error flag
!INTEGER     :: i, j

filename = fname
nm = (Nmax+1)**2-1
dims = (/nm,nm/)
!Taa = T_mat(1:nm,1:nm)
!Tab = T_mat(1:nm, nm+1:2*nm)
!Tba = T_mat(nm+1:2*nm,1:nm)
!Tbb = T_mat(nm+1:2*nm,nm+1:2*nm)
     
CALL h5open_f(error)

! Create a new file using default properties.
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
! Create the dataspace.
CALL h5screate_simple_f(rank, dims, dspace_id, error)
! Create and write dataset using default properties.
CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(T_mat(1:nm,1:nm)), dims, error)
CALL h5dclose_f(dset_id1, error)

CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)
CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(T_mat(1:nm,1:nm)), dims, error)
CALL h5dclose_f(dset_id2, error)

!******************************************************************

CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id3, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)
CALL h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, real(T_mat(1:nm, nm+1:2*nm)), dims, error)
CALL h5dclose_f(dset_id3, error)

CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id4, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)
CALL h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, imag(T_mat(1:nm, nm+1:2*nm)), dims, error)
CALL h5dclose_f(dset_id4, error)

!******************************************************************

CALL h5dcreate_f(file_id, dsetname5, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id5, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)
CALL h5dwrite_f(dset_id5, H5T_NATIVE_DOUBLE, real(T_mat(nm+1:2*nm,1:nm)), dims, error)
CALL h5dclose_f(dset_id5, error)

CALL h5dcreate_f(file_id, dsetname6, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id6, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)
CALL h5dwrite_f(dset_id6, H5T_NATIVE_DOUBLE, imag(T_mat(nm+1:2*nm,1:nm)), dims, error)
CALL h5dclose_f(dset_id6, error)

!******************************************************************

CALL h5dcreate_f(file_id, dsetname7, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id7, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)
CALL h5dwrite_f(dset_id7, H5T_NATIVE_DOUBLE, real(T_mat(nm+1:2*nm,nm+1:2*nm)), dims, error)
CALL h5dclose_f(dset_id7, error)

CALL h5dcreate_f(file_id, dsetname8, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id8, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)
CALL h5dwrite_f(dset_id8, H5T_NATIVE_DOUBLE, imag(T_mat(nm+1:2*nm,nm+1:2*nm)), dims, error)
CALL h5dclose_f(dset_id8, error)

!******************************************************************
! Terminate access to the data space.
CALL h5sclose_f(dspace_id, error)
! Close the file.
CALL h5fclose_f(file_id, error)
! Close FORTRAN interface.
CALL h5close_f(error)

end subroutine T_write2file

end module 

