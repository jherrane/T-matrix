module io
use hdf5
use common
implicit none

contains 

subroutine write2file(A,fname)
double complex, intent(in) :: A(:,:)
CHARACTER(LEN=38), intent(in) :: fname
CHARACTER(LEN=38) :: filename
!CHARACTER(LEN=7), PARAMETER :: filename = "A.h5" ! File name
CHARACTER(LEN=3), PARAMETER :: dsetname1 = "A_r" ! Dataset name
CHARACTER(LEN=3), PARAMETER :: dsetname2 = "A_i" ! Dataset name

INTEGER(HID_T) :: file_id       ! File identifier
INTEGER(HID_T) :: dset_id1       ! Dataset identifier
INTEGER(HID_T) :: dset_id2       ! Dataset identifier
INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


INTEGER(HSIZE_T), DIMENSION(2) :: dims  ! Dataset dimensions
INTEGER     ::    rank = 2                       ! Dataset rank

INTEGER     ::   error ! Error flag
!INTEGER     :: i, j

filename = fname
dims = (/size(A,1),size(A,2)/)
     
CALL h5open_f(error)

!
! Create a new file using default properties.
!
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

!
! Create the dataspace.
!
CALL h5screate_simple_f(rank, dims, dspace_id, error)

!
! Create and write dataset using default properties.
!
CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(A), dims, error)

CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(A), dims, error)

!
! End access to the dataset and release resources used by it.
!
CALL h5dclose_f(dset_id1, error)
CALL h5dclose_f(dset_id2, error)
!
! Terminate access to the data space.
!
CALL h5sclose_f(dspace_id, error)

!
! Close the file.
!
CALL h5fclose_f(file_id, error)

!
! Close FORTRAN interface.
!
CALL h5close_f(error)

end subroutine write2file

!_______________________________________________________________________

subroutine real_write2file(A,fname)
double precision, intent(in) :: A(:,:)
CHARACTER(LEN=38), intent(in) :: fname
CHARACTER(LEN=38) :: filename 

CHARACTER(LEN=7), PARAMETER :: dsetname1 = "mueller" ! Dataset name


INTEGER(HID_T) :: file_id       ! File identifier
INTEGER(HID_T) :: dset_id1       ! Dataset identifier
INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


INTEGER(HSIZE_T), DIMENSION(2) :: dims  ! Dataset dimensions
INTEGER     ::    rank = 2                       ! Dataset rank

INTEGER     ::   error ! Error flag
!INTEGER     :: i, j

dims = (/size(A,1),size(A,2)/)
filename = fname  
 
CALL h5open_f(error)

!
! Create a new file using default properties.
!
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

!
! Create the dataspace.
!
CALL h5screate_simple_f(rank, dims, dspace_id, error)

!
! Create and write dataset using default properties.
!
CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, A, dims, error)


!
! End access to the dataset and release resources used by it.
!
CALL h5dclose_f(dset_id1, error)

!
! Terminate access to the data space.
!
CALL h5sclose_f(dspace_id, error)

!
! Close the file.
!
CALL h5fclose_f(file_id, error)

!
! Close FORTRAN interface.
!
CALL h5close_f(error)

end subroutine real_write2file
!_---______________________________________________________________

subroutine cmplx_vec_write2file(vec,fname)
double complex, intent(in) :: vec(:)
CHARACTER(LEN=28), intent(in) :: fname
CHARACTER(LEN=28) :: filename

!CHARACTER(LEN=7), PARAMETER :: filename = "J.h5" ! File name
CHARACTER(LEN=3), PARAMETER :: dsetname1 = "J_r" ! Dataset name
CHARACTER(LEN=3), PARAMETER :: dsetname2 = "J_i" ! Dataset name
!INTEGER, PARAMETER :: NX = 1128

INTEGER(HID_T) :: file_id       ! File identifier
INTEGER(HID_T) :: dset_id1       ! Dataset identifier
INTEGER(HID_T) :: dset_id2       ! Dataset identifier
INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


INTEGER(HSIZE_T), DIMENSION(1) :: dims           ! Dataset dimensions
INTEGER     ::    rank = 1                       ! Dataset rank

INTEGER     ::   error ! Error flag
!INTEGER     :: i, j

filename = fname
dims = (/size(vec)/)
     
CALL h5open_f(error)

! Create a new file using default properties.
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

! Create the dataspace.
CALL h5screate_simple_f(rank, dims, dspace_id, error)

! Create and write dataset using default properties.

CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(vec), dims, error)

CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(vec), dims, error)


! End access to the dataset and release resources used by it.
CALL h5dclose_f(dset_id1, error)
CALL h5dclose_f(dset_id2, error)

! Terminate access to the data space.
CALL h5sclose_f(dspace_id, error)


! Close the file.
CALL h5fclose_f(file_id, error)

! Close FORTRAN interface.
CALL h5close_f(error)

end subroutine cmplx_vec_write2file

!________________________________________________________________________

subroutine real_vec_write2file(vec)
double precision, intent(in) :: vec(:)

CHARACTER(LEN=7), PARAMETER :: filename = "rcs.h5" ! File name
CHARACTER(LEN=3), PARAMETER :: dsetname1 = "rcs" ! Dataset name

INTEGER(HID_T) :: file_id       ! File identifier
INTEGER(HID_T) :: dset_id1       ! Dataset identifier
INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


INTEGER(HSIZE_T), DIMENSION(1) :: dims     ! Dataset dimensions
INTEGER     ::    rank = 1                       ! Dataset rank

INTEGER     ::   error ! Error flag

dims  =  (/size(vec)/)

CALL h5open_f(error)

!
! Create a new file using default properties.
!
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

!
! Create the dataspace.
!
CALL h5screate_simple_f(rank, dims, dspace_id, error)

!
! Create and write dataset using default properties.
!
CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, vec, dims, error)

!
! End access to the dataset and release resources used by it.
!
CALL h5dclose_f(dset_id1, error)

!
! Terminate access to the data space.
!
CALL h5sclose_f(dspace_id, error)

!
! Close the file.
!
CALL h5fclose_f(file_id, error)

!
! Close FORTRAN interface.
!
CALL h5close_f(error)

end subroutine real_vec_write2file

subroutine read_mesh(mesh,file)
type (mesh_struct) :: mesh 
CHARACTER(LEN=28) :: file

!CHARACTER(LEN=8), PARAMETER :: file = "mesh.h5" ! File name
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

end subroutine read_mesh



subroutine read_field_points(matrices)
type (data) :: matrices 

CHARACTER(LEN=18), PARAMETER :: file = "field_points.h5" ! File name
CHARACTER(LEN=16), PARAMETER :: dataset1 = "coord"     ! Dataset name


integer(HID_T) :: file_id 
integer(HID_T) :: dataset1_id
integer(HID_T) :: dataspace_id
integer(HSIZE_T), dimension(2) :: coord_dims, dims_out

integer :: error

double precision, dimension(:,:), allocatable :: coord

call h5open_f(error) ! initialize interface
call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error) ! open file


call h5dopen_f(file_id, dataset1, dataset1_id, error) ! open dataset

call h5dget_space_f(dataset1_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, coord_dims, error)

allocate(coord(coord_dims(1), coord_dims(2)))
call h5dread_f(dataset1_id, H5T_NATIVE_DOUBLE, coord, coord_dims, error)

call h5dclose_f(dataset1_id, error) ! close dataset

call h5fclose_f(file_id, error) ! close file
call h5close_f(error) ! close inteface

matrices%field_points = coord

end subroutine read_field_points



subroutine T_write2file(Taa, Tab, Tba, Tbb, fname)
double complex, intent(in) :: Taa(:,:)
double complex, intent(in) :: Tab(:,:)
double complex, intent(in) :: Tba(:,:)
double complex, intent(in) :: Tbb(:,:)

CHARACTER(LEN=38), intent(in) :: fname
CHARACTER(LEN=38) :: filename

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

INTEGER     ::   error ! Error flag
!INTEGER     :: i, j

filename = fname
dims = (/size(Taa,1),size(Taa,2)/)
     
CALL h5open_f(error)

!
! Create a new file using default properties.
!
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

!
! Create the dataspace.
!
CALL h5screate_simple_f(rank, dims, dspace_id, error)

!
! Create and write dataset using default properties.
!
CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(Taa), dims, error)

CALL h5dclose_f(dset_id1, error)
!____

CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(Taa), dims, error)

CALL h5dclose_f(dset_id2, error)

!******************************************************************

CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id3, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, real(Tab), dims, error)

CALL h5dclose_f(dset_id3, error)
!____

CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id4, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, imag(Tab), dims, error)

CALL h5dclose_f(dset_id4, error)

!******************************************************************

CALL h5dcreate_f(file_id, dsetname5, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id5, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id5, H5T_NATIVE_DOUBLE, real(Tba), dims, error)

CALL h5dclose_f(dset_id5, error)
!____

CALL h5dcreate_f(file_id, dsetname6, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id6, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id6, H5T_NATIVE_DOUBLE, imag(Tba), dims, error)

CALL h5dclose_f(dset_id6, error)

!******************************************************************


CALL h5dcreate_f(file_id, dsetname7, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id7, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id7, H5T_NATIVE_DOUBLE, real(Tbb), dims, error)

CALL h5dclose_f(dset_id7, error)
!____

CALL h5dcreate_f(file_id, dsetname8, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id8, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id8, H5T_NATIVE_DOUBLE, imag(Tbb), dims, error)

CALL h5dclose_f(dset_id8, error)

!******************************************************************
! Terminate access to the data space.
!
CALL h5sclose_f(dspace_id, error)

!
! Close the file.
!
CALL h5fclose_f(file_id, error)

!
! Close FORTRAN interface.
!
CALL h5close_f(error)

end subroutine T_write2file





subroutine read_T(file, Taa, Tab, Tba, Tbb)
!CHARACTER(LEN=28), intent(in) :: fname

CHARACTER(LEN=38) :: file  ! File name
CHARACTER(LEN=16), PARAMETER :: dataset1 = "Taa_r"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset2 = "Taa_i"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset3 = "Tab_r"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset4 = "Tab_i"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset5 = "Tba_r"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset6 = "Tba_i"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset7 = "Tbb_r"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset8 = "Tbb_i"     ! Dataset name

integer(HID_T) :: file_id 
integer(HID_T) :: dataset1_id, dataset2_id, dataset3_id, dataset4_id
integer(HID_T) :: dataset5_id, dataset6_id, dataset7_id, dataset8_id
integer(HID_T) :: dataspace_id

integer(HSIZE_T), dimension(2) :: dims_out, dims

integer :: error

double precision, dimension(:,:), allocatable :: Taa_r, Taa_i
double precision, dimension(:,:), allocatable :: Tab_r, Tab_i
double precision, dimension(:,:), allocatable :: Tba_r, Tba_i
double precision, dimension(:,:), allocatable :: Tbb_r, Tbb_i
double complex, dimension(:,:), allocatable :: Taa, Tab, Tba, Tbb 


call h5open_f(error) ! initialize interface
call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error) ! open file

!__________ _____________________________
call h5dopen_f(file_id, dataset1, dataset1_id, error) ! open dataset
call h5dget_space_f(dataset1_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Taa_r(dims(1), dims(2)))
call h5dread_f(dataset1_id, H5T_NATIVE_DOUBLE, Taa_r,dims, error)
call h5dclose_f(dataset1_id, error) ! close dataset
!__________________________________________
call h5dopen_f(file_id, dataset2, dataset2_id, error) ! open dataset
call h5dget_space_f(dataset2_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Taa_i(dims(1), dims(2)))
call h5dread_f(dataset2_id, H5T_NATIVE_DOUBLE, Taa_i,dims, error)
call h5dclose_f(dataset2_id, error) ! close dataset
!__________________________________________
allocate(Taa(size(Taa_r,1),size(Taa_r,1)))
Taa = dcmplx(Taa_r, Taa_i)


!__________ _____________________________
call h5dopen_f(file_id, dataset3, dataset3_id, error) ! open dataset
call h5dget_space_f(dataset3_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tab_r(dims(1), dims(2)))
call h5dread_f(dataset3_id, H5T_NATIVE_DOUBLE, Tab_r,dims, error)
call h5dclose_f(dataset3_id, error) ! close dataset
!__________________________________________
call h5dopen_f(file_id, dataset4, dataset4_id, error) ! open dataset
call h5dget_space_f(dataset4_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tab_i(dims(1), dims(2)))
call h5dread_f(dataset4_id, H5T_NATIVE_DOUBLE, Tab_i,dims, error)
call h5dclose_f(dataset4_id, error) ! close dataset
!__________________________________________
allocate(Tab(size(Tab_r,1),size(Tab_r,1)))
Tab = dcmplx(Tab_r, Tab_i)


!__________ _____________________________
call h5dopen_f(file_id, dataset5, dataset5_id, error) ! open dataset
call h5dget_space_f(dataset5_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tba_r(dims(1), dims(2)))
call h5dread_f(dataset5_id, H5T_NATIVE_DOUBLE, Tba_r,dims, error)
call h5dclose_f(dataset5_id, error) ! close dataset
!__________________________________________
call h5dopen_f(file_id, dataset6, dataset6_id, error) ! open dataset
call h5dget_space_f(dataset6_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tba_i(dims(1), dims(2)))
call h5dread_f(dataset6_id, H5T_NATIVE_DOUBLE, Tba_i,dims, error)
call h5dclose_f(dataset6_id, error) ! close dataset
!__________________________________________
allocate(Tba(size(Tba_r,1),size(Tba_r,1)))
Tba = dcmplx(Tba_r, Tba_i)




!__________ _____________________________
call h5dopen_f(file_id, dataset7, dataset7_id, error) ! open dataset
call h5dget_space_f(dataset7_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tbb_r(dims(1), dims(2)))
call h5dread_f(dataset7_id, H5T_NATIVE_DOUBLE, Tbb_r,dims, error)
call h5dclose_f(dataset7_id, error) ! close dataset
!__________________________________________
call h5dopen_f(file_id, dataset8, dataset8_id, error) ! open dataset
call h5dget_space_f(dataset8_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tbb_i(dims(1), dims(2)))
call h5dread_f(dataset8_id, H5T_NATIVE_DOUBLE, Tbb_i,dims, error)
call h5dclose_f(dataset8_id, error) ! close dataset
!__________________________________________
allocate(Tbb(size(Tbb_r,1),size(Tbb_r,1)))
Tbb = dcmplx(Tbb_r, Tbb_i)



!_______________________________________________________________

call h5fclose_f(file_id, error) ! close file
call h5close_f(error) ! close inteface


end subroutine read_T





subroutine T_write2file2(Taa, Tab, Tba, Tbb, fname)
double complex, intent(in) :: Taa(:,:,:)
double complex, intent(in) :: Tab(:,:,:)
double complex, intent(in) :: Tba(:,:,:)
double complex, intent(in) :: Tbb(:,:,:)

CHARACTER(LEN=38), intent(in) :: fname
CHARACTER(LEN=38) :: filename

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


INTEGER(HSIZE_T), DIMENSION(3) :: dims  ! Dataset dimensions
INTEGER     ::    rank = 3                       ! Dataset rank

INTEGER     ::   error ! Error flag
!INTEGER     :: i, j

filename = fname
dims = (/size(Taa,1),size(Taa,2),size(Taa,3)/)
     
CALL h5open_f(error)

!
! Create a new file using default properties.
!
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

!
! Create the dataspace.
!
CALL h5screate_simple_f(rank, dims, dspace_id, error)

!
! Create and write dataset using default properties.
!
CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(Taa), dims, error)

CALL h5dclose_f(dset_id1, error)
!____

CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(Taa), dims, error)

CALL h5dclose_f(dset_id2, error)

!******************************************************************

CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id3, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, real(Tab), dims, error)

CALL h5dclose_f(dset_id3, error)
!____

CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id4, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, imag(Tab), dims, error)

CALL h5dclose_f(dset_id4, error)

!******************************************************************

CALL h5dcreate_f(file_id, dsetname5, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id5, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id5, H5T_NATIVE_DOUBLE, real(Tba), dims, error)

CALL h5dclose_f(dset_id5, error)
!____

CALL h5dcreate_f(file_id, dsetname6, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id6, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id6, H5T_NATIVE_DOUBLE, imag(Tba), dims, error)

CALL h5dclose_f(dset_id6, error)

!******************************************************************


CALL h5dcreate_f(file_id, dsetname7, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id7, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id7, H5T_NATIVE_DOUBLE, real(Tbb), dims, error)

CALL h5dclose_f(dset_id7, error)
!____

CALL h5dcreate_f(file_id, dsetname8, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id8, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id8, H5T_NATIVE_DOUBLE, imag(Tbb), dims, error)

CALL h5dclose_f(dset_id8, error)

!******************************************************************
! Terminate access to the data space.
!
CALL h5sclose_f(dspace_id, error)

!
! Close the file.
!
CALL h5fclose_f(file_id, error)

!
! Close FORTRAN interface.
!
CALL h5close_f(error)

end subroutine T_write2file2



end module 

