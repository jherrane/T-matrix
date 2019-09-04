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

subroutine allocate_Ti(matrices)
   type (data) :: matrices
   integer :: Nmax, i, ind1, ind2, nm

   Nmax = maxval(matrices%Nmaxs)

   allocate (matrices%Taai((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1, matrices%bars))
   allocate (matrices%Tabi((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1, matrices%bars))
   allocate (matrices%Tbai((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1, matrices%bars))
   allocate (matrices%Tbbi((Nmax + 1)**2 - 1, (Nmax + 1)**2 - 1, matrices%bars))

   matrices%Taai(:, :, :) = dcmplx(0.0, 0.0)
   matrices%Tabi(:, :, :) = dcmplx(0.0, 0.0)
   matrices%Tbai(:, :, :) = dcmplx(0.0, 0.0)
   matrices%Tbbi(:, :, :) = dcmplx(0.0, 0.0)

end subroutine allocate_Ti

!******************************************************************************

subroutine T_empty(matrices, mesh)
      type (mesh_struct) :: mesh
      type (data) :: matrices
      character(len=80) :: fname
      character(len=80) :: filename

      character(len=6), PARAMETER :: dsetname1 = "Taai_r"
      character(len=6), PARAMETER :: dsetname2 = "Taai_i"
      character(len=6), PARAMETER :: dsetname3 = "Tabi_r"
      character(len=6), PARAMETER :: dsetname4 = "Tabi_i"
      character(len=6), PARAMETER :: dsetname5 = "Tbai_r"
      character(len=6), PARAMETER :: dsetname6 = "Tbai_i"
      character(len=6), PARAMETER :: dsetname7 = "Tbbi_r"
      character(len=6), PARAMETER :: dsetname8 = "Tbbi_i"
      character(len=16), PARAMETER :: dsetname9 = "T-ref-a"
      character(len=16), PARAMETER :: dsetname10 = "T-wavlens"

      integer(HID_T) :: file_id
      integer(HID_T) :: dset_id1
      integer(HID_T) :: dset_id2
      integer(HID_T) :: dset_id3
      integer(HID_T) :: dset_id4
      integer(HID_T) :: dset_id5
      integer(HID_T) :: dset_id6
      integer(HID_T) :: dset_id7
      integer(HID_T) :: dset_id8
      integer(HID_T) :: dset_id9
      integer(HID_T) :: dset_id10

      integer(HID_T) :: dspace_id, dspace_id2, dspace_id3

      integer(HSIZE_T), dimension(1) :: dims, dimswl, dimsinfo
      integer     ::    rank = 1
      integer     ::   error ! Error flag

      complex(dp), dimension(:), allocatable    ::   emptyT

      fname = matrices%tname
      if(file_exists(fname)) return

      filename = fname
      dims = int8((/T_size/))
      dimswl = int8((/matrices%bars/))
      dimsinfo = int8((/3/))
      allocate (emptyT(T_size))
      emptyT = dcmplx(0d0, 0d0)

      CALL h5open_f(error)
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
      CALL h5screate_simple_f(1, dimsinfo, dspace_id2, error)
      CALL h5screate_simple_f(1, dimswl, dspace_id3, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(emptyT), dims, error)
      CALL h5dclose_f(dset_id1, error)
      CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(emptyT), dims, error)
      CALL h5dclose_f(dset_id2, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id3, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, real(emptyT), dims, error)
      CALL h5dclose_f(dset_id3, error)
      CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id4, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, imag(emptyT), dims, error)
      CALL h5dclose_f(dset_id4, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname5, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id5, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id5, H5T_NATIVE_DOUBLE, real(emptyT), dims, error)
      CALL h5dclose_f(dset_id5, error)
      CALL h5dcreate_f(file_id, dsetname6, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id6, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id6, H5T_NATIVE_DOUBLE, imag(emptyT), dims, error)
      CALL h5dclose_f(dset_id6, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname7, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id7, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id7, H5T_NATIVE_DOUBLE, real(emptyT), dims, error)
      CALL h5dclose_f(dset_id7, error)
      CALL h5dcreate_f(file_id, dsetname8, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id8, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id8, H5T_NATIVE_DOUBLE, imag(emptyT), dims, error)
      CALL h5dclose_f(dset_id8, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname9, H5T_NATIVE_DOUBLE, dspace_id2, &
                       dset_id9, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id9, H5T_NATIVE_DOUBLE, [real(mesh%param(1)), &
                                                    imag(mesh%param(1)), mesh%a], dimsinfo, error)
      CALL h5dclose_f(dset_id9, error)

!****************************************************************************80

      CALL h5dcreate_f(file_id, dsetname10, H5T_NATIVE_DOUBLE, dspace_id3, &
                       dset_id10, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5dwrite_f(dset_id10, H5T_NATIVE_DOUBLE, 2d0*pi/mesh%ki, dimswl, error)
      CALL h5dclose_f(dset_id10, error)

!****************************************************************************80

      CALL h5sclose_f(dspace_id, error)
      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error)

   end subroutine T_empty

!****************************************************************************80

   subroutine singleT_write2file(matrices, mesh, which_i)
      type (mesh_struct) :: mesh
      type (data) :: matrices
      character(len=80) :: filename

      character(len=6), PARAMETER :: dsetname1 = "Taai_r"
      character(len=6), PARAMETER :: dsetname2 = "Taai_i"
      character(len=6), PARAMETER :: dsetname3 = "Tabi_r"
      character(len=6), PARAMETER :: dsetname4 = "Tabi_i"
      character(len=6), PARAMETER :: dsetname5 = "Tbai_r"
      character(len=6), PARAMETER :: dsetname6 = "Tbai_i"
      character(len=6), PARAMETER :: dsetname7 = "Tbbi_r"
      character(len=6), PARAMETER :: dsetname8 = "Tbbi_i"

      integer(HID_T) :: file_id
      integer(HID_T) :: dset_id1
      integer(HID_T) :: dset_id2
      integer(HID_T) :: dset_id3
      integer(HID_T) :: dset_id4
      integer(HID_T) :: dset_id5
      integer(HID_T) :: dset_id6
      integer(HID_T) :: dset_id7
      integer(HID_T) :: dset_id8

      integer(HID_T) :: dspace_id

      integer(HSIZE_T), dimension(1) :: dims, dims_out
      real(dp), dimension(:), allocatable :: Ti_r, Ti_i
      integer     ::   error ! Error flag
      integer :: i, ind1, ind2, nm, choosebar
      integer, optional :: which_i

      filename = matrices%tname
      dims = int8((/T_size/))
      allocate(Ti_r(T_size), Ti_i(T_size))

      ind1 = 1
      if(present(which_i)) then
         choosebar = which_i
         do i = 1, matrices%bars
            nm = (matrices%Nmaxs(i)+1)**2-1
            if(i==which_i)then
               ind2 = ind1+nm**2-1
               exit
            else
               ind1 = ind1+nm**2
            end if
         end do
      else
         choosebar = 1
         do i = 1, matrices%bars
            nm = (matrices%Nmaxs(i)+1)**2-1
            if(i==1)then
               ind2 = ind1+nm**2-1
               exit
            else
               ind1 = ind1+nm**2
            end if
         end do
      end if

      call h5open_f(error)
      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

!****************************************************************************80

      call h5dopen_f(file_id, dsetname1, dset_id1, error)
      call h5dget_space_f(dset_id1, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id1, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      Ti_r(ind1:ind2) = reshape(real(matrices%Taai(1:nm, 1:nm, &
         choosebar)), [nm**2])
      CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      call h5dclose_f(dset_id1, error)

      call h5dopen_f(file_id, dsetname2, dset_id2, error)
      call h5dget_space_f(dset_id2, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id2, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      Ti_i(ind1:ind2) = reshape(imag(matrices%Taai(1:nm, 1:nm, &
         choosebar)), [nm**2])
      CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      call h5dclose_f(dset_id2, error)

!****************************************************************************80

      call h5dopen_f(file_id, dsetname3, dset_id3, error)
      call h5dget_space_f(dset_id3, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id3, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      Ti_r(ind1:ind2) = reshape(real(matrices%Tabi(1:nm, 1:nm, &
         choosebar)), [nm**2])
      CALL h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      call h5dclose_f(dset_id3, error)

      call h5dopen_f(file_id, dsetname4, dset_id4, error)
      call h5dget_space_f(dset_id4, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id4, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      Ti_i(ind1:ind2) = reshape(imag(matrices%Tabi(1:nm, 1:nm, &
         choosebar)), [nm**2])
      CALL h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      call h5dclose_f(dset_id4, error)

!****************************************************************************80

      call h5dopen_f(file_id, dsetname5, dset_id5, error)
      call h5dget_space_f(dset_id5, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id5, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      Ti_r(ind1:ind2) = reshape(real(matrices%Tbai(1:nm, 1:nm, &
         choosebar)), [nm**2])
      CALL h5dwrite_f(dset_id5, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      call h5dclose_f(dset_id5, error)

      call h5dopen_f(file_id, dsetname6, dset_id6, error)
      call h5dget_space_f(dset_id6, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id6, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      Ti_i(ind1:ind2) = reshape(imag(matrices%Tbai(1:nm, 1:nm, &
         choosebar)), [nm**2])
      CALL h5dwrite_f(dset_id6, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      call h5dclose_f(dset_id6, error)

!****************************************************************************80

      call h5dopen_f(file_id, dsetname7, dset_id7, error)
      call h5dget_space_f(dset_id7, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id7, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      Ti_r(ind1:ind2) = reshape(real(matrices%Tbbi(1:nm, 1:nm, &
         choosebar)), [nm**2])
      CALL h5dwrite_f(dset_id7, H5T_NATIVE_DOUBLE, Ti_r, dims, error)
      call h5dclose_f(dset_id7, error)

      call h5dopen_f(file_id, dsetname8, dset_id8, error)
      call h5dget_space_f(dset_id8, dspace_id, error)
      call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
      call h5dread_f(dset_id8, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      Ti_i(ind1:ind2) = reshape(imag(matrices%Tbbi(1:nm, 1:nm, &
         choosebar)), [nm**2])
      CALL h5dwrite_f(dset_id8, H5T_NATIVE_DOUBLE, Ti_i, dims, error)
      call h5dclose_f(dset_id8, error)

!****************************************************************************80

      CALL h5sclose_f(dspace_id, error)
      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error)

   end subroutine singleT_write2file


end module 

