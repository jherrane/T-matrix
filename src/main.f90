program main
use geometry
use io
use field
use solver
use precorrection
use projection
use rhs
use transformation_matrices
use T_matrix
use mie
use translations
use mpi

implicit none

type (mesh_struct) :: mesh
type (data) :: matrices

integer ::  T1, T2, rate, ios, restart, maxit, wt1, wt2, wrate, &
 order, near_zone, expansion_order, mueller, ave, halton_init
double precision :: k, tol, khat(3), E0(3), cell_size, phi, theta, ka, elem_ka
character(5) :: projector

double precision, dimension(:,:), allocatable :: mueller_mat, mueller_mat_ave

integer, dimension(:), allocatable :: T_cubes
integer :: ierr, N_procs, my_id, rc, M1, M2, i, las, i1, locs, Nbasis, i2, las2
integer, dimension(:), allocatable :: sp_size, disp
complex, dimension(:,:,:), allocatable :: sp_mat
integer, dimension(:,:), allocatable :: sp_ind, M
complex(dp), dimension(:,:), allocatable :: mat, T_mat, Taa, Tab, Tba, Tbb 
complex(dp), dimension(:), allocatable :: a_in, b_in, a_in2, b_in2
complex(dp), dimension(:), allocatable :: a_nm, b_nm, a_nm2, b_nm2
complex(8), dimension(:), allocatable :: rotD
integer, dimension(:,:), allocatable :: indD
real(dp), dimension(:,:), allocatable :: SS

CHARACTER(LEN=38) :: meshname, fname, mueller_out, J_out, arg_name, arg
CHARACTER(LEN=38) :: mueller_out2
CHARACTER(LEN=38) :: tname
integer :: cont, tet, num_args, i_arg, Nmax, n, Tmat

integer :: status(MPI_STATUS_SIZE)
call MPI_INIT(ierr)

call system_clock(wt1,wrate)

if (ierr .ne. MPI_SUCCESS) then
   print *,'Error starting MPI program. Terminating.'
   call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
end if

call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
call MPI_COMM_SIZE (MPI_COMM_WORLD, N_procs, ierr)

if(my_id == 0) then 
   print *,'*************************************************************'
   print *,'**                                                         **'
   print *,'**               JVIE T-matrix v. 0.1                       **'
   print *,'**                                                         **'
   print *,'*************************************************************'
  
   ! Default arguments
   meshname = 'mesh.h5'
   k = 1
   khat = [0,0,1] 
   E0 = [1,0,0]
   mueller = 1 
   phi = 0.0
   theta = 0.0
   ave = 0
   halton_init = 0
   tol = 1e-3 
   maxit = 50  
   restart = 4
   projector = 'pfft'
   order = 2
   cell_size = 1.8
   near_zone  = 1
   expansion_order = 0
   tname = 'T.h5'
   Tmat = 1
   
   num_args = command_argument_count()
   do i_arg = 1,num_args,2
      call get_command_argument(i_arg,arg_name)
  
      select case(arg_name)
      
      case('-mesh')
         call get_command_argument(i_arg+1,arg)
         meshname = arg
         print*, 'mesh file is:', meshname
      case('-T_out')
         call get_command_argument(i_arg+1,arg)
         tname = arg
         print*, 'T-matrix is written in the file:', tname
      case('-k')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) k 
      case('-Tmatrix')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) Tmat

      case('-help')
         print*, 'Command line parameters' 
         print*, '-mesh mesh.h5      "Read mesh from file"' 
         print*, '-k 1.0             "Wavenumber"'
         stop
      case default 
         print '(a,a,/)', 'Unrecognized command-line option: ', arg_name
         stop
      end select
   end do

   mesh%tol = tol
   mesh%restart = restart
   mesh%maxit = maxit
   matrices%khat = khat
   matrices%E0 = dcmplx(E0, 0.0)
   mesh%grid_size = cell_size
   mesh%near_zone = near_zone
   mesh%order = expansion_order
   mesh%M_ex = order
   mesh%k= k

   print*,'Reading mesh...'
   call read_mesh(mesh, meshname) ! mesh.h5

   print*,'   Wavenumber           =', real(k)
   print*,'   Wavelength           = ', real(2*pi / k)
   print*, 'Done'
end if

call MPI_Bcast(mesh%N_tet,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(mesh%tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(mesh%restart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(mesh%maxit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(mesh%grid_size,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(mesh%near_zone,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(mesh%order,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(mesh%M_ex,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(mesh%k,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(mesh%N_node,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_Bcast(matrices%khat,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(matrices%E0,3,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

if(my_id .ne. 0) then
   allocate(mesh%coord(3,mesh%N_node))
   allocate(mesh%etopol(4,mesh%N_tet))
   allocate(mesh%param(mesh%N_tet))
end if

call MPI_Bcast(mesh%coord,size(mesh%coord),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(mesh%param,size(mesh%param),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(mesh%etopol,size(mesh%etopol),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

if(my_id == 0) then 
   print*,'Initializing FFT... '
end if

   call build_grid2(mesh) 
   call build_box(mesh)
   call tetras_in_cubes(mesh)
   call build_G(matrices, mesh, my_id, N_procs)

if(my_id == 0) then 
   print *,'   Grid size            =', (/mesh%Nx,mesh%Ny,mesh%Nz/)
   print *,'   Delta grid           =', real(mesh%delta) 
   print *,'   Number of cubes      =', mesh%N_cubes
   print *,'   Delta cubes          =', real(mesh%box_delta)
   print *,'   Elems. in cube (max) =', mesh%N_tet_cube
   print*,'Done'
end if

!________________________________________________________________________

allocate(matrices%rhs(3*mesh%N_tet))
allocate(matrices%x(3*mesh%N_tet))
allocate(matrices%Ax(3*mesh%N_tet))

if(my_id == 0) then
   print*,'Constructing ', projector,'-projectors...'
end if

call system_clock(T1,rate)
call pfft_projection_const(matrices,mesh,my_id)
call system_clock(T2)

if(my_id == 0) then
   print*,'Done in', real(T2-T1)/real(rate), 'seconds'
end if

call MPI_BARRIER(MPI_COMM_WORLD,ierr)

call MPI_Bcast(matrices%indS,size(matrices%indS),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(matrices%S,size(matrices%S),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(matrices%Sx,size(matrices%Sx),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(matrices%Sy,size(matrices%Sy),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(matrices%Sz,size(matrices%Sz),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

!_________________________________________________________!

if(my_id == 0) then
   print*,'Constructing the sparse part of the system matrix'
end if
call system_clock(T1,rate)
call compute_near_zone_interactions_const(matrices,mesh)
call system_clock(T2)

if(my_id == 0) print*,'Done in', real(T2-T1)/real(rate), 'seconds'

!**********************************************

Nmax = truncation_order(mesh%k*(dble(maxval([mesh%Nx, mesh%Ny, mesh%Nz])) &
                    *mesh%delta)/2.0d0)

allocate(T_mat(2*((Nmax+1)**2-1), 2*((Nmax+1)**2-1)))
allocate(Taa((Nmax+1)**2-1, (Nmax+1)**2-1))
allocate(Tab((Nmax+1)**2-1, (Nmax+1)**2-1))
allocate(Tba((Nmax+1)**2-1, (Nmax+1)**2-1))
allocate(Tbb((Nmax+1)**2-1, (Nmax+1)**2-1))
Taa(:,:) = dcmplx(0.0)
Tab(:,:) = dcmplx(0.0)
Tba(:,:) = dcmplx(0.0)
Tbb(:,:) = dcmplx(0.0)

if(my_id == 0) print*, 'Tmatrix size:', size(Taa,1)*2, size(Taa,2)*2

call compute_T_matrix(matrices, mesh, mesh%k, Nmax, Taa, Tab, Tba, Tbb, my_id, N_procs)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)

call MPI_ALLREDUCE(MPI_IN_PLACE,Taa,size(Taa), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD,ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE,Tab,size(Tab), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD,ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE,Tba,size(Tba), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD,ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE,Tbb,size(Tbb), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD,ierr)

if(my_id == 0) call T_write2file(Taa, Tab, Tba, Tbb, tname)

call system_clock(wt2)
if(my_id == 0) print*,'Total wall-time ', real(wt2-wt1)/real(wrate), 'seconds'

call MPI_FINALIZE(ierr)

end program main
