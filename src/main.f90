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

call system_clock(wt1,wrate)

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
   print*,'Initializing FFT... '

   call build_grid2(mesh) 
   call build_box(mesh)
   call tetras_in_cubes(mesh)
   call build_G(matrices, mesh, my_id, N_procs)

   print *,'   Grid size            =', (/mesh%Nx,mesh%Ny,mesh%Nz/)
   print *,'   Delta grid           =', real(mesh%delta) 
   print *,'   Number of cubes      =', mesh%N_cubes
   print *,'   Delta cubes          =', real(mesh%box_delta)
   print *,'   Elems. in cube (max) =', mesh%N_tet_cube
print*,'Done'

!________________________________________________________________________

allocate(matrices%rhs(3*mesh%N_tet))
allocate(matrices%x(3*mesh%N_tet))
allocate(matrices%Ax(3*mesh%N_tet))

print*,'Constructing ', projector,'-projectors...'
call system_clock(T1,rate)
call pfft_projection_const(matrices,mesh)
call system_clock(T2)
print*,'Done in', real(T2-T1)/real(rate), 'seconds'
!_________________________________________________________!

print*,'Constructing the sparse part of the system matrix'
call system_clock(T1,rate)
call compute_near_zone_interactions_const(matrices,mesh)
call system_clock(T2)
print*,'Done in', real(T2-T1)/real(rate), 'seconds'

!**********************************************

Nmax = truncation_order(k) 

allocate(T_mat(2*((Nmax+1)**2-1), 2*((Nmax+1)**2-1)))
allocate(Taa((Nmax+1)**2-1, (Nmax+1)**2-1))
allocate(Tab((Nmax+1)**2-1, (Nmax+1)**2-1))
allocate(Tba((Nmax+1)**2-1, (Nmax+1)**2-1))
allocate(Tbb((Nmax+1)**2-1, (Nmax+1)**2-1))

print*, 'Tmatrix size:', size(Taa,1)*2, size(Taa,2)*2

call compute_T_matrix(matrices, mesh, k, Nmax, Taa, Tab, Tba, Tbb)

call T_write2file(Taa, Tab, Tba, Tbb, tname)

call system_clock(wt2)
print*,'Total wall-time ', real(wt2-wt1)/real(wrate), 'seconds'

end program main
