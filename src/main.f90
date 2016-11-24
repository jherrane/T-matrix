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
double precision :: k, tol, khat(3), E0(3), cell_size, phi, theta, ka
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

CHARACTER(LEN=28) :: meshname, fname, mueller_out, J_out, arg_name, arg
CHARACTER(LEN=28) :: mueller_out2
CHARACTER(LEN=38) :: tname
integer :: cont, tet, num_args, i_arg, Nmax, n, Tmat

!integer :: status(MPI_STATUS_SIZE)
!call MPI_INIT(ierr)

call system_clock(wt1,wrate)


!if (ierr .ne. MPI_SUCCESS) then
!   print *,'Error starting MPI program. Terminating.'
!   call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
!end if

!call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
!call MPI_COMM_SIZE (MPI_COMM_WORLD, N_procs, ierr)

!if(my_id == 0) then 

   print *,'*************************************************************'
   print *,'**                                                         **'
   print *,'**               JVIE T-matrix v. 0.1                       **'
   print *,'**                                                         **'
   print *,'*************************************************************'

  
   ! Default arguments
   meshname = 'mesh.h5'
   mueller_out = 'mueller.h5'
   J_out = 'J.h5'
   k = 1
   khat = [0,0,1] 
   E0 = [1,0,0]
   mueller = 1 
   phi = 0.0
   theta = 0.0
   ave = 0
   halton_init = 0
   tol = 1e-5 
   maxit = 50  
   restart = 4
   projector = 'pfft'
   order = 2
   cell_size = 1.8
   near_zone  = 1
   expansion_order = 0
   tname = 'T.h5'
   Tmat = 0

   num_args = command_argument_count()
   do i_arg = 1,num_args,2
      call get_command_argument(i_arg,arg_name)
  
      select case(arg_name)
      
      case('-mesh')
         call get_command_argument(i_arg+1,arg)
         meshname = arg
         print*, 'mesh file is:', meshname
      case('-S_out')
         call get_command_argument(i_arg+1,arg)
         mueller_out = arg
         print*, 'Mueller matrix is written in the file:', meshname
      case('-J_out')
         call get_command_argument(i_arg+1,arg)
         J_out = arg
         print*, 'Solution is written in the file:', J_out
      case('-T_out')
         call get_command_argument(i_arg+1,arg)
         tname = arg
         print*, 'T-matrix is written in the file:', tname
      case('-k')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) k 
      case('-mueller')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) mueller
      case('-phi')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) phi
      case('-theta')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) theta
      case('-ave')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) ave
      case('-halton_int')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) halton_init
      case('-tol')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) tol
      case('-maxit')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) maxit
      case('-restart')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) restart
      case('-projector')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) projector
      case('-order')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) order
      case('-cell_size')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) cell_size
      case('-near_zone')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) near_zone
      case('-expansion_order')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) expansion_order
      case('-Tmatrix')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) Tmat
      case('-help')
         print*, 'Command line parameters' 
         print*, '-mesh mesh.h5      "Read mesh from file"' 
         print*, '-S_out mueller.h5  "Output file: Mueller matrix"'
         print*, '-J_out J.h5        "Output file: Solution coefficients"'
         print*, '-k 1.0             "Wavenumber"'
         print*, '-mueller 1         "Compute mueller matrix (1) yes (0) no"'
         print*, '-phi 0.0           "Incident angel phi"'
         print*, '-theta 0.0         "Incident angel theta"'
         print*, '-ave 0             "Orientation averaging (number of orientations)"'
         print*, '-halton_init 0     "Orientation averaging (beginning of Halton sequence)"'
         print*, '-tol 1e-5          "GMRES tolerance"'
         print*, '-restart 4         "GMRES restart"'
         print*, '-maxit 50          "Maximum number of GMRES iterations"'
         print*, '-cell_size 1.8     "Cell size of the auxiliary grid"'
         print*, '-near_zone 1       "Near_zone distance (in aux-cells)"'
         print*, '-expansion_order 0 "Order of basis and testing functions"'
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

   print*, 'Order of basis functions  =', expansion_order

   if(order > 3 .or. order < 1) then
      print*, 'ERROR: order should be 1, 2 or 3' 
      stop   
   end if
   mesh%M_ex = order
   mesh%k= k
   if(expansion_order > 1 .or. order < 0) then
      print*, 'ERROR: Expansion order should be 0 or 1 ' 
      stop   
   end if

   print*,'Reading mesh...'
   call read_mesh(mesh, meshname) ! mesh.h5
   !call read_field_points(matrices)
   !call reader(mesh)   ! .txt file

   print*,'   Wavenumber           =', real(k)
   print*,'   Wavelength           = ', real(2*pi / k)
   Print*, 'Done'
 
 
   print*,'Intializing FFT... '

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
!
!________________________________________________________________________

if(mesh%order == 0) then 
   allocate(matrices%rhs(3*mesh%N_tet))
   allocate(matrices%x(3*mesh%N_tet))
   allocate(matrices%Ax(3*mesh%N_tet))
end if
if(mesh%order == 1) then 
   allocate(matrices%rhs(4*3*mesh%N_tet))
   allocate(matrices%x(4*3*mesh%N_tet))
   allocate(matrices%Ax(4*3*mesh%N_tet))
end if



print*,'Constructing ', projector,'-projectors...'


call system_clock(T1,rate)

if(mesh%order == 1) then
   call pfft_projection_lin(matrices,mesh)
end if
if(mesh%order == 0) then
   call pfft_projection_const(matrices,mesh)
end if


call system_clock(T2)
print*,'Done in', real(T2-T1)/real(rate), 'seconds'

!__________________________________________________________!


print*,'Constructing the sparse part of the system matrix'


call system_clock(T1,rate)


if(mesh%order == 0) then
   call compute_near_zone_interactions_const(matrices,mesh)
end if
if(mesh%order == 1) then 
   call compute_near_zone_interactions_lin(matrices,mesh)
end if

call system_clock(T2)
print*,'Done in', real(T2-T1)/real(rate), 'seconds'

!**********************************************

if(Tmat == 1) then


   ka = k * (maxval([mesh%Nx,mesh%Ny,mesh%Nz]) *mesh%delta)/2.0

   Nmax = truncation_order(ka) 
   allocate(T_mat(2*((Nmax+1)**2-1), 2*((Nmax+1)**2-1)))

   allocate(Taa((Nmax+1)**2-1, (Nmax+1)**2-1))
   allocate(Tab((Nmax+1)**2-1, (Nmax+1)**2-1))
   allocate(Tba((Nmax+1)**2-1, (Nmax+1)**2-1))
   allocate(Tbb((Nmax+1)**2-1, (Nmax+1)**2-1))
   
   call compute_T_matrix(matrices, mesh, k, Nmax, Taa, Tab, Tba, Tbb)
   
   call T_write2file(Taa, Tab, Tba, Tbb, tname)

   allocate(a_in((Nmax+1)**2-1))
   allocate(b_in((Nmax+1)**2-1))
   call planewave(Nmax, dble(k), a_in, b_in)

   las = 0
   do n = 1,Nmax
      las = las + (2*n+1)**2
   end do

   allocate(rotD(las))
   allocate(indD(las,2))

   call sph_rotation_sparse(0.0d0, -pi/2.0d0, Nmax, rotD, indD)
   a_in2 = sparse_matmul(rotD,indD,a_in,(Nmax+1)**2-1)
   b_in2 = sparse_matmul(rotD,indD,b_in,(Nmax+1)**2-1)
   
   a_nm = matmul(Taa,a_in) + matmul(Tab,b_in)
   b_nm = matmul(Tbb,b_in) + matmul(Tba,a_in)
   
   a_nm2 = matmul(Taa,a_in2) + matmul(Tab,b_in2)
   b_nm2 = matmul(Tbb,b_in2) + matmul(Tba,a_in2)
   
   call mueller_matrix_coeff(a_nm, b_nm, a_nm2, b_nm2, dcmplx(k), 180, 1, Nmax, SS)
   mueller_out2 = 'mueller2.h5'
   call real_write2file(SS,mueller_out2)
!**********************************************

end if

if(mueller == 0) then


call rhs_planewave(matrices, mesh)

print*,'Solving...'
call system_clock(T1,rate)
call gmres(matrices,mesh)
call system_clock(T2)
print*,'Done in', real(T2-T1)/real(rate), 'seconds'



print*,'Computing scatterered fields...'
call compute_rcs(matrices,mesh)

call real_vec_write2file(matrices%rcs)
call cmplx_vec_write2file(matrices%x,J_out)

call compute_fields(matrices,mesh)
fname = "E_fields.h5"
call write2file(matrices%E_field,fname) 


end if
!_______________ Post prosessing __________________________!



if (mueller==1) then

   if(ave == 0) then
      mueller_mat = compute_mueller(matrices, mesh, pi/180.0*phi, pi/180.0*theta, 360, pi/180.0 * 0.0)

   else
      
      mueller_mat = orientation_ave(mesh,matrices,181,ave,halton_init)
   end if
  
   call real_write2file(mueller_mat,mueller_out)



end if



call system_clock(wt2)
print*,'Total wall-time ', real(wt2-wt1)/real(wrate), 'seconds'

end program main
