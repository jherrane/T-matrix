module T_matrix
use common
use transformation_matrices
use gmres_module 
implicit none 

contains

subroutine compute_T_matrix(matrices, mesh, k, Nmax, T_mat, my_id, N_procs)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: k, Cabs
integer :: Nmax, ierr, N_procs, my_id, block_size

complex(dp) :: mat(3*mesh%N_tet,2*((Nmax+1)**2-1))
integer :: nm, T1, T2, rate, N, N1, N2
complex(dp) :: T_mat(2*((Nmax+1)**2-1), 2*((Nmax+1)**2-1))
complex(dp) :: sc

if(my_id==0) print*, 'Compute transformations...'
call system_clock(T1,rate)
call vswf2constant(mesh, dcmplx(k), Nmax, mat)
call system_clock(T2)
if(my_id==0) print*,'Done in', real(T2-T1)/real(rate), 'seconds'

N = size(mat,2)
block_size = int(ceiling(dble(N)/dble(N_procs)))
N1 = 1 + my_id * block_size
N2 =  (my_id + 1) * block_size

if((my_id)*block_size > N) then
  N1 = N
end if

if((my_id+1)*block_size > N) then 
  N2 = N
end if

!$omp parallel num_threads(nthreads) default(private) &
!$omp shared(mesh, matrices, T_mat)
!$omp do schedule(guided,32)
do nm = N1,N2
   matrices%rhs = mat(:,nm) 
   call gmres(matrices,mesh, nm, size(mat,2))
   !$omp critical
   T_mat(:,nm) = dcmplx(0.0,k**3.0)*matmul(transpose(conjg(mat)),matrices%x)
   !$omp end critical   
end do
!$omp end do
!$omp end parallel

end subroutine compute_T_matrix

end module T_matrix
