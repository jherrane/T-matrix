module T_matrix
use common
use transformation_matrices
use gmres_module 
use field
implicit none 

contains

subroutine compute_T_matrix(matrices, mesh, k, Nmax, Taa, Tab, Tba, Tbb, my_id, N_procs)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: k, Cabs
integer :: Nmax, ierr, N_procs, my_id, block_size

complex(dp) :: mat(3*mesh%N_tet,2*((Nmax+1)**2-1))
integer :: nm, T1, T2, rate, N, N1, N2
complex(dp) :: T_mat(2*((Nmax+1)**2-1), 2*((Nmax+1)**2-1))
complex(dp), dimension((Nmax+1)**2-1,(Nmax+1)**2-1) :: Taa, Tab, Tba, Tbb
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
   call gmres(matrices,mesh)
   !$omp critical
   T_mat(:,nm) = matmul(transpose(conjg(mat)),matrices%x)
   !$omp end critical   
   write(*,'(2(A,I0))') 'Iteration ', nm, '/', size(mat,2)
end do
!$omp end do
!$omp end parallel

nm = (Nmax+1)**2-1

sc = dcmplx(0.0,k**3.0)

Taa = T_mat(1:nm,1:nm) * sc
Tab = T_mat(1:nm, nm+1:2*nm) * sc
Tba = T_mat(nm+1:2*nm,1:nm) * sc
Tbb = T_mat(nm+1:2*nm,nm+1:2*nm) *sc

end subroutine compute_T_matrix

end module T_matrix
