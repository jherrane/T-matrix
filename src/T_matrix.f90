module T_matrix
use common
use transformation_matrices
use gmres_module 

implicit none 

contains

subroutine compute_T_matrix(matrices, mesh, k, Nmax, Taa, Tab, Tba, Tbb)
type (mesh_struct) :: mesh
type (data) :: matrices
real(dp) :: k
integer :: Nmax

complex(dp) :: mat(3*mesh%N_tet,2*((Nmax+1)**2-1))
integer :: nm, T1, T2, rate 
complex(dp) :: T_mat(2*((Nmax+1)**2-1), 2*((Nmax+1)**2-1))
complex(dp), dimension((Nmax+1)**2-1,(Nmax+1)**2-1) :: Taa, Tab, Tba, Tbb
complex(dp) :: sc

print*, 'Compute transformations...'
call system_clock(T1,rate)
call vswf2constant(mesh, dcmplx(k), Nmax, mat)
call system_clock(T2)
print*,'Done in', real(T2-T1)/real(rate), 'seconds'


do nm = 1,size(mat,2)

   matrices%rhs = mat(:,nm)
   call gmres(matrices,mesh)
   T_mat(:,nm) = matmul(transpose(conjg(mat)),matrices%x)

   print*, nm, '/',size(mat,2)
end do

nm = (Nmax+1)**2-1

sc = dcmplx(0.0,k**3.0)

Taa = T_mat(1:nm,1:nm) * sc
Tab = T_mat(1:nm, nm+1:2*nm) * sc
Tba = T_mat(nm+1:2*nm,1:nm) * sc
Tbb = T_mat(nm+1:2*nm,nm+1:2*nm) *sc

end subroutine compute_T_matrix

end module T_matrix
