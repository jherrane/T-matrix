module rhs
use common
use integration_points

implicit none

contains

!**************************************************

subroutine rhs_planewave(matrices, mesh)

type (mesh_struct) :: mesh
type (data) :: matrices


integer :: ele, i1
double precision :: pt(3,4), v, r(3), khat(3), omega
double precision, dimension(:), allocatable :: W, w0_tet
double precision, dimension(:,:), allocatable :: P, P0_tet, shape_tet
double complex :: eri, E(3), E0(3)
double complex, dimension(:), allocatable :: b, I

if(mesh%order == 0) then
   allocate(b(3*mesh%N_tet), I(3))
   !allocate(matrices%rhs(3*mesh%N_tet))
end if

if(mesh%order == 1) then
   allocate(b(12*mesh%N_tet), I(12))
   !allocate(matrices%rhs(12*mesh%N_tet))
end if

omega = mesh%k * 299792458.0

khat = matrices%khat 
E0 = matrices%E0 

call inttetra(P0_tet, w0_tet,3)

allocate(shape_tet(4,size(w0_tet)))

shape_tet(1,:) = 1-P0_tet(1,:)-P0_tet(2,:) -P0_tet(3,:) 
shape_tet(2,:) = P0_tet(1,:)
shape_tet(3,:) = P0_tet(2,:)
shape_tet(4,:) = P0_tet(3,:)

allocate(P(3,size(w0_tet)))
allocate(w(size(w0_tet)))

do ele = 1, mesh%N_tet

   pt = mesh%coord(:,mesh%etopol(:,ele))
   call linmap_tet(P,W, pt, P0_tet, W0_tet)
   v = tetra_volume(pt)

   eri = mesh%param(ele) - dcmplx(1.0,0.0)

   I(:) = dcmplx(0.0, 0.0)

   do i1 = 1,size(w)

      r = P(:,i1)
      !E = eri * E0 * exp(dcmplx(0.0, dot_product(mesh%k*khat,r)))
      E = E0 * exp(dcmplx(0.0, dot_product(mesh%k*khat,r)))
      if(mesh%order == 0) then
         I = I + E * W(i1)
      end if
      if(mesh%order == 1) then
         I(1:4) = I(1:4) + E(1) * shape_tet(:,i1) * W(i1)
         I(5:8) = I(5:8) + E(2) * shape_tet(:,i1) * W(i1)
         I(9:12) = I(9:12) + E(3) * shape_tet(:,i1) * W(i1)
      end if

   end do
  

 if(mesh%order == 0) then
    b(3*(ele-1)+1 : 3*(ele-1)+3) = -dcmplx(0.0, 1.0) * omega*epsilon*I / sqrt(v)
 end if

 if(mesh%order == 1) then
    b(12*(ele-1)+1 : 12*(ele-1)+12) = -dcmplx(0.0, 1.0) * omega*epsilon*I / sqrt(v)
 end if


end do

matrices%rhs = b

end subroutine rhs_planewave

subroutine rhs_vswf(matrices, mesh)

type (mesh_struct) :: mesh
type (data) :: matrices


integer :: ele, i1
double precision :: pt(3,4), v, r(3), khat(3), omega
double precision, dimension(:), allocatable :: W, w0_tet
double precision, dimension(:,:), allocatable :: P, P0_tet, shape_tet
double complex :: eri, E(3), E0(3)
double complex, dimension(:), allocatable :: b, I

if(mesh%order == 0) then
   allocate(b(3*mesh%N_tet), I(3))
   !allocate(matrices%rhs(3*mesh%N_tet))
end if

if(mesh%order == 1) then
   allocate(b(12*mesh%N_tet), I(12))
   !allocate(matrices%rhs(12*mesh%N_tet))
end if

omega = mesh%k * 299792458.0

khat = matrices%khat 
E0 = matrices%E0 

call inttetra(P0_tet, w0_tet,3)

allocate(shape_tet(4,size(w0_tet)))

shape_tet(1,:) = 1-P0_tet(1,:)-P0_tet(2,:) -P0_tet(3,:) 
shape_tet(2,:) = P0_tet(1,:)
shape_tet(3,:) = P0_tet(2,:)
shape_tet(4,:) = P0_tet(3,:)

allocate(P(3,size(w0_tet)))
allocate(w(size(w0_tet)))

do ele = 1, mesh%N_tet

   pt = mesh%coord(:,mesh%etopol(:,ele))
   call linmap_tet(P,W, pt, P0_tet, W0_tet)
   v = tetra_volume(pt)

   eri = mesh%param(ele) - dcmplx(1.0,0.0)

   I(:) = dcmplx(0.0, 0.0)

   do i1 = 1,size(w)

      r = P(:,i1)
      !E = eri * E0 * exp(dcmplx(0.0, dot_product(mesh%k*khat,r)))
      E = E0 * exp(dcmplx(0.0, dot_product(mesh%k*khat,r)))

      if(mesh%order == 0) then
         I = I + E * W(i1)
      end if
      if(mesh%order == 1) then
         I(1:4) = I(1:4) + E(1) * shape_tet(:,i1) * W(i1)
         I(5:8) = I(5:8) + E(2) * shape_tet(:,i1) * W(i1)
         I(9:12) = I(9:12) + E(3) * shape_tet(:,i1) * W(i1)
      end if

   end do
  

 if(mesh%order == 0) then
    b(3*(ele-1)+1 : 3*(ele-1)+3) = -dcmplx(0.0, 1.0) * omega*epsilon*I / sqrt(v)
 end if

 if(mesh%order == 1) then
    b(12*(ele-1)+1 : 12*(ele-1)+12) = -dcmplx(0.0, 1.0) * omega*epsilon*I / sqrt(v)
 end if


end do

matrices%rhs = b

end subroutine rhs_vswf

end module rhs
