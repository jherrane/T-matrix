module field
use common
use integration_points
use io
use omp_lib

implicit none

contains 


function compute_far_field_const(matrices, mesh, r) result(E)

type (data), intent(in) :: matrices  
type (mesh_struct), intent(in) :: mesh
double precision, intent(in) :: r(3)
double complex :: E(3)

integer :: Nb, ele, nodes(4), i1
double precision :: B_coord(3,4), vol_B, rp(3), omega, RR
double complex :: I,  G, Ex(3), Ey(3), Ez(3), I1x(3),I1y(3),I1z(3), ikhat(3)

double precision, dimension(:,:), allocatable :: P0, Pb
double precision, dimension(:), allocatable :: w0, wb

RR = (sqrt(dot_product(r,r)))
ikhat = dcmplx(0.0, r/(sqrt(dot_product(r,r))))

omega = mesh%k*299792458.0

call inttetra(P0, w0, 3)


Nb = size(w0)

allocate(Pb(3,Nb), wb(Nb))
E(:) = dcmplx(0.0, 0.0)
Ex(:) = dcmplx(0.0, 0.0)
Ey(:) = dcmplx(0.0, 0.0)
Ez(:) = dcmplx(0.0, 0.0)

do ele = 1, mesh%N_tet

   nodes = mesh%etopol(:,ele)
   B_coord = mesh%coord(:,nodes)   
   vol_B = tetra_volume(B_coord)

   call linmap_tet(Pb, wb, B_coord, P0, w0)

   I = dcmplx(0.0, 0.0)
   I1x(:) = dcmplx(0.0, 0.0)
   I1y(:) = dcmplx(0.0, 0.0)
   I1z(:) = dcmplx(0.0, 0.0)
 
   do i1 = 1, Nb
      rp = Pb(:,i1)
      G = Gr(r, rp,mesh%k)
      !G = cdexp(dcmplx(0.0,mesh%k*vlen(rp)))
      I = I + G*wb(i1)
   end do

   I1x(1) = I
   I1y(2) = I
   I1z(3) = I

   
   Ex = Ex + mesh%k**2 * crossCC(ikhat,crossCC(ikhat, I1x)) * matrices%x(3*(ele-1)+1) / sqrt(vol_B)
   Ey = Ey + mesh%k**2 * crossCC(ikhat,crossCC(ikhat, I1y)) * matrices%x(3*(ele-1)+2) / sqrt(vol_B) 
   Ez = Ez + mesh%k**2 * crossCC(ikhat,crossCC(ikhat, I1z)) * matrices%x(3*(ele-1)+3) / sqrt(vol_B) 

end do

E = -(Ex + Ey + Ez) / dcmplx(0.0, omega*epsilon)

!E = E*cdexp(dcmplx(0.0,mesh%k*RR))/(4*pi*RR)
end function compute_far_field_const

!__________________________________________________________________________

function compute_far_field_lin(matrices, mesh, r) result(E)

type (data), intent(in) :: matrices  
type (mesh_struct), intent(in) :: mesh
double precision, intent(in) :: r(3)
double complex :: E(3)

integer :: Nb, ele, nodes(4), i1, i2
double precision :: B_coord(3,4), vol_B, rp(3), omega
double complex :: I(4),  G, Ex(3), Ey(3), Ez(3), I1x(3),I1y(3),I1z(3), ikhat(3)

double precision, dimension(:,:), allocatable :: P0, Pb, shape_tet
double precision, dimension(:), allocatable :: w0, wb


ikhat = dcmplx(0.0, r/(sqrt(dot_product(r,r))))

omega = mesh%k*299792458.0

call inttetra(P0, w0, 3)

allocate(shape_tet(4,size(w0)))

shape_tet(1,:) = 1-P0(1,:)-P0(2,:) -P0(3,:) 
shape_tet(2,:) = P0(1,:)
shape_tet(3,:) = P0(2,:)
shape_tet(4,:) = P0(3,:)

Nb = size(w0)

allocate(Pb(3,Nb), wb(Nb))
E(:) = dcmplx(0.0, 0.0)
Ex(:) = dcmplx(0.0, 0.0)
Ey(:) = dcmplx(0.0, 0.0)
Ez(:) = dcmplx(0.0, 0.0)
do ele = 1, mesh%N_tet

   nodes = mesh%etopol(:,ele)
   B_coord = mesh%coord(:,nodes)   
   vol_B = tetra_volume(B_coord)

   call linmap_tet(Pb, wb, B_coord, P0, w0)

   I(:) = dcmplx(0.0, 0.0)
   I1x(:) = dcmplx(0.0, 0.0)
   I1y(:) = dcmplx(0.0, 0.0)
   I1z(:) = dcmplx(0.0, 0.0)
 
   do i1 = 1, Nb
      rp = Pb(:,i1)
      G = Gr(r, rp,mesh%k)
      I = I + G*wb(i1)*shape_tet(:,i1)
   end do


   do i2 = 1,4

      I1x(:) = dcmplx(0.0, 0.0)
      I1y(:) = dcmplx(0.0, 0.0)
      I1z(:) = dcmplx(0.0, 0.0)

      I1x(1) = I(i2)
      I1y(2) = I(i2)
      I1z(3) = I(i2)

   
      Ex = Ex + mesh%k**2 * crossCC(ikhat,crossCC(ikhat, I1x)) * matrices%x(12*(ele-1)+i2) / sqrt(vol_B)
      Ey = Ey + mesh%k**2 * crossCC(ikhat,crossCC(ikhat, I1y)) * matrices%x(12*(ele-1)+i2+4) / sqrt(vol_B) 
      Ez = Ez + mesh%k**2 * crossCC(ikhat,crossCC(ikhat, I1z)) * matrices%x(12*(ele-1)+i2+8) / sqrt(vol_B) 


   end do
end do

E = -(Ex + Ey + Ez) / dcmplx(0.0, omega*epsilon)

end function compute_far_field_lin

!__________________________________________________________________________

function compute_near_field(matrices, mesh, r) result(E)

type (data), intent(in) :: matrices  
type (mesh_struct), intent(in) :: mesh
double precision, intent(in) :: r(3)
double complex :: E(3)

integer :: Nb, Nb_tri, ele,face, nodes(4), i1, i2
double precision :: B_coord(3,4), vol_B, B_tri_coord(3,3), B_nvec(3), rp(3), rp2(3), omega, khat(3)
double complex :: I, Ix(3),Iy(3),Iz(3), G, gradG(3), Ex(3), Ey(3), Ez(3), I1x(3),I1y(3),I1z(3), ikhat(3)

double precision, dimension(:,:), allocatable :: P0, Pb, P0_tri, Pb_tri
double precision, dimension(:), allocatable :: w0, wb, w0_tri, wb_tri

khat=[0,0,1]
ikhat(1) = dcmplx(0.0,0.0)
ikhat(2) = dcmplx(0.0,0.0)
ikhat(3) = dcmplx(0.0,1.0)

omega = mesh%k*299792458.0

call inttetra(P0, w0, 5)
call inttri(P0_tri,w0_tri,5)

Nb = size(w0)
Nb_tri = size(w0_tri)
allocate(Pb(3,Nb), wb(Nb), Pb_tri(3,Nb_tri), wb_tri(Nb_tri))
E(:) = dcmplx(0.0, 0.0)
Ex(:) = dcmplx(0.0, 0.0)
Ey(:) = dcmplx(0.0, 0.0)
Ez(:) = dcmplx(0.0, 0.0)
do ele = 1, mesh%N_tet

   nodes = mesh%etopol(:,ele)
   B_coord = mesh%coord(:,nodes)   
   vol_B = tetra_volume(B_coord)

   call linmap_tet(Pb, wb, B_coord, P0, w0)


   I = dcmplx(0.0, 0.0)
   I1x(:) = dcmplx(0.0, 0.0)
   I1y(:) = dcmplx(0.0, 0.0)
   I1z(:) = dcmplx(0.0, 0.0)
 
   do i1 = 1, Nb
      rp = Pb(:,i1)
      G = Gr(r, rp,mesh%k)
      I = I + G*wb(i1)
   end do

   I1x(1) = I
   I1y(2) = I
   I1z(3) = I

   Ix(:) = dcmplx(0.0,0.0)
   Iy(:) = dcmplx(0.0,0.0)
   Iz(:) = dcmplx(0.0,0.0)

   do face = 1,4
      B_tri_coord = tetra_face(B_coord,face)
      B_nvec = tri_n_vectors(B_tri_coord)
      call linmap_tri(Pb_tri, wb_tri, B_tri_coord, P0_tri, W0_tri)

      do i2 = 1,Nb_tri
         rp2 = Pb_tri(:,i2)
         gradG = grad_G(r,rp2,mesh%k)

         Ix = Ix + B_nvec(1) * gradG * wb_tri(i2)
         Iy = Iy + B_nvec(2) * gradG * wb_tri(i2)
         Iz = Iz + B_nvec(3) * gradG * wb_tri(i2)
      end do

   end do
Ex = Ex + (Ix - mesh%k**2 * I1x) * matrices%x(3*(ele-1)+1) / sqrt(vol_B)
Ey = Ey + (Iy - mesh%k**2 * I1y) * matrices%x(3*(ele-1)+2) / sqrt(vol_B) 
Ez = Ez + (Iz - mesh%k**2 * I1z) * matrices%x(3*(ele-1)+3) / sqrt(vol_B) 

end do

E = (Ex + Ey + Ez) / dcmplx(0.0, omega*epsilon)

end function compute_near_field




!_______________________________________________________________________
!
!
!
!

subroutine compute_rcs(matrices,mesh)
type (data) :: matrices
type (mesh_struct), intent(in) :: mesh
integer :: N_theta, i2
double precision :: r(3), RR, lambda, theta, phi
double complex :: E(3)
double precision, dimension(:), allocatable :: rcs


N_theta = 180

allocate(rcs(N_theta+1))

lambda = (2*pi)/mesh%k
RR = 1000*lambda

phi = 0.0!pi/2.0

do i2 = 0, N_theta

   theta = pi * (i2) / (N_theta) 
  
   r(1) = RR*sin(theta)*cos(phi)
   r(2) = RR*sin(theta)*sin(phi)
   r(3) = RR*cos(theta)

   if(mesh%order == 0)then
      E = compute_far_field_const(matrices, mesh, r)
   end if
   if(mesh%order == 1)then
      E = compute_far_field_lin(matrices, mesh, r)
   end if
   
   rcs(i2+1) = 4*pi*RR**2/(lambda**2) * dble(dot_product(E,E))
      
end do
matrices%rcs = rcs
end subroutine compute_rcs

!___________________________________________________________________

subroutine compute_fields(matrices, mesh)
type (data) :: matrices
type (mesh_struct) :: mesh
integer :: N, i1 
double precision :: r(3)
double complex, dimension(:,:), allocatable :: E  
N = size(matrices%field_points,2)

allocate(E(3,N))
do i1 = 1, N
   r = matrices%field_points(:,i1)

   if(mesh%order == 0)then
      E(:,i1) = compute_far_field_const(matrices, mesh, r) 
      !E(:,i1) = compute_near_field(matrices, mesh, r)
   end if
   if(mesh%order == 1)then
      E(:,i1) = compute_far_field_lin(matrices, mesh, r) 
   end if

end do

!call write2file(E) 
matrices%E_field = E

end subroutine compute_fields

end module field
