module projection
use common
use possu
use integration_points
use integrals
use geometry
use sparse

implicit none

contains

subroutine pfft_projection_const(matrices,mesh)
type (mesh_struct), intent(in) :: mesh
type (data) :: matrices

double complex ::  erI
integer, dimension(:), allocatable :: box_nodes
integer ::  i1, N_basis, M_box, M, tet_B, cube,  m1, Nb, face, Nt, i2, box_node, u

double precision, dimension(:,:), allocatable :: Pb, P_tri, P0_tet, P0_tri, P_sph, P0_sph
double precision, dimension(:), allocatable :: wb, W_tri, w0_tet, w0_tri
double precision ::  cp(3), vol_B, tri_coord(3,3), nvec(3), P_node(3),  B_coord(3,4), node1(3), cp_cube(3)
double complex, dimension(:), allocatable :: Q, Q2x, Q2y, Q2z, sigma, sigma2x, sigma2y, sigma2z
double complex, dimension(:,:), allocatable :: W, WW

double complex :: Q2
integer :: B_nodes(4), Ns, ierr

N_basis = 3 * mesh%N_tet
M_box = size(mesh%etopol_box,1) 
M = int(M_box**(1.0/3.0) - 1)
!print*,'   Expansion order      =',M

allocate(box_nodes(M_box))

allocate(sigma(M_box), sigma2x(M_box), sigma2y(M_box), sigma2z(M_box))

allocate(matrices%S(M_box,mesh%N_tet))
allocate(matrices%Sx(M_box,mesh%N_tet))
allocate(matrices%Sy(M_box,mesh%N_tet))
allocate(matrices%Sz(M_box,mesh%N_tet))
allocate(matrices%indS(M_box,mesh%N_tet))

call inttri(P0_tri, w0_tri,1)
call inttetra(P0_tet, w0_tet,1)

call sample_points_cube(P0_sph)

Nb = size(w0_tet)
allocate(Pb(3,Nb))
allocate(wb(Nb))

Nt = size(w0_tri)
allocate(P_tri(3,Nt))
allocate(W_tri(Nt))

Ns = size(P0_sph,2)

allocate(P_sph(3,Ns))

allocate(Q(Ns),Q2x(Ns),Q2y(Ns),Q2z(Ns))
allocate(W(Ns,M_box))
allocate(WW(M_box,Ns))

P_sph = 2.0*mesh%box_delta * P0_sph 
 
!--------------- Compute W-----------------------------------------!
tet_B = 1

B_nodes =  mesh%etopol(:,tet_B)
B_coord = mesh%coord(:,B_nodes)
vol_B = tetra_volume(B_coord)

call linmap_tet(Pb,Wb, B_coord, P0_tet, W0_tet)
cp = (B_coord(:,1) + B_coord(:,2) + B_coord(:,3) +  B_coord(:,4))/4

cube = point_in_cube(mesh,cp)
box_nodes = mesh%etopol_box(:,cube)
node1 = mesh%nodes(:,box_nodes(1))
!print*, node1

cp_cube(1) = node1(1) + mesh%box_delta/2.0
cp_cube(2) = node1(2) + mesh%box_delta/2.0
cp_cube(3) = node1(3) + mesh%box_delta/2.0
  
W(:,:) = dcmplx(0.0, 0.0)

do m1 = 1,Ns
   do u = 1,M_box
      box_node = box_nodes(u)
      P_node = mesh%nodes(:,box_node)
      W(m1,u) = Gr(P_sph(:,m1)+cp_cube, P_node,mesh%k)
   end do
end do

WW=pinv(W)

do tet_B = 1,mesh%N_tet
   B_nodes =  mesh%etopol(:,tet_B)
   B_coord = mesh%coord(:,B_nodes)
   vol_B = tetra_volume(B_coord)

   erI = mesh%param(tet_B) - dcmplx(1.0, 0.0)

   call linmap_tet(Pb,Wb, B_coord, P0_tet, W0_tet)

   cp = (B_coord(:,1) + B_coord(:,2) + B_coord(:,3) +  B_coord(:,4))/4

   cube = point_in_cube(mesh,cp)
   box_nodes = mesh%etopol_box(:,cube)
   node1 = mesh%nodes(:,box_nodes(1))
   
   cp_cube(1) = node1(1) + mesh%box_delta/2.0
   cp_cube(2) = node1(2) + mesh%box_delta/2.0
   cp_cube(3) = node1(3) + mesh%box_delta/2.0

   Q(:) = dcmplx(0.0, 0.0)
   Q2x(:) = dcmplx(0.0, 0.0)
   Q2y(:) = dcmplx(0.0, 0.0)
   Q2z(:) = dcmplx(0.0, 0.0)

   do m1 = 1,Ns
      do i1 = 1,Nb
         Q(m1) = Q(m1) + Gr(P_sph(:,m1)+cp_cube ,Pb(:,i1), mesh%k) * wb(i1)
      end do
      
      do face =1,4
         tri_coord = tetra_face(B_coord,face)
         call linmap_tri(P_tri, W_tri, tri_coord, P0_tri, w0_tri)
         nvec = tri_n_vectors(tri_coord)
         Q2 = dcmplx(0.0, 0.0)

         if(dot_product(nvec,tri_coord(:,1)-B_coord(:,face))<0) then
            print*, 'Something is wrong... (normal vectors / node ordering)'
            stop
         end if

         do i2 = 1, Nt
            Q2 = Q2 + W_tri(i2) * Gr(P_sph(:,m1)+cp_cube, P_tri(:,i2), mesh%k)   
         end do
         Q2x(m1) = Q2x(m1) + nvec(1) * Q2
         Q2y(m1) = Q2y(m1) + nvec(2) * Q2
         Q2z(m1) = Q2z(m1) + nvec(3) * Q2               
      end do                                     
   end do
  
   sigma = matmul(WW,Q)/sqrt(vol_B)
   sigma2x = matmul(WW,Q2x)/sqrt(vol_B) 
   sigma2y = matmul(WW,Q2y)/sqrt(vol_B)
   sigma2z = matmul(WW,Q2z)/sqrt(vol_B) 

   matrices%indS(:,tet_B) = box_nodes
   matrices%S(:,tet_B) = sigma
   matrices%Sx(:,tet_B) = sigma2x
   matrices%Sy(:,tet_B) = sigma2y
   matrices%Sz(:,tet_B) = sigma2z  
end do

end subroutine pfft_projection_const

end module projection
