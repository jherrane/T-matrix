module common

implicit none

integer, parameter :: dp = selected_real_kind(15, 307)
double precision, parameter :: pi = 3.141592653589793
double precision, parameter :: epsilon = 8.854187817*(10**(-12.0))
double precision, parameter :: mu = 4.0*pi*10.0**(-7.0)
integer, parameter :: nthreads = 24


type mesh_struct
  integer, dimension(:,:), allocatable :: etopol, etopol_box, tetras, etopol_edges, edges
  double precision, dimension(:,:), allocatable :: coord, nodes
  double complex, dimension(:), allocatable :: param 
  integer :: Nx, Ny, Nz, N_cubes, N_tet, N_tet_cube, M_ex, Nx_cube, Ny_cube,&
       Nz_cube, restart, maxit, near_zone, order, N_node
  double precision, dimension(3) :: min_coord  
  double precision :: delta, k, box_delta, tol, grid_size
  integer :: M1_loc, M2_loc, N1_loc, N2_loc
end type mesh_struct


type data
   double complex, dimension(:,:,:), allocatable :: Fg
   double complex, dimension(:,:), allocatable :: S, Sx, Sy, Sz, E_field 
   double complex, dimension(:,:), allocatable :: S_loc, Sx_loc, Sy_loc, Sz_loc
   integer, dimension(:), allocatable :: S_tet_loc
   integer, dimension(:,:), allocatable :: indS, T_ind, indS_loc 
   double complex, dimension(:), allocatable ::  x,Ax,rhs
   double precision, dimension(:), allocatable :: rcs
   double precision, dimension(3) :: khat
   double complex, dimension(3) :: E0
   double precision, dimension(:,:), allocatable :: field_points, T
   complex, dimension(:,:,:), allocatable :: sp_mat
   integer, dimension(:,:), allocatable :: sp_ind
end type data

contains

!***********************************************************

function crossRR(a,b) result(c)
double precision, intent(in) :: a(3), b(3)
double precision :: c(3)
c(1) = a(2)*b(3) - a(3)*b(2)
c(2) = -a(1)*b(3) + a(3)*b(1)
c(3) = a(1)*b(2) - a(2)*b(1)

end function crossRR

!***********************************************************

function crossCC(a,b) result(c)
double complex, intent(in) :: a(3), b(3)
double complex :: c(3)
c(1) = a(2)*b(3) - a(3)*b(2)
c(2) = -a(1)*b(3) + a(3)*b(1)
c(3) = a(1)*b(2) - a(2)*b(1)

end function crossCC

!*************************************************************************

function norm(r, rp) result(c)
implicit none
double precision, dimension(3), intent(in) :: r, rp
double precision :: c
c = sqrt((r(1)-rp(1))**2 + (r(2)-rp(2))**2 + (r(3)-rp(3))**2)
!norm=sqrt(dot_product(r-rp,r-rp))
end function norm

!*************************************************************************
function vlen(r) result(c)
implicit none
double precision, dimension(3), intent(in) :: r
double precision :: c
c = sqrt((r(1))**2 + (r(2))**2 + (r(3))**2)
!norm=sqrt(dot_product(r-rp,r-rp))
end function vlen

!**********************************************************************

function tet_area(coord) result(c)
implicit none
double precision, dimension(3,4), intent(in) :: coord
double precision :: c, tri1, tri2, tri3, tri4


tri1 = vlen(crossRR(coord(:,2)-coord(:,1),coord(:,3)-coord(:,1))) / 2
tri2 = vlen(crossRR(coord(:,3)-coord(:,4),coord(:,2)-coord(:,4))) / 2
tri3 = vlen(crossRR(coord(:,4)-coord(:,3),coord(:,1)-coord(:,3))) / 2
tri4 = vlen(crossRR(coord(:,1)-coord(:,2),coord(:,4)-coord(:,2))) / 2

c = tri1+tri2+tri3+tri4
!norm=sqrt(dot_product(r-rp,r-rp))
end function tet_area

!*************************************************************************


subroutine linmap_tri(P, W, tri_coord, P0, W0)
implicit none
double precision, intent(in) :: tri_coord(3,3)
double precision, intent(in) :: P0(:,:)
double precision, intent(in) :: W0(:)

double precision :: P(3,size(W0)), W(size(W0))

double precision, dimension(3) :: a, b, ww
double precision :: ww1
integer i1, n

n=size(W0)

a=tri_coord(:,2)-tri_coord(:,1);
b=tri_coord(:,3)-tri_coord(:,1);

ww(1) = (a(1)*a(1) +  a(2)*a(2) + a(3)*a(3)) * b(1) 
ww(2) = (a(1)*a(1) +  a(2)*a(2) + a(3)*a(3)) * b(2) 
ww(3) = (a(1)*a(1) +  a(2)*a(2) + a(3)*a(3)) * b(3) 

ww1=sqrt((ww(1)*b(1)+ww(2)*b(2)+ww(3)*b(3))-((a(1)*b(1)+a(2)*b(2)+a(3)*b(3))*(a(1)*b(1)+a(2)*b(2)+a(3)*b(3))));

do i1 = 1,n 
P(:,i1) = tri_coord(:,1) + a * P0(1,i1) + b * P0(2,i1)
end do

W = W0 * ww1 

end subroutine linmap_tri

!*****************************************************************************
function tetra_volume(tet_coord) result(vol)
double precision, intent(in) :: tet_coord(:,:)
double precision :: p41(3), p21(3), p31(3), cr(3)
double precision :: dt, vol
p41 = tet_coord(:,4) - tet_coord(:,1)
p21 = tet_coord(:,2) - tet_coord(:,1)
p31 = tet_coord(:,3) - tet_coord(:,1)

cr = crossRR(p21,p31)
dt = dot_product(p41,cr)
vol = sqrt(dt*dt)/6

end function tetra_volume

!************************************************************************

subroutine linmap_tet(P,W, tet_coord, P0, W0)
double precision, intent(in) :: tet_coord(:,:)
double precision, intent(in) :: P0(:,:)
double precision,  intent(in) :: W0(:)

double precision :: P(3,size(W0)), W(size(W0))
integer :: n, i1 
double precision :: N1, N2, N3, N4
n=size(W0)

do i1 = 1,n
  N1 = 1 - P0(1,i1) - P0(2,i1) - P0(3,i1)
  N2 = P0(1,i1)
  N3 = P0(2,i1)
  N4 = P0(3,i1)
   
  P(:,i1) = tet_coord(:,1) * N1 + tet_coord(:,2) * N2 + tet_coord(:,3) * N3 + tet_coord(:,4) * N4  

end do

W=W0 * 6 * tetra_volume(tet_coord)

end subroutine linmap_tet

!*********************************************************************************

function Gr(rf,rp,k) result(Green)
double precision, intent(in) :: rf(3), rp(3)
double precision, intent(in) :: k
double precision :: R
double complex :: Green
R=norm(rf,rp)
Green = cdexp(dcmplx(0.0,k*R)) / (4*pi*R)
end function Gr

!************************************************************************

function grad_G(rf,rp,k) result(gradG) 
double precision, intent(in) :: rf(3), rp(3), k
double complex :: gradG(3)  

double precision :: R
R=norm(rf,rp)
gradG = cdexp(dcmplx(0.0,k*R)) / (4*pi*R**3) * dcmplx(-1.0,k*R) * (rf-rp)

end function grad_G

!********************************************************************

function tetra_face(tet_coord, face) result(tri_coord)
double precision, intent(in) :: tet_coord(3,4)
integer, intent(in) :: face
double precision :: tri_coord(3,3)

if(face == 1) then
 tri_coord(:,1) = tet_coord(:,2) 
 tri_coord(:,2) = tet_coord(:,3) 
 tri_coord(:,3) = tet_coord(:,4) 
end if
if(face == 2) then
   tri_coord(:,1) = tet_coord(:,1) 
   tri_coord(:,2) = tet_coord(:,4) 
   tri_coord(:,3) = tet_coord(:,3) 
end if
if(face == 3) then
 tri_coord(:,1) = tet_coord(:,1) 
 tri_coord(:,2) = tet_coord(:,2) 
 tri_coord(:,3) = tet_coord(:,4)

end if
if(face == 4) then
 tri_coord(:,1) = tet_coord(:,1) 
 tri_coord(:,2) = tet_coord(:,3) 
 tri_coord(:,3) = tet_coord(:,2)
end if
 
end function tetra_face

!******************************************************************************

function tri_n_vectors(tri_coord) result(nvec)
double precision, intent(in) :: tri_coord(3,3)
double precision :: nvec(3)

double precision :: a1(3), a2(3), n(3), len_n

a1 = tri_coord(:,2) - tri_coord(:,1)
a2 = tri_coord(:,3) - tri_coord(:,1)

n=crossRR(a1,a2)
len_n = sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3)) 

nvec = n/len_n

end function tri_n_vectors


!*********************************************************************************

function tetra_n_vectors(tet_coord) result(n_vectors)
double precision, intent(in) :: tet_coord(3,4)
double precision :: n_vectors(3,4), tri_coord(3,3), n_vec(3), test(3), tt

integer :: n

do n = 1,4
  tri_coord = tetra_face(tet_coord,n)
  n_vec = tri_n_vectors(tri_coord)

  test = tri_coord(:,1) - tet_coord(:,n)
  tt = dot_product(n_vec,test)
  if(tt<0) then
    print*,'Something is wrong'
  end if

  n_vectors(:,n) = n_vec

end do

end function tetra_n_vectors

!**************************************************************************

function vec2arr(mesh,vec) result(arr)
type (mesh_struct) :: mesh 
double complex, dimension(:), intent(in) :: vec
double complex, dimension(:,:,:), allocatable :: arr
integer :: las, m1, l1, k1, las2,t1,t2,rate

las = 1

allocate(arr(2*mesh%Nx, 2*mesh%Ny, 2*mesh%Nz))

arr(:,:,:) = dcmplx(0.0, 0.0) 

do m1 = 1, mesh%Nz
   do l1 = 1, mesh%Ny 
      do k1 = 1, mesh%Nx        
         arr(k1, l1, m1) = vec(las)       
         las = las + 1
      end do
   end do
end do

end function vec2arr

!****************************************************************************

function arr2vec(mesh,arr) result(vec)
type (mesh_struct) :: mesh 
double complex, dimension(:), allocatable :: vec
double complex, dimension(:,:,:), intent(in) :: arr
integer :: las, m1, l1, k1

las = 1
allocate(vec(mesh%Nx*mesh%Ny*mesh%Nz))

do m1 = 1, mesh%Nz
   do l1 = 1, mesh%Ny 
      do k1 = 1, mesh%Nx
         vec(las) = arr(k1,l1,m1)
         las = las + 1

      end do
   end do
end do

!scale IFFT
vec = vec / (8*mesh%Nx*mesh%Ny*mesh%Nz)

end function arr2vec
!*******************************************************************************

subroutine vec2arr2(arr,mesh,vec) 
type (mesh_struct), intent(in) :: mesh 
double complex, dimension(:), intent(in) :: vec
double complex, dimension(:,:,:) :: arr
integer :: las, m1, l1, k1, las2,t1,t2,rate, Nx, Ny, Nz

Nx = mesh%Nx
Ny = mesh%Ny
Nz = mesh%Nz

arr = dcmplx(0.0, 0.0) 


!!$omp parallel num_threads(nthreads) default(none) &
!!$omp firstprivate(Nx, Ny, Nz)  &
!!$omp private(m1,l1,k1,las) &
!!$omp shared(arr, vec)
!!$omp do 
!!do m1 = 1, 2*Nz
!!   do l1 = 1, 2*Ny 
!!      do k1 = 1, 2*Nx   
!!         las = k1 + (l1-1)*Nx + (m1-1)*Ny*Nx
!!         arr(k1, l1, m1) = dcmplx(0.0, 0.0)   
!         !las = las + 1
!!      end do
!!   end do
!!end do
!!$omp end do
!!$omp end parallel



!$omp parallel num_threads(nthreads) default(none) &
!$omp firstprivate(Nx, Ny, Nz)  &
!$omp private(m1,l1,k1,las) &
!$omp shared(arr, vec)
!$omp do 
do m1 = 1, Nz
   do l1 = 1, Ny 
      do k1 = 1, Nx   
         las = k1 + (l1-1)*Nx + (m1-1)*Ny*Nx
         arr(k1, l1, m1) = vec(las)    
         !las = las + 1
      end do
   end do
end do
!$omp end do
!$omp end parallel


end subroutine vec2arr2

function calc_det3(mat) result(det)
double precision, dimension(3,3), intent(in) :: mat
double precision :: det

det = mat(1,1) * (mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)) &
-  mat(1,2) * (mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1)) &
+  mat(1,3) * (mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))

end function calc_det3

subroutine gradshape(rot,dN,T_coord)
double precision, dimension(3,4), intent(in) :: T_coord
double precision, dimension(3,6) :: rot
double precision, dimension(3,4) :: dN

double precision :: J(3,3), detJ, J0(3,3)

J(:,1) = T_coord(:,2) - T_coord(:,1)
J(:,2) = T_coord(:,3) - T_coord(:,1)
J(:,3) = T_coord(:,4) - T_coord(:,1)

detJ = calc_det3(J)

if(detJ <= 0) then 
   print*,'Something is wrong (detJ <= 0)'
   stop
end if 

J0(1,1) = J(2,2)*J(3,3) - J(2,3)*J(3,2);
J0(2,1) = -(J(1,2)*J(3,3) - J(1,3)*J(3,2));
J0(3,1) = (J(1,2)*J(2,3) - J(1,3)*J(2,2));

J0(1,2) = -(J(2,1)*J(3,3) - J(2,3)*J(3,1));
J0(2,2) = (J(1,1)*J(3,3) - J(1,3)*J(3,1));
J0(3,2) = -(J(1,1)*J(2,3) - J(1,3)*J(2,1));

J0(1,3) = J(2,1)*J(3,2) - J(2,2)*J(3,1);
J0(2,3) = -(J(1,1)*J(3,2) - J(1,2)*J(3,1));
J0(3,3) = (J(1,1)*J(2,2) - J(1,2)*J(2,1));

dN(:,2) = J0(:,1) / detJ
dN(:,3) = J0(:,2) / detJ
dN(:,4) = J0(:,3) / detJ
dN(:,1) = - dN(:,2) - dN(:,3) - dN(:,4)

rot(:,1) = 2*crossRR(dN(:,1),dN(:,2))
rot(:,2) = 2*crossRR(dN(:,1),dN(:,3))
rot(:,3) = 2*crossRR(dN(:,1),dN(:,4))
rot(:,4) = 2*crossRR(dN(:,2),dN(:,3))
rot(:,5) = 2*crossRR(dN(:,3),dN(:,4))
rot(:,6) = 2*crossRR(dN(:,2),dN(:,4))

end subroutine gradshape

function sph_unit_vectors(theta,phi) result(vec)
double precision, intent(in) :: theta, phi
double precision :: vec(3,3)

vec(:,1) = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)] ! r_hat
vec(:,2) = [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)] ! theta_hat
vec(:,3) = [-sin(phi), cos(phi), dble(0.0)] ! phi_hat

end function sph_unit_vectors


function real_outer_product(a,b) result(vec)
double precision, dimension(:) :: a, b 
double precision, dimension(:,:), allocatable :: vec
integer :: i1, Na, Nb

Na = size(a)
Nb = size(b)

allocate(vec(Na,Nb))

do i1 = 1,Na
   vec(i1,:) = b*a(i1)
end do

end function real_outer_product


function rotation(x, axis, angle) result(xp)
double precision :: x(3), angle, xp(3)
integer :: axis

if(axis == 1) then
   xp(1) = x(1)
   xp(2) = cos(angle)*x(2) - sin(angle)*x(3)
   xp(3) = sin(angle)*x(2) + cos(angle)*x(3)
end if

if(axis == 2) then
   xp(1) = cos(angle)*x(1) + sin(angle)*x(3)
   xp(2) = x(2)
   xp(3) = -sin(angle)*x(1) + cos(angle)*x(3)
end if

if(axis == 3) then
   xp(1) = cos(angle)*x(1) - sin(angle)*x(2)
   xp(2) = sin(angle)*x(1) + cos(angle)*x(2)
   xp(3) = x(3)
end if

end function rotation


function rotation_matrix(a,b,g) result(rot)
real(dp) :: a, b, g
real(dp) :: rot(3,3)

rot(1,1) = cos(b)*cos(g);
rot(1,2) = cos(a)*sin(g) + sin(a)*sin(b)*cos(g);
rot(1,3) = sin(a)*sin(g) - cos(a)*sin(b)*cos(g);

rot(2,1) = -cos(b)*sin(g);
rot(2,2) = cos(a)*cos(g) - sin(a)*sin(b)*sin(g);
rot(2,3) = sin(a)*cos(g) + cos(a)*sin(b)*sin(g);

rot(3,1) = sin(b);
rot(3,2) = -sin(a)*cos(b);
rot(3,3) = cos(a)*cos(b);

end function rotation_matrix



function halton_seq(index, base) result(num)
integer :: base, index, i
double precision :: num, f
num = 0.0

f = 1.0/dble(base)

i = index

do while(i>0) 
   num = num + f * mod(i,base)
   i = floor(dble(i)/dble(base))
   f = f / dble(base)
end do

end function halton_seq


function neighbour_tet(T_nodes, B_nodes) result(num)
integer :: T_nodes(4), B_nodes(4)
integer :: i, j, num
num = 0
do i=1,4
   do j=1,4
      if(T_nodes(i) == B_nodes(j)) then 
         num = 1
      end if
   end do
end do

end function neighbour_tet


function factorial(a) result(c)
integer :: i1, a
real(dp) :: c, b

b=1.0
do i1 = 1,a
   b = b * i1
end do
c=dble(b)

end function factorial


function cart2sph(x) result(vec)
real(dp), intent(in) :: x(3)
real(dp) :: vec(3)

vec(1) = sqrt(x(1)**2 + x(2)**2 + x(3)**2) ! r
vec(2) = acos(x(3)/vec(1)) ! theta
vec(3) = atan2(x(2),x(1))   

end function cart2sph


function sph2cart(r,theta,phi) result(x)
real(dp), intent(in) :: r, theta, phi
real(dp) :: x(3)

x(1) = r*sin(theta)*cos(phi)
x(2) = r*sin(theta)*sin(phi)
x(3) = r*cos(theta)

end function sph2cart

function sph2cart_vec(theta, phi, vec) result(vec2)
complex(dp), intent(in) :: vec(3)
real(dp), intent(in) :: theta, phi
complex(dp) :: vec2(3)
real(dp) :: H(3,3)

H(1,:) = [sin(theta)*cos(phi), cos(theta)*cos(phi), -sin(phi)]
H(2,:) = [sin(theta)*sin(phi), cos(theta)*sin(phi), cos(phi)]
H(3,:) = [cos(theta), -sin(theta), dble(0.0)];

vec2 = matmul(H,vec)

end function sph2cart_vec


function binomial(n,k) result(c)
integer :: n, k, i1
real(dp) :: c

c = 1.0d0
do i1 = 1,k
   c = c * (n + 1.0d0 - i1)/i1
end do

end function binomial


function rotation_angles(x) result(vec)
real(dp), intent(in) :: x(3)
real(dp) :: vec(3)

vec(1) = sqrt(x(1)**2 + x(2)**2 + x(3)**2) 
vec(3) = atan2(x(2),x(1))     ! phi      
vec(2) = atan2(x(1),x(3))   

end function rotation_angles


function truncation_order(ka) result(Nmax)
real(dp) :: ka 
integer :: Nmax

if(ka > 1) then
   Nmax = floor(ka + 4.0* (ka)**(1.0/3.0)) 
else
   Nmax = 4
end if

end function truncation_order



end module common
