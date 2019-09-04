module transformation_matrices
use common
use sfunctions

implicit none

contains

!******************************************************************************

subroutine vswf2constant(mesh, k, Nmax, mat)
type (mesh_struct) :: mesh
integer :: Nmax, inout
complex(dp), dimension(3,(Nmax+1)**2-1) :: MM_nm, NN_nm
complex(dp), dimension(3*mesh%N_tet, 2*((Nmax+1)**2 -1)) :: mat
complex(dp) :: F(3), G(3), k
integer :: tet, ind1(3), nm
real(dp) :: T_coord(3,4), r(3), vol

nm = (Nmax+1)**2 -1

do tet = 1, mesh%N_tet
   T_coord = mesh%coord(:,mesh%etopol(:,tet))
   r = (T_coord(:,1) +T_coord(:,2) + T_coord(:,3) + T_coord(:,4))/4.0

   call calc_MN(MM_nm, NN_nm, Nmax, k, r, 0)

   vol= tetra_volume(T_coord)

   ind1 = [3*(tet-1)+1,3*(tet-1)+2,3*(tet-1)+3]

   mat(ind1,1:nm) = MM_nm * sqrt(vol)
   mat(ind1,nm+1:2*nm) = NN_nm * sqrt(vol)
end do

end subroutine vswf2constant

!******************************************************************************

subroutine calc_MN(MM_nm, NN_nm, Nmax, k, Po, inou)
real(dp) :: Po(3)
complex(dp) :: k
integer :: inou
integer :: Nmax, n, m, ind, mm
real(dp) :: r, theta, phi, vec(3), q, omega
complex(dp) :: kr, alpha, beta, gamma, cc
complex(dp), dimension(:), allocatable :: sphj, sphy, sphh
complex(dp) :: P(3), B(3), C(3), Y, Y1, Y2, M_nm(3), N_nm(3)
real(dp), dimension(:), allocatable :: L, L1, L2
complex(dp), dimension(3,(Nmax+1)**2-1) :: MM_nm, NN_nm

vec = cart2sph(Po)

r = vec(1)
theta = vec(2)
phi = vec(3)

kr = k*r

omega = real(k)*299792458.0

allocate(sphj(Nmax+2), sphy(Nmax+2), sphh(Nmax+2))

! spherical bessel functions at kr
if(inou == 0) then
   call cspherebessel(Nmax+1,kr, sphj, sphy)
   sphh = sphj
else    
   sphh = sphankel(Nmax+1, kr)
end if

ind = 0

do n = 1, Nmax
   alpha = sphh(n+1)
   beta = sqrt(dble(n*(n+1)))/kr * sphh(n+1)
   gamma = (n+1.0d0)/kr * sphh(n+1) - sphh(n+2)
 
   allocate(L(n+1),L1(n+2),L2(n))

   call legendre2(n,cos(theta),L)  
   call legendre2(n+1,cos(theta),L1)
   call legendre2(n-1,cos(theta),L2)

   q=(sqrt(n*(n+1.0d0)))/((n*2d0+1.0d0)*sin(theta));

   do m = -n, n
      ind = ind+1
      mm = abs(m)
      
      cc = sqrt((2d0*n+1.0d0)*factorial(n-mm)/factorial(n+mm)/(4d0*pi));
     
      ! Unnormalized complex scalar spherical harmonics
      Y=L(mm+1)*exp(dcmplx(0.0, m*phi));
      Y1=L1(mm+1)*exp(dcmplx(0.0, m*phi));
     
      if(mm == n) then
         Y2 = dcmplx(0.0,0.0)
      else 
         Y2 = L2(mm+1)*exp(dcmplx(0.0, m*phi)) 
      end if

      ! vector spherical harmonics
      P(:) = dcmplx(0.0,0.0)
      P(1) = Y

      Y1=Y1*((n-mm+1.0d0)/(n+1.0d0))

      Y2=Y2*(dble(n+mm)/dble(n))

      B(:) = dcmplx(0.0,0.0)
      B(2) = Y1-Y2
      B(3)=((dcmplx(0.0, m*(2*n+1.0)))/(n*(n+1.0)))*Y

      B = B*q

      C(:) = dcmplx(0.0,0.0)
      C(2) = B(3) 
      C(3) = -B(2)

      ! Spherical vector wave functions
      M_nm = cc * alpha * C
      N_nm = cc * (beta*P + gamma * B)
      
      MM_nm(:,ind) = sph2cart_vec(theta, phi, M_nm)
      NN_nm(:,ind) = sph2cart_vec(theta, phi, N_nm)
   end do
   
   deallocate(L,L1,L2)
end do

end subroutine calc_MN


end module transformation_matrices
