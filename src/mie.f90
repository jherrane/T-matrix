module mie
use sfunctions
use common
use translations


implicit none

contains




!*****************************************************************
!
! Calculates fields from coefficients
!
! F = electric field
! G = Magnetic field
!
!*******************************************************************

subroutine calc_fields(a_nm, b_nm, k, Po, F, G, inou)
complex(dp), dimension(:) :: a_nm, b_nm
real(dp) :: Po(3)
complex(dp) :: k, F(3), G(3) 
integer :: inou
integer :: Nmax, n, m, ind, mm
real(dp) :: r, theta, phi, vec(3), q, omega
complex(dp) :: kr, alpha, beta, gamma, cc
complex(dp), dimension(:), allocatable :: sphj, sphy, sphh
complex(dp) :: P(3), B(3), C(3), Y, Y1, Y2, M_nm(3), N_nm(3)
real(dp), dimension(:), allocatable :: L, L1, L2

Nmax = int(sqrt(dble(1+size(a_nm))))-1

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
F(:) = dcmplx(0.0,0.0)
G(:) = dcmplx(0.0,0.0)

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
      
      F = F + a_nm(ind) * M_nm + b_nm(ind) * N_nm 
      G = G + k * (a_nm(ind) * N_nm + b_nm(ind) * M_nm)
     
   end do
   
   deallocate(L,L1,L2)
  
end do

F = sph2cart_vec(theta, phi, F)

G = sph2cart_vec(theta, phi, G)
G = G/(dcmplx(0.0,1.0) * omega * mu)

end subroutine calc_fields


!*****************************************************************
!
! Calculates fields from coefficients
!
! F = electric field
! G = Magnetic field
!
!*******************************************************************

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


end module mie

