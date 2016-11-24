module solver
use common
use field
use io
use integration_points
use gmres_module
use rhs

implicit none

contains
function compute_mueller_ave(matrices, mesh, phi, theta, np, mp) result(mueller)
type (data) :: matrices
type (mesh_struct) :: mesh
double precision :: phi, theta, phi_sca
integer :: np, mp
double precision :: inc_vec(3,3), theta_sca, R , abcd2(2,2)
double precision :: obs_point(3), sca_vec(3,3), ksca(3), abcd(2,2), mat(2,2) 
integer :: T1, T2, rate, n, m, las, ind1, ind2
double precision, dimension(:,:), allocatable :: mueller
double complex :: E(3), S(mp*np,5), A, F(mp*np,5)
character(len=18) :: fname

allocate(mueller(mp*np,17))
mueller(:,:) = 0.0
R  = 100000.0 * 2*pi/mesh%k 

inc_vec = sph_unit_vectors(theta,phi)

abcd(1,1) = 1.0!inc_vec(:,2)
abcd(1,2) = 0.0!inc_vec(:,2)
abcd(2,1) = 0.0!inc_vec(:,3)
abcd(2,2) = 1.0!inc_vec(:,3)


matrices%khat = inc_vec(:,1)
matrices%E0 = inc_vec(:,2)

call rhs_planewave(matrices, mesh)

print*,'Solving...'
call system_clock(T1,rate)
call gmres(matrices,mesh)
call system_clock(T2)
print*,'Done in', real(T2-T1)/real(rate), 'seconds'

print*,'Compute fields...'
call system_clock(T1,rate)
!$omp parallel num_threads(nthreads) default(private) &
!$omp firstprivate(mp, np, theta, phi,R, abcd)  &
!$omp shared(matrices, mesh, S, F)
!$omp do  

do m = 0,mp-1

   phi_sca = 2*pi*m/(mp)
   mat(1,1) = cos(phi_sca)
   mat(1,2) = sin(phi_sca)
   mat(2,1) = sin(phi_sca)
   mat(2,2) = -cos(phi_sca)
   abcd2 = matmul(abcd,mat)

  

   do n = 0,np-1
         
      theta_sca = pi*n/(np-1)
      obs_point = sph2cart(R,theta_sca,phi_sca)
      obs_point = rotation(obs_point,2,theta)
      obs_point = rotation(obs_point,3,phi)
     
      if(mesh%order == 0) then
         E = compute_far_field_const(matrices, mesh, obs_point)
      end if
      if(mesh%order == 1) then
         E = compute_far_field_lin(matrices, mesh, obs_point)
      end if

      sca_vec = sph_unit_vectors(theta_sca,phi_sca)
      sca_vec(:,2) =  rotation(sca_vec(:,2),2,theta)
      sca_vec(:,2) =  rotation(sca_vec(:,2),3,phi)

      sca_vec(:,3) =  rotation(sca_vec(:,3),2,theta)
      sca_vec(:,3) =  rotation(sca_vec(:,3),3,phi)

      A = cdexp(dcmplx(0.0, mesh%k*(R))) / (mesh%k*R) 

      las = m*(np) + n+1
      F(las,1) = 1/A  * sum(E * sca_vec(:,2))
      F(las,3) = 1/A  * sum(E * sca_vec(:,3))
     

   end do
end do

!$omp end do
!$omp end parallel
call system_clock(T2)
print*,'Done in', real(T2-T1)/real(rate), 'seconds'

matrices%E0 = inc_vec(:,3)

call rhs_planewave(matrices, mesh)


print*,'Solving...'
call system_clock(T1,rate)
call gmres(matrices,mesh)
call system_clock(T2)
print*,'Done in', real(T2-T1)/real(rate), 'seconds'

print*,'Compute fields...'
call system_clock(T1,rate)



!$omp parallel num_threads(nthreads) default(private) &
!$omp firstprivate(mp, np, theta, phi,R, abcd)  &
!$omp shared(matrices, mesh, S, F)
!$omp do 

do m = 0,mp-1
   phi_sca = 2*pi*m/(mp)
   mat(1,1) = cos(phi_sca)
   mat(1,2) = sin(phi_sca)
   mat(2,1) = sin(phi_sca)
   mat(2,2) = -cos(phi_sca)
   abcd2 = matmul(abcd,mat)

   do n = 0,np-1
      
     
      theta_sca = pi*n/(np-1)
     
      obs_point = sph2cart(R,theta_sca,phi_sca)
      obs_point = rotation(obs_point,2,theta)
      obs_point = rotation(obs_point,3,phi)
    

      if(mesh%order == 0) then
         E = compute_far_field_const(matrices, mesh, obs_point)
      end if
      if(mesh%order == 1) then
         E = compute_far_field_lin(matrices, mesh, obs_point)
      end if
  

      sca_vec = sph_unit_vectors(theta_sca,phi_sca)
      sca_vec(:,2) =  rotation(sca_vec(:,2),2,theta)
      sca_vec(:,2) =  rotation(sca_vec(:,2),3,phi)

      sca_vec(:,3) =  rotation(sca_vec(:,3),2,theta)
      sca_vec(:,3) =  rotation(sca_vec(:,3),3,phi)


      ksca = obs_point/R * mesh%k

      A = cdexp(dcmplx(0.0, mesh%k*(R) )) / (mesh%k*R) 

      las = m*(np) + n+1
      F(las,2) = 1/A  * sum(E * sca_vec(:,2))
      F(las,4) = 1/A  * sum(E * sca_vec(:,3))
      F(las,5) = theta_sca

   end do
end do
!$omp end do
!$omp end parallel

call system_clock(T2)
print*,'Done in', real(T2-T1)/real(rate), 'seconds'

do m = 0,mp-1
   phi_sca = 2*pi*m/(mp)
   mat(1,1) = cos(phi_sca)
   mat(1,2) = sin(phi_sca)
   mat(2,1) = sin(phi_sca)
   mat(2,2) = -cos(phi_sca)
   abcd2 = matmul(abcd,mat)

   ind1 = m*np + 1
   ind2 = (m+1)*np

   S(ind1:ind2,2) =  -dcmplx(0.0, 1.0) * (F(ind1:ind2,1) * abcd2(1,1) + F(ind1:ind2,2) * abcd2(2,1))
   S(ind1:ind2,3) =  dcmplx(0.0, 1.0) * (F(ind1:ind2,1) * abcd2(1,2) + F(ind1:ind2,2) * abcd2(2,2))
   S(ind1:ind2,4) =  dcmplx(0.0, 1.0) * (F(ind1:ind2,3) * abcd2(1,1) + F(ind1:ind2,4) * abcd2(2,1))
   S(ind1:ind2,1) =  -dcmplx(0.0, 1.0) * (F(ind1:ind2,3) * abcd2(1,2) + F(ind1:ind2,4) * abcd2(2,2)) 

   S(ind1:ind2,5) = F(ind1:ind2,5)

end do


mueller(:,1) = real(S(:,5))*180/pi

mueller(:,2) = 0.5 * (abs(S(:,1))**2 + abs(S(:,2))**2 + abs(S(:,3))**2 + abs(S(:,4))**2)
mueller(:,3) = 0.5 * (abs(S(:,2))**2 - abs(S(:,1))**2 + abs(S(:,4))**2 - abs(S(:,3))**2)
mueller(:,4) = -(real(S(:,2)*conjg(S(:,3)) + S(:,1)*conjg(S(:,4))))
mueller(:,5) = imag(S(:,2)*conjg(S(:,3)) - S(:,1)*conjg(S(:,4)))

mueller(:,6) = 0.5 * (abs(S(:,2))**2 - abs(S(:,1))**2 + abs(S(:,3))**2 - abs(S(:,4))**2)
mueller(:,7) = 0.5 * (abs(S(:,1))**2 + abs(S(:,2))**2 - abs(S(:,3))**2 - abs(S(:,3))**2)
mueller(:,8) = real(S(:,2)*conjg(S(:,3)) - S(:,1)*conjg(S(:,4)))
mueller(:,9) = imag(S(:,2)*conjg(S(:,3)) + S(:,1)*conjg(S(:,4)))

mueller(:,10) = real(S(:,2)*conjg(S(:,4)) + S(:,1)*conjg(S(:,3)))
mueller(:,11) = -(real(S(:,2)*conjg(S(:,4)) - S(:,1)*conjg(S(:,3))))
mueller(:,12) = -(real(S(:,1)*conjg(S(:,2)) + S(:,3)*conjg(S(:,4))))
mueller(:,13) = -(imag(S(:,2)*conjg(S(:,1)) + S(:,4)*conjg(S(:,3))))

mueller(:,14) = -(imag(S(:,4)*conjg(S(:,2)) + S(:,1)*conjg(S(:,3))))
mueller(:,15) = -(imag(S(:,4)*conjg(S(:,2)) - S(:,1)*conjg(S(:,3))))
mueller(:,16) = -(imag(S(:,1)*conjg(S(:,2)) - S(:,3)*conjg(S(:,4))))
mueller(:,17) = -(real(S(:,1)*conjg(S(:,2)) - S(:,3)*conjg(S(:,4))))

fname = "mueller.h5"
!call write2file(S,fname)
end function compute_mueller_ave

!____________________________________________________________________________
!
!
!
!____________________________________________________________________________-

function compute_mueller(matrices, mesh, phi, theta, np, phi_sca) result(mueller)
type (data) :: matrices
type (mesh_struct) :: mesh
double precision :: phi, theta, phi_sca
integer :: np
double precision :: inc_vec(3,3), theta_sca, R, abcd(2,2), mat(2,2) 
double precision :: obs_point(3), sca_vec(3,3), abcd2(2,2)
integer :: T1, T2, rate, n
double precision, dimension(:,:), allocatable :: mueller
double complex :: E(3), S(np,5), A, F(np,5)
character(len = 18) :: fname

allocate(mueller(np,17))
mueller(:,:) = 0.0
R  = 100000.0 * 2*pi/mesh%k 

inc_vec = sph_unit_vectors(theta,phi)

abcd(1,1) = 1.0!inc_vec(:,2)
abcd(1,2) = 0.0!inc_vec(:,2)
abcd(2,1) = 0.0!inc_vec(:,3)
abcd(2,2) = 1.0!inc_vec(:,3)


mat(1,1) = cos(phi_sca)
mat(1,2) = sin(phi_sca)
mat(2,1) = sin(phi_sca)
mat(2,2) = -cos(phi_sca)
abcd2 = matmul(abcd,mat)



matrices%khat = inc_vec(:,1)
matrices%E0 = inc_vec(:,2)

call rhs_planewave(matrices, mesh)

print*,'Solving...'
call system_clock(T1,rate)
call gmres(matrices,mesh)
call system_clock(T2)
print*,'Done in', real(T2-T1)/real(rate), 'seconds'



do n = 0,np-1

   theta_sca = pi*n/(np-1)
 
   obs_point = sph2cart(R,theta_sca,phi_sca)
   obs_point = rotation(obs_point,2,theta)
   obs_point = rotation(obs_point,3,phi)

   if(mesh%order == 0) then
      E = compute_far_field_const(matrices, mesh, obs_point)
   end if
   if(mesh%order == 1) then
      E = compute_far_field_lin(matrices, mesh, obs_point)
   end if
   
   sca_vec = sph_unit_vectors(theta_sca,phi_sca)
   sca_vec(:,2) =  rotation(sca_vec(:,2),2,theta)
   sca_vec(:,2) =  rotation(sca_vec(:,2),3,phi)

   sca_vec(:,3) =  rotation(sca_vec(:,3),2,theta)
   sca_vec(:,3) =  rotation(sca_vec(:,3),3,phi)

   A = cdexp(dcmplx(0.0, mesh%k*(R) )) / (mesh%k*R) 

   F(n+1,1) = 1/A  * sum(E * sca_vec(:,2)) !f11
   F(n+1,3) = 1/A  * sum(E * sca_vec(:,3)) !ff21 
end do

matrices%E0 = inc_vec(:,3)

call rhs_planewave(matrices, mesh)


print*,'Solving...'
call system_clock(T1,rate)
call gmres(matrices,mesh)
call system_clock(T2)
print*,'Done in', real(T2-T1)/real(rate), 'seconds'

do n = 0,np-1

   theta_sca = pi*n/(np-1)
  
   obs_point = sph2cart(R,theta_sca,phi_sca)
   obs_point = rotation(obs_point,2,theta)
   obs_point = rotation(obs_point,3,phi)

   if(mesh%order == 0) then
      E = compute_far_field_const(matrices, mesh, obs_point)
   end if
   if(mesh%order == 1) then
      E = compute_far_field_lin(matrices, mesh, obs_point)
   end if
  
   sca_vec = sph_unit_vectors(theta_sca,phi_sca)
   sca_vec(:,2) =  rotation(sca_vec(:,2),2,theta)
   sca_vec(:,2) =  rotation(sca_vec(:,2),3,phi)

   sca_vec(:,3) =  rotation(sca_vec(:,3),2,theta)
   sca_vec(:,3) =  rotation(sca_vec(:,3),3,phi)

   A = cdexp(dcmplx(0.0, mesh%k*(R) )) / (mesh%k*R) 

   F(n+1,2) = 1/A  * sum(E * sca_vec(:,2)) !f12
   F(n+1,4) = 1/A  * sum(E * sca_vec(:,3)) !f22
   S(n+1,5) = theta_sca
end do



S(:,2) =  -dcmplx(0.0, 1.0) * (F(:,1) * abcd2(1,1) + F(:,2) * abcd2(2,1))
S(:,3) =  dcmplx(0.0, 1.0) * (F(:,1) * abcd2(1,2) + F(:,2) * abcd2(2,2))
S(:,4) =  dcmplx(0.0, 1.0) * (F(:,3) * abcd2(1,1) + F(:,4) * abcd2(2,1))
S(:,1) =  -dcmplx(0.0, 1.0) * (F(:,3) * abcd2(1,2) + F(:,4) * abcd2(2,2))

mueller(:,1) = real(S(:,5))*180/pi

mueller(:,2) = 0.5 * (abs(S(:,1))**2 + abs(S(:,2))**2 + abs(S(:,3))**2 + abs(S(:,4))**2)
mueller(:,3) = 0.5 * (abs(S(:,2))**2 - abs(S(:,1))**2 + abs(S(:,4))**2 - abs(S(:,3))**2)
mueller(:,4) = -(real(S(:,2)*conjg(S(:,3)) + S(:,1)*conjg(S(:,4))))
mueller(:,5) = imag(S(:,2)*conjg(S(:,3)) - S(:,1)*conjg(S(:,4)))

mueller(:,6) = 0.5 * (abs(S(:,2))**2 - abs(S(:,1))**2 + abs(S(:,3))**2 - abs(S(:,4))**2)
mueller(:,7) = 0.5 * (abs(S(:,1))**2 + abs(S(:,2))**2 - abs(S(:,3))**2 - abs(S(:,3))**2)
mueller(:,8) = real(S(:,2)*conjg(S(:,3)) - S(:,1)*conjg(S(:,4)))
mueller(:,9) = imag(S(:,2)*conjg(S(:,3)) + S(:,1)*conjg(S(:,4)))

mueller(:,10) = real(S(:,2)*conjg(S(:,4)) + S(:,1)*conjg(S(:,3)))
mueller(:,11) = -(real(S(:,2)*conjg(S(:,4)) - S(:,1)*conjg(S(:,3))))
mueller(:,12) = -(real(S(:,1)*conjg(S(:,2)) + S(:,3)*conjg(S(:,4))))
mueller(:,13) = -(imag(S(:,2)*conjg(S(:,1)) + S(:,4)*conjg(S(:,3))))

mueller(:,14) = -(imag(S(:,4)*conjg(S(:,2)) + S(:,1)*conjg(S(:,3))))
mueller(:,15) = -(imag(S(:,4)*conjg(S(:,2)) - S(:,1)*conjg(S(:,3))))
mueller(:,16) = -(imag(S(:,1)*conjg(S(:,2)) - S(:,3)*conjg(S(:,4))))
mueller(:,17) = -(real(S(:,1)*conjg(S(:,2)) - S(:,3)*conjg(S(:,4))))

fname = "mueller.h5"
!call write2file(S,fname)

end function compute_mueller

!______ Orientation averaging routine_____________1

function orientation_ave(mesh,matrices, N_theta,M,halton_init) result(mueller_ave)
type (data) :: matrices
type (mesh_struct) :: mesh

double precision, dimension(:,:), allocatable :: mueller, mueller_ave
integer :: M,N, i, N_theta, i1, halton_init
double precision :: vec(3)
 
N = 16!int(floor(sqrt(dble(M))*2))

allocate(mueller_ave(N_theta,17))
mueller_ave(:,:) = 0.0 

do i = 1,M

   vec(2) = acos(2*halton_seq(halton_init+i, 2)-1)
   vec(3) = halton_seq(halton_init+i, 3)*2*pi

   print*,'theta=',180/pi*vec(2), 'phi=',180/pi*vec(3)  
   
   mueller = compute_mueller_ave(matrices, mesh, vec(3), vec(2), N_theta, N) 

   do i1 = 1,N
      mueller_ave = mueller_ave + mueller(N_theta*(i1-1)+1:N_theta*i1,:)/dble(N*M)
   end do
   
   print*,'Orientation averaging', i,'/',M 
end do


end function orientation_ave

end module solver
