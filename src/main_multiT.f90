program main_multiT
use io
use mie
use translations

implicit none


CHARACTER(LEN=38) :: arg_name, arg
CHARACTER(LEN=38) :: T_in, T_multi, T_coh, x1, Tname, mueller_file, T_multi_tot

integer :: num_args, i_arg, Nmax, N, nm, i1, N_theta, N_phi, phi, th
real(8) :: elem_ka, crs(4), Csca_ave, k, Csca_ic_ave, Cabs_ave, albedo, Cabs_ic_ave, Cext_ic_ave 
complex(8), dimension(:,:,:), allocatable :: multiT, multi_T_tot
complex(8), dimension(:,:,:), allocatable :: Taa, Tab, Tba, Tbb
complex(8), dimension(:,:,:), allocatable :: Taa2, Tab2, Tba2, Tbb2
complex(8), dimension(:,:), allocatable :: T2aa, T2ab, T2ba, T2bb
complex(8), dimension(:,:), allocatable :: Taa_c, Tab_c, Tba_c, Tbb_c
complex(8), dimension(:,:), allocatable :: Taa_ic, Tab_ic, Tba_ic, Tbb_ic
complex(8), dimension(:), allocatable :: a_nm, b_nm, a_nm2, b_nm2
complex(8), dimension(:), allocatable :: a_in, b_in, a_in2, b_in2
real(8), dimension(:,:), allocatable :: S_out, S_ave, S_out_ave, Cexts
real(8), dimension(:), allocatable :: Cabs_vec
real(8) :: Csca, Cabs_ic, Cextu, Cscau, Cabsu

   ! Default arguments
   T_in = 'T'
   T_multi = 'T_multi_ic.h5'
   T_multi_tot = 'T_multi_tot.h5'
   T_coh = 'T_coh.h5'
   elem_ka = 10.0
   k = 1
   N=1
   N_theta = 181
   N_phi = 64
   mueller_file = 'mueller_elem.h5'

   num_args = command_argument_count()
   do i_arg = 1,num_args,2
      call get_command_argument(i_arg,arg_name)
  
      select case(arg_name)
      
      case('-T_in')
         call get_command_argument(i_arg+1,arg)
          T_in = arg        
      case('-T_multi_ic')
         call get_command_argument(i_arg+1,arg)
         T_multi = arg
      case('-T_multi')
         call get_command_argument(i_arg+1,arg)
         T_multi_tot = arg
         
      case('-T_coh')
         call get_command_argument(i_arg+1,arg)
         T_coh = arg
      case('-elem_ka')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) elem_ka
      case('-N_Tin')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) N
      case('-k')
         call get_command_argument(i_arg+1,arg)
         read(arg,*) k
      case default 
         print '(a,a,/)', 'Unrecognized command-line option: ', arg_name
         stop
      end select
   end do

Nmax = truncation_order(elem_ka) 

nm = (Nmax+1)**2-1

allocate(Taa(nm,nm,N))
allocate(Tab(nm,nm,N))
allocate(Tba(nm,nm,N))
allocate(Tbb(nm,nm,N))

allocate(Taa2(nm,nm,N))
allocate(Tab2(nm,nm,N))
allocate(Tba2(nm,nm,N))
allocate(Tbb2(nm,nm,N))

allocate(Taa_c(nm,nm))
allocate(Tab_c(nm,nm))
allocate(Tba_c(nm,nm))
allocate(Tbb_c(nm,nm))

allocate(Taa_ic(nm,nm))
allocate(Tab_ic(nm,nm))
allocate(Tba_ic(nm,nm))
allocate(Tbb_ic(nm,nm))

allocate(Cabs_vec(N))
allocate(Cexts(2,N))

allocate(a_nm(nm), b_nm(nm))
allocate(a_in(nm), b_in(nm))

Taa_c(:,:) = dcmplx(0.0,0.0)
Tab_c(:,:) = dcmplx(0.0,0.0)
Tba_c(:,:) = dcmplx(0.0,0.0)
Tbb_c(:,:) = dcmplx(0.0,0.0)

Taa_ic(:,:) = dcmplx(0.0,0.0)
Tab_ic(:,:) = dcmplx(0.0,0.0)
Tba_ic(:,:) = dcmplx(0.0,0.0)
Tbb_ic(:,:) = dcmplx(0.0,0.0)

call planewave(Nmax, k, a_in, b_in)


Csca_ave = 0.0
do i1 = 1,N

   write(x1,'(i10)') i1
   Tname = trim(T_in)//'_'//trim(adjustl(x1))//'.h5'

   print*, 'Read T-matrix from file:',trim(Tname)

   call read_T(Tname, T2aa, T2ab, T2ba, T2bb)
   Taa(:,:,i1) = T2aa
   Tab(:,:,i1) = T2ab
   Tba(:,:,i1) = T2ba
   Tbb(:,:,i1) = T2bb

   !print*, 'Tmat size', sqrt(dble(size(Taa(:,1,i1)))+1)-1

   !Coherent T-matrix
   Taa_c = Taa_c + Taa(:,:,i1)/N
   Tab_c = Tab_c + Tab(:,:,i1)/N
   Tba_c = Tba_c + Tba(:,:,i1)/N
   Tbb_c = Tbb_c + Tbb(:,:,i1)/N

   call tr_T(T2aa, T2bb, T2ab, T2ba, k, crs)

   !print*, size(a_nm), size(a_in), size(T2aa,1)
   a_nm = matmul(T2aa,a_in) + matmul(T2ab,b_in)
   b_nm = matmul(T2bb,b_in) + matmul(T2ba,a_in)

   call cross_sections(a_nm, b_nm, a_in, b_in, dcmplx(k,0.0), Nmax, Cextu, Cscau, Cabsu)

   print*, 'Cabs', crs(3), 'Csca', crs(2)

   Csca_ave = Csca_ave + crs(2)/N
   Cabs_ave = Cabs_ave + crs(3)/N

   if(crs(3)<0) print*, 'negative absorption'
   Cabs_vec(i1) = abs(crs(3))

   

   deallocate(T2aa, T2ab, T2ba, T2bb)
end do



call tr_T(Taa_c, Tbb_c, Tab_c, Tba_c, k, crs)

!print*, 'Cext_c ave =', crs(1)
!print*, 'Csca_c ave =', crs(2)

print*, 'Csca ave =', Csca_ave
print*, 'Cabs ave =', Cabs_ave

deallocate(a_in,b_in,a_nm,b_nm)

allocate(a_in(nm), b_in(nm),a_in2(nm),b_in2(nm))
allocate(a_nm(nm), b_nm(nm),a_nm2(nm),b_nm2(nm))
call planewave2(Nmax, k, a_in, b_in, a_in2, b_in2)

allocate(S_out_ave(N_theta*N_phi,18))
Cext_ic_ave = 0.0
Csca_ic_ave = 0.0
S_out_ave(:,:) = 0.0


do i1 = 1,N

   Taa_ic = Taa(:,:,i1)-Taa_c
   Tab_ic = Tab(:,:,i1)-Tab_c
   Tba_ic = Tba(:,:,i1)-Tba_c
   Tbb_ic = Tbb(:,:,i1)-Tbb_c

   call tr_T(Taa_ic, Tbb_ic, Tab_ic, Tba_ic, k, crs)

   Cext_ic_ave = Cext_ic_ave + crs(1)/N
   Csca_ic_ave = Csca_ic_ave + crs(2)/N

   Cexts(1,i1) = crs(2) + Cabs_vec(i1)
   Cexts(2,i1) = crs(2)/(crs(2) + Cabs_vec(i1))


   print*, 'Cext=', Cexts(1,i1), 'al', Cexts(2,i1)

   ! solve scattering 
   a_nm = matmul(Taa_ic,a_in) + matmul(Tab_ic,b_in)
   b_nm = matmul(Tbb_ic,b_in) + matmul(Tba_ic,a_in)

   a_nm2 = matmul(Taa_ic,a_in2) + matmul(Tab_ic,b_in2)
   b_nm2 = matmul(Tbb_ic,b_in2) + matmul(Tba_ic,a_in2)

   call mueller_matrix_coeff(a_nm, b_nm, a_nm2, b_nm2, dcmplx(k,0.0), N_theta, N_phi, Nmax, S_out)

   S_out_ave = S_out_ave + S_out/N


   Taa2(:,:,i1) = Taa(:,:,i1)
   Tab2(:,:,i1) = Tab(:,:,i1)
   Tba2(:,:,i1) = Tba(:,:,i1)
   Tbb2(:,:,i1) = Tbb(:,:,i1)


   Taa(:,:,i1) = Taa(:,:,i1)-Taa_c
   Tab(:,:,i1) = Tab(:,:,i1)-Tab_c
   Tba(:,:,i1) = Tba(:,:,i1)-Tba_c
   Tbb(:,:,i1) = Tbb(:,:,i1)-Tbb_c


   
   deallocate(S_out)
end do

call T_write2file2(Taa, Tab, Tba, Tbb, Cexts, T_multi)
call T_write2file2(Taa2, Tab2, Tba2, Tbb2, Cexts, T_multi_tot)

allocate(S_ave(N_theta,17))

S_ave(:,:) = 0.0
do th = 1,N_theta
   do phi = 1,N_phi
      S_ave(th,1:17) = S_ave(th,1:17) + S_out_ave(th+N_theta*(phi-1),2:18)/N_phi;           
   end do
   
end do

albedo = Csca_ave/(Csca_ave+Cabs_ave)
Cabs_ic_ave = Csca_ic_ave * (1.0/albedo - 1)


S_ave(:,1) = S_ave(:,1)*180/pi;

call real_write2file(S_ave,mueller_file)

print*, 'Csca ave =', Csca_ave
print*, 'Cabs ave =', Cabs_ave
print*, 'Cext ave =', Cabs_ave + Csca_ave
print*, 'Albedo  =', Csca_ic_ave/(Csca_ic_ave+Cabs_ave)
!print*, 'Albedo (old) =', albedo
print*, 'Csca_ic ave =', Csca_ic_ave
print*, 'Cabs_ic ave =', Cabs_ave
print*, 'Cext_ic ave =', Csca_ic_ave + Cabs_ave

print*, 'kappa_s_ic =', Csca_ic_ave / (4.0/3.0*pi*((elem_ka/k)**3.0))

print*, 'mfp_ic =', (4.0/3.0*pi*((elem_ka/k)**3.0)) / (Csca_ic_ave + Cabs_ave)

!print*, 'mfp_ic(old) =', (4.0/3.0*pi*((elem_ka/k)**3.0)) / (Csca_ic_ave + Cabs_ic_ave)




end program main_multiT
