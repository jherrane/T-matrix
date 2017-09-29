program main_multiT
use io
use mie


implicit none


CHARACTER(LEN=38) :: arg_name, arg
CHARACTER(LEN=38) :: T_in, T_multi, T_coh, x1, Tname

integer :: num_args, i_arg, Nmax, N, nm, i1
real(8) :: elem_ka, crs(4), Csca_ave, k, Csca_ic_ave
complex(8), dimension(:,:,:), allocatable :: multiT
complex(8), dimension(:,:,:), allocatable :: Taa, Tab, Tba, Tbb
complex(8), dimension(:,:), allocatable :: T2aa, T2ab, T2ba, T2bb
complex(8), dimension(:,:), allocatable :: Taa_c, Tab_c, Tba_c, Tbb_c
complex(8), dimension(:,:), allocatable :: Taa_ic, Tab_ic, Tba_ic, Tbb_ic

   ! Default arguments
   T_in = 'T'
   T_multi = 'T_multi.h5'
   T_coh = 'T_coh.h5'
   elem_ka = 10.0
   k = 1
   N=1
   num_args = command_argument_count()
   do i_arg = 1,num_args,2
      call get_command_argument(i_arg,arg_name)
  
      select case(arg_name)
      
      case('-T_in')
         call get_command_argument(i_arg+1,arg)
          T_in = arg        
      case('-T_multi')
         call get_command_argument(i_arg+1,arg)
         T_multi = arg       
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

allocate(Taa_c(nm,nm))
allocate(Tab_c(nm,nm))
allocate(Tba_c(nm,nm))
allocate(Tbb_c(nm,nm))

allocate(Taa_ic(nm,nm))
allocate(Tab_ic(nm,nm))
allocate(Tba_ic(nm,nm))
allocate(Tbb_ic(nm,nm))


Taa_c(:,:) = dcmplx(0.0,0.0)
Tab_c(:,:) = dcmplx(0.0,0.0)
Tba_c(:,:) = dcmplx(0.0,0.0)
Tbb_c(:,:) = dcmplx(0.0,0.0)

Taa_ic(:,:) = dcmplx(0.0,0.0)
Tab_ic(:,:) = dcmplx(0.0,0.0)
Tba_ic(:,:) = dcmplx(0.0,0.0)
Tbb_ic(:,:) = dcmplx(0.0,0.0)


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

   !Coherent T-matrix
   Taa_c = Taa_c + Taa(:,:,i1)/N
   Tab_c = Tab_c + Tab(:,:,i1)/N
   Tba_c = Tba_c + Tba(:,:,i1)/N
   Tbb_c = Tbb_c + Tbb(:,:,i1)/N

   call tr_T(T2aa, T2bb, T2ab, T2ba, k, crs)

   Csca_ave = Csca_ave + crs(2)/N

   deallocate(T2aa, T2ab, T2ba, T2bb)
end do

call T_write2file2(Taa, Tab, Tba, Tbb, T_multi)
call T_write2file(Taa_c, Tab_c, Tba_c, Tbb_c, T_coh)

print*, 'Csca ave =', Csca_ave

Csca_ic_ave = 0.0
do i1 = 1,N

   Taa_ic = Taa(:,:,i1)-Taa_c
   Tab_ic = Tab(:,:,i1)-Tab_c
   Tba_ic = Tba(:,:,i1)-Tba_c
   Tbb_ic = Tbb(:,:,i1)-Tbb_c

   call tr_T(Taa_ic, Tbb_ic, Tab_ic, Tba_ic, k, crs)

   Csca_ic_ave = Csca_ic_ave + crs(2)/N
end do

print*, 'Csca_ic ave =', Csca_ic_ave
print*, 'kappa_s_ic =', Csca_ic_ave / (4.0/3.0*pi*((elem_ka/k)**3.0))
print*, 'mfp_ic =', (4.0/3.0*pi*((elem_ka/k)**3.0)) / Csca_ic_ave

end program main_multiT
