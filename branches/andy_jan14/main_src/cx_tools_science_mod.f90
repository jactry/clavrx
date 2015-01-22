! $Id:$
!
!

module cx_tools_science_mod


contains
   !
   !   This computes a pseudo toa emissivity
   !       by 11um channel, for which we assume a emissivity of 1
   !
   subroutine compute_emis ( bt_11um , rad,  sensor , ch_num &
            & , emis )
      use cx_tools_science_planck_mod
      implicit none
      
      real, intent(in) :: bt_11um (:,:)
      real, intent(in) :: rad (:,:)
      character(len=10) , intent(in) :: sensor
      integer, intent(in) :: ch_num
      
      real,  intent(out) :: emis (:,:)
      
      real, allocatable :: rad_bt11um_in_chn20 (:,:)
      
      ! --- executable
      allocate ( rad_bt11um_in_chn20 ( size(rad,1), size(rad,2) ) )
      
      call planck_bt2rad ( bt_11um , sensor, ch_num , rad_bt11um_in_chn20)
      
      emis = rad / rad_bt11um_in_chn20
      
      deallocate ( rad_bt11um_in_chn20 )
   
   end subroutine
   
   !  ===========================================
   !
   ! ===========================================
   subroutine compute_solar_refl20 ( &
         rad20 &
      , bt_ch31 &
      , solar_f0_ch20 &
      , sol_zen &
      , sun_earth_distance &
      , sensor &
      ,  refl_sol_ch20) 
      
      use cx_tools_science_planck_mod
      implicit none
      
      real, intent(in) :: rad20(:,:)
      real, intent(in) :: bt_ch31(:,:)
      real, intent(in) :: solar_f0_ch20
      real, intent(in) :: sol_zen (:,:)
      real, intent(in) :: sun_earth_distance
      character(len=10) , intent(in) :: sensor
      real, intent(out) :: refl_sol_ch20(:,:)
      real, allocatable :: rad_bt11um_in_chn20 (:,:)
     
      real , allocatable :: rad_to_refl(:,:)
      integer :: dum(2)
      
      dum = shape ( rad20)
      allocate ( rad_to_refl(dum(1),dum(2)))
      rad_to_refl = -999.
      
      allocate ( rad_bt11um_in_chn20(dum(1),dum(2)), source=rad20 )
      rad_bt11um_in_chn20 = -999.
      
      call planck_bt2rad ( bt_ch31 , sensor, 20 , rad_bt11um_in_chn20)
      
      
      call get_rad_refl_ch20_factor ( sensor , sol_zen  , rad_to_refl )
      
      
               ! --->    AW 10/20/2014
         ! ch(20) %ref_toa is the pseudo solar reflectance in 3.9 channels
         !  Rad_obs = Rad_sol + ( 1 - R ) Rad_ch20_ems
         !  Rad_obs = (R * F_0 * mu) / PI + ( 1 - R ) Rad_ch20_ems
         !  == >   R = ( PI (Rad_obs - Rad_ch20_ems )) / ( F_o * mu - PI * Rad_ch20_ems )
         !     see Kaufman and Remer IEEE 1994:
         !   "Detection of  Forests Using Mid-IR Reflectance: An  Application for Aerosol Studies"
         !   http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=297984
         !
      refl_sol_ch20 = 100. * (rad20 - rad_bt11um_in_chn20) / &
          ( 1 / rad_to_refl -  rad_bt11um_in_chn20)
     
      deallocate ( rad_bt11um_in_chn20 )
      deallocate ( rad_to_refl)
   
   end subroutine compute_solar_refl20 
   
   !
   !
   !
   subroutine get_rad_refl_ch20_factor ( sensor , sol_zen  , rad_to_refl )
      implicit none
      character (len =  * )  :: sensor
      real, intent ( in ) :: sol_zen(:,:)
  
      real :: sun_earth_distance
      real, parameter :: PI = 3.14159265
      real, parameter :: DTOR = PI / 180.
      real , intent(out):: rad_to_refl(:,:)
      real :: solar, ew , solar_rad_20
      integer :: start_day
    
      start_day = 100

      !- these are all sensor specifiv values and van be removed an computed in dcomp box!
      sun_earth_distance = 1.0 - 0.016729 * cos ( 0.9856 * ( start_day - 4.0) * DTOR )
   
      select case ( sensor )
      case ( 'GOES-08')
         solar = 2.1734486
         ew = 149.27568
      case ( 'GOES-09')
         solar = 2.2212146
         ew = 152.97360
      case ( 'GOES-10')
         solar = 2.2627592
         ew = 156.03971
      case ( 'GOES-11')
         solar = 2.2609792
         ew = 155.07054   
      case ( 'GOES-12')
         solar = 4.204
         ew = 270.43962  
      case ( 'GOES-13')
         solar = 3.2502353
         ew = 222.80349
      case ( 'GOES-14')
         solar = 3.4876949 
         ew = 237.43905 
      case ( 'GOES-15')
         solar = 3.5111240
         ew = 241.02064
      case ( 'MODIS-AQUA' , 'MODIS-TERRA')
         solar = 1.9553725
         ew = 127.12892
      case ( 'VIIRS')
         solar = 2.2671891
         ew = 140.86442
      case ('NOAA-05')
         solar =  4.1604741 
         ew = 267.14124
      case ( 'NOAA-06')
         solar = 3.9980827
         ew = 254.00470
      case ('NOAA-07')
         solar =  4.5634907 
         ew = 287.01810
      case ( 'NOAA-08')
         solar = 4.0794900
         ew = 262.90717
      case ( 'NOAA-09')
         solar = 4.611
         ew = 288.84289   
      case ('NOAA-10')
         solar =  4.2889941 
         ew = 272.12618
      case ( 'NOAA-11')
         solar = 4.448
         ew = 278.85792 
      case ( 'NOAA-12')
         solar = 4.204
         ew = 270.43962
      case ('NOAA-14')
         solar =  4.448 
         ew = 284.69366
      case ( 'NOAA-15')
         solar = 3.781
         ew = 236.53016    
      case ( 'NOAA-16')
         solar = 3.7372757
         ew = 236.38144
      case ( 'NOAA-17')
         solar = 4.2443861
         ew = 269.79606
      case ( 'NOAA-18')
         solar = 4.0620453
         ew = 259.57508 
      case ( 'NOAA-19')
         solar = 4.1725582
         ew = 265.07816   
      case ('MTSAT-2') 
         solar = 5.1190589
         ew = 322.06623 
      case ('MTSAT-1R')
         solar = 4.3693553
         ew = 280.93528
      case('METOP-A')
         solar = 4.642
         ew = 291.06
      case('METOP-B')
         solar = 3.9731433
         ew = 253.16404
      case ('Meteosat-8')
         solar = 5.3444818
         ew = 365.59826
      case('Meteosat-9')
         solar = 5.4832609
         ew = 374.42649
      case('Meteosat-10')
         solar = 5.4694165
         ew = 374.16782
      case('GOES-16','ABI')
         solar = 3.5111240 ! faked from goes-15
         ew = 241.02064 ! faked
      case('COMS-1')
         solar = 4.8461549
         ew = 306.29122      
      case default
         print*,'missing sensor calibration'
         print*, 'add to get_rad_refl_factor.f90'
         stop
      end select

      solar_rad_20 =  1000.* solar / ew
      rad_to_refl    = PI / cos (sol_zen * DTOR )/ solar_rad_20 / sun_earth_distance ** 2

end subroutine get_rad_refl_ch20_factor
   

   !-------------------------------------------------------------------------
   ! subroutine LOCATE(xx, n, x, j)
   ! Numerical recipes bisection search - x will be between xx(j) and xx(j+1)
   !--------------------------------------------------------------------------
  subroutine LOCATE(xx, n, x, j)

!   Arguments
    integer,                        intent(in)  :: n
    integer,                        intent(out) :: j
    real,               intent(in)  :: x
    real , dimension(:), intent(in)  :: xx

!   Local variables
    integer :: i, jl, jm, ju

    jl = 0
    ju = n + 1
    do i = 1, 2*n
       if (ju-jl <= 1) then
          exit
       endif
       jm = (ju + jl) / 2
       if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if (x == xx(1)) then
       j=1
    else if (x == xx(n)) then
       j = n - 1
    else
       j = jl
    endif

  end subroutine LOCATE
   

end module cx_tools_science_mod
