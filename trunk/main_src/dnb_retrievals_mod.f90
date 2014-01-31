!  $Header: /cvsroot/akh_clavrx_src/dnb_retrievals_mod.f90,v 1.2.2.2 2014/01/26 04:48:33 heidinger Exp $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: dnb_retrievals_mod.f90 (src)
!       dnb_retrievals_mod (program)
!
! PURPOSE: this program computes lunar eflectance
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
!--------------------------------------------------------------------------------------
module dnb_retrievals_mod
   use FILE_UTILITY

   private
   public :: compute_lunar_reflectance


contains
   subroutine compute_lunar_reflectance ( rad_chdnb_input &
             & , solzen, lunzen &
             & , start_year, month,day_of_month,start_time &
             & , moon_phase_angle , moon_illum_frac &
             & , ancil_data_dir &
             & , ref_chdnb_lunar)

      implicit none
  
      integer, parameter :: int2 = selected_int_kind(3)
      integer, parameter :: int4 = selected_int_kind(8)
      integer, parameter:: real4 = selected_real_kind(6,37)

      real , intent(in) , dimension(:,:) :: rad_chdnb_input
      real , intent(in) , dimension(:,:) :: solzen, lunzen
      real ( kind = real4) , intent(in) :: moon_phase_angle , moon_illum_frac
      integer ( kind = int4 ) :: start_time
      integer ( kind = int2) , intent(in) ::  start_year , month , day_of_month
      character ( len = * ), intent(in) :: ancil_data_dir
  
      real, intent ( out ) , dimension(:,:),allocatable :: ref_chdnb_lunar
      real, dimension(:,:) , allocatable :: rad_chdnb
  
      double precision , parameter :: MEAN_EARTH_SUN_DIST = 149598022.6071
      double precision , parameter :: MEAN_EARTH_MOON_DIST = 384400.0
      real , parameter :: EARTH_RADIUS_KM = 6378.140

      double precision :: yyyymmddhh
      double precision :: curr_phase_angle, cos_phase_angle, curr_earth_sun_dist
      double precision :: curr_earth_moon_dist,curr_mean_irrad,cos_weighted_irrad,dum
  
      character( len = 128 ) ::dnb_lut_file
      integer( kind = 4 ) :: num_irrad_tabvals = 181
      integer( kind = 4 ) :: num_dist_tabvals = 184080
      real, dimension(:,:), allocatable :: lunar_irrad_lut
      double precision, dimension( : , : ), allocatable :: dist_phase_lut
      double precision, dimension( : ), allocatable :: phase_array
  
      integer :: float_size = 4
      integer :: double_size = 8
      real :: hour_fraction , phase_fraction
      real :: minute , hour
      real :: lunar_irrad_dnb 
      integer :: i , j , dtg_index , irrad_index , allocatestatus
      double precision :: denorm1 , denorm2 , denorm3 , denorm_factor
      real :: SRF_INTEG = 0.32560294 ! integral of the DNB sensor response function (micron)
      real :: ASTRO_DARK_THRESH = 109.0 ! Sun 19 degrees or more below horizon
      real :: MIN_LUNAR_IRRAD_DNB = 1.0e-5 ! W/m^2 = 1.0e-09 W/cm^2, threshold for doing calcs
 
      real , parameter:: PI = 3.14159
      real ::  DTOR = PI / 180.
  
      logical :: dnb_verbose = .false.
    ! logical :: dnb_verbose = .true.
      character(len=128) :: lunar_irrad_file 
      character(len=128) :: distance_table_file
      integer :: num_pix , num_elem
 
      ! --- executable
      lunar_irrad_file = trim(ancil_data_dir)//'dnb_ancils/lunar_irrad_Mean_DNB.bin'
      distance_table_file=trim(ancil_data_dir)//'dnb_ancils/DIST_2010-2030_double.bin'
      
      num_pix = ubound(rad_chdnb_input,1)
      num_elem = ubound(rad_chdnb_input,2)
      
      allocate ( lunar_irrad_LUT(2,num_irrad_tabvals),dist_phase_LUT(4,num_dist_tabvals))
      allocate ( phase_array(num_irrad_tabvals) )
      allocate ( ref_chdnb_lunar(num_pix,num_elem),rad_chdnb(num_pix, num_elem))
  
      ! check if files are exist
      if (file_exists(trim(lunar_irrad_file)) .eqv. .false. &
         .or. file_exists(trim(distance_table_file)) .eqv. .false.) then
         print *,'DNB Ancilary data are missing, stopping'
         stop 88
      endif
         
      open (unit=1,file=trim(lunar_irrad_file),status="old",action="read",&
            access="direct",form="unformatted",recl=float_size*2*num_irrad_tabvals)
      read (unit=1,rec=1) lunar_irrad_lut
      close (1)
 
      open (unit=1,file=trim(distance_table_file),status="old",action="read",&
             access="direct",form="unformatted",recl=double_size*4*num_dist_tabvals)
      read (unit=1,rec=1) dist_phase_lut
      close (1)

      !***************************************************************
      ! 3. compute toa downwelling lunar irradiance (lunar_irrad_dnb) for current date/time
      !***************************************************************
 
      minute = mod(start_time/1000./60., 60.)
      hour_fraction = minute / 60.0
      hour = floor(start_time/1000./60./60.)
      yyyymmddhh = start_year * 1000000 + month*10000+ day_of_month *100 + hour
  
      dtg_index = index_in_vector(dist_phase_lut(1,:),num_dist_tabvals,yyyymmddhh)
      dtg_index = min(num_dist_tabvals,dtg_index)
      dtg_index = max(1,dtg_index)

      curr_phase_angle    = dist_phase_lut(2,dtg_index) + &
                       hour_fraction*(dist_phase_lut(2,dtg_index+1)-dist_phase_lut(2,dtg_index))
      curr_earth_sun_dist  = dist_phase_lut(3,dtg_index) + &
                       hour_fraction*(dist_phase_lut(3,dtg_index+1)-dist_phase_lut(3,dtg_index))
      curr_earth_moon_dist = dist_phase_lut(4,dtg_index) + &
                       hour_fraction*(dist_phase_lut(4,dtg_index+1)-dist_phase_lut(4,dtg_index))

      if (dnb_verbose) then

         print *,''
         print *,'compare (these two values should be about the same):'
         print *, 'lunar_phase (from viirs granule) = ',moon_phase_angle
         print *, 'curr_phase_angle (from lut) = ',curr_phase_angle
         print *,''
         print *,''
         print *,'dist_phase_lut(:,dtg_index) = ',dist_phase_lut(:,dtg_index)
         print *,'dist_phase_lut(:,dtg_index+1) = ',dist_phase_lut(:,dtg_index+1)
         print *,'hour_fraction = ',hour_fraction
         print *,'curr_phase_angle = ',curr_phase_angle
         print *,'curr_earthsun_dist = ',curr_earth_sun_dist
         print *,'curr_earthmoon_dist = ',curr_earth_moon_dist
         print *,''
         print *,'solar angle min max ', minval(solzen),maxval(solzen) 
      end if
 
      ! b) interpolate lunar_irrad_lut() to get current mean-geometry lunar irradiance pre-convolved to dnb srf 
      phase_fraction = curr_phase_angle - int(curr_phase_angle)
      phase_array = lunar_irrad_lut(1,:)
      irrad_index =  index_in_vector(phase_array,num_irrad_tabvals,curr_phase_angle)
      curr_mean_irrad = lunar_irrad_lut(2,irrad_index) + &
                 &   phase_fraction*(lunar_irrad_lut(2,irrad_index+1)-lunar_irrad_lut(2,irrad_index))

      if (dnb_verbose) then

         print *,'phase_fraction = ',phase_fraction
         print *,'irrad_index = ',irrad_index
         print *,'lunar_irrad_lut(:,irrad_index) = ',lunar_irrad_lut(:,irrad_index)
         print *,'lunar_irrad_lut(:,irrad_index+1) = ',lunar_irrad_lut(:,irrad_index+1)
         print *,'curr_mean_irrad = ',curr_mean_irrad

      end if
 
      !  c) define denormalization parameters to scale irradiance to current sun/earth/moon geometry
      cos_phase_angle = cos(curr_phase_angle * DTOR)
      denorm1 = mean_earth_sun_dist ** 2 + mean_earth_moon_dist ** 2 + &
           2.0 * mean_earth_moon_dist * mean_earth_sun_dist * cos_phase_angle
      denorm2 = curr_earth_sun_dist ** 2 + curr_earth_moon_dist ** 2 + &
           2.0 * curr_earth_moon_dist * curr_earth_sun_dist * cos_phase_angle
      denorm3 = ((mean_earth_moon_dist - earth_radius_km ) / (curr_earth_moon_dist - earth_radius_km)) ** 2.0
      denorm_factor = (denorm1 / denorm2 ) * denorm3

 
      !  d) denormalize mean-geometry irradiance to current geometry
      !     also, convert from mw/m^2-um to w/m^2 (divide by 1000 and multiply by srf_integ)
      lunar_irrad_dnb = curr_mean_irrad * denorm_factor * (srf_integ * 1.0e-03)

      if (dnb_verbose) then
         print *,'curr_mean_irrad (mw/m^2-micron)= ',curr_mean_irrad
         print *,'denorm_factor = ',denorm_factor
         print *,'--> dnb band-integrated lunar irradiance (w/m^2)= ',lunar_irrad_dnb
      end if

      rad_chdnb= rad_chdnb_input  * 1.0e+04
      ref_chdnb_lunar = -999.0
      do i=1,num_pix
         do j=1,num_elem
          
            if (rad_chdnb ( i , j ) > -1.0) ref_chdnb_lunar( i , j ) = 0.0
          
            if (rad_chdnb( i , j ) < 0 .or. &
                  & solzen( i , j ) < astro_dark_thresh .or. &
                  & lunzen( i , j ) > 90.0) cycle 
             
            cos_weighted_irrad = cos(lunzen( i , j ) * DTOR) * lunar_irrad_dnb
             
            if (cos_weighted_irrad > MIN_LUNAR_IRRAD_DNB) then 
               ref_chdnb_lunar( i , j ) = 100.0 * ( PI * rad_chdnb( i , j ) ) / ( cos_weighted_irrad )
 
               if (ref_chdnb_lunar( i , j ) > 150.0) then
                  ref_chdnb_lunar( i , j ) = 150.0
               end if

            else
               ref_chdnb_lunar( i , j ) = 0.0
            end if 

         end do !j
      end do !i
       

   contains
      !
      !
      !
      integer function index_in_vector(xx,n,x)
         integer, intent(in)::n
         double precision, intent(in)::x
         double precision,dimension(:),intent(in)::xx
         integer:: i,jl,jm,ju
         jl=0
         ju=n+1

         do i=1,2*n
            if (ju-jl <= 1) then
               exit
            end if
            jm=(ju+jl)/2
            if ( (xx(n) >= xx(1)) .eqv. (x >= xx(jm)) ) then
               jl=jm
            else
               ju=jm
            end if
         end do

         if (x == xx(1)) then
            index_in_vector =1
         else if (x == xx(n)) then
            index_in_vector = n-1
         else
            index_in_vector = jl
         end if
  
      end function index_in_vector
  
   end subroutine compute_lunar_reflectance
 
end module dnb_retrievals_mod
