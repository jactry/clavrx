!  $Header: https://svn.ssec.wisc.edu/repos/aw_clavrx/trunk/cr_physics_mod.f90 4 2014-03-20 21:14:10Z awalther $

! name:                      cr_physics
! function:                   module  houses physical routines 
! description:
! reference:
! calling sequence:
! inputs:
! outputs:
! dependencies:
! restrictions:
! history:                              added  Jan 2013 (AW)
!-----------------------------------------------------------------------------------------------------------------------

module cr_physics_mod
  use cr_constants_module
  implicit none
  private
  public :: wvmr_from_rh_temp_press
  
contains
  ! - transforms relative humidity to water vapor mixing ratio [g/kg  or ppmv]
   elemental real function wvmr_from_rh_temp_press  ( rh , temp, press )
      implicit none
	  real(kind = REAL4), intent(in) :: rh
	  real(kind = REAL4), intent(in) :: temp
	  real(kind = REAL4), intent(in) :: press
	  real(kind = REAL4) :: sat_wvp
	  real(kind = REAL4) :: wvp
	
	
	  sat_wvp = saturation_vapor_pressure ( temp )
	  wvp = sat_wvp * rh /100.
 	  wvmr_from_rh_temp_press = 1000.0 * 0.622 *  ( wvp / ( press - wvp ) )
	
   end function wvmr_from_rh_temp_press
  
   !-----------------------------------------------------------------
   !
   !
   !
   !-----------------------------------------------------------------
   elemental real function saturation_vapor_pressure ( temp , ice_flag )
      implicit none
	  real(kind = REAL4), intent(in) :: temp
	  logical, intent(in) , optional :: ice_flag
	  saturation_vapor_pressure = 6.112 * exp(17.67 * (temp - 273.16) / (temp - 29.66))
	  if ( present ( ice_flag )) saturation_vapor_pressure = &
	           & 6.1078 * exp(21.8745584 * (temp-273.16) / (temp - 7.66))
   end function saturation_vapor_pressure 
  
  !-----------------------------------------------------------------
  !
  !
  !
  !-----------------------------------------------------------------
  elemental real function droplet_number_concentration ( cod , sigma_ext , dz_meters )
    real(kind = REAL4) , intent(in) :: cod
	real(kind = REAL4) , intent(in) :: sigma_ext
	real(kind = REAL4) , intent(in) :: dz_meters
    
	real(kind = REAL4) :: dz_microns
	
	dz_microns = dz_meters * 1.0e6
	
	droplet_number_concentration = cod / ( sigma_ext * dz_microns )
	
  
  end function droplet_number_concentration
  


end module cr_physics_mod
