program test_dcomp_lut

use dcomp_lut_mod

type(lut_type) :: d
type( lut_output) :: out
character(len=1020) :: ancil_path


ancil_path = '/DATA/Ancil_Data/clavrx_ancil_data/luts/cld/'
call d.initialize('VIIRS',trim(ancil_path))
call d.set_angles(33.,29.,120.)
call d.get_data(1, 1, 0.23,0.1, out)

print*,out%refl, out%trn_sol
end program test_dcomp_lut
