! $Id: constant.f90,v 1.1.1.1 2012/05/18 20:18:40 heidinger Exp $
module CONSTANTS
  implicit none
  integer, parameter, public:: int1 = selected_int_kind(1)
  integer, parameter, public:: int2 = selected_int_kind(3)
  integer, parameter, public:: int4 = selected_int_kind(8)
  integer, parameter, public:: real4 = selected_real_kind(6,37)
  integer, parameter, public:: real8 = selected_real_kind(15,307)
  integer, parameter, public:: ipre = real4
  real (kind=real4), parameter, public:: pi = 3.14159265
  real (kind=real4), parameter, public:: dtor = 3.14159265/180.0
  real(kind=real4), parameter, public:: missing_value_real4 = -99.0
end module CONSTANTS
