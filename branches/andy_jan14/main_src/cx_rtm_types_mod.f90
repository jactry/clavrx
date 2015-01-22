! $Id:$
!
!   defines rtm structure and substructures with allocation/deallocation routines
!
!     rtm
!        % sgrid
!        % ngrid
!

module cx_rtm_types_mod
   
   
  

   integer, parameter :: real4 = selected_real_kind(6,37)
   integer, parameter :: int1 = selected_int_kind(1)
   
   integer, parameter, public:: RTM_NVZEN = 50
   real, parameter, public::  RTM_VZA_BINSIZE = 1./ RTM_NVZEN
   
  
 

end module cx_rtm_types_mod
