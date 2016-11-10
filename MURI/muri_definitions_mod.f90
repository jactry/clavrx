! $Id$

module muri_definitions_mod

    type muri_input_type
      real :: rfl(6) 
      real :: sol
      real :: sat
      real :: azi  
   end type muri_input_type
   
   type muri_output_type
      real :: aot
   
   end type  muri_output_type
   

end module muri_definitions_mod
