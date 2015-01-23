! $Id:$

module cx_prd_sfc_mod


   type prd_sfc_main_type
      real, allocatable :: ndvi (:,:)
      real, allocatable :: sst_masked (:,:)
      real, allocatable :: sst_unmasked (:,:)
   end type prd_sfc_main_type

contains

   subroutine populate 
   
   
   end subroutine populate

end module cx_prd_sfc_mod
