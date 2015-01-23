module cx_sfc_tools_mod
   public :: update_sfc

contains

   subroutine update_sfc ( sfc_obj , sat_obj)
      use cx_sfc_mod
      use cx_sat_mod
      
      type ( sfc_main_type ) :: sfc_obj
      type ( sat_main_type ) :: sat_obj
     
		if ( allocated ( sat_obj % chn(2) % ref )) then
      	sfc_obj % ndvi % data = ndvi (sat_obj % chn(2) % ref , sat_obj % chn(1) % ref) 
		end if	
     
   
   end subroutine
   
   elemental real function ndvi ( rfl_1, rfl_2 )
      real , intent(in) :: rfl_1, rfl_2
      ndvi = (rfl_2-rfl_1 ) / ( rfl_2+rfl_1) 
   end function ndvi
   
   
   !====================================================================
   !
   ! Attempt to fix the land classification based on observed ndvi
   !
   ! if the ndvi is high and the land class is not land, this pixel should be land
   ! if the ndvi is low and the land class is land, this pixel should be water
   !
   !====================================================================
   subroutine MODIFY_LAND_CLASS_WITH_NDVI( land_class, ndvi , solzen, bad_mask )
      
      integer , intent(inout) :: land_class (:,:)
      real, intent(in) :: ndvi (:,:)
      real, intent(in) :: solzen (:,:)
      logical , intent(in) :: bad_mask(:,:)
      
      real, parameter:: Ndvi_Land_Threshold = 0.25
      real, parameter:: Ndvi_Water_Threshold = -0.25
      real, parameter:: Solzen_Threshold = 60.0
      integer ::LAND=1
      integer :: SHALLOW_INLAND_WATER = 7
      
      logical, allocatable :: to_apply_land (:,:) , to_apply_water (:,:)
      

    !  if (Sensor%Chan_On_Flag_Default(1) == sym%NO) return
     ! if (Sensor%Chan_On_Flag_Default(2) == sym%NO) return
     ! if (index(Sensor%Sensor_Name,'MODIS') > 0) return
  
       !modis ch2 saturates, need to modify for MODIS
       
       
      
      to_apply_land = .not. Bad_Mask  .and. &
                  Solzen <= Solzen_Threshold .and. &
                  ndvi > Ndvi_Land_Threshold
                                    
      to_apply_water = .not. Bad_Mask .and. &
                  Solzen <= Solzen_Threshold .and. &
                  ndvi < Ndvi_Water_Threshold    .and. &
                   land_class  ==  LAND         
      
      where ( to_apply_land )
         land_class  =  LAND      
      end where  
      
      where ( to_apply_water )
          land_class  = SHALLOW_INLAND_WATER      
      end where     
   
      deallocate ( to_apply_land)
      deallocate ( to_apply_water )
      

   end subroutine MODIFY_LAND_CLASS_WITH_NDVI   
   

end module 
