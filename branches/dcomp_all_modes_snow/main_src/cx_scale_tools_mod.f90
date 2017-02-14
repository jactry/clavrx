! $Id$

module cx_scale_tools_mod
  
   
   use cx_constants_mod, only: &
      int4, real4, int2, int1 &
      , missing_value_int2 &
      , missing_value_int1
      
    implicit none   
   
   private
   integer(kind=int4), parameter :: TWO_BYTE_MAX = 32767, & !(2**15)/2 - 1
                                        TWO_BYTE_MIN = -32767   !-(2**15)/2
                                        
   integer(kind=int4), parameter :: ONE_BYTE_MAX = 127, & !(2**8)/2 - 1
                                        ONE_BYTE_MIN = -127   !-(2**8)/2
    
    
   public :: scale_i2_rank2   
   public :: scale_i1_rank2
                                    
contains
   !
   !
   !
    subroutine scale_i2_rank2 ( data_r4 &
            , unscaled_min &
            , unscaled_max &
            , unscaled_missing &
            , out_i2) 
            
      implicit none
      
      real(kind=real4), intent(in) :: data_r4(:,:)
      
      real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
      integer :: size_arrays(2)
      integer (kind = int2), allocatable, intent(out) ::out_i2(:,:)
      
      
      size_arrays = shape ( data_r4)
      allocate (out_i2(size_arrays(1),size_arrays(2)))
      out_i2 = two_byte_min + min(1.0,max(0.0,(data_r4 - unscaled_min)/(unscaled_max - unscaled_min))) &
         * (TWO_BYTE_MAX - TWO_BYTE_MIN)
           
      !--- set scaled missing values
      where (data_r4 == unscaled_missing)
         out_i2 = missing_value_int2
      end where
   
   end subroutine scale_i2_rank2
   
   subroutine scale_i1_rank2 ( data_r4 &
            , unscaled_min &
            , unscaled_max &
            , unscaled_missing &
            , out_i1) 
            
      implicit none
      
      real(kind=real4), intent(in) :: data_r4(:,:)
      
      real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
      integer :: size_arrays(2)
      integer (kind = int1), allocatable, intent(out) ::out_i1(:,:)
      
      
      size_arrays = shape ( data_r4)
      allocate (out_i1(size_arrays(1),size_arrays(2)))
      out_i1 = one_byte_min + min(1.0,max(0.0,(data_r4 - unscaled_min)/(unscaled_max - unscaled_min))) &
         * (ONE_BYTE_MAX - ONE_BYTE_MIN)
           
      !--- set scaled missing values
      where (data_r4 == unscaled_missing)
         out_i1 = missing_value_int1
      end where
      
   
   end subroutine scale_i1_rank2
   
   
end module cx_scale_tools_mod
