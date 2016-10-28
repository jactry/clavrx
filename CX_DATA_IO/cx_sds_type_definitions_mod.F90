module cx_sds_type_definitions_mod

   integer, parameter, public :: MAXNCDIM = 32
   integer, parameter, public :: MAXNCNAM = 128
   
    include 'cx_sds_constants.F90'
   
   type cx_sds_data_type
      integer :: type
      integer :: datasize
      integer :: size
      integer :: nval
      integer :: rank
      integer :: dimsize (MAXNCDIM)
      integer :: utype
      integer :: calbrtd
        
      real ( kind = 8 ) :: calibr(4)
      character, allocatable :: c1values(:)
      integer(kind=1), allocatable :: i1values(:)
      integer(kind=2), allocatable :: i2values(:)
      integer(kind=4), allocatable :: i4values(:)
      real(kind=4), allocatable :: r4values(:)
      real(kind=8), allocatable :: r8values(:)
      
      contains
      procedure :: info=>cx_sds_data_type__info  
      procedure :: transform_to_real
      procedure :: transform_to_dbl
   end type cx_sds_data_type
   
   type cx_att_type
      character ( len = MAXNCNAM) :: name
      type ( cx_sds_data_type ) :: data
      
      
      
   end type cx_att_type
   
   type cx_sds_type
      character( len = MAXNCNAM) :: name
      integer :: nattr
      type ( cx_sds_data_type) :: data
      type(cx_att_type), dimension(:),allocatable :: attr
      contains
      procedure :: info=>cx_sds_type__info 
      procedure :: get_att => cx_sds_type__get_att 
   end type cx_sds_type
   
   
    contains 
      !  -------------------------------------------------------------------------------
      !
      !-------------------------------------------------------------------------------
   subroutine cx_sds_data_type__info (self)
      class(cx_sds_data_type):: self
      
      print*,'data size: ', self.size
      print*,'data rank: ',self.rank
      print*,'dimsize: ',self.dimsize(1:self.rank)
      print*,'Calibr: ',self.calibr
      print*,'data_type: ',self.type
      select case ( self.type)
      case (DFNT_CHAR8)
         print*,'CHAR8 '
      case (DFNT_UCHAR8, DFNT_UINT8, DFNT_INT8)
         print*,'UINT8 '
         print*,'shape data: ',shape(self % i1values)
      case (DFNT_UINT16, DFNT_INT16)  
         print*,'UINT16 '
      case (DFNT_UINT32, DFNT_INT32)  
         print*,'UINT32 '
      case (DFNT_FLOAT32)
         print*,'FLOAT32 '
      case (DFNT_FLOAT64)
         print*,'FLOAT64 '
      end select
      
      
      print*, (allocated( self% i1values))
      print*, (allocated(self% i2values))
      print*, (allocated( self%i4values))
      print*, (allocated( self%r4values))
      print*, (allocated( self%r8values))
      
      
   end subroutine cx_sds_data_type__info
   
   
   !-------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------
   subroutine transform_to_real (self, data_real)
      class(cx_sds_data_type), intent(in) :: self
      real, intent(out) :: data_real(:)
         
      select case ( self.type)
      case (DFNT_CHAR8)
        
      case (DFNT_UCHAR8, DFNT_UINT8, DFNT_INT8)
         data_real = real (self % i1values)
      case (DFNT_UINT16, DFNT_INT16)  
         data_real = real (self % i2values)
      case (DFNT_UINT32, DFNT_INT32)  
         data_real = real (self % i4values)
      case (DFNT_FLOAT32)
         !print*,self % r4values
         data_real = real (self % r4values)
         
      case (DFNT_FLOAT64)
         data_real = real (self % r8values)
      end select
   
   end subroutine transform_to_real  

   !-------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------
   subroutine transform_to_dbl (self, data_dbl)
      class(cx_sds_data_type), intent(in) :: self
      real(kind = 8), intent(out) :: data_dbl(:)
   
      select case ( self.type)
      case (DFNT_CHAR8)
        
      case (DFNT_UCHAR8, DFNT_UINT8, DFNT_INT8)
         data_dbl = real (self % i1values)
      case (DFNT_UINT16, DFNT_INT16)  
         data_dbl = real (self % i2values)
      case (DFNT_UINT32, DFNT_INT32)  
         data_dbl = real (self % i4values)
      case (DFNT_FLOAT32)
         data_dbl = real (self % r4values)
      case (DFNT_FLOAT64)
         data_dbl = real (self % r8values)
      end select
   
   end subroutine transform_to_dbl     
   
   !-------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------      
   subroutine cx_sds_type__info(self)
      class(cx_sds_type):: self
      integer :: i
      print*,self.name
      print*,'info SDS TYPE'
      print*,'number attributes: ', self.nattr
      
      do i=1, self.nattr
         print*,self.attr(i).name
      
      end do
      
   end subroutine cx_sds_type__info
   !-------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------    
   function cx_sds_type__get_att (self, att_name)
      class(cx_sds_type) , target :: self
      character(len=*) :: att_name
      real :: cx_sds_type__get_att(1)
      type ( cx_sds_data_type), pointer :: pd
      integer :: i
      
      do i =1, self.nattr

         if (self.attr(i).name .EQ. trim(att_name)) then
            pd=>self.attr(i).data
            call pd.transform_to_real(cx_sds_type__get_att)
         end if
               
      end do
      
      
   
   end function cx_sds_type__get_att


end module cx_sds_type_definitions_mod
