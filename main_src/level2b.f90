!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: level2b.f90 (src)
!       LEVEL2B_ROUTINES (program)
!
! PURPOSE: Routines for creating, writing and closing  L2b output files
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
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
! REVISON HISTORY:
!    Creation Date May 2009
!--------------------------------------------------------------------------------------
module LEVEL2B_ROUTINES
   use CONSTANTS, only:  &
   int4 &
   , real4 &
   , int2 &
   , int1 &
   , sym &
   , EXE_PROMPT &
   , missing_value_int1 &
   , missing_value_int2 &
   , missing_value_int4 &
   , Missing_Value_Real4 &
   , No_Attribute_Missing_Value
   

   use FILE_UTILITY,only: &
    get_lun

   implicit none
   include 'hdf.f90' 
   private
   
   public:: INIT_RANDOM_SEED
  

 

!interface DEFINE_SDS
!    module procedure  &
!        DEFINE_SDS_RANK1,  &
!        DEFINE_SDS_RANK2,  &
!        DEFINE_SDS_RANK3
!end interface

 type, public :: Sds_Struct
  character(len=100):: Variable_Name
  integer(kind=int4):: Data_Type
  integer(kind=int4):: Rank
  integer(kind=int1):: Scaling_Type
  character(len=100):: Units
  character(len=300):: Long_Name
  character(len=100):: Standard_Name
  real(kind=real4):: Unscaled_Min
  real(kind=real4):: Unscaled_Max
  real(kind=real4):: Unscaled_Missing
  real(kind=real4):: Scale_Factor
  real(kind=real4):: Add_Offset
  real(kind=real4), dimension(2):: Actual_Range
  real(kind=real4):: Fill_Value
  integer(kind=int4), dimension(2):: Valid_Range
  integer(kind=int4):: Status_Flag
  integer(kind=int4):: Id_Input
  integer(kind=int4):: Id_Output
  integer(kind=int4):: Compression_Type
  integer(kind=int4):: Num_Attrs
 end type Sds_Struct

 character(len=30), parameter, private :: MOD_PROMPT = " PATMOSx_LEVEL2B_ROUTINES: "
 character(len=18), private, parameter:: Coordinates_String = "longitude latitude"

 real, dimension(:,:), allocatable, public, save:: Lat_Input
 real, dimension(:,:), allocatable, public, save:: Lon_Input
 logical, dimension(:,:), allocatable, public, save:: Gap_Pixel_Mask_Input
 real, dimension(:,:), allocatable, public, save:: Lat_Output
 real, dimension(:,:), allocatable, public, save:: Lon_Output
 integer(kind=int4), dimension(:,:), allocatable, public, save:: Ielem_Output
 integer(kind=int4), dimension(:,:), allocatable, public, save:: Iline_Output

 !--- no compression
 !integer, parameter, private:: Comp_Type = 0
 !integer, dimension(2), parameter, private:: Comp_Prm = /(0,0)/
 !--- gzip compression
 integer, parameter, private:: Comp_Type = 4
 integer, dimension(4), parameter, private:: Comp_Prm = (/6,0,0,0/)
 !--- szip compression
 !integer, parameter, private:: Comp_Type = 5
 !integer, dimension(2), parameter, private:: Comp_Prm = /(32,0)/


 contains







!------------------------------------------------------------------------------------------------
! SUBROUTINE Name: READ_SDS_RANK3
!
! Function:
!     Reads 3-dimensional SDSs from Level 2b files and reshape to 2 dimensions
!     (remove time coordinate). 
!-----------------------------------------------------------------------------------------------
  subroutine READ_SDS_RANK3(Sd_Id,     &
                            Scaled_Sds_Data,   &
                            Sds)

    integer, intent(in):: Sd_Id
    type(Sds_Struct), intent(inout):: Sds
    real(kind=real4), dimension(:,:,:), intent(out):: Scaled_Sds_Data
    integer, dimension(3):: Sds_Dims
    integer, dimension(3):: Sds_Start
    integer, dimension(3):: Sds_Stride
    integer, dimension(3):: Sds_Edges
    integer(kind=int1), dimension(:,:,:), allocatable:: Temp_I1
    integer(kind=int2), dimension(:,:,:), allocatable:: Temp_I2
    integer(kind=int4), dimension(:,:,:), allocatable:: Temp_I4
    real(kind=real4), dimension(:,:,:), allocatable:: Temp_R4
    integer:: Num_Attrs
    integer:: sfn2index
    integer:: sfselect
    
    integer:: sfginfo
    integer:: sfrdata
    integer:: sfendacc
    integer:: Istatus
    character(100):: Sds_Name_Temp

    Istatus = 0

    !--- open Sds for reading
    Sds%Id_Input = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds%Variable_Name)))
    if ( Sds%Id_Input < 0 ) then
       Sds%Data_Type = -999
       return
    endif
    Istatus = sfginfo(Sds%Id_Input, Sds_Name_Temp, Sds%rank, Sds_Dims, Sds%Data_Type, Num_Attrs) + Istatus
    Sds_Start = (/ 0, 0, 0 /)
    Sds_Stride = (/ 1, 1, 1 /)
    Sds_Edges = Sds_Dims

    !--- read sds scaling attributes
    call READ_SCALING_ATTRIBUTES(Sds,Istatus)

    !--- allocate arrays for holding data, read data and store in output array
    if (Sds%Data_Type == DFNT_INT8) then
       allocate(Temp_I1(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_I1) + Istatus
       Scaled_Sds_Data = real(Temp_I1)
    elseif (Sds%Data_Type == DFNT_INT16) then
       allocate(Temp_I2(Sds_Dims(1), Sds_Dims(2),Sds_Dims(3)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_I2) + Istatus
       Scaled_Sds_Data = real(Temp_I2)
    elseif (Sds%Data_Type == DFNT_INT32) then
       allocate(Temp_I4(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_I4) + Istatus
       Scaled_Sds_Data = real(Temp_I4)
    elseif (Sds%Data_Type == DFNT_FLOAT32) then
       allocate(Temp_R4(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
       Istatus = sfrdata(Sds%Id_Input, Sds_Start, Sds_Stride, Sds_Edges, Temp_R4) + Istatus
       Scaled_Sds_Data = real(Temp_R4)
    else
       print *, "Possibly fatal error at location 2:"
       print *, "data type was ", Sds%Data_Type
       print *, "sds complete value is "
       print *, Sds
       print *, MOD_PROMPT, "attempt to read unsupported data type, stopping"
       Scaled_Sds_Data = Missing_Value_Real4
       Sds%Data_Type = -999
       return
    endif

    !---deallocate temp arrays
    if (allocated(Temp_I1)) deallocate(Temp_I1)
    if (allocated(Temp_I2)) deallocate(Temp_I2)
    if (allocated(Temp_I4)) deallocate(Temp_I4)
    if (allocated(Temp_R4)) deallocate(Temp_R4)

    !--- close Sds
    Istatus = sfendacc(Sds%Id_Input) + Istatus

   end subroutine READ_SDS_RANK3
!------------------------------------------------------------------------------------------------
! SUBROUTINE Name: UNSCALE_SDS_RANK1
!
! Function:
!     Unscales scaled SDS data from level2b files
!-----------------------------------------------------------------------------------------------
   subroutine UNSCALE_SDS_RANK1(Sds, Scaled_Data, Unscaled_Data)
      type(Sds_Struct), intent(in):: Sds
      real(kind=real4), dimension(:), intent(in):: Scaled_Data
      real(kind=real4), dimension(:), intent(out):: Unscaled_Data

      !--- unscale Sds
      
      if (Sds%Scaling_Type /= sym%NO_SCALING) then

        !---- linear
        if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
           Unscaled_Data = Sds%Scale_Factor * real(Scaled_Data) + Sds%Add_Offset
        endif

        where (Scaled_Data == Sds%Fill_Value)
         Unscaled_Data = Missing_Value_Real4
        endwhere

      else

        Unscaled_Data = Scaled_Data

      endif

   end subroutine UNSCALE_SDS_RANK1
!------------------------------------------------------------------------------------------------
! SUBROUTINE Name: UNSCALE_SDS_RANK2
!
! Function:
!     Unscales scaled SDS data from level2b files
!-----------------------------------------------------------------------------------------------
   subroutine UNSCALE_SDS_RANK2(Sds, Scaled_Data, Unscaled_Data)
      type(Sds_Struct), intent(in):: Sds
      real(kind=real4), dimension(:,:), intent(in):: Scaled_Data
      real(kind=real4), dimension(:,:), intent(out):: Unscaled_Data

      !--- unscale Sds

      if (Sds%Scaling_Type /= sym%NO_SCALING) then

        !---- linear
        if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
           Unscaled_Data = Sds%Scale_Factor * real(Scaled_Data) + Sds%Add_Offset
        endif

        where (Scaled_Data == Sds%Fill_Value)
         Unscaled_Data = Missing_Value_Real4
        endwhere

      else

        Unscaled_Data = Scaled_Data

      endif

   end subroutine UNSCALE_SDS_RANK2

!------------------------------------------------------------------------------------------------
! SUBROUTINE Name: UNSCALE_SDS_RANK3
!
! Function:
!     Unscales scaled SDS data from level2b files
!-----------------------------------------------------------------------------------------------
   subroutine UNSCALE_SDS_RANK3(Sds, Scaled_Data, Unscaled_Data)
      type(Sds_Struct), intent(in):: Sds
      real(kind=real4), dimension(:,:,:), intent(in):: Scaled_Data
      real(kind=real4), dimension(:,:,:), intent(out):: Unscaled_Data

      !--- unscale Sds
      if (Sds%Scaling_Type /= sym%NO_SCALING) then

        !---- linear
        if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
           Unscaled_Data = Sds%Scale_Factor * real(Scaled_Data) + Sds%Add_Offset
        endif

        !--- handle missing
        where (Scaled_Data == Sds%Fill_Value)
         Unscaled_Data = Sds%Unscaled_Missing
        endwhere
        
      else
        Unscaled_Data = Scaled_Data
      endif

   end subroutine UNSCALE_SDS_RANK3
   !------------------------------------------------------------------------------------------------
   ! SUBROUTINE Name: SCALE_SDS_RANK1
   !
   ! Function:
   !     Scales unscaled SDS data from level2b files
   !-----------------------------------------------------------------------------------------------
   subroutine SCALE_SDS_RANK1(Sds, Unscaled_Data, Scaled_Data)
      type(Sds_Struct), intent(in):: Sds
      real(kind=real4), dimension(:), intent(in):: Unscaled_Data
      real(kind=real4), dimension(:), intent(out):: Scaled_Data
      
      !--- scale Sds
      if (Sds%Scaling_Type /= sym%NO_SCALING) then
      
        !---- linear
        if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
           Scaled_Data = (Unscaled_Data - Sds%Add_Offset)/(Sds%Scale_Factor)
           where (Scaled_Data < Sds%Valid_Range(1))
             Scaled_Data = Sds%Valid_Range(1)
           endwhere
           where (Scaled_Data > Sds%Valid_Range(2))
             Scaled_Data = Sds%Valid_Range(2)
           endwhere
        endif

        !--- set scaled missing values
        where (Unscaled_Data == Sds%Unscaled_Missing)
         Scaled_Data = Sds%Fill_Value
        end where

      else
         Scaled_Data = Unscaled_Data
      end if

   end subroutine SCALE_SDS_RANK1
!------------------------------------------------------------------------------------------------
! SUBROUTINE Name: SCALE_SDS_RANK2
!
! Function:
!     Scales unscaled SDS data from level2b files
!-----------------------------------------------------------------------------------------------
   subroutine SCALE_SDS_RANK2(Sds, Unscaled_Data, Scaled_Data)
      type(Sds_Struct), intent(in):: Sds
      real(kind=real4), dimension(:,:), intent(in):: Unscaled_Data
      real(kind=real4), dimension(:,:), intent(out):: Scaled_Data

      !--- scale Sds
      if (Sds%Scaling_Type /= sym%NO_SCALING) then

        !---- linear
        if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
           Scaled_Data = (Unscaled_Data - Sds%Add_Offset)/(Sds%Scale_Factor)
           where (Scaled_Data < Sds%Valid_Range(1))
             Scaled_Data = Sds%Valid_Range(1)
           endwhere
           where (Scaled_Data > Sds%Valid_Range(2))
             Scaled_Data = Sds%Valid_Range(2)
           endwhere
        endif

        !--- set scaled missing values
        where (Unscaled_Data == Sds%Unscaled_Missing)
         Scaled_Data = Sds%Fill_Value
        endwhere

     else
       Scaled_Data = Unscaled_Data
     endif

   end subroutine SCALE_SDS_RANK2

!------------------------------------------------------------------------------------------------
! SUBROUTINE Name: SCALE_SDS_RANK3
!
! Function:
!     Scales unscaled SDS data from level2b files
!-----------------------------------------------------------------------------------------------
   subroutine SCALE_SDS_RANK3(Sds, Unscaled_Data, Scaled_Data)
      type(Sds_Struct), intent(in):: Sds
      real(kind=real4), dimension(:,:,:), intent(in):: Unscaled_Data
      real(kind=real4), dimension(:,:,:), intent(out):: Scaled_Data

      !--- scale Sds
      if (Sds%Scaling_Type /= sym%NO_SCALING) then

        !---- linear
        if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
           Scaled_Data = (Unscaled_Data - Sds%Add_Offset)/(Sds%Scale_Factor)
           where (Scaled_Data < Sds%Valid_Range(1))
             Scaled_Data = Sds%Valid_Range(1)
           endwhere
           where (Scaled_Data > Sds%Valid_Range(2))
             Scaled_Data = Sds%Valid_Range(2)
           endwhere
        endif

        !--- set scaled missing values
        where (Unscaled_Data == Sds%Unscaled_Missing)
         Scaled_Data = Sds%Fill_Value
        endwhere

     else
       Scaled_Data = Unscaled_Data
     endif

   end subroutine SCALE_SDS_RANK3


!-----------------------------------------------------------------------------
! SUBROUTINE Name: WRITE_SDS_RANK1
!
! Function:
!     Write 1 dimensional data to a HDF file
!-----------------------------------------------------------------------------
  subroutine WRITE_SDS_RANK1(Scaled_Sds_Data,   &
                             Sds)

    type(Sds_Struct), intent(in):: Sds
    real(kind=real4), dimension(:), intent(in):: Scaled_Sds_Data
    integer, dimension(1):: Sds_Dims
    integer, dimension(1):: Sds_Start
    integer, dimension(1):: Sds_Stride
    integer, dimension(1):: Sds_Edges
    integer(kind=int1), dimension(:), allocatable:: Temp_I1
    integer(kind=int2), dimension(:), allocatable:: Temp_I2
    integer(kind=int4), dimension(:), allocatable:: Temp_I4
    real(kind=real4), dimension(:), allocatable:: Temp_R4
    integer(kind=int4):: Istatus
    integer:: sfwdata
    
    
    print*,'ffgggg 1d'

    Sds_Dims(1) = size(Scaled_Sds_Data,1)
    Sds_Start = (/ 0 /)
    Sds_Stride = (/ 1 /)
    Sds_Edges = Sds_Dims


    Istatus = 0

    Sds_Dims(1) = size(Scaled_Sds_Data,1)

    if (Sds%Data_Type == DFNT_INT8) then
       allocate(Temp_I1(Sds_Dims(1)))
       Temp_I1 = int(Scaled_Sds_Data,kind=int1)
       Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I1) + Istatus
    elseif (Sds%Data_Type == DFNT_INT16) then
       allocate(Temp_I2(Sds_Dims(1)))
       Temp_I2 = int(Scaled_Sds_Data,kind=int2)
       Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I2) + Istatus
    elseif (Sds%Data_Type == DFNT_INT32) then
       allocate(Temp_I4(Sds_Dims(1)))
       Temp_I4 = int(Scaled_Sds_Data,kind=int4)
       Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I4) + Istatus
    elseif (Sds%Data_Type == DFNT_FLOAT32) then
       allocate(Temp_R4(Sds_Dims(1)))
       Temp_R4 = real(Scaled_Sds_Data,kind=real4)
       Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_R4) + Istatus
    endif


   !---deallocate temp arrays
   if (allocated(Temp_I1)) deallocate(Temp_I1)
   if (allocated(Temp_I2)) deallocate(Temp_I2)
   if (allocated(Temp_I4)) deallocate(Temp_I4)
   if (allocated(Temp_R4)) deallocate(Temp_R4)

   end subroutine WRITE_SDS_RANK1

!------------------------------------------------------------------------------------------------
! SUBROUTINE Name: WRITE_SDS_RANK2
!
! Function:
!     Write 2 dimensional data to a HDF file
!-----------------------------------------------------------------------------------------------
    subroutine WRITE_SDS_RANK2(Scaled_Sds_Data,Sds)

     type(Sds_Struct), intent(in):: Sds
     real(kind=real4), dimension(:,:), intent(in):: Scaled_Sds_Data
     integer, dimension(2):: Sds_Dims
     integer, dimension(2):: Sds_Start
     integer, dimension(2):: Sds_Stride
     integer, dimension(2):: Sds_Edges
     integer(kind=int1), dimension(:,:), allocatable:: Temp_I1
     integer(kind=int2), dimension(:,:), allocatable:: Temp_I2
     integer(kind=int4), dimension(:,:), allocatable:: Temp_I4
     real(kind=real4), dimension(:,:), allocatable:: Temp_R4
     integer(kind=int4):: Istatus
     integer:: sfwdata
       print*,'ffgggg 2d'
     Sds_Dims(1) = size(Scaled_Sds_Data,1)
     Sds_Dims(2) = size(Scaled_Sds_Data,2)
     Sds_Start = (/ 0, 0 /)
     Sds_Stride = (/ 1, 1 /)
     Sds_Edges = Sds_Dims

     Istatus = 0
 print*,'kk', Sds%Data_Type, Sds_Dims(1),Sds_Dims(2)
     if (Sds%Data_Type == DFNT_INT8) then
        allocate(Temp_I1(Sds_Dims(1),Sds_Dims(2)))
        Temp_I1 = int(Scaled_Sds_Data,kind=int1)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I1) + Istatus
     elseif (Sds%Data_Type == DFNT_INT16) then
        allocate(Temp_I2(Sds_Dims(1),Sds_Dims(2)))
        Temp_I2 = int(Scaled_Sds_Data,kind=int2)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I2) + Istatus
     elseif (Sds%Data_Type == DFNT_INT32) then
        allocate(Temp_I4(Sds_Dims(1),Sds_Dims(2)))
        Temp_I4 = int(Scaled_Sds_Data,kind=int4)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I4) + Istatus
     elseif (Sds%Data_Type == DFNT_FLOAT32) then
         print*,'w1'
        allocate(Temp_R4(Sds_Dims(1),Sds_Dims(2)))
        print*,'w2',size(temp_r4)
        Temp_R4 = real(Scaled_Sds_Data,kind=real4)
        print*,'w3'
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_R4) + Istatus
        print*,'w4'
     endif
 print*,'fff'
   !---deallocate temp arrays
   if (allocated(Temp_I1)) deallocate(Temp_I1)
   if (allocated(Temp_I2)) deallocate(Temp_I2)
   if (allocated(Temp_I4)) deallocate(Temp_I4)
   if (allocated(Temp_R4)) deallocate(Temp_R4)

end subroutine WRITE_SDS_RANK2

!------------------------------------------------------------------------------------------------
! SUBROUTINE Name: WRITE_SDS_RANK3
!
! Function:
!     Write 3 dimensional data to a HDF file
!-----------------------------------------------------------------------------------------------
    subroutine WRITE_SDS_RANK3(Scaled_Sds_Data,Sds)

     type(Sds_Struct), intent(in):: Sds
     real(kind=real4), dimension(:,:,:), intent(in):: Scaled_Sds_Data
     integer, dimension(3):: Sds_Dims
     integer, dimension(3):: Sds_Start
     integer, dimension(3):: Sds_Stride
     integer, dimension(3):: Sds_Edges
     integer(kind=int1), dimension(:,:,:), allocatable:: Temp_I1
     integer(kind=int2), dimension(:,:,:), allocatable:: Temp_I2
     integer(kind=int4), dimension(:,:,:), allocatable:: Temp_I4
     real(kind=real4), dimension(:,:,:), allocatable:: Temp_R4
     integer(kind=int4):: Istatus
     integer:: sfwdata
       print*,'ffgggg 3d ',sds % variable_name
     Sds_Dims(1) = size(Scaled_Sds_Data,1)
     Sds_Dims(2) = size(Scaled_Sds_Data,2)
     Sds_Dims(3) = size(Scaled_Sds_Data,3)
     Sds_Start = (/ 0, 0, 0 /)
     Sds_Stride = (/ 1, 1, 1 /)
     Sds_Edges = Sds_Dims

     Istatus = 0
      print*,'kk', Sds%Data_Type
     if (Sds%Data_Type == DFNT_INT8) then
        allocate(Temp_I1(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
        Temp_I1 = int(Scaled_Sds_Data,kind=int1)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I1) + Istatus
     elseif (Sds%Data_Type == DFNT_INT16) then
        allocate(Temp_I2(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
        Temp_I2 = int(Scaled_Sds_Data,kind=int2)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I2) + Istatus
     elseif (Sds%Data_Type == DFNT_INT32) then
        allocate(Temp_I4(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
        Temp_I4 = int(Scaled_Sds_Data,kind=int4)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_I4) + Istatus
     elseif (Sds%Data_Type == DFNT_FLOAT32) then
        allocate(Temp_R4(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
        Temp_R4 = real(Scaled_Sds_Data,kind=real4)
        Istatus = sfwdata(Sds%Id_Output, Sds_Start, Sds_Stride, Sds_Edges, Temp_R4) + Istatus
     endif
     
   !---deallocate temp arrays
   if (allocated(Temp_I1)) deallocate(Temp_I1)
   if (allocated(Temp_I2)) deallocate(Temp_I2)
   if (allocated(Temp_I4)) deallocate(Temp_I4)
   if (allocated(Temp_R4)) deallocate(Temp_R4)

end subroutine WRITE_SDS_RANK3

!------------------------------------------------------------------------------------------------
! SUBROUTINE Name: COPY_GLOBAL_ATTRIBUTES
!
! Function:
!    Copies the global attributes from one file and puts them in a different file
!-----------------------------------------------------------------------------------------------

subroutine COPY_GLOBAL_ATTRIBUTES(Sd_Id_Input,Sd_Id_Output)
   
   integer(kind=int4):: Sd_Id_Input
   integer(kind=int4):: Sd_Id_Output
   integer:: Istatus
   integer(kind=int4):: num_global_attrs
   integer(kind=int4):: Num_Sds_Input
   integer(kind=int4):: iattr
   character(len=100):: Attr_Name
   integer(kind=int4):: Data_Type
   integer(kind=int4):: count
   integer(kind=int4):: Attr_Buffer_i4
   integer(kind=int2):: Attr_Buffer_i2
   integer(kind=int1):: Attr_Buffer_i1
   real(kind=real4):: Attr_Buffer_r4
   character(len=500):: Attr_Buffer_char

   integer:: sffinfo
   integer:: sfgainfo
   integer:: sfrnatt
   integer:: sfsnatt
   
   
   Istatus = sffinfo(Sd_Id_Input,Num_Sds_Input,num_global_attrs)

   do iattr = 0, num_global_attrs-1

       Istatus = sfgainfo(Sd_Id_Input,iattr,Attr_Name,Data_Type,count)

       !--- skip certain level2 attributes in level2b
       if (Attr_Name == "FILENAME" .or. &
           Attr_Name == "L1B" .or. &
           Attr_Name == "NUMBER_OF_SCANS_LEVEL1B" .or. &
           Attr_Name == "NUMBER_OF_SCANS_LEVEL2" .or. &
           Attr_Name == "START_YEAR" .or.  &
           Attr_Name == "END_YEAR" .or.  &
           Attr_Name == "START_DAY" .or.  &
           Attr_Name == "END_DAY" .or.  &
           Attr_Name == "START_TIME" .or.  &
           Attr_Name == "END_TIME") then
           cycle
        endif

       if (Data_Type == DFNT_INT8) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_i1)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,Attr_Buffer_i1)
       endif

       if (Data_Type == DFNT_INT16) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_i2)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,Attr_Buffer_i2)
       endif

       if (Data_Type == DFNT_INT32) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_i4)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,Attr_Buffer_i4)
       endif

       if (Data_Type == DFNT_FLOAT32) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_r4)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,Attr_Buffer_r4)
       endif

       if (Data_Type == DFNT_CHAR) then
           Istatus = sfrnatt(Sd_Id_Input,iattr,Attr_Buffer_char)
           Istatus = sfsnatt(Sd_Id_Output,trim(Attr_Name),Data_Type,count,trim(Attr_Buffer_char))
       endif
   enddo

   end subroutine 

!------------------------------------------------------------------------------------------------
! SUBROUTINE Name: SUBSET_LEVEL2b
!
! Function:
!    Reads in a level-2b file and then write out a subset of it
!
! Notes:
!    This routine assumes there are only 1d and 2d arrays
!    This routine assumes the only 1d arrays are latitude and longitude
!   
!-----------------------------------------------------------------------------------------------
  subroutine SUBSET_LEVEL2b(File_Input, &
                            Path_Input, & 
                            File_Output, &
                            Path_Output, & 
                            Lon_West_Out_Temp, &
                            Lon_East_Out_Temp, &
                            Dlon_Out, &
                            Lat_South_Out_Temp, &
                            Lat_North_Out_Temp, &
                            DLat_Out, &
                            Data_Description_String, &
                            Sds_Output_Names)

  character(len=*), intent(in):: File_Input
  character(len=*), intent(in):: Path_Input
  character(len=*), intent(in):: File_Output
  character(len=*), intent(in):: Path_Output
  character(len=*), intent(in):: Data_Description_String
  character(len=*), dimension(:), intent(in):: Sds_Output_Names
  real(kind=real4), intent(in):: Lon_West_Out_Temp
  real(kind=real4), intent(in):: Lon_East_Out_Temp
  real(kind=real4), intent(in):: Dlon_Out
  real(kind=real4), intent(in):: Lat_South_Out_Temp
  real(kind=real4), intent(in):: Lat_North_Out_Temp
  real(kind=real4), intent(in):: DLat_Out

  real(kind=real4) :: Lon_West_Out
  real(kind=real4) :: Lon_East_Out
  real(kind=real4) :: Lat_South_Out
  real(kind=real4) :: Lat_North_Out

  integer(kind=int4):: Sd_Id_Input
  integer(kind=int4):: Sd_Id_Output

  type(Sds_Struct) :: Sds_Input
  type(Sds_Struct) :: Sds_Output

  integer(kind=int4) :: Nlon_Out
  integer(kind=int4) :: NLat_Out
  real(kind=real4):: Lon_West_In
  real(kind=real4):: Lon_East_In
  real(kind=real4):: Dlon_In
  integer(kind=int4):: Nlon_In
  real(kind=real4):: Lat_South_In
  real(kind=real4):: Lat_North_In
  real(kind=real4):: Dlat_In
  integer(kind=int4):: Nlat_In

  integer(kind=int4):: Lon_Stride
  integer(kind=int4):: Lat_Stride
  integer(kind=int4):: Lon_Start
  integer(kind=int4):: Lat_Start

  real(kind=real4), dimension(:), allocatable:: Scaled_Sds_Longitude_Input
  real(kind=real4), dimension(:), allocatable:: Scaled_Sds_Latitude_Input
  real(kind=real4), dimension(:,:), allocatable:: Scaled_Sds_Data_Input
  real(kind=real4), dimension(:), allocatable:: Scaled_Sds_Longitude_Output
  real(kind=real4), dimension(:), allocatable:: Scaled_Sds_Latitude_Output
  real(kind=real4), dimension(:,:), allocatable:: Scaled_Sds_Data_Output

  integer:: sfstart, sfend, sfscatt, sfsnatt, sfrnatt,  &
            sfginfo, sfselect, sffattr, sffinfo
  integer:: Isds, Isds_Inner, Num_Sds_Input, Num_Global_Attrs, Num_Sds_Output
  integer:: Istatus
  integer:: Ifound
  integer(kind=int4), dimension(2):: Sds_Dims

  

  end subroutine SUBSET_LEVEL2b
!----------------------------------------------------------------------
! An example of how to reset a random seed based on system time
!
! taken from:
! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
!----------------------------------------------------------------------
  SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
  END SUBROUTINE init_random_seed
!--------------------------------------------------------------------
! File Search - similar to IDL
!--------------------------------------------------------------------
  subroutine FILE_SEARCH(Path,Search_String,Num_Files,Files)
        character(len=*), intent(in):: Path
        character(len=*), intent(in):: Search_String
        character(len=1020), dimension(:), allocatable, intent(out):: Files
        integer(kind=int4), intent(out):: Num_Files
        integer(kind=int4):: Lun
        integer(kind=int4):: Idx

        Lun = GET_LUN()

        !-- first determine number of files
!       call system('ls -1 -p '//trim(Path)//' | grep level2b.hdf | wc -l > subset_file_list')
        call system('ls -1 -p '//trim(Path)//' | grep '//trim(Search_String)//' | wc -l > subset_file_list')

        !-- next, list the files
!       call system('ls -1 -p '//trim(Path)//' | grep level2b.hdf >> subset_file_list')
        call system('ls -1 -p '//trim(Path)//' | grep '//trim(Search_String)//' >> subset_file_list')

        open (unit=Lun, file = 'subset_file_list', status="old",action="read")

        read(unit=Lun,fmt=*) Num_Files

        if (Num_Files == 0) then
           print *, "FILE_SEARCH found no files with ",trim(Search_String)," in  ", trim(Path)
           return
        endif

        allocate(Files(Num_Files))

        do Idx = 1, Num_Files
          read(unit=Lun,fmt=*) Files(Idx)
        enddo 

        close(unit=Lun) 
        
  end subroutine
!---------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------
subroutine COMPUTE_WMO_ID_KNOWING_SENSOR_NAME(Sensor_Name,WMO_Id)

character(len=*), intent(in):: Sensor_Name
integer(kind=int2), intent(out):: WMO_Id

select case (trim(Sensor_Name))


case("zen")
  WMO_Id = 0
case("metop-02")
  WMO_Id = 4
case("metop-01")
  WMO_Id = 3
case("metop-03")
  WMO_Id = 5
case("met8")
  WMO_Id = 55
case("met9")
  WMO_Id = 56
case("met10")
  WMO_Id = 57
case("met11")
  WMO_Id = 70
case("mtsat-1r")
  WMO_Id = 171 
case("mtsat-2")
  WMO_Id = 172 
case("noaa-06")
  WMO_Id = 706
case("noaa-07")
  WMO_Id = 707
case("noaa-05")
  WMO_Id = 708
case("noaa-08")
  WMO_Id = 200
case("noaa-09")
  WMO_Id = 201
case("noaa-10")
  WMO_Id = 202
case("noaa-11")
  WMO_Id = 203
case("noaa-12")
  WMO_Id = 204
case("noaa-14")
  WMO_Id = 205
case("noaa-15")
  WMO_Id = 206
case("noaa-16")
  WMO_Id = 207
case("noaa-17")
  WMO_Id = 208
case("noaa-18")
  WMO_Id = 209
case("noaa-19")
  WMO_Id = 223
case("npp")
  WMO_Id = 224
case("iff_npp")
  WMO_Id = 224
case("npp-nasa")
  WMO_Id = 224
case("goes-08")
  WMO_Id = 252
case("goes-09")
  WMO_Id = 253
case("goes-10")
  WMO_Id = 254
case("goes-11")
  WMO_Id = 255
case("goes-12")
  WMO_Id = 256
case("goes-13")
  WMO_Id = 257
case("goes-14")
  WMO_Id = 258
case("goes-15")
  WMO_Id = 259
case("terra-modis")
  WMO_Id = 783
case("aqua-modis")
  WMO_Id = 784
case("coms")
  WMO_Id = 810

case default

print *, "Unknown Sensor Name, Can not assign Sc_Id"
stop

end select

end subroutine COMPUTE_WMO_ID_KNOWING_SENSOR_NAME

!--------------------------------------------------------------------
!  Write the flag_values and flag_meanings attributes for
!  relevant variables
!--------------------------------------------------------------------
subroutine WRITE_FLAG_ATTRIBUTES(Sds_Id,Sds_Name,IStatus_Sum)
   integer, intent(in):: Sds_Id
   integer, intent(inout):: Istatus_Sum
   character(len=*), intent(in):: Sds_Name
   character(len=500):: Temp_Name
   integer(kind=int1), dimension(2):: yes_no_flag_values
   integer(kind=int1), dimension(4):: qf_flag_values

   integer:: sfsnatt
   integer:: sfscatt

   yes_no_flag_values = int((/0,1/),kind=int1)
   qf_flag_values = int((/0,1,2,3/),kind=int1)

     if (Sds_Name == "acha_quality") then
       Temp_Name = "Processed "// &
                   "valid_Tc_retrieval "// &
                   "valid_ec_retrieval "// &
                   "valid_beta_retrieval "// &
                   "degraded_Tc_retrieval "// &
                   "degraded_ec_retrieval "// &
                   "degraded_beta_retrieval "
       Istatus_Sum = sfsnatt(Sds_Id, "flag_masks", DFNT_INT8, 7, int((/1,2,4,8,16,32,64/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum
     end if

     if (Sds_Name == "acha_info") then
       Temp_Name = "Cloud_Height_Attempted "// &
                   "Bias_Correction_Employed "// &
                   "Ice_Cloud_Retrieval "// &
                   "Local_Radiative_Center_Processing_Used "// &
                   "Multi_Layer_Retrieval "// &
                   "Lower_Cloud_Interpolation_Used "// &
                   "Boundary_Layer_Inversion_Assumed "
       Istatus_Sum = sfsnatt(Sds_Id, "flag_masks", DFNT_INT8, 7, int((/1,2,4,8,16,32,64/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum
     end if

     if (Sds_Name == "dcomp_quality") then
       Temp_Name = "Processed "// &
                   "valid_COD_retrieval "// &
                   "valid_REF_retrieval "// &
                   "degraded_COD_retrieval "// &
                   "degraded_REF_retrieval "// &
                   "convergency "// &
                   "glint "
       Istatus_Sum = sfsnatt(Sds_Id, "flag_masks", DFNT_INT8, 7, int((/1,2,4,8,16,32,64/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum
     end if

     if (Sds_Name == "dcomp_info") then
       Temp_Name = "land_or_sea_mask "// &
                   "day_or_night_mask "// &
                   "twilight_@65-82_solar_zenith "// &
                   "snow "// &
                   "sea_ice "// &
                   "liquid_or_ice_phase "// &
                   "thick_cloud "// &
                   "thin_cloud "
       Istatus_Sum = sfsnatt(Sds_Id, "flag_masks", DFNT_INT16, 8, int((/1,2,4,8,16,32,64,128/),kind=int2)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum
     end if

     if (Sds_Name == "cld_opd_dcomp_qf") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, qf_flag_values) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                                " failed "// &
                                " low_quality "// &
                                " high_quality ") + Istatus_Sum
     end if

     if (Sds_Name == "cld_reff_dcomp_qf") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, qf_flag_values) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                                " failed "// &
                                " low_quality "// &
                                " high_quality ") + Istatus_Sum
     end if

     if (Sds_Name == "cld_temp_acha_qf") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, qf_flag_values) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                                " failed "// &
                                " low_quality "// &
                                " high_quality ") + Istatus_Sum
     end if

     if (Sds_Name == "cld_emiss_acha_qf") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, qf_flag_values) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                                " failed "// &
                                " low_quality "// &
                                " high_quality ") + Istatus_Sum
     end if

     if (Sds_Name == "cld_beta_acha_qf") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, qf_flag_values) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                                " failed "// &
                                " low_quality "// &
                                " high_quality ") + Istatus_Sum
     end if

     if (Sds_Name == "cloud_mask") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, qf_flag_values) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 47, "clear "// &
                                " probably_clear "// &
                                " probably_cloudy "// &
                                " cloudy ") + Istatus_Sum
     end if

     if (Sds_Name == "cloud_type") then
       Temp_Name = "clear "// &
                   "probably_clear "// &
                   "fog "// &
                   "water "// &
                   "supercooled_water "// &
                   "mixed "// &
                   "opaque_ice "// &
                   "cirrus "// &
                   "overlapping "// &
                   "overshooting "// &
                   "unknown "// &
                   "dust "// &
                   "smoke "
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 13, int((/0,1,2,3,4,5,6,7,8,9,10,11,12/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum
     end if

     if (Sds_Name == "surface_type") then
       Temp_Name = "water "// &
                   "evergreen_needle "// &
                   "evergreen_broad "// &
                   "deciduous_needle "// &
                   "deciduous_broad "// &
                   "mixed_forest "// &
                   "woodlands "// &
                   "wooded_grass "// &
                   "closed_shrubs "// &
                   "open_shrubs "// &
                   "grasses "// &
                   "croplands "// &
                   "bare "// &
                   "urban "
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 14, int((/0,1,2,3,4,5,6,7,8,9,10,11,12,13/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum
     end if

     if (Sds_Name == "land_class") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 8, int((/0,1,2,3,4,5,6,7/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 109, "ocean "// &
                                " land "// &
                                " coastline "// &
                                " shallow_inland_water "// &
                                " ephemeral_water "// &
                                " deep_inland_water "// &
                                " moderate_ocean "// &
                                " deep_ocean ") + Istatus_Sum
     end if

     if (Sds_Name == "snow_class") then
       Istatus_Sum = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 3, int((/1,2,3/),kind=int1)) + Istatus_Sum
       Istatus_Sum = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 30, "no_snow_or_ice "// &
                                " sea_ice "// &
                                " snow ") + Istatus_Sum
     end if

end subroutine WRITE_FLAG_ATTRIBUTES

!------------------------------------------------------------------------------------------------
!  write the sds attributes that communicate how data is scaled
!------------------------------------------------------------------------------------------------
subroutine WRITE_SCALING_ATTRIBUTES(Sds,Istatus_Sum)
    type(Sds_Struct), intent(in):: Sds
    integer, intent(inout):: Istatus_Sum
    integer:: sfsnatt


     !--- determined if a scaled Sds, if so write needed attributes for scaling
     if (Sds%Scaling_Type > 0) then

!      Istatus_Sum = sfsnatt(Sds%Id_Output, "actual_missing", DFNT_FLOAT32, 1, Sds%Unscaled_Missing) + Istatus_Sum
      Istatus_Sum = sfsnatt(Sds%Id_Output, "actual_range", DFNT_FLOAT32, 2, Sds%Actual_Range) + Istatus_Sum
      if (Sds%Data_Type == DFNT_INT8) then
         Istatus_Sum = sfsnatt(Sds%Id_Output, "valid_range", Sds%Data_Type, 2, int(Sds%valid_range,kind=int1)) + Istatus_Sum
      elseif (Sds%Data_Type == DFNT_INT16) then
         Istatus_Sum = sfsnatt(Sds%Id_Output, "valid_range", Sds%Data_Type, 2, int(Sds%valid_range,kind=int2)) + Istatus_Sum
      elseif (Sds%Data_Type == DFNT_FLOAT32) then
         Istatus_Sum = sfsnatt(Sds%Id_Output, "valid_range", Sds%Data_Type, 2, Sds%valid_range) + Istatus_Sum
     endif



      !--- write remaining attributes
      if (Sds%Scaling_Type == 1) then
       Istatus_Sum = sfsnatt(Sds%Id_Output, "scale_factor", DFNT_FLOAT32, 1, Sds%Scale_Factor) + Istatus_Sum
       Istatus_Sum = sfsnatt(Sds%Id_Output, "add_offset", DFNT_FLOAT32, 1, Sds%Add_Offset) + Istatus_Sum
!     else  !deprecated
!      Istatus_Sum = sfsnatt(Sds%Id_Output, "RANGE_MIN", DFNT_FLOAT32, 1, Sds%Unscaled_Min) + Istatus_Sum
!      Istatus_Sum = sfsnatt(Sds%Id_Output, "RANGE_MAX", DFNT_FLOAT32, 1, Sds%Unscaled_Max) + Istatus_Sum
!      Istatus_Sum = sfsnatt(Sds%Id_Output, "SCALED_MIN", DFNT_INT32, 1, Sds%Scaled_Min) + Istatus_Sum
!      Istatus_Sum = sfsnatt(Sds%Id_Output, "SCALED_MAX", DFNT_INT32, 1, Sds%Scaled_Max) + Istatus_Sum
!      Istatus_Sum = sfsnatt(Sds%Id_Output, "SCALED_MISSING", DFNT_INT32, 1, Sds%Scaled_Missing) + Istatus_Sum
      endif

     endif

     !--- write out fill value for all but packed variables
     if (Sds%Unscaled_Missing /= No_Attribute_Missing_Value) then
         if (Sds%Data_Type == DFNT_INT8) then
           Istatus_Sum = sfsnatt(Sds%Id_Output, "_FillValue", Sds%Data_Type, 1, int(Sds%Fill_Value,kind=int1)) + Istatus_Sum
         elseif (Sds%Data_Type == DFNT_INT16) then
           Istatus_Sum = sfsnatt(Sds%Id_Output, "_FillValue", Sds%Data_Type, 1, int(Sds%Fill_Value,kind=int2)) + Istatus_Sum
         elseif (Sds%Data_Type == DFNT_FLOAT32) then
           Istatus_Sum = sfsnatt(Sds%Id_Output, "_FillValue", Sds%Data_Type, 1, Sds%Fill_Value) + Istatus_Sum
         endif
     endif

end subroutine WRITE_SCALING_ATTRIBUTES
   !------------------------------------------------------------------------------------------------
   !  write the sds attributes that communicate how data is scaled
   !------------------------------------------------------------------------------------------------
   subroutine READ_SCALING_ATTRIBUTES(Sds,Istatus)
      type(Sds_Struct), intent(inout):: Sds
      integer, intent(inout):: Istatus

      integer(kind=int1):: fill_value_i1
      integer(kind=int2):: fill_value_i2
      integer(kind=int4):: fill_value_i4
      real(kind=real4):: fill_value_r4
      integer(kind=int1), dimension(2):: valid_range_i1
      integer(kind=int2), dimension(2):: valid_range_i2
      integer(kind=int4), dimension(2):: valid_range_i4
      real(kind=real4), dimension(2):: valid_range_r4

      integer:: sffattr, sfrnatt, sfrcatt
      
      
      !--- read Sds attributes
      Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"SCALED"), Sds%Scaling_Type)
      
      Istatus = sfrcatt(Sds%Id_Input, sffattr(Sds%Id_Input,"units"), Sds%Units)
      if (Istatus /= 0) Sds%Units = "none"
      Istatus = sfrcatt(Sds%Id_Input, sffattr(Sds%Id_Input,"standard_name"), Sds%Standard_Name)
      if (Istatus /= 0) Sds%Standard_Name = "none"
      Istatus = sfrcatt(Sds%Id_Input, sffattr(Sds%Id_Input,"long_name"), Sds%Long_Name)
      if (Istatus /= 0) Sds%Long_Name = "none"
       
      if (Sds%Scaling_Type /= sym%NO_SCALING) then
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"actual_missing"), Sds%Unscaled_Missing)
         print*,'a22',istatus
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"actual_range"), Sds%Actual_Range)
      end if
       print*,'a2',istatus,Sds%Actual_Range
      if (Sds%Scaling_Type /= sym%NO_SCALING) then
         if (Sds%Data_Type == DFNT_INT8) then
            Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"valid_range"), valid_range_i1)
            Sds%Valid_Range = int(valid_range_i1,kind=int4)
         elseif (Sds%Data_Type == DFNT_INT16) then
            Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"valid_range"), valid_range_i2)
            Sds%Valid_Range = int(valid_range_i2,kind=int4)
         elseif (Sds%Data_Type == DFNT_INT32) then
            Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"valid_range"), valid_range_i4)
            Sds%Valid_Range = int(valid_range_i4,kind=int4)
         elseif (Sds%Data_Type == DFNT_FLOAT32) then
            Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"valid_range"), valid_range_r4)
            Sds%Valid_Range = int(valid_range_r4,kind=real4)
         endif
      endif
 print*,'a2',istatus
      !--- fill value (if absent, assume this is packed variable and set
      !--- Unscaled_Missing to -888
      if (Sds%Data_Type == DFNT_INT8) then
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"_FillValue"), fill_value_i1)
         if (Istatus == 0) then
            Sds%Fill_Value = int(fill_value_i1,kind=int4)
         else
            Sds%Fill_Value = Missing_Value_int4
            Sds%Unscaled_Missing = No_Attribute_Missing_Value
         endif
      elseif (Sds%Data_Type == DFNT_INT16) then
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"_FillValue"), fill_value_i2)
         if (Istatus == 0) then
            Sds%Fill_Value = int(fill_value_i2,kind=int4)
         else
           Sds%Fill_Value = Missing_Value_int4
           Sds%Unscaled_Missing = No_Attribute_Missing_Value
         endif
      elseif (Sds%Data_Type == DFNT_INT32) then
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"_FillValue"), fill_value_i4)
         if (Istatus == 0) then
            Sds%Fill_Value = int(fill_value_i4,kind=int4)
         else
           Sds%Fill_Value = Missing_Value_int4
           Sds%Unscaled_Missing = No_Attribute_Missing_Value
         endif
      elseif (Sds%Data_Type == DFNT_FLOAT32) then
         Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"_FillValue"), fill_value_r4)
         if (Istatus == 0) then
            Sds%Fill_Value = int(fill_value_r4,kind=int4)
         else
           Sds%Fill_Value = Missing_Value_int4
           Sds%Unscaled_Missing = No_Attribute_Missing_Value
         endif
      endif
 print*,'a4',istatus
      !-- if scaled, read attributes that allow unscaling
      if (Sds%Scaling_Type > 0) then
         if (Sds%Scaling_Type == sym%LINEAR_SCALING) then
            Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"scale_factor"), Sds%Scale_Factor) + Istatus
            Istatus = sfrnatt(Sds%Id_Input, sffattr(Sds%Id_Input,"add_offset"), Sds%Add_Offset) + Istatus
         endif
      endif
print*,'a5',istatus
   end subroutine READ_SCALING_ATTRIBUTES


end module LEVEL2B_ROUTINES
