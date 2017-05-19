!$Id$
!==============================================================================
! Module for a simple cloud optical depth for use in masking and typing
!==============================================================================
module SIMPLE_COD

   use CONSTANTS
   use PIXEL_COMMON
   use NUMERICAL_ROUTINES
   use FILE_UTILITY
   use SURFACE_PROPERTIES

   implicit none

   public:: COMPUTE_SIMPLE_SOLAR_COD
   public:: COMPUTE_SIMPLE_LUNAR_COD

   private:: READ_LUT
 
   integer, private, save:: Table_Read_Flag = 0 
   real, private, parameter:: SOLZEN_LIMIT = 80.0
   integer, private, save:: Number_Solzen
   integer, private, save:: Number_Senzen
   integer, private, save:: Number_Relaz
   integer, private, save:: Number_Opd

   real, private, save, dimension(:),allocatable:: Opd_Lut
   real, private, save, dimension(:),allocatable:: Solzen_Lut
   real, private, save, dimension(:),allocatable:: Senzen_Lut
   real, private, save, dimension(:),allocatable:: Relaz_Lut
   real, private, save, dimension(:,:,:,:),allocatable:: Ref_Lut
   real, private, save, dimension(:),allocatable:: Spherical_Albedo_Lut
   real, private, save, dimension(:,:),allocatable:: Albedo_Lut
   real, private, save, dimension(:,:),allocatable:: Transmission_Lut
   real, private, save, dimension(:),allocatable:: Ref_Vector
   real, private, save, dimension(:),allocatable:: Ref_Toa_Vector
   real, private, save, dimension(:),allocatable:: Temp_Vector
   real, private, save:: Solzen_Delta
   real, private, save:: Senzen_Delta
   real, private, save:: Relaz_Delta
   real, private, save:: Opd_Delta
   character(len=100):: filename = 'simple_063um_ref_lut_wat_cld.bin'

   contains

   !===========================================================================
   ! simplistic optical depth assuming a liquid water cloud with Reff = 10 um
   ! This uses the solar source
   !===========================================================================
   subroutine COMPUTE_SIMPLE_SOLAR_COD(Number_Elements, Number_Lines)

      integer, intent(in):: Number_Lines
      integer, intent(in):: Number_Elements
      integer:: Solzen_Idx, Senzen_Idx, Relaz_Idx, Opd_Idx
      integer:: Read_Table_Error
      real:: Ref_Toa, Opd, dRef, dOpd_dRef, Alb_Sfc, Ref_Toa_Temp
      integer:: Elem_Idx, Line_Idx
      logical:: Negative_Opd

      if (Sensor%Chan_On_Flag_Default(1) == sym%NO)  return

      if (Table_Read_Flag == 0) then
         call READ_LUT(Read_Table_Error)
         if (Read_Table_Error /= 0) return
         Table_Read_Flag = 1
      endif

      Element_Loop: do Elem_Idx = 1, Number_Elements
         Line_Loop: do Line_Idx = 1, Number_Lines

            if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == 1) cycle
            if (Geo%Solzen(Elem_Idx,Line_Idx) > SOLZEN_LIMIT) cycle

            Ref_Toa = ch(1)%Ref_Toa(Elem_Idx,Line_Idx) / 100.0

            Solzen_Idx = int((Geo%Solzen(Elem_Idx,Line_Idx) - Solzen_Lut(1))/Solzen_Delta) + 1
            Solzen_Idx = min(max(1,Solzen_Idx),Number_Solzen)

            Senzen_Idx = int((Geo%Satzen(Elem_Idx,Line_Idx) - Senzen_Lut(1))/Senzen_Delta) + 1
            Senzen_Idx = min(max(1,Senzen_Idx),Number_Senzen)

            Relaz_Idx = int((Geo%Relaz(Elem_Idx,Line_Idx) - Relaz_Lut(1))/Relaz_Delta) + 1
            Relaz_Idx = min(max(1,Relaz_Idx),Number_Relaz)

            Ref_Vector = Ref_Lut(:,Solzen_Idx,Senzen_Idx,Relaz_Idx)

            !--- Account for surface reflection and atmospheric transmission
            Alb_Sfc = ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx)
            if (Alb_Sfc < 0.0) then
               Alb_Sfc = Ch1_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx))
            endif
            if (Sfc%Snow(Elem_Idx,Line_Idx) /=sym%NO_SNOW) then
               Alb_Sfc = Ch1_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx))
            endif
            Alb_Sfc = Alb_Sfc / 100.0
            Temp_Vector = Alb_Sfc / (1.0 - Alb_Sfc * Spherical_Albedo_Lut)
            Ref_Toa_Vector = (Ref_Vector + Temp_Vector * Transmission_Lut(:,Solzen_Idx)*Transmission_Lut(:,Senzen_Idx))

            !--------------------------------------------------------
            !handle negative optical depths
            !--------------------------------------------------------
            Ref_Toa_Temp = Ref_Toa
            Negative_Opd = .false.
            if (Ref_Toa < Ref_Toa_Vector(1)) then
               Ref_Toa_Temp = Ref_Toa_Vector(1) + (Ref_Toa_Vector(1)-Ref_Toa) 
               Negative_Opd = .true.
            endif 

            call LOCATE(Ref_Toa_Vector, Number_Opd, Ref_Toa_Temp, Opd_Idx)
            Opd_Idx = min(Number_Opd-1,max(1,Opd_Idx))

            dRef = Ref_Toa_Vector(Opd_Idx+1)-Ref_Toa_Vector(Opd_Idx)
            dOpd_dRef = 0.0
            if (dRef > 0) then
               dOpd_dRef = (Opd_Lut(Opd_Idx+1) - Opd_Lut(Opd_Idx))/dRef
            endif
            Opd = Opd_Lut(Opd_Idx) + dOpd_dRef * (Ref_Toa_Temp - Ref_Toa_Vector(Opd_Idx))
            Opd = min(max(Opd_Lut(1),Opd),Opd_Lut(Number_Opd))
            ch(1)%Opd(Elem_Idx,Line_Idx) = 10.0**Opd

            if (Negative_Opd) then
               ch(1)%Opd(Elem_Idx,Line_Idx) = Opd_Lut(1) - ch(1)%Opd(Elem_Idx,Line_Idx) 
           endif

         enddo Line_Loop
      enddo Element_Loop

   end subroutine COMPUTE_SIMPLE_SOLAR_COD
   !===========================================================================
   ! simplistic optical depth assuming a liquid water cloud with Reff = 10 um
   ! This uses the lunar source
   ! note, we should combine this with solar version 
   !===========================================================================
   subroutine COMPUTE_SIMPLE_LUNAR_COD(Number_Elements, Number_Lines)

      integer, intent(in):: Number_Lines
      integer, intent(in):: Number_Elements
      integer:: Solzen_Idx, Senzen_Idx, Relaz_Idx, Opd_Idx
      integer:: Read_Table_Error
      real:: Ref_Toa, Opd, dRef, dOpd_dRef, Alb_Sfc, Ref_Toa_Temp
      integer:: Elem_Idx, Line_Idx
      logical:: Negative_Opd

      if (Sensor%Chan_On_Flag_Default(44) == sym%NO)  return
      if (Sensor%Chan_On_Flag_Default(1) == sym%NO)  return     !issue - we need ch1 and ch2 white sky
      if (Sensor%Chan_On_Flag_Default(2) == sym%NO) return 

      if (Table_Read_Flag == 0) then
         call READ_LUT(Read_Table_Error)
         if (Read_Table_Error /= 0) return
         Table_Read_Flag = 1
      endif

      Element_Loop: do Elem_Idx = 1, Number_Elements
         Line_Loop: do Line_Idx = 1, Number_Lines

            if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == 1) cycle
            if (Geo%Lunzen(Elem_Idx,Line_Idx) > SOLZEN_LIMIT) cycle

            Ref_Toa = ch(44)%Ref_Lunar_Toa(Elem_Idx,Line_Idx) / 100.0

            Solzen_Idx = int((Geo%Lunzen(Elem_Idx,Line_Idx) - Solzen_Lut(1))/Solzen_Delta) + 1
            Solzen_Idx = min(max(1,Solzen_Idx),Number_Solzen)

            Senzen_Idx = int((Geo%Satzen(Elem_Idx,Line_Idx) - Senzen_Lut(1))/Senzen_Delta) + 1
            Senzen_Idx = min(max(1,Senzen_Idx),Number_Senzen)

            Relaz_Idx = int((Geo%LunRelaz(Elem_Idx,Line_Idx) - Relaz_Lut(1))/Relaz_Delta) + 1
            Relaz_Idx = min(max(1,Relaz_Idx),Number_Relaz)

            Ref_Vector = Ref_Lut(:,Solzen_Idx,Senzen_Idx,Relaz_Idx)

            !--- Account for surface reflection and atmospheric transmission
            if (ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4 .and.  &
                ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Alb_Sfc = 0.5*(ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx)+ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx))
            else
              Alb_Sfc = 0.5*(Ch1_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) + Ch2_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)))
            endif

            if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then
              Alb_Sfc = 0.5*(Ch1_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) &
                           + Ch2_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)))
            endif


            Alb_Sfc = Alb_Sfc / 100.0
            Temp_Vector = Alb_Sfc / (1.0 - Alb_Sfc * Spherical_Albedo_Lut)
            Ref_Toa_Vector = (Ref_Vector + Temp_Vector * Transmission_Lut(:,Solzen_Idx)*Transmission_Lut(:,Senzen_Idx))

            !--------------------------------------------------------
            !handle negative optical depths
            !--------------------------------------------------------
            Ref_Toa_Temp = Ref_Toa
            Negative_Opd = .false.
            if (Ref_Toa < Ref_Toa_Vector(1)) then
               Ref_Toa_Temp = Ref_Toa_Vector(1) + (Ref_Toa_Vector(1)-Ref_Toa) 
               Negative_Opd = .true.
            endif 

            call LOCATE(Ref_Toa_Vector, Number_Opd, Ref_Toa_Temp, Opd_Idx)
            Opd_Idx = min(Number_Opd-1,max(1,Opd_Idx))

            dRef = Ref_Toa_Vector(Opd_Idx+1)-Ref_Toa_Vector(Opd_Idx)
            dOpd_dRef = 0.0
            if (dRef > 0) then
               dOpd_dRef = (Opd_Lut(Opd_Idx+1) - Opd_Lut(Opd_Idx))/dRef
            endif
            Opd = Opd_Lut(Opd_Idx) + dOpd_dRef * (Ref_Toa_Temp - Ref_Toa_Vector(Opd_Idx))
            Opd = min(max(Opd_Lut(1),Opd),Opd_Lut(Number_Opd))
            ch(44)%Opd(Elem_Idx,Line_Idx) = 10.0**Opd

            if (Negative_Opd) then
               ch(44)%Opd(Elem_Idx,Line_Idx) = Opd_Lut(1) - ch(44)%Opd(Elem_Idx,Line_Idx) 
            endif

         enddo Line_Loop
      enddo Element_Loop

   end subroutine COMPUTE_SIMPLE_LUNAR_COD

   !===========================================================================
   ! Read Routine for the LUT
   !===========================================================================
   subroutine READ_LUT(ierror)

       integer, intent(out):: ierror
       integer:: lun

       lun = get_lun()

       open(unit=lun,file=trim(ancil_data_dir)//"static/luts/simple_cod_table/"//trim(filename), &
            action="read",form="unformatted",access="stream",iostat=ierror)

       if (ierror /= 0) return

       read(unit=lun) Number_Opd, Number_Solzen, Number_Senzen, Number_Relaz

       allocate(Opd_Lut(Number_Opd))
       allocate(Solzen_Lut(Number_Solzen))
       allocate(Senzen_Lut(Number_Senzen))
       allocate(Relaz_Lut(Number_Relaz))
       allocate(Ref_Lut(Number_Opd,Number_Solzen,Number_Senzen,Number_Relaz))
       allocate(Spherical_Albedo_Lut(Number_Opd))
       allocate(Albedo_Lut(Number_Opd,Number_Solzen))
       allocate(Transmission_Lut(Number_Opd,Number_Solzen))
       allocate(Ref_Vector(Number_Opd))
       allocate(Ref_Toa_Vector(Number_Opd))
       allocate(Temp_Vector(Number_Opd))

       read(unit=lun) Opd_Lut, Solzen_Lut, Senzen_Lut, Relaz_Lut,  &
                      Spherical_Albedo_Lut, Albedo_Lut, Transmission_Lut, Ref_Lut

       close(unit=lun)

       Solzen_Delta = Solzen_Lut(2) - Solzen_Lut(1)
       Senzen_Delta = Senzen_Lut(2) - Senzen_Lut(1)
       Relaz_Delta = Relaz_Lut(2) - Relaz_Lut(1)
       Opd_Delta = Opd_Lut(2) - Opd_Lut(1)

!      print *, "Table Read ", Number_Opd, Number_Solzen, Number_Senzen, Number_Relaz
!      print *, "Opd_Lut = ", Opd_Lut
!      print *, "Solzen_Lut = ", Solzen_Lut
!      print *, "Relaz_Lut = ", Relaz_Lut
   end subroutine READ_LUT
  
end module SIMPLE_COD
