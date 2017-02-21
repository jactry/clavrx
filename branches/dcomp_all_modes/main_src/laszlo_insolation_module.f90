!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: laszlo_insolation_module.f90 (src)
!       laszlo_insolation (program)
!
! PURPOSE: Insolation Module
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
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
!--------------------------------------------------------------------------------------
module LASZLO_INSOLATION
      use CX_CONSTANTS_MOD
      use PIXEL_COMMON
      use NWP_COMMON
      use FILE_TOOLS, only: &
         get_lun

   implicit none
   
   private
 
   public:: INSOLATION

   integer, parameter:: Maxfpr = 5
   integer, parameter:: Maxint = 3
   integer, parameter:: Nfpar = 5
   integer, parameter:: Nintv = 3
   real, parameter:: Solzen_Limit_Sasrab = 78.4

contains

   subroutine INSOLATION(Line_Idx_Min_Segment,Number_Of_Lines_Read_This_Segment)
      integer, intent(in):: Line_Idx_Min_Segment
      integer, intent(in):: Number_Of_Lines_Read_This_Segment
      
      character(len=1020):: Auxpath
      real(kind=real4):: Glat
      real(kind=real4):: Glon
      integer(kind=int4):: Year
      integer(kind=int4):: Month_Local
      integer(kind=int4):: Jday
      real(kind=real4):: Gmtime
      real(kind=real4):: Snowfr
      real(kind=real4):: Cdfrac
      real(kind=real4):: Cldref
      real(kind=real4):: Clrref
      real(kind=real4):: Cmpref
      real(kind=real4):: Ozone
      real(kind=real4):: Pwater
      real(kind=real4):: Solmu
      real(kind=real4):: Satmu
      real(kind=real4):: Relaz_Local
      real(kind=real4):: Glint
      integer(kind=int4):: Isatid
      real(kind=real4):: Misval
      real(kind=real4):: Cldtau
      real(kind=real4):: Aertau
      real (kind=real4), dimension(Maxfpr,Maxint) :: flxall
      real (kind=real4), dimension(Maxfpr,Maxint) :: flxclr

      integer:: Elem_Idx
      integer:: Line_Idx
      integer:: First_Element
      integer:: Number_of_Elements
      integer:: Number_of_Lines

      First_Element = Line_Idx_Min_Segment
      Number_of_Elements = Image%Number_Of_Elements
      Number_of_Lines = Number_Of_Lines_Read_This_Segment

      !----------------------------------------------------------
      ! loop over each pixel in segment
      !----------------------------------------------------------
      Insolation_All_Sky = Missing_Value_Real4
      Insolation_All_Sky_Diffuse = Missing_Value_Real4
      Insolation_Clear_Sky = Missing_Value_Real4
      Insolation_Cld_Opd = Missing_Value_Real4
      Insolation_Aer_Opd = Missing_Value_Real4

      !--- check for required channels
      if ( .not. Sensor%Chan_On_Flag_Default(1) ) return
      
      Year = Image%Start_Year
      Month_Local = month
      Jday = Image%Start_Doy
      

      element_loop: do Elem_Idx = First_Element, First_Element + Number_of_Elements - 1
         line_loop: do Line_Idx = 1, Number_of_Lines

            !--- check for a bad pixel
            if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) ) then
               cycle
            endif
            
            if (Space_Mask(Elem_Idx,Line_Idx) == sym%YES) then
               cycle
            end if
            
            if (Geo%Solzen(Elem_Idx,Line_Idx) > Solzen_Limit_SASRAB) then
               Insolation_All_Sky(Elem_Idx,Line_Idx) = 0.0
               Insolation_All_Sky_Diffuse(Elem_Idx,Line_Idx) = 0.0
               Insolation_Clear_Sky(Elem_Idx,Line_Idx) = 0.0
               Insolation_Cld_Opd(Elem_Idx,Line_Idx) = 0.0
               Insolation_Aer_Opd(Elem_Idx,Line_Idx) = 0.0
               cycle
            end if

            Glat = Nav%Lat(Elem_Idx,Line_Idx)
            Glon = Nav%Lon(Elem_Idx,Line_Idx)
            
            Gmtime = Utc_Scan_Time_Hours(Line_Idx)
            if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then
               Snowfr = 1.0
            else
               Snowfr = 0.0
            end if
            
     
            if ((CLDMASK%Cld_Mask(Elem_Idx,Line_Idx) == sym%CLOUDY) .or. &
                (CLDMASK%Cld_Mask(Elem_Idx,Line_Idx) == sym%PROB_CLOUDY)) then
               Cdfrac = 1.0
               Cldref = ch(1)%Ref_Toa(Elem_Idx,Line_Idx)/100.0
               Clrref = Missing_Value_Real4
            else
               Cdfrac = 0.0
               Clrref = ch(1)%Ref_Toa(Elem_Idx,Line_Idx)/100.0
               Cldref = Missing_Value_Real4
            end if

            if ((trim(Dark_Composite_Name) /= "no_file") .and. &
                  (Ref_Ch1_Dark_Composite(Elem_Idx,Line_Idx) /= Missing_Value_Real4)) then
               Cmpref = Ref_Ch1_Dark_Composite(Elem_Idx,Line_Idx)/100.0
            else
               Cmpref = ch(1)%Ref_Toa_Clear(Elem_Idx,Line_Idx)/100.0
            end if

            Ozone = Ozone_Nwp_Pix(Elem_Idx,Line_Idx)  !Dobson
            !--- convert Ozone into atm-cm
            if (Ozone /= Missing_Value_Real4) Ozone = Ozone / 1000.0

            Pwater = Tpw_Nwp_Pix(Elem_Idx,Line_Idx)   !g/cm^2 or cm
            Solmu = Geo%Cossolzen(Elem_Idx,Line_Idx)
            Satmu = Geo%Coszen(Elem_Idx,Line_Idx)
            Relaz_local = Geo%Relaz(Elem_Idx,Line_Idx) !convention?
            IF ( Relaz_Local /= Missing_Value_Real4 ) Relaz_Local = 180.0 - Relaz_Local
            Glint = Geo%Glintzen(Elem_Idx,Line_Idx)

            !---------------------------------------------------
            ! Isatid : satellite id.
            !                            1 - noaa-7/ch1
            !                            2 - goes-6/ch1
            !                            3 - goes-5/ch1
            !                            4 - meteosat-2/ch1
            !                            5 - gms-2/ch1
            !                            6 - goes-8/ch1
            !---------------------------------------------------
            Isatid = 6                 !Default to GOES

            select case (Sensor%WMO_Id)

            case(3:5,200:209,223,706:708)       !AVHRR
               Isatid = 1

            case(252:259)       !GOES
               Isatid = 6

            case(55:57)         !MSG
               Isatid = 4

            case(171:172)       !GMS (MTSAT?)
               Isatid = 5

            case(224)           !VIIRS
               Isatid = 1

            case(783:784)       !MODIS
               Isatid = 1

            case(810)           !COMS
               Isatid = 6

            end select
 
            Misval = Missing_Value_Real4

            !------------------------------------------------------
            !  Qc
            !------------------------------------------------------
            if  (Cdfrac == 0.0 .and. Clrref <= 0.0) cycle 
            if  (Cdfrac == 1.0 .and. Cldref <= 0.0) cycle 
            if  (Cmpref <= 0.0) cycle 
     

            !------------------------------------------------------
            !  call the main insolation algorithm routine
            !------------------------------------------------------
            Auxpath = trim(Ancil_Data_Dir)//"static/sasrab/"


            call Sasrab(Auxpath, &
                 Glat, Glon, Year, Month_Local, Jday, Gmtime,  &
                 Snowfr, Cdfrac, Cldref, Clrref, Cmpref, Ozone,        &
                 Pwater, Solmu, Satmu, Relaz_Local, Glint, Isatid,     &
                 Misval, Maxfpr, Maxint, Flxall, Flxclr, Aertau, &
                 Cldtau, Nfpar, Nintv )

            !---- store sasrab output into correct output variables
            Insolation_All_Sky(Elem_Idx,Line_Idx) = Flxall(3,1)
            Insolation_All_Sky_Diffuse(Elem_Idx,Line_Idx) = Flxall(5,1)
            Insolation_Clear_Sky(Elem_Idx,Line_Idx) = Flxclr(3,1)
            Insolation_Cld_Opd(Elem_Idx,Line_Idx) = Cldtau
            Insolation_Aer_Opd(Elem_Idx,Line_Idx) = Aertau

         end do line_loop
      end do element_loop

    
   end subroutine INSOLATION

end module LASZLO_INSOLATION
