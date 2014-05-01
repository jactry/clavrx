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
     use CONSTANTS
     use PIXEL_COMMON
     use NWP_COMMON
     use FILE_UTILITY

 implicit none
 
 public:: INSOLATION

 integer, private, parameter:: Maxfpr = 5
 integer, private, parameter:: Maxint = 3
 integer, private, parameter:: Nfpar = 5
 integer, private, parameter:: Nintv = 3
 real, private, parameter:: Solzen_Limit_Sasrab = 78.4

 contains

  subroutine INSOLATION(Line_Idx_Min_Segment,Num_Scans_Read)
    integer, intent(in):: Line_Idx_Min_Segment
    integer, intent(in):: Num_Scans_Read
    character(len=120):: Auxpath
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
    Number_of_Elements = Num_Pix
    Number_of_Lines = Num_Scans_Read

!----------------------------------------------------------
! loop over each pixel in segment
!----------------------------------------------------------
Insolation_All_Sky = Missing_Value_Real4
Insolation_All_Sky_Diffuse = Missing_Value_Real4
Insolation_Clear_Sky = Missing_Value_Real4
Insolation_Cld_Opd = Missing_Value_Real4
Insolation_Aer_Opd = Missing_Value_Real4

!--- check for required channels
if (Chan_On_Flag_Default(1) == sym%NO) return

element_loop: do Elem_Idx = First_Element, First_Element + Number_of_Elements - 1
line_loop: do Line_Idx = 1, Number_of_Lines

     !--- check for a bad pixel
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
             cycle
     endif
     if (Space_Mask(Elem_Idx,Line_Idx) == sym%YES) then
             cycle
     endif
     if (Solzen(Elem_Idx,Line_Idx) > Solzen_Limit_SASRAB) then
             Insolation_All_Sky(Elem_Idx,Line_Idx) = 0.0
             Insolation_All_Sky_Diffuse(Elem_Idx,Line_Idx) = 0.0
             Insolation_Clear_Sky(Elem_Idx,Line_Idx) = 0.0
             Insolation_Cld_Opd(Elem_Idx,Line_Idx) = 0.0
             Insolation_Aer_Opd(Elem_Idx,Line_Idx) = 0.0
             cycle
     endif

     Glat = Lat(Elem_Idx,Line_Idx)
     Glon = Lon(Elem_Idx,Line_Idx)
     Year = Start_Year
     Month_Local = month
     Jday = start_day
     Gmtime = Utc_Scan_Time_Hours(Line_Idx)
     if (Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then
             Snowfr = 1.0
     else
             Snowfr = 0.0
     endif
     !if (Posterior_Cld_Probability(Elem_Idx,Line_Idx) >= 0.5) then
     if ((Cld_Mask(Elem_Idx,Line_Idx) == sym%CLOUDY) .or. (Cld_Mask(Elem_Idx,Line_Idx) == sym%PROB_CLOUDY)) then
             Cdfrac = 1.0
             Cldref = ch(1)%Ref_Toa(Elem_Idx,Line_Idx)/100.0
             Clrref = Missing_Value_Real4
     else
             Cdfrac = 0.0
             Clrref = ch(1)%Ref_Toa(Elem_Idx,Line_Idx)/100.0
             Cldref = Missing_Value_Real4
     endif

     if ((trim(Dark_Composite_Name) /= "no_file") .and. &
         (Ref_Ch1_Dark_Composite(Elem_Idx,Line_Idx) /= Missing_Value_Real4)) then
       Cmpref = Ref_Ch1_Dark_Composite(Elem_Idx,Line_Idx)/100.0
     else
       Cmpref = ch(1)%Ref_Toa_Clear(Elem_Idx,Line_Idx)/100.0
     endif

     Ozone = Ozone_Nwp_Pix(Elem_Idx,Line_Idx)  !Dobson
     !--- convert Ozone into atm-cm
     if (Ozone /= Missing_Value_Real4) Ozone = Ozone / 1000.0

     Pwater = Tpw_Nwp_Pix(Elem_Idx,Line_Idx)   !g/cm^2 or cm
     Solmu = Cossolzen(Elem_Idx,Line_Idx)
     Satmu = Coszen(Elem_Idx,Line_Idx)
     Relaz_local = Relaz(Elem_Idx,Line_Idx) !convention?
     IF ( Relaz_Local /= Missing_Value_Real4 ) Relaz_Local = 180.0 - Relaz_Local
     Glint = Glintzen(Elem_Idx,Line_Idx)

     Isatid = 6
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

     select case (Sc_Id_WMO)

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
     Auxpath = trim(Ancil_Data_Dir)//"insolation/"

!    Glat = 60.0
!    Glon = -163.9
!    Year = 2009
!    Month_Local = 8
!    Jday = 213
!    Gmtime = 22.0
!    Snowfr = 0.0
!    Cdfrac = 1.0
!    Cldref = 0.56
!    Clrref = missing_value_real4  !0.05
!    Cmpref = 0.10
!    Ozone = 0.3
!    Pwater = 0.2
!    Solmu = 0.8
!    Satmu = 0.3
!    Relaz_Local = 11.0
!    Glint = 115.0
!    Isatid = 6

!    print *, "calling sasrab with ", &
!                Glat, Glon, Year, Month_Local, Jday, Gmtime,  &
!                Snowfr, Cdfrac, Cldref, Clrref, Cmpref, Ozone,        &
!                Pwater, Solmu, Satmu, Relaz_Local, Glint, Isatid,     &
!                Misval, Maxfpr, Maxint

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
