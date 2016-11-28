!$Id: goes_area_io_mod.f90,v 1.1.1.1 2012/05/18 20:18:40 heidinger Exp $
module GOES_AREA_IO_MODULE
use CONSTANTS
implicit none

public:: GET_HEADERS,GET_IMAGE
private:: READ_IMGR_RP,READ_IMGR_SIN,READ_IMGR_MON,FIX_VAR,FIX_VAR2
private

!-----------------------------------------------------------------------------
! define derived data types used for holding GVAR parameters
!-----------------------------------------------------------------------------
!
!C	Imager ONA repeat sinusoid T.
!C

        type, public:: IMGR_SIN 
          integer :: mag_sinu
          integer :: phase_ang_sinu
        end type IMGR_SIN
!C
!C	Imager repeat monomial T.
!C
        type, public:: IMGR_MON
          integer :: order_appl_sinu
          integer :: order_mono_sinu
          integer :: mag_mono_sinu
          integer :: phase_ang_sinu
          integer :: ang_from_epoch
        end type IMGR_MON
!C
!C	Imager repeat T.
!C
        type, public:: IMGR_RP 
           integer :: exp_mag
           integer :: exp_time_const
           integer :: mean_att_ang_const
           integer :: num_sinu_per_angle
           type(IMGR_SIN), dimension(15):: sinusoid
           integer :: num_mono_sinu
           type(IMGR_MON), dimension(4) :: monomial
        end type IMGR_RP
!
!     Define McIDAS area struct%re
!
      type, public:: AREA_STRUCT
	integer ::area_status           ! Area status
	integer ::version_num       ! Area version number
	integer ::sat_id_num        ! Satellite ID (SSS)
	integer ::img_date          ! Year and Julian day (YYDDD)
	integer ::img_time          ! Time of image (HHMMSS)
	integer ::north_bound       ! Upper left line in Sat coords (Y-Coord)
	integer ::west_vis_pixel    ! Upper left element in Sat coords (X-Coord)
	integer ::z_coor            ! Upper z-coord (Z-Coord)
	integer ::num_line          ! Number of lines in image (Y-SIZE)
	integer ::num_elem          ! Number of elememts in image (X-SIZE)
	integer ::bytes_per_pixel   ! Number of bytes per data element
	integer ::line_res          ! Line resolution (Y-RES) 
	integer ::elem_res          ! Element resolution (X-RES)
	integer ::num_chan          ! Number of bands (Z-RES)
	integer ::num_byte_ln_prefix ! Number of bytes in line prefix (multiple of 4)
	integer ::proj_num	        ! Project number
	integer ::creation_date	        ! Creation date (YYDDD)
	integer ::creation_time         ! Creation time (HHMMSS)
	integer ::sndr_filter_map       ! Filter map for soundings
	integer ::img_id_num            ! Image ID number
	integer, dimension(4) :: id     ! Reserved for radar appications
	character(len=32):: comment	! 32 char comments
	integer ::pri_key_calib	        ! Calibration colicil (area number)
	integer ::pri_key_nav	        ! Primary navigation codicil  (data)
	integer ::sec_key_nav	        ! Secondary navigation codicil (nav)
	integer ::val_code              ! Validity code
	integer, dimension(8) :: pdl    ! PDL in packed-byte format
	integer ::band8	                ! Where band-8 came from
	integer ::act_img_date	        ! Actual image start day (YYDDD)
	integer ::act_img_time	        ! Actual image start time (HHMMSS)
	integer ::act_start_scan        ! Actual start scan
	integer ::len_prefix_doc        ! Length of prefix documentation (DOC)
	integer ::len_prefix_calib      ! Length of prefix calibration (CAL)
	integer ::len_prefix_lev        ! Length of prefix level (LEV)
	character(len=4)::src_type      ! Source type
	character(len=4)::calib_type    ! Calibration type
	integer ::avg_or_sample	        ! Data was averaged (0) or sampled (1)
	integer ::poes_signal	        ! LAC, GAC, HRPT
	integer ::poes_up_down	        ! POES ascending/descending
	character(len=4) :: orig_src_type   ! Original source type of data
	integer, dimension(7) ::reserved    ! Reserved (6 calib pointer)
      end type AREA_STRUCT
!
!	Define GVAR Navigation block
!
        type, public:: GVAR_NAV
          character(len=4):: nav_type
          character(len=4):: IMC_status
          integer, dimension(12) ::  spare1
          integer ::    stat
          integer ::    ref_long
          integer ::    ref_rad_dist
          integer ::    ref_lat
          integer ::    ref_orb_yaw
          integer ::    ref_att_roll
          integer ::    ref_att_pitch
          integer ::    ref_att_yaw
          integer(kind=int1), dimension(8)::  epoch_time  ! BCD_TIME
          integer ::    start_time
          integer ::    IMC_corr_roll
          integer ::    IMC_corr_pitch
          integer ::    IMC_corr_yaw
          integer, dimension(13) ::    ref_long_change
          integer, dimension(11) ::    ref_rad_dist_change
          integer, dimension(9) ::    sine_lat
          integer, dimension(9) ::    sine_orb_yaw
          integer ::    solar_rate
          integer ::    exp_start_time
          type(IMGR_RP) :: roll_att
          integer, dimension(10) ::     spare2
          character(len=4):: more1
          character(len=4):: gvar1
          type(IMGR_RP) :: pitch_att
          type(IMGR_RP) :: yaw_att
          integer, dimension(16) ::    spare3
          character(len=4):: more2
          character(len=4):: gvar2
          type(IMGR_RP) :: roll_misalgn
          type(IMGR_RP) :: pitch_misalgn
          integer ::    img_date
          integer ::    img_time
          integer ::    instr
          integer, dimension(9) ::    spare4
          integer ::    NS_CYL
          integer ::    EW_CYL
          integer ::    NS_INC
          integer ::    EW_INC
          character(len=4):: more3
          character(len=4):: gvar3
          integer, dimension(126) ::     spare5
          character(len=4):: more4
          character(len=4):: gvar4
          integer, dimension(127) ::     spare6
        end type GVAR_NAV
!
!     Calibration table. 
!    
      type, public:: MC_IMGR_CAL_VIS_DET_T
        integer(kind=int4), dimension(8) :: det
      end type MC_IMGR_CAL_VIS_DET_T

      type, public:: MC_IMGR_CAL_IR_DET_T
        integer(kind=int4) :: det0
        integer(kind=int4) :: det2
        integer(kind=int4) :: det4
        integer(kind=int4) :: det6
      end type MC_IMGR_CAL_IR_DET_T

      type, public:: MC_IMGR_CAL_ALL_IR_DET_T 
        type(MC_IMGR_CAL_IR_DET_T):: PRI
        type(MC_IMGR_CAL_IR_DET_T):: SEC
      end type MC_IMGR_CAL_ALL_IR_DET_T

      type, public:: CALIB_STRUCT
	     type(MC_IMGR_CAL_VIS_DET_T):: VIS_BIAS
	     type(MC_IMGR_CAL_VIS_DET_T):: VIS_GAIN1
	     type(MC_IMGR_CAL_VIS_DET_T):: VIS_GAIN2
	  integer(kind=int4) :: vis_albedo
	  type (MC_IMGR_CAL_ALL_IR_DET_T):: IR_SCL_BIAS
	  type (MC_IMGR_CAL_ALL_IR_DET_T):: IR_SCL_GAIN1
	  integer(kind=int4), dimension(8) :: unknown
	  integer(kind=int4), dimension(79) :: z
      end type CALIB_STRUCT


contains

!-----------------------------------------------------------------------------------------------------
! routine to read to AREA and NAVIGATION headers
!------------------------------------------------------------------------------------------------------
 subroutine GET_HEADERS(filename,AREAstr,NAVStr)
   character(len=*), intent(in):: filename
   type(AREA_STRUCT), intent(out):: AREAstr
   type(GVAR_NAV), intent(out):: NAVstr
   integer:: i,swap,recnum,nav_offset
   integer:: a
   character(len=4):: in

   open(unit=1,file=trim(filename),form="unformatted",access="direct",recl=4,status="old",action="read")

   swap = 0
        read (unit=1, rec = 1) AREAstr%area_status
        read (unit=1, rec = 2) AREAstr%version_num
        if (AREAstr%version_num .ne. 4) then
           swap = 1
           print *, "byte swapping is necessary stopping"
           stop
        endif
        read (unit=1, rec = 3) AREAstr%sat_id_num
        read (unit=1, rec = 4) AREAstr%img_date
        read (unit=1, rec = 5) AREAstr%img_time
        read (unit=1, rec = 6) AREAstr%north_bound
        read (unit=1, rec = 7) AREAstr%west_vis_pixel
        read (unit=1, rec = 8) AREAstr%z_coor
        read (unit=1, rec = 9) AREAstr%num_line
        read (unit=1, rec = 10) AREAstr%num_elem
        read (unit=1, rec = 11) AREAstr%bytes_per_pixel
        read (unit=1, rec = 12) AREAstr%line_res
        read (unit=1, rec = 13) AREAstr%elem_res
        read (unit=1, rec = 14) AREAstr%num_chan
        read (unit=1, rec = 15) AREAstr%num_byte_ln_prefix
        read (unit=1, rec = 16) AREAstr%proj_num
        read (unit=1, rec = 17) AREAstr%creation_date
        read (unit=1, rec = 18) AREAstr%creation_time
        read (unit=1, rec = 19) AREAstr%sndr_filter_map
        read (unit=1, rec = 20) AREAstr%img_id_num
        do i=1,4
           read (unit=1, rec = i + 20) AREAstr%id(i)
        enddo
        do i=1,8
           read (unit=1, rec = i + 24) in
           AREAstr%comment = AREAstr%comment // in
        enddo
        read (unit=1, rec = 33) AREAstr%pri_key_calib
        read (unit=1, rec = 34) AREAstr%pri_key_nav
        read (unit=1, rec = 35) AREAstr%sec_key_nav
        read (unit=1, rec = 36) AREAstr%val_code
        do i=1,8
           read (unit=1, rec=i+36) AREAstr%pdl(i)
        enddo
        read (unit=1, rec = 45) AREAstr%band8
        read (unit=1, rec = 46) AREAstr%act_img_date
        read (unit=1, rec = 47) AREAstr%act_img_time
        read (unit=1, rec = 48) AREAstr%act_start_scan
        read (unit=1, rec = 49) AREAstr%len_prefix_doc
        read (unit=1, rec = 50) AREAstr%len_prefix_calib
        read (unit=1, rec = 51) AREAstr%len_prefix_lev
        do i=1,1
           read (unit=1, rec = i + 51) in
           AREAstr%src_type = AREAstr%src_type//in
        enddo
        do i=1,1
           read (unit=1, rec = i + 52) in
           AREAstr%calib_type = AREAstr%calib_type//in
        enddo
        read (unit=1, rec = 54) AREAstr%avg_or_sample
                read (unit=1, rec = 55) AREAstr%poes_signal
        read (unit=1, rec= 56) AREAstr%poes_up_down
        do i=1,4
           read (unit=1, rec = i + 56) in
           AREAstr%orig_src_type = AREAstr%orig_src_type//in
        enddo
        do i=1,4
           read (unit=1, rec = i + 60) AREAstr%reserved(i)
        enddo

!------------------------------------------------------------
! read in navigation data
!------------------------------------------------------------
        nav_offset = AREAstr%sec_key_nav/4
        read (unit=1, rec = nav_offset + 1) NAVstr%nav_type
        read (unit=1, rec = nav_offset + 2) NAVstr%IMC_status
        read (unit=1, rec = nav_offset + 3) NAVstr%spare1(1)
        read (unit=1, rec = nav_offset + 4) NAVstr%spare1(2)
        read (unit=1, rec = nav_offset + 5) NAVstr%stat
        read (unit=1, rec = nav_offset + 6) NAVstr%ref_long
        read (unit=1, rec = nav_offset + 7) NAVstr%ref_rad_dist
        read (unit=1, rec = nav_offset + 8) NAVstr%ref_lat
        read (unit=1, rec = nav_offset + 9) NAVstr%ref_orb_yaw
        read (unit=1, rec = nav_offset + 10) NAVstr%ref_att_roll
        read (unit=1, rec = nav_offset + 11) NAVstr%ref_att_pitch
        read (unit=1, rec = nav_offset + 12) NAVstr%ref_att_yaw
        read (unit=1, rec = nav_offset + 13) a !NAVstr%epoch_time(1:4)
        read (unit=1, rec = nav_offset + 14) a !NAVstr%epoch_time(5:8)
        read (unit=1, rec = nav_offset + 15) NAVstr%start_time
        read (unit=1, rec = nav_offset + 16) NAVstr%imc_corr_roll
        read (unit=1, rec = nav_offset + 17) NAVstr%imc_corr_pitch
        read (unit=1, rec = nav_offset + 18) NAVstr%imc_corr_yaw
        do i=1,13
           read (unit=1, rec = i + nav_offset + 18) NAVstr%ref_long_change(i)
        enddo
        do i=1,11
           read (unit=1, rec = i + nav_offset + 31) NAVstr%ref_rad_dist_change(i)
        enddo
        do i=1,9
           read (unit=1, rec = i + nav_offset + 42) NAVstr%sine_lat(i)
        enddo
        do i=1,9
           read (unit=1, rec = i + nav_offset + 51) NAVstr%sine_orb_yaw(i)
        enddo
        read (unit=1, rec = nav_offset + 61) NAVstr%solar_rate
        read (unit=1, rec = nav_offset + 62) NAVstr%exp_start_time
        recnum = nav_offset + 62
        call READ_IMGR_RP(recnum, NAVstr%roll_att)
        do i=1,10
           recnum = recnum + 1
           read (unit=1, rec = recnum) NAVstr%spare2(i)
        enddo
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%more1
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%gvar1
        call READ_IMGR_RP(recnum, NAVstr%pitch_att)
        call READ_IMGR_RP(recnum, NAVstr%yaw_att)
        do i=1,16
           recnum = recnum + 1
           read (unit=1, rec = recnum) NAVstr%spare3(i)
        enddo
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%more2
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%gvar2
        call READ_IMGR_RP(recnum, NAVstr%roll_misalgn)
        call READ_IMGR_RP(recnum, NAVstr%pitch_misalgn)
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%img_date
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%img_time
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%instr
        do i=1,9
           recnum = recnum + 1
          read (unit=1, rec = recnum) NAVstr%spare4(i)
        enddo
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%ns_cyl
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%ew_cyl
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%ns_inc
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%ew_inc
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%more3
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%gvar3
        do i=1,126
           recnum = recnum + 1
           read (unit=1, rec = recnum) NAVstr%spare5(i)
        enddo
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%more4
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%gvar4
        do i=1,127
           recnum = recnum + 1
           read (unit=1, rec = recnum) NAVstr%spare6(i)
        enddo


   close(unit=1)

 end subroutine GET_HEADERS

!-----------------------------------------------------------------------
subroutine GET_IMAGE(filename,imgstart,prefix,nx,ny,image,byte_shift,error_status)
  character(len=*), intent(in):: filename
  integer, intent(in):: imgstart, prefix, nx, ny, byte_shift
  integer(kind=int2), dimension(:,:), intent(out):: image
  integer, intent(out):: error_status
  integer:: iline, mypos
  integer:: ielem_start, i
  integer(kind=int2), dimension(nx):: i2_buffer

  integer(kind=int2):: i2_dummy
  integer(kind=int4):: i4_dummy


  error_status = 0

  !--- open area file
  open (unit=1,file=filename,form="unformatted",access="stream", &
        status="old",action="read",iostat=error_status)

  if (error_status /= 0) return

  do iline=1,ny !loop over lines

    !--- determine start of data for this line (in bytes)
    ielem_start = imgstart + prefix + (prefix + nx*2)*(iline-1) + 1

!   print *, 'starting pos = ',ielem_start

    !--- read in one line of data
    read (unit=1,pos = ielem_start, iostat=error_status) i2_buffer

    if (error_status /= 0) then
      close(unit=1)
      return
    endif

!   inquire(unit=1,POS=mypos)
!   print *, "mypos after read = ", mypos

    i2_buffer = ishft(i2_buffer,byte_shift)

    image(:,iline) = i2_buffer

  enddo

  !--- close area file
  close(unit=1)

 end subroutine GET_IMAGE


!-----------------------------------------------------------------------
 subroutine READ_IMGR_RP(recnum, struct)
        integer, intent(inout):: recnum
        type(IMGR_RP), intent(out):: struct
        integer:: i

        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%exp_mag
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%exp_time_const
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%mean_att_ang_const
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%num_sinu_per_angle
        do i=1,15
           call READ_IMGR_SIN(recnum, struct%sinusoid(i))
        enddo
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%num_mono_sinu
        do i=1,4
           call READ_IMGR_MON(recnum, struct%monomial(i))
        enddo
 end subroutine READ_IMGR_RP

!-----------------------------------------------------------------------
  subroutine READ_IMGR_SIN(recnum, struct)
    integer, intent(inout):: recnum
    type(imgr_sin), intent(out):: struct
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%mag_sinu
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%phase_ang_sinu
  end subroutine READ_IMGR_SIN
!-----------------------------------------------------------------------
  subroutine READ_IMGR_MON(recnum, struct)
    integer, intent(inout):: recnum
     type(imgr_mon), intent(out):: struct
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%order_appl_sinu
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%order_mono_sinu
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%mag_mono_sinu
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%phase_ang_sinu
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%ang_from_epoch
   end subroutine READ_IMGR_MON

!-----------------------------------------------------------------------
        subroutine FIX_VAR(var)
        integer, intent(inout):: var

        var=ishftc(ishftc(ishftc(var,-8,32),-8,24),-8,16)

        end subroutine FIX_VAR
!-----------------------------------------------------------------------
        subroutine FIX_VAR2(var)
!       integer(kind=int2), intent(inout):: var
        integer*2, intent(inout):: var

        var=ishftc(var,-8,16)

        end subroutine FIX_VAR2

end module GOES_AREA_IO_MODULE
