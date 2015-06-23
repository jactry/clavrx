!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: geo_dark_sky_vis.f90 (src)
!       GEO_DARK_SKY_VIS (program)
!
! PURPOSE: a code to  make ch1 dakr-sky composites for geostationary data
!
! DESCRIPTION: reads in AREA files in raw format and does a simplistic 
!              renavigation based on a reference line and element origin
!
!              this assumes that the image may shift around but always has
!              the same number of pixels (nx,ny)
!
!              this version takes the values of beg_sc and beg_el from the
!              first image and uses them as the reference for subsequent images.
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
! REVISION HISTORY:
!   2000: Andrew K Heidinger, NOAA/NESDIS Created
!   2006: modified to use stream i/o for faster reading
!   2012: Andrew K Heidinger included MTSAT and SEVIRI
!   2013: Incorporated in the clavrx src and modified to 
!         use clavrx routines instead of separate ones

!  Geostationary WMO ID's
!
!select case (Sc_Id_WMO)
!   case(252)
!        satname = "goes08"
!   case(253)
!        satname = "goes09"
!   case(254)
!        satname = "goes10"
!   case(255)
!        satname = "goes11"
!   case(256)
!        satname = "goes12"
!   case(257)
!        satname = "goes13"
!   case(258)
!        satname = "goes14"
!   case(259)
!        satname = "goes15"
!   case(55)
!        satname = "met8"
!   case(56)
!        satname = "met9"
!   case(57)
!        satname = "met10"
!   case(171)
!        satname = "mtsat-1r"
!   case(172)
!        satname = "mtsat-2"
!   case(810)
!        satname = "coms-1"
!   case default
!        return
!end select
!--------------------------------------------------------------------------------------
program GEO_DARK_SKY_VIS
  use CONSTANTS
  use GOES_MODULE

  implicit none

  type(AREA_STRUCT)::  AREAstr
  type(GVAR_NAV)::     NAVstr
  type(AREA_STRUCT)::  AREAstr_ref
  type(GVAR_NAV)::   NAVstr_ref


  integer(kind=int2), allocatable, dimension(:,:) :: image,temp_image
  integer(kind=int2), allocatable, dimension(:,:) :: dark_image
  integer(kind=int2), allocatable, dimension(:,:) :: second_dark_image
  integer(kind=int2), allocatable, dimension(:) :: tempx
  integer(kind=int2), allocatable, dimension(:) :: tempy
  integer(kind=int1):: space_mask
  integer(kind=int4):: wmo_id
  integer:: byte_shift
  integer:: nx_expected
  integer:: ny_expected_min
  integer:: ny_expected_max
  integer:: nx, ny, ny_read, nx_ref,ny_ref,imgstart,prefix , irow,icol
  integer:: i,j,beg_sc, beg_el, el_res, ln_res, length,ierr,jj,ii,j1,j2,irec, &
                  iscan,file_number,ios,beg_sc_ref,beg_el_ref, &
                  nx_offset,ny_offset, ilen

  integer:: error_status
  integer:: nday_input
  integer:: nday
  integer:: jday_temp
  integer:: year_temp

  integer:: jday,isat,itime,year,ileap,big_count,space_count,jday_first,jday_last
  character(len=7):: domain
  character(len=8):: satname
  real(kind=real4):: time
  real(kind=real4):: lat,lon,rlat,rlon

  character(len=120):: system_string
  character(len=120):: input_name, output_name, data_path, output_path, &
                      home_path,areafile_compressed,data_path_temp
  character(len=4):: year_string
  character(len=4):: year_temp_string
  character(len=1):: chan_string
  character(len=3):: jday_string
  character(len=3):: jday_temp_string
  character(len=3):: nday_string
  character(len=4):: time_string
  character(len=3):: comp_string
  character(len=10):: extra_path

  integer, parameter:: file_number_thresh = 1

 !----------------------------------------------------------------
 ! set reference values of beg_sc, and beg_el
 !----------------------------------------------------------------
 big_count = 1000
 space_count = 25

 !----- read in gsip startup file
 home_path = "./"
 open(unit=8,file=trim(home_path)//"dark_ch1_start_file",status="old",action="read")
 read(unit=8,fmt="(a)") data_path
 read(unit=8,fmt="(a)") output_path
 read(unit=8,fmt=*) year
 read(unit=8,fmt=*) jday
 read(unit=8,fmt=*) nday_input              !number of days to use in composite
 read(unit=8,fmt=*) itime
 read(unit=8,fmt="(a)") domain
 read(unit=8,fmt=*) wmo_id 

 !---- if nday_input is less than zero, 
 nday = nday_input
 if (nday_input < 0) then
     nday = -1*nday_input
 endif

 !----  based on GOES satellite id, choose expected satellite position
 select case (Sc_Id_WMO)
    case(252)
      satname = "goes08"
    case(253)
      satname = "goes09"
    case(254)
      satname = "goes10"
    case(255)
      satname = "goes11"
    case(256)
      satname = "goes12"
    case(257)
      satname = "goes13"
    case(258)
      satname = "goes14"
    case(259)
      satname = "goes15"
    case(55)
      satname = "met8"
    case(56)
      satname = "met9"
    case(57)
      satname = "met10"
    case(171)
      satname = "mtsat-1r"
    case(172)
      satname = "mtsat-2"
    case(810)
      satname = "coms-1"
    case default
      print *,"Unable to process this GOES satellite, stopping"
      stop 
 end select

 select case (Sc_Id_WMO)
    case(252:259)
      space_count = 25
      byte_shift = -5
    case default
      space_count = 0
      byte_shift = 0
 end select



 !----  based on GOES satellite id, choose expected satellite position
 select case (trim(satname))
    case("goes08","goes12","goes13")
      ny_conus =  883
      ny_nhem =  1827
      ny_fulldisk = 2705 
    case("goes09","goes10","goes11","goes15")
      ny_conus =  983
      ny_nhem = 1355 
      ny_fulldisk = 2705 
    case("met8","met9","met10")
      ny_fulldisk = 3712
    case("mtsat-1r","mtsat-2")
      ny_nhem = 1375 
      ny_fulldisk = 2750
    case("coms-1")
      ny_nhem = 1544
      ny_fulldisk = 2750
 end select


 !--- based on julian day and year, derive other date values

 !--- set leap year flag
 ileap = 0
 if (mod(year,4) == 0) ileap = 1   
 if (mod(year,100) == 0) ileap = 0 
 if (mod(year,400) == 0) ileap = 1

 !--- make some strings
 write (year_string,  '(I4.4)') year
 write (time_string,  '(I4.4)') itime
 write (jday_string,  '(I3.3)') jday
 write (nday_string,  '(I3.3)') nday

 !--- compute image time
 time = itime / 100
 time = time + (itime - time*100)/60.0

 !----------- compute range of days
 if (nday_input >= 0) then 
     jday_first =jday - nday
     jday_last = jday 
   else
     jday_first =jday - nday
     jday_last = jday + nday
 endif

 !--- add year to output path
 output_path = trim(output_path)//trim(satname)//"/"//trim(year_string)//"/"


!   if (nday_input >= 0) then 
!      output_name = trim(satname)//"_"//trim(domain)//"_"//trim(year_string)//"_"// &
!                    trim(jday_string)//"_"//trim(nday_string)//"_"//trim(time_string)//"_drk_ch1_pix.dat"
!   else
!      output_name = trim(satname)//"_"//trim(domain)//"_"//trim(year_string)//"_"// &
!                    trim(jday_string)//"_pm"//trim(nday_string)//"_"//trim(time_string)//"_drk_ch1_pix.dat"
!   endif

  output_name = trim(satname)//"_"//trim(year_string)//"_"// &
                  trim(jday_string)//"_"//trim(time_string)//"_drk_ch1_pix.dat"

!------------------------------------------------------------
! make a list of files to be composited
!------------------------------------------------------------
  jday_temp = jday
  year_temp = year
  write (jday_temp_string,  '(I3.3)') jday_temp
  write (year_temp_string,  '(I4.4)') year_temp

  !--- for data as stored on saga
  extra_path = trim(year_temp_string)//"_"//trim(jday_temp_string)//"/"

  !---- tiros
  extra_path = "/"


  data_path_temp  = trim(data_path)//trim(extra_path)


  system_string =  "ls "//trim(data_path_temp)//"*_1_*"// &
                   time_string//".area* | grep "// &
                   year_temp_string//"_"//jday_temp_string//" > ch1_list"

  call system(system_string)

!------------------------------------------------------------
! loop over other days make a list of files to be composited
!------------------------------------------------------------
do j = jday_first,jday_last

    if (j == jday) cycle   !skip the actual day (it is first)

    jday_temp_string = char(jday_temp/100 + 48) //  &
                       char(jday_temp/10 - 10*(jday_temp/100) + 48) //  &
                       char(jday_temp - 10*(jday_temp/10) +48)


    jday_temp = j
    year_temp = year

    if (j <= 0) then
       jday_temp = 365 + j
       year_temp = year - 1
    endif

    if (j > 365) then
       jday_temp = j - 365
       year_temp = year + 1
    endif


    jday_temp_string = char(jday_temp/100 + 48) //  &
                       char(jday_temp/10 - 10*(jday_temp/100) + 48) //  &
                       char(jday_temp - 10*(jday_temp/10) +48)

    year_temp_string = char(year_temp/1000 + 48) //  &
              char(year_temp/100 - 10*(year_temp/1000) + 48) //  &
              char(year_temp/10 - 10*(year_temp/100) + 48) //  &
              char(year_temp - 10*(year_temp/10) +48)

    !---- this works for the real-time tiros archive - note only selects compressed files
    data_path_temp  = trim(data_path)//trim(extra_path)
    system_string =  "ls "//trim(data_path_temp)//"*_1_*"// &
                     time_string//".area* | grep "// &
                     year_temp_string//"_"//jday_temp_string//" >> ch1_list"

    !--- make file list 
    call system(system_string)

enddo

!--- determine expected size
if (domain == "nhem") then
   ny_expected_min = 1247
   ny_expected_max = 1827
endif
if (domain == "disk") then
   !---stw ny_expected_min = 2705
   !---stw ny_expected_max = 2707
   ny_expected_min = 2705
   ny_expected_max = 2750
endif
if (domain == "conus") then 
   ny_expected_min = 00
   ny_expected_max = 1000
endif
if (domain == "east1km") then 
   ny_expected_min = 1850
   ny_expected_max = 1950
endif
if (ssec_id == 84 .or. ssec_id == 85) then
   ny_expected_min = 2725
   ny_expected_max = 2775
endif
if (ssec_id == 51 .or. ssec_id == 52 .or. ssec_id == 53 ) then
   ny_expected_min = 3710
   ny_expected_max = 3715
endif

!--------------- loop over files in list
open(unit=8,file="ch1_list",form="formatted",status="old",action="read")
file_number = 1
file_loop: do

   read(unit=8,fmt="(a)", iostat = ios) areafile_compressed

   if (ios /= 0) then
    exit
   endif

!-----------------------------------------------------------
! uncompress file
!-----------------------------------------------------------
ilen = len_trim(areafile_compressed)
!areafile_compressed = adjustr(areafile_compressed)
!ilen = len(areafile_compressed)
comp_string = areafile_compressed(ilen-2:ilen+1)


if (comp_string == ".gz") then 
  system_string =  "gunzip -c "//trim(areafile_compressed)//" > temp_dark_area_file"
elseif (comp_string == "bz2") then 
  system_string =  "bunzip2 -c "//trim(areafile_compressed)//" > temp_dark_area_file"
else
  system_string =  "cp "//trim(areafile_compressed)//"  temp_dark_area_file"
endif

call system(system_string)

!print *, "dark_ch1 processing this file: ", areafile_compressed

!------------------------------------------------------------
! read in header derive some needed parameters
!------------------------------------------------------------
call GET_GOES_HEADERS("temp_dark_area_file",AREAstr,NAVstr)
nx = AREAstr.num_elem
ny = AREAstr.num_line

!--- skip this image if it is not the correct size - only check y
if (ny < ny_expected_min .or. ny > ny_expected_max) then
  print *, "skipping because image size differs from expectations ", ny, ny_expected_min,ny_expected_max
  cycle
endif


imgstart = AREAstr.pri_key_nav
prefix = AREAstr.num_byte_ln_prefix
beg_sc = AREAstr.NORTH_BOUND
beg_el = AREAstr.WEST_VIS_PIXEL
el_res = AREAstr.elem_res
ln_res = AREAstr.line_res

!------------------------------------------------------------
! take navigation and reference point from first file
!------------------------------------------------------------
if (file_number == 1) then
   beg_sc_ref = beg_sc
   beg_el_ref = beg_el
   AREAstr_ref = AREAstr
   NAVstr_ref = NAVstr
   nx_ref = nx
   ny_ref = ny
   allocate (image(nx_ref,ny_ref),dark_image(nx_ref,ny_ref),second_dark_image(nx_ref,ny_ref))
   image = 0
   image = 0
   dark_image = big_count
   second_dark_image = big_count
!  print *,"reference image set"
!  print *,"reference image size = ", nx_ref, ny_ref
!  print *,"reference origin  = ", beg_sc_ref, beg_el_ref
endif
!--------------------------------------------------
! allocate memory for current image
!---------------------------------------------------
allocate(temp_image(nx,ny),tempx(nx),tempy(ny))
temp_image = 0
tempx = 0
tempy = 0

call GET_IMAGE_FROM_AREAFILE("temp_dark_area_file", &
                             byte_shift, &
                             AREAstr, &
                             1, &
                             1, &
                             ny, &
                             ny_read, &
                             temp_image(:,:))

if (error_status /= 0) then
   print *, "skipping file because of read error"
   if (allocated(temp_image)) deallocate(temp_image)
   if (allocated(tempx)) deallocate(tempx)
   if (allocated(tempy)) deallocate(tempy)
   cycle
endif

!------------ east/west offsets
if (beg_el /= beg_el_ref) then

nx_offset = abs(beg_el - beg_el_ref)/el_res

if (nx_offset < 0) then
       print *,"shifting to west ", nx_offset
else
       print *,"shifting to east", nx_offset
endif

do j = 1, ny
   tempx = temp_image(:,j)
   if (beg_el < beg_el_ref) then    !this is a shift to the west
       temp_image(1:nx - nx_offset, j) = tempx(nx_offset+1:nx)
       temp_image(nx-nx_offset+1:nx,j) = big_count
   endif
   if (beg_el > beg_el_ref) then    !this is a shift to the east
       temp_image(1:nx_offset, j) = big_count
       temp_image(nx_offset+1:nx,j) = tempx(1:nx - nx_offset)
   endif
enddo
endif

!------- north/south offsets
 if (beg_sc /= beg_sc_ref) then
   ny_offset = abs(beg_sc - beg_sc_ref)/ln_res
   do i = 1, nx
    tempy = temp_image(i,:)
    if (beg_sc < beg_sc_ref) then    !this is a shift to the north
       temp_image(i,1:ny - ny_offset) = tempy(ny_offset+1:ny)
       temp_image(i,ny-ny_offset+1:ny) = big_count
    endif
    if (beg_sc > beg_sc_ref) then    !this is a shift to the south
       temp_image(i,1:ny_offset) = big_count
       temp_image(i,ny_offset+1:ny) = tempy(1:ny - ny_offset)
    endif
   enddo
  endif 
!--------------------------------------------------------------------
! trim or pad ccr image to match current data
!--------------------------------------------------------------------
   image = 0
   image(1:min(nx,nx_ref),1:min(ny,ny_ref)) = &
      temp_image(1:min(nx,nx_ref),1:min(ny,ny_ref))
   deallocate(temp_image,tempx,tempy)

!-------------------------------------------------------------------
! construct dark composite
!--------------------------------------------------------------------

!---update second darkest when darkest is updated
!where (image < dark_image .and. image > space_count)
!   second_dark_image = image
!   dark_image = image
!endwhere

 do i = 1, nx_ref
   do j = 1, ny_ref
     if ((image(i,j) < dark_image(i,j)).and.(image(i,j) > space_count)) then
       second_dark_image(i,j) = image(i,j)
       dark_image(i,j) = image(i,j)
     endif
   enddo
 enddo

!---update second darkest when darkest is not updated
!where (image > dark_image .and. image < second_dark_image .and. image > space_count)
!   second_dark_image = image
!endwhere

 do i = 1, nx_ref
   do j = 1, ny_ref
     if ((image(i,j) > dark_image(i,j)) .and. &
         (image(i,j) < second_dark_image(i,j)) .and. &
         (image(i,j) > space_count)) then
         second_dark_image(i,j) = image(i,j)
     endif
   enddo
 enddo

file_number = file_number + 1

enddo file_loop

print *, "Dark Comp Includes ", file_number, " files"

if (file_number > file_number_thresh) then

!------------------------- reset remaining high values to 0
!where (dark_image > 999)
!   dark_image = 0
!endwhere
 do i = 1, nx_ref
   do j = 1, ny_ref
    if (dark_image(i,j) > 999) then
     dark_image(i,j) = 0
    endif
   end do
 end do


!where (second_dark_image > 999)
!   second_dark_image = 0
!endwhere

 do i = 1, nx_ref
  do j = 1, ny_ref
   if (second_dark_image(i,j) > 999) then
    second_dark_image(i,j) = 0
   endif
  end do
 end do

!------------------------------------------------------
! output
!------------------------------------------------------
print *,"writing out to file: ", trim(output_name), ' ', nx_ref, ny_ref
open(unit=11,file=trim(output_path)//trim(output_name), &
              access="direct",recl=2*nx_ref, &
              action="write", status="replace",form="unformatted")
write(unit=11,rec=1) AREAstr_ref
write(unit=11,rec=2) NAVstr_ref
do j = 1, ny_ref
   write(unit=11, rec=j+2) second_dark_image(:,j)
enddo
close(unit=11)

endif

if (allocated(image)) deallocate(image)
if (allocated(second_dark_image)) deallocate(second_dark_image)
if (allocated(dark_image)) deallocate(dark_image)

end program GEO_DARK_SKY_VIS
