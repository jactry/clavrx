!$Id: batch_dark_ch1.f90,v 1.3 2013/11/12 16:28:01 stevew Exp $
!------------------------------------------------------------
! a fortran! code to run gsip dark ch1 compositing multiple times
!
!------------------------------------------------------------
 program batch_gsip

  implicit none

  character(len=80):: data_path,dark_ch1_path,avn_path
  character(len=7):: domain
  character(len=6):: sat_id
  integer:: ssec_id
  integer:: year,hour,minute,ix_skip,ios
  integer:: jday_first
  integer:: jday_last

  integer:: jday, nday
  integer:: first_hour
  integer:: first_minute
  integer:: last_hour
  integer:: last_minute
  integer:: minute_spacing
  integer:: first_itime
  integer:: last_itime
  integer:: itime_delta
  integer:: itime


!-------------------------------------------------
! read in input from batch_dark_ch1_input
!-------------------------------------------------
open(unit=7,file="batch_dark_ch1_input",status="old",action="read")
read(unit=7,fmt="(a)") data_path
read(unit=7,fmt="(a)") dark_ch1_path
read(unit=7,fmt=*) year
read(unit=7,fmt=*) first_hour
read(unit=7,fmt=*) first_minute
read(unit=7,fmt=*) last_hour
read(unit=7,fmt=*) last_minute
read(unit=7,fmt=*) minute_spacing
read(unit=7,fmt="(a7)") domain
read(unit=7,fmt="(a6)") sat_id 
read(unit=7,fmt=*) nday
read(unit=7,fmt=*) jday_first
read(unit=7,fmt=*) jday_last
close(unit=7)

if (trim(sat_id) == "goes08") ssec_id = 70
if (trim(sat_id) == "goes09") ssec_id = 72
if (trim(sat_id) == "goes10") ssec_id = 74
if (trim(sat_id) == "goes11") ssec_id = 76
if (trim(sat_id) == "goes12") ssec_id = 78
if (trim(sat_id) == "goes13") ssec_id = 180 
if (trim(sat_id) == "goes14") ssec_id = 182 
if (trim(sat_id) == "goes15") ssec_id = 184 
if (trim(sat_id) == "mtsat-1r") ssec_id = 84
if (trim(sat_id) == "mtsat-2") ssec_id = 85
if (trim(sat_id) == "met8") ssec_id = 51
if (trim(sat_id) == "met9") ssec_id = 52
if (trim(sat_id) == "coms-1") ssec_id = 250

ios = 0
first_itime = first_hour*100 + first_minute
last_itime = last_hour*100 + last_minute

jday_loop: do jday = jday_first, jday_last

itime = -999
  do

  if (itime == -999) then
     itime = first_itime
  else
   hour = itime/100
   minute = itime - hour*100
   minute = minute + minute_spacing
   if (minute >= 180) then 
      minute = minute - 180
      hour = hour + 3
   elseif (minute >= 120) then 
      minute = minute - 120
      hour = hour + 2
   elseif (minute >= 60) then 
      minute = minute - 60
      hour = hour + 1
   endif
   itime = hour*100 + minute 
  endif

  if (itime > last_itime)  then
    print *, "exiting time loop"
    exit 
  endif
  
!------------------------------------------------------------
! write out startfile
!------------------------------------------------------------
 open(unit=8,file="dark_ch1_start_file",&
      status="replace",action="write")
 write(unit=8,fmt="(a)") data_path
 write(unit=8,fmt="(a)") dark_ch1_path
 write(unit=8,fmt=*) year
 write(unit=8,fmt=*) jday
 write(unit=8,fmt=*) nday
 write(unit=8,fmt=*) itime
 write(unit=8,fmt="(a7)") domain
 write(unit=8,fmt=*) ssec_id
 close(unit=8)

 write (unit=6,fmt=*)


  call system("date")
  print *, "Calling DARK_CH1 for jday = ", jday," and time = ", itime
  call system('./dark_ch1')

  enddo

end do jday_loop

end program batch_gsip 
