! $Id$

module date_tools_mod
   
   implicit none
   integer, parameter :: r15 = selected_real_kind(7)
    
   ! interface operator(+)
   !   procedure add_day
   ! end interface
    
   type :: date_type 
      integer :: year = 1970
      integer :: month = 1
      integer :: day = 1
      integer :: hour = 0
      integer :: minute = 0 
      integer :: second = 0
      integer :: dayOfYear
      real :: hour_frac
      integer :: msec_of_day
      character(6) :: yymmdd
      character(8) :: yymmddhh
      real (kind = r15) :: julday
      real :: mod_julday
           
   contains
      procedure :: set_date 
      
      procedure :: set_date_with_doy
      procedure :: set_jday
      procedure :: print_data 
      procedure :: date_string 
      procedure :: add_time 
      procedure :: get_date
      procedure , private :: update 
      procedure , private :: set_julday
      procedure , private::  add_day
      procedure , private :: update_from_jd
      procedure ::  period_16    
   
   end type date_type
   
  
contains   
   ! --
   !
   ! ---
   subroutine set_date (this , year , month, day, hour ,minute, second)
      class ( date_type) :: this
      integer , optional :: year , month ,day, hour , minute, second

       if ( present(year) ) this % year = year
       if ( present(month) ) this % month = month
       if ( present(day) ) this % day = day
       if ( present ( hour ) ) this % hour = hour
       if ( present ( minute ) ) this % minute = minute
       if ( present ( second ) ) this % second = second
  
       call this % update()
       
   end subroutine set_date
   

   ! ----------------------------
   !
   ! ----------------------------
   subroutine get_date ( this , year , month, day, hour ,minute, second, doy , hour_frac, msec_of_day)
      class ( date_type) :: this
      integer , optional :: year , month ,day, hour , minute, second, doy
      real, optional :: hour_frac
      integer, optional :: msec_of_day
      
      if ( present(year) )  year = this % year
      if ( present(month) ) month = this % month 
      if ( present(day) ) day = this % day 
      if ( present ( hour ) ) hour = this % hour
      if ( present ( minute ) ) minute =  this % minute 
      if ( present ( second ) ) second = this % second 
      if ( present ( doy ) ) doy = this % dayOfYear
      if ( present ( hour_frac )) hour_frac = this % hour_frac
      if ( present ( msec_of_day )) msec_of_day = this % msec_of_day
   
   end subroutine get_date
   
   
   
   
   
   !
   !
   !
   function period_16 (this ) result ( out_string )
      class ( date_type) :: this
      character ( len = 3 ) :: out_string
      integer :: iperiod16
      
       iperiod16 = 16 * (( this % dayofyear-1) / 16  ) + 1
     
      write (  out_string  , fmt = '(i3.3)') iperiod16
   
   end function period_16
   !
   !
   !
   
   function first_doy_month ( this ) result ( out_string)
      type ( date_type ) :: this
      character ( len = 3 ) :: out_string
      character(len=3), dimension(12) :: first_doy
     
     
      first_doy = ['001','032','060','091','121','152','182','213','244','274','305','335']
      out_string = first_doy ( this % month)
   
   end function  first_doy_month

   !
   !
   !
   subroutine set_date_with_doy (this , year , doy, hour , minute, second)
      class ( date_type) :: this
      integer , optional :: year , doy , hour , minute, second
      integer , dimension(12) :: jmonth 
      integer, dimension(12) :: last_day_month
      integer :: i
      integer , dimension(12) :: day_dum
      integer :: month(1)
       
       if ( present(year) ) this % year = year
      
      if ( present(doy) ) this % dayOfYear = doy
      if ( present ( hour ) ) this % hour = hour
      if ( present ( minute ) ) this % minute = minute
      if ( present ( second ) ) this % second = second
       
      jmonth = [31,28,31,30,31,30,31,31,30,31,30,31]
      if  ( modulo ( this % year , 4) == 0 ) jmonth(2) = 29 
      last_day_month(1) = 31
       
      do i = 2 , 12
         last_day_month (i) = last_day_month (i-1) + jmonth(i)  
      end do
      day_dum = doy - last_day_month
      month = maxloc ( day_dum , mask = day_dum <= 0 )
      
      this % month = month(1)
      this % day = jmonth(month(1)) + day_dum(month(1))
      
      call this % update()
       
   end subroutine set_date_with_doy
   
   ! --
   
   subroutine set_jday ( this, jday)
      class ( date_type) :: this
      integer ( kind = r15 ) :: jday
      this % julday = jday
      call this % update_from_jd()
   end subroutine set_jday
   !
   !
   !
   function add_day ( this, day_to_add) result ( out )
      class ( date_type)  :: this
      type ( date_type ) :: out
  
      integer , intent(in):: day_to_add
      
      select type (this)
         type is (date_type)
             out = this
      
      end select
      out % day = out % day + day_to_add
      call out % update ()
      
   end function add_day
   !
   !
   !
   function next_6h ( date_str , count ) result ( out)
      type ( date_type) , intent (in) :: date_str
      type (date_type)  :: out
      integer , intent(in), optional :: count
      integer :: cnt
      cnt = 1
      if ( present ( count )) cnt = count
      if (cnt <= 0) cnt = cnt + 1
      
      out = date_str
      
      out % hour = 6 * ((date_str % hour )/ 6 + cnt )
      out % minute = 0
      out % second = 0
      call out % update()
   
   end function next_6h
   
   
  
   
   ! --
   subroutine update ( this )
      class ( date_type ) :: this
        character ( len = 2 ) :: year_s2d
      character ( len = 4 ) :: year_s
      character ( len = 2 ) :: month_s
      character ( len = 2 ) :: day_s
      character ( len = 2 ) :: hour_s
      character ( len = 2 ) :: minute_s
      
      this % julday = jd (this % year , this % month &
                           , this % day ,this % hour &
                           , this % minute )
       this % dayofyear = 1+  this % julday -  jd (this % year , 1 &
                           , 1 ,this % hour &
                           , this % minute )                   
      call this % update_from_jd()
      
      write ( year_s2d, fmt ='(i2.2)') mod(this % year , 100)
      write ( year_s, fmt = '(i4.4)') this % year
      write ( month_s, fmt = '(i2.2)') this % month
      write ( day_s, fmt = '(i2.2)') this % day
      write ( hour_s , fmt = '(i2.2)') this % hour
      write ( minute_s , fmt = '(i2.2)') this % minute
      this % yymmdd =  year_s2d//month_s//day_s
      this % yymmddhh =  year_s2d//month_s//day_s//hour_s
   end subroutine update
   
   
  subroutine update_from_jd ( this )
      class ( date_type ) :: this
      call cdate ( this % julday , this % year , this % month , this % day ) 
      this % hour =  int (24.* ( this % julday - int(this % julday) ))
      this % minute = int((60)  &
             * ((24. * ( this % julday - int(this % julday) )) - this % hour))
      this % hour_frac = this % hour + this % minute / 60.
     
      this % msec_of_day =  60* 60* 1000 * this % hour &
                           + 60* 1000 * this % minute &
                           + 1000 * this % second
                           
  
   end subroutine update_from_jd
   
   
   subroutine set_julday ( this , julday)
      class ( date_type ) :: this
      real ( kind= r15 ) :: julday
      call cdate ( julday , this % year , this % month , this % day ) 
      this % hour =  int (24.* ( julday - int(julday) ))
      this % minute = int((60)  &
             * ((24. * ( julday - int(julday) )) - this % hour))
      call this % update()
   
   end subroutine set_julday
   
   
   ! -- 
   subroutine print_data( this )
      class (date_type) :: this
      
      print*,'year: ', this % year
      print*,'month: ', this % month
      print*,'day: ', this % day
      print*,'hour: ', this % hour
      print*,'minute: ', this % minute
      print*,'second: ', this % second
      print*,'jday: ' , this % julday
      print*, 'frac of hour: ',this % hour_frac
   
   end subroutine print_data
   !
   !
   !
   function date_string ( this , fmt ) result(out)
      class ( date_type) :: this
      character ( len = * ) , intent (in) :: fmt 
      character (len=20) , allocatable :: out
      character ( len = 2 ) :: year_s2d
      character ( len = 4 ) :: year_s
      character ( len = 2 ) :: month_s
      character ( len = 2 ) :: day_s
      character ( len = 2 ) :: hour_s
      character ( len = 2 ) :: minute_s
      character ( len =3 ) :: doy_s
      integer :: len_fmt
      
      len_fmt = len(fmt) 
      ! not supported by gfortran 4.7!
      allocate ( character(len = len_fmt ) :: out  )
      
      write ( year_s2d, fmt ='(i2.2)') mod(this % year , 100)
      write ( year_s, fmt = '(i4.4)') this % year
      write ( month_s, fmt = '(i2.2)') this % month
      write ( day_s, fmt = '(i2.2)') this % day
      write ( hour_s , fmt = '(i2.2)') this % hour
      write ( minute_s , fmt = '(i2.2)') this % minute
      write ( doy_s , fmt = '(i3.3)') this % dayofyear
      
      out='start'  
      
      
      select case (fmt)
         case ('yymmdd')
            out = year_s2d//month_s//day_s
         case ('yymmddhhmm')
           out = year_s2d//month_s//day_s//hour_s//minute_s  
          case ('yymmddhh')
            out = year_s2d//month_s//day_s//hour_s     
         case ('yy/mm/dd')
           out = year_s2d//'/'//month_s//'/'//day_s 
         case ('yy/mm/dd/hh')
           out = year_s2d//'/'//month_s//'/'//day_s//'/'//hour_s   
         case ('yy/mm/hhmm')
           out = year_s2d//'/'//month_s//'/'//hour_s//minute_s     
         case ('yyyy')
           out = year_s  
         case ('yyyy_doy')
           out = year_s//'_'//doy_s    
         case default
            out='format not set'
            print*,'WARNING: ',out, '> ',fmt
      end select
      
      
   end function date_string
   !
   !
   !
   subroutine add_time ( this  , day , hour , minute)
      class ( date_type ) :: this
      integer , optional ::   day , hour , minute
      integer ::   day_add, hour_add , minute_add
      real ( kind = r15) :: julday_add
      
      day_add = 0
      hour_add = 0
      minute_add = 0
      
      if ( present ( day) )      day_add = day
      if ( present ( hour ) )    hour_add = hour
      if ( present ( minute ) )  minute_add = minute
      
      julday_add = this % julday + day_add + hour_add/24. + minute_add/(24.*60.)
      call cdate ( julday_add , this % year , this % month , this % day ) 
      this % hour =  int (24.* ( julday_add - int(julday_add) ))
      this % minute = int((60)  &
             * ((24. * ( julday_add - int(julday_add) )) - this % hour))
      call this % update()
   end subroutine add_time
   
   function time_diff_weight ( time, time0, time1 ) result(wgt)
      type (date_type) , intent(in) :: time, time0 , time1
      real (kind = r15) :: wgt
      
      wgt = (time%julday - time0%julday)/ (time1%julday - time0 %julday)
      
   
   end function time_diff_weight
   
   function time_is_in_window ( time, time0, time1 ) 
      type (date_type) , intent(in) :: time, time0 , time1
      logical :: time_is_in_window
      real :: diff
      time_is_in_window = .false.
      
      diff = time_diff_weight ( time, time0, time1 )
      if ( diff .GE. 0. .and.  diff .le. 1) time_is_in_window = .true.
   
   end function time_is_in_window
   
   
   !           arithmetic functions 'izlr' and 'iday' are taken from remark on
   !        algorithm 398, by j. douglas robertson, cacm 15(10):918.
   function iday(yyyy, mm, dd) result(ival)
      !------iday is a companion to calend; given a calendar date, yyyy, mm,
      !           dd, iday is returned as the day of the year.
      !           example: iday(1984, 4, 22) = 113

      integer, intent(in) :: yyyy, mm, dd
      integer             :: ival

      ival = 3055*(mm+2)/100 - (mm+10)/13*2 -91 + (1-(mod(yyyy, 4)+3)/4 +  &
         (mod(yyyy, 100) + 99)/100 - (mod(yyyy, 400)+399)/400)*(mm+10)/13 + dd

      return
   end function iday


   FUNCTION izlr(yyyy, mm, dd) RESULT(ival)
      !------IZLR(YYYY, MM, DD) GIVES THE WEEKDAY NUMBER 0 = SUNDAY, 1 = MONDAY,
      !      ... 6 = SATURDAY.  EXAMPLE: IZLR(1970, 1, 1) = 4 = THURSDAY

      INTEGER, INTENT(IN) :: yyyy, mm, dd
      INTEGER             :: ival

      ival = MOD((13*(mm+10-(mm+10)/13*12)-1)/5 + dd + 77 + 5*(yyyy+(mm-14)/12 -  &
           (yyyy+(mm-14)/12)/100*100)/4 + (yyyy+(mm-14)/12)/400 -  &
           (yyyy+(mm-14)/12)/100*2, 7)

      RETURN
   END FUNCTION izlr


   SUBROUTINE cdate(jd, yyyy, mm, dd)
      !=======GIVEN A JULIAN DAY NUMBER, NNNNNNNN, YYYY,MM,DD ARE RETURNED AS THE
      !              CALENDAR DATE. JD = NNNNNNNN IS THE JULIAN DATE FROM AN EPOCH
      !              IN THE VERY DISTANT PAST.  SEE CACM 1968 11(10):657,
      !              LETTER TO THE EDITOR BY FLIEGEL AND VAN FLANDERN.
      !    EXAMPLE CALL CDATE(2440588, YYYY, MM, DD) RETURNS 1970 1 1 .

      real ( kind = r15), INTENT(IN)   :: jd
      INTEGER, INTENT(OUT)  :: yyyy
      INTEGER, INTENT(OUT)  :: mm
      INTEGER, INTENT(OUT)  :: dd
      INTEGER :: l, n

      l = jd + 68569
      n = 4*l/146097
      l = l - (146097*n + 3)/4
      yyyy = 4000*(l+1)/1461001
      l = l - 1461*yyyy/4 + 31
      mm = 80*l/2447
      dd = l - 2447*mm/80
      l = mm/11
      mm = mm + 2 - 12*l
      yyyy = 100*(n-49) + yyyy + l
      RETURN
   END SUBROUTINE cdate


   SUBROUTINE daysub(jd, yyyy, mm, dd, wd, ddd)
      !========GIVEN JD, A JULIAN DAY # (SEE ASF JD), THIS ROUTINE CALCULATES DD,
      !        THE DAY NUMBER OF THE MONTH; MM, THE MONTH NUMBER; YYYY THE YEAR;
      !        WD THE WEEKDAY NUMBER, AND DDD THE DAY NUMBER OF THE YEAR.

      !   EXAMPLE: CALL DAYSUB(2440588, YYYY, MM, DD, WD, DDD) YIELDS 1970 1 1 4 1.

      real ( kind = r15) , INTENT(IN)   :: jd
      INTEGER, INTENT(OUT)  :: yyyy
      INTEGER, INTENT(OUT)  :: mm
      INTEGER, INTENT(OUT)  :: dd
      INTEGER, INTENT(OUT)  :: wd
      INTEGER, INTENT(OUT)  :: ddd

      CALL cdate(jd, yyyy, mm, dd)
      wd = izlr(yyyy, mm, dd)
      ddd = iday(yyyy, mm, dd)

      RETURN
   END SUBROUTINE daysub


   function jd(yyyy, mm, dd, hh , uu ) result(ival)
     
      integer, intent(in)  :: yyyy
      integer, intent(in)  :: mm
      integer, intent(in)  :: dd
      integer, intent(in) , optional :: hh , uu 
      real ( kind = r15 )   :: ival

      ival = dd - 32075 + 1461*(yyyy+4800+(mm-14)/12)/4 +  &
       367*(mm-2-((mm-14)/12)*12)/12 - 3*((yyyy+4900+(mm-14)/12)/100)/4
      ival = ival + hh/24. + uu/(24.* 60) 
      return
   end function jd


   function ndays(mm1, dd1, yyyy1, mm2, dd2, yyyy2) result(ival)

      integer, intent(in)  :: mm1
      integer, intent(in)  :: dd1
      integer, intent(in)  :: yyyy1
      integer, intent(in)  :: mm2
      integer, intent(in)  :: dd2
      integer, intent(in)  :: yyyy2
      integer              :: ival

      !==============NDAYS IS RETURNED AS THE NUMBER OF DAYS BETWEEN TWO
      !              DATES; THAT IS  MM1/DD1/YYYY1 MINUS MM2/DD2/YYYY2,
      !              WHERE DATEI AND DATEJ HAVE ELEMENTS MM, DD, YYYY.
      !-------NDAYS WILL BE POSITIVE IFF DATE1 IS MORE RECENT THAN DATE2.

      ival = jd(yyyy1, mm1, dd1) - jd(yyyy2, mm2, dd2)

      RETURN
   END FUNCTION ndays

   !-------------------------------------------------
! subroutine JULIAN(iday,imonth,iyear,jday)
! compute julian day
! input:
!         iday - integer day
!         imonth - integer month
!         iyear - integer year (2 or four digits)
!         
! output : jday - julian day
!--------------------------------------------------
 subroutine JULIAN(iday,imonth,iyear,jday)

!-- Computes julian day (1-365/366)
        integer, intent(in)::  iday,imonth,iyear
        integer, intent(out):: jday
        integer::  j
        integer, dimension(12)::  jmonth

        jmonth = reshape ((/31,28,31,30,31,30,31,31,30,31,30,31/),(/12/))

        jday = iday
        if (modulo(iyear,4) == 0) then
            jmonth(2)=29
        endif

        do j = 1,imonth-1
           jday = jday + jmonth(j)
        end do

   end subroutine JULIAN

!--------------------------------------------
! compute the month
!---------------------------------------------
 function COMPUTE_MONTH(jday,ileap) result(month)
   integer, intent(in):: ileap
   integer, intent(in):: jday
   integer:: month

   month = 0
   if (jday < 32) then
     month = 1
   elseif (jday < 60+ileap) then
     month = 2
   elseif (jday < 91+ileap) then
     month = 3
   elseif (jday < 121+ileap) then
     month = 4
   elseif (jday < 152+ileap) then
     month = 5
   elseif (jday < 182+ileap) then
     month = 6
   elseif (jday < 213+ileap) then
     month = 7
   elseif (jday < 244+ileap) then
     month = 8
   elseif (jday < 274+ileap) then
     month = 9
   elseif (jday < 305+ileap) then
     month = 10
   elseif (jday < 335+ileap) then
     month = 11
   else
     month = 12
   endif

 end function COMPUTE_MONTH
!--------------------------------------------
! compute the day
!---------------------------------------------
 function COMPUTE_DAY(jday,ileap) result(day)
   integer, intent(in):: ileap
   integer, intent(in):: jday
   integer:: day

   if (jday < 32) then
     day = jday
   elseif (jday < 60) then
     day = jday - 31
   elseif ((jday == 60).and.(ileap == 1)) then
     day = jday - 31
   elseif ((jday == 60).and.(ileap == 0)) then
     day = jday - 59
   elseif (jday < 91+ileap) then
     day = jday - (59 + ileap)
   elseif (jday < 121+ileap) then
     day = jday - (90 + ileap)
   elseif (jday < 152+ileap) then
     day = jday - (120 + ileap)
   elseif (jday < 182+ileap) then
     day = jday - (151 + ileap)
   elseif (jday < 213+ileap) then
     day = jday - (181 + ileap)
   elseif (jday < 244+ileap) then
     day = jday - (212 + ileap)
   elseif (jday < 274+ileap) then
     day = jday - (243 + ileap)
   elseif (jday < 305+ileap) then
     day = jday - (273 + ileap)
   elseif (jday < 335+ileap) then
     day = jday - (304 + ileap)
   else
     day = jday - (334 + ileap)
   endif

   return

 end function COMPUTE_DAY

 !---- function to return the system in fractional hours
 function COMPUTE_TIME_HOURS()   result(time_hours)
   character(len=8):: system_date
   character(len=10):: system_time
   character(len=5):: system_time_zone
   integer, dimension(8):: system_time_value

   real:: time_hours

   call DATE_AND_TIME(system_date,system_time,system_time_zone, system_time_value)

   time_hours = real(system_time_value(5)) +  &
                     (real(system_time_value(6)) + &
                      real(system_time_value(7) + &
                      real(system_time_value(8))/1000.0)/60.0)/60.0
   return

 end function COMPUTE_TIME_HOURS
 
 !---------------------------------------------------------------------
! function leap_year_fct(yyyy) result(leap_flg)
!
! Function to determine if an input year is a leap year.
!---------------------------------------------------------------------

function leap_year_fct(yyyy) result(leap_flg)
  integer, intent(in) :: yyyy
  integer :: leap_flg

  leap_flg = 0

  if ((modulo(yyyy,4) == 0 .and. modulo(yyyy,100) /= 0) .or. &
       modulo(yyyy,400) == 0) leap_flg = 1

  return

end function leap_year_fct


   
   
end module date_tools_mod

!program test_it
!  use date_tools_mod
!   type ( date_type) :: t , q
!   character ( len = 26 ) :: pp
!   character ( len =4) :: yy
   
   !
  ! call t % set_date ( year=2003,month=12, hour=19 , minute = 23, second=33)
  ! call t % print_data ()
   !
  ! print*, t % date_string('yyyy')
  ! yy = t % date_string('yyyy') 
 !  print*,yy//'pp'
   
 !  print*, t % date_string('yymmdd')
  ! print*, t % date_string('yy/mm/dd')
  ! print*, t % date_string('yy/mm/hhmm')
  ! print*, jd ( 2003,12,31, 22,21)
   !call t % add_time (hour = 99 )
   
  ! call t % print_data ()
  ! call t % add_time ( minute = 47)
  ! call t % print_data()
   
  ! print*,'======='
  ! call t % set_date ( year=1900,month=2,day=28, hour=19 , minute = 23, second=33)
   
  ! q = t % next_6h()
  ! call q % print_data()
   
  
  ! call q % print_data()
  ! print*,t % date_string('yy/mm/dd')
   
  ! print*,'==========='
   
  ! print*,'test'
   
  ! q = t % next_6h()
  ! call q % print_data
  ! print*, q % date_string ( 'yy/mm/dd')
   
  ! q = t % next_6h(-1)
  ! call q % print_data
  ! print*, '-1:',q % date_string ( 'yy/mm/dd')
   
  ! q = t % next_6h(2)
  ! call q % print_data
  ! print*, '2:',q % date_string ( 'yy/mm/dd')
   
   
  ! print*,'========'
   
   
  ! call t % set_jday ( 2499993 )
  ! call t % print_data

!end program test_it
