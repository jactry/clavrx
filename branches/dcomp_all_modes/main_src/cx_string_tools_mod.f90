! $Id$
!>
!!
!!
module cx_string_tools_mod

   use precision

   private :: value_dr,value_sr,value_di,value_si
   private :: write_dr,write_sr,write_di,write_si
   private :: writeq_dr,writeq_sr,writeq_di,writeq_si

   interface value  ! Generic operator for converting a number string to a 
                 ! number. Calling syntax is 'call value(numstring,number,ios)' 
                 ! where 'numstring' is a number string and 'number' is a 
                 ! real number or an integer (single or double precision).         
      module procedure value_dr
      module procedure value_sr
      module procedure value_di
      module procedure value_si
   end interface

   interface writenum  ! Generic  interface for writing a number to a string. The 
                    ! number is left justified in the string. The calling syntax
                    ! is 'call writenum(number,string,format)' where 'number' is
                    ! a real number or an integer, 'string' is a character string
                    ! containing the result, and 'format' is the format desired, 
                    ! e.g., 'e15.6' or 'i5'.
      module procedure write_dr
      module procedure write_sr
      module procedure write_di
      module procedure write_si
   end interface

   interface writeq  ! Generic interface equating a name to a numerical value. The
                  ! calling syntax is 'call writeq(unit,name,value,format)' where
                  ! unit is the integer output unit number, 'name' is the variable
                  ! name, 'value' is the real or integer value of the variable, 
                  ! and 'format' is the format of the value. The result written to
                  ! the output unit has the form <name> = <value>.
   module procedure writeq_dr
   module procedure writeq_sr
   module procedure writeq_di
   module procedure writeq_si
end interface

INTERFACE Copy    ! generic
   MODULE PROCEDURE copy_a2s, copy_s2a
END INTERFACE Copy


!**********************************************************************

contains

!**********************************************************************

subroutine parse(str,delims,args,nargs)

! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.

character(len=*) :: str,delims
character(len=len_trim(str)) :: strsav
character(len=*),dimension(:) :: args
integer :: nargs, na, i, lenstr, k

strsav=str
call compact(str)
na=size(args)
do i=1,na
  args(i)=' '
end do  
nargs=0
lenstr=len_trim(str)
if(lenstr==0) return
k=0

do
   if(len_trim(str) == 0) exit
   nargs=nargs+1
   call split(str,delims,args(nargs))
   call removebksl(args(nargs))
end do   
str=strsav

end subroutine parse

!**********************************************************************

subroutine compact(str)

! Converts multiple spaces and tabs to single spaces; deletes control characters;
! removes initial spaces.

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str)):: outstr
integer :: k,isp,lenstr,i,ich

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
isp=0
k=0

do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)
  
  select case(ich)
  
    case(9,32)     ! space or tab character
      if(isp==0) then
        k=k+1
        outstr(k:k)=' '
      end if
      isp=1
      
    case(33:)      ! not a space, quote, or control character
      k=k+1
      outstr(k:k)=ch
      isp=0
      
  end select
  
end do

str=adjustl(outstr)

end subroutine compact

!**********************************************************************

subroutine removesp(str)

! Removes spaces, tabs, and control characters in string str

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str))::outstr
integer :: k,lenstr,i,ich

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0

do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)
  select case(ich)    
    case(0:32)  ! space, tab, or control character
         cycle       
    case(33:)  
      k=k+1
      outstr(k:k)=ch
  end select
end do

str=adjustl(outstr)

end subroutine removesp

!**********************************************************************

subroutine value_dr(str,rnum,ios)

! Converts number string to a double precision real number

character(len=*)::str
real(kr8)::rnum
integer :: ios,ilen,ipos

ilen=len_trim(str)
ipos=scan(str,'Ee')
if(.not.is_digit(str(ilen:ilen)) .and. ipos/=0) then
   ios=3
   return
end if
read(str,*,iostat=ios) rnum

end subroutine value_dr

!**********************************************************************

subroutine value_sr(str,rnum,ios)

! Converts number string to a single precision real number

character(len=*)::str
real(kr4) :: rnum
real(kr8) :: rnumd 
integer :: ios

call value_dr(str,rnumd,ios)
if( abs(rnumd) > huge(rnum) ) then
  ios=15
  return
end if
if( abs(rnumd) < tiny(rnum) ) rnum=0.0_kr4
rnum=rnumd

end subroutine value_sr

!**********************************************************************

subroutine value_di(str,inum,ios)

! Converts number string to a double precision integer value

character(len=*)::str
integer(ki8) :: inum
real(kr8) :: rnum
integer :: ios

call value_dr(str,rnum,ios)
if(abs(rnum)>huge(inum)) then
  ios=15
  return
end if
inum=nint(rnum,ki8)

end subroutine value_di

!**********************************************************************

subroutine value_si(str,inum,ios)

! Converts number string to a single precision integer value

character(len=*)::str
integer(ki4) :: inum
real(kr8) :: rnum
integer :: ios

call value_dr(str,rnum,ios)
if(abs(rnum)>huge(inum)) then
  ios=15
  return
end if
inum=nint(rnum,ki4)

end subroutine value_si

!**********************************************************************

subroutine shiftstr(str,n)
 
! Shifts characters in in the string 'str' n positions (positive values
! denote a right shift and negative values denote a left shift). Characters
! that are shifted off the end are lost. Positions opened up by the shift 
! are replaced by spaces.

character(len=*):: str
integer :: n,lenstr,nabs

lenstr=len(str)
nabs=iabs(n)
if(nabs>=lenstr) then
  str=repeat(' ',lenstr)
  return
end if
if(n<0) str=str(nabs+1:)//repeat(' ',nabs)  ! shift left
if(n>0) str=repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right 
return

end subroutine shiftstr

!**********************************************************************

subroutine insertstr(str,strins,loc)

! Inserts the string 'strins' into the string 'str' at position 'loc'. 
! Characters in 'str' starting at position 'loc' are shifted right to
! make room for the inserted string. Trailing spaces of 'strins' are 
! removed prior to insertion

character(len=*):: str,strins
character(len=len(str))::tempstr
integer :: loc,lenstrins

lenstrins=len_trim(strins)
tempstr=str(loc:)
call shiftstr(tempstr,lenstrins)
tempstr(1:lenstrins)=strins(1:lenstrins)
str(loc:)=tempstr
return

end subroutine insertstr

!**********************************************************************

subroutine delsubstr(str,substr)

! Deletes first occurrence of substring 'substr' from string 'str' and
! shifts characters left to fill hole. Trailing spaces or blanks are
! not considered part of 'substr'.

character(len=*):: str,substr
integer :: lensubstr,ipos

lensubstr=len_trim(substr)
ipos=index(str,substr)
if(ipos==0) return
if(ipos == 1) then
   str=str(lensubstr+1:)
else
   str=str(:ipos-1)//str(ipos+lensubstr:)
end if   
return

end subroutine delsubstr

!**********************************************************************

subroutine delall(str,substr)

! Deletes all occurrences of substring 'substr' from string 'str' and
! shifts characters left to fill holes.

character(len=*):: str,substr
integer :: lensubstr,ipos

lensubstr=len_trim(substr)
do
   ipos=index(str,substr)
   if(ipos == 0) exit
   if(ipos == 1) then
      str=str(lensubstr+1:)
   else
      str=str(:ipos-1)//str(ipos+lensubstr:)
   end if
end do   
return

end subroutine delall

!**********************************************************************

function uppercase(str) result(ucstr)

! convert string to upper case

character (len=*):: str
character (len=len_trim(str)):: ucstr
integer :: ilen,ioffset,iquote,i,iav,iqc

ilen=len_trim(str)
ioffset=iachar('A')-iachar('a')     
iquote=0
ucstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('a') .and. iav <= iachar('z')) then
    ucstr(i:i)=achar(iav+ioffset)
  else
    ucstr(i:i)=str(i:i)
  end if
end do
return

end function uppercase

!**********************************************************************

function lowercase(str) result(lcstr)

! convert string to lower case

character (len=*):: str
character (len=len_trim(str)):: lcstr
integer :: ilen,ioffset,iquote,i,iqc,iav

ilen=len_trim(str)
ioffset=iachar('A')-iachar('a')
iquote=0
lcstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('A') .and. iav <= iachar('Z')) then
    lcstr(i:i)=achar(iav-ioffset)
  else
    lcstr(i:i)=str(i:i)
  end if
end do
return

end function lowercase

!**********************************************************************

subroutine readline(nunitr,line,ios)

! Reads line from unit=nunitr, ignoring blank lines
! and deleting comments beginning with an exclamation point(!)

character (len=*):: line
integer :: nunitr,ios,ipos

do  
  read(nunitr,'(a)', iostat=ios) line      ! read input line
  if(ios /= 0) return
  line=adjustl(line)
  ipos=index(line,'!')
  if(ipos == 1) cycle
  if(ipos /= 0) line=line(:ipos-1)
  if(len_trim(line) /= 0) exit
end do
return

end subroutine readline

!**********************************************************************

subroutine match(str,ipos,imatch)

! Sets imatch to the position in string of the delimiter matching the delimiter
! in position ipos. Allowable delimiters are (), [], {}, <>.

character(len=*) :: str
character :: delim1,delim2,ch
integer :: ipos,imatch,lenstr,idelim2,istart,iend,inc,isum,i

lenstr=len_trim(str)
delim1=str(ipos:ipos)
select case(delim1)
   case('(')
      idelim2=iachar(delim1)+1
      istart=ipos+1
      iend=lenstr
      inc=1
   case(')')
      idelim2=iachar(delim1)-1
      istart=ipos-1
      iend=1
      inc=-1
   case('[','{','<')
      idelim2=iachar(delim1)+2
      istart=ipos+1
      iend=lenstr
      inc=1
   case(']','}','>')
      idelim2=iachar(delim1)-2
      istart=ipos-1
      iend=1
      inc=-1
   case default
      write(*,*) delim1,' is not a valid delimiter'
      return
end select
if(istart < 1 .or. istart > lenstr) then
   write(*,*) delim1,' has no matching delimiter'
   return
end if
delim2=achar(idelim2) ! matching delimiter

isum=1
do i=istart,iend,inc
   ch=str(i:i)
   if(ch /= delim1 .and. ch /= delim2) cycle
   if(ch == delim1) isum=isum+1
   if(ch == delim2) isum=isum-1
   if(isum == 0) exit
end do
if(isum /= 0) then
   write(*,*) delim1,' has no matching delimiter'
   return
end if   
imatch=i

return

end subroutine match

!**********************************************************************

subroutine write_dr(rnum,str,fmt)

! Writes double precision real number rnum to string str using format fmt

real(kr8) :: rnum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) rnum
str=adjustl(str)

end subroutine write_dr

!***********************************************************************

subroutine write_sr(rnum,str,fmt)

! Writes single precision real number rnum to string str using format fmt

real(kr4) :: rnum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) rnum
str=adjustl(str)

end subroutine write_sr

!***********************************************************************

subroutine write_di(inum,str,fmt)

! Writes double precision integer inum to string str using format fmt

integer(ki8) :: inum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) inum
str=adjustl(str)

end subroutine write_di

!***********************************************************************

subroutine write_si(inum,str,fmt)

! Writes single precision integer inum to string str using format fmt

integer(ki4) :: inum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) inum
str=adjustl(str)

end subroutine write_si

!***********************************************************************

subroutine trimzero(str)

! Deletes nonsignificant trailing zeroes from number string str. If number
! string ends in a decimal point, one trailing zero is added.

character(len=*) :: str
character :: ch
character(len=10) :: exp
integer :: ipos,lstr,i

ipos=scan(str,'eE')
if(ipos>0) then
   exp=str(ipos:)
   str=str(1:ipos-1)
endif
lstr=len_trim(str)
do i=lstr,1,-1
   ch=str(i:i)
   if(ch=='0') cycle          
   if(ch=='.') then
      str=str(1:i)//'0'
      if(ipos>0) str=trim(str)//trim(exp)
      exit
   endif
   str=str(1:i)
   exit
end do
if(ipos>0) str=trim(str)//trim(exp)

end subroutine trimzero

!**********************************************************************

subroutine writeq_dr(unit,namestr,value,fmt)

! Writes a string of the form <name> = value to unit

real(kr8) :: value
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr

call writenum(value,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_dr

!**********************************************************************

subroutine writeq_sr(unit,namestr,value,fmt)

! Writes a string of the form <name> = value to unit

real(kr4) :: value
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr

call writenum(value,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_sr

!**********************************************************************

subroutine writeq_di(unit,namestr,ivalue,fmt)

! Writes a string of the form <name> = ivalue to unit

integer(ki8) :: ivalue
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr
call writenum(ivalue,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_di

!**********************************************************************

subroutine writeq_si(unit,namestr,ivalue,fmt)

! Writes a string of the form <name> = ivalue to unit

integer(ki4) :: ivalue
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr
call writenum(ivalue,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_si

!**********************************************************************

function is_letter(ch) result(res)

! Returns .true. if ch is a letter and .false. otherwise

character :: ch
logical :: res

select case(ch)
case('A':'Z','a':'z')
  res=.true.
case default
  res=.false.
end select
return

end function is_letter

!**********************************************************************

function is_digit(ch) result(res)

! Returns .true. if ch is a digit (0,1,...,9) and .false. otherwise

character :: ch
logical :: res

select case(ch)
case('0':'9')
  res=.true.
case default
  res=.false.
end select
return

end function is_digit

!**********************************************************************

subroutine split(str,delims,before,sep)

! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the 
! found delimiter. A delimiter in 'str' is treated like an ordinary 
! character if it is preceded by a backslash (\). If the backslash 
! character is desired in 'str', then precede it with another backslash.

character(len=*) :: str,delims,before
character,optional :: sep
logical :: pres
character :: ch,cha
integer :: lenstr,k,ibsl,i,ipos,iposa

pres=present(sep)
str=adjustl(str)
call compact(str)
lenstr=len_trim(str)
if(lenstr == 0) return        ! string str is empty
k=0
ibsl=0                        ! backslash initially inactive
before=' '
do i=1,lenstr
   ch=str(i:i)
   if(ibsl == 1) then          ! backslash active
      k=k+1
      before(k:k)=ch
      ibsl=0
      cycle
   end if
   if(ch == '\') then          ! backslash with backslash inactive
      k=k+1
      before(k:k)=ch
      ibsl=1
      cycle
   end if
   ipos=index(delims,ch)         
   if(ipos == 0) then          ! character is not a delimiter
      k=k+1
      before(k:k)=ch
      cycle
   end if
   if(ch /= ' ') then          ! character is a delimiter that is not a space
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
   cha=str(i+1:i+1)            ! character is a space delimiter
   iposa=index(delims,cha)
   if(iposa > 0) then          ! next character is a delimiter
      str=str(i+2:)
      if(pres) sep=cha
      exit
   else
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
end do
if(i >= lenstr) str=''
str=adjustl(str)              ! remove initial spaces
return

end subroutine split

!**********************************************************************

subroutine removebksl(str)

! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str))::outstr
integer :: lenstr,k,ibsl,i

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0
ibsl=0                        ! backslash initially inactive

do i=1,lenstr
  ch=str(i:i)
  if(ibsl == 1) then          ! backslash active
   k=k+1
   outstr(k:k)=ch
   ibsl=0
   cycle
  end if
  if(ch == '\') then          ! backslash with backslash inactive
   ibsl=1
   cycle
  end if
  k=k+1
  outstr(k:k)=ch              ! non-backslash with backslash inactive
end do

str=adjustl(outstr)

end subroutine removebksl

!**********************************************************************
FUNCTION is_numeric(string)
  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: string
  LOGICAL :: is_numeric
  REAL :: x
  INTEGER :: e
  READ(string,*,IOSTAT=e) x
  is_numeric = e == 0
END FUNCTION is_numeric


! ------------------------
PURE FUNCTION Copy_a2s(a)  RESULT (s)    ! copy char array to string
CHARACTER,INTENT(IN) :: a(:)
CHARACTER(SIZE(a)) :: s
INTEGER :: i
DO i = 1,SIZE(a)
   s(i:i) = a(i)
END DO
END FUNCTION Copy_a2s

! ------------------------
PURE FUNCTION Copy_s2a(s)  RESULT (a)   ! copy s(1:Clen(s)) to char array
CHARACTER(*),INTENT(IN) :: s
CHARACTER :: a(LEN(s))
INTEGER :: i
DO i = 1,LEN(s)
   a(i) = s(i:i)
END DO
END FUNCTION Copy_s2a

! ------------------------
PURE INTEGER FUNCTION Clen(s)      ! returns same result as LEN unless:
CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
INTEGER :: i
Clen = LEN(s)
i = LEN_TRIM(s)
IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
END FUNCTION Clen

! ------------------------
PURE INTEGER FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
INTEGER :: i                       ! then len of C string is returned, note:
                                   ! Ctrim is only user of this function
i = LEN_TRIM(s) ; Clen_trim = i
IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
END FUNCTION Clen_trim

! ----------------
FUNCTION Ctrim(s1)  RESULT(s2)     ! returns same result as TRIM unless:
CHARACTER(*),INTENT(IN)  :: s1     ! last non-blank char is null in which
CHARACTER(Clen_trim(s1)) :: s2     ! case trailing blanks prior to null
s2 = s1                            ! are output
END FUNCTION Ctrim

! --------------------
INTEGER FUNCTION Count_Items(s1)  ! in string or C string that are blank or comma separated
CHARACTER(*) :: s1
CHARACTER(Clen(s1)) :: s
INTEGER :: i, k

s = s1                            ! remove possible last char null
k = 0  ; IF (s /= ' ') k = 1      ! string has at least 1 item
DO i = 1,LEN_TRIM(s)-1
   IF (s(i:i) /= ' '.AND.s(i:i) /= ',' &
                    .AND.s(i+1:i+1) == ' '.OR.s(i+1:i+1) == ',') k = k+1
END DO
Count_Items = k
END FUNCTION Count_Items

! --------------------
FUNCTION Reduce_Blanks(s)  RESULT (outs)
CHARACTER(*)      :: s
CHARACTER(LEN_TRIM(s)) :: outs
INTEGER           :: i, k, n

n = 0  ; k = LEN_TRIM(s)          ! k=index last non-blank (may be null)
DO i = 1,k-1                      ! dont process last char yet
   n = n+1 ; outs(n:n) = s(i:i)
   IF (s(i:i+1) == '  ') n = n-1  ! backup/discard consecutive output blank
END DO
n = n+1  ; outs(n:n)  = s(k:k)    ! last non-blank char output (may be null)
IF (n < k) outs(n+1:) = ' '       ! pad trailing blanks
END FUNCTION Reduce_Blanks

! ------------------
FUNCTION Replace_Text (s,text,rep)  RESULT(outs)
CHARACTER(*)            :: s,text,rep
CHARACTER(LEN(s)+100)   :: outs     ! provide outs with extra 100 char len
INTEGER                 :: i, nt, nr
INTEGER                 :: last_i  ! to avoid overwriting 
outs = s 
nt = LEN_TRIM(text) 
nr = LEN_TRIM(rep)

last_i = 0

DO
   i = INDEX(outs,text(:nt)) 
   if (i == last_i ) exit
   IF (i == 0) EXIT
   last_i = i
   outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
END DO
END FUNCTION Replace_Text

! ---------------------------------
FUNCTION Spack (s,ex)  RESULT (outs)
CHARACTER(*) :: s,ex
CHARACTER(LEN(s)) :: outs
CHARACTER :: aex(LEN(ex))   ! array of ex chars to extract
INTEGER   :: i, n

n = 0  ;  aex = Copy(ex)
DO i = 1,LEN(s)
   IF (.NOT.ANY(s(i:i) == aex)) CYCLE   ! dont pack char
   n = n+1 ; outs(n:n) = s(i:i)
END DO
outs(n+1:) = ' '     ! pad with trailing blanks
END FUNCTION Spack

! --------------------
INTEGER FUNCTION Tally (s,text)
CHARACTER(*) :: s, text
INTEGER :: i, nt

Tally = 0 ; nt = LEN_TRIM(text)
DO i = 1,LEN(s)-nt+1
   IF (s(i:i+nt-1) == text(:nt)) Tally = Tally+1
END DO
END FUNCTION Tally

! ---------------------------------
FUNCTION Translate(s1,codes)  RESULT (s2)
CHARACTER(*)       :: s1, codes(2)
CHARACTER(LEN(s1)) :: s2
CHARACTER          :: ch
INTEGER            :: i, j

DO i = 1,LEN(s1)
   ch = s1(i:i)
   j = INDEX(codes(1),ch) ; IF (j > 0) ch = codes(2)(j:j)
   s2(i:i) = ch
END DO
END FUNCTION Translate

! ---------------------------------
FUNCTION Upper(s1)  RESULT (s2)
CHARACTER(*)       :: s1
CHARACTER(LEN(s1)) :: s2
CHARACTER          :: ch
INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
INTEGER            :: i

DO i = 1,LEN(s1)
   ch = s1(i:i)
   IF (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
   s2(i:i) = ch
END DO
END FUNCTION Upper

! ---------------------------------
FUNCTION Lower(s1)  RESULT (s2)
CHARACTER(*)       :: s1
CHARACTER(LEN(s1)) :: s2
CHARACTER          :: ch
INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
INTEGER            :: i

DO i = 1,LEN(s1)
   ch = s1(i:i)
   IF (ch >= 'A'.AND.ch <= 'Z') ch = CHAR(ICHAR(ch)-DUC)
   s2(i:i) = ch
END DO
END FUNCTION Lower

!------------------------------------------------------------------------------------- 

!-------------------------------------------------------------------------------------                                                   
! count how many characters in the string 
!-------------------------------------------------------------------------------------                                                                 
function COUNTSUBSTRING (s1, s2) result(c)
  character(*), intent(in) :: s1, s2
  integer :: c, p, posn

  c = 0
  if(len(s2) == 0) return
  p = 1
  do
    posn = index(s1(p:), s2)
    if(posn == 0) return
    c = c + 1
    p = p + posn + len(s2)
  end do

end function

!-------------------------------------------------------------------------------------                                                   
! splits a string to substrings, returns array
!-------------------------------------------------------------------------------------

function SPLIT_STRING (str, separ, dims, word) result(error_status)
  character(*), intent(in) :: str
  character(*), intent(in) :: separ
  integer, intent(in) :: dims
  integer :: error_status
  character(100), dimension(:), allocatable, intent (out) :: word
  integer :: pos1, pos2, i

  error_status = 0
  pos1 = 1
  if (allocated (word) ) deallocate (word)
  allocate (word (dims))
  do i = 1, dims
    pos2 = INDEX(str(pos1:), separ)
    if (pos2 == 0) THEN
       word(i) = str(pos1:)
       EXIT
    endif
    word(i) = str(pos1:pos1+pos2-2)
    pos1 = pos2+pos1
 enddo

end function

!------------------------------------------------------------------------------------- 

function REPLACE_CHAR_IN_STRG (string_inout,target_char,substring_char, &
                    what_del) result(error_status)
 character(*), intent(inout) :: string_inout
 character(*), intent(in) :: target_char, substring_char, what_del
 integer :: indx, error_status

 error_status = 0
 indx = index(string_inout,target_char)
 if (indx > 0) then
    if (trim(what_del) == 'before') then
       string_inout = trim(substring_char)//trim(string_inout(indx+1:))
    elseif (trim(what_del) == 'after') then
       string_inout = trim(string_inout(:indx-1))//trim(substring_char)
    else
       error_status = 1
    endif
 else
    error_status = 1
    return
 endif

endfunction 



end module cx_string_tools_mod


