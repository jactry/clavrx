!
program TEST_CHARS

  implicit none

  character(len=10):: satname

  print *,"Enter sat name (up to 10 chars)"
  print *, "(Enter goes13 for this test)"
  read(unit=*,fmt="(a)") satname

  if (satname == "goes13") then
     print *, "it is a match"
  else
     print *, "it is not a match"
  endif

end program TEST_CHARS

