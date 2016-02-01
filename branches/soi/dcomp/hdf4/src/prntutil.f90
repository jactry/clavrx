subroutine prntvals(data)

  include 'icahdfsdsinc.f90'

  type(hdf_data), intent(in)      :: &
       data

  call prnt__(data, 1, data%dimsize(1))

end subroutine prntvals

subroutine prntdimvals(data)

  include 'icahdfsdsinc.f90'

  type(hdf_data), intent(in)      :: &
       data

  integer, dimension(1:data%rank) :: &
       dind

  integer id, i, j, d

  if (data%rank == 1) then
     call prntvals(data)
     return
  endif

  dind = 1

250 continue

  i = dind(1)
  d = data%dimsize(1)
  do j = 2, data%rank
     i = i + (dind(j)-1)*d
     d = d*data%dimsize(j)
  enddo

  if (data%rank > 2) then
     print *
     print *, '(', (dind(j),',', j = data%rank, 3, -1), ' *, * )'
     print *, '--------------------------------'
  endif

  do j = 1, data%dimsize(2)
     call prnt__(data, i, i+data%dimsize(1)-1)
     i = i + data%dimsize(1)
  enddo
  
  do j = 3, data%rank
     if (dind(j).ge.data%dimsize(j)) then
        if ((j.eq.data%rank)) return
        dind(j) = 1
     else
        dind(j) = dind(j) + 1
        goto 250
     endif
  enddo
  
end subroutine prntdimvals

subroutine prnt__(data, i0, i1)
  
  include 'icahdfsdsinc.f90'
  
  type(hdf_data), intent(in) :: &
       data
  
  integer, intent(in)        :: &
       i0,                      &
       i1
  
  if (allocated(data%c1values)) print *, data%c1values(i0:i1)
  if (allocated(data%i1values)) print *, data%i1values(i0:i1)
  if (allocated(data%i2values)) print *, data%i2values(i0:i1)
  if (allocated(data%i4values)) print *, data%i4values(i0:i1)
! En prévision...
!  if (allocated(data%i8values)) print *, data%i8values(i0:i1)
  if (allocated(data%r4values)) print *, data%r4values(i0:i1)
  if (allocated(data%r8values)) print *, data%r8values(i0:i1)
  
end subroutine prnt__

