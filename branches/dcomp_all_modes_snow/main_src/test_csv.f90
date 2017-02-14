program test_csv

use strings

character ( len=200) ::csv_file_name

integer   ( kind = 4 ) csv_file_status
  integer   ( kind = 4 ) csv_file_unit
  integer   ( kind = 4 ) csv_record_status
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) line_num
  character ( len = 400 ) record
  integer   ( kind = 4 ) value_count
  character (len =300) :: before
   character ( len = 300) :: rec_arr ( 13 )
  integer :: j
csv_file_name='clavrx_prd.csv'

call csv_file_line_count ( csv_file_name, line_num )
print*,line_num

 write ( *, '(a,i8,a)' ) '  File contains ', line_num, ' lines.'

  call csv_file_open_read ( csv_file_name, csv_file_unit )


  do i = 1, line_num
    read ( csv_file_unit, '(a)', iostat = csv_file_status ) record
    write ( *, '(a)' ) i, trim ( record )
    call csv_value_count ( record, csv_record_status, value_count )
    write ( *, * ) i, value_count
    
    rec_arr = extract_single ( trim ( record ) )
    do j =1, 13
         print*,trim(rec_arr(j))
    end do     
    write (*,*) '========'
  end do

  call csv_file_close_read ( csv_file_name, csv_file_unit )
  
  
contains

   function extract_single ( record) result  (record_single)
      character ( len = * ) :: record
      character ( len = 300) :: record_single ( 13 )
      character (len = 1000 ) :: record_local
      character ( len =300) :: before
      integer :: i
      
      record_local = trim (record)
      n_val = size ( record_single)
      record_single(:) = 'not_set' 
      do i = 1, n_val 
         call split ( record_local , ',', before)
         record_single ( i ) = trim(before)
      
      end do
      
      
      
   
   
   end function extract_single
  
end program test_csv
