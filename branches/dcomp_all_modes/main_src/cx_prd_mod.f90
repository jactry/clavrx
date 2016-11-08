! $Id$
module cx_prd_mod
   
      use csv_mod, only: &
      csv_file_close_read &
      , csv_file_line_count &
      , csv_file_open_read &
      , csv_value_count
   
    use cx_string_tools_mod, only: &
      split
   
    use CONSTANTS, only: &
      int4   
   
   type prd_individual_dtype
      logical :: switch
      integer :: dim
      character(len=50) :: name
      character(len=50) :: name_clavrx
      integer :: dtype
      integer :: scaling
      real :: act_min
      real :: act_max
      integer :: val_min
      integer :: val_max
      integer :: val_fill
      real :: add_offset
      real :: slope
      character(len=300) :: standard_name
      character(len=30) :: unit
      character (len=300):: long_name
      character (len=300):: flagvalues
      
      integer(kind=int4) :: Sds_Id
      
   end type prd_individual_dtype
   
   type prd_dtype
      logical :: is_set = .false.
      character(len=200) :: csv_filename ='clavrx_level2_products.csv'
      integer :: num_products
      integer :: num_products_on
      type (prd_individual_dtype), allocatable :: product (:) 
      contains
      procedure :: read_products
   end type prd_dtype
   
   public:: prd_dtype
   
   contains
   !
   !   read csv file into product structure
   !
   subroutine read_products ( this)
      implicit none
      class (prd_dtype) :: this
      
      integer :: line_num
      integer  (kind = int4) :: csv_fle_unit
      integer  (kind = int4) :: csv_file_status
      integer  (kind = int4) :: csv_record_status
      integer  (kind = int4) :: csv_file_unit
      character ( len = 1000 ) record
      integer :: i_prd
      integer   ( kind = 4 ) :: value_count
      integer, parameter :: N_CSV = 12
      character ( len = 300) :: rec_arr ( N_CSV )
      
      call csv_file_line_count ( this % csv_filename, line_num )
      this % num_products = line_num - 1
      allocate (this % product (this % num_products))
      call csv_file_open_read ( this % csv_filename, csv_file_unit )
      
      ! - read header which is on first line
      read ( csv_file_unit, '(a)', iostat = csv_file_status ) record
      
      do i_prd = 1, this % num_products
        
         read ( csv_file_unit, '(a)', iostat = csv_file_status ) record
         call csv_value_count ( record, csv_record_status, value_count )
         rec_arr = extract_single ( trim ( record ) )
         
         this % product(i_prd) % switch = trim(rec_arr(1))  .eq. "1"
              
         read ( rec_arr(2), * ) this%product(i_prd)%name 
         read( rec_arr(3), '(a)' ) this%product(i_prd)%name_clavrx
         read ( rec_arr(4), * ) this%product(i_prd)%dim
         read ( rec_arr(5), * ) this%product(i_prd)%dtype
         read ( rec_arr(6), * ) this%product(i_prd)%scaling
         read ( rec_arr(7), * ) this%product(i_prd)%act_min
         read ( rec_arr(8), * ) this%product(i_prd)%act_max
         !read ( rec_arr(9), '(a)' ) this%product(i_prd)%add_offset
         !read ( rec_arr(10), '(a)' ) this%product(i_prd)%slope
         read ( rec_arr(9), '(a)' ) this%product(i_prd)%standard_name
         read ( rec_arr(10), '(a)' ) this%product(i_prd)%unit
         read ( rec_arr(11), '(a)' ) this%product(i_prd)%long_name
         read ( rec_arr(12), '(a)' ) this%product(i_prd)%flagvalues
         
      end do
       
      call csv_file_close_read ( this % csv_filename, csv_file_unit )
      this % is_set = .true.
      
      contains
      
      !
      !
      !
      function extract_single ( record) result  (record_single)
         implicit none
         character ( len = * ) :: record
         character ( len = 300) :: record_single ( N_CSV )
         character (len = 1000 ) :: record_local
         character ( len =300) :: before
         integer :: i
         integer :: n_val
      
         record_local = trim (record)
         n_val = size ( record_single)
         record_single(:) = 'not_set' 
         do i = 1, n_val 
            call split ( record_local , ',', before)
            record_single ( i ) = trim(before)
         end do
      
      end function extract_single
      
      
   end subroutine 



end module cx_prd_mod
