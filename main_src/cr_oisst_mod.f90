!$Id:$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: oisst_analysis.f90 (src)
!       OISST_ANALYSIS (program)
!
! PURPOSE: Routine to handle the Reynolds OISST analysis
!
! DESCRIPTION: see below
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
! REVISION HISTORY:
!  July 2004 - modified to look for previous file if nearest is not present
!  Apr 2007 - Lahey fortran does not all byte swapping at compilation so it requires
!   the CONVERT="BIG_ENDIAN" in the open statement.  This is not standard.
!  October 2009 - moved to daily 0.25 degree Reynolds SST
!
!
! data is one degree resolution
! 
! http://www.emc.ncep.noaa.gov/research/cmb/sst_analysis/
!
! global (89.875S - 89.875N)  (1440 x 720) starts at 89.5S and GM
!
!
!  DESCRIPTION OF THE DAILY OI SEA SURFACE TEMPERATURE (SST) ANALYSIS Version2
!  
!  The SST analysis is computed daily on a 0.25 degree latitude/longitude grid. 
!  This is version 1.0. There are two products with different satellite data. 
!  Both products use in situ data from ships and buoys. Also SSTs are generated 
!  for sea-ice concentrations above 50%. The sea ice for 1981-2004 is from 
!  http://nsidc.org/data/nsidc-0051.html 
!  (Cavalieri D., C. Parkinson, P. Gloerson, and H.J. Zwally. 1997, updated 2005. 
!  Sea ice concentrations from Nimbus-7 SMMR and DMSP SSM/I passive microwave data, 
!  June to September 2001. Boulder, CO, USA). The sea ice from 2005 to present is from 
!  http://polar.ncep.noaa.gov/seaice/ 
!  (Grumbine, R. W., 1996: Automated passive microwave sea ice concentration 
!  analysis at NCEP, 13pp. Unpublished manuscript available from NCEP/NWS/NOAA, 
!  5200 Auth Road, Camp Springs, MD, 20746, USA.)
!  
!  The first product uses NODC's AVHRR Pathfinder Version 5 
!  http://pathfinder.nodc.noaa.gov 
!  for September 1, 1981 though December 31, 2005 and the operational US Navy AVHRR data 
!  (May, D.A., M. M. Parmeter, D. S. Olszewski and B. D. McKenzie, 1998: Operational 
!  processing of satellite sea surface temperature retrievals at the Naval Oceanographic 
!  Office, Bull. Amer. Met. Soc., 79, 397-407) from January 1, 2006, through present. 
!  This product will henceforth be referred to as the AVHRR product.
!  
!  The second product adds AMSR-E version 5 data obtained from 
!  http://www.remss.com/ 
!  along with the AVHRR data used in version 1a and is available from June 1, 2002, 
!  (the start of AMSR-E) through present. This product will henceforth be 
!  termed the AVHRR + AMSR product.
!  
!  Both analyses include a bias correction of the satellite data with respect to 
!  in situ data using an empirical orthogonal teleconnection (EOT) algorithm. 
!  A short description of the complete analysis procedure can be found in the 
!  AMS extended abstract file (Reynolds-reviewed-rev.pdf).
!  
!  The SST analyses are available in individual daily files. The AVHRR product 
!  is named avhrr-only-v2.YYYYMMDD where YYYY is the year, MM is the month, 
!  and DD is the day. The files can be found on 
!  ftp://eclipse.ncdc.noaa.gov/pub/OI-daily-v2/IEEE/YYYY/AVHRR 
!  where YYYY is the year: 1981 to present. The files were written in IEEE 
!  binary (big-endian) and must be decompressed using gunzip. 
!  The AVHRR + AMSR-E product is written with the same format as the AVHRR product.
!  However, the file names are avhrr-only-v2.YYYYMMDD. The files can be found on
!  ftp://eclipse.ncdc.noaa.gov/pub/OI-daily-v2/IEEE/YYYY/AVHRR-AMSR 
!  where YYYY is the year: 2002 to present.
!  
!  Each file contains 4 records with integer*4 year, month, day, 
!  followed by a gridded integer*2 array. The first array is SST. 
!  The second array is the SST anomaly with respect to a 1971-2000 
!  base period. The third array is the sea ice concentration. 
!  The fourth array is the standard deviation of the analysis error 
!  which includes sampling, random and bias error. 
!  
!  Note: The SST, SST ANOMALY AND ERROR ARRAYS MUST BE MULTIPLIED BY 0.01 
!  TO CONVERT THE VALUES TO DEGREE C. The sea ice concentration array is in per cent (0-100). 
!  Missing values are -999.
!  
!  All arrays consist of 1440 spatial points in longitude from 0.125E to 359.875E 
!  in intervals of 0.25 increasing eastward, and 720 spatial points in latitude 
!  from 89.875S to 89.875N in intervals of 0.25 increasing northward.
!  
!  Each day consists of four FORTRAN records:
!  1. Three 4-byte integers for the year, month and day followed by 
!     1440*720 2-byte integer SST values.
!  2. Three 4-byte integers for the year, month and day followed by 
!     1440*720 2-byte integer SST anomaly values.
!  3. Three 4-byte integers for the year, month and day followed by 
!     1440*720 2-byte integer error values.
!  4. Three 4-byte integers for the year, month and day followed by 
!     1440*720 2-byte integer ice concentration values.
!  
!  Each record is written with a FORTRAN unformatted write which adds 
!  an extra 4 byte header and trailer word to the total record.
!--------------------------------------------------------------------------------------
module cr_oisst_mod



      
   use FILE_tools,only: &
    file_test &
    , getlun
   
   implicit none
   private
   
   
   integer, parameter :: int1 = selected_int_kind(1)
   integer, parameter :: int2 = selected_int_kind(3)
   integer, parameter :: int4 = selected_int_kind(8)
   integer, parameter :: int8 = selected_int_kind(10)
   integer, parameter :: real4 = selected_real_kind(6,37)
   integer, parameter :: real8 = selected_real_kind(15,307)
   integer, parameter :: ipre = real4
   real, parameter :: missing_value_real4 = -999.
   
   character(len=14) :: EXE_PROMPT = 'OISST READER> '
   
   !  public:: READ_OISST_ANALYSIS_MAP, GET_PIXEL_SST_ANALYSIS
   public:: GET_PIXEL_SST_ANALYSIS, get_oisst_map_filename, READ_OISST_ANALYSIS_MAP,get_oisst_data
   integer, parameter, public:: num_lon_sst_anal =1440, num_lat_sst_anal= 720 
   
   
   real, parameter, private:: min_sst_anal = -4.0, max_sst_anal=36.0
   integer(kind=int2), dimension(num_lon_sst_anal,num_lat_sst_anal):: temp_i2_buffer
   real(kind=real4), dimension(num_lon_sst_anal,num_lat_sst_anal), public, save:: oisst_anal_map
   
   real(kind=real4), dimension(num_lon_sst_anal,num_lat_sst_anal), public, save:: oisst_anal_map_uni
   real(kind=real4), dimension(num_lon_sst_anal,num_lat_sst_anal), public, save:: oisst_cice_map

   

contains

   !-------------------------------------------------------------------
   ! Function to find the oisst map name.
   !-------------------------------------------------------------------
   FUNCTION get_oisst_map_filename(date , oisst_path ) result(oisst_filename)
      
      use date_tools_mod, only: &
         date_type
      
      CHARACTER(*), intent(in) :: oisst_path
      type ( date_type ) , intent (in) :: date
      
      type ( date_type ) ::date_local
      CHARACTER(len=256) :: oisst_filename
      CHARACTER(len=256) :: oisst_filename_tmp
      CHARACTER(len=256) :: oisst_filename_full
      CHARACTER(len=256) :: oisst_filename_full_prel
      INTEGER(kind=int4) :: iday, year, month, day, jday, ileap
      CHARACTER (len=2)   :: day_string, month_string
      CHARACTER (len=4)   :: year_string
      
      INTEGER(kind=int4), parameter :: MAX_OISST_LATENCY = 5 !including current day

      oisst_filename = "no_file"
      
      date_local = date
      
      do iday = 0, MAX_OISST_LATENCY - 1
         call date_local % add_time ( day = -1 )

         oisst_filename_tmp = 'avhrr-only-v2.'//date_local % yyyymmdd
         oisst_filename_full = trim(oisst_path)//date_local % yyyy//'/'//trim(oisst_filename_tmp)//'.gz'
         oisst_filename_full_prel = trim(oisst_path)//date_local % yyyy//'/'//trim(oisst_filename_tmp)//'_preliminary.gz'
         
         !--- check for regular file
         if (file_test(trim(oisst_filename_full)) ) then
            oisst_filename = oisst_filename_full
            
            exit
         end if

         !--- check for preliminary file (true of recent files)
         if (file_test(trim(oisst_filename_full_prel)) ) then
            oisst_filename = oisst_filename_full_prel
            print *, EXE_PROMPT, "Found ", trim(oisst_filename)
            exit
         end if

      end do
      
      return

   END FUNCTION get_oisst_map_filename
   
   
   subroutine get_oisst_data ( oisst_filename, temp_path, lat , lon , sst ,  sea_ice , sst_uni)
      character ( len = * ) , intent(in) :: oisst_filename
      character ( len = * ) , intent(in) :: temp_path
      real , intent(in) :: lat ( :,:)
      real , intent(in) :: lon ( :,:)
      real , intent(out) :: sst ( :,:)
      real , intent(out) :: sea_ice ( :,:)
      real , intent(out) :: sst_uni ( :,:)
      character (len = 200 ) , save :: identifier_last_file
      
      
      sst = -999.
      sea_ice = -999. 
      if (trim(identifier_last_file) /= trim(oisst_filename)) then
         print *, EXE_PROMPT, "Found ", trim(oisst_filename)
         call READ_OISST_ANALYSIS_MAP(oisst_filename , temp_path)
         identifier_last_file =  oisst_filename 
      end if   
      call GET_PIXEL_SST_ANALYSIS( lat, lon , sst, sea_ice, sst_uni)
   
   end subroutine get_oisst_data

   !----------------------------------------------------------------------
   ! Read the Data
   !----------------------------------------------------------------------
   subroutine READ_OISST_ANALYSIS_MAP(sst_anal_name , temp_path)
      character(len=*), intent(in):: sst_anal_name
      character ( len = * ) ,intent(in) :: temp_path
      character(len=256):: file_temp
      character(len=256):: system_string
      integer(kind=int4):: iyrst,imst,idst,i,j
      integer:: oisst_lun
      integer:: ios0,ios1,ios2,ios3,ios4,ios5,erstat
      integer:: ii,jj

   
      !--- uncompress
      system_string = "gunzip -c "//trim(sst_anal_name)// &
        " > "//trim(temp_path)//"temp_oisst_file"
        
      
      call system(system_string)
  
      file_temp = trim(temp_path)//"temp_oisst_file"
      !------------------------------------------------------------------------------
      ! read in data
      !------------------------------------------------------------------------------
      ios0 = 0
      ios1 = 0
      ios2 = 0
      ios3 = 0
      ios4 = 0
      ios5 = 0
  
      oisst_lun = GETLUN()
      open(unit=oisst_lun,file=trim(file_temp),access="sequential",status="old", &
            action="read",form="unformatted", iostat = ios0,CONVERT="BIG_ENDIAN") !use this line for Lahey or gfortran

      if (ios0 /= 0) then
         erstat = 66
         print *,EXE_PROMPT, "Error opening oisst file, ",trim(sst_anal_name)," ios0 = ", ios0
         stop 66
      endif

      read(unit=oisst_lun,iostat=ios1) iyrst,imst,idst,((temp_i2_buffer(i,j),i=1,num_lon_sst_anal),j=1,num_lat_sst_anal)
      oisst_anal_map = temp_i2_buffer * 0.01

      read(unit=oisst_lun,iostat=ios2) iyrst,imst,idst,((temp_i2_buffer(i,j),i=1,num_lon_sst_anal),j=1,num_lat_sst_anal)

      read(unit=oisst_lun,iostat=ios3) iyrst,imst,idst,((temp_i2_buffer(i,j),i=1,num_lon_sst_anal),j=1,num_lat_sst_anal)
   

      read(unit=oisst_lun,iostat=ios4) iyrst,imst,idst,((temp_i2_buffer(i,j),i=1,num_lon_sst_anal),j=1,num_lat_sst_anal)
      oisst_cice_map = temp_i2_buffer * 0.01

      close(unit=oisst_lun,iostat=ios5)

      if (ios0+ios1+ios2+ios3+ios4+ios5 /= 0) then 

         print *, EXE_PROMPT,"Error reading OISST Analysis file, file = ", trim(sst_anal_name)

        

         return
         
      end if 

      print *, EXE_PROMPT,"OISST analysis read in successfully"


      !-------- convert from Celsius to Kelvin
      oisst_anal_map = oisst_anal_map + 273.15

      where(oisst_anal_map < 270.0)
         oisst_anal_map = missing_value_real4
      endwhere

      !---- compute the 3x3 uniformity of sst_anal
      do i = 1, num_lon_sst_anal
         ii = max(2,min(num_lon_sst_anal-1,i))
         do j = 1, num_lat_sst_anal
            jj = max(2,min(num_lat_sst_anal-1,j))
            oisst_anal_map_uni(i,j) = (maxval(oisst_anal_map(ii-1:ii+1,jj-1:jj+1)) - &
                        minval(oisst_anal_map(ii-1:ii+1,jj-1:jj+1)) ) / 3.0
         end do
      end do

      where (oisst_anal_map_uni > 8.0)
         oisst_anal_map_uni = 8.0
      end where

   end subroutine READ_OISST_ANALYSIS_MAP

   !----------------------------------------------------------------------------------------
   ! Interpolate SST Analysis to each AVHRR pixel
   !
   !  input: j1 - the first scan index of this segment
   !         j2 - the last scan index of this segment
   !----------------------------------------------------------------------------------------
   subroutine GET_PIXEL_SST_ANALYSIS( lat, lon , sst, sea_ice, sst_uni)
      real, intent(in) :: lat(:,:) , lon (:,:) 
      real, intent(out) ::sst(:,:) , sea_ice (:,:), sst_uni(:,:)
      integer :: j1,j2
      integer:: i
      integer:: j
      integer :: nx , ny
      integer:: ilon_sst_anal
      integer:: ilat_sst_anal
      integer:: ilon_sst_anal_x
      integer:: ilat_sst_anal_x
      real:: xlon
      real:: lon_sst_anal
      real:: lat_sst_anal
      real:: lon_weight
      real:: lat_weight
      
      real, parameter:: first_lon_sst_anal = 0.125, first_lat_sst_anal = -89.875, & 
                            last_lon_sst_anal = 359.875, last_lat_sst_anal = 89.875
      
      real, parameter:: del_lon_sst_anal = (last_lon_sst_anal - first_lon_sst_anal)/(num_lon_sst_anal-1)
      real, parameter:: del_lat_sst_anal = (last_lat_sst_anal - first_lat_sst_anal)/(num_lat_sst_anal-1)
      
      
       nx = size ( lat , 1)
       ny = size ( lat , 2)  
      do j = 1 , ny 
         do i = 1, nx

            if ( lon(i,j) < 0.0) then
               xlon =  lon(i,j) + 360.0
            else
               xlon =  lon(i,j)
            endif

            ilon_sst_anal = max(1,min(num_lon_sst_anal,int((xlon - first_lon_sst_anal)/del_lon_sst_anal+1)))
            ilat_sst_anal = max(1,min(num_lat_sst_anal,int(( lat(i,j) - first_lat_sst_anal)/del_lat_sst_anal+1)))
            ilon_sst_anal_x = ilon_sst_anal+1
            ilat_sst_anal_x = ilat_sst_anal+1

            ilon_sst_anal_x = max(1,min(num_lon_sst_anal,ilon_sst_anal_x))
            ilat_sst_anal_x = max(1,min(num_lat_sst_anal,ilat_sst_anal_x))

            lon_sst_anal = first_lon_sst_anal + (ilon_sst_anal-1) * del_lon_sst_anal
            lat_sst_anal = first_lat_sst_anal + (ilat_sst_anal-1) * del_lat_sst_anal

            lon_weight = max(0.0,min(1.0,(xlon - lon_sst_anal) / del_lon_sst_anal))
            lat_weight = max(0.0,min(1.0,( lat(i,j) - lat_sst_anal) / del_lat_sst_anal))

            if (minval(oisst_anal_map(ilon_sst_anal:ilon_sst_anal_x,ilat_sst_anal:ilat_sst_anal_x))  &
                     & /= missing_value_real4) then
               sst(i,j) = (1.0-lon_weight)*(1.0-lat_weight)*oisst_anal_map(ilon_sst_anal,ilat_sst_anal) + &
                     (lon_weight)*(1.0-lat_weight)*oisst_anal_map(ilon_sst_anal_x,ilat_sst_anal) + &
                     (1.0-lon_weight)*(lat_weight)*oisst_anal_map(ilon_sst_anal,ilat_sst_anal_x) + &
                     (lon_weight)*(lat_weight)*oisst_anal_map(ilon_sst_anal_x,ilat_sst_anal_x) 
            else
               
               sst(i,j) = oisst_anal_map(ilon_sst_anal,ilat_sst_anal)

            endif

    
            sea_ice(i,j) = oisst_cice_map(ilon_sst_anal,ilat_sst_anal)
            sst_uni(i,j) = oisst_anal_map_uni(ilon_sst_anal,ilat_sst_anal)

         end do
      end do

   end subroutine GET_PIXEL_SST_ANALYSIS


end module cr_oisst_mod
