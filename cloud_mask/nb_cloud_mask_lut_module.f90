!$Id$
!----------------------------------------------------------------------
! MODULE name: NB_CLOUD_MASK
! 
! Routines for the reading of LUTS of the naive Bayesian cloud mask
!
! Authors: Andrew Heidinger, NOAA/NESDIS
!          Andi Walther, CIMSS
!          Denis Botambekov, CIMSS
!          William Straka, CIMSS
!
! DEPENDENCIES: NETCDF library
!
!----------------------------------------------------------------------
module NB_CLOUD_MASK_LUT_MODULE

 use NB_CLOUD_MASK_NETCDF_READ
 use NB_CLOUD_MASK_SERVICES, only: int1, int2, int4, real4, real8, MISSING_VALUE_REAL4

 implicit none

 public:: READ_NAIVE_BAYES_LUT
 public:: RESET_NB_CLOUD_MASK_LUT
 public:: READ_PRIOR
 public:: COMPUTE_PRIOR
 public:: RESET_NB_CLOUD_MASK_PRIOR_LUT
 public:: SET_CLOUD_MASK_LUT_VERSION

 character(len=*), parameter, private :: EXE_PROMPT_CM = "Naive Bayesian Cloud Mask LUT Module >> "
 character(len=100), private:: Cloud_Mask_Lut_Version

 !--- set thresholds and algorithm specific constants
 include 'nb_cloud_mask.inc'

 !--------------------------------------------------------------------------------
 ! store table values as module-wide arrays
 !--------------------------------------------------------------------------------
 integer, public, save:: N_class
 integer, public, save:: N_Sfc_Bayes 
 integer, public, save:: N_bounds
 real, dimension(:), allocatable, public, save:: Prior_Yes
 real, dimension(:), allocatable, public, save:: Prior_No
 real, dimension(:), allocatable, public, save:: Optimal_Posterior_Prob
 integer, dimension(:), allocatable, public, save:: Skip_Sfc_Type_Flag
 real, dimension(:,:), allocatable, public, save:: Classifier_Bounds_Min
 real, dimension(:,:), allocatable, public, save:: Classifier_Bounds_Max
 real, dimension(:,:), allocatable, public, save:: Delta_Classifier_Bounds
 integer, dimension(:,:), allocatable, public, save:: First_valid_Classifier_Bounds
 integer, dimension(:,:), allocatable, public, save:: Last_valid_Classifier_Bounds
 real, dimension(:,:,:), allocatable, public, save:: Class_Cond_Ratio
 character (len=30), dimension(:,:), allocatable, public, save:: Classifier_Value_Name
 integer, dimension(:), allocatable, public, save:: Class_To_Test_Idx
 integer, public, save:: Class_Idx_Etropo

 logical, public, save:: Is_Classifiers_Read = .false.

 ! netCDF parameters
 integer, parameter, private :: sds_rank_1d = 1
 integer, dimension(sds_rank_1d), private :: sds_start_1d, sds_edge_1d, &
          sds_stride_1d

 integer, parameter, private :: sds_rank_2d = 2
 integer, dimension(sds_rank_2d), private :: sds_start_2d, sds_edge_2d, &
          sds_stride_2d, sds_dims_2d

 integer, parameter, private :: sds_rank_3d = 3
 integer, dimension(sds_rank_3d), private :: sds_start_3d, sds_edge_3d, &
          sds_stride_3d, sds_dims_3d

 !--- Prior Variables
 logical, public, save:: Is_Prior_Read
 integer(kind=int4), private, save:: nlon_prior, nlat_prior, nmonths_prior, ndiurnal_prior
 real(kind=real4), private, save:: dlon_prior, dlat_prior
 real(kind=real4), private, save:: lon_min_prior, lon_max_prior
 real(kind=real4), private, save:: lat_min_prior, lat_max_prior
 real(kind=real4), dimension(:,:,:,:), allocatable, private, save:: Prior_Table

 !-------------------------------------------------------------------
 integer(kind=int1), private, parameter:: symbol_no = 0
 integer(kind=int1), private, parameter:: symbol_yes = 1

  contains

!-------------------------------------------------------------------------------
! NetCDF LUT routines:
!-------------------------------------------------------------------------------
subroutine SET_CLOUD_MASK_LUT_VERSION(Cloud_Mask_Lut_Version_Out)
  character(len=*), intent(out):: Cloud_Mask_Lut_Version_Out
  integer:: n_in, n_out
  n_in = len(Cloud_Mask_Lut_Version)
  n_out = len(Cloud_Mask_Lut_Version_Out)
  Cloud_Mask_Lut_Version_Out = Cloud_Mask_Lut_Version(1:min(n_in,n_out))
end subroutine SET_CLOUD_MASK_LUT_VERSION

!====================================================================
! SUBROUTINE Name: READ_NAIVE_BAYES_LUT
!
! Function:
!   Allocate and Read in the LUT needed for Bayesian cloud mask tables
!
!====================================================================
 subroutine READ_NAIVE_BAYES_LUT(Naive_Bayes_File_Name_Full_Path, &
                                 Cloud_Mask_Bayesian_Flag)

   character(len=*), intent(in):: Naive_Bayes_File_Name_Full_Path

   !--- Need a method to flag things
   integer, intent(out):: Cloud_Mask_Bayesian_Flag

   !local variables
   integer:: ncid
   integer:: status
   character(30):: var_name
   integer:: Sfc_Idx
   integer:: Class_Idx

   Is_Classifiers_Read = .FALSE.

   status = nf90_open(Naive_Bayes_File_Name_Full_Path, mode = nf90_nowrite, ncid = ncid)
   if (status /= nf90_noerr) then
      print *, EXE_PROMPT_CM , 'ERROR: Bayesian Cloud Mask Classifier Open Failed '
      print *, EXE_PROMPT_CM , 'Bayesian Cloud Mask Turned Off'
      cloud_mask_bayesian_flag = symbol_NO
      print*, EXE_PROMPT_CM ,' filename is: ', Naive_Bayes_File_Name_Full_Path
      return
   endif

   status = nf90_get_att(ncid, nf90_global, "data_file", Cloud_Mask_Lut_Version)
   if (status /= nf90_noerr) then
      print *, EXE_PROMPT_CM , 'ERROR: Bayesian Cloud Mask Version Read Failed'
      return
   endif

   status = nf90_get_att(ncid, nf90_global, "n_class", N_Class)
   status = nf90_get_att(ncid, nf90_global, "n_bounds_reg", N_bounds)
   status = nf90_get_att(ncid, nf90_global, "n_sfc_type", N_sfc_bayes)

   !--- allocate
   allocate(Prior_Yes(N_sfc_bayes))
   allocate(Prior_No(N_sfc_bayes))
   allocate(Optimal_Posterior_Prob(N_sfc_bayes))
   allocate(Skip_Sfc_Type_Flag(N_sfc_bayes))
   allocate(Classifier_Bounds_Min(N_class,N_sfc_bayes))
   allocate(Classifier_Bounds_Max(N_class,N_sfc_bayes))
   allocate(Delta_Classifier_Bounds(N_class,N_sfc_bayes))
   allocate(First_valid_Classifier_Bounds(N_class,N_sfc_bayes))
   allocate(Last_valid_Classifier_Bounds(N_class,N_sfc_bayes))
   allocate(Class_Cond_Ratio(N_bounds-1,N_class,N_sfc_bayes))
   allocate(Classifier_Value_Name(N_class,N_sfc_bayes))
   allocate(Class_To_Test_Idx(N_class))

   !--- initialize
   Prior_Yes = Missing_Value_Real4
   Prior_No = Missing_Value_Real4
   Optimal_Posterior_Prob = Missing_Value_Real4
   First_valid_Classifier_Bounds = 0
   Last_valid_Classifier_Bounds = 0
   Skip_Sfc_Type_Flag = symbol_NO


   !Now read in to the 1D variables

   var_name="prior_yes"
   call read_netcdf_1d_real(ncid,N_sfc_bayes,var_name,Prior_Yes)

   var_name="prior_no"
   call read_netcdf_1d_real(ncid,N_sfc_bayes,var_name,Prior_No)

   var_name="optimal_posterior_prob"
   call read_netcdf_1d_real(ncid, N_sfc_bayes, var_name, Optimal_Posterior_Prob)

   !--- if data are missing, skip this surface type
   do Sfc_Idx = 1, N_sfc_bayes
      if (Optimal_Posterior_Prob(Sfc_Idx) == Missing_Value_Real4) then
            Skip_Sfc_Type_Flag(Sfc_Idx) = symbol_YES
      endif
   end do

   !Now the 2D variables

   sds_start_2d = 1
   sds_edge_2d(1) = N_class
   sds_edge_2d(2) = N_sfc_bayes

   var_name="bin_start"
   call read_netcdf_2d_real(ncid, sds_start_2d, sds_edge_2d, &
                              var_name,Classifier_Bounds_Min)

   var_name="bin_end" !real
   call read_netcdf_2d_real(ncid, sds_start_2d, sds_edge_2d, & 
                              var_name,Classifier_Bounds_Max)

   var_name="delta_bin" !real
   call read_netcdf_2d_real(ncid, sds_start_2d, sds_edge_2d, &
                              var_name,Delta_Classifier_Bounds)

   var_name="first_valid_bounds" !integer
   call read_netcdf_2d_int(ncid, sds_start_2d, sds_edge_2d, &
                              var_name, First_valid_Classifier_Bounds)

   var_name="last_valid_bounds" !integer
   call read_netcdf_2d_int(ncid, sds_start_2d, sds_edge_2d, &
                              var_name, Last_valid_Classifier_Bounds)

   var_name="classifier_names" !character
   call read_netcdf_2d_char(ncid, sds_start_2d, sds_edge_2d, var_name, Classifier_Value_Name)

   !finally 3D variables
   sds_start_3d = 1
   sds_edge_3d(1) = N_bounds-1
   sds_edge_3d(2) = N_class
   sds_edge_3d(3) = N_sfc_bayes

   var_name="class_cond_ratio_reg" !real
   call read_netcdf_3d_real(ncid, sds_start_3d, sds_edge_3d, &
                            var_name,Class_Cond_Ratio)

   status = nf90_close(ncid)

   Is_Classifiers_Read = .true.

   !---set up Classifier to Test Mapping
   do Class_Idx = 1, N_Class

          select case (trim(Classifier_Value_Name(Class_Idx,1)))
                    case("T_11")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 1
                    case("T_Max-T")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 2
                    case("T_Std")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 3
                    case("Emiss_Tropo")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 4
                       Class_Idx_Etropo = Class_Idx
                    case("FMFT")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 5
                    case("Btd_11_67")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 6
                    case("Bt_11_67_Covar")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 7
                    case("Btd_11_85")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 8
                    case("Emiss_375")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 9
                    case("Btd_375_11_All")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 9
                    case("Btd_375_11_Day")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 10
                    case("Emiss_375_Day")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 10
                    case("Btd_375_11_Night")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 11
                    case("Emiss_375_Night")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 11
                    case("Spare")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 12
                    case("Ref_063_Day")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 13
                    case("Ref_Std")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 14
                    case("Ref_063_Min_3x3_Day")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 15
                    case("Ref_Ratio_Day")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 16
                    case("Ref_138_Day")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 17
                    case("Ndsi_Day")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 18
                    case default
                       print *, "Unknown Classifier Naive Bayesian Cloud Mask, returning"
                       print *, "Name = ",trim(Classifier_Value_Name(Class_Idx,1))
                       return
           end select

         end do


 end subroutine READ_NAIVE_BAYES_LUT


!-----------------------------------------------------------------------------------------
! This routine deallocates all LUT arrays and resets Is_Classifiers_Read to be  false
!
! This should be called from the bridge code using whatever mechanism exists to 
! identify the appropriate time to reset the LUT.  This could be
! 1. end of a file
! 2. change from one sensor to another
!-----------------------------------------------------------------------------------------
subroutine RESET_NB_CLOUD_MASK_LUT()
        if (allocated(Prior_Yes)) deallocate(Prior_Yes)
        if (allocated(Prior_No)) deallocate(Prior_No)
        if (allocated(Optimal_Posterior_Prob)) deallocate(Optimal_Posterior_Prob)
        if (allocated(Skip_Sfc_Type_Flag)) deallocate(Skip_Sfc_Type_Flag)
        if (allocated(Classifier_Bounds_Min)) deallocate(Classifier_Bounds_Min)
        if (allocated(Classifier_Bounds_Max)) deallocate(Classifier_Bounds_Max)
        if (allocated(Delta_Classifier_Bounds)) deallocate(Delta_Classifier_Bounds)
        if (allocated(First_Valid_Classifier_Bounds)) deallocate(First_Valid_Classifier_Bounds)
        if (allocated(Last_Valid_Classifier_Bounds)) deallocate(Last_Valid_Classifier_Bounds)
        if (allocated(Class_Cond_Ratio)) deallocate(Class_Cond_Ratio)
        if (allocated(Classifier_Value_Name)) deallocate(Classifier_Value_Name)
        if (allocated(Class_To_Test_Idx)) deallocate(Class_To_Test_Idx)
        Is_Classifiers_Read = .false.
end subroutine RESET_NB_CLOUD_MASK_LUT 


!-----------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------
subroutine RESET_NB_CLOUD_MASK_PRIOR_LUT()

 if (allocated(Prior_Table)) deallocate(Prior_Table)

 Is_Prior_Read = .false.

end subroutine RESET_NB_CLOUD_MASK_PRIOR_LUT
!================D============================================================
!dimensions:
!       longitude = 36 ;
!       latitude = 18 ;
!       month = 12 ;
!       diurnal = 3 ;
!Variables:
!       float cloud_fraction_table(diurnal, month, latitude, longitude) ;
!       float cloud_fraction_table_smoothed(diurnal, month, latitude, longitude)
!
!
! global attributes:
!               :description = "Andy Heidinger made this" ;
!               :number_longitudes = 36.f ;
!               :number_latitudes = 18.f ;
!               :number_months = 12.f ;
!               :number_times = 3.f ;
!               :latitude_min = -90.f ;
!               :latitude_max = 90.f ;
!               :longitude_min = -180.f ;
!               :longitude_max = 180.f ;
!               :longitude_spacing = 10.f ;
!               :latitude_spacing = 10.f ;
!============================================================================
subroutine READ_PRIOR(File_Name)

  character(len=*), intent(in):: File_Name

  integer:: status
  integer:: ncid

  Is_Prior_Read = .FALSE.

  status = nf90_open(File_Name, mode = nf90_nowrite, ncid = ncid)

  if (status /= nf90_noerr) then
      print *, EXE_PROMPT_CM , 'ERROR: NB Prior File Read Failed'
      print *, EXE_PROMPT_CM , 'Default Priors Used'
      print*, EXE_PROMPT_CM ,' filename is: ', File_Name
      return
  endif

  status = nf90_get_att(ncid, nf90_global, "number_longitudes", nlon_prior)
  status = nf90_get_att(ncid, nf90_global, "number_latitudes", nlat_prior)
  status = nf90_get_att(ncid, nf90_global, "number_months", nmonths_prior)
  status = nf90_get_att(ncid, nf90_global, "number_times", ndiurnal_prior)
  status = nf90_get_att(ncid, nf90_global, "longitude_spacing", dlon_prior)
  status = nf90_get_att(ncid, nf90_global, "latitude_spacing", dlat_prior)
  status = nf90_get_att(ncid, nf90_global, "longitude_min", lon_min_prior)
  status = nf90_get_att(ncid, nf90_global, "longitude_max", lon_max_prior)
  status = nf90_get_att(ncid, nf90_global, "latitude_min", lat_min_prior)
  status = nf90_get_att(ncid, nf90_global, "latitude_max", lat_max_prior)

  allocate(Prior_Table(nlon_prior,nlat_prior,nmonths_prior,ndiurnal_prior))

  Prior_Table = MISSING_VALUE_REAL4

  call read_netcdf_4d_real(ncid, (/1,1,1,1/), (/nlon_prior,nlat_prior,nmonths_prior, ndiurnal_prior/), &
                           "cloud_fraction_table_smoothed",Prior_Table)

  status = nf90_close(ncid)

  Is_Prior_Read = .true.

end subroutine READ_PRIOR
!============================================================================
subroutine COMPUTE_PRIOR(Lon,Lat,Month,Prior)
    real, dimension(:,:), intent(in):: Lon
    real, dimension(:,:), intent(in):: Lat
    integer(kind=int2), intent(in):: Month
    real, dimension(:,:), intent(out):: Prior

    integer:: Nx, Ny
    integer, parameter:: idiurnal = 1    !1 = daily averaged
    integer:: i,j,ilon,ilat,imonth

    Prior = MISSING_VALUE_REAL4

    Nx = size(Lon,1)
    Ny = size(Lon,2)

    imonth = Month

    do i = 1, Nx
       do j = 1, Ny

        ilon = min(Nlon_Prior,max(1,int((Lon(i,j) - Lon_Min_Prior) / (dLon_Prior))+1))
        ilat = min(Nlat_Prior,max(1,int((Lat(i,j) - Lat_Min_Prior) / (dLat_Prior))+1))

        Prior(i,j) = Prior_Table(ilon,ilat,imonth,idiurnal)

       enddo
    enddo

end subroutine COMPUTE_PRIOR

!-----------------------------------------------------------------------------------

end module NB_CLOUD_MASK_LUT_MODULE

