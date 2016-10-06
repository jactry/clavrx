module cx_grid_tools_mod
   real, parameter :: real4 = selected_real_kind(6,37)
   integer, parameter :: int4 = selected_int_kind(8)
contains

   !====================================================================
   ! SUBROUTINE Name: INDEX_IN_REGULAR_GRID
   !
   ! Function: 
   !       routines to find indicies in regular grid
   !
   ! Description:
   !   This subroutine is a routine to remap data from lat/lon to a 
   !   static grid.
   !
   ! Inputs:
   !       Lon_Output = vector of the longitudes of the grid
   !       Lat_Output = vector of the latitudes of the grid
   !       Ielem_Output = array of element indices for each grid point
   !       Iline_Output = array of line indices for each grid point
   !       Random_Flag = set true to select a random point in the grid, set
   !                     set false to select nearest point to grid lattice
   !
   ! Output:
   !
   !  2016/09/30: Rewrite (AW)
   !====================================================================
   subroutine INDEX_IN_REGULAR_GRID( &
        Lat_South, &
        dlat,      &
        Nlat,      &
        Lon_West,  &
        dlon,      &
        Nlon,      &
        Random_Flag &
        , lon_inp &
        , lat_inp &
        , gap_inp &
        , Ielem_out &
        , Iline_out )
      
      real(kind=real4), intent(in):: Lat_South
      real(kind=real4), intent(in):: dlat
      integer(kind=int4), intent(in):: Nlat
      real(kind=real4), intent(in):: Lon_West
      real(kind=real4), intent(in):: dlon
      integer(kind=int4), intent(in):: Nlon
      logical, intent(in):: Random_Flag
      real, intent(in) :: lon_inp (:,:)
      real, intent(in) :: lat_inp (:,:)
      logical, intent(in) :: gap_inp(:,:)
      integer, intent(inout), allocatable :: Ielem_out(:,:)
      integer, intent(inout), allocatable :: Iline_out(:,:)

      real(kind=real4):: big_number
      real(kind=real4):: Max_Allowable_Distance
      real(kind=real4), dimension(:,:),allocatable:: Min_Dist_Grid
      integer:: Ielem
      integer:: Iline
      integer:: Ilon
      integer:: Ilat
      integer:: nx
      integer:: ny
      real:: Lat_North
      real:: Lon_East
      real:: xdist
      real:: xdist_01
      real:: xdist_02
      real:: xdist_10
      real:: xdist_20
      real:: lat_width
      real:: lon_width
      integer:: Ilon_1, Ilon_2
      integer:: Ilon_min, Ilon_max
      real:: lon_1, lon_2
      integer:: Ilat_1, Ilat_2
      integer:: Ilat_min, Ilat_max
      real:: lat_1, lat_2
      integer:: i1,i2,j1,j2
      logical :: dateline_flag
      real:: Lon_East_abs
      real:: lon_inp_abs
      real:: xrand

      big_number = 999999.9

      !--- construct boundaries
      Lon_East_abs = Lon_West + dlon * (Nlon-1)
      Lon_East = Lon_East_abs
      if (Lon_East > 180.0) Lon_East = Lon_East - 360.0

      dateline_flag = .false.
      if (Lon_West > 0 .and. Lon_East < 0.0) dateline_flag = .true.

      Lat_North = Lat_South + dlat * (Nlat-1)

      nx = size(Lat_Inp,1)
      ny = size(Lat_Inp,2)
      
      print*,nx,ny,nlon,nlat
      print*,lat_south,dlat,nlat,lat_north
      

      Ielem_Output = Missing_Value_Int4
      Iline_Output = Missing_Value_Int4

      Max_Allowable_Distance = 1.0*sqrt(dlat**2 + dlon**2)

      ALLOCATE(Min_Dist_Grid(Nlon,Nlat))
      
      Min_Dist_Grid = big_number

      do Ielem = 1, nx
         do Iline = 1, ny
            
            if ( lat_inp(Ielem,Iline) .LT. lat_south .OR. lat_inp(Ielem,Iline) .GT. lat_north ) cycle
            
            if ( .NOT. dateline_flag .AND. ( Lon_Inp(Ielem,Iline) .LT. Lon_West) ) cycle
            if ( .NOT. dateline_flag .AND. ( Lon_Inp(Ielem,Iline) .GT. Lon_East) ) cycle
                       
            lon_inp_abs = lon_inp(Ielem,Iline)
            if (lon_inp_abs < 0.0) lon_inp_abs = lon_inp_abs + 360.0
            
            if ( dateline_flag .AND.  lon_inp_abs .LT. Lon_West  ) cycle
            if ( dateline_flag .AND.  lon_inp_abs .GT. Lon_East_abs  ) cycle
            
       
            i1 = max(1,Ielem-1)
            i2 = min(nx,Ielem+1)
            j1 = max(1,Iline-1)
            j2 = min(ny,Iline+1)

            lon_width = abs(lon_inp(i1,Iline) - lon_inp(i2,Iline))
            lon_1 = lon_inp(Ielem,Iline) + lon_width/2.0
            lon_2 = lon_inp(Ielem,Iline) - lon_width/2.0

            !-- lon width constraints prevents large values in presence of scanline jumps
            lat_width = min(lon_width,abs(Lat_Inp(Ielem,j1) - Lat_Inp(Ielem,j2)))
            lat_1 = lat_inp(Ielem,Iline) + lat_width/2.0
            lat_2 = lat_inp(Ielem,Iline) - lat_width/2.0

            call POSITION_TO_INDICES (Ilat, Ilon, Lat_Inp(Ielem,Iline), Lon_Inp(Ielem,Iline),  &
               Lat_South, Lat_North, dlat, Nlat, Lon_West, Lon_East,dlon, Nlon, xdist)
               
            if ( Ilon .LT.1 ) cycle
            if ( Ilat .LT.1 ) cycle   
            if ( xdist .GT. Max_Allowable_Distance) cycle   
            
            !----- if random, replace true xdist with a random number
            if (Random_Flag ) then
               call random_number(xrand)
               xdist = xrand * Max_Allowable_Distance
            endif
            
            !update closest grid point
            if ( .NOT. Gap_Inp(ielem,iline) ) then
               
               
               Ielem_Out(Ilon,Ilat)  = Ielem
              
               Iline_Out(Ilon,Ilat)  = Iline
               
               Min_Dist_Grid(Ilon,Ilat) = xdist
            endif
            

            if (xdist .GE. Min_Dist_Grid(Ilon,Ilat)) cycle  !check for pixels that are closer 
            
            
            !update surrounding grid points that fall within the pixel's footprint but
            !have not already been filled in by a closest pixel
            call POSITION_TO_INDICES (Ilat, Ilon_1, lat_inp(Ielem,Iline), lon_1,  &
               Lat_South, Lat_North, dlat, Nlat, Lon_West, Lon_East,dlon, Nlon, xdist_10)

            call POSITION_TO_INDICES (Ilat, Ilon_2, lat_inp(Ielem,Iline), lon_2,  &
               Lat_South, Lat_North, dlat, Nlat, Lon_West, Lon_East,dlon, Nlon, xdist_20)

            call POSITION_TO_INDICES (Ilat_1, Ilon, lat_1, lon_inp(Ielem,Iline),  &
               Lat_South, Lat_North, dlat, Nlat, Lon_West, Lon_East,dlon, Nlon, xdist_01)

            call POSITION_TO_INDICES (Ilat_2, Ilon, lat_2, lon_inp(Ielem,Iline),  &
               Lat_South, Lat_North, dlat, Nlat, Lon_West, Lon_East,dlon, Nlon, xdist_02)
       
            Ilon_min = min(Ilon_1,Ilon_2)
            Ilon_max = max(Ilon_1,Ilon_2)
            Ilat_min = min(Ilat_1,Ilat_2)
            Ilat_max = max(Ilat_1,Ilat_2)
            
            if (Ilon_min .LE. 1) cycle
            if (Ilat_min .LE. 1) cycle
            
            do Ilon_1 = Ilon_min, Ilon_max
               do Ilat_1 = Ilat_min, Ilat_max
                  !-- do not reset grid point closest to pixel
                  if ((Ilon_1 == Ilon) .and. (Ilat_1 == Ilat)) cycle
                  !-- fill in values for all grid points within pixel footprint
                  !-- ignore those already set previously
                  if (Min_Dist_Grid(Ilon_1,Ilat_1) == big_number .and.  &
                           .not. Gap_Inp(Ielem,Iline)) then
                     Ielem_Out(Ilon_1,Ilat_1)  = Ielem
                     Iline_Out(Ilon_1,Ilat_1)  = Iline
                  end if
               end do
            end do
            
           !  print*,'L2: ',Ielem, Iline,lon_inp(Ielem,Iline),lat_inp(Ielem,Iline), Ielem_Out(Ilon,Ilat),Iline_Out(Ilon,Ilat)
         end do
      end do

      DEALLOCATE(Min_Dist_Grid)
      !--- make grid lat and lons
 

   end subroutine INDEX_IN_REGULAR_GRID
   
   ! ----------------------- !
   !                         ! 
   ! ----------------------- ! 
   subroutine lonlat_edge_to_array (  &
         Lat_South &
        , dlat      &
        , Nlat      &
        , Lon_West  &
        , dlon      &
        , Nlon      &
        , lat_out &
        , lon_out )
      
      real(kind=real4), intent(in):: Lat_South
      real(kind=real4), intent(in):: dlat
      integer(kind=int4), intent(in):: Nlat
      real(kind=real4), intent(in):: Lon_West
      real(kind=real4), intent(in):: dlon
      integer(kind=int4), intent(in):: Nlon
      real, intent(out) :: lat_out(:,:)
      real, intent(out) :: lon_out(:,:)
      
      
      integer :: ilat
      integer :: ilon
      
      do Ilat = 1, Nlat
         do Ilon = 1, Nlon
            Lat_Out(Ilon,Ilat) = Lat_South + (Ilat-1)*dlat
            Lon_Out(Ilon,Ilat) = Lon_West + (Ilon-1)*dlon
         end do
      end do

      where(Lon_Out > 180.0)
         Lon_Out = Lon_Out - 360.0
      end where
   
   
   end subroutine lonlat_edge_to_array 
   


!====================================================================
! SUBROUTINE Name: POSITION_TO_INDICES
!
! Function: 
!       Converts lat/lon position to grid indicies 
!
! Description:
!   This subroutine is a routine to convert lat/lon position to grid index
!
!====================================================================

subroutine POSITION_TO_INDICES(Ilat,Ilon, lat, lon,  &
                               Lat_South, Lat_North, &
                               dlat, Nlat, Lon_West, Lon_East, &
                               dlon, Nlon, distance)

       real(kind=real4), intent(in):: lat
       real(kind=real4), intent(in):: lon
       real(kind=real4), intent(in):: Lat_South
       real(kind=real4), intent(in):: Lat_North
       real(kind=real4), intent(in):: dlat
       integer(kind=int4), intent(in):: Nlat
       real(kind=real4), intent(in):: Lon_West
       real(kind=real4), intent(in):: Lon_East
       real(kind=real4), intent(in):: dlon
       integer(kind=int4), intent(in):: Nlon
       integer(kind=int4), intent(out):: Ilat
       integer(kind=int4), intent(out):: Ilon
       real(kind=real4), intent(out):: distance
       real(kind=real4):: lon_abs
       real(kind=real4):: Lon_East_abs
       real(kind=real4):: lat_abs

       Ilat = -999
       Ilon = -999

      Lon_East_abs = Lon_East
      if (Lon_East_abs < 0.0) Lon_East_abs = Lon_East_abs + 360.0
      lon_abs = lon
      if (lon_abs < 0.0) lon_abs = lon_abs + 360.0

      if ((lat >= Lat_South) .and. (lat <= Lat_North)) then
         Ilat = nint((lat - Lat_South) / dlat) + 1
         Ilat = min(Nlat,max(1,Ilat))
      endif

     !-- -no date line
     if (Lon_West  < Lon_East) then
      if ((lon >= Lon_West) .and. (lon <= Lon_East)) then
        Ilon = nint((lon - Lon_West) / dlon) + 1
        Ilon = min(Nlon,max(1,Ilon))
      end if
     else 
      if (lon < 0.0) lon_abs = lon + 360.0
      if ((lon_abs <= Lon_East_abs) .and. (lon_abs >= Lon_West)) then
        Ilon = nint((lon_abs - Lon_West) / dlon) + 1
        Ilon = min(Nlon,max(1,Ilon))
      endif
     end if

     distance = 999.9
     if ((Ilon >= 1) .and. (Ilat >= 1)) then
             lat_abs = Lat_South + (Ilat-1)*dlat
             Lon_East_abs = Lon_West + (Ilon-1)*dlon
             lon_abs = lon
             if (Lon_West > 0.0 .and. lon < 0.0) lon_abs = lon + 360.0
             distance = (lat-lat_abs)**2 + (lon_abs-Lon_East_abs)**2
     end if

end subroutine POSITION_TO_INDICES




end module cx_grid_tools_mod
