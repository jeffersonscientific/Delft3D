!----- AGPL --------------------------------------------------------------------
!
!  Copyright (C)  Stichting Deltares, 2017-2024.
!
!  This file is part of Delft3D (D-Flow Flexible Mesh component).
!
!  Delft3D is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Affero General Public License as
!  published by the Free Software Foundation version 3.
!
!  Delft3D  is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Affero General Public License for more details.
!
!  You should have received a copy of the GNU Affero General Public License
!  along with Delft3D.  If not, see <http://www.gnu.org/licenses/>.
!
!  contact: delft3d.support@deltares.nl
!  Stichting Deltares
!  P.O. Box 177
!  2600 MH Delft, The Netherlands
!
!  All indications and logos of, and references to, "Delft3D",
!  "D-Flow Flexible Mesh" and "Deltares" are registered trademarks of Stichting
!  Deltares, and remain the property of Stichting Deltares. All rights reserved.
!
!-------------------------------------------------------------------------------

module m_his_file

   use stdlib_kinds, only: dp
   use MessageHandling, only: mess, err, LEVEL_ERROR, LEVEL_DEBUG, LEVEL_WARN
   use dfm_error, only: DFM_NOERR
   use fm_statistical_output, only: model_is_3D
   use unstruc_netcdf, only: ihisfile, definencvar
   use m_flowtimes, only: handle_extra, it_his, ti_his, ti_hiss, ti_hise, time_his
   use netcdf_utils, only: ncu_set_att, check_netcdf_error
   use timers, only: timon, timstrt, timstop
   use m_missing, only: dmiss
   use m_his_file_netcdf_ids

   implicit none
   
   private
   
   public :: def_and_put_his_file
   
   integer              :: strlen_netcdf         !< string length definition for (station) names on history file
   integer              :: nc_precision          !< Precision of NetCDF variables (e.g. nf90_int, nf90_double, etc.)
   
   real(dp)             :: curtime_split = 0d0   !< Current time-partition that the file writer has open.
   real(dp)             :: time_his_prev         !< Last time that data was written to the his file
   
   logical              :: add_latlon = .false.  !< Whether or not to include station lat/lon coordinates in the his file
   integer              :: jawrizc    = 0        !< Whether or not to write zc coordinates to the his file
   integer              :: jawrizw    = 0        !< Whether or not to write zw coordinates to the his file
   
   character(len=1024)  :: statcoordstring       !< String listing the coordinates for each variables
   
   integer, allocatable :: id_hwq(:)
   
contains

!> Create and write to the his file
subroutine def_and_put_his_file(time)
   use fm_statistical_output, only: model_has_obs_stations
   use m_flowparameters, only: jahiszcor
   use m_map_his_precision, only: md_nc_his_precision, netcdf_data_type
   use netcdf, only: nf90_inquire, nf90_sync, nf90_enddef
   use unstruc_netcdf, only: unc_close, unc_writeopts, UG_WRITE_LATLON, unc_noforcedflush
   use m_monitoring_crosssections, only: ncrs
   use m_monitoring_runupgauges, only: num_rugs
   use m_flowparameters, only: jahiscgen, jahisbal
   use m_sferic, only: jsferic
   use m_his_file_statistical, only: def_his_file_statoutput, put_his_file_statoutput
   use m_his_file_structures, only: def_his_file_time_independent_structures, put_his_file_time_independent_structures

   double precision, intent(in) :: time !< Current time, should in fact be time1, since the data written is always s1, ucx, etc.

   integer                      :: ierr, ndims

   if (jahiszcor > 0) then
      jawrizc = 1
      jawrizw = 1
   end if

   nc_precision = netcdf_data_type(md_nc_his_precision)

   if (timon) then 
      call timstrt ( "def_and_put_his_file", handle_extra(54))
   end if

   ! Close/reset any previous hisfile.
   if (ihisfile /= 0) then  ! reset stord ncid to zero if file not open
      ierr = nf90_inquire(ihisfile, ndims)
      if (ierr /= 0) ihisfile = 0
   end if

   if (ihisfile > 0 .and. it_his == 0) then
      ierr = unc_close(ihisfile)
      ihisfile = 0
   end if

   ! When no crs/obs present, return immediately.
   if (model_has_obs_stations() == .false. .and. ncrs <= 0 .and. jahisbal <= 0 .and. jahiscgen <= 0 .and. num_rugs <= 0) then
      if (ihisfile == 0) then
         call mess(LEVEL_WARN, 'No observations nor cross sections defined. Will not produce a history file.')
      end if
      ihisfile = -1 ! -1 stands for: no file open, no obs/crs defined.
      return
   end if

   ! Only add auto-tranformed lat/lon coordinates if model is Cartesian and user has requested extra latlon output.
#ifdef HAVE_PROJ
   add_latlon = jsferic == 0 .and. iand(unc_writeopts, UG_WRITE_LATLON) == UG_WRITE_LATLON
#else
   add_latlon = .false.
#endif

   if (ihisfile == 0) then
   
      if (timon) then
         call timstrt ( "def_and_put_his_file INIT/DEF", handle_extra(61))
      end if
   
      ! Create the his file and define all its basic time-independent dimensions and variables
      call def_his_file_time_independent_basic()

      ! Define all the time-independent dimensions and variables for structures
      if (timon) then
         call timstrt ( "def_and_put_his_file DEF structures", handle_extra(60))
      end if
      call def_his_file_time_independent_structures(nc_precision, add_latlon, jawrizc, jawrizw, statcoordstring)
      if (timon) then
         call timstop (handle_extra(60))
      end if
   
      ! Define all user-requested statistical output variables
      call def_his_file_statoutput(statcoordstring)

      ! Switch from define mode to data mode
      call check_netcdf_error(nf90_enddef(ihisfile))
       
      ! Write all time-independent data to the his file
      if (it_his == 0) then
         if (timon) then
            call timstrt('def_and_put_his_file timeindep data', handle_extra(63))
         end if
         call put_his_file_time_independent_structures()
         if (timon) then
            call timstop( handle_extra(63))
         end if
      end if
   
      if (timon) then
         call timstop(handle_extra(61))
      end if
      
   end if
    
   ! Increment output counters in m_flowtimes.
   if (it_his == 0) then
      time_his_prev = time
   end if
   time_his = time
   it_his   = it_his + 1
    
   ! Write all time-independent data to the his file
   if (timon) then
      call timstrt('def_and_put_his_file time data', handle_extra(64))
   end if
   call put_his_file_time_dependent_basic()
   if (timon) then
      call timstop( handle_extra(64))
   end if
   
   ! Write all user-requested statistical output variables
   call put_his_file_statoutput()

   if (unc_noforcedflush == 0) then
      call check_netcdf_error( nf90_sync(ihisfile)) ! Flush file
   end if

   if (timon) then
      call timstop (handle_extra(54))
   end if

end subroutine def_and_put_his_file

!> Create the his file and define all its basic time-independent dimensions and variables
subroutine def_his_file_time_independent_basic()
   use m_flowtimes, only: ti_split, time_split0, Tudunitstr
   use unstruc_files, only: defaultFilename
   use unstruc_netcdf, only: unc_create, unc_meta_add_user_defined, unc_add_time_coverage, &
                             unc_def_var_nonspatial
   use netcdf, only: nf90_noerr, nf90_def_dim, nf90_unlimited, nf90_char, nf90_double, nf90_int
   use MessageHandling, only: idlen
   use m_flow, only: kmx
   use m_sediment, only: stm_included, stmpar, jased
   use m_flowparameters, only: jahissed, jacheckmonitor, jahisbedlev
   use m_transport, only: ISED1
   use m_ug_nc_attribute, only: ug_nc_attribute
   use m_fm_wq_processes, only: jawaqproc
   use fm_statistical_output, only: model_has_obs_stations

   character(len=255)    :: filename
   integer               :: ierr
   type(ug_nc_attribute) :: attributes(4)
   
   ! Another time-partitioned file needs to start, reset iteration count (and file).
   if (ti_split > 0d0 .and. curtime_split /= time_split0) then
      it_his        = 0
      curtime_split = time_split0
   end if

   if (ti_split > 0d0) then
      filename = defaultFilename('his', timestamp=time_split0)
   else
      filename = defaultFilename('his')
   end if

   ierr = unc_create(filename, 0, ihisfile)
   if (ierr /= nf90_noerr) then
      call mess(LEVEL_WARN, 'Could not create history file.')
   end if

   ierr = unc_meta_add_user_defined(ihisfile)
   ierr = unc_add_time_coverage(ihisfile, ti_hiss, ti_hise, ti_his)
   
   ! General purpose dimensions
   !if (unc_nounlimited > 0) then ! UNST-4764: His file has never shown good results with NcNoUnlimited option on.
   call check_netcdf_error( nf90_def_dim(ihisfile, 'time', nf90_unlimited, id_timedim))
   call check_netcdf_error( nf90_def_dim(ihisfile, 'two', 2, id_twodim))

   strlen_netcdf = idlen  !< Max string length of Ids.
   call check_netcdf_error( nf90_def_dim(ihisfile, 'name_len', strlen_netcdf, id_strlendim))

   if (model_is_3D()) then
      call check_netcdf_error( nf90_def_dim(ihisfile, 'laydim', kmx, id_laydim))
      call check_netcdf_error( nf90_def_dim(ihisfile, 'laydimw', kmx+1, id_laydimw))
   end if

   if (stm_included .and. jahissed > 0) then
      if (ISED1 > 0) then
         ! New implementation, sedsus fraction is additional dimension
         call check_netcdf_error( nf90_def_dim(ihisfile, 'nSedTot', stmpar%lsedtot, id_sedtotdim))
         call check_netcdf_error( nf90_def_dim(ihisfile, 'nSedSus', stmpar%lsedsus, id_sedsusdim))
         call definencvar(ihisfile, id_frac_name, nf90_char, (/ id_strlendim, id_sedtotdim /), 'sedfrac_name', 'sediment fraction identifier')
      end if
      if (jased > 0 .and. stmpar%morlyr%settings%iunderlyr==2) then
         call check_netcdf_error( nf90_def_dim(ihisfile, 'nBedLayers', stmpar%morlyr%settings%nlyr, id_nlyrdim))
      end if
   end if
   
   ! Time
   call ncu_set_att(attributes(1), 'standard_name', 'time')
   call ncu_set_att(attributes(2), 'bounds', 'time_bds')
   call definencvar(ihisfile, id_time, nf90_double, (/ id_timedim /), 'time', unit=trim(Tudunitstr),extra_attributes=attributes(1:2))
   call definencvar(ihisfile, id_timebds, nf90_double,(/ id_twodim, id_timedim /), 'time_bds', 'Time interval for each point in time.', unit=trim(Tudunitstr),extra_attributes=attributes(1:1))

   ! Size of latest timestep
   ierr = unc_def_var_nonspatial(ihisfile, id_timestep, nf90_double, (/ id_timedim /), 'timestep', '',     'latest computational timestep size in each output interval', 's')

   if (jacheckmonitor == 1) then
      call definencvar(ihisfile, id_checkmon, nc_precision, (/ id_laydim, id_timedim /), 'checkerboard_monitor','Checkerboard mode monitor', unit='m s-1')
      call definencvar(ihisfile, id_num_timesteps, nf90_int, (/ id_timedim /), 'num_timesteps')
      call definencvar(ihisfile, id_comp_time, nc_precision, (/ id_timedim /), 'comp_time')
   end if

   ! WAQ statistic outputs are kept outside of the statistical output framework
   if (jawaqproc <= 0) then
      ierr = def_his_file_station_waq_statistic_outputs(statcoordstring, nc_precision, id_hwq)
   end if

   ! Bottom level is written separately from statout if it is static
   if (jahisbedlev > 0 .and. model_has_obs_stations() .and. .not. stm_included ) then
      call definencvar(ihisfile, id_varb, nc_precision, [id_statdim], 'bedlevel', 'bottom level', unit='m', namecoord=statcoordstring)
   end if
       
end subroutine def_his_file_time_independent_basic

!> Write all time-dependent data to the his file
subroutine put_his_file_time_dependent_basic()
   use netcdf, only: nf90_put_var
   use m_flowtimes, only: dts, dnt, handle_steps
   use m_flowexternalforcings, only: numsrc, qsrc, qstss
   use m_transport, only: numconst
   use fm_statistical_output, only: model_has_obs_stations
   use m_sediment, only: stm_included
   use m_flowparameters, only: jahisbedlev, jacheckmonitor
   use m_observations, only: valobs, IPNT_BL
   use m_fm_wq_processes, only: jawaqproc
   use m_globalparameters, only: ST_MAX_TYPE
   use unstruc_netcdf, only: unc_write_struc_input_coordinates
   use m_filter, only: checkmonitor
   use timers, only: tim_get_wallclock
   use m_his_file_structures, only: put_his_file_station_coord_vars
   
   integer :: i, ierr, n
   
   call check_netcdf_error( nf90_put_var(ihisfile, id_time, time_his, (/ it_his /)))
   call check_netcdf_error( nf90_put_var(ihisfile, id_timebds, (/ time_his_prev, time_his /), (/ 1, it_his /)))
   time_his_prev = time_his
   call check_netcdf_error( nf90_put_var(ihisfile, id_timestep, dts, (/ it_his /)))

   ! Fill average source-sink discharge with different array on first timestep
   if (it_his == 1) then
      do i = 1, numsrc
         qsrc(i) = qstss((numconst+1)*(i-1)+1)
      end do
   end if
   
   ! Bottom level is written separately from statout if it is static
   if (jahisbedlev > 0 .and. model_has_obs_stations() .and. .not. stm_included ) then
      call check_netcdf_error( nf90_put_var(ihisfile, id_varb, valobs(:,IPNT_BL), start = (/ 1 /) ))
   end if

   ! WAQ statistic outputs are kept outside of the statistical output framework
   if (.not. model_has_obs_stations() .or. jawaqproc <= 0) then
      ierr = put_his_file_station_waq_statistic_outputs(id_hwq)
   end if

   ! Write x/y-, lat/lon- and z-coordinates for the observation stations every time (needed for moving observation stations)
   ierr = put_his_file_station_coord_vars(add_latlon, jawrizc, jawrizw, &
                                       id_statx, id_staty, id_statlat, id_statlon, &
                                       id_zcs, id_zws, id_zwu, it_his, &
                                       id_statgeom_node_count, id_statgeom_node_coordx, id_statgeom_node_coordy, &
                                       id_statgeom_node_lon, id_statgeom_node_lat)

   if (it_his == 1) then
      do n = 1, ST_MAX_TYPE
         call unc_write_struc_input_coordinates(ihisfile, n)
      end do
   end if

   if (jacheckmonitor == 1) then
      call check_netcdf_error( nf90_put_var(ihisfile, id_checkmon, checkmonitor, start = [1, it_his]))
      call check_netcdf_error( nf90_put_var(ihisfile, id_num_timesteps, int(dnt), start = [it_his]))
      call check_netcdf_error( nf90_put_var(ihisfile, id_comp_time, tim_get_wallclock(handle_steps), start = [it_his]))
   end if
    
end subroutine put_his_file_time_dependent_basic

!> Define variables for WAQ statistic outputs (not to be confused with the general statistical output framework).
!! They are not part of the statistical output framework, because it is useless to allow statistics about statistics,
!! because writing output only in the final time step is currently not supported in statistical output,
!! and because statistic outputs may be redundant in the future when the statistical output framework is feature complete.
function def_his_file_station_waq_statistic_outputs(statcoordstring, nc_precision, waq_statistics_ids) result(ierr)
   use m_fm_wq_processes, only: wq_user_outputs => outputs, noout_statt, noout_state, noout_user
   use netcdf, only: nf90_put_att
   use string_module, only: replace_multiple_spaces_by_single_spaces

   character(len=1024),  intent(in   ) :: statcoordstring       !< String listing the coordinate variables associated with the stations
   integer,              intent(in   ) :: nc_precision          !< Precision of NetCDF variables (e.g. nf90_int, nf90_double, etc.)
   integer, allocatable, intent(  out) :: waq_statistics_ids(:) !< NetCDF ids for the water quality statistic output variables
   integer                             :: ierr                  !< D-Flow FM error code

   character(len = 255)  :: variable_name, description
   character(len = 1024) :: station_coordinate_string
   integer               :: statistics_index
   integer, allocatable  :: nc_dimensions(:), specific_nc_dimensions(:)

   ierr = DFM_NOERR

   allocate(waq_statistics_ids(noout_statt + noout_state))

   if (model_is_3D()) then
      nc_dimensions = [id_laydim, id_statdim, id_timedim]
      station_coordinate_string = trim(statcoordstring) // ' zcoordinate_c'
   else
      nc_dimensions = [id_statdim, id_timedim]
      station_coordinate_string = statcoordstring
   end if

   do statistics_index = 1, noout_statt + noout_state
      if (statistics_index > noout_statt) then
         specific_nc_dimensions = nc_dimensions(1 : size(nc_dimensions) - 1) ! Drop time dimension for end statistics (stat-e) variables
      else
         specific_nc_dimensions = nc_dimensions
      end if
      variable_name = ' '
      write (variable_name, "('water_quality_stat_',I0)") statistics_index
      call definencvar(ihisfile, waq_statistics_ids(statistics_index), nc_precision, specific_nc_dimensions, &
                        trim(variable_name), trim(wq_user_outputs%names(statistics_index)), &
                        trim(wq_user_outputs%units(statistics_index)), trim(station_coordinate_string), 'station_geom', fillVal = dmiss)
      description = trim(wq_user_outputs%names(statistics_index))//' - '//trim(wq_user_outputs%description(statistics_index))//' in flow element'
      call replace_multiple_spaces_by_single_spaces(description)
      call check_netcdf_error(nf90_put_att(ihisfile, waq_statistics_ids(statistics_index), 'description', description))
   end do
end function def_his_file_station_waq_statistic_outputs

!> Write data to WAQ statistic output variables (not to be confused with the general statistical output framework).
function put_his_file_station_waq_statistic_outputs(waq_statistics_ids) result(ierr)
   use m_observations, only: numobs, nummovobs, IPNT_HWQ1, valobs
   use m_flow, only: kmx
   use m_fm_wq_processes, only: noout_statt, noout_state, noout_user
   use precision_basics, only: comparereal
   use m_flowparameters, only: eps10
   use netcdf, only: nf90_put_var

   integer,          intent(in) :: waq_statistics_ids(:) !< NetCDF ids for the water quality statistic output variables
   integer                      :: ierr                  !< D-Flow FM error code

   integer               :: ntot
   integer, allocatable  :: nc_start(:), nc_count(:)
   integer               :: start_index_valobs, statistics_index, num_layers

   ierr = DFM_NOERR

   ! Default start and count for stat-t variables
   ntot = numobs + nummovobs
   if (model_is_3D()) then
      nc_start = [1, 1, it_his]
      nc_count = [kmx, ntot, 1]
   else
      nc_start = [1, it_his]
      nc_count = [ntot, 1]
   end if

   do statistics_index = 1, noout_statt + noout_state
      if (statistics_index == noout_statt + 1 ) then
         if (comparereal(time_his, ti_hise, eps10) < 0) then
            return ! The end statistic outputs (stat-e) are only written in the last his time step
         end if
         if (model_is_3D()) then
            nc_start = [1, 1]
            nc_count = [kmx, ntot]
         else
            nc_start = [1]
            nc_count = [ntot]
         end if
      end if

      num_layers = max(kmx, 1)

      start_index_valobs = IPNT_HWQ1 - 1 + (noout_user + statistics_index - 1) * num_layers + 1
      call check_netcdf_error(nf90_put_var(ihisfile, waq_statistics_ids(statistics_index), &
                                           transpose(valobs(:, start_index_valobs : start_index_valobs + num_layers - 1)), &
                                           start = nc_start, count = nc_count))
   end do
end function put_his_file_station_waq_statistic_outputs

end module m_his_file