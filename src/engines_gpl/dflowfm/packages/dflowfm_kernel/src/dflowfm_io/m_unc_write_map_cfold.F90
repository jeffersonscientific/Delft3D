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
   
module m_unc_write_map_cfold
   
   use stdlib_kinds, only: dp
   
   implicit none
   
   private
   
   public :: unc_write_map_filepointer
   
   contains
   
!> Writes map/flow data to an already opened netCDF dataset.
!! The netnode and -links have been written already.
subroutine unc_write_map_filepointer(imapfile, tim, jaseparate) ! wrimap
    use m_flow
    use m_flowtimes
    use m_flowgeom
    use m_sobekdfm
    use m_heatfluxes
    use m_sferic
    use network_data
    use m_sediment
    use m_bedform
    use m_wind
    use m_flowparameters, only: jatrt, jacali
    use m_mass_balance_areas
    use m_fm_wq_processes
    use m_xbeach_data
    use m_transportdata
    use bedcomposition_module, only: bedcomp_getpointer_integer
    use m_alloc
    use m_missing
    use string_module, only: replace_multiple_spaces_by_single_spaces
    use netcdf_utils, only: ncu_append_atts
    use m_fm_icecover, only: ice_mapout, ice_af, ice_h, ice_p, ice_t, snow_h, snow_t, ja_icecover, ICECOVER_SEMTNER
    use netcdf, only: nf90_inquire, nf90_inq_dimid, nf90_def_dim, nf90_unlimited, nf90_def_var, nf90_double, nf90_put_att, &
                      nf90_def_var, nf90_int, nf90_noerr, nf90_global, nf90_inq_varid, nf90_inquire_variable, nf90_enddef, &
                      nf90_put_var
    use unstruc_netcdf, only: unc_write_flowgeom_filepointer, unc_add_gridmapping_att, unc_append_3dflowgeom_put, &
                              unc_write_net_filepointer, definencvar, unc_append_3dflowgeom_def, check_error
    use string_module, only: replace_char

    implicit none

    integer,           intent(in) :: imapfile
    real(kind=hp),     intent(in) :: tim
    integer, optional, intent(in) :: jaseparate   !< Whether this save is manual/by user (not part of the standard map write series)

    integer                       :: jaseparate_, idims(2)

    logical, dimension(2), save   :: firststep = .true.

    integer, save                 :: ierr, ndim
    integer, dimension(2), save   :: &
    !id_netcelldim, id_netcellmaxnodedim, id_netcellcontourptsdim, &
    id_laydim, id_wdim, &
        id_flowelemdim, &
    id_maxfracdim,  &
    id_erolaydim,   &
    id_flowlinkdim, &
    id_netlinkdim,  &
    id_1d2ddim,     &
    id_timedim,     &
    id_time, id_timestep, &
    id_sbcx, id_sbcy, id_sbcx_reconstructed, id_sbcy_reconstructed, &
    id_sbwx, id_sbwy, id_sbwx_reconstructed, id_sbwy_reconstructed, &
    id_sswx, id_sswy, id_sswx_reconstructed, id_sswy_reconstructed, &
    id_sourse, id_sinkse, id_ws, &
    id_sxtot, id_sytot, id_rsedeq, id_umod, id_zumod, id_ustar, id_dzdn, id_dzdt, id_morbl, id_aks, id_rca, &
    id_bodsed, id_dpsed, id_msed, id_lyrfrac, id_thlyr, id_poros, id_nlyrdim, &
    id_sedtotdim, id_sedsusdim, id_rho, id_rhop, id_viu, id_diu, id_q1, id_spircrv, id_spirint, &
    id_q1main, &
    id_s1, id_taus, id_ucx, id_ucy, id_ucz, id_ucxa, id_ucya, id_unorm, id_ww1, id_sa1, id_tem1, id_sed, id_ero, id_s0, id_u0, id_cfcl, id_cftrt, id_czs, id_czu, &
    id_qsun, id_qeva, id_qcon, id_qlong, id_qfreva, id_qfrcon, id_qtot, &
    id_patm, id_ice_af, id_ice_h, id_ice_p, id_ice_t, id_snow_h, id_snow_t, id_tair, id_rhum, id_clou, id_E, id_R, id_H, id_D, id_DR, id_urms, id_thetamean, &
    id_cwav, id_cgwav, id_sigmwav, &
    id_ust, id_vst, id_windx, id_windy, id_windxu, id_windyu, id_numlimdt, id_hs, id_bl, id_zk, &
    id_1d2d_edges, id_1d2d_zeta1d, id_1d2d_crest_level, id_1d2d_b_2di, id_1d2d_b_2dv, id_1d2d_d_2dv, id_1d2d_q_zeta, id_1d2d_q_lat, &
    id_1d2d_cfl, id_1d2d_flow_cond, id_1d2d_sb, id_1d2d_s1_2d, id_1d2d_s0_2d, id_tidep, id_salp, id_inttidesdiss, &
    id_duneheight, id_dunelength, id_ksd, id_ksr, id_ksmr, id_ks, &
    id_taurat, id_dm, id_dg, id_dgsd, id_frac, id_mudfrac, id_sandfrac, id_fixfac, id_hidexp, id_mfluff, id_scrn, id_urmscc, id_Fxcc, id_Fycc, &
    id_sscx, id_sscy, id_sscx_reconstructed, id_sscy_reconstructed, &
    id_turkin1, id_tureps1, id_vicwwu, id_vicwws, id_swanbl, &
    id_rnveg, id_diaveg, id_veg_stemheight

    integer,          dimension(:,:),   allocatable, save :: id_dxx                     ! fractions
    double precision, dimension(:),     allocatable       :: dum
    double precision, dimension(:,:),   allocatable       :: poros
    double precision, dimension(:,:,:), allocatable       :: frac
    double precision, dimension(:),     allocatable       :: toutput
    double precision, dimension(:,:),   allocatable       :: toutputx, toutputy
    double precision, dimension(:),     allocatable       :: rks

    integer,          dimension(:), allocatable :: idum

    integer :: iid, i, j, jj, itim, k, kb, kt, kk, n, LL, Ltx, Lb, L, nm, nlayb,nrlay, nlaybL, nrlayLx, varid, ndims
    integer :: ndxndxi ! Either ndx or ndxi, depending on whether boundary nodes also need to be written.
    double precision, dimension(:), allocatable :: windx, windy
    double precision, dimension(:), allocatable :: numlimdtdbl ! TODO: WO/AvD: remove this once integer version of unc_def_map_var is available
    double precision :: vicc, dicc
    integer :: jaeulerloc

    double precision   :: rhol
    character(16)      :: dxname, zw_elem, zcc_elem, zwu_link, zu_link
    character(64)      :: dxdescr
    character(10)      :: transpunit
    character(len=255) :: tmpstr

    integer, dimension(:), allocatable :: flag_val
    character(len=10000)               :: flag_mean

    double precision, dimension(:), pointer :: dens

    if (.not. allocated(id_dxx) .and. stm_included) allocate(id_dxx(1:stmpar%morpar%nxx,1:2))

    ! If jaseparate_==1 or this map file was just opened for the first time:
    ! only write net+vardefs first time, and write subsequent flow snapshots in later calls.
    ! jaseparate_==2: write com file
    if (present(jaseparate)) then
        jaseparate_ = jaseparate
    else
        jaseparate_ = 0
    endif

    if (jaseparate_ == 0 .or. jaseparate_ == 1) then
       ! mapfile, store/use ids number 1
       iid = 1
       ndxndxi = ndxi
    elseif (jaseparate_ == 2) then
       ! comfile, store/use ids number 2
       iid = 2
       ndxndxi = ndx ! Com file, include boundary nodes
    else
       ! error
       iid = 0
    endif

    ! Use nr of dimensions in netCDF file a quick check whether vardefs were written
    ! before in previous calls.
    ndim = 0
    ierr = nf90_inquire(imapfile, nDimensions=ndim)

    ! Only write net and flow geometry data the first time, or for a separate map file.
    if (ndim == 0) then

        call unc_write_net_filepointer(imapfile)      ! Write standard net data as well

        if (jaseparate_ == 2) then
           call unc_write_flowgeom_filepointer(imapfile, jabndnd = 1) ! Write time-independent flow geometry data, with boundary nodes
           ierr = nf90_inq_dimid(imapfile, 'nFlowElemWithBnd', id_flowelemdim(iid))
        else
           call unc_write_flowgeom_filepointer(imapfile) ! Write time-independent flow geometry data
           ierr = nf90_inq_dimid(imapfile, 'nFlowElem', id_flowelemdim(iid))
        endif

        ierr = nf90_inq_dimid(imapfile, 'nFlowLink', id_flowlinkdim(iid))
        ierr = nf90_inq_dimid(imapfile, 'nNetLink' , id_netlinkdim(iid))

        if (nbnd1d2d > 0) then
           ierr = nf90_def_dim(imapfile, 'nBnd1d2d', nbnd1d2d, id_1d2ddim(iid))
        endif

        ! Time
        ierr = nf90_def_dim(imapfile, 'time', nf90_unlimited, id_timedim(iid))
        call check_error(ierr, 'def time dim')
        ierr = nf90_def_var(imapfile, 'time', nf90_double, id_timedim(iid),  id_time(iid))
        ierr = nf90_put_att(imapfile, id_time(iid),  'units'        , trim(Tudunitstr))
        ierr = nf90_put_att(imapfile, id_time(iid),  'standard_name', 'time')

        ! 3D
        if ( kmx > 0 ) then
           call unc_append_3dflowgeom_def(imapfile)              ! Append definition of time-independent 3d flow geometry data
           ierr = nf90_inq_dimid(imapfile, 'laydim', id_laydim(iid))
           ierr = nf90_inq_dimid(imapfile, 'wdim', id_wdim(iid))
        endif

        ! Size of latest timestep
        ierr = nf90_def_var(imapfile, 'timestep', nf90_double, id_timedim(iid),  id_timestep(iid))
        ierr = nf90_put_att(imapfile, id_timestep(iid),  'units'        , 'seconds')
        ierr = nf90_put_att(imapfile, id_timestep(iid),  'standard_name', 'timestep')

        if (jamaps1 > 0 .or. jaseparate_==2) then
            ! Flow data on centres: water level at latest timestep
            ierr = nf90_def_var(imapfile, 's1',  nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_s1(iid))
            ierr = nf90_put_att(imapfile, id_s1(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
            ierr = nf90_put_att(imapfile, id_s1(iid),   'standard_name', 'sea_surface_height') ! sorry for inland water people
            ierr = nf90_put_att(imapfile, id_s1(iid),   'long_name'    , 'water level')
            ierr = nf90_put_att(imapfile, id_s1(iid),   'units'        , 'm')
            ierr = unc_add_gridmapping_att(imapfile, (/ id_s1(iid) /), jsferic)
        endif

        if (jaseparate_ == 0 .or. jaseparate_ == 1) then ! to mapfile
            ! Flow data on centres: water level timestep before the latest timestep

            if (jamaps0 > 0) then
                ierr = nf90_def_var(imapfile, 's0',  nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_s0(iid))
                ierr = nf90_put_att(imapfile, id_s0(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                ierr = nf90_put_att(imapfile, id_s0(iid),   'standard_name', 'sea_surface_height') ! sorry for inland water people
                ierr = nf90_put_att(imapfile, id_s0(iid),   'long_name'    , 'water level at previous timestep')
                ierr = nf90_put_att(imapfile, id_s0(iid),   'units'        , 'm')
                ierr = unc_add_gridmapping_att(imapfile, (/ id_s0(iid) /), jsferic)
            endif

            idims(1) = id_flowelemdim(iid)
            idims(2) = id_timedim(iid)

            if (jamaphs > 0) then
                call definencvar(imapfile,id_hs(iid)   ,nf90_double,idims, 'waterdepth'  , 'water depth', 'm', 'FlowElem_xcc FlowElem_ycc')
            endif

            if (jamapheatflux > 0 .and. jatem > 1) then ! Heat modelling only
               call definencvar(imapfile,id_tair(iid)   ,nf90_double,idims, 'Tair'  , 'air temperature', 'degC', 'FlowElem_xcc FlowElem_ycc')
               call definencvar(imapfile,id_rhum(iid)   ,nf90_double,idims, 'rhum'  , 'Relative humidity', ' ','FlowElem_xcc FlowElem_ycc')
               call definencvar(imapfile,id_clou(iid)   ,nf90_double,idims, 'clou'  , 'cloudiness', ' ', 'FlowElem_xcc FlowElem_ycc')

               if (jatem == 5) then
                  call definencvar(imapfile,id_qsun(iid)   ,nf90_double,idims, 'Qsun'  , 'solar influx', 'W m-2', 'FlowElem_xcc FlowElem_ycc')
                  call definencvar(imapfile,id_Qeva(iid)   ,nf90_double,idims, 'Qeva'  , 'evaporative heat flux', 'W m-2', 'FlowElem_xcc FlowElem_ycc')
                  call definencvar(imapfile,id_Qcon(iid)   ,nf90_double,idims, 'Qcon'  , 'sensible heat flux', 'W m-2', 'FlowElem_xcc FlowElem_ycc')
                  call definencvar(imapfile,id_Qlong(iid)  ,nf90_double,idims, 'Qlong' , 'long wave back radiation', 'W m-2', 'FlowElem_xcc FlowElem_ycc')
                  call definencvar(imapfile,id_Qfreva(iid) ,nf90_double,idims, 'Qfreva', 'free convection evaporative heat flux', 'W m-2', 'FlowElem_xcc FlowElem_ycc')
                  call definencvar(imapfile,id_Qfrcon(iid) ,nf90_double,idims, 'Qfrcon', 'free convection sensible heat flux', 'W m-2', 'FlowElem_xcc FlowElem_ycc')
               endif

               call definencvar(imapfile,id_Qtot(iid)   ,nf90_double,idims, 'Qtot'  , 'total heat flux', 'W m-2', 'FlowElem_xcc FlowElem_ycc')
            endif

            if (jamapnumlimdt > 0) then
                call definencvar(imapfile,id_numlimdt(iid)  ,nf90_double,idims, 'numlimdt' , 'number of times flow element was Courant limiting', '1', 'FlowElem_xcc FlowElem_ycc')
            endif

            if (jamaptaucurrent > 0) then
                ! Flow data on centres
                ierr = nf90_def_var(imapfile, 'taus' ,  nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_taus(iid))
                ierr = nf90_put_att(imapfile, id_taus(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                ierr = nf90_put_att(imapfile, id_taus(iid),  'standard_name', 'taucurrent')
                ierr = nf90_put_att(imapfile, id_taus(iid),  'long_name'    , 'taucurrent in flow element')
                ierr = nf90_put_att(imapfile, id_taus(iid),  'units'        , 'N m-2')
            endif

            if (jamaptidep >0 .and. jatidep >0) then
               ierr = nf90_def_var(imapfile, 'TidalPotential', nf90_double, (/ id_flowelemdim(iid), id_timedim(iid)/), id_tidep(iid))
               ierr = nf90_put_att(imapfile, id_tidep(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
               ierr = nf90_put_att(imapfile, id_tidep(iid),  'standard_name', 'TidalPotential')
               ierr = nf90_put_att(imapfile, id_tidep(iid),  'long_name'    , 'Tidal Potential generated by celestial forces in flow element center')
               ierr = nf90_put_att(imapfile, id_tidep(iid),  'units'        , 'm2 s-2')
            endif
            if (jamapselfal > 0) then
               if ( jaselfal.gt.0 ) then
                  ierr = nf90_def_var(imapfile, 'SALPotential', nf90_double, (/ id_flowelemdim(iid), id_timedim(iid)/), id_salp(iid))
                  ierr = nf90_put_att(imapfile, id_salp(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                  ierr = nf90_put_att(imapfile, id_salp(iid),  'standard_name', 'SALPotential')
                  ierr = nf90_put_att(imapfile, id_salp(iid),  'long_name'    , 'Self-attraction and loading Potential in flow element center')
                  ierr = nf90_put_att(imapfile, id_salp(iid),  'units'        , 'm2 s-2')
               endif
            endif

            if (jaFrcInternalTides2D >0 .and. jamapIntTidesDiss >0) then
               ierr = nf90_def_var(imapfile, 'internal_tides_dissipation', nf90_double, (/ id_flowelemdim(iid), id_timedim(iid)/), id_IntTidesDiss(iid))
               ierr = nf90_put_att(imapfile, id_inttidesdiss(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
               ierr = nf90_put_att(imapfile, id_inttidesdiss(iid),  'standard_name', 'internal_tides_dissipation')
               ierr = nf90_put_att(imapfile, id_inttidesdiss(iid),  'long_name'    , 'internal tides dissipation in flow element center')
               ierr = nf90_put_att(imapfile, id_inttidesdiss(iid),  'units'        , 'J s-1 m-2')
            endif

            if (kmx > 0) then
                !     3D
                if (jamapu1 > 0) then
                    ierr = nf90_def_var(imapfile, 'unorm', nf90_double, (/ id_laydim(iid), id_flowlinkdim(iid), id_timedim (iid)/) , id_unorm(iid))
                endif
                if (jamapu0 > 0) then
                    ierr = nf90_def_var(imapfile, 'u0'   , nf90_double, (/ id_laydim(iid), id_flowlinkdim(iid), id_timedim (iid)/) , id_u0(iid)   )
                endif
                if (jamapq1 > 0) then
                    ierr = nf90_def_var(imapfile, 'q1'   , nf90_double, (/ id_laydim(iid), id_flowlinkdim(iid), id_timedim (iid)/) , id_q1(iid)   )
                endif
                if (jamapq1main > 0 .and. allocated(q1_main)) then
                    ierr = nf90_def_var(imapfile, 'q1main', nf90_double, (/ id_laydim(iid), id_flowlinkdim(iid), id_timedim (iid)/) , id_q1main(iid)   )
                endif
                if (jamapviu > 0) then
                    ierr = nf90_def_var(imapfile, 'viu'   , nf90_double, (/ id_laydim(iid), id_flowlinkdim(iid), id_timedim (iid)/) , id_viu(iid)   )
                endif
                if (jamapdiu > 0) then
                    ierr = nf90_def_var(imapfile, 'diu'   , nf90_double, (/ id_laydim(iid), id_flowlinkdim(iid), id_timedim (iid)/) , id_diu(iid)   )
                endif

                if (jamapucvec > 0) then
                   ! JRE Velocity vector needs to be written, irrespective of kmx, also for com file. Statements moved down outside if-clause
                   !    ierr = nf90_def_var(imapfile, 'ucx'   , nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_ucx(iid)  )
                   !    ierr = nf90_def_var(imapfile, 'ucy'   , nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_ucy(iid)  )
                    ierr = nf90_def_var(imapfile, 'ucz'  , nf90_double, (/ id_laydim(iid), id_flowelemdim(iid), id_timedim (iid)/) , id_ucz(iid)  )

                   ! Depth-averaged cell-center velocities in 3D:
                   ierr = nf90_def_var(imapfile, 'ucxa' , nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_ucxa(iid)  )
                   ierr = nf90_def_var(imapfile, 'ucya' , nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_ucya(iid)  )

                endif
                if (jamapww1 > 0) then
                    ierr = nf90_def_var(imapfile, 'ww1'  , nf90_double, (/ id_wdim(iid), id_flowelemdim(iid), id_timedim (iid)/) , id_ww1(iid))
                endif
                if (jamaprho > 0) then
                    if ( density_is_pressure_dependent() ) then
                        ierr = nf90_def_var(imapfile, 'density', nf90_double, (/ id_laydim(iid), id_flowelemdim(iid), id_timedim (iid)/) , id_rho(iid))
                    else
                        ierr = nf90_def_var(imapfile, 'rho',     nf90_double, (/ id_laydim(iid), id_flowelemdim(iid), id_timedim (iid)/) , id_rhop(iid))
                    endif
                endif
              !
                if (jamapucvec > 0) then
                  ierr = nf90_put_att(imapfile, id_ucz(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                  ierr = nf90_put_att(imapfile, id_ucz(iid),  'standard_name', 'upward_sea_water_velocity')
                  ierr = nf90_put_att(imapfile, id_ucz(iid),  'long_name'    , 'upward velocity on flow element center')
                  ierr = nf90_put_att(imapfile, id_ucz(iid),  'units'        , 'm s-1')
                  ierr = nf90_put_att(imapfile, id_ucz(iid),  '_FillValue'   , dmiss)

                  ierr = nf90_put_att(imapfile, id_ucxa(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                  if (jsferic == 0) then
                     ierr = nf90_put_att(imapfile, id_ucxa(iid),  'standard_name', 'sea_water_x_velocity')
                  else
                     ierr = nf90_put_att(imapfile, id_ucxa(iid),  'standard_name', 'eastward_sea_water_velocity')
                  endif

                  ierr = nf90_put_att(imapfile, id_ucxa(iid),  'long_name'    , 'depth-averaged velocity on flow element center, x-component')
                  ierr = nf90_put_att(imapfile, id_ucxa(iid),  'units'        , 'm s-1')

                  ierr = nf90_put_att(imapfile, id_ucya(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                  if (jsferic == 0) then
                     ierr = nf90_put_att(imapfile, id_ucya(iid),  'standard_name', 'sea_water_y_velocity')
                  else
                     ierr = nf90_put_att(imapfile, id_ucya(iid),  'standard_name', 'northward_sea_water_velocity')
                  endif
                  ierr = nf90_put_att(imapfile, id_ucya(iid),  'long_name'    , 'depth-averaged velocity on flow element center, y-component')
                  ierr = nf90_put_att(imapfile, id_ucya(iid),  'units'        , 'm s-1')
                endif
                if (jamapww1 > 0) then
                  ierr = nf90_put_att(imapfile, id_ww1(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                  ierr = nf90_put_att(imapfile, id_ww1(iid),  'standard_name', 'upward_sea_water_velocity')              ! same standard name allowed?
                  ierr = nf90_put_att(imapfile, id_ww1(iid),  'long_name'    , 'upward velocity on vertical interface')  ! (upward normal or upward)?
                  ierr = nf90_put_att(imapfile, id_ww1(iid),  'units'        , 'm s-1')
                  ierr = nf90_put_att(imapfile, id_ww1(iid),  '_FillValue'   , dmiss)
                  !?elevation
                endif
                if (jamaprho > 0) then
                  if ( density_is_pressure_dependent() ) then
                    ierr = nf90_put_att(imapfile, id_rho(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_rho(iid),  'standard_name', 'sea_water_density')
                    ierr = nf90_put_att(imapfile, id_rho(iid),  'long_name'    , 'flow mass density')
                    ierr = nf90_put_att(imapfile, id_rho(iid),  'units'        , 'kg m-3')
                    ierr = nf90_put_att(imapfile, id_rho(iid),  '_FillValue'   , dmiss)
                  else
                    ierr = nf90_put_att(imapfile, id_rhop(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_rhop(iid),  'standard_name', 'sea_water_potential_density')
                    ierr = nf90_put_att(imapfile, id_rhop(iid),  'long_name'    , 'flow mass potential density')
                    ierr = nf90_put_att(imapfile, id_rhop(iid),  'units'        , 'kg m-3')
                    ierr = nf90_put_att(imapfile, id_rhop(iid),  '_FillValue'   , dmiss)
                  endif
                endif
            endif  ! kmx>0

            if (kmx == 0) then
               if (jamapu1 > 0) then
                  ierr = nf90_def_var(imapfile, 'unorm' , nf90_double, (/ id_flowlinkdim(iid), id_timedim (iid)/) , id_unorm(iid))
               endif
               if (jamapu0 > 0) then
                  ierr = nf90_def_var(imapfile, 'u0'    , nf90_double, (/ id_flowlinkdim(iid), id_timedim (iid)/) , id_u0(iid)   )
               endif
               if (jamapq1 > 0) then
                  ierr = nf90_def_var(imapfile, 'q1'    , nf90_double, (/ id_flowlinkdim(iid), id_timedim (iid)/) , id_q1(iid)   )
               endif
               if (jamapq1main > 0 .and. allocated(q1_main)) then
                  ierr = nf90_def_var(imapfile, 'q1main', nf90_double, (/ id_flowlinkdim(iid), id_timedim (iid)/) , id_q1main(iid)   )
               endif
               if (jamapviu > 0) then
                  ierr = nf90_def_var(imapfile, 'viu'    , nf90_double, (/ id_flowlinkdim(iid), id_timedim (iid)/) , id_viu(iid)   )
               endif
               if (jamapdiu > 0) then
                  ierr = nf90_def_var(imapfile, 'diu'    , nf90_double, (/ id_flowlinkdim(iid), id_timedim (iid)/) , id_diu(iid)   )
               endif
            endif

            if (jamapu1 > 0) then
               ierr = nf90_put_att(imapfile, id_unorm(iid),'coordinates'  , 'FlowLink_xu FlowLink_yu')
               ierr = nf90_put_att(imapfile, id_unorm(iid),'long_name', 'normal component of sea_water_speed')
               ierr = nf90_put_att(imapfile, id_unorm(iid),'units'        , 'm s-1')
               ierr = nf90_put_att(imapfile, id_unorm(iid),'_FillValue'   , dmiss)
            endif

            if (jamapu0 > 0) then
               ierr = nf90_put_att(imapfile, id_u0(iid)   ,'coordinates'  , 'FlowLink_xu FlowLink_yu')
               ierr = nf90_put_att(imapfile, id_u0(iid)   ,'long_name',     'normal component of sea_water_speed at previous timestep')
               ierr = nf90_put_att(imapfile, id_u0(iid)   ,'units'        , 'm s-1')
               ierr = nf90_put_att(imapfile, id_u0(iid)   ,'_FillValue'   , dmiss)
            endif
            if (jamapq1 > 0) then
               ierr = nf90_put_att(imapfile, id_q1(iid)   ,'coordinates'  , 'FlowLink_xu FlowLink_yu')
               !ierr = nf90_put_att(imapfile, id_q1(iid)   ,'standard_name', 'discharge') ! not CF
               ierr = nf90_put_att(imapfile, id_q1(iid)   ,'long_name'    , 'flow flux')
               ierr = nf90_put_att(imapfile, id_q1(iid)   ,'units'        , 'm3 s-1')
               ierr = nf90_put_att(imapfile, id_q1(iid)   ,'_FillValue'   , dmiss)
            endif
            if (jamapq1main > 0 .and. allocated(q1_main)) then
               ierr = nf90_put_att(imapfile, id_q1main(iid)   ,'coordinates'  , 'FlowLink_xu FlowLink_yu')
               !ierr = nf90_put_att(imapfile, id_q1main(iid)   ,'standard_name', 'discharge') ! not CF
               ierr = nf90_put_att(imapfile, id_q1main(iid)   ,'long_name'    , 'flow flux in main channel')
               ierr = nf90_put_att(imapfile, id_q1main(iid)   ,'units'        , 'm3 s-1')
               ierr = nf90_put_att(imapfile, id_q1main(iid)   ,'_FillValue'   , dmiss)
            endif

            if (jamapviu > 0) then
               ierr = nf90_put_att(imapfile, id_viu(iid)   ,'coordinates'  , 'FlowLink_xu FlowLink_yu')
               ierr = nf90_put_att(imapfile, id_viu(iid)   ,'long_name',     'horizontal viscosity')
               ierr = nf90_put_att(imapfile, id_viu(iid)   ,'units'        , 'm2 s-1')
               ierr = nf90_put_att(imapfile, id_viu(iid)   ,'_FillValue'   , dmiss)
            endif
            if (jamapdiu > 0) then
               ierr = nf90_put_att(imapfile, id_diu(iid)   ,'coordinates'  , 'FlowLink_xu FlowLink_yu')
               ierr = nf90_put_att(imapfile, id_diu(iid)   ,'long_name',     'horizontal diffusivity')
               ierr = nf90_put_att(imapfile, id_diu(iid)   ,'units'        , 'm2 s-1')
               ierr = nf90_put_att(imapfile, id_diu(iid)   ,'_FillValue'   , dmiss)
            endif
        endif   ! jaseparate =/ 2
        !
        if (kmx==0) then
           if (jamapucvec > 0 .or. jaseparate_==2) then
               ierr = nf90_def_var(imapfile, 'ucx'   , nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_ucx(iid)  )
               ierr = nf90_def_var(imapfile, 'ucy'   , nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_ucy(iid)  )
           endif
        else
           if (jamapucvec > 0 .or. jaseparate_==2) then
              ierr = nf90_def_var(imapfile, 'ucx'  , nf90_double, (/ id_laydim(iid), id_flowelemdim(iid), id_timedim (iid)/) , id_ucx(iid)  )
              ierr = nf90_def_var(imapfile, 'ucy'  , nf90_double, (/ id_laydim(iid), id_flowelemdim(iid), id_timedim (iid)/) , id_ucy(iid)  )
           endif
        endif

        if (jamapucvec > 0 .or. jaseparate_==2) then
            ierr = nf90_put_att(imapfile, id_ucx(iid)  ,'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
            if (jsferic == 0) then
               ierr = nf90_put_att(imapfile, id_ucx(iid)  ,'standard_name', 'sea_water_x_velocity')
            else
               ierr = nf90_put_att(imapfile, id_ucx(iid)  ,'standard_name', 'eastward_sea_water_velocity')
            endif

            if (jaeulervel==0 .or. jaseparate_==2) then
               ierr = nf90_put_att(imapfile, id_ucx(iid)  ,'long_name'    , 'velocity on flow element center, x-component')
            else
               ierr = nf90_put_att(imapfile, id_ucx(iid)  ,'long_name'    , 'Eulerian velocity on flow element center, x-component')
            endif
            ierr = nf90_put_att(imapfile, id_ucx(iid)  ,'units'        , 'm s-1')
            ierr = nf90_put_att(imapfile, id_ucx(iid)  ,'_FillValue'   , dmiss)

            ierr = nf90_put_att(imapfile, id_ucy(iid)  ,'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
            if (jsferic == 0) then
               ierr = nf90_put_att(imapfile, id_ucy(iid)  ,'standard_name', 'sea_water_y_velocity')
            else
               ierr = nf90_put_att(imapfile, id_ucy(iid)  ,'standard_name', 'northward_sea_water_velocity')
            endif

            if (jaeulervel==0 .or. jaseparate_==2) then
               ierr = nf90_put_att(imapfile, id_ucy(iid)  ,'long_name'    , 'velocity on flow element center, y-component')
            else
               ierr = nf90_put_att(imapfile, id_ucy(iid)  ,'long_name'    , 'Eulerian velocity on flow element center, y-component')
            endif
            ierr = nf90_put_att(imapfile, id_ucy(iid)  ,'units'        , 'm s-1')
            ierr = nf90_put_att(imapfile, id_ucy(iid)  ,'_FillValue'   , dmiss)
        endif

        if (jaseparate_ /= 2) then
           if (jamapsal > 0 .and. jasal > 0) then
              if ( kmx > 0 ) then  !        3D
                 ierr = nf90_def_var(imapfile, 'sa1' , nf90_double, (/ id_laydim(iid), id_flowelemdim (iid), id_timedim (iid)/) , id_sa1(iid))
              else
                 ierr = nf90_def_var(imapfile, 'sa1' , nf90_double, (/ id_flowelemdim (iid), id_timedim (iid)/) , id_sa1(iid))
              endif
              ierr = nf90_put_att(imapfile, id_sa1(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
              ierr = nf90_put_att(imapfile, id_sa1(iid),  'standard_name', 'sea_water_salinity')
              ierr = nf90_put_att(imapfile, id_sa1(iid),  'long_name'    , 'salinity')
              ierr = nf90_put_att(imapfile, id_sa1(iid),  'units'        , '1e-3')
              ierr = nf90_put_att(imapfile, id_sa1(iid),  '_FillValue'   , dmiss)
           endif

           if (jamaptem > 0 .and. jatem > 0) then
              if ( kmx > 0 ) then !        3D
                ierr = nf90_def_var(imapfile, 'tem1' , nf90_double, (/ id_laydim(iid), id_flowelemdim(iid) , id_timedim(iid) /) , id_tem1(iid))
              else
                ierr = nf90_def_var(imapfile, 'tem1' , nf90_double, (/ id_flowelemdim(iid) , id_timedim(iid) /) , id_tem1(iid))
              endif
              ierr = nf90_put_att(imapfile, id_tem1(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
              ierr = nf90_put_att(imapfile, id_tem1(iid),  'standard_name', 'sea_water_temperature')
              ierr = nf90_put_att(imapfile, id_tem1(iid),  'long_name'    , 'temperature')
              ierr = nf90_put_att(imapfile, id_tem1(iid),  'units'        , 'degC')
              ierr = nf90_put_att(imapfile, id_tem1(iid),  '_FillValue'   , dmiss)
           endif

!          tracers
           if (jamapconst > 0 .and. ITRA1 > 0) then
              do j=ITRA1,ITRAN
                 tmpstr = const_names(j)
                 ! Forbidden chars in NetCDF names: space, /, and more.
                 call replace_char(tmpstr,32,95)
                 call replace_char(tmpstr,47,95)
                 if ( kmx > 0 ) then  !        3D
                    ierr = nf90_def_var(imapfile, trim(tmpstr), nf90_double, (/ id_laydim(iid), id_flowelemdim (iid), id_timedim (iid)/) , id_const(iid,j))
                 else
                    ierr = nf90_def_var(imapfile, trim(tmpstr), nf90_double, (/ id_flowelemdim (iid), id_timedim (iid)/) , id_const(iid,j))
                 endif
                 ierr = nf90_put_att(imapfile, id_const(iid,j),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_const(iid,j),  'standard_name', trim(tmpstr))
                 ierr = nf90_put_att(imapfile, id_const(iid,j),  'long_name'    , trim(tmpstr))
                 if (const_units(j).ne.' ') then
                    tmpstr = const_units(j)
                 else
                    tmpstr = '-'
                 endif
                 ierr = nf90_put_att(imapfile, id_const(iid,j),  'units'        , tmpstr)
                 ierr = nf90_put_att(imapfile, id_const(iid,j),  '_FillValue'   , dmiss)
              enddo
           endif

!          water quality bottom variables
           if (numwqbots > 0) then
              call realloc(id_wqb, (/ 3, numwqbots /), keepExisting=.false., fill = 0)
              do j=1,numwqbots
                 tmpstr = wqbotnames(j)
                 ! Forbidden chars in NetCDF names: space, /, and more.
                 call replace_char(tmpstr,32,95)
                 call replace_char(tmpstr,47,95)
                 ierr = nf90_def_var(imapfile, trim(tmpstr), nf90_double, (/ id_flowelemdim (iid), id_timedim (iid)/) , id_wqb(iid,j))
                 ierr = nf90_put_att(imapfile, id_wqb(iid,j),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_wqb(iid,j),  'standard_name', trim(tmpstr))
                 ierr = nf90_put_att(imapfile, id_wqb(iid,j),  'long_name'    , trim(tmpstr))
                 tmpstr = wqbotunits(j)
                 ierr = nf90_put_att(imapfile, id_wqb(iid,j),  'units'        , tmpstr)
                 ierr = nf90_put_att(imapfile, id_wqb(iid,j),  '_FillValue'   , dmiss)
              enddo
              if (wqbot3D_output == 1) then
                 call realloc(id_wqb3d, (/ 3, numwqbots /), keepExisting=.false., fill = 0)
                 do j=1,numwqbots
                    tmpstr = wqbotnames(j)
                    ! Forbidden chars in NetCDF names: space, /, and more.
                    call replace_char(tmpstr,32,95)
                    call replace_char(tmpstr,47,95)
                    ierr = nf90_def_var(imapfile, trim(tmpstr)//'_3D', nf90_double, (/ id_laydim(iid), id_flowelemdim (iid), id_timedim (iid)/) , id_wqb3d(iid,j))
                    ierr = nf90_put_att(imapfile, id_wqb3d(iid,j),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_wqb3d(iid,j),  'standard_name', trim(tmpstr))
                    ierr = nf90_put_att(imapfile, id_wqb3d(iid,j),  'long_name'    , trim(tmpstr))
                    tmpstr = wqbotunits(j)
                    ierr = nf90_put_att(imapfile, id_wqb3d(iid,j),  'units'        , tmpstr)
                    ierr = nf90_put_att(imapfile, id_wqb3d(iid,j),  '_FillValue'   , dmiss)
                 enddo
              endif
           endif

!          waq output
           if (jawaqproc > 0) then
              if (noout_map > 0) then
                 call realloc(id_waq, (/ 3, noout_map /), keepExisting=.false., fill = 0)
                 do j=1,noout_map
                    tmpstr = ' '
                    write (tmpstr, "('water_quality_output_',I0)") j
                    if ( kmx > 0 ) then  !        3D
                       ierr = nf90_def_var(imapfile, tmpstr, nf90_double, (/ id_laydim(iid), id_flowelemdim (iid), id_timedim (iid)/) , id_waq(iid,j))
                    else
                       ierr = nf90_def_var(imapfile, tmpstr, nf90_double, (/ id_flowelemdim (iid), id_timedim (iid)/) , id_waq(iid,j))
                    endif
                    tmpstr = trim(outputs%names(j))//' - '//trim(outputs%description(j))//' in flow element'
                    call replace_multiple_spaces_by_single_spaces(tmpstr)
                    ierr = nf90_put_att(imapfile, id_waq(iid,j),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_waq(iid,j),  'long_name'    , trim(outputs%names(j)))
                    ierr = nf90_put_att(imapfile, id_waq(iid,j),  'units'        , trim(outputs%units(j)))
                    ierr = nf90_put_att(imapfile, id_waq(iid,j),  'description'  , tmpstr)
                    ierr = nf90_put_att(imapfile, id_waq(iid,j),  '_FillValue'   , dmiss)
                 enddo
              endif
              if (noout_statt > 0) then
                 call realloc(id_wqst, (/ 3, noout_statt /), keepExisting=.false., fill = 0)
                 do j=1,noout_statt
                    jj = noout_user + j
                    tmpstr = ' '
                    write (tmpstr, "('water_quality_stat_',I0)") j
                    if ( kmx > 0 ) then  !        3D
                       ierr = nf90_def_var(imapfile, tmpstr, nf90_double, (/ id_laydim(iid), id_flowelemdim (iid), id_timedim (iid)/) , id_wqst(iid,j))
                    else
                       ierr = nf90_def_var(imapfile, tmpstr, nf90_double, (/ id_flowelemdim (iid), id_timedim (iid)/) , id_wqst(iid,j))
                    endif
                    tmpstr = trim(outputs%names(jj))//' - '//trim(outputs%description(jj))//' in flow element'
                    call replace_multiple_spaces_by_single_spaces(tmpstr)
                    ierr = nf90_put_att(imapfile, id_wqst(iid,j),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_wqst(iid,j),  'long_name'    , trim(outputs%names(jj)))
                    ierr = nf90_put_att(imapfile, id_wqst(iid,j),  'units'        , trim(outputs%units(jj)))
                    ierr = nf90_put_att(imapfile, id_wqst(iid,j),  'description'  , tmpstr)
                    ierr = nf90_put_att(imapfile, id_wqst(iid,j),  '_FillValue'   , dmiss)
                 enddo
              endif
              if (noout_state > 0) then
                 call realloc(id_wqse, (/ 3, noout_state /), keepExisting=.false., fill = 0)
                 do j=1,noout_state
                    jj = noout_user + noout_statt + j
                    tmpstr = ' '
                    write (tmpstr, "('water_quality_stat_',I0)") noout_statt + j
                    if ( kmx > 0 ) then  !        3D
                       ierr = nf90_def_var(imapfile, tmpstr, nf90_double, (/ id_laydim(iid), id_flowelemdim (iid)/) , id_wqse(iid,j))
                    else
                       ierr = nf90_def_var(imapfile, tmpstr, nf90_double, (/ id_flowelemdim (iid)/) , id_wqse(iid,j))
                    endif
                    tmpstr = trim(outputs%names(jj))//' - '//trim(outputs%description(jj))//' in flow element'
                    call replace_multiple_spaces_by_single_spaces(tmpstr)
                    ierr = nf90_put_att(imapfile, id_wqse(iid,j),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_wqse(iid,j),  'long_name'    , trim(outputs%names(jj)))
                    ierr = nf90_put_att(imapfile, id_wqse(iid,j),  'units'        , trim(outputs%units(jj)))
                    ierr = nf90_put_att(imapfile, id_wqse(iid,j),  'description'  , tmpstr)
                    ierr = nf90_put_att(imapfile, id_wqse(iid,j),  '_FillValue'   , dmiss)
                 enddo
              endif
           endif

           ! water quality mass balance areas
           if (nomba > 0) then
              ierr = nf90_def_var(imapfile,  'water_quality_mba', nf90_int, (/ id_flowelemdim (iid) /) , id_mba(iid))
              ierr = nf90_put_att(imapfile, id_mba(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
              ierr = nf90_put_att(imapfile, id_mba(iid),  'long_name'    , 'Water quality mass balance areas')
              ierr = unc_add_gridmapping_att(imapfile, (/ id_mba(iid) /), jsferic)
              call realloc(flag_val, nomba, keepExisting = .false., fill = 0)
              flag_mean = ' '
              do j=nomba,1,-1
                 flag_val(j) = j
                 flag_mean = trim(mbaname(j))//' '//flag_mean
              enddo
              ierr = nf90_put_att(imapfile, id_mba(iid), 'flag_values', flag_val)
              ierr = nf90_put_att(imapfile, id_mba(iid), 'flag_meanings', flag_mean)
           endif

           if ( jasecflow > 0 .and. jamapspir > 0) then
              if (kmx < 2) then
                 ierr = nf90_def_var(imapfile, 'spircrv' , nf90_double, (/ id_flowelemdim (iid), id_timedim (iid) /) , id_spircrv(iid))
                 ierr = nf90_put_att(imapfile, id_spircrv(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_spircrv(iid),  'long_name'    , 'streamline curvature')
                 ierr = nf90_put_att(imapfile, id_spircrv(iid),  'units'        , 'm-1')
              endif
              ierr = nf90_def_var(imapfile, 'spirint' , nf90_double, (/ id_flowelemdim (iid), id_timedim (iid) /) , id_spirint(iid))
              ierr = nf90_put_att(imapfile, id_spirint(iid)  ,'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
              ierr = nf90_put_att(imapfile, id_spirint(iid)  ,'long_name'    , 'Spiral flow intensity')
              ierr = nf90_put_att(imapfile, id_spirint(iid)  ,'units'        , 'm/s')
           endif


           if (jamaptur > 0 .and. kmx > 0) then
              if ( iturbulencemodel >= 3 ) then
                 ierr = nf90_def_var(imapfile, 'turkin1' , nf90_double, (/ id_wdim(iid), id_flowlinkdim(iid) , id_timedim(iid) /) , id_turkin1(iid))
                 ierr = nf90_put_att(imapfile, id_turkin1(iid),  'coordinates'  , 'FlowLink_xu FlowLink_yu')
                 ierr = nf90_put_att(imapfile, id_turkin1(iid),  'standard_name', 'specific_turbulent_kinetic_energy_of_sea_water')
                 ierr = nf90_put_att(imapfile, id_turkin1(iid),  'long_name'    , 'turbulent kinetic energy')
                 ierr = nf90_put_att(imapfile, id_turkin1(iid),  'units'        , 'm2 s-2')
                 ierr = nf90_put_att(imapfile, id_turkin1(iid),  '_FillValue'   , dmiss)

                 ierr = nf90_def_var(imapfile, 'vicwwu' , nf90_double, (/ id_wdim(iid), id_flowlinkdim(iid) , id_timedim(iid) /) , id_vicwwu(iid))
                 ierr = nf90_put_att(imapfile, id_vicwwu(iid),  'coordinates'  , 'FlowLink_xu FlowLink_yu')
                 ierr = nf90_put_att(imapfile, id_vicwwu(iid),  'long_name'    , 'turbulent vertical eddy viscosity at velocity points')
                 ierr = nf90_put_att(imapfile, id_vicwwu(iid),  'units'        , 'm2 s-1')
                 ierr = nf90_put_att(imapfile, id_vicwwu(iid),  '_FillValue'   , dmiss)

                 ierr = nf90_def_var(imapfile, 'vicwws' , nf90_double, (/ id_wdim(iid), id_flowelemdim(iid) , id_timedim(iid) /) , id_vicwws(iid))
                 ierr = nf90_put_att(imapfile, id_vicwws(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_vicwws(iid),  'long_name'    , 'turbulent vertical eddy viscosity at pressure points')
                 ierr = nf90_put_att(imapfile, id_vicwws(iid),  'units'        , 'm2 s-1')
                 ierr = nf90_put_att(imapfile, id_vicwws(iid),  '_FillValue'   , dmiss)

                 ierr = nf90_def_var(imapfile, 'tureps1' , nf90_double, (/ id_wdim(iid), id_flowlinkdim(iid) , id_timedim(iid) /) , id_tureps1(iid))
                 ierr = nf90_put_att(imapfile, id_tureps1(iid),  'coordinates'  , 'FlowLink_xu FlowLink_yu')
                 ierr = nf90_put_att(imapfile, id_tureps1(iid),  '_FillValue'   , dmiss)

                 if ( iturbulencemodel == 3 ) then
                    ierr = nf90_put_att(imapfile, id_tureps1(iid),  'standard_name', 'specific_turbulent_kinetic_energy_dissipation_in_sea_water')
                    ierr = nf90_put_att(imapfile, id_tureps1(iid),  'long_name'    , 'turbulent energy dissipation')
                    ierr = nf90_put_att(imapfile, id_tureps1(iid),  'units'        , 'm2 s-3')
                 else if ( iturbulencemodel == 4 ) then
                    ierr = nf90_put_att(imapfile, id_tureps1(iid),  'long_name'    , 'turbulent time scale')
                    ierr = nf90_put_att(imapfile, id_tureps1(iid),  'units'        , 's-1')
                 endif
              endif
           endif

           if (jamapsed > 0 .and. stm_included) then
              ierr = nf90_def_dim(imapfile, 'nSedTot', stmpar%lsedtot, id_sedtotdim(iid))
              ierr = nf90_def_dim(imapfile, 'nSedSus', stmpar%lsedsus, id_sedsusdim(iid))
              ierr = nf90_def_dim(imapfile, 'nBedLayers', stmpar%morlyr%settings%nlyr, id_nlyrdim(iid))
              !
              select case(stmpar%morpar%moroutput%transptype)
                 case (0)
                    transpunit = 'kg/(s m)'
                 case (1)
                    transpunit = 'm3/(s m)'
                 case (2)
                    transpunit = 'm3/(s m)'
              end select
              !
              ! fall velocity
              if (stmpar%lsedsus > 0) then
                 if (kmx > 0) then
                    ierr = nf90_def_var(imapfile, 'ws', nf90_double, (/ id_laydim(iid), id_flowelemdim(iid) , id_sedsusdim(iid) , id_timedim(iid) /), id_ws(iid))
                 else ! '2D' fall velocity, ref fm_erosed(), to check...
                    ierr = nf90_def_var(imapfile, 'ws', nf90_double, (/ id_flowelemdim(iid) , id_sedsusdim(iid) , id_timedim(iid) /), id_ws(iid))
                 endif
                 ierr = nf90_put_att(imapfile, id_ws(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_ws(iid) ,  'long_name'    , 'Sediment settling velocity')
                 ierr = nf90_put_att(imapfile, id_ws(iid) ,  'units'        , 'm s-1')
                 !
                 ! equilibrium concentration, 2D only
                 if (kmx == 0) then
                    ierr = nf90_def_var(imapfile, 'rsedeq', nf90_double, (/ id_flowelemdim(iid) , id_sedsusdim(iid) , id_timedim(iid) /), id_rsedeq(iid))
                    ierr = nf90_put_att(imapfile, id_rsedeq(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_rsedeq(iid) ,  'long_name'    , 'Equilibrium sediment concentration')
                    ierr = nf90_put_att(imapfile, id_rsedeq(iid) ,  'units'        , 'kg m-3')
                 endif
                 !
                 ! reference height
                 if (stmpar%morpar%moroutput%aks) then
                    ierr = nf90_def_var(imapfile, 'aks', nf90_double, (/ id_flowelemdim(iid) , id_sedsusdim(iid) , id_timedim(iid) /), id_aks(iid))
                    ierr = nf90_put_att(imapfile, id_aks(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_aks(iid) ,  'long_name'    , 'Near-bed reference concentration height')
                    ierr = nf90_put_att(imapfile, id_aks(iid) ,  'units'        , 'm')

                    ierr = nf90_def_var(imapfile, 'rca', nf90_double, (/ id_flowelemdim(iid) , id_sedsusdim(iid) , id_timedim(iid) /), id_rca(iid))
                    ierr = nf90_put_att(imapfile, id_rca(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_rca(iid) ,  'long_name'    , 'Near-bed reference concentration')
                    ierr = nf90_put_att(imapfile, id_rca(iid) ,  'units'        , 'kg m-3')
                 endif

                 if (stmpar%morpar%moroutput%sourcesink) then
                    ierr = nf90_def_var(imapfile, 'sourse' , nf90_double, (/ id_flowelemdim(iid) , id_sedsusdim(iid) , id_timedim(iid) /) , id_sourse(iid))
                    ierr = nf90_put_att(imapfile, id_sourse(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_sourse(iid) ,  'long_name'    , 'Source term suspended sediment fractions')
                    ierr = nf90_put_att(imapfile, id_sourse(iid) ,  'units'        , 'kg/(m3 s)')

                    ierr = nf90_def_var(imapfile, 'sinkse' , nf90_double, (/ id_flowelemdim(iid) , id_sedsusdim(iid) , id_timedim(iid) /) , id_sinkse(iid))
                    ierr = nf90_put_att(imapfile, id_sinkse(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_sinkse(iid) ,  'long_name'    , 'Sink term suspended sediment fractions')
                    ierr = nf90_put_att(imapfile, id_sinkse(iid) ,  'units'        , 's-1')
                 endif

                 if (stmpar%morpar%moroutput%suvcor) then
                    ierr = nf90_def_var(imapfile, 'e_scrn' , nf90_double, (/ id_flowlinkdim(iid) , id_sedsusdim(iid) , id_timedim(iid) /) , id_scrn(iid))
                    ierr = nf90_put_att(imapfile, id_scrn(iid) ,  'coordinates'  , 'FlowLink_xu FlowLink_yu')
                    ierr = nf90_put_att(imapfile, id_scrn(iid) ,  'long_name'    , 'Near-bed transport correction in face-normal direction')
                    ierr = nf90_put_att(imapfile, id_scrn(iid) ,  'units'        , transpunit)

                    !ierr = nf90_def_var(imapfile, 'e_scrt' , nf90_double, (/ id_flowlinkdim(iid) , id_sedsusdim(iid) , id_timedim(iid) /) , id_scrt(iid))
                    !ierr = nf90_put_att(imapfile, id_scrt(iid) ,  'coordinates'  , 'FlowLink_xu FlowLink_yu')
                    !ierr = nf90_put_att(imapfile, id_scrt(iid) ,  'long_name'    , 'Near-bed transport correction face-tangential direction')
                    !ierr = nf90_put_att(imapfile, id_scrt(iid) ,  'units'        , transpunit)
                 endif
                 !
                 ! Suspended fractions
                 !
                 do j=ISED1,ISEDN
                    tmpstr = const_names(j)
                    ! Forbidden chars in NetCDF names: space, /, and more.
                    call replace_char(tmpstr,32,95)
                    call replace_char(tmpstr,47,95)
                    if ( kmx > 0 ) then  !        3D
                       ierr = nf90_def_var(imapfile, trim(tmpstr), nf90_double, (/ id_laydim(iid), id_flowelemdim (iid), id_timedim (iid)/) , id_const(iid,j))
                    else
                       ierr = nf90_def_var(imapfile, trim(tmpstr), nf90_double, (/ id_flowelemdim (iid), id_timedim (iid)/) , id_const(iid,j))
                    endif
                    ierr = nf90_put_att(imapfile, id_const(iid,j),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_const(iid,j),  'standard_name', trim(tmpstr)//' concentration')
                    ierr = nf90_put_att(imapfile, id_const(iid,j),  'long_name'    , trim(tmpstr)//' concentration')
                    ierr = nf90_put_att(imapfile, id_const(iid,j),  'units'        , 'kg m-3')
                 enddo
              endif

              if (stmpar%morpar%moroutput%dzduuvv) then ! bedslope
                 ierr = nf90_def_var(imapfile, 'e_dzdn', nf90_double, (/ id_flowlinkdim(iid) , id_timedim(iid) /), id_dzdn(iid))
                 ierr = nf90_put_att(imapfile, id_dzdn(iid) ,  'coordinates'  , 'FlowLink_xu FlowLink_yu')
                 ierr = nf90_put_att(imapfile, id_dzdn(iid) ,  'long_name'    , 'Bed slope, n-component')
                 ierr = nf90_put_att(imapfile, id_dzdn(iid) ,  'units'        , '-')

                 ierr = nf90_def_var(imapfile, 'e_dzdt', nf90_double, (/ id_flowlinkdim(iid) , id_timedim(iid) /), id_dzdt(iid))
                 ierr = nf90_put_att(imapfile, id_dzdt(iid) ,  'coordinates'  , 'FlowLink_xu FlowLink_yu')
                 ierr = nf90_put_att(imapfile, id_dzdt(iid) ,  'long_name'    , 'Bed slope, t-component')
                 ierr = nf90_put_att(imapfile, id_dzdt(iid) ,  'units'        , '-')
              endif

              if (stmpar%morpar%moroutput%umod) then
                 ierr = nf90_def_var(imapfile, 'umod', nf90_double, (/ id_flowelemdim(iid) , id_timedim(iid) /), id_umod(iid))
                 ierr = nf90_put_att(imapfile, id_umod(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_umod(iid) ,  'long_name'    , 'Characteristic velocity magnitude in cell centre')
                 ierr = nf90_put_att(imapfile, id_umod(iid) ,  'units'        , 'm s-1')
              endif

              if (stmpar%morpar%moroutput%zumod) then
                 ierr = nf90_def_var(imapfile, 'zumod', nf90_double, (/ id_flowelemdim(iid) , id_timedim(iid) /), id_zumod(iid))
                 ierr = nf90_put_att(imapfile, id_zumod(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_zumod(iid) ,  'long_name'    , 'Height above bed for characteristic velocity in cell centre')
                 ierr = nf90_put_att(imapfile, id_zumod(iid) ,  'units'        , 'm')
              endif

              if (stmpar%morpar%moroutput%ustar) then
                 ierr = nf90_def_var(imapfile, 'ustar', nf90_double, (/ id_flowelemdim(iid) , id_timedim(iid) /), id_ustar(iid))
                 ierr = nf90_put_att(imapfile, id_ustar(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_ustar(iid) ,  'long_name'    , 'Bed shear velocity u* in cell centre')
                 ierr = nf90_put_att(imapfile, id_ustar(iid) ,  'units'        , 'm s-1')
              endif

              if (stmpar%morpar%moroutput%sbcuv) then
                 ierr = nf90_def_var(imapfile, 'sbcx' , nf90_double, (/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /) , id_sbcx(iid))
                 ierr = nf90_put_att(imapfile, id_sbcx(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sbcx(iid) ,  'long_name'    , 'bed load transport due to currents, x-component')
                 ierr = nf90_put_att(imapfile, id_sbcx(iid) ,  'units'        , transpunit)

                 ierr = nf90_def_var(imapfile, 'sbcy' , nf90_double, (/ id_flowelemdim (iid), id_sedtotdim(iid), id_timedim (iid)/) , id_sbcy(iid))
                 ierr = nf90_put_att(imapfile, id_sbcy (iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sbcy (iid),  'long_name'    , 'bed load transport due to currents, y-component')
                 ierr = nf90_put_att(imapfile, id_sbcy (iid),  'units'        , transpunit)

                 ierr = nf90_def_var(imapfile, 'sbcx_reconstructed' , nf90_double, (/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /) , id_sbcx_reconstructed(iid))
                 ierr = nf90_put_att(imapfile, id_sbcx(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sbcx(iid) ,  'long_name'    , 'bed load transport due to currents (reconstructed), x-component')
                 ierr = nf90_put_att(imapfile, id_sbcx(iid) ,  'units'        , transpunit)

                 ierr = nf90_def_var(imapfile, 'sbcy_reconstructed' , nf90_double, (/ id_flowelemdim (iid), id_sedtotdim(iid), id_timedim (iid)/) , id_sbcy_reconstructed(iid))
                 ierr = nf90_put_att(imapfile, id_sbcy (iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sbcy (iid),  'long_name'    , 'bed load transport due to currents (reconstructed), y-component')
                 ierr = nf90_put_att(imapfile, id_sbcy (iid),  'units'        , transpunit)
              endif

              if (stmpar%morpar%moroutput%sbwuv) then
                 ierr = nf90_def_var(imapfile, 'sbwx' , nf90_double, (/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /) , id_sbwx(iid))
                 ierr = nf90_put_att(imapfile, id_sbwx(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sbwx(iid) ,  'long_name'    , 'bed load transport due to waves, x-component')
                 ierr = nf90_put_att(imapfile, id_sbwx(iid) ,  'units'        , transpunit)

                 ierr = nf90_def_var(imapfile, 'sbwy' , nf90_double, (/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /) , id_sbwy(iid))
                 ierr = nf90_put_att(imapfile, id_sbwy(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sbwy(iid) ,  'long_name'    , 'bed load transport due to waves, y-component')
                 ierr = nf90_put_att(imapfile, id_sbwy(iid) ,  'units'        , transpunit)

                 ierr = nf90_def_var(imapfile, 'sbwx_reconstructed' , nf90_double, (/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /) , id_sbwx_reconstructed(iid))
                 ierr = nf90_put_att(imapfile, id_sbwx(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sbwx(iid) ,  'long_name'    , 'bed load transport due to waves (reconstructed), x-component')
                 ierr = nf90_put_att(imapfile, id_sbwx(iid) ,  'units'        , transpunit)

                 ierr = nf90_def_var(imapfile, 'sbwy_reconstructed' , nf90_double, (/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /) , id_sbwy_reconstructed(iid))
                 ierr = nf90_put_att(imapfile, id_sbwy(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sbwy(iid) ,  'long_name'    , 'bed load transport due to waves (reconstructed), y-component')
                 ierr = nf90_put_att(imapfile, id_sbwy(iid) ,  'units'        , transpunit)
              endif

              if (stmpar%morpar%moroutput%sswuv) then
                 ierr = nf90_def_var(imapfile, 'sswx' , nf90_double, (/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /) , id_sswx(iid))
                 ierr = nf90_put_att(imapfile, id_sswx(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sswx(iid) ,  'long_name'    , 'suspended load transport due to waves, x-component')
                 ierr = nf90_put_att(imapfile, id_sswx(iid) ,  'units'        , transpunit)

                 ierr = nf90_def_var(imapfile, 'sswy' , nf90_double, (/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /) , id_sswy(iid))
                 ierr = nf90_put_att(imapfile, id_sswy(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sswy(iid) ,  'long_name'    , 'suspended load transport due to waves, y-component')
                 ierr = nf90_put_att(imapfile, id_sswy(iid) ,  'units'        , transpunit)

                 ierr = nf90_def_var(imapfile, 'sswx_reconstructed' , nf90_double, (/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /) , id_sswx_reconstructed(iid))
                 ierr = nf90_put_att(imapfile, id_sswx(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sswx(iid) ,  'long_name'    , 'suspended load transport due to waves (reconstructed), x-component')
                 ierr = nf90_put_att(imapfile, id_sswx(iid) ,  'units'        , transpunit)

                 ierr = nf90_def_var(imapfile, 'sswy_reconstructed' , nf90_double, (/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /) , id_sswy_reconstructed(iid))
                 ierr = nf90_put_att(imapfile, id_sswy(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sswy(iid) ,  'long_name'    , 'suspended load transport due to waves (reconstructed), y-component')
                 ierr = nf90_put_att(imapfile, id_sswy(iid) ,  'units'        , transpunit)
              endif

              if (stmpar%morpar%moroutput%sscuv) then
                 ierr = nf90_def_var(imapfile, 'sscx' , nf90_double, (/ id_flowelemdim(iid) , id_sedsusdim(iid) , id_timedim(iid) /) , id_sscx(iid))
                 ierr = nf90_put_att(imapfile, id_sscx(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sscx(iid) ,  'long_name'    , 'suspended load transport due to currents, x-component')
                 ierr = nf90_put_att(imapfile, id_sscx(iid) ,  'units'        , transpunit)

                 ierr = nf90_def_var(imapfile, 'sscy' , nf90_double, (/ id_flowelemdim(iid) , id_sedsusdim(iid) , id_timedim(iid) /) , id_sscy(iid))
                 ierr = nf90_put_att(imapfile, id_sscy(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sscy(iid) ,  'long_name'    , 'suspended load transport due to currents, y-component')
                 ierr = nf90_put_att(imapfile, id_sscy(iid) ,  'units'        , transpunit)

                 ierr = nf90_def_var(imapfile, 'sscx_reconstructed' , nf90_double, (/ id_flowelemdim(iid) , id_sedsusdim(iid) , id_timedim(iid) /) , id_sscx_reconstructed(iid))
                 ierr = nf90_put_att(imapfile, id_sscx(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sscx(iid) ,  'long_name'    , 'suspended load transport due to currents (reconstructed), x-component')
                 ierr = nf90_put_att(imapfile, id_sscx(iid) ,  'units'        , transpunit)

                 ierr = nf90_def_var(imapfile, 'sscy_reconstructed' , nf90_double, (/ id_flowelemdim(iid) , id_sedsusdim(iid) , id_timedim(iid) /) , id_sscy_reconstructed(iid))
                 ierr = nf90_put_att(imapfile, id_sscy(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sscy(iid) ,  'long_name'    , 'suspended load transport due to currents (reconstructed), y-component')
                 ierr = nf90_put_att(imapfile, id_sscy(iid) ,  'units'        , transpunit)
              endif

              ierr = nf90_def_var(imapfile, 'sxtot' , nf90_double, (/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /) , id_sxtot(iid))
              ierr = nf90_put_att(imapfile, id_sxtot(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
              ierr = nf90_put_att(imapfile, id_sxtot(iid) ,  'long_name'    , 'total sediment transport in flow cell center, x-component')
              ierr = nf90_put_att(imapfile, id_sxtot(iid) ,  'units'        , transpunit)

              ierr = nf90_def_var(imapfile, 'sytot' , nf90_double, (/ id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /) , id_sytot(iid))
              ierr = nf90_put_att(imapfile, id_sytot(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
              ierr = nf90_put_att(imapfile, id_sytot(iid) ,  'long_name'    , 'total sediment transport in flow cell center, y-component')
              ierr = nf90_put_att(imapfile, id_sytot(iid) ,  'units'        , transpunit)

              ierr = nf90_def_var(imapfile, 'mor_bl' , nf90_double, (/ id_flowelemdim(iid) , id_timedim(iid) /) , id_morbl(iid))
              ierr = nf90_put_att(imapfile, id_morbl(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
              ierr = nf90_put_att(imapfile, id_morbl(iid) ,  'long_name'    , 'Time-varying bottom level in flow cell center')
              ierr = nf90_put_att(imapfile, id_morbl(iid) ,  'units'        , 'm')


              select case (stmpar%morlyr%settings%iunderlyr)
              case (1)
                 ierr = nf90_def_var(imapfile, 'bodsed' , nf90_double, (/ id_sedtotdim(iid) , id_flowelemdim(iid) , id_timedim(iid) /) , id_bodsed(iid))
                 ierr = nf90_put_att(imapfile, id_bodsed(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_bodsed(iid) ,  'long_name'    , 'available sediment in the bed in flow cell center')
                 ierr = nf90_put_att(imapfile, id_bodsed(iid) ,  'units'        , 'kg m-2')

                 ierr = nf90_def_var(imapfile, 'dpsed' , nf90_double, (/ id_flowelemdim(iid) , id_timedim(iid) /) , id_dpsed(iid))
                 ierr = nf90_put_att(imapfile, id_dpsed(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_dpsed(iid) ,  'long_name'    , 'sediment thickness in the bed in flow cell center')
                 ierr = nf90_put_att(imapfile, id_dpsed(iid) ,  'units'        , 'm')
              case (2)
                 ierr = nf90_def_var(imapfile, 'msed' , nf90_double, (/ id_sedtotdim(iid) , id_nlyrdim(iid) , id_flowelemdim(iid) , id_timedim(iid) /) , id_msed(iid))
                 ierr = nf90_put_att(imapfile, id_msed(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_msed(iid) ,  'long_name'    , 'available sediment in a layer of the bed in flow cell center')
                 ierr = nf90_put_att(imapfile, id_msed(iid) ,  'units'        , 'kg m-2')

                 ierr = nf90_def_var(imapfile, 'lyrfrac' , nf90_double, (/ id_sedtotdim(iid) , id_nlyrdim(iid) , id_flowelemdim(iid) , id_timedim(iid) /) , id_lyrfrac(iid))
                 ierr = nf90_put_att(imapfile, id_lyrfrac(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_lyrfrac(iid) ,  'long_name'    , 'volume fraction in a layer of the bed in flow cell center')
                 ierr = nf90_put_att(imapfile, id_lyrfrac(iid) ,  'units'        , '-')

                 ierr = nf90_def_var(imapfile, 'thlyr' , nf90_double, (/ id_nlyrdim(iid) , id_flowelemdim(iid) , id_timedim(iid) /) , id_thlyr(iid))
                 ierr = nf90_put_att(imapfile, id_thlyr(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_thlyr(iid) ,  'long_name'    , 'thickness of a layer of the bed in flow cell center')
                 ierr = nf90_put_att(imapfile, id_thlyr(iid) ,  'units'        , 'm')

                 if (stmpar%morlyr%settings%iporosity>0) then
                    ierr = nf90_def_var(imapfile, 'poros' , nf90_double, (/ id_nlyrdim(iid) , id_flowelemdim(iid) , id_timedim(iid) /) , id_poros(iid))
                    ierr = nf90_put_att(imapfile, id_poros(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_poros(iid) ,  'long_name'    , 'porosity of a layer of the bed in flow cell center')
                    ierr = nf90_put_att(imapfile, id_poros(iid) ,  'units'        , '-')
                 endif
              end select

              if (stmpar%morpar%moroutput%taurat) then
                 ierr = nf90_def_var(imapfile, 'taurat' , nf90_double, (/id_flowelemdim(iid) , id_sedtotdim(iid) ,id_timedim(iid) /) , id_taurat(iid))
                 ierr = nf90_put_att(imapfile, id_taurat(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_taurat(iid) ,  'long_name'    , 'Excess bed shear ratio')
                 ierr = nf90_put_att(imapfile, id_taurat(iid) ,  'units'        , '-')
              endif
              if (stmpar%morpar%moroutput%dm) then
                 ierr = nf90_def_var(imapfile, 'dm' , nf90_double, (/id_flowelemdim(iid) , id_timedim(iid) /) , id_dm(iid))
                 ierr = nf90_put_att(imapfile, id_dm(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_dm(iid) ,  'long_name'    , 'Arithmetic mean sediment diameter')
                 ierr = nf90_put_att(imapfile, id_dm(iid) ,  'units'        , 'm')
              endif
              if (stmpar%morpar%moroutput%dg) then
                 ierr = nf90_def_var(imapfile, 'dg' , nf90_double, (/id_flowelemdim(iid) , id_timedim(iid) /) , id_dg(iid))
                 ierr = nf90_put_att(imapfile, id_dg(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_dg(iid) ,  'long_name'    , 'Geometric mean sediment diameter')
                 ierr = nf90_put_att(imapfile, id_dg(iid) ,  'units'        , 'm')
              endif
              if (stmpar%morpar%moroutput%dgsd) then
                 ierr = nf90_def_var(imapfile, 'dgsd' , nf90_double, (/id_flowelemdim(iid) , id_timedim(iid) /) , id_dgsd(iid))
                 ierr = nf90_put_att(imapfile, id_dgsd(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_dgsd(iid) ,  'long_name'    , 'Geometric standard deviation of particle size mix')
                 ierr = nf90_put_att(imapfile, id_dgsd(iid) ,  'units'        , 'm')
              endif
              if (stmpar%morpar%moroutput%percentiles) then
                 do l = 1, stmpar%morpar%nxx
                    write(dxname,'(A,I2.2)') 'DXX',l
                    write(dxdescr,'(A,F4.1,A)') 'Sediment diameter percentile '    , stmpar%morpar%xx(l)*100d0,' %'
                    ierr = nf90_def_var(imapfile, dxname , nf90_double, (/id_flowelemdim(iid) , id_timedim(iid) /) , id_dxx(l,iid))
                    ierr = nf90_put_att(imapfile, id_dxx(l,iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                    ierr = nf90_put_att(imapfile, id_dxx(l,iid) ,  'long_name'    , dxdescr)
                    ierr = nf90_put_att(imapfile, id_dxx(l,iid) ,  'units'        , 'm')
                 enddo
              endif
              if (stmpar%morpar%moroutput%frac) then
                 ierr = nf90_def_var(imapfile, 'frac' , nf90_double, (/id_flowelemdim(iid) , id_sedtotdim(iid) , id_timedim(iid) /) , id_frac(iid))
                 ierr = nf90_put_att(imapfile, id_frac(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_frac(iid) ,  'long_name'    , 'Availability fraction in top layer')
                 ierr = nf90_put_att(imapfile, id_frac(iid) ,  'units'        , '-')
              endif
              if (stmpar%morpar%moroutput%mudfrac) then
                 ierr = nf90_def_var(imapfile, 'mudfrac' , nf90_double, (/id_flowelemdim(iid) , id_timedim(iid) /) , id_mudfrac(iid))
                 ierr = nf90_put_att(imapfile, id_mudfrac(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_mudfrac(iid) ,  'long_name'    , 'Mud fraction in top layer')
                 ierr = nf90_put_att(imapfile, id_mudfrac(iid) ,  'units'        , '-')
              endif
              if (stmpar%morpar%moroutput%sandfrac) then
                 ierr = nf90_def_var(imapfile, 'sandfrac' , nf90_double, (/id_flowelemdim(iid) , id_timedim(iid) /) , id_sandfrac(iid))
                 ierr = nf90_put_att(imapfile, id_sandfrac(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_sandfrac(iid) ,  'long_name'    , 'Sand fraction in top layer')
                 ierr = nf90_put_att(imapfile, id_sandfrac(iid) ,  'units'        , '-')
              endif
              if (stmpar%morpar%moroutput%fixfac) then
                 ierr = nf90_def_var(imapfile, 'fixfac' , nf90_double, (/id_flowelemdim(iid) , id_sedtotdim(iid), id_timedim(iid) /) , id_fixfac(iid))
                 ierr = nf90_put_att(imapfile, id_fixfac(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_fixfac(iid) ,  'long_name'    , 'Reduction factor due to limited sediment thickness')
                 ierr = nf90_put_att(imapfile, id_fixfac(iid) ,  'units'        , '-')
              endif
              if (stmpar%morpar%moroutput%hidexp) then
                 ierr = nf90_def_var(imapfile, 'hidexp' , nf90_double, (/id_flowelemdim(iid) , id_sedtotdim(iid), id_timedim(iid) /) , id_hidexp(iid))
                 ierr = nf90_put_att(imapfile, id_hidexp(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_hidexp(iid) ,  'long_name'    , 'Hiding and exposure factor')
                 ierr = nf90_put_att(imapfile, id_hidexp(iid) ,  'units'        , '-')
              endif
              ! Fluff layers
              if (stmpar%morpar%flufflyr%iflufflyr>0 .and. stmpar%lsedsus>0) then
                 ierr = nf90_def_var(imapfile, 'mfluff' , nf90_double, (/id_flowelemdim(iid) , id_sedsusdim(iid), id_timedim(iid) /) , id_mfluff(iid))
                 ierr = nf90_put_att(imapfile, id_mfluff(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_mfluff(iid) ,  'long_name'    , 'Sediment mass in fluff layer')
                 ierr = nf90_put_att(imapfile, id_mfluff(iid) ,  'units'        , 'kg m-2 ')
              endif
           endif

           if (bfmpar%lfbedfrmout) then
              if (bfmpar%lfbedfrm) then
                 ! DUNEHEIGHT
                 ierr = nf90_def_var(imapfile, 'duneheight' , nf90_double, (/ id_flowelemdim(iid) , id_timedim(iid) /) , id_duneheight(iid))
                 ierr = nf90_put_att(imapfile, id_duneheight(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_duneheight(iid) ,  'long_name'    , 'Time-varying dune height in flow cell centers')
                 ierr = nf90_put_att(imapfile, id_duneheight(iid) ,  'units'        , 'm')
                 ! DUNELENGTH
                 ierr = nf90_def_var(imapfile, 'dunelength' , nf90_double, (/ id_flowelemdim(iid) , id_timedim(iid) /) , id_dunelength(iid))
                 ierr = nf90_put_att(imapfile, id_dunelength(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_dunelength(iid) ,  'long_name'    , 'Time-varying dune length in flow cell centers')
                 ierr = nf90_put_att(imapfile, id_dunelength(iid) ,  'units'        , 'm')
              endif
              if (bfmpar%lfbedfrmrou) then
                 call realloc(rks,ndx, keepExisting=.false.,fill=0d0)
                 ! KSR
                 ierr = nf90_def_var(imapfile, 'ksr' , nf90_double, (/ id_flowelemdim(iid) , id_timedim(iid) /) , id_ksr(iid))
                 ierr = nf90_put_att(imapfile, id_ksr(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_ksr(iid) ,  'long_name'    , 'Ripple roughness height in flow cell center')
                 ierr = nf90_put_att(imapfile, id_ksr(iid) ,  'units'        , 'm')
                 ! KSMR
                 ierr = nf90_def_var(imapfile, 'ksmr' , nf90_double, (/ id_flowelemdim(iid) , id_timedim(iid) /) , id_ksmr(iid))
                 ierr = nf90_put_att(imapfile, id_ksmr(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_ksmr(iid) ,  'long_name'    , 'Mega-ripple roughness height in flow cell center')
                 ierr = nf90_put_att(imapfile, id_ksmr(iid) ,  'units'        , 'm')
                 ! KSD
                 ierr = nf90_def_var(imapfile, 'ksd' , nf90_double, (/ id_flowelemdim(iid) , id_timedim(iid) /) , id_ksd(iid))
                 ierr = nf90_put_att(imapfile, id_ksd(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_ksd(iid) ,  'long_name'    , 'Dune roughness height in flow cell center')
                 ierr = nf90_put_att(imapfile, id_ksd(iid) ,  'units'        , 'm')
                 ! KS
                 ierr = nf90_def_var(imapfile, 'ks' , nf90_double, (/ id_flowelemdim(iid) , id_timedim(iid) /) , id_ks(iid))
                 ierr = nf90_put_att(imapfile, id_ks(iid) ,  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                 ierr = nf90_put_att(imapfile, id_ks(iid) ,  'long_name'    , 'Bedform roughness height in flow cell center')
                 ierr = nf90_put_att(imapfile, id_ks(iid) ,  'units'        , 'm')
              endif
           endif
           if (jased > 0 .and. .not.stm_included) then
              ierr = nf90_def_dim(imapfile, 'nFrac', mxgr, id_maxfracdim(iid))

              if (jaceneqtr == 1) then
                  ierr = nf90_inq_dimid(imapfile, 'nFlowElem', id_erolaydim(iid)) ! Note: points to an existing dimension (either nNetNode, or nFlowElem)
                  if (ierr /= nf90_noerr) then
                     ierr = nf90_inq_dimid(imapfile, 'nFlowElemWithBnd', id_erolaydim(iid))
                  endif
              else
                  ierr = nf90_inq_dimid(imapfile, 'nNetNode' , id_erolaydim(iid)) ! Note: points to an existing dimension (either nNetNode, or nFlowElem)
              endif

              ierr = nf90_def_var(imapfile, 'sed'  , nf90_double, (/ id_maxfracdim  (iid), id_flowelemdim(iid), id_timedim (iid)/) , id_sed(iid))
              ierr = nf90_put_att(imapfile, id_sed(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
              ierr = nf90_put_att(imapfile, id_sed(iid),  'long_name'    , 'sediment concentration')
              ierr = nf90_put_att(imapfile, id_sed(iid),  'units'        , 'kg m-3')
              ierr = nf90_def_var(imapfile, 'ero' , nf90_double, (/ id_maxfracdim  (iid), id_erolaydim(iid), id_timedim (iid)/) , id_ero(iid))
              if (jaceneqtr == 1) then
                  ierr = nf90_put_att(imapfile, id_ero(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
                  ierr = nf90_put_att(imapfile, id_ero(iid),  'long_name', 'erodable layer thickness per size fraction in flow element center')
              else
                  ierr = nf90_put_att(imapfile, id_ero(iid),  'coordinates'  , 'NetNode_x NetNode_y')
                  ierr = nf90_put_att(imapfile, id_ero(iid),  'long_name', 'erodable layer thickness per size fraction at flow element corners')
              endif
              ierr = nf90_put_att(imapfile, id_ero(iid),  'standard_name'    , 'Erodable layer thickness') ! Not CF
              ierr = nf90_put_att(imapfile, id_ero(iid),  'units'        , 'm')

              if (jaceneqtr .ne. 1) then
                 idims(1) = id_erolaydim(iid)
                 call definencvar(imapfile,id_zk(iid)   ,nf90_double,idims, 'netnode_bedlevel_zk'  , 'Flow element corner bedlevel (zk)', 'm', 'NetNode_x NetNode_y')
              endif
              idims(1) = id_flowelemdim(iid)
              call definencvar(imapfile,id_bl(iid)   ,nf90_double,idims, 'flowelem_bedlevel_bl'  , 'Flow element center bedlevel (bl)', 'm', 'FlowElem_xcc FlowElem_ycc')

           endif

           ! JRE waves
           if (jawave .eq. 4) then
             ierr = nf90_def_var(imapfile, 'E',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_E(iid))
             ierr = nf90_put_att(imapfile, id_E(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             ierr = nf90_put_att(imapfile, id_E(iid),   'standard_name', 'sea_surface_bulk_wave_energy')                          ! not CF
             ierr = nf90_put_att(imapfile, id_E(iid),   'long_name'    , 'wave energy per square meter')
             ierr = nf90_put_att(imapfile, id_E(iid),   'units'        , 'J m-2')

             ierr = nf90_def_var(imapfile, 'R',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_R(iid))
             ierr = nf90_put_att(imapfile, id_R(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             ierr = nf90_put_att(imapfile, id_R(iid),   'standard_name', 'sea_surface_bulk_roller_energy')                          ! not CF
             ierr = nf90_put_att(imapfile, id_R(iid),   'long_name'    , 'roller energy per square meter')
             ierr = nf90_put_att(imapfile, id_R(iid),   'units'        , 'J m-2')

             ierr = nf90_def_var(imapfile, 'DR',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_DR(iid))
             ierr = nf90_put_att(imapfile, id_DR(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             ierr = nf90_put_att(imapfile, id_DR(iid),   'standard_name', 'sea_surface_bulk_roller_dissipation')                          ! not CF
             ierr = nf90_put_att(imapfile, id_DR(iid),   'long_name'    , 'roller energy dissipation per square meter')
             ierr = nf90_put_att(imapfile, id_DR(iid),   'units'        , 'W m-2')

             ierr = nf90_def_var(imapfile, 'D',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_D(iid))
             ierr = nf90_put_att(imapfile, id_D(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             ierr = nf90_put_att(imapfile, id_D(iid),   'standard_name', 'sea_surface_wave_breaking_dissipation')                          ! not CF
             ierr = nf90_put_att(imapfile, id_D(iid),   'long_name'    , 'wave breaking energy dissipation per square meter')
             ierr = nf90_put_att(imapfile, id_D(iid),   'units'        , 'W m-2')
             ! JRE TO DO: change definition in unc file to correct one
             ierr = nf90_def_var(imapfile, 'H',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_H(iid))
             ierr = nf90_put_att(imapfile, id_H(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             ierr = nf90_put_att(imapfile, id_H(iid),   'standard_name', 'sea_surface_wave_rms_height')
             ierr = nf90_put_att(imapfile, id_H(iid),   'long_name'    , 'Root mean square wave height based on wave energy')
             ierr = nf90_put_att(imapfile, id_H(iid),   'units'        , 'm')

             ierr = nf90_def_var(imapfile, 'urms_cc',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_urmscc(iid))
             ierr = nf90_put_att(imapfile, id_urmscc(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             ierr = nf90_put_att(imapfile, id_urmscc(iid),   'standard_name', 'sea_surface_wave_orbital_velocity')
             ierr = nf90_put_att(imapfile, id_urmscc(iid),   'long_name'    , 'Root mean square orbital velocity on flow centers')
             ierr = nf90_put_att(imapfile, id_urmscc(iid),   'units'        , 'm/s')

             ierr = nf90_def_var(imapfile, 'Fx_cc',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_Fxcc(iid))
             ierr = nf90_put_att(imapfile, id_Fxcc(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             ierr = nf90_put_att(imapfile, id_Fxcc(iid),   'standard_name', 'sea_surface_wave_force_east')
             ierr = nf90_put_att(imapfile, id_Fxcc(iid),   'long_name'    , 'Wave induced flow forcing in cell centre, east component')
             ierr = nf90_put_att(imapfile, id_Fxcc(iid),   'units'        , 'kg m s-2')

             ierr = nf90_def_var(imapfile, 'Fy_cc',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_Fycc(iid))
             ierr = nf90_put_att(imapfile, id_Fycc(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             ierr = nf90_put_att(imapfile, id_Fycc(iid),   'standard_name', 'sea_surface_wave_force_north')
             ierr = nf90_put_att(imapfile, id_Fycc(iid),   'long_name'    , 'Wave induced flow forcing in cell centre, north component')
             ierr = nf90_put_att(imapfile, id_Fycc(iid),   'units'        , 'kg m s-2')

             ierr = nf90_def_var(imapfile, 'thetamean',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid)/) , id_thetamean(iid))
             ierr = nf90_put_att(imapfile, id_thetamean(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             ierr = nf90_put_att(imapfile, id_thetamean(iid),   'standard_name', 'sea_surface_wave_from_direction')                          ! not CF
             ierr = nf90_put_att(imapfile, id_thetamean(iid),   'long_name'    , 'mean wave angle')
             ierr = nf90_put_att(imapfile, id_thetamean(iid),   'units'        , 'deg')

             ierr = nf90_def_var(imapfile, 'cwav',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_cwav(iid))
             ierr = nf90_put_att(imapfile, id_cwav(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             ierr = nf90_put_att(imapfile, id_cwav(iid),   'standard_name', 'sea_surface_wave_phase_celerity')                          ! not CF
             ierr = nf90_put_att(imapfile, id_cwav(iid),   'long_name'    , 'phase celerity')
             ierr = nf90_put_att(imapfile, id_cwav(iid),   'units'        , 'm s-1')

             ierr = nf90_def_var(imapfile, 'cgwav',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_cgwav(iid))
             ierr = nf90_put_att(imapfile, id_cgwav(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             ierr = nf90_put_att(imapfile, id_cgwav(iid),   'standard_name', 'sea_surface_wave_group_celerity')                          ! not CF
             ierr = nf90_put_att(imapfile, id_cgwav(iid),   'long_name'    , 'group celerity')
             ierr = nf90_put_att(imapfile, id_cgwav(iid),   'units'        , 'm s-1')

             ierr = nf90_def_var(imapfile, 'sigmwav',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_sigmwav(iid))
             ierr = nf90_put_att(imapfile, id_sigmwav(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             ierr = nf90_put_att(imapfile, id_sigmwav(iid),   'standard_name', 'sea_surface_wave_mean_frequency')                          ! not CF
             ierr = nf90_put_att(imapfile, id_sigmwav(iid),   'long_name'    , 'mean wave frequency')
             ierr = nf90_put_att(imapfile, id_sigmwav(iid),   'units'        , 'rad s-1')

             !if ( (windmodel.eq.1) .and. (jawsource.eq.1) ) then
             !
             !   ierr = nf90_def_var(imapfile, 'SwE',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_SwE(iid))
             !   ierr = nf90_put_att(imapfile, id_SwE(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             !   ierr = nf90_put_att(imapfile, id_SwE(iid),   'standard_name', 'source_term_wind_on_E')                          ! not CF
             !   ierr = nf90_put_att(imapfile, id_SwE(iid),   'long_name'    , 'source term wind on wave energy')
             !   ierr = nf90_put_att(imapfile, id_SwE(iid),   'units'        , 'J m-2 s-1')
             !
             !   ierr = nf90_def_var(imapfile, 'SwT',  nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_SwT(iid))
             !   ierr = nf90_put_att(imapfile, id_SwT(iid),   'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
             !   ierr = nf90_put_att(imapfile, id_SwT(iid),   'standard_name', 'source_term_wind_on_T')                          ! not CF
             !   ierr = nf90_put_att(imapfile, id_SwT(iid),   'long_name'    , 'source term wind on wave period')
             !   ierr = nf90_put_att(imapfile, id_SwT(iid),   'units'        , 's s-1')
             !
             !endif
           endif

           if ( NUMCONST.eq.0 ) then
              ierr = unc_add_gridmapping_att(imapfile, (/ id_s1(iid), id_taus(iid), id_ucx(iid), id_ucy(iid), id_unorm(iid), id_sa1(iid), id_sed(iid) /), jsferic)   ! add id_ucz(iid)?
           else
              if (allocated(idum)) then
                 deallocate(idum)
              endif
              allocate(idum(7+NUMCONST))
              idum(1:7) = (/ id_s1(iid), id_taus(iid), id_ucx(iid), id_ucy(iid), id_unorm(iid), id_sa1(iid), id_sed(iid) /)
              do j=1,NUMCONST
                 idum(7+j) = id_const(iid,j)
              enddo
              ierr = unc_add_gridmapping_att(imapfile, idum, jsferic)
           endif
           if (kmx > 0) then
              if ( density_is_pressure_dependent() ) then
                  ierr = unc_add_gridmapping_att(imapfile, (/ id_ucz(iid), id_ucxa(iid), id_ucya(iid), id_ww1(iid), id_rho(iid) /), jsferic)
              else
                  ierr = unc_add_gridmapping_att(imapfile, (/ id_ucz(iid), id_ucxa(iid), id_ucya(iid), id_ww1(iid), id_rhop(iid) /), jsferic)
              endif
           endif

           if (jamaptrachy > 0 .and. jatrt == 1) then
               ! Roughness data on net-links
               ierr = nf90_def_var(imapfile, 'cftrt' , nf90_double, (/ id_netlinkdim(iid), id_timedim(iid) /) , id_cftrt(iid))
               if (ifrctypuni == 0) then
                   ierr = nf90_put_att(imapfile, id_cftrt(iid),'long_name'    , 'Chezy roughness from trachytopes')
                   ierr = nf90_put_att(imapfile, id_cftrt(iid),'units'        , 'm0.5s-1')                ! WO: does not follow standard ? (which accepts only integral powers?)
               elseif (ifrctypuni == 1) then
                   ierr = nf90_put_att(imapfile, id_cftrt(iid),'long_name'    , 'Manning roughness from trachytopes')
                   ierr = nf90_put_att(imapfile, id_cftrt(iid),'units'        , 'sm-0.333')               ! WO: does not follow standard ? (which accepts only integral powers?)
               elseif ((ifrctypuni == 2) .or. (ifrctypuni == 3)) then
                   ierr = nf90_put_att(imapfile, id_cftrt(iid),'long_name'    , 'White-Colebrook roughness from trachytopes')
                   ierr = nf90_put_att(imapfile, id_cftrt(iid),'units'        , 'm')
               else
                   ierr = nf90_put_att(imapfile, id_cftrt(iid),'long_name'    , 'roughness from trachytopes')
                   ierr = nf90_put_att(imapfile, id_cftrt(iid),'units'        , ' ')
               endif
           endif

           if (jamapcali > 0 .and. jacali == 1) then
               ! Calibration factor for roughness data on net-links
               ierr = nf90_def_var(imapfile, 'cfcl' , nf90_double, (/ id_netlinkdim(iid), id_timedim(iid) /) , id_cfcl(iid))
               ierr = nf90_put_att(imapfile, id_cfcl(iid),'long_name'    , 'Calibration factor for roughness')
               ierr = nf90_put_att(imapfile, id_cfcl(iid),'units'        , ' ')
           endif

           if (jamap_chezy_elements > 0) then
               ! Chezy data on flow-nodes
               ierr = nf90_def_var(imapfile, 'czs' , nf90_double, (/ id_flowelemdim(iid), id_timedim(iid) /) , id_czs(iid))
               ierr = nf90_put_att(imapfile, id_czs(iid),'long_name'    , 'Chezy roughness')
               ierr = nf90_put_att(imapfile, id_czs(iid),'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
               ierr = nf90_put_att(imapfile, id_czs(iid),'units'        , 'm0.5s-1')                ! WO: does not follow standard ? (which accepts only integral powers?)
           endif
           if (jamap_chezy_links > 0) then
               ! Chezy data on flow-links
               ierr = nf90_def_var(imapfile, 'czu' , nf90_double, (/ id_flowlinkdim(iid), id_timedim(iid) /) , id_czu(iid))
               ierr = nf90_put_att(imapfile, id_czu(iid),'long_name'    , 'Chezy roughness on flow links')
               ierr = nf90_put_att(imapfile, id_czu(iid),'coordinates'  , 'FlowLink_xu FlowLink_yu')
               ierr = nf90_put_att(imapfile, id_czu(iid),'units'        , 'm0.5s-1')
           endif

           ! 1D2D boundaries
           if (nbnd1d2d > 0) then
               ierr = nf90_def_var(imapfile, '1d2d_flowlinknrs' , nf90_int, (/ id_1d2ddim(iid) /) , id_1d2d_edges(iid))
               ierr = nf90_put_att(imapfile, id_czs(iid),'long_name'    , 'flow link numbers of the open 1D2D boundary links')

               ierr = nf90_def_var(imapfile, '1d2d_zeta' , nf90_double, (/ id_1d2ddim(iid), id_timedim(iid) /) , id_1d2d_zeta1d(iid))
               ierr = nf90_put_att(imapfile, id_1d2d_zeta1d(iid),'standard_name', 'sea_surface_height_above_geoid')
               ierr = nf90_put_att(imapfile, id_1d2d_zeta1d(iid),'long_name'    , '1D water level next to each 1d2d boundary link')
               ierr = nf90_put_att(imapfile, id_1d2d_zeta1d(iid),'units'        , 'm')

               ierr = nf90_def_var(imapfile, '1d2d_crest_level' , nf90_double, (/ id_1d2ddim(iid), id_timedim(iid) /) , id_1d2d_crest_level(iid))
               ierr = nf90_put_att(imapfile, id_1d2d_crest_level(iid),'standard_name', 'sea_surface_height_above_geoid')
               ierr = nf90_put_att(imapfile, id_1d2d_crest_level(iid),'long_name'    , 'crest level of 1d2d boundary link')
               ierr = nf90_put_att(imapfile, id_1d2d_crest_level(iid),'units'        , 'm')

               ierr = nf90_def_var(imapfile, '1d2d_b_2di' , nf90_double, (/ id_1d2ddim(iid), id_timedim(iid) /) , id_1d2d_b_2di(iid))
               ierr = nf90_put_att(imapfile, id_1d2d_b_2di(iid),'standard_name', 'b_2di')
               ierr = nf90_put_att(imapfile, id_1d2d_b_2di(iid),'long_name'    , 'coefficient for 1d2d interface b_2di')
               ierr = nf90_put_att(imapfile, id_1d2d_b_2di(iid),'units'        , '-')

               ierr = nf90_def_var(imapfile, '1d2d_b_2dv' , nf90_double, (/ id_1d2ddim(iid), id_timedim(iid) /) , id_1d2d_b_2dv(iid))
               ierr = nf90_put_att(imapfile, id_1d2d_b_2dv(iid),'standard_name', 'b_2dv')
               ierr = nf90_put_att(imapfile, id_1d2d_b_2dv(iid),'long_name'    , 'coefficient for 1d2d interface b_2di')
               ierr = nf90_put_att(imapfile, id_1d2d_b_2dv(iid),'units'        , '-')

               ierr = nf90_def_var(imapfile, '1d2d_d_2dv' , nf90_double, (/ id_1d2ddim(iid), id_timedim(iid) /) , id_1d2d_d_2dv(iid))
               ierr = nf90_put_att(imapfile, id_1d2d_d_2dv(iid),'standard_name', 'd_2dv')
               ierr = nf90_put_att(imapfile, id_1d2d_d_2dv(iid),'long_name'    , 'coefficient for 1d2d interface d_2dv')
               ierr = nf90_put_att(imapfile, id_1d2d_d_2dv(iid),'units'        , '-')

               ierr = nf90_def_var(imapfile, '1d2d_qzeta' , nf90_double, (/ id_1d2ddim(iid), id_timedim(iid) /) , id_1d2d_q_zeta(iid))
               ierr = nf90_put_att(imapfile, id_1d2d_q_zeta(iid),'standard_name', 'q_zeta_1d2d')
               ierr = nf90_put_att(imapfile, id_1d2d_q_zeta(iid),'long_name'    , 'q_zeta_1d2d')
               ierr = nf90_put_att(imapfile, id_1d2d_q_zeta(iid),'units'        , 'm2 s-1')

               ierr = nf90_def_var(imapfile, '1d2d_q_lat' , nf90_double, (/ id_1d2ddim(iid), id_timedim(iid) /) , id_1d2d_q_lat(iid))
               ierr = nf90_put_att(imapfile, id_1d2d_q_lat(iid),'standard_name', 'q_lat')
               ierr = nf90_put_att(imapfile, id_1d2d_q_lat(iid),'long_name'    , 'q_lat')
               ierr = nf90_put_att(imapfile, id_1d2d_q_lat(iid),'units'        , 'm3 s-1')

               ierr = nf90_def_var(imapfile, '1d2d_cfl' , nf90_double, (/ id_1d2ddim(iid), id_timedim(iid) /) , id_1d2d_cfl(iid))
               ierr = nf90_put_att(imapfile, id_1d2d_cfl(iid),'standard_name', 'cfl')
               ierr = nf90_put_att(imapfile, id_1d2d_cfl(iid),'long_name'    , 'wave flow courant')
               ierr = nf90_put_att(imapfile, id_1d2d_cfl(iid),'units'        , '-')

               ierr = nf90_def_var(imapfile, '1d2d_sb' , nf90_double, (/ id_1d2ddim(iid), id_timedim(iid) /) , id_1d2d_sb(iid))
               ierr = nf90_put_att(imapfile, id_1d2d_sb(iid),'standard_name', '1d2d_sb')
               ierr = nf90_put_att(imapfile, id_1d2d_sb(iid),'long_name'    , 'water levels in boundary points')
               ierr = nf90_put_att(imapfile, id_1d2d_sb(iid),'units'        , 'm')

               ierr = nf90_def_var(imapfile, '1d2d_s0_2d' , nf90_double, (/ id_1d2ddim(iid), id_timedim(iid) /) , id_1d2d_s0_2d(iid))
               ierr = nf90_put_att(imapfile, id_1d2d_s0_2d(iid),'standard_name', '1d2d_s0_2d')
               ierr = nf90_put_att(imapfile, id_1d2d_s0_2d(iid),'long_name'    , 'water levels on interface at previous time step')
               ierr = nf90_put_att(imapfile, id_1d2d_s0_2d(iid),'units'        , 'm')

               ierr = nf90_def_var(imapfile, '1d2d_s1_2d' , nf90_double, (/ id_1d2ddim(iid), id_timedim(iid) /) , id_1d2d_s1_2d(iid))
               ierr = nf90_put_att(imapfile, id_1d2d_s1_2d(iid),'standard_name', '1d2d_s1_2d')
               ierr = nf90_put_att(imapfile, id_1d2d_s1_2d(iid),'long_name'    , 'water levels on interface at current time step')
               ierr = nf90_put_att(imapfile, id_1d2d_s1_2d(iid),'units'        , 'm')

               ierr = nf90_def_var(imapfile, '1d2d_flow_cond' , nf90_int, (/ id_1d2ddim(iid), id_timedim(iid) /) , id_1d2d_flow_cond(iid))
               ierr = nf90_put_att(imapfile, id_1d2d_flow_cond(iid),'standard_name', 'flow_condition')
               ierr = nf90_put_att(imapfile, id_1d2d_flow_cond(iid),'long_name'    , 'flow Condition 0: closed, 1: free 1d to 2d, 2: free 2d to 1d, 3: submerged')
               ierr = nf90_put_att(imapfile, id_1d2d_flow_cond(iid),'units'        , '-')

           endif
        endif

        if (jamapwind > 0 .and. japatm > 0) then
            call definencvar(imapfile,id_patm(iid)   ,nf90_double,idims, 'Patm'  , 'Atmospheric Pressure', 'N m-2', 'FlowElem_xcc FlowElem_ycc')
        endif

        if (ice_mapout) then
            call definencvar(imapfile,id_ice_af(iid)  ,nf90_double,idims, 'ice_af' , 'Fraction of the surface area covered by floating ice', '1', 'FlowElem_xcc FlowElem_ycc')
            call definencvar(imapfile,id_ice_h(iid)   ,nf90_double,idims, 'ice_h'  , 'Thickness of floating ice cover', 'm', 'FlowElem_xcc FlowElem_ycc')
            call definencvar(imapfile,id_ice_p(iid)   ,nf90_double,idims, 'ice_p'  , 'Pressure exerted by the floating ice cover', 'N m-2', 'FlowElem_xcc FlowElem_ycc')
            if (ja_icecover == ICECOVER_SEMTNER) then
               call definencvar(imapfile,id_ice_t(iid)   ,nf90_double,idims, 'ice_t'  , 'Temperature of the floating ice cover', 'degC', 'FlowElem_xcc FlowElem_ycc')
               call definencvar(imapfile,id_snow_h(iid)  ,nf90_double,idims, 'snow_h'  , 'Thickness of the snow layer', 'm', 'FlowElem_xcc FlowElem_ycc')
               call definencvar(imapfile,id_snow_t(iid)  ,nf90_double,idims, 'snow_t'  , 'Temperature of the snow layer', 'degC', 'FlowElem_xcc FlowElem_ycc')
            endif
        endif

        if ((jamapwind > 0 .or. jamapwindstress > 0 .or. jaseparate_==2) .and. jawind /= 0) then
           if (jawindstressgiven == 0 .or. jaseparate_==2) then
              ierr = nf90_def_var(imapfile, 'windx', nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_windx(iid))
              ierr = nf90_def_var(imapfile, 'windy', nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_windy(iid))
           else
              ierr = nf90_def_var(imapfile, 'windstressx', nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_windx(iid))
              ierr = nf90_def_var(imapfile, 'windstressy', nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_windy(iid))
           endif

           ierr = nf90_put_att(imapfile, id_windx(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
           if (jawindstressgiven == 0 .or. jaseparate_==2) then
              if (jsferic == 0 ) then
                 ierr = nf90_put_att(imapfile, id_windx(iid),  'standard_name', 'x_wind')
                 ierr = nf90_put_att(imapfile, id_windx(iid),  'long_name'    , 'velocity of air on flow element center, x-component')
              else
                 ierr = nf90_put_att(imapfile, id_windx(iid),  'standard_name', 'eastward_wind')
                 ierr = nf90_put_att(imapfile, id_windx(iid),  'long_name'    , 'eastward air velocity on flow element center, x-component')
              endif
              ierr = nf90_put_att(imapfile, id_windx(iid),  'units'        , 'm s-1')
           else
              if (jsferic == 0 ) then
                 ierr = nf90_put_att(imapfile, id_windx(iid),  'standard_name', 'x_windstress')
                 ierr = nf90_put_att(imapfile, id_windx(iid),  'long_name'    , 'windstress on flow element center, x-component')
              else
                 ierr = nf90_put_att(imapfile, id_windx(iid),  'standard_name', 'eastward_windstress')
                 ierr = nf90_put_att(imapfile, id_windx(iid),  'long_name'    , 'eastward windstress on flow element center, x-component')
              endif
              ierr = nf90_put_att(imapfile, id_windx(iid),  'units'        , 'N m-2')
           endif

           ierr = nf90_put_att(imapfile, id_windy(iid),  'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
           if (jawindstressgiven == 0 .or. jaseparate_==2) then
              if (jsferic == 0 ) then
                 ierr = nf90_put_att(imapfile, id_windy(iid),  'standard_name', 'y_wind')
                 ierr = nf90_put_att(imapfile, id_windy(iid),  'long_name'    , 'velocity of air on flow element center, y-component')
              else
                 ierr = nf90_put_att(imapfile, id_windy(iid),  'standard_name', 'northward_wind')
                 ierr = nf90_put_att(imapfile, id_windy(iid),  'long_name'    , 'northward air velocity on flow element center, y-component')
              endif
              ierr = nf90_put_att(imapfile, id_windy(iid),  'units'        , 'm s-1')
           else
              if (jsferic == 0 ) then
                 ierr = nf90_put_att(imapfile, id_windy(iid),  'standard_name', 'y_windstress')
                 ierr = nf90_put_att(imapfile, id_windy(iid),  'long_name'    , 'windstress air on flow element center, y-component')
              else
                 ierr = nf90_put_att(imapfile, id_windy(iid),  'standard_name', 'northward_windstress')
                 ierr = nf90_put_att(imapfile, id_windy(iid),  'long_name'    , 'northward windstress on flow element center, y-component')
              endif
              ierr = nf90_put_att(imapfile, id_windy(iid),  'units'        , 'N m-2')
           endif
        endif

        if (jamapwind > 0 .and. jawind /= 0 .and. jawindstressgiven == 0) then
           ! Also wind on flow links
           ierr = nf90_def_var(imapfile, 'windxu', nf90_double, (/ id_flowlinkdim(iid), id_timedim (iid)/) , id_windxu(iid))
           ierr = nf90_def_var(imapfile, 'windyu', nf90_double, (/ id_flowlinkdim(iid), id_timedim (iid)/) , id_windyu(iid))

           ierr = nf90_put_att(imapfile, id_windxu(iid),  'coordinates'  , 'FlowLink_xu FlowLink_yu')
           if (jsferic == 0) then
              ierr = nf90_put_att(imapfile, id_windxu(iid),  'long_name'    , 'velocity of air on flow links, x-component')
              ierr = nf90_put_att(imapfile, id_windxu(iid),  'standard_name', 'x_velocity_wind')
           else
              ierr = nf90_put_att(imapfile, id_windxu(iid),  'long_name'    , 'eastward air velocity on flow links, x-component')
              ierr = nf90_put_att(imapfile, id_windxu(iid),  'standard_name', 'eastward_wind')
           endif
           ierr = nf90_put_att(imapfile, id_windxu(iid),  'units'        , 'm s-1')

           ierr = nf90_put_att(imapfile, id_windyu(iid),  'coordinates'  , 'FlowLink_xu FlowLink_yu')
           if (jsferic == 0) then
              ierr = nf90_put_att(imapfile, id_windyu(iid),  'long_name'    , 'velocity of air on flow links, y-component')
              ierr = nf90_put_att(imapfile, id_windyu(iid),  'standard_name', 'y_velocity_wind')
           else
              ierr = nf90_put_att(imapfile, id_windyu(iid),  'long_name'    , 'northward air velocity on flow links, y-component')
              ierr = nf90_put_att(imapfile, id_windyu(iid),  'standard_name', 'northward_wind')
           endif
           ierr = nf90_put_att(imapfile, id_windyu(iid),  'units'        , 'm s-1')
        endif
        !
        ierr = unc_add_gridmapping_att(imapfile, (/ id_windx(iid), id_windy(iid), id_windxu(iid), id_windyu(iid),  nf90_global /), jsferic)

        if (javeg > 0) then
           ierr = nf90_def_var(imapfile, 'rnveg', nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_rnveg(iid))
           ierr = nf90_put_att(imapfile, id_rnveg(iid),'long_name'    , 'Stem density of vegetation')
           ierr = nf90_put_att(imapfile, id_rnveg(iid),'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
           ierr = nf90_put_att(imapfile, id_rnveg(iid),'units'        , 'm-2')

           ierr = nf90_def_var(imapfile, 'diaveg', nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_diaveg(iid))
           ierr = nf90_put_att(imapfile, id_diaveg(iid),'long_name'    , 'Stem diameter of vegetation')
           ierr = nf90_put_att(imapfile, id_diaveg(iid),'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
           ierr = nf90_put_att(imapfile, id_diaveg(iid),'units'        , 'm')

           ierr = nf90_def_var(imapfile, 'veg_stemheight', nf90_double, (/ id_flowelemdim(iid), id_timedim (iid)/) , id_veg_stemheight(iid))
           ierr = nf90_put_att(imapfile, id_veg_stemheight(iid),'long_name'    , 'Stem height of vegetation')
           ierr = nf90_put_att(imapfile, id_veg_stemheight(iid),'coordinates'  , 'FlowElem_xcc FlowElem_ycc')
           ierr = nf90_put_att(imapfile, id_veg_stemheight(iid),'units'        , 'm')
        endif

        ! For all 3D variables, expand the coordinate attribute with a vertical coordinate
        ierr = nf90_inq_varid( imapfile, 'LayCoord_cc', varid)
        if (ierr==NF90_NOERR) then
           zcc_elem = 'LayCoord_cc'
           zw_elem = 'LayCoord_w'
           zu_link = 'LayCoord_cc'   ! z/sigma coords are the same for u-positions and cc-positions.
           zwu_link = 'LayCoord_w'
        else
           zcc_elem = 'FlowElem_zcc'
           zw_elem = 'FlowElem_zw'
           zu_link = ''                        ! To be added, Issue UNST-4880
           zwu_link = ''
        endif
        if (nf90_inquire(imapfile, nVariables=varid)==NF90_NOERR) then
           do while (varid>0)
              if (nf90_inquire_variable(imapfile, varid, ndims=ndims)==NF90_NOERR) then
                 call realloc(idum,ndims,keepexisting=.False.)
                 if (nf90_inquire_variable(imapfile, varid, dimids=idum)==NF90_NOERR) then
                    if (any(idum==id_wdim(iid)) .and. any(idum==id_flowelemdim(iid))) then
                        ierr = ncu_append_atts( imapfile, varid, 'coordinates', trim(zw_elem))
                    endif
                    if (any(idum==id_laydim(iid)) .and. any(idum==id_flowelemdim(iid))) then
                        ierr = ncu_append_atts( imapfile, varid, 'coordinates', trim(zcc_elem))
                    endif
                    if (any(idum==id_wdim(iid)) .and. any(idum==id_flowlinkdim(iid))) then
                        ierr = ncu_append_atts( imapfile, varid, 'coordinates', trim(zwu_link))
                    endif
                    if (any(idum==id_laydim(iid)) .and. any(idum==id_flowelemdim(iid))) then
                        ierr = ncu_append_atts( imapfile, varid, 'coordinates', trim(zu_link))
                    endif
                 endif
              endif
              varid = varid - 1
           enddo
        endif


        ierr = nf90_enddef(imapfile)

        ! 1D2D boundaries
        if (nbnd1d2d > 0 .and. jaseparate_ /= 2) then
           if (allocated(idum)) then 
              deallocate(idum)
           endif
           allocate(idum(nbnd1d2d))
           do i=1,nbnd1d2d
              idum(i) = kbnd1d2d(3, i) ! Flow link nrs
           enddo
           ierr = nf90_put_var(imapfile, id_1d2d_edges(iid), idum)
           deallocate(idum)
        endif

        if (nomba > 0) then
           ierr = nf90_put_var(imapfile, id_mba(iid), mbadef(1:NdxNdxi))
        endif

        firststep(iid) = .false.

    endif

    ! End of writing time-independent flow geometry data.
    ! -- Inquire id's belonging to map file ------------------------
    if (firststep(iid) .and. ndim>0) then ! TODO: AvD: UNST-530
       !
       !
       ! this step is necessary because if a snapshot_map.nc file is written
       ! in between two map file outputs the saved id's may have changed
       !
       firststep(iid) = .false.
       !
       ierr = nf90_inq_dimid(imapfile, 'nFlowElem', id_flowelemdim(iid))
       if (ierr /= nf90_noerr) then
          ierr = nf90_inq_dimid(imapfile, 'nFlowElemWithBnd', id_flowelemdim(iid))
       endif

       ierr = nf90_inq_dimid(imapfile, 'nFlowLink', id_flowlinkdim(iid))
       !
       ! Time
       ierr = nf90_inq_dimid(imapfile, 'time', id_timedim(iid))
       ierr = nf90_inq_varid(imapfile, 'time', id_time(iid))
       !
       if ( kmx>0 ) then
          ierr = nf90_inq_dimid(imapfile, 'laydim', id_laydim(iid))
          ierr = nf90_inq_dimid(imapfile, 'wdim', id_wdim(iid))
       endif
       !
       ! Size of latest timestep

       ! Why ask for id_*, they are in a save statement no?

       ierr = nf90_inq_varid(imapfile, 'timestep', id_timestep(iid))
       ierr = nf90_inq_varid(imapfile, 'taus' ,  id_taus(iid))
       !
       if ( kmx>0 ) then     !  3D
          ierr = nf90_inq_varid(imapfile, 'ucx', id_ucx(iid))
          ierr = nf90_inq_varid(imapfile, 'ucy', id_ucy(iid))
          ierr = nf90_inq_varid(imapfile, 'ucz', id_ucz(iid))
          ierr = nf90_inq_varid(imapfile, 'ucxa', id_ucxa(iid))
          ierr = nf90_inq_varid(imapfile, 'ucya', id_ucya(iid))
          ierr = nf90_inq_varid(imapfile, 'ww1', id_ww1(iid))
          if ( density_is_pressure_dependent() ) then
            ierr = nf90_inq_varid(imapfile, 'density', id_rho(iid))
          else
            ierr = nf90_inq_varid(imapfile, 'rho', id_rhop(iid))
          endif
          if ( iturbulencemodel >= 3 ) then
             ierr = nf90_inq_varid(imapfile, 'turkin1', id_turkin1(iid))
             ierr = nf90_inq_varid(imapfile, 'tureps1', id_tureps1(iid))
             ierr = nf90_inq_varid(imapfile, 'vicwwu' , id_vicwwu(iid) )
          endif
        else
          ierr = nf90_inq_varid(imapfile, 'ucx', id_ucx(iid))
          ierr = nf90_inq_varid(imapfile, 'ucy', id_ucy(iid))
          ierr = nf90_inq_varid(imapfile, 'spircrv', id_spircrv(iid))
          ierr = nf90_inq_varid(imapfile, 'spirint', id_spirint(iid))
        endif
        !
        if (jasal > 0) then
           ierr = nf90_inq_varid(imapfile, 'sa1', id_sa1(iid))
        endif

        if (jatem > 0) then
           ierr = nf90_inq_varid(imapfile, 'tem1', id_tem1(iid))
        endif

        if (ITRA1 > 0) then
           do j=ITRA1,ITRAN
              tmpstr = const_names(j)
              ! Forbidden chars in NetCDF names: space, /, and more.
              call replace_char(tmpstr,32,95)
              call replace_char(tmpstr,47,95)
              ierr = nf90_inq_varid(imapfile, trim(tmpstr), id_const(iid,j))
           enddo
        endif

        !
        if (stm_included) then
           ierr = nf90_inq_varid(imapfile, 'nSedTot', id_sedtotdim(iid))
           ierr = nf90_inq_varid(imapfile, 'nSedSus', id_sedsusdim(iid))
           ierr = nf90_inq_varid(imapfile, 'nBedLayers', id_nlyrdim(iid))

           if (stmpar%lsedsus > 0) then
              ierr = nf90_inq_varid(imapfile, 'ws', id_ws(iid))
              !
              ! equilibrium concentration, 2D only
              if (kmx == 0) then
                 ierr = nf90_inq_varid(imapfile, 'rsedeq', id_rsedeq(iid))
              endif

              if (stmpar%morpar%moroutput%sourcesink) then
                 ierr = nf90_inq_varid(imapfile, 'sourse', id_sourse(iid))
                 ierr = nf90_inq_varid(imapfile, 'sinkse', id_sinkse(iid))
              endif

              if (stmpar%morpar%moroutput%suvcor) then
                 ierr = nf90_inq_varid(imapfile, 'e_scrn', id_scrn(iid))
                 !ierr = nf90_inq_varid(imapfile, 'e_scrt', id_scrt(iid))
              endif

              if (stmpar%morpar%moroutput%aks) then
                 ierr = nf90_inq_varid(imapfile, 'aks', id_aks(iid))
                 ierr = nf90_inq_varid(imapfile, 'rca', id_rca(iid))
              endif
              !
              ! Suspended fractions
              if (stmpar%lsedsus .gt. 0) then
                 do j=ISED1,ISEDN
                    tmpstr = const_names(j)
                    ! Forbidden chars in NetCDF names: space, /, and more.
                    call replace_char(tmpstr,32,95)
                    call replace_char(tmpstr,47,95)
                    ierr = nf90_inq_varid(imapfile, trim(tmpstr), id_const(iid,j))
                 enddo
              endif
           endif

           if (stmpar%morpar%moroutput%dzduuvv) then ! bedslope
              ierr = nf90_inq_varid(imapfile, 'e_dzdn', id_dzdn(iid))
              ierr = nf90_inq_varid(imapfile, 'e_dzdt', id_dzdt(iid))
           endif

           if (stmpar%morpar%moroutput%umod) then
              ierr = nf90_inq_varid(imapfile, 'umod', id_umod(iid))
           endif

           if (stmpar%morpar%moroutput%zumod) then
              ierr = nf90_inq_varid(imapfile, 'zumod', id_zumod(iid))
           endif

           if (stmpar%morpar%moroutput%ustar) then
              ierr = nf90_inq_varid(imapfile, 'ustar', id_ustar(iid))
           endif

           if (stmpar%morpar%moroutput%sbcuv) then
              ierr = nf90_inq_varid(imapfile, 'sbcx', id_sbcx(iid))
              ierr = nf90_inq_varid(imapfile, 'sbcy', id_sbcy(iid))
              ierr = nf90_inq_varid(imapfile, 'sbcx_reconstructed', id_sbcx_reconstructed(iid))
              ierr = nf90_inq_varid(imapfile, 'sbcy_reconstructed', id_sbcy_reconstructed(iid))
           endif

           if (stmpar%morpar%moroutput%sbwuv) then
              ierr = nf90_inq_varid(imapfile, 'sbwx', id_sbwx(iid))
              ierr = nf90_inq_varid(imapfile, 'sbwy', id_sbwy(iid))
              ierr = nf90_inq_varid(imapfile, 'sbwx_reconstructed', id_sbwx_reconstructed(iid))
              ierr = nf90_inq_varid(imapfile, 'sbwy_reconstructed', id_sbwy_reconstructed(iid))
           endif

           if (stmpar%morpar%moroutput%sswuv) then
              ierr = nf90_inq_varid(imapfile, 'sswx', id_sswx(iid))
              ierr = nf90_inq_varid(imapfile, 'sswy', id_sswy(iid))
              ierr = nf90_inq_varid(imapfile, 'sswx_reconstructed', id_sswx_reconstructed(iid))
              ierr = nf90_inq_varid(imapfile, 'sswy_reconstructed', id_sswy_reconstructed(iid))
           endif

           if (stmpar%morpar%moroutput%sscuv) then
              ierr = nf90_inq_varid(imapfile, 'sscx', id_sscx(iid))
              ierr = nf90_inq_varid(imapfile, 'sscy', id_sscy(iid))
              ierr = nf90_inq_varid(imapfile, 'sscx_reconstructed', id_sscx_reconstructed(iid))
              ierr = nf90_inq_varid(imapfile, 'sscy_reconstructed', id_sscy_reconstructed(iid))
           endif

           ierr = nf90_inq_varid(imapfile, 'sxtot', id_sxtot(iid))
           ierr = nf90_inq_varid(imapfile, 'sytot', id_sytot(iid))

           ierr = nf90_inq_varid(imapfile, 'mor_bl', id_morbl(iid))

           select case (stmpar%morlyr%settings%iunderlyr)
           case (1)
              ierr = nf90_inq_varid(imapfile, 'bodsed', id_bodsed(iid))
              ierr = nf90_inq_varid(imapfile, 'dpsed', id_dpsed(iid))
           case (2)
              ierr = nf90_inq_varid(imapfile, 'msed', id_msed(iid))
              ierr = nf90_inq_varid(imapfile, 'lyrfrac', id_lyrfrac(iid))
              ierr = nf90_inq_varid(imapfile, 'thlyr', id_thlyr(iid))
              if (stmpar%morlyr%settings%iporosity>0) then
                 ierr = nf90_inq_varid(imapfile, 'poros', id_poros(iid))
              endif
           end select
        !
           if (stmpar%morpar%moroutput%taurat) then
              ierr = nf90_inq_varid(imapfile, 'taurat' ,id_taurat(iid))
           endif
           if (stmpar%morpar%moroutput%dm) then
              ierr = nf90_inq_varid(imapfile, 'dm' ,id_dm(iid))
           endif
           if (stmpar%morpar%moroutput%dg) then
              ierr = nf90_inq_varid(imapfile, 'dg' ,id_dg(iid))
           endif
           if (stmpar%morpar%moroutput%dgsd) then
              ierr = nf90_inq_varid(imapfile, 'dgsd' ,id_dgsd(iid))
           endif
           if (stmpar%morpar%moroutput%percentiles) then
              do l = 1, stmpar%morpar%nxx
                 write(dxname,'(A,I2.2)') 'DXX',l
                 ierr = nf90_inq_varid(imapfile, dxname ,id_dxx(l,iid))
              enddo
           endif
           if (stmpar%morpar%moroutput%frac) then
              ierr = nf90_inq_varid(imapfile, 'frac' ,id_frac(iid))
           endif
           if (stmpar%morpar%moroutput%mudfrac) then
              ierr = nf90_inq_varid(imapfile, 'mudfrac' ,id_mudfrac(iid))
           endif
           if (stmpar%morpar%moroutput%sandfrac) then
              ierr = nf90_inq_varid(imapfile, 'sandfrac' ,id_sandfrac(iid))
           endif
           if (stmpar%morpar%moroutput%fixfac) then
              ierr = nf90_inq_varid(imapfile, 'fixfac' ,id_fixfac(iid))
           endif
           if (stmpar%morpar%moroutput%hidexp) then
              ierr = nf90_inq_varid(imapfile, 'hidexp' ,id_hidexp(iid))
           endif
           ! Fluff layers
           if (stmpar%morpar%flufflyr%iflufflyr>0 .and. stmpar%lsedsus>0) then
              ierr = nf90_inq_varid(imapfile, 'mfluff' ,id_mfluff(iid))
           endif
        endif

        if (bfmpar%lfbedfrmout) then
           if (bfmpar%lfbedfrm) then
              ierr = nf90_inq_varid(imapfile, 'duneheight' ,id_duneheight(iid))
              ierr = nf90_inq_varid(imapfile, 'dunelength' ,id_dunelength(iid))
           endif
           if (bfmpar%lfbedfrmrou) then
              ierr = nf90_inq_varid(imapfile, 'ksr' ,id_ksr(iid))
              ierr = nf90_inq_varid(imapfile, 'ksmr' ,id_ksmr(iid))
              ierr = nf90_inq_varid(imapfile, 'ksd' ,id_ksd(iid))
              ierr = nf90_inq_varid(imapfile, 'ks' ,id_ks(iid))
           endif
        endif
        !
        if (jased > 0 .and. .not.stm_included) then
           ierr = nf90_inq_dimid(imapfile, 'nFrac', id_maxfracdim(iid))
           if (jaceneqtr == 1) then
              ierr = nf90_inq_dimid(imapfile, 'nFlowElem', id_erolaydim(iid)) ! Note: points to an existing dimension (either nNetNode, or nFlowElem)
              if (ierr /= nf90_noerr) then
                 ierr = nf90_inq_dimid(imapfile, 'nFlowElemWithBnd', id_erolaydim(iid))
              endif
           else
              ierr = nf90_inq_dimid(imapfile, 'nNetNode', id_erolaydim(iid)) ! Note: points to an existing dimension (either nNetNode, or nFlowElem)
           endif
           !
           ierr = nf90_inq_varid(imapfile, 'sed', id_sed(iid))
           !
           ierr = nf90_inq_varid(imapfile, 'ero', id_ero(iid))
        endif
        !
        ! JRE - XBeach
        if (jawave .eq. 4) then
           ierr = nf90_inq_varid(imapfile, 'E'        , id_E(iid))
           ierr = nf90_inq_varid(imapfile, 'R'        , id_R(iid))
           ierr = nf90_inq_varid(imapfile, 'H'        , id_H(iid))
           ierr = nf90_inq_varid(imapfile, 'D'        , id_D(iid))
           ierr = nf90_inq_varid(imapfile, 'DR'       , id_DR(iid))
           ierr = nf90_inq_varid(imapfile, 'urms'     , id_urms(iid))
           ierr = nf90_inq_varid(imapfile, 'urms_cc'  , id_urmscc(iid))
           ierr = nf90_inq_varid(imapfile, 'ust'      , id_ust(iid))
           ierr = nf90_inq_varid(imapfile, 'vst'      , id_vst(iid))
           ierr = nf90_inq_varid(imapfile, 'Fx_cc'    , id_Fxcc(iid))
           ierr = nf90_inq_varid(imapfile, 'Fy_cc'    , id_Fycc(iid))

           ierr = nf90_inq_varid(imapfile, 'thetamean', id_thetamean(iid))
           ierr = nf90_inq_varid(imapfile, 'cwav'     , id_cwav(iid))
           ierr = nf90_inq_varid(imapfile, 'cgwav'    , id_cgwav(iid))
           ierr = nf90_inq_varid(imapfile, 'sigmwav'  , id_sigmwav(iid))

           !if ( (windmodel .eq. 1) .and. (jawsource .eq. 1) ) then
           !   ierr = nf90_inq_varid(imapfile, 'SwE'  , id_SwE(iid))
           !   ierr = nf90_inq_varid(imapfile, 'SwT'  , id_SwT(iid))
           !endif

        endif

        ! 1D2D boundaries
        if (nbnd1d2d > 0) then
           ierr = nf90_inq_varid(imapfile, '1d2d_flowlinknrs' , id_1d2d_edges(iid))
           ierr = nf90_inq_varid(imapfile, '1d2d_zeta'        , id_1d2d_zeta1d(iid))
           ierr = nf90_inq_varid(imapfile, '1d2d_crest_level' , id_1d2d_crest_level(iid))
           ierr = nf90_inq_varid(imapfile, '1d2d_b_2di'       , id_1d2d_b_2di(iid))
           ierr = nf90_inq_varid(imapfile, '1d2d_b_2dv'       , id_1d2d_b_2dv(iid))
           ierr = nf90_inq_varid(imapfile, '1d2d_d_2dv'       , id_1d2d_d_2dv(iid))
           ierr = nf90_inq_varid(imapfile, '1d2d_q_zeta'      , id_1d2d_q_zeta(iid))
           ierr = nf90_inq_varid(imapfile, '1d2d_q_lat'       , id_1d2d_q_lat(iid))
           ierr = nf90_inq_varid(imapfile, '1d2d_cfl'         , id_1d2d_cfl(iid))
           ierr = nf90_inq_varid(imapfile, '1d2d_sb'          , id_1d2d_sb(iid))
           ierr = nf90_inq_varid(imapfile, 'id_1d2d_s0_2d'    , id_1d2d_s0_2d(iid))
           ierr = nf90_inq_varid(imapfile, 'id_1d2d_s1_2d'    , id_1d2d_s1_2d(iid))
           ierr = nf90_inq_varid(imapfile, '1d2d_flow_cond'   , id_1d2d_flow_cond(iid))
        endif

        if ( jamaptidep.eq.1 .and. jatidep > 0 ) then
           if ( jaselfal.eq.0 ) then
              ierr = nf90_inq_varid(imapfile, 'TidalPotential', id_tidep(iid))
           else
              ierr = nf90_inq_varid(imapfile, 'TidalPotential_no_SAL', id_tidep(iid))
           endif
           if ( jaselfal.gt.0 ) then
              ierr = nf90_inq_varid(imapfile, 'SALPotential', id_salp(iid))
           endif
        endif

        if ( jamapIntTidesDiss.eq.1 .and. jaFrcInternalTides2D.gt.0 ) then
           ierr = nf90_inq_varid(imapfile, 'internal_tides_dissipation', id_inttidesdiss(iid))
        endif

        !
        ! Flow data on edges
        ierr = nf90_inq_varid(imapfile, 'unorm' , id_unorm(iid))
        !
        ! Flow data on edges
        ierr = nf90_inq_varid(imapfile, 'u0'    , id_u0(iid))
        ierr = nf90_inq_varid(imapfile, 'q1'    , id_q1(iid))
        ierr = nf90_inq_varid(imapfile, 'viu'   , id_viu(iid))
        ierr = nf90_inq_varid(imapfile, 'diu'   , id_diu(iid))
        !
        if (jawind/=0) then
            ierr = nf90_inq_varid(imapfile, 'windx', id_windx(iid))
            ierr = nf90_inq_varid(imapfile, 'windy', id_windy(iid))
        endif

        if (jaseparate_==2 .and. javeg > 0) then
           ierr = nf90_inq_varid(imapfile, 'rnveg', id_rnveg(iid))
           ierr = nf90_inq_varid(imapfile, 'diaveg', id_diaveg(iid))
           ierr = nf90_inq_varid(imapfile, 'veg_stemheight', id_veg_stemheight(iid))
        endif
    endif

    ! -- Start data writing (flow data) ------------------------
    if (jaseparate_ == 1) then
        itim = 1
        firststep(iid) = .true.
    elseif (jaseparate_ == 2) then
        itim = 1
    else
        it_map   = it_map+1
        itim     = it_map ! Increment time dimension index
    endif

    ! Time
    ierr = nf90_put_var(imapfile, id_time    (iid), tim, (/ itim /))
    ierr = nf90_put_var(imapfile, id_timestep(iid), dts, (/ itim /))

    !
    ! Transform uxy/ucy into Eulerian velocities,
    ! only when the user asks for it and only if we are not writing to com-file
    !
    jaeulerloc = 0
    if (jaeulervel==1 .and. jaseparate_/=2 .and. jawave.gt.0) then
       jaeulerloc = 1
    endif
    !
    call getucxucyeulmag(ndkx, workx, worky, ucmag, jaeulerloc, 0)
    !
    !  Hack to pass time varying bottom levels to SWAN
    !  Also needed for morphostatic runs in 3D
    !
    if (jaseparate_==2) then
       ! JRE: was _zcc, but this has laydim included as dimension, which does not work in 3D
       ierr = nf90_inq_varid(imapfile, 'FlowElem_bl', id_swanbl(iid))
       ierr = nf90_put_var(imapfile, id_swanbl(iid),  -bl,   (/ 1, itim /), (/ ndxndxi, 1 /))
    endif
    !
    ! Water level
    if (jamaps1>0 .or. jaseparate_==2) then
       ierr = nf90_put_var(imapfile, id_s1(iid),  s1,   (/ 1, itim /), (/ ndxndxi, 1 /))
    endif
    !
    if (jamapucvec>0 .or. jaseparate_==2) then
       if ( kmx==0 ) then
          ierr = nf90_put_var(imapfile, id_ucx  (iid), workx,  (/ 1, itim /), (/ ndxndxi, 1 /))
          ierr = nf90_put_var(imapfile, id_ucy  (iid), worky,  (/ 1, itim /), (/ ndxndxi, 1 /))
       endif
    endif

    if ( kmx>0 ) then
       call unc_append_3dflowgeom_put(imapfile, jaseparate_, itim) ! needed for 3D wave coupling on comfile: Flowelem_zw
       if (jamapucvec>0 .or. jaseparate_==2) then
          do kk=1,ndxndxi
             work1(:, kk) = dmiss ! For proper fill values in z-model runs.
             call getkbotktop(kk,kb,kt)
             call getlayerindices(kk, nlayb, nrlay)
             do k = kb,kt
                work1(k-kb+nlayb,kk) = workx(k)
             enddo
          enddo
          ierr = nf90_put_var(imapfile, id_ucx(iid), work1(1:kmx,1:ndxndxi), start=(/ 1, 1, itim /), count=(/ kmx, ndxndxi, 1 /))

          do kk=1,ndxndxi
             work1(:, kk) = dmiss ! For proper fill values in z-model runs.
             call getkbotktop(kk,kb,kt)
             call getlayerindices(kk, nlayb, nrlay)
             do k = kb,kt
                work1(k-kb+nlayb,kk) = worky(k)
             enddo
          enddo
          ierr = nf90_put_var(imapfile, id_ucy(iid), work1(1:kmx,1:ndxndxi), start=(/ 1, 1, itim /), count=(/ kmx, ndxndxi, 1 /))
       endif
    endif

    if (jaseparate_ /= 2) then
        if (jamaps0>0) then
           ierr = nf90_put_var(imapfile, id_s0(iid),  s0,   (/ 1, itim /), (/ ndxndxi, 1 /))
        endif

        if (jamaphs>0) then
           ierr = nf90_put_var(imapfile, id_hs(iid),  hs,   (/ 1, itim /), (/ ndxndxi, 1 /))
        endif
       ! Tau current and chezy roughness
       if (jamaptaucurrent > 0 .or. jamap_chezy_elements > 0 .or. jamap_chezy_links > 0) then
          if (jawave==0) then       ! Else, get taus from subroutine tauwave (taus = f(taucur,tauwave))
             call gettaus(1,1)       ! Update taus and czs
          else if (jamap_chezy_links > 0) then
             call gettaus(2,1)       ! Only update czs
          endif
          if (jawave>0 .and. .not. flowWithoutWaves) then
             call gettauswave(jawaveswartdelwaq)
          endif
       endif
       !
       if (jamap_chezy_links > 0) then
          do LL = 1,lnx
             if (frcu(LL) > 0d0) then
                call getcz (hu(LL), frcu(LL), ifrcutp(LL), czu(LL), LL)
             endif
          enddo
       endif
       !
       if (jamaptaucurrent > 0) then
           ierr = nf90_put_var(imapfile, id_taus(iid), taus,  (/ 1, itim /), (/ ndxndxi, 1 /))
       endif
       !
       if (jamap_chezy_elements > 0) then
           ierr = nf90_put_var(imapfile, id_czs(iid), czs,  (/ 1, itim /), (/ ndxndxi, 1 /))
       endif
       if (jamap_chezy_links > 0) then
           ierr = nf90_put_var(imapfile, id_czu(iid), czu,  (/ 1, itim /), (/ lnx, 1 /))
       endif

       ! Velocities
       if ( kmx>0 ) then
!         3D
          if (jamapucvec>0) then
             call reconstructucz(0)
             !
             do kk=1,ndxndxi
                work1(:, kk) = dmiss ! For proper fill values in z-model runs.
                call getkbotktop(kk,kb,kt)
                call getlayerindices(kk, nlayb, nrlay)
                do k = kb,kt
                   work1(k-kb+nlayb,kk) = ucz(k)
                enddo
             enddo
             ierr = nf90_put_var(imapfile, id_ucz(iid), work1(1:kmx,1:ndxndxi), start=(/ 1, 1, itim /), count=(/ kmx, ndxndxi, 1 /))

             ierr = nf90_put_var(imapfile, id_ucxa(iid), ucxq(1:ndxndxi), start=(/ 1, itim /), count=(/ ndxndxi, 1 /))
             ierr = nf90_put_var(imapfile, id_ucya(iid), ucyq(1:ndxndxi), start=(/ 1, itim /), count=(/ ndxndxi, 1 /))
          endif

          if (jamapww1 > 0) then
             do kk=1,ndxndxi
                work0(:, kk) = dmiss ! For proper fill values in z-model runs.
                call getkbotktop(kk,kb,kt)
                call getlayerindices(kk, nlayb, nrlay)
                do k = kb-1,kt
                   work0(k-kb+nlayb,kk) = ww1(k)
                enddo
             enddo
             ierr = nf90_put_var(imapfile, id_ww1(iid), work0(0:kmx,1:ndxndxi), start=(/ 1, 1, itim /), count=(/ kmx+1, ndxndxi, 1 /))
          endif

          if (jamapu1>0) then
             do LL=1,lnx
                work1(:, LL) = dmiss ! For proper fill values in z-model runs.
                call getLbotLtopmax(LL,Lb,Ltx)
                call getlayerindicesLmax(LL, nlaybL, nrlayLx)
                do L = Lb,Ltx
                    work1(L-Lb+nlaybL,LL) = u1(L)
                enddo
             enddo
             ierr = nf90_put_var(imapfile, id_unorm(iid), work1(1:kmx,1:lnx), start=(/ 1, 1, itim /), count=(/ kmx, lnx, 1 /))
          endif

          if (jamapu0>0) then
             do LL=1,lnx
                work1(:, LL) = dmiss ! For proper fill values in z-model runs.
                call getLbotLtopmax(LL,Lb,Ltx)
                call getlayerindicesLmax(LL, nlaybL, nrlayLx)
                do L = Lb,Ltx
                    work1(L-Lb+nlaybL,LL) = u0(L)
                enddo
             enddo
             ierr = nf90_put_var(imapfile, id_u0(iid)   , work1(1:kmx,1:lnx), start=(/ 1, 1, itim /), count=(/ kmx, lnx, 1 /))
          endif

          if (jamapq1 > 0) then
             do LL=1,lnx
                work1(:, LL) = dmiss ! For proper fill values in z-model runs.
                call getLbotLtopmax(LL,Lb,Ltx)
                call getlayerindicesLmax(LL, nlaybL, nrlayLx)
                do L = Lb,Ltx
                    work1(L-Lb+nlaybL,LL) = q1(L)
                enddo
             enddo
             ierr = nf90_put_var(imapfile, id_q1(iid)   , work1(1:kmx,1:lnx), start=(/ 1, 1, itim /), count=(/ kmx, lnx, 1 /))
          endif

          if (jamapviu>0) then
             do LL=1,lnx
                work1(:, LL) = dmiss ! For proper fill values in z-model runs.
                call getLbotLtopmax(LL,Lb,Ltx)
                call getlayerindicesLmax(LL, nlaybL, nrlayLx)
                if (javiusp == 1) then       ! user specified part
                    vicc = viusp(LL)
                else
                    vicc = vicouv
                endif
                do L = Lb,Ltx
                    work1(L-Lb+nlaybL,LL) = viu(L) + vicc
                enddo
             enddo
             ierr = nf90_put_var(imapfile, id_viu(iid)   , work1(1:kmx,1:lnx), start=(/ 1, 1, itim /), count=(/ kmx, lnx, 1 /))
          endif

          if (jamapdiu>0) then
             do LL=1,lnx
                work1(:, LL) = dmiss ! For proper fill values in z-model runs.
                call getLbotLtopmax(LL,Lb,Ltx)
                call getlayerindicesLmax(LL, nlaybL, nrlayLx)
                if (jadiusp == 1) then
                    dicc = diusp(LL)
                else
                    dicc = dicouv
                endif
                do L = Lb,Ltx
                    work1(L-Lb+nlaybL,LL) = viu(L) / 0.7 + dicc
                enddo
             enddo
             ierr = nf90_put_var(imapfile, id_diu(iid)   , work1(1:kmx,1:lnx), start=(/ 1, 1, itim /), count=(/ kmx, lnx, 1 /))
          endif

          if (jamaprho>0) then
             do kk=1,ndxndxi
                work1(:, kk) = dmiss ! For proper fill values in z-model runs.
                call getkbotktop(kk,kb,kt)
                call getlayerindices(kk, nlayb, nrlay)
                do k = kb,kt
                   work1(k-kb+nlayb, kk) = rho(k)
                enddo
             enddo
             if ( density_is_pressure_dependent() ) then
                 ierr = nf90_put_var(imapfile, id_rho(iid),  work1(1:kmx,1:ndxndxi), start=(/ 1, 1, itim /), count=(/ kmx, ndxndxi, 1 /))
             else
                 ierr = nf90_put_var(imapfile, id_rhop(iid), work1(1:kmx,1:ndxndxi), start=(/ 1, 1, itim /), count=(/ kmx, ndxndxi, 1 /))
             endif
          endif

          if (jamaptur > 0 .and. iturbulencemodel >= 3) then
             do LL=1,lnx
                work0(:, LL) = dmiss ! For proper fill values in z-model runs.
                call getLbotLtopmax(LL,Lb,Ltx)
                call getlayerindicesLmax(LL, nlaybL, nrlayLx)
                do L = Lb-1,Ltx
                   work0(L-Lb+nlaybL,LL) = turkin1(L)
                enddo
             enddo
             ierr = nf90_put_var(imapfile, id_turkin1(iid)   , work0(0:kmx,1:lnx), start=(/ 1, 1, itim /), count=(/ kmx+1, lnx, 1 /))
             do LL=1,lnx
                work0(:, LL) = dmiss ! For proper fill values in z-model runs.
                call getLbotLtopmax(LL,Lb,Ltx)
                call getlayerindicesLmax(LL, nlaybL, nrlayLx)
                do L = Lb-1,Ltx
                   work0(L-Lb+nlaybL,LL) = tureps1(L)
                enddo
             enddo
             ierr = nf90_put_var(imapfile, id_tureps1(iid)   , work0(0:kmx,1:lnx), start=(/ 1, 1, itim /), count=(/ kmx+1, lnx, 1 /))
             do LL=1,lnx
                work0(:, LL) = dmiss ! For proper fill values in z-model runs.
                call getLbotLtopmax(LL,Lb,Ltx)
                call getlayerindicesLmax(LL, nlaybL, nrlayLx)
                do L = Lb-1,Ltx
                   work0(L-Lb+nlaybL,LL) = vicwwu(L)
                enddo
             enddo
             ierr = nf90_put_var(imapfile, id_vicwwu(iid)   , work0(0:kmx,1:lnx), start=(/ 1, 1, itim /), count=(/ kmx+1, lnx, 1 /))
             work0 = dmiss
             do kk = 1,ndxi
             call getkbotktop(kk,kb,kt)
             call getlayerindices(kk, nlayb, nrlay)
             do k = kb-1,kt
                 work0(k-kb+nlayb, kk) = vicwws(k)
             enddo
             enddo
             ierr = nf90_put_var(imapfile, id_vicwws(iid), work0(0:kmx,1:ndxi), (/ 1, 1, itim /), (/ kmx+1, ndxi, 1 /))
          endif

       endif

       if ( jasecflow > 0 .and. jamapspir > 0) then
          ierr = nf90_put_var(imapfile, id_spirint(iid), spirint, (/ 1, itim /), (/ ndxndxi, 1 /))
          if ( kmx == 0 ) then
             ierr = nf90_put_var(imapfile, id_spircrv(iid), spircrv, (/ 1, itim /), (/ ndxndxi, 1 /))
          endif
       endif

       if ( kmx == 0 ) then
          if (jamapu1>0) then
             ierr = nf90_put_var(imapfile, id_unorm(iid), u1 ,  (/ 1, itim /), (/ lnx , 1 /))
          endif

          if (jamapu0>0) then
             ierr = nf90_put_var(imapfile, id_u0   (iid), u0 ,  (/ 1, itim /), (/ lnx , 1 /))
          endif

          if (jamapq1>0) then
             ierr = nf90_put_var(imapfile, id_q1 (iid)  , q1     , (/ 1, itim /), (/ lnx    , 1 /))
          endif

          if (jamapviu>0) then
             do LL=1,lnx
                work1(:,LL) = dmiss
                if (javiusp == 1) then       ! user specified part
                   vicc = viusp(LL)
                else
                   vicc = vicouv
                endif
                work1(1,LL) = viu(LL) + vicc
             enddo
             ierr = nf90_put_var(imapfile, id_viu (iid), work1(1:1,1:lnx) ,  (/ 1, itim /), (/ lnx , 1 /))
          endif

          if (jamapdiu>0) then
             do LL=1,lnx
                work1(:,LL) = dmiss
                if (jadiusp == 1) then
                   dicc = diusp(LL)
                else
                   dicc = dicouv
                endif
                work1(1,LL) = viu(LL) * 0.7 + dicc
             enddo
             ierr = nf90_put_var(imapfile, id_diu (iid), work1(1:1,1:lnx) ,  (/ 1, itim /), (/ lnx , 1 /))
           endif
       endif

    endif

    if (jaseparate_ /= 2) then


       ! Salinity
       if (jamapsal > 0 .and. jasal > 0) then
          if ( kmx>0 ) then
!            3D
             !do kk=1,ndxndxi
             !   call getkbotktop(kk,kb,kt)
             !   ierr = nf90_put_var(imapfile, id_sa1(iid), sa1(kb:kt), (/ 1, kk, itim /), (/ kt-kb+1, 1, 1 /))
             !enddo
             do kk=1,ndxndxi
                 work1(:,kk) = dmiss ! For proper fill values in z-model runs.
                call getkbotktop(kk,kb,kt)
                call getlayerindices(kk, nlayb, nrlay)
                do k = kb,kt
                   work1(k-kb+nlayb, kk) = constituents(isalt, k)
                enddo
             enddo
             ierr = nf90_put_var(imapfile, id_sa1(iid), work1(1:kmx,1:ndxndxi), (/ 1, 1, itim /), (/ kmx, ndxndxi, 1 /))
          else
             do k = 1, ndxndxi
                sa1(k) = constituents(isalt, k)
            enddo
             ierr = nf90_put_var(imapfile, id_sa1(iid), sa1, (/ 1, itim /), (/ ndxndxi, 1 /))
          endif
       endif

       if (jamaptem > 0 .and. jatem > 0) then
          if ( kmx>0 ) then ! 3D
             !do kk=1,ndxndxi
             !   call getkbotktop(kk,kb,kt)
             !   ierr = nf90_put_var(imapfile, id_tem1(iid), tem1(kb:kt), (/ 1, kk, itim /), (/ kt-kb+1, 1, 1 /))
             !enddo
             do kk=1,ndxndxi
                work1(:,kk) = dmiss ! For proper fill values in z-model runs.
                call getkbotktop(kk,kb,kt)
                call getlayerindices(kk, nlayb, nrlay)
                do k = kb,kt
                   work1(k-kb+nlayb, kk) = constituents(itemp,k)
                enddo
             enddo
             ierr = nf90_put_var(imapfile, id_tem1(iid), work1(1:kmx,1:ndxndxi), (/ 1, 1, itim /), (/ kmx, ndxndxi, 1 /))
          else
             do k = 1, ndxndxi
                tem1(k) = constituents(itemp, k)
             enddo
             ierr = nf90_put_var(imapfile, id_tem1(iid), tem1, (/ 1, itim /), (/ ndxndxi, 1 /))
          endif
       endif

!      tracers
       if (jamapconst > 0 .and. ITRA1 > 0) then ! Note: numtracers is only counting tracer boundaries. SPvdP: now also includes tracers with initial conditions only
          allocate(dum(NdxNdxi))

          do j=ITRA1,ITRAN
             if ( kmx>0 ) then
!               3D
                do kk=1,ndxndxi
                   work1(:, kk) = dmiss ! For proper fill values in z-model runs.
                   call getkbotktop(kk,kb,kt)
                   call getlayerindices(kk, nlayb, nrlay)
                   do k = kb,kt
                      work1(k-kb+nlayb, kk) = constituents(j,k)
                   enddo
                enddo
                ierr = nf90_put_var(imapfile, id_const(iid,j), work1(1:kmx,1:ndxndxi), (/ 1, 1, itim /), (/ kmx, ndxndxi, 1 /))
                !   if ( ierr.ne.0 ) exit  ! probably newly added tracer in the GUI
             else
                do kk=1,NdxNdxi
                   dum(kk) = constituents(j,kk)
                enddo
                ierr = nf90_put_var(imapfile, id_const(iid,j), dum, (/ 1, itim /), (/ NdxNdxi, 1 /) )
             endif
          enddo

          if ( allocated(dum) ) deallocate(dum)
       endif

       ! water quality bottom variables outputs
       if (numwqbots > 0) then
          allocate(dum(NdxNdxi))
          do j=1,numwqbots
             do kk=1,NdxNdxi
                call getkbotktop(kk,kb,kt)
                dum(kk) = wqbot(j,kb)
             enddo
             ierr = nf90_put_var(imapfile, id_wqb(iid,j), dum, (/ 1, itim /), (/ NdxNdxi, 1 /) )
          enddo
          if (wqbot3D_output == 1) then
             do j=1,numwqbots
                do kk=1,ndxndxi
                   work1(:, kk) = dmiss ! For proper fill values in z-model runs.
                   call getkbotktop(kk,kb,kt)
                   call getlayerindices(kk, nlayb, nrlay)
                   do k = kb,kt
                      work1(k-kb+nlayb, kk) = wqbot(j,k)
                   enddo
                enddo
                ierr = nf90_put_var(imapfile, id_wqb3d(iid,j), work1(1:kmx,1:ndxndxi), (/ 1, 1, itim /), (/ kmx, ndxndxi, 1 /))
             enddo
          endif
          if ( allocated(dum) ) deallocate(dum)
       endif

       ! WAQ extra outputs
       if (jawaqproc > 0) then
          do j=1,noout_map
             if (outvar(j)>0)then
                work1 = DMISS ! For proper fill values in z-model runs.
                if ( kmx>0 ) then
!                  3D
                   do kk=1,ndxndxi
                      work1(:, kk) = dmiss ! For proper fill values in z-model runs.
                      call getkbotktop(kk,kb,kt)
                      call getlayerindices(kk, nlayb, nrlay)
                      do k = kb,kt
                         work1(k-kb+nlayb, kk) = waqoutputs(j,k-kbx+1)
                      enddo
                   enddo
                   ierr = nf90_put_var(imapfile, id_waq(iid,j), work1(1:kmx,1:ndxndxi), (/ 1, 1, itim /), (/ kmx, ndxndxi, 1 /))
                else
                   call realloc(dum,NdxNdxi, keepExisting=.false.)
                   do kk=1,NdxNdxi
                      dum(kk) = waqoutputs(j,kk)
                   enddo
                   ierr = nf90_put_var(imapfile, id_waq(iid,j), dum, (/ 1, itim /), (/ NdxNdxi, 1 /) )
                   if (allocated(dum)) then
                      deallocate(dum)
                   endif
                endif
             endif
          enddo
          do j=1,noout_statt
             jj = noout_user + j
             if (outvar(jj)>0)then
                work1 = DMISS ! For proper fill values in z-model runs.
                if ( kmx>0 ) then
!                  3D
                   do kk=1,ndxndxi
                      work1(:, kk) = dmiss ! For proper fill values in z-model runs.
                      call getkbotktop(kk,kb,kt)
                      call getlayerindices(kk, nlayb, nrlay)
                      do k = kb,kt
                         work1(k-kb+nlayb, kk) = waqoutputs(jj,k-kbx+1)
                      enddo
                   enddo
                   ierr = nf90_put_var(imapfile, id_wqst(iid,j), work1(1:kmx,1:ndxndxi), (/ 1, 1, itim /), (/ kmx, ndxndxi, 1 /))
                else
                   call realloc(dum,NdxNdxi, keepExisting=.false.)
                   do kk=1,NdxNdxi
                      dum(kk) = waqoutputs(jj,kk)
                   enddo
                   ierr = nf90_put_var(imapfile, id_wqst(iid,j), dum, (/ 1, itim /), (/ NdxNdxi, 1 /) )
                   if (allocated(dum)) then 
                      deallocate(dum)
                   endif
                endif
             endif
          enddo
          if (comparereal(tim, ti_mape, eps10) == 0) then
             do j=1,noout_state
                jj = noout_user + noout_statt + j
                if (outvar(jj)>0)then
                   work1 = DMISS ! For proper fill values in z-model runs.
                   if ( kmx>0 ) then
!                     3D
                      do kk=1,ndxndxi
                         work1(:, kk) = dmiss ! For proper fill values in z-model runs.
                         call getkbotktop(kk,kb,kt)
                         call getlayerindices(kk, nlayb, nrlay)
                         do k = kb,kt
                            work1(k-kb+nlayb, kk) = waqoutputs(jj,k-kbx+1)
                         enddo
                      enddo
                      ierr = nf90_put_var(imapfile, id_wqse(iid,j), work1(1:kmx,1:ndxndxi), (/ 1, 1 /), (/ kmx, ndxndxi, 1 /))
                   else
                      call realloc(dum,NdxNdxi, keepExisting=.false.)
                      do kk=1,NdxNdxi
                         dum(kk) = waqoutputs(jj,kk)
                      enddo
                      ierr = nf90_put_var(imapfile, id_wqse(iid,j), dum, (/ 1 /), (/ NdxNdxi, 1 /) )
                      if (allocated(dum)) then 
                         deallocate(dum)
                      endif
                   endif
                endif
             enddo
          endif
       endif

       if (jased>0 .and. stm_included) then
          if (stmpar%lsedsus > 0) then
             if (kmx > 0) then
                do kk = 1, ndxndxi
                   call getkbotktop(kk, kb, kt)
                   ierr = nf90_put_var(imapfile, id_ws(iid), mtd%ws(kb:kt,:), (/ 1, kk , 1 , itim /), (/ kt-kb+1, 1 , stmpar%lsedsus , 1 /))
                enddo
             else
                ierr = nf90_put_var(imapfile, id_ws(iid),mtd%ws, (/ 1 , 1 , itim /), (/ ndxndxi , stmpar%lsedsus , 1 /))
             endif
             !
             ! equilibrium concentration, 2D only
             if (kmx == 0) then
                ierr = nf90_put_var(imapfile, id_rsedeq(iid), sedtra%rsedeq(1:ndxndxi, :), (/ 1 , 1 , itim /), (/ ndxndxi , stmpar%lsedsus , 1 /))
             endif

             if (stmpar%morpar%moroutput%sourcesink) then
                ierr = nf90_put_var(imapfile, id_sourse(iid) , sedtra%sourse(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedsus, 1 /))
                ierr = nf90_put_var(imapfile, id_sinkse(iid) , sedtra%sinkse(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedsus, 1 /))
             endif

             if (stmpar%morpar%moroutput%suvcor) then
                ierr = nf90_put_var(imapfile, id_scrn(iid) , sedtra%e_scrn(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedsus, 1 /))
                !ierr = nf90_put_var(imapfile, id_scrt(iid) , sedtra%e_scrt(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedsus, 1 /))
             endif

             if (stmpar%morpar%moroutput%aks) then
                ierr = nf90_put_var(imapfile, id_aks(iid) , sedtra%aks(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedsus, 1 /))
                ierr = nf90_put_var(imapfile, id_rca(iid) , sedtra%rca(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedsus, 1 /))
             endif
             !
             ! Suspended fractions
             call realloc(dum,NdxNdxi, keepExisting=.false.)
             do j=ISED1,ISEDN
                if ( kmx>0 ) then
      !            3D
                   do kk=1,ndxndxi
                      call getkbotktop(kk,kb,kt)
                      do k = kb,kt
                         ! TODO: UNST-976, incorrect for Z-layers:
                         work1(k-kb+1,kk) = constituents(j,k)
                      enddo
                   enddo
                   ierr = nf90_put_var(imapfile, id_const(iid,j), work1(1:kmx,1:ndxndxi), (/ 1, 1, itim /), (/ kmx, ndxndxi, 1 /))
                else
                   do kk=1,NdxNdxi
                      dum(kk) = constituents(j,kk)
                   enddo
                   ierr = nf90_put_var(imapfile, id_const(iid,j), dum, (/ 1, itim /), (/ NdxNdxi, 1 /) )
                endif
             enddo
             if ( allocated(dum) ) deallocate(dum)
          endif

          if (stmpar%morpar%moroutput%dzduuvv) then ! bedslope
             ierr = nf90_put_var(imapfile, id_dzdn(iid), sedtra%e_dzdn, (/ 1, itim /), (/ lnxi , 1 /))
             ierr = nf90_put_var(imapfile, id_dzdt(iid), sedtra%e_dzdt, (/ 1, itim /), (/ lnxi , 1 /))
          endif

          if (stmpar%morpar%moroutput%umod) then
             ierr = nf90_put_var(imapfile, id_umod(iid), sedtra%umod, (/ 1, itim /), (/ ndxndxi, 1 /))
          endif

          if (stmpar%morpar%moroutput%zumod) then
             ierr = nf90_put_var(imapfile, id_zumod(iid), sedtra%zumod, (/ 1, itim /), (/ ndxndxi, 1 /))
          endif

          if (stmpar%morpar%moroutput%ustar) then
             ierr = nf90_put_var(imapfile, id_ustar(iid), sqrt(sedtra%ust2), (/ 1, itim /), (/ ndxndxi, 1 /))
          endif

          if (stmpar%morpar%moroutput%sbcuv) then
             do l = 1, stmpar%lsedtot
                select case(stmpar%morpar%moroutput%transptype)
                case (0)
                   rhol = 1d0
                case (1)
                   rhol = stmpar%sedpar%cdryb(l)
                case (2)
                   rhol = stmpar%sedpar%rhosol(l)
                end select
                sedtra%sbcx(:,l) = sedtra%sbcx(:,l)/rhol
                sedtra%sbcy(:,l) = sedtra%sbcy(:,l)/rhol
             enddo
             ierr = nf90_put_var(imapfile, id_sbcx(iid) , sedtra%sbcx(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
             ierr = nf90_put_var(imapfile, id_sbcy(iid) , sedtra%sbcy(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
          endif

          if (stmpar%morpar%moroutput%sbwuv) then
             do l = 1, stmpar%lsedtot
                select case(stmpar%morpar%moroutput%transptype)
                case (0)
                   rhol = 1d0
                case (1)
                   rhol = stmpar%sedpar%cdryb(l)
                case (2)
                   rhol = stmpar%sedpar%rhosol(l)
                end select
                sedtra%sbwx(:,l) = sedtra%sbwx(:,l)/rhol
                sedtra%sbwy(:,l) = sedtra%sbwy(:,l)/rhol
             enddo
             ierr = nf90_put_var(imapfile, id_sbwx(iid) , sedtra%sbwx(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
             ierr = nf90_put_var(imapfile, id_sbwy(iid) , sedtra%sbwy(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
          endif

          if (stmpar%morpar%moroutput%sswuv) then
             do l = 1, stmpar%lsedtot
                select case(stmpar%morpar%moroutput%transptype)
                case (0)
                   rhol = 1d0
                case (1)
                   rhol = stmpar%sedpar%cdryb(l)
                case (2)
                   rhol = stmpar%sedpar%rhosol(l)
                end select
                sedtra%sswx(:,l) = sedtra%sswx(:,l)/rhol
                sedtra%sswy(:,l) = sedtra%sswy(:,l)/rhol
             enddo
             ierr = nf90_put_var(imapfile, id_sswx(iid) , sedtra%sswx(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
             ierr = nf90_put_var(imapfile, id_sswy(iid) , sedtra%sswy(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
          endif

          if (stmpar%morpar%moroutput%sscuv) then
             call realloc(toutputx, (/ndx, stmpar%lsedsus /), keepExisting=.false., fill = -999d0)
             call realloc(toutputy, (/ndx, stmpar%lsedsus /), keepExisting=.false., fill = -999d0)
             do l = 1, stmpar%lsedsus
                select case(stmpar%morpar%moroutput%transptype)
                case (0)
                   rhol = 1d0
                case (1)
                   rhol = stmpar%sedpar%cdryb(sedtot2sedsus(sedtot2sedsus(l)))
                case (2)
                   rhol = stmpar%sedpar%rhosol(sedtot2sedsus(sedtot2sedsus(l)))
                end select
                toutputx(:,l) = sedtra%sscx(:,sedtot2sedsus(l))/rhol         ! mapping necessary because dim(sscx)=lsedtot
                toutputy(:,l) = sedtra%sscy(:,sedtot2sedsus(l))/rhol
             enddo
             ierr = nf90_put_var(imapfile, id_sscx(iid) , toutputx(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedsus, 1 /))
             ierr = nf90_put_var(imapfile, id_sscy(iid) , toutputy(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedsus, 1 /))
          endif

          ! Get cell centre transport values, removed from fm_erosed and fm_bott3d, and calculated here and stored in sscx, sscy, sbcx, sbcy, sbwx, sbwy, sswx, sswy
          call reconstructsedtransports()

          if (stmpar%morpar%moroutput%sbcuv) then
             do l = 1, stmpar%lsedtot
                select case(stmpar%morpar%moroutput%transptype)
                case (0)
                   rhol = 1d0
                case (1)
                   rhol = stmpar%sedpar%cdryb(l)
                case (2)
                   rhol = stmpar%sedpar%rhosol(l)
                end select
                sedtra%sbcx(:,l) = sedtra%sbcx(:,l)/rhol
                sedtra%sbcy(:,l) = sedtra%sbcy(:,l)/rhol
             enddo
             ierr = nf90_put_var(imapfile, id_sbcx_reconstructed(iid) , sedtra%sbcx(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
             ierr = nf90_put_var(imapfile, id_sbcy_reconstructed(iid) , sedtra%sbcy(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
          endif

          if (stmpar%morpar%moroutput%sbwuv) then
             do l = 1, stmpar%lsedtot
                select case(stmpar%morpar%moroutput%transptype)
                case (0)
                   rhol = 1d0
                case (1)
                   rhol = stmpar%sedpar%cdryb(l)
                case (2)
                   rhol = stmpar%sedpar%rhosol(l)
                end select
                sedtra%sbwx(:,l) = sedtra%sbwx(:,l)/rhol
                sedtra%sbwy(:,l) = sedtra%sbwy(:,l)/rhol
             enddo
             ierr = nf90_put_var(imapfile, id_sbwx_reconstructed(iid) , sedtra%sbwx(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
             ierr = nf90_put_var(imapfile, id_sbwy_reconstructed(iid) , sedtra%sbwy(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
          endif

          if (stmpar%morpar%moroutput%sswuv) then
             do l = 1, stmpar%lsedtot
                select case(stmpar%morpar%moroutput%transptype)
                case (0)
                   rhol = 1d0
                case (1)
                   rhol = stmpar%sedpar%cdryb(l)
                case (2)
                   rhol = stmpar%sedpar%rhosol(l)
                end select
                sedtra%sswx(:,l) = sedtra%sswx(:,l)/rhol
                sedtra%sswy(:,l) = sedtra%sswy(:,l)/rhol
             enddo
             ierr = nf90_put_var(imapfile, id_sswx_reconstructed(iid) , sedtra%sswx(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
             ierr = nf90_put_var(imapfile, id_sswy_reconstructed(iid) , sedtra%sswy(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
          endif

          if (stmpar%morpar%moroutput%sscuv) then
             call realloc(toutputx, (/ndx, stmpar%lsedsus /), keepExisting=.false., fill = -999d0)
             call realloc(toutputy, (/ndx, stmpar%lsedsus /), keepExisting=.false., fill = -999d0)
             do l = 1, stmpar%lsedsus
                select case(stmpar%morpar%moroutput%transptype)
                case (0)
                   rhol = 1d0
                case (1)
                   rhol = stmpar%sedpar%cdryb(sedtot2sedsus(sedtot2sedsus(l)))
                case (2)
                   rhol = stmpar%sedpar%rhosol(sedtot2sedsus(sedtot2sedsus(l)))
                end select
                toutputx(:,l) = sedtra%sscx(:,sedtot2sedsus(l))/rhol         ! mapping necessary because dim(sscx)=lsedtot
                toutputy(:,l) = sedtra%sscy(:,sedtot2sedsus(l))/rhol
             enddo
             ierr = nf90_put_var(imapfile, id_sscx_reconstructed(iid) , toutputx(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedsus, 1 /))
             ierr = nf90_put_var(imapfile, id_sscy_reconstructed(iid) , toutputy(1:ndxndxi, :),  (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedsus, 1 /))
          endif

          do l = 1, stmpar%lsedtot
             call realloc(toutputx, (/ndx, stmpar%lsedtot /), keepExisting=.false., fill = -999d0)
             call realloc(toutputy, (/ndx, stmpar%lsedtot /), keepExisting=.false., fill = -999d0)
             select case(stmpar%morpar%moroutput%transptype)
             case (0)
                rhol = 1d0
             case (1)
                rhol = stmpar%sedpar%cdryb(l)
             case (2)
                rhol = stmpar%sedpar%rhosol(l)
             end select
             toutputx(:,l) = sedtra%sxtot(:,l)/rhol
             toutputy(:,l) = sedtra%sytot(:,l)/rhol
          enddo
          ierr = nf90_put_var(imapfile, id_sxtot(iid), toutputx(1:ndxndxi, :), (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
          ierr = nf90_put_var(imapfile, id_sytot(iid), toutputy(1:ndxndxi, :), (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))

          if (stmpar%morpar%bedupd) then
             ierr = nf90_put_var(imapfile, id_morbl(iid), bl(1:ndxndxi), (/ 1, itim /), (/ ndxndxi, 1 /))
          endif

          select case (stmpar%morlyr%settings%iunderlyr)
          case (1)
             ierr = nf90_put_var(imapfile, id_bodsed(iid), stmpar%morlyr%state%bodsed(:, 1:ndxndxi), (/ 1, 1, itim /), (/ stmpar%lsedtot, ndxndxi, 1 /))
             ierr = nf90_put_var(imapfile, id_dpsed(iid), stmpar%morlyr%state%dpsed(1:ndxndxi), (/ 1, itim /), (/ ndxndxi, 1 /))
          case (2)
             !
             ! Calculate values for lyrfrac and porosity
             !
             ! lyrfrac
             if (.not. allocated(frac) ) allocate( frac(stmpar%lsedtot,1:stmpar%morlyr%settings%nlyr,1:ndx  ) )
             frac = -999d0
             if (stmpar%morlyr%settings%iporosity==0) then
                dens => stmpar%sedpar%cdryb
             else
                dens => stmpar%sedpar%rhosol
             endif
             do k = 1, stmpar%morlyr%settings%nlyr
                do nm = 1, ndxndxi
                   if (stmpar%morlyr%state%thlyr(k,nm)>0.0_fp) then
                      do l = 1, stmpar%lsedtot
                           frac(l, k, nm) = stmpar%morlyr%state%msed(l, k, nm)/(dens(l)*stmpar%morlyr%state%svfrac(k, nm) * &
                                            stmpar%morlyr%state%thlyr(k, nm))
                      enddo
                   else
                      frac(:, k, nm) = 0d0
                   endif
                enddo
             enddo
             !
             if (stmpar%morlyr%settings%iporosity>0) then
                if (.not. allocated(poros) ) allocate( poros(1:stmpar%morlyr%settings%nlyr, 1:ndx ) )
                poros = 1d0-stmpar%morlyr%state%svfrac
             endif
             !
             ! Avoid stack overflow issues with large models
             do l = 1, stmpar%lsedtot
                ierr = nf90_put_var(imapfile, id_msed(iid), stmpar%morlyr%state%msed(l,:,1:ndxndxi), (/ l, 1, 1, itim /), (/ 1, stmpar%morlyr%settings%nlyr, ndxndxi, 1 /))
                ierr = nf90_put_var(imapfile, id_lyrfrac(iid), frac(l,:,1:ndxndxi), (/ l, 1, 1, itim /), (/ 1, stmpar%morlyr%settings%nlyr, ndxndxi, 1 /))
             enddo
             ierr = nf90_put_var(imapfile, id_thlyr(iid), stmpar%morlyr%state%thlyr(:,1:ndxndxi), (/ 1, 1, itim /), (/ stmpar%morlyr%settings%nlyr, ndxndxi, 1 /))
             if (stmpar%morlyr%settings%iporosity>0) then
                ierr = nf90_put_var(imapfile, id_poros(iid), poros(:,1:ndxndxi), (/ 1, 1, itim /), (/ stmpar%morlyr%settings%nlyr, ndxndxi, 1 /))
             endif
          end select

          if (stmpar%morpar%moroutput%taurat) then
             ierr = nf90_put_var(imapfile, id_taurat(iid), sedtra%taurat(1:ndxndxi,:), (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
          endif
          if (stmpar%morpar%moroutput%dm) then
             ierr = nf90_put_var(imapfile, id_dm(iid), sedtra%dm(1:ndxndxi), (/ 1, itim /), (/ ndxndxi, 1 /))

          endif
          if (stmpar%morpar%moroutput%dg) then
             ierr = nf90_put_var(imapfile, id_dg(iid), sedtra%dg(1:ndxndxi), (/ 1, itim /), (/ ndxndxi, 1 /))
          endif
          if (stmpar%morpar%moroutput%dgsd) then
             ierr = nf90_put_var(imapfile, id_dgsd(iid), sedtra%dgsd(1:ndxndxi), (/ 1, itim /), (/ ndxndxi, 1 /))
          endif
          if (stmpar%morpar%moroutput%percentiles) then    ! JRE to do: check with Arthur
             call realloc(dum,ndxndxi, keepExisting=.false.)
             do l = 1, stmpar%morpar%nxx
                do kk=1,NdxNdxi
                   dum(kk) = sedtra%dxx(kk, l)
                enddo
                ierr = nf90_put_var(imapfile, id_dxx(l,iid), dum, (/ 1, itim /), (/ ndxndxi, 1 /))
             enddo
          endif
          if (stmpar%morpar%moroutput%frac) then
             ierr = nf90_put_var(imapfile, id_frac(iid), sedtra%frac(1:ndxndxi, :), (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
          endif
          if (stmpar%morpar%moroutput%mudfrac) then
             ierr = nf90_put_var(imapfile, id_mudfrac(iid), sedtra%mudfrac(1:ndxndxi), (/ 1, itim /), (/ ndxndxi,  1 /))
          endif
          if (stmpar%morpar%moroutput%sandfrac) then
             ierr = nf90_put_var(imapfile, id_sandfrac(iid), sedtra%sandfrac(1:ndxndxi), (/ 1, itim /), (/ ndxndxi,  1 /))
          endif
          if (stmpar%morpar%moroutput%fixfac) then
             ierr = nf90_put_var(imapfile, id_fixfac(iid), sedtra%fixfac(1:ndxndxi,:), (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
          endif
          if (stmpar%morpar%moroutput%hidexp) then
             ierr = nf90_put_var(imapfile, id_hidexp(iid), sedtra%hidexp(1:ndxndxi,:), (/ 1, 1, itim /), (/ ndxndxi, stmpar%lsedtot, 1 /))
          endif
          ! Fluff layers
          if (stmpar%morpar%flufflyr%iflufflyr>0 .and. stmpar%lsedsus>0) then
             do l = 1, stmpar%lsedsus
                call realloc(toutput, ndx, keepExisting=.false., fill=-999d0)
                toutput = stmpar%morpar%flufflyr%mfluff(l,1:ndx)
                ierr = nf90_put_var(imapfile, id_mfluff(iid), toutput(1:ndxndxi), (/ 1, l, itim /), (/ ndxndxi, 1, 1 /))
             enddo
          endif
       endif ! stm

       ! Bedform pars
       if (bfmpar%lfbedfrmout) then
          if (bfmpar%lfbedfrm) then
             ierr = nf90_put_var(imapfile, id_duneheight(iid), bfmpar%duneheight(1:ndxndxi), (/ 1, itim /), (/ ndxndxi, 1 /))
             ierr = nf90_put_var(imapfile, id_dunelength(iid), bfmpar%dunelength(1:ndxndxi), (/ 1, itim /), (/ ndxndxi, 1 /))
          endif
          if (bfmpar%lfbedfrmrou) then
             ierr = nf90_put_var(imapfile, id_ksr(iid), bfmpar%rksr(1:ndxndxi), (/ 1, itim /), (/ ndxndxi, 1 /))
             ierr = nf90_put_var(imapfile, id_ksmr(iid), bfmpar%rksmr(1:ndxndxi), (/ 1, itim /), (/ ndxndxi, 1 /))
             ierr = nf90_put_var(imapfile, id_ksd(iid), bfmpar%rksd(1:ndxndxi), (/ 1, itim /), (/ ndxndxi, 1 /))

             do k = 1,ndxndxi
                rks(k) = sqrt(bfmpar%rksr(k)**2 + bfmpar%rksmr(k)**2 + bfmpar%rksd(k)**2)
             enddo
             ierr = nf90_put_var(imapfile, id_ks(iid), rks(1:ndxndxi), (/ 1, itim /), (/ ndxndxi, 1 /))
          endif
       endif
       ! Sediment Herman
       if (jased > 0 .and. .not.stm_included) then
          ierr = nf90_put_var(imapfile, id_sed(iid), sed, (/ 1, 1, itim /), (/ mxgr, ndxndxi, 1 /))
          ierr = nf90_put_var(imapfile, id_ero(iid), grainlay, (/ 1, 1, itim /), (/ mxgr, size(grainlay,2) , 1 /))

          ierr = nf90_put_var(imapfile, id_bl(iid), bl, (/ 1, itim /), (/ ndxndxi , 1 /))
          if (jaceneqtr .ne. 1) then
              ierr = nf90_put_var(imapfile, id_zk(iid), zk, (/ 1, itim /), (/ numk , 1 /))
          endif


       ! TODO: AvD: size(grainlay,2) is always correct (mxn), but we have a problem if jaceneqtr==2 and mxn/=numk,
       ! because then the dimension for ero is set to nNetNode, and coordinate attribute refers to NetNode_x
       ! (both length numk), whereas ero itself is shorter than numk.
       endif

       ! 1D2D boundaries
       if (nbnd1d2d > 0) then
          ierr = nf90_put_var(imapfile, id_1d2d_zeta1d(iid),      zbnd1d2d1,  (/ 1, itim /), (/ nbnd1d2d, 1 /))
          ierr = nf90_put_var(imapfile, id_1d2d_crest_level(iid), zcrest1d2d, (/ 1, itim /), (/ nbnd1d2d, 1 /))
          ierr = nf90_put_var(imapfile, id_1d2d_b_2di(iid),       b_2di,      (/ 1, itim /), (/ nbnd1d2d, 1 /))
          ierr = nf90_put_var(imapfile, id_1d2d_b_2dv(iid),       b_2dv,      (/ 1, itim /), (/ nbnd1d2d, 1 /))
          ierr = nf90_put_var(imapfile, id_1d2d_d_2dv(iid),       d_2dv,      (/ 1, itim /), (/ nbnd1d2d, 1 /))
          ierr = nf90_put_var(imapfile, id_1d2d_q_zeta(iid),      qzeta_1d2d, (/ 1, itim /), (/ nbnd1d2d, 1 /))
          ierr = nf90_put_var(imapfile, id_1d2d_q_lat(iid),       qlat_1d2d,  (/ 1, itim /), (/ nbnd1d2d, 1 /))
          ierr = nf90_put_var(imapfile, id_1d2d_cfl(iid),         cfl,        (/ 1, itim /), (/ nbnd1d2d, 1 /))
          ierr = nf90_put_var(imapfile, id_1d2d_sb(iid),          sb_1d2d,    (/ 1, itim /), (/ nbnd1d2d, 1 /))
          ierr = nf90_put_var(imapfile, id_1d2d_s1_2d(iid),       s1_2d,      (/ 1, itim /), (/ nbnd1d2d, 1 /))
          ierr = nf90_put_var(imapfile, id_1d2d_s0_2d(iid),       s0_2d,      (/ 1, itim /), (/ nbnd1d2d, 1 /))
          ierr = nf90_put_var(imapfile, id_1d2d_flow_cond(iid),   FlowCond,   (/ 1, itim /), (/ nbnd1d2d, 1 /))
       endif

       if ( jatidep > 0 .and. jamaptidep.eq.1 ) then
          if ( jaselfal.eq.0 ) then
             do k=1,Ndx
                workx(k) = tidep(1,k)
             enddo
             ierr = nf90_put_var(imapfile, id_tidep(iid), workx,  (/ 1, itim /), (/ ndxndxi, 1 /))
          else ! write potential without SAL and SAL potential
             do k=1,Ndx
                workx(k) = tidep(1,k) - tidep(2,k)
!                worky(k) = tidep(2,k)
             enddo
             ierr = nf90_put_var(imapfile, id_tidep(iid), workx,  (/ 1, itim /), (/ ndxndxi, 1 /))
!             ierr = nf90_put_var(imapfile, id_salp(iid),  worky,  (/ 1, itim /), (/ ndxndxi, 1 /))
          endif
       endif
       if ( jaselfal.gt.0 .and. jamapselfal.eq.1 ) then
          do k=1,Ndx
             worky(k) = tidep(2,k)
          enddo
          ierr = nf90_put_var(imapfile, id_salp(iid),  worky,  (/ 1, itim /), (/ ndxndxi, 1 /))
       endif

       if ( jaFrcInternalTides2D.gt.0 .and. jamapIntTidesDiss.eq.1 ) then
          ierr = nf90_put_var(imapfile, id_inttidesdiss(iid), DissInternalTidesPerArea,  (/ 1, itim /), (/ ndxndxi, 1 /))
       endif
    endif

    if (jawind > 0 .and. ((jamapwind > 0 .and. jawindstressgiven == 0) .or. (jaseparate_==2))) then
       allocate (windx(ndxndxi), windy(ndxndxi), stat=ierr)
       if (ierr /= 0) call aerr( 'windx/windy', ierr, ndxndxi)
       !windx/y is not set to 0.0 for flownodes without links !
       windx = 0.0d0
       windy = 0.0d0
       do n = 1,ndxndxi
          !
          ! Currently, wx/y is defined on the links
          ! TO DO: EC-module should not be asked for wind components on the links but on the cells
          !
          if (nd(n)%lnx > 0) then
             do i = 1,nd(n)%lnx
                windx(n) = windx(n) + wx(iabs(nd(n)%ln(i)))
                windy(n) = windy(n) + wy(iabs(nd(n)%ln(i)))
             enddo
             windx(n) = windx(n) / nd(n)%lnx
             windy(n) = windy(n) / nd(n)%lnx
          else
             j=1
          endif
       enddo
       ierr = nf90_put_var(imapfile, id_windx  (iid), windx,  (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_windy  (iid), windy,  (/ 1, itim /), (/ ndxndxi, 1 /))
       deallocate (windx, stat=ierr)
       ierr = nf90_put_var(imapfile, id_windxu  (iid), wx,  (/ 1, itim /), (/ lnx, 1 /))
       ierr = nf90_put_var(imapfile, id_windyu  (iid), wy,  (/ 1, itim /), (/ lnx, 1 /))
    endif

    if (jamapwind > 0 .and. japatm > 0) then
       ierr = nf90_put_var(imapfile, id_patm(iid)  , Patm, (/ 1, itim /), (/ ndxndxi, 1 /))
    endif

    if (ice_mapout) then
       ierr = nf90_put_var(imapfile, id_ice_af(iid) , ice_af, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_ice_h(iid)  , ice_h , (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_ice_p(iid)  , ice_p , (/ 1, itim /), (/ ndxndxi, 1 /))
       if (ja_icecover == ICECOVER_SEMTNER) then
          ierr = nf90_put_var(imapfile, id_ice_t(iid)  , ice_t , (/ 1, itim /), (/ ndxndxi, 1 /))
          ierr = nf90_put_var(imapfile, id_snow_h(iid) , snow_h, (/ 1, itim /), (/ ndxndxi, 1 /))
          ierr = nf90_put_var(imapfile, id_snow_t(iid) , snow_t, (/ 1, itim /), (/ ndxndxi, 1 /))
       endif
    endif

    if (jamapheatflux > 0 .and. jatem > 1) then    ! Heat modelling only
       ierr = nf90_put_var(imapfile, id_tair(iid)  , Tair, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_rhum(iid)  , Rhum, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_clou(iid)  , Clou, (/ 1, itim /), (/ ndxndxi, 1 /))

        if (jatem == 5) then
           ierr = nf90_put_var(imapfile, id_qsun(iid)  , Qsunmap  , (/ 1, itim /), (/ ndxndxi, 1 /))
           ierr = nf90_put_var(imapfile, id_qeva(iid)  , Qevamap  , (/ 1, itim /), (/ ndxndxi, 1 /))
           ierr = nf90_put_var(imapfile, id_qcon(iid)  , Qconmap  , (/ 1, itim /), (/ ndxndxi, 1 /))
           ierr = nf90_put_var(imapfile, id_qlong(iid) , Qlongmap , (/ 1, itim /), (/ ndxndxi, 1 /))
           ierr = nf90_put_var(imapfile, id_qfreva(iid), Qfrevamap, (/ 1, itim /), (/ ndxndxi, 1 /))
           ierr = nf90_put_var(imapfile, id_qfrcon(iid), Qfrconmap, (/ 1, itim /), (/ ndxndxi, 1 /))
        endif
        ierr = nf90_put_var(imapfile, id_qtot(iid)  , Qtotmap  , (/ 1, itim /), (/ ndxndxi, 1 /))
    endif
    call realloc(numlimdtdbl, ndxndxi, keepExisting=.false.)
    numlimdtdbl = dble(numlimdt) ! To prevent stack overflow. TODO: remove once integer version is available.
    ierr = nf90_put_var(imapfile, id_numlimdt(iid)  , numlimdtdbl, (/ 1, itim /), (/ ndxndxi, 1 /))
    deallocate(numlimdtdbl)

    ! Roughness from trachytopes
    if (jatrt == 1) then
       ierr = nf90_put_var(imapfile, id_cftrt(iid),  cftrt(:,2),   (/ 1, itim /), (/ numl, 1 /))
    endif

    ! Roughness calibration factors
    if (jacali == 1) then
       ierr = nf90_put_var(imapfile, id_cfcl(iid),  cfclval,   (/ 1, itim /), (/ numl, 1 /))
    endif

    ! JRE - XBeach
    if (jawave .eq. 4) then
       ierr = nf90_put_var(imapfile, id_E(iid), E, (/ 1, itim /), (/ ndxndxi, 1 /)) ! direction integrated
       ierr = nf90_put_var(imapfile, id_R(iid), R, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_H(iid), H, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_urmscc(iid), uorb, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_Fxcc(iid), Fx_cc, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_Fycc(iid), Fy_cc, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_D(iid), D, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_DR(iid), DR, (/ 1, itim /), (/ ndxndxi, 1 /))

       ierr = nf90_put_var(imapfile, id_sigmwav(iid), sigmwav, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_cwav(iid), cwav, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_cgwav(iid), cgwav, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_thetamean(iid), 270d0 - thetamean*180d0/pi, (/ 1, itim /), (/ ndxndxi, 1 /))
       !if ( (windmodel .eq. 1) .and. (jawsource .eq. 1) ) then
       !   ierr = nf90_put_var(imapfile, id_SwE(iid), SwE, (/ 1, itim /), (/ ndxndxi, 1 /))
       !   ierr = nf90_put_var(imapfile, id_SwT(iid), SwT, (/ 1, itim /), (/ ndxndxi, 1 /))
       !endif
    endif

!   deallocate
    if ( NUMCONST.gt.0 ) then
       if ( allocated(idum)     ) deallocate(idum)
    endif

    if (jaseparate_==2 .and. javeg > 0) then
       ierr = nf90_put_var(imapfile, id_rnveg(iid), rnveg, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_diaveg(iid), diaveg, (/ 1, itim /), (/ ndxndxi, 1 /))
       ierr = nf90_put_var(imapfile, id_veg_stemheight(iid), stemheight, (/ 1, itim /), (/ ndxndxi, 1 /))
    endif

end subroutine unc_write_map_filepointer

end module m_unc_write_map_cfold