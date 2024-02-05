subroutine upwbed(su        ,sv        ,suu       ,svv       ,kfu       , &
                & dsbcudu   ,dsbcudv   ,dsbcvdu   ,dsbcvdv   , &
                & kfv       ,kcs       ,kfsed     ,lsedtot   , &
                & nmmax     ,icx       ,icy       ,sutot     ,svtot     , &
                & gvu       ,guv       ,hu        ,hv        ,umean     , &
                & vmean     ,cdryb     ,dtmor     ,ag        ,dps       , &
                & gdp       )
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2024.                                
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!  
!  
!!--description-----------------------------------------------------------------
!
!    Function: Copy transport rate from cell centres to velocity points
!              using an upwind or central approach.
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use sediment_basics_module, only: has_bedload
    use morphology_data_module, only: BL_SCHEME_UPWSB, BL_SCHEME_CENTRAL
    use flux_limiters
    !
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    !
    integer, dimension(:)                , pointer :: tratyp
    real(fp)                             , pointer :: bed
    type (mornumericstype)               , pointer :: mornum
!
! Global variables
!
    integer                                            , intent(in)  :: lsedtot
    integer                                            , intent(in)  :: icx
    integer                                            , intent(in)  :: icy
    integer                                            , intent(in)  :: nmmax
    integer , dimension(gdp%d%nmlb:gdp%d%nmub)         , intent(in)  :: kcs
    integer , dimension(gdp%d%nmlb:gdp%d%nmub)         , intent(in)  :: kfsed
    integer , dimension(gdp%d%nmlb:gdp%d%nmub)         , intent(in)  :: kfu
    integer , dimension(gdp%d%nmlb:gdp%d%nmub)         , intent(in)  :: kfv
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, lsedtot), intent(in)  :: su
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, lsedtot), intent(in)  :: sv
    real(fp), dimension(lsedtot, gdp%d%nmlb:gdp%d%nmub), intent(in)  :: dsbcudu
    real(fp), dimension(lsedtot, gdp%d%nmlb:gdp%d%nmub), intent(in)  :: dsbcudv
    real(fp), dimension(lsedtot, gdp%d%nmlb:gdp%d%nmub), intent(in)  :: dsbcvdu
    real(fp), dimension(lsedtot, gdp%d%nmlb:gdp%d%nmub), intent(in)  :: dsbcvdv
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, lsedtot), intent(in)  :: sutot
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, lsedtot), intent(in)  :: svtot
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, lsedtot), intent(out) :: suu
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, lsedtot), intent(out) :: svv
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)         , intent(in)  :: gvu   !< distance between cell centres in u direction (m)
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)         , intent(in)  :: guv   !< distance between cell centres in v direction (m)
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)         , intent(in)  :: hu    !< water depth at u point (m)
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)         , intent(in)  :: hv    !< water depth at v point (m)
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)         , intent(in)  :: umean !< u component of depth averaged velocity at u point (m/s)
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)         , intent(in)  :: vmean !< v component of depth averaged velocity at v point (m/s)
    real(fp), dimension(lsedtot)                       , intent(in)  :: cdryb !< dry bed density (kg/m3)
    real(fp)                                           , intent(in)  :: dtmor !< morphological time step (s)
    real(fp)                                           , intent(in)  :: ag    !< gravitational acceleration (m/s2)
    real(prec), dimension(gdp%d%nmlb:gdp%d%nmub)         , intent(in)  :: dps   !< bed level depth relative to reference plane bed in cell centres (m)
!
! Local variables
!
    integer  :: l
    integer  :: nm
    integer  :: nmu
    integer  :: nmuu
    integer  :: num
    integer  :: nmd
    integer  :: ndm
    integer  :: ndmu
    integer  :: numd
    integer  :: numu
    integer  :: nuum
    logical  :: laterallyaveragedbedload
    integer  :: bedload_scheme
    real(fp) :: cfl                      !< local Courant number (-)
    real(fp) :: dsbcuduu                 !< dsbcudu in u point (kg/m2)
    real(fp) :: dsbcudvu                 !< dsbcudv in u point (kg/m2)
    real(fp) :: dsbcvduv                 !< dsbcvdu in v point (kg/m2)
    real(fp) :: dsbcvdvv                 !< dsbcvdv in v point (kg/m2)
    real(fp) :: h                        !< local water depth (m)
    real(fp) :: phi                      !< result of limiter function (-)
    real(fp) :: rhoceler                 !< bed celerity multiplied by bed density (kg/m2/s)
    real(fp) :: slpd                     !< bed slope at lower M/N index (-)
    real(fp) :: slpu                     !< bed slope at higher M/N index (-)
    real(fp) :: suv1                     !< transport capacity normal to edge in cel 1 (kg/m/s)
    real(fp) :: suv2                     !< transport capacity normal to edge in cel 2 (kg/m/s)
    real(fp) :: u                        !< u component of the velocity (m/s)
    real(fp) :: um2                      !< velocity magnitude squared (m/s)
    real(fp) :: v                        !< v component of the velocity (m/s)
!
    abstract interface
       function fluxlim (slp1, slp2) result (phi)
          use precision
          real(fp), intent (in) :: slp1
          real(fp), intent (in) :: slp2
          real(fp)              :: phi
       end function fluxlim
    end interface

    procedure (fluxlim), pointer :: fluxlimiter => null ()
!
!! executable statements -------------------------------------------------------
!
    tratyp              => gdp%gdsedpar%tratyp
    bed                 => gdp%gdmorpar%bed
    mornum              => gdp%gdmorpar%mornum
    !
    bedload_scheme           = mornum%bedload_scheme
    laterallyaveragedbedload = mornum%laterallyaveragedbedload
    !
    fluxlimiter => fluxlim_koren
    !
    do l = 1, lsedtot
       ! if the transport of the fraction may include a bedload component
       if (has_bedload(tratyp(l))) then
          do nm = 1, nmmax
             !
             ! set bedload transport at u points
             !
             nmd  = nm - icx
             nmu  = nm + icx
             ndm  = nm - icy
             num  = nm + icy
             !
             nmuu = nm + icx + icx
             nuum = nm + icy + icy
             !
             ndmu = nm - icy + icx
             numd = nm + icy - icx
             numu = nm + icy + icx
             !
             ! if active velocity point with two adjacent active sediment cells
             ! (done to prevent bed-load transport into kfsed=0 cells)
             !
             if ((kfu(nm)*kfsed(nm)*kfsed(nmu)) /= 0) then
                if (laterallyaveragedbedload) then
                   suv1 = (4.0_fp * su(nm, l) &
                         & + kfv(nm) * (su(num, l) - su(nm,l)) &
                         & + kfv(ndm) * (su(ndm, l) - su(nm,l)) &
                         & ) / 4.0_fp
                   suv2 = (4.0_fp*su(nmu,l) &
                         & + kfv(nmu) * (su(numu,l) - su(nmu,l)) &
                         & + kfv(ndmu) * (su(ndmu,l) - su(nmu,l)) &
                         & ) / 4.0_fp
                else
                   suv1 = su(nm, l)
                   suv2 = su(nmu, l)
                endif
                !
                ! At a domain decomposition boundary, always use the internal value.
                ! In the domain decomposition mapper, the fluxes will be made consistent.
                !
                if (kcs(nmu) == 3) then
                   suu(nm, l) = suv1
                elseif (kcs(nm) == 3) then
                   suu(nm, l) = suv2
                !
                ! At other locations use the selected numerical scheme.
                !
                elseif (bedload_scheme == BL_SCHEME_UPWSB) then
                   if (sutot (nm, l)>0.0_fp .and. sutot(nmu,l)>0.0_fp) then
                      suu(nm, l) = suv1
                   elseif (sutot(nm, l)<0.0_fp .and. sutot(nmu, l)<0.0_fp) then
                      suu(nm, l) = suv2
                   else
                      suu(nm, l) = (suv1 + suv2)/2.0_fp
                   endif
                elseif (bedload_scheme == BL_SCHEME_CENTRAL) then
                   suu(nm, l) = (suv1 + suv2)/2.0_fp
                else ! bedload_scheme == BL_SCHEME_UPWIND
                   u   = umean(nm)
                   v   = 0.25_fp * (vmean(nm) + vmean(nmu) + vmean(ndm) + vmean(ndmu))
                   um2 = u**2 + v**2
                   h   = hu(nm)
                   dsbcuduu = 0.5_fp * (dsbcudu(l,nm) + dsbcudu(l,nmu))
                   dsbcudvu = 0.5_fp * (dsbcudv(l,nm) + dsbcudv(l,nmu))
                   rhoceler = ( dsbcuduu * (ag * u * (ag * h - v**2)) &
                       &      + dsbcudvu * (ag * v * (ag * h - u**2)) &
                       &      ) / (ag * h * (ag * h - um2))
                   cfl = (abs(rhoceler) * dtmor) / (cdryb(l) * gvu(nm))
                   if (rhoceler >= 0.0_fp) then
                      if (kcs(nm) == 2) then
                         phi = 0.0_fp
                      else
                         slpd = - real(dps(nm) - dps(nmd), fp) / gvu(nmd)  ! use dzduu instead?
                         slpu = - real(dps(nmu) - dps(nm), fp) / gvu(nm)
                         phi = fluxlimiter(slpd, slpu)
                      endif
                      suu(nm, l) = (suv1 + suv2)/2.0_fp + &
                          & - 0.5_fp * rhoceler * ((1.0_fp - cfl) * phi - 1.0_fp) * real(dps(nmu) - dps(nm), fp)
                   else ! celerity < 0.0_fp
                      if (kcs(nmu) == 2) then
                         phi = 0.0_fp
                      else
                         slpu = - real(dps(nmuu) - dps(nmu), fp) / gvu(nmu)
                         slpd = - real(dps(nmu) - dps(nm), fp) / gvu(nm)
                         phi = fluxlimiter(slpu, slpd)
                      endif
                      suu(nm, l) = (suv1 + suv2)/2.0_fp + &
                          & - 0.5_fp * rhoceler * ((1.0_fp - cfl) * phi - 1.0_fp) * real(dps(nm) - dps(nmu), fp)
                   endif
                endif
             else
                suu(nm, l) = 0.0_fp
             endif
             !
             ! set bedload transport at v points
             !
             ! if active velocity point with two adjacent active sediment cells
             ! (done to prevent bed-load transport into kfsed=0 cells)
             !
             if ((kfv(nm)*kfsed(nm)*kfsed(num)) /= 0) then
                if (laterallyaveragedbedload) then
                   suv1 = (4.0_fp * sv(nm, l) &
                         & + kfu(nm) * (sv(nmu, l) - sv(nm,l)) &
                         & + kfu(nmd) * (sv(nmd, l) - sv(nm,l)) &
                         & ) / 4.0_fp
                   suv2 = (4.0_fp*sv(num,l) &
                         & + kfu(num) * (sv(numu,l) - sv(num,l)) &
                         & + kfu(numd) * (sv(numd,l) - sv(num,l)) &
                         & ) / 4.0_fp
                else
                   suv1 = sv(nm, l)
                   suv2 = sv(num, l)
                endif
                !
                ! At a domain decomposition boundary, always use the internal value.
                ! In the domain decomposition mapper, the fluxes will be made consistent.
                !
                if (kcs(num) == 3) then
                   svv(nm, l) = suv1
                elseif (kcs(nm) == 3) then
                   svv(nm, l) = suv2
                !
                ! At other locations use the selected numerical scheme.
                !
                elseif (bedload_scheme == BL_SCHEME_UPWSB) then
                   if (svtot(nm, l)>0.0_fp .and. svtot(num, l)>0.0_fp) then
                      svv(nm, l) = suv1
                   elseif (svtot(nm, l)<0.0_fp .and. svtot(num, l)<0.0_fp) then
                      svv(nm, l) = suv2
                   else
                      svv(nm, l) = (suv1 + suv2)/2.0_fp
                   endif
                elseif (bedload_scheme == BL_SCHEME_CENTRAL) then
                   svv(nm, l) = (suv1 + suv2)/2.0_fp
                else ! bedload_scheme == BL_SCHEME_UPWIND
                   u   = 0.25_fp * (umean(nm) + umean(num) + umean(nmd) + umean(numd))
                   v   = vmean(nm)
                   um2 = u**2 + v**2
                   h   = hv(nm)
                   dsbcvduv = 0.5_fp * (dsbcvdu(l,nm) + dsbcvdu(l,num))
                   dsbcvdvv = 0.5_fp * (dsbcvdv(l,nm) + dsbcvdv(l,num))
                   rhoceler = ( dsbcvduv * (ag * u * (ag * h - v**2)) &
                       &      + dsbcvdvv * (ag * v * (ag * h - u**2)) &
                       &      ) / (ag * h * (ag * h - um2))
                   cfl = (abs(rhoceler) * dtmor) / (cdryb(l) * guv(nm))
                   if (rhoceler >= 0.0_fp) then
                      if (kcs(nm) == 2) then
                         phi = 0.0_fp
                      else
                         slpd = - real(dps(nm) - dps(ndm), fp) / guv(ndm)  ! use dzdvv instead?
                         slpu = - real(dps(num) - dps(nm), fp) / guv(nm)
                         phi = fluxlimiter(slpd, slpu)
                      endif
                      svv(nm, l) = (suv1 + suv2)/2.0_fp + &
                          & - 0.5_fp * rhoceler * ((1.0_fp - cfl) * phi - 1.0_fp) * real(dps(num) - dps(nm), fp)
                   else ! celerity < 0.0_fp
                      if (kcs(num) == 2) then
                         phi = 0.0_fp
                      else
                         slpu = - real(dps(nuum) - dps(num), fp) / guv(num)
                         slpd = - real(dps(num) - dps(nm), fp) / guv(nm)
                         phi = fluxlimiter(slpu, slpd)
                      endif
                      svv(nm, l) = (suv1 + suv2)/2.0_fp + &
                          & - 0.5_fp * rhoceler * ((1.0_fp - cfl) * phi - 1.0_fp) * real(dps(nm) - dps(num), fp)
                   endif
                endif
             else
                svv(nm, l) = 0.0_fp
             endif
          enddo
       endif
    enddo
end subroutine upwbed
