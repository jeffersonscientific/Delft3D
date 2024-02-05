!! \section eqtrander calls eqtran for the sediment transport formula
!! multiple times to return the transport quantities and the derivates with
!! respect to the effective depth-averaged velocity components u and v.

subroutine eqtrander(sig    ,thick     ,kmax      ,ws        ,ltur      , &
                & frac      ,sigmol    ,dicww     ,lundia    ,taucr0    , &
                & rksrs     ,i2d3d     ,lsecfl    ,spirint   ,suspfrac  , &
                & tetacr    ,concin    , &
                & dzduu     ,dzdvv     ,ubot      ,tauadd    ,sus       , &
                & bed       ,susw      ,bedw      ,espir     ,wave      , &
                & scour     ,ubot_from_com        ,camax     ,eps       , &
                & iform     ,npar      ,par       ,numintpar ,numrealpar, &
                & numstrpar ,dllfunc   ,dllhandle ,intpar    ,realpar   , &
                & strpar    , &
!output:
                & aks       ,caks      ,taurat    ,seddif    ,rsedeq    , &
                & kmaxsd    ,conc2d    ,sbcu      ,sbcv      ,sbwu      , &
                & sbwv      ,sswu      ,sswv      ,dss       ,caks_ss3d , &
                & aks_ss3d  ,ust2      ,t_relax   ,error     , &
                & dsbcudu   ,dsbcudv   ,dsbcvdu   ,dsbcvdv   )

!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2021.                                
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
!  $Id: eqtrander.f90 69757 2021-09-15 21:28:51Z jagers $
!  $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/branches/research/Deltares/20210729_changes_implementation_exner/src/utils_gpl/morphology/packages/morphology_kernel/src/eqtrander.f90 $
!!--description-----------------------------------------------------------------
!
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use mathconsts, only: ee_hp
    use morphology_data_module
    !
    implicit none
!
! Call variables
!
    integer(pntrsize)                   , intent(in)    :: dllhandle
    integer                             , intent(in)    :: i2d3d
    integer                             , intent(in)    :: iform
    integer                             , intent(in)    :: kmax     !  Description and declaration in esm_alloc_int.f90
    integer                             , intent(in)    :: lsecfl   !  Description and declaration in esm_alloc_int.f90
    integer                             , intent(in)    :: ltur     !  Description and declaration in esm_alloc_int.f90
    integer                             , intent(in)    :: lundia   !  Description and declaration in inout.igs
    integer                             , intent(in)    :: npar
    integer                             , intent(in)    :: numintpar
    integer                             , intent(in)    :: numrealpar
    integer                             , intent(in)    :: numstrpar
    integer      , dimension(numintpar) , intent(inout) :: intpar
    real(fp)                            , intent(in)    :: bed
    real(fp)                            , intent(in)    :: bedw
    real(fp)                            , intent(in)    :: camax
    real(fp)     , dimension(kmax)      , intent(inout) :: concin
    real(fp)     , dimension(0:kmax)    , intent(in)    :: dicww    !  Description and declaration in esm_alloc_real.f90
    real(fp)                            , intent(in)    :: dzduu     !  Description and declaration in esm_alloc_real.f90
    real(fp)                            , intent(in)    :: dzdvv     !  Description and declaration in esm_alloc_real.f90
    real(fp)                            , intent(in)    :: eps
    real(fp)                            , intent(in)    :: espir
    real(fp)                            , intent(in)    :: frac     !  Description and declaration in esm_alloc_real.f90
    real(fp)     , dimension(npar)      , intent(inout) :: par
    real(fp)                            , intent(in)    :: rksrs    !  Description and declaration in esm_alloc_real.f90
    real(fp)     , dimension(kmax)      , intent(in)    :: sig      !  Description and declaration in esm_alloc_real.f90
    real(fp)                            , intent(in)    :: sigmol   !  Description and declaration in esm_alloc_real.f90
    real(fp)                            , intent(in)    :: spirint  !  Spiral flow intensity
    real(fp)                            , intent(in)    :: sus
    real(fp)                            , intent(in)    :: susw
    real(fp)                            , intent(in)    :: tauadd
    real(fp)                            , intent(in)    :: taucr0
    real(fp)                            , intent(in)    :: tetacr
    real(fp)     , dimension(kmax)      , intent(in)    :: thick    !  Description and declaration in esm_alloc_real.f90
    real(fp)                            , intent(in)    :: ubot     !  Description and declaration in esm_alloc_real.f90
    real(fp)     , dimension(0:kmax)    , intent(in)    :: ws       !  Description and declaration in esm_alloc_real.f90
    real(hp)     , dimension(numrealpar), intent(inout) :: realpar
    logical                             , intent(in)    :: scour
    logical                             , intent(in)    :: suspfrac !  suspended sediment fraction
    logical                             , intent(in)    :: ubot_from_com
    logical                             , intent(in)    :: wave
    character(256)                      , intent(in)    :: dllfunc
    character(256), dimension(numstrpar), intent(inout) :: strpar
!
    logical                         , intent(out)  :: error     !< error flag
    integer                         , intent(out)  :: kmaxsd    !< layer number interacting with the bed (-)
    real(fp)                        , intent(inout):: aks       !< out parameter for Van Rijn, in parameter for others (m)
    real(fp)                        , intent(out)  :: aks_ss3d  !< Van Rijn reference height (m)
    real(fp)                        , intent(out)  :: caks      !< reference concentration (kg m-3)
    real(fp)                        , intent(out)  :: caks_ss3d !< reference concentration (kg m-3)
    real(fp)                        , intent(out)  :: conc2d    !< equilibrium concentration (kg m-3)
    real(fp)                        , intent(out)  :: dss       !< suspended sediment diameter (m)
    real(fp), dimension(kmax)       , intent(out)  :: rsedeq    !< equilibrium concentration (kg m-3)
    real(fp)                        , intent(out)  :: sbcu      !< bed load in u direction due to currents (kg m-1 s-1)
    real(fp)                        , intent(out)  :: sbcv      !< bed load in v direction due to currents (kg m-1 s-1)
    real(fp)                        , intent(out)  :: sbwu      !< bed load in u direction due to waves (kg m-1 s-1)
    real(fp)                        , intent(out)  :: sbwv      !< bed load in v direction due to waves (kg m-1 s-1)
    real(fp), dimension(0:kmax)     , intent(out)  :: seddif    !< vertical profile for the sediment diffusion 
    real(fp)                        , intent(out)  :: sswu      !< suspended load in u direction due to waves (kg m-1 s-1)
    real(fp)                        , intent(out)  :: sswv      !< suspended load in v direction due to waves (kg m-1 s-1)
    real(fp)                        , intent(out)  :: t_relax   !< relaxation time scale (s)
    real(fp)                        , intent(out)  :: taurat    !< 
    real(fp)                        , intent(out)  :: ust2      !< square of shear velocity (m2 s-2)
    real(fp)                        , intent(out)  :: dsbcudu   !< derivative of sbcu with respect to u (kg m-2)
    real(fp)                        , intent(out)  :: dsbcudv   !< derivative of sbcu with respect to v (kg m-2)
    real(fp)                        , intent(out)  :: dsbcvdu   !< derivative of sbcv with respect to u (kg m-2)
    real(fp)                        , intent(out)  :: dsbcvdv   !< derivative of sbcv with respect to v (kg m-2)
!
! Local variables
!
    integer                :: icall    !< loop index for eqtran call
    real(fp)               :: sbcu_du  !< bed load due to currents in u direction when u is increased (kg m-1 s-1)
    real(fp)               :: sbcu_dv  !< bed load due to currents in u direction when v is increased (kg m-1 s-1)
    real(fp)               :: sbcv_du  !< bed load due to currents in v direction when u is increased (kg m-1 s-1)
    real(fp)               :: sbcv_dv  !< bed load due to currents in v direction when v is increased (kg m-1 s-1)
    real(hp)               :: du       !< adjustment of effective depth-averaged velocity in u direction (m s-1)
    real(hp)               :: dv       !< adjustment of effective depth-averaged velocity in v direction (m s-1)
    real(hp)               :: r        !< ratio of modified and original effective depth-averaged velocities (-)
    real(hp)               :: u        !< modified effective depth-averaged velocity in u direction (m s-1)
    real(hp)               :: utot     !< modified effective depth-averaged velocity magnitude (m s-1)
    real(hp)               :: utot0    !< original effective depth-averaged velocity magnitude (m s-1)
    real(hp)               :: v        !< modified effective depth-averaged velocity in v direction (m s-1)
!
! backup arrays for inout arguments
!
    integer       , dimension(numintpar)  :: intpar0
    real(fp)      , dimension(kmax)       :: concin0
    real(fp)      , dimension(npar)       :: par0
    real(hp)      , dimension(numrealpar) :: realpar0
    character(256), dimension(numstrpar)  :: strpar0
    real(fp)                              :: aks0
!
!! executable statements -------------------------------------------------------
!
    ! first: backup all intent(inout) parameters such that we can reset them
    aks0      = aks
    intpar0   = intpar
    concin0   = concin
    par0      = par
    realpar0  = realpar
    strpar0   = strpar
    !
    du = 0.001_hp
    dv = 0.001_hp
    !
    ! call eqtran three times
    ! the first two times it is called with modified velocities
    ! final time it is with the original velocities
    do icall = 1, 5
       if (icall > 1) then
          ! reset intent(inout) arguments
          aks      = aks0
          intpar   = intpar0
          concin   = concin0
          par      = par0
          realpar  = realpar0
          strpar   = strpar0
       endif
       !
       if (icall == 5) then
          ! use all the original values
       else
          ! modify the effective depth-averaged velocities
          if (icall == 1) then
             u = realpar(RP_EFUMN) + du
             v = realpar(RP_EFVMN)
          elseif (icall == 2) then
             u = realpar(RP_EFUMN)
             v = realpar(RP_EFVMN) + dv
          elseif (icall == 3) then
             u = realpar(RP_EFUMN) - du
             v = realpar(RP_EFVMN)
          elseif (icall == 4) then
             u = realpar(RP_EFUMN)
             v = realpar(RP_EFVMN) - dv
          endif
          utot  = sqrt(u**2 + v**2)
          utot0 = realpar(RP_EFVLM)
          !
          realpar(RP_EFUMN) = u
          realpar(RP_EFVMN) = v
          realpar(RP_EFVLM) = utot
          !
          ! modify the characteristic velocities
          if (comparereal(utot0, 0.0_hp) > 0) then
             r = utot/utot0
             realpar(RP_VELCH) = r * realpar(RP_VELCH)
             realpar(RP_UCHAR) = r * realpar(RP_UCHAR)
             realpar(RP_VCHAR) = r * realpar(RP_VCHAR)
             realpar(RP_USTAR) = r * realpar(RP_USTAR)
          else
             realpar(RP_VELCH) = utot
             realpar(RP_UCHAR) = u
             realpar(RP_VCHAR) = v
             realpar(RP_ZVLCH) = realpar(RP_DEPTH)/ee_hp
             realpar(RP_USTAR) = utot * realpar(RP_VNKAR) / log(1.0_hp + realpar(RP_ZVLCH)/realpar(RP_Z0ROU))
          endif
       endif
       !
       call eqtran(sig       ,thick     ,kmax      ,ws        ,ltur      , &
                 & frac      ,sigmol    ,dicww     ,lundia    ,taucr0    , &
                 & rksrs     ,i2d3d     ,lsecfl    ,spirint   ,suspfrac  , &
                 & tetacr    ,concin    , &
                 & dzduu     ,dzdvv     ,ubot      ,tauadd    ,sus       , &
                 & bed       ,susw      ,bedw      ,espir     ,wave      , &
                 & scour     ,ubot_from_com        ,camax     ,eps       , &
                 & iform     ,npar      ,par       ,numintpar ,numrealpar, &
                 & numstrpar ,dllfunc   ,dllhandle ,intpar    ,realpar   , &
                 & strpar    , &
!output:
                 & aks       ,caks      ,taurat    ,seddif    ,rsedeq    , &
                 & kmaxsd    ,conc2d    ,sbcu      ,sbcv      ,sbwu      , &
                 & sbwv      ,sswu      ,sswv      ,dss       ,caks_ss3d , &
                 & aks_ss3d  ,ust2      ,t_relax   ,error     )
       !
       if (icall == 1) then
          ! store the bed load vector due to currents
          sbcu_du = sbcu
          sbcv_du = sbcv
       elseif (icall == 2) then
          ! store the bed load vector due to currents
          sbcu_dv = sbcu
          sbcv_dv = sbcv
       elseif (icall == 3) then
          ! compute the derivatives for u change
          ! from first call
          dsbcudu = 0.5_fp * (sbcu_du - sbcu)/du
          dsbcvdu = 0.5_fp * (sbcv_du - sbcv)/du
       elseif (icall == 4) then
          ! compute the derivatives for v change
          ! from second call
          dsbcudv = 0.5_fp * (sbcu_dv - sbcu)/dv
          dsbcvdv = 0.5_fp * (sbcv_dv - sbcv)/dv
       else
          ! final call; nothing to do
       endif
    enddo ! icall
end subroutine eqtrander
