!!  Copyright (C)  Stichting Deltares, 2012-2023.
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License version 3,
!!  as published by the Free Software Foundation.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program. If not, see <http://www.gnu.org/licenses/>.
!!
!!  contact: delft3d.support@deltares.nl
!!  Stichting Deltares
!!  P.O. Box 177
!!  2600 MH Delft, The Netherlands
!!
!!  All indications and logos of, and references to registered trademarks
!!  of Stichting Deltares remain the property of Stichting Deltares. All
!!  rights reserved.
module m_flocsd
use m_waq_type_definitions


implicit none

contains


      subroutine flocsd     ( pmsa   , fl     , ipoint , increm, noseg , &
                              noflux , iexpnt , iknmrk , noq1  , noq2  , &
                              noq3   , noq4   )
!
!*******************************************************************************
!
      use m_evaluate_waq_attribute
      use flocculation_dwq

      implicit none
!
!     type    name         i/o description
!
      integer(kind=int_32),  parameter :: nopmsa = 17

      real(kind=sp)  ::pmsa(*)        !i/o process manager system array, window of routine to process library
      real(kind=sp)  ::fl(*)          ! o  array of fluxes made by this process in mass/volume/time
      integer(kind=int_32)  ::ipoint(nopmsa) ! i  array of pointers in pmsa to get and store the data
      integer(kind=int_32)  ::increm(nopmsa) ! i  increments in ipoint for segment loop, 0=constant, 1=spatially varying
      integer(kind=int_32)  ::noseg          ! i  number of computational elements in the whole model schematisation
      integer(kind=int_32)  ::noflux         ! i  number of fluxes, increment in the fl array
      integer(kind=int_32)  ::iexpnt(4,*)    ! i  from, to, from-1 and to+1 segment numbers of the exchange surfaces
      integer(kind=int_32)  ::iknmrk(*)      ! i  active-inactive, surface-water-bottom, see manual for use
      integer(kind=int_32)  ::noq1           ! i  nr of exchanges in 1st direction (the horizontal dir if irregular mesh)
      integer(kind=int_32)  ::noq2           ! i  nr of exchanges in 2nd direction, noq1+noq2 gives hor. dir. reg. grid
      integer(kind=int_32)  ::noq3           ! i  nr of exchanges in 3rd direction, vertical direction, pos. downward
      integer(kind=int_32)  ::noq4           ! i  nr of exchanges in the bottom (bottom layers, specialist use only)
      integer(kind=int_32)  ::ipnt(nopmsa)   ! local work array for the pointering
      integer(kind=int_32)  ::iseg           ! local loop counter for computational element loop
!
!*******************************************************************************
!
!     type    name         i/o description                                        unit
!
      integer(kind=int_32)  ::idflocim1
      integer(kind=int_32)  ::idflocim2

      integer(kind=int_32)  ::ip8, in8, ip9, in9, ipwmac, inwmac, ipwmic, inwmic, iq, noq, ivan

      real(kind=sp)  ::cmacro      ! i  inorganic matter (im1; macro flocs)                (gdm/m3)
      real(kind=sp)  ::cmicro      ! i  inorganic matter (im2; micro flocs)                (gdm/m3)
      real(kind=sp)  ::tpm         ! i  total particulate matter (including algae)         (gdw/m3)
      real(kind=sp)  ::tau         ! i  bottom shear stress                                (Pa)
      real(kind=sp)  ::tke         ! i  turbulent kinetic energy                           (J)
      integer(kind=int_32)  ::swfloform   ! i  1=Manning/Dyer, 2=Chassagne/Safar, 3=NA            (-)
      real(kind=sp)  ::rcfloc      ! i  flocculation rate                                  (1/d)
      real(kind=sp)  ::rcbreakup   ! i  floc break-up rate                                 (1/d)
      real(kind=sp)  ::delt        ! i  timestep for processes                             (d)
      real(kind=sp)  ::total_depth ! i  total depth for segment (bottom to surface)        (m)
      real(kind=sp)  ::local_depth ! i  local depth for segment (segment to surface)       (m)
      real(kind=sp)  ::rho_water   ! i  density of water                                   (kg/m3)
      real(kind=sp)  ::viscosity   ! i  kinetic viscosity                                  (m/s2)
      real(kind=sp)  ::spmratioem  ! o  flocculation ratio macro:micro empirical model     (-)
      real(kind=sp)  ::dflocim1    ! f  flocculation or break-up flux im1                  (g/m3/d)
      real(kind=sp)  ::dflocim2    ! f  flocculation or break-up flux im2                  (g/m3/d)
      integer(kind=int_32)  ::ikmrk1      !    first segment attribute
      logical active      !    active segment
      logical bottom      !    sediment bed segment
      real(kind=sp)  ::macro       !    concentration macro flocs                            (g/m3)
      real(kind=sp)  ::micro       !    concentration micro flocs                            (g/m3)
      real(kind=sp)  ::tim         !    total concentration flocs                            (g/m3)
      real(kind=sp)  ::macroeq     !    concentration macro flocs in equilibrium             (g/m3)
      real(kind=sp)  ::dfloc       !    flocculation or break-up flux                      (g/m3/d)
      real(kind=sp)  ::ws_macro    !    fall velocity of macro flocs                       (m/d)
      real(kind=sp)  ::ws_micro    !    fall velocity of micro flocs                       (m/d)
!
!*******************************************************************************
!
      ipnt      = ipoint
      idflocim1 = 1
      idflocim2 = 2

      do iseg = 1 , noseg

         cmacro      = pmsa( ipnt(  1) )    ! IM1
         cmicro      = pmsa( ipnt(  2) )    ! IM2
         tpm         = pmsa( ipnt(  3) )
         tke         = pmsa( ipnt(  4) )
         swfloform   = pmsa( ipnt(  5) )
         rcfloc      = pmsa( ipnt(  6) )
         rcbreakup   = pmsa( ipnt(  7) )
         ws_macro    = pmsa( ipnt(  8) )
         ws_micro    = pmsa( ipnt(  9) )
         rho_water   = pmsa( ipnt( 10) )
         viscosity   = pmsa( ipnt( 11) )
         delt        = pmsa( ipnt( 12) )
         total_depth = pmsa( ipnt( 13) )
         local_depth = pmsa( ipnt( 14) )

         ! only for active water segments

         active = btest(iknmrk(iseg),0)
         call evaluate_waq_attribute(1,iknmrk(iseg),ikmrk1)
         bottom = ikmrk1.eq.3
         if ( active .and. .not. bottom ) then
            tke = tau / param_soulsby ! Very coarse estimate!
            call flocculate_dwq( swfloform, cmacro, cmicro, tpm, tke, tau, total_depth, local_depth, viscosity, rho_water, &
                                 spmratioem, ws_macro, ws_micro )


            ! calculate flocculatio/break-up flux and restrict flux to 50% in one timestep for stability

            tim     = cmacro+cmicro
            macroeq = spmratioem*tim/(1.+spmratioem)
            if ( macroeq .gt. cmacro ) then
               dfloc = (macroeq-cmacro)*rcfloc
               dfloc = min(dfloc,0.5*cmicro/delt)
            else
               dfloc = (macroeq-cmacro)*rcbreakup
               dfloc = max(dfloc,-0.5*cmacro/delt)
            endif

         else

            spmratioem = -999.
            dfloc      = 0.0

         endif

         ! pass values back to the system

         fl  ( idflocim1   ) =  dfloc
         fl  ( idflocim2   ) = -dfloc
         pmsa( ipnt( 15)   ) =  spmratioem
         pmsa( ipnt(  8)   ) =  ws_macro
         pmsa( ipnt(  9)   ) =  ws_micro

         idflocim1   = idflocim1   + noflux
         idflocim2   = idflocim2   + noflux
         ipnt        = ipnt        + increm

      enddo

      !
      ! Now fill in the fall velocities
      !
      noq = noq1 + noq2 + noq3

      ipwmac = ipoint(16)
      inwmac = increm(16)
      ipwmic = ipoint(17)
      inwmic = increm(17)

      !
      ! Horizontal exchanges - set to zero
      !
      do iq = 1,noq1+noq2

         pmsa(ipwmac) = 0.0
         pmsa(ipwmic) = 0.0

         ipwmac = ipwmac + inwmac
         ipwmic = ipwmic + inwmic

      enddo

      do iq=noq1+noq2+1,noq

         ivan = iexpnt(1,iq)
!
!        sedimentation velocity from segment to exchange-area
!
         if ( ivan .gt. 0 ) then
            ip8 = ipoint(8) + (ivan-1) * in8
            ip9 = ipoint(9) + (ivan-1) * in9
            pmsa(ipwmac) = pmsa( ip8 )           ! Correct entries? I am in doubt because of the workarray aspect
            pmsa(ipwmic) = pmsa( ip9 )
         endif

         ipwmac = ipwmac + inwmac
         ipwmic = ipwmic + inwmic

      enddo

      end subroutine flocsd

end module m_flocsd
