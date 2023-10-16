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
module m_macnut
use m_waq_type_definitions


implicit none

contains


      subroutine MACNUT     ( pmsa   , fl     , ipoint , increm, noseg , &
                              noflux , iexpnt , iknmrk , noq1  , noq2  , &
                              noq3   , noq4   )
      use m_evaluate_waq_attribute

!
!*******************************************************************************
!
      implicit none
!
!     type    name         i/o description
!
      real(kind=sp)  ::pmsa(*)     !i/o process manager system array, window of routine to process library
      real(kind=sp)  ::fl(*)       ! o  array of fluxes made by this process in mass/volume/time
      integer(kind=int_32)  ::ipoint( 45) ! i  array of pointers in pmsa to get and store the data
      integer(kind=int_32)  ::increm( 45) ! i  increments in ipoint for segment loop, 0=constant, 1=spatially varying
      integer(kind=int_32)  ::noseg       ! i  number of computational elements in the whole model schematisation
      integer(kind=int_32)  ::noflux      ! i  number of fluxes, increment in the fl array
      integer(kind=int_32)  ::iexpnt(4,*) ! i  from, to, from-1 and to+1 segment numbers of the exchange surfaces
      integer(kind=int_32)  ::iknmrk(*)   ! i  active-inactive, surface-water-bottom, see manual for use
      integer(kind=int_32)  ::noq1        ! i  nr of exchanges in 1st direction (the horizontal dir if irregular mesh)
      integer(kind=int_32)  ::noq2        ! i  nr of exchanges in 2nd direction, noq1+noq2 gives hor. dir. reg. grid
      integer(kind=int_32)  ::noq3        ! i  nr of exchanges in 3rd direction, vertical direction, pos. downward
      integer(kind=int_32)  ::noq4        ! i  nr of exchanges in the bottom (bottom layers, specialist use only)
      integer(kind=int_32)  ::ipnt( 45)   !    local work array for the pointering
      integer(kind=int_32)  ::iseg        !    local loop counter for computational element loop
!
!*******************************************************************************
!
!     type    name         i/o description                                        unit
!
      real(kind=sp)  ::depth       ! i  depth of segment                                   (m)
      real(kind=sp)  ::totaldepth  ! i  total depth water column                           (m)
      real(kind=sp)  ::locseddept  ! i  sediment layer depth to bottom of segment          (m)
      real(kind=sp)  ::po4s12      ! i  concentration of PO4 in S12 bottom (pores)         (g/m3)
      real(kind=sp)  ::nh4s12      ! i  concentration of NH4 in S12 bottom (pores)         (g/m3)
      integer(kind=int_32)  ::ibotseg     ! i  bottom segment number                              (-)
      real(kind=sp)  ::frbmlay     ! i  fraction in layer of SM biomass in column          (-)
      real(kind=sp)  ::rootdesm01  ! i  rooting depth sm01                                 (m)
      real(kind=sp)  ::poros       ! i  volumetric porosity                                (-)
      real(kind=sp)  ::nh4         ! i  ammonium (nh4)                                     (gn/m3)
      real(kind=sp)  ::no3         ! i  nitrate (no3)                                      (gn/m3)
      real(kind=sp)  ::po4         ! i  ortho-phosphate (po4)                              (gp/m3)
      real(kind=sp)  ::disco2      ! i  concentration of dissolved carbon dioxide          (g/m3)
      real(kind=sp)  ::dish2co3    ! i  concentration of dissolved true h2co3              (gc/m3)
      real(kind=sp)  ::dishco3     ! i  concentration of dissolved hco3(-)                 (gc/m3)
      real(kind=sp)  ::prfnh4sm01  ! i  ammonium preferency over nitrate sm01              (-)
      real(kind=sp)  ::kmdinsm01w  ! i  half-saturation value n sm01 in water              (gn/m3)
      real(kind=sp)  ::kmpsm01w    ! i  half-saturation value p sm01 in water              (gp/m3)
      real(kind=sp)  ::kmco2sm01   ! i  half-saturation value co2 + h2c03 sm01             (gc/m3)
      real(kind=sp)  ::kmhco3sm01  ! i  half-saturation value hc03 sm01                    (gc/m3)
      real(kind=sp)  ::kmdinsm01b  ! i  half-saturation value n sm01 in bottom             (gn/m3)
      real(kind=sp)  ::kmpsm01b    ! i  half-saturation value p sm01 in bottom             (gp/m3)
      real(kind=sp)  ::temp        ! i  ambient water temperature                         (oC)
      real(kind=sp)  ::limnsm01w   ! o  nitrogen limitation function sm01 <0-1> water      (-)
      real(kind=sp)  ::limpsm01w   ! o  phosphorus limitation function sm01 <0-1> water    (-)
      real(kind=sp)  ::lco2sm01    ! o  co2+h2c03 limitation function sm01 <0-1>           (-)
      real(kind=sp)  ::limnsm01b   ! o  nitrogen limitation function sm01 <0-1> bottom     (-)
      real(kind=sp)  ::limpsm01b   ! o  phosphorus limitation function sm01 <0-1>bottom    (-)
      real(kind=sp)  ::limnutsm01  ! o  nutrient limitation function sm01 <0-1>            (-)
      real(kind=sp)  ::frootuptn   ! o  fraction root uptake nitrogen sm01                 (-)
      real(kind=sp)  ::frootuptp   ! o  fraction root uptake phosphorus                    (-)
      real(kind=sp)  ::cdinsm01w   ! o  average water concentration din for sm01           (m)
      real(kind=sp)  ::cpo4sm01w   ! o  average water concentration po4 for sm01           (m)
      real(kind=sp)  ::cco2sm01    ! o  average concentration co2+h2co3 for sm01           (m)
      real(kind=sp)  ::chco3sm01   ! o  average concentration hco3 for sm01                (m)
      real(kind=sp)  ::cdinsm01b   ! o  average sediment concentration din for sm01        (m)
      real(kind=sp)  ::cpo4sm01b   ! o  average sediment concentration po4 for sm01        (m)

      ! local

      integer(kind=int_32)  ::ikmrk1      !    first attribute
      integer(kind=int_32)  ::ikmrk2      !    second attribute
      real(kind=sp)  ::din         ! l  dissolved inorganic nitrogen, corrected for preference
      real(kind=sp)  ::z1          ! l  z1
      real(kind=sp)  ::fr_avg      ! l  fr_avg
      real(kind=sp)  ::hroot       ! l  hroot
      real(kind=sp)  ::limn        ! l  limn
      real(kind=sp)  ::limp        ! l  limp

      ! zero the average concentrations for all segments

      ipnt  = ipoint
      do iseg = 1 , noseg
         pmsa(ipnt(39)) = 0.0
         pmsa(ipnt(40)) = 0.0
         pmsa(ipnt(41)) = 0.0
         pmsa(ipnt(42)) = 0.0
         pmsa(ipnt(43)) = 0.0
         pmsa(ipnt(44)) = 0.0
         pmsa(ipnt(45)) = 0.0
         ipnt  = ipnt + increm
      enddo

      ! first loop average concentration over height and rooting depth

      ipnt  = ipoint
      do iseg = 1 , noseg

         depth      = pmsa( ipnt(  1) )
         totaldepth = pmsa( ipnt(  2) )
         locseddept = pmsa( ipnt(  4) )
         nh4s12     = pmsa( ipnt(  5) )
         po4s12     = pmsa( ipnt(  6) )
         ibotseg    = nint(pmsa( ipnt(  7) ))
         FrBmLay    = pmsa( ipnt(  8) )
         rootdesm01 = pmsa( ipnt(  9) )
         poros      = pmsa( ipnt( 10) )
         nh4        = max(pmsa( ipnt( 11) ),0.0)
         no3        = max(pmsa( ipnt( 12) ),0.0)
         po4        = max(pmsa( ipnt( 13) ),0.0)
         disco2     = max(pmsa( ipnt( 14) ),0.0)
         dish2co3   = max(pmsa( ipnt( 15) ),0.0)
         dishco3    = max(pmsa( ipnt( 16) ),0.0)
         prfnh4sm01 = pmsa( ipnt( 17) )
         temp       = pmsa( ipnt( 30) )

         ! convert co2 to carbon and add h2co3, calculated din with preference
         ! adjust for porosity

         din        = (nh4 + no3/prfnh4sm01)/poros
         po4        =  po4/poros
         disco2     = (disco2*12./44. + dish2co3)/poros
         dishco3    =  dishco3/poros

         call evaluate_waq_attribute(1,iknmrk(iseg),ikmrk1)
         if (ikmrk1.eq.1) then

            ! active water segment

            fr_avg = FrBmLay

            pmsa(botidx(39)) = pmsa(botidx(39)) + din*fr_avg
            pmsa(botidx(40)) = pmsa(botidx(40)) + po4*fr_avg
            pmsa(botidx(41)) = pmsa(botidx(41)) + disco2*fr_avg
            pmsa(botidx(42)) = pmsa(botidx(42)) + dishco3*fr_avg
            pmsa(botidx(45)) = pmsa(botidx(45)) + temp*fr_avg

            ! S12 sediment concentration

            call evaluate_waq_attribute(2,iknmrk(iseg),ikmrk2)
            if ((ikmrk2.eq.0).or.(ikmrk2.eq.3)) then
              if (nh4s12.gt.0.0) pmsa(botidx(43)) = nh4s12
              if (po4s12.gt.0.0) pmsa(botidx(44)) = po4s12
            endif

         elseif (ikmrk1.eq.3) then

            ! sediment bed segment, distribution of roots in bed

            hroot = min(rootdesm01,totaldepth)
            z1 = locseddept - depth

            if (hroot .gt. locseddept) then
               ! completely in segment:
               fr_avg = min(1.0,depth/hroot)
            elseif (hroot .gt. z1 ) then
               ! partialy in segment:
               fr_avg = (hroot-z1)/hroot
            else
               ! not in segment:
               fr_avg = 0.0
            endif

            pmsa(botidx(43)) = pmsa(botidx(43)) + din*fr_avg
            pmsa(botidx(44)) = pmsa(botidx(44)) + po4*fr_avg
            pmsa(botidx(45)) = pmsa(botidx(45)) + temp*fr_avg

         endif

         ipnt = ipnt + increm

      enddo

      ! second loop limiting factors and water / sediment uptake ratio

      ipnt  = ipoint
      do iseg = 1 , noseg

         call evaluate_waq_attribute(1,iknmrk(iseg),ikmrk1)
         if (ikmrk1.eq.1) then
            call evaluate_waq_attribute(2,iknmrk(iseg),ikmrk2)
            if ((ikmrk2.eq.0).or.(ikmrk2.eq.3)) then

               kmdinsm01w = pmsa( ipnt( 18) )
               kmpsm01w   = pmsa( ipnt( 19) )
               kmco2sm01  = pmsa( ipnt( 20) )
               kmhco3sm01 = pmsa( ipnt( 21) )
               kmdinsm01b = pmsa( ipnt( 22) )
               kmpsm01b   = pmsa( ipnt( 23) )
               cdinsm01w  = pmsa( ipnt( 24) )
               cpo4sm01w  = pmsa( ipnt( 25) )
               cco2sm01   = pmsa( ipnt( 26) )
               chco3sm01  = pmsa( ipnt( 27) )
               cdinsm01b  = pmsa( ipnt( 28) )
               cpo4sm01b  = pmsa( ipnt( 29) )

               ! n limitation

               if ( kmdinsm01b .lt. 1e-20 ) then
                  ! only in water
                  limnsm01w  = cdinsm01w/(cdinsm01w+kmdinsm01w)
                  limnsm01b  = -1.
               else
                  if ( kmdinsm01w .lt. 1e-20 ) then
                     ! only in bottom
                     limnsm01b  = cdinsm01b/(cdinsm01b+kmdinsm01b)
                     limnsm01w  = -1.
                  else
                     ! both
                     limnsm01w  = cdinsm01w/(cdinsm01w+kmdinsm01w)
                     limnsm01b  = cdinsm01b/(cdinsm01b+kmdinsm01b)
                  endif
               endif
               limn = max(limnsm01w,limnsm01b)
               if ( cdinsm01w .gt. 1e-10 ) then
                  if ( cdinsm01b .gt. 1e-10 ) then
                     frootuptn = .998/(1.+2.66*(cdinsm01b/cdinsm01w)**(-0.83))
                  else
                     frootuptn = 0.0
                  endif
               else
                  frootuptn = 1.0
               endif

               ! p limitation

               if ( kmpsm01b .lt. 1e-20 ) then
                  ! only in water
                  limpsm01w  = cpo4sm01w/(cpo4sm01w+kmpsm01w)
                  limpsm01b  = -1.
               else
                  if ( kmpsm01w .lt. 1e-20 ) then
                     ! only in bottom
                     limpsm01b  = cpo4sm01b/(cpo4sm01b+kmpsm01b)
                     limpsm01w  = -1.
                  else
                     ! both
                     limpsm01w  = cpo4sm01w/(cpo4sm01w+kmpsm01w)
                     limpsm01b  = cpo4sm01b/(cpo4sm01b+kmpsm01b)
                  endif
               endif
               limp = max(limpsm01w,limpsm01b)
               if ( cpo4sm01w .gt. 1e-10 ) then
                  if ( cpo4sm01b .gt. 1e-10 ) then
                     frootuptp = .998/(1.+2.66*(cpo4sm01b/cpo4sm01w)**(-0.83))
                  else
                     frootuptp = 0.0
                  endif
               else
                  frootuptp = 1.0
               endif

               ! c limitation

               if ( kmco2sm01 .lt. 1e-20 ) then
                  if ( kmhco3sm01 .lt. 1e-20 ) then
                     ! no c limitation calculation
                     lco2sm01 = 1.0
                  else
                     ! only hco3 limitation calculatation
                     lco2sm01 = chco3sm01/(chco3sm01+kmhco3sm01)
                  endif
               else
                  if ( kmhco3sm01 .lt. 1e-20 ) then
                     ! only co2 limitation calculatation
                     lco2sm01 = cco2sm01/(cco2sm01+kmco2sm01)
                  else
                     ! both co2 and hco3 limitation calculatation, take the max
                     lco2sm01 = max(cco2sm01/(cco2sm01+kmco2sm01),chco3sm01/(chco3sm01+kmhco3sm01))
                  endif
               endif

               ! overall limitation is minimum of n, p and c

               limnutsm01 = min(limn,limp,lco2sm01)

               pmsa( ipnt( 31)   ) = limnsm01w
               pmsa( ipnt( 32)   ) = limpsm01w
               pmsa( ipnt( 33)   ) = lco2sm01
               pmsa( ipnt( 34)   ) = limnsm01b
               pmsa( ipnt( 35)   ) = limpsm01b
               pmsa( ipnt( 36)   ) = limnutsm01
               pmsa( ipnt( 37)   ) = frootuptn
               pmsa( ipnt( 38)   ) = frootuptp

            endif
         endif
         ipnt        = ipnt        + increm

      enddo

      return
      contains

      ! Auxiliary function to compute the index for the associated bottom segment
      integer function botidx(number)
         integer, intent(in) :: number

         botidx = ipoint(number)+(ibotseg-1)*increm(number)
      end function botidx

      end subroutine

end module m_macnut
