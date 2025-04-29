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

!
!

module m_addbaroclinicpressure

   implicit none

   private

   public :: addbaroclinicpressure

contains

   !> Computes and adds the baroclinic pressure gradient contributions to the momentum equations
   subroutine addbaroclinicpressure()
      use precision, only: dp, comparereal
      use m_addbarocl, only: addbarocL, addbarocLrho_w, addbarocL_use_rho_directly
      use m_addbarocn, only: addbarocn, addbarocnrho_w, addbarocn_use_rho_directly
      use m_addbaroc, only: addbaroc
      use m_flowgeom, only: lnxi, lnx, ndx
      use m_flow, only: hu, kmx
      use m_turbulence, only: rvdn, grn
      use m_physcoef, only: jabarocponbnd, rhointerfaces

      implicit none

      integer :: LL, l_bot, l_top, cell_index_2d, lnxbc

      if (jabarocponbnd == 0) then
         lnxbc = lnxi
      else
         lnxbc = lnx
      end if

      if (kmx == 0) then
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(LL)
         do LL = 1, lnxbc
            if (comparereal(hu(LL), 0.0_dp) == 0) then
               cycle
            end if
            call addbaroc(LL)
         end do
         !$OMP END PARALLEL DO
      else

         rvdn = 0.0_dp
         grn = 0.0_dp

         if (rhointerfaces == 0) then

            !$OMP PARALLEL DO &
            !$OMP PRIVATE(cell_index_2d)
            do cell_index_2d = 1, ndx
               call addbarocn(cell_index_2d)
            end do
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO &
            !$OMP PRIVATE(LL,l_bot,l_top)
            do LL = 1, lnxbc
               if (comparereal(hu(LL), 0.0_dp) == 0) then
                  cycle
               end if
               call getLbotLtop(LL, l_bot, l_top)
               if (l_top < l_bot) then
                  cycle
               end if
               call addbarocL(LL, l_bot, l_top)
            end do
            !$OMP END PARALLEL DO

         elseif (rhointerfaces == 1) then

            !$OMP PARALLEL DO &
            !$OMP PRIVATE(cell_index_2d)
            do cell_index_2d = 1, ndx
               call addbarocnrho_w(cell_index_2d)
            end do
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO &
            !$OMP PRIVATE(LL,l_bot,l_top)
            do LL = 1, lnxbc
               if (comparereal(hu(LL), 0.0_dp) == 0) then
                  cycle
               end if
               call getLbotLtop(LL, l_bot, l_top)
               if (l_top < l_bot) then
                  cycle
               end if
               call addbarocLrho_w(LL, l_bot, l_top)
            end do
            !$OMP END PARALLEL DO

         elseif (rhointerfaces == 2) then

            !$OMP PARALLEL DO &
            !$OMP PRIVATE(cell_index_2d)
            do cell_index_2d = 1, ndx
               call addbarocn_use_rho_directly(cell_index_2d)
            end do
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO &
            !$OMP PRIVATE(LL,l_bot,l_top)
            do LL = 1, lnxbc
               if (comparereal(hu(LL), 0.0_dp) == 0) then
                  cycle
               end if
               call getLbotLtop(LL, l_bot, l_top)
               if (l_top < l_bot) then
                  cycle
               end if
               call addbarocL_use_rho_directly(LL, l_bot, l_top)
            end do
            !$OMP END PARALLEL DO

         end if
      end if
   end subroutine addbaroclinicpressure
end module m_addbaroclinicpressure
