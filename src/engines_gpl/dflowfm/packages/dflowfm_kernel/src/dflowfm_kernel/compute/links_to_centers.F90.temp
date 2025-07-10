!----- AGPL --------------------------------------------------------------------
!
!  Copyright (C)  Stichting Deltares, 2017-2025.
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
module m_links_to_centers
   use precision, only: dp, sp
   use m_flow, only: lnkx, ndkx, kmx, kmxn, ktop, kbot, lbot, ltop
   use m_flowgeom, only: lnx, ln, wcL, ndx

   implicit none
   private

   public :: links_to_centers

   interface links_to_centers
      module procedure links_to_centers_dp
      module procedure links_to_centers_sp
   end interface links_to_centers
contains
   !> Set flow node value based on flow link values, where vlin is real(kind=dp)
   subroutine links_to_centers_dp(vnod, vlin)
      real(kind=dp), intent(in) :: vlin(lnkx)
      include 'links_to_centers_body.inc'
   end subroutine links_to_centers_dp

   subroutine links_to_centers_sp(vnod, vlin)
      real(kind=sp), intent(in) :: vlin(lnkx)
      include 'links_to_centers_body.inc'
   end subroutine links_to_centers_sp

end module m_links_to_centers
