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

module m_intrinsic_replacements
   use precision, only: dp

   implicit none
   private
   public :: cosd, sind

   real(kind=dp), parameter, public :: pi = 4.0_dp * atan(1.0_dp)
   real(kind=dp), parameter, public :: degrees_to_radians = (pi / 180.0)

contains

   pure elemental function cosd(degrees) result(res)
      real(kind=dp), intent(in) :: degrees
      real(kind=dp) :: res

      res = cos(degrees * degrees_to_radians)
   end function cosd

   pure elemental function sind(degrees) result(res)
      real(kind=dp), intent(in) :: degrees
      real(kind=dp) :: res

      res = sin(degrees * degrees_to_radians)
   end function sind

end module m_intrinsic_replacements
