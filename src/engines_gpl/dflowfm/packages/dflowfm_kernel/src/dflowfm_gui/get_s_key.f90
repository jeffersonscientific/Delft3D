!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2023.                                
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

      !> Detects S-keypress or left mouse, indicating user-interrupt in the GUI.
      SUBROUTINE get_s_key(was_pressed)
      implicit none
      integer, intent(  out) :: was_pressed !< Whether or not (1/0) S-key or left-mouse button was pressed since previous time an Interacter event was checked.

      integer :: numkey

      was_pressed = 0

      CALL INKEYEVENTIMM(NUMKEY)

      ! Check lower- and uppercase 's' and left mouse click
      IF (NUMKEY == 115 .or. NUMKEY == 115-32 .or. NUMKEY == 251) then
         was_pressed = 1
         call inflush()
      endif
      RETURN
      END SUBROUTINE get_s_key
