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

! $Id$
! $HeadURL$

!> Checks whether or not a user-interrupt was recently entered in the Interacter GUI.
!! Non-GUI runs will always return .false..
!! User-interrupt is typically the S-key pressed, or the left mouse button.
function check_gui_interrupt() result(was_interrupted)
   use unstruc_display, only: jaGUI
   use m_partitioninfo, only: jampi, reduce_key
   use MessageHandling
   logical :: was_interrupted !< Whether or not a user-interrupt was recently entered in the GUI.

   was_interrupted = .false.

   if (jaGUI == 1) then

      call get_s_key(key)

      if (jampi == 1) then
         call reduce_key(key)
      end if

      if (key == 1) then
         call mess(LEVEL_INFO, 'User interrupt')
         was_interrupted = .true.
      end if
   end if
end function check_gui_interrupt