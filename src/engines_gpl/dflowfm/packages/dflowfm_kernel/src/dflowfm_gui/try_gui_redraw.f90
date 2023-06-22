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

!> Does a redraw of the Interacter GUI at the current time, if that is needed.
!! Non-GUI runs will directly return, do nothing extra.
!! A redraw is considered needed if current user timestep number is a new integer multiple
!! of NTEK refresh interval.
subroutine try_gui_redraw()
   use unstruc_display, only: jaGUI, ntek, plottofile
   use m_flowtimes, only: dnt_user

   integer          :: ndraw
   COMMON /DRAWTHIS/   ndraw(50)

   integer          :: key

   if (jaGUI == 0) then
      return
   end if
   
   key = 3                                          ! this part is for online visualisation
   if (ntek > 0) then
      if (mod(int(dnt_user),ntek) .eq. 0) then
            if (plottofile == 1) then
               ndraw(10) = plottofile 
            endif      
            call drawnu(key)
            if (key == 1) then
               continue !goto 1234
            endif
      endif
   endif
end subroutine try_gui_redraw
   