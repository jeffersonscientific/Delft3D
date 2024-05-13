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

!> Finds the flow nodes/cell numbers for each given x,y point (e.g., an observation station)
subroutine find_flownode(n, xobs, yobs, namobs, kobs, jakdtree, jaoutside, iLocTp)
   use unstruc_messages
   use m_partitioninfo
   use m_flowgeom
   use m_GlobalParameters, only: INDTP_1D, INDTP_2D, INDTP_ALL
   use kdtree2Factory
   use geometry_module, only: dbdistance
   use m_missing, only: dmiss
   use m_sferic, only: jsferic, jasfer3D

   implicit none

   integer,                            intent(in   )  :: n           !< number of points
   double precision,     dimension(n), intent(in   )  :: xobs, yobs  !< points coordinates
   character(len=IdLen), dimension(n), intent(in   )  :: namobs      !< names of points
   integer,              dimension(n), intent(inout)  :: kobs        !< associated flow nodes, if found.
   integer,                            intent(inout)  :: jakdtree    !< use kdtree (1) or not (other)
   integer,                            intent(in   )  :: jaoutside   !< allow outside cells (for 1D) (1) or not (0)
   integer,                            intent(in   )  :: iLocTp      !< Node type, one of INDTP_1D/2D/ALL.
   integer                                            :: ierror      !< error (1) or not (0)
   integer                                            :: i, k, k1b
   integer,           dimension(1)                    :: idum
   double precision                                   :: d1, d2

   ierror = 1

   if (jakdtree == 1) then
      call find_flowcells_kdtree(treeglob, n, xobs, yobs, kobs, jaoutside, iLocTp, ierror)

      if (jampi == 1) then
         ! globally reduce ierror
         idum(1) = ierror
         call reduce_int_max(1, idum)
         ierror = idum(1)
      end if

      if (ierror /= 0) then
         jakdtree = 0   ! retry without kdtree
      end if

      ! disable observation stations without attached flowlinks
      do i = 1, n
         k = kobs(i)
         if (k > 0) then
            if (nd(k)%lnx < 1) then
               kobs(i) = 0
            end if
         end if
      end do
   
   else
      
      do i = 1, n
         call inflowcell(xobs(i), yobs(i), k, jaoutside, iLocTp)
         if (jaoutside == 1 .and. (iLocTp == INDTP_1D .or. iLocTp == INDTP_ALL)) then
            call closeto1Dorbnd(xobs(i), yobs(i), k1b)
            if (k /= 0 .and. k1b /= 0) then
                d1 = dbdistance(xz(k1b), yz(k1b), xobs(i), yobs(i), jsferic, jasfer3D, dmiss)
                d2 = dbdistance(xz(k  ), yz(k  ), xobs(i), yobs(i), jsferic, jasfer3D, dmiss)
                if (d1 < d2) then
                   k = k1b
                end if
            else if (k1b /= 0) then
                k = k1b
            end if
         end if
         kobs(i) = 0
         if (k /= 0) then
            if ( nd(k)%lnx > 0 ) then
               kobs(i) = k
            end if
         end if
      end do
   end if

   if (jampi == 1 .and. n > 0) then
      ! check if this subdomain owns the observation station
      call reduce_kobs(n, kobs, xobs, yobs, jaoutside)
   end if

   do i = 1, n
      if (kobs(i) == 0) then
          write(msgbuf, '(a,i0,a,a,a)') 'Could not find flowcell for point #', i, ' (', trim(namobs(i)), '). Discarding.'
          call msg_flush()
      end if
   end do

   ierror = 0
1234 continue

   end subroutine find_flownode
