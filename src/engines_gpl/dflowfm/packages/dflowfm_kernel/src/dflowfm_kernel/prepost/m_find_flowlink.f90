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
   
module m_find_flowlink
   
   implicit none
   
   private
   
   public :: find_1D_or_boundary_flowlink_bruteforce
   
   contains
   
   !> Find the flowlinks with the shortest perpendicular distances to the points [xx,yy]
   !! Uses the k-d tree routines
   subroutine find_flowlinks_kdtree(xx, yy, link_ids_closest)
      use stdlib_kinds, only: dp
      use MessageHandling, only: mess, LEVEL_ERROR
      use m_flowgeom, only: lnx, ln, xz, yz
      use m_alloc, only: realloc
      use kdtree2Factory, only: kdtree_instance, find_nearest_sample_kdtree
      use geometry_module, only: dlinedis
      use m_sferic, only: jsferic, jasfer3D
      use m_missing, only: dmiss
      
      real(dp), dimension(:), intent(in   ) :: xx, yy           !< x,y-coordinates
      integer,  dimension(:), intent(  out) :: link_ids_closest !< ids of the flowlinks whose midpoints lie closest to the points [xx,yy]
      
      real(dp), dimension(lnx)             :: flowlink_midpoints_x
      real(dp), dimension(lnx)             :: flowlink_midpoints_y
      integer                              :: link_id, ka, kb
      integer                              :: number_of_points, i_point
      real(dp), dimension(lnx)             :: zs_dummy
      type(kdtree_instance)                :: treeinst
      integer, parameter                   :: n_nearest_kdtree = 100
      integer, dimension(n_nearest_kdtree) :: link_ids_closest_midpoint
      integer                              :: ierror
      integer                              :: i_sample
      real(dp)                             :: dist_perp, dist_perp_min
      
      if (size(xx) /= size(yy) .or. size(xx) /= size(link_ids_closest)) then
         call mess(LEVEL_ERROR,'find_flowlinks: unmatched input array size')
      end if
      
      link_ids_closest = 0
      
      ! Calculate the x,y-coordinates of the midpoints of all the flowlinks
      do link_id = 1, lnx
         ka = ln(1,link_id)
         kb = ln(2,link_id)
         flowlink_midpoints_x(link_id) = (xz(ka) + xz(kb)) / 2.0_dp
         flowlink_midpoints_y(link_id) = (yz(ka) + yz(kb)) / 2.0_dp
      end do
      
      number_of_points = size(xx)
      do i_point = 1, number_of_points
         
         ! The k-d tree uses the distance to the midpoints of the flowlinks
         ! which is not actually what we want (namely the perpendicular distance).
         ! To solve this, we query the nearest n points (default 10) and then
         ! do a brute-force search on this (very) small subset
         call find_nearest_sample_kdtree(treeinst, lnx, 2, flowlink_midpoints_x, flowlink_midpoints_y, zs_dummy, &
                                         xx(i_point), yy(i_point), n_nearest_kdtree, link_ids_closest_midpoint, &
                                         ierror, jsferic, dmiss)
         
         ! Now loop over the n nearest points to find the one with the shortest perpendicular distance
         dist_perp_min   = huge(dist_perp_min)
         do i_sample = 1, n_nearest_kdtree
            link_id = link_ids_closest_midpoint(i_sample)
            call perpendicular_distance_to_flowlink(xx(i_point), yy(i_point), link_id, dist_perp)
            if (dist_perp < dist_perp_min) then
               link_ids_closest(i_point) = link_id
               dist_perp_min = dist_perp
            end if
         end do
         
      end do
      
   end subroutine find_flowlinks_kdtree
   
   !> Find the flowlinks with the shortest perpendicular distances to the points [xx,yy]
   !! Brute-force approach: simply check all flowlinks in the entire grid
   subroutine find_flowlinks_bruteforce(xx, yy, link_ids_closest)
      use stdlib_kinds, only: dp
      use MessageHandling, only: mess, LEVEL_ERROR
      use m_flowgeom, only: lnx, ln, xz, yz
      use m_sferic, only: jsferic, jasfer3D
      use geometry_module, only: dlinedis
      use m_missing, only: dmiss
      
      real(dp), dimension(:), intent(in   ) :: xx, yy           !< x,y-coordinates
      integer,  dimension(:), intent(  out) :: link_ids_closest !< ids of the flowlinks whose midpoints lie closest to the points [xx,yy]
      
      integer  :: number_of_points, i_point
      integer  :: link_id
      real(dp) :: dist_perp, dist_perp_min
      
      if (size(xx) /= size(yy) .or. size(xx) /= size(link_ids_closest)) then
         call mess(LEVEL_ERROR,'find_flowlinks: unmatched input array size')
      end if
      
      link_ids_closest = 0
      
      number_of_points = size(xx)
      do i_point = 1, number_of_points
         dist_perp_min   = huge(dist_perp_min)
         do link_id = 1, lnx
            call perpendicular_distance_to_flowlink(xx(i_point), yy(i_point), link_id, dist_perp)
            if (dist_perp < dist_perp_min) then
               link_ids_closest(i_point) = link_id
               dist_perp_min = dist_perp
            end if
         end do
      end do
      
   end subroutine find_flowlinks_bruteforce
   
   !> Find the 1-D or boundary flowlink with the shortest perpendicular distance to the point [x,y]
   !! Brute-force approach: simply check all flowlinks in the entire grid
   subroutine find_1D_or_boundary_flowlink_bruteforce(x, y, link_id_closest)
      use stdlib_kinds, only: dp
      use m_flowgeom, only: lnx, lnx1D, lnxi, ln, xz, yz, ndx
      use m_sferic, only: jsferic, jasfer3D
      use geometry_module, only: dlinedis
      use m_missing, only: dmiss
      
      real(dp), intent(in   ) :: x, y            !< x,y-coordinates
      integer,  intent(  out) :: link_id_closest !< id of the flowlink whose midpoint lies closest to the point [x,y]
      
      integer  :: link_id
      real(dp) :: dist_perp, dist_perp_min
      integer  :: ja
      
      link_id_closest = 0
      dist_perp_min   = huge(dist_perp_min)
   
      do link_id = 1, lnx
         if (link_id <= lnx1D .or. link_id > lnxi) then ! Only check 1-D and boundary flowlinks
            call perpendicular_distance_to_flowlink(x, y, link_id, dist_perp, ja)
            if (ja == 1) then
               if (dist_perp < dist_perp_min) then
                  link_id_closest = link_id
                  dist_perp_min = dist_perp
               end if
            end if
         end if
      end do
      
   end subroutine find_1D_or_boundary_flowlink_bruteforce
   
   !> Calculate the perpendicular distance from the point [x,y] to a flowlink
   subroutine perpendicular_distance_to_flowlink(x, y, link_id, perpendicular_distance, ja)
      use stdlib_kinds, only: dp
      use m_flowgeom, only: lnx, ln, xz, yz
      use m_sferic, only: jsferic, jasfer3D
      use geometry_module, only: dlinedis
      use m_missing, only: dmiss
      
      real(dp),          intent(in   ) :: x, y                   !< Coordinates of the point
      integer,           intent(in   ) :: link_id                !< id of the flowlink
      real(dp),          intent(  out) :: perpendicular_distance !< Perpendicular distance from the point [x,y] to a flowlink
      integer, optional, intent(  out) :: ja                     !< Whether or not (1/0) the computation was possible. If line points 1 and 2
                                                                 !! coincide, ja==0, and distance is just Euclidean distance between 3 and 1. 
      integer  :: ka, kb
      real(dp) :: xa, ya, xb, yb
      integer  :: ja_
      real(dp) :: xn, yn
      
      ka = ln(1,link_id)
      kb = ln(2,link_id)
      xa = xz(ka)
      ya = yz(ka)
      xb = xz(kb)
      yb = yz(kb)
      
      call dlinedis(x, y, xa, ya, xb, yb, ja_, perpendicular_distance, xn, yn, jsferic, jasfer3D, dmiss)
      
      if (present(ja)) then
         ja = ja_
      end if
            
   end subroutine perpendicular_distance_to_flowlink
   
end module m_find_flowlink