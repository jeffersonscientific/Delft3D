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
module m_output_to_polygon
   use precision, only: dp
   use messagehandling, only: idLen
   implicit none

   type, public :: t_variables_inside_polygon
      character(len=idLen) :: filename
      logical :: is_defined
      integer, allocatable, dimension(:) :: cell_indices
      integer, allocatable, dimension(:) :: link_indices
      integer, allocatable, dimension(:, :) :: link_to_nodes
      integer, allocatable, dimension(:) :: links_mask
      integer, allocatable, dimension(:,:) :: netlink_to_netnodes
      integer, allocatable, dimension(:) :: netnode_indices
      integer, allocatable, dimension(:) :: netlink_indices
      integer :: ndx
      integer :: ndxi
      integer :: ndx2d
      integer :: ndx1db
      integer :: lnx
      integer :: lnx1d
      integer :: lnx1db  
      integer :: lnxi
      integer :: numk
      integer :: numl
      integer :: numl1d
      real(kind=dp), pointer :: cell_values(:)
   contains
   procedure :: create_mask_arrays => create_mask_arrays_impl
   procedure :: findcells_inside_polygon => findcells_inside_polygon_impl
   procedure :: remap_cells => remap_cells_impl
   end type t_variables_inside_polygon

   interface 
      module subroutine create_mask_arrays_impl(this)
         implicit none
         class(t_variables_inside_polygon), intent(inout) :: this
      end subroutine create_mask_arrays_impl
   end interface 

   interface 
      module subroutine findcells_inside_polygon_impl(this, cells_mask)
         implicit none
         class(t_variables_inside_polygon), intent(inout) :: this
         integer, allocatable, dimension(:), intent(inout) :: cells_mask
      end subroutine findcells_inside_polygon_impl
   end interface 

   interface 
      module function remap_cells_impl(this, var, start_index, end_index) result(mapped_var)
         implicit none
         class(t_variables_inside_polygon), intent(inout) :: this
         real(kind=dp), dimension(:), intent(in) :: var
         real(kind=dp), dimension(:), pointer :: mapped_var
         integer, intent(in) :: start_index  
         integer, intent(in) :: end_index
      end function remap_cells_impl
      
   end interface 
   private

contains
end module m_output_to_polygon
