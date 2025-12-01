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
submodule(m_output_to_polygon) m_output_to_polygon_sub
   use precision, only: dp
   use messagehandling, only: idLen
   implicit none

contains

   ! Create a mask array indicating which elements are inside the polygon
   module subroutine create_mask_arrays_impl(this)
      use messagehandling, only: msgbuf, fatal_flush
      use m_reapol, only: reapol
      use geometry_module, only: dbpinpol
      use m_alloc, only: realloc
      use m_flowgeom, only: ndx, lnx, ln, ndxi, ndx2d, ndx1db, lnx1d, lnx1db, lnxi

      class(t_variables_inside_polygon), intent(inout) :: this

      integer :: k
      integer, allocatable, dimension(:) :: cells_mask, links_mask

      if (allocated(this%cell_indices)) then
         return
      end if
      allocate (this%cell_indices(ndx))
      allocate (cells_mask(ndx))
      allocate (this%link_indices(lnx))
      allocate (this%link_to_nodes(2, lnx))
      allocate (links_mask(lnx))

      if (.not. this%is_defined) then
         this%ndx = ndx
         this%lnx = lnx
         this%ndxi = ndxi
         this%ndx2d = ndx2d
         this%ndx1db = ndx1db
         this%cell_indices = [(k, k=1, ndx)]
         this%link_indices = [(k, k=1, lnx)]
         this%link_to_nodes = ln
         return
      end if

      associate (cell_indices => this%cell_indices, &
                 link_indices => this%link_indices, &
                 num_vertices_inside => this%lnx, &
                 link_to_nodes => this%link_to_nodes)

         call this%findcells_inside_polygon(cells_mask)

         this%ndxi = get_newindex(cells_mask, ndxi)
         this%ndx2d = get_newindex(cells_mask, ndx2d)
         this%ndx1db = get_newindex(cells_mask, ndx1db)

         call find_vertices_inside_polygon(this, cells_mask)
         this%lnx1d = get_newindex(links_mask, lnx1d)
         this%lnx1db = get_newindex(links_mask, lnx1db)
         this%lnxi = get_newindex(links_mask, lnxi)

         call find_netnodes_inside_polygon(this, cells_mask)
      end associate

      call realloc(this%cell_indices, this%ndx, &
                   keepExisting=.true.)
      call realloc(this%link_indices, this%lnx, &
                   keepExisting=.true.)
      call realloc(this%link_to_nodes, [2, this%lnx], &
                   keepExisting=.true.)
      allocate (this%cell_values(this%ndx))
   end subroutine create_mask_arrays_impl

   function get_newindex(cells_mask, old_ndx) result(new_ndx)
      integer, intent(in) :: cells_mask(:)
      integer, intent(in) :: old_ndx
      integer :: new_ndx
      integer :: k

      new_ndx = 0
      do k = old_ndx, 1, -1
         if (cells_mask(k) > 0) then
            new_ndx = cells_mask(k)
            exit
         end if
      end do
   end function get_newindex

   module function remap_cells_impl(this, var, start_index, end_index) result(mapped_var)
      implicit none
      class(t_variables_inside_polygon), intent(inout) :: this
      real(kind=dp), dimension(:), intent(in) :: var
      real(kind=dp), dimension(:), pointer :: mapped_var
      integer, intent(in) :: start_index
      integer, intent(in) :: end_index

      integer :: i

      do i = start_index, end_index
         this%cell_values(i - start_index + 1) = var(this%cell_indices(i))
      end do
      mapped_var => this%cell_values(1:end_index - start_index + 1)
   end function remap_cells_impl

   module subroutine findcells_inside_polygon_impl(this, cells_mask)
      use m_reapol, only: reapol
      use m_polygon, only: npl, xpl, ypl, zpl
      use geometry_module, only: dbpinpol
      use m_missing, only: dmiss, jins
      use m_cell_geometry, only: xz, yz
      use m_flowgeom, only: ndx
      use m_polygon, only: npl, xpl, ypl, zpl, savepol, restorepol
      use m_filez, only: oldfil, doclose

      class(t_variables_inside_polygon), intent(inout) :: this
      integer, allocatable, dimension(:), intent(inout) :: cells_mask
      integer :: k, inside, minp
      call savepol()
      call oldfil(minp, this%filename)
      call reapol(minp, 0)

      associate (cell_indices => this%cell_indices, num_cells_inside => this%ndx)

         inside = -1
         num_cells_inside = 0
         do k = 1, ndx
            call dbpinpol(xz(k), yz(k), cells_mask(k), dmiss, jins, npl, xpl, ypl, zpl)
            if (cells_mask(k) == 1) then
               num_cells_inside = num_cells_inside + 1
               cells_mask(k) = num_cells_inside
               cell_indices(num_cells_inside) = k
            end if
         end do

      end associate

      call restorepol()
      call doclose(minp)

   end subroutine findcells_inside_polygon_impl

   subroutine find_vertices_inside_polygon(this, cells_mask)
      use m_flowgeom, only: lnx, ln

      class(t_variables_inside_polygon), intent(inout) :: this
      integer, dimension(:), intent(in) :: cells_mask

      integer :: L

      associate (link_indices => this%link_indices, &
                 link_to_nodes => this%link_to_nodes, &
                 num_vertices_inside => this%lnx, &
                 links_mask => this%links_mask)
         num_vertices_inside = 0
         links_mask = 0
         do L = 1, lnx
            if (cells_mask(ln(1, L)) > 0 .and. cells_mask(ln(2, L)) > 0) then
               num_vertices_inside = num_vertices_inside + 1
               link_indices(num_vertices_inside) = L
               links_mask(L) = num_vertices_inside
               link_to_nodes(1, num_vertices_inside) = cells_mask(ln(1, L))
               link_to_nodes(2, num_vertices_inside) = cells_mask(ln(2, L))
            end if
         end do
      end associate

   end subroutine find_vertices_inside_polygon

   subroutine find_netnodes_inside_polygon(this, cells_mask)
      use m_flowgeom, only: nd
      use network_data, only: numk, numL
      use m_alloc, only: realloc

      class(t_variables_inside_polygon), intent(inout) :: this
      integer, dimension(:), intent(in) :: cells_mask

      integer, dimension(:), allocatable :: netnodes_mask

      integer :: i, k, L, node1, node2, max_netnodes

      allocate (netnodes_mask(numk), this%netnode_indices(numk))
      allocate (this%netnode_indices(numk))
      allocate (this%netlink_indices(numL))
      netnodes_mask = 0

      associate (ndx => this%ndx, &
                 netnode_indices => this%netnode_indices)
         do k = 1, ndx
            do i = 1, size(nd(k)%nod)
               netnodes_mask(nd(k)%nod(i)) = 1
            end do
         end do

         max_netnodes = 0
         do k = 1, numk
            if (netnodes_mask(k) == 1) then
               max_netnodes = max_netnodes + 1
               netnode_indices(max_netnodes) = k
               netnodes_mask(k) = max_netnodes
            end if
         end do
      end associate
      call realloc(this%netnode_indices, max_netnodes, &
                   keepExisting=.true.)


   end subroutine find_netnodes_inside_polygon

end submodule m_output_to_polygon_sub
