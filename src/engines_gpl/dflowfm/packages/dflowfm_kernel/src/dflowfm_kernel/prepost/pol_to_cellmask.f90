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

!> update cellmask from samples
!> a cell is dry when it is:
!>   1) inside ANY "1"-polygon (drypnts), OR
!>   2) outside ALL "-1"-polygons (enclosures)
module m_pol_to_cellmask

   use precision, only: dp
   implicit none

   private

   public :: pol_to_cellmask

   ! Module-level variables for cellmask polygon checking
   real(kind=dp), allocatable :: xpmin_cellmask(:), ypmin_cellmask(:)
   real(kind=dp), allocatable :: xpmax_cellmask(:), ypmax_cellmask(:)
   real(kind=dp), allocatable :: zpl_start_cellmask(:)
   real(kind=dp), allocatable :: xpl_cellmask(:), ypl_cellmask(:)
   integer, allocatable :: iistart_cellmask(:), iiend_cellmask(:)
   integer :: Npoly_cellmask = 0
   logical :: cellmask_initialized = .false.

contains

   subroutine pol_to_cellmask()
   use network_data, only: cellmask, nump1d2d, npl, nump, xzw, yzw, xpl, ypl, zpl

   if (allocated(cellmask)) deallocate(cellmask)
   allocate(cellmask(nump1d2d))
   cellmask = 0

   if (NPL == 0) return

   ! Initialize once
   call dbpinpol_cellmask_init(NPL, xpl, ypl, zpl)

   ! Vectorized assignment - beautiful and simple!
   cellmask(1:nump) = dbpinpol_cellmask(xzw(1:nump), yzw(1:nump))

   ! Cleanup
   call dbpinpol_cellmask_cleanup()

end subroutine pol_to_cellmask

  !> Initialize polygon data structures for cellmask checking.
   !! Must be called once before using dbpinpol_cellmask.
   subroutine dbpinpol_cellmask_init(NPL, xpl, ypl, zpl)
      use m_alloc
      use m_missing, only: dp, dmiss

      implicit none

      integer, intent(in) :: NPL
      real(kind=dp), intent(in) :: xpl(NPL), ypl(NPL), zpl(NPL)

      integer :: MAXPOLY
      integer :: ipoint, istart, iend, ipoly

      if (cellmask_initialized) then
         call dbpinpol_cellmask_cleanup()
      end if

      if (NPL == 0) then
         cellmask_initialized = .true.
         return
      end if

      MAXPOLY = 1000
      call realloc(xpmin_cellmask, maxpoly, keepExisting=.false.)
      call realloc(xpmax_cellmask, maxpoly, keepExisting=.false.)
      call realloc(ypmin_cellmask, maxpoly, keepExisting=.false.)
      call realloc(ypmax_cellmask, maxpoly, keepExisting=.false.)
      call realloc(iistart_cellmask, maxpoly, keepExisting=.false.)
      call realloc(iiend_cellmask, maxpoly, keepExisting=.false.)
      call realloc(zpl_start_cellmask, maxpoly, keepExisting=.false.)

      ipoint = 1
      ipoly = 0

      do while (ipoint < NPL)
         ipoly = ipoly + 1
         if (ipoly > maxpoly) then
            maxpoly = ceiling(maxpoly * 1.1)
            call realloc(xpmin_cellmask, maxpoly, keepExisting=.true.)
            call realloc(xpmax_cellmask, maxpoly, keepExisting=.true.)
            call realloc(ypmin_cellmask, maxpoly, keepExisting=.true.)
            call realloc(ypmax_cellmask, maxpoly, keepExisting=.true.)
            call realloc(iistart_cellmask, maxpoly, keepExisting=.true.)
            call realloc(iiend_cellmask, maxpoly, keepExisting=.true.)
            call realloc(zpl_start_cellmask, maxpoly, keepExisting=.true.)
         end if

         ! Get polygon start and end pointer
         call get_startend(NPL - ipoint + 1, xpl(ipoint:NPL), ypl(ipoint:NPL), istart, iend, dmiss)
         istart = istart + ipoint - 1
         iend = iend + ipoint - 1

         if (istart >= iend .or. iend > NPL) exit

         ! Store bounding box
         xpmin_cellmask(ipoly) = minval(xpl(istart:iend))
         xpmax_cellmask(ipoly) = maxval(xpl(istart:iend))
         ypmin_cellmask(ipoly) = minval(ypl(istart:iend))
         ypmax_cellmask(ipoly) = maxval(ypl(istart:iend))

         iistart_cellmask(ipoly) = istart
         iiend_cellmask(ipoly) = iend
         
         ! Store polygon type (sign of zpl)
         if (zpl(istart) /= dmiss) then
            zpl_start_cellmask(ipoly) = zpl(istart)
         else
            zpl_start_cellmask(ipoly) = dmiss
         end if

         ! Advance pointer
         ipoint = iend + 2
      end do
      
      Npoly_cellmask = ipoly
      cellmask_initialized = .true.

   end subroutine dbpinpol_cellmask_init

  !> Check if point should be masked.
   elemental function dbpinpol_cellmask(xp, yp) result(mask)
      use m_missing, only: dmiss  !< Use globals directly
      implicit none

      integer :: mask
      real(kind=dp), intent(in) :: xp, yp

      integer :: istart, iend, ipoly, in_test
      logical :: found_inside_enclosure
      integer :: num_enclosures

      mask = 0
      if (.not. cellmask_initialized) return

      found_inside_enclosure = .false.
      num_enclosures = 0

      do ipoly = 1, Npoly_cellmask
         istart = iistart_cellmask(ipoly)
         iend = iiend_cellmask(ipoly)

         if (zpl_start_cellmask(ipoly) == dmiss) cycle

         ! Bounding box check
         if (xp < xpmin_cellmask(ipoly) .or. xp > xpmax_cellmask(ipoly) .or. &
             yp < ypmin_cellmask(ipoly) .or. yp > ypmax_cellmask(ipoly)) cycle

         ! Point-in-polygon test using globals dmiss and JINS
         !call pinpok(xp, yp, iend - istart + 1, &
         !            xpl_cellmask(istart:iend), ypl_cellmask(istart:iend), &
         !            in_test, JINS, dmiss)
         in_test = pinpok_elemental(xp, yp, iend - istart + 1)

         if (zpl_start_cellmask(ipoly) > 0.0_dp) then
            ! Dry point polygon
            if (in_test == 1) then
               mask = 1
               return
            end if
         else if (zpl_start_cellmask(ipoly) < 0.0_dp) then
            ! Enclosure polygon
            num_enclosures = num_enclosures + 1
            if (in_test == 1) found_inside_enclosure = .true.
         end if
      end do

      ! Outside all enclosures = dry
      if (num_enclosures > 0 .and. .not. found_inside_enclosure) mask = 1

   end function dbpinpol_cellmask

   !> Clean up module-level cellmask polygon data structures.
   subroutine dbpinpol_cellmask_cleanup()
      implicit none

      if (allocated(xpmin_cellmask)) deallocate(xpmin_cellmask)
      if (allocated(xpmax_cellmask)) deallocate(xpmax_cellmask)
      if (allocated(ypmin_cellmask)) deallocate(ypmin_cellmask)
      if (allocated(ypmax_cellmask)) deallocate(ypmax_cellmask)
      if (allocated(zpl_start_cellmask)) deallocate(zpl_start_cellmask)
      if (allocated(iistart_cellmask)) deallocate(iistart_cellmask)
      if (allocated(iiend_cellmask)) deallocate(iiend_cellmask)
      
      Npoly_cellmask = 0
      cellmask_initialized = .false.

      return
   end subroutine dbpinpol_cellmask_cleanup

   !> Optimized elemental point-in-polygon test using ray casting algorithm.
   !! Accesses polygon data via module arrays.
   elemental function pinpok_elemental(xl, yl, ipoly) result(inside)
      use m_missing, only: dmiss, JINS  !< Use globals directly

      implicit none

      real(kind=dp), intent(in) :: xl, yl     !< Point coordinates (scalar)
      integer, intent(in) :: ipoly            !< Polygon index
      integer :: inside                       !< Result: 1=inside, 0=outside     

      integer :: i, j, istart, iend, crossings
      real(kind=dp) :: x1, x2, y1, y2, xinters

      inside = 0
      
      ! Get polygon bounds from module variables
      istart = iistart_cellmask(ipoly)
      iend = iiend_cellmask(ipoly)
      
      if (iend - istart + 1 <= 2) then
         inside = 1
         goto 999
      end if

      ! Ray-casting algorithm
      crossings = 0
      j = iend
      
      do i = istart, iend
         if (xpl_cellmask(i) == dmiss) exit
         
         x1 = xpl_cellmask(j)
         y1 = ypl_cellmask(j)
         x2 = xpl_cellmask(i)
         y2 = ypl_cellmask(i)
         
         ! Check if point is on vertex
         if (xl == x1 .and. yl == y1) then
            inside = 1
            goto 999
         end if
         
         ! Check if ray crosses edge
         if ((y1 > yl) .neqv. (y2 > yl)) then
            xinters = x1 + (yl - y1) * (x2 - x1) / (y2 - y1)
            
            if (xl < xinters) then
               crossings = crossings + 1
            else if (xl == xinters) then
               inside = 1
               goto 999
            end if
         end if
         
         j = i
      end do

      inside = mod(crossings, 2)
      999   continue
      if (jins == 0) inside = 1 - inside

   end function pinpok_elemental

end module m_pol_to_cellmask
