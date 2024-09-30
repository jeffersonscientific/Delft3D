!!  Copyright (C)  Stichting Deltares, 2024-2024.
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License version 3,
!!  as published by the Free Software Foundation.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program. If not, see <http://www.gnu.org/licenses/>.
!!
!!  contact: delft3d.support@deltares.nl
!!  Stichting Deltares
!!  P.O. Box 177
!!  2600 MH Delft, The Netherlands
!!
!!  All indications and logos of, and references to registered trademarks
!!  of Stichting Deltares remain the property of Stichting Deltares. All
!!  rights reserved.

module test_netcdf_define_mode
   use ftnunit
   use netcdf
   implicit none
   private

   public :: tests_netcdf_define_mode

contains

   subroutine tests_netcdf_define_mode()
      call test(test_define_mode_after_create, "Verify that NetCDF is in define mode after nf90_create is called")
   end subroutine tests_netcdf_define_mode

   subroutine test_define_mode_after_create()
      integer :: nc_id, status
      character(len=:), allocatable :: file_name
      integer :: dim_id, var_id
      integer, dimension(1) :: dim_ids
      integer, dimension(1) :: var_data

      file_name = "test_define_mode_after_create.nc"
      status = nf90_create(file_name, NF90_CLOBBER, nc_id)
      call assert_equal(status, NF90_NOERR, "nf90_create")

      status = nf90_enddef(nc_id)
      call assert_equal(status, NF90_NOERR, "nf90_enddef")

      status = nf90_close(nc_id)
      call assert_equal(status, NF90_NOERR, "nf90_close")
   end subroutine test_define_mode_after_create
end module test_netcdf_define_mode
