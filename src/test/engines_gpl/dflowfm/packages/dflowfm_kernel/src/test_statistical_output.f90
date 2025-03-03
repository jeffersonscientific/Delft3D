!!  Copyright (C)  Stichting Deltares, 2025-2025.
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

module test_statistical_output
   use m_statistical_output
   use ftnunit
   use precision, only: dp
   implicit none
   private

   real(kind=dp), parameter :: test_tolerance = 1e-3_dp

   public :: tests_statistical_output

contains
   subroutine tests_statistical_output
      call test(test_realloc_unallocated, 'Tests realloc function with unallocated array')
      call test(test_realloc_allocated, 'Tests realloc function with allocated array')
      call test(test_realloc_crop_no_elements, 'Tests realloc function with crop and no elements')
      call test(test_realloc_crop_one_element, 'Tests realloc function with crop and an element')
      call test(test_dealloc, 'Tests dealloc function')
   end subroutine tests_statistical_output

   subroutine test_realloc_unallocated()
      type(t_output_variable_set) :: output_set
      call assert_equal(allocated(output_set%statout), .false., '')
      call assert_equal(output_set%capacity, 0, '')

      call realloc(output_set)
      call assert_equal(allocated(output_set%statout), .true., '')
      call assert_equal(output_set%capacity > 0, .true., '')
   end subroutine test_realloc_unallocated

   subroutine test_realloc_allocated()
      type(t_output_variable_set) :: output_set
      integer, parameter :: value_to_test = 496

      output_set%count = 3
      call realloc(output_set)
      output_set%statout(3)%id_var = value_to_test

      output_set%count = 300
      call realloc(output_set)
      call assert_equal(allocated(output_set%statout), .true., '')
      call assert_equal(output_set%capacity >= 300, .true., '')
      call assert_equal(output_set%statout(3)%id_var, value_to_test, '')
   end subroutine test_realloc_allocated

   subroutine test_realloc_crop_no_elements()
      type(t_output_variable_set) :: output_set
      call realloc(output_set, .true.)
      call assert_equal(allocated(output_set%statout), .true., '')
      call assert_equal(output_set%capacity, output_set%count, '')
   end subroutine test_realloc_crop_no_elements

   subroutine test_realloc_crop_one_element()
      type(t_output_variable_set) :: output_set
      output_set%count = 1
      call realloc(output_set, .true.)
      call assert_equal(allocated(output_set%statout), .true., '')
      call assert_equal(output_set%capacity, output_set%count, '')
   end subroutine test_realloc_crop_one_element

   subroutine test_dealloc()
      type(t_output_variable_set) :: output_set
      call realloc(output_set)
      call assert_equal(allocated(output_set%statout), .true., '')

      call dealloc(output_set)
      call assert_equal(allocated(output_set%statout), .false., '')
   end subroutine test_dealloc
end module test_statistical_output
