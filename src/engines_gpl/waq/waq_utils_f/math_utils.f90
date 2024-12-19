!!  Copyright (C)  Stichting Deltares, 2012-2023.
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
module math_utils
    use m_waq_precision

    implicit none

    private
    public :: greatest_common_divisor, salinity_from_chloride, chlorinity_from_sal

    real(kind = real_wp), parameter  :: sal0 = 0.03  ! g/kg
    real(kind = real_wp), parameter  :: gtcl = 1.805 !

contains


    subroutine greatest_common_divisor(num_elements, numbers, gcd_result)

        ! calculates the greatest common divisor (gcd) (largest common denominator) for a set of numbers.
        integer(kind = int_wp) :: num_elements, gcd_result, i, inner_index
        integer(kind = int_wp) :: numbers(num_elements)
        integer(kind = int_wp) :: min_number
        logical :: divisor_found

        min_number = numbers(1)
        do i = 2, num_elements
            min_number = MIN(numbers(i), min_number)
        enddo

        do i = min_number, 1, -1
            divisor_found = .TRUE.
            do inner_index = 1, num_elements
                if (MOD(numbers(inner_index), i) /= 0) then
                    divisor_found = .FALSE.
                    exit
                end if
            end do
            if (divisor_found) then
                gcd_result = i
                exit
            end if
        end do

    end subroutine greatest_common_divisor

    ! Auxiliary function: convert chlorinity to salinity
    ! based on temperature
    !
    ! It might be useful to have the converse available,
    ! so here is a copy of the original code from SALCHL:
    !
    ! DENS = 1000. + 0.7 * SAL / (1 - SAL / 1000.) &
    !         - 0.0061 * (TEMP - 4.0) * (TEMP - 4.0)
    ! !
    ! IF (SAL <= SAL0) THEN
    !     SAL = 0.0
    ! ELSE
    !     SAL = SAL - SAL0
    ! ENDIF
    ! !
    ! !     g/m3 = (g/kg)*(kg/m3)/(g/g)
    ! !
    ! CL = SAL * DENS / GTCL
    !
    subroutine salinity_from_chloride( cl, temp, sal, density )
        real(kind = real_wp), intent(in)  :: cl
        real(kind = real_wp), intent(in)  :: temp
        real(kind = real_wp), intent(out) :: sal
        real(kind = real_wp), intent(out) :: density

        density = 1000.0 + 0.7 * cl / 1000.0 * gtcl &
                  - 0.0061 * (temp - 4.0) * (temp - 4.0)
        sal = cl * gtcl / density + sal0

    end subroutine salinity_from_chloride

    function chlorinity_from_sal( sal, temp ) result(cl)
        real(kind = real_wp), intent(in)  :: sal
        real(kind = real_wp), intent(in)  :: temp
        real(kind = real_wp)              :: cl

        real(kind = real_wp)              :: density

        density = 1000. + 0.7 * sal / (1 - sal / 1000.) &
                  - 0.0061 * (temp - 4.0) * (temp - 4.0)

        cl = max( 0.0_real_wp, sal - sal0 ) * density / gtcl

    end function chlorinity_from_sal

end module math_utils
