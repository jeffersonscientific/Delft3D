!!  Copyright (C)  Stichting Deltares, 2012-2024.
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
module m_dlwqj3
    use m_waq_precision

    implicit none

    private
    public :: dlwqj3

contains

    SUBROUTINE DLWQJ3 (binary_work_file, ascii_output_file_unit, output_file_width, num_blocks, integer_array, &
            real_array, values_arr, num_items, num_dims, iorder, &
            has_scale_factors, convert_breakpoint, info_from_binary_file, time_function_type, memory_type, &
            time_factor, is_date_format, is_yyddhh_format, cum_items, cum_dims, &
            flow_ignore_string, item_string, value_conc_string, breakpoint_harm_string, output_file_option)

        !! Prints and writes blocks of data
        !     output_file_option  INTEGER     1       INPUT   output file option

        use m_dlwqj2
        use date_time_utils, only : convert_time_format
        use timers       !   performance timers

        logical :: has_scale_factors, info_from_binary_file, defaults_on
        logical, intent(in) :: convert_breakpoint   !! T = Breakpoints are converted
        logical, intent(in) :: is_date_format       !! 'date'-format 1st time scale
        logical, intent(in) :: is_yyddhh_format     !! 'date'-format (F;ddmmhhss,T;yydddhh)
        character(len = *), intent(in)  :: item_string, value_conc_string, breakpoint_harm_string, flow_ignore_string(:)
        integer(kind = int_wp), intent(in) :: num_dims, num_items        !! number of subs, items to write
        integer(kind = int_wp), intent(in) :: iorder                 !! 1 = groups of subs per item
        integer(kind = int_wp), intent(in) :: ascii_output_file_unit !!Unit of ASCII output file(formatted output file)
        integer(kind = int_wp), intent(in) :: binary_work_file          !! Unit of binary work file
        integer(kind = int_wp), intent(inout) :: integer_array(:)
        integer(kind = int_wp), intent(in) :: time_function_type !!1 is block, 2 is linear, 3 is harmonics, 4 is fourier
        integer(kind = int_wp), intent(in) :: memory_type                           !! 0 is non permanent memory
        integer(kind = int_wp), intent(inout) :: cum_items, cum_dims    !! cumulative integer/real space count
        integer(kind = int_wp), intent(in) :: num_blocks                !! number of blocks to write
        integer(kind = int_wp) :: output_file_option
        integer(kind = int_wp) :: output_file_width         !! Width of the output file
        integer(kind = int_wp) :: time_factor               !! factor between clocks
        real(kind = real_wp), intent(inout) :: real_array(:), values_arr(:)

        integer(kind = int_wp) :: itel2, i1dum, i2dum, nodi2
        integer(kind = int_wp) :: k, ie, ie2, i1, i2
        integer(kind = int_wp) :: itels, itel, i3
        integer(kind = int_wp) :: ioffb, ioffi, ioffs, iskip, iskp2, notot, iss
        integer(kind = int_wp) :: ithndl = 0
        if (timon) call timstrt("dlwqj3", ithndl)

        ! write headers
        defaults_on = .false.
        if (num_dims < 0) defaults_on = .true.
        nodi2 = num_dims
        if (num_dims <= 0) nodi2 = 1
        if (iorder == 1) then
            write (ascii_output_file_unit, 1000) num_items, nodi2, value_conc_string
            write (binary_work_file) iorder, &
                    num_items, (integer_array(k), k = 1, num_items), &
                    num_dims, (integer_array(k), k = num_items + 1, num_items + num_dims), &
                    time_function_type, memory_type
        elseif (iorder == 2) then
            write (ascii_output_file_unit, 1000) nodi2, num_items, item_string
            if (binary_work_file > 0) &
                    write (binary_work_file) iorder, &
                            num_dims, (integer_array(k), k = 1, num_dims), &
                            num_items, (integer_array(k), k = num_dims + 1, num_dims + num_items), &
                            time_function_type, memory_type
        endif

        cum_items = cum_items + num_items + max(0, num_dims) + 5

        ! just declare array space for binary files and return
        if (info_from_binary_file) then
            write (ascii_output_file_unit, 1130) memory_type
            cum_items = cum_items + 3
            cum_dims = cum_dims + max(1, num_dims) * max(1, num_items) * 3
            goto 70
        endif

        if (num_blocks == 0) then
            has_scale_factors = .false.
            goto 70
        endif
        ioffb = num_items + nodi2 + 1
        ioffi = 0
        ioffs = num_items
        iskip = 1
        iskp2 = nodi2
        notot = num_items * nodi2
        if (time_function_type == 3 .or. time_function_type == 4) notot = notot + 1
        if (iorder == 2) then
            ioffi = max(num_dims, 0)
            ioffs = 0
            iskip = num_items
            iskp2 = 1
        endif

        ! scale factors
        iss = 1
        if (has_scale_factors) then
            has_scale_factors = .false.
            iss = 1
            if (output_file_option >= 4) then
                write (ascii_output_file_unit, 1010)
                do i2 = 1, num_dims, output_file_width
                    ie = min(i2 + output_file_width - 1, num_dims)
                    write (ascii_output_file_unit, 1020) (integer_array(ioffs + k), k = i2, ie)
                    write (ascii_output_file_unit, 1025) (flow_or_ignore(flow_ignore_string, integer_array(ioffs + k)), k = i2, ie)
                    write (ascii_output_file_unit, 1030) (real_array(k), k = i2, ie)
                end do
            endif
            do i1 = 1, num_blocks
                do i2 = 0, notot - 1
                    if (iorder == 1) itel2 = mod(i2, num_dims) + 1
                    if (iorder == 2) itel2 = i2 / num_dims + 1
                    values_arr(iss + i2) = values_arr(iss + i2) * real_array(itel2)
                end do
                iss = iss + notot
            end do
        endif

        ! convert breakpoints
        if (num_blocks > 1) then
            if (output_file_option >= 4) write (ascii_output_file_unit, 1040) breakpoint_harm_string, num_blocks
            if (.not. convert_breakpoint) &
                    call convert_time_format(integer_array(ioffb:), num_blocks, time_factor, is_date_format, is_yyddhh_format)
            if (defaults_on .and. output_file_option >= 4) write (ascii_output_file_unit, 1050)
        else
            if (defaults_on) then
                if (output_file_option >= 4) write (ascii_output_file_unit, 1050)
            else
                if (output_file_option >= 4) write (ascii_output_file_unit, 1060)
            endif
        endif

        ! write binary file
        if (binary_work_file > 0) then
            i1dum = 0
            i2dum = 0
            ! write table in binary format to wrk file.
            call write_breakpoint_data_blocks(binary_work_file, num_blocks, notot, 1, integer_array(ioffb:), values_arr, i1dum, i2dum)
            cum_items = cum_items + i1dum
            cum_dims = cum_dims + i2dum
        endif

        ! write formatted output
        if (output_file_option >= 4) then
            itels = 0
            do i1 = 1, num_blocks
                if (num_blocks > 1) then
                    if (time_function_type == 1) &
                            write (ascii_output_file_unit, 1070) breakpoint_harm_string, i1, integer_array(ioffb + i1 - 1)
                    if (time_function_type == 2) &
                            write (ascii_output_file_unit, 1070) breakpoint_harm_string, i1, integer_array(ioffb + i1 - 1)
                    if (time_function_type == 3) then
                        itels = itels + 1
                        write (ascii_output_file_unit, 1080) i1, integer_array(ioffb + i1 - 1), values_arr(itels)
                    endif
                    if (time_function_type == 4) then
                        itels = itels + 1
                        write (ascii_output_file_unit, 1090) i1, integer_array(ioffb + i1 - 1), values_arr(itels)
                    endif
                endif
                do i2 = 1, nodi2, output_file_width
                    ie2 = min(i2 + output_file_width - 1, nodi2)
                    if (num_dims > 0) then
                        write (ascii_output_file_unit, 1100) value_conc_string, (integer_array(ioffs + k), k = i2, ie2)
                        write (ascii_output_file_unit, 1150) item_string, &
                                (flow_or_ignore(flow_ignore_string, integer_array(ioffs + k)), k = i2, ie2)
                    endif
                    itel = itels
                    do i3 = 1, num_items
                        write (ascii_output_file_unit, 1120)  abs(integer_array(ioffi + i3)), &
                                (values_arr(itel + k), k = (i2 - 1) * iskip + 1, (ie2 - 1) * iskip + 1, iskip)
                        itel = itel + iskp2
                    end do
                end do
                itels = itels + nodi2 * num_items
            end do
        endif

        70 write (ascii_output_file_unit, 1140)
        if (timon) call timstop(ithndl)
        return

        1000 format (/' DATA grouped in', I10, ' blocks of', I10, ' ', A)
        1010 format (' Scale factors for this block of data: ')
        1020 format (' Scale    :', I6, 9I12)
        1025 format (' Substance:', 10('  ', A10))
        1030 format (' values   :', 10E12.4)
        1040 format (/' Number of ', A, 's with full data:', I5)
        1050 format (' Default values in this block.')
        1060 format (' Constant values in this block.')
        1070 format (' ', A, ' ', I7, ' :', I10)
        1080 format (' Harmonic: ', I3, ' :', I10, ' Phase: ', 10E12.4)
        1090 format (' Fourier : ', I3, ' :', I10, ' Phase: ', 10E12.4)
        1100 format (' ', A, I20, 9I12)  ! ( ' ',A,I6,9I12)
        1150 format (' ', A, ' ', 10('  ', A10))
        1120 format (I10, 2X, 1P, 10E12.4)
        1130 format (' Info comes at runtime from binary file at unit: ', I3)
        1140 format(/' ====> input item completed <==== '//)

    end subroutine dlwqj3

    character*20 function flow_or_ignore(string, i)
        integer(kind = int_wp) :: i
        character(len = *) :: string(*)
        if (i > 0) then
            flow_or_ignore = string(i)
        elseif (i == 0) then
            flow_or_ignore = 'FLOW'
        else
            flow_or_ignore = 'ignored'
        endif
    end function flow_or_ignore

end module m_dlwqj3
