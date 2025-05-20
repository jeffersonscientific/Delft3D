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
module m_averad
   use m_waq_precision

   implicit none

contains

   subroutine AVERAD(process_space_real, FL, IPOINT, INCREM, num_cells, &
                     NOFLUX, IEXPNT, IKNMRK, num_exchanges_u_dir, num_exchanges_v_dir, &
                     num_exchanges_z_dir, num_exchanges_bottom_dir)
      use m_extract_waq_attribute
      use m_logger_helper, only: get_log_unit_number, write_error_message, stop_with_error

      implicit none

      !     arguments

      real(kind=real_wp) :: process_space_real(*) ! in/out input-output array space to be adressed with IPOINT/INCREM
      real(kind=real_wp) :: FL(*) ! in/out flux array
      integer(kind=int_wp) :: IPOINT(*) ! in     start index input-output parameters in the process_space_real array (segment or exchange number 1)
      integer(kind=int_wp) :: INCREM(*) ! in     increment for each segment-exchange for the input-output parameters in the process_space_real array
      integer(kind=int_wp) :: num_cells ! in     number of segments
      integer(kind=int_wp) :: NOFLUX ! in     total number of fluxes (increment in FL array)
      integer(kind=int_wp) :: IEXPNT(4, *) ! in     exchange pointer table
      integer(kind=int_wp) :: IKNMRK(*) ! in     segment features array
      integer(kind=int_wp) :: num_exchanges_u_dir ! in     number of exchanges in first direction
      integer(kind=int_wp) :: num_exchanges_v_dir ! in     number of exchanges in second direction
      integer(kind=int_wp) :: num_exchanges_z_dir ! in     number of exchanges in third direction
      integer(kind=int_wp) :: num_exchanges_bottom_dir ! in     number of exchanges in fourth direction

      !     from process_space_real array

      !     8 inputs, 3 outputs.

      integer(kind=int_wp) :: IP1 !<REAL RADSURF        input 1, actual irradiation at the water surface            (W/m2)
      integer(kind=int_wp) :: IP2 !<REAL AveRadTIni     input 2, starting time for average irradiance      TINIT       (d)
      integer(kind=int_wp) :: IP3 !<REAL AveRadPeri     input 3, averaging period of irradiance    PERIOD           (d)
      integer(kind=int_wp) :: IP4 !<REAL ITIME          input 4, DELWAQ time                                        (s)
      integer(kind=int_wp) :: IP5 !<REAL DELT           input 5, timestep for processes          (d)
      integer(kind=int_wp) :: IP6 !<REAL AuxSys         input 6, ratio between days and system clock   (scu/d, i.e. system clock unit/day)
      integer(kind=int_wp) :: IP7 !<REAL SumAveRad      input 7,   Work array for summing over time  (W/m2)
      integer(kind=int_wp) :: IP8 !<REAL SumAveRadT     input 8,   Count of times   TCOUNT  (d)
      integer(kind=int_wp) :: IP9 !<REAL  SumAveRad     output 9,   Work array for summing over time  (W/m2)
      integer(kind=int_wp) :: IP10 !<REAL SumAveRadT    output 10,  Count of times   TCOUNT  (d)
      integer(kind=int_wp) :: IP11 !<REAL RadSurfAve    output 11,  average irradiance over the period              (W/m2)

      integer(kind=int_wp) :: IN1, IN2, IN3, IN4, IN5, &
                              IN6, IN7, IN8, IN9, IN10, IN11
      integer(kind=int_wp) :: IKMRK, ISEG
      integer(kind=int_wp) :: IACTION, lunrep
      integer(kind=int_wp) :: ATTRIB
      real(kind=real_wp) :: TINIT, PERIOD, TIME, DELT, TCOUNT, AuxSys

      integer(kind=int_wp), parameter :: MAXWARN = 50
      integer(kind=int_wp), save :: NOWARN = 0

      call get_log_unit_number(lunrep)

      !     IACTION is in 3 parts. 0, 2, 3.

      IP1 = IPOINT(1)
      IP2 = IPOINT(2)
      IP3 = IPOINT(3)
      IP4 = IPOINT(4)
      IP5 = IPOINT(5)
      IP6 = IPOINT(6)
      IP7 = IPOINT(7)
      IP8 = IPOINT(8)
      IP9 = IPOINT(9)
      IP10 = IPOINT(10)
      IP11 = IPOINT(11)

      IN1 = INCREM(1)
      IN2 = INCREM(2)
      IN3 = INCREM(3)
      IN4 = INCREM(4)
      IN5 = INCREM(5)
      IN6 = INCREM(6)
      IN7 = INCREM(7)
      IN8 = INCREM(8)
      IN9 = INCREM(9)
      IN10 = INCREM(10)
      IN11 = INCREM(11)

      TINIT = process_space_real(IP2)
      PERIOD = process_space_real(IP3)
      TIME = process_space_real(IP4)
      DELT = process_space_real(IP5)
      AuxSys = process_space_real(IP6)

      if (PERIOD < DELT) then
         call write_error_message('AveRadSurf: Period of averaging should be larger than DELWAQ time step.')
      end if

      !
      ! TIME in second should be changed to TIME in days.
      !
      TIME = TIME / AuxSys
      !
      !      Start and stop criteria are somewhat involved:
      !      - The first time for the first period is special, as this
      !        is the only time there is no previous period.
      !      - If there is a previous period, update the averages
      !        for that period and reset the accumulative values
      !        for the next
      !
      !      To formulate the ideas more clearly:
      !      - The first period is a closed interval
      !      - All other periods are half-open intervals (the last time
      !        of the previous period should not be reused.)
      !
      IACTION = 0
      if (TIME >= TINIT - 0.5 * DELT) then
         IACTION = 2
         if (TIME <= TINIT + 0.5 * DELT) then
            do ISEG = 1, num_cells
               IP9 = IPOINT(9) + (ISEG - 1) * INCREM(9)
               IP10 = IPOINT(10) + (ISEG - 1) * INCREM(10)
               process_space_real(IP9) = 0.0
               process_space_real(IP10) = 0.0
            end do
         end if
      end if

      if (TIME >= TINIT + PERIOD - 0.5 * DELT .and. TIME <= TINIT + PERIOD + 0.5 * DELT) then
         IACTION = 3
      end if

      if (IACTION == 0) return

      do ISEG = 1, num_cells
         !
         !           Keep track of the time within the current quantile specification
         !           that each segment is active
         !
         TCOUNT = process_space_real(IP8) + DELT
         process_space_real(IP10) = TCOUNT

         process_space_real(IP9) = process_space_real(IP7) + process_space_real(IP1) * DELT

         !
         !        Always do the final processing whether the segment is active at this moment or not
         !

         if (IACTION == 3) then
            if (TCOUNT > 0.0) then
               process_space_real(IP11) = process_space_real(IP9) / TCOUNT
            else
               process_space_real(IP11) = 0.0

               if (NOWARN < MAXWARN) then
                  call extract_waq_attribute(3, IKNMRK(ISEG), ATTRIB)
                  if (ATTRIB /= 0) then
                     NOWARN = NOWARN + 1
                     write (lunrep, '(a,i0)') 'Periodic average of RadSurf could not be determined for segment ', ISEG
                     write (lunrep, '(a)') '    - division by zero. Average set to zero'

                     if (NOWARN == MAXWARN) then
                        write (lunrep, '(a)') '(Further messages suppressed)'
                     end if
                  end if
               end if
            end if

            !
            !           Reset for the next round
            !

            process_space_real(IP9) = 0.0
            process_space_real(IP10) = 0.0

         end if

         IP1 = IP1 + IN1
         IP7 = IP7 + IN7
         IP8 = IP8 + IN8
         IP9 = IP9 + IN9
         IP10 = IP10 + IN10
         IP11 = IP11 + IN11
      end do
      !
      !     Be sure to also reset the initial time, so that we can restart the
      !     averaging for the next period
      !
      if (IACTION == 3) then
         process_space_real(IP2) = TINIT + PERIOD
      end if

      return
   end
end module m_averad
