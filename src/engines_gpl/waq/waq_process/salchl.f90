module m_salchl
    use m_waq_precision
    use math_utils, only: salinity_from_chloride
    use m_logger_helper, only : stop_with_error, get_log_unit_number

    implicit none
    private
    public :: salchl

contains

    subroutine salchl (process_space_real, fl, ipoint, increm, num_cells, &
            noflux, iexpnt, iknmrk, num_exchanges_u_dir, num_exchanges_v_dir, &
            num_exchanges_z_dir, num_exchanges_bottom_dir)
        !>\file
        !>       Converts salinity into chloride or vice versa (Aquatic Chemistry 2nd ed 1981 p567)

        !----- GPL ---------------------------------------------------------------------
        !
        !  Copyright (C)  Stichting Deltares, 2011-2024.
        !
        !  This program is free software: you can redistribute it and/or modify
        !  it under the terms of the GNU General Public License as published by
        !  the Free Software Foundation version 3.
        !
        !  This program is distributed in the hope that it will be useful,
        !  but WITHOUT ANY WARRANTY; without even the implied warranty of
        !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        !  GNU General Public License for more details.
        !
        !  You should have received a copy of the GNU General Public License
        !  along with this program.  If not, see <http://www.gnu.org/licenses/>.
        !
        !  contact: delft3d.support@deltares.nl
        !  Stichting Deltares
        !  P.O. Box 177
        !  2600 MH Delft, The Netherlands
        !
        !  All indications and logos of, and references to, "Delft3D" and "Deltares"
        !  are registered trademarks of Stichting Deltares, and remain the property of
        !  Stichting Deltares. All rights reserved.
        !
        !-------------------------------------------------------------------------------
        !
        !
        !-------------------------------------------------------------------------------
        !
        !     Description of the module :
        !
        ! Name    T   L I/O   Description                                  Units
        ! ----    --- -  -    -------------------                           ----
        ! CL      R*4 1 I/O  chloride concentration                         [g/m3]
        ! SAL     R*4 1 I/O  salinity                                       [g/kg]
        ! SAL0    R*4 1 I    salinity at zero chloride                      [g/kg]
        ! GTCL    R*4 1 I    ratio of salinity and chloride                 [g/g]
        ! TEMP    R*4 1 I    ambient temperature                            [oC]
        ! DENS    R*4 1 -    densioty of water with dissolved salt          [kg/m3]
        ! SWSALCL R*4 1 I    option: 0 SAL simulated, 1 CL simulated
        !
        !   Logical Units : -
        !   Modules called : -

        !     Name     Type   Library
        !     ------   -----  ------------
        !
!!        IMPLICIT REAL    (A-H, J-Z)
!!        IMPLICIT INTEGER (I)
        !
        REAL(kind = real_wp) :: process_space_real  (*), FL    (*)
        INTEGER(kind = int_wp) :: IPOINT(*), INCREM(*), num_cells, NOFLUX, &
                IEXPNT(4, *), IKNMRK(*), num_exchanges_u_dir, num_exchanges_v_dir, num_exchanges_z_dir, num_exchanges_bottom_dir
        !
        REAL(kind = real_wp) :: CL, SAL, SAL0, GTCL, TEMP, DENS, SWSALCL
        integer(kind = int_wp) :: iseg, ip1, ip2, ip3, ip4, ip5, iflux
        integer(kind = int_wp) :: lunrep

        !
        IP1 = IPOINT(1)
        IP2 = IPOINT(2)
        IP3 = IPOINT(3)
        IP4 = IPOINT(4)
        IP5 = IPOINT(5)
        !
        IFLUX = 0
        DO ISEG = 1, num_cells

            IF (BTEST(IKNMRK(ISEG), 0)) THEN
                !
                CL = process_space_real(IP1)
                TEMP = process_space_real(IP2)
                SWSALCL = process_space_real(IP3)
                !
                !***********************************************************************
                !**** Processes connected to the normalization RIZA method
                !***********************************************************************
                !
                !     factor 0.7 in density correction was derived empirically from RIZA Standard methods
                !     table 210 on p 109 is repoduced within 0.15%
                !     basic relation sal-chlorinity: sal = 0.03 +1.805*chlor/density
                !     density = f(temp and salt concentration)
                !
                !     Note: gtcl and sal0 are fixed to 1.805 and 0.03 respectively
                !           Chlorinity expressed in g/m3, temperature in degrees C
                !
                !     Note: the switch is maintained to make sure that incorrect use is caught
                !           Of old the routine allowed the reverse conversion via this switch.
                !
                IF (NINT(SWSALCL) == 1) THEN
                    call salinity_from_chloride( cl, temp, sal, dens )
                ELSE
                    CALL get_log_unit_number(LUNREP)
                    WRITE(LUNREP, *) 'ERROR in SALCHL'
                    WRITE(LUNREP, *) 'Obsolete option for conversion - only the value 1 is allowed'
                    WRITE(LUNREP, *) 'Option in input:', SWSALCL
                    WRITE(*, *) 'ERROR in SALCHL'
                    WRITE(*, *) 'Obsolete option for conversion - only the value 1 is allowed'
                    WRITE(*, *) 'Option in input:', SWSALCL
                    CALL stop_with_error()
                ENDIF

                !
                process_space_real (IP4) = DENS
                process_space_real (IP5) = SAL
                !
            ENDIF
            !
            IFLUX = IFLUX + NOFLUX
            IP1 = IP1 + INCREM (1)
            IP2 = IP2 + INCREM (2)
            IP3 = IP3 + INCREM (3)
            IP4 = IP4 + INCREM (4)
            IP5 = IP5 + INCREM (5)
            !
        end do
        !
        RETURN
    END

end module m_salchl
