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
module m_write_balance_output
    use m_waq_precision

    implicit none

    private
    public :: write_balance_output

contains


    SUBROUTINE write_balance_output(LUBAL, FILBAL, ITIME, MONAME, NOTOT, &
            NOFLUX, SYNAME, NDMPAR, DANAME, ASMASS, &
            FLXINT, NOTOT2, CONC2, INIT)
        ! Writes balance output

        !
        !     NAME    KIND     LENGTH     FUNCT.  DESCRIPTION
        !     ----    -----    ------     ------- -----------
        !     LUBAL   INTEGER       1     INPUT   Logical unit balance file
        !     FILBAL  CHAR*(*)      1     INPUT   Name balance file
        !     ITIME   INTEGER       1     INPUT   Simulation time ( scu )
        !     MONAME  CHAR*40       4     INPUT   model identification
        !     NOTOT   INTEGER       1     INPUT   Total number of substances
        !     NOFLUX  INTEGER       1     INPUT   Nr. of fluxes
        !     SYNAME  CHAR*20   NOTOT     INPUT   names of substances
        !     NDMPAR  INTEGER       1     INPUT   Number of dump segments
        !     DANAME  CHAR*20   NDMPAR    INPUT   names of monitoring stations
        !     ASMASS  REAL NOTOT*NDMPAR*6 INPUT   Mass balance terms
        !     FLXINT  REAL  NOFLUX*NDMPAR INPUT   Integrated fluxes
        !     NOTOT2  REAL          1     INPUT   Number of extra variables
        !     CONC2   REAL NOTOT2*NDMPAR  INPUT   Extra variables
        !     INIT    INTEGER       1     IN/OUT  Init flag (1=yes,!1=no)

        use timers

        integer(kind = int_wp) :: lubal, itime, init, notot, noflux, &
                ndmpar, notot2, nopout
        real(kind = real_wp) :: asmass(notot, ndmpar, 6), flxint(noflux, ndmpar), &
                conc2(notot2, ndmpar)
        character(len = 20)  syname(*), daname(*)
        character(len = 40)  moname(4)
        character(len = *) filbal
        integer(kind = int_wp) :: j, i, k, isys, iflx, ihlp
        integer(kind = int_wp) :: ithandl = 0
        if (timon) call timstrt ("write_balance_output", ithandl)

        ! Initialize file
        if (init == 1) then
            init = 0

            ! write header
            write (lubal) (moname(i), i = 1, 4)
            nopout = 6 * notot + noflux + 2
            write (lubal) nopout, ndmpar, notot
            write (lubal) (syname(i), i = 1, notot)
            write (lubal) (daname(i), i = 1, ndmpar)
        endif

        ! Perform output
        write (lubal) itime, (&
                ((asmass(isys, j, k), k = 1, 6), isys = 1, notot), &
                (flxint(iflx, j), iflx = 1, noflux), &
                (conc2(ihlp, j), ihlp = 1, 2), &
                j = 1, ndmpar)

        if (timon) call timstop (ithandl)

    end subroutine write_balance_output

end module m_write_balance_output
