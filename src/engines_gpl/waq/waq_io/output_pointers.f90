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
module output_pointers
    use m_waq_precision

    implicit none

    private
    public :: get_output_pointers, set_default_output

contains


    subroutine get_output_pointers(noutp, nrvar, nrvarm, dlwnam, iopoin, &
            nmis, notot, syname, nocons, coname, &
            nopa, paname, nofun, funame, nosfun, &
            sfname, lurep)

        !! Sets the pointers for all extra vars

        use timers       !   performance timers

        integer(kind = int_wp), intent(in) :: noutp                  !< Number of output files
        integer(kind = int_wp), intent(in) :: nrvar (noutp)         !< No of output vars per file
        integer(kind = int_wp), intent(in) :: nrvarm                 !< Maximum of output variables p.p.
        integer(kind = int_wp), intent(out) :: nmis                   !< Number of missing input vars
        character*(*), intent(in) :: dlwnam(nrvarm, noutp)  !< Name of input variables
        integer(kind = int_wp), intent(out) :: iopoin(nrvarm, noutp)   !< Number of missing input vars
        integer(kind = int_wp), intent(in) :: notot                  !< Total number of substances
        integer(kind = int_wp), intent(in) :: nopa                   !< Number of parameters
        integer(kind = int_wp), intent(in) :: nosfun                 !< Number of segment functions
        character(20), intent(in) :: syname(notot)         !< Names of systems
        integer(kind = int_wp), intent(in) :: nocons                 !< Number of constants used
        integer(kind = int_wp), intent(in) :: nofun                  !< Number of functions ( user )
        character(20), intent(in) :: coname(nocons)        !< Constant names
        character(20), intent(in) :: paname(nopa)        !< Parameter names
        character(20), intent(in) :: funame(nofun)        !< Function names
        character(20), intent(in) :: sfname(nosfun)        !< Segment function names
        integer(kind = int_wp), intent(in) :: lurep                  !< Unit nr. report file

        character(20) varnam            ! Name of variable to be identified
        integer(kind = int_wp) :: ivarip             ! Pointer in the SSA
        integer(kind = int_wp) :: iout               ! loop variable of outputs
        integer(kind = int_wp) :: inrv               ! loop variable output number
        integer(kind = int_wp) :: ithndl = 0
        if (timon) call timstrt("get_output_pointers", ithndl)

        write(lurep, *)
        write(lurep, *) ' Determining the place of the output variables'
        write(lurep, *)

        nmis = 0

        do iout = 1, noutp
            do inrv = 1, nrvar(iout)
                varnam = dlwnam(inrv, iout)
                if (varnam == ' ') then
                    ivarip = 0
                else
                    call set_output_pointers (notot, nopa, nosfun, syname, nocons, &
                            nofun, coname, paname, funame, sfname, &
                            varnam, ivarip, lurep)
                    if (ivarip == -1) then
                        nmis = nmis + 1
                        write(lurep, '(3a)') '   INFO:', varnam, &
                                '; NOT FOUND, delwaq will detect variables from process library'
                    endif
                endif
                iopoin(inrv, iout) = ivarip
            end do
        end do

        if (timon) call timstop(ithndl)
        return
    end subroutine get_output_pointers

    subroutine set_output_pointers(notot, nopa, nosfun, syname, nocons, &
            nofun, coname, paname, funame, sfname, &
            varnam, ivarip, lurep)

        !! sets pointers for output variables
        !!     Logical unitnumbers : lurep   - report file

        use timers       !   performance timers
        use m_string_utils

        integer(kind = int_wp), intent(in) :: notot              !< Total number of substances
        integer(kind = int_wp), intent(in) :: nopa               !< Number of parameters
        integer(kind = int_wp), intent(in) :: nosfun             !< Number of segment functions
        character(20), intent(in) :: syname(notot)     !< Names of systems
        integer(kind = int_wp), intent(in) :: nocons             !< Number of constants used
        integer(kind = int_wp), intent(in) :: nofun              !< Number of functions ( user )
        character(20), intent(in) :: coname(nocons)    !< Constant names
        character(20), intent(in) :: paname(nopa)    !< Parameter names
        character(20), intent(in) :: funame(nofun)    !< Function names
        character(20), intent(in) :: sfname(nosfun)    !< Segment function names
        character(20), intent(in) :: varnam            !< Name of variable to be identified
        integer(kind = int_wp), intent(out) :: ivarip             !< Pointer in the SSA
        integer(kind = int_wp), intent(in) :: lurep              !< Unit nr. report file

        !     Local

        integer(kind = int_wp), parameter :: nopred = 6
        integer(kind = int_wp) :: indx            !  index in array of names
        character(20) predef(3)
        integer(kind = int_wp) :: ithndl = 0
        if (timon) call timstrt("set_output_pointers", ithndl)

        predef(1) = 'volume'
        predef(2) = 'itime'
        predef(3) = 'idt'

        !     determine how VAL is modelled

        indx = index_in_array(varnam, predef)
        if (indx == 1) then
            write(lurep, *) '       ', varnam, '; Using DELWAQ VOLUME'
            ivarip = 1
            goto 800
        endif
        if (indx == 2) then
            write(lurep, *) '       ', varnam, '; Using DELWAQ ITIME'
            ivarip = 2
            goto 800
        endif
        if (indx == 3) then
            write(lurep, *) '       ', varnam, '; Using DELWAQ IDT'
            ivarip = 3
            goto 800
        endif

        !     as model variable ?

        indx = index_in_array(varnam, syname)
        if (indx > 0) then
            write(lurep, *) '       ', varnam, '; Using substance nr', indx
            ivarip = nopred + nocons + nopa + nofun + nosfun + indx
            goto 800
        endif

        !     as segment function ?

        indx = index_in_array(varnam, sfname)
        if (indx > 0) then
            write(lurep, *) '       ', varnam, '; Using segment function nr', indx
            ivarip = nopred + nocons + nopa + nofun + indx
            goto 800
        endif

        !     as function ?

        indx = index_in_array(varnam, funame)
        if (indx > 0) then
            write(lurep, *) '       ', varnam, '; Using function nr', indx
            ivarip = nopred + nocons + nopa + indx
            goto 800
        endif

        !     as parameter ?

        indx = index_in_array(varnam, paname)
        if (indx > 0) then
            write(lurep, *) '       ', varnam, '; Using parameter nr', indx
            ivarip = nopred + nocons + indx
            goto 800
        endif

        !     as constant ?

        indx = index_in_array(varnam, coname)
        if (indx > 0) then
            write(lurep, *) '       ', varnam, '; Using constant nr', indx
            ivarip = nopred + indx
            goto 800
        endif

        !     not found

        ivarip = -1

        800 continue

        if (timon) call timstop(ithndl)
        return
    end subroutine set_output_pointers

    subroutine set_default_output(noutp, nrvar, iostrt, iostop, iostep, &
            isrtou, igrdou)
        !! Sets default output behavior

        use timers       !   performance timers
        use results, only : ncopt, idmp, ihi3, imap, iba3, iba2, ibal, ima2, imo4, imo3
        use m_sysi          ! Timer characteristics

        integer(kind = int_wp), intent(in) :: noutp          !< Number of output files
        integer(kind = int_wp), intent(out) :: nrvar (noutp)  !< Number of extra output vars
        integer(kind = int_wp), intent(out) :: iostrt(noutp)  !< Output start time (scu)
        integer(kind = int_wp), intent(out) :: iostop(noutp)  !< Output stop time (scu)
        integer(kind = int_wp), intent(out) :: iostep(noutp)  !< Output step time (scu)
        integer(kind = int_wp), intent(out) :: isrtou(noutp)  !< Sort output indication
        integer(kind = int_wp), intent(out) :: igrdou(noutp)  !< Output grid indication

        integer(kind = int_wp), parameter :: igseg = 1, igmon = 2, iggrd = 3, igsub = 4
        integer(kind = int_wp) :: ithndl = 0
        if (timon) call timstrt("set_default_output", ithndl)

        ! set output system default action
        ! MONitor file
        iostrt(1) = imstrt
        iostop(1) = imstop
        iostep(1) = imstep
        nrvar (1) = 0
        isrtou(1) = imo3
        igrdou(1) = igsub

        ! GRID file
        iostrt(2) = idstrt
        iostop(2) = idstop
        iostep(2) = idstep
        nrvar (2) = 0
        isrtou(2) = idmp
        igrdou(2) = iggrd

        ! HIStory file
        iostrt(3) = ihstrt
        iostop(3) = ihstop
        iostep(3) = ihstep
        nrvar (3) = 0
        isrtou(3) = ihi3
        igrdou(3) = igsub

        ! MAP file
        iostrt(4) = idstrt
        iostop(4) = idstop
        iostep(4) = idstep
        nrvar (4) = 0
        isrtou(4) = imap
        igrdou(4) = igseg

        ! BAL file
        if (noutp >= 5) then
            iostrt(5) = imstrt
            iostop(5) = imstop
            iostep(5) = imstep
            nrvar (5) = 0
            if (mod(intopt, 64) >= 32) then
                isrtou(5) = iba3
            elseif (mod(intopt, 32) >= 16) then
                isrtou(5) = iba2
            else
                isrtou(5) = ibal
            endif
            igrdou(5) = igsub
        endif

        ! NEFIS HIS file
        if (noutp >= 6) then
            iostrt(6) = ihstrt
            iostop(6) = ihstop
            iostep(6) = ihstep
            nrvar (6) = 0
            isrtou(6) = 0
            igrdou(6) = igsub
        endif

        ! NEFIS MAP file
        if (noutp >= 7) then
            iostrt(7) = idstrt
            iostop(7) = idstop
            iostep(7) = idstep
            nrvar (7) = 0
            isrtou(7) = 0
            igrdou(7) = igseg
        endif

        ! STATistical output MAP
        if (noutp >= 8) then
            iostrt(8) = itstop
            iostop(8) = itstop
            iostep(8) = idt
            nrvar (8) = 0
            isrtou(8) = ima2
            igrdou(8) = igseg
        endif

        ! STATistical output MON
        if (noutp >= 9) then
            iostrt(9) = itstop
            iostop(9) = itstop
            iostep(9) = idt
            nrvar (9) = 0
            isrtou(9) = imo4
            igrdou(9) = igsub
        endif

        if (timon) call timstop(ithndl)
        return
    end subroutine set_default_output

end module output_pointers
