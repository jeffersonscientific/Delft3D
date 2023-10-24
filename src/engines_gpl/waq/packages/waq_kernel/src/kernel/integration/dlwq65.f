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
      module m_dlwq65
      use m_waq_type_definitions


      implicit none

      contains


      SUBROUTINE DLWQ65 ( AMASS  , VOLUME , IDT    , NOSEG  )
!
!     Deltares     SECTOR WATERRESOURCES AND ENVIRONMENT
!
!     CREATED             : june 1988 by L.Postma
!
!     FUNCTION            : makes a closure error correction
!                           for steady state computations
!
!     LOGICAL UNITNUMBERS : none
!
!     SUBROUTINES CALLED  : none
!
!     PARAMETERS          :
!
!     NAME    KIND     LENGTH     FUNCT.  DESCRIPTION
!     ----    -----    ------     ------- -----------
!     AMASS   REAL(kind=sp) ::NOSEG     IN/OUT  closure error correction
!     VOLUME  REAL(kind=sp) ::NOSEG     INPUT   volume
!     IDT     INTEGER(kind=int_32) ::1       INPUT   Time between AMASS and VOLUME
!     NOSEG   INTEGER(kind=int_32) ::1       INPUT   number of segments
!
      use timers

      real(kind=sp) ::AMASS (  *) , VOLUME(*)
      integer(kind=int_32) ::idt, noseg
      integer(kind=int_32) ::i
      integer(kind=int_32) ::ithandl = 0
      if ( timon ) call timstrt ( "dlwq65", ithandl )
!
!         loop oversehe number of segments
!
      DO 10 I = 1 , NOSEG
      AMASS(I) = ( AMASS(I)-VOLUME(I) )/IDT
   10 CONTINUE
!
      if ( timon ) call timstop ( ithandl )
      RETURN
      END

      end module m_dlwq65
