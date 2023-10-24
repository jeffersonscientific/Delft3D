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
      module m_dlwq70
      use m_waq_type_definitions


      implicit none

      contains


      SUBROUTINE DLWQ70 ( DISP   , DISPER , AREA   , FLOW   , ALENG  ,
     *                    VELO   , BOUND  , IPOINT , NOTOT  , ISYS   ,
     *                    NSYS   , NOQ1   , NOQ2   , NOQ    , NODISP ,
     *                    NOVELO , IDPNT  , IVPNT  , DERIV  , AMAT   ,
     *                                      JTRACK , IOPT   , ILFLAG )
!
!     Deltares     SECTOR WATERRESOURCES AND ENVIRONMENT
!
!     CREATED             : june 1988 by L.Postma
!
!     FUNCTION            : Fills matrix according to central
!                                       differencing in space.
!
!     LOGICAL UNITNUMBERS : none
!
!     SUBROUTINES CALLED  : none
!
!     PARAMETERS          :
!
!     NAME    KIND     LENGTH     FUNCT.  DESCRIPTION
!     ----    -----    ------     ------- -----------
!     DISP    REAL(kind=sp) ::3       INPUT   dispersion in 3 directions
!     DISPER  REAL(kind=sp) ::NODISP*NOQ   INPUT   additional dispersion array
!     AREA    REAL(kind=sp) ::NOQ      INPUT   exchange surface area
!     FLOW    REAL(kind=sp) ::NOQ      INPUT   flows accross exchange surfs
!     ALENG   REAL(kind=sp) ::2*NOQ     INPUT   from- and to lengthes
!     VELO    REAL(kind=sp) ::NOVELO*NOQ   INPUT   additional velocity array
!     BOUND   REAL(kind=sp) ::NOTOT*?    INPUT   boundary concentrations
!     IPOINT  INTEGER(kind=int_32) ::4*NOQ     INPUT   exchange pointers
!     NOTOT   INTEGER(kind=int_32) ::1       INPUT   number of total substances
!     ISYS    INTEGER(kind=int_32) ::1       INPUT   system number considered
!     NSYS    INTEGER(kind=int_32) ::1       INPUT   number of systems considered
!     NOQ1    INTEGER(kind=int_32) ::1       INPUT   nr of exchanges in first dir.
!     NOQ2    INTEGER(kind=int_32) ::1       INPUT   nr of exchanges in second dir.
!     NOQ     INTEGER(kind=int_32) ::1       INPUT   total number of exchanges
!     NODISP  INTEGER(kind=int_32) ::1       INPUT   number of additional dispers.
!     NOVELO  INTEGER(kind=int_32) ::1       INPUT   number of additional velos.
!     IDPNT   INTEGER(kind=int_32) ::NOSYS     INPUT   pointer systems to dispersions
!     IVPNT   INTEGER(kind=int_32) ::NOSYS     INPUT   pointer systems to velocities
!     DERIV   REAL(kind=sp) ::NOTOT*NOSEG  OUTPUT  derivatives
!     AMAT    REAL(kind=sp) ::large     IN/OUT  matrix to be updated
!     JTRACK  INTEGER(kind=int_32) ::1       INPUT   number of codiagonals of AMAT
!     IOPT    INTEGER(kind=int_32) ::1       INPUT   = 0 or 2 DISP at zero flow
!                                         = 1 or 3 no DISP at zero flow
!                                         = 0 or 1 DISP over boundary
!                                         = 2 or 3 no DISP over boundary
!     ILFLAG  INTEGER(kind=int_32) ::1       INPUT   if 0 then 3 length values
!
      use timers

      real(kind=sp) ::DISP  (  3) , DISPER(*) , AREA (*) , FLOW (*) ,
     *           ALENG (  *) , VELO  (*) , BOUND(*) , AMAT (*) ,
     *           DERIV(*)
      integer(kind=int_32) ::IPOINT(4,*) , IDPNT(*)  , IVPNT(*)

      integer(kind=int_32) ::iband, iq, i, it, i3, i4, iopt, ilflag, isys
      integer(kind=int_32) ::j, jtrack, jt
      integer(kind=int_32) ::noq, noq1, noq2, nodisp, novelo, notot, nsys
      integer(kind=int_32) ::kt, k1, k2

      real(kind=sp) ::a, q, q1, q2, e, al, f1, f2, dl

      integer(kind=int_32) ::ithandl = 0
      if ( timon ) call timstrt ( "dlwq70", ithandl )
!
      IBAND = 2*JTRACK + 1
      DO 50 IQ = 1 , NOQ
!
!         initialisations , check for transport anyhow
!
      I    = IPOINT(1,IQ)
      J    = IPOINT(2,IQ)
      IF ( I .EQ. 0 .OR. J .EQ. 0 ) GOTO 50
      A    = AREA(IQ)
      Q    = FLOW(IQ)
      IF ( MOD(IOPT,2) .EQ. 1 .AND. ABS(Q) .LT. 10.0E-25 ) GOTO 50
           IF ( A .LT. 1.0E-25 )  A = 1.0
      E  = DISP(1)
      AL = ALENG(1)
      IF ( IQ .GT. NOQ1      ) THEN
           E  = DISP (2)
           AL = ALENG(2)
      ENDIF
      IF ( IQ .GT. NOQ1+NOQ2 ) THEN
           E  = DISP (3)
           AL = ALENG(3)
      ENDIF
      IF ( ILFLAG .EQ. 1 ) THEN
           DL = A/(ALENG(2*IQ-1) + ALENG(2*IQ))
           F1 = ALENG(2*IQ  )*DL/A
           F2 = ALENG(2*IQ-1)*DL/A
      ELSE
           DL = A/AL
           F1 = 0.5
           F2 = 0.5
      ENDIF
      E  = E*DL
      IF (IDPNT(ISYS).GT.0) E = E + DISPER((IQ-1)*NODISP+IDPNT(ISYS))*DL
      IF (IVPNT(ISYS).GT.0) Q = Q + VELO  ((IQ-1)*NOVELO+IVPNT(ISYS))*A
      Q1 = F1*Q
      Q2 = F2*Q
      IF ( I .LT. 0 ) GOTO 10
      IF ( J .LT. 0 ) GOTO 30
!
!        the regular case
!
      JT = (I-1)*IBAND + JTRACK + 1
      KT = JT + (J-I)
      AMAT(JT) = AMAT(JT) + Q1 + E
      AMAT(KT) = AMAT(KT) + Q2 - E
      IT = (J-1)*IBAND + JTRACK + 1
      KT = IT + (I-J)
      AMAT(IT) = AMAT(IT) - Q2 + E
      AMAT(KT) = AMAT(KT) - Q1 - E
      GOTO 50
!
!        The 'from' segment is a boundary
!
   10 IF ( J    .LT. 0 ) GOTO 50
      IF ( MOD(IOPT,4) .GT. 1 ) E = 0.0
      IF ( MOD(IOPT,8) .GE. 4 ) THEN
           IF ( Q .GT. 0.0 ) THEN
                Q1 = Q
                Q2 = 0.0
           ELSE
                Q1 = 0.0
                Q2 = Q
           ENDIF
      ENDIF
      K1 = (-I-1)*NOTOT
      I4 = ( J-1)*NSYS  + 1
      DO 20 I3=ISYS,ISYS+NSYS-1
      DERIV(I4) = DERIV(I4) + ( Q1+E) * BOUND(K1+I3)
   20 I4=I4+1
      IT = (J-1)*IBAND + JTRACK + 1
      AMAT(IT) = AMAT(IT) - Q2 + E
      GOTO 50
!
!        The 'to' element was a boundary.
!
   30 IF ( MOD(IOPT,4) .GT. 1 ) E = 0.0
      IF ( MOD(IOPT,8) .GE. 4 ) THEN
           IF ( Q .GT. 0.0 ) THEN
                Q1 = Q
                Q2 = 0.0
           ELSE
                Q1 = 0.0
                Q2 = Q
           ENDIF
      ENDIF
      K2 = (-J-1)*NOTOT
      I4 = ( I-1)*NSYS  + 1
      DO 40 I3=ISYS,ISYS+NSYS-1
      DERIV(I4) = DERIV(I4) + (-Q2+E) * BOUND(K2+I3)
   40 I4=I4+1
      JT = (I-1)*IBAND + JTRACK + 1
      AMAT(JT) = AMAT(JT) + Q1 + E
!
!        end of the loop over exchanges
!
   50 CONTINUE
!
      if ( timon ) call timstop ( ithandl )
      RETURN
      END

      end module m_dlwq70
