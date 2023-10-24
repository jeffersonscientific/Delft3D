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
      module m_dlwqb5
      use m_waq_type_definitions


      implicit none

      contains


      SUBROUTINE DLWQB5 ( DISP   , DISPER , AREA   , FLOW   , ALENG  ,
     *                    VELO   , CONC   , BOUND  , IPOINT , NOSYS  ,
     *                    NOTOT  , NOQ1   , NOQ2   , NOQ    , NODISP ,
     *                    NOVELO , IDPNT  , IVPNT  , IOPT   , AMASS2 ,
     *                    ILFLAG , DMPQ   , NDMPQ  , IDT    , IQDMP  )
!
!     Deltares     SECTOR WATERRESOURCES AND ENVIRONMENT
!
!     CREATED             : march 1988 by L.Postma
!                           may 1992 by J.v.Gils
!
!     FUNCTION            : Makes a mass balance final to
!                           implicit integration methods.
!                           Identical to DLWQ64, but with dimension
!                           BOUND(NOSYS,*), with multiplication
!                           factor IDT on mass balances and with
!                           loop over active substances only
!                           (loops 10, 30 and 50).
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
!     CONC    REAL(kind=sp) ::NOTOT*NOSEG  INPUT   concentrations
!     BOUND   REAL(kind=sp) ::NOSYS*?    INPUT   boundary concentrations
!     IPOINT  INTEGER(kind=int_32) ::4*NOQ     INPUT   exchange pointers
!     NOSYS   INTEGER(kind=int_32) ::1       INPUT   number  of active substances
!     NOTOT   INTEGER(kind=int_32) ::1       INPUT   number  of total substances
!     NOQ1    INTEGER(kind=int_32) ::1       INPUT   nr of exchanges in first dir.
!     NOQ2    INTEGER(kind=int_32) ::1       INPUT   nr of exchanges in second dir.
!     NOQ3    INTEGER(kind=int_32) ::1       INPUT   nr of exchanges in third dir.
!     NOQ     INTEGER(kind=int_32) ::1       INPUT   total number of exchanges
!     NODISP  INTEGER(kind=int_32) ::1       INPUT   number  of additional dispers.
!     NOVELO  INTEGER(kind=int_32) ::1       INPUT   number  of additional velos.
!     IDPNT   INTEGER(kind=int_32) ::NOSYS     INPUT   pointer systems to dispersions
!     IVPNT   INTEGER(kind=int_32) ::NOSYS     INPUT   pointer systems to velocities
!     IOPT    INTEGER(kind=int_32) ::1       INPUT   = 0 or 2 DISP at zero flow
!                                         = 1 or 3 no DISP at zero flow
!                                         = 0 or 1 DISP over boundary
!                                         = 2 or 3 no DISP over boundary
!     AMASS2  REAL(kind=sp) ::NOTOT*5    IN/OUT  mass balance array
!     ILFLAG  INTEGER(kind=int_32) ::1       INPUT   if 0 then 3 length values
!     DMPQ    REAL(kind=sp) ::NOTOT*NDMPQ*? IN/OUT  mass balance dumped exchange
!                                         if INTOPT > 7
!     NDMPQ   INTEGER(kind=int_32) ::1       INPUT   number of dumped exchanges
!     IDT     INTEGER(kind=int_32) ::1       INPUT   timestep (or 1 for steady state)
!     IQDMP   INTEGER(kind=int_32) ::*       INPUT   pointer dumped exchanges
!
      use timers

      INTEGER(kind=int_32) ::NDMPQ
      INTEGER(kind=int_32) ::IQDMP   (*) , IPOINT(4,*) , IDPNT(*) , IVPNT(*)
      real(kind=sp) ::DISP  (  3) , DISPER(*)   , AREA (*) , FLOW  (*) ,
     *           ALENG (  *) , VELO  (*)   , CONC (*) , BOUND (*) ,
     *           AMASS2(*)   , DMPQ(*)
      integer(kind=int_32) ::NOSYS, NOTOT, NOQ1, NOQ2, NOQ, NODISP, NOVELO
      integer(kind=int_32) ::IOPT, ILFLAG, IDT

      integer(kind=int_32) ::i, i3, i4, i5, i6, iq, ipq, is
      integer(kind=int_32) ::noq12, ibflag, j, ipb, k1, k2

      real(kind=sp) ::a, q, e, al, dl, d, v, dq

      integer(kind=int_32) ::ithandl = 0
      if ( timon ) call timstrt ( "dlwqb5", ithandl )
!
!         loop accross the number of exchanges
!
      I4 = 3*NOTOT
      I5 = 4*NOTOT
      I6 = NOSYS*NDMPQ
      NOQ12 = NOQ1 + NOQ2
      IF ( MOD(IOPT,16) .GE. 8  ) THEN
         IBFLAG = 1
      ELSE
         IBFLAG = 0
      ENDIF
!
      DO 60 IQ = 1 , NOQ
!
!         initialistations, check for transport anyhow
!
      I    = IPOINT(1,IQ)
      J    = IPOINT(2,IQ)
      IF ( I .EQ. 0 .OR. J .EQ. 0 ) GOTO 60
!
!     Check if exchange is dump exchange, set IPB
!
      IF ( IBFLAG .EQ. 1 ) THEN
         IF ( IQDMP(IQ) .GT. 0 ) THEN
            IPB = IQDMP(IQ)
            IPQ = (IQDMP(IQ)-1)*NOSYS
         ELSE
            IPB = 0
         ENDIF
      ELSE
         IPB = 0
      ENDIF
      IF ( I .GT. 0 .AND. J .GT. 0 .AND. IPB .EQ. 0 ) GOTO 60
      A    = AREA(IQ)
      Q    = FLOW(IQ)
      IF ( MOD(IOPT,2) .EQ. 1 .AND. IQ .LE. NOQ12) THEN
           IF ( ABS(Q) .LT. 10.0E-25 )  GOTO 60
      ENDIF
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
           DL = A / (ALENG(2*IQ-1) + ALENG(2*IQ))
      ELSE
           DL = A / AL
      ENDIF
      E  = E*DL
      IF ( I .LT. 0 ) GOTO 20
      IF ( J .LT. 0 ) GOTO 40
!
!         The regular case
!
      K1 = (I-1)*NOTOT
      K2 = (J-1)*NOTOT
      DO 10 I3=1,NOSYS
         IS = MIN ( I3 , NOSYS )
!
!        dispersion
!
         IF ( IDPNT(IS) .GT. 0 ) THEN
            D = E + DISPER((IQ-1)*NODISP+IDPNT(IS))*DL
         ELSE
            D = E
         ENDIF
!
!        flow
!
         IF ( IVPNT(IS) .GT. 0 ) THEN
            V = Q + VELO  ((IQ-1)*NOVELO+IVPNT(IS))*A
         ELSE
            V = Q
         ENDIF
!
!        transport
!
         IF ( V .GT. 0.0 ) THEN
            DQ = ((V+D)*CONC(K1+I3) - D*CONC(K2+I3))*IDT
         ELSE
            DQ = ((V-D)*CONC(K2+I3) + D*CONC(K1+I3))*IDT
         ENDIF
!
!        mass balance
!
         IF ( DQ .GT. 0.0 ) THEN
            DMPQ(IPQ+I3)=DMPQ(IPQ+I3) + DQ
         ELSE
            DMPQ(IPQ+I3+I6)=DMPQ(IPQ+I3+I6) - DQ
         ENDIF
!
   10 CONTINUE
      GOTO 60
!
!        The 'from' element was a boundary. Note the 2 options.
!
   20 IF ( J .LT. 0 ) GOTO 60
      K1 = (-I-1)*NOSYS
      K2 = ( J-1)*NOTOT
      DO 30 I3=1,NOSYS
      IS = MIN ( I3 , NOSYS )
      V  = Q
      D  = 0.0
      IF ( IVPNT(IS) .GT. 0 ) V = V + VELO  ((IQ-1)*NOVELO+IVPNT(IS))*A
      IF ( MOD(IOPT,4) .LT.  2 ) THEN
           D  = E
           IF ( IDPNT(IS).GT.0 ) D=D+ DISPER((IQ-1)*NODISP+IDPNT(IS))*DL
      ENDIF
      IF ( V .GT. 0.0 ) THEN
           DQ = ((V+D)*BOUND(K1+I3) - D*CONC (K2+I3))*IDT
      ELSE
           DQ = ((V-D)*CONC (K2+I3) + D*BOUND(K1+I3))*IDT
      ENDIF
      IF ( DQ .GT. 0.0 ) THEN
           AMASS2(I3+I4) = AMASS2(I3+I4) + DQ
      ELSE
           AMASS2(I3+I5) = AMASS2(I3+I5) - DQ
      ENDIF
      IF ( IPB .GT. 0 ) THEN
         IF ( DQ .GT. 0.0 ) THEN
            DMPQ(IPQ+I3)=DMPQ(IPQ+I3) + DQ
         ELSE
            DMPQ(IPQ+I3+I6)=DMPQ(IPQ+I3+I6) - DQ
         ENDIF
      ENDIF
   30 CONTINUE
      GOTO 60
!
!        The 'to' element was a boundary.
!
   40 K1 = ( I-1)*NOTOT
      K2 = (-J-1)*NOSYS
      DO 50 I3=1,NOSYS
      IS = MIN ( I3 , NOSYS )
      V  = Q
      D  = 0.0
      IF ( IVPNT(IS) .GT. 0 ) V = V + VELO  ((IQ-1)*NOVELO+IVPNT(IS))*A
      IF ( MOD(IOPT,4)  .LT.  2 ) THEN
           D  = E
           IF ( IDPNT(IS).GT.0 ) D=D+ DISPER((IQ-1)*NODISP+IDPNT(IS))*DL
      ENDIF
      IF ( V .GT. 0.0 ) THEN
           DQ = ((V+D)*CONC (K1+I3) - D*BOUND(K2+I3))*IDT
      ELSE
           DQ = ((V-D)*BOUND(K2+I3) + D*CONC (K1+I3))*IDT
      ENDIF
      IF ( DQ .GT. 0.0 ) THEN
           AMASS2(I3+I5) = AMASS2(I3+I5) + DQ
      ELSE
           AMASS2(I3+I4) = AMASS2(I3+I4) - DQ
      ENDIF
      IF ( IPB .GT. 0 ) THEN
         IF ( DQ .GT. 0.0 ) THEN
            DMPQ(IPQ+I3)=DMPQ(IPQ+I3) + DQ
         ELSE
            DMPQ(IPQ+I3+I6)=DMPQ(IPQ+I3+I6) - DQ
         ENDIF
      ENDIF
   50 CONTINUE
!
!        end of the loop over exchanges
!
   60 CONTINUE
!
      if ( timon ) call timstop ( ithandl )
      RETURN
      END

      end module m_dlwqb5
