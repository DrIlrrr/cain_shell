      SUBROUTINE COPYBM(N1,N2)
	USE BEAMCM
      IMPLICIT NONE
      INTEGER N1,N2
C      INCLUDE 'include/beamcm.h'
      INTEGER I
      IF(NP.EQ.0.OR.N1.LE.0.OR.N1.GT.NP.OR.N2.LE.0.OR.N2.GT.NP
     %     .OR.N1.EQ.N2) RETURN
	LOST(N2)=LOST(N1)
      PNAME(N2)=PNAME(N1)
      WGT(N2)=WGT(N1)
      KIND(N2)=KIND(N1)
      GEN(N2)=GEN(N1)
      ISBIN(N2)=ISBIN(N1)
      IADRS(1,N2)=IADRS(1,N1)
      IADRS(2,N2)=IADRS(2,N1)
      LBBFIN(N2)=LBBFIN(N1)
      DO 200 I=0,3
        TXYS(I,N2)=TXYS(I,N1)
        EP(I,N2)=EP(I,N1)
 200  CONTINUE
      DO 220 I=1,3
        SPIN(I,N2)=SPIN(I,N1)
        FLD(I,1,N2)=FLD(I,1,N1)
        FLD(I,2,N2)=FLD(I,2,N1)
 220  CONTINUE
      RETURN
      END
