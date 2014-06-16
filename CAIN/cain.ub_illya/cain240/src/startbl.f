	SUBROUTINE STARTBL(NAM,ID,APERT)
C  Start beamline specifiend by the name.
C    ID>0   beamline ID
C      <0   error
	USE BEAMLN
	IMPLICIT NONE
	INTEGER ID
	CHARACTER(*) NAM
	REAL(8) APERT(2)
      INCLUDE 'include/ctrlcm.h'
	INCLUDE 'include/blexpand.h'
	INTEGER I

	IF(NBEAMLINE.LE.0) THEN
	  ID=-1
	  RETURN
	ENDIF
	BLID(1)=0
	DO I=1,NBEAMLINE
	  IF(NAM.EQ.BL(I)%NAME) THEN
	    BLID(1)=I
	    EXIT
	  ENDIF
	ENDDO
	IF(BLID(1).EQ.0) THEN
	  ID=-2
	  IF(MSGLVL.GE.0) THEN
	    WRITE(MSGFL,120) NAM
120       FORMAT(' (SUBR.STARTBL) Beamline ',A,' undefined.')
        ENDIF
	  RETURN
	ENDIF
	LVL=1
	ELN(1)=0
	DIR(1)=1
	APERT=BL(BLID(1))%APERT
	APERT0(1:2,LVL)=APERT
	ID=BLID(1)
	RETURN
	END