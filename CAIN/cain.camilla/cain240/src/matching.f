      SUBROUTINE MATCHING(MSGLVL,MSGFL,IRTN)
	USE FLCHTYP
c	USE BEAMLN
	USE ARRAYMOD
	USE MATCHMOD
	IMPLICIT NONE
	INTEGER MSGLVL,MSGFL,IRTN
	INCLUDE 'include/nameleng.h'
	INCLUDE 'include/evparc.h'
	INTEGER MSGLVL1

	INTEGER I,ITMAX,IT,ISTAT,IRTN1
	REAL(8) OBJFUN,DDER
	LOGICAL SCALE/.FALSE./,GRAD/.FALSE./
	REAL(8) BAI(2,2),BAF(2,2),EPI(4),EPF(4)
	TYPE(FLCHTYPE) FC
	CHARACTER(80) ERR
	EXTERNAL MATCHFUN
	REAL(8), ALLOCATABLE:: CONFUN(:)

	ALLOCATE(CONFUN(NCOND),STAT=ISTAT)
	IF(ISTAT.NE.0) GOTO 900
	DO I=1,NVAR
	  IF(KVAR(1,I).GE.0) THEN
	    XVAR(I)=VPAR(KVAR(1,I))
	  ELSE
	    XVAR(I)=ARR(-KVAR(1,I))%VAL(KVAR(2,I))
	  ENDIF
	  XVAR0(I)=XVAR(I)
	ENDDO
	IF(MSGLVL.GE.1) CALL PRMATCHCOND(0,MSGFL)
	DDER=0.001D0
	MSGLVL1=0
	CALL SQP(MATCHFUN,NVAR,NCOND-NPOS,NCOND,XVAR,OBJFUN,CONFUN,
     %       SCALE,GRAD,DDER,ITMAX,IT,MSGLVL1,MSGFL,IRTN1)
	IF(IRTN1.NE.0) THEN
	  IF(MSGLVL.GE.0) THEN
C              repeat for error message
          WRITE(MSGFL,200)
200       FORMAT(' +++ Matching failed. Try once more for message +++')
	    DO I=1,NVAR
	      XVAR(I)=XVAR0(I)
	    ENDDO
	    MSGLVL1=1
	    CALL SQP(MATCHFUN,NVAR,NCOND-NPOS,NCOND,XVAR,OBJFUN,CONFUN,
     %         SCALE,GRAD,DDER,ITMAX,IT,MSGLVL1,MSGFL,IRTN1)
	  ENDIF
	  CALL EVDEFP('Convergence',1D10,IRTN1)
C          Failure of matching does not cause CAIN abnormal termination.
C          User must check the 'Convergence' variables to stop or to continue.
C          It is set to a big value here.
	ELSE
	  CALL MATCHFUN(NVAR,NCOND,XVAR,OBJFUN,CONFUN)
	  IF(MSGLVL.GE.1) CALL PRMATCHCOND(1,MSGFL)
	  CALL EVDEFP('Convergence',CONV,IRTN1)
	ENDIF
	IRTN=0
	GOTO 1000
900   IRTN=1000
	WRITE(MSGFL,910)
910   FORMAT(' (SUBR.MATCHING) Allocation error.')
      GOTO 1000
940   IRTN=100
	WRITE(MSGFL,945)
945   FORMAT(' (SUBR.MATCHING) Matching failed.')
      GOTO 1000
1000  DEALLOCATE(CONFUN,STAT=ISTAT)
	RETURN
	END

	SUBROUTINE PRMATCHCOND(KK,IFL)
	USE FLCHTYP
	USE BEAMLN
	USE ARRAYMOD
	USE MATCHMOD
	IMPLICIT NONE
	INTEGER KK,IFL
	INCLUDE 'include/nameleng.h'
	INCLUDE 'include/evparc.h'
	CHARACTER(3) CTWIS(8)/'Bx','By','Ax','Ay','Ex','Ey','Epx','Epy'/
	CHARACTER(1) CH
	INTEGER I,NC,IPR
	REAL(8) X1
	CHARACTER(30) FMT

	IF(KK.EQ.0) THEN
	  WRITE(IFL,200) 'condition',BLNAM(1:NCBLNAM)
	ELSE
	  WRITE(IFL,200) 'results',BLNAM(1:NCBLNAM)
	ENDIF
200   FORMAT(' +++ Matching ',A,' for beamline ',A,' +++')
	IF(LPER.NE.0) THEN
	  IF(KK.EQ.0) WRITE(IFL,220)
220     FORMAT(' Periodic')
      ELSE
	  IF(KK.EQ.0) THEN
	    IPR=1
	  ELSE
	    IPR=0
	    DO I=1,8
	      IF(TWISSIN(I)%L.EQ.2) THEN
	        IPR=1
	        EXIT
	      ENDIF
	    ENDDO
	  ENDIF
	  IF(IPR.NE.0) WRITE(IFL,240)
240     FORMAT(' Entrance parameters')
	  DO I=1,8
	    IF(TWISSIN(I)%L.EQ.1) THEN
	      IF(KK.EQ.0) WRITE(IFL,260) CTWIS(I),TWISSIN(I)%X
260         FORMAT(5X,A,1PD13.5)
          ELSE
	      WRITE(IFL,280) CTWIS(I),TWISSIN(I)%X,
     %        GSTRMATCH(TWISSIN(I)%C(1):TWISSIN(I)%C(2))
280         FORMAT(5X,A,1PD13.5,'  "',A,'"')
          ENDIF
	  ENDDO
	ENDIF
	WRITE(IFL,300) NVAR
300   FORMAT(' Variables   number of variables=',I3)
      DO I=1,NVAR
	  WRITE(IFL,'(3X,A,1PD16.7)') VARNAM(I)(1:NCVARNAMMAX),XVAR(I)
	ENDDO
	WRITE(IFL,320) NCOND-NPOS,NPOS
320   FORMAT(' Conditions:  Zero:',I3,'  Positive:',I3)
	WRITE(FMT,340) NCCONDMAX+3+4+3
340   FORMAT("(3X,A,' ',A,' 0',T",I3,",1PD12.5)")
	DO I=1,NCOND
	  CH='='
	  IF(I.GT.NCOND-NPOS) CH='>'
	  WRITE(IFL,FMT) GSTRMATCH(COND(I)%C(1):COND(I)%C(2)),CH,
     %    COND(I)%X
      ENDDO
	IF(KK.NE.0) WRITE(IFL,360) CONV
360   FORMAT(' Convergence ',1PD14.6)
	RETURN
	END