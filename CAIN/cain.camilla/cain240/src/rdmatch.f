      SUBROUTINE RDMATCH(LN,LINE,NCHAR,IRTN)
	USE FLCHTYP
	USE READMOD
	USE BEAMLN
	USE ARRAYMOD
      IMPLICIT NONE
      INTEGER LN(2,3),NCHAR(*),IRTN
      CHARACTER*(*) LINE(*)
      INTEGER MOP
      PARAMETER (MOP=10)
      CHARACTER*13 OP(MOP)/'BEAMLINE','PERIODIC',
     %   'BETA','ALPHA','ETA','ETAPRIME',   ! must not change the order of these four
     %   'VARIABLE','ZERO','POSITIVE','MINFUN'/
      INTEGER NFF(MOP)/1,0,2,2,2,2,1,1,1,1/
      INTEGER IDBMLN,IDBETA
	INTEGER IPPER,IPVAR,IPZERO,IPPOS
      INTEGER IOP(MOP)
      INCLUDE 'include/ctrlcm.h'
C      INCLUDE 'include/readcm.h'
      INTEGER J,I,NF,NC,K,NCTXT2(0:3)
	INTEGER LPER,IDBL
	TYPE(FLCHTYPE) TWISSIN1(8),FC

	REAL(8) CONV,APERT1(2)
	INTEGER NCBLNAM
      CHARACTER(16) BLNAM
	CHARACTER(MCTEXT) TEXT2(0:3)
C
      IRTN=0
      IF(LN(1,2).EQ.0) RETURN
      CALL CMDBLK('MATCHING',MOP,OP,0,NFF,MBL,NBL,KOP,KBL,REL,LNKW,LNBL,
     %    LN,LINE,NCHAR,IRTN)
      IF(IRTN.GT.1) GOTO 990
      IF(IRTN.EQ.1.OR.NBL.LE.0) RETURN
      DO I=1,MOP
        IOP(I)=0
      ENDDO
      J=1
      DO I=1,MOP
        ID(I)=J
        IF(OP(I).EQ.'BEAMLINE') IDBMLN=ID(I)
	  IF(OP(I).EQ.'BETA') IDBETA=ID(I)
	  IF(OP(I).EQ.'PERIODIC') IPPER=I
	  IF(OP(I).EQ.'VARIABLE') IPVAR=I
	  IF(OP(I).EQ.'ZERO') IPZERO=I
	  IF(OP(I).EQ.'POSITIVE') IPPOS=I
        J=J+MAX(0,NFF(I))
      ENDDO
      DO I=1,J-1
        PAR(I)=UNDEF
      ENDDO
	NGSTRRD=0
	NCTXT2=0
C
      DO 400 J=1,NBL
        I=KBL(J)
	  IOP(I)=1
        IF(NFF(I).EQ.0) THEN
        ELSEIF(NFF(I).GE.1) THEN
          CALL BLKREC(LNBL(1,1,J),LINE,NCHAR,TEXT,NC,IRTN,MSGFL)
          IF(IRTN.NE.0) GOTO 990
          IF(OP(I).EQ.'VARIABLE'.OR.OP(I).EQ.'ZERO'
     %       .OR.OP(I).EQ.'POSITIVE'.OR.OP(I).EQ.'MINFUN') THEN
	      IF(OP(I).EQ.'VARIABLE') THEN
	        K=0
	      ELSEIF(OP(I).EQ.'ZERO') THEN
	        K=1
	      ELSEIF(OP(I).EQ.'POSITIVE') THEN
	        K=2
	      ELSEIF(OP(I).EQ.'MINFUN') THEN
	        K=3
	      ENDIF
	      NCTXT2(K)=NC
	      TEXT2(K)=TEXT(1:NC)
          ELSE
            DO K=1,NFF(I)
	        FFF(K)=UNDEF
            ENDDO
            CALL BLKRHS(TEXT(1:NC),FFF,MFFF,NF,GSTRRD,NGSTRRD,IRTN)
            IF(IRTN.NE.0) GOTO 990
            IF(NF.GT.NFF(I)) GOTO 910
            IF(NF.GE.1) THEN
              DO K=1,NF
                IF(FFF(K).NE.UNDEF) THEN
                  PAR(ID(I)+K-1)=FFF(K)
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDIF
 400  CONTINUE
 	IF(PAR(IDBMLN).EQ.UNDEF) GOTO 920
	IF(PAR(IDBMLN)%L.NE.2) GOTO 920
	BLNAM=GSTRRD(PAR(IDBMLN)%C(1):PAR(IDBMLN)%C(2))
	NCBLNAM=PAR(IDBMLN)%C(2)-PAR(IDBMLN)%C(1)+1
	CALL STARTBL(BLNAM(1:NCBLNAM),IDBL,APERT1)
	IF(IDBL.LE.0) GOTO 925
	LPER=IOP(IPPER)
	IF(LPER.EQ.0) THEN
	  DO I=1,8
	    IF(PAR(IDBETA+I-1).EQ.UNDEF) THEN
	      TWISSIN1(I)=ZERO
	    ELSE
	      TWISSIN1(I)=PAR(IDBETA+I-1)
	    ENDIF
	  ENDDO
	ENDIF
	
	CALL RDMATCH0(IDBL,BLNAM(1:NCBLNAM),LPER,TWISSIN1,TEXT2,NCTXT2,
     %   CONV,GSTRRD,NGSTRRD,IRTN)
	
	IF(IRTN.NE.0) GOTO 1000

 800  IRTN=0
      GOTO 1000
C---
C 900  IRTN=1000
C      WRITE(MSGFL,905) TEXT(1:NC)
C 905  FORMAT(' (SUBR.RDMATCH) Invalid operand "',A,'".')
C      GOTO 1000
 910  IRTN=1010
      WRITE(MSGFL,915) OP(I)
 915  FORMAT(' (SUBR.RDMATCH) Too many numbers for ',
     %  'operand "',A,'".')
      GOTO 1000
 920  IRTN=1020
      WRITE(MSGFL,921)
 921  FORMAT(' (SUBR.RDMATCH) Beamline name not specified. ')
      GOTO 1000
 925  IRTN=1025
      WRITE(MSGFL,926) BLNAM(1:NCBLNAM)
 926  FORMAT(' (SUBR.RDMATCH) Beamline "',A,'" does not exist.')
      GOTO 1000
 990  IRTN=1090
      GOTO 1000
 1000 RETURN
      END

	SUBROUTINE RDMATCH0(IDBL0,BLNAM0,LPER0,TWISSIN1,TXT,NCTXT,CONV0,
     %   GSTRRD,NGSTRRD,IRTN)
	USE FLCHTYP
	USE BEAMLN
	USE ARRAYMOD
	USE MATCHMOD
	IMPLICIT NONE
      INCLUDE 'include/nameleng.h'
      INCLUDE 'include/evparc.h'
	INCLUDE 'include/ctrlcm.h'
	INTEGER IDBL0,LPER0,NCTXT(0:3),NGSTRRD,IRTN
	TYPE(FLCHTYPE) TWISSIN1(8)
	CHARACTER(*) BLNAM0,TXT(0:3),GSTRRD
	REAL(8) CONV0
	CHARACTER(8) KWD(0:2)/'VARIABLE','ZERO','POSITIVE'/
	
	INTEGER, PARAMETER:: MITEM=100  !  max number of var, zero, pos
	INTEGER NITEM,ITEM(2,MITEM)
      INTEGER, PARAMETER:: MIND=50
	INTEGER NNIND(MIND),IND(MIND),NIND

	CHARACTER(16) NAM
	INTEGER I,J,K,ISTAT,N,N1,II,KK,K1,IDARR,NHIT,NCNAM,NGSTRMATCH0,
     %    IREP,KVAR2
	REAL(8) BAI(2,2),BAF(2,2),EPI(4),EPF(4)
	TYPE(FLCHTYPE) FC
	CHARACTER(80) ERR
	INTEGER MF,NF,NP
	INTEGER, PARAMETER:: MP=200
	CHARACTER(16) NAMP(MP),NAMF(1)

      LPER=LPER0
	IDBLMATCH=IDBL0
	BLNAM=BLNAM0
	NCBLNAM=LEN(BLNAM0)
	CONV=CONV0
	NGSTRMATCH=0
C--- Variables
	CALL BLKRHS0(TXT(0)(1:NCTXT(0)),MITEM,NITEM,ITEM,IRTN)
	IF(IRTN.NE.0) GOTO 990
	I=0
	DO IREP=1,2
C         repeat twice. First to count the number of variables.
        NVAR=0
	  DO K=1,NITEM
	    CALL EVLHS(TXT(I)(ITEM(1,K):ITEM(2,K)),KK,NAM,NCNAM,IDARR,
     %              MIND,NIND,NNIND,IND,ERR)
	    IF(ERR.NE.' ') GOTO 914
	    IF(KK.EQ.1) THEN    !  scalar
	      NVAR=NVAR+1
	      IF(IREP.EQ.2) KVAR(1,NVAR)=IDARR
	    ELSEIF(KK.EQ.2) THEN   !  array element
	      IF(NIND.EQ.0) THEN
	        N=1
	        DO N1=1,ARR(IDARR)%RANK
	          N=N*ARR(IDARR)%DIM(3,N1)
	        ENDDO
	        IF(IREP.EQ.2) THEN
	          DO N1=1,N
	            KVAR(1,NVAR+N1)=-IDARR
	            KVAR(2,NVAR+N1)=N1
	          ENDDO
	        ENDIF
	        NVAR=NVAR+N
	      ELSE
	        IF(NIND.NE.ARR(IDARR)%RANK) GOTO 916
	        DO K1=1,NIND
	          IF(NNIND(K1).NE.1) GOTO 916
	        ENDDO
	        CALL ARRIND2N(IND,NIND,ARR(IDARR)%DIM,KVAR2,IRTN)
	        IF(IRTN.NE.0) GOTO 916
	        NVAR=NVAR+1
	        IF(IREP.EQ.2) THEN
					  KVAR(1,NVAR)=-IDARR
	          KVAR(2,NVAR)=KVAR2
	        ENDIF
	      ENDIF
	    ELSEIF(KK.EQ.3) THEN
	      GOTO 918  !  character var
	    ELSE
	      GOTO 920  !  undefined var
	    ENDIF
	  ENDDO
	  IF(IREP.EQ.1) THEN
	    IF(NVAR.LE.0) GOTO 910
	    ALLOCATE(KVAR(2,NVAR),LVAR(NVAR),XVAR(NVAR),XVAR0(NVAR),
     %       VARNAM(NVAR),NCVARNAM(NVAR),STAT=ISTAT)
	    IF(ISTAT.NE.0) GOTO 900
	  ENDIF
	ENDDO
C  Store Twiss at entrance
	IF(LPER.EQ.0) THEN
	  TWISSIN=TWISSIN1
	  DO I=1,8
	    IF(TWISSIN1(I)%L.EQ.2) THEN
	      CALL FLCHSET3(GSTRRD(TWISSIN1(I)%C(1):TWISSIN1(I)%C(2)),
     %        TWISSIN(I),GSTRMATCH,NGSTRMATCH,ERR)
	      IF(ERR.NE.' ') GOTO 965
	    ENDIF
	  ENDDO
	ENDIF
C  Check Twiss at entrance
	CALL MATCHSET(IRTN)
	IF(IRTN.NE.0) GOTO 990
C       This will cause a stop if the initial values of variables
C       are not good (no stable periodic solution, etc)
C       Call of BLOPTICS is needed because Twiss functions will
C       be called in checking the matching conditions.
C       In the future, an error of BLOPTICS at this stage should
C       be ignored.

C  Count the number of conditions
	NCOND=0
	DO I=1,2
	  IF(NCTXT(I).LE.0) CYCLE
	  CALL BLKRHS0(TXT(I)(1:NCTXT(I)),MITEM,NITEM,ITEM,IRTN)
	  IF(IRTN.NE.0) GOTO 990
	  NCOND=NCOND+NITEM
	ENDDO
	IF(NCOND.EQ.0) GOTO 935
	ALLOCATE(COND(NCOND),STAT=ISTAT)
	IF(ISTAT.NE.0) GOTO 900
C  Check conditions
	NCOND=0
	NPOS=0
	DO I=1,2
	  IF(NCTXT(I).LE.0) CYCLE
	  CALL BLKRHS0(TXT(I)(1:NCTXT(I)),MITEM,NITEM,ITEM,IRTN)
	  DO K=1,NITEM
	    NCOND=NCOND+1
	    IF(I.EQ.2) NPOS=NPOS+1
	    CALL EVAL0(TXT(I)(ITEM(1,K):ITEM(2,K)),FC,ERR)
	    IF(ERR.NE.' ') GOTO 914
C        a little complicated to allow expression both enclosed
C        and non-enclosed by apostrophes
	    IF(FC%L.EQ.1) THEN
	      CALL FLCHSET3(TXT(I)(ITEM(1,K):ITEM(2,K)),COND(NCOND),
     %        GSTRMATCH,NGSTRMATCH,ERR)
	      IF(ERR.NE.' ') GOTO 965
	      COND(NCOND)%X=FC%X
	    ELSEIF(FC%L.EQ.2) THEN
	      NGSTRMATCH0=NGSTRMATCH
C               save the string temporarily in GSTRMATCH
	      CALL FLCHSET3(GSTR2(EVALLAST)(FC%C(1):FC%C(2)),FC,
     %        GSTRMATCH,NGSTRMATCH0,ERR)
	      IF(ERR.NE.' ') GOTO 960
	      CALL EVAL0(GSTRMATCH(FC%C(1):FC%C(2)),FC,ERR)
	      IF(ERR.NE.' ') GOTO 950
	      IF(FC%L.NE.1) GOTO 950
	      CALL FLCHSET3(GSTR2(EVALLAST)(FC%C(1):FC%C(2)),COND(NCOND),
     %        GSTRMATCH,NGSTRMATCH,ERR)
	      IF(ERR.NE.' ') GOTO 960
	      COND(NCOND)%X=FC%X
	    ELSE
	      GOTO 914
	    ENDIF
	  ENDDO
	ENDDO
C  Store MINFUN
      IF(NCTXT(3).EQ.0) THEN
	  MINFUN%L=1
	  MINFUN%X=0
	ELSE
	  CALL EVAL0(TXT(3)(1:NCTXT(3)),FC,ERR)
	  IF(ERR.NE.' ') GOTO 912
	  IF(FC%L.EQ.1) THEN
	    CALL FLCHSET3(TXT(3)(1:NCTXT(3)),MINFUN,
     %        GSTRMATCH,NGSTRMATCH,ERR)
	    IF(ERR.NE.' ') GOTO 960
	    MINFUN%X=FC%X
	  ELSE
	    NGSTRMATCH0=NGSTRMATCH
C             save the string temporarily in GSTRMATCH
	    CALL FLCHSET3(GSTR2(EVALLAST)(FC%C(1):FC%C(2)),FC,
     %      GSTRMATCH,NGSTRMATCH0,ERR)
	    IF(ERR.NE.' ') GOTO 960
	    CALL EVAL0(GSTRMATCH(FC%C(1):FC%C(2)),FC,ERR)
	    IF(ERR.NE.' ') GOTO 950
	    IF(FC%L.NE.1) GOTO 950
	    CALL FLCHSET3(GSTR2(EVALLAST)(FC%C(1):FC%C(2)),MINFUN,
     %      GSTRMATCH,NGSTRMATCH,ERR)
	    IF(ERR.NE.' ') GOTO 960
	    MINFUN%X=FC%X
	  ENDIF
	ENDIF

C  Check if variables appear
	DO N=1,NVAR
	  IF(KVAR(1,N).GT.0) THEN
	    VARNAM(N)=NAMPAR(KVAR(1,N))
	  ELSE
	    VARNAM(N)=ARR(-KVAR(1,N))%NAME
	  ENDIF
	ENDDO
	LVAR=0
	DO I=1,BL(IDBLMATCH)%NEXP
	  IF(BL(IDBLMATCH)%MAGNID(I).NE.1) CYCLE
	  J=BL(IDBLMATCH)%MAGID(I)
	  DO K=1,2
	    IF(K.EQ.1) THEN
	      FC=MAG(J)%LENGTH
	    ELSEIF(K.EQ.2) THEN
	      FC=MAG(J)%K1
	    ENDIF
	    IF(FC%L.EQ.2) THEN
	      CALL EVCHKVAR(GSTRMG(FC%C(1):FC%C(2)),MP,NP,NAMP,
     %        MF,NF,NAMF,FC,IRTN)
	      IF(IRTN.NE.0) GOTO 922
	      CALL COUNTVAR(NVAR,VARNAM,LVAR,NP,NAMP,NHIT)
	    ENDIF
	  ENDDO
	ENDDO
	IF(LPER.EQ.0) THEN
	  DO I=1,8
	    IF(TWISSIN(I)%L.EQ.2) THEN
	      CALL EVCHKVAR(GSTRMATCH(TWISSIN(I)%C(1):TWISSIN(I)%C(2)),
     %        MP,NP,NAMP,MF,NF,NAMF,FC,IRTN)
	      IF(IRTN.NE.0) GOTO 924
	      CALL COUNTVAR(NVAR,VARNAM,LVAR,NP,NAMP,NHIT)
	    ENDIF
	  ENDDO
	ENDIF
	DO I=1,NCOND
	  IF(COND(I)%L.EQ.2) THEN
	    CALL EVCHKVAR(GSTRMATCH(COND(I)%C(1):COND(I)%C(2)),
     %      MP,NP,NAMP,MF,NF,NAMF,FC,IRTN)
	    IF(IRTN.NE.0) GOTO 927
	    CALL COUNTVAR(NVAR,VARNAM,LVAR,NP,NAMP,NHIT)
c	    IF(NHIT.EQ.0) THEN
c	      IF(MSGLVL.GE.0) WRITE(MSGFL,400)
c400         FORMAT(' (SUBR.RDMATCH0) The matching condition "',A,'" ',/,
c     %        '    does not contain matching variables. Ignored.')
c	      COND(I)%L=0
c	    ENDIF
C            (turned out hard to check because the optics functions
C             may contain variables.)
	  ENDIF
	ENDDO
	N=0
	N1=0
	NCCONDMAX=0
	DO I=1,NCOND
	  IF(COND(I)%L.EQ.2) THEN
	    N=N+1
	    IF(I.NE.N) COND(N)=COND(I)
	    IF(I.GT.NCOND-NPOS) N1=N1+1
	    NCCONDMAX=MAX(NCCONDMAX,COND(N)%C(2)-COND(N)%C(1)+1)
	  ENDIF
	ENDDO
	NCOND=N
	NPOS=N1
	IF(NCOND.EQ.0) GOTO 935
	DO I=1,NVAR
	  DO J=MCHAR,1,-1
	    IF(VARNAM(I)(J:J).NE.' ') THEN
	      NCVARNAM(I)=J
	      EXIT
	    ENDIF
	  ENDDO
C          Add subscripts for arrays (for printing only)
	  IF(KVAR(1,I).LT.0) THEN
	    IF(ARR(-KVAR(1,I))%RANK.GE.1) THEN
	      CALL ARRN2IND(IND,ARR(-KVAR(1,I))%RANK,
     %        ARR(-KVAR(1,I))%DIM,KVAR(2,I),IRTN)
	      CALL ARRAYSTR(VARNAM(I)(1:NCVARNAM(I)),ARR(-KVAR(1,I))%RANK,
     %        1,1,IND,ERR,NCVARNAM(I))
C                   ERR used as temp work
	      NCVARNAM(I)=MIN(LEN(VARNAM),NCVARNAM(I))
	      VARNAM(I)=ERR
	    ENDIF
	  ENDIF
	ENDDO
	N=0
	NCVARNAMMAX=0
	DO 500 I=1,NVAR
	  IF(LVAR(I).NE.0) THEN
C         Check if the same variable already appears
          IF(N.GE.1) THEN
	      DO J=1,N
	        IF(KVAR(1,J).EQ.KVAR(1,I).AND.KVAR(2,J).EQ.KVAR(2,I)) THEN
	        IF(MSGLVL.GE.0) WRITE(MSGFL,410) VARNAM(I)(1:NCVARNAM(I))
410             FORMAT(' (SUBR.RDMATCH0) The variable "',A,'" ',
     %               'duplicates. Ignored.')
	          GOTO 500
	        ENDIF
	      ENDDO
	    ENDIF
	    N=N+1
	    IF(I.NE.N) THEN
	      KVAR(1:2,N)=KVAR(1:2,I)
	      VARNAM(N)=VARNAM(I)
	      NCVARNAM(N)=NCVARNAM(I)
	    ENDIF
	    NCVARNAMMAX=MAX(NCVARNAMMAX,NCVARNAM(N))
	  ELSE
	    IF(MSGLVL.GE.0) WRITE(MSGFL,420) VARNAM(I)(1:NCVARNAM(I))
420       FORMAT(' (SUBR.RDMATCH0) The variable "',A,'" does not ',
     %      'appear in conditions/magnet parameters. Ignored.')
	  ENDIF
500   CONTINUE
	NVAR=N
	IF(NVAR.LE.0) GOTO 930
	IF(NVAR.LT.NCOND-NPOS) GOTO 932
      CALL MATCHING(MSGLVL,MSGFL,IRTN)
	IF(IRTN.NE.0) GOTO 1000

	IRTN=0
	GOTO 1000

900   IRTN=1000
	IF(MSGLVL.GE.0) WRITE(MSGFL,902)
902   FORMAT(' (SUBR.RDMATCH0) Allocation error.')
	GOTO 1000
910   IRTN=1010
	IF(MSGLVL.GE.0) WRITE(MSGFL,911)
911   FORMAT(' (SUBR.RDMATCH0) No variable specified.')
	GOTO 1000
912   IRTN=1012
      WRITE(MSGFL,913) TXT(3)(1:NCTXT(3)),ERR
913   FORMAT(' (SUBR.RDMATCH0) Invalid expression "',A,'"',
     %  ' for MINFUN',/,3X,A)
      GOTO 1000
914   IRTN=1014
      WRITE(MSGFL,915) TXT(I)(ITEM(1,K):ITEM(2,K)),KWD(I),ERR
915   FORMAT(' (SUBR.RDMATCH0) Invalid expression "',A,'"',
     %  ' for ',A,/,3X,A)
      GOTO 1000

916   IRTN=1016
      WRITE(MSGFL,917) NAM(1:NCNAM)
917   FORMAT(' (SUBR.RDMATCH0) Invalid array subscript for "',A,'"')
      GOTO 1000
918   IRTN=1018
      WRITE(MSGFL,919) NAM(1:NCNAM)
919   FORMAT(' (SUBR.RDMATCH0) Character variable "',A,'" cannot be ',
     %  'a matching variable.')
      GOTO 1000
920   IRTN=1020
      WRITE(MSGFL,921) NAM(1:NCNAM)
921   FORMAT(' (SUBR.RDMATCH0) Undefined variable name "',A,'"')
      GOTO 1000
922   IRTN=1022
	IF(MSGLVL.GE.0) WRITE(MSGFL,923)
     %   GSTRMG(FC%C(1):FC%C(2)),MAG(J)%NAME
923   FORMAT(' (SUBR.RDMATCH0) Invalid expression ',A,' for magnet ',A)
	GOTO 1000
924   IRTN=1024
	IF(MSGLVL.GE.0) WRITE(MSGFL,926)
     %   GSTRMG(TWISSIN1(I)%C(1):TWISSIN1(I)%C(2))
926   FORMAT(' (SUBR.RDMATCH0) Invalid Twiss param expression ',A)
	GOTO 1000
927   IRTN=1026
	IF(MSGLVL.GE.0) WRITE(MSGFL,928)
     %   GSTRMG(COND(I)%C(1):COND(I)%C(2))
928   FORMAT(' (SUBR.RDMATCH0) Invalid matching condition ',A)
	GOTO 1000
930   IRTN=1030
	IF(MSGLVL.GE.0) WRITE(MSGFL,931)
931   FORMAT(' (SUBR.RDMATCH0) No actually used variables.')
	GOTO 1000
932   IRTN=1032
	IF(MSGLVL.GE.0) WRITE(MSGFL,933) NVAR,NCOND
933   FORMAT(' (SUBR.RDMATCH0) No. of variables=',I2,' less than ',
     %   'no. of equality constraints=',I2)
	GOTO 1000
935   IRTN=1035
      WRITE(MSGFL,936)
936   FORMAT(' (SUBR.RDMATCH0) No conditions specified.')
      GOTO 1000
950   IRTN=1050
      IF(MSGLVL.GE.0) WRITE(MSGFL,951) GSTRRD(FC%C(1):FC%C(2))
951   FORMAT(' (SUBR.RDMATCH0) Invalid expression ',A)
      GOTO 1000
960   IRTN=1060
      WRITE(MSGFL,961)
961   FORMAT(' (SUBR.RDMATCH0) Character buffer GSTRRD full.')
      GOTO 1000
965   IRTN=1065
      WRITE(MSGFL,966)
966   FORMAT(' (SUBR.RDMATCH0) Character buffer GSTRMATCH full.')
      GOTO 1000
990   IRTN=1090
	IF(MSGLVL.GE.0) WRITE(MSGFL,992) ERR
992   FORMAT(' (SUBR.RDMATCH0) ',A)
	GOTO 1000

1000  IDBLMATCH=0
      DEALLOCATE(KVAR,LVAR,XVAR,XVAR0,VARNAM,NCVARNAM,COND,STAT=ISTAT)
	RETURN
	END

	
	SUBROUTINE COUNTVAR(NVAR,VARNAM,LVAR,NP,NAMP,NHIT)
	IMPLICIT NONE
	INTEGER NVAR,LVAR(NVAR),NP,NHIT
	CHARACTER(*) VARNAM(NVAR),NAMP(NP)
	INTEGER I,N
	NHIT=0
	IF(NP.NE.0) THEN
	  DO I=1,NP
	    DO N=1,NVAR
	      IF(NAMP(I).EQ.VARNAM(N)) THEN
	        NHIT=NHIT+1
	        LVAR(N)=LVAR(N)+1
C                Loop still continues because VARNAM(N) does not contain
C                the array element index. Other VARNAM may be the same.
	      ENDIF
	    ENDDO
	  ENDDO
	ENDIF
	RETURN
	END

