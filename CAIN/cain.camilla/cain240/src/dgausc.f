      SUBROUTINE DGAUSC(X,N,GCUT)
      IMPLICIT NONE
      INTEGER N
      REAL*8 X(N),GCUT
      INTEGER I
C
      CALL DGAUSN(X,N)
      IF(GCUT.LE.0) RETURN
      DO 200 I=1,N
 160    IF(ABS(X(I)).GT.GCUT) THEN
          CALL DGAUSN(X(I),1)
          GOTO 160
        ENDIF
 200  CONTINUE
      RETURN
      END
