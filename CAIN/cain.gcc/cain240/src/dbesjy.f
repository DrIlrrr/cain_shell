C******************* DBESJY **************************************
      FUNCTION DBESJY(NORDER,JY,X)
C  REAL*8 BESSEL FUNCTION I0,I1,K0,K1.   CERNLIB C312
C  DBESJ0(X)=DBESJY(0,1,X)=J0(X),
C  DBESJ1(X)=DBESJY(1,1,X)=J1(X),
C  DBESY0(X)=DBESJY(0,2,X)=Y0(X),
C  DBESY1(X)=DBESJY(1,2,X)=Y1(X),
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL L
	INCLUDE 'include/ctrlcm.h'
C
      IENT=0
      IF(NORDER.LT.0.OR.NORDER.GE.2) GOTO 920
      IF(JY.LE.0.OR.JY.GE.3) GOTO 940
      IF(JY.EQ.1) THEN
          IF(NORDER.EQ.0) GOTO 100
          GOTO 200
        ELSE
          IF(NORDER.EQ.0) GOTO 300
          GOTO 400
        ENDIF
C
      ENTRY DBESJ0(X)
      IENT=1
  100 L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 4
    8 F=0.0625D0*X**2-2.0
      A =           - 0.00000 00000 000008D0
      B = F * A     + 0.00000 00000 000413D0
      A = F * B - A - 0.00000 00000 019438D0
      B = F * A - B + 0.00000 00000 784870D0
      A = F * B - A - 0.00000 00026 792535D0
      B = F * A - B + 0.00000 00760 816359D0
      A = F * B - A - 0.00000 17619 469078D0
      B = F * A - B + 0.00003 24603 288210D0
      A = F * B - A - 0.00046 06261 662063D0
      B = F * A - B + 0.00481 91800 694676D0
      A = F * B - A - 0.03489 37694 114089D0
      B = F * A - B + 0.15806 71023 320973D0
      A = F * B - A - 0.37009 49938 726498D0
      B = F * A - B + 0.26517 86132 033368D0
      A = F * B - A - 0.00872 34423 528522D0
      A = F * A - B + 0.31545 59429 497802D0
      DBESJY=0.5*(A-B)
      IF(L) GOTO 500

      A =           + 0.00000 00000 000016D0
      B = F * A     - 0.00000 00000 000875D0
      A = F * B - A + 0.00000 00000 040263D0
      B = F * A - B - 0.00000 00001 583755D0
      A = F * B - A + 0.00000 00052 487948D0
      B = F * A - B - 0.00000 01440 723327D0
      A = F * B - A + 0.00000 32065 325377D0
      B = F * A - B - 0.00005 63207 914106D0
      A = F * B - A + 0.00075 31135 932578D0
      B = F * A - B - 0.00728 79624 795521D0
      A = F * B - A + 0.04719 66895 957634D0
      B = F * A - B - 0.17730 20127 811436D0
      A = F * B - A + 0.26156 73462 550466D0
      B = F * A - B + 0.17903 43140 771827D0
      A = F * B - A - 0.27447 43055 297453D0
      A = F * A - B - 0.06629 22264 065699D0
      DBESJY=0.636619772367581D0*LOG(X)*DBESJY+0.5*(A-B)
      GOTO 500

    4 F=256.D0/X**2-2.0
      B =           + 0.00000 00000 000007D0
      A = F * B     - 0.00000 00000 000051D0
      B = F * A - B + 0.00000 00000 000433D0
      A = F * B - A - 0.00000 00000 004305D0
      B = F * A - B + 0.00000 00000 051683D0
      A = F * B - A - 0.00000 00000 786409D0
      B = F * A - B + 0.00000 00016 306465D0
      A = F * B - A - 0.00000 00517 059454D0
      B = F * A - B + 0.00000 30751 847875D0
      A = F * B - A - 0.00053 65220 468132D0
      A = F * A - B + 1.99892 06986 950373D0
      P=A-B
      B =           - 0.00000 00000 000006D0
      A = F * B     + 0.00000 00000 000043D0
      B = F * A - B - 0.00000 00000 000334D0
      A = F * B - A + 0.00000 00000 003006D0
      B = F * A - B - 0.00000 00000 032067D0
      A = F * B - A + 0.00000 00000 422012D0
      B = F * A - B - 0.00000 00007 271916D0
      A = F * B - A + 0.00000 00179 724572D0
      B = F * A - B - 0.00000 07414 498411D0
      A = F * B - A + 0.00006 83851 994261D0
      A = F * A - B - 0.03111 17092 106740D0
      Q=8.0*(A-B)/V
      F=V-0.785398163397448D0
      A=COS(F)
      B=SIN(F)
      F=0.398942280401432D0/SQRT(V)
      IF(L) GO TO 6
      DBESJY=F*(Q*A+P*B)
      GOTO 500
    6 DBESJY=F*(P*A-Q*B)
      GOTO 500

      ENTRY DBESJ1(X)
      IENT=2

  200 L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 5
    3 F=0.0625D0*X**2-2.0
      B =           + 0.00000 00000 000114D0
      A = F * B     - 0.00000 00000 005777D0
      B = F * A - B + 0.00000 00000 252812D0
      A = F * B - A - 0.00000 00009 424213D0
      B = F * A - B + 0.00000 00294 970701D0
      A = F * B - A - 0.00000 07617 587805D0
      B = F * A - B + 0.00001 58870 192399D0
      A = F * B - A - 0.00026 04443 893486D0
      B = F * A - B + 0.00324 02701 826839D0
      A = F * B - A - 0.02917 55248 061542D0
      B = F * A - B + 0.17770 91172 397283D0
      A = F * B - A - 0.66144 39341 345433D0
      B = F * A - B + 1.28799 40988 576776D0
      A = F * B - A - 1.19180 11605 412169D0
      A = F * A - B + 1.29671 75412 105298D0
      DBESJY=0.0625D0*(A-B)*X
      IF(L) GOTO 500

      B =           - 0.00000 00000 000244D0
      A = F * B     + 0.00000 00000 012114D0
      B = F * A - B - 0.00000 00000 517212D0
      A = F * B - A + 0.00000 00018 754703D0
      B = F * A - B - 0.00000 00568 844004D0
      A = F * B - A + 0.00000 14166 243645D0
      B = F * A - B - 0.00002 83046 401495D0
      A = F * B - A + 0.00044 04786 298671D0
      B = F * A - B - 0.00513 16411 610611D0
      A = F * B - A + 0.04231 91803 533369D0
      B = F * A - B - 0.22662 49915 567549D0
      A = F * B - A + 0.67561 57807 721877D0
      B = F * A - B - 0.76729 63628 866459D0
      A = F * B - A - 0.12869 73843 813500D0
      A = F * A - B + 0.04060 82117 718685D0
      DBESJY=0.636619772367581D0*LOG(X)*DBESJY-0.636619772367581D0/X
     1     +0.0625D0*(A-B)*X
      GOTO 500

    5 F=256.D0/X**2-2.0
      B =           - 0.00000 00000 000007D0
      A = F * B     + 0.00000 00000 000055D0
      B = F * A - B - 0.00000 00000 000468D0
      A = F * B - A + 0.00000 00000 004699D0
      B = F * A - B - 0.00000 00000 057049D0
      A = F * B - A + 0.00000 00000 881690D0
      B = F * A - B - 0.00000 00018 718907D0
      A = F * B - A + 0.00000 00617 763396D0
      B = F * A - B - 0.00000 39872 843005D0
      A = F * B - A + 0.00089 89898 330859D0
      A = F * A - B + 2.00180 60817 200274D0
      P=A-B
      B =           + 0.00000 00000 000007D0
      A = F * B     - 0.00000 00000 000046D0
      B = F * A - B + 0.00000 00000 000360D0
      A = F * B - A - 0.00000 00000 003264D0
      B = F * A - B + 0.00000 00000 035152D0
      A = F * B - A - 0.00000 00000 468636D0
      B = F * A - B + 0.00000 00008 229193D0
      A = F * B - A - 0.00000 00209 597814D0
      B = F * A - B + 0.00000 09138 615258D0
      A = F * B - A - 0.00009 62772 354916D0
      A = F * A - B + 0.09355 55741 390707D0
      Q=8.0*(A-B)/V
      F=V-2.356194490192345D0
      A=COS(F)
      B=SIN(F)
      F=0.398942280401432D0/SQRT(V)
      IF(L) GO TO 7
      DBESJY=F*(Q*A+P*B)
      GOTO 500
    7 DBESJY=F*(P*A-Q*B)
      IF(X .LT. 0.0) DBESJY=-DBESJY
      GOTO 500

      ENTRY DBESY0(X)
      IENT=3

  300 IF(X .LE. 0.0) GO TO 900
      L=.FALSE.
      V=X
      IF(V .GE. 8.0) GO TO 4
      GO TO 8

      ENTRY DBESY1(X)
      IENT=4

  400 IF(X .LE. 0.0) GO TO 900
      L=.FALSE.
      V=X
      IF(V .GE. 8.0) GO TO 5
      GO TO 3
C-- NORMAL RETURN
 500  GOTO (510,520,530,540,550), IENT+1
 510  RETURN
 520  DBESJ0=DBESJY
      RETURN
 530  DBESJ1=DBESJY
      RETURN
 540  DBESY0=DBESJY
      RETURN
 550  DBESY1=DBESJY
      RETURN
C-- ERROR RETURN
  900 DBESJY=0.
      WRITE(MSGFL,910) X
  910 FORMAT(' (FUNC.DBESJY) NON-POSITIVE ARGUMENT X =',1PD14.7)
      RETURN
  920 WRITE(MSGFL,930) NORDER
  930 FORMAT(' (FUNC.DBESJY) NORDER=',I5,' MUST BE 0 OR 1.')
      RETURN
  940 WRITE(MSGFL,950) JY
  950 FORMAT(' (FUNC.DBESJY) JY=',I5,' MUST BE 1 OR 2.')
      RETURN
      END
