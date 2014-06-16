C******************** DBESIK *******************************************
      FUNCTION DBESIK(N,IK,X)
C  REAL*8 BESSEL FUNCTION I0,I1,K0,K1.   CERNLIB
C  BESI0(X)=DBESIK(0,1,X)=I0(X), EBESI0(X)=DBESIK(0,-1,X)=EXP(-X)*I0(X)
C  BESI1(X)=DBESIK(1,1,X)=I1(X), EBESI1(X)=DBESIK(1,-1,X)=EXP(-X)*I1(X)
C  BESK0(X)=DBESIK(0,2,X)=K0(X), EBESK0(X)=DBESIK(0,-2,X)=EXP(+X)*K0(X)
C  BESK1(X)=DBESIK(1,2,X)=K1(X), EBESK1(X)=DBESIK(1,-2,X)=EXP(+X)*K1(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL L,E
      IENT=0
      IF(N.LT.0.OR.N.GT.1) GOTO 900
      IK1=ABS(IK)
      IF(IK1.LE.0.OR.IK1.GE.3) GOTO 900
      IF(IK.LT.0) IK1=IK1+2
      GOTO (210,410,200,400,310,510,300,500), 4*N+IK1
      ENTRY EBESI0(X)
      IENT=1
  200 E=.TRUE.
      GO TO 220
      ENTRY BESI0(X)
      IENT=2
  210 E=.FALSE.
  220 L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 260
C I0(K0) SMALL X
  240 F=0.0625D0*X**2-2D0
      A =               0.00000 00000 00002D0
      B = F * A     +   0.00000 00000 00120D0
      A = F * B - A +   0.00000 00000 06097D0
      B = F * A - B +   0.00000 00002 68828D0
      A = F * B - A +   0.00000 00101 69727D0
      B = F * A - B +   0.00000 03260 91051D0
      A = F * B - A +   0.00000 87383 15497D0
      B = F * A - B +   0.00019 24693 59688D0
      A = F * B - A +   0.00341 63317 66012D0
      B = F * A - B +   0.04771 87487 98174D0
      A = F * B - A +   0.50949 33654 39983D0
      B = F * A - B +   4.01167 37601 79349D0
      A = F * B - A +  22.27481 92424 62231D0
      B = F * A - B +  82.48903 27440 24100D0
      A = F * B - A + 190.49432 01727 42844D0
      A = F * A - B + 255.46687 96243 62167D0
      DBESIK=0.5D0*(A-B)
      IF(L .AND. E) DBESIK=EXP(-V)*DBESIK
      IF(L) GOTO 600
C K0 SMALL X
      A =           +   0.00000 00000 00003D0
      B = F * A     +   0.00000 00000 00159D0
      A = F * B - A +   0.00000 00000 07658D0
      B = F * A - B +   0.00000 00003 18588D0
      A = F * B - A +   0.00000 00112 81211D0
      B = F * A - B +   0.00000 03351 95256D0
      A = F * B - A +   0.00000 82160 25940D0
      B = F * A - B +   0.00016 27083 79043D0
      A = F * B - A +   0.00253 63081 88086D0
      B = F * A - B +   0.03008 07224 20512D0
      A = F * B - A +   0.25908 44324 34900D0
      B = F * A - B +   1.51153 56760 29228D0
      A = F * B - A +   5.28363 28668 73920D0
      B = F * A - B +   8.00536 88687 00334D0
      A = F * B - A -   4.56343 35864 48395D0
      A = F * A - B -  21.05766 01774 02440D0
      DBESIK=-LOG(0.125D0*X)*DBESIK+0.5D0*(A-B)
      IF(E) DBESIK=EXP(X)*DBESIK
      GOTO 600
C I0 LARGE X
  260 F=32D0/V-2D0
      B =           - 0.00000 00000 00001D0
      A = F * B     - 0.00000 00000 00001D0
      B = F * A - B + 0.00000 00000 00004D0
      A = F * B - A + 0.00000 00000 00010D0
      B = F * A - B - 0.00000 00000 00024D0
      A = F * B - A - 0.00000 00000 00104D0
      B = F * A - B + 0.00000 00000 00039D0
      A = F * B - A + 0.00000 00000 00966D0
      B = F * A - B + 0.00000 00000 01800D0
      A = F * B - A - 0.00000 00000 04497D0
      B = F * A - B - 0.00000 00000 33127D0
      A = F * B - A - 0.00000 00000 78957D0
      B = F * A - B + 0.00000 00000 29802D0
      A = F * B - A + 0.00000 00012 38425D0
      B = F * A - B + 0.00000 00085 13091D0
      A = F * B - A + 0.00000 00568 16966D0
      B = F * A - B + 0.00000 05135 87727D0
      A = F * B - A + 0.00000 72475 91100D0
      B = F * A - B + 0.00017 27006 30778D0
      A = F * B - A + 0.00844 51226 24921D0
      A = F * A - B + 2.01655 84109 17480D0
      DBESIK=0.199471140200717D0*(A-B)/SQRT(V)
      IF(E) GOTO 600
      DBESIK=EXP(V)*DBESIK
      GOTO 600
      ENTRY EBESI1(X)
      IENT=3
  300 E=.TRUE.
      GO TO 320
      ENTRY BESI1(X)
      IENT=4
  310 E=.FALSE.
  320 L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 360
C I1(K1) SMALL X
  340 F=0.0625D0*X**2-2D0
      A =           +   0.00000 00000 00001D0
      B = F * A     +   0.00000 00000 00031D0
      A = F * B - A +   0.00000 00000 01679D0
      B = F * A - B +   0.00000 00000 79291D0
      A = F * B - A +   0.00000 00032 27617D0
      B = F * A - B +   0.00000 01119 46285D0
      A = F * B - A +   0.00000 32641 38122D0
      B = F * A - B +   0.00007 87567 85754D0
      A = F * B - A +   0.00154 30190 15627D0
      B = F * A - B +   0.02399 30791 47841D0
      A = F * B - A +   0.28785 55118 04672D0
      B = F * A - B +   2.57145 99063 47755D0
      A = F * B - A +  16.33455 05525 22066D0
      B = F * A - B +  69.39591 76337 34448D0
      A = F * B - A + 181.31261 60405 70265D0
      A = F * A - B + 259.89023 78064 77292D0
      DBESIK=0.0625D0*(A-B)*X
      IF(L .AND. E) DBESIK=EXP(-V)*DBESIK
      IF(L) GOTO 600
C K1 SMALL X
      A =           +   0.00000 00000 00001D0
      B = F * A     +   0.00000 00000 00042D0
      A = F * B - A +   0.00000 00000 02163D0
      B = F * A - B +   0.00000 00000 96660D0
      A = F * B - A +   0.00000 00036 96783D0
      B = F * A - B +   0.00000 01193 67971D0
      A = F * B - A +   0.00000 32025 10692D0
      B = F * A - B +   0.00007 00106 27855D0
      A = F * B - A +   0.00121 70569 94516D0
      B = F * A - B +   0.01630 00492 89816D0
      A = F * B - A +   0.16107 43016 56148D0
      B = F * A - B +   1.10146 19930 04852D0
      A = F * B - A +   4.66638 70268 62842D0
      B = F * A - B +   9.36161 78313 95389D0
      A = F * B - A -   1.83923 92242 86199D0
      A = F * A - B -  26.68809 54808 62668D0
      DBESIK=LOG(0.125D0*X)*DBESIK+1D0/X-0.0625D0*(A-B)*X
      IF(E) DBESIK=EXP(X)*DBESIK
      GOTO 600
C I1 LARGE X
  360 F=32D0/V-2D0
      B =           + 0.00000 00000 00001D0
      A = F * B     + 0.00000 00000 00001D0
      B = F * A - B - 0.00000 00000 00005D0
      A = F * B - A - 0.00000 00000 00010D0
      B = F * A - B + 0.00000 00000 00026D0
      A = F * B - A + 0.00000 00000 00107D0
      B = F * A - B - 0.00000 00000 00053D0
      A = F * B - A - 0.00000 00000 01024D0
      B = F * A - B - 0.00000 00000 01804D0
      A = F * B - A + 0.00000 00000 05103D0
      B = F * A - B + 0.00000 00000 35408D0
      A = F * B - A + 0.00000 00000 81531D0
      B = F * A - B - 0.00000 00000 47563D0
      A = F * B - A - 0.00000 00014 01141D0
      B = F * A - B - 0.00000 00096 13873D0
      A = F * B - A - 0.00000 00659 61142D0
      B = F * A - B - 0.00000 06297 24239D0
      A = F * B - A - 0.00000 97321 46728D0
      B = F * A - B - 0.00027 72053 60764D0
      A = F * B - A - 0.02446 74429 63276D0
      A = F * A - B + 1.95160 12046 52572D0
      DBESIK=0.199471140200717D0*(A-B)/SQRT(V)
      IF(X .LT. 0.0) DBESIK=-DBESIK
      IF(E) GOTO 600
      DBESIK=EXP(V)*DBESIK
      GOTO 600
      ENTRY EBESK0(X)
      IENT=5
  400 E=.TRUE.
      GO TO 420
      ENTRY BESK0(X)
      IENT=6
  410 E=.FALSE.
  420 IF(X .LE. 0.0) GO TO 920
      L=.FALSE.
      V=X
      IF(X .LT. 5.0) GO TO 240
C K0 LARGE X
      F=20D0/X-2D0
      A =           - 0.00000 00000 00002D0
      B = F * A     + 0.00000 00000 00011D0
      A = F * B - A - 0.00000 00000 00079D0
      B = F * A - B + 0.00000 00000 00581D0
      A = F * B - A - 0.00000 00000 04580D0
      B = F * A - B + 0.00000 00000 39044D0
      A = F * B - A - 0.00000 00003 64547D0
      B = F * A - B + 0.00000 00037 92996D0
      A = F * B - A - 0.00000 00450 47338D0
      B = F * A - B + 0.00000 06325 75109D0
      A = F * B - A - 0.00001 11066 85197D0
      B = F * A - B + 0.00026 95326 12763D0
      A = F * B - A - 0.01131 05046 46928D0
      A = F * A - B + 1.97681 63484 61652D0
      DBESIK=0.626657068657750D0*(A-B)/SQRT(X)
      IF(E) GOTO 600
      Z=DBESIK
      DBESIK=0.
      IF(X.LT.180.) DBESIK=EXP(-X)*Z
      GOTO 600
      ENTRY EBESK1(X)
      IENT=7
  500 E=.TRUE.
      GO TO 520
      ENTRY BESK1(X)
      IENT=8
  510 E=.FALSE.
  520 IF(X .LE. 0.0) GO TO 920
      L=.FALSE.
      V=X
      IF(X .LT. 5.0) GO TO 340
C K1 LARGE X
      F=20D0/X-2D0
      A =           + 0.00000 00000 00002D0
      B = F * A     - 0.00000 00000 00013D0
      A = F * B - A + 0.00000 00000 00089D0
      B = F * A - B - 0.00000 00000 00663D0
      A = F * B - A + 0.00000 00000 05288D0
      B = F * A - B - 0.00000 00000 45757D0
      A = F * B - A + 0.00000 00004 35417D0
      B = F * A - B - 0.00000 00046 45555D0
      A = F * B - A + 0.00000 00571 32218D0
      B = F * A - B - 0.00000 08451 72048D0
      A = F * B - A + 0.00001 61850 63810D0
      B = F * A - B - 0.00046 84750 28167D0
      A = F * B - A + 0.03546 52912 43331D0
      A = F * A - B + 2.07190 17175 44716D0
      DBESIK=0.626657068657750D0*(A-B)/SQRT(X)
      IF(E) GOTO 600
      Z=DBESIK
      DBESIK=0.
      IF(X.LT.180.) DBESIK=EXP(-X)*Z
      GOTO 600
C-- NORMAL RETURN
 600  GOTO (700,610,620,630,640,650,660,670,680), IENT+1
 610  EBESI0=DBESIK
      RETURN
 620  BESI0=DBESIK
      RETURN
 630  EBESI1=DBESIK
      RETURN
 640  BESI1=DBESIK
      RETURN
 650  EBESK0=DBESIK
      RETURN
 660  BESK0=DBESIK
      RETURN
 670  EBESK1=DBESIK
      RETURN
 680  BESK1=DBESIK
      RETURN
 700  RETURN
C-- ERROR RETURN
  900 WRITE(6,910) N,IK,X
  910 FORMAT(' (FUNC.DBESIK) INVALID ARGUMENT. N=',I5,' IK=',I5,
     % ' X=',1PD13.6)
      RETURN
  920 DBESIK=0.
      WRITE(6,930) X
  930 FORMAT(' (FUNC.DBESIK) NON-POSITIVE ARGUMENT X =',1PD13.6)
      RETURN
      END
