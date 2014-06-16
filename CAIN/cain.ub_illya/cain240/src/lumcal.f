      SUBROUTINE LUMCAL(IADRS,NP,KIND,EP,WGT,SPIN,IRTN)
	USE LUMCOM
      IMPLICIT NONE
      INTEGER NP,IADRS(2,NP),KIND(NP),IRTN
      REAL*8 EP(0:3,NP),WGT(NP),SPIN(3,NP)
      INCLUDE 'include/lumcom2.h'
      INTEGER M,K,L,K0,KIN,KIN1,KIN2,IL,NMDIM
      REAL*8 WBIN,EBIN(2),EMIN(2)
C
      CALL CPUTIM('LUMCAL',1) 
      DO 280 L=1,2
        DO 260 KIN=1,3
          IF(NPKINL(KIN,L).NE.0) THEN
            DO 240 M=MM-1,0,-1
              DO 220 K=0,NDD2(M)-1
                K0=IPDD(M+1)+4*K
                DIST(IPDD(M)+K,KIN,L)=DIST(K0,KIN,L)
     %              +DIST(K0+1,KIN,L)+DIST(K0+2,KIN,L)+DIST(K0+3,KIN,L)
                VDIST(IPDD(M)+K,KIN,L)=VDIST(K0,KIN,L)
     %           +VDIST(K0+1,KIN,L)+VDIST(K0+2,KIN,L)+VDIST(K0+3,KIN,L)
 220          CONTINUE
 240        CONTINUE
          ENDIF
 260    CONTINUE
 280  CONTINUE
C
      NMDIM=(4*NDD2(MM)-1)/3
      DO 400 IL=1,NLUM
        KIN1=KLUM(1,IL)
        KIN2=KLUM(2,IL)
        IF(NPKINL(KIN1,1).EQ.0.OR.NPKINL(KIN2,2).EQ.0) GOTO 400
        IF(LBINW(IL).EQ.1) THEN
          WBIN=(WMMLUM(2,IL)-WMMLUM(1,IL))/NBNWLM(IL)
        ELSE
          WBIN=0
        ENDIF
        DO 320 K=1,2
          IF(LBINE(1,IL).NE.0) THEN
            EMIN(K)=EMMLUM(1,K,IL)
            IF(LBINE(K,IL).EQ.1) THEN
              EBIN(K)=(EMMLUM(2,K,IL)-EMMLUM(1,K,IL))/NBNELM(K,IL)
            ENDIF
          ENDIF
 320    CONTINUE
        CALL LUMCAL0(NMDIM,DIST(0,KIN1,1),DIST(0,KIN2,2),
     %     VDIST(0,KIN1,1),VDIST(0,KIN2,2),LUMFAC,LHEL(0,IL),
     %     LUM(0,IL),VLUM(IL),
     %     DLUM(IPWLUM(0,IL)),DLUM(IPWLUM(1,IL)),DLUM(IPWLUM(2,IL)),
     %     LBINW(IL),NBNWLM(IL),WMMLUM(1,IL),WBIN,DLUM(IPBINW(IL)),
     %     DLUM(IPELUM(0,IL)),DLUM(IPELUM(1,IL)),DLUM(IPELUM(2,IL)),
     %     LBINE(1,IL),NBNELM(1,IL),NBNELM(2,IL),EMIN,EBIN,
     %     DLUM(IPBINE(1,IL)),DLUM(IPBINE(2,IL)),
     %     IPPFLG(IL),
     %     MM,IPDD,NDD,NDD2,
     %     NDD2(MM),IPBN(0,KIN1,1),IPBN(0,KIN2,2),
     %     IADRS,NP,KIND,EP,WGT,SPIN,IRTN)
        IF(IRTN.NE.0) GOTO 1000
 400  CONTINUE
C
      IRTN=0
 1000 CALL CPUTIM('LUMCAL',2) 
      RETURN
      END
