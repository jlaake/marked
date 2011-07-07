      SUBROUTINE CJS(CH,PHI,P,FRST,LST,FRQ,LOC,PHIF,PF,TINT,N,M,K,L,
     $                     LNL,P0)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER LOC(N), IFRST, ILST
       DOUBLE PRECISION CH(N,M),PHI(N,M-1),P(N,M-1),FRST(N)
       DOUBLE PRECISION LNL, LST(N), CUMPHI(N),P0(N)
       DOUBLE PRECISION PHICUMPROD(N,M), CUMP(N), PLAST(N,M)
       DOUBLE PRECISION FRQ(N), TINT(N,M-1), PHIF(K,3), PF(L,3)
       DOUBLE PRECISION X,Y
       LNL=0D0
       DO 2 I=1,N
       DO 2 J=INT(FRST(I))+1,M
         P(I,J-1)=1/(1+EXP(-P(I,J-1)))
	   PHI(I,J-1)=1/(1+EXP(-PHI(I,J-1)))
  2    CONTINUE
       IF(K.GT.1 .OR. PHIF(1,1).GT.0)THEN
	   DO 3 I=1,K
           PHI(INT(PHIF(I,1)),INT(PHIF(I,2)))=PHIF(I,3)
  3      CONTINUE
       ENDIF
       IF(L.GT.1 .OR. PF(1,1).GT.0)THEN
	   DO 4 I=1,L
           P(INT(PF(I,1)),INT(PF(I,2))-1)=PF(I,3)
  4      CONTINUE
       ENDIF
 	 DO 100 I=1,N
	   DO 5 J=1,INT(FRST(I))
           PHICUMPROD(I,J)=1D0
  5      CONTINUE
	   IF(INT(FRST(I)).LT.M) THEN
	      DO 10 J=INT(FRST(I))+1,M
			IF(TINT(I,J-1).NE.1D0)PHI(I,J-1)=PHI(I,J-1)**TINT(I,J-1)
	        PHICUMPROD(I,J)=PHICUMPROD(I,J-1)*PHI(I,J-1)
  10        CONTINUE
         ENDIF
	   CUMPHI(I)=PHICUMPROD(I,INT(LST(I)))
	   CUMP(I)=1D0
            IF(INT(FRST(I)).LT.INT(LST(I)))THEN
            DO 20 J=INT(FRST(I))+1,INT(LST(I))
	         IF(CH(I,J)<.5) THEN
		        CUMP(I)=CUMP(I)*(1-P(I,J-1))
	         ELSE
		        CUMP(I)=CUMP(I)*P(I,J-1)
	         ENDIF
  20        CONTINUE
         ENDIF
         PLAST(I,INT(LST(I)))=1
	   IF(LST(I).GT.1) THEN
	      DO 30 J=1,INT(LST(I))-1
	         PLAST(I,J)=0D0
  30        CONTINUE
         ENDIF
	   IF(INT(LST(I)).LE.M-1) THEN
	      DO 40 J=INT(LST(I))+1,M
	        PLAST(I,J)=PLAST(I,J-1)*(1-P(I,J-1))
  40         CONTINUE
         ENDIF
         IF(CUMP(I).LE.0) THEN
	      CUMP(I)=1D-8
	   ENDIF
         IF(CUMPHI(I).LE.0) THEN
	      CUMPHI(I)=1D-8
	   ENDIF
         P0(I)=0D0
	   IF(LOC(I).EQ.0.AND.CUMPHI(I).GT.0D0) THEN
	    DO 50 J=INT(LST(I)),M
	    IF(J.EQ.M) THEN
	    P0(I)=P0(I)+PHICUMPROD(I,J)/CUMPHI(I)*PLAST(I,J)
	    ELSE
	    P0(I)=P0(I)+PHICUMPROD(I,J)/CUMPHI(I)*PLAST(I,J)*(1-PHI(I,J))
          ENDIF
  50       CONTINUE
          IF(P0(I).LE. 0D0) THEN
		   P0(I)=1D-8
	    ENDIF
          LNL=LNL-2*FRQ(I)*(LOG(P0(I))+LOG(CUMP(I))+LOG(CUMPHI(I)))
         ELSE
          LNL=LNL-2*FRQ(I)*(LOG(CUMP(I))+LOG(CUMPHI(I)))
	   ENDIF
  100  CONTINUE
       RETURN
       END
 