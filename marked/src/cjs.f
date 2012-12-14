       SUBROUTINE CJS(CH,PHI,P,FRST,LST,FRQ,LOC,PHIF,PF,TINT,N,M,K,L,
     $                     LNL,P0)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER LOC(N), IFRST, ILST, KI
       DOUBLE PRECISION CH(N,M),PHI(N,M-1),P(N,M-1),FRST(N)
       DOUBLE PRECISION LNL, LST(N),PCH(N),P0(N)
       DOUBLE PRECISION PHICUMPROD(M), CUMP(M)
       DOUBLE PRECISION FRQ(N), TINT(N,M-1), PHIF(K,3), PF(L,3)
       DOUBLE PRECISION X,Y
       LNL=0D0
       IFLAG=0
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
	     IF(INT(FRST(I)).LT.INT(M)) THEN
	        PHICUMPROD(INT(FRST(I)))=1D0
	        DO 10 J=INT(FRST(I))+1,INT(M)
		    	IF(TINT(I,J-1).NE.1D0)PHI(I,J-1)=PHI(I,J-1)**TINT(I,J-1)
	            PHICUMPROD(J)=PHICUMPROD(J-1)*PHI(I,J-1)
  10        CONTINUE
         ENDIF
         IF(INT(FRST(I)).LT.INT(M))THEN
           CUMP(INT(FRST(I)))=1D0
           DO 20 J=INT(FRST(I))+1,INT(M)
 	         IF(CH(I,J)<.5) THEN
		        CUMP(J)=CUMP(J-1)*(1-P(I,J-1))
	         ELSE
		        CUMP(J)=CUMP(J-1)*P(I,J-1)
	         ENDIF
  20       CONTINUE
         ENDIF
         IF(FRST(I).EQ.M) THEN
            P0(I)=1
            PCH(I)=1D0
         ELSE
            PCH(I)=0
            KI=INT((1-LOC(I))*M+LOC(I)*LST(I))
	        DO 30 J=INT(LST(I)),KI
	          IF(LOC(I).EQ.1.OR.J.EQ.M)THEN
	             PCH(I)=PCH(I)+CUMP(J)*PHICUMPROD(J)
	          ELSE
	             PCH(I)=PCH(I)+CUMP(J)*PHICUMPROD(J)*(1-PHI(I,J))
	          ENDIF
	          IF(PCH(I)<1E-15.AND.FRQ(I).GT.0) PCH(I)=1D-307
  30        CONTINUE
            P0(I)=PCH(I)/(CUMP( INT(LST(I)))*PHICUMPROD(INT(LST(I))))
         ENDIF
         LNL=LNL-FRQ(I)*LOG(PCH(I))
  100  CONTINUE
       RETURN
       END
  