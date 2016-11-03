       SUBROUTINE MVMSP(P,DELTA,N,M,F,T,NP,NOBS,PCOUNTS,INDICES,
     $                      UNK,UNKINIT,PMAT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER NOBS,N,F(N),T,INDEX,I,J,K,L,M,NP,PCOUNTS(NP),UNK,UNKINIT
	   INTEGER INDICES(NP,2), ICOL
       DOUBLE PRECISION PMAT(N,T,NOBS,M), P(N,(M-1)*(T-1))
       DOUBLE PRECISION DELTA(N,NP*T), DMAT(NOBS,M), COLSUM(M)
C      Zero out values       
       DO 1 I=1,N
       DO 1 J=1,T
       DO 1 K=1,NOBS
       DO 1 L=1,M
          PMAT(I,J,K,L)=0.0D0
  1    CONTINUE
       DO 2 K=1,NOBS
       DO 2 L=1,M
          DMAT(K,L)=0.0D0
  2    CONTINUE
C      Loop over each capture history and from F(I) to T
       DO 20 I=1,N
       DO 20 J=F(I),T
C         For first release occasion, p=1 and delta only applied if for initial
C         release state could be uncertainty
          IF(J.EQ.F(I)) THEN
                DO 5 K=1,NP
				  ICOL=INDICES(K,2)
                  PMAT(I,J,INDICES(K,1),ICOL)=1.0D0 
  5             CONTINUE
                PMAT(I,J,1,M)=1.0D0
                IF(UNKINIT.EQ.1) THEN
                   DO 9 L=1,M-1
		              COLSUM(L)=0.0D0
  9                CONTINUE
                   INDEX=(J-1)*NP
                   DO 10 K=1,NP
                      INDEX=INDEX+1
                      ICOL=INDICES(K,2)
                      DMAT(INDICES(K,1),ICOL)=DELTA(I,INDEX)
			          COLSUM(ICOL)=COLSUM(ICOL)+DELTA(I,INDEX)
 10                CONTINUE
                   DO 12 L=1,M-1
  		           DO 12 K=2,NOBS
			           DMAT(K,L)=DMAT(K,L)/COLSUM(L)
			           PMAT(I,J,K,L)=PMAT(I,J,K,L)*DMAT(K,L)
 12                CONTINUE
                ENDIF             
           ELSE
C               For remaining occasions compute matrix with p values
C               Uses M-1 to exclude death state
                INDEX=(J-2)*(M-1)
                DO 13 K=1,NP
                   ICOL=INDICES(K,2)
                   PMAT(I,J,INDICES(K,1),ICOL)=P(I,INDEX+ICOL)
                   PMAT(I,J,1,ICOL)=1-P(I,INDEX+ICOL)
  13             CONTINUE
                PMAT(I,J,1,M)=1.0D0
C               Compute values of delta IF there is uncertainty (UNK=1)
C               Uses M-1 to exclude death state
                IF(UNK.EQ.1) THEN
		           DO 14 L=1,M-1
		              COLSUM(L)=0.0D0
 14                CONTINUE
                   INDEX=(J-1)*NP
                   DO 15 K=1,NP
                      INDEX=INDEX+1
                      ICOL=INDICES(K,2)
                      DMAT(INDICES(K,1),ICOL)=DELTA(I,INDEX)
			          COLSUM(ICOL)=COLSUM(ICOL)+DELTA(I,INDEX)
 15                CONTINUE
                   DO 17 L=1,M-1
  		           DO 17 K=2,NOBS
			           DMAT(K,L)=DMAT(K,L)/COLSUM(L)
			           PMAT(I,J,K,L)=PMAT(I,J,K,L)*DMAT(K,L)
 17                CONTINUE
                ENDIF
           ENDIF
 20    CONTINUE
       RETURN
       END
