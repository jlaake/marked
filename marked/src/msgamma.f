       SUBROUTINE MSGAM(S,PSI,N,M,F,T,TMAT)
       INTEGER N,F(N),T,I,J,K,L,INDEX,INDEX1
       DOUBLE PRECISION S(N,(M-1)*(T-1)),SURV
       DOUBLE PRECISION PSI(N,(M-1)*(M-1)*(T-1)),TMAT(N,T-1,M,M)
       DOUBLE PRECISION PSIMAT(M,M),PSIX,PSISUM
C      Loop over each capture history and T-1 intervals
       DO 30 I=1,N
	   IF(F(I).LT.T) THEN
       DO 27 J=1,T-1
C        Zero out values for intervals before first
         IF(J.LT.F(I)) THEN
           DO 1 K=1,M
           DO 1 L=1,M
             TMAT(I,J,K,L)=0.0D0
 1         CONTINUE
        ELSE   
C          For remaining intervals compute Psi matrix
           INDEX=(J-1)*(M-1)
           DO 4 K=1,M-1
              INDEX=INDEX+1
              SURV=S(I,INDEX)
              DO 3 L=1,M-1
                 TMAT(I,J,K,L)=SURV
  3           CONTINUE
              TMAT(I,J,K,M)=1.0D0-SURV
  4        CONTINUE         
           DO 7 L=1,M-1
              TMAT(I,J,M,L)=0.0D0
  7        CONTINUE
           TMAT(I,J,M,M)=1.0D0
           INDEX1=(J-1)*(M-1)*(M-1)
           DO 17 K=1,M-1
              PSISUM=0.0D0
              DO 13 L=1,M-1
                 INDEX1=INDEX1+1
                 PSIX=PSI(I,INDEX1)
                 PSISUM=PSISUM+PSIX
                 PSIMAT(K,L)=PSIX
  13          CONTINUE
              PSIMAT(K,M)=1.0D0
              DO 14 L=1,M-1
                 PSIMAT(K,L)=PSIMAT(K,L)/PSISUM
  14          CONTINUE
  17       CONTINUE         
           DO 19 L=1,M
              PSIMAT(M,L)=1.0D0
  19       CONTINUE
C          Multiply survival matrix in TMAT and Psi matrix and store in TMAT
           DO 25 K=1,M
           DO 25 L=1,M
              TMAT(I,J,K,L)=PSIMAT(K,L)*TMAT(I,J,K,L)
  25       CONTINUE
       ENDIF
  27   CONTINUE
       ENDIF
  30   CONTINUE
       RETURN
       END
  