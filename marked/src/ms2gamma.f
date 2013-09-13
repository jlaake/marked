       SUBROUTINE MS2GAM(S,PSI,ALPHA,N,NS,MA,MS,F,T,TMAT)
       INTEGER N,F(N),T,I,J,K,L,INDEX
       DOUBLE PRECISION SMAT(N,T-1,NS,NS), S(N,MA*MS*(T-1)),SURV
       DOUBLE PRECISION PSI(N,MS*MS*(T-1)),TMAT(N,T-1,NS,NS)
       DOUBLE PRECISION PSIMAT(N,T-1,NS,NS),PSIX,ALPHAX
       DOUBLE PRECISION AMAT(MA,MA),ALPHA(N,NA*NA*(T-1))
       DO 10 I=1,N
       DO 10 J=F(I),T-1
           INDEX=(J-1)*MS*MA
           DO 4 K=1,MS*MA
              INDEX=INDEX+1
              SURV=S(I,INDEX)
              DO 3 L=1,MS*MA
                 SMAT(I,J,K,L)=SURV
  3           CONTINUE
              SMAT(I,J,K,NS)=1.0D0-SURV
  4        CONTINUE         
           DO 7 L=1,NS-1
              SMAT(I,J,NS,L)=0.0D0
  7        CONTINUE
           SMAT(I,J,NS,NS)=1.0D0
  10   CONTINUE
       DO 20 I=1,N
       DO 20 J=F(I),T-1
           INDEX=(J-1)*MS*MS
           DO 17 K=1,NS-1,MS
              PSISUM=0.0D0
              DO 13 L=1,M-1
                 INDEX=INDEX+1
                 PSIX=PSI(I,INDEX)
                 PSISUM=PSISUM+PSIX
                 PSIMAT(I,J,K,L)=PSIX
  13          CONTINUE
              PSIMAT(I,J,K,M)=1.0D0
              DO 14 L=1,M-1
                 PSIMAT(I,J,K,L)=PSIMAT(I,J,K,L)/PSISUM
  14          CONTINUE
  17       CONTINUE         
           DO 19 L=1,M
              PSIMAT(I,J,M,L)=1.0D0
  19       CONTINUE
  20   CONTINUE
       DO 30 I=1,N
       DO 30 J=1,T-1
       IF(J.LT.F(I)) THEN
           DO 22 K=1,M
           DO 22 L=1,M
             TMAT(I,J,K,L)=0.0D0
  22       CONTINUE
       ELSE
           DO 25 K=1,M
           DO 25 L=1,M
              TMAT(I,J,K,L)=PSIMAT(I,J,K,L)*SMAT(I,J,K,L)
  25       CONTINUE
       ENDIF
  30   CONTINUE
       RETURN
       END
  