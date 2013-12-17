       SUBROUTINE MS2GAM(S,PSI,ALPHA,N,NS,MA,MS,F,T,TMAT)
       INTEGER N,F(N),T,I,J,K,L,INDEX
       DOUBLE PRECISION S(N,MA*MS*(T-1)),SURV
       DOUBLE PRECISION PSI(N,MS*MS*(T-1)),TMAT(N,T-1,NS,NS)
       DOUBLE PRECISION PSIMAT(MS,MS),PSIX,ALPHAX,ROWSUM
       DOUBLE PRECISION AMAT(MA,MA),ALPHA(N,MA*MA*(T-1))
C      Loop over each capture history and T-1 occasions       
       DO 10 I=1,N
       DO 10 J=F(I),T-1
           INDEX=(J-1)*MS*MA
           DO 4 K=1,MS*MA
              INDEX=INDEX+1
              SURV=S(I,INDEX)
              DO 3 L=1,MS*MA
                 TMAT(I,J,K,L)=SURV
  3           CONTINUE
              TMAT(I,J,K,NS)=1.0D0-SURV
  4        CONTINUE         
           DO 7 L=1,NS-1
              TMAT(I,J,NS,L)=0.0D0
  7        CONTINUE
           TMAT(I,J,NS,NS)=1.0D0
  10   CONTINUE
       DO 20 I=1,N
       DO 20 J=F(I),T-1
           INDEX=(J-1)*MS*MS
           DO 13 K=1,MS
              ROWSUM=0.0D0
              DO 11 L=1,MS
                 INDEX=INDEX+1
                 PSIX=PSI(I,INDEX)
                 ROWSUM=ROWSUM+PSIX
                 PSIMAT(K,L)=PSIX
  11          CONTINUE
              DO 12 L=1,MS
                 PSIMAT(K,L)=PSIMAT(K,L)/ROWSUM
  12          CONTINUE    
  13       CONTINUE         
           INDEX=(J-1)*MA*MA
           DO 17 K=1,MA
              ROWSUM=0.0D0
              DO 15 L=1,MA
                 INDEX=INDEX+1
                 ALPHAX=ALPHA(I,INDEX)
                 ROWSUM=ROWSUM+ALPHAX
                 AMAT(K,L)=ALPHAX
  15          CONTINUE
              DO 16 L=1,MA
                 AMAT(K,L)=AMAT(K,L)/ROWSUM
  16          CONTINUE    
  17       CONTINUE         
           DO 18 K=1,NS-1
           DO 18 L=1,NS-1
              INDEXA1=INT((K-1)/MS)+1
              INDEXA2=INT((L-1)/MS)+1
              INDEXS1=K-INT((K-1)/MS)*MS
              INDEXS2=L-INT((L-1)/MS)*MS
              TMAT(I,J,K,L)=TMAT(I,J,K,L)*
     *          AMAT(INDEXA1,INDEXA2)*PSIMAT(INDEXS1,INDEXS2)
  18       CONTINUE
  20   CONTINUE
       RETURN
       END
  