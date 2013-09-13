       SUBROUTINE UMS2P(P,DELTA,N,MA,MS,F,T,PMAT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER N,F(N),T,INDEX,I,J,K,L,MA,MS,M
       DOUBLE PRECISION PMAT(N,T,MA*(MS+1)+1,MS*MA+1)
       DOUBLE PRECISION P(N,MA*MS*(T-1))
       DOUBLE PRECISION DELTA(N,MS*(T-1))
       DO 1 I=1,N
       DO 1 J=1,T
       DO 1 K=1,MA*(MS+1)+1
       DO 1 L=1,MS*MA+1
          PMAT(I,J,K,L)=0.0D0
  1    CONTINUE
       DO 20 I=1,N
       DO 20 J=F(I),T-1
          IF(J.EQ.F(I)) THEN
              IAREA=0
              DO 5 K=2,MA*(MS+1)+1,MS+1
                 IAREA=IAREA+1
                 DO 5 L=1,MS
                     M=(IAREA-1)*MS
                     PMAT(I,J,K+L-1,M+L)=1.0D0 
  5           CONTINUE
              PMAT(I,J,1,MA*MS+1)=1.0D0
           ENDIF
           INDEX=(J-1)*MS*MA
           IAREA=0
           DO 9 K=2,MA*(MS+1)+1,MS+1
              IAREA=IAREA+1
              DO 9 L=1,MS
                  INDEX=INDEX+1
                  M=(IAREA-1)*MS
                  PMAT(I,J+1,K+L-1,M+L)=P(I,INDEX)*DELTA(I,INDEX) 
  9        CONTINUE
           INDEX=(J-1)*MS*MA
           IAREA=0
           DO 10 K=MS+2,MA*(MS+1)+1,MS+1
           IAREA=IAREA+1
           DO 10 L=(IAREA-1)*MS+1,IAREA*MS
              INDEX=INDEX+1
              PMAT(I,J+1,K,L)=P(I,INDEX)*(1-DELTA(I,INDEX))
  10       CONTINUE
           INDEX=(J-1)*MS*MA
           DO 11 L=1,MA*MS
              INDEX=INDEX+1
              PMAT(I,J+1,1,L)=1-P(I,INDEX) 
 11        CONTINUE
           PMAT(I,J+1,1,MS*MA+1)=1.0D0
 20    CONTINUE
       RETURN
       END
  