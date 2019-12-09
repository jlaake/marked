       SUBROUTINE UMS2P(P,DELTA,N,MA,MS,F,T,PMAT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER N,F(N),T,INDEX,I,J,K,L,MA,MS,M
       DOUBLE PRECISION PMAT(N,T,MA*(MS+1)+1,MS*MA+1)
       DOUBLE PRECISION P(N,MA*MS*(T-1))
       DOUBLE PRECISION DELTA(N,MS*(T-1))
C      Zero out values       
       DO 4 I=1,N
       DO 3 J=1,T
       DO 2 K=1,MA*(MS+1)+1
       DO 1 L=1,MS*MA+1
          PMAT(I,J,K,L)=0.0D0
 1    CONTINUE
 2    CONTINUE
 3    CONTINUE
 4    CONTINUE
C      Loop over each capture history and F(I) to T-1 occasions
       DO 25 I=1,N
       DO 20 J=F(I),T-1
C         For first release occasion p=1
          IF(J.EQ.F(I)) THEN
              IAREA=0
              DO 6 K=2,MA*(MS+1)+1,MS+1
                 IAREA=IAREA+1
                 DO 5 L=1,MS
                     M=(IAREA-1)*MS
                     PMAT(I,J,K+L-1,M+L)=1.0D0 
  5           CONTINUE
  6           CONTINUE
              PMAT(I,J,1,MA*MS+1)=1.0D0
           ENDIF
C          For each remaining occasion, compute matrix with p and delta values
C          Each area value has an uncertain state observation.
           INDEX=(J-1)*MS*MA
           IAREA=0
           DO 9 K=2,MA*(MS+1)+1,MS+1
              IAREA=IAREA+1
              DO 8 L=1,MS
                  INDEX=INDEX+1
                  M=(IAREA-1)*MS
                  PMAT(I,J+1,K+L-1,M+L)=P(I,INDEX)*DELTA(I,INDEX) 
  8        CONTINUE
  9        CONTINUE
           INDEX=(J-1)*MS*MA
           IAREA=0
           DO 11 K=MS+2,MA*(MS+1)+1,MS+1
           IAREA=IAREA+1
           DO 10 L=(IAREA-1)*MS+1,IAREA*MS
              INDEX=INDEX+1
              PMAT(I,J+1,K,L)=P(I,INDEX)*(1-DELTA(I,INDEX))
  10       CONTINUE
  11       CONTINUE
           INDEX=(J-1)*MS*MA
           DO 12 L=1,MA*MS
              INDEX=INDEX+1
              PMAT(I,J+1,1,L)=1-P(I,INDEX) 
 12        CONTINUE
           PMAT(I,J+1,1,MS*MA+1)=1.0D0
 20    CONTINUE
 25    CONTINUE
       RETURN
       END
  