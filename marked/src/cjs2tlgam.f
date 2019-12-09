       SUBROUTINE CJS2TLGAM(PHI,TAU,N,F,T,TMAT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER N,F(N),T,K,L,I,FIRST
       DOUBLE PRECISION TMAT(N,T-1,5,5), PHI(N,T-1)
       DOUBLE PRECISION TAUVEC(4), TAU0(2), TAU(N,4*(T-1)) 
C      Loop over each capture history and from F(I) to T-1 occasions   
       DO 10 I=1,N
       DO 5 J=F(I),T-1
           FIRST=(J-1)*4+1
           TAUSUM=TAU(I,FIRST)+TAU(I,FIRST+1)+TAU(I,FIRST+2)
           TAUSUM=TAUSUM+TAU(I,FIRST+3)
           TAUVEC(1)=TAU(I,FIRST)/TAUSUM
           TAUVEC(2)=TAU(I,FIRST+1)/TAUSUM
           TAUVEC(3)=TAU(I,FIRST+2)/TAUSUM
           TAUVEC(4)=TAU(I,FIRST+3)/TAUSUM
           IF((TAUVEC(2)+TAUVEC(4))>0) THEN
              TAU0(1)=TAUVEC(4)/(TAUVEC(2)+TAUVEC(4))
           ELSE
              TAU0(1)=0.0D0
           ENDIF
           IF((TAUVEC(3)+TAUVEC(4))>0) THEN
              TAU0(2)=TAUVEC(4)/(TAUVEC(3)+TAUVEC(4))
           ELSE
              TAU0(2)=0.0D0
           ENDIF
           DO 2 K=1,5
           DO 1 L=1,5
              TMAT(I,J,K,L)=0.0D0
  1        CONTINUE
  2        CONTINUE
           DO 3 L=1,4
              TMAT(I,J,1,L)=PHI(I,J)*TAUVEC(L)
  3        CONTINUE
           DO 4 K=1,4
              TMAT(I,J,K,5)=1-PHI(I,J)
  4        CONTINUE
           TMAT(I,J,4,4)=PHI(I,J)
           TMAT(I,J,5,5)=1.0D0
           TMAT(I,J,2,2)=PHI(I,J)*(1-TAU0(1))
           TMAT(I,J,2,4)=PHI(I,J)*TAU0(1)
           TMAT(I,J,3,3)=PHI(I,J)*(1-TAU0(2))
           TMAT(I,J,3,4)=PHI(I,J)*TAU0(2)
  5    CONTINUE
 10    CONTINUE
       RETURN
       END
  