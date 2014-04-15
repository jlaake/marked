       SUBROUTINE CJS1TLGAM(PHI,TAU,N,F,T,TMAT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER N,F(N),T,K,L,I,FIRST
       DOUBLE PRECISION TMAT(N,T-1,3,3), PHI(N,T-1)
       DOUBLE PRECISION TAUVEC(2), TAU(N,T-1)
C      Loop over each capture history and from F(I) to T-1 occasions   
       DO 5 I=1,N
       DO 5 J=F(I),T-1
           TAUVEC(1)=1.0D0-TAU(I,J)
           TAUVEC(2)=1.0D0-TAUVEC(1)
           DO 1 K=1,3
           DO 1 L=1,3
              TMAT(I,J,K,L)=0.0D0
  1        CONTINUE
           DO 2 L=1,2
              TMAT(I,J,1,L)=PHI(I,J)*TAUVEC(L)
  2        CONTINUE
           DO 3 K=1,2
              TMAT(I,J,K,3)=1-PHI(I,J)
  3        CONTINUE
           TMAT(I,J,2,2)=PHI(I,J)
           TMAT(I,J,3,3)=1.0D0
  5    CONTINUE
       RETURN
       END
  