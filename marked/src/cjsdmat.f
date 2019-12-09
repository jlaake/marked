       SUBROUTINE CJSP(P,N,F,T,PMAT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER N,F(N),T
       DOUBLE PRECISION PMAT(N,T,2,2), P(N,T-1)
C      Loop over each capture history and T-1 occasions
       DO 3 I=1,N
       DO 2 J=F(I),T-1
C      First occasion is a release occasion, so p=1
       IF(J.EQ.F(I)) THEN
            PMAT(I,J,1,1)=0.0D0
            PMAT(I,J,1,2)=1.0D0
            PMAT(I,J,2,1)=1.0D0
            PMAT(I,J,2,2)=0.0D0
       ENDIF
C      For each possible recapture occasion create matrix
C      using p and 1-p
       PMAT(I,J+1,1,1)=1-P(I,J)
       PMAT(I,J+1,1,2)=1.0D0
       PMAT(I,J+1,2,1)=P(I,J)
       PMAT(I,J+1,2,2)=0.0D0
  2    CONTINUE
  3    CONTINUE
       RETURN
       END
  