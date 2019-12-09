       SUBROUTINE CJS1TLP(P,N,F,T,PMAT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER N,F(N),T,I,J,K,L
       DOUBLE PRECISION PMAT(N,T,3,3), P(N,2*(T-1))
C      Loop over each capture history and T-1 occasions
       DO 10 I=1,N
       DO 5 J=F(I),T-1
C           First occasion is a release occasion, so p=1
            DO 2 K=1,3
            DO 1 L=1,3
              PMAT(I,J+1,K,L)=0.0D0
              IF(J.EQ.F(I)) THEN
                  PMAT(I,J,K,L)=0.0D0
                  IF(K.EQ.L) PMAT(I,J,K,L)=1.0D0
              ENDIF
  1         CONTINUE
  2         CONTINUE
C        For each possible recapture occasion create matrix
C        using p and 1-p
         DO 3 L=1,2
            PMAT(I,J+1,3,L)=1-P(I,2*(J-1)+L)
  3      CONTINUE      
         PMAT(I,J+1,3,3)=1.0D0
         DO 4 K=1,2
            PMAT(I,J+1,K,K)=P(I,2*(J-1)+K)
  4      CONTINUE       
  5    CONTINUE
  10   CONTINUE
       RETURN
       END
      