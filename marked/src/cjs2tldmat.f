       SUBROUTINE CJS2TLP(P,N,F,T,PMAT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER N,F(N),T,I,J,K,L
       DOUBLE PRECISION PMAT(N,T,5,5), P(N,4*(T-1))
C      Loop over each capture history and T-1 occasions
       DO 5 I=1,N
       DO 5 J=F(I),T-1
C           First occasion is a release occasion, so p=1
            DO 1 K=1,5
            DO 1 L=1,5
              PMAT(I,J+1,K,L)=0.0D0
              IF(J.EQ.F(I)) THEN
                  PMAT(I,J,K,L)=0.0D0
                  IF(K.EQ.L) PMAT(I,J,K,L)=1.0D0
              ENDIF
  1         CONTINUE
C        For each possible recapture occasion create matrix
C        using p and 1-p
         DO 2 L=1,4
            PMAT(I,J+1,5,L)=1-P(I,4*(J-1)+L)
  2      CONTINUE      
         PMAT(I,J+1,5,5)=1.0D0
         DO 3 K=1,4
            PMAT(I,J+1,K,K)=P(I,4*(J-1)+K)
  3      CONTINUE       
  5    CONTINUE
       RETURN
       END
