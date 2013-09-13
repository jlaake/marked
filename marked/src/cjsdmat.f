       SUBROUTINE CJSP(P,N,F,T,PMAT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER N,F(N),T
       DOUBLE PRECISION PMAT(N,T,2,2), P(N,T-1)
       DO 2 I=1,N
       DO 2 J=1,T-1
       IF(J.LE.F(I)) THEN
            IF(J.EQ.F(I)) THEN
                PMAT(I,J,1,1)=0.0D0
                PMAT(I,J,1,2)=1.0D0
                PMAT(I,J,2,1)=1.0D0
                PMAT(I,J,2,2)=0.0D0
            ELSE
                PMAT(I,J,1,1)=0.0D0
                PMAT(I,J,1,2)=0.0D0
                PMAT(I,J,2,1)=0.0D0
                PMAT(I,J,2,2)=0.0D0            
            ENDIF
       ENDIF
       IF(J.GE.F(I)) THEN
            PMAT(I,J+1,1,1)=1-P(I,J)
            PMAT(I,J+1,1,2)=1.0D0
            PMAT(I,J+1,2,1)=P(I,J)
            PMAT(I,J+1,2,2)=0.0D0
       ENDIF
  2    CONTINUE
       RETURN
       END
  