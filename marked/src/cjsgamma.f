       SUBROUTINE CJSGAM(PHI,N,F,T,PHIMAT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER N,F(N),T
       DOUBLE PRECISION PHIMAT(N,T-1,2,2), PHI(N,T-1)
       DO 2 I=1,N
       DO 2 J=F(I),T-1
	   IF(J.LT.F(I)) THEN
            PHIMAT(I,J,1,1)=0.0D0
            PHIMAT(I,J,1,2)=0.0D0
            PHIMAT(I,J,2,1)=0.0D0
            PHIMAT(I,J,2,1)=0.0D0
       ELSE
            PHIMAT(I,J,1,1)=PHI(I,J)
            PHIMAT(I,J,1,2)=1-PHI(I,J)
            PHIMAT(I,J,2,1)=0.0D0
            PHIMAT(I,J,2,2)=1.0D0
       ENDIF
  2    CONTINUE
       RETURN
       END
  