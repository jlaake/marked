       SUBROUTINE CJSGAM(PHI,N,F,T,PHIMAT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER N,F(N),T
       DOUBLE PRECISION PHIMAT(N,T-1,2,2), PHI(N,T-1)
C      Loop over each capture history and from F(I) to T-1 occasions   
       DO 2 I=1,N
       DO 2 J=F(I),T-1
C      For intervals past initial release occasion, create
C      matrix using Phi and 1-Phi values
          PHIMAT(I,J,1,1)=PHI(I,J)
          PHIMAT(I,J,1,2)=1-PHI(I,J)
          PHIMAT(I,J,2,1)=0.0D0
          PHIMAT(I,J,2,2)=1.0D0
  2    CONTINUE
       RETURN
       END
  