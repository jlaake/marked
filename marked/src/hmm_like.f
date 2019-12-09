       SUBROUTINE HMMLIKE(X,N,M,T,NO,START,FREQ,DMAT,GAMMA,DELTA,LNL)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER N,M,T,NO,I,J,FIRST
       INTEGER X(N,T),START(N)
       DOUBLE PRECISION V(M),U,PHI(M), FREQ(N)
       DOUBLE PRECISION DMAT(N,T,NO,M), GAMMA(N,T-1,M,M), DELTA(N,M)
       DOUBLE PRECISION LNL(N) 
C      Loop over each capture history
       DO 20 I=1,N
           LNL(I)=0D0
           U=0D0
C          Compute initial state vector phi and likelihood 
C          contribution (if any) for first capture history position
         FIRST=START(I)
         DO 1 K=1,M
             V(K)=DELTA(I,K)*DMAT(I,FIRST,X(I,FIRST),K)
             U=U+V(K)
 1       CONTINUE
         DO 2 K=1,M
             PHI(K)=V(K)/U
 2       CONTINUE
         LNL(I)=DLOG(U)*FREQ(I)
C          Loop over remaining capture history positions, computing
C          log-likelihood value and updating state vector(phi)
         DO 10 J=FIRST+1,T
               DO 6 K=1,M
                  V(K)=0D0
                  DO 5 L=1,M
                    V(K)=V(K)+PHI(L)*GAMMA(I,J-1,L,K)
  5            CONTINUE
  6            CONTINUE
               U=0D0
               DO 7 K=1,M 
                  V(K)=V(K)*DMAT(I,J,X(I,J),K)
                  U=U+V(K)
  7            CONTINUE
               LNL(I)=LNL(I)+DLOG(U)*FREQ(I)
               DO 9 K=1,M
                 PHI(K)=V(K)/U
  9            CONTINUE
  10       CONTINUE
  20   CONTINUE
       RETURN
       END
