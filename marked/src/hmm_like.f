       SUBROUTINE hmmlike(x,N,M,T,no,start,dmat,gamma,delta,lnl,P)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER N,M,T,no,I,J,first
	   INTEGER x(N,T),start(N,2)
       DOUBLE PRECISION P(M,M)
       DOUBLE PRECISION dmat(N,T,no,M), gamma(N,T-1,M,M), delta(N,M)
       DOUBLE PRECISION LNL 
       LNL=0D0
       DO 10 I=1,N
	       first=start(I,2)
	       DO 1 K=1,M
		   DO 1 L=1,M
		      P(K,L)=0.0D0
		      IF(K.eq.L) P(K,K)=dmat(I,first,x(I,first),K)	
   1       CONTINUE			  
           DO 5 J=start(I,2)+1,M
              LNL=1D0
  5        CONTINUE		    
  10   CONTINUE
       RETURN
       END
