// TMB Version: Fixed-effect Multivariate Multi-State Cormack-Jolly-Seber model with state uncertainty based on work of Johnson et al 2015.
// Fixed-effect Multi-State Cormack-Jolly-Seber model with unobservable states
// Jeff Laake; 9 Aug 2018

#include <TMB.hpp>                              // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_INTEGER(n);                            // number of capture histories
    DATA_INTEGER(m);                            // number of capture occasions
	  DATA_INTEGER(nS);                           // number of states excluding death state
	  DATA_INTEGER(nobs);                         // number of possible observations including 0
    DATA_IMATRIX(ch);                           // capture history matrix); uses numeric values for observations (1=0 not seen)
    DATA_IVECTOR(frst);                         // occasion first seen for each history
	  DATA_VECTOR(freq);                          // frequency of each capture history
	  DATA_MATRIX(tint);                          // time interval between occasions for each history-interval
 
    DATA_INTEGER(nrowphi);                      // number of rows in the simplified design matrix for Phi - survival
    DATA_MATRIX(phidm);                         // design matrix for Phi
    DATA_VECTOR(phifix);                        // phi fixed values
    DATA_IVECTOR(phiindex);                     // phi indices
    
    DATA_INTEGER(nrowp);                        // number of rows in the simplified design matrix for p
    DATA_MATRIX(pdm);                           // design matrix for p
    DATA_VECTOR(pfix);                          // p fixed values
    DATA_IVECTOR(pindex);                       // p indices

    DATA_INTEGER(npos);                         // number of positions in dmat containing p and delta
    DATA_IMATRIX(ipos);                         // 2 column matrix of row and col values of positions for p and delta    
	
    DATA_INTEGER(nrowpsi);                      // number of rows in the simplified design matrix for psi
    DATA_MATRIX(psidm);                         // design matrix for psi
    DATA_VECTOR(psifix);                        // psi fixed values
    DATA_IVECTOR(psiindex);                     // psi indices

    DATA_INTEGER(nrowd);                        // number of rows in the simplified design matrix for delta
    DATA_MATRIX(ddm);                           // design matrix for delta
    DATA_VECTOR(dfix);                          // delta fixed values
    DATA_IVECTOR(dindex);                       // delta indices
   
    DATA_INTEGER(nrowpi);                       // number of rows in the simplified design matrix for pi
    DATA_MATRIX(pidm);                          // design matrix for pi
    DATA_VECTOR(pifix);                         // pi fixed values
    DATA_IVECTOR(piindex);                      // pi indices

    DATA_INTEGER(initknown);                    // if 1 then delta not used on release occasion
    DATA_INTEGER(debug);                        // if 1 then write out parameters and likelihood at each iteration.

    PARAMETER_VECTOR(phibeta);                  // parameter vector for Phi
    PARAMETER_VECTOR(pbeta);                    // parameter vector for p
    PARAMETER_VECTOR(dbeta);                    // parameter vector for delta
    PARAMETER_VECTOR(psibeta);                  // parameter vector for psi
    PARAMETER_VECTOR(pibeta);                   // parameter vector for pi

    int nrows;                                  // number of entries in design matrix m-1 values for each of nS states
    nrows=nS*(m-1);   

    int nT;                                     // number of transitions excluding death
    nT=nS*nS*(m-1);                       
    int dcells;                                 // number of cells in delta to sum over
    dcells=npos/nS;              
    int i,j,k,bindex,bindex2,k2,k3;             // indices and counters
    int bindex3,bindex4,index,irow,icol;

    vector<Type> uniquephi(nrowphi);   // all unique phi values    
    vector<Type> phi(nrows);           // temp vector for Phis for an individual
    vector<Type> uniquep(nrowp);       // all unique p values    
    vector<Type> p(nrows);             // temp vector for ps for an individual
    vector<Type> uniquepsi(nrowpsi);   // temp vector for psis 
    Type psisum;                       // sum of psi for each state to normalize with
    array<Type> psi(m-1,nS,nS);        // matrix for psis for each occasion 
    vector<Type> uniquepi(nrowpi);     // temp vector for pi 
    vector<Type> pi(nS);               // temp vector for pi for an individual 
    Type pisum;                        // sum of pi across states for an individual to normalize with
    vector<Type> uniquedelta(nrowd);   // temp vector for delta 
    matrix<Type> delta(m,npos);        // matrix of delta values by occasion for an individual 
    Type deltasum;                     // sum of delta across possible observations within a state	      
    array<Type> gamma(m-1,nS+1,nS+1);  // transition probability matrices for individual i
    array<Type> pmat(m,nobs,nS+1);     // observation probability matrices for individual i
    array<Type> dmat(m,nobs,nS+1);     // state dependent observation matrices for individual i
    Type u;                            // sum of state probabilities
    vector<Type> pS(nS+1);             // update vector for prob of being in state j=1,nS + death       
    vector<Type> S(nS+1);              // prob of being in state j=1,nS + death for each occasion
    vector<Type> v(nS+1);              // temporary update vector
    Type Lglki;                        // log-likelihood accumulator
    Type g=0;                          // initialize objective function value 
    
    uniquephi=phidm*phibeta;                          // compute phi
    uniquephi=1/(1+exp(-uniquephi));                  
    for(j=0;j<=nrowphi-1;j++)                         
       if(phifix(j)> -0.5)
          uniquephi(j)=phifix(j);                     // fixed value
       
    uniquep=pdm*pbeta;                                // compute p 
    uniquep=1/(1+exp(-uniquep));                      // compute p 
    for(j=0;j<=nrowp-1;j++)                           
       if(pfix(j)> -0.5)
          uniquep(j)=pfix(j);                         // fixed value

    uniquepsi=exp(psidm*psibeta);                     // compute exp of psidm*psibeta; these are components of psi
    for(j=0;j<=nrowpsi-1;j++)
    {
       if(psifix(j) > -0.5)
           uniquepsi(j)=psifix(j);                    // fixed exp Psi value     
    }  

    if(!initknown)uniquepi=exp(pidm*pibeta);          // compute exp of pidm*pibeta; these are components of pi
    for(j=0;j<=nrowpi-1;j++)
    {
       if(pifix(j) > -0.5)
           uniquepi(j)=pifix(j);                      // fixed exp pi value     
    }  

    uniquedelta=exp(ddm*dbeta);                       // compute exp of ddm*dbeta; these are components of delta
    for(j=0;j<=nrowd-1;j++)
    {
       if(dfix(j) > -0.5)
           uniquedelta(j)=dfix(j);                   // fixed exp delta value     
    }  

    for(i=0;i<=n-1;i++)                             // loop over capture histories - one per capture history
    {
        bindex=i*nrows;                             // initialize indices into index values for the ith history
        bindex2=i*nT;
        bindex3=i*nS;    
        bindex4=i*m*npos;
        pisum=0;                                    // compute initial state distribution
        for (k=0;k<=nS-1;k++)                        
          pisum=pisum+uniquepi(piindex(bindex3+k)-1);
        for (k=0;k<=nS-1;k++)                        
          pi(k)=uniquepi(piindex(bindex3+k)-1)/pisum;
                                                          
	      for (j=0;j<=m-2;j++)                       // loop over occasions for the individual computing phi,psi,p
	      {                       
  	       for (k=0;k<=nS-1;k++)                        
           {   
	           p(j*nS+k)=uniquep(pindex(bindex+j*nS+k)-1);
             phi(j*nS+k)=pow(uniquephi(phiindex(bindex+j*nS+k)-1),tint(i,j));     
	           psisum=0;
	           for(k2=0;k2<=nS-1;k2++)
                psisum=psisum+uniquepsi(psiindex(bindex2+k*nS+k2)-1);
  	         for(k2=0;k2<=nS-1;k2++)	             
	             psi(j,k,k2)=uniquepsi(psiindex(bindex2+k*nS+k2)-1)/psisum;
	          }
	          bindex2=bindex2+nS*nS;   
        }        
	      for (j=0;j<=m-1;j++)                       // loop over occasions for delta
	      {                       
  	       for (k=0;k<=nS-1;k++)                        
           {   
	           deltasum=0;
	           for(k2=0;k2<=dcells-1;k2++)
	               deltasum=deltasum+uniquedelta(dindex(bindex4+k*dcells+k2)-1);
	           for(k2=0;k2<=dcells-1;k2++)              
	               delta(j,k*dcells+k2)=uniquedelta(dindex(bindex4+k*dcells+k2)-1)/deltasum;
	        }
	        bindex4=bindex4+npos;      
        }        

	    
	      // ********************************  gamma *************************** 
	      //  compute transition matrices for each occasion

	      gamma.setZero();                               // initialize all transitions to zero
	      index=0;
	      for(j=0;j<=m-2;j++)                            // loop over intervals 
	      {
	         for(k=0;k<=nS-1;k++)                        // loop over states creating p and gamma values
	         {
	           for(k2=0;k2<=nS-1;k2++)
	             gamma(j,k,k2)=psi(j,k,k2)*phi(index);   // adjust psi for survival
	           gamma(j,k,nS)=1-phi(index);               // add death state value for each state
	           index++;
	         }
	         gamma(j,nS,nS)=1;                           // death is an absorbing state
	      }
	    
	    
	      // ********************************  dmat *************************** 
	      //  compute state dependent observation probability matrices for each occasion
	      dmat.setZero();
	      pmat.setZero();
	      for(j=frst(i);j<=m;j++)
	      {
	        bindex=(j-2)*nS;
	        for(k=0;k<=npos-1;k++)
	        {
	          irow=ipos(k,0);
	          icol=ipos(k,1);
	          if(j==frst(i)) 
	          {
	            if(initknown>0)
	              dmat(j-1,irow-1,icol-1)=1;  
	            else
	              dmat(j-1,irow-1,icol-1)=delta(j-1,k);  
	            dmat(j-1,0,icol-1)=0;  
	          } else
	          {	    
	            dmat(j-1,irow-1,icol-1)=p(bindex+icol-1)*delta(j-1,k);  
	            dmat(j-1,0,icol-1)=1-p(bindex+icol-1);
	          }  
	        }
	        dmat(j-1,0,nS)=1;
	      }
	      
	      //  *****************  HMM algorithm ***********************************************
	      pS.setZero();                                    // initialize values to 0
	      Lglki=0;	
	      S.setZero();                                     // initialize values to 0
	      for(j=0;j<=nS-1;j++)
	        S(j)=pi(j);                                   // set state prob to pi for first observation
	      for(j=frst(i);j<=m;j++)                          // loop over possible occasions from first(i)+1 to m
	      {
	        if(j==frst(i))
	          pS=S;
	        else
	        {
	          for(k=0;k<=nS;k++)
	          {
	            pS(k)=0;
	            for(k2=0;k2<=nS;k2++)
	              pS(k)+= S(k2)*gamma(j-2,k2,k);
	          }
	        }
	        for(k3=0;k3<=nS;k3++)
	        {
	          v(k3)=pS(k3)*dmat(j-1,ch(i,j-1)-1,k3);         // v is temp state vector alpha in Z&M
	        }
	        u=sum(v);                                     // sum across states
	        if(u==0) {
	          Rcout << "\n Check Psi or p values set to 0";
	          Rcout << "\n i = " << i << " ch = " << ch(i);
	        }
	        S=v/u;                                        // update S;S is normalized alpha(j) (phi) in Z&M
	        Lglki+=log(u);    	                           // accumulate log-likelihood value
	      }
	      g-=freq(i)*Lglki;
	      
         if(debug>0)
         {
           Rcout <<  "\nphi " <<phibeta;
           Rcout <<  "\np " <<pbeta;
           Rcout <<  "\ndelta " <<dbeta;
           Rcout << "\npsi " <<psibeta;
           Rcout <<  "\npi " <<pibeta;
           Rcout <<  "\n-lnl = " << g << "\n";
         }
    }
    return g;      
}
