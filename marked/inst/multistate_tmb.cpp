// TMB Version: Fixed-effect Multi-State Cormack-Jolly-Seber model with unobservable states
// Jeff Laake; 22 Feb 2017

#include <TMB.hpp>                              // Links in the TMB libraries
 
template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_INTEGER(n);                            // number of capture histories
    DATA_INTEGER(m);                            // number of capture occasions
    DATA_INTEGER(nS);                           // number of states excluding death state
    DATA_IMATRIX(ch);                           // capture history matrix; uses numeric values for states
    DATA_IVECTOR(frst);                         // occasion first seen for each history
    DATA_VECTOR(freq);                          // frequency of each capture history
    DATA_MATRIX(tint);                          // time interval between occasions for each history-interval
    DATA_INTEGER(kphi);                         // number of columns in the design matrix for Phi - survival
    DATA_INTEGER(nrowphi);                      // number of rows in the simplified design matrix for Phi - survival
                                                // last column in each design matrix is either -1 (estimated) or a fixed value
    DATA_MATRIX(phidm);                         // design matrix for Phi
    DATA_VECTOR(phifix);                        // phi fixed values
    DATA_IVECTOR(phiindex);                     // phi indices
    
    DATA_INTEGER(kp);                           // number of columns in the design matrix for p - capture probability
    DATA_INTEGER(nrowp);                        // number of rows in the simplified design matrix for p
    DATA_MATRIX(pdm);                           // design matrix for p
    DATA_VECTOR(pfix);                          // p fixed values
    DATA_IVECTOR(pindex);                       // p indices

    DATA_INTEGER(kpsi);                         // number of columns in Psi design matrix
    DATA_INTEGER(nrowpsi);                       // number of rows in the simplified design matrix for psi
    DATA_MATRIX(psidm);                         // design matrix for psi
    DATA_VECTOR(psifix);                        // psi fixed values
    DATA_IVECTOR(psiindex);                     // psi indices
		
    PARAMETER_VECTOR(phibeta);                  // parameter vector for Phi
    PARAMETER_VECTOR(pbeta);                    // parameter vector for p
    PARAMETER_VECTOR(psibeta);                  // parameter vector for Psi
	
    Type g=0; 

	
    int nrows;                           // number of entries in design matrix m-1 values for each of nS states
    int all_nrows;                       // number of rows in design matrix for all histories
    nrows=nS*(m-1);                
    all_nrows=n*nrows;                
    int nT;                             // number of transitions excluding death
    int all_nT;                         // number of transitions for all histories
    nT=nS*nS*(m-1);                       
    all_nT=n*nT;                       

    int i,j,k,bindex,bindex2,k2;        // indices and counters
    vector<Type> uniquephi(nrowphi);    // all unique phi values    
    vector<Type> phi(nrows);            // temp vector for Phis for an individual
    vector<Type> uniquep(nrowp);        // all unique p values    
    vector<Type> p(nrows);              // temp vector for ps for an individual
    vector<Type> uniquepsi(nrowpsi);    // temp vector for psis 
    Type psisum;                        // sum of psi for each state to normalize with

    array<Type> psi(m-1,nS,nS);         // matrix for psis for each occasion 
	array<Type> gamma(m-1,nS+1,nS+1);   // transition probability matrices for individual i
	array<Type> dmat(m-1,nS+1,nS+1);    // observation probability matrices for individual i
    Type u;                             // sum of state probabilities
	vector<Type> pS(nS+1);              // update vector for prob of being in state j=1,nS + death       
	vector<Type> S(nS+1);               // prob of being in state j=1,nS + death for each occasion
	vector<Type> v(nS+1);               // temporary update vector
	vector<Type> vec;                   // temporary vector
	Type Lglki=0;                       // log-likelihood accumulator
 

    for(j=1;j<=nrowphi;j++)                                // compute all unique phi values
       if(phifix(j-1)< -0.5)
	   {
	      vec=phidm.row(j-1);
          uniquephi(j-1)=1/(1+exp(-((vec*phibeta).sum())));   // compute phi
	   }
       else
          uniquephi(j-1)=phifix(j-1);                      // fixed value

    for(j=1;j<=nrowp;j++)                                  // compute all unique p values
       if(pfix(j-1)< -0.5)
       {
	      vec=pdm.row(j-1);
	      uniquep(j-1)=1/(1+exp(-((vec*pbeta).sum())));         // compute p 
       }
	   else
          uniquep(j-1)=pfix(j-1);                          // fixed value
          
    for(j=1;j<=nrowpsi;j++)
    {
       if(psifix(j-1) < -0.5)
	   {
 	       vec=psidm.row(j-1);
           uniquepsi(j-1)=exp((vec*psibeta).sum());        // compute exp of psidm*psibeta; these are components of psi
	   }
       else
           uniquepsi(j-1)=psifix(j-1);                    // fixed exp Psi value     
    }  
          
    for(i=1;i<=n;i++)                                     // loop over capture histories - one per capture history
    {
                                                          //  compute phi and p values for the individual   
        bindex=(i-1)*nrows;                               // initialize indices into index values for the ith history
        bindex2=(i-1)*nT;    
	    for (j=1;j<=m-1;j++)
	    {                       
  	       for (k=1;k<=nS;k++)                        
               {   
	           p((j-1)*nS+k-1)=uniquep(pindex(bindex+(j-1)*nS+k-1)-1);
                   phi((j-1)*nS+k-1)=pow(uniquephi(phiindex(bindex+(j-1)*nS+k-1)-1),tint(i-1,j-1));           
	           psisum=0;
	           for(k2=1;k2<=nS;k2++)              
		           psisum=psisum+uniquepsi(psiindex(bindex2+(k-1)*nS+k2-1)-1);
  	           for(k2=1;k2<=nS;k2++)
	               psi(j-1,k-1,k2-1)=uniquepsi(psiindex(bindex2+(k-1)*nS+k2-1)-1)/psisum;
	       }
	       bindex2=bindex2+nS*nS;         
        }

        //  compute transition matrices for each occasion
   	    gamma.setZero();                        // initialize all transitions to zero
        bindex=1;
        for(j=1;j<=m-1;j++)                        // loop over intervals 
	    {
	        for(k=1;k<=nS;k++)                     // loop over states creating p and gamma values
		   {
	          for(k2=1;k2<=nS;k2++)
		         gamma(j-1,k-1,k2-1)=psi(j-1,k-1,k2-1)*phi(bindex-1);    // adjust psi for survival
		      gamma(j-1,k-1,nS)=1-phi(bindex-1);              // add death state value for each state
		      bindex++;
		    }
		    gamma(j-1,nS,nS)=1;                             // death is an absorbing state
        }
		

        //  compute state dependent observation matrices for each occasion
	     dmat.setZero();
	     bindex=1;
         for(j=1;j<=m-1;j++)
        {
	       for(k=1;k<=nS;k++)
	       {
             dmat(j-1,k,k-1)=p(bindex-1);  
             dmat(j-1,0,k-1)=1-dmat(j-1,k,k-1);  
		     bindex++;
	       }
	       dmat(j-1,0,nS)=1;
	    }
        //  HMM algorithm
	    pS.setZero();                                    // initialize values to 0
        Lglki=0;	
	    S.setZero();                                     
	    S(ch(i-1,frst(i-1)-1)-1)=1;                        // set state prob to 1 for observed state at first observation
		//if(i==3) 
		//for(k=1;k<=nS+1;k++)
		//  Rcout << "k= " <<k << " S = " << S(k-1);
        for(j=frst(i-1)+1;j<=m;j++)                      // loop over possible occasions from first(i)+1 to m
        {
		    for(k=1;k<=nS+1;k++)
			{
			   pS(k-1)=0;
			   for(k2=1;k2<=nS+1;k2++)
			     pS(k-1)+= S(k2-1)*gamma(j-2,k2-1,k-1);
			}
            //pS=gamma(j-2)*S;        	                 // previous scaled state probability*transition matrix
			//if(i==3) 
			//{
			  //Rcout << frst(i-1);
			  //for(k=1;k<=nS+1;k++)
			  //{
			  //Rcout << "\n";
			  //for(k2=1;k2<=nS+1;k2++)
              //  Rcout << gamma(j-2,k-1,k2-1) << " ";
			  //}
			  //Rcout << "j = " << j << "  " << pS << "\n";
			  
			//}
			for(k=1;k<=nS+1;k++)
			{
              v(k-1)=pS(k-1)*dmat(j-2,ch(i-1,j-1),k-1); // v is temp state vector alpha in Z&M
			 // if(i==3)
             // {
 			 //    Rcout <<  ch(i-1,j-1) << "\n";
 			 //    Rcout <<  "k = " << k << "dmat " << dmat(j-2,ch(i-1,j-1),k-1) << "\n";
 			 //    Rcout <<  "pS " << pS(k-1) << "\n";
			 // }
			}
            if(v.sum()==0) {
              //   Rcout << "\n Check Psi or p values set to 0";
              //   Rcout << "\n i = " << i << " ch = " << ch(i-1);
            }
    	    u=v.sum();                                       // sum across states
            S=v/u;                                          // update S;S is normalized alpha(j) (phi) in Z&M
	        Lglki+=log(u);    	                            // accumulate log-likelihood value
			//if(i==3) Rcout << " u " << u << " i = " << i << "  " << Lglki << "\n";

	    }
	    g-=freq(i-1)*Lglki;
     }
     return g;
}
