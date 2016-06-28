// Fixed-effect Multi-State Cormack-Jolly-Seber model with unobservable states based on work of Jessica Ford in MEE Dec 2012
// Jeff Laake; 23 June 2016
// Modified from original version to split off code to compute dmat and gamma calculations
// Need to add fixed values to handle unobservable state
DATA_SECTION 
    init_int n;                            // number of capture histories
    init_int m;                            // number of capture occasions
	init_int nS;                           // number of states excluding death state
	init_int uS;                           // number of unobserved states
    init_imatrix ch(1,n,1,m);              // capture history matrix; uses numeric values for states
    init_ivector frst(1,n);                // occasion first seen for each history
	init_vector freq(1,n);                 // frequency of each capture history
    init_matrix tint(1,n,1,m-1);           // time interval between occasions for each history-interval
    init_int kphi;                         // number of columns in the design matrix for Phi - survival
    int nrows;                             // number of entries in design matrix m-1 values for each of nS states
    int all_nrows;                         // number of rows in design matrix for all histories
    !! nrows=nS*(m-1);                
    !! all_nrows=n*nrows;                
                                             // last column in each design matrix is either -1 (estimated) or a fixed value
    init_matrix phidm(1,all_nrows,1,kphi);   // design matrix for Phi
    init_vector phifix(1,all_nrows);         // phi fixed values
    init_int kp;                             // number of columns in the design matrix for p - capture probability
    init_matrix pdm(1,all_nrows,1,kp);       // design matrix for p
    init_vector pfix(1,all_nrows);           // p fixed values
	int nT;                                  // number of transitions excluding death
	int all_nT;                              // number of transitions for all histories
	!! nT=nS*nS*(m-1);                       
	!! all_nT=n*nT;                       
	init_int kpsi;                           // number of columns in Psi design matrix
	init_matrix psidm(1,all_nT,1,kpsi);      // design matrix for transition probs
    init_vector psifix(1,all_nT);            // Psi fixed values
		
PARAMETER_SECTION
    init_vector phibeta(1,kphi);       // parameter vector for Phi
    init_vector pbeta(1,kp);           // parameter vector for p
    init_vector psibeta(1,kpsi);       // parameter vector for p
	objective_function_value g; 

PROCEDURE_SECTION
    int i;                             // index over observations
    for(i=1;i<=n;i++)                  // loop over capture histories - one per capture history
        ll_i(i,phibeta,pbeta,psibeta);

FUNCTION dvar3_array get_dmat(const int i, const dvar_vector& pbeta )
	dvar3_array dmat(1,m-1,1,nS+1,1,nS+1); // observation probability matrices for individual i
	int j,bindex,k;
	dmat.initialize();
	bindex=(i-1)*nrows+1;                  // initialize index into p for ith history
    for(j=1;j<=m-1;j++)
	{
	  for(k=1;k<=nS;k++)
	  {
	     if(pfix(bindex)==-1)
	        dmat(j,k+1,k)=1/(1+exp(-pdm(bindex)*pbeta));  
	     else
	        dmat(j,k+1,k)=pfix(bindex);  	     
	     dmat(j,1,k)=1-dmat(j,k+1,k);  
		 bindex++;
	  }
	  dmat(j,1,nS+1)=1;
	}
	return dmat;


FUNCTION dvar3_array get_gamma(const int i,  const dvar_vector& phibeta, const dvar_vector& psibeta )
	int j,bindex,k,index,index2,k2;
    dvar_vector phi(1,nrows);                  // temp vector for Phis for an occasion
    dvar_vector psiexp(1,nT);                  // temp vector for psis 
    dvariable psisum;                          // sum of psiexp for each state
    dvariable u;                               // sum of state probabilities
    dvar_matrix psi(1,nS,1,nS);                // matrix for psis 
	dvar3_array gamma(1,m-1,1,nS+1,1,nS+1);    // transition probability matrices for individual i
	int ni=m-frst(i)+1;                	       // length of capture history for individual i
    
	bindex=(i-1)*nrows;                        // initialize index into phi for ith history
    for(j=1;j<=nrows;j++)
	   	if(phifix(bindex+j)==-1)
	      phi(j)=1/(1+exp(-phidm(bindex+j)*phibeta));   // compute phi for the interval
		else
		  phi(j)=phifix(bindex+j);                      // fixed real phi  
	for (j=1;j<=m-1;j++)                                
	   if(phifix(bindex+j)==-1)
   	   for (k=1;k<=nS;k++)
          phi((j-1)*nS+k)=pow(phi((j-1)*nS+k),tint(i,j)); // adjust phi for the time interval length if needed
   
    bindex=(i-1)*nT;                                    // initialize index into psi for ith history
    for(j=1;j<=nT;j++)
    {
       if(psifix(bindex+j)==-1)
           psiexp(j)=exp(psidm(bindex+j)*psibeta);      // compute exp of psidm*psibeta; these are components of psi
       else
           psiexp(j)=psifix(bindex+j);                  // fixed exp Psi value     
    }  
	gamma.initialize();                                 // initialize all transitions to zero
    for(j=2;j<=ni;j++)                                  // loop over remaining possible occasions 
	{
 	    int j2=j+frst(i)-1;
 	    index=(j2-2)*nS+1;                              // index into phi vector for each occasion
 	    index2=(j2-2)*nS*nS;                            // index in psi vector
	    for(k=1;k<=nS;k++)                              // loop over states computing psi matrix 
	    {
	       psisum=0;
	       for(k2=(k-1)*nS+1;k2<=k*nS;k2++)              
		     psisum=psisum+psiexp(index2+k2);
  	       for(k2=1;k2<=nS;k2++)
	          psi(k,k2)=psiexp((k-1)*nS+k2+index2)/psisum;
	    }
	    for(k=1;k<=nS;k++)                              // loop over states creating p and gamma values
		{
	       for(k2=1;k2<=nS;k2++)
		     gamma(j-1,k,k2)=psi(k,k2)*phi(index);          // adjust psi for survival
		   gamma(j-1,k,nS+1)=1-phi(index);                  // add death state value for each state
		   index++;
		}
		gamma(j-1,nS+1,nS+1)=1;                             // death is an absorbing state
    }
	return gamma;

FUNCTION void ll_i(const int i, const dvar_vector& phibeta, const dvar_vector& pbeta, const dvar_vector& psibeta )
//  Code implements an algorithm described on pg 47 in Zucchini and MacDonald (Z&M) for computing likelihood value for a 
//  Hidden Markov Model. The algorithm is recursive over the encounter occasions. The algorithm involves a simpler recursion 
//  shown on pages 38 and 45 but to avoid numerical underflows it is scaled by the sum of the state probabilities which makes
//  it appear a little more complex.
	dvar3_array gamma(1,m-1,1,nS+1,1,nS+1); // transition probability matrices for individual i
	dvar3_array dmat(1,m-1,1,nS+1,1,nS+1);  // observation probability matrices for individual i
    dvariable u;                            // sum of state probabilities
	int ni=m-frst(i)+1;                     // length of capture history for individual i
	dvar_vector pS(1,nS+1);                 // update vector for prob of being in state j=1,nS + death       
	dvar_vector S(1,nS+1);                  // prob of being in state j=1,nS + death for each occasion
	dvar_vector v(1,nS+1);                  // temporary update vector
	dvariable Lglki;                        // log-likelihood accumulator
    int j, j2;                              // loop variables
//  compute transition matrices for each occasion
	gamma=get_gamma(i,phibeta,psibeta);
//  compute state dependent observation matrices for each occasion
	dmat=get_dmat(i,pbeta);
//  HMM algorithm
	pS.initialize();                                    // initialize values to 0
    Lglki.initialize();	
	S.initialize();                                     
	S(ch(i,frst(i)))=1;                                 // set state prob to 1 for observed state at first observation
    for(j=2;j<=ni;j++)                                  // loop over remaining possible occasions from 2 to ni
    {
 	    j2=j+frst(i)-1; 
        pS=S*gamma(j-1);        	                    // previous scaled state probability*transition matrix
        v=elem_prod(pS,dmat(j2-1,ch(i,j2)+1));          // v is temp state vector alpha in Z&M
    	u=sum(v);                                       // sum across states
        S=v/u;                                          // update S;S is normalized alpha(j) (phi) in Z&M
	    Lglki+=log(u);    	                            // accumulate log-likelihood value
	}
	g-=freq(i)*Lglki;
