// Fixed-effect Multi-State Cormack-Jolly-Seber model with unobservable states based on work of Jessica Ford in MEE Dec 2012
// Jeff Laake; 29 June 2016
// Modified from original version to split off code to compute dmat and gamma calculations
// Now allows fixed real parameters and uses simplification
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
    init_int nrowphi;                      // number of rows in the simplified design matrix for Phi - survival
    int nrows;                             // number of entries in design matrix m-1 values for each of nS states
    int all_nrows;                         // number of rows in design matrix for all histories
    !! nrows=nS*(m-1);                
    !! all_nrows=n*nrows;                
                                             // last column in each design matrix is either -1 (estimated) or a fixed value
    init_matrix phidm(1,nrowphi,1,kphi);     // design matrix for Phi
    init_vector phifix(1,nrowphi);           // phi fixed values
    init_ivector phiindex(1,all_nrows);      // phi indices
    
    init_int kp;                             // number of columns in the design matrix for p - capture probability
    init_int nrowp;                          // number of rows in the simplified design matrix for p
    init_matrix pdm(1,nrowp,1,kp);           // design matrix for p
    init_vector pfix(1,nrowp);               // p fixed values
    init_ivector pindex(1,all_nrows);        // p indices

	int nT;                                  // number of transitions excluding death
	int all_nT;                              // number of transitions for all histories
	!! nT=nS*nS*(m-1);                       
	!! all_nT=n*nT;                       
	init_int kpsi;                           // number of columns in Psi design matrix
    init_int nrowpsi;                        // number of rows in the simplified design matrix for psi
    init_matrix psidm(1,nrowpsi,1,kpsi);     // design matrix for psi
    init_vector psifix(1,nrowpsi);           // psi fixed values
    init_ivector psiindex(1,all_nT);         // psi indices
		
PARAMETER_SECTION
    init_vector phibeta(1,kphi);       // parameter vector for Phi
    init_vector pbeta(1,kp);           // parameter vector for p
    init_vector psibeta(1,kpsi);       // parameter vector for p
	objective_function_value g; 

PROCEDURE_SECTION
    int i,j,k,bindex,bindex2,k2;       // indices and counters
    dvar_vector uniquephi(1,nrowphi);  // all unique phi values    
    dvar_vector phi(1,nrows);          // temp vector for Phis for an individual
    dvar_vector uniquep(1,nrowp);      // all unique p values    
    dvar_vector p(1,nrows);            // temp vector for ps for an individual
    dvar_vector uniquepsi(1,nrowpsi);  // temp vector for psis 
    dvariable psisum;                  // sum of psi for each state to normalize with
    dvar3_array psi(1,m-1,1,nS,1,nS);  // matrix for psis for each occasion 

    for(j=1;j<=nrowphi;j++)                           // compute all unique phi values
       if(phifix(j)< -0.5)
          uniquephi(j)=1/(1+exp(-phidm(j)*phibeta));  // compute phi
       else
          uniquephi(j)=phifix(j);                     // fixed value

    for(j=1;j<=nrowp;j++)                             // compute all unique p values
       if(pfix(j)< -0.5)
          uniquep(j)=1/(1+exp(-pdm(j)*pbeta));        // compute p 
       else
          uniquep(j)=pfix(j);                         // fixed value
          
    for(j=1;j<=nrowpsi;j++)
    {
       if(psifix(j) < -0.5)
           uniquepsi(j)=exp(psidm(j)*psibeta);        // compute exp of psidm*psibeta; these are components of psi
       else
           uniquepsi(j)=psifix(j);                    // fixed exp Psi value     
    }  
          
    for(i=1;i<=n;i++)                              // loop over capture histories - one per capture history
    {
                                                   //  compute phi and p values for the individual   
        bindex=(i-1)*nrows;                        // initialize indices into index values for the ith history
        bindex2=(i-1)*nT;    
	    for (j=1;j<=m-1;j++)
	    {                       
  	       for (k=1;k<=nS;k++)                        
           {   
	           p((j-1)*nS+k)=uniquep(pindex(bindex+(j-1)*nS+k));
               phi((j-1)*nS+k)=pow(uniquephi(phiindex(bindex+(j-1)*nS+k)),tint(i,j));           
	           psisum=0;
	           for(k2=1;k2<=nS;k2++)              
		           psisum=psisum+uniquepsi(psiindex(bindex2+(k-1)*nS+k2));
  	           for(k2=1;k2<=nS;k2++)
	               psi(j,k,k2)=uniquepsi(psiindex(bindex2+(k-1)*nS+k2))/psisum;
	       }
	       bindex2=bindex2+nS*nS;         
        }
           
        ll_i(i,phi,p,psi);                 // compute neg log likelihod and increment
     }

FUNCTION dvar3_array get_dmat(const int i, const dvar_vector& p )
	dvar3_array dmat(1,m-1,1,nS+1,1,nS+1); // observation probability matrices for individual i
	int j,k,bindex;
	dmat.initialize();
	bindex=1;
    for(j=1;j<=m-1;j++)
    {
	  for(k=1;k<=nS;k++)
	  {
        dmat(j,k+1,k)=p(bindex);  
        dmat(j,1,k)=1-dmat(j,k+1,k);  
		bindex++;
	  }
	  dmat(j,1,nS+1)=1;
	}
	return dmat;


FUNCTION dvar3_array get_gamma(const int i, const dvar_vector& phi, const dvar3_array& psi )
	int j,k,index,k2;
	dvar3_array gamma(1,m-1,1,nS+1,1,nS+1);    // transition probability matrices for individual i

	gamma.initialize();                        // initialize all transitions to zero
    index=1;
    for(j=1;j<=m-1;j++)                        // loop over intervals 
	{
	    for(k=1;k<=nS;k++)                     // loop over states creating p and gamma values
		{
	       for(k2=1;k2<=nS;k2++)
		     gamma(j,k,k2)=psi(j,k,k2)*phi(index);    // adjust psi for survival
		   gamma(j,k,nS+1)=1-phi(index);              // add death state value for each state
		   index++;
		}
		gamma(j,nS+1,nS+1)=1;                             // death is an absorbing state
    }
	return gamma;

FUNCTION void ll_i(const int i, const dvar_vector& phi, const dvar_vector& p, const dvar3_array& psi )
//  Code implements an algorithm described on pg 47 in Zucchini and MacDonald (Z&M) for computing likelihood value for a 
//  Hidden Markov Model. The algorithm is recursive over the encounter occasions. The algorithm involves a simpler recursion 
//  shown on pages 38 and 45 but to avoid numerical underflows it is scaled by the sum of the state probabilities which makes
//  it appear a little more complex.
	dvar3_array gamma(1,m-1,1,nS+1,1,nS+1); // transition probability matrices for individual i
	dvar3_array dmat(1,m-1,1,nS+1,1,nS+1);  // observation probability matrices for individual i
    dvariable u;                            // sum of state probabilities
	dvar_vector pS(1,nS+1);                 // update vector for prob of being in state j=1,nS + death       
	dvar_vector S(1,nS+1);                  // prob of being in state j=1,nS + death for each occasion
	dvar_vector v(1,nS+1);                  // temporary update vector
	dvariable Lglki;                        // log-likelihood accumulator
    int j, j2;                              // loop variables
   
//  compute transition matrices for each occasion
	gamma=get_gamma(i,phi,psi);
//  compute state dependent observation matrices for each occasion
	dmat=get_dmat(i,p);
//  HMM algorithm
	pS.initialize();                                    // initialize values to 0
    Lglki.initialize();	
	S.initialize();                                     
	S(ch(i,frst(i)))=1;                                 // set state prob to 1 for observed state at first observation
    for(j=frst(i)+1;j<=m;j++)                           // loop over possible occasions from first(i)+1 to m
    {
// 	    j2=j+frst(i)-1; 
        pS=S*gamma(j-1);        	                    // previous scaled state probability*transition matrix
        v=elem_prod(pS,dmat(j-1,ch(i,j)+1));            // v is temp state vector alpha in Z&M
    	u=sum(v);                                       // sum across states
        S=v/u;                                          // update S;S is normalized alpha(j) (phi) in Z&M
	    Lglki+=log(u);    	                            // accumulate log-likelihood value
	}
	g-=freq(i)*Lglki;
