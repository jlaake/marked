// Fixed-effect Multivariate Multi-State Cormack-Jolly-Seber model with state uncertainty based on work of Johnson et al 2015.
// Jeff Laake; 29 June 2016
DATA_SECTION 
    init_int n;                            // number of capture histories
    init_int m;                            // number of capture occasions
	init_int nS;                           // number of states excluding death state
	init_int nobs;                         // number of possible observations including 0
	init_int nd;                           // number of records in delta design data per id-occ
    init_imatrix ch(1,n,1,m);              // capture history matrix; uses numeric values for observations (1=0 not seen)
    init_ivector frst(1,n);                // occasion first seen for each history
	init_vector freq(1,n);                 // frequency of each capture history
    init_matrix tint(1,n,1,m-1);           // time interval between occasions for each history-interval
    init_int kphi;                         // number of columns in the design matrix for Phi - survival
    init_int nrowphi;                      // number of rows in the simplified design matrix for Phi - survival
    int nrows;                             // number of entries in design matrix m-1 values for each of nS states
    int all_nrows;                         // number of rows in design matrix for all histories
    int ndrows;                            // number of entries in delta design matrix per id
    int all_ndrows;                        // number of rows in delta design matrix for all histories
    !! nrows=nS*(m-1);   
    !! ndrows=m*nd;   
    !! all_nrows=n*nrows;                
    !! all_ndrows=n*ndrows;                
    init_matrix phidm(1,nrowphi,1,kphi);     // design matrix for Phi
    init_vector phifix(1,nrowphi);           // phi fixed values
    init_ivector phiindex(1,all_nrows);      // phi indices
    
    init_int kp;                             // number of columns in the design matrix for p - capture probability
    init_int nrowp;                          // number of rows in the simplified design matrix for p
    init_matrix pdm(1,nrowp,1,kp);           // design matrix for p
    init_vector pfix(1,nrowp);               // p fixed values
    init_ivector pindex(1,all_nrows);        // p indices

    init_int npos;                           // number of positions in dmat containing p and delta
    init_imatrix ipos(1,npos,1,2);           // 2 column matrix of row and col values of positions for p and delta    
	
 	int nT;                                  // number of transitions excluding death
	int all_nT;                              // number of transitions for all histories
	!! nT=nS*nS*(m-1);                       
	!! all_nT=n*nT;                       
	init_int kpsi;                           // number of columns in Psi design matrix
    init_int nrowpsi;                        // number of rows in the simplified design matrix for psi
    init_matrix psidm(1,nrowpsi,1,kpsi);     // design matrix for psi
    init_vector psifix(1,nrowpsi);           // psi fixed values
    init_ivector psiindex(1,all_nT);         // psi indices

    int dcells;                              // number of cells in delta to sum over
	!! dcells=npos/nS;              
    init_int kd;                             // number of columns in the design matrix for delta
    init_int nrowd;                          // number of rows in the simplified design matrix for delta

    init_matrix ddm(1,nrowd,1,kd);           // design matrix for delta
    init_vector dfix(1,nrowd);               // delta fixed values

    init_ivector dindex(1,all_ndrows);       // delta indices
   
    init_int kpi;                             // number of columns in the design matrix for pi - initial state distribution if unknown
    init_int nrowpi;                          // number of rows in the simplified design matrix for pi
    init_matrix pidm(1,nrowpi,1,kpi);         // design matrix for pi
    init_vector pifix(1,nrowpi);              // pi fixed values
    init_ivector piindex(1,n*nS);             // pi indices
    init_int initknown;                       // if 1 then delta not used on release occasion
    init_int debug;                           // if 1 then write out parameters and likelihood at each iteration.

		
PARAMETER_SECTION
    init_vector phibeta(1,kphi);       // parameter vector for Phi
    init_vector pbeta(1,kp);           // parameter vector for p
    init_vector psibeta(1,kpsi);       // parameter vector for psi
    init_vector pibeta(1,kpi);         // parameter vector for pi
    init_vector dbeta(1,kd);           // parameter vector for delta
    
	objective_function_value g; 

PROCEDURE_SECTION
    int i,j,k,bindex,bindex2,k2;       // indices and counters
    int bindex3,bindex4;
    dvar_vector uniquephi(1,nrowphi);  // all unique phi values    
    dvar_vector phi(1,nrows);          // temp vector for Phis for an individual
    dvar_vector uniquep(1,nrowp);      // all unique p values    
    dvar_vector p(1,nrows);            // temp vector for ps for an individual
    dvar_vector uniquepsi(1,nrowpsi);  // temp vector for psis 
    dvariable psisum;                  // sum of psi for each state to normalize with
    dvar3_array psi(1,m-1,1,nS,1,nS);  // matrix for psis for each occasion 
    dvar_vector uniquepi(1,nrowpi);    // temp vector for pi 
    dvar_vector pi(1,nS);              // temp vector for pi for an individual 
    dvariable pisum;                   // sum of pi across states for an individual to normalize with
    dvar_vector uniquedelta(1,nrowd);  // temp vector for delta 
    dvar_matrix delta(1,m,1,npos);     // matrix of delta values by occasion for an individual 
    dvariable deltasum;                // sum of delta across possible observations within a state

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
          
    for(j=1;j<=nrowpi;j++)
    {
       if(pifix(j) < -0.5)
           uniquepi(j)=exp(pidm(j)*pibeta);           // compute exp of pidm*pibeta; these are components of pi
       else
           uniquepi(j)=pifix(j);                      // fixed exp pi value     
    }  
 
    for(j=1;j<=nrowd;j++)
    {
       if(dfix(j) < -0.5)
           uniquedelta(j)=exp(ddm(j)*dbeta);         // compute exp of ddm*dbeta; these are components of delta
       else
           uniquedelta(j)=dfix(j);                   // fixed exp delta value     
    }  

  

    for(i=1;i<=n;i++)                              // loop over capture histories - one per capture history
    {
        bindex=(i-1)*nrows;                        // initialize indices into index values for the ith history
        bindex2=(i-1)*nT;
        bindex3=(i-1)*nS;    
        bindex4=(i-1)*m*npos;
        pisum=0;                                   // compute initial state distribution
        for (k=1;k<=nS;k++)                        
          pisum=pisum+uniquepi(piindex(bindex3+k));
        for (k=1;k<=nS;k++)                        
          pi(k)=uniquepi(piindex(bindex3+k))/pisum;
                                                          
	    for (j=1;j<=m-1;j++)                       // loop over occasions for the individual computing phi,psi,p
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
	    for (j=1;j<=m;j++)                       // loop over occasions for delta
	    {                       
  	       for (k=1;k<=nS;k++)                        
           {   
	           deltasum=0;
	           for(k2=1;k2<=dcells;k2++)              
	               deltasum=deltasum+uniquedelta(dindex(bindex4+(k-1)*dcells+k2));
	           for(k2=1;k2<=dcells;k2++)              
	               delta(j,(k-1)*dcells+k2)=uniquedelta(dindex(bindex4+(k-1)*dcells+k2))/deltasum;
	       }
	       bindex4=bindex4+npos;      
        }        

        ll_i(i,phi,p,psi,pi,delta,ipos);                 // compute neg log likelihod and increment
     }
     if(debug>0)
     {
         cout << "\nphi " <<phibeta;
         cout << "\np " <<pbeta;
         cout << "\npsi " <<psibeta;
         cout << "\npi " <<pibeta;
         cout << "\ndelta " <<dbeta;
         cout << "\n-lnl = " << g << "\n";
     }

FUNCTION dvar3_array get_dmat(const int i, const dvar_vector& p, const dvar_matrix& delta, const imatrix ipos)
	dvar3_array dmat(1,m,1,nobs,1,nS+1);   // state dependent observation matrices for individual i
	dvar3_array pmat(1,m,1,nobs,1,nS+1);   // observation probability matrices for individual i
	int j,k,bindex,irow,icol;
	dmat.initialize();
	pmat.initialize();
    for(j=frst(i);j<=m;j++)
    {
      bindex=(j-2)*nS;
	  for(k=1;k<=npos;k++)
	  {
	    irow=ipos(k,1);
	    icol=ipos(k,2);
	    if(j==frst(i)) 
        {
           if(initknown>0)
              dmat(j,irow,icol)=1;  
           else
              dmat(j,irow,icol)=delta(j,k);  
           dmat(j,1,icol)=0;  
        } else
        {	    
	       dmat(j,irow,icol)=p(bindex+icol)*delta(j,k);  
           dmat(j,1,icol)=1-p(bindex+icol);
        }  
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

FUNCTION void ll_i(const int i, const dvar_vector& phi, const dvar_vector& p, const dvar3_array& psi, const dvar_vector& pi, const dvar_matrix& delta, const imatrix ipos)
//  Code implements an algorithm described on pg 47 in Zucchini and MacDonald (Z&M) for computing likelihood value for a 
//  Hidden Markov Model. The algorithm is recursive over the encounter occasions. The algorithm involves a simpler recursion 
//  shown on pages 38 and 45 but to avoid numerical underflows it is scaled by the sum of the state probabilities which makes
//  it appear a little more complex.
	dvar3_array gamma(1,m-1,1,nS+1,1,nS+1); // transition probability matrices for individual i
	dvar3_array dmat(1,m,1,nobs,1,nS+1);    // observation probability matrices for individual i
    dvariable u;                            // sum of state probabilities
	dvar_vector pS(1,nS+1);                 // update vector for prob of being in state j=1,nS + death       
	dvar_vector S(1,nS+1);                  // prob of being in state j=1,nS + death for each occasion
	dvar_vector v(1,nS+1);                  // temporary update vector
	dvariable Lglki;                        // log-likelihood accumulator
    int j;                                  // loop variable
     
//  compute transition matrices for each occasion
	gamma=get_gamma(i,phi,psi);
//  compute state dependent observation matrices for each occasion
	dmat=get_dmat(i,p,delta,ipos);
//  HMM algorithm
	pS.initialize();                                    // initialize values to 0
    Lglki.initialize();	
	S.initialize();                                     // initialize values to 0
	for(j=1;j<=nS;j++)
	  S(j)=pi(j);                                       // set state prob to pi for first observation
    for(j=frst(i);j<=m;j++)                             // loop over possible occasions from first(i)+1 to m
    {
        if(j==frst(i))
          pS=S;
        else
          pS=S*gamma(j-1);        	                    // previous scaled state probability*transition matrix
        v=elem_prod(pS,dmat(j,ch(i,j)));                // v is temp state vector alpha in Z&M
    	u=sum(v);                                       // sum across states
        S=v/u;                                          // update S;S is normalized alpha(j) (phi) in Z&M
	    Lglki+=log(u);    	                            // accumulate log-likelihood value
	}
	g-=freq(i)*Lglki;
