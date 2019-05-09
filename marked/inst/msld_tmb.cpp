// TMB Version: Mixed-effect Multi-State Cormack-Jolly-Seber + recovery model with unobservable states
// Jeff Laake; 18 Mar 2019

#include <TMB.hpp>                              // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n);                            // number of capture histories
  DATA_INTEGER(m);                            // number of capture occasions
  DATA_INTEGER(nS);                           // number of states excluding death states (newly dead and dead)
  DATA_IMATRIX(ch);                           // capture history matrix; uses numeric values for states
  DATA_IVECTOR(frst);                         // occasion first seen for each history
  DATA_VECTOR(freq);                          // frequency of each capture history
  DATA_MATRIX(tint);                          // time interval between occasions for each history-interval
  
  DATA_INTEGER(nrowphi);                      // number of rows in the simplified design matrix for Phi - survival
                                              // last column in each design matrix is either -1 (estimated) or a fixed value
  DATA_MATRIX(phidm);                         // design matrix for Phi
  DATA_VECTOR(phifix);                        // phi fixed values
  DATA_IVECTOR(phiindex);                     // phi indices
  DATA_INTEGER(phi_nre);                      // number of random effects for phi
  DATA_INTEGER(phi_krand);                    // number of columns in phi random effect DM
  DATA_MATRIX(phi_randDM);                    // phi random effect DM
  DATA_IMATRIX(phi_randIndex);                // phi random effect indices for DM
  DATA_IVECTOR(phi_counts);                   // count of phi random effect indices by id
  DATA_IMATRIX(phi_idIndex);                  // phi random effect indices by id

  DATA_INTEGER(nrowr);                        // number of rows in the simplified design matrix for r - recovery
                                              // last column in each design matrix is either -1 (estimated) or a fixed value
  DATA_MATRIX(rdm);                           // design matrix for r
  DATA_VECTOR(rfix);                          // r fixed values
  DATA_IVECTOR(rindex);                       // r indices
  DATA_INTEGER(r_nre);                        // number of random effects for r
  DATA_INTEGER(r_krand);                      // number of columns in r random effect DM
  DATA_MATRIX(r_randDM);                      // r random effect DM
  DATA_IMATRIX(r_randIndex);                  // r random effect indices for DM
  DATA_IVECTOR(r_counts);                     // count of r random effect indices by id
  DATA_IMATRIX(r_idIndex);                    // r random effect indices by id
  
  DATA_INTEGER(nrowp);                        // number of rows in the simplified design matrix for p
  DATA_MATRIX(pdm);                           // design matrix for p
  DATA_VECTOR(pfix);                          // p fixed values
  DATA_IVECTOR(pindex);                       // p indices
  DATA_INTEGER(p_nre);                        // number of random effects for p
  DATA_INTEGER(p_krand);                      // number of columns in p random effect DM
  DATA_MATRIX(p_randDM);                      // p random effect DM
  DATA_IMATRIX(p_randIndex);                  // p random effect indices for DM; index into p_u
  DATA_IVECTOR(p_counts);                     // count of p random effect indices by id
  DATA_IMATRIX(p_idIndex);                    // p random effect indices by id; index into u_phi to construct phi_u
  
  DATA_INTEGER(nrowpsi);                      // number of rows in the simplified design matrix for psi
  DATA_MATRIX(psidm);                         // design matrix for psi
  DATA_VECTOR(psifix);                        // psi fixed values
  DATA_IVECTOR(psiindex);                     // psi indices
  DATA_INTEGER(psi_nre);                      // number of random effects for psi
  DATA_INTEGER(psi_krand);                    // number of columns in psi random effect DM
  DATA_MATRIX(psi_randDM);                    // psi random effect DM
  DATA_IMATRIX(psi_randIndex);                // psi random effect indices for DM; index into phi_u
  DATA_IVECTOR(psi_counts);                   // count of psi random effect indices by id
  DATA_IMATRIX(psi_idIndex);                  // psi random effect indices by id; index into u_phi to construct phi_u
  
  DATA_INTEGER(getreals);                     // if 1, report reals and std errors
  
  PARAMETER_VECTOR(phibeta);                  // parameter vector for Phi
  PARAMETER_VECTOR(rbeta);                    // parameter vector for r
  PARAMETER_VECTOR(pbeta);                    // parameter vector for p
  PARAMETER_VECTOR(psibeta);                  // parameter vector for Psi
  PARAMETER_VECTOR(log_sigma_phi);
  PARAMETER_VECTOR(log_sigma_r);
  PARAMETER_VECTOR(log_sigma_p);
  PARAMETER_VECTOR(log_sigma_psi);
  PARAMETER_VECTOR(u_phi);
  PARAMETER_VECTOR(u_r);
  PARAMETER_VECTOR(u_p);
  PARAMETER_VECTOR(u_psi);
  
  Type g=0; 
  
  int nrows;                           // number of entries in design matrix m-1 values for each of nS states
  nrows=nS*(m-1);                
  int nT;                              // number of transitions excluding death
  nT=nS*nS*(m-1);                       
  
  int i,j,k,bindex,bindex2,k2,idx,i2;  // indices and counters
  int L; 
  vector<Type> uniquephi(nrowphi);     // all unique phi values    
  vector<Type> phi(nrows);             // temp vector for Phis for an individual
  vector<Type> uniquer(nrowr);         // all unique r values    
  vector<Type> r(nrows);               // temp vector for rs for an individual
  vector<Type> uniquep(nrowp);         // all unique p values    
  vector<Type> p(nrows);               // temp vector for ps for an individual
  vector<Type> uniquepsi(nrowpsi);     // temp vector for psis 
  Type psisum;                         // sum of psi for each state to normalize with
  
  array<Type> psi(m-1,nS,nS);             // matrix for psis for each occasion 
  array<Type> gamma(m-1,2*nS+1,2*nS+1);   // transition probability matrices for individual i
  array<Type> dmat(m-1,nS+2,2*nS+1);      // observation probability matrices for individual i
  Type u;                                 // sum of state probabilities
  vector<Type> pS(2*nS+1);              // update vector for prob of being in state j=1,2*nS + 1       
  vector<Type> S(2*nS+1);               // prob of being in state j=1,2*nS + 1 for each occasion
  vector<Type> v(2*nS+1);               // temporary update vector
  vector<Type> vec;                   // temporary vector
  Type Lglki=0;                       // log-likelihood accumulator
  Type mu;
  int nphicounts=n;                   // number of counts for phi random effects by id
  if(phi_nre==0)nphicounts=0;
  int nrcounts=n;                      // number of counts for r random effects by id
  if(r_nre==0)nrcounts=0;
  int npcounts=n;                     // number of counts for p random effects by id
  if(p_nre==0)npcounts=0;
  int npsicounts=n;                   // number of counts for psi random effects by id
  if(psi_nre==0)npsicounts=0;
  
  if(phi_krand>0)	                                     // likelihood contribution for n(0,1) re for phi
    for (int i=0;i<=phi_nre-1;i++)	   	     
      g-= dnorm(u_phi(i),Type(0),Type(1),true);
  
  if(r_krand>0)	                                        // likelihood contribution for n(0,1) re for r
    for (int i=0;i<=r_nre-1;i++)	   	     
      g-= dnorm(u_r(i),Type(0),Type(1),true);
  
  if(p_krand>0)	                                        // likelihood contribution for n(0,1) re for p
    for (int i=0;i<=p_nre-1;i++)
      g-= dnorm(u_p(i),Type(0),Type(1),true);
  
  if(psi_krand>0)	                                        // likelihood contribution for n(0,1) re for psi
    for (int i=0;i<=psi_nre-1;i++)
      g-= dnorm(u_psi(i),Type(0),Type(1),true);
  
  uniquephi=phidm*phibeta;                              // compute unique parameter sets on link scale
  uniquer=rdm*rbeta;    
  uniquep=pdm*pbeta;
  uniquepsi=psidm*psibeta;
  
  for(i=1;i<=n;i++)                                     // loop over capture histories - one per capture history
  {
    vector<Type> p_u(p_idIndex.cols());        // define random effects vector for p, Phi,r and psi used                  
    vector<Type> phi_u(phi_idIndex.cols());    // just for this capture history copied from full vectors
    vector<Type> r_u(r_idIndex.cols());
    vector<Type> psi_u(psi_idIndex.cols());
    p_u.setZero();	   
    phi_u.setZero();	   
    r_u.setZero();	   
    if(nphicounts >0)                          // if any random effects for phi, copy values from u_phi to phi_u
    {
      if(phi_counts(i-1)==0)
        phi_u(0)=0;
      else
        for(j=0;j<=phi_counts(i-1)-1;j++)
          phi_u(j)=u_phi(phi_idIndex(i-1,j)-1);
    } 
    
    if(nrcounts >0)                          // if any random effects for r, copy values from u_r to r_u
    {
      if(r_counts(i-1)==0)
        r_u(0)=0;
      else
        for(j=0;j<=r_counts(i-1)-1;j++)
          r_u(j)=u_r(r_idIndex(i-1,j)-1);
    } 
    
    if(npcounts >0)                           // if any random effects for p, copy values from u_p to p_u
    {
      if(p_counts(i-1)==0)
        p_u(0)=0;
      else
        for(j=0;j<=p_counts(i-1)-1;j++)
          p_u(j)=u_p(p_idIndex(i-1,j)-1);
    } 
    
    if(npsicounts >0)                           // if any random effects for psi, copy values from u_psi to psi_u
    {
      if(psi_counts(i-1)==0)
        psi_u(0)=0;
      else
        for(j=0;j<=psi_counts(i-1)-1;j++)
          psi_u(j)=u_psi(psi_idIndex(i-1,j)-1);
    } 
    //  compute phi and p values for the individual   
    bindex=(i-1)*nrows;                               // initialize indices into index values for the ith history
    bindex2=(i-1)*nT;    
    for (j=1;j<=m-1;j++)
    {                       
      for (k=1;k<=nS;k++)                        
      {   
        i2=bindex+(j-1)*nS+k;
        idx=pindex(i2-1)-1;
        if(pfix(idx)< -0.5)
        {
          mu=0;
          if(npcounts>0)
            if(p_counts(i-1) > 0)	                        // random portion of mean if any
            {
              for(L=1;L<=p_krand;L++)
                if(p_randIndex(i2-1,L-1)>0)
                  mu+=p_randDM(i2-1,L-1)*p_u(p_randIndex(i2-1,L-1)-1)*exp(log_sigma_p(L-1));
            }	           
            p((j-1)*nS+k-1)=1/(1+exp(-(uniquep(idx)+mu)));
        }
        else
          p((j-1)*nS+k-1)=pfix(idx);
        
        idx=phiindex(i2-1)-1;
        if(phifix(idx)< -0.5)
        {
          mu=0;
          if(nphicounts>0)
            if(phi_counts(i-1) > 0)	                        // random portion of mean if any
            {
              for(L=1;L<=phi_krand;L++)
                if(phi_randIndex(i2-1,L-1)>0)
                  mu+=phi_randDM(i2-1,L-1)*phi_u(phi_randIndex(i2-1,L-1)-1)*exp(log_sigma_phi(L-1));
            }	
            phi((j-1)*nS+k-1)=pow(1/(1+exp(-(uniquephi(idx)+mu))),tint(i-1,j-1)); 
        }
        else
          phi((j-1)*nS+k-1)=phifix(idx);     

        idx=rindex(i2-1)-1;
        if(rfix(idx)< -0.5)
        {
          mu=0;
          if(nrcounts>0)
            if(r_counts(i-1) > 0)	                        // random portion of mean if any
            {
              for(L=1;L<=r_krand;L++)
                if(r_randIndex(i2-1,L-1)>0)
                  mu+=r_randDM(i2-1,L-1)*r_u(r_randIndex(i2-1,L-1)-1)*exp(log_sigma_r(L-1));
            }	
            r((j-1)*nS+k-1)=(1/(1+exp(-(uniquer(idx)+mu)))); 
        }
        else
          r((j-1)*nS+k-1)=rfix(idx);     
        
        
        psisum=0;
        for(k2=1;k2<=nS;k2++)
        {
          i2=bindex2+(k-1)*nS+k2;
          idx=psiindex(i2-1)-1;              
          if(psifix(idx)< -0.5)
          {
            mu=0;
            if(npsicounts>0)
              if(psi_counts(i-1) > 0)	                 
              {
                for(L=1;L<=psi_krand;L++)
                  if(psi_randIndex(i2-1,L-1)>0)
                    mu+=psi_randDM(i2-1,L-1)*psi_u(psi_randIndex(i2-1,L-1)-1)*exp(log_sigma_psi(L-1));
              }       
              psi(j-1,k-1,k2-1)=exp(uniquepsi(idx)+mu);
          }	
          else
            psi(j-1,k-1,k2-1)=psifix(idx);
          psisum+=psi(j-1,k-1,k2-1);
        }
        for(k2=1;k2<=nS;k2++)
          psi(j-1,k-1,k2-1)=psi(j-1,k-1,k2-1)/psisum;
      }
      bindex2=bindex2+nS*nS;         
    }
    if(getreals>0)                                               // if requested report phi, r, p and psi values
    {
      ADREPORT(phi);
      ADREPORT(r);
      ADREPORT(p);
      ADREPORT(psi);
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
        gamma(j-1,k-1,nS+k-1)=1-phi(bindex-1);                    // add newly dead state value for each state
        gamma(j-1,nS+k-1,2*nS)=1;                                 // newly dead to permanently dead transition
        bindex++;
      }
      gamma(j-1,2*nS,2*nS)=1;                                    // permanently dead is an absorbing state
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
        dmat(j-1,0,nS+k-1)=1-r(bindex-1);  //did not recover newly dead
        dmat(j-1,nS+1,nS+k-1)=r(bindex-1); // recovered newly dead
        bindex++;
      }
      dmat(j-1,0,2*nS)=1;                   //cannot observe permanently dead
    }
    //  HMM algorithm
    pS.setZero();                                      // initialize values to 0
    Lglki=0;	
    S.setZero();                                     
    S(ch(i-1,frst(i-1)-1)-1)=1;                        // set state prob to 1 for observed state at first observation
    for(j=frst(i-1)+1;j<=m;j++)                       // loop over possible occasions from first(i)+1 to m
    {
      for(k=1;k<=2*nS+1;k++)
      {
        pS(k-1)=0;
        for(k2=1;k2<=2*nS+1;k2++)
          pS(k-1)+= S(k2-1)*gamma(j-2,k2-1,k-1);
      }
      for(k=1;k<=2*nS+1;k++)
      {
        v(k-1)=pS(k-1)*dmat(j-2,ch(i-1,j-1),k-1); // v is temp state vector alpha in Z&M
      }
      if(v.sum()==0) {
        Rcout << "\n Check Psi or p values set to 0";
        Rcout << "\n i = " << i << " ch = " << ch(i-1);
      }
      u=v.sum();                                       // sum across states
      S=v/u;                                          // update S;S is normalized alpha(j) (phi) in Z&M
      Lglki+=log(u);    	                            // accumulate log-likelihood value
      
    }
    g-=freq(i-1)*Lglki;
    
  }
  return g;
}
