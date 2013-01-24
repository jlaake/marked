// Cormack-Jolly-Seber model for fixed, random or mixed effect models for survival (Phi) and capture probability (p)
// Replaces admbcjs.tpl and admbcjsre.tpl
// Each capture history must represent a single animal
// Jeff Laake; 19 Jan 2013
DATA_SECTION 
    init_int n;                                       // number of capture histories
    init_int m;                                       // number of capture occasions
    init_imatrix ch(1,n,1,m);                         // capture history matrix
    init_ivector frst(1,n);                           // occasion first seen for each history
    init_ivector lst(1,n);                            // occasion last seen for each history
    init_ivector loc(1,n);                            // 0 or 1, 1 if lost on capture at last event
    init_matrix tint(1,n,1,m-1);                      // time interval between occasions for each history-interval
    int nrows;                                        // number of entries in design matrix m-1 values for each of n histories
    !! nrows=n*(m-1);
    init_int kphi;                                    // number of fixed effect columns in the design matrix for Phi - survival
    init_matrix phi_fixedDM(1,nrows,1,kphi);          // phi fixed effect DM
    init_int phi_krand;                               // number of columns in phi random effect DM
	int phi_phase;                                    // phase for re for phi; -1 if no random effects
	!! phi_phase=2;
	!! if(phi_krand==0)phi_phase=-1;
    init_int phi_nre;                                 // number of random effects for phi
    init_matrix phi_randDM(1,nrows,1,phi_krand);      // phi random effect DM
    init_imatrix phi_randIndex(1,nrows,1,phi_krand);  // phi random effect indices for DM
    init_int kp;                                      // number of fixed effect columns in the design matrix for p - capture probability
    init_matrix p_fixedDM(1,nrows,1,kp);              // p fixed effect DM
    init_int p_krand;                                  // number of columns in p random effect DM
	int p_phase;                                      // phase for re for p; -1 if no random effects
	!! p_phase=2;
	!! if(p_krand==0)p_phase=-1;
    init_int p_nre;                                   // number of random effects for p
    init_matrix p_randDM(1,nrows,1,p_krand);          // p random effect DM
    init_imatrix p_randIndex(1,nrows,1,p_krand);      // p random effect indices for DM
    init_int K;                                       // number of fixed Phi values 
    init_matrix PhiF(1,K,1,2);                        // Phi fixed matrix with index in first column and value in second column
    init_int L;                                       // number of fixed p values
    init_matrix pF(1,L,1,2);                          // p fixed matrix with index in first column and value in second column
        
PARAMETER_SECTION
    init_vector phi_beta(1,kphi,1);                   // fixed parameter vector for Phi
    init_vector p_beta(1,kp,1);                       // fixed parameter vector for p
    init_vector phi_sigma(1,phi_krand,phi_phase);     // log-sigma for random effects in phi;             
    init_vector p_sigma(1,p_krand,p_phase);           // log-sigma for random effects in p;             
    objective_function_value f;                       // objective function - negative log-likelihood
    random_effects_vector phi_u(1,phi_nre,phi_phase); // random effect vector for phi
    random_effects_vector p_u(1,p_nre,p_phase);       // random effect vector for p

TOP_OF_MAIN_SECTION
  arrmblsize=5000000;                                 // not sure how to set this as a function of problem size
  
PROCEDURE_SECTION
    int i;                                                         // index over observations
	f=0;
    if(phi_krand>0)	                                               // likelihood contribution for n(0,1) re for Phi
        for (i=1;i<=phi_nre;i++)
            n01_prior(phi_u(i));         

    if(p_krand>0)	                                               // likelihood contribution for n(0,1) re for p
        for (i=1;i<=p_nre;i++)
            n01_prior(p_u(i));         

    for(i=1;i<=n;i++)                                              // loop over capture histories - one per animal
    {
          ll_i(i,phi_sigma,p_sigma,phi_u,p_u,phi_beta,p_beta);
    }

SEPARABLE_FUNCTION void n01_prior(const prevariable&  u)           // taken from glmmadmb.tpl; uses PI 
    f -= -0.5*log(2.0*PI) - 0.5*square(u);

SEPARABLE_FUNCTION void ll_i(const int i, const dvar_vector& phi_sigma,const dvar_vector& p_sigma,const dvar_vector& phi_u,const dvar_vector& p_u,const dvar_vector& phi_beta, const dvar_vector& p_beta )
    dvar_vector phi(1,m);                                          // temp vector for Phis for each occasion for a single history
    dvar_vector p(1,m-1);                                          // temp vector for ps for each occasion for a single history
    dvar_vector phicumprod(1,m);                                   // cummulative survival probability across occasions
    dvar_vector cump(1,m);                                         // cummulative probability of being seen across occasions
    dvariable pch;                                                 // probability of capture history
    int i1,i2,j,L;                                                   // miscellaneous ints
    dvariable mu;                                                  // link function value
	    
    phi=0;                                                         // set all phi values to 0
    phicumprod=1.0;                                                // set cummulative survival to 1
    cump=1.0;                                                      // set cummulative p to 1
    i1=(m-1)*(i-1);                                                // compute beginning index in design matrix
    for(j=frst(i)+1;j<=m;j++)                                      // loop over occasions from frst to m
    {
       i2=i1+j;                                                    // increment index in design matrix
       // phi computation // 
       mu=phi_fixedDM(i)*phi_beta;                                 // fixed portion of mean   
       if(phi_krand > 0)                                           // random portion of mean if any
	      for(L=1;L<=phi_krand;L++)                                 
	         mu+=phi_randDM(i2-1,L)*phi_u(phi_randIndex(i2-1,L))*mfexp(phi_sigma(L));
       phi(j-1)=1/(1+exp(-mu));                                    // compute phi for the interval
       if(tint(i,j-1)!=1)
          phi(j-1)=pow(phi(j-1),tint(i,j-1));                      // adjust phi for the time interval length
       // p computation //
  	   mu=p_fixedDM(i)*p_beta;                                     // fixed portion of mean   
       if(p_krand > 0)                                             // random portion of mean if any
	      for(L=1;L<=p_krand;L++)                                 
	         mu+=p_randDM(i2-1,L)*p_u(p_randIndex(i2-1,L))*mfexp(p_sigma(L));
       p(j-1)=1/(1+exp(-mu));                                     
       phicumprod(j)=phicumprod(j-1)*phi(j-1);                     // compute cummulative survival
       cump(j)=cump(j-1)*((1-p(j-1))*(1-ch(i,j))+p(j-1)*ch(i,j));  // compute cummulative capture probability
    }   
    pch=0.0;                                                       // initialize prob of the capture history to 0
    for(j=lst(i);j<=((1-loc(i))*m+loc(i)*lst(i));j++)              // loop over last occasion to m unless loss on capture
    {                                                              // to compute prob of the capture history 
       i2=i1+j;                                                    // index occasion
       if(loc(i)==1)
          pch=pch+cump(j)*phicumprod(j);                           // probability of history given possible last occasion alive
       else       
          pch=pch+cump(j)*phicumprod(j)*(1-phi(j));                // probability of history given possible last occasion alive
    }   
    f-= log(pch+.0000000000000001);                                // sum log-likelihood log(pr(ch))
