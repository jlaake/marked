// Cormack-Jolly-Seber model with individual random effects for survival (Phi) and capture probability (p)
// Each capture history must represent a single animal
// Jeff Laake; 15 Nov 2012
DATA_SECTION 
    init_int n;                        // number of capture histories
    init_int m;                        // number of capture occasions
    init_imatrix ch(1,n,1,m);          // capture history matrix
    init_ivector frst(1,n);            // occasion first seen for each history
    init_ivector lst(1,n);             // occasion last seen for each history
    init_ivector loc(1,n);             // 0 or 1, 1 if lost on capture at last event
    init_matrix tint(1,n,1,m-1);       // time interval between occasions for each history-interval
    init_int kphi;                     // number of columns in the design matrix for Phi - survival
    int nrows;                         // number of entries in design matrix m-1 values for each of n histories
    !! nrows=n*(m-1);
    init_matrix phidm(1,nrows,1,kphi); // design matrix for Phi
    init_int kp;                       // number of columns in the design matrix for p - capture probability
    init_matrix pdm(1,nrows,1,kp);     // design matrix for p
    init_int K;                        // number of fixed Phi values 
    init_matrix PhiF(1,K,1,2);         // Phi fixed matrix with index in first column and value in second column
    init_int L;                        // number of fixed p values
    init_matrix pF(1,L,1,2);           // p fixed matrix with index in first column and value in second column
        
PARAMETER_SECTION
    init_vector phibeta(1,kphi,1);       // parameter vector for Phi
    init_vector pbeta(1,kp,1);           // parameter vector for p
    init_bounded_number phisigeps(0.000001,6,2);   
                                       // sigma for random effect in phi;             
    init_bounded_number psigeps(0.000001,6,2);   
                                       // sigma for random effect in p;             
    objective_function_value f;        // objective function - negative log-likelihood
    random_effects_vector phiu(1,n,2);   // random effect for scale
    random_effects_vector pu(1,n,2);     // random effect for scale

TOP_OF_MAIN_SECTION
  arrmblsize=5000000;                  // not sure how to set this as a function of problem size
  
PROCEDURE_SECTION
    int i;                             // index over observations
    for(i=1;i<=n;i++)                  // loop over capture histories - one per animal
    {
          ll_i(i,phisigeps,psigeps,phiu(i),pu(i),phibeta,pbeta);
    }

SEPARABLE_FUNCTION void ll_i(const int i, const dvariable& phisigeps,const dvariable& psigeps,const dvariable& phiu,const dvariable& pu,const dvar_vector& phibeta, const dvar_vector& pbeta )
    dvar_vector phi(1,m);              // temp vector for Phis for each occasion for a single history
    dvar_vector p(1,m-1);              // temp vector for ps for each occasion for a single history
    dvar_vector phicumprod(1,m);       // cummulative survival probability across occasions
    dvar_vector cump(1,m);             // cummulative probability of being seen across occasions
    dvariable pch;                     // probability of capture history
    int i1,i2,j;                       // miscellaneous ints
    dvariable phieps=phiu*phisigeps;   // scaled re for phi
    dvariable peps=pu*psigeps;         // scaled re for p 
                                                                   //    
    phi=0;                                                         // set all phi values to 0
    phicumprod=1.0;                                                // set cummulative survival to 1
    cump=1.0;                                                      // set cummulative p to 1
    i1=(m-1)*(i-1);                                                // compute beginning index in design matrix
    for(j=frst(i)+1;j<=m;j++)                                      // loop over occasions from frst to m
    {
       i2=i1+j;                                                    // increment index in design matrix
       phi(j-1)=1/(1+exp(-phidm(i2-1)*phibeta-phieps));            // compute phi for the interval
       if(tint(i,j-1)!=1)
          phi(j-1)=pow(phi(j-1),tint(i,j-1));                      // adjust phi for the time interval length
       p(j-1)=1/(1+exp(-pdm(i2-1)*pbeta-peps));                    // compute p for the occasion
       phicumprod(j)=phicumprod(j-1)*phi(j-1);                     // compute cummulative survival
       cump(j)=cump(j-1)*((1-p(j-1))*(1-ch(i,j))+p(j-1)*ch(i,j));  // compute cummulative capture probability
    }   
    pch=0.0;                                                       // initialize prob of the capture history to 0
    for(j=lst(i);j<=((1-loc(i))*m+loc(i)*lst(i));j++)              // loop over last occasion to m unless loss on capture
    {                                                              //  to compute prob of the capture history 
       i2=i1+j;                                                    // index occasion
       if(loc(i)==1)
          pch=pch+cump(j)*phicumprod(j);                           // probability of history given possible last occasion alive
       else       
          pch=pch+cump(j)*phicumprod(j)*(1-phi(j));                // probability of history given possible last occasion alive
    }   
    f-= log(pch+.0000000000000001);                                // sum log-likelihood log(pr(ch))
    f -= -0.5*square(pu)-0.9189385332046727;                       // normal random effect distr; constant is log(sqrt(2*pi))
    f -= -0.5*square(phiu)-0.9189385332046727;                     // normal random effect distr; constant is log(sqrt(2*pi))

