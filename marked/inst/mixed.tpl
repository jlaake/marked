// Skeleton template for a mixed effects model structure; with logistic regression likelihood
// Jeff Laake; 18 Jan 2013
DATA_SECTION 
    init_int n;                           // number of rows in data
	init_vector y(1,n);                   // vector of responses
    init_int kfixed;                      // number of columns in design matrix(DM) for fixed effects
    init_matrix fixedDM(1,n,1,kfixed);    // fixed effect DM
	init_int rand;                        // 0 or 1 to set; phase either -1 (fixed only) or 2 (fixed+random)
	int phase;   
	!! phase=2;
	!! if(rand==0)phase=-1;
    init_int krand;                       // number of columns in random effect DM
    init_int nre;                         // number of random effects
    init_matrix randDM(1,n,1,krand);      // random effect DM
    init_imatrix randIndex(1,n,1,krand);  // random effect indices for DM
        
PARAMETER_SECTION
    init_vector fixedBeta(1,kfixed,1);    // parameter vector for fixed effects
    init_vector randBeta(1,krand,phase);  // parameter vector for log(sigma)
    objective_function_value g;           // objective function - negative log-likelihood
    random_effects_vector u(1,nre,phase); // random effects vector
 
PROCEDURE_SECTION
    int i;                                // index variable
	g=0;
    if(krand>0)	
        for (i=1;i<=nre;i++)
            n01_prior(u(i));              // u's are N(0,1) distributed

	for(i=1;i<=n;i++)                     // loop over rows in data
    {
          ll_i(i,fixedBeta,randBeta,u);   //compute negative log-likelihood for each row
    }

SEPARABLE_FUNCTION void n01_prior(const prevariable&  u)  // taken from glmmadmb.tpl; uses PI 
 g -= -0.5*log(2.0*PI) - 0.5*square(u);
 
SEPARABLE_FUNCTION void ll_i(const int i, const dvar_vector& fixedBeta,const dvar_vector& randBeta,const dvar_vector& u)
    int j;
    dvariable mu;                                                    // mean
    mu=fixedDM(i)*fixedBeta;                                         // fixed portion of mean   
    if(krand > 0)
	   for(j=1;j<=krand;j++)                                         // random portion of mean
	      mu+=randDM(i,j)*u(randIndex(i,j))*mfexp(randBeta(j));
	
	g-= y(i)*mu - log(1+exp(mu));

TOP_OF_MAIN_SECTION
