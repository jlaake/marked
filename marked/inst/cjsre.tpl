// Cormack-Jolly-Seber model with individual random effects for survival (Phi) and capture probability (p)
// Jeff Laake; 6 March 2013
DATA_SECTION 
    init_int n;                                   // number of capture histories
	init_int m;                                   // number of capture occasions
	init_imatrix ch(1,n,1,m);                     // capture history matrix
	init_ivector frst(1,n);                       // occasion first seen for each history
	init_ivector lst(1,n);                        // occasion last seen for each history
    init_vector frq(1,n);                         // frequency of each history
	init_ivector loc(1,n);                        // 0 or 1, 1 if lost on capture at last event
	init_matrix tint(1,n,1,m-1);                  // time interval between occasions for each history-interval
    init_int kphi;                                // number of columns in the design matrix for Phi - survival
	int nrows;                                    // number of entries in design matrix m-1 values for each of n histories
	!! nrows=n*(m-1);
    init_matrix phi_fixedDM(1,nrows,1,kphi+1);    // design matrix for Phi; final column can contain fixed values
    init_int phi_nre;                             // number of random effects for phi
	int phi_phase;                                // phase for re for phi; -1 if no random effects                                   
	!! phi_phase=2;                                   
	!! if(phi_nre==0)phi_phase=-1;
	int phi_nre_rows;
	!! phi_nre_rows=0;
	!! if(phi_nre>0)phi_nre_rows=phi_nre;
    init_int phi_krand;                           // number of columns in phi random effect DM
    init_matrix phi_randDM(1,nrows,1,phi_krand);  // phi random effect DM
    init_int kp;                                  // number of columns in the design matrix for p - capture probability
    init_matrix p_fixedDM(1,nrows,1,kp+1);        // design matrix for p; final column can contain fixed values

    init_int p_nre;                               // number of random effects for p
	int p_phase;                                  // phase for re for p; -1 if no random effects                                   
	!! p_phase=2;                                   
	!! if(p_nre==0)p_phase=-1;
	int p_nre_rows;
	!! p_nre_rows=0;
	!! if(p_nre>0)p_nre_rows=p_nre;
    init_int p_krand;                             // number of columns in phi random effect DM
    init_matrix p_randDM(1,nrows,1,p_krand);      // p random effect DM
 
    // Weights applied to the likelihood to obtain normalizing probablity
    // Note that (by mistake) weights run 2,4,6,....
    vector w(1,2*n)
	int ii;
    !! w=0.0;       
	!! for(ii=1;ii<=n;ii++)
	!!    w(2*ii)=frq(ii);
	
PARAMETER_SECTION
	init_vector phi_beta(1,kphi,1);                                      // parameter vector for Phi
	init_vector p_beta(1,kp,1);                                          // parameter vector for p
    init_vector phi_sigma(1,phi_krand,phi_phase);                        // log-sigma for random effects in phi;             
    init_vector p_sigma(1,p_krand,p_phase);                              // log-sigma for random effects in p;             
	objective_function_value f;                                          // objective function - negative log-likelihood
    !!set_multinomial_weights(w);                                        // weights to handle frequency
    random_effects_matrix phi_u(1,phi_nre_rows,1,phi_krand,phi_phase);   // random effect for scale
    random_effects_matrix p_u(1,p_nre_rows,1,p_krand,p_phase);           // random effect for scale

GLOBALS_SECTION
// Modification of files needed to make set_multinomial_weights work
   #include "minfil.cpp"	    // Laplace approx (too inaccurate here); f1b2fnl3.cpp in ADMB code base
   #include "df1b2gh.cpp"	    // Gauss-Hermite integration
   #include "xmodelm5.cpp"	    // Gauss-Hermite integration
   #include "df1b2ghmult.cpp"   // multiple G-H random effects
   
TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(100000); 
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(100000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  arrmblsize=5000000;

PROCEDURE_SECTION
    int i,i1,i2,j,L;	                   // miscellaneous ints
  //  int i,j;	                           // index over observations
	f=0;

	for(i=1;i<=n;i++)                      // loop over capture histories - one per animal
    {
        if(phi_nre_rows>0)
		{
		   if(p_nre_rows>0)
		      ll_i(i,phi_sigma,p_sigma,phi_u(i),p_u(i),phi_beta,p_beta);   // phi and p random effects
		   else
  		      nopll_i(i,phi_sigma,phi_u(i),phi_beta,p_beta);               // phi random effects only   	
        } else
		{
   		   if(p_nre_rows>0)
		      nophill_i(i,p_sigma,p_u(i),phi_beta,p_beta);                 // p random effects only
		   else
  		      norell_i(i,phi_beta,p_beta);                                 // no random effects 
		}
    }



SEPARABLE_FUNCTION void ll_i(const int i, const dvar_vector& phi_sigma,const dvar_vector& p_sigma,const dvar_vector& phi_u,const dvar_vector& p_u,const dvar_vector& phi_beta, const dvar_vector& p_beta )
    dvar_vector phi(1,m);              // temp vector for Phis for each occasion for a single history
    dvar_vector p(1,m-1);              // temp vector for Phis for each occasion for a single history
	dvar_vector phicumprod(1,m);       // cummulative survival probability across occasions
	dvar_vector cump(1,m);             // cummulative probability of being seen across occasions
	dvariable pch;                     // probability of capture history
    int i1,i2,j,L;	                   // miscellaneous ints
    dvariable mu;                      // link function value
    phi=0;                                                         // set all phi values to 0
    p=0;
   	phicumprod=1.0;                                                // set cummulative survival to 1
	cump=1.0;                                                      // set cummulative p to 1
	i1=(m-1)*(i-1);                                                // compute beginning index in design matrix
	for(j=frst(i)+1;j<=m;j++)                                      // loop over occasions from frst to m
	{
	   i2=i1+j-1;                                                  // increment index in design matrix
       /////// phi computation //////////
	   if(phi_fixedDM(i2,kphi+1)== -1)
	   {
	      mu=0;
	      for(L=1;L<=kphi;L++)
            mu+=phi_fixedDM(i2,L)*phi_beta(L);                     // fixed portion of mean 
          mu+=phi_randDM(i2)*elem_prod(phi_u,mfexp(phi_sigma));    // random portion
          phi(j-1)=1/(1+exp(-mu));                                 // compute phi for the interval
          if(tint(i,j-1)!=1)
             phi(j-1)=pow(phi(j-1),tint(i,j-1));                   // adjust phi for the time interval length
	   } else
	      phi(j-1)=phi_fixedDM(i2,kphi+1);
	   /////// p computation /////////
	   if(p_fixedDM(i2,kp+1)== -1)
	   {
   	      mu=0;
	      for(L=1;L<=kp;L++)
            mu+=p_fixedDM(i2,L)*p_beta(L);                         // fixed portion of mean 
          mu+=p_randDM(i2)*elem_prod(p_u,mfexp(p_sigma));          // random portion
          p(j-1)=1/(1+exp(-mu));                                     
       } else
	     p(j-1)=p_fixedDM(i2,kp+1);

	   phicumprod(j)=phicumprod(j-1)*phi(j-1);                     // compute cummulative survival
	   cump(j)=cump(j-1)*((1-p(j-1))*(1-ch(i,j))+p(j-1)*ch(i,j));  // compute cummulative capture probability
	}   
    pch=0.0;                                                       // initialize prob of the capture history to 0
    for(j=lst(i);j<=((1-loc(i))*m+loc(i)*lst(i));j++)              // loop over last occasion to m unless loss on capture
    {                                                              // to compute prob of the capture history 
       if(loc(i)==1)
          pch=pch+cump(j)*phicumprod(j);                           // probability of history given possible last occasion alive
       else       
          pch=pch+cump(j)*phicumprod(j)*(1-phi(j));                // probability of history given possible last occasion alive
    }   
    f-= log(pch+0.000000000000000000000001);                                // sum log-likelihood log(pr(ch))
	f-= sum(-0.5*log(2.0*M_PI) - 0.5*square(phi_u));
	f-= sum(-0.5*log(2.0*M_PI) - 0.5*square(p_u));

SEPARABLE_FUNCTION void nopll_i(const int i, const dvar_vector& phi_sigma,const dvar_vector& phi_u,const dvar_vector& phi_beta, const dvar_vector& p_beta )
    dvar_vector phi(1,m);              // temp vector for Phis for each occasion for a single history
    dvar_vector p(1,m-1);              // temp vector for Phis for each occasion for a single history
	dvar_vector phicumprod(1,m);       // cummulative survival probability across occasions
	dvar_vector cump(1,m);             // cummulative probability of being seen across occasions
	dvariable pch;                     // probability of capture history
    int i1,i2,j,L;	                   // miscellaneous ints
    dvariable mu;                      // link function value
    phi=0;                                                         // set all phi values to 0
    p=0;
   	phicumprod=1.0;                                                // set cummulative survival to 1
	cump=1.0;                                                      // set cummulative p to 1
	i1=(m-1)*(i-1);                                                // compute beginning index in design matrix
	for(j=frst(i)+1;j<=m;j++)                                      // loop over occasions from frst to m
	{
	   i2=i1+j-1;                                                  // increment index in design matrix
       /////// phi computation //////////
	   if(phi_fixedDM(i2,kphi+1)== -1)
	   {
	      mu=0;
	      for(L=1;L<=kphi;L++)
            mu+=phi_fixedDM(i2,L)*phi_beta(L);                     // fixed portion of mean 
          mu+=phi_randDM(i2)*elem_prod(phi_u,mfexp(phi_sigma));    // random portion
          phi(j-1)=1/(1+exp(-mu));                                 // compute phi for the interval
          if(tint(i,j-1)!=1)
             phi(j-1)=pow(phi(j-1),tint(i,j-1));                   // adjust phi for the time interval length
	   } else
	      phi(j-1)=phi_fixedDM(i2,kphi+1);
	   /////// p computation /////////
	   if(p_fixedDM(i2,kp+1)== -1)
	   {
   	      mu=0;
	      for(L=1;L<=kp;L++)
            mu+=p_fixedDM(i2,L)*p_beta(L);                       // fixed portion of mean 
          p(j-1)=1/(1+exp(-mu));                                     
       } else
	     p(j-1)=p_fixedDM(i2,kp+1);
	   
       phicumprod(j)=phicumprod(j-1)*phi(j-1);                     // compute cummulative survival
	   cump(j)=cump(j-1)*((1-p(j-1))*(1-ch(i,j))+p(j-1)*ch(i,j));  // compute cummulative capture probability
	}   
    pch=0.0;                                                       // initialize prob of the capture history to 0
    for(j=lst(i);j<=((1-loc(i))*m+loc(i)*lst(i));j++)              // loop over last occasion to m unless loss on capture
    {                                                              // to compute prob of the capture history 
       if(loc(i)==1)
          pch=pch+cump(j)*phicumprod(j);                           // probability of history given possible last occasion alive
       else       
          pch=pch+cump(j)*phicumprod(j)*(1-phi(j));                // probability of history given possible last occasion alive
    }   
    f-= log(pch+0.000000000000000000000001);                                // sum log-likelihood log(pr(ch))
	f-= sum(-0.5*log(2.0*M_PI) - 0.5*square(phi_u));

SEPARABLE_FUNCTION void nophill_i(const int i, const dvar_vector& p_sigma,const dvar_vector& p_u,const dvar_vector& phi_beta, const dvar_vector& p_beta )
    dvar_vector phi(1,m);              // temp vector for Phis for each occasion for a single history
    dvar_vector p(1,m-1);              // temp vector for Phis for each occasion for a single history
	dvar_vector phicumprod(1,m);       // cummulative survival probability across occasions
	dvar_vector cump(1,m);             // cummulative probability of being seen across occasions
	dvariable pch;                     // probability of capture history
    int i1,i2,j,L;	                   // miscellaneous ints
    dvariable mu;                      // link function value
    phi=0;                                                         // set all phi values to 0
    p=0;
   	phicumprod=1.0;                                                // set cummulative survival to 1
	cump=1.0;                                                      // set cummulative p to 1
	i1=(m-1)*(i-1);                                                // compute beginning index in design matrix
	for(j=frst(i)+1;j<=m;j++)                                      // loop over occasions from frst to m
	{
	   i2=i1+j-1;                                                  // increment index in design matrix
       /////// phi computation //////////
	   if(phi_fixedDM(i2,kphi+1)== -1)
	   {
	      mu=0;
	      for(L=1;L<=kphi;L++)
            mu+=phi_fixedDM(i2,L)*phi_beta(L);                     // fixed portion of mean 
          phi(j-1)=1/(1+exp(-mu));                                 // compute phi for the interval
          if(tint(i,j-1)!=1)
             phi(j-1)=pow(phi(j-1),tint(i,j-1));                   // adjust phi for the time interval length
	   } else
	     phi(j-1)=phi_fixedDM(i2,kphi+1);	   
	   /////// p computation /////////
	   if(p_fixedDM(i2,kp+1)== -1)
	   {
   	      mu=0;
	      for(L=1;L<=kp;L++)
            mu+=p_fixedDM(i2,L)*p_beta(L);                         // fixed portion of mean 
          mu+=p_randDM(i2)*elem_prod(p_u,mfexp(p_sigma));          // random portion
          p(j-1)=1/(1+exp(-mu));                                     
       } else
	     p(j-1)=p_fixedDM(i2,kp+1);

       phicumprod(j)=phicumprod(j-1)*phi(j-1);                     // compute cummulative survival
	   cump(j)=cump(j-1)*((1-p(j-1))*(1-ch(i,j))+p(j-1)*ch(i,j));  // compute cummulative capture probability
	}   
    pch=0.0;                                                       // initialize prob of the capture history to 0
    for(j=lst(i);j<=((1-loc(i))*m+loc(i)*lst(i));j++)              // loop over last occasion to m unless loss on capture
    {                                                              // to compute prob of the capture history 
       if(loc(i)==1)
          pch=pch+cump(j)*phicumprod(j);                           // probability of history given possible last occasion alive
       else       
          pch=pch+cump(j)*phicumprod(j)*(1-phi(j));                // probability of history given possible last occasion alive
    }   
    f-= log(pch+0.000000000000000000000001);                                // sum log-likelihood log(pr(ch))
	f-= sum(-0.5*log(2.0*M_PI) - 0.5*square(p_u));

SEPARABLE_FUNCTION void norell_i(const int i,const dvar_vector& phi_beta, const dvar_vector& p_beta )
    dvar_vector phi(1,m);              // temp vector for Phis for each occasion for a single history
    dvar_vector p(1,m-1);              // temp vector for Phis for each occasion for a single history
	dvar_vector phicumprod(1,m);       // cummulative survival probability across occasions
	dvar_vector cump(1,m);             // cummulative probability of being seen across occasions
	dvariable pch;                     // probability of capture history
    int i1,i2,j,L;	                   // miscellaneous ints
    dvariable mu;                      // link function value
    phi=0;                                                         // set all phi values to 0
    p=0;
   	phicumprod=1.0;                                                // set cummulative survival to 1
	cump=1.0;                                                      // set cummulative p to 1
	i1=(m-1)*(i-1);                                                // compute beginning index in design matrix
	for(j=frst(i)+1;j<=m;j++)                                      // loop over occasions from frst to m
	{
	   i2=i1+j-1;                                                  // increment index in design matrix
       /////// phi computation //////////
	   if(phi_fixedDM(i2,kphi+1)== -1)
	   {
	      mu=0;
	      for(L=1;L<=kphi;L++)
            mu+=phi_fixedDM(i2,L)*phi_beta(L);                     // fixed portion of mean 
          phi(j-1)=1/(1+exp(-mu));                                 // compute phi for the interval
          if(tint(i,j-1)!=1)
             phi(j-1)=pow(phi(j-1),tint(i,j-1));                   // adjust phi for the time interval length
	   } else
	     phi(j-1)=phi_fixedDM(i2,kphi+1);
	   /////// p computation /////////
	   if(p_fixedDM(i2,kp+1)== -1)
	   {
   	      mu=0;
	      for(L=1;L<=kp;L++)
            mu+=p_fixedDM(i2,L)*p_beta(L);                         // fixed portion of mean 
          p(j-1)=1/(1+exp(-mu));                                     
       } else
	     p(j-1)=p_fixedDM(i2,kp+1);
	   
	   phicumprod(j)=phicumprod(j-1)*phi(j-1);                     // compute cummulative survival
	   cump(j)=cump(j-1)*((1-p(j-1))*(1-ch(i,j))+p(j-1)*ch(i,j));  // compute cummulative capture probability
	}   
    pch=0.0;                                                       // initialize prob of the capture history to 0
    for(j=lst(i);j<=((1-loc(i))*m+loc(i)*lst(i));j++)              // loop over last occasion to m unless loss on capture
    {                                                              // to compute prob of the capture history 
       if(loc(i)==1)
          pch=pch+cump(j)*phicumprod(j);                           // probability of history given possible last occasion alive
       else       
          pch=pch+cump(j)*phicumprod(j)*(1-phi(j));                // probability of history given possible last occasion alive
    }   
       f-= log(pch+0.000000000000000000000001)*frq(i);                 // sum log-likelihood log(pr(ch))
		
REPORT_SECTION
    dvar_vector phi(1,m-1);              // temp vector for Phis for each occasion for a single history
    dvar_vector p(1,m-1);                // temp vector for Phis for each occasion for a single history
    int i,i1,i2,j,L;	                   // miscellaneous ints
    dvariable mu;                      // link function value
    for(i=1;i<=n;i++)
	{
    phi=0;                                                         // set all phi values to 0
	i1=(m-1)*(i-1);                                                // compute beginning index in design matrix
	for(j=frst(i)+1;j<=m;j++)                                      // loop over occasions from frst to m
	{
	   i2=i1+j-1;                                                  // increment index in design matrix
       /////// phi computation //////////
	   if(phi_fixedDM(i2,kphi+1)== -1)
	   {
	      mu=0;
	      for(L=1;L<=kphi;L++)
            mu+=phi_fixedDM(i2,L)*phi_beta(L);                     // fixed portion of mean 
		  if(phi_nre_rows>0)
		    mu+=phi_randDM(i2)*elem_prod(phi_u(i),mfexp(phi_sigma));
          phi(j-1)=1/(1+exp(-mu));                                 // compute phi for the interval
          if(tint(i,j-1)!=1)
             phi(j-1)=pow(phi(j-1),tint(i,j-1));                   // adjust phi for the time interval length
	   } else
	     phi(j-1)=phi_fixedDM(i2,kphi+1);
	}
	report << phi << endl;
	}
    for(i=1;i<=n;i++)
	{
	p=0;
	i1=(m-1)*(i-1);                                                // compute beginning index in design matrix
	for(j=frst(i)+1;j<=m;j++)                                      // loop over occasions from frst to m
	{
	   i2=i1+j-1;                                                  // increment index in design matrix
	   /////// p computation /////////
	   if(p_fixedDM(i2,kp+1)== -1)
	   {
   	      mu=0;
	      for(L=1;L<=kp;L++)
            mu+=p_fixedDM(i2,L)*p_beta(L);                         // fixed portion of mean 
		  if(p_nre_rows>0)
		    mu+=p_randDM(i2)*elem_prod(p_u(i),mfexp(p_sigma));
          p(j-1)=1/(1+exp(-mu));                                     
       } else
	     p(j-1)=p_fixedDM(i2,kp+1);
	   }
	report << p << endl;
	}


