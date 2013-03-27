// Cormack-Jolly-Seber model for fixed, random or mixed effect models for survival (Phi) and capture probability (p)
// Same as cjs.tpl except that it uses reml estimation. Each capture history must represent a single animal.
// Jeff Laake; 6 Mar 2013
DATA_SECTION 
    init_int n;                                           // number of capture histories
    init_int m;                                           // number of capture occasions
    init_imatrix ch(1,n,1,m);                             // capture history matrix
    init_ivector frst(1,n);                               // occasion first seen for each history
    init_ivector lst(1,n);                                // occasion last seen for each history
    init_ivector loc(1,n);                                // 0 or 1, 1 if lost on capture at last event
    init_matrix tint(1,n,1,m-1);                          // time interval between occasions for each history-interval
    int nrows;                                            // number of entries in design matrix m-1 values for each of n histories
    !! nrows=n*(m-1);

    init_int kphi;                                        // number of fixed effect columns in the design matrix for Phi - survival
    init_matrix phi_fixedDM(1,nrows,1,kphi+1);            // phi fixed effect DM
    init_int phi_nre;                                     // number of random effects for phi
	int phi_phase;                                        // phase for re for phi; -1 if no random effects                                   
	!! phi_phase=2;                                   
	!! if(phi_nre==0)phi_phase=-1;
    init_int phi_krand;                                   // number of columns in phi random effect DM
    init_matrix phi_randDM(1,nrows,1,phi_krand);          // phi random effect DM
    init_imatrix phi_randIndex(1,nrows,1,phi_krand);      // phi random effect indices for DM
    int nphicounts;                                       // number of counts for phi random effects by id
	!! nphicounts=n;
    !! if(phi_nre==0)nphicounts=0;
	int phi_nre_rows;
	!! phi_nre_rows=0;
	!! if(phi_nre>0)phi_nre_rows=phi_nre-1;
    init_ivector phi_counts(1,nphicounts);                // count of phi random effect indices by id
    init_imatrix phi_idIndex(1,nphicounts,1,phi_counts);  // phi random effect indices by id

    init_int kp;                                          // number of fixed effect columns in the design matrix for p - capture probability
    init_matrix p_fixedDM(1,nrows,1,kp+1);                // p fixed effect DM
    init_int p_nre;                                       // number of random effects for p
	int p_phase;                                          // phase for re for p; -1 if no random effects
	!! p_phase=2;                                     
	!! if(p_nre==0)p_phase=-1;
	int p_nre_rows;
	!! p_nre_rows=0;
	!! if(p_nre>0)p_nre_rows=p_nre-1;
    init_int p_krand;                                     // number of columns in p random effect DM
    init_matrix p_randDM(1,nrows,1,p_krand);              // p random effect DM
    init_imatrix p_randIndex(1,nrows,1,p_krand);          // p random effect indices for DM
    int npcounts;                                         // number of counts for p random effects by id
	!! npcounts=n;
    !! if(p_nre==0)npcounts=0;
    init_ivector p_counts(1,npcounts);                    // count of p random effect indices by id
    init_imatrix p_idIndex(1,npcounts,1,p_counts);        // p random effect indices by id

    init_int K;                                           // number of fixed Phi values 
    init_matrix PhiF(1,K,1,2);                            // Phi fixed matrix with index in first column and value in second column
    init_int L;                                           // number of fixed p values
    init_matrix pF(1,L,1,2);                              // p fixed matrix with index in first column and value in second column
        
PARAMETER_SECTION
    init_vector phi_sigma(1,phi_krand);                   // log-sigma for random effects in phi;             
    init_vector p_sigma(1,p_krand);                       // log-sigma for random effects in p;             
    objective_function_value f;                           // objective function - negative log-likelihood
    random_effects_vector phi_beta(1,kphi);               // re parameter vector for Phi
    random_effects_vector p_beta(1,kp);                   // re parameter vector for p
    random_effects_vector phi_u(0,phi_nre_rows);          // random effect vector for phi
    random_effects_vector p_u(0,p_nre_rows);              // random effect vector for p

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(10000); 
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(10000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(10000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  arrmblsize=5000000;

PROCEDURE_SECTION
    int i;                                                         // index over observations
	f=0;
    if(phi_krand>0)	                                               // likelihood contribution for n(0,1) re for Phi
        for (i=0;i<=phi_nre-1;i++)
            n01_prior(phi_u(i));         
	
    if(p_krand>0)	                                               // likelihood contribution for n(0,1) re for p
        for (i=0;i<=p_nre-1;i++)
            n01_prior(p_u(i));         

	for(i=1;i<=n;i++)                                              // loop over capture histories - one per animal
    {
	      int phi_upper,p_upper;                                   // begin block of code to handle 0 indices
		  if(nphicounts >0)
	      {
		     if(phi_counts(i)==0)
		        phi_upper=1;
		     else 
		        phi_upper=phi_counts(i);
		  } else
		      phi_upper=1;
		  if(npcounts >0)
	      {
    	      if(p_counts(i)==0)
	    	     p_upper=1;
		      else 
		         p_upper=p_counts(i);
		  } else
		      p_upper=1;
          ivector phi_indices(1,phi_upper);                          
          ivector p_indices(1,p_upper);                              
		  if(nphicounts >0)
	      {
   		     if(phi_counts(i)==0)
		        phi_indices(1)=0;
		     else
		        phi_indices=phi_idIndex(i)-1;
		  } else
         	  phi_indices(1)=0;	  
		  if(npcounts >0)
	      {
   		      if(p_counts(i)==0)
		         p_indices(1)=0;
		      else
		         p_indices=p_idIndex(i)-1;				 
		  } else
		      p_indices(1)=0;                                      // end block of code to handle 0 indices
			  
 	      if(phi_krand==0 & p_krand>0)                             
		     nophill_i(i,p_sigma,p_u(p_indices),phi_beta,p_beta);                              // p random effects; phi fixed
		  else
		     if(p_krand==0 & phi_krand>0)                         
  		         nopll_i(i,phi_sigma,phi_u(phi_indices),phi_beta,p_beta);                      // phi random effects; p fixed
			 else 			 
		  	     ll_i(i,phi_sigma,p_sigma,phi_u(phi_indices),p_u(p_indices),phi_beta,p_beta);  // phi and p random effects		 
    }
	

SEPARABLE_FUNCTION void n01_prior(const prevariable&  u)           // taken from glmmadmb.tpl; uses PI 
    f -= -0.5*log(2.0*M_PI) - 0.5*square(u);

SEPARABLE_FUNCTION void ll_i(const int i, const dvar_vector& phi_sigma,const dvar_vector& p_sigma,const dvar_vector& phi_u,const dvar_vector& p_u,const dvar_vector& phi_beta, const dvar_vector& p_beta )
    dvar_vector phi(1,m);                                          // temp vector for Phis for each occasion for a single history
    dvar_vector p(1,m-1);                                          // temp vector for ps for each occasion for a single history
    dvar_vector phicumprod(1,m);                                   // cummulative survival probability across occasions
    dvar_vector cump(1,m);                                         // cummulative probability of being seen across occasions
    dvariable pch;                                                 // probability of capture history
    int i1,i2,j,L;                                                 // miscellaneous ints
    dvariable mu;                                                  // link function value
    phi=0;                                                         // set all phi and p values to 0
	p=0;
    phicumprod=1.0;                                                // set cummulative survival to 1
    cump=1.0;                                                      // set cummulative p to 1
    i1=(m-1)*(i-1);  	                                           // compute beginning index in design matrix
    for(j=frst(i)+1;j<=m;j++)                                      // loop over occasions from frst to m
    {
       i2=i1+j-1;                                                  // increment index in design matrix
       /////// phi computation //////////

	   if(phi_fixedDM(i2,kphi+1)== -1)
	   {
	      mu=0;
	      for(L=1;L<=kphi;L++)
            mu+=phi_fixedDM(i2,L)*phi_beta(L);                     // fixed portion of mean 
       if(phi_phase>0 & current_phase()==phi_phase)
	   {
           if(phi_counts(i) > 0)	                               // random portion of mean if any
	       {
	          for(L=1;L<=phi_krand;L++)
		         if(phi_randIndex(i2,L)>0)
	                mu+=phi_randDM(i2,L)*phi_u(phi_randIndex(i2,L))*mfexp(phi_sigma(L));
	       }
	   }
       phi(j-1)=1/(1+exp(-mu));                                    // compute phi for the interval
       if(tint(i,j-1)!=1)
          phi(j-1)=pow(phi(j-1),tint(i,j-1));                      // adjust phi for the time interval length
	   } else
	      phi(j-1)=phi_fixedDM(i2,kphi+1);
	   /////// p computation /////////
	   if(p_fixedDM(i2,kp+1)== -1)
	   {
	      mu=0;
	      for(L=1;L<=kp;L++)
            mu+=p_fixedDM(i2,L)*p_beta(L);                         // fixed portion of mean 
	   if(p_phase>0 & current_phase()==p_phase)
	   {
           if(p_counts(i) > 0)                                     // random portion of mean if any
	       {
	         for(L=1;L<=p_krand;L++)                                 
		         if(p_randIndex(i2,L)>0)
	                mu+=p_randDM(i2,L)*p_u(p_randIndex(i2,L))*mfexp(p_sigma(L));
	       }
	   }
          p(j-1)=1/(1+exp(-mu));                                    // compute p for the occasion
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

SEPARABLE_FUNCTION void nophill_i(const int i, const dvar_vector& p_sigma,const dvar_vector& p_u,const dvar_vector& phi_beta, const dvar_vector& p_beta )
    dvar_vector phi(1,m);                                          // temp vector for Phis for each occasion for a single history
    dvar_vector p(1,m-1);                                          // temp vector for ps for each occasion for a single history
    dvar_vector phicumprod(1,m);                                   // cummulative survival probability across occasions
    dvar_vector cump(1,m);                                         // cummulative probability of being seen across occasions
    dvariable pch;                                                 // probability of capture history
    int i1,i2,j,L;                                                 // miscellaneous ints
    dvariable mu;                                                  // link function value
    phi=0;                                                         // set all phi and p values to 0
	p=0;
    phicumprod=1.0;                                                // set cummulative survival to 1
    cump=1.0;                                                      // set cummulative p to 1
    i1=(m-1)*(i-1);  	                                           // compute beginning index in design matrix
    for(j=frst(i)+1;j<=m;j++)                                      // loop over occasions from frst to m
    {
       i2=i1+j-1;                                                  // increment index in design matrix
       /////// phi computation //////////
	   if(phi_fixedDM(i2,kphi+1)== -1)
	   {
	      mu=0;
	      for(L=1;L<=kphi;L++)
            mu+=phi_fixedDM(i2,L)*phi_beta(L);                     // fixed portion of mean 
       phi(j-1)=1/(1+exp(-mu));                                    // compute phi for the interval
       if(tint(i,j-1)!=1)
          phi(j-1)=pow(phi(j-1),tint(i,j-1));                      // adjust phi for the time interval length
	   } else
	      phi(j-1)=phi_fixedDM(i2,kphi+1);
	   /////// p computation /////////
	   if(p_fixedDM(i2,kp+1)== -1)
	   {
	      mu=0;
	      for(L=1;L<=kp;L++)
            mu+=p_fixedDM(i2,L)*p_beta(L);                         // fixed portion of mean 
	   if(p_phase>0 & current_phase()==p_phase)
	   {
           if(p_counts(i) > 0)                                     // random portion of mean if any
	       {
	         for(L=1;L<=p_krand;L++)                                 
		         if(p_randIndex(i2,L)>0)
	                mu+=p_randDM(i2,L)*p_u(p_randIndex(i2,L))*mfexp(p_sigma(L));
	       }
	   }
          p(j-1)=1/(1+exp(-mu));                                    // compute p for the occasion
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
	
SEPARABLE_FUNCTION void nopll_i(const int i, const dvar_vector& phi_sigma,const dvar_vector& phi_u,const dvar_vector& phi_beta, const dvar_vector& p_beta )
    dvar_vector phi(1,m);                                          // temp vector for Phis for each occasion for a single history
    dvar_vector p(1,m-1);                                          // temp vector for ps for each occasion for a single history
    dvar_vector phicumprod(1,m);                                   // cummulative survival probability across occasions
    dvar_vector cump(1,m);                                         // cummulative probability of being seen across occasions
    dvariable pch;                                                 // probability of capture history
    int i1,i2,j,L;                                                 // miscellaneous ints
    dvariable mu;                                                  // link function value
    phi=0;                                                         // set all phi and p values to 0
	p=0;
    phicumprod=1.0;                                                // set cummulative survival to 1
    cump=1.0;                                                      // set cummulative p to 1
    i1=(m-1)*(i-1);  	                                           // compute beginning index in design matrix
    for(j=frst(i)+1;j<=m;j++)                                      // loop over occasions from frst to m
    {
       i2=i1+j-1;                                                  // increment index in design matrix
       /////// phi computation //////////
	   if(phi_fixedDM(i2,kphi+1)== -1)
	   {
	      mu=0;
	      for(L=1;L<=kphi;L++)
            mu+=phi_fixedDM(i2,L)*phi_beta(L);                     // fixed portion of mean 
          if(phi_phase>0 & current_phase()==phi_phase)
	   {
           if(phi_counts(i) > 0)	                               // random portion of mean if any
	       {
	          for(L=1;L<=phi_krand;L++)
		         if(phi_randIndex(i2,L)>0)
	                mu+=phi_randDM(i2,L)*phi_u(phi_randIndex(i2,L))*mfexp(phi_sigma(L));
	       }
	   }
       phi(j-1)=1/(1+exp(-mu));                                    // compute phi for the interval
       if(tint(i,j-1)!=1)
          phi(j-1)=pow(phi(j-1),tint(i,j-1));                      // adjust phi for the time interval length
	   } else
	      phi(j-1)=phi_fixedDM(i2,kphi+1);
	   /////// p computation /////////
	   if(p_fixedDM(i2,kp+1)== -1)
	   {
	      mu=0;
	      for(L=1;L<=kp;L++)
            mu+=p_fixedDM(i2,L)*p_beta(L);                         // fixed portion of mean 
          p(j-1)=1/(1+exp(-mu));                                    // compute p for the occasion
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
          if(phi_phase>0)
	      {
              if(phi_counts(i) > 0)	                               // random portion of mean if any
	          {
	             for(L=1;L<=phi_krand;L++)
		            if(phi_randIndex(i2,L)>0)
	                   mu+=phi_randDM(i2,L)*phi_u(phi_randIndex(i2,L))*mfexp(phi_sigma(L));
	          }  
	      }
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
          if(p_phase>0)
	      {
              if(p_counts(i) > 0)	                               // random portion of mean if any
	          {
	             for(L=1;L<=p_krand;L++)
		            if(p_randIndex(i2,L)>0)
	                   mu+=p_randDM(i2,L)*p_u(p_randIndex(i2,L))*mfexp(p_sigma(L));
	          }  
	      }
          p(j-1)=1/(1+exp(-mu));                                     
       } else
	     p(j-1)=p_fixedDM(i2,kp+1);
	   }
	report << p << endl;
	}



	

