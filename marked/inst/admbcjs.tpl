// Cormack-Jolly-Seber model; fixed effects only for survival (Phi) and capture probability (p)
// Capture history can represent more than one animal; frq is the frequency for the capture history
// Jeff Laake; 15 Nov 2012
DATA_SECTION 
    init_int n;                        // number of capture histories
    init_int m;                        // number of occasions
    init_imatrix ch(1,n,1,m);          // capture history matrix
    init_ivector frst(1,n);            // occasion first seen for each history
    init_ivector lst(1,n);             // occasion last seen for each history
    init_vector frq(1,n);              // frequency of each history
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
    init_vector phibeta(1,kphi);       // parameter vector for Phi
    init_vector pbeta(1,kp);           // parameter vector for p
    vector phi(1,m);                   // temp vector for Phis for each occasion for a single history
    vector p(1,m-1);                   // temp vector for Phis for each occasion for a single history
    number pch;                        // probability of capture history
    vector phicumprod(1,m);            // cummulative survival probability across occasions
    vector cump(1,m);                  // cummulative probability of being seen across occasions
    vector phix(1,nrows);              // vector of all real survival values
    vector px(1,nrows);                // vector of all real capture probability values    
    objective_function_value f;        // objective function - negative log-likelihood 
PROCEDURE_SECTION
    int i,i1,i2,j;                                         // miscellaneous ints 
    f=0.0;                                                 // initialize to 0     
    phix=1/(1+exp(-phidm*phibeta));                        // compute all phi values using inverse logit
    for(i=1;i<=K;i++)                                      // assign any fixed real phi values
       phix(PhiF(i,1))=PhiF(i,2);
    px=1/(1+exp(-pdm*pbeta));                              // compute all p values using inverse logit
    for(i=1;i<=L;i++)                                      // assign any fixed real p values 
       px(pF(i,1))=pF(i,2); 
    phi=0;                                                 // set all phi values to 0
    for(i=1;i<=n;i++)                                      // loop over each history
    {
        phicumprod=1.0;                                    // set cummulative survival to 1
        cump=1.0;                                          // set cummulative capture prob to 1
        i1=(m-1)*(i-1);                                    // compute beginning index in design matrix
        for(j=frst(i)+1;j<=m;j++)                          // loop over occasions from frst to m
        {
           i2=i1+j;                                        // compute index in design matrix for this occasion
           phi(j-1)=pow(phix(i2-1),tint(i,j-1));           // get phi for the interval adjusted for time length
           p(j-1)=px(i2-1);                                // get p for the occasion
           phicumprod(j)=phicumprod(j-1)*phi(j-1);         // compute cummulative survival
           cump(j)=cump(j-1)*((1-p(j-1))*
                  (1-ch(i,j))+p(j-1)*ch(i,j));             // compute cummulative capture probability
        }   
        pch=0.0;                                           // initialize capture history probability 
        for(j=lst(i);j<=((1-loc(i))*m+loc(i)*lst(i));j++)  // loop over last occasion to m unless loss on capture
        {                                                  //  to compute prob of the capture history 
           i2=i1+j;                                        // index occasion
           if(loc(i)==1)
              pch=pch+cump(j)*phicumprod(j);               // probability of history for loss on captures
           else
              pch=pch+cump(j)*phicumprod(j)*(1-phi(j));    // probability of history given possible last occasion alive
        }   
        if(pch < 1E-15 & frq(i)> 0) pch=1E-307;            // avoid log(0)
        f-=frq(i)*log(pch);                                // sum -log-likelihood frequency*log(pr(ch))
    }
