#' Fitting function for Multivariate Multistate CJS with uncertainty models
#' 
#' A function for computing MLEs for MVMS models following Johnson et al (2015)  via ADMB. It works very much like mscjs but with more parameters. While this function is for fitting the model via ADMB, the
#' documentation contained here is also relevant when this type of model is fitted with optimx and FORTRAN code (use.admb=F).
#'  
#' It is easiest to call this function or any model built in \code{marked} through the function \code{\link{crm}}.
#' This function has been exported but to fit a model it should be called through the crm function
#' which does some testing of the arguments to avoid errors and sets up the calling arguments appropriately.
#' 
#' The mvmscjs model is quite flexible but also quite complex. Prior to using this model, read Johnson et al. (2015) 
#' and particularly Supplement B which goes through the code used in the analysis of Johnson et al. (2015).
#' Supplement B is useful but there have been a number of changes to the code since the paper was written and it is far from complete.
#' The documentation here will fill in some of those blanks.
#'   
#' There are 5 classes of parameters in an mvms model. They are 1)Phi - survival probability, 2) p-capture probability, 3) Psi - state transition
#' probabilities, 4) delta - certainty/uncertainty state probabilities, 5) pi - initial state probability. The final
#' class pi was described in Johnson et al. (2015) but was not implemented in \code{marked} at the time the paper was written. However, with 
#' the sealions example, the initial state is known for all releases and pi is not relevant.
#' 
#' Before I jump into a description of each parameter, I'll describe some characteristics that are common to all of the 
#' parameters. First, all of the parameters are probabilities and thus must be in the unit interval. To accomplish this
#' each of the parameters uses either a logit or mlogit (multinomial logit) link function.
#' 
#' A logit link is log(theta/(1-theta))= a+bx... where theta is the probability (called a real parameter in \code{MARK} terminology) 
#' and a+bx... is the set of covariates (e.g. x) and working parameters (a,b) that are unbounded.  Working parameters are called beta parameters in 
#' \code{MARK} terminology. The inverse logit link is  1/(1+exp(-a-bx)) for this simple example and it can be computed from
#' the plogis function in R (plogis(a+bx)). A largish negative value of the link (eg -5) mean theta will be nearly 0 and
#' and a largish positive value (+5) means theta is nearly 1.
#' 
#' The mlogit link is simply an extension of the logit link for more than 2 classes.  For example, the logit link for survival is
#' for the 2 classes alive/dead.  I'll use the logit link as a simple example to describe the way mlogit links are formed in \code{marked}.
#' If survival is Phi=1/(1+exp(-a))=exp(a)/(1+exp(a)) then mortality is 1-Phi = 1-exp(a)/(1+exp(a))= 1/(1+exp(a)). Imagine that we
#' have 2 link values a - for survival and b - for mortality.  Then we can rewrite the inverse logit function for Phi as exp(a)/(exp(b)+exp(a))
#' and then 1-Phi = exp(b)/(exp(b)+exp(a)). It is obvious that Phi+1-Phi=1 which must be the case.  Now, if you actually tried to fit a model
#' with both a and b the parameters a and b are not uniquely identifiable because there is only one parameter Phi. One of the values must be a 
#' reference cell which means the value is fixed.  Typically we assume b=0 which yields exp(0)=1 and we get the logit link.
#' 
#' Now consider 3 classes (A,B,C) and we'll use the link values a,b,c for a generic parameter theta. As with the 2 class example, the
#' probabilities are theta(A)=exp(a)/(exp(a)+exp(b)+exp(c)),  theta(B)=exp(b)/(exp(a)+exp(b)+exp(c)) and theta(C)=exp(c)/(exp(a)+exp(b)+exp(c)).
#' As before theta(A)+theta(B)+theta(C)=1 but for a multinomial with 3 classes there are only 2 free estimable parameters due to the constraint that they have to sum to 1.
#' Just as with the example above, a reference cell must be selected.  If we chose B as the reference cell then b=0 and the probabilities become
#' theta(A)=exp(a)/(exp(a)+1+exp(c)),  theta(B)=1/(exp(a)+1+exp(c)) and theta(C)=exp(c)/(exp(a)+1+exp(c)). This extends to an k-class multinomial where
#' k-1 parameters are freely estimable and 1 must be fixed.
#' 
#' For parameters with logit links there is an obvious choice for the reference cell.  For example, for survival probability we are interested in
#' describing the parameters that affect survival, so we use death as the reference cell which is computed by subtraction (1-Phi). Likewise, for
#' capture (sighting) probability we use not caught (1-p) as the reference cell.  For mlogit links, there is not always an obvious choice for a
#' reference cell and it can be an advantageous to allow the reference cell to be flexible and user-defined.
#' 
#' To understand how to set the reference cell, you should review the material in the help file for create.dmdf and Supplement B of Johnson et al. (2015) 
#' which describes the design data but I'll give a brief description here. All model building in the \code{marked} package revolves around the design data. 
#' If the design data are incorrect the model will be incorrect. Briefly, the design data contain the data for the parameter model, 
#' providing the maximum structure from which specific models can be built.. For example, for p there are M-1 records for each of the I individuals and M occasions
#' (note that if the data are accumulated then an "individual" is actually a unique set of capture history and covariate data that represents "freq" animals). 
#' Thus, capture probability could be allowed to differ for each occasion, individual, age, strata, or any combination thereof.  In practice, such a model would
#'  be overparameterized, but this structure allows maximum flexibility; users still need to use formula to provide specific constraints to reduce the number of 
#' parameters to a reasonable level. With this design data approach essentially there is a place holder for each possible p and any or all of those values 
#' can be fixed to a real value. 
#' 
#' "Real" parameters can be fixed by adding a field called fix to the design data for a parameter. The value of fix should be either NA if 
#' the parameter is to be estimated or it should be the value of the real parameter (inverse link value). For example, we could fix the capture probability 
#' for time 2 to be 0 for all individuals at time 2 or at time 2 in strata B because there was no capture effort.  Or you could fix p=0 for an individual at a time
#' because it was temporarily removed. Or you could fix p=1 for a set of individuals that were always recaptured via telemetry.
#' 
#' Now let's go back to the mlogit parameters. There is an important difference between the way \code{marked} and \code{RMark/MARK} work 
#' with regard to mlogit parameters like Psi for state transition probabilities. An mlogit parameter is one in which the 
#' sum of the probabilities is 1. For Psi, if I'm in stratum A and can go to B or C or remain in A, the probabilities A to A, A to B and A to C must sum to 1 because that is
#' all of the possibilities. In \code{RMark/MARK} the design data would only contain 2 records which are determined 
#' based on what you select as subtract.stratum. If subtract.stratum was set as A for the A stratum, the design 
#' data for Psi would only contain records for A to B and A to C. The value for A to A would be computed by 
#' subtraction and is the reference cell.
#' 
#' By contrast, in \code{marked}, all 3 records are in the design data and to set the transition from A to A as the reference cell (computed by subtraction) 
#' the value of fix is set to 1. So if fix is the inverse logit value why is it being
#' set to 1 for the mlogit?  Because the link for these parameters is actually a log-link and the inverse link is exp(link).
#' So the link value is fixed to 0 and the inverse link ("real value") is exp(0)=1. The mlogit parameter is constructed
#' by dividing each exp() value by the sum of the exp() values across the set defining the possible classes for the multinomial. 
#' Following along with the example, 
#' 
#' \preformatted{
#' Psi_AtoB=exp(beta_AtoB)/(1+exp(beta_AtoB)+exp(beta_AtoC))
#' Psi_AtoC=exp(beta_AtoC)/(1+exp(beta_AtoB)+exp(beta_AtoC)) and
#' Psi_AtoA=1-Psi_AtoB-Psi_AtoC = 1/(1+exp(beta_AtoB)+exp(beta_AtoC)). 
#' }
#' The "appropriate" set will depend on the mlogit parameter and the data/model structure.
#' 
#' I used this structure for several reasons. First, you get a real parameter estimate for the subtracted stratum which you don't get in \code{RMark/MARK}. 
#' Secondly, you can change the value to be subtracted at will and it is not fixed across the entire model fit,
#' but you do have to be careful when specifying the model when you do that because the formula specifies the 
#' parameters for those that are not fixed. Finally, you can change the reference cell simply by modifying the value of fix, 
#' whereas with the subtract.stratum approach the design data would have to be recreated which can take some time with large data sets.
#'  
#' Next I'll briefly describe each parameter in the model and its link function structure.
#' 
#' Phi
#' -----
#' Phi is survival probability in the model. Survival probability is also called S in the straight multi-state model. In retrospect I wish I had stayed with 
#' S but I didn't.  In capture-recapture Phi is usually reserved for apparent survival rate which might include permanent emigration. 
#' With a multi-state model, the presumption is that "all" states are covered and animals can't emigrate permanently to some 
#' unobservable state.  Whether that is true or not is always a question. Anyhow I use Phi for survival in \code{mvmscjs} models.
#' 
#' Phi uses a logit link (log(Phi/(1-Phi))= a+bx...) to relate survival to any covariates which have parameter values
#' across the entire real line.  Phi is what I call an interval parameter in that it describes the probability of surviving 
#' during an interval of time between sampling occasions. The design data for interval parameters are describe in terms of
#' the time and age of the animal at the beginning of the interval.  For example, if occasions are at time 1,2,...J there
#' will be J-1 data records for each individual with time values 1,...,J-1.
#'  
#' p
#' -----
#' p is recapture/resighting probability.  While I'll often fall back and call it capture or sighting probability, technically it is
#' recapture or resighting because this model is based on a Cormack-Jolly-Seber (CJS) formulation which does not model initial "capture"
#' probability and only models events after they are marked and released. p also uses the logit link and it is an occasion parameter because
#' it describes an event associated with an occasion. Because they are "recaptures", there are J-1 design data records for
#' each individual at times 2,...,J. 
#' 
#' An example of using fix for p is given in the sealions example.  For the data, there was no survey effort in stratum A for the last occasion.
#' \preformatted{
#' ddl$p$fix = ifelse(ddl$p$Time==17 & ddl$p$area=="A", 0, ddl$p$fix)
#' }
#' 
#' Psi
#' -----
#' Psi describes the transition probabilities of moving between states. For each individual at each occasion for each of the possible states,
#' it is formulated as a multinomial with a cell probability for movement from the state to all of the other possible states and remaining in the
#' in the same state. Thus if there are M states, then there are M*M records for each individual for each of the J-1 intervals. It is treated as
#' interval parameter with movement occurring between occasions but after survival.  If an animal is in state A it survives and moves to state B with
#' probability Phi_A*Psi_AtoB. If those were to be reversed you could use a trick proposed by Bill Kendall of adding dummy occasions between each real set
#' of occasions and fixing real parameters. For example using the first 2 occasions with data A0, the data would become A00 with 3 occasions. For the
#' first interval all the Phis would be fixed to 1 and Psi_AtoB would be estimated. For occasion 2, all the ps would be fixed to 0 because there are no observations and
#' for the second interval all the Psi's should be fixed so it will remain in the same state (eg. Psi_AtoA=1 and Psi_Atox=0).
#' 
#' This can get rather involved when you have more than 1 variable describing the states which is the whole purpose of the
#' multi-variate multi-state model. It just needs some careful thought.  As an example, here are some of the columns of the design data for the first interval for
#' the first individual for stratum A-- which for this example is area A with missing left and right tags. There are 8 records because
#' there are 8 possible strata - 2 areas *2 left tag states * 2 right tag states.
#' 
#' \preformatted{
#' > head(ddl$Psi[,c(1:10,20),],8)
#'  id occ stratum tostratum area ltag rtag toarea toltag tortag fix
#'   1   1     A--       A--    A    -    -      A      -      -   1
#'   1   1     A--       A-+    A    -    -      A      -      +   0
#'   1   1     A--       A+-    A    -    -      A      +      -   0
#'   1   1     A--       A++    A    -    -      A      +      +   0
#'   1   1     A--       S--    A    -    -      S      -      -  NA
#'   1   1     A--       S-+    A    -    -      S      -      +   0
#'   1   1     A--       S+-    A    -    -      S      +      -   0
#'   1   1     A--       S++    A    -    -      S      +      +   0
#' }
#'  
#' When the design data are created, it creates fields based on the definition of strata.labels in the call to process.data. For this 
#' example it was  strata.labels=list(area=c("A","S"),ltag=c("+","-","u"),rtag=c("+","-","u"))).  The design data include a
#' stratum and tostratum field which are the concatenated values of all the variables defining the states. In addition field and tofield
#' variables are created for each strata variables.  Here the fields are area,ltag,rtag and the tofields are toarea,toltag and tortag.
#' Each simply partitions the stratum and tostratum fields.  Note that if only a single variable is used to define the strata.labels the
#' the design data will only contain a stratum and tostratum fields.
#' 
#' Now, when the design data are initially created for Psi, it creates the field fix and assigns the value NA to all values except for
#' the record in which stratum=tostratum which gets the value of fix=1 making it the reference cell. Those are the default values and it is up to you to change as needed.
#' In this case almost all of the Psi values are fixed to 0 because they cannot occur. The only possible movement is from A-- to S-- because
#' all of the other transitions would have it gaining tags that have already been lost.  As another example, the stratum A++ would have all of
#' the values of fix=NA except for tostratum=A++ which would be fixed to 1.
#' 
#' delta
#' -----
#' delta is the probability of being certain or uncertain about each state variable for a given stratum. The definition can change based on 
#' how you define the reference cell. For delta, the number of records depends on the structure of strata.labels and the number of 
#' state variables and how many state variables can be uncertain.  Below are the design data for a single stratum for a single individual and occasion from the sealions
#' example.
#' 
#' \preformatted{
#' 
#'id occ stratum area ltag rtag obs.area obs.ltag obs.rtag obs.ltag.u obs.rtag.u
#'1   1     A--    A    -    -        A        -        -          0          0
#'1   1     A--    A    -    -        A        -        u          0          1
#'1   1     A--    A    -    -        A        u        -          1          0
#'1   1     A--    A    -    -        A        u        u          1          1
#'}
#' For this example, the area is certain but the left tag and right tag status can be uncertain. The above data are for the stratum
#' A--. Each field used to define the states has a field value and an obs.field value (eg., ltag, obs.ltag). The fields obs.ltag.u and obs.rtag.u were created in code 
#' after the call to make.design.data (see ?sealions). There would be a set like these for each of the 8 possible stratum value for id=1 and occ=1. Because each tag status can be either 
#' known or uncertain there are 2*2 records and the sum of the 4 cell probabilities must sum to 1. The easiest way to do that is
#' by specifying fix=1 for one of the cells and assigning the other 3 fix=NA so they are estimated. If you set fix=1 for the first record (observation A--),
#' then the parameters for delta are describing the probability of being uncertain. Whereas, if you set fix=1 for the fourth record (observation Auu) then
#' the parameters for delta are describing the probability of being certain.
#' 
#' Now the use of fix is not required to set a reference cell because you can accomplish the same thing by defining the formula for
#' delta such that the link value is 0 for one of the cells.  For example, the formula for delta in the example is
#'  \preformatted{
#'   ~ -1 + obs.ltag.u + obs.rtag.u + obs.ltag.u:obs.rtag.u)
#' }
#' The numeric value of obs.ltag.u + obs.rtag.u + obs.ltag.u*obs.rtag.u (: equivalent to mulitplication here) is 0 for the first record and -1 in the
#' formula specifies that there shouldn't be an intercept. Thus regardless of the 3 beta values for delta, the link value for the 
#' first record is 0 and the inverse link value is exp(0)=1 which is equivalent to setting fix=1 for the first record to set it as the reference cell.
#' If you run the sealions example with the current version of the software it will give the message:
#' \preformatted{
#'  No values provided for fix for delta. At least need to set a reference cell
#' }
#' Messages like this are warnings and do not stop the model from being fitted because as shown above it is possible to set the reference cell with an appropriate formula.
#' If a reference cell is not set the beta estimates are not all identifiable.  Symptoms include troubles with optimization, different starting values giving different
#' beta estimates for the parameter with the same likelihood value, and most easily seen is large standard errors for all of the beta estimates involved.
#' 
#' pi
#' -----
#' pi is the probability of being in a state for the initial occasion. It is only used when there is one or more individual capture histories have an uncertain (u) 
#' state variable value for the initial occasion of the capture history. For the sealions example, all are released in the state S++ so there was no uncertainty.
#' 
#' The design data for pi includes a record for each possible stratum for each individual.  For the sealions example, the following are the
#' design data for the first individual:
#' 
#' \preformatted{
#'  id occ stratum area ltag rtag fix
#'   1   1     A--    A    -    -   0
#'   1   1     A-+    A    -    +   0
#'   1   1     A+-    A    +    -   0
#'   1   1     A++    A    +    +   0
#'   1   1     S--    S    -    -   0
#'   1   1     S-+    S    -    +   0
#'   1   1     S+-    S    +    -   0
#'   1   1     S++    S    +    +   1
#' }
#' 
#' When the design data are created the code recognizes that the initial value for this individual capture history is known for all of
#' the strata values so it assigns 0 to the strata that are not the initial value and 1 for the strata that is the initial value.
#' 
#' To demonstrate what is obtained when there is uncertainty I modified one of the S++ values to Suu and this is the resulting design data 
#'\preformatted{
#'  id occ stratum area ltag rtag fix
#'  335   1     A--    A    -    -   0
#'  335   1     A-+    A    -    +   0
#'  335   1     A+-    A    +    -   0
#'  335   1     A++    A    +    +   0
#'  335   1     S--    S    -    -   1
#'  335   1     S-+    S    -    +  NA
#'  335   1     S+-    S    +    -  NA
#'  335   1     S++    S    +    +  NA
#' }
#' 
#' All of the fix values for area A are still 0 because area is certain and it is in area S.  However, each tag status is uncertain
#' so there are 4 possible true values. S-- is set to the reference cell by default and the other 3 would be estimated.
#' 
#' There is a subtlety to the definition of pi that should be understood. The structure of the default design data for pi means that pi is the conditional 
#' probability of being in a state given there is uncertainty and that an animal is first observed in that period. In other words, it is the probability of 
#' being in states just for the records that have uncertainty and are first observed at a given time. For the sealions, consider having released 1/4 of the 
#' animals in A with all certain tag status and for the 3/4 released in S neither tag status was known. In that case pi would only relate to S and the proportion
#'  given both tags will possibly be less in S than in both areas combined. However, if you defined fix for all records as above the estimated pi would apply to 
#' the whole set. 
#'  
#' The mvmscjs model is constructed as a hidden Markov model as described in Johnson et al. (2015).  As such the 5 parameter classes (Phi,p,Psi,delta and pi)
#' are used to construct the matrices gamma (eq 2.13 in Johnson et al. (2015)),  dmat (eq 2.14 in Johnson et al. (2015)) and delta which is called pi in Johnson et al. (2015).
#' In the code I used the name delta for the matrix which has I rows and M colums because the symbol delta is used in Zucchini et al's book on HMMs. So think of delta
#' as the matrix and pi as the values in each row of the matrix. I should have used pi as the matrix as well to avoid confusion between the delta parameters and the
#' matrix delta for initial values.
#' 
#' One other difference between the coding and Johnson et al. (2015) is that
#' the 0 observation in dmat is the last row in Jonhson et al (2015) and it is the first row in the code.  Likewise, for gamma the death state
#' is the last row in Johnson et al. (2015) and the first row in the code. In Johnson et al. (2015) dmat and gamma are matrices but in the code they are
#' 4-d arrays where the first 2 dimensions are individuals and occasions and the last 2 dimensions are the matrices defined in Johnson et al. (2015).
#' 
#' Johnson et al. (2015) show delta being used for each occasion including the initial release occasion.  In fact, design data for delta are created for
#' each occasion. But in constructing the formula there is one exception that you need to understand.  If the initial values are either all certain or
#' all uncertain then formula for delta is not used for that initial release because all values of dmat are either 0 or 1.
#' 
#' If these matrices are setup incorrectly the optimization will likely fail. In general, that is all handled by the code but it is
#' possible to mess up the values in those matrices if the values of fix are incorrectly defined.  One place where this can occur is if
#' the value of fix=0 for all records defining a multinomial because when it normalizes to sum to 1, it will produce a NaN due to 0/0.
#' To prevent problems the function HMMLikelihood with the argument check=TRUE is run before attempting the optimization. With the starting 
#' parameter values it checks to make sure that the columns of dmat sum to 1 and the rows of gamma sum to 1 and all values of pi in delta sum to 1.
#' If any of these fail outside of rounding error, the program will issue errors showing the problematic records and values and will stop.
#' 
#' Unobservable states can be included in the data structure and model but you must be careful in constructing formulas. Keep in
#' mind that there are no observations for unobservable states. That states the obvious but it can be easily forgotten when constructing
#' formulas.  For an easy example, consider a univariate example with area defining states.  Will assume that we have areas A,B and C that
#' are observable and X which is unobservable.  For Psi you can have transitions to and from the X and the other states. Obviously, you'll need to
#' set fix=0 for p for X because it can't be seen. Also, if you were to set the formula for Phi to be ~stratum there are no data in X so it will not be able to
#' estimate survival in X. So don't do that.  Likewise, the unobservable state value X should not be in a formula for delta because it represents
#' certainty/uncertainty of state variables for observations and there are no observations in X. Also, because this is a CJS formulation that is
#' conditional on the first release (observation), the formula for pi cannot include X, again because there are no observations. 
#'
#' To avoid including parameters for levels of a factor for unobservable state values, the easiest approach is to set fix values for the
#' design data rows for the unobservable states.  For example, if fix=0 for p for stratum X, then ~stratum could be used in the formula
#' and it would have an intercept and parameters for B and C.  The stratum level X would be excluded because all rows of the resulting
#' design matrix with fixed values are set to 0 and then any columns which have all 0 values are dropped which would drop the X factor level.
#' Another approach is to define a different variable that excludes the factor level. Note that this approach would not work to set p=0 for 
#' the unobservable state. Assume in addition to area we had another variable to define states and it could be uncertain. If we thought
#' that the uncertainty in that variable was related to area, we could use a formula for delta like ~-1+notX:area where notX was a variable
#' with value 1 when the area was not X and was 0 when area was X.  This would result in parameters for A,B and C but not X. A similar approach
#' would be to create 0/1 dummy variables A,B,C which were 1 when the record was for area=A,B,C respectively and 0 otherwise. Then you could
#' use the delta formula ~-1+A+B+C specifically excluding X.  It is not important for X to have a value for delta because there are no observations
#' for X and if you look at the last row (observation =0) in eq 2.14 in Johnson et al. (2015) you'll see is 1-p and doesn't use delta.
#' 
#' \code{marked} doesn't include any specific robust design models because they can be relatively easily handled by fixing real parameters.
#' In a standard closed multi-state robust design, during the secondary periods Phi=1 and Psi_AtoA= 1 for all states (no transitions).
#' Any variation can be created simply by fixing the appropriate parameters. For example, with the sealions example, you might allow
#' transitions between areas but assume no transitions for tag status.
#'
#' 
#' @param x processed dataframe created by process.data
#' @param ddl list of dataframes for design data; created by call to
#' \code{\link{make.design.data}} and then simplified
#' @param fullddl list of dataframes for design data prior to simplification
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param model_data a list of all the relevant data for fitting the model including
#' imat, Phi.dm,p.dm,Psi.dm,Phi.fixed,p.fixed,Psi.fixed and time.intervals. It is used to save values
#' and avoid accumulation again if the model was re-rerun with an additional call to cjs when
#' using autoscale or re-starting with initial values.  It is stored with returned model object.
#' @param parameters equivalent to \code{model.parameters} in \code{\link{crm}}
#' @param accumulate if TRUE will accumulate capture histories with common
#' value and with a common design matrix for S and p to speed up execution
#' @param initial list of initial values for parameters if desired; if each is a named vector
#' from previous run it will match to columns with same name
#' @param method method to use for optimization; see \code{optim}
#' @param hessian if TRUE will compute and return the hessian
#' @param debug if TRUE will print out information for each iteration
#' @param chunk_size specifies amount of memory to use in accumulating capture
#' histories; amount used is 8*chunk_size/1e6 MB (default 80MB)
#' @param refit non-zero entry to refit
#' @param itnmax maximum number of iterations
#' @param control control string for optimization functions
#' @param scale vector of scale values for parameters
#' @param re if TRUE creates random effect model admbcjsre.tpl and runs admb optimizer
#' @param compile if TRUE forces re-compilation of tpl file
#' @param extra.args optional character string that is passed to admb if use.admb==TRUE
#' @param clean if TRUE, deletes the tpl and executable files for amdb if use.admb=T
#' @param sup supplemental index values for constructing mvms model
#' @param ... not currently used
#' @export
#' @import R2admb
#' @return The resulting value of the function is a list with the class of
#' crm,cjs such that the generic functions print and coef can be used.
#' \item{beta}{named vector of parameter estimates} \item{lnl}{-2*log
#' likelihood} \item{AIC}{lnl + 2* number of parameters}
#' \item{convergence}{result from \code{optim}; if 0 \code{optim} thinks it
#' converged} \item{count}{\code{optim} results of number of function
#' evaluations} \item{reals}{dataframe of data and real S and p estimates for
#' each animal-occasion excluding those that occurred before release}
#' \item{vcv}{var-cov matrix of betas if hessian=TRUE was set}
#' @author Jeff Laake 
#' @references Johnson, D. S., J. L. Laake, S. R. Melin, and DeLong, R.L. 2015. Multivariate State Hidden Markov Models for Mark-Recapture Data. 31:233-244.
mvmscjs=function(x,ddl,fullddl,dml,model_data=NULL,parameters,accumulate=TRUE,initial=NULL,method,
		hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale,
		re=FALSE,compile=FALSE,extra.args="",clean=TRUE,sup,...)
{
	accumulate=FALSE
	nocc=x$nocc
	xstart=x$start
#  Time intervals has been changed to a matrix (columns=intervals,rows=animals)
#  so that the initial time interval can vary by animal; use x$intervals if none are in ddl$Phi
	if(!is.null(ddl$Phi$time.interval))
	  time.intervals=matrix(fullddl$Phi$time.interval[fullddl$Phi$stratum==x$strata.labels[1]],nrow(x$data),ncol=nocc-1,byrow=TRUE)
	else
	if(is.vector(x$time.intervals))
		time.intervals=matrix(x$time.intervals,nrow=nrow(x$data),ncol=nocc-1,byrow=TRUE)
	else
		time.intervals=x$time.intervals
	chmat=x$ehmat	
#  Store data from x$data into x
	strata.labels=x$strata.labels
	uS=x$unobserved
	x=x$data
#  set default frequencies if not used
	freq=NULL
	if(!is.null(x$freq))freq=x$freq
#  get first and last vectors, loc and chmat with process.ch and store in imat
	ch=x$ch
	imat=process.ch(ch,freq,all=FALSE)
#	chmat=matrix((unlist(strsplit(ch,","))),byrow=TRUE,ncol=nocc,nrow=length(ch))
#	for(nlabel in 1:length(strata.labels))
#		chmat=t(apply(chmat,1,sub,pattern=strata.labels[nlabel],replacement=nlabel))
#  Use specified initial values or create if null
	if(is.null(initial))
	{
		par=list(Psi=rep(0,ncol(dml$Psi$fe)),
				p=rep(0,ncol(dml$p$fe)),
				Phi=rep(0,ncol(dml$Phi$fe)))
		if(ncol(dml$pi$fe)>0) par$pi=rep(0,ncol(dml$pi$fe))
		if(ncol(dml$delta$fe)>0) par$delta=rep(0,ncol(dml$delta$fe))
	}
	else
		par=set.initial(names(dml),dml,initial)$par
#  Create list of model data for optimization
	model_data=list(Phi.dm=dml$Phi$fe,p.dm=dml$p$fe,Psi.dm=dml$Psi$fe,delta.dm=dml$delta$fe,pi.dm=dml$pi$fe,imat=imat,Phi.fixed=parameters$Phi$fixed,
			p.fixed=parameters$p$fixed,Psi.fixed=parameters$Psi$fixed,delta.fixed=parameters$delta$fixed,
			pi.fixed=parameters$pi$fixed,time.intervals=time.intervals)
#   If data are to be accumulated based on ch and design matrices do so here;
	if(accumulate)
	{
		cat("Accumulating capture histories based on design. This can take awhile.\n")
		flush.console()
		model_data.save=model_data   
		#model_data=mscjs.accumulate(x,model_data,nocc,freq,chunk_size=chunk_size)
	}else
		model_data.save=model_data
#   Create links  -- not used at present; idea here is to use sin links for parameters where you can   
#   Phi.links=create.links(Phi.dm)
#   Phi.links=which(Phi.links==1)
#   p.links=create.links(p.dm)
#   p.links=which(p.links==1)
#  Scale the design matrices and parameters with either input scale or computed scale
	scale=1
	scale=set_scale(names(dml),model_data,scale)
	model_data=scale_dm(model_data,scale)
	# setup tpl to be multistate.tpl 
	if(!re)
		tpl="mvms"
	else
		stop("random effect portion not completed for this model")
	# setup admb exe and cleanup old files and previous tpl; checks for exe 
	# in package directory and if found uses it otherwise compiles tpl file
	setup_admb(tpl,compile,clean,re=FALSE)
	# create .dat file to write its contents 
	con=file(paste(tpl,".dat",sep=""),open="wt")
	# Number of observations
	n=length(model_data$imat$freq)
	write(n,con,append=FALSE)
	# Number of occasions
	nocc=model_data$imat$nocc
	write(nocc,con,append=TRUE)
	# Number of states
	nS=length(strata.labels)
	write(nS,con,append=TRUE)
	# Number of possible observations
    nobs=length(sup$obslevels)
	write(nobs,con,append=TRUE)
	# Number of delta data records per id-occasion
	write(sup$np,con,append=TRUE)
	
	# capture history matrix
	write(t(chmat),con,ncolumns=nocc,append=TRUE)
	# first occasions seen 
	write(model_data$imat$first,con,ncolumns=n,append=TRUE)
	# frequency of capture history 
	if(!re)
	{
		write(model_data$imat$freq,con,ncolumns=n,append=TRUE)
	} else
	{
		if(any(model_data$imat$freq!=1))stop("\n cannot use random effects with frequency >1")
	}
	write(t(model_data$time.intervals),con,ncolumns=nocc-1,append=TRUE)

	# Phi design matrix
    # zero out rows with fixed parameters and remove any unneeded columns
    if(!is.null(ddl$Phi$fix))
        model_data$Phi.dm[!is.na(ddl$Phi$fix),]=0
    if(ncol(model_data$Phi.dm)!=0)
    {
	    select=vector("logical",length=ncol(model_data$Phi.dm))
	    for (i in 1:ncol(model_data$Phi.dm))
	    	select[i]=any(model_data$Phi.dm[,i]!=0)
	    model_data$Phi.dm=model_data$Phi.dm[,select,drop=FALSE]
    }
	phidm=model_data$Phi.dm
	phifix=rep(-1,nrow(phidm))
	if(!is.null(ddl$Phi$fix))
		phifix[!is.na(ddl$Phi$fix)]=ddl$Phi$fix[!is.na(ddl$Phi$fix)]
	slist=simplify_indices(cbind(phidm,phifix))
	write(ncol(phidm),con,append=TRUE)
	write(length(slist$set),con,append=TRUE)
	if(ncol(phidm)>0)
		write(t(phidm[slist$set,,drop=FALSE]),con,ncolumns=ncol(phidm),append=TRUE)
	write(phifix[slist$set],con,append=TRUE)
	write(slist$indices[ddl$Phi.indices],con,append=TRUE)

	# p design matrix
    # zero out rows with fixed parameters and remove any unneeded columns
    if(!is.null(ddl$p$fix))
       model_data$p.dm[!is.na(ddl$p$fix),]=0
    if(ncol(model_data$p.dm)!=0)
    {
	    select=vector("logical",length=ncol(model_data$p.dm))
	    for (i in 1:ncol(model_data$p.dm))
	    	select[i]=any(model_data$p.dm[,i]!=0)
	    model_data$p.dm=model_data$p.dm[,select,drop=FALSE]
    }
	pdm=model_data$p.dm
	pfix=rep(-1,nrow(pdm))
	if(!is.null(ddl$p$fix))
		pfix[!is.na(ddl$p$fix)]=ddl$p$fix[!is.na(ddl$p$fix)]
	slist=simplify_indices(cbind(pdm,pfix))
	write(ncol(pdm),con,append=TRUE)
	write(length(slist$set),con,append=TRUE)
	if(ncol(pdm)>0)
		write(t(pdm[slist$set,,drop=FALSE]),con,ncolumns=ncol(pdm),append=TRUE)
	write(pfix[slist$set],con,append=TRUE)
	write(slist$indices[ddl$p.indices],con,append=TRUE)
    # write out indices for completing dmat
    write(nrow(sup$indices_forp),con,append=TRUE)
    write(t(sup$indices_forp),con,append=TRUE)

	# Psi design matrix
	# zero out subtracted stratum (fixed) and remove any unneeded columns
    if(!is.null(ddl$Psi$fix))
	   model_data$Psi.dm[!is.na(ddl$Psi$fix),]=0
	if(ncol(model_data$Psi.dm)!=0)
	{
		select=vector("logical",length=ncol(model_data$Psi.dm))
		for (i in 1:ncol(model_data$Psi.dm))
			select[i]=any(model_data$Psi.dm[,i]!=0)
		model_data$Psi.dm=model_data$Psi.dm[,select,drop=FALSE]
	}
	psidm=model_data$Psi.dm
	psifix=rep(-1,nrow(psidm))
	if(!is.null(ddl$Psi$fix))
		psifix[!is.na(ddl$Psi$fix)]=ddl$Psi$fix[!is.na(ddl$Psi$fix)]
	slist=simplify_indices(cbind(psidm,psifix))
	write(ncol(psidm),con,append=TRUE)
	write(length(slist$set),con,append=TRUE)
	if(ncol(psidm)>0)
	   write(t(psidm[slist$set,,drop=FALSE]),con,ncolumns=ncol(psidm),append=TRUE)
	write(psifix[slist$set],con,append=TRUE)
	write(slist$indices[ddl$Psi.indices],con,append=TRUE)

	# delta design matrix
    if(!is.null(ddl$delta$fix))
	    model_data$delta.dm[!is.na(ddl$delta$fix),]=0
	if(ncol(model_data$delta.dm)!=0)
	{
		select=vector("logical",length=ncol(model_data$delta.dm))
		for (i in 1:ncol(model_data$delta.dm))
			select[i]=any(model_data$delta.dm[,i]!=0)
		model_data$delta.dm=model_data$delta.dm[,select,drop=FALSE]
	}
    deltadm=model_data$delta.dm
    deltafix=rep(-1,nrow(deltadm))
    if(!is.null(ddl$delta$fix))
	   deltafix[!is.na(ddl$delta$fix)]=ddl$delta$fix[!is.na(ddl$delta$fix)]
    slist=simplify_indices(cbind(deltadm,deltafix))
    write(ncol(deltadm),con,append=TRUE)
    write(length(slist$set),con,append=TRUE)
	if(ncol(deltadm)>0)
       write(t(deltadm[slist$set,,drop=FALSE]),con,ncolumns=ncol(deltadm),append=TRUE)
    write(deltafix[slist$set],con,append=TRUE)
    write(slist$indices[ddl$delta.indices],con,append=TRUE)

	# pi design matrix
    if(!is.null(ddl$pi$fix))
	    model_data$pi.dm[!is.na(ddl$pi$fix),]=0
	if(ncol(model_data$pi.dm)!=0)
	{
		select=vector("logical",length=ncol(model_data$pi.dm))
		for (i in 1:ncol(model_data$pi.dm))
			select[i]=any(model_data$pi.dm[,i]!=0)
		model_data$pi.dm=model_data$pi.dm[,select,drop=FALSE]
	}
    pidm=model_data$pi.dm
    pifix=rep(-1,nrow(pidm))
    if(!is.null(ddl$pi$fix))
	   pifix[!is.na(ddl$pi$fix)]=ddl$pi$fix[!is.na(ddl$pi$fix)]
    slist=simplify_indices(cbind(pidm,pifix))
    write(ncol(pidm),con,append=TRUE)
    write(length(slist$set),con,append=TRUE)
	if(ncol(pidm)>0)
		write(t(pidm[slist$set,,drop=FALSE]),con,ncolumns=ncol(pidm),append=TRUE)
    write(pifix[slist$set],con,append=TRUE)
    write(slist$indices[ddl$pi.indices],con,append=TRUE)

	# unkinit set to 0 unless all unknown or all known at initial release;
	# when unkinit=0, delta is applied in dmat
	unkinit=as.numeric(all(is.na(xstart[,1])) | all(!is.na(xstart[,1])))
    write(unkinit,con,append=TRUE)

   if(!debug)
	   write(0,con,append=TRUE)
   else
	   write(1,con,append=TRUE)
   
	close(con)
#   write out initial values for betas
	con=file(paste(tpl,".pin",sep=""),open="wt")
	if(ncol(dml$Phi$fe)>0) 
		write(par$Phi,con,ncolumns=length(par$Phi),append=FALSE)
	if(ncol(dml$p$fe)>0) 
		write(par$p,con,ncolumns=length(par$p),append=TRUE)
	if(ncol(dml$delta$fe)>0) 
	  write(par$delta,con,ncolumns=length(par$delta),append=TRUE)
	if(ncol(dml$Psi$fe)>0) 
		write(par$Psi,con,ncolumns=length(par$Psi),append=TRUE)
	if(ncol(dml$pi$fe)>0) 
		write(par$pi,con,ncolumns=length(par$pi),append=TRUE)
	close(con)   
	if(hessian)
		xx=run_admb(tpl,extra.args=extra.args)
	else
		xx=run_admb(tpl,extra.args=paste(extra.args,"-nohess"))
	convergence=attr(xx,"status")
	if(is.null(convergence))convergence=0
	res=read_admb(tpl)
	beta=list(unscale_par(c(res$coeflist$phibeta,res$coeflist$pbeta,res$coeflist$dbeta,res$coeflist$psibeta,res$coeflist$pi),scale))
	parnames=names(unlist(beta))
	fixed.npar=length(unlist(beta))
	if(!is.null(res$hes))
	{
		beta.vcv=solvecov(res$hes)$inv
		rownames(res$hes)=parnames
		colnames(res$hes)=rownames(res$hes)
		if(all(diag(beta.vcv>0))) 
			res$cor=beta.vcv/outer(sqrt(diag(beta.vcv)),sqrt(diag(beta.vcv)))
	}  else
		beta.vcv=res$vcov
	rownames(beta.vcv)=parnames
	colnames(beta.vcv)=rownames(beta.vcv)
	rownames(res$cor)=rownames(beta.vcv)
	colnames(res$cor)=rownames(beta.vcv)
	res$vcov=NULL
	optim.details=c(fn=res$fn,maxgrad=res$maxgrad,eratio=res$eratio)
	options=list(extra.args=extra.args)
	res$cor=NULL
	res$maxgrad=NULL
	results=c(beta=beta,neg2lnl=-2*res$loglik,AIC=-2*res$loglik+2*res$npar,convergence=convergence)
	results$hessian=res$hes
	results$optim.details=optim.details
	results$options=options
	results$coeflist=res$coeflist
	results$npar=list(npar=res$npar,npar_sdrpt=res$npar_sdrpt,npar_total=res$npar_total)
	results$beta.vcv=beta.vcv
	res=results
#  Restore non-accumulated, non-scaled dm's etc
	res$model_data=model_data.save
#  Assign S3 class values and return
	class(res)=c("crm","admb","mle","mvmscjs")
	return(res)
}



