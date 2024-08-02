# Document required functions from R packages
#' @importFrom stats as.formula median predict
#' @importFrom survival coxph survfit Surv
#' @importFrom rpart rpart
#' @importFrom randomForest randomForest importance


GfuncSurvivalTree=function(obs,delta,dtype,xx)
{

  num = length(obs)
  nu = num

  # Fitting the Cox model
  # Creating the data frame
  # Ensure 'xx' is a data frame
  xx <- if(is.data.frame(xx)) xx else as.data.frame(xx)

  # Creating the data frame for the Cox model, adjusting for censoring
  data.used <- cbind(obs = obs, delta.cens = 1 - delta, xx)

  # Ensure the data frame for modeling includes obs and delta, and uses correct column names
  names(data.used)[-c(1, 2)] <- paste0("xx", 1:ncol(xx))

  # Fitting the Survival Tree
  surv.tree = rpart(Surv(obs,delta.cens)~., data = data.used, minbucket = 30)

  # Getting the Survival Curves.
  pred.surv.tree <- predict(surv.tree, proximity = FALSE)

  # Finding the terminal nodes
  sett=unique(surv.tree[['where']])
  nset=length(sett)

  cens.est = matrix(0, nrow = num, ncol = num)
  obs.used = rep(NA, num)
  delta.used = rep(NA, num)

  for (i in 1:nset){
    # Finding the subset corresponding the ith node of the tree
    subset=(1:nu)[surv.tree[['where']]==sett[i]]
    nlen=length(subset)
    # Observed time within the node
    sobs=obs[subset]
    # Failure indicators within each node.
    sdelta=delta[subset]

    # Doing truncation within that subset
    # Changing the dataset to account for truncation.
    aa=datach(sobs,sdelta,dtype = dtype)
    # Observed time after truncation
    sobs=aa[[1]]
    # Failure indicator after truncation.
    sdelta = aa[[2]]

    obs.used[subset] = sobs
    delta.used[subset] = sdelta

    # Calculating the KM estimator censoring curve within a node
    # Calculating the jumps in the KM estimator
    hazC=mapply(function(xx,dd){dd/sum((sobs>=xx))},xx=sobs,dd=1-sdelta)
    surC_km=mapply(function(xx){prod(1-hazC[sobs<=xx])},xx=obs)
    cens.est[subset, ] = matrix(surC_km,nrow=length(subset),ncol=length(surC_km),byrow=TRUE)
  }

  return(list(cens.est, obs.used, delta.used, surv.tree[['where']]))
}



externalRegularTreeKM = function(obs,delta,x_list, mtype,dtype)
{
  n = length(obs)
  nu = length(obs)


  # Calculating the conditional expectation
  xx <- if(is.data.frame(x_list)) x_list else as.data.frame(x_list)
  m1 <- mfunc(obs, delta, xx, mtype)


  # Calculating the conditional censoring distribution.
  tem=GfuncSurvivalTree(obs,delta,dtype,xx)
  # Calculating the censoring distribution
  surC_rf=tem[[1]]
  # Observed event times for adjusted for truncation
  obs=tem[[2]]
  # Failure indicator adjusted for truncation
  delta=tem[[3]]
  # Finding which observations fall in which terminal node
  term.node = tem[[4]]

  # Calculating a0, a1, b0, b1, c0, c1
  a0=delta/diag(surC_rf)
  a1=a0*log(obs)

  b0=(1-delta)/diag(surC_rf)
  b1=b0 * diag(m1)

  c0 <- rep(NA, n)
  c1 <- rep(NA, n)

  # Creating the ordered data
  ord.used = order(obs)
  obs.order = obs[ord.used]
  delta.order = delta[ord.used]

  # Finding the terminal nodes
  sett=unique(term.node)
  nset=length(sett)

  for (i in 1:nset){
    # Finding the subset corresponding the ith node of the tree
    subset=(1:nu)[term.node==sett[i]]
    nlen=length(subset)
    # Observed time within the node
    sobs=obs[subset]
    # Failure indicators within each node.
    sdelta=delta[subset]


    kk=mapply(function(tt){sum((tt<=sobs))},tt=sobs)
    c0[subset]=mapply(function(tt){sum(b0[subset]*(sobs<=tt)/kk)},tt=sobs)
    c1[subset]=mapply(function(tt,i){sum(b0[subset]*(sobs<=tt)*m1[subset,i]/kk)},tt=sobs,i=1:nlen)
  }

  parms = c(a0,a1,b0,b1,c0,c1,obs,delta,diag(m1),nu)

  return(parms)
}


BJexternal=function(obs,delta,x_list, mtype,dtype)
{
  n = length(obs)
  nu = n
  # Creating the new T(t) dataset
  # Calculating the conditional expectation
  m1 <- mfunc(obs, delta, x_list, mtype)
  a1 = delta *  log(obs) + (1 - delta) * diag(m1)

  parms = a1

  return(parms)
}


tree_km=function(obs,delta,xx)
{
  n <- length(obs)
  nu <- length(obs)
  xx <- as.matrix(xx)
  # Fitting an exponential tree.
  fitEXPt <- rpart(Surv(obs,delta)~xx,method="exp", cp = 0)
  nl=length(fitEXPt[['cptable']][,1])
  nnum=(1:nl)[fitEXPt[['cptable']][,4]==min(fitEXPt[['cptable']][,4])]
  CPEXP=mean(fitEXPt[['cptable']][ ,1][c(nnum-1,nnum)])

  # Pruning the tree
  fitEXP=rpart(Surv(obs,delta)~xx,method="exp",cp=CPEXP)

  # Finding the terminal nodes
  sett=unique(fitEXP[['where']])
  nset=length(sett)
  m1=matrix(0,nu, nu)

  for (i in 1:nset)
  {
    # Finding the subset corresponding the ith node of the tree
    subset=(1:nu)[fitEXP[['where']]==sett[i]]
    nlen=length(subset)
    # Observed time within the node
    sobs=obs[subset]
    # Failure indicators within each node.
    sdelta=delta[subset]

    # Calculating the KM estimator survival curve within a node
    a <- survfit(Surv(sobs, sdelta)~1)
    time.used <- a[['time']]
    surv.used <- a[['surv']]

    # Calculating the jumps in the KM estimator
    surv.diff <- c(1, a[['surv']][-length(a[['surv']])]) - a[['surv']]

    for (j in 1:nu)
    {
      if(delta[j]==FALSE){
        # Calculating the conditional expectation as the log(obs) at jump points
        if(obs[j] < max(sobs[sdelta == 1])){
          m1[j,subset] <- sum(log(time.used[time.used > obs[j]]) * surv.diff[time.used > obs[j]])/sum(surv.diff[time.used > obs[j]])
        }

        # Taking care of problem observations.
        if (obs[j] > max(sobs[sdelta == 1])){
          m1[j,subset]=log(obs[j])}
      }
    }
  }
  # Returns a matrix of m_{1i}(tilde T_j) for i,j computed only if Delta_j = 0
  return(m1)
}


coxxRegular = function(obs, delta, x_list) {
  # Check if x_list is properly named
  if (is.null(names(x_list)) || any(names(x_list) == "")) {
    stop("Covariate names are missing, cannot construct the formula.")
  }


  # Convert x_list to a data frame and add obs and delta
  data.used <- data.frame(obs = obs, delta = delta, x_list)

  # Dynamically construct the Cox formula
  # The Cox model formula is dynamically constructed using the names in x_list.
  # This makes the formula construction flexible to changes in the number of covariates or their names.
  covariate_names <- names(x_list)
  formula_rhs <- paste(covariate_names, collapse = " + ")
  cox_formula <- as.formula(paste("Surv(obs, delta) ~", formula_rhs))

  # Fit the Cox model
  cox.mod <- coxph(cox_formula, data = data.used)


  # Getting the Survival Curves.
  cox.surv <- survfit(cox.mod, newdata= data.used)
  nu <- length(obs)

  # Calculation of the survivor function
  m1=matrix(0,nu, nu)
  time.used <- matrix(0, ncol = length(unique(obs)), nrow = nu)
  surv.used <- matrix(0, ncol = length(unique(obs)), nrow = nu)
  surv.diff <- matrix(0, ncol = length(unique(obs)), nrow = nu)
  for(i in 1:nu){
    # Getting the properties of the survival cox function for observation i.
    time.used[i, ] <- cox.surv[i][['time']]
    surv.used[i, ] <- cox.surv[i][['surv']]

    # Calculating the jumps in the Cox model survival curve estimator
    surv.diff[i, ] <- c(1, cox.surv[i][['surv']][-length(cox.surv[i][['surv']])]) - cox.surv[i][['surv']]
  }

  for (j in 1:nu)
  {
    if(delta[j]==FALSE){

      for(i in 1:nu)
      {
        if(obs[j]<=obs[i]){

          if(max(obs[delta==1]) > obs[j]){
            # Calculating the conditional expectation
            m1[j,i]= sum(log(time.used[i, ][time.used[i, ] > obs[j]]) * surv.diff[i, ][time.used[i, ] > obs[j]])/sum(surv.diff[i, ][time.used[i, ] > obs[j]])
          }
        }
      }
      if (max(obs[delta==1]) <= obs[j]){
        m1[j,]=log(obs[j])}

    }
  }

  return(m1)
}


mfunc=function(obs,delta,x_list, mtype)
{
  num=length(obs)
  xx <- as.matrix(data.frame(x_list))

  if (mtype=="tree_km")
  {
    m1=tree_km(obs,delta,xx)
  }

  if (mtype=="cox")
  {

    m1=coxxRegular(obs,delta, x_list) # Currently this function works
  }

  if (mtype=="rand.for")
  {
    m1=randomForest(obs,delta, xx)
  }
  return(m1)
}


datach=function(obs,delta,dtype)
{
  nu = length(obs)
  if(dtype=="b0")
  {
    delta[obs==max(obs)]=TRUE
  }



  if(dtype=="b5")
  {
    delta[order(obs)][floor(nu-0.05*nu):nu]=TRUE
  }


  if(dtype=="b10")
  {
    delta[order(obs)][floor(nu-0.10*nu):nu]=TRUE
  }


  if(dtype=="b15")
  {
    delta[order(obs)][floor(nu-0.15*nu):nu]=TRUE
  }

  if(dtype=="a5")
  {
    delta[order(obs)][floor(nu-0.05*nu):nu]=T
    obs[order(obs)][floor(nu-0.05*nu):nu]=obs[order(obs)][floor(nu-0.05*nu)]
  }
  if(dtype=="a10")
  {
    delta[order(obs)][floor(nu-0.10*nu):nu]=T
    obs[order(obs)][floor(nu-0.10*nu):nu]=obs[order(obs)][floor(nu-0.10*nu)]
  }


  if(dtype=="a15")
  {
    delta[order(obs)][floor(nu-0.15*nu):nu]=T
    obs[order(obs)][floor(nu-0.15*nu):nu]=obs[order(obs)][floor(nu-0.15*nu)]
  }

  if(dtype=="none")
  {
    obs=obs
    delta=delta
  }
  return(list(obs,delta))
}


simL2NP <- function(obs,delta,x_list, n.tree, tau1, data.test, mtype, dtype, mtry){

  p = length(x_list)
  n <- length(obs)
  num = n
  mtry = mtry

  # surv.est[i,j] = P(T >tau1|W_i, tree j)
  surv.est.sim <- matrix(NA, nrow = nrow(data.test), ncol = n.tree)


  parms <- externalRegularTreeKM(obs = obs, delta = delta, x_list = x_list, mtype = mtype, dtype = dtype)


  a1 <- parms[(num+1):(2*num)]
  b1 <- parms[(3*num+1):(4*num)]
  c1 <- parms[(5*num+1):(6*num)]
  y.imp.all <- a1 + b1 - c1

  # Imputing the observations
  # Computing the parameter vector for the simualated martingale
  bs <- sample(1:num, size = num, replace = TRUE)

  # Use 'covariates' x_list directly if it's a data frame, otherwise create one from the list
  xx <- data.frame(x_list)

  # Fitting an external tree
  exttree.one = randomForest(xx, y.imp.all, n.tree=n.tree, mtry=mtry, importance=TRUE)

  importance.results = importance(exttree.one, type=1)

  # This for loop fits one tree in the forest at each iteration
  for(j in 1:n.tree){

    # Imputing the observations
    # Computing the parameter vector for the simualated martingale
    bs <- sample(1:num, size = num, replace = TRUE)

    # Use 'covariates' x_list directly if it's a data frame, otherwise create one from the list
    xx <- data.frame(x_list)

    # Fitting a single tree
    tree.one = randomForest(xx[bs, ], y.imp.all[bs],wt = rep(1,n), ntree=1, mtry=mtry, sampsize=num, replace=F, nodesize=5)

    # Predicting the value of the observations in the test set in order to know in what
    # terminal node they fall.
    predict.test <- predict(tree.one, newdata = data.test)
    predict.train <- predict(tree.one, newdata = xx)

    # Computing the Kaplan-Meier estimator for the simulated martinagel imputed tree
    for(i in 1:length(unique(predict.train))){
      # Finding observations in the bootstrap sample that fall in that node
      obs.term.node <- which(predict.train == unique(predict.train)[i])

      # Finding the observations falling in that terminal node
      delt <- delta[obs.term.node]
      obst <- obs[obs.term.node]

      # Finding the observations in the test set that fall in that node
      obs.orig <- which(predict.test == unique(predict.train)[i])
      # Calculating the probability P(T > t|W_test) using the Kaplan Meier estimator
      #for that bootstrap tree.

      # Finding P(T >tau1|W) for that terminal node using the Kaplan Meier product limit estimator on the dataset
      hazT=mapply(function(xx,dd){dd/sum(obst>=xx)},dd=delt,xx=obst)
      surv.est.sim[obs.orig, j] = prod(1-hazT[obst <= tau1])
    } #end i loop
  } # End j loop

  # Predicting P(T > tau1|W) by averaging over the KM estimator predictions
  pred.surv.sim <- apply(surv.est.sim, 1, mean)

  # Prepare to return both results, each wrapped in their own list
  results <- list(
    predictions = pred.surv.sim,
    importance = importance.results
  )

  return(results)
}

BJsimL2NP <- function(obs,delta,x_list, n.tree, tau1, data.test, mtype, dtype, mtry){

  p = length(x_list)
  n <- length(obs)
  nu <- n
  num = n
  mtry = mtry

  # surv.est[i,j] = P(T >tau1|W_i, tree j)
  surv.est.sim <- matrix(1,nrow = nrow(data.test), ncol = n.tree)

  parms = BJexternal(obs = obs, delta = delta, x_list = x_list, mtype = mtype, dtype = dtype)

  y.imp.all <- parms

  # Imputing the observations
  # Computing the parameter vector for the simualated martingale
  bs <- sample(1:num, size = num, replace = TRUE)

  # Creating the imputed dataset
  xx <- data.frame(x_list)

  # Fitting an external tree
  exttree.one = randomForest(xx, y.imp.all, n.tree=n.tree, mtry=mtry, importance=TRUE)

  importance.results = importance(exttree.one, type=1)

  # This for loop fits one tree in the forest at each iteration
  for(j in 1:n.tree){

    # Imputing the observations
    # Computing the parameter vector for the simualated martingale
    bs <- sample(1:num, size = num, replace = TRUE)

    # Creating the imputed dataset
    xx <- data.frame(x_list)

    # Fitting a single tree
    tree.one = randomForest(xx[bs, ], y.imp.all[bs], wt = rep(1, n), ntree=1, mtry=mtry, sampsize=num, replace=F, nodesize=5)

    # Predicting the value of the observations in the test set in order to know in what
    # terminal node they fall.
    predict.test <- predict(tree.one, newdata = data.test)
    predict.train <- predict(tree.one, newdata = xx)

    # Computing the Kaplan-Meier estimator for the simulated martinagel imputed tree
    for(i in 1:length(unique(predict.train))){
      # Finding observations in the bootstrap sample that fall in that node
      obs.term.node <- which(predict.train == unique(predict.train)[i])

      # Finding the observations falling in that terminal node
      delt <- delta[obs.term.node]
      obst <- obs[obs.term.node]

      # Finding the observations in the test set that fall in that node
      obs.orig <- which(predict.test == unique(predict.train)[i])
      # Calculating the probability P(T > t|W_test) using the Kaplan Meier estimator
      #for that bootstrap tree.

      # Finding P(T >tau1|W) for that terminal node using the Kaplan Meier product limit estimator on the dataset
      hazT=mapply(function(xx,dd){dd/sum(obst>=xx)},dd=delt,xx=obst)
      surv.est.sim[obs.orig, j] = prod(1-hazT[obst <= tau1])
    } #end i loop
  } # End j loop

  # Predicting P(T > tau1|W) by averaging over the KM estimator predictions
  pred.surv.sim <- apply(surv.est.sim, 1, mean)

  # Prepare to return both results, each wrapped in their own list
  results <- list(
    predictions = pred.surv.sim,
    importance = importance.results
  )

  return(results)

}

one.sim <- function(obs, delta, xx, xx.test, n.tree, tau1, type, mtry) {
  tryCatch({
    # Convert training and testing data to data frames
    data.train <- as.data.frame(xx)
    #names(data.train) <- paste0("x", 1:ncol(xx))

    data.test <- as.data.frame(xx.test)
    #names(data.test) <- paste0("x", 1:ncol(xx.test))

    # Create a list of training data covariates to pass as x_list
    x_list <- lapply(data.train, function(column) column)

    # Prepare the rest of the arguments as before
    additional_args <- list(
      obs = obs,
      delta = delta,
      n.tree = n.tree,
      tau1 = tau1,
      data.test = data.test,
      mtype = "cox",
      dtype = "b10",
      mtry = mtry
    )

    # Combine the x_list with the additional arguments
    args_list <- c(list(x_list = x_list), additional_args)

    # Call appropriate function based on the type parameter
    if (type == "drl") {
      pred.new.imp <- do.call(simL2NP, args_list)
    } else if (type == "bjl") {
      pred.new.imp <- do.call(BJsimL2NP, args_list)
    } else {
      stop("Invalid type specified")
    }

    return(pred.new.imp)
  }, error = function(e) {
    cat("Error in one.sim: ", e$message, "\n")
    return(NULL)  # Return NULL
  })
}


#' @title Survival analysis using censoring unbiased random forests
#'
#' @description This function calculates survival predictions (e.g., one year survival probability) and variable importance measures using censoring unbiased random forests for both Buckley-James and doubly robust loss functions.
#'
#' @param formula an object of class "formula"
#'
#' @param train_data A dataframe containing the training dataset.
#'
#' @param test_data A dataframe containing the test dataset for which survival predictions will be made for.
#'
#' @param n.tree Number of trees used in the forest
#
#' @param type The type of loss function to be used in the random forest algorithm. Type = "bjl" uses the Buckley-James loss function and is the default and type= "drl" uses the doubly robust loss function.
#'
#' @param mtry Number of variables randomly sampled as candidates at each split.
#'
#' @param tau1 At what time-point survival predictions should be calculated (e.g, if time is measured in years and tau=1 the algorithm would predict one year survival probabilities).
#'
#' @return A list containing survival predictions on the test set and variable importance measures
#'
#' @examples
#' train_data_path <- system.file("extdata", "train_data.csv", package = "CURT")
#' test_data_path <- system.file("extdata", "test_data.csv", package = "CURT")
#' train_data <- read.csv(train_data_path)
#' test_data <- read.csv(test_data_path)
#' result <- curt(Surv(obs, delta) ~ . -x2-x3, train_data, test_data, n.tree=500, tau1=0.6, type="drl")
#' @export
#'
curt <- function(formula, train_data, test_data, n.tree = 1000, tau1 = NULL, type='drl', mtry=NULL) {
  tryCatch({
    # Extracting the survival object names from the formula
    surv_obj <- as.character(formula[[2]])
    obs_name <- surv_obj[2]
    delta_name <- surv_obj[3]

    # Retrieving survival data from the training dataset
    obs <- train_data[[obs_name]]
    delta <- train_data[[delta_name]]

    # Set default tau1 to the median of obs if not provided
    if (is.null(tau1)) {
      tau1 <- median(obs, na.rm = TRUE)
    }

    # Extracting all covariates
    all_covariates <- names(train_data)[!(names(train_data) %in% c(obs_name, delta_name))]
    formula_str <- deparse(formula)

    # Identifying explicit inclusions and exclusions
    if (grepl("~\\s*\\.", formula_str)) {
      inclusions <- all_covariates
    } else {
      inclusions <- regmatches(formula_str, gregexpr("(?<=~|\\+)[^\\+\\-]+", formula_str, perl = TRUE))[[1]]
    }
    exclusions <- regmatches(formula_str, gregexpr("(?<=\\-)[^\\+\\-]+", formula_str, perl = TRUE))[[1]]
    inclusions <- trimws(inclusions)
    exclusions <- trimws(exclusions)

    # Validate exclusions and inclusions
    invalid_exclusions <- exclusions[!exclusions %in% all_covariates]
    invalid_inclusions <- inclusions[!inclusions %in% all_covariates]
    if (length(invalid_exclusions) > 0 || length(invalid_inclusions) > 0) {
      stop(paste("Invalid covariate modifications:", toString(c(invalid_exclusions, invalid_inclusions))))
    }

    # Determining final covariates based on specified inclusions or exclusions
    if (length(inclusions) > 0) {
      valid_covariates <- inclusions
    } else {
      valid_covariates <- setdiff(all_covariates, exclusions)
    }
    valid_covariates <- setdiff(valid_covariates, exclusions)  # Ensure exclusions are applied

    # Selecting valid covariates and forming the covariate matrix
    xx <- as.matrix(train_data[valid_covariates])
    xx.test <- as.matrix(test_data[valid_covariates])

    # Calculate the default mtry value if not provided by the user
    if (is.null(mtry)) {
      p <- ncol(xx)  # Number of predictor variables
      mtry <- ceiling(sqrt(p))
    }

    # Running the survival analysis prediction model
    result <- one.sim(
      obs = obs,
      delta = delta,
      xx = xx,
      xx.test = xx.test,
      n.tree = n.tree,
      tau1 = tau1,
      type = type,
      mtry = mtry
    )

    return(result)
  }, error = function(e) {
    cat("Error in curt: ", e$message, "\n")
    return(NULL)  # Return NULL in case of an error
  })
}
