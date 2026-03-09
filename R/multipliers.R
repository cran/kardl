

#' Calculating multipliers of estimated model
#'
#'
#' By running \code{m<-mplier(kardl_model)}, the object \code{m$mpsi} contains the final multiplier values, \code{m$omega} holds the omega values, and  \code{m$lambda} returns the Lambda values of the model.
#'
#' @param kmodel  The model produced by kardl function. This is the model object from which the dynamic multipliers will be calculated. The function expects an object of class \code{kardl_lm}, which is the standard output class for linear models estimated using the kardl package. If the input object does not inherit from this class, the function will throw an error, ensuring that it operates on compatible model objects.
#' The function extracts necessary information from the model, such as coefficients, lag structure, and variable names, to compute the dynamic multipliers. It calculates the short-run coefficients, Lambda values, and omega values based on the model's parameters and lag structure. The output includes a matrix of dynamic multipliers (mpsi), which can be used for further analysis or visualization.
#' @param horizon  The horizon over which multipliers will be computed. This parameter defines the time frame for the analysis, allowing users to specify how many periods into the future they want to calculate the multipliers for. For example, setting \code{horizon = 40} will compute the multipliers for 40 periods ahead. The function uses this horizon to determine the number of rows in the output matrix of multipliers and to structure the calculations accordingly.
#' @param minProb  The minimum p-value threshold for including coefficients in the multipliers calculation. Coefficients with p-values above this threshold will be set to zero in the multipliers calculation. This parameter allows users to control the inclusion of coefficients based on their statistical significance. Setting a threshold can help focus the analysis on more relevant variables, but it may also exclude potentially important effects if set too stringently. The default value is \code{0}, which means that all coefficients will be included regardless of their p-values.
#'
#' @return
#' A list containing the following elements:
#' \itemize{
#'  \item \strong{mpsi}: A matrix of dynamic multipliers for each
#'  variable across the specified horizon. The first column typically represents the time horizon (e.g., 0 to 40), while the subsequent columns contain the multiplier values for each variable.
#'  \item \strong{omega}: A vector of omega values calculated based on the  model's parameters and lag structure. These values are used in the calculation of the dynamic multipliers and provide insight into the persistence of effects over time.
#'  \item \strong{lambda}: A matrix of Lambda values that represent the short run coefficients and their transformations based on the model's parameters. These values are essential for understanding the immediate effects of changes in the independent variables on the dependent variable.
#'  \item \strong{horizon}: The horizon over which the multipliers were computed, as specified in the function's input. This value indicates the time frame for which the dynamic multipliers are relevant and can be used for interpretation and analysis of the results.
#'  \item \strong{model}: The original model object that was passed to the function. This allows users to reference the model used for the calculations and can be useful for further analysis or validation of the results.
#'  }
#' @export
#' @import stats msm ggplot2
#' @seealso    \code{\link{bootstrap}}
#' @details
#' The \code{mplier} function computes dynamic multipliers based on the coefficients and lag structure of a model estimated using the \code{kardl} package. The function extracts necessary
#' information from the model, such as coefficients, lag structure, and variable names, to compute the dynamic multipliers. It calculates the short-run coefficients, Lambda values, and omega values based on the model's parameters and lag structure. The output includes a matrix of dynamic multipliers (mpsi), which can be used for further analysis or visualization.
#' The dynamic multipliers provide insight into how changes in the independent variables affect the dependent variable over time, allowing for a deeper understanding of the relationships captured by the model. The function also allows users to set a minimum p-value threshold for including coefficients in the calculation, providing flexibility in focusing on statistically significant effects.
#' The function is designed to work specifically with models estimated using the \code{kardl} package, and it ensures that the input model is of the correct class before proceeding with the calculations. The resulting multipliers can be used for various applications, including forecasting, policy analysis, and understanding the dynamic relationships between variables in time series data.
#'
#' \deqn{
#' {{\psi_{h}^{+} = {\sum_{i = 0}^{h}\frac{\partial y_{t + i}}{\partial x_{t}^{+}}}} ; {\psi_{h}^{-} = {\sum_{i = 0}^{h}\frac{\partial y_{t + i}}{\partial x_{t}^{-}}}}}
#' }
#' The above equation represents the cumulative dynamic multipliers for positive and negative changes in an independent variable (x) on the dependent variable (y) over a specified horizon (h). The multipliers are calculated as the sum of the partial derivatives of the dependent variable with respect to the positive and negative changes in the independent variable across the time horizon. These multipliers provide insight into how changes in the independent variable affect the dependent variable over time, allowing for a deeper understanding of the relationships captured by the model.
#'
#'
#' @examples
#'
#' # This example demonstrates how to use the mplier function to calculate dynamic multipliers
#' # from a model estimated using the kardl package. The example includes fitting a model with
#' # the kardl function, calculating the multipliers, and visualizing the results using both
#' # base R plotting and ggplot2.
#'
#'  kardl_model<-kardl(imf_example_data, CPI~ER+PPI+asym(ER)+deterministic(covid)+trend,
#'  mode=c(1,2,3,0))
#'  m<-mplier(kardl_model,40)
#'  head(m$mpsi)
#'  head(m$omega)
#'  head(m$lambda)
#'
#'  # Displaying the summary of the multipliers object
#'  summary(m)
#'
#'  # Visualize the dynamic multipliers
#'
#' plot(m)
#'
#'  # For plotting specific variables, you can specify them in the plot function. For example,
#'  # to plot the multipliers for the variable "ER":
#'
#'  plot(m, variable = "ER")
#'
#'
mplier<-function(kmodel,horizon=80,minProb=0){

  if(! inherits(kmodel, "kardl_lm")){

    stop("Input object must be of class 'kardl_lm'.",call. = FALSE)
  }
   var_names <-names(coefficients(kmodel))
  properLag<-kmodel$lagInfo$OptLag
  h<-horizon +1
  vars<-kmodel$extractedInfo

  depvar<- replace_lag_var(.kardl_Settings_env$LongCoef ,vars$dependentVar,1)
  a_index <- which(var_names ==depvar) # Index the first lagged dependent in the matrix

 my_alpha_lags <- properLag[[vars$dependentVar]]
  my_alpha_index<-which(var_names %in% unlist( lapply(1:my_alpha_lags,function(i){
    replace_lag_var(.kardl_Settings_env$ShortCoef ,vars$dependentVar,i)
    })))



  model_params <- kmodel$coefficients

  if(minProb > 0){
    tempModel<-summary(kmodel)$coefficients
    for (v in 1:nrow(tempModel)) {
      if(tempModel[v,4]>minProb){
        model_params[v]<-0
      }
    }
  }

  my_theta <- model_params[a_index] # alpha is reserved functipn in R

  my_alpha <- model_params[my_alpha_index] # beta1


  p<-length(my_alpha_index) # the length of the short-run lags of the dependent variable

  my_omega<-numeric(p)
  my_omega[1] = 1 + model_params[a_index] + my_alpha[1]
  if (p >= 2) {
    for (i in 2:p) {
      my_omega[i] = my_alpha[i] - my_alpha[i-1]
    }
  }
  my_omega[p+1] = -my_alpha[p]

  ColNum<-length(vars$independentVars)

  ColName<- ColName_F<-c()
  for(x in vars$independentVars){
    posName<-paste0(.kardl_Settings_env$AsymPrefix[1],x,.kardl_Settings_env$AsymSuffix[1]) ### posName<- paste0(x,".POS")
    negName<-paste0(.kardl_Settings_env$AsymPrefix[2],x,.kardl_Settings_env$AsymSuffix[2]) ### negName<- paste0(x,".NEG")
    diff_name<-paste0(x,"_dif")
    ColName_F<- c(ColName_F,posName,negName,diff_name)
    ColName<- c(ColName,posName,negName)
  }


  shortRun_coefs<-matrix(0, nrow = h, ncol = ColNum*2)
  lambda_mtx<-matrix(0, nrow = h, ncol = ColNum*2)
  mpsi_mtx<-matrix(0, nrow = h, ncol = ColNum*2)
  colnames(shortRun_coefs)<-colnames(lambda_mtx)<-colnames(mpsi_mtx)<-ColName

  for (v in 1:ColNum) {
    mtx_varMultiplier<-vars$independentVars[v]
    mtx_POS<-mtx_NEG<-mtx_varMultiplier
    if(mtx_varMultiplier %in% vars$ALvars){
      mtx_POS <-paste0(.kardl_Settings_env$AsymPrefix[1],mtx_varMultiplier,.kardl_Settings_env$AsymSuffix[1]) ### mtx_POS<- paste0(mtx_varMultiplier,".POS")
     mtx_NEG <-paste0(.kardl_Settings_env$AsymPrefix[2],mtx_varMultiplier,.kardl_Settings_env$AsymSuffix[2]) ### mtx_NEG<- paste0(mtx_varMultiplier,".NEG")
      # mtx_POS<-paste0(vars$AsymPrefix[1],mtx_varMultiplier,vars$AsymSuffix[1])
     # mtx_NEG<-paste0(vars$AsymPrefix[2],mtx_varMultiplier,vars$AsymSuffix[2])
     mtx_index_POS<-which(var_names == replace_lag_var(.kardl_Settings_env$LongCoef ,mtx_POS,1))
      mtx_index_NEG<-which(var_names == replace_lag_var(.kardl_Settings_env$LongCoef ,mtx_NEG,1))

     # mtx_index_POS<-which(var_names == paste0("L1.",mtx_POS))
     # mtx_index_NEG<-which(var_names == paste0("L1.",mtx_NEG))
    }else{
      mtx_index_NEG<-mtx_index_POS <- which(var_names == replace_lag_var(.kardl_Settings_env$LongCoef ,mtx_varMultiplier,1))
                                            #paste0("L1.", mtx_varMultiplier ))
    }
    mtx_beta_index_p <- which(var_names %in% unlist( lapply(0:properLag[mtx_POS],function(i){
      replace_lag_var(.kardl_Settings_env$ShortCoef ,mtx_POS,i)
      # paste0("L",i,".d.",mtx_POS)

      })))
    mtx_beta_index_n <- which(var_names %in% unlist( lapply(0:properLag[mtx_NEG],function(i){
      replace_lag_var(.kardl_Settings_env$ShortCoef ,mtx_NEG,i)
      #paste0("L",i,".d.",mtx_NEG)
      })))

    mtx_beta_POS<- model_params[mtx_beta_index_p]
    mtx_beta_NEG<-model_params[mtx_beta_index_n] #beta2 beta3

    tempColNo<-(v-1)*2+1
    shortRun_coefs[1:length(mtx_beta_POS),tempColNo]<-mtx_beta_POS
    shortRun_coefs[1:length(mtx_beta_NEG),tempColNo+1]<-mtx_beta_NEG
    lambda_mtx[1,tempColNo]<-shortRun_coefs[1,tempColNo]
    lambda_mtx[1,tempColNo+1]<-shortRun_coefs[1,tempColNo+1]
    lambda_mtx[2,tempColNo]<- model_params[mtx_index_POS] - shortRun_coefs[1,tempColNo]+shortRun_coefs[2,tempColNo]
    lambda_mtx[2,tempColNo+1]<-model_params[mtx_index_NEG] - shortRun_coefs[1,tempColNo+1]+shortRun_coefs[2,tempColNo+1]
    my_beta_lags_p <-  properLag[[mtx_POS]]
    if(my_beta_lags_p>1){
      for (i in 2:length(mtx_beta_POS)) {
        lambda_mtx[i+1,tempColNo] <- shortRun_coefs[i+1,tempColNo] - shortRun_coefs[i ,tempColNo]
      }
    }
    if(my_beta_lags_p>0){
      lambda_mtx[my_beta_lags_p+2,tempColNo]<- -shortRun_coefs[my_beta_lags_p+1,tempColNo]
    }

    if( mtx_POS==mtx_NEG){
      lambda_mtx[,tempColNo+1]<-lambda_mtx[,tempColNo]
    }else{
      my_beta_lags_n <-  properLag[[mtx_NEG]]
      if(my_beta_lags_n>1){
        for (i in 2:length(mtx_beta_NEG)) {
          lambda_mtx[i+1,tempColNo+1] <- shortRun_coefs[i+1,tempColNo+1] - shortRun_coefs[i,tempColNo+1 ]
        }
      }
      if(my_beta_lags_n>0){
        lambda_mtx[my_beta_lags_n+2,tempColNo+1]<- -shortRun_coefs[my_beta_lags_n+1,tempColNo+1]
      }
    }
  }

  mpsi_mtx[1,]<-lambda_mtx[1,]
  p2<-p+1
  for (v in 2:h) {
    n_say<-ifelse(v>p2,p2,v-1)
    g_say<-ifelse(v-p2<=0,1,v-p2)
    n_omega <- my_omega[1:n_say]
    mpsi_mtx[v,] <-lambda_mtx[v,] + n_omega %*% mpsi_mtx[(v-1):g_say,]

  }
  # even_columns <- seq(2, ncol(mpsi_mtx), by = 2)
  # mpsi_mtx[, even_columns] <- mpsi_mtx[, even_columns] * -1
  mpsi_mtx2 <-apply(mpsi_mtx, 2, cumsum)

  mpsi<-matrix(0, nrow = h, ncol = ColNum*3)
  for (v in 1:ColNum) {
    colNo2<-v*2-1
    colNo3<-v*3-2
    mpsi[,colNo3]<-mpsi_mtx2[,colNo2]
    mpsi[,colNo3+1]<--mpsi_mtx2[,colNo2+1]
    mpsi[,colNo3+2]<-mpsi_mtx2[,colNo2]-mpsi_mtx2[,colNo2+1]
  }

  mpsi<-cbind(0:horizon,mpsi)
  colnames(mpsi)<- c("h",ColName_F)

  structure(
    list(
      call = match.call(),
      mpsi = mpsi,
      omega = my_omega,
      lambda = lambda_mtx,
      horizon = horizon,
      vars=vars
    ),
    class = "kardl_mplier"
  )
}

#' Produce Bootstrap Confidence Intervals for Dynamic Multipliers
#'
#' This function computes bootstrap confidence intervals (CI) for dynamic multipliers
#' of a specified variable in a model estimated using the \code{kardl} package. The bootstrap method
#' generates resampled datasets to estimate the variability of the dynamic multipliers,
#' providing upper and lower bounds for the confidence interval.
#'
#' @param kmodel The model produced by the \code{\link{kardl}} function. This is the model object
#'        from which the dynamic multipliers are calculated.
#' @param horizon An integer specifying the horizon over which dynamic multipliers will be computed.
#'        The horizon defines the time frame for the analysis (e.g., 40 periods).
#' @param replications An integer indicating the number of bootstrap replications to perform.
#'        Higher values increase accuracy but also computational time. Default is \code{100}.
#' @param level A numeric value specifying the confidence level for the intervals (e.g., 95 for 95% confidence).
#'        Default is \code{90}.
#' @param minProb A numeric value specifying the minimum p-value threshold for including coefficients in the bootstrap.
#'       Coefficients with p-values above this threshold will be set to zero in the bootstrap samples. Default is \code{0} (no threshold).
#'       This parameter allows users to control the inclusion of coefficients in the bootstrap process based on their statistical significance. Setting a threshold can help focus the analysis on more relevant variables, but it may also exclude potentially important effects if set too stringently.
#' @seealso
#' \code{\link{mplier}} for calculating dynamic multipliers
#'
#' @return A list containing the following elements:
#' \itemize{
#' \item \strong{mpsi}: A data frame containing the dynamic multiplier estimates along with their upper and lower confidence intervals for each variable and time horizon.
#' \item \strong{level}: The confidence level used for the intervals (e.g 95).
#' \item \strong{horizon}: The horizon over which the multipliers were computed (e.g., 40).
#' \item \strong{vars}: A list of variable information extracted from the model, including dependent variable, independent variables, asymmetric variables, and deterministic terms.
#' \item \strong{replications}: The number of bootstrap replications performed.
#' \item \strong{type}: A character string indicating the type of analysis, in this case "bootstrap".
#' }
#'
#'
#' @details
#'
#' The \code{mpsi} component of the output contains the dynamic multiplier estimates along with their upper
#' and lower confidence intervals. These values are provided for each variable and at each time horizon.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#'
#' library(dplyr)
#'
#'   # Example usage of the bootstrap function
#'
#'
#'  # Fit a model using kardl
#'  kardl_model <- kardl(imf_example_data,
#'                       CPI ~ ER + PPI + asy(ER) +
#'                        det(covid) + trend,
#'                       mode = c(1, 2, 3, 0))
#'
#'  # Perform bootstrap with specific variables for plotting
#'  boot <-
#'    bootstrap(kardl_model,   replications=5)
#'  # The boot object will include all plots for the specified variables
#'  # Displaying the boot object provides an overview of its components
#'  names(boot)
#'
#'
#'  # Inspect the first few rows of the dynamic multiplier estimates
#'   head(boot$mpsi)
#'
#'
#'   summary(boot)
#'
#'  # Retrieve plots generated during the bootstrap process
#'  # Accessing all plots
#'   plot(boot)
#'
#'  # Accessing the plot for a specific variable by its name
#'   plot(boot, variable = "PPI")
#'  plot(boot, variable = "ER")
#'
#'
#'
#'  # Using dplyr
#'  library(dplyr)
#'
#'    imf_example_data %>% kardl( CPI ~ PPI + asym(ER) +trend, maxlag=2) %>%
#'    bootstrap(replications=5) %>% plot(variable = "ER")
#'
#'
bootstrap<-function(kmodel,
                    horizon=80,
                    replications=100,
                    level=95,
                    minProb=0){

  if(! inherits(kmodel, "kardl_lm")){

    stop("Input object must be of class 'kardl_lm'.",call. = FALSE)
  }

  My_mp<-mplier(kmodel,horizon, minProb)
  mpsi<-as.data.frame(My_mp[["mpsi"]])
  vars<-kmodel$extractedInfo

  if(length( vars$ASvars)>0){
    b0 <- coef(kmodel )
    tmin <- kmodel$estInfo$start
    tmax <-OBS <- kmodel$estInfo$end
    OptLag<-kmodel$lagInfo$OptLag
    nardlres<-residuals(kmodel)
    NewDiff<-matrix(NA,nrow = 1,ncol = length(vars$ASvars)+1)
    colnames(NewDiff)<-c("h",  vars$ASvars )
    phi <-  My_mp$omega
    phi.length<-length(phi)

    p<-OptLag[[vars$dependentVar]]+1

    q<-max(OptLag[names(OptLag)!= vars$dependentVar])+2
    # q<-max(OptLag[-1])+2
    theta<-My_mp$lambda[1:q,]


    # excluding symmetric independent vars
    symPosNames<- unlist(lapply(vars$indepASexcluded ,function(x){paste0(.kardl_Settings_env$AsymPrefix[1],x,.kardl_Settings_env$AsymSuffix[1])}) )
    symNegNames<- unlist(lapply(vars$indepASexcluded ,function(x){paste0(.kardl_Settings_env$AsymPrefix[2],x,.kardl_Settings_env$AsymSuffix[2])}) )

     if(!is.null(symPosNames)){
      theta<-theta[, colnames(theta) != symPosNames]
      colnames(theta)[colnames(theta) %in% symNegNames]<-vars$indepASexcluded
    }

    ## names of the indep. variables that matter for multipliers
    inddColnames <- colnames(theta)

    ## build the column list; add deterministic terms only if they exist
    reqCols <- c(vars$dependentVar,
                 inddColnames,
                 if (!is.null(vars$deterministic)) vars$deterministic)

    vars$data <- as.data.frame(vars$data)

    ## --- NEW: make sure every requested column exists --------------------
    missingCols <- setdiff(reqCols, colnames(vars$data))
    if (length(missingCols)) {
      vars$data[ , missingCols] <- 0          # fill with zeros
    }
    ## --------------------------------------------------------------------

    ## subset the data safely
    my_newData <- as.data.frame(vars$data)[ , reqCols, drop = FALSE]


    # Addin dif suffix to asymmetric short-run var names


    # Replace leading NA values with 0 for each column
    for(col in colnames(my_newData)) {
      first_non_na_index <- which(!is.na(my_newData[,col]))[1] # Find the index of first non-NA value
      if (!is.na(first_non_na_index) && first_non_na_index>1) { # Check if non-NA value exists
        my_newData[1:(first_non_na_index - 1), col] <- 0 # Replace leading NA values with 0
      }
    }



    newFormula <- kmodel$argsInfo$formula
    newDepData<- my_newData[,vars$dependentVar] #my_newData[,1]

    for (r in 1:replications) {
        for (t in tmin:tmax) {
          newy1<-sum(phi*newDepData[(t-1):(t-phi.length)]) # calculating dependent
          newy2<- sum(theta*my_newData[(t-0):(t-q+1),inddColnames]) #calculating independent vars
          newy3<-sum(my_newData[t,vars$deterministic]*b0[vars$deterministic]) #calculating deterministic vars
          newy <- newy1+ newy2+newy3+ b0[[1]] # adding constant
          newDepData[t]  <-   newy  +   nardlres[sample(1:NROW(nardlres), 1)] #adding random residual
        }
        if(r==1){
          BSnewdata<-cbind(newDepData,vars$data) # adding to original data
          colnames(BSnewdata)<-c("BStmp",colnames(vars$data))
          tempFormula<-BStmp~x+ff
          newFormula[[2]]<-tempFormula[[2]]
          names(OptLag)[1]<-"BStmp"
          OptLag<-as.double(OptLag)
        }else{
          BSnewdata[,"BStmp"]<-newDepData
        }

        My_model<-kardl(BSnewdata,newFormula,mode=OptLag)
        My_mpNew<- mplier(My_model,horizon ,minProb)
        NewDiff<-rbind(NewDiff, cbind(horizon=0:horizon, val=My_mpNew$mpsi[,paste0( vars$ASvars,"_dif")]))
    }

    NewDiff<-as.data.frame( NewDiff[-1,])
    rownames(NewDiff)<-1:nrow(NewDiff)

    lowerZ <- 0.5 * (1 - level / 100)
    rZ<-replications * lowerZ
    lowerIndex1 <- ceiling(rZ)
    lowerIndex2 <- ifelse(lowerIndex1 - rZ != 0,lowerIndex1,rZ + 1)

    medianZ<-0.5
    rZ<-replications * medianZ
    medianIndex1 <- ceiling(rZ)
    medianIndex2 <- ifelse(medianIndex1 - rZ != 0,medianIndex1,rZ + 1)

    upperZ<-1-lowerZ
    rZ<-replications * upperZ
    upperIndex1 <- ceiling(rZ)
    upperIndex2 <- ifelse(upperIndex1 - rZ != 0,upperIndex1,rZ + 1)

    for (variable in vars$ASvars) {
      lowerCI<-medianCI<-upperCI<-c()
      for (jh in 0:horizon) {
        nardlbs_sub <- NewDiff[NewDiff$h==jh,]
        nardlbs_sub <- nardlbs_sub[order(nardlbs_sub[,variable]), variable]
        # nardlbs_sub <- unlist( nardlbs_sub[2])
        lowerCI[jh+1]<-(nardlbs_sub[lowerIndex1]+nardlbs_sub[lowerIndex2])/2
        medianCI[jh+1]<-(nardlbs_sub[medianIndex1]+nardlbs_sub[medianIndex2])/2
        upperCI[jh+1]<-(nardlbs_sub[upperIndex1]+nardlbs_sub[upperIndex2])/2
      }
      dynNames<-colnames(mpsi)
      mpsi<-cbind(mpsi,upperCI,lowerCI)
      colnames(mpsi)<-c(dynNames,paste0(variable,c("_uCI","_lCI")))

    }
  }


  my_output<-  list(
    mpsi=mpsi,
    level=level,
    horizon=horizon,
    vars=vars,
    replications=replications,
    type="bootstrap"
  )
  class(my_output)<-c("kardl_boot","kardl_mplier")
  my_output

}

mplierggplot<-function(mpsi,vars, varName  ){
  x_p <- paste0(.kardl_Settings_env$AsymPrefix[1], varName, .kardl_Settings_env$AsymSuffix[1]) # paste0(varName,".POS")
  x_n <- paste0(.kardl_Settings_env$AsymPrefix[2], varName, .kardl_Settings_env$AsymSuffix[2]) # paste0(varName,".NEG")

  x_diff<-paste0(varName,"_dif")
  upperCI<-paste0(varName,"_uCI")
  lowerCI<-paste0(varName,"_lCI")

  p<- ggplot(mpsi, aes(x=.data$h)) +
    geom_line(aes(y = .data[[x_p]], color = "Positive")) +
    geom_line(aes(y = .data[[x_n]], color = "Negative"))
    # geom_line(aes(y=!!rlang::sym(x_p),color="Positive")) +
    # geom_line(aes(y=!!rlang::sym(x_n),color="Negative"))


  myBreaks<- c( "Positive", "Negative")
  labels <- c(
    "Positive" = "Positive Change",
    "Negative" = "Negative Change")
  values <- c( # Specify colors for each legend item
    "Positive" =  alpha("blue",1),
    "Negative" =  alpha("red",1) )
  color<-c("blue", "red")
  fill<-c("transparent","transparent")
  linetype <-c(1,1) # c("dashed", "solid","solid","solid") ,
  size<-c(8,8)  # Increases the width of the legend key

  if(upperCI  %in% colnames(mpsi) ){ # !is.null(mpsi[[upperCI]])){

    # p <- p+    geom_line(aes(y=!!rlang::sym(x_diff),color="Difference"),linetype = "dashed") +
    #   geom_ribbon(aes(ymin = !!rlang::sym(lowerCI), ymax = !!rlang::sym(upperCI),color="CI") ,fill=  "blue"  ,   alpha = 0.2 )
     p <- p+    geom_line(aes(y = .data[[x_diff]], color = "Difference"), linetype = "dashed") +
      geom_ribbon(aes(ymin = .data[[lowerCI]], ymax = .data[[upperCI]], color = "CI"), fill = "blue", alpha = 0.2)

    myBreaks<- c("Difference",myBreaks, "CI")
    labels <- c("Difference" = "Asymmetry", labels  , "CI" = "Confidence Interval")
    values <- c("Difference" = alpha("black",1),  values  ,       "CI" =  alpha("blue",0.2) )
    color<-c("black",color,"transparent")
    fill<-c("transparent", fill,"blue")
    linetype <-c(2,linetype,0)
    size<-c( 8,size,8)
  }

  p <- p+     scale_colour_manual(name="",
                                  breaks =myBreaks,  # Include "Confidence Interval" in breaks
                                  labels = labels,  # Include label for "Confidence Interval"
                                  values = values


  ) +  # Use the same color for "Confidence Interval" as the fill color
    #scale_fill_manual(name="",values=  "#84CA72"                   ,labels="CI")+
    scale_x_continuous(limits = c(0,nrow(mpsi)-1),  expand = c(0, 0))+
    labs(x = "time",y="") + ggtitle(paste0("Cumulative effect of ",varName," on ",vars$dependentVar))+
    guides( color = guide_legend(
      override.aes = list(
        color = color,
        fill  = fill,
        linetype =linetype, # c("dashed", "solid","solid","solid") ,

        size=size  # Increases the width of the legend key
      ) )
    )+  theme_bw()  +
    theme( legend.key.width = unit(0, "cm") ,
           legend.position = "bottom"  ,
           legend.margin = margin(t = 0, unit='cm')

    )

  p
}


