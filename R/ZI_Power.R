#' @title ZI_Power
#' @description FUNCTION_DESCRIPTION
#' @param model PARAM_DESCRIPTION
#' @param family PARAM_DESCRIPTION
#' @param cov_interest PARAM_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param nsim PARAM_DESCRIPTION
#' @param grid PARAM_DESCRIPTION, Default: seq(100, 1000, 100)
#' @param alpha the specified significance level, Default: .05
#' @param padj specification for an option p-value adjustment using options in stats::p.adjust documentation. Value of 0 gives no correct and a string of any method in p.adjust() documentation can be used.  Default: 0
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @references Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics, 29, 1165–1188. doi:10.1214/aos/1013699998.
#' @references 
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #ZI_Power(model=Depression~Sex+EOD_total,cov_interest = "EOD_total",family = "poisson",data = dat,nsim = 500,grid = seq(200,4000,100),cort = 0)
#'  #ZI_Power(model=Depression~Sex+EOD_total,cov_interest = "EOD_total",family = "poisson",data = dat,nsim = 500,grid = seq(200,4000,100),cort = "BY")
#'  }
#' }
#' @seealso 
#'  \code{\link[pacman]{p_load}}
#'  \code{\link[pscl]{zeroinfl}}
#'  \code{\link[parameters]{p_value}}
#' @rdname ZI_Power2
#' @export 
#' @importFrom pacman p_load
#' @importFrom pscl zeroinfl
#' @importFrom parameters p_value


ZI_Power <- function(model,family,cov_interest,data,nsim,grid=seq(100,1000,100),alpha,padj=0){
  if (!require("pacman")) install.packages("pacman") 
  pacman::p_load(tidyverse,performance,pscl,stringr,stats,ggplot2,plyr)

  ## Extract vector of covariate names
  cov <- attr(terms(formula(model)), which = "term.labels")
  ## Calculate number of covariates
  ncov <- length(cov)
  ## Extract response variable name
  response <- as.character(attr(terms(formula(model)), which = "variables")[[2]])
  
  ## Logic to determine correct index of coefficients 
  ##  in both components of mixture model for any number of covariates
  covi1 <- 2:(ncov+1)
  covi2 <- (3+ncov):(ncov+4)
  
  ## Create vector of x1,x2,...,xn where n is number of covariates
  xs <- paste0("x",covi1-1)
  
  ## Extract model formula
  model <- as.formula(model)  
  
  ## Fit initial model to extract coefficients
  init <- pscl::zeroinfl(model,data,dist=family)
  coefs <- coef(init)
  
  ## Below create string representations of the linear combinations for 
  ##  each coefficients and predictor to define helper functions for relating
  ##  simulated covariates to the response mean through link function
  
  ## Note, this is only necessary to add flexibility for any number of
  ##  covariates
  c1s <- paste0("(coefs[",covi1,"]*",xs,")")
  c1s <- paste0(c1s,collapse = "+")
  
  c2s <- paste0("(coefs[",covi2,"]*",xs,")")
  c2s <- paste0(c2s,collapse = "+")
  
  xc <- paste0(xs,collapse = ",")
  
  lc1 <- paste0("coefs[1]+",c1s,collapse = "")
  lc2 <- paste0("coefs[",2+ncov,"]+",c2s,collapse = "")
  
  ## Actual string representations of the functions creating linear combination
  ##  of predictors and coefficients
  xbcp <- paste0("XBC <- function(",xc,"){",lc1,"}")
  xbbp <- paste0("XBB <- function(",xc,"){",lc2,"}")
  
  ## Parse and evaluate these strings to create functions in environment
  eval(parse(text=xbcp))
  eval(parse(text=xbbp))
  
  ## Define sigmoid function
  logistic <- function(xb){
    exp(xb)/(1+exp(xb))
  }
  
  
  ## Function to simulate all covariates by bootstrapping rows
  Sim_Cov <- function(ss,cov_names){
    
    
    i <- sample(size = ss,x = 1:nrow(data),replace = T)
    return(data[i,c(cov_names)])
    
  } 
  
  
  ## Extract all coefficient names for mixture (ZI) GLM model
  cn <- names(coefs)
  cn <- str_split(cn,"_",simplify = T,n = 2)
  cn <- cn[,2]
  
  ## Identify which coefficients relate to the covariate of interest for power
  pind <- cn %in% cov_interest
  
  
  ## Case when using ZIP model
  if(family=="poisson"){
    
    ## Functions which simulate the individual components necessary to simulate our response,
    ##  which is a mixture of poisson and bernoulli random variables with parameters determined
    ##  by link function, model coefficients, and simulated covariates.
    rpp <- paste0("rp <- function(",xc,"){rpois(length(",xs[1],"),lambda=exp(XBC(",xc,")))","}")
    rbnp <- paste0("rbn <- function(",xc,"){rbinom(length(",xs[1],"),1,prob=1-logistic(XBB(",xc,")))","}")
    
    ## Create functions within environment
    eval(parse(text=rpp))
    eval(parse(text=rbnp))
    
    
    
    
    
    ## Function to perform simulation and model fitting for single iteration
    Sim <- function(ss){
      
      
      ## Simulate covariates
      COV <- Sim_Cov(ss,cov)
      ## String so we can create data.frame of both simulated covariates and response
      COVc <- paste0("COV[,",1:ncov,"]",collapse = ",")
      
      ## Code to simulate response based on simulated covariates
      srp <- paste0("SimResp <- rbn(",COVc,")*rp(",COVc,")",collapse = "")
      eval(parse(text=srp))
      
      ## Create data.frame of all simulated predictors and responses
      sdp <- paste0("SimData <- data.table::data.table(Response=SimResp,",COVc,")",collapse = "")
      eval(parse(text=sdp))
      colnames(SimData) <- c(response,cov)
      
      ## Fit new model on simulated data
      m <- pscl::zeroinfl(model,data=SimData,dist=family)
      
      ## Determine index of p-values we need for the covariate of interest
      pindt <- c(which(pind==TRUE))
      ## Extract p-values for cov_interest coefficient in both components of ZIP model
      p <- c(parameters::p_value(m)$p[pindt[1]],parameters::p_value(m)$p[[pindt[2]]])
      ## Return the p-values in vector
      return(c(Count=p[1],Zero=p[2]))
      
      
    }
    
    
    
  } else {
    
    ## Functions which simulate the individual components necessary to simulate our response,
    ##  which is a mixture of poisson and bernoulli random variables with parameters determined
    ##  by link function, model coefficients, and simulated covariates.
    ## Also, extract estimated theta from initial model
    theta <- init$theta  
    rnbp <- paste0("rnb <- function(",xc,"){rnbinom(length(",xs[1],"),mu=exp(XBC(",xc,")),size=",theta,")}")
    rbnp <- paste0("rbn <- function(",xc,"){rbinom(length(",xs[1],"),1,prob=1-logistic(XBB(",xc,")))","}")
    ## Create functions within environment
    eval(parse(text=rnbp))
    eval(parse(text=rbnp))
    
    
    ## Function to perform simulation and model fitting for single iteration
    ##  All code is same except we simulated response with different distributions
    Sim <- function(ss){
      
      
      COV <- Sim_Cov(ss,cov)
      COVc <- paste0("COV[,",1:ncov,"]",collapse = ",")
      
      srp <- paste0("SimResp <- rbn(",COVc,")*rnb(",COVc,")",collapse = "")
      
      eval(parse(text=srp))
      
      
      sdp <- paste0("SimData <- data.table::data.table(Response=SimResp,",COVc,")",collapse = "")
      eval(parse(text=sdp))
      
      colnames(SimData) <- c(response,cov)
      
      m <- pscl::zeroinfl(model,data=SimData,dist=family)
      
      pindt <- c(which(pind==TRUE))
      p <- c(parameters::p_value(m)$p[pindt[1]],parameters::p_value(m)$p[[pindt[2]]])
      
      return(c(Count=p[1],Zero=p[2]))
      
      
    }
    
  }
  
  
  
  
  ## Function to performed Power Analysis by replicating operation many times for any value
  ##  of sample size
  PowerAnalysis <- function(ss,nsim,cor=0){
    
    ## Replicate single iteration nsim times for given sample size
    sims <- replicate(nsim,expr = Sim(ss=ss),simplify = F)
    ## Combine all results into a matrix
    ps <- do.call(rbind,sims)
    
    ## Based on p-value adjustment/correction setting,
    ##  we take the proportion of times both unadjusted/adjusted p-value
    ##  for the coefficients in each model component are less than alpha
    if(cor==0){
      return(mean(ps[,1]<alpha|ps[,2]<alpha))
    } else if(is.character(cor)){
      ps <- apply(ps,MARGIN = 1,FUN=p.adjust,method=cor) %>% t()
      return(mean(ps[,1]<alpha|ps[,2]<alpha))
    }
    
    return(mean(ps[,1]<alpha | ps[,2]<alpha))  
    
    
  }
  
  
  system.time(calcs <- lapply(X = grid,FUN = PowerAnalysis,nsim=nsim,cor=padj))
  
  out <- data.frame(SampleSize=grid,Power=unlist(calcs))
  
  
  plot <- ggplot(data=out)+
    geom_smooth(aes(x=SampleSize,y=Power),color="black",se = FALSE,method = "loess")+
    theme_bw()+
    scale_x_continuous(breaks = grid)+
    scale_y_continuous(breaks = seq(0,1,.05))+
    labs(x="Sample Size",y="Power",title="Power Curve")+
    theme(axis.text.x = element_text(angle=60,vjust=.5),axis.title  = element_text(hjust = .5),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust=0.5))
  return(list(Results=out,Plot=plot))
  
  
  
  
  
}
