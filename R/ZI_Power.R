#' @title ZI_Power
#' @description This function calculates power with Monte-Carlo simulation when using zero-inflated count GLM models.  
#' @param model An R model formula specifying the regression relationship to be used.
#' @param family A character value of either “poisson” or “negbin”, which specifies the family that the count distribution in the zero-inflated GLM model should follow.
#' @param cov_interest A character string for the name of the predictor you are interested in calculating power for.
#' @param data A data frame containing the data set to be used in fitting the initial model.
#' @param nsim An integer value indicating the number of simulation iterations to perform for each different sample size value.
#' @param grid A vector which specifies the sequence of sample size values to perform the calculation over, Default: seq(100, 1000, 100)
#' @param alpha The significance level you would like to use when calculating the proportion of times p-values were less than a certain size, Default: .05
#' @param padj Character string specifying a p-value adjustment from the options available for the stats::p.adjust function.  Default: 0
#' @return A list object with two components: Results and Plot. Results is a data.frame object containing the sample size values used and the corresponding power calculation for each sample size. Plot contains a ggplot object which create a power curve showing the power as a function of sample size.
#' @details For the padj argument, a value of 0 gives no correction and a string of any method in p.adjust() documentation will perform the corresponding adjustment in each power calculation.
#' @references  
#' @examples 
#' \dontrun{
#' Below code corresponds to usage based on the model outlined in our research paper:
#'  ZI_Power(model=Depression~Sex+EOD_total,cov_interest = "EOD_total",family = "poisson",data = dat,nsim = 500,grid = seq(200,4000,100),padj = 0)
#'  ZI_Power(model=Depression~Sex+EOD_total,cov_interest = "EOD_total",family = "poisson",data = dat,nsim = 500,grid = seq(200,4000,100),padj = "BY")
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


ZI_Power <- function(model,family,cov_interest,data,nsim,grid=seq(100,1000,100),alpha=0.05,padj=0){
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
  
  
  system.time(calcs <- llply(as.list(grid),PowerAnalysis,nsim=nsim,cor=padj,.progress = "text"))
    
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
