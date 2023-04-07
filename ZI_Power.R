## Not sure any of these are necessary since packages are referenced with :: in function
library(pscl)
library(tidyverse)
library(performance)
library(dplyr)
library(SIBERG)
library(MASS)
library(lmtest)


# dat <- read.csv("~/Downloads/DataSubset_102922.csv")
# dat <- dat %>% mutate(Sex=Sex-1)


ZI_Power2 <- function(model,family,cov_interest,data,nsim,grid=seq(100,1000,100),cort=0){
  require(tidyverse)
  require(performance)
  
  # model <- Depression ~ EOD_total + Sex
  # data <- dat
  # cov_interest <- "EOD_total"
  # nsim <- 10
  # grid <- seq(200,500,100)
  # cort <- 0
  # family <- "poisson"
  
    ## Determine the response, number of covariates, and name of covariates
  
  
  
  cov <- attr(terms(formula(model)), which = "term.labels")
  ncov <- length(cov)
  
  response <- as.character(attr(terms(formula(model)), which = "variables")[[2]])
  ## Logic to determine correct index of coefficients 
  ##  in both components of mixture model for any number of covariates
  covi1 <- 2:(ncov+1)
  covi2 <- (3+ncov):(ncov+4)
  
  xs <- paste0("x",covi1-1)
  
  model <- as.formula(model)  
  
  init <- pscl::zeroinfl(model,data,dist=family)
  
  coefs <- coef(init)
  
  ## Below are all for purpose of allowing any given number of covariate values
  c1s <- paste0("(coefs[",covi1,"]*",xs,")")
  c1s <- paste0(c1s,collapse = "+")
  
  c2s <- paste0("(coefs[",covi2,"]*",xs,")")
  c2s <- paste0(c2s,collapse = "+")
  
  xc <- paste0(xs,collapse = ",")
  
  lc1 <- paste0("coefs[1]+",c1s,collapse = "")
  lc2 <- paste0("coefs[",2+ncov,"]+",c2s,collapse = "")
  
  ## Functions create linear combinations of covariates and coefficients 
  ## necessary to simulate response
  xbcp <- paste0("XBC <- function(",xc,"){",lc1,"}")
  xbbp <- paste0("XBB <- function(",xc,"){",lc2,"}")
  
  eval(parse(text=xbcp))
  eval(parse(text=xbbp))
  
  logistic <- function(xb){
    exp(xb)/(1+exp(xb))
  }
  
  
  ## Function to simulate a single covariate
  Sim_Cov <- function(ss,cov_name){
    
    i <- sample(size = ss,x = 1:nrow(data),replace = T)
    return(data[i,cov_name])
    
  }
  
  ## Function to simulate every covariate and return a data.frame
  Sim_All_Cov <- function(ss){
  
      
    sims <- list()
    for(i in 1:length(cov)){
      sims[[i]] <- Sim_Cov(ss,cov[i])
      
    }
    return(do.call(cbind,sims))
  }
  
  
  cn <- names(coefs)
  cn <- str_split(cn,"_",simplify = T,n = 2)
  cn <- cn[,2]
  
  pind <- cn %in% cov_interest
  
  if(family=="poisson"){
    
    
    rpp <- paste0("rp <- function(",xc,"){rpois(length(",xs[1],"),lambda=exp(XBC(",xc,")))","}")
    rbnp <- paste0("rbn <- function(",xc,"){rbinom(length(",xs[1],"),1,prob=1-logistic(XBB(",xc,")))","}")
    
    eval(parse(text=rpp))
    eval(parse(text=rbnp))
    
    
    
    
    
    ## Function to perform simulation and model fitting for single iteration
    Sim <- function(ss){
      
      
      ## Simulate covariates
      COV <- Sim_All_Cov(ss)
      COVc <- paste0("COV[,",1:ncov,"]",collapse = ",")
      
      
      srp <- paste0("SimResp <- rbn(",COVc,")*rp(",COVc,")",collapse = "")
      
      eval(parse(text=srp))
      
      
      sdp <- paste0("SimData <- data.table::data.table(Response=SimResp,",COVc,")",collapse = "")
      eval(parse(text=sdp))
      
      colnames(SimData) <- c(response,cov)
      
      m <- pscl::zeroinfl(model,data=SimData,dist=family)
      
      
      
      
      pindt <- c(which(pind==TRUE))
      p <- c(parameters::p_value(m)$p[pindt[1]],parameters::p_value(m)$p[[pindt[2]]])
      
      return(c(Count=p[1],Zero=p[2]))
      
      
    }
    
    
    
  } else {
    
    theta <- init$theta  
    rnbp <- paste0("rnb <- function(",xc,"){rnbinom(length(",xs[1],"),mu=exp(XBC(",xc,")),size=",theta,")}")
    rbnp <- paste0("rbn <- function(",xc,"){rbinom(length(",xs[1],"),1,prob=1-logistic(XBB(",xc,")))","}")
    
    eval(parse(text=rnbp))
    eval(parse(text=rbnp))
    
    
    
    Sim <- function(ss){
      
      
      
      ## Simulate covariates
      COV <- Sim_All_Cov(ss)
      
      
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
  
  
  
  
  
  PowerAnalysis <- function(ss,nsim,cor=0){
    
    sims <- replicate(nsim,expr = Sim(ss=ss),simplify = F)
    ps <- do.call(rbind,sims)
    
    if(cor==0){
      return(mean(ps[,1]<.05|ps[,2]<.05))
    } else if(cor==1){
      ps <- apply(ps,MARGIN = 1,FUN=p.adjust,method="bonferroni") %>% t()
      return(mean(ps[,1]<.05|ps[,2]<.05))
    } else if(cor==2){
      ps <- apply(ps,MARGIN = 1,FUN=p.adjust,method="BY") %>% t()
      return(mean(ps[,1]<.05|ps[,2]<.05))
    }
    
    return(mean(ps[,1]<.05 | ps[,2]<.05))  
    
    
  }
  
  
  system.time(calcs <- lapply(X = grid,FUN = PowerAnalysis,nsim=nsim,cor=cort))
  
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

## Note we changed theta from 1/theta, need to confirm this relation with Dr. O'Connell but I think this is correct now
##  ZINB should give almost same results as ZIP since ZINB is reducing to ZIP

# Original 

o1 <- ZI_Power(model="Depression~Sex+EOD_total",cov_interest = "EOD_total",family = "poisson",data = dat,nsim = 500,grid = seq(200,4000,100),cort = 0)

o2 <- ZI_Power2(model=Depression~Sex+EOD_total,cov_interest = "EOD_total",family = "poisson",data = dat,nsim = 500,grid = seq(200,4000,100),cort = 0)

o3 <- ZI_Power2(model=Depression~Sex+EOD_total,cov_interest = "EOD_total",family = "negbin",data = dat,nsim = 100,grid = seq(200,2500,500),cort = 0)

