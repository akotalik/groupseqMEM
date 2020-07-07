#R code used to yield the simulation results in the manuscript "A group-sequential randomized trial design utilizing supplemental trial data"
#Ales Kotalik

###########################################################################################
################################### Binary Outcome ##########################################
###########################################################################################

library(parallel)
library(gsDesign)



cores=8

#function that draws from the mixture for beta
rmixture = function(wts, betas){
  ##random generation of the indices
  id = sample(1:length(wts),prob=wts,size=nrow(betas),replace=TRUE)  
  id = cbind(1:nrow(betas),id)
  betas[id]
}


doSim <- function(seed, max, p0, looks, p, beta_t, a=NULL, b=NULL, bdryMEM=NULL, pi_omega=NULL, constrained=NULL, method) {
  ls = list()  
  for(k in 1:looks) {
    #gather 1/looks fraction of data, run the method, either stop early or loop again
    n=max*(1/looks)
    n_sources = length(n)
    
    #simulate data:
    #primary
    trt<- rep(0, n[1])
    trt[sample(1:n[1],n[1]*p[1])] <- 1
    Y <- rbinom(n[1],1,p0[1] + beta_t[1]*trt)
    primary <- data.frame(Y,trt)
    primary$df <- 1
    primary$secondary <- 0
    
    #secondary
    secondary <- NULL
    for(i in 2:n_sources) {
      trt<- rep(0, n[i])
      trt[sample(1:n[i],n[i]*p[i])] <- 1
      Y <- rbinom(n[i],1,p0[i] + beta_t[i]*trt)
      df <- i
      secondary <- rbind(secondary, data.frame(Y,trt,df))
    }
    secondary$secondary <- 1
    
    #data frame together
    dt <- data.frame(rbind(primary, secondary))
    
    names(dt) <- c("Y", "trt", "df", "secondary")
    ls[[k]] <- dt
    combined <- data.frame(do.call(rbind, ls))
    
    if(method == "MEM") {
      samples <- NULL
      for(i in 0:1) {
        y_1t <- sum(combined$Y[combined$df==1 & combined$trt==i])
        y_2t <- sum(combined$Y[combined$df==2 & combined$trt==i])
        n_1t <- sum(combined$df==1 & combined$trt==i)
        n_2t <- sum(combined$df==2 & combined$trt==i)
        
        #model 1: no borrowing
        betas <- NULL
        marg <- NULL
        betas <- cbind(betas, rbeta(100000, y_1t+a, b+n_1t-y_1t))
        marg <- c(marg, (1/(beta(a, b))^2)*beta(y_1t+a, n_1t-y_1t+b)*beta(y_2t+a, n_2t-y_2t+b))
        #model 2: borrowing
        betas <- cbind(betas, rbeta(100000, y_1t+y_2t+a, b+n_1t+n_2t-y_1t-y_2t))
        marg <- c(marg, (1/(beta(a, b)))*beta(y_1t+y_2t+a, b+n_1t+n_2t-y_1t-y_2t))
        
        #weights
        wts <- pi_omega*marg/sum(pi_omega*marg)
        
        #constraint on weights
        if(constrained==1){
          ess <- c(a + b + n_1t, a + b + n_1t+ n_2t)
          totesss <- sum(wts*(ess-n_1t))
          
          if(k<looks){
            if(totesss > (n[1]*p[1])) {
              s <- (n[1]*p[1])/(totesss)
              wts <- c(wts[1]+(1-s)*(1-wts[1]), wts[2:length(wts)]*s)
            } }
        }
        
        
        #posterior: draw from individual posteriors according to weights
        samples <- cbind(samples, rmixture(wts=wts, betas=betas))
        
      }
      
      
      final <- samples[,2] - samples[,1]
      CI <- c(quantile(final, 0.025), quantile(final, 0.975))
      pointest <- mean(final)
      stopped <- quantile(final, bdryMEM) > 0
    }
    else if(method == "OBF") {
      
      bdry <- gsDesign(k=looks, test.type=1, alpha=0.025, beta=0.1, delta=0.5, timing=seq(0,1,1/looks)[-1], sfu="Pocock")
      bdry <- bdry$upper$bound
      fit <- prop.test(x = c(sum(combined$Y[combined$df==1 & combined$trt==1]), sum(combined$Y[combined$df==1 & combined$trt==0])), n = c(sum(combined$df==1 & combined$trt==1), sum(combined$df==1 & combined$trt==0)), alternative = "two.sided", correct = F)
      CI <- fit$conf.int
      pointest <- mean(combined$Y[combined$df==1 & combined$trt==1]) - mean(combined$Y[combined$df==1 & combined$trt==0])
      wts=rep(0, 2^(n_sources-1))
      fit <- prop.test(x = c(sum(combined$Y[combined$df==1 & combined$trt==1]), sum(combined$Y[combined$df==1 & combined$trt==0])), n = c(sum(combined$df==1 & combined$trt==1), sum(combined$df==1 & combined$trt==0)), alternative = "greater", correct = F)
      stopped <- fit$p.value < (1-pnorm(bdry[k]))
      
    }
    
    if (stopped==TRUE) break
  }
  
  
  return(list(CI=CI, mean= pointest, weights= wts, stoplook=k, sig=stopped))
}



nsim=10000
RNGkind("L'Ecuyer-CMRG")
pars=list(
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0.2, 0.2), a=1, b=1, pi_omega=c(0.9, 0.1), constrained=0, bdryMEM=0.0092, method="OBF"),
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0.2, 0.2), a=1, b=1, pi_omega=c(0.95, 0.05), constrained=0, bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0.2, 0.2), a=1, b=1, pi_omega=c(0.9, 0.1), constrained=0, bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0.2, 0.2), a=1, b=1, pi_omega=c(0.8, 0.2), constrained=1, bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0.2, 0.2), a=1, b=1, pi_omega=c(0.5, 0.5), constrained=1, bdryMEM=0.0092, method="MEM"),
  
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.2), a=1, b=1, pi_omega=c(0.9, 0.1), constrained=0, bdryMEM=0.0092, method="OBF"),
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.2), a=1, b=1, pi_omega=c(0.95, 0.05), constrained=0, bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.2), a=1, b=1, pi_omega=c(0.9, 0.1), constrained=0, bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.2), a=1, b=1, pi_omega=c(0.8, 0.2), constrained=1, bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.2), a=1, b=1, pi_omega=c(0.5, 0.5), constrained=1, bdryMEM=0.0092, method="MEM"),
  
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.1), a=1, b=1, pi_omega=c(0.9, 0.1), constrained=0, bdryMEM=0.0092, method="OBF"),
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.1), a=1, b=1, pi_omega=c(0.95, 0.05), constrained=0, bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.1), a=1, b=1, pi_omega=c(0.9, 0.1), constrained=0, bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.1), a=1, b=1, pi_omega=c(0.8, 0.2), constrained=1, bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), p0=c(0.4, 0.4), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.1), a=1, b=1, pi_omega=c(0.5, 0.5), constrained=1, bdryMEM=0.0092, method="MEM")
 )


#summary matrix
summ <- matrix(nrow=length(pars),ncol=9)

is.between <- function(x, a, b) {
  (x - a)  *  (b - x) >= 0
}


wt <- NULL
res <- list()
pb= txtProgressBar(min = 0, max = length(pars), style = 3, char=":)")
system.time(for(i in 1:length(pars)) {
  ## Define the parameters
  par <- pars[[i]] ####
  
  #### Do the simulation with these settings
  set.seed(3)
  sim.results <- t(mclapply((1:nsim),doSim,max=par$max, p0=par$p0, looks=par$looks, p=par$p, beta_t=par$beta_t, a=par$a, b=par$b, method = par$method, bdryMEM=par$bdryMEM, constrained=par$constrained, pi_omega=par$pi_omega, mc.cores=cores, mc.silent=T))
  result <- unlist(sim.results)
  result <- matrix(result, ncol=5+2^(length(par$max)-1), byrow=T)
  res[[i]] <- result
  #coverage: target of estimation is beta_t[1] in primary study
  cvrg <- is.between(as.numeric(par$beta_t[1]), result[,1], result[,2])
  
  ## Summarize the simulation results
  summ[i,1] <- sum(cvrg==T)/nsim
  summ[i,2] <- mean(result[,3] - as.numeric(par$beta_t[1]))
  summ[i,3] <- sum((result[,3] - as.numeric(par$beta_t[1]))**2)/nsim
  summ[i,4] <- mean(result[,ncol(result)-1])
  summ[i,5] <- sum(result[,ncol(result)-1]==1)/nsim
  summ[i,6] <- sum(result[,ncol(result)-1]==2)/nsim
  summ[i,7] <- sum(result[,ncol(result)-1]==3)/nsim
  summ[i,8] <- sum(result[,ncol(result)-1]==4)/nsim
  summ[i,9] <- mean(result[,ncol(result)])
  
  wt <- cbind(wt, result[,4])
  #sink('output2not.txt')
  print(summ[1:i,])
  setTxtProgressBar(pb, i)
  #sink()
})
close(pb)

summ

###########################################################################################
################################### Normal Outcome ##########################################
###########################################################################################


library(parallel)
library(gsDesign)

cores=8

#function that draws from the mixture for beta
rmixture = function(wts, betas){
  ##random generation of the indices
  id = sample(1:length(wts),prob=wts,size=nrow(betas),replace=TRUE)  
  id = cbind(1:nrow(betas),id)
  betas[id]
}


doSim <- function(seed, max, mu0, looks, p, beta_t, sd, method, bdryMEM=NULL, constrained=NULL, pi_omega=NULL) {
  ls = list()  
  for(k in 1:looks) {
    #gather 1/looks fraction of data, run the method, either stop early or loop again
    n=max*(1/looks)
    n_sources = length(n)
    
    #simulate data:
    #primary
    trt<- rep(0, n[1])
    trt[sample(1:n[1],n[1]*p[1])] <- 1
    Y <- rnorm(n[1], mean= mu0[1] + beta_t[1]*trt, sd=sd[1])
    primary <- data.frame(Y,trt)
    primary$df <- 1
    primary$secondary <- 0
    
    #secondary
    secondary <- NULL
    for(i in 2:n_sources) {
      trt<- rep(0, n[i])
      trt[sample(1:n[i],n[i]*p[i])] <- 1
      Y <- rnorm(n[i], mean= mu0[i] + beta_t[i]*trt, sd=sd[i])
      df <- i
      secondary <- rbind(secondary, data.frame(Y,trt,df))
    }
    secondary$secondary <- 1
    
    #data frame together
    dt <- data.frame(rbind(primary, secondary))
    
    names(dt) <- c("Y", "trt", "df", "secondary")
    #combined <- rbind(combined, dt)
    ls[[k]] <- dt
    combined <- data.frame(do.call(rbind, ls))
    
    if(method == "MEM") {
      samples <- NULL
      for(i in 0:1) {
        ybar_pt <- mean(combined$Y[combined$df==1 & combined$trt==i])
        ybar_1t <- mean(combined$Y[combined$df==2 & combined$trt==i])
        n_pt <- sum(combined$df==1 & combined$trt==i)
        n_1t <- sum(combined$df==2 & combined$trt==i)
        var_pt <- var(combined$Y[combined$df==1 & combined$trt==i])
        var_1t <- var(combined$Y[combined$df==2 & combined$trt==i])
        
        #model 1: no borrowing
        betas <- NULL
        marg <- NULL
        betas <- cbind(betas, rnorm(10000, ybar_pt, sqrt(var_pt/n_pt)))
        #marg <- c(marg, log(sqrt(2*pi)^2/sqrt((n_pt * n_1t)/(var_pt*var_1t))) +0.5*(n_pt*ybar_pt^2/var_pt + n_1t*ybar_1t^2/var_1t)) 
        marg <- log(2*pi/sqrt(1/((var_pt/n_pt)*(var_1t/n_1t))))
        
        #model 2: borrowing
        betas <- cbind(betas, rnorm(10000, (ybar_pt* n_pt/var_pt)/(n_pt/var_pt + n_1t/var_1t) + (ybar_1t* n_1t/var_1t)/(n_pt/var_pt + n_1t/var_1t), sqrt(1/(n_pt/var_pt + n_1t/var_1t))))
        #marg <- c(marg, log(sqrt(2*pi)/sqrt(n_pt/var_pt + n_1t/var_1t)) +0.5*(n_pt*ybar_pt/var_pt + n_1t*ybar_1t/var_1t)^2/(n_pt/var_pt + n_1t/var_1t)) 
        marg <- c(marg, log(sqrt(2*pi)/sqrt(1/((n_pt/var_pt) + (n_1t/var_1t)))) - 0.5*( (ybar_pt-ybar_1t)^2/((var_pt/n_pt) + (var_1t/n_1t)) ))
        
        #weights
        wts <- pi_omega*exp(marg-min(marg))/sum(pi_omega*exp(marg-min(marg)))
        #wts <- pi_omega*exp(marg)/sum(pi_omega*exp(marg))
        
        # #constraint on weights 
        if(constrained==1){
          rats <- c((n_pt/var_pt)/(n_pt/var_pt) -1, (n_pt/var_pt + n_1t/var_1t)/(n_pt/var_pt) -1)
          totesss <- n_pt *(sum(wts*rats))
          if(k<looks){
            if(totesss > (n[1]*p[1])) {
              s <- (n[1]*p[1])/(totesss)
              wts <- c(wts[1]+(1-s)*(1-wts[1]), wts[2:length(wts)]*s)
            } }
        }

        #posterior: draw from individual posteriors according to weights
        samples <- cbind(samples, rmixture(wts=wts, betas=betas))
        
      }
      
      
      final <- samples[,2] - samples[,1]
      CI <- c(quantile(final, 0.025), quantile(final, 0.975))
      pointest <- mean(final)
      stopped <- quantile(final, bdryMEM) > 0
    }
    else if(method == "OBF") {
      
      bdry <- gsDesign(k=looks, test.type=1, alpha=0.025, beta=0.1, delta=0.5, timing=seq(0,1,1/looks)[-1], sfu="Pocock")
      bdry <- bdry$upper$bound
      
      #2 sample t test
      fit <- t.test(x = combined$Y[combined$df==1 & combined$trt==1], y= combined$Y[combined$df==1 & combined$trt==0], alternative = "greater")
      CI <- fit$conf.int
      pointest <- mean(combined$Y[combined$df==1 & combined$trt==1]) - mean(combined$Y[combined$df==1 & combined$trt==0])
      wts=rep(0, 2^(n_sources-1))
      stopped <- fit$p.value < (1-pnorm(bdry[k]))
      
    }
    
    if (stopped==TRUE) break
  }
  
  
  return(list(CI=CI, mean= pointest, weights= wts, stoplook=k, sig=stopped))
}


nsim=10000
RNGkind("L'Ecuyer-CMRG")
pars=list(
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(1, 1), sd=c(3,4), bdryMEM=0.0092, constrained=0, pi_omega=c(0.95, 0.05), method="OBF"),
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(1, 1), sd=c(3,4), bdryMEM=0.0092, constrained=0, pi_omega=c(0.95, 0.05), method="MEM"),
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(1, 1), sd=c(3,4), bdryMEM=0.0092, constrained=0, pi_omega=c(0.9, 0.1), method="MEM"),
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(1, 1), sd=c(3,4), bdryMEM=0.0092, constrained=1, pi_omega=c(0.8, 0.2), method="MEM"),
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(1, 1), sd=c(3,4), bdryMEM=0.0092, constrained=1, pi_omega=c(0.5, 0.5), method="MEM"),
  
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(0, 1), sd=c(3,4), bdryMEM=0.0092, constrained=0, pi_omega=c(0.95, 0.05), method="OBF"),
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(0, 1), sd=c(3,4), bdryMEM=0.0092, constrained=0, pi_omega=c(0.95, 0.05), method="MEM"),
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(0, 1), sd=c(3,4), bdryMEM=0.0092, constrained=0, pi_omega=c(0.9, 0.1), method="MEM"),
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(0, 1), sd=c(3,4), bdryMEM=0.0092, constrained=1, pi_omega=c(0.8, 0.2), method="MEM"),
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(0, 1), sd=c(3,4), bdryMEM=0.0092, constrained=1, pi_omega=c(0.5, 0.5), method="MEM"),
  
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.5), sd=c(3,4), bdryMEM=0.0092, constrained=0, pi_omega=c(0.95, 0.05), method="OBF"),
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.5), sd=c(3,4), bdryMEM=0.0092, constrained=0, pi_omega=c(0.95, 0.05), method="MEM"),
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.5), sd=c(3,4), bdryMEM=0.0092, constrained=0, pi_omega=c(0.9, 0.1), method="MEM"),
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.5), sd=c(3,4), bdryMEM=0.0092, constrained=1, pi_omega=c(0.8, 0.2), method="MEM"),
  list(max=c(200, 400), mu0=c(5, 5), looks=4, p=c(0.5, 0.5), beta_t=c(0, 0.5), sd=c(3,4), bdryMEM=0.0092, constrained=1, pi_omega=c(0.5, 0.5), method="MEM")
  )

summ <- matrix(nrow=length(pars),ncol=9)

is.between <- function(x, a, b) {
  (x - a)  *  (b - x) >= 0
}

wt <- NULL
pb= txtProgressBar(min = 0, max = length(pars), style = 3, char=":)")
system.time(for(i in 1:length(pars)) {
  ## Define the parameters
  par <- pars[[i]] ####
  
  #### Do the simulation with these settings
  set.seed(1)
  sim.results <- t(mclapply(1:nsim,doSim,max=par$max, mu0=par$mu0, looks=par$looks, p=par$p, beta_t=par$beta_t, sd=par$sd, method = par$method, bdryMEM=par$bdryMEM, constrained=par$constrained, pi_omega=par$pi_omega, mc.cores=cores, mc.silent=T))
  result <- unlist(sim.results)
  result <- matrix(result, ncol=5+2^(length(par$max)-1), byrow=T)
  #coverage: target of estimation is beta_t[1] in primary study
  cvrg <- is.between(as.numeric(par$beta_t[1]), result[,1], result[,2])
  
  ## Summarize the simulation results
  summ[i,1] <- sum(cvrg==T)/nsim
  summ[i,2] <- mean(result[,3] - as.numeric(par$beta_t[1]))
  summ[i,3] <- sum((result[,3] - as.numeric(par$beta_t[1]))**2)/nsim
  summ[i,4] <- mean(result[,ncol(result)-1])
  summ[i,5] <- sum(result[,ncol(result)-1]==1)/nsim
  summ[i,6] <- sum(result[,ncol(result)-1]==2)/nsim
  summ[i,7] <- sum(result[,ncol(result)-1]==3)/nsim
  summ[i,8] <- sum(result[,ncol(result)-1]==4)/nsim
  summ[i,9] <- mean(result[,ncol(result)])
  
  wt <- cbind(wt, result[,4])
  #sink('output2not.txt')
  print(summ[1:i,])
  setTxtProgressBar(pb, i)
  #sink()
})
close(pb)


summ


###########################################################################################
################################### Linear Regression ##########################################
###########################################################################################



library(parallel)
library(R2jags)
library(mvtnorm)
library(gsDesign)

cores=10

#function that draws from the mixture for beta
rmixture = function(wts, betas){
  ##random generation of the indices
  id = sample(1:length(wts),prob=wts,size=nrow(betas),replace=TRUE)  
  id = cbind(1:nrow(betas),id)
  betas[id]
}

doSim <- function(seed, max, mu0, looks, beta_z, p, beta_t, beta_m, muz, varz, sd, var_b=NULL, a=NULL, b=NULL, pi_omega=NULL, method, burnin=NULL, ndraw=NULL, BICiterations=NULL, nint=NULL, marginal=NULL, constrained, bdryMEM) {
  
  lst = list()  
  for(k in 1:looks) {
    
    n=max*(1/looks)
    n_cur <- max*(k/looks)
    n_sources = length(n)
    n_z = length(beta_z)
    n_tot=sum(n_cur)
    n_m = length(beta_m)
    #simulate data:
    #primary
    trt <- rbinom(n[1], 1, p[1])
    Z <- rmvnorm(n[1], mean = muz[[1]], sigma=as.matrix(varz[[1]]))
    Y <- rnorm(n[1], mean= mu0[1] +Z%*%as.matrix(beta_z) + beta_t[1]*trt + (Z[,1:n_m]%*%as.matrix(beta_m))*trt, sd=sd[1])
    primary <- data.frame(Y,Z,trt)
    primary$df <- 1
    primary$secondary <- 0
    
    #secondary
    secondary <- NULL
    for(i in 2:n_sources) {
      trt <- rbinom(n[i], 1, p[i])
      Z <- rmvnorm(n[i], mean = muz[[i]], sigma=as.matrix(varz[[i]]))
      Y <- rnorm(n[i], mean= mu0[i] +Z%*%as.matrix(beta_z) + beta_t[i]*trt + (Z[,1:n_m]%*%as.matrix(beta_m))*trt, sd=sd[i])
      df <- i
      secondary <- rbind(secondary, data.frame(Y,Z,trt,df))
    }
    secondary$secondary <- 1
    
    #data frame together
    data <- data.frame(rbind(primary, secondary))
    
    Z<- as.matrix(data[,2:(length(beta_z)+1)])
    for(i in 1:n_sources) {
      for(j in 1:length(beta_z)) {
        x <- 1*(data$df==i) * (Z[,j])
        data <- cbind(data, x)
      }
    }
    
    Z<- as.matrix(data[,2:(length(beta_z)+1)])
    for(i in 1:n_sources) {
      for(j in 1:length(beta_z)) {
        x <- 1*(data$df==i) * (Z[,j] - mean(Z[data$secondary==0,j]))
        data <- cbind(data, x)
      }
    }
    
    #dummies for source
    for(i in 2:n_sources) {
      x <- 1*(data$df==i)
      data <- cbind(data, x)
    }
    
    #source*trt*z interactions
    for(i in 1:n_sources) {
      for(j in 1:n_m) {
        x <- 1*(data$df==i) * (Z[,j])*data$trt
        data <- cbind(data, x)
      }
      
    }
    
    for(i in 1:n_sources) {
      for(j in 1:n_m) {
        x <- 1*(data$df==i) * (Z[,j] - mean(Z[data$secondary==0,j]))*data$trt
        data <- cbind(data, x)
      }
      
    }
    
    
    #create interactions for source*trt which will decide whether I borrow or not on the trt effect
    for(i in 2:n_sources) {
      x <- 1*(data$df==i) * data$trt
      data <- cbind(data, x)
    }
    names(data) <- c("Y", paste("Z", 1:length(beta_z), sep=""), "trt", "df", "secondary", paste("Z", rep(1:n_sources, each=length(beta_z)), rep(1:length(beta_z), n_sources), sep=""), paste("Z_centered", rep(1:n_sources, each=length(beta_z)), rep(1:length(beta_z), n_sources), sep="") , paste("S", 2:n_sources, sep=""),paste("Z_trt", rep(1:n_sources, each=n_m), rep(1:n_m, n_sources), sep=""), paste("Z_centered_trt", rep(1:n_sources, each=n_m), rep(1:n_m, n_sources), sep=""),paste("trt_S", 2:n_sources, sep=""))
    
    data$type <- 1*(data$df==1 & data$trt==0)+ 2*(data$df==1 & data$trt==1) + 3*(data$df==2 & data$trt==0)+ 4*(data$df==2 & data$trt==1)+ 5*(data$df==3 & data$trt==0)+ 6*(data$df==3 & data$trt==1)
    
    
    lst[[k]] <- data
    combined <- data.frame(do.call(rbind, lst))
    combined <- combined[order(combined$df),]
    
    if(method == "MEM") {
      
      ls <- vector("list", n_sources-1)
      for(i in 1:(n_sources-1)) {
        x<-c(TRUE, FALSE)
        ls[[i]] <- x
      }
      
      regMat <- expand.grid(ls)
      regressors <- paste("trt_S", 2:n_sources, sep="")
      
      allModelsList <- apply(regMat, 1, function(x) as.formula(
        paste(c("Y ~ 1 + trt",paste("Z", rep(1:n_sources, each=length(beta_z)), rep(1:length(beta_z), n_sources), sep=""),paste("S", 2:n_sources, sep=""), regressors[x]),
              collapse=" + ")) )
      
      allModelsResults <- lapply(allModelsList,
                                 function(x) lm(x, data=combined))
      
      bic <- NULL
      aic <- NULL
      betas <- NULL
      for(i in 1:length(allModelsResults)) {
        fit <- allModelsResults[[i]]
        X <- as.matrix(model.matrix(fit))
        mub <- rep(0, ncol(X))
        VB <- var_b*diag(ncol(X))
        jags.inits <- function(){
          list("inv.var"=rep(1, n_sources))
        }
        
        #JAGS model:
        M <- function() {
          # Likelihood
          for(i in 1:n_tot){
            Y[i]   ~ dnorm(mu[i],prec[i])
            mu[i] <- X[i, ] %*% beta
          }
          
          # Prior for beta
          for(j in 1:J){
            beta[j] ~ dnorm(mub[j],1/(var_b))
          }
          
          # Prior for the inverse variance
          for(j in 1:nsource){
            inv.var[j]   ~ dgamma(a[j], b[j])
            sigma[j] <- 1/inv.var[j]
          }
          for(i in 1:n_tot){
            prec[i] <- inv.var[df[i]]
          }
        }
        
        
        out <- jags(list(Y=combined$Y, X=X, n_tot=nrow(combined), df=combined$df, nsource=n_sources, mub = mub, var_b= var_b, a = a, b =b,
                         J=ncol(X)),  inits = jags.inits, n.chains=1, n.thin=2, n.iter=ndraw, n.burnin = burnin, DIC=F, parameters.to.save =  c('beta'), M)
        
        out_mcmc <- as.mcmc(out)
        temp <- data.frame(unlist(out_mcmc[[1]]))
        names(temp) <- gsub("\\.", "", names(temp))
        names(temp) <- gsub("beta", "", names(temp))
        temp <- temp[ , order(as.numeric(names(temp)))]
        betas <- cbind(betas, temp[,2])
        
        if(marginal == "BIC" | marginal == "AIC") {
          #BIC weights:
          #iterated weighted least squares
          cf <- coef(fit)
          s <- rep(NA, n_sources)
          for(e in 1:BICiterations) {
            old <- cf
            for(j in 1:n_sources){
              s[j] <- sum((combined$Y[combined$df==j] - X[combined$df==j,] %*% cf)^2)/(length(combined$Y[combined$df==j]))
            }
            #W <- diag(rep(s, times=n))
            #cf <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% combined$Y
            
            cf <- coef(lm(combined$Y ~ X + 0, weights = sqrt(rep(s, times=n_cur))))
            if (sum(abs(old - cf)) < 0.0000000000001)
            {
              break;
            }
          }
          
          #BIC:
          bic <- c(bic, - log(nrow(combined))*(ncol(X)+n_sources) + 2*dmvnorm(combined$Y, mean= X %*% cf, sigma = diag(rep(s, times=n_cur)), log=TRUE))
          aic <- c(aic, - 2*(ncol(X)+n_sources) + 2*dmvnorm(combined$Y, mean= X %*% cf, sigma = diag(rep(s, times=n_cur)), log=TRUE))
          
        }
        
      }
      
      if(marginal == "BIC") {
        wts <- pi_omega*exp(0.5*(bic-max(bic)))/sum(pi_omega*exp(0.5*(bic-max(bic))))
      }else if(marginal == "AIC") {
        wts <- pi_omega*exp(0.5*(aic-max(aic)))/sum(pi_omega*exp(0.5*(aic-max(aic))))
      }
      
      if(constrained==1){
        #ESS calculation: use empirical estimate of posterior precision
        rats <- (1/diag(var(betas)))/(1/var(betas[,1])) -1
        totesss <- n[1] *(sum(wts*rats))
        if(k<looks){
          if(totesss > n[1]) {
            s <- n[1]/(totesss)
            wts <- c(wts[1]+(1-s)*(1-wts[1]), wts[2:length(wts)]*s)
          } }
      }
      
      #posterior: draw from individual posteriors according to weights
      samples <- rmixture(wts=wts, betas=betas)
      CI <- c(quantile(samples, 0.025), quantile(samples, 0.975))
      pointest <- mean(samples)
      stopped <- quantile(samples, bdryMEM) > 0
    } 
    else if(method == "OBF") {
      bdry <- gsDesign(k=looks, test.type=1, alpha=0.025, beta=0.1, delta=0.5, timing=seq(0,1,1/looks)[-1], sfu="Pocock")
      bdry <- bdry$upper$bound
      
      fmla <- as.formula(
        paste(c("Y ~ 1 + trt",paste("Z", 1, 1:n_z, sep="")),
              collapse=" + "))
      fit <- lm(fmla, data=combined[combined$df==1,])
      CI <- confint(fit)[2,]
      pointest <- fit$coef[2]
      wts=c(0,0)
      stopped <- pt(summary(fit)$coefficients[2,3], df=fit$df.residual, lower.tail = F) < (1-pnorm(bdry[1]))
      
    }
    if (stopped==TRUE) break
  }
  
  
  return(list(CI=CI, mean= pointest, weights= wts, stoplook=k, sig=stopped))
  
}


nsim=10000
RNGkind("L'Ecuyer-CMRG")
pars=list(
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(0,2), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.95, 0.05), constrained=0,marginal="BIC", bdryMEM=0.0092, method="OBF"),
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(0,1), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.95, 0.05), constrained=0,marginal="BIC", bdryMEM=0.0092, method="OBF"),
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(2,2), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.95, 0.05), constrained=0,marginal="BIC", bdryMEM=0.0092, method="OBF"),
  
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(0,2), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.95, 0.05), constrained=0,marginal="BIC", bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(0,1), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.95, 0.05), constrained=0,marginal="BIC", bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(2,2), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.95, 0.05), constrained=0,marginal="BIC", bdryMEM=0.0092, method="MEM"),
  
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(0,2), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.9, 0.1), constrained=0,marginal="BIC", bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(0,1), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.9, 0.1), constrained=0,marginal="BIC", bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(2,2), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.9, 0.1), constrained=0,marginal="BIC", bdryMEM=0.0092, method="MEM"),
  
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(0,2), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.5, 0.5), constrained=1, marginal="BIC", bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(0,1), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.5, 0.5), constrained=1,marginal="BIC", bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(2,2), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.5, 0.5), constrained=1,marginal="BIC", bdryMEM=0.0092, method="MEM"),
  
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(0,2), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.8, 0.2), constrained=1, marginal="BIC", bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(0,1), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.8, 0.2), constrained=1,marginal="BIC", bdryMEM=0.0092, method="MEM"),
  list(max=c(200, 400), mu0=c(2, 2), looks=4, beta_t=c(2,2), beta_z=c(0.8), p=c(0.5, 0.5), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=8000, BICiterations=1000, pi_omega=c(0.8, 0.2), constrained=1,marginal="BIC", bdryMEM=0.0092, method="MEM")
  
)


#summary matrix:
summ <- matrix(nrow=length(pars),ncol=9)

is.between <- function(x, a, b) {
  (x - a)  *  (b - x) >= 0
}

wt <- NULL
res <- list()
pb= txtProgressBar(min = 0, max = length(pars), style = 3, char=":)")
system.time(for(i in 1:length(pars)) {
  ## Define the parameters
  par <- pars[[i]] ####
  
  #### Do the simulation with these settings
  set.seed(3)
  sim.results <- t(mclapply((1:nsim),doSim,max=par$max, mu0=par$mu0, looks=par$looks, beta_z=par$beta_z, p=par$p, beta_t=par$beta_t, beta_m=c(0), muz=par$muz, varz=par$varz, var_b=par$var_b, sd=par$sd, a=par$a, b=par$b, BICiterations=par$BICiterations, burnin=par$burnin, ndraw=par$ndraw, method = par$method, marginal=par$marginal, constrained=par$constrained, bdryMEM=par$bdryMEM, pi_omega=par$pi_omega, mc.cores=cores, mc.silent=T))
  result <- unlist(sim.results)
  result <- matrix(result, ncol=7, byrow=T)
  res[[i]] <- result
  
  #coverage: target of estimation is beta_t[1] in primary study
  cvrg <- is.between(as.numeric(par$beta_t[1]), result[,1], result[,2])
  
  ## Summarize the simulation results
  summ[i,1] <- sum(cvrg==T)/nsim
  summ[i,2] <- mean(result[,3] - as.numeric(par$beta_t[1]))
  summ[i,3] <- sum((result[,3] - as.numeric(par$beta_t[1]))**2)/nsim
  summ[i,4] <- mean(result[,ncol(result)-1])
  summ[i,5] <- sum(result[,ncol(result)-1]==1)/nsim
  summ[i,6] <- sum(result[,ncol(result)-1]==2)/nsim
  summ[i,7] <- sum(result[,ncol(result)-1]==3)/nsim
  summ[i,8] <- sum(result[,ncol(result)-1]==4)/nsim
  summ[i,9] <- mean(result[,ncol(result)])
  
  wt <- cbind(wt, result[,4])
  #sink('output2not.txt')
  print(summ[1:i,])
  setTxtProgressBar(pb, i)
  #sink()
})
close(pb)

summ
