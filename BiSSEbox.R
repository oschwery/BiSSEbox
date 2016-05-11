# load packages
library(ape)
library(geiger)
library(diversitree)
#Code to run BiSSEbox
BiSSEbox <- function(tree, data, Nsteps=5000) {  # later: samp.freq, taxofile
  # generate output objects
  traitMLfits <- list()
  MCMCs <- list()
  AICs <- list()
  logLiks <- list()
  coefs <- list()
  # vector of trait names
  traits <- colnames(data)

  for (i in 1:length(traits)) {
    ##chose desired setup (no missing sampling, richness  file, sampling frequency)
    print(paste("Current trait:", traits[i], sep=" "))
    current.trait <- data[, i]
    names(current.trait) <- row.names(data)
    print("Running ML fits.")
    lik <- make.bisse(tree, current.trait)
    #lik <- make.bisse(tree, SLA.v, unresolved=richnessSLA)
    #sampling.f_SLA<-c(0.051095,0.046405)
    #lik <- make.bisse(tree, SLA.v, unresolved=sampling.f_SLA)
    #sampling.f_SLAdistr<-c(0.047726,0.044845)
    #lik <- make.bisse(tree, SLA.v, unresolved=sampling.f_SLAdistr)

    ##unconstrained run
    p <- starting.point.bisse(tree)
    fit <- find.mle(lik, p)
    #fit$logLik
    #round(coef(fit), 3)
    ##constrained run ALL MODELS
    #Constrain lambda, and/or mu and/or q
    lik.l <- constrain(lik, lambda0 ~ lambda1)
    lik.m<- constrain(lik, mu0 ~ mu1)
    lik.q <- constrain(lik, q01 ~ q10)
    lik.lm<- constrain(lik, mu0 ~ mu1, lambda0 ~ lambda1)
    lik.lq<- constrain(lik, q01 ~ q10, lambda0 ~ lambda1)
    lik.mq<- constrain(lik, q01 ~ q10,  mu0 ~ mu1)
    lik.lmq<- constrain(lik, mu0 ~ mu1, lambda0 ~ lambda1, q01 ~ q10)
    # run fits for liks
    fit.l <- find.mle(lik.l, p[argnames(lik.l)])
    fit.m<- find.mle(lik.m, p[argnames(lik.m)])
    fit.q <- find.mle(lik.q, p[argnames(lik.q)])
    fit.lm<- find.mle(lik.lm, p[argnames(lik.lm)])
    fit.lq <- find.mle(lik.lq, p[argnames(lik.lq)])
    fit.mq<- find.mle(lik.mq, p[argnames(lik.mq)])
    fit.lmq<- find.mle(lik.lmq, p[argnames(lik.lmq)])
    #save the fits
    traitMLfits[[i]] <- list(unconstrained=fit, equal.lambda=fit.l, equal.mu=fit.m, equal.q=fit.q, equal.lambda.mu=fit.lm, equal.lambda.q=fit.lq, equal.mu.q=fit.mq, all.equal=fit.lmq)
    names(traitMLfits)[[i]] <- paste(traits[i])
    #round(rbind(full=coef(fit), equal.l=coef(fit.l, TRUE)), 3)
    #round(rbind(full=coef(fit), equal.m=coef(fit.m, TRUE)), 3)
    #round(rbind(full=coef(fit), equal.q=coef(fit.q, TRUE)), 3)
    #round(rbind(full=coef(fit), equal.lm=coef(fit.lm, TRUE)), 3)
    #round(rbind(full=coef(fit), equal.lq=coef(fit.lq, TRUE)), 3)
    #round(rbind(full=coef(fit), equal.mq=coef(fit.mq, TRUE)), 3)
    #round(rbind(full=coef(fit), equal.lmq=coef(fit.lmq, TRUE)), 3)

    ##ANOVA to test if constrained and unconstrained model perform different
    #anova(fit, equal.l=fit.l, equal.m=fit.m, equal.q=fit.q, equal.lm=fit.lm, equal.lq=fit.lq, equal.mq=fit.mq, equal.lmq=fit.lmq )
    AllLiks <- c(lik, lik.l, lik.m, lik.q, lik.lm, lik.lq, lik.mq, lik.lmq)
    AllFits <- list(fit, fit.l, fit.m, fit.q, fit.lm, fit.lq, fit.mq, fit.lmq)
    Modelnames <- c("unconstrained", "equal.lambda", "equal.mu", "equal.q", "equal.lambda.mu", "equal.lambda.q", "equal.mu.q", "all.equal")
    # add AIC's hethere...
    currentAICs <- c()
    for (j in 1:length(Modelnames)) {
      currentAICs <- c(currentAICs, AIC(AllFits[[j]]))
    }
    names(currentAICs) <- Modelnames
    # Likelihoods
    currentlogLiks <- c()
    for (j in 1:length(Modelnames)) {
      currentlogLiks <- c(currentlogLiks, logLik(AllFits[[j]]))
    }
    names(currentlogLiks) <- Modelnames
    # estimated rates
    currentcoefs <- list()
    for (j in 1:length(Modelnames)) {
      currentcoefs[[j]] <- coef(AllFits[[j]])
    }
    names(currentcoefs) <- Modelnames
    # add all of these to list and name by current trait
    AICs[[i]] <- currentAICs
    logLiks[[i]] <- currentlogLiks
    coefs[[i]] <- currentcoefs
    names(AICs)[[i]] <- paste(traits[i])
    names(logLiks)[[i]] <- paste(traits[i])
    names(coefs)[[i]] <- paste(traits[i])
#    suffix <- c("l", "m", "ml")
#    for (i in sequence(lenght))

    #selected.model.lik <- AllLiks[max(which(currentAICs == min(currentAICs)))]
    #selected.model.fit <- AllFits[[max(which(currentAICs == min(currentAICs)))]]
    #selected.model.name <- Modelnames[max(which(currentAICs == min(currentAICs)))]
    # find model with smalles AIC and get corresponding likelihood function, fit and name
    selected.model.lik <- AllLiks[[which(currentAICs == min(currentAICs))]]
    selected.model.fit <- AllFits[[which(currentAICs == min(currentAICs))]]
    selected.model.name <- Modelnames[which(currentAICs == min(currentAICs))]
    # display name of chosen model
    print(selected.model.name)
    # run MCMC
    # make prior
    prior <- make.prior.exponential(2*(log(length(tree$tip.label))/(max(branching.times(tree)))))
    # run MCMC short to determine w
    print("Start MCMC calibration")
    mcmc.bisse<-mcmc(selected.model.lik,selected.model.fit$par,nsteps=100,prior=prior,w=0.1)
    mcmc.bisse
    w=diff(sapply(mcmc.bisse[2:(ncol(mcmc.bisse)-1)],quantile,c(0.05,0.95)))
    # run actual MCMC
    print("Start MCMC")
    #For real...
    mcmc.bisse2<-mcmc(selected.model.lik,selected.model.fit$par,nsteps=Nsteps,w=w,prior=prior)
#    traitMLfits <- c(traitMLfits, MLfits)
    # save MCMC result and name by trait and model
    MCMCs[[i]] <- mcmc.bisse2
    names(MCMCs)[[i]] <- paste(traits[i], selected.model.name, sep="_")
  }
  return(list(fits=traitMLfits, AICs=AICs, MCMCs=MCMCs))
}
