# load packages
library(ape)
library(geiger)
library(diversitree)
#Code to run BiSSEbox
BiSSEbox <- function(tree, data) {
  # generate output objects
  traitMLfits <- list()
  MCMCs <- list()

  for (i in length(traits) {
    ##chose desired setup (no missing sampling, richness  file, sampling frequency)
    current.trait <- traits[, i]
    names(current.trait) <- row.names(traits)
    lik <- make.bisse(tree, current.trait)
    #lik <- make.bisse(tree, SLA.v, unresolved=richnessSLA)
    #sampling.f_SLA<-c(0.051095,0.046405)
    #lik <- make.bisse(tree, SLA.v, unresolved=sampling.f_SLA)
    #sampling.f_SLAdistr<-c(0.047726,0.044845)
    #lik <- make.bisse(tree, SLA.v, unresolved=sampling.f_SLAdistr)

    ##unconstrained run
    p <- starting.point.bisse(tree)
    fit <- find.mle(lik, p)
    #fit$lnLik
    #round(coef(fit), 3)

    ##constrained run ALL MODELS
    #Contrain lambda, and/or mu and/or q

    lik.l <- constrain(lik, lambda0 ~ lambda1)
    lik.m<- constrain(lik, mu0 ~ mu1)
    lik.q <- constrain(lik, q01 ~ q10)
    lik.lm<- constrain(lik, mu0 ~ mu1, lambda0 ~ lambda1)
    lik.lq<- constrain(lik, q01 ~ q10, lambda0 ~ lambda1)
    lik.mq<- constrain(lik, q01 ~ q10,  mu0 ~ mu1)
    lik.lmq<- constrain(lik, mu0 ~ mu1, lambda0 ~ lambda1, q01 ~ q10)

    fit.l <- find.mle(lik.l, p[argnames(lik.l)])
    fit.m<- find.mle(lik.m, p[argnames(lik.m)])
    fit.q <- find.mle(lik.q, p[argnames(lik.q)])
    fit.lm<- find.mle(lik.lm, p[argnames(lik.lm)])
    fit.lq <- find.mle(lik.lq, p[argnames(lik.lq)])
    fit.mq<- find.mle(lik.mq, p[argnames(lik.mq)])
    fit.lmq<- find.mle(lik.lmq, p[argnames(lik.lmq)])
    #save the fits
    MLfits <- list(unconstrained=fit, equal.lambda=fit.l, equal.mu=fit.m, equal.q=fit.q, equal.lambda.mu=fit.lm, equal.lambda.q=fit.lq, equal.mu.q=fit.mq, all.equal=fit.lmq)
    #round(rbind(full=coef(fit), equal.l=coef(fit.l, TRUE)), 3)
    #round(rbind(full=coef(fit), equal.m=coef(fit.m, TRUE)), 3)
    #round(rbind(full=coef(fit), equal.q=coef(fit.q, TRUE)), 3)
    #round(rbind(full=coef(fit), equal.lm=coef(fit.lm, TRUE)), 3)
    #round(rbind(full=coef(fit), equal.lq=coef(fit.lq, TRUE)), 3)
    #round(rbind(full=coef(fit), equal.mq=coef(fit.mq, TRUE)), 3)
    #round(rbind(full=coef(fit), equal.lmq=coef(fit.lmq, TRUE)), 3)

    ##ANOVA to test if constrained and unconstrained model perform different
    anova(fit, equal.l=fit.l, equal.m=fit.m, equal.q=fit.q, equal.lm=fit.lm, equal.lq=fit.lq, equal.mq=fit.mq, equal.lmq=fit.lmq )

# add AIC's hethere...
AIC(fit)
logLik(fit)
coef(fit)


    AIC.uncon <- (-2*(ln(likelihood))+(2*6))
    AIC.l <-
    AIC.m <-
    AIC.q <-
    AIC.lm <-
    AIC.lq <-
    AIC.mq <-
    AIC.lmq <-

    suffix <- c("l", "m", "ml")
    for (i in sequence(lenght))

    selected.model.lik <- ...
    selected.model.fit <- ...
    # run MCMC of that one
    prior <- make.prior.exponential(2*(log(length(tree$tip.label))/(max(branching.times(tree)))))

    mcmc.bisse<-mcmc(selected.model.lik,selected.model.fit$par,nsteps=100,prior=prior,w=0.1)
    mcmc.bisse
    w=diff(sapply(mcmc.bisse[2:5],quantile,c(0.05,0.95)))

    #For real...
    mcmc.bisse2<-mcmc(selected.model.lik,selected.model.fit$par,nsteps=5000,w=w,prior=prior)

    #Now we can examine the 95% credible intervals of the posterior samples for each parameter.  If the intervals do not overlap, then we have posterior Bayesian support for a difference in rates.
    sapply(mcmc.bisse2[,2:6],quantile,c(0.025,0.975))


    save outputs (all the fit, round table of fits, anova, aic, chosen model, mcmc.bisse2)

    traitMLfits <- c(traitMLfits, MLfits)

    MCMCs <- c(MCMCs, mcmc.bisse2)


    clear all the objects?
    MLfits <- list()
  }

  return(list(fits=traitMLfits, tablesofrates=roundtableofallthefits, ANOVA=anova, AICs=AICs, ChosenModels=ChosenModels, MCMCs=MCMCs))

  Save output to allow seeing all and plotting all and having a summary of the results
}
