### BiSSE ERICACEAE

### variables to replace here (apart from paths, but change those AFTER):
### trait file, richness file, trait vector...

#######################################################################
###   preparation                      ################################
#######################################################################


##load packages
library(ape)
library(geiger)
library(diversitree)

##read tree
tree<-read.nexus("C:\\Users\\Orlando\\Documents\\UZH\\Master\\Publication\\Data\\Phylogeny\\BEAST\\Final_files_superlast\\Ericac_rbcL_matK_Ranunc1_R1_f3_superlast_mean_noOG.tre")

##read trait states (0,1,NA)
SLA<-read.table("C:\\Users\\Orlando\\Documents\\UZH\\Master\\Publication\\Analyses\\BiSSE\\BiSSE_Paper\\BiSSE_SLA_fullR1.txt", row.names=1)

##read richness file for unresolved sampling
#richnessSLA<- read.csv("C:\\Users\\Orlando\\Documents\\UZH\\Master\\Thesis\\Analyses\\BiSSE\\Rich_eric_SLA.csv", as.is=TRUE)

##this is just to check
##head(richnessSLA)
head(SLA)

##make vector out of trait and assign names to it
SLA.v<-SLA[,1]
names(SLA.v)<-row.names(SLA)

#######################################################################
###   Model Choice                     ################################
#######################################################################


##chose desired setup (no missing sampling, richness file, sampling frequency)
lik <- make.bisse(tree, SLA.v)
#lik <- make.bisse(tree, SLA.v, unresolved=richnessSLA)
#sampling.f_SLA<-c(0.051095,0.046405)
#lik <- make.bisse(tree, SLA.v, unresolved=sampling.f_SLA)
#sampling.f_SLAdistr<-c(0.047726,0.044845)
#lik <- make.bisse(tree, SLA.v, unresolved=sampling.f_SLAdistr)

##unconstrained run
p <- starting.point.bisse(tree)
fit <- find.mle(lik, p)
fit$lnLik
round(coef(fit), 3)

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
fit.m<- find.mle (lik.m, p[argnames(lik.m)])
fit.q <- find.mle(lik.q, p[argnames(lik.q)])
fit.lm<- find.mle (lik.lm, p[argnames(lik.lm)])
fit.lq <- find.mle(lik.lq, p[argnames(lik.lq)])
fit.mq<- find.mle (lik.mq, p[argnames(lik.mq)])
fit.lmq<- find.mle (lik.lmq, p[argnames(lik.lmq)])

round(rbind(full=coef(fit), equal.l=coef(fit.l, TRUE)), 3)
round(rbind(full=coef(fit), equal.m=coef(fit.m, TRUE)), 3)
round(rbind(full=coef(fit), equal.q=coef(fit.q, TRUE)), 3)
round(rbind(full=coef(fit), equal.lm=coef(fit.lm, TRUE)), 3)
round(rbind(full=coef(fit), equal.lq=coef(fit.lq, TRUE)), 3)
round(rbind(full=coef(fit), equal.mq=coef(fit.mq, TRUE)), 3)
round(rbind(full=coef(fit), equal.lmq=coef(fit.lmq, TRUE)), 3)

##ANOVA to test if constrained and unconstrained model perform different
anova(fit, equal.l=fit.l, equal.m=fit.m, equal.q=fit.q, equal.lm=fit.lm, equal.lq=fit.lq, equal.mq=fit.mq, equal.lmq=fit.lmq )

### NOW DECIDE FOR SIMPLEST SIGNIFICANT MODEL ###



#######################################################################
###   MCMC  FULL                       ################################
#######################################################################
#
##MCMC syndrome FULL MODEL
#lik <- make.bisse(tree, MtAss.v)
#
#p <- starting.point.bisse(tree)
#prior <- make.prior.exponential(2*(log(450)/112.995))
#
#mcmc.bisse<-mcmc(lik,fit$par,nsteps=100,prior=prior,w=0.1)
#mcmc.bisse
#w=diff(sapply(mcmc.bisse[2:7],quantile,c(0.05,0.95)))
#
##For real...
#mcmc.bisse2<-mcmc(lik,fit$par,nsteps=2000,w=w,prior=prior)
#
##Now we can examine the 95% credible intervals of the posterior samples for each parameter.  If the intervals do not overlap, then we have posterior Bayesian support for a difference in rates.
#sapply(mcmc.bisse2[,2:7],quantile,c(0.025,0.975))
#
##To plot this
#col <- c("#004165", "#eaab00")
#windows()
#profiles.plot(mcmc.bisse2[c("lambda1", "lambda0")], col.line=col, las=1,xlab="Speciation rate", legend="topright")
#windows()
#profiles.plot(mcmc.bisse2[c("mu1", "mu0")], col.line=col, las=1,xlab="Extinction rate", legend="topright")
#windows()
#profiles.plot(mcmc.bisse2[c("q01", "q10")], col.line=col, las=1,xlab="Transition rate", legend="topright")
#
#mcmc.bisse2$r0=mcmc.bisse2$lambda0-mcmc.bisse2$mu0
#mcmc.bisse2$r1=mcmc.bisse2$lambda1-mcmc.bisse2$mu1
#
#windows()
#profiles.plot(mcmc.bisse2[c("r0", "r1")], col.line=col, las=1,xlab="Diversification rate", legend="topright")
#
#write.csv(mcmc.bisse2, file="syndrome_all_full_model.csv")
#


#######################################################################
###   MCMC  CONSTRAINED                       #########################
#######################################################################

#################MCMC syndrome CONSTRAINED MODEL (lambda, q)#####################################################

prior <- make.prior.exponential(2*(log(450)/112.995))

mcmc.bisse<-mcmc(lik.mq,fit.mq$par,nsteps=100,prior=prior,w=0.1)
mcmc.bisse
w=diff(sapply(mcmc.bisse[2:5],quantile,c(0.05,0.95)))

#For real...
mcmc.bisse2<-mcmc(lik.mq,fit.mq$par,nsteps=5000,w=w,prior=prior)

#Now we can examine the 95% credible intervals of the posterior samples for each parameter.  If the intervals do not overlap, then we have posterior Bayesian support for a difference in rates.
sapply(mcmc.bisse2[,2:6],quantile,c(0.025,0.975))

#To plot this
col <- c("#004165", "#eaab00")
windows()
profiles.plot(mcmc.bisse2[c("lambda1", "lambda0")], col.line=col, las=1,xlab="Speciation rate", legend="topright")
windows()
profiles.plot(mcmc.bisse2[c("mu1", "mu0")], col.line=col, las=1,xlab="Extinction rate", legend="topright")
windows()
profiles.plot(mcmc.bisse2[c("q01", "q10")], col.line=col, las=1,xlab="Transition rate", legend="topright")

mcmc.bisse2$r0=mcmc.bisse2$lambda0-mcmc.bisse2$mu1
mcmc.bisse2$r1=mcmc.bisse2$lambda1-mcmc.bisse2$mu1

windows()
profiles.plot(mcmc.bisse2[c("r1", "r0")], col.line=col, las=1,xlab="Diversification rate", legend="topright")

write.csv(mcmc.bisse2, file="SLA_full_constrained_model_superlast_poster.csv")

##############################################################################
