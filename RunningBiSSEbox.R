# Settings code for BiSSEbox
#office paths

source("/Users/orlandoschwery/Documents/UT/Courses&Workshops/Spring16/Phylometh/BiSSEbox/BiSSEbox.R")

tree <- read.tree("/Users/orlandoschwery/Documents/UT/Courses&Workshops/Spring16/Phylometh/BiSSEbox/BiSSEbox_testtree.tre")

data <- read.csv("/Users/orlandoschwery/Documents/UT/Courses&Workshops/Spring16/Phylometh/BiSSEbox/BiSSEbox_testdata.csv", row.names=1, header=TRUE)

#parameters <-

ALLtheBiSSE <- BiSSEbox(tree, data, 200)


get the summary
get the individual results
get the plots


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
