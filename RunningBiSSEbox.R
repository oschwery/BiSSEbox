## Settings code for BiSSEbox
#
source("/Users/orlandoschwery/Documents/UT/Courses&Workshops/Spring16/Phylometh/BiSSEbox/BiSSEbox.R")

tree <- read.tree("/Users/orlandoschwery/Documents/UT/Courses&Workshops/Spring16/Phylometh/BiSSEbox/BiSSEbox_testtree.tre")

data <- read.csv("/Users/orlandoschwery/Documents/UT/Courses&Workshops/Spring16/Phylometh/BiSSEbox/BiSSEbox_testdata.csv", row.names=1, header=TRUE)

ALLtheBiSSE <- BiSSEbox(tree, data, 200)

#Now we can examine the 95% credible intervals of the posterior samples for each parameter.  If the intervals do not overlap, then we have posterior Bayesian support for a difference in rates.
sapply(as.data.frame(ALLtheBiSSE$MCMCs[[1]])[,2:(ncol(as.data.frame(ALLtheBiSSE$MCMCs[1]))-1)],quantile,c(0.025,0.975))

#To plot this
col <- c("#004165", "#eaab00")
quartz()
profiles.plot(as.data.frame(ALLtheBiSSE$MCMCs[[1]])[c("lambda1", "lambda0")], col.line=col, las=1,xlab="Speciation rate", legend="topright")
