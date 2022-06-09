###
# Reversible Jump MCMC algorithm:
# applied to estimate the diversification rates over time.
#
# author: Sebastian Hoehna
#
#library(TESS) loads packages ape and coda and desolve as dependencies
library(TESS)

treefile <- "pruned"
analysis_name <- sprintf("CoMET_%s",treefile)
CDT <- "survival"
EXPECTED_NUM_EVENTS <- 2
MCMC_ITERATIONS <- 1000000000
ALLOW_MASS_EXTINCTION <- TRUE
if ( ALLOW_MASS_EXTINCTION == TRUE ) {
   analysis_name <- sprintf("CoMET_%s_ME",treefile)
} else {
   analysis_name <- sprintf("CoMET_%s",treefile)
}

tree <- read.tree(file=sprintf("data/%s.tre",treefile) ) 

#Assign sampling value at present

rho <- 1.0

#If there is incomplete taxon sampling, estimate as:

total <- 28
rho <- (tree$Nnode+1)/total


# Assign priorForms; in this case, a lognormal prior is used for speciation and extinction rates, defined by mean=0.0 and stdev=1.0. Number of rate shifts in extinction and speciation and number MEEs are modeled by three independent Compound Poisson Process with a distribution hyperprior with lambda=2 (0.5 probability given to 0 rate shifts)

#If you do not want to fix the prior distribution but let TESS estimate it with empiricalHyperPriors=TRUE, then leave open the brackets as in TESS manual.
#priorForms <- c()

priorForms <- c("lognormal","normal","gamma")
priorForms <- c("lognormal")

#################################################### RUN ANALYSIS ##############################################################
#To run an analysis with no MEEs, use the code [1] below. This will produce an output called "CoMET_pruned" containing the files with parameter values. It will also produce at the end a PDF figure named "CoMET_pruned.pdf" containing vignettes for 4 figures. This function will create a directory with the name of "analysis_name" given above in the current working directory getwd()

[1]

tess.analysis(tree=tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, numExpectedMassExtinctions=EXPECTED_NUM_EVENTS, initialSpeciationRate=2.0, initialExtinctionRate=1.0, empiricalHyperPriorInflation = 10.0, empiricalHyperPriorForm = priorForms, samplingProbability=rho, estimateMassExtinctionTimes = ALLOW_MASS_EXTINCTION, estimateNumberMassExtinctions = ALLOW_MASS_EXTINCTION, MAX_ITERATIONS = MCMC_ITERATIONS, THINNING = 10,  MAX_TIME = Inf, MIN_ESS = 1000, CONDITION=CDT, dir = analysis_name)  


#To run an analysis with MEEs, use the code [2] below. This will produce an output folder called "CoMET_ME_pruned" containing the files with parameter posterior MCMC values. It will also produce at the end a PDF figure named "CoMET_ME_pruned.pdf" containing vignettes for six figures.

[2]

ALLOW_MASS_EXTINCTION <- TRUE

tess.analysis(tree=tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, numExpectedMassExtinctions=EXPECTED_NUM_EVENTS, initialSpeciationRate=2.0, initialExtinctionRate=1.0, empiricalHyperPriorInflation = 10.0, empiricalHyperPriorForm = priorForms, samplingProbability=rho, estimateMassExtinctionTimes = ALLOW_MASS_EXTINCTION, estimateNumberMassExtinctions = ALLOW_MASS_EXTINCTION, MAX_ITERATIONS = MCMC_ITERATIONS, THINNING = 10,  MAX_TIME = Inf, MIN_ESS = 1000, CONDITION=CDT, dir = analysis_name)
          
######## CHANGE SAMPLING STRATEGY FROM UNIFORM TO DIVERSIFIED ##########################
# INCOMPLETE TAXON SAMPLING strategy. Default in tess.analysis function is "samplingStrategy=uniform" (species are sampled uniformly at random), which means that each species has the same probability (rho) of being sampled in the phylogeny. 
# To apply alternative: "samplingStrategy=diversified" (species are sampled to maximize the diversity sampled in the phylogeny: i.e., only the oldest 25% of divergence events are included in the reconstructed phylogeny with sampling probability rho, and all later divergence events are excluded), we can use the function "tess.analysis.diversified.R" created by modifying the argument "samplingStrategy=diversified" in the "tess.likelihood" function part of the code. First, source it:

source("tess.analysis.diversified.R")

tess.analysis.diversified(tree=tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, numExpectedMassExtinctions=EXPECTED_NUM_EVENTS, initialSpeciationRate=2.0, initialExtinctionRate=1.0, empiricalHyperPriorInflation = 10.0, empiricalHyperPriorForm = priorForms, samplingProbability=rho, estimateMassExtinctionTimes = ALLOW_MASS_EXTINCTION, estimateNumberMassExtinctions = ALLOW_MASS_EXTINCTION, MAX_ITERATIONS = MCMC_ITERATIONS, THINNING = 100,  MAX_TIME = Inf, MIN_ESS = 1000, CONDITION=CDT, dir = analysis_name)   

############################################# PROCESS ANALYSIS OUTPUT ########################################################
### NO MEEs
out <- tess.process.output(analysis_name, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)

#If there are two analyses run: run1 and run2, then we can use the following commands to construct a combined sample to plot
output1 <- tess.process.output(analysis_name1, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)
output2 <- tess.process.output(analysis_name2, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)

outputs <- list(output1, output2)
#And use this for any analysis
test.plot.multichain.diagnostics(outputs)


### WITH MEEs

out2 <- tess.process.output(analysis_name, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, numExpectedMassExtinctions=EXPECTED_NUM_EVENTS)

#Arguments: tree = NULL, numExpectedRateChanges = numExpectedMassExtinctions = EXPECTED_NUM_EVENTS = 2, burnin = 0.25, numIntervals = 100, criticalBayesFactors = c(2,6,10)

############################################ CHECK MCMC DIAGNOSTICS WITH CODA ################################################
#If only one chain and with EES diagnostic
#Remember to give only one parameter or divide into vignettes with "par" function. Otherwise, plots will be overlaid. Similarly, choose between "EES" and "geweke" to avoid overlaying of values. There are other arguments, which we do not consider here: (col=NULL,xaxt="n",yaxt="s")
pdf(file = "mass extinction times" )
tess.plot.singlechain.diagnostics(out2,parameters=c("mass extinction times"),diagnostics=c("ESS"),ess.crit=c(100,200),correction="bonferroni",xlab="million years ago",pch=19)
dev.off()
pdf(file = "speciation rates" )
tess.plot.singlechain.diagnostics(out2,parameters=c("speciation rates"),diagnostics=c("ESS"),ess.crit=c(100,200),correction="bonferroni",xlab="million years ago",pch=19)
dev.off()
pdf(file = "extinction rates" )
tess.plot.singlechain.diagnostics(out2,parameters=c("extinction rates"),diagnostics=c("ESS"),ess.crit=c(100,200),correction="bonferroni",xlab="million years ago",pch=19)
dev.off()
pdf(file = "speciation shift times")
tess.plot.singlechain.diagnostics(out2,parameters=c("speciation shift times"),diagnostics=c("ESS"),ess.crit=c(100,200),correction="bonferroni",xlab="million years ago",pch=19)
dev.off()
pdf(file = "extinction shift times")
tess.plot.singlechain.diagnostics(out2,parameters=c("extinction shift times"),diagnostics=c("ESS"),ess.crit=c(100,200),correction="bonferroni",xlab="million years ago",pch=19)
dev.off()
pdf(file = "net-diversification rates")
tess.plot.singlechain.diagnostics(out2,parameters=c("net-diversification rates"),diagnostics=c("ESS"),ess.crit=c(100,200),correction="bonferroni",xlab="million years ago",pch=19)
dev.off()
pdf(file = "relative-extinction rates")
tess.plot.singlechain.diagnostics(out2,parameters=c("relative-extinction rates"),diagnostics=c("ESS"),ess.crit=c(100,200),correction="bonferroni",xlab="million years ago",pch=19)
dev.off()
#tess.plot.singlechain.diagnostics(out,parameters=c("speciation rates", "extinction rates", "speciation shift times", "extinction shift times", "net-diversification rates", "relative-extinction rates", "mass extinction times"),diagnostics=c("ESS"),ess.crit=c(100,200),correction="bonferroni",xlab="million years ago",pch=19)

#If multiple chains and with geweke diagnostics

tess.plot.multiplechain.diagnostics(outputs,parameters=c("speciation rates", "extinction rates", "speciation shift times", "extinction shift times", "net-diversification rates", "relative-extinction rates", "mass extinction times"),diagnostics=c("geweke"),geweke.crit=0.05,correction="bonferroni",xlab="million years ago",pch=19)

############################################ PLOT RESULTS IN VIGNETTES #######################################################
#Now plot the results. Results will be plotted as vignettes within a PDF: nrows = 2, ncolumns= FIG_TYPES/2
#If no MEEs are allowed as in the code [1] above, only 4 figures will be plotted; if ALLOW_MASS_EXTINCTION = TRUE as in the code [2] above, 6 figures will be plotted as vignettes. Alternatively, 8 vignettes including everything can be plotted. Rate shift times and mass extinction times are plotted together with Bayes Factors (2lnBF) significance levels at 2,6,10 on the right y axis (this can be assigned a name with argument "yaxt")

NUM_FIGS <- 6
FIG_TYPES <- c("speciation rates","speciation shift times","extinction rates","extinction shift times", "net-diversification rates", "relative-extinction rates", "mass extinction times")
FIG_TYPES <- c("speciation rates","speciation shift times","extinction rates","extinction shift times","mass extinction Bayes factors","mass extinction times")
#if ( ALLOW_MASS_EXTINCTION == FALSE ) {
#   NUM_FIGS <- 4
#   FIG_TYPES <- c("speciation rates","extinction rates", "net-diversification rates", "relative-extinction rates")
#}
FIG_TYPES <- c("speciation rates", "net-diversification rates","relative-extinction rates","speciation shift times", "mass extinction Bayes factors","mass extinction times")
pdf(sprintf("%s.pdf",analysis_name))
layout.mat <- matrix(1:NUM_FIGS,nrow=2,ncol=NUM_FIGS / 2)
layout(layout.mat)
tess.plot.output(out2,fig.types=FIG_TYPES,las=1)
dev.off()

