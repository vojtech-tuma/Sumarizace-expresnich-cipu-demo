# Rscript mcmc.R core D KID GCRMA D1100KID.CEL snp
message("source functions")
library(dplyr)
source("BEST.R")
source('exon_assignment_new.R')

message("parse args")
args<-commandArgs(TRUE)
level <- args[2]
genotype <- args[3]
tissue <- args[4]
method <- args[5]
summary <- paste(genotype, tissue, method, sep="_")

data <- args[6]
treatment <- args[7]

message("run BEST")

if(level == "core"){
	load("core.RData")
	if(treatment == "SNP"){
		mcmc <- BESTarray(core_snp[[summary]], core_nosnp[[summary]], data, "SNP", numSavedSteps = 5000)
	} else if(treatment == "random"){
		mcmc <- BESTarray(core_snp[[summary]], core_random[[summary]], data, "random", numSavedSteps = 5000)
	} else {
		stop("Wrong treatment specified")
	}

} else if(level == "mps") {
	load("mps.RData")
	if(treatment == "SNP"){
		mcmc <- BESTarray(mps_snp[[summary]], mps_nosnp[[summary]], data, "SNP", numSavedSteps = 100000)
	} else if(treatment == "random"){
		mcmc <- BESTarray(mps_snp[[summary]], mps_random[[summary]], data, "random", numSavedSteps = 100000)
	} else {
		stop("Wrong treatment specified")
	}

} else {
	stop("Wrong data specified")
}



message("saving data")
cwd <- getwd()
setwd(file.path(cwd,"results"))
save(mcmc, file=paste(paste(level,data,method,treatment,sep="_"),"RData",sep="."))
setwd(cwd)
message("done")