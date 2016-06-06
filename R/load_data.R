#### Initialize ####
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(magrittr)

source('exon_assignment.R')

#### Read common files ####
cwd = getwd()
setwd(file.path(cwd,"out"))
ids<-read.table("probeset_exon_transcript", header=T)
ids_probesets<-unique(ids[,c("probeset_id","exon_cluster_id","transcript_cluster_id")])
snp<-unique(read.delim("probesets.snp.affected"))
snp_kill <- read.table("B:\\Diplomka\\release\\out\\probes.kill", header=TRUE)
snp_random <- read.table("B:\\Diplomka\\release\\out\\random.kill", header=TRUE)


#### Read exons ####
core_snp<-read.exons("core_snp", genotypes=c("D", "M"), tissues=c("KID"), methods=c("PLIER", "GCRMA"))
core_nosnp<-read.exons("core_nosnp", genotypes=c("D", "M"), tissues=c("KID"), methods=c("PLIER", "GCRMA")) #, "SPN", "TES"
core_random<-read.exons("core_random", genotypes=c("D", "M"), tissues=c("KID"), methods=c("PLIER", "GCRMA"))

#### Read genes ####
mps_snp <- read.genes("core_mps_snp", genotypes=c("D", "M"), tissues=c("KID"), methods=c("ITER_PLIER", "PLIER", "GCRMA"))
mps_nosnp <- read.genes("core_mps_nosnp", genotypes=c("D", "M"), tissues=c("KID"), methods=c("ITER_PLIER", "PLIER", "GCRMA"))
mps_random <- read.genes("core_mps_random", genotypes=c("D", "M"), tissues=c("KID"), methods=c("ITER_PLIER", "PLIER", "GCRMA"))

#### Prepare and tidy ####
setwd(cwd)

# treated <- core_nosnp$D_KID_PLIER
# original <- core_snp$D_KID_PLIER
# random <- core_random$D_KID_PLIER
# 
# dmt_signal <- treated %>% select(-median,-ends_with(".dabg")) %>% gather(chip, signal, c(ends_with(".CEL"))) 
# 
# dmt_dabg <- treated %>% select(-median,-ends_with(".CEL"))
# names(dmt_dabg) <- gsub(".dabg","",names(dmt_dabg))
# dmt_dabg <- dmt_dabg %>% gather(chip, dabg, ends_with(".CEL")) 
# 
# dmt <- merge(dmt_signal, dmt_dabg)
# 
# 
# dmo_signal <- original %>% select(-median,-ends_with(".dabg")) %>% gather(chip, signal, c(ends_with(".CEL"))) 
# 
# dmo_dabg <- original %>% select(-median,-ends_with(".CEL"))
# names(dmo_dabg) <- gsub(".dabg","",names(dmo_dabg))
# dmo_dabg <- dmo_dabg %>% gather(chip, dabg, ends_with(".CEL")) 
# 
# dmo <- merge(dmo_signal, dmo_dabg)
# 
# dmr_signal <- random %>% select(-median,-ends_with(".dabg")) %>% gather(chip, signal, c(ends_with(".CEL"))) 
# 
# dmr_dabg <- random %>% select(-median,-ends_with(".CEL"))
# names(dmr_dabg) <- gsub(".dabg","",names(dmr_dabg))
# dmr_dabg <- dmr_dabg %>% gather(chip, dabg, ends_with(".CEL")) 
# 
# dmr <- merge(dmr_signal, dmr_dabg)
# 
# 
# dm <- bind_rows(data.frame(treat='nosnp', dmt), data.frame(treat='original', dmo))