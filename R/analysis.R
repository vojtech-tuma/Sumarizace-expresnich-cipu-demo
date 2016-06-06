#### Initialize ####
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(Cairo)
source("BEST.R")
source('exon_assignment_new.R')

############################################################################################
theme_set(theme_bw())

#### Effect size example ####
a <- rnorm(5000, 0, 1)
b <- rnorm(5000, 1, 1)
c <- rnorm(5000, 0, 3)
d <- rnorm(5000, 1, 3)
e <- as.data.frame(a)
e$b <- b
e$c <- c
e$d <- d
a_name <- ("\u03BC = 0, \u03BC = 1") ### needs manual replacement  
b_name <- ("\u03BC = 1, \u03BC = 1")
c_name <- ("\u03BC = 0, \u03BC = 3")
d_name <- ("\u03BC = 1, \u03BC = 3")
names(e) <- c(a_name,b_name,c_name,d_name)
e <- gather(e)
e$key <- head(replace(e$key, "<U+03C3>", "\u03C3"),n=-1)
e$group <- ifelse(e$key == a_name | e$key == b_name,
                  paste("velikost účinku = ", round((mean(a) - mean(b)) / (sqrt((sd(a)+sd(b))/2)),3)),
                  paste("velikost účinku = ", round((mean(c) - mean(d)) / (sqrt((sd(c)+sd(d))/2)),3))
            ) 
e$group_f <- factor(e$group, levels = c(paste("velikost účinku = ", round((mean(a) - mean(b)) / (sqrt((sd(a)+sd(b))/2)),3)),paste("velikost účinku = ", round((mean(c) - mean(d)) / (sqrt((sd(c)+sd(d))/2)),3))))

e %>% ggplot(aes(x=value, fill=key)) + geom_density(alpha=.5) + 
  xlab("hodnota")+
  ylab("hustota")+
  scale_fill_discrete(name="parametry\nrozdělení") +
  theme(legend.position = c(0.85, 0.8)) +
  scale_x_continuous(breaks=c(-10,-5,0,1,5,10)) +
  facet_wrap(~group_f,ncol=1) -> p_example
ggsave(file.path("results","effect_size_example.pdf"), p_example, width = 15, height = 15, units = "cm", device = cairo_pdf)

############################################################################################
#### Read preprocessed files (directory ./out) ####
cwd = getwd()
setwd(file.path(cwd,"out"))
ids<-read.table("probeset_exon_transcript", header=T)
ids_probesets<-unique(ids[,c("probeset_id","exon_cluster_id","transcript_cluster_id")])
snp<-unique(read.delim("probesets.snp.affected"))
snp_kill <- read.table("probes.kill", header=TRUE)
snp_random <- read.table("random.kill", header=TRUE)

mps_gene <- read.table('mps.gene', header = T) %>% 
  merge(snp_kill %>% 
          group_by(probeset_id) %>% 
          summarize(snp_kill=n_distinct(probe_id)),all.x = T) %>%
  merge(snp_random %>% group_by(probeset_id) %>%
          summarize(random_kill=n_distinct(probe_id)),all.x = T)


#### Read exons ####
core_snp<-read.exons("core_snp", genotypes=c("D", "M"), tissues=c("KID", "SPN", "TES"), methods=c("ITER_PLIER", "PLIER", "GCRMA"))
core_nosnp<-read.exons("core_nosnp", genotypes=c("D", "M"), tissues=c("KID", "SPN", "TES"), methods=c("ITER_PLIER", "PLIER", "GCRMA"))
core_random<-read.exons("core_random", genotypes=c("D", "M"), tissues=c("KID", "SPN", "TES"), methods=c("ITER_PLIER", "PLIER", "GCRMA"))

# apply log to PLIER summaries (function side-effect "<<-" )
lapply(names(core_snp), function(name){ if(length(grep("PLIER",name,value=TRUE)) > 0) {core_snp[[name]]$signal <<- log2(core_snp[[name]]$signal + 16); return(name)}})
lapply(names(core_nosnp), function(name){ if(length(grep("PLIER",name,value=TRUE)) > 0) {core_nosnp[[name]]$signal <<- log2(core_nosnp[[name]]$signal + 16); return(name)}})
lapply(names(core_random), function(name){ if(length(grep("PLIER",name,value=TRUE)) > 0) {core_random[[name]]$signal <<- log2(core_random[[name]]$signal + 16); return(name)}})

#### Read genes ####
mps_snp <- read.genes("core_mps_snp", genotypes=c("D", "M"), tissues=c("KID", "SPN", "TES"), methods=c("ITER_PLIER", "PLIER", "GCRMA"))
mps_nosnp <- read.genes("core_mps_nosnp", genotypes=c("D", "M"), tissues=c("KID", "SPN", "TES"), methods=c("ITER_PLIER", "PLIER", "GCRMA"))
mps_random <- read.genes("core_mps_random", genotypes=c("D", "M"), tissues=c("KID", "SPN", "TES"), methods=c("ITER_PLIER", "PLIER", "GCRMA"))

# apply log to PLIER summaries (function side-effect "<<-" )
lapply(names(mps_snp), function(name){ if(length(grep("PLIER",name,value=TRUE)) > 0) {mps_snp[[name]]$signal <<- log2(mps_snp[[name]]$signal + 16); return(name)}})
lapply(names(mps_nosnp), function(name){ if(length(grep("PLIER",name,value=TRUE)) > 0) {mps_nosnp[[name]]$signal <<- log2(mps_nosnp[[name]]$signal + 16); return(name)}})
lapply(names(mps_random), function(name){ if(length(grep("PLIER",name,value=TRUE)) > 0) {mps_random[[name]]$signal <<- log2(mps_random[[name]]$signal + 16); return(name)}})


setwd(cwd)

# save data for parallel analysis
save(core_snp, core_nosnp, core_random, file="core.RData")
save(mps_snp, mps_nosnp, mps_random, file="mps.RData")

############################################################################################
theme_set(theme_bw())

#### DABG analysis ####

#### + DABG density ####
signal_dabg <-rbind(core_snp$D_KID_PLIER %>% 
               filter(chip=="D1100KID.CEL") %>% 
               mutate(treat="original"),
             core_nosnp$D_KID_PLIER %>% 
               filter(chip=="D1100KID.CEL") %>% 
               mutate(treat="nosnp")) 
signal_dabg %>%
  ggplot(aes(dabg, colour=treat)) +
  geom_line(stat="density") +
  labs(x="DABG",y="hustota") +
  scale_color_discrete(name="ošetření",labels=c("bez SNP","původní"))+
  scale_y_log10() + 
  scale_x_continuous(breaks=c(0,0.05,0.25,0.50,0.75,1)) -> p_dabg_density_line
ggsave(file.path("results","dabg-density-line.svg"), p_dabg_density_line, width = 16, height = 12, units = "cm")

signal_dabg %>%
  filter(dabg < 0.05) %>%
  ggplot(aes(dabg, colour=treat)) +
  labs(x="DABG",y="hustota") +
  scale_color_discrete(name="ošetření",labels=c("bez SNP","původní"))+
  geom_line(stat="density") +
  scale_y_log10() -> p_dabg_density_005 
ggsave(file.path("results","dabg-density-005.svg"), p_dabg_density_005, width = 16, height = 12, units = "cm")

#### + DABG ANOVA ####

dabg_extract <- merge(merge(core_snp$D_KID_PLIER %>% 
                              filter(chip=="D1100KID.CEL") %>% 
                              select(probeset_id,dabg),
                            core_nosnp$D_KID_PLIER %>% 
                              filter(chip=="D1100KID.CEL") %>% 
                              select(probeset_id, dabg), by="probeset_id",suffixes = c(".o","")),
                      core_random$D_KID_PLIER %>% 
                        filter(chip=="D1100KID.CEL") %>% 
                        select(probeset_id,dabg), by="probeset_id",suffixes=c(".t",".r")
                ) %>% gather(treat,dabg,-probeset_id)

aov.dabg=aov(dabg ~ treat + Error(probeset_id/treat), data=dabg_extract)
summary(aov.dabg)

##### BEST mcmc #####
# example run of BEST mcmc:
# ### + D1100TES.CEL_GCRMA ####
# # 
# mcmc_D1100TES_GCRMA_SNP <- BESTarray(core_snp$D_KID_GCRMA, core_nosnp$D_KID_GCRMA, "D1100TES.CEL", "SNP", numSavedSteps = 2000)
# save(mcmc_D1100TES_GCRMA_SNP,file="D1100TES_GCRMA_SNP.RData")
# 
# mcmc_D1100TES_GCRMA_random <- BESTarray(core_snp$D_KID_GCRMA, core_random$D_KID_GCRMA, "D1100TES.CEL", "random", numSavedSteps = 2000)
# save(mcmc_D1100TES_GCRMA_random,file="D1100TES_GCRMA_random.RData")
# 
# #### + D1100KID.CEL_PLIER ####
# mcmc_D1100KID_PLIER_SNP <- BESTarray(mps_snp$D_KID_PLIER, mps_nosnp$D_KID_PLIER, "D1100KID.CEL", "SNP", numSavedSteps = 100000)
# save(mcmc_D1100KID_PLIER_SNP,file="D1100KID_PLIER_SNP.RData")
# 
# mcmc_D1100KID_PLIER_random <- BESTarray(core_snp$D_KID_PLIER, core_random$D_KID_PLIER, "D1100KID.CEL", "random", numSavedSteps = 2000, log=T)
# save(mcmc_D1100KID_PLIER_random,file="D1100KID_PLIER_random.RData")
# 
# #### + D1100KID.CEL_ITER_PLIER ####
# mcmc_D1100KID_ITER_PLIER_SNP <- BESTarray(mps_snp$D_KID_ITER_PLIER, mps_nosnp$D_KID_ITER_PLIER, "D1100KID.CEL", "SNP", numSavedSteps = 100000)
# save(mcmc_D1100KID_ITER_PLIER_SNP,file="D1100KID_ITER_PLIER_SNP.RData")
# 
# mcmc_D1100KID_ITER_PLIER_random <- BESTarray(mps_snp$D_KID_ITER_PLIER, mps_random$D_KID_ITER_PLIER, "D1100KID.CEL", "random", numSavedSteps = 100000, log=T)
# save(mcmc_D1100KID_ITER_PLIER_random,file="D1100KID_ITER_PLIER_random.RData")

#### + load data from files ####

cwd = getwd()
setwd(file.path(cwd,"out/BEST"))
core.mcmc.list <- list()
invisible(
  lapply(list.files(pattern="core_.*.RData"), 
       function(name) {
          tmp.env <- new.env()
          load(name,envir=tmp.env)
          # assign(mcmc.list[gsub("(.CEL|.RData)","",name)],tmp.env$mcmc, pos=.GlobalEnv)
          .GlobalEnv$core.mcmc.list[[gsub("ITER_PLIER","IterPLIER",gsub("(.CEL|.RData)","",name))]] <- tmp.env$mcmc
       }))
mps.mcmc.list <- list()
invisible(
  lapply(list.files(pattern="mps_.*.RData"), 
         function(name) {
           tmp.env <- new.env()
           load(name,envir=tmp.env)
           # assign(mcmc.list[gsub("(.CEL|.RData)","",name)],tmp.env$mcmc, pos=.GlobalEnv)
           .GlobalEnv$mps.mcmc.list[[gsub("ITER_PLIER","IterPLIER",gsub("(.CEL|.RData)","",name))]] <- tmp.env$mcmc
         }))
setwd(cwd)

#### GCRMA Analysis ####

## violin plot for method comparison ##
BEST_violin <- function(mcmc.list, level, met, param, limits, breaks, x){
  # extract results
  mcmc.list %>% 
    .[grep(paste(level,"_","(D|M)","[0-9]*","(KID|SPN|TES)","_",met,sep=""),names(.))] -> selected
  selected %>%  #select only some results defined by variables above
    lapply(function(x) { 
      x <- as.data.frame(x$mcmc);                           #extract mcmc table as data frame
      names(x) <- sub("\\[","",sub("\\]","",names(x)));     #remove "[]" characters from names for easier manipulation
      mutate(x, muDiff=mu1-mu2, sigmaDiff=sigma1-sigma2, effSz=(mu1-mu2)/sqrt((sigma1^2+sigma2^2)/2)) %>%   #compute parameter differences and effect size
        return
    }) -> results
  
  # form a single long- table
  lapply(names(results), function(name){
    df <- as.data.frame(results[[name]])
    case <- strsplit(name,"[0-9_]+")                        # split name to spread case description i.e. core_D1302KID_GCRMA_SNP => D | KID | GCRMA | SNP | D1302KID | D1302KID_GCRMA_SNP
    df$level <- rep(case[[1]][1],dim(df)[1])
    df$genome <- rep(case[[1]][2],dim(df)[1])
    df$tissue <- rep(case[[1]][3],dim(df)[1])
    df$method <- rep(case[[1]][4],dim(df)[1])
    df$treatment <- rep(case[[1]][5],dim(df)[1])
    df$chip <- rep(gsub("_.*", "", name))
    df$case <- rep(name,dim(df)[1])
    return(df)
  }) -> spreaded
  bind_rows(spreaded) -> test
  
  test %>%
    group_by(method,genome,tissue, treatment, chip) %>%
    mutate(med = median(get(param, pos=.))) -> test2
  
  # plot the result
  test2 %>%
    ggplot(aes(y=get(param),x=factor(interaction(method,treatment, tissue, genome)),color=factor(interaction(method,treatment)))) + 
    geom_violin() + 
    stat_boxplot(geom ='errorbar') + 
    geom_boxplot(width=.5, outlier.shape = NA) +
    theme(panel.grid.major.y = element_blank(), 
          legend.position="none",
          axis.text.x = element_text(angle=90, vjust=0.5))+
    ggtitle(paste("Porovnání parametrů BEST mezi ošetřeními pro", met,"-",x,sep = " ")) +
    labs(y=x,x="sumarizace")  -> p_BEST_param
  ggsave(file.path("results","BEST_summary",paste("BEST",level,met,param,".pdf",sep = "_")), p_BEST_param, width = 20, height = 25, units = "cm", device = cairo_pdf)
}

## violin plot for chip comparison ##
BEST_violin_chip <- function(mcmc.list, level, met, gen, tis, param, limits, breaks, x){
  # extract results
  mcmc.list %>% 
    .[grep(paste(level,"_",gen,"[0-9]*",tis,"_",met,sep=""),names(.))] -> selected
  selected %>%  #select only some results defined by variables above
    lapply(function(x) { 
      x <- as.data.frame(x$mcmc);                           #extract mcmc table as data frame
      names(x) <- sub("\\[","",sub("\\]","",names(x)));     #remove "[]" characters from names for easier manipulation
      mutate(x, muDiff=mu1-mu2, sigmaDiff=sigma1-sigma2, effSz=(mu1-mu2)/sqrt((sigma1^2+sigma2^2)/2)) %>%   #compute parameter differences and effect size
        return
    }) -> results
  
  # form a single long- table
  lapply(names(results), function(name){
    df <- as.data.frame(results[[name]])
    case <- strsplit(name,"[0-9_]+")                        # split name to spread case description i.e. core_D1302KID_GCRMA_SNP => core | D | KID | GCRMA | SNP | D1302KID | D1302KID_GCRMA_SNP
    df$level <- rep(case[[1]][1],dim(df)[1])
    df$genome <- rep(case[[1]][2],dim(df)[1])
    df$tissue <- rep(case[[1]][3],dim(df)[1])
    df$method <- rep(case[[1]][4],dim(df)[1])
    df$treatment <- rep(case[[1]][5],dim(df)[1])
    df$chip <- rep(gsub("[^0-9]","",name))
    df$case <- rep(name,dim(df)[1])
    return(df)
  }) -> spreaded
  bind_rows(spreaded) -> test
  
  test %>%
    group_by(method,genome,tissue, treatment) %>%
    mutate(med = median(get(param, pos=.))) -> test2
  
  # plot the result
  test2 %>%
    ggplot(aes(y=get(param),x=factor(interaction(treatment, chip)),color=factor(treatment))) + 
    geom_violin() + 
    stat_boxplot(geom ='errorbar') + 
    geom_boxplot(width=.5, outlier.shape = NA) +
    
    theme(panel.grid.major.y = element_blank(), 
          legend.position="none",
          axis.text.x = element_text(angle=90, vjust=0.5))+
    ggtitle(paste("Porovnání parametrů BEST mezi čipy pro", gen, tis, met, "-",x,sep = " ")) +
    labs(y=x,x="sumarizace")  -> p_BEST_param
  ggsave(file.path("results","BEST_summary",paste("BEST",level,gen,tis,met,param,".pdf",sep = "_")), p_BEST_param, width = 20, height = 25, units = "cm", device = cairo_pdf)
}

####################################################################
theme_set(theme_bw())
# create all the plots
for(level in c("mps","core")){ #
  
  if(level == "core") {mcmc.list <- core.mcmc.list}
  if(level == "mps") {mcmc.list <- mps.mcmc.list}
  
  for(met in c( "GCRMA","IterPLIER", "PLIER")){ # 
    
    for(param in c("sigmaDiff", "muDiff", "effSz")){#
      
      if(param == "sigmaDiff") {limits = c(-.02, 0.07); breaks=seq(-1,1,.01); x="rozdíl odchylek"}
      if(param == "muDiff") { limits = c(-.40, -.12); breaks=seq(-1,0,.02); x="rozdíl středních hodnot" }
      if(param == "effSz") {limits <- c(-.20, -.04); breaks=seq(-1,0,.02); x="velikost účinku"}
      
      
      
      for(gen in c("D", "M")){
        for(tis in c("KID", "SPN", "TES")){
          print(paste(level,met,gen,tis,param))
          BEST_violin_chip(mcmc.list,level,met,gen,tis,param,limits,breaks,x)
        }
      }
      BEST_violin(mcmc.list, level, met, param, limits, breaks, x)
    }
  }
}













