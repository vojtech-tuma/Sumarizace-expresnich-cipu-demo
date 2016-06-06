read.genes<-function(version, genotypes = c("D"), tissues = c("KID"), methods = c("PLIER")){
  data<-list()
  for (genotype in genotypes) {
    for (tissue in tissues) {
      print(paste("reading: ",file.path(version,"DABG",paste(genotype,tissue,"DABG",sep="_"),"summary.txt")))
      mps_dabg <-read.table(file.path(version,"DABG",paste(genotype,tissue,"DABG",sep="_"),"summary.txt"), header=T)
      mps_dabg_long <- (gather(mps_dabg,chip,dabg,c(ends_with(".CEL"))))
      for (method in methods) {
        print(paste("gene", version, genotype, tissue, method, sep=" "))
        mps_signal<-read.table(file.path(version,method,paste(genotype,tissue,method,sep="_"),"summary.txt"), header=T) %>% 
          gather(chip, signal, c(ends_with(".CEL"))) %>%
          merge(mps_gene, all.x=T)
        mps <- merge(mps_signal,mps_dabg_long, by=c("probeset_id","chip"))
        
        name<-paste(genotype, tissue, method, sep="_")
        
        data[[name]]<-mps
      }
    }
  }
  return(data)
}

read.exons<-function(version, genotypes = c("D"), tissues = c("KID"), methods = c("PLIER")){
  data<-list()
  kill <- merge(snp_kill %>% group_by(probeset_id) %>% summarise(snp_kill=n()), 
                snp_random %>% group_by(probeset_id) %>% summarise(random_kill=n()), 
                by="probeset_id", all=T)
  kill[is.na(kill)] <- 0
  for (genotype in genotypes) {
    for (tissue in tissues) {
      print(paste("reading: ",file.path(version,"DABG",paste(genotype,tissue,"DABG",sep="_"),"summary.txt")))
      dabg <-read.table(file.path(version,"DABG",paste(genotype,tissue,"DABG",sep="_"),"summary.txt"), header=T)
      dabg_long <- (gather(dabg,chip,dabg,c(ends_with(".CEL"))))
      for (method in methods) {
        print(paste("exon", version, genotype, tissue, method, sep=" "))
        signal <- read.table(file.path(version,method,paste(genotype,tissue,method,sep="_"),"summary.txt"), header=T) %>% 
                  merge(ids_probesets, by="probeset_id", all.x = TRUE) %>%
                  merge(snp, by="probeset_id", all.x = TRUE) %>%
                  merge(kill, by="probeset_id") %>%
                  gather(chip, signal, c(ends_with(".CEL"))) %>% 
                  merge(dabg_long, by=c("probeset_id","chip"))
        name<-paste(genotype, tissue, method, sep="_")
        
        data[[name]]<-signal
      }
    }
  }
  return(data)
}

exon.summary<-function(data, genotype = "D", logscale = F){
  if(genotype == "D"){
    if(logscale){
      data <- data %>% mutate(signal = log2(signal))
    }
    result <- data %>% 
              #filter(dabg<0.05) %>%
              group_by(exon_cluster_id) %>% 
                               summarise(transcript_cluster_id = max(transcript_cluster_id),
                                         mean=mean(signal), 
                                         median=median(signal), 
                                         sd=sd(signal), 
                                         var=var(signal), 
                                         dabg=max(dabg), 
                                         SNP_hits=sum(WSB_sum)/n(), 
                                         SNP_probes=sum(WSB_probes)/n(), 
                                         probes=sum(total)/n(),
                                         probesets=n_distinct(probeset_id),
                                         snp_kill=sum(snp_kill)/n(),
                                         random_kill=sum(random_kill)/n(),
                                         sd.mean=sd/mean)
  } else {
    result <- data %>% 
              #filter(dabg<0.05) %>%
              group_by(exon_cluster_id) %>% 
                               summarise(transcript_cluster_id = max(transcript_cluster_id),
                                         mean=mean(signal), 
                                         median=median(signal), 
                                         sd=sd(signal), 
                                         var=var(signal), 
                                         dabg=max(dabg), 
                                         SNP_hits=sum(PWK_sum)/n(), 
                                         SNP_probes=sum(PWK_probes)/n(), 
                                         probes=sum(total)/n(),
                                         probesets=n_distinct(probeset_id),
                                         snp_kill=sum(snp_kill)/n(),
                                         random_kill=sum(random_kill)/n(),
                                         sd.mean=sd/mean)
  }
  return(result)
}

#### BEST ####
BESTarray <- function(x,y,chip,treat,numSavedSteps=10000,log=F){
  if(treat=="SNP"){
    data1 <- x %>% filter(chip == chip) %>% filter(dabg < 0.05) %>% filter(snp_kill > 0) %>% select(signal) %>% as.list
    data2 <- y %>% filter(chip == chip) %>% filter(dabg < 0.05) %>% filter(snp_kill > 0) %>% select(signal) %>% as.list
  } else {
    data1 <- x %>% filter(chip == chip) %>% filter(dabg < 0.05) %>% filter(random_kill > 0) %>% select(signal) %>% as.list
    data2 <- y %>% filter(chip == chip) %>% filter(dabg < 0.05) %>% filter(random_kill > 0) %>% select(signal) %>% as.list
  }
  if(log){
    data1$signal <- log2(data1$signal)
    data2$signal <- log2(data2$signal)
  }
  mcmc = BESTmcmc(data1$signal,data2$signal,numSavedSteps = numSavedSteps)
  return(list(x=data1,y=data2,mcmc=mcmc))
}

BESTarray.plot <- function(result){
  BESTplot(result$x$signal, result$y$signal, result$mcmc, pairsPlot = F)
}

BESTsummary.extract <- function(mcmc.list){
   summaries <- sapply(mcmc.list, 
           function(result){ 
             ret<- as.data.frame(t(BESTsummary(result$x,result$y,result$mcmc)))
             ret$stat <- rownames(ret)
             return(ret %>% gather(parameter,value,-stat))
             },
           simplify = FALSE, 
           USE.NAMES = TRUE)
    names(summaries) <- sub("mcmc_","",names(mcmc.list)) %>% sub("ITER_PLIER", "IterPLIER",.)
    
    return(lapply(seq_along(summaries),function(i){summaries[[i]] <- summaries[[i]] %>% mutate(chip=names(summaries)[[i]])}) %>% 
             bind_rows() %>%
             separate(col = chip, into = c("level","genome","tissue","method","treatment"),sep = "[0-9_]+",remove=F, convert=T) %>%
             spread(stat,value)
           )
}
