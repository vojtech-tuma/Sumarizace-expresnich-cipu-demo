t_M <- exon.summary(core_nosnp$M_KID_GCRMA,"M", logscale=F)
o_M <- exon.summary(core_snp$M_KID_GCRMA,"M", logscale=F)
r_M <- exon.summary(core_random$M_KID_GCRMA,"M", logscale=F)

t_D <- exon.summary(core_nosnp$D_KID_GCRMA,"D", logscale=F)
o_D <- exon.summary(core_snp$D_KID_GCRMA,"D", logscale=F)
r_D <- exon.summary(core_random$D_KID_GCRMA,"D", logscale=F)

M <- merge(o_M,t_M,by="exon_cluster_id", suffixes = c(".o",""), all=F) %>% 
  merge(r_M,by="exon_cluster_id",suffixes = c(".t",".r"), all=F) %>%
  select(exon_cluster_id, starts_with("mean"), starts_with("median"), starts_with("sd"), starts_with("var"), starts_with("dabg"), 
         SNP_hits = SNP_hits.o, SNP_probes = SNP_probes.o, probes = probes.o, snp_kill=snp_kill.o, random_kill=random_kill.o,
         -SNP_hits.t,-SNP_probes.t,-random_kill.t,-snp_kill.t,-SNP_hits.r,-SNP_probes.r,-snp_kill.r,-random_kill.r) 
  

D <- merge(o_D,t_D,by="exon_cluster_id", suffixes = c(".o",""), all=F) %>% 
  merge(r_D,by="exon_cluster_id",suffixes = c(".t",".r"), all=F) %>%
  select(exon_cluster_id, starts_with("mean"), starts_with("median"), starts_with("sd"), starts_with("var"), starts_with("dabg"), 
         SNP_hits = SNP_hits.o, SNP_probes = SNP_probes.o, probes = probes.o, snp_kill=snp_kill.o, random_kill=random_kill.o,
         -SNP_hits.t,-SNP_probes.t,-random_kill.t,-snp_kill.t,-SNP_hits.r,-SNP_probes.r,-snp_kill.r,-random_kill.r) 

M %>%
  filter(SNP_hits > 0) %>%
  ggplot(aes(var.o,var.t, alpha=0.1)) + 
  geom_point(color='black',alpha=0.1,fill=NA, na.value=NA) +
  stat_smooth()

M %>%
  filter(SNP_hits > 0) %>%
  ggplot() +
  #  geom_histogram(aes(var.o, alpha=0.3), fill='red', binwidth=0.01) +
  #  geom_histogram(aes(var.t, alpha=0.3), fill='blue', binwidth=0.01) +
  geom_density(aes(var.o, alpha=0.3), fill='red') +
  geom_density(aes(var.t, alpha=0.3), fill='blue') +
  xlim(0,3)

M %>% 
  filter(SNP_hits > 0) %>%
  ggplot(aes(sd.mean.o,sd.mean.r)) +
  geom_point() +
  geom_smooth()+
  geom_abline()

M %>% 
  filter(snp_kill > 0 | random_kill > 0) %>%
  ggplot(aes(log2(sd.mean.t),log2(sd.mean.r))) +
  geom_point(alpha=0.01) +
  geom_smooth()+
  geom_abline() 

D %>%
  filter(snp_kill > 0) %>%
  ggplot(aes(sd.mean.o, sd.mean.t)) +
  geom_point(alpha=0.1) +
  geom_smooth(colour="red")

D %>%
  filter(random_kill > 0) %>%
  ggplot(aes(sd.mean.o, sd.mean.r)) +
  geom_point(alpha=0.1, fill='black') +
  geom_smooth(colour="red")


# o_M %>% 
#   #filter(probes <= 40) %>%
#   ggplot(aes(probes)) + 
#   geom_histogram(aes(fill=factor(probesets)), binwidth=1) + 
#   #xlab("prób v exonu") + 
#   #ylab("počet") + 
#   coord_trans(ytrans="sqrt") + 
#   scale_fill_discrete(name="probesetů") +
#   scale_x_continuous(breaks=seq(0,66,4)) + 
#   #scale_y_continuous(breaks=(c(1:400)*c(1:400))[seq(1,400,10)]) #TODO y scale
#   scale_y_continuous(breaks=c(1,50,200,500,1000,2000,5000,10000,20000,35000,55000,80000,110000))

#### Normal density comparison ####
M %>% filter(snp_kill > 0 | random_kill > 0) -> M_kill
M_kill %>% ggplot() + 
      stat_function(fun = dnorm, arg=list(mean = mean(M_kill$median.o),sd = sd(M_kill$median.o)),aes(colour = 'Normal')) + 
      stat_function(fun = dt, arg=list(df = 10^2.82, ncp=length(M_kill$mean.o)),aes(colour = 't')) + 
      geom_density(aes(x=median.o,colour='original')) +
      geom_density(aes(x=median.t,colour='treat'))

test <- merge(core_snp$D_KID_GCRMA, core_nosnp$D_KID_GCRMA %>% select(-PWK_probes,-PWK_sum,-WSB_probes,-WSB_sum,-total,-snp_kill,-random_kill), by=c("probeset_id","chip","exon_cluster_id","transcript_cluster_id"))
test2 <- test %>% group_by(exon_cluster_id) %>% do(t=t.test(.$signal.x,.$signal.y,paired=TRUE))

t.test(M$mean.o,M$mean.t,paired=T,alt="less")

##### probe - exon - gene ####
#### M ####
M_chip_exon <- core_snp$M_KID_GCRMA %>% 
  filter(chip=="M1006KID.CEL") %>% 
  group_by(exon_cluster_id) %>% 
  summarise(transcript_cluster_id = max(transcript_cluster_id),
            hits=sum(PWK_sum),
            probes=sum(total), 
            probesets=n_distinct(probeset_id),
            snp_kill=sum(snp_kill),
            random_kill=sum(random_kill),
            dabg=max(dabg))
M_chip_gene <- M_chip_exon %>%
  group_by(transcript_cluster_id) %>% 
  summarise(exons=n_distinct(exon_cluster_id), 
            probes=sum(probes), 
            probesets=sum(probesets), 
            snp_kill=sum(snp_kill), 
            random_kill=sum(random_kill), 
            dabg=max(dabg))
M_chip_exon %>% ggplot(aes(probes)) + geom_histogram(aes(fill=factor(probesets)),binwidth=1) + scale_y_sqrt()

#### D ####
D_chip_exon <- core_nosnp$D_KID_GCRMA %>% 
  filter(chip=="D1062KID.CEL") %>% 
  group_by(exon_cluster_id) %>% 
  summarise(transcript_cluster_id = max(transcript_cluster_id),
            hits=sum(WSB_sum),
            probes=sum(total), 
            probesets=n_distinct(probeset_id),
            snp_kill=sum(snp_kill),
            random_kill=sum(random_kill),
            dabg=max(dabg))
D_chip_gene <- D_chip_exon %>%
  group_by(transcript_cluster_id) %>% 
  summarise(exons=n_distinct(exon_cluster_id), 
            probes=sum(probes), 
            probesets=sum(probesets), 
            snp_kill=sum(snp_kill), 
            random_kill=sum(random_kill), 
            dabg=max(dabg))

#### SNP exploration ####
#### + SNP hits per exon
theme_set(theme_bw())
M_chip_exon %>% 
  filter(probes < 26) %>%
  ggplot(aes(probes, snp_kill)) + 
  geom_jitter(color='black',alpha=0.1, position=position_jitter(width=0.4, height=0.4), shape=20) +
  scale_x_continuous(breaks=seq(0,26,1)) +
  scale_y_continuous(breaks=seq(0,11,1)) +
  theme(panel.grid.major = element_line(color=NA),
  panel.grid.minor = element_line(color='gray')) +
  xlab("počet prób v exonu") +
  ylab("počet prób zasažených SNP") -> p_probe_hit_exon
ggsave(file.path("results","probe-hit-exon.png"), p_probe_hit_exon, width = 25, height = 15, units="cm" )

#### affected probes histogram ####
D_chip_exon %>%
  ggplot(aes(factor(hits))) +
  geom_histogram(binwidth=1) + 
  labs(x="WSB - počet zásahů SNP v exonu",y="četnost") +
  stat_bin(binwidth=1, origin = -0.7, geom="text", aes(label=..count..), vjust=-1, hjust=-.0001, size=2) +
  scale_y_sqrt(breaks=c(500,2000,5000,10000,20000,35000,60000,80000), limits = c(0,95000), expand = c(0,0)) +
  scale_x_discrete(breaks=c(0:17)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) -> p_wsb_snp
ggsave(file.path("results","wsb_snp.svg"), p_wsb_snp, width = 10, height = 10, units="cm" )

M_chip_exon %>%
  ggplot(aes(factor(hits))) +
  geom_histogram(binwidth=1) + 
  stat_bin(binwidth=1, origin = -0.7, geom="text", aes(label=..count..), vjust=-1, hjust=-.0001, size=2) +
  scale_y_sqrt(breaks=c(500,2000,5000,10000,20000,35000,60000,80000), limits = c(0,95000), expand = c(0,0)) +
  scale_x_discrete(breaks=c(0:17)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  xlab("PWK - počet zásahů SNP v exonu") +
  ylab("četnost") -> p_pwk_snp
ggsave(file.path("results","pwk_snp.svg"), p_pwk_snp, width = 10, height = 10, units="cm" )


#### #####
# plot probes per exon #
D_chip_exon %>% 
  filter(probes < 30) %>%
  ggplot(aes(factor(probes))) + 
  geom_histogram(aes(fill=factor(probesets)),binwidth=1)  + 
  coord_trans(ytrans="sqrt") + 
  labs(x="počet prób v exonu",y="četnost") +
  scale_fill_discrete(name="počet \nprobesetů") +
  scale_y_continuous(breaks=c(1,50,200,500,1000,2000,5000,10000,20000,33000,50000,70000,110000)) + 
  scale_x_discrete(breaks=seq(0,40,by=4)) -> p_probes_exon
ggsave(file.path("results","probes_per_exon.svg"), p_probes_exon, width = 20, height = 15, units="cm" )



# D_chip_gene %>%
#   ggplot(aes(factor(probes))) +
#   geom_histogram(aes(fill=factor(exons)), binwidth=1) + 
#   scale_y_continuous(breaks=c(1,10,50,200,500,1000))
D_chip_gene %>%
  filter(exons <= 30) %>%
  ggplot(aes(factor(exons))) +
  geom_histogram(aes(fill=cut(probes, c(1,2,5,10,20,50,100,200))), binwidth=1) + 
  labs(x="počet exonů v genu",y="četnost") +
  scale_fill_discrete(name="počet prób") +
  scale_y_continuous(breaks=seq(0,3000,300))+#c(10,100,300,600,1000,1400,1900,2400)) +
  scale_x_discrete(breaks=seq(0,40,by=4)) -> p_probes_exon_gene
ggsave(file.path("results","exons_per_gene.svg"), p_probes_exon_gene, width = 20, height = 15, units="cm" )



##### presentation ####
model_original <- 1
names(model_original) <- 'x'
model_original$x <- c(1,2,3,4)
model_original$y <- c(2.9,2.1,2.4,.4)

model_original <- data.frame(model_original$x, model_original$y)
names(model_original)<-c("x","y")
model_reduced <- model_original[-4,]
model_nohit <- model_original
model_nohit$y <- c(2.9,2.1,2.4,2.7)

theme_set(theme_bw())
p <- ggplot(model_nohit,aes(x, y, size=5)) + 
  geom_point(aes(size=5)) +
  #geom_point(data=model_reduced,aes(size=5)) + 
  geom_smooth(colour="blue", method="lm", se=FALSE, fullrange=TRUE,aes(size=5)) +
  geom_point(aes(4,2.7),color='green',shape=1) +
  xlab("próba") + ylab("signál") +
  coord_fixed(ratio = 1.3, xlim = NULL, ylim = c(0,3.3))
p + theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          legend.position="none",
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) -> p_snp_nohit
ggsave(file.path("results","snp_nohit.svg"), p_snp_nohit, width = 8, height = 10, units="cm" )


theme_set(theme_bw())
p <- ggplot(model_original,aes(x, y, size=5)) + 
  geom_point(aes(size=5)) +
  #geom_point(data=model_reduced,aes(size=5)) + 
  geom_smooth(colour="blue", method="lm", se=FALSE, fullrange=TRUE,aes(size=5)) +
  geom_point(aes(4,.4,color='green'),shape=1) +
  xlab("próba") + ylab("signál") +
  coord_fixed(ratio = 1.3, xlim = NULL, ylim = c(0,3.3))
p + theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          legend.position="none",
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) -> p_snp_original
ggsave(file.path("results","snp_original.svg"), p_snp_original, width = 8, height = 10, units="cm" )


theme_set(theme_bw())
p <- ggplot(model_original,aes(x, y, size=5)) + 
  geom_point(aes(size=5)) +
  #geom_point(data=model_reduced,aes(size=5)) + 
  geom_smooth(data=model_reduced,colour="blue", method="lm", se=FALSE, fullrange=TRUE,aes(size=5)) +
  geom_point(aes(4,.4,color='red'),shape=13) +
  xlab("próba") + ylab("signál") +
  coord_fixed(ratio = 1.3, xlim = NULL, ylim = c(0,3.3))
p + theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          legend.position="none",
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) -> p_snp_treated
ggsave(file.path("results","snp_treated.svg"), p_snp_treated, width = 8, height = 10, units="cm" )


#### DABG density ####
test <-rbind(core_snp$D_KID_PLIER %>% 
               filter(chip=="D1100KID.CEL") %>% 
               mutate(treat="original"),
             core_nosnp$D_KID_PLIER %>% 
               filter(chip=="D1100KID.CEL") %>% 
               mutate(treat="nosnp")) 
test %>%
  ggplot(aes(dabg, colour=treat)) +
  geom_line(stat="density") +
  scale_y_log10() + 
  scale_x_continuous(breaks=c(0,0.05,0.25,0.50,0.75,1)) -> p_dabg_density_line
ggsave(file.path("results","dabg-density-line.pdf"), p_dabg_density_line, width = 16, height = 12, units = "cm")

test %>%
  filter(dabg < 0.05) %>%
  ggplot(aes(dabg, colour=treat)) +
  geom_line(stat="density") +
  scale_y_log10() -> p_dabg_density_005 
ggsave(file.path("results","dabg-density-005.pdf"), p_dabg_density_005, width = 16, height = 12, units = "cm")

test2 <- merge(merge(core_snp$D_KID_PLIER %>% 
                 filter(chip=="D1100KID.CEL") %>% 
                 select(probeset_id,dabg),
               core_nosnp$D_KID_PLIER %>% 
                 filter(chip=="D1100KID.CEL") %>% 
                 select(probeset_id, dabg), by="probeset_id",suffixes = c(".o","")),
               core_random$D_KID_PLIER %>% 
                 filter(chip=="D1100KID.CEL") %>% 
                 select(probeset_id,dabg), by="probeset_id",suffixes=c(".t",".r")
              )

ks.test(test2$dabg.o,test2$dabg.t)

test2 %>%
  filter(dabg.o < 0.05 & dabg.t < 0.05) %>%
  mutate(diff=dabg.o-dabg.t) %>%
  ggplot(aes(dabg.o,dabg.t)) +
  geom_point() +
  geom_smooth()
  geom_line(stat="density") +
  scale_y_log10()
  

test <- gather(test2,treat,dabg,-probeset_id)
head(test)
aov.out=aov(dabg ~ treat + Error(probeset_id/treat), data=test)
summary(aov.out)
  
  
#### BEST #####
#### + M1006TES_GCRMA ####
mcmc_M1006TES_GCRMA_SNP <- BESTarray(core_snp$M_KID_GCRMA, core_nosnp$M_KID_GCRMA, "M1006TES.CEL", "SNP", numSavedSteps = 2000)
save(mcmc_M1006TES_GCRMA_SNP,file="M1006TES_GCRMA_SNP.RData")

mcmc_M1006TES_GCRMA_random <- BESTarray(core_snp$M_KID_GCRMA, core_random$M_KID_GCRMA, "M1006TES.CEL", "random", numSavedSteps = 2000)
save(mcmc_M1006TES_GCRMA_random,file="M1006TES_GCRMA_random.RData")

#### + M1024TES_GCRMA ####
mcmc_M1024TES_GCRMA_SNP <- BESTarray(core_snp$M_KID_GCRMA, core_nosnp$M_KID_GCRMA, "M1024TES.CEL", "SNP", numSavedSteps = 2000)
save(mcmc_M1024TES_GCRMA_SNP,file="M1024TES_GCRMA_SNP.RData")

mcmc_M1024TES_GCRMA_random <- BESTarray(core_snp$M_KID_GCRMA, core_random$M_KID_GCRMA, "M1024TES.CEL", "random", numSavedSteps = 2000)
save(mcmc_M1024TES_GCRMA_random,file="M1024TES_GCRMA_random.RData")

#### + D1062TES.CEL_GCRMA ####
mcmc_D1062TES_GCRMA_SNP <- BESTarray(core_snp$D_KID_GCRMA, core_nosnp$D_KID_GCRMA, "D1062TES.CEL", "SNP", numSavedSteps = 2000)
save(mcmc_D1062TES_GCRMA_SNP,file="D1062TES_GCRMA_SNP.RData")

mcmc_D1062TES_GCRMA_random <- BESTarray(core_snp$D_KID_GCRMA, core_random$D_KID_GCRMA, "D1062TES.CEL", "random", numSavedSteps = 2000)
save(mcmc_D1062TES_GCRMA_random,file="D1062TES_GCRMA_random.RData")

#### + D1100TES.CEL_GCRMA ####
mcmc_D1100TES_GCRMA_SNP <- BESTarray(core_snp$D_KID_GCRMA, core_nosnp$D_KID_GCRMA, "D1100TES.CEL", "SNP", numSavedSteps = 2000)
save(mcmc_D1100TES_GCRMA_SNP,file="D1100TES_GCRMA_SNP.RData")

mcmc_D1100TES_GCRMA_random <- BESTarray(core_snp$D_KID_GCRMA, core_random$D_KID_GCRMA, "D1100TES.CEL", "random", numSavedSteps = 2000)
save(mcmc_D1100TES_GCRMA_random,file="D1100TES_GCRMA_random.RData")

##### BEST plots #####
theme_set(theme_bw())
#### + mcmc.list ####
mcmc.list <- list(mcmc_D1100TES_GCRMA_random=mcmc_D1100TES_GCRMA_random,
                  mcmc_D1062TES_GCRMA_random=mcmc_D1062TES_GCRMA_random,
                  mcmc_D1100KID_GCRMA_random=mcmc_D1100KID_GCRMA_random,
                  mcmc_D1062KID_GCRMA_random=mcmc_D1062KID_GCRMA_random,
                  mcmc_D1100SPN_GCRMA_random=mcmc_D1100SPN_GCRMA_random,
                  #mcmc_D1062SPN_GCRMA_random=mcmc_D1062SPN_GCRMA_random,
                  mcmc_M1024TES_GCRMA_random=mcmc_M1024TES_GCRMA_random,
                  mcmc_M1006TES_GCRMA_random=mcmc_M1006TES_GCRMA_random,
                  mcmc_M1024KID_GCRMA_random=mcmc_M1024KID_GCRMA_random,
                  mcmc_M1006KID_GCRMA_random=mcmc_M1006KID_GCRMA_random,
                  mcmc_M1024SPN_GCRMA_random=mcmc_M1024SPN_GCRMA_random,
                  mcmc_M1006SPN_GCRMA_random=mcmc_M1006SPN_GCRMA_random,
                  mcmc_D1062KID_GCRMA_SNP=mcmc_D1062KID_GCRMA_SNP,
                  mcmc_D1100KID_GCRMA_SNP=mcmc_D1100KID_GCRMA_SNP,
                  mcmc_D1062TES_GCRMA_SNP=mcmc_D1062TES_GCRMA_SNP,
                  mcmc_D1100TES_GCRMA_SNP=mcmc_D1100TES_GCRMA_SNP,
                  #mcmc_D1062SPN_GCRMA_SNP=mcmc_D1062SPN_GCRMA_SNP,
                  mcmc_D1100SPN_GCRMA_SNP=mcmc_D1100SPN_GCRMA_SNP,
                  mcmc_M1006KID_GCRMA_SNP=mcmc_M1006KID_GCRMA_SNP,
                  mcmc_M1024KID_GCRMA_SNP=mcmc_M1024KID_GCRMA_SNP,
                  mcmc_M1006TES_GCRMA_SNP=mcmc_M1006TES_GCRMA_SNP,
                  mcmc_M1024TES_GCRMA_SNP=mcmc_M1024TES_GCRMA_SNP,
                  mcmc_M1006SPN_GCRMA_SNP=mcmc_M1006SPN_GCRMA_SNP,
                  mcmc_M1024SPN_GCRMA_SNP=mcmc_M1024SPN_GCRMA_SNP,
                  mcmc_D1100TES_PLIER_random=mcmc_D1100TES_PLIER_random,
                  mcmc_D1062TES_PLIER_random=mcmc_D1062TES_PLIER_random,
                  mcmc_D1100KID_PLIER_random=mcmc_D1100KID_PLIER_random,
                  mcmc_D1062KID_PLIER_random=mcmc_D1062KID_PLIER_random,
                  mcmc_D1100SPN_PLIER_random=mcmc_D1100SPN_PLIER_random,
                  #mcmc_D1062SPN_PLIER_random=#mcmc_D1062SPN_PLIER_random,
                  mcmc_M1024TES_PLIER_random=mcmc_M1024TES_PLIER_random,
                  mcmc_M1006TES_PLIER_random=mcmc_M1006TES_PLIER_random,
                  mcmc_M1024KID_PLIER_random=mcmc_M1024KID_PLIER_random,
                  mcmc_M1006KID_PLIER_random=mcmc_M1006KID_PLIER_random,
                  mcmc_M1024SPN_PLIER_random=mcmc_M1024SPN_PLIER_random,
                  mcmc_M1006SPN_PLIER_random=mcmc_M1006SPN_PLIER_random,
                  mcmc_D1062KID_PLIER_SNP=mcmc_D1062KID_PLIER_SNP,
                  mcmc_D1100KID_PLIER_SNP=mcmc_D1100KID_PLIER_SNP,
                  mcmc_D1062TES_PLIER_SNP=mcmc_D1062TES_PLIER_SNP,
                  mcmc_D1100TES_PLIER_SNP=mcmc_D1100TES_PLIER_SNP,
                  #mcmc_D1062SPN_PLIER_SNP=mcmc_D1062SPN_PLIER_SNP,
                  mcmc_D1100SPN_PLIER_SNP=mcmc_D1100SPN_PLIER_SNP,
                  mcmc_M1006KID_PLIER_SNP=mcmc_M1006KID_PLIER_SNP,
                  mcmc_M1024KID_PLIER_SNP=mcmc_M1024KID_PLIER_SNP,
                  mcmc_M1006TES_PLIER_SNP=mcmc_M1006TES_PLIER_SNP,
                  mcmc_M1024TES_PLIER_SNP=mcmc_M1024TES_PLIER_SNP,
                  mcmc_M1006SPN_PLIER_SNP=mcmc_M1006SPN_PLIER_SNP,
                  mcmc_M1024SPN_PLIER_SNP=mcmc_M1024SPN_PLIER_SNP,
                  mcmc_D1100TES_ITER_PLIER_random=mcmc_D1100TES_ITER_PLIER_random,
                  mcmc_D1062TES_ITER_PLIER_random=mcmc_D1062TES_ITER_PLIER_random,
                  mcmc_D1100KID_ITER_PLIER_random=mcmc_D1100KID_ITER_PLIER_random,
                  mcmc_D1062KID_ITER_PLIER_random=mcmc_D1062KID_ITER_PLIER_random,
                  mcmc_D1100SPN_ITER_PLIER_random=mcmc_D1100SPN_ITER_PLIER_random,
                  #mcmc_D1062SPN_ITER_PLIER_random=#mcmc_D1062SPN_ITER_PLIER_random,
                  mcmc_M1024TES_ITER_PLIER_random=mcmc_M1024TES_ITER_PLIER_random,
                  mcmc_M1006TES_ITER_PLIER_random=mcmc_M1006TES_ITER_PLIER_random,
                  mcmc_M1024KID_ITER_PLIER_random=mcmc_M1024KID_ITER_PLIER_random,
                  mcmc_M1006KID_ITER_PLIER_random=mcmc_M1006KID_ITER_PLIER_random,
                  mcmc_M1024SPN_ITER_PLIER_random=mcmc_M1024SPN_ITER_PLIER_random,
                  mcmc_M1006SPN_ITER_PLIER_random=mcmc_M1006SPN_ITER_PLIER_random,
                  mcmc_D1062KID_ITER_PLIER_SNP=mcmc_D1062KID_ITER_PLIER_SNP,
                  mcmc_D1100KID_ITER_PLIER_SNP=mcmc_D1100KID_ITER_PLIER_SNP,
                  mcmc_D1062TES_ITER_PLIER_SNP=mcmc_D1062TES_ITER_PLIER_SNP,
                  mcmc_D1100TES_ITER_PLIER_SNP=mcmc_D1100TES_ITER_PLIER_SNP,
                  #mcmc_D1062SPN_ITER_PLIER_SNP=mcmc_D1062SPN_ITER_PLIER_SNP,
                  mcmc_D1100SPN_ITER_PLIER_SNP=mcmc_D1100SPN_ITER_PLIER_SNP,
                  mcmc_M1006KID_ITER_PLIER_SNP=mcmc_M1006KID_ITER_PLIER_SNP,
                  mcmc_M1024KID_ITER_PLIER_SNP=mcmc_M1024KID_ITER_PLIER_SNP,
                  mcmc_M1006TES_ITER_PLIER_SNP=mcmc_M1006TES_ITER_PLIER_SNP,
                  mcmc_M1024TES_ITER_PLIER_SNP=mcmc_M1024TES_ITER_PLIER_SNP,
                  mcmc_M1006SPN_ITER_PLIER_SNP=mcmc_M1006SPN_ITER_PLIER_SNP,
                  mcmc_M1024SPN_ITER_PLIER_SNP=mcmc_M1024SPN_ITER_PLIER_SNP
              )


test <- as.data.frame(t(BESTsummary(mcmc_D1062KID_GCRMA_random$x,mcmc_D1062KID_GCRMA_random$y,mcmc_D1062KID_GCRMA_random$mcmc)))
test$stat <- rownames(test)
test <- test %>% gather(parameter,value,-stat)




#### + plot ####
summaries <- BESTsummary.extract(mcmc.list = mcmc.list)


met<-"GCRMA" # "PLIER" "IterPLIER"
#param <-"sigmaDiff" # "muDiff" "effSz"
param <-"muDiff"
#param <-"effSz"
summaries %>% 
  filter(parameter==param) %>%
  #filter(method==met) %>%
  group_by(method,genome,tissue,treatment) %>%
  summarize(low=max(HDIlow),hi=min(HDIhigh),value=mean(mean),chip=paste(unique(method),unique(genome),unique(tissue),unique(treatment))) %>%
  ggplot(aes(value,chip,color=factor(interaction(method,genome,tissue),levels=sort(levels(interaction(method,genome,tissue)))))) + 
  #geom_point() +
  geom_errorbarh(aes(xmin=low,xmax=hi,y=chip), height = .2) +
  geom_text(aes(x=value,y=chip,label=paste(genome, tissue, method, treatment)), vjust=-.7) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position="none", 
        axis.ticks.y = element_blank(), 
        axis.text.y=element_blank())+
  labs(x=param,y="sumarizace")
  #labs(x="síla efektu",y="sumarizace")         -> p_BEST_param
  #scale_x_continuous(limits = c(-.23, -.03),breaks=seq(-10,10,.02)) -> p_BEST_param
ggsave(file.path("results",paste("BEST_",param,"_",met,".svg",sep = "")), p_BEST_param, width = 25, height = 20, units = "cm")


###################################################################
random_mu <- cbind(
  arrange(data.frame(mu.x=mcmc_M1024TES_GCRMA_random$mcmc[,"mu[1]"] - mcmc_M1024TES_GCRMA_random$mcmc[,"mu[2]"]),mu.x),
  arrange(data.frame(mu.y=mcmc_M1006TES_GCRMA_random$mcmc[,"mu[1]"] - mcmc_M1006TES_GCRMA_random$mcmc[,"mu[2]"]),mu.y)
)
random_mu %>% ggplot(aes(mu.x,mu.y)) + geom_point()

snp_mu <- cbind(
  arrange(data.frame(mu.x=mcmc_M1024TES_GCRMA_SNP$mcmc[,"mu[1]"] - mcmc_M1024TES_GCRMA_SNP$mcmc[,"mu[2]"]),mu.x),
  arrange(data.frame(mu.y=mcmc_M1006TES_GCRMA_SNP$mcmc[,"mu[1]"] - mcmc_M1006TES_GCRMA_SNP$mcmc[,"mu[2]"]),mu.y)
)
snp_mu %>% ggplot(aes(mu.x,mu.y)) + geom_point()

