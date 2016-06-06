#### Initialization #####
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(magrittr)
library(scales)

#library(extrafont)


#### + font import ####
font_import()
#library(grid)
#library(gridExtra)

#### + set basic variables ####
treated <- core_nosnp$D_KID_PLIER
original <- core_snp$D_KID_PLIER
random <- core_random$D_KID_PLIER
#output <- "D_KID_PLIER_random"
output <- "D_KID_PLIER"

#treated <- core_nosnp$M_KID_PLIER
#original <- core_snp$M_KID_PLIER
#output <- "M_KID_PLIER"


dir.create(file.path("results","exon",output), recursive=TRUE)

#theme_set(theme_bw() + theme(text=element_text(size = 12, 
#                                               family = "Arial", 
#                                               face = "plain", 
#                                               colour = "black", 
#                                               hjust = 0.5, 
#                                               vjust = 0.5,
#                                               angle = 0, 
#                                               lineheight = 0.9)))

#theme_set(theme_bw(base_family = "Arial"))
theme_set(theme_bw())

names(original)[names(original)=="WSB_probes"] <- "SNP_probes"
names(original)[names(original)=="PWK_probes"] <- "SNP_probes"
names(original)[names(original)=="WSB_sum"] <- "SNP_sum"
names(original)[names(original)=="PWK_sum"] <- "SNP_sum"
names(treated)[names(treated)=="WSB_probes"] <- "SNP_probes"
names(treated)[names(treated)=="PWK_probes"] <- "SNP_probes"
names(treated)[names(treated)=="WSB_sum"] <- "SNP_sum"
names(treated)[names(treated)=="PWK_sum"] <- "SNP_sum"
names(random)[names(random)=="WSB_probes"] <- "SNP_probes"
names(random)[names(random)=="PWK_probes"] <- "SNP_probes"
names(random)[names(random)=="WSB_sum"] <- "SNP_sum"
names(random)[names(random)=="PWK_sum"] <- "SNP_sum"


#### + Reshape the data into tidy format ####
#### + + treated ####
dmt_signal <- treated %>% select(-median,-ends_with(".dabg")) %>% gather(chip, signal, c(ends_with(".CEL"))) 

dmt_dabg <- treated %>% select(-median,-ends_with(".CEL"))
names(dmt_dabg) <- gsub(".dabg","",names(dmt_dabg))
dmt_dabg <- dmt_dabg %>% gather(chip, dabg, ends_with(".CEL")) 

dmt <- merge(dmt_signal, dmt_dabg)

#### + + original ####
dmo_signal <- original %>% select(-median,-ends_with(".dabg")) %>% gather(chip, signal, c(ends_with(".CEL"))) 

dmo_dabg <- original %>% select(-median,-ends_with(".CEL"))
names(dmo_dabg) <- gsub(".dabg","",names(dmo_dabg))
dmo_dabg <- dmo_dabg %>% gather(chip, dabg, ends_with(".CEL")) 

dmo <- merge(dmo_signal, dmo_dabg)

#### + + random ####
dmr_signal <- random %>% select(-median,-ends_with(".dabg")) %>% gather(chip, signal, c(ends_with(".CEL"))) 

dmr_dabg <- random %>% select(-median,-ends_with(".CEL"))
names(dmr_dabg) <- gsub(".dabg","",names(dmr_dabg))
dmr_dabg <- dmr_dabg %>% gather(chip, dabg, ends_with(".CEL")) 

dmr <- merge(dmr_signal, dmr_dabg)


#### + + bind together ####
dm <- bind_rows(data.frame(treat='nosnp', dmt), data.frame(treat='original', dmo), data.frame(treat='random', dmr))


#### + List common exons ####
common_exons <- intersect(dmt$exon_cluster_id, dmo$exon_cluster_id)
common_exons_random <- intersect(common_exons, dmr$exon_cluster_id)

#### + List common exons hit by snp ####
dmo %>% 
  filter(SNP_probes > 0) %>%
  filter(exon_cluster_id %in% common_exons) -> dmo_filter

common_hits <- intersect(dmo_filter$exon_cluster_id, dmt$exon_cluster_id)
common_hits_random <- intersect(common_hits, dmr$exon_cluster_id)

#### DABG exploration ####
#### + dabg_signal variable####
merge(dmt %>% select(exon_cluster_id, chip, signal,dabg), 
      dmo %>% select(exon_cluster_id, chip,signal,dabg), 
      by=c("exon_cluster_id", "chip"), suffixes = c(".nosnp",".original")) %>%
  rowwise() %>%
  mutate(dabg = max(dabg.original,dabg.nosnp)) %>% 
  mutate(diff = signal.original - signal.nosnp) %>% 
  mutate(diff_log = log2(signal.original) - log2(signal.nosnp)) -> dabg_signal

dmr %>% 
  select(exon_cluster_id, chip, signal, dabg) %>%
  merge(dabg_signal, by=c("exon_cluster_id", "chip"), suffixes = c(".random","")) -> dabg_signal_random
names(dabg_signal_random)[names(dabg_signal_random) == "signal"] <- "signal.random"

# now some simple plots
#
# question: is dabg different between the treatments

#### + dabg density ####
#dm %>%
#  ggplot(aes(dabg, fill=treat)) +
#  geom_density(colour=NA, alpha=0.3) +  
#  scale_y_log10()
#ggsave(file.path("results",output,"dabg-density.pdf"))

dm %>%
  ggplot(aes(dabg, colour=treat)) +
  geom_line(stat="density") +
  scale_y_log10() + 
  scale_x_continuous(breaks=c(0,0.05,0.25,0.50,0.75,1)) -> p_dabg_density_line
ggsave(file.path("results",output,"dabg-density-line.pdf"), p_dabg_density_line, width = 16, height = 12, units = "cm")

dm %>%
  filter(dabg < 0.05) %>%
  ggplot(aes(dabg, colour=treat)) +
  geom_line(stat="density") +
  scale_y_log10() -> p_dabg_density_005 
ggsave(file.path("results",output,"dabg-density-005.pdf"), p_dabg_density_005, width = 16, height = 12, units = "cm")

#grid.arrange(p_dabg_density_line, p_dabg_density_005, ncol = 2)
#arrangeGrob(p_dabg_density_line, p_dabg_density_line, ncol = 2) -> p_dabg_density_multiplot
#multiplot(p_dabg_density_line, p_dabg_density_005, cols=2) 
#ggsave(file.path("results",output,"dabg-density-multiplot.pdf"), p_dabg_density_multiplot, width = 30, height = 15, units = "cm")


#### + dabg plot ####
dabg_signal %>%
  ggplot(aes(dabg.original, dabg.nosnp, colour=chip)) + 
  geom_point( aes(alpha = 1 - dabg)) -> p_dabg_plot
ggsave(file.path("results",output,"dabg-plot.png"), p_dabg_plot, width = 40, height = 30, units = "cm")

dabg_signal %>% 
  #filter(dabg.original <= 0.05) %>%
  #filter(dabg.nosnp <= 0.05) %>%
  filter(dabg < 0.05) %>%
  ggplot(aes(dabg.original, dabg.nosnp, colour=chip)) + 
  geom_point(aes(alpha = 1 - dabg)) -> p_dabg_005
ggsave(file.path("results",output,"dabg-plot-005.png"), p_dabg_005, width = 40, height = 30, units = "cm")

dabg_signal %>%
  #filter(dabg.original <= 0.05) %>%
  #filter(dabg.nosnp <= 0.05) %>%
  filter(dabg < 0.05) %>%
  ggplot(aes(dabg.original, dabg.nosnp)) + 
  stat_bin2d(binwidth=c(0.005,0.005)) + 
  scale_fill_gradient2(breaks=c(100,1000,10000,20000,40000,80000,240000)) # -> p_dabg_bin
ggsave(file.path("results",output,"dabg-bin.png"), p_dabg_bin, width = 40, height = 30, units = "cm")

#### + dabg qqplot ####
# http://stats.stackexchange.com/a/12394

as.data.frame(qqplot(dmo$dabg, dmt$dabg, plot.it=FALSE)) %>%
  ggplot() + 
  geom_point(aes(x=x, y=y)) +
  xlab("dabg orignal") +
  ylab("dabg no snp") -> p_dabg_qqplot
ggsave(file.path("results",output,"dabg-qqplot.png"), p_dabg_qqplot, width = 20, height = 15, units = "cm")

#### Signal exploration ####
# comparison of expresion values
# question: Does the method change the overall distribution?

#### + signal qqplot < 0.05 ####
dmt %>% filter(dabg<0.05) -> dmt_005
dmo %>% filter(dabg<0.05) -> dmo_005
as.data.frame(qqplot(dmo_005$signal, dmt_005$signal, plot.it=FALSE)) %>%
  ggplot() + 
  geom_point(aes(x=x, y=y)) +
  xlab("signal orignal (DABG < 0.05)") +
  ylab("signal no snp (DABG < 0.05)") -> p_signal_qqplot
ggsave(file.path("results",output,"signal-qqplot.pdf"), p_signal_qqplot, width = 20, height = 15, units = "cm")


#### + signal boxplot ####
# there is something weird, try boxplot to see the outliers
# looks like an outlier issue in the qqplot
dm %>%
  filter(dabg < 0.05) %>%
  ggplot(aes(treat, signal)) +
  geom_boxplot() + 
  scale_y_log10() +
  xlab("treatment") -> p_signal_boxplot
ggsave(file.path("results",output,"signal-boxplot.pdf"), p_signal_boxplot, width = 20, height = 15, units = "cm")

#### + signal plot ####
dabg_signal %>% 
  #filter(dabg.original <= 0.05) %>%
  #filter(dabg.nosnp <= 0.05) %>%
  filter(dabg < 0.05) %>%
  ggplot(aes(log2(signal.original), log2(signal.nosnp))) + 
  #geom_smooth() +
  geom_point(alpha = 0.005) -> p_signal_005
ggsave(file.path("results",output,"signal-plot-005.png"), p_signal_005, width = 40, height = 30, units = "cm")

dabg_signal %>% 
  filter(exon_cluster_id %in% common_hits) %>%
  filter(dabg.original <= 0.05) %>%
  filter(dabg.nosnp <= 0.05) %>%
  ggplot(aes(log2(signal.original), log2(signal.nosnp))) + 
  #geom_smooth() +
  geom_point(alpha = 0.05) -> p_signal_hits
ggsave(file.path("results",output,"signal-plot-hits.png"), p_signal_hits, width = 40, height = 30, units = "cm")


#### + signal density ####
# just to be sure check density plots
# - looks like the original data have some outlier at the very end
dm %>%
  filter(dabg < 0.05) %>%
  ggplot(aes(log2(signal), colour=treat)) +
  geom_line(stat="density",alpha=0.7) -> p_signal_005_density
ggsave(file.path("results",output,"signal-005-density.pdf"), p_signal_005_density, width = 20, height = 10, units = "cm")

#### Exon exploration ####
#### + probes per exon ####
original %>% 
  filter(total <= 40) %>%
  ggplot(aes(total)) + 
  geom_histogram(aes(fill=factor(probesets)), binwidth=1) + 
  xlab("prób v exonu") + 
  ylab("počet") + 
  coord_trans(ytrans="sqrt") + 
  scale_fill_discrete(name="probesetů") +
  scale_x_continuous(breaks=seq(0,66,4)) + 
  #scale_y_continuous(breaks=(c(1:400)*c(1:400))[seq(1,400,10)]) #TODO y scale
  scale_y_continuous(breaks=c(1,50,200,500,1000,2000,5000,10000,20000,35000,55000,80000,110000)) -> p_probes_per_exon_original
ggsave(file.path("results",output,"probes_exon-hist-original.pdf"), p_probes_per_exon_original, width = 20, height = 15, units = "cm")

treated %>%
  filter(total <= 40) %>%
  ggplot(aes(total)) + 
  geom_histogram(aes(fill=factor(probesets)), binwidth=1) + 
  xlab("prób v exonu") + 
  ylab("počet") + 
  coord_trans(ytrans="sqrt") + 
  scale_fill_discrete(name="probesetů") +
  scale_x_continuous(breaks=seq(0,66,4)) +
  scale_y_continuous(breaks=c(1,50,200,500,1000,2000,5000,10000,20000,35000,55000,80000,110000)) -> p_probes_per_exon_treated
ggsave(file.path("results",output,"probes_exon-hist-trated.pdf"), p_probes_per_exon_treated, width = 20, height = 15, units = "cm")

#### + + aggregate counts of exons and probes per gene ####

original %>% 
  select(exon_cluster_id) %>% 
  aggregate(by=list(transcript_cluster_id=original$transcript_cluster_id), length) -> exon_count

original %>% 
  select(total) %>% 
  aggregate(by=list(transcript_cluster_id=original$transcript_cluster_id), sum) -> probe_count

counts <- merge(exon_count, probe_count)
names(counts) <- c("transcript_cluster_id", "exons", "probes")

#### + Probes per gene ####
counts %>%
  filter(exons <= 35) %>%
  ggplot(aes(factor(exons), probes)) +
  #geom_point(position=position_jitter(w = 0.3, h=0.1)) + 
  #geom_boxplot(aes(fill=probes)) +
  geom_violin(aes(fill=probes)) +
  scale_y_continuous(breaks=seq(0,300,10)) +
  xlab("exonů v genu") +
  ylab("prób v genu") -> p_probes_exon_gene
ggsave(file.path("results",output,"probes_exon_gene.pdf"), p_probes_exon_gene, width = 20, height = 15, units = "cm")

  
  

# now go one level deeper - check only the affected exon clusters
# question: How many probes are the affected?
# the 'treated' should contain correct data on probe hits by snps
#treated %>% 
#  mutate(hit_probes = PWK_probes + WSB_probes, any_hit = hit_probes > 0) %>%
#  group_by(any_hit) %>%
#  summarise(n=n())

# looks like 1/3 of the exon clusters has at least one hit
#any_hit     n
#1   FALSE 65406
#2    TRUE 31921

# now check the 'incidence' 
# question: How is the total number of probes in exon cluster related
#           to the number of probes with hits.
# one would expect some pretty good correlation
#  a chance and biological preference for mutations in some genes
#  being parts of the model of the probablity
#treated %>% 
#  mutate(hit_probes = PWK_probes + WSB_probes) %>%
#  ggplot(aes(probes, hit_probes)) +
#  geom_point()

#### SNP exploration ####
#### + SNP hits per exon
original %>% 
  filter(total < 35) %>%
  ggplot(aes(total, SNP_probes)) + 
  geom_jitter(alpha=0.1, position=position_jitter(width=0.3, height=0.3)) +
  xlab("probes per exon") +
  ylab("probes hit by SNP ") -> p_probe_hit_exon
ggsave(file.path("results",output,"probe-hit-exon.png"), p_probe_hit_exon, width = 20, height = 10, units="cm" )
  
  
#### + SNP histogram ####
core_snp$D_KID_PLIER %>%
  ggplot(aes(WSB_probes)) +
  geom_histogram(binwidth=1) + 
  xlab("WSB probes affected") +
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1, hjust=-.0001, size=2) +
  scale_y_sqrt(breaks=c(10,50,200,500,1000,2000,5000,10000,20000,50000,10000)) + 
  scale_x_continuous(breaks=c(0:12)) -> p_wsb_snp
ggsave(file.path("results",output,"wsb_snp.pdf"), p_wsb_snp, width = 10, height = 10, units="cm" )

core_snp$M_KID_PLIER %>%
  ggplot(aes(PWK_probes)) +
  geom_bar(binwidth=1) + 
  xlab("PWK probes affected") +
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1, hjust=-.0001, size=2) +
  scale_x_continuous(breaks=c(0:12)) +
  scale_y_sqrt(breaks=c(10,50,200,500,1000,2000,5000,10000,20000,50000,10000)) -> p_pwk_snp
ggsave(file.path("results",output,"pwk_snp.pdf"), p_pwk_snp, width = 10, height = 10, units="cm" )




# now check only the affected probes
#dm %>%
#  filter(PWK_probes + WSB_probes > 0) %>%
#  ggplot(aes(treat, signal)) +
#  geom_boxplot() +
#  scale_y_log10()

#### Affected probes only ####
#### + signal boxplot ####
dm %>%
  filter(exon_cluster_id %in% common_hits) %>%
  ggplot(aes(treat, signal)) +
  geom_boxplot() +
  scale_y_continuous(trans=log2_trans()) -> p_signal_hits_boxplot_log2
ggsave(file.path("results",output,"affected-probes-boxplot.pdf"),p_signal_hits_boxplot_log2, width=20, height=10, units="cm")
  
# find a reasonable limit for the outliers
quantile(dm$signal, .99)
# 6210

#### + signal density ####
dm %>%
  filter(exon_cluster_id %in% common_hits) %>%
  filter(dabg < 0.05) %>%
  ggplot(aes(log2(signal), colour=treat)) +
  geom_line(stat="density") -> p_signal_hits_005_density
  #scale_y_continuous(trans=log2_trans()) + 
  #xlim(0, 6243)
ggsave(file.path("results",output,"density-hits.pdf"),p_signal_hits_005_density, width=20, height=10, units="cm")

# in general the distribution does not seem to change much
### ... or does it? Slight increase in signal for treated data


#### Pairwise comparison ####
# we've got paired data - two treatments of one chip
# so we can check the difference

### old: dsigdabg (dabg_signal) generation

# check the differenes in general
# question: Is there any difference in the pairwise values?
# - there seems to be a difference
# question: What is the trend in the differences?
# - it is original minus treated, and the mean is ~ -20.64477 (!)# -7.2
#   the treated values tend to be higher, that's expectable, because 
#   the signal is not spoiled by badly matching probes in the set

#### Hits only -- dsigdabg var (dabg_signal %in% common_hits) ####
dsigdabg <- dabg_signal %>% filter(exon_cluster_id %in% common_hits)

#### + diff histogram ####
dsigdabg %>%
  filter(abs(diff) < 500) %>%
  ggplot(aes(diff)) +
  geom_histogram(aes(alpha=1),binwidth=1) +
  geom_histogram(aes(alpha=0.3),binwidth=2) +
  geom_histogram(aes(alpha=0.1),binwidth=10) -> p_diff_hist
ggsave(file.path("results",output,"difference-hist.pdf"), p_diff_hist, width = 20, height = 10, units = "cm")

#### + diff histogram with dabg fill ####
dsigdabg %>%
  filter(abs(diff) < 500) %>%
  mutate(dabg_treshold = dabg >= 0.05) %>%
  ggplot(aes(diff)) +
  #geom_histogram(aes(alpha=1, fill=factor(dabg_treshold)), binwidth=1) +
  geom_histogram(aes(alpha=0.3, fill=factor(dabg_treshold)), binwidth=2) + 
  geom_histogram(aes(alpha=0.1, fill=factor(dabg_treshold)), binwidth=10) + 
  scale_fill_discrete(name="dabg",labels=c("< 0.05", ">= 0.05")) -> p_diff_hist_dabg
ggsave(file.path("results",output,"diff-hist-dabg.pdf"),p_diff_hist_dabg, width=20, height=10, units="cm")


# find the mean without the outliers
dsigdabg %>% 
  filter(abs(diff) < 1000) %>% 
  {mean(.$diff)}

#### + binplot ratio - diff ####
library(hexbin)
dsigdabg %>%
  filter(dabg < 0.05) %>%
  mutate(ratio=(log2(signal.nosnp) / log2(signal.original))) %>%
  filter(abs(ratio) < 5) %>%
  filter(abs(diff_log) < 7) %>%
  ggplot(aes((ratio), diff_log)) + 
  stat_binhex(binwidth=c(0.033,0.1)) +
  scale_fill_gradient(trans="log10") +
  scale_x_continuous(breaks=-10:15) -> p_ratio_diff
  #geom_point(aes(alpha=0.2)) 
ggsave(file.path("results",output,"ratio_diff-binplot.pdf"), p_ratio_diff, width = 20, height = 10, units = "cm")

# what's up with Domesticus' exon 817526?

#### + dabg - diff ####
# now we can check how the difference relates to dabg
# outliers spoil the image again, filter them out

dsigdabg %>%
  filter(abs(diff) < 500) %>%
  filter(dabg < 0.05) %>%
  ggplot(aes(dabg.nosnp + 0.00001, diff)) +
  geom_point(alpha=0.02) +
  geom_smooth() +
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.05,0.1,0.5,1)) -> p_dabg_diff
  #scale_x_continuous(breaks=c(0.0001,0.001,0.01,0.05,0.1,0.5,1))
ggsave(file.path("results",output,"dabg-diff.png"),p_dabg_diff, width = 20, height = 10, units="cm")
ggsave(file.path("results",output,"dabg-diff.pdf"),p_dabg_diff, width = 20, height = 10, units="cm")

#### + dabg - signal ####
# this makes an interestign plot when not filtering for dabg (I guess it doesn't matter since it's on 1 side of dabg range)
dsigdabg %>%
  filter(abs(diff) < 500) %>%
  filter(dabg < 0.05) %>%
  ggplot(aes(dabg.nosnp, signal.nosnp)) +
  geom_point(alpha=0.02) +
  geom_smooth() + 
  scale_y_log10() -> p_dabg_signal
ggsave(file.path("results",output,"dabg-signal.png"), p_dabg_signal, width = 20, height = 10, units="cm")
ggsave(file.path("results",output,"dabg-signal.pdf"), p_dabg_signal, width = 20, height = 10, units="cm")

#### Treatment exploration ####
# now check if the treatment helped to reduce the variability
# of the affected probes
# need to use reshape2 cast, tidyr's spread is not versatile enough here

#### + dvar means and SDs ####
### could be computed from dsigdabg...
dm %>%
  filter(exon_cluster_id %in% common_hits) %>%
  ### WARNING: ONLY for Plier as gcRMA is already log-scaled
  mutate(signal=log10(signal)) %>%  # log2 would be more consistent with gcRMA, but absolute numbers (fold change) don't really matter in this case
  group_by(treat, exon_cluster_id) %>%
  summarise(mean=mean(signal), sd=sd(signal)) %>% 
  gather(measure, value, mean:sd) %>% 
  dcast(exon_cluster_id ~ treat + measure) %>% 
  mutate(sd_diff=original_sd - nosnp_sd, mean_diff=original_mean - nosnp_mean, sd_diff_random=original_sd - random_sd, mean_diff_random=original_mean - random_mean) ->
  dvar

dm %>%
  filter(exon_cluster_id %in% common_hits_random) %>%
  filter(treat %in% c("original","nosnp")) %>%
  select(-transcript_cluster_id,-SNP_dum,-SNP_probes,-total,-probesets) %>% 
  group_by(exon_cluster_id, chip) %>%
  summarise(diff=signal-lag(signal))
  mutate(signal=log2(signal)) %>%
  mutate(o_t_diff = )
  
  

#### + SD hist ####
# question: Is the standard deviation (as a measure of variance)
#           lower for the corrected data?
# distribution is centered around 0, but seems to be a bit lower
dvar %>%
  filter(abs(sd_diff) < .5) %>% 
  ggplot(aes(sd_diff)) + 
  geom_histogram() -> p_sd_diff_hist
ggsave(file.path("results",output,"sd_diff-hist.pdf"), p_sd_diff_hist, width = 20, height=10, units="cm")

dvar %>% {mean(.$sd_diff, na.rm=T)}

dvar %>%
  filter(abs(sd_diff_random) < .5) %>% 
  ggplot(aes(sd_diff_random)) + 
  geom_histogram() -> p_sd_diff_hist_random
ggsave(file.path("results",output,"sd_diff-hist_random.pdf"), p_sd_diff_hist_random, width = 20, height=10, units="cm")

dvar %>% {mean(.$sd_diff_random, na.rm=T)}

#### + mean hist ####
# once more the same for mean
dvar %>%
  ggplot(aes(mean_diff)) + 
  geom_histogram(binwidth=0.05) +
  xlim(-0.5, 1) -> p_mean_diff_hist
ggsave(file.path("results",output,"mean_diff-hist.pdf"), p_mean_diff_hist, width = 20, height = 10, units = "cm")

dvar %>%
  ggplot(aes(mean_diff_random)) + 
  geom_histogram(binwidth=0.05) +
  xlim(-0.5, 1) -> p_mean_diff_hist_random
ggsave(file.path("results",output,"mean_diff-hist_random.pdf"), p_mean_diff_hist_random, width = 20, height = 10, units = "cm")

# the expresion values got higher 
# (no surprise, if they got higher pairwise)
# this is fold change, if we're in log space (log(a/b) == log a - log b)
dvar %>% {mean(.$mean_diff, na.rm=T)}

#### + sd - mean ####
# is the change in sd dependent on mean?
# - looks like most of the smoother is below 0, 
#   meaning we get a bit more variable data
dvar %>%
  filter(original_mean > 0.1) %>%
  ggplot(aes(original_mean, sd_diff)) +
  geom_point(alpha=0.1) +
  geom_smooth(colour="red") +
  ylim(-0.4, 0.4) +
  scale_x_log10() -> p_mean_sd_diff_log
ggsave(file.path("results",output,"mean-sd_diff.pdf"), p_mean_sd_diff, width = 20, height = 10, units = "cm")

# TODO: are those ↑↓ plots identical???? 

# a special variant for log space
dvar %>%
  filter(original_mean > -2) %>%
  ggplot(aes(original_mean, sd_diff)) +
  geom_point(alpha=0.1) +
  geom_smooth(colour="red") -> p_mean_sd_diff
ggsave(file.path("results",output,"mean-sd_diff-log.pdf"), p_mean_sd_diff, width=20, height=10, units= "cm")

#### + diff to diff ####
### XXX: is this really log space??
# diff to diff plot;)
# no surprise again, quite nice linear fit 
# in log-log scale, the higher the mean difference, the 
# higher the variance difference
# possible explanation: one of the value got significantly higner
#  so it increased the mean and the sd at the same time
# beware for log scale for the diff - there is a lot of negative values
# taking absolute value should not harm the distribution that much..;)
# - a dirty trick of folding the axes, and distinguishing the quadrants based on shapes;)
dvar %>%
  filter(mean_diff + sd_diff > 0) %>%
  ggplot(aes(abs(mean_diff), abs(sd_diff), 
             shape=paste(mean_diff < 0, sd_diff < 0),
             size=original_mean)) +
  geom_point(alpha=0.2) +
  scale_x_log10() +
  scale_y_log10() -> p_sd_diff_mean_diff
ggsave(file.path("results",output,"mean_diff-sd_diff.png"), p_sd_diff_mean_diff, width = 30, height = 30, units = "cm")

# things get easier in the log space
# - most of the poins seem to be centered around 0,0
#   and some show an inverse dependency - if the mean value of exon
#   grew after treatment, the sd grew as well.. noting surprising
dvar %>%
  ggplot(aes(sd_diff, mean_diff,
             size=nosnp_mean)) +
  geom_point(alpha=0.05) -> p_sd_diff_mean_diff_log
ggsave(file.path("results",output,"sd_diff-mean_diff-log.png"), p_sd_diff_mean_diff_log, width = 30, height = 30, units = "cm")

dvar %>% filter(mean_diff == 0) %>% summarise(n())
dvar %>% filter(sd_diff == 0) %>% summarise(n())

##FIXME: the question is - should i calculate sd on raw or log-scaled values?
#  i guess that without log scaling - the higher the values, the higher the sd..
# - using sd(log10(signal)) removes the correlation in the diff-diff plot

# what if we normalize the sd with the mean - even better, 
# to stay in the same units, we'd rather use var, that after dividing
# with mean should give the same units as signal
# we're in log space, sd and mean have the same units - use difference instead of fraction
# - looks like some (few) results got more noisy, none got less noisy
dvar %>%
  ggplot(aes(original_sd - original_mean, nosnp_sd - nosnp_mean)) +
  geom_point(alpha=0.1) +
  geom_smooth(colour="red") +
  ggtitle("comparison of sds normalized by mean") -> p_sd_mean_original_treated
ggsave(file.path("results",output,"sd_mean-original-treated.png"), p_sd_mean_original_treated, width = 30, height = 30, units = "cm")

# this will make more sense in the light of different comparisons
# domesticus is actually supposed to change a little
# musculus results shold differ more dramatically

# the basic idea of the work is to leave out the probes where the animals 
# are genetically different from the reference used to design the chip 
# from the summarization to correct for the genetic differences between the subspecies
# and the reference genome
#
# to get some sense of validity of the approach, we need to check 
# the robustness of the summarization methods to leaving out probes
# robustness means giving the same results for treated and untreated data
# in the case where the animal's genome is the same as the reference one

# the best way to compare this is to make a model of probe hits based on the 
# real world data, and then randomly leave out probes so it fits the model 
# withou the biological factors (i.e. leave out probes from the exon clusters 
# based on some conditional proability distribution - would this be the same
# as leaving out probes completely at random?)
# and then compare this to the original results and the biology-based ones
# the expectation is that after leaving out probes with poora affinity
# the signal will go up
# in short, the signal level increase should ordered like this: 
# random treated D (~ 0) < bio-treated D < bio-treated M
# where for the same-species the comparison can be made on the chip level, 
# for the D-M the comparison has to be made on the exon cluster level (there is no pairs)

#### difference test -- variance exploration ####
dsigdabg %>%
  filter(dabg<0.05) %>%
  filter(exon_cluster_id %in% common_hits) %>%
  mutate(signal.original=log2(signal.original)) %>% 
  mutate(signal.nosnp=log2(signal.nosnp)) %>% 
  ggplot() +
  geom_histogram(aes(signal.nosnp), alpha=0.2, binwidth=0.3) +
  geom_histogram(aes(signal.original), alpha=0.2, fill="red", binwidth=0.3)

### https://en.wikipedia.org/wiki/F-test
### https://en.wikipedia.org/wiki/Paired_difference_test
### https://en.wikipedia.org/wiki/Student%27s_t-test
# within-group variability : ((nosnp_chip - exon_mean)^2 + (snp_chip - exon_mean)^2) / 
dsigdabg %>% 
  filter(dabg<0.05) %>% 
  filter(exon_cluster_id %in% common_hits) %>%
  mutate(signal.original=log2(signal.original)) %>% 
  mutate(signal.nosnp=log2(signal.nosnp)) %>% 
  ggplot() + 
  stat_qq(aes(sample=signal.original, colour="red")) + 
  stat_qq(aes(sample=signal.nosnp))


# is this considered a normal distribution? 
# https://explorable.com/dependent-t-test-for-paired-samples
# "It has however been shown that minor departures from normality do not affect this test - this is indeed an advantage."

#### + t stastistic ####
# http://www.statstutor.ac.uk/resources/uploaded/paired-t-test.pdf
# n ... number of chips
# d_i = y_i - x_i ==> diff = signal.treated - signal.original
# d_mean = mean(diff)
# s_d = sd(diff)
# SE(d_mean) = s_d/sqrt(n)  -- not really necessary as it can be expressed directly in the T computation
#
# T = d_mean/SE(d_mean) = (d_mean * sqrt(n))/s_d
# t = T on n-1 df
# p based on t distribution

dm %>% 
  filter(treat %in% c("nosnp", "original")) %>%             # analysis for SNP treated
  filter(exon_cluster_id %in% common_hits) %>%
  mutate(signal = log2(signal)) %>%                         # WARNING: log transformation
  group_by(transcript_cluster_id,exon_cluster_id,treat) %>% 
  mutate(n = length(chip)) %>%                              # count chips in group
  mutate(dabg_max = max(dabg)) %>%                          # extract max dabg between chips of the same exon / treatment
  ungroup %>%
  group_by(transcript_cluster_id, exon_cluster_id) %>% 
  filter(dabg_max=all(dabg_max<0.05)) %>%                   # drop exons that exceed dabg treshold in any treatment
  ungroup %>% 
  group_by(transcript_cluster_id,exon_cluster_id, chip) %>% # there are only two values when grouping by chip now
  arrange %>%
  mutate(diff = ifelse(is.na(lead(signal)),signal - lag(signal), signal - lead(signal))) %>%
  ungroup %>%
  group_by(transcript_cluster_id, exon_cluster_id) %>%
  filter(treat == "nosnp") %>%
  mutate(d_mean = mean(diff)) %>%
  mutate(s_d = sd(diff)) %>%
  mutate(T = d_mean * sqrt(n)/s_d) %>%
  mutate(t = pt(T,n-1)) -> t_stat

summary(t_stat)
  
t_stat %>%
  ggplot() +
  geom_histogram(aes(t), binwidth=0.05)

t_stat %>%
  ggplot() +
  geom_point(aes(t, signal), alpha=0.05) +
  geom_vline(xintercept = 0.05, colour="blue")
  


  
  
  
  
  
  
  
  















