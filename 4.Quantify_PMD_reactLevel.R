library(tidyverse)

highfreq.pmd.global = read_tsv("3.Pruned_highfreq_PMD_threshold_8.tsv")
highfreq.pmd.pairwise = read_tsv("2.Highfreq.PMD.pairwise.Aver.distance_ppm10corrected.tsv") %>%
  select(c(mz1,mz2,I1,I2,T1,T2,Prefix))
intensity.sum = read_tsv("Intensity_sum.tsv")
intensity.sum1 = intensity.sum %>%
  filter(str_detect(Name,"IN")) %>%
  mutate(Name=str_remove(Name,"_.*")) %>%
  rename(Prefix=Name, Isum1=Isum)
intensity.sum2 = intensity.sum %>%
  filter(str_detect(Name,"SED")) %>%
  mutate(Name=str_remove(Name,"_.*")) %>%
  rename(Prefix=Name, Isum2=Isum)

# Determine the degree of each node
WWTPs = unique(highfreq.pmd.pairwise$Prefix)
for (WWTP in WWTPs) {
  dfx = highfreq.pmd.pairwise %>%
    filter(Prefix == WWTP) %>%
    left_join(highfreq.pmd.global) %>%
    na.omit()
  net = igraph::graph_from_data_frame(dfx[,c(1,2)], directed = F)
  topology = data.frame(Prefix = WWTP,
                        Degree = igraph::degree(net, mode = "all",
                                                loops = F, normalized = F))
  topology$mz = as.double(row.names(topology))
  if (WWTP == WWTPs[1]) {
    topology.res = topology
  }else{
    topology.res = rbind(topology.res,topology)
  }
}
topology.res1 = topology.res %>%
  rename(mz1 = mz, Degree1 = Degree)
topology.res2 = topology.res %>%
  rename(mz2 = mz, Degree2 = Degree)

# Join absolute abundance, relative abundance and degree
dat = left_join(highfreq.pmd.global,highfreq.pmd.pairwise) %>%
  left_join(intensity.sum1) %>%
  left_join(intensity.sum2) %>%
  mutate(Ir1 = I1/Isum1,
         Ir2 = I2/Isum2) %>%
  select(-c(Isum1,Isum2)) %>%
  left_join(topology.res1) %>%
  left_join(topology.res2) %>%
  mutate(I1.weighted = round(I1/Degree1,0),
         I2.weighted = round(I2/Degree2,0),
         Ir1.weighted = Ir1/Degree1,
         Ir2.weighted = Ir2/Degree2) %>%
  mutate(ReactLevel.pairwise = (I1.weighted+I2.weighted)/2,
         ReactLevel.pairwise.r = (Ir1.weighted+Ir2.weighted)/2,
         difference.pairwise = abs(I1-I2)/max(I1,I2),
         difference.pairwise.r = abs(Ir1-Ir2)/max(Ir1,Ir2))
write_tsv(dat,"4.Pruned_highfreq_PMD_with_quant_global.tsv")

# Get sample-wise PMD traits quantification (using weighted direct)
dat.stat = dat %>%
  group_by(Prefix,PMD.abs,PMD.direct,PMD.direct.weighted) %>%
  summarise(Freq = n(),
            Ir1.sum = sum(Ir1),
            Ir2.sum = sum(Ir2),
            React.level = sum(ReactLevel.pairwise),
            React.level.r = sum(ReactLevel.pairwise.r),
            Difference.sum = sum(difference.pairwise),
            Difference = Difference.sum/Freq,
            Difference.r = sum(difference.pairwise.r)/Freq)
write_tsv(dat.stat,"4.Pruned_highfreq_PMD_with_quant_summary_global.tsv")

# Get sample-wise PMD traits quantification (using direct only)
dat.stat = dat %>%
  group_by(Prefix,PMD.abs,PMD.direct) %>%
  summarise(Freq = n(),
            Ir1.sum = sum(Ir1),
            Ir2.sum = sum(Ir2),
            React.level = sum(ReactLevel.pairwise),
            React.level.r = sum(ReactLevel.pairwise.r),
            Difference.sum = sum(difference.pairwise),
            Difference = Difference.sum/Freq,
            Difference.r = sum(difference.pairwise.r)/Freq)
write_tsv(dat.stat,"4.Pruned_highfreq_PMD_with_quant_summary_global_onlyDirect.tsv")
