library(tidyverse)
#library(xcms)

#files = dir(path="../FT/mz/",pattern = "*.txt",full.names = T) # Raw data
#files = dir(path="mz_deiso/",pattern = "*.txt",full.names = T) # De-Iso data
files = dir(path="mz_deiso_recal/",pattern = "*.txt",full.names = T) 
for(file in files){
  if(file==files[1]){
    i=1
    j=1
    mz.dat = read_tsv(file) %>%
      mutate(Sample=str_remove(str_remove(file,"\\.txt"),".*\\/"),
             Index=i,
             Location = str_remove(Sample,".*_"),
             Location.index = case_when(Location=="IN"~1,
                                        Location=="SED"~2,
                                        Location=="OUT"~3))
  }else{
    i=i+1
    tmp = read_tsv(file) %>%
      mutate(Sample=str_remove(str_remove(file,"\\.txt"),".*\\/"),
             Index=i,
             Location = str_remove(Sample,".*_"),
             Location.index = case_when(Location=="IN"~1,
                                        Location=="SED"~2,
                                        Location=="OUT"~3))
    mz.dat = rbind(mz.dat,tmp)
  }
}
mz.dat = mz.dat %>%
  mutate(Peak.index = row.names(mz.dat))

# Cluster mz
ppm.threshold=10

mz.peak = data.frame(mz=mz.dat$`m/z`,rt=-1,sample=mz.dat$Index)
sampleGroup0 = data.frame(sample=mz.dat$Index,sample.group=mz.dat$Location) %>%
  unique()
sampleGroup = sampleGroup0$sample.group

mzcluster = xcms::do_groupPeaks_mzClust(peaks = mz.peak,
                                        sampleGroups = sampleGroup,
                                        ppm=ppm.threshold,
                                        absMz = 0,
                                        minFraction = 2/12,
                                        minSamples = 2)
head(mzcluster[[1]]) #featureDefinitions
head(mzcluster[[2]]) #peakIndex

# Post processsing
aligned.peak = data.frame(mzcluster[[1]])
aligned.peak = aligned.peak %>%
  mutate(Aligned.peak.index = row.names(aligned.peak))

aligned.peak.group = data.frame(Aligned.peak.index=aligned.peak$Aligned.peak.index,Peak.list=0)
for(i in 1:nrow(aligned.peak)){
  j=1
  for(j in 1:length(unlist(mzcluster[[2]][i]))){
    if(j==1){
      peak.list=as.character(unlist(mzcluster[[2]][i])[j])
    }else{
      peak.list = paste(peak.list,unlist(mzcluster[[2]][i])[j],sep=",")
    }
  }
  aligned.peak.group$Peak.list[i] = peak.list
}

col.num = max(str_count(aligned.peak.group$Peak.list,","))+1 # determine the col number needed
aligned.peak.group.1to1 = aligned.peak.group %>%
  separate(Peak.list, paste0("X", 1:col.num, seq = ""), sep = ",") %>% 
  pivot_longer(-Aligned.peak.index,values_to = "Peak.index") %>%
  dplyr::select(Aligned.peak.index, Peak.index) %>% 
  na.omit() %>%
  base::unique()

# Merge data
mz.dat.aligned = left_join(mz.dat,aligned.peak.group.1to1)
## Found dup
dup = which(duplicated(mz.dat.aligned[,c("Index","Aligned.peak.index")])&!is.na(mz.dat.aligned$Aligned.peak.index))
mz.dat.aligned = mz.dat.aligned %>%
  mutate(Presence = case_when(is.na(Aligned.peak.index)~"Singleton",
                              Peak.index %in% dup~"Duplicate",
                              T~"Peak"))

## Statics
statics = mz.dat.aligned %>%
  group_by(Sample,Presence) %>%
  summarise(Count=n(),I=sum(I))

library(paletteer)
ggplot(statics, aes(x=Sample,y=Count, fill = Presence)) +
  geom_col(position = "fill",width = .8) +
  scale_fill_paletteer_d("rcartocolor::Bold") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  ylab("Proportion") +
  theme(axis.text.x = element_text(size = 8,angle = 60,
                                   hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        axis.ticks = element_line(size = .5),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.height = unit(6,"mm"),
        legend.key.width = unit(6,"mm"),
        legend.text = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.background = element_blank(),
        panel.border = element_rect(size = .5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        aspect.ratio = 1)
ggsave(paste("1.Sample.mz.statics_ppm",ppm.threshold,".pdf",sep=""), 
       wi = 6, he = 4.5, useDingbats = F)

## Discard singleton and sum up duplicate
mz.dat.aligned2 = mz.dat.aligned %>%
  filter(Presence != "Singleton") %>%
  left_join(aligned.peak[,-c(4:6)]) %>%
  group_by(Sample,Aligned.peak.index) %>%
  summarise(mz=mzmed,I=sum(I)) %>%
  base::unique()
### Wider
mz.dat.aligned3 = mz.dat.aligned2 %>%
  pivot_wider(names_from = "Sample", values_from = "I") %>%
  arrange(as.numeric(Aligned.peak.index))
write_tsv(mz.dat.aligned3,
          paste("1.Aligned_mz_ppm",ppm.threshold,".tsv",sep=""),
          na = "")
