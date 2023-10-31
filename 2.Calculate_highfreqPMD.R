setwd("F:/OneDrive/Study/paper/WWTP2021/directedPMD2")
# Load codes in R/
source("R/highfreq_pmd.R")
source("R/highfreq_pmd_averdis.R")

library(tidyverse)
digit = 6
accuracy = 3
is.knownPMD = F
Literature = T
filter.method = c("Aver.distance","Power.law","RMT")[1]
topN = 50
#dedupPMD = T
#correct.highfreqPMD = T
allow.error = 0.001 # 0.001 Da
ppm.threshold=10 # detemine the cluster data at specific ppm threshold

# Load aligned mz data
mz.dat = read_tsv(paste("1.Aligned_mz_ppm",ppm.threshold,".tsv",sep = "")) %>%
  mutate(mz = round(mz,digit))

if(Literature){
  # Load KEGG db
  literature.pmd = read_tsv("pmdDB_literature.tsv")
  kegg.pmd = read_tsv("pmdDB_KEGG.tsv")
  known.pmd = unique(c(round(kegg.pmd$pmd,accuracy),round(literature.pmd$PMD,accuracy)))
  rm(literature.pmd,kegg.pmd)
}else{
  kegg.pmd = read_tsv("pmdDB_KEGG.tsv")
  known.pmd = unique(c(round(kegg.pmd$pmd,accuracy)))
  rm(kegg.pmd)
}
known.pmd.upper = round(known.pmd + allow.error, accuracy)
known.pmd.lower = round(known.pmd - allow.error, accuracy)
known.pmd.all = data.frame(upper = known.pmd.upper,
                           lower = known.pmd.lower,
                           real.pmd = known.pmd)

WWTPs = unique(str_remove(names(mz.dat)[-c(1:2)],"_.*"))
now.sample = 1
pmd.result = list()
pmd.result.valid = list()

# Main
for(WWTP in WWTPs){
  message(paste("Now procssing ",WWTP," (",now.sample,"/",length(WWTPs),")...",sep=""))
  mz.sample1 = str_c(WWTP,"IN",sep="_")
  mz.sample2 = str_c(WWTP,"SED",sep="_")
  mz.list1 = mz.dat[,c("mz",mz.sample1)] %>%
    na.omit()
  mz.list2 = mz.dat[,c("mz",mz.sample2)] %>%
    na.omit()
  
  # Determine transformed, produced, and persist
  message("Determining the transformation status of mz")
  mz.list = full_join(mz.list1,mz.list2) %>%
    rename(I1=!!as.name(mz.sample1),I2=!!as.name(mz.sample2)) %>%
    mutate(Type = case_when(!is.na(I1)&is.na(I2)~"Completely transformed",
                            !is.na(I1)&!is.na(I2)&I1>5*I2~"Transformed",
                            is.na(I1)&!is.na(I2)~"Newly produced",
                            !is.na(I1)&!is.na(I2)&5*I1<I2~"Produced",
                            T~"Persisted")) %>%
    mutate(Type2 = case_when(Type=="Completely transformed"|Type=="Transformed"~"Transformed",
                             Type=="Newly produced"|Type=="Produced"~"Produced",
                             T~"Persisted"))
  transf.statics = mz.list %>%
    group_by(Type) %>%
    summarise(Count=n())
  #print(transf.statics)
  #write_tsv(transf.statics,paste("2.",str_remove(mz.sample1,"_.*"),".tsv",sep=""))
  
  t1 = mz.list %>%
    select(c(mz,Type2)) %>%
    rename(mz1=mz,T1=Type2) %>%
    filter(T1!="Produced") # Non-sense for mz1 as "Produced"
  t2 = mz.list %>%
    select(c(mz,Type2)) %>%
    rename(mz2=mz,T2=Type2)  %>%
    filter(T2!="Transformed") # Non-sense for mz2 as "Transformed"
  
  # Calculate PMDs
  message("Calculating all possible PMDs...")
  
  mz.list1 = mz.list1[mz.list1$mz %in% t1$mz1,]
  mz.list2 = mz.list2[mz.list2$mz %in% t2$mz2,]
  
  pmd.matrix = data.frame(outer(mz.list2$mz,mz.list1$mz,"-")) %>%
    round(accuracy)
  names(pmd.matrix) = mz.list1$mz
  pmd.matrix$X1 = mz.list2$mz
  
  pmd.3col = pivot_longer(pmd.matrix,-X1,names_to = "X2",values_to = "PMD") %>%
    filter(PMD!=0) %>%
    rename(mz1=X2,mz2=X1) %>%
    mutate(PMD.abs = abs(PMD),
           PMD.direct = case_when(PMD<0~"-",
                                  PMD>0~"+")) %>%
    mutate(mz1=as.numeric(mz1)) %>%
	left_join(t1) %>%
	left_join(t2) %>%
    select(c(mz1,mz2,PMD,PMD.abs,T1,T2,PMD.direct))
  
  # Find high-frequency PMDs based on mean distance of the network
  message(paste("Filtering high-frequency PMDs by",filter.method,"..."))
  
  highFreqPMD = highfreq_pmd(pmd.3col = pmd.3col,
                             topN = topN,
                             prefix = WWTP,
                             filter.method = filter.method)
   highFreqPMD.abs = unique(highFreqPMD$PMD.abs)
    
   highFreqPMD = pmd.3col %>%
      filter(PMD.abs %in% highFreqPMD.abs) %>%
      group_by(PMD.abs) %>%
      summarise(n=n()) %>%
      arrange(-n) %>%
      mutate(Name=WWTP) %>%
      select(Name,PMD.abs,n)
    message(paste(length(highFreqPMD$PMD.abs), "groups were finally found as high frequency PMD group."))
    message(paste(highFreqPMD$PMD.abs, "was finally found as high frequency PMD.", 
                  "\n"))
  }
  
  if(WWTP==WWTPs[1]){
    highFreqPMD.statics = highFreqPMD
  }else{
    highFreqPMD.statics = rbind(highFreqPMD.statics,highFreqPMD)
  }
  
  # Merge info
  message("Assessing direction strength...")
  
  pmd = pmd.3col %>%
    filter(PMD.abs %in% highFreqPMD$PMD.abs) %>%
    left_join(rename(mz.list1,mz1=mz,I1=!!as.name(mz.sample1))) %>%
    left_join(rename(mz.list2,mz2=mz,I2=!!as.name(mz.sample2))) %>%
    #left_join(t1) %>%
    #left_join(t2) %>%
    filter(T1!="Produced"&T2!="Transformed") %>%
    mutate(PMD.direct.strength = case_when(T1=="Transformed"&T2=="Produced"~"Strong",
                                           T1=="Persisted"&T2=="Produced"~"Weak",
                                           T1=="Transformed"&T2=="Persisted"~"Weak",
                                           T~"Undirected")) %>%
    mutate(PMD.direct.weighted = case_when(PMD.direct.strength=="Strong"~str_c(PMD.direct,PMD.direct,PMD.direct,sep = ""),
                                           PMD.direct.strength=="Weak"~str_c(PMD.direct,PMD.direct,sep = ""),
                                           T~PMD.direct))  %>%
    mutate(PMD = case_when(PMD.direct.weighted=="---"|
	                           PMD.direct.weighted=="--"|
	                           PMD.direct.weighted=="-"~PMD.abs*-1,
	                         PMD.direct.weighted=="+++"|
	                           PMD.direct.weighted=="++"|
	                           PMD.direct.weighted=="+"~PMD.abs))
  
  rm(highFreqPMD)
  pmd.result[[now.sample]] = pmd
  #pmd %>% group_by(PMD.direct,T1,T2) %>% summarise(n=n())
  #pmd %>% group_by(PMD.direct.weighted) %>% summarise(n=n())
  
  # Valid by known PMDs if set is.knownPMD = T
  if(is.knownPMD){
    message("Filtering valid reactions by KEGG database and literature...")
    pmd.valid = pmd %>%
      filter(PMD.abs%in%known.pmd)
    pmd.result.valid[[now.sample]] = pmd.valid
  }
  
  now.sample=now.sample+1
}

# Aggregate results
for(i in 1:length(pmd.result)){
  if(i==1){
    agg=mutate(pmd.result[[i]], Prefix=WWTPs[i])
  }else{
    agg=rbind(agg,mutate(pmd.result[[i]], Prefix=WWTPs[i]))
  }
}

suffix = ""
agg = agg %>% select(mz1,mz2,PMD,PMD.abs,I1,I2,T1,T2,PMD.direct,PMD.direct.strength,PMD.direct.weighted,Prefix)
write_tsv(agg,paste("2.Highfreq.PMD.pairwise.",filter.method,"_ppm",ppm.threshold,suffix,".tsv",sep=""))
#write_tsv(highFreqPMD.statics,paste("2.Highfreq.PMD.statics.",filter.method,"_ppm",ppm.threshold,suffix,".tsv",sep=""))

# Get global pmd network
#highFreqPMD.statics = read_tsv(paste("2.Highfreq.PMD.statics.",filter.method,
#                                     "_ppm",ppm.threshold,"corrected.tsv",sep=""))
agg = read_tsv(paste("2.Highfreq.PMD.pairwise.",filter.method,
                     "_ppm",ppm.threshold,"corrected.tsv",sep="")) %>%
  group_by(mz1,mz2) %>%
  summarise(PMD = names(which.max(table(PMD))),
            PMD.abs = names(which.max(table(PMD.abs))),
            PMD.direct = unique(PMD.direct),
            PMD.direct.weighted = names(which.max(table(PMD.direct.weighted))),
            Presence=n()) %>%
  arrange(-Presence)
write_tsv(agg,paste("2.Highfreq.PMD.global.",filter.method,"_ppm",ppm.threshold,suffix,".tsv",sep=""))
