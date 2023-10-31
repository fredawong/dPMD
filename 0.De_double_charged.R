library(tidyverse)

if(!dir.exists("mz_de_double_charged")){
  dir.create("mz_de_double_charged")
}

files = dir(path="F:/OneDrive/Study/paper/WWTP2021/FT/mz",pattern = "*.txt",full.names = T)
now.sample = 1
for(file in files){
  message(paste("Now procssing ",stringr::str_remove(file,".*\\/"),
                " (",now.sample,"/",length(files),")...",sep=""))
  data = read_tsv(file)
  masses = data$`m/z`
  intensities = data$I
  mass_defects=masses-floor(masses)
  c13d2=(13.00335-12)/2
  ppm_tol=1
  flags=rep(0,length(masses))
  
  for(m in 1:length(masses)){
    if (mass_defects[m]>0.4&mass_defects[m]<0.8) {
      flags[m]=1
      mono_mass=masses[m]-c13d2
      m_diff=abs(masses-mono_mass)
      mono_pos=which(m_diff==min(m_diff))
      if (min(m_diff)/m*1e6<ppm_tol) {
        if (any(intensities[mono_pos]/intensities[m]<10)){
          flags[mono_pos[which(intensities[mono_pos]/intensities[m]<10)]]=2
        }
      }
    }
  }
  message(paste("Retain ",length(which(flags==0))," masses in ",stringr::str_remove(file,".*\\/"),
                " (",round(length(which(flags==0))/length(masses),4)*100,"%)",sep=""))
  # Export
  write_tsv(data[which(flags==0),],paste("mz_de_double_charged/",
                                       stringr::str_remove(file,".*\\/"),sep=""))
  now.sample=now.sample+1
}

