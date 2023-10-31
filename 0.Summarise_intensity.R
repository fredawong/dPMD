library(tidyverse)

files = dir(path="F:/OneDrive/Study/paper/WWTP2021/directedPMD2/mz_de_double_charged/",pattern = "*.txt",full.names = T)

for (file in files) {
  dat  = read_tsv(file)
  tmp = data.frame(Name=str_remove(str_remove(file,".*\\/"),"\\.txt"),
                   Isum = sum(dat$I))
  if (file == files[1]) {
    res = tmp
  }else{
    res = rbind(res,tmp)
  }
}
write_tsv(res,"Intensity_sum.tsv")
