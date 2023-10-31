library(MFAssignR)
source("MFAssignCHO_RMD_modified.R")

files = dir(path="F:/OneDrive/Study/paper/WWTP2021/directedPMD2/mz_de_double_charged",pattern = "*.txt",full.names = T)
set_SN = 0

if(!dir.exists("mz_deiso_recal")){
  dir.create("mz_deiso_recal")
}

now.sample = 1
for(file in files){
  message(paste("Now procssing ",stringr::str_remove(file,".*\\/"),
                " (",now.sample,"/",length(files),")...",sep=""))
  Data = read.table(paste(file), sep = "\t", header = T)
  # Isotope Filtering
  Mono_Iso = MFAssignR::IsoFiltR(Data, SN = 0)
  Mono = Mono_Iso[["Mono"]]
  Iso = Mono_Iso[["Iso"]]
  
  # Calibrate
  Assign = MFAssignCHO_RMD2(peaks = Mono, isopeaks = Iso,
                                      ionMode = "neg", lowMW =100,
                                      highMW = 800,  Mx = 1, Ex = 1, Ox = 27,
                                      ppm_err = 3, H_Cmin = 0.3, Omin = 1,
                                      SN = set_SN, HetCut = "off", NMScut = "on")
  Unambig = Assign[["Unambig"]]
  Ambig = Assign[["Ambig"]]
  Unassigned = Assign[["None"]]
  RecalList = MFAssignR::RecalList(df = Unambig)
  callist = dplyr::arrange(RecalList, -RecalList[,2], -RecalList[,6])[,c(1,2,6)]
  callist = dplyr::filter(callist,`Abundance Score`>100&`Number Observed`>5)
  if(nrow(callist)<3){
    warning(paste("Insufficient calibrants: use original Mono and Iso for",
                  stringr::str_remove(file,".*\\/")))
    Recal_Mono = Mono
    Recal_Iso = Iso
  }else{
    Recalcheck = try(MFAssignR::Recal(df = Unambig, peaks = Mono,
                                  isopeaks = Iso, mode = "neg", 
                                  SN = set_SN, mzRange = 70, 
                                  series1 = callist$Series[1], series2 = callist$Series[2], 
                                  series3 = callist$Series[3], series4 = callist$Series[4],
                                  series5 = callist$Series[5], series6 = callist$Series[6],
                                  series7 = callist$Series[7], series8 = callist$Series[8],
                                  series9 = callist$Series[9], series10 = callist$Series[10]
                                  ))
    if(any(stringr::str_detect(Recalcheck[[1]],"Error"))){
      message(paste("Unable to calibrate: use original Mono and Iso for",
                    stringr::str_remove(file,".*\\/")))
      Recal_Mono = Mono
      Recal_Iso = Iso
    }else{
      Recal_Mono = Recalcheck[["Mono"]]
      Recal_Iso = Recalcheck[["Iso"]]
      List = Recalcheck[["RecalList"]]
    }
  }
  
  # Export
  names(Recal_Mono) = c("m/z","I")
  readr::write_tsv(Recal_Mono,paste("mz_deiso_recal/",
                                    stringr::str_remove(file,".*\\/"),sep=""))
  
  now.sample=now.sample+1
}
