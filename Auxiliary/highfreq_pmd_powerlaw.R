highfreq_pmd_powerlaw = function(pmd.3col,topN,prefix="Undefined"){
  pmd.freq = pmd.3col %>%
    group_by(PMD.abs) %>%
    summarise(n=n()) %>%
    arrange(-n) %>%
    mutate(Index = 1:length(unique(pmd.3col$PMD.abs)))
  
  powerLawP = c()
  for(i in c(1:ifelse(nrow(pmd.freq) > topN, topN, nrow(pmd.freq)))){
    pmd = as.numeric(pmd.freq$PMD.abs)[1:i]
    dfx = pmd.3col[pmd.3col$PMD.abs %in% pmd, c(1, 2)]
    net = igraph::graph_from_data_frame(dfx, directed = F)
    powerLawP[i] = igraph::fit_power_law(igraph::degree(net))$KS.p
  }
  n = min(which(powerLawP<0.05))
  if (is.infinite(n)&filter.method == "Power.law") {
    warning("Unable to fit power law distribution.")
    message("Using topN instead,but interpret results with caution!")
    n = topN
  }
  freqt = pmd.freq$n[n]
  sda = pmd.freq$PMD.abs[pmd.freq$n >= freqt]
  message(paste("PMD frequency cutoff is ", freqt,", by power law distribution.",sep=""))
  message(paste(length(sda), "groups were found as high frequency PMD group."))
  message(paste(sda, "was found as high frequency PMD.", 
                "\n"))
  #return(data.frame(Name=prefix,PMD.abs=sda,Method="Power law",n=pmd.freq$n[which(pmd.freq$PMD.abs==sda)]))
  return(data.frame(Name=prefix,PMD.abs=sda,n=pmd.freq$n[which(pmd.freq$PMD.abs==sda)]))
}
