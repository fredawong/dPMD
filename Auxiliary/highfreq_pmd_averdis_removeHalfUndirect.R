highfreq_pmd_averdis = function(pmd.3col,topN,prefix="Undefined"){
  pmd.freq = pmd.3col %>%
	filter(!(T1=="Persisted"&T2=="Persisted"&PMD<0)) %>%
    group_by(PMD.abs) %>%
    summarise(n=n()) %>%
    arrange(-n) %>%
    mutate(Index = 1:length(unique(pmd.3col$PMD.abs)))
  
  dis = c()
    for(i in c(1:ifelse(nrow(pmd.freq) > topN, topN, nrow(pmd.freq)))){
    pmd = as.numeric(pmd.freq$PMD.abs)[1:i]
    dfx = pmd.3col[pmd.3col$PMD.abs %in% pmd, c(1, 2)]
    net = igraph::graph_from_data_frame(dfx, directed = F)
    dis[i] = igraph::mean_distance(net)
    }
  n = which.max(dis)
  freqt = pmd.freq$n[n]
  sda = pmd.freq$PMD.abs[pmd.freq$n >= freqt]
  if (n == topN) {
    warning("Average distance is still increasing, increase topN.")
  }
  message(paste("PMD frequency cutoff is", freqt, 
                "by PMD network analysis with largest network average distance", 
                round(max(dis), 2), "."))
  message(paste(length(sda), "groups were found as high frequency PMD group."))
  message(paste(sda, "was found as high frequency PMD.", 
                "\n"))
  #return(data.frame(Name=prefix,PMD.abs=sda,Method="Average distance",n=pmd.freq$n[which(pmd.freq$PMD.abs==sda)]))
  return(data.frame(Name=prefix,PMD.abs=sda,n=pmd.freq$n[which(pmd.freq$PMD.abs==sda)]))
}

