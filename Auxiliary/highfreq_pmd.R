highfreq_pmd = function(pmd.3col,topN,prefix="Undefined",filter.method="Aver.distance"){
  if(!filter.method%in%c("Aver.distance","Power.law","RMT")){
    stop("\n\n  highfreq_pmd: 'method' must be one of 'Aver.distance', 'Power.law' or 'RMT'.\n\n")
  }
  if(filter.method == "Aver.distance"){
    res = highfreq_pmd_averdis(pmd.3col = pmd.3col,
                               topN = topN,
                               prefix = WWTP)
  }
  ## 2.High-freq pmd by power law (scale-free)
  if(filter.method == "Power.law"){
    res = highfreq_pmd_powerlaw(pmd.3col = pmd.3col,
                                topN = topN,
                                prefix = WWTP)
  }
  ## 3.High-freq pmd by RMT
  if(filter.method == "RMT"){
    res = highfreq_pmd_rmt(pmd.3col = pmd.3col,
                           topN = topN,
                           prefix = WWTP,
                           discard.outliers = T,
                           unfold.method = "gaussian")
  }
  return(res)
}
