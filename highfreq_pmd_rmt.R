highfreq_pmd_rmt = function(pmd.3col,topN,prefix="Undefined",
                            discard.outliers=T,
                            unfold.method = "gaussian",
                            max.ev.spacing = 3,
                            dist.method = "LL",
                            nr.breaks = 101){
  pmd.freq = pmd.3col %>%
    group_by(PMD.abs) %>%
    summarise(n=n()) %>%
    arrange(-n) %>%
    mutate(Index = 1:length(unique(pmd.3col$PMD.abs)))
  
  results = list()
  for(i in c(1:ifelse(nrow(pmd.freq) > topN, topN, nrow(pmd.freq)))){
    pmd = as.numeric(pmd.freq$PMD.abs)[1:i]
    dfx = pmd.3col[pmd.3col$PMD.abs %in% pmd, c(1, 2)]
    net = igraph::graph_from_data_frame(dfx, directed = F)
    
    # Remove outlier
    remove.outliers = function(x, factor = 1.5) {
      q25 = quantile(x, probs = 0.25)
      q75 = quantile(x, probs = 0.75)
      iqr = unname(q75 - q25)
      lower.threshold = q25 - (iqr * factor)
      upper.threshold = q75 + (iqr * factor)
      return(x[(x >= lower.threshold) & (x <= upper.threshold)])
    }
    
    # Calculate eigenvalue
    adj.m = igraph::get.adjacency(net, names = T, sparse = F)
    eigenvalues = eigen(adj.m, only.values = T)$values
    eigenvalues = eigenvalues[order(eigenvalues)]/max(abs(eigenvalues))
    if (discard.outliers) {
      orig.number.ev = length(eigenvalues)
      eigenvalues = remove.outliers(unique(eigenvalues))
      new.number.ev = length(eigenvalues)
      nr.outliers.removed = orig.number.ev - new.number.ev
    }
    else {
      nr.outliers.removed = 0
    }
    cat(paste("  Number of discarded outlier eigenvalues:", 
              nr.outliers.removed, "\n"))
    
    # Unfold
    if (unfold.method == "gaussian") {
      uf = RMThreshold::rm.unfold.gauss(eigenvalues)
    }
    if (unfold.method == "spline") {
      nr.fit.points = min(251, floor(0.75 * length(eigenvalues)))
      uf = RMThreshold::rm.unfold.spline(eigenvalues, 
                                         nr.fit.points = nr.fit.points)
    }
    ev.spacing = uf$ev.spacing
    ev = uf$eigenvalues
    l1 = length(ev.spacing)
    if (!is.null(max.ev.spacing)) {
      ev.spacing = ev.spacing[ev.spacing <= max.ev.spacing]
    }
    l2 = length(ev.spacing)
    cat(paste("      Number of large spacings not considered ( larger than", 
              max.ev.spacing, ") :", l1 - l2, "\n"))
    epsilon = max.ev.spacing/1000
    nr.small.spacings = sum(ev.spacing < epsilon)
    perc.small.spacings = nr.small.spacings/length(ev.spacing) * 
      100
    cat(paste("      Percentage of small spacings ( <", 
              epsilon, ") =", round(perc.small.spacings, 2), "\n"))
    ev.spacing = ev.spacing[ev.spacing > epsilon]
    p.val.ks.test = ks.test(unique(ev.spacing), "pexp", 
                            1)$p.value
    sse.exp = RMThreshold::rm.sse(ev.spacing)
    # Get results
    results[["nr.small.spacings"]][i] = nr.small.spacings
    results[["perc.small.spacings"]][i] = perc.small.spacings
    results[["nr.spacings"]][i] = length(ev.spacing)
    results[["tested.thresholds"]][i] = i
    results[["nr.uniq.ev"]][i] = length(unique(ev))
    results[["max.ev.mult"]][i] = max(table(ev))
    if (discard.outliers) {
      results[["nr.outliers.removed"]][i] = nr.outliers.removed
    }
    results[["p.ks"]][i] = p.val.ks.test
    results[["sse.exp"]][i] = sse.exp
    dres = RMThreshold::rm.get.distance(ev.spacing, dist.method = dist.method,
                                        nr.breaks = nr.breaks)
    results[["dist.Wigner"]][i] = dres$dist.Wigner
    results[["dist.Expon"]][i] = dres$dist.Expon
  }
  thresholds = results[["tested.thresholds"]]
  nr.unique.ev = results[["nr.uniq.ev"]]
  max.ev.mult = results[["max.ev.mult"]]
  nr.zeros = results[["nr.zeros"]]
  nr.spacings = results[["nr.spacings"]]
  dist.Wigner = results[["dist.Wigner"]]
  dist.Expon = results[["dist.Expon"]]
  sum.sq.err = results[["sse.exp"]]
  
  #fn = paste("RMT.num.uniq.ev", "png", sep = ".")
  #png(fn)
  #mtxt = "Number of unique eigenvalues vs. threshold"
  #plot(thresholds, nr.unique.ev, col = "blue", main = mtxt, 
  #     font.main = 1, xlab = "threshold", ylab = "nr unique ev")
  #dev.off()
  
  #fn = paste("RMT.max.ev.mult", "png", sep = ".")
  #png(fn)
  #mtxt = "Maximum eigenvalue multiplicity vs. threshold"
  #plot(thresholds, max.ev.mult, col = "blue", main = mtxt, 
  #     font.main = 1, xlab = "threshold", ylab = "max. ev multiplicity")
  #dev.off()
  
  #fn = paste("RMT.num.ev.spacings", "png", sep = ".")
  #png(fn)
  #mtxt = "Number of ev spacings vs. threshold"
  #plot(thresholds, nr.spacings, col = "blue", main = mtxt, 
  #     font.main = 1, xlab = "threshold", ylab = "nr. ev spacings")
  #dev.off()
  
  fn = paste("RMT.Dist.vs.Thres", "png", sep = ".")
  if(dist.method == "LL"){
    main.res = RMThreshold::rm.likelihood.plot(thresholds, log.le = dist.Expon, 
                                               log.lw = dist.Wigner, fn = fn, 
                                               interactive = F)
  }
  
  if (dist.method == "KLD") {
    main.res = RMThreshold::rm.distance.plot(thresholds, dist.Expon = dist.Expon, 
                                dist.Wigner = dist.Wigner,  
                                fn = fn, interactive = F)
  }
  
  cat(paste("\n  Distance plot saved to '", fn, "' \n\n"))
  results[["distance.plot"]] = fn
  results[["chosen.thresholds"]] = main.res$chosen.thresholds
  fn = "RMT.pks.vs.Thres.png"
  maintxt = "p-values for KS-test"
  tres <- rm.show.test(thresholds = thresholds, p.values = results[["p.ks"]], 
                       main = maintxt, fn = fn)
  cat(paste("\n  Plot with KS-test saved to '", fn, "' \n\n"))
  results[["p.ks.plot"]] = fn
  results[["p.ks.chosen"]] = tres$chosen
  fn = "RMT.SSE.vs.Thres.png"
  maintxt = "SSE for NNSD <--> Exponential"
  sse.res <- rm.sse.plot(thresholds = thresholds, sse.values = sum.sq.err, 
                         main = maintxt, fn = fn)
  cat(paste("\n  SSE plot saved to '", fn, "' \n\n"))
  results[["sse.plot"]] = fn
  results[["sse.chosen"]] = sse.res$chosen
  return(results)
  
  n = rmt.threshold
  freqt = pmd.freq$n[n]
  sda = pmd.freq$PMD.abs[pmd.freq$n >= freqt]
  if (n == topN) {
    warning("Transition form Possion to GOE not observed, increase topN.")
  }
  message(paste("PMD frequency cutoff is ", freqt,", by RMT.",sep=""))
  message(paste(length(sda), "groups were found as high frequency PMD group."))
  message(paste(sda, "was found as high frequency PMD.", 
                "\n"))
  #return(data.frame(Name=prefix,PMD.abs=sda,Method="RMT",n=pmd.freq$n[which(pmd.freq$PMD.abs==sda)]))
  return(data.frame(Name=prefix,PMD.abs=sda,n=pmd.freq$n[which(pmd.freq$PMD.abs==sda)]))
}