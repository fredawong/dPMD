library(tidyverse)
source("Auxiliary/highfreq_pmd.R")
source("Auxiliary/highfreq_pmd_averdis.R")

# Parameters
## Global
digit = 6
accuracy = 3
ppm.threshold = 10
suffix = "corrected"
filter.method = c("Aver.distance","Power.law","RMT")[1]
## RMT
discard.outliers=T
unfold.method = "gaussian"
max.ev.spacing = 3
dist.method = "LL"
nr.breaks = 101

# Load global highfreq pmd network
pmd.highfreq = read_tsv(paste("2.Highfreq.PMD.global.",filter.method,"_ppm",
                              ppm.threshold,suffix,".tsv",sep=""))

## Define remove outlier function
remove.outliers = function(x, factor = 1.5) {
  q25 = quantile(x, probs = 0.25)
  q75 = quantile(x, probs = 0.75)
  iqr = unname(q75 - q25)
  lower.threshold = q25 - (iqr * factor)
  upper.threshold = q75 + (iqr * factor)
  return(x[(x >= lower.threshold) & (x <= upper.threshold)])
}

# Build network and test by RMT
dist.Wigner = 0
dist.Expon = 0
dist.Wigner.list = c()
dist.Expon.list = c()
threshold.now = max(pmd.highfreq$Presence)
i=1
while (!(dist.method=="LL"&dist.Expon<dist.Wigner|
       dist.method=="KLD"&dist.Expon>dist.Wigner|
       threshold.now==1)) {
  
  message(paste("Run ",i,": PMD presence threshold:",threshold.now,sep=""))
  dfx = pmd.highfreq[pmd.highfreq$Presence>=threshold.now,c("mz1","mz2")]
  net = igraph::graph_from_data_frame(dfx, directed = F)
  
  ## Calculate eigenvalue
  adj.m = igraph::get.adjacency(net, names = T, sparse = T)
  eigenvalues = eigen(adj.m, only.values = T)$values
  eigenvalues = eigenvalues[order(eigenvalues)]/max(abs(eigenvalues))
  if (discard.outliers) {
    orig.number.ev = length(eigenvalues)
    eigenvalues = remove.outliers(unique(eigenvalues))
    new.number.ev = length(eigenvalues)
    nr.outliers.removed = orig.number.ev - new.number.ev
  }else {
    nr.outliers.removed = 0
  }
  cat(paste("Number of discarded outlier eigenvalues:", 
            nr.outliers.removed, "\n"))
  ## Unfold
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
  cat(paste("Number of large spacings not considered ( larger than", 
            max.ev.spacing, ") :", l1 - l2, "\n"))
  epsilon = max.ev.spacing/1000
  nr.small.spacings = sum(ev.spacing < epsilon)
  perc.small.spacings = nr.small.spacings/length(ev.spacing) * 100
  cat(paste("Percentage of small spacings ( <", 
            epsilon, ") =", round(perc.small.spacings, 2), "\n"))
  ev.spacing = ev.spacing[ev.spacing > epsilon]
  p.val.ks.test = ks.test(unique(ev.spacing), "pexp", 
                          1)$p.value
  sse.exp = RMThreshold::rm.sse(ev.spacing)
  ## Get likelihood or distance
  dres = RMThreshold::rm.get.distance(ev.spacing, dist.method = dist.method,
                                      nr.breaks = nr.breaks)
  dist.Wigner = dres$dist.Wigner
  dist.Expon = dres$dist.Expon
  dist.Wigner.list[i] = dist.Wigner
  dist.Expon.list[i] = dist.Expon
  names(dist.Wigner.list)[i] = threshold.now
  names(dist.Expon.list)[i] = threshold.now
  message(paste(dist.method,": Wigner=",dist.Wigner," dist.Expon=",dist.Expon,sep=""))
  
  threshold.now = threshold.now - 1
  i=i+1
}

if(threshold.now == 1){
  message("No threshold could be achieved!")
}else{
  threshold.sel = threshold.now + 1
  message(paste("Get threshold:",threshold.sel))
}

# Get highfreq PMDs above threshold
pmd.highfreq.sel = pmd.highfreq[pmd.highfreq$Presence>=threshold.sel,]
write_tsv(pmd.highfreq.sel,
          paste("3.Pruned_highfreq_PMD_threshold_",threshold.sel,".tsv",sep=""))
