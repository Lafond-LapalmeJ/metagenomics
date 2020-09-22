plot_taxa_summary = function(physeq, Rank, GroupBy = NULL, TopN=-1, showLegend = FALSE){
  # Get taxa summary table  
  dt1 = sum_taxo(physeq, Rank = Rank, GroupBy = GroupBy)
  dt2 = sum_taxo(physeq, Rank = Rank)
  
  if (TopN>-1) {
    if (TopN>nrow(dt2)) TopN=nrow(dt2)
    indices=unlist(dt1[,1]) %in% unlist(dt2[1:TopN,1])
    dt1=dt1[indices,]
  }
  
  # Set factor appropriately for plotting
  RankCol = which(colnames(dt1) == Rank)
  setorder(dt1, -meanRA)
  dt1[, RankFac := factor(dt1[[Rank]],levels = rev(unique(dt1[[Rank]])))]
  dt1[, ebarMax := max(c(0, min(meanRA + sdRA))), by = eval(Rank)]
  dt1[, ebarMin := max(c(0, min(meanRA - sdRA))), by = eval(Rank)]
  # Set zeroes to one-tenth the smallest value
  ebarMinFloor = dt1[(ebarMin > 0), min(ebarMin)]
  ebarMinFloor <- ebarMinFloor / 10
  dt1[(ebarMin == 0), ebarMin := ebarMinFloor]
  pRank = ggplot(dt1, aes(x = meanRA, y = RankFac)) +
    scale_x_log10() +
    xlab("Mean Relative Abundance (log10)") +
    ylab(Rank) +
    theme_bw()
  if(!is.null(GroupBy)){
    # pRank <- pRank + facet_wrap(facets = as.formula(paste("~", GroupBy)))
    pRank <- pRank + geom_point(mapping = aes_string(colour = GroupBy),
                                size = 5, show.legend = showLegend)
  } else {
    # Don't include error bars for faceted version
    pRank <- pRank + geom_errorbarh(aes(xmax = ebarMax,
                                        xmin = ebarMin), show.legend = showLegend)
  }
  return(pRank)
}  
