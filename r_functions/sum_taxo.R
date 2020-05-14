sum_taxo <- function(physeq, Rank, GroupBy = NULL){
  # check if rank exist
  if(!Rank %in% rank_names(physeq)){
    message("The argument to `Rank` was:\n", Rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(physeq), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  # check if group exist 
  if(!is.null(GroupBy)){
    if(!GroupBy %in% sample_variables(physeq)){
      message("The argument to `GroupBy` was:\n", GroupBy,
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  
  # Merge any NA value
  taxo <-  data.frame(tax_table(physeq))
  indexNA1 <- which(is.na(taxo[[Rank]]))
  indexNA2 <- which(grepl(pattern = "__NA$",x = taxo[[Rank]]))
  indexNA3 <- which(taxo[[Rank]] == "NA")
  indexNA <- unique(c(indexNA1,indexNA2,indexNA3))
  if(length(unique(taxo[[Rank]][indexNA]))>1){
    warning("merging all NA values to 'NA'")
    print(unique(taxo[[Rank]][indexNA]))
    physeq <- merge_taxa(x = physeq, eqtaxa = indexNA,archetype = 1)
  }
  # Remove upper taxonomic level
  taxo <- data.frame(tax_table(physeq))
  taxo_lvl <- colnames(taxo)
  ind_rank = which(taxo_lvl == Rank)
  for(lvl in 1:(ind_rank-1)){
    t <- taxo_lvl[lvl]
    taxo[t] <- rep(t, length(taxo[[t]]))
  }
  tax_table(physeq) <- as.matrix(taxo)
  
  # group by selected taxonomy rank
  phy_rank <- tax_glom(physeq,Rank,NArm = FALSE)
  
  # get taxo and otu table
  taxo <- as.data.frame(phy_rank@tax_table@.Data)
  otu <-  data.frame(otu_table(phy_rank))
  # compute rel. abundance 
  rel_otu <- data.frame(apply(otu,MARGIN = 2, FUN = function(x){return(x/sum(x))}))
  
  # add taxonomy to otu table
  df <- rel_otu %>% mutate(temp = taxo[[Rank]])
  names(df)[names(df) == 'temp'] <- Rank
  
  # case no group
  if(is.null(GroupBy)){
    sumdf <- apply(df, 1, FUN = function(x){
      name = x[Rank]
      data = as.numeric(x[-length(x)])
      meanRA = mean(data)
      sdRA = sd(data)
      minRA = min(data)
      maxRA = max(data)
      return(list(temptaxo = name, meanRA = meanRA, sdRA = sdRA, minRA = minRA, maxRA = maxRA ))
        })
  
  # case with group
  }else{
    meta <- sample_data(physeq)
    group <- meta[[GroupBy]]
    sumdf <- apply(df, 1, FUN = function(x){
      name <- vector()
      meanRA <- vector()
      temp <- vector()
      sdRA <- vector()
      minRA <- vector()
      maxRA <- vector()
      for(g in unique(group)){
        name <-  append(x = name,x[Rank])
        ind_group <- which(meta[[GroupBy]] == g)
        data <- as.numeric(x[-length(x)])
        data <- data[ind_group]
        meanRA = append(meanRA,mean(data))
        temp <- append(temp, g)
        sdRA = append(sdRA, sd(data))
        minRA = append(minRA, min(data))
        maxRA = append(maxRA, max(data))
      }
    
      return(list(temptaxo = name, temp = temp, meanRA = meanRA, sdRA = sdRA, minRA = minRA, maxRA = maxRA ))
    })
  }
  
  # convert list to a data.frame
  sumdf = do.call(rbind, sumdf)
  sumdf = data.frame(apply(sumdf, 2, unlist), stringsAsFactors = FALSE)
  
  # fix column type 
  taxo_name <- as.character(sumdf$temptaxo)
  taxo_name[is.na(taxo_name)] <- "NA"
  sumdf$temptaxo <- factor(taxo_name)
  sumdf$meanRA <- as.numeric(sumdf$meanRA)
  sumdf$sdRA <- as.numeric(sumdf$sdRA)
  sumdf$minRA <- as.numeric(sumdf$minRA)
  sumdf$maxRA <- as.numeric(sumdf$maxRA)
  if(!is.null(GroupBy)){
    sumdf$temp = factor(sumdf$temp)
    names(sumdf)[names(sumdf)=='temp'] = GroupBy
  }
  names(sumdf)[names(sumdf)=='temptaxo'] = Rank
  sumdf <- sumdf %>% arrange(-meanRA)
  return(data.table(sumdf))
}