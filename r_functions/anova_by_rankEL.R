anova_by_rankEL = function(physeq, Rank="Family", Traitements="",topN=25, filename="tmp.txt"){ 
  physeq <- physeq %>%  scale_reads(round = "round")
  dat=sum_taxo(physeq,Rank) #summarize taxa choisi (ex:family) dans otu table   
  meta=data.frame(t(data.frame(t(sample_data(physeq)))))	
  
  colTrt=which(colnames(meta)==Traitements)
  if (nrow(dat)<topN) {
    topN=nrow(dat)
  } 
  dat=dat[order(-dat$meanRA),][1:topN] 
  meta[,colTrt]=as.character(meta[,colTrt])
  otu=data.frame(otu_table(physeq))
  taxa=data.frame(physeq@tax_table@.Data)
  bigt=cbind(otu,taxa) 
  pos_col=(which(colnames(bigt)==Rank)) 
  sp=as.character(unlist(dat[,1])) 
  
  
  meta[,colTrt]=as.factor(meta[,colTrt])
  lTreatment=levels(meta[,colTrt]) #liste mois dans otu table		  
  cat("\n",Rank,"Total_meanRA","sdRA","minRA","maxRA","n","\n",sep="\t")
  cat("\n",Rank,"Total_meanRA","sdRA","minRA","maxRA","n","\n",file=filename, append=F,sep="\t")
  
  # Create Table Norm_Count_OTU x Traitement 
  count_traitement_level=array(0,c(topN,length(lTreatment)))
  for (species_i in 1:length(sp)) {
    species=sp[species_i]
    #print(species)
    v=as.character(bigt[,pos_col])
    
    indexrow=which(v==species) 
    if (length(indexrow)>0) {
      
      dat2=otu[indexrow,] #lignes avec la bonne esp?ce pour chaque species
      # Creation des niveaux pour l'anova 
      #browser()
      lrow=nrow(dat2) #nb otu
      lcol=ncol(dat2) #nb rep
      fTraitement=array("",c(lrow,lcol))			
      for (i in 1:lcol) {
        # get the Treatment and Event associated with this col				
        if (all(colnames(dat2)[i]==rownames(meta))==FALSE) {
          #cat("Error. Colnames and meta-data id are not identical.\n")
          rownames(meta)=gsub("-","\\.",rownames(meta))
          #if (all(colnames(dat2)[i]==rownames(meta))==FALSE) return()
        }
        cTraitement=as.character(unlist(meta[colnames(dat2)[i]==rownames(meta),colTrt]))								
        lcTraitement=which(lTreatment==cTraitement)
        #browser()
        fTraitement[,i]=cTraitement				
        for (j in 1:lrow) {									
          count_traitement_level[species_i,lcTraitement]=count_traitement_level[species_i,lcTraitement]+dat2[j,i]
        }
        
      }
      
      
      fTraitement=as.factor(fTraitement)
      
      #print(wilcox.test(lm(unlist(dat2) ~ fTraitement*fEvent, paired=F)))
      #print(anova(lm(unlist(dat2) ~ fTraitement*fEvent)))
      normalized_otu_counts=unlist(dat2)
      cat(unlist(dat[species_i,]),length(normalized_otu_counts),"\n",sep="\t")
      cat(unlist(dat[species_i,]),length(normalized_otu_counts),"\n",sep="\t",file=filename, append=T)
      #cat(species,"\n",file=filename, append=T)
      
      df=data.frame(normalized_otu_counts,fTraitement)
      capture.output(df,file=filename,append=T)
      #capture.output(shapiro.test(normalized_otu_counts))
      if (length(normalized_otu_counts)>3&&length(normalized_otu_counts)<5000) {
        s=shapiro.test(normalized_otu_counts)
      } else {
        s=list(p.value=0.00) #Assume data is not normal
      }
      #browser()
      if (s$p.value>0.05) { 
        cat("Normal: true (shapiro) ",s$p.value,"\n",file=filename,append=T) 
        capture.output(anova(lm(normalized_otu_counts ~ fTraitement, data=df)),file=filename,append=T)
      } else {
        cat("Normal: false (shapiro) ",s$p.value,"\n",file=filename,append=T) 
        capture.output(kruskal.test(normalized_otu_counts ~ fTraitement, data=df),file=filename,append=T)
        #capture.output(adonis(normalized_otu_counts ~ fTraitement, data=df),file=filename,append=T)
      }
      
    }
  } #End top species
  ll=count_traitement_level
  #browser()
  count_traitement_level=t(apply(ll,1,rbind)) #Collapse the Event dim
  count_traitement_level=apply(count_traitement_level,2,function(x) {(x/sum(x))}) #Relative abundance
  rownames(count_traitement_level)=sp
  #Create colonne names
  #vet=c()
  #for (t in lTreatment) {
  #		vet=c(vet,paste0(as.character(e),"x",as.character(t)))
  #}
  
  colnames(count_traitement_level)=lTreatment
  
  #print(count_traitement_level)
  return (count_traitement_level)
}
