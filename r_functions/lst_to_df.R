lst_to_df=function(ll) {
  levels=max(unlist(lapply(ll,length)))
  df=array(dim = c(1,levels),0)
  for (l in ll) {
    df2=array(dim=c(1,levels),NA)
    l2=unlist(l)
    for (j in 1:length(l)) df2[1,j]=l2[j]
    df=rbind(df,df2)
  }
  df=df[-1,]
  return(df)
}