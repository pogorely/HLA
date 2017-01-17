#artifdf<-data.frame(seq=c("AAAA","GGGG","ACAA","GGTT"),reads=c(1000,500,1,2),stringsAsFactors = F)
#res<-matrix(0,ncol=nrow(artifdf),nrow=nrow(artifdf));for (i in 1:nrow(artifdf))for (j in 1:nrow(artifdf)){res[i,j]<-dbinom(x = artifdf$reads[j],size = artifdf$reads[j]+artifdf$reads[i],prob = 0.001^stringdist(artifdf$seq[i],artifdf$seq[j],method = "hamming"))}
#colSums(res)#gives the prob, that given sequence is an error. 

add_error_col<-function(amplist,err_prob=0.01){
  res<-matrix(0,ncol=nrow(amplist),nrow=nrow(amplist));
  for (i in 1:nrow(amplist))
    for (j in 1:nrow(amplist))if (i!=j){
      res[i,j]<-dbinom(x = amplist$readnumber[j],size = amplist$readnumber[j]+amplist$readnumber[i],prob = 0.1^stringdist(amplist$assembled[i],amplist$assembled[j],method = "lv"))
      }
  amplist$error_prob<-colSums(res)
  amplist
}
