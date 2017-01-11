load("HLA_base.rda")
library(igraph)
library(stringdist)
library(Biostrings)
library(data.table)
HLA_set<-DNAStringSet(HLA_base$Sequence)
#from tcr
revcomp<-function (.seq) 
{
  rc.table <- c(A = "T", T = "A", C = "G", G = "C",N="N")
  sapply(strsplit(.seq, "", T, F, T), function(l) paste0(rc.table[l][length(l):1], 
                                                         collapse = ""), USE.NAMES = F)
}
#get two fastq files as text files! update with connection - save some memory here. 
reader<-function(read1,read2){
  #supernaive reader:
  Fastq1<-readLines(read1);
  Fastq2<-readLines(read2);
  res<-data.frame(read1=Fastq1[(1:length(Fastq1))%%4==2],read2=Fastq2[(1:length(Fastq2))%%4==2],stringsAsFactors = F)
  rm(Fastq1)
  rm(Fastq2)
  res
}


gzreader<-function(read1,read2){
  #supernaive reader:
  gzread1 <- gzfile(read1, open = "r")
  gzread2 <- gzfile(read2, open = "r")
  Fastq1<-readLines(read1);
  Fastq2<-readLines(read2);
  res<-data.frame(read1=Fastq1[(1:length(Fastq1))%%4==2],read2=Fastq2[(1:length(Fastq2))%%4==2],stringsAsFactors = F)
  rm(Fastq1)
  rm(Fastq2)
  close(gzread1)
  close(gzread2)
  res
}

#merge them, sort them, filter unique in crosstab fashion. trim them. 
get_groups<-function(readed,threshold=10,read_length=301,primer_length=16){
  readed$read1<-substr(readed$read1,primer_length,read_length)#trim them
  readed$read2<-substr(readed$read2,primer_length,read_length)#trim them
  readed<-readed[duplicated(readed$read1),]
  readed<-split(readed,f = readed$read1) #merge them
  if (length(readed)>0)
  readed<-readed[(rank(-sapply(readed,nrow),ties.method = "first")<=threshold)&(sapply(readed,nrow)>2)]
  #readed<-readed[sapply(readed,nrow)>threshold] #leave only huge groups. 
  readed
}

merge_cons<-function(readed){
  consensus<-character(length(readed))
  if(length(readed)>0)
  for (i in 1:length(readed))
  {#print(i)
    consensus[i]<-docons(readed[[i]]$read2)
  }
  consensus
  data.frame(read1=names(readed),read2cons=consensus,readnumber=sapply(readed,nrow),stringsAsFactors = F)
}
#replace with graph approach? possibly
get_exact<-function(merged){
  if (nrow(merged)==0)return(merged)
  resM<-character(nrow(merged))
  for (i in 1:length(resM))
  resM[i]<-paste0(HLA_base[grepl(merged$read1[i],x = HLA_base$Sequence,fixed=TRUE),]$Allele,collapse = " ")
  merged$Exact<-resM
  merged
}
make_overlaps<-function(seq1,seq2,posthres=20){
  over1<-LCS_comp2(seq1,seq2)
  #print(over1)
  over2<-LCS_comp2(seq2,seq1)
  if (over1[1]>=over2[1]){#first overlap is higher
    if(over1[2]>posthres)paste0(substr(seq1,1,over1[2]-1),seq2)
  }
  else if(over2[2]>posthres){paste0(substr(seq2,1,over2[2]-1),seq1)}
}

make_overlaps_one_side<-function(seq1,seq2,posthres=20){
  over1<-LCS_comp2(seq1,seq2)
    if(over1[2]>posthres)
      {
      data.frame(seq=paste0(substr(seq1,1,over1[2]-1),seq2),
                 Overlap_length=over1[1],Overlap_start=over1[2], 
                 Overlap_mismatches=stringdist(substr(seq1,over1[2],nchar(seq1)),substr(seq2,1,nchar(seq1)-over1[2]+1),method = "hamming")-max(N_content(substr(seq1,over1[2],nchar(seq1))),N_content(substr(seq2,1,nchar(seq1)-over1[2]+1))),stringsAsFactors = F)
      }
}

make_overlaps_one_side_grep<-function(seq1,seq2,length=50){
  Overlapstart<-regexpr(gsub("N","[AGCT]",substr(seq2,1,length) ,fixed=T),seq1)
  Overlaplength<-nchar(seq1)-Overlapstart+1
  over1<-c(Overlaplength,Overlapstart)
  if(over1[2]!=-1)
  {
    data.frame(seq=paste0(substr(seq1,1,over1[2]-1),seq2),
               Overlap_length=over1[1],Overlap_start=over1[2], 
               Overlap_mismatches=stringdist(substr(seq1,over1[2],nchar(seq1)),substr(seq2,1,nchar(seq1)-over1[2]+1),method = "hamming")-max(N_content(substr(seq1,over1[2],nchar(seq1))),N_content(substr(seq2,1,nchar(seq1)-over1[2]+1))),stringsAsFactors = F)
  }
}

get_overlap_merged<-function(merged){
  temp<-list();
  if(nrow(merged)==0)return(merged)
  
  for (i in 1:nrow(merged))temp[[i]]<-LCS_comp2(merged$read1[i],revcomp(merged$read2cons[i]))
  temp<-do.call(rbind,temp)
  merged$Overlap_max_aligned<-temp[,1]
  merged$Overlap_start<-temp[,2]
  merged$assembled<-paste0(merged$read1,substr(revcomp(merged$read2cons),nchar(merged$read1)-merged$Overlap_start+2,nchar(merged$read2cons)))
  merged
}

get_overlap_merged_fix<-function(merged,shift=1){
  temp<-list();
  if(nrow(merged)==0)return(merged)
  #for (i in 1:nrow(merged))temp[[i]]<-LCS_comp2(merged$read1[i],revcomp(merged$read2cons[i]))
  #temp<-do.call(rbind,temp)
  merged$Overlap_max_aligned<-shift
  merged$Overlap_start<-nchar(merged$read1[1])-shift+1
  merged$assembled<-paste0(merged$read1,substr(revcomp(merged$read2cons),nchar(merged$read1)-merged$Overlap_start+2,nchar(merged$read2cons)))
  merged
}

alphabet<-c(A=1,G=2,C=3,T=4,N=5)
#HLA file Make it faster with substr
docons<-function(vec,thres=0.7){
  paste0(sapply(apply(do.call(rbind,(strsplit(vec,""))),MARGIN = 2,table),function(x){if(max(prop.table(x))>thres){names(x)[x==max(x)][1]}else{"N"}}),collapse = "")
}

docons_exact<-function(vec){
  paste0(sapply(apply(do.call(rbind,(strsplit(vec,""))),MARGIN = 2,table),function(x){names(x)[x==max(x)][1]}),collapse = "")
}

docons2<-function(vec){
  res<-(sapply(apply(do.call(rbind,(strsplit(vec,""))),MARGIN = 2,table),function(x){x[c('A','G','C','T','N','.','-','*')]}))
  row.names(res)<-c("A","G","C","T","N",".","-","*")
  res[is.na(res)]<-0
  res
}

#"SEQUENCE_ID=example","PRIMER_TASK=check_primers","PRIMER_MIN_SIZE=8","PRIMER_MAX_SIZE=36",
check_primers<-function(forw,rev){
  #setwd("primer3")
  write(paste("SEQUENCE_ID=example","PRIMER_TASK=check_primers","PRIMER_MIN_SIZE=8","PRIMER_MAX_SIZE=36","PRIMER_MIN_TM=40",paste0("SEQUENCE_PRIMER=",forw),paste0("SEQUENCE_PRIMER_REVCOMP=",rev),"PRIMER_EXPLAIN_FLAG=1","=", sep="\n"),file = "primer.txt")
  system("./primer3_core < primer.txt -output=out.txt")
  out<-do.call(rbind,strsplit(readLines("out.txt"),"="))
  out<-data.frame(Par=out[,1],res=out[,2],stringsAsFactors = F)
  row.names(out)<-out[,1]
  out
  #setwd("..")
}  
docons3<-function(vec,thres=0.8){#make 2 cons anyway... 
  constbl<-sapply(apply(do.call(rbind,(strsplit(vec,"")))[1:4,],MARGIN = 2,function(x)prop.table(table(x))),function(x){x[c('A','G','C','T','N','.','-','*')]})
  conf<-which(apply(constbl,MARGIN = 2,function(x){max(x,na.rm=T)<thres}))#get.pos with confidence<0.8?
  print(conf)
  #conf<-which(constbl2)
  print(apply(do.call(rbind,(strsplit(vec,"")))[,conf],MARGIN=1,paste0,collapse=""))
  champs<-(sort(table(apply(do.call(rbind,(strsplit(vec,"")))[,conf],MARGIN=1,paste0,collapse="")),decreasing=T))[1:20]
  cons<-docons(vec)
  champs
  #paste them and get 2 champions
  #return champions with read_numbers
  #2cons with read numbers. 
}

intersect_amplicones<-function(exact_list){
  res<-list()
  for (i in 1:nrow(exact_list)){res[[i]]<-list();#print(i)
    for (j in 1:nrow(exact_list))
      if (i<j)
        if (sum(duplicated(c(unlist(strsplit(as.character(exact_list[i,]$Exact)," ")),unlist(strsplit(as.character(exact_list[j,]$Exact)," ")))))!=0)
          {
           res[[i]][[j]]<-make_overlaps(exact_list[i,]$assembled,exact_list[j,]$assembled)
        }
  }
  res
}

intersect_amplicones2<-function(exact_list){
  res<-list()
  for (i in 1:nrow(exact_list)){res[[i]]<-list();#print(i)
  for (j in 1:nrow(exact_list))
    if (i<j)
      if (sum(duplicated(c(unlist(strsplit(as.character(exact_list[i,]$Exact)," ")),unlist(strsplit(as.character(exact_list[j,]$Exact)," ")))))!=0)
      {
        over<-make_overlaps(exact_list[i,]$assembled,exact_list[j,]$assembled)
        if(length(over)!=0)
        res[[i]][[j]]<-data.frame(over,first_reads=exact_list[i,]$readnumber,second_reads=exact_list[j,]$readnumber,first_assembled=exact_list[i,]$assembled,second_assembled=exact_list[j,]$assembled)
      }
  }
  res
}
intersect_amplicones3<-function(first_amp_list,second_amp_list){
  #res<-list()
  res<-data.frame()
  if ((nrow(first_amp_list)==0)|(nrow(second_amp_list)==0)) return (res)
  for (i in 1:nrow(first_amp_list)){#res[[i]]<-list();#print(i)
  for (j in 1:nrow(second_amp_list))
   # if (i<j)
      #if (sum(duplicated(c(unlist(strsplit(as.character(exact_list[i,]$Exact)," ")),unlist(strsplit(as.character(exact_list[j,]$Exact)," ")))))!=0)
      {
        over<-make_overlaps_one_side_grep(first_amp_list[i,]$assembled,second_amp_list[j,]$assembled)
        if(length(over)!=0){
          over$first_reads=first_amp_list[i,]$readnumber;
          over$second_reads=second_amp_list[j,]$readnumber;
          over$first_reads_freq=first_amp_list[i,]$freq;
          over$second_reads_freq=second_amp_list[j,]$freq;
          over$mean_geom_freq=sqrt(first_amp_list[i,]$freq*second_amp_list[j,]$freq)
          over$geom_mean_reads=sqrt(first_amp_list[i,]$readnumber*second_amp_list[j,]$readnumber)
          over$geom_min_coverage=min(first_amp_list[i,]$readnumber,second_amp_list[j,]$readnumber)
          over$first_assembled=first_amp_list[i,]$assembled
          over$second_assembled=second_amp_list[j,]$assembled
          over$Ns=sapply(over$seq,N_content)
          over$Length=nchar(over$seq)
          over$Exact_match<-paste0(HLA_base[grepl(gsub("N","[AGCT]",over$seq,fixed=T),x = HLA_base$Sequence),]$Allele,collapse = " ")
          over$One_mismatch<-paste0(HLA_base[vcountPattern(pattern = over$seq,subject = HLA_set,algorithm = "auto",fixed = F,max.mismatch = 1,min.mismatch = 1,with.indels = F)!=0,]$Allele,collapse = " ")
          over$Two_mismatch<-paste0(HLA_base[vcountPattern(pattern = over$seq,subject = HLA_set,algorithm = "auto",fixed = F,max.mismatch = 2,min.mismatch = 2,with.indels = F)!=0,]$Allele,collapse = " ")
          res<-rbind(res,over)
          #res[[i]][[j]]<-over
        }
      }
  }
  res
  #do.call(rbind,do.call(rbind,res))
}

mismatch_vector<-function(seq,err_prob=0.001,max_mismatch=3){ #returns a vector???
  mism<-rep(0.0000000001,nrow(HLA_base)) #set very low prob 
  mism[grepl(gsub("N","[AGCT]",seq,fixed=T),x = HLA_base$Sequence)]<-1 #if something is exactly matched it gets weight 1
  for (i in 1:max_mismatch)
    mism[vcountPattern(pattern = seq,subject = HLA_set,algorithm = "auto",fixed = F,max.mismatch = i,min.mismatch = i,with.indels = F)!=0]<-err_prob**i #if something is not exactly matched it gets weight err_prob**i
  #mism<-mism/sum(mism,na.rm = T)
  prop.table(mism)
}

mismatch_matrix<-function(amplist,err_prob=0.001,max_mismatch=3){
  do.call(cbind,lapply(amplist$assembled,mismatch_vector,err_prob=0.001,max_mismatch=3))
}

Alleles_bayes<-function(amp1list, amp2list){
#get mismat.

#fi get fi
amp1mism<-mismatch_matrix(amp1list)
amp2mism<-mismatch_matrix(amp2list)
amp1maps<-(apply(amp1mism,MARGIN = 1,function(x){paste0(which(x>0.001),collapse=" ")}))
amp2maps<-(apply(amp2mism,MARGIN = 1,function(x){paste0(which(x>0.001),collapse=" ")}))
fi<-amp1list$freq%*%t(amp1mism) #why is that? it could be replaced with some prob, that given seq is error from bigger seq (matrix n_seq by n_seq), ratios? and then we simply multiply each by rowsums?== integrate over n2
bi<-amp2list$freq%*%t(amp2mism) #PROB THAT IT IS TRUE SEQUENCE. 

#list(as.vector(fi),as.vector(bi),as.vector(fi*bi))
as.data.table(data.frame(fi=as.vector(fi),bi=as.vector(bi),fibi=as.vector(fi)*as.vector(bi),allele=HLA_base$Allele,class=HLA_base$HLA_class,a1=amp1maps,a2=amp2maps,a1_a2=paste(amp1maps,amp2maps,sep="_")))
#get mi.
#tstluk[class=="C",,][prop.table(fibi)>0.01,,][,list(amps=paste0(allele,collapse=" "),fibi=sum(fibi)),a1_a2]
}

bayesian_alleles_report<-function(allelelist){
  a<-allelelist[class=="A",,][,,][,list(amps=paste0(allele,collapse=" "),fibi=sum(fibi)),a1_a2][,rel_fibi:=prop.table(fibi),][order(-rel_fibi)]
  b<-allelelist[class=="B",,][,,][,list(amps=paste0(allele,collapse=" "),fibi=sum(fibi)),a1_a2][,rel_fibi:=prop.table(fibi),][order(-rel_fibi)]
  c<-allelelist[class=="C",,][,,][,list(amps=paste0(allele,collapse=" "),fibi=sum(fibi)),a1_a2][,rel_fibi:=prop.table(fibi),][order(-rel_fibi)]
  DRB<-allelelist[grepl("DRB",class),,][,,][,list(amps=paste0(allele,collapse=" "),fibi=sum(fibi)),a1_a2][,rel_fibi:=prop.table(fibi),][order(-rel_fibi)]
  DQB<-allelelist[class=="DQB1",,][,,][,list(amps=paste0(allele,collapse=" "),fibi=sum(fibi)),a1_a2][,rel_fibi:=prop.table(fibi),][order(-rel_fibi)]
  list(A=a,B=b,C=c,DRB=DRB,DQB=DQB)
}

LCS_comp2<-function(a,b){#one side overlap a b
  a<-alphabet[unlist(strsplit(a,""))]
  b<-alphabet[unlist(strsplit(b,""))]
  res<-a
  for (i in 1:length(a))
  {
    res[i]<-suppressWarnings(max(rle(a[i:length(a)]-b[1:(length(a)-i+1)])$lengths[rle(a[i:length(a)]-b[1:(length(a)-i+1)])$values==0],na.rm=T))
   }
  c(max(res),which(res==max(res)))
}
#LCS_comp<-function(a,b)paste0(LCS(unlist(strsplit(a,"")),unlist(strsplit(b,"")))$va,collapse="")

Run_pipeline<-function(read1,read2,threshold=10,read_length=301){
  print("File read start")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  #readed<-reader(read1,read2)
  readed<-gzreader(read1,read2)
  #readed<-readed[sample(1:nrow(readed),10000),]
  print("File read end. Looking for champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-merge_cons(get_groups(readed,threshold = threshold,read_length = read_length))
  print("We are the champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-get_overlap_merged(readed)
  print("Reads_assmbled")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  #write.table
  readed<-get_exact(readed)
  print("Exact matches read1 found")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed_exact<-readed[readed$Exact!="",]
  Amp<-intersect_amplicones(readed_exact)
  print("Amplicones intersected")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  Amplist<-unlist(Amp)
  resM<-character(length(Amplist))
  
    if(length(Amplist>0)){
  for (i in 1:length(resM))
    resM[i]<-paste0(HLA_base[grepl(gsub("N","[AGCT]",Amplist[i],fixed=T),x = HLA_base$Sequence),]$Allele,collapse = " ")
  }
  print("Amplicones exact found")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  list(AMP=Amp,safety4=table(resM),safety3=data.frame(Full_seq=Amplist,Exact_search=resM,stringsAsFactors = F),safety2=readed_exact,safety1=readed)
}

Run_pipeline_first_class<-function(read1,read2,threshold=10,read_length=301){
  print("File read start")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  #readed<-reader(read1,read2)
  readed<-gzreader(read1,read2)
  readed<-readed[grepl("CCCTGACC[GC]AGACCTG|CGACGGCAA[AG]GATTAC",readed$read1),]
  #readed<-readed[sample(1:nrow(readed),10000),]
  print("File read end. Looking for champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-merge_cons(get_groups(readed,threshold = threshold,read_length = read_length))
  print("We are the champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-get_overlap_merged(readed)
  print("Reads_assmbled")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  #write.table
  readed<-get_exact(readed)
  print("Exact matches read1 found")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed_exact<-readed[readed$Exact!="",]
  Amp<-intersect_amplicones(readed_exact)
  print("Amplicones intersected")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  Amplist<-unlist(Amp)
  resM<-character(length(Amplist))
  
  if(length(Amplist>0)){
    for (i in 1:length(resM))
      resM[i]<-paste0(HLA_base[grepl(gsub("N","[AGCT]",Amplist[i],fixed=T),x = HLA_base$Sequence),]$Allele,collapse = " ")
  }
  print("Amplicones exact found")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  list(AMP=Amp,safety4=table(resM),safety3=data.frame(Full_seq=Amplist,Exact_search=resM,stringsAsFactors = F),safety2=readed_exact,safety1=readed)
}

make.sequence.graph <- function (.data, .name = '',max_errs=1) {
  G <- graph.empty(n = length(.data), directed=F)
  tmp<-stringdistmatrix(.data,.data,method="hamming")
  G <- add.edges(G, t(which(tmp==max_errs,arr.ind=T)))
  G <- igraph::simplify(G)
  G <- set.vertex.attribute(G, 'label', V(G), .data)
  #print(G)
  G
}
N_content<-function(str){
  table(c("N",unlist(strsplit(str,""))))["N"]-1
}
fathers_and_children<-function(readed){#graph and cluster analysis
  if (nrow(readed)==0)return(readed)
  Gr<-make.sequence.graph(readed$read1)
  readed$neighbours_read1<-degree(Gr)
  cl<-clusters(Gr)
  readed$clusters<-cl$membership
  readed$freq<-prop.table(readed$readnumber)
  parents<-numeric(nrow(readed))
  for (clust in unique(cl$membership))
    parents[readed$clusters==clust]<-rank(-readed[readed$clusters==clust,]$neighbours_read1,ties.method = "min")
  readed$parents<-parents
  readed
}
primer_list_class1_reg<-c(HLA_I_amp1_next_For="CCCTGACC[GC]AGACCTG",HLA_I_amp1_next_Rev2="AGAG[CT]CTACCTGGA[GT]GG",HLA_I_amp2_next_For="CGACGGCAA[AG]GATTAC",HLA_I_amp2_next_Rev="TCTTCCCAG[AG]CCACC[AG]")#,Iamp1rev="TCTTCCCAGACCACCA",Iamp1rev="TCTTCCCAGGCCACCG",Iamp1rev="TCTTCCCAGACCACCG",Iamp1rev="TCTTCCCAGGCCACCA")

Run_pipeline_novel<-function(read1,read2,threshold=10,read_length=301){
  print("File read start")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  #readed<-reader(read1,read2)
  readed<-gzreader(read1,read2)
  #Iamp1log<-grepl("CCCTGACC[GC]AGACCTG",readed$read1)
  inv<-grepl("[CT]GGTGG[AG]CTGGGAAGA|[CT]CAGCAGGTTGTGGTG|CCAG[GC]AGGTT[AG]TGGTG|CCAC[GT]TGGCAGGTGTA|CCACTTGGCAAGTGTA",substr(readed$read1,1,20))
  readed[inv,c(1,2)]<-readed[inv,c(2,1)]
  readed<-list(Iamp1=readed[grepl("CCCTGACC[GC]AGACCTG",substr(readed$read1,1,20)),],
               Iamp2=readed[grepl("CGACGGCAA[AG]GATTAC",substr(readed$read1,1,20)),],
               IIamp1_DQB=readed[grepl("AG[GT]CTTTGCGGATCCC",substr(readed$read1,1,20)),],
               IIamp1_Others=readed[grepl("CTGAGCTCCC[GC]ACTGG",substr(readed$read1,1,20)),],
               IIamp2=readed[grepl("GGAACAGCCAGAAGGA",substr(readed$read1,1,20)),])
  print(sapply(readed,nrow))
  if (sum(sapply(readed,nrow)==0)>0)return(list(safety1="fail", safety2="fail", Iclass="fail", IIclassDQB="fail", IIclassOthers="fail", 
                                                ClassI_report="fail", ClassII_DQB_report="fail", ClassII_Others_report="fail"))
  #readed<-readed[sample(1:nrow(readed),10000),]
  print("File read end. Looking for champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-lapply(lapply(readed,get_groups,threshold = threshold,read_length = read_length),merge_cons)
  print("We are the champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  print(sapply(readed,nrow))
  if (sum(sapply(readed,nrow)==0)>0)return(list(safety1=readed, safety2="fail", Iclass="fail", IIclassDQB="fail", IIclassOthers="fail", 
                                                 ClassI_report="fail", ClassII_DQB_report="fail", ClassII_Others_report="fail"))
  readed<-lapply(readed,get_overlap_merged)
  print("Reads_assembled")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-lapply(readed,get_exact)
  print("Exact matches read1 found")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-lapply(readed,fathers_and_children)
  #readed_intersect<-lapply(readed,function(x){x[x$parents==1|x$parents==2,]})#||rank(-x$readnumber,ties.method="first")<15|x$Exact!=""
  readed_intersect<-lapply(readed,function(x){x[x$parents==1,]})
  res<-list(safety1=readed,safety2=readed_intersect,Iclass=do.call(rbind,do.call(rbind,(intersect_amplicones3(readed_intersect$Iamp1,readed_intersect$Iamp2)))),
                                              IIclassDQB=do.call(rbind,do.call(rbind,(intersect_amplicones3(readed_intersect$IIamp1_DQB,readed_intersect$IIamp2)))),
                                              IIclassOthers=do.call(rbind,do.call(rbind,(intersect_amplicones3(readed_intersect$IIamp1_Others,readed_intersect$IIamp2)))))
  res$ClassI_report<-data.frame(Iclass_exact=names(table(res$Iclass$Exact_match)),Iclass_exact_num=table(res$Iclass$Exact_match))
  res$ClassII_DQB_report<-data.frame(IIclassDQB_exact=names(table(res$IIclassDQB$Exact_match)),IIclassDQB_exact_num=table(res$IIclassDQB$Exact_match))
  res$ClassII_Others_report<-data.frame(IIclassOthersclass_exact=names(table(res$IIclassOthers$Exact_match)),IIclassOthers_exact_num=table(res$IIclassOthers$Exact_match))
  print("Amplicones intersected. Report generated.")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  res
}

Run_pipeline_novel_alternative<-function(read1,read2,threshold=10,read_length=301){
  print("File read start")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  #readed<-reader(read1,read2)
  readed<-gzreader(read1,read2)
  #Iamp1log<-grepl("CCCTGACC[GC]AGACCTG",readed$read1)
  #inv<-grepl("[CT]GGTGG[AG]CTGGGAAGA|[CT]CAGCAGGTTGTGGTG|CCAG[GC]AGGTT[AG]TGGTG|CCAC[GT]TGGCAGGTGTA|CCACTTGGCAAGTGTA",substr(readed$read1,1,20))
  #readed[inv,c(1,2)]<-readed[inv,c(2,1)]
  readed<-list(Iamp1alt=readed[grepl("[CT]CACTCCATGAGGTATTTC|CCACTCCATGAAGTATTTC",substr(readed$read1,1,22)),],
               Iamp2alt=readed[grepl("GGCAA[AG]GATTACATCGCC|GGCAAGGATTACATCGCT",substr(readed$read1,1,22)),],
               IIamp1_DQBalt=readed[grepl("TGAGGGCAGAGAC[CT]CTCC",substr(readed$read1,1,22)),],
 #              IIamp1_DRBalt=readed[grepl("TGACAGTGACACTGATGG|TGACAGTGACATTGACGG",substr(readed$read1,1,22)),],
 #              IIamp2_DRBalt=readed[grepl("GGTTTCTATCCAGGCAGC",substr(readed$read1,1,22))&grepl("CCAGAGTGTCCTTTCTGA",substr(readed$read2,1,22)),],
               IIamp2_DQBalt=readed[grepl("ACCATCTCCCCATCCAG",substr(readed$read1,1,22))&grepl("TG[CT]TCTGGGCAGATTCAG",substr(readed$read2,1,22)),]) 

  print(sapply(readed,nrow))
  if (sum(sapply(readed,nrow)==0)>0)return(list(safety1="fail", safety2="fail", Iclass="fail", IIclassDQB="fail", IIclassOthers="fail", 
                                                ClassI_report="fail", ClassII_DQB_report="fail", ClassII_Others_report="fail"))
  #readed<-readed[sample(1:nrow(readed),10000),]
  print("File read end. Looking for champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-lapply(lapply(readed,get_groups,threshold = threshold,read_length = read_length,primer_length=22),merge_cons)
  print("We are the champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  print(sapply(readed,nrow))
  #if (sum(sapply(readed,nrow)==0)>0)return(list(safety1=readed, safety2="fail", Iclass="fail", IIclassDQB="fail", IIclassOthers="fail", 
  #                                              ClassI_report="fail", ClassII_DQB_report="fail", ClassII_Others_report="fail"))
  readed<-lapply(readed,get_overlap_merged)
  readed[[1]]<-get_overlap_merged_fix(readed[[1]])
  print("Reads_assembled")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-lapply(readed,get_exact)
  print("Exact matches read1 found")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-lapply(readed,fathers_and_children)
  #readed_intersect<-lapply(readed,function(x){x[x$parents==1|x$parents==2,]})#||rank(-x$readnumber,ties.method="first")<15|x$Exact!=""
  readed_intersect<-lapply(readed,function(x){x[x$parents==1,]})
  res<-list(safety1=readed,safety2=readed_intersect,Iclass=do.call(rbind,do.call(rbind,(intersect_amplicones3(readed_intersect$Iamp1alt,readed_intersect$Iamp2alt)))),
            IIclassDQB=do.call(rbind,do.call(rbind,(intersect_amplicones3(readed_intersect$IIamp1_DQBalt,readed_intersect$IIamp2_DQBalt))))  ) #,
         #   IIclassOthers=do.call(rbind,do.call(rbind,(intersect_amplicones3(readed_intersect$IIamp1_DRBalt,readed_intersect$IIamp2_DRBalt)))))
  res$ClassI_report<-data.frame(Iclass_exact=names(table(res$Iclass$Exact_match)),Iclass_exact_num=table(res$Iclass$Exact_match))
  res$ClassII_DQB_report<-data.frame(IIclassDQB_exact=names(table(res$IIclassDQB$Exact_match)),IIclassDQB_exact_num=table(res$IIclassDQB$Exact_match))
 # res$ClassII_Others_report<-data.frame(IIclassOthersclass_exact=names(table(res$IIclassOthers$Exact_match)),IIclassOthers_exact_num=table(res$IIclassOthers$Exact_match))
  print("Amplicones intersected. Report generated.")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  res
}

File_list_pipeline<-function(filelist,read_length=300,threshold=100){
  resHLA<-lapply(filelist[,1],"[",1)
  names(resHLA)<-filelist[,1]
  for (i in 1:nrow(filelist)){
  print(names(resHLA)[i])
  resHLA[[i]]<-Run_pipeline_novel(filelist[i,2],filelist[i,3],read_length = read_length,threshold=threshold)
  if(resHLA[[i]]$safety2!="fail"){
  write.table(resHLA[[i]]$safety1$Iamp1[order(-resHLA[[i]]$safety1$Iamp1$readnumber),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_Iclass_amp1_safety1.csv",collapse = ""))
  write.table(resHLA[[i]]$safety1$Iamp2[order(-resHLA[[i]]$safety1$Iamp2$readnumber),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_Iclass_amp2_safety1.csv",collapse = ""))
  write.table(resHLA[[i]]$safety1$IIamp1_DQB[order(-resHLA[[i]]$safety1$IIamp1_DQB$readnumber),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_IIclass_DQB_amp1_safety1.csv",collapse = ""))
  write.table(resHLA[[i]]$safety1$IIamp1_Others[order(-resHLA[[i]]$safety1$IIamp1_Others$readnumber),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_IIclass_Others_safety1.csv",collapse = ""))
  write.table(resHLA[[i]]$safety1$IIamp2[order(-resHLA[[i]]$safety1$IIamp2$readnumber),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_IIclass_amp2_safety1.csv",collapse = ""))
  }
  if(resHLA[[i]]$Iclass!="fail"){
  write.table(resHLA[[i]]$Iclass[order(resHLA[[i]]$Iclass$Overlap_mismatches,-resHLA[[i]]$Iclass$geom_min_coverage),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_Iclass_Full.csv",collapse = ""))
  write.table(resHLA[[i]]$IIclassDQB[order(resHLA[[i]]$IIclassDQB$Overlap_mismatches,-resHLA[[i]]$IIclassDQB$geom_min_coverage),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_IIclass_Full_DQB.csv",collapse = ""))
  write.table(resHLA[[i]]$IIclassOthers[order(resHLA[[i]]$IIclassOthers$Overlap_mismatches,-resHLA[[i]]$IIclassOthers$geom_min_coverage),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_IIclass_Full_Others.csv",collapse = ""))
  }
  }
  resHLA
}

File_list_pipeline_alternative<-function(filelist,read_length=300,threshold=100){
  resHLA<-lapply(filelist[,1],"[",1)
  names(resHLA)<-filelist[,1]
  for (i in 1:nrow(filelist)){
    print(names(resHLA)[i])
    resHLA[[i]]<-Run_pipeline_novel_alternative(filelist[i,2],filelist[i,3],read_length = read_length,threshold=threshold)
    if(resHLA[[i]]$safety2!="fail"){
      write.table(resHLA[[i]]$safety1$Iamp1alt[order(-resHLA[[i]]$safety1$Iamp1alt$readnumber),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_Iclass_amp1_safety1_alt.csv",collapse = ""))
      write.table(resHLA[[i]]$safety1$Iamp2alt[order(-resHLA[[i]]$safety1$Iamp2alt$readnumber),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_Iclass_amp2_safety1_alt.csv",collapse = ""))
      #write.table(resHLA[[i]]$safety1$IIamp1_DQB[order(-resHLA[[i]]$safety1$IIamp1_DQB$readnumber),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_IIclass_DQB_amp1_safety1.csv",collapse = ""))
      #write.table(resHLA[[i]]$safety1$IIamp1_Others[order(-resHLA[[i]]$safety1$IIamp1_Others$readnumber),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_IIclass_Others_safety1.csv",collapse = ""))
      #write.table(resHLA[[i]]$safety1$IIamp2[order(-resHLA[[i]]$safety1$IIamp2$readnumber),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_IIclass_amp2_safety1.csv",collapse = ""))
    }
    if(resHLA[[i]]$Iclass!="fail"){
      write.table(resHLA[[i]]$Iclass[order(resHLA[[i]]$Iclass$Overlap_mismatches,-resHLA[[i]]$Iclass$geom_min_coverage),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_Iclass_Full_alt.csv",collapse = ""))
      #write.table(resHLA[[i]]$IIclassDQB[order(resHLA[[i]]$IIclassDQB$Overlap_mismatches,-resHLA[[i]]$IIclassDQB$geom_min_coverage),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_IIclass_Full_DQB.csv",collapse = ""))
      #write.table(resHLA[[i]]$IIclassOthers[order(resHLA[[i]]$IIclassOthers$Overlap_mismatches,-resHLA[[i]]$IIclassOthers$geom_min_coverage),],sep=";",quote = F,row.names = F,file = paste0(names(resHLA)[i],"_IIclass_Full_Others.csv",collapse = ""))
    }
  }
  resHLA
}

exclusive<-function(Zero_list){
  
  #for (i in 1:nrow(Zero_list))
  Zero_list<-Zero_list[!duplicated(Zero_list$first_assembled),];
  Zero_list<-Zero_list[!duplicated(Zero_list$second_assembled),];
  Zero_list
}
report_producer<-function(res_list,name){
  #sort all by mismatches and frequency.
  classI<-res_list$Iclass[order(res_list$Iclass$Overlap_mismatches,-res_list$Iclass$mean_geom_freq),]
  DQB<-res_list$IIclassDQB[order(res_list$IIclassDQB$Overlap_mismatches,-res_list$IIclassDQB$mean_geom_freq),]
  DRB<-res_list$IIclassOthers[order(res_list$IIclassOthers$Overlap_mismatches,-res_list$IIclassOthers$mean_geom_freq),]
  #
  classI<-exclusive(classI[classI$Overlap_mismatches==0,])
  DRB<-exclusive(DRB[DRB$Overlap_mismatches==0,])
  DQB<-exclusive(DQB[DQB$Overlap_mismatches==0,])
  #
  classIgr<-paste0(gsub(" ","|",classI$Exact_match[classI$Exact_match!=""],fixed=T),collapse = "|")
  classIgr<-gsub("*","\\*",classIgr,fixed = T)
  
  Iclass_Iamp<-res_list$safety1$Iamp1[!grepl(classIgr,res_list$safety1$Iamp1$Exact),]
  Iclass_Iamp<-Iclass_Iamp[Iclass_Iamp$Exact!="",]
  Iclass_IIamp<-res_list$safety1$Iamp2[!grepl(classIgr,res_list$safety1$Iamp2$Exact),]
  Iclass_IIamp<-Iclass_IIamp[Iclass_IIamp$Exact!="",]
  
  classI<-classI[order(classI$Exact_match),]
  if(nrow(Iclass_Iamp)!=0)classI<-rbind(classI,data.frame(seq="1 AMP ONLY",Overlap_length=Iclass_Iamp$parent,Overlap_start="",
                          Overlap_mismatches="",first_reads=Iclass_Iamp$readnumber,
                          second_reads="",first_reads_freq=Iclass_Iamp$freq,second_reads_freq="",
                          mean_geom_freq="",geom_mean_reads="",geom_min_coverage=Iclass_Iamp$readnumber,
                          first_assembled=Iclass_Iamp$assembled,second_assembled="",Ns="",Length="",Exact_match=Iclass_Iamp$Exact,
                          One_mismatch="",Two_mismatch=""))
  if(nrow(Iclass_IIamp)!=0)
    classI<-rbind(classI,data.frame(seq="2 AMP ONLY",Overlap_length=Iclass_IIamp$parent,Overlap_start="",
                   Overlap_mismatches="",first_reads="",
                   second_reads=Iclass_IIamp$readnumber,first_reads_freq="",second_reads_freq=Iclass_IIamp$freq,
                   mean_geom_freq="",geom_mean_reads="",geom_min_coverage=Iclass_IIamp$readnumber,
                   first_assembled="",second_assembled=Iclass_IIamp$assembled,Ns="",Length="",Exact_match=Iclass_IIamp$Exact,
                   One_mismatch="",Two_mismatch=""))
  #Go from top to bottom, read_exact, write results
  classI$guy=name
  classI<-cbind(classI[,c(19,16)],classI[,-c(16,19)])
  
  DRBgr<-paste0(gsub(" ","|",DRB$Exact_match[DRB$Exact_match!=""],fixed=T),collapse = "|")
  DRBgr<-gsub("*","\\*",DRBgr,fixed = T)
  #print(DRBgr)
  DRB_Iamp<-res_list$safety1$IIamp1_Others[!grepl(DRBgr,res_list$safety1$IIamp1_Others$Exact),]
  DRB_Iamp<-DRB_Iamp[DRB_Iamp$Exact!="",]
  DRB_IIamp<-res_list$safety1$IIamp2[!grepl(DRBgr,res_list$safety1$IIamp2$Exact),]
  DRB_IIamp<-DRB_IIamp[DRB_IIamp$Exact!="",]
  #print(DRB_Iamp)
  DRB<-DRB[order(DRB$Exact_match),]
  if(nrow(DRB_Iamp)!=0)DRB<-rbind(DRB,data.frame(seq="1 AMP ONLY",Overlap_length=DRB_Iamp$parent,Overlap_start="",
                                                 Overlap_mismatches="",first_reads=DRB_Iamp$readnumber,
                                                 second_reads="",first_reads_freq=DRB_Iamp$freq,second_reads_freq="",
                                                 mean_geom_freq="",geom_mean_reads="",geom_min_coverage=DRB_Iamp$readnumber,
                                                 first_assembled=DRB_Iamp$assembled,second_assembled="",Ns="",Length="",Exact_match=DRB_Iamp$Exact,
                                                 One_mismatch="",Two_mismatch=""))
  if(nrow(DRB_IIamp)!=0)
    DRB<-rbind(DRB,data.frame(seq="2 AMP ONLY",Overlap_length=DRB_IIamp$parent,Overlap_start="",
                              Overlap_mismatches="",first_reads="",
                              second_reads=DRB_IIamp$readnumber,first_reads_freq="",second_reads_freq=DRB_IIamp$freq,
                              mean_geom_freq="",geom_mean_reads="",geom_min_coverage=DRB_IIamp$readnumber,
                              first_assembled="",second_assembled=DRB_IIamp$assembled,Ns="",Length="",Exact_match=DRB_IIamp$Exact,
                              One_mismatch="",Two_mismatch=""))
  #Go from top to bottom, read_exact, write results
  DRB$guy=name
  DRB<-cbind(DRB[,c(19,16)],DRB[,-c(16,19)])
  
  DQBgr<-paste0(gsub(" ","|",DQB$Exact_match[DQB$Exact_match!=""],fixed=T),collapse = "|")
  DQBgr<-gsub("*","\\*",DQBgr,fixed = T)
  #print(DQBgr)
  DQB_Iamp<-res_list$safety1$IIamp1_DQB[!grepl(DQBgr,res_list$safety1$IIamp1_DQB$Exact),]
  DQB_Iamp<-DQB_Iamp[DQB_Iamp$Exact!="",]
  DQB_IIamp<-res_list$safety1$IIamp2[!grepl(DQBgr,res_list$safety1$IIamp2$Exact),]
  DQB_IIamp<-DQB_IIamp[DQB_IIamp$Exact!="",]
  #print(DQB_Iamp)
  DQB<-DQB[order(DQB$Exact_match),]
  if(nrow(DQB_Iamp)!=0)DQB<-rbind(DQB,data.frame(seq="1 AMP ONLY",Overlap_length=DQB_Iamp$parent,Overlap_start="",
                                                 Overlap_mismatches="",first_reads=DQB_Iamp$readnumber,
                                                 second_reads="",first_reads_freq=DQB_Iamp$freq,second_reads_freq="",
                                                 mean_geom_freq="",geom_mean_reads="",geom_min_coverage=DQB_Iamp$readnumber,
                                                 first_assembled=DQB_Iamp$assembled,second_assembled="",Ns="",Length="",Exact_match=DQB_Iamp$Exact,
                                                 One_mismatch="",Two_mismatch=""))
  if(nrow(DQB_IIamp)!=0)
    DQB<-rbind(DQB,data.frame(seq="2 AMP ONLY",Overlap_length=DQB_IIamp$parent,Overlap_start="",
                              Overlap_mismatches="",first_reads="",
                              second_reads=DQB_IIamp$readnumber,first_reads_freq="",second_reads_freq=DQB_IIamp$freq,
                              mean_geom_freq="",geom_mean_reads="",geom_min_coverage=DQB_IIamp$readnumber,
                              first_assembled="",second_assembled=DQB_IIamp$assembled,Ns="",Length="",Exact_match=DQB_IIamp$Exact,
                              One_mismatch="",Two_mismatch=""))
  #Go from top to bottom, read_exact, write results
  DQB$guy=name
  DQB<-cbind(DQB[,c(19,16)],DQB[,-c(16,19)])
  
  list(classI,DRB,DQB,Iclass_Iamp,Iclass_IIamp)

}

get_forwards<-function(Alignment, pos,length=16){
  pr<-substr(Alignment,pos,pos+length-1)
  pr<-pr[!grepl("*",pr,fixed=T)]
  unique(pr)
}
get_forwards_tbl<-function(Alignment, pos,length=16){
  pr<-substr(Alignment,pos,pos+length-1)
  pr<-pr[!grepl("*",pr,fixed=T)]
  table(pr)
}

get_unmatched_first_class<-function(pos,length=16){
  primersinpos<-unique(unlist(lapply(ClassIlist,"[[",pos)))
  c(AlignHLAspf$A[!is.element(substr(AlignHLAspf$A$Alignment,pos,pos+length-1),primersinpos),1],
  AlignHLAspf$B[!is.element(substr(AlignHLAspf$B$Alignment,pos,pos+length-1),primersinpos),1],
  AlignHLAspf$C[!is.element(substr(AlignHLAspf$C$Alignment,pos,pos+length-1),primersinpos),1])
}

check_primers_list<-function(forwards,reverses){
  TMfor<-numeric(length(forwards))
  TMrev<-numeric(length(reverses))
  DIMERselffor<-numeric(length(forwards))
  DIMERselfrev<-numeric(length(reverses))
  HAIRselffor<-numeric(length(forwards))
  HAIRselfrev<-numeric(length(reverses))
  DIMER_other<-matrix(0,nrow=length(forwards),ncol=length(reverses))
  for (i in 1:length(forwards))
    for (j in 1:length(reverses))
    {
      report<-check_primers(forwards[i],reverses[j])
     # print(report)
      
      DIMER_other[i,j]<-as.numeric(report["PRIMER_PAIR_0_COMPL_ANY_TH",2])
      TMfor[i]<-as.numeric(report["PRIMER_LEFT_0_TM",2])
      TMrev[j]<-as.numeric(report["PRIMER_RIGHT_0_TM",2])
      DIMERselffor[i]<-as.numeric(report["PRIMER_LEFT_0_SELF_ANY_TH",2])
      DIMERselfrev[j]<-as.numeric(report["PRIMER_RIGHT_0_SELF_ANY_TH",2])
      HAIRselffor[i]<-as.numeric(report["PRIMER_LEFT_0_HAIRPIN_TH",2])
      HAIRselfrev[j]<-as.numeric(report["PRIMER_RIGHT_0_HAIRPIN_TH",2])
    }
  list(TMfor=TMfor,TMrev=TMrev,DIMERselffor=DIMERselffor,DIMERselfrev=DIMERselfrev,HAIRselffor=HAIRselffor,HAIRselfrev=HAIRselfrev,DIMER_other=DIMER_other)
}

get_forwards_prop_tbl<-function(Alignment, pos,length=16){
  pr<-substr(Alignment,pos,pos+length-1)
  pr<-pr[!grepl("*",pr,fixed=T)]
  prop.table(table(pr))
}
get_reverse<-function(){
  
}

Run_pipeline_short<-function(read1,read2,threshold=300,read_length=250){
  print("File read start")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  #readed<-reader(read1,read2)
  readed<-gzreader(read1,read2)
  #readed<-readed[sample(1:nrow(readed),10000),]
  print("File read end. Looking for champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-merge_cons(get_groups(readed,threshold = threshold,read_length = read_length,primer_length = 17))
  print("We are the champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  #readed<-get_overlap_merged(readed)
  #print("Reads_assmbled")
  #print(format(Sys.time(), "%a %b %d %X %Y"))
  #write.table
  readed<-get_exact(readed)
  readed<-get_exact_rev(readed)
  for (i in 1:nrow(primers)){
    readed$primer[grepl(primers[i,2],substr(readed$read1,1,20))]<-primers[i,1]
    }
  print("Exact matches read1 found")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed_exact<-readed[readed$Exact!=""|readed$Exact_rev!="",]
  intersected_amp<-intersect_for_inv(readed_exact)
  #make_amplicones table here with assembled and ready to intersect amplicones. 
  amplicones<-NULL
  #for (i in 1:nrow(readed_exact))
  #  for (j in 1:nrow(readed_exact))
  #  {
  #    print(paste0(i,j,collapse="_"))
  #    if (is.null(intersected_amp[[i]][[j]]))
  #    amplicones<-rbind(amplicones, data.frame(assembled=intersected_amp[[i]][[j]],firstreads=readed_exact[i,]$readnumber,secondreads=readed_exact[j,]$readnumber,stringsAsFactors = F))
  #  }
  amplicones<-data.frame(assembled=unlist(intersected_amp),stringsAsFactors = F)
  resM<-character(nrow(amplicones))
  if(length(amplicones>0)){
    for (i in 1:length(resM))
      resM[i]<-paste0(HLA_base[grepl(gsub("N","[AGCT]",amplicones$assembled[i],fixed=T),x = HLA_base$Sequence),]$Allele,collapse = " ")
  }
  amplicones$Exact<-resM
  Amp<-intersect_amplicones(amplicones)
  print("Amplicones intersected")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  Amplist<-unlist(Amp)
  resM<-character(length(Amplist))
 # print(amplicones)
    if(length(Amplist>0)){
      for (i in 1:length(resM))
  resM[i]<-paste0(HLA_base[grepl(gsub("N","[AGCT]",Amplist[i],fixed=T),x = HLA_base$Sequence),]$Allele,collapse = " ")
   }
  print("Amplicones exact found")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  list(safety2=readed_exact,safety1=readed,answer=cbind(resM,Amplist))
  #list(AMP=Amp,safety4=table(resM),safety3=data.frame(Full_seq=Amplist,Exact_search=resM,stringsAsFactors = F),safety2=readed_exact,safety1=readed)
}
get_exact_rev<-function(merged){
  if (nrow(merged)==0)return(merged)
  
  resM<-character(nrow(merged))
  for (i in 1:length(resM))
    resM[i]<-paste0(HLA_base[grepl(revcomp(merged$read1[i]),x = HLA_base$Sequence,fixed=TRUE),]$Allele,collapse = " ")
  merged$Exact_rev<-resM
  merged
}

intersect_for_inv<-function(exact_list){
  res<-list()
  for (i in 1:nrow(exact_list)){res[[i]]<-list();#print(i)
  for (j in 1:nrow(exact_list))
  #  if (i<j)
      if (sum(duplicated(c(unlist(strsplit(as.character(exact_list[i,]$Exact)," ")),unlist(strsplit(as.character(exact_list[j,]$Exact_rev)," ")))))!=0)
      {
        res[[i]][[j]]<-make_overlaps(exact_list[i,]$read1,revcomp(exact_list[j,]$read1))
      }
  }
  res
}

get_primer_lists<-function(readed,forwards,reverses){
  tblpr_R1<-lapply(primers_all$V2,function(x){sum(grepl(x,substr(readed$read1,1,20)))})
  tblpr_R2<-lapply(primers_all$V2,function(x){sum(grepl(x,substr(readed$read2,1,20)))})
  primers_allc<-primers_all
  primers_allc$R1_pr<-unlist(tblpr_R1)
  primers_allc$R2_pr<-unlist(tblpr_R2)
  primers_allc
}

get_primer_lists_allele<-function(allele,forwards,reverses){
  tblpr_R1<-lapply(primers_all$V2,function(x){(regexpr(x,(allele)))})
  tblpr_R2<-lapply(primers_all$V2,function(x){(regexpr(x,revcomp(allele)))})
  primers_allc<-primers_all
  primers_allc$R1_pr<-unlist(tblpr_R1)
  primers_allc$R2_pr<-nchar(allele)-unlist(tblpr_R2)
  primers_allc
}

get_primer_lists_2D_MRD<-function(readed,primerlist){
  tblpr_R1<-lapply(primerlist$V2,function(x){(grepl(x,substr(readed$read1,1,30)))})
  tblpr_R2<-lapply(primerlist$V2,function(x){(grepl(x,substr(readed$read2,1,30)))})
  tblpr_R1<-do.call(cbind,tblpr_R1)
  tblpr_R2<-do.call(cbind,tblpr_R2)
  res<-t(tblpr_R1)%*%tblpr_R2
  rownames(res)<-primerlist$V1
  colnames(res)<-primerlist$V1
  list(read1=tblpr_R1,read2=tblpr_R2,combinations=res)
}

get_primer_lists_2D<-function(readed){
  tblpr_R1<-lapply(primers_all$V2,function(x){(grepl(x,substr(readed$read1,1,20)))})
  tblpr_R2<-lapply(primers_all$V2,function(x){(grepl(x,substr(readed$read2,1,20)))})
  tblpr_R1<-do.call(cbind,tblpr_R1)
  tblpr_R2<-do.call(cbind,tblpr_R2)
  t(tblpr_R1)%*%tblpr_R2
}

get_stats_mat<-function(prmat){
  tmp<-prmat[ampmat]
  tmp<-c(tmp,sum(prmat)-sum(prmat[ampmat]))
  names(tmp)<-c(row.names(ampmat),"Other")
  tmp
}

get_stats<-function(primers_mat_list){
  lapply(primers_mat_list,get_stats_mat)
}

File_list_pipeline_primer_stats<-function(filelist){
  resHLA<-lapply(filelist[,1],"[",1)
  names(resHLA)<-filelist[,1]
  for (i in 1:nrow(filelist)){
    print(names(resHLA)[i])
    print("File read start")
    print(format(Sys.time(), "%a %b %d %X %Y"))
    resHLA[[i]]<-get_primer_lists_2D(gzreader(filelist[i,2],filelist[i,3]))
  }
  resHLA
}

HLA_amplicones_full<-function(read1,read2,threshold=100,read_length=250){
  print("File read start")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  #readed<-reader(read1,read2)
  readed<-gzreader(read1,read2)
  #Iamp1log<-grepl("CCCTGACC[GC]AGACCTG",readed$read1)
  #inv<-grepl("[CT]GGTGG[AG]CTGGGAAGA|[CT]CAGCAGGTTGTGGTG|CCAG[GC]AGGTT[AG]TGGTG|CCAC[GT]TGGCAGGTGTA|CCACTTGGCAAGTGTA",substr(readed$read1,1,20))
  #readed[inv,c(1,2)]<-readed[inv,c(2,1)]
  readed<-list(Iamp1=readed[grepl("CCCTGACC[GC]AGACCTG",substr(readed$read1,1,20)),],
               Iamp2=readed[grepl("CGACGGCAA[AG]GATTAC",substr(readed$read1,1,20)),],
               IIamp1_DQB=readed[grepl("AG[GT]CTTTGCGGATCCC",substr(readed$read1,1,20)),],
               IIamp1_Others=readed[grepl("CTGAGCTCCC[GC]ACTGG",substr(readed$read1,1,20)),],
               IIamp2=readed[grepl("GGAACAGCCAGAAGGA",substr(readed$read1,1,20)),],
               Iamp1alt=readed[grepl("[CT]CACTCCATGAGGTATTTC|CCACTCCATGAAGTATTTC",substr(readed$read1,1,22)),],
               Iamp2alt=readed[grepl("GGCAA[AG]GATTACATCGCC|GGCAAGGATTACATCGCT",substr(readed$read1,1,22)),],
               IIamp1_DQBalt=readed[grepl("TGAGGGCAGAGAC[CT]CTCC",substr(readed$read1,1,22)),],
               IIamp1_DRBalt=readed[grepl("TGACAGTGACACTGATGG|TGACAGTGACATTGACGG",substr(readed$read1,1,22)),],
               IIamp2_DRBalt=readed[grepl("GGTTTCTATCCAGGCAGC",substr(readed$read1,1,22))&grepl("TG[CT]TCTGGGCAGATTCAG",substr(readed$read2,1,22)),],
               IIamp2_DQBalt=readed[grepl("ACCATCTCCCCATCCAG",substr(readed$read1,1,22))&grepl("TG[CT]TCTGGGCAGATTCAG",substr(readed$read2,1,22)),],
               Iamp1_inv=readed[grepl("GGGCCGCCTCC[AC]ACTTG|GGGCCGTCTCCCACTTG|GGACCGCCTCCCACTTG",substr(readed$read1,1,20)),],
               Iamp2_inv=readed[grepl("[CT]GGTGG[AG]CTGGGAAGA",substr(readed$read1,1,20)),],
               IIamp1_DQB_inv=readed[grepl("[CT]CAGCAGGTTGTGGTG|CCAG[GC]AGGTT[AG]TGGTG",substr(readed$read1,1,20))&grepl("AG[GT]CTTTGCGGATCCC",substr(readed$read2,1,20)),],
               IIamp1_Others_inv=readed[grepl("[CT]CAGCAGGTTGTGGTG|CCAG[GC]AGGTT[AG]TGGTG",substr(readed$read1,1,20))&grepl("CTGAGCTCCC[GC]ACTGG",substr(readed$read2,1,20)),],
               IIamp2_inv=readed[grepl("CCAC[GT]TGGCAGGTGTA|CCACTTGGCAAGTGTA",substr(readed$read1,1,20)),],
               Iamp1alt_inv=readed[grepl("GAGC[GC]ACTCCACGCAC|GAGCCCGTCCACGCAC",substr(readed$read1,1,22)),],
               Iamp2alt_inv=readed[grepl("TCAGGGTGAGGGGCT[CT]|TCAGGGTGCAGGGCTC",substr(readed$read1,1,22)),],
               IIamp1_DQBalt_inv=readed[grepl("GTCCAGTCACC[AG]TTCCTA",substr(readed$read1,1,22)),],
               IIamp1_DRBalt_inv=readed[grepl("CAG[CT]CTTCTCTTCCTGGC",substr(readed$read1,1,22)),],
               IIamp2_DRBalt_inv=readed[grepl("GGTTTCTATCCAGGCAGC",substr(readed$read2,1,22))&grepl("TG[CT]TCTGGGCAGATTCAG",substr(readed$read1,1,22)),],
               IIamp2_DQBalt_inv=readed[grepl("ACCATCTCCCCATCCAG",substr(readed$read2,1,22))&grepl("TG[CT]TCTGGGCAGATTCAG",substr(readed$read1,1,22)),]
               )
  print(sapply(readed,nrow))
  print("File read end. Looking for champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-lapply(lapply(readed,get_groups,threshold = threshold,read_length = read_length, primer_length=20),merge_cons)
  print("We are the champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  print(sapply(readed,nrow))
  readed<-lapply(readed,get_overlap_merged)
  readed$Iamp1alt<-get_overlap_merged_fix(readed$Iamp1alt)
  readed$Iamp1alt_inv<-get_overlap_merged_fix(readed$Iamp1alt_inv)
  print("Reads_assembled")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-lapply(readed,get_exact)
  readed<-lapply(readed,get_exact_rev)
  print("Exact matches read1 found")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-lapply(readed,fathers_and_children)
  readed[12:22]<-lapply(readed[12:22],straight_inverse)
 # readed<-amplist
  readed_intersect<-lapply(readed,function(x){x[(x$parents==1|(x$parents==2&x$freq>0.005)),]})
  print("Graph made. Ready for papas intersection")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  print(sapply(readed_intersect,nrow))
  res<-list(safety1=readed,
            safety2=readed_intersect,
            Iclass=intersect_amplicones3(rbind(readed_intersect$Iamp1,readed_intersect$Iamp1_inv),rbind(readed_intersect$Iamp2,readed_intersect$Iamp2_inv)),
            IIclassDQB=intersect_amplicones3(rbind(readed_intersect$IIamp1_DQB,readed_intersect$IIamp1_DQB_inv),rbind(readed_intersect$IIamp2,readed_intersect$IIamp2_inv)),
            IIclassOthers=intersect_amplicones3(rbind(readed_intersect$IIamp1_Others,readed_intersect$IIamp1_Others_inv),rbind(readed_intersect$IIamp2,readed_intersect$IIamp2_inv)),
            #Iclass_inv=intersect_amplicones3(readed_intersect$Iamp1_inv,readed_intersect$Iamp2_inv),
            #IIclassDQB_inv=intersect_amplicones3(readed_intersect$IIamp1_DQB_inv,readed_intersect$IIamp2_inv),
            #IIclassOthers_inv=intersect_amplicones3(readed_intersect$IIamp1_Others_inv,readed_intersect$IIamp2_inv),
            
            Iclass_alt=intersect_amplicones3(rbind(readed_intersect$Iamp1alt,readed_intersect$Iamp1alt_inv),rbind(readed_intersect$Iamp2alt,readed_intersect$Iamp2alt_inv)),
            IIclassDQB_alt=intersect_amplicones3(rbind(readed_intersect$IIamp1_DQBalt,readed_intersect$IIamp1_DQBalt_inv),rbind(readed_intersect$IIamp2_DQBalt,readed_intersect$IIamp2_DQBalt_inv)),
            IIclassOthers_alt=intersect_amplicones3(rbind(readed_intersect$IIamp1_DRBalt,readed_intersect$IIamp1_DRBalt_inv),rbind(readed_intersect$IIamp2_DRBalt,readed_intersect$IIamp2_DRBalt_inv))
            # Iclass_inv_alt=intersect_amplicones3(readed_intersect$Iamp1alt_inv,readed_intersect$Iamp2alt_inv),
            # IIclassDQB_inv_alt=intersect_amplicones3(readed_intersect$IIamp1_DQBalt_inv,readed_intersect$IIamp2_DQBalt_inv),
            #IIclassOthers_inv_alt=intersect_amplicones3(readed_intersect$IIamp1_DRBalt_inv,readed_intersect$IIamp2_DRBalt_inv)
  )
  print("Amplicones intersected. Report generated.")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  res
}  
  

File_list_pipeline_amplicones<-function(filelist,read_length=250,threshold=100){
  resHLA<-lapply(filelist[,1],"[",1)
  names(resHLA)<-filelist[,1]
  for (i in 1:nrow(filelist)){
    print(names(resHLA)[i])
    resHLA[[i]]<-HLA_amplicones_full(filelist[i,2],filelist[i,3],read_length = read_length,threshold=threshold)
  }
  resHLA
}

straight_inverse<-function(readed){
  if (nrow(readed)==0)return(readed)
  readed$Exact<-readed$Exact_rev
  readed$assembled<-revcomp(readed$assembled)
  readed
}  

add_grep<-function(readed,length=50){
  readed$grepleft<-substr(readed$assembled,1,length)  
  readed$grepright<-substr(readed$assembled,nchar(readed$assembled)-length+1,nchar(readed$assembled))
  readed$grepleft<-gsub("N","[AGCT]",readed$grepleft,fixed=T)
  readed$grepright<-gsub("N","[AGCT]",readed$grepright,fixed=T)
}

HLA_amplicones_intersect<-function(amplist){
#add grep to each amplicone(right for right, left fro left)  
#use amplicone intersect matrix(with names?)
  amplist[12:22]<-lapply(amplist[12:22],straight_inverse)
  readed<-amplist
  readed_intersect<-lapply(readed,function(x){x[(x$parents==1|x$parents==2)&(x$freq>0.01),]})
  print("Amplicones intersected. Report generated.")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  print(sapply(readed_intersect,nrow))
  res<-list(safety1=readed,
            safety2=readed_intersect,
            Iclass=intersect_amplicones3(readed_intersect$Iamp1,readed_intersect$Iamp2),
            IIclassDQB=intersect_amplicones3(readed_intersect$IIamp1_DQB,readed_intersect$IIamp2),
            IIclassOthers=intersect_amplicones3(readed_intersect$IIamp1_Others,readed_intersect$IIamp2),
            Iclass_inv=intersect_amplicones3(readed_intersect$Iamp1_inv,readed_intersect$Iamp2_inv),
            IIclassDQB_inv=intersect_amplicones3(readed_intersect$IIamp1_DQB_inv,readed_intersect$IIamp2_inv),
            IIclassOthers_inv=intersect_amplicones3(readed_intersect$IIamp1_Others_inv,readed_intersect$IIamp2_inv),
  
            Iclass_alt=intersect_amplicones3(readed_intersect$Iamp1alt,readed_intersect$Iamp2alt),
            IIclassDQB_alt=intersect_amplicones3(readed_intersect$IIamp1_DQBalt,readed_intersect$IIamp2_DQBalt),
            IIclassOthers_alt=intersect_amplicones3(readed_intersect$IIamp1_DRBalt,readed_intersect$IIamp2_DRBalt),
            Iclass_inv_alt=intersect_amplicones3(readed_intersect$Iamp1alt_inv,readed_intersect$Iamp2alt_inv),
            IIclassDQB_inv_alt=intersect_amplicones3(readed_intersect$IIamp1_DQBalt_inv,readed_intersect$IIamp2_DQBalt_inv),
            IIclassOthers_inv_alt=intersect_amplicones3(readed_intersect$IIamp1_DRBalt_inv,readed_intersect$IIamp2_DRBalt_inv)
            )
  print("Amplicones intersected. Report generated.")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  res
  
}

HLA_amplicones_intersect_all<-function(amplist){
  #add grep to each amplicone(right for right, left fro left)  
  #use amplicone intersect matrix(with names?)
  amplist[12:22]<-lapply(amplist[12:22],straight_inverse)
  readed<-amplist
  readed_intersect<-lapply(readed,function(x){x[(x$parents==1|(x$parents==2&x$freq>0.005)),]})
  print("Amplicones intersected. Report generated.")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  print(sapply(readed_intersect,nrow))
  res<-list(safety1=readed,
            safety2=readed_intersect,
            Iclass=intersect_amplicones3(rbind(readed_intersect$Iamp1,readed_intersect$Iamp1_inv),rbind(readed_intersect$Iamp2,readed_intersect$Iamp2_inv)),
            IIclassDQB=intersect_amplicones3(rbind(readed_intersect$IIamp1_DQB,readed_intersect$IIamp1_DQB_inv),rbind(readed_intersect$IIamp2,readed_intersect$IIamp2_inv)),
            IIclassOthers=intersect_amplicones3(rbind(readed_intersect$IIamp1_Others,readed_intersect$IIamp1_Others_inv),rbind(readed_intersect$IIamp2,readed_intersect$IIamp2_inv)),
            #Iclass_inv=intersect_amplicones3(readed_intersect$Iamp1_inv,readed_intersect$Iamp2_inv),
            #IIclassDQB_inv=intersect_amplicones3(readed_intersect$IIamp1_DQB_inv,readed_intersect$IIamp2_inv),
            #IIclassOthers_inv=intersect_amplicones3(readed_intersect$IIamp1_Others_inv,readed_intersect$IIamp2_inv),
            
            Iclass_alt=intersect_amplicones3(rbind(readed_intersect$Iamp1alt,readed_intersect$Iamp1alt_inv),rbind(readed_intersect$Iamp2alt,readed_intersect$Iamp2alt_inv)),
            IIclassDQB_alt=intersect_amplicones3(rbind(readed_intersect$IIamp1_DQBalt,readed_intersect$IIamp1_DQBalt_inv),rbind(readed_intersect$IIamp2_DQBalt,readed_intersect$IIamp2_DQBalt_inv)),
            IIclassOthers_alt=intersect_amplicones3(rbind(readed_intersect$IIamp1_DRBalt,readed_intersect$IIamp1_DRBalt_inv),rbind(readed_intersect$IIamp2_DRBalt,readed_intersect$IIamp2_DRBalt_inv))
           # Iclass_inv_alt=intersect_amplicones3(readed_intersect$Iamp1alt_inv,readed_intersect$Iamp2alt_inv),
           # IIclassDQB_inv_alt=intersect_amplicones3(readed_intersect$IIamp1_DQBalt_inv,readed_intersect$IIamp2_DQBalt_inv),
            #IIclassOthers_inv_alt=intersect_amplicones3(readed_intersect$IIamp1_DRBalt_inv,readed_intersect$IIamp2_DRBalt_inv)
  )
  print("Amplicones intersected. Report generated.")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  res
  
}
Aggregate<-function(amplist){
  amplist<-lapply(names(amplist[-c(1,2)]),function(x){amplist[[x]]$Name=x;amplist[[x]]})
  tmp<-as.data.table(do.call(rbind,amplist[(sapply(amplist,function(x)!is.null(ncol(x))))]))
 
  if (nrow(tmp)!=0)
  tmp[,list(coverage=sum(geom_min_coverage),Unique_names=paste0(unique(Name),collapse=" "),Lengths=paste0(unique(Length),collapse=" "),two_signs=two_symbols(Exact_match)),Exact_match][order(Exact_match)][-1,,]
}
Aggregate_list<-function(amplist_list){
  tmp<-names(amplist_list)
  amplist_list<-lapply(names(amplist_list),function(x){Aggregate(amplist_list[[x]]);})
  names(amplist_list)<-tmp
  amplist_list<-lapply(names(amplist_list),function(x){if (is.data.table(amplist_list[[x]]))amplist_list[[x]][,Donor:=x,];amplist_list[[x]]})
  do.call(rbind,amplist_list)
}
two_symbols<-function(x){
  if (length(x)>1)sapply(x,two_symbols)else{
  allelelist<-unlist(strsplit(x," "))
  paste0(unique(sapply(allelelist,function(x){paste0(unlist(strsplit(x,":"))[c(1,2)],collapse=":")})),collapse=" ")}
}
remove_duplicates<-function(tbl){tbl[!duplicated(tbl),]}
write_safety1<-function(tbl,fname){
  if (nrow(tbl)!=0)
  write.table(tbl[order(-tbl$readnumber),],sep=";",quote = F,row.names = F,file = fname)
}

write_full<-function(tbl,fname){
  if (!is.null(nrow(tbl)))
  {
  tbl<-remove_duplicates(tbl)
  write.table(tbl[order(tbl$Overlap_mismatches,-tbl$geom_min_coverage),],sep=";",quote = F,row.names = F,file = fname)}
}

deserialise_HLA<-function(amplist,prefix)
{
for (i in 1:length(amplist$safety1))
  write_safety1(amplist$safety1[[i]],paste0(prefix,"_",names(amplist$safety1)[i],"_","safety1.csv",collapse=""))

for (i in 1:length(amplist[-c(1,2)]))
  write_full(amplist[-c(1,2)][[i]],paste0(prefix,"_",names(amplist[-c(1,2)])[i],"_","full.csv",collapse=""))
}

deserialise_HLA_list<-function(HLA_list){
  for (i in 1:length(HLA_list))
    deserialise_HLA(HLA_list[[i]],names(HLA_list[i]))
}


