gene_to_contig=function(x,y){
  for (i in 1:nrow(x)){
  tmp=strsplit(x[i,y],split = '_')
  tmp=as.data.frame(do.call(rbind,tmp))
  tmp=paste(tmp[,1:ncol(tmp)-1],collapse = '_')
  x$contigs[i]=tmp
  }
  return(x$contigs)
}
match=read.csv("Total_GFVD_jumbophage_95-80_drep.ffn.fna.mmseqs.out_filtering",header = F,sep = '\t')
geographic=read.table("vOTUs_geo.tsv",sep = '\t',header = 1,row.names = 1)

match$contigs_q=gene_to_contig(match,1)
match$contigs_t=gene_to_contig(match,2)
match$geo_q=geographic[as.character(match$contigs_q),1]
match$geo_t=geographic[as.character(match$contigs_t),1]
library(stringr)
for (i in 1:nrow(match)){
  if(str_detect(match$geo_q[i],"later")==TRUE && str_detect(match$geo_t[i],"later")==TRUE){
    match$group[i]="intra_late"
  }
  else if(str_detect(match$geo_q[i],"early")==TRUE && str_detect(match$geo_t[i],"early")==TRUE){
    print("intra_early")
    match$group[i]="intra_early"
  }
  else if(str_detect(match$geo_q[i],"early")==TRUE && str_detect(match$geo_t[i],"later")==TRUE){
    match$group[i]="inter_stage"
  }
  else if(str_detect(match$geo_q[i],"later")==TRUE && str_detect(match$geo_t[i],"early")==TRUE){
    match$group[i]="inter_stage"
  }
  else if(str_detect(match$geo_q[i],"inter")==TRUE || str_detect(match$geo_t[i],"inter")==TRUE){
    match$group[i]="inter_stage"
  }
}
table(match$group)
nrow(unique(rbind(match[grep("late",match$group),1],match[grep("late",match$group),2])))
tmp=rbind(as.matrix(match[grep("early",match$group),1]),as.matrix(match[grep("early",match$group),2]))
tmp=as.data.frame(tmp)
percent=data.frame(group=c("inter_stage","intra_late","intra_early"),percent=c(3,1.09,0.6))
percent$group=factor(percent$group,levels = c("inter_stage","intra_late","intra_early"))
library(ggplot2)
library(ggsci)
ggplot(percent,aes(group,percent,fill=group))+
  geom_bar(stat = 'identity',width = 0.8)+theme_bw()+
  scale_fill_aaas()
  


tmp=tmp[!duplicated(tmp$V1),]
tmp=rbind(as.matrix(match$contigs_q),as.matrix(match$contigs_t))%>%as.data.frame()
tmp=tmp[!duplicated(tmp$V1),]%>%as.data.frame()


network=match[,c(15,16)]
network$value=1
library(reshape2)
network=dcast(data = network,formula = contigs_q~contigs_t)
network=melt(network)
network=network[network$value>0,]
write.table(x = network,file = "lateral_gene_transfer.net",sep = '\t',row.names = F,quote = F)
