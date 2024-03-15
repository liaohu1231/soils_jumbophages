library(picante)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(dplyr)
library(ggpubr)
instrain=read.table("all_PJPD_scaffold_nucl.tsv",header = 1,row.names = 1,check.names = F)
instrain_shannon=as.data.frame(diversity(t(instrain), "shannon"))
data=read.table("../Macrodiversity/all.PJPD.rpkm",header = 1,row.names = 1,check.names = F)
relative_abundance=function(d){
  dta_sum=apply(d,2,function(x){x/sum(x)})
}
data_rel=relative_abundance(data)
abu_shannon=as.data.frame(diversity(t(data_rel), "shannon"))
colnames(abu_shannon)="abu_shannon"
colnames(instrain_shannon)="instrain_shannon"

instrain_shannon$abu_shannon=abu_shannon[rownames(instrain_shannon),1]

summary(lm(data = instrain_shannon,formula = abu_shannon~instrain_shannon))

ggplot(instrain_shannon,aes(instrain_shannon,abu_shannon))+
  #geom_boxplot()+
  geom_point(alpha=0.5,size=4)+
  theme_classic()+
  labs(x="shannon index of Microdiversity",y="shannon index of Macrodiversity")+
  geom_smooth(aes(instrain_shannon,abu_shannon),method='lm',
              formula=y~x,level=0.95,color='blue')

ggsave("shannon_index_microdiversity.pdf",device = 'pdf',width = 6,height = 5)

##instrain_melt
library(reshape2)
instrain_melt=melt(as.matrix(instrain))
instrain_melt=instrain_melt[instrain_melt$value>0,]
instrain_melt=na.omit(instrain_melt)
length=read.table("quality_summary_putative_jumbo.tsv",sep = '\t',header = F,row.names = 1,check.names = F)
colnames(length)='length'
group=read.csv("../Macrodiversity//metadata.txt",sep = '\t')
group$stage=rep(c("early_stage","later_stage"),c(12,12))
Freq=data.frame(Freq=apply(data_rel,1,function(x){sum(x>0)}))
instrain_melt$Freq=as.factor(Freq[as.character(instrain_melt$Var1),1])
instrain_melt=instrain_melt[instrain_melt$Freq!='0',]
instrain_melt$length=length[as.character(instrain_melt$Var1),1]

for (i in 1:nrow(data_rel)){
  for(j in 1:ncol(data_rel)){
    if (data_rel[i,j]>0){
      data_rel[i,j]=1
    }
  }
}
data_rel=as.data.frame(t(data_rel))
rownames(group)=group$sample
data_rel$site=group[rownames(data_rel),2]
data_rel$stage=group[rownames(data_rel),4]

data_sites=as.data.frame(matrix(NA,8,nrow(data)))
colnames(data_sites)=rownames(data)
for (i in 1:nrow(data)) {
  print(i)
  data_sites[,i]=tapply(data_rel[,i],INDEX = data_rel$site,sum)
  for(j in 1:nrow(data_sites)){
    if(data_sites[j,i]>0){
      data_sites[j,i]=1
    }
  }
}
tmp=as.data.frame(tapply(data_rel[,1],INDEX = data_rel$site,sum))
rownames(data_sites)=rownames(tmp)
data_stages=as.data.frame(matrix(NA,2,nrow(data)))
colnames(data_stages)=rownames(data)

for (i in 1:nrow(data)){
  print(i)
  data_stages[,i]=tapply(data_rel[,i],INDEX = data_rel$stage,sum)
  for(j in 1:nrow(data_stages)){
    if(data_stages[j,i]>0){
      data_stages[j,i]=1
    }
  }
}
tmp=as.data.frame(tapply(data_rel[,1],INDEX = data_rel$stage,sum))
rownames(data_stages)=rownames(tmp)

vOTUs_geo=data.frame(geogra_range=rep(c(NA),c(nrow(data))))
rownames(vOTUs_geo)=rownames(data)
for (i in 1:nrow(data)){
  if (sum(data_rel[,i])==1 && data_stages[1,i]==1){
    vOTUs_geo[i,]='sample_specific (early)'
  }
  else if(sum(data_rel[,i])==1 && data_stages[2,i]==1){
    vOTUs_geo[i,]='sample_specific (later)'
  }else{
    if(sum(data_sites[,i])==1 && data_stages[1,i]==1){
      vOTUs_geo[i,]='site_specific (early)'
    }
    else if(sum(data_sites[,i])==1 && data_stages[2,i]==1){
      vOTUs_geo[i,]='site_specific (later)'
    }else{
    if(sum(data_stages[,i])==2){
      vOTUs_geo[i,]="inter_stage"
    }else{
      if(data_stages[1,i]==1){
        vOTUs_geo[i,]='intra_early_stage'
      }
      if(data_stages[2,i]==1){
        vOTUs_geo[i,]='intra_later_stage'
      }
    }
  }
  }
}
data_sites_freq=data.frame(site_freq=apply(data_sites,2,sum))

instrain_melt$geographic_range=vOTUs_geo[as.character(instrain_melt$Var1),1]
write.table(vOTUs_geo,"puative_vOTUs_geo.tsv",sep = '\t',row.names = T,quote = F)

ggplot(instrain_melt,aes(geographic_range,log10(value)))+
  geom_boxplot(aes(fill=geographic_range))+
  #geom_jitter(height = 0,width = 0.3,alpha=0.5,size=4,shape=21)+
  theme_classic()+
  labs(x="Geographic range of jumbo vOTUs",y="Log10(Microdiversity (pi))")+
  scale_fill_d3()+
  stat_compare_means(comparisons = list(c("intra_later_stage","intra_early_stage"),
                                        c("intra_later_stage","inter_stage"),
                                        c("sample_specific (early)","sample_specific (later)"),
                                        c('site_specific (early)','site_specific (later)'),
                                        c("intra_later_stage",'sample_specific (early)')))+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,size = 10,face = "bold"),
        axis.text.y = element_text(size = 13,face = "bold"),
        axis.title = element_text(size=18))

ggsave("Putative_geographic_range_microdiveristy.pdf",device = 'pdf',width = 10,height = 5)
###inter_stage
instrain_inter=instrain_melt[grep("inter",instrain_melt$geographic_range),]
instrain_inter$stage=group[as.character(instrain_inter$Var2),4]
ggplot(instrain_inter,aes(x=stage,y=log10(value),group=stage))+
  #geom_density(mapping = aes(fill=stage,color=stage),alpha=0.6)+
  geom_jitter(height = 0,width = 0.28,alpha=0.6,size=2.7,shape=21,fill="#F2BD14")+
  geom_boxplot(fill="#90AFDA",alpha=0.9)+
  theme_classic()+
  labs(y="Log10(Microdiversity (pi))")+
  scale_color_aaas()+
  stat_compare_means(comparisons = list(c("early_stage","later_stage")),
                     method = "t.test")+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,size = 10,face = "bold"),
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title = element_text(size=15))
shapiro.test(x = instrain_inter[grep("later",instrain_inter$stage),3])
shapiro.test(x = instrain_inter[grep("early",instrain_inter$stage),3])
wilcox.test(instrain_inter[grep("later",instrain_inter$stage),3],
       instrain_inter[grep("early",instrain_inter$stage),3],alternative = "greater")
median(instrain_melt$value)
ggsave("putative_inter_stage_boxplot.pdf",device = 'pdf',width = 3,height = 5)

for(i in 1:nrow(instrain_melt)){
  instrain_melt$relative_abundance[i]=data_rel[which(rownames(data_rel)==instrain_melt$Var1[i]),
                                            grep(instrain_melt$Var2[i],colnames(data_rel))]
}
summary(lm(formula = value~relative_abundance,data = instrain_melt))
cor(x = instrain_melt$relative_abundance,y = instrain_melt$value,method = 'pearson')
ggplot(instrain_melt,aes(relative_abundance,log10(value)))+
  geom_point()+theme_bw()+
  geom_smooth(formula = y~x,method = 'lm')+
  labs(y="log10(Microdiversity (pi))")
ggsave("relative_abundance_microdiveristy.pdf",device = 'pdf',width = 6,heigh=5)

#host
host=read.csv("../tree/all_GFJPD_host.csv",header = F,row.names = 1)
table_host=as.data.frame(table(host$V4))
table_host=table_host[order(table_host$Freq,decreasing = T),]

for (i in 1:nrow(instrain_melt)){
  tmp=strsplit(as.character(instrain_melt$Var1)[i],split = '.fna')
  tmp=as.data.frame(do.call(rbind,tmp))
  instrain_melt$candidate_contig[i]=tmp
}
instrain_melt$host=host[as.character(instrain_melt$candidate_contig),3]

instrain_host=rbind(instrain_melt[which(as.character(instrain_melt$host)=="Firmicutes"),],
                    instrain_melt[which(as.character(instrain_melt$host)=="Proteobacteria"),],
                    instrain_melt[which(as.character(instrain_melt$host)=="Bacteroidetes"),],
                    instrain_melt[which(as.character(instrain_melt$host)=="Actinobacteria"),])


ggplot(instrain_host,aes(host,log10(value),fill=stages))+
  geom_boxplot()+
  #geom_jitter(height = 0,width = 0.3,alpha=0.5,size=4,shape=21)+
  labs(y="Microdiversity (pi)")+
  theme_classic()+
  scale_fill_aaas()+
  stat_compare_means(comparisons = list(c("Firmicutes","Proteobacteria"),
                                           c("Bacteroidetes","Proteobacteria"),
                                          c("Actinobacteria","Proteobacteria"),
                                          c("Actinobacteria","Firmicutes"),
                                          c("Actinobacteria","Bacteroidetes")),
                                     method = "wilcox.test")+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))
wilcox.test(instrain_host[grep(Actinobacteria),])
summary(lm(data = instrain_melt,formula = log(value)~length))
cor.test(x = instrain_melt$value,y = instrain_melt$length/1000,method="spearman")
ggplot(instrain_melt,aes(length,log10(value)))+
  #geom_boxplot()+
  geom_point(alpha=0.5,size=2)+
  theme_classic()+
  labs(x="length",y="log(Microdiversity (pi))")+
  geom_smooth(aes(as.numeric(length),log10(value)),method='lm',
              formula=y~x,level=0.95,color='blue')+
  scale_fill_d3(palette = "category10")

instrain_melt$length=length[instrain_melt$Var1,1]
instrain_melt$stages=group[instrain_melt$Var2,4]
instrain_melt$quality=length[instrain_melt$Var1,7]
instrain_melt_length=instrain_melt[-grep("Medium-quality",instrain_melt$quality),]
instrain_melt_length=instrain_melt_length[-grep("Low-quality",instrain_melt_length$quality),]

for (i in 1:nrow(instrain_melt)){
  if(instrain_melt$length[i] > 330000){
    instrain_melt$phages[i] = ">330kbp"
  }else{
    instrain_melt$phages[i] = "<330kbp"
  }
}
library(ggpointdensity)
ggplot(instrain_melt,aes(phages,log10(value),fill=stages))+
  #geom_point(alpha=0.8,size=2)+
  geom_boxplot()+
  #geom_pointdensity() +
  theme_classic()+
  labs(x="length",y="log(Microdiversity (pi))")+
  #geom_smooth(aes(as.numeric(length),log10(value)),method='glm',
   #           formula=y~x,level=0.95,color='blue')+
  scale_fill_d3(palette = "category10")+
  stat_compare_means(comparisons = list(c(">330kbp","<330kbp")))
summary(glm(data = instrain_melt,
           formula = log(value)~length))


###gene_microdiversity——pNpS
pnps=read.csv("all.gene.pNpS.tsv",sep='\t',row.names = 1,header = 1,check.names = F)
for (i in 1:nrow(pnps)){
  tmp=strsplit(rownames(pnps)[i],split = '_')
  tmp=as.data.frame(do.call(rbind,tmp))
  tmp=paste(tmp[,1:ncol(tmp)-1],collapse = '_')
  pnps$contigs[i]=tmp
}

library(reshape2)
pnps$contigs=as.factor(pnps$contigs)
pnps$genes=rownames(pnps)

pnps$geo=vOTUs_geo[as.character(pnps$contigs),1]
pnps_geo=as.data.frame(matrix(data = NA,nrow = 7,ncol = 1))
rownames(pnps_geo)=unique(vOTUs_geo$geogra_range)[1:7]
pnps_geo=pnps_geo[,-1]
source("new.cbind.R")
for (i in 1:24) {
  tmp=pnps[,c(i,27)]
  tmp=na.omit(tmp)
  tmp=as.matrix(tapply(X = tmp[,1],tmp$geo,mean))
  colnames(tmp)=colnames(pnps)[i]
  tmp=as.data.frame(tmp)
  pnps_geo=new.cbind(pnps_geo,tmp)
}

pnps_melt=melt(as.matrix(pnps_geo))
pnps_melt=na.omit(pnps_melt)
#pnps_melt$geographic_range=vOTUs_geo[as.character(pnps_melt$contigs),1]
#pnps_melt=na.omit(pnps_melt)

ggplot(pnps_melt,aes(Var1,value))+
  geom_violin(aes(color=Var1),fill="gray90")+
  geom_boxplot(aes(fill=Var1),width=0.5)+
  #geom_jitter(height = 0,width = 0.3,alpha=0.5,size=4,shape=21)+
  theme_classic()+
  labs(x="Geographic range of jumbo vOTUs",y="log(pNpS)")+
  scale_fill_d3()+
  stat_compare_means(comparisons = list(c("intra_later_stage","intra_early_stage"),
                                        c("intra_later_stage","inter_stage"),
                                        c("sample_specific (later)","site_specific (later)"),
                                        c("intra_later_stage","site_specific (later)")))+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,size = 10,face = "bold"),
        axis.text.y = element_text(size = 13,face = "bold"),
        axis.title = element_text(size=18))

###pnps_host
for (i in 1:nrow(pnps_melt)){
  tmp=strsplit(as.character(pnps_melt$contigs)[i],split = '.fna')
  tmp=as.data.frame(do.call(rbind,tmp))
  pnps_melt$candidate_contig[i]=tmp
}
pnps_melt$host=host[as.character(pnps_melt$candidate_contig),3]

pnps_host=rbind(pnps_melt[which(as.character(pnps_melt$host)=="Firmicutes"),],
                    pnps_melt[which(as.character(pnps_melt$host)=="Proteobacteria"),],
                    pnps_melt[which(as.character(pnps_melt$host)=="Bacteroidetes"),],
                    pnps_melt[which(as.character(pnps_melt$host)=="Actinobacteria"),])

ggplot(pnps_host,aes(host,log10(value),fill=host))+
  geom_boxplot()+
  #geom_jitter(height = 0,width = 0.3,alpha=0.5,size=4,shape=21)+
  labs(y="log10(pNpS)")+
  theme_classic()+
  scale_fill_aaas()+
  stat_compare_means(comparisons = list(c("Firmicutes","Proteobacteria"),
                                        c("Bacteroidetes","Proteobacteria"),
                                        c("Actinobacteria","Proteobacteria"),
                                        c("Actinobacteria","Firmicutes"),
                                        c("Actinobacteria","Bacteroidetes")),
                     method = "wilcox.test")+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))

kegg=read.csv("../huge_phage_merge_derep_95_80_g16.bestexec_anno2",header = F,row.names = 1,sep='\t')
pNpS$ko=as.factor(kegg[pNpS$contig_gene,1])
pNpS_positive=pNpS[pNpS$pNpS_ratio>1,]
nrow(pNpS_positive[!duplicated(pNpS_positive$contig_gene),])
pNpS_positive=na.omit(pNpS_positive)
nrow(pNpS_positive[!duplicated(pNpS_positive$contig_gene),])

nd=matrix(NA,nrow(pNpS_positive),1)
for(i in 1:nrow(nd)){
  name=pNpS_positive$ko[i]
  map=kegg[grep(name,x = kegg$V2),5]%>%as.data.frame()
  nd[i,]=map[1,1]
}

pNpS_positive$enzyme=nd
pNpS_positive$flag=paste(pNpS_positive$enzyme,'(',pNpS_positive$ko,')')
pNpS_positive_1=pNpS_positive[!duplicated(pNpS_positive$contig_gene),]
tmp=table(as.character(pNpS_positive_1$flag))%>%as.data.frame()

library(ggplot2)
ggplot(pNpS_positive,aes(x = ko,y=pNpS_ratio))+
  geom_boxplot()+
  coord_flip()+theme_classic()

ggplot(pnps_melt,mapping =aes(x=log10(value)))+
  geom_density(fill="lightblue",alpha=0.7)+
  theme_classic()+
  geom_vline(xintercept = -0.6828)
10^-0.6828
write.table(pNpS_positive,"pNpS_positive.txt",sep='\t',quote = F,row.names = F)

##positive_genes_annotation
positive_genes=pnps_melt[pnps_melt$value>0.2,]
positive_genes_drep=positive_genes[!duplicated(positive_genes$genes),]
annotation=read.table("../functions/Total_GFVD_jumbophage_95-80_drep.best_kofam_annotation",row.names = 1,header = F,sep='\t')
positive_genes_drep$KEGG=annotation[as.character(positive_genes_drep$genes),1]
positive_genes_drep$annotaion=annotation[as.character(positive_genes_drep$genes),5]
anno_table=as.data.frame(table(positive_genes_drep$annotaion))
positive_genes_drep=as.data.frame(positive_genes_drep)

write.table(x = positive_genes_drep[,c(1:5,7:9)],file = "Table_S6_positive_genes.tsv",quote = F,
          row.names = F,sep='\t',na = "NA")

sum(anno_table$Freq)
###conANI
instrain_compare=read.table("compare_comparisonsTable.tsv",sep='\t',
                            header = 1)
instrain_compare$group="NA"
for (i in 1:nrow(instrain_compare)){
  stage_2=group[grep(instrain_compare[i,2],group$sample),4]
  stage_3=group[grep(instrain_compare[i,3],group$sample),4]
  sites_2=group[grep(instrain_compare[i,2],group$sample),2]
  sites_3=group[grep(instrain_compare[i,3],group$sample),2]
  if(stage_2 == "early_stage" && stage_3 == "early_stage")
  {
    if(sites_2==sites_3){
      instrain_compare$group[i]="intra_sites (early)"
    }else{
      instrain_compare$group[i] = "intra_early"
    }
  }
  else if(stage_2 == "late_stage" && stage_3 == "late_stage")
  {
    if(sites_2==sites_3){
    instrain_compare$group[i]="intra_sites (late)"
    }else{
      instrain_compare$group[i] = "intra_late"
      }
  }
  else{
    instrain_compare$group[i]= "inter_stage"
  }
}
instrain_compare$group=factor(instrain_compare$group,
                              levels = c("intra_sites (early)","intra_sites (later)",
                                         "intra_early","intra_later",
                                         "inter_stage"))
library(ggplot2)
library(ggpubr)
library(ggsci)
instrain_compare1=instrain_compare[instrain_compare$coverage_overlap>0.5,]
p=ggplot(instrain_compare1,aes(group,conANI,fill=group))+
  geom_boxplot()+theme_bw()+
  stat_compare_means(
                     comparisons = list(c("inter_stage","intra_later"),
                                        c("intra_early","intra_later"),
                                        c("intra_early","inter_stage"),
                                        c("intra_sites (early)","intra_later")))+
  theme(axis.text.x = element_text(angle = 40,hjust = 1))+
  scale_fill_lancet()

p
ggsave("conANI_boxplot.pdf",plot = p,device = "pdf",width = 5,height = 4)
instrain_compare1$year_distance=NA
for (i in 1:nrow(instrain_compare1)){
  stage_2=group[grep(instrain_compare1[i,2],group$sample),3]
  stage_3=group[grep(instrain_compare1[i,3],group$sample),3]
  instrain_compare1$year_distance[i]=abs(as.numeric(stage_2) - as.numeric(stage_3))
}
intra_later=instrain_compare1[grep("intra_later",instrain_compare1$group),]
intra_early=instrain_compare1[grep("intra_early",instrain_compare1$group),]
inter_stage=instrain_compare1[grep("inter_stage",instrain_compare1$group),]
lm(data = intra_later,formula = conANI~year_distance)%>%summary()
p2=ggplot(instrain_compare1,aes(year_distance,conANI,color=group))+
  geom_point()+
  theme_bw()+geom_smooth(data = instrain_compare1,
                         aes(year_distance,conANI,group=group),
                         formula=y~x,level=0,method='lm')
ggsave("conANI_lm.pdf",plot = p2,device = "pdf",width = 5,height = 3)
