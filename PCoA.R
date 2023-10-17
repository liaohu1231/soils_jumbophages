data=read.csv("all.jumphages.rpkm_drep.tsv",header = 1,row.names = 1,check.names = F,sep = '\t')
group=read.table("metadata.txt",header = 1,stringsAsFactors = T)

library(ggplot2)
relative_abundance=function(d){
  dta_sum=apply(d,2,function(x){x/sum(x)})
}
data_rel=as.data.frame(relative_abundance(data))
##PCoA
library(dplyr)
library(vegan)
results5=t(data_rel)
results5=as.matrix(results5)
library(vegan)
distance <- as.matrix(vegdist(results5, method= "bray",na.rm = T))
colnames(distance)=rownames(results5)
rownames(distance)=rownames(results5)
distance2=matrix(NA,nrow(distance),ncol(distance))
for(i in 1:nrow(distance)){
  for(j in 1:ncol(distance)){
    if ( i != j){
      distance2[i,j]=distance[i,j]
    }
  }
}
rownames(distance2)=rownames(distance)
colnames(distance2)=colnames(distance)

library(pheatmap)
p=pheatmap(distance2,cluster_rows = F,cluster_cols = F,clustering_method = "average",
         na_col = "lightblue",border_color = F,cellwidth = 10,cellheight = 8,
         color = colorRampPalette(c("white","white","white","white",
                                    "white","white","white","white",
                                    "white","white","white","white","white","white",
                                    "white","white","white","#40CCA7"))(100))
p
ggsave("bray_distance_pheatmap.pdf",plot = p,device = "pdf",width = 6,height = 5)

pcoa_v <- cmdscale(distance, k = (nrow(results5) - 1), eig = TRUE)
point_v <- data.frame(pcoa_v$point)
#species <- wascores(pcoa$points[,1:2], results5)
pcoa_eig_v <- (pcoa_v$eig)[1:2] / sum(pcoa_v$eig)
sample_site_v <- data.frame ({pcoa_v$point})[1:3]
sample_site_v$names <- rownames(sample_site_v)
names(sample_site_v)[1:3] <- c('PCoA1', 'PCoA2','PCoA3')

write.table(sample_site_v,"viral_pcoa123.txt",sep = '\t',quote = F)

nd=sample_site_v[,1:3]
group$stage=rep(c("early_stage","later_stage"),c(12,12))

p_value=anosim(t(data_rel),group$stage,permutations = 9999)
p_value=anosim(t(data_rel),group$sites,permutations = 9999)
summary(p_value)

dist = vegdist(t(data_rel), method = 'bray')

p_value_1=adonis2(formula = dist~sites,group,permutations = 9999)
p_value_1
p_value_2=adonis2(formula = dist~year,group,permutations = 9999)
p_value_2
p_value_3=adonis2(formula = dist~stage,group,permutations = 9999)
p_value_3
dist_year=vegdist(group$year)
mantel(distance,dist_year)

##sum_abu
sum_abu=data.frame(sum_abu=apply(data, 2, sum))
rownames(group)=group$sample
sum_abu$sites=group[rownames(sum_abu)]
abu_sites=data.frame(average_abu=tapply(sum_abu$sum_abu,sum_abu$site,mean))
abu_sites$sd=tapply(sum_abu$sum_abu,sum_abu$site,sd)
abu_sites$sites=rownames(abu_sites)
library(ggplot2)
ggplot(abu_sites,aes(sites,average_abu))+
  geom_bar(stat = "identity",fill="#E377C2",color="black",
           alpha=0.7,width = 0.8)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=average_abu-sd,ymax=average_abu+sd),
                size=0.5,width=0.4)+
  theme_classic()

abu2=abu
abu2=as.data.frame(relative_abundance(abu2))
abu2$distribution=length[rownames(abu2),2]
source("sum_sam_colnames.R")
abu2=sum_same_col(x1 = abu2,25)
library(reshape2)
library(ggpubr)
abu2=melt(abu2)
ggplot(abu2,aes(Var1,value*100))+
  geom_boxplot(outlier.size = 0,outlier.colour = "white",fill="#C4E1E1")+
  geom_jitter(height = 0,width = 0.2,size=5,alpha=0.4)+
  theme_classic()+stat_compare_means(paired = T,label = "p.signif",
                                     comparisons = list(c("local","global")))+
  labs(y="Relative abundance",x="")

abu3=as.data.frame(t(data[,1:24]))
abu3$site=group[,2]
abu3=sum_same_col(abu3,ncol(abu3))
abu3=as.data.frame(t(abu3))
number=as.data.frame(apply(abu3, 1, function(x){sum(x>50)}))
nrow(as.data.frame(number[number[,1]>1,]))

##PCoA
me = unique(nd$type)
result = data.frame()
otu=t(data)
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij = subset(nd, type %in% c(as.character(group_name[i]), as.character(group_name[j])))
    otu_ij = otu[rownames(group_ij), ]
    adonis_result_otu_ij = adonis(otu_ij~group, group_ij, permutations = 999, distance = 'bray')
    res.temp = as.data.frame(adonis_result_otu_ij$aov.tab)[1,]
    rownames(res.temp) = paste(as.character(group_name[i]),'/',as.character(group_name[j]))
    result = rbind(result,res.temp)
  }
}
head(result,nrow(result))

nd=matrix(data = NA,nrow = nrow(sample_site_v),ncol = 1)
summary(p_value)
sample_site_v$site=group[,2]
sample_site_v$stage=group[rownames(sample_site_v),4]
library(ggplot2)
library(ggsci)
p <- ggplot(sample_site_v, aes(PCoA1, PCoA2)) +
  theme(panel.grid = element_line(color = 'gray', linetype = 2), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(aes(color=site),size = 6, alpha = 0.5) + 
  scale_color_d3() + 
  scale_shape_manual(values = c(21,24))+
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig_v[1], 2), '%'), 
       y = paste('PCoA axis2: ', round(100 * pcoa_eig_v[2], 2), '%'))+
  #stat_ellipse(data = sample_site_v,mapping = aes(PCoA1, PCoA2,group = site),level = 0.95, show.legend = TRUE,inherit.aes = F)+
  #stat_ellipse(data = sample_site_v,mapping = aes(PCoA1, PCoA2,group = site),level = 0.975, show.legend = TRUE,inherit.aes = F)+
  annotate('text', colour="#8766d")
p

ggsave("PCoA.pdf",plot = p,device = "pdf",width = 6,height = 4.5)


###shared populations
source('new.cbind.R')
nd=matrix(NA,0,0)
abu_vc_2=data_rel
for (i in 1:ncol(abu_vc_2)){
  col=as.matrix(abu_vc_2[,i])
  rownames(col)=rownames(abu_vc_2)
  col=as.matrix(col[col[,1]>0,])
  colnames(col)=colnames(abu_vc_2)[i]
  nd=as.matrix(new.cbind(nd,col))
}

for(i in 1:nrow(nd)){
  for (j in 1:ncol(nd)){
    if(is.na(nd[i,j])==FALSE){
      nd[i,j]=1
    }
    else{
      nd[i,j]=0
    }
  }
}

share=matrix(NA,ncol(nd),ncol(nd))
for (i in 1:ncol(nd)){
  for (j in 1:ncol(nd)){
    print(j)
    sum=as.matrix(apply(nd[,c(i,j)],1,sum))
    shared=as.data.frame(nrow(nd[sum==2,]))
    if(nrow(shared)>0){
      a=sum(nd[,i])
      b=sum(nd[,j])
      viral_shared_content=((shared/a)+(shared/b))/2
      share[i,j]=viral_shared_content[1,1]
    }else{
      share[i,j]==0
    }
  }
}

colnames(share)=colnames(nd)
rownames(share)=colnames(nd)
share[is.na(share)]=0

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
share1=get_lower_tri(share)

group$stage=rep(x = c("early","later"),c(12,12))

library(reshape2)
share1=melt(share1)
share1=na.omit(share1)
share1=share1[share1$value!=1,]
share1=as.data.frame(share1)

nd=matrix(NA,nrow(share1),1)
for (i in 1:nrow(share1)){
  print(i)
  name=as.matrix(group[grep(share1[i,1],group$sample),4])
  name_1=as.matrix(group[grep(share1[i,2],group$sample),4])
  if ( name_1[1,1]==name[1,1]){
    if (name_1[1,1]=='early'){
     nd[i,1]='intra_early_stage' 
    }
    if (name_1[1,1]=="later"){
      nd[i,1]='intra_later_stage'
    }
  }else{
    nd[i,1]='inter_stage'
  }

}

share1$zone=nd
library(ggplot2)
library(ggpubr)
share1=as.data.frame(share1)
p=ggplot(share1,aes(zone,value*100))+
  geom_boxplot(fill="lightblue",outlier.size = 0,outlier.colour = 'white')+
  geom_jitter(width = 0.1,alpha=0.7,size=2)+
  theme(axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60,size = 10,hjust = 1),
        legend.key.size = unit(5, 'mm'),
        legend.text = element_text(size = 10))+
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.2), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4)
p
p+stat_compare_means(value ~ group, data = share1,
                     method = "wilcox.test",label="p.signif",paired = F,
                     comparisons =list(c("intra_early_stage","inter_stage"
                     ),c("intra_later_stage","intra_early_stage")))
ggsave("share_boxplot.pdf",device = "pdf",width = 3.5,height = 5)

library(pheatmap)
annotation_row=as.data.frame(group_non[,c(1,2)])
rownames(annotation_row)=group_non[,4]
colnames(annotation_row)=c("type","zone")

p=pheatmap(share,
         border_color = "grey",
         color = colorRampPalette(c("white","blue","yellow","lightblue","black"))(100),
         cluster_rows = F,cluster_cols = F,
         clustering_method = "average")
p
rownames(group)=group$sample
ggsave(filename = "share_heatmap.pdf",plot = p,device = 'pdf',width = 6,height = 5)

shannon=data.frame(shannon_index=diversity(t(data_rel), "shannon"))
shannon$site=group[rownames(shannon),2]
shannon_average=data.frame(shannon_average=tapply(shannon$shannon_index,shannon$site,mean))
shannon_average$sd=tapply(shannon$shannon_index,shannon$site,sd)
group2=group[!duplicated(group$sites),]
rownames(group2)=group2$sites
shannon_average$year=group2[rownames(shannon_average),3]
model <- loess(shannon_average~year, data = shannon_average)
summary(model)
sum_abu$year=group[rownames(sum_abu),3]
shannon_average$site=group2[rownames(shannon_average),2]
ggplot(data = shannon_average,mapping = aes(x = year,y=shannon_average))+
  theme_classic()+
  labs(y="shannon index")+
  geom_smooth(method = 'loess', formula = y~x, span = 1,color="blue")+
  geom_point(stat = 'identity',size=4,shape=21,fill="lightblue")+
  geom_errorbar(aes(ymin=shannon_average-sd,ymax=shannon_average+sd))+
  theme_classic()+geom_line()+
  geom_text(aes(label=site),size=4,vjust=1.7)+
  labs(y="shannon index")
ggsave("shannon_index.pdf",device = 'pdf',width = 6,height = 4.5)

count=data.frame(count=apply(data,2,function(x){sum(x>0)}))
count$site=group[rownames(count),2]
count_average=data.frame(count_average=tapply(count$count,count$site,mean))
count_average$sd=tapply(count$count,count$site,sd)
count_average$site=group2[rownames(count_average),2]
count_average$year=group2[rownames(count_average),3]
ggplot(count)
ggplot(data = count_average,mapping = aes(x = year,y=count_average))+
  theme_classic()+
  labs(y="shannon index")+
  geom_smooth(method = 'loess', formula = y~x, span = 1,color="blue")+
  geom_point(stat = 'identity',size=4,shape=21,fill="lightblue")+
  geom_errorbar(aes(ymin=count_average-sd,ymax=count_average+sd))+
  theme_classic()+geom_line()+
  geom_text(aes(label=site),size=4,vjust=1.7)+
  labs(y="shannon index")
ggsave("jumbo_count.pdf",device = 'pdf',width = 6,height = 4.5)


###abu_frequency
data_rel=as.data.frame(relative_abundance(data))
data_rel$frequency=apply(data_rel,1,function(x){sum(x>0)})
data_rel$frequency=as.factor(data_rel$frequency)
data_rel$contigs=rownames(data_rel)
library(reshape2)
data_rel=melt(data_rel)
data_rel=data_rel[data_rel$value>0,]
data_rel$site=group[data_rel$variable,2]
library(ggplot2)
ggplot(data_rel,aes(frequency,log(value),fill=site))+
  geom_point(stat = 'identity',shape=21)



