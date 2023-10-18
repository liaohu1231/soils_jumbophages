library(picante)
library(ggplot2)
library(ggsci)
tree=read.tree("concat.faa.treefile")
data=read.table("all.jumphages.rpkm_drep.tsv",header = 1,row.names = 1,check.names = F)
relative_abundance=function(d){
  dta_sum=apply(d,2,function(x){x/sum(x)})
}
data=relative_abundance(data)
sum=apply(data,1,sum)
data1=data[sum==0,]
prune_tree1<-prune.sample(t(data1),tree)
dist_sum1=sum(as.data.frame(cophenetic(prune_tree1)))
dist_sum=sum(as.data.frame(cophenetic(tree)))
dist_sum/dist_sum1

abu=as.data.frame(t(data))
group=read.csv("../metadata.txt",sep = '\t',header = 1,row.names = 1)
abu$sites=group[,1]
source("sum_sam_colnames.R")
sites_abu=sum_same_col(abu,1844)


prune_tree<-prune.sample(sites_abu,tree)
dist <- as.data.frame(cophenetic(prune_tree))
sites_abu=sites_abu[,colnames(dist)]
library(dplyr)
mpd=ses.mpd(sites_abu, dist, null.model = "richness",abundance.weighted = F, runs = 999, iterations = 1000)

ses.mntd=ses.mntd(sites_abu, dist, null.model = "richness",abundance.weighted = T, runs = 999, iterations = 1000)

ses.pd=ses.pd(sites_abu,tree,include.root = F, null.model = "richness",
             runs = 999, iterations = 1000)
pd=pd(samp = sites_abu,tree = tree,include.root = F)

mpd$sites=rownames(mpd)
ses.mntd$sites=rownames(ses.mntd)

ggplot(mpd,aes(sites,mpd.obs.z,fill=sites))+
  geom_bar(stat = "identity")+
  theme_classic()+
  scale_fill_d3()

ggplot(ses.mntd,aes(sites,mntd.obs.z,fill=sites))+
  geom_bar(stat = "identity")+
  theme_classic()+
  scale_fill_d3()

ses.pd$sites=rownames(ses.pd)
ggplot(ses.pd,aes(sites,pd.obs,fill=sites))+
  geom_bar(stat = "identity")+
  theme_classic()+
  scale_fill_d3()

group1=read.csv("group1.txt",sep = '\t',header = F,row.names = 1)
ses.mntd$year=group1[rownames(ses.mntd),1]
mpd$year=group1[rownames(mpd),1]
summary(lm(year~mntd.obs.z,ses.mntd))
summary(lm(year~mpd.obs.z,mpd))

linearplot=data.frame(year=mpd$year,sites=mpd$sites,mpd=mpd$mpd.obs.z,mntd=ses.mntd$mntd.obs.z)
ggplot(linearplot)+
  geom_point(aes(x = year,y = mpd),shape=21,size=5,fill="#89B510")+
  geom_point(aes(x = year,y = mntd),shape=21,size=5,fill="red")+
  labs(y="Phylogenetic diversity")+
  theme_classic()


library(phylocomr)
library(vegan)
library(reshape2)

inter_mpd=as.matrix(comdistnt(comm = sites_abu,dis = dist,abundance.weighted = T))
inter_mpd1=melt(inter_mpd)
inter_mpd1=inter_mpd1[inter_mpd1$value>0,]
ggplot(inter_mpd1,aes(Var1,value))+geom_boxplot()+
  theme_classic()+labs(x="sites",y="MNTD")
inter_mpd=as.dist(inter_mpd)
mpd_hclust=hclust(inter_mpd,method="complete")

plot(mpd_hclust)

plot(ses.mntd$mntd.obs.z ~ group1$V2, xlab = "Sites", ylab = "SES(MNTD)")
abline(h = 0, col = "gray")

