length=read.csv("Contigs_quality_summary_filtering4.tsv",sep='\t',header = F)
length=length[-grep("high kmer",length$V14),]
length2=as.data.frame(table(cut(length$V2,breaks = c(200000,250000,300000,350000,400000,450000,500000,600000,700000,1000000,1500000),labels = c(
                                       "200","250","300","350","400","450","500","600","700","1000"),ordered_result = T)))

mag=read.csv("MAGs_quality_summary.tsv",sep='\t',header = T)
mag=mag[-grep("high kmer",mag$warnings),]
mag=mag[mag$contamination<2,]
mag=mag[-grep("contig",mag$warnings),]
completeness=read.table("completeness.tsv",header = T)
completeness$aai_completeness
mag2=as.data.frame(table(cut(mag$contig_length,breaks = c(200000,250000,
                                                    300000,350000,400000,450000,500000,600000,700000,1000000,1500000),labels = c(
                                                      "200","250",
                                                      "300","350","400","450","500","600","700","1000"
                                                    ),ordered_result = T)))
length2$mag_Freq=mag2$Freq
write.table(mag,"MAGs_quality_summary_filtering.tsv",quote = F,sep='\t',row.names = F)
library(ggplot2)
library(reshape2)
length2=melt(length2)
ggplot(length2)+
  geom_bar(stat="identity",aes(Var1,value,fill=variable),alpha=0.8,width=0.8)+
  theme_classic()+labs(y="Count",x="Scaffold length (kbp)")
ggsave("figure1_length_count.pdf",width = 7,height = 4)
