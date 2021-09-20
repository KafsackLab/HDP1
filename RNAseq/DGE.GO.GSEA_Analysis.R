library(ggplot2); library(reshape2); library(dplyr); library(ggrepel); library(fgsea);library(data.table);

setwd("/Users/bkafsack/kafsacklab@gmail.com - Google Drive/My Drive/Projects/Riward/PfSneezy/RNASeq")
dge<-read.delim("/Users/bkafsack/kafsacklab@gmail.com - Google Drive/My Drive/Projects/Riward/PfSneezy/RNASeq/HDP1_DGE.tab",header=T,sep="\t")
het<-read.csv("/Users/bkafsack/kafsacklab@gmail.com - Google Drive/My Drive/Knowledge Base/GeneLists/GeneInfo_Heterochromatin.csv")
geneInfo<-read.csv("/Users/bkafsack/kafsacklab@gmail.com - Google Drive/My Drive/Knowledge Base/GeneLists/GeneInfoTable/GeneInfoTable_expanded_May2017_PlasmoDB32.csv")
all_imc<-c("PF3D7_0109000","PF3D7_0303400","PF3D7_0304000","PF3D7_0304100","PF3D7_0423500","PF3D7_0510300","PF3D7_0515700","PF3D7_0518900","PF3D7_0522600","PF3D7_0525800","PF3D7_0611600","PF3D7_0621400","PF3D7_0708600","PF3D7_0708800","PF3D7_0722900","PF3D7_0723300","PF3D7_0805200","PF3D7_0816500","PF3D7_0823500","PF3D7_0918000","PF3D7_1003600","PF3D7_1011000","PF3D7_1017500","PF3D7_1028900","PF3D7_1031200","PF3D7_1115900","PF3D7_1141900","PF3D7_1221400","PF3D7_1222700","PF3D7_1246400","PF3D7_1323700","PF3D7_1341500","PF3D7_1341800","PF3D7_1342600","PF3D7_1345600","PF3D7_1351700","PF3D7_1355600","PF3D7_1406800","PF3D7_1407700","PF3D7_1417000","PF3D7_1430800","PF3D7_1431100","PF3D7_1431900","PF3D7_1445700","PF3D7_1447500","PF3D7_1460600")

rownames(dge)<-dge$gene

dge$fpkm_KO[dge$fpkm_KO==0]=min(dge$fpkm_KO[dge$fpkm_KO!=0])
dge$fpkm_KO[dge$fpkm_WT==0]=min(dge$fpkm_KO[dge$fpkm_WT!=0])
dge$logX<-log2(dge$fpkm_KO/dge$fpkm_WT)

dge$fam<-geneInfo$family[match(rownames(dge),geneInfo$geneID)];
dge$CVfam<-"other"
dge$CVfam[dge$fam %in% c("var","rifin","stevor","surfin","PHISTa","PHISTb","PHISTc")]<-as.character(dge$fam[dge$fam %in% c("var","rifin","stevor","surfin","PHISTa","PHISTb","PHISTc")])
dge$CVfam[dge$fam %in% c("PHISTa","PHISTb","PHISTc")]<-"PHIST"

dge$hc<-as.logical(geneInfo$heterochrom[match(rownames(dge),geneInfo$geneID)])
head(dge)
dge$pseudo[dge$description=="erythrocyte membrane protein 1 (PfEMP1), exon 2"]<-"Yes"
dge$fam[dge$description=="erythrocyte membrane protein 1-like"]<-"other"
dge<-dge[dge$pseudo=="No",]

dge2<-dge[dge$fpkm_KO>=10 |dge$fpkm_WT>=10,] # at least one sample must have mean of FPKM>10
rownames(het)<-het$geneID

plotme<-dge2
plotme$q_value[plotme$q_value==min(plotme$q_value)]<-plotme$q_value[plotme$q_value==min(plotme$q_value)]*runif(length(plotme$q_value[plotme$q_value==min(plotme$q_value)]),0.95,1.05) # add 5% jitter to lowest qval
plotme$logX[plotme$logX > 4]<- 4
plotme$logX[plotme$logX < -4]<- -4

plotme$het<-0
plotme$het[plotme$gene %in% het$geneID] = het[plotme$gene[plotme$gene %in% het$geneID],]$H3K9me3_Bunnik_Asex+het[plotme$gene[plotme$gene %in% het$geneID],]$H3K9me3_Bunnik_GC*2
plotme$het2<-0
plotme$het2[plotme$gene %in% het$geneID] = het[plotme$gene[plotme$gene %in% het$geneID],]$HP1_Fraschka_Asex +het[plotme$gene[plotme$gene %in% het$geneID],]$HP1_Fraschka_eGC*2
plotme$het3<-0
plotme$het3[plotme$gene %in% het$geneID] = het[plotme$gene[plotme$gene %in% het$geneID],]$HP1_Flueck_Asex
rownames(plotme)<-plotme$gene

plotme$het[is.na(plotme$het)]<-0;plotme$het2[is.na(plotme$het2)]<-0;plotme$het3[is.na(plotme$het3)]<-0
hetgenes<-((plotme$het>0)+(plotme$het2>0)+(plotme$het3>0))>1
names(hetgenes)<-rownames(plotme)
hetgenes<-names(hetgenes[hetgenes==T])
plotme$allhet<-0; plotme$allhet[plotme$gene %in% hetgenes]<-1

logX<-plotme$logX; names(logX)<-plotme$gene
logX[logX==Inf]=max(logX[logX!=Inf])
logX[logX==-Inf]=min(logX[logX!=-Inf])
logX<-logX[order(logX)]
ranks<-rank(logX,ties.method = "random")
plotme2<-data.frame(rank=max(ranks)-ranks,logX=logX,row.names = names(logX))

################ Volcano Plot
ggplot()+
  geom_hline(yintercept = -log10(0.05),color="black",linetype="dashed")+
  geom_point(data=plotme,aes(x=logX,y=-log10(q_value)),size=1,fill="white",shape=21)+
  geom_point(data=plotme[plotme$gene%in% hetgenes,],aes(x=logX,y=-log10(q_value)),size=2,fill="cyan",shape=21)+
  geom_point(data=plotme[plotme$gene%in% all_imc,],aes(x=logX,y=-log10(q_value)),size=2,fill="orange",shape=21)+
  geom_text_repel(data=plotme[plotme$gene=="PF3D7_1466200",],aes(x=logX,y=-log10(q_value)),label="hdp1", nudge_y=-0.13, nudge_x=-0.13)+
  geom_point(data=plotme[plotme$gene=="PF3D7_1466200",],aes(x=logX,y=-log10(q_value)),size=2,fill="black",shape=21)+
  geom_text_repel(data=plotme[plotme$gene=="PF3D7_0936600",],aes(x=logX,y=-log10(q_value)),label=c("gexp5"), nudge_x=+0.6, nudge_y = -.3, segment.color = "red", color="red")+
  geom_text_repel(data=plotme[plotme$gene =="PF3D7_0406200",],aes(x=logX,y=-log10(q_value)),label=c("pfs16"), nudge_x=-0.5, nudge_y = -.18, segment.color = "red", color="red")+
  geom_point(data=plotme[plotme$gene %in% c("PF3D7_0936600","PF3D7_0406200"),],aes(x=logX,y=-log10(q_value)),size=2,fill="red",shape=21)+
  scale_x_continuous(limits = c(-4.1, 4.1),breaks=-4:4,labels =c("<1/16","1/8","1/4","1/2",1,2,4,8,">16"))+
  scale_y_continuous(breaks=-log10(c(1,0.05,0.01,0.002)),labels =c(1,0.05,0.01," <0.002"))+
  labs(x="KO vs Parental Fold-Change", y="False Discovery Rate")+
  theme(
    text = element_text(colour = "black"),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
    axis.text = element_text(colour = "black", size=10)
  )
ggsave("/Volumes/GoogleDrive-100638174929531249371/My Drive/Manuscripts+Abstracts/Manuscripts/Submitted/Riward_HDP1 paper/DGE&GSEA/VolcanoPlot_all.pdf",height=5,width=10)


####### subtelomeric vs centromeric var/rif

z<-dge[dge$fam %in% c("var","rifin","stevor","PHISTa","PHISTb","PHISTc"),]
z$subtelo.HP1<- factor(z$subtelo.HP1,levels=c(TRUE,FALSE))
z$logX[z$logX>4]<-4; z$logX[z$logX<(-4)]<-(-4);

ggplot(z[z$fam %in% c("var","rifin"),], aes(x=fam,y=logX,fill=fam, lty=subtelo.HP1))+
  geom_hline(yintercept = 0, lty="dotted")+
  geom_violin(aes(color=fam), position=position_dodge(width = 1), alpha=0.3, color="black")+
  geom_boxplot(aes(color=fam), position=position_dodge(width = 1), alpha=0.3, color="black", width=0.3)+
  geom_point(shape=21,position=position_jitterdodge(jitter.width =0.2,dodge.width =1))+
  theme_classic()+
  scale_fill_manual(values=c("mediumorchid","chartreuse"),labels=c("rifins", "var"))+
  scale_linetype_manual(values=c("solid","dashed"),labels=c("subtelomeric", "non-subtelomeric"))+
  scale_y_continuous(breaks=-2:4,labels =c("<1/4","1/2",1,2,4,8,">=16"), limits=c(-2.2,4.5))+
  scale_x_discrete(labels=c("rifins","vars"))+
  labs(y="HDP1 KO vs Parental fold-change", x="", fill="",lty="")+
  theme(legend.position = "none")
ggsave("/Volumes/GoogleDrive-100638174929531249371/My Drive/Manuscripts+Abstracts/Manuscripts/Submitted/Riward_HDP1 paper/DGE&GSEA/Chr.Position.pdf", height= 4, width=5)

ggplot(z, aes(x=fam,y=logX,fill=fam))+
  geom_hline(yintercept = 0, lty="dotted")+
  geom_violin(aes(color=fam), position=position_dodge(width = 1), alpha=0.3, color="black")+
  geom_point(shape=21,position=position_jitterdodge(jitter.width =0.2,dodge.width =1))+
  theme_classic()+
  #scale_fill_manual(values=c("mediumorchid","chartreuse"),labels=c("subtelomeric", "non-subtelomeric"))+
  scale_y_continuous(breaks=-3:4,labels =c("<1/8","1/4","1/2",1,2,4,8,">=16"), limits=c(-3.1,4.5))+
  #scale_x_discrete(labels=c("inner\nmembrane\ncomplex","rifins","stevors","vars"))+
  labs(y="HDP1 KO vs Parental fold-change", x="", fill="chromosomal\nposition")+
  theme(legend.position = "none")
ggsave("/Volumes/GoogleDrive-100638174929531249371/My Drive/Manuscripts+Abstracts/Manuscripts/Submitted/Riward_HDP1 paper/DGE&GSEA/gene_fam.pdf", height= 2.5, width=8)

t.test(dge$logX[dge$fam=="var"])
t.test(dge$logX[dge$fam=="rifin"])
t.test(dge$logX[dge$fam=="stevor"])
t.test(dge$logX[dge$fam=="Pfmc-2TM"])
t.test(z$logX[z$fam=="PHISTa"])
t.test(dge$logX[dge$fam=="PHISTb"])
t.test(dge$logX[dge$fam=="PHISTc"])
t.test(dge$logX[dge$fam=="Pfmc-2TM"])

####### Downregulated IMC proteins      #########
imc_down<-c("PF3D7_1430800","PF3D7_1431100","PF3D7_1345600","PF3D7_1323700","PF3D7_1221400","PF3D7_0109000","PF3D7_0717600","PF3D7_1011000","PF3D7_1222700","PF3D7_0621400","PF3D7_1406800","PF3D7_0918000","PF3D7_0423500")

z<-plotme[plotme$gene %in% imc_down,]
z=z[,c("gene","fpkm_KO","fpkm_WT","fam","name")]
z=melt(z,id.vars=c("gene","fam","name"))
z$variable=factor(z$variable,levels=c("fpkm_WT","fpkm_KO"))
z$name[z$name==""]<-z$gene[z$name==""]

ggplot(z,aes(x=variable, fill=variable,color=name, y=value,group=gene))+
  geom_line(color="black",size=0.5)+
  geom_point(shape=1,size=3,color="black",fill="white")+
  scale_y_continuous(limits=c(30,500))+
  theme_classic()+
  xlab("")+
  geom_text_repel(aes(label=ifelse(variable=="fpkm_WT",paste0(name),"")),nudge_x = -0.2,color="black")+
  scale_x_discrete(labels=c("WT","HDP1 KO"))+
  ylab("fpkm")+
  theme(axis.text = element_text(colour = "black",size=12))

ggsave("/Volumes/GoogleDrive-100638174929531249371/My Drive/Manuscripts+Abstracts/Manuscripts/Submitted/Riward_HDP1 paper/DGE&GSEA/IMC_down.pdf",height=4,width=4)

######  GO Analysis with fgsea
set.seed(42)
pathways<-gmtPathways("/Volumes/GoogleDrive-100638174929531249371/My Drive/Manuscripts+Abstracts/Manuscripts/Submitted/Riward_HDP1 paper/DGE&GSEA/Pf3D7.gmt")
pathways$hetgenes<-hetgenes

fgseaRes <- fgsea(pathways, logX, minSize=15, maxSize=200,eps=0.0, scoreType="std")
fgseaRes<-fgseaRes[fgseaRes$padj<=0.05,]
fgseaRes<-fgseaRes[order(fgseaRes$padj),]

fgseaRes.df<-data.frame(fgseaRes[,1:7])
fgseaRes.df$set.genes<-unlist(lapply(fgseaRes$leadingEdge,paste,collapse=", "))

write.csv(fgseaRes.df,"/Volumes/GoogleDrive-100638174929531249371/My Drive/Manuscripts+Abstracts/Manuscripts/Submitted/Riward_HDP1 paper/DGE&GSEA/GSEA.csv")

### IMC genes GO
Ecurve<-plotEnrichment(pathways[["GO:0070258"]],stats=logX)$data

# ggplot()+
#   geom_tile(data=plotme2,aes(y=mean(c(max(Ecurve$y),min(Ecurve$y))),x=rank,fill=logX,height=max(Ecurve$y)-min(Ecurve$y)))+
#   scale_fill_gradient2(high="red",low="blue",breaks=c(4,2,0,-2,-4), labels=c(">4",2,0,-2,"<-4"))+
#   geom_hline(yintercept = 0, lty="dotted")+
#   geom_line(data=Ecurve,aes(x,y))+
#   theme_classic()+
#   scale_x_continuous(expand = c(0,0), labels=NULL, breaks = c(Ecurve$x[2:(length(Ecurve$x)-1)]), sec.axis = sec_axis( trans=~.*1, name="",breaks=c(320,4620),labels=c("HDP1 KO","WT")))+
#   scale_y_continuous(expand = c(0,0))+
#   theme(
#     axis.ticks.x.bottom = element_line(color = "black",size=0.2), axis.ticks.length.x.bottom = unit(0.2,"cm"),
#     axis.ticks.x.top = element_blank(),
#     axis.text.x.top = element_text(color = "black",size=12),
#     axis.title.x.bottom = element_text(color = "black",size=14, vjust=-1, hjust=0),
#     axis.text.y = element_text(color = "black",size=10),
#     legend.title.align = 0.5, legend.title=element_text(size=8), legend.position = "none"
#   )+
#   labs(y="Enrichment score",x="inner membrane complex genes, NES= -1.69, FDR= 0.02")

ggplot()+
  geom_tile(data=plotme2,aes(y=-mean(c(max(Ecurve$y),min(Ecurve$y))),x=-rank,fill=logX,height=max(Ecurve$y)-min(Ecurve$y)))+
  scale_fill_gradient2(high="red",low="blue",breaks=c(4,2,0,-2,-4), labels=c(">4",2,0,-2,"<-4"))+
  geom_hline(yintercept = 0, lty="dotted")+
  geom_line(data=Ecurve,aes(-x,-y))+
  theme_classic()+
  scale_x_continuous(expand = c(0,0), labels=NULL, breaks = -c(Ecurve$x[2:(length(Ecurve$x)-1)]), sec.axis = sec_axis( trans=~.*1, name="",breaks=c(-4120,-350),labels=c("WT","HDP1 KO")))+
  scale_y_continuous(expand = c(0,0))+
  theme(
    axis.ticks.x.bottom = element_line(color = "black",size=0.2), axis.ticks.length.x.bottom = unit(0.2,"cm"),
    axis.ticks.x.top = element_blank(),
    axis.text.x.top = element_text(color = "black",size=12),
    axis.title.x.bottom = element_text(color = "black",size=14, vjust=-1, hjust=0),
    axis.text.y = element_text(color = "black",size=10),
    legend.title.align = 0.5, legend.title=element_text(size=8), legend.position = "none"
  )+
  labs(y="Enrichment score",x="gene rank\ninner membrane complex genes, NES= 1.69, FDR= 0.02")

ggsave("/Volumes/GoogleDrive-100638174929531249371/My Drive/Manuscripts+Abstracts/Manuscripts/Submitted/Riward_HDP1 paper/DGE&GSEA/IMC_GSEAplot.pdf", width=6, height=2)

######## Heterochromatin ############
Ecurve<-plotEnrichment(hetgenes,stats=logX)$data

ggplot()+
  geom_tile(data=plotme2,aes(y=mean(c(max(Ecurve$y),min(Ecurve$y))),x=-rank,fill=logX,height=max(Ecurve$y)-min(Ecurve$y)))+
  scale_fill_gradient2(high="red",low="blue",breaks=c(4,2,0,-2,-4), labels=c(">4",2,0,-2,"<-4"))+
  geom_hline(yintercept = 0, lty="dotted")+
  geom_line(data=Ecurve,aes(-x,y))+
  theme_classic()+
  scale_x_continuous(expand = c(0,0), labels=NULL, breaks = -c(Ecurve$x[2:(length(Ecurve$x)-1)]), sec.axis = sec_axis( trans=~.*1, name="",breaks=c(-4120,-350),labels=c("WT","HDP1 KO")))+
  scale_y_continuous(expand = c(0,0))+
  theme(
    axis.ticks.x.bottom = element_line(color = "black",size=0.2), axis.ticks.length.x.bottom = unit(0.2,"cm"),
    axis.ticks.x.top = element_blank(),
    axis.text.x.top = element_text(color = "black",size=12),
    axis.title.x.bottom = element_text(color = "black",size=14, vjust=-1, hjust=0),
    axis.text.y = element_text(color = "black",size=10),
    legend.title.align = 0.5, legend.title=element_text(size=8), legend.position = "none"
  )+
  labs(y="Enrichment score",x="gene rank\nHeterchromatinized genes, NES = 2.63, FDR< 0.001")
  ggsave("/Volumes/GoogleDrive-100638174929531249371/My Drive/Manuscripts+Abstracts/Manuscripts/Submitted/Riward_HDP1 paper/DGE&GSEA/Hetchrom_GSEAplot.pdf", width=6, height=2)

  