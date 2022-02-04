### SETUP ####
setwd("/mnt/raid5/ChIPseqArchive/LeRoch_HDP1_GFP/")
require(rtracklayer);require(plyr);require(ggplot2); library(GenomicRanges); require(Gviz);require(ChIPpeakAnno);require(reshape2);require(motifmatchr);require(plyranges); require(scales)

# get 3D7 genome info
library(BSGenome.Pf3D7.PlasmoDB.51)
Pf3D7.seqinfo<-data.frame(cbind(Chr=seqlevels(BSGenome.Pf3D7.PlasmoDB.51),chr=c(paste0("Chr",substring(names(seqinfo(BSGenome.Pf3D7.PlasmoDB.51))[1:14],7,8)),"Mito","Apico"), data.frame(seqinfo(BSGenome.Pf3D7.PlasmoDB.51))))

# function switches between PlasmoDB and ChIPpeakAnno-compatible ASM9120v1 seqlevelsStyle, see seqlevelsStyle() 
renameChrs<-function(gr){  
  gr<-data.frame(gr)
  gr$seqnames<-as.character(gr$seqnames)
  if(grepl("Pf3D7_",gr$seqnames[1],fixed=T)){
    gr$seqnames<-Pf3D7.seqinfo$chr[match(gr$seqnames,Pf3D7.seqinfo$Chr)]
    gr<-makeGRangesFromDataFrame(data.frame(gr),keep.extra.columns = T)
    seqlengths(gr)<-Pf3D7.seqinfo$seqlengths[match(names(seqlengths(gr)),Pf3D7.seqinfo$chr)]
    isCircular(gr)<-Pf3D7.seqinfo$isCircular[match(names(seqlengths(gr)),Pf3D7.seqinfo$chr)]
    genome(gr)<-"51"
  }else{
    gr$seqnames<-Pf3D7.seqinfo$Chr[match(gr$seqnames,Pf3D7.seqinfo$chr)]
    gr<-makeGRangesFromDataFrame(data.frame(gr),keep.extra.columns = T)
    seqlengths(gr)<-Pf3D7.seqinfo$seqlengths[match(names(seqlengths(gr)),Pf3D7.seqinfo$Chr)]
    isCircular(gr)<-Pf3D7.seqinfo$isCircular[match(names(seqlengths(gr)),Pf3D7.seqinfo$Chr)]
    genome(gr)<-"51"
  }
  return(gr)
}

# import GFF
GFF<-sort(import.gff3("./genome/PlasmoDB-51_Pfalciparum3D7.gff"))
seqlengths(GFF)[1:16]<-seqlengths(Pf3D7)[c(1:14,16,15)]
isCircular(GFF)<-c(rep(F,14),T,T)
genome(GFF)<-genome(Pf3D7)
seqinfo(GFF)

# modify GFF for compatibility with later functions
gff<-renameChrs(GFF)
gff$feature<-gff$type
names(gff)<-gff$ID

gff$gene<-substr(gff$ID,regexpr("PF3D7_",gff$ID),regexpr("PF3D7_",gff$ID)+12)
gff[gff@seqnames %in% c("Apico","Mito")]$gene<-substr(gff[gff@seqnames %in% c("Apico","Mito")]$ID,regexpr("PF3D7_",gff[gff@seqnames %in% c("Apico","Mito")]$ID),regexpr("PF3D7_",gff[gff@seqnames %in% c("Apico","Mito")]$ID)+13)
gff$transcript<-gff$ID
gff[gff$type %in% c("protein_coding_gene","pseudogene","ncRNA_gene")]$transcript<-NA
gff[gff$type %in% c("exon","CDS","five_prime_UTR","three_prime_UTR")]$transcript<-gff[gff$type %in% c("exon","CDS","five_prime_UTR","three_prime_UTR")]$Parent
gff$exon<-NA; gff$exon[gff$feature=="exon"]<-gff$ID[gff$feature=="exon"]
gff$source<-NULL; gff$phase<-NULL; gff$Note<-NULL; gff$type<-NULL; gff$protein_source_id<-NULL; gff$score<-NULL

seqlevelsStyle(gff)

# Subsets of GFF
genes<-gff[gff$feature%in% c("protein_coding_gene","pseudogene","ncRNA_gene"),]
genes$transcript<-NULL; genes$Parent<-NULL; genes$description<-NULL; genes$exon<-NULL
genes<-unlist(range(split(genes,genes$ID)),use.names = T)
genes$gene<-names(genes)

TSS<-genes
TSS[TSS@strand=="+"]@ranges=IRanges(start = TSS[TSS@strand=="+"]@ranges@start, width=1)
TSS[TSS@strand=="-"]@ranges=IRanges(start = TSS[TSS@strand=="-"]@ranges@start + TSS[TSS@strand=="-"]@ranges@width, width=1)

CDS<-gff[gff$feature %in% c("CDS"),]
CDS<-unlist(range(split(CDS,CDS$gene)),use.names = T)
CDS$gene<-names(CDS)

fiveUTR<-sort(c(
  GRanges(
    seqnames=genes[genes@strand=="+"& genes$gene %in% CDS$gene]@seqnames,
    strand="+",
    ranges=IRanges(
      start=genes[genes@strand=="+"& genes$gene %in% CDS$gene]@ranges@start,
      end=CDS[CDS@strand=="+"]@ranges@start-1
    ),
    gene=genes[genes@strand=="+" & genes$gene %in% CDS$gene]$gene),
  GRanges(
    seqnames=genes[genes@strand=="-"& genes$gene %in% CDS$gene]@seqnames,strand="-",
    ranges=IRanges(
      start=CDS[CDS@strand=="-"]@ranges@start+CDS[CDS@strand=="-"]@ranges@width,
      end=genes[genes@strand=="-" & genes$gene %in% CDS$gene]@ranges@start + genes[genes@strand=="-" & genes$gene %in% CDS$gene]@ranges@width -1
    ),
    gene=genes[genes@strand=="-" &genes$gene %in% CDS$gene]$gene
  )
)
)
fiveUTR<-renameChrs(renameChrs(fiveUTR[fiveUTR@ranges@width>0]))

fiveUTR$feature<-"five_prime_UTR"


threeUTR<-sort(c(
  GRanges(
    seqnames=genes[genes@strand=="+"& genes$gene %in% CDS$gene]@seqnames,strand="+",
    ranges=IRanges(
      start=CDS[CDS@strand=="+"]@ranges@start + CDS[CDS@strand=="+"]@ranges@width,
      end=genes[genes@strand=="+"& genes$gene %in% CDS$gene]@ranges@start + genes[genes@strand=="+"& genes$gene %in% CDS$gene]@ranges@width-1,
    ),
    gene=genes[genes@strand=="+" & genes$gene %in% CDS$gene]$gene
  )
  ,
  GRanges(
    seqnames=genes[genes@strand=="-"& genes$gene %in% CDS$gene]@seqnames,
    strand="-",
    ranges=IRanges(
      start=genes[genes@strand=="-" & genes$gene %in% CDS$gene]@ranges@start,
      end=CDS[CDS@strand=="-"]@ranges@start-1
    ),
    gene=genes[genes@strand=="-" &genes$gene %in% CDS$gene]$gene
  )
))
threeUTR<-renameChrs(renameChrs(threeUTR[threeUTR@ranges@width>0]))
threeUTR$feature<-"three_prime_UTR"

# start peak calling analysis
if(!dir.exists("./bams/"))dir.create("./bams")
if(!dir.exists("./standard/"))dir.create("./standard/")
if(!dir.exists("./DirectGFPvsTy1/"))dir.create("./DirectGFPvsTy1/")
if(!dir.exists("./reads/"))dir.create("./reads/")

R1<-list.files(path="./reads",pattern="R1.fastq.gz")
R2<-list.files(path="./reads",pattern="R2.fastq.gz")

sample<-ldply(strsplit(R1,"R"))[,1]

#### TRIM READS ####
SysCommands<-paste(
  "trimmomatic PE -threads 4",
  paste0("./reads/",R1),
  paste0("./reads/",R2),
  paste0("./reads/",sample,"R1.paired.fq.gz"),
  paste0("./reads/",sample,"R1.unpaired.fastq.gz"),
  paste0("./reads/",sample,"R2.paired.fq.gz"),
  paste0("./reads/",sample,"R2.unpaired.fastq.gz"),
  "ILLUMINACLIP:/home/bkafsack/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:30 &"
)
#		run SysCommands in terminal
#   system("rm ./reads/*unpaired.fastq.gz",ignore.stdout=T,ignore.stderr=F)

#### ALIGN READS WITH BWA ####
SysCommands<-paste0(
  'bwa mem -M -t 30 -R "@RG\\tID:kafsacklab\\tSM:',
  sample,
  '\\tPL:Illumina\\tLB:lib\\tPU:unit"',
  ' ./genome/PlasmoDB-51_Pfalciparum3D7_Genome.fasta ',
  paste0("./reads/",sample,"R1.paired.fq.gz "),
  paste0("./reads/",sample,"R2.paired.fq.gz "),
  '| samtools view -hbf 2 - | samtools sort -n --threads 30 --reference ./genome/PlasmoDB-51_Pfalciparum3D7_Genome.fasta -O BAM -o ',
  paste0("./bams/",sample,"_sorted.bam"),'\n')

system("rm runme.sh") 
cat(c("#!/bin/bash\nsource ~/.profile\n",SysCommands[c(2,4,6,8)]),file="runme.sh")
system('chmod u+x runme.sh')
# system("./runme.sh&") 

##### SORT BAMS BY NAME ##### ##### 
SysCommands<-paste("samtools sort -n -l 9 --threads 30 -O BAM --reference ./genome/PlasmoDB-51_Pfalciparum3D7_Genome.fasta -o", paste0("./bams/",sample,"sorted.bam"),paste0("./bams/",sample,"bam"))
#	run SysCommands in terminal

#### CALL PEAKS with MACS2 use Ty1 Ab directly as control calling for both replicates combined ############################

# call peaks with FC>2 & q< 0.001
SysCommands<-paste(
  "macs2 callpeak",
  "-t","./bams/Hdp1_GFP_rep1.sorted.bam","./bams/Hdp1-GFP_rep2.sorted.bam", 
  "-c", "./bams/Hdp1_Ty1_rep1.sorted.bam", "./bams/Hdp1-Ty1_rep2.sorted.bam",
  "-f BAMPE -B -g 2.3e7 -q 0.001 --cutoff-analysis --fe-cutoff 2 --call-summits --outdir ./DirectGFPvsTy1 -n", 
  "Hdp1_GFPvsTy1_both_2FC_q0.001",
  "&") 

#		run SysCommands in terminal

SysCommands<-paste(
  "macs2 bdgcmp",
  "-t","./DirectGFPvsTy1/Hdp1_GFPvsTy1_both_2FC_q0.001_treat_pileup.bdg",
  "-c","./DirectGFPvsTy1/Hdp1_GFPvsTy1_both_2FC_q0.001_control_lambda.bdg",
  "-m FE",
  "-o ./DirectGFPvsTy1/Hdp1_GFPvsTy1_both_FE.bdg",
  "&")
#		run SysCommands in terminal
SysCommands<-paste(
  "macs2 bdgcmp",
  "-t","./DirectGFPvsTy1/Hdp1_GFPvsTy1_both_2FC_q0.001_treat_pileup.bdg",
  "-c","./DirectGFPvsTy1/Hdp1_GFPvsTy1_both_2FC_q0.001_control_lambda.bdg",
  "-m qpois",
  "-o ./DirectGFPvsTy1/Hdp1_GFPvsTy1_both_qpois.bdg",
  "&")

#		run SysCommands in terminal

########################################################
#### Analyzing Called Peaks ####
########################################################

summits<-read.delim("DirectGFPvsTy1/Hdp1_GFPvsTy1_both_2FC_q0.001_peaks.narrowPeak",header = F)
colnames(summits)<-c("chr","peak.start","peak.end","PeakID","score","strand","signal","pval","fdr","start")
summits$PeakID<-sub("Hdp1_GFPvsTy1_both_2FC_q0.001_","",summits$PeakID)
summits$start<-summits$peak.start+summits$start
summits$strand<-"*"
summits$end<-summits$start
summits$summitID<-summits$PeakID
summits$PeakID[!grepl("\\d",substring(summits$PeakID,nchar(summits$PeakID),nchar(summits$PeakID)))]<-   substring(summits$PeakID[!grepl("\\d",substring(summits$PeakID,nchar(summits$PeakID),nchar(summits$PeakID)))],1,nchar(summits$PeakID[!grepl("\\d",substring(summits$PeakID,nchar(summits$PeakID),nchar(summits$PeakID)))])-1)
summits<-summits[,c("chr","start","end","strand","summitID","PeakID","score","signal","pval","fdr","peak.start","peak.end")]
summits$summitID<-sub("peak","summit",summits$summitID)

peaks<-summits[,c("chr","peak.start","peak.end","summitID","PeakID","score","strand","signal","pval","fdr","start")]
peaks$start<-peaks$start-peaks$peak.start
peaks$strand="."
write.table(peaks,file="./DirectGFPvsTy1/peaks.narrowPeak",quote=F,row.names = F,col.names = F,sep="\t")
write.table(summits[,c("chr","start","end","PeakID","score","signal","fdr")],file="./DirectGFPvsTy1/summits.bed",quote=F,row.names = F,col.names = F,sep="\t")
summits<-makeGRangesFromDataFrame(summits,keep.extra.columns = T)
summits<-renameChrs(summits)

summits.genic<-annoPeaks(summits, annoData=annoGR(genes), bindingType="fullRange", bindingRegion=c(-2001,1), ignore.peak.strand=T, select="all")
summits.promoter<-annoPeaks(summits, annoData=annoGR(genes), bindingType="startSite", bindingRegion=c(-2001,1), ignore.peak.strand=T, select="all")
summits.promoter$insideFeature<-"promoter"
length(unique(summits.promoter$summitID)) # 904/1188 76.1%
summits.5UTR<-annoPeaks(summits, annoData=annoGR(fiveUTR), bindingType="fullRange", bindingRegion=c(-1,1), ignore.peak.strand=T, select="all")
summits.5UTR$insideFeature<-"5'UTR"
length(unique(summits.5UTR$summitID)) # 233/1188  19.6%
summits.CDS<-annoPeaks(summits, annoData=annoGR(CDS), bindingType="fullRange", bindingRegion=c(-1,1), ignore.peak.strand=T, select="all")
summits.CDS$insideFeature<-"CDS"
length(unique(summits.CDS$summitID)) # 171/1188  14.4%
summits.3UTR<-annoPeaks(summits, annoData=annoGR(threeUTR), bindingType="fullRange", bindingRegion=c(-1,1), ignore.peak.strand=T, select="all")
summits.3UTR$insideFeature<-"3'UTR"
length(unique(summits.3UTR$summitID)) # 123/1188 (10.4)

summits.upstream<-c(summits.5UTR,summits.promoter)
names(summits.upstream)<-paste0("X",1:length(summits.upstream))
length(unique(summits.upstream$summitID)) # 1011/1188 85.1%

summits.intergenic<-summits[!(summits$summitID %in% c(summits.genic$summitID,summits.promoter$summitID)),]
length(summits.intergenic)

summits$gene<-"intergenic"
temp<-split(summits.upstream,summits.upstream$summitID)
temp2<-unlist(lapply(temp, function(gr)return(paste(paste(gr$gene,gr$insideFeature),collapse=";"))))
summits$gene[match(names(temp2),summits$summitID)]<-paste0(temp2,";")

write.table(data.frame(renameChrs(summits))[,c("summitID","seqnames","start","score","signal","fdr","PeakID","peak.start","peak.end","gene")],file="./DirectGFPvsTy1/summits.tab",quote=F,row.names = F,col.names = T,sep="\t")

# Plot summits with respect to gene regions 
reldists.5UTR<-c(
  (summits.5UTR[summits.5UTR$feature.strand=="+"]@ranges@start-summits.5UTR[summits.5UTR$feature.strand=="+"]$feature.ranges@start)/summits.5UTR[summits.5UTR$feature.strand=="+"]$feature.ranges@width,
  (summits.5UTR[summits.5UTR$feature.strand=="-"]$feature.ranges@start+summits.5UTR[summits.5UTR$feature.strand=="-"]$feature.ranges@width-summits.5UTR[summits.5UTR$feature.strand=="-"]@ranges@start)/summits.5UTR[summits.5UTR$feature.strand=="-"]$feature.ranges@width
)
reldists.3UTR<-c((summits.3UTR[summits.3UTR$feature.strand=="+"]@ranges@start-summits.3UTR[summits.3UTR$feature.strand=="+"]$feature.ranges@start)/summits.3UTR[summits.3UTR$feature.strand=="+"]$feature.ranges@width,(summits.3UTR[summits.3UTR$feature.strand=="-"]$feature.ranges@start+summits.3UTR[summits.3UTR$feature.strand=="-"]$feature.ranges@width-summits.3UTR[summits.3UTR$feature.strand=="-"]@ranges@start)/summits.3UTR[summits.3UTR$feature.strand=="-"]$feature.ranges@width)
reldists.CDS<-c((summits.CDS[summits.CDS$feature.strand=="+"]@ranges@start-summits.CDS[summits.CDS$feature.strand=="+"]$feature.ranges@start)/summits.CDS[summits.CDS$feature.strand=="+"]$feature.ranges@width,(summits.CDS[summits.CDS$feature.strand=="-"]$feature.ranges@start+summits.CDS[summits.CDS$feature.strand=="-"]$feature.ranges@width-summits.CDS[summits.CDS$feature.strand=="-"]@ranges@start)/summits.CDS[summits.CDS$feature.strand=="-"]$feature.ranges@width)
reldists.upstream<-summits.promoter$distance*-1/2000

reldists=data.frame(
  region=rep(c("upstream","5'UTR","CDS","3'UTR"),times=c(length(reldists.upstream),length(reldists.5UTR),length(reldists.CDS),length(reldists.3UTR))),
  dist=c(
    reldists.upstream*2000,
    reldists.5UTR*median(fiveUTR@ranges@width),
    reldists.CDS*median(CDS@ranges@width)+median(fiveUTR@ranges@width),
    reldists.3UTR*median(threeUTR@ranges@width)+median(fiveUTR@ranges@width)+median(CDS@ranges@width)
  )
)

median(fiveUTR@ranges@width)
median(CDS[seqnames(CDS) %in% seqlevels(CDS)[3:16]]@ranges@width)
median(threeUTR@ranges@width)

ggplot(data=reldists,aes(x=dist))+
  geom_density(adjust= 0.5,color="black",size=1,trim=F)+
  geom_area(aes(x = stage(dist, after_stat = oob_censor(x, c(-2002, 5)))),stat = "density", adjust= 0.5, fill="purple")+
  geom_area(aes(x = stage(dist, after_stat = oob_censor(x, c(-5, 510)))),stat = "density", adjust= 0.5, fill="darkgreen")+
  geom_area(aes(x = stage(dist, after_stat = oob_censor(x, c(500, 2190)))),stat = "density", adjust= 0.5, fill="darkorange")+
  geom_area(aes(x = stage(dist, after_stat = oob_censor(x, c(2179, 2587)))),stat = "density", adjust= 0.5, fill="darkred")+
  scale_x_continuous(breaks=seq(from=-2000,to=2500,by=500),labels=NULL)+
  scale_y_continuous(breaks=seq(from=0,to=8e-4,by=1e-4),labels=NULL)+
  theme_classic()+
  theme(axis.text=element_text(colour = "black"))+
  labs(y="",x="")

ggsave("summit.positions.pdf",width=6,height=3,unit="in")

#### #### #### #### #### #### #### #### #### #### #### 
############  CALLING MOTIFS NEAR SUMMITS
#### #### #### #### #### #### #### #### #### #### #### 

para.summits<-renameChrs(flank(summits,50,both = T))
names(para.summits)<-para.summits$summitID
summary(para.summits$signal)
saveme<-data.frame(para.summits)
saveme<-saveme[,c("seqnames","start","end","summitID","score","strand","signal","pval","fdr")]
saveme$strand="."
saveme$dist=51
write.table(saveme,"./DirectGFPvsTy1/100bpFlankingSummits.narrowPeak",sep="\t",quote=F,row.names = F,col.names = F)

SysCommands<-paste(
  "findMotifsGenome.pl",
  "./DirectGFPvsTy1/100bpFlankingSummits.narrowPeak",
  "./genome/PlasmoDB-51_Pfalciparum3D7_Genome.fasta",
  "./Hdp1_GFPvsTy1_both_2FC_q0.001_summit_100bp_len5-12",
  "-size given -len 5,6,7,8,9,10,11,12 -p 32 -S 5",
  "&")

########### Read in motifs

homer.motifA<-read.table(textConnection(readLines("Hdp1_GFPvsTy1_both_2FC_q0.001_summit_100bp_len5-12/homerResults/motif1.motif")[-1]), header = F, sep="\t")
colnames(homer.motifA)<-c("A","C","G","T")
homer.motifB<-read.table(textConnection(readLines("Hdp1_GFPvsTy1_both_2FC_q0.001_summit_100bp_len5-12/homerResults/motif2.motif")[-1]), header = F, sep="\t")
colnames(homer.motifB)<-c("A","C","G","T")

hdp1.mtxA <- PFMatrix(ID="HDP1_motif1", name="HDP1_motif1", 
                      matrixClass="homeo", strand="+",
                      bg=c(A=0.41, C=0.09, G=0.09, T=0.41),
                      profileMatrix=t(homer.motifA*1000)
)

hdp1.mtxB <- PFMatrix(ID="HDP1_motif2", name="HDP1_motif2", 
                      matrixClass="homeo", strand="+",
                      bg=c(A=0.41, C=0.09, G=0.09, T=0.41),
                      profileMatrix=t(homer.motifB*1000)
)
# 
peaks<- summits[!duplicated(summits$PeakID)]
peaks@ranges <-IRanges(start=peaks$peak.start,end=peaks$peak.end)
peaks$summitID<-NULL;peaks$peak.start<-NULL;peaks$peak.end<-NULL;
peaks2<-renameChrs(peaks)

#### HOMER MOTIF A
motifAinPeaks <- matchMotifs(hdp1.mtxA, peaks2, genome = BSGenome.Pf3D7.PlasmoDB.51,out = "positions",bg="genome",p.cutoff=1e-4)[[1]]
motifAinPeaks$motifID<-paste0("motifA.",1:length(motifAinPeaks))
seqinfo(motifAinPeaks)<-seqinfo(peaks2)

peaksWithMotifA<-subsetByOverlaps(peaks2,motifAinPeaks)
write.table(data.frame(peaksWithMotifA),file="./DirectGFPvsTy1/peaksWithMotifA.tab",quote=F,row.names = F,col.names = T,sep="\t")

fun<-function(gr1) gr1$PeakID<-subsetByOverlaps(peaksWithMotifA,gr1)$PeakID #adds peakID for each motif

tmp<-sapply(split(motifAinPeaks,motifAinPeaks$motifID),fun)
motifAinPeaks$PeakID[match(names(tmp),motifAinPeaks$motifID)]<-tmp

#### HOMER MOTIF B
motifBinPeaks <- matchMotifs(hdp1.mtxB, peaks2, genome = BSGenome.Pf3D7.PlasmoDB.51,out = "positions",bg="genome",p.cutoff=1e-4)[[1]]
motifBinPeaks$motifID<-paste0("motifB.",1:length(motifBinPeaks))
seqinfo(motifBinPeaks)<-seqinfo(peaks2)

peaksWithMotifB<-subsetByOverlaps(peaks2,motifBinPeaks)
write.table(data.frame(peaksWithMotifB),file="./DirectGFPvsTy1/peaksWithmotifB.tab",quote=F,row.names = F,col.names = T,sep="\t")

fun<-function(gr1) gr1$PeakID<-subsetByOverlaps(peaksWithMotifB,gr1)$PeakID #adds peakID for each motif

tmp<-sapply(split(motifBinPeaks,motifBinPeaks$motifID),fun)
motifBinPeaks$PeakID[match(names(tmp),motifBinPeaks$motifID)]<-tmp

# 1e-4
length(peaksWithMotifA) # 615/1003 (61.3%)
length(peaksWithMotifB) # 478/1003 (47.6%)
length(unique(c(motifAinPeaks$PeakID,motifBinPeaks$PeakID)))  # 784/1003 (78.1%) has at least one motif
length(intersect(motifAinPeaks$PeakID,motifBinPeaks$PeakID)) # 309/1003 (30.8%) have both
length(setdiff(motifAinPeaks$PeakID,motifBinPeaks$PeakID))   # 306 have motifA only
length(setdiff(motifBinPeaks$PeakID,motifAinPeaks$PeakID))   # 169 have motifB only

######## Distance of Motifs to nearest summit
fun2<-function(mip){  #find closest summit to motif
  smts<-summits[summits$PeakID==mip$PeakID]
  mip$closestSummmit<-smts[abs(mip@ranges@start+3-smts@ranges@start)==min(abs(mip@ranges@start+3-smts@ranges@start))]$summitID
  mip$distToClosestSummit<-mip@ranges@start+3-smts[abs(mip@ranges@start+3-smts@ranges@start)==min(abs(mip@ranges@start+3-smts@ranges@start))]@ranges@start
  return(mip)
}

tmp<-lapply(split(motifAinPeaks,motifAinPeaks$motifID),fun2)
motifAinPeaks<-unlist(GRangesList(tmp))
upstreamPeaksWithMotifA<-peaksWithMotifA[peaksWithMotifA$PeakID %in% summits.upstream$PeakID]

mean(motifAinPeaks$distToClosestSummit[abs(motifAinPeaks$distToClosestSummit)<=400])
sd(motifAinPeaks$distToClosestSummit[abs(motifAinPeaks$distToClosestSummit)<=400])/sqrt(length(motifAinPeaks$distToClosestSummit[abs(motifAinPeaks$distToClosestSummit)<=400]))

ggplot(data.frame(motifAinPeaks),aes(x=distToClosestSummit))+
  geom_density(adjust= 0.5,color="black",size=1,trim=F,fill="black")+
  theme_classic()+
  scale_x_continuous(breaks=seq(-400,400,100),limits=c(-400,400))+
  scale_y_continuous(breaks=NULL,labels=NULL,guide = NULL)+
  theme(axis.text=element_text(colour = "black"))+
  labs(y="",x="",title="motif A")
ggsave("motifADist2Summit.pdf",width=6,height=2)

tmp<-lapply(split(motifBinPeaks,motifBinPeaks$motifID),fun2)
motifBinPeaks<-unlist(GRangesList(tmp))
upstreamPeaksWithMotifB<-peaksWithMotifB[peaksWithMotifB$PeakID %in% summits.upstream$PeakID]

mean(motifBinPeaks$distToClosestSummit[abs(motifBinPeaks$distToClosestSummit)<=400])
sd(motifBinPeaks$distToClosestSummit[abs(motifBinPeaks$distToClosestSummit)<=400])/sqrt(length(motifBinPeaks$distToClosestSummit[abs(motifBinPeaks$distToClosestSummit)<=400]))

ggplot(data.frame(motifBinPeaks),aes(x=distToClosestSummit))+
  geom_density(adjust= 0.5,color="black",size=1,trim=F,fill="black")+
  theme_classic()+
  scale_x_continuous(breaks=seq(-400,400,100),limits=c(-400,400))+
  scale_y_continuous(breaks=NULL,labels=NULL,guide = NULL)+
  theme(axis.text=element_text(colour = "black"))+
  labs(y="",x="",title="motif B")
ggsave("motifBDist2Summit.pdf",width=6,height=2)

write.table(data.frame(motifBinPeaks),file="./DirectGFPvsTy1/MotifBInPeaks.tab",quote=F,row.names = F,col.names = T,sep="\t")

length(upstreamPeaksWithMotifA$PeakID) #526/808 (65.1%) at 1e-4
length(upstreamPeaksWithMotifB$PeakID) #399/808 (49.4%) at 1e-4
length(unique(c(upstreamPeaksWithMotifA$PeakID,upstreamPeaksWithMotifB$PeakID))) #669/808 (82.8%) at 1e-4

####### ####### ####### ####### ####### ####### ####### 
####### Peak/Motif co-occurrence with DGE
####### ####### ####### ####### ####### ####### ####### 
imc.down<-c("PF3D7_1430800","PF3D7_1431100","PF3D7_1345600","PF3D7_1323700","PF3D7_1221400","PF3D7_0109000","PF3D7_0717600","PF3D7_1011000","PF3D7_1222700","PF3D7_0621400","PF3D7_1406800","PF3D7_0918000","PF3D7_0423500")
genesWithUpstreamPeaks<-unique(summits.upstream$gene) 
length(imc.down[imc.down %in% genesWithUpstreamPeaks]) # 11/13

DGE<-read.delim("all_DGE.tab")

DGE$upstream.peaks<-""

temp<-unlist(lapply(split(summits.upstream[summits.upstream$gene %in% DGE$gene.ID],summits.upstream[summits.upstream$gene %in% DGE$gene.ID]$gene), function(gr) return(paste0(gr$summitID,collapse=";"))))
DGE$upstream.peaks[match(names(temp),DGE$gene.ID)]<-temp

tmp<-peaks[peaks$PeakID %in% summits.upstream$PeakID & peaks$PeakID %in% upstreamPeaksWithMotifA$PeakID]$gene
tmp<-gsub(" inside;","",tmp)
tmp<-gsub(" upstream;","",tmp)
tmp.ids<-substr(tmp,1,13)
tmp<-sub("PF3D7_\\d+","",tmp,perl=T)
tmp.ids<-c(tmp.ids, substr(tmp,1,13))
tmp<-sub("PF3D7_\\d+","",tmp,perl=T)
tmp.ids<-c(tmp.ids, substr(tmp,1,13))

genesWithMotifAInUpstreamPeaks<-tmp.ids[nchar(tmp.ids)>0]
DGE$motifAInPeak<-DGE$gene.ID %in% genesWithMotifAInUpstreamPeaks

tmp<-peaks[peaks$PeakID %in% summits.upstream$PeakID & peaks$PeakID %in% upstreamPeaksWithMotifB$PeakID]$gene
tmp<-gsub(" inside;","",tmp)
tmp<-gsub(" upstream;","",tmp)
tmp.ids<-substr(tmp,1,13)
tmp<-sub("PF3D7_\\d+","",tmp,perl=T)
tmp.ids<-c(tmp.ids, substr(tmp,1,13))
tmp<-sub("PF3D7_\\d+","",tmp,perl=T)
tmp.ids<-c(tmp.ids, substr(tmp,1,13))

genesWithMotifBInUpstreamPeaks<-tmp.ids[nchar(tmp.ids)>0]
DGE$motifBInPeak<-DGE$gene.ID %in% genesWithMotifBInUpstreamPeaks

# calculate distance of TSS to nearest summit

fun3<-function(smt.us,tss=TSS[TSS$gene %in% genesWithUpstreamPeaks]){  #find closest summit to motif
  if(as.logical(tss[tss$gene==smt.us$gene]@strand=="+")) smt.us$dist2tss<-smt.us@ranges@start-tss[tss$gene==smt.us$gene]@ranges@start+1
  if(as.logical(tss[tss$gene==smt.us$gene]@strand=="-")) smt.us$dist2tss<-tss[tss$gene==smt.us$gene]@ranges@start-smt.us@ranges@start
  return(smt.us)
} 
summits.upstream<-unlist(GRangesList(lapply(split(summits.upstream,names(summits.upstream)),fun3)))
temp<-ddply(data.frame(summits.upstream),.(gene),summarize,min=min(dist2tss))

DGE$DistOfnearestSummit2TSS<-temp$min[match(DGE$gene.ID,temp$gene)]

write.table(DGE,file="all_DGE+summits.tab",sep="\t",quote = F,row.names = F,col.names = T)
write.table(DGE[DGE$q_value<=0.05 & abs(DGE$log2.KO.WT.)>=0.5 & (DGE$fpkm_KO>=10|DGE$fpkm_WT>=10),],file="sign_DGE+summits.tab",sep="\t",quote = F,row.names = F,col.names = T)

DOWN <- DGE[DGE$q_value<=0.05 & round(DGE$log2.KO.WT.,2)<=-0.5 & (DGE$fpkm_KO>=10|DGE$fpkm_WT>=10),]
UP <- DGE[DGE$q_value<=0.05 & round(DGE$log2.KO.WT.,2)>=0.5 & (DGE$fpkm_KO>=10|DGE$fpkm_WT>=10),]

nrow(DOWN[DOWN$upstream.peaks!="",])/nrow(DOWN) # 86/155 (55.4%) CORRECT
nrow(UP[UP$upstream.peaks!="",])/nrow(UP) # 18/103 (17.4%)  CORRECT

plotme<-DGE
plotme$peakUpstream<-factor("all",levels=c("all","upstream","upstream+motif"))
plotme$peakUpstream[plotme$gene.ID %in% genesWithUpstreamPeaks] <-"upstream"
plotme$fpkm_WT[plotme$fpkm_WT>10000]=10000
plotme$fpkm_WT[plotme$fpkm_WT<1]=1
plotme<-plotme[order(plotme$fpkm_WT),]

ggplot(plotme,aes(x=log10(fpkm_WT),fill=peakUpstream))+
  geom_histogram(bins=35)+
  theme_classic()+
  scale_x_continuous(labels=NULL)+# labels=c("<=1",10,100,1000,">=10000"))+
  scale_y_continuous(breaks=seq(0,500,100),labels=NULL)+# labels=c("<=1",10,100,1000,">=10000"))+ 
  labs(x="", y="", fill="")#+theme(legend.position = "NULL")
ggsave("upstreamMotifs_vsFPKM.pdf",width=6,height=4,unit="in")

wilcox.test(1:5607,which(plotme$peakUpstream=="upstream"))
mean(plotme$fpkm_WT[plotme$peakUpstream=="upstream"])
mean(plotme$fpkm_WT)
mean(log(plotme$fpkm_WT[plotme$peakUpstream=="upstream"]))
mean(log(plotme$fpkm_WT))
median(plotme$fpkm_WT[plotme$peakUpstream=="upstream"])
median(plotme$fpkm_WT)


#UP<-DGE$gene.ID[DGE$q_value<=0.05 & DGE$log2.KO.WT.>=0.5 & (DGE$fpkm_WT >=10 | DGE$fpkm_KO >=10)]
#DOWN<-DGE$gene.ID[DGE$q_value<=0.05 & DGE$log2.KO.WT.<=0.5 & (DGE$fpkm_WT >=10 | DGE$fpkm_KO >=10)]

DOWN.pos<-c(promoter=sum(DOWN$gene.ID %in% summits.promoter$gene[summits.promoter$gene %in% DOWN$gene.ID]),
            fiveUTR=sum(DOWN$gene.ID %in% summits.5UTR$gene[summits.5UTR$gene %in% DOWN$gene.ID]),
            CDS=sum(DOWN$gene.ID %in% summits.CDS$gene[summits.CDS$gene %in% DOWN$gene.ID]),
            threeUTR=sum(DOWN$gene.ID %in% summits.3UTR$gene[summits.3UTR$gene %in% DOWN$gene.ID])
)
DOWN.pos["no peak"]<-length(DOWN$gene.ID)-sum(DOWN.pos)
round(DOWN.pos/sum(DOWN.pos),3)

UP.pos<-c(promoter=sum(UP$gene.ID %in% summits.promoter$gene[summits.promoter$gene %in% UP$gene.ID]),
          fiveUTR=sum(UP$gene.ID %in% summits.5UTR$gene[summits.5UTR$gene %in% UP$gene.ID]),
          CDS=sum(UP$gene.ID %in% summits.CDS$gene[summits.CDS$gene %in% UP$gene.ID]),
          threeUTR=sum(UP$gene.ID %in% summits.3UTR$gene[summits.3UTR$gene %in% UP$gene.ID])
)
UP.pos["no peak"]<-length(UP$gene.ID)-sum(UP.pos)
round(UP.pos/sum(UP.pos),3)

plotme<-melt(cbind(UP.pos,DOWN.pos))
plotme$Var2<-rep(c("up","down"),each=5)
plotme$pct<-NA
plotme$pct[1:5]<-plotme$value[1:5]/sum(plotme$value[1:5])
plotme$pct[6:10]<-plotme$value[6:10]/sum(plotme$value[6:10])

ggplot(plotme,aes(x="",y=pct, fill=Var1))+
  geom_bar(stat="identity",color="black")+
  coord_polar("y", start=0)+facet_wrap(~Var2)+
  theme_classic()+scale_y_continuous(breaks=NULL,name="")+scale_x_discrete(breaks=NULL,name=NULL)+
  scale_fill_manual(values=c("purple","darkgreen","darkorange","darkred","grey80"),name=NULL)+
  theme(legend.position = "bottom")

# what fraction of DOWN/UP regulated genes have motifs upstream?
summitsInPeaksWithMotifA<-summits[summits$PeakID %in% motifAinPeaks$PeakID]
summitsInPeaksWithMotifB<-summits[summits$PeakID %in% motifBinPeaks$PeakID]

length(DOWN[DOWN %in% genesWithUpstreamPeaks])/length(DOWN) #80/156 (51.2%) of DOWN genes have peaks upstream
length(DOWN[DOWN %in% c(genesWithMotifAInUpstreamPeaks,genesWithMotifBInUpstreamPeaks)])/length(DOWN[DOWN %in% genesWithUpstreamPeaks]) #74/80 (92.5%) of DOWN genes with peaks upstream have motifs 

sort(DOWN[DOWN %in% genesWithUpstreamPeaks])
sort(DOWN[DOWN %in% genesWithMotifAInUpstreamPeaks])
sort(DOWN[DOWN %in% genesWithMotifBInUpstreamPeaks])

length(UP[UP %in% genesWithUpstreamPeaks])/length(UP) # 15/103 14.6% of UP genes have peaks & motifs upstream
length(UP[UP %in% c(genesWithMotifAInUpstreamPeaks,genesWithMotifBInUpstreamPeaks)])/length(UP[UP %in% genesWithUpstreamPeaks]) # 12/15 (80%) of UP genes with peaks upstream have motifs
sort(UP[UP %in% genesWithUpstreamPeaks])
sort(UP[UP %in% genesWithMotifAInUpstreamPeaks])
sort(UP[UP %in% genesWithMotifBInUpstreamPeaks])

##################################
write.table(data.frame(motifAinPeaks)[,c("seqnames","start","end","strand","score")],file="./DirectGFPvsTy1/MotifAInPeaks.bed",quote=F,row.names = F,col.names = F,sep="\t")
write.table(data.frame(motifBinPeaks)[,c("seqnames","start","end","strand","score")],file="./DirectGFPvsTy1/MotifBInPeaks.bed",quote=F,row.names = F,col.names = F,sep="\t")
# Plot Gviz Tracks

qpois<-data.frame(import.bedGraph("DirectGFPvsTy1/Hdp1_GFPvsTy1_both_qpois.bdg"))
qpois$seqnames<-as.character(paste0("Chr",substring(qpois$seqnames,7,8)))
qpois<-qpois[qpois$seqnames %in% seqlevels(gff)[3:16],]
qpois<-sort(makeGRangesFromDataFrame(qpois,keep.extra.columns = T))
seqlengths(qpois)<-seqlengths(gff)[3:16]
isCircular(qpois)<-isCircular(gff)[3:16]
genome(qpois)<-"51"
seqinfo(qpois)

FE<-data.frame(import.bedGraph("DirectGFPvsTy1/Hdp1_GFPvsTy1_both_FE.bdg"))
FE$seqnames<-as.character(paste0("Chr",substring(FE$seqnames,7,8)))
FE<-FE[FE$seqnames %in% seqlevels(gff)[3:16],]
FE<-sort(makeGRangesFromDataFrame(FE,keep.extra.columns = T))
seqlengths(FE)<-seqlengths(gff)[3:16]
isCircular(FE)<-isCircular(gff)[3:16]
genome(FE)<-"51"
seqinfo(FE)
FE$log2FE<-log(FE$score,2)

saveme<-data.frame(summits)

saveme$seqnames<-paste0("PF3D7_",substr(saveme$seqnames,4,5),"_v3")
saveme<-saveme[,c(1,2,6,8,10:12,13)]
colnames(saveme)<-c("chr","pos","summit.id","FE","fdr","peak.start","peak.end","gene.associations")
write.table(saveme,file="./summits.tab",sep="\t",col.names = T,row.names = F,quote = F)

temp<-peaks[grepl("\\d\\D",peaks$PeakID,perl=TRUE)]
temp$PeakID<-substr(temp$PeakID,1,nchar(temp$PeakID)-1)
peaks<-sort(c(peaks[!grepl("\\d\\D",peaks$PeakID,perl=TRUE)],temp[!duplicated(temp$PeakID)]))

plotTrax<-function(roi,file.name=NA){
  
  # Genome Axis Track
  axis.track <- GenomeAxisTrack(
    range=roi,
    col="black",
    name=as.character(seqnames(roi)[1]),
    showTitle=T,rotation.title=0,
    fill.range="transparent",
    fontcolor="black",
    labelPos="below",
    cex.title=0.75,
    col.border.title="white"
  )
  #displayPars(axis.track)
  
  # Data Tracks
  qpois.track<-DataTrack(subsetByOverlaps(qpois,roi), name="-log10 FDR", col="black",col.histogram="black",cex.title=0.75)
  FE.track<-DataTrack(subsetByOverlaps(FE,roi), name="log2 GFP/Ty1 FE", col="black",col.histogram="black",cex.title=0.75)
  
  # Combined Peaks & summits Track
  peaks.track<-AnnotationTrack(
    subsetByOverlaps(peaks,roi),
    name="peaks",
    id=subsetByOverlaps(peaks,roi)$PeakID,
    fill="transparent",
    col="black",
    rotation.title=0,
    just.group = "above",
    stacking="dense",
    cex=.75,
    cex.title=0.75
  )
  
  summits.track<-AnnotationTrack(
    subsetByOverlaps(summits,roi),
    id=" ",
    name="peaks",
    fill="black",
    col="black",
    rotation.title=0,
    just.group = "above",
    stacking="dense",
    cex=.75,
    cex.title=0.75
  )
  
  peaks.summits.track<-OverlayTrack(trackList=list(summits.track,peaks.track))
  
  motifA.track<-AnnotationTrack(
    subsetByOverlaps(renameChrs(motifAinPeaks),roi),
    name="motifs",
    id=" ",
    fill="red",
    col="transparent",
    rotation.title=0,
    cex=.75,
    just.group = "below",
    stacking = "dense",
    cex.title=0.75
  )
  
  motifB.track<-AnnotationTrack(
    subsetByOverlaps(renameChrs(motifBinPeaks),roi),
    name="motifs",
    id=" ",
    fill="blue",
    col="transparent",
    rotation.title=0,
    cex=.75,
    just.group = "below",
    stacking = "dense",
    cex.title=0.75
  )
  
  motifAB.track<-AnnotationTrack(
    subsetByOverlaps(renameChrs(motifAinPeaks),subsetByOverlaps(renameChrs(motifBinPeaks),roi)),
    name="motifs",
    id=" ",
    fill="purple",
    col="transparent",
    rotation.title=0,
    cex=.75,
    just.group = "below",
    stacking = "dense",
    cex.title=0.75
  )
  motif.track<-OverlayTrack(trackList=list(motifA.track,motifB.track,motifAB.track))
  
  # GeneRegionTrack
  plotme<-data.frame(sort(subsetByOverlaps(gff[gff$feature%in% c("CDS","ncRNA"),],roi)))
  plotme$STRAND<-"plus"
  plotme$STRAND[plotme$strand=="-"]<-"minus"
  # 
  CDS.track<-GeneRegionTrack(
    plotme,
    chr=as.character(seqnames(roi)[1]),
    stacking = "pack",
    feature=plotme$STRAND,
    col="black",
    name="genes",
    rotation.title=0,
    transcriptAnnotation = "gene",
    plus="blue",
    minus="red",
    cex.group=0.75,
    cex.title=0.75,
    just.group="left"
  )
  
  # Plot Tracks
  if(!is.na(file.name)) pdf(file.name,width = 8 ,height=4)
  
  plotTracks(
    c(qpois.track,peaks.summits.track,motif.track,CDS.track,axis.track),
    sizes = c(1,0.3,0.3,0.5,0.45)*1.38,
    type="histogram",
    background.panel = "white",
    background.title = "black",
    panel.only=F,
    groupAnnotation="id",
    #just.group = "below",
    labelPos="below",
    col.line="black",
    fontcolor.group="black",
    cex.id=0.7,
    cex.group=0.7,
    cex=1
  )
  
  if(!is.na(file.name)) dev.off()
}

Phil1<-GRanges(seqnames="Chr01",range=IRanges(start=357500,end=369492))
plotTrax(Phil1,'./peaks/Phil1_PF3D7_0109000.pdf')

gap45<-GRanges(seqnames="Chr12",range=IRanges(start=907000,end=922500))
plotTrax(gap45,'./peaks/gap45_PF3D7_1222700.pdf')

PF3D7_1345600<-GRanges(seqnames="Chr13",range=IRanges(start=1815596,end=1830119))
plotTrax(PF3D7_1345600,'./peaks/PF3D7_1345600.pdf')

AP2<-GRanges(seqnames="Chr08",range=IRanges(start=144673,end=171353))
plotTrax(AP2,'./peaks/AP2_PF3D7_0802100.pdf')

AP2.O2<-GRanges(seqnames="Chr05",range=IRanges(start=699300,end=716800))
plotTrax(AP2.O2,'./peaks/AP2-O2_PF3D7_0516800.pdf')

pip3<-GRanges(seqnames="Chr14",range=IRanges(start=1210500,end=1217500))
plotTrax(pip3,'./peaks/pip3_PF3D7_1430800.pdf')

gapm1<-GRanges(seqnames="Chr13",range=IRanges(start=978572,end=989034))
plotTrax(gapm1,'./peaks/gapm1_PF3D7_1323700.pdf')

IMC1h<-GRanges(seqnames="Chr12",range=IRanges(start=849498,end=863100))
plotTrax(IMC1h,'./peaks/IMC1h_PF3D7_1221400.pdf')

gapm3<-GRanges(seqnames="Chr14",range=IRanges(start=246342,end=257342))
plotTrax(gapm3,'./peaks/gapm3_PF3D7_1406800.pdf')

gap50<-GRanges(seqnames="Chr09",range=IRanges(start=736484,end=750644))
plotTrax(gap50,'./peaks/gap50_PF3D7_0918000.pdf')

gapm2<-GRanges(seqnames="Chr04",range=IRanges(start=1051500,end=1067000))
plotTrax(gapm2,'./peaks/gapm2_PF3D7_0423500.pdf')

mdv1<-GRanges(seqnames="Chr12",range=IRanges(start=652725,end=666382))
plotTrax(mdv1,'./peaks/mdv1_PF3D7_1216500.pdf')

rnf1<-GRanges(seqnames="Chr03",range=IRanges(start=586119,end=607815))
plotTrax(rnf1,'./peaks/rnf1_PF3D7_0314700.pdf')

hdp1<-GRanges(seqnames="Chr14",range=IRanges(start=2694014,end=2718476))
plotTrax(hdp1,'./peaks/hdp1_PF3D7_1466200.pdf')

