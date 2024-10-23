
#############################Fig.S18_Traits density of differential OTUs##########################

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(splitstackshape)

### 1 Bacteria differential OTUs
rm(list = setdiff(ls(), c("p_amf","p_bac1")))
setwd("E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Bacteria/otu")

# read group information
group = read.csv("sample_metadata.csv", header=T, row.names=1) 

group_MC <- group %>% filter(country == "Uzbekistan",plants=="Maize"|plants=="Cotton")
group_MC$type <- rep("MC",nrow(group_MC))

group_TS <- group %>% filter(country == "Uzbekistan",plants=="Tamarix chinensis"|plants=="Suaeda")
group_TS$type <- rep("TS",nrow(group_TS))

group2 <- rbind(group_MC,group_TS)
group2$type <- as.factor(group2$type)

# read otu table
otu_table = read.csv("otu_bacteriaFlattening.csv", header=T, row.names=1, check.names=F)
otu_table <- otu_table[,rownames(group2)]
otu_table <- otu_table[rowSums(otu_table) > 0,]

# filter otu
otu_relative <- apply(otu_table, 2, function(x){x*100/sum(x)})
threshold = 0.0001
idx <- rowSums(otu_relative > threshold) >= 1
otu <- as.data.frame(otu_table[idx, ])
otu_relative <- as.data.frame(otu_relative[idx, ])

# read taxonomy
taxonomy = read.csv("HMME.bact.OTU.ID1.csv", row.names= 1,header=T)
taxonomy <- taxonomy[match(rownames(otu),taxonomy$OTU.ID),]
table(taxonomy$OTU.ID ==rownames(otu))

# differential analysis
dds <- DESeqDataSetFromMatrix(countData = otu, colData = group2, design = ~type) 

# normaliz
dds <- DESeq(dds)
dds

# extract results
group='type'
control = 'MC'
treatment ='TS'
res <- results(dds, contrast=c(group, treatment, control))

# order results by P value
res = res[order(res$pvalue),]
res
summary(res)
table(res$padj<0.05)

# extract differential OTU
diff_OTU_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_OTU_deseq2)
head(diff_OTU_deseq2)
#write.csv(diff_OTU_deseq2, file= paste("DEOTU_",control,"_vs_",treatment,".csv"))

# calculate average relative abundance
abundance<-aggregate(t(otu_relative),by=list(group2$type),FUN=mean)
abundance<-as.data.frame(t(abundance))
colnames(abundance)<-abundance[1,]
abundance<-abundance[-1,]
abundance$MC <- as.numeric(abundance$MC)
abundance$TS <- as.numeric(abundance$TS)

#abundance<-as.data.frame(lapply(abundance, as.numeric))
abundance2 <- as.data.frame(apply(abundance, 1, function(x){mean(x)}))
names(abundance2) <- "abundance"; abundance2$OTU.ID <- rownames(abundance2)
taxonomy <- merge(taxonomy,abundance2,by="OTU.ID")

# merge data
test<- as.data.frame(res)
test$OTU.ID <- row.names(test)
data <- merge(test, taxonomy,by="OTU.ID",sort=FALSE,all=F)
data<-data[(data$abundance>0.0005),]

# transform p value
data$neglogp = -log10(data$padj)

data<-as.data.frame(cbind(data$OTU.ID, data$log2FoldChange, data$padj, data$V5_2, data$abundance, data$neglogp))
colnames(data)<-c("otu","log2FoldChange","pvalue","Phylum","abundance","neglogp")

#data$Genus <- substr(data$Genus,4,nchar(data$Genus))

# change data type
data<-transform(data, log2FoldChange = as.numeric(log2FoldChange),abundance = as.numeric(abundance), neglogp = as.numeric(neglogp), pvalue= as.numeric(pvalue))

# label differential OTU 
data$level = as.factor(ifelse(data$log2FoldChange>0, "Wildland > Farmland","Farmland > Wildland"))
data$level <- factor(data$level, levels=c("Wildland > Farmland","Farmland > Wildland"))
#write.csv(data, file= paste("OTU",control,"_vs_",treatment,".csv"))

# plot Manhattan
# data[data$neglogp>30,]$neglogp  = 30
Title=paste("Prokaryotes")

#col<-c("#0000FF","#FF3399","#FFCC33","#ff00ff","#00ff00", "deepskyblue", "darkgreen","black", "maroon3",   "grey")
col<-c("green","#FF3399","blue","#FFCC33","purple","black","darkgreen","chocolate", "red","deepskyblue","grey")

select <- as.data.frame(table(data[data$pvalue < 0.05,]$Phylum))[,1]
data <- data[data$Phylum %in% select,]
#data <- data[data$log2FoldChange > 2.5|data$log2FoldChange < -2.5,]
data <- na.omit(data)

blast2 <- readRDS("E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Figures/traits/blast_metadata.rda")

data2 <- data %>% separate(otu, into = c("otu"), "[.]")

data3 <- merge(data2, blast2[,c("V1", "genome_size")], by.x="otu",by.y="V1", all.x=TRUE)

set.seed(123)
print(data3[data3$neglogp >-log10(0.05),]$otu)



# 2. traits

setwd("E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Bacteria/otu/gtdb")

# read result blast
blast <- read.table("blast.uotus.bacteria.gtdb.txt",comment.char = "")  
blast <- blast[!duplicated(blast$V1),]

# read metadata of GTDB
metadata_bac <- read.delim("bac120_metadata_r214.tsv",sep = "\t",header = T)
metadata_arc <- read.delim("ar53_metadata_r214.tsv",sep = "\t",header = T)
metadata <- rbind(metadata_bac,metadata_arc)

blast1 <- blast %>% separate(V5, into = c("First", "Second"), "~")
blast1 <- blast1[!duplicated(blast1$V1),]  # delete duolicate OTUs

# merge results of blast and GTDB metadata
blast2 <- merge(blast1,metadata,by.x = "First",by.y="accession")
rownames(blast2) <- blast2$V1



# rrna 
setwd("E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Bacteria/otu/rrndb/")

blast_rrna <- read.delim("blast.uotus.bacteria.rrndb.all.txt", check.names = F,header = F)
blast_rrna <- blast_rrna[!duplicated(blast_rrna$V1),]

blast_rrna1 <-
  concat.split(
    data = blast_rrna,
    split.col = "V7",
    sep = "[|]",
    drop = TRUE
  )  # use ??/|?? for ??|??

rrndb <- read.delim("rrnDB-5.8.tsv", check.names = F)
names(rrndb)[1] <- "id"
names(rrndb)[12] <- "16S"
blast_rrna2 <- merge(blast_rrna1[, c("V1","V7_2")], rrndb[, c("id","16S")] ,by.x="V7_2", by.y="id")
blast_rrna2 <- as.data.frame(blast_rrna2)
rownames(blast_rrna2) <- blast_rrna2$V1


# plot

df_gs <- data.frame(id=data3[data3$neglogp >(-log10(0.05)),]$otu, value = blast2[data3[data3$neglogp > (-log10(0.05)),]$otu,]$genome_size)

df_gs$value <- df_gs$value/1000000

p_gs <- ggplot(df_gs, aes(x=value)) +
  geom_density(aes(value,..density..), size = 1) +
  geom_vline(aes(xintercept = median(value)), linetype = 5, size = 1) +
  geom_vline(aes(xintercept = 6.3), linetype = 5, size = 1,color="red") +
  annotate("text", x=4,y=0.08,label="Median", size=3.5, color="black") +
  annotate("text", x=8.1,y=0.235,label="OTU5~(italic(Acinetobacter))", parse=T, size=3.5, color="red") +
  annotate("text", x=8.1,y=0.219,label="Genome size = 6.3 Mb", size=3.5, color="red") +
  theme_bw() +
  xlim(2,10) +
  labs(x="Genome size", y="Density",title = "") +
  theme(
    plot.title = element_text(size = 20,color = "black",face = "bold", hjust = 0.5), 
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(colour = "black",size = 12,face = "bold"),
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(colour = "black", size = 12, face = "bold"),
    legend.spacing = unit(0.1, "cm"))
p_gs


df_gc <- data.frame(id=data3[data3$neglogp >(-log10(0.05)),]$otu, value = blast2[data3[data3$neglogp > (-log10(0.05)),]$otu,]$gc_percentage)
df_gc$value <- df_gc$value/100

p_gc <- ggplot(df_gc, aes(x=value)) +
  geom_density(aes(value,..density..), size = 1) +
  geom_vline(aes(xintercept = median(value)), linetype = 5, size = 1) +
  geom_vline(aes(xintercept = 0.6641), linetype = 5, size = 1,color="red") +
  scale_x_continuous(labels = scales::percent,limits = c(0.3, 0.9)) +
  annotate("text", x=0.51,y=0.8,label="Median", size=3.5, color="black") +
  annotate("text", x=0.81,y=3.7,label="OTU5~(italic(Acinetobacter))", parse=T, size=3.5, color="red") +
  annotate("text", x=0.81,y=3.45,label="GC content = 66.4%", size=3.5, color="red") +
  theme_bw() +
  labs(x="GC content", y="Density",title = "") +
  theme(
    plot.title = element_text(size = 20,color = "black",face = "bold", hjust = 0.5), 
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(colour = "black",size = 12,face = "bold"),
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(colour = "black", size = 12, face = "bold"),
    legend.spacing = unit(0.1, "cm"))
p_gc



df_rrna <- data.frame(id=data3[data3$neglogp > (-log10(0.05)),]$otu, value = blast_rrna2[data3[data3$neglogp > (-log10(0.05)),]$otu,"16S"])
df_rrna$value <- df_rrna$value

p_rrna <- ggplot(df_rrna, aes(x=value)) +
  geom_density(aes(value,..density..), size = 1) +
  geom_vline(aes(xintercept = median(value)), linetype = 5, size = 1) +
  geom_vline(aes(xintercept = 7), linetype = 5, size = 1,color="red") +
  annotate("text", x=3.5,y=0.02,label="Median", size=3.5, color="black") +
  annotate("text", x=10.5,y=0.08,label="OTU5~(italic(Acinetobacter))", parse=T, size=3.5, color="red") +
  annotate("text", x=10.5,y=0.075,label="rDNA copy number = 7", size=3.5, color="red") +
  theme_bw() +
  labs(x="rDNA copy number", y="Density",title = "") +
  theme(
    plot.title = element_text(size = 20,color = "black",face = "bold", hjust = 0.5), 
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(colour = "black",size = 12,face = "bold"),
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(colour = "black", size = 12, face = "bold"),
    legend.spacing = unit(0.1, "cm"))
p_rrna


saveRDS(p_gs,"E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Figures/traits/gs_dengsity.rda")
saveRDS(p_gc,"E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Figures/traits/gc_dengsity.rda")
saveRDS(p_rrna,"E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Figures/traits/rrna_dengsity.rda")


# h5.5 w12