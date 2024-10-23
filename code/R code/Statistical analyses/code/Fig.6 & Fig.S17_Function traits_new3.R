
#############################Fig.6 & Fig.S17_Function traits_new3##########################

library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(splitstackshape)
library(patchwork)
library(tibble)
library(psych)
library(reshape2)
library(ggpubr)
library(linkET)
library(agricolae) # kruskal for nonparametric test
library(ggh4x)  # set individual y axis for facet

rm(list=ls())
setwd("E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Bacteria/otu/gtdb")

# read rarefied otu table
otu_table <- otu_table <- read.csv("../otu_bacteriaFlattening.csv")
otu_table1 <- otu_table %>% separate(X, into = c("First", "Second"), "[.]")

# read result blast
blast <- read.table("blast.uotus.bacteria.gtdb.txt",comment.char = "")  
blast <- blast[!duplicated(blast$V1),]

# read metadata of GTDB
metadata_bac <- read.delim("bac120_metadata_r214.tsv",sep = "\t",header = T)
metadata_arc <- read.delim("ar53_metadata_r214.tsv",sep = "\t",header = T)
metadata <- rbind(metadata_bac,metadata_arc)

# filter the results according to otu table
blast1 <- blast[blast[,1] %in% otu_table1[,1],]

blast1 <- blast1 %>% separate(V5, into = c("First", "Second"), "~")
blast1 <- blast1[!duplicated(blast1$V1),]  # delete duolicate OTUs

# merge results of blast and GTDB metadata
blast2 <- merge(blast1,metadata,by.x = "First",by.y="accession")
otu_table2 <- merge(otu_table1,blast2,by.x = "First",by.y = "V1")
saveRDS(blast2,"E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Figures/traits/blast_metadata.rda")

### 1.1 calculate avergae genome size 
my_gs <- function(name){
  
  otu_table3 <- otu_table2
  otu_table3[,3:27] <- apply(otu_table3[,3:27],2,function(x) x/sum(x))  # percentage by column
  colSums(otu_table3[,3:27])  # check if the sum of each column is equal to 1
  
  # average genome size
  for(i in 3:27){  
    otu_table3[,i] <- otu_table3[,i]*otu_table3[,name]
  }
  
  genomesize_ave <- as.data.frame(colSums(otu_table3[,3:27]))
  
  genomesize_ave3 <- genomesize_ave
  genomesize_ave3$sample.id <- rownames(genomesize_ave3)
  names(genomesize_ave3)[1] <- "genome_size"
  genomesize_ave3 <<- genomesize_ave3
}

### 1.2 calculate avergae genome size Remove OTU5 
my_gs2 <- function(name){
  
  otu_table3 <- otu_table2[otu_table2$First != "OTU5",]
  otu_table3[,3:27] <- apply(otu_table3[,3:27],2,function(x) x/sum(x))  # percentage by column
  colSums(otu_table3[,3:27])  # check if the sum of each column is equal to 1
  
  # average genome size
  for(i in 3:27){  
    otu_table3[,i] <- otu_table3[,i]*otu_table3[,name]
  }
  
  genomesize_ave <- as.data.frame(colSums(otu_table3[,3:27]))
  
  genomesize_ave3 <- genomesize_ave
  genomesize_ave3$sample.id <- rownames(genomesize_ave3)
  names(genomesize_ave3)[1] <- "genome_size"
  genomesize_ave3_notu5 <<- genomesize_ave3
}

### 2.1 calculate avergae gc content
my_gc <- function(name){
  
  otu_table3 <- otu_table2
  otu_table3[,3:27] <- apply(otu_table3[,3:27],2,function(x) x/sum(x))  # percentage by column
  colSums(otu_table3[,3:27])  # check if the sum of each column is equal to 1
  
  # gc content
  for(i in 3:27){ 
    otu_table3[,i] <- otu_table3[,i]*otu_table3[,name]
  }
  
  gc_ave <- as.data.frame(colSums(otu_table3[,3:27]))
  
  gc_ave3 <- gc_ave
  gc_ave3$sample.id <- rownames(gc_ave3)
  
  names(gc_ave3)[1] <- "gc_content"
  gc_ave3 <<- gc_ave3
}


### 2.2 calculate avergae gc content Remove OTU5 
my_gc2 <- function(name){
  
  otu_table3 <- otu_table2[otu_table2$First != "OTU5",]
  otu_table3[,3:27] <- apply(otu_table3[,3:27],2,function(x) x/sum(x))  # percentage by column
  colSums(otu_table3[,3:27])  # check if the sum of each column is equal to 1
  
  # gc content
  for(i in 3:27){ 
    otu_table3[,i] <- otu_table3[,i]*otu_table3[,name]
  }
  
  gc_ave <- as.data.frame(colSums(otu_table3[,3:27]))
  
  gc_ave3 <- gc_ave
  gc_ave3$sample.id <- rownames(gc_ave3)
  
  names(gc_ave3)[1] <- "gc_content"
  gc_ave3_notu5 <<- gc_ave3
}

my_gs("genome_size")
my_gc("gc_percentage")

my_gs2("genome_size")
my_gc2("gc_percentage")


### 3. rRNA copy number from rrndb

#rm(list = setdiff(ls(), c("p1", "p2")))
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

otu_all0 <- read.csv("../otu_bacteriaFlattening.csv", row.names = 1)
otu_all0 <- as.data.frame(apply(otu_all0,2, function(x) x/sum(x)))
colSums(otu_all0)

otu_all0$id <- rownames(otu_all0)

otu_all1 <-
  concat.split(
    data = otu_all0,
    split.col = "id",
    sep = "\\.",
    drop = TRUE
  )

otu_all1 <- otu_all1[,1:26]

otu_all2 <- as.data.frame(merge(otu_all1,blast_rrna2, by.x="id_1", by.y="V1", all.x = TRUE))


my_rrna <- function(name){
  for(i in 2:26){ 
    otu_all2[,i] <- otu_all2[,i] * otu_all2[,name]
  }
  otu_all2 <- na.omit(otu_all2)
  
  rrna_ave <- as.data.frame(colSums(otu_all2[,2:26]))
  rrna_ave3 <- rrna_ave
  rrna_ave3$sample.id <- rownames(rrna_ave3)
  names(rrna_ave3)[1] <- "rrna"
  rrna_ave3 <<- rrna_ave3
}

my_rrna2 <- function(name){
  for(i in 2:26){ 
    otu_all2 <- otu_all2[otu_all2$id_1 != "OTU5",]
    otu_all2[,i] <- otu_all2[,i] * otu_all2[,name]
  }
  otu_all2 <- na.omit(otu_all2)
  
  rrna_ave <- as.data.frame(colSums(otu_all2[,2:26]))
  rrna_ave3 <- rrna_ave
  rrna_ave3$sample.id <- rownames(rrna_ave3)
  names(rrna_ave3)[1] <- "rrna"
  rrna_ave3_notu5 <<- rrna_ave3
}

my_rrna("16S")
my_rrna2("16S")

group <- read.csv("sample_metadata.csv")
group <- na.omit(group)
group$cnratio <- group$tc/group$tn
group$npration <- group$tn/group$tp
group$cpration <- group$tc/group$tp


df_list <- list(genomesize_ave3, gc_ave3, rrna_ave3, genomesize_ave3_notu5, gc_ave3_notu5, rrna_ave3_notu5)
df_list2 <- list(genomesize_ave3, gc_ave3, rrna_ave3,genomesize_ave3_notu5, gc_ave3_notu5, rrna_ave3_notu5,group[,c(1,8:16)])

# use Reduce() to megge multiple dataframe one time
traits <- Reduce(function(x, y) merge(x, y, by = "sample.id"), df_list)
traits <- column_to_rownames(traits,"sample.id")
group2 <- group[,c(1,9:17)]
group2 <- column_to_rownames(group2,"sample.id")

traits <- traits[rownames(group2),]

traits2 <- Reduce(function(x, y) merge(x, y, by = "sample.id"), df_list2)
traits2 <- column_to_rownames(traits2,"sample.id")

corr_matrix <- corr.test(traits,group2, method = 'spearman', adjust = 'fdr')
corr_matrix$r
corr_matrix$p

p1_cor <- qcorrplot(corr_matrix) +
  geom_square() +
  geom_mark(sep = '\n',sig_thres = 0.05, size = 2.5, color = "white") +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_x_discrete(labels = c(
    pH = "pH", 
    salt_content = "Salt content", 
    tc = "Total carbon", 
    tn = "Total nitrogen", 
    tp = "Total phosphorus", 
    avail_p = "Available phosphorus",
    cnratio = "Total carbon:Total nitrogen",
    npration = "Total nitrogen:Total phosphorus",
    cpration = "Total carbon:Total phosphorus"), 
    limit = c("pH","salt_content","tc","tn","tp","avail_p","cnratio","npration","cpration"), 
    position = "bottom") +
  scale_y_discrete(labels = c(
    genome_size.x = "Genome size (All OTUs)", 
    gc_content.x = "GC content (All OTUs)",
    rrna.x = "rDNA copy number (All OTUs)",
    genome_size.y = "Genome size (OTU5 removed)", 
    gc_content.y = "GC content (OTU5 removed)",
    rrna.y = "rDNA copy number (OTU5 removed)"
  ), limit = c("rrna.y","gc_content.y","genome_size.y","rrna.x","gc_content.x","genome_size.x")) +
  labs(title = "") +
  guides(fill = guide_colorbar(title = "Spearman's Rho")) +
  theme(plot.margin = margin(0,0,0,0, unit = "cm"), 
        #plot.title = element_text(size=20,face = "bold", hjust = 0.5),
        plot.title = element_blank(),
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 15,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"))

p1_cor


# data preperation

traits$id <- substr(rownames(traits), 1, 2)

## genome size_all
test_gs_all <- with(traits,kruskal(genome_size.x, id, p.adj="fdr"))
test_gc_all <- with(traits,kruskal(gc_content.x, id, p.adj="fdr"))
test_rrna_all <- with(traits,kruskal(rrna.x, id, p.adj="fdr"))

test_gs_all2 <- with(traits,kruskal(genome_size.y, id, p.adj="fdr"))
test_gc_all2 <- with(traits,kruskal(gc_content.y, id, p.adj="fdr"))
test_rrna_all2 <- with(traits,kruskal(rrna.y, id, p.adj="fdr"))


multcomp_gs <- data.frame(
  id = c(rownames(test_gs_all$groups),rownames(test_gs_all2$groups)),
  x = c(rep("genome_size.x",5), rep("genome_size.y",5)),
  y = c(test_gs_all$means[rownames(test_gs_all$groups),]$Max/1000000, 
        test_gs_all2$means[rownames(test_gs_all2$groups),]$Max/1000000),
  groups = c(test_gs_all$groups$groups,test_gs_all2$groups$groups))
multcomp_gs$id <- factor(multcomp_gs$id, levels = c("CK","JP","CL","MH","YM"))


traits_genomesize <- melt(traits[,c(1,4,7)],id="id")

p_gs_all <- ggplot(data = traits_genomesize,aes(y=value/1000000,color=factor(id, levels = c("CK","JP","CL","MH","YM")), x=factor(variable)))+
  stat_boxplot(geom = "errorbar",
               width=0.7)+
  #  ymaxdata <- layer_data(p)
  geom_boxplot(width=0.7)+
  #  geom_jitter(size=0.8)+
  scale_color_manual(name = "Habitat", 
                     values = c("red","blue","purple","darkgreen","deepskyblue"),
                     labels = c( "CK"="Bareland", 
                                 "JP"="Suaeda",
                                 "CL"="Tamarix",
                                 "MH"="Cotton", 
                                 "YM"="Maize")) +
  scale_x_discrete(labels=c(
    "genome_size.x" = "All OTUs",
    "genome_size.y" = "OTU5 removed"
  )) +
  labs(x="", y="Mb", title = "Average genome size")+
  geom_text(data=multcomp_gs,
            aes(x=x,y=y,label=groups, group = id), color="black", size = 4, position=position_dodge2(width = .7), 
            vjust=-1) +
  
  annotate("text",x="genome_size.x", y=5.8, 
           label=c(sprintf("italic(P)  ==  %.3g ~ (italic(chi)^2  ==  %.3g)", 
                   test_gs_all$statistics[2],as.numeric(test_gs_all$statistic[1]))), parse=T,color="black",size=4) +
  annotate("text",x="genome_size.y", y=5.7,
           label=c(sprintf("italic(P)  ==  %.3g ~ (italic(chi)^2  ==  %.3g)", 
                   test_gs_all2$statistics[2],as.numeric(test_gs_all2$statistics[1]))), parse=T,color="black",size=4) + 
  theme_bw() +
  ylim(4.9,5.9) +
  theme(panel.grid=element_blank(), 
        axis.ticks.y=element_line(color="black",linewidth=0.5))+
  theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),
        axis.title.x = element_blank(),
        axis.text=element_text(size=12,face="bold",color = "black"),
        legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"))
p_gs_all


### rrna_all
multcomp_rrna <- data.frame(
  id = c(rownames(test_rrna_all$groups),rownames(test_rrna_all2$groups)),
  x = c(rep("rrna.x",5), rep("rrna.y",5)),
  y = c(test_rrna_all$means[rownames(test_rrna_all$groups),]$Max, 
        test_rrna_all2$means[rownames(test_rrna_all2$groups),]$Max),
  groups = c(test_rrna_all$groups$groups,test_rrna_all2$groups$groups))
multcomp_rrna$id <- factor(multcomp_rrna$id, levels = c("CK","JP","CL","MH","YM"))


traits_rrna <- melt(traits[,c(3,6,7)],id="id")


#df3 <- left_join(df,multcomp.data,by=c("id"="x"))

p_rrna_all <- ggplot(data = traits_rrna,aes(y=value,color=factor(id, levels = c("CK","JP","CL","MH","YM")), x=factor(variable)))+
  stat_boxplot(geom = "errorbar",
               width=0.7)+
  #  ymaxdata <- layer_data(p)
  geom_boxplot(width=0.7)+
  #  geom_jitter(size=0.8)+
  scale_color_manual(name = "Habitat", 
                     values = c("red","blue","purple","darkgreen","deepskyblue"),
                     labels = c( "CK"="Bareland", 
                                 "JP"="Suaeda",
                                 "CL"="Tamarix",
                                 "MH"="Cotton", 
                                 "YM"="Maize")) +
  scale_x_discrete(labels=c(
    "rrna.x" = "All OTUs",
    "rrna.y" = "OTU5 removed"
  )) +
  labs(x="", y="", title = "rDNA copy number")+
  geom_text(data=multcomp_rrna,
            aes(x=x,y=y,label=groups, group = id), color="black", size = 4, position=position_dodge2(width = .7), 
            vjust=-1) +
  
  annotate("text",x="rrna.x", y=4.8, 
           label=c(sprintf("italic(P)  ==  %.3g ~ (italic(chi)^2  ==  %.3g)", 
                           test_rrna_all$statistics[2],as.numeric(test_rrna_all$statistic[1]))), parse=T,color="black",size=4) +
  annotate("text",x="rrna.y", y=3.8,
           label=c(sprintf("italic(P)  ==  %.3g ~ (italic(chi)^2  ==  %.3g)", 
                           test_rrna_all2$statistics[2],as.numeric(test_rrna_all2$statistics[1]))), parse=T,color="black",size=4) + 
  theme_bw() +
  ylim(2,5) +
  theme(panel.grid=element_blank(), 
        axis.ticks.y=element_line(color="black",linewidth=0.5))+
  theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),
        axis.title.x = element_blank(),
        axis.text=element_text(size=12,face="bold",color = "black"),
        legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"))
p_rrna_all


### gc_all
multcomp_gc <- data.frame(
  id = c(rownames(test_gc_all$groups),rownames(test_gc_all2$groups)),
  x = c(rep("gc_content.x",5), rep("gc_content.y",5)),
  y = c(test_gc_all$means[rownames(test_gc_all$groups),]$Max, 
        test_gc_all2$means[rownames(test_gc_all2$groups),]$Max),
  groups = c(test_gc_all$groups$groups,test_gc_all2$groups$groups))
multcomp_gc$id <- factor(multcomp_gc$id, levels = c("CK","JP","CL","MH","YM"))


traits_gc <- melt(traits[,c(2,5,7)],id="id")


#df3 <- left_join(df,multcomp.data,by=c("id"="x"))

p_gc_all <- ggplot(data = traits_gc,aes(y=value,color=factor(id, levels = c("CK","JP","CL","MH","YM")), x=factor(variable)))+
  stat_boxplot(geom = "errorbar",
               width=0.7)+
  #  ymaxdata <- layer_data(p)
  geom_boxplot(width=0.7)+
  #  geom_jitter(size=0.8)+
  scale_color_manual(name = "Habitat", 
                     values = c("red","blue","purple","darkgreen","deepskyblue"),
                     labels = c( "CK"="Bareland", 
                                 "JP"="Suaeda",
                                 "CL"="Tamarix",
                                 "MH"="Cotton", 
                                 "YM"="Maize")) +
  scale_x_discrete(labels=c(
    "gc_content.x" = "All OTUs",
    "gc_content.y" = "OTU5 removed"
  )) +
  labs(x="", y="", title = "GC content (%)")+
  geom_text(data=multcomp_gc,
            aes(x=x,y=y,label=groups, group = id), color="black", size = 4, position=position_dodge2(width = .7), 
            vjust=-1) +
  
  annotate("text",x="gc_content.x", y=65, 
           label=c(sprintf("italic(P)  ==  %.3g ~ (italic(chi)^2  ==  %.3g)", 
                           test_gc_all$statistics[2],as.numeric(test_gc_all$statistic[1]))), parse=T,color="black",size=4) +
  annotate("text",x="gc_content.y", y=64.8,
           label=c(sprintf("italic(P)  ==  %.3g ~ (italic(chi)^2  ==  %.3g)", 
                           test_gc_all2$statistics[2],as.numeric(test_gc_all2$statistics[1]))), parse=T,color="black",size=4) + 
  theme_bw() +
  ylim(59,66) +
  theme(panel.grid=element_blank(), 
        axis.ticks.y=element_line(color="black",linewidth=0.5))+
  theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),
        axis.title.x = element_blank(),
        axis.text=element_text(size=12,face="bold",color = "black"),
        legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"))
p_gc_all



## data praperation for two groups
# group
traits3 <- traits[traits$id != "CK",]
traits3[traits3$id == "JP" | traits3$id == "CL",]$id <- "Wildland"
traits3[traits3$id == "YM" | traits3$id == "MH",]$id <- "Farmland"

# t test
test_gs <- wilcox.test(genome_size.x~id,traits3)
test_gc <- wilcox.test(gc_content.x~id,traits3)
test_rrna <- wilcox.test(rrna.x~id,traits3)

test_gs2 <- wilcox.test(genome_size.y~id,traits3)
test_gc2 <- wilcox.test(gc_content.y~id,traits3)
test_rrna2 <- wilcox.test(rrna.y~id,traits3)


traits4_gs <- melt(traits3[,c(1,4,7)],id="id")

p_gs <- ggplot(data = traits4_gs,aes(y=value/1000000,color=factor(id, levels = c("Wildland", "Farmland")), x=factor(variable)))+
  stat_boxplot(geom = "errorbar",
               width=0.4)+
  #  ymaxdata <- layer_data(p)
  geom_boxplot(width=0.4)+
  #  geom_jitter(size=0.8)+
  scale_x_discrete(labels=c(
    'genome_size.x' = "All OTUs",
    'genome_size.y' = "OTU5 Removed"
  )) +
  scale_color_manual(name="Type", values = c("red", "blue"), labels=c(
    "Farmland" = "Farmland (Cotton & Maize)",
    "Wildland" = "Wildland (Tamarix & Suaeda)"
  )) +
  labs(x="", y="", title = "Genome size (Mb)")+
  annotate("text",x="genome_size.x", y=5.65,label=c(sprintf("italic(P)  ==  %.3g ~ (italic(W)  ==  %.3g)", 
                                                            test_gs$p.value,as.numeric(test_gs$statistic))), parse=T,color="black",size=3.5) +
  annotate("text",x="genome_size.y", y=5.45,label=c(sprintf("italic(P)  ==  %.3g ~ (italic(W)  ==  %.3g)", 
                                                              test_gs2$p.value,as.numeric(test_gs2$statistic))), parse=T,color="black",size=3.5) + 
  theme_bw()+
  #ylim(4.9,5.8)+
  theme(panel.grid=element_blank(), 
        axis.ticks.y=element_line(color="black",linewidth=0.5))+
  theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
        #axis.title=element_text(size=16,face="bold"),
        axis.title = element_blank(),
        axis.text=element_text(size=12,face="bold",color = "black"),
        legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"))
p_gs



## gc_two

traits4_gc <- melt(traits3[,c(2,5,7)],id="id")

p_gc <- ggplot(data = traits4_gc,aes(y=value,color=factor(id,levels = c("Wildland", "Farmland")), x=factor(variable)))+
  stat_boxplot(geom = "errorbar",
               width=0.4)+
  #  ymaxdata <- layer_data(p)
  geom_boxplot(width=0.4)+
  #  geom_jitter(size=0.8)+
  scale_x_discrete(labels=c(
    'gc_content.x' = "All OTUs",
    'gc_content.y' = "OTU5 Removed"
  )) +
  scale_color_manual(name="Type", values = c("red", "blue"), labels=c(
    "Farmland" = "Farmland (Cotton & Maize)",
    "Wildland" = "Wildland (Tamarix & Suaeda)"
  )) +
  labs(x="", y="", title = "GC content (%)")+
  annotate("text",x="gc_content.x", y=61.9,label=c(sprintf("italic(P)  ==  %.3g ~ (italic(W)  ==  %.3g)", 
                                                            test_gc$p.value,as.numeric(test_gc$statistic))), parse=T,color="black",size=3.5) +
  annotate("text",x="gc_content.y", y=61.7,label=c(sprintf("italic(P)  ==  %.3g ~ (italic(W)  ==  %.3g)", 
                                                            test_gc2$p.value,as.numeric(test_gc2$statistic))), parse=T,color="black",size=3.5) + 
  theme_bw()+
  #ylim(4.9,5.8)+
  theme(panel.grid=element_blank(), 
        axis.ticks.y=element_line(color="black",linewidth=0.5))+
  theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
        #axis.title=element_text(size=16,face="bold"),
        axis.title = element_blank(),
        axis.text=element_text(size=12,face="bold",color = "black"),
        legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"))
p_gc



## rrna_two

traits4_rrna <- melt(traits3[,c(3,6,7)],id="id")

p_rrna <- ggplot(data = traits4_rrna,aes(y=value,color=factor(id,levels = c("Wildland", "Farmland")), x=factor(variable)))+
  stat_boxplot(geom = "errorbar",
               width=0.4)+
  #  ymaxdata <- layer_data(p)
  geom_boxplot(width=0.4)+
  #  geom_jitter(size=0.8)+
  scale_x_discrete(labels=c(
    'rrna.x' = "All OTUs",
    'rrna.y' = "OTU5 Removed"
  )) +
  scale_color_manual(name="Type", values = c("red", "blue"), labels=c(
    "Farmland" = "Farmland (Cotton & Maize)",
    "Wildland" = "Wildland (Tamarix & Suaeda)"
  )) +
  labs(x="", y="", title = "rDNA copy number")+
  annotate("text",x="rrna.x", y=4.4,label=c(sprintf("italic(P)  ==  %.3g ~ (italic(W)  ==  %.3g)", 
                                                           test_rrna$p.value,as.numeric(test_rrna$statistic))), parse=T,color="black",size=3.5) +
  annotate("text",x="rrna.y", y=3.4,label=c(sprintf("italic(P)  ==  %.3g ~ (italic(W)  ==  %.3g)", 
                                                           test_rrna2$p.value,as.numeric(test_rrna2$statistic))), parse=T,color="black",size=3.5) + 
  theme_bw()+
  #ylim(4.9,5.8)+
  theme(panel.grid=element_blank(), 
        axis.ticks.y=element_line(color="black",linewidth=0.5))+
  theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
        #axis.title=element_text(size=16,face="bold"),
        axis.title = element_blank(),
        axis.text=element_text(size=12,face="bold",color = "black"),
        legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"))
p_rrna



saveRDS(p1_cor,"E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Figures/traits/Correlation.rda")
saveRDS(p_gs_all,"E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Figures/traits/Genome_size_five.rda")
saveRDS(p_gc_all,"E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Figures/traits/GC_content_five.rda")
saveRDS(p_rrna_all,"E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Figures/traits/rRNA_five.rda")
saveRDS(p_gs,"E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Figures/traits/GS_two.rda")
saveRDS(p_gc,"E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Figures/traits/GC_two.rda")
saveRDS(p_rrna,"E:/Saline-alkali soil/Saline-alkali/amplicon/result-U/Figures/traits/rRNA_two.rda")

