# Import Read Counts
setwd("")

#Packages
library("tidyverse")
library("reshape2")
library("DESeq2")
library("gplots")
library("MASS")
library("mitch")
library("limma")
library("kableExtra")
library(readxl)
library(edgeR)
library(org.Hs.eg.db)
library(RColorBrewer)
library(FactoMineR)
library(factoextra)
library("EnhancedVolcano")
library(cowplot)

# PRE-PROCESSING
###########
#Import gene names
ENSG <- read.table("Homo_sapiens.GRCh38.tx2gene.tsv", header = F, fill = T)
head(ENSG)

#Import in read counts
tmp <- read.table("3col_ALL.tsv",header=F)
head(tmp)

#Merge df so that GeneName and geneID are included in dataframe.
tmpALL<-merge(tmp, ENSG, by.x="V2", by.y="V1")
tmpALL= tmpALL %>%
  rename("Transcript" = V2,
         "ID" = V1,
         "GeneCount"= V3.x,
         "GeneID"=V2.y,
         "GeneSYM"=V3.y)
tmpALL$GeneID = sub("\\..*", "", tmpALL$GeneID)
head(tmpALL)

#Create Matrix and aggregate according to ENSG
x <- as.matrix(acast(tmpALL, GeneID + GeneSYM  ~ ID, value.var= "GeneCount", fun.aggregate = sum))
x <- as.data.frame(x)
xx <- round(x)
head(xx)
write.csv(xx, "GeneCounts_RAW DATA.csv")

###########
xx=read_excel("GeneCounts_ RAW DATA_ALL.xlsx", sheet = 1)
xx=xx%>%
  tibble::column_to_rownames(var = "ENSEMBL")
xx= xx[,-1]

#Import in sample sheet
ss <- read_excel("Pheno.xlsx", sheet = 2)
ss = data.frame(ss, row.names=1)
ss$timepoint <- factor(ss$timepoint, levels= c("PRE", "POST"))
ss$ID <- factor(ss$ID)
ss$batch <- factor(ss$batch)

## QC analysis
par(mar=c(5,8,3,1))
sums <- colSums(xx)
sums <- sums[order(sums)]
barplot(sums,horiz=TRUE,las=1,xlab="num reads",cex.names=0.8)
abline(v=20000000,col="red")

#PCA PLOT
#Is there a batch effect? ID is nested within batch therefore we need to have ID in model and this will take care of the batch effect
M_norm_t=t(scale(xx))
res.pca = PCA(M_norm_t, graph = FALSE)
Time=ss$timepoint
Batch = ss$batch
ID= ss$ID

fviz_pca_ind(res.pca, geom = "point",
                 pointsize = 2,
                 habillage=ss$batch,
                 addEllipses=TRUE, 
                 ellipse.type = "confidence")+
  geom_point(aes(shape=Batch, colour= Batch))


####### Differential expression analysis ##########
xx <- xx[,match(rownames(ss) , colnames(xx))] 
all(rownames(ss) == colnames(xx)) 
xx1 <- xx[which(rowMeans(xx)>=10),] #Remove poorly detected genes 
dim(xx1)

ss$timepoint <- factor(ss$timepoint , levels = c("PRE","POST")) 
ss$ID <- factor(ss$ID) 

dds <- DESeqDataSetFromMatrix(countData = xx1, 
                              colData = ss, 
                              design = ~ ID + timepoint)
res <- DESeq(dds)
z <- results(res) # In this case we don't need to specify coefficient in results(dds) as it will be  based on the last variable in the design formula in this case is timepoint. 

zz <-cbind(as.data.frame(z),xx1) #Merge results and gene counts matrix
def <-as.data.frame(zz[order(zz$pvalue),])

#Volcano Plot
sig <- subset(def, padj < 0.05)
N_SIG=nrow(sig)
N_UP=nrow(subset(sig,log2FoldChange>0))
N_DN=nrow(subset(sig,log2FoldChange<0))

def = mutate(def, sig=ifelse(def$padj <0.05, "FDR<0.05", "Not Sig"))
def = mutate(def, coef=ifelse(def$stat <0, "neg", "pos"))
def = mutate(def, color=ifelse(def$padj>0.05,"black",ifelse(def$padj <0.05 & def$coef=="neg","blue","red")))

ggplot(def, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=color), size=1.5)+
  scale_color_manual(values=c("black","blue", "red"))+
  labs(title = "Differential gene expression with resistance training", 
       subtitle = paste(N_UP,"up-regulated genes &", N_DN, "down-regulated genes"),
       x="LogFC", y="-log10(p-value)")+
  coord_cartesian(clip = 'off')+#for the one infinite point
  labs(text=element_text(family="Times", size = 16))+
  theme_classic()+
  theme(legend.position = "none")

#ggpubr::ggarrange(VC)+
#tiff('Volcano Plot.tiff', width = 10, height = 7, units = 'in', res=600)
#dev.off()


# Visualisation
vsd <- vst(dds, blind=FALSE) #Apply variance stabilizing transformation
# In the above function calls, we specified blind = FALSE, which means that differences between ID and timepoint (the variables in the design) will not contribute to the expected variance-mean trend of the experiment.
plotPCA(vsd, "batch")

#visualize the transformed data with batch variation removed
design0 <- model.matrix(~ timepoint, data = ss)
assay(vsd) <- removeBatchEffect(assay(vsd), vsd$batch, design=design0)
plotPCA(vsd, "batch")


#P Value Histogram
hist(def$pvalue[def$baseMean > 1], breaks = 0:20/20,
     col = "blue", border = "white")


