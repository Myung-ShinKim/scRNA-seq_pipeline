library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(svglite)
library(DoubletFinder)

options(scipen=999)

HER2_data <- Read10X(data.dir = "filtered_feature_bc_matrix")
HER2.all <- CreateSeuratObject(counts = HER2_data, min.cells = 3, min.features = 200)
prot_genes <- read.table("Ntab.Prot_gene.list", header=T)
HER2 <- HER2.all[rownames(HER2.all) %in% prot_genes$geneID, ] ### 39482 protein coding genes + HER2 + mito + pltd

# Add number of genes per UMI for each cell to metadata
HER2$log10GenesPerUMI <- log10(HER2$nFeature_RNA) / log10(HER2$nCount_RNA)

# Compute percent mito and pltd ratio
HER2$mitoRatio <- PercentageFeatureSet(object = HER2, pattern = "^NitaM")
HER2$mitoRatio <- HER2@meta.data$mitoRatio / 100
HER2$pltdRatio <- PercentageFeatureSet(object = HER2, pattern = "^NitaC")
HER2$pltdRatio <- HER2@meta.data$pltdRatio / 100

# Compute percent HER2 ratio
HER2$HER2Ratio <- PercentageFeatureSet(object = HER2, pattern = "^anti-HER2-VHH")
HER2$HER2Ratio <- HER2@meta.data$HER2Ratio / 100

# Create metadata dataframe
metadata <- HER2@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)

# Create sample column
metadata$sample <- "Unfiltered_HER2"

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	facet_wrap(~sample)
ggsave("Figure3a.svg", width= 30, height=20, units="cm")

metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=pltdRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	facet_wrap(~sample)
ggsave("Figure3b.svg", width= 30, height=20, units="cm")

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  	ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = c(300,3000))
ggsave("Supplementary_Fig2a.svg", width= 20, height=20, units="cm")

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  	ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  	geom_boxplot() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells vs NGenes")
ggsave("Supplementary_Fig2b.svg", width= 20, height=20, units="cm")

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  	geom_density(alpha = 0.2) +
  	theme_classic()
ggsave("Supplementary_Fig2c.svg", width= 20, height=20, units="cm")

# Add metadata back to Seurat object
HER2@meta.data <- metadata
                           
# Filter out low quality reads and genes
# Comparison between scRNAbulk and bulkRNAseq
HER2_count <- GetAssayData(object = HER2, slot = "counts") %>% Matrix::rowSums()
scRNAbulk <- data.frame(cbind(HER2_count))
scRNAbulk$geneID <- rownames(scRNAbulk)

BulkRNAseq <- read.table("BulkRNAseq.count.txt", header=F, col.names=c("geneID","Bulk_count"))
rownames(BulkRNAseq) <- BulkRNAseq$geneID
Bulk_count <- BulkRNAseq[rownames(BulkRNAseq) %in% scRNAbulk$geneID,]
PseudoBulk <- merge(scRNAbulk, Bulk_count, by='geneID')

# Plot correlation between scRNAseq and bulkRNAseq
svglite("Supplementary_Fig2d.svg", width= 8, height=8)
plot(log(PseudoBulk$HER2_count+1,2),log(PseudoBulk$Bulk_count+1,2))
abline(lm(log(PseudoBulk$Bulk_count+1,2) ~ log(PseudoBulk$HER2_count+1,2)), col = "red", lwd = 3)
text(paste("Correlation:", round(cor(log(PseudoBulk$HER2_count+1,2), log(PseudoBulk$Bulk_count+1,2)), 2)), x = 3, y = 18)
dev.off()

# Make list of highly expressed genes (log2 fold change > 5)
filter1.HER2 <- subset(PseudoBulk, subset=log(PseudoBulk$HER2_count/PseudoBulk$Bulk_count,2) > 5)
filter2.HER2 <- filter1.HER2[!grepl("HER2",filter1.HER2$geneID),]
write.table(filter2.HER2$geneID, "Supplementary_Table3.txt", quote=F, row.names=F, col.names=F)

# Filtering low quality cells
boxplot(metadata$nGene)$stats
filtered_HER2 <- subset(x = HER2, 
                         subset= (nUMI >= 500) & (nUMI <= 10000) &
                           (nGene >= 300) & (nGene <= 3000) &
                           (log10GenesPerUMI > 0.90) & 
                           (mitoRatio < 0.05) &
                           (pltdRatio < 0.25))
VlnPlot(filtered_HER2,features=c("nUMI","nGene","mitoRatio","pltdRatio"), group.by="sample", ncol=2, pt.size=0)
ggsave("Figure3c.svg", width= 20, height=20, units="cm")

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_HER2, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 50 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 50  & !names(Matrix::rowSums(nonzero)) %in% filter2.HER2$geneID
summary(keep_genes)

# Keeping those genes expressed in more than 50 cells and removing stress-responsive genes during protoplasting
filtered_counts <- counts[keep_genes, ] 
dim(counts)
dim(filtered_counts)

# Reassign to filtered Seurat object and add metadata back to Seurat object
filtered_HER2 <- CreateSeuratObject(filtered_counts, meta.data = filtered_HER2@meta.data)
metadata1 <- filtered_HER2@meta.data
metadata1$sample <- "HER2"
filtered_HER2@meta.data <- metadata1

# cell cycle scoring (https://rnabio.org/module-08-scrna/0008/02/01/scRNA/)
cell.cycle.Ntab <- read.table("Cell_Cycle.Nta_Ath.txt", header=TRUE)
s.genes = unique(cell.cycle.Ntab$Ntab.Gene.ID[which(cell.cycle.Ntab$List == "G1/S")])
g2m.genes = unique(cell.cycle.Ntab$Ntab.Gene.ID[which(cell.cycle.Ntab$List == "G2/M")])
s_g2m_genes <- intersect(s.genes, g2m.genes)
s.genes <- setdiff(s.genes, s_g2m_genes)
g2m.genes <- setdiff(g2m.genes, s_g2m_genes)

# normalize, transformation, PCA, and UMAP
filtered_HER2 <- NormalizeData(filtered_HER2, verbose = TRUE)
filtered_HER2 <- CellCycleScoring(filtered_HER2, g2m.features=g2m.genes, s.features=s.genes)
filtered_HER2 <- SCTransform(filtered_HER2, vst.flavor = "v2", vars.to.regress = c("mitoRatio", "pltdRatio","S.Score","G2M.Score"))
filtered_HER2 <- RunPCA(filtered_HER2, npcs = 100, verbose = FALSE)
PCAPlot(filtered_HER2, group.by="Phase")
ggsave("Figure4a.svg", width= 25, height=20, units="cm")
filtered_HER2 <- RunUMAP(filtered_HER2, dims = 1:100, reduction = "pca")                        
DimPlot(filtered_HER2, reduction="umap", group.by="Phase")
ggsave("Figure4b.svg", width= 25, height=20, units="cm")

# Find doublet using DoubletFinder v3
sweep.res <- paramSweep_v3(filtered_HER2, sct=TRUE, num.cores=20) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
nExp <- round(ncol(filtered_HER2) * 0.04)
nExp
filtered_HER2 <- doubletFinder_v3(filtered_HER2, pN = 0.25, pK = pK_choose, nExp = nExp, PCs = 1:100, sct=TRUE)
DF.name = colnames(filtered_HER2@meta.data)[grepl("DF.classification", colnames(filtered_HER2@meta.data))]
filtered_HER2 = filtered_HER2[, filtered_HER2@meta.data[, DF.name] == "Singlet"]

# Clustering
filtered_HER2 <- FindNeighbors(filtered_HER2, reduction="pca", dims = 1:100, verbose = FALSE) %>%
    FindClusters(resolution = c(0.6,0.8,1.0), verbose = FALSE)
DimPlot(filtered_HER2, label = T, repel = T) # res 1
ggsave("Figure2a.svg", width = 20, height = 20, units = "cm")
DimPlot(filtered_HER2, split.by="Phase", ncol=2)
ggsave("Figure4c.svg", width = 20, height = 20, units = "cm")

# Feature plot
metrics <-  c("nCount_RNA", "nFeature_RNA", "mitoRatio", "pltdRatio")
FeaturePlot(filtered_HER2, reduction = "umap", features = metrics, pt.size = 0.4, order=TRUE, min.cutoff = 'q10',label = TRUE)
ggsave("Figure5.svg", width = 16, height = 16, units = "cm")

# Marker identification
options(scipen=1)
DefaultAssay(filtered_HER2) <- "RNA"
Idents(object = filtered_HER2) <- "SCT_snn_res.1"
Marker.res10 <- FindAllMarkers(filtered_HER2, only.pos=TRUE, min.pct=0.30, logfc.threshold=0.25)
Marker.res10 <- Marker.res10[Marker.res10$p_val_adj < 0.05,]

Known_Marker <- read.table("PlantscRNAdb/Leaf.Markers.txt",header=F,sep="\t")
s1 <- intersect(rownames(filtered_HER2),Known_Marker$V1)
s2 <- Known_Marker %>% filter(V1 %in% s1) %>% select(V1,V2,V3,V4)
s2 <- rename(s2, "gene"="V1", "AthID"="V2", "Type"="V3", "Degree"="V4")
res10 <- inner_join(Marker.res10, s2, by='gene')

EP <- res10[grep("Epidermal", res10$Type),]
EP1 <- EP[grep("_1", EP$Degree),]
DotPlot(filtered_HER2, features=c(unique(EP$gene)))+RotatedAxis()
ggsave("Supplementary_Fig1a.svg", width = 30, height = 20, units = "cm")

PC <- res10[grep("pavement", res10$Type),] 
DotPlot(filtered_HER2, features=c(unique(PC$gene)))+RotatedAxis() +
ggsave("Supplementary_Fig1b.svg", width = 8, height = 8, units = "cm")

VS <- res10[grep("Vasculature", res10$Type),]
DotPlot(filtered_HER2, features=c(unique(VS$gene)))+RotatedAxis()
ggsave("Supplementary_Fig1c.svg", width = 8, height = 8, units = "cm")

VlnPlot(filtered_HER2, features=c("LOC107768773","LOC107790609","LOC107774993","LOC107826207","LOC107831628","anti-HER2-VHH"), ncol=1)
ggsave("Figure2b.svg", width = 16, height = 30, units = "cm")

#Mesophyll LOC107768773     AT3G01500  CA1
#Epidermis LOC107790609     AT1G66400
#Pavement LOC107774993     AT1G13210
#Vasculature LOC107826207   AT2G16850
#Guard cell  LOC107831628   AT3G26744   SCRM




