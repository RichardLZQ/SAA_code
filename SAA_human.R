## This script included human data import and analysis
## Karmath etc and Semajic etc two datasets were included here

library(Seurat)
library(tidyverse)
library(harmony)
library(homologene)
library(biomaRt)

#### Karmath etc ####

#1. Import data
karmath <- Read10X("./Kamath/RAW/") # Load paper's data

karmath <- CreateSeuratObject(counts = count)
ac.label <- read_tsv("./Kamath/astro_UMAP.tsv") # Extract the AC barcode only keep AC cells 
kar.ast <- karmath[,(colnames(karmath) %in% ac.label$NAME)]
kar.ast <- NormalizeData(kar.ast) 
kar.ast <- FindVariableFeatures(kar.ast) 
kar.ast <- ScaleData(kar.ast,vars.to.regress = "donor_id")
kar.ast <- RunPCA(kar.ast)
kar.ast <- RunUMAP(kar.ast,dims = 1:20,umap.method = "umap-learn")

ac.label2 <- as.matrix(ac.label[,c(2,3)])
ac.label2 <- ac.label2[-1,]
ac.label2[,1] <- as.numeric(ac.label2[,1])
ac.label2[,2] <- as.numeric(ac.label2[,2])
rownames(ac.label2) <- ac.label$NAME[-1]
colnames(ac.label2) <- c(1,2)

kar.ast[["umap2"]] <- CreateDimReducObject(embeddings = ac.label2,key = "UMAP2_", assay = DefaultAssay(kar.ast))

kar.ast <- RunHarmony(kar.ast,group.by.vars = "donor_id",dims.use = 1:30)

#2. Add meta data
kar.meta <- read_tsv('./Kamath/METADATA_PD.tsv')
kar.meta <- kar.meta[kar.meta$NAME %in% colnames(kar.ast),]
kar.meta <- as.data.frame(kar.meta)
rownames(kar.meta) <- kar.meta$NAME
kar.meta <- kar.meta[,-1]
kar.ast <- AddMetaData(kar.ast,kar.meta)

# 3. Inferring the SAA
AC.fc.dat2 <- CreateSeuratObject(assay = "SCT", AC.fc.dat@assays$SCT@data) # Here we used frontal cortex data as reference dataset
AC.fc.dat2 <- FindVariableFeatures(AC.fc.dat2)
AC.fc.dat2 <- ScaleData(AC.fc.dat2)
AC.fc.dat2 <- RunPCA(AC.fc.dat2)
AC.fc.dat2$def_1 <- AC.fc.dat$def_1

kar.hu.mtx <- kar.ast@assays$RNA@data
kar.hu.mtx <- kar.hu.mtx[kar.ast@assays$RNA@var.features,]
kar.hu.gene <- rownames(kar.hu.mtx)
kar.hu.gene <- human2mouse(kar.hu.gene)
kar.hu.gene <- kar.hu.gene[!duplicated(kar.hu.gene$humanGene),]
kar.hu.mtx <- kar.hu.mtx[rownames(kar.hu.mtx) %in% kar.hu.gene$humanGene,]
rownames(kar.hu.mtx) <- kar.hu.gene$mouseGene[match(rownames(kar.hu.mtx),kar.hu.gene$humanGene)]
kar.hu.mtx <- kar.hu.mtx[!duplicated(rownames(kar.hu.mtx)),]

kar.hu.ast <- CreateSeuratObject(assay = "integrated",kar.hu.mtx)
kar.hu.ast <- FindVariableFeatures(kar.hu.ast)
kar.hu.ast <- ScaleData(kar.hu.ast)
kar.hu.ast <- RunPCA(kar.hu.ast)

spe.anchor <- FindTransferAnchors(reference = AC.fc.dat2,query = kar.hu.ast,project.query = TRUE)
predictions <- TransferData(anchorset = spe.anchor, refdata = AC.fc.dat2$def_1)

def1 <- predictions$predicted.id
names(def1) <- colnames(kar.ast)
kar.ast$def_mo <- def1

##### Semajc etc ####

# 1. Load the data

path2data <- "./"

count <- read_tsv("./matrix.tsv") # Convert the tsv to mtx
gene <- read_tsv("./genes.tsv")
cell <- read_tsv("./barcodes.tsv")

cell <- column_to_rownames(cell,"barcode")

dataset="hsapiens_gene_ensembl"  # Prepare a dic for gene name converting
mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset, 
                        host = paste0("feb2021", ".archive.ensembl.org"), path = "/biomart/martservice", archive = FALSE) # Convert the gene name to symbol name
annotations <- biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id", "external_gene_name"))
annotations <- annotations %>% filter(ensembl_gene_id %in% gene$gene)



count <- as.data.frame(count)
rownames(count) <- gene$gene
colnames(count) <- rownames(cell)

prada <- CreateSeuratObject(counts = count ,project = "Prada",meta.data = cell)

prada.list <- SplitObject(prada, split.by = "patient")

# Normalize and remove cell cycle + batch effect
s.genes <- cc.genes$s.genes
s.genes <- annotations$ensembl_gene_id[match(s.genes,annotations$external_gene_name)]
s.genes[4] <- "ENSG00000168496"
s.genes[15] <- "ENSG00000151725"
s.genes[29] <- "ENSG00000051180"

g2m.genes <- cc.genes$g2m.genes
g2m.genes <- annotations$ensembl_gene_id[match(g2m.genes,annotations$external_gene_name)]
g2m.genes[16] <- "ENSG00000129195"
g2m.genes[30] <- "ENSG00000189159"

for (i in 1:length(prada.list)) {
  prada.list[[i]] <- SCTransform(prada.list[[i]],variable.features.n = 4000)
  prada.list[[i]] <- CellCycleScoring(prada.list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
}

# find anchors
anchors <- FindIntegrationAnchors(object.list = prada.list,anchor.features = 4000)

# integrate data
prada <- IntegrateData(anchorset = anchors)
prada.bp <- prada # Backup the data
# Standard scRNA-seq data analysis pipeline
prada <- ScaleData(prada)
prada <- RunPCA(prada)
prada <- RunUMAP(prada,dims = 1:25,umap.method = "umap-learn")
prada <- FindNeighbors(prada,dims = 1:25)
prada <- FindClusters(prada,resolution = 1.5)


# Extract the AC cells data and redo above steps
Idents(prada) <- "cell_ontology"
ac.prada <- subset(prada,idents = "Astrocytes")


ac.prada.list <- SplitObject(ac.prada, split.by = "patient")
for (i in 1:length(ac.prada.list)) {
  ac.prada.list[[i]] <- SCTransform(ac.prada.list[[i]],variable.features.n = 4000)
  #ac.prada.list[[i]] <- CellCycleScoring(ac.prada.list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
}
anchors <- FindIntegrationAnchors(object.list = ac.prada.list,anchor.features = 4000)
ac.prada <- IntegrateData(anchorset = anchors)
ac.prada <- ScaleData(ac.prada)
ac.prada <- RunPCA(ac.prada)
ac.prada <- RunUMAP(ac.prada,dims = 1:25,umap.method = "umap-learn")
ac.prada <- FindNeighbors(ac.prada,dims = 1:25)
ac.prada <- FindClusters(ac.prada,resolution = 1.0)

old.gene <- ac.prada@assays$integrated@data@Dimnames[[1]]
new.gene <- annotations[match(old.gene,annotations$ensembl_gene_id),2]
new.gene[duplicated(new.gene)] <- paste0(new.gene[duplicated(new.gene)],1:10)
new.gene <- str_split_fixed(string = new.gene,"\\.",2)[,1]
new.gene.mtx <- human2mouse(new.gene)
new.gene.mtx <- new.gene.mtx[!duplicated(new.gene.mtx$mouseGene),]
ac.prada.mo <- ac.prada[new.gene %in% new.gene.mtx$humanGene,]
ac.prada.mo@assays$integrated@data@Dimnames[[1]] <- new.gene.mtx$mouseGene

new.gene.hu <- human2mouse(rownames(humanAC.mtx))
new.gene.hu <- new.gene.hu[!duplicated(new.gene.hu),]
humanAC2.mtx <- humanAC.mtx[rownames(humanAC.mtx)%in% new.gene.hu$humanGene,]
rownames(humanAC2.mtx) <- new.gene.hu$mouseGene[match(rownames(humanAC2.mtx),new.gene.hu$humanGene)]
humanAC2.mtx <- humanAC2.mtx[!duplicated(rownames(humanAC2.mtx)),]

ac.prada.new <- CreateSeuratObject(assay = "integrated",humanAC2.mtx)
ac.prada.new <- FindVariableFeatures(ac.prada.new)
ac.prada.new <- ScaleData(ac.prada.new)
ac.prada.new <- RunPCA(ac.prada.new)

new.name <- str_split_fixed(colnames(AC.fc.dat),"_",2)[,2]
AC.fc.dat2 <- CreateSeuratObject(assay = "SCT", AC.fc.dat@assays$SCT@data)
AC.fc.dat2 <- FindVariableFeatures(AC.fc.dat2)
AC.fc.dat2 <- ScaleData(AC.fc.dat2)
AC.fc.dat2 <- RunPCA(AC.fc.dat2)
AC.fc.dat2$def_1 <- AC.fc.dat$def_1

spe.anchor <- FindTransferAnchors(reference = AC.fc.dat2,query = ac.prada.new,project.query = TRUE)
predictions <- TransferData(anchorset = spe.anchor, refdata = AC.fc.dat2$def_1)

def1 <- predictions$predicted.id
names(def1) <- colnames(ac.prada.new)
ac.prada$def_mo <- def1

ac.prada$SAAscore <- predictions$prediction.score.SAA

Idents(ac.prada) <- "def_mo"
hu.acmrks <- FindAllMarkers(ac.prada,logfc.threshold = 0.05)
hu.acmrks2 <- hu.acmrks %>% filter(p_val_adj <0.05) %>% group_by(cluster) %>% top_n(50,avg_logFC)
hu.acmrks2$Symbol <- annotations$external_gene_name[match(hu.acmrks2$gene,annotations$ensembl_gene_id)]