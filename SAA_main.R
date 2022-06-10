# The mian code for SAA project analysis
# Based on new cellranger and seurat(3.1.1) version
# Ko lab 2021

# Start from here
# Load library
library(Seurat)
library(ggplot2)
library(sctransform)
library(tidyverse)
library(venn)
library(EnhancedVolcano)
library(ggpmisc)

#### 1. Loading raw data and merge by regions ####
tg.fc.data <- Read10X(data.dir = "~/bioinfo/Caro/Astro/Aligned/Batch5_SNCA/5Tg1/")
tg.mb.data <- Read10X(data.dir = "~/bioinfo/Caro/Astro/Aligned/Batch5_SNCA/5Tg2/")
wt.fc.data <- Read10X(data.dir = "~/bioinfo/Caro/Astro/Aligned/Batch5_SNCA/5Lm1/")
wt.mb.data <- Read10X(data.dir = "~/bioinfo/Caro/Astro/Aligned/Batch5_SNCA/5Lm2/")

tg.fc <- CreateSeuratObject(counts=tg.fc.data,min.cells = 5,project = "tg.fc")
tg.mb <- CreateSeuratObject(counts=tg.mb.data,min.cells = 5,project = "tg.mb")
wt.fc <- CreateSeuratObject(counts=wt.fc.data,min.cells = 5,project = "wt.fc")
wt.mb <- CreateSeuratObject(counts=wt.mb.data,min.cells = 5,project = "wt.mb")

tg.fc <- RenameCells(tg.fc, add.cell.id = "tgfc")
tg.mb <- RenameCells(tg.mb, add.cell.id = "tgmb")
wt.fc <- RenameCells(wt.fc, add.cell.id = "wtfc")
wt.mb <- RenameCells(wt.mb, add.cell.id = "wtmb")

ast.fc <- merge(x = tg.fc,y = wt.fc ,project = "b5_fc")
ast.mb <- merge(x = tg.mb,y = wt.mb, project = "b5_mb")

#### 2. Data normalization and dimention reduction####
ast.fc <- PercentageFeatureSet(ast.fc, pattern = "^mt-", col.name = "percent.mt")
ast.mb <- PercentageFeatureSet(ast.mb, pattern = "^mt-", col.name = "percent.mt")

ast.fc <- SCTransform(ast.fc, verbose = T)
ast.mb <- SCTransform(ast.mb, verbose = T)

ast.fc <- RunPCA(ast.fc)
ast.fc <- RunUMAP(ast.fc,dims = 1:30,umap.method = "umap-learn")
ast.fc <- FindNeighbors(ast.fc,dims = 1:30)
ast.fc <- FindClusters(ast.fc,resolution = 1.0)

ast.mb <- RunPCA(ast.mb)
ast.mb <- RunUMAP(ast.mb,dims = 1:30,umap.method = "umap-learn")
ast.mb <- FindNeighbors(ast.mb,dims = 1:30)
ast.mb <- FindClusters(ast.mb,resolution = 1.0)

#### 3. Clustering and cell type identification and quality control ####
VlnPlot(ast.fc,features = 'percent.mt') # 19 removed due to abnormal mito gene percent
VlnPlot(ast.mb,features = 'percent.mt') # 23 removed due to abnormal mito gene percent

Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("12")))<-"PC"
Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("9")))<-"SMC"
Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("0","5","6","14")))<-"EC"
Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("2","3","7","8")))<-"MG"
Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("1","4","11","22")))<-"AC"
Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("23","27")))<-"OPC"
Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("21")))<-"OLG"
Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("16")))<-"imNeur"
Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("13")))<-"mNeur"
Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("29")))<-"EPC"
Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("20")))<-"CPC"
Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("17")))<-"Hb_EC"
Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("15")))<-"MAC"
Idents(ast.fc, cells = WhichCells(object = ast.fc, idents = c("25")))<-"MNC"

Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("11")))<-"PC"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("7","26")))<-"SMC"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("0","9","6","14")))<-"EC"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("2","16","15","21")))<-"MG"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("1","3","4","10","22","5")))<-"AC"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("25")))<-"OPC"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("30")))<-"OLG"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("29")))<-"imNeur"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("20")))<-"mNeur"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("19")))<-"EPC"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("8","12","13")))<-"CPC"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("24")))<-"Hb_EC"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("18")))<-"MAC"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("28")))<-"TNC"
Idents(ast.mb, cells = WhichCells(object = ast.mb, idents = c("27")))<-"MNC"

ast.fc$Celltype <- Idents(ast.fc)
ast.mb$Celltype <- Idents(ast.mb)

ast.fc <- subset(ast.fc, idents =c("PC","SMC","EC","MG","AC","OPC","OLG","imNeur","mNeur","EPC","CPC","Hb_EC","MAC","MNC"))
ast.mb <- subset(ast.mb, idents =c("PC","SMC","EC","MG","AC","OPC","OLG","imNeur","mNeur","EPC","CPC","Hb_EC","MAC","TNC","MNC"))

#### 4. The cell type prioritization calculate ####
library(Augur)
augur.fc = calculate_auc(ast.fc,  cell_type_col = "Celltype", label_col = "orig.ident",n_threads = 16)
augur.mb = calculate_auc(ast.mb,  cell_type_col = "Celltype", label_col = "orig.ident",n_threads = 16)

#### 5. Sub-cluster analysis in AC data ####

AC.fc.dat <- SCTransform(AC.fc.dat,verbose = T)
AC.mb.dat <- SCTransform(AC.mb.dat,verbose = T)

AC.fc.dat <- RunPCA(AC.fc.dat)
AC.fc.dat <- RunUMAP(AC.fc.dat,dims = 1:30,umap.method = "umap-learn")
AC.fc.dat <- FindNeighbors(AC.fc.dat,dims = 1:30)
AC.fc.dat <- FindClusters(AC.fc.dat,resolution = 0.8)

AC.mb.dat <- RunPCA(AC.mb.dat)
AC.mb.dat <- RunUMAP(AC.mb.dat,dims = 1:30,umap.method = "umap-learn")
AC.mb.dat <- FindNeighbors(AC.mb.dat,dims = 1:30)
AC.mb.dat <- FindClusters(AC.mb.dat,resolution = 0.8)
AC.mb.dat <- FindClusters(AC.mb.dat,resolution = 1)

#### 5.1 Re-cluster based on signature expression level ####

Idents(AC.fc.dat) <- "SCT_snn_res.0.8"
Idents(AC.fc.dat, cells = WhichCells(object = AC.fc.dat, idents = c("1","7")))<-"SAA"
Idents(AC.fc.dat, cells = WhichCells(object = AC.fc.dat, idents = c("4","2")))<-"Gfap_high" #Add 8
Idents(AC.fc.dat, cells = WhichCells(object = AC.fc.dat, idents = c("0","3","5","6","8")))<-"Gfap_low"
AC.fc.dat$def_1 <- Idents(AC.fc.dat)

Idents(AC.fc.dat) <- "def_1"
AC.fc.dat <- subset(AC.fc.dat,idents = c("SAA","Gfap_high","Gfap_low"))
AC.fc.dat$def_1 <- droplevels(AC.fc.dat$def_1)

Idents(AC.mb.dat) <- "SCT_snn_res.1"
Idents(AC.mb.dat, cells = WhichCells(object = AC.mb.dat, idents = c("10")))<-"SAA"
Idents(AC.mb.dat, cells = WhichCells(object = AC.mb.dat, idents = c("6","11","12")))<-"Gfap_high"
Idents(AC.mb.dat, cells = WhichCells(object = AC.mb.dat, idents = c("0","1","2","3","4","5","7","8","9")))<-"Gfap_low"
AC.mb.dat$def_1 <- Idents(AC.mb.dat)

Idents(AC.mb.dat) <- "def_1"
AC.mb.dat <- subset(AC.mb.dat,idents = c("SAA","Gfap_high","Gfap_low"))
AC.mb.dat$def_1 <- droplevels(AC.mb.dat$def_1)

#### 5.2 Calculate the cluster marker genes and inter-cluster DEGs ####
Idents(AC.fc.dat) <- "def_1"

FC.SvGh.diff <- FindMarkers(AC.fc.dat,ident.1 = "SAA",ident.2 = "Gfap_high",logfc.threshold = 0,test.use = "MAST")
FC.SvGl.diff <- FindMarkers(AC.fc.dat,ident.1 = "SAA",ident.2 = "Gfap_low",logfc.threshold = 0,test.use = "MAST")

SvGh.diff <- SvGh.diff %>% as_tibble(rownames = "Gene")
SvGl.diff <- SvGl.diff %>% as_tibble(rownames = "Gene")

FC.SvAll.dif <- FindMarkers(AC.fc.dat,ident.1 = "SAA",test.use = "MAST",logfc.threshold = 0)
FC.SvAll.diff <- as_tibble(FC.SvAll.dif,rownames = "Gene")
FC.SvAll.diff <- FC.SvAll.diff %>% filter(abs(avg_logFC)>0.3&p_val_adj<0.005) %>% arrange(rev(avg_logFC))

MB.SvAll.dif <- FindMarkers(AC.mb.dat,ident.1 = "SAA",group.by = "def_1",test.use = "MAST",logfc.threshold = 0)
MB.SvGh.dif <- FindMarkers(AC.mb.dat,ident.1 = "SAA",ident.2 = "Gfap_high",test.use = "MAST",logfc.threshold = 0)
MB.SvGl.dif <- FindMarkers(AC.mb.dat,ident.1 = "SAA",ident.2 = "Gfap_low",test.use = "MAST",logfc.threshold = 0)

MB.SvAll.diff <- as_tibble(MB.SvAll.dif,rownames = "Gene")
MB.SvAll.diff <- MB.SvAll.diff %>% filter(abs(avg_logFC)>0.3&p_val_adj<0.005) %>% arrange(rev(avg_logFC))

#### 5.3 Output the DEGs for pathway analysis (Gene Analytic platform) ####
FC.SvGh.diff = SvGh.diff %>% filter(p_val_adj<0.005 & abs(avg_logFC) >0.3)
FC.SvGl.diff = SvGl.diff %>% filter(p_val_adj<0.005 & abs(avg_logFC) >0.3)

write_csv(as_tibble(FC.SvGh.diff,rownames = "Gene"),path = "Pathway/Input/FC_SvGh.csv")
write_csv(as_tibble(FC.SvGl.diff,rownames = "Gene"),path = "Pathway/Input/FC_SvGl.csv")

MB.SvGh.diff <- MB.SvGh.dif %>% filter(p_val_adj<0.005 & abs(avg_logFC) >0.3)
MB.SvGl.diff <- MB.SvGl.dif %>% filter(p_val_adj<0.005 & abs(avg_logFC) >0.3)

write_csv(as_tibble(MB.SvGh.diff,rownames = "Gene"),path = "Pathway/Input/MB_SvGh.csv")
write_csv(as_tibble(MB.SvGl.diff,rownames = "Gene"),path = "Pathway/Input/MB_SvGl.csv")

#### 5.4 Merged astrocyte, only for show UMAP ####
ast.merge <- merge(x = tg.fc,y = list(wt.fc,tg.mb,wt.mb) ,project = "b5_merge")

fc.bc <- AC.fc.dat$def_1 %>% as.character() 
fc.bc <- paste0("FC_",fc.bc)
mb.bc <- AC.mb.dat$def_1 %>% as.character()
mb.bc <- paste0("MB_",mb.bc)
ast.merge.bc <- c(fc.bc,mb.bc)
names(ast.merge.bc) <- c(names(AC.fc.dat$def_1),names(AC.mb.dat$def_1))

ast.merge <- SubsetData(ast.merge,cells = names(ast.merge.bc))
ast.merge <- SCTransform(ast.merge, verbose = T)

fc.area <- rep("FC",dim(AC.fc.dat)[2])
mb.area <- rep("MB",dim(AC.mb.dat)[2])

ast.area <- c(fc.area,mb.area)

names(ast.area) <- c(colnames(AC.fc.dat),colnames(AC.mb.dat))
ast.merge$area <- ast.area
ast.merge$def_1 <-  ast.merge.bc

ast.merge <- RunPCA(ast.merge)
ast.merge <- RunUMAP(ast.merge,dims = 1:30,umap.method = "umap-learn")

#### 6. Microglia inter-region analysis ####

Idents(ast.fc) <- "Celltype"
MG.fc.dat <- subset(ast.fc,idents = "MG")

MG.fc.dat <- SCTransform(MG.fc.dat,verbose = T)

MG.fc.dat <- RunPCA(MG.fc.dat)
MG.fc.dat <- RunUMAP(MG.fc.dat,dims = 1:30,umap.method = "umap-learn")
MG.fc.dat <- FindNeighbors(MG.fc.dat,dims = 1:30)
MG.fc.dat <- FindClusters(MG.fc.dat,resolution = 0.8)

Idents(ast.mb) <- "Celltype"
MG.mb.dat <- subset(ast.mb,idents = "MG")

MG.mb.dat <- SCTransform(MG.mb.dat,verbose = T)

MG.mb.dat <- RunPCA(MG.mb.dat)
MG.mb.dat <- RunUMAP(MG.mb.dat,dims = 1:30,umap.method = "umap-learn")
MG.mb.dat <- FindNeighbors(MG.mb.dat,dims = 1:30)
MG.mb.dat <- FindClusters(MG.mb.dat,resolution = 0.8)

MG.fc.diff <- FindMarkers(MG.fc.dat,ident.1 = "tg.fc",ident.2 = "wt.fc",group.by = "orig.ident",logfc.threshold = 0,test.use = "MAST")
MG.mb.diff <- FindMarkers(MG.mb.dat,ident.1 = "tg.mb",ident.2 = "wt.mb",group.by = "orig.ident",logfc.threshold = 0,test.use = "MAST")

#### 6.1 Output the MG AC expression matrix and finish cell cell interaction alaysis in CellPhoneDB ####

Idents(ast.fc) <- "Celltype"
CCI.dat <- subset(ast.fc,idents = c("AC","MG"))
Idents(CCI.dat) <- "orig.ident"
CCI.dat <- subset(CCI.dat,idents = "tg.fc")

WT.bc <- CCI.dat$orig.ident %>% enframe(name = "Cellbc",value = "Strain")
AC.sub.idents <- AC.fc.dat$def_1 %>% enframe(name = "Cellbc",value = "Type") %>% filter(str_detect(Cellbc,"^tg"))
MG.idents <- MG.fc.dat$SCT_snn_res.0.6 %>% enframe(name = "Cellbc",value = "Type") %>% filter(str_detect(Cellbc,"^tg"))
CCI.dat <- SubsetData(CCI.dat,cells = c(AC.sub.idents$Cellbc,MG.idents$Cellbc))
MG.idents$Type <- "MG"

CCI.meta <- bind_rows(AC.sub.idents,MG.idents)
colnames(CCI.meta) <- c("Cell","cell_type")
CCI.dat <- CCI.dat@assays$SCT@data

allmmu.genes <- rownames(CCI.dat)
allhsa.genes <- lapply(allmmu.genes,mmu2hsa)
CCI.dat <- CCI.dat[!is.na(allhsa.genes),]
CCI.dat@Dimnames[[1]] <- unlist(allhsa.genes[!is.na(allhsa.genes)])

write.csv(CCI.dat, "~/bioinfo/Caro/Astro/Astro/DEGs/CellphoneDB/Input/FC_MGACtg.csv")
write.csv(CCI.meta,"~/bioinfo/Caro/Astro/Astro/DEGs/CellphoneDB/Input/FC_MGACtgmeta.csv")

Idents(ast.mb) <- "Celltype"
CCI.dat <- subset(ast.mb,idents = c("AC","MG"))
AC.sub.idents <- AC.mb.dat$def_1 %>% enframe(name = "Cellbc",value = "Type")
MG.idents <- MG.mb.dat$SCT_snn_res.1 %>% enframe(name = "Cellbc",value = "Type")
CCI.dat <- SubsetData(CCI.dat,cells = c(AC.sub.idents$Cellbc,MG.idents$Cellbc))
MG.idents$Type <- "MG"

CCI.meta <- bind_rows(AC.sub.idents,MG.idents)
colnames(CCI.meta) <- c("Cell","cell_type")
CCI.dat <- CCI.dat@assays$SCT@data


allmmu.genes <- rownames(CCI.dat)
allhsa.genes <- lapply(allmmu.genes,mmu2hsa)
CCI.dat <- CCI.dat[!is.na(allhsa.genes),]
CCI.dat@Dimnames[[1]] <- unlist(allhsa.genes[!is.na(allhsa.genes)])

write.csv(CCI.dat, "~/bioinfo/Caro/Astro/Astro/DEGs/CellphoneDB/Input/MB_MGAC.csv")
write.csv(CCI.meta,"~/bioinfo/Caro/Astro/Astro/DEGs/CellphoneDB/Input/MB_MGACmeta.csv")

