# This script is aiming to reproduce the paper's figure (only sequencing data sourced plot)
# Some figures were undergo post-process in Adobe AI

library(ggpubr)
library(ggpmisc)
library(ggplot2)
library(viridisLite)
library(ggrepel)

#### Figure 1 ####
# Figure 1f and 1g

AC.fc.dat$def_1 <- factor(AC.fc.dat$def_1,levels = c("SAA", "Gfap_high", "Gfap_low"  ))
ast.fc.less <- subset(ast.fc,cells = sample(colnames(ast.fc),6000))
DimPlot(ast.fc.less,group.by = "orig.ident",label = F, pt.size = 0.1,cols = c("#bf2928","#4db8c2"))+ NoLegend()

AC.mb.dat$def_1 <- factor(AC.mb.dat$def_1,levels = c("SAA", "Gfap_high", "Gfap_low"  ))
ast.mb.less <- subset(ast.mb,cells = sample(colnames(ast.mb),6000))
DimPlot(ast.mb.less,group.by = "orig.ident",label = F, pt.size = 0.1,cols = c("#bf2928","#4db8c2"))

ast.fc.less$Celltype <- droplevels(ast.fc.less$Celltype)
ast.fc.less$Celltype <- factor(ast.fc.less$Celltype,levels = levels(ast.fc.less$Celltype))
DimPlot(ast.fc.less,group.by = "Celltype",label = T, pt.size = 0.1,label.size = 6)

ast.mb.less$Celltype <- droplevels(ast.mb.less$Celltype)
ast.mb.less$Celltype <- factor(ast.mb.less$Celltype,levels = c(levels(ast.fc.less$Celltype),"TNC"))
DimPlot(ast.mb.less,group.by = "Celltype",label = T, pt.size = 0.1,label.size = 6)

#### Figure 2 ####
# Figure 2a~2e

DimPlot(ast.merge,group.by = "ori",label = F, pt.size = 0.4)+ 
  ggplot2::coord_cartesian(xlim=c(-4,12 ),ylim=c(0,15))+ NoLegend()

FeaturePlot(ast.merge,features = "Gfap",pt.size = 0.4,max.cutoff = 'q90')+
  ggplot2::coord_cartesian(xlim=c(-4,12 ),ylim=c(0,15))+
  theme(legend.position = "none", plot.title = element_blank())

DimPlot(ast.merge,group.by = "orig.ident",label = F, pt.size = 0.4)+
  ggplot2::coord_cartesian(xlim=c(-4,12 ),ylim=c(0,15))+ NoLegend()

DimPlot(ast.merge,group.by = "def_1",label = F, pt.size = 0.4)+
  ggplot2::coord_cartesian(xlim=c(-4,12 ),ylim=c(0,15))+ NoLegend()

# Figure 2h

greendot.gene2 <- read_csv("~/bioinfo/Caro/Astro/seq.csv") # A marker gene list from published paper 

green_dot_ac2(greendot.gene2$gene3)

green_dot_ac2 <- function(tar.genes){
  # Prepare the expression data from FC
  greendot.fcdat <- subset(AC.fc.dat,features = unique(tar.genes))
  
  # Prepare the expression data from MB
  greendot.mbdat <- subset(AC.mb.dat,features = unique(tar.genes))
  
  # Combine two expression matrix and normalzie
  greendot.dat <- merge(greendot.fcdat, greendot.mbdat,merge.data = T)
  greendot.dat <- SCTransform(greendot.dat)
  greendot.mat <- greendot.dat@assays$SCT@scale.data
  greendot.mat <- as_tibble(t(as.matrix(greendot.mat)),rownames="Cellbc")
  
  
  greendot.fcclust <- greendot.fcdat$def_1
  greendot.fcclust <- recode(greendot.fcclust, SAA = "SAA_FC",Gfap_high = "Gfap_high_FC", Gfap_low = "Gfap_low_FC")
  greendot.fcclust <- enframe(greendot.fcclust,name = "Cellbc",value = "Clusters")
  
  greendot.mbclust <- greendot.mbdat$def_1
  greendot.mbclust <- recode(greendot.mbclust, SAA = "SAA_MB",Gfap_high = "Gfap_high_MB", Gfap_low = "Gfap_low_MB")
  greendot.mbclust <- enframe(greendot.mbclust,name = "Cellbc",value = "Clusters")
  
  greendot.clust <- rbind(greendot.fcclust,greendot.mbclust)
  greendot.mat <- left_join(greendot.clust,greendot.mat,by="Cellbc")
  
  greendot.pct <- greendot.mat %>% group_by(Clusters) %>% summarise_all(funs(count.pct))
  greendot.pct <- greendot.pct[,-2]
  greendot.pct <- greendot.pct %>% gather(key = "Gene", value = "pct",-Clusters)
  
  greendot.mat <- greendot.mat %>% group_by(Clusters) %>% summarise_all(funs(mean))
  greendot.mat <- greendot.mat[,-2]
  greendot.mat <- greendot.mat %>% gather(key = "Gene", value = "avg_expr",-Clusters)
  
  greendot.dat <- inner_join(greendot.pct,greendot.mat,by=c("Clusters","Gene"))
  greendot.dat$Gene <- factor(greendot.dat$Gene,levels = unique(tar.genes))
  greendot.dat$Clusters <- factor(greendot.dat$Clusters,levels = rev(levels(greendot.dat$Clusters)))
  
  
  greendot.dat$avg_expr <- ifelse(greendot.dat$avg_expr>4,4,greendot.dat$avg_expr)
  greendot.dat$avg_expr <- ifelse(greendot.dat$avg_expr< -4, -4,greendot.dat$avg_expr)
  
  ggplot(greendot.dat, aes(x=Gene, y=Clusters)) + 
    geom_point(aes(color=avg_expr,size=pct))+scale_size_continuous(range = c(0,8))+
    theme(axis.text.x = element_text(angle = 45,hjust=1))+scale_color_viridis(option = "D")+
    expand_limits(size = seq(0.05, 1.1, by = 0.25))
  
}

# Figure 2I Comparision between DAA and 

# Prepare the DEGlist first
DAA.diff <- read_csv("~/bioinfo/Caro/Astro/Out_source_data/DAA_DEGs.csv") # From published paper
DAAvsGl <- DAA.diff %>% filter(`Cluster ID #1`==1 & `Cluster ID #2` ==4)
DAAvsGh <- DAA.diff %>% filter(`Cluster ID #1`==4 & `Cluster ID #2` ==6)

sig.list <- read_csv("~/bioinfo/Caro/Astro/Out_source_data/Signature_top50.csv")
MCAO.list <- sig.list %>% filter(Source=="MACO") 
LPS.list <- sig.list %>% filter(Source=="LPS")

Idents(AC.fc.dat) <- "def_1"
SvAll.dif.fc <- FindMarkers(AC.fc.dat,ident.1 = "SAA",test.use = "MAST",logfc.threshold = 0)
SvAll.diff.fc <- as_tibble(SvAll.dif.fc,rownames = "Gene")
SvAll.diff.fc <- SvAll.diff.fc %>% filter(abs(avg_logFC)>0.3&p_val_adj<0.005) %>% arrange(rev(avg_logFC))

Idents(AC.mb.dat) <- "def_1"
SvAll.dif.mb <- FindMarkers(AC.mb.dat,ident.1 = "SAA",test.use = "MAST",logfc.threshold = 0)
SvAll.dif.mb <- arrange(SvAll.dif.mb,desc(abs(avg_logFC)))

FC.DAA.cpr.list <- list(
  #SvGl = SvGl.diff %>% filter(p_val_adj<0.005 & abs(avg_logFC) >0.3) %>% pull(Gene),
  SvAll = SvAll.diff.fc %>% pull(Gene),
  DAAvGl = DAAvsGl %>% pull(Gene),
  DAAvsGh = DAAvsGh %>%  pull(Gene),
  MCAO = MCAO.list %>% pull(`Gene symbol`),
  LPS = LPS.list %>% pull(`Gene symbol`)
)
venn(FC.DAA.cpr.list, zcolor = c("#bf2928","#4d089a","#323edd","#dc2ade","#E91E63") ,cexil = 1, cexsn = 0.8,borders = T,opacity = 0.6) #DEG of SvAll FC overlapped with others 

MB.SvAll.diff <- as_tibble(MB.SvAll.dif,rownames = "Gene") %>% filter(abs(avg_logFC)>0.3&p_val_adj<0.005) %>% arrange(rev(avg_logFC))

MB.DAA.cpr.list <- list(
  #SvGl = SvGl.diff %>% filter(p_val_adj<0.005 & abs(avg_logFC) >0.3) %>% pull(Gene),
  SvAll = MB.SvAll.diff %>% pull(Gene),
  DAAvGl = DAAvsGl %>% pull(Gene),
  DAAvsGh = DAAvsGh %>%  pull(Gene),
  MCAO = MCAO.list %>% pull(`Gene symbol`),
  LPS = LPS.list %>% pull(`Gene symbol`)
)

venn(MB.DAA.cpr.list, zcolor = c("#bf2928","#4d089a","#323edd","#dc2ade","#E91E63") ,cexil = 1, cexsn = 0.8,borders = T,opacity = 0.6) #DEG of SvAll MB overlapped with others MB.SvAll.dif

#Figure 2J Volcano plot 
inter20 <- intersect(rownames(SvAll.dif)[1:55],rownames(SvAll.dif.mb)[1:55])

vol.genes <- c(inter20[1:20], greendot.gene2$gene3[1:20]) %>% unique()

vol.genes <- c(vol.genes,"Slc6a9","Slc6a11","Agt","Bdnf","Plat")

SvAll.dif <- arrange(SvAll.dif,desc(abs(avg_logFC)))

EnhancedVolcano(toptable = SvAll.dif, #FC volcano plot
                lab = rownames(SvAll.dif),
                selectLab = vol.genes,
                x = 'avg_logFC',
                y = 'p_val_adj',
                pCutoff = 0.005,
                FCcutoff = 0.3,
                drawConnectors = T,
                arrowheads = FALSE)

EnhancedVolcano(toptable = SvAll.dif.mb,
                lab = rownames(SvAll.dif.mb),
                selectLab = vol.genes,
                x = 'avg_logFC',
                y = 'p_val_adj',
                pCutoff = 0.005,
                FCcutoff = 0.3,
                drawConnectors = T,
                arrowheads = FALSE)

#### Figure 3 ####

# Figure 3a, 3b
DimPlot(MG.fc.dat,group.by = "orig.ident",label = F, pt.size = 0.4,cols = c("#bf2928","#4db8c2"))+
  ggplot2::coord_cartesian(xlim=c(-11, 6),ylim=c(-5,9))

DimPlot(MG.mb.dat,group.by = "orig.ident",label = F, pt.size = 0.5,cols = c("#bf2928","#4db8c2"))

# Figure 3a, 3b
EnhancedVolcano(toptable = MG.fc.diff,
                lab = rownames(MG.fc.diff),
                selectLab = c(rownames(MG.fc.diff)[which(abs(MG.fc.diff$avg_logFC)>0.3&MG.fc.diff$p_val_adj<0.005)],
                              "Slc6a9","Slc6a11","Agt"),
                x = 'avg_logFC',
                y = 'p_val_adj',
                pCutoff = 0.005,
                FCcutoff = 0.3,
                drawConnectors = T,
                arrowheads = FALSE)

EnhancedVolcano(toptable = MG.mb.diff,
                lab = rownames(MG.mb.diff),
                selectLab = rownames(MG.mb.diff)[which(abs(MG.mb.diff$avg_logFC)>0.3&MG.mb.diff$p_val_adj<0.005)],
                x = 'avg_logFC',
                y = 'p_val_adj',
                pCutoff = 0.005,
                FCcutoff = 0.3,
                drawConnectors = T,
                arrowheads = FALSE)

# Figure 3e

library(devtools)
source_gist("524eade46135f6348140")

MG.lable.add <- read_csv("../MG_DEG/MG_DEG.csv")

MG.diff.cpr <- MG.diff.cpr %>% mutate(FCweight=(avg_logFC_FC)+(avg_logFC_MB)) 

MG.lable.top5 <-top_n(MG.diff.cpr,5,FCweight) %>% select(Gene)
MG.lable.bottom5 <- top_n(MG.diff.cpr,-5,FCweight) %>% select(Gene)

MG.label.upright <- top_n(MG.diff.cpr,15,avg_logFC_MB) %>% select(Gene)

MG.diff.cpr <- MG.diff.cpr %>% mutate(Label2 = ifelse(Gene %in% c(MG.lable.add$MG_gene,MG.lable.top5$Gene,
                                                                  MG.lable.bottom5$Gene,MG.label.upright$Gene),Gene,""))

ggplot(MG.diff.cpr, aes(x=avg_logFC_FC,y=avg_logFC_MB))+ geom_point(size=3,alpha = 0.8,aes(colour=Type2))+
  stat_smooth_func(geom="text",method="lm",hjust=0,vjust=2,parse=T,ypos = 1.0,xpos = -0.75)+
  stat_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8,fullrange = T)+theme_bw()+
  geom_text_repel(aes(label=Label2),min.segment.length=0.1)+
  stat_fit_glance(method = "lm", 
                  method.args = list(formula = fml),
                  label.x = "right",
                  label.y = "bottom",
                  aes(label = paste("italic(P)*\"-value = \"*", 
                                    signif(..p.value.., digits = 4), sep = "")),
                  parse = TRUE)+
  scale_color_manual(values=c("#CD534CFF","#868686FF"))+
  geom_hline(yintercept = 0) +geom_vline(xintercept = 0)

#### Figure 5 ####
# Figure 5a

DimPlot(kar.ast,label = F,group.by = "def_mo", reduction = "harmony")+theme(legend.position = "none")+ ggtitle("")+
  scale_y_continuous(limits = c(5.689721, 16.175358),breaks=seq(7.5,15,2.5))+
  scale_x_continuous(limits = c(0.7222885, 12.4235525),breaks=seq(2.5,12.5,2.5))

Idents(kar.ast) <-  "Status"
DimPlot(kar.ast,label = F,group.by = "Status", reduction = "harmony")+theme(legend.position = "none")+ ggtitle("")+
  scale_y_continuous(limits = c(5.689721, 16.175358),breaks=seq(7.5,15,2.5))+
  scale_x_continuous(limits = c(0.7222885, 12.4235525),breaks=seq(2.5,12.5,2.5))

# Figure 5b and 5d

greendot.gene2 <- read_csv("~/bioinfo/Caro/Astro/seq.csv") # Load the genes appeared in figure 2h
engene.name <- mouse2human(greendot.gene2$gene3)
engene.name <- engene.name[engene.name$humanGene %in% rownames(wang.ast),]

pdf("./Kamath/Output/F5_Vlnplot_grouped_by_SAA.pdf")
for (i in 1:length(engene.name$humanGene)) {
  en.name <- engene.name$humanGene[i]
  print(VlnPlot(object = kar.ast, features = en.name,pt.size = 0)+ggtitle(en.name))
}
dev.off()

# Figure 5c

kar.ast$SAAscore <- predictions$prediction.score.SAA

FeaturePlot(kar.ast, features = "SAAscore",reduction = "harmony") +
  scale_colour_continuous(limits = c(0,1))+
  ggtitle("") & scale_color_viridis_c() 

Idents(kar.ast) <- factor(kar.ast$def_mo,levels = c("SAA","Gfap_high","Gfap_low"))
VlnPlot(kar.ast,features = "SAAscore",pt.size = 0) +  NoLegend() +ggtitle("")

# Figure 5e

FeaturePlot(kar.ast, features = "SLC1A2",reduction = "harmony") + NoLegend() + ggtitle("") & scale_color_viridis_c() 
FeaturePlot(kar.ast, features = "AGT",reduction = "harmony")  & scale_color_viridis_c() 

VlnPlot(kar.ast,features = "SLC1A2",pt.size = 0 ) + NoLegend() + ggtitle("")
VlnPlot(kar.ast,features = "AGT",pt.size = 0 ) + NoLegend() + ggtitle("")


##### Supplyment figure 1####

# Sf 1c and 1d
VlnPlot(object = ast.fc, features = c("Kcnj8"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.fc, features = c("Acta2"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.fc, features = c("Cldn5"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.fc, features = c("Ctss"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.fc, features = c("Ntsr2"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.fc, features = c("Syt1"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.fc, features = c("Pdgfra"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.fc, features = c("Cldn11"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.fc, features = c("Sox11"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.fc, features = c("Ccdc153"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.fc, features = c("Ttr"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.fc, features = c("Alas2"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.fc, features = c("Pf4"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.fc, features = c("Plac8"),pt.size = 0,ncol = 1,y.max = 8)

VlnPlot(object = ast.mb, features = c("Kcnj8"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Acta2"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Cldn5"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Ctss"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Ntsr2"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Syt1"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Pdgfra"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Cldn11"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Sox11"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Ccdc153"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Ttr"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Alas2"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Pf4"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Rax"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ast.mb, features = c("Plac8"),pt.size = 0,ncol = 1,y.max = 8)

##### Supplyment figure 2####

#Sf 2a and 2b

DimPlot(ast.merge,group.by = "cluster",label = F, pt.size = 0.4)+
  ggplot2::coord_cartesian(xlim=c(-4,12 ),ylim=c(0,15))+ NoLegend()

#Sf 2c

VlnPlot(ast.merge,features = "Gfap",idents = c(paste0("FC_",0:8)),pt.size = 0,sort = T)
VlnPlot(ast.merge,features = "Gfap",idents = c(paste0("MB_",0:12)),pt.size = 0,sort = T)

#Sf 2d

VlnPlot(AC.fc.dat,features = "Gfap")
VlnPlot(AC.mb.dat,features = "Gfap")

##### Supplyment figure 3####

#Sf 3a and 3b

AC.fc.dat$def_1 <- factor(AC.fc.dat$def_1,levels = c("SAA","Gfap_high","Gfap_low"))
Idents(AC.fc.dat) <- "def_1"
AC.fc.top20 <- FindAllMarkers(AC.fc.dat,test.use = "MAST")
AC.fc.top20 <- AC.fc.top20 %>% as_tibble(rownames = "Gene") %>% group_by(cluster) %>% top_n(20,avg_logFC)
DoHeatmap(AC.fc.dat,features = rev(AC.fc.top20$gene))

AC.mb.dat$def_1 <- factor(AC.mb.dat$def_1,levels = c("SAA","Gfap_high","Gfap_low"))
Idents(AC.mb.dat) <- "def_1"
AC.mb.top20 <- FindAllMarkers(AC.mb.dat,test.use = "MAST")
AC.mb.top20 <- AC.mb.top20 %>% as_tibble(rownames = "Gene") %>% group_by(cluster) %>% top_n(20,avg_logFC)
DoHeatmap(AC.mb.dat,features = rev(AC.mb.top20$gene))

##### Supplyment figure 8####

#Sf 8a, 8b and 8c

Idents(kar.ast) <- "Status"
DimPlot(kar.ast,cells = WhichCells(kar.ast,idents = "PD"),group.by = "def_mo",reduction = "harmony")+theme(legend.position = "none")+ ggtitle("")+
  scale_x_continuous(limits = c(-22, 21),breaks=seq(-20,20,10))+
  scale_y_continuous(limits = c(-20, 12),breaks=seq(-10,10,10))

DimPlot(kar.ast,cells = WhichCells(kar.ast,idents = "LBD"),group.by = "def_mo",reduction = "harmony")+theme(legend.position = "none")+ ggtitle("")+
  scale_x_continuous(limits = c(-22, 21),breaks=seq(-20,20,10))+
  scale_y_continuous(limits = c(-20, 12),breaks=seq(-10,10,10))

DimPlot(kar.ast,cells = WhichCells(kar.ast,idents = "Ctrl"),group.by = "def_mo",reduction = "harmony")+theme(legend.position = "none")+ ggtitle("")+
  scale_x_continuous(limits = c(-22, 21),breaks=seq(-20,20,10))+
  scale_y_continuous(limits = c(-20, 12),breaks=seq(-10,10,10))

#Sf 8d

Idents(kar.ast) <- "def_mo"
kar.acmrks <- FindAllMarkers(kar.ast,logfc.threshold = 0.05)
kar.acmrks2 <- kar.acmrks %>% filter(p_val_adj <0.05) %>% group_by(cluster) %>% top_n(50,avg_logFC)

saa.heatmap.gene.top20 <- kar.acmrks2 %>% group_by(cluster) %>% top_n(20,avg_logFC) %>% arrange(desc(avg_logFC))
saa.heatmap.gene.top20 <- saa.heatmap.gene.top20 %>% group_by(cluster) %>% arrange(desc(avg_logFC),.by_group=T)
Idents(kar.ast) <- factor(kar.ast$def_mo,levels = c("SAA","Gfap_high","Gfap_low"))

DoHeatmap(kar.ast,cells = sample(colnames(kar.ast),15000),features = (saa.heatmap.gene.top20$gene))

##### Supplyment figure 9####

#Sf 9a and 9c

DimPlot(ac.prada,label = F)+theme(legend.position = "none")+ ggtitle("")+
  scale_y_continuous(limits = c(5.689721, 16.175358),breaks=seq(7.5,15,2.5))+
  scale_x_continuous(limits = c(0.7222885, 12.4235525),breaks=seq(2.5,12.5,2.5))

Idents(ac.prada) <- "PD"
DimPlot(ac.prada,label = F) +theme(legend.position = "none")+ ggtitle("")+
  scale_y_continuous(limits = c(5.689721, 16.175358),breaks=seq(7.5,15,2.5))+
  scale_x_continuous(limits = c(0.7222885, 12.4235525),breaks=seq(2.5,12.5,2.5))

FeaturePlot(ac.prada, features = "SAAscore") + NoLegend()+ expand_limits(colour = seq(0, 1, by = 0.25)) +
  ggtitle("") + expand_limits(colour = seq(0, 1, by = 0.25))+
  scale_y_continuous(limits = c(5.689721, 16.175358),breaks=seq(7.5,15,2.5))+
  scale_x_continuous(limits = c(0.7222885, 12.4235525),breaks=seq(2.5,12.5,2.5)) & scale_color_viridis_c() 

VlnPlot(ac.prada, features = "SAAscore",group.by = "def_mo",pt.size = 0)

#Sf 9b, 9d and 9g

Idents(ac.prada) <- factor(ac.prada$def_mo,levels = c("SAA","Gfap_high","Gfap_low"))

pdf("../Parada/Figures/Vlnplot_SAAmkrs_human_max4.pdf")
for (i in 1:length(engene.name3)) {
  en.name <- engene.name3[i]
  syb.name <- engene.name$external_gene_name[match(en.name, engene.name$ensembl_gene_id)]
  print(VlnPlot(object = ac.prada, features = en.name ,pt.size = 0,y.max = 4)+ggtitle(syb.name))
}
dev.off()

#Sf 9f

Idents(ac.prada) <- "PD"
DimPlot(ac.prada,cells = WhichCells(ac.prada,idents = "PD"),group.by = "def_mo")+theme(legend.position = "none")+ ggtitle("")+
  scale_y_continuous(limits = c(5.689721, 16.175358),breaks=seq(7.5,15,2.5))+
  scale_x_continuous(limits = c(0.7222885, 12.4235525),breaks=seq(2.5,12.5,2.5))

DimPlot(ac.prada,cells = WhichCells(ac.prada,idents = "Control"),group.by = "def_mo")+theme(legend.position = "none")+ ggtitle("")+
  scale_y_continuous(limits = c(5.689721, 16.175358),breaks=seq(7.5,15,2.5))+
  scale_x_continuous(limits = c(0.7222885, 12.4235525),breaks=seq(2.5,12.5,2.5))

#Sf 9g

FeaturePlot(ac.prada, features = "ENSG00000110436" ) +ggtitle("") + NoLegend()+
  scale_y_continuous(limits = c(5.689721, 16.175358),breaks=seq(7.5,15,2.5))+
  scale_x_continuous(limits = c(0.7222885, 12.4235525),breaks=seq(2.5,12.5,2.5)) & scale_color_viridis_c() 

FeaturePlot(ac.prada, features = "ENSG00000135744" ) +ggtitle("") + NoLegend()+
  scale_y_continuous(limits = c(5.689721, 16.175358),breaks=seq(7.5,15,2.5))+
  scale_x_continuous(limits = c(0.7222885, 12.4235525),breaks=seq(2.5,12.5,2.5))  & scale_color_viridis_c() 

#Sf 9e

hu.acmrks2 <- hu.acmrks %>% filter(p_val_adj <0.05) %>% group_by(cluster) %>% top_n(50,avg_logFC)
hu.acmrks2$Symbol <- annotations$external_gene_name[match(hu.acmrks2$gene,annotations$ensembl_gene_id)]

saa.heatmap.gene.top20 <- hu.acmrks2 %>% group_by(cluster) %>% top_n(20,avg_logFC) %>% arrange(desc(avg_logFC))
saa.heatmap.gene.top20 <- saa.heatmap.gene.top20[!duplicated(saa.heatmap.gene.top20$Symbol),]
saa.heatmap.gene.top20 <- saa.heatmap.gene.top20 %>% group_by(cluster) %>% arrange(desc(avg_logFC),.by_group=T)
Idents(ac.prada) <- factor(ac.prada$def_mo,levels = c("SAA","Gfap_high","Gfap_low"))

DoHeatmap(ac.prada,features = (saa.heatmap.gene.top20$gene),group.by = "def_mo")

