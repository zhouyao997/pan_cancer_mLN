seurat_qc@meta.data$cancer1 <- seurat_qc@meta.data$cancer
seurat_qc@meta.data$cancer1[which(seurat_qc@meta.data$cancer%in%c('BRCA1','BRCA2'))] <- 'BRCA'
seurat_qc@meta.data$cancer1[which(seurat_qc@meta.data$cancer%in%c('prostate'))] <- 'Prostate'
seurat_qc@meta.data$cancer1[which(seurat_qc@meta.data$cancer%in%c('Lung'))] <- 'LUNG'
yr <- c(
  "#fafac0", "#f5eca6", "#fee391", "#fec44f",
  "#fe9929", "#ec7014", "#cc4c02", "#8c2d04",
  "#611f03"
)

gr <- c(
  "#b9e9e0", "#7dc9b7", "#59b898", "#41ae76",
  "#16d355", "#238b45", "#116d37", "#025826",
  "#003516"
)

bl <- c(
  "#82cbf5", "#7ba7e0", "#5199cc", "#488dbe",
  "#3690c0", "#0570b0", "#0d71aa", "#045a8d",
  "#023858"
)

pur <- c(
  "#8c97c6", "#a28abd", "#997abd", "#9362cc",
  "#88419d", "#810f7c", "#4d004b"
)

bro <- c(
  "#8c510a", "#995401", "#be7816", "#be9430",
  "#ad8d36", "#a07540"
)

color1 <- c(bl[1:8], yr, pur, gr[c(3, 5, 6, 7, 9)], bro, "#A9A9A9")
color2 <- c(
  "#E0D39B", "#D05146", "#748EAE", "#3377a9",
  "#574F84", "#a35c4e","#a28abd","#997abd"
)
# main cell annotation ----------------------------------------------------

genes_to_check <- c('PTPRC',
  'CD3D','CD3E','CD3G',##T
  'NKG7','GNLY','KLRD1','KLRB1',###NK
  'MS4A1','CD79B',##B
  'CD68','CD14','TPSAB1' , 'TPSB2',##MYOLD
  'EPCAM','KRT8','KRT18','KRT19',##epi
  'DCN','COL1A1',##fib
  'VWF','PECAM1',###ENDO
  'IGHA1','IGHG1','JCHAIN'##PLASMA
  )
seurat_qc@meta.data$RNA_snn_res.0.3 <- paste0('c',seurat_qc@meta.data$RNA_snn_res.0.3)
seurat_qc@meta.data$RNA_snn_res.0.3 <- factor(seurat_qc@meta.data$RNA_snn_res.0.3,levels=paste0('c',0:20))

# annotation --------------------------------------------------------------

#####添加注释
B = paste0('c',c(4))
T =  paste0('c',c(1,13,14,20))
Myeloid=paste0('c',c(3,12,16))
NK=paste0('c',c(2))
Epithelium=paste0('c',c(0,7,8,9,11,15,17,19))
Endothelial=paste0('c',c(6))
Fibroblast=paste0('c',c(5,18))
Plasma=paste0('c',c(10))

#Doublet=c(15)

current.cluster.ids <- c(B,
                         T,
                         Myeloid,
                         NK,
                         Epithelium,
                         Endothelial,
                         Fibroblast,
                         Plasma
                         )

new.cluster.ids <- c(rep("B",length(B)),
                     rep("T",length(T)),
                     rep("Myeloid",length(Myeloid)),
                     rep("NK/NKT",length(NK)),
                     rep("Epithelium",length(Epithelium)),
                     rep("Endothelial",length(Endothelial)),
                     rep("Fibroblast",length(Fibroblast)),
                     rep("Plasma",length(Plasma))
                     )
                    

seurat_qc@meta.data$celltype <- plyr::mapvalues(x = as.character(seurat_qc@meta.data$RNA_snn_res.0.3), from = current.cluster.ids, to = new.cluster.ids)

# immune cell annotation --------------------------------------------------

gc()###清理内存
seurat_immune <- subset(seurat_qc, cells=row.names(seurat_qc@meta.data)
                                     [which(seurat_qc@meta.data$celltype%in%c('T','NK/NKT','B','Myeloid'))])

# immune <- list(immune = c('PTPRC'),
#                
#                    T =  c('CD3D','CD3E','CD3G'),##T
#                    CD4 =    c('CD4'),
#                CD8_Ex =c('CD8A','LAG3','PDCD1'),##CD8 Ex
#                CD8_EFF =c('IFNG','GZMK'),#CD8 EFF
#                Treg=c('FOXP3','IL2RA'),##Treg
#                Memory_T=c('CD44','IL7R','LTB'),##Memory T
#                Naive_T=c('CCR7','TCF7','LEF1','SELL'),###Naive T
#                   NK=c('NKG7','GNLY','KLRD1','KLRB1'),###NK
#                    B=c('CD79A','CD79B','MS4A1'),##B
#                #Myeloid=c('CD68','CD14','TPSAB1','TPSB2')
#                Mac=c("C1QA","C1QB","C1QC","CD68","CD163",'APOE'),
#                 #Mono_CD14 =c('CD14',"S100A8",'S100A9',"FCN1"),###经典
#                # Mono_CD16 = c('FCGR3A','LST1','LILRB2'),###非经典
#                Monocyte =c('CD14','FCGR3A',"S100A8",'S100A9',"FCN1",'LST1','LILRB2'),
#                Neutrophils =  c("FCGR3B",'CSF3R','FPR1',"CXCR2","SLC25A37"),
#                # cDC1 = c("CLEC9A",'FLT3','IDO1'),
#                # cDC2=c('CD1C','FCER1A','CLEC10A'),
#                # cDC3 =  c("LAMP3",'FSCN1'),
#                cDC=c('CD1C','FCER1A'),
#                Mast= c('TPSAB1','TPSB2')
# )

immune <- list(immune = c('PTPRC'),
               Epi=c('EPCAM'),
               T =  c('CD3D','CD3E','CD3G'),##T
               CD4 = c('CD4'),
               #T_h =c('IL23R','RORC','IL17A','FURIN'),
               CD8_Ex =c('PDCD1','CXCL13','LAYN'),##CD8 Ex
               CD8_EFF =c('CD8A','IFNG','GZMK'),#CD8 EFF
               Treg=c('FOXP3','IL2RA','IKZF2','TNFRSF9'),##Treg,IKZF2,CD4
               Memory_T=c('IL7R','GPR183','ZFP36L2','CXCR4','ZNF683'),##Memory T
               Naive_T=c('TCF7','CCR7','LEF1','SELL'),###Naive T
               NK=c('NKG7','GNLY','KLRD1','KLRB1'),###NK
               B=c('CD79A','CD79B','MS4A1'),##B
               #Myeloid=c('CD68','CD14','TPSAB1','TPSB2')
               Mac=c("C1QA","C1QB","C1QC","CD68","CD163",'APOE'),
               #Mono_CD14 =c('CD14',"S100A8",'S100A9',"FCN1"),###经典
               # Mono_CD16 = c('FCGR3A','LST1','LILRB2'),###非经典
               Monocyte =c('CD14','FCGR3A',"S100A8",'S100A9',"FCN1",'LST1','LILRB2'),
               Neutrophils =  c("FCGR3B",'CSF3R','FPR1',"CXCR2","SLC25A37"),
               # cDC1 = c("CLEC9A",'FLT3','IDO1'),
               # cDC2=c('CD1C','FCER1A','CLEC10A'),
               # cDC3 =  c("LAMP3",'FSCN1'),
               cDC=c('CD1C','FCER1A'),
               Mast= c('TPSAB1','TPSB2')
)
DotPlot(seurat_immune, features = immune, cluster.idents = TRUE, 
        assay='RNA' ,group.by = 'celltype_NEW' ) + theme_bw()+
  theme( panel.grid=element_blank(),axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5)) 
                                                                                                                            
DimPlot(
  seurat_immune,reduction = "tsne",pt.size = 1,group.by = "tissue"#label.color = color_cluster[unique(allseu@meta.data$cluster)]
  ,label=TRUE
)
T=c(0,1,10,12,17,20,21,25,26,3,6,7)
CD8_Ex =c()
CD8_EFF =c()
Treg=c()
Memory_T=c()
Naive_T=c()
NKT=c()
NK=c(8)
B=c(2,5)               
Mac=c(13,27,29,30,33,4)               
Monocyte =c(11)
Neutrophils = c(24)         
cDC=c(9)
Mast= c(15)
undefined=c(14,16,18,19,22,23,28,31,32)



current.cluster.ids <- c(#CD8_Ex, 
                         #CD8_EFF,
                         #Treg,
                         #Memory_T,
                         #Naive_T,
                         #NKT,
                         T,
                         NK,
                         B ,              
                         Mac,            
                         Monocyte, 
                         Neutrophils,    
                         cDC,
                         Mast,
                         undefined
                        )

new.cluster.ids <- c(rep("T",length(T)),
                     rep("NK",length(NK)),
                    rep("B",length(B)),
                    rep("Mac",length(Mac)),
                    rep("Monocyte",length(Monocyte)) ,
                   rep("Neutrophils",length(Neutrophils)),
                   rep("cDC",length(cDC)),
                    rep("Mast",length(Mast)) ,
                   rep("undefined",length(undefined))
                   )


seurat_immune@meta.data$celltype_NEW <- plyr::mapvalues(x = seurat_immune@meta.data$RNA_snn_res.1, from = current.cluster.ids, to = new.cluster.ids)

markers_immue <- FindAllMarkers(seurat_immune, only.pos = TRUE, min.pct = 0.25)

# T annotation ------------------------------------------------------------
seurat_T <- readRDS('/home/zhouyao/Pan_cancer/data/raw_data/RDS/last/seurat_T.Rds')
pos <- match(row.names(seurat_T@meta.data),row.names(seurat_immune@meta.data))
seurat_immune@meta.data$celltype_NEW[pos] <-  seurat_T@meta.data$celltype_T
seurat_T <- subset(seurat_immune, cells=row.names(seurat_immune@meta.data)
                        [which(seurat_immune@meta.data$celltype_NEW%in%c('T'))])
seurat_T %<>% NormalizeData(normalization.method = "LogNormalize") %>%
  FindVariableFeatures(selection.method = "vst")%>%
  ScaleData(vars.to.regress = c("nCount_RNA","percent.mt"))%>%
  RunPCA()%>%
  RunHarmony("patient", plot_convergence = TRUE)%>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = seq(0.1,1,0.1))%>%
  RunUMAP(reduction = "harmony", dims = 1:20)


immune_T <- list(immune = c('PTPRC'),
               Epi=c('EPCAM'),
               T =  c('CD3D','CD3E','CD3G'),##T
               CD4 =    c('CD4'),
               CD8_Ex =c('CD8A','LAG3','PDCD1'),##CD8 Ex
               CD8_EFF =c('IFNG','GZMK'),#CD8 EFF
               Treg=c('FOXP3','IL2RA'),##Treg
               Memory_T=c('CD44','IL7R','LTB'),##Memory T
               Naive_T=c('CCR7','TCF7','LEF1','SELL')###Naive T
)

DotPlot(seurat_T, features =immune_T,
        assay='RNA' ,group.by = 'RNA_snn_res.1' ) +theme_bw()+
  theme( panel.grid=element_blank(),axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5)) 

DimPlot(
  seurat_T,reduction = "tsne",pt.size = 1,group.by = "tissue"#label.color = color_cluster[unique(allseu@meta.data$cluster)]
  ,label=TRUE
)
saveRDS(seurat_T, paste0('/home/zhouyao/Pan_cancer/data/raw_data/RDS/last/', paste0('seurat_T',".Rds")))
saveRDS(seurat_immune, paste0('/home/zhouyao/Pan_cancer/data/raw_data/RDS/last/', paste0('seurat_immune',".Rds")))
saveRDS(seurat_qc, paste0('/home/zhouyao/Pan_cancer/data/raw_data/RDS/last/', paste0('seurat_qc',".Rds")))






# Myeloid annotation ------------------------------------------------------

seurat_Myeloid <- subset(seurat_immune, cells=row.names(seurat_immune@meta.data)
                   [which(seurat_immune@meta.data$celltype_NEW%in%c('Myeloid'))])
seurat_Myeloid %<>% NormalizeData(normalization.method = "LogNormalize") %>%
  FindVariableFeatures(selection.method = "vst")%>%
  ScaleData(vars.to.regress = c("nCount_RNA","percent.mt"))%>%
  RunPCA()%>%
  RunHarmony("patient", plot_convergence = TRUE)%>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = c(0.3,1))%>%
  RunUMAP(reduction = "harmony", dims = 1:20)


genes_to_myeloid = c(
'FCGR3B','FPR1',###Neutrophil
'S100A9', 'S100A8',# monocyte 单核细胞
'ITGAX','ITGAM','CD86','IL1B',###M1
'CD163','MRC1','MSR1',#M2
'CD1C','FCER1A','CLEC10A',#cDC
'LILRA4','CCDC50','IL3RA',##pDC
'TPSAB1' , 'TPSB2'# mast cells
)
genes_to_myeloid = c(
  'FCGR3B', 'CSF3R','FPR1',###Neutrophil,CSF3R,CD14-,9
  'CD14','FCGR3A',# monocyte 单核细胞
  'CD68','CD86','CCL4',###M1,'FCGR3A'+
  'CD163','MRC1','MSR1',#M2
  'CD1C','FCER1A',#cDC
  #'LILRA4',##pDC
  'TPSAB1' , 'TPSB2'# mast cells
)
genes_to_myeloid = c('CD3D','CD3E','GNLY',
  'EPCAM','KRT19','MS4A1','NKG7','IL7R',
  'FCGR3B', 'CSF3R',###Neutrophil,CSF3R,CD14-,9
  'CD14','LYZ','FCGR3A',"MS4A7",# monocyte 单核细胞0,1,
  'INHBA','NLRP3','C1QC','CD68','CD163',##mac
  'CD1C','FCER1A',#cDC
  #'LILRA4',##pDC
  'TPSAB1' , 'TPSB2'# mast cells
)
DotPlot(seurat_Myeloid, features = unique(genes_to_myeloid),
        group.by = "RNA_snn_res.1")  + coord_flip()
DimPlot(
  seurat_Myeloid,reduction = "umap",pt.size = 1,
  group.by = "RNA_snn_res.1"
  ,label=TRUE
)
DimPlot(
  seurat_Myeloid,reduction = "umap",pt.size = 1,
  group.by = "tissue"
  ,label=TRUE
)
 FeaturePlot(seurat_qc,c('LILRA4'))
# markers_myeloid <- FindAllMarkers(seurat_Myeloid, only.pos = TRUE, min.pct = 0.25)
# seurat_qc$celltype_1 <- seurat_qc$celltype
# cells <- row.names(seurat_Myeloid@meta.data)[which(seurat_Myeloid$seurat_clusters%in%c(5))]
# seurat_qc$celltype_1[which(row.names(seurat_qc@meta.data)%in%cells)] <- 'NK/NKT' 

th=theme( panel.grid=element_blank(),axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))  
myeloids = list(
  # NK <- c('CD3D','CD3E','GNLY','NKG7','IL7R'),
  # epi <- c('EPCAM','KRT19'),
  # B <- c('MS4A1'),
  Mac=c("C1QA","C1QB","C1QC","CD68","CD163",'APOE'),
  #   M1=c('CD68','CD86'),###M1,'FCGR3A'+
  #   M2=c('CD163','MRC1','MSR1'),#M2
  # Mac_INHBA=c('INHBA','IL1RN','CCL4'),
  # Mac_NLRP3=c('NLRP3','EREG','IL1B'),
  # Mac_LYVE1=c('LYVE1','PLTP','SEPP1'),
  # Mac_C1QC=c('C1QC','C1QA','APOE'),
  #mono=c('CD14','FCGR3A',"VCAN","FCN1","CD300E","S100A12"),
   Mono_CD14 =c('CD14',"S100A8",'S100A9',"FCN1"),###经典
   Mono_CD16 = c('FCGR3A','LST1','LILRB2'),###非经典
  neutrophils =  c("FCGR3B",'CSF3R','FPR1',"CXCR2","SLC25A37"),
  # pDC =c('LILRA4','GZMB','IL3RA'),
   cDC1 = c("CLEC9A",'FLT3','IDO1'),
   cDC2=c('CD1C','FCER1A','CLEC10A'),
   cDC3 =  c("LAMP3","CCR7",'FSCN1'),
   mast= c('TPSAB1','TPSB2')
  )

DotPlot(seurat_Myeloid, features = myeloids,
             assay='RNA' ,group.by = 'RNA_snn_res.1' )  +theme_bw()+theme( panel.grid=element_blank(),axis.text.x = element_text(angle = 45, 
                                                                                                                                 vjust = 0.5, hjust=0.5))
View(filter(markers_myeloid,cluster==24))

# T annotation ------------------------------------------------------------

# saveRDS(seurat_T,file='D:\\博士课题\\淋巴\\seurat_T.Rds')
# seurat_T <- readRDS('D:\\博士课题\\淋巴\\seurat_T.Rds')
# seurat_T %<>% NormalizeData(normalization.method = "LogNormalize") %>%
#   FindVariableFeatures(selection.method = "vst")%>%
#   ScaleData(vars.to.regress = c("nCount_RNA","percent.mt"))%>%
#   RunPCA()%>%
#   RunHarmony("patient", plot_convergence = TRUE)%>%
#   FindNeighbors(reduction = "harmony", dims = 1:30) %>%
#   FindClusters(resolution = c(0.3,1))%>%
#   RunUMAP(reduction = "harmony", dims = 1:20)
# 
# immune_T <- list(immune = c('PTPRC'),
#                  Epi=c('EPCAM'),
#                  
#                  T =  c('CD3D','CD3E','CD3G'),##T
#                  NK=c('NKG7','GNLY','KLRD1','KLRB1'),
#                  CD4 = c('CD4'),
#                  #T_h =c('IL23R','RORC','IL17A','FURIN'),
#                  CD8_Ex =c('CD8A','PDCD1','CXCL13','LAYN'),##CD8 Ex
#                  CD8_EFF =c('IFNG','GZMK'),#CD8 EFF
#                  Treg=c('FOXP3','IL2RA','IKZF2','TNFRSF9'),##Treg,IKZF2,CD4
#                  Memory_T=c('IL7R','GPR183','ZFP36L2','CXCR4','ZNF683'),##Memory T
#                  Naive_T=c('TCF7','CCR7','LEF1','SELL')###Naive T
# )
# #FeaturePlot(seurat_T,c('CD4','FOXP3','IL2RA','LTB'))
# DotPlot(seurat_T, features =immune_T,
#         assay='RNA' ,group.by = 'RNA_snn_res.0.3' ) +theme_bw()+
#   theme(panel.grid=element_blank(),axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5)) 
# DimPlot(
#   seurat_T,reduction = "umap",pt.size = 1,group.by = "RNA_snn_res.0.3"#label.color = color_cluster[unique(allseu@meta.data$cluster)]
#   ,label=TRUE
# )
# DimPlot(
#   seurat_T,reduction = "umap",pt.size = 1,group.by = "tissue"#label.color = color_cluster[unique(allseu@meta.data$cluster)]
#   ,label=TRUE
# )
# FeaturePlot(seurat_T,c('CD4','FOXP3','IL2RA','IL7R'))                                                                                                                                    
# T_ex =c(8)##CD8 Ex
# CD8_EFF =c(1)#CD8 EFF
# T_reg=c(3)##Treg
# Memory_T=c(10,17,2,4,7,9)##Memory T
# Naive_T=c(6)###Naive T
# NKT=c(5)
# undefined=c(12,16,0,11,13,14,15)
# current.cluster.ids <- c(T_ex,##CD8 Ex
#                          CD8_EFF,#CD8 EFF
#                          T_reg,##Treg
#                          Memory_T,##Memory T
#                          Naive_T,###Naive T
#                          NKT,
#                          undefined
# )
# 
# new.cluster.ids <- c(rep("T_ex",length(T_ex)),
#                      rep("CD8_EFF",length(CD8_EFF)),
#                      rep("T_reg",length(T_reg)),
#                      rep("Memory_T",length(Memory_T)),
#                      rep("Naive_T",length(Naive_T)),
#                      rep("NKT",length(NKT)),
#                      rep("undefined",length(undefined))
# )
# 
# #FeaturePlot(seurat_T,c('CD4','FOXP3','IL2RA','LTB'))
# seurat_T@meta.data$celltype_T <- plyr::mapvalues(x = seurat_T@meta.data$RNA_snn_res.0.3, from = current.cluster.ids, to = new.cluster.ids)
# Idents(seurat_T) <- seurat_T$RNA_snn_res.0.3
# markers_T <- FindAllMarkers(seurat_T, only.pos = TRUE, min.pct = 0.25)
# View(filter(markers_T,cluster==0))
# 
