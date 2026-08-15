###############################################BRCA-COMMOT#####################################
tissue_positions<-read.csv(file = "/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-COMMOT/outs/spatial/tissue_positions_list.csv",row.names = 1,header = F)
colSums(tissue_positions)

library(STdeconvolve)
library(SpatialExperiment)

set.seed(2022)
matrix_dir="~/Documents/BulkToST/Dataset/BRCA-COMMOT/outs"
#matrix_dir = "~/Downloads/A1_spaceranger_output/outs/"
#barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
#features.path <- paste0(matrix_dir, "features.tsv.gz")
#matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

se <- SpatialExperiment::read10xVisium(samples = matrix_dir,
                                       type = "HDF5",
                                       data = "filtered",load = T)

se
##Change gene names from ENSEMBLE IDs to Symbols:
#se@assays@data@listData$counts@seed@dimnames[[1]]<-se@rowRanges@elementMetadata@listData$symbol
## this is the genes x barcode sparse count matrix
cd <- se@assays@data@listData$counts

### Convert the TENxMatrix to dgCMatrix so that the cleanCounts() works correctly
names_dims<-cd@seed@dimnames
cd<-Matrix::Matrix(cd, sparse = TRUE,nrow = 36601, ncol = 3798, dimnames = names_dims)


pos <- SpatialExperiment::spatialCoords(se)

## change column names to x and y
## for this dataset, we will visualize barcodes using "pxl_col_in_fullres" = "y" coordinates, and "pxl_row_in_fullres" = "x" coordinates
colnames(pos) <- c("y", "x")

cd<-as.matrix(cd)
##Change gene names from ENSEMBLE IDs to Symbols:
rownames(cd)<-se@rowRanges@elementMetadata@listData$symbol
cd<-t(cd) ### To be Spots X Genes

pos_cor<-tissue_positions[intersect(rownames(pos),rownames(tissue_positions)),2:3]
colnames(pos_cor) <- c("y", "x")

write.csv(cd,file="/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-COMMOT/processed/STcount.csv",quote = F,row.names = T)
write.csv(pos_cor,file = "/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-COMMOT/processed/Pos.csv",quote = F,row.names = T)


tissue_mask<-matrix(0,nrow(pos_cor),1)
for(i in 1:nrow(pos_cor)){
  tissue_mask[i,1]<-c("True")
}
write.csv(tissue_mask,file="/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-COMMOT/processed/tissue_mask.csv",quote = F,row.names = F)
write.csv(colnames(cd),file="/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-COMMOT/processed/Gene_names.csv",quote = F,row.names = F)
write.csv(rownames(cd),file="/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-COMMOT/processed/Barcodes.csv",quote = F,row.names = F)


#######################Bulk-BRCA-TCGA (Downloaded from COMICS paper repository)
GSE1456<-read.csv(file="/Users/naminiyakan/Documents/BulkToST/Dataset/Bulk-BRCA-TCGA/GSE1456.csv",header = T)

library(GEOquery)
gset <- getGEO("GSE1456", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

meta_data<-cbind(as.numeric(gset$`RELAPSE:ch1`),as.numeric(gset$`SURV_RELAPSE:ch1`),as.numeric(gset$`DEATH:ch1`),as.numeric(gset$`DEATH_BC:ch1`),as.numeric(gset$`SURV_DEATH:ch1`),gset$`SUBTYPE:ch1`,as.numeric(gset$`ELSTON:ch1`))
colnames(meta_data)<-c("Relapse","Surv_relapse","Death","Death_BC","Surv_death","Subtype","Elston")

###### Total Death Occurence 
colSums(as.matrix(as.numeric(meta_data[,3])))


###### Total Death Occurence Due to BRCA
colSums(as.matrix(as.numeric(meta_data[,4])))


#######Network_Construction
neighbors<-read.csv(file = "/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-COMMOT/processed/neighbors.csv",header = T,row.names = 1)
#tot_ST<-read.csv(file = "/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/Total/tot_ST.csv",header = T,row.names = 1)

sim_graph<-matrix(0,nrow(cd),nrow(cd))
colnames(sim_graph)<-rownames(cd)
rownames(sim_graph)<-rownames(cd)
for(i in 1:nrow(neighbors)){
  sim_graph[(neighbors$X0[i])+1,(neighbors$X1[i])+1]<-1
  sim_graph[(neighbors$X1[i])+1,(neighbors$X0[i])+1]<-1
}

library(dplyr)
# Remove Duplicate columns
cd <- cd[, !colnames(cd) %in% c("ARMCX5-GPRASP2","CYB561D2","GGT1","GOLGA8M","HSPA14","LINC01238",
                                "LINC01505","MATR3","TBCE","TMSB15B")]

library(Seurat)
sc_exprs <- as.matrix(t(cd))
Seurat_tmp <- CreateSeuratObject(counts=t(cd))
Seurat_tmp <- NormalizeData(object = Seurat_tmp,normalization.method = "LogNormalize",scale.factor = 10000,verbose = T)
Seurat_tmp <- FindVariableFeatures(Seurat_tmp, selection.method = "vst", verbose = F)
Seurat_tmp <- ScaleData(Seurat_tmp, verbose = F)
Seurat_tmp <- RunPCA(Seurat_tmp, features = VariableFeatures(Seurat_tmp), verbose = F)
Seurat_tmp <- FindNeighbors(Seurat_tmp, dims = 1:10, verbose = F)
Seurat_tmp <- FindClusters( object = Seurat_tmp, resolution = 0.6)
Seurat_tmp <- RunTSNE(object = Seurat_tmp, dims = 1:10)
Seurat_tmp <- RunUMAP(object = Seurat_tmp, dims = 1:10)
Seurat_tmp@graphs$RNA_snn<-as.Graph(sim_graph)
#SetAssayData(Seurat_tmp,layer = "data",new.data = as.matrix(t(tot_ST)), assay = "RNA")
#Seurat_tmp$RNA@data<-as.matrix(t(tot_ST))
network  <- as.matrix(Seurat_tmp@graphs$RNA_snn)


all(colnames(GSE1456) == rownames(meta_data))
phenotype <- cbind(meta_data[,5],meta_data[,3])
colnames(phenotype) <- c("time", "status")
test1<-as.numeric(phenotype[,1])
test1<-as.matrix(test1)
test2<-as.numeric(phenotype[,2])
test2<-as.matrix(test2)
phenotype<-cbind(test1,test2)
colnames(phenotype) <- c("time", "status")
head(phenotype)


cutoff = 0.2
alpha = 0.03
#alpha = seq(1,10,2)/1000
Save_file = c("/Users/naminiyakan/Documents/BulkToST/Res-BRCA-COMMOT/Scissor_BRCA-COMMOT_DEATH_SURVDEATHPheno.RData")
family = c("cox")
common <- intersect(rownames(GSE1456), rownames(Seurat_tmp))

if (class(Seurat_tmp) == "Seurat"){
  sc_exprs <- as.matrix(Seurat_tmp@assays$RNA$data)
  network  <- as.matrix(Seurat_tmp@graphs$RNA_snn)
}

#scanpy_norm<-read.csv(file='/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-COMMOT/processed/Normalized_5k_BRCA_COMMOT.csv',header = T,row.names = 1)
#scanpy_norm<-as.matrix(t(scanpy_norm))
#sc_exprs<-scanpy_norm
#common <- intersect(rownames(GSE1456), rownames(sc_exprs))

library(preprocessCore)
dataset0 <- cbind(GSE1456[common,], sc_exprs[common,])         # Dataset before quantile normalization.
dataset1 <- normalize.quantiles(as.matrix(dataset0))                           # Dataset after  quantile normalization.
rownames(dataset1) <- rownames(dataset0)
colnames(dataset1) <- colnames(dataset0)

Expression_bulk <- dataset1[,1:ncol(GSE1456)]
Expression_cell <- dataset1[,(ncol(GSE1456) + 1):ncol(dataset1)]
X <- cor(Expression_bulk, Expression_cell)


quality_check <- quantile(X)
print("|**************************************************|")
print("Performing quality-check for the correlations")
print("The five-number summary of correlations:")
print(quality_check)
print("|**************************************************|")
if (quality_check[3] < 0.01){
  warning("The median correlation between the single-cell and bulk samples is relatively low.")
}


if (family == "cox"){
  Y <- as.matrix(phenotype)
  if (ncol(Y) != 2){
    stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
  }else{
    print("Perform cox regression on the given clinical outcomes:")
  }
}

save(X, Y, network, Expression_bulk, Expression_cell, file = Save_file)

for (i in 1:length(alpha)){
  set.seed(123)
  fit0 <- APML1(X, Y, family = family, penalty = "Net", alpha = alpha[i], Omega = network, nlambda = 100, nfolds = min(10,nrow(X)))
  fit1 <- APML1(X, Y, family = family, penalty = "Net", alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)
  if (family == "binomial"){
    Coefs <- as.numeric(fit1$Beta[2:(ncol(X)+1)])
  }else{
    Coefs <- as.numeric(fit1$Beta)
  }
  Cell1 <- colnames(X)[which(Coefs > 0)]
  Cell2 <- colnames(X)[which(Coefs < 0)]
  percentage <- (length(Cell1) + length(Cell2)) / ncol(X)
  print(i)
  print(sprintf("alpha = %s", alpha[i]))
  print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", length(Cell1), length(Cell2)))
  print(sprintf("The percentage of selected cell is: %s%%", formatC(percentage*100, format = 'f', digits = 3)))
  
  if (percentage < cutoff){
    break
  }
  cat("\n")
}
print("|**************************************************|")

res<-list(para = list(alpha = alpha[i], lambda = fit0$lambda.min, family = family),
          Coefs = Coefs,
          Scissor_pos = Cell1,
          Scissor_neg = Cell2)

Scissor_select <- rep(0, ncol(sc_exprs))
names(Scissor_select) <- colnames(sc_exprs)
Scissor_select[res$Scissor_pos] <- 1
Scissor_select[res$Scissor_neg] <- 2

Scissor_select<-as.factor(as.matrix(Scissor_select))


dat <- data.frame("emb1" = pos[,1],
                  "emb2" = pos[,2],
                  "Cluster" = Scissor_select)


plt<-ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(alpha=1,ggplot2::aes(x = emb2, y = emb1,
                                           color = Cluster), size = 0.9) +
  
  ggplot2::scale_color_manual(values = rainbow(n = length(levels(Scissor_select)))) +
  
  # ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +
  
  ggplot2::labs(title = "",
                x = "X",
                y = "Y") +
  
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text( color = "black"),
                 axis.text.y = ggplot2::element_text( color = "black"),
                 axis.title.y = ggplot2::element_text(),
                 axis.title.x = ggplot2::element_text(),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text( colour = "black"),
                 legend.title = ggplot2::element_text( colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1.5, colour = "black")
                 # legend.position="none"
  ) +
  
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +
  
  ggplot2::coord_equal()

plt



