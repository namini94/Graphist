#library(BayesSpace)
#melanoma <- getRDS("2018_thrane_melanoma", "ST_mel1_rep2")

mel_ST<-read.table(file="/Users/naminiyakan/Documents/BulkToST/Dataset/Melanoma-ST/ST-Melanoma-Datasets_1/ST_mel1_rep2_counts.tsv",sep='\t',header = T,row.names = 1)
colnames(mel_ST) <- gsub("^X", "", colnames(mel_ST))

pos <- do.call(rbind, strsplit(colnames(mel_ST), "x"))
pos <- as.data.frame(pos)

# Rename columns of the new data frame
colnames(pos) <- c("y", "x")

pos$y <- as.numeric(as.character(pos$y))
pos$x <- as.numeric(as.character(pos$x))


Genes <- do.call(rbind, strsplit(rownames(mel_ST), " "))
Genes <- as.data.frame(Genes)

# Rename columns of the new data frame
colnames(Genes) <- c("x", "ENSMBLID")
Genes <- as.data.frame(Genes$x)
colnames(Genes)<- c('x')

# Make row names unique by adding a suffix
rownames(mel_ST) <- make.unique(as.character(Genes$x))

colnames(mel_ST) <- paste0("S", seq_along(colnames(mel_ST)))
rownames(pos)<-colnames(mel_ST)

Barcodes<-as.data.frame(colnames(mel_ST))
colnames(Barcodes)<-c('x')


tissue_mask<-matrix(0,nrow(pos),1)
for(i in 1:nrow(pos)){
  tissue_mask[i,1]<-c("True")
}


write.csv(t(mel_ST),file="/Users/naminiyakan/Documents/BulkToST/Dataset/Melanoma-ST/processed/STcount.csv",quote = F,row.names = T)
write.csv(pos,file = "/Users/naminiyakan/Documents/BulkToST/Dataset/Melanoma-ST/processed/Pos.csv",quote = F,row.names = T)

write.csv(tissue_mask,file="/Users/naminiyakan/Documents/BulkToST/Dataset/Melanoma-ST/processed/tissue_mask.csv",quote = F,row.names = F)
write.csv(Genes,file="/Users/naminiyakan/Documents/BulkToST/Dataset/Melanoma-ST/processed/Gene_names.csv",quote = F,row.names = F)
write.csv(Barcodes,file="/Users/naminiyakan/Documents/BulkToST/Dataset/Melanoma-ST/processed/Barcodes.csv",quote = F,row.names = F)



#######################Bulk-BRCA-TCGA (Downloaded from Scissor paper)
GSE78220<-read.table(file="/Users/naminiyakan/Documents/BulkToST/Dataset/Melanoma-Bulk/GSE78220_norm_counts_TPM_GRCh38.p13_NCBI.tsv",sep='\t',header = T)
GSE78220 <- GSE78220[!GSE78220$Symbol %in% c("TRNAV-CAC"),]
rownames(GSE78220)<-GSE78220$Symbol
GSE78220<-GSE78220[,-1]

library(GEOquery)
gset <- getGEO("GSE78220", GSEMatrix =TRUE, getGPL=FALSE)

meta_data<-cbind(gset$GSE78220_series_matrix.txt.gz$source_name_ch1,gset$GSE78220_series_matrix.txt.gz$`overall survival (days):ch1`)
colnames(meta_data)<-c("Response","survival")

meta_res_binary<-read.csv(file = "/Users/naminiyakan/Documents/BulkToST/Dataset/Melanoma-Bulk/meta_response_binary.csv",row.names = 1,header = F)

meta_tot<-cbind(meta_data,meta_res_binary)
colnames(meta_tot)<-c("Response","survival","ResBin")

#######Network_Construction
neighbors<-read.csv(file = "/Users/naminiyakan/Documents/BulkToST/Dataset/Melanoma-ST/processed/neighbors.csv",header = T,row.names = 1)
#tot_ST<-read.csv(file = "/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/Total/tot_ST.csv",header = T,row.names = 1)
cd<-t(mel_ST)

sim_graph<-matrix(0,nrow(cd),nrow(cd))
colnames(sim_graph)<-rownames(cd)
rownames(sim_graph)<-rownames(cd)
for(i in 1:nrow(neighbors)){
  sim_graph[(neighbors$X0[i])+1,(neighbors$X1[i])+1]<-1
  sim_graph[(neighbors$X1[i])+1,(neighbors$X0[i])+1]<-1
}


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
network  <- as.matrix(Seurat_tmp@graphs$RNA_snn)


all(colnames(GSE78220) == rownames(meta_tot))

phenotype <-meta_tot$ResBin
tag <- c('Responders', 'Non-Responders')


cutoff = 0.2
alpha = 0.005
#alpha = seq(1,10,2)/1000
Save_file = c("/Users/naminiyakan/Documents/BulkToST/Res-Melanoma/Scissor_Melanoma_ResponsePheno.RData")
family = c("binomial")
common <- intersect(rownames(GSE78220), rownames(Seurat_tmp))

if (class(Seurat_tmp) == "Seurat"){
  sc_exprs <- as.matrix(Seurat_tmp@assays$RNA$data)
  network  <- as.matrix(Seurat_tmp@graphs$RNA_snn)
}

library(preprocessCore)
dataset0 <- cbind(GSE78220[common,], sc_exprs[common,])         # Dataset before quantile normalization.
dataset1 <- normalize.quantiles(as.matrix(dataset0))                           # Dataset after  quantile normalization.
rownames(dataset1) <- rownames(dataset0)
colnames(dataset1) <- colnames(dataset0)

Expression_bulk <- dataset1[,1:ncol(GSE78220)]
Expression_cell <- dataset1[,(ncol(GSE78220) + 1):ncol(dataset1)]
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

if (family == "binomial"){
  Y <- as.numeric(phenotype)
  z <- table(Y)
  if (length(z) != length(tag)){
    stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
  }else{
    print(sprintf("Current phenotype contains %d %s and %d %s samples.", z[1], tag[1], z[2], tag[2]))
    print("Perform logistic regression on the given phenotypes:")
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
#tot_pos<-read.csv(file="/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/Total/tot_pos.csv",row.names = 1)

dat <- data.frame("emb1" = pos$y,
                  "emb2" = pos$x,
                  "Cluster" = Scissor_select)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_tile(alpha=1,ggplot2::aes(x = emb1, y = emb2,
                                           fill = Cluster), color = "black",width = 1, height = 1) +
  
  ggplot2::scale_fill_manual(values = rainbow(n = length(levels(Scissor_select)))) +
  
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






