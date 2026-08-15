library(pheatmap)    # for heatmap visualization
library(RColorBrewer) # for color schemes
library(grid)
##############################*****KNOCKDOWN##################################
#raw_bulk<-read.table(file="/Users/naminiyakan/Documents/BulkToST/Dataset/Bulk-PDAC-Wu/GSE171485_raw_counts_GRCh38.p13_NCBI.tsv",sep="\t",header = T)
raw_bulk<-read.table(file="/Users/naminiyakan/Documents/BulkToST/Dataset/Bulk-PDAC-Wu/Knockdown-Bulk/GSE171486_norm_counts_TPM_GRCh38.p13_NCBI.tsv",sep="\t",header = T)
raw_bulk<-raw_bulk[-2152,]
raw_bulk<-raw_bulk[-2032,]
raw_bulk<-raw_bulk[-2031,]
rownames(raw_bulk)<-raw_bulk$Symbol
raw_bulk<-raw_bulk[,-1]
#raw_bulk<-raw_bulk[,16:21]
control<-raw_bulk[,16:18]
knockdown_DX<-raw_bulk[,19:21]
knockdown_MY<-raw_bulk[,25:27]
raw_bulk<-cbind(control,knockdown_DX,knockdown_MY)

bulk_meta<-read.table(file="/Users/naminiyakan/Documents/BulkToST/Dataset/Bulk-PDAC-Wu/Knockdown-Bulk/Metadata_Panc1-DDXandMYE.txt",header = F,row.names = 1,sep = '\t')



scanpy_norm<-read.csv(file='/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/Processed-A/Normalized_PDAC.csv',header = T,row.names = 1)
scanpy_norm<-as.matrix(t(scanpy_norm))
sc_exprs<-scanpy_norm
common <- intersect(rownames(raw_bulk), rownames(sc_exprs))

raw_bulk<-raw_bulk[common,]


# Define colors for the metadata (sample conditions)
annotation_colors <- list(Condition = c(Control = "gray", DDX60L_Knockdown = "purple", MYEOV_Knockdown = "gold"))


# Create a row annotation dataframe for the heatmap
sample_annotation <- data.frame(Condition = bulk_meta$V2)
rownames(sample_annotation) <- rownames(bulk_meta)


z_score_loc<-matrix(0,nrow(raw_bulk),ncol(raw_bulk))
for(i in 1:nrow(raw_bulk)){
  for(j in 1:ncol(raw_bulk)){
    z_score_loc[i,j]<-(as.numeric(raw_bulk[i,j])-mean(as.numeric(raw_bulk[i,])))/sd(as.numeric(raw_bulk[i,]))
  }
}
colnames(z_score_loc)<-colnames(raw_bulk)
rownames(z_score_loc)<-rownames(raw_bulk)

z_score_loc<-na.omit(z_score_loc)

# Plot the heatmap
pheatmap(z_score_loc, 
         annotation_col = sample_annotation, 
         annotation_colors = annotation_colors, 
         color = colorRampPalette(c("navy","white","firebrick3"))(100), # Color scheme
         cluster_rows = TRUE, cluster_cols = TRUE, # Clustering
         show_rownames = FALSE, show_colnames = FALSE, # Show labels
         main = "GSE171486", treeheight_row = 0,     # Hide row dendrogram
         treeheight_col = 0, border_color = "black")     # Hide column dendrogram) # Title

