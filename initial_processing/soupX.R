" Goal: Cleanup gex data for etx scRNAseq 
Date:230129
Author:Carsten Knutsen
"
#install.packages('SoupX')
#install.packages('installr')
#library(installr)
#install.Rtools()
#install.packages('Matrix')
#BiocManager::install("DropletUtils")
library(DropletUtils)
library(Matrix)
library(SoupX)
library(Seurat)
data_dir <- '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/cellranger_output/230609_aggregate/outs/'
output_dir <-'/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files/soupx'

cellnames <- read.csv(sprintf('%s/filtered_feature_bc_matrix/barcodes.tsv.gz', data_dir),header =FALSE)
filt.matrix <- Read10X_h5(sprintf("%s/filtered_feature_bc_matrix.h5",data_dir),use.names = T)
raw.matrix <- Read10X_h5(sprintf("%s/raw_feature_bc_matrix.h5",data_dir),use.names = T)
filt.matrix.gene <- filt.matrix[["Gene Expression"]]
raw.matrix.gene <- raw.matrix[["Gene Expression"]]
soup.channel = SoupChannel(raw.matrix.gene, filt.matrix.gene)
srat <- CreateSeuratObject(counts = filt.matrix.gene)
srat <-RenameCells(srat, new.names = cellnames$V1)
srat    <- SCTransform(srat, verbose = F)
srat    <- RunPCA(srat, verbose = F)
srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat    <- FindClusters(srat, verbose = T)
meta    <- srat@meta.data
umap    <- srat@reductions$umap@cell.embeddings
soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)
soup.channel  <- autoEstCont(soup.channel)
adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
DropletUtils:::write10xCounts(output_dir, adj.matrix)

