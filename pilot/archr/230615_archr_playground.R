" Goal: archR on P7 lung multiome test ground
Date:230615
Author:Carsten Knutsen
"


output_fol <- '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/pilot/230615_archr_pilot'
dir.create(output_fol, recursive = TRUE,showWarnings = FALSE)
setwd(output_fol)

library(ArchR)
library(Seurat)
library(anndata)
arch_proj <- loadArchRProject(path = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/pilot/230615_archr_pilot')

metadata = read.table('/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files/share/p7_multiome_metadata.csv',
                      fill=TRUE,
                      colClasses = 'character',
                      sep = ",",
                      header = TRUE,
                      row.names = 1)
metadata$cellname <-paste('P7_lung_multiome#', row.names(metadata),sep = "")
subset_cells <- Reduce(intersect,list(metadata$cellname,arch_proj$cellNames))
metadata <-metadata[metadata$cellname %in% subset_cells, ]
metadata$celltype[metadata$celltype == 'Ecm1+ epi cell'] <- 'Basal'
as.data.frame(sapply(metadata,class))


subset_fol <- '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/pilot/230615_archr_trimmed_'
dir.create(subset_fol, recursive = TRUE,showWarnings = FALSE)
subset_proj <- subsetArchRProject(ArchRProj = arch_proj, outputDirectory = subset_fol,cells = subset_cells, force=true())

rna_fn <-"/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files/share/p7_multiome_rna.gz.h5ad"
data <- read_h5ad(rna_fn)
data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
data <- data[rownames(data) %in% rownames(metadata),]
data <- RenameCells(data, new.names =paste("P7_lung_multiome#",Cells(data), sep=""))
scRNA <- as.SingleCellExperiment(data)
addGeneExpressionMatrix(input = subset_proj, seRNA = scRNA, force = TRUE)

paste("P7_lung_multiome#",Cells(object = data))
Cells(data)
for (col in colnames(metadata)){
  subset_proj <- addCellColData(ArchRProj = subset_proj, data = metadata[,col], cells = subset_cells, name = col  )
}

subset_proj <- addIterativeLSI(
  ArchRProj = subset_proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 

)

subset_proj <- addHarmony(
  ArchRProj = subset_proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "mouse"
)

subset_proj <- addClusters(
  input = subset_proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "cluster",
  resolution = 0.8
)

subset_proj <- addUMAP(
  ArchRProj = subset_proj, 
  reducedDims = "Harmony", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
p1 <- plotEmbedding(ArchRProj = subset_proj, colorBy = "cellColData", name = "mouse", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = subset_proj, colorBy = "cellColData", name = "cluster", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = subset_proj, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = subset_proj, colorBy = "cellColData", name = "treatment", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = subset_proj, colorBy = "cellColData", name = "lineage", embedding = "UMAP")
plotPDF(list(p1, p2,p3,p4,p5), name = "archr_atac_umaps.pdf", ArchRProj = subset_proj, addDOC = FALSE, width = 5, height = 5)
help(plotPDF)

subset_proj <-addGroupCoverages(
  ArchRProj = subset_proj,
  groupBy = 'celltype',
)


pathToMacs2 <- findMacs2()
subset_proj <- addReproduciblePeakSet(
  ArchRProj = subset_proj, 
  groupBy = "celltype", 
  pathToMacs2 = pathToMacs2
)
subset_proj <- addPeakMatrix(subset_proj)
markersPeaks <- getMarkerFeatures(
  ArchRProj = subset_proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.5",
  limits = c(-2, 2),
  nLabel = 1,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  )
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 10, height = 10, ArchRProj = subset_proj, addDOC = FALSE)

subset_proj$celltype_treatment <- paste(subset_proj$celltype,"_" ,subset_proj$treatment)
ArchRBrowser(subset_proj)
p <- plotBrowserTrack(
  ArchRProj = subset_proj, 
  groupBy = "celltype", 
  useGroups = c('Arterial EC','Cap1','Cap2','Lymphatic EC','Proliferating EC','Venous EC'),
  geneSymbol = c("Slc6a2"),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Venous EC"],
  upstream = 50000,
  downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$Slc6a2)
plotPDF(p, name = "Plot-Tracks-With-Features-slc6a2", width = 5, height = 5, ArchRProj = subset_proj, addDOC = FALSE)

p <- plotBrowserTrack(
  ArchRProj = subset_proj, 
  groupBy = "celltype", 
  geneSymbol = c("Acta1"),
  upstream = 50000,
  downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$Slc6a2)
plotPDF(p, name = "Plot-Tracks-With-Features-acta1", width = 5, height = 10, ArchRProj = subset_proj, addDOC = FALSE)


