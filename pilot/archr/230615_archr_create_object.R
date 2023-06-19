" Goal: Create archR object on P7 lung multiome 
Date:230615
Author:Carsten Knutsen
"
#BiocManager::install("TFBSTools")
#BiocManager::install("BSgenome")
#devtools::install_github("GreenleafLab/motifmatchr")
#BiocManager::install(c("Bioconductor/GenomeInfoDb", "Bioconductor/BSgenome")) # recieved error this solved https://github.com/Bioconductor/BSgenome/issues/9


###Installation was a big ole headache
#Sys.setenv(CONDA_BUILD_SYSROOT="/")
#devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

output_fol <- '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/pilot/230615_archr_pilot'
dir.create(output_fol, recursive = TRUE,showWarnings = FALSE)
setwd(output_fol)

library(ArchR)
#ArchR::installExtraPackages()
addArchRThreads(threads = 6)
addArchRGenome("mm10")

inputFiles <- c('/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/cellranger_output/230609_aggregate/outs/atac_fragments.tsv.gz')
type(inputFiles)
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = 'P7_lung_multiome',
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = output_fol,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
saveArchRProject(ArchRProj = projHeme1, outputDirectory = output_fol, load = FALSE)