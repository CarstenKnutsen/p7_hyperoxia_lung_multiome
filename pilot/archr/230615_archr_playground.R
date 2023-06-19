" Goal: archR on P7 lung multiome test ground
Date:230615
Author:Carsten Knutsen
"

output_fol <- '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/pilot/230615_archr_pilot'
dir.create(output_fol, recursive = TRUE,showWarnings = FALSE)
setwd(output_fol)

library(ArchR)
arch_proj <- loadArchRProject(path = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/pilot/230615_archr_pilot')

metadata = read.csv('/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files/share/p7_multiome_metadata.csv')
metadata$X
metadata$cellname <-paste('P7_lung_multiome#', metadata$X,sep = "")
metadata$cellname
subset_cells <- Reduce(intersect,list(metadata$cellname,arch_proj$cellNames))
subset_fol <- '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/pilot/230615_archr_trimmed_'
dir.create(subset_fol, recursive = TRUE,showWarnings = FALSE)

subset_proj <- subsetArchRProject(ArchRProj = arch_proj, outputDirectory = subset_fol,cells = subset_cells, force=true())

#metadata$cellname
#metadata[cellname]
#
#head(projHeme1$cellNames)
#saveArchRProject(ArchRProj = projHeme1, outputDirectory = output_fol, load = FALSE)