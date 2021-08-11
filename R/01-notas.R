# Se requiere la instalación de últiples paquetes
# BiocManager::install(c("batchelor", "BiocFileCache", "BiocSingular", "bluster",
#"cowplot", "dplyr", "DropletUtils", "EnsDb.Hsapiens.v86", "ExperimentHub", "fossil",
#"gert", "gh", "here", "iSEE", "kableExtra", "org.Mm.eg.db", "patchwork", "PCAtools", "pheatmap",
#"plotly", "pryr", "RColorBrewer", "rsthemes", "Rtsne", "scater", "scPipe", "scran", "scRNAseq", "sessioninfo",
#"Seurat", "SingleCellExperiment", "suncalc", "TENxPBMCData", "usethis", "uwot"))


#Tengo problemas con la instalación de los siguientes
##BiocManager::install(c("batchelor", "DropletUtils", "EnsDb.Hsapiens.v86", "iSEE", "scater", "scPipe", "scran", "scRNAseq"))

#El problema provenía de mi versión de R, al haber instalado la versión para M1 nativa, no permitía instalar
#paquetes que no estuvieran diseñados para M1, para lograrlo, tuve que reinstalar R,
#Usando la versión diseñada para la arquitectura intel, de manera que utilice Rosetta para traducir
#entre arquitecturas y así sea posible instalar librerías diseñadas para intel.


#Uso del booleano %in% int_gen
