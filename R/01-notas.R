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


library("gitcreds", "gert")
library('gh')

#Revisar la configuración de Git
gh::gh_whoami() # para checar cómo quedó la configuración

gert::git_commit_all("Nuevos")
gert::git_push()

usethis::use_r('08_reduccion_de_dimensiones.Rmd')

##Change the error
#options(error= function() starwarssay::say("Do or do not... There is no try", by = "baby_yoda"))

options(prompt = "🦖 ~ 🧬 >")

starwarssay::say("Join the dark side.... #CODE \n Welcome back Pablo", by = "darth_vader")

usethis::edit_r_profile()
starwarssay::say("Be a Jedi", by = "yoda")




if (interactive()) {
  .First <- function() {
    # aquí van los mensajes o las funciones para la bienvenida
  }

  .Last <- function() {
    # aquí van los mensajes o las funciones para la despedida
  }
}
