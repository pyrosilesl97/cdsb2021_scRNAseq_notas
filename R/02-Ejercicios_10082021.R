library("scRNAseq")
sce.416b <- LunSpikeInData(which = "416b")

# Carga el paquete SingleCellExperiment
library("SingleCellExperiment")

"Primera parte
aquí checamos el slot assays"

# Extrae la matriz de cuentas del set de datos de 416b
counts.416b <- counts(sce.416b)
# CHEQUEMOS clase y dimensiones
class(counts.416b) # es matriz
dim(counts.416b) # indicará genes y células


# CONSTRUIR un nuevo objeto SCE de la matriz de cuentas !!!!!!
sce <- SingleCellExperiment(assays = list(counts = counts.416b))

# Revisa el objeto que acabamos de crear
sce

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce)

# Accesa la matriz de cuenta del compartimento (slot) "assays"
# assays(sce, "counts")
# OJO: ¡esto puede inundar tu sesión de R!

# 1. El método general
assay(sce, "counts")[110:115, 1:3] # gene, cell

# 2. El método específico para accesar la matriz de cuentas "counts"
counts(sce)[110:115, 1:3]

# AGREGAR MAS ASSAYS
sce <- scater::logNormCounts(sce)
# Revisa el objeto que acabamos de actualizar
sce
## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce)

# 1. El método general
assay(sce, "logcounts")[110:115, 1:3]

# 2. El método específico para accesar la matriz de cuentas
# transformadas "logcounts"
logcounts(sce)[110:115, 1:3]

# agregemos un assay mas, esta vez de manera manual
assay(sce, "counts_100") <- assay(sce, "counts") + 100 # suma 100 a counts assay
# Enumera los "assays" en el objeto
assays(sce) # indica num y nombre de assays


assayNames(sce) # solo nos dará los nombres de los assays

#      assay(sce, "counts_100")[110:115, 1:3]
## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce)

############"segunda parte:aquí checaremos metadata de las células"
# Extrae la información de las muestras (metadata) del set de datos de 416b
colData.416b <- colData(sce.416b) # podemos checar objeto en la cajita de environment de RStudio!!

# explorar datooos
table(colData.416b$phenotype)


table(colData.416b$block) # fue en varios dias?


# Agrega algo de esa información a nuestro objeto de SCE
colData(sce) <- colData.416b[, c("phenotype", "block")]
# Revisa el objeto que acabamos de actualizar
sce

# Accesa a la información de las muestras (metadata) en nuestro SCE
colData(sce) # usar head?

# Accesa una columna específica de la información de las muestras (metadata)
table(sce$block)

# Ejemplo de una función que agrega columnas nuevas al colData
sce <- scater::addPerCellQC(sce.416b) # añade datos de control de calidad
# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
colData(sce)

# Revisa el objeto que acabamos de actualizar
sce

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce)

# Ejemplo: obtén el subconjunto de células de fenotipo "wild type"
# Acuérdate que las células son columnas del SCE !!!!
sce[, sce$phenotype == "wild type phenotype"]


###################################################
##"Tercera parte: examinaremos metadata de features (rowData)"
# Accesa la información de los genes de nuestro SCE
# ¡Está vació actualmente!
rowData(sce)

# Ejemplo de una función que agrega campos nuevos en el rowData
sce <- scater::addPerFeatureQC(sce)
# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
rowData(sce)

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce)

# Descarga los archivos de anotación de la base de datos de Ensembl
# correspondientes usando los recursos disponibles vía AnnotationHub
library("AnnotationHub")
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))



# Obtén la posición del cromosoma para cada gen
ensdb <- ah[["AH73905"]]
chromosome <- mapIds(ensdb,
                     keys = rownames(sce),
                     keytype = "GENEID",
                     column = "SEQNAME"
)
rowData(sce)$chromosome <- chromosome

# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
rowData(sce)

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce)
