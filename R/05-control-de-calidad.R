## ----cargar_paquetes, message = FALSE---------------------------------------------------------------------------------
## Paquetes de este capítulo
library("scRNAseq") ## para descargar datos de ejemplo
library("AnnotationHub") ## para obtener información de genes
library("scater") ## para gráficas y control de calidad
library("BiocFileCache") ## para descargar datos
library("DropletUtils") ## para detectar droplets
library("Matrix") ## para leer datos en formatos comprimidos


## ----datos_ejercicio_1------------------------------------------------------------------------------------------------
## Datos
library("scRNAseq")
sce.416b <- LunSpikeInData(which = "416b")
sce.416b$block <- factor(sce.416b$block)

# Descarga los archivos de anotación de la base de datos de Ensembl
# correspondientes usando los recursos disponibles vía AnnotationHub
library("AnnotationHub")
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))
# Obtén la posición del cromosoma para cada gen
ens.mm.v97 <- ah[["AH73905"]]
location <- mapIds(
  ens.mm.v97,
  keys = rownames(sce.416b),
  keytype = "GENEID",
  column = "SEQNAME"
)

# Identifica los genes mitocondriales
is.mito <- which(location == "MT")
library("scater")
sce.416b <- addPerCellQC(sce.416b,
                         subsets = list(Mito = is.mito)
)
## Si quieres guarda los resultados de addPerCellQC() para responder
## las preguntas del ejercicio. Eventualmente si necesitaremos los
## resultados de addPerCellQC() para las secciones posteriores a este
## ejercicio.


## ----visualizar_qc----------------------------------------------------------------------------------------------------
plotColData(sce.416b, x = "block", y = "detected")
plotColData(sce.416b, x = "block", y = "detected") +
  scale_y_log10()

#Ejercicio
#En colData
boxplot(detected ~ block, data=colData(sce.416b))

###PROBLEMA
plotColData(sce.416b,
            x = "block",
            y = "detected",
            other_fields = "phenotype"
) +
  scale_y_log10() +
  facet_wrap(~phenotype)
###BETTER

plotColData(sce.416b,
            x = "phenotype",
            y = "detected",
            other_fields = "block"
) +
  scale_y_log10() +
  facet_wrap(~block)

## ----qc_altexps_ERCC_percent, echo = FALSE----------------------------------------------------------------------------
plotColData(sce.416b, x = "block", y = "altexps_ERCC_sum")

plotColData(sce.416b,
            x = "block",
            y = "altexps_ERCC_percent",
            other_fields = "phenotype", colour_by = 'phenotype'
) + facet_wrap(~phenotype)













## ----valores_qc-------------------------------------------------------------------------------------------------------
# Valores de límite ejemplo
qc.lib <- sce.416b$sum < 100000
qc.nexprs <- sce.416b$detected < 5000
qc.spike <- sce.416b$altexps_ERCC_percent > 10
qc.mito <- sce.416b$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito

# Obtenemos un resumen del número de células
# eliminadas por cada filtro
DataFrame(
  LibSize = sum(qc.lib),
  NExprs = sum(qc.nexprs),
  SpikeProp = sum(qc.spike),
  MitoProp = sum(qc.mito),
  Total = sum(discard)
)
#-----------------------------------VALORES DE CORTE MÓVILES--------------------------
## Usando isOutlier() para determinar los valores de corte
qc.lib2 <- isOutlier(sce.416b$sum, log = TRUE, type = "lower")
qc.nexprs2 <- isOutlier(sce.416b$detected,
                        log = TRUE,
                        type = "lower"
)
qc.spike2 <- isOutlier(sce.416b$altexps_ERCC_percent,
                       type = "higher"
)
qc.mito2 <- isOutlier(sce.416b$subsets_Mito_percent,
                      type = "higher"
)
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2

# Extraemos los límites de valores (thresholds)
attr(qc.lib2, "thresholds")
attr(qc.nexprs2, "thresholds")

# Obtenemos un resumen del número de células
# eliminadas por cada filtro
DataFrame(
  LibSize = sum(qc.lib2),
  NExprs = sum(qc.nexprs2),
  SpikeProp = sum(qc.spike2),
  MitoProp = sum(qc.mito2),
  Total = sum(discard2)
)

## Más pruebas
plotColData(sce.416b,
            x = "block",
            y = "detected",
            other_fields = "phenotype"
) +
  scale_y_log10() +
  facet_wrap(~phenotype)

## Determino el bloque (batch) de muestras
batch <- paste0(sce.416b$phenotype, "-", sce.416b$block)

## Versión de isOutlier() que toma en cuenta los bloques de muestras
qc.lib3 <- isOutlier(sce.416b$sum,
                     log = TRUE,
                     type = "lower",
                     batch = batch
)
qc.nexprs3 <- isOutlier(sce.416b$detected,
                        log = TRUE,
                        type = "lower",
                        batch = batch
)
qc.spike3 <- isOutlier(sce.416b$altexps_ERCC_percent,
                       type = "higher",
                       batch = batch
)
qc.mito3 <- isOutlier(sce.416b$subsets_Mito_percent,
                      type = "higher",
                      batch = batch
)
discard3 <- qc.lib3 | qc.nexprs3 | qc.spike3 | qc.mito3

# Extraemos los límites de valores (thresholds)
attr(qc.lib3, "thresholds")
attr(qc.nexprs3, "thresholds")

# Obtenemos un resumen del número de células
# eliminadas por cada filtro
DataFrame(
  LibSize = sum(qc.lib3),
  NExprs = sum(qc.nexprs3),
  SpikeProp = sum(qc.spike3),
  MitoProp = sum(qc.mito3),
  Total = sum(discard3)
)


###FUE NECESARIO UTILIZAR QC.LIB
table(qc.lib, qc.nexprs | qc.spike | qc.mito)

#SÍ, PORQUE OBTENEMOS DOS


#¿CUÁL ES MEJOR, DISCARD O DISCARD2?
addmargins(table(discard,discard2))

#USAR BATCH MEJORA LA CALIDAD
table(discard,discard2,discard3) #recuerda que TRUE es baja calidad


## ----grun_problema----------------------------------------------------------------------------------------------------
sce.grun <- GrunPancreasData()
sce.grun <- addPerCellQC(sce.grun)

## ¿Qué patrón revela esta gráfica?
plotColData(sce.grun, x = "donor", y = "altexps_ERCC_percent")



###Hay un problema grave en los datos de D10, D3

## ----grun_isOutlier---------------------------------------------------------------------------------------------------
## isOutlier() puede ayudarnos cuando un grupo de muestras
## tuvo más problemas que el resto
discard.ercc <- isOutlier(sce.grun$altexps_ERCC_percent,
                          type = "higher",
                          batch = sce.grun$donor
)
discard.ercc2 <- isOutlier(
  sce.grun$altexps_ERCC_percent,
  type = "higher",
  batch = sce.grun$donor,
  subset = sce.grun$donor %in% c("D17", "D2", "D7")
)# Solo las células de los donadores chidos toma en cuenta para el valor de corte, y lo usa para los demás

## isOutlier() tomando en cuenta el batch
plotColData(
  sce.grun,
  x = "donor",
  y = "altexps_ERCC_percent",
  colour_by = data.frame(discard = discard.ercc)
) #Poner colores en las que eliminariamos


## isOutlier() tomando en cuenta batch y muestras que fallaron
plotColData(
  sce.grun,
  x = "donor",
  y = "altexps_ERCC_percent",
  colour_by = data.frame(discard = discard.ercc2)
)




## ----qc_extra_416b----------------------------------------------------------------------------------------------------
# Agregamos información sobre que células
# tienen valores extremos
sce.416b$discard <- discard2

# Haz esta gráfica para cada medida de
# control de calidad (QC)
plotColData(
  sce.416b,
  x = "block",
  y = "sum",
  colour_by = "discard",
  other_fields = "phenotype"
) +
  facet_wrap(~phenotype) +
  scale_y_log10()

# Otra gráfica de diagnóstico útil
plotColData(
  sce.416b,
  x = "sum",
  y = "subsets_Mito_percent",
  colour_by = "discard",
  other_fields = c("block", "phenotype")
) +
  facet_grid(block ~ phenotype)
#plotly::ggplotly(p)
#Gráfica interactiva cuando guardamos esto

# -------------------------------EJERCICIO

sce.grun$discard <- discard.ercc2


plotColData(
  sce.grun,
  x = "altexps_ERCC_detected",
  y = "altexps_ERCC_sum",
  colour_by = "discard",
  other_fields = c("donor")
) +
  facet_grid(~donor)



#
#
#
#
#













## ----qc_extra_grun, echo = FALSE--------------------------------------------------------------------------------------
sce.grun$discard <- discard.ercc2

plotColData(
  sce.grun,
  x = "altexps_ERCC_detected",
  y = "altexps_ERCC_sum",
  colour_by = "discard",
  other_fields = c("donor")
) +
  facet_grid(~donor)





######IDENTIFICACIÓN DE DROPLETS VACIOS

###SISTEMA DE CHROMIUM
## Descarguemos los datos
library("BiocFileCache")
bfc <- BiocFileCache()
raw.path <-
  bfcrpath(
    bfc,
    file.path(
      "http://cf.10xgenomics.com/samples",
      "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
    )
  )
untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

## Leamos los datos en R
library("DropletUtils")
library("Matrix")
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = TRUE)
bcrank <- barcodeRanks(counts(sce.pbmc))

# Mostremos solo los puntos únicos para acelerar
# el proceso de hacer esta gráfica
uniq <- !duplicated(bcrank$rank)
plot(
  bcrank$rank[uniq],
  bcrank$total[uniq],
  log = "xy",
  xlab = "Rank",
  ylab = "Total UMI count",
  cex.lab = 1.2
)
abline(
  h = metadata(bcrank)$inflection,
  col = "darkgreen",
  lty = 2
)
abline(
  h = metadata(bcrank)$knee,
  col = "dodgerblue",
  lty = 2
)
legend(
  "bottomleft",
  legend = c("Inflection", "Knee"),
  col = c("darkgreen", "dodgerblue"),
  lty = 2,
  cex = 1.2
)


## ----emptyDrops-------------------------------------------------------------------------------------------------------
## Usemos DropletUtils para encontrar los droplets
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))


# Revisa ?emptyDrops para una explicación de porque hay valores NA
summary(e.out$FDR <= 0.001)
set.seed(100)
limit <- 100
all.out <-
  emptyDrops(counts(sce.pbmc), lower = limit, test.ambient = TRUE) #TEST AMBIENT
#Significa que calculara los resultados para todos los barcodes

# Idealmente, este histograma debería verse uniforme.
# Picos grandes cerca de cero indican que los _barcodes_
# con un número total de cuentas menor a "lower" no son
# de origen ambiental.
hist(all.out$PValue[all.out$Total <= limit &
                      all.out$Total > 0],
     xlab = "P-value",
     main = "",
     col = "grey80"
)
summary(all.out$PValue<= 0.001)
summary(e.out$PValue<=0.001)



###Porcentaje mitocondrial
## ----pbmc_chrMT_ayuda-------------------------------------------------------------------------------------------------
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]
is.mito <- grep("^MT-", rowData(sce.pbmc)$Symbol)
sce.pmbc <- addPerCellQC(sce.pbmc, subsets = list(MT = is.mito))
discard.mito <-
  isOutlier(sce.pmbc$subsets_MT_percent, type = "higher")
plot(
  sce.pmbc$sum,
  sce.pmbc$subsets_MT_percent,
  log = "x",
  xlab = "Total count",
  ylab = "Mitochondrial %"
)
abline(h = attr(discard.mito, "thresholds")["higher"], col = "red")


## ----pbmc_combined, echo = FALSE--------------------------------------------------------------------------------------
## Leer los datos crudos de nuevo
sce.pbmc <- read10xCounts(fname, col.names = TRUE)

## Almacenar si es o no es célula
sce.pbmc$is_cell <- e.out$FDR <= 0.001

## Filtar los casos para los cuales no tuvimos pruebas estadísticas
sce.pbmc <- sce.pbmc[, which(!is.na(sce.pbmc$is_cell))]

## Identicar chr mitocondrial
is.mito <- grep("^MT-", rowData(sce.pbmc)$Symbol)

## Calcular medidas de QC
sce.pbmc <- addPerCellQC(sce.pbmc, subsets = list(MT = is.mito))

## Almacenar cuales hay que descartar por MT
sce.pbmc$discard_mito <-
  isOutlier(sce.pbmc$subsets_MT_percent, type = "higher")

## Visualizar resultados
plotColData(
  sce.pbmc,
  x = "total",
  y = "subsets_MT_percent",
  colour_by = "discard_mito"
) + facet_grid(~ ifelse(sce.pbmc$is_cell, "Célula", "Droplet vacío"))


## ----filtrar_o_no-----------------------------------------------------------------------------------------------------
# Eliminemos las células de calidad baja
# al quedarnos con las columnas del objeto sce que NO
# queremos descartar (eso hace el !)
filtered <- sce.416b[, !discard2]
# Alternativamente, podemos marcar
# las células de baja calidad
marked <- sce.416b
marked$discard <- discard2



## ----isee_basic-------------------------------------------------------------------------------------------------------
## Hagamos un objeto sencillo de tipo RangedSummarizedExperiment
library("SummarizedExperiment")
## ?SummarizedExperiment

## De los ejemplos en la ayuda oficial
## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6

## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)

## Información de nuestros genes
rowRanges <- GRanges(
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))

## Información de nuestras muestras
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)

## Juntamos ahora toda la información en un solo objeto de R
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

## Exploremos el objeto resultante
rse

## Explora el objeto rse de forma interactiva
library("iSEE")
if (interactive()) {
  iSEE::iSEE(rse)
}


## ----isee_416b--------------------------------------------------------------------------------------------------------
## Explora el objeto sce.416b de forma interactiva
if (interactive()) {
  iSEE::iSEE(sce.416b, appTitle = "sce.416b")
}







## ---------------------------------------------------------------------------------------------------------------------
## Información de la sesión de R
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
