library(Seurat)
library(SoupX)
library(DropletUtils)
library(ggplot2)
library(knitr)

### MCAI016
run_soupx <- function(toc,tod,rho=NULL) {
  all <- toc
  all <- CreateSeuratObject(all)
  all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
  all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(all)
  all <- ScaleData(all, features = all.genes)
  all <- RunPCA(all, features = VariableFeatures(all), npcs = 30, verbose = F)
  all <- FindNeighbors(all, dims = 1:30)
  all <- FindClusters(all, resolution = 0.5)
  all <- RunUMAP(all, dims = 1:30)
  matx <- all@meta.data

  sc = SoupChannel(tod, toc)
  sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))

  if (is.null(rho)) {
   tryCatch(
    {sc = autoEstCont(sc)},
     error=function(e) {
    print("autoEstCont Error !")
    sc = setContaminationFraction(sc, 0.2)}
   )
   }else{
   sc = setContaminationFraction(sc, rho)
 }
  out = adjustCounts(sc, roundToInt = TRUE)
  saveRDS(sc,"Young/MCAI016/sc.rds")
  DropletUtils:::write10xCounts("Young/MCAI016/soupX_matrix", out,version="3")
}


data_dir<-"Young/MCAI016/RawMatrix"
list.files(data_dir)
tod <- Read10X(data.dir = data_dir, gene.column = 1)

data_dir<-"Young/MCAI016/04.Matrix/FilterMatrix"
list.files(data_dir)
toc <- Read10X(data.dir = data_dir, gene.column = 1)

tod <- tod[rownames(toc),]
run_soupx(toc,tod, rho = 0.2)