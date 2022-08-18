# Run RePACT on snATAC Seurat object
## Prepare the object for RePACT
This step is to make sure your Seurat obj contains LSI info
```
OBJ <- RunTFIDF(OBJ)
OBJ <- FindTopFeatures(OBJ, min.cutoff = "q5")
OBJ <- RunSVD(OBJ)
OBJ <- RunUMAP(OBJ, reduction = "lsi", dims = 1:50)
```
Load OBJ
```
OBJ <- readRDS("data/Beta.snATAC.rds")

```
