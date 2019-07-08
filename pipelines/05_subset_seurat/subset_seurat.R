# import libraries
library(Seurat)

seurat.obj.addr = "../../data/test_yolk_sac_subset.RDS"
save.at = "../../data/test_yolk_sac_subset_cluster_1.RDS"
process = T
add.dr  = T
filter.args = list("cell.labels" = c("HSC/MPP" , "Neutrophil-myeloid progenitor", "Monocyte-DC precursor", "DC2"), "gender" = c("male"))


#######################################################################################################
#######################################################################################################
#######################################################################################################

# load the seurat object
print("Loading the data ... ")
seurat.obj = readRDS(seurat.obj.addr)
print("Data loaded")

source("../../tools/bunddle_utils.R")

cells.to.keep = rep(T, length(seurat.obj@ident))
for(k in 1:length(filter.args)){
  cat = filter.args[k]
  conditions = as.vector(unlist(cat))
  cat = as.vector(names(cat))
  satisfy = seurat.obj@meta.data[, cat] %in% conditions
  cells.to.keep = cells.to.keep & satisfy
}

if(!any(cells.to.keep)){
  print("No cells have been selected. Relax the conditions.")
}else{
  seurat.obj = SubsetData(object=seurat.obj, cells.use=names(seurat.obj@ident)[cells.to.keep], subset.raw=T, do.clean=T)
  # add processing
  
  if (process){
    # normaliza data
    print("Normalizing data ...")
    seurat.obj = NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
    
    print("Computing variable genes ...")
    # find variable genes
    seurat.obj = FindVariableGenes(object = seurat.obj, mean.function = ExpMean, 
                                   dispersion.function = LogVMR, x.low.cutoff = .0125, 
                                   x.high.cutoff = 3, y.cutoff = .625)
    # calculate percentage of variable genes
    print(paste("Percentage of variable genes:", round(100 * length(seurat.obj@var.genes) / dim(seurat.obj@data)[1], digits = 2), sep = " "))
    
    # scale data in variable genes, otherwise pca is not possible
    print("Scaling data ...")
    seurat.obj = ScaleData(object=seurat.obj)
    
    # run PCA
    print("Performing PCA ...")
    seurat.obj = RunPCA(object = seurat.obj, pc.genes = seurat.obj@var.genes, do.print = TRUE, pcs.print = 1:20, genes.print = 10)
  }
  
  if(add.dr){
    # run TSNE
    print("Performing TSNE")
    seurat.obj = RunTSNE(object=seurat.obj, dims.use=1:20, seed.use=42, do.fast=TRUE)
    
    # run umap
    print("running UMAP")
    umap.coordinates = RunUMAP(pca.df=seurat.obj@dr$pca@cell.embeddings, tool_addr=tool_addr, python.addr=python.addr)
    rownames(umap.coordinates) = names(seurat.obj@ident)
    seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="umap", slot="cell.embeddings", new.data=as.matrix(umap.coordinates))
    seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="umap", slot="key", new.data="umap")
    
    # run force-directed graph
    print("Running force directed graph")
    seurat.obj = BuildSNN(object=seurat.obj, reduction.type="pca", dims.use=1:20, plot.SNN=F, force.recalc=TRUE, prune.SNN=.1)
    fdg_coordinates = runFDG(pca.df=seurat.obj@dr$pca@cell.embeddings, snn=seurat.obj@snn, iterations=2000, tool_addr=tool_addr, python.addr=python.addr)
    seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="fdg", slot="cell.embeddings", new.data=as.matrix(fdg_coordinates))
    seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="fdg", slot = "key", new.data = "fdg")
  }
  
  # save seurat object
  print(sprintf("saving data at: %s", save.at))
  saveRDS(seurat.obj, save.at)
}

file.remove("Rplots.pdf")

print("Ended beautifully ... ")
