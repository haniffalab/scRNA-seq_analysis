args = commandArgs(trailingOnly=T)
seurat.addr = args[1]
set.ident = "cell.labels"

passed.threshold         = .5
competition.threshold    = .8
expansion                = 4
    
library(Seurat)
library(plyr)
library(dplyr)

source("../../tools/bunddle_utils.R")

seurat.addr = file.path("../../data", seurat.addr)

gene_mean_expression = function(seurat.obj){
  no.genes = nrow(seurat.obj@data)
  start_index = 1
  while (start_index < no.genes){
    end_index = start_index + 999
    end_index = min(end_index, no.genes)
    expression.data_ = data.matrix(seurat.obj@data[start_index:end_index, ])
    expression.data_ = t(expression.data_)
    expression.data_ = as.data.frame(expression.data_)
    expression.data_ = cbind(data.frame(CellLabels = as.vector(seurat.obj@ident)), expression.data_)
    expression.data_ = aggregate(expression.data_[2:dim(expression.data_)[2]], list(expression.data_$CellLabels), mean)
    expression.data_ = cbind(data.frame(CellType = expression.data_$Group.1), expression.data_[, 2:dim(expression.data_)[2]])
    rownames(expression.data_) = expression.data_$CellType
    expression.data_ = expression.data_[, 2:ncol(expression.data_)]
    print(start_index)
    if (start_index == 1){
      expression.data = expression.data_
    }else{
      expression.data = cbind(expression.data, expression.data_)
    }
    start_index = start_index + 1000
  }
  expression.data
}

print("Loading data ... ")
seurat.obj = readRDS(seurat.addr)
seurat.obj = SetAllIdent(seurat.obj, id=set.ident)
print("Data loaded.")

# get aggregated expression matrix
expression.data = gene_mean_expression(seurat.obj)
saveRDS(expression.data, "gene_mean_expression.RDS")

#
gene.names = colnames(expression.data)
cell.types = rownames(expression.data)
n.cell.types = length(cell.types)
super.markers = rep("", length(cell.types))
names(super.markers) = cell.types

gene.max.expression  = apply(expression.data, 2, max)
expression.data.norm = t(t(expression.data) / gene.max.expression)

for(j in seq_along(gene.names)){
  gene.name = gene.names[j]
  super.vector = expression.data.norm[, gene.name]
  super.order  = order(super.vector, decreasing=T)
  cell.type    = names(super.vector[super.order[1]])
  expr.val     = as.vector(expression.data[cell.type, gene.name] )
  competition = all(super.vector[super.order[expansion:length(super.order)]] < competition.threshold)
  if(competition & expr.val > passed.threshold){
    super.markers[cell.type] = paste(super.markers[cell.type], gene.name, sep = ", ")
  }
}

# order markers
for (i in seq_along(super.markers)){
  ms        = super.markers[i]
  cell.type = names(ms)[1]
  ms = unlist(strsplit(as.character(ms), split=", "))
  ms = ms[ms != ""]
  expression.row = expression.data[cell.type, ms]
  ms = ms[rev(order(expression.row))]
  ms = paste(ms, collapse = ", ")
  super.markers[i] = ms
}

# make classification markers
classification.markers = super.markers

seurat.obj = FindVariableGenes(object = seurat.obj, mean.function = ExpMean, 
                               dispersion.function = LogVMR, x.low.cutoff = .0125, 
                               x.high.cutoff = 9, y.cutoff = .2)

for (k in seq_along(cell.types)){
  cell.type = cell.types[k]
  print("printing cell.type")
  print(cell.type)
  eval(parse(text = sprintf("cells.markers = super.markers[['%s']]", cell.type)))
  cells.markers = unlist(strsplit(cells.markers, split=", "))
  print("printing cell.markers")
  print(cells.markers)
  cor.expression.data = expression.data[, cells.markers]
  print("printing cor.expression.data")
  print(cor.expression.data)
  compare.to = cell.types[cell.types != cell.type]
  print("printing compare.to")
  print(compare.to)
  t.cor.expression.data = t(cor.expression.data)
  colnames(t.cor.expression.data) <- cell.types
  rownames(t.cor.expression.data) <- cells.markers
  print("printing t.cor.expression.data")
  print(t.cor.expression.data)
  cor.dist = cor(t.cor.expression.data, method="spearman")
  rownames(cor.dist) <- cell.types
  colnames(cor.dist) <- cell.types
  print("printing cor.dist")
  print(cor.dist)
  cor.dist = cor.dist[cell.type, compare.to]
  print("printing cor.dist[cell.type, compare.to]")
  print(cor.dist)
  compare.to = names(cor.dist)[cor.dist > .3]
  print("printing length of compare.to")
  print(length(compare.to))
  if(length(compare.to) > .45){
    cor.dist = cor.dist[compare.to]
    cor.dist = cor.dist[order(cor.dist, decreasing=T)]
    compare.to = names(cor.dist)[1:3]
  }
  if (length(compare.to) > 0 & !any(is.na(compare.to))){
    print(sprintf("Comparing %s to: %s", cell.type, paste(compare.to, collapse = ", ")))
    DEGs = list()
    for (m in seq_along(compare.to)){
      DEG = FindMarkers(object=seurat.obj, ident.2=compare.to[m], ident.1=cell.type,max.cells.per.ident=500, only.pos=T, min.pct=.5, genes.use=seurat.obj@var.genes, )
      DEG$cluster = compare.to[m]
      DEG$gene = rownames(DEG)
      DEGs[[m]] = DEG
    }
    DEGs = Reduce(f=rbind, x=DEGs)
    additional.markers = (DEGs %>% group_by(cluster) %>% top_n(5, avg_logFC))$gene
    new.markers = unique(c(cells.markers, additional.markers))
    new.markers = paste(new.markers, collapse = ", ")
    eval(parse(text = sprintf("classification.markers[['%s']] = new.markers", cell.type)))
  }
}


# order markers
for (i in seq_along(classification.markers)){
  ms        = classification.markers[i]
  cell.type = names(ms)[1]
  ms = unlist(strsplit(as.character(ms), split=", "))
  ms = ms[ms != ""]
  expression.row = expression.data[cell.type, ms]
  ms = ms[rev(order(expression.row))]
  ms = paste(ms, collapse = ", ")
  classification.markers[i] = ms
}

saveRDS(classification.markers, "signatures.RDS")

sig.df = data.frame(CellTypes = names(classification.markers), Signatures = classification.markers)
write.csv(sig.df, "Signatures.csv", row.names = F)

print("Ended beautifully ... ")
