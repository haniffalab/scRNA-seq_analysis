args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))
args = gsub(pattern = '@@', replacement = ' ', x = args)

arguments.list = "
seurat.addr.arg     = args[1]
set.ident.arg       = args[2]
cell.types.arg      = args[3]
root_cell_type.arg  = args[4]
type.to.colours.arg = args[5]
lineage.name.arg    = args[6]
"

expected_arguments = unlist(strsplit(arguments.list, "\n"))
expected_arguments = expected_arguments[!(expected_arguments == "")]

if(length(args) != length(expected_arguments)){
  error.msg = sprintf('This pipeline requires %s parameters', as.character(length(expected_arguments)))
  expected_arguments = paste(unlist(lapply(strsplit(expected_arguments, ".arg"), "[", 1)), collapse = "\n")
  stop(sprintf('This pipeline requires %s parameters: ', length(expected_arguments)))
}

eval(parse(text = arguments.list))

for(n in 1:length(expected_arguments)){
  argument = expected_arguments[n]
  argument.name = unlist(strsplit(argument, "="))[1]
  variable.name = gsub(pattern=".arg", replacement="", argument.name)
  variable.name = gsub(pattern=" ", replacement="", argument.name)
  argument.content = eval(parse(text = argument.name))
  eval(parse(text = argument.content))
  if (!exists(variable.name)){
    stop(sprintf("Argument %s not passed. Stopping ... ", variable.name))
  }
}

python.addr = 'python'

# create required folders for output and work material
output_folder = gsub(pattern="^\\d+_", replacement="", x=basename(getwd()))
output_folder = paste(output_folder, seurat.addr, sep = "_")
c.time = Sys.time()
c.time = gsub(pattern=" BST", replacement="", x=c.time)
c.time = gsub(pattern=":", replacement="", x=c.time)
c.time = gsub(pattern=" ", replacement="", x=c.time)
c.time = gsub(pattern="-", replacement="", x=c.time)
c.time = substr(x=c.time, start=3, stop=nchar(c.time))
output_folder = paste(output_folder, c.time, sep = "_")
output_folder = file.path("../../output", output_folder)
dir.create(output_folder)

output_folder_material = file.path(output_folder, "material")
dir.create(output_folder_material)

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(monocle)
library(dplyr)
library(reshape2)

#######################################################################################################
###########

ma = function(arr, kernel = 50){
  res = arr
  n = 2 * kernel
  for(i in 1:length(arr)){
    start_index = max(1, i - kernel)
    stop_index = min(length(arr), i + kernel)
    res[i] = mean(arr[start_index:stop_index])
  }
  res
}

adaptive.moving_average = function(arr, kernel = 10, minim_kernel = 10, range.factor = 5){
  res = arr
  n = 2 * kernel
  for(i in 1:length(arr)){
    start_index = max(1, i - kernel)
    stop_index = min(length(arr), i + kernel)
    local_sd = sd(arr[start_index:stop_index])
    local_kernel = minim_kernel + round(range.factor / (local_sd + .1))
    start_index = max(1, i - local_kernel)
    stop_index = min(length(arr), i + local_kernel)
    res[i] = mean(arr[start_index:stop_index])
  }
  res
}

###########
#######################################################################################################

print("Loading data ...")
seurat.addr = file.path("../../data", seurat.addr)
seurat.obj = readRDS(seurat.addr)
seurat.obj = SetAllIdent(object=seurat.obj, id=set.ident)
print("Data loaded.")

print("Subseting data on singlets and required cell populations")
if(cell.types == "all"){
  cell.types = as.vector(unique(seurat.obj@ident))
}

print("Subseting data ...")
to.keep = names(seurat.obj@ident)[as.vector(seurat.obj@ident) %in% cell.types]
seurat.obj = SubsetData(object=seurat.obj, cells.use=to.keep)
seurat.obj@ident = factor(seurat.obj@ident, levels = cell.types)

print("Writing data to disk ...")
# save raw data to disk
raw_data = seurat.obj@raw.data
raw_data = raw_data[rownames(seurat.obj@data), colnames(seurat.obj@data)]

to_exclude    = readRDS('cellcycle_genes.RDS')
genes_to_keep = rownames(raw_data)
genes_to_keep = genes_to_keep[!(genes_to_keep %in% to_exclude)]
raw_data = raw_data[genes_to_keep, colnames(seurat.obj@data)]

writeMM(raw_data, file.path(output_folder_material, "raw_data.mtx"))
# save gene names
gene_names = rownames(raw_data)
write.csv(data.frame(Genes = gene_names), file.path(output_folder_material, "genenames.csv"))
# save cell names
cell_names = colnames(raw_data)
write.csv(data.frame(Cells = cell_names), file.path(output_folder_material, "cellnames.csv"))
# write cell labels to disk
write.csv(data.frame(Cells = names(seurat.obj@ident), Labels = seurat.obj@ident), file.path(output_folder_material, "cell_labels.csv"), row.names = F)

print("Computing pseudotime...")
# compute pseudotime in python scanpy
command = sprintf("%s pdt_scanpy.py %s %s %s", python.addr, root_cell_type, output_folder, lineage.name)
system(command, wait=T)

# get cell labels and colours
if (!is.na(type.to.colours)){
  type.to.colours = file.path("../../resources", type.to.colours)
  type.to.colour = read.csv(type.to.colours)
  filter.key   =  type.to.colour$CellTypes %in% as.vector(unique(seurat.obj@ident))
  cell.labels  = as.vector(type.to.colour$CellTypes[filter.key])
  cell.colours = as.vector(type.to.colour$Colours[filter.key])
}else{
  cell.labels  = sort(as.vector(unique(seurat.obj@ident)))
  cell.colours = sample(colorRampPalette(brewer.pal(12, "Paired"))(length(cell.labels)))
}

# load pseudotime
print('reading pseudotime values')
pseudotime = read.csv(file.path(output_folder_material, "pseudotime.csv"), row.names = 1, header = F)
print("Are the cells in the same order in both pseudotime and seurat object? ")
print(all(rownames(pseudotime) == names(seurat.obj@ident)))
pseudotime$CellTypes = seurat.obj@ident
colnames(pseudotime) = c("Pseudotime", "CellType")

pseudotime$Color = mapvalues(x=pseudotime$CellType, from=cell.labels, to=cell.colours)
pseudotime$Color = factor(as.vector(pseudotime$Color), levels = cell.colours)
pseudotime$CellType = factor(as.vector(pseudotime$CellType), levels = cell.labels)
colnames(pseudotime) = c("Pseudotime", "Cell Type", "Color")

# compute diff genes
print("Computing var genes by cell type...")
cds                   = newCellDataSet(cellData = as.matrix(raw_data), phenoData=NULL, featureData=NULL, expressionFamily = negbinomial.size())
pData(cds)$Cluster    = as.vector(seurat.obj@ident)
cds                   = estimateSizeFactors(cds)
pData(cds)$Pseudotime = pseudotime$Pseudotime

var.genes.total = c()
print('Computing variable genes ... ')
for (j in 1:length(cell.labels)){
  print(sprintf("Choice %s out of %s ... ", as.character(j), as.character(length(cell.labels))))
  choices = pseudotime$`Cell Type` == cell.labels[j]
  var.genes = differentialGeneTest(cds[, choices], fullModelFormulaStr = "~sm.ns(Pseudotime)")
  var.genes = cbind(var.genes, data.frame(gene_id = rownames(var.genes)))
  var.genes = var.genes[var.genes$qval < .0001, ]
  var.genes.ch = var.genes %>% arrange(qval)
  var.genes.ch = as.vector(var.genes.ch$gene_id)
  var.genes.total = union(var.genes.total, var.genes.ch)
}

print("Computing var genes globally...")
var.genes = differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)")
var.genes = cbind(var.genes, data.frame(gene_id = rownames(var.genes)))
var.genes = var.genes[var.genes$qval < .0001, ]
var.genes.ch = as.vector(var.genes$gene_id)
var.genes.total = union(var.genes.total, var.genes.ch)
MT_genes = var.genes.total[grep("^MT-", x=var.genes.total, ignore.case=T)]
var.genes.total = setdiff(var.genes.total, MT_genes)
print(sprintf("Number of var genes total is : %d", length(var.genes.total)))
var.genes.total = sort(var.genes.total)

# min-max normalized expression
###################################################################################################
raw_data_genes = as.matrix(seurat.obj@data[var.genes.total, order(pseudotime$Pseudotime)])
raw_data_genes = t(apply(raw_data_genes, 1, adaptive.moving_average, kernel = 15, minim_kernel = 1, range.factor=15))

# min-max normalization
raw_data_genes_min = apply(raw_data_genes, 1, min)
raw_data_genes = raw_data_genes - raw_data_genes_min
raw_data_genes_max = apply(raw_data_genes, 1, max)
raw_data_genes = raw_data_genes / raw_data_genes_max

# non-normalized expression
###################################################################################################
raw_data_genes = as.matrix(seurat.obj@data[var.genes.total, order(pseudotime$Pseudotime)])
raw_data_genes = t(apply(raw_data_genes, 1, adaptive.moving_average, kernel = 15, minim_kernel = 1, range.factor=15))

# save diffusion map coordinates and expression data for found genes
by.pdt.order = order(pseudotime$Pseudotime)

dm.df = read.csv(file.path(output_folder_material, "dm.csv"), row.names = 1, header = F)
dm.df = as.data.frame(dm.df[, 2:4])

dm.df$Labels = factor(seurat.obj@ident, levels = cell.labels)
dm.df$Colours = mapvalues(x = dm.df$Labels, from = cell.labels, to = cell.colours)
dm.df = dm.df[by.pdt.order, ]
colnames(dm.df) = c("DM1", "DM2", "DM3", "Labels", "Colours")
print(head(dm.df))

expression_data_and_pdt = as.data.frame(t(as.matrix(seurat.obj@data[var.genes.total, by.pdt.order])))
pdt.data = data.frame(Pseudotime = pseudotime[by.pdt.order, c(1)])
pdt.data = cbind(dm.df, pdt.data, expression_data_and_pdt)
pdt.data.fp = file.path(output_folder, "pdt_and_expression.csv")

if (nrow(pdt.data) > 10000){
  sample.token = sort(sample(1:nrow(pdt.data), size=10000, replace=F))
  pdt.data = pdt.data[sample.token, ]
}

write.csv(pdt.data, pdt.data.fp, row.names = F)

# make interactive diffusion map
dir.create(file.path(output_folder, "genes"))
command = sprintf("%s pdt_3D_webportal.py %s %s %s", python.addr, output_folder, pdt.data.fp, lineage.name)
system(command, wait = T)

print("Ended beautifully ... ")
