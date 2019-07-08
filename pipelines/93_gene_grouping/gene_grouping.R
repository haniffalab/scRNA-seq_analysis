args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.addr.arg = args[1]
no_clusters.arg = args[2]
"

python.addr = "python"

expected_arguments = unlist(strsplit(arguments.list, "\n"))
expected_arguments = expected_arguments[!(expected_arguments == "")]

if(length(args) != length(expected_arguments)){
  error.msg = sprintf('This pipeline requires %s parameters', as.character(length(expected_arguments)))
  expected_arguments = paste(unlist(lapply(strsplit(expected_arguments, ".arg"), "[", 1)), collapse = "\n")
  stop(sprintf('This pipeline requires %s parameters: '))
}

eval(parse(text = arguments.list))

for(n in 1:length(expected_arguments)){
  argument = expected_arguments[n]
  argument = gsub(pattern=" ", replacement="", x=argument)
  argument.name = unlist(strsplit(argument, "="))[1]
  variable.name = gsub(pattern=".arg", replacement="", argument.name)
  argument.content = eval(parse(text = argument.name))
  eval(parse(text = argument.content))
  if (!exists(variable.name)){
    stop(sprintf("Argument %s not passed. Stopping ... ", variable.name))
  }
}

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

library(Seurat)
library(RColorBrewer)
library(dplyr)
library(plyr)

#######################################################################################################

# load data
print("loading data ... ")
seurat.obj = readRDS(seurat.addr)
print("Data loaded.")

# check if LouvainClustering is present
if ("LouvainClustering" %in% colnames(seurat.obj@meta.data)){
  print("Identifying gene outliers but first need to aggregate gene expression by clusters")
  seurat.obj = SetAllIdent(object=seurat.obj, id="LouvainClustering")
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
  # saving the expression matrix
  write.csv(expression.data, file.path(output_folder, "expression.csv"))
  # run python script to identify outliers
  command = sprintf("%s clustering.py %s %s", python.addr, output_folder, no_clusters)
  system(command, wait = T)
  # remove the expression csv file
  file.remove(file.path(output_folder, "expression.csv"))
  # load gene clustering
  gene_clustering = read.csv(file.path(output_folder, "clustering.csv"), row.names = 1)
  # save feature plots
  gene_names = as.vector(unique(gene_clustering$GeneNames))
  features_folder = file.path(output_folder, "features")
  dir.create(features_folder)
  dr_coordinates = seurat.obj@dr$umap@cell.embeddings
  for (i in seq_along(gene_names)){
    gene_name = gene_names[i]
    png_name = paste(file.path(features_folder, gene_name), "png", sep = ".")
    dframe = data.frame(X = dr_coordinates[, 1], Y = dr_coordinates[, 2], Expression = seurat.obj@data[gene_name, ])
    plot.obj = ggplot(dframe, aes(x = X, y = Y, color = Expression))
    plot.obj = plot.obj + geom_point(size = .5)
    plot.obj = plot.obj + theme_void() + theme(panel.background = element_rect(fill = 'black', colour = 'black'))
    plot.obj = plot.obj + scale_colour_gradient(low = "blue", high = "red")
    png(png_name, width = 500, height = 500)
    print(plot.obj)
    dev.off()
    if (i %% 10 == 0){
      print(sprintf("%s / %s", i, length(gene_names)))
    }
  }
}else{
  print("Data needs to be clustered first")
}

print("Ended beautifully ... ")