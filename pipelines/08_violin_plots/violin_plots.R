args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))
if(length(args) != 7){
  stop('This pipeline requires 7 parameters: seurat.addr, set.ident (from meta.data on which to set identity of cells,
       type.to.colours (file holding the type to colour key), cell.labels (list of categories or all), plot.width, 
       plot.height, features.file (name of file that holds the list of genes, should be placed in the resource folder)')
}

arguments.list = "
seurat.addr.arg     = args[1]
set.ident.arg       = args[2]
type.to.colours.arg = args[3]
cell.labels.arg     = args[4]
plot.width.arg      = args[5]   
plot.height.arg     = args[6]
features.file.arg   = args[7]
"
eval(parse(text = arguments.list))

arguments.list = unlist(strsplit(arguments.list, "\n"))
arguments.list = arguments.list[!(arguments.list == "")]

for(n in 1:length(arguments.list)){
  argument = arguments.list[n]
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
output_folder = paste("08_violin_plots", seurat.addr, sep = "_")
c.time = Sys.time()
c.time = gsub(pattern=" BST", replacement="", x=c.time)
c.time = gsub(pattern=":", replacement="", x=c.time)
c.time = gsub(pattern=" ", replacement="", x=c.time)
c.time = gsub(pattern="-", replacement="", x=c.time)
c.time = substr(x=c.time, start=3, stop=nchar(c.time))
output_folder = paste(output_folder, c.time, sep = "_")
output_folder = file.path("../../output", output_folder)
dir.create(output_folder)

seurat.addr = file.path("../../data", seurat.addr)
features.file   = file.path("../../resources", features.file)

features.file = file(features.file, "r")
features      = readLines(features.file)
close(features.file)

source("../../tools/bunddle_utils.R")

library(Seurat)
library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)
library(RColorBrewer)

#######################################################################################################

print("Loading data ...")
seurat.obj <- readRDS(seurat.addr)
print("Data loaded.")
seurat.obj <- SetAllIdent(object=seurat.obj, id=set.ident)
print("Data loaded.")
print("Features present in data:")
print(table(features %in% rownames(seurat.obj@data)))
features = intersect(features, rownames(seurat.obj@data))

if(cell.labels == "all"){
  cell.labels=as.vector(unique(seurat.obj@ident))
}else{
  cell.labels = file.path("../../resources", cell.labels)
  cell.labels.file = file(cell.labels, "r")
  cell.labels = readLines(cell.labels.file)
  close(cell.labels.file)
}

if (!is.na(type.to.colours)){
  type.to.colours = file.path("../../resources", type.to.colours)
  type.to.colour <- read.csv(type.to.colours)
  match.key <- match(cell.labels, type.to.colour$CellTypes)
  cell.colours <- as.vector(type.to.colour$Colours[match.key])
}else{
  cell.colours <- sample(colorRampPalette(brewer.pal(12, "Paired"))(length(cell.labels)))
}

to.keep = names(seurat.obj@ident)[as.vector(seurat.obj@ident) %in% cell.labels]
seurat.obj = SubsetData(object=seurat.obj, cells.use=to.keep)

plot_violins <- function(seurat.obj, class.order, class.colours, features){
  expression.data <- t(as.data.frame(as.matrix(seurat.obj@data[features, ])))
  expression.data <- melt(expression.data)
  colnames(expression.data) <- c("CellLabels", "Gene", "Expression")
  expression.data$CellLabels <- mapvalues(x=expression.data$CellLabels, from=names(seurat.obj@ident), to=as.vector(seurat.obj@ident))
  expression.data$CellLabels <- factor(as.vector(expression.data$CellLabels), levels = class.order)
  plot.obj <- ggplot(data=expression.data, aes(x=CellLabels, y = Expression, fill = CellLabels))
  plot.obj <- plot.obj + geom_violin(scale = 'width')
  plot.obj <- plot.obj + facet_wrap(~Gene, ncol=1)
  plot.obj <- plot.obj + theme(axis.text.x = element_text(angle = 90), legend.position = "none")
  plot.obj <- plot.obj +  scale_fill_manual(values=class.colours)
  plot.obj
}

f.name = file.path(output_folder, "features.pdf")
pdf(f.name, width = plot.width, height = plot.height)
print(plot_violins(seurat.obj, class.order=cell.labels, class.colours=cell.colours, features))
dev.off()

print("Ended beautifully ... ")
