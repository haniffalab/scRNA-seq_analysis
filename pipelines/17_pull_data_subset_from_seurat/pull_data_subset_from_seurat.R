args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.addr.arg   = args[1]
genes.file.arg    = args[2]
set.ident.arg     = args[3]
cell.types.arg    = args[4]
save.meta.data    = args[5]
save.raw.data.arg = args[6]
save.data.arg     = args[7]
save.dr.arg       = args[8]
"

expected_arguments = unlist(strsplit(arguments.list, "\n"))
expected_arguments = expected_arguments[!(expected_arguments == "")]

if(length(args) != length(expected_arguments)){
  error.msg = sprintf('This pipeline requires %s parameters', as.character(length(expected_arguments)))
  expected_arguments = paste(unlist(lappy(strsplit(expected_arguments, ".arg"), "[", 1)), collapse = "\n")
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

seurat.addr = file.path("../../data", seurat.addr)

library(Seurat)

#######################################################################################################

# load data
print("loading data ... ")
seurat.obj = readRDS(seurat.addr)
print("Data loaded.")
seurat.obj = SetAllIdent(object=seurat.obj, id=set.ident)

if(!is.na(genes.file)){
  genes.file = file.path("../../resources", genes.file)
  genes.file = file(genes.file, "r")
  genes      = readLines(genes.file)
  genes.raw = genes
  close(genes.file)
}else{
  genes = rownames(seurat.obj@data)
  genes.raw = rownames(seurat.obj@raw.data)
}

if(cell.types == 'all'){
  cell.types = as.vector(unique(seurat.obj@ident))
}
cells.to.keep = names(seurat.obj@ident)[seurat.obj@ident %in% cell.types]

if(save.meta){
  meta.data = as.data.frame(seurat.obj@meta.data)[cells.to.keep, ]
  write.csv(meta.data, file.path(output_folder, "meta_data.csv"))
}

if(save.raw.data){
  saveRDS(seurat.obj@raw.data[genes.raw, cells.to.keep], file.path(output_folder, "raw_data.RDS"))
}

if(save.data){
  saveRDS(as.data.frame(as.matrix(seurat.obj@data[genes, cells.to.keep])), file.path(output_folder, "expression_data.RDS"))
}

if (save.dr){
  dr.list = list()
  for (i in seq_along(names(seurat.obj@dr))){
    dr.type = names(seurat.obj@dr)[i]
    eval(parse(text = sprintf("dr.data = seurat.obj@dr$%s@cell.embeddings", dr.type)))
    dr.list[[dr.type]] = dr.data
  }
  saveRDS(dr.list, file.path(output_folder, "dr_data.RDS"))
}

#######################################################################################################

print("Ended beautifully ... ")
