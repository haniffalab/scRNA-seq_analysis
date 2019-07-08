args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

print("The scanpy object will be converted to a seurat object")
print("Please keep in mind that the resulted seurat object will have some limitations.")
print("Any dimensionality reduction computed for the scanpy object will not be transfered to the seurat object.")
print("Furthermore the raw data will not be the original count data so you should be very very carefull subsetting or merging of the resulted seurat obj or just avoid doing these.")
print("All seurat objects saved from this pipeline will have a tag appended saying 'converted_from_scanpy' to raise awarness about the processing history of the data.")

arguments.list = "
scanpy.addr.arg = args[1]
save.to.arg     = args[2]
"

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
output_folder = "scanpy_to_seurat_"
c.time = Sys.time()
c.time = gsub(pattern=" BST", replacement="", x=c.time)
c.time = gsub(pattern=":", replacement="", x=c.time)
c.time = gsub(pattern=" ", replacement="", x=c.time)
c.time = gsub(pattern="-", replacement="", x=c.time)
c.time = substr(x=c.time, start=3, stop=nchar(c.time))
output_folder = paste(output_folder, c.time, sep = "")
output_folder = file.path("../../output", output_folder)
output_folder = paste(output_folder, paste(sample(c(letters, LETTERS), 10, F), collapse = ""), sep = "_")
dir.create(output_folder)

scanpy.addr = file.path("../../data", scanpy.addr)
save.to = paste(unlist(strsplit(save.to, "\\.RDS")), "converted_from_scanpy.RDS", sep = "_")

source("../../tools/bunddle_utils.R")

library(Seurat)
library(methods)
library(Matrix)

#######################################################################################################

# load data
command = sprintf("%s scanpy_to_seurat.py %s %s", python.addr, scanpy.addr, output_folder)
print("deconstructing the scanpy object in Python ...")
system(command, wait = T)

# load gene names
print("Loading gene names ...")
gene.names = read.csv(file.path(output_folder, "gene_names.csv"), header = F, row.names=1, stringsAsFactors = F)$V2

# load cell names
print("Loading cell names ...")
cell.names = read.csv(file.path(output_folder, "cell_names.csv"), header = F, row.names=1, stringsAsFactors = F)$V2

# load meta data
print("Loading meta data ... ")
meta.data = read.csv(file.path(output_folder, "meta_data.csv"), stringsAsFactors = F)

# load the expression matrix
print("Loading the expression matrix ...")
expression.data = readMM(file.path(output_folder, 'expression.mtx'))
expression.data = t(expression.data)

# make the seurat object
print("Constructing the Seurat object")
colnames(expression.data) = cell.names
rownames(expression.data) = gene.names
seurat.obj = CreateSeuratObject(raw.data = expression.data, min.cells = 0, min.genes = 0, project = "")
seurat.obj@meta.data = cbind(seurat.obj@meta.data, meta.data)
seurat.obj@data = expression.data

print("Saving the seurat object ... ")
saveRDS(seurat.obj, save.to)

unlink(output_folder, recursive=T, force=T)

print("Ended beautifully ... ")
