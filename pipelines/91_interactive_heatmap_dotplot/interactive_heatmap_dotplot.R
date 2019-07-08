library(Seurat)
library(methods)

python.addr = 'python'

args = commandArgs(trailingOnly=T)
options_file = args[1]

options_fobj   = file(options_file, 'r')
options_fields = readLines(options_fobj)
close(options_fobj)

file_name      = options_fields[1]
set.ident      = options_fields[2]
output_folder  = options_fields[3]
save_to        = options_fields[4]
data_name      = options_fields[5]

dir.create(output_folder)

# check file_mame extension
# if is RDS assume this is a Seurat object and go on
# if is h5ad then assume it is a scanpy object and get help form Python to extract the required data
file_name_extension = unlist(strsplit(file_name, "\\."))
file_name_extension = file_name_extension[length(file_name_extension)]

if (file_name_extension == 'h5ad'){ # handle a scanpy object
  # must obtained an expression matrix aggregated by genes
  # rownames are cell types
  # colnames are genes
  print('Handling a Scanpy object')
  command = sprintf("%s aggregate_matrix_from_seurat.py %s %s %s", python.addr, file_name, set.ident, output_folder)
  system(command, wait = T)
  # read the expression sparse matrix from disk
  expression.data = readMM(file.path(output_folder, 'expression.mtx'))
  # convert the expression data to dgCMatrix so aggregation will be faster
  expression.data = as(expression.data, "dgCMatrix")
  # reading gene names from disk
  input_file = file(file.path(output_folder, 'gene_names.txt'))
  gene_names = readLines(input_file)
  close(input_file)
  # reading cell types from disk
  input_file = file(file.path(output_folder, "cell_types.txt"))
  cell_types = readLines(input_file)
  close(input_file)
  # update colnames and rownames of expression data
  rownames(expression.data) = cell_types
  colnames(expression.data) = gene_names
  # aggregate the expression matrix by gene
  no.genes = ncol(expression.data)
  start_index = 1
  while (start_index < no.genes){
    end_index = start_index + 999
    end_index = min(end_index, no.genes)
    expression.data_ = data.matrix(expression.data[, start_index:end_index])
    expression.data_ = as.data.frame(expression.data_)
    expression.data_ = cbind(data.frame(CellLabels = cell_types), expression.data_)
    expression.data_ = aggregate(expression.data_[2:dim(expression.data_)[2]], list(expression.data_$CellLabels), mean)
    expression.data_ = cbind(data.frame(CellType = expression.data_$Group.1), expression.data_[, 2:dim(expression.data_)[2]])
    rownames(expression.data_) = expression.data_$CellType
    expression.data_ = expression.data_[, 2:ncol(expression.data_)]
    print(start_index)
    if (start_index == 1){
      gene.expression.data = expression.data_
    }else{
      gene.expression.data = cbind(gene.expression.data, expression.data_)
    }
    start_index = start_index + 1000
  }
  # Save the expression matrix (aggregated by cell type using median) to the output folder
  write.csv(gene.expression.data, file.path(output_folder, 'expression.csv'))
  file.remove(file.path(output_folder, 'gene_names.txt'))
  file.remove(file.path(output_folder, 'cell_types.txt'))
}else{

  print("Loading data ... ")
  seurat.obj = readRDS(file_name)
  seurat.obj = SetAllIdent(object=seurat.obj, id=set.ident)

  # create expression data aggregated by the median expression of each gene in each cell population
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

  # Save the expression matrix (aggregated by cell type using median) to the output folder
  write.csv(expression.data, file.path(output_folder, 'expression.csv'))
}

# start the python script
command = sprintf('%s compile_app.py %s', python.addr, options_file)
system(command, wait = T)

# clean-up
file.remove(file.path(output_folder, 'expression.csv'))
file.remove(file.path(output_folder, 'expression.mtx'))

# end
print('Ended beautifully')
