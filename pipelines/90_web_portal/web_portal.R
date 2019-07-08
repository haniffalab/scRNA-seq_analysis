args = commandArgs(trailingOnly=T)

python.addr = 'python'

# load the required libraries
library(Seurat)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(methods)

options_file = args[1]

##########################################################################################
##########################################################################################
##########################################################################################
# FUNCTIONS

# a function to convert hex colours to subunit rgb values - that are red by WebGL
hex_to_floats= function(hex_str){
  red = round(strtoi(paste("0x", substr(x=hex_str, start=2, stop=3), sep = "")) / 255, digits=2)
  green = round(strtoi(paste("0x", substr(x=hex_str, start=4, stop=5), sep = "")) / 255, digits=2)
  blue = round(strtoi(paste("0x", substr(x=hex_str, start=6, stop=7), sep = "")) / 255, digits=2)
  paste(c(red, green, blue), collapse = ",")
}

##########################################################################################
##########################################################################################
##########################################################################################

# extract options
options_fobj = file(options_file)
options_val  = readLines(options_fobj)
close(options_fobj)
# make file_name variable
file_name = gsub(" ", "", unlist(strsplit(options_val[1], ':'))[2])
# make output_folder variable
output_folder = gsub(" ", "", unlist(strsplit(options_val[2], ':'))[2])
# make dr coordinates pointers
dr_coordinates = gsub(" ", "", unlist(strsplit(options_val[3], ':'))[2])
dr_coordinates = unlist(lapply(unlist(strsplit(dr_coordinates, ';')), strsplit, "->"))
dr_names = dr_coordinates[seq(1, length(dr_coordinates), 2)]
dr_slots = dr_coordinates[seq(2, length(dr_coordinates), 2)]
# make categories pointers
categories       = unlist(lapply(unlist(strsplit(options_val[4], ';')), strsplit, '->'))
categories_names = categories[seq(1, length(categories), 3)]
categories_slot  = categories[seq(2, length(categories), 3)]
categories_cols  = categories[seq(3, length(categories), 3)]

# make the output folder and the other folders
dir.create(output_folder)
gene_output_folder = file.path(output_folder, "genes")
categories_output_folder = file.path(output_folder, "categories")
dr_output_folder = file.path(output_folder, "dr")
dir.create(gene_output_folder)
dir.create(categories_output_folder)
dir.create(dr_output_folder)

# check file_mame extension
# if is RDS assume this is a Seurat object and go on
# if is h5ad then assume it is a scanpy object and get help form Python to extract the required data
file_name_extension = unlist(strsplit(file_name, "\\."))
file_name_extension = file_name_extension[length(file_name_extension)]
if (file_name_extension == 'h5ad'){ # handle a scanpy object
  print('Handling a Scanpy object')
  command = sprintf("%s scanpy_to_seurat.py %s", python.addr, options_file)
  system(command, wait = T)
  # read expression matrix - should include nUMI and nGene
  expression_matrix = readMM(file.path(output_folder, 'expression.mtx'))
  # read the gene names
  gene_name_fobj = file(file.path(output_folder, 'gene_names.txt'))
  gene_names     = readLines(gene_name_fobj)
  close(gene_name_fobj)
  # add gene names as colnames to gene_expression
  colnames(expression_matrix) = gene_names
  # transpose the expression matrix
  expression_matrix = t(expression_matrix)
  # read dr data
  dr_matrix = read.csv(file.path(output_folder, 'dr.csv'), header = F)
  # read the categories
  categories_matrix = read.csv(file.path(output_folder, 'categories.csv'), row.names = 1, check.names=FALSE)
  # remove the files
  file.remove(file.path(output_folder, 'expression.mtx'))
  file.remove(file.path(output_folder, 'gene_names.txt'))
  file.remove(file.path(output_folder, 'dr.csv'))
  file.remove(file.path(output_folder, 'categories.csv'))
  # extract nUMI and nGENE
  nUMI  = expression_matrix['nUMI', ]
  nGene = expression_matrix['nGene', ]
  # remove nUMI and nGene rows from the expression matrix
  indices = which(rownames(expression_matrix) %in% c('nUMI', 'nGene'))
  expression_matrix = expression_matrix[-c(indices), ]
}else{ # handle a seurat object
  print('Handling a Seurat object')
  print('Loading data ...')
  seurat.obj = readRDS(file_name)
  print('Data loaded')
  # extract expression_matrix
  expression_matrix = seurat.obj@data
  # make the dr matrix
  for(i in 1:length(dr_slots)){
    dr_slot = dr_slots[i]
    command = paste("dr_slot_data = seurat.obj@dr$", dr_slot, "@cell.embeddings", sep = "")
    eval(parse(text = command))
    if (i == 1){
      dr_matrix = dr_slot_data
    }else{
      dr_matrix = cbind(dr_matrix, dr_slot_data)
    }
  }
  # make the categories_matrix
  for (i in 1:length(categories_slot)){
    cat_slot = categories_slot[i]
    command = paste("category_data = seurat.obj@meta.data$", cat_slot, sep = "")
    eval(parse(text = command))
    category_data = as.vector(category_data)
    if(categories_cols[i] == "null"){
      # randomly generate the colors
      unique_cats  = sort(as.vector(unique(category_data)))
      cat_colours = sample(colorRampPalette(brewer.pal(12, "Paired"))(length(unique_cats)))
    }else{
      type.to.colour = read.csv(categories_cols[i])
      unique_cats  = as.vector(type.to.colour$CellTypes)
      cat_colours = as.vector(type.to.colour$Colours)
    }
    cat_colours = mapvalues(x=category_data, from=unique_cats, to=cat_colours)
    cat_data = data.frame(Labels = category_data, Colours = cat_colours)
    colnames(cat_data) = c(categories_names[i], paste(categories_names[i], "_Colour", sep = ""))
    if ( i == 1 ){
      categories_matrix = cat_data
    }else{
      categories_matrix = cbind(categories_matrix, cat_data)
    }
  }
  # extract nUMI and nGENE
  nGene  = seurat.obj@meta.data$nGene
  nUMI = seurat.obj@meta.data$nUMI
}

# save color key for each category
# categories_names
# categories_slot
# categories_cols
for (i in 1:length(categories_names)){
  category_name = categories_names[i]
  category_col  = paste(category_name, 'Colour', sep = '_')
  category_slot = categories_slot[i]
  category_tags = as.vector(categories_matrix[, category_name])
  category_cols = as.vector(categories_matrix[, category_col])
  category_unique = unique(category_tags)
  cols_unique     = mapvalues(x=category_unique, from = category_tags, to = category_cols, warn_missing=F)
  # caol keys
  color.key = data.frame(Categories = category_unique, Colours = cols_unique)
  color_key_fname = paste(category_slot, "_color_key.csv", sep = '')
  write.csv(color.key, file.path(categories_output_folder, color_key_fname), row.names=FALSE)
  # color floats webgl-readable
  float_colors = unlist(lapply(category_cols, hex_to_floats))
  float_colors = paste(float_colors, collapse = ",")
  float_colors_file = paste(category_slot, 'colors', sep = '_')
  float_colors_file = file(file.path(categories_output_folder, float_colors_file), 'w')
  writeLines(float_colors, float_colors_file)
  close(float_colors_file)
  # save category indices
  cell.types.indices = ""
  for(i in 1:length(category_unique)){
    cell.label = category_unique[i]
    cell.label.indices  = which(category_tags %in% c(cell.label)) - 1
    cell.label.indices = paste(cell.label.indices, collapse = ",")
    assignment = paste(cell.label, cell.label.indices, sep = "->")
    cell.types.indices = c(cell.types.indices, assignment)
  }
  cell.types.indices = paste(cell.types.indices, collapse = ";")
  type_fname = paste(category_slot, 'indices', sep = '_')
  cell.types.indices.file = file.path(categories_output_folder, type_fname)
  cell.types.indices.file = file(cell.types.indices.file, "w")
  writeLines(cell.types.indices, cell.types.indices.file)
  close(cell.types.indices.file)
}

# save categories types key
# categories_names
# categories_slot
categories.menu.content = ""
for(i in 1:length(categories_names)){
  category_name = categories_names[i]
  category_slot = categories_slot[i]
  category_arrow = paste(category_name, "->", category_slot, ";", sep = '')
  categories.menu.content = paste(categories.menu.content, category_arrow, sep = "")
}
categories.menu.file = file.path(categories_output_folder, "categories_key")
categories.menu.file = file(categories.menu.file, "w")
writeLines(categories.menu.content, categories.menu.file)
close(categories.menu.file)

# save key to dr
#dr_names
#dr_slots
dr.menu.content = ""
for(i in 1:length(dr_names)){
  dr_name = dr_names[i]
  dr_slot = dr_slots[i]
  dr_arrow = paste(dr_name, "->", dr_slot, ";", sep = "")
  dr.menu.content = paste(dr.menu.content, dr_arrow, sep = "")
}
dr.menu.file = file.path(dr_output_folder, "dr_key")
dr.menu.file = file(dr.menu.file, "w")
writeLines(dr.menu.content, dr.menu.file)
close(dr.menu.file)

# save dr coordinates
#dr_names
#dr_slots
#dr_matrix
for(i  in 0:(length(dr_names) - 1)){
  dr_name = dr_names[i + 1]
  start_ = 2*i + 1
  end_   = 2*i + 2
  coordinates = dr_matrix[, c(start_, end_)]
  coordinates.x.q = mean(quantile(coordinates[, 1], c(.01, .99)))
  coordinates.y.q = mean(quantile(coordinates[, 2], c(.01, .99)))
  coordinates[, 1] = coordinates[, 1] - coordinates.x.q
  coordinates[, 2] = coordinates[, 2] - coordinates.y.q
  divide.x.at = quantile(coordinates[, 1], c(.99)) * 1.2
  divide.y.at = quantile(coordinates[, 2], c(.99)) * 1.2
  coordinates[, 1] = coordinates[, 1] / divide.x.at
  coordinates[, 2] = coordinates[, 2] / divide.y.at
  coordinates_ = c()
  for(i in 1:dim(coordinates)[1]){
    varx = round(coordinates[i, 1], digits = 3)
    vary = round(coordinates[i, 2], digits = 3)
    coordinates_ = c(coordinates_ ,paste(varx, vary, sep = ","))
  }
  coordinates = paste(coordinates_, collapse = ",")
  coordinates.file = paste(dr_name, "coordinates", sep = "_")
  coordinates.file = file.path(dr_output_folder, coordinates.file)
  coordinates.file = file(coordinates.file, "w")
  writeLines(coordinates, coordinates.file)
  close(coordinates.file)
}

# save gene expression
# get all genes names in alphabetical order
gene_names = read.csv("genes.tsv", sep = "\t", header = F)
gene_names = as.vector(unique(gene_names$V2))
gene_names = sort(gene_names)
# read the hgnc data frame and keep just the columns relevant for gene symbol and gene family
hgnc = read.csv("./hgnc.csv", sep = '\t', row.names =NULL)
hgnc = hgnc[, c("symbol", "gene_family")]
# filter out any rows where the gene symbol is not found in the seurat genes or where the gene symbol is not part of any family
hgnc = hgnc[hgnc$symbol %in% gene_names, ]
hgnc = hgnc[!(hgnc$gene_family == ""), ]
# extract all gene families and arrange them in alphabetical order
gene.families = sort(unique(unlist(strsplit(as.vector(hgnc$gene_family), split="\\|"))))
gene.families.contents = rep("", length(gene.families))
names(gene.families.contents) = gene.families
# population an empty vector which has names set to gene families, with all the genes part of that family
for(j in 1:dim(hgnc)[1]){
  gene.fam.concat =unlist(strsplit( as.vector(hgnc$gene_family)[j], split = "\\|"))
  quoted.gene.symbol = as.vector(hgnc$symbol)[j]
  for (k in 1:length(gene.fam.concat)){
    gene.fam = gene.fam.concat[k]
    if(gene.families.contents[gene.fam] != ""){
      gene.families.contents[gene.fam] = paste(gene.families.contents[gene.fam], paste("'", quoted.gene.symbol, "'", sep = ""), sep = ",")
    }else{
      gene.families.contents[gene.fam] = as.character(paste("'", quoted.gene.symbol, "'", sep = ""))
    }
  }
}
# insert each gene family into double quotes - for js compatibility
gene.families = unlist(lapply(gene.families, function(gene.fam){paste('\"', gene.fam, '\"', sep = "")}))
# insert each gene name in single quotes - for js compatibility
gene_names = unlist(lapply(gene_names, function(gene.n){paste("'", gene.n, "'", sep = "")}))
gene_names = paste(gene_names, collapse = ",")
# create a javacript point to gene_list and populate it with all gene names
gene.list = sprintf("gene_list = [%s];", gene_names)
# create the javacript named array gene_familie and populate it with indices for all genes for each gene family
gene.list = paste(gene.list, "gene_families = [];", sep = "")
# loop over each gene family name and insert its gene indices in gene.list
for(k in 1:length(gene.families)){
  gene.fam = gene.families[k]
  gene.indices = as.vector(gene.families.contents[k])
  gene.list = paste(gene.list, sprintf("gene_families[%s] = [%s];", gene.fam, gene.indices), sep = "")
}
gene.list.file = file.path(gene_output_folder, "gene_lists")
gene.list.file = file(gene.list.file, "w")
writeLines(gene.list, gene.list.file)
close(gene.list.file)

####### save gene expression files
gene_names = read.csv("genes.tsv", sep = "\t", header = F)
gene_names = as.vector(unique(gene_names$V2))
gene_names = sort(gene_names)
gene.expression.scale.limit = 1;
# loop though each gene name, get its expression and convert it to a color scale
# if gene name does not exists in the seurat object (has not passed filtering) save its expression as zeros
expression_matrix = as(expression_matrix, 'dgCMatrix')
for (i in 1:length(gene_names)){
  gene_name = gene_names[i]
  if (gene_name %in% rownames(expression_matrix)){
    gene_expression = expression_matrix[gene_name, ]
  }else{
    gene_expression = rep(0, ncol(expression_matrix))
  }
  gene_expression.upper.limit = max(max(gene_expression), gene.expression.scale.limit)
  gene_expression = gene_expression / gene_expression.upper.limit
  gene_expression = lapply(gene_expression, function(val){round(val, digits = 2)})
  gene_expression = unlist(gene_expression)
  gene_expression = c(round(gene_expression.upper.limit, digits=1), gene_expression)
  gene_expression = paste(gene_expression, collapse = ",")
  gene_expression.fname = file.path(gene_output_folder, gsub(pattern="/", replacement="__", x=gene_name));
  gene_expression.file = file(gene_expression.fname)
  writeLines(gene_expression, gene_expression.file)
  close(gene_expression.file)
  print(sprintf("%s - %s", gene_name, i))
}

# add the nUMI and nGene
for (j in 1:2){
  gene_name = c('nGene', 'nUMI')[j]
  if (j == 'nGene'){
    gene_expression = nGene
  }else{
    gene_expression = nUMI
  }
  gene_expression.upper.limit = max(gene_expression)
  gene_expression = gene_expression / gene_expression.upper.limit
  gene_expression = lapply(gene_expression, function(val){round(val, digits = 2)})
  gene_expression = unlist(gene_expression)
  gene_expression = c(round(gene_expression.upper.limit, digits=1), gene_expression)
  gene_expression = paste(gene_expression, collapse = ",")
  gene_expression.fname = file.path(gene_output_folder, gsub(pattern="/", replacement="__", x=gene_name));
  gene_expression.file = file(gene_expression.fname)
  writeLines(gene_expression, gene_expression.file)
  close(gene_expression.file)
  print(sprintf("%s - %s", gene_name, j))
}

# tar zip the output so it will be faster to download/upload between servers
print("Tar-zipping the output so it will be much faster for downloading/uploading operations between servers ...")
output_folder.arh = paste(output_folder, ".tar.gz", sep = "")
tar(output_folder.arh, files = output_folder, compression = 'gzip')
print("(The output is the tar.zip file which is much faster to download/upload than a folder with tens of thousands of files.)")

# remove the out folder
print("Deleting the output folder ...")
unlink(x=output_folder, recursive=T, force=T)

print("Ended beautifully ... ")