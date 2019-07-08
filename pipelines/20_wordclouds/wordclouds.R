args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.addr.arg = args[1]
set.ident.arg   = args[2]
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
dir.create(file.path(output_folder, "gene_word_clouds"))
dir.create(file.path(output_folder, "celltype_word_clouds"))

seurat.addr = file.path("../../data", seurat.addr)

source("../../tools/bunddle_utils.R")

library(Seurat)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(wordcloud)


gene_to_weighted_cell_mention = function(gene.expr){
  idx = which(as.vector(gene_to_pop$V1) %in% names(gene.expr))
  gene.expr = gene.expr[as.vector(gene_to_pop$V1)[idx]]
  pop.expr = c()
  pop.names = c()
  for (k in 1:length(idx)){
    gene.name = names(gene.expr)[k]
    gene.value = gene.expr[k]
    pop.flags = as.vector(gene_to_pop$V2)[as.vector(gene_to_pop$V1) == gene.name]
    pop.flags = unlist(strsplit(pop.flags, ", "))
    for (p in 1:length(pop.flags)){
      pop.flag = pop.flags[p]
      gene.v = 100 * gene.value / populations.weight[pop.flag]
      if (pop.flag %in% pop.names){
        pop.expr[pop.flag] = pop.expr[pop.flag] + gene.v
      }else{
        pop.names = c(pop.names, pop.flag)
        pop.expr = c(pop.expr, gene.v)
        names(pop.expr) = pop.names
      }
    }
  }
  pop.expr
}

#######################################################################################################

# load data
print("loading data ... ")
seurat.obj = readRDS(seurat.addr)
print("Data loaded.")

seurat.obj = SetAllIdent(object=seurat.obj, id=set.ident)

print('Making gene clouds')
idents = as.vector(unique(seurat.obj@ident))
for (i  in 1:length(idents)){
  ident = idents[i]
  ident = names(seurat.obj@ident)[seurat.obj@ident == ident]
  expression.data = as.matrix(seurat.obj@data[,ident])
  expression.data = rowMeans(expression.data)
  genes = names(expression.data)
  genes = genes[!(genes %in% genes[grep(pattern='^MT-', x=genes)])]
  expression.data = expression.data[genes]
  fname = sprintf('%s/%s.pdf', file.path(output_folder, "gene_word_clouds"), gsub(pattern="/", replacement="-", x=idents[i]))
  pdf(fname, width = 15, height = 15)
  freq.weight = round(expression.data * 100)
  wordcloud(words=names(expression.data), freq.weight, min.freq = 1, max.words=500, 
            random.order=FALSE, rot.per=0.0, colors=brewer.pal(8, "Dark2"),
            order.color = T)
  dev.off()
}

print('Making cell type flags clouds')


expression.data = seurat.obj@data
mito.genes = grep(pattern="^MT-", x=rownames(expression.data))
expression.data = expression.data[-c(mito.genes), ]

gene_to_pop = read.csv("./gene_to_pop.tsv", sep = '\t', header = F)
populations = paste(as.vector(gene_to_pop$V2), collapse = ", ")
populations = unlist(strsplit(populations, ", "))
populations.table = table(populations)
populations.weight = as.vector(populations.table)
names(populations.weight) = names(populations.table)

idents = as.vector(unique(seurat.obj@ident))
for (i  in 1:length(idents)){
  ident = idents[i]
  print(ident)
  ident = names(seurat.obj@ident)[seurat.obj@ident == ident]
  expression.data = as.matrix(seurat.obj@data[,ident])
  expression.data = rowMeans(expression.data)
  genes = names(expression.data)
  genes = genes[!(genes %in% genes[grep(pattern='^MT-', x=genes)])]
  expression.data = expression.data[genes]
  pop.expr = gene_to_weighted_cell_mention(expression.data)
  clouder = round(100 * pop.expr)
  fname = sprintf('%s.pdf', idents[i])
  fname = gsub(pattern="/", replacement="-", x=fname)
  fname = file.path(file.path(output_folder, "celltype_word_clouds"), fname)
  pdf(fname, width = 10, height = 10)
  wordcloud(words=names(clouder), clouder, min.freq = 1, max.words=500, 
            random.order=FALSE, rot.per=0.0, colors=brewer.pal(8, "Dark2"),
            order.color = T)
  dev.off()
}

print("Ended beautifully ... ")
