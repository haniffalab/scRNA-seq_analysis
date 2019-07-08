args = commandArgs(trailingOnly=T)
option.file = args[1]

# create required folders for output and work material
output_folder = gsub(pattern="^\\d+_", replacement="", x=basename(getwd()))
c.time = Sys.time()
c.time = gsub(pattern=" BST", replacement="", x=c.time)
c.time = gsub(pattern=":", replacement="", x=c.time)
c.time = gsub(pattern=" ", replacement="", x=c.time)
c.time = gsub(pattern="-", replacement="", x=c.time)
c.time = substr(x=c.time, start=3, stop=nchar(c.time))
output_folder = paste(output_folder, c.time, sep = "_")
output_folder = file.path("../../output", output_folder)
dir.create(output_folder)

source("../../tools/bunddle_utils.R")

library(Seurat)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(magrittr)

#######################################################################################################

# parse options
option.file = file(option.file, "r")
option_lines = readLines(option.file)
close(option.file)
#remove comments
option_lines = option_lines[-grep(pattern="^#", x=option_lines)]
# split options into blocks
option_lines =paste(option_lines, collapse = "@@@")
option_blocks = option_lines %>% strsplit(split="name: ") %>% unlist
option_blocks = option_blocks[option_blocks != ""]

# loop through each block of options, load the data sets and make the AGAs
for(k in seq_along(option_blocks)){
  option_block = option_blocks[k] %>% strsplit(split="@@@") %>% unlist
  print(sprintf("Making AGA for %s", option_block[1]))
  AGA_save_to = file.path(output_folder, option_block[1])
  dir.create(AGA_save_to)
  output_folder_material = file.path(AGA_save_to, "material")
  AGA_folder = file.path(AGA_save_to, "AGA_folder")
  dir.create(output_folder_material)
  dir.create(AGA_folder)
  # load data sets
  data_files = option_block[grep(pattern="^data", x=option_block)] %>% gsub(pattern="^data.: ", replacement="")
  set.idents = option_block[grep(pattern="^set.ident", x=option_block)] %>% gsub(pattern="^set.ident..: ", replacement="")
  label.tags = option_block[grep(pattern="^tag", x=option_block)] %>% gsub(pattern="^tag.: ", replacement="")
  categories = option_block[grep(pattern="^categories", x=option_block)]
  
  data.list = list()
  for(i in seq_along(data_files)){
    data.file = file.path("../../data", data_files[i])
    print(sprintf("Loading %s", data.file))
    seurat.obj = readRDS(data.file)
    seurat.obj %<>% SetAllIdent(id=set.idents[i])
    cell.labels = categories[i] %>% strsplit(split=": ") %>% unlist %>% (function(x)x[2]) %>% strsplit(split=", ") %>% unlist
    seurat.obj %<>% SubsetData(ident.use=cell.labels)
    seurat.obj@meta.data$AGA_labels = paste(paste(label.tags[i], "::", sep=""), as.vector(seurat.obj@ident), sep="")
    eval(parse(text=sprintf("data.list$data%s = seurat.obj", i)))
  }
  print("Subsetting and merging datasets ...")
  seurat.obj = Reduce(f=MergeSeurat, x=data.list)
  seurat.obj %<>% SetAllIdent(id="AGA_labels")
  
  write.csv(data.frame(Cells = names(seurat.obj@ident), Labels = seurat.obj@ident), file.path(output_folder_material, "cell_labels.csv"), row.names = F)
  
  # save raw data to disk
  raw_data = seurat.obj@raw.data
  raw_data = raw_data[rownames(seurat.obj@data), colnames(seurat.obj@data)]
  writeMM(raw_data, file.path(output_folder_material, "raw_data.mtx"))
  # save gene names
  gene_names = rownames(raw_data)
  write.csv(data.frame(Genes = gene_names), file.path(output_folder_material, "genenames.csv"))
  # save cell names
  cell_names = colnames(raw_data)
  write.csv(data.frame(Cells = cell_names), file.path(output_folder_material, "cellnames.csv"))
  # write cell labels to disk
  write.csv(data.frame(Cells = names(seurat.obj@ident), Labels = seurat.obj@ident), file.path(output_folder_material, "cell_labels.csv"), row.names = F)
  
  # running AGA
  command = file.path(tool_addr, "AGA/AGA_from_Seurat.py")
  command = paste(paste(python.addr, command, sep = " "), AGA_save_to, sep = " ")
  command = paste(command, option_block[1], sep =" ")
  system(command, wait = T)
  
  # read the AGA output
  coordinates = read.csv(file.path(AGA_folder, "coordinates.csv"), row.names = 1)
  connectivities = read.csv(file.path(AGA_folder, "connectivities.csv"), row.names = 1)
  colnames(connectivities) = rownames(connectivities)
  cell.labels = rownames(coordinates)
  cell.colours = sample(colorRampPalette(brewer.pal(12, "Paired"))(length(cell.labels)))
  
  ######## now make the interactive AGA app 
  #########################################
  print("Making the AGA app ... ")
  # save colours
  colours.df = data.frame(CellTypes = cell.labels, Colours = cell.colours)
  write.csv(colours.df, file.path(AGA_folder, "colours.csv"), row.names = F)
  
  # run python to built the AGA app
  command = sprintf("%s make_AGA_app.py %s %s", python.addr, AGA_save_to, option_block[1])
  system(command, wait = T)
  
  AGA.file = file.path(AGA_save_to, paste(c("AGAlinkage_map_", option_block[1], ".html"), collapse = ""))
  AGA.final.destination = file.path(output_folder, paste(c("AGAlinkage_map_", option_block[1], ".html"), collapse = ""))
  file.rename(from=AGA.file, to=AGA.final.destination)
  
  # make FDG
  print("Making force directed graph interactive app ...")
  seurat.obj = FindVariableGenes(object = seurat.obj, mean.function = ExpMean, 
                                 dispersion.function = LogVMR, x.low.cutoff = .0125, 
                                 x.high.cutoff = 3, y.cutoff = .625, do.plot=F)
  seurat.obj = ScaleData(object=seurat.obj)
  seurat.obj = RunPCA(object = seurat.obj, pc.genes = seurat.obj@var.genes, do.print = FALSE)
  seurat.obj = BuildSNN(object=seurat.obj, reduction.type="pca", dims.use=1:20, plot.SNN=F, force.recalc=T)
  fdg_coordinates = runFDG(pca.df=seurat.obj@dr$pca@cell.embeddings, snn=seurat.obj@snn, iterations=2000, tool_addr=tool_addr, python.addr=python.addr)
  seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="fdg", slot="cell.embeddings", new.data=as.matrix(fdg_coordinates))
  seurat.obj = SetDimReduction(object=seurat.obj, reduction.type="fdg", slot = "key", new.data = "fdg")
  
  interactive_plot_df = data.frame(X = seurat.obj@dr$fdg@cell.embeddings[, 1],
                                   Y = seurat.obj@dr$fdg@cell.embeddings[, 2])
  interactive_plot_df$Labels = factor(seurat.obj@ident, levels = cell.labels)
  interactive_plot_df$Colours = mapvalues(x = interactive_plot_df$Labels, from = cell.labels, to = cell.colours)
  
  interactive_fdg_filename = file.path(output_folder, paste(paste("Interactive_FDG", option_block[1], sep = "_"), "html", sep = "."))
  make_2D_interactive_page(data_frame_2D=interactive_plot_df, tool_addr=tool_addr, python.addr=python.addr, save.to=interactive_fdg_filename)
  
  unlink(AGA_save_to, recursive=T, force=T)
}

print("Ended beautifully ... ")
