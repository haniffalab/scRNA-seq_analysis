library(Seurat)
library(methods)

python.addr = 'python'

args = commandArgs(trailingOnly=T)
options_file = args[1]

options_fobj   = file(options_file, 'r')
options_fields = readLines(options_fobj)
close(options_fobj)

file_name      = options_fields[1]
output_folder  = options_fields[2]
save_to        = options_fields[3]
data_name      = options_fields[4]

dir.create(output_folder)

# start the python script
command = sprintf('%s compile_app.py %s', python.addr, options_file)
system(command, wait = T)

# end
print('Ended beautifully')
