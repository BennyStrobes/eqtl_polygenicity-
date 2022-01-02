




args <- commandArgs(TRUE)
input_file <- args[1]  # Input file
output_file <- args[2]  # Output file




load(input_file)

write.table(wgt.matrix[,1], file=output_file, quote=FALSE,sep='\t', col.names=FALSE)