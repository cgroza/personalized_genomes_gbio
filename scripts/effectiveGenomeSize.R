# formula found in chipseq.py of mugqic_pipelines
args = commandArgs(TRUE)
file = args[1]
data = read.csv(file, sep = "\t", header = FALSE)
cat(sum(as.numeric(data$V2)) * 0.8)
