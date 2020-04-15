# script to load the metagenomics R functions

# Load the required packages
require(ggplot2)
require(dplyr)
require(data.table)
require(reshape2)
require(plyr)
require(phyloseq)

git_path <- dirname(sys.frame(1)$ofile)
git_path <-  paste(git_path, 'r_functions', sep = '/')

fcts <- list.files(path = git_path, pattern = "*.R")

message(sprintf("Source the following files from directory %s \n%s",git_path, paste(fcts,collapse = '\n')))

for(f in fcts){
  source(paste(git_path, f, sep = "/"))
}