#!/usr/bin/env Rscript

###### Bradley Nov 2023
library(optparse)
library(dplyr)


# Inherit options
option_list = list(
    make_option(c("-c", "--catdir"), action = "store", default = NA, type ="character",
        help="Where are the files located?")
    )
opt = parse_args(OptionParser(option_list = option_list))
catdir = opt$c

# Load files and process
nPCs = c("5","10","15","20")
res = data.frame("nGenotype_PCs"=nPCs, "neGenes"="")
cols = c("gene", "ACAT_p", "top_MarkerID", "top_pval")
for(f in nPCs){
    temp = read.delim(paste0(catdir, "/chr1_nPC_", f, "_ACAT_all.txt"), header=F)
    temp = temp[-1,]
    colnames(temp) = cols
    nsig = nrow(temp[temp$top_pval < 5e-8,])
    res$neGenes[which(nPCs == f)] = nsig
    rm(temp)
}

# Save the value of the nPCs at which the maximum number of 
max_res = res[res$neGenes == max(res$neGenes),]$nGenotype_PCs
if(length(max_res) > 1){
    max_res = max_res[1]
}
write.table(max_res, paste0(catdir, "/optim_nPCs_chr1.txt"), row.names=F, col.names=F, quote=F)
