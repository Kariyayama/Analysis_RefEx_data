library(tidyverse)
source("211024_extact_tissue_expression.human_function.R")

## Backgroud
est_list <- "~/Labolatory/db/RefEx/211020/RefEx_expression_CAGE40_human_PRJDB1099.tsv"
est <- read.table(est_list, sep="\t", h=T, row.names=1)
est <- est[, !apply(est, 2, function(x) sum(as.numeric(x >= 0)) < length(x))]
print("EST dimension")
print(dim(est))

# vertebrate
vert_table <- preprocess("H_sapiens_vertebrate_class2.Refseq_id.release92.txt", est)
v_id<-read.table("../H_sapiens_vertebrate_class2.symbol_NCBI_Refseq_id.release92.txt", sep="\t", h=T)[,3:4]
v_id <- unique(v_id)
rownames(v_id) <- v_id[,2]
rownames(vert_table) <- v_id[rownames(vert_table),1]
vert_table.l <- table_longer(vert_table)
vert_table.c <- count_data(vert_table)

# gnathostome
gnath_table <- preprocess("H_sapiens_gnathostome_class2.Refseq_id.release92.txt", est)
g_id<-read.table("../H_sapiens_gnathostome_class2.symbol_NCBI_Refseq_id.release92.txt", sep="\t", h=T)[,3:4]
g_id <- unique(g_id)
rownames(g_id) <- g_id[,2]
rownames(gnath_table) <- g_id[rownames(gnath_table),1]
gnath_table.l <- table_longer(gnath_table)
gnath_table.c <- count_data(gnath_table)

# Highest expression
gnath_high <- rbind(rownames(gnath_table), colnames(gnath_table)[apply(gnath_table, 1, which.max)])

# Order of high expression
order_out <- c()
for(i in 1:nrow(gnath_table)){
    order_out <- rbind(order_out, colnames(sort(as.vector(gnath_table[i,]), decreasing=T)))
}
rownames(order_out) <- rownames(gnath_table)

# More than log2(TPM) > 1
morelog2 <- gnath_table > 1


