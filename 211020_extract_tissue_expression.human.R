library(tidyverse)
source("211024_extact_tissue_expression.human_function.R")

## CAGE40
# background data
est_list <- "~/Labolatory/db/RefEx/211020/RefEx_expression_CAGE40_human_PRJDB1099.tsv"
est <- read.table(est_list, sep="\t", h=T, row.names=1)
est <- est[, !apply(est, 2, function(x) sum(as.numeric(x >= 0)) < length(x))]
print("EST dimension")
print(dim(est))

# vertebrate
vert_table <- preprocess("H_sapiens_vertebrate_class2.Refseq_id.release92.txt", est)
v_id<-read.table("../H_sapiens_vertebrate_class2.symbol_NCBI_Refseq_id.release92.txt",
                 sep="\t", h=T)[,3:4]
v_id <- unique(v_id)
rownames(v_id) <- v_id[,2]
rownames(vert_table) <- v_id[rownames(vert_table),1]
vert_test <- test_data(vert_table, est, method=fisher.test)
vert_table.l <- expression_box("vertebrate_name", vert_table)
vert_count <- count_bar("vertebrate_name", vert_table)
heatmap_data(vert_table, "211027_heatmap_vertebrate_name_CAGE40.png")

# gnathostome
gnath_table <- preprocess("H_sapiens_gnathostome_class2.Refseq_id.release92.txt", est)
g_id<-read.table("../H_sapiens_gnathostome_class2.symbol_NCBI_Refseq_id.release92.txt",
                 sep="\t", h=T)[,3:4]
g_id <- unique(g_id)
rownames(g_id) <- g_id[,2]
rownames(gnath_table) <- g_id[rownames(gnath_table),1]
gnath_test <- test_data(gnath_table, est, method=fisher.test)
gnath_table.l <- expression_box("gnathostome_name", gnath_table)
gnath_count <- count_bar("gnathostome_name", gnath_table)
heatmap_data(gnath_table, "211027_heatmap_gnathostome_name_CAGE40.png")


# ## all_human
# # background
# all_list <- "~/Labolatory/db/RefEx/211020/RefEx_expression_CAGE_all_human_PRJDB1099.tsv"
# all_human <- as.matrix(readr::read_tsv(all_list))
# gene <- all_human[,1]
# all_human <- all_human[,c(2:ncol(all_human))] %>%
#     as.data.frame
# rownames(all_human) <- gene
# 
# vert_table.all <- preprocess("H_sapiens_vertebrate_class2.Refseq_id.release92.txt", all_human)
# vert_test.all <- test_data(vert_table.all, all_human, method=fisher.test)
# vert_table.l.all <- expression_box("vertebrate", vert_table.all, output=F)
# vert_count.all <- count_bar("vertebrate", vert_table.all, output=F)
# 
# # gnathostome
# gnath_table.all <- preprocess("H_sapiens_gnathostome_class2.Refseq_id.release92.txt", all_human)
# gnath_test.all <- test_data(gnath_table.all, all_human, method=fisher.test)
# gnath_table.l.all <- expression_box("gnathostome", gnath_table.all)
# gnath_count.all <- count_bar("gnathostome", gnath_table.all)
