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
v_id<-read.table("../H_sapiens_vertebrate_class2.symbol_NCBI_Refseq_id.release92.txt", sep="\t", h=T)[,3:4]
v_id <- unique(v_id)
rownames(v_id) <- v_id[,2]
rownames(vert_table) <- v_id[rownames(vert_table),1]

vert_test <- test_data(vert_table, est, method=fisher.test)
vert_table.l <- table_longer(vert_table)
vert_table.c <- count_data(vert_table)

vert_box <- expression_box(vert_table)
ggsave("H_sapiens_vertebrate_class2.CAGE_40.expression_boxplot.png",
       vert_box, units="px", width=4200, height=2100)
vert_bar <- count_bar("vertebrate_name", vert_table)
ggsave("H_sapiens_vertebrate_class2.CAGE_40.expression_barplot.png",
       vert_bar, units="px", width=4200, height=2100)
vert_heat <- heatmap_data(vert_table)
ggsave("211027_heatmap_vertebrate_name_CAGE40.png", vert_heat)

# gnathostome
gnath_table <- preprocess("H_sapiens_gnathostome_class2.Refseq_id.release92.txt", est)
g_id<-read.table("../H_sapiens_gnathostome_class2.symbol_NCBI_Refseq_id.release92.txt",
                 sep="\t", h=T)[,3:4]
g_id <- unique(g_id)
rownames(g_id) <- g_id[,2]
rownames(gnath_table) <- g_id[rownames(gnath_table),1]

gnath_test <- test_data(gnath_table, est, method=fisher.test)
gnath_table.l <- table_longer(gnath_table)
gnath_table.c <- count_data(gnath_table)

gnath_box <- expression_box(gnath_table)
ggsave("H_sapiens_gnathostome_class2.CAGE_40.expression_boxplot.png",
       gnath_box, units="px", width=4200, height=2100)
gnath_bar <- count_bar("vertebrate_name", vert_table)
ggsave("H_sapiens_gnathostome_class2.CAGE_40.expression_barplot.png",
       gnath_bar, units="px", width=4200, height=2100)
gnath_heat <- heatmap_data(gnath_table)
ggsave("211027_heatmap_gnathostome_name_CAGE40.png", gnath_heat)
