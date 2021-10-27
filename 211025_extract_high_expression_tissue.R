source("211024_extact_tissue_expression.human_function.R")

## all_human
# background
all_list <- "~/Labolatory/db/RefEx/211020/RefEx_expression_CAGE_all_human_PRJDB1099.tsv"
all_human <- as.matrix(readr::read_tsv(all_list))
gene <- all_human[,1]
all_human <- all_human[,c(2:ncol(all_human))] %>%
    as.data.frame
rownames(all_human) <- gene

each_cell_median <- apply(all_human, 2, median)
vert_table.all <- preprocess("H_sapiens_vertebrate_class2.Refseq_id.release92.txt", all_human)
# vert_test.all <- test_data(vert_table.all, all_human, method=fisher.test)
# vert_table.l.all <- expression_box("vertebrate", vert_table.all, output=F)
# vert_count.all <- count_bar("vertebrate", vert_table.all, output=F)
# 
# gnathostome
gnath_table.all <- preprocess("H_sapiens_gnathostome_class2.Refseq_id.release92.txt", all_human)
# gnath_test.all <- test_data(gnath_table.all, all_human, method=fisher.test)
# gnath_table.l.all <- expression_box("gnathostome", gnath_table.all)
# gnath_count.all <- count_bar("gnathostome", gnath_table.all)
