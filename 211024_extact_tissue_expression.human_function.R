library(tidyverse)
library(RColorBrewer)

preprocess <- function(input, background){
    id<-read.table(input , h=T, sep="\t")
    print("ID dimension")
    print(dim(id))

   target <- background[as.character(unique(id[,4])),]
   id_na  <- target %>%
        apply(1, anyNA)

   target_table <- target[!id_na,]
   return(target_table)
}

table_longer <- function(tpm){
    # Expression
    tpm.r <- cbind(rownames(tpm), tpm)
    colnames(tpm.r)[1] <- "gene"
    target.l <- tpm.r %>%
        pivot_longer(cols=c(2:ncol(tpm.r)), names_to="Tissue", values_to="tpm")
    return(target.l)
}

expression_box <- function(type, tpm_data, output=TRUE){
    target_table.l <- table_longer(tpm_data)
    p1 <- ggplot(target_table.l, aes(x=Tissue, y=tpm)) +
            geom_boxplot() +
            coord_flip() +
            theme_light() +
            theme(panel.grid=element_blank(),
                  axis.text = element_text(size=16),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank())
    return(p1)
}

count_data <- function(tpm){
    # Absolute value
    count.tpm <- apply(tpm, 2, function(x) sum(as.numeric(x > 1)))
    count.tpm <- data.frame(count=count.tpm, Tissue=names(count.tpm))
    count.tpm <- count.tpm %>%
        as.data.frame
    return(count.tpm)
}

count_bar <- function(tpm_data){
    table_count <- count_data(tpm_data)
    p2 <- ggplot(table_count, aes(x=Tissue, y=count)) +
        geom_bar(stat="identity") +
        coord_flip() +
        theme_light() +
        theme(panel.grid=element_blank(),
              axis.text = element_text(size=16),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
    return(p2)
}

test_data <- function(target, bkground, method=chisq.test, alternative="t"){
    target_gene_number <- dim(target)[1]
    bkground_number <- dim(bkground)[1]
    target.count <- count_data(target)
    bkground.count <- count_data(bkground)

    pvalue_result <- c()
    for(i in 1:dim(target.count)[1]){
        sample_data1 <- matrix(c(target.count[i,1],    target_gene_number - target.count[i,1],
                                 bkground.count[i,1] - target.count[i,1],  bkground_number - target_gene_number - bkground.count[i,1]),
                               ncol=2, byrow=2)
        result <- method(sample_data1, alternative=alternative)
        pvalue_result <- c(pvalue_result, result$p.value)
    }
    names(pvalue_result) <- target.count[,2]
    return(pvalue_result)
}

utest_data <- function(target, bkground, method=wilcox.test, alternative="t"){
    pvalue_result <- c()
    for(i in 1:ncol(target)){
        result <- method(target[,i], bkground[,i], alternative=alternative)
        pvalue_result <- c(pvalue_result, result$p.value)
    }
    names(pvalue_result) <- colnames(target)
    return(pvalue_result)
}


## 21.10.27 追記
heatmap_data <- function(tpm_table, distfun="spearman", hclustfun="word.D2"){
    ## https://stats.biopapyrus.jp/r/ggplot/geom-tile.html
    tpm_table.l <- table_longer(tpm_table)

    # clustering
    clr <- heatmap(as.matrix(tpm_table), Rowv=T)

    gene.idx  <- rownames(data)[clr$rowInd]
    group.idx <- colnames(data)[clr$colInd]

    tpm_table.l$Gene  <- factor(tpm_table.l$gene, levels = gene.idx)
    tpm_table.l$Group <- factor(tpm_table.l$Tissue, levels = group.idx)

    ghm <- ggplot(tpm_table.l, aes(x = gene, y = Tissue, fill = tpm)) +
        geom_tile() +
        theme_bw() +
        theme(plot.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              strip.background = element_rect(fill = "white", colour = "white"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        xlab("Tissue") + ylab("Gene")
#        scale_fill_gradientn("tpm", colours = rev(brewer.pal(9, "Spectral")), na.value = "white") +
    return(ghm)
}

