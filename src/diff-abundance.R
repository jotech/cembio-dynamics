library(DESeq2)
library(phyloseq)
library(ggplot2)

ps2 <- readRDS("../dat/phyloseq_ps2.RDS")

deseq_diffspec <- function(ps2){
  ps2dds <- phyloseq_to_deseq2(ps2, ~ time + source)
  gm_mean = function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
  geoMeans = apply(counts(ps2dds), 1, gm_mean)
  ps2dds <- estimateSizeFactors(ps2dds, geoMeans=geoMeans); rownames(ps2dds) <- ps2@tax_table[,8]
  ps2dds <- estimateDispersions(ps2dds)
  ps2dds <- DESeq(ps2dds, fitType="local")
  ps2dds_res <- as.data.table(results(ps2dds)); ps2dds_res[,id:=ps2@tax_table[,8]]
	return(ps2dds_res)
}




ps2.asso.worm  <- subset_samples(ps2, source %in% c("host","associated"))
ps2dds.asso.worm  <- deseq_diffspec(ps2.asso.worm)
ggplot(ps2dds.asso.worm[padj<=0.05], aes(x=as.character(id), y=log2FoldChange, fill=log(baseMean))) + geom_bar(stat="identity") + coord_flip() + theme_minimal(base_size=14) + ylab("log2FC associated vs. host") + xlab("") + scale_fill_gradient(low="gray", high="steelblue")
ggsave("~/uni/cembio.ext/img/agnes/dds-diff-otu-asso.worm.pdf", height=2, width=6)


ps2.alone.asso <- subset_samples(ps2, source %in% c("alone","associated") & !time %in% c(0, 2)) # 2 more time points in cembio.alone
ps2dds.alone.asso <- deseq_diffspec(ps2.alone.asso)
ggplot(ps2dds.alone.asso[padj<=0.05], aes(x=as.character(id), y=log2FoldChange, fill=log(baseMean))) + geom_bar(stat="identity") + coord_flip() + theme_minimal(base_size=14) + ylab("log2FC alone vs. associated") + xlab("") + scale_fill_gradient(low="gray", high="steelblue")
ggsave("~/uni/cembio.ext/img/agnes/dds-diff-otu-alone.asso.pdf", height=2, width=6)

ps2.alone.worm <- subset_samples(ps2, source %in% c("alone","host") & !time %in% c(0, 2)) # 2 more time points in cembio.alone
ps2dds.alone.worm <- deseq_diffspec(ps2.alone.worm)
ggplot(ps2dds.alone.worm[padj<=0.05], aes(x=as.character(id), y=log2FoldChange, fill=log(baseMean))) + geom_bar(stat="identity") + coord_flip() + theme_minimal(base_size=14) + ylab("log2FC alone vs. host") + xlab("") + scale_fill_gradient(low="gray", high="steelblue")
ggsave("~/uni/cembio.ext/img/agnes/dds-diff-otu-alone.worm.pdf", height=2, width=6)
