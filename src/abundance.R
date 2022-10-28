library(phyloseq)
library(ggplot2)

ps2 <- readRDS("../dat/phyloseq_ps2.RDS")

ps2B <- merge_samples(ps2, "condition")
ps2B@sam_data$condition <- rownames(ps2B@sam_data)
ps2B@sam_data$source <- str_extract(rownames(ps2B@sam_data), "host|associated|alone")
ps2B@sam_data$time <- factor(as.numeric(str_extract(rownames(ps2B@sam_data), "(?<=_)[0-9]+$")))
ps2Brel <- transform_sample_counts(ps2B, function(x){x / sum(x)})
ps2Brel@tax_table <- tax_table(cbind(ps2Brel@tax_table, matrix(ifelse(ps2Brel@tax_table[,8]=="unknown","x_other",ifelse(apply(ps2Brel@otu_table,2,max)>0.05, as.vector(ps2Brel@tax_table[,8]), "x_cembio (low)")), dimnames=list(NULL, "id2"))))		 
plot_bar(ps2Brel, x="time", y="Abundance", fill="id2") + facet_wrap(~source, scales="free") + scale_fill_manual(values=BacArena::colpal3) + xlab("Time [h]") + ylab("Relative abundance") + theme_minimal(base_size=14) + labs(fill="")
ggsave("../img/ASV-barplot.pdf", width=12, height=5)

