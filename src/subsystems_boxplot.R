library(ggplot2)

ps2.asso.worm.feat.dds <- readRDS("../dat/dds-diff-asso.worm.feat.RDS")
ps2.alone.asso.feat.dds <- readRDS("../dat/dds-diff-alone.asso.feat.RDS")
ps2.alone.worm.feat.dds <- readRDS("../dat/dds-diff-alone.worm.feat.RDS")

ps2.combined.feat.dds <- rbind(data.table(cmp="associated vs. host", ps2.asso.worm.feat.dds[padj<=0.05]),
	  data.table(cmp="alone vs. associated", ps2.alone.asso.feat.dds[padj<=0.05]),
	  data.table(cmp="alone vs. host",ps2.alone.worm.feat.dds[padj<=0.05]))
ps2.combined.feat.dds[!is.na(hierarchy), subsystem:="metabolism"]
ggplot(ps2.combined.feat.dds, aes(x=subsystem, y=log2FoldChange)) + geom_boxplot() + coord_flip() + theme_minimal(base_size=14) + ylab("log2 fold change") + xlab("") + geom_hline(yintercept=0, linetype="dashed", color = "red") + scale_x_discrete(limits=rev) + facet_wrap(~cmp)
