library(ggplot2)
library(ggrepel)
library(phyloseq)

ps2 <- readRDS("../dat/phyloseq_ps2.RDS")
ps2.sync <- subset_samples(ps2, !time %in% c(0,2)) # sync time
ps2rel.sync <- transform_sample_counts(ps2.sync, function(x){x / sum(x)})

ps2rel.sync.ord <- ordinate(ps2rel.sync, "MDS", "bray")
plot_ordination(ps2rel.sync, ps2rel.sync.ord, type="samples", color="time") + geom_point(size=3, shape=21, color="black", aes(fill=time)) + facet_wrap(~source) + scale_fill_brewer(type="qual", palette="Oranges") + theme_bw(base_size=14)
ggsave("../img/samples-mds_bray.pdf", width=8, height=2.5)
