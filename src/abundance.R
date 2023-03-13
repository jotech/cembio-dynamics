library(phyloseq)
library(ggplot2)
library(ggrepel)

ps2 <- readRDS("../dat/phyloseq_ps2.RDS")

# barplot
ps2B <- merge_samples(ps2, "condition")
ps2B@sam_data$condition <- rownames(ps2B@sam_data)
ps2B@sam_data$source <- factor(str_extract(rownames(ps2B@sam_data), "host|associated|alone"))
levels(ps2B@sam_data$source) <- c("control", "substrate", "host") # rename sample sources
ps2B@sam_data$time <- factor(as.numeric(str_extract(rownames(ps2B@sam_data), "(?<=_)[0-9]+$")))
ps2Brel <- transform_sample_counts(ps2B, function(x){x / sum(x)})
ps2Brel@tax_table <- tax_table(cbind(ps2Brel@tax_table, matrix(ifelse(ps2Brel@tax_table[,8]=="unknown","x_other",ifelse(apply(ps2Brel@otu_table,2,max)>0.05, as.vector(ps2Brel@tax_table[,8]), "x_cembio (low)")), dimnames=list(NULL, "id2"))))		 
plot_bar(ps2Brel, x="time", y="Abundance", fill="id2") + facet_wrap(~source, scales="free") + scale_fill_manual(values=BacArena::colpal2) + xlab("Time [h]") + ylab("Relative abundance") + theme_minimal(base_size=14) + labs(fill="")
ggsave("../img/ASV-barplot.pdf", width=12, height=5)


# time series
colors.taxa <- fread("~/uni/cembio.ext/src/agnes/dat/colors_taxa.csv")
ps2rel <- transform_sample_counts(ps2, function(x){x / sum(x)})
levels(ps2rel@sam_data$source) <- c("control", "substrate", "host") # rename sample sources
ps2rel.dt <- melt(data.table(ps2rel@otu_table, keep.rownames=T), id.vars="rn", variable.name="sample_id", value.name="abundance")
setnames(ps2rel.dt, old="rn", new="org")
ps2rel.dt <- merge(ps2rel.dt, data.frame(ps2rel@sam_data), by="sample_id", all.x=T)
ps2rel.dt <- merge(ps2rel.dt, data.frame(ps2rel@tax_table[,6:8]), by.x="org", by.y="id", all.x=T)
ps2rel.dt <- merge(ps2rel.dt, colors.taxa, by.x="Genus", by.y="group", all.x=T)
ps2rel.sel.dt <- ps2rel.dt[!time%in%c(0,2) & org%in%c("MYb71,MYb49", "MYb10", "MYb388", "MYb328", "JUb44", "JUb134", "CEent1", "MYb396,MYb69,MYb21", "JUb19")]
ps2rel.sel.dt[time==42,label:=org]
ggplot(ps2rel.sel.dt, aes(x=time, y=abundance, color=org, group=org)) + stat_summary(geom="point", fun=mean, size=2) + stat_summary(geom="line", linewidth=1.1, fun=mean) + stat_summary(geom="errorbar", fun.data=mean_se, alpha=0.5, width=0.2) + stat_summary(geom="label_repel", fun=mean, aes(label = label), size=2.5, na.rm = TRUE, max.overlaps=15) + facet_wrap(~source) + scale_colour_manual(values=setNames(ps2rel.sel.dt$color, ps2rel.sel.dt$org)) + theme_bw(base_size=14) + theme(legend.position="none") + ylab("Relative abundance") + xlab("Time [h]")
ggsave("../img/ASV-timeseries.pdf", width=12, height=4)


