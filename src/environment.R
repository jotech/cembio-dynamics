library(phyloseq)
library(ggplot2)
library(data.table)
library(stringr)
library(ggrepel)
library(ggpubr)

ps2 <- readRDS("../dat/phyloseq_ps2.RDS")
ps2 <- subset_samples(ps2, !time %in% c(0,2)) # sync time
levels(ps2@sam_data$source) <- c("control", "substrate", "host") # rename sample sources
ps2 <- subset_taxa(ps2, !id %in% c("unknown","OP50")) # remove non-cembio members
#ps2 <- transform_sample_counts(ps2, function(x){x / sum(x)})
ps2 <- transform_sample_counts(ps2, function(x){log(x+1)})

otu.melt.dt <- data.table(reshape2:::melt.matrix(ps2@otu_table, varnames=c("org","sample"), value.name="abundance"))
otu.melt.dt[,source:=ps2@sam_data$source[match(sample,ps2@sam_data$sample_id)]]
otu.melt.dt[,time:=ps2@sam_data$time[match(sample,ps2@sam_data$sample_id)]]

otu.melt.sum.dt <- otu.melt.dt[,list(abundance.med=median(abundance),abundance.mad=mad(abundance)),by=.(source,org)]

otu.wide.sum.dt <- dcast(otu.melt.sum.dt,org~source, value.var=c("abundance.med","abundance.mad"))
otu.wide.sum.dt <- otu.wide.sum.dt[rowSums(otu.wide.sum.dt[,.SD,.SDcols=is.numeric])>0] # remove zero abundant otus

ggplot(otu.wide.sum.dt, aes(x=abundance.med_host, y=abundance.med_control)) + geom_point(size=3) + geom_errorbar(alpha=0.3, aes(xmin=abundance.med_host-abundance.mad_host, xmax=abundance.med_host+abundance.mad_host)) + geom_errorbar(alpha=0.3, aes(ymin=abundance.med_control-abundance.mad_control, ymax=abundance.med_control+abundance.mad_control)) + theme_minimal(base_size=14) + geom_text_repel(aes(label=org)) + stat_cor(method="spearman") + geom_abline(intercept=0,slope=1, color="red", linetype="dashed") + xlab("Abundance (host)") + ylab("Abundance (control)")
ggsave("../img/abundance-cor_host-control.pdf", width=6,height=6)

ggplot(otu.wide.sum.dt, aes(x=abundance.med_host, y=abundance.med_substrate)) + geom_point(size=3) + geom_errorbar(alpha=0.3, aes(xmin=abundance.med_host-abundance.mad_host, xmax=abundance.med_host+abundance.mad_host)) + geom_errorbar(alpha=0.3, aes(ymin=abundance.med_substrate-abundance.mad_substrate, ymax=abundance.med_substrate+abundance.mad_substrate)) + theme_minimal(base_size=14) + geom_text_repel(aes(label=org)) + stat_cor(method="spearman") + geom_abline(intercept=0,slope=1, color="red", linetype="dashed") + xlab("Abundance (host)") + ylab("Abundance (substrate)")
ggsave("../img/abundance-cor_host-substrate.pdf", width=6,height=6)

ggplot(otu.wide.sum.dt, aes(x=abundance.med_substrate, y=abundance.med_control)) + geom_point(size=3) + geom_errorbar(alpha=0.3, aes(xmin=abundance.med_substrate-abundance.mad_substrate, xmax=abundance.med_substrate+abundance.mad_substrate)) + geom_errorbar(alpha=0.3, aes(ymin=abundance.med_control-abundance.mad_control, ymax=abundance.med_control+abundance.mad_control)) + theme_minimal(base_size=14) + geom_text_repel(aes(label=org)) + stat_cor(method="spearman") + geom_abline(intercept=0,slope=1, color="red", linetype="dashed") + xlab("Abundance (substrate)") + ylab("Abundance (control)")
ggsave("../img/abundance-cor_substrate-control.pdf", width=6,height=6)
