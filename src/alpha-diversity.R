library(ggplot2)
library(ggrepel)
library(phyloseq)
library(ggpubr)

ps2 <- readRDS("../dat/phyloseq_ps2.RDS")
ps2.sync <- subset_samples(ps2, !time %in% c(0,2)) # sync time


plot_richness(ps2.sync, x="time", measures=c("Shannon")) + geom_boxplot() + facet_wrap(~source) + theme_bw(base_size=14) + ylab("Alpha diversity (Shannon)")  + stat_compare_means(method="anova", label.x.npc="center") + stat_cor(method="spearman") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + geom_smooth()
ggsave("../img/alpha-diversity_shannon.pdf", width=8, height=2.5)

plot_richness(ps2.sync, x="time", measures=c("Chao1")) + geom_boxplot() + facet_wrap(~source) + theme_bw(base_size=14) + ylab("Alpha diversity (Chao1)") + stat_compare_means(method="anova", label.x.npc="center") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave("../img/alpha-diversity_chao1.pdf", width=8, height=2.5)

plot_richness(ps2.sync, x="time", measures=c("Simpson")) + geom_boxplot() + facet_wrap(~source) + theme_bw(base_size=14) + ylab("Alpha diversity (Simpson)") + stat_compare_means(method="anova", label.x.npc="center") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave("../img/alpha-diversity_simpson.pdf", width=8, height=2.5)


alpha.div.dt <- data.table(estimate_richness(ps2.sync, measures=c("Shannon", "Chao1", "Simpson")), source=ps2.sync@sam_data$source, time=ps2.sync@sam_data$time, sample=ps2.sync@sam_data$sample_name)

summary(aov(Shannon ~ source*time, data = alpha.div.dt))
summary(aov(Simpson ~ source*time, data = alpha.div.dt))
summary(aov(Chao1 ~ source*time, data = alpha.div.dt))

TukeyHSD(aov(Shannon ~ time, data = alpha.div.dt))
TukeyHSD(aov(Simpson ~ time, data = alpha.div.dt))
TukeyHSD(aov(Chao1 ~ time, data = alpha.div.dt))

TukeyHSD(aov(Shannon ~ source, data = alpha.div.dt))
TukeyHSD(aov(Simpson ~ source, data = alpha.div.dt))
TukeyHSD(aov(Chao1 ~ source, data = alpha.div.dt))

alpha.div.melt <- melt(alpha.div.dt, id.vars=c("source","time","sample"))[variable=="Shannon"]; alpha.div.dt[,time:=as.numeric(as.character(time))]
ggscatter(alpha.div.melt, x="time",y="value", add="reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int=T, facet.by="source") + stat_cor(method="spearman") + theme_minimal(base_size=14) + ylab("Alpha diversity (Shannon)") + scale_x_continuous(breaks=unique(alpha.div.melt$time))
ggsave("../img/alpha-diversity_shannon-trend.pdf", width=8, height=3)
