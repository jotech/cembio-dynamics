library(ggplot2)
library(ggrepel)
library(phyloseq)
library(ggpubr)
library(foreach)
library(doParallel)
library(vegan)

ps2.tree <- readRDS("../dat/phyloseq_ps2.tree.RDS")
ps2.tree.sync <- subset_samples(ps2.tree, !time %in% c(0,2)) # sync time
ps2rel.tree.sync <- transform_sample_counts(ps2.tree.sync, function(x){x / sum(x)})
cores <- 2
dist.bray    = phyloseq::distance(ps2rel.tree.sync, method="bray")


# Ordination plots
ord.mds.bray <- ordinate(ps2rel.tree.sync, "MDS", "bray")

plot_ordination(ps2rel.tree.sync, ord.mds.bray, type="samples", color="time") + geom_point(size=3, shape=21, color="black", aes(fill=time)) + facet_wrap(~source) + scale_fill_brewer(type="qual", palette="Oranges") + theme_bw(base_size=14)
#ggsave("../img/ordination-mds_bray.pdf", width=8, height=2.5)


# PERMANOVA

anova(vegan::betadisper(dist.bray, sample_data(ps2rel.tree.sync)$source)) 
adonis2(dist.bray ~ sample_data(ps2rel.tree.sync)$source)

permanova.src.dt <- data.table()
for(src in levels(ps2.tree.sync@sam_data$source)){
    ps2.tmp <- subset_samples(ps2rel.tree.sync, source==src)
    bray.tmp    <- phyloseq::distance(ps2.tmp, method="bray")
    disp.test1 <- anova(vegan::betadisper(bray.tmp, sample_data(ps2.tmp)$time))
    permanova1 <- adonis2(bray.tmp ~ sample_data(ps2.tmp)$time)
    permanova.src.dt <- rbind(permanova.src.dt, data.table(src, dispersion.bray=disp.test1$`Pr(>F)`[1], permanova.bray=permanova1$`Pr(>F)`[1]))
}

permanova.time.dt <- data.table()
for(ti in levels(ps2.tree.sync@sam_data$time)){
    ps2.tmp <- subset_samples(ps2rel.tree.sync, time==ti)
    bray.tmp    <- phyloseq::distance(ps2.tmp, method="bray")
    disp.test1 <- anova(vegan::betadisper(bray.tmp, sample_data(ps2.tmp)$source))
    permanova1 <- adonis2(bray.tmp ~ sample_data(ps2.tmp)$source)
    permanova.time.dt <- rbind(permanova.time.dt, data.table(time=ti, dispersion.bray=disp.test1$`Pr(>F)`[1], permanova.bray=permanova1$`Pr(>F)`[1]))
}


# Compare beta diversity

dist_melt <- function(mydist){
	dist.melt <- data.table(reshape2::melt(as.matrix(mydist), value.name="dist"))
	dist.melt[,source1:=sample_data(ps2rel.tree.sync)$source[match(Var1, sample_data(ps2rel.tree.sync)$sample_id)]]
	dist.melt[,source2:=sample_data(ps2rel.tree.sync)$source[match(Var2, sample_data(ps2rel.tree.sync)$sample_id)]]
	dist.melt[,time1:=sample_data(ps2rel.tree.sync)$time[match(Var1, sample_data(ps2rel.tree.sync)$sample_id)]]
	dist.melt[,time2:=sample_data(ps2rel.tree.sync)$time[match(Var2, sample_data(ps2rel.tree.sync)$sample_id)]]
	dist.melt[,code1:=sample_data(ps2rel.tree.sync)$sample_code[match(Var1, sample_data(ps2rel.tree.sync)$sample_id)]]
	dist.melt[,code2:=sample_data(ps2rel.tree.sync)$sample_code[match(Var2, sample_data(ps2rel.tree.sync)$sample_id)]]
	dist.melt <- dist.melt[Var1!=Var2]
	dist.melt <- dist.melt[!duplicated(data.table(pmin(as.character(Var1),as.character(Var2)), pmax(as.character(Var1),as.character(Var2))))]
	return(dist.melt)
}
dist.bray.melt <- dist_melt(dist.bray)

dist.bray.source.time <-dist.bray.melt[time1==time2 & source1==source2]
ggplot(dist.bray.source.time, aes(x=time1, y=dist)) + geom_boxplot() + facet_wrap(~source1) + theme_bw(base_size=14) + xlab("Time [h]") + ylab("Beta-diversity (bray)") + stat_compare_means(method="anova", label.x.npc="center") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#ggsave("../img/beta-diversity_bray.pdf", width=8, height=2.5)



# Comparison pairs

dist.bray.pairs <- dist.bray.melt[code1==code2]
dist.bray.asso.alone <- dist.bray.melt[time1==time2 & ((source1=="associated" & source2=="alone") | (source1=="alone" & source2=="associated"))]
dist.bray.host.alone <- dist.bray.melt[time1==time2 & ((source1=="host" & source2=="alone") | (source1=="alone" & source2=="host"))]
#

registerDoParallel(cores)
rand_asso <- function(dist.pairs, dist.asso.alone, dist.host.alone){
	mycombine <- function(x, ...) { mapply(rbind,x,...,SIMPLIFY=FALSE) }
	rand.lst <- foreach(i=1:100, .combine="mycombine", .multicombine=TRUE) %dopar%{
		rand.asso.dt <- data.table()
		for(j in 1:nrow(dist.pairs)){
			time <- dist.pairs[j,time1]
			pair <- dist.pairs[j,code1]
			asso.alone <- dist.asso.alone[time1==time][sample(nrow(dist.asso.alone[time1==time]),1),]
			host.alone <- dist.host.alone[time1==time][sample(nrow(dist.host.alone[time1==time]),1),]
			rand.asso.dt <- rbind(rand.asso.dt, data.table(run=i,time, pair, dist.host.asso=dist.pairs[j,dist], dist.asso.alone=asso.alone$dist, dist.host.alone=host.alone$dist))
		}
		stat.host.asso <- compare_means(dist.host.asso~time, data=rand.asso.dt[run==i])
		stat.asso.alone <- compare_means(dist.asso.alone~time, data=rand.asso.dt[run==i])
		stat.host.alone <- compare_means(dist.host.alone~time, data=rand.asso.dt[run==i])

		list(data.table(run=i, rbind(stat.host.asso, stat.asso.alone, stat.host.alone)), rand.asso.dt)
	}
	return(rand.lst)
}
bray.rand.lst <- rand_asso(dist.bray.pairs, dist.bray.asso.alone, dist.bray.host.alone)
#saveRDS(bray.rand.lst, "../dat/beta-div_bray.rand-asso.RDS")

#bray.rand.lst <- readRDS("../dat/beta-div_bray.rand-asso.RDS")
bray.rand.asso.stat.dt <- bray.rand.lst[[1]]; bray.rand.asso.dt <- bray.rand.lst[[2]]
bray.rand.asso.stat.melt <- melt(bray.rand.asso.dt[,lapply(.SD,mean),by=.(pair,time),.SDcols=!grep("run",colnames(bray.rand.asso.dt))], id.vars=c("pair","time"))
bray.rand.asso.stat.melt[, cmp:=str_remove(variable, "^dist\\.")]
ggplot(bray.rand.asso.stat.melt, aes(x=time, y=value)) + geom_boxplot() + facet_wrap(~cmp) + theme_bw(base_size=14) + ylab("Beta-diversity (bray)") + stat_compare_means(method="anova", label.x.npc="center") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#ggsave("../img/beta-diversity_bray-pairs.pdf", width=8, height=2.5)

