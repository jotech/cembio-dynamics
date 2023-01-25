library(ggplot2)
library(ggrepel)
library(phyloseq)
library(ggpubr)
library(foreach)
library(doParallel)
library(vegan)
library(microViz)

cores <- 2
ps2.tree <- readRDS("../dat/phyloseq_ps2.tree.RDS")
ps2.tree.sync <- subset_samples(ps2.tree, !time %in% c(0,2)) # sync time
levels(ps2.tree.sync@sam_data$source) <- c("control", "substrate", "host") # rename sample sources
ps2rel.tree.sync <- transform_sample_counts(ps2.tree.sync, function(x){x / sum(x)})

dist.bray      <- microViz::dist_calc(ps2rel.tree.sync, dist="bray")@dist
dist.aitchison <- microViz::dist_calc(ps2rel.tree.sync, dist="aitchison")@dist
dist.unifrac   <- microViz::dist_calc(ps2rel.tree.sync, dist="wunifrac")@dist


# Ordination plots
ord.mds.bray <- ordinate(ps2rel.tree.sync, method="PCoA", distance=dist.bray)
plot_ordination(ps2rel.tree.sync, ord.mds.bray, type="samples", color="time") + geom_point(size=3, shape=21, color="black", aes(fill=time)) + facet_wrap(~source) + scale_fill_brewer(type="qual", palette="Oranges") + theme_bw(base_size=14)
#ggsave("../img/ordination-pcoa_bray.pdf", width=8, height=2.5)

ord.mds.aitchison <- ordinate(ps2rel.tree.sync, method="PCoA", distance=dist.aitchison)
plot_ordination(ps2rel.tree.sync, ord.mds.aitchison, type="samples", color="time") + geom_point(size=3, shape=21, color="black", aes(fill=time)) + facet_wrap(~source) + scale_fill_brewer(type="qual", palette="Oranges") + theme_bw(base_size=14)
#ggsave("../img/ordination-pcoa_aitchison.pdf", width=8, height=2.5)

ord.mds.unifrac <- ordinate(ps2rel.tree.sync, method="PCoA", distance=dist.unifrac)
plot_ordination(ps2rel.tree.sync, ord.mds.unifrac, type="samples", color="time") + geom_point(size=3, shape=21, color="black", aes(fill=time)) + facet_wrap(~source) + scale_fill_brewer(type="qual", palette="Oranges") + theme_bw(base_size=14)
#ggsave("../img/ordination-pcoa_unifrac.pdf", width=8, height=2.5)


# PERMANOVA

anova(vegan::betadisper(dist.bray, sample_data(ps2rel.tree.sync)$source)) 
adonis2(dist.bray ~ sample_data(ps2rel.tree.sync)$source, permutations=10000)

anova(vegan::betadisper(dist.aitchison, sample_data(ps2rel.tree.sync)$source)) 
adonis2(dist.aitchison ~ sample_data(ps2rel.tree.sync)$source, permutations=10000)

anova(vegan::betadisper(dist.unifrac, sample_data(ps2rel.tree.sync)$source)) 
adonis2(dist.unifrac ~ sample_data(ps2rel.tree.sync)$source, permutations=10000)

permanova.src.dt <- data.table()
for(src in levels(ps2.tree.sync@sam_data$source)){
    ps2.tmp <- subset_samples(ps2rel.tree.sync, source==src)
    bray.tmp    <- phyloseq::distance(ps2.tmp, method="bray")
    disp.test1 <- anova(vegan::betadisper(bray.tmp, sample_data(ps2.tmp)$time))
    permanova1 <- adonis2(bray.tmp ~ sample_data(ps2.tmp)$time, permutations=10000)
    aitchison.tmp <- microViz::dist_calc(ps2.tmp, dist="aitchison")@dist
    disp.test2 <- anova(vegan::betadisper(aitchison.tmp, sample_data(ps2.tmp)$time))
    permanova2 <- adonis2(aitchison.tmp ~ sample_data(ps2.tmp)$time, permutations=10000)
    unifrac.tmp <- microViz::dist_calc(ps2.tmp, dist="wunifrac")@dist
    disp.test3 <- anova(vegan::betadisper(unifrac.tmp, sample_data(ps2.tmp)$time))
    permanova3 <- adonis2(unifrac.tmp ~ sample_data(ps2.tmp)$time, permutations=10000)
    permanova.src.dt <- rbind(permanova.src.dt, data.table(src, dispersion.bray=disp.test1$`Pr(>F)`[1], permanova.bray=permanova1$`Pr(>F)`[1], dispersion.aitchison=disp.test2$`Pr(>F)`[1], permanova.aitchison=permanova2$`Pr(>F)`[1], dispersion.unifrac=disp.test3$`Pr(>F)`[1], permanova.unifrac=permanova3$`Pr(>F)`[1]))
}

permanova.time.dt <- data.table()
for(ti in levels(ps2.tree.sync@sam_data$time)){
    ps2.tmp <- subset_samples(ps2rel.tree.sync, time==ti)
    bray.tmp    <- phyloseq::distance(ps2.tmp, method="bray")
    disp.test1 <- anova(vegan::betadisper(bray.tmp, sample_data(ps2.tmp)$source))
    permanova1 <- adonis2(bray.tmp ~ sample_data(ps2.tmp)$source, permutations=10000)
    aitchison.tmp    <- microViz::dist_calc(ps2.tmp, dist="aitchison")@dist
    disp.test2 <- anova(vegan::betadisper(aitchison.tmp, sample_data(ps2.tmp)$source))
    permanova2 <- adonis2(aitchison.tmp ~ sample_data(ps2.tmp)$source, permutations=10000)
    unifrac.tmp    <- microViz::dist_calc(ps2.tmp, dist="wunifrac")@dist
    disp.test3 <- anova(vegan::betadisper(unifrac.tmp, sample_data(ps2.tmp)$source))
    permanova3 <- adonis2(unifrac.tmp ~ sample_data(ps2.tmp)$source, permutations=10000)
    permanova.time.dt <- rbind(permanova.time.dt, data.table(time=ti, dispersion.bray=disp.test1$`Pr(>F)`[1], permanova.bray=permanova1$`Pr(>F)`[1], dispersion.aitchison=disp.test2$`Pr(>F)`[1], permanova.aitchison=permanova2$`Pr(>F)`[1], dispersion.unifrac=disp.test3$`Pr(>F)`[1], permanova.unifrac=permanova3$`Pr(>F)`[1]))
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
dist.aitchison.melt <- dist_melt(dist.aitchison)
dist.unifrac.melt <- dist_melt(dist.unifrac)
permanova.src.dt[,source1:=src] # dummy column for facet plotting

dist.bray.source.time <-dist.bray.melt[time1==time2 & source1==source2]
ggplot(dist.bray.source.time, aes(x=time1, y=dist)) + geom_boxplot() + facet_wrap(~factor(source1)) + theme_bw(base_size=14) + xlab("Time [h]") + ylab("Beta-diversity (Bray-Curtis)") + geom_text(data=permanova.src.dt, mapping = aes(x = -Inf, y = Inf, label=paste("PERMANOVA, p =",round(permanova.bray,3))), hjust=-0.1, vjust=1.5) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#ggsave("../img/beta-diversity_bray.pdf", width=8, height=2.5)

dist.aitchison.source.time <-dist.aitchison.melt[time1==time2 & source1==source2]
ggplot(dist.aitchison.source.time, aes(x=time1, y=dist)) + geom_boxplot() + facet_wrap(~factor(source1)) + theme_bw(base_size=14) + xlab("Time [h]") + ylab("Beta-diversity (Aitchison)") + geom_text(data=permanova.src.dt, mapping = aes(x = -Inf, y = Inf, label=paste("PERMANOVA, p =",round(permanova.aitchison,3))), hjust=-0.1, vjust=1.5) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#ggsave("../img/beta-diversity_aitchison.pdf", width=8, height=2.5)

dist.unifrac.source.time <-dist.unifrac.melt[time1==time2 & source1==source2]
ggplot(dist.unifrac.source.time, aes(x=time1, y=dist)) + geom_boxplot() + facet_wrap(~source1) + theme_bw(base_size=14) + xlab("Time [h]") + ylab("Beta-diversity (Unifrac)") + geom_text(data=permanova.src.dt, mapping = aes(x = -Inf, y = Inf, label=paste("PERMANOVA, p =",round(permanova.unifrac,3))), hjust=-0.1, vjust=1.5) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#ggsave("../img/beta-diversity_unifrac.pdf", width=8, height=2.5)

# Comparison pairs

dist.bray.pairs <- dist.bray.melt[code1==code2]
dist.bray.asso.alone <- dist.bray.melt[time1==time2 & ((source1=="substrate" & source2=="control") | (source1=="control" & source2=="substrate"))]
dist.bray.host.alone <- dist.bray.melt[time1==time2 & ((source1=="host" & source2=="control") | (source1=="control" & source2=="host"))]

dist.aitchison.pairs <- dist.aitchison.melt[code1==code2]
dist.aitchison.asso.alone <- dist.aitchison.melt[time1==time2 & ((source1=="substrate" & source2=="control") | (source1=="control" & source2=="substrate"))]
dist.aitchison.host.alone <- dist.aitchison.melt[time1==time2 & ((source1=="host" & source2=="control") | (source1=="control" & source2=="host"))]

dist.unifrac.pairs <- dist.unifrac.melt[code1==code2]
dist.unifrac.asso.alone <- dist.unifrac.melt[time1==time2 & ((source1=="substrate" & source2=="control") | (source1=="control" & source2=="substrate"))]
dist.unifrac.host.alone <- dist.unifrac.melt[time1==time2 & ((source1=="host" & source2=="control") | (source1=="control" & source2=="host"))]
#

registerDoParallel(cores)
rand_asso <- function(dist.pairs, dist.asso.alone, dist.host.alone){
	mycombine <- function(x, ...) { mapply(rbind,x,...,SIMPLIFY=FALSE) }
	stat.host.asso  <- adonis2(dist.pairs$dist~dist.pairs$time1, permutation=10000) # pairs are constant for host-substrate
	rand.lst <- foreach(i=1:100, .combine="mycombine", .multicombine=TRUE) %dopar%{
		rand.asso.dt <- data.table()
		for(j in 1:nrow(dist.pairs)){
			time <- dist.pairs[j,time1]
			pair <- dist.pairs[j,code1]
			asso.alone <- dist.asso.alone[time1==time][sample(nrow(dist.asso.alone[time1==time]),1),]
			host.alone <- dist.host.alone[time1==time][sample(nrow(dist.host.alone[time1==time]),1),]
			rand.asso.dt <- rbind(rand.asso.dt, data.table(run=i,time, pair, dist.host.asso=dist.pairs[j,dist], dist.asso.alone=asso.alone$dist, dist.host.alone=host.alone$dist))
		}
		stat.asso.alone <- adonis2(rand.asso.dt$dist.asso.alone~rand.asso.dt$time, permutation=10000)
		stat.host.alone <- adonis2(rand.asso.dt$dist.host.alone~rand.asso.dt$time, permutation=10000)

		list(data.table(run=i, rbind(data.table(cmp="host.substrate",r2=stat.host.asso$R2[1], pval=stat.host.asso$`Pr(>F)`[1]), data.table(cmp="substrate.control",r2=stat.asso.alone$R2[1], pval=stat.asso.alone$`Pr(>F)`[1]), data.table(cmp="host.control",r2=stat.host.alone$R2[1], pval=stat.host.alone$`Pr(>F)`[1]))), rand.asso.dt)
	}
	return(rand.lst)
}
bray.rand.lst <- rand_asso(dist.bray.pairs, dist.bray.asso.alone, dist.bray.host.alone)
aitchison.rand.lst <- rand_asso(dist.aitchison.pairs, dist.aitchison.asso.alone, dist.aitchison.host.alone)
unifrac.rand.lst <- rand_asso(dist.unifrac.pairs, dist.unifrac.asso.alone, dist.unifrac.host.alone)
#saveRDS(bray.rand.lst, "../dat/beta-div_bray.rand-asso.RDS")
#saveRDS(aitchison.rand.lst, "../dat/beta-div_aitchison.rand-asso.RDS")
#saveRDS(unifrac.rand.lst, "../dat/beta-div_unifrac.rand-asso.RDS")

bray.rand.lst[[1]][,list(mean=mean(pval), median=median(pval), sd=sd(pval)),by=cmp]
aitchison.rand.lst[[1]][,list(mean=mean(pval), median=median(pval), sd=sd(pval)),by=cmp]
unifrac.rand.lst[[1]][,list(mean=mean(pval), median=median(pval), sd=sd(pval)),by=cmp]

#bray.rand.lst <- readRDS("../dat/beta-div_bray.rand-asso.RDS")
bray.rand.asso.dt <- bray.rand.lst[[2]]
permanova.bray.rand.asso.dt <- bray.rand.lst[[1]][,list(pval_mean=mean(pval)),by=cmp] 
bray.rand.asso.melt <- melt(bray.rand.asso.dt, id.vars=c("pair","time","run"))
bray.rand.asso.melt[, cmp:=str_remove(gsub("asso","substrate",gsub("alone","control",variable)), "^dist\\.")]
bray.rand.asso.melt <- bray.rand.asso.melt[!(run>1 & cmp=="host.substrate")]
ggplot(bray.rand.asso.melt, aes(x=time, y=value)) + geom_boxplot() + facet_wrap(~cmp) + theme_bw(base_size=14) + ylab("Beta-diversity (Bray-Curtis)") + geom_text(data=permanova.bray.rand.asso.dt, mapping = aes(x = -Inf, y = Inf, label=paste("PERMANOVA, p =",round(pval_mean,3))), hjust=-0.1, vjust=1.5) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#ggsave("../img/beta-diversity_bray-pairs.pdf", width=8, height=2.5)

#aitchison.rand.lst <- readRDS("../dat/beta-div_aitchison.rand-asso.RDS")
aitchison.rand.asso.dt <- aitchison.rand.lst[[2]]
permanova.aitchison.rand.asso.dt <- aitchison.rand.lst[[1]][,list(pval_mean=mean(pval)),by=cmp] 
aitchison.rand.asso.melt <- melt(aitchison.rand.asso.dt, id.vars=c("pair","time","run"))
aitchison.rand.asso.melt[, cmp:=str_remove(gsub("asso","substrate",gsub("alone","control",variable)), "^dist\\.")]
aitchison.rand.asso.melt <- aitchison.rand.asso.melt[!(run>1 & cmp=="host.substrate")]
ggplot(aitchison.rand.asso.melt, aes(x=time, y=value)) + geom_boxplot() + facet_wrap(~cmp) + theme_bw(base_size=14) + ylab("Beta-diversity (Aitchison)") + geom_text(data=permanova.aitchison.rand.asso.dt, mapping = aes(x = -Inf, y = Inf, label=paste("PERMANOVA, p =",round(pval_mean,3))), hjust=-0.1, vjust=1.5) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#ggsave("../img/beta-diversity_aitchison-pairs.pdf", width=8, height=2.5)

#unifrac.rand.lst <- readRDS("../dat/beta-div_unifrac.rand-asso.RDS")
unifrac.rand.asso.dt <- unifrac.rand.lst[[2]]
permanova.unifrac.rand.asso.dt <- unifrac.rand.lst[[1]][,list(pval_mean=mean(pval)),by=cmp] 
unifrac.rand.asso.melt <- melt(unifrac.rand.asso.dt, id.vars=c("pair","time","run"))
unifrac.rand.asso.melt[, cmp:=str_remove(gsub("asso","substrate",gsub("alone","control",variable)), "^dist\\.")]
unifrac.rand.asso.melt <- unifrac.rand.asso.melt[!(run>1 & cmp=="host.substrate")]
ggplot(unifrac.rand.asso.melt, aes(x=time, y=value)) + geom_boxplot() + facet_wrap(~cmp) + theme_bw(base_size=14) + ylab("Beta-diversity (Unifrac)") + geom_text(data=permanova.unifrac.rand.asso.dt, mapping = aes(x = -Inf, y = Inf, label=paste("PERMANOVA, p =",round(pval_mean,3))), hjust=-0.1, vjust=1.5) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#ggsave("../img/beta-diversity_unifrac-pairs.pdf", width=8, height=2.5)
