library(NST)
library(phyloseq)
library(ggpubr)
library(foreach)
library(doParallel)

cores <- 3
rand.time <- 1000

ps2 <- readRDS("../dat/phyloseq_ps2.RDS")
ps2 <- subset_taxa(ps2, id!="unknown")
ps2 <- subset_samples(ps2, !time %in% c(0,2)) # sync time

comm <- data.frame(t(otu_table(ps2)))
treat <- data.frame(sample_data(ps2))[,c("time","source"), drop=F]
treat$source <- as.character(treat$source)

tax.ck=NST::match.name(cn.list = list(comm=comm))
comm=tax.ck$comm
all(rownames(comm) == rownames(treat))

treat.use=treat[,2,drop=F] # source
#tnst=NST::tNST(comm=comm, group=treat.use, meta.group=treat.use, dist.method="bray", abundance.weighted=TRUE, rand=rand.time, nworker=cores, null.model="PF", output.rand = TRUE)
tnst=NST::tNST(comm=comm, group=treat.use, rand=rand.time, nworker=cores, output.rand = TRUE)
tnst$index.grp
tnst.bt=NST::nst.boot(nst.result=tnst, rand=rand.time, nworker=cores, out.detail=TRUE, between.group=TRUE)
data.table(tnst.bt$compare)[Index=="NST",.(Index,group1,group2,group1.obs,group2.obs,p.count)]
#save(pnst, pnst.bt, tnst, tnst.bt, file="../dat/NST.RData", compress="xz")

tnst.bt$summary[,1:8]
tnst.dt <- tnst.bt$summary[tnst.bt$summary$Index=="NST" & tnst.bt$summary$Group %in% c("host", "associated", "alone"),c("Group","mean","stdev")]
stat.dt <- tnst.bt$compare[tnst.bt$compare$Index=="NST" ,c("group1","group2","p.count")]
stat.dt$y.position <- c(.72,.67,.62)
stat.dt$p.value <- round(as.numeric(stat.dt$p.count),3); stat.dt$p.signif <- ifelse(stat.dt$p.value>=0.05, "ns", "*")
ggplot(tnst.dt) + geom_point(aes(x=Group, y=mean), size=1.5) + geom_errorbar(aes(x=Group, ymin=ifelse(mean-stdev<0,0,mean-stdev), ymax=mean+stdev),width=0.2) + theme_bw(base_size=14) + stat_pvalue_manual(stat.dt, label = "p.signif", hide.ns=T) + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + xlab("") + ylab("Stochasticity [tNST]") + geom_hline(yintercept=0.5, linetype="dashed", color = "red")
ggsave("../img/nst-tnst.pdf", width=3, height=2.5)


#
# per source
#
nst.s.boot.dt <- data.table()
nst.s.pval.dt <- data.table()
for(i in c("host", "associated", "alone")){
	ps2.s   <- subset_samples(ps2, source==i & time %in% c(16, 42, 66, 90, 138, 186))
	comm.s  <- data.frame(t(otu_table(ps2.s)))
    treat.s <- data.frame(sample_data(ps2.s))[,c("time"), drop=F]; treat.s$time <- as.numeric(as.character(treat.s$time))
    #tnst.s=NST::tNST(comm=comm.s, group=treat.s, meta.group=treat.s, dist.method="bray", abundance.weighted=TRUE, rand=rand.time, nworker=cores, null.model="PF", output.rand = TRUE)
    tnst.s=NST::tNST(comm=comm.s, group=treat.s, rand=rand.time, nworker=cores, output.rand = TRUE)
    tnst.s.boot=NST::nst.boot(nst.result=tnst.s, rand=rand.time, nworker=cores)
    nst.s.boot.dt <- rbind(nst.s.boot.dt, data.table(source=i, algo="tnst", tnst.s.boot$summary[,1:8]))
    nst.s.pval.dt <- rbind(nst.s.pval.dt, data.table(source=i, algo="tnst", tnst.s.boot$compare))
}
#save(nst.s.boot.dt,nst.s.pval.dt, file="~/uni/cembio.ext/dat/NST_time.RData")
#load("~/uni/cembio.ext/dat/NST_time.RData")
stat.dt <- nst.s.pval.dt[Index=="NST" ,c("source","algo","group1","group2","p.count")]
stat.dt[,group1:=factor(group1)]; stat.dt[,group2:=factor(group2)]
stat.dt$p.value <- round(as.numeric(stat.dt$p.count),3); stat.dt$p.signif <- ifelse(stat.dt$p.value>=0.05, "ns", "*")
stat.dt <- stat.dt[p.signif!="ns"]
y.pos.N <- max(stat.dt[,.N,by=.(source,algo)]$N)
stat.dt[algo=="tnst", y.position:=rep(seq(1,1.5,by=0.5/y.pos.N),3)[1:.N]]
ggplot(nst.s.boot.dt[algo=="tnst" & Index=="NST"]) + geom_point(aes(x=factor(as.numeric(Group)), y=mean), size=1.5) + geom_line(aes(x=factor(as.numeric(Group)), y=mean, group=source)) + geom_errorbar(aes(x=factor(as.numeric(Group)), ymin=ifelse(mean-stdev<0,0,mean-stdev), ymax=mean+stdev), alpha=0.5, width=0.2) + stat_pvalue_manual(stat.dt[algo=="tnst"], label = "p.signif", hide.ns=T) + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + facet_wrap(~source) + theme_bw(base_size=14) + ylab("Stochasticity [tNST]") + xlab("Time [h]") + geom_hline(yintercept=0.5, linetype="dashed", color = "red")
ggsave("../img/nst-tnst_source.pdf", width=9, height=2.5)


#
# per time
#
nst.t.boot.dt <- data.table()
nst.t.pval.dt <- data.table()
for(i in c(16, 42, 66, 90, 138, 186)){
	ps2.t   <- subset_samples(ps2, time==i)
	comm.t  <- data.frame(t(otu_table(ps2.t)))
    treat.t <- data.frame(sample_data(ps2.t))[,c("source"), drop=F]; treat.t$source <- as.character(treat.t$source)
    #tnst.t=NST::tNST(comm=comm.t, group=treat.t, meta.group=treat.t, dist.method="bray", abundance.weighted=TRUE, rand=rand.time, nworker=cores, null.model="PF", output.rand = TRUE)
    tnst.t=NST::tNST(comm=comm.t, group=treat.t, rand=rand.time, nworker=cores, output.rand = TRUE)
    tnst.t.boot=NST::nst.boot(nst.result=tnst.t, rand=rand.time, nworker=cores)
    nst.t.boot.dt <- rbind(nst.t.boot.dt, data.table(t=i, algo="tnst", tnst.t.boot$summary[,1:8]))
    nst.t.pval.dt <- rbind(nst.t.pval.dt, data.table(t=i, algo="tnst", tnst.t.boot$compare))
}
#save(nst.t.boot.dt,nst.t.pval.dt, file="~/uni/cembio.ext/dat/NST_source.RData", compress="xz")
stat.dt <- nst.t.pval.dt[Index=="NST" ,c("t","algo","group1","group2","p.count")]
stat.dt$p.value <- round(as.numeric(stat.dt$p.count),3); stat.dt$p.signif <- ifelse(stat.dt$p.value>=0.05, "ns", "*")
stat.dt <- stat.dt[p.signif!="ns"]
stat.dt[algo=="tnst", y.position:=c(1,1.1,1.2)[1:.N]]
stat.dt[,time:=factor(t)]; nst.t.boot.dt[,time:=factor(t)]
ggplot(nst.t.boot.dt[algo=="tnst" & Index=="NST"]) + geom_point(aes(x=Group, y=mean), size=1.5) + geom_errorbar(aes(x=Group, ymin=ifelse(mean-stdev<0,0,mean-stdev), ymax=ifelse(mean+stdev>1,1,mean+stdev)), width=0.2) + stat_pvalue_manual(stat.dt[algo=="tnst"], label = "p.signif", hide.ns=T) + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + facet_wrap(~time, labeller=label_both) + theme_bw(base_size=14) + ylab("Stochasticity [tNST]") + geom_hline(yintercept=0.5, linetype="dashed", color = "red")+ xlab("")
ggsave("../img/nst-tnst_time.pdf", width=7, height=5)


