library(DESeq2)
library(phyloseq)
library(ggplot2)
library(data.table)
library(stringr)

ps2 <- readRDS("../dat/phyloseq_ps2.RDS")
ps2 <- subset_samples(ps2, !time %in% c(0,2)) # sync time
levels(ps2@sam_data$source) <- c("control", "substrate", "host") # rename sample sources

ps2.asso.worm  <- subset_samples(ps2, source %in% c("host","substrate"))
ps2.alone.asso <- subset_samples(ps2, source %in% c("control","substrate")) # 2 more time points in cembio.alone
ps2.alone.worm <- subset_samples(ps2, source %in% c("control","host")) # 2 more time points in cembio.alone


meta.pwy <- fread("../dat/meta_pwy.tbl")
custom.pwy <- fread("../dat/custom_pwy.tbl")
pwy.desc <- rbind(meta.pwy, custom.pwy, fill=T)
abricate.desc.dt <- fread("../dat/abricate-description.csv"); abricate.desc.dt[,GENE:=gsub("'","",GENE)]
gutsmash.desc.dt <- fread("../dat/gutsmash-description.csv")
dbcan.desc.dt <- fread("../dat/CAZyDB.07302020.fam-activities.txt", skip=
1, fill=T, sep="\t", col.names=c("id","desc","V3"))
seed.dt <- fread("../dat/seed_metabolites_edited.tsv")

pwy.abundances.count <- read.table("../dat/pwy-table_count.tbl")
abricate.abundances.count <- read.table("../dat/abricate-table_count.tbl")
gutsmash.abundances.count <- read.table("../dat/gutsmash-table_count.tbl")
dbcan.abundances.count <- read.table("../dat/dbcan-table_count.tbl")	
uast.abundances.count <- read.table("../dat/uast-table_count.tbl")	
medium.abundances.count <- read.table("../dat/medium-table_count.tbl")	
csferm.abundances.count <- read.table("../dat/fermcs-table_count.tbl")	
interact.abundances.count <- read.table("../dat/interact-table_count.tbl")

feat.abundances.count <- rbind(pwy.abundances.count, abricate.abundances.count, gutsmash.abundances.count, dbcan.abundances.count, uast.abundances.count, medium.abundances.count, csferm.abundances.count, interact.abundances.count)
dim(feat.abundances.count)
feat.abundances.count <- feat.abundances.count[which(rowSums(feat.abundances.count)!=0),] # remove zero count feat
dim(feat.abundances.count)

ps2.feat   <- phyloseq(sample_data(ps2), otu_table(feat.abundances.count, taxa_are_rows = TRUE))
ps2.feat.old <- readRDS("../dat/phyloseq_ps2-feat.RDS")
if(!is.logical(all.equal(ps2.feat,ps2.feat.old))){
    warning("phyloseq object mismatch")
    saveRDS(ps2.feat, "../dat/phyloseq_ps2-feat.RDS")
}
ps2.asso.worm.feat   <- phyloseq(sample_data(ps2.asso.worm), otu_table(feat.abundances.count, taxa_are_rows = TRUE))
ps2.alone.asso.feat   <- phyloseq(sample_data(ps2.alone.asso), otu_table(feat.abundances.count, taxa_are_rows = TRUE))
ps2.alone.worm.feat   <- phyloseq(sample_data(ps2.alone.worm), otu_table(feat.abundances.count, taxa_are_rows = TRUE))

ddseq2_analysis <- function(ps.pwy){
	ps.pwy.dds <- phyloseq_to_deseq2(ps.pwy, ~ time + source)
	gm_mean = function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
	ps.pwy.dds.geoMeans = apply(counts(ps.pwy.dds), 1, gm_mean)
	ps.pwy.dds <- estimateSizeFactors(ps.pwy.dds, geoMeans=ps.pwy.dds.geoMeans)
	ps.pwy.dds <- estimateDispersions(ps.pwy.dds)
	ps.pwy.dds <- DESeq(ps.pwy.dds, fitType="local")
	ps.pwy.dds_res <- as.data.table(results(ps.pwy.dds)); ps.pwy.dds_res[,id:=taxa_names(ps.pwy)]
	ps.pwy.dds_res[,name:=pwy.desc$name[match(id, pwy.desc$id)]]
	ps.pwy.dds_res[,hierarchy:=pwy.desc$hierarchy[match(id, pwy.desc$id)]]
	ps.pwy.dds_res[,subsystem:=str_extract(hierarchy, ("(Degradation|Energy-Metabolism|Biosynthesis|Detoxification|Glycan-Pathways|Macromolecule-Modification|Metabolic-Clusters|Signaling-Pathways|Transport-Pathways|Transport|Activation-Inactivation-Interconversion|Bioluminescence|Enzyme-Test)"))]
	ps.pwy.dds_res[is.na(name), name:=abricate.desc.dt$PRODUCT[match(id, abricate.desc.dt$GENE)]]
	ps.pwy.dds_res[is.na(name), name:=gutsmash.desc.dt$`Most similar known cluster`[match(id, gutsmash.desc.dt$Type)]]
	ps.pwy.dds_res[is.na(name), name:=dbcan.desc.dt$desc[match(id, dbcan.desc.dt$id)]]
	ps.pwy.dds_res[is.na(name), name:=seed.dt$name[match(str_extract(id,"cpd[0-9]+"),seed.dt$id)]]
	ps.pwy.dds_res[is.na(name) | name=="" | length(name)==0, name:=id]
	ps.pwy.dds_res[is.na(subsystem), subsystem:=ifelse(id %in% gsub("'","",rownames(abricate.abundances.count)), "virulence", ifelse(id %in% rownames(gutsmash.abundances.count), "gut", ifelse(id %in% rownames(dbcan.abundances.count), "cazyme", ifelse(id %in% rownames(uast.abundances.count), "uast", ifelse(id %in% rownames(medium.abundances.count), "medium", ifelse(id %in% rownames(csferm.abundances.count), "exchange", ifelse(id %in% rownames(interact.abundances.count),"interactions",NA)))))))]
	return(ps.pwy.dds_res)
}

ps2.asso.worm.feat.dds <- ddseq2_analysis(ps2.asso.worm.feat)
ps2.alone.asso.feat.dds <- ddseq2_analysis(ps2.alone.asso.feat)
ps2.alone.worm.feat.dds <- ddseq2_analysis(ps2.alone.worm.feat)

# diversity

# subsystem comparison
ps2.combined.feat.dds <- rbind(data.table(cmp="substrate vs. host", ps2.asso.worm.feat.dds[padj<=0.05]),
	  data.table(cmp="control vs. substrate", ps2.alone.asso.feat.dds[padj<=0.05]),
	  data.table(cmp="control vs. host",ps2.alone.worm.feat.dds[padj<=0.05]))
ps2.combined.feat.dds[!is.na(hierarchy), subsystem:="metabolism"]
ps2.combined.feat.dds[,cmp:=factor(cmp, levels=c("control vs. substrate","control vs. host","substrate vs. host"))]

#saveRDS(ps2.combined.feat.dds, "../dat/diff-functions.RDS")
# ps2.combined.feat.dds <- readRDS("../dat/diff-functions.RDS")

# plotting

ggplot(ps2.combined.feat.dds, aes(x=subsystem, y=log2FoldChange)) + geom_boxplot() + coord_flip() + theme_minimal(base_size=14) + ylab("log2 fold change") + xlab("Subsystem") + geom_hline(yintercept=0, linetype="dashed", color = "red") + scale_x_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_subsystem.pdf", height=2.5, width=7)

ggplot(ps2.combined.feat.dds[subsystem=="uast"], aes(y=name, x=log2FoldChange)) +  geom_segment(aes(yend=name), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("Strategies (UAST)") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_uast.pdf", height=2, width=7.5)

ps2.combined.feat.dds[,id2:=gsub("_"," ",gsub("2"," to ",id))]
ggplot(ps2.combined.feat.dds[subsystem=="gut"], aes(y=id2, x=log2FoldChange)) +  geom_segment(aes(yend=id2), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("Gut-related gene cluster") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_gut.pdf", height=6, width=9.5)

ggplot(ps2.combined.feat.dds[subsystem=="medium"], aes(y=name, x=log2FoldChange)) +  geom_segment(aes(yend=name), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("Growth medium") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_medium.pdf", height=2.5, width=6.7)

ps2.combined.feat.dds[subsystem=="exchange", name:=paste0(name, " (",ifelse(str_extract(id,"(cs|ferm)")=="ferm","pro","up"),")")]
ggplot(ps2.combined.feat.dds[subsystem=="exchange"], aes(y=name, x=log2FoldChange)) +  geom_segment(aes(yend=name), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("Uptake and production") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_exchanges.pdf", height=7, width=8)

ggplot(ps2.combined.feat.dds[subsystem=="virulence"], aes(y=id, x=log2FoldChange)) +  geom_segment(aes(yend=id), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("Virulence genes") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_virulence.pdf", height=20, width=12)

ggplot(ps2.combined.feat.dds[subsystem=="cazyme"], aes(y=str_trunc(name, 55, "right"), x=log2FoldChange)) +  geom_segment(aes(yend=str_trunc(name, 55, "right")), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("carbohydrate-activae enzymes (cazymes)") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_cazyme.pdf", height=15, width=12)

ggplot(ps2.combined.feat.dds[subsystem=="interactions"], aes(y=str_trunc(name, 55, "right"), x=log2FoldChange)) +  geom_segment(aes(yend=str_trunc(name, 55, "right")), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("interactions") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_interactions.pdf", height=1.7, width=7)

ps2.combined.feat.dds[subsystem=="metabolism", subsystem2:=str_extract(hierarchy, "Activation-Inactivation-Interconversion|Bioluminescence|Biosynthesis|Degradation|Detoxification|Energy-Metabolism|Glycan-Pathways|Macromolecule-Modification|Metabolic-Clusters|Signaling-Pathways|Transport-Pathways")]
ps2.combined.feat.dds[subsystem=="metabolism" & is.na(subsystem2), subsystem2:="Other"]
ggplot(ps2.combined.feat.dds[subsystem=="metabolism"], aes(y=str_trunc(subsystem2, 55, "right"), x=log2FoldChange)) +  geom_segment(aes(yend=str_trunc(subsystem2, 55, "right")), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("MetaCyc subsystems") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_metabolism.pdf", height=3.5, width=9)


# dbcan substrates
dbcan.long <- fread("~/uni/cembio.ext/src/agnes/dat/dbcan_substrates.csv")
dbcan.long[,id2:=str_remove(id,"_[0-9]+$")]
 ps2.combined.feat.dds[subsystem=="cazyme", substrate:=paste0(unique(dbcan.long$substrate[grep(id, dbcan.long$id)]), collapse=","), by=id]

dbcan.sub.dds <- merge(ps2.combined.feat.dds[subsystem=="cazyme"], dbcan.long[,-"id"], by.x="id", by.y="id2")

ggplot(dbcan.sub.dds, aes(x=substrate, y=log2FoldChange)) + geom_boxplot() + coord_flip() + theme_minimal(base_size=14) + ylab("log2 fold change") + xlab("Cazymes predicted substrate") + geom_hline(yintercept=0, linetype="dashed", color = "red") + scale_x_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_cazyme-substrates.pdf", height=12, width=12)


# selecting features with FC in same direction for host vs. control/substrate
feat.2cmp <- ps2.combined.feat.dds[cmp%in%c("control vs. host", "substrate vs. host"),.N,by=id][N==2,id]
ps2.2cmp.sign <- ps2.combined.feat.dds[id %in%feat.2cmp, prod(log2FoldChange)>0, by=id]$id
ps2.combined.feat.dds[id%in%ps2.2cmp.sign,-c(hierarchy)][order(id)]
ps2.combined.feat.dds[subsystem=="cazyme", substrate:=paste0(unique(dbcan.long$substrate[grep(id, dbcan.long$id)]), collapse=","), by=id]

# PCA
library(sparsepca)
spca.dat <- spca(t(ps2.feat@otu_table), k=3)
spca.loadings.dt <- data.table(id=rownames(ps2.feat@otu_table), spca.dat$loadings)
spca.dt <- data.table(sample=rownames(spca.dat$scores),spca.dat$scores))
spca.dt[,time:=ps2@sam_data$time[match(sample,ps2@sam_data$sample_id)]]
spca.dt[,source:=ps2@sam_data$source[match(sample,ps2@sam_data$sample_id)]]
ggplot(spca.dt, aes(x=V1,y=V2)) + geom_point(size=3, aes(color=source)) + theme_minimal(base_size=14)

dist.aitchison <- microViz::dist_calc(ps2.feat, dist="aitchison")@dist
ord.mds.aitchison <- ordinate(ps2.feat, method="PCoA", distance=dist.aitchison)
plot_ordination(ps2.feat, ord.mds.aitchison, type="samples", color="time") + geom_point(size=3, shape=21, color="black", aes(fill=time)) + facet_wrap(~source) + scale_fill_brewer(type="qual", palette="Oranges") + theme_bw(base_size=14)
#ggsave("../img/ordination-pcoa_aitchison.pdf", width=8, height=2.5)
