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
#saveRDS(list(ps2.asso.worm.feat.dds, ps2.alone.asso.feat.dds, ps2.alone.worm.feat.dds), "../dat/diff-functions_all-cmp.RDS")

# subsystem comparison
if(nrow(ps2.alone.asso.feat.dds[padj<=0.05])==0) print("No differential function between control and substrate")
ps2.combined.feat.dds <- rbind(data.table(cmp="substrate vs. host", ps2.asso.worm.feat.dds[padj<=0.05]),
	  #data.table(cmp="control vs. substrate", ps2.alone.asso.feat.dds[padj<=0.05]),
	  data.table(cmp="control vs. host",ps2.alone.worm.feat.dds[padj<=0.05]))
ps2.combined.feat.dds[!is.na(hierarchy), subsystem:="metabolism"]
ps2.combined.feat.dds[,cmp:=factor(cmp, levels=c("control vs. substrate","control vs. host","substrate vs. host"))]

# saveRDS(ps2.combined.feat.dds, "../dat/diff-functions.RDS")
# fwrite(ps2.combined.feat.dds, "../dat/diff-functions.csv")
# ps2.combined.feat.dds <- readRDS("../dat/diff-functions.RDS")

# plotting

ggplot(ps2.combined.feat.dds, aes(x=subsystem, y=log2FoldChange)) + geom_boxplot() + coord_flip() + theme_minimal(base_size=14) + ylab("log2 fold change") + xlab("Subsystem") + geom_hline(yintercept=0, linetype="dashed", color = "red") + scale_x_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_subsystem.pdf", height=2.5, width=5)

ggplot(ps2.combined.feat.dds[subsystem=="uast"], aes(y=name, x=log2FoldChange)) +  geom_segment(aes(yend=name), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("Strategies (UAST)") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_uast.pdf", height=2.5, width=7)

ps2.combined.feat.dds[,id2:=gsub(",",", ",gsub("_"," ",gsub("2"," to ",id)))]
substr(ps2.combined.feat.dds$id2, 1, 1) <- toupper(substr(ps2.combined.feat.dds$id2, 1, 1))
ggplot(ps2.combined.feat.dds[subsystem=="gut"], aes(y=id2, x=log2FoldChange)) +  geom_segment(aes(yend=id2), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("Gut-related gene cluster") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_gut.pdf", height=6, width=9.7)

ggplot(ps2.combined.feat.dds[subsystem=="medium"], aes(y=name, x=log2FoldChange)) +  geom_segment(aes(yend=name), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("Growth medium") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_medium.pdf", height=4.5, width=6.7)

ps2.combined.feat.dds[subsystem=="exchange", name:=paste0(name, " (",ifelse(str_extract(id,"(cs|ferm)")=="ferm","pro","up"),")")]
ggplot(ps2.combined.feat.dds[subsystem=="exchange"], aes(y=name, x=log2FoldChange)) +  geom_segment(aes(yend=name), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("Uptake and production") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_exchanges.pdf", height=7, width=6.5)

vfdb.categories <- c("Adherence", "Antimicrobial activity/Competitive advantage", "Biofilm", "Effector delivery system", "Exotoxin", "Exoenzyme", "Immune modulation", "Invasion", "Motility", "Nutritional/Metabolic factor", "Regulation", "Stress survival", "Post-translational modification", "Others")
ggplot(ps2.combined.feat.dds[subsystem=="virulence" & id %in% vfdb.categories], aes(y=id, x=log2FoldChange)) +  geom_segment(aes(yend=id), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("Virulence genes") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_virulence.pdf", height=3, width=9)

ggplot(ps2.combined.feat.dds[subsystem=="cazyme"], aes(y=str_trunc(name, 55, "right"), x=log2FoldChange)) +  geom_segment(aes(yend=str_trunc(name, 55, "right")), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("carbohydrate-activae enzymes (cazymes)") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_cazyme.pdf", height=15, width=10)

ggplot(ps2.combined.feat.dds[subsystem=="interactions"], aes(y=str_trunc(name, 55, "right"), x=log2FoldChange)) +  geom_segment(aes(yend=str_trunc(name, 55, "right")), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("interactions") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_interactions.pdf", height=2, width=7)

ps2.combined.feat.dds[subsystem=="metabolism", subsystem2:=str_extract(hierarchy, "Activation-Inactivation-Interconversion|Bioluminescence|Biosynthesis|Degradation|Detoxification|Energy-Metabolism|Glycan-Pathways|Macromolecule-Modification|Metabolic-Clusters|Signaling-Pathways|Transport-Pathways")]
ps2.combined.feat.dds[subsystem=="metabolism" & is.na(subsystem2), subsystem2:="Other"]
ggplot(ps2.combined.feat.dds[subsystem=="metabolism"], aes(y=str_trunc(subsystem2, 55, "right"), x=log2FoldChange)) +  geom_segment(aes(yend=str_trunc(subsystem2, 55, "right")), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("MetaCyc subsystems") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_metabolism.pdf", height=3.5, width=7)


# hydroxy proline/acetoin
ggplot(ps2.combined.feat.dds[grepl("acetoin", name) | (grepl("proline", name) & grepl("hydroxy",name))], aes(y=str_trunc(name, 55, "right"), x=log2FoldChange)) +  geom_segment(aes(yend=str_trunc(name, 55, "right")), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_hydroxyproline-acetoin.pdf", height=2.5, width=9)

# BCAA
ps2.combined.feat.dds[grepl("leucine|valine", name),-"hierarchy"]
ggplot(ps2.combined.feat.dds[grepl("leucine|valine", name)], aes(y=str_trunc(name, 55, "right"), x=log2FoldChange)) +  geom_segment(aes(yend=str_trunc(name, 55, "right")), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_bcaa.pdf", height=2, width=5)

# known microbial metabolites influencing CE
ggplot(ps2.combined.feat.dds[grepl("\\bethanol\\b|\\bindole\\b|dopamine|spermidin|quinolin|propionate|\\bacetate\\b", name)], aes(y=str_trunc(name, 55, "right"), x=log2FoldChange)) +  geom_segment(aes(yend=str_trunc(name, 55, "right")), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_knownmicmet.pdf", height=6, width=9)

# b12
ggplot(ps2.combined.feat.dds[grepl("cobalamin|b12|cobalt", name,ignore.case=T)], aes(y=str_trunc(name, 55, "right"), x=log2FoldChange)) +  geom_segment(aes(yend=str_trunc(name, 55, "right")), xend=0, colour="grey50") + geom_point(size=4, aes(color=baseMean)) + theme_minimal(base_size=14) + xlab("log2 fold change") + ylab("") + geom_vline(xintercept=0, linetype="dashed", color = "red") + scale_y_discrete(limits=rev) + facet_wrap(~cmp)
ggsave("../img/diff-func_b12.pdf", height=2, width=7)

