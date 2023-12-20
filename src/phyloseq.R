library(phyloseq)
library(DESeq2)
library(readxl)
library(ggplot2)
library(ggpubr)
library(ape)

asv.table <- fread("../dat/ASV_table.tsv")[,-1]
sample.table <- fread("../dat/sample.sheet.tsv")
taxa.table <- fread("../dat/taxa.own.tbl")
meta.dat <- data.table(read_excel("../dat/CE-CeMbio_V1_revised_20220420.xlsx"))

sample.dt <- sample.table[,list(sample_id=sampleID, sample_name=name)]
sample.dt[,time:=factor(as.numeric(meta.dat$`sampling_time_(h)`[match(sample.dt$sample_name, meta.dat$BARCODE_sample)]))]
sample.dt[,sample_code:=meta.dat$sample_code[match(sample.dt$sample_name, meta.dat$BARCODE_sample)]]
sample.dt[,source:=meta.dat$Association[match(sample.dt$sample_name, meta.dat$BARCODE_sample)]]
sample.dt[,inclusion:=meta.dat$Inclusion_main_analysis[match(sample.dt$sample_name, meta.dat$BARCODE_sample)]]
sample.dt[,condition:=paste0(source,"_",sprintf("%03d",time))]
sample.df <- sample_data(sample.dt); sample_names(sample.df) <- sample.dt$sample_id
tax.mat <- tax_table(as.matrix(taxa.table))

ps.org <- phyloseq(tax.mat, sample.df, otu_table(asv.table, taxa_are_rows = TRUE))
taxa_names(ps.org) <- paste0("ASV",1:length(taxa_names(ps.org))) # rename to ASV

# remove samples (according to meta data)
ps <- subset_samples(ps.org, inclusion==1)

# read count histogram per source
read.count.dt <- data.table(data.frame(id=sample_data(ps)$sample_id, name=sample_data(ps)$sample_name, source=sample_data(ps)$source, time=sample_data(ps)$time, reads=colSums(otu_table(ps))))
g1 <- ggplot(read.count.dt, aes(x=reads, fill=source)) + geom_histogram() + theme_bw(base_size=14) + ggtitle("All samples")
g2 <- ggplot(read.count.dt[reads<=1000], aes(x=reads, fill=source)) + geom_histogram() + theme_bw(base_size=14)  + ggtitle("Samples with reads <= 1000")
ggpubr::ggarrange(g1, g2, labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave("../img/sample-reads_hist.pdf", width=8, height=3)

# read count per source and time (cap for reads > 20000)
read.count.dt[,reads2:=ifelse(reads>20000,20000,reads)]
ggplot(read.count.dt[!time%in%c(0,2)], aes(x=reads2, color=source)) + geom_freqpoly(breaks=c(seq(0, 20000, by=1000))) + theme_bw(base_size=14) + facet_wrap(~time, ncol=2,nrow=4) #+ scale_x_continuous(limits=c(0, 20000), breaks=c(seq(0, 20000, by=5000)), labels=c(seq(0,19000, by=5000), "20000+"))
ggsave("../img/sample-reads_hist2.pdf", width=8, height=6)

# check for bias in sample count per source and time
compare_means(reads~source, data=read.count.dt, group.by="time")
compare_means(reads~source, data=read.count.dt[reads>=12000], group.by="time")

read.count.dt[,source:=factor(source)]
levels(read.count.dt$source) <- c("control","substrate", "host") # rename sample sources
ggplot(read.count.dt, aes(x=reads2, color=source)) + geom_freqpoly(breaks=c(seq(0, 20000, by=1000))) + theme_bw(base_size=14) + xlab("Number of reads per sample") + ylab("Count")
ggsave("../img/sample-reads_hist3.pdf", width=5, height=3)

#############
# filtering #
#############

# 1) read count
ps <- prune_samples(!sample_names(ps) %in% names(which(colSums(asv.table)<12000)), ps) # Remove low read samples (142 samples remain)

# 2) prevalence
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) # remove uncharacterized 
prev = unname(apply(X = otu_table(ps0), 1, FUN = function(x){sum(x > 0)}))
prev.dt = data.table(Prevalence = prev,TotalAbundance = unname(taxa_sums(ps0)),data.frame(tax_table(ps0),row.names=NULL))
prev.phylum.dt <- prev.dt[,list(Prevalence.mean=mean(as.double(Prevalence)), Prevalence.sum=sum(Prevalence), Abundance=round(sum(TotalAbundance)/sum(prev.dt$TotalAbundance),4)),by=Phylum]
#
prev.phylum.dt[order(Abundance, decreasing=T)]
ps1 = subset_taxa(ps0, !Phylum %in% prev.phylum.dt[Prevalence.sum<5,Phylum]) # cyanobacteria for prevalence cutoff


#################
# agglomerating #
#################
ps2 <- ps1
for(tax in unique(unname(ps1@tax_table[,8]))){
  idx <- which(ps2@tax_table[,8]==as.character(tax))
  if(length(idx) >= 2) {
    ps2 <- merge_taxa(ps2, idx, 2) # sums up reads
    if(any(is.na(ps2@tax_table[,8]))) ps2@tax_table[,8][which(is.na(ps2@tax_table[,8]))] <- as.character(tax) # group entry can be missing for merged taxa
  }
}
ps2@tax_table[,8][which(ps2@tax_table[,8]=="MYb371,MYb331,MYb330,MYb177")] <- "MYb371,MYb331,MYb330" # MYb177 is Acinetobacter not Pseudomonas (new sequencing 07/2023)
ps2 <- prune_taxa(! taxa_names(ps2) %in% taxa_names(ps2)[ps2@tax_table[,8] == "JUb66,MYb186"],ps2) # remove ambigious entry
ps2 <- filter_taxa(ps2, function(x) mean(x) > 0, TRUE) # remove low abundant
# check mean reads per species
ps2.tax.dt <- data.table(id=as.character(ps2@tax_table[,8]), mean.reads=colMeans(ps2@otu_table))


##########
# export #
##########

ps2@sam_data$condition <- paste0(ps2@sam_data$source, "_", ps2@sam_data$time)
sample_data(ps2)$source <- factor(sample_data(ps2)$source)

# fix name for merged species
ps2@tax_table[,7][which(ps2@tax_table[,8]=="MYb71,MYb49")] <- "Ochrobactrum sp."
ps2@tax_table[,7][which(ps2@tax_table[,8]=="MYb176,MYb174")] <- "Enterobacter sp."
ps2@tax_table[,7][which(ps2@tax_table[,8]=="MYb191,MYb177")] <- "Acinetobacter sp."
ps2@tax_table[,7][which(ps2@tax_table[,8]=="MYb371,MYb331,MYb330")] <- "Pseudomonas sp."
# Ochrobactrums had wrong order, family, genus classified by DADA2 taxonomy
ps2@tax_table[,4][which(ps2@tax_table[,8]=="MYb71,MYb49")] <- as.character(ps2@tax_table[,4][which(ps2@tax_table[,8]=="MYb58")]) 
ps2@tax_table[,5][which(ps2@tax_table[,8]=="MYb71,MYb49")] <- as.character(ps2@tax_table[,5][which(ps2@tax_table[,8]=="MYb58")]) 
ps2@tax_table[,6][which(ps2@tax_table[,8]=="MYb71,MYb49")] <- as.character(ps2@tax_table[,6][which(ps2@tax_table[,8]=="MYb58")]) 
# missing genus for Comamonas
ps2@tax_table[,6][which(ps2@tax_table[,8]=="MYb396,MYb69,MYb21")] <- "Comamonas"
ps2@tax_table[,7][which(ps2@tax_table[,8]=="MYb396,MYb69,MYb21")] <- "Comamonas sp."
# MYb174,176 are Enterobacter, wrong Genus classified by DADA2 taxonomy
ps2@tax_table[,6][which(ps2@tax_table[,8]=="MYb176,MYb174")] <- "Enterobacter"
# set taxa names to ids
taxa_names(ps2) <- ps2@tax_table[,8]

# export phyloseq object
ps2.old <- readRDS("../dat/phyloseq_ps2.RDS")
if(!is.logical(all.equal(ps2,ps2.old))){
    warning("phyloseq object mismatch")
    saveRDS(ps2, "../dat/phyloseq_ps2.RDS")
}

asv.old <- read.table("../dat/ASV_table-phyloseq.tsv")
if(!is.logical(all.equal(asv.old,as.data.frame(ps2@otu_table)))){
    warning("filtered ASV table mismatch")
    write.table(ps2@otu_table, "../dat/ASV_table-phyloseq.tsv", sep="\t")
}

ps2.rel <- transform_sample_counts(ps2, function(x){x / sum(x)})
asv.rel.old <- read.table("../dat/ASVrel_table-phyloseq.tsv")
if(!is.logical(all.equal(asv.rel.old,as.data.frame(ps2.rel@otu_table)))){
    warning("filtered relative ASV table mismatch")
    write.table(ps2.rel@otu_table, "../dat/ASVrel_table-phyloseq.tsv", sep="\t")
}

# tree

tree <- read.tree("../dat/gtdbtk.bac120.user_msa.fasta.gz_aggregated-rooted.treefile")
tree$tip.label <- gsub("\\.",",", tree$tip.label)
tree$tip.label <- gsub("\\'","", tree$tip.label)
ps2a <- subset_taxa(ps2, id!="unknown") # remove unknown group
ps2.tree <- phyloseq(otu_table(ps2a), tax_table(ps2a), sample_data(ps2a), phy_tree(tree))
ps2.tree.old <- readRDS("../dat/phyloseq_ps2.tree.RDS")
if(!is.logical(all.equal(ps2.tree,ps2.tree.old))){
    warning("phyloseq object mismatch")
    saveRDS(ps2.tree, "../dat/phyloseq_ps2.tree.RDS")
}
