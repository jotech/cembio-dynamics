library(phyloseq)
library(DESeq2)
library(readxl)
library(ggplot2)
library(ggpubr)
library(decontam)

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
sample.dt[,control:=meta.dat$Comparison_controls[match(sample.dt$sample_name, meta.dat$BARCODE_sample)]]
sample.dt[,neg.control:=control %in% c("NGM_only","Non-template-control")]
sample.df <- sample_data(sample.dt); sample_names(sample.df) <- sample.dt$sample_id
tax.mat <- tax_table(as.matrix(taxa.table))

ps.org <- phyloseq(tax.mat, sample.df, otu_table(asv.table, taxa_are_rows = TRUE))
taxa_names(ps.org) <- paste0("ASV",1:length(taxa_names(ps.org))) # rename to ASV

# remove samples (according to meta data)
ps <- subset_samples(ps.org, inclusion==1 | neg.control | control%in%c("Mock","Zymo_mock"))
ps <- prune_samples(!sample_names(ps) %in% names(which(colSums(asv.table)<12000)), ps) # Remove low read samples (142 samples remain)


df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=neg.control)) + geom_point() + theme_bw(base_size=14)
ggsave("../img/library_sizes.pdf", width=6, height=4)


contamdf.prev <- data.table(isContaminant(ps, method="prevalence", neg="neg.control"), keep.rownames=T)

idx.contam <- match(contamdf.prev[contaminant==T,rn], rownames(ps@tax_table))
ps@tax_table[idx.contam,] # taxonomy of contaminations
rowSums(ps@otu_table[idx.contam,])
