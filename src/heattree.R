library(DESeq2)
library(metacoder)
library(phyloseq)
library(ggplot2)
#devtools::install_github("grunwaldlab/metacoder")

ps2 <- readRDS("../dat/phyloseq_ps2.RDS")
ps2 <- subset_taxa(ps2, !id %in% c("unknown","OP50")) # remove non-cembio members
ps2@tax_table[,7] <- ps2@tax_table[,8]; ps2@tax_table <- ps2@tax_table[,-8] # remove id column for accumulated abundance calculation (otherwise species will be treated twice)
ps2.taxmap <- parse_phyloseq(ps2)
ps2.taxmap$data$tax_abund <- calc_taxon_abund(ps2.taxmap, "otu_table")

deseq_diffspec_taxmap <- function(ps2.tmp){
  ps2.tmp.taxmap <- parse_phyloseq(ps2.tmp)
  ps2.tmp.taxmap$data$tax_abund <- calc_taxon_abund(ps2.tmp.taxmap, "otu_table")
  ps2.tmp.tax <- as_phyloseq(ps2.tmp.taxmap, otu_table="tax_abund", otu_id_col="taxon_id")
	#ps2.tmp.tax <- as_phyloseq(ps2.tmp.taxmap, otu_table="otu_table", otu_id_col="otu_id")
  ps2.tmp.tax@tax_table[,7] <- apply(ps2.tmp.tax@tax_table,1,function(x)na.omit(x)[length(na.omit(x))]) # use non-na taxon as id
	
	ps2.tmpdds <- phyloseq_to_deseq2(ps2.tmp.tax, ~ time + source)
  gm_mean = function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
  geoMeans = apply(counts(ps2.tmpdds), 1, gm_mean)
  ps2.tmpdds <- estimateSizeFactors(ps2.tmpdds, geoMeans=geoMeans); rownames(ps2.tmpdds) <- ps2.tmp.tax@tax_table[,7]
  ps2.tmpdds <- estimateDispersions(ps2.tmpdds)
  ps2.tmpdds <- DESeq(ps2.tmpdds, fitType="local")
  ps2.tmpdds_res <- as.data.table(results(ps2.tmpdds)); ps2.tmpdds_res[,id:=ps2.tmp.tax@tax_table[,7]]
	ps2.tmpdds_res <- merge(ps2.tmpdds_res, data.table(id=ps2.taxmap$taxon_names(), taxon.id=names(ps2.taxmap$taxon_names())), by="id") # add taxmap ids
	return(ps2.tmpdds_res)
}

ps2.deseq2.difftable <- data.table()
ps2.asso.worm  <- subset_samples(ps2, source %in% c("host","associated"))
ps2dds.asso.worm.tax  <- deseq_diffspec_taxmap(ps2.asso.worm)
ps2.deseq2.difftable <- rbind(ps2.deseq2.difftable, data.table(taxon_id=ps2dds.asso.worm.tax$taxon.id, treatment_2="associated", treatment_1="host", log2FoldChange=ps2dds.asso.worm.tax$log2FoldChange, baseMean=ps2dds.asso.worm.tax$baseMean, padj=ps2dds.asso.worm.tax$padj))

ps2.alone.asso <- subset_samples(ps2, source %in% c("alone","associated") & !time %in% c(0, 2)) # 2 more time points in cembio.alone
ps2dds.alone.asso <- deseq_diffspec_taxmap(ps2.alone.asso)	
if(nrow(ps2dds.alone.asso[padj<=0.05]) > 0){ # alone vs. asso has now differential abundany species
	ps2.deseq2.difftable <- rbind(ps2.deseq2.difftable, data.table(taxon_id=ps2dds.alone.asso$taxon.id, treatment_2="alone", treatment_1="associated", log2FoldChange=ps2dds.alone.asso$log2FoldChange, baseMean=ps2dds.alone.asso$baseMean, padj=ps2dds.alone.asso$padj))
}
								 
ps2.alone.worm <- subset_samples(ps2, source %in% c("alone","host") & !time %in% c(0, 2)) # 2 more time points in cembio.alone
ps2dds.alone.worm <- deseq_diffspec_taxmap(ps2.alone.worm)								 
ps2.deseq2.difftable <- rbind(ps2.deseq2.difftable, data.table(taxon_id=ps2dds.alone.worm$taxon.id, treatment_2="alone", treatment_1="host", log2FoldChange=ps2dds.alone.worm$log2FoldChange, baseMean=ps2dds.alone.worm$baseMean, padj=ps2dds.alone.worm$padj))	

# diff abundance test from metacoder package does not account for 'time' and 'source', use deseq2 instead								 
ps2.taxmap$data$diff_table <- tibble::tibble(ps2.deseq2.difftable)
ps2.taxmap$data$diff_table$log2FoldChange[ps2.taxmap$data$diff_table$padj > 0.05 | is.na(ps2.taxmap$data$diff_table$log2FoldChange)] <- 0
ps2.taxmap$data$diff_table$baseMean[ps2.taxmap$data$diff_table$padj > 0.05 | is.na(ps2.taxmap$data$diff_table$baseMean)] <- 0
								 
ps2.taxmap %>%
  heat_tree_matrix(data = "diff_table",
                   #node_label = cleaned_names,
                   node_label = taxon_names,
                   node_size = baseMean, #n_obs, # number of OTUs
                   node_color = log2FoldChange, #log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "base mean",
                   node_color_axis_label = "Log2 fold change",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2) 
ggplot2::ggsave("../img/heat_tree-differential.pdf", width=5.5, height=5)

