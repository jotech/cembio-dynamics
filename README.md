Script and data for C. elegenas mirobiome time series analysis

# Files
## src
- `abundances.R`: abundances plots (barplots, time series)
- `alpha-diversity.R`: within-sample diversity analysis
- `attribution.R`: contribution of species to functions
- `beta-diversity.R`: between-sample diversity analysis
- `diff-abundance.R`: differentially abundant species (DESeq2)
- `diff-functions.R`: differentially abundant functions (DESeq2)
- `heattree.R`: heattree showing differentially abundant species
- `mca.R`: mutiple correspondence analysis of metabolic reactions
- `nst.R`: normalized stochasticity ratio (null model)
- `over-representation-analysis.R`: Over-representation analysis of differentially abundant metabolic pathways 
- `phyloseq.R`: analysis of microbial census using phyloseq (filtering, agglomeration)
- `subsystems_boxplot.R`: Differentially abundant function summarized in subsystems determined by DESeq2


## dat
- `CE-CeMbio_V1_revised_20220420.xlsx`: Meta data experiment
- `sample.sheet.tsv`: Meta data sequencing
- `ASV_table.tsv`: ASV/OTU table with read counts
- `tax.own.tbl`: Taxonomic classification of ASVs
- `phyloseq_ps2.RDS`: phyloseq object used for further analysis
- `gtdbtk.bac120.user_msa.fasta.gz_aggregated-rooted.treefile`: phylogenetic tree file
- `diff-functions.csv`: differentially abundant functions
- `cembio.ext_20230818_premed.RDS`: metabolic models produced by gapseq
- `NST.RData`: normalized stochasticity ratio results (`NST_time.RData`, `NST_source.RData`)
- `colors_taxa.csv`: color codes
- `uast-*.tbl`, `gutsmash-*.tbl`, `medium-*.tbl`, `fermcs-*.tbl`, `meta_pwy.tbl`, `interact-tabale_count.tbl`, `dbcan-*.tbl`: predicted functions for each species
