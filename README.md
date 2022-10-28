Script and data for C. elegenas mirobiome time series analysis

# Files
## src
- `phyloseq.R`: analysis of microbial census using phyloseq (filtering, agglomeration)
- `subsystems_boxplot.R`: Differentially abundant function summarized in subsystems determined by DESeq2
- `over-representation-analysis.R`: Over-representation analysis of differentially abundant metabolic pathways 

## dat
- `CE-CeMbio_V1_revised_20220420.xlsx`: Meta data experiment
- `sample.sheet.tsv`: Meta data sequencing
- `ASV_table.tsv`: ASV/OTU table with read counts
- `tax.own.tbl`: Taxonomic classification of ASVs
- `dds-diff-alone.asso.feat.RDS`: Differentially abundant features from DESeq2 for alone vs. associated samples
- `dds-diff-alone.worm.feat.RDS`: Differentially abundant features from DESeq2 for alone vs. host samples

- `dds-diff-asso.worm.feat.RDS`:  Differentially abundant features from DESeq2 for associated vs. host samples

