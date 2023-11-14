library(ggplot2)
library(clusterProfiler)

ps2.feat.dds <- readRDS("../dat/diff-functions_all-cmp.RDS")
ps2.asso.worm.feat.dds  <- ps2.feat.dds[[1]] 
ps2.alone.asso.feat.dds <- ps2.feat.dds[[2]]
ps2.alone.worm.feat.dds <- ps2.feat.dds[[3]]

biosynthesis <- c("Aminoacyl-tRNAs-Charging","AROMATIC-COMPOUNDS-BIOSYN","Carbohydrates-Biosynthesis","Cell-Structure-Biosynthesis","Cofactor-Biosynthesis","Lipid-Biosynthesis","Metabolic-Regulators","Nucleotide-Biosynthesis","Other-biosynthesis","Polyprenyl-Biosynthesis","SECONDARY-METABOLITE-BIOSYNTHESIS","Storage-Compounds-Biosynthesis","Tetrapyrrole-Biosynthesis")
degradation <- c("Alcohol-Degradation","Aldehyde-Degradation","AMINE-DEG","Amino-Acid-Degradation","AROMATIC-COMPOUNDS-DEGRADATION","C1-COMPOUNDS","Carbohydrates-Degradation","CARBOXYLATES-DEG","CHLORINATED-COMPOUNDS-DEG","COFACTOR-DEGRADATION","Other-Degradation","Fatty-Acid-and-Lipid-Degradation","HORMONE-DEG","Noncarbon-Nutrients","NUCLEO-DEG","Polymer-Degradation", "Protein-Degradation","SECONDARY-METABOLITE-DEGRADATION")
energy <- c("Acetyl-CoA-Biosynthesis","CHEMOAUTOTROPHIC-ENERGY-METABOLISM","Electron-Transfer","Entner-Duodoroff-Pathways","Fermentation","GLYCOLYSIS-VARIANTS","Hydrogen-Production","OTHER-ENERGY", "Pentose-Phosphate-Cycle","Photosynthesis","Respiration","TCA-VARIANTS")
detox <- c("8-Oxo-GTP-Detoxification","Acid-Resistance","Antibiotic-Resistance","Arsenic-Detoxification","Cyanide-Detoxification","Mercury-Detoxification","Methylglyoxal-Detoxification","REACTIVE-OXYGEN-SPECIES-DEGRADATION","Seleno-Amino-Acid-Detoxification")
glycan <- c("Glycan-Biosynthesis","Glycan-Degradation")
macro <- c("Nucleic-Acid-Processing","Protein-Modification")
activation <- c("Activation","Inactivation","Interconversion")

ora_fun <- function(feat.dds, dir){
	meta.sub.db.dt <- feat.dds[!is.na(hierarchy), lapply(.SD, function(x) unlist(tstrsplit(x, ",", fixed=TRUE))), by = .(id, name), .SDcols = "hierarchy"]; meta.sub.db.dt[,`:=`(hierarchy=gsub("\\|","",hierarchy))]
	meta.sub.db.dt <- meta.sub.db.dt[hierarchy %in% c(biosynthesis, degradation,energy,detox,glycan,macro,activation)]

	if(dir=="up"){
		feat.lst <- feat.dds[padj<=0.05 & !is.na(hierarchy) & log2FoldChange > 0, id]
	} else	feat.lst <- feat.dds[padj<=0.05 & !is.na(hierarchy) & log2FoldChange < 0, id]
	feat.ora <- enricher(feat.lst, TERM2GENE = meta.sub.db.dt[,.(hierarchy,id)])
	return(feat.ora)																	  
}

														 
ora.dt <- rbind(data.table(id="host vs. substrate", data.frame(ora_fun(ps2.asso.worm.feat.dds, dir="up"))),
	  data.table(id="substrate vs. host", data.frame(ora_fun(ps2.asso.worm.feat.dds, dir="down"))),
	  data.table(id="host vs. control", data.frame(ora_fun(ps2.alone.worm.feat.dds, dir="up"))),
	  data.table(id="control vs. host", data.frame(ora_fun(ps2.alone.worm.feat.dds, dir="down"))))
	  #data.table(id="substrate vs. control", data.frame(ora_fun(ps2.alone.asso.feat.dds, dir="up"))),
	  #data.table(id="control vs. substrate", data.frame(ora_fun(ps2.alone.asso.feat.dds, dir="down"))))

print(ora.dt[!is.na(ID), -"geneID"])
