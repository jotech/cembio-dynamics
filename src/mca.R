library(sybil)
library(FactoMineR)
library(ggplot2)
library(ggrepel)

models <- readRDS("../dat/cembio.ext_20230818_premed.RDS")
taxonomy.dt <- fread("../dat/taxonomy.csv")
colors.tax <- fread("../dat/colors_taxa.csv")
taxonomy.dt <- merge(taxonomy.dt, colors.tax, by.x="genus", by.y="group", all.x=T)
taxonomy.dt[is.na(color), color:=setdiff(BacArena::colpal3,color)[1:.N]]

rxn.dt <- data.table()
for(m in models){
	rxn.dt <- rbind(rxn.dt, data.table(id=m@mod_id,rxn=m@react_id))
}
rxn.dt <- rxn.dt[!id %in% c("MYbb2","MYbb4", "OP50", "MYb45")]
rxn.large.dt <- data.table(dcast(rxn.dt, id~rxn, fun.aggregate=length))
rxn.mat <- as.matrix(rxn.large.dt[,-"id"]); rownames(rxn.mat) <- rxn.large.dt$id
rxn.mat <- rxn.mat[,apply(rxn.mat, 2, function(x){(length(unique(x))>1)})]

rxn.mat.bool <- ifelse(rxn.mat==0,F,T)
pca <- MCA(rxn.mat.bool, graph = F)
pca.dt <- data.table(dim1=pca$ind$coord[,1], dim2=pca$ind$coord[,2], id=rownames(pca$ind$coord), group=taxonomy.dt$genus[match(rownames(pca$ind$coord), taxonomy.dt$id)])
ggplot(pca.dt, aes(x=dim1, y=dim2, color=group)) + geom_point(size=3) + geom_text_repel(size=3, aes(label=id), bg.color = "black", bg.r = .01, max.overlaps=15) + xlab(paste0("Dim 1 (", round(pca$eig[1,2],1),"%)")) + ylab(paste0("Dim 2 (", round(pca$eig[2,2],1),"%)")) + theme_minimal(base_size = 14) + scale_color_manual(values=setNames(taxonomy.dt$color,taxonomy.dt$genus)) + labs(color="")
#ggsave("../img/cembio-mca_models-rxn.pdf", width=7,height=5)

