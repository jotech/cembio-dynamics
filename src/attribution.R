library(phyloseq)
library(ggplot2)
library(data.table)
library(stringr)
library(glmnet)
library(Boruta)
library(foreach)
library(doParallel)

cores <- 16
registerDoParallel(cores)

ps2 <- readRDS("../dat/phyloseq_ps2.RDS")
ps2 <- subset_samples(ps2, !time %in% c(0,2)) # sync time
levels(ps2@sam_data$source) <- c("control", "substrate", "host") # rename sample sources
ps2 <- subset_taxa(ps2, !id %in% c("unknown","OP50")) # remove non-cembio members
ps2.rel <- transform_sample_counts(ps2, function(x){x / sum(x)})
rownames(ps2@otu_table) <- make.names(rownames(ps2@otu_table)) # harmonize otu names
rownames(ps2.rel@otu_table) <- make.names(rownames(ps2.rel@otu_table)) # harmonize otu names

ps2.combined.feat.dds <- readRDS("../dat/diff-functions.RDS")
ps2.feat <- readRDS("../dat/phyloseq_ps2-feat.RDS")

group.dt <- data.table(org=factor(make.names(ps2@tax_table[,"id"])), group=factor(ps2@tax_table[,"Genus"]))[order(group)]
colors.taxa <- fread("../dat/colors_taxa.csv")
group.dt <- merge(group.dt,colors.taxa)


#
# cumulative attribution
#
dat.pwy.mat <- read.table("../dat/pwy-org-table.tbl")
dat.abricate.mat <- read.table("../dat/abricate-org-table.tbl")
dat.gut.mat <- read.table("../dat/gutsmash-org-table.tbl")
dat.dbcan.mat <- read.table("../dat/dbcan-org-table.tbl")
dat.medium.mat <- read.table("../dat/medium-org-table.tbl")
dat.uast.mat <- read.table("../dat/uast-org-table.tbl")
dat.fermcs.mat <- read.table("../dat/fermcs-org-table.tbl")
all.equal(sort(rownames(dat.pwy.mat)), sort(rownames(dat.abricate.mat)),sort(rownames(dat.gut.mat)),sort(rownames(dat.dbcan.mat)), sort(row.names(dat.medium.mat)), sort(rownames(dat.uast.mat)), sort(rownames(dat.fermcs.mat)))

dat.pwy.mat <- dat.pwy.mat[order(rownames(dat.pwy.mat)),]
dat.abricate.mat <- dat.abricate.mat[order(rownames(dat.abricate.mat)),]
dat.gut.mat <- dat.gut.mat[order(rownames(dat.gut.mat)),]
dat.dbcan.mat <- dat.dbcan.mat[order(rownames(dat.dbcan.mat)),]
dat.medium.mat <- dat.medium.mat[order(rownames(dat.medium.mat)),]
dat.uast.mat <- dat.uast.mat[order(rownames(dat.uast.mat)),]
dat.fermcs.mat <- dat.fermcs.mat[order(rownames(dat.fermcs.mat)),]

all.equal(rownames(dat.pwy.mat), rownames(dat.abricate.mat),rownames(dat.gut.mat),rownames(dat.dbcan.mat), row.names(dat.medium.mat), rownames(dat.uast.mat), rownames(dat.fermcs.mat))
dat.feat.mat <- cbind(dat.pwy.mat, dat.abricate.mat[match(rownames(dat.abricate.mat),rownames(dat.pwy.mat)),],dat.gut.mat[match(rownames(dat.gut.mat),rownames(dat.pwy.mat)),],dat.dbcan.mat[match(rownames(dat.dbcan.mat),rownames(dat.pwy.mat)),], dat.medium.mat[match(rownames(dat.medium.mat),rownames(dat.pwy.mat)),], dat.uast.mat[match(rownames(dat.uast.mat),rownames(dat.pwy.mat)),], dat.fermcs.mat[match(rownames(dat.fermcs.mat),rownames(dat.pwy.mat)),])

merge_pwy <- function(merge.lst, dat.mat){
	new.rows <- matrix(nrow=0, ncol=ncol(dat.mat), dimnames=list(NULL,colnames(dat.mat)))
	idx.lst <- c()
	for(pair in merge.lst){
	  idx <- match(unlist(str_split(pair,"\\.")), rownames(dat.mat))
		print(pair)
		print(idx)
		col.merged <- as.logical(apply(dat.mat[idx,],2,max))
		new.rows <- rbind(new.rows, col.merged)
		rownames(new.rows)[nrow(new.rows)] <- pair
		idx.lst <- c(idx.lst, idx)
	}
	return(rbind(dat.mat[-unique(idx.lst),], new.rows))
}

feat.mat <- merge_pwy(c("MYb191.MYb177", "MYb71.MYb49", "MYb396.MYb69.MYb21", "MYb371.MYb331.MYb330", "MYb176.MYb174"), dat.feat.mat)
cum.attr.mat <- matrix(0, nrow=nrow(ps2.rel@otu_table), ncol=ncol(feat.mat), dimnames=list(rownames(ps2.rel@otu_table),colnames(feat.mat)))
org.idx <- match(rownames(ps2.rel@otu_table), rownames(feat.mat))
for(i in 1:ncol(feat.mat)){
    cum.attr.mat[,i] <- (ps2.rel@otu_table * feat.mat[org.idx,i])[,1]
}

cum.attr.dt <- data.table(reshape2:::melt.matrix(cum.attr.mat, varnames=c("org","id"),value.name="cumulative"))

cum.attr.dt <- cum.attr.dt[id %in% ps2.combined.feat.dds$id]
cum.attr.dt <- cbind(cum.attr.dt, ps2.combined.feat.dds[,.(name,subsystem)][match(cum.attr.dt$id, make.names(ps2.combined.feat.dds$id)),])
cum.attr.dt[is.na(cumulative),cumulative:=0]

cum.attr.dt[,org:=factor(org, levels=group.dt$org)]
levels(cum.attr.dt$org) <- gsub("\\.",",",levels(cum.attr.dt$org))
ggplot(cum.attr.dt, aes(x=org,y=cumulative)) + geom_boxplot() + facet_wrap(~subsystem,nrow=1,scales="free_x") + coord_flip() + theme_minimal(base_size=14) + xlab("") + ylab("Cumulative abundance") + theme(axis.text.y = element_text(size=10,color=group.dt$color))
ggsave("../img/attr-cumulative.pdf", width=12, height=7)


#
# regression/random forest attribution
#
ps2.feat.uni.id <- unique(ps2.combined.feat.dds[,.(id,name,subsystem)])
attr.dt <- foreach(i=1:nrow(ps2.feat.uni.id), .combine=rbind) %dopar%{
    cat("",i,"/",ncol(ps2.feat.uni.id))
    response.id=ps2.feat.uni.id[i,id]
    response.name=ps2.feat.uni.id[i,name]
    response.subsystem=ps2.feat.uni.id[i,subsystem]
    response  <- as.vector(ps2.feat@otu_table[match(response.id, rownames(ps2.feat@otu_table)),])
    covariats <- t(ps2@otu_table)
    covariats <- cbind(data.table(time=as.numeric(ps2@sam_data$time)),covariats)
    covariats <- cbind(data.table(source=as.numeric(ps2@sam_data$source)),covariats)

    cor.dt <- rbindlist(apply(covariats, 2, function(x) {cor <- cor.test(x, response, method="spearman"); data.table(cor.pval=cor$p.value, cor.rho=cor$estimate)}), idcol="org")
    cor.dt[,cor.pval.adj:=p.adjust(cor.pval)]
    #cor.dt[cor.pval.adj<=0.05][order(cor.rho)]

    cv <- cv.glmnet(x=as.matrix(covariats), y=response, alpha=1) # estimate shrinkage parameter
    lasso <- glmnet(x=as.matrix(covariats), y=response, alpha = 1, lambda = cv$lambda.min)
    lasso.dt <- data.table(org=rownames(coef(lasso))[-1], lasso.coef=coef(lasso)[,1][-1])

    boruta <- Boruta(x=covariats, y=response)
    boruta.dt <- data.table(attStats(boruta), keep.rownames=T)
    boruta.dt2 <- data.table(org=boruta.dt$rn, boruta.imp=ifelse(boruta.dt$decision=="Confirmed", boruta.dt$meanImp, 0))

    data.table(feat=response.id, name=response.name, subsystem=response.subsystem, merge(cor.dt, merge(lasso.dt,boruta.dt2,by="org", all=T), by="org", all=T))
}
for (i in seq_along(attr.dt)) set(attr.dt, i=which(is.na(attr.dt[[i]])), j=i, value=0)
#saveRDS(attr.dt, "../dat/attribution.RDS", compress="xz")


#
# Evaluation
#
attr.dt <- readRDS("../dat/attribution.RDS")
attr.dt[,org2:=factor(org, levels=group.dt$org)]
levels(attr.dt$org2) <- gsub("\\.",",",levels(attr.dt$org2))

ggplot(attr.dt[!is.na(org2) & subsystem!="interactions"], aes(x=org2,y=lasso.coef)) + geom_boxplot() + facet_wrap(~subsystem,nrow=1,scales="free_x") + coord_flip() + theme_minimal(base_size=14) + xlab("") + ylab("Regression coefficient (lasso)") + theme(axis.text.y = element_text(size=10, color=group.dt$color)) + scale_y_continuous(n.breaks = 3)
ggsave("../img/attr-lasso.pdf", width=12, height=7)

ggplot(attr.dt[!is.na(org2)], aes(x=org2,y=boruta.imp)) + geom_boxplot() + facet_wrap(~subsystem,nrow=1,scales="free_x") + coord_flip() + theme_minimal(base_size=14) + xlab("") + ylab("Importance (random forest)") + theme(axis.text.y = element_text(size=10, color=group.dt$color))
#ggplot(attr.dt, aes(x=org2,y=boruta.imp)) + geom_boxplot() + facet_wrap(~subsystem,ncol=1,scales="free_y") + theme_minimal(base_size=14) + xlab("") + ylab("Importance score (Boruta)") + theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1, color=group.dt$color))
ggsave("../img/attr-boruta.pdf", width=12, height=7)


ggplot(attr.dt[subsystem=="uast" & lasso.coef!=0], aes(y=org2,x=feat,fill=lasso.coef)) + geom_raster() + theme_minimal(base_size=14) + xlab("") + ylab("") + labs(fill="regression\n coeficient") + theme(axis.text.x=element_text(angle=45, hjust=1), strip.background = element_blank(),axis.ticks = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + scale_fill_gradient(low="white", high="black")
ggsave("../img/attr_uast-lasso.pdf", width=4.5, height=6)

ggplot(attr.dt[subsystem=="uast" & boruta.imp!=0], aes(y=org2,x=feat,fill=boruta.imp)) + geom_raster() + theme_minimal(base_size=14) + xlab("") + ylab("") + labs(fill="Boruta\nimportance") + theme(axis.text.x=element_text(angle=45, hjust=1), strip.background = element_blank(),axis.ticks = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + scale_fill_gradient(low="white", high="black")
ggsave("../img/attr_uast-boruta.pdf", width=4.5, height=6)

# b12 & propionate
attr.dt[grepl("cobalamin|cobalt", name) & lasso.coef!=0, .(name,org,lasso.coef)][order(rank(name),lasso.coef)]
attr.dt[grepl("propionate", name) & (lasso.coef!=0 & boruta.imp!=0)]

# bcaa
attr.dt[grepl("leucine|valine", name) & (lasso.coef!=0 & boruta.imp!=0)][order(rank(name),lasso.coef)]
