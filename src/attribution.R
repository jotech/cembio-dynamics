library(phyloseq)
library(ggplot2)
library(data.table)
library(stringr)
library(glmnet)
library(Boruta)
library(foreach)
library(doParallel)

cores <- 8
registerDoParallel(cores)
ps2 <- readRDS("../dat/phyloseq_ps2.RDS")
ps2 <- subset_samples(ps2, !time %in% c(0,2)) # sync time
levels(ps2@sam_data$source) <- c("control", "substrate", "host") # rename sample sources
ps2 <- subset_taxa(ps2, !id %in% c("unknown","OP50")) # remove non-cembio members
rownames(ps2@otu_table) <- make.names(rownames(ps2@otu_table)) # harmonize otu names

ps2.combined.feat.dds <- readRDS("../dat/diff-functions.RDS")
ps2.feat <- readRDS("../dat/phyloseq_ps2-feat.RDS")

dat.pwy.mat <- read.table("~/uni/cembio.ext/dat/agnes/pwy-org-table.tbl")
dat.pwy.mat[,1]
ps2.rel <- transform_sample_counts(ps2, function(x){x / sum(x)})

feat.mat <- merge_pwy(c("MYb191.MYb177", "MYb71.MYb49", "MYb396.MYb69.MYb21", "MYb371.MYb331.MYb330.MYb177", "MYb176.MYb174"), dat.pwy.mat[,-1])
feat.attr.mat <- matrix(0, nrow=nrow(ps2.rel@otu_table), ncol=ncol(feat.mat), dimnames=list(rownames(ps2.rel@otu_table),colnames(feat.mat)))
org.idx <- match(rownames(ps2.rel@otu_table), rownames(feat.mat))
for(i in 1:ncol(feat.mat)){
    feat.attr.mat[,i] <- (ps2.rel@otu_table * feat.mat[org.idx,i])[,1]
}


# regression

ps2.feat.uni.id <- unique(ps2.combined.feat.dds[,.(id,name,subsystem)])
attr.dt <- foreach(i=1:nrow(ps2.feat.uni.id), .combine=rbind) %dopar%{
    cat("",i,"/",ncol(ps2.feat.uni.id))
    response.id=ps2.feat.uni.id[i,id]
    response.name=ps2.feat.uni.id[i,name]
    response.subsystem=ps2.feat.uni.id[i,subsystem]
    response  <- as.vector(ps2.feat@otu_table[match(response.id, rownames(ps2.feat@otu_table)),])
    covariats <- t(ps2@otu_table)

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
group.dt <- data.table(org=factor(make.names(ps2@tax_table[,"id"])), group=factor(ps2@tax_table[,"Genus"]))[order(group)]
attr.dt[,org:=factor(org, levels=group.dt$org)]
org.col <- BacArena::colpal3[group.dt$group]

ggplot(attr.dt, aes(x=org,y=lasso.coef)) + geom_boxplot() + facet_wrap(~subsystem,ncol=1,scales="free_y") + theme_minimal(base_size=14) + xlab("") + ylab("Regression coefficient (lasso)") + theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1, color=org.col))
ggsave("../img/attr-lasso.pdf", width=7, height=12)

ggplot(attr.dt, aes(x=org,y=boruta.imp)) + geom_boxplot() + facet_wrap(~subsystem,ncol=1,scales="free_y") + theme_minimal(base_size=14) + xlab("") + ylab("Importance score (Boruta)") + theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1, color=org.col))
ggsave("../img/attr-boruta.pdf", width=7, height=12)


ggplot(attr.dt[subsystem=="uast" & lasso.coef!=0], aes(y=org,x=feat,fill=lasso.coef)) + geom_raster() + theme_minimal(base_size=14) + xlab("") + ylab("") + labs(fill="regression\n coeficient") + theme(axis.text.x=element_text(angle=45, hjust=1), strip.background = element_blank(),axis.ticks = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + scale_fill_gradient(low="white", high="black")
ggsave("../img/attr_uast-lasso.pdf", width=4.5, height=6)

ggplot(attr.dt[subsystem=="uast" & boruta.imp!=0], aes(y=org,x=feat,fill=boruta.imp)) + geom_raster() + theme_minimal(base_size=14) + xlab("") + ylab("") + labs(fill="Boruta\nimportance") + theme(axis.text.x=element_text(angle=45, hjust=1), strip.background = element_blank(),axis.ticks = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + scale_fill_gradient(low="white", high="black")
ggsave("../img/attr_uast-boruta.pdf", width=4.5, height=6)


