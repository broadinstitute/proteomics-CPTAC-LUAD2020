options(stringsAsFactors = F)
rm(list=ls())
library(pacman)
p_load(cmapR)
p_load(dplyr)
p_load(glue)

setwd('c:/Users/karsten/Dropbox/Devel/CPTAC-LUAD2019/')

## import gct files
data.dir <- 'g:/Shared drives/CPTAC3.0_LUAD_Data/luad-v3.1-data-freeze/'
label <- 'v3.1-tumor-over-nat'

ord.column <- 'NMF.consensus'


gct.str <- c(glue("{data.dir}luad-v3.1-acetylome-ratio-norm-NArm.gct"),
             glue("{data.dir}luad-v3.1-cnv-gene-LR.gct"),
             glue("{data.dir}luad-v3.1-phosphoproteome-ratio-norm-NArm.gct"),
             glue("{data.dir}luad-v3.1-proteome-ratio-norm-NArm.gct"),
             glue("{data.dir}luad-v3.1-rnaseq-prot-uq-rpkm-log2-NArm-row-norm.gct")
             )

gct <- lapply(gct.str, parse.gctx)
names(gct) <- c('5_acK', '1_CNA', '4_pSTY', '3_Protein', '2_RNAseq')


###########################################
## calculate tumor-normal ratios
gct.expr <- vector('list', length(gct))
names(gct.expr) <- names(gct)
gct.rdesc <- gct.cdesc <- gct.rid <- gct.cid <- gct.expr
  
for(i in names(gct)){
  cdesc <- gct[[i]]@cdesc
  rdesc <- gct[[i]]@rdesc
  rid <- gct[[i]]@rid 
  cid <- gct[[i]]@cid
  mat <- gct[[i]]@mat
  
  if('NAT' %in% cdesc$Type){
    nat.idx <- which(cdesc$Type == 'NAT')
    tumor.idx <- which(cdesc$Type == 'Tumor')
    
    ## separate tumor and nat
    mat.n <- mat[, nat.idx]
    mat.t <- mat[, tumor.idx]
    
    ## - exclude tumors without nat
    ## - bring into same order
    n.cid <- colnames(mat.n)
    t.cid <- colnames(mat.t)
    t.cid <- t.cid[ match( sub('\\.N$','', n.cid), t.cid ) ]
    mat.t <- mat.t[, t.cid]
    
    ## tumor-nat
    mat <- mat.t - mat.n
    cid <- t.cid
    cdesc <- cdesc[cid, ]
  }
  
  gct.expr[[i]] <- mat
  gct.rdesc[[i]] <- rdesc
  gct.cdesc[[i]] <- cdesc
  gct.rid[[i]] <- rid
  gct.cid[[i]] <- cid
}

## rdesc
#gct.rdesc <- lapply(gct, function(x)x@rdesc)

## expression
#gct.expr <- lapply(gct, function(x) x@mat)

## column ids
#gct.cid <-lapply(gct, function(x) x@cid)

## row ids
#gct.rid <-lapply(gct, function(x) x@rid)





## common sample names
samp.common <- Reduce(f = intersect, gct.cid)
gct.expr <- lapply(gct.expr, function(x)x[, samp.common])

###################
## rdesc 
gct.gene <- lapply(gct.rdesc, function(x) try(x[ ,grep('geneSymbol', colnames(x), value=T, ignore.case = T)[1]], silent=T))
#gct.gene <- lapply(gct.gene.id, )
if(sum(sapply(gct.gene, length) < 2) > 0){
  idx <- which( sapply(gct.gene, length) < 2 )
  for(i in idx)
    gct.gene[[i]] <- gct.rid[[i]]
}

## data type
gct.data.type <- lapply(names(gct.gene), function(x) rep(x, nrow(gct.expr[[x]])))
names(gct.data.type) <- names(gct.gene)

## ids
gct.id <- lapply(gct.rdesc, function(x)x[, 'id'])
#gct.id <- lapply(gct.rdesc, function(x) try(grep('id', colnames(x), value=T), silent=T))

if(sum(sapply(gct.id, length) < 2) > 0){
  idx <- which( sapply(gct.id, length) < 2 )
  for(i in idx)
    gct.id[[i]] <- gct.rid[[i]]
}

## combine
gct.rdesc <- lapply(names(gct.gene), function(x) data.frame(ID=gct.id[[x]], geneSymbol=gct.gene[[x]], DataType=gct.data.type[[x]]))

rdesc <- Reduce(f = rbind, gct.rdesc)

###############################
## cdesc from proteome
cdesc <- gct[['3_Protein']]@cdesc

#################################
## single data frames
column.anno <- cdesc[samp.common, ]
row.anno <- rdesc
tab.expr.all <- Reduce(f=rbind, gct.expr)
rownames(tab.expr.all) <- make.unique(rownames(tab.expr.all)) 
rownames(row.anno) <- rownames(tab.expr.all)

## reorder
ord.idx <- order(column.anno[, ord.column])
column.anno <- column.anno[ord.idx,]
tab.expr.all <- tab.expr.all[, ord.idx]

## reorder rows
row.idx <- with(row.anno, order(geneSymbol, DataType))
tab.expr.all <- tab.expr.all[row.idx, ]
row.anno <- row.anno[row.idx, ]

## export
save(column.anno, row.anno, tab.expr.all, file = glue('data/data-luad-{label}.RData'))

