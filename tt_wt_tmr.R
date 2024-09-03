library(Seurat)
library(devtools)
library(harmony)
library(qusage)
library(GSVA)
library(ggplot2)

# loading data
tiantan.wt.all = readRDS('./tiantan.wt.all.anno.RDS')
tiantan.wt.tmr = readRDS('./tiantan.wt.tmr.anno.RDS')

# overview, cell types - tumor subpopulations
tmp = tiantan.wt.all@meta.data[which(tiantan.wt.all$redefine_type!='malignant-like'),]
tmp$newtype = tmp$redefine_type
tmp$cnt = 1
stat1 = aggregate(cnt ~ orig.ident+newtype,tmp,sum)
tmp = tiantan.wt.tmr@meta.data
tmp$newtype = tmp$redefine_cluster
tmp$cnt = 1
stat2 = aggregate(cnt ~ orig.ident+newtype,tmp,sum)
stat = rbind(stat1,stat2)
ggplot(stat, aes(fill=newtype, y=cnt, x=orig.ident)) + 
  geom_bar(position="fill", stat="identity")

# overview, tumor subpopulations
DimPlot(tiantan.wt.tmr,reduction='umap',group.by='redefine_cluster',cols=c('#46AEA0','#FABF7B','#3C93C2','#F0746E','#9CCB86'))
DotPlot(tiantan.wt.tmr,group.by='redefine_cluster',features=c('CHI3L1','DLL3','CENPF','VEGFA','RND3'))

# PCA
hmks = qusage::read.gmt('./h.all.v6.1.symbols.gmt')
mat = as.matrix(tiantan.wt.tmr[['RNA']]@data)
ssgsea.obj = ssgseaParam(mat,hmks)
tiantan.wt.tmr.hmks <- gsva(ssgsea.obj)
subcluster.hmks = apply(tiantan.wt.tmr.hmks,1,function(x,y){tapply(x,y,mean)},y=tiantan.wt.tmr$redefine_cluster)
pca = prcomp(subcluster.hmks)
plot(pca$x[,1],pca$x[,2],cex=3)
text(pca$x[,1],pca$x[,2],label=rownames(subcluster.hmks))

# ICMs
persample_stat = function(gene_vec,smpl,grp,mtd='frac'){
  grp = as.factor(grp)
  sample_vec = unique(smpl)
  stat = sapply(sample_vec,function(x){
    idx = which(smpl==x)
    tapply(gene_vec[idx],grp[idx],function(a){
      if(mtd=='frac'){
        length(which(a>0))
      }else if(mtd=='avg'){
        mean(a,na.rm=T)
      }
    })
  })
  stat[is.na(stat)] = 0
  stat
}
mat = as.matrix(tiantan.wt.tmr[['RNA']]@data)
pool = sapply(c('PVR','CD274','PDCD1LG2','LGALS9','CD80','CD86','CD200'),function(g){
  stat = persample_stat(unlist(mat[g,]),tiantan.wt.tmr$orig.ident,tiantan.wt.tmr$redefine_cluster)
  apply(t(stat)/colSums(stat),2,function(x){mean(x[!is.nan(x)])})
})
cols = colorRampPalette(c('#FCDE9C','#FFFFFF','#F0746E'))
pheatmap::pheatmap(pool,scale='column',color=cols(20))

# focus on tmr_chi3l1
pool = apply(mat[c('PVR','CD274','PDCD1LG2','LGALS9','CD80','CD86'),which(tiantan.wt.tmr$subcluster=='tmr_CHI3L1')],1,function(x){
  ifelse(x>0,1,0)
})
pool.f = pool[which(rowSums(pool)>0),]; dim(pool.f)
head(pool.f)
stat = table(apply(pool.f,1,function(x){paste0(x,collapse='')}))
stat = stat/sum(stat)
barplot(stat,las=2)
stat

# cell communications
cccdir = grep('out$',list.dirs('./ccc'),value=T)
pool = matrix(NA,nrow=length(cccdir),ncol=7)
pool = as.data.frame(pool)
colnames(pool) = c('tmr_CHI3L1.CD4_CCR7','tmr_CHI3L1.CD4_CD40LG','tmr_CHI3L1.CD4_FOXP3','tmr_CHI3L1.CD8_FGFBP2','tmr_CHI3L1.CD8_GZMK','tmr_CHI3L1.CD8_MKI67','tmr_CHI3L1.NK_FGFBP2')
for(i in 1:length(cccdir)){
  cat(i)
  demo = read.table(paste0(cccdir[i],'/pvalues.txt'),header=T,sep='\t')
  demo.f = demo[grep('PVR',demo$gene_a),c(5,6,grep('^tmr_CHI3L1.(CD|NK)',colnames(demo)))]
  idx.name = colnames(demo.f)[-c(1:2)]
  if(length(which(demo.f$gene_b=='TIGIT'))>0){
    pool[i,idx.name] = unlist(demo.f[which(demo.f$gene_a=='PVR' & demo.f$gene_b=='TIGIT'),-c(1:2)])
  }
}
rownames(pool) = basename(dirname(cccdir))
pool.f = pool
pool.f[pool.f==0] = 0.001
cols = colorRampPalette(c('#FCDE9C','#FFFFFF','#F0746E'))
pheatmap::pheatmap(as.matrix(-log(pool.f)),cluster_row=F,cluster_col=F,color=cols(50))




