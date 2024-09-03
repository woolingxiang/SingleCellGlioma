library(Seurat)
library(devtools)
library(harmony)
library(qusage)
library(GSVA)
library(ggplot2)

# loading data
tiantan.wt.all = readRDS('./tiantan.wt.all.anno.RDS')
tiantan.wt.all = readRDS('./tiantan.wt.tmr.anno.RDS')

# overview, cell types
DimPlot(tiantan.wt.all,reduction='umap',group.by='redefine_type',cols=c('#46AEA0','#FABF7B','#3C93C2','#F0746E','#DDDDDD'))
tmp = tiantan.wt.all@meta.data
tmp$newtype = tmp$redefine_type
tmp$cnt = 1
stat = aggregate(cnt ~ orig.ident+redefine_type,tmp,sum)
ggplot(stat, aes(fill=redefine_type, y=cnt, x=orig.ident)) + 
  geom_bar(position="fill", stat="identity")

# overview, features
FeaturePlot(tiantan.wt.all,reduction='umap',features=c('PTPRZ1','PTPRC','CD68','CD3D'),raster=F);
FeaturePlot(tiantan.wt.all,reduction='umap',features=c('PTPRZ1','PTPRC','MOG','CD3D'),raster=F);
DotPlot(tiantan.wt.all,features=c('PTPRZ1','PTPRC','MAG','CD248','CD68','CDH5','CHI3L1'),group.by = 'redefine_type')

# overview, others
ggplot(tiantan.wt.all@meta.data,aes(y=nCount_RNA, x=orig.ident)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",fill='#F0746E') +
  theme_bw()
ggplot(tiantan.wt.all@meta.data,aes(y=nFeature_RNA, x=orig.ident)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent",fill='#FABF7B') +
  theme_bw()




