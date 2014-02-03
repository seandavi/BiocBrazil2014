
## ----setup,results='hide',echo=FALSE-------------------------------------
library(knitr)
opts_chunk$set(cache=TRUE,out.width='0.5\\textwidth')


## ----geoquery,echo=FALSE,results='hide',message=FALSE,warnings=FALSE-----
library(GEOquery)
eset = getGEO('GSE1577')[[1]]


## ----featureData,echo=TRUE,results='hide'--------------------------------
head(fData(eset))


## ----featureData2,echo=FALSE---------------------------------------------
head(fData(eset)[,c(11,12,13,10)])


## ----phenoData,echo=TRUE,results='hide'----------------------------------
head(pData(eset))


## ----phenoData2,echo=FALSE-----------------------------------------------
pData(eset)[c(1:3,10:12),c(1:2,4,6:8,10:11)]


## ----assayData,echo=TRUE,results='hide'----------------------------------
head(assayDataElement(eset,'exprs')) 
# OR 
head(exprs(eset))


## ----assayData2,echo=FALSE-----------------------------------------------
head(assayDataElement(eset,'exprs')) # head(exprs(eset))


## ----getgse,messages=FALSE,warnings=FALSE--------------------------------
browseURL("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1577")
library(GEOquery)
eset = getGEO('GSE1577')[[1]]


## ----exprs1--------------------------------------------------------------
class(eset)


## ----exprs2,results='hide'-----------------------------------------------
help(ExpressionSet)
help('ExpressionSet-class')


## ----exprs10-------------------------------------------------------------
show(eset)


## ----exprs20-------------------------------------------------------------
# dimensions of the ExpressionSet
dim(eset)
# Number of features
nrow(eset) 
# number of samples
ncol(eset) 


## ----exprs30-------------------------------------------------------------
# "names" of samples--must be unique
sampleNames(eset)[1:8]
# "names" of features--must be unique
featureNames(eset)[1:8]


## ----phenoData22,echo=TRUE,results='hide'--------------------------------
head(pData(eset))
class(pData(eset))
colnames(pData(eset))
all.equal(sampleNames(eset),rownames(pData(eset)))
summary(pData(eset))


## ----featureData22,echo=TRUE,results='hide'------------------------------
head(fData(eset))
class(fData(eset))
colnames(fData(eset))
all.equal(featureNames(eset),rownames(fData(eset)))
summary(fData(eset))


## ----assayData22,echo=TRUE,results='hide'--------------------------------
head(assayDataElement(eset,'exprs')) 
# OR 
head(exprs(eset))
class(exprs(eset))
summary(exprs(eset))
all.equal(colnames(exprs(eset)),sampleNames(eset))
all.equal(rownames(exprs(eset)),featureNames(eset))


## ----subsetting10,results='hide'-----------------------------------------
eset[1:10,]
eset[,1:10]


## ----subsetting20--------------------------------------------------------
levels(pData(eset)$source_name_ch1)


## ----subsetting30--------------------------------------------------------
eset[,pData(eset)$source_name_ch1=="Bone marrow sample"]


## ----heatmapprelim,results='hide'----------------------------------------
library(gplots)
?heatmap.2


## ----heatmap10,results='hide'--------------------------------------------
summary(exprs(eset))
hist(exprs(eset))


## ----heatmap20,results='hide'--------------------------------------------
summary(log2(exprs(eset)))
hist(log2(exprs(eset)))


## ----heatmap30,results='hide'--------------------------------------------
exprs(eset) = log2(exprs(eset))
hist(exprs(eset))


## ----heatmap40,results='hide'--------------------------------------------
stdDev = apply(exprs(eset),1,sd)
hist(stdDev)


## ----heatmap50,results='hide'--------------------------------------------
eset2 = eset[order(stdDev,decreasing=TRUE)[1:20],]


## ----heatmap60,results='hide',fig.show='hide'----------------------------
heatmap.2(exprs(eset2),trace='none')


## ----heatmap70,results='hide',echo=FALSE,fig.show='hide'-----------------
heatmap.2(exprs(eset2),trace='none',labRow=fData(eset)$'Gene Symbol',
  labCol=pData(eset)$source_name_ch1)


