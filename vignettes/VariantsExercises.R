
## ----setup,results='hide',echo=FALSE-------------------------------------
library(knitr)
opts_chunk$set(cache=TRUE,out.width='0.5\\textwidth',fig.align='center',warnings=FALSE,messages=FALSE)


## ----locateVariants,results="hide"---------------------------------------
library(VariantAnnotation)
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(fl, genome="hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
codingvar <- locateVariants(vcf, txdb, CodingVariants())
seqlevels(vcf)=paste0('chr',seqlevels(vcf))
codingvar <- locateVariants(vcf, txdb, CodingVariants())
head(codingvar, 3)
allvar <- locateVariants(vcf, txdb, AllVariants())
head(allvar)


## ----predictCoding,results="hide"----------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
coding <- predictCoding(vcf, txdb, seqSource=Hsapiens)
coding[5:7]


## ----predictCoding_frameshift,results="hide"-----------------------------
## CONSEQUENCE is 'frameshift' where translation is not possible
coding[mcols(coding)$CONSEQUENCE == "frameshift"]


## ----gviz10,fig.show=TRUE,results='hide',fig.align='center',error=FALSE,warning=FALSE----
library(Gviz)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
geneTrack = GeneRegionTrack(txdb,chromosome='chr22',start=50e6,end=51e6)
variantTrack = AnnotationTrack(rowData(vcf))
plotTracks(list(geneTrack,variantTrack),from=50.5e6,to=50.6e6)


## ----readvcf,results='hide'----------------------------------------------
# you may need to change the file name to load in the data
vcfMel = readVcf('../data/SK-MEL-5.vcf.gz',genome='hg19')
vcfLung = readVcf('../data/A549_ATCC.vcf.gz',genome='hg19')
nrow(vcfMel)
nrow(vcfLung)
vcfMel
vcfLung
rowData(vcfMel)
info(vcfMel)
rowData(vcfMel)
header(vcfMel)
ref(vcfMel)
alt(vcfMel)


## ----functiondefs,results='hide'-----------------------------------------
vcf2snpDF = function(vcf) {
    refall = as.character(ref(vcf))
    altall = as.character(unlist(alt(vcf))[start(PartitioningByEnd(alt(vcf)))])
    tmpDF = data.frame(ref=refall,alt=altall,stringsAsFactors=FALSE)
    # snps only
    tmpDF = tmpDF[nchar(tmpDF$ref)==1 & nchar(tmpDF$alt)==1,]
    return(tmpDF)
}


## ----titvfunction,results='hide'-----------------------------------------
titv = function(vcf) {
    variantDF = vcf2snpDF(vcf)
    tbl = table(variantDF)
    ti = tbl['C','T']+tbl['T','C']+tbl['A','G']+tbl['G','A']
    tv = tbl['A','C']+tbl['C','A']+tbl['A','T']+tbl['T','A']+tbl['C','G']+tbl['G','C']+tbl['G','T']+tbl['T','G']
    return(list(ti=ti,tv=tv,tbl=tbl,titv=ti/tv))
}


## ----subsetvcfs----------------------------------------------------------
vcfMelT2 = vcfMel[grep('Type2',rowData(vcfMel)$FILTER)]
vcfLungT2 = vcfLung[grep('Type2',rowData(vcfLung)$FILTER)]


## ----titvcald,results='hide'---------------------------------------------
titv(vcfMelT2)
titv(vcfLungT2)


## ----ggplot--------------------------------------------------------------
library(ggplot2)
ggplot(vcf2snpDF(vcfMelT2),aes(x=ref,fill=alt))+geom_bar()
ggplot(vcf2snpDF(vcfLungT2),aes(x=ref,fill=alt))+geom_bar()


