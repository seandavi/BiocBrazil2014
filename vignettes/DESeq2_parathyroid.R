
## ----knitrsetup,echo=FALSE,hide=TRUE-------------------------------------
library(knitr)
opts_chunk$set(fig.lp='',size='footnotesize',dev=c('pdf','png'),fig.align='center',fig.width=4,fig.height=4.5,cache=TRUE)


## ----options, echo=FALSE-------------------------------------------------
options(digits=3, width=80)
#library("cacheSweave")
#setCacheDir("cachedir")


## ----libraries-----------------------------------------------------------
library( "DESeq2" )
library( "parathyroidSE" )


## ----loadEcs-------------------------------------------------------------
data("parathyroidGenesSE")


## ----countTable----------------------------------------------------------
head( assay( parathyroidGenesSE ) )


## ----nrowSE--------------------------------------------------------------
nrow(parathyroidGenesSE)


## ----colData-------------------------------------------------------------
colData( parathyroidGenesSE )


## ----rowData-------------------------------------------------------------
rowData( parathyroidGenesSE )


## ----columnData----------------------------------------------------------
as.data.frame( colData(parathyroidGenesSE)[,c("sample","patient","treatment","time")] )


## ----split---------------------------------------------------------------
allColSamples <- colData(parathyroidGenesSE)$sample
sp <- split( seq(along=allColSamples), allColSamples )


## ----addSamples----------------------------------------------------------
countdata <- sapply(sp, function(columns) 
   rowSums( assay(parathyroidGenesSE)[,columns,drop=FALSE] ) )
head(countdata)


## ----poormanssum---------------------------------------------------------
a <- assay(parathyroidGenesSE)
countdata2 <- cbind( a[,1:8], a[,9]+a[,10], a[,11], 
   a[,12]+a[,13], a[,14:22], a[,23]+a[,24], a[,25], a[,26]+a[,27] )
all( countdata == countdata2 )


## ----subsetMetaData------------------------------------------------------
coldata <- colData(parathyroidGenesSE)[sapply(sp, `[`, 1),]
rownames(coldata) <- coldata$sample
coldata


## ----colDataSubsetCols---------------------------------------------------
coldata <- coldata[ , c( "patient", "treatment", "time" ) ]
head( coldata )


## ----rowdata-------------------------------------------------------------
rowdata <- rowData(parathyroidGenesSE)
rowdata


## ----makeddsfull---------------------------------------------------------
ddsFull <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~ patient + treatment,
  rowData = rowdata)
ddsFull  


## ----testCountSum, echo=FALSE--------------------------------------------
stopifnot(sum(assay(parathyroidGenesSE)) == sum(counts(ddsFull)))


## ----subset--------------------------------------------------------------
dds <- ddsFull[ , colData(ddsFull)$treatment %in% c("Control","DPN") & 
                  colData(ddsFull)$time == "48h" ]


## ----refactor------------------------------------------------------------
dds$patient <- factor(dds$patient)
dds$treatment <- factor(dds$treatment)


## ----relevel-------------------------------------------------------------
dds$treatment <- relevel( dds$treatment, "Control" )


## ----multifactorColData--------------------------------------------------
colData(dds)


## ----runDESeq, cache=TRUE------------------------------------------------
dds <- DESeq(dds)


## ----extractResults------------------------------------------------------
res <- results(dds)
res


## ----resCols-------------------------------------------------------------
mcols(res)


## ----checkBaseMean-------------------------------------------------------
all.equal(res$baseMean, rowMeans(counts(dds)))
all.equal(res$baseMean, rowMeans(counts(dds,normalized=TRUE)))


## ----rawpvalue-----------------------------------------------------------
sum( res$pvalue < 0.05, na.rm=TRUE )
table( is.na(res$pvalue) )


## ----adjpvalue-----------------------------------------------------------
sum( res$padj < 0.1, na.rm=TRUE )


## ----strongGenes---------------------------------------------------------
resSig <- res[ which(res$padj < 0.1 ), ]
head( resSig[ order( resSig$log2FoldChange ), ] )


## ----strongGenesUp-------------------------------------------------------
tail( resSig[ order( resSig$log2FoldChange ), ] )


## ----proportionSigGenes--------------------------------------------------
table(sign(resSig$log2FoldChange))


## ----plotMA1,fig.cap="The MA-plot shows the $\\log_2$ fold changes from the treatment over the mean of normalized counts, i.e. the average of counts normalized by size factor. The \\Rpackage{DESeq2} package incorporates a prior on $\\log_2$ fold changes, resulting in moderated estimates from genes with low counts and highly variable counts, as can be seen by the narrowing of spread of points on the left side of the plot."----
plotMA(dds, ylim = c( -1.5, 1.5 ) )


## ----dispPlot,fig.cap='Plot of dispersion estimates.  See text for details.'----
plotDispEsts( dds )


## ----plotMApadjchange,fig.cap="The MA-plot with red points indicating adjusted p value less than 0.5."----
plotMA(dds, pvalCutoff = 0.5, ylim = c( -1.5, 1.5) )


## ----pvalHist,fig.cap="Histogram of the p values returned by the test for differential expression."----
hist( res$pvalue, breaks=100     )


## ----filter--------------------------------------------------------------
filterThreshold <- 2.0
keep <- rowMeans( counts( dds, normalized=TRUE ) ) > filterThreshold
table( keep )


## ----adjFilter-----------------------------------------------------------
min( res$padj[!keep], na.rm=TRUE )


## ----padjFilterCmp-------------------------------------------------------
table( p.adjust( res$pvalue, method="BH" ) < .1 )
table( p.adjust( res$pvalue[keep], method="BH" ) < .1 )


## ----pvalHistFilt,fig.cap="Histogram of the p values returned by the test for differential expression."----
hist( res$pvalue[keep], breaks=100 )


## ----loadOrg-------------------------------------------------------------
library( "org.Hs.eg.db" )


## ----keyType-------------------------------------------------------------
cols(org.Hs.eg.db)


## ----convertIDs----------------------------------------------------------
convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
   stopifnot( inherits( db, "AnnotationDb" ) )
   ifMultiple <- match.arg( ifMultiple )
   suppressWarnings( selRes <- AnnotationDbi::select( 
      db, keys=ids, keytype=fromKey, cols=c(fromKey,toKey) ) )
   if( ifMultiple == "putNA" ) {
      duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
      selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
   return( selRes[ match( ids, selRes[,1] ), 2 ] )
}


## ----addSymbols----------------------------------------------------------
res$symbol <- convertIDs( row.names(res), "ENSEMBL", "SYMBOL", org.Hs.eg.db )
res


## ----writeCSV------------------------------------------------------------
write.csv( as.data.frame(res), file="results.csv" )


## ----loadReactome--------------------------------------------------------
library( "reactome.db" )


## ----addEntrez-----------------------------------------------------------
res$entrez <- convertIDs( row.names(res), "ENSEMBL", "ENTREZID", org.Hs.eg.db )


## ----subsetTores2--------------------------------------------------------
res2 <- res[ res$entrez %in% keys( reactome.db, "ENTREZID" ) & 
   !is.na( res$pvalue) , ]
head(res2)


## ----queryReactome-------------------------------------------------------
reactomeTable <- AnnotationDbi::select( reactome.db, keys=res2$entrez, 
   keytype="ENTREZID", columns=c("ENTREZID","REACTOMEID") )
head(reactomeTable)


## ----CreateIncidenceMatrix-----------------------------------------------
incm <- do.call( rbind, with(reactomeTable, tapply( 
   ENTREZID, factor(REACTOMEID), function(x) res2$entrez %in% x ) ))
colnames(incm) <- res2$entrez

str(incm)


## ----PruneIncm-----------------------------------------------------------
incm <- incm[ rowSums(incm) >= 5, ]


## ----testFun-------------------------------------------------------------
testCategory <- function( reactomeID ) {
  isMember <- incm[ reactomeID, ]
  data.frame( 
     reactomeID = reactomeID,
     numGenes = sum( isMember ),
     avgLFC = mean( res2$log2FoldChange[isMember] ),
     strength = sum( res2$log2FoldChange[isMember] ) / sqrt(sum(isMember)),
     pvalue = t.test( res2$log2FoldChange[ isMember ] )$p.value,
     reactomeName = reactomePATHID2NAME[[reactomeID]] ) }


## ----testTestFun---------------------------------------------------------
testCategory("109581")


## ----runGSEA,cache=TRUE--------------------------------------------------
reactomeResult <- do.call( rbind, lapply( rownames(incm), testCategory ) )


## ----padjGSEA------------------------------------------------------------
reactomeResult$padjust <- p.adjust( reactomeResult$pvalue, "BH" )


## ----GSEAresult----------------------------------------------------------
reactomeResultSignif <- reactomeResult[ reactomeResult$padjust < 0.05, ]
reactomeResultSignif[ order(reactomeResultSignif$strength), ]


## ----topDEGene-----------------------------------------------------------
deGeneID <- "ENSG00000099194"
res[deGeneID,]
deGene <- range(rowData(dds[deGeneID,])[[1]])
names(deGene) <- deGeneID
deGene


## ----defineRange---------------------------------------------------------
as.character(seqnames(deGene))
ucscChrom <- paste0("chr",as.character(seqnames(deGene)))
ucscRanges <- ranges(flank(deGene,width=10e6,both=TRUE))
subsetRange <- GRanges(ucscChrom, ucscRanges)
subsetRange


## ----downloadERBS, eval=FALSE--------------------------------------------
## ##
## ## Please do not run this code if you do not have an internet connection,
## ## alternatively use the local file import in the next code chunk.
## ##
## library( "rtracklayer" )
## trackName <- "wgEncodeRegTfbsClusteredV2"
## tableName <- "wgEncodeRegTfbsClusteredV2"
## trFactor <- "ERalpha_a"
## mySession <- browserSession()
## ucscTable <- getTable(ucscTableQuery(mySession, track=trackName,
##                                      range=subsetRange, table=tableName,
##                                      name=trFactor))


## ----localERBS-----------------------------------------------------------
ucscTableFile <- system.file('extdata/localUcscTable.csv.gz',package='BiocBrazil2014')
ucscTable <- read.csv(gzfile(ucscTableFile), stringsAsFactors=FALSE)


## ----ERBS2Peaks----------------------------------------------------------
peaks <- with(ucscTable, GRanges(chrom, IRanges(chromStart, chromEnd), 
                                 score=score))
seqlevels(peaks) <- gsub("chr(.+)","\\1",seqlevels(peaks))
seqlevels(peaks) <- seqlevels(deGene)


## ----nearestPeakDo, cache=TRUE, echo=FALSE-------------------------------
# suppress warning which exists for BioC 2.12 of changes to distance()
suppressWarnings({d2nearest <- distanceToNearest(deGene, peaks)})


## ----nearestPeakShow, eval=FALSE-----------------------------------------
## d2nearest <- distanceToNearest(deGene, peaks)


## ----distanceShow, eval=FALSE--------------------------------------------
## distance(deGene, peaks)


## ----distanceDo, echo=FALSE----------------------------------------------
suppressWarnings({distance(deGene, peaks)})


## ----distanceToNearest---------------------------------------------------
d2nearest


## ----nearestPeakOut------------------------------------------------------
deGene
peaks[subjectHits(d2nearest)]


## ----plotPeaksAndGene, fig.width=6, fig.height=4, fig.cap="A 2 Mb genomic range showing the location of the differentially expressedgene (labelled 'g'), and the peaks (labelled 'p'). "----
plotRange <- start(deGene) + 1e6 * c(-1,1)
peakNearest <- ( seq_along(peaks) == subjectHits(d2nearest) )
plot(x=start(peaks), y=ifelse(peakNearest,.3,.2),
     ylim=c(0,1), xlim=plotRange, pch='p',
     col=ifelse(peakNearest,"red","grey60"),
     yaxt="n", ylab="",
     xlab=paste("2 Mb on chromosome",as.character(seqnames(deGene))))
points(x=start(deGene),y=.8,pch='g')


## ----calculateSumForCaption, echo=FALSE----------------------------------
plotGRange <- GRanges(seqnames(peaks[1]),IRanges(plotRange[1],plotRange[2]))
numPeaksInPlotRange <- sum(peaks %over% plotGRange)


## ----peakDists-----------------------------------------------------------
  peakDists <- diff(sort(start(peaks)))
  summary(peakDists)
  mean(peakDists)


## ----rld, cache=TRUE-----------------------------------------------------
rld <- rlogTransformation(dds)
head( assay(rld) )


## ----rldPlot,fig.width=10,fig.height=5,fig.cap="Scatter plot of sample 2 versus sample 1. Left: using an ordinary $\\log_2$ transformation. Right: Using the rlog transformation."----
par( mfrow = c( 1, 2 ) )
plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )


## ----euclDist------------------------------------------------------------
sampleDists <- dist( t( assay(rld) ) )
sampleDists


## ----sampleDistHeatmap,fig.cap="Heatmap of Euclidean sample distances after rlog transformation."----
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( colData(rld)$treatment, 
   colData(rld)$patient, sep="-" )
   
library( "gplots" )   
heatmap.2( sampleDistMatrix, trace="none" )


## ----betterheatmap,fig.cap="Using RColorBrewer to improve our heatmap"----
library("RColorBrewer")  
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)


## ----samplePCA,fig.cap="The same heatmap as in Figure \\ref{sampleDistHeatmap} but with better colours"----
print( plotPCA( rld, intgroup = c( "patient", "treatment") ) )


## ----topVarGenes---------------------------------------------------------
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )


## ----geneHeatmap, fig.width=9, fig.height=9, fig.cap="Heatmap with gene clustering."----
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
     trace="none", dendrogram="column", 
     col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))


## ----sessInfo, echo=FALSE------------------------------------------------
sessionInfo()


## ----relocateAnswers, echo=FALSE-----------------------------------------
relocateAnswers = function(con="DESeq2_parathyroid.tex"){
  txt = readLines(con=con)
  i1 = grep("\\begin{answer}", txt, fixed=TRUE)
  i2 = grep("\\end{answer}",  txt, fixed=TRUE)
  if(!((length(i1)==length(i2))&&all(i2>i1))) stop("\\begin{answer}/\\end{answer} macros are unbalanced.")
  i = unlist(mapply(`:`, i1, i2))
  res = txt[-i]
  sol = txt[i]
  dest = grep("\\section{Solutions}", res, fixed=TRUE)
  if(length(dest)!=1) stop("\\section{Solutions} not found exactly once.")
  res = c(res[1:dest], sol, res[(dest+1):length(res)])
  writeLines(res, con=con)
}


