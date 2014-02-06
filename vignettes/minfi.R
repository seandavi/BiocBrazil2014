
## ----options,echo=FALSE,results='hide'-----------------------------------
library(knitr)
opts_chunk$set(fig.show='asis',size='footnotesize',dev=c('pdf','png'),fig.align='center',fig.pos='!htbp',cache=TRUE)
options(width=70)


## ----load,results='hide',messages=FALSE----------------------------------
require(minfi)
require(minfiData)


## ----baseDir-------------------------------------------------------------
baseDir <- system.file("extdata", package = "minfiData")
list.files(baseDir)


## ----baseDir2------------------------------------------------------------
list.files(file.path(baseDir, "5723646052"))


## ----sampleSheet---------------------------------------------------------
targets <- read.450k.sheet(baseDir)
targets


## ----BasenameColumn------------------------------------------------------
sub(baseDir, "", targets$Basename)


## ----paths---------------------------------------------------------------
RGset <- read.450k.exp(base = baseDir, targets = targets)


## ----pData---------------------------------------------------------------
show(RGset)
pd <- pData(RGset)
pd


## ----read2---------------------------------------------------------------
RGset2 = read.450k.exp(file.path(baseDir, "5723646052"))
RGset3 = read.450k.exp(baseDir, recursive = TRUE)


## ----detp----------------------------------------------------------------
detP <- detectionP(RGset)
failed <- detP > 0.01
head(failed, n=3)


## ----colmeans1,cache=FALSE-----------------------------------------------
colMeans(failed)


## ----rowmeans1,cache=FALSE-----------------------------------------------
sum(rowMeans(failed)>0.5)


## ------------------------------------------------------------------------
MSet <- preprocessRaw(RGset)
MSet


## ------------------------------------------------------------------------
head(getMeth(MSet),n=3)
head(getUnmeth(MSet),n=3)


## ----qc1,fig.cap='A QC plot can show outlier samples that may be of poor quality'----
qc <- getQC(MSet)
head(qc)
plotQC(qc)


## ----qcReport-quick,eval=FALSE-------------------------------------------
## qcReport(RGset, sampNames = pd$Sample_Name,
##          sampGroups = pd$Sample_Group, pdf = "qcReport.pdf")


## ----densityPlot,fig.cap="Beta density plots"----------------------------
densityPlot(RGset, sampGroups = pd$Sample_Group, 
            main = "Beta", xlab = "Beta")


## ----densityBeanPlot,fig.cap="Bean plot of beta values"------------------
par(oma=c(2,10,1,1))
densityBeanPlot(RGset, sampGroups = pd$Sample_Group, 
                sampNames = pd$Sample_Name)


## ----controlStripPlot,fig.cap="Beta stripplot"---------------------------
controlStripPlot(RGset, controls="BISULFITE CONVERSION II", 
                 sampNames = pd$Sample_Name)


## ----Msetraw-------------------------------------------------------------
MSet.raw <- preprocessRaw(RGset)


## ----allMsets------------------------------------------------------------
MSet.norm <- preprocessIllumina(RGset, bg.correct = TRUE,
                                 normalize = "controls", reference = 2)


## ----MSet----------------------------------------------------------------
getMeth(MSet.raw)[1:4,1:3]
getUnmeth(MSet.raw)[1:4,1:3]
getBeta(MSet.raw, type = "Illumina")[1:4,1:3]
getM(MSet.raw)[1:4,1:3]


## ----MDS,fig.cap="Multi-dimensional scaling plot"------------------------
mdsPlot(MSet.norm, numPositions = 1000, sampGroups = pd$Sample_Group, 
	sampNames = pd$Sample_Name)


## ----preprocessSwan------------------------------------------------------
Mset.swan <- preprocessSWAN(RGsetEx, MsetEx)


## ----plotBetaTypes,fig.width=8,fig.height=4,fig.cap="The effect of normalizing using SWAN."----
par(mfrow=c(1,2))
plotBetasByType(MsetEx[,1], main = "Raw")
plotBetasByType(Mset.swan[,1], main = "SWAN")


## ----subset-mset---------------------------------------------------------
mset <- MSet.norm[1:20000,]


## ----dmpFinder-categorical-----------------------------------------------
table(pd$Sample_Group)
M <- getM(mset, type = "beta", betaThreshold = 0.001)
dmp <- dmpFinder(M, pheno=pd$Sample_Group, type="categorical")
head(dmp)


## ----plot-dmps-categorical-----------------------------------------------
cpgs <- rownames(dmp)[1:4]
par(mfrow=c(2,2))
plotCpg(mset, cpg=cpgs, pheno=pd$Sample_Group)


## ----set-seed, echo=FALSE------------------------------------------------
set.seed(123)


## ----sim-pheno-----------------------------------------------------------
continuousPheno <- rnorm(nrow(pd))


## ----dmpFinder-continuous------------------------------------------------
dmp <- dmpFinder(mset, pheno=continuousPheno, type="continuous")
dmp[1:3,]


## ----filter-dmp----------------------------------------------------------
dmp <- subset(dmp, abs(beta)>1)


## ----plot-dmps-continuous------------------------------------------------
cpgs <- rownames(dmp)[1:4]
par(mfrow=c(2,2))
plotCpg(mset, cpg=cpgs, type="continuous",
        pheno=continuousPheno, xlab="Phenotype 1")


## ----sessionInfo,results='asis',echo=FALSE-------------------------------
toLatex(sessionInfo())


