# Setting workpath
setwd("D:\\Ziyi\\School\\PMO\\metanalysis\\Simulation\\combine")

library(stringi)
library(GEOquery)
library(simpleaffy)
library(RColorBrewer)
library(affyPLM)
library(R.utils)
library(affy)
library(annotate)

# download dataset from GEO dataset, download "affy_geo" format data to the "data" fold, if download "affy_series_matrix"
# data directly from GEO dataset 

#getGEOSuppFiles("GSE10072")
#untar("GSE10072/GSE10072_RAW.tar", exdir="data/GSE10072")
#cels <- list.files("data/GSE10072", pattern = "[gz]")
#sapply(paste("data/GSE10072", cels, sep="/"), gunzip)

#getGEOSuppFiles("GSE37745")
#untar("GSE37745/GSE37745_RAW.tar", exdir="data/GSE37745")
#cels <- list.files("data/GSE37745", pattern = "[gz]")
#sapply(paste("data/GSE37745", cels, sep="/"), gunzip)
#celfiles <- read.affy(covdesc="phenodata.txt", path="data")

# install all packages
source('metaAnalyze.R')
source('metaAnalyzePamr.R')

# Check that the corresponding platformInfo is in the list of supported platforms.
#getSupportedPlatforms() 

## according to the research to modify parameters
metaAnalysisName = 'luca_main' # name of file of processed expression data and prefix for output files of meta-analysis
studyMetadataFileName = 'luca_main_study_metadata.csv' # table of study metadata
sampleMetadataFileName = 'luca_sample_metadata.csv' # table of sample metadata
parentFolderPath = 'data' # path to folder that contains the expression data
denovo = TRUE # process the raw data de novo or load processed data
className = 'class' # column name in table of sample metadata that contains data to predict
classesTrain = c('Tumor','Normal') # values for multinomial classification


# load the table of study metadata
studyMetadata = read.csv(studyMetadataFileName, header=TRUE, stringsAsFactors=FALSE)
rownames(studyMetadata) = studyMetadata[,'study']
studyMetadata[,'discovery'] = studyMetadata[,'discovery']==1
studyMetadata[,'validation'] = studyMetadata[,'validation']==1

# load the table of sample metadata
sampleMetadataTmp = read.csv(sampleMetadataFileName, header=TRUE, stringsAsFactors=FALSE)
rownames(sampleMetadataTmp) = sampleMetadataTmp[,'sample']
sampleMetadata = sampleMetadataTmp[sampleMetadataTmp[,'study'] %in% studyMetadata[,'study'],]

# load the expression data for all studies
if (denovo) {
  esetList = getStudyDataList(parentFolderPath, studyMetadata) ############
  
  saveRDS(esetList, file=paste0(metaAnalysisName, '.rds'))
} else {
  esetList = readRDS(paste0(metaAnalysisName, '.rds'))}

# select samples that are in sampleMetadata, extract the expression matrix
ematList = cleanStudyData(esetList, sampleMetadata)

# merge expression data and perform cross-study normalization of discovery datasets
ematMergedDiscoveryAllClasses = mergeStudyData(ematList[studyMetadata[studyMetadata[,'discovery'], 'study']],
                                               sampleMetadata, covariateName=NA)
geneIds = Reduce(intersect, lapply(ematList, function(x) rownames(x)))
geneSymbols = getSYMBOL(geneIds, 'org.Hs.eg')

# save the expression data and gene name
write.table(ematMergedDiscoveryAllClasses,file="luca_alldata.txt",sep="\t")
write.table(geneSymbols,file="luca_gene_symbols.txt",sep="\t")
