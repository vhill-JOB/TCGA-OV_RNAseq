library(TCGAbiolinks)
library(dplyr)
library(DT)
#query GDC Data portal
query2 <- GDCquery(project = "TCGA-OV", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  file.type = "htseq.counts")
#retrieve barcodes for sample identification
samplesDown <- getResults(query2,cols=c("cases"))
#Download barcodes
GDCdownload(query = query2,directory = "RES_DEA/")
#Prepare expression matrice
dataPrep <- GDCprepare(query = query2, save = TRUE, directory =  "RES_DEA/",save.filename = "RES_DEA/counts.rda", summarizedExperiment = TRUE)
#Gene expression and visualization
dataPrepro <- TCGAanalyze_Preprocessing(object = dataPrep,cor.cut = 0.6,datatype = "HTSeq - Counts")
# normalization of genes. Other method: geneLength
dataNormalized <- TCGAanalyze_Normalization(tabDF = dataPrepro,geneInfo=TCGAbiolinks::geneInfoHT,method="gcContent")
# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNormalized,method = "quantile",qnt.cut =  0.25)
#see size of data
dim(dataFilt)
# selection of normal samples "NT"
samplesNormalTissues <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),typesample = c("NT"))
# selection of tumor samples "TP"
samplesPrimaryTumors <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),typesample = c("TP"))
# selection of normal samples "TR"
samplesRecurrentTumors <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),typesample = c("TR"))
dataFiltTP <- dataFilt[,samplesPrimaryTumors]
dataFiltTR <- dataFilt[,samplesRecurrentTumors]
dataFiltNT <- dataFilt[,samplesNormalTissues]
#Differential Expression Analysis to characterize Primary tumor
dataDEGsPrimaryTumors <- TCGAanalyze_DEA(mat1 = dataFiltNT,
                            mat2 = dataFiltTP,
                            Cond1type = "Normal",
                            Cond2type = "Tumor",fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")
                            
#Differential Expression Analysis to characterize Recurrent tumor                            
dataDEGsRecurrentTumors <- TCGAanalyze_DEA(mat1 = dataFiltNT,
                            mat2 = dataFiltTR,
                            Cond1type = "Normal",
                            Cond2type = "Tumor",fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")


