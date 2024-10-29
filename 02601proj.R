library(WGCNA)          ###used for topological overlap calculation and clustering steps
library(RColorBrewer)   ###used to create nicer colour palettes
library(preprocessCore) ###used by the quantile normalization function

#Note: the data can be downloaded from the Gene Expression Omnibus
# http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS2901
data<-as.matrix(read.csv(file="GDS2901.soft",skip=166,row.names=1,sep="\t",header=T))
data<-data[-15924,]
rawData<-matrix(as.numeric(data[,-1]),nrow=15923)
dimnames(rawData)<-dimnames(data[,-1])
#we create an annotation matrix containing the matches between probesets and gene names
anno<-as.matrix(data[-2475,1]) 
normData<-normalize.quantiles(log2(rawData))
dimnames(normData)<-dimnames(rawData)

#we remove the probeset at index 2475 because
#after quantile normalization it has zero variance
#(the probeset has the highest signal of all samples)
normData<-normData[-2475,]  

datC1<-t(normData[,c(1:12,25:36,37:48)]) ### these samples correspond to the Eker mutants.
# Note that since the Eker mutants have two sets of 12 control samples (13:24 and 37:48)
# we discard one to have a symmetric perturbation (carcinogenic vs control) between the two conditions (Eker mutants vs wild-types)
datC2<-t(normData[,49:84]) ###those samples correspond to the wild-types