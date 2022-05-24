library(cBioPortalData)
library(survival)
get_GBM_data = function(d=1){
  if (!d%in%c(1,2,3)){
    stop('degree of mechanistic model should be 1, 2 or 3')
  }
  mydat = cBioDataPack("gbm_tcga_pub")
  #preprocess data from cbioportal
  #caveat: the data are updated irregularly, modification might be needed
  #get expression
  expression_med = mydat[[4L]]
  exp_med = assay(expression_med)
  exp_med = exp_med[complete.cases(exp_med),]
  #get miRNA profile
  expression_miRNA = mydat[[6L]]
  exp_miRNA = assay(expression_miRNA)
  exp_miRNA = exp_miRNA[complete.cases(exp_miRNA),]
  #get methylation profile
  expression_methy = mydat[[9L]]
  exp_methy = assay(expression_methy)
  exp_methy = exp_methy[complete.cases(exp_methy),-c(1,2)]
  #expression_med@NAMES
  (cbio <- cBioPortal())
  mols <- molecularProfiles(cbio, "gbm_tcga_pub")
  #mols[["molecularProfileId"]]
  #get clinical data
  clinicalDat = clinicalData(cbio,"gbm_tcga_pub")
  clinicalDat2013 = clinicalData(cbio,"gbm_tcga_pub2013")
  clinicalDat$AGE = clinicalDat2013[clinicalDat2013$patientId%in%clinicalDat$patientId,]$AGE
  clinicalDF = data.frame(clinicalDat[,c("AGE","OS_MONTHS","OS_STATUS","PRIOR_GLIOMA",
                                         "SEX","KARNOFSKY_PERFORMANCE_SCORE","PRETREATMENT_HISTORY")])
  clinicalDF$AGE = as.vector(clinicalDat$AGE)
  rownames(clinicalDF) = clinicalDat$patientId
  clinicalIDs = clinicalDat$patientId[complete.cases(clinicalDF)]
  myIDs = intersect(paste0(clinicalIDs,'-01'),colnames(exp_med))
  finalIDs = intersect(myIDs, colnames(exp_methy))
  myexp_med = exp_med[,finalIDs]
  #myexp_miRNA = exp_miRNA[,finalIDs]
  myexp_methy = exp_methy[,finalIDs]
  myclinical = clinicalDF[paste0(rownames(clinicalDF),'-01')%in%finalIDs,]
  myclinical$AGE = as.numeric(myclinical$AGE)
  myclinical$OS_MONTHS = as.numeric(myclinical$OS_MONTHS)
  myclinical$OS_STATUS = ifelse(myclinical$OS_STATUS=="1:DECEASED",1,0)
  myclinical$PRIOR_GLIOMA = ifelse(myclinical$PRIOR_GLIOMA=="NO",0,1)
  myclinical$SEX = ifelse(myclinical$SEX=="Male",0,1)
  myclinical$KARNOFSKY_PERFORMANCE_SCORE = as.numeric(myclinical$KARNOFSKY_PERFORMANCE_SCORE)
  myclinical$PRETREATMENT_HISTORY = ifelse(myclinical$PRETREATMENT_HISTORY=="NO",0,1)
  Y=myclinical$OS_MONTHS
  Delta = myclinical$OS_STATUS
  #use AFT model to select 1000 genes having strongest association with survival
  p.surv = matrix(nrow=nrow(myexp_med),ncol=1)
  for (i in 1:nrow(myexp_med)){
    survregWeibull <- survreg(Surv(OS_MONTHS, OS_STATUS) ~ AGE+myexp_med[i,],
                              myclinical, dist = "weibull")
    s = summary(survregWeibull)
    p.surv[i] = s$table[3,4]
  }
  G=t(myexp_med[order(p.surv)[1:1000],])
  C = as.matrix(myclinical[,c('AGE','PRIOR_GLIOMA','SEX','PRETREATMENT_HISTORY')])
  M.raw=t(myexp_methy[colnames(G)[colnames(G)%in%rownames(myexp_methy)],])
  N=nrow(M.raw)
  if (d==3){
    M = matrix(0,nrow=N,ncol = 3*ncol(M.raw))
    for (j in 1:ncol(M.raw)) {
      M[,(3*(j-1)+1):(3*j)]=ns(M.raw[,j],df=3)
    }
  }else if (d==2){
    M = matrix(0,nrow=N,ncol = 2*ncol(M.raw))
    for (j in 1:ncol(M.raw)) {
      M[,(2*(j-1)+1)]=M.raw[,j]
      M[,(2*j)]=M.raw[,j]^2
    }
  }else{
    M=M.raw}
  return(list(G=G,C=C,M=M,Y=Y,Delta=Delta))
}

get_grouping = function(full_gene_list,file_dir="functional_classification.txt"){
  K0 = length(full_gene_list)
  functionalClassification <- read.delim(file_dir, header=FALSE, comment.char="#") 
  breakPoint = which(startsWith(as.character(functionalClassification[,1]),'Gene Group'))
  gene_names = unique(functionalClassification[-c(breakPoint,breakPoint+1),1])
  grouping = matrix(0,length(gene_names),length(breakPoint))
  rownames(grouping) = gene_names
  for (i in 1:(length(breakPoint)-1)){
    genes = functionalClassification[(breakPoint[i]+2):(breakPoint[i+1]-1),1]
    grouping[as.character(genes),i]=1
  }
  genes = functionalClassification[(tail(breakPoint,1)+2):nrow(functionalClassification),1]
  grouping[as.character(genes),i+1] = 1
  R = ncol(grouping)+1
  gene_group = rep(R,K0)
  for (k in 1:K0){
    genename = full_gene_list[k]
    if (genename %in% rownames(grouping)){
      gene_group[k] = which(grouping[genename,]>0)[1]
    }
  }
  names(gene_group) = full_gene_list
  return(gene_group)
}

GBM_data = get_GBM_data(d=2)
gene_group = get_grouping(colnames(GBM_data$G), "functional_classification.txt")
save(GBM_data, gene_group, file='GBM_data2.RData')
