library(gridExtra)
library(ggplot2)
library(glmnet)
library(MASS)
library(survival)
library(MASS)
library(coda)
library(caret)
library(truncnorm)
# library(clusterProfiler)
# library(org.Hs.eg.db)
set.seed(2020)

#EMVS for the first stage
#setting up known parameters:
get_GBM_data = function(mechanistic_mod='linear'){
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
  mydict <- read.delim("~/Dropbox/harvard/emvs/data/gbm_tcga_pub/data_methylation_hm27.txt")
  id2sym = as.list(mydict[,1])
  names(id2sym) <- mydict[,2]
  # sym2id = as.list(mydict$Entrez_Gene_Id)
  # names(sym2id) = mydict$Hugo_Symbol
  exp_methy = assay(expression_methy)
  #rownames(exp_methy) = unlist(id2sym[rownames(exp_methy)],use.names = FALSE)
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
  if (mechanistic_mod=='cubic'){
    M = matrix(0,nrow=N,ncol = 3*ncol(M.raw))
    for (j in 1:ncol(M.raw)) {
      M[,(3*(j-1)+1):(3*j)]=ns(M.raw[,j],df=3)
    }
  }else if (mechanistic_mod=='quadratic'){
    M = matrix(0,nrow=N,ncol = 2*ncol(M.raw))
    for (j in 1:ncol(M.raw)) {
      M[,(2*(j-1)+1)]=M.raw[,j]
      M[,(2*j)]=M.raw[,j]^2
    }
  }else{
      M=M.raw}
  return(list(G=G,C=C,M=M,Y=Y,Delta=Delta))
}

get_grouping = function(){
  #subdata
  #functionalClassification <- read.delim("~/Dropbox/harvard/emvs/code/results/clustering/functionalClassification.txt", header=FALSE, comment.char="#")
  #full data
  functionalClassification <- read.delim("~/Dropbox/harvard/emvs/code/results/clustering/functional_classification.txt", header=FALSE, comment.char="#") 
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
    genename = colnames(G)[k]
    if (genename %in% rownames(grouping)){
      gene_group[k] = which(grouping[genename,]>0)[1]
    }
  }
  return(gene_group)
  
  # sym2id = mapIds(org.Hs.eg.db,colnames(G),'ENTREZID','SYMBOL')
  # kk <- enrichKEGG(gene         = sym2id,
  #                  organism     = 'hsa',
  #                  pvalueCutoff = 0.05)
  # seen_gene = list()
  # gene2pw = vector(mode = "list", length = K0)
  # names(gene2pw) = sym2id[colnames(G[,1:K0])]
  # for (i in 1:nrow(kk@result)){
  #   pwid = rownames(kk@result)[i]
  #   genes = strsplit(kk@result[pwid,'geneID'], split = "/")[[1]]
  #   for (gene in genes){
  #     if (!gene %in% seen_gene){
  #       gene2pw[gene] = pwid
  #       seen_gene = append(seen_gene,gene)
  #     }
  #   }
  # }
  # pwid2idx = vector(mode = "list", length = length(unique(unlist(gene2pw))))
  # names(pwid2idx) = unique(unlist(gene2pw))
  # idx=1
  # for (pw in unique(unlist(gene2pw))){
  #   pwid2idx[pw] = idx
  #   idx=idx+1
  # }
  # R = length(pwid2idx)+1
  # gene_group = rep(R,K0)
  # for (k in 1:K0){
  #   genename = colnames(G)[k]
  #   id = sym2id[genename]
  #   if (!is.null(unlist(gene2pw[unlist(id)]))){
  #     gene_group[k]=unlist(pwid2idx[unlist(gene2pw[unlist(id)])])
  #   }
  # }
  # return(gene_group)
}

EMVS = function(M,G,grouping,nu0=0.5,nu1=10^3,nu=1,lambda=1,a=1,b=1){
  #get dimension of data
  N=nrow(M)
  K=ncol(G)
  J=ncol(M)
  #get functional cluster information
  R = length(unique(grouping))
  Kr = list()
  for (r in unique(grouping)){
    Kr[[r]]=which(grouping==r)
  }
  #initialize parameters
  estalpha = numeric(R)
  esttau2 = numeric(R)
  estOme=matrix(0,nrow=J,ncol=K)
  estsig2=numeric(K)
  esttheta=numeric(K)
  iteration=numeric(K)
  #estimate random intercept of each functional cluster
  for (r in 1:R){
    print(r)
    I=100
    #initialize estimation
    ome = thresh_ome = array(0,dim=c(I,J,length(Kr[[r]])))
    sig2 = array(0,dim=c(I,length(Kr[[r]])))
    theta = array(0,dim=c(I,length(Kr[[r]])))
    alpha = array(0,I)
    tau2 = array(0,I)
    ome[1,,]=numeric(J)
    sig2[1,]=1
    theta[1,]=0.1
    alpha[1]=0.1
    tau2[1]=0.1
    #E-step
    for (i in 2:I){
      p=array(0,dim=c(J,length(Kr[[r]])))
      d=array(0,dim=c(J,length(Kr[[r]])))
      for (k in 1:length(Kr[[r]])){
        #initial value of parameters:
        for (j in 1:J){
          p1=dnorm(ome[(i-1),j,k],mean=0,sd=sqrt(sig2[i-1,k]*nu1))*theta[i-1,k]
          p0=dnorm(ome[(i-1),j,k],mean=0,sd=sqrt(sig2[i-1,k]*nu0))*(1-theta[i-1,k])
          p[j,k]=p1/(p1+p0)
          d[j,k]=(1-p[j,k])/nu0+p[j,k]/nu1
        }
        D=diag(d[,k])
        #M-step
        ome[i,,k]=coef(glmnet(x=M,y=(G[,Kr[[r]][k]]-alpha[i-1]), intercept = F, standardize = F,
                              lambda=(sum(d[,k])/ncol(M))/N,alpha=0,penalty.factor = sqrt(d[,k])))[-1,1]
        # ome[i,,k]=(solve(D)-solve(D)%*%t(M)%*%solve(diag(nrow=N,ncol=N)+M%*%solve(D)%*%t(M))%*%M%*%solve(D))%*%t(M)%*%(G[,Kr[[r]][k]]-alpha[i-1])
        sig2[i,k]=(t(G[,Kr[[r]][k]]-M%*%ome[i,,k]-alpha[i-1])%*%(G[,Kr[[r]][k]]-M%*%ome[i,,k]-alpha[i-1])+t(sqrt(D)%*%ome[i,,k])%*%(sqrt(D)%*%ome[i,,k])+nu*lambda)/(N+J+nu+2)
        theta[i,k]=(sum(p[,k])+a-1)/(a+b+J-2)
      }
      #thresholding
      model=matrix(0,nrow=J,ncol=length(Kr[[r]]))
      for (k in 1:length(Kr[[r]])){
        postp=numeric(J)
        for (j in 1:J){
          p1=dnorm(ome[i,j,k],mean=0,sd=sqrt(sig2[i,k]*nu1))*theta[i,k]
          p0=dnorm(ome[i,j,k],mean=0,sd=sqrt(sig2[i,k]*nu0))*(1-theta[i,k])
          postp[j]=p1/(p1+p0)
          model[j,k]=ifelse(postp[j]>=0.5,1,0)
          thresh_ome[i,j,k]=ifelse(model[j,k]==1,ome[i,j,k],0)
        }
      }
      tau2[i]=(alpha[i]^2+1)/4
      alpha[i]=((tau2[i]/sig2[i,])%*%colSums(G[,Kr[[r]]]-M%*%thresh_ome[i,,]))/(1+N*sum((tau2[i]/sig2[i,])))
      # print('alpha')
      print(alpha[i])
      # print('ome:')
      # print(sum(abs(thresh_ome[i,,])))
      # print('----------------')
      if (sum(abs(ome[i,,]-ome[(i-1),,]))+sum(abs(sig2[i,]-sig2[i-1,]))+sum(abs(theta[i,]-theta[i-1,]))+abs(alpha[i]-alpha[i-1])+abs(tau2[i]-tau2[i-1])<0.0001) {
        #print('converge')
        estOme[,Kr[[r]]]=thresh_ome[i,,] 
        estsig2[Kr[[r]]]=sig2[i,] 
        esttheta[Kr[[r]]]=theta[i,] 
        esttau2[r]=tau2[i]
        estalpha[r]=alpha[i]
        iteration[Kr[[r]]]=i
        #print(paste0('alpha',estalpha[r]))
        break}
    }
    estOme[,Kr[[r]]]=thresh_ome[i,,] 
    estsig2[Kr[[r]]]=sig2[i,] 
    esttheta[Kr[[r]]]=theta[i,] 
    esttau2[r]=tau2[i]
    estalpha[r]=alpha[i]
    iteration[Kr[[r]]]=i
  }
  #compute R2 for each gene site
  R2=numeric(K)
  for (kk in 1:K){
    ind=(estOme[,kk]!=0)
    R2[kk]=ifelse(sum(ind)==0,0,(summary(lm((G[,kk]-estalpha[grouping[kk]])~M[,ind]-1))$r.squared))
  }
  return(list(estOme=estOme,estsig2=estsig2,esttheta=esttheta,
              esttau2=esttau2,estalpha=estalpha,
              iteration=iteration,R2=R2))
}

NEG=function(G,Y,Zmatrix,mypath='/Users/xuehao/Dropbox/harvard/emvs/code/results/mixed_effect/',
             a=0.1,gstr=1,alpha=1,gam=1,s1=1,s2=2,c=1,d=1){
  #Define LOG Parabolic Cylinder Function (using "fAsianOptions" package with Whittaker W function representation)
  write.table(" ",paste0(mypath,"out1",a,"_",gstr,".txt"),append=FALSE,row.names=FALSE,col.names=FALSE)
  lpcf=function(k,v,z){
    vz=as.matrix(t(c(k,v,z)))
    write.table(vz,paste0(mypath,"zv",a,"_",gstr,".txt"),sep=',',append=FALSE,row.names=FALSE,col.names=FALSE)
    #evoke lpcf.py for sampling from LOG Parabolic Cylinder Function
    system(paste0("python3 ", mypath, "lpcf.py", " ", mypath, " ", a," ",gstr))
    A=read.table(paste0(mypath,"out1",a,"_",gstr,".txt"))
    return(A)
  }
  
  #E-M loop
  #hyperparameter
  K=ncol(G)
  N=nrow(Y)
  I=100
  #Store the values of expectations in E-step
  Elambdaj =as.vector(rep(0,K))
  Elambdaj2=as.vector(rep(0,K)) 
  nz=length(Zmatrix[1,])
  #initial values
  sigmak=rep(1,I)
  betak=matrix(1,nrow=I,ncol=K)
  bbk=matrix(1,nrow=I,ncol=nz)
  v1=as.vector(rep((-(2*a+1+s1)),K))
  v2=as.vector(rep((-(2*a+1+s2)),K))
  vb=as.vector(rep(-2*(a+0.5),K))
  if (gstr=='scale'){
    g=N^(-2)
  }else if (is.numeric(as.numeric(gstr))){
    #g=1
    g=as.numeric(as.numeric(gstr))
  }else{
    print('g is supposed to be a number or string indicating its value, not a valide g detected, using defalt value g=1')
  }
  #starting iteration
  for(kkk in 2:I){
    #E-step
    hf=1/((Zmatrix)%*%bbk[kkk-1,])
    ZZ=as.vector(abs(betak[kkk-1,])*sqrt(hf)/(sqrt(g)*sigmak[kkk-1]))
    lPCFb=lpcf(K,vb,ZZ)
    lPCF1=lpcf(K,v1,ZZ)
    lPCF2=lpcf(K,v2,ZZ)
    logmar=as.vector(rep(0,K))
    for (h in 1:K){
      logmar[h]=as.numeric(log(a*(2^a)*sqrt(hf[h])/(sqrt(pi)*sqrt(g)*sigmak[kkk-1]))+lgamma(a+0.5)+(betak[kkk-1,h])^2*hf[h]/(4*(g)*sigmak[kkk-1]^2)+lPCFb[h])
    }
    for(h in 1:K){
      Elambdaj[h]=as.numeric(exp(lgamma(2*a+s1+1)+(0.5*s1+0.5)*log(hf[h])-log(sqrt(g)*sigmak[kkk-1]*gamma(a)*2^(a+0.5*s1))+betak[kkk-1,h]^2*(hf[h])/(4*(g)*sigmak[kkk-1]^2)-logmar[h]+lPCF1[h]))
      Elambdaj2[h]=as.numeric(exp(lgamma(2*a+s2+1)+(0.5*s2+0.5)*log(hf[h])-log(sqrt(g)*sigmak[kkk-1]*gamma(a)*2^(a+0.5*s2))+betak[kkk-1,h]^2*(hf[h])/(4*(g)*sigmak[kkk-1]^2)-logmar[h]+lPCF2[h]))
    }
    #M-step
    #update Beta through Adaptive LASSO, Zou (2006) (using "lqa" package)
    Bayelsso = glmnet(x=G,y=Y,alpha=1,intercept = F,
                      family="gaussian", lambda=2*sqrt(2)/sqrt(g)*sigmak[kkk-1]*(sum(Elambdaj)/N)/(2*N),
                      penalty.factor=as.vector(Elambdaj)/(sum(Elambdaj)/N))
    betak[kkk,]=coef(Bayelsso)[-1]
    # Bayelsso=lqa.default(x=G,y=as.vector(Y),intercept=FALSE,
    #                      family=gaussian(),
    #                      penalty=adaptive.lasso(lambda=2*sqrt(2)/sqrt(g)*sigmak[kkk-1],
    #                                             al.weights=as.vector(Elambdaj)))
    # betak[kkk,]=Bayelsso$coefficients
    plot(betak[kkk,],main=kkk)
    #update sigma (correct the typo on denominator in the paper)
    sigmak[kkk]=as.numeric(((sqrt(2)*sum(abs(Elambdaj/sqrt(g)*betak[kkk,]))+sqrt(2*(sum(abs(Elambdaj/sqrt(g)*betak[kkk,])))^2+4*((t(Y-G%*%betak[kkk,])%*%(Y-G%*%betak[kkk,]))^1+2*d)*(N+K+2*c+2)))/(2*(N+K+2*c+2))))
    #update b
    Q2=function(bk){
      sum(-a*log(1/(t(bk)%*%t(Zmatrix)))-Elambdaj2*(t(bk)%*%t(Zmatrix)))+sum((alpha-1)*log(bk)-gam*bk)
    }
    bbk[kkk,]=as.vector(optim(rep(1,nz),Q2,lower=rep(1e-7,nz),
                              upper=rep(1000,nz),method="L-BFGS-B",
                              gr = NULL,control=list(fnscale=-1))$par)
    #convergence criterion  
    betat=betak[kkk,]-betak[kkk-1,]
    bt=bbk[kkk,]-bbk[kkk-1,]
    if(sum(abs(betat))+sum(abs(bt))<10^(-3)){break}
  }
  lst=list(b=bbk,beta=betak,k=kkk,a=a,g=g)
  return(lst)
}

NEG_fun=function(Y,G,C,a0,gstr,Zmatrix,mypath){
  a=a0
  #Define LOG Parabolic Cylinder Function (using "fAsianOptions" package with Whittaker W function representation)
  write.table(" ",paste0(mypath,"out1",a,"_",gstr,".txt"),append=FALSE,row.names=FALSE,col.names=FALSE)
  lpcf=function(k,v,z){
    vz=as.matrix(t(c(k,v,z)))
    write.table(vz,paste0(mypath,"zv",a,"_",gstr,".txt"),sep=',',append=FALSE,row.names=FALSE,col.names=FALSE)
    system(paste0("python3 ", mypath, "lpcf.py", " ", mypath, " ", a," ",gstr))
    A=read.table(paste0(mypath,"out1",a,"_",gstr,".txt"))
    return(A)
  }
  alpha=1
  gam=1
  s1=1
  s2=2
  c=1
  d=1
  N = nrow(G)
  K0 = ncol(G)
  L = ncol(C)
  K = K0+L
  I=300
  #Store the values of expectations in E-step
  Elambdaj =as.vector(rep(0,K0))
  Elambdaj2=as.vector(rep(0,K0)) 
  nz=length(Zmatrix[1,])
  #initial values
  sigmak=rep(1,I)
  betak=matrix(0,nrow=I,ncol=K)
  bbk=matrix(1,nrow=I,ncol=nz)
  v1=as.vector(rep((-(2*a+1+s1)),K0))
  v2=as.vector(rep((-(2*a+1+s2)),K0))
  vb=as.vector(rep(-2*(a+0.5),K0))
  #g=N^(-2)
  if (gstr=='scale'){
    g0=1/(N^2)
  }else{
    g0 = as.numeric(gstr)
  }
  g=g0
  #starting iteration
  for(kkk in 2:I){
    #E-step
    print(bbk[kkk-1,])
    hf=1/((Zmatrix)%*%bbk[kkk-1,])
    ZZ=as.vector(abs(betak[kkk-1,1:K0])*sqrt(hf)/(sqrt(g)*sigmak[kkk-1]))
    lPCFb=lpcf(K0,vb,ZZ)
    lPCF1=lpcf(K0,v1,ZZ)
    lPCF2=lpcf(K0,v2,ZZ)
    logmar=as.vector(rep(0,K0))
    for (h in 1:K0){
      logmar[h]=as.numeric(log(a*(2^a)*sqrt(hf[h])/(sqrt(pi)*sqrt(g)*sigmak[kkk-1]))+lgamma(a+0.5)+(betak[kkk-1,h])^2*hf[h]/(4*(g)*sigmak[kkk-1]^2)+lPCFb[h])
    }
    for(h in 1:K0){
      Elambdaj[h]=as.numeric(exp(lgamma(2*a+s1+1)+(0.5*s1+0.5)*log(hf[h])-log(sqrt(g)*sigmak[kkk-1]*gamma(a)*2^(a+0.5*s1))+betak[kkk-1,h]^2*(hf[h])/(4*(g)*sigmak[kkk-1]^2)-logmar[h]+lPCF1[h]))
      Elambdaj2[h]=as.numeric(exp(lgamma(2*a+s2+1)+(0.5*s2+0.5)*log(hf[h])-log(sqrt(g)*sigmak[kkk-1]*gamma(a)*2^(a+0.5*s2))+betak[kkk-1,h]^2*(hf[h])/(4*(g)*sigmak[kkk-1]^2)-logmar[h]+lPCF2[h]))
    }
    
    #M-step
    #update Beta through Adaptive LASSO, Zou (2006) (using "lqa" package)
    Bayelsso_gen = glmnet(x=as.matrix(G),
                          y=Y-as.vector(as.matrix(C)%*%betak[kkk-1,(K0+1):K]),
                      intercept=FALSE, alpha=1,
                      family="gaussian", lambda=2*sqrt(2)/sqrt(g)*sigmak[kkk-1]*(sum(Elambdaj[1:K0])/N)/(2*N),
                      penalty.factor =as.vector(Elambdaj[1:K0]))
    betak[kkk,1:K0]=coef(Bayelsso_gen)[-1]
    # Bayelsso_gen=lqa.default(x=as.matrix(G),
    #                          y=Y-as.vector(as.matrix(C)%*%betak[kkk-1,(K0+1):K]),
    #                          intercept=FALSE,
    #                          family=gaussian(),
    #                          penalty=adaptive.lasso(lambda=2*sqrt(2)/sqrt(g)*sigmak[kkk-1],
    #                                                 al.weights=as.vector(Elambdaj[1:K0])))
    # betak[kkk,1:K0]=Bayelsso_gen$coefficients
    Bayelsso_clinical = glmnet(x=as.matrix(C),
                               y=Y-as.vector(as.matrix(G)%*%betak[kkk,1:K0]),
                               intercept=FALSE, alpha=0,
                               family="gaussian", lambda=1/N)
    betak[kkk,(K0+1):K]=coef(Bayelsso_clinical)[-1]
    
    plot(betak[kkk,],main=kkk)
    #update sigma (correct the typo on denominator in the paper)
    sigmak[kkk]=as.numeric(((sqrt(2)*sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0]))+sqrt(2*(sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0])))^2+4*((t(Y-cbind(G,C)%*%betak[kkk,])%*%(Y-cbind(G,C)%*%betak[kkk,]))^1+betak[kkk,(K0:K)]%*%betak[kkk,(K0:K)]+2*d)*(N+K+2*c+2+L)))/(2*(N+K+2*c+2+L))))
    #sigmak[kkk]=as.numeric(((sqrt(2)*sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0]))+sqrt(2*(sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0])))^2+4*((t(Y_full[train_indx]-as.matrix(train[,-1])%*%betak[kkk,])%*%(Y_full[train_indx]-as.matrix(train[,-1])%*%betak[kkk,]))^1+betak[kkk,(K0:K)]%*%betak[kkk,(K0:K)]+2*d)*(N+K+2*c+2+L)))/(2*(N+K+2*c+2+L))))
    #update b
    Q2=function(bk){
      sum(-a*log(1/(t(bk)%*%t(Zmatrix)))-Elambdaj2[1:K0]*(t(bk)%*%t(Zmatrix)))+sum((alpha-1)*log(bk)-gam*bk)
    }
    # bbk[kkk,]=as.vector(optim(rep(1,nz),Q2,lower=rep(10^-10,nz),upper=rep(10,nz),method="L-BFGS-B",gr = NULL,control=list(fnscale=-1))$par)
    bbk[kkk,]=as.vector(optim(rep(1,nz),Q2,lower=rep(1e-7,nz),
                              upper=rep(1000,nz),method="L-BFGS-B",
                              gr = NULL,control=list(fnscale=-1))$par)
    #convergence criterion  
    betat=betak[kkk,]-betak[kkk-1,]
    bt=bbk[kkk,]-bbk[kkk-1,]
    if(sum(abs(betat))+sum(abs(bt))<10^(-3)){break}
    print(sum(abs(betat))+sum(abs(bt)))
  }
  lst=list(b=bbk,beta=betak,k=kkk,a=a,g=g)
  return(lst)
}

NEG_mcmc=function(Y,G,C,a0,gstr,Zmatrix,mypath){
  a=a0
  #Define LOG Parabolic Cylinder Function (using "fAsianOptions" package with Whittaker W function representation)
  write.table(" ",paste0(mypath,"out1",a,"_",gstr,".txt"),append=FALSE,row.names=FALSE,col.names=FALSE)
  lpcf=function(k,v,z){
    vz=as.matrix(t(c(k,v,z)))
    write.table(vz,paste0(mypath,"zv",a,"_",gstr,".txt"),sep=',',append=FALSE,row.names=FALSE,col.names=FALSE)
    system(paste0("python3 ", mypath, "lpcf.py", " ", mypath, " ", a," ",gstr))
    A=read.table(paste0(mypath,"out1",a,"_",gstr,".txt"))
    return(A)
  }
  alpha=1
  gam=1
  s1=1
  s2=2
  c=1
  d=1
  N = nrow(G)
  K0 = ncol(G)
  L = ncol(C)
  K = K0+L
  I=300
  #Store the values of expectations in E-step
  Elambdaj =as.vector(rep(0,K0))
  Elambdaj2=as.vector(rep(0,K0)) 
  nz=length(Zmatrix[1,])
  #initial values
  sigmak=rep(1,I)
  betak=matrix(0,nrow=I,ncol=K)
  bbk=matrix(1,nrow=I,ncol=nz)
  v1=as.vector(rep((-(2*a+1+s1)),K0))
  v2=as.vector(rep((-(2*a+1+s2)),K0))
  vb=as.vector(rep(-2*(a+0.5),K0))
  #g=N^(-2)
  if (gstr=='scale'){
    g0=1/(N^2)
  }else{
    g0 = as.numeric(gstr)
  }
  g=g0
  #starting iteration
  for(kkk in 2:I){
    #E-step
    print(bbk[kkk-1,])
    hf=1/((Zmatrix)%*%bbk[kkk-1,])
    ZZ=as.vector(abs(betak[kkk-1,1:K0])*sqrt(hf)/(sqrt(g)*sigmak[kkk-1]))
    lPCFb=lpcf(K0,vb,ZZ)
    lPCF1=lpcf(K0,v1,ZZ)
    lPCF2=lpcf(K0,v2,ZZ)
    logmar=as.vector(rep(0,K0))
    for (h in 1:K0){
      logmar[h]=as.numeric(log(a*(2^a)*sqrt(hf[h])/(sqrt(pi)*sqrt(g)*sigmak[kkk-1]))+lgamma(a+0.5)+(betak[kkk-1,h])^2*hf[h]/(4*(g)*sigmak[kkk-1]^2)+lPCFb[h])
    }
    for(h in 1:K0){
      Elambdaj[h]=as.numeric(exp(lgamma(2*a+s1+1)+(0.5*s1+0.5)*log(hf[h])-log(sqrt(g)*sigmak[kkk-1]*gamma(a)*2^(a+0.5*s1))+betak[kkk-1,h]^2*(hf[h])/(4*(g)*sigmak[kkk-1]^2)-logmar[h]+lPCF1[h]))
      Elambdaj2[h]=as.numeric(exp(lgamma(2*a+s2+1)+(0.5*s2+0.5)*log(hf[h])-log(sqrt(g)*sigmak[kkk-1]*gamma(a)*2^(a+0.5*s2))+betak[kkk-1,h]^2*(hf[h])/(4*(g)*sigmak[kkk-1]^2)-logmar[h]+lPCF2[h]))
    }
    
    #M-step
    #update Beta through Adaptive LASSO, Zou (2006) (using "lqa" package)
    Bayelsso_gen = glmnet(x=as.matrix(G),
                          y=Y-as.vector(as.matrix(C)%*%betak[kkk-1,(K0+1):K]),
                          intercept=FALSE, alpha=1,
                          family="gaussian", lambda=2*sqrt(2)/sqrt(g)*sigmak[kkk-1]*(sum(Elambdaj[1:K0])/N)/(2*N),
                          penalty.factor =as.vector(Elambdaj[1:K0]))
    betak[kkk,1:K0]=coef(Bayelsso_gen)[-1]
    # Bayelsso_gen=lqa.default(x=as.matrix(G),
    #                          y=Y-as.vector(as.matrix(C)%*%betak[kkk-1,(K0+1):K]),
    #                          intercept=FALSE,
    #                          family=gaussian(),
    #                          penalty=adaptive.lasso(lambda=2*sqrt(2)/sqrt(g)*sigmak[kkk-1],
    #                                                 al.weights=as.vector(Elambdaj[1:K0])))
    # betak[kkk,1:K0]=Bayelsso_gen$coefficients
    Bstar = solve(t(C)%*%C+diag(1,ncol(C))/(sigmak[kkk-1])^2)
    bstar = Bstar%*%t(C)%*%(Y-as.vector(as.matrix(G)%*%betak[kkk,1:K0]))/(sigmak[kkk-1]^2)
    betak[kkk,(K0+1):K]=mvrnorm(n=1,mu=bstar,Sigma=Bstar)
    
    plot(betak[kkk,],main=kkk)
    #update sigma (correct the typo on denominator in the paper)
    sigmak[kkk]=as.numeric(((sqrt(2)*sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0]))+sqrt(2*(sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0])))^2+4*((t(Y-cbind(G,C)%*%betak[kkk,])%*%(Y-cbind(G,C)%*%betak[kkk,]))^1+betak[kkk,(K0:K)]%*%betak[kkk,(K0:K)]+2*d)*(N+K+2*c+2+L)))/(2*(N+K+2*c+2+L))))
    #sigmak[kkk]=as.numeric(((sqrt(2)*sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0]))+sqrt(2*(sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0])))^2+4*((t(Y_full[train_indx]-as.matrix(train[,-1])%*%betak[kkk,])%*%(Y_full[train_indx]-as.matrix(train[,-1])%*%betak[kkk,]))^1+betak[kkk,(K0:K)]%*%betak[kkk,(K0:K)]+2*d)*(N+K+2*c+2+L)))/(2*(N+K+2*c+2+L))))
    #update b
    Q2=function(bk){
      sum(-a*log(1/(t(bk)%*%t(Zmatrix)))-Elambdaj2[1:K0]*(t(bk)%*%t(Zmatrix)))+sum((alpha-1)*log(bk)-gam*bk)
    }
    # bbk[kkk,]=as.vector(optim(rep(1,nz),Q2,lower=rep(10^-10,nz),upper=rep(10,nz),method="L-BFGS-B",gr = NULL,control=list(fnscale=-1))$par)
    bbk[kkk,]=as.vector(optim(rep(1,nz),Q2,lower=rep(1e-7,nz),
                              upper=rep(1000,nz),method="L-BFGS-B",
                              gr = NULL,control=list(fnscale=-1))$par)
    #convergence criterion  
    betat=betak[kkk,]-betak[kkk-1,]
    bt=bbk[kkk,]-bbk[kkk-1,]
    if(sum(abs(betat))+sum(abs(bt))<10^(-3)){break}
    print(sum(abs(betat))+sum(abs(bt)))
  }
  lst=list(b=bbk,beta=betak,k=kkk,a=a,g=g)
  return(lst)
}

NEG_censor=function(Y,G,C,a0,gstr,Zmatrix,mypath,Delta){
  a=a0
  #Define LOG Parabolic Cylinder Function (using "fAsianOptions" package with Whittaker W function representation)
  write.table(" ",paste0(mypath,"out1",a,"_",gstr,".txt"),append=FALSE,row.names=FALSE,col.names=FALSE)
  lpcf=function(k,v,z){
    vz=as.matrix(t(c(k,v,z)))
    write.table(vz,paste0(mypath,"zv",a,"_",gstr,".txt"),sep=',',append=FALSE,row.names=FALSE,col.names=FALSE)
    system(paste0("python3 ", mypath, "lpcf.py", " ", mypath, " ", a," ",gstr))
    A=read.table(paste0(mypath,"out1",a,"_",gstr,".txt"))
    return(A)
  }
  alpha=1
  gam=1
  s1=1
  s2=2
  c=1
  d=1
  N = nrow(G)
  K0 = ncol(G)
  L = ncol(C)
  K = K0+L
  I=300
  #Store the values of expectations in E-step
  Elambdaj =as.vector(rep(0,K0))
  Elambdaj2=as.vector(rep(0,K0)) 
  nz=length(Zmatrix[1,])
  #initial values
  sigmak=rep(1,I)
  betak=matrix(0,nrow=I,ncol=K)
  bbk=matrix(1,nrow=I,ncol=nz)
  v1=as.vector(rep((-(2*a+1+s1)),K0))
  v2=as.vector(rep((-(2*a+1+s2)),K0))
  vb=as.vector(rep(-2*(a+0.5),K0))
  Z = matrix(0,nrow=N,ncol=I)
  #g=N^(-2)
  if (gstr=='scale'){
    g0=1/(N^2)
  }else{
    g0 = as.numeric(gstr)
  }
  g=g0
  #starting iteration
  for(kkk in 2:I){
    #E-step
    print(bbk[kkk-1,])
    hf=1/((Zmatrix)%*%bbk[kkk-1,])
    ZZ=as.vector(abs(betak[kkk-1,1:K0])*sqrt(hf)/(sqrt(g)*sigmak[kkk-1]))
    lPCFb=lpcf(K0,vb,ZZ)
    lPCF1=lpcf(K0,v1,ZZ)
    lPCF2=lpcf(K0,v2,ZZ)
    logmar=as.vector(rep(0,K0))
    for (h in 1:K0){
      logmar[h]=as.numeric(log(a*(2^a)*sqrt(hf[h])/(sqrt(pi)*sqrt(g)*sigmak[kkk-1]))+lgamma(a+0.5)+(betak[kkk-1,h])^2*hf[h]/(4*(g)*sigmak[kkk-1]^2)+lPCFb[h])
    }
    for(h in 1:K0){
      Elambdaj[h]=as.numeric(exp(lgamma(2*a+s1+1)+(0.5*s1+0.5)*log(hf[h])-log(sqrt(g)*sigmak[kkk-1]*gamma(a)*2^(a+0.5*s1))+betak[kkk-1,h]^2*(hf[h])/(4*(g)*sigmak[kkk-1]^2)-logmar[h]+lPCF1[h]))
      Elambdaj2[h]=as.numeric(exp(lgamma(2*a+s2+1)+(0.5*s2+0.5)*log(hf[h])-log(sqrt(g)*sigmak[kkk-1]*gamma(a)*2^(a+0.5*s2))+betak[kkk-1,h]^2*(hf[h])/(4*(g)*sigmak[kkk-1]^2)-logmar[h]+lPCF2[h]))
    }
    
    #M-step
    #update Beta through Adaptive LASSO, Zou (2006) (using "lqa" package)
    Bayelsso_gen = glmnet(x=as.matrix(G),
                          y=Y-as.vector(as.matrix(C)%*%betak[kkk-1,(K0+1):K]),
                          intercept=FALSE, alpha=1,
                          family="gaussian", lambda=2*sqrt(2)/sqrt(g)*sigmak[kkk-1]*(sum(Elambdaj[1:K0])/N)/(2*N),
                          penalty.factor =as.vector(Elambdaj[1:K0]))
    betak[kkk,1:K0]=coef(Bayelsso_gen)[-1]
    # Bayelsso_gen=lqa.default(x=as.matrix(G),
    #                          y=Y-as.vector(as.matrix(C)%*%betak[kkk-1,(K0+1):K]),
    #                          intercept=FALSE,
    #                          family=gaussian(),
    #                          penalty=adaptive.lasso(lambda=2*sqrt(2)/sqrt(g)*sigmak[kkk-1],
    #                                                 al.weights=as.vector(Elambdaj[1:K0])))
    # betak[kkk,1:K0]=Bayelsso_gen$coefficients
    Bayelsso_clinical = glmnet(x=as.matrix(C),
                               y=Y-as.vector(as.matrix(G)%*%betak[kkk,1:K0]),
                               intercept=FALSE, alpha=0,
                               family="gaussian", lambda=1/N)
    betak[kkk,(K0+1):K]=coef(Bayelsso_clinical)[-1]
    
    plot(betak[kkk,],main=kkk)
    #update sigma (correct the typo on denominator in the paper)
    sigmak[kkk]=as.numeric(((sqrt(2)*sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0]))+sqrt(2*(sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0])))^2+4*((t(Y-cbind(G,C)%*%betak[kkk,])%*%(Y-cbind(G,C)%*%betak[kkk,]))^1+betak[kkk,(K0:K)]%*%betak[kkk,(K0:K)]+2*d)*(N+K+2*c+2+L)))/(2*(N+K+2*c+2+L))))
    #sigmak[kkk]=as.numeric(((sqrt(2)*sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0]))+sqrt(2*(sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0])))^2+4*((t(Y_full[train_indx]-as.matrix(train[,-1])%*%betak[kkk,])%*%(Y_full[train_indx]-as.matrix(train[,-1])%*%betak[kkk,]))^1+betak[kkk,(K0:K)]%*%betak[kkk,(K0:K)]+2*d)*(N+K+2*c+2+L)))/(2*(N+K+2*c+2+L))))
    #update b
    Q2=function(bk){
      sum(-a*log(1/(t(bk)%*%t(Zmatrix)))-Elambdaj2[1:K0]*(t(bk)%*%t(Zmatrix)))+sum((alpha-1)*log(bk)-gam*bk)
    }
    # bbk[kkk,]=as.vector(optim(rep(1,nz),Q2,lower=rep(10^-10,nz),upper=rep(10,nz),method="L-BFGS-B",gr = NULL,control=list(fnscale=-1))$par)
    bbk[kkk,]=as.vector(optim(rep(1,nz),Q2,lower=rep(1e-7,nz),
                              upper=rep(1000,nz),method="L-BFGS-B",
                              gr = NULL,control=list(fnscale=-1))$par)
    #convergence criterion  
    betat=betak[kkk,]-betak[kkk-1,]
    bt=bbk[kkk,]-bbk[kkk-1,]
    Z[,kkk]=rtruncnorm(n=length(Y),a=Y,b=Inf,
                 mean=as.vector(as.matrix(G)%*%betak[kkk,1:K0])+as.vector(as.matrix(C)%*%betak[kkk-1,(K0+1):K]),
                 sd=sigmak[kkk])
    Y[Delta==0]=Z[Delta==0,kkk]
    if(sum(abs(betat))+sum(abs(bt))<10^(-3)){break}
    print(sum(abs(betat))+sum(abs(bt)))
  }
  lst=list(b=bbk,beta=betak,k=kkk,a=a,g=g)
  return(lst)
}

rsq <- function(x, y) summary(lm(y~x))$r.squared
cindx = function(pred,actual){
  phi = 0
  phi_pred = 0
  n_input = length(actual)
  for (ci in 1:n_input){
    for (cj in 1:n_input){
      if (actual[ci]>actual[cj]){
        phi=phi+1
        if (pred[ci]>pred[cj]){
          phi_pred=phi_pred+1
        }
      }
    }
  }
  return(phi_pred/phi)
}