rm(list=ls())
library(gridExtra)
library(ggplot2)
library(glmnet)
library(MASS)
library(cBioPortalData)
library(dplyr)
library(survival)
library(coda)
library(caret)
set.seed(2021)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  source('functions.R')
  mechanistic_mod = 'quadratic'
  clinic_mod = 'em'
  # parameters for clinic model
  local = T
  dirr = getwd()
}else{
  source('functions_final.R')
  mechanistic_mod = args[1]
  clinic_mod = args[2]
  local = F
  dirr = args[5]
}

if (local==T){
  parentPath = getwd()
  if (mechanistic_mod=='cubic'){
    data_path = paste0(parentPath,"/GBM_data3.RData")
  }else if (mechanistic_mod=='quadratic'){
    data_path = paste0(parentPath,"/GBM_data2.RData")
  }else{
    data_path = paste0(parentPath,"/GBM_data1.RData")
  }
}else{
  parentPath = "/home/hx222/emvs/"
  if (mechanistic_mod=='cubic'){
    data_path = paste0(parentPath,"/GBM_data3.RData")
  }else if (mechanistic_mod=='quadratic'){
    data_path = paste0(parentPath,"/GBM_data2.RData")
  }else{
    data_path = paste0(parentPath,"/GBM_data1.RData")
  }
}
# if load subdata, uncomment this trunk
# load(data_path)
# Y = log(Y)
# Delta = myclinical$OS_STATUS

# if load full data, uncomment this trunk
# The loaded data includes gene expression, clinical data, overall survival time, sensorship and gene groupings
load(data_path)
G = GBM_data$G
C = GBM_data$C
Y = log(GBM_data$Y)
M = GBM_data$M
Delta = GBM_data$Delta

#
#
#
#
# ################################################################################
# #first-stage
# Making Grouping matrix Zmatrix
# run emvs algorithm
EMVS_res = EMVS(M,G,grouping=gene_group,nu0=0.5,nu1=10^3,nu=1,lambda=1,a=1,b=1)
save(EMVS_res,file=paste0('GBM_data2_EMVS_res.Rdata'))
load("GBM_data2_EMVS_res.Rdata")
N=nrow(G)
L=ncol(C)
J=ncol(M)
K=ncol(G)
if (length(args)==0){
  a0 = .1
  g0 = gstr = 1/(N^2)
}else{
  a0 = as.numeric(args[3])
  gstr = args[4]
}
mypath = paste0(getwd(),'/')
if (!file.exists(mypath)){
  dir.create(file.path(mypath))
}
# plot(EMVS_res$R2)
################################################################################
#second-stage

Zmatrix=cbind(matrix(1,nrow=K,ncol=1),matrix(0,nrow=K,ncol=3))
for(ww in 1:K){
  if(EMVS_res$R2[ww]>=0.8){Zmatrix[ww,2]=1}
  else{ if(EMVS_res$R2[ww]>=0.2){Zmatrix[ww,3]=1} else{Zmatrix[ww,4]=1}}
}

###  NEG   ###
###  NEG Group ###
# colnames(G)[abs(NEG_res$beta[NEG_res$k,])>1e-5]
# NEG_res$beta[NEG_res$k,abs(NEG_res$beta[NEG_res$k,])>1e-5]

#E-M loop
#hyperparameter

n_fold =10
set.seed(123)
folds <- createFolds(1:N, k = n_fold)
final_list = list()
print("a0")
print(a0)
print('g')
print(gstr)
for (jjj in 1:n_fold){
  test_indx = folds[[jjj]]
  train_delta = Delta[-test_indx]
  test_delta = Delta[test_indx]
  if (clinic_mod=='mcmc'){
    lst = NEG_mcmc(Y=Y[-test_indx,],G=as.matrix(G[-test_indx,]),
                   C=as.matrix(C[-test_indx,]),
                   a0=a0,gstr=gstr,Zmatrix=Zmatrix,mypath=mypath)
  }else if (clinic_mod=='em'){
    lst = NEG_em(Y=Y[-test_indx],G=as.matrix(G[-test_indx,]),
                 C=as.matrix(C[-test_indx,]),
                 a0=a0,gstr=gstr,Zmatrix=Zmatrix,mypath=mypath)
  }else if (clinic_mod=='censor'){
    lst = NEG_censor(Y=train[,1],G=as.matrix(train[,2:(K0+1)]),
                     C=as.matrix(train[,(K0+2):(K+1)]),a0=a0,gstr=gstr,
                     Zmatrix=Zmatrix,mypath=mypath,Delta=train_delta)
  }
  #final_list[[indx_g]][[indx_a]][[jjj]]=lst
  final_list[[jjj]]=lst
  save(final_list,file=paste0(mypath,"10fold_",a0,"_",gstr,"_",clinic_mod,"_censor.RData"))
}
foreach(a0=c(0.1,1,10,50)) %:% 
  foreach(gstr=c('scale')) %:%
    foreach(jjj=1:n_fold) %do%{
      test_indx = folds[[jjj]]
      train_delta = Delta[-test_indx]
      test_delta = Delta[test_indx]
      if (clinic_mod=='mcmc'){
        lst = NEG_mcmc(Y=Y[-test_indx,],G=as.matrix(G[-test_indx,]),
                       C=as.matrix(C[-test_indx,]),
                       a0=a0,gstr=gstr,Zmatrix=Zmatrix,mypath=mypath)
      }else if (clinic_mod=='em'){
        lst = NEG_em(Y=Y[-test_indx],G=as.matrix(G[-test_indx,]),
                     C=as.matrix(C[-test_indx,]),
                     a0=a0,gstr=gstr,Zmatrix=Zmatrix,mypath=mypath)
      }else if (clinic_mod=='censor'){
        lst = NEG_censor(Y=train[,1],G=as.matrix(train[,2:(K0+1)]),
                         C=as.matrix(train[,(K0+2):(K+1)]),a0=a0,gstr=gstr,
                         Zmatrix=Zmatrix,mypath=mypath,Delta=train_delta)
      }
      #final_list[[indx_g]][[indx_a]][[jjj]]=lst
      final_list[[jjj]]=lst
      save(final_list,file=paste0(mypath,"10fold_",a0,"_",gstr,"_",clinic_mod,".RData"))
    }

# #load(paste0(mypath,"10fold_",a0,"_",gstr,"_em.RData"))
# # load("~/Dropbox/harvard/emvs/code/results/mcmc_dist/10fold_0.5_1.RData")
# 
# rsq <- function(x, y) summary(lm(y~x))$r.squared
# res_table = data.frame(matrix(0,n_fold,19))
# colnames(res_table) = c('size support','r2_train','r2_test',
#                         'cindex_train','cindex_test','mse_train','mse_test',
#                         'beta_AGE_lower','beta_AGE_mean','beta_AGE_upper',
#                         'beta_PRIOR_GLIOMA_lower','beta_PRIOR_GLIOMA_mean','beta_PRIOR_GLIOMA_upper',
#                         'beta_SEX_lower','beta_SEX_mean','beta_SEX_upper',
#                         'beta_PRETREATMENT_HISTORY_lower','beta_PRETREATMENT_HISTORY_mean','beta_PRETREATMENT_HISTORY_upper')
# for (jjj in 1:n_fold){
#   test_indx = folds[[jjj]]
#   output = final_list[[jjj]]
#   output$beta = data.frame(output$beta)
#   colnames(output$beta) = colnames(G)
#   estbeta = output$beta[output$k,]
#   estbeta[(K0+1):K] =colMeans(output$beta[round(1/5*output$k):output$k,(K0+1):K])
#   pred_train = as.matrix(GenDat)[-test_indx,-1]%*%as.numeric(estbeta)
#   pred_test = as.matrix(GenDat)[test_indx,-1]%*%as.numeric(estbeta)
#   res_table[jjj,'size support']=sum(abs(estbeta)>1e-5)
#   res_table[jjj,'r2_train']=rsq(pred_train,Y[-test_indx])
#   res_table[jjj,'r2_test']=rsq(pred_test,Y[test_indx])
#   res_uotable[jjj,'cindex_train']=cindx(pred_train,Y[-test_indx])
#   res_table[jjj,'cindex_test']=cindx(pred_test,Y[test_indx])
#   res_table[jjj,'mse_train']=mean((pred_train-Y[-test_indx])^2)
#   res_table[jjj,'mse_test']=mean((pred_test-Y[test_indx])^2)
#   beta.mcmc = as.mcmc(output$beta[round(1/5*output$k):output$k,(K0+1):K])
# 
#   # find HPD interval using the CODA function
#   for (clifac in c("AGE","PRIOR_GLIOMA","SEX","PRETREATMENT_HISTORY")){
#     res_table[jjj,c(paste0('beta_',clifac,'_lower'),paste0('beta_',clifac,'_upper'))]=HPDinterval(beta.mcmc)[clifac,]
#     res_table[jjj,paste0('beta_',clifac,'_mean')]=estbeta[clifac]
#   }
# }
# #
#summarize cross-validation results
final_table = data.frame(matrix(0,nrow=8,ncol=13))
box_tables = data.frame()
colnames(final_table)[1:2] =c('a','g')
colnames(final_table)[3:13] = c('size support',
                                'r2_train',
                                'r2_test',
                                'cindex_train',
                                'cindex_test',
                                'mse_train',
                                'mse_test',
                                'AGE',
                                'PRIOR_GLIOMA',
                                'SEX',
                                'PRETREATMENT_HISTORY')

ppp=1
for (gstr in c('scale',1)){
  for (a0 in c(0.1,1,10,50)){
    fpath = paste0(mypath,"10fold_",
                   a0,"_",gstr,"_",clinic_mod,".RData")
    load(fpath)
    print(length(final_list))
    res_table = data.frame(matrix(0,n_fold,11))
    colnames(res_table) = c('size support','r2_train','r2_test',
                            'cindex_train','cindex_test','mse_train','mse_test',
                            'AGE','PRIOR_GLIOMA','SEX','PRETREATMENT_HISTORY')
    selected_biomarkers = c()
    for (jjj in 1:n_fold){
      test_indx = folds[[jjj]]
      output = final_list[[jjj]]
      output$beta = data.frame(output$beta)
      colnames(output$beta) = colnames(G)
      estbeta = output$beta[output$k,]
      selected_biomarkers = unique(c(selected_biomarkers,names(estbeta)[abs(estbeta)>1e-5]))
      pred_train = cbind(G[-test_indx,],C[-test_indx,])%*%as.numeric(estbeta)
      pred_test = cbind(G[test_indx,],C[test_indx,])%*%as.numeric(estbeta)
      res_table[jjj,'size support']=sum(abs(estbeta)>1e-5)
      res_table[jjj,'r2_train']=rsq(pred_train,Y[-test_indx])
      res_table[jjj,'r2_test']=rsq(pred_test,Y[test_indx])
      res_table[jjj,'cindex_train']=cindx(pred_train,Y[-test_indx])
      res_table[jjj,'cindex_test']=cindx(pred_test,Y[test_indx])
      res_table[jjj,'mse_train']=mean((pred_train-Y[-test_indx])^2)
      res_table[jjj,'mse_test']=mean((pred_test-Y[test_indx])^2)
      res_table[jjj,c("AGE","PRIOR_GLIOMA","SEX","PRETREATMENT_HISTORY")]=estbeta[1,1001:1004]
    }
    final_table[ppp,1]=a0
    final_table[ppp,2]=ifelse(gstr==1/(N^2),'scale',gstr)
    final_table[ppp,3:13]=round(colMeans(res_table),3)
    box_table=data.frame(a=a0,g=ifelse(gstr==1/(N^2),'scale',gstr),
                         AGE=res_table$AGE,
                         PRIOR_GLIOMA=res_table$PRIOR_GLIOMA,
                         SEX=res_table$SEX,
                         PRETREATMENT_HISTORY=res_table$PRETREATMENT_HISTORY)
    box_tables = rbind(box_tables,box_table)
    ppp=ppp+1
  }
}
box_tables$group = paste0('g=',box_tables$g,',a=',box_tables$a)
g1=ggplot(data=box_tables[1:50,],aes(x=as.factor(group),y=AGE))+
  geom_boxplot()+xlab(NULL)+ylab('beta age')+theme(axis.text.x = element_text(angle = 15))
g2=ggplot(data=box_tables[1:50,],aes(x=as.factor(group),y=PRIOR_GLIOMA))+
  geom_boxplot()+xlab(NULL)+ylab('beta prior glioma')+theme(axis.text.x = element_text(angle = 15))
g3=ggplot(data=box_tables[1:50,],aes(x=as.factor(group),y=SEX))+
  geom_boxplot()+xlab(NULL)+ylab('beta sex')+theme(axis.text.x = element_text(angle = 15))
g4=ggplot(data=box_tables[1:50,],aes(x=as.factor(group),y=PRETREATMENT_HISTORY))+
  geom_boxplot()+xlab(NULL)+ylab('beta pretreatment history')+theme(axis.text.x = element_text(angle = 15))
grid.arrange(g1,g2,g3,g4, ncol = 2)
#
final_table
paste0(final_table[1,1:7],collapse = '&')
