library(foreach)
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
  foreach(r=1:R) %do%{
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
  #estimate random intercept of each functional cluster
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

NEG_em=function(Y,G,C,a0,gstr,Zmatrix,mypath){
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
                          intercept=FALSE, standardize=F, alpha=1,
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
                               intercept=FALSE, standardize=F, alpha=0,
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