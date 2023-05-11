library(tibble)
NN<-c(2000,5000,10000)
for(nn in 1:length(NN)){


h=0.01
n=NN[nn]
T=n*h
beta=1
alpha=3
sigma=0.2
x0=0
N=1000

print(n)
time=vector()
hatalpha=vector()
hatbeta=vector()
deltaL=vector()
for(a in 1:N){
  w=rnorm(n,0,sqrt(h))
  W=cumsum(w)
  
  X=vector()
  X[1]=x0
  for(k in 1:(n-1)){
    dl=max(0,-beta*h+alpha*X[k]*h-sigma*w[k])
    deltaL[k]=max(0,dl-X[k])
    X[k+1]=X[k]+beta*h-alpha*X[k]*h+sigma*w[k]+deltaL[k]
    if(X[k+1]<0){
      X[k+1]=0
    }else{
      X[k+1]=X[k+1]
    }
  }
  
  
  
  s_time <- Sys.time()
  deltaX=X[2:n]-X[1:n-1]
  X2=sum(X^2)
  X1=sum(X)
  X3=t(X[1:n-1])%*%(deltaX-deltaL)
  atop=(X[n]-X[1])*X1-n*X3
 
  btop=(X[n]-X[1])*X2-X3*X1
  bbottom=n*h*X2-X1^2*h
  
  hatalpha[a]=atop/bbottom
  hatbeta[a]=btop/bbottom
  e_time <- Sys.time()
  time[a]=e_time-s_time

  }


alphaBias=mean(hatalpha)-alpha
alphaStd.dev=sd(hatalpha)
alphaAsy.var=var(sqrt(T)*hatalpha)
alphaMSE=var(hatalpha)+alphaBias^2
betaBias=mean(hatbeta)-beta
betaStd.dev=sd(hatbeta)
betaAsy.var=var(sqrt(T)*hatbeta)
betaMSE=var(hatbeta)+betaBias^2


alpharesults=tibble(alphaBias,alphaStd.dev,alphaAsy.var,alphaMSE)
betaresults=tibble(betaBias,betaStd.dev,betaAsy.var,betaMSE)
print(alpharesults)
print(betaresults)

totaltime=sum(time)
print(totaltime)


EX1=hatbeta/hatalpha+sigma*dnorm(sqrt(2*hatalpha)*hatbeta/alpha/sigma)/sqrt(2*hatalpha)/(1-pnorm(-sqrt(2*hatalpha)*hatbeta/hatalpha/sigma))
EX2=sigma^2/2/hatalpha+hatbeta^2/hatalpha^2+hatbeta*sigma*dnorm(sqrt(2*hatalpha)*hatbeta/alpha/sigma)/sqrt(2*hatalpha)/(1-pnorm(-sqrt(2*hatalpha)*hatbeta/hatalpha/sigma))/hatalpha
sdalpha=sigma^2/(EX2-EX1^2)
sdbeta=sigma^2*(EX2^2-EX2*EX1^2)/(EX2-EX1^2)^2
  
Ialpha=(alpha<hatalpha+1.96*sqrt(sdalpha)/sqrt(T))*(alpha>hatalpha-1.96*sqrt(sdalpha)/sqrt(T))
CRalpha=sum(Ialpha)/N
print(CRalpha)

Ibeta=(beta<hatbeta+1.96*sqrt(sdbeta)/sqrt(T))*(beta>hatbeta-1.96*sqrt(sdbeta)/sqrt(T))
CRbeta=sum(Ibeta)/N
print(CRbeta)
}
qqnorm(hatalpha)
qqline(hatalpha)
qqnorm(hatbeta)
qqline(hatbeta)