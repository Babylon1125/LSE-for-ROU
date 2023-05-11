library(tibble)
NN<-c(2000,5000,10000)
for(nn in 1:length(NN)){
  
start_time <- Sys.time()
h=0.01
n=NN[nn]
T=n*h
theta=1
sigma=0.5
x0=0
N=1000

print(n)
time=vector()
hattheta=vector()
for(a in 1:N){
w=rnorm(n,0,sqrt(h))
W=cumsum(w)

X=vector()
deltaL=vector()
X[1]=x0
for(k in 1:(n-1)){
  dl=max(0,theta*X[k]*h-sigma*w[k])
  deltaL[k]=max(0,dl-X[k])
  X[k+1]=X[k]-theta*X[k]*h+sigma*w[k]+deltaL[k]
  if(X[k+1]<0){
    X[k+1]=0
  }else{
    X[k+1]=X[k+1]
  }
}

s_time<- Sys.time()
deltaX=X[2:n]-X[1:n-1]
top=t(X[1:n-1])%*%(deltaX-deltaL)


bottom=t(X)%*%X*h


hattheta[a]=-top/bottom

e_time <- Sys.time()
time[a]=e_time-s_time
}

mean(hattheta)

Bias=mean(hattheta)-theta
Std.dev=sd(hattheta)
Asy.var=var(sqrt(T)*hattheta)
MSE=var(hattheta)+Bias^2

results=tibble(Bias,Std.dev,Asy.var,MSE)
print(results)

totaltime=sum(time)
print(totaltime)


I=(theta<hattheta+1.96*sqrt(2*hattheta)/sqrt(T))*(theta>hattheta-1.96*sqrt(2*hattheta)/sqrt(T))
CR=sum(I)/N
print(CR)
}
qqnorm(hattheta)
qqline(hattheta)
