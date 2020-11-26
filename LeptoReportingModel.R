library(MCMCpack)
# Load data file containing:
#   years: vector of study years
#   yrs: length of vector of study years (=18)
#   tps: total number of timepoints (=216)
#   occupations: vector of occupation names (including unspecified)
#   occs: vector of cocupation names (not including unspecified)
#   occ: length of occupations vector (=6)
#   serovars: vector of serovar names (including unknown)
#   srv: length of serovars (=5)
#   cases: an array of dimension occ x srv x tps with entry i,j,t containing the number of lepto cases from occupation i, of serovar j, at timepoint t.
load('SummaryData.Rdata')
getB<-function(eta) {
  return(cbind(diag(rep(eta,srv-1)),rep(1-eta,srv-1)))
}
# log density function of the Dirichlet distribution
lddirichlet<-function(x,a) {
  wh<-which(a!=0)
  return(lgamma(sum(a[wh]))+sum((a[wh]-1)*log(x[wh])-lgamma(a[wh])))
}
# function y(t): gives the (index of the) year of timepoint t.
getYear<-function(t) {return(1+(t-1)%/%12)}
# function m(t): gives the (index of the) month of timepoint t.
getMonth<-function(t) {return(1+(t-1)%%12)}
# Build demographic data matrix.
n<-matrix(NA,4,yrs)
tot.raw<-c(3820749,4027947,4242048)
n.raw<-matrix(c(26331,24795,26577,46983,46194,44976,20157,18858,14997,NA,NA,NA),3,4)
n.raw[,4]<-tot.raw-apply(n.raw[,1:3],1,sum)
y.raw<-c(2001,2006,2013)
for (i in 1:4) {
  m<-lm(log(n.raw[,i]) ~ y.raw)
  n[i,]<-round(exp(predict(m,list(y.raw=years))),0)
}
# summarise case data for quick reference
cases.by.month<-apply(array(cases,c(occ,srv,12,yrs)),c(1,3),sum)
cases.by.sv.occ<-apply(cases,1:2,sum)
cases.by.year<-array(NA,c(occ,srv,yrs))
for (k in 1:yrs) {
  cases.by.year[,,k]<-apply(cases[,,1:12+(k-1)*12],1:2,sum)
}
# Priors
alpha<-rep(0.8/(srv-1),srv-1) # for p
a<-1 # for lambda
b<-1
a.eta<-0.5 # for eta
b.eta<-0.5
delta<-matrix(0,occ,occ-2) # for phi
prior<-0.8/3
delta[1,1]<-prior
delta[2,2]<-prior
delta[3,3]<-prior
delta[4,4]<-prior
delta[5:6,]<-prior
mu.alpha<-rep(0.5,12)
# MCMC parameters
thinning<-10
iters<-2000
burnin<-500
# Starting values
lambda<-matrix(NA,occ-2,yrs)
mu<-matrix(NA,occ-2,12) 
p<-matrix(NA,occ-2,srv-1) 
phi<-matrix(NA,occ,occ-2) 
for (i in 1:(occ-2)) {# add proper starting values
  phi[,i]<-rdirichlet(1,delta[,i])
  mu[i,]<-rdirichlet(1,mu.alpha+cases.by.month[i,])
  for (k in 1:yrs) {
    lambda[i,k]<-rgamma(1,a+sum(cases[i,,1:12+(k-1)*12]),b+n[i,k]*phi[i,i])
  }
  p[i,]<-rdirichlet(1,alpha+cases.by.sv.occ[i,1:(srv-1)])
}
eta<-rbeta(1,a.eta+sum(cases[,-srv,]),b.eta+sum(cases[,srv,]))
# Monitors
lambda.accept<-matrix(0,occ-2,yrs)
lambda.reject<-matrix(0,occ-2,yrs)
p.accept<-rep(0,occ-2)
p.reject<-rep(0,occ-2)
phi.accept<-rep(0,occ-2)
phi.reject<-rep(0,occ-2)
mu.accept<-rep(0,occ-2)
mu.reject<-rep(0,occ-2)
# adaptive parameters
phi.beta<-c(39,26,66,37)#rep(10,occ-2)
mu.beta<-rep(10,occ-2)
adapt<-TRUE
# Storage
lambda.stored<-array(NA,c(occ-2,yrs,iters))
p.stored<-array(NA,c(occ-2,srv-1,iters))
phi.stored<-array(NA,c(occ,occ-2,iters))
mu.stored<-array(NA,c(occ-2,12,iters))
eta.stored<-rep(NA,iters)
# precomputation
B<-getB(eta)
# MCMC iterations start here
for (it in 1:(burnin+iters)) {
  for (th in 1:thinning) {
    for (i in 1:(occ-2)) {
      # update lambda
      for (k in 1:yrs) {
        proposal<-lambda
        proposal[i,k]<-rgamma(1,a+sum(cases[i,,1:12+(k-1)*12]),b+n[i,k]*phi[i,i])
        lar<-sum(dpois(cases[(occ-1):occ,1:srv,1:12+(k-1)*12],
             rep(phi[(occ-1):occ,]%*%diag((n[,k]*proposal[,k]))%*%p%*%B,12)*mu[i,],log=TRUE))-
             sum(dpois(cases[(occ-1):occ,1:srv,1:12+(k-1)*12],
             rep(phi[(occ-1):occ,]%*%diag((n[,k]*  lambda[,k]))%*%p%*%B,12)*mu[i,],log=TRUE))
        u<-runif(1)
        if (log(u)<=lar) {
          lambda[i,k]<-proposal[i,k]
          lambda.accept[i,k]<-lambda.accept[i,k]+1
        } else {
          lambda.reject[i,k]<-lambda.reject[i,k]+1
        }
      }
      # update p
      proposal<-p
      proposal[i,]<-rdirichlet(1,alpha+cases.by.sv.occ[i,1:(srv-1)])
      lar<-0
      for (k in 1:tps) {
        lar<-lar+sum(dpois(cases[(occ-1):occ,1:(srv-1),k],
                 phi[(occ-1):occ,]%*%diag(n[,getYear(k)]*lambda[,getYear(k)]*mu[,getMonth(k)])%*%proposal*eta,log=TRUE))-
                 sum(dpois(cases[(occ-1):occ,1:(srv-1),k],
                 phi[(occ-1):occ,]%*%diag(n[,getYear(k)]*lambda[,getYear(k)]*mu[,getMonth(k)])%*%p       *eta,log=TRUE))
      }
      u<-runif(1)
      if (log(u)<=lar) {
        p[i,]<-proposal[i,]
        p.accept[i]<-p.accept[i]+1
      } else {
        p.reject[i]<-p.reject[i]+1
      }
      # update phi
      proposal<-phi
      proposal[,i]<-rdirichlet(1,delta[,i]+phi.beta[i]*phi[,i])
      lar<-lddirichlet(phi[,i]     ,delta[,i]+phi.beta[i]*proposal[,i])-
           lddirichlet(proposal[,i],delta[,i]+phi.beta[i]*phi[,i]     )
      for (k in 1:tps) {
        lar<-lar+sum(dpois(cases[,,k],proposal%*%diag(n[,getYear(k)]*lambda[,getYear(k)]*mu[,getMonth(k)])%*%p%*%B,log=TRUE))-
                 sum(dpois(cases[,,k],phi     %*%diag(n[,getYear(k)]*lambda[,getYear(k)]*mu[,getMonth(k)])%*%p%*%B,log=TRUE))
      }       
      u<-runif(1)
      if (log(u)<=lar) {
        phi[,i]<-proposal[,i]
        phi.accept[i]<-phi.accept[i]+1
        if (adapt) {phi.beta[i]<-max(0,phi.beta[i]-3)}
      } else {
        phi.reject[i]<-phi.reject[i]+1
        if (adapt) {phi.beta[i]<-phi.beta[i]+1}
      }
      # update mu
      proposal<-mu
      proposal[i,]<-rdirichlet(1,mu.alpha+cases.by.month[i,]+mu.beta[i]*mu[i,])
      lar<-lddirichlet(proposal[i,],mu.alpha+cases.by.month[i,])+lddirichlet(      mu[i,],mu.alpha+cases.by.month[i,]+mu.beta[i]*proposal[i,])-
           lddirichlet(      mu[i,],mu.alpha+cases.by.month[i,])-lddirichlet(proposal[i,],mu.alpha+cases.by.month[i,]+mu.beta[i]*mu[i,])
      for (k in 1:tps) {
        lar<-lar+sum(dpois(cases[(occ-1):occ,,k],phi[(occ-1):occ,]%*%diag(n[,getYear(k)]*lambda[,getYear(k)]*proposal[,getMonth(k)])%*%p%*%B,log=TRUE))-
                 sum(dpois(cases[(occ-1):occ,,k],phi[(occ-1):occ,]%*%diag(n[,getYear(k)]*lambda[,getYear(k)]*      mu[,getMonth(k)])%*%p%*%B,log=TRUE))
      }
      u<-runif(1)
      if (log(u)<=lar) {
        mu[i,]<-proposal[i,]
        mu.accept[i]<-mu.accept[i]+1
        if (it<burnin) {mu.beta[i]<-max(0,mu.beta[i]-3)} 
      } else {
        mu.reject[i]<-mu.reject[i]+1
        if (it<burnin) {mu.beta[i]<-mu.beta[i]+1}
      }      
    }
    # update eta
    eta<-rbeta(1,a.eta+sum(cases[,-srv,]),b.eta+sum(cases[,srv,]))
    B<-getB(eta)
  }
  if (it>burnin) {
    lambda.stored[,,it-burnin]<-lambda
    p.stored[,,it-burnin]<-p
    phi.stored[,,it-burnin]<-phi
    mu.stored[,,it-burnin]<-mu
    eta.stored[it-burnin]<-eta
  } else if (it==burnin) {
    adapt<-FALSE
  }
}
save.image(file='LeptoOutput.Rdata')
par(mfcol=c(2,5))
for (i in 1:(occ-2)) {
  sv.cols<-c("red","blue","violet","darkorange")
  plot(1:iters,ylim=c(0,0.8),t="n",main=occs[i],ylab=round(p.accept[i]/(p.accept[i]+p.reject[i]),3),xlab="p")
  for (j in 1:(srv-1)) {
    lines(1:iters,p.stored[i,j,],col=sv.cols[j])
  }
  points(rep(-10,srv-1),apply(cases[i,1:(srv-1),],1,sum)/sum(cases[i,1:(srv-1),]),col=sv.cols)
  points(rep(iters+10,srv-1),apply(cases[i,1:(srv-1),],1,sum)/sum(cases[i,1:(srv-1),]),col=sv.cols)
  legend("topright",serovars[-srv],lty=1,col=sv.cols)
  cols<-rainbow(occ)
  plot(1:iters,ylim=c(0,1),t="n",main=occs[i],ylab=round(phi.accept[i]/(phi.accept[i]+phi.reject[i]),3),xlab=phi.beta[i])
  for (j in 1:occ) {
    lines(1:iters,phi.stored[j,i,],col=cols[j])
  }
  #legend("right",occupations,lty=1,col=cols)
}
#
lambda.quantiles<-apply(lambda.stored,1:2,quantile,probs=c(0.05,0.5,0.95))
plot(years,t="n",ylim=c(0,max(lambda.quantiles[3,,]*n)),xlim=c(years[1],years[yrs]),main="Annual cases",ylab="90% CI",xlab="study year")
for (i in 1:(occ-2)) {
  segments(years+0.05*(i-2),lambda.quantiles[1,i,]*n[i,],years+0.05*(i-2),lambda.quantiles[3,i,]*n[i,],col=cols[i])
  points(years+0.05*(i-2),lambda.quantiles[2,i,]*n[i,],col=cols[i])
}
#
plot(1:iters,eta.stored,t="l",ylim=c(0,1),main="Proportion with serovar",ylab="")
legend("bottom",occupations,lty=1,col=cols)
#
e.cases<-array(NA,c(occ,occ-2,iters))
for (it in 1:iters) {
  e.cases[,,it]<-phi.stored[,,it]%*%diag(n*lambda.stored[,,it])
}
mean.cases<-apply(e.cases,1:2,mean)
prop.cases<-mean.cases/apply(mean.cases,1,sum)
#
par(mfrow=c(2,6))
for (k in 1:12) {
  plot(1:iters,t="n",ylim=c(0,0.2),main=month.abb[k],xlab="iteration",ylab="")
  for (i in 1:(occ-2)) {
    lines(1:iters,mu.stored[i,k,],col=cols[i])
    points(c(-10,iters+10),rep(cases.by.month[i,k]/sum(cases.by.month[i,]),2),col=cols[i])
  }
}
#
# Figure plots
#
cols<-c("#FF0000","#009900","#0000FF","#FF00FF")
shade.cols<-c("#FF000066","#00990066","#0000FF66","#FF00FF66")
pdf("Years.pdf",width=6,height=6,pointsize=8)
par(mfrow=c(2,2))
for (i in 1:4) { # max boundary might not work now n is matrix
  plot(years,t="n",xlab="Year",ylab="Expected number of cases",ylim=c(0,max(rep(n,each=3)*lambda.quantiles)),main=occupations[i],xlim=c(min(years),max(years)))
  polygon(c(years,rev(years)),c(n[i,]*lambda.quantiles[1,i,],rev(n[i,]*lambda.quantiles[3,i,])),col=shade.cols[i],border=NA)
  lines(years,n[i,]*lambda.quantiles[2,i,],col=cols[i],lwd=2)
  points(years,apply(cases.by.year[i,,],2,sum),pch=4,col=cols[i])
}
dev.off()
#
mu.quantiles<-apply(mu.stored,1:2,quantile,probs=c(0.05,0.5,0.95))
pdf("Months.pdf",width=6,height=6,pointsize=8)
par(mfrow=c(2,2))
for (i in 1:4) {
  plot(1:12,t="n",xlab="Month",ylab="Proportion of cases",ylim=c(0,max(mu.quantiles)),main=occupations[i],xlim=c(1,12),xaxt="n")
  axis(1,1:12,substr(month.abb,1,1))
  lines(c(1,12),rep(1/12,2),col="grey")
  polygon(c(1:12,12:1),c(mu.quantiles[1,i,],rev(mu.quantiles[3,i,])),col=shade.cols[i],border=NA)
  lines(1:12,mu.quantiles[2,i,],col=cols[i],lwd=2)
  #cis<-MultinomCI(cases.by.month[i,],conf.level=0.9,method="goodman")
  #segments(1:12,cis[,2],1:12,cis[,3],col=cols[i])
  points(1:12,cases.by.month[i,]/sum(cases.by.month[i,]),pch=4,col=cols[i])
}
dev.off()
#
attributed.cases<-array(NA,c(2,occ-2,iters-burnin))
attributed.props<-array(NA,c(2,occ-2,iters-burnin))
for (it in 1:(iters-burnin)) {
  for (k in 1:2) {
    attributed.cases[k,,it]<-phi.stored[4+k,,it]*apply(n*lambda.stored[,,it],1,sum)
    attributed.props[k,,it]<-attributed.cases[k,,it]/sum(attributed.cases[k,,it])
  }
}
#
epsilon<-0.25
occupations.short<-occupations
occupations.short[c(1:2,5)]<-c("Dairy\nfarmer","Dry stock\nfarmer","Farmer not\nspecified")
pdf("Attribution.pdf",width=6,height=3,pointsize=8)
par(mfrow=c(1,2))
attribution<-apply(attributed.props,1:2,quantile,probs=c(0.05,0.5,0.95))
for (i in 1:2) {
  plot(1:4,t="n",xaxt="n",ylab="Proportion of cases",ylim=c(0,max(attribution)),main=occupations[4+i],xlim=c(1-epsilon,4+epsilon),xlab="")
  axis(1,1:4,occupations.short[1:4],mgp=c(3,1.5,0),cex.axis=0.93)
  segments(1:4,attribution[1,i,],1:4,attribution[3,i,],col=cols[1:4])
  points(1:4,attribution[2,i,],col=cols[1:4])
}
dev.off()
#
library(DescTools)
p.quantiles<-apply(p.stored,1:2,quantile,probs=c(0.05,0.5,0.95))
pdf("Serovar.pdf",width=6,height=6,pointsize=8)
par(mfrow=c(2,2))
for (i in 1:(occ-2)) {
  plot(1:4,t="n",xaxt="n",ylab="Proportion of cases",ylim=c(0,max(p.quantiles)),main=occupations[i],xlim=c(1-epsilon,4+epsilon),xlab="")
  axis(1,1:4,serovars[1:4])
  segments(1:4-epsilon/2,p.quantiles[1,i,],1:4-epsilon/2,p.quantiles[3,i,],col=cols[i])
  points(1:4-epsilon/2,p.quantiles[2,i,],col=cols[i])
  cis<-MultinomCI(cases.by.sv.occ[i,1:4],conf.level=0.9,method="goodman")
  segments(1:4+epsilon/2,cis[,2],1:4+epsilon/2,cis[,3],col=cols[i])
  points(1:4+epsilon/2,cases.by.sv.occ[i,1:4]/sum(cases.by.sv.occ[i,1:4]),pch=4,col=cols[i])
}
dev.off()
#
phi.quantiles<-apply(phi.stored,1:2,quantile,probs=c(0.05,0.5,0.95))
pdf("Reporting.pdf",width=6,height=6,pointsize=8)
par(mfrow=c(2,2),mar=c(5.5,4,4,2)+0.1)
for (i in 1:(occ-2)) {
  plot(1:occ,t="n",xaxt="n",ylab="Proportion of cases",ylim=c(0,max(phi.quantiles)),main=occupations[i],xlim=c(1-epsilon,6+epsilon),xlab="")
  axis(1,1:occ,occupations.short,las=2)#mgp=c(3,1.5,0),cex.axis=0.9)
  segments(1:occ,phi.quantiles[1,,i],1:occ,phi.quantiles[3,,i],col=cols[i])
  points(1:occ,phi.quantiles[2,,i],col=cols[i])
}
dev.off()
#
pdf("Theta.pdf",width=4,height=4,pointsize=8)
par(mar=c(5,4,4,2)+0.1)
hist(eta.stored,col="lightgrey",prob=TRUE,main="Probability of obtaining serovar",xlab="",xlim=c(0.6,0.8),breaks=40)
cis<-prop.test(sum(cases.by.sv.occ[,1:4]), sum(cases.by.sv.occ), conf.level=0.90, correct = FALSE)$conf.int
segments(cis[1],-1,cis[2],-1)
points(sum(cases.by.sv.occ[,1:4])/sum(cases.by.sv.occ),-1,pch=4)
dev.off()
#
# Output Table of lambdas
#
lambda.med<-apply(lambda.stored,1:2,median)
file<-"lambdaTable.txt"
cat("& Dairy farmer & Dry stock farmer & Meatworker & Other\\\\\\hline\n",file=file)
for (k in 1:yrs) {
  cat(years[k]," & ",paste(round(10^5*lambda.med[,k],2),collapse=" & "),"\\\\\n",file=file,append=TRUE)
}

  