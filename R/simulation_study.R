
n = 1e4
simdat<-rnorm(n,0,.25)
#bound the data below by -1
for (i in 1:n){
  simdat[i]<-max(-1,simdat[i])
}

#transform data to be in terms of delta
delta.dat<-(2^(simdat+1))-1

#transform data to be in terms of r (assume M on eh1)
r.dat<-delta.dat/(1+delta.dat)

cutoffs<-quantile(r.dat,probs=c(.025,.975))
hist(r.dat)
abline(v=cutoffs[1], col = "red")
abline(v=cutoffs[2], col = "red")

#transform data to be in terms of r (assume M on eh2)
r.dat.eh2<-1/(1+delta.dat)

cutoffs2<-quantile(r.dat.eh2,probs=c(.025,.975))
hist(r.dat.eh2)
abline(v=cutoffs2[1], col = "red")
abline(v=cutoffs2[2], col = "red")

