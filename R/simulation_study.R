Vg = .25
sim_n = 1e4
N = 100


simlxn<-function(Vg, sim_n, N){
  simdat<-rnorm(sim_n, 0, Vg)
  #bound the data below by -1
  for (i in 1:n){
    simdat[i]<-max(-1,simdat[i])
  }
  #transform data to be in terms of delta
  delta.dat<-(2^(simdat+1))-1

  #transform data to be in terms of r (assume M on eh1)
  r.dat<-delta.dat/(1+delta.dat)

  eh1<-numeric(sim_n)
  for (i in 1:sim_n){
    eh1[i]<-rbinom(1,1000,r.dat[i])
  }
  hist(eh1, probability = TRUE)
}

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

eh1<-numeric(length(r.dat))
for (i in (1:length(r.dat))){
  eh1[i]<-rbinom(1,1000,r.dat[i])
}

hist(eh1)


#transform data to be in terms of r (assume M on eh2)
r.dat.eh2<-1/(1+delta.dat)

eh2<-numeric(length(r.dat.eh2))
for (i in (1:length(r.dat.eh2))){
  eh2[i]<-rbinom(1,1000,r.dat.eh2[i])
}

hist(eh2)

###############################DON'T NEED########################################

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

