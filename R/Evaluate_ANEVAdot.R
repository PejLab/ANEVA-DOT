#In this script, we produce simulated data and verify that the ANEVAdot test is calibrated properly

#Generate simulated data
simN  = 10000
Sg = runif(simN)*.25
N  = round(10^(3 + rnorm(simN,0,1)*.25)) #total count in binomial
N[N>1e4]<-1e4 #cap at 1e4
dE = rnorm(simN,0,1)*Sg
dE[dE<(-1)]<-(-1)
k  = 2^(dE+1)-1
r_h1 = k/(k+1)
#eh1  = rbinom(1,N, r_h1)
eh1 = numeric(length(N))
for (i in 1:length(N)){
  eh1[i]<-rbinom(1,N[i],r_h1[i])
}

test.data<-cbind(eh1,N-eh1)
test.data<-t(test.data)
write(test.data, file = "testdata.txt",ncolumns = 2, sep = "\t")
test.filepath<-"testdata.txt"

####add column titles manually####

start_time <- Sys.time()
test.result<-ANEVAdot(filepath = test.filepath, output_columns = c("eh1","eh2"),eh1 = "eh1", eh2 = "eh2",Eg_std=Sg)
end_time <- Sys.time()
runtime<-end_time - start_time
qqplot(runif(10000),test.result$p.val)

