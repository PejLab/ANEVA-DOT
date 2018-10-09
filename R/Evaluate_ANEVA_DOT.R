#In this script, we produce simulated data and verify that the ANEVAdot test is calibrated properly

#Generate simulated data
simN  = 1000
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
#test.data<-t(test.data)
colnames(test.data) <- c("eh1", "eh2")
test.data<-as.data.frame(test.data)

start_time <- Sys.time()
test.result<-ANEVA_DOT(test.data, output_columns = c("eh1","eh2"),
                       eh1 = "eh1", eh2 = "eh2",Eg_std=Sg)
end_time <- Sys.time()
runtime<-end_time - start_time #Time difference of 24.45896 secs
qqplot(runif(1000),test.result$p.val)
abline(a = 0, b = 1, col = "blue")



