
#Goal
#p.val<-integrate(integrand,-rad,rad,eh1,eh2, Eg_std, r0, p0, log_BinCoeff, abs.tol = 0, rel.tol = 1e-4)$value

rad=1.3
eh1=2150
eh2=2760
r0=.34889
p0=0.0003
Eg_std=.83
coeffx<-log_BinCoeffs(eh1+eh2)

p.val<-integrate(integrand,-rad,rad,eh1,eh2, Eg_std, r0, p0, log_BinCoeff, abs.tol = 0, rel.tol = 1e-4)$value
p.val<-integrate(integrand,-rad,rad,eh1,eh2, Eg_std, r0=r0, p0=p0, log_BinCoeff=coeffx, abs.tol = 0, rel.tol = 1e-4)$value


