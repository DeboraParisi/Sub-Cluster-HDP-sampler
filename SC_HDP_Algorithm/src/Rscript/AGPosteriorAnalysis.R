sink(file = "RChainOutput.txt", append = TRUE)

write("", stdout())


if(alpha_try == 'y'){
write("#####################################################", stdout())
write("##################                ###################",stdout())
write("################## ALPHA ANALYSIS ###################", stdout())
write("##################                ###################",stdout())
write("#####################################################", stdout())
write("Length of Alpha chain before burnin/thinning: ", stdout())
print(N_alpha)
write("Alpha Burnin: ",stdout())
print(alphaburnin)
write("Alpha Thinning: ",stdout())
print(alphathinning)
write("Length of Alpha chain after burnin/thinning",stdout())
print(after_alphaburnin)


alphaburnin <- 0
alphathinning <- 1

write("Alpha plot",stdout())
png(filename = "Alpha_plot.png")
plot(AllAlpha)
title(main="Alpha")
write("Done",stdout())
dev.off()

alpha_mean <- mean(AllAlpha)
alpha_var <- var(AllAlpha)


write("",stdout())
write("I'm creating Alpha MCMC object... ",stdout())
Alpha_mcmc <- mcmc(AllAlpha,start=alphaburnin + 1,end=after_alphaburnin,thin=alphathinning)
write("Alpha Summary",stdout())
print(summary(Alpha_mcmc))
write("",stdout())
write("I'm creating Alpha density...",stdout())
png(filename="Alpha_density.png")
densplot(Alpha_mcmc,show.obs=TRUE,main="Alpha Density")
lines(density(Alpha_mcmc),col="red",lwd=5)
#val <-seq(0,max(AllAlpha),0.1)

#density(AllAlpha,given.Rkern=false,type="l")
write("Done..",stdout())
dev.off()

write("",stdout())
write("I'm creating Alpha Traceplot... ",stdout())
png(filename="Alpha_traceplot.png")
plot(Alpha_mcmc,trace=TRUE,density=FALSE)
write("Done",stdout())
dev.off()


write("Autocorrelation",stdout())
png(filename = "Autocorrelation_Alpha.png")
acf(Alpha_mcmc,main="Autocorrelation Alpha")
write("Done",stdout())
dev.off()

write("Heidel diagnostic",stdout())
print(heidel.diag(Alpha_mcmc))
write("Done",stdout())

write("Effective sample size",stdout())
print(effectiveSize(Alpha_mcmc))
write("Done",stdout())

write("Done",stdout())

}

if(gamma_try == 'y'){
write("#####################################################", stdout())
write("##################                ###################",stdout())
write("################## GAMMA ANALYSIS ###################", stdout())
write("##################                ###################",stdout())
write("#####################################################", stdout())
write("Length of Gamma chain before burnin/thinning: ", stdout())
print(N_gamma)
write("Gamma Burnin: ",stdout())
print(gammaburnin)
write("Gamma Thinning: ",stdout())
print(gammathinning)
write("Length of Gamma chain after burnin/thinning",stdout())
print(after_gammaburnin)


gammaburnin <- 0
gammathinning <- 1

write("Gamma plot",stdout())
png(filename = "Gamma_plot.png")
plot(AllGamma)
title(main="Gamma")
write("Done",stdout())
dev.off()

write("",stdout())
write("I'm creating Gamma MCMC object... ",stdout())
Gamma_mcmc <- mcmc(AllGamma,start=gammaburnin + 1,end=after_gammaburnin,thin=gammathinning)
write("Alpha Summary",stdout())
print(summary(Gamma_mcmc))

write("",stdout())
write("I'm creating Gamma density... ",stdout())
png(filename="Gamma_density.png")
densplot(Gamma_mcmc,show.obs=TRUE,main="Gamma Density")
lines(density(Gamma_mcmc),col="red",lwd=5)
write("Done",stdout())
dev.off()


write("",stdout())
write("I'm creating Gamma traceplot... ",stdout())
png(filename="Gamma_traceplot.png")
plot(Gamma_mcmc,trace=TRUE,density=FALSE)
write("Done",stdout())
dev.off()


write("Autocorrelation",stdout())
png(filename = "Autocorrelation_Gamma.png")
acf(Gamma_mcmc,main="Autocorrelation Gamma")
write("Done",stdout())
dev.off()

write("Heidel diagnostic",stdout())
print(heidel.diag(Gamma_mcmc))
write("Done",stdout())

write("Effective sample size",stdout())
print(effectiveSize(Gamma_mcmc))
write("Done",stdout())

}


sink()
