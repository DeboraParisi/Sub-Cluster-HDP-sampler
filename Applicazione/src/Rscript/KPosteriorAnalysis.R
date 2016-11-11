#This is an R Script, whit command to analyze and print
#the analysis of K in the algorithm

#Setwd
#setwd("~/Progetto_Bayesiana_PACS") 
#Apro il file per la scrittura dell'output
sink(file = "RChainOutput.txt", append = FALSE)
write("#####################################################", stdout())
write("####################            #####################",stdout())
write("#################### K ANALYSIS #####################", stdout())
write("####################            #####################",stdout())
write("#####################################################", stdout())
write("Length of chain before burnin/thinning: ", stdout())
print(N)
write("Burnin: ",stdout())
print(burnin)
write("Thinning: ",stdout())
print(thinning)
write("Length of chain after burnin/thinning",stdout())
print(after_burnin)

burnin = 0
thinning = 1

write("Loading AllK...",stdout())
png(filename = "AllK.png")
plot(AllK)
title(main="All K")
write("Done",stdout())
dev.off()

write("I'm creating mcmc object...",stdout())
k.mc <- mcmc(AllK,start=burnin+1,end=after_burnin,thin=thinning)

write("Summary of the chain",stdout())
print(summary(k.mc))

write("TracePlot",stdout())
png(filename = "Traceplot_AllK.png")
plot(k.mc,trace=TRUE,density=FALSE,main="AllK")
write("Done",stdout())
dev.off()

write("Density",stdout())
png(filename = "Density_AllK.png")
plot(k.mc,trace=FALSE,density=TRUE,main="AllK")
write("Done",stdout())
dev.off()

write("Autocorrelation",stdout())
png(filename = "Autocorrelation_AllK.png")
acf(k.mc,main = "AC all K")
write("Done",stdout())
dev.off()

write("Heidel diagnostic",stdout())
print(heidel.diag(k.mc))
write("Done",stdout())

write("Effective sample size",stdout())
print(effectiveSize(k.mc))
write("Done",stdout())

sink()
