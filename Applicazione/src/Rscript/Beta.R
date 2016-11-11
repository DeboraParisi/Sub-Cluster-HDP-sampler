#This is an R Script, with command to visualize Beta
pdf("Beta.pdf")
V<-Beta[[1]]
b_dim <- length(V)
tau<- seq(0,b_dim-1,1)
oth <- order(tau)
plot(c(min(tau)-1,tau[oth],max(tau)+1),c(0,cumsum(V[oth]),1),type="s",col="blue", ylab = "Topic weight", xlab = "Topic label",xaxt="n")

for(i in 2:Dim){
V<-Beta[[i]]
b_dim <- length(V)
tau<- seq(0,b_dim-1,1)
oth <- order(tau)
lines(c(min(tau)-1,tau[oth],max(tau)+1),c(0,cumsum(V[oth]),1),type="s",col="blue")
}

V<- BetaMean
b_dim <- length(V)
tau<- seq(0,b_dim-1,1)
oth <- order(tau)
lines(c(min(tau)-1,tau[oth],max(tau)+1),c(0,cumsum(V[oth]),1),type="s",col="black",lwd=5)

if(BestClustering > 0){
	V <- Beta[[BestClustering]]
	b_dim <- length(V)
	tau<- seq(0,b_dim-1,1)
	oth <- order(tau)
	lines(c(min(tau)-1,tau[oth],max(tau)+1),c(0,cumsum(V[oth]),1),type="s",col="red",lwd = 1.5)
	legend(x = "topleft", c("BestClustering"),lty=1,lwd=2.5,col="red")
}


title("Traiettorie DP globale")
axis(1,at=at ,labels=TopicLabels,las=0)

dev.off()



