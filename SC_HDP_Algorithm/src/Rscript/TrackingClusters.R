
Matrix<-cbind(Beta[[1]])

if(Dim>1){
	for(i in 2:Dim){
		V<-Beta[[i]]
		Matrix <- cbind(Matrix,V)
	}
}

x<-rep(0,Dim)
y<-Matrix[1,]
Lab<-as.character(Labels)
pdf("Cluster.pdf")
matplot(Matrix,type="l",lty=1,lwd=1.5,col=rainbow(Dim),main="Topic during last 100 iterations",ylab = "Topic weight")
legend("topleft", title="Topics",Lab, col=rainbow(Dim),lty=1)
dev.off()

