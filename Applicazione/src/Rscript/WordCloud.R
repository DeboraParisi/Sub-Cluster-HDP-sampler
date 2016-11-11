sink(file = "RWordcloud.txt")
#lista
#print(Weights)
#scalare
#print(Dim)
#vettore di numeri
#print(Labels)
#scalare
#print(Dim_Voc)
#vettore di parole
#print(Voc)
#lista
#print(Words)
#numero di parole da visualizzare
#print(Dim_w)

for(i in 1:Dim){

	
	idx_word <- Words[[i]]
	#print("Indice parole")
	#print(idx_word)
	#print("Frequenze")
	freq <- Weights[[i]]
	#print(freq)
	cloud <- cbind(Voc[idx_word[1]])
	for(y in 2:Dim_w){
		cloud <- cbind(cloud,Voc[idx_word[y]])			
	}
	cloud <- as.character(cloud)
	cat("Words inside topic ",Labels[i],": \n")
	print(cloud,quote=FALSE)
	write(" ",stdout())
	filename<-file.path( paste("WordCloud of Topic ",Labels[i],".pdf",sep=""))
	pdf(filename)
	wordcloud(cloud,freq,scale=c(4,0.5),min.freq=1,random.order=TRUE,colors=brewer.pal(Dim_w, "Dark2"))
	dev.off()
	
}
sink()



