data<-read.table(file="ecumor.str",header=FALSE) #a modified pyRAD structure output
#basically, I manually enter the population information in the 2nd column for each line

n.ind <- nrow(data);#count the number of samples
n.loci <- ncol(data) - 2;#count the number of loci

#identify loci with shared by less than 20% of the samples
low.cov.loci <- NULL;#create an empty vector
cov.loci <- NULL;#don't really need this; a vector to store the number of coverage per locus
for(iter in 3:ncol(data)){
	no.data <- sum(data[,iter] == -9);
	if(no.data > n.ind*0.25){
		low.cov.loci <- c(low.cov.loci, iter);
	}
	cov.loci <- c(cov.loci, no.data)
}

newdata <- data[,-low.cov.loci];#a new file with loci less than 75% coverage deleted
write.table(newdata,file='ecumor_c75.str',quote=FALSE,sep= '\t', row.names=FALSE, col.names=FALSE);

