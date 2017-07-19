#histograms for theta from number of segregating sites
#calculating number of variable sites
library(pegas);
input.data <- scan('Dy.AmaSP_all_loci.txt', what = 'character', sep = '\n');#change file name
id.lines <- grep('^\\d+[[:space:]]\\d+', input.data); 
working.lines <- c(id.lines, length(input.data) + 1);

n.var.sites <- NULL;
theta.vec <- NULL;

for(i in 1:length(id.lines)){
	
	start.line <- working.lines[i] + 1;
	end.line <- working.lines[i + 1] - 1;
	work.range <- c(start.line : end.line);
	
	s.pattern <- '\t.+[[:space:]]+';
	rep.pattern <- '';
	
	var.sites.vec <- NULL;
	
	if(start.line != end.line){
			
			for(ii in 1 : length(work.range) - 1){
		
			seq.1 <- sub(s.pattern, rep.pattern, input.data[work.range[ii]]);
			seq.2 <- sub(s.pattern, rep.pattern, input.data[work.range[ii + 1]]);		
			str.1 <- unlist(strsplit(seq.1, ''));
			str.2 <- unlist(strsplit(seq.2, ''));
			temp.seq <- which(str.1 != str.2 & str.1 != '-' & str.2 != '-' & str.1 != 'N' & str.2 != 'N');
			var.sites.vec <- c(var.sites.vec, temp.seq);
		
		}
				
		temp.var.site <- length(unique(var.sites.vec));
		n.var.sites <- c(n.var.sites, temp.var.site);
		n.ind.seq <- length(work.range);
		temp.theta <- (theta.s(temp.var.site, n.ind.seq))/120;
		theta.vec <- c(theta.vec, temp.theta)
	
	}	
}


hist(theta.vec, xlim = c(-0.01, 0.05), main = 'Estimated THETA\nThe Number of Segregating Sites Per Site', xlab = 'Theta.s/Sites');


#subset bpp; in the following case I chose 100 loci that have 2 to 10 var.sites; only an example.
#you can specify whatever you want
lower.limit <- 0; #at least how many var sites in a alignment#change this according to your need
upper.limit <- 6; #at most how many var sites are allowed in an alignment
selected.loci <- which(n.var.sites > lower.limit & n.var.sites < upper.limit);
chose.loci <- sample(selected.loci, 50); #randomly choose 100 loci from the selected loci that have desired number of var sites

for(i in 1:length(chose.loci)){
	temp.start.line <- working.lines[chose.loci[i]];
	temp.end.line <- working.lines[chose.loci[i] + 1] - 1;
	temp.wanted.lines <- temp.start.line:temp.end.line;
	temp.bpp.file <- input.data[temp.wanted.lines];
	write.table(temp.bpp.file, append = TRUE, file='L50_10.txt', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE);#change the file name here
}
