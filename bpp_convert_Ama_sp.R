#from unlinked SNPs to bpp inputs
data.loci <- scan('GamaAllSP.loci', what = 'character', sep = '\n');#used the edited .loci file
species.vec <- c('GecuC', 'GecuP', 'GecuE', 'GlicC', 'GlicP', 'GlicE', 'GmorB', 'GpasB');#enter here the taxa


break.lines <- grep('//', data.loci);
index.lines <- c(0, break.lines);
loci.count <- 0;#number of loci counter. will print the number of loci at the end

for(ii in 1:956){
	
	num.1 <- index.lines[ii] + 1;
	num.2 <- index.lines[ii + 1] - 1;
	
	temp.text <- data.loci[c(num.1:num.2)];

	num.sp <- 0;
	ind.vec <- NULL;
	for(iter in 1:length(species.vec)){
		temp.ind.vec <- grep(species.vec[iter], temp.text);
		ind.vec <- c(ind.vec, temp.ind.vec);
		num.sp <- num.sp + (length(temp.ind.vec) > 0);
		
		if(num.sp == length(species.vec)){
			loci.count <- loci.count + 1;
			
			subset.temp <- temp.text[ind.vec];
			
			string.ed <- subset.temp[1];
			find.p <- '>GecuC_[[:alnum:]]+[[:space:]]+'#for cal length of the loci; data here are all 110
			sub.pa <- '';
			temp.dna.seq <- sub(find.p,sub.pa,string.ed);
			dna.seq <- strsplit(temp.dna.seq, '');
			loci.length <- length(dna.seq[[1]])
			n.ind <- length(subset.temp);
			first.line <- paste(n.ind, loci.length, sep = ' ');#first line for bpp each locus
			
			new.string <- NULL;
			Gblu.vec <- grep(species.vec[1], subset.temp);
			count <- 1;
			for(iii in Gblu.vec){
				s.pattern <- '>GecuC_[[:alnum:]]+';
				bpp.tail <- paste('^', count, sep = '');
				rep.pattern <- paste('	GecuC', bpp.tail, sep = '');
				new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
				count <- count + 1;
			}
			GherD.vec <- grep(species.vec[2], subset.temp);
			count <- 1;
			for(iii in GherD.vec){
				s.pattern <- '>GecuP_[[:alnum:]]+';
				bpp.tail <- paste('^', count+3, sep = '');#count + 3 because there are 3 GecuC samples
				rep.pattern <- paste('	GecuP', bpp.tail, sep = '');
				new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
				count <- count + 1;
			}			
			GherG.vec <- grep(species.vec[3], subset.temp);
			count <- 1;
			for(iii in GherG.vec){
				s.pattern <- '>GecuE_[[:alnum:]]+';
				bpp.tail <- paste('^', count+7, sep = '');
				rep.pattern <- paste('	GecuE', bpp.tail, sep = '');
				new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
				count <- count + 1;
			}							
			GreiM.vec <- grep(species.vec[4], subset.temp);
			count <- 1;
			for(iii in GreiM.vec){
				s.pattern <- '>GlicC_[[:alnum:]]+';
				bpp.tail <- paste('^', count+14, sep = '');
				rep.pattern <- paste('	GlicC', bpp.tail, sep = '');
				new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
				count <- count + 1;
			}
			GreiS.vec <- grep(species.vec[5], subset.temp);
			count <- 1;
			for(iii in GreiS.vec){
				s.pattern <- '>GlicP_[[:alnum:]]+';
				bpp.tail <- paste('^', count+16, sep = '');
				rep.pattern <- paste('	GlicP', bpp.tail, sep = '');
				new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
				count <- count + 1;
			}							
			Gtri.vec <- grep(species.vec[6], subset.temp);
			count <- 1;
			for(iii in Gtri.vec){
				s.pattern <- '>GlicE_[[:alnum:]]+';
				bpp.tail <- paste('^', count+19, sep = '');#two samples of Grei
				rep.pattern <- paste('	GlicE', bpp.tail, sep = '');
				new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
				count <- count + 1;
			}
			GreiS.vec <- grep(species.vec[7], subset.temp);
			count <- 1;
			for(iii in GreiS.vec){
				s.pattern <- '>GmorB_[[:alnum:]]+';
				bpp.tail <- paste('^', count+26, sep = '');
				rep.pattern <- paste('	GmorB', bpp.tail, sep = '');
				new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
				count <- count + 1;
			}							
			Gtri.vec <- grep(species.vec[8], subset.temp);
			count <- 1;
			for(iii in Gtri.vec){
				s.pattern <- '>GpasB_[[:alnum:]]+';
				bpp.tail <- paste('^', count+27, sep = '');#two samples of Grei
				rep.pattern <- paste('	GpasB', bpp.tail, sep = '');
				new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
				count <- count + 1;
			}
			
			
							
			temp.bpp.file <- c(first.line, new.string, '\n');
			write.table(temp.bpp.file, append = TRUE, file='Dy.AmaSP_all_loci.txt', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE);
			
		}
	}
	
}
cat('a total of', loci.count, 'loci');

#a total of 956 loci