args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1){
	stop("One argument must be supplied")
}
rem <- as.integer(args[1])
get_counts <- function(i,rem=NA){
	d <- file.path(paste0('Sample_',i),'Output')
	h.f <- file.path(d,list.files(d)[grepl("(?=.*CpG)(?=.*hist)",list.files(d),perl=T)])
	h <- read.table(h.f,sep='\t',skip=1,header=T,check.names=F,comment.char = "")
	seq <- readLines(paste0('sample',i,'.fa'))[2]
	n <- length(gregexpr('CG',seq)[[1]])
	h <- h[h[,3]+h[,5]==n,]
	if(!is.na(rem)){
		n <- n-1
		remT <- h[,rem+6] =='T'
		remC <- h[,rem+6] =='C'
		h[remC,3] <- h[remC,3]-1
		h[remT,5] <- h[remT,5]-1
	}
	s <- integer(n+1)
	t <- n:0
	for (i in 1:(n+1)){
		s[i] <- sum(h[h[,5]==t[i],1])
	}
	s <- c(sum(h[,1]),s)
	return(s)
}

ss <- read.table('sample_sheet.txt',sep='\t',header=F,stringsAsFactors=F)
n.s <- nrow(ss)
S <- sapply(1:n.s,get_counts,rem)
if (class(S)=="matrix"){
	S <- as.list(data.frame(S))
}
# print to file
N <- unlist(lapply(S,length))-2
con <- file('Summary_results_select.txt',open='w')
title <- paste(c('Sample #','Gene','Sample','CpGs',
		'All reads','All T',paste0('All T - ',1:max(N))),collapse='\t')
writeLines(title,con)
for (i in 1:n.s){
	line <- paste(c(i,ss[i,1],ss[i,2],N[i],S[[i]]),collapse='\t')
	writeLines(line,con)
}
close(con)
