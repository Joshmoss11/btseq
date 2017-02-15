rev.comp<-function(x,rev=TRUE)
{
x<-toupper(x)
y<-rep("N",nchar(x))
xx<-unlist(strsplit(x,NULL))
for (bbb in 1:nchar(x))
    {
        if(xx[bbb]=="A") y[bbb]<-"T"    
        if(xx[bbb]=="C") y[bbb]<-"G"    
        if(xx[bbb]=="G") y[bbb]<-"C"    
        if(xx[bbb]=="T") y[bbb]<-"A"
    }
if(rev==FALSE) 
    {
    for(ccc in (1:nchar(x)))
        {
        if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
        }
    }
if(rev==T)
    {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
        {
        zz[ccc]<-y[nchar(x)+1-ccc]
        if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
        }
    }
    return(yy)  
}


d <- getwd()
ss <- read.table(file.path(d,'sample_sheet.txt'),sep='\t',header=F,stringsAsFactors=F)
#library(Biostrings)
#rev_comp <- function(x){ return(as.character(reverseComplement(DNAString(x))))}
bc <- sapply(ss[,3],rev.comp)
bc.unique <- unique(bc)
bc.n <- length(bc.unique)
#sequences <- paste(ss[,4],'NN')
#sequences.id <- sapply(ss[,1],function(x) {strsplit(x,'-')[[1]][1]})

#seq.pair <- cbind(sequences.id,sequences)

## write fasta files
for (i in 1:nrow(ss)){
	#fileConn<-file(paste0(ss[i,1],".fa"),open='w')
	fileConn<-file(paste0('sample',i,".fa"),open='w')
	writeLines(c(paste0(">",ss[i,1]),paste(c(rep('N',20),ss[i,4],rep('N',20)),sep='',collapse='')), fileConn)
	close(fileConn)
}

## write sample sheet for bcl2fastq
bc.name <- paste0('A00',1:bc.n)
if(bc.n>9){
	for (i in 10:bc.n){
		bc.name[i] <- paste0('A0',i)
	}
}

bc_ss <- cbind(1:bc.n,bc.name,bc.unique)
colnames(bc_ss) <- c('Sample_ID','Sample_Name','index')
fileConn<-file('bc_ss.csv',open='w')
writeLines('[Data]', fileConn)
writeLines(paste(colnames(bc_ss),collapse=','),fileConn)
write.table(bc_ss,fileConn,append=T,sep=',',row.names=F,quote=F,col.names=F)
close(fileConn)


# write input file for BSTarget analysis
bc.id <- match(bc,bc.unique)
#bs_input <- cbind(1:nrow(ss),paste(ss[1,],ss[,2],sep='_'),file.path(d,'trim_galore',paste0(

fileConn<-file('BSTarget_input.txt',open='w')
for (i in 1:nrow(ss)){
	#writeLines(paste(i,paste(ss[i,1],ss[i,2],sep='_'),file.path(d,'trim_galore',paste0(bc.name[bc.id[i]],'_S',bc.id[i],'_R1_001_trimmed.fq.gz')),file.path(d,paste0('sample',i,".fa")),sep='\t'), fileConn)
	    writeLines(paste(i,'Sample',file.path(d,'trim_galore',paste0(bc.name[bc.id[i]],'_S',bc.id[i],'_R1_001_trimmed.fq.gz')),file.path(d,paste0('sample',i,".fa")),sep='\t'), fileConn)
}
close(fileConn)
