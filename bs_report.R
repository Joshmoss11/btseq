bs_report <- function(summary_file){
  f <- summary_file
  x.l <- readLines(f); x.l1 <- which(x.l=="All C sites")
  x  <- read.table(f,sep='\t',skip=x.l1,header=T)
  CH <- which(x[,2]!='CG')
  CH.me <- 100-round(mean(x[CH,5]),1)
  return(CH.me)
}

sf <- '~/Downloads/19_6_19_2/summaries/Sample_31.summary'
print(bs_report(sf))
