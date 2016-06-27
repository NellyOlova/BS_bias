setwd("/bi/group/bioinf/Nelly_Olova/GC_Bias/")

window.size <- "20kb"

read.delim(paste("FOR_reads_count_",window.size,"_filtered_composition.txt",sep="")) -> for.counts

# Remove windows with N in
for.counts[for.counts$N == 0,] -> for.counts

# Make a %G and %C value
for.counts$G <- 100 * ((for.counts$G) / ((for.counts$End - for.counts$Start)+1))
for.counts$C <- 100 * ((for.counts$C) / ((for.counts$End - for.counts$Start)+1))

# Get just the count and G+C data
for.counts[,c(13:(ncol(for.counts)-5),ncol(for.counts)-4,ncol(for.counts)-1)] -> all.counts

# Remove any samples which don't contain any data
all.counts[,colSums(all.counts) != 0] -> all.counts

gsub("lane._","",colnames(all.counts)) -> colnames(all.counts)
gsub("_L00.*","",colnames(all.counts)) -> colnames(all.counts)

colSums(all.counts)/1000000 -> total.counts

as.data.frame(t(apply(all.counts, 1, function(x) x/total.counts))) -> all.counts

all.counts$G <- for.counts$G
all.counts$C <- for.counts$C


plot.relationship <- function(letter) {
  # Collate the counts for each methylation level
  by(all.counts[,1:(ncol(all.counts)-2)],as.integer(all.counts[[letter]]),function(x){
    if (nrow(x) > 100) {
      return(colMeans(x))
    }
    return (NA)
  }) -> split.all.sums
  
  do.call(rbind,split.all.sums) -> split.all.sums
  
  split.all.sums[!is.na(split.all.sums[,1]),] -> split.all.sums
  
  plot(density(as.integer(all.counts[[letter]]),adjust=3))
  
  library(RColorBrewer)
  
  brewer.pal(9,"Set1") -> colours
  
  plot(NA,
       xlim=c(0,100),
       ylim=c(0,max(split.all.sums)),
       xlab=paste("%",letter),
       ylab="Coverage")
  
  sapply(1:ncol(split.all.sums),
         function(x) {
           lines(as.integer(rownames(split.all.sums)),
                 split.all.sums[,x],
                 col=colours[x%%length(colours)]
           )
         }
  )
  
  legend("topright",colnames(split.all.sums),fill=colours)
}

plot.relationship("C")
plot.relationship("G")