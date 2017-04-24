setwd("/bi/group/reik/Nelly/Degradation/New_relative_analysis/")
read.delim("Probe_report_methylation_PBAT_vs_BS-seq.txt") -> meth.data
head(meth.data)

par(mar=c(5,8,4,2))
boxplot(meth.data[,13:ncol(meth.data)],horizontal = TRUE,las=1)

bs.diffs <- (meth.data$Sperm.BS.seq - meth.data$ESC.BS.seq)-(meth.data$Sperm.PBAT - meth.data$ESC.PBAT)

map.colours <- function (value,high.low,palette) {
  
  value[value<high.low[1]] <- high.low[1]
  value[value>high.low[2]] <- high.low[2]
  
  proportion <- ((value-high.low[1])/(high.low[2]-high.low[1]))
  index <- round ((length(palette)-1)*proportion)+1
  return (palette[index])
}

colorRampPalette(c("red3","grey","mediumblue"))(100) -> my.palette

range(na.omit(bs.diffs))

diff.colours <- map.colours(bs.diffs,c(-30,30),my.palette)

plot((meth.data$Sperm.PBAT+meth.data$ESC.PBAT)/2,
     meth.data$Sperm.PBAT-meth.data$ESC.PBAT,
     pch=19,
     cex=0.5,
     col=diff.colours,
     xlab="Mean methylation level between Sperm and ESC",
     ylab="Methylation difference (Sperm - ESC) in PBAT")

legend("topleft",
       c("BS-Seq diff 30% larger than PBAT","BS-Seq diff same as PBAT","BS-Seq diff 30% smaller than PBAT"),
       fill=c("red3","grey","mediumblue"))
