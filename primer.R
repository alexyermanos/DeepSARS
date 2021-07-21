library(Rsamtools)
library(tidyr)
library(rtracklayer)


source("/Users/karenhung/Desktop/lssi/SARS-CoV-2-NGS/r32_clean/mapRead.R")

setwd("/Users/karenhung/Desktop/lssi/SARS-CoV-2-NGS/r32_clean/data/primer")

############################    PLOT READS PER PRIMER   #######################

# input: bam file of primers
# output: df with a row per primer: (name, n.read, start, end)
# create a dataframe to count the number of reads per primer
createPrimerInterval <- function(primer.bam){
  Primer <- loadBam(primer.bam)
  primer.interval <- data.frame(name = Primer[[1]]$qname, n.read = 0, start = Primer[[1]]$pos,
                                end = Primer[[1]]$pos + 150)
  # primer.interval$name <- sub("P_PCR_","", primer.interval$name)
  primer.interval$name <- sub("extraction","ext", primer.interval$name)
  primer.interval$name <- sub("PCR_","", primer.interval$name)
  primer.interval$name <- factor(primer.interval$name, levels=primer.interval$name)
  return(primer.interval)
}


# input: list of dfs
# plot, for each file, the number of reads per primer
plotReadPerPrimer <- function(df){
  
  setwd("/Users/karenhung/Desktop/lssi/SARS-CoV-2-NGS/r32_clean/plot/primer")
  list_plot <- list()
  for (i in seq_along(df)) {
    p <- ggplot(df[[i]], aes(x = name, y = n.read)) + geom_col() +
      labs(title = paste0("Read per primer, file ", i), x = "primer", y = "number of read") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(plot.title = element_text(hjust = 0.5)) + scale_y_log10()+
      geom_text(aes(label=n.read,alpha=as.factor(n.read)),vjust=-0.3)+
      scale_alpha_manual(values=c(0,rep(1,length(df[[i]]$name-1))))+ 
      theme(legend.position = "none")
    ggsave(paste0("ReadPerPrimer_file", i,".png"))
    list_plot[[i]] <- p
  }
  do.call(grid.arrange,list_plot)
  
}
plotReadPerPrimer2 <- function(df){
  
  setwd("/Users/karenhung/Desktop/lssi/SARS-CoV-2-NGS/r32_clean/plot/primer")
  list_plot <- list()
  for (j in seq_along(df)) {
    p <- ggplot(df[[j]], aes(x = name, y = n.read)) + geom_col() +
      labs(title = paste0("Read per primer, file ", i), x = "primer", y = "number of read") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(plot.title = element_text(hjust = 0.5)) + scale_y_log10()+
      geom_text(aes(label=n.read,alpha=as.factor(n.read)),vjust=-0.3)+
      scale_alpha_manual(values=c(0,rep(1,length(df[[j]]$name-1))))+ 
      theme(legend.position = "none")
    ggsave(paste0("ReadPerPrimer_file", i,".png"))
    list_plot[[j]] <- p
  }
  do.call(grid.arrange,list_plot)
  
}

# createPrimerInterval("Primer_mapping.bam") %>% getCountPerIntervalFast(H, .) %>% plotReadPerPrimer()

###### comparison number of reads map to primer (diff: no primer for RNAP, multiple primer with same pos) ######

# input: list of reads df and primer count df
# output: df with the ratio of reads
# compare the number of read map to primer vs the total number of read
quantifyPrimerRead <- function(bam, primerInterval) {
  n.reads.primer <- c()
  n.reads.h <- c()
  for(i in seq_along(primerInterval)){
    n.reads.primer[i] <- sum(primerInterval[[i]]$n.read)
    n.reads.h[i] <- length(bam[[i]]$pos)
  }
  n.reads.primer
  n.reads.h
  quantify.n.primer.read = data.frame(n.read = n.reads.h, n.read.primer = n.reads.primer)
  quantify.n.primer.read %<>% mutate(ratio = n.read.primer / n.read)
  return(quantify.n.primer.read)
}
# h.primer.count <- createPrimerInterval("Primer_mapping.bam") %>% getCountPerIntervalFast(H, .)
# quantify.n.primer.read <- quantifyPrimerRead(H, h.primer.count)


############# CREATE CSV WITH PRIMER READS PER SAMPLE ###################

setwd("/Users/karenhung/Desktop/lssi/SARS-CoV-2-NGS/r32_clean/data/primer")


# input: bam primer file, list of df with the number of reads per primer
# output: df with the number of reads per primer for each sample
# create a single df with the primer name, sequence, and reads for each sample
createPrimerDf <- function(primerFile, primerCountDf){
  Primer <- loadBam(primerFile)
  primer.df = data.frame(name =Primer[[1]]$qname, sequence = Primer[[1]]$seq)
  
  for(i in seq_along(primerCountDf)){
    primer.df$col = primerCountDf[[i]]$n.read
    names(primer.df)[i+2] <- unlist(mapply(function(x,y) paste(x, y, sep="."), "sample", i))
  }
  return(primer.df)
}
# p.df <- createPrimerDf("Primer_mapping.bam", h.primer.count)
# 
# write.csv(p.df,"primer_map.csv", row.names=FALSE)

