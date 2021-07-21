library(tidyverse)
library(magrittr)
library(Rsamtools)
library(tidyr)
library(ggplot2)
library(data.table)
library(compare)
library(gridExtra)



###############################   COUNT THE NUMBER OF READS PER GENE ################


# input: list containing the path/names of the .bam files
# output: a list of dataframe
# put the bam in a list of dataframes and removes unaligned reads from each file
loadBam <- function(files_path){
  # load the bam file
  BAM <- list()
  for (i in seq_along(files_path)){
    BAM[i] <-scanBam(files_path[[i]])
  }
  # put the bam files into dataframes and remove the unaligned reads
  BAM.df <- list()
  for (b in seq_along(BAM)){
    BAM.df[[b]] <- data.frame(BAM[[b]])
    BAM.df[[b]] <- BAM.df[[b]] %>% drop_na(pos)
  }
  
  return(BAM.df)
}


# input: bedfile
# output: a dataframe with each gene on a row
# create a dataframe to count the number of reads per gene for CoV-2 with RNAP (name, n.read, start, end)
createGeneInterval <- function(bed){
  
  #get bed as df
  bed <- as.data.frame(bed)
  bed
  #get only genes
  bed_gene <-bed[grepl("gene", bed$name),]
  bed_gene
  
  #create new df with gene names, n.reads, start and end position
  reads_per_genes <- data.frame("name" = bed_gene$name, n.read = 0)
  reads_per_genes$start = bed_gene$start
  reads_per_genes$end = bed_gene$end
  
  reads_per_genes$name=as.character(reads_per_genes$name)
  #reads_per_genes <- rbind(reads_per_genes, list("RNAP",0,427,2758))# add rnap
  #reads_per_genes <- rbind(reads_per_genes, list("GAPDH",0,1,426))# add GAPDH
  reads_per_genes %<>% filter(name != "orf1ab gene")# split gene orf1ab
  reads_per_genes <- rbind(reads_per_genes, list("ORF1a", 0,3024,16229))
  reads_per_genes <- rbind(reads_per_genes, list("ORF1b", 0,16230,24313))
  
  #reads_per_genes <- rbind(reads_per_genes, list("RNAP", 0,1,2332)) # add rnap
  #reads_per_genes %<>% filter(name != "orf1ab gene")# split gene orf1ab
  #reads_per_genes <- rbind(reads_per_genes, list("ORF1A", 0,2598,15803))
  #reads_per_genes <- rbind(reads_per_genes, list("ORF1B", 0,15804,23887))
  
  reads_per_genes$name <- sub(" .*","", reads_per_genes$name)
  reads_per_genes$name <- factor(reads_per_genes$name, levels = reads_per_genes$name[order(reads_per_genes$start)])
  
  reads_per_genes<- reads_per_genes[order(reads_per_genes$start),]
  return(reads_per_genes)
}


# # input: list of bam files, a dataframe containg the different intervals and their start and end position: (name, n.read, start, end)
# # output: list of dataframes with the number of reads per interval
# # count the number of reads per interval
# #TODO: optimize
# getCountPerInterval <- function(bam, intervals){
#   
#   df <- list()
#   for (i in seq_along(bam)) {
#     # loop through the bam position, add a read in the right gene
#     print(paste0("bam: ", i))
#     intervals$n.read <- 0
#     for(pos in bam[[i]]$pos){
#       if(!is.na(pos)){
#         for(r in 1:nrow(intervals)) {
#           row <- intervals[r,]
#           if(pos >= row$start && pos < row$end){
#             intervals[r,2] <- intervals[r,2] +1 #TODO use name
#           }
#         }
#       }
#     }
#     df[[i]] <- intervals
#   }
#   return(df)
# }

# input: list of bam files, a dataframe containg the different intervals and their start and end position: (name, n.read, start, end)
# output: list of dataframes with the number of reads per interval
# count the number of reads per interval
getCountPerIntervalFast <- function(bam, intervals){
  
  inter <- data.table(intervals)
  df <- list()
  for (i in seq_along(bam)) {
    tryCatch({
    # loop through the bam position, add a read in the right gene
    # print(paste0("bam: ", i))
    inter$n.read <- 0
    # create a table with read position and number of occurence
    t <- table(bam[[i]]$pos)
    t <- data.frame(t, stringsAsFactors=F)
    t %<>% mutate(Position = strtoi(Var1, base = 0L))
    for (pos in t$Position){
      # for each row of interval, if pos is in the interval, add the number of occurence
      for(r in 1:nrow(inter)) {
        row <- inter[r,]
        if(pos >= row$start && pos < row$end){
          inter[r,2] <- inter[r,2] + t[t$Position==pos,]$Freq
        }
      }
    }
    df[[i]] <- data.frame(inter)
    }, error=function(e){cat("Number",i,"ERROR :",conditionMessage(e), "\n")})
    }
  return(df)
}
 

# input: a dataframe with the number of reads (n.read) per gene (name)
# plot the number of reads per gene for each file
plotReadPerInterval <- function(df){

  list_plot <- list()
  for (i in seq_along(df)) {
    p <- ggplot(df[[i]], aes(x = name, y = n.read)) + geom_col() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste0("Number of Read, file ", i), x = "gene", y = "number of read") +
      theme(plot.title = element_text(hjust = 0.5)) + scale_y_log10()
    list_plot[[i]] <- p
  }
  do.call(grid.arrange,list_plot)
}

# input: a dataframe with the number of reads (n.read) per gene (name), and a title for the plot
# plot the number of reads per gene for each file and save the plots
plotReadPerIntervalAndSave <- function(df, title){
  setwd("~/lssi/SARS-CoV-2-NGS/r32_clean/plot/gene")
  
  list_plot <- list()
  for (i in seq_along(df)) {
    p <- ggplot(df[[i]], aes(x = name, y = n.read)) + geom_col() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste0(), x = "gene", y = "number of read") +
      theme(plot.title = element_text(hjust = 0.5)) + scale_y_log10()
    ggsave(paste0(gsub(" ","",title),"_file", i,".png"))
    list_plot[[i]] <- p
  }
  do.call(grid.arrange,list_plot)
}


# setwd("~/lssi/SARS-CoV-2-NGS/data")
# bamFiles <- c("H1.bam","H2.bam","H3.bam","H4.bam","H5.bam","H6.bam","H7.bam","H8.bam","H9.bam","H10.bam","H11.bam","H12.bam")
# bedFile <- import("Reference_RNAP_COV.bed", format="bed")
# 
# H <- loadBam(bamFiles)
# bedFile %>% createGeneInterval() %>% getCountPerIntervalFast(H,.) %>% plotReadPerInterval()

