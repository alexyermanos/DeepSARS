library(Rsubread)
library(Rsamtools)
library(stringr)
library(ggplot2)
library(reshape2)
library(ShortRead)
library(seqinr) 


rm(list=ls())
wd="/Users/karenhung/Desktop/lssi/SARS-CoV-2-NGS/r32_clean/data"
setwd(wd)


#Barcode

# r_barcode_generate<-function(x){
#   bases=c("A","C","G","T")
#   xx=unlist(strsplit(toupper(x),NULL))
#   substring(reverse(
#   paste(unlist(lapply(xx,function(bbb){
#     if(bbb=="A") compString<-"T"
#     if(bbb=="C") compString<-"G"
#     if(bbb=="G") compString<-"C"
#     if(bbb=="T") compString<-"A"
#     if(!bbb %in% bases) compString<-"N"
#     return(compString)
#   })),collapse="")),2)
# }


f_barcodes = c("ATTAGATTAG","TCTACAACTG","ATGGCACTT","CCGCCAACAA","GGTTGTTGTT","TCCTGGTATT","ACCATGTTGT","CAACAATGGT"
               ,"CCAAGCCAGT","ATGTGGTTCT","CGAGTTGCGT","GACAGTAGCT")

#r_barcodes= sapply(f_barcodes, r_barcode_generate)

#Input Data file CSV
file=as.vector(read.csv("r32.csv",header = F)[,1])


#Seperate each patient barcode in to different fq file
barcode_extract=function(x){
  
  if(file.exists("fq")){
    print("fq folder exist")
  }else{
    dir.create(file.path(wd,"fq"))
  }
  
  
experiment=paste(unlist(strsplit(x,split="_"))[7],unlist(strsplit(x,split="_"))[8],sep="-")
seq_direct=unlist(strsplit(x,split="_"))[12]
reads=readFastq(x)
  
  if(seq_direct=="R1"){
    barcodes=f_barcodes
    i=c(110,151)
  }else if (seq_direct=="R2"){
    barcodes=r_barcodes
    i=c(1,15)
  }else{
    print("error")
  }

outputcsv=data.frame()
for(j in c(1:length(barcodes))) {
  seqs <- sread(reads) # get sequence list
  qual <- quality(reads) # get quality score list
  qual <- quality(qual) # strip quality score type
  seqs2 <- narrow(seqs, start=i[1], end=i[2])

  trimCoords=vmatchPattern(barcodes[j],seqs2,
                           max.mismatch=1, min.mismatch=0,
                           with.indels=FALSE, fixed=TRUE,
                           algorithm="auto")
  pattern=lengths(trimCoords)
  pattern=pattern>0
  seqs=seqs[(pattern)]
  qual=qual[(pattern)]
  
  qual <- SFastqQuality(qual)
  trimmed <- ShortReadQ(sread=seqs, quality=qual, id=id(reads[(pattern)]))
  outputFileName <- paste(experiment,"_",seq_direct,"_",f_barcodes[j], ".fq", sep="")
  writeFastq(trimmed, file.path(wd,"fq",outputFileName))
  print(paste("wrote", length(trimmed),"/",length(trimCoords),
              "reads to file", outputFileName))
  output=data.frame(outputFileName)
  outputcsv=rbind(outputcsv,output)
}
return(outputcsv)
}



barcode_result=sapply(file,barcode_extract)
fqfilelist=as.character(unlist(barcode_result)) 
write.csv(fqfilelist,"fqfilelist.csv")

#Map each fq file into reference
ref <- system.file("extdata","MN908947_master_RNAP-GAPDH.fasta",package="Rsubread")
buildindex(basename="my_index",reference="MN908947_master_RNAP-GAPDH.fasta")


if(file.exists("bam")){
  print("bam folder exist")
}else{
  dir.create(file.path(wd,"bam"))
}



for (i in c(1:length(fqfilelist))) {
  

  align(index="my_index",readfile1=file.path(wd,"fq",fqfilelist[i]),type="DNA",
        output_file=file.path(wd,"bam",paste0(gsub(".fq","",fqfilelist[i]),".bam")),nTrim5 = 30,nTrim3 = 30,minFragLength=10,maxFragLength=90)
  
}






