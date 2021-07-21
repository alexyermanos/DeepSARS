library(Rsamtools)
library("data.table")


source('/Users/karenhung/Desktop/lssi/SARS-CoV-2-NGS/r32_clean/mapRead.R')
source('/Users/karenhung/Desktop/lssi/SARS-CoV-2-NGS/r32_clean/primer.R')

wd="/Users/karenhung/Desktop/lssi/SARS-CoV-2-NGS/r32_clean/data"
setwd(wd)


#Input bam file
bamFiles=read.csv("fqfilelist.csv")
bamFiles=gsub(".fq",".bam",bamFiles$x)
setwd(file.path(wd,"bam"))
H <- loadBam(bamFiles) 

setwd(wd)
bedFile <- import("MN908947_master_RNAP_GAPDH.bed", format="bed")


setwd("/Users/karenhung/Desktop/lssi/SARS-CoV-2-NGS/r32_clean/data/primer")
df=createPrimerInterval("pool of mutation_RNAP_GAPDH.bam") %>% getCountPerIntervalFast(H, .) 
df2=data.table()
for (i in c(1:length(df))){
if (length(df[[i]])!=0){
df3=data.frame(df[[i]],file=bamFiles[i])
df2=rbind(df2,df3)
}
}
df2=separate(df2,file,c("Exp","Seqdir","barcode"),sep="_")
df2$n.read[df2$n.read==0]=NA


#plot each index with different barcode together
exp=unique(df2$Exp)
for (i in c(1:6)) {
ggplot(data=df2[(df2$Exp==exp[i]),],aes(x=name,y=n.read,group=barcode,fill=barcode))+
  geom_col(stat="identity",position = position_dodge(preserve = "single"),width = 0.8,color="black")+
  geom_text(aes(label=n.read),position = position_dodge(0.8),vjust=-0.8)+
  scale_y_log10()+
  scale_alpha_manual(values=c(1,0.5))+
  facet_grid(Seqdir~.)+
  labs(title = exp[i], x = "primer", y = "number of read")
  ggsave(paste0(exp[i],".png"),width = 15,height = 10)
  }


#plot each index with different barcode seperate
df4=df2
df4$barcode_number[(df4$barcode%in%c("ATTAGATTAG.bam"))]="barcode01"
df4$barcode_number[(df4$barcode%in%c("TCTACAACTG.bam"))]="barcode02"
df4$barcode_number[(df4$barcode%in%c("ATGGCACTT.bam"))]="barcode03"
df4$barcode_number[(df4$barcode%in%c("CCGCCAACAA.bam"))]="barcode04"
df4$barcode_number[(df4$barcode%in%c("GGTTGTTGTT.bam"))]="barcode05"
df4$barcode_number[(df4$barcode%in%c("TCCTGGTATT.bam"))]="barcode06"
df4$barcode_number[(df4$barcode%in%c("ACCATGTTGT.bam"))]="barcode07"
df4$barcode_number[(df4$barcode%in%c("CAACAATGGT.bam"))]="barcode08"
df4$barcode_number[(df4$barcode%in%c("CCAAGCCAGT.bam"))]="barcode09"
df4$barcode_number[(df4$barcode%in%c("ATGTGGTTCT.bam"))]="barcode10"
df4$barcode_number[(df4$barcode%in%c("CGAGTTGCGT.bam"))]="barcode11"
df4$barcode_number[(df4$barcode%in%c("GACAGTAGCT.bam"))]="barcode12"


for (i in c(1:6)) {
  ggplot(data = df4[(df4$Exp==exp[i]),], mapping = aes(x = name, fill = name, y = n.read))+ 
    geom_col() +
    facet_wrap(~ `barcode_number`,ncol=5)+
    scale_y_log10()+
    geom_text(size=2.5,aes(label=n.read),position = position_dodge(0.8),vjust=-0.5)+
    theme(axis.text.x = element_text(size=9,angle = 90, hjust = 1, vjust = 1))+
    labs(title = exp[i])
  ggsave(paste0(exp[i],"_group.png"),width = 15,height = 5)
}



#plot single barcode (control)
exp3=subset(df4,Exp=="r32-3"&barcode_number=="barcode01")
exp4=subset(df4,Exp=="r32-4"&barcode_number=="barcode01")
exp5=subset(df4,Exp=="r32-5"&barcode_number=="barcode01")
exp6=subset(df4,Exp=="r32-6"&barcode_number=="barcode01")

exp=list(exp3,exp4,exp5,exp6)

exp[[1]]


for(i in c(1:length(exp))){
ggplot(data=exp[[i]],aes(x=name,y=n.read,group=barcode,fill=name))+
 geom_col(stat="identity",position = position_dodge(preserve = "single"),width = 0.8,color="black")+
   geom_text(size=5,aes(label=n.read),position = position_dodge(0.8),vjust=-0.8)+
   scale_y_log10()+
   scale_alpha_manual(values=c(1,0.5))+
   facet_grid(Seqdir~.)+
   theme(axis.text.x = element_text(size=18,angle = 45, hjust = 1, vjust = 1))+
   labs(title = unique(exp[[i]]$Exp), x = "primer", y = "number of read")
   ggsave(paste0(unique(exp[[i]]$Exp),"_s.png"),width = 15,height = 10)
}

  
 