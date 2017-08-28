# This script reads in the data files and caluclates chromosomal intermingling degree (IMD) and protein enrichment

dira="insert path to datafiles folder in sample1"
setwd(dira)
binary_nucleus_ChrAChrB_S1 <- read.csv("binary_nucleus.csv")
binary_chromsome_A_ChrAChrB_S1 <- read.csv("binary_chromsome_2.csv")
binary_chromsome_B_ChrAChrB_S1 <- read.csv("binary_chromsome_3.csv")
binary_IMR_ChrAChrB_S1 <- read.csv("binary_IMR.csv")
pol2_IMR_ChrAChrB_S1 <- read.csv("pol2_IMR_binary.csv")
pol2_nucelus_ChrAChrB_S1 <- read.csv("pol2_nucelus_binary.csv")

dira="insert path to datafiles folder in sample2"
setwd(dira)
binary_nucleus_ChrAChrB_S2 <- read.csv("binary_nucleus.csv")
binary_chromsome_A_ChrAChrB_S2 <- read.csv("binary_chromsome_2.csv")
binary_chromsome_B_ChrAChrB_S2 <- read.csv("binary_chromsome_3.csv")
binary_IMR_ChrAChrB_S2 <- read.csv("binary_IMR.csv")
pol2_IMR_ChrAChrB_S2 <- read.csv("pol2_IMR_binary.csv")
pol2_nucelus_ChrAChrB_S2 <- read.csv("pol2_nucelus_binary.csv")

# insert names of the conditions 
cond<-c("S1","S2")

# The intermingling degree was calculated by dividing the volume of the IMR  between two chromosomes by the sum total volume of the two chromosomes.
IMD_S1<-binary_IMR_ChrAChrB_S1[,3]/(binary_chromsome_A_ChrAChrB_S1[,3]+binary_chromsome_B_ChrAChrB_S1[,3])
IMD_S2<-binary_IMR_ChrAChrB_S2[,3]/(binary_chromsome_A_ChrAChrB_S2[,3]+binary_chromsome_B_ChrAChrB_S2[,3])


#The enrichment of active RNAPII in the intermingling regions was obtained by dividing the mean intensity of active RNAPII in the IMR by the mean intensity of the active RNAPII in the entire nucleus.
Enrichment_pol2_S1<-pol2_IMR_ChrAChrB_S1[,4]/pol2_nucelus_ChrAChrB_S1[,4]
Enrichment_pol2_S2<-pol2_IMR_ChrAChrB_S2[,4]/pol2_nucelus_ChrAChrB_S2[,4]


# datavisualisation and statistical tests 

dira="insert path to R project folder"
setwd(dira)

png(filename=paste("IMD.png",sep=""), units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis = 2,font.lab=2,family="serif")
boxplot(IMD_S1,IMD_S2,las=1,col="gray",pch=18,cex=0.6,names=cond,lty=1,ylab="Intermingling Degree")
dev.off()

t.test(IMD_S1,IMD_S2)

a<-nrow(subset(binary_IMR_ChrAChrB_S1,binary_IMR_ChrAChrB_S1[,3]==0))/nrow(binary_IMR_ChrAChrB_S1)
b<-(nrow(binary_IMR_ChrAChrB_S1)-a)/nrow(binary_IMR_ChrAChrB_S1)
c<-nrow(subset(binary_IMR_ChrAChrB_S2,binary_IMR_ChrAChrB_S2[,3]==0))/nrow(binary_IMR_ChrAChrB_S2)
d<-(nrow(binary_IMR_ChrAChrB_S2)-c)/nrow(binary_IMR_ChrAChrB_S2)

Values <- matrix(c(a,b,c,d),nrow = 2,ncol = 2,byrow = FALSE)

png(filename=paste("percentage.png",sep=""), units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis = 2,font.lab=2,family="serif")
barplot(Values,cex.names=0.8,names=cond,ylab = "Percentage of nuclei",xlim=c(0,4),las=1,col=c("blue","red"),cex.axis=0.8)
box()
legend("topright",c("Intermingling","Separate"), cex = 0.7,fill=c("red","blue"))
dev.off()


png(filename=paste("Pol2.png",sep=""), units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis = 2,font.lab=2,family="serif")
boxplot(Enrichment_pol2_S1,Enrichment_pol2_S2,las=1,col="gray",pch=18,names=cond,lty=1,ylab="Active RNAPII Enrichment Factor",outline=T)
dev.off()



