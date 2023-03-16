library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)
library(SummarizedExperiment)
library(mzR)
library(MSnbase)
library(Cairo)
library(readr)

###############################################################################
 # Exponential phase
###############################################################################

# Load files
Directory <- "W:/users/Abraham/Exp_004_version_2/Data/Exponential"

setwd(Directory)

Directory <- "W:/users/Abraham/Exp_004_version_2/Data/Exponential/"

files<-list.files()


files_CID <- grep("CID",files,value=T)
files_EAD <- grep("EAD",files,value=T)

samples_CID <- grep("Exp",files_CID,value = T)
blank_CID <- grep("Exp",files_CID,value = T,invert=T)


samples_EAD <- grep("Exp",files_EAD,value = T)
blank_EAD <- grep("Exp",files_EAD,value = T,invert=T)


s_CID <- sort(rep(paste0("CID ",c("WT",1:15,"Con")),each=3 ))
s__CID <- paste0(s_CID,"_",c(1:3))

s_EAD <- sort(rep(paste0("EAD ",c("WT",1:15,"Con")),each=3 ))
s__EAD <- paste0(s_EAD,"_",c(1:3))

b_EAD <- rep("Blank EAD",5)
b__EAD <- paste0(b_EAD," ",c(1:5))
b_CID <- rep("Blank CID",6)
b__CID <- paste0(b_CID," ",c(1:6))


sample_ID <-c(s__CID,s__EAD,b__CID,b__EAD)
condition <- c(s_CID,s_EAD,b_CID,b_EAD)
metadata <- data.frame(Samples=sample_ID,Condition=condition)

files <- c(samples_CID,samples_EAD,blank_CID,blank_EAD)

files <- paste0(Directory,files)

data<-vector(mode="list",length = length(files))
names(data) <-  sample_ID
for(i in 1:length(files)){
  data[[i]] <- readMSData(files = files[i] ,mode="onDisk",verbose=FALSE) #openMSfile(files[1])  
}


scan<-table(msLevel(data1))

scans <- matrix(ncol=2,nrow=113)
colnames(scans) <- c("ms1","ms2")

for(i in 1:length(files)){
scans[i,] <-unlist(table(msLevel(data[[i]])))
}

color<-c(rep(c("coral1","palegreen3"),each=51),rep("grey",11))

label<-c(paste0("CID ",c(1:15,"Con","WT")),paste0("EAD ",c(1:15,"Con","WT")),"CID blanks","EAD blanks")
sequence <- seq(1,102,3)+3
space<-rep(0.01,length(files))
space[sequence]<-0.8
space[103+6]<-0.8

sequence <- c(seq(1,102,3)+1,106,110)
name<-rep(NA,length(files))
name[sequence] <- label

SAVE<-"W:/users/Abraham/Exp_004_version_2/Plots/1) Raw data plots/"
Cairo(file=paste0(SAVE,"barplot_EP_ms1.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
par(mar=c(9,6,5,5),mfrow=c(1,1))
barplot(scans[,1],col = color, main="Total Number of MS1 Scans per Sample (EP)",ylab="# of scans",
        space=space,names.arg =name ,las=3,cex.names  = 1)
dev.off()


Cairo(file=paste0(SAVE,"barplot_EP_ms2.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
par(mar=c(9,6,5,5),mfrow=c(1,1))
barplot(scans[,2],col = color, main="Total Number of MS2 Scans per Sample (EP)",ylab="# of scans",
        space=space,names.arg =name ,las=3,cex.names  = 1)
dev.off()








CID<-c(1:45,49:51)
EAD<-c(52:96,100:102)

Cairo(file=paste0(SAVE,"boxplot_EP_ms1_2.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
par(cex.main=1.5,cex.lab=1.6,cex.axis=1.4,mar=c(5,6,5,2),mfrow=c(1,2))
boxplot(cbind(scans[CID,1],scans[EAD,1]),col=c("coral1","palegreen3"),names=c("CID","EAD"),
        main="Total Number of MS1 Scans/Sample (EP)",ylab="# of scans",xaxt="n")
axis(1,at=c(1,2),labels = c("CID","EAD"),cex.axis=2)
boxplot(cbind(scans[CID,2],scans[EAD,2]),col=c("coral1","palegreen3"),names=c("CID","EAD"),
        main="Total Number of MS2 Scans/Sample (EP)",ylab="# of scans",xaxt="n")
axis(1,at=c(1,2),labels = c("CID","EAD"),cex.axis=2)
dev.off()


mean(scans[EAD,1])
mean(scans[EAD,2])

mean(scans[CID,1])
mean(scans[CID,2])


rt=c(5.714-0.2,5.714+0.1)*60
mz=c(704.521-0.01,704.521+0.01)
CID <- data[[51]]
EAD <- data[[102]]




Cairo(file=paste0(SAVE,"XIC_EAD_CID_PE_33_1.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
par(mfrow=c(1,2),cex.lab=1.5)
chr_raw <- MSnbase::chromatogram(CID, mz = mz, rt = rt)
par(mar=c(5,5,5,2))
plot(chr_raw,type="o",pch=16,cex=2,lwd=3,col="coral1",main="CID: PE 33:1  (mz = 704.511-704.531)")
text(333,1445000,labels="XIC",cex=2)



chr_raw <- MSnbase::chromatogram(EAD, mz = mz, rt = rt)
plot(chr_raw,type="o",pch=16,cex=2,lwd=3,col="palegreen3",main="EAD: PE 33:1  (mz = 704.511-704.531)")
text(333,1288000,labels="XIC",cex=2)
dev.off()



###############################################################################
# Stationary phase
###############################################################################

# Load files
Directory <- "W:/users/Abraham/Exp_004_version_2/Data/Stationary"

setwd(Directory)

Directory <- "W:/users/Abraham/Exp_004_version_2/Data/Stationary/"

files<-list.files()


files_CID <- grep("CID",files,value=T)
files_EAD <- grep("EAD",files,value=T)

samples_CID <- grep("Exp",files_CID,value = T)
blank_CID <- grep("Exp",files_CID,value = T,invert=T)


samples_EAD <- grep("Exp",files_EAD,value = T)
blank_EAD <- grep("Exp",files_EAD,value = T,invert=T)


s_CID <- sort(rep(paste0("CID ",c("WT",1:15,"Con")),each=3 ))
s__CID <- paste0(s_CID,"_",c(1:3))

s_EAD <- sort(rep(paste0("EAD ",c("WT",1:15,"Con")),each=3 ))
s__EAD <- paste0(s_EAD,"_",c(1:3))

b_EAD <- rep("Blank EAD",6)
b__EAD <- paste0(b_EAD," ",c(1:6))
b_CID <- rep("Blank CID",6)
b__CID <- paste0(b_CID," ",c(1:6))


sample_ID <-c(s__CID,s__EAD,b__CID,b__EAD)
condition <- c(s_CID,s_EAD,b_CID,b_EAD)
metadata <- data.frame(Samples=sample_ID,Condition=condition)

files <- c(samples_CID,samples_EAD,blank_CID,blank_EAD)

files <- paste0(Directory,files)

data<-vector(mode="list",length = length(files))
names(data) <-  sample_ID
for(i in 1:length(files)){
    data[[i]] <- readMSData(files = files[i] ,mode="onDisk",verbose=FALSE) #openMSfile(files[1])  
}

table(msLevel(data[[2]]))

scans <- matrix(ncol=2,nrow=114)
colnames(scans) <- c("ms1","ms2")

for(i in 1:length(files)){
  scans[i,] <-unlist(table(msLevel(data[[i]])))
}

color<-c(rep(c("coral1","palegreen3"),each=51),rep("grey",12))

label<-c(paste0("CID ",c(1:15,"Con","WT")),paste0("EAD ",c(1:15,"Con","WT")),"CID blanks","EAD blanks")
sequence <- seq(1,102,3)+3
space<-rep(0.01,length(files))
space[sequence]<-0.8
space[103+6]<-0.8

sequence <- c(seq(1,102,3)+1,106,112)
name<-rep(NA,length(files))
name[sequence] <- label

SAVE<-"W:/users/Abraham/Exp_004_version_2/Plots/1) Raw data plots/"
Cairo(file=paste0(SAVE,"barplot_SP_ms1.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
par(mar=c(9,6,5,5),mfrow=c(1,1))
barplot(scans[,1],col = color, main="Total Number of MS1 Scans per Sample (SP)",ylab="# of scans",
        space=space,names.arg =name ,las=3,cex.names  = 1)
dev.off()


Cairo(file=paste0(SAVE,"barplot_SP_ms2.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
par(mar=c(9,6,5,5),mfrow=c(1,1))
barplot(scans[,2],col = color, main="Total Number of MS2 Scans per Sample (SP)",ylab="# of scans",
        space=space,names.arg =name ,las=3,cex.names  = 1)
dev.off()










CID<-c(1:45,49:51)
EAD<-c(52:96,100:102)

Cairo(file=paste0(SAVE,"boxplot_SP_ms1_2.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
par(cex.main=1.5,cex.lab=1.6,cex.axis=1.4,mar=c(5,6,5,2),mfrow=c(1,2))
boxplot(cbind(scans[CID,1],scans[EAD,1]),col=c("coral1","palegreen3"),names=c("CID","EAD"),
        main="Total Number of MS1 Scans/Sample (SP)",ylab="# of scans",xaxt="n")
axis(1,at=c(1,2),labels = c("CID","EAD"),cex.axis=2)
boxplot(cbind(scans[CID,2],scans[EAD,2]),col=c("coral1","palegreen3"),names=c("CID","EAD"),
        main="Total Number of MS2 Scans/Sample (SP)",ylab="# of scans",xaxt="n")
axis(1,at=c(1,2),labels = c("CID","EAD"),cex.axis=2)
dev.off()


mean(scans[EAD,1])
mean(scans[EAD,2])

mean(scans[CID,1])
mean(scans[CID,2])





EAD<- read_table("C:/Users/amoyal/Downloads/EAD_PE_33_1_Na.msp", 
                             na = "empty")
EAD<-EAD[-c(1:9),-3]
colnames(EAD)<- c("mz","intensity")

EAD$intensity<-as.numeric(EAD$intensity)/max(as.numeric(EAD$intensity))*100

CID<-read_table("C:/Users/amoyal/Downloads/CID-PE_33_1_Na.msp", 
                na = "empty")
CID<-CID[-c(1:9),-3]
colnames(CID)<- c("mz","intensity")
CID$intensity<-as.numeric(CID$intensity)/max(as.numeric(CID$intensity))*100


Cairo(file=paste0(SAVE,"MS2_spectra.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
par(mfrow=c(2,1))
par(mar=c(0,5,4,4))
plot(CID,type="h",yaxs="i",lwd=3,col="coral1",ylim=c(0,110),xaxt="n",cex.axis=1,
     ylab="rel. Intensity",main=expression("Example MS2 spectra: PE 33:1 [M+NA]"^"+")
       )
text(90,100,labels="CID",cex=2,col="coral1")
par(mar=c(4,5,0,4))
plot(EAD,type="h",yaxs="i",lwd=3,col="palegreen3",ylim=c(110,0),cex.axis=1,ylab="rel. Intensity")
text(149,100,labels="EAD",cex=2,col="palegreen3")
dev.off()
