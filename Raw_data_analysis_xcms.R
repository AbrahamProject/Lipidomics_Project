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
      type="png",bg="white",dpi=300,units = "cm", height =18,width =35 )
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





ms2_cid_ead<-function(type=T){
  library("grid")
  library("cowplot")
  setwd("W:/users/Abraham/Exp_004_version_2/Data/msp file standards/CID vs EAD")
  files <- list.files(pattern = "msp")
  
  CID <- files[grep("CID",files)]
  EAD <- files[grep("EAD",files)]
  
  C<-vector(mode="list",5) 
  names(C)<-gsub("\\..*","",CID)
  
  for(i in 1:length(CID)){
    
    library <- read.delim(CID[i], header=FALSE, comment.char="#",sep="\t", dec=",", quote="")
    
    
    library1 <- library[11:nrow(library),1]
    
    library2 <- matrix(library1,ncol=2,byrow=T)
    
    
    C[[i]] <- library2 
    
    
  }
  
  
  
  E<-vector(mode="list",5) 
  names(E)<-gsub("\\..*","",EAD)
  
  for(i in 1:length(EAD)){
    
    library <- read.delim(EAD[i], header=FALSE, comment.char="#",sep="\t", dec=",", quote="")
    
    
    library1 <- library[11:nrow(library),1]
    
    library2 <- matrix(library1,ncol=2,byrow=T)
    
    
    E[[i]] <- library2 
    
    
  }
  
  
  
  
  
  
  title<-c(
    "PE 36:2",
    "PE 33:2",
    "PE 34:2")
  
  
  
  test<-matrix(ncol=3,nrow=3)
  
  test[1,]<-1
  test[2,]<-c(2,3,4)
  test[3,]<-c(5,6,7)
  
  layout(test,heights = c(1,8,8))
  par(mar=c(0,2,0,2))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')  
  mtext(3,text="MS2 Spectra: CID vs EAD fragmentation",line=-2.5,cex=1.5)
  
  
  
  
  for(i in 1:3){
    d2 <- C[[i]]
    d4 <- E[[i]]
    
    
    x2 <- as.numeric(d2[,1])
    y2 <- as.numeric(d2[,2])/max(as.numeric(d2[,2]))*100
    
    
    x4 <- as.numeric(d4[,1])
    y4 <- as.numeric(d4[,2])/max(as.numeric(d4[,2]))*100
    
    
    xmin=min(as.numeric(c(x2,x4)))
    
    xmax=max(as.numeric(c(x2,x4)))
    
    y2max=max(y2)
    
    y4max=max(y4)
    
    par(mar=c(0,2,5,2),cex.lab=1,cex.axis=1,cex.main=1)
    
    plot(x2,y2,type= "h",xaxt="n",xlab="",ylab="",main=paste0(title[i]," "  ),col="coral1",lwd=2,xaxs = "i",
         yaxs = "i",xlim=c(xmin,xmax+10),ylim=c(0,150),yaxt="n",cex.main=1.8)
    mtext(3,text="CID ms2 Spectra",adj=0.1,col="coral1",line=-2,cex=1.2)
    
    
    
    
    
  }
  
  
  
  
  for(i in 1:3){
    d2 <- C[[i]]
    d4 <- E[[i]]
    
    
    x2 <- as.numeric(d2[,1])
    y2 <- as.numeric(d2[,2])/max(as.numeric(d2[,2]))*100
    
    
    x4 <- as.numeric(d4[,1])
    y4 <- as.numeric(d4[,2])/max(as.numeric(d4[,2]))*100
    
    
    xmin=min(as.numeric(c(x2,x4)))
    
    xmax=max(as.numeric(c(x2,x4)))
    
    y2max=max(y2)
    
    y4max=max(y4)
    
    
    
    par(mar=c(5,2,0,2))
    plot(x4,y4,type= "h",xlab="m/z",ylab="",col="palegreen3",lwd=2,xaxs = "i",
         yaxs = "i",xlim=c(xmin,xmax+10) ,ylim=c(150,0),yaxt="n",cex.axis=1.4,cex.lab=1.4)  
    mtext(1,text="EAD ms2 Spectra",adj=0.1,col="palegreen3",line=-1.2,cex=1.2)
    
    
    
    
  }
  p<- recordPlot()
  par(mar=c(0,4,5,4),mfrow=c(1,1))
  plot.new()
  
  return(p)
  
  
  
  
}

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/1) Raw data plots/"


while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_CID_MS2.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
ms2_cid_ead()
while (dev.cur()>1) dev.off()


############################################################################################

setwd("C:/FUNCTION_R_DATA")

layout(matrix(c(1,2,3,4),ncol=2,byrow=T))

CID_S <- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_CID.mzML",mode = "onDisk") 

EAD_S<- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_EAD.mzML",mode = "onDisk")

ch_CID <- chromatogram(CID_S, aggregationFun = "max")

ch_EAD<- chromatogram(EAD_S, aggregationFun = "max")




ms1 <- CID_S %>%
  filterMsLevel(1L) 

rt1<-c(as.numeric((rtime(ms1))),100000)


ms2 <- CID_S %>%
  filterMsLevel(2L) 

rt2<-as.numeric((rtime(ms2)))

C<-matrix(nrow=length(rtime(ms1)),ncol=2)
C[,1]<-1:length(rtime(ms1))
colnames(C)<-c("MS1","MS2")
i=1
for(i in 1:(length(rt1)-1)){
  index<-which(rt2>rt1[i] & rt2<rt1[i+1])
  
  C[i,2]<-length(index)
  
  
  
}


C_percentage<- 100-( length(which(C[,2]==0)))/nrow(C)*100
C_dist<-C[which(C[,2]!=0),2]
CC<-nrow(C)
par(mar=c(5,5,5,5),cex.lab=1.5,cex.axis=1.2,cex.main=1.5)
plot(rt1[C[,1]],C[,2],type="h",col="lightpink",ylim=c(0,30),las=1,xlab="mz",ylab="# of ms2 scans per Cycle",main="CID: Number of MS2 scans per Cycle ")

points(rtime(ch_CID[1,1]),intensity(ch_CID[1,1])/max(intensity(ch_CID[1,1]))*30,type="l",col="red",lwd=2)
axis(4,at=c(7.5,15,22.5,30),labels = c(25,50,75,100),las=2)
mtext(4,text="Rel. Intensity",line=3,cex=1.5)
legend("topright",legend = c("Base Peak Chromatogram","# MS2 scans"), lwd=3, col=c("red","lightpink"),bty="n" )

ms1 <- EAD_S %>%
  filterMsLevel(1L) 

rt1<-c(as.numeric((rtime(ms1))),100000)



ms2 <- EAD_S %>%
  filterMsLevel(2L) 

rt2<-as.numeric((rtime(ms2)))

C<-matrix(nrow=length(rtime(ms1)),ncol=2)
C[,1]<-1:length(rtime(ms1))
colnames(C)<-c("MS1","MS2")
i=1
for(i in 1:(length(rt1)-1)){
  index<-which(rt2>rt1[i] & rt2<rt1[i+1])
  
  C[i,2]<-length(index)
  
  
  
}
E_percentage<- 100-( length(which(C[,2]==0)))/nrow(C)*100
E_dist<-C[which(C[,2]!=0),2]
EE<-nrow(C)


par(mar=c(5,5,5,5),cex.lab=1.5,cex.axis=1.2,cex.main=1.5)
plot(rt1[C[,1]],C[,2],type="h",col="#93F1AC",ylim=c(0,30),las=1,xlab="mz",ylab="# of ms2 scans per Cycle",main="EAD: Number of MS2 scans per Cycle ")

points(rtime(ch_EAD[1,1]),intensity(ch_EAD[1,1])/max(intensity(ch_CID[1,1]))*30,type="l",col="#3D954C",lwd=2)
axis(4,at=c(7.5,15,22.5,30),labels = c(25,50,75,100),las=2)
mtext(4,text="Rel. Intensity",line=3,cex=1.5)
legend("topright",legend = c("Base Peak Chromatogram","# MS2 scans"), lwd=3, col=c("#3D954C","#93F1AC"),bty="n" )


a<-barplot(c(C_percentage,E_percentage),col=c("coral1","palegreen3"),names.arg = c("CID","EAD"),ylab="Percentage [%]",
           main="Percentage of MS1 scans with MS2 Spectra",cex.names = 2)
text(a,c(C_percentage,E_percentage),c(paste0(length(C_dist),"/",(CC)),paste0(length(E_dist),"/",(EE))),cex=2,pos=1)

boxplot(C_dist,E_dist,col=c("coral1","palegreen3"),names=c("CID","EAD"),cex.axis=2,yaxt="n",
        ylab="# of ms2 scans per Cycle",main="Distribution: Number of MS2 scans per Cycle ")
axis(2,cex.axis=1.2)

p<-recordPlot()
plot.new()

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/1) Raw data plots/"


while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"ms2_scans_per_cycle.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 30,width =37 )
p
while (dev.cur()>1) dev.off()


#####################################################################################

# Number of data point per peak


setwd("C:/FUNCTION_R_DATA")


CID_S <- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_CID.mzML",mode = "onDisk") 

EAD_S<- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_EAD.mzML",mode = "onDisk")


C<-read_delim("W:/users/Abraham/Exp 004/Wild Type 3/CID.txt", delim = "\t", escape_double = FALSE,trim_ws = TRUE)


E<-read_delim("W:/users/Abraham/Exp 004/Wild Type 3/EAD.txt", delim = "\t", escape_double = FALSE,trim_ws = TRUE)

data<-C
data_1 <- data[which(data$Ontology !="Unknown"),]
data_2 <- data_1[which(data_1$Ontology !="null"),]
data_3 <- data_2[which(data_2$Ontology !="Others"),]

C_lipid <- 
  
E_lipid <-



data_C<-matrix(ncol=2,nrow=nrow(C))
data_C[,1] <- c(1:nrow(C))
for(i in 1:nrow(data_C)){
  rt1<-C$`RT left(min)`[i]*60
  rt2<-C$`RT right (min)`[i]*60
  rtr<-c(rt1,rt2)
  mz<-C$`Precursor m/z`[i]
  mzr <-c(mz-0.05,mz+0.05)
  
  chr_raw <- chromatogram(CID_S, mz = mzr, rt = rtr)
  
  
  int<-intensity(chr_raw[1,1])
  
  int<- na.omit(int)
  
  l<-as.numeric(length(int))
  
  data_C[i,2]<-l
  
}



data_E<-matrix(ncol=2,nrow=nrow(E))
data_E[,1] <- c(1:nrow(E))

for(i in 1:nrow(data_E)){
  rt1<-E$`RT left(min)`[i]*60
  rt2<-E$`RT right (min)`[i]*60
  rtr<-c(rt1,rt2)
  mz<-E$`Precursor m/z`[i]
  mzr <-c(mz-0.05,mz+0.05)
  
  chr_raw <- chromatogram(EAD_S, mz = mzr, rt = rtr)
  
  
  int<-intensity(chr_raw[1,1])
  
  int<- na.omit(int)
  
  l<-as.numeric(length(int))
  
  data_E[i,2] <-l
  
}








