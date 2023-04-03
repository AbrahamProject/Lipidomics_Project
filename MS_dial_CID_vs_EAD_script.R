################################################################################
#
# Abraham Moyal
# 16.02.23
# Data Analysis
#
################################################################################
################################################################################
# Libraries

install.packages("devtools")
library(devtools)

################################################################################
# source of all used functions

source_url("https://raw.githubusercontent.com/AbrahamProject/Lipidomics_Project/main/Data_analysis_function.R" )
source_url("https://raw.githubusercontent.com/AbrahamProject/Lipidomics_Project/main/Processing_function.R" )

################################################################################
#-------------------------------------------------------------------------------
################################################################################
# data import: use Rscript_data_processing.R


CID_EP <- DATA$CID$`CID-EP-data`
CID_SP <- DATA$CID$`CID-SP-data`

EAD_EP <- DATA$EAD$`EAD-EP-data`
EAD_SP <- DATA$EAD$`EAD-SP-data`



# CID/EAD data location and name  for functions
#----------------------------------------------

CID_wd <- "W:/users/Abraham/Exp_004_version_2/MS dial Output/CID"
CID_name <-"Area_CID15min.txt"

EAD_wd <- "W:/users/Abraham/Exp_004_version_2/MS dial Output/EAD"
EAD_name<-"EAD_15min_all.txt"

################################################################################

# Feature Overlap



# Venn Diagram
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Feature Overlap Venn Diagram"

S<-feature_overlap(CID_wd,CID_name,EAD_wd,EAD_name,directory_venn_diagram = SAVE)



# Barplot number of Features 

while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Feature Barplot number of feature per method/"

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 18.6,width =20 )
S$Barplot
while (dev.cur()>1) dev.off()



while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Feature RT_MZ map/"

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 22,width =32 )
S$`MZ-RT plot`
while (dev.cur()>1) dev.off()






################################################################################

# Lipid Overlap



# Venn Diagram
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Lipid Overlap Venn Diagram"

S<-lipid_overlap(CID_wd,CID_name,EAD_wd,EAD_name,directory_venn_diagram = SAVE)


# Barplot number of Features 

while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Lipid Barplot number of lipid per method/"

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 18.6,width =20 )
S$Barplot
while (dev.cur()>1) dev.off()



################################################################################

# Lipid Overlap per class



# Absolute
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Overlap per Class Barplot/Absolute/"

S<-overlap_barplot(CID_wd,CID_name,EAD_wd,EAD_name,percentage = F)


while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =28 )
S[[1]]
while (dev.cur()>1) dev.off()


# Percentage
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Overlap per Class Barplot/Percentage/"

S<-overlap_barplot(CID_wd,CID_name,EAD_wd,EAD_name,percentage = T)


while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =28 )
S[[1]]
while (dev.cur()>1) dev.off()









################################################################################

#   Coefficient of Variation


# CV all
#-------

S<-cv_all(CID_wd,CID_name,EAD_wd,EAD_name)

# Barplot:

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/CV all/Barplot/"

while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 28,width =20 )
S[[1]]
while (dev.cur()>1) dev.off()


# Pie Chart
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/CV all/Pie Chart/"

while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 28,width =20 )
S[[2]]
while (dev.cur()>1) dev.off()







# CV Standard
#------------

P<-cv_standard(CID_wd,CID_name,EAD_wd,EAD_name)
# Barplot:

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/CV standard/Barplot/"

while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 28,width =20 )
P[[1]]
while (dev.cur()>1) dev.off()


# Pie Chart
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/CV standard/Pie Chart/"

while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 28,width =20 )
P[[2]]
while (dev.cur()>1) dev.off()




# CV blanks
#------------

S<-cv_blanks(CID_wd,CID_name,EAD_wd,EAD_name)
# Barplot:

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/CV blanks/Barplot/"

while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 28,width =20 )
S[[1]]
while (dev.cur()>1) dev.off()


# Pie Chart
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/CV blanks/Pie Chart/"

while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 28,width =20 )
S[[2]]
while (dev.cur()>1) dev.off()



################################################################################

#   CV vs Intensity

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/CV vs Intensity/"


cv_vs_intensity(CID_wd,CID_name,EAD_wd,EAD_name)

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =28 )
cv_vs_intensity(CID_wd,CID_name,EAD_wd,EAD_name)
while (dev.cur()>1) dev.off()

################################################################################

#   log EAD vs log CID

# EP
#---


   # Venn Diagramm
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Venn Diagram Enriched Lipids/EP"
S<-CID_vs_EAD(CID_EP,EAD_EP,phase = "EP",fragmentation = "CID" ,Directory = SAVE)

# log vs log
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/log vs log/EP/"

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
S$`log vs. log`
while (dev.cur()>1) dev.off()

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"CID_EAD_paper.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
S$`log vs. log`
while (dev.cur()>1) dev.off()

# Pie Chart Enriched
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Pie Chart Enriched Lipids/EP/"

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 25,width =35 )
S$`Piechart: IDs per class all`
while (dev.cur()>1) dev.off()


# Table
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Table IDs Enriched Lipids/EP/"

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29,width = 40 )
S$`Venn Diagram unique IDs`
while (dev.cur()>1) dev.off()




# SP
#---


# Venn Diagramm
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Venn Diagram Enriched Lipids/SP"
S<-CID_vs_EAD(CID_SP,EAD_SP,phase = "SP",fragmentation = "CID" ,Directory = SAVE)

# log vs log
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/log vs log/SP/"

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
S$`log vs. log`
while (dev.cur()>1) dev.off()


# Pie Chart Enriched
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Pie Chart Enriched Lipids/SP/"

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 25,width =35 )
S$`Piechart: IDs per class all`
while (dev.cur()>1) dev.off()


# Table
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Table IDs Enriched Lipids/SP/"

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29,width = 40 )
S$`Venn Diagram unique IDs`
while (dev.cur()>1) dev.off()




################################################################################

#   signal to noise EAD vs CID

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Signal to Noise/"

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width = 30 )
SN(CID_SP,EAD_SP)
while (dev.cur()>1) dev.off()




################################################################################

#   Gaussian similarity EAD vs CID

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Gaussian Similarity/"

while (dev.cur()>1) dev.off()

S<-gaussian(CID_EP,EAD_EP)

while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_T1.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 25,width = 25 )
S$`CID TEIL I`
while (dev.cur()>1) dev.off()


while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_T2.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 25,width = 25 )
S$`CID TEIL II`
while (dev.cur()>1) dev.off()



while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"EAD_T1.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 25,width = 25 )
S$`EAD TEIL I`
while (dev.cur()>1) dev.off()



while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"EAD_T2.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 25,width = 25 )
S$`EAD TEIL II`
while (dev.cur()>1) dev.off()





#################################################################################################
# Signal to Noise ratio
setwd("C:/FUNCTION_R_DATA")


CID_S <- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_CID.mzML",mode = "onDisk") 

EAD_S<- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_EAD.mzML",mode = "onDisk")


C<-read_delim("W:/users/Abraham/Exp 004/Wild Type 3/CID.txt", delim = "\t", escape_double = FALSE,trim_ws = TRUE)


E<-read_delim("W:/users/Abraham/Exp 004/Wild Type 3/EAD.txt", delim = "\t", escape_double = FALSE,trim_ws = TRUE)





data<-C
data_1 <- data[which(data$Ontology !="Unknown"),]
data_2 <- data_1[which(data_1$Ontology !="null"),]
data_3 <- data_2[which(data_2$Ontology !="Others"),]

C_lipid <- data_3


data<-E
data_1 <- data[which(data$Ontology !="Unknown"),]
data_2 <- data_1[which(data_1$Ontology !="null"),]
data_3 <- data_2[which(data_2$Ontology !="Others"),]

E_lipid <- data_3



data_C<-matrix(ncol=2,nrow=nrow(C))
data_C[,1] <- c(1:nrow(C))

for(i in 1:nrow(data_C)){
  sn<-C$`S/N`[i]
  
  data_C[i,2]<-sn
  
}



data_E<-matrix(ncol=2,nrow=nrow(E))
data_E[,1] <- c(1:nrow(E))

for(i in 1:nrow(data_E)){
  sn<-E$`S/N`[i]
  
  data_E[i,2]<-sn
  
}


data_C_lipid<-matrix(ncol=2,nrow=nrow(C_lipid))
data_C_lipid[,1] <- c(1:nrow(C_lipid))

for(i in 1:nrow(data_C_lipid)){
  sn<-C_lipid$`S/N`[i]
  
  data_C_lipid[i,2]<-sn
  
}



data_E_lipid<-matrix(ncol=2,nrow=nrow(E_lipid))
data_E_lipid[,1] <- c(1:nrow(E_lipid))

for(i in 1:nrow(data_E_lipid)){
  sn<-E_lipid$`S/N`[i]
  
  data_E_lipid[i,2]<-sn
  
}
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Signal to Noise/"

dev.off()

Cairo(file=paste0(SAVE,"SN.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 14,width =30 )

layout(matrix(c(1,2),ncol=2))
par(mar=c(3,5,5,2),cex.lab=1.5,cex.axis=1.4,cex.main=1.5)

boxplot(data_C[,2],data_E[,2],log="y",col=c("coral1","palegreen3"),names=c("CID","EAD"),
        yaxt="n",ylab="Signal to Noise Ratio (S/N)",
        main="S/N: All Features")
axis(2,cex=1.2)
boxplot(data_C_lipid[,2],data_E_lipid[,2],log="y",col=c("coral1","palegreen3"),names=c("CID","EAD"),
        yaxt="n",ylab="Signal to Noise Ratio (S/N)",
        main="S/N: All Lipids")
axis(2,cex=1.2)
dev.off()
mean(data_C[,2])
mean(data_E[,2])

mean(data_C_lipid[,2])
mean(data_E_lipid[,2])
#####################################################################################

# Number of data point per peak


setwd("C:/FUNCTION_R_DATA")


CID_S <- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_CID.mzML",mode = "onDisk") 

EAD_S<- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_EAD.mzML",mode = "onDisk")


C<-read_delim("W:/users/Abraham/Exp 004/Wild Type 3/CID.txt", delim = "\t", escape_double = FALSE,trim_ws = TRUE)


E<-read_delim("W:/users/Abraham/Exp 004/Wild Type 3/EAD.txt", delim = "\t", escape_double = FALSE,trim_ws = TRUE)




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
  rm(chr_raw)
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
  rm(chr_raw)
  
}



data<-C
data_1 <- data[which(data$Ontology !="Unknown"),]
data_2 <- data_1[which(data_1$Ontology !="null"),]
data_3 <- data_2[which(data_2$Ontology !="Others"),]

C_lipid <- data_3


data<-E
data_1 <- data[which(data$Ontology !="Unknown"),]
data_2 <- data_1[which(data_1$Ontology !="null"),]
data_3 <- data_2[which(data_2$Ontology !="Others"),]

E_lipid <- data_3




data_C_lipid<-matrix(ncol=2,nrow=nrow(C_lipid))
data_C_lipid[,1] <- c(1:nrow(C_lipid))
i=36
for(i in 1:nrow(data_C_lipid)){
  rt1<-C_lipid$`RT left(min)`[i]*60
  rt2<-C_lipid$`RT right (min)`[i]*60
  rtr<-c(rt1,rt2)
  mz<-C_lipid$`Precursor m/z`[i]
  mzr <-c(mz-0.05,mz+0.05)
  
  chr_raw <- chromatogram(CID_S, mz = mzr, rt = rtr)
  
  plot(chr_raw)
  int<-intensity(chr_raw[1,1])
  
  int<- na.omit(int)
  
  l<-as.numeric(length(int))
  
  data_C_lipid[i,2]<-l
  rm(chr_raw)
}



data_E_lipid<-matrix(ncol=2,nrow=nrow(E_lipid))
data_E_lipid[,1] <- c(1:nrow(E_lipid))

for(i in 1:nrow(data_E_lipid)){
  rt1<-E_lipid$`RT left(min)`[i]*60
  rt2<-E_lipid$`RT right (min)`[i]*60
  rtr<-c(rt1,rt2)
  mz<-E_lipid$`Precursor m/z`[i]
  mzr <-c(mz-0.05,mz+0.05)
  
  chr_raw <- chromatogram(EAD_S, mz = mzr, rt = rtr)
  
  
  int<-intensity(chr_raw[1,1])
  
  int<- na.omit(int)
  
  l<-as.numeric(length(int))
  
  data_E_lipid[i,2] <-l
  rm(chr_raw)
  
}

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Number of data points/"

dev.off()

Cairo(file=paste0(SAVE,"data_points.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 15,width =30 )



layout(matrix(c(1,2),ncol=2))
par(mar=c(3,5,5,2),cex.lab=1.5,cex.axis=1.4,cex.main=1.5)

boxplot(data_C[ind_C,2],data_E[ind_E,2],col=c("coral1","palegreen3"),names=c("CID","EAD"),
        yaxt="n",ylab="# data points per peak",
        main="# Data Points Per Peak: All Features")
axis(2,cex=1.2)
boxplot(data_C_lipid[,2],data_E_lipid[,2],col=c("coral1","palegreen3"),names=c("CID","EAD"),
        yaxt="n",ylab="# data points per peak",
        main="# Data Points Per Peak: All Lipids")
axis(2,cex=1.2)


dev.off()
dev.off()





#################################################################################################
  # Gaussian similarity
setwd("C:/FUNCTION_R_DATA")


CID_S <- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_CID.mzML",mode = "onDisk") 

EAD_S<- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_EAD.mzML",mode = "onDisk")

setwd("W:/users/Abraham/Exp_004_version_2/Tables")
d<-read.csv("gaussian.csv",header = F)


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Gaussian Similarity/"

dev.off()

Cairo(file=paste0(SAVE,"GASI.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 14,width =30 )

layout(matrix(c(1,2),ncol=2))
par(mar=c(3,5,5,2),cex.lab=1.5,cex.axis=1.4,cex.main=1.5)

boxplot(d[,1],d[,2],col=c("coral1","palegreen3"),names=c("CID","EAD"),
        yaxt="n",ylab="Dot Product Score",
        main="Gaussian Similarity: All Features",ylim=c(0.4,1))
axis(2,cex=1.2)
boxplot(d[,3],d[,4],col=c("coral1","palegreen3"),names=c("CID","EAD"),
        yaxt="n",ylab="Dot Product Score",
        main="Gaussian Similarity: All Lipids",ylim=c(0.4,1))
axis(2,cex=1.2)
dev.off()
mean(d[,1],na.rm=T)
mean(d[,2],na.rm=T)

mean(d[,3],na.rm=T)
mean(d[,4],na.rm=T)

var(d[,1],na.rm=T)
var(d[,2],na.rm=T)

var(d[,3],na.rm=T)
var(d[,4],na.rm=T)


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Barplot overlap/"

dev.off()

Cairo(file=paste0(SAVE,"Barplot.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 12,width =28 )

Overlap_enriched_per_ko(EAD,CID,phase="EP",fragmentation="CID")
dev.off()



#######################################################################################################
# score distribution

setwd("W:/users/Abraham/Exp_004_version_2/Tables")

score<-read.csv("scores.csv")

par(cex.main=1.5,cex.lab=1.6,cex.axis=1.4,mar=c(5,5,5,1))


c_min<-min(unlist(score[,1:3]))-0.05
c_max<-max(unlist(score[,1:3]))+0.05
layout(matrix(c(1,2),ncol=2))
plot(1:10,score[,1],type="l",col="coral1",lwd=2,ylab="MS-dial score",xlab="Candidate Rank",main="",ylim=c(c_min,c_max),xaxt="n")
axis(1,at=1:10,labels = 1:10)
points(1:10,score[,1],pch=16,col="coral1",cex=2)
points(1:10,score[,2],type="l",col="coral1",lwd=2)
points(1:10,score[,2],pch=16,col="coral1")
points(1:10,score[,3],type="l",col="coral1",lwd=2)
points(1:10,score[,3],pch=16,col="coral1",cex=2)
mtext(side=3,"CID data",cex=1.4)

c_min<-min(unlist(score[,4:6]))-0.05
c_max<-max(unlist(score[,4:6]))+0.05

plot(1:10,score[,4],type="l",col="palegreen3",lwd=2,ylab="MS-dial score",xlab="Candidate Rank",main="",ylim=c(c_min,c_max),xaxt="n")
axis(1,at=1:10,labels = 1:10)

points(1:10,score[,4],pch=16,col="palegreen3",cex=2)
points(1:10,score[,5],type="l",col="palegreen3",lwd=2)
points(1:10,score[,5],pch=16,col="palegreen3",cex=2)
points(1:10,score[,6],type="l",col="palegreen3",lwd=2)
points(1:10,score[,6],pch=16,col="palegreen3",cex=2)
mtext(side=3,"EAD data",cex=1.4)
mtext(side=3,"Matching Score Distribution (Top 10 Candidates)",outer = T,line=-2,cex=1.8)

p<-recordPlot()
plot.new()



SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/Score distribution/"

dev.off()

Cairo(file=paste0(SAVE,"Barplot.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 15,width =28 )

p
dev.off()

