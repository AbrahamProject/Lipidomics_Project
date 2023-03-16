################################################################################
#
# Abraham Moyal
# 16.02.23
# Data Analysis MS1 Feature level  of MS-dial output
# 
################################################################################
# Libraries

install.packages("devtools")
library(devtools)

################################################################################
# source of all used functions

source_url( "https://raw.githubusercontent.com/AbrahamProject/data-processing/main/data_processing.R" )


################################################################################
#-------------------------------------------------------------------------------
################################################################################

# Load CID, EAD MS-dial output
#------------------------------


Directory <- "W:/users/Abraham/Exp_004_version_2/MS dial Output/CID"
file <-"Area_CID15min.txt"

dataframe <- data_import(Directory,file)
data_CID <-dataframe$data

Directory <- "W:/users/Abraham/Exp_004_version_2/MS dial Output/EAD"
file <-"EAD_15min_all.txt"

dataframe <- data_import(Directory,file)
data_EAD <-dataframe$data


par(cex.main=1.5,cex.lab=1.6,cex.axis=1.4,mar=c(5,6,5,5))

SAVE<-"W:/users/Abraham/Exp_004_version_2/Plots/2) MS-dial Features/"

Cairo(file=paste0(SAVE,"Barplot_Features_CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
par(cex.main=1.5,cex.lab=1.6,cex.axis=1.4,mar=c(5,6,5,5))
barplot(c(nrow(data_CID),nrow(data_EAD)),col=c("coral1","palegreen3"),names=c("CID","EAD"),cex.names = 2,
        main="Number of Total MS1 Features per Fragmentation Method",ylab="# Features")
dev.off()

Cairo(file=paste0(SAVE,"Feature_overlap_map.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
par(cex.main=1.5,cex.lab=1.6,cex.axis=1.4,mar=c(5,6,5,5))
plot(data_CID$`Average Rt(min)`,data_CID$`Average Mz`,col="coral1",xlab="RT",ylab="mz",pch=16,
     main="EAD/CID Feature Map ")
points(data_EAD$`Average Rt(min)`,data_EAD$`Average Mz`,col="palegreen3",xlab="RT",ylab="mz",pch=16)
legend("topright",pch=16,pt.cex =3,legend = c("CID","EAD"),col=c("coral1","palegreen3"),cex=1.8,bty="n")
dev.off()




Directory<-"W:/users/Abraham/Exp_004_version_2/Plots/2) MS-dial Features"

overlap(data_CID,data_EAD,Directory,mz_tol = 0.008,rt_tol = 0.1)




