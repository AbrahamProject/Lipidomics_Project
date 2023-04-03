################################################################################
#
# Abraham Moyal
# 16.02.23
# Data Processing
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

# CID data processing

# Input for all data processing


Directory <- "W:/users/Abraham/Exp 004/MS data processing/MS dial output/CID Alignment"
file <-"Area_CID15min.txt"
FC<-10
n_CID <- c(1,2,5,6,8,10)
fix <- 4
################################################################################

# Data import 

dataframe <- data_import(Directory,file)

data <-dataframe$data

standard <- dataframe$standard[n_CID,]

Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Raw data/CID_EP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(data)
dev.off()

Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Raw data/CID_SP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(data,phase = "SP")
dev.off()
################################################################################

# Filtering for ms2 identified lipids

data <- filtering(data)


Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Filtered/CID_EP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(data)
dev.off()

Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Filtered/CID_SP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(data,phase = "SP")
dev.off()
################################################################################

# Remove features based on blank samples

type="EP"

dataframe_EP <- blanking(data,type,FC=10)
Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Blanked/CID_EP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_EP)
dev.off()

type <- "SP"

dataframe_SP <- blanking(data,type,FC)
Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Blanked/CID_SP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_SP,phase="SP")
dev.off()
################################################################################

# Normalizaion using internal standard

type="EP"

dataframe_EP <- normalization(dataframe_EP,type,standard,fix)

Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Normalized/CID_EP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_EP,type = "Concentration")
dev.off()

type="SP"

dataframe_SP <- normalization(dataframe_SP,type,standard,fix)

Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Normalized/CID_SP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_SP,phase="SP",type = "Concentration")
dev.off()
################################################################################

#   OD correction using OD data
################################################################################

dataframe_EP <- OD_correction(dataframe_EP)
Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/OD corrected/CID_EP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_EP,type = "Concentration")
dev.off()


dataframe_SP <- OD_correction(dataframe_SP)

Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/OD corrected/CID_SP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_SP,phase="SP",type = "Concentration")
dev.off()
################################################################################

# CID data storage

CID <- list(dataframe_EP,
            dataframe_SP
)

names(CID) <-c("CID-EP-data","CID-SP-data")
################################################################################
#-------------------------------------------------------------------------------
################################################################################

# EAD data processing
#--------------------
# Input for all data processing


Directory <- "W:/users/Abraham/Exp 004/MS data processing/MS dial output/EAD  Alignment"
file <-"EAD_15min_all.txt"
FC<-10
n_EAD <- c(1,3)
fix <- 1
################################################################################

# Data import 

dataframe <- data_import(Directory,file)

data <-dataframe$data

standard <- dataframe$standard[n_EAD,]


Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Raw data/EAD_EP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(data,fragmentation = "EAD")
dev.off()

Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Raw data/EAD_SP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(data,phase = "SP",fragmentation = "EAD")
dev.off()
################################################################################

# Filtering for ms2 identified lipids

data <- filtering(data)


Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Filtered/EAD_EP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(data,fragmentation = "EAD")
dev.off()

Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Filtered/EAD_SP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(data,phase = "SP",fragmentation = "EAD")
dev.off()
################################################################################

# Remove features based on blank samples

type="EP"

dataframe_EP <- blanking(data,type,FC=10)
Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Blanked/EAD_EP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_EP,fragmentation = "EAD")
dev.off()

type <- "SP"

dataframe_SP <- blanking(data,type,FC)
Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Blanked/EAD_SP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_SP,phase="SP",fragmentation = "EAD")
dev.off()
################################################################################

# Normalizaion using internal standard

type="EP"

dataframe_EP <- normalization(dataframe_EP,type,standard,fix)

Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Normalized/EAD_EP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_EP,fragmentation = "EAD",type = "Concentration")
dev.off()

type="SP"

dataframe_SP <- normalization(dataframe_SP,type,standard,fix)

Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/Normalized/EAD_SP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_SP,phase="SP",fragmentation = "EAD",type = "Concentration")
dev.off()
################################################################################

#   OD correction using OD data
################################################################################

dataframe_EP <- OD_correction(dataframe_EP)
Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/OD corrected/EAD_EP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_EP,fragmentation = "EAD",type = "Concentration")
dev.off()


dataframe_SP <- OD_correction(dataframe_SP)

Cairo(file="//pasteur/SysBC-Home/amoyal/Desktop/Report/6) Plots/1) Data Pre processing/OD corrected/EAD_SP.png",
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_SP,phase="SP",fragmentation = "EAD",type = "Concentration")
dev.off()
################################################################################
# EAD data storage

EAD <- list(dataframe_EP,
            dataframe_SP
)

names(EAD) <-c("EAD-EP-data","EAD-SP-data")


################################################################################

# Data storage


DATA <- list(CID=CID,EAD=EAD)

################################################################################
