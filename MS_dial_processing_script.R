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


Directory <- "W:/users/Abraham/Exp_004_version_2/MS dial Output/CID"
file <-"Area_CID15min.txt"
FC<-10
n_CID <- c(1,2,5,6,8,10)
fix <- 4
################################################################################

# Data import 

dataframe <- data_import(Directory,file)

data <-dataframe$data

standard <- dataframe$standard[n_CID,]

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Pre-processing/Intensity distribution/"

Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(data,fragmentation = "CID",phase = "EP")
dev.off()

Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(data,phase = "SP",fragmentation = "CID")
dev.off()
################################################################################

# Filtering for ms2 identified lipids

data <- filtering(data)



distribution(data)



distribution(data,phase = "SP")

################################################################################

# Remove features based on blank samples

type="EP"

dataframe_EP <- blanking(data,type,FC=10)

distribution(dataframe_EP)

type <- "SP"

dataframe_SP <- blanking(data,type,FC=10)

distribution(dataframe_SP,phase="SP")

################################################################################

# Normalizaion using internal standard

type="EP"

dataframe_EP <- normalization(dataframe_EP,type,standard,fix)


distribution(dataframe_EP,type = "Concentration")


type="SP"

dataframe_SP <- normalization(dataframe_SP,type,standard,fix)


distribution(dataframe_SP,phase="SP",type = "Concentration")

################################################################################

#   OD correction using OD data
################################################################################

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Pre-processing/Concentration distribution/"

dataframe_EP <- OD_correction(dataframe_EP)
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_EP,type = "Concentration",phase = "EP",fragmentation = "CID")
dev.off()


dataframe_SP <- OD_correction(dataframe_SP)

Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_SP,phase="SP",type = "Concentration",fragmentation = "CID")
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

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Pre-processing/Intensity distribution/"


Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(data,fragmentation = "EAD",phase = "EP")
dev.off()

Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(data,phase = "SP",fragmentation = "EAD")
dev.off()
################################################################################

# Filtering for ms2 identified lipids

data <- filtering(data)



distribution(data,fragmentation = "EAD")



distribution(data,phase = "SP",fragmentation = "EAD")
################################################################################

# Remove features based on blank samples

type="EP"

dataframe_EP <- blanking(data,type,FC=10)

distribution(dataframe_EP,fragmentation = "EAD")

type <- "SP"

dataframe_SP <- blanking(data,type,FC)

distribution(dataframe_SP,phase="SP",fragmentation = "EAD")
################################################################################

# Normalizaion using internal standard

type="EP"

dataframe_EP <- normalization(dataframe_EP,type,standard,fix)


distribution(dataframe_EP,fragmentation = "EAD",type = "Concentration")

type="SP"

dataframe_SP <- normalization(dataframe_SP,type,standard,fix)


distribution(dataframe_SP,phase="SP",fragmentation = "EAD",type = "Concentration")
################################################################################

#   OD correction using OD data
################################################################################

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Pre-processing/Concentration distribution/"

dataframe_EP <- OD_correction(dataframe_EP)
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
distribution(dataframe_EP,fragmentation = "EAD",type = "Concentration",phase = "EP")
dev.off()


dataframe_SP <- OD_correction(dataframe_SP)

Cairo(file=paste0(SAVE,"EAD_SP.png"),
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
