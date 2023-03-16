################################################################################
#
# Abraham Moyal
# 16.02.23
# Data Analysis
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
      type="png",bg="white",dpi=300,units = "cm", height = 15,width =21 )
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
      type="png",bg="white",dpi=300,units = "cm", height = 15,width =21 )
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
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =28 )
S[[1]]
while (dev.cur()>1) dev.off()


# Pie Chart
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/CV all/Pie Chart/"

while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =28 )
S[[2]]
while (dev.cur()>1) dev.off()







# CV Standard
#------------

S<-cv_standard(CID_wd,CID_name,EAD_wd,EAD_name)
# Barplot:

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/CV standard/Barplot/"

while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =28 )
S[[1]]
while (dev.cur()>1) dev.off()


# Pie Chart
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/CV standard/Pie Chart/"

while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =28 )
S[[2]]
while (dev.cur()>1) dev.off()






# CV blanks
#------------

S<-cv_blanks(CID_wd,CID_name,EAD_wd,EAD_name)
# Barplot:

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/CV blanks/Barplot/"

while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =28 )
S[[1]]
while (dev.cur()>1) dev.off()


# Pie Chart
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/CID vs EAD/CV blanks/Pie Chart/"

while (dev.cur()>1) dev.off()

Cairo(file=paste0(SAVE,"CID_EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =28 )
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

