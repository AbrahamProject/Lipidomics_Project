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

################################################################################

# lipid_class distribution

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Summary plots/Lipid class distribution/"


Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
class_distribution(CID_EP,fragmentation = "CID",phase="EP")
dev.off()

Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
class_distribution(CID_SP,fragmentation = "CID",phase="SP")
dev.off()

Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
class_distribution(EAD_EP,fragmentation = "EAD",phase="EP")
dev.off()

Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
class_distribution(EAD_SP,fragmentation = "EAD",phase="SP")
dev.off()

################################################################################

# lipid_count number of identified in EAD/CID

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Summary plots/Lipid counts/"


Cairo(file=paste0(SAVE,"EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
lipid_count(CID$`CID-EP-data`,EAD$`EAD-EP-data`,phase = "EP")
dev.off()


Cairo(file=paste0(SAVE,"SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
lipid_count(CID$`CID-SP-data`,EAD$`EAD-SP-data`,phase = "SP")
dev.off()
################################################################################

# lipid_depth percentage number of identified in EAD/CID

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Summary plots/Lipid depth percentage/"


Cairo(file=paste0(SAVE,"EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =30 )
lipid_depth(CID$`CID-EP-data`,EAD$`EAD-EP-data`,phase = "EP")
dev.off()


Cairo(file=paste0(SAVE,"SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
lipid_depth(CID$`CID-SP-data`,EAD$`EAD-SP-data`,phase = "SP")
dev.off()
################################################################################

# lipid_isomerers percentage for sum composition isomers

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Summary plots/Lipid isomer distribution/"

Cairo(file= paste0(SAVE,"EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
lipid_isomers(CID$`CID-EP-data`,EAD$`EAD-EP-data`,phase = "EP")
dev.off()


Cairo(file=paste0(SAVE,"SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
lipid_isomers(CID$`CID-SP-data`,EAD$`EAD-SP-data`,phase = "SP")
dev.off()

################################################################################

# lipid_isomers_identified for highest level available

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Summary plots/Lipid isomer highest level identified/"


Cairo(file=paste0(SAVE,"EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =30 )
lipid_isomers_identified(CID$`CID-EP-data`,EAD$`EAD-EP-data`,phase="EP")
dev.off()


Cairo(file=paste0(SAVE,"SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
lipid_isomers_identified(CID$`CID-SP-data`,EAD$`EAD-SP-data`,phase="SP")
dev.off()

################################################################################

# csv file for resolved isomers by EAD

SAVE <-  "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Summary plots/Resolved isomers/"
Directory <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Summary plots/Resolved isomers"

Cairo(file=paste0(SAVE,"EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
resolved_isomers(CID$`CID-EP-data`,EAD$`EAD-EP-data`,Directory,phase = "EP")
dev.off()


Cairo(file=paste0(SAVE,"SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
resolved_isomers(CID$`CID-SP-data`,EAD$`EAD-SP-data`,Directory,phase = "SP")
dev.off()


DF<-data.table(CID=c(37,33),EAD=c(16,14))
rownames(DF)<-c("EP","SP")
png("Resolved_isomers_percentage.png")
p <- tableGrob(DF)
grid.arrange(p)
dev.off()

################################################################################

# venn() overlap of sum product id EAD,CID
Directory <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Summary plots/Overlap of sum composition"

venn(CID$`CID-EP-data`,EAD$`EAD-EP-data`,Directory)


################################################################################

################################################################################

# DG 16:1 18:1 Isomers
par(mfrow=c(1,1))
while (dev.cur()>1) dev.off()
p<-plot_DG_16_1_18_1()

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Summary plots/DG 16_1_18_1/"

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"CID.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 18,width =26 )
p[[1]]
while (dev.cur()>1) dev.off()


Cairo(file=paste0(SAVE,"EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 18,width =26 )
p[[2]]
while (dev.cur()>1) dev.off()



################################################################################

# Chromatogram

while (dev.cur()>1) dev.off()
p<-chrom(CID_EP,EAD_EP)

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Summary plots/chromatogram/"

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"CID.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 18,width =26 )
p[[1]]
while (dev.cur()>1) dev.off()


Cairo(file=paste0(SAVE,"EAD.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 18,width =26 )
p[[2]]
while (dev.cur()>1) dev.off()


################################################################################

# Growth Curve

while (dev.cur()>1) dev.off()

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Summary plots/Growth Curve/"

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"Growth_Curve.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 18,width =33 )
growth()
while (dev.cur()>1) dev.off()


################################################################################

# MS2 library match

while (dev.cur()>1) dev.off()

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Summary plots/MS2 library match/"

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"test.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )

ms2()
while (dev.cur()>1) dev.off()


