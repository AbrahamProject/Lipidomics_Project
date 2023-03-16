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

################################################################################
# PCA analysis


# CID EP

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/PCA plot/"

P<-PCA(CID_EP,phase="EP",fragmentation = "CID")


Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$PCA
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Variance/"

Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$Variance
dev.off()



SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Loading 1/"

Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`PC1 loadings`
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Loading 2/"

Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`PC2 loadings`
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Loadings/"

Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$loadings
dev.off()





# CID SP

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/PCA plot/"

P<-PCA(CID_SP,phase="SP",fragmentation = "CID")


Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$PCA
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Variance/"

Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$Variance
dev.off()



SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Loading 1/"

Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`PC1 loadings`
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Loading 2/"

Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`PC2 loadings`
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Loadings/"

Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$loadings
dev.off()




# EAD EP

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/PCA plot/"

P<-PCA(EAD_EP,phase="EP",fragmentation = "EAD")


Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$PCA
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Variance/"

Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$Variance
dev.off()



SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Loading 1/"

Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`PC1 loadings`
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Loading 2/"

Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`PC2 loadings`
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Loadings/"

Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$loadings
dev.off()



# EAD SP

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/PCA plot/"

P<-PCA(EAD_SP,phase="SP",fragmentation = "EAD")


Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$PCA
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Variance/"

Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$Variance
dev.off()



SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Loading 1/"

Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`PC1 loadings`
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Loading 2/"

Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`PC2 loadings`
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PCA/Loadings/"

Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$loadings
dev.off()


################################################################################
# PLS-DA analysis

# CID EP

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/PLSDA plot/"

P<-PLSDA(CID_EP,phase="EP",fragmentation = "CID")


Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$PLSDA
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/Variance/"

Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$Variance
dev.off()



SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/VIP 1/"

Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`VIP C1`
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/VIP 2/"

Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`VIP C2`
dev.off()






# CID SP

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/PLSDA plot/"

P<-PLSDA(CID_SP,phase="SP",fragmentation = "CID")


Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$PLSDA
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/Variance/"

Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$Variance
dev.off()



SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/VIP 1/"

Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`VIP C1`
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/VIP 2/"

Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`VIP C2`
dev.off()



# EAD EP

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/PLSDA plot/"

P<-PLSDA(EAD_EP,phase="EP",fragmentation = "EAD")


Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$PLSDA
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/Variance/"

Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$Variance
dev.off()



SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/VIP 1/"

Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`VIP C1`
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/VIP 2/"

Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`VIP C2`
dev.off()




# EAD SP

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/PLSDA plot/"

P<-PLSDA(EAD_SP,phase="SP",fragmentation = "EAD")


Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$PLSDA
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/Variance/"

Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$Variance
dev.off()



SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/VIP 1/"

Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =30 )
P$`VIP C1`
dev.off()


SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/PLSDA/VIP 2/"

Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 20,width =21 )
P$`VIP C2`
dev.off()



################################################################################
# Volcano Plot


P<-Volcano_plot(CID_EP,phase = "EP",fragmentation = "CID")

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Volcano Plots/Plot/"

Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
P$`Volcano plots`
dev.off()


P<-Volcano_plot(CID_SP,phase = "SP",fragmentation = "CID")

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Volcano Plots/Plot/"

Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
P$`Volcano plots`
dev.off()


P<-Volcano_plot(EAD_EP,phase = "EP",fragmentation = "EAD")

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Volcano Plots/Plot/"

Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
P$`Volcano plots`
dev.off()


P<-Volcano_plot(EAD_SP,phase = "SP",fragmentation = "EAD")

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Volcano Plots/Plot/"

Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
P$`Volcano plots`
dev.off()


################################################################################
# Sample clustering "Euclidean"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Sample Clustering/Unscaled/Euclidean/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Euclidean")
dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Euclidean")
dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Euclidean")
dev.off()

# CID SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Euclidean")
dev.off()



# Sample clustering "Poiss"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Sample Clustering/Unscaled/Poisson/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Poiss")
dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Poiss")
dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Poiss")
dev.off()

# CID SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Poiss")
dev.off()


### scaled data

# Sample clustering "Euclidean"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Sample Clustering/Scaled/Euclidean/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Euclidean",scale=T)
dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Euclidean",scale=T)
dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Euclidean",scale=T)
dev.off()

# CID SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Euclidean",scale=T)
dev.off()



# Sample clustering "Poiss"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Sample Clustering/Scaled/Poisson/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Poiss",scale=T)
dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Poiss",scale=T)
dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Poiss",scale=T)
dev.off()

# CID SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =29 )
sample_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Poiss",scale=T)
dev.off()

################################################################################
################################################################################
# Lipids clustering "Euclidean"
while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Lipids Clustering/Unscaled/Euclidean/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Euclidean")
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Euclidean")
while (dev.cur()>1) dev.off()
# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Euclidean")
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Euclidean")
while (dev.cur()>1) dev.off()


# Lipid clustering "Poiss"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Lipids Clustering/Unscaled/Poisson/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Poiss")
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Poiss")
while (dev.cur()>1) dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Poiss")
while (dev.cur()>1) dev.off()

# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Poiss")
while (dev.cur()>1) dev.off()



### scaled data

# Lipid clustering "Euclidean"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Lipids Clustering/Scaled/Euclidean/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Euclidean",scale=T)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Euclidean",scale=T)
while (dev.cur()>1) dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Euclidean",scale=T)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Euclidean",scale=T)
while (dev.cur()>1) dev.off()



# Lipids clustering "Poiss"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Lipids Clustering/Scaled/Poisson/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Poiss",scale=T)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Poiss",scale=T)
while (dev.cur()>1) dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Poiss",scale=T)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Poiss",scale=T)
while (dev.cur()>1) dev.off()








################################################################################
# Top N lipid clustering
################################################################################
# Lipids clustering "Euclidean"
while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Lipids Clustering/top N/Unscaled/Euclidean/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Euclidean",topN = 50)
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Euclidean",topN = 50)
while (dev.cur()>1) dev.off()
# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Euclidean",topN = 50)
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Euclidean",topN = 50)
while (dev.cur()>1) dev.off()


# Lipid clustering "Poiss"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Lipids Clustering/top N/Unscaled/Poisson/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Poiss",topN = 50)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Poiss",topN = 50)
while (dev.cur()>1) dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Poiss",topN = 50)
while (dev.cur()>1) dev.off()

# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Poiss",topN = 50)
while (dev.cur()>1) dev.off()



### scaled data

# Lipid clustering "Euclidean"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Lipids Clustering/top N/Scaled/Euclidean/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Euclidean",scale=T,topN = 50)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Euclidean",scale=T,topN = 50)
while (dev.cur()>1) dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Euclidean",scale=T,topN = 50)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Euclidean",scale=T,topN = 50)
while (dev.cur()>1) dev.off()



# Lipids clustering "Poiss"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Lipids Clustering/top N/Scaled/Poisson/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Poiss",scale=T,topN = 50)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Poiss",scale=T,topN = 50)
while (dev.cur()>1) dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Poiss",scale=T,topN = 50)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
lipid_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Poiss",scale=T,topN = 50)
while (dev.cur()>1) dev.off()













################################################################################
################################################################################
# Bi clustering "Euclidean"
while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Bi Clustering/Unscaled/Euclidean/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Euclidean")
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Euclidean")
while (dev.cur()>1) dev.off()
# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Euclidean")
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Euclidean")
while (dev.cur()>1) dev.off()


# Lipid clustering "Poiss"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Bi Clustering/Unscaled/Poisson/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Poiss")
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Poiss")
while (dev.cur()>1) dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Poiss")
while (dev.cur()>1) dev.off()

# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Poiss")
while (dev.cur()>1) dev.off()



### scaled data

# Lipid clustering "Euclidean"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Bi Clustering/Scaled/Euclidean/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Euclidean",scale=T)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Euclidean",scale=T)
while (dev.cur()>1) dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Euclidean",scale=T)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Euclidean",scale=T)
while (dev.cur()>1) dev.off()



# Bi Clustering "Poiss"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Bi Clustering/Scaled/Poisson/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Poiss",scale=T)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Poiss",scale=T)
while (dev.cur()>1) dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Poiss",scale=T)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Bi_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Poiss",scale=T)
while (dev.cur()>1) dev.off()








################################################################################
# Top N lipid clustering
################################################################################
# Bi Clustering "Euclidean"
while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Bi Clustering/top N/Unscaled/Euclidean/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Euclidean",topN = 50)
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Euclidean",topN = 50)
while (dev.cur()>1) dev.off()
# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Euclidean",topN = 50)
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Euclidean",topN = 50)
while (dev.cur()>1) dev.off()


# Lipid clustering "Poiss"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Bi Clustering/top N/Unscaled/Poisson/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Poiss",topN = 50)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Poiss",topN = 50)
while (dev.cur()>1) dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Poiss",topN = 50)
while (dev.cur()>1) dev.off()

# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Poiss",topN = 50)
while (dev.cur()>1) dev.off()



### scaled data

# Lipid clustering "Euclidean"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Bi Clustering/top N/Scaled/Euclidean/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Euclidean",scale=T,topN = 50)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Euclidean",scale=T,topN = 50)
while (dev.cur()>1) dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Euclidean",scale=T,topN = 50)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Euclidean",scale=T,topN = 50)
while (dev.cur()>1) dev.off()



# Bi Clustering "Poiss"

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Bi Clustering/top N/Scaled/Poisson/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(CID_EP,phase="EP",fragmentation="CID",mode="Poiss",scale=T,topN = 50)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(CID_SP,phase="SP",fragmentation="CID",mode="Poiss",scale=T,topN = 50)
while (dev.cur()>1) dev.off()

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(EAD_EP,phase="EP",fragmentation="EAD",mode="Poiss",scale=T,topN = 50)
while (dev.cur()>1) dev.off()

# CID SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =27 )
Bi_cluster(EAD_SP,phase="SP",fragmentation="EAD",mode="Poiss",scale=T,topN = 50)
while (dev.cur()>1) dev.off()







################################################################################
# Lipid species fold change 
################################################################################
while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Fold Change HeatMap/All/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Fold_change(CID_EP,phase="EP",fragmentation="CID")
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Fold_change(CID_SP,phase="SP",fragmentation="CID")
while (dev.cur()>1) dev.off()
# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Fold_change(EAD_EP,phase="EP",fragmentation="EAD")
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Fold_change(EAD_SP,phase="SP",fragmentation="EAD")
while (dev.cur()>1) dev.off()


# top N

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Fold Change HeatMap/top N/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Fold_change(CID_EP,phase="EP",fragmentation="CID",topN = 50)
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Fold_change(CID_SP,phase="SP",fragmentation="CID",topN = 50)
while (dev.cur()>1) dev.off()
# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Fold_change(EAD_EP,phase="EP",fragmentation="EAD",topN = 50)
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
Fold_change(EAD_SP,phase="SP",fragmentation="EAD",topN = 50)
while (dev.cur()>1) dev.off()





################################################################################
# Chain Length fold change 
################################################################################
while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Chain Length HeatMap/All/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
chain_length(CID_EP,phase="EP",fragmentation="CID")
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
chain_length(CID_SP,phase="SP",fragmentation="CID")
while (dev.cur()>1) dev.off()
# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
chain_length(EAD_EP,phase="EP",fragmentation="EAD")
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
chain_length(EAD_SP,phase="SP",fragmentation="EAD")
while (dev.cur()>1) dev.off()


# per class

while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Chain Length HeatMap/Per Class/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
chain_length(CID_EP,phase="EP",fragmentation="CID",class = T)
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
chain_length(CID_SP,phase="SP",fragmentation="CID",class = T)
while (dev.cur()>1) dev.off()
# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
chain_length(EAD_EP,phase="EP",fragmentation="EAD",class = T)
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
chain_length(EAD_SP,phase="SP",fragmentation="EAD",class = T)
while (dev.cur()>1) dev.off()




################################################################################
# Number of DB fold change 
################################################################################
while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Number of DB HeatMap/All/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
db_number(CID_EP,phase="EP",fragmentation="CID")
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
db_number(CID_SP,phase="SP",fragmentation="CID")
while (dev.cur()>1) dev.off()
# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
db_number(EAD_EP,phase="EP",fragmentation="EAD")
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
db_number(EAD_SP,phase="SP",fragmentation="EAD")
while (dev.cur()>1) dev.off()


# per class

while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Number of DB HeatMap/Per Class/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
db_number(CID_EP,phase="EP",fragmentation="CID",class = T)
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
db_number(CID_SP,phase="SP",fragmentation="CID",class = T)
while (dev.cur()>1) dev.off()
# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
db_number(EAD_EP,phase="EP",fragmentation="EAD",class = T)
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
db_number(EAD_SP,phase="SP",fragmentation="EAD",class = T)
while (dev.cur()>1) dev.off()





################################################################################
# Position of DB fold change 
################################################################################
while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/DB position HeatMap/All/"


Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
db_position(EAD_EP,phase="EP",fragmentation="EAD")
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
db_position(EAD_SP,phase="SP",fragmentation="EAD")
while (dev.cur()>1) dev.off()


# per class

while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/DB position HeatMap/Per Class/"


Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
db_position(EAD_EP,phase="EP",fragmentation="EAD",class = T)
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
db_number(EAD_SP,phase="SP",fragmentation="EAD",class = T)
while (dev.cur()>1) dev.off()



################################################################################
# Acyl Chain fold change 
################################################################################
while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Acyl Chain HeatMap/All/"


# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
sn(CID_EP,phase="EP",fragmentation="CID",size=6)
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
sn(CID_SP,phase="SP",fragmentation="CID",size=6)
while (dev.cur()>1) dev.off()


# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
sn(EAD_EP,phase="EP",fragmentation="EAD",size=6)
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
sn(EAD_SP,phase="SP",fragmentation="EAD",size=6)
while (dev.cur()>1) dev.off()


# per class

while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Acyl Chain HeatMap/Per Class/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
sn(CID_EP,phase="EP",fragmentation="CID",size=6,class=T)
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
sn(CID_SP,phase="SP",fragmentation="CID",size=6,class=T)
while (dev.cur()>1) dev.off()




# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
sn(EAD_EP,phase="EP",fragmentation="EAD",size=4,class=T)
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
sn(EAD_SP,phase="SP",fragmentation="EAD",size=4,class=T)
while (dev.cur()>1) dev.off()




################################################################################
# Cyclopropane
################################################################################
while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Cyclopropane/All/"


# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
cyclo_number(CID_EP,phase="EP",fragmentation="CID",size=14)
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
cyclo_number(CID_SP,phase="SP",fragmentation="CID",size=14)
while (dev.cur()>1) dev.off()


# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
cyclo_number(EAD_EP,phase="EP",fragmentation="EAD",size=14)
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
cyclo_number(EAD_SP,phase="SP",fragmentation="EAD",size=14)
while (dev.cur()>1) dev.off()


# per class

while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Cyclopropane/Per Class/"

# CID EP
Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
cyclo_number(CID_EP,phase="EP",fragmentation="CID",size=10,class=T)
while (dev.cur()>1) dev.off()
# CID SP
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
cyclo_number(CID_SP,phase="SP",fragmentation="CID",size=10,class=T)
while (dev.cur()>1) dev.off()




# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
cyclo_number(EAD_EP,phase="EP",fragmentation="EAD",size=10,class=T)
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =25 )
cyclo_number(EAD_SP,phase="SP",fragmentation="EAD",size=10,class=T)
while (dev.cur()>1) dev.off()


################################################################################
# SN attachement of cyclopropane
################################################################################

# Fold Change (ALL)

while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Cyclopropane sn position/FC/All/"
# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =32 )
cyclo_sn(EAD_EP,phase="EP",fragmentation="EAD",size=14,FC=T,class=F)
while (dev.cur()>1) dev.off()
# EAD SP  
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =32 )
cyclo_sn(EAD_SP,phase="SP",fragmentation="EAD",size=14,FC=T,class=F)
while (dev.cur()>1) dev.off()

# Fold Change (Per Class)

while (dev.cur()>1) dev.off()

SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Cyclopropane sn position/FC/Per Class/"

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =32 )
cyclo_sn(EAD_EP,phase="EP",fragmentation="EAD",size=14,FC=T,class=T)
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =32 )
cyclo_sn(EAD_SP,phase="SP",fragmentation="EAD",size=14,FC=T,class=T)
while (dev.cur()>1) dev.off()


# Absolute (ALL)

while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Cyclopropane sn position/Absolute/All/"

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =32 )
cyclo_sn(EAD_EP,phase="EP",fragmentation="EAD",size=14,FC=F,class=F)
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =32 )
cyclo_sn(EAD_SP,phase="SP",fragmentation="EAD",size=14,FC=F,class=F)
while (dev.cur()>1) dev.off()


# Absolute (Per class)

while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Cyclopropane sn position/Absolute/Per Class/"

# EAD EP
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =32 )
cyclo_sn(EAD_EP,phase="EP",fragmentation="EAD",size=14,FC=F,class=T)
while (dev.cur()>1) dev.off()
# EAD SP
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 16,width =32 )
cyclo_sn(EAD_SP,phase="SP",fragmentation="EAD",size=14,FC=F,class=T)
while (dev.cur()>1) dev.off()









################################################################################
# Barplot FC per class
################################################################################

# CID EP

while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Barplot FC per Class/"

Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
FC_class_bp(CID_EP,phase="EP",fragmentation = "CID")
while (dev.cur()>1) dev.off()


# CID SP

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
FC_class_bp(CID_SP,phase="SP",fragmentation = "CID")
while (dev.cur()>1) dev.off()

# EAD EP

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
FC_class_bp(EAD_EP,phase="EP",fragmentation = "EAD")
while (dev.cur()>1) dev.off()

# EAD SP

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 29.7,width =21 )
FC_class_bp(EAD_SP,phase="SP",fragmentation = "EAD")
while (dev.cur()>1) dev.off()





################################################################################
# Chord Diagram
################################################################################

# CID EP

while (dev.cur()>1) dev.off()
SAVE <- "W:/users/Abraham/Exp_004_version_2/Plots/3) MS-dial data analysis/Differential Analysis/Chord Diagram/"

Cairo(file=paste0(SAVE,"CID_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 21,width =21 )
Chord_plot(CID_EP,phase="EP",fragmentation = "CID")
while (dev.cur()>1) dev.off()


# CID SP

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"CID_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 21,width =21 )
Chord_plot(CID_SP,phase="SP",fragmentation = "CID")
while (dev.cur()>1) dev.off()

# EAD EP

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"EAD_EP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 21,width =21 )
Chord_plot(EAD_EP,phase="EP",fragmentation = "EAD")
while (dev.cur()>1) dev.off()

# EAD SP

while (dev.cur()>1) dev.off()
Cairo(file=paste0(SAVE,"EAD_SP.png"),
      type="png",bg="white",dpi=300,units = "cm", height = 21,width =21 )
Chord_plot(EAD_SP,phase="SP",fragmentation = "EAD")
while (dev.cur()>1) dev.off()



