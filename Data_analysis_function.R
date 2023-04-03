################################################################################
#
# Abraham Moyal
# 16.02.23
# Data Analysis
# 1) Distribution of peaks
# 2) Lipid class distribution
# 3) Number of lipids identified with EAD and CID
# 4) Percetage of identificitaion for each level
# 5) Number of isomers in EAD and CID
# 6) Number of resolved isomers by EAD
# 7)
################################################################################
#
# Load libraries
library(circlize)
library(huxtable)
library(matrixStats)
library("mixOmics")
library(pheatmap)
library(PoiClaClu)
library(factoextra)
library(ggplot2)
library(dplyr)
library(ggVennDiagram)
library(tidyverse)
library(rcdk)
library(cluster)
library(svDialogs)
library(plotly)
library(rsconnect)
library(stringi)
library(insight)
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library(VennDiagram)
library(RColorBrewer)
library(Cairo)
library(ggplot2)
library(gridExtra)
library(data.table)
library(dendsort)
library(PoiClaClu)

################################################################################

# function 1: distribution of lipids
distribution<-function(dataframe,type="intensity",fragmentation="CID",phase="EP"){
  conditions <- c(expression(Delta*pgpB),
                  expression(Delta*cdh),
                  expression(Delta*cfa),
                  expression(Delta*plsX),
                  expression(Delta*opgB),
                  expression(Delta*clsC),
                  expression(Delta*pgpC),
                  expression(Delta*pldC),
                  expression(Delta*clsB),
                  expression(Delta*clsA),
                  expression(Delta*pldB),
                  expression(Delta*pldA),
                  expression(Delta*pgpA),
                  expression(Delta*aas),
                  expression(Delta*dgkA),
                  "Wilt Type"
  )
  if(phase=="EP"){
    Intensity_data<-dataframe[,c(36:80,84:86)]
    
  }else if(phase=="SP"){
    Intensity_data<-dataframe[,c(87:131,135:137)]
  }
  if(fragmentation=="CID"){
    color<-rep(c(rep("coral1",3),rep("coral3",3)),8)
  }else if(fragmentation=="EAD"){
    color<-rep(c(rep("palegreen3",3),rep("palegreen1",3)),8)
    
  }
  main= paste0(fragmentation," ",phase,": ",type," ","distribution")
  yy <- bquote(log[10](.(type) ))
  
  
  
  Intensity_data[Intensity_data<=0]<- NA
  
  Intensity_data <- log10(matrix(as.numeric(as.matrix(Intensity_data)), ncol=ncol(Intensity_data)))
  
  layout(matrix(c(1,1),ncol=2),width=c(8))
  par(yaxs="r")
  par(mar=c(6.5,5,5,2))
  boxplot(Intensity_data,col=color,las=2,xaxt="n" ,main=main,cex.main=1.8)
  axis(1,at=seq(1,48,3)+1,labels = conditions,las=2,cex.axis=1.4 )
  mtext(2,line=2,text = yy,cex=1.5)
  
  
}



################################################################################

# function 2: Distribution of lipid class in each sample


class_distribution<-function(dataframe,fragmentation="CID",phase="EP"){
  color <-c("azure4","blue4","cyan3","darkmagenta",
            "darkorange","coral1","cyan4","tomato3",
            "chartreuse4","red1","yellow3","black",
            "tan","wheat2","chartreuse1","steelblue1")
  
  l_class<-c("LPE","LPC","DG","Cer_NDS",
             "PE","CoQ","EtherPC","EtherDG",
             "PG", "Cer_NS","PS","PC",
             "TG","NAE","EtherTG","AHexCer")
  conditions <- c(expression(Delta*pgpB),
                  expression(Delta*cdh),
                  expression(Delta*cfa),
                  expression(Delta*plsX),
                  expression(Delta*opgB),
                  expression(Delta*clsC),
                  expression(Delta*pgpC),
                  expression(Delta*pldC),
                  expression(Delta*clsB),
                  expression(Delta*clsA),
                  expression(Delta*pldB),
                  expression(Delta*pldA),
                  expression(Delta*pgpA),
                  expression(Delta*aas),
                  expression(Delta*dgkA),
                  "Wilt Type"
  )
  if(phase=="EP"){
    Intensity_data<-dataframe[,c(36:80,84:86)]
  }else if(phase=="SP"){
    Intensity_data<-dataframe[,c(87:131,135:137)]
  }
  
  Intensity_data[Intensity_data<0]<-0
  class<-sort(unique(dataframe$Ontology))
  
  count<-matrix(ncol=48,nrow=length(class))
  rownames(count) <- class
  for(i in 1:length(class)){
    name<-class[i]
    index <- which(dataframe$Ontology==name)
    
    for(j in 1:48){
      count[i,j]<- sum(as.numeric(Intensity_data[index,j]))
    }
  }
  
  for(i in 1:48){
    count[,i] <- count[,i]/sum(count[,i])*100
  }
  
  color_s<-rep(NA,nrow(count))
  for(i in 1:nrow(count)){
    name<-rownames(count)[i]
    index<-which(l_class==name)
    color_s[i] <- color[index]
  }
  main= paste0(fragmentation," ",phase,": ","Lipid Class Distribution")
  yy <- "Percentage [%]"
  
  space<-rep(0,48)
  space[seq(1,45,3)+3]<-0.4
  
  color<-color_s
  layout(matrix(c(1,2),ncol=2),width=c(8,2))
  par(yaxs="r")
  par(mar=c(6.5,5,5,0))
  x<- barplot(count,col=color,las=2,xaxt="n" ,main=main,cex.main=1.8,xaxt="n",space=space)
  axis(1,at=x[seq(1,48,3)+1],labels = conditions,las=2,cex.axis=1.4 )
  mtext(2,line=2.4,text = yy,cex=1.5)
  par(mar=c(9,0,5,0))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')  
  legend("topleft",col=color,pch=15,legend=rownames(count),bty="n",pt.cex = 1.5)
}


################################################################################

# function 3: Amount of lipids identified


lipid_count<-function(CID,EAD,phase){
  c<-unique(CID$`Metabolite name`)
  e<-unique(EAD$`Metabolite name`)
  y<-c(nrow(CID),nrow(EAD),length(c),length(e))
  par(mfrow=c(1,1),mar=c(5,5,5,5))
  x<-barplot(y, 
             names.arg = c("CID","EAD","CID","EAD"),
             main=paste0("Number of Lipids identified (",phase,")"),
             col = c("coral1","palegreen3"),cex.main=1.8,cex.names  =1.6,
             ylab = "# of Lipids",cex.lab =1.5,ylim = c(0,max(y)+5))
  text(x,y,y,pos=1,cex=1.5)
  abline(v=2.5,lwd=1.5,lty=c("19"))
  mtext(3,at=c(1.3,3.7),text=c("Feature Level","Lipid ID level"),cex = 1.4)
  
}

################################################################################

# function 4: Lipid characteristic

#   EAD characteristics 
#----------------------



lipid_EAD <- function(dataframe){
  
  lipid_names <- dataframe$`Metabolite name`
  
  caracteristics<- matrix(ncol=6,nrow=nrow(dataframe))
  
  
  colnames(caracteristics) <- c("Lipid class","Sum composition","Fatty acid identification","Fatty acid sn position","Double bond position","Stereochemistry")
  for(i in 1:length(lipid_names)){
    
    x<-lipid_names[i]
    if(grepl("|", x)==T){
      x <- gsub(".*\\|", "",x )
    }else{
      x <- x  
    }
    
    
    caracteristics[i,1] <-gsub("\\ .*", "",x )
    
    y <- gsub(".*\\ ", "",x )
    y <- gsub("O-","",y)
    
    
    if(grepl(":",y)==T){
      caracteristics[i,2] <- x
    }
    
    if(grepl("/",y)==T | grepl("_",y)==T){
      caracteristics[i,3] <- x
    }
    
    if(grepl("/",y)==T){
      caracteristics[i,4] <- x
    }
    
    if(grepl('\\(',y)==T){
      caracteristics[i,5] <- x
    }
    
    
    if(grepl("Z",y)==T | grepl("E",y)==T ){
      caracteristics[i,6] <- x
    }
    
    
  }
  caracteristics <- as.data.frame(caracteristics)
  
  return(caracteristics)
}






lipid_CID <- function(dataframe){
  
  lipid_names <- dataframe$`Metabolite name`
  
  caracteristics<- matrix(ncol=6,nrow=nrow(dataframe))
  
  
  colnames(caracteristics) <- c("Lipid class","Sum composition","Fatty acid identification","Fatty acid sn position","Double bond position","Stereochemistry")
  i=9
  for(i in 1:length(lipid_names)){
    
    x<-lipid_names[i]
    if(grepl("|", x)==T){
      x <- gsub(".*\\|", "",x )
    }else{
      x <- x  
    }
    
    
    caracteristics[i,1] <-gsub("\\ .*", "",x )
    
    y <- gsub(".*\\ ", "",x )
    y <- gsub("O-","",y)
    
    
    if(grepl(":",y)==T){
      caracteristics[i,2] <- x
    }
    
    if(grepl("/",y)==T | grepl("_",y)==T){
      caracteristics[i,3] <- x
    }
    
    if(grepl("/",y)==T){
      caracteristics[i,4] <- x
    }
    
    if(grepl('\\(',y)==T){
      caracteristics[i,5] <- x
    }
    
    
    if(grepl("Z",y)==T | grepl("E",y)==T ){
      caracteristics[i,6] <- x
    }
    
    
  }
  caracteristics <- as.data.frame(caracteristics)
  
  return(caracteristics)
}

lipid_depth  <- function(CID,EAD,phase){
  l_c <- lipid_CID(CID)
  l_e <- lipid_EAD(EAD)
  
  
  count<-matrix(nrow=6,ncol=2)
  rownames(count)<- c("Lipid class","Sum composition","Fatty acid identification","Fatty acid sn position","Double bond position","Stereochemistry")
  colnames(count)<-c("CID","EAD")
  
  for(i in 1:6){
    count[i,1] <- length(which(is.na(l_c[,i])==F))
    count[i,2] <- length(which(is.na(l_e[,i])==F))
  }
  
  
  
  count[,1]<- count[,1]/count[1,1]*100
  count[,2]<- count[,2]/count[1,2]*100
  
  color<-c("forestgreen","springgreen4","palegreen4","palegreen3","palegreen2","palegreen1",
           "red4","red2","red1","tomato3","tomato1","coral1") 
  
  par(mfrow=c(1,1),mar=c(11,5,5,5))
  x<-barplot(count,beside=T,ylab="percentage [%]",cex.lab=1.5,names.arg =c("",""),col=color,ylim=c(0,110)
  )
  abline(v=7.5,lwd=2,lty=2)
  mtext(3,at=c(3.5,10.5),text = c("CID","EAD"),line = 0.5,cex = 1.2)
  mtext(3,at=c(7.5),text = paste0("Lipid Level Identification (",phase,")"),line = 1.6,cex = 1.5)
  axis(1,at=unlist(x),rep(rownames(count),2),las=2)
  text(unlist(x),
       unlist(count),
       paste0(round(unlist(count),digits=1)," %"),
       pos=3,offset = 0.1,cex=0.9
  )
}

################################################################################

# function 5: Lipid isomers

#   EAD characteristics 
#----------------------
characteristics_EAD <- function(dataframe){
  
  lipid_names <- dataframe$`Metabolite name`
  
  caracteristics<- matrix(ncol=11,nrow=nrow(dataframe))
  
  
  colnames(caracteristics) <- c("Subclass level","Species level", "Molecular Species level","sn position level"," DBE position", "Full structure","Chain length","# DB","cyclopropane","DB positions","Chains")
  
  caracteristics[,9]<- rep("NO",nrow(dataframe))
  
  for(i in 1:length(lipid_names)){
    
    x<-lipid_names[i]
    if(grepl("|", x)==T){
      x <- gsub(".*\\|", "",x )
    }else{
      x <- x  
    }
    
    
    caracteristics[i,1] <- l1<-gsub("\\ .*", "",x )
    
    y <- gsub(".*\\ ", "",x )
    y <- gsub("O-","",y)
    if(grepl("/",y)==T ){
      
      FA<-strsplit(y,split = '/')[[1]]
      
      db<-strsplit(y,split = '/')[[1]]
      a<-grepl('\\(',db)
      db<-db[a]
      chain<-gsub("\\(.*", "",db )
      
      db<-gsub(".*\\(", "",db )
      
      db<-gsub(")","",db)
      db<-gsub("Z","",db)
      db<-gsub("E","",db)
      db<-paste(db,collapse = ",")
      
      
      caracteristics[i,10] <- db
      FA <- gsub("\\(.*", "",FA )
      
      caracteristics[i,11] <- paste(FA,collapse = "_")
      c<-paste(FA,collapse = "/")
      
      caracteristics[i,3]<-paste0(l1," ",c) 
      
      a<-as.numeric(gsub("\\:.*", "",FA ))
      b<-as.numeric(gsub(".*\\:", "",FA ) )
      
      caracteristics[i,2] <- paste0(l1," ",sum(a),":",sum(b))
      caracteristics[i,7] <- sum(a)
      caracteristics[i,8] <- sum(b)
      
      test <- as.numeric(a)  %% 2 == 0 
      if(all(test)!=T){
        caracteristics[i,9] <- "Yes"
      }
      
    }else if(grepl("/",y)==F & grepl("_",y)== F & grepl(" ",x)==T){
      
      FA<-strsplit(y,split = '/')[[1]]
      
      db<-strsplit(y,split = '/')[[1]]
      a<-grepl('\\(',db)
      db<-db[a]
      gsub("\\(.*", "",db )
      
      db<-gsub(".*\\(", "",db )
      
      db<-gsub(")","",db)
      db<-gsub("Z","",db)
      db<-gsub("E","",db)
      db<-paste(db,collapse = ",")
      
      
      caracteristics[i,10] <- db
      FA <- gsub("\\(.*", "",FA )
      
      caracteristics[i,11] <- paste(FA,collapse = "_")
      c<-paste(FA,collapse = "/")
      
      caracteristics[i,3]<-paste0(l1," ",c) 
      
      a<-as.numeric(gsub("\\:.*", "",FA ))
      b<-as.numeric(gsub(".*\\:", "",FA ) )
      
      caracteristics[i,2] <- paste0(l1," ",sum(a),":",sum(b))
      caracteristics[i,7] <- sum(a)
      caracteristics[i,8] <- sum(b)
      
      test <- as.numeric(a)  %% 2 == 0 
      if(all(test)!=T){
        caracteristics[i,9] <- "Yes"
      }
      
    }else if(grepl("_",y)==T ){
      
      FA<-strsplit(y,split = '_')[[1]]
      
      db<-strsplit(y,split = '_')[[1]]
      a<-grepl('\\(',db)
      db<-db[a]
      gsub("\\(.*", "",db )
      
      db<-gsub(".*\\(", "",db )
      
      db<-gsub(")","",db)
      db<-gsub("Z","",db)
      db<-gsub("E","",db)
      db<-paste(db,collapse = ",")
      
      
      caracteristics[i,10] <- db
      FA <- gsub("\\(.*", "",FA )
      caracteristics[i,11] <- paste(FA,collapse = "/")
      c<-paste(FA,collapse = "_")
      
      caracteristics[i,3]<-paste0(l1," ",c) 
      
      a<-as.numeric(gsub("\\:.*", "",FA ))
      b<-as.numeric(gsub(".*\\:", "",FA ) )
      
      caracteristics[i,2] <- paste0(l1," ",sum(a),":",sum(b))
      caracteristics[i,7] <- sum(a)
      caracteristics[i,8] <- sum(b)
      
      test <- as.numeric(a)  %% 2 == 0 
      if(all(test)!=T){
        caracteristics[i,9] <- "Yes"
      }
    }
    
    
    if(grepl("/",x)==T){
      caracteristics[i,4] <- x
    } 
    
    if(grepl('\\(',x)==T){
      caracteristics[i,5] <- x
    }
    
    if(grepl("E)",x)==T | grepl("Z)",x)==T){
      caracteristics[i,6] <- x
    }
    
    
    
  }
  caracteristics <- as.data.frame(caracteristics)
  
  return(caracteristics)
}




#   CID characteristics 
#----------------------

characteristics_CID <- function(dataframe){
  
  lipid_names <- dataframe$`Metabolite name`
  
  caracteristics<- matrix(ncol=11,nrow=nrow(dataframe))
  
  
  colnames(caracteristics) <- c("Subclass level","Species level", "Molecular Species level","sn position level"," DBE position", "Full structure","Chain length","# DB","cyclopropane","DB positions","Chains")
  
  caracteristics[,9]<- rep("NO",nrow(dataframe))
  
  for(i in 1:length(lipid_names)){
    
    x<-lipid_names[i]
    if(grepl("|", x)==T){
      x <- gsub(".*\\|", "",x )
    }else{
      x <- x  
    }
    
    
    caracteristics[i,1] <- l1<-gsub("\\ .*", "",x )
    
    y <- gsub(".*\\ ", "",x )
    y <- gsub("O-","",y)
    if(grepl("/",y)==T ){
      
      FA<-strsplit(y,split = '/')[[1]]
      
      db<-strsplit(y,split = '/')[[1]]
      a<-grepl('\\(',db)
      db<-db[a]
      chain<-gsub("\\(.*", "",db )
      
      db<-gsub(".*\\(", "",db )
      
      db<-gsub(")","",db)
      db<-gsub("Z","",db)
      db<-gsub("E","",db)
      db<-paste(db,collapse = ",")
      
      
      caracteristics[i,10] <- db
      FA <- gsub("\\(.*", "",FA )
      
      caracteristics[i,11] <- paste(FA,collapse = "_")
      c<-paste(FA,collapse = "/")
      
      caracteristics[i,3]<-paste0(l1," ",c) 
      
      a<-as.numeric(gsub("\\:.*", "",FA ))
      b<-as.numeric(gsub(".*\\:", "",FA ) )
      
      caracteristics[i,2] <- paste0(l1," ",sum(a),":",sum(b))
      caracteristics[i,7] <- sum(a)
      caracteristics[i,8] <- sum(b)
      
      test <- as.numeric(a)  %% 2 == 0 
      if(all(test)!=T){
        caracteristics[i,9] <- "Yes"
      }
      
    }else if(grepl("/",y)==F & grepl("_",y)== F & grepl(" ",x)==T){
      
      FA<-strsplit(y,split = '/')[[1]]
      
      db<-strsplit(y,split = '/')[[1]]
      a<-grepl('\\(',db)
      db<-db[a]
      gsub("\\(.*", "",db )
      
      db<-gsub(".*\\(", "",db )
      
      db<-gsub(")","",db)
      db<-gsub("Z","",db)
      db<-gsub("E","",db)
      db<-paste(db,collapse = ",")
      
      
      caracteristics[i,10] <- db
      FA <- gsub("\\(.*", "",FA )
      
      caracteristics[i,11] <- paste(FA,collapse = "_")
      c<-paste(FA,collapse = "/")
      
      caracteristics[i,3]<-paste0(l1," ",c) 
      
      a<-as.numeric(gsub("\\:.*", "",FA ))
      b<-as.numeric(gsub(".*\\:", "",FA ) )
      
      caracteristics[i,2] <- paste0(l1," ",sum(a),":",sum(b))
      caracteristics[i,7] <- sum(a)
      caracteristics[i,8] <- sum(b)
      
      test <- as.numeric(a)  %% 2 == 0 
      if(all(test)!=T){
        caracteristics[i,9] <- "Yes"
      }
      
    }else if(grepl("_",y)==T ){
      
      FA<-strsplit(y,split = '_')[[1]]
      
      db<-strsplit(y,split = '_')[[1]]
      a<-grepl('\\(',db)
      db<-db[a]
      gsub("\\(.*", "",db )
      
      db<-gsub(".*\\(", "",db )
      
      db<-gsub(")","",db)
      db<-gsub("Z","",db)
      db<-gsub("E","",db)
      db<-paste(db,collapse = ",")
      
      
      caracteristics[i,10] <- db
      FA <- gsub("\\(.*", "",FA )
      caracteristics[i,11] <- paste(FA,collapse = "/")
      c<-paste(FA,collapse = "_")
      
      caracteristics[i,3]<-paste0(l1," ",c) 
      
      a<-as.numeric(gsub("\\:.*", "",FA ))
      b<-as.numeric(gsub(".*\\:", "",FA ) )
      
      caracteristics[i,2] <- paste0(l1," ",sum(a),":",sum(b))
      caracteristics[i,7] <- sum(a)
      caracteristics[i,8] <- sum(b)
      
      test <- as.numeric(a)  %% 2 == 0 
      if(all(test)!=T){
        caracteristics[i,9] <- "Yes"
      }
    }
    
    
    if(grepl("/",x)==T){
      caracteristics[i,4] <- x
    } 
    
    if(grepl('\\(',x)==T){
      caracteristics[i,5] <- x
    }
    
    if(grepl("E)",x)==T | grepl("Z)",x)==T){
      caracteristics[i,6] <- x
    }
    
    
    
  }
  caracteristics <- as.data.frame(caracteristics)
  
  return(caracteristics)
}

lipid_isomers<-function(CID,EAD,phase){
  l_c<-characteristics_CID(CID)
  l_e<-characteristics_EAD(EAD)
  
  iso_CID <- unique(l_c$`Species level` )
  iso_EAD <- unique(l_e$`Species level`)
  
  i_c <- matrix(ncol=1,nrow=length(iso_CID))
  rownames(i_c) <- iso_CID
  
  i_e <- matrix(ncol=1,nrow=length(iso_EAD))
  rownames(i_e) <- iso_EAD
  
  
  for(i in 1:length(iso_CID)){
    name<-iso_CID[i]
    index<-which(l_c$`Species level`==name)
    i_c[i,1] <- length(index)
  }
  
  for(i in 1:length(iso_EAD)){
    name<-iso_EAD[i]
    index<-which(l_e$`Species level`==name)
    i_e[i,1] <- length(index)
  }
  
  u_c<- as.matrix(table(i_c)[-1])
  u_e<- as.matrix(table(i_e)[-1])
  
  u <- sort(unique(as.numeric(c(rownames(u_c),rownames(u_e)))))
  
  
  bar<-matrix(0,ncol=2,nrow=length(u))
  colnames(bar)<-c("CID","EAD")
  rownames(bar)<-u
  
  
  for(i in 1:length(u)){
    name<-as.character(u[i])
    index_c<-which(rownames(u_c)==name)
    if(length(index_c)>0){
      bar[i,1]<-u_c[index_c]
    }
    
    index_e<-which(rownames(u_e)==name) 
    
    if(length(index_e)>0){
      bar[i,2]<-u_e[index_e]
    }
    
  }
  
  bar[,1]<-bar[,1]/sum(bar[,1])*100
  bar[,2]<-bar[,2]/sum(bar[,2])*100
  
  x<-barplot(bar,beside = T,col=c(rep("palegreen3",length(u)),rep("coral1",length(u))),
             names.arg = c("",""),ylim=c(0,60),ylab="Percentage [%]",cex.lab=1.5)
  abline(v=9.5,lwd=2,lty=2)
  mtext(3,at=c(5,14),text = c("CID","EAD"),line = 0,cex = 1.2)
  mtext(3,at=c(9.5),text = paste0("Number of Sum Composition Isomers (",phase,")"),line = 1.6,cex = 1.5)
  axis(1,at=unlist(x),rep(u,2),las=1)
  text(unlist(x)+0.1,
       unlist(bar),
       paste0(round(unlist(bar),digits=1)," %"),
       pos=3,offset = 0.1,cex=0.9
  )
  mtext(1,at=9.5,text="Number of isomeric compounts for each sum composition",cex=1.6,line=2.5)
  
}




lipid_isomers_identified<-function(CID,EAD,phase){
  l_c<-CID$`Metabolite name`
  l_e<-EAD$`Metabolite name`
  
  iso_CID <- unique(l_c)
  iso_EAD <- unique(l_e)
  
  i_c <- matrix(ncol=1,nrow=length(iso_CID))
  rownames(i_c) <- iso_CID
  
  i_e <- matrix(ncol=1,nrow=length(iso_EAD))
  rownames(i_e) <- iso_EAD
  
  
  
  for(i in 1:length(iso_CID)){
    name<-iso_CID[i]
    index<-which(l_c==name)
    i_c[i,1] <- length(index)
  }
  
  
  for(i in 1:length(iso_EAD)){
    name<-iso_EAD[i]
    index<-which(l_e==name)
    i_e[i,1] <- length(index)
  }
  
  
  u_c<- as.matrix(table(i_c))
  u_e<- as.matrix(table(i_e))
  
  u <- sort(unique(as.numeric(c(rownames(u_c),rownames(u_e)))))
  
  
  bar<-matrix(0,ncol=2,nrow=length(u))
  colnames(bar)<-c("CID","EAD")
  rownames(bar)<-u
  
  
  for(i in 1:length(u)){
    name<-as.character(u[i])
    index_c<-which(rownames(u_c)==name)
    if(length(index_c)>0){
      bar[i,1]<-u_c[index_c]
    }
    
    index_e<-which(rownames(u_e)==name) 
    
    if(length(index_e)>0){
      bar[i,2]<-u_e[index_e]
    }
    
  }
  
  bar[,1]<-bar[,1]/sum(bar[,1])*100
  bar[,2]<-bar[,2]/sum(bar[,2])*100
  x<-barplot(bar,beside = T,col=c(rep("palegreen3",length(u)),rep("coral1",length(u))),
             names.arg = c("",""),ylim=c(0,100),ylab="Percentage [%]",cex.lab=1.5)
  abline(v=7.5,lwd=2,lty=2)
  mtext(3,at=c(3.5,10.5),text = c("CID","EAD"),line = 0,cex = 1.2)
  mtext(3,at=c(7.5),text = paste0("Number of Isomers at Highest Resolved Structure Level (",phase,")"),line = 1.6,cex = 1.5)
  axis(1,at=unlist(x),rep(u,2),las=1)
  text(unlist(x)+0.1,
       unlist(bar),
       paste0(round(unlist(bar),digits=1)," %"),
       pos=3,offset = 0.1,cex=0.9
  )
  mtext(1,at=7.5,text="Number of isomeric compounts for each Lipid species",cex=1.6,line=2.5)
  
}


resolved_isomers <-function(CID,EAD,Directory,phase){
  setwd(Directory)
  
  
  
  c<- characteristics_CID(CID)
  e<- characteristics_EAD(EAD)
  
  CID$'sum product'<-c[,2]
  EAD$'sum product'<-e[,2]
  d_c<-CID[,c(4,221,2,3,5)]
  d_e<-EAD[,c(4,219,2,3,5)]
  
  
  ###### sum composition match
  
  sum_composition<- intersect(d_c$`sum product`,d_e$`sum product` )
  
  d_c_match<-d_c[which(is.na(match(d_c$`sum product`,sum_composition))==F),]
  d_e_match<-d_e[which(is.na(match(d_e$`sum product`,sum_composition))==F),]
  
  
  ###### Only lipids with different isomers
  j=c()
  for(i in 1:nrow(d_c_match)){
    name <- d_c_match$`sum product`[i]
    index<-length(which(d_c_match$`sum product`==name))
    
    if(index==1){
      j<-c(j,i)
    }
  }
  d_c_match<-d_c_match[-j,]
  
  
  
  j=c()
  for(i in 1:nrow(d_e_match)){
    name <- d_e_match$`sum product`[i]
    index<-length(which(d_e_match$`sum product`==name))
    
    if(index==1){
      j<-c(j,i)
    }
  }
  d_e_match<-d_e_match[-j,]
  
  
  d_e_match<-d_e_match[-which(is.na(d_e_match$`sum product`) ),]
  d_c_match<-d_c_match[-which(is.na(d_c_match$`sum product`) ),]
  
  
  match<-matrix(ncol=2,nrow=0)
  colnames(match)<-c("CID","EAD")
  par(mar=c(5,5,5,5))
  plot(d_c_match$`Average Rt(min)`,d_c_match$`Average Mz`,main=paste0("Features Overlap (",phase,")"),
       pch=16,cex=1.3,col="palegreen3",ylab="mz",xlab="RT (min)",cex.lab=1.5,cex.main=1.8)
  points(d_e_match$`Average Rt(min)`,d_e_match$`Average Mz`,col="coral1",cex=1.4,lwd=2)
  legend("topleft",legend=c("CID","EAD"),col=c("palegreen3","coral1"),pch=c(16,1),bty="n",pt.cex = c(1.8,1.8),cex = 1.3)
  
  delta_mz<-0.01
  delta_rt<-0.01
  for(i in 1:nrow(d_c_match) ){
    ref <- d_c_match[i,]
    for(j in 1:nrow(d_e_match)){
      test <- d_e_match[j,]
      
      if(ref$`sum product`== test$`sum product` ){
        N<-T
      }else{N<-F}
      
      if(abs(as.numeric(ref$`Average Rt(min)`)- as.numeric(test$`Average Rt(min)`))<=delta_rt  ){
        RT<-T
      }else{RT<-F}
      
      if(abs(as.numeric(ref$`Average Mz`)- as.numeric(test$`Average Mz`))<=delta_mz  ){
        MZ<-T
      }else{MZ<-F}
      
      if((ref$`Adduct type`== test$`Adduct type`)==T){
        ad<-T
      }else{ad<-F}
      
      if(all(c(N,RT,MZ,ad))==T){
        add<-matrix(c(i,j),byrow=T,ncol=2)
        colnames(add)<-c("CID","EAD")
        match<-rbind(match,add)
      }
    }
  }
  
  
  
  
  
  match_id<-cbind(d_c_match$`sum product`[match[,1]],d_c_match$`Metabolite name`[match[,1]],d_e_match$`Metabolite name`[match[,2]])
  
  match_id<-as.data.frame(match_id)
  m<-(unique(match_id[duplicated(match_id[,2]),2]))
  m_CID<-match_id %>% filter(V2 %in% m)
  
  resolved<-matrix(ncol=(1+max(table(match_id$V2))),nrow=length(m))
  colnames(resolved)<-c("CID",rep("EAD",max(table(match_id$V2))))
  resolved[,1]<-m
  for(i in 1:length(m)){
    name<-m[i]
    index<-which(m_CID[,2]==name)
    
    eee<- length(unique(m_CID[index,3]))
    if(eee>1){
      resolved[i,2:(eee+1)]<-unique(m_CID[index,3])
    }else{resolved[i,2]<-unique(m_CID[index,3])}
    
  }
  
  na_cols <- which(apply(resolved, 2, function(x) all(is.na(x))))
  resolved<-resolved[,-na_cols]
  
  write.csv(resolved,paste0("Isomers_",phase,".csv"))
  
  return(c(nrow(resolved),nrow(m_CID)))
}


venn<-function(CID,EAD,Directory){
  c<- characteristics_CID(CID)
  e<- characteristics_EAD(EAD)
  
  l_c <- na.omit(c$`Species level`)
  l_e <-na.omit(e$`Species level`)
  
  myCol <- brewer.pal(3, "Pastel2")[1:2]
  setwd(Directory)
  venn.diagram(
    x=list(l_c,l_e),
    category.names=c("CID","EAD"),filename = "venn_diagram.png",
    output=T,
    
    wd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = 1.4,
    fontface = "bold",
    fontfamily = "sans",
    
    # # Set names
    cat.cex = 1.8,
    cat.fontface = "bold",
    
    cat.fontfamily = "sans",
  )  
  
  
  
}
################################################################################
overlap<-function(CID,EAD,Directory,mz_tol=0.05,rt_tol=0.1){
  EAD<-EAD[,c(1,2,3,5,8)]
  
  master_peak_list <- CID[,c(1,2,3,5,8)]
  master_peak_list$CID<-T
  master_peak_list$EAD<-F
  mz_tol<-mz_tol
  rt_tol<-rt_tol #min
  i=5
  for(i in 1:nrow(EAD)){
    RT<- as.numeric(EAD[i,2]) - as.numeric(master_peak_list[,2])
    RT<- abs(RT) < rt_tol
    MZ<- as.numeric(EAD[i,3]) - as.numeric(master_peak_list[,3])
    MZ <- abs(MZ) < mz_tol
    Adduct <- EAD[i,4]==master_peak_list[,4]
    
    score<-data.frame(RT=RT,MZ=MZ,Adduct=Adduct)
    score[score==T]<-1
    
    index<-which(rowSums(score)==3)
    
    if(length(index)>1){
      RT<- as.numeric(EAD[i,2]) - as.numeric(master_peak_list[,2])
      minimum <- which(RT==min(RT[index]))
      index <- minimum
    }
    
    if(length(index)!=0){
      master_peak_list[index,7]<-T
    }else if(length(index)==0){
      add<-matrix(c(EAD[i,],FALSE,TRUE),nrow=1)
      colnames(add) <- colnames(master_peak_list)
      master_peak_list <- rbind(master_peak_list,add)
      master_peak_list <- data.frame(master_peak_list)
    }
    
  }
  
  
  
  master_peak_list$CID <- paste0(master_peak_list$CID, c(1:nrow(master_peak_list)))
  master_peak_list$EAD <- paste0(master_peak_list$EAD, c(1:nrow(master_peak_list)))
  
  
  C<-master_peak_list$CID
  E<-master_peak_list$EAD
  
  C<-C[grep("TRUE",C)]
  E<-E[grep("TRUE",E)]

    myCol <- brewer.pal(3, "Pastel2")[1:2]
    venn.diagram(
      x=list(C,E),
      category.names=c("CID","EAD"),filename = "venn_diagram.png",
      output=T,
      
      wd = 2,
      lty = 'blank',
      fill = myCol,
      
      # Numbers
      cex = 1.4,
      fontface = "bold",
      fontfamily = "sans",
      
      # # Set names
      cat.cex = 1.8,
      cat.fontface = "bold",
      
      cat.fontfamily = "sans",
    )  
    
  
}



################################################################################
################################################################################
# Differential Analysis
################################################################################
# function scale

scaling<-function(data){
  S<-data
  for(i in 1:nrow(data)){
    d<-as.numeric(data[i,])
    min_i<-min(as.numeric(data[i,]))
    max_i<-max(as.numeric(data[i,]))
    
    ss <- (d - min_i)/(max_i - min_i)
    
    S[i,]<- as.numeric( ss)
    
  }
  return(S)
}


################################################################################
# function 1: PCA Analysis
PCA<- function(data,phase,fragmentation){
  title <- paste0(fragmentation," data (",phase,")")
  conditions <- c("Wilt Type",
                  expression(Delta*pgpB),
                  expression(Delta*cdh),
                  expression(Delta*cfa),
                  expression(Delta*plsX),
                  expression(Delta*opgB),
                  expression(Delta*clsC),
                  expression(Delta*pgpC),
                  expression(Delta*pldC),
                  expression(Delta*clsB),
                  expression(Delta*clsA),
                  expression(Delta*pldB),
                  expression(Delta*pldA),
                  expression(Delta*pgpA),
                  expression(Delta*aas),
                  expression(Delta*dgkA)
                  
  )
  output<- vector(mode = "list",6)
  names(output)<-c("PCA","Variance","PC1 loadings","PC2 loadings","loadings","data")
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }

  d<-data[,s]
  d<-data.frame(d)
  colnames(d)<-samples
  d<-t(d)
  colnames(d)<- data$`Metabolite name`
  
  coul<-c("azure4","blue4","red1","cyan3","chartreuse2","deeppink","cyan4",
          "chartreuse4","orange4","darkorchid3","brown4",
          "darkorange","darkmagenta","yellow3","deeppink4","black")
  
  color<-rep(coul,each=3)
  
  
  d[d<0]<-0
  index <- which(colSums(d)==0)
  if(length(index)!=0){
    d <- d[,-index]
  }
  
  d <- log10(d+1)
  
  d_p <-d
  for(p in 1:ncol(d)){
    d_p[,p] <- d[,p]/sqrt(sd(d[,p]))
  }
  
  PC<-prcomp(d_p,scale = F,center = T)
  
  layout(matrix(c(1,2),ncol=2),widths = c(8,2))
  par(mar=c(5,5,5,0))
  plot(PC$x[,1:2],col=color,pch=rep(c(16,17),each = 48),cex=2, 
       main=paste0("Principle Component Analysis: ",title),cex.lab=1.5,cex.main=1.5,cex.axis=1.2)   
  par(mar=c(5,0,5,0))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  
  color_matrix <- data.frame(Sample=unique(gsub("EP_","",samples)),color=coul)
  color_matrix <- color_matrix[c(16,c(1:15)),]
  
  legend("topleft",legend = conditions ,col =color_matrix$color ,cex=1.5,pch=15,bty="n" ,pt.cex = 2)
  p <- recordPlot()
  
  output[[1]]<- p
  plot.new() 
  
  
  VAR<-summary(PC)
  VAR<-VAR$importance[2,1:10]*100
  par(mfrow=c(1,1),mar=c(5,5,5,5))
  barplot(VAR,col="tomato3",ylab="Variance Explained [%]", main=paste0("Variance explained: ",title),cex.lab=1.5,cex.main=1.5,cex.axis=1.2,cex.names = 1.2)
  
  v<-recordPlot()
  output[[2]] <- v
  plot.new()
  
  loadings1<-PC$rotation
  d<-sort(abs(loadings1[,1]),decreasing = T)[1:100]
  D<-loadings1[order(abs(loadings1[,1]),decreasing = T),1]
  if(length(D)>20){
    D<-D[1:20]
  }
  
  par(mfrow=c(1,1),mar=c(14,6,3,5))
  barplot(D,ylab = "Eigenvalues / Weights",
          cex.axis = 1.2,cex.lab=1.5,cex.names = 1.1,main=paste0("Loadings PC 1: ",title),las=2,col="tomato3" ,cex.main=1.5)
  
  p <- recordPlot()
  
  output[[3]] <- p
  plot.new() 
  
  
  loadings2<-PC$rotation
  d<-sort(abs(loadings2[,2]),decreasing = T)[1:100]
  D<-loadings2[order(abs(loadings2[,2]),decreasing = T),2]
  if(length(D)>20){
    D<-D[1:20]
  }
  
  par(mfrow=c(1,1),mar=c(14,6,3,5))
  barplot(D,ylab = "Eigenvalues / Weights",
          cex.axis = 1.2,cex.lab=1.5,cex.names = 1.1,main=paste0("Loadings PC 2: ",title),las=2,col="tomato3",cex.main=1.5 )
  
  p <- recordPlot()
  
  output[[4]] <- p
  plot.new() 
  
  x<-PC$rotation[,1]
  y<-PC$rotation[,2]
  
  lab<- gsub("\\|.*", "",colnames(d_p) )
  class<-data$Ontology
  class_unique<-unique(class)
  
  col_load <- coul[1:length(class_unique)]
  clrs <- col_load[factor(class)]
  
  par(mfrow=c(1,1),mar=c(5,5,5,5))
  plot(x,y,pch=16,xlab="PC1",ylab="PC2",main=paste0("PCA Loadings Plot: ",title),cex=1.4,
       col=clrs,ylim=c(min(y)-0.05,max(y)+0.05),xlim=c(min(x)-0.01,max(x)+0.05),cex.main=1.5,cex.lab = 1.5)
  load<-data.frame(x=x,y=y,label=lab)
  
  load1<-load[which(abs(load$x)>0.1),]
  text(load1$x,load1$y,load1$label, pos=4)
  
  load<-load[which(abs(load$x)<=0.1),]
  load2<-load[which(abs(load$y)>0.1),]
  text(load2$x,load2$y,load2$label, pos=2)
  
  v<-recordPlot()
  output[[5]] <- v
  plot.new()
  
  output[[6]]<- PC
  
  return(output)
}

################################################################################

# function 1: PLS-DA Analysis
PLSDA<- function(data,phase,fragmentation){
  title <- paste0(fragmentation," data (",phase,")")
  conditions <- c("Wilt Type",
                  expression(Delta*pgpB),
                  expression(Delta*cdh),
                  expression(Delta*cfa),
                  expression(Delta*plsX),
                  expression(Delta*opgB),
                  expression(Delta*clsC),
                  expression(Delta*pgpC),
                  expression(Delta*pldC),
                  expression(Delta*clsB),
                  expression(Delta*clsA),
                  expression(Delta*pldB),
                  expression(Delta*pldA),
                  expression(Delta*pgpA),
                  expression(Delta*aas),
                  expression(Delta*dgkA)
                  
  )
  coul<-c("azure4","blue4","red1","cyan3","chartreuse2","deeppink","cyan4",
          "chartreuse4","orange4","darkorchid3","brown4",
          "darkorange","darkmagenta","yellow3","deeppink4","black")
  
  output<- vector(mode = "list",4)
  names(output)<-c("PLSDA","Variance","VIP C1","VIP C2")
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  
  d<-data[,s]
  
  d<-t(d)

  X <- d
  colnames(X)<- paste0(c(1:ncol(X)),"|",data$`Metabolite name`)
  
  Y <-  as.factor(samples)

  samples <- paste0(samples,".",c(1:3))
  rownames(X)<-samples  
  
  color<-rep(coul,each=3)
  
  
  X[X<0]<-0
  index <- which(colSums(X)==0)
  if(length(index)!=0){
    X <- X[,-index]
  }
  
  X <- log10(X+1)
  
  X_p <-X
  for(p in 1:ncol(X)){
    X_p[,p] <- X[,p]/sqrt(sd(X[,p]))
  }
  
  plsda.d <- plsda(X,Y, ncomp = 10)
  
  
  
  VAR<-explained_variance(plsda.d$X, variates=plsda.d$variates$X, ncomp=10)
  
  
  par(mfrow=c(1,1),mar=c(5,5,5,5))
  barplot(VAR*100,col="tomato3",ylab="Variance Explained [%]", main=paste0("Variance explained: ",title),cex.lab=1.5,cex.main=1.5,cex.axis=1.2,cex.names = 1.2)
  
  v<-recordPlot()
  output[[2]] <- v
  plot.new()
  
  final.plsda.d <- plsda(X,Y, ncomp = 2)
  
  z<-plotIndiv(final.plsda.d, ind.names = FALSE, legend=TRUE,
               comp=c(1,2), ellipse = TRUE, 
               title = 'PLS-DA on SRBCT comp 1-2',
               X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2')
  
  
  layout(matrix(c(1,2),ncol=2),widths = c(8,2))
  par(mar=c(5,5,5,0))
  plot(z$df$x,z$df$y,col=color,pch=rep(c(16,17),each = 48),cex=2, 
       main=paste0("Partial Least Squares Analysis: ",title),cex.lab=1.5,cex.main=1.5,cex.axis=1.2,
       xlab="Component 1",ylab="Component 2")   
  par(mar=c(5,0,5,0))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  
  color_matrix <- data.frame(Sample=unique(gsub("EP_","",samples)),color=coul)
  color_matrix <- color_matrix[c(16,c(1:15)),]
  
  legend("topleft",legend = conditions ,col =color_matrix$color ,cex=1.5,pch=15,bty="n" ,pt.cex = 2)
  
  p <- recordPlot()
  
  output[[1]]<- p
  plot.new() 
  par(mar=c(5,14,5,5),mfrow=c(1,1))
  
  q<-vip(final.plsda.d)
  q<-q[order(q[,1]),]
  i<-c((nrow(q)-39):nrow(q))
  q1 <- q[i,]
  
  plot(q1[,1],1:nrow(q1),xlab="VIP scores",ylab="",yaxt="n",
       main=paste0("Variable Importance Plot Component 1 (",fragmentation,"/",phase,")"),
       cex.lab=1.5,cex.main=1.5,cex.axis=1.2,pch=16,cex=2)
  
  l<-gsub(".*\\|", "",rownames(q1) )
  
  axis(2,at=(1:nrow(q1)),labels = l,las=2,cex.axis=1)
  abline(h=1:40,lwd=0.1,lty=2)
  
  p <- recordPlot()
  
  output[[3]] <- p
  plot.new() 
  
  
  
  q<-vip(final.plsda.d)
  q<-q[order(q[,2]),]
  i<-c((nrow(q)-39):nrow(q))
  q2 <- q[i,]
  
  plot(q2[,2],1:nrow(q2),xlab="VIP scores",ylab="",yaxt="n",
       main=paste0("Variable Importance Plot Component 2 (",fragmentation,"/",phase,")"),
       cex.lab=1.5,cex.main=1.5,cex.axis=1.2,pch=16,cex=2)
  
  l<-gsub(".*\\|", "",rownames(q2) )
  
  axis(2,at=(1:nrow(q2)),labels = l,las=2,cex.axis=1)
  abline(h=1:40,lwd=0.1,lty=2)
  p <- recordPlot()
  
  output[[4]] <- p
  plot.new() 
  

  return(output)
}

################################################################################
# function 2: Volcano Plots
Volcano_plot<-function(data,phase,fragmentation){
  title <- paste0("Volcano Plots WT vs. KO (",fragmentation,"-",phase,")")
  conditions <- c("Wilt Type",
                  expression(Delta*pgpB),
                  expression(Delta*cdh),
                  expression(Delta*cfa),
                  expression(Delta*plsX),
                  expression(Delta*opgB),
                  expression(Delta*clsC),
                  expression(Delta*pgpC),
                  expression(Delta*pldC),
                  expression(Delta*clsB),
                  expression(Delta*clsA),
                  expression(Delta*pldB),
                  expression(Delta*pldA),
                  expression(Delta*pgpA),
                  expression(Delta*aas),
                  expression(Delta*dgkA)
                  
  )
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  d<-data[,s]
  lipids<-data$`Metabolite name`
  d[d<0]<-0
  sequence<-seq(1,45,3)

  df<-vector(mode="list",15)
  test<-matrix(ncol=4,nrow=8)
  test[,1]<-2
  test[1,]<-1
  test[7,]<-18
  test[8,]<-19
  test[2:6,2:4]<-matrix(c(3:17),ncol=3)
  layout(test,widths=c(1,10,10,10),heights=c(1,8,8,8,8,8,1,2))
  
  par(mar=c(0,2,0,4))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center",legend=title,cex=2.5,bty="n")
  par(mar=c(2,0,1,0))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  mtext(2,text = expression(-log[10](p-value)),line=-2.5,cex=1.3)
  
  for(i in 1:15){
    name<-conditions[i+1]
    wt<-c(46,47,48)
    k<-sequence[i]
    ko<-c(k,k+1,k+2)
    dataframe <- d[,c(wt,ko)]
    dataframe1<-matrix(as.numeric(as.matrix(dataframe)),ncol=6)
    colnames(dataframe1)<-colnames(dataframe)
    dataframe<-dataframe1
    colSums(dataframe)
    
    w<-rowMeans(dataframe[,c(1:3)])
    k<-rowMeans(dataframe[,c(4:6)])
    
    FC<-k/w
    
    FC_log<-log2(FC)
    
    FC_log[FC_log==-Inf]<- min(FC_log)-1
    FC_log[FC_log==Inf]<- max(FC_log)+1
    
    pvalue <- rep(NA,nrow(dataframe))  
    
    for(j in 1:nrow(dataframe)){
      ww<-dataframe[j,c(1:3)]
      kk<-dataframe[j,c(4:6)]
      
      t<-t.test(ww,kk,alternative = "two.sided",conf.level = 0.95)
      value<- -log10(t$p.value)
      
      pvalue[j]<-value
    }
    
    col<-rep("grey",nrow(dataframe))
    m<-cbind(FC_log,pvalue)
    colnames(m)<-c("FC","pvalue")
    
    p<--log10(0.05)
    f_neg<-log2(0.5)
    f_pos<-log2(2)
    
    index<-which(m[,1]>=f_pos & m[,2]>=p)
    col[index] <-"palegreen3"
    
    index<-which(m[,1]<=f_neg & m[,2]>=p)
    col[index] <-"coral1"
    
    par(mar=c(2,2,1,2),mgp=c(2,0.4,0),tcl=-0.2)
    plot(FC_log,pvalue,pch=16,cex=1.5,ylab="",xlab="",col=col,cex.axis=1)
    abline(h=p,col="red",lty=c("19"),lwd=1)
    abline(v=f_neg,col="red",lty=c("19"),lwd=1)
    abline(v=f_pos,col="red",lty=c("19"),lwd=1)
    points(FC_log,pvalue,pch=16,cex=1.5,col=col)
    
    legend("topright",legend = name,bty="n",cex=1.2, adj = 0.5 )
    
    m<-data.frame(m)
    m$lipid<-lipids 
    se<-which(abs(m[,1])>=1 & m[,2]>p)
    
    df[[i]] <- se
  }
  par(mar=c(0,2,0,4))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  mtext(3,text = expression(log[2](FC)),line=-2.5,cex=1.3)
  par(mar=c(0,0,0,0))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center",legend = c("WT s. enriched","KO s. enriched"),col=c("coral1","palegreen3"),horiz = T,bty="n",cex=2,pch=16,pt.cex = 4,
         text.width = 0.2)
  
  p <- recordPlot()
  
  output<-vector(mode="list",2)
  names(output) <- c("Volcano plots","list of enriched lipids")
  
  output[[1]]<-p
  output[[2]]<-df
  return(output)
  plot.new()
  
  
  
}


################################################################################
# function 2: Sample Clustering

sample_cluster<-function(data,phase,fragmentation,mode,log=F,scale=F){
  title <- paste0("Sample clustering: ",phase," (",fragmentation,")")
  conditions <- c(
                  expression(Delta*pgpB),
                  expression(Delta*cdh),
                  expression(Delta*cfa),
                  expression(Delta*plsX),
                  expression(Delta*opgB),
                  expression(Delta*clsC),
                  expression(Delta*pgpC),
                  expression(Delta*pldC),
                  expression(Delta*clsB),
                  expression(Delta*clsA),
                  expression(Delta*pldB),
                  expression(Delta*pldA),
                  expression(Delta*pgpA),
                  expression(Delta*aas),
                  expression(Delta*dgkA,
                             "Wilt Type")
                  
  )
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP ",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP ",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  
  
  d<-data[,s]
  lipids<-data$`Metabolite name`
  d[d<0]<-0
  
  
  if(log==TRUE){
    d<-log10(d+1)
  }
  
  # Count matrix
  
  x<-samples
  x<-paste0(x," (",c(1:3),")")
  colnames(d) <- x
  rownames(d) <- paste0(1:length(lipids) ," ",lipids)
  
  
  nrow(d)
  
  keep <- rowSums(d,na.rm=T) > 0
  d <- d[keep,]
  lipids<-lipids[keep]
  nrow(d)
  
  if(scale==T){
    d<-scaling(d)
  }
  
  #poison distance
  if(mode=="Poiss"){
  poisd <- PoissonDistance(t((d)))

  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <- rep(conditions,each=3)
  colnames(samplePoisDistMatrix) <- NULL
  pheatmap(samplePoisDistMatrix,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,labels_row = rep(conditions,each=3)
           ,main=title,fontsize = 15,
           fontsize_row = 10
  )
  
  }else if(mode=="Euclidean"){
  
  # euclidean distance
  
  sampleDists <- dist(t((d)))
  
  #sampleDists <- dist(t(assay(vsd)))
  
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- colnames(d)
  colnames(sampleDistMatrix) <- NULL
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           labels_row = rep(conditions,each=3)
           ,main=title,fontsize = 15,
           fontsize_row = 10
  )
  
  }
  
  p<- recordPlot()
  
  plot.new()
  
  return(p)
}

################################################################################

# function 3: Lipid Clustering


lipid_cluster<-function(data,phase,fragmentation,mode,log=F,topN=F,scale=F){
  
  
  title <- paste0("Lipid Species clustering: ",phase," (",fragmentation,")")
  
   lipid<-data$`Metabolite name`
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP ",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP ",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  
  
  d<-data[,s]
  lipids<-data$`Metabolite name`
  d[d<0]<-0
  
  
  if(log==TRUE){
    d<-log10(d+1)
  }
  
  if(scale==T){
    d<-scaling(d)
  }
  # Count matrix
  
  x<-samples
  x<-paste0(x," (",c(1:3),")")
  colnames(d) <- x
  rownames(d) <- paste0(1:length(lipids) ," ",lipids)
  
  
  nrow(d)
  
  keep <- rowSums(d,na.rm=T) > 0
  d <- d[keep,]
lipids<-lipids[keep]
  nrow(d)
  
  size=5
  if(is.numeric(topN)==T){
    variation <- rowVars(as.matrix(d),na.rm=T)
    variation <- cbind(1:length(variation),variation)
    
    variation <- variation[order(variation[,2],decreasing = T),]
    
    order <- variation[,1]
    
    d<-d[order,]
    
    lipids<-lipids[order]
    
    index<-c(1:topN)
    
    d<-d[index,]
    lipids<-lipids[index]
    size=9
  }
  
  
  
  #poison distance
  if(mode=="Poiss"){
    poisd <- PoissonDistance(d)
    
    samplePoisDistMatrix <- as.matrix( poisd$dd )
    rownames(samplePoisDistMatrix) <- lipids
    colnames(samplePoisDistMatrix) <- NULL
    
    pheatmap(samplePoisDistMatrix,
             clustering_distance_rows = poisd$dd,
             clustering_distance_cols = poisd$dd,labels_row = lipids,
             fontsize_row = size,
             fontsize=15,main = title
    )
  }else if(mode=="Euclidean"){
    
    # euclidean distance
    
    sampleDists <- dist(((d)))
    

    sampleDistMatrix <- as.matrix( sampleDists )
    rownames(sampleDistMatrix) <- lipids
    colnames(sampleDistMatrix) <- NULL
    pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             labels_row = lipids,
             fontsize_row = size,
             fontsize = 15,main=title
    )
    
  }
  
  p<- recordPlot()
  
  plot.new()
  
  return(p)
}

################################################################################

# function 4: 


Bi_cluster<-function(data,phase,fragmentation,mode,log=F,topN=F,scale=F){
  title <- paste0("Lipid-Sample Bi-clustering: ",phase," (",fragmentation,")")
  
  conditions <- c(
    expression(Delta*pgpB),
    expression(Delta*cdh),
    expression(Delta*cfa),
    expression(Delta*plsX),
    expression(Delta*opgB),
    expression(Delta*clsC),
    expression(Delta*pgpC),
    expression(Delta*pldC),
    expression(Delta*clsB),
    expression(Delta*clsA),
    expression(Delta*pldB),
    expression(Delta*pldA),
    expression(Delta*pgpA),
    expression(Delta*aas),
    expression(Delta*dgkA,
               "Wilt Type")
    
  )
  lipids<-data$`Metabolite name`
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP ",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP ",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  
  
  d<-data[,s]
  lipids<-data$`Metabolite name`
  d[d<0]<-0
  
  
  if(log==TRUE){
    d<-log10(d+1)
  }
  
  if(scale==T){
    d<-scaling(d)
  }
  
  
  # Count matrix
  
  x<-samples
  x<-paste0(x," (",c(1:3),")")
  colnames(d) <- x
  rownames(d) <- paste0(1:length(lipids) ," ",lipids)
  
  
  nrow(d)
  
  keep <- rowSums(d,na.rm=T) > 0
  d <- d[keep,]
  nrow(d)
  
  size=5
  if(is.numeric(topN)==T){
    variation <- rowVars(as.matrix(d),na.rm=T)
    variation <- cbind(1:length(variation),variation)
    
    variation <- variation[order(variation[,2],decreasing = T),]
    
    order <- variation[,1]
    
    d<-d[order,]
    
    lipids<-lipids[order]
    
    index<-c(1:topN)
    
    d<-d[index,]
    lipids<-lipids[index]
    size=9
  }
  
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  
  
  
  #poison distance
  if(mode=="Poiss"){
    mat_cluster_cols <- sort_hclust(hclust(PoissonDistance(t(d))$dd))
    
    mat_cluster_rows <- sort_hclust(hclust(PoissonDistance(d)$dd))
    
    pheatmap(d,
                cluster_cols      = mat_cluster_cols,
                cluster_rows      = mat_cluster_rows,
                fontsize_col = 10,
                fontsize = 15,
                fontsize_row = size,
                labels_row = lipids,
                labels_col = rep(conditions,each=3),angle_col="90",
             main = title
             
    )
    
    
  }else if(mode=="Euclidean"){
    
    # euclidean distance
    
    sampleDists <- dist(((d)))
    
    mat_cluster_cols <- sort_hclust(hclust(dist(t(d))))
    
    mat_cluster_rows <- sort_hclust(hclust(dist(d)))
    
    pheatmap(d,
             cluster_cols      = mat_cluster_cols,
             cluster_rows      = mat_cluster_rows,
             fontsize = 15,
             fontsize_col = 10,
             main=title,
             fontsize_row = size,
             labels_row = lipids,
             labels_col = rep(conditions,each=3),angle_col="90"
             
    )
    
  }
  
  p<- recordPlot()
  
  plot.new()
  
  return(p)
}

Bi_cluster<-function(data,phase,fragmentation,mode,log=F,topN=F,scale=F){
  title <- paste0("Lipid-Sample Bi-clustering: ",phase," (",fragmentation,")")
  conditions <- c(
    expression(Delta*pgpB),
    expression(Delta*cdh),
    expression(Delta*cfa),
    expression(Delta*plsX),
    expression(Delta*opgB),
    expression(Delta*clsC),
    expression(Delta*pgpC),
    expression(Delta*pldC),
    expression(Delta*clsB),
    expression(Delta*clsA),
    expression(Delta*pldB),
    expression(Delta*pldA),
    expression(Delta*pgpA),
    expression(Delta*aas),
    expression(Delta*dgkA,
               "Wilt Type")
    
  )
  lipids<-data$`Metabolite name`
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP ",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP ",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  
  
  d<-data[,s]
  lipids<-data$`Metabolite name`
  d[d<0]<-0
  
  
  if(log==TRUE){
    d<-log10(d+1)
  }
  
  if(scale==T){
    d<-scaling(d)
    
  }
  
  
  # Count matrix
  lipids[which(lipids=="DG 25:7(4,7,10,13,16,19,22)/27:7(6,9,12,15,18,21,24)")]<-"DG 25:7/27:7"
  x<-samples
  x<-paste0(x," (",c(1:3),")")
  colnames(d) <- x
  rownames(d) <- paste0(1:length(lipids) ," ",lipids)
  
  
  nrow(d)
  
  keep <- rowSums(d,na.rm=T) > 0
  d <- d[keep,]
  nrow(d)
  
  size=5
  if(is.numeric(topN)==T){
    
    wt<-d[,46:48]
    dgka<-d[,43:45]
    
    
    variation<-rowMeans(dgka,na.rm = T)/rowMeans(wt,na.rm = T)
    variation[variation==Inf]<-0
    variation <- cbind(1:length(variation),variation)
    
    variation <- variation[order(variation[,2],decreasing = T),]
    
    order <- variation[,1]
    
    d<-d[order,]
    
    lipids<-lipids[order]
    
    index<-c(1:(topN-10))
    index_not1<-c(109,129,133,136,147)
    index_not2<-c(89,78,65,146,112)
    
    index<-c(index_not1,index,index_not2)
    d<-d[index,]
    lipids<-lipids[index]
    size=9
  }
  
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  
  
  
  #poison distance
  if(mode=="Poiss"){
    mat_cluster_cols <- sort_hclust(hclust(PoissonDistance(t(d))$dd))
    
    mat_cluster_rows <- sort_hclust(hclust(PoissonDistance(d)$dd))
    
    pheatmap(d,
             cluster_cols      = mat_cluster_cols,
             cluster_rows      = mat_cluster_rows,
             fontsize_col = 10,
             fontsize = 15,
             fontsize_row = size,
             labels_row = lipids,
             labels_col = rep(conditions,each=3),angle_col="90",
             main = title
             
             
    )
    
    
  }else if(mode=="Euclidean"){
    
    # euclidean distance
    
    sampleDists <- dist(((d)))
    
    mat_cluster_cols <- sort_hclust(hclust(dist(t(d))))
    
    mat_cluster_rows <- sort_hclust(hclust(dist(d)))
    
    pheatmap(d,
             cluster_cols      = mat_cluster_cols,
             cluster_rows      = mat_cluster_rows,
             fontsize = 15,
             fontsize_col = 10,
             main=title,
             fontsize_row = size,
             labels_row = lipids,
             labels_col = rep(conditions,each=3),angle_col="90"
             
    )
    
  }
  
  p<- recordPlot()
  
  plot.new()
  
  return(p)
}
###############################################################################
# FC 
# function 2: Volcano Plots
Fold_change<-function(data,phase,fragmentation,topN=F){
  title <- paste0("Fold Change: WT vs. KO (",fragmentation,"-",phase,")")
  conditions <- c(
    expression(Delta*pgpB),
    expression(Delta*cdh),
    expression(Delta*cfa),
    expression(Delta*plsX),
    expression(Delta*opgB),
    expression(Delta*clsC),
    expression(Delta*pgpC),
    expression(Delta*pldC),
    expression(Delta*clsB),
    expression(Delta*clsA),
    expression(Delta*pldB),
    expression(Delta*pldA),
    expression(Delta*pgpA),
    expression(Delta*aas),
    expression(Delta*dgkA)
    
  )
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  d<-data[,s]
  lipids<-data$`Metabolite name`
  d[d<0]<-0
  sequence<-seq(1,45,3)
  
  keep <- rowSums(d[,c(46:48)],na.rm=T) > 0
  d <- d[keep,]
  nrow(d)
  lipids<-lipids[keep]
  size=5
  if(is.numeric(topN)==T){
    variation <- rowVars(as.matrix(d),na.rm=T)
    variation <- cbind(1:length(variation),variation)
    
    variation <- variation[order(variation[,2],decreasing = T),]
    
    order <- variation[,1]
    
    d<-d[order,]
    
    lipids<-lipids[order]
    
    index<-c(1:topN)
    
    d<-d[index,]
    lipids<-lipids[index]
    size=9
  }
  
  
  
  
  
  FoldChange<-p_value<-matrix(ncol=15,nrow=nrow(d))
  
  i=1
  for(i in 1:15){
    wt<-c(46,47,48)
    k<-sequence[i]
    ko<-c(k,k+1,k+2)
    dataframe <- d[,c(wt,ko)]
    dataframe1<-matrix(as.numeric(as.matrix(dataframe)),ncol=6)
    colnames(dataframe1)<-colnames(dataframe)
    dataframe<-dataframe1
    
    w<-rowMeans(dataframe[,c(1:3)])
    k<-rowMeans(dataframe[,c(4:6)])
    
    FC<-k/w
    
    FC_log<-log2(FC)
    
    FC_log[FC_log==-Inf]<- min(FC_log)-1
    FC_log[FC_log==Inf]<- max(FC_log)+1
    
    pvalue <- rep(NA,nrow(dataframe))  
    
    for(j in 1:nrow(dataframe)){
      ww<-dataframe[j,c(1:3)]
      kk<-dataframe[j,c(4:6)]
      
      t<-t.test(ww,kk,alternative = "two.sided",conf.level = 0.95)
      value<- (t$p.value)
      if(value > 0.05){
        value<-""
      }else if(value <= 0.05 & value > 0.01){
        value<-"*"
      }else if(value <= 0.01 & value > 0.001){
        value <- value<-"**"
      }else if(value <= 0.001){
        value="***"
      }else{value<-""}
      pvalue[j]<-value
    }
    
   
    
    
    FoldChange[,i]<-FC_log
    p_value[,i]<-pvalue
    
    FoldChange[FoldChange=="NaN"] <- NA
    
    p_value[p_value=="NaN"] <- NA
    
    
  }
  
  FoldChange[FoldChange== -Inf] <- min(is.finite(unlist(FoldChange)))
  FoldChange[FoldChange== Inf] <-max(is.finite(unlist(FoldChange)))
  pheatmap(FoldChange,display_numbers = p_value,cluster_rows=FALSE, cluster_cols=FALSE,
           labels_row = lipids,fontsize = 15,main=title,
           labels_col = (conditions),fontsize_row = size,angle_col="90",fontsize_number = 7,na_col="black",border_color = NA)
  
  p<-recordPlot()
  plot.new()  
  
  return(p)
}

##################################################################################
FC_f<-function(data,topN=F,size,title,WT,log=T){
  conditions <- c(
    expression(Delta*pgpB),
    expression(Delta*cdh),
    expression(Delta*cfa),
    expression(Delta*plsX),
    expression(Delta*opgB),
    expression(Delta*clsC),
    expression(Delta*pgpC),
    expression(Delta*pldC),
    expression(Delta*clsB),
    expression(Delta*clsA),
    expression(Delta*pldB),
    expression(Delta*pldA),
    expression(Delta*pgpA),
    expression(Delta*aas),
    expression(Delta*dgkA)
    
  )
  
  if(log==T){
    lab<-"log(Concentration) [WT]"
    WT<-log2(WT)
  }else if(log==F){
    lab="Concentration [WT]"
  }
  
  d<-data
  lipids<-rownames(data)
  
  sequence<-seq(1,45,3)
  
  keep <- rowSums(d[,c(46:48)],na.rm=T) > 0
  d <- d[keep,]
  nrow(d)
  WT<-WT[keep]
  lipids<-lipids[keep]
  if(is.numeric(topN)==T){
    variation <- rowVars(as.matrix(d),na.rm=T)
    variation <- cbind(1:length(variation),variation)
    
    variation <- variation[order(variation[,2],decreasing = T),]
    
    order <- variation[,1]
    
    d<-d[order,]
    
    lipids<-lipids[order]
    
    index<-c(1:topN)
    
    d<-d[index,]
    lipids<-lipids[index]
  }
  
  
  
  
  
  FoldChange<-p_value<-matrix(ncol=15,nrow=nrow(d))
  
  for(i in 1:15){
    wt<-c(46,47,48)
    k<-sequence[i]
    ko<-c(k,k+1,k+2)
    dataframe <- d[,c(wt,ko)]
    dataframe1<-matrix(as.numeric(as.matrix(dataframe)),ncol=6)
    colnames(dataframe1)<-colnames(dataframe)
    dataframe<-dataframe1
    
    w<-rowMeans(dataframe[,c(1:3)])
    k<-rowMeans(dataframe[,c(4:6)])
    
    FC<-k/w
    
    FC_log<-log2(FC)
    
    FC_log[FC_log==-Inf]<- min(FC_log)-1
    FC_log[FC_log==Inf]<- max(FC_log)+1
    
    pvalue <- rep(NA,nrow(dataframe))  
    
    for(j in 1:nrow(dataframe)){
      ww<-dataframe[j,c(1:3)]
      kk<-dataframe[j,c(4:6)]
      
      t<-t.test(ww,kk,alternative = "two.sided",conf.level = 0.95)
      value<- (t$p.value)
      
      if(value > 0.05){
        value<-""
      }else if(value <= 0.05 & value > 0.01){
        value<-"*"
      }else if(value <= 0.01 & value > 0.001){
        value <- value<-"**"
      }else if(value <= 0.001){
        value="***"
      }

      
      pvalue[j]<-value
    }
    
    
    
    FoldChange[,i]<-FC_log
    p_value[,i]<-pvalue
    
    FoldChange[FoldChange=="NaN"] <- NA
    
    p_value[p_value=="NaN"] <- NA
    
    
  }
  
  FoldChange[FoldChange== -Inf] <- min(is.finite(unlist(FoldChange)))
  FoldChange[FoldChange== Inf] <-max(is.finite(unlist(FoldChange)))
  
  x<-pheatmap(FoldChange,display_numbers = p_value,cluster_rows=FALSE, cluster_cols=FALSE,
              labels_row = lipids,fontsize_col = 15,
              labels_col = conditions ,fontsize_row = size,main=title,angle_col="90",na_col="black" )#labels_col = (conditions)
  
  
  y<-ggplot(data.frame(WT),aes(seq_along(WT),WT,width=0.9))+geom_bar(stat="identity",fill="steelblue")+ 
    theme_minimal() +
    ylab(lab) + 
    xlab("Concentration Distribution in Wild Type Strain") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    coord_flip()+ scale_y_reverse()+scale_x_reverse(expand = c(0,0))+
    theme(axis.title.y=element_text(size=11),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x = element_text( size=11))+
          theme(plot.margin = margin(20, 0, 30, 5 ))
  
  grid.arrange(y,x[[4]], ncol=2, nrow=1, widths=c(1,4)) 
  
  p<-recordPlot()
  
  return(p)
  
  
  
  
  
  

}
##################################################################################
# Fold Change overview class



FC_class_bp <- function(data,phase,fragmentation){
  conditions <- c(
    expression(Delta*pgpB),
    expression(Delta*cdh),
    expression(Delta*cfa),
    expression(Delta*plsX),
    expression(Delta*opgB),
    expression(Delta*clsC),
    expression(Delta*pgpC),
    expression(Delta*pldC),
    expression(Delta*clsB),
    expression(Delta*clsA),
    expression(Delta*pldB),
    expression(Delta*pldA),
    expression(Delta*pgpA),
    expression(Delta*aas),
    expression(Delta*dgkA)
    
  )
  
  title=paste0("Fold Change per Lipid Class (KO/WT) ", "[",fragmentation,"-",phase,"]")
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  
  if(fragmentation=="CID"){
    color<-"coral1"
  }else if(fragmentation=="EAD"){
    color <- "palegreen3"
  }
  class <- data$Ontology
  d<-data[,s]
  d[d<0]<-0
  x<-unique(class)
  
  
  
  count<-matrix(ncol=ncol(d),nrow=length(x))
  
  rownames(count)<-x
  
  for(i in 1:length(x)){
    name<-x[i]
    index<-which(class==name)
    for(j in 1:48){
      count[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
    }
  }
  
  
  sequence<-seq(1,45,3)
  
  
  test<-matrix(ncol=4,nrow=7)
  test[,1]<-2
  test[1,]<-1
  test[7,]<-c(18,19,20,21)
  test[2:6,2:4]<-matrix(c(3:17),ncol=3)
  layout(test,widths=c(1,10,10,10),heights=c(2,8,8,8,8,8,7))
  
  par(mar=c(0,2,0,4))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center",legend=title,cex=2.5,bty="n")
  par(mar=c(2,0,1,0))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  mtext(2,text = expression(log[2](FC)),line=-2.5,cex=1.3)
  i=1
  for(i in 1:15){
    name<-conditions[i]
    wt<-c(46,47,48)
    k<-sequence[i]
    ko<-c(k,k+1,k+2)
    dataframe <- count[,c(wt,ko)]
    dataframe1<-matrix(as.numeric(as.matrix(dataframe)),ncol=6)
    colnames(dataframe1)<-colnames(dataframe)
    dataframe<-dataframe1
    colSums(dataframe)
    
    w<-rowMeans(dataframe[,c(1:3)])
    k<-rowMeans(dataframe[,c(4:6)])
    
    null <- which(w+k==0)
    
    FC<-k/w
    
    FC_log<-log2(FC)
    FC_log[null]<- 0
    
    
    
    filter<-!is.infinite(FC_log)
    FC_log[FC_log==-Inf]<- min(FC_log[filter])-1
    FC_log[FC_log==Inf]<- max(FC_log[filter])+1
    
    pvalue <- rep(NA,nrow(dataframe))  
    j=3
    for(j in 1:nrow(dataframe)){
      ww<-dataframe[j,c(1:3)]
      kk<-dataframe[j,c(4:6)]
      
      t<-t.test(ww,kk,alternative = "two.sided",conf.level = 0.95)
      value<- (t$p.value)
      
      
      if(sum(ww+kk)==0){
        value<-0.0001
      }
      
      if(value > 0.05){
        value<-""
      }else if(value <= 0.05 & value > 0.01){
        value<-"*"
      }else if(value <= 0.01 & value > 0.001){
        value <- value<-"**"
      }else if(value <= 0.001){
        value="***"
      }else{value<-""}
      pvalue[j]<-value
    }
    
    par(mar=c(2,2,1,2),mgp=c(2,0.4,0),tcl=-0.2)
    q<-barplot(FC_log,ylim=c(-max(abs(FC_log)),max(abs(FC_log))),col=color)
    
    y<-FC_log
    y[y>=0]<-y[y>=0]+0.3
    y[y<0]<-y[y<0]-0.3
    text(q,y,pvalue)
    abline(h=1,col="red",lty=c("19"),lwd=1)
    abline(h=-1,col="red",lty=c("19"),lwd=1)
    mtext(3,text=name,adj=0,line=-0.4)

    
    
    
  }
  par(mar=c(0,0,0,0),mgp=c(2,0.4,0),tcl=-0.2)
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  
  
  par(mar=c(2,2,0,2),mgp=c(2,0.4,0),tcl=-0.2)
  q<-barplot(FC_log,ylim=c(0,10),xaxt="n",yaxt="n",bty="n",type="n",col = NA, border = NA)
  text(x=q,y=rep(10,length(q)),x,srt = 90,cex=1.2,pos=2)
  
  
  par(mar=c(2,2,0,2),mgp=c(2,0.4,0),tcl=-0.2)
  q<-barplot(FC_log,ylim=c(0,10),xaxt="n",yaxt="n",bty="n",type="n",col = NA, border = NA)
  text(x=q,y=rep(10,length(q)),x,srt = 90,cex=1.2,pos=2)
  
  par(mar=c(2,2,0,2),mgp=c(2,0.4,0),tcl=-0.2)
  q<-barplot(FC_log,ylim=c(0,10),xaxt="n",yaxt="n",bty="n",type="n",col = NA, border = NA)
  text(x=q,y=rep(10,length(q)),x,srt = 90,cex=1.2,pos=2)
  
  p <- recordPlot()
  plot.new()
  return(p)
  
  
  
  
  
}

##################################################################################

# Chain Length

chain_length<-function(data,phase,fragmentation,topN=F,class=F,size=9){
  title=paste0("Chain Length (KO/WT) FC ", "[",fragmentation,"-",phase,"]")
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  
  if(fragmentation=="CID"){
    ch<-characteristics_CID(data)
  }else if(fragmentation=="EAD"){
    ch<- characteristics_EAD(data)
  }
  
  d<-data[,s]
  d[d<0]<-0
  
  unique <- sort(unique(ch$`Chain length`))
  db<-matrix(ncol=48,nrow=length(unique))
  x<-rep(sort(paste0( "EP ",c("WT",1:15)) ),each=3)
  x<-paste0(x," (",c(1:3),")")
  colnames(db) <- x
  rownames(db) <- paste0("Chain length = ",unique)
  
  for(i in 1:length(unique)){
    n<-unique[i]
    index <- which(ch$`Chain length`==n)
    for(j in 1:48){
      db[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
    }
  }
  
  
  WT<- rowSums(db[,46:48])
  FC_f(db,topN = F,size=9,title,WT,log = F)
  
  
  if(class==T){
    title=paste0("Per Lipid Class: Chain Length (KO/WT) FC ", "[",fragmentation,"-",phase,"]")
    
    ch$feature <- paste0(ch$`Subclass level`, ": Chain length = ",ch$`Chain length` )
    
    unique <- sort(unique(ch$feature))
    
    db<-matrix(ncol=48,nrow=length(unique))
    x<-rep(sort(paste0( "EP ",c("WT",1:15)) ),each=3)
    x<-paste0(x," (",c(1:3),")")
    colnames(db) <- x
    rownames(db) <- unique
    i=1
    for(i in 1:length(unique)){
      n<-unique[i]
      index <- which( ch$feature == n)
      for(j in 1:48){
        db[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
      }
    }
    
    
    WT<- rowSums(db[,46:48])
    FC_f(db,topN = F,size=9,title,WT,log = F)    
  }
  
}


##################################################################################

# Double Bonds
db_number<-function(data,phase,fragmentation,topN=F,class=F,size=9){
  title=paste0("Number of DB (KO/WT) FC ", "[",fragmentation,"-",phase,"]")
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  
  if(fragmentation=="CID"){
    ch<-characteristics_CID(data)
  }else if(fragmentation=="EAD"){
    ch<- characteristics_EAD(data)
  }
  
  d<-data[,s]
  d[d<0]<-0
  
  unique <- sort(unique( ch$`# DB`))
  db<-matrix(ncol=48,nrow=length(unique))
  x<-rep(sort(paste0( "EP ",c("WT",1:15)) ),each=3)
  x<-paste0(x," (",c(1:3),")")
  colnames(db) <- x
  rownames(db) <- paste0("Number of DB = ",unique)
  
  for(i in 1:length(unique)){
    n<-unique[i]
    index <- which(ch$`# DB`==n)
    for(j in 1:48){
      db[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
    }
  }
  index<- grep("NA",rownames(db) )
  if(length(index)>0){
  db<-db[-index,]
  }
  
  WT<- rowSums(db[,46:48])
  FC_f(db,topN = F,size=9,title,WT,log = F)

  
  if(class==T){
    title=paste0("Per Lipid Class: Number of DB (KO/WT) FC ", "[",fragmentation,"-",phase,"]")
    
    ch$feature <- paste0(ch$`Subclass level`, ": # of DB = ",ch$`# DB` )
    
    unique <- sort(unique(ch$feature))
    
    db<-matrix(ncol=48,nrow=length(unique))
    x<-rep(sort(paste0( "EP ",c("WT",1:15)) ),each=3)
    x<-paste0(x," (",c(1:3),")")
    colnames(db) <- x
    rownames(db) <- unique
    i=1
    for(i in 1:length(unique)){
      n<-unique[i]
      index <- which( ch$feature == n)
      for(j in 1:48){
        db[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
      }
    }
    
    index<- grep("NA",rownames(db) )
    db<-db[-index,]
    
    WT<- rowSums(db[,46:48])
    FC_f(db,topN = F,size=9,title,WT,log = F)
    

  }
  
}



##################################################################################

# DB position

db_position<-function(data,phase,fragmentation,topN=F,class=F,size=9){
  title=paste0("DB position (KO/WT) FC ", "[",fragmentation,"-",phase,"]")
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  
  if(fragmentation=="CID"){
    print("No positional information for CID data")
  }else if(fragmentation=="EAD"){
    ch<- characteristics_EAD(data)
    
  
  
  d<-data[,s]
  d[d<0]<-0
  
  unique <- sort(unique( ch$`DB positions`))
  index<-which(unique=="")
  unique<-unique[-index]
  
  unique<-str_split(unique,",")
  unique<-unlist(unique)
  unique<-sort(as.numeric(unique(unique)))
  position <- matrix(ncol = length(unique),nrow=nrow(data))
  colnames(position)<-unique
  
  for(i in 1:nrow(d)){
    pos<- str_split(ch$`DB positions`[i],",")[[1]]
    test<-as.numeric(pos)
    if(all(is.na(test))==F){
    for(j in 1:length(pos)){
      p<-pos[j]
      index<-which(colnames(position)==p)
      position[i,index] <- 1
    }
    }  
  }
  
  
  
  
  
  db<-matrix(ncol=48,nrow=ncol(position))
  x<-rep(sort(paste0( "EP ",c("WT",1:15)) ),each=3)
  x<-paste0(x," (",c(1:3),")")
  colnames(db) <- x
  rownames(db) <- colnames(position)
  
  for(i in 1:length(colnames(position))){
    n<-colnames(position)[i]
    index <- which(position[,i] == 1)
    for(j in 1:48){
      db[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
    }
  }
  rownames(db) <- paste0("DB position = ",colnames(position))
  WT<- rowSums(db[,46:48])
  FC_f(db,topN = F,size=9,title,WT,log=F)
  
  
  
  
  
  
  
  if(class==T){
    title=paste0("Per Lipid Class: DB position (KO/WT) FC ", "[",fragmentation,"-",phase,"]")
    
    class <- unique(data$Ontology )
    name<-c()
    for(i in 1:nrow(d)){
      pos<- str_split(ch$`DB positions`[i],",")[[1]]
      test<-as.numeric(pos)
      if(all(is.na(test))==F){
        l<-data$Ontology[i]
        
        new <- paste0(l,":",pos)
        
        name <- c(name,new)
      }
    }
    
    unique <- sort(unique(name))
    
    p<-gsub(".*\\:", "",unique )
    
    c<-gsub("\\:.*", "",unique )
    
    name_s<-unique
    for(i in 1:length(unique(c))){
      n<-unique(c)[i]
      index<-which(c==n)
      
      name_s[index]<-paste0(n,":",sort(as.numeric(p[index])))
      
    }
    
    
    position <- matrix(ncol = length(name_s),nrow=nrow(data))
    colnames(position)<-name_s
    
    for(i in 1:nrow(d)){
      pos<- str_split(ch$`DB positions`[i],",")[[1]]
      l<-data$Ontology
      test<-as.numeric(pos)
      if(all(is.na(test))==F){
        pos<-paste0(l,":",pos)
        for(j in 1:length(pos)){
          p<-pos[j]
          index<-which(colnames(position)==p)
          position[i,index] <- 1
        }
      }  
    }
    
  
    db<-matrix(ncol=48,nrow=length(name_s))
    x<-rep(sort(paste0( "EP ",c("WT",1:15)) ),each=3)
    x<-paste0(x," (",c(1:3),")")
    colnames(db) <- x
    rownames(db) <- name_s

    
    
    for(i in 1:length(name_s)){
      n<-name_s[i]
      index <- which( position[,i] == 1)
      for(j in 1:48){
        db[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
      }
    }
    
    p<-gsub(".*\\:", "",name_s )
    
    c<-gsub("\\:.*", "",name_s )
    
    
    rownames(db) <- paste0(c,": DB position = ",p)
    WT<- rowSums(db[,46:48])
    FC_f(db,topN = F,size=9,title,WT,log=F)
    


  }
  
}
}
#################################################################################
# SN position
sn<-function(data,phase,fragmentation,topN=F,class=F,size=9){
  title=paste0("Acyl Chain species (KO/WT) FC ", "[",fragmentation,"-",phase,"]")
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  
  if(fragmentation=="CID"){
    ch<-characteristics_CID(data)
  }else if(fragmentation=="EAD"){
    ch<- characteristics_EAD(data)
  }
  
  d<-data[,s]
  d[d<0]<-0
  
  index1<-grep("_",ch$Chains)
  index2<-grep("/",ch$Chains)
  index<-unique(c(index1,index2))
  ind<-index
  
  unique <- sort(unique( ch$Chains[index]))
  ####hier i am
  
  
  db<-matrix(ncol=48,nrow=length(unique))
  x<-rep(sort(paste0( "EP ",c("WT",1:15)) ),each=3)
  x<-paste0(x," (",c(1:3),")")
  colnames(db) <- x
  rownames(db) <- paste0("Chain species = ",unique)
  
  for(i in 1:length(unique)){
    n<-unique[i]
    index <- which(ch$Chains == n)
    for(j in 1:48){
      db[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
    }
  }
  WT<- rowSums(db[,46:48])
  FC_f(db,topN = F,size,title,WT,log=F)
  
  
  if(class==T){
    title=paste0("Per Lipid Class: Acyl Chain species (KO/WT) FC ", "[",fragmentation,"-",phase,"]")
    
    ch$feature <- paste0(ch$`Subclass level`, ": Chain species = ",ch$Chains )
    
    
    unique <- sort(unique(ch$feature[ind]))
    
    db<-matrix(ncol=48,nrow=length(unique))
    x<-rep(sort(paste0( "EP ",c("WT",1:15)) ),each=3)
    x<-paste0(x," (",c(1:3),")")
    colnames(db) <- x
    rownames(db) <- unique
    i=1
    for(i in 1:length(unique)){
      n<-unique[i]
      index <- which( ch$feature == n)
      for(j in 1:48){
        db[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
      }
    }
    
    WT<- rowSums(db[,46:48])
    FC_f(db,topN = F,size,title,WT,log=F)
    

  }
  
}



#################################################################################
# Number of cyclopropane
cyclo_number<-function(data,phase,fragmentation,topN=F,class=F,size=9){
  title=paste0("Lipids containing Cyclopropanes (KO/WT) FC ", "[",fragmentation,"-",phase,"]")
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  
  if(fragmentation=="CID"){
    ch<-characteristics_CID(data)
  }else if(fragmentation=="EAD"){
    ch<- characteristics_EAD(data)
  }
  
  d<-data[,s]
  d[d<0]<-0
  

  
  unique <- sort(unique( ch$cyclopropane))
  ####hier i am
  
  
  db<-matrix(ncol=48,nrow=length(unique))
  x<-rep(sort(paste0( "EP ",c("WT",1:15)) ),each=3)
  x<-paste0(x," (",c(1:3),")")
  colnames(db) <- x
  rownames(db) <- paste0("Lipids containing Cyclopropane= ",unique)
  
  for(i in 1:length(unique)){
    n<-unique[i]
    index <- which(ch$cyclopropane == n)
    for(j in 1:48){
      db[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
    }
  }
  
  WT<- rowSums(db[,46:48])
  FC_f(db,topN = F,size,title,WT,log=F)

  
  if(class==T){
    title=paste0("Per Lipid Class: Lipids containing Cyclopropanes  (KO/WT) FC ", "[",fragmentation,"-",phase,"]")
    
    ch$feature <- paste0(ch$`Subclass level`, ": Lipids contaiting Cyclopropanes  = ",ch$cyclopropane )
    
    
    unique <- sort(unique(ch$feature))
    
    db<-matrix(ncol=48,nrow=length(unique))
    x<-rep(sort(paste0( "EP ",c("WT",1:15)) ),each=3)
    x<-paste0(x," (",c(1:3),")")
    colnames(db) <- x
    rownames(db) <- unique
    i=1
    for(i in 1:length(unique)){
      n<-unique[i]
      index <- which( ch$feature == n)
      for(j in 1:48){
        db[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
      }
    }
    
    
    WT<- rowSums(db[,46:48])
    FC_f(db,topN = F,size,title,WT,log=F)

  }
  
}

#################################################################################
# Effect of Sn position and cyclopropane
# Effect of Sn position and cyclopropane
cyclo_sn<-function(data,phase,fragmentation,topN=F,class=F,size=9,FC=T){
  if(FC==T){
    title=paste0("Cyclopropane contaiting Acyl Chain FC ", "[",fragmentation,"-",phase,"]")
    if(phase=="EP"){
      s<-c(c(36:80),c(84:86))
      name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
      samples <- name[order(name)]
    }else if(phase=="SP"){
      s<-c(c(87:131),c(135:137))
      name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
      samples <- name[order(name)]
    }
    
    if(fragmentation=="CID"){
      
      stop("No sn attachement information for CID data!")
      
      
    }else if(fragmentation=="EAD"){
      ch<- characteristics_EAD(data)
    }
    
    d<-data[,s]
    d[d<0]<-0
    
    ch$sn1<-NA
    ch$sn2<-NA
    
    index1<-grep("_",ch$Chains)
    index2<-grep("/",ch$Chains)
    ind<-unique(c(index1,index2))
    
    names<-ch$Chains[ind]
    
    ch$Chains<-NA
    
    ch$Chains[ind]<-names
    
    store<-matrix(ncol=3,nrow=nrow(d))
    
    i=1
    for(i in 1:nrow(d)){
      
      if(grepl("_",ch$Chains[i])==T ){
        value<-str_split(ch$Chains[i],"_")[[1]]
      }else if(grepl("/",ch$Chains[i])==T){
        value<-str_split(ch$Chains[i],"/")[[1]]
      }else{
        value<-c("","")
      }
      
      j<-c(1:length(value))
      store[i,j] <-value
      
    }
    
    ch$sn1 <- gsub("\\:.*", "",store[,1] )
    ch$sn2 <- gsub("\\:.*", "",store[,2] )
    
    
    t1<-as.numeric(ch$sn1) %% 2
    t1[t1==1]<-"Yes"
    t1[t1==0]<-"No"
    t2<-as.numeric(ch$sn2) %% 2
    t2[t2==1]<-"Yes"
    t2[t2==0]<-"No"
    
    ch$sn1 <- t1
    ch$sn2 <- t2
    
    unique <- c("Yes","No")
    
    d<-d[c(7:9,46:48)]
    db1<-db2<-matrix(ncol=6,nrow=length(unique))
    x<-rep(sort(paste0( "EP ",c("WT","KO")) ),each=3)
    x<-paste0(x," (",c(1:3),")")
    colnames(db1)<-colnames(db2) <- x
    rownames(db1) <- paste0("SN1 containing Cyclopropane= ",unique)
    rownames(db2) <- paste0("SN1 containing Cyclopropane= ",unique)
    for(i in 1:length(unique)){
      n<-unique[i]
      index <- which(ch$sn1 == n)
      for(j in 1:6){
        db1[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
      }
    }
    
    FC <-pval <- matrix(nrow=nrow(db1),ncol=2)
    colnames(FC)<-colnames(pval) <- c("sn1","sn2")
    for(i in 1:nrow(db1)){
      KO<-as.numeric(db1[i,c(1,2,3)])
      WT<-as.numeric(db1[i,c(4,5,6)])
      
      FC[i,1]<- mean(KO,na.rm=T)/mean(WT,na.rm=T)
      t<-t.test(KO,WT,alternative = "two.sided",conf.level = 0.95)
      value<- (t$p.value)
      pval[i,1]<-value
      
    }
    
    
    
    for(i in 1:length(unique)){
      n<-unique[i]
      index <- which(ch$sn2 == n)
      for(j in 1:6){
        db2[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
      }
    }
    
    
    for(i in 1:nrow(db2)){
      KO<-as.numeric(db2[i,c(1,2,3)])
      WT<-as.numeric(db2[i,c(4,5,6)])
      
      FC[i,2]<- mean(KO,na.rm=T)/mean(WT,na.rm=T)
      t<-t.test(KO,WT,alternative = "two.sided",conf.level = 0.95)
      value<- (t$p.value)
      pval[i,2]<-value
      
    }
    
    pval[pval<=0.001]<-"***"
    pval[pval<=0.01]<-"**"
    pval[pval>0.01]<-"*"
    
    if(class==F){
      par(mar=c(2,2,2,2))
      
      layout(matrix(c(1,2),ncol=2),widths=c(8,4))
      par(mar=c(5,5,5,0))
      x<-barplot(FC,beside = T,col = c("coral1","palegreen3"),main=title,ylab = "Fold Change (KO/WT)",cex.lab=1.5,cex.names = 1.6,cex.main=1.5)
      text(unlist(x),unlist(FC)-0.01,unlist(pval),pos=1,cex=2)
      par(mar=c(5,0,5,0))
      plot(1:10,1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      legend(x=1,y=9.4,legend=c("Yes","No"),cex=2,bty="n",col=c("coral1","palegreen3"),pch=15,pt.cex = 5)
      text(x=5,y=9.5,"Contains Cyclopropane:",cex=2)
    }
    
    if(class==T){
      title=paste0("Per Lipid Class: Cyclopropane contaiting Acyl Chain FC ", "[",fragmentation,"-",phase,"]")
      
      ch$feature_sn1 <- paste0(ch$`Subclass level`, ": Cyclopropane contaiting Acyl Chain  = ",ch$sn1 )
      ch$feature_sn2 <- paste0(ch$`Subclass level`, ": Cyclopropane contaiting Acyl Chain  = ",ch$sn2 )
      
      
      unique1 <- sort(unique(ch$feature_sn1))
      index<-!grepl("= NA",unique1)
      unique1<-unique1[which(index==T)]
      unique2 <- sort(unique(ch$feature_sn2))
      index<-!grepl("= NA",unique2)
      unique2<-unique2[which(index==T)]
      
      
      db1<-matrix(ncol=6,nrow=length(unique1))
      x<-rep(sort(paste0( "EP ",c("WT","KO")) ),each=3)
      x<-paste0(x," (",c(1:3),")")
      colnames(db1) <- x
      rownames(db1) <- unique1
      for(i in 1:length(unique1)){
        n<-unique1[i]
        index <- which( ch$feature_sn1 == n)
        for(j in 1:6){
          db1[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
        }
      }
      
      index<-which(rowSums(db1)==0)
      if(length(index)>0){
        db1<-db1[-index,]
      }
      
      
      
      
      db2<-matrix(ncol=6,nrow=length(unique2))
      x<-rep(sort(paste0( "EP ",c("WT","KO")) ),each=3)
      x<-paste0(x," (",c(1:3),")")
      colnames(db2) <- x
      rownames(db2) <- unique2
      for(i in 1:length(unique2)){
        n<-unique2[i]
        index <- which( ch$feature_sn2 == n)
        for(j in 1:6){
          db2[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
        }
      }
      
      index<-which(rowSums(db2[,4:6])==0)
      if(length(index)>0){
        db2<-db2[-index,]
      }
      
      name<-sort(unique(c(rownames(db1),rownames(db2))))
      
      FC <-pval <- matrix(nrow=length(name),ncol=2)
      colnames(FC)<-colnames(pval) <- c("sn1","sn2")
      for(i in 1:length(name)){
        n<-name[i]
        
        j<-which(rownames(db1)==n)
        if(length(j>0)){   
          KO<-as.numeric(db1[j,c(1,2,3)])
          WT<-as.numeric(db1[j,c(4,5,6)])
          
          FC[i,1]<- mean(KO,na.rm=T)/mean(WT,na.rm=T)
          t<-t.test(KO,WT,alternative = "two.sided",conf.level = 0.95)
          value<- (t$p.value)
          pval[i,1]<-value
        }
        
      }
      
      
      
      for(i in 1:length(name)){
        n<-name[i]
        j<-which(rownames(db2)==n)
        KO<-as.numeric(db2[j,c(1,2,3)])
        WT<-as.numeric(db2[j,c(4,5,6)])
        
        FC[i,2]<- mean(KO,na.rm=T)/mean(WT,na.rm=T)
        t<-t.test(KO,WT,alternative = "two.sided",conf.level = 0.95)
        value<- (t$p.value)
        pval[i,2]<-value
        
        
      }
      
      pval[pval<=0.001]<-"***"
      pval[pval<=0.01]<-"**"
      pval[pval>0.01]<-"*"
      
      
      max_tot<-max(max(FC[,1],na.rm = T),max(FC[,2],na.rm = T))
      label<-gsub("\\:.*", "",name )
      cyclo<-gsub(".*\\=", "",name )
      
      col=rep("coral1",length(name))
      col[which(cyclo==" No")]<-"palegreen3"
      label <- paste0(label,": ",cyclo)
      
      plot.new()
      par(mar=c(5,5,5,5))
      layout(matrix(c(1,1,1,2,3,4),ncol=3,byrow=T),widths=c(8,8,8),heights=c(1,9))
      par(mar=c(0,5,0,5))
      plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      mtext(3,text=title,cex=1.8,line=-4)
      par(mar=c(12,5,2,0))
      x<-barplot((FC[,1]),col = col,main="sn1 Acyl Chain",ylab = "Fold Change (KO/WT)",cex.lab=1.8,cex.names = 1.6,cex.main=1.5,
                 names.arg = rep(label,1),las=2,ylim = c(0,max_tot+1))
      text(unlist(x),unlist(FC[,1])-0.01,unlist(pval[,1]),pos=1,cex=2)
      par(mar=c(12,5,2,0))
      
      x<-barplot((FC[,2]),col = col,main="sn2 Acyl Chain",cex.lab=1.4,cex.names = 1.6,cex.main=1.5,
                 names.arg = rep(label,1),las=2,ylim = c(0,max_tot+1),yaxt="n",ylab="")
      text(unlist(x),unlist(FC[,2])-0.01,unlist(pval[,2]),pos=1,cex=2)
      
      par(mar=c(5,0,12,0))
      plot(1:10,1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      legend(x=1.5,y=9.4,legend=c("Yes","No"),cex=2,bty="n",col=c("coral1","palegreen3"),pch=15,pt.cex = 5)
      text(x=4,y=9.5,"Contains Cyclopropane:",cex=2)
      
    }
  }else if(FC==F){ ###############################################################################################
    
    
    
    
    
    
    title=paste0("Cyclopropane contaiting Acyl Chain Abundance ", "[",fragmentation,"-",phase,"]")
    if(phase=="EP"){
      s<-c(c(36:80),c(84:86))
      name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
      samples <- name[order(name)]
    }else if(phase=="SP"){
      s<-c(c(87:131),c(135:137))
      name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
      samples <- name[order(name)]
    }
    # Hier
    if(fragmentation=="CID"){
      stop("No sn attachement information for CID data!")
    }else if(fragmentation=="EAD"){
      ch<- characteristics_EAD(data)
    }
    
    d<-data[,s]
    d[d<0]<-0
    
    ch$sn1<-NA
    ch$sn2<-NA
    
    index1<-grep("_",ch$Chains)
    index2<-grep("/",ch$Chains)
    ind<-unique(c(index1,index2))
    
    names<-ch$Chains[ind]
    
    ch$Chains<-NA
    
    ch$Chains[ind]<-names
    
    store<-matrix(ncol=3,nrow=nrow(d))
    
    i=1
    for(i in 1:nrow(d)){
      
      if(grepl("_",ch$Chains[i])==T ){
        value<-str_split(ch$Chains[i],"_")[[1]]
      }else if(grepl("/",ch$Chains[i])==T){
        value<-str_split(ch$Chains[i],"/")[[1]]
      }else{
        value<-c("","")
      }
      
      j<-c(1:length(value))
      store[i,j] <-value
      
    }
    
    ch$sn1 <- gsub("\\:.*", "",store[,1] )
    ch$sn2 <- gsub("\\:.*", "",store[,2] )
    
    
    t1<-as.numeric(ch$sn1) %% 2
    t1[t1==1]<-"Yes"
    t1[t1==0]<-"No"
    t2<-as.numeric(ch$sn2) %% 2
    t2[t2==1]<-"Yes"
    t2[t2==0]<-"No"
    
    ch$sn1 <- t1
    ch$sn2 <- t2
    
    unique <- c("Yes","No")
    
    d<-d[c(7:9,46:48)]
    db1<-db2<-matrix(ncol=6,nrow=length(unique))
    x<-rep(sort(paste0( "EP ",c("WT","KO")) ),each=3)
    x<-paste0(x," (",c(1:3),")")
    colnames(db1)<-colnames(db2) <- x
    rownames(db1) <- paste0("SN1 containing Cyclopropane= ",unique)
    rownames(db2) <- paste0("SN1 containing Cyclopropane= ",unique)
    for(i in 1:length(unique)){
      n<-unique[i]
      index <- which(ch$sn1 == n)
      for(j in 1:6){
        db1[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
      }
    }
    
    
    for(i in 1:length(unique)){
      n<-unique[i]
      index <- which(ch$sn2 == n)
      for(j in 1:6){
        db2[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
      }
    }
    
    FC1 <-matrix(nrow=nrow(db1),ncol=2)
    pval1 <- matrix(nrow=nrow(db1),ncol=1)
    rownames(FC1)<-rownames(pval1) <- rownames(db1)
    for(i in 1:nrow(db1)){
      KO<-as.numeric(db1[i,c(1,2,3)])
      WT<-as.numeric(db1[i,c(4,5,6)])
      
      FC1[i,1]<- mean(KO,na.rm=T)
      FC1[i,2]<- mean(WT,na.rm=T)
      
      t<-t.test(KO,WT,alternative = "two.sided",conf.level = 0.95)
      value<- (t$p.value)
      pval1[i,1]<-value
      
    }
    
    max1<-max(unlist(FC1))
    
    for(i in 1:length(unique)){
      n<-unique[i]
      index <- which(ch$sn2 == n)
      for(j in 1:6){
        db2[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
      }
    }
    
    
    FC2 <-matrix(nrow=nrow(db2),ncol=2)
    pval2 <- matrix(nrow=nrow(db2),ncol=1)
    rownames(FC2)<-rownames(pval2) <- rownames(db2)
    for(i in 1:nrow(db2)){
      KO<-as.numeric(db2[i,c(1,2,3)])
      WT<-as.numeric(db2[i,c(4,5,6)])
      
      FC2[i,1]<- mean(KO,na.rm=T)
      FC2[i,2]<- mean(WT,na.rm=T)
      
      t<-t.test(KO,WT,alternative = "two.sided",conf.level = 0.95)
      value<- (t$p.value)
      pval2[i,1]<-value
      
    }
    
    p<-function(x,y,pval){
      yy<-max(y)+20
      segments(x[1],y[1],x[1],yy,lwd=3)
      segments(x[2],y[2],x[2],yy,lwd=3)
      segments(x[1],yy,x[2],yy,lwd=3)
      x1<-x[2]-( (x[2]-x[1])/2)
      y1<-yy
      text(x1,y1,pval,pos=3,cex=2)
    }
    
    max2<-max(unlist(FC2))
    
    max_t<-max(c(max1,max2))+40
    
    pval1[pval1<=0.001]<-"***"
    pval1[pval1<=0.01]<-"**"
    pval1[pval1>0.01]<-"*"
    
    pval2[pval2<=0.001]<-"***"
    pval2[pval2<=0.01]<-"**"
    pval2[pval2>0.01]<-"*"
    
    
    if(class==F){
      
      plot.new()
      par(mar=c(5,5,5,5))
      layout(matrix(c(1,1,1,2,3,4),ncol=3,byrow=T),widths=c(8,8,8),heights=c(1,9))
      par(mar=c(0,5,0,5))
      plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      mtext(3,text=title,cex=1.8,line=-4)
      par(mar=c(12,5,2,0))
      label=c("Yes","No")
      x<-barplot(t(FC1),beside=T,col = c("coral1","palegreen3"),main="sn1 Acyl Chain",ylab = "Abundance",cex.lab=1.8,cex.names = 1.6,cex.main=1.5,
                 names.arg = rep(label,1),las=2,ylim = c(0,max_t))
      p(x=x[,1],y=FC1[1,],pval=pval1[1,1])
      p(x=x[,2],y=FC1[2,],pval=pval1[2,1])
      par(mar=c(12,5,2,0))
      
      x<-barplot(t(FC2),beside=T,col = c("coral1","palegreen3"),main="sn2 Acyl Chain",cex.lab=1.4,cex.names = 1.6,cex.main=1.5,
                 names.arg = rep(label,1),las=2,ylim = c(0,max_t),yaxt="n",ylab="")
      p(x=x[,1],y=FC2[1,],pval=pval2[1,1])
      p(x=x[,2],y=FC2[2,],pval=pval2[2,1])
      par(mar=c(5,0,5,0))
      plot(1:10,1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      legend(x=1.5,y=9.9,legend=c(expression(Delta*cfa),"Wild Type"),cex=2,bty="n",col=c("coral1","palegreen3"),pch=15,pt.cex = 5)
      text(x=3,y=9.9,"Conditions:",cex=2)
      
    }
    
    
    
    
    
    if(class==T){
      title=paste0("Per Lipid Class: Cyclopropane contaiting Acyl Chain Abundance ", "[",fragmentation,"-",phase,"]")
      
      ch$feature_sn1 <- paste0(ch$`Subclass level`, ": Cyclopropane contaiting Acyl Chain  = ",ch$sn1 )
      ch$feature_sn2 <- paste0(ch$`Subclass level`, ": Cyclopropane contaiting Acyl Chain  = ",ch$sn2 )
      
      
      unique1 <- sort(unique(ch$feature_sn1))
      index<-!grepl("= NA",unique1)
      unique1<-unique1[which(index==T)]
      unique2 <- sort(unique(ch$feature_sn2))
      index<-!grepl("= NA",unique2)
      unique2<-unique2[which(index==T)]
      
      
      db1<-matrix(ncol=6,nrow=length(unique1))
      x<-rep(sort(paste0( "EP ",c("WT","KO")) ),each=3)
      x<-paste0(x," (",c(1:3),")")
      colnames(db1) <- x
      rownames(db1) <- unique1
      for(i in 1:length(unique1)){
        n<-unique1[i]
        index <- which( ch$feature_sn1 == n)
        for(j in 1:6){
          db1[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
        }
      }
      
      
      
      
      
      
      db2<-matrix(ncol=6,nrow=length(unique2))
      x<-rep(sort(paste0( "EP ",c("WT","KO")) ),each=3)
      x<-paste0(x," (",c(1:3),")")
      colnames(db2) <- x
      rownames(db2) <- unique2
      for(i in 1:length(unique2)){
        n<-unique2[i]
        index <- which( ch$feature_sn2 == n)
        for(j in 1:6){
          db2[i,j] <- sum(as.numeric(d[index,j]),na.rm=T)
        }
      }
      
      
      
      name<-sort(unique(c(rownames(db1),rownames(db2))))
      
      
      FC1 <-matrix(nrow=nrow(db1),ncol=2)
      pval1 <- matrix(nrow=nrow(db1),ncol=1)
      rownames(FC1)<-rownames(pval1) <- rownames(db1)
      for(i in 1:nrow(db1)){
        KO<-as.numeric(db1[i,c(1,2,3)])
        WT<-as.numeric(db1[i,c(4,5,6)])
        
        FC1[i,1]<- mean(KO,na.rm=T)
        FC1[i,2]<- mean(WT,na.rm=T)
        
        t<-t.test(KO,WT,alternative = "two.sided",conf.level = 0.95)
        value<- (t$p.value)
        pval1[i,1]<-value
        
      }
      
      max1<-max(unlist(FC1))
      
      
      
      
      FC2 <-matrix(nrow=nrow(db2),ncol=2)
      pval2 <- matrix(nrow=nrow(db2),ncol=1)
      rownames(FC2)<-rownames(pval2) <- rownames(db2)
      for(i in 1:nrow(db2)){
        KO<-as.numeric(db2[i,c(1,2,3)])
        WT<-as.numeric(db2[i,c(4,5,6)])
        
        FC2[i,1]<- mean(KO,na.rm=T)
        FC2[i,2]<- mean(WT,na.rm=T)
        
        t<-t.test(KO,WT,alternative = "two.sided",conf.level = 0.95)
        value<- (t$p.value)
        pval2[i,1]<-value
        
      }
      
      p<-function(x,y,pval){
        yy<-max(y)+10
        segments(x[1],y[1],x[1],yy,lwd=3)
        segments(x[2],y[2],x[2],yy,lwd=3)
        segments(x[1],yy,x[2],yy,lwd=3)
        x1<-x[2]-( (x[2]-x[1])/2)
        y1<-yy
        text(x1,y1,pval,pos=3,cex=2)
      }
      
      max2<-max(unlist(FC2))
      
      max_t<-max(c(max1,max2))+30
      
      pval1[pval1<=0.001]<-"***"
      pval1[pval1<=0.01]<-"**"
      pval1[pval1>0.01]<-"*"
      
      pval2[pval2<=0.001]<-"***"
      pval2[pval2<=0.01]<-"**"
      pval2[pval2>0.01]<-"*"
      
      
      
      label<-gsub("\\:.*", "",name )
      cyclo<-gsub(".*\\=", "",name )
      
      col=c("coral1","palegreen3")
      label <- paste0(label,": ",cyclo)
      plot.new()
      par(mar=c(5,5,5,5))
      layout(matrix(c(1,1,1,2,3,4),ncol=3,byrow=T),widths=c(8,8,4),heights=c(1,9))
      par(mar=c(0,5,0,5))
      plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      mtext(3,text=title,cex=1.8,line=-4)
      par(mar=c(12,5,2,0))
      x<-barplot(t(FC1),beside=T,col = col,main="sn1 Acyl Chain",ylab = "Abundance",cex.lab=1.8,cex.names = 1.6,cex.main=1.5,
                 names.arg = rep(label,1),las=2,ylim = c(0,max_t))
      for(i in 1:length(name)){
        p(x=x[,i],y=FC1[i,],pval=pval1[i,1])
        
      }
      par(mar=c(12,5,2,0))
      
      x<-barplot(t(FC2),beside=T,col = col,main="sn2 Acyl Chain",cex.lab=1.4,cex.names = 1.6,cex.main=1.5,
                 names.arg = rep(label,1),las=2,ylim = c(0,max_t),yaxt="n",ylab="")
      for(i in 1:length(name)){
        p(x=x[,i],y=FC2[i,],pval=pval2[i,1])
        
      }
      
      par(mar=c(5,0,12,0))
      plot(1:10,1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      legend(x=0.5,y=9.9,legend=c(expression(Delta*cfa),"Wild Type"),cex=2,bty="n",col=c("coral1","palegreen3"),pch=15,pt.cex = 5)
      text(x=2.7,y=10,"Conditions:",cex=2)
      
      
      
      
      
      
    }
    
    
  }
  p<-recordPlot()
  return(p)
  par(mar=c(6,6,6,6))
  plot.new()
  
}




################################################################################
#
# Abraham Moyal
# 16.02.23
# Data Analysis CID vs EAD

################################################################################
#
# Load libraries
library(vctrs)
library(ggplot2)
library(factoextra)
library(dplyr)
library(ggVennDiagram)
library(tidyverse)
library(rcdk)
library(cluster)
library(svDialogs)
library(plotly)
library(rsconnect)
library(stringi)
library(insight)
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library(Cairo)
#library(ggpubr)
library(xcms)
library(mzR)
library(MSnbase)
library(pracma)
################################################################################



feature_overlap<-function(CID_wd,CID_name,EAD_wd,EAD_name,mz_tol = 0.05,rt_tol = 0.1,directory_venn_diagram){
  
  
  Output<-vector(mode="list",2)
  names(Output)<-c("Barplot","MZ-RT plot")
  Directory <- CID_wd
  file <-CID_name
  
  dataframe <- data_import(Directory,file)
  data_CID <-dataframe$data
  
  Directory <- EAD_wd
  file <- EAD_name
  
  dataframe <- data_import(Directory,file)
  data_EAD <-dataframe$data
  
  par(cex.main=1.5,cex.lab=1.6,cex.axis=1.4,mar=c(5,6,5,5))
  x<-barplot(c(nrow(data_CID),nrow(data_EAD)),col=c("coral1","palegreen3"),names=c("CID","EAD"),cex.names = 2.5,
             main="Number of Total MS1 Features per Fragmentation Method",ylab="# Features")
  text(x,c(nrow(data_CID),nrow(data_EAD)),c(nrow(data_CID),nrow(data_EAD)),cex=1.8,pos=1)
  p<-recordPlot()
  Output[[1]] <- p
  plot.new()
  
  
  
  
  
  
  
  c1 <-col2rgb("coral1")
  c1 <- rgb(c1[1],c1[2],c1[3], max = 255, alpha = 125)
  
  c2 <-col2rgb("palegreen3")
  c2 <- rgb(c2[1],c2[2],c2[3], max = 255, alpha = 125)
  
  par(cex.main=1.5,cex.lab=1.5,cex.axis=1.4,mar=c(5,6,5,5))
  
  plot(data_CID$`Average Rt(min)`,data_CID$`Average Mz`,col=c1,xlab="RT",ylab="mz",pch=16,main="2D EAD/CID Feature Map ",cex=1.5)
  
  
  points(data_EAD$`Average Rt(min)`,data_EAD$`Average Mz`,col=c2,xlab="RT",ylab="mz",pch=16,cex=1.5)
  legend("topright",pch=16,pt.cex =2,legend = c("CID","EAD"),col=c(c1,c2),cex=2,bty="n")
  
  p<-recordPlot()
  Output[[2]] <- p
  plot.new()
  
  
  EAD<-data_EAD[,c(1,2,3,5,8)]
  
  master_peak_list <- data_CID[,c(1,2,3,5,8)]
  master_peak_list$CID<-T
  master_peak_list$EAD<-F
  mz_tol<-mz_tol
  rt_tol<-rt_tol #min
  i=5
  for(i in 1:nrow(EAD)){
    RT<- as.numeric(EAD[i,2]) - as.numeric(master_peak_list[,2])
    RT<- abs(RT) < rt_tol
    MZ<- as.numeric(EAD[i,3]) - as.numeric(master_peak_list[,3])
    MZ <- abs(MZ) < mz_tol
    Adduct <- EAD[i,4]==master_peak_list[,4]
    
    score<-data.frame(RT=RT,MZ=MZ,Adduct=Adduct)
    score[score==T]<-1
    
    index<-which(rowSums(score)==3)
    
    if(length(index)>1){
      RT<- as.numeric(EAD[i,2]) - as.numeric(master_peak_list[,2])
      minimum <- which(RT==min(RT[index]))
      index <- minimum
    }
    
    if(length(index)!=0){
      master_peak_list[index,7]<-T
    }else if(length(index)==0){
      add<-matrix(c(EAD[i,],FALSE,TRUE),nrow=1)
      colnames(add) <- colnames(master_peak_list)
      master_peak_list <- rbind(master_peak_list,add)
      master_peak_list <- data.frame(master_peak_list)
    }
    
  }
  
  master_peak_list$CID <- paste0(master_peak_list$CID, c(1:nrow(master_peak_list)))
  master_peak_list$EAD <- paste0(master_peak_list$EAD, c(1:nrow(master_peak_list)))
  
  
  C<-master_peak_list$CID
  E<-master_peak_list$EAD
  
  C<-C[grep("TRUE",C)]
  E<-E[grep("TRUE",E)]
  
  
  wd<-getwd()
  setwd(directory_venn_diagram)
  myCol <- c("coral1","palegreen3")
  venn.diagram(
    x=list(C,E),
    category.names=c("CID","EAD"),filename = "Feature_venn_diagram.png",
    output=T,
    main="Feature based Venn Diagram",
    main.cex = 2,
    wd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = 1.4,
    fontface = "bold",
    fontfamily = "sans",
    main.fontfamily =  'sans',
    # # Set names
    cat.cex = 1.8,
    cat.fontface = "bold",
    
    cat.fontfamily = "sans",
  ) 
  
  return(Output)
  setwd(wd)
}











lipid_overlap<-function(CID_wd,CID_name,EAD_wd,EAD_name,directory_venn_diagram){
  
  
  Output<-vector(mode="list",1)
  names(Output)<-c("Barplot")
  Directory <- CID_wd
  file <-CID_name
  
  dataframe <- data_import(Directory,file)
  data_CID <-dataframe$data
  
  Directory <- EAD_wd
  file <- EAD_name
  
  dataframe <- data_import(Directory,file)
  data_EAD <-dataframe$data
  
  
  data_CID<-filtering(data_CID)
  data_EAD<-filtering(data_EAD)
  
  
  CID<-characteristics_CID(data_CID)
  
  
  EAD<-characteristics_EAD(data_EAD)
  
  
  CID_unique <- unique(CID$`Species level`)
  
  EAD_unique <- unique(EAD$`Species level`)
  
  CID_unique<-CID_unique[!is.na(CID_unique)]
  EAD_unique<-EAD_unique[!is.na(EAD_unique)]
  
  par(cex.main=1.5,cex.lab=1.6,cex.axis=1.4,mar=c(5,6,5,5))
  
  d<-c(nrow(data_CID),nrow(data_EAD),length(CID_unique),length(EAD_unique))
  x<-barplot(d,col=c("coral1","palegreen3"),names=c("CID","EAD","CID","EAD"),cex.names = 2.5,
             main="Number of Lipid Species per Fragmentation Method",ylab="# Features",ylim = c(0,max(d)+20))
  text(x,d,d,cex=1.8,pos=1)
  
  v<-x[2]+((x[3]-x[2])/2)
  
  abline(v=v,lwd=2,lty=c("19"))
  x1<-x[1]+((x[2]-x[1])/2)
  x2<-x[3]+((x[4]-x[3])/2)
  text(x1,max(d)+10,"All Lipid Species",cex=1.5)
  text(x2,max(d)+10,"Unique Lipid Species",cex=1.5)
  
  
  p<-recordPlot()
  Output[[1]] <- p
  plot.new()
  
  
  C<-CID_unique[!is.na(CID_unique)]
  E<-EAD_unique[!is.na(EAD_unique)]
  
  
  
  C_store <- setdiff(C,E)
  E_store <- setdiff(E,C)
  UU<-matrix(ncol=2,nrow=max(c(length(C_store),length(E_store))))
  colnames(UU)<-c("CID","EAD")
  UU[1:length(C_store),1]<-C_store
  UU[1:length(E_store),2]<-E_store
  
  wd<-getwd()
  setwd(directory_venn_diagram)
  myCol <- c("coral1","palegreen3")
  venn.diagram(
    x=list(C,E),
    category.names=c("EAD","CID"),filename = "Lipid_venn_diagram.png",
    output=T,
    main="Lipid Species based Venn Diagram",
    main.cex = 2,
    wd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = 1.4,
    fontface = "bold",
    fontfamily = "sans",
    rotation.degree =180,
    # # Set names
    cat.cex = 1.8,
    cat.fontface = "bold",
    main.fontfamily =  'sans',
    cat.fontfamily = "sans",
    cat.dist = c(0.1, 0.08)
  ) 
  
  write.csv(UU,paste0(directory_venn_diagram,"/unique_lipids.csv"))
  setwd(wd)
  
  return(Output)
}


overlap_barplot<-function(CID_wd,CID_name,EAD_wd,EAD_name,percentage=T){
  Output<-vector(mode="list",1)
  names(Output)<-c("Barplot")
  Directory <- CID_wd
  file <-CID_name
  
  dataframe <- data_import(Directory,file)
  data_CID <-dataframe$data
  
  Directory <- EAD_wd
  file <- EAD_name
  
  dataframe <- data_import(Directory,file)
  data_EAD <-dataframe$data
  
  
  data_CID<-filtering(data_CID)
  data_EAD<-filtering(data_EAD)
  
  
  CID<-characteristics_CID(data_CID)
  
  
  EAD<-characteristics_EAD(data_EAD)
  
  
  CID_f <- CID %>% 
    distinct(`Species level`, .keep_all = TRUE)
  
  CID_f<-CID_f[!is.na(CID_f$`Species level`),]
  
  
  
  EAD_f <- EAD %>%
    distinct(`Species level`, .keep_all = TRUE)
  
  EAD_f<-EAD_f[!is.na(EAD_f$`Species level`),]
  
  
  E<-unique(EAD_f$`Subclass level`)  
  C<-unique(CID_f$`Subclass level`)
  
  Tot<-unique(E,C)
  
  
  count<-matrix(ncol=3,nrow=length(Tot))
  rownames(count)<-Tot
  colnames(count)<- c("CID","EID","Overlap")
  i=1
  
  if(percentage==T){
    
    
    for(i in 1:length(Tot)){
      name<-Tot[i]
      
      indexC <- which(CID_f$`Subclass level`==name)
      indexE <- which(EAD_f$`Subclass level`==name)  
      
      cc<-CID_f$`Species level`[indexC]
      ee<-EAD_f$`Species level`[indexE]
      
      overlap<-intersect(cc,ee)
      union<-union(cc,ee)
      cc_u<-setdiff(cc,ee)
      ee_u<-setdiff(ee,cc)
      
      
      count[i,1] <- length(cc_u)/length(union)*100
      count[i,2] <- length(ee_u)/length(union)*100
      count[i,3] <- length(overlap)/length(union)*100
    }
    
    layout(matrix(c(1,2),ncol=2),widths=c(8,2))
    par(cex.main=1.8,cex.lab=1.6,cex.axis=1.4,mar=c(5,6,5,0))
    
    barplot(t(count),col=c("coral1","palegreen3","grey"),ylab ="Percentage of Identified Lipid Species",main="Per Class Species Overlap")
    par(mar=c(5,0,5,0))
    plot(1:10,1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
    legend("topleft",legend=c("CID: unique","EAD unique","Overlap"),cex=1.7,bty="n",col=c("coral1","palegreen3","grey"),pch=15,pt.cex = 3.5,
           y.intersp=2)
    
  }
  
  
  if(percentage==F){
    
    
    for(i in 1:length(Tot)){
      name<-Tot[i]
      
      indexC <- which(CID_f$`Subclass level`==name)
      indexE <- which(EAD_f$`Subclass level`==name)  
      
      cc<-CID_f$`Species level`[indexC]
      ee<-EAD_f$`Species level`[indexE]
      
      overlap<-intersect(cc,ee)
      union<-union(cc,ee)
      cc_u<-setdiff(cc,ee)
      ee_u<-setdiff(ee,cc)
      
      
      count[i,1] <- length(cc_u)
      count[i,2] <- length(ee_u)
      count[i,3] <- length(overlap)
    }
    
    layout(matrix(c(1,2),ncol=2),widths=c(8,2))
    par(cex.main=1.8,cex.lab=1.6,cex.axis=1.4,mar=c(5,6,5,0))
    
    barplot(t(count),col=c("coral1","palegreen3","grey"),ylab ="Number of Identified Lipid Species",main="Per Class Species Overlap")
    par(mar=c(5,0,5,0))
    plot(1:10,1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
    legend("topleft",legend=c("CID: unique","EAD unique","Overlap"),cex=1.7,bty="n",col=c("coral1","palegreen3","grey"),pch=15,pt.cex = 3.5,
           y.intersp=2)
    
  }
  
  P<-recordPlot()
  Output[[1]]<-P
  plot.new()
  return(Output)
}




cv_all<-function(CID_wd,CID_name,EAD_wd,EAD_name){
  conditions <- c(
    expression(Delta*pgpB),
    expression(Delta*cdh),
    expression(Delta*cfa),
    expression(Delta*plsX),
    expression(Delta*opgB),
    expression(Delta*clsC),
    expression(Delta*pgpC),
    expression(Delta*pldC),
    expression(Delta*clsB),
    expression(Delta*clsA),
    expression(Delta*pldB),
    expression(Delta*pldA),
    expression(Delta*pgpA),
    expression(Delta*aas),
    expression(Delta*dgkA),
    "Wild Type"
    
  )
  Output<-vector(mode="list",2)
  names(Output)<-c("Barplot","Pie Chart")
  
  Directory <- CID_wd
  file <-CID_name
  
  dataframe <- data_import(Directory,file)
  data_CID <-dataframe$data
  
  Directory <- EAD_wd
  file <- EAD_name
  
  dataframe <- data_import(Directory,file)
  data_EAD <-dataframe$data
  
  
  data_CID<-filtering(data_CID)
  data_EAD<-filtering(data_EAD)
  
  CID<-matrix(ncol=16,nrow=2*nrow(data_CID))
  
  EAD<-matrix(ncol=16,nrow=2*nrow(data_EAD))
  ep<-c(c(36:80),c(84:86))
  sp<-c(c(87:131),c(135:137))
  
  
  sequence<-seq(1,48,3)
  i=1
  for(i in 1:16){
    j<-ep[i]
    ind_ep<-c(j,j+1,j+2)
    
    k<-sp[i]
    ind_sp<-c(k,k+1,k+2)
    
    for(h in 1:nrow(data_CID)){
      m_1 <- mean(as.numeric(data_CID[h,ind_ep]),na.rm=T)
      sd_1 <-sd(as.numeric(data_CID[h,ind_ep]),na.rm=T)
      cv_1 <- sd_1/m_1
      CID[h,i] <-cv_1
      
      m_2 <- mean(as.numeric(data_CID[h,ind_sp]),na.rm=T)
      sd_2 <-sd(as.numeric(data_CID[h,ind_sp]),na.rm=T)
      cv_2 <- sd_2/m_2
      CID[(h+nrow(data_CID)),i] <-cv_2
      
    }
    
    
  } 
  CID[CID=="NaN"]<-NA
  CID[CID==Inf]<-NA
  CID[abs(CID)>1.50]<-1.50
  par(mar=c(5,5,5,5))
  
  
  
  for(i in 1:16){
    j<-ep[i]
    ind_ep<-c(j,j+1,j+2)
    
    k<-sp[i]
    ind_sp<-c(k,k+1,k+2)
    
    for(h in 1:nrow(data_EAD)){
      m_1 <- mean(as.numeric(data_EAD[h,ind_ep]),na.rm=T)
      sd_1 <-sd(as.numeric(data_EAD[h,ind_ep]),na.rm=T)
      cv_1 <- sd_1/m_1
      EAD[h,i] <-cv_1
      
      m_2 <- mean(as.numeric(data_EAD[h,ind_sp]),na.rm=T)
      sd_2 <-sd(as.numeric(data_EAD[h,ind_sp]),na.rm=T)
      cv_2 <- sd_2/m_2
      EAD[(h+nrow(data_EAD)),i] <-cv_2
      
    }
    
    
  } 
  EAD[EAD=="NaN"]<-NA
  EAD[EAD==Inf]<-NA
  EAD[abs(EAD)>1.5]<-1.50
  
  
  
  
  d<-rbind(colMeans(CID,na.rm = T),colMeans(EAD,na.rm = T))*100
  
  layout(matrix(c(1,2),ncol=2),widths=c(8,2))
  par(cex.main=1.8,cex.lab=1.6,cex.axis=1.4,mar=c(8,6,5,0))
  barplot(d,beside = T,col=c("coral1","palegreen3"),names.arg = conditions,las=2,
          ylab = "Coefficient of Variation [%]",
          main="Mean CV: Biological Samples")
  par(mar=c(8,0,5,0))
  plot(1:10,1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("topleft",legend=c("CID","EAD"),cex=2,bty="n",col=c("coral1","palegreen3"),pch=15,pt.cex = 4,y.intersp=2)
  
  p<-recordPlot()
  Output[[1]]<-p
  plot.new()
  
  CID <- unlist(data.frame(CID))
  
  EAD<-unlist(data.frame(EAD))
  
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  
  x<-matrix(ncol=2,nrow=3)
  x[,1] <- c("CV < 20%","20% < CV < 100%", "CV > 100%")
  x[1,2] <- length(which(abs(CID)<0.2))
  x[2,2] <- length(which(abs(CID)>0.2 & abs(CID)<1 ))
  x[3,2] <- length(which(abs(CID)>1))
  
  # Create Data
  data <- data.frame(
    group=x[,1],
    value=x[,2]
  )
  
  # Compute the position of labels
  data <- data %>% 
    arrange(desc(as.numeric(value))) %>%
    mutate(prop = as.numeric(value) / sum(as.numeric(data$value)) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  
  color<-c("azure4",
           "chartreuse4","tomato3")
  
  data$group <- factor(data$group, levels = rev(as.character(data$group)))
  
  percentage <- round(data$prop,digits = 0)
  percentage[which(as.numeric(percentage)<5)]<- 0
  percentage[which(percentage>1)]<- paste0(percentage[which(percentage>1)]," %")
  percentage[which(percentage=="0")]<-""
  
  color_index <- nrow(data)-length(unique(color))
  if(color_index>0){
    col=c(rep("black",(color_index) ),unique(color) )
  }else{
    col=unique(color)
  }
  
  # Basic piechart
  
  
  
  
  pie1<-ggplot(data, aes(x="", y=prop, fill=group)) +
    labs(title = "Biological Samples")+
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="right") +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_text(aes(y = ypos, label = percentage), color = "white", size=6) +
    scale_fill_manual(values= col )+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(size=22))+
    theme(legend.key.size = unit(1.3, 'cm'))+
    theme(legend.title = element_text(size=20))+
    theme(legend.text = element_text(size=20))+
    theme(legend.direction="horizontal")+
    guides(fill=guide_legend(title=""))
  
  
  y<-matrix(ncol=2,nrow=3)
  y[,1] <- c("CV < 20%","20% < CV < 100%", "CV > 100%")
  y[1,2] <- length(which(abs(EAD)<0.2))
  y[2,2] <- length(which(abs(EAD)>0.2 & abs(EAD)<1 ))
  y[3,2] <- length(which(abs(EAD)>1))
  
  # Create Data
  data2 <- data.frame(
    group=y[,1],
    value1=y[,2]
  )
  
  # Compute the position of labels
  data2 <- data2 %>% 
    arrange(desc(as.numeric(value1))) %>%
    mutate(prop1 = as.numeric(value1) / sum(as.numeric(data2$value1)) *100) %>%
    mutate(ypos = cumsum(prop1)- 0.5*prop1 )
  
  
  color<-c("azure4",
           "chartreuse4","tomato3")
  
  data2$group <- factor(data2$group, levels = rev(as.character(data2$group)))
  
  percentage1 <- round(data2$prop,digits = 0)
  percentage1[which(as.numeric(percentage1)<5)]<- 0
  percentage1[which(percentage1>1)]<- paste0(percentage1[which(percentage1>1)]," %")
  percentage1[which(percentage1=="0")]<-""
  
  color_index <- nrow(data2)-length(unique(color))
  if(color_index>0){
    col=c(rep("black",(color_index) ),unique(color) )
  }else{
    col=unique(color)
  }
  
  # Basic piechart
  
  
  
  
  pie2<-ggplot(data2, aes(x="", y=prop1, fill=group)) +
    labs(title = "")+
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="right") +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_text(aes(y = ypos, label = percentage1), color = "white", size=6) +
    scale_fill_manual(values= col )+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(size=22))+
    theme(legend.key.size = unit(2, 'cm'))+
    theme(legend.title = element_text(size=20))+
    theme(legend.text = element_text(size=20))
  
  mylegend<-g_legend(pie1)
  grid.arrange(arrangeGrob(pie1+theme(legend.position = "none"),
                           pie2+theme(legend.position = "none"),nrow=2),
               mylegend, nrow=2,heights=c(8,1)) 
  
  p<-recordPlot()
  Output[[2]]<-p
  plot.new()
  
  return(Output)
  
}







cv_standard<-function(CID_wd,CID_name,EAD_wd,EAD_name){
  conditions <- c(
    expression(Delta*pgpB),
    expression(Delta*cdh),
    expression(Delta*cfa),
    expression(Delta*plsX),
    expression(Delta*opgB),
    expression(Delta*clsC),
    expression(Delta*pgpC),
    expression(Delta*pldC),
    expression(Delta*clsB),
    expression(Delta*clsA),
    expression(Delta*pldB),
    expression(Delta*pldA),
    expression(Delta*pgpA),
    expression(Delta*aas),
    expression(Delta*dgkA),
    "Wild Type"
    
  )
  Output<-vector(mode="list",2)
  names(Output)<-c("Barplot","Pie Chart")
  
  Directory <- CID_wd
  file <-CID_name
  
  dataframe <- data_import(Directory,file)
  data_CID <-dataframe$standard
  
  Directory <- EAD_wd
  file <- EAD_name
  
  dataframe <- data_import(Directory,file)
  data_EAD <-dataframe$standard
  
  
  
  CID<-matrix(ncol=16,nrow=2*nrow(data_CID))
  
  EAD<-matrix(ncol=16,nrow=2*nrow(data_EAD))
  ep<-c(c(36:80),c(84:86))
  sp<-c(c(87:131),c(135:137))
  
  
  sequence<-seq(1,48,3)
  i=1
  for(i in 1:16){
    j<-ep[i]
    ind_ep<-c(j,j+1,j+2)
    
    k<-sp[i]
    ind_sp<-c(k,k+1,k+2)
    
    for(h in 1:nrow(data_CID)){
      m_1 <- mean(as.numeric(data_CID[h,ind_ep]),na.rm=T)
      sd_1 <-sd(as.numeric(data_CID[h,ind_ep]),na.rm=T)
      cv_1 <- sd_1/m_1
      CID[h,i] <-cv_1
      
      m_2 <- mean(as.numeric(data_CID[h,ind_sp]),na.rm=T)
      sd_2 <-sd(as.numeric(data_CID[h,ind_sp]),na.rm=T)
      cv_2 <- sd_2/m_2
      CID[(h+nrow(data_CID)),i] <-cv_2
      
    }
    
    
  } 
  CID[CID=="NaN"]<-NA
  CID[CID==Inf]<-NA
  CID[abs(CID)>1.5]<-1.50
  
  par(mar=c(5,5,5,5))
  
  
  
  for(i in 1:16){
    j<-ep[i]
    ind_ep<-c(j,j+1,j+2)
    
    k<-sp[i]
    ind_sp<-c(k,k+1,k+2)
    
    for(h in 1:nrow(data_EAD)){
      m_1 <- mean(as.numeric(data_EAD[h,ind_ep]),na.rm=T)
      sd_1 <-sd(as.numeric(data_EAD[h,ind_ep]),na.rm=T)
      cv_1 <- sd_1/m_1
      EAD[h,i] <-cv_1
      
      m_2 <- mean(as.numeric(data_EAD[h,ind_sp]),na.rm=T)
      sd_2 <-sd(as.numeric(data_EAD[h,ind_sp]),na.rm=T)
      cv_2 <- sd_2/m_2
      EAD[(h+nrow(data_EAD)),i] <-cv_2
      
    }
    
    
  } 
  EAD[EAD=="NaN"]<-NA
  EAD[EAD==Inf]<-NA
  EAD[abs(EAD)>1.5]<-1.50
  
  
  
  
  
  d<-rbind(colMeans(CID,na.rm = T),colMeans(EAD,na.rm = T))*100
  
  layout(matrix(c(1,2),ncol=2),widths=c(8,2))
  par(cex.main=1.8,cex.lab=1.6,cex.axis=1.4,mar=c(8,6,5,0))
  barplot(d,beside = T,col=c("coral1","palegreen3"),names.arg = conditions,las=2,
          ylab = "Coefficient of Variation [%]",
          main="Mean CV: Internal Standards")
  par(mar=c(8,0,5,0))
  plot(1:10,1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("topleft",legend=c("CID","EAD"),cex=2,bty="n",col=c("coral1","palegreen3"),pch=15,pt.cex = 4,y.intersp=2)
  
  p<-recordPlot()
  Output[[1]]<-p
  plot.new()
  
  CID <- unlist(data.frame(CID))
  
  EAD<-unlist(data.frame(EAD))
  
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  
  x<-matrix(ncol=2,nrow=3)
  x[,1] <- c("CV < 20%","20% < CV < 100%", "CV > 100%")
  x[1,2] <- length(which(abs(CID)<0.2))
  x[2,2] <- length(which(abs(CID)>0.2 & abs(CID)<1 ))
  x[3,2] <- length(which(abs(CID)>1))
  
  # Create Data
  data <- data.frame(
    group=x[,1],
    value=x[,2]
  )
  
  # Compute the position of labels
  data <- data %>% 
    arrange(desc(as.numeric(value))) %>%
    mutate(prop = as.numeric(value) / sum(as.numeric(data$value)) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  
  color<-c("azure4",
           "chartreuse4","tomato3")
  
  data$group <- factor(data$group, levels = rev(as.character(data$group)))
  
  percentage <- round(data$prop,digits = 0)
  percentage[which(as.numeric(percentage)<5)]<- 0
  percentage[which(percentage>1)]<- paste0(percentage[which(percentage>1)]," %")
  percentage[which(percentage=="0")]<-""
  
  color_index <- nrow(data)-length(unique(color))
  if(color_index>0){
    col=c(rep("black",(color_index) ),unique(color) )
  }else{
    col=unique(color)
  }
  
  # Basic piechart
  
  
  
  
  pie1<-ggplot(data, aes(x="", y=prop, fill=group)) +
    labs(title = "Internal Standards")+
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="right") +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_text(aes(y = ypos, label = percentage), color = "white", size=6) +
    scale_fill_manual(values= col )+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(size=22))+
    theme(legend.key.size = unit(1.3, 'cm'))+
    theme(legend.title = element_text(size=20))+
    theme(legend.text = element_text(size=20))+
    theme(legend.direction="horizontal")+
    guides(fill=guide_legend(title=""))
  
  
  y<-matrix(ncol=2,nrow=3)
  y[,1] <- c("CV < 20%","20% < CV < 100%", "CV > 100%")
  y[1,2] <- length(which(abs(EAD)<0.2))
  y[2,2] <- length(which(abs(EAD)>0.2 & abs(EAD)<1 ))
  y[3,2] <- length(which(abs(EAD)>1))
  
  # Create Data
  data2 <- data.frame(
    group=y[,1],
    value1=y[,2]
  )
  
  # Compute the position of labels
  data2 <- data2 %>% 
    arrange(desc(as.numeric(value1))) %>%
    mutate(prop1 = as.numeric(value1) / sum(as.numeric(data2$value1)) *100) %>%
    mutate(ypos = cumsum(prop1)- 0.5*prop1 )
  
  
  color<-c("azure4",
           "chartreuse4","tomato3")
  
  data2$group <- factor(data2$group, levels = rev(as.character(data2$group)))
  
  percentage1 <- round(data2$prop,digits = 0)
  percentage1[which(as.numeric(percentage1)<5)]<- 0
  percentage1[which(percentage1>1)]<- paste0(percentage1[which(percentage1>1)]," %")
  percentage1[which(percentage1=="0")]<-""
  
  
  color_index <- nrow(data2)-length(unique(color))
  if(color_index>0){
    col=c(rep("black",(color_index) ),unique(color) )
  }else{
    col=unique(color)
  }
  
  # Basic piechart
  
  
  
  
  pie2<-ggplot(data2, aes(x="", y=prop1, fill=group)) +
    labs(title = "")+
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="right") +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_text(aes(y = ypos, label = percentage1), color = "white", size=6) +
    scale_fill_manual(values= col )+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(size=22))+
    theme(legend.key.size = unit(2, 'cm'))+
    theme(legend.title = element_text(size=20))+
    theme(legend.text = element_text(size=20))
  
  mylegend<-g_legend(pie1)
  grid.arrange(arrangeGrob(pie1+theme(legend.position = "none"),
                           pie2+theme(legend.position = "none"),nrow=2),
               mylegend, nrow=2,heights=c(8,1)) 
  
  p<-recordPlot()
  Output[[2]]<-p
  plot.new()
  
  return(Output)
  
}








cv_blanks<-function(CID_wd,CID_name,EAD_wd,EAD_name){
  conditions <- c(
    "EP Blanks","SP Blanks"
    
  )
  Output<-vector(mode="list",2)
  names(Output)<-c("Barplot","Pie Chart")
  
  Directory <- CID_wd
  file <-CID_name
  
  dataframe <- data_import(Directory,file)
  data_CID <-dataframe$data
  
  Directory <- EAD_wd
  file <- EAD_name
  
  dataframe <- data_import(Directory,file)
  data_EAD <-dataframe$data
  
  
  
  CID<-matrix(ncol=2,nrow=nrow(data_CID))
  
  EAD<-matrix(ncol=2,nrow=nrow(data_EAD))
  ep_CID<-c(138:143)
  sp_CID<-c(144:149)
  
  ep_EAD<-c(138:142)
  sp_EAD<-c(143:147)
  
  
  
  
  
  
  ind_ep<-ep_CID
  ind_sp<-sp_CID
  h<-1
  for(h in 1:nrow(data_CID)){
    m_1 <- mean(as.numeric(data_CID[h,ind_ep]),na.rm=T)
    sd_1 <-sd(as.numeric(data_CID[h,ind_ep]),na.rm=T)
    cv_1 <- sd_1/m_1
    CID[h,1] <-cv_1
    
    m_2 <- mean(as.numeric(data_CID[h,ind_sp]),na.rm=T)
    sd_2 <-sd(as.numeric(data_CID[h,ind_sp]),na.rm=T)
    cv_2 <- sd_2/m_2
    CID[h,2] <-cv_2
    
  }
  
  
  
  CID[CID=="NaN"]<-NA
  CID[CID==Inf]<-NA
  CID[abs(CID)>1.5]<-1.50
  
  par(mar=c(5,5,5,5))
  
  
  
  ind_ep<-ep_EAD
  ind_sp<-sp_EAD
  
  for(h in 1:nrow(data_EAD)){
    m_1 <- mean(as.numeric(data_EAD[h,ind_ep]),na.rm=T)
    sd_1 <-sd(as.numeric(data_EAD[h,ind_ep]),na.rm=T)
    cv_1 <- sd_1/m_1
    EAD[h,1] <-cv_1
    
    m_2 <- mean(as.numeric(data_EAD[h,ind_sp]),na.rm=T)
    sd_2 <-sd(as.numeric(data_EAD[h,ind_sp]),na.rm=T)
    cv_2 <- sd_2/m_2
    EAD[h,2] <-cv_2
    
  }
  
  
  
  EAD[EAD=="NaN"]<-NA
  EAD[EAD==Inf]<-NA
  EAD[abs(EAD)>1.5]<-1.50
  
  
  
  
  
  
  d<-rbind(colMeans(CID,na.rm = T),colMeans(EAD,na.rm = T))*100
  
  layout(matrix(c(1,2),ncol=2),widths=c(8,2))
  par(cex.main=1.8,cex.lab=1.6,cex.axis=1.4,mar=c(8,6,5,0))
  barplot(d,beside = T,col=c("coral1","palegreen3"),names.arg = conditions,las=2,
          ylab = "Coefficient of Variation [%]",
          main="Mean CV: Blank Samples")
  par(mar=c(8,0,5,0))
  plot(1:10,1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("topleft",legend=c("CID","EAD"),cex=2,bty="n",col=c("coral1","palegreen3"),pch=15,pt.cex = 4,y.intersp=2)
  
  p<-recordPlot()
  Output[[1]]<-p
  plot.new()
  
  CID <- unlist(data.frame(CID))
  
  EAD<-unlist(data.frame(EAD))
  
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  
  x<-matrix(ncol=2,nrow=3)
  x[,1] <- c("CV < 20%","20% < CV < 100%", "CV > 100%")
  x[1,2] <- length(which(abs(CID)<0.2))
  x[2,2] <- length(which(abs(CID)>0.2 & abs(CID)<1 ))
  x[3,2] <- length(which(abs(CID)>1))
  
  # Create Data
  data <- data.frame(
    group=x[,1],
    value=x[,2]
  )
  
  # Compute the position of labels
  data <- data %>% 
    arrange(desc(as.numeric(value))) %>%
    mutate(prop = as.numeric(value) / sum(as.numeric(data$value)) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  
  color<-c("azure4",
           "chartreuse4","tomato3")
  
  data$group <- factor(data$group, levels = rev(as.character(data$group)))
  
  percentage <- round(data$prop,digits = 0)
  percentage[which(percentage<5)]<-0
  percentage[which(percentage>1)]<- paste0(percentage[which(percentage>1)]," %")
  percentage[which(percentage=="0")]<-""
  
  
  color_index <- nrow(data)-length(unique(color))
  if(color_index>0){
    col=c(rep("black",(color_index) ),unique(color) )
  }else{
    col=unique(color)
  }
  
  # Basic piechart
  
  
  
  
  pie1<-ggplot(data, aes(x="", y=prop, fill=group)) +
    labs(title = "Blank Samples")+
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="right") +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_text(aes(y = ypos, label = percentage), color = "white", size=6) +
    scale_fill_manual(values= col )+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(size=22))+
    theme(legend.key.size = unit(1.3, 'cm'))+
    theme(legend.title = element_text(size=20))+
    theme(legend.text = element_text(size=20))+
    theme(legend.direction="horizontal")+
    guides(fill=guide_legend(title=""))
  
  
  y<-matrix(ncol=2,nrow=3)
  y[,1] <- c("CV < 20%","20% < CV < 100%", "CV > 100%")
  y[1,2] <- length(which(abs(EAD)<0.2))
  y[2,2] <- length(which(abs(EAD)>0.2 & abs(EAD)<1 ))
  y[3,2] <- length(which(abs(EAD)>1))
  
  # Create Data
  data2 <- data.frame(
    group=y[,1],
    value1=y[,2]
  )
  
  # Compute the position of labels
  data2 <- data2 %>% 
    arrange(desc(as.numeric(value1))) %>%
    mutate(prop1 = as.numeric(value1) / sum(as.numeric(data2$value1)) *100) %>%
    mutate(ypos = cumsum(prop1)- 0.5*prop1 )
  
  
  color<-c("azure4",
           "chartreuse4","tomato3")
  
  data2$group <- factor(data2$group, levels = rev(as.character(data2$group)))
  
  percentage1 <- round(data2$prop,digits = 0)
  percentage1[which(percentage1<5)]<-0
  percentage1[which(percentage1>1)]<- paste0(percentage1[which(percentage1>1)]," %")
  percentage1[which(percentage1=="0")]<-""
  
  
  color_index <- nrow(data2)-length(unique(color))
  if(color_index>0){
    col=c(rep("black",(color_index) ),unique(color) )
  }else{
    col=unique(color)
  }
  
  # Basic piechart
  
  
  
  
  pie2<-ggplot(data2, aes(x="", y=prop1, fill=group)) +
    labs(title = "")+
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="right") +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_text(aes(y = ypos, label = percentage1), color = "white", size=6) +
    scale_fill_manual(values= col )+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(size=22))+
    theme(legend.key.size = unit(2, 'cm'))+
    theme(legend.title = element_text(size=20))+
    theme(legend.text = element_text(size=20))
  
  mylegend<-g_legend(pie1)
  grid.arrange(arrangeGrob(pie1+theme(legend.position = "none"),
                           pie2+theme(legend.position = "none"),nrow=2),
               mylegend, nrow=2,heights=c(8,1)) 
  
  p<-recordPlot()
  Output[[2]]<-p
  plot.new()
  
  return(Output)
  
}



cv_vs_intensity<-function(CID_wd,CID_name,EAD_wd,EAD_name){
  Output<-vector(mode="list",2)
  names(Output)<-c("Barplot","Pie Chart")
  
  Directory <- CID_wd
  file <-CID_name
  
  dataframe <- data_import(Directory,file)
  data_CID <-dataframe$data
  
  Directory <- EAD_wd
  file <- EAD_name
  
  dataframe <- data_import(Directory,file)
  data_EAD <-dataframe$data
  
  
  data_CID<-filtering(data_CID)
  data_EAD<-filtering(data_EAD)
  
  CID<-matrix(ncol=16,nrow=2*nrow(data_CID))
  CID_int<-CID
  EAD<-matrix(ncol=16,nrow=2*nrow(data_EAD))
  EAD_int<-EAD
  ep<-c(c(36:80),c(84:86))
  sp<-c(c(87:131),c(135:137))
  
  
  sequence<-seq(1,48,3)
  i=1
  for(i in 1:16){
    j<-ep[i]
    ind_ep<-c(j,j+1,j+2)
    
    k<-sp[i]
    ind_sp<-c(k,k+1,k+2)
    
    for(h in 1:nrow(data_CID)){
      m_1 <- mean(as.numeric(data_CID[h,ind_ep]),na.rm=T)
      sd_1 <-sd(as.numeric(data_CID[h,ind_ep]),na.rm=T)
      cv_1 <- sd_1/m_1
      CID[h,i] <-cv_1
      CID_int[h,i]<-m_1
      
      m_2 <- mean(as.numeric(data_CID[h,ind_sp]),na.rm=T)
      sd_2 <-sd(as.numeric(data_CID[h,ind_sp]),na.rm=T)
      cv_2 <- sd_2/m_2
      CID[(h+nrow(data_CID)),i] <-cv_2
      CID_int[(h+nrow(data_CID)),i] <-m_2
      
    }
    
    
  } 
  
  dc<-cbind(unlist(data.frame(CID_int)),unlist(data.frame(CID)))
  
  dc[which(dc[,2]=="NaN"),]<-NA
  dc[which(abs(dc[,2])>1.50),]<-NA
  
  
  
  
  
  
  
  for(i in 1:16){
    j<-ep[i]
    ind_ep<-c(j,j+1,j+2)
    
    k<-sp[i]
    ind_sp<-c(k,k+1,k+2)
    
    for(h in 1:nrow(data_EAD)){
      m_1 <- mean(as.numeric(data_EAD[h,ind_ep]),na.rm=T)
      sd_1 <-sd(as.numeric(data_EAD[h,ind_ep]),na.rm=T)
      cv_1 <- sd_1/m_1
      EAD[h,i] <-cv_1
      EAD_int[h,i]<-m_1
      m_2 <- mean(as.numeric(data_EAD[h,ind_sp]),na.rm=T)
      sd_2 <-sd(as.numeric(data_EAD[h,ind_sp]),na.rm=T)
      cv_2 <- sd_2/m_2
      EAD[(h+nrow(data_EAD)),i] <-cv_2
      EAD[(h+nrow(data_EAD)),i] <-m_2
      
    }
    
    
  } 
  
  
  de<-cbind(unlist(data.frame(EAD_int)),unlist(data.frame(EAD)))
  
  de[which(de[,2]=="NaN"),]<-NA
  de[which(abs(de[,2])>1.5),]<-NA
  
  layout(matrix(c(1,2),ncol=2),widths=c(8,4))
  par(mar=c(5,5,5,1))
  plot(log(dc[,1]),abs((dc[,2]*100))^2,col="coral1",pch=16,cex=1.6,main="Avg Squared CV vs log(mean Intensity) ",
       ylab=expression(CV^2),
       xlab="log(mean Intensity)",
       cex.main=1.5,
       cex.lab=1.5,
       cex.axis=1.4)
  points(log(de[,1]),abs((de[,2]*100))^2,col="palegreen3",pch=16,cex=1.6)
  
  
  model_c<-lm(abs((dc[,2]*100))^2~log(dc[,1]))
  abline(model_c,lwd=3,col="coral3")
  
  model_e<-lm(abs((de[,2]*100))^2~log(de[,1]))
  abline(model_e,lwd=3,col="palegreen4")
  par(mar=c(5,0,5,0))
  plot(1:10,1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("topleft",legend=c("CID data","EAD data","CID linear fit","EAD linear fit"),cex=1.5,bty="n",col=c("coral1","palegreen3","coral3","palegreen4")
         ,pch=c(16,16,NA,NA),lwd=c(NA,NA,4,4),pt.cex = 3,y.intersp=1.8)
  
  p<-recordPlot()
  return(p)
  plot.new()
  return(Output)
}





###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

CID_vs_EAD<-function(EAD,CID,phase,fragmentation,Directory){
  Output<-vector(mode="list",3)
  names(Output) <- c("log vs. log",
                     "Piechart: IDs per class all",
                     "Venn Diagram unique IDs") 
  conditions <- c("Wilt Type",
                  expression(Delta*pgpB),
                  expression(Delta*cdh),
                  expression(Delta*cfa),
                  expression(Delta*plsX),
                  expression(Delta*opgB),
                  expression(Delta*clsC),
                  expression(Delta*pgpC),
                  expression(Delta*pldC),
                  expression(Delta*clsB),
                  expression(Delta*clsA),
                  expression(Delta*pldB),
                  expression(Delta*pldA),
                  expression(Delta*pgpA),
                  expression(Delta*aas),
                  expression(Delta*dgkA)
                  
  )
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  
  
  ######################################################################
  E<-EAD[,c(1,2,3,5,8)]
  
  master_peak_list <- CID[,c(1,2,3,5,8)]
  master_peak_list$CID<-c(1:nrow(CID))
  master_peak_list$EAD<-F
  mz_tol<-0.05
  rt_tol<-0.1 #min
  i=5
  for(i in 1:nrow(EAD)){
    RT<- as.numeric(E[i,2]) - as.numeric(master_peak_list[,2])
    RT<- abs(RT) < rt_tol
    MZ<- as.numeric(E[i,3]) - as.numeric(master_peak_list[,3])
    MZ <- abs(MZ) < mz_tol
    Adduct <- E[i,4]==master_peak_list[,4]
    
    score<-data.frame(RT=RT,MZ=MZ,Adduct=Adduct)
    score[score==T]<-1
    
    index<-which(rowSums(score)==3)
    
    if(length(index)>1){
      RT<- as.numeric(E[i,2]) - as.numeric(master_peak_list[,2])
      minimum <- which(RT==min(RT[index]))
      index <- minimum
    }
    
    if(length(index)!=0){
      master_peak_list[index,7]<-i
    }else if(length(index)==0){
      add<-matrix(c(E[i,],FALSE,i),nrow=1)
      colnames(add) <- colnames(master_peak_list)
      master_peak_list <- rbind(master_peak_list,add)
      master_peak_list <- data.frame(master_peak_list)
    }
    
  }
  
  
  
  
  master_peak_list <- master_peak_list[-which(master_peak_list$CID==F),]
  master_peak_list <- master_peak_list[-which(master_peak_list$EAD==F),]
  
  
  master<-cbind(EAD[unlist(master_peak_list$EAD),c(2,3,4,5,12)],CID[unlist(master_peak_list$CID),c(2,3,4,5,12)])
  master$`CID id`<- unlist(master_peak_list$CID)
  master$`EAD id`<- unlist(master_peak_list$EAD) 
  
  ######################################################################
  
  
  dc<-CID[unlist(master_peak_list$CID),s]
  lipids<-CID$`Metabolite name`[unlist(master_peak_list$CID)]
  dc[dc<0]<-0
  sequence<-seq(1,45,3)
  
  
  de<-EAD[unlist(master_peak_list$EAD),s]
  lipids<-EAD$`Metabolite name`[unlist(master_peak_list$EAD)]
  de[de<0]<-0
  sequence<-seq(1,45,3)
  
  
  
  
  
  ######################################################################
  
  
  
  title <- bquote(paste(.(phase),": ",log[2](FC[CID])*" vs."*log[2](FC[EAD]))) 
  df<-vector(mode="list",15)
  test<-matrix(ncol=4,nrow=8)
  test[,1]<-2
  test[1,]<-1
  test[7,]<-18
  test[8,]<-19
  test[2:6,2:4]<-matrix(c(3:17),ncol=3)
  layout(test,widths=c(1,10,10,10),heights=c(1.5,8,8,8,8,8,1,2))
  
  par(mar=c(0,2,0,4))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center",legend=title,cex=2.5,bty="n")
  par(mar=c(2,0,1,0))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  mtext(2,text = expression(log[2](FC[EAD])),line=-2.5,cex=1.3)
  i=1
  for(i in 1:15){
    name<-conditions[i+1]
    wt<-c(46,47,48)
    k<-sequence[i]
    ko<-c(k,k+1,k+2)
    dataframe_C <- dc[,c(wt,ko)]
    dataframeC<-matrix(as.numeric(as.matrix(dataframe_C)),ncol=6)
    colnames(dataframeC)<-colnames(dataframe_C)
    dataframe_C<-dataframeC
    
    
    w<-rowMeans(dataframe_C[,c(1:3)])
    k<-rowMeans(dataframe_C[,c(4:6)])
    
    FC<-k/w
    
    FC_log_C<-log2(FC)
    
    FC_log_C[FC_log_C==-Inf]<- min(FC_log_C)-1
    FC_log_C[FC_log_C==Inf]<- max(FC_log_C)+1
    
    
    
    pvalue <- rep(NA,nrow(dataframe_C))  
    
    for(j in 1:nrow(dataframe_C)){
      ww<-dataframe_C[j,c(1:3)]
      kk<-dataframe_C[j,c(4:6)]
      
      t<-t.test(ww,kk,alternative = "two.sided",conf.level = 0.95)
      value<- -log10(t$p.value)
      
      pvalue[j]<-value
    }
    
    P<-rep("16",nrow(dataframe_C))
    P[which(pvalue>0.05)]<-1

    dataframe_E <- de[,c(wt,ko)]
    dataframeE<-matrix(as.numeric(as.matrix(dataframe_E)),ncol=6)
    colnames(dataframeE)<-colnames(dataframe_E)
    dataframe_E<-dataframeE
    
    
    w<-rowMeans(dataframe_E[,c(1:3)])
    k<-rowMeans(dataframe_E[,c(4:6)])
    
    FC<-k/w
    
    FC_log_E<-log2(FC)
    
    FC_log_E[FC_log_E==-Inf]<- min(FC_log_E)-1
    FC_log_E[FC_log_E==Inf]<- max(FC_log_E)+1
    
    
    
    
    FF<-FC_log_E/FC_log_C
    
    index<-which(FF>0.7 & FF < 1.3)
    P[index]<-16
    
    par(mar=c(2,2,1,2),mgp=c(2,0.4,0),tcl=-0.2)
    plot(FC_log_C,FC_log_E,cex=1.5,ylab="",xlab="",col="black",cex.axis=1,pch=as.numeric(P))
    abline(coef = c(0,1),col="red")
    legend("topleft",legend = name,bty="n",cex=1.2, adj = 0.2 )
    
    
  }
  par(mar=c(0,2,0,4))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  mtext(3,text = expression(log[2](FC[CID])),line=-2.5,cex=1.3)
  par(mar=c(0,5,0,0))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center",legend = c("p-value < 0.05","p-value > 0.05","x = y line"),col=c("black","black","red"),horiz = T,bty="n",cex=2,pch=c(16,1,NA)
         ,pt.cex = 4,text.width = 0.2,lwd=c(NA,NA,2))
  
  
  
  p <- recordPlot()
  Output[[1]]<-p
  plot.new()
  
  ##################################################################################################
  # Piechart ID's
  
  
  # CID all:
  
  unique_CID <- unique(CID$Ontology)
  
  count_C <- data.frame(matrix(ncol=2,nrow=length(unique_CID)))
  colnames(count_C)<-c("group","value")
  for(i in 1:length(unique_CID)){
    name<-unique_CID[i]
    count_C[i,1] <- name
    index<-which(CID$Ontology==name)
    count_C[i,2]<-length(index)
    
  }
  
  
  
  
  
  
  # CID significantly enriched:
  
  SE<-c()
  for(i in 1:15){
    name<-conditions[i+1]
    wt<-c(46,47,48)
    k<-sequence[i]
    ko<-c(k,k+1,k+2)
    dataframe <- dc[,c(wt,ko)]
    dataframe1<-matrix(as.numeric(as.matrix(dataframe)),ncol=6)
    colnames(dataframe1)<-colnames(dataframe)
    dataframe<-dataframe1
    colSums(dataframe)
    
    w<-rowMeans(dataframe[,c(1:3)])
    k<-rowMeans(dataframe[,c(4:6)])
    
    FC<-k/w
    
    FC_log<-log2(FC)
    
    FC_log[FC_log==-Inf]<- min(FC_log)-1
    FC_log[FC_log==Inf]<- max(FC_log)+1
    
    pvalue <- rep(NA,nrow(dataframe))  
    
    for(j in 1:nrow(dataframe)){
      ww<-dataframe[j,c(1:3)]
      kk<-dataframe[j,c(4:6)]
      
      t<-t.test(ww,kk,alternative = "two.sided",conf.level = 0.95)
      value<- -log10(t$p.value)
      
      pvalue[j]<-value
    }
    
    m<-cbind(FC_log,pvalue)
    colnames(m)<-c("FC","pvalue")
    
    
    
    index<-which(abs(m[,1])>=1 & m[,2]>=0.05)
    
    SE<-c(SE,index)
    SE<-unique(SE)
    
  }
  
  CC<-CID[SE,]
  
  unique_CC <- unique(CC$Ontology)
  
  count_CC <- data.frame(matrix(ncol=2,nrow=length(unique_CC)))
  colnames(count_CC)<-c("group","value")
  for(i in 1:length(unique_CC)){
    name<-unique_CC[i]
    count_CC[i,1] <- name
    index<-which(CC$Ontology==name)
    count_CC[i,2]<-length(index)
    
  }
  
  
  # EAD all:
  
  unique_EAD <- unique(EAD$Ontology)
  
  count_E <- data.frame(matrix(ncol=2,nrow=length(unique_EAD)))
  colnames(count_E)<-c("group","value")
  for(i in 1:length(unique_EAD)){
    name<-unique_EAD[i]
    count_E[i,1] <- name
    index<-which(EAD$Ontology==name)
    count_E[i,2]<-length(index)
    
  }
  # EAD significantly enriched:
  
  SE<-c()
  for(i in 1:15){
    name<-conditions[i+1]
    wt<-c(46,47,48)
    k<-sequence[i]
    ko<-c(k,k+1,k+2)
    dataframe <- de[,c(wt,ko)]
    dataframe1<-matrix(as.numeric(as.matrix(dataframe)),ncol=6)
    colnames(dataframe1)<-colnames(dataframe)
    dataframe<-dataframe1
    colSums(dataframe)
    
    w<-rowMeans(dataframe[,c(1:3)])
    k<-rowMeans(dataframe[,c(4:6)])
    
    FC<-k/w
    
    FC_log<-log2(FC)
    
    FC_log[FC_log==-Inf]<- min(FC_log)-1
    FC_log[FC_log==Inf]<- max(FC_log)+1
    
    pvalue <- rep(NA,nrow(dataframe))  
    
    for(j in 1:nrow(dataframe)){
      ww<-dataframe[j,c(1:3)]
      kk<-dataframe[j,c(4:6)]
      
      t<-t.test(ww,kk,alternative = "two.sided",conf.level = 0.95)
      value<- -log10(t$p.value)
      
      pvalue[j]<-value
    }
    
    m<-cbind(FC_log,pvalue)
    colnames(m)<-c("FC","pvalue")
    
    
    
    index<-which(abs(m[,1])>=1 & m[,2]>=0.05)
    
    SE<-c(SE,index)
    SE<-unique(SE)
    
  }
  
  EE<-EAD[SE,]
  
  unique_EE <- unique(EE$Ontology)
  
  count_EE <- data.frame(matrix(ncol=2,nrow=length(unique_EE)))
  colnames(count_EE)<-c("group","value")
  for(i in 1:length(unique_EE)){
    name<-unique_EE[i]
    count_EE[i,1] <- name
    index<-which(EE$Ontology==name)
    count_EE[i,2]<-length(index)
    
  }
  
  
  U<-unique(c(unique_CID,unique_EAD))  
  
  count<-matrix(ncol=4,nrow=length(U))
  colnames(count) <- c("CID","EAD","CID s.","EAD s.")  
  rownames(count) <- U    
  i=2
  for(i in 1:length(U)){
    name <- U[i]
    index<-which(count_C$group==name)
    if(length(index)>0){
      count[i,1]<- count_C$value[index]
    }
    
    index<-which(count_E$group==name)
    if(length(index)>0){
      count[i,2]<- count_E$value[index]
    }
    
    index<-which(count_CC$group==name)
    if(length(index)>0){
      count[i,3]<- count_CC$value[index]
    }
    
    index<-which(count_EE$group==name)
    if(length(index)>0){
      count[i,4]<- count_EE$value[index]
    }
    
    
    
  }
  
  count[is.na(count)]<-0
  
  color<- color<-c("azure4",
                   "chartreuse4",
                   "tomato3",
                   "black",
                   "red",
                   "deeppink1",
                   "goldenrod4",
                   "blue3",
                   "gold1",
                   "steelblue3",
                   "yellowgreen",
                   "chartreuse1",
                   "lightblue4",
                   "firebrick4",
                   "magenta4",
                   "orange3")
  
  color<-color[1:nrow(count)]
  
  
  
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  
  
  
  # Create Data
  data1 <- data.frame(
    group=rownames(count),
    value=as.numeric(count[,1])
  )
  
  # Compute the position of labels
  data1 <- data1 %>% 
    arrange(desc(as.numeric(value))) %>%
    mutate(prop = as.numeric(value) / sum(as.numeric(data1$value)) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  
  
  
  data1$group <- factor(data1$group, levels = rev(as.character(data1$group)))
  
  percentage <- round(data1$prop,digits = 0)
  percentage[which(percentage<5)]<-0
  percentage[which(percentage>1)]<- paste0(percentage[which(percentage>1)]," %")
  percentage[which(percentage=="0")]<-""
  
  
  color_index <- nrow(data1)-length(unique(color))
  if(color_index>0){
    col=c(rep("black",(color_index) ),unique(color) )
  }else{
    col=unique(color)
  }
  
  # Basic piechart
  
  
  
  
  pie1<-ggplot(data1, aes(x="", y=prop, fill=group)) +
    labs(title = "CID All Lipids")+
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="right") +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_text(aes(y = ypos, label = percentage), color = "white", size=6) +
    scale_fill_manual(values= col )+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(size=22))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.title = element_text(size=20))+
    theme(legend.text = element_text(size=15))+
    theme(legend.direction="horizontal")+
    guides(fill=guide_legend(title=""))+
    theme(plot.margin = margin(0, 0, 0, 6, "cm"))
  
  
  # Create Data
  data2 <- data.frame(
    group2=rownames(count),
    value2=as.numeric(count[,2])
  )
  
  # Compute the position of labels
  data2 <- data2 %>% 
    arrange(desc(as.numeric(value2))) %>%
    mutate(prop2 = as.numeric(value2) / sum(as.numeric(data2$value2)) *100) %>%
    mutate(ypos = cumsum(prop2)- 0.5*prop2 )
  
  
  
  
  data2$group2 <- factor(data2$group2, levels = rev(as.character(data2$group2)))
  
  percentage2 <- round(data2$prop2,digits = 0)
  percentage2[which(percentage2<5)]<-0
  percentage2[which(percentage2>1)]<- paste0(percentage2[which(percentage2>1)]," %")
  percentage2[which(percentage2=="0")]<-""
  
  
  color_index <- nrow(data2)-length(unique(color))
  if(color_index>0){
    col=c(rep("black",(color_index) ),unique(color) )
  }else{
    col=unique(color)
  }
  
  # Basic piechart
  
  
  
  
  pie2<-ggplot(data2, aes(x="", y=prop2, fill=group2)) +
    labs(title = "EAD All Lipids")+
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="right") +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_text(aes(y = ypos, label = percentage2), color = "white", size=6) +
    scale_fill_manual(values= col )+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(size=22))+
    theme(legend.key.size = unit(2, 'cm'))+
    theme(legend.title = element_text(size=20))+
    theme(legend.text = element_text(size=20))+
    theme(plot.margin = margin(0, 6, 0, 0, "cm"))
  
  
  
  
  # Create Data
  data3 <- data.frame(
    group3=rownames(count),
    value3=as.numeric(count[,3])
  )
  
  # Compute the position of labels
  data3 <- data3 %>% 
    arrange(desc(as.numeric(value3))) %>%
    mutate(prop = as.numeric(value3) / sum(as.numeric(data3$value3)) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  
  
  
  data3$group3 <- factor(data3$group3, levels = rev(as.character(data3$group3)))
  
  percentage3 <- round(data3$prop,digits = 0)
  percentage3[which(percentage3<5)]<-0
  percentage3[which(percentage3>1)]<- paste0(percentage3[which(percentage3>1)]," %")
  percentage3[which(percentage3=="0")]<-""
  
  
  color_index <- nrow(data3)-length(unique(color))
  if(color_index>0){
    col=c(rep("black",(color_index) ),unique(color) )
  }else{
    col=unique(color)
  }
  
  # Basic piechart
  
  
  
  
  pie3<-ggplot(data3, aes(x="", y=prop, fill=group3)) +
    labs(title = "CID Significantly Enrich. Lipids")+
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="right") +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_text(aes(y = ypos, label = percentage3), color = "white", size=6) +
    scale_fill_manual(values= col )+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(size=22))+
    theme(legend.key.size = unit(2, 'cm'))+
    theme(legend.title = element_text(size=20))+
    theme(legend.text = element_text(size=20))+
    theme(plot.margin = margin(0, 0, 0, 6, "cm"))
  
  
  
  # Create Data
  data4 <- data.frame(
    group4=rownames(count),
    value4=as.numeric(count[,4])
  )
  
  # Compute the position of labels
  data4 <- data4 %>% 
    arrange(desc(as.numeric(value4))) %>%
    mutate(prop = as.numeric(value4) / sum(as.numeric(data4$value4)) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  
  
  
  data4$group4 <- factor(data4$group4, levels = rev(as.character(data4$group4)))
  
  percentage4 <- round(data4$prop,digits = 0)
  percentage4[which(percentage4<5)]<-0
  percentage4[which(percentage4>1)]<- paste0(percentage4[which(percentage4>1)]," %")
  percentage4[which(percentage4=="0")]<-""
  
  
  color_index <- nrow(data4)-length(unique(color))
  if(color_index>0){
    col=c(rep("black",(color_index) ),unique(color) )
  }else{
    col=unique(color)
  }
  
  # Basic piechart
  
  
  
  
  pie4<-ggplot(data4, aes(x="", y=prop, fill=group4)) +
    labs(title = "EAD Significantly Enrich. Lipids")+
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="right") +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_text(aes(y = ypos, label = percentage4), color = "white", size=6) +
    scale_fill_manual(values= col )+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(size=22))+
    theme(legend.key.size = unit(2, 'cm'))+
    theme(legend.title = element_text(size=20))+
    theme(legend.text = element_text(size=20))+
    theme(plot.margin = margin(0, 6, 0, 0, "cm"))
  
  
  
  mylegend<-g_legend(pie1)
  
  layout<-matrix(c(1,2,3,4,5,5),ncol=2,byrow=T)
  grid.arrange(pie1+theme(legend.position = "none"),
               pie2+theme(legend.position = "none"),
               pie3+theme(legend.position = "none"),
               pie4+theme(legend.position = "none"),
               mylegend,layout_matrix=layout ,heights=c(8,8,2),widths=c(8,8),top=textGrob("ID distribution",gp=gpar(fontsize=25))) 
  
  
  p <- recordPlot()
  Output[[2]]<-p
  par(mar=c(1,1,1,1))
  plot.new()
  #####################################################################################
  # overlap all
  
  venn<-function(CID,EAD,Directory,name,title){
    c<- characteristics_CID(CID)
    e<- characteristics_EAD(EAD)
    
    l_c <- na.omit(c$`Species level`)
    l_e <-na.omit(e$`Species level`)
    
    
    myCol =c("coral1","palegreen3")
    setwd(Directory)
    venn.diagram(
      x=list(l_c,l_e),
      category.names=c("CID","EAD"),filename = name,
      output=T,
      main=title,
      main.cex = 1.5,
      wd = 2,
      lty = 'blank',
      fill = myCol,
      
      # Numbers
      cex = 1.4,
      fontface = "bold",
      fontfamily = "sans",
      main.fontfamily = "sans" ,
      # # Set names
      cat.cex = 1.8,
      cat.fontface = "bold",
      
      cat.fontfamily = "sans"
    )  
    
    
    
  }
  
  name<-paste0(phase,"_","overlap_all.png")
  title<- paste0(phase,": ","Overlap all Lipid ID's")
  venn(CID,EAD,Directory,name,title)
  
  
  
  
  name<-paste0(phase,"_","overlap_enriched.png")
  title<- paste0(phase,": ","Overlap Enriched Lipid ID's")
  venn(CC,EE,Directory,name,title)
  
  
  c<- characteristics_CID(CID)
  e<- characteristics_EAD(EAD)
  
  l_c <- na.omit(c$`Species level`)
  l_e <-na.omit(e$`Species level`)
  
  myCol =c("coral1","palegreen3")
  
  C_C<-setdiff(l_c,l_e)
  E_E<-setdiff(l_e,l_c)
  
  myCol =c("coral1","palegreen3")
  
  
  m<-max(c(length(C_C),length(E_E)))
  x<-seq(1,3,length=m)
  y<-seq(1,m,length=m)
  layout(matrix(c(1,2),ncol=2),widths=c(6,6))
  par(mar=c(5,15,5,0))
  plot( x,y ,xlim=c(0.5,2.5),ylim=rev(range(y)),xaxt="n",yaxt="n",ylab="",xlab="",type="n")
  text(x=rep(1,length(C_C)),y=c(length(C_C):1),C_C,cex = 1.7,col=myCol[1])
  text(x=rep(2,length(E_E)),y=c(length(E_E):1),E_E,cex = 1.7,col=myCol[2])
  mtext(3,text = "Unique ID's all Lipids",cex=2)
  c<- characteristics_CID(CC)
  e<- characteristics_EAD(EE)
  
  l_c <- na.omit(c$`Species level`)
  l_e <-na.omit(e$`Species level`)
  
  myCol =c("coral1","palegreen3")
  
  C_C<-setdiff(l_c,l_e)
  E_E<-setdiff(l_e,l_c)
  
  
  
  
  par(mar=c(5,0,5,15))
  plot( x,y ,xlim=c(0.5,2.5),ylim=rev(range(y)),xaxt="n",yaxt="n",ylab="",xlab="",type="n")
  text(x=rep(1,length(C_C)),y=c(length(C_C):1),C_C,cex = 1.7,col=myCol[1])
  text(x=rep(2,length(E_E)),y=c(length(E_E):1),E_E,cex = 1.7,col=myCol[2])
  mtext(3,text = "Unique ID's enriched Lipids",cex=2)
  
  
  legend("topright",legend = c("CID","EAD") ,col = myCol ,cex=2.5,pch=15,bty="n" ,pt.cex = 3,xpd=TRUE,inset=c(-0.5,0))
  mtext(paste0(phase,": ","Unique identified Lipid species"), side = 3, line = -2, outer = TRUE,cex = 2.6)
  
  p <- recordPlot()
  Output[[3]]<-p
  plot.new()
  
  return(Output)
}



SN<-function(EAD,CID){
  sn_c<-CID$`S/N average`
  sn_e<-EAD$`S/N average`
  
  sn_c<-log10(as.numeric(sn_c))
  sn_e<-log10(as.numeric(sn_e))
  par(mfrow=c(1,1),mar=c(6,6,6,8),cex.lab=1.5,cex.axis=1.4,cex.main=1.6)
  boxplot(sn_c,sn_e,col=c("coral1","palegreen3"),names=c("",""))
  axis(1,at=c(1,2),label=c("CID","EAD"),cex.axis=3,padj=1)
  mtext(2,text=expression(log[10](Signal/Noise)),line=2,cex = 2)
  mtext(3,text = "Distribution of Signal/Noise",cex=2.5,line=1)
}




gaussian <- function(CID,EAD){
  E<-EAD[,c(1,2,3,5,8)]
  
  master_peak_list <- CID[,c(1,2,3,5,8)]
  master_peak_list$CID<-c(1:nrow(CID))
  master_peak_list$EAD<-F
  mz_tol<-0.01
  rt_tol<-0.1 #min
  i=5
  for(i in 1:nrow(EAD)){
    RT<- as.numeric(E[i,2]) - as.numeric(master_peak_list[,2])
    RT<- abs(RT) < rt_tol
    MZ<- as.numeric(E[i,3]) - as.numeric(master_peak_list[,3])
    MZ <- abs(MZ) < mz_tol
    Adduct <- E[i,4]==master_peak_list[,4]
    
    score<-data.frame(RT=RT,MZ=MZ,Adduct=Adduct)
    score[score==T]<-1
    
    index<-which(rowSums(score)==3)
    
    if(length(index)>1){
      RT<- as.numeric(E[i,2]) - as.numeric(master_peak_list[,2])
      minimum <- which(RT==min(RT[index]))
      index <- minimum
    }
    
    if(length(index)!=0){
      master_peak_list[index,7]<-i
    }else if(length(index)==0){
      add<-matrix(c(E[i,],FALSE,i),nrow=1)
      colnames(add) <- colnames(master_peak_list)
      master_peak_list <- rbind(master_peak_list,add)
      master_peak_list <- data.frame(master_peak_list)
    }
    
  }
  
  
  
  
  master_peak_list <- master_peak_list[-which(master_peak_list$CID==F),]
  master_peak_list <- master_peak_list[-which(master_peak_list$EAD==F),]
  
  
  master<-cbind(CID[unlist(master_peak_list$CID),c(2,3,4,5,12,86,137)],
                EAD[unlist(master_peak_list$EAD),c(2,3,4,5,12,86,137)])
  master$`CID id`<- unlist(master_peak_list$CID)
  master$`EAD id`<- unlist(master_peak_list$EAD) 
  
  class <- unique(master$Ontology)
  count_1<-matrix(ncol=4,nrow=length(class))
  count_2 <- matrix(ncol=4,nrow=length(class))
  
  colnames(count_1) <-colnames(count_2) <- c("CID (row)","CID (col)","EAD (row)","EAD (col)")
  rownames(count_1) <-rownames(count_2) <- class
  i=5
  for(i in 1:length(class)){
    name<-class[i]
    index<-which(master$Ontology==name)
    
    
    
    d <- as.matrix(data.frame(EP=as.numeric(master[,6]),SP=as.numeric(master[,7])))
    d<-d[index,]
    if(length(index)==1){
      d<-matrix(d,nrow=1)
    }
    
    ind_c<-master$`CID id`[index]
    if(nrow(d)>=2){
      r <- which(d == max(d), arr.ind = TRUE)
      d[r[1],r[2]] <- 0
      
      count_1[i,1] <- ind_c[r[1]]
      count_1[i,2] <- r[2]
      
      
      r <- which(d == max(d), arr.ind = TRUE)
      
      
      count_2[i,1] <- ind_c[r[1]]
      count_2[i,2] <- r[2]
    }else if(nrow(d)==1){
      r<-which(d==max(d))
      count_1[i,1]<-ind_c
      count_1[i,2] <- r
    }
    
    
    d <- as.matrix(data.frame(EP=as.numeric(master[,13]),SP=as.numeric(master[,14])))
    d<-d[index,]
    if(length(index)==1){
      d<-matrix(d,nrow=1)
    }
    ind_e<-master$`EAD id`[index]
    
    if(nrow(d)>=2){
      r <- which(d == max(d), arr.ind = TRUE)
      d[r[1],r[2]] <- 0
      
      count_1[i,3] <- ind_e[r[1]]
      count_1[i,4] <- r[2]
      
      
      
      r <- which(d == max(d), arr.ind = TRUE)
      
      
      count_2[i,3] <- ind_e[r[1]]
      count_2[i,4] <- r[2]
    }else if(nrow(d)==1){
      r<-which(d==max(d))
      count_1[i,3]<-ind_e
      count_1[i,4] <- r
    }
    
  }
  
  count<-rbind(count_1,count_2)
  count<-count[-which(is.na(count[,1])),]
  count<-count[order(rownames(count)),]
  wd<-getwd()
  
  
  
  Output<-vector(mode="list",4)
  names(Output) <-c("CID TEIL I",
                    "CID TEIL II",
                    "EAD TEIL I",
                    "EAD TEIL II")
  
  fit<-function(x,y,pp,method,rt){
    
    param<-pp
    objective<-function(param,x,y){
      h<-param[1]
      tau<-param[2]
      sigma<-param[3]
      mu<-param[4]
      
      y_pred = (h*sigma/tau) * (sqrt(pi/2))*exp( ((mu-x)/tau) + ((sigma^2)/(2*(tau^2) ) ) ) * erfc(  (1/sqrt(2))*( ((mu-x)/sigma) +(sigma/tau)    )  )
      
      sum((y-y_pred)^2)
    }
    
    
    
    f <- optim(par=pp,fn=objective,x=x,y=y,method = method,control = list(maxit = 1000))
    
    h<-f$par[1]
    tau<-f$par[2]
    sigma<-f$par[3]
    mu<-f$par[4]
    
    x_new<-seq(0.5*min(x),1.5*max(x),length=10000*length(x))
    y_new = (h*sigma/tau) * (sqrt(pi/2))*exp( ((mu-x_new)/tau) + ((sigma^2)/(2*(tau^2) ) ) ) * erfc(  (1/sqrt(2))*( ((mu-x_new)/sigma) +(sigma/tau)    )  )
    
    xx<-rt
    yy <- (h*sigma/tau) * (sqrt(pi/2))*exp( ((mu-xx)/tau) + ((sigma^2)/(2*(tau^2) ) ) ) * erfc(  (1/sqrt(2))*( ((mu-xx)/sigma) +(sigma/tau)    )  )
    points(x_new,y_new,col="red",type="l")
    
    return(yy)
    
  }
  
  
  dotscore<-function(rt,ref,gaussian,col){
    par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
    plot(0,0,xlim=c(min(rt),max(rt)+5),ylim=c(-1.05*max(gaussian),1.05*max(ref)),xaxt="n",yaxt="n",
         main = "Dot Product Score")
    abline(h=0)
    points(rt,ref,type="h",col=col,lwd=2)
    points(rt,-gaussian,type="h",col="blue",lwd=2)
    
    dot_product <- sum(ref * gaussian)
    
    # Calculate the magnitude of each vector
    magnitude1 <- norm(ref, "2")
    magnitude2 <- norm(gaussian, "2")
    
    # Calculate the spectral similarity normalized dot product score
    spectral_score <- dot_product / (magnitude1 * magnitude2)
    legend("topright",legend = paste0("DPS: ",round(spectral_score,digits = 3)),bty="n",cex=1.2)
    
    
  }
  
  setwd("C:/FUNCTION_R_DATA")
  CID_S <- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_CID.mzML",mode = "onDisk") 
  
  EAD_S<- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_EAD.mzML",mode = "onDisk")
  
  
  
  test<-matrix(ncol=5,nrow=7)
  test[,1]<-2
  test[1,]<-1
  test[7,]<-23
  test[2:6,2:5]<- matrix(c(3:22),ncol=4,byrow = T)
  layout(test,widths=c(1,10,10,10,10),heights=c(1,8,8,8,8,8,2))
  
  par(mar=c(0,2,0,4))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center",legend="CID: Exponential Modiefied Gaussian Similarity",cex=2,bty="n")
  par(mar=c(2,0,1,0))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  mtext(2,text = expression(intensity),line=-1.8,cex=1.3)
  
  
  # PE 36:2
  title <- "PE 36:2"
  title <- paste0("CID: ",title)
  XIC<-CID[114,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(CID_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  pp<-c(max(smoothed$y),1,9,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  
  R<-as.numeric(rt)
  #fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"coral1")
  
  
  
  
  
  
  # PE 35:1
  title <- "PE 35:1"
  title <- paste0("CID: ",title)
  XIC<-CID[104,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(CID_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,9,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  #fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"coral1")  
  
  
  
  
  
  
  
  
  
  ###############################################################################################
  
  
  # CoQ7
  
  title <- "CoQ7"
  title <- paste0("CID: ",title)
  XIC<-CID[55,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(CID_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  pp<-c(max(smoothed$y),1,5,smoothed$x[which(smoothed$y==max(smoothed$y))])      
  
  R<-as.numeric(rt)
  #fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"coral1")
  
  
  # CoQ8
  title <- "CoQ8"
  title <- paste0("CID: ",title)
  XIC<-CID[101,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(CID_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  plot(smoothed$x,smoothed$y,col="coral1",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,5,smoothed$x[which(smoothed$y==max(smoothed$y))]+1)
  
  R<-as.numeric(rt)
  #fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"coral1")  
  
  
  
  
  ###############################################################################################
  rm(CID_S)
  rm(EAD_S)
  CID_S <- readMSData(files = "22_12_22_Exp_004_SP_Sam_4_2_CID.mzML",mode = "onDisk") 
  
  EAD_S<- readMSData(files = "22_12_22_Exp_004_SP_Sam_4_2_EAD.mzML",mode = "onDisk")
  
  # PG 33:1
  title <- "PG 33:1"
  title <- paste0("CID: ",title)
  XIC<-CID[119,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(CID_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,col="coral1",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  
  plot(smoothed$x,smoothed$y,col="coral1",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,9,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  #fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"coral1")  
  
  
  par(mar=c(0,2,2,4))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  mtext(3,text = "RT [s]",line=0,cex=1.3)
  
  p<-recordPlot()
  Output[[1]]<-p
  par(mar=c(0,0,0,0))
  plot.new()
  
  while (dev.cur()>1) dev.off()
  
  
  
  
  
  test<-matrix(ncol=5,nrow=7)
  test[,1]<-2
  test[1,]<-1
  test[7,]<-23
  test[2:6,2:5]<- matrix(c(3:22),ncol=4,byrow = T)
  layout(test,widths=c(1,10,10,10,10),heights=c(1,8,8,8,8,8,2))
  
  par(mar=c(0,2,0,4))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center",legend="CID: Exponential Modiefied Gaussian Similarity",cex=2,bty="n")
  par(mar=c(2,0,1,0))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  mtext(2,text = expression(intensity),line=-1.8,cex=1.3)
  
  
  
  
  
  # PG 35:1 
  title <- "PG 35:1"
  title <- paste0("CID: ",title)
  XIC<-CID[138,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(CID_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,col="coral1",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  
  plot(smoothed$x,smoothed$y,col="coral1",cex=2,pch=16,yaxt="n",ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,9,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"coral1")
  
  
  
  ###############################################################################################
  rm(CID_S)
  rm(EAD_S)
  CID_S <- readMSData(files = "22_12_22_Exp_004_SP_Sam_9_2_CID.mzML",mode = "onDisk") 
  
  EAD_S<- readMSData(files = "22_12_22_Exp_004_SP_Sam_9_2_EAD.mzML",mode = "onDisk")
  
  # DG 33:1
  title <- "DG 33:1"
  title <- paste0("CID: ",title)
  XIC<-CID[25,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(CID_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,col="coral1",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  plot(smoothed$x,smoothed$y,col="coral1",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,9,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"coral1")
  
  
  
  # DG 34:2
  title <- "DG 34:2"
  title <- paste0("CID: ",title)
  XIC<-CID[33,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(CID_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  
  plot(smoothed$x,smoothed$y,col="coral1",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,9,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"coral1")
  
  
  
  # TG 48:1
  title <- "TG 48:1"
  title <- paste0("CID: ",title)
  XIC<-CID[149,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(CID_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,col="coral1",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  plot(smoothed$x,smoothed$y,col="coral1",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y)+1000,1,10,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"coral1")  
  
  
  # TG 50:2
  title <- "TG 50:2"
  title <- paste0("CID: ",title)
  XIC<-CID[160,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(CID_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="coral1",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,3,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"coral1")  
  
  
  par(mar=c(0,2,2,4))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  mtext(3,text = "RT [s]",line=0,cex=1.3)
  
  p<-recordPlot()
  Output[[2]]<-p
  par(mar=c(0,0,0,0))
  plot.new()
  
  while (dev.cur()>1) dev.off()
  
  ###########################################################
  ###########################################################
  ###########################################################
  ###########################################################
  
  
  setwd("C:/FUNCTION_R_DATA")
  CID_S <- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_CID.mzML",mode = "onDisk") 
  
  EAD_S<- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_EAD.mzML",mode = "onDisk")
  
  test<-matrix(ncol=5,nrow=7)
  test[,1]<-2
  test[1,]<-1
  test[7,]<-23
  test[2:6,2:5]<- matrix(c(3:22),ncol=4,byrow = T)
  layout(test,widths=c(1,10,10,10,10),heights=c(1,8,8,8,8,8,2))
  
  par(mar=c(0,2,0,4))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center",legend="EAD: Exponential Modiefied Gaussian Similarity",cex=2,bty="n")
  par(mar=c(2,0,1,0))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  mtext(2,text = expression(intensity),line=-1.8,cex=1.3)
  
  
  
  # PE 36:2
  title <- "PE 36:2"
  title <- paste0("EAD: ",title)
  n<-master$`EAD id`[which(master$`CID id`==114)]
  XIC<-EAD[n,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(EAD_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  pp<-c(max(smoothed$y),1,9,smoothed$x[which(smoothed$y==max(smoothed$y))])      
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"palegreen3")
  
  
  
  
  # PE 35:1
  title <- "PE 35:1"
  title <- paste0("EAD: ",title)
  n<-master$`EAD id`[which(master$`CID id`==104)]
  XIC<-EAD[n,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(EAD_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,9,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"palegreen3")  
  
  
  
  
  
  
  
  
  
  ###############################################################################################
  
  
  # CoQ7
  
  title <- "CoQ7"
  title <- paste0("EAD: ",title)
  n<-master$`EAD id`[which(master$`CID id`==55)]
  
  XIC<-EAD[n,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.01),(as.numeric(XIC[1,2])+0.01))
  chr_raw <- chromatogram(EAD_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  pp<-c(max(smoothed$y),0.5,5,smoothed$x[which(smoothed$y==max(smoothed$y))])      
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"palegreen3")  
  
  
  
  
  # CoQ8
  title <- "CoQ8"
  title <- paste0("EAD: ",title)
  n<-master$`EAD id`[which(master$`CID id`==101)]
  XIC<-EAD[n,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(EAD_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  plot(smoothed$x,smoothed$y,col="palegreen3",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,5,smoothed$x[which(smoothed$y==max(smoothed$y))]+1)
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"palegreen3")  
  
  
  
  
  
  
  
  ###############################################################################################
  rm(CID_S)
  rm(EAD_S)
  CID_S <- readMSData(files = "22_12_22_Exp_004_SP_Sam_4_2_CID.mzML",mode = "onDisk") 
  
  EAD_S<- readMSData(files = "22_12_22_Exp_004_SP_Sam_4_2_EAD.mzML",mode = "onDisk")
  
  # PG 33:1
  title <- "PG 33:1"
  title <- paste0("EAD: ",title)
  n<-master$`EAD id`[which(master$`CID id`==119)]
  XIC<-EAD[n,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(EAD_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,col="palegreen3",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  
  plot(smoothed$x,smoothed$y,col="palegreen3",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,5,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"palegreen3")  
  
  
  par(mar=c(0,2,2,4))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  mtext(3,text = "RT [s]",line=0,cex=1.3)
  
  p<-recordPlot()
  Output[[3]]<-p
  par(mar=c(0,0,0,0))
  plot.new()
  
  while (dev.cur()>1) dev.off()
  
  
  
  
  test<-matrix(ncol=5,nrow=7)
  test[,1]<-2
  test[1,]<-1
  test[7,]<-23
  test[2:6,2:5]<- matrix(c(3:22),ncol=4,byrow = T)
  layout(test,widths=c(1,10,10,10,10),heights=c(1,8,8,8,8,8,2))
  
  par(mar=c(0,2,0,4))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center",legend="EAD: Exponential Modiefied Gaussian Similarity",cex=2,bty="n")
  par(mar=c(2,0,1,0))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  mtext(2,text = expression(intensity),line=-1.8,cex=1.3)
  
  
  
  
  
  # PG 35:1 
  title <- "PG 35:1"
  title <- paste0("EAD: ",title)
  n<-master$`EAD id`[which(master$`CID id`==138)]
  XIC<-EAD[n,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(EAD_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,col="palegreen3",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  
  plot(smoothed$x,smoothed$y,col="palegreen3",cex=2,pch=16,yaxt="n",ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,9,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"palegreen3")  
  
  
  
  ###############################################################################################
  rm(CID_S)
  rm(EAD_S)
  CID_S <- readMSData(files = "22_12_22_Exp_004_SP_Sam_9_2_CID.mzML",mode = "onDisk") 
  
  EAD_S<- readMSData(files = "22_12_22_Exp_004_SP_Sam_9_2_EAD.mzML",mode = "onDisk")
  
  # DG 33:1
  title <- "DG 33:1"
  title <- paste0("EAD: ",title)
  n<-master$`EAD id`[which(master$`CID id`==25)]
  
  XIC<-EAD[n,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(EAD_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,col="palegreen3",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  plot(smoothed$x,smoothed$y,col="palegreen3",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,9,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"palegreen3")  
  
  
  
  # DG 34:2
  title <- "DG 34:2"
  title <- paste0("EAD: ",title)
  n<-master$`EAD id`[which(master$`CID id`==33)]
  
  XIC<-EAD[n,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(EAD_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  
  plot(smoothed$x,smoothed$y,col="palegreen3",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,9,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"palegreen3")  
  
  
  
  # TG 48:1
  title <- "TG 48:1"
  title <- paste0("EAD: ",title)
  n<-master$`EAD id`[which(master$`CID id`==149)]
  
  XIC<-EAD[n,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(EAD_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,col="palegreen3",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  plot(smoothed$x,smoothed$y,col="palegreen3",yaxt="n",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,10,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"palegreen3")  
  
  
  # TG 50:2
  title <- "TG 50:2"
  title <- paste0("EAD: ",title)
  n<-master$`EAD id`[which(master$`CID id`==160)]
  XIC<-EAD[n,c(2,3,4,137)]
  rtr<-c( (as.numeric(XIC[1,1])-0.2)*60,(as.numeric(XIC[1,1])+0.2)*60)
  
  mzr<-c( (as.numeric(XIC[1,2])-0.08),(as.numeric(XIC[1,2])+0.08))
  chr_raw <- chromatogram(EAD_S, mz = mzr, rt = rtr)
  
  rt<-rtime(chr_raw[1,1])
  int<-intensity(chr_raw[1,1])
  par(mar=c(2,2,2,2),cex.lab=1,cex.axis=1,cex.main=1)
  plot(rt,int,col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Raw Data"),xlim=rtr)
  
  rt<-rt[!is.na(int)]
  int<-int[!is.na(int)]
  smoothed <- smooth.spline(rt,int,spar=0.5)
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Savitzky Golay filter"),xlim=rtr)
  
  
  
  
  
  plot(smoothed$x,smoothed$y,yaxt="n",col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",main=paste0(title," Exp. modified Gaussian"),xlim=rtr)
  
  
  
  pp<-c(max(smoothed$y),1,3,smoothed$x[which(smoothed$y==max(smoothed$y))])
  
  R<-as.numeric(rt)
  G<-  fit(smoothed$x,smoothed$y,pp,method="SANN",rt=R)
  
  dotscore(R,int,G,"palegreen3")  
  
  
  par(mar=c(0,2,2,4))
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  mtext(3,text = "RT [s]",line=0,cex=1.3)
  
  p<-recordPlot()
  Output[[4]]<-p
  par(mar=c(0,0,0,0))
  plot.new()
  
  while (dev.cur()>1) dev.off()
  
  return(Output)
  
  
}



Chord_plot<-function(data,phase,fragmentation){
  title <- paste0("Chord Diagram:  (",fragmentation,"-",phase,")")
  conditions <- c(
    expression(Delta*pgpB),
    expression(Delta*cdh),
    expression(Delta*cfa),
    expression(Delta*plsX),
    expression(Delta*opgB),
    expression(Delta*clsC),
    expression(Delta*pgpC),
    expression(Delta*pldC),
    expression(Delta*clsB),
    expression(Delta*clsA),
    expression(Delta*pldB),
    expression(Delta*pldA),
    expression(Delta*pgpA),
    expression(Delta*aas),
    expression(Delta*dgkA)
    
  )
  
  color<-c("azure4","blue4","red1","cyan3","chartreuse2","deeppink","cyan4",
           "chartreuse4","orange4","darkorchid3","brown4",
           "darkorange","darkmagenta","yellow3","deeppink4")
  
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  
  if(fragmentation=="CID"){
    ch<-characteristics_CID(data)
  }else if(fragmentation=="EAD"){
    ch<- characteristics_EAD(data)
  }
  
  
  d<-data[,s]
  lipids<-ch[,2]
  class<-ch[,1]
  lipids[which(is.na(lipids)==T)]<- "Other"
  d[d<0]<-0
  sequence<-seq(1,45,3)
  nrow(data)
  keep <- rowSums(d[,c(46:48)],na.rm=T) > 0
  d <- d[keep,]
  nrow(d)
  lipids<-lipids[keep]
  class<-class[keep]
  
  
  
  
  
  
  FoldChange<-p_value<-Result<-matrix(ncol=15,nrow=nrow(d))
  
  i=1
  for(i in 1:15){
    wt<-c(46,47,48)
    k<-sequence[i]
    ko<-c(k,k+1,k+2)
    dataframe <- d[,c(wt,ko)]
    dataframe1<-matrix(as.numeric(as.matrix(dataframe)),ncol=6)
    colnames(dataframe1)<-colnames(dataframe)
    dataframe<-dataframe1
    
    w<-rowMeans(dataframe[,c(1:3)])
    k<-rowMeans(dataframe[,c(4:6)])
    
    FC<-k/w
    
    FC_log<-log2(FC)
    
    FC_log[FC_log==-Inf]<- min(FC_log)-1
    FC_log[FC_log==Inf]<- max(FC_log)+1
    
    pvalue <- rep(NA,nrow(dataframe))  
    
    for(j in 1:nrow(dataframe)){
      ww<-dataframe[j,c(1:3)]
      kk<-dataframe[j,c(4:6)]
      
      t<-t.test(ww,kk,alternative = "two.sided",conf.level = 0.95)
      value<- (t$p.value)
      
      
      
      pvalue[j] <-value
    }
    
    FoldChange[,i]<-FC_log
    p_value[,i]<-pvalue
    
    FoldChange[FoldChange=="NaN"] <- NA
    
    p_value[p_value=="NaN"] <- NA
    
    
  }
  
  FC_final<- rep(0,length(unlist(FoldChange)))
  
  index_FC<-which(abs(unlist(FoldChange))>1) 
  index_pval<-which(unlist(p_value)<0.05)
  index<-intersect(index_FC,index_pval)
  FC_final[index]<-1
  result<-matrix(unlist(FC_final),ncol=15)
  
  
  result<-t(result)
  result[result==1]<-10000
  rownames(result)<-sort(as.character(c(1:15)))
  
  colnames(result) <-lipids
  
  
  
  nm = (unlist(dimnames(result)))
  group = structure(c(1:15,class), names = nm)
  group
  
  l<-unique(group)
  n<-c(1:length(l))
  
  c<-rep(NA,length(group))
  for(i in 1:length(l)){
    save<-l[i]
    index<-which(group==save)
    c[index]<-n[i]
  }
  
  
  grid.col=structure(c,names=nm)
  grid.col[1:15] <-color
  par(mar=c(0,0,0,0))
  circos.clear()
  circos.par(gap.degree = 0.01)
  chordDiagram(result, group = group,grid.col = grid.col,annotationTrack = c("grid"),big.gap = 2,small.gap = 0,
               annotationTrackHeight = c(0.1),preAllocateTracks = list(track.height = mm_h(20)),link.lwd = 20,directional = -1,
               link.arr.width=2)
  
  
  
  
  index<-data.frame(samples=c(1:15),
                    real=sort(as.character(1:15)))
  label<-get.all.sector.index()
  for(i in 1:15){
    test<-label[i]
    j<-which(index$real==test)
    j<-as.numeric(j)
    label[i]<-conditions[j]
    
  }
  
  
  label[1:15] <-conditions
  
  for(i in 1:15) {
    si<-get.all.sector.index()[i]
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), label[i], sector.index = si, track.index = 1, 
                facing = "clockwise", niceFacing = TRUE, col = "black",adj = c(0.5, 0.5),cex=1.5)
  }
  
  
  for(i in 16:length(label)) {
    si<-get.all.sector.index()[i]
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), label[i], sector.index = si, track.index = 1, 
                facing = "clockwise", niceFacing = TRUE, col = "black",adj = c(0.5, 0),cex=.8)
  }
  
  mtext(3,text=title,line=-1.5,cex=1.5,font=2)
  
  
  p<-recordPlot()
  plot.new()  
  
  return(p)
}







plot_DG_16_1_18_1<-function(type=T){
  Output<-vector(mode="list",2)
  setwd("C:/FUNCTION_R_DATA")
  
  CID_15 <- readMSData(files = "22_12_22_Exp_004_SP_Sam_9_2_CID.mzML",mode = "onDisk") 
  
  EAD_15<- readMSData(files = "22_12_22_Exp_004_SP_Sam_9_2_EAD.mzML",mode = "onDisk")
  
  
  
  CID_12 <- readMSData(files = "22_12_22_Exp_004_EP_Sam_6_1_CID.mzML",mode = "onDisk") 
  
  EAD_12<- readMSData(files = "22_12_22_Exp_004_EP_Sam_6_1_EAD.mzML",mode = "onDisk")
  
  
  
  CID_9 <- readMSData(files = "22_12_22_Exp_004_EP_Sam_9_2_CID.mzML",mode = "onDisk") 
  
  EAD_9<- readMSData(files = "22_12_22_Exp_004_EP_Sam_9_2_EAD.mzML",mode = "onDisk")
  
  
  
  stoch_line<-function(x1,y1,x2,y2,color,N=40){
    a<-(y2-y1)/(x2-x1)
    b<-y1-a*x1
    x<-seq(x1,x2,length=N)
    y<-a*x+b
    
    
    for(i in 1:length(y)){
      random<-runif(1,min=0,max=2)
      r<-runif(1,min=0,max=1)
      
      if(r>=0.4){
        d<-1
      }else{ d<- -1}
      
      value<-y[i]+(random*d)
      
      if(value<1){
        value<- -1*value
      }
      y[i]<-value
    }
    
    y[1]<-y1
    y[length(y)] <- y2
    
    points(x,y,col=color,lwd=2,type="l")
    
  }
  
  # DG 16:1(9Z)/18:1(9)
  
  # CID
  rt <- 4.65
  delta <- 10
  rtr <- c((rt*60)-delta,(rt*60)+delta)
  
  mzr<-c( (615.49522-0.08),(615.49522+0.08))
  
  c9 <- chromatogram(CID_9, mz = mzr, rt = rtr)
  
  rt<-rtime(c9[1,1])
  int<-intensity(c9[1,1])
  
  
  x<-rt[!is.na(int)]
  y<-int[!is.na(int)]
  smoothed <- smooth.spline(x,y,spar=0.5)
  
  rt<-smoothed$x
  
  int<-smoothed$y
  
  int <- int/max(int,na.rm=T)*100
  
  

  
  
  
  
  rt_c9<-rt
  int_c9<-int
  
  
  # EAD
  rt <- 4.65
  delta <- 10
  rtr <- c((rt*60)-delta,(rt*60)+delta)
  
  mzr<-c( (615.49522-0.08),(615.49522+0.08))
  
  e9 <- chromatogram(EAD_9, mz = mzr, rt = rtr)
  
  rt<-rtime(e9[1,1])
  int<-intensity(e9[1,1])
  
  
  x<-rt[!is.na(int)]
  y<-int[!is.na(int)]
  smoothed <- smooth.spline(x,y,spar=0.5)
  
  rt<-smoothed$x
  
  int<-smoothed$y
  
  int <- int/max(int,na.rm=T)*100
  
  
  
  
  rt_e9<-rt
  int_e9<-int
  
  
  #####################################################################################
  
  # DG 16:1(9Z)/18:1(12)
  
  # CID
  rt <- 10
  delta <- 12
  rtr <- c((rt*60)-delta,(rt*60)+delta)
  
  mzr<-c( (615.49522-0.08),(615.49522+0.08))
  
  c12 <- chromatogram(CID_12, mz = mzr, rt = rtr)
  
  rt<-rtime(c12[1,1])
  int<-intensity(c12[1,1])
  
  
  x<-rt[!is.na(int)]
  y<-int[!is.na(int)]
  smoothed <- smooth.spline(x,y,spar=0.8)
  
  rt<-smoothed$x
  
  int<-smoothed$y
  
  int <- int/max(int,na.rm=T)*100
  

  rt_c12<-rt
  int_c12<-int
  
  
  # EAD
  rt <- 10
  delta <- 12
  rtr <- c((rt*60)-delta,(rt*60)+delta)
  
  mzr<-c( (615.49522-0.08),(615.49522+0.08))
  
  e12 <- chromatogram(EAD_12, mz = mzr, rt = rtr)
  
  rt<-rtime(e12[1,1])
  int<-intensity(e12[1,1])
  
  
  x<-rt[!is.na(int)]
  y<-int[!is.na(int)]
  smoothed <- smooth.spline(x,y,spar=0.8)
  
  rt<-smoothed$x
  
  int<-smoothed$y
  
  int <- int/max(int,na.rm=T)*100
  

  
  rt_e12<-rt
  int_e12<-int
  
  ##############################################################################
  
  
  # DG 16:1(9Z)/18:1(15)
  
  # CID
  rt <- 6.25
  delta <- 12
  rtr <- c((rt*60)-delta,(rt*60)+delta)
  
  mzr<-c( (615.49522-0.08),(615.49522+0.08))
  
  c15 <- chromatogram(CID_15, mz = mzr, rt = rtr)
  
  rt<-rtime(c15[1,1])
  int<-intensity(c15[1,1])
  
  
  
  x<-rt[!is.na(int)]
  y<-int[!is.na(int)]
  smoothed <- smooth.spline(x,y,spar=0.5)
  
  rt<-smoothed$x
  
  int<-smoothed$y
  
  int <- int/max(int,na.rm=T)*100

  
  rt_c15<-rt
  int_c15<-int
  
  
  # EAD
  rt <- 6.25
  delta <- 12
  rtr <- c((rt*60)-delta,(rt*60)+delta)
  
  mzr<-c( (615.49522-0.08),(615.49522+0.08))
  
  e15 <- chromatogram(EAD_15, mz = mzr, rt = rtr)
  
  rt<-rtime(e15[1,1])
  int<-intensity(e15[1,1])
  
  
  x<-rt[!is.na(int)]
  y<-int[!is.na(int)]
  smoothed <- smooth.spline(x,y,spar=0.4)
  
  rt<-smoothed$x
  
  int<-smoothed$y
  
  int <- int/max(int,na.rm=T)*100

  
  rt_e15<-rt
  int_e15<-int
  
  
  ##########################################################################################
  # CID
  par(yaxs="i")
  par(mar=c(5,5,5,5),cex.lab=1.4,cex.axis=1,cex.main=1.5)
  plot(rt_c9,int_c9,col="coral1",cex=2,pch=16,ylab="Rel. Intensity",xlab="Retention Time [s]",
       main="XID DG 16:1/18:1 CID data",type="l",lwd=2,xlim=c(0,900),ylim=c(0,110))
  points(rt_c12,int_c12,col="coral1",lwd=2,type="l")
  points(rt_c15,int_c15,col="coral1",lwd=2,type="l")
  
  x1<-5
  y1<-3
  
  x2<-rt_c9[1]
  y2<-int_c9[1]
  
  
  stoch_line(x1,y1,x2,y2,color="coral1")
  
  
  x1<-rt_c9[length(rt_c9)]
  y1<-int_c9[length(rt_c9)]
  
  x2<-rt_c15[1]
  y2<-int_c15[1]
  
  
  stoch_line(x1,y1,x2,y2,color="coral1",N=20)
  
  
  x1<-rt_c15[length(rt_c15)]
  y1<-int_c15[length(rt_c15)]
  
  x2<-rt_c12[1]
  y2<-int_c12[1]
  
  
  stoch_line(x1,y1,x2,y2,color="coral1")
  
  
  x1<-rt_c12[length(rt_c12)]
  y1<-int_c12[length(rt_c12)]
  
  x2<-897
  y2<-6
  
  
  stoch_line(x1,y1,x2,y2,color="coral1")
  
  p<-recordPlot()
  
  Output[[1]]<-p
  plot.new()
  # EAD
  par(yaxs="i")
  plot(rt_e9,int_e9,col="palegreen3",cex=2,pch=16,ylab="Intensity",xlab="Retention Time [s]",
       main="XID DG 16:1/18:1 CID data",type="l",lwd=2,xlim=c(0,900),ylim=c(0,110))
  points(rt_e12,int_e12,col="palegreen3",lwd=2,type="l")
  points(rt_e15,int_e15,col="palegreen3",lwd=2,type="l")
  
  x1<-5
  y1<-3
  
  x2<-rt_e9[1]
  y2<-int_e9[1]
  
  
  stoch_line(x1,y1,x2,y2,color="palegreen3")
  
  
  x1<-rt_e9[length(rt_e9)]
  y1<-int_e9[length(rt_e9)]
  
  x2<-rt_e15[1]
  y2<-int_e15[1]
  
  
  stoch_line(x1,y1,x2,y2,color="palegreen3",N=20)
  
  
  x1<-rt_e15[length(rt_e15)]
  y1<-int_e15[length(rt_e15)]
  
  x2<-rt_e12[1]
  y2<-int_e12[1]
  
  
  stoch_line(x1,y1,x2,y2,color="palegreen3")
  
  
  x1<-rt_e12[length(rt_e12)]
  y1<-int_e12[length(rt_e12)]
  
  x2<-897
  y2<-6
  
  
  stoch_line(x1,y1,x2,y2,color="palegreen3")
  
  p<-recordPlot()
  Output[[2]]<-p
  plot.new()
  
  return(Output)
}


chrom<-function(CID,EAD){
  Output<-vector(mode="list",2)
  setwd("C:/FUNCTION_R_DATA")
  
  CID_WT <- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_CID.mzML",mode = "onDisk") 
  
  EAD_WT<- readMSData(files = "22_12_22_Exp_004_SP_Sam_WT_3_EAD.mzML",mode = "onDisk")
  
  WT_C <- chromatogram(CID_WT, aggregationFun = "max")
  
  WT_E <- chromatogram(EAD_WT, aggregationFun = "max")
  
  rt_C<-rtime(WT_C[1,1])
  int_C<-intensity(WT_C[1,1])
  int_C<-int_C/max(int_C,na.rm = T)*700
  
  rt_E<-rtime(WT_E[1,1])
  int_E<-intensity(WT_E[1,1])
  int_E <- int_E/max(int_E,na.rm = T)*700
  
  
  
  #####
  color<-c("azure4","blue4","red1","cyan3","chartreuse2","deeppink","cyan4",
           "chartreuse4","orange4","darkorchid3","brown4",
           "darkorange","darkmagenta","yellow3","deeppink4")
  
  
  
  D_CID <- data.frame(rt=CID$`Average Rt(min)`,
                      mz=CID$`Average Mz`,
                      class=CID$Ontology,
                      col=NA)
  
  D_EAD <- data.frame(rt=EAD$`Average Rt(min)`,
                      mz=EAD$`Average Mz`,
                      class=EAD$Ontology,
                      col=NA)
  
  class<-unique(c(CID$Ontology,EAD$Ontology))
  
  
  for(i in 1:length(class)){
    name<-class[i]
    index_c<- which(D_CID$class==name)
    index_e<-which(D_EAD$class==name )
    D_CID$col[index_c] <- color[i]
    D_EAD$col[index_e] <-color[i]
  }
  
  
  layout(matrix(c(1,2),ncol=2),width=c(8,2))
  par(yaxs="r")
  par(mar=c(5,5,5,0),cex.lab=1.4,cex.axis=1,cex.main=1.5)
  plot(rt_C,int_C,col="coral1",cex=2,pch=16,ylab="mz",xlab="Retention Time [s]",
       main="CID data: Total Ion Chromatogram",type="l",lwd=2,xlim=c(0,900),ylim=c(0,1000))
  points(as.numeric(D_CID$rt)*60,D_CID$mz,pch=16,cex=1.3,col=D_CID$col)
  par(mar=c(5,0,5,0))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("topleft",legend=class,col=color[1:length(class)],pch=16,pt.cex=2,cex=1.4,bty="n")
  
  
  p<-recordPlot()
  plot.new()
  Output[[1]]<-p
  
  
  
  layout(matrix(c(1,2),ncol=2),width=c(8,2))
  par(yaxs="r")
  par(mar=c(5,5,5,0),cex.lab=1.4,cex.axis=1,cex.main=1.5)
  plot(rt_E,int_E,col="palegreen3",cex=2,pch=16,ylab="mz",xlab="Retention Time [s]",
       main="EAD data: Total Ion Chromatogram",type="l",lwd=2,xlim=c(0,900),ylim=c(0,1000))
  points(as.numeric(D_EAD$rt)*60,D_EAD$mz,pch=16,cex=1.3,col=D_EAD$col)
  par(mar=c(5,0,5,0))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("topleft",legend=class,col=color[1:length(class)],pch=16,pt.cex=2,cex=1.4,bty="n")
  
  p<-recordPlot()
  plot.new()
  Output[[2]]<-p
  names(Output)<-c("CID","EAD")
  
  return(Output)
  
}



growth <- function(type=T){
  setwd("W:/users/Abraham/od measurement")
  conditions <- c("Wilt Type",
                  expression(Delta*pgpB),
                  expression(Delta*cdh),
                  expression(Delta*cfa),
                  expression(Delta*plsX),
                  expression(Delta*opgB),
                  expression(Delta*clsC),
                  expression(Delta*pgpC),
                  expression(Delta*pldC),
                  expression(Delta*clsB),
                  expression(Delta*clsA),
                  expression(Delta*pldB),
                  expression(Delta*pldA),
                  expression(Delta*pgpA),
                  expression(Delta*aas),
                  expression(Delta*dgkA)
  )
  
  coul<-c("black","azure4","blue4","red1","cyan3","chartreuse2","deeppink","cyan4",
          "chartreuse4","orange4","darkorchid3","brown4",
          "darkorange","darkmagenta","yellow3","deeppink4")
  
  
  ############################################################################
  # Exponential and Stationary Growth phase:
  ############################################################################
  data<-read.csv("Summary_SP2_timeline.csv",row.names = 1,header=F)
  data<-as.matrix(data)
  data<-data[,-12]
  x<-data[1,]
  data<-data[-1,]
  y<-data[1,]
  
  a<-x
  b<-y
  
  model1 <- nls(b ~ t/ ( 1 +   (( (t-b[1])/b[1])*exp(-r*a)) )   ,start = list( t = 20,r = 0.1))
  
  new_a <- seq(x[1],x[length(x)]+100,length=3000) # creates new sequence of x
  
  
  new_b <- predict(model1,list(a=new_a)) # predicted seuqnce of y 
  
  
  
  
  par(mfrow=c(1,1),mar=c(5,5,3,0))
  layout(matrix(c(1,2),ncol=2),widths = c(8,2))
  plot(x,y,ylab="OD",xlab="Time [min]",col="black",pch=16,cex=1.6,ylim=c(0,3),main="Growth: Individual Mutants",cex.main=1.5,cex.lab=1.6,cex.axis=1.5)
  lines(new_a,new_b,lty = c("79") ,col="black",lwd=3,lend="butt") # add line
  points(a,b,cex=1.6,pch=16)
  
  points(x,y,pch=16,col="black")
  
  
  
  for(i in 2:16){
    
    y<-data[i,]
    
    
    model1 <- nls(y ~ a/ ( 1 +   (( (a-y[1])/y[1])*exp(-r*x)) )   ,start = list( a = 20,r = 0.01))
    
    new_x <- seq(x[1],x[length(x)]+100,length=3000) # creates new sequence of x
    
    
    new_y <- predict(model1,list(x=new_x)) # predicted seuqnce of y 
    
    lines(new_x,new_y,lty = c("79"),col=coul[i],lwd=2,lend="butt") # add line
    
    points(x,y,col=coul[i],pch=16,cex=1.6)
    
  }
  
  points(a,b,cex=1.6,pch=16)
  
  
  par(mar=c(0,0,3,0))
  plot(x,y, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  
  legend("topleft",pch=16,lty=1,col=c(coul),legend=conditions,bty="n",cex=1.5)
  
  
}





##################################################################################

# ms2 library match

ms2<-function(type=T){
  library("grid")
  library("cowplot")
  setwd("W:/users/Abraham/Exp_004_version_2/Data/msp file standards/Experiment")
  files <- list.files(pattern = "msp")
  
  Exp<-vector(mode="list",8) 
  names(Exp)<-gsub("\\..*","",files)
  
  for(i in 1:length(files)){
    
    library <- read.delim(files[i], header=FALSE, comment.char="#",sep="\t", dec=",", quote="")
    
    
    library1 <- library[11:nrow(library),1]
    
    library2 <- matrix(library1,ncol=2,byrow=T)
    
    
    Exp[[i]] <- library2 
    
    
  }
  
  
  
  
  setwd("W:/users/Abraham/Exp_004_version_2/Data/msp file standards/Library")
  files <- list.files(pattern = "txt")
  
  Lib<-vector(mode="list",8) 
  names(Lib)<-gsub("\\..*","",files)
  
  for(i in 1:length(files)){
    
    lib <- read.delim(files[i], header=FALSE, comment.char="#",sep="\t", dec=",", quote="")
    
    
    lib1 <- lib[13:nrow(lib),1]
    
    library2 <- matrix(lib1,ncol=2,byrow=T)
    
    
    Lib[[i]] <- library2 
    
    
  }
  
  
  
  title<-c("DG 15:0-18:1(d7)",
           "LPC 18:1(d7)",
           "LPE 18:1(d7)",
           "PC 15:0-18:1(d7)",
           "PE 15:0-18:1(d7)",
           "PG 15:0-18:1(d7)",
           "PS 15:0-18:1(d7)",
           "TG 15:0-18:1(d7)-15:0")
  
  
  
  test<-matrix(ncol=2,nrow=9)
  
  test[1,]<-1
  test[2,]<-c(2,3)
  test[3,]<-c(10,11)
  test[4,]<-c(4,5)
  test[5,]<-c(12,13)
  test[6,]<-c(6,7)
  test[7,]<-c(14,15)
  test[8,]<-c(8,9)
  test[9,]<-c(16,17)
  layout(test,heights = c(1,8,8,8,8,8,8,8,8))
  par(mar=c(0,2,0,2))
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')  
  mtext(3,text="Standard Library Spectral Match",line=-2,cex=1.5)
  
  
  
  
  
  for(i in 1:8){
    d2 <- Exp[[i]]
    d4 <- Lib[[i]]
    
    
    x2 <- as.numeric(d2[,1])
    y2 <- as.numeric(d2[,2])/max(as.numeric(d2[,2]))*100
    
    
    x4 <- as.numeric(d4[,1])
    y4 <- as.numeric(d4[,2])/max(as.numeric(d4[,2]))*100
    
    
    xmin=min(as.numeric(c(x2,x4)))
    
    xmax=max(as.numeric(c(x2,x4)))
    
    y2max=max(y2)
    
    y4max=max(y4)
    
    
    par(mar=c(0,4,5,4),cex.lab=1,cex.axis=1,cex.main=1)
    
    plot(x2,y2,type= "h",xaxt="n",xlab="",ylab="",main=paste0(title[i]," Standard Lipid"  ),col="red",lwd=3,xaxs = "i",
         yaxs = "i",xlim=c(xmin,xmax+10),ylim=c(0,150),yaxt="n")
    mtext(3,text="Experimental Spectra",adj=0,col="red",line=-1.2)
    
    
    
    
  }
  
  
  
  
  
  for(i in 1:8){
    d2 <- Exp[[i]]
    d4 <- Lib[[i]]
    
    
    x2 <- as.numeric(d2[,1])
    y2 <- as.numeric(d2[,2])/max(as.numeric(d2[,2]))*100
    
    
    x4 <- as.numeric(d4[,1])
    y4 <- as.numeric(d4[,2])/max(as.numeric(d4[,2]))*100
    
    
    xmin=min(as.numeric(c(x2,x4)))
    
    xmax=max(as.numeric(c(x2,x4)))
    
    y2max=max(y2)
    
    y4max=max(y4)
    
    
    
    par(mar=c(5,4,0,4))
    plot(x4,y4,type= "h",xlab="m/z",ylab="",col="blue",lwd=3,xaxs = "i",
         yaxs = "i",xlim=c(xmin,xmax+10) ,ylim=c(150,0),yaxt="n")  
    mtext(1,text="Reference Spectra",adj=0,col="blue",line=-1.2)
    
    
    
    
  }
  p<- recordPlot()
  par(mar=c(0,4,5,4),mfrow=c(1,1))
  plot.new()
  
  return(p)
  
  
  
  
}



Overlap_enriched_per_ko<-function(EAD,CID,phase,fragmentation,Directory){
  
  conditions <- c("Wilt Type",
                  expression(Delta*pgpB),
                  expression(Delta*cdh),
                  expression(Delta*cfa),
                  expression(Delta*plsX),
                  expression(Delta*opgB),
                  expression(Delta*clsC),
                  expression(Delta*pgpC),
                  expression(Delta*pldC),
                  expression(Delta*clsB),
                  expression(Delta*clsA),
                  expression(Delta*pldB),
                  expression(Delta*pldA),
                  expression(Delta*pgpA),
                  expression(Delta*aas),
                  expression(Delta*dgkA)
                  
  )
  if(phase=="EP"){
    s<-c(c(36:80),c(84:86))
    name <- rep(c(paste0("EP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }else if(phase=="SP"){
    s<-c(c(87:131),c(135:137))
    name <- rep(c(paste0("SP_",c("WT",c(1:15)))),each=3)
    samples <- name[order(name)]
  }
  sequence<-seq(1,45,3)
  
  
  cc_C<-characteristics_CID(CID)
  cc_E<-characteristics_EAD(EAD)
  
  ######################################################################
  
  
  Table_C<-matrix(ncol=15,nrow=nrow(CID))
  Table_E<-matrix(ncol=15,nrow=nrow(EAD))
  
  
  
  
  dc<-CID[,s]
  de<-EAD[,s]
  i=1
  for(i in 1:15){
    name<-conditions[i+1]
    
    wt<-c(46,47,48)
    k<-sequence[i]
    ko<-c(k,k+1,k+2)
    dataframe_C <- dc[,c(wt,ko)]
    dataframeC<-matrix(as.numeric(as.matrix(dataframe_C)),ncol=6)
    colnames(dataframeC)<-colnames(dataframe_C)
    dataframe_C<-dataframeC
    
    
    w<-rowMeans(dataframe_C[,c(1:3)])
    k<-rowMeans(dataframe_C[,c(4:6)])
    
    FC<-k/w
    
    FC_log_C<-log2(FC)
    
    FC_log_C[FC_log_C==-Inf]<- min(FC_log_C)-1
    FC_log_C[FC_log_C==Inf]<- max(FC_log_C)+1
    
    
    
    pvalue <- rep(NA,nrow(dataframe_C))  
    
    for(j in 1:nrow(dataframe_C)){
      ww<-dataframe_C[j,c(1:3)]
      kk<-dataframe_C[j,c(4:6)]
      
      t<-t.test(ww,kk,alternative = "two.sided",conf.level = 0.95)
      value<- (t$p.value)
      
      pvalue[j]<-value
    }
    
    index <- which(abs(FC_log_C)>=1 & pvalue <= 0.05)
    
    Table_C[index,i] <- cc_C$`Species level`[index]
    
    dataframe_E <- de[,c(wt,ko)]
    dataframeE<-matrix(as.numeric(as.matrix(dataframe_E)),ncol=6)
    colnames(dataframeE)<-colnames(dataframe_E)
    dataframe_E<-dataframeE
    
    
    w<-rowMeans(dataframe_E[,c(1:3)])
    k<-rowMeans(dataframe_E[,c(4:6)])
    
    FC<-k/w
    
    FC_log_E<-log2(FC)
    
    FC_log_E[FC_log_E==-Inf]<- min(FC_log_E)-1
    FC_log_E[FC_log_E==Inf]<- max(FC_log_E)+1
    
    pvalue <- rep(NA,nrow(dataframe_E))  
    
    for(j in 1:nrow(dataframe_E)){
      ww<-dataframe_E[j,c(1:3)]
      kk<-dataframe_E[j,c(4:6)]
      
      t<-t.test(ww,kk,alternative = "two.sided",conf.level = 0.95)
      value<- (t$p.value)
      
      pvalue[j]<-value
    }
    
    index <- which(abs(FC_log_E)>=1 & pvalue <= 0.05)
    
    Table_E[index,i] <- cc_E$`Species level`[index]
    
    
    
    
    
    
    
    
  }
  overlap<-matrix(ncol=15,nrow=3)
  for(i in 1:15){
    c<-Table_C[,i]
    e<-Table_E[,i]
    
    c<-na.omit(c)
    e<-na.omit(e)
    overlap[3,i]<-length(which(c %in% e==T))
    overlap[2,i]<-length(which(e %in% c==F))
    
    overlap[1,i]<-length(which(c %in% e==F))
    
  }
  
  
  ov_norm<-overlap
  for(i in 1:15){
    ov_norm[,i] <- ov_norm[,i]/sum(ov_norm[,i])*100
  }
  layout(matrix(c(1,2,3),ncol=3),width=c(8,8,2))
  par(mar=c(6,5,5,2),cex.lab=1.7,cex.axis=1.5,cex.main=1.3)
  barplot(overlap/1.5,col=c("coral1","palegreen3","grey"),names.arg = conditions[2:16],las=2,
          ylab = "# of Enriched Lipids",main="",cex.names=1.7)
  mtext(side=3,"Absolute Values",cex=1.4,line=0.1)
  x<-barplot(ov_norm,col=c("coral1","palegreen3","grey"),names.arg = conditions[2:16],las=2,
             ylab = "Percentage of Enriched Lipids",main="",cex.names=1.7)
  mtext(side=3,"Percentage",cex=1.4,line=0.1)
  
  text(x,ov_norm[1,]-7,cex=1.2,las=1,paste0(round(ov_norm[1,],digits = 0),"%"),srt = 90)
  
  text(x[-3],colSums(ov_norm[1:2,-3])-c(rep(4.7,13),3.1),cex=1.2,las=1,paste0(round((ov_norm[2,-3]),digits = 0),"%"),srt = 90)
  
  
  text(x,colSums(ov_norm[1:3,])-15,cex=1.2,las=1,paste0(round((ov_norm[3,]),digits = 0),"%"),srt = 90)
  
  
  mtext(side=3,text="Enriched Lipids Overlap across Conditions",outer = T,line=-1.8,cex=1.5)
  
  
  par(mar=c(0,0,5,0))
  plot(1:10,1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')  
  legend("topleft",col=c("coral1","palegreen3","grey"),pch=15,legend=c("CID","EAD","Overlap"),bty="n",pt.cex = 2.6,cex=2)
  #####################################################################################
  
}





