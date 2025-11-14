library(tidyverse)
library(scales)
library(caret)
library(pROC)
library(ROCR)

Clustering=function(Seq,Res) {
  CL=0 #Cluster Length
  CC=0 #Current Cluster
  for (i in 1:nchar(Seq)) {
    if (i==nchar(Seq) & grepl(Res,substr(Seq,i,i))) {
      CC=CC+1
      if (CC>CL) {
        CL=CC
      }
    } else if (grepl(Res,substr(Seq,i,i))) {
      CC=CC+1
    } else if (!grepl(Res,substr(Seq,i,i)) & CC>CL) {
      CL=CC
      CC=0
    } else {
      CC=0
    }
  }
  CL
}

CalcRules=function(Seqs,
                   Slopes,
                   BinaryFunction=NULL,
                   RedundancyLimit=2,
                   ChargeLimit=-2,
                   BalanceDevLimit=3,
                   PositionLimit=12,
                   ClusterRatioLimit=0.5,
                   ClusterRedundancy=3){
  ##Parameters
  #Seqs is a vector containing sequences of a consistent length
  #Slopes is a vector of the same length as Seqs where 0 is a cutoff for functional sequences
  #BinaryFunction is a vector of the same length as Seqs where 1 is functional 0 is nonfunctional
  #Redundancy rule followed if AroCount and AcidicCount are greater than RedundancyLimit
  #Charge rule followed if net charge is less than ChargeLimit
  #Balance rule followed if abs(AroCount-AcidicCount) is less than BalanceDevLimit
  #Position rule followed if average position from C-terminus is less than PositionLimit
  #Cluster rule followed if the cluster neighbor ratio is less than ClusterRatioLimit
  
  ##Initial Sequence Features and Function
  if (length(BinaryFunction)==0) {
    BinaryFunction=ifelse(Slopes>0,1,0)
  }
  MainDF=data.frame(Sequence=Seqs,
                    Slope=Slopes,
                    SlopeBinary=factor(BinaryFunction))
  for (i in c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")) {
    MainDF[[i]]=nchar(Seqs)-nchar(gsub(paste("[",i,"]",sep = ""),"",Seqs))
  }
  MainDF$AroCount=nchar(Seqs)-nchar(gsub("[WF]","",Seqs)) #Previously "[WFY]"
  MainDF$AcidicCount=nchar(Seqs)-nchar(gsub("[DE]","",Seqs))
  MainDF$BasicCount=nchar(Seqs)-nchar(gsub("[RK]","",Seqs))
  MainDF$SeqLength=nchar(Seqs)
  
  ##Calculate Additional Features
  #Redundancy
  Redundancy=c()
  for (i in 1:nrow(MainDF)) {
    Redundancy=c(Redundancy,min(c(MainDF[i,]$AroCount,MainDF[i,]$AcidicCount)))
  }
  MainDF$Redundancy=Redundancy #The smaller of AroCount and AcidicCount
  
  #Charge
  MainDF$Charge=MainDF$BasicCount-MainDF$AcidicCount
  
  #Balance
  MainDF$BalDev=abs(MainDF$AroCount-MainDF$AcidicCount)  
  
  #Position
  AvgPos=c()
  for (Seq in MainDF$Sequence) {
    SeqVec=str_split(Seq,"")%>%unlist
    i=0
    j=0
    Temp=c()
    for (position in SeqVec) {
      i=i+1
      if (grepl("[WFY]",position)) {
        Temp=c(Temp,i)
        j=j+1
      }
    }
    if (j==0) {
      Temp=0
    }
    AvgPos=c(AvgPos,round(mean(Temp)))
  }
  MainDF$AvgPosition=AvgPos
  MainDF$PositionDev=MainDF$SeqLength-MainDF$AvgPosition
  
  #Clustering
  NeighborRatio=c()
  for (Seq in MainDF$Sequence) {
    Wgroups=unlist(strsplit(gsub("[^WF]","X",Seq),"X+")) #Previously [^WFY]
    Temp=c()
    for (group in 1:length(Wgroups)) {
      if (nchar(Wgroups[group])==1) {
        Temp=c(Temp,0)
      }
      if (nchar(Wgroups[group])>1) {
        Temp=c(Temp,rep(0.5,2),rep(1,nchar(Wgroups[group])-2))
      }
    }
    if (length(Temp)==0){
      NeighborRatio=c(NeighborRatio,0)
    }
    if (length(Temp)>0){
      NeighborRatio=c(NeighborRatio,mean(Temp))
    }
  }
  MainDF$NeighborRatio=NeighborRatio
  AvgNeighborRatio=1-mean(MainDF[MainDF$Redundancy>=ClusterRedundancy,]$NeighborRatio)
  MainDF$AroCluster=sapply(MainDF$Sequence, Clustering, Res="[WF]") #Integer of length of largest cluster of aromatics
  LeftFlankAcidic=c()
  RightFlankAcidic=c()
  IsFlanked=c()
  for (i in 1:nrow(MainDF)) {
    MaskedSeq=gsub("[WF]","X",MainDF[i,]$Sequence)
    if (grepl("X",MaskedSeq)) {
      Flanks=unlist(strsplit(MaskedSeq,paste("X{",MainDF[i,]$AroCluster,"}",sep = "")))
      LeftFlankAcidic=c(LeftFlankAcidic,grepl("[DE]",str_sub(Flanks[1],-1,-1)))
      RightFlankAcidic=c(RightFlankAcidic,grepl("[DE]",str_sub(Flanks[2],1,1)))
      IsFlanked=c(IsFlanked,ifelse((LeftFlankAcidic[i]+RightFlankAcidic[i])==2,1,0))
    } else {
      LeftFlankAcidic=c(LeftFlankAcidic,0)
      RightFlankAcidic=c(RightFlankAcidic,0)
      IsFlanked=c(IsFlanked,0)
    }
  }
  MainDF$Flanked=IsFlanked #1 for flanked, 0 for not flanked or no aromatics to flank
  
  ##Redundancy
  MainDF$RuleRedundancyOld=with(MainDF,ifelse(Redundancy>2,1,0))
  MainDF$RuleRedundancyCutoff=with(MainDF,ifelse(Redundancy>RedundancyLimit,1,0))
  MainDF$RuleRedundancyNew=rescale(with(MainDF,ifelse(Redundancy>5,1,
                                                      ifelse(Redundancy<3,-3,(Redundancy-5)))),c(0,1),c(-3,1))
  
  ##Charge
  MainDF$RuleChargeOld=with(MainDF,ifelse(BasicCount==0,1,0))
  MainDF$RuleChargeCutoff=with(MainDF,ifelse(Charge<ChargeLimit,1,0))
  MainDF$RuleChargeNew=rescale(with(MainDF,ifelse(Charge<(-5),1,
                                                  ifelse(Charge>(-3),-3,(-5-Charge)))),c(0,1),c(-3,1))
  
  ##Balance
  MainDF$RuleBalanceOld=with(MainDF,ifelse(BalDev==0,1,0))
  MainDF$RuleBalanceCutoff=with(MainDF,ifelse(BalDev<BalanceDevLimit,1,0)) 
  MainDF$RuleBalanceNew=rescale(with(MainDF,ifelse(BalDev<3,1,
                                                   ifelse(BalDev>7,-5,(3-BalDev)))),c(0,1),c(-5,1))
  
  ##Position
  MainDF$RulePositionOld=with(MainDF,ifelse(PositionDev<9,1,0))
  MainDF$RulePositionCutoff=with(MainDF,ifelse(PositionDev<PositionLimit,1,0))
  MainDF$RulePositionNew=rescale(with(MainDF,ifelse(PositionDev<12,1,
                                                    ifelse(PositionDev>15,-4,(12-PositionDev)))),c(0,1),c(-4,1))
  
  ##Cluster (Aromatics WF)
  MainDF$RuleClusterOld=with(MainDF,ifelse(NeighborRatio<0.5,1,0))
  MainDF$RuleClusterCutoff=with(MainDF,ifelse(NeighborRatio<ClusterRatioLimit,1,0))
  MainDF$RuleClusterNew=with(MainDF,ifelse(Flanked==1,1,ifelse(Redundancy<ClusterRedundancy,AvgNeighborRatio,(1-NeighborRatio))))
  
  return(MainDF)
}
