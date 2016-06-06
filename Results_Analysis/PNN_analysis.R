setwd("~/Documentos/Omar/Tesis/Taxa/Results/May18/Raw_IndexR/")

PNN <- data.frame(Area=read.csv("DT.PNN")$area,
                  TD=read.csv("DT.PNN")$W,
                  PD=read.csv("PD.PNN")$PD,
                  AvTD=read.csv("AvTD.PNN")$Dplus)
PNN

## Tres calculos
##             /PNNout+PNNin  /PNNout+PNN+ /NABin
## PNNin
## PNN+
## PNNout-NAB


## For TD

PNNin <- PNN$TD[grep("PNNin",PNN$Area)]
PNN1 <- sum(PNN$TD[1:117])
PNNNAB <- PNN$TD[grep("NABin",PNN$Area)] - PNN$TD[grep("PNNout",PNN$Area)] 
NAB <- PNN$TD[grep("NABin",PNN$Area)]

PNN.td <- matrix(0, ncol = 3, nrow = 3)

PNN.td[1,1] <- PNNin/(PNN$TD[grep("PNNout",PNN$Area)]+PNNin)
PNN.td[2,1] <- PNN1/(PNN$TD[grep("PNNout",PNN$Area)]+PNNin)
PNN.td[3,1] <- PNNNAB/(PNN$TD[grep("PNNout",PNN$Area)]+PNNin)

PNN.td[1,2] <- PNNin/(PNN$TD[grep("PNNout",PNN$Area)]+PNN1)
PNN.td[2,2] <- PNN1/(PNN$TD[grep("PNNout",PNN$Area)]+PNN1)
PNN.td[3,2] <- PNNNAB/(PNN$TD[grep("PNNout",PNN$Area)]+PNN1)

PNN.td[1,3] <- PNNin/NAB
PNN.td[2,3] <- PNN1/NAB
PNN.td[3,3] <- PNNNAB/NAB

colnames(PNN.td) <- c("/PNNout+PNNin","/PNNout+PNN+","/NAB")
rownames(PNN.td) <- c("PNNin","PNN+","PNNout-NAB")

PNN.td

## For PD

PNNin <- PNN$PD[grep("PNNin",PNN$Area)]
PNN1 <- sum(PNN$PD[1:117])
PNNNAB <- PNN$PD[grep("NABin",PNN$Area)] - PNN$PD[grep("PNNout",PNN$Area)] 
NAB <- PNN$PD[grep("NABin",PNN$Area)]

PNN.pd <- matrix(0, ncol = 3, nrow = 3)

PNN.pd[1,1] <- PNNin/(PNN$PD[grep("PNNout",PNN$Area)]+PNNin)
PNN.pd[2,1] <- PNN1/(PNN$PD[grep("PNNout",PNN$Area)]+PNNin)
PNN.pd[3,1] <- PNNNAB/(PNN$PD[grep("PNNout",PNN$Area)]+PNNin)

PNN.pd[1,2] <- PNNin/(PNN$PD[grep("PNNout",PNN$Area)]+PNN1)
PNN.pd[2,2] <- PNN1/(PNN$PD[grep("PNNout",PNN$Area)]+PNN1)
PNN.pd[3,2] <- PNNNAB/(PNN$PD[grep("PNNout",PNN$Area)]+PNN1)

PNN.pd[1,3] <- PNNin/NAB
PNN.pd[2,3] <- PNN1/NAB
PNN.pd[3,3] <- PNNNAB/NAB

colnames(PNN.pd) <- c("/PNNout+PNNin","/PNNout+PNN+","/NAB")
rownames(PNN.pd) <- c("PNNin","PNN+","PNNout-NAB")

PNN.pd


## For AvTD

PNNin <- PNN$AvTD[grep("PNNin",PNN$Area)]
PNN1 <- sum(PNN$AvTD[1:117])
PNNNAB <- (PNN$AvTD[grep("NABin",PNN$Area)] - PNN$AvTD[grep("PNNout",PNN$Area)])*-1 
NAB <- PNN$AvTD[grep("NABin",PNN$Area)]

PNN.avtd <- matrix(0, ncol = 3, nrow = 3)

PNN.avtd[1,1] <- PNNin/(PNN$AvTD[grep("PNNout",PNN$Area)]+PNNin)
PNN.avtd[2,1] <- PNN1/(PNN$AvTD[grep("PNNout",PNN$Area)]+PNNin)
PNN.avtd[3,1] <- PNNNAB/(PNN$AvTD[grep("PNNout",PNN$Area)]+PNNin)

PNN.avtd[1,2] <- PNNin/(PNN$AvTD[grep("PNNout",PNN$Area)]+PNN1)
PNN.avtd[2,2] <- PNN1/(PNN$AvTD[grep("PNNout",PNN$Area)]+PNN1)
PNN.avtd[3,2] <- PNNNAB/(PNN$AvTD[grep("PNNout",PNN$Area)]+PNN1)

PNN.avtd[1,3] <- PNNin/NAB
PNN.avtd[2,3] <- PNN1/NAB
PNN.avtd[3,3] <- PNNNAB/NAB

colnames(PNN.avtd) <- c("/PNNout+PNNin","/PNNout+PNN+","/NAB")
rownames(PNN.avtd) <- c("PNNin","PNN+","PNNout-NAB")

PNN.avtd



setwd("~/Documentos/Omar/Tesis/Taxa/Results/May18/Raw_IndexR/")

write.table(PNN.td,
            file = "TD_PNN.r",
            col.names = T, row.names = T, quote = F, sep = ",")

write.table(PNN.pd,
            file = "PD_PNN.r",
            col.names = T, row.names = T, quote = F, sep = ",")

write.table(PNN.avtd,
            file = "AvTD_PNN.r",
            col.names = T, row.names = T, quote = F, sep = ",")
