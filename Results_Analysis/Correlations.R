library(ggplot2)

## Set working directory
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/May18/")

## Read the absence/presence tables for the total richness

area.rich <- read.csv("Area.dist.matrix")
grid.rich <- read.csv("grid_25g.dist.matrix")

## Create two vectors where the total richness will are

area.t.rich <- c()
grid.t.rich <- c()

## Make the summatory of the richness (This could be made by a apply function, check it later)

for(i in 2:length(colnames(area.rich))){
  total <- sum(area.rich[,i])
  area.t.rich[i] <- total
}

for(i in 2:length(colnames(grid.rich))){
  total <- sum(grid.rich[,i])
  grid.t.rich[i] <- total
}


area.t.rich <- area.t.rich[-1]

grid.t.rich <- grid.t.rich[-1]

# Create the two data frames with the index values and richness values
# f.Rich are the the filogenetic richness
# t.Rich are the whole richness

setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/May18/Raw_IndexR/")

index.area <- data.frame(Area=read.csv("DT.Area")$area,
                         DT=read.csv("DT.Area")$W,
                         PD=read.csv("PD.Area")$PD,
                         AvDT=read.csv("AvTD.Area")$Dplus,
                         f.Rich=read.csv("DT.Area")$rich,
                         t.Rich=area.t.rich)

index.grid <- data.frame(Cell.Grid=read.csv("DT.grid25")$area,
                         DT=read.csv("DT.grid25")$W,
                         PD=read.csv("PD.grid25")$PD,
                         AvDT=read.csv("AvTD.grid25")$Dplus,
                         f.Rich=read.csv("DT.grid25")$rich,
                         t.Rich=grid.t.rich)

## Check first f.Rich vs t.Rich

f.t  <- ggplot()
f.t <- f.t + geom_point(aes(x=index.grid$t.Rich,y=index.grid$f.Rich))
f.t

## f.Rich vs Index

fi <- ggplot()
fi <- fi + geom_point(aes(x=index.grid$f.Rich,y = log(index.grid$DT)),colour="red",shape=1)
fi <- fi + geom_point(aes(x=index.grid$f.Rich,y = log(index.grid$PD)),colour="green",shape=1)
fi <- fi + geom_point(aes(x=index.grid$f.Rich,y = index.grid$AvDT),colour="orange",shape=1)
fi

## t.Rich vs Index

fi <- ggplot()
fi <- fi + geom_point(aes(x=index.grid$t.Rich,y = log(index.grid$DT)),colour="red",shape=2)
fi <- fi + geom_point(aes(x=index.grid$t.Rich,y = log(index.grid$PD)),colour="green",shape=2)
fi <- fi + geom_point(aes(x=index.grid$t.Rich,y = index.grid$AvDT),colour="orange",shape=2)
fi

# t.Rich+f.rich vs Index

fi <- ggplot()
fi <- fi + geom_point(aes(x=index.grid$f.Rich,y = log(index.grid$DT)),colour="red",shape=1)
fi <- fi + geom_point(aes(x=index.grid$f.Rich,y = log(index.grid$PD)),colour="green",shape=1)
fi <- fi + geom_point(aes(x=index.grid$f.Rich,y = index.grid$AvDT),colour="orange",shape=1)
fi <- fi + geom_point(aes(x=index.grid$t.Rich,y = log(index.grid$DT)),colour="red",shape=2)
fi <- fi + geom_point(aes(x=index.grid$t.Rich,y = log(index.grid$PD)),colour="green",shape=2)
fi <- fi + geom_point(aes(x=index.grid$t.Rich,y = index.grid$AvDT),colour="orange",shape=2)
fi

## Create outcome matrix

correlations <- matrix(NA, nrow=2,ncol=3)
rownames(correlations) <- c("Total.Rich","Filo.Rich")
colnames(correlations) <- c("DT", "PD", "AVDT")

## Correlations t.rich
##DT vs Total Richness
dt <- cor.test(index.grid$t.Rich,index.grid$DT,method = "spearman")
correlations[1,1] <- dt$estimate
## PD vs Total Richness
pd <- cor.test(index.grid$t.Rich,index.grid$PD,method = "spearman")
correlations[1,2] <- pd$estimate
## AVDT vs Total Richness
avdt <- cor.test(index.grid$t.Rich,index.grid$AvDT,method = "spearman")
correlations[1,3] <- avdt$estimate

## Correlations f.rich
## DT vs Filogenetic Richness
dt <- cor.test(index.grid$f.Rich,index.grid$DT,method = "spearman")
correlations[2,1] <- dt$estimate
## PD vs Filogenetic Richness
pd <- cor.test(index.grid$f.Rich,index.grid$PD,method = "spearman")
correlations[2,2] <- pd$estimate
## AVDT vs Filogenetic Richness
avdt <- cor.test(index.grid$f.Rich,index.grid$AvDT,method = "spearman")
correlations[2,3] <- avdt$estimate

correlations

write.table(correlations,file="~/Documentos/Omar/Tesis/Taxa/Results/May18/Corr.RichvIndex",
            quote = F, row.names = T,col.names = T)

ind.cor <- cor(index.grid[,2:4],method = "spearman")

write.table(ind.cor,file = "~/Documentos/Omar/Tesis/Taxa/Results/May18/Corr.IndvInd", 
            quote = F, row.names = T,col.names = T)
