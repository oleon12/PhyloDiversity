setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final/Endemic/")

#################################################################################################
##                                                                                             ##
##                                   Grido of 0.25° Cells                                      ##
##                                                                                             ##
#################################################################################################

############################################################
# TD
#

end<- read.csv("Grid.end")

head(end)

class <- readShapePoly("~/Documentos/Omar/Tesis/Taxa/Results/Final/Raw_IndexR/Grid25_TD.shp")

ind <- class$Index

ind <- rep(ind,(length(end$Area)/length(unique(end$Area))))
           
end <- cbind(end,ind)
     
q5 <- grep("Q5",end$ind)
           
q4 <- grep("Q4",end$ind)

end2 <- end[-which(end$TD==0), ]

kruskal.test(end2$TD[c(q5,q4)],as.factor(end2$Percentage[c(q5,q4)]))

############################################################
# PD
#

end<- read.csv("Grid.end")

head(end)

class <- readShapePoly("~/Documentos/Omar/Tesis/Taxa/Results/Final/Raw_IndexR/Grid25_PD.shp")

ind <- class$Index

ind <- rep(ind,(length(end$Area)/length(unique(end$Area))))

end <- cbind(end,ind)

q5 <- grep("Q5",end$ind)

q4 <- grep("Q4",end$ind)

end2 <- end[-which(end$PD==0), ]

kruskal.test(end2$PD[c(q5,q4)],as.factor(end2$Percentage[c(q5,q4)]))

############################################################
# AvTD
#

end<- read.csv("Grid.end")

head(end)

class <- readShapePoly("~/Documentos/Omar/Tesis/Taxa/Results/Final/Raw_IndexR/Grid25_AvTD.shp")

ind <- class$Index

ind <- rep(ind,(length(end$Area)/length(unique(end$Area))))

end <- cbind(end,ind)

q5 <- grep("Q5",end$ind)

q4 <- grep("Q4",end$ind)

end2 <- end[-which(end$AvTD==0), ]

kruskal.test(end2$AvTD[c(q5,q4)],as.factor(end2$Percentage[c(q5,q4)]))

#################################################################################################
##                                                                                             ##
##                                   Areas of endemism                                         ##
##                                                                                             ##
#################################################################################################

############################################################
# TD
#

end<- read.csv("Area.end")

kruskal.test(end$TD,as.factor(end$Percentage))

############################################################
# PD
#
end<- read.csv("Area.end")

kruskal.test(end$PD,as.factor(end$Percentage))

############################################################
# AvTD
#
end<- read.csv("Area.end")

kruskal.test(end$AvTD,as.factor(end$Percentage))
