setwd("~/Documentos/Omar/Tesis/Taxa/Results/Julio1/Endemic/")

end.table <- read.csv("Grid.end")

head(end.table)

############################################################
# AvTD
#
which(end.table$TD==0)

kruskal.test(end.table$TD[-which(end.table$TD==0)],as.factor(end.table$Percentage[-which(end.table$TD==0)]))

############################################################
# AvTD
#
which(end.table$PD==0)

kruskal.test(end.table$PD[-which(end.table$PD==0)],as.factor(end.table$Percentage[-which(end.table$PD==0)]))

############################################################
# AvTD
#
which(end.table$AvTD==0)

kruskal.test(end.table$AvTD[-which(end.table$AvTD==0)],as.factor(end.table$Percentage[-which(end.table$AvTD==0)]))
