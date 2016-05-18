library(rgeos)
library(rgbif)
library(maptools)

polygon<-'POLYGON((-80.3311270045821 -3.04137854749849,-81.0029497681885 -2.23389926431788,-80.9900300996576 -1.03237009094502,-80.1502516451497 0.731164663521601,-77.5146392648479 8.93515418063738,-74.8531875474844 11.1573371679506,-71.6620294203543 12.4622236895706,-69.9953921798694 12.2555089930764,-68.4191926191007 11.2477748476669,-67.8765665408032 10.4725947358135,-69.2202120680159 9.36150324215685,-71.6749490888851 7.39771362546142,-71.6232704147616 6.44165815417549,-72.0366998077501 5.64063870526025,-75.4474922999053 1.48050543831336,-77.2045672201065 -1.31014296435909,-78.3156587137631 -2.36955578389215,-80.3311270045821 -3.04137854749849))'

occ <- occ_search(scientificName = 'Viperidae',geometry = polygon,limit = 200000 )


occ.dist <- data.frame(Sp=occ$data$species,
                       Lon=occ$data$decimalLongitude,
                       Lat=occ$data$decimalLatitude)

if(TRUE%in%is.na(occ.dist$Sp)==T){

  occ.dist <- occ.dist[-which(is.na(occ.dist$Sp)==T),]

}

occ.dist

unique(occ.dist$Sp)

setwd("~/Documentos/Omar/Tesis/Taxa/Squamata//Viperidae")

write.table(occ.dist,"Viperidae.occ",
            quote = F,row.names = F,col.names = T,sep = ",")

setwd("~/Documentos/Omar/Tesis/Taxa/Richness/")

write.table(occ.dist,file="Richness.occ",
            quote = F,row.names = F,col.names = F,
            append = T,sep=",")
