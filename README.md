# PhyloDiversity
___

Hi, this is the repository for: **Conserving Northern Andean Block: What Phylogenies Say - Leon-Alvarado & Miranda-Esquivel (Submitted).**

Here, all the data and R-script used for the analyses are stored. Everything here is under the GPL V2 lincense. So, you can use our data if you desire, but don't forget the cite. For a detailed explanaiton about the methods that we used, please carefully read the following document. 

___

**Data Information**

For phylogenetic diversity index (or at least with the approcha that we implemented) is just neccesary two types of data:

  1. Species distributions (occurrences)
  2. Phylogenies
  
  
First, you must choose your taxon or taxa, so you can use your selection rules, or if you desire you can follow our [rules](https://github.com/oleon12/PhyloDiversity/blob/master/Img/Diagrama_Flujo.png).

Once you had been selected your taxon or taxa you sould create the data for the distribution and phylogenies. Basically for phylogenies you can use the algorithm of predilection, but remeber that only TD index do not take into account branch lengths. Our code was build for trees in [Newick format](http://evolution.genetics.washington.edu/phylip/newicktree.html), so we stongly recommended to fit your trees in this format.

___

**Distributions**


1. First, we need to download the occurrences for each taxon, here we use the [Occ_Download.R](https://github.com/oleon12/PhyloDiversity/blob/master/R/Occ_Download.R) script. For this step you need the polygon of the study area (Northern Andean Block here) and the name of the taxon. You can use any taxonomic level in the search (Orden, Family, Genera or Species), but given our experience we recommended use only Genera or Species level. At the end, the results is a CSV file with three columns: Species, Longitud and Latitude.

> *_Please, be careful with the append parameter in the last write.table command, it only works properly if you add informtion over a existing file, but if you create a file with this paramater set TRUE will generated a Warning message_

2. Now, with the occurrences downloades it is time to create the absence-presence matrix. This matrix is just a table with Species as rows and Areas as columns, where, if a species is present in an area, the intersection will be filled with 1, and 0 if the species is absent. For this task, we use the [Absence.Presence_Matrix.R](https://github.com/oleon12/PhyloDiversity/blob/master/R/Absence.Presence_Matrix.R) script.
