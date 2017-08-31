# PhyloDiversity
___

Hi, this is the repository for: **Conserving Northern Andean Block: What Phylogenies Say - Leon-Alvarado & Miranda-Esquivel (Submitted).**

Here, all the data and R-script used for the analyses are stored. Everything here is under the [GPL V2 lincense](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html). So, you can use our data if you desire, but don't forget the cite. For a detailed explanaiton about the methods that we used, please carefully read the following document. 

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

> _Please, be careful with the append parameter in the last write.table command, it only works properly if you add informtion over a existing file, but if you create a file with this paramater set TRUE, it will generate a Warning message_

2. Due to some occurrences could be inside the sea, and are not of interest in this study it is necesary to remove it. For this task we use the [CleanOcc.R](https://github.com/oleon12/PhyloDiversity/blob/master/R/CleanOcc.R) script. This script use a polygon in Shape File format, here we use the polygon of the Northern Andean Block which are just continental land, and from the occurrences downloaded the script found the points outside of the polygon and remove it from the table.

3. Now, with the occurrences downloades it is time to create the absence-presence matrix. This matrix is just a table with Species as rows and Areas as columns, where, if a species is present in an area, the intersection will be filled with 1, and 0 if the species is absent. For this task, we use the [Absence.Presence_Matrix.R](https://github.com/oleon12/PhyloDiversity/blob/master/R/Absence.Presence_Matrix.R) script. To use this script you need the occurences (previously downloaded) and polygons in Shape Files (Areas). In this case we use three different areas: Areas of endemism, Grid cells and Protected Areas. All of these areas are in Shape Files. We strongly recommended that all areas be inside in one file, for example, we use more than 100 Protected Areas, however, all areas are inside in a single Shape File. This recomendation is just for save memory.

___

**Phylogenetic Diversity Index**

Now, we already have the two necesary data (Species Distribution and Phylogenies) for calculate the Phylogenetic Diversity indices. Here we implemented three indices: Taxonomic Distinctness (TD), a topological based index ([Vane-Wright _et al_., 1991](http://www.sciencedirect.com/science/article/pii/000632079190030D)); Phylogenetic Diversity (PD), a minimum spanning tree based index ([Faith, 1992](http://www.sciencedirect.com/science/article/pii/0006320792912013)) and Average Taxonomic Distinctness (AvTD), a pair-wise distance based index ([Clarke & Warwick, 1998](http://onlinelibrary.wiley.com/doi/10.1046/j.1365-2664.1998.3540523.x/full)).

Here we use the [Multi.Index_Calculation.R](https://github.com/oleon12/PhyloDiversity/blob/master/R/Multi.Index_Caculation.R) script, which use multiple phylogenies in a **MultiPhylo** object and multiple areas in a **List** object. This script is a little complex because use loops inside loops, so if you want to modify or understand what it does, please read it with a lot of patience.

The final result are three objects: 1. A **List** object which contain the three index values for all set of areas; 2. A table with information about the species such Branch Length of each one or if is Endemic or not; and 3. A Table also with information such number of nodes of the phylogeny, the distance between the tip and the root, number of tips in the phylogeny. Please, denote that at the end, we save the index values separately by Area and Index. At the end of the script you will find a code to plot some information obtained from the data, that code was just for our personal propose, so you cand avoid it if you want.

___
