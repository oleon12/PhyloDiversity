# PhyloDiversity
___

Hi, this is the repository for:

<p align="center">
 <b>Conserving Northern Andean Block: What Phylogenies Say</b><br>
 <b>Leon-Alvarado & Miranda-Esquivel (Submitted)</b><br>
</p>



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

The final result are three objects: 

- A **List** object which contain the three index values for all set of areas.
- A table with information about the species such Branch Length of each one or if is Endemic or not.
- A Table also with information such number of nodes of the phylogeny, the distance between the tip and the root, number of tips in the phylogeny. 

Please, denote that we save the index values separately by Area and Index. 

> _At the end of the script you will find a code to plot some information obtained from the data, that code was just for our personal purpose, so you cand avoid it if you want._

___

**Prioritization**

Once we have the index values it is time to prioritiza which areas contain more evolutive information and are the most importan for conservation purpose. Now we use the [RtoGIS.R](https://github.com/oleon12/PhyloDiversity/blob/master/R/RtoGIS.R) script. We decide to use the AvTD index using the grid cell size of 0.25°.

This script take all cells values and divide it in Quantiles, then the cells whith index values within the last Quantile (Q5) are identified, then the cells with index values within the Q4 and so on until the Q2. Once all grid cell were identified to which Quantile belong given their index value, a new object was created labeling each grid cell with their Quantile, then, this object was added to the attibute table in the Shape File of the grid cell and save again as a Shape File.

The Shape File generated can be reade in any GIS software, in this case we use [QGIS](https://www.qgis.org/en/site/index.html) for this task. In summary, we load the Shape File in QGIS and each cell received a color given their Quantile label, so at the end we obtained cells with 4 colours corresponding to 4 Quantile labels (Q5, Q4, Q3 and Q2), so in this way we was able to identify the Q5 cells, our cells of interest.

___

**Correlations**

Also, besides the prioritization we want to know the dependence of the Phylogenetic Diversity indices used with the species richness. For this we decided to implement a Bayesian Approach following [Kruschke. 2014](https://books.google.es/books?hl=es&lr=&id=FzvLAwAAQBAJ&oi=fnd&pg=PP1&dq=Doing+Bayesian+data+analysis:+A+tutorial+with+R,+JAGS,+and+Stan,+2nd+edn&ots=CfpkO0ydXE&sig=be7YfMbeKSV3RIttp6r-xjCfbxs#v=onepage&q=Doing%20Bayesian%20data%20analysis%3A%20A%20tutorial%20with%20R%2C%20JAGS%2C%20and%20Stan%2C%202nd%20edn&f=false) using a Bayesian Simple Linear Regression. Here, we modified the Script from Kruschke ([SimpleLinearRegressionJags.R](http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/Programs/SimpleLinearRegressionJags.R)), and created the function [BayesSlope](https://github.com/oleon12/PhyloDiversity/blob/master/R/BayesSlope.R), yet, if you use this function please cite Kruschke (2014). All correlations and plots are found in the [BayesianSLR.R](https://github.com/oleon12/PhyloDiversity/blob/master/R/BayesianSLR.R) script. For detailed information about the Bayesian Simple Linear Regression and Bayesian Statistics we suggest you to read [Kruschke. 2014](https://books.google.es/books?hl=es&lr=&id=FzvLAwAAQBAJ&oi=fnd&pg=PP1&dq=Doing+Bayesian+data+analysis:+A+tutorial+with+R,+JAGS,+and+Stan,+2nd+edn&ots=CfpkO0ydXE&sig=be7YfMbeKSV3RIttp6r-xjCfbxs#v=onepage&q=Doing%20Bayesian%20data%20analysis%3A%20A%20tutorial%20with%20R%2C%20JAGS%2C%20and%20Stan%2C%202nd%20edn&f=false).

> _Warning, to use the BayesianSlope function, you also require the functions [openGraphSaveGraph](http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/Programs/openGraphSaveGraph.R), [plotPost](http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/Programs/plotPost.R) and the [rjags package](https://cran.r-project.org/web/packages/rjags/index.html)._

___

**Endemic species**

Another question that we wanted to answer was the contribution of Endemic Species to the Phylogenetic Diversity, here we considere a species with their distribution restricted to an area of endemism as endemic. 

1. First, we recalculated the three Phylogenetic Diversity indices removing the 25%, 50%, 75% and 100% of endemic species from the distribution data with the [Endemic_MultindexCalc.R](https://github.com/oleon12/PhyloDiversity/blob/master/R/Endemic_MultindexCalc.R).

2. Then, with the results of the previous step we build the [Endemic_Plot-1.0.R](https://github.com/oleon12/PhyloDiversity/blob/master/R/Endemic_Plot-1.0.R) script to plot the results. This script fit perfectly with our data, so if you want to used for your data you can modified, however, results with more than 10 areas are very difficult to visualized.

___

**Protected Areas**

Finally, our final goal with this work was to calculated how much Phylogenetic Diversity had the protected areas within the NAB and what is their contribution to the total of the Phylogenetic Diversity. Here we calculate the Phylogenetic Diversity indices with the absence/presence matrix for tha PA areas. Here, we calculate the Phylogenetic Diversity of each PA, the Phylogenetic Diversity of all PAs and finally the Phylogenetic Diversity of the whole NAB. Once all calculations are made, we compute the contribution of the PAs to the total Phylogenetic Diversity. Everything was mde with the [PNN_analysis.R](https://github.com/oleon12/PhyloDiversity/blob/master/R/PNN_analysis.R) script.

> _Here, we avoid the contribution calculated with the PNN+ label, which correspond to the PNN1 object__

___
