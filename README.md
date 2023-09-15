# Phylogenetic Diversity and North Andean Block Conservation

## [Omar Daniel Leon-Alvarado](https://leon-alvarado.weebly.com/) and [Daniel Rafael Miranda-Esquivel](https://www.researchgate.net/profile/Daniel-Miranda-Esquivel)

</br>

<b>*Abstract*</b>
<p align="justify">

**Background.** The Northern Andean Block (NAB) harbors high biodiversity, being one of the most important areas of the Neotropics, nevertheless, the settlement of several human populations has triggered rapid transformations of the ecosystems, leading to the extinction or endangerment of many species. As phylogenetic diversity indices quantify the distinctness between species, they are an adequate tool to evaluate conservation priority areas. 

**Methods.** We reconstructed 93 phylogenies from 1255 species, and using their occurrence from GBIF, we calculated the Average Taxonomic Distinctness Index for a grid cell of 0.25° in NAB. With the index values of each cell, we classified them into quantiles, selecting those cells with values in the upper quantile (Q5) as the most important cells. We also calculated the contribution of endemic species to NAB’s phylogenetic diversity, and how much Phylogenetic Diversity is preserved within the Protected Areas (PAs).

**Results.** We found that the Andes region has the highest phylogenetic diversity values across the NAB, especially the middle and south regions of Colombia's Cordilleras. We found that endemic species contribute a small percentage of the NAB’s Phylogenetic Diversity (1.2%). The protected areas, although representing a small geographic area, harbor 40% of the total phylogenetic diversity.

**Discussion.** Although the NAB's Andean region has been identified as the most crucial area in terms of Phylogenetic Diversity, there are also regions in the Amazonian Piedemonte and Pacific lowlands that exhibit high levels of Phylogenetic Diversity. Interestingly, some protected areas (PAs) have been found to harbor a higher amount of Phylogenetic Diversity than expected, given their smaller size. While the delimitation of new PAs and species richness has been the primary factor driving the expansion of PAs, it's also essential to consider the evolutionary information of the species if we want to conserve all aspects of biodiversity, or at least cover most of them. Therefore, using Phylogenetic Diversity measures and the results of this study can contribute to expanding the PAs network and improving connectivity between PAs. This will help conserve different aspects of biodiversity and ensure that evolutionary relationships between species are also preserved.

An interactive **map** with all results is available [here](https://rpubs.com/oleon12/PhyloDiv)
</br>
</p>

---

### Some of the methods employed in the study are more effectively elucidated through visual representations, as the adage goes, 'a picture is worth a thousand words.'

## Groups selection rules

<img src="https://github.com/oleon12/PhyloDiversity/blob/master/Supplementary_Material/Supplemental_Figure_S1.png">

## Input data

The basic input data we utilized consists of an uncalibrated phylogeny featuring branch lengths, coupled with data representing species distribution. Specifically, these [distributions](https://github.com/oleon12/PhyloDiversity/tree/master/Data/Distributions) are CSV matrices that denote the presence or absence of species in individual cells.

<img src="https://github.com/oleon12/PhyloDiversity/blob/master/Supplementary_Material/Fig_1Mesa%20de%20trabajo%201.png">

## Prioritization scheme

We worked with a total of 93 taxonomic groups, each of which was equipped with its own associated phylogeny and distribution data. Our approach involved calculating the Average Taxonomic Distinctiveness (AvTD) for each of these taxonomic groups. Subsequently, we summed the AvTD values across all cells. To establish a prioritization scheme, we employed quantiles.

<img src="https://github.com/oleon12/PhyloDiversity/blob/master/Supplementary_Material/Fig_2Mesa%20de%20trabajo%201.png">

## Endemic species

We began by identifying endemic species and excluding them from our analysis. Following this step, we randomly selected and removed 25%, 50%, 75%, and 100% of these endemic species. After each removal, we calculated the Average Taxonomic Distinctiveness (AvTD) index and conducted the prioritization process. It's worth noting that we repeated the removal process 100 times for 25%, 50%, and 75% removal scenarios, while the 100% removal was performed only once.

<img src="https://github.com/oleon12/PhyloDiversity/blob/master/Supplementary_Material/Fig_3Mesa%20de%20trabajo%201.png">

## Data robustness

We used a jackknife approach to evaluate the robustness of the phylogenetic data and the resulting prioritization.

<img src="https://github.com/oleon12/PhyloDiversity/blob/master/Supplementary_Material/Fig_4Mesa%20de%20trabajo%201.png">

