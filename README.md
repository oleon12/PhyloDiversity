# Phylogenetic Diversity and North Andean Block Conservation

## [Omar Daniel Leon-Alvarado](https://leon-alvarado.weebly.com/) and [Daniel Rafael Miranda-Esquivel](https://www.researchgate.net/profile/Daniel-Miranda-Esquivel)

</br>

<b>*Abstract*</b>
<p align="justify">

**Background.** The Northern Andean Block (NAB) harbors a high biodiversity; therefore, it is one of the most important areas in the Neotropics. Nevertheless, the settlement of several human populations has triggered rapid transformation of ecosystems, leading to the extinction or endangerment of many species.

**Methods.** As phylogenetic diversity indices quantify the distinctness between species, they are adequate tools for evaluating priority conservation areas. We reconstructed 93 phylogenies encompassing 1252 species, and utilizing their occurrence data sourced from the Global Biodiversity Information Facility (GBIF), we computed the Average Taxonomic Distinctness Index (AvTD) for each grid cell with a spatial resolution of 0.25° within the North Andean Block (NAB). The index values for each grid cell were categorized into quantiles, and grid cells displaying values falling within the upper quantile (Q5) were identified as the most significant in terms of phylogenetic diversity. We also calculated the contribution of endemic species to the overall phylogenetic diversity within the North Andean Block, specifically focusing on areas preserved within Protected Areas.

**Results.** The NAB Andean region exhibited the highest Average Taxonomic Distinctness Index (AvTD), with a concentration of high AvTD values observed in the middle and southern regions of Colombia's Cordilleras. The endemic species made a relatively modest contribution to the overall phylogenetic diversity of NAB, accounting for only 1.2% of the total. Despite their relatively small geographical footprint, Protected Areas within the NAB have emerged as crucial repositories of biodiversity, encompassing a substantial 40% of the total phylogenetic diversity in the region.

**Discussion.** Although the NAB Andean region has been identified as the most crucial area in terms of AvTD, some regions in the Amazonian Piedemonte and Pacific lowlands exhibit high AvTD levels. Interestingly, some protected areas have been found to harbor higher AvTDs than expected, given their smaller size. Although the delimitation of new PAs and species richness have been the primary factors driving the expansion of PAs, it is also essential to consider the evolutionary information of species if we want to conserve all aspects of biodiversity or at least cover most of them. Therefore, using phylogenetic diversity measures and the results of this study can contribute to expanding the PAs network and improving connectivity between PAs. This will help conserve different aspects of biodiversity and ensure that evolutionary relationships between species are preserved.


An interactive **map** with all results is available [here](https://rpubs.com/oleon12/PhyloDiv)
</br>
</p>

---

### Some of the methods employed in the study are more effectively elucidated through visual representations, as the adage goes, 'a picture is worth a thousand words.'

## Groups selection rules

<img src="https://github.com/oleon12/PhyloDiversity/blob/master/Supplementary_Material/Supplemental_Figure_S1_.png">

## Input data

The basic input data we utilized consists of an uncalibrated phylogeny featuring branch lengths, coupled with data representing species distribution. Specifically, these [distributions](https://github.com/oleon12/PhyloDiversity/tree/master/Data/Distributions) are CSV matrices that denote the presence or absence of species in individual cells.

<img src="https://github.com/oleon12/PhyloDiversity/blob/master/Supplementary_Material/Methods_Fig_1.png">

## Prioritization scheme

We worked with a total of 93 taxonomic groups, each of which was equipped with its own associated phylogeny and distribution data. Our approach involved calculating the Average Taxonomic Distinctiveness (AvTD) for each of these taxonomic groups. Subsequently, we summed the AvTD values across all cells. To establish a prioritization scheme, we employed quantiles.

<img src="https://github.com/oleon12/PhyloDiversity/blob/master/Supplementary_Material/Methods_Fig_2.png">

## Endemic species

We began by identifying endemic species and excluding them from our analysis. Following this step, we randomly selected and removed 25%, 50%, 75%, and 100% of these endemic species. After each removal, we calculated the Average Taxonomic Distinctiveness (AvTD) index and conducted the prioritization process. It's worth noting that we repeated the removal process 100 times for 25%, 50%, and 75% removal scenarios, while the 100% removal was performed only once.

<img src="https://github.com/oleon12/PhyloDiversity/blob/master/Supplementary_Material/Methods_Fig_3.png">

## Data robustness

We used a jackknife approach to evaluate the robustness of the phylogenetic data and the resulting prioritization.

<img src="https://github.com/oleon12/PhyloDiversity/blob/master/Supplementary_Material/Methods_Fig_4.png">

