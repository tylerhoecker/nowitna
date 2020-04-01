This repository contains data and code associated with:

Hoecker TJ, Higuera PE. 2019. Forest succession and climate variability interacted to control fire activity over the last four centuries in an Alaskan boreal forest. *Landscape Ecology* DOI: doi.org/10.1007/s10980-018-00766-8

Data are stored in this repository (in `data/`), and are also in a [Dryad Data Repository: doi.org/10.5061/dryad.hg19c6n](doi.org/10.5061/dryad.hg19c6n).

To view the code and results of the analysis, open `analysis_main.html` in any web browser. 

To manipulate the code, open `analysis_main.Rmd` in RStudio. The default workflow loads all data at once from the provided `.Rdata` file (`compelete_data.RData`), and then performs the analysis. 

The code reproduces simplified versions of the published figures, and complete versions of results tables (including non-significant results, which were not included in the publication). 

Alternatively, individual datasets can be loaded using the `load_indv_data.R` script before proceeding with the analysis in `analysis_main`. `load_indv_data.R` performs some minor data wrangling procedures, and calls upon `analysis_composite.R` for developing the composite biomass burning record.

[![DOI](https://zenodo.org/badge/104701659.svg)](https://zenodo.org/badge/latestdoi/104701659)
