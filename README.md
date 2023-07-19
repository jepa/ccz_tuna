# Climate change to drive increasing overlap between deep-sea mining and Pacific tuna fisheries

This repository is intended to support the project *Climate change to drive increasing overlap between deep-sea mining and Pacific tuna fisheries* published in 2023 in the journal *npj sustainabillity*.

**Authors**: Diva J. Amon^1,2^, Juliano Palacios-Abrantes^3^, Jeffrey C. Drazen^4^, Hannah Lily^5^, Neil Nathan^2^, Jesse van der Grient^6^, Douglas McCauley^2^

1. SpeSeas, D’Abadie, Trinidad and Tobago
2. Marine Science Institute, University of California, Santa Barbara, Santa Barbara, CA, USA
3. The Institute for the Oceans and Fisheries, University of British Columbia, Canada
4. University of Hawaii at Manoa, Honolulu, HI, USA
5. Independent Consultant, UK
6. South Atlantic Environmental Research Institute, Stanley, Falkland Islands

Citation: Amon, D.J., Palacios-Abrantes, J., Drazen, J.C. et al. Climate change to drive increasing overlap between Pacific tuna fisheries and emerging deep-sea mining industry. *npj Ocean Sustain* 2, 9 (2023)

[![DOI](https://zenodo.org/badge/29789533.svg)](https://doi.org/10.1038/s44183-023-00016-8)


# Files and folders organization:

In this repository you will find all of the code related to the data analysis. Note that all data is available from Bell 
*et al* (2020), see *Data* section below.

## Folders

- `results`: folder containing the results of the manuscript. It is sub-divided in `tbl` containing tables and `fig` containing the figures 

- `scripts`: this folder has the scripts used to do the analysis
  - `functions.R`: script that contains all of the functions needed for the analysis. Note that most functions are adopted from Bell *et al* (2020) 
  - `main_script.R`: this script contains the step-by-step analysis of the data

# Data

The raw data used in this analysis comes from the manuscript lead by J.D. Bell published in the journal *nature sustainability* in 2021 with the title "Pathways to sustaining tuna-dependent Pacific Island economies during climate change". The data is publicly available and can be found [here](https://osf.io/qa8w4/). For access to the shapefiles used please contact [The Sea Around Us](http://www.seaaroundus.org/) for the RFMO shapefile and [Jesse van der Grient](https://www.south-atlantic-research.org/2022/07/06/dr-jesse-van-der-grient-joins-the-saeri-team/) for the CCZ shapefiles.

# Relevant References

- Bell, J. D., Senina, I., Adams, T., Aumont, O., Calmettes, B., Clark, S., Dessert, M., Gehlen, M., Gorgues, T., Hampton, J., Hanich, Q., Harden-Davies, H., Hare, S. R., Holmes, G., Lehodey, P., Lengaigne, M., Mansfield, W., Menkes, C., Nicol, S., … Williams, P. (2021). Pathways to sustaining tuna-dependent Pacific Island economies during climate change. Nature Sustainability, 4(10), 900–910. https://doi.org/10.1038/s41893-021-00745-z

- Nicol, S., Ghergariu, M., Lehodey, P., Senina, I., Bell, J., & Potts, J. (2023, February 13). Tuna_Redistribution. Retrieved from osf.io/qa8w4

- Grient, J. M. A. van der, & Drazen, J. C. (2021). Potential spatial intersection between high-seas fisheries and deep-sea mining in international waters. Marine Policy, 129, 104564. https://doi.org/10.1016/j.marpol.2021.104564

***As of July of 2023, the paper has been published and this repository has been archived***
