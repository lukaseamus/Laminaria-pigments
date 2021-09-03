# Photosynthetic pigments of co-occurring Northeast Atlantic *Laminaria* spp. are unaffected by decomposition
This repository contains data and R code accompanying article 10.3354/meps13886 in *Marine Ecology Progress Series*.

The repository is split into two folders: **Analysis** and **Deconvolution**. The former contains all files to perform the statitsical analysis. The latter contains all files to perform the spectral deconvolution adapted from 10.1371/journal.pone.0137645. Below is a description of each file within those folders.

**Analysis**
1. `Environmental.csv`: *In situ* light and temperature data.
    - *d* = day from start of experiment as an integer
                        *day* = day from start of experiment as a numeric
                        *lux* = irradiance in lumens per square metre
   `Environmental.R`: R code to analyse these data.
2. `Decomposition.csv`: Biomass loss data.
                        *species* = *Laminaria digitata* (d), *Laminaria hyperborea* (h) or *Laminaria ochroleuca* (o)
                        *bag* = ID of galvanised steel mesh bag
                        *age* = detrital age in days
                        *loss* = biomass loss per day (%)
   `Decomposition.R`: R code to analyse these data.
3. `Extraction.csv`: Comparison of pigment extraction from fresh and lyophilised tissue.
                     *id* = sample ID
                     *pigment* = Chlorophll *a*, Chlorophyll *c* or Carotenoids
                     *wet* = fresh extract pigment concentration
                     *dry* = lyophilised extract pigment concentration
   `Extraction.R`: R code to analyse these data.
4. `Pigments.csv`: Main pigment data.
                     *id* = sample ID
                     *site* = Mount Batten (MB) or West Hoe (WH)
                     *species* = *Laminaria digitata* (d), *Laminaria hyperborea* (h) or *Laminaria ochroleuca* (o)
                     *bag* = ID of plant (P) or galvanised steel mesh bag (B)
                     *age* = detrital age in days
                     *bb.Car* = β,β-carotene
                     *Chl.a* = Chlorophll *a*
                     *Chl.c1* = Chlorophll *c*<sub>1</sub>
                     *Chl.c2* = Chlorophll *c*<sub>2</sub>
                     *Fuco* = Fucoxanthin
                     *Phe.a* = Pheophytin *a*
                     *Viola* = Violaxanthin
                     *Zea* = Zeaxanthin
   `Pigments.R`: R code to analyse these data.

**Deconvolution**

Luka Seamus Wright, 3 September 2021
