# INTACT-Updraft-modelling

This ESRI ArcGIS toolbox enables high-resolution updraft modelling for bird-friendly micro-siting of wind-turbines. The toolbox is spatially explicit,  utilize high resolution spatial data and remote sensing thermal imagery (Landsat 8), and is easy to implement cost-effectively in the pre-construction phase of a wind-farm development project.

Orographic updraft is estimated from a digital elevation model and proxy wind direction and wind speed.

Thermal updraft is estimated from LandSat 8 (band 3-4-10) and some proxy climate and atmospheric constants. The Landsat Thermal Remote Sensing ArcGIS algorithms for calculation of Land Surface Temperature (LST) developed by Wlawender et.al. 2012 are incorporated in the thermal updraft tool.

All Landsat 8 imagery for the validation are available from USGS Earth Explorer (https://earthexplorer.usgs.gov/). The DTM10 for the validation is supplied as supplementary materials at NINAs geodata portal (https://geodata.nina.no/layers/geonode:dem10utm). All image acquisition dates, and atmospheric and climatic constants are listed in Table 1 in the submitted manuscript.The R scripts and data files to perform the white-tailed eagle validation are supplied as supplementary materials at the Mendeley Data repository (http://dx.doi.org/10.17632/pkbsnftv5h.1).

Refererences:

Hanssen, F., May, R. and Nygård, T. (2020). High-Resolution Modeling of Uplift Landscapes can Inform Micrositing of Wind Turbines for Soaring Raptors. Journal of Environmental Management, http://dx.doi.org/10.1007/s00267-020-01318-0. 

Marques, Ana T.; Santos, Carlos D.; Hanssen, Frank Ole; Muñoz, Antonio-Román; Onrubia, Alejandro; Wikelski, Martin; Moreira, Francisco; Palmeirim, Jorge Manuel; Silva, João P. Wind turbines cause functional habitat loss for migratory soaring birds. Journal of Animal Ecology 2019 s. 1-11. https://doi.org/10.1111/1365-2656.12961

Santos, CD, Hanssen, F, Muñoz, A-R, Onrubia, A, Wikelski, M, May, R, Silva, JP (2017) Match between soaring modes of black kites and the fine-scale distribution of updrafts. Scientific Reports 7:6421. https://doi.org/10.1038/s41598-017-05319-8

Walawender, JP, Hajto, MJ, Iwaniuk, P (2012) A new ArcGIS toolset for automated mapping of land surface temperature with the use of LANDSAT satellite data, in: IEEE (Ed.), International Geoscience and Remote Sensing Symposium (IGARSS). Institute of Electrical and Electronics Engineers (IEEE), Munich, Germany, pp. 4371-4374. DOI: 10.1109/IGARSS.2012.6350405 
