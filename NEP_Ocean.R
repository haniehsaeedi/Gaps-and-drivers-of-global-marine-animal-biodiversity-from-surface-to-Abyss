# Retrieve occurrence data from OBIS
library("dplyr") 
library("robis") 
obisNEPacific <- robis::occurrence(c('Animalia'), geometry='POLYGON((-76 0,-66 7,-103 28,-123 57,-180 67,-180 0,-76 0)))')
write.csv(obisNEPacific,'obisNEPacific.csv')

obisNEPacific <- read.csv('obisNEPacific.csv')

# get citation for OBIS data
datasetids <- unique(obisNEPacific$dataset_id)

citations <- data.frame(id = datasetids, citation = NA)

for (i in 1:nrow(citations)) {
  message(citations$id[i])
  d <- dataset(datasetid = citations$id[i])
  citations$citation[i] <- d$citation
}
write.csv(citations, file = "citations_obisNEPacific.csv", row.names = FALSE)

#Trim OBIS data
obisNEPacific_t <- subset(obisNEPacific, absence!="TRUE")
obisNEPacific_t <- subset(obisNEPacific_t, decimalLatitude!="NA", decimalLongitude!="NA")
obisNEPacific_t <- subset(obisNEPacific_t, basisOfRecord!="FossilSpecimen")
obisNEPacific_t <- subset(obisNEPacific_t, coordinateUncertaintyInMeters<= 100000 
                          | is.na(coordinateUncertaintyInMeters)) # remove rows with "Coordinate Uncertainty" exceeding 100 km 

# https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html
# Retrieve occurrence data from GBIF
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
install.packages("remotes")
#remotes::install_github("ropensci/taxize")
library(taxize) # for get_gbifid
library(dplyr)

# fill in your gbif.org credentials 
user <- "GBIF_USER" 
pwd <- "GBIF_PWD" 
email <- "GBIF_EMAIL"

name_backbone("Animalia") # get the taxon id information

large_wkt <- "POLYGON ((-76 0,-66 7,-103 28,-123 57,-180 67,-180 0,-76 0))"
occ_download(
  pred_within(large_wkt),
  pred_gte("depth", 0),
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"), 
  pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN"))),
  pred("taxonKey", 1),
  pred_lt("coordinateUncertaintyInMeters",100000),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)
# Check status using the download key
occ_download_wait('')
# After it finishes, use
gbifSEPacific <- occ_download_get('') %>%
  occ_download_import()
write.csv(gbifNEPacific, file = "gbifNEPacific.csv", row.names = FALSE)

# obisSEPacific <- read.csv('obisSEPacific.csv')
# gbifSEPacific <- read.csv('gbifSEPacific.csv')

# Merge the OBIS and GBIF data
# Filter OBIS data
obisSEPacific_fil <- obisSEPacific_t %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# Filter GBIF data
gbifSEPacific_fil <- gbifSEPacific %>%
  dplyr::select(scientificName, datasetKey, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# rename columns dataset_id to merge the data 
colnames(gbifSEPacific_fil)[2] <- "dataset_id" 

# make the integer columns to character to merge the data 
gbifSEPacific_fil <- gbifSEPacific_fil %>% mutate(individualCount = as.character(individualCount))
gbifSEPacific_fil <- gbifSEPacific_fil %>% mutate(year = as.character(year))
gbifSEPacific_fil <- gbifSEPacific_fil %>% mutate(coordinatePrecision = as.character(coordinatePrecision))
gbifSEPacific_fil <- gbifSEPacific_fil %>% mutate(day = as.character(day))
gbifSEPacific_fil <- gbifSEPacific_fil %>% mutate(month = as.character(month))
gbifSEPacific_fil <- gbifSEPacific_fil %>% mutate(coordinateUncertaintyInMeters = as.character(coordinateUncertaintyInMeters))
gbifSEPacific_fil <- gbifSEPacific_fil %>% mutate(dateIdentified = as.character(dateIdentified))

# Merge OBIS, GBIF, and integrated datasets and remove the duplicates by function distinct 
SEPacificdata <- bind_rows(obisSEPacific_fil, gbifSEPacific_fil) %>% 
  distinct()
write.csv(SEPacificdata, 'SEPacificdata.csv')

# data cleaning using obistools
#install.packages("devtools")
#devtools::install_github("iobis/obistools")
library("obistools")

#Plot points on a map
#plot_map(SEPacificdata, zoom = TRUE)

#Check points on land
check_onland(SEPacificdata)
data_on_land <- check_onland(SEPacificdata, report = TRUE, buffer = 100) # plot records on land with 100 meter buffer
#print(data_on_land, n=95857)
plot_map_leaflet(SEPacificdata[data_on_land$row,], popup = "id")
SEPacificdata_clean <- SEPacificdata[-1 * data_on_land$row,] #Remove the points on land

#Check depth using obistools
plot_map(check_depth(SEPacificdata_clean, depthmargin = 50), zoom = TRUE)
report <- check_depth(SEPacificdata_clean, report=T, depthmargin = 50)
head(report)# as only min and max depth are missing we do not need to do anything

#remove the duplicates for taxon matching as it took so ling for taxon matching the whole dataset
library(tidyverse)
taxmatch_SEP <- SEPacificdata_clean[!duplicated(SEPacificdata_clean$scientificName), ]

#taxonmatch with WoRMS
names <- (taxmatch_SEP$scientificName)
match_taxa(names)# all taxa matched with worms, click the info to get the unmatched list and save it as Taxmatch, then match this with WoRMS and delete the unaccepted species
SEPacificdata_clean_taxmatch <- subset(SEPacificdata_clean, !(scientificName %in% c("Neolabrus Steindachner, 1875", "Paguridea"))) # remove the unaccepted species
write.csv(SEPacificdata_clean_taxmatch,'SEPacificdata_clean_taxmatch.csv')

# Remove the point data out of the shape file polygon of the study area # this took so long in R, I have done it in QGIS via clip

