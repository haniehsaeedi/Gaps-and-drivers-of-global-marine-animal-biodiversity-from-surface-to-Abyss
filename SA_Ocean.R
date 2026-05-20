# Retrieve occurrence data from OBIS

library("dplyr") 
library("robis") 
obisSAtlantic <- robis::occurrence(c('Animalia'), geometry='POLYGON((12 0,-57 0,-39 -6,-69 -45,-69 -63,24 -63,24 -31,18 -14,12 0)))')
write.csv(obisSAtlantic,'obisSAtlantic.csv')

# get citation for OBIS data
datasetids <- unique(obisSAtlantic$dataset_id)

citations <- data.frame(id = datasetids, citation = NA)

for (i in 1:nrow(citations)) {
  message(citations$id[i])
  d <- dataset(datasetid = citations$id[i])
  citations$citation[i] <- d$citation
}
write.csv(citations, file = "citations_obisSAtlnatic.csv", row.names = FALSE)

#Trim OBIS data
obisSAtlantic_t <- subset(obisSAtlantic, absence!="TRUE")
obisSAtlantic_t <- subset(obisSAtlantic_t, decimalLatitude!="NA", decimalLongitude!="NA")
obisSAtlantic_t <- subset(obisSAtlantic_t, basisOfRecord!="FossilSpecimen")
obisSAtlantic_t <- subset(obisSAtlantic_t, coordinateUncertaintyInMeters<= 100000 
                          | is.na(coordinateUncertaintyInMeters)) # remove rows with "Coordinate Uncertainty" exceeding 100 km 

# https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html
# Retrieve occurrence data from GBIF
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(taxize) # for get_gbifid
library(dplyr)

# fill in your gbif.org credentials 
user <- "GBIF_USER" 
pwd <- "GBIF_PWD" 
email <- "GBIF_EMAIL"

name_backbone("Animalia") # get the taxon id information

large_wkt <- "POLYGON ((12 0,-57 0,-39 -6,-69 -45,-69 -63,24 -63,24 -31,18 -14,12 0))"
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
occ_download_wait('0110275-240321170329656')
# After it finishes, use
gbifSAtlantic <- occ_download_get('0110275-240321170329656') %>%
  occ_download_import()
write.csv(gbifSAtlantic, file = "gbifSAtlantic.csv", row.names = FALSE)

# obisSAtlantic <- read.csv('obisSAtlantic.csv')
# gbifSAtlantic <- read.csv('gbifSAtlantic.csv')

# Merge the OBIS and GBIF data
# Filter OBIS data
obisSAtlantic_fil <- obisSAtlantic_t %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# Filter GBIF data
gbifSAtlantic_fil <- gbifSAtlantic %>%
  dplyr::select(scientificName, datasetKey, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# rename columns dataset_id to merge the data 
colnames(gbifSAtlantic_fil)[2] <- "dataset_id" 

# make the integer columns to character to merge the data 
gbifSAtlantic_fil <- gbifSAtlantic_fil %>% mutate(individualCount = as.character(individualCount))
gbifSAtlantic_fil <- gbifSAtlantic_fil %>% mutate(year = as.character(year))
gbifSAtlantic_fil <- gbifSAtlantic_fil %>% mutate(coordinatePrecision = as.character(coordinatePrecision))
gbifSAtlantic_fil <- gbifSAtlantic_fil %>% mutate(day = as.character(day))
gbifSAtlantic_fil <- gbifSAtlantic_fil %>% mutate(dateIdentified = as.character(dateIdentified))

# Merge OBIS, GBIF, and integrated datasets and remove the duplicates by function distinct 
SAtlanticdata <- bind_rows(obisSAtlantic_fil, gbifSAtlantic_fil) %>% 
  distinct()
write.csv(SAtlanticdata, 'SAtlanticdata.csv')

# data cleaning using obistools
#install.packages("devtools")
#devtools::install_github("iobis/obistools")
library("obistools")

#Plot points on a map
#plot_map(SEPacificdata, zoom = TRUE)

#Check points on land
check_onland(SAtlanticdata)
data_on_land <- check_onland(SAtlanticdata, report = TRUE, buffer = 100) # plot records on land with 100 meter buffer
#print(data_on_land, n=95857)
plot_map_leaflet(SAtlanticdata[data_on_land$row,], popup = "id")
SAtlanticdata_clean <- SAtlanticdata[-1 * data_on_land$row,] #Remove the points on land

#Check depth using obistools
plot_map(check_depth(SAtlanticdata_clean, depthmargin = 50), zoom = TRUE)
report <- check_depth(SAtlanticdata_clean, report=T, depthmargin = 50)
head(report)# as only min and max depth are missing we do not need to do anything

#remove the duplicates for taxon matching as it took so ling for taxon matching the whole dataset
library(tidyverse)
taxmatch_SA <- SAtlanticdata_clean[!duplicated(SAtlanticdata_clean$scientificName), ]

#taxonmatch with WoRMS
names <- (taxmatch_SA$scientificName)
match_taxa(names)# all taxa matched with worms, click the info to get the unmatched list and save it as Taxmatch, then match this with WoRMS and delete the unaccepted species
SAtlanticdata_clean_taxmatch <- subset(SAtlanticdata_clean, !(scientificName %in% c("Soestella De Weerdt, 2000"))) # remove the unaccepted species
write.csv(SAtlanticdata_clean_taxmatch,'SAtlanticdata_clean_taxmatch.csv')

# Remove the point data out of the shape file polygon of the study area # this took so long in R, I have done it in QGIS via clip

