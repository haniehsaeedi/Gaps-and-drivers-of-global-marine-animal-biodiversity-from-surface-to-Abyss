# Retrieve occurrence data from OBIS

library("dplyr") 
library("robis") 
obisIndia <- robis::occurrence(c('Animalia'), geometry='POLYGON((148.5 -57,147.5 -35,124.8 -26,117.1 -27,130.6 -17,130.6 -17,144.5 -22,146 -7,100.5 27,80 20.5,59.1 33.2,38.5 33.5,26.4 28.8,44.3 6,34.7 0,28.5 -23,15.2 -32,15 -58,148.5 -57)))')
write.csv(obisIndia,'obisIndia.csv')

# get citation for OBIS data
datasetids <- unique(obisIndia$dataset_id)

citations <- data.frame(id = datasetids, citation = NA)

for (i in 1:nrow(citations)) {
  message(citations$id[i])
  d <- dataset(datasetid = citations$id[i])
  citations$citation[i] <- d$citation
}
write.csv(citations, file = "citations_obisIndia.csv", row.names = FALSE)

#Trim OBIS data
obisIndia_t <- subset(obisIndia, absence!="TRUE")
obisIndia_t <- subset(obisIndia_t, decimalLatitude!="NA", decimalLongitude!="NA")
obisIndia_t <- subset(obisIndia_t, basisOfRecord!="FossilSpecimen")
obisIndia_t <- subset(obisIndia_t, coordinateUncertaintyInMeters<= 100000 
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

large_wkt <- "POLYGON ((148.5 -57,147.5 -35,124.8 -26,117.1 -27,130.6 -17,130.6 -17,144.5 -22,146 -7,100.5 27,80 20.5,59.1 33.2,38.5 33.5,26.4 28.8,44.3 6,34.7 0,28.5 -23,15.2 -32,15 -58,148.5 -57))"
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
occ_download_wait('0145596-240321170329656')
# After it finishes, use
gbifIndia <- occ_download_get('0145596-240321170329656') %>%
  occ_download_import()
write.csv(gbifIndia, file = "gbifIndia.csv", row.names = FALSE)

# Indiadata <- read.csv('Indiadata.csv')
# gbifIndia <- read.csv('gbifIndia.csv')

# Merge the OBIS and GBIF data
# Filter OBIS data
obisIndia_fil <- obisIndia_t %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# Filter GBIF data
gbifIndia_fil <- gbifIndia %>%
  dplyr::select(scientificName, datasetKey, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# rename columns dataset_id to merge the data 
colnames(gbifIndia_fil)[2] <- "dataset_id" 

# make the integer columns to character to merge the data 
gbifIndia_fil <- gbifIndia_fil %>% mutate(individualCount = as.character(individualCount))
gbifIndia_fil <- gbifIndia_fil %>% mutate(year = as.character(year))
gbifIndia_fil <- gbifIndia_fil %>% mutate(coordinatePrecision = as.character(coordinatePrecision))
gbifIndia_fil <- gbifIndia_fil %>% mutate(day = as.character(day))
gbifIndia_fil <- gbifIndia_fil %>% mutate(month = as.character(month))
gbifIndia_fil <- gbifIndia_fil %>% mutate(coordinateUncertaintyInMeters = as.character(coordinateUncertaintyInMeters))
gbifIndia_fil <- gbifIndia_fil %>% mutate(dateIdentified = as.character(dateIdentified))

# Merge OBIS, GBIF, and integrated datasets and remove the duplicates by function distinct 
Indiadata <- bind_rows(obisIndia_fil, gbifIndia_fil) %>% 
  distinct()  
write.csv(Indiadata, 'Indiadata.csv')

# data cleaning using obistools
#install.packages("devtools")
#devtools::install_github("iobis/obistools")
library("obistools")

#Plot points on a map
#plot_map(Mediterdata, zoom = TRUE)

#Check points on land
check_onland(Indiadata)
data_on_land <- check_onland(Indiadata, report = TRUE, buffer = 100) # plot records on land with 100 meter buffer
#print(data_on_land, n=95857)
plot_map_leaflet(Indiadata[data_on_land$row,], popup = "id")
Indiadata_clean <- Indiadata[-1 * data_on_land$row,] #Remove the points on land

#Check depth using obistools
plot_map(check_depth(Indiadata_clean, depthmargin = 50), zoom = TRUE)
report <- check_depth(Indiadata_clean, report=T, depthmargin = 50)
head(report)# as only min and max depth are missing we do not need to do anything

#remove the duplicates for taxon matching as it took so ling for taxon matching the whole dataset
library(tidyverse)
taxmatch_India <- Indiadata_clean[!duplicated(Indiadata_clean$scientificName), ]
write.csv(taxmatch_India,'taxmatch_India.csv')
#options(max.print=1000)

#taxonmatch with WoRMS
names <- (taxmatch_India$scientificName)
match_taxa(names)# all taxa matched with worms, click the info to get the unmatched list and save it as Taxmatch, then match this with WoRMS and delete the unaccepted species
Indiadata_clean_taxmatch <- subset(Indiadata_clean, !(scientificName %in% c("Berbyce", "Caridae", "Centrogobius Bleeker, 1874", "Delibus Andronov, 2007", "Devisina Fowler, 1931", "Leucosidae", "Myrmilla Wesmael, 1852", "Periaulax Cossmann, 1888", "Sea Hayward, 1950"))) # remove the unaccepted species
write.csv(Indiadata_clean_taxmatch,'Indiadata_clean_taxmatch.csv')

# Remove the point data out of the shape file polygon of the study area # this took so long in R, I have done it in QGIS via clip

