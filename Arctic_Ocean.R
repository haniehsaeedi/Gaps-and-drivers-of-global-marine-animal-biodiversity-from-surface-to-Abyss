# Retrieve occurrence data from OBIS

library("dplyr") 
library("robis") 
obisarctic <- robis::occurrence(c('Animalia'), geometry='POLYGON((180 66, 180 90, -180 90, -180 66, 180 66))')
write.csv(obisarctic,'obisarctic.csv')

# get citation for OBIS data
datasetids <- unique(obisarctic$dataset_id)

citations <- data.frame(id = datasetids, citation = NA)

for (i in 1:nrow(citations)) {
  message(citations$id[i])
  d <- dataset(datasetid = citations$id[i])
  citations$citation[i] <- d$citation
}
write.csv(citations, file = "citations_obisarctic.csv", row.names = FALSE)

#Trim OBIS data
obisarctic_t <- subset(obisarctic, absence!="TRUE")
obisarctic_t <- subset(obisarctic_t, decimalLatitude!="NA", decimalLongitude!="NA")
obisarctic_t <- subset(obisarctic_t, basisOfRecord!="FossilSpecimen")
obisarctic_t <- subset(obisarctic_t, coordinateUncertaintyInMeters<= 100000 
                         | is.na(coordinateUncertaintyInMeters)) # remove rows with "Coordinate Uncertainty" exceeding 100 km 

# https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html
# Retrieve occurrence data from GBIF
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(taxize) # for get_gbifid
library(dplyr)

# fill in your gbif.org credentials 
user <- "" 
pwd <- "" 
email <- ""

name_backbone("Animalia") # get the taxon id information

large_wkt <- "POLYGON ((180 66, 180 90, -180 90, -180 66, 180 66))"
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
occ_download_wait('0042850-240314170635999')
# After it finishes, use
gbifarctic <- occ_download_get('0042850-240314170635999') %>%
  occ_download_import()
write.csv(gbifarctic, file = "gbifarctic.csv", row.names = FALSE)

# obisarctic <- read.csv('obisarctic.csv')
# gbifarctic <- read.csv('gbifarctic.csv')

# Merge the OBIS and GBIF data
# Filter OBIS data
obisarctic_fil <- obisarctic_t %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# Filter GBIF data
gbifarctic_fil <- gbifarctic %>%
  dplyr::select(scientificName, datasetKey, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# rename columns dataset_id to merge the data 
colnames(gbifarctic_fil)[2] <- "dataset_id" 

# make the integer columns to character to merge the data 
gbifarctic_fil <- gbifarctic_fil %>% mutate(individualCount = as.character(individualCount))
gbifarctic_fil <- gbifarctic_fil %>% mutate(year = as.character(year))

# Merge OBIS, GBIF, and integrated datasets and remove the duplicates by function distinct 
arcticdata <- bind_rows(obisarctic_fil, gbifarctic_fil) %>% 
  distinct()
write.csv(arcticdata, 'arcticdata.csv')

# data cleaning using obistools
#install.packages("devtools")
#devtools::install_github("iobis/obistools")
library("obistools")

#Plot points on a map
plot_map(arcticdata, zoom = TRUE)

#Check points on land
check_onland(arcticdata)
data_on_land <- check_onland(arcticdata, report = TRUE, buffer = 100) # plot records on land with 100 meter buffer
#print(data_on_land, n=95857)
plot_map_leaflet(arcticdata[data_on_land$row,], popup = "id")
arcticdata_clean <- arcticdata[-1 * data_on_land$row,] #Remove the points on land

#Check depth using obistools
plot_map(check_depth(arcticdata_clean, depthmargin = 50), zoom = TRUE)
report <- check_depth(arcticdata_clean, report=T, depthmargin = 50)
head(report)# as only min and max depth are missing we do not need to do anything

#taxonmatch with WoRMS
names <- (arcticdata_clean$scientificName)
match_taxa(names)# all taxa matched with worms, click the info to get the unmatched list and save it as Taxmatch, then match this with WoRMS and delete the unaccepted species
arcticdata_clean_taxmatch <- subset(arcticdata_clean, !(scientificName %in% c("Glycera Lamarck, 1818", "Leuctra Stephens, 1836", "Notoproctus oculatus arcticus"))) # remove the unaccepted species
write.csv(arcticdata_clean_taxmatch,'arcticdata_clean_taxmatch.csv')
#taxmatch_matched <- read.csv('taxmatch_matched.csv')
#data_GBIF_taxmatch <- anti_join(dfg, taxmatch_matched, by = "scientificName") # Removed the unaccepted non matched WoRMS species



