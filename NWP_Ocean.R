# Retrieve occurrence data from OBIS

library("dplyr") 
library("robis") 
obisNWPacific <- robis::occurrence(c('Animalia'), geometry='POLYGON((95 0,180 0,180 66,95 66,95 0)))')
write.csv(obisNWPacific,'obisNWPacific.csv')

# get citation for OBIS data
datasetids <- unique(obisNWPacific$dataset_id)

citations <- data.frame(id = datasetids, citation = NA)

for (i in 1:nrow(citations)) {
  message(citations$id[i])
  d <- dataset(datasetid = citations$id[i])
  citations$citation[i] <- d$citation
}
write.csv(citations, file = "citations_obisNWPacific.csv", row.names = FALSE)

#Trim OBIS data
obisNWPacific_t <- subset(obisNWPacific, absence!="TRUE")
obisNWPacific_t <- subset(obisNWPacific_t, decimalLatitude!="NA", decimalLongitude!="NA")
obisNWPacific_t <- subset(obisNWPacific_t, basisOfRecord!="FossilSpecimen")
obisNWPacific_t <- subset(obisNWPacific_t, coordinateUncertaintyInMeters<= 100000 
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

large_wkt <- "POLYGON ((95 0,180 0,180 66,95 66,95 0))"
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
occ_download_wait('0042922-240321170329656')
# After it finishes, use
gbifNWPacific <- occ_download_get('0042922-240321170329656') %>%
  occ_download_import()
write.csv(gbifNWPacific, file = "gbifNWPacific.csv", row.names = FALSE)

# obisNWPacific <- read.csv('obissubarctic.csv')
# gbifNWPacific <- read.csv('gbifsubarctic.csv')

# Merge the OBIS and GBIF data
# Filter OBIS data
obisNWPacific_fil <- obisNWPacific_t %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# Filter GBIF data
gbifNWPacific_fil <- gbifNWPacific %>%
  dplyr::select(scientificName, datasetKey, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# rename columns dataset_id to merge the data 
colnames(gbifNWPacific_fil)[2] <- "dataset_id" 

# make the integer columns to character to merge the data 
gbifNWPacific_fil <- gbifNWPacific_fil %>% mutate(individualCount = as.character(individualCount))
gbifNWPacific_fil <- gbifNWPacific_fil %>% mutate(year = as.character(year))
gbifNWPacific_fil <- gbifNWPacific_fil %>% mutate(coordinatePrecision = as.character(coordinatePrecision))
gbifNWPacific_fil <- gbifNWPacific_fil %>% mutate(day = as.character(day))
gbifNWPacific_fil <- gbifNWPacific_fil %>% mutate(month = as.character(month))
gbifNWPacific_fil <- gbifNWPacific_fil %>% mutate(coordinateUncertaintyInMeters = as.character(coordinateUncertaintyInMeters))
gbifNWPacific_fil <- gbifNWPacific_fil %>% mutate(dateIdentified = as.character(dateIdentified))

# Merge OBIS, GBIF, and integrated datasets and remove the duplicates by function distinct 
NWPacificdata <- bind_rows(obisNWPacific_fil, gbifNWPacific_fil) %>% 
  distinct()
write.csv(NWPacificdata, 'NWPacificdata.csv')

# data cleaning using obistools
#install.packages("devtools")
#devtools::install_github("iobis/obistools")
library("obistools")

#Plot points on a map
plot_map(NWPacificdata, zoom = TRUE)

#Check points on land
check_onland(NWPacificdata)
data_on_land <- check_onland(NWPacificdata, report = TRUE, buffer = 100) # plot records on land with 100 meter buffer
#print(data_on_land, n=95857)
plot_map_leaflet(NWPacificdata[data_on_land$row,], popup = "id")
NWPacificdata_clean <- NWPacificdata[-1 * data_on_land$row,] #Remove the points on land

#Check depth using obistools
plot_map(check_depth(NWPacificdata_clean, depthmargin = 50), zoom = TRUE)
report <- check_depth(NWPacificdata_clean, report=T, depthmargin = 50)
head(report)# as only min and max depth are missing we do not need to do anything

#taxonmatch with WoRMS
names <- (NWPacificdata_clean$scientificName)
match_taxa(names)# all taxa matched with worms, click the info to get the unmatched list and save it as Taxmatch, then match this with WoRMS and delete the unaccepted species
NWPacificdata_clean_taxmatch <- subset(NWPacificdata_clean, !(scientificName %in% c("Caridae", "Delibus", "Golfingia murinae murinae", "Ostraea imbricata (Lamarck, 1819)", "Paguridea", "Ringiculina doliaris Gould, 1860", "Talostolida teres teres", "Typosyllis aciculata orientalis"))) # remove the unaccepted species
write.csv(NWPacificdata_clean_taxmatch,'NWPacificdata_clean_taxmatch.csv')

# Remove the point data out of the shape file polygon of the study area # this took so long in R, I have done it in QGIS via clip
NWPacificdata_clean_taxmatch <- read.csv('NWPacificdata_clean_taxmatch.csv')
library(sf)
library(sp)
study_area <- st_read("NWP_Ocean_diss.shp") %>% #add the shape file polygon
  sf::st_as_sf()
#plot(study_area)
data_sf <- NWPacificdata_clean_taxmatch[, c("scientificName", "dataset_id", "decimalLatitude", "decimalLongitude", "depth", "occurrenceID", "basisOfRecord", "kingdom", "order", "species", "countryCode", "occurrenceStatus", "coordinatePrecision", "day", "institutionCode", "recordNumber", "license", "typeStatus", "phylum", "family", "infraspecificEpithet", "locality", "individualCount", "month", "collectionCode", "identifiedBy", "rightsHolder", "establishmentMeans", "class", "genus", "taxonRank", "stateProvince", "coordinateUncertaintyInMeters", "eventDate", "year","catalogNumber", "dateIdentified", "recordedBy")]
data_sf <- data_sf %>% st_as_sf(coords = c('decimalLongitude','decimalLatitude'))
st_crs(data_sf) = 4326
sf_use_s2(FALSE) # to switch of Spherical geometry (s2)
out <- st_intersection(data_sf, study_area) #here remeber that the shape file must be dissolved as one polygon
library(magrittr) #assign coordinates to the dataframe
df <- out %>%
  dplyr::mutate(decimalLongitude = sf::st_coordinates(.)[,1],
                decimalLatitude = sf::st_coordinates(.)[,2])
NWPacificdata_final <- st_drop_geometry(df) # this is to drop the geometry column to be able to run the data_on_land
write.csv(NWPacificdata_final,'NWPacificdata_final.csv')
