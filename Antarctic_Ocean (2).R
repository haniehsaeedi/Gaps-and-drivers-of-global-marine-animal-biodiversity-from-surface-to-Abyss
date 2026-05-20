# Retrieve occurrence data from OBIS

library("dplyr") 
library("robis") 
obisantarctic <- robis::occurrence(c('Animalia'), geometry='POLYGON((180 -60, 180 -90, -180 -90, -180 -60, 180 -60))')
write.csv(obisantarctic,'obisantarctic.csv')

# get citation for OBIS data
datasetids <- unique(obisantarctic$dataset_id)

citations <- data.frame(id = datasetids, citation = NA)

for (i in 1:nrow(citations)) {
  message(citations$id[i])
  d <- dataset(datasetid = citations$id[i])
  citations$citation[i] <- d$citation
}
write.csv(citations, file = "citations_obisantarctic.csv", row.names = FALSE)

#Trim OBIS data
obisantarctic_t <- subset(obisantarctic, absence!="TRUE")
obisantarctic_t <- subset(obisantarctic_t, decimalLatitude!="NA", decimalLongitude!="NA")
obisantarctic_t <- subset(obisantarctic_t, basisOfRecord!="FossilSpecimen")
obisantarctic_t <- subset(obisantarctic_t, coordinateUncertaintyInMeters<= 100000 
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

large_wkt <- "POLYGON ((180 -60, 180 -90, -180 -90, -180 -60, 180 -60))"
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
occ_download_wait('0004300-240321170329656')
# After it finishes, use
gbifantarctic <- occ_download_get('0004300-240321170329656') %>%
  occ_download_import()
write.csv(gbifantarctic, file = "gbifantarctic.csv", row.names = FALSE)# the rGBIF was retrieving only 5o records, so I have downloaded it from GBIF manually 

# obisantarctic <- read.csv('obisantarctic.csv')
# gbifantarctic <- read.csv('gbifantarctic.csv', sep='\t')

# Merge the OBIS and GBIF data
# Filter OBIS data
obisantarctic_fil <- obisantarctic_t %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# Filter GBIF data
gbifantarctic_fil <- gbifantarctic %>%
  dplyr::select(scientificName, datasetKey, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# rename columns dataset_id to merge the data 
colnames(gbifantarctic_fil)[2] <- "dataset_id" 

# make the integer columns to character to merge the data 
gbifantarctic_fil <- gbifantarctic_fil %>% mutate(individualCount = as.character(individualCount))
gbifantarctic_fil <- gbifantarctic_fil %>% mutate(year = as.character(year))
gbifantarctic_fil <- gbifantarctic_fil %>% mutate(coordinatePrecision = as.character(coordinatePrecision))
gbifantarctic_fil <- gbifantarctic_fil %>% mutate(day = as.character(day))
gbifantarctic_fil <- gbifantarctic_fil %>% mutate(month = as.character(month))
gbifantarctic_fil <- gbifantarctic_fil %>% mutate(coordinateUncertaintyInMeters = as.character(coordinateUncertaintyInMeters))

# Merge OBIS, GBIF, and integrated datasets and remove the duplicates by function distinct 
antarcticdata <- bind_rows(obisantarctic_fil, gbifantarctic_fil) %>% 
  distinct()

# data cleaning using obistools
#install.packages("devtools")
#devtools::install_github("iobis/obistools")
library("obistools")

#Plot points on a map
plot_map(antarcticdata, zoom = TRUE)

#Check points on land
check_onland(antarcticdata)
data_on_land <- check_onland(antarcticdata, report = TRUE, buffer = 100) # plot records on land with 100 meter buffer
#print(data_on_land, n=95857)
plot_map_leaflet(antarcticdata[data_on_land$row,], popup = "id")
antarcticdata_clean <- antarcticdata[-1 * data_on_land$row,] #Remove the points on land

#Check depth using obistools
plot_map(check_depth(antarcticdata_clean, depthmargin = 50), zoom = TRUE)
report <- check_depth(antarcticdata_clean, report=T, depthmargin = 50)
head(report)# as only min and max depth are missing we do not need to do anything

#taxonmatch with WoRMS
#matched <- obistools::match_taxa(antarcticdata_clean$scientificName, ask = FALSE) %>%
 # select(scientificName, scientificNameID)
#matched

#antarcticdata_clean_taxmatch <- bind_cols(antarcticdata_clean_taxmatch, matched)
#antarcticdata_clean_taxmatch

#non_matches <- antarcticdata_clean_taxmatch %>%
 # filter(is.na(scientificNameID)) %>%
 # group_by(scientificName) %>%
 # summarize(n = n()) %>%
#  arrange(desc(n))

#write.table(non_matches, file = file.path(output_dir, "nonmatches.txt"), sep = "\t", row.names = FALSE, na = "", quote = FALSE)
#non_matches

names <- (antarcticdata_clean$scientificName)
match_taxa(names)# all taxa matched with worms, click the info to get the unmatched list and save it as Taxmatch, then match this with WoRMS and delete the unaccepted species
antarcticdata_clean_taxmatch <- subset(antarcticdata_clean, !(scientificName %in% c("Soestella De Weerdt, 2000", "Oswaldaria Stechow, 1921"))) # remove the unaccepted species
write.csv(antarcticdata_clean_taxmatch,'antarcticdata_clean_taxmatch.csv')

# Remove the point data out of the shape file polygon of the study area
#antarcticdata_clean_taxmatch <- read.csv('antarcticdata_clean_taxmatch.csv')
library(sf)
library(sp)
study_area <- st_read("antarctic.shp") %>% #add the shape file polygon
  sf::st_as_sf()
#plot(study_area)
data_sf <- antarcticdata_clean_taxmatch[, c("scientificName", "dataset_id", "decimalLatitude", "decimalLongitude", "depth", "occurrenceID", "basisOfRecord", "kingdom", "order", "species", "countryCode", "occurrenceStatus", "coordinatePrecision", "day", "institutionCode", "recordNumber", "license", "typeStatus", "phylum", "family", "infraspecificEpithet", "locality", "individualCount", "month", "collectionCode", "identifiedBy", "rightsHolder", "establishmentMeans", "class", "genus", "taxonRank", "stateProvince", "coordinateUncertaintyInMeters", "eventDate", "year","catalogNumber", "dateIdentified", "recordedBy")]
data_sf <- data_sf %>% st_as_sf(coords = c('decimalLongitude','decimalLatitude'))
st_crs(data_sf) = 4326
sf_use_s2(FALSE) # to switch of Spherical geometry (s2)
out <- st_intersection(data_sf, study_area)
library(magrittr) #assign coordinates to the dataframe
df <- out %>%
  dplyr::mutate(decimalLongitude = sf::st_coordinates(.)[,1],
                decimalLatitude = sf::st_coordinates(.)[,2])
antarcticdata_final <- st_drop_geometry(df) # this is to drop the geometry column to be able to run the data_on_land
write.csv(antarcticdata_final,'antarcticdata_final.csv')

