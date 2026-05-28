# Retrieve occurrence data from OBIS

library("dplyr") 
library("robis") 
obis <- robis::occurrence(c('Animalia'), geometry='POLYGON((9.5 64,-28 71,-46 63.2,-73 64.5,-84 64.5,-64.3 53.8,-86 35,-102 31,-101 21.5,-94 15.5,-72 3,-41 -9,19.5 -9,13 6,-6 13.5,1 40,15.7 59,9.5 64))', startdepth = 50, enddepth = 11000)
write.csv(obisNAtlantic_1,'obisNAtlantic_1.csv')

# get citation for OBIS data
datasetids <- unique(obisNAtlantic_1$dataset_id)

citations <- data.frame(id = datasetids, citation = NA)

for (i in 1:nrow(citations)) {
  message(citations$id[i])
  d <- dataset(datasetid = citations$id[i])
  citations$citation[i] <- d$citation
}
write.csv(citations, file = "citations_obisNAtlantic_1.csv", row.names = FALSE)


#Trim OBIS data
obisNAtlantic_1_t <- subset(obisNAtlantic_1, absence!="TRUE")
obisNAtlantic_1_t <- subset(obisNAtlantic_1_t, decimalLatitude!="NA", decimalLongitude!="NA")
obisNAtlantic_1_t <- subset(obisNAtlantic_1_t, basisOfRecord!="FossilSpecimen")
obisNAtlantic_1_t <- subset(obisNAtlantic_1_t, coordinateUncertaintyInMeters<= 100000 
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

large_wkt <- "POLYGON ((9.5 64,-28 71,-46 63.2,-73 64.5,-84 64.5,-64.3 53.8,-86 35,-102 31,-101 21.5,-94 15.5,-72 3,-41 -9,19.5 -9,13 6,-6 13.5,1 40,15.7 59,9.5 64))"
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
occ_download_wait('0158882-240321170329656')
# After it finishes, use
gbifNAtlantic <- occ_download_get('0158882-240321170329656') %>%
  occ_download_import()
write.csv(gbifNAtlantic, file = "gbifNAtlantic.csv", row.names = FALSE)

# obisCaspian <- read.csv('obisCaspian.csv')
# gbifCaspian <- read.csv('gbifCaspian.csv')

# Merge the OBIS and GBIF data
# Filter OBIS data
obisNAtlantic_1_fil <- obisNAtlantic_1_t %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# Filter GBIF data
gbifNAtlantic_fil <- gbifNAtlantic %>%
  dplyr::select(scientificName, datasetKey, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, infraspecificEpithet, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

# rename columns dataset_id to merge the data 
colnames(gbifNAtlantic_fil)[2] <- "dataset_id" 

# make the integer columns to character to merge the data 
gbifNAtlantic_fil <- gbifNAtlantic_fil %>% mutate(individualCount = as.character(individualCount))
gbifNAtlantic_fil <- gbifNAtlantic_fil %>% mutate(year = as.character(year))
gbifNAtlantic_fil <- gbifNAtlantic_fil %>% mutate(coordinatePrecision = as.character(coordinatePrecision))
gbifNAtlantic_fil <- gbifNAtlantic_fil %>% mutate(day = as.character(day))
gbifNAtlantic_fil <- gbifNAtlantic_fil %>% mutate(month = as.character(month))
gbifNAtlantic_fil <- gbifNAtlantic_fil %>% mutate(coordinateUncertaintyInMeters = as.character(coordinateUncertaintyInMeters))
gbifNAtlantic_fil <- gbifNAtlantic_fil %>% mutate(dateIdentified = as.character(dateIdentified))

# Merge OBIS, GBIF, and integrated datasets and remove the duplicates by function distinct 
NAtlanticdata <- bind_rows(obisNAtlantic_1_fil, gbifNAtlantic_fil) %>% 
  distinct()
write.csv(NAtlanticdata, 'NAtlanticdata.csv')

# data cleaning using obistools
#install.packages("devtools")
#devtools::install_github("iobis/obistools")
library("obistools")

#Plot points on a map
plot_map(NAtlanticdata, zoom = TRUE)

#Check points on land
check_onland(NAtlanticdata)
data_on_land <- check_onland(NAtlanticdata, report = TRUE, buffer = 100) # plot records on land with 100 meter buffer
#print(data_on_land, n=95857)
plot_map_leaflet(NAtlanticdata[data_on_land$row,], popup = "id")
NAtlanticdata_clean <- NAtlanticdata[-1 * data_on_land$row,] #Remove the points on land
write.csv(NAtlanticdata_clean,'NAtlanticdata_clean.csv')

#Check depth using obistools
plot_map(check_depth(NAtlanticdata_clean, depthmargin = 50), zoom = TRUE)
report <- check_depth(NAtlanticdata_clean, report=T, depthmargin = 50)
head(report)# as only min and max depth are missing we do not need to do anything

#remove the duplicates for taxon matching as it took so ling for taxon matching the whole dataset
library(tidyverse)
taxmatch_NAtlantic <- NAtlanticdata_clean[!duplicated(NAtlanticdata_clean$scientificName), ]
options(max.print=1500)

#taxonmatch with WoRMS
names <- (taxmatch_NAtlantic$scientificName)
match_taxa(names)# all taxa matched exact with WoRMS, no action needed
NAtlanticdata_clean_taxmatch <- subset(NAtlanticdata_clean, !(scientificName %in% c("Acanthopria Ashmead, 1895", "Afrida Möschler, 1886", "Anochetus Mayr, 1861", "Anochetus Mayr, 1861", "Anthrax Scopoli, 1763", "Atlides Hübner, 1819", "Bracon Fabricius, 1804", "Calycopis Scudder, 1876", "Chalinea Grentzenberg, 1891", "Columbina minuta (Linnaeus, 1766)", "Cryptanura Brullé, 1846", "Delibus", "Delibus Andronov, 2007", "Emesis Fabricius, 1807", "Euderus Haliday, 1844", "Eupeodes Osten-Sacken, 1877", "Eurytoma Illiger, 1807", "Gastrops Williston, 1897", "Heterospilus Haliday, 1836", "Leptacis Förster, 1856", "Leptogaster Meigen, 1803", "Lespesia Robineau-Desvoidy, 1863", "Ligyra Newman, 1844", "Microdon Meigen, 1803", "Mixogaster Kertész, 1910", "Monohelea Kieffer, 1917", "Novakia Strobl, 1893", "Orfelia Costa, 1857", "Paguridea", "Palorus depressus (Fabricius, 1790)", "Palpada Macquart, 1834", "Paraneda Timberlake, 1943", "Pentapria Kieffer, 1905", "Pericoma Walker, 1856", "Physiphora Fallén, 1810", "Pipunculus Latreille, 1802", "Setacera Cresson, 1930", "Spioniades Hübner, 1819", "Stratiomys Geoffroy, 1762", "Synapte pecta Evans, 1955", "Theope Moore, 1857", "Torymus Dalman, 1820", "Triatoma Laporte, 1832", "Trichacis Förster, 1856", "Trichopoda Berthold, 1827", "Triorla Parks, 1968", "Villa Lioy, 1864"))) # remove the unaccepted species
write.csv(NAtlanticdata_clean_taxmatch,'NAtlanticdata_clean_taxmatch.csv')

# Remove the point data out of the shape file polygon of the study area # this took so long in R, I have done it in QGIS via clip

