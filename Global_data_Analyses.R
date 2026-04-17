Sys.setenv(LANG = "en")
library(readxl)
library(openxlsx)
library(tidyverse)
library(sf)
library(vegan)
library(pvclust)
library(dplyr) 
library(ggplot2)
library(robis)
library(obistools)
library(nortest) # for Anderson-Darling test
library(stringi) # for encoding UTF-8
library(corrplot)
library(brms)           # For Bayesian Methods
library(randomForest)   # For Machine Learning
library(tidyr)          # For data manipulation
library(purrr)          # For functional programming
library(diptest)  # For Hartigan's Dip Test

# Merge all the OBIS citations
# Get all CSV file names in the working directory
file_list <- list.files(pattern = "\\.csv$", full.names = TRUE)

# Read and merge all CSV files
merged_OBIS_Citations <- map_dfr(file_list, read.csv)
write.csv(merged_OBIS_Citations, file = "merged_OBIS_Citations.csv")

# Import and Formatting
# 1 Caspian Sea
Caspian_data <- read_csv("Caspiandata_clean_taxmatch_clip.csv")
#Add a column named Antarctic Ocean
Caspian_data$Ocean_Sea <- c("Caspian_Sea")

# 2 Antarctic Ocean
Antarctic_data <- read_csv("antarcticdata_clean_taxmatch_clip.csv")
#Add a column named Antarctic Ocean
Antarctic_data$Ocean_Sea <- c("Antarctic_Ocean")

# 3 NA Ocean
NA_data <- read_csv("entireNAtlanticdata_clean_taxmatch_clip.csv")
#Add a column named Antarctic Ocean
NA_data$Ocean_Sea <- c("NA_Ocean")

# 4 Baltic Sea
Baltic_data <- read_csv("Balticdata_clean_taxmatch_clip.csv")
#Add a column named Antarctic Ocean
Baltic_data$Ocean_Sea <- c("Baltic_Sea")

# 5 Arctic Ocean
Arctic_data <- read_csv("entirearcticdata_clean_taxmatch_clip.csv")
#Add a column named Antarctic Ocean
Arctic_data$Ocean_Sea <- c("Arctic_Ocean")

# 6 Indian Ocean
Indian_data <- read_csv("Indiandata_clean_taxmatch_clip.csv")
#Add a column named Antarctic Ocean
Indian_data$Ocean_Sea <- c("Indian_Ocean")

# 7 Mediterranean Sea
Mediter_data <- read_csv("Mediterdata_clean_taxmatch_clip.csv")
#Add a column named Antarctic Ocean
Mediter_data$Ocean_Sea <- c("Mediterranean_Sea")

# 8 NWP Ocean
NWP_data <- read_csv("NWPacificdata_clean_taxmatch_clip.csv")
#Add a column named Antarctic Ocean
NWP_data$Ocean_Sea <- c("NWP_Oean")

# 9 NEP Ocean
NEP_data <- read_csv("NEPacificdata_clean_taxmatch_clip.csv")
#Add a column named Antarctic Ocean
NEP_data$Ocean_Sea <- c("NEP_Oean")

# 10 SA Ocean
SA_data <- read_csv("SAtlanticdata_clean_taxmatch_clip.csv")
#Add a column named Antarctic Ocean
SA_data$Ocean_Sea <- c("SA_Oean")

# 11 SEP Ocean
SEP_data <- read_csv("SEPacificdata_clean_taxmatch_clip.csv")
#Add a column named Antarctic Ocean
SEP_data$Ocean_Sea <- c("SEP_Oean")

# 12 SWP Ocean
SWP_data <- read_csv("SWPacificdata_clean_taxmatch_clip.csv")
#Add a column named Antarctic Ocean
SWP_data$Ocean_Sea <- c("SWP_Oean")

# 1 Filter Caspian_data
Caspian_data_fil <- Caspian_data %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy)

# 2 Filter Antarctic data
Antarctic_data_fil <- Antarctic_data %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy)

# 3 Filter Arctic data
Arctic_data_fil <- Arctic_data %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) 

# 4 Filter NA data
NA_data_fil <- NA_data %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) 

# 5 Filter Baltic data
Baltic_data_fil <- Baltic_data %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) 

# 6 Filter Indian data
Indian_data_fil <- Indian_data %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy)

# 7 Filter NWP data
NWP_data_fil <- NWP_data %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) 

# 8 Filter SWP data
SWP_data_fil <- SWP_data %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) 

# 9 Filter SEP data
SEP_data_fil <- SEP_data %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) 

# 10 Filter Mediter data
Mediter_data_fil <- Mediter_data %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) 

# 11 Filter SA data
SA_data_fil <- SA_data %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) 

# 12 Filter NEP data
NEP_data_fil <- NEP_data %>%
  dplyr::select(scientificName, dataset_id, decimalLatitude, decimalLongitude, depth, occurrenceID, basisOfRecord, kingdom, order, species, countryCode, occurrenceStatus, coordinatePrecision, day, institutionCode, recordNumber, license, typeStatus, phylum, family, locality, individualCount, month, collectionCode, identifiedBy, rightsHolder, establishmentMeans, class, genus, taxonRank, stateProvince, coordinateUncertaintyInMeters, eventDate, year,catalogNumber, dateIdentified, recordedBy) 

Caspian_data_fil <- Caspian_data_fil %>% mutate(recordNumber = as.character(recordNumber))
Antarctic_data_fil <- Antarctic_data_fil %>% mutate(recordNumber = as.character(recordNumber))
Arctic_data_fil <- Antarctic_data_fil %>% mutate(recordNumber = as.character(recordNumber))
NA_data_fil <- NA_data_fil %>% mutate(recordNumber = as.character(recordNumber))
Baltic_data_fil <- Baltic_data_fil %>% mutate(recordNumber = as.character(recordNumber))
Indian_data_fil <- Indian_data_fil %>% mutate(recordNumber = as.character(recordNumber))
NWP_data_fil <- NWP_data_fil %>% mutate(recordNumber = as.character(recordNumber))
SWP_data_fil <- SWP_data_fil %>% mutate(recordNumber = as.character(recordNumber))
SEP_data_fil <- SEP_data_fil %>% mutate(recordNumber = as.character(recordNumber))
Mediter_data_fil <- Mediter_data_fil %>% mutate(recordNumber = as.character(recordNumber))
SA_data_fil <- SA_data_fil %>% mutate(recordNumber = as.character(recordNumber))
NEP_data_fil <- NEP_data_fil %>% mutate(recordNumber = as.character(recordNumber))

Caspian_data_fil <- Caspian_data_fil %>% mutate(day = as.character(day))
Antarctic_data_fil <- Antarctic_data_fil %>% mutate(day = as.character(day))
Arctic_data_fil <- Antarctic_data_fil %>% mutate(day = as.character(day))
NA_data_fil <- NA_data_fil %>% mutate(day = as.character(day))
Baltic_data_fil <- Baltic_data_fil %>% mutate(day = as.character(day))
Indian_data_fil <- Indian_data_fil %>% mutate(day = as.character(day))
NWP_data_fil <- NWP_data_fil %>% mutate(day = as.character(day))
SWP_data_fil <- SWP_data_fil %>% mutate(day = as.character(day))
SEP_data_fil <- SEP_data_fil %>% mutate(day = as.character(day))
Mediter_data_fil <- Mediter_data_fil %>% mutate(day = as.character(day))
SA_data_fil <- SA_data_fil %>% mutate(day = as.character(day))
NEP_data_fil <- NEP_data_fil %>% mutate(day = as.character(day))

Caspian_data_fil <- Caspian_data_fil %>% mutate(month = as.character(month))
Antarctic_data_fil <- Antarctic_data_fil %>% mutate(month = as.character(month))
Arctic_data_fil <- Antarctic_data_fil %>% mutate(month = as.character(month))
NA_data_fil <- NA_data_fil %>% mutate(month = as.character(month))
Baltic_data_fil <- Baltic_data_fil %>% mutate(month = as.character(month))
Indian_data_fil <- Indian_data_fil %>% mutate(month = as.character(month))
NWP_data_fil <- NWP_data_fil %>% mutate(month = as.character(month))
SWP_data_fil <- SWP_data_fil %>% mutate(month = as.character(month))
SEP_data_fil <- SEP_data_fil %>% mutate(month = as.character(month))
Mediter_data_fil <- Mediter_data_fil %>% mutate(month = as.character(month))
SA_data_fil <- SA_data_fil %>% mutate(month = as.character(month))
NEP_data_fil <- NEP_data_fil %>% mutate(month = as.character(month))

Caspian_data_fil <- Caspian_data_fil %>% mutate(dateIdentified = as.character(dateIdentified))
Antarctic_data_fil <- Antarctic_data_fil %>% mutate(dateIdentified = as.character(dateIdentified))
Arctic_data_fil <- Antarctic_data_fil %>% mutate(dateIdentified = as.character(dateIdentified))
NA_data_fil <- NA_data_fil %>% mutate(dateIdentified = as.character(dateIdentified))
Baltic_data_fil <- Baltic_data_fil %>% mutate(dateIdentified = as.character(dateIdentified))
Indian_data_fil <- Indian_data_fil %>% mutate(dateIdentified = as.character(dateIdentified))
NWP_data_fil <- NWP_data_fil %>% mutate(dateIdentified = as.character(dateIdentified))
SWP_data_fil <- SWP_data_fil %>% mutate(dateIdentified = as.character(dateIdentified))
SEP_data_fil <- SEP_data_fil %>% mutate(dateIdentified = as.character(dateIdentified))
Mediter_data_fil <- Mediter_data_fil %>% mutate(dateIdentified = as.character(dateIdentified))
SA_data_fil <- SA_data_fil %>% mutate(dateIdentified = as.character(dateIdentified))
NEP_data_fil <- NEP_data_fil %>% mutate(dateIdentified = as.character(dateIdentified))

Caspian_data_fil <- Caspian_data_fil %>% mutate(year = as.character(year))
Antarctic_data_fil <- Antarctic_data_fil %>% mutate(year = as.character(year))
Arctic_data_fil <- Antarctic_data_fil %>% mutate(year = as.character(year))
NA_data_fil <- NA_data_fil %>% mutate(year = as.character(year))
Baltic_data_fil <- Baltic_data_fil %>% mutate(year = as.character(year))
Indian_data_fil <- Indian_data_fil %>% mutate(year = as.character(year))
NWP_data_fil <- NWP_data_fil %>% mutate(year = as.character(year))
SWP_data_fil <- SWP_data_fil %>% mutate(year = as.character(year))
SEP_data_fil <- SEP_data_fil %>% mutate(year = as.character(year))
Mediter_data_fil <- Mediter_data_fil %>% mutate(year = as.character(year))
SA_data_fil <- SA_data_fil %>% mutate(year = as.character(year))
NEP_data_fil <- NEP_data_fil %>% mutate(year = as.character(year))

# Merge all dataset
Global_data <- bind_rows(Caspian_data_fil, Antarctic_data_fil, Arctic_data_fil, NA_data_fil, Baltic_data_fil, Indian_data_fil, NWP_data_fil, SWP_data_fil, SEP_data_fil, Mediter_data_fil, SA_data_fil, NEP_data_fil) 
save(Global_data, file="Global_data.rda") #to save as .rda
load("Global_data.rda") #to load the .rda
#write.csv(Global_data, 'Global_data.csv')

#Data formatting
Global_data <- Global_data %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) # remove NAs
Global_data$scientificName <- stri_encode(Global_data$scientificName, "", "UTF-8") # re-mark encodings

Global_data$phylum <- as.factor(Global_data$phylum)
Global_data <- mutate(Global_data, lat_5 = ceiling(decimalLatitude/5) * 5)
Global_data <- mutate(Global_data, dep_rnd = floor(depth/100) * 100)

#Get summary of the dataset
summary(Global_data)
Global_data_table <- summary(Global_data)
write.csv(Global_data_table, "Global_data_table.csv")
subset(Global_data, select = scientificName) %>% n_distinct() # get number of species

#Get the number of distinct scientific names
distinct_species_count <- subset(Global_data, select = scientificName) %>% n_distinct()
print(distinct_species_count)

# Count the total number of scientific names per latitude bin
total_counts_per_lat <- Global_data %>%
  group_by(lat_5) %>%
  summarise(total_scientific_names = n())

# Calculate max, min, and mean of the total number of scientific names
summary_stats <- total_counts_per_lat %>%
  summarise(
    max_count = max(total_scientific_names),
    min_count = min(total_scientific_names),
    mean_count = mean(total_scientific_names)
  )
print(summary_stats)

# Count the number of unique scientific names per latitude bin
unique_counts_per_lat <- Global_data %>%
  group_by(lat_5) %>%
  summarise(unique_scientific_names = n_distinct(scientificName))

# Calculate max, min, and mean of the unique scientific name counts
summary_stats <- unique_counts_per_lat %>%
  summarise(
    max_count = max(unique_scientific_names),
    min_count = min(unique_scientific_names),
    mean_count = mean(unique_scientific_names)
  )
print(summary_stats)

# Count the number of unique phylum names per latitude bin
unique_phylum_per_lat <- Global_data %>%
  group_by(lat_5) %>%
  summarise(unique_phylum_names = n_distinct(phylum))

# Count the total number of scientific names per depth bin
total_counts_per_dep <- Global_data %>%
  group_by(dep_rnd) %>%
  summarise(total_scientific_names = n())

# Calculate max, min, and mean of the total number of scientific names
summary_stats <- total_counts_per_dep %>%
  summarise(
    max_count = max(total_scientific_names),
    min_count = min(total_scientific_names),
    mean_count = mean(total_scientific_names)
  )
print(summary_stats)

# Count the number of unique scientific names per depth bin
unique_counts_per_dep <- Global_data %>%
  group_by(dep_rnd) %>%
  summarise(unique_scientific_names = n_distinct(scientificName))

# Calculate max, min, and mean of the unique scientific name counts
summary_stats <- unique_counts_per_dep %>%
  summarise(
    max_count = max(unique_scientific_names),
    min_count = min(unique_scientific_names),
    mean_count = mean(unique_scientific_names)
  )
print(summary_stats)

#Get summary of the phylum column
Global_Phylum_data_table <- Global_data %>%
  group_by(phylum) %>%
  summarise_each(funs(n_distinct))
write.csv(Global_Phylum_data_table, "Global_Phylum_data_table.csv")

#divide to shallow, meso, and deep
#shallow
Global_data_shallow <- filter(Global_data, depth<=200)
summary(Global_data_shallow)
subset(Global_data_shallow, select = scientificName) %>% n_distinct() # get number of species
save(Global_data_shallow, file="Global_data_shallow.rda") #to save as .rda
plot_map(Global_data_shallow, zoom = TRUE) +
  theme_bw()

#meso
Global_data_meso <- filter(Global_data, depth>200 & depth<=500)
summary(Global_data_meso)
subset(Global_data_meso, select = scientificName) %>% n_distinct() # get number of species
save(Global_data_meso, file="Global_data_meso.rda") #to save as .rda
plot_map(Global_data_meso, zoom = TRUE) +
  theme_bw()

#deep
Global_data_deep <- filter(Global_data, depth>500)
summary(Global_data_deep)
subset(Global_data_deep, select = scientificName) %>% n_distinct() # get number of species
save(Global_data_deep, file="Global_data_deep.rda") #to save as .rda
plot_map(Global_data_deep, zoom = TRUE) +
  theme_bw()

#Plot points on a map
plot_map(Global_data, zoom = TRUE) +
  theme_bw()

#Format the dataset
Global_data$phylum <- as.factor(Global_data$phylum)
Global_data <- mutate(Global_data, lat_5 = ceiling(decimalLatitude/5) * 5)
Global_data <- mutate(Global_data, dep_rnd = floor(depth/100) * 100)
Global_data <- mutate(Global_data, dep_rat = factor(case_when(depth <= 200 ~ "shallow", depth>200 & depth<=500 ~ "meso", depth > 500 ~ "deep")))

#divide data to seven taxa
#Annelida
Global_data_Annelida <- filter(Global_data, phylum=="Annelida")
save(Global_data_Annelida, file="Global_data_Annelida.rda") #to save as .rda
summary(Global_data_Annelida)
Global_data_Annelida_table <- summary(Global_data_Annelida)
write.csv(Global_data_Annelida_table, "Global_data_Annelida_table.csv")
subset(Global_data_Annelida, select = scientificName) %>% n_distinct() # get number of species

#Arthropoda
Global_data_Arthropoda <- filter(Global_data, phylum=="Arthropoda")
save(Global_data_Arthropoda, file="Global_data_Arthropoda.rda") #to save as .rda
summary(Global_data_Arthropoda)
Global_data_Arthropoda_table <- summary(Global_data_Arthropoda)
write.csv(Global_data_Arthropoda_table, "Global_data_Arthropoda_table.csv")
subset(Global_data_Arthropoda, select = scientificName) %>% n_distinct() # get number of species

#Chordata
Global_data_Chordata <- filter(Global_data, phylum=="Chordata")
save(Global_data_Chordata, file="Global_data_Chordata.rda") #to save as .rda
summary(Global_data_Chordata)
Global_data_Chordata_table <- summary(Global_data_Chordata)
write.csv(Global_data_Chordata_table, "Global_data_Chordata_table.csv")
subset(Global_data_Chordata, select = scientificName) %>% n_distinct() # get number of species

#Cnidaria
Global_data_Cnidaria <- filter(Global_data, phylum=="Cnidaria")
save(Global_data_Cnidaria, file="Global_data_Cnidaria.rda") #to save as .rda
summary(Global_data_Cnidaria)
Global_data_Cnidaria_table <- summary(Global_data_Cnidaria)
write.csv(Global_data_Cnidaria_table, "Global_data_Cnidaria_table.csv")
subset(Global_data_Cnidaria, select = scientificName) %>% n_distinct() # get number of species

#Echinodermata
Global_data_Echinodermata <- filter(Global_data, phylum=="Echinodermata")
save(Global_data_Echinodermata, file="Global_data_Echinodermata.rda") #to save as .rda
summary(Global_data_Echinodermata)
Global_data_Echinodermata_table <- summary(Global_data_Echinodermata)
write.csv(Global_data_Echinodermata_table, "Global_data_Echinodermata_table.csv")
subset(Global_data_Echinodermata, select = scientificName) %>% n_distinct() # get number of species

#Mollusca
Global_data_Mollusca <- filter(Global_data, phylum=="Mollusca")
save(Global_data_Mollusca, file="Global_data_Mollusca.rda") #to save as .rda
summary(Global_data_Mollusca)
Global_data_Mollusca_table <- summary(Global_data_Mollusca)
write.csv(Global_data_Mollusca_table, "Global_data_Mollusca_table.csv")
subset(Global_data_Mollusca, select = scientificName) %>% n_distinct() # get number of species

#Porifera
Global_data_Porifera <- filter(Global_data, phylum=="Porifera")
save(Global_data_Porifera, file="Global_data_Porifera.rda") #to save as .rda
summary(Global_data_Porifera)
Global_data_Porifera_table <- summary(Global_data_Porifera)
write.csv(Global_data_Porifera_table, "Global_data_Porifera_table.csv")
subset(Global_data_Porifera, select = scientificName) %>% n_distinct() # get number of species

# Global_data has a column 'phylum' with specific phylum names
phylum_colors <- c(
  "Acanthocephala" = "tan3", 
  "Annelida" = "purple", 
  "Arthropoda" = "red2",
  "Brachiopoda" = "cornflowerblue", 
  "Bryozoa" = "magenta", 
  "Cephalorhyncha" = "darkolivegreen4",
  "Chaetognatha" = "wheat4",
  "Chordata" = "black",
  "Cnidaria" = "tan4", 
  "Coelenterata" = "indianred1", 
  "Ctenophora" = "firebrick4",
  "Cycliophora" = "yellowgreen",
  "Dicyemida" = "lightsalmon",
  "Echinodermata" = "orange",
  "Entoprocta" = "gold",
  "Gastrotricha" = "violet",
  "Gnathostomulida" = "#DDAD4B", 
  "Hemichordata" = "darkgreen", 
  "Kinorhyncha" = "seagreen2",
  "Loricifera" = "moccasin", 
  "Mollusca" = "blue", 
  "Nematoda" = "yellow",
  "Nematomorpha" = "cadetblue1",
  "Nemertea" = "lightgrey",
  "Orthonectida" = "tomato3", 
  "Phoronida" = "lightgreen", 
  "Platyhelminthes" = "lightyellow",
  "Porifera" = "chartreuse", 
  "Priapulida" = "firebrick", 
  "Rotifera" = "wheat3",
  "Sipuncula" = "lightblue",
  "Tardigrada" = "wheat", 
  "Xenacoelomorpha" = "pink"
)

# Plot records per latitude per phylum
Lat_Num_Rec_Phy <- ggplot(Global_data, aes(x = lat_5, fill = phylum)) + 
  geom_bar(position="dodge") +
  scale_y_continuous(limits=c(0, 4500000), breaks=seq(0, 4500000, by=1000000), expand = c(0, 0)) +
  scale_x_continuous(limits=c(-90,90), breaks=seq(-90, 90, by=20), expand = c(0, 0)) +
  xlab("Latitude (degree)") + ylab("Number of Records") +
  scale_fill_manual(values = phylum_colors) +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),    
    legend.text = element_text(size = 8),          # Adjust legend text size
    legend.title = element_text(size = 10),        # Adjust legend title size
    legend.key.size = unit(0.5, "cm")               # Adjust legend key size
  )
print (Lat_Num_Rec_Phy)

ggsave("Lat_Num_Rec_Phy.tiff", plot = Lat_Num_Rec_Phy, width = 12, height = 6, dpi = 600)

# Plot records per latitude (black and white)
ggplot(Global_data, aes(x = lat_5)) + 
  geom_bar(position="dodge") +
  scale_y_continuous(limits=c(0, 8000000), breaks=seq(0, 8000000, by=1000000), expand = c(0, 0)) +
  scale_x_continuous(limits=c(-90,90), breaks=seq(-90, 90, by=20), expand = c(0, 0)) +
  xlab("Latitude (degree)") + ylab("Number of Records") +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = 22, color = "black"),     
    axis.text = element_text(size = 18, color = "black"),
    axis.line = element_line(color = "black"),  # Change axis lines to black
    axis.ticks = element_line(color = "black"),  # Change axis ticks to black
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust the plot margins to add space around the plot
  )

# Plot records per latitude (with Kernel Estimation)
# Step 1: Calculate the histogram (number of records per 5-degree latitude bin)
hist_data_lat <- hist(Global_data$lat_5, breaks = seq(-90, 90, by = 5), plot = FALSE)

# Step 2: Find the maximum count
max_count_lat <- max(hist_data_lat$counts)
print(max_count_lat)
min_count_lat <- min(hist_data_lat$counts)
print(min_count_lat)
mean_count_lat <- mean(hist_data_lat$counts)
print(mean_count_lat)

# Perform the Anderson-Darling test
ad_test <- ad.test(hist_data_lat$counts)
print(ad_test) # A = 3.1139, p-value = 5.611e-08

# Perform the Dip Test for unimodality
dip_test <- dip.test(hist_data_lat$counts)
print(dip_test) # D = 0.074077, p-value < 2.2e-16

# Define axis limits for consistency
x_limits <- c(-90, 90) # Latitude range
y_limits <- c(0, 10000000) # Fixed y-axis limits

# Calculate relative annotation positions
x_annotate <- x_limits[1] + 0.29 * diff(x_limits)  # 25% from the left
y_annotate_ad <- y_limits[1] + 0.92 * diff(y_limits)  # 90% up from the bottom
y_annotate_dip <- y_limits[1] + 0.82 * diff(y_limits)  # 80% up from the bottom

# Create the plot
Lat_Num_Rec_Kernel_col <- ggplot(Global_data, aes(x = lat_5)) +
  geom_histogram(binwidth = 5, fill = "#1a80bb", color = "black", position = "identity") +
  geom_density(aes(y = after_stat(count * (8000000 / max(after_stat(count))))), bw = 10,
               alpha = 0.5, fill = "NA", color = "#ea801c", size = 1.5) +
  scale_y_continuous(limits = y_limits, breaks = seq(0, 10000000, by = 2000000), expand = c(0, 0), labels = scales::label_comma()) +
  scale_x_continuous(limits = x_limits, breaks = seq(-90, 90, by = 20), expand = c(0, 0)) +
  xlab("Latitude (degree)") + ylab("Number of Records") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = 24, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  # Add annotations for test results
  annotate("text", x = x_annotate, y = y_annotate_ad, 
           label = paste("AD-p:", formatC(ad_test$p.value, format = "e", digits = 2)), 
           size = 8, color = "black") +
  annotate("text", x = x_annotate, y = y_annotate_dip, 
           label = paste("DIP-p:", formatC(dip_test$p.value, format = "e", digits = 2)), 
           size = 8, color = "black")

# Save the plot
ggsave("Lat_Num_Rec_Kernel_col.tiff", plot = Lat_Num_Rec_Kernel_col, width = 8, height = 5, dpi = 600)

# Add kernel diversity on the graph
# Filter distinct values
distinct_data <- distinct(Global_data, scientificName, lat_5, .keep_all = TRUE)

# Remove any non-finite values
distinct_data <- distinct_data %>% filter(is.finite(lat_5))

# Count distinct scientific names in each 5-degree latitude bin
name_counts <- distinct_data %>%
  group_by(lat_5) %>%
  summarise(distinct_names = n_distinct(scientificName)) %>%
  ungroup()

# Perform the Anderson-Darling test
ad_test_1 <- ad.test(name_counts$distinct_names)
print(ad_test_1) # 0.6432, p-value = 0.08557

# Perform the Dip Test for unimodality
dip_test_1 <- dip.test(name_counts$distinct_names)
print(dip_test_1) # D = 0.058996, p-value = 0.4704

# Define axis limits for consistency
x_limits <- c(-90, 90) # Latitude range
y_limits <- c(0, 60000) # Fixed y-axis limits

# Calculate relative annotation positions
x_annotate <- x_limits[1] + 0.26 * diff(x_limits)  # 30% from the left
y_annotate_ad <- y_limits[1] + 0.92 * diff(y_limits)  # 95% up from the bottom
y_annotate_dip <- y_limits[1] + 0.82 * diff(y_limits)  # 85% up from the bottom

# Create the plot
Lat_Num_Spe_Kernel_col <- ggplot(distinct_data, aes(x = lat_5)) +
  geom_histogram(binwidth = 5, fill = "#1a80bb", color = "black", position = "identity") +
  geom_density(aes(y = after_stat(count * (40000 / max(after_stat(count))))), bw = 10,
               alpha = 0.5, fill = "NA", color = "#ea801c", size = 1.5) +
  scale_y_continuous(limits = y_limits, breaks = seq(0, 60000, by = 12000), expand = c(0, 0), labels = scales::label_comma()) +
  scale_x_continuous(limits = x_limits, breaks = seq(-90, 90, by = 20), expand = c(0, 0)) +
  xlab("Latitude (degree)") + ylab("Number of Species") +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = 24, color = "black"),     
    axis.text = element_text(size = 20, color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust the plot margins to add space around the plot
  ) +
  # Add annotations for test results
  annotate("text", x = x_annotate, y = y_annotate_ad, 
           label = paste("AD-p:", formatC(ad_test_1$p.value, format = "e", digits = 2)), 
           size = 8, color = "black") +
  annotate("text", x = x_annotate, y = y_annotate_dip, 
           label = paste("Dip-p:", formatC(dip_test_1$p.value, format = "e", digits = 2)), 
           size = 8, color = "black")

# Save the plot
ggsave("Lat_Num_Spe_Kernel_col.tiff", plot = Lat_Num_Spe_Kernel_col, width = 7, height = 5, dpi = 600)

# Plot records per depth with pylum
Dep_Num_Rec_Phy <- ggplot(Global_data, aes(x = dep_rnd, fill = phylum)) + 
  geom_bar(position="dodge") +
  scale_y_continuous(limits=c(0, 3500000), breaks=seq(0, 3500000, by=1000000), expand = c(0, 0)) +
  scale_x_continuous(limits=c(0,11000), breaks=seq(0, 11000, by=2000), expand = c(0, 0)) +
  xlab("Depth (m)") + ylab("Number of Records") +
  scale_fill_manual(values = phylum_colors) +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),    
    legend.text = element_text(size = 8),          # Adjust legend text size
    legend.title = element_text(size = 10),        # Adjust legend title size
    legend.key.size = unit(0.5, "cm")               # Adjust legend key size
  )
Print (Dep_Num_Rec_Phy)
ggsave("Dep_Num_Rec_Phy.tiff", plot = Dep_Num_Rec_Phy, width = 12, height = 6, dpi = 600)

# Plot records per depth (with a maximum of 4000)
Dep_Num_Rec_Phy_max4000 <- ggplot(Global_data, aes(x = dep_rnd, fill = phylum)) + 
  geom_bar(position="dodge") +
  scale_y_continuous(limits=c(0, 4000), breaks=seq(0, 4000, by=1000), expand = c(0, 0)) +
  scale_x_continuous(limits=c(0,11000), breaks=seq(0, 11000, by=2000), expand = c(0, 0)) +
  xlab("Depth (m)") + ylab("Number of Records") +
  scale_fill_manual(values = phylum_colors) +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),    
    legend.text = element_text(size = 8),          # Adjust legend text size
    legend.title = element_text(size = 10),        # Adjust legend title size
    legend.key.size = unit(0.5, "cm")               # Adjust legend key size
  )

print(Dep_Num_Rec_Phy_max4000)
ggsave("Dep_Num_Rec_Phy_max4000.tiff", plot = Dep_Num_Rec_Phy_max4000, width = 12, height = 6, dpi = 600)

# Plot records per depth (black and white and a maximum of 500000 (average of number of records per band))
Dep_Num_Rec <- ggplot(Global_data, aes(x = dep_rnd)) + 
  geom_bar(position="dodge") +
  scale_y_continuous(limits=c(0, 500000), breaks=seq(0, 500000, by=1000), expand = c(0, 0)) +
  scale_x_continuous(limits=c(0,11000), breaks=seq(0, 11000, by=2000), expand = c(0, 0)) +
  xlab("Depth (m)") + ylab("Number of Records") +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),    
    legend.text = element_text(size = 8),          # Adjust legend text size
    legend.title = element_text(size = 10),        # Adjust legend title size
    legend.key.size = unit(0.5, "cm")               # Adjust legend key size
  )
print (Dep_Num_Rec)
ggsave("Dep_Num_Rec_500000.tiff", plot = Dep_Num_Rec, width = 12, height = 6, dpi = 600)

# Plot records per depth (with Kernel Estimation)
# Step 1: Calculate the histogram
hist_data_dep <- hist(Global_data$dep_rnd, breaks=seq(0, 11000, by=100), plot=FALSE)

# Step 2: Find the maximum count
max_count_dep <- max(hist_data_dep$counts)
print(max_count_dep)
min_count_dep <- min(hist_data_dep$counts)
print(min_count_dep)
mean_count_dep <- mean(hist_data_dep$counts)
print(mean_count_dep)

# Perform the Anderson-Darling test
ad_test_2 <- ad.test(hist_data_dep$counts)
print(ad_test_2) # A = 39.72, Extremely Small p-value (< 2.2e-16) than significance level of 0.05. You can conclude with high confidence that your data does not follow a normal distribution.

# Perform the Dip Test for unimodality
dip_test_2 <- dip.test(hist_data_dep$counts)
print(dip_test_2) # D = 0.025735, p-value = 0.9133

# Define axis limits for consistency
x_limits <- c(0, 11000) # Latitude range
y_limits <- c(0, 40000000) # Fixed y-axis limits

# Calculate relative annotation positions
x_annotate <- x_limits[1] + 0.29 * diff(x_limits)  # 30% from the left
y_annotate_ad <- y_limits[1] + 0.92 * diff(y_limits)  # 95% up from the bottom
y_annotate_dip <- y_limits[1] + 0.82 * diff(y_limits)  # 85% up from the bottom

# Create the plot (the max is 38000000)
Dep_Num_Rec_Kernel_col <- ggplot(Global_data, aes(x = dep_rnd)) +
  geom_histogram(binwidth = 100, fill = "#1a80bb", color = "black", position = "identity") +
  geom_density(aes(y = after_stat(count * (38000000 / max(after_stat(count))))), bw = 20,
               alpha = 0.5, fill = "NA", color = "#ea801c", size = 1.5) +
  scale_y_continuous(limits=c(0, 40000000), breaks=seq(0, 40000000, by=8000000), expand = c(0, 0), labels = scales::label_comma()) +
  scale_x_continuous(limits=c(0, 11000), breaks=seq(0, 11000, by=2000), expand = c(0, 0), labels = scales::label_comma()) +
  xlab("Depth (m)") + ylab("Number of Records") +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = 24, color = "black"),     
    axis.text = element_text(size = 20, color = "black"),
    axis.line = element_line(color = "black"),  # Change axis lines to black
    axis.ticks = element_line(color = "black"),  # Change axis ticks to black
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust the plot margins to add space around the plot
  ) +
  # Add annotations for test results
  annotate("text", x = x_annotate, y = y_annotate_ad, 
           label = paste("AD-p:", formatC(ad_test_2$p.value, format = "e", digits = 2)), 
           size = 8, color = "black") +
  annotate("text", x = x_annotate, y = y_annotate_dip, 
           label = paste("Dip-p:", formatC(dip_test_2$p.value, format = "e", digits = 2)), 
           size = 8, color = "black")

ggsave("Dep_Num_Rec_Kernel_col_test.tiff", plot = Dep_Num_Rec_Kernel_col, width = 7, height = 5, dpi = 600)

# Plot species per depth
Dep_Num_Spe_Phy <- distinct(Global_data, scientificName, dep_rnd, .keep_all = TRUE) %>% 
  ggplot(aes(x = dep_rnd, fill = phylum)) + 
  geom_bar(position="dodge") +
  scale_y_continuous(limits=c(0, 30000), breaks=seq(0, 30000, by=5000), expand = c(0, 0)) +
  scale_x_continuous(limits=c(0,11000), breaks=seq(0, 11000, by=2000), expand = c(0, 0)) +
  xlab("Depth (m)") + ylab("Number of Species")+
  scale_fill_manual(values = phylum_colors) +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),    
    legend.text = element_text(size = 8),          # Adjust legend text size
    legend.title = element_text(size = 10),        # Adjust legend title size
    legend.key.size = unit(0.5, "cm")               # Adjust legend key size
  )

print(Dep_Num_Spe_Phy)
ggsave("Dep_Num_Spe_Phy.tiff", plot = Dep_Num_Spe_Phy, width = 12, height = 6, dpi = 600)

# Plot species per depth (with a maximum of 4000)
Dep_Num_Spe_Phy_max4000 <- distinct(Global_data, scientificName, dep_rnd, .keep_all = TRUE) %>% 
  ggplot(aes(x = dep_rnd, fill = phylum)) + 
  geom_bar(position="dodge") +
  scale_y_continuous(limits=c(0, 4000), breaks=seq(0, 4000, by=1000), expand = c(0, 0)) +
  scale_x_continuous(limits=c(0,11000), breaks=seq(0, 11000, by=2000), expand = c(0, 0)) +
  xlab("Depth (m)") + ylab("Number of Species")+
  scale_fill_manual(values = phylum_colors) +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),    
    legend.text = element_text(size = 8),          # Adjust legend text size
    legend.title = element_text(size = 10),        # Adjust legend title size
    legend.key.size = unit(0.5, "cm")               # Adjust legend key size
  )

print(Dep_Num_Spe_Phy_max4000)
ggsave("Dep_Num_Spe_Phy_max4000.tiff", plot = Dep_Num_Spe_Phy_max4000, width = 12, height = 6, dpi = 600)

# Add kernel density
# Filter distinct values
distinct_data_2 <- distinct(Global_data, scientificName, dep_rnd, .keep_all = TRUE)

# Remove any non-finite values
distinct_data_2 <- distinct_data_2 %>% filter(is.finite(dep_rnd))

# Count distinct scientific names in each 5-degree latitude bin
name_counts_2 <- distinct_data_2 %>%
  group_by(dep_rnd) %>%
  summarise(distinct_names_2 = n_distinct(scientificName)) %>%
  ungroup()

# Perform the Anderson-Darling test
ad_test_3 <- ad.test(name_counts_2$distinct_names_2)
print(ad_test_3) # A = 23.238, p-value < 2.2e-16

# Perform the Dip Test for unimodality
dip_test_3 <- dip.test(name_counts_2$distinct_names_2)
print(dip_test_3) # D = 0.036067, p-value = 0.4029

# Define axis limits for consistency
x_limits <- c(0, 11000) # Depth range
y_limits <- c(0, 70000) # Fixed y-axis limits

# Calculate relative annotation positions
x_annotate <- x_limits[1] + 0.26 * diff(x_limits)  # 30% from the left
y_annotate_ad <- y_limits[1] + 0.92 * diff(y_limits)  # 95% up from the bottom
y_annotate_dip <- y_limits[1] + 0.82 * diff(y_limits)  # 85% up from the bottom

#Create the plot
Dep_Num_Spe_Kernel_col <- ggplot(distinct_data_2, aes(x = dep_rnd)) +
  geom_histogram(binwidth = 100, fill = "#1a80bb", color = "black", position = "identity") +
  geom_density(aes(y = after_stat(density * nrow(distinct_data_2) * 100)), bw = 90,
               alpha = 0.5, fill = "NA", color = "#ea801c", size = 1.5) +
  scale_y_continuous(limits = c(0, 70000), breaks = seq(0, 70000, by = 10000), expand = c(0, 0), labels = scales::label_comma()) +
  scale_x_continuous(limits = c(0, 11000), breaks = seq(0, 11000, by = 2000), expand = c(0, 0), labels = scales::label_comma()) +
  xlab("Depth (m)") + ylab("Number of Species") +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = 24, color = "black"),     
    axis.text = element_text(size = 20, color = "black"),
    axis.line = element_line(color = "black"),  # Change axis lines to black
    axis.ticks = element_line(color = "black"),  # Change axis ticks to black
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust the plot margins to add space around the plot
  ) +
  # Add annotations for test results
  annotate("text", x = x_annotate, y = y_annotate_ad, 
           label = paste("AD-p:", formatC(ad_test_3$p.value, format = "e", digits = 2)), 
           size = 8, color = "black") +
  annotate("text", x = x_annotate, y = y_annotate_dip, 
           label = paste("Dip-p:", formatC(dip_test_3$p.value, format = "e", digits = 2)), 
           size = 8, color = "black")

ggsave("Dep_Num_Spe_Kernel_col.tiff", plot = Dep_Num_Spe_Kernel_col, width = 7, height = 5, dpi = 600)

# Plot (records) latitude against depth
ggplot(Global_data, aes(x = decimalLatitude, y = depth, colour = phylum)) +
  geom_jitter(alpha = 1) +
  xlab("Latitude (degree)") + ylab("Depth (m)") +
  scale_colour_manual(values = phylum_colors) +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank()
  )

# Plot (records) latitude / longitude
ggplot(Global_data, aes(x = decimalLongitude, y = decimalLatitude, colour = phylum)) +
  geom_jitter(alpha = 1) +
  xlab("Longitude (degree)") + ylab("Latitude (degree)") +
  scale_colour_manual(values = phylum_colors) +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank()
  )
  
  # Matrix Function
data <- Global_data
data <- data %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) # remove NAs
data$scientificName <- stri_encode(data$scientificName, "", "UTF-8") # re-mark encodings

# remember to format the dataset
#data$phylum <- as.factor(data$phylum)
#data <- mutate(data, lat_5 = ceiling(decimalLatitude/5) * 5)
#data <- mutate(data, dep_rnd = floor(depth/100) * 100)
#data <- mutate(data, dep_rat = factor(case_when(depth <= 200 ~ "shallow", depth>200 & depth<=500 ~ "meso", depth > 500 ~ "deep")))

# Define the function with a proper name
presence_absence_matrix <- function(data, sites.col, sp.col, keep.n = TRUE) {
  stopifnot(
    length(sites.col) == 1,
    length(sp.col) == 1,
    sites.col != sp.col,
    (sites.col %in% 1:ncol(data) || sites.col %in% names(data)),
    (sp.col %in% 1:ncol(data) || sp.col %in% names(data)),
    is.logical(keep.n)
  )
  
  if (is.character(sites.col)) sites.col <- match(sites.col, names(data))
  if (is.character(sp.col)) sp.col <- match(sp.col, names(data))
  
  presabs <- table(data[, c(sites.col, sp.col)])
  presabs <- as.data.frame.matrix(unclass(presabs))
  
  if (!keep.n) presabs[presabs > 1] <- 1
  
  presabs <- data.frame(row.names(presabs), presabs)
  names(presabs)[1] <- names(data)[sites.col]
  rownames(presabs) <- NULL
  
  return(presabs)
}

# Matrix
lat <- presence_absence_matrix(subset(data), "lat_5", "scientificName")%>%  
  column_to_rownames(var = "lat_5")

dep <- presence_absence_matrix(subset(data), "dep_rnd", "scientificName")%>%  
  column_to_rownames(var = "dep_rnd")

# Abundance 
abu_lat <- apply(lat, MARGIN = 1, sum)
#write.table(abu_lat,file="adundance.txt")

# Gamma richness: per 5-degree Latitudinal bands
sp_rich_lat <- specnumber(lat)

# Rarefaction
rarecurve(lat, step = 20, sample = 50, col = "black", cex = 0.6, label = TRUE, main = "Global")

# Calculate ES50 for each 5-degree latitude bin and store it as a table
ES50_table <- lat %>%
  rarefy(sample = 50, MARGIN = 1) %>%
  as.data.frame() %>%
  rownames_to_column("Latitude") %>%
  mutate(Latitude = as.numeric(Latitude))

# Rarefaction plots
lat %>%
  rarefy(sample = 50, MARGIN = 1) %>%
  as.data.frame() %>%
  rownames_to_column("Latitude") %>%
  mutate(Latitude = as.numeric(Latitude)) %>%
  ggplot(aes(x = Latitude, y = .)) +
  geom_col(fill = "lightgrey", color = "black") +  # Set bar color to light grey
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-90, 90), breaks = seq(-90, 90, by = 20), expand = c(0, 0)) +
  xlab("Latitude (degree)") + ylab("ES 50") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = 22, color = "black"),     
    axis.text = element_text(size = 18, color = "black"),
    axis.line = element_line(color = "black"),  # Change axis lines to black
    axis.ticks = element_line(color = "black"),  # Change axis ticks to black
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust the plot margins to add space around the plot
  )

# Add Kernel Density
# Rarefaction process Latitude
lat_rarefied <- lat %>%
  rarefy(sample = 50, MARGIN = 1) %>%
  as.data.frame() %>%
  rownames_to_column("Latitude") %>%
  mutate(Latitude = as.numeric(Latitude))

# Perform the Anderson-Darling test
ad_test_4 <- ad.test(lat_rarefied$.)
print(ad_test_4) # A = 1.3559, p-value = 0.001364

# Perform the Dip Test for unimodality
dip_test_4 <- dip.test(lat_rarefied$.)
print(dip_test_4) # D = 0.064924, p-value = 0.3042

# Define axis limits for consistency
x_limits <- c(-90, 90) # Latitude range
y_limits <- c(0, 70) # Fixed y-axis limits

# Calculate relative annotation positions
x_annotate <- x_limits[1] + 0.24 * diff(x_limits)  # 30% from the left
y_annotate_ad <- y_limits[1] + 0.92 * diff(y_limits)  # 95% up from the bottom
y_annotate_dip <- y_limits[1] + 0.82 * diff(y_limits)  # 85% up from the bottom

# Create the plot
Lat_ES50_Kernel_col_2 <- ggplot(lat_rarefied, aes(x = Latitude, y = .)) +
  geom_col(fill = "#1a88bb", color = "black", width = 5) +  # Set bar color to light grey
  geom_density(aes(y = after_stat(density * nrow(lat_rarefied) * 245)), bw = 30,
  #geom_density(data = lat_rarefied, aes(x = Latitude, y = ..density.. * 245 * nrow(lat_rarefied)), 
               alpha = 0.5, fill = "NA", color = "#ea801c", size = 1.5) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, by = 10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-90, 90), breaks = seq(-90, 90, by = 20), expand = c(0, 0)) +
  xlab("Latitude (degree)") + ylab("ES 50") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = 24, color = "black"),     
    axis.text = element_text(size = 20, color = "black"),
    axis.line = element_line(color = "black"),  # Change axis lines to black
    axis.ticks = element_line(color = "black"),  # Change axis ticks to black
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust the plot margins to add space around the plot
  ) +
  # Add annotations for test results
  annotate("text", x = x_annotate, y = y_annotate_ad, 
           label = paste("AD-p:", formatC(ad_test_4$p.value, format = "e", digits = 2)), 
           size = 8, color = "black") +
  annotate("text", x = x_annotate, y = y_annotate_dip, 
           label = paste("Dip-p:", formatC(dip_test_4$p.value, format = "e", digits = 2)), 
           size = 8, color = "black")

ggsave("Lat_ES50_Kernel_col.tiff", plot = Lat_ES50_Kernel_col_2, width = 7, height = 5, dpi = 600)


# Rarefaction process Depth
# Calculate row totals
row_totals <- rowSums(dep)

# Print row totals for verification
print(row_totals)

# Filter out rows with fewer than 50 samples
dep_filtered <- dep[row_totals >= 50, ]

# Check if any rows remain after filtering
if (nrow(dep_filtered) == 0) {
  stop("No rows with at least 50 samples available for rarefaction.")
}

# Perform rarefaction
dep_rarefied <- dep_filtered %>%
  rarefy(sample = 50, MARGIN = 1) %>%
  as.data.frame() %>%
  rownames_to_column("Depth") %>%
  mutate(Depth = as.numeric(Depth))

# Perform the Anderson-Darling test
ad_test_5 <- ad.test(dep_rarefied$.)
print(ad_test_5) # A = 4.0619, p-value = 3.369e-10

# Perform the Dip Test for unimodality
dip_test_5 <- dip.test(dep_rarefied$.)
print(dip_test_5) # D = 0.037655, p-value = 0.5876

# Define axis limits for consistency
x_limits <- c(0, 11000) # Latitude range
y_limits <- c(0, 70) # Fixed y-axis limits

# Calculate relative annotation positions
x_annotate <- x_limits[1] + 0.24 * diff(x_limits)  # 30% from the left
y_annotate_ad <- y_limits[1] + 0.92 * diff(y_limits)  # 95% up from the bottom
y_annotate_dip <- y_limits[1] + 0.82 * diff(y_limits)  # 85% up from the bottom

# Create the plot
Dep_ES50_Kernel_col_2 <- ggplot(dep_rarefied, aes(x = Depth, y = .)) +
  geom_col(fill = "#1a88bb", color = "black", width = 100) +  # Set bar color to light grey
  geom_density(aes(y = after_stat(density * nrow(dep_rarefied) * 4900)), bw = 300, # to create col_2
  #geom_density(data = dep_rarefied, aes(x = Depth, y = ..density.. * 4900 * nrow(dep_rarefied)), 
               alpha = 0.5, fill = "NA", color = "#ea801c", size = 1.5) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, by = 10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 11000), breaks = seq(0, 11000, by = 2000), expand = c(0, 0), labels = scales::label_comma()) +
  xlab("Depth (m)") + ylab("ES 50") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = 24, color = "black"),     
    axis.text = element_text(size = 20, color = "black"),
    axis.line = element_line(color = "black"),  # Change axis lines to black
    axis.ticks = element_line(color = "black"),  # Change axis ticks to black
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust the plot margins to add space around the plot
  ) +
  # Add annotations for test results
  annotate("text", x = x_annotate, y = y_annotate_ad, 
           label = paste("AD-p:", formatC(ad_test_5$p.value, format = "e", digits = 2)), 
           size = 8, color = "black") +
  annotate("text", x = x_annotate, y = y_annotate_dip, 
           label = paste("Dip-p:", formatC(dip_test_5$p.value, format = "e", digits = 2)), 
           size = 8, color = "black")

ggsave("Dep_ES50_Kernel_col_2.tiff", plot = Dep_ES50_Kernel_col_2, width = 7, height = 5, dpi = 600)


# Calculate the data gaps fro the centrl tropical areas
# Load data
load("Global_data.rda")  # assumes object = Global_data

library(ggplot2)
library(dplyr)
library(maps)

# Filter tropical region
tropical_data <- Global_data %>%
  filter(decimalLatitude >= -5 & decimalLatitude <= 5)

tropical_data <- tropical_data %>%
  filter(!is.na(year) & year > 1800 & year <= 2025)

year_counts <- tropical_data %>%
  group_by(year) %>%
  summarise(n_records = n())

# World map
world <- map_data("world")

# Plot
Tropical_Occurrence_Data_Year <- ggplot(year_counts, aes(x = year, y = n_records)) +
  geom_col(fill = "#1a80bb") +
  theme_bw() +
  labs(
    x = "Year",
    y = "Tropical occurrence records (n)"
  ) +
  scale_x_continuous(breaks = seq(min(year_counts$year, na.rm = TRUE),
                                  max(year_counts$year, na.rm = TRUE),
                                  by = 20)) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

ggsave("Tropical_Occurrence_Data_Year.tiff", plot = Tropical_Occurrence_Data_Year, width = 7, height = 5, dpi = 600)

Global_data_clean <- Global_data %>%
  filter(!is.na(year) & year >= 2013 & year <= 2023)

tropical_subset <- Global_data_clean %>%
  filter(decimalLatitude >= -5 & decimalLatitude <= 5)

total_global <- nrow(Global_data_clean)
total_tropical <- nrow(tropical_subset)

percentage_tropical <- (total_tropical / total_global) * 100

percentage_tropical

yearly_percent <- Global_data %>%
  filter(!is.na(year) & year >= 2013 & year <= 2023) %>%
  mutate(region = ifelse(decimalLatitude >= -5 & decimalLatitude <= 5, "Tropical", "Global")) %>%
  group_by(year) %>%
  summarise(
    total = n(),
    tropical = sum(decimalLatitude >= -5 & decimalLatitude <= 5),
    percent_tropical = (tropical / total) * 100
  )

Tropical_Global_Percent_Occurrence_Data_Year <- ggplot(yearly_percent, aes(x = year, y = percent_tropical)) +
  geom_col(fill = "#1a80bb") +
  geom_text(aes(label = round(percent_tropical, 1)),
            vjust = -0.5, size = 5) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  ) +
  labs(
    x = "Year",
    y = "Tropical occurrence records (%)"
  ) +
  expand_limits(y = max(yearly_percent$percent_tropical, na.rm = TRUE) * 1.1)

ggsave("Tropical_Global_Percent_Occurrence_Data_Year.tiff", plot = Tropical_Global_Percent_Occurrence_Data_Year, width = 7, height = 5, dpi = 600)

