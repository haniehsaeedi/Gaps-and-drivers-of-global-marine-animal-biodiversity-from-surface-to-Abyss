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
library(purrr)          # For functional programming
library(diptest)  # For Hartigan's Dip Test

# load data
load("Global_data_Cnidaria.rda") #to load the .rda

#Summary
summary(Global_data_Cnidaria)
Global_data_Cnidaria_table <- summary(Global_data_Cnidaria)
subset(Global_data_Cnidaria, select = scientificName) %>% n_distinct() # get number of species

# Count the total number of scientific names per latitude bin
total_counts_per_lat <- Global_data_Cnidaria %>%
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
unique_counts_per_lat <- Global_data_Cnidaria %>%
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

# Count the total number of scientific names per depth bin
total_counts_per_dep <- Global_data_Cnidaria %>%
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
unique_counts_per_dep <- Global_data_Cnidaria %>%
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
..............................................
# Plot records per latitude (with Kernel Estimation)
# Step 1: Calculate the histogram
hist_data_lat <- hist(Global_data_Cnidaria$lat_5, breaks=seq(-90, 90, by=5), plot=FALSE)

# Step 2: Find the maximum count
print(c(max_count_lat = max(hist_data_lat$counts),
        min_count_lat = min(hist_data_lat$counts),
        mean_count_lat = mean(hist_data_lat$counts)))

# Perform the Anderson-Darling test
ad_test <- ad.test(hist_data_lat$counts)
print(ad_test) # A = 2.5414, p-value = 1.489e-06

# Perform the Dip Test for unimodality
dip_test <- dip.test(hist_data_lat$counts)
print(dip_test) # D = 0.060672, p-value = 0.3745

# Define axis limits for consistency
x_limits <- c(-90, 90) # Latitude range
y_limits <- c(0, 500000) # Fixed y-axis limits

# Calculate relative annotation positions
x_annotate <- x_limits[1] + 0.26 * diff(x_limits)  # 30% from the left
y_annotate_ad <- y_limits[1] + 0.92 * diff(y_limits)  # 95% up from the bottom
y_annotate_dip <- y_limits[1] + 0.82 * diff(y_limits)  # 85% up from the bottom

# Create the plot
Lat_Num_Rec_Kernel_Cnidaria_col <- ggplot(Global_data_Cnidaria, aes(x = lat_5)) +
  geom_histogram(binwidth = 5, fill = "#1a80bb", color = "black", position = "identity") +
  geom_density(aes(y = after_stat(count * (380000 / max(after_stat(count))))), bw = 10, 
               alpha = 0.5, fill = "NA", color = "#ea801c", size = 1.5) +
  scale_y_continuous(limits=c(0, 500000), breaks=seq(0, 500000, by=100000), expand = c(0, 0), labels = scales::label_comma()) +
  scale_x_continuous(limits=c(-90, 90), breaks=seq(-90, 90, by=20), expand = c(0, 0)) +
  xlab("Latitude (degree)") + ylab("Number of Records") +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    text = element_text(size = 20),           
    axis.title = element_text(size = 24, color = "black"),     
    axis.text = element_text(size = 20, color = "black"),
    axis.line = element_line(color = "black"),  # Change axis lines to black
    axis.ticks = element_line(color = "black"),  # Change axis ticks to black
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust the plot margins to add space around the plot
  ) +
  # Add annotations for test results
  annotate("text", x = x_annotate, y = y_annotate_ad, 
           label = paste("AD/p:", formatC(ad_test$p.value, format = "e", digits = 2)), 
           size = 8, color = "black") +
  annotate("text", x = x_annotate, y = y_annotate_dip, 
           label = paste("Dip/p:", formatC(dip_test$p.value, format = "e", digits = 2)), 
           size = 8, color = "black")

ggsave("Lat_Num_Rec_Kernel_Cnidaria_col.tiff", plot = Lat_Num_Rec_Kernel_Cnidaria_col, width = 7, height = 5, dpi = 600)

# Plot species numbers per latitude (with Kernel Estimation)
# Filter distinct values
distinct_data <- distinct(Global_data_Cnidaria, scientificName, lat_5, .keep_all = TRUE)

# Remove any non-finite values
distinct_data <- distinct_data %>% filter(is.finite(lat_5))

# Count distinct scientific names in each 5-degree latitude bin
name_counts <- distinct_data %>%
  group_by(lat_5) %>%
  summarise(distinct_names = n_distinct(scientificName)) %>%
  ungroup()

# Perform the Anderson-Darling test
ad_test_1 <- ad.test(name_counts$distinct_names)
print(ad_test_1) # A = 0.57407, p-value = 0.1256

# Perform the Dip Test for unimodality
dip_test_1 <- dip.test(name_counts$distinct_names)
print(dip_test_1) # D = 0.060182, p-value = 0.4342

# Define axis limits for consistency
x_limits <- c(-90, 90) # Latitude range
y_limits <- c(0, 4500) # Fixed y-axis limits

# Calculate relative annotation positions
x_annotate <- x_limits[1] + 0.25 * diff(x_limits)  # 30% from the left
y_annotate_ad <- y_limits[1] + 0.92 * diff(y_limits)  # 95% up from the bottom
y_annotate_dip <- y_limits[1] + 0.82 * diff(y_limits)  # 85% up from the bottom

# Create the plot
Lat_Num_Spe_Kernel_Cnidaria_col <- ggplot(distinct_data, aes(x = lat_5)) +
  geom_histogram(binwidth = 5, fill = "#1a80bb", color = "black", position = "identity") +
  geom_density(aes(y = after_stat(count * (3500 / max(after_stat(count))))), bw = 10,
               alpha = 0.5, fill = "NA", color = "#ea801c", size = 1.5) +
  scale_y_continuous(limits=c(0, 4500), breaks=seq(0, 4500, by=900), expand = c(0, 0), labels = scales::label_comma()) +
  scale_x_continuous(limits=c(-90, 90), breaks=seq(-90, 90, by=20), expand = c(0, 0)) +
  xlab("Latitude (degree)") + ylab("Number of Species") +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    text = element_text(size = 20),           
    axis.title = element_text(size = 24, color = "black"),     
    axis.text = element_text(size = 20, color = "black"),
    axis.line = element_line(color = "black"),  # Change axis lines to black
    axis.ticks = element_line(color = "black"),  # Change axis ticks to black
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust the plot margins to add space around the plot
  ) +
  # Add annotations for test results
  annotate("text", x = x_annotate, y = y_annotate_ad, 
           label = paste("AD-p:", formatC(ad_test_1$p.value, format = "e", digits = 2)), 
           size = 8, color = "black") +
  annotate("text", x = x_annotate, y = y_annotate_dip, 
           label = paste("Dip-p:", formatC(dip_test_1$p.value, format = "e", digits = 2)), 
           size = 8, color = "black")


ggsave("Lat_Num_Spe_Kernel_Cnidaria_col.tiff", plot = Lat_Num_Spe_Kernel_Cnidaria_col, width = 7, height = 5, dpi = 600)

...................................
# Plot records per depth (with Kernel Estimation)
# Step 1: Calculate the histogram
hist_data_dep <- hist(Global_data_Cnidaria$dep_rnd, breaks=seq(0, 11000, by=100), plot=FALSE)

# Step 2: Find the maximum count
print(c(max_count_dep = max(hist_data_dep$counts),
        min_count_dep = min(hist_data_dep$counts),
        mean_count_dep = mean(hist_data_dep$counts)))

# Perform the Anderson-Darling test
ad_test_2 <- ad.test(hist_data_dep$counts)
print(ad_test_2) # A = 34.829, p-value < 2.2e-16

# Perform the Dip Test for unimodality
dip_test_2 <- dip.test(hist_data_dep$counts)
print(dip_test_2) # D = 0.027273, p-value = 0.8526

# Define axis limits for consistency
x_limits <- c(0, 11000) # Latitude range
y_limits <- c(0, 250000) # Fixed y-axis limits

# Calculate relative annotation positions
x_annotate <- x_limits[1] + 0.26 * diff(x_limits)  # 30% from the left
y_annotate_ad <- y_limits[1] + 0.92 * diff(y_limits)  # 95% up from the bottom
y_annotate_dip <- y_limits[1] + 0.82 * diff(y_limits)  # 85% up from the bottom

# Create the plot 
Dep_Num_Rec_Kernel_Cnidaria_col <- ggplot(Global_data_Cnidaria, aes(x = dep_rnd)) +
  geom_histogram(binwidth = 100, fill = "#1a80bb", color = "black", position = "identity") +
  geom_density(aes(y = after_stat(count * (190000 / max(after_stat(count))))), bw = 350,
               alpha = 0.5, fill = "NA", color = "#ea801c", size = 1.5) +
  scale_y_continuous(limits=c(0, 250000), breaks=seq(0, 250000, by=50000), expand = c(0, 0), labels = scales::label_comma()) +
  scale_x_continuous(limits=c(0, 11000), breaks=seq(0, 11000, by=2000), expand = c(0, 0), labels = scales::label_comma()) +
  xlab("Depth (m)") + ylab("Number of Records") +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    text = element_text(size = 20),           
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

ggsave("Dep_Num_Rec_Kernel_Cnidaria_col.tiff", plot = Dep_Num_Rec_Kernel_Cnidaria_col, width = 7, height = 5, dpi = 600)

.........................................
# Plot species per depth (with Kernel Estimation)
# Filter distinct values
distinct_data_2 <- distinct(Global_data_Cnidaria, scientificName, dep_rnd, .keep_all = TRUE)

# Remove any non-finite values
distinct_data_2 <- distinct_data_2 %>% filter(is.finite(dep_rnd))

# Count distinct scientific names in each 100 m depth bin
name_counts_2 <- distinct_data_2 %>%
  group_by(dep_rnd) %>%
  summarise(distinct_names_2 = n_distinct(scientificName)) %>%
  ungroup()

# Step 2: Find the maximum count
print(c(max_count_dep = max(name_counts_2$distinct_names_2),
        min_count_dep = min(name_counts_2$distinct_names_2),
        mean_count_dep = mean(name_counts_2$distinct_names_2)))

# Perform the Anderson-Darling test
ad_test_3 <- ad.test(name_counts_2$distinct_names_2)
print(ad_test_3) # A = 14.689, p-value < 2.2e-16

# Perform the Dip Test for unimodality
dip_test_3 <- dip.test(name_counts_2$distinct_names_2)
print(dip_test_3) # D = 0.04023, p-value = 0.387

# Define axis limits for consistency
x_limits <- c(0, 11000) # Latitude range
y_limits <- c(0, 6000) # Fixed y-axis limits

# Calculate relative annotation positions
x_annotate <- x_limits[1] + 0.25 * diff(x_limits)  # 30% from the left
y_annotate_ad <- y_limits[1] + 0.92 * diff(y_limits)  # 95% up from the bottom
y_annotate_dip <- y_limits[1] + 0.82 * diff(y_limits)  # 85% up from the bott

#Create the plot
Dep_Num_Spe_Kernel_Cnidaria_col <- ggplot(distinct_data_2, aes(x = dep_rnd)) +
  geom_histogram(binwidth = 100, fill = "#1a88bb", color = "black", position = "identity") +
  geom_density(aes(y = after_stat(density * nrow(distinct_data_2) * 130)), bw = 250, 
               alpha = 0.5, fill = "NA", color = "#ea801c", size = 1.5) +
  scale_y_continuous(limits = c(0, 6000), breaks = seq(0, 6000, by = 1200), expand = c(0, 0), labels = scales::label_comma()) +
  scale_x_continuous(limits = c(0, 11000), breaks = seq(0, 11000, by = 2000), expand = c(0, 0), labels = scales::label_comma()) +
  xlab("Depth (m)") + ylab("Number of Species") +
  theme_bw() +   
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    text = element_text(size = 20),           
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

ggsave("Dep_Num_Spe_Kernel_Cnidaria_col.tiff", plot = Dep_Num_Spe_Kernel_Cnidaria_col, width = 7, height = 5, dpi = 600)

.........................................
# Matrix Function
data <- Global_data_Cnidaria 
data <- data %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) # remove NAs
data$scientificName <- stri_encode(data$scientificName, "", "UTF-8") # re-mark encodings
remove(Global_data_Cnidaria)

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
tiff("Rarefaction_Cnidaria.tiff", width = 6, height = 6, units = "in", res = 600)
rarecurve(lat, step = 20, sample = 50, col = "black", cex = 0.6, label = TRUE, main = "Cnidaria")
dev.off()

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
    text = element_text(size = 20),           
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
print(ad_test_4) # A = 0.87585, p-value = 0.02211

# Perform the Dip Test for unimodality
dip_test_4 <- dip.test(lat_rarefied$.)
print(dip_test_4) # D = 0.042485, p-value = 0.9407

# Define axis limits for consistency
x_limits <- c(-90, 90) # Latitude range
y_limits <- c(0, 70) # Fixed y-axis limits

# Calculate relative annotation positions
x_annotate <- x_limits[1] + 0.23 * diff(x_limits)  # 30% from the left
y_annotate_ad <- y_limits[1] + 0.92 * diff(y_limits)  # 95% up from the bottom
y_annotate_dip <- y_limits[1] + 0.82 * diff(y_limits)  # 85% up from the bottom

# Create the plot
Lat_ES50_Kernel_Cnidaria_col <- ggplot(lat_rarefied, aes(x = Latitude, y = .)) +
  geom_col(fill = "#1a80bb", color = "black", width = 5) +  # Set bar color to light grey
  geom_density(aes(y = after_stat(density * nrow(lat_rarefied) * 250)), bw = 20,
  #geom_density(data = lat_rarefied, aes(x = Latitude, y = ..density.. * 245 * nrow(lat_rarefied)), 
               alpha = 0.5, fill = "NA", color = "#ea801c", size = 1.5) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, by = 10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-90, 90), breaks = seq(-90, 90, by = 20), expand = c(0, 0)) +
  xlab("Latitude (degree)") + ylab("ES 50") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    text = element_text(size = 20),           
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

ggsave("Lat_ES50_Kernel_Cnidaria_col.tiff", plot = Lat_ES50_Kernel_Cnidaria_col, width = 7, height = 5, dpi = 600)

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
print(ad_test_5) # A = 1.4823, p-value = 0.000716

# Perform the Dip Test for unimodality
dip_test_5 <- dip.test(dep_rarefied$.)
print(dip_test_5) # D = 0.031946, p-value = 0.9631

# Define axis limits for consistency
x_limits <- c(0, 11000) # Latitude range
y_limits <- c(0, 70) # Fixed y-axis limits

# Calculate relative annotation positions
x_annotate <- x_limits[1] + 0.23 * diff(x_limits)  # 30% from the left
y_annotate_ad <- y_limits[1] + 0.92 * diff(y_limits)  # 95% up from the bottom
y_annotate_dip <- y_limits[1] + 0.82 * diff(y_limits)  # 85% up from the bottom

# Create the plot
Dep_ES50_Kernel_Cnidaria_col <- ggplot(dep_rarefied, aes(x = Depth, y = .)) +
  geom_col(fill = "#1a80bb", color = "black", width = 100) +  # Set bar color to light grey
  geom_density(aes(y = after_stat(density * nrow(dep_rarefied) * 4500)), bw = 300, # to create col_2
  #geom_density(data = dep_rarefied, aes(x = Depth, y = ..density.. * 4900 * nrow(dep_rarefied)), 
               alpha = 0.5, fill = "NA", color = "#ea801c", size = 1.5) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, by = 10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 11000), breaks = seq(0, 11000, by = 2000), expand = c(0, 0), labels = scales::label_comma()) +
  xlab("Depth (m)") + ylab("ES 50") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    text = element_text(size = 20),           
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

ggsave("Dep_ES50_Kernel_Cnidaria_col.tiff", plot = Dep_ES50_Kernel_Cnidaria_col, width = 7, height = 5, dpi = 600)

# Define a 5-color gradient
#custom_colors <- c("#aadee4", "#009999", "#ffdd3c", "#fdae61", "#d73027")

custom_colors <- c( "#b3e1e8", "#66cccc", "olivedrab3", "orange", "orangered3")

# Define a function to create a presence-absence matrix
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

# read the shapefiles
sf_land <- st_read("ne_110m_land.shp")
sf_cat <- st_read("hexgrid3.shp")
sf_data <- data[, c("scientificName", "decimalLatitude", "decimalLongitude")]
sf_data <- sf_data %>% st_as_sf(coords = c('decimalLongitude','decimalLatitude'))
st_crs(sf_data) = 4326

# Define the spatial function with adjusted rarefaction calculation
spatial_fun <- function(sf_land, sf_cat, sf_data)
  
# Step 1: Ensure data is in WGS84
sf_land <- st_transform(sf_land, crs = 4326)
sf_cat <- st_transform(sf_cat, crs = 4326)
sf_data <- st_transform(sf_data, crs = 4326)

# Step 2: Perform spatial joins and summarizations
spatial_data <- st_join(sf_cat, sf_data, join = st_contains) %>% 
  group_by(ID) %>% 
  dplyr::summarize(
    Num_Records = length(scientificName[!is.na(scientificName)]), 
    Num_Species = n_distinct(scientificName[!is.na(scientificName)])
  )

spatial_data <- spatial_data %>% as.data.frame()
spatial_data <- spatial_data[rowSums(spatial_data[, c("Num_Records", "Num_Species")]) != 0, ]


write.csv(spatial_data, "spatial_data.csv")
write.table(spatial_data,file="spatial_data.txt")
library(writexl)
write_xlsx(spatial_data, path = "spatial_data.xlsx")

# Step 3: Create presence-absence matrix
spatial_data2 <- st_join(sf_data, sf_cat, join = st_intersects) %>% 
  as.data.frame() %>% 
  presence_absence_matrix("ID", "scientificName")

# Step 4: Filter out rows with less than 50 samples
spatial_data2 <- spatial_data2[rowSums(spatial_data2[,-1]) >= 50, ]

# Step 5: Adjust rarefaction sample size to a reasonable threshold
min_sample_size <- min(rowSums(spatial_data2[,-1]))

# Step 6: Perform rarefaction with adjusted sample size
spatial_data_r <- rarefy(spatial_data2[,-1], sample = min_sample_size, MARGIN = 1) %>% 
  as.data.frame()
colnames(spatial_data_r) <- "ES50"
spatial_data_r$ID <- spatial_data2$ID

# Step 7: Merge data
spatial_data <- spatial_data %>% 
  merge(spatial_data_r, by = "ID", all.x = TRUE) %>% 
  st_as_sf()

# Calculate other Indices
calculate_indices <- function(spatial_data2) {
  chao1 <- estimateR(spatial_data2[,-1])["S.chao1", ]
  ace <- estimateR(spatial_data2[,-1])["S.ACE", ]
  weight_adjustment <- diversity(spatial_data2[,-1], index = "invsimpson")
  data.frame(ID = spatial_data2$ID, Chao1 = chao1, ACE = ace, Weighted = weight_adjustment)
}

add_indices_to_data <- function(spatial_data, estimate_results) {
  merge(spatial_data, estimate_results, by = "ID", all.x = TRUE)
}

# Add the Indices to the data
estimate_results <- calculate_indices(spatial_data2)
spatial_data <- add_indices_to_data(spatial_data, estimate_results)

summary(spatial_data)

#Chao2 needs a new dataframe using only presence absence data 

# Convert the count matrix to a binary presence/absence matrix
spatial_data2_binary <- spatial_data2 %>%
  mutate(across(-ID, ~ ifelse(. > 0, 1, 0)))

# Initialize a vector to store Chao2 values
chao2_values <- numeric(nrow(spatial_data2_binary))

# Loop through each row in the binary matrix
for (i in 1:nrow(spatial_data2_binary)) {
  # Select the single sample (excluding the ID column)
  single_sample <- as.matrix(spatial_data2_binary[i, -1])
  
  # Check if the single sample has more than one species presence
  if (sum(single_sample) > 1) {
    # Calculate Chao2 using specpool
    spec_pool <- specpool(single_sample)
    if (!is.null(spec_pool$chao) && length(spec_pool$chao) > 0) {
      chao2_values[i] <- spec_pool$chao[1]
    } else {
      chao2_values[i] <- NA  # Assign NA if Chao2 is not calculated
    }
  } else {
    chao2_values[i] <- NA  # Assign NA if there is not more than one species presence
  }
}

# Create a data frame to store Chao2 results along with the ID
chao2_results <- data.frame(
  ID = spatial_data2_binary$ID,  # Include the ID column from the original data
  Chao2 = chao2_values
)

# Add the Indices to the data for Chao2
add_indices_to_data_2 <- function(spatial_data, chao2_results) {
  merge(spatial_data, chao2_results, by = "ID", all.x = TRUE)
}
spatial_data_Chao2 <- add_indices_to_data_2(spatial_data, chao2_results)

summary(spatial_data_Chao2)

# Step 8: Create bins for each variable with labels showing the ranges
breaks_Num_Records <- c(1, 100000, 200000, 300000, 400000, 413067)
labels_Num_Records <- c("1-100000", "100000-200000", "200000-300000", "300000-400000", "400000-413067")
spatial_data$Num_Records_Binned <- cut(spatial_data$Num_Records, 
                                       breaks = breaks_Num_Records, 
                                       labels = labels_Num_Records, 
                                       include.lowest = TRUE)

breaks_Num_Species <- c(1, 100, 300, 500, 700, 1772)
labels_Num_Species <- c("1-100", "100-300", "300-500", "500-700", "700-1772")
spatial_data$Num_Species_Binned <- cut(spatial_data$Num_Species, 
                                       breaks = breaks_Num_Species, 
                                       labels = labels_Num_Species, 
                                       include.lowest = TRUE)

breaks_ES50 <- c(1, 10, 20, 30, 40, 50)
labels_ES50 <- c("1-10", "10-20", "20-30", "30-40", "40-50")
spatial_data$ES50_Binned <- cut(spatial_data$ES50, 
                                breaks = breaks_ES50, 
                                labels = labels_ES50, 
                                include.lowest = TRUE)

breaks_Chao1 <- c(1, 500, 1000, 1500, 2000, 2343)
labels_Chao1 <- c("1-500", "500-1000", "1000-1500", "1500-2000", "2000-2343")
spatial_data$Chao1_Binned <- cut(spatial_data$Chao1, 
                                 breaks = breaks_Chao1, 
                                 labels = labels_Chao1, 
                                 include.lowest = TRUE)

breaks_ACE <- c(1, 500, 1000, 1500, 2000, 2893)
labels_ACE <- c("1-500", "500-1000", "1000-1500", "1500-2000", "2000-2893")
spatial_data$ACE_Binned <- cut(spatial_data$ACE, 
                               breaks = breaks_ACE, 
                               labels = labels_ACE, 
                               include.lowest = TRUE)

breaks_Weighted  <- c(1, 50, 100, 150, 200, 361)
labels_Weighted  <- c("1-50", "50-100", "100-150", "150-200", "200-361")
spatial_data$Weighted_Binned <- cut(spatial_data$Weighted, 
                                    breaks = breaks_Weighted, 
                                    labels = labels_Weighted, 
                                    include.lowest = TRUE)

breaks_Chao2 <- c(1, 300, 600, 900, 1000, 1772)
labels_Chao2 <- c("1-300", "300-600", "600-900", "900-1000", "1000-1772")
spatial_data_Chao2$Chao2_Binned <- cut(spatial_data_Chao2$Chao2, 
                                 breaks = breaks_Chao2, 
                                 labels = labels_Chao2, 
                                 include.lowest = TRUE)
# Step 9: Generate plots
abbreviated_labels_Num_Records <- c("1-100k", "100-200k", "200-300k", "300-400K", "400-413K")

plot1 <- ggplot() +
  geom_sf(data = spatial_data, aes(fill = Num_Records_Binned)) +
  scale_fill_manual(
    values = custom_colors,
    drop = FALSE,
    name = "Number of Records",
    labels = abbreviated_labels_Num_Records  # Use abbreviated labels
  ) +
  geom_sf(data = sf_land, fill = "black") +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme_bw() +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.key.width = unit(0.5, "cm"),  # Width of legend keys
    legend.key.height = unit(0.1, "cm"),  # Height of legend keys
    legend.box.spacing = unit(0.02, "cm"),  # Space around the legend box
    legend.text = element_text(size = 16),  # Font size for legend text
    legend.title = element_text(size = 16)  # Font size for legend title
  ) +
  guides(
    fill = guide_legend(
      keywidth = 0.5,  # Width of legend keys
      keyheight = 0.5,  # Height of legend keys
      title.position = "top",  # Position of the legend title
      label.position = "bottom"  # Position of the legend labels
    )
  )

plot2 <- ggplot() +
  geom_sf(data = spatial_data, aes(fill = Num_Species_Binned)) +
  scale_fill_manual(values = custom_colors, name = "Number of Species") +
  geom_sf(data = sf_land, fill = "black") +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme_bw() +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.key.width = unit(0.5, "cm"),  # Width of legend keys
    legend.key.height = unit(0.1, "cm"),  # Height of legend keys
    legend.box.spacing = unit(0.02, "cm"),  # Space around the legend box
    legend.text = element_text(size = 16),  # Font size for legend text
    legend.title = element_text(size = 16)  # Font size for legend title
  ) +
  guides(
    fill = guide_legend(
      keywidth = 0.5,  # Width of legend keys
      keyheight = 0.5,  # Height of legend keys
      title.position = "top",  # Position of the legend title
      label.position = "bottom"  # Position of the legend labels
    )
  )

# Filter out NA values in the data
spatial_data_filtered <- spatial_data %>% 
  dplyr::filter(!is.na(ES50_Binned))

plot3 <- ggplot() +
  geom_sf(data = spatial_data_filtered, aes(fill = ES50_Binned)) + # Use filtered data here
  scale_fill_manual(values = custom_colors, name = "ES50", na.value = "white") + 
  geom_sf(data = sf_land, fill = "black") +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme_bw() +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.key.width = unit(0.5, "cm"),  # Width of legend keys
    legend.key.height = unit(0.1, "cm"),  # Height of legend keys
    legend.box.spacing = unit(0.02, "cm"),  # Space around the legend box
    legend.text = element_text(size = 16),  # Font size for legend text
    legend.title = element_text(size = 16)  # Font size for legend title
  ) +
  guides(
    fill = guide_legend(
      keywidth = 0.5,  # Width of legend keys
      keyheight = 0.5,  # Height of legend keys
      title.position = "top",  # Position of the legend title
      label.position = "bottom"  # Position of the legend labels
    )
  )


spatial_data_filtered_2 <- spatial_data %>% 
  dplyr::filter(!is.na(Chao1_Binned))

plot4 <- ggplot() +
  geom_sf(data = spatial_data_filtered_2, aes(fill = Chao1_Binned)) + # Use filtered data here
  scale_fill_manual(values = custom_colors, name = "Chao1", na.value = "white", drop = FALSE) + 
  geom_sf(data = sf_land, fill = "black") +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme_bw() +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.key.width = unit(0.5, "cm"),  # Width of legend keys
    legend.key.height = unit(0.1, "cm"),  # Height of legend keys
    legend.box.spacing = unit(0.02, "cm"),  # Space around the legend box
    legend.text = element_text(size = 16),  # Font size for legend text
    legend.title = element_text(size = 16)  # Font size for legend title
  ) +
  guides(
    fill = guide_legend(
      keywidth = 0.5,  # Width of legend keys
      keyheight = 0.5,  # Height of legend keys
      title.position = "top",  # Position of the legend title
      label.position = "bottom"  # Position of the legend labels
    )
  )

spatial_data_filtered_3 <- spatial_data %>% 
  dplyr::filter(!is.na(ACE_Binned))

plot5 <- ggplot() +
  geom_sf(data = spatial_data_filtered_3, aes(fill = ACE_Binned)) + # Use filtered data here
  scale_fill_manual(values = custom_colors, name = "ACE", na.value = "white", drop = FALSE) + 
  geom_sf(data = sf_land, fill = "black") +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme_bw() +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.key.width = unit(0.5, "cm"),  # Width of legend keys
    legend.key.height = unit(0.1, "cm"),  # Height of legend keys
    legend.box.spacing = unit(0.02, "cm"),  # Space around the legend box
    legend.text = element_text(size = 16),  # Font size for legend text
    legend.title = element_text(size = 16)  # Font size for legend title
  ) +
  guides(
    fill = guide_legend(
      keywidth = 0.5,  # Width of legend keys
      keyheight = 0.5,  # Height of legend keys
      title.position = "top",  # Position of the legend title
      label.position = "bottom"  # Position of the legend labels
    )
  )


spatial_data_filtered_4 <- spatial_data %>% 
  dplyr::filter(!is.na(Weighted_Binned))

plot6 <- ggplot() +
  geom_sf(data = spatial_data_filtered_4, aes(fill = Weighted_Binned)) + # Use filtered data here
  scale_fill_manual(values = custom_colors, name = "Inverse Simpson", na.value = "white", drop = FALSE) + 
  geom_sf(data = sf_land, fill = "black") +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme_bw() +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.key.width = unit(0.5, "cm"),  # Width of legend keys
    legend.key.height = unit(0.1, "cm"),  # Height of legend keys
    legend.box.spacing = unit(0.02, "cm"),  # Space around the legend box
    legend.text = element_text(size = 16),  # Font size for legend text
    legend.title = element_text(size = 16)  # Font size for legend title
  ) +
  guides(
    fill = guide_legend(
      keywidth = 0.5,  # Width of legend keys
      keyheight = 0.5,  # Height of legend keys
      title.position = "top",  # Position of the legend title
      label.position = "bottom"  # Position of the legend labels
    )
  )

spatial_data_filtered_5 <- spatial_data_Chao2 %>% 
  dplyr::filter(!is.na(Chao2_Binned))

plot7 <- ggplot() +
  geom_sf(data = spatial_data_filtered_5, aes(fill = Chao2_Binned)) + # Use filtered data here
  scale_fill_manual(values = custom_colors, name = "Chao2", na.value = "white", drop = FALSE) + 
  geom_sf(data = sf_land, fill = "black") +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme_bw() +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.key.width = unit(0.5, "cm"),  # Width of legend keys
    legend.key.height = unit(0.1, "cm"),  # Height of legend keys
    legend.box.spacing = unit(0.02, "cm"),  # Space around the legend box
    legend.text = element_text(size = 16),  # Font size for legend text
    legend.title = element_text(size = 16)  # Font size for legend title
  ) +
  guides(
    fill = guide_legend(
      keywidth = 0.5,  # Width of legend keys
      keyheight = 0.5,  # Height of legend keys
      title.position = "top",  # Position of the legend title
      label.position = "bottom"  # Position of the legend labels
    )
  )

# Display the plots
plot1 # Plot number of occurrences
ggsave("Num_Rec_Cnidaria.tiff", plot = plot1, width = 6, height = 4, dpi = 600)
plot2 # Plot number of species
ggsave("Num_Spe_Cnidaria.tiff", plot = plot2, width = 6, height = 4, dpi = 600)
plot3 # Plot ES50
ggsave("ES50_Cnidaria.tiff", plot = plot3, width = 6, height = 4, dpi = 600)
plot4 # Plot Chao1
ggsave("Chao1_Cnidaria.tiff", plot = plot4, width = 6, height = 4, dpi = 600)
plot5 # Plot ACE
ggsave("ACE_Cnidaria.tiff", plot = plot5, width = 6, height = 4, dpi = 600)
plot6 # Plot Weighted
ggsave("Weighted_Cnidaria.tiff", plot = plot6, width = 6, height = 4, dpi = 600)
plot7 # Plot Chao2
ggsave("Chao2_Cnidaria.tiff", plot = plot7, width = 6, height = 4, dpi = 600)

save.image(file = "Cnidaria_Analyses.RData")
load("Cnidaria_Analyses.RData")
