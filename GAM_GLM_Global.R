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
library(mgcv)
library(ggeffects)
library(DHARMa) #simulations package gam- 
library(knitr)
library(qpcR)
library(scales)


# ===================================
#                GLOBAL-GLM
# ===================================
#Species Counts and Environment, 5 degree bands
Ecological_Data_Global <- read.csv("Ecological_Data_Global.csv", sep = ";")
summary(Ecological_Data_Global)

#First we’re going to load in our data and then trim the data frame down to just the columns we need.
analysis.cols <- c("Margin_Sum", "Shelf_Sum", "CurVel_Mean", "Depth_Mean", "HumImp_Mean", "IceCov_Mean", 
                   "PhotoActi_Mean", "PrimProd_Mean", "Area_Sum", "O2_Mean", "Nitrate_Mean","ThemM_Mean", "ThemR_Mean", 
                   "NumSpe","NumRec", "ES50")
Ecological_Data_Global <- Ecological_Data_Global[,analysis.cols]
Ecological_Data_Global  <- Ecological_Data_Global [complete.cases(Ecological_Data_Global ),]

# Calculate the correlation matrix
corr_matrix <- cor(Ecological_Data_Global)

# Create the correlation plot with black font for text
corrplot(corr_matrix, tl.col = "black")

# Function to calculate significance (p-values) matrix
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Calculate the p-value matrix
p_matrix <- cor.mtest(Ecological_Data_Global)

# Set correlation values less than 0.6 to NA
corr_matrix[abs(corr_matrix) < 0.6] <- NA

# Set the tiff device to save the plot
tiff("Correlation_Plot_Global.tiff", width = 6, height = 4, units = "in", res = 600)

# Create the correlation plot with significant values and formatted coefficients
corrplot(corr_matrix, 
         tl.col = "black",        # Set font color to black
         addCoef.col = "white",   # Add correlation coefficients in white
         number.cex = 0.5,        # Set size of the correlation coefficients
         p.mat = p_matrix,        # Use p-value matrix to highlight significance
         sig.level = 0.05,        # Only show significant correlations (p < 0.05)
         insig = "blank",         # Leave insignificant correlations blank
         na.label = " ",          # Leave NA cells blank (no text or label)
         tl.cex = 0.7,            # Set smaller font size for row/column labels
         number.digits = 1)       # Display correlation coefficients with 1 digit after decimal

# Close the tiff device
dev.off()

#GLMs for number of species, Global data
global.numsp.intercept <- glm(NumSpe ~ 1, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.intercept)

global.numsp.numrec <- glm(NumSpe ~ NumRec, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.numrec)

global.numsp.margin <- glm(NumSpe ~ NumRec + Margin_Sum, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.margin)

global.numsp.shelf <- glm(NumSpe ~ NumRec + Shelf_Sum, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.shelf)

global.numsp.current <- glm(NumSpe ~ NumRec + CurVel_Mean, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.current)

global.numsp.depth <- glm(NumSpe ~ NumRec + Depth_Mean, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.depth)

global.numsp.humimp <- glm(NumSpe ~ NumRec + HumImp_Mean, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.humimp)

global.numsp.icecov <- glm(NumSpe ~ NumRec + IceCov_Mean, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.icecov)

global.numsp.phoact <- glm(NumSpe ~ NumRec + PhotoActi_Mean, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.phoact)

global.numsp.primprod <- glm(NumSpe ~ NumRec + PrimProd_Mean, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.primprod)

global.numsp.area <- glm(NumSpe ~ NumRec + Area_Sum, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.area)

global.numsp.o2 <- glm(NumSpe ~ NumRec + O2_Mean, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.o2)

global.numsp.nitrate <- glm(NumSpe ~ NumRec + Nitrate_Mean, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.nitrate)

global.numsp.themmean <- glm(NumSpe ~ NumRec + ThemM_Mean, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.themmean)

global.numsp.themrange <- glm(NumSpe ~ NumRec + ThemR_Mean, family = "poisson", data = Ecological_Data_Global)
summary(global.numsp.themrange)


#Model selection for number of species, Global
global.numsp.models <- list(global.numsp.intercept = global.numsp.intercept,
                            global.numsp.numrec = global.numsp.numrec,
                            global.numsp.margin = global.numsp.margin,
                            global.numsp.shelf = global.numsp.shelf,
                            global.numsp.current = global.numsp.current,
                            global.numsp.depth = global.numsp.depth,
                            global.numsp.humimp = global.numsp.humimp,
                            global.numsp.icecov = global.numsp.icecov,
                            global.numsp.phoact = global.numsp.phoact,
                            global.numsp.primprod = global.numsp.primprod,
                            global.numsp.area = global.numsp.area,
                            global.numsp.o2 = global.numsp.o2,
                            global.numsp.nitrate = global.numsp.nitrate,
                            global.numsp.themmean = global.numsp.themmean,
                            global.numsp.themrange = global.numsp.themrange)
global.numsp.aic.df <- data.frame(Model = names(global.numsp.models),
                                   AIC = sapply(global.numsp.models, function(x) AICc(x)),
                                   akaike.weights(sapply(global.numsp.models, function(x) AICc(x))))


global.numsp.aic.df <- global.numsp.aic.df[order(global.numsp.aic.df$AIC),]
global.numsp.aic.df$Cumulative.Weight <- cumsum(global.numsp.aic.df$weights)

kable(global.numsp.aic.df, row.names = FALSE)

#write.csv(global.numsp.aic.df, file = "global.5.degree.numsp.aic.csv")

Predicted_Res_Global <- ggplot(global.numsp.aic.df, aes(x = Model, y = AIC)) +
  geom_bar(stat = "identity", fill = "#1a80bb") +  # Create bars with blue color
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Variables",  
    y = "AIC"  
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),  # Increase x-axis text size and rotate labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )

ggsave("Predicted_Res_Global.tiff", plot = Predicted_Res_Global, width = 6, height = 4, dpi = 600)

#Plots for number of species, Global
Predicted_NumSp_NumRe_col <- ggplot(Ecological_Data_Global, aes(x = NumRec, y = predict(global.numsp.numrec, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Number of Records",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_NumRe_col.tiff", plot = Predicted_NumSp_NumRe_col, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Margin_col <- ggplot(Ecological_Data_Global, aes(x = Margin_Sum, y = predict(global.numsp.margin, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Continental Margin (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Margin_col.tiff", plot = Predicted_NumSp_Margin_col, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Shelf_col <- ggplot(Ecological_Data_Global, aes(x = Shelf_Sum, y = predict(global.numsp.shelf, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Continental Shelf (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Shelf_col.tiff", plot = Predicted_NumSp_Shelf_col, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Current_col <- ggplot(Ecological_Data_Global, aes(x = CurVel_Mean, y = predict(global.numsp.current, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Current Velocity (m.s-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Current_col.tiff", plot = Predicted_NumSp_Current_col, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Depth_col <- ggplot(Ecological_Data_Global, aes(x = Depth_Mean, y = predict(global.numsp.depth, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Depth (m)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Depth_col.tiff", plot = Predicted_NumSp_Depth_col, width = 6, height = 4, dpi = 600)

Predicted_NumSp_humimp_col <- ggplot(Ecological_Data_Global, aes(x = HumImp_Mean, y = predict(global.numsp.humimp, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Human Impact",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_humipm_col.tiff", plot = Predicted_NumSp_humimp_col, width = 6, height = 4, dpi = 600)

Predicted_NumSp_icecover_col <- ggplot(Ecological_Data_Global, aes(x = IceCov_Mean, y = predict(global.numsp.icecov, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Ice Cover (fraction)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
    theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_icecover_col.tiff", plot = Predicted_NumSp_icecover_col, width = 6, height = 4, dpi = 600)

Predicted_NumSp_phoact_col <- ggplot(Ecological_Data_Global, aes(x = PhotoActi_Mean, y = predict(global.numsp.phoact, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Phot. Avai. Radi (E.m-2.day-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_phoact_col.tiff", plot = Predicted_NumSp_phoact_col, width = 6, height = 4, dpi = 600)

Predicted_NumSp_primprod_col <- ggplot(Ecological_Data_Global, aes(x = PrimProd_Mean, y = predict(global.numsp.primprod, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Primary Productivity (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_primprod_col.tiff", plot = Predicted_NumSp_primprod_col, width = 6, height = 4, dpi = 600)

Predicted_NumSp_area_col <- ggplot(Ecological_Data_Global, aes(x = Area_Sum, y = predict(global.numsp.area, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  #geom_smooth(method = "glm", formula = y ~ x, color = "black") +  # Add a smooth black line using GLM
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Area (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_area_col.tiff", plot = Predicted_NumSp_area_col, width = 6, height = 4, dpi = 600)

Predicted_NumSp_O2_col <- ggplot(Ecological_Data_Global, aes(x = O2_Mean, y = predict(global.numsp.o2, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "O2 (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_O2_col.tiff", plot = Predicted_NumSp_O2_col, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Nitrate_col <- ggplot(Ecological_Data_Global, aes(x = Nitrate_Mean, y = predict(global.numsp.nitrate, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Nitrate (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Nitrate_col.tiff", plot = Predicted_NumSp_Nitrate_col, width = 6, height = 4, dpi = 600)

Predicted_NumSp_themmean_col <- ggplot(Ecological_Data_Global, aes(x = ThemM_Mean, y = predict(global.numsp.themmean, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  #geom_smooth(method = "glm", formula = y ~ x, color = "black") +  # Add a smooth black line using GLM
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Mean (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_themmean_col.tiff", plot = Predicted_NumSp_themmean_col, width = 6, height = 4, dpi = 600)

Predicted_NumSp_themrange_col <- ggplot(Ecological_Data_Global, aes(x = ThemR_Mean, y = predict(global.numsp.themrange, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  #geom_smooth(method = "glm", formula = y ~ x, color = "black") +  # Add a smooth black line using GLM
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Range (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_themrange_col.tiff", plot = Predicted_NumSp_themrange_col, width = 6, height = 4, dpi = 600)


# ===================================
#                SHALLOW-NUMSPe-GLM
# ===================================

#Species Counts and Environment, 5 degree bands
Ecological_Data_Global <- read.csv("Ecological_Data_Global.csv", sep = ";")
summary(Ecological_Data_Global)

#First we’re going to load in our data and then trim the data frame down to just the columns we need.
analysis.cols <- c("Shelf_Sum", "CurVel_Mean", "Depth_Mean", "HumImp_Mean", "IceCov_Mean", "PrimProd_Mean", "Area_Sum", "Nitrate_Mean", 
                   "ThemM_Mean", "NumSpe_Sha","NumRec_Sha", "NumPhy_Sha", "ES50_Sha")
Ecological_Data_Global <- Ecological_Data_Global[,analysis.cols]
Ecological_Data_Global  <- Ecological_Data_Global [complete.cases(Ecological_Data_Global ),]

# Calculate the correlation matrix
corr_matrix <- cor(Ecological_Data_Global)

# Create the correlation plot with black font for text
corrplot(corr_matrix, tl.col = "black")

# Function to calculate significance (p-values) matrix
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Calculate the p-value matrix
p_matrix <- cor.mtest(Ecological_Data_Global)

# Set correlation values less than 0.6 to NA
corr_matrix[abs(corr_matrix) < 0.6] <- NA

# Set the tiff device to save the plot
tiff("Correlation_Plot_Shallow.tiff", width = 6, height = 4, units = "in", res = 600)

# Create the correlation plot with significant values and formatted coefficients
corrplot(corr_matrix, 
         tl.col = "black",        # Set font color to black
         addCoef.col = "white",   # Add correlation coefficients in white
         number.cex = 0.5,        # Set size of the correlation coefficients
         p.mat = p_matrix,        # Use p-value matrix to highlight significance
         sig.level = 0.05,        # Only show significant correlations (p < 0.05)
         insig = "blank",         # Leave insignificant correlations blank
         na.label = " ",          # Leave NA cells blank (no text or label)
         tl.cex = 0.7,            # Set smaller font size for row/column labels
         number.digits = 1)       # Display correlation coefficients with 1 digit after decimal

# Close the tiff device
dev.off()

#GLMs for number of species, Global data, Shallow
shallow.numsp.intercept <- glm(NumSpe_Sha ~ 1, family = "poisson", data = Ecological_Data_Global)
summary(shallow.numsp.intercept)

shallow.numsp.numrec <- glm(NumSpe_Sha ~ NumRec_Sha, family = "poisson", data = Ecological_Data_Global)
summary(shallow.numsp.numrec)

shallow.numsp.shelfsum <- glm(NumSpe_Sha ~ NumRec_Sha + Shelf_Sum, family = "poisson", data = Ecological_Data_Global)
summary(shallow.numsp.shelfsum)

shallow.numsp.current <- glm(NumSpe_Sha ~ NumRec_Sha + CurVel_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.numsp.current)

shallow.numsp.depth <- glm(NumSpe_Sha ~ NumRec_Sha + Depth_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.numsp.depth)

shallow.numsp.humimp <- glm(NumSpe_Sha ~ NumRec_Sha + HumImp_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.numsp.humimp)

shallow.numsp.icecov <- glm(NumSpe_Sha ~ NumRec_Sha + IceCov_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.numsp.icecov)

shallow.numsp.primprod <- glm(NumSpe_Sha ~ NumRec_Sha + PrimProd_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.numsp.primprod)

shallow.numsp.themmean <- glm(NumSpe_Sha ~ NumRec_Sha + ThemM_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.numsp.themmean)

shallow.numsp.area <- glm(NumSpe_Sha ~ NumRec_Sha + Area_Sum, family = "poisson", data = Ecological_Data_Global)
summary(shallow.numsp.area)

shallow.numsp.nitrate <- glm(NumSpe_Sha ~ NumRec_Sha + Nitrate_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.numsp.nitrate)

#Model selection for number of species, Global
shallow.numsp.models <- list(Intercept = shallow.numsp.intercept,
                            NumRec = shallow.numsp.numrec,
                            ConShe = shallow.numsp.shelfsum,
                            CurVel = shallow.numsp.current,
                            Depth = shallow.numsp.depth,
                            HumImp = shallow.numsp.humimp,
                            IceCov = shallow.numsp.icecov,
                            PriPro = shallow.numsp.primprod,
                            TheMea = shallow.numsp.themmean,
                            Area = shallow.numsp.area,
                            Nitrate = shallow.numsp.nitrate)
shallow.numsp.aic.df <- data.frame(Model = names(shallow.numsp.models),
                                  AIC = sapply(shallow.numsp.models, function(x) AICc(x)),
                                  akaike.weights(sapply(shallow.numsp.models, function(x) AICc(x))))


shallow.numsp.aic.df <- shallow.numsp.aic.df[order(shallow.numsp.aic.df$AIC),]
shallow.numsp.aic.df$Cumulative.Weight <- cumsum(shallow.numsp.aic.df$weights)

kable(shallow.numsp.aic.df, row.names = FALSE)

#write.csv(shallow.numsp.aic.df, file = "global.5.degree.numsp.shallow.aic.csv")

#Plots for number of species, Shallow
Predicted_NumSp_NumRe_sha <- ggplot(Ecological_Data_Global, aes(x = NumRec_Sha, y = predict(shallow.numsp.numrec, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  labs(
    x = "Number of Records",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_NumRe_sha.tiff", plot = Predicted_NumSp_NumRe_sha, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Shelf_sha <- ggplot(Ecological_Data_Global, aes(x = Shelf_Sum, y = predict(shallow.numsp.shelfsum, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  labs(
    x = "Continental Shelf (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Shelf_sha.tiff", plot = Predicted_NumSp_Shelf_sha, width = 6, height = 4, dpi = 600)

Predicted_NumSp_CurVel_sha <- ggplot(Ecological_Data_Global, aes(x = CurVel_Mean, y = predict(shallow.numsp.current, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Current Velocity (m.s-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_CurVel_sha.tiff", plot = Predicted_NumSp_CurVel_sha, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Depth_sha <- ggplot(Ecological_Data_Global, aes(x = Depth_Mean, y = predict(shallow.numsp.depth, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  labs(
    x = "Depth (m)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Depth_sha.tiff", plot = Predicted_NumSp_Depth_sha, width = 6, height = 4, dpi = 600)

Predicted_NumSp_humimp_sha <- ggplot(Ecological_Data_Global, aes(x = HumImp_Mean, y = predict(shallow.numsp.humimp, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Human Impact",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_humipm_sha.tiff", plot = Predicted_NumSp_humimp_sha, width = 6, height = 4, dpi = 600)

Predicted_NumSp_icecover_sha <- ggplot(Ecological_Data_Global, aes(x = IceCov_Mean, y = predict(shallow.numsp.icecov, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Ice Cover (fraction)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_icecover_sha.tiff", plot = Predicted_NumSp_icecover_sha, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_phoact_sha <- ggplot(Ecological_Data_Global, aes(x = PhotoActi_Mean, y = predict(shallow.numsp.phoact, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Phot. Avai. Radi (E.m-2.day-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_phoact_sha.tiff", plot = Predicted_NumSp_phoact_sha, width = 6, height = 4, dpi = 600)

Predicted_NumSp_primprod_sha <- ggplot(Ecological_Data_Global, aes(x = PrimProd_Mean, y = predict(shallow.numsp.primprod, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Primary Productivity (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_primprod_sha.tiff", plot = Predicted_NumSp_primprod_sha, width = 6, height = 4, dpi = 600)

Predicted_NumSp_themmean_sha <- ggplot(Ecological_Data_Global, aes(x = ThemM_Mean, y = predict(shallow.numsp.themmean, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Mean (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_themmean_sha.tiff", plot = Predicted_NumSp_themmean_sha, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_themrange_sha <- ggplot(Ecological_Data_Global, aes(x = ThemR_Mean, y = predict(shallow.numsp.themrange, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Range (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_themrange_sha.tiff", plot = Predicted_NumSp_themrange_sha, width = 6, height = 4, dpi = 600)

Predicted_NumSp_area_sha <- ggplot(Ecological_Data_Global, aes(x = Area_Sum, y = predict(shallow.numsp.area, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Area (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_area_sha.tiff", plot = Predicted_NumSp_area_sha, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_oxygen_sha <- ggplot(Ecological_Data_Global, aes(x = O2_Mean, y = predict(shallow.numsp.oxygen, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "O2 (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_oxygen_sha.tiff", plot = Predicted_NumSp_oxygen_sha, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Nitrate_sha <- ggplot(Ecological_Data_Global, aes(x = Nitrate_Mean, y = predict(shallow.numsp.nitrate, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Nitrate (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Nitrate_sha.tiff", plot = Predicted_NumSp_Nitrate_sha, width = 6, height = 4, dpi = 600)

# ===================================
#                SHALLOW-ES50-GLM
# ===================================
shallow.es50.intercept <- glm(ES50_Sha ~ 1, family = "poisson", data = Ecological_Data_Global)
summary(shallow.es50.intercept)

shallow.es50.themmean <- glm(ES50_Sha ~ ThemM_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.es50.themmean)

shallow.es50.numrec <- glm(ES50_Sha ~ NumRec_Sha, family = "poisson", data = Ecological_Data_Global)
summary(shallow.es50.numrec)

shallow.es50.shelfsum <- glm(ES50_Sha ~ NumRec_Sha + Shelf_Sum, family = "poisson", data = Ecological_Data_Global)
summary(shallow.es50.shelfsum)

shallow.es50.current <- glm(ES50_Sha ~ NumRec_Sha + CurVel_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.es50.current)

shallow.es50.depth <- glm(ES50_Sha ~ Depth_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.es50.depth)

shallow.es50.humimp <- glm(ES50_Sha ~ HumImp_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.es50.humimp)

shallow.es50.icecov <- glm(ES50_Sha ~ IceCov_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.es50.icecov)

#shallow.es50.phoact <- glm(ES50_Sha ~ PhotoActi_Mean, family = "poisson", data = Ecological_Data_Global)
#summary(shallow.es50.phoact)

shallow.es50.primprod <- glm(ES50_Sha ~ PrimProd_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.es50.primprod)

#shallow.es50.themrange <- glm(ES50_Sha ~ ThemR_Mean, family = "poisson", data = Ecological_Data_Global)
#summary(shallow.es50.themrange)

shallow.es50.area <- glm(ES50_Sha ~ NumRec_Sha + Area_Sum, family = "poisson", data = Ecological_Data_Global)
summary(shallow.es50.area)

#shallow.es50.oxygen <- glm(ES50_Sha ~ O2_Mean, family = "poisson", data = Ecological_Data_Global)
#summary(shallow.es50.oxygen)

shallow.es50.nitrate <- glm(ES50_Sha ~ Nitrate_Mean, family = "poisson", data = Ecological_Data_Global)
summary(shallow.es50.nitrate)

#Model selection for number of species, Global
shallow.es50.models <- list(Intercept = shallow.es50.intercept,
                           NumRec = shallow.es50.numrec,
                           ConShe = shallow.es50.shelfsum,
                           CurVel = shallow.es50.current,
                           Depth = shallow.es50.depth,
                           HumImp = shallow.es50.humimp,
                           IceCov = shallow.es50.icecov,
                           #PhoAct = shallow.es50.phoact,
                           PriPro = shallow.es50.primprod,
                           TheMea = shallow.es50.themmean,
                           #TheRan = shallow.es50.themrange,
                           Area = shallow.es50.area,
                           #O2 = shallow.es50.oxygen,
                           Nitrate = shallow.es50.nitrate)
shallow.es50.aic.df <- data.frame(Model = names(shallow.es50.models),
                                 AIC = sapply(shallow.es50.models, function(x) AICc(x)),
                                 akaike.weights(sapply(shallow.es50.models, function(x) AICc(x))))


shallow.es50.aic.df <- shallow.es50.aic.df[order(shallow.es50.aic.df$AIC),]
shallow.es50.aic.df$Cumulative.Weight <- cumsum(shallow.es50.aic.df$weights)

kable(shallow.es50.aic.df, row.names = FALSE)

#write.csv(shallow.es50.aic.df, file = "global.5.degree.es50.shallow.aic.csv")

#Plots for ES50, shallow
Predicted_ES50_NumRe_sha <- ggplot(Ecological_Data_Global, aes(x = NumRec_Sha, y = ES50_Sha)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(limits=c(18, 52), breaks=seq(20, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Number of Records",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_NumRe_sha.tiff", plot = Predicted_ES50_NumRe_sha, width = 6, height = 4, dpi = 600)

Predicted_ES50_Shelf_sha <- ggplot(Ecological_Data_Global, aes(x = Shelf_Sum, y = ES50_Sha)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(limits=c(18, 52), breaks=seq(20, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Continental Shelf (km2)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_Shelf_sha.tiff", plot = Predicted_ES50_Shelf_sha, width = 6, height = 4, dpi = 600)

Predicted_ES50_CurVel_sha <- ggplot(Ecological_Data_Global, aes(x = CurVel_Mean, y = ES50_Sha)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(18, 52), breaks=seq(20, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Current Velocity (m.s-1)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_CurVel_sha.tiff", plot = Predicted_ES50_CurVel_sha, width = 6, height = 4, dpi = 600)

Predicted_ES50_Depth_sha <- ggplot(Ecological_Data_Global, aes(x = Depth_Mean, y = ES50_Sha)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(limits=c(18, 52), breaks=seq(20, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Depth (m)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_Depth_sha.tiff", plot = Predicted_ES50_Depth_sha, width = 6, height = 4, dpi = 600)

Predicted_ES50_humimp_sha <- ggplot(Ecological_Data_Global, aes(x = HumImp_Mean, y = ES50_Sha)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(18, 52), breaks=seq(20, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Human Impact",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_humipm_sha.tiff", plot = Predicted_ES50_humimp_sha, width = 6, height = 4, dpi = 600)

Predicted_ES50_icecover_sha <- ggplot(Ecological_Data_Global, aes(x = IceCov_Mean, y = ES50_Sha)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(18, 52), breaks=seq(20, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Ice Cover (fraction)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_icecover_sha.tiff", plot = Predicted_ES50_icecover_sha, width = 6, height = 4, dpi = 600)

#Predicted_ES50_phoact_sha <- ggplot(Ecological_Data_Global, aes(x = PhotoActi_Mean, y = ES50_Sha)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(18, 52), breaks=seq(20, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Phot. Avai. Radi (E.m-2.day-1)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_phoact_sha.tiff", plot = Predicted_ES50_phoact_sha, width = 6, height = 4, dpi = 600)

Predicted_ES50_primprod_sha <- ggplot(Ecological_Data_Global, aes(x = PrimProd_Mean, y = ES50_Sha)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(18, 52), breaks=seq(20, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Primary Productivity (mmol . m-3)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_primprod_sha.tiff", plot = Predicted_ES50_primprod_sha, width = 6, height = 4, dpi = 600)

Predicted_ES50_themmean_sha <- ggplot(Ecological_Data_Global, aes(x = ThemM_Mean, y = ES50_Sha)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(18, 52), breaks=seq(20, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Temperature Mean (ºC)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_themmean_sha.tiff", plot = Predicted_ES50_themmean_sha, width = 6, height = 4, dpi = 600)

#Predicted_ES50_themrange_sha <- ggplot(Ecological_Data_Global, aes(x = ThemR_Mean, y = ES50_Sha)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(18, 52), breaks=seq(20, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Temperature Range (ºC)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_themrange_sha.tiff", plot = Predicted_ES50_themrange_sha, width = 6, height = 4, dpi = 600)

Predicted_ES50_area_sha <- ggplot(Ecological_Data_Global, aes(x = Area_Sum, y = ES50_Sha)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(limits=c(18, 52), breaks=seq(20, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Area (km2)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_area_sha.tiff", plot = Predicted_ES50_area_sha, width = 6, height = 4, dpi = 600)

#Predicted_ES50_oxygen_sha <- ggplot(Ecological_Data_Global, aes(x = O2_Mean, y = ES50_Sha)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(18, 52), breaks=seq(20, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "O2 (mmol . m-3)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_oxygen_sha.tiff", plot = Predicted_ES50_oxygen_sha, width = 6, height = 4, dpi = 600)

Predicted_ES50_Nitrate_sha <- ggplot(Ecological_Data_Global, aes(x = Nitrate_Mean, y = ES50_Sha)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(18, 52), breaks=seq(20, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Nitrate (mmol . m-3)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_Nitrate_sha.tiff", plot = Predicted_ES50_Nitrate_sha, width = 6, height = 4, dpi = 600)


# ===================================
#                MESO-NUMSPe-GLM
# ===================================

#Species Counts and Environment, 5 degree bands
Ecological_Data_Global <- read.csv("Ecological_Data_Global.csv", sep = ";")
summary(Ecological_Data_Global)

#First we’re going to load in our data and then trim the data frame down to just the columns we need.
analysis.cols <- c("Margin_Sum","CurVel_Mean", "Depth_Mean", "HumImp_Mean", 
                   "PrimProd_Mean", "ThemM_Mean", "Area_Sum", "Nitrate_Mean", "NumSpe_Mes",
                   "NumRec_Mes", "NumPhy_Mes", "ES50_Mes")
Ecological_Data_Global <- Ecological_Data_Global[,analysis.cols]
Ecological_Data_Global  <- Ecological_Data_Global [complete.cases(Ecological_Data_Global ),]

# Calculate the correlation matrix
corr_matrix <- cor(Ecological_Data_Global)

# Create the correlation plot with black font for text
corrplot(corr_matrix, tl.col = "black")

# Function to calculate significance (p-values) matrix
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Calculate the p-value matrix
p_matrix <- cor.mtest(Ecological_Data_Global)

# Set correlation values less than 0.6 to NA
corr_matrix[abs(corr_matrix) < 0.6] <- NA

# Set the tiff device to save the plot
tiff("Correlation_Plot_Meso.tiff", width = 6, height = 4, units = "in", res = 600)

# Create the correlation plot with significant values and formatted coefficients
corrplot(corr_matrix, 
         tl.col = "black",        # Set font color to black
         addCoef.col = "white",   # Add correlation coefficients in white
         number.cex = 0.5,        # Set size of the correlation coefficients
         p.mat = p_matrix,        # Use p-value matrix to highlight significance
         sig.level = 0.05,        # Only show significant correlations (p < 0.05)
         insig = "blank",         # Leave insignificant correlations blank
         na.label = " ",          # Leave NA cells blank (no text or label)
         tl.cex = 0.7,            # Set smaller font size for row/column labels
         number.digits = 1)       # Display correlation coefficients with 1 digit after decimal

# Close the tiff device
dev.off()

#GLMs for number of species, Global data, Meso
meso.numsp.intercept <- glm(NumSpe_Mes ~ 1, family = "poisson", data = Ecological_Data_Global)
summary(meso.numsp.intercept)

meso.numsp.numrec <- glm(NumSpe_Mes ~ NumRec_Mes, family = "poisson", data = Ecological_Data_Global)
summary(meso.numsp.numrec)

meso.numsp.marginsum <- glm(NumSpe_Mes ~ NumRec_Mes + Margin_Sum, family = "poisson", data = Ecological_Data_Global)
summary(meso.numsp.marginsum)

meso.numsp.current <- glm(NumSpe_Mes ~ NumRec_Mes + CurVel_Mean, family = "poisson", data = Ecological_Data_Global)
summary(meso.numsp.current)

meso.numsp.depth <- glm(NumSpe_Mes ~ NumRec_Mes + Depth_Mean, family = "poisson", data = Ecological_Data_Global)
summary(meso.numsp.depth)

meso.numsp.humimp <- glm(NumSpe_Mes ~ NumRec_Mes + HumImp_Mean, family = "poisson", data = Ecological_Data_Global)
summary(meso.numsp.humimp)

#meso.numsp.icecov <- glm(NumSpe_Mes ~ NumRec_Mes + IceCov_Mean, family = "poisson", data = Ecological_Data_Global)
#summary(meso.numsp.icecov)

#meso.numsp.phoact <- glm(NumSpe_Mes ~ NumRec_Mes + PhotoActi_Mean, family = "poisson", data = Ecological_Data_Global)
#summary(meso.numsp.phoact)

meso.numsp.primprod <- glm(NumSpe_Mes ~ NumRec_Mes + PrimProd_Mean, family = "poisson", data = Ecological_Data_Global)
summary(meso.numsp.primprod)

meso.numsp.themmean <- glm(NumSpe_Mes ~ NumRec_Mes + ThemM_Mean, family = "poisson", data = Ecological_Data_Global)
summary(meso.numsp.themmean)

#meso.numsp.themrange <- glm(NumSpe_Mes ~ NumRec_Mes + ThemR_Mean, family = "poisson", data = Ecological_Data_Global)
#summary(meso.numsp.themrange)

meso.numsp.area <- glm(NumSpe_Mes ~ NumRec_Mes + Area_Sum, family = "poisson", data = Ecological_Data_Global)
summary(meso.numsp.area)

#meso.numsp.oxygen <- glm(NumSpe_Mes ~ NumRec_Mes + O2_Mean, family = "poisson", data = Ecological_Data_Global)
#summary(meso.numsp.oxygen)

meso.numsp.nitrate <- glm(NumSpe_Mes ~ NumRec_Mes + Nitrate_Mean, family = "poisson", data = Ecological_Data_Global)
summary(meso.numsp.nitrate)


#Model selection for number of species, Global
meso.numsp.models <- list(Intercept = meso.numsp.intercept,
                            NumRec = meso.numsp.numrec,
                            ConMar = meso.numsp.marginsum,
                            CurVel = meso.numsp.current,
                            Depth = meso.numsp.depth,
                            HumImp = meso.numsp.humimp,
                            #IceCov = meso.numsp.icecov,
                            #PhoAct = meso.numsp.phoact,
                            PriPro = meso.numsp.primprod,
                            TheMea = meso.numsp.themmean,
                            #TheRan = meso.numsp.themrange,
                            Area = meso.numsp.area,
                            #O2 = meso.numsp.oxygen,
                            Nitrate = meso.numsp.nitrate)
meso.numsp.aic.df <- data.frame(Model = names(meso.numsp.models),
                                  AIC = sapply(meso.numsp.models, function(x) AICc(x)),
                                  akaike.weights(sapply(meso.numsp.models, function(x) AICc(x))))


meso.numsp.aic.df <- meso.numsp.aic.df[order(meso.numsp.aic.df$AIC),]
meso.numsp.aic.df$Cumulative.Weight <- cumsum(meso.numsp.aic.df$weights)

kable(meso.numsp.aic.df, row.names = FALSE)

#write.csv(meso.numsp.aic.df, file = "global.5.degree.numsp.meso.aic.csv")

#Plots for number of species, Global
Predicted_NumSp_NumRe_mes <- ggplot(Ecological_Data_Global, aes(x = NumRec_Mes, y = predict(meso.numsp.numrec, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  labs(
    x = "Number of Records",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_NumRe_mes.tiff", plot = Predicted_NumSp_NumRe_mes, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Margin_mes <- ggplot(Ecological_Data_Global, aes(x = Margin_Sum, y = predict(meso.numsp.marginsum, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  labs(
    x = "Continental Margin (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Margin_mes.tiff", plot = Predicted_NumSp_Margin_mes, width = 6, height = 4, dpi = 600)

Predicted_NumSp_CurVel_mes <- ggplot(Ecological_Data_Global, aes(x = CurVel_Mean, y = predict(meso.numsp.marginsum, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Current Velocity (m.s-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_CurVel_mes.tiff", plot = Predicted_NumSp_CurVel_mes, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Depth_mes <- ggplot(Ecological_Data_Global, aes(x = Depth_Mean, y = predict(meso.numsp.depth, Ecological_Data_Global))) +
  geom_point(size = 3) +  # Add scatter plot points
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  labs(
    x = "Depth (m)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Depth_mes.tiff", plot = Predicted_NumSp_Depth_mes, width = 6, height = 4, dpi = 600)

Predicted_NumSp_humimp_mes <- ggplot(Ecological_Data_Global, aes(x = HumImp_Mean, y = predict(meso.numsp.humimp, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Human Impact",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_humipm_mes.tiff", plot = Predicted_NumSp_humimp_mes, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_icecover_mes <- ggplot(Ecological_Data_Global, aes(x = IceCov_Mean, y = predict(meso.numsp.icecov, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Ice Cover (fraction)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_icecover_mes.tiff", plot = Predicted_NumSp_icecover_mes, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_phoact_mes <- ggplot(Ecological_Data_Global, aes(x = PhotoActi_Mean, y = predict(meso.numsp.phoact, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Phot. Avai. Radi (E.m-2.day-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_phoact_mes.tiff", plot = Predicted_NumSp_phoact_mes, width = 6, height = 4, dpi = 600)


Predicted_NumSp_primprod_mes <- ggplot(Ecological_Data_Global, aes(x = PrimProd_Mean, y = predict(meso.numsp.primprod, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Primary Productivity (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_primprod_mes.tiff", plot = Predicted_NumSp_primprod_mes, width = 6, height = 4, dpi = 600)

Predicted_NumSp_themmean_mes <- ggplot(Ecological_Data_Global, aes(x = ThemM_Mean, y = predict(meso.numsp.themmean, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Mean (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_themmean_mes.tiff", plot = Predicted_NumSp_themmean_mes, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_themrange_mes <- ggplot(Ecological_Data_Global, aes(x = ThemR_Mean, y = predict(meso.numsp.themrange, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Range (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_themrange_mes.tiff", plot = Predicted_NumSp_themrange_mes, width = 6, height = 4, dpi = 600)

Predicted_NumSp_area_mes <- ggplot(Ecological_Data_Global, aes(x = Area_Sum, y = predict(meso.numsp.area, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  labs(
    x = "Area (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_area_mes.tiff", plot = Predicted_NumSp_area_mes, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_oxygen_mes <- ggplot(Ecological_Data_Global, aes(x = O2_Mean, y = predict(meso.numsp.oxygen, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "O2 (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_oxygen_mes.tiff", plot = Predicted_NumSp_oxygen_mes, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Nitrate_mes <- ggplot(Ecological_Data_Global, aes(x = Nitrate_Mean, y = predict(meso.numsp.nitrate, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Nitrate (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Nitrate_mes.tiff", plot = Predicted_NumSp_Nitrate_mes, width = 6, height = 4, dpi = 600)

# ===================================
#                MESO-ES50-GLM
# ===================================

meso.es50.intercept <- glm(ES50_Mes ~ 1, family = "poisson", data = Ecological_Data_Global)
summary(meso.es50.intercept)

meso.es50.numrec <- glm(ES50_Mes ~ NumRec_Mes, family = "poisson", data = Ecological_Data_Global)
summary(meso.es50.numrec)

meso.es50.marginsum <- glm(ES50_Mes ~ Margin_Sum, family = "poisson", data = Ecological_Data_Global)
summary(meso.es50.marginsum)

meso.es50.current <- glm(ES50_Mes ~ CurVel_Mean, family = "poisson", data = Ecological_Data_Global)
summary(meso.es50.current)

meso.es50.depth <- glm(ES50_Mes ~ Depth_Mean, family = "poisson", data = Ecological_Data_Global)
summary(meso.es50.depth)

meso.es50.humimp <- glm(ES50_Mes ~ HumImp_Mean, family = "poisson", data = Ecological_Data_Global)
summary(meso.es50.humimp)

#meso.es50.icecov <- glm(ES50_Mes ~ IceCov_Mean, family = "poisson", data = Ecological_Data_Global)
#summary(meso.es50.icecov)

#meso.es50.phoact <- glm(ES50_Mes ~ PhotoActi_Mean, family = "poisson", data = Ecological_Data_Global)
#summary(meso.es50.phoact)

meso.es50.primprod <- glm(ES50_Mes ~ PrimProd_Mean, family = "poisson", data = Ecological_Data_Global)
summary(meso.es50.primprod)

meso.es50.themmean <- glm(ES50_Mes ~ ThemM_Mean, family = "poisson", data = Ecological_Data_Global)
summary(meso.es50.themmean)

#meso.es50.themrange <- glm(ES50_Mes ~ ThemR_Mean, family = "poisson", data = Ecological_Data_Global)
#summary(meso.es50.themrange)

meso.es50.area <- glm(ES50_Mes ~ Area_Sum, family = "poisson", data = Ecological_Data_Global)
summary(meso.es50.area)

#meso.es50.oxygen <- glm(ES50_Mes ~ O2_Mean, family = "poisson", data = Ecological_Data_Global)
#summary(meso.es50.oxygen)

meso.es50.nitrate <- glm(ES50_Mes ~ Nitrate_Mean, family = "poisson", data = Ecological_Data_Global)
summary(meso.es50.nitrate)

#Model selection for number of species, Global
meso.es50.models <- list(Intercept = meso.es50.intercept,
                            NumRec = meso.es50.numrec,
                            ConMar = meso.es50.marginsum,
                            CurVel = meso.es50.current,
                            Depth = meso.es50.depth,
                            HumIpm = meso.es50.humimp,
                            #IceCov = meso.es50.icecov,
                            #PhoAct = meso.es50.phoact,
                            PriPro = meso.es50.primprod,
                            TheMea = meso.es50.themmean,
                            #TheRan = meso.es50.themrange,
                            Area = meso.es50.area,
                            #O2 = meso.es50.oxygen,
                            Nitrate = meso.es50.nitrate)
meso.es50.aic.df <- data.frame(Model = names(meso.es50.models),
                                  AIC = sapply(meso.es50.models, function(x) AICc(x)),
                                  akaike.weights(sapply(meso.es50.models, function(x) AICc(x))))


meso.es50.aic.df <-meso.es50.aic.df[order(meso.es50.aic.df$AIC),]
meso.es50.aic.df$Cumulative.Weight <- cumsum(meso.es50.aic.df$weights)

kable(meso.es50.aic.df, row.names = FALSE)

#write.csv(meso.es50.aic.df, file = "global.5.degree.es50.meso.aic.csv")

#Plots for number of species, Global
Predicted_ES50_NumRe_mes <- ggplot(Ecological_Data_Global, aes(x = NumRec_Mes, y = ES50_Mes)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Number of Records",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_NumRe_mes.tiff", plot = Predicted_ES50_NumRe_mes, width = 6, height = 4, dpi = 600)

Predicted_ES50_Margin_mes <- ggplot(Ecological_Data_Global, aes(x = Margin_Sum, y = ES50_Mes)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Continental Margin (km2)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_Margin_mes.tiff", plot = Predicted_ES50_Margin_mes, width = 6, height = 4, dpi = 600)

Predicted_ES50_CurVel_mes <- ggplot(Ecological_Data_Global, aes(x = CurVel_Mean, y = ES50_Mes)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Current Velocity (m.s-1)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_CurVel_mes.tiff", plot = Predicted_ES50_CurVel_mes, width = 6, height = 4, dpi = 600)

Predicted_ES50_Depth_mes <- ggplot(Ecological_Data_Global, aes(x = Depth_Mean, y = ES50_Mes)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Depth (m)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_Depth_mes.tiff", plot = Predicted_ES50_Depth_mes, width = 6, height = 4, dpi = 600)

Predicted_ES50_humimp_mes <- ggplot(Ecological_Data_Global, aes(x = HumImp_Mean, y = ES50_Mes)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Human Impact",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_humipm_mes.tiff", plot = Predicted_ES50_humimp_mes, width = 6, height = 4, dpi = 600)

#Predicted_ES50_icecover_mes <- ggplot(Ecological_Data_Global, aes(x = IceCov_Mean, y = ES50_Mes)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Ice Cover (fraction)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_icecover_mes.tiff", plot = Predicted_ES50_icecover_mes, width = 6, height = 4, dpi = 600)

#Predicted_ES50_phoact_mes <- ggplot(Ecological_Data_Global, aes(x = PhotoActi_Mean, y = ES50_Mes)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Phot. Avai. Radi (E.m-2.day-1)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_phoact_mes.tiff", plot = Predicted_ES50_phoact_mes, width = 6, height = 4, dpi = 600)

Predicted_ES50_primprod_mes <- ggplot(Ecological_Data_Global, aes(x = PrimProd_Mean, y = ES50_Mes)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Primary Productivity (mmol . m-3)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_primprod_mes.tiff", plot = Predicted_ES50_primprod_mes, width = 6, height = 4, dpi = 600)

Predicted_ES50_themmean_mes <- ggplot(Ecological_Data_Global, aes(x = ThemM_Mean, y = ES50_Mes)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Temperature Mean (ºC)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_themmean_mes.tiff", plot = Predicted_ES50_themmean_mes, width = 6, height = 4, dpi = 600)

#Predicted_ES50_themrange_mes <- ggplot(Ecological_Data_Global, aes(x = ThemR_Mean, y = ES50_Mes)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Temperature Range (ºC)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_themrange_mes.tiff", plot = Predicted_ES50_themrange_mes, width = 6, height = 4, dpi = 600)

Predicted_ES50_area_mes <- ggplot(Ecological_Data_Global, aes(x = Area_Sum, y = ES50_Mes)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Area (km2)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_area_mes.tiff", plot = Predicted_ES50_area_mes, width = 6, height = 4, dpi = 600)

#Predicted_ES50_oxygen_mes <- ggplot(Ecological_Data_Global, aes(x = O2_Mean, y = ES50_Mes)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "O2 (mmol . m-3)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_oxygen_mes.tiff", plot = Predicted_ES50_oxygen_mes, width = 6, height = 4, dpi = 600)

Predicted_ES50_Nitrate_mes <- ggplot(Ecological_Data_Global, aes(x = Nitrate_Mean, y = ES50_Mes)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Nitrate (mmol . m-3)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_Nitrate_mes.tiff", plot = Predicted_ES50_Nitrate_mes, width = 6, height = 4, dpi = 600)


# ===================================
#                DEEP-NUMSPe-GLM
# ===================================

#Species Counts and Environment, 5 degree bands
Ecological_Data_Global <- read.csv("Ecological_Data_Global.csv", sep = ";")
summary(Ecological_Data_Global)

#First we’re going to load in our data and then trim the data frame down to just the columns we need.
analysis.cols <- c("Margin_Sum", "Depth_Mean", "HumImp_Mean", "CurVel_Bot_Mean", "Nitrate_Bot_Mean",	
                   "O2_Bot_Mean","Area_Sum","PrimProd_Bot_Mean", "ThemM_Bot_Mean", "NumSpe_Dee",	
                   "NumRec_Dee",	"NumPhy_Dee",	"ES50_Dee")
Ecological_Data_Global <- Ecological_Data_Global[,analysis.cols]
Ecological_Data_Global  <- Ecological_Data_Global [complete.cases(Ecological_Data_Global ),]

# Calculate the correlation matrix
corr_matrix <- cor(Ecological_Data_Global)

# Create the correlation plot with black font for text
corrplot(corr_matrix, tl.col = "black")

# Function to calculate significance (p-values) matrix
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Calculate the p-value matrix
p_matrix <- cor.mtest(Ecological_Data_Global)

# Set correlation values less than 0.6 to NA
corr_matrix[abs(corr_matrix) < 0.6] <- NA

# Set the tiff device to save the plot
tiff("Correlation_Plot_Deep.tiff", width = 6, height = 4, units = "in", res = 600)

# Create the correlation plot with significant values and formatted coefficients
corrplot(corr_matrix, 
         tl.col = "black",        # Set font color to black
         addCoef.col = "white",   # Add correlation coefficients in white
         number.cex = 0.5,        # Set size of the correlation coefficients
         p.mat = p_matrix,        # Use p-value matrix to highlight significance
         sig.level = 0.05,        # Only show significant correlations (p < 0.05)
         insig = "blank",         # Leave insignificant correlations blank
         na.label = " ",          # Leave NA cells blank (no text or label)
         tl.cex = 0.7,            # Set smaller font size for row/column labels
         number.digits = 1)       # Display correlation coefficients with 1 digit after decimal

# Close the tiff device
dev.off()

#GLMs for number of species, Global data, Deep
deep.numsp.intercept <- glm(NumSpe_Dee ~ 1, family = "poisson", data = Ecological_Data_Global)
summary(deep.numsp.intercept)

deep.numsp.numrec <- glm(NumSpe_Dee ~ NumRec_Dee, family = "poisson", data = Ecological_Data_Global)
summary(deep.numsp.numrec)

deep.numsp.marginsum <- glm(NumSpe_Dee ~ NumRec_Dee + Margin_Sum, family = "poisson", data = Ecological_Data_Global)
summary(deep.numsp.marginsum)

deep.numsp.depth <- glm(NumSpe_Dee ~ NumRec_Dee + Depth_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.numsp.depth)

deep.numsp.humimp <- glm(NumSpe_Dee ~ NumRec_Dee + HumImp_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.numsp.humimp)

deep.numsp.curvel <- glm(NumSpe_Dee ~ NumRec_Dee + CurVel_Bot_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.numsp.curvel)

deep.numsp.nitrate <- glm(NumSpe_Dee ~ NumRec_Dee + Nitrate_Bot_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.numsp.nitrate)

deep.numsp.oxygen <- glm(NumSpe_Dee ~ NumRec_Dee + O2_Bot_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.numsp.oxygen)

deep.numsp.area <- glm(NumSpe_Dee ~ NumRec_Dee + Area_Sum, family = "poisson", data = Ecological_Data_Global)
summary(deep.numsp.area)

deep.numsp.primprod <- glm(NumSpe_Dee ~ NumRec_Dee + PrimProd_Bot_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.numsp.primprod)

deep.numsp.themmean <- glm(NumSpe_Dee ~ NumRec_Dee + ThemM_Bot_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.numsp.themmean)

#deep.numsp.themrange <- glm(NumSpe_Dee ~ NumRec_Dee + ThemR_Bot_Mean, family = "poisson", data = Ecological_Data_Global)
#summary(deep.numsp.themrange)

#Model selection for number of species, Global
deep.numsp.models <- list(Intercept = deep.numsp.intercept,
                            NumRec = deep.numsp.numrec,
                            ConMar = deep.numsp.marginsum,
                            Depth = deep.numsp.depth,
                            HumImp = deep.numsp.humimp,
                            CurVel = deep.numsp.curvel,
                            Nitrate = deep.numsp.nitrate,
                            O2 = deep.numsp.oxygen,
                            Area = deep.numsp.area,
                            PriPro = deep.numsp.primprod,
                            #TheRan = deep.numsp.themrange,
                            TheMean = deep.numsp.themmean)
deep.numsp.aic.df <- data.frame(Model = names(deep.numsp.models),
                                  AIC = sapply(deep.numsp.models, function(x) AICc(x)),
                                  akaike.weights(sapply(deep.numsp.models, function(x) AICc(x))))


deep.numsp.aic.df <- deep.numsp.aic.df[order(deep.numsp.aic.df$AIC),]
deep.numsp.aic.df$Cumulative.Weight <- cumsum(deep.numsp.aic.df$weights)

kable(deep.numsp.aic.df, row.names = FALSE)

#write.csv(deep.numsp.aic.df, file = "global.5.degree.numsp.deep.aic.csv")

#Plots for number of species, Deep
Predicted_NumSp_NumRe_dee <- ggplot(Ecological_Data_Global, aes(x = NumRec_Dee, y = predict(deep.numsp.numrec, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  labs(
    x = "Number of Records",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_NumRe_dee.tiff", plot = Predicted_NumSp_NumRe_dee, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Margin_dee <- ggplot(Ecological_Data_Global, aes(x = Margin_Sum, y = predict(deep.numsp.marginsum, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  labs(
    x = "Continental Margin (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Margin_dee.tiff", plot = Predicted_NumSp_Margin_dee, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Depth_dee <- ggplot(Ecological_Data_Global, aes(x = Depth_Mean, y = predict(deep.numsp.depth, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  labs(
    x = "Depth (m)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Depth_dee.tiff", plot = Predicted_NumSp_Depth_dee, width = 6, height = 4, dpi = 600)

Predicted_NumSp_humimp_dee <- ggplot(Ecological_Data_Global, aes(x = HumImp_Mean, y = predict(deep.numsp.humimp, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  #geom_smooth(method = "glm", formula = y ~ x, color = "black") +  # Add a smooth black line using GLM
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Human Impact",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_humipm_dee.tiff", plot = Predicted_NumSp_humimp_dee, width = 6, height = 4, dpi = 600)

Predicted_NumSp_curvel_dee <- ggplot(Ecological_Data_Global, aes(x = CurVel_Bot_Mean, y = predict(deep.numsp.curvel, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  #geom_smooth(method = "glm", formula = y ~ x, color = "black") +  # Add a smooth black line using GLM
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Current Velocity (m.s-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_curvel_dee.tiff", plot = Predicted_NumSp_curvel_dee, width = 6, height = 4, dpi = 600)

Predicted_NumSp_Nitrate_dee <- ggplot(Ecological_Data_Global, aes(x = Nitrate_Bot_Mean, y = predict(deep.numsp.nitrate, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  #geom_smooth(method = "glm", formula = y ~ x, color = "black") +  # Add a smooth black line using GLM
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Nitrate (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_Nitrate_dee.tiff", plot = Predicted_NumSp_Nitrate_dee, width = 6, height = 4, dpi = 600)

Predicted_NumSp_oxygen_dee <- ggplot(Ecological_Data_Global, aes(x = O2_Bot_Mean, y = predict(deep.numsp.oxygen, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  #geom_smooth(method = "glm", formula = y ~ x, color = "black") +  # Add a smooth black line using GLM
  theme_bw() +  # Use the black and white theme
  labs(
    x = "O2 (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_oxygen_dee.tiff", plot = Predicted_NumSp_oxygen_dee, width = 6, height = 4, dpi = 600)

Predicted_NumSp_area_dee <- ggplot(Ecological_Data_Global, aes(x = Area_Sum, y = predict(deep.numsp.area, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  labs(
    x = "Area (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_area_dee.tiff", plot = Predicted_NumSp_area_dee, width = 6, height = 4, dpi = 600)

Predicted_NumSp_primprod_dee <- ggplot(Ecological_Data_Global, aes(x = PrimProd_Bot_Mean, y = predict(deep.numsp.primprod, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Primary Productivity (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_primprod_dee.tiff", plot = Predicted_NumSp_primprod_dee, width = 6, height = 4, dpi = 600)

Predicted_NumSp_themmean_dee <- ggplot(Ecological_Data_Global, aes(x = ThemM_Bot_Mean, y = predict(deep.numsp.themmean, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Mean (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_themmean_dee.tiff", plot = Predicted_NumSp_themmean_dee, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_themrange_dee <- ggplot(Ecological_Data_Global, aes(x = ThemR_Bot_Mean, y = predict(deep.numsp.themrange, Ecological_Data_Global))) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Range (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_themrange_dee.tiff", plot = Predicted_NumSp_themrange_dee, width = 6, height = 4, dpi = 600)

# ===================================
#                DEEP-ES50-GLM
# ===================================

deep.es50.intercept <- glm(ES50_Dee ~ 1, family = "poisson", data = Ecological_Data_Global)
summary(deep.es50.intercept)

deep.es50.numrec <- glm(ES50_Dee ~ NumRec_Dee, family = "poisson", data = Ecological_Data_Global)
summary(deep.es50.numrec)

deep.es50.marginsum <- glm(ES50_Dee ~ Margin_Sum, family = "poisson", data = Ecological_Data_Global)
summary(deep.es50.marginsum)

deep.es50.depth <- glm(ES50_Dee ~ Depth_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.es50.depth)

deep.es50.humimp <- glm(ES50_Dee ~ HumImp_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.es50.humimp)

deep.es50.curvel <- glm(ES50_Dee ~ CurVel_Bot_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.es50.curvel)

deep.es50.nitrate <- glm(ES50_Dee ~ Nitrate_Bot_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.es50.nitrate)

deep.es50.oxygen <- glm(ES50_Dee ~ O2_Bot_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.es50.oxygen)

deep.es50.area <- glm(ES50_Dee ~ Area_Sum, family = "poisson", data = Ecological_Data_Global)
summary(deep.es50.area)

deep.es50.primprod <- glm(ES50_Dee ~ PrimProd_Bot_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.es50.primprod)

deep.es50.themmean <- glm(ES50_Dee ~ ThemM_Bot_Mean, family = "poisson", data = Ecological_Data_Global)
summary(deep.es50.themmean)

#deep.es50.themrange <- glm(ES50_Dee ~ ThemR_Bot_Mean, family = "poisson", data = Ecological_Data_Global)
#summary(deep.es50.themrange)

#Model selection for number of species, Global
deep.es50.models <- list(Intercept = deep.es50.intercept,
                           NumRec = deep.es50.numrec,
                           ConMar = deep.es50.marginsum,
                           Depth = deep.es50.depth,
                           HumImp = deep.es50.humimp,
                           CurVel = deep.es50.curvel,
                           Nitrate = deep.es50.nitrate,
                           O2 = deep.es50.oxygen,
                           Area = deep.es50.area,
                           PriPro = deep.es50.primprod,
                           #TheRan = deep.es50.themrange,
                           TheMea = deep.es50.themmean)
deep.es50.aic.df <- data.frame(Model = names(deep.es50.models),
                                 AIC = sapply(deep.es50.models, function(x) AICc(x)),
                                 akaike.weights(sapply(deep.es50.models, function(x) AICc(x))))


deep.es50.aic.df <- deep.es50.aic.df[order(deep.es50.aic.df$AIC),]
deep.es50.aic.df$Cumulative.Weight <- cumsum(deep.es50.aic.df$weights)

kable(deep.es50.aic.df, row.names = FALSE)

#write.csv(deep.es50.aic.df, file = "deep.5.degree.es50.deep.aic.csv")

#Plots for ES50, Deep
Predicted_ES50_NumRe_dee <- ggplot(Ecological_Data_Global, aes(x = NumRec_Dee, y = ES50_Dee)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Number of Records",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_NumRe_dee.tiff", plot = Predicted_ES50_NumRe_dee, width = 6, height = 4, dpi = 600)

Predicted_ES50_Margin_dee <- ggplot(Ecological_Data_Global, aes(x = Margin_Sum, y = ES50_Dee)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Continental Margin (km2)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_Margin_dee.tiff", plot = Predicted_ES50_Margin_dee, width = 6, height = 4, dpi = 600)

Predicted_ES50_Depth_dee <- ggplot(Ecological_Data_Global, aes(x = Depth_Mean, y = ES50_Dee)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Depth (m)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_Depth_dee.tiff", plot = Predicted_ES50_Depth_dee, width = 6, height = 4, dpi = 600)

Predicted_ES50_humimp_dee <- ggplot(Ecological_Data_Global, aes(x = HumImp_Mean, y = ES50_Dee)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Human Impact",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_humipm_dee.tiff", plot = Predicted_ES50_humimp_dee, width = 6, height = 4, dpi = 600)

Predicted_ES50_curvel_dee <- ggplot(Ecological_Data_Global, aes(x = CurVel_Bot_Mean, y = ES50_Dee)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Current Velocity (m.s-1)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_curvel_dee.tiff", plot = Predicted_ES50_curvel_dee, width = 6, height = 4, dpi = 600)

Predicted_ES50_Nitrate_dee <- ggplot(Ecological_Data_Global, aes(x = Nitrate_Bot_Mean, y = ES50_Dee)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Nitrate (mmol . m-3)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_Nitrate_dee.tiff", plot = Predicted_ES50_Nitrate_dee, width = 6, height = 4, dpi = 600)

Predicted_ES50_oxygen_dee <- ggplot(Ecological_Data_Global, aes(x = O2_Bot_Mean, y = ES50_Dee)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "O2 (mmol . m-3)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_oxygen_dee.tiff", plot = Predicted_ES50_oxygen_dee, width = 6, height = 4, dpi = 600)

Predicted_ES50_area_dee <- ggplot(Ecological_Data_Global, aes(x = Area_Sum, y = ES50_Dee)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Area (km2)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_area_dee.tiff", plot = Predicted_ES50_area_dee, width = 6, height = 4, dpi = 600)

Predicted_ES50_primprod_dee <- ggplot(Ecological_Data_Global, aes(x = PrimProd_Bot_Mean, y = ES50_Dee)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Primary Productivity (mmol . m-3)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_primprod_dee.tiff", plot = Predicted_ES50_primprod_dee, width = 6, height = 4, dpi = 600)

Predicted_ES50_themmean_dee <- ggplot(Ecological_Data_Global, aes(x = ThemM_Bot_Mean, y = ES50_Dee)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Temperature Mean (ºC)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_themmean_dee.tiff", plot = Predicted_ES50_themmean_dee, width = 6, height = 4, dpi = 600)

#Predicted_ES50_themrange_dee <- ggplot(Ecological_Data_Global, aes(x = ThemR_Bot_Mean, y = ES50_Dee)) +
  geom_smooth(method = "glm", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  scale_y_continuous(limits=c(29, 52), breaks=seq(30, 50, by=5), expand = c(0, 0)) +
  labs(
    x = "Temperature Range (ºC)",  # Shorten the x-axis title
    y = "ES50"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_themrange_dee.tiff", plot = Predicted_ES50_themrange_dee, width = 6, height = 4, dpi = 600)

# ===================================
#                SHALLOW-NUMSPe-GAM
# ===================================

#Species Counts and Environment, hexagons
Ecological_Data_Global_Hex_sp <- read.csv("Ecological_Data_Global_Hex_sp.csv", sep = ";")

summary(Ecological_Data_Global_Hex_sp)
#Ecological_Data_Global_Hex_sp <- write.csv(Ecological_Data_Global_Hex_sp, file="Ecological_Data_Global_Hex_sp.csv")

#Correlation analyses for variables
#select only numeric values
numeric_data <- Ecological_Data_Global_Hex_sp %>% select_if(is.numeric)

#select for all columns
correlation_matrix <- cor(numeric_data, use = "complete.obs", method = "spearman")

# Function to calculate p-values for each pairwise correlation
cor_test_pval <- function(x) {
  mat <- numeric(ncol(x) * ncol(x))
  dim(mat) <- c(ncol(x), ncol(x))
  
  for (i in 1:ncol(x)) {
    for (j in i:ncol(x)) {
      test <- cor.test(x[, i], x[, j], method = "spearman", use = "complete.obs")
      mat[i, j] <- test$p.value
      mat[j, i] <- mat[i, j]
    }
  }
  colnames(mat) <- colnames(x)
  rownames(mat) <- colnames(x)
  mat
}

# Calculate p-values
p_values <- cor_test_pval(numeric_data)

# Set correlations to NA where p > 0.05 (not significant)
correlation_matrix[p_values > 0.05 | abs(correlation_matrix) < 0.7] <- NA

# Plot only significant correlations
corrplot(correlation_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, na.label = " ", tl.cex = 0.4, # Adjust font size with tl.cex
         title = "Significant Correlations (Spearman, p < 0.05)", mar = c(0, 0, 1, 0))

#First we’re going to load in our data and then trim the data frame down to just the columns we need.
analysis.cols <- c("Latitude", "Longitude", "Depth_Mean", "Shelf_Mean", "CurVel_Mean", "HumImp_Mean", "IceCov_Mean", 
                   "Nitrate_Mean", "PrimProd_Mean", "ThemM_mean", "ThemR_mean", "NumRec_sha", "NumSpe_sha", "ES50_sha")
Ecological_Data_Global_Hex_sp <- Ecological_Data_Global_Hex_sp [,analysis.cols]
Ecological_Data_Global_Hex_sp  <- Ecological_Data_Global_Hex_sp [complete.cases(Ecological_Data_Global_Hex_sp),]

# Calculate the correlation matrix
corr_matrix <- cor(Ecological_Data_Global_Hex_sp)

# Create the correlation plot with black font for text
corrplot(corr_matrix, tl.col = "black")

# Function to calculate significance (p-values) matrix
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Calculate the p-value matrix
p_matrix <- cor.mtest(Ecological_Data_Global_Hex_sp)

# Set correlation values less than 0.6 to NA
corr_matrix[abs(corr_matrix) < 0.4] <- NA
p_matrix[abs(corr_matrix) < 0.4] <- NA

# Set the tiff device to save the plot
tiff("Correlation_Plot_Global_GAM_Shallo.tiff", width = 6, height = 4, units = "in", res = 600)

# Create the correlation plot with significant values and formatted coefficients
corrplot(corr_matrix, 
         tl.col = "black",        # Set font color to black
         addCoef.col = "white",   # Add correlation coefficients in white
         number.cex = 0.5,        # Set size of the correlation coefficients
         p.mat = p_matrix,        # Use p-value matrix to highlight significance
         sig.level = 0.05,        # Only show significant correlations (p < 0.05)
         insig = "blank",         # Leave insignificant correlations blank
         na.label = " ",          # Leave NA cells blank (no text or label)
         tl.cex = 0.7,            # Set smaller font size for row/column labels
         number.digits = 1)       # Display correlation coefficients with 1 digit after decimal

# Close the tiff device
dev.off()

#GAMs for number of species, shallow water
shallow.numsp.intercept <- gam(NumSpe_sha ~ 1, data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.numsp.intercept)
summary(shallow.numsp.intercept)

shallow.numsp.latlon <- gam(NumSpe_sha ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_sha), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.numsp.latlon)
summary(shallow.numsp.latlon)

shallow.numsp.depth <- gam(NumSpe_sha ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_sha) + s(Depth_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.numsp.depth)
summary(shallow.numsp.depth)

shallow.numsp.shelf <- gam(NumSpe_sha ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_sha) + s(Shelf_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.numsp.shelf)
summary(shallow.numsp.shelf)

shallow.numsp.current <- gam(NumSpe_sha ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_sha) + s(CurVel_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.numsp.current)
summary(shallow.numsp.current)

shallow.numsp.humimp <- gam(NumSpe_sha ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_sha) + s(HumImp_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.numsp.humimp)
summary(shallow.numsp.humimp)

shallow.numsp.icecover <- gam(NumSpe_sha ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_sha) + s(IceCov_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.numsp.icecover)
summary(shallow.numsp.icecover)

shallow.numsp.nitrate <- gam(NumSpe_sha ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_sha) + s(Nitrate_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.numsp.nitrate)
summary(shallow.numsp.nitrate)

#shallow.numsp.o2 <- gam(NumSpe_sha ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_sha) + s(O2_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(shallow.numsp.o2)
#summary(shallow.numsp.o2)

#shallow.numsp.PhotoActi <- gam(NumSpe_sha ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_sha) + s(PhotoActi_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(shallow.numsp.PhotoActi)
#summary(shallow.numsp.PhotoActi)

shallow.numsp.PrimProd <- gam(NumSpe_sha ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_sha) + s(PrimProd_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.numsp.PrimProd)
summary(shallow.numsp.PrimProd)

shallow.numsp.ThemM <- gam(NumSpe_sha ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_sha) + s(ThemM_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.numsp.ThemM)
summary(shallow.numsp.ThemM)

shallow.numsp.ThemR <- gam(NumSpe_sha ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_sha) + s(ThemR_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.numsp.ThemR)
summary(shallow.numsp.ThemR)

shallow.numsp.env <- gam(NumSpe_sha ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_sha) + s(Depth_Mean) +  s(Shelf_Mean) + s(CurVel_Mean) + s(HumImp_Mean) + s(IceCov_Mean) + s(Nitrate_Mean) + s(PrimProd_Mean) + s(ThemM_mean) + s(ThemR_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.numsp.env)
summary(shallow.numsp.env)

shallow.numsp.models <- list(Intercept = shallow.numsp.intercept,
                             LatLon = shallow.numsp.latlon,
                             Depth = shallow.numsp.depth,
                             ConShe = shallow.numsp.shelf,
                             CurVel = shallow.numsp.current,
                             HumImp = shallow.numsp.humimp,
                             IceCov = shallow.numsp.icecover,
                             Nitrate = shallow.numsp.nitrate,
                             #O2 = shallow.numsp.o2,
                             #PhoAct = shallow.numsp.PhotoActi,
                             PriPro = shallow.numsp.PrimProd,
                             TheMea = shallow.numsp.ThemM,
                             TheRan = shallow.numsp.ThemR,
                             Environment = shallow.numsp.env)
shallow.numsp.aic.df <- data.frame(Model = names(shallow.numsp.models),
                                   AIC = sapply(shallow.numsp.models, function(x) x$aic),
                                   akaike.weights(sapply(shallow.numsp.models, function(x) x$aic)))

shallow.numsp.aic.df <- shallow.numsp.aic.df[order(shallow.numsp.aic.df$AIC),]
shallow.numsp.aic.df$Cumulative.Weight <- cumsum(shallow.numsp.aic.df$weights)

kable(shallow.numsp.aic.df, row.names = FALSE)

#write.csv(shallow.numsp.aic.df, file = "shallow.numsp.aic.GAM.csv")

#Plots for number of species, shallow water
Predicted_NumSp_depth_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Depth_Mean, y = predict(shallow.numsp.depth, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Depth (m)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_depth_sha_gam.tiff", plot = Predicted_NumSp_depth_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_shelf_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Shelf_Mean, y = predict(shallow.numsp.shelf, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Continental Shelf (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_shelf_sha_gam.tiff", plot = Predicted_NumSp_shelf_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_current_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = CurVel_Mean, y = predict(shallow.numsp.current, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Current Velocity (m.s-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_current_sha_gam.tiff", plot = Predicted_NumSp_current_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_humimp_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = HumImp_Mean, y = predict(shallow.numsp.humimp, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Human Impact",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_humimp_sha_gam.tiff", plot = Predicted_NumSp_humimp_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_icecover_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = IceCov_Mean, y = predict(shallow.numsp.icecover, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Ice Cover (fraction)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_icecover_sha_gam.tiff", plot = Predicted_NumSp_icecover_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_nitrate_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Nitrate_Mean, y = predict(shallow.numsp.nitrate, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Nitrate (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_nitrate_sha_gam.tiff", plot = Predicted_NumSp_nitrate_sha_gam, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_o2_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = O2_Mean, y = predict(shallow.numsp.o2, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "O2 (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_o2_sha_gam.tiff", plot = Predicted_NumSp_o2_sha_gam, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_PhotoActi_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = PhotoActi_Mean, y = predict(shallow.numsp.PhotoActi, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Phot. Acti. Radi (E.m-2.day-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_PhotoActi_sha_gam.tiff", plot = Predicted_NumSp_PhotoActi_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_PrimProd_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = PrimProd_Mean, y = predict(shallow.numsp.PrimProd, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Primary Productivity (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_PrimProd_sha_gam.tiff", plot = Predicted_NumSp_PrimProd_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_ThemM_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = ThemM_mean, y = predict(shallow.numsp.ThemM, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Mean (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_ThemM_sha_gam.tiff", plot = Predicted_NumSp_ThemM_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_ThemR_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = ThemR_mean, y = predict(shallow.numsp.ThemR, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Range (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_ThemR_sha_gam.tiff", plot = Predicted_NumSp_ThemR_sha_gam, width = 6, height = 4, dpi = 600)

# Generate predictions for the model
Ecological_Data_Global_Hex_sp$Prediction <- predict(shallow.numsp.env, type = "response")

# ===================================
#                SHALLOW-ES50-GAM
# ===================================

shallow.ES50.intercept <- gam(ES50_sha ~ 1, data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.ES50.intercept)
summary(shallow.ES50.intercept)

shallow.ES50.latlon <- gam(ES50_sha ~ s(Latitude, Longitude, bs = "sos"), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.ES50.latlon)
summary(shallow.ES50.latlon)

shallow.ES50.depth <- gam(ES50_sha ~ s(Latitude, Longitude, bs = "sos") + s(Depth_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.ES50.depth)
summary(shallow.ES50.depth)

shallow.ES50.shelf <- gam(ES50_sha ~ s(Latitude, Longitude, bs = "sos")  + s(Shelf_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.ES50.shelf)
summary(shallow.ES50.shelf)

shallow.ES50.current <- gam(ES50_sha ~ s(Latitude, Longitude, bs = "sos") + s(CurVel_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.ES50.current)
summary(shallow.ES50.current)

shallow.ES50.humimp <- gam(ES50_sha ~ s(Latitude, Longitude, bs = "sos") + s(HumImp_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.ES50.humimp)
summary(shallow.ES50.humimp)

shallow.ES50.icecover <- gam(ES50_sha ~ s(Latitude, Longitude, bs = "sos") + s(IceCov_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.ES50.icecover)
summary(shallow.ES50.icecover)

shallow.ES50.nitrate <- gam(ES50_sha ~ s(Latitude, Longitude, bs = "sos") + s(Nitrate_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.ES50.nitrate)
summary(shallow.ES50.nitrate)

#shallow.ES50.o2 <- gam(ES50_sha ~ s(Latitude, Longitude, bs = "sos") + s(O2_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(shallow.ES50.o2)
#summary(shallow.ES50.o2)

#shallow.ES50.PhotoActi <- gam(ES50_sha ~ s(Latitude, Longitude, bs = "sos") + s(PhotoActi_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(shallow.ES50.PhotoActi)
#summary(shallow.ES50.PhotoActi)

shallow.ES50.PrimProd <- gam(ES50_sha ~ s(Latitude, Longitude, bs = "sos") + s(PrimProd_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.ES50.PrimProd)
summary(shallow.ES50.PrimProd)

shallow.ES50.ThemM <- gam(ES50_sha ~ s(Latitude, Longitude, bs = "sos") + s(ThemM_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.ES50.ThemM)
summary(shallow.ES50.ThemM)

shallow.ES50.ThemR <- gam(ES50_sha ~ s(Latitude, Longitude, bs = "sos") + s(ThemR_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.ES50.ThemR)
summary(shallow.ES50.ThemR)

shallow.ES50.env <- gam(ES50_sha ~ s(Latitude, Longitude, bs = "sos") + s(Depth_Mean) + s(Shelf_Mean) + s(CurVel_Mean) + s(HumImp_Mean) + s(IceCov_Mean) + s(Nitrate_Mean) + s(PrimProd_Mean) + s(ThemM_mean) + s(ThemR_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(shallow.ES50.env)
summary(shallow.ES50.env)

shallow.ES50.models <- list(Intercept = shallow.ES50.intercept,
                             LatLon = shallow.ES50.latlon,
                             Depth = shallow.ES50.depth,
                             ConShe= shallow.ES50.shelf,
                             CurVel = shallow.ES50.current,
                             HumImp = shallow.ES50.humimp,
                             IceCov = shallow.ES50.icecover,
                             Nitrate = shallow.ES50.nitrate,
                             #O2 = shallow.ES50.o2,
                             #PhoAct = shallow.ES50.PhotoActi,
                             PriPro = shallow.ES50.PrimProd,
                             TheMea = shallow.ES50.ThemM,
                             TheRan = shallow.ES50.ThemR,
                             Environment = shallow.numsp.env)
shallow.ES50.aic.df <- data.frame(Model = names(shallow.ES50.models),
                                   AIC = sapply(shallow.ES50.models, function(x) x$aic),
                                   akaike.weights(sapply(shallow.ES50.models, function(x) x$aic)))

shallow.ES50.aic.df <- shallow.ES50.aic.df[order(shallow.ES50.aic.df$AIC),]
shallow.ES50.aic.df$Cumulative.Weight <- cumsum(shallow.ES50.aic.df$weights)

kable(shallow.ES50.aic.df, row.names = FALSE)

#write.csv(shallow.ES50.aic.df, file = "shallow.ES50.aic.GAM.csv")

#Plots for ES50, shallow water
Predicted_ES50_depth_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Depth_Mean, y = predict(shallow.ES50.depth, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Depth (m)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_depth_sha_gam.tiff", plot = Predicted_ES50_depth_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_shelf_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Shelf_Mean, y = predict(shallow.ES50.shelf, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Continental Shelf (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_shelf_sha_gam.tiff", plot = Predicted_ES50_shelf_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_current_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = CurVel_Mean, y = predict(shallow.ES50.current, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Current Velocity (m.s-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_current_sha_gam.tiff", plot = Predicted_ES50_current_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_humimp_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = HumImp_Mean, y = predict(shallow.ES50.humimp, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Human Impact",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_humimp_sha_gam.tiff", plot = Predicted_ES50_humimp_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_icecover_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = IceCov_Mean, y = predict(shallow.ES50.icecover, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Ice Cover (fraction)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_icecover_sha_gam.tiff", plot = Predicted_ES50_icecover_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_nitrate_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Nitrate_Mean, y = predict(shallow.ES50.nitrate, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Nitrate (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_nitrate_sha_gam.tiff", plot = Predicted_ES50_nitrate_sha_gam, width = 6, height = 4, dpi = 600)

#Predicted_ES50_o2_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = O2_Mean, y = predict(shallow.ES50.o2, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "O2 (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_o2_sha_gam.tiff", plot = Predicted_ES50_o2_sha_gam, width = 6, height = 4, dpi = 600)

#Predicted_ES50_PhotoActi_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = PhotoActi_Mean, y = predict(shallow.ES50.PhotoActi, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Phot. Acti. Radi (E.m-2.day-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_PhotoActi_sha_gam.tiff", plot = Predicted_ES50_PhotoActi_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_PrimProd_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = PrimProd_Mean, y = predict(shallow.ES50.PrimProd, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Primary Productivity (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_PrimProd_sha_gam.tiff", plot = Predicted_ES50_PrimProd_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_ThemM_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = ThemM_mean, y = predict(shallow.ES50.ThemM, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Mean (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_ThemM_sha_gam.tiff", plot = Predicted_ES50_ThemM_sha_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_ThemR_sha_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = ThemR_mean, y = predict(shallow.ES50.ThemR, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Range (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_ThemR_sha_gam.tiff", plot = Predicted_ES50_ThemR_sha_gam, width = 6, height = 4, dpi = 600)

# ===================================
#                MESO-NUMSPe-GAM
# ===================================

#Species Counts and Environment, hexagons
Ecological_Data_Global_Hex_sp <- read.csv("Ecological_Data_Global_Hex_sp.csv", sep = ";")

#select only numeric values
numeric_data <- Ecological_Data_Global_Hex_sp %>% select_if(is.numeric)

#First we’re going to load in our data and then trim the data frame down to just the columns we need.
analysis.cols <- c("Latitude", "Longitude", "Depth_Mean", "Margin_Sum", "CurVel_Mean", "HumImp_Mean", 
                   "Nitrate_Mean", "PrimProd_Mean", "ThemM_mean", "ThemR_mean", "NumRec_mes", "NumSpe_mes", "ES50_mes")
Ecological_Data_Global_Hex_sp <- Ecological_Data_Global_Hex_sp [,analysis.cols]
Ecological_Data_Global_Hex_sp  <- Ecological_Data_Global_Hex_sp [complete.cases(Ecological_Data_Global_Hex_sp),]

# Calculate the correlation matrix
corr_matrix <- cor(Ecological_Data_Global_Hex_sp)

# Create the correlation plot with black font for text
corrplot(corr_matrix, tl.col = "black")

# Function to calculate significance (p-values) matrix
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Calculate the p-value matrix
p_matrix <- cor.mtest(Ecological_Data_Global_Hex_sp)

# Set correlation values less than 0.6 to NA
corr_matrix[abs(corr_matrix) < 0.4] <- NA
p_matrix[abs(corr_matrix) < 0.4] <- NA

# Set the tiff device to save the plot
tiff("Correlation_Plot_Global_GAM_Mes.tiff", width = 6, height = 4, units = "in", res = 600)

# Create the correlation plot with significant values and formatted coefficients
corrplot(corr_matrix, 
         tl.col = "black",        # Set font color to black
         addCoef.col = "white",   # Add correlation coefficients in white
         number.cex = 0.5,        # Set size of the correlation coefficients
         p.mat = p_matrix,        # Use p-value matrix to highlight significance
         sig.level = 0.05,        # Only show significant correlations (p < 0.05)
         insig = "blank",         # Leave insignificant correlations blank
         na.label = " ",          # Leave NA cells blank (no text or label)
         tl.cex = 0.7,            # Set smaller font size for row/column labels
         number.digits = 1)       # Display correlation coefficients with 1 digit after decimal

# Close the tiff device
dev.off()

#GAMs for number of species, shallow water
meso.numsp.intercept <- gam(NumSpe_mes ~ 1, data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.numsp.intercept)
summary(meso.numsp.intercept)

meso.numsp.latlon <- gam(NumSpe_mes ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_mes), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.numsp.latlon)
summary(meso.numsp.latlon)

meso.numsp.depth <- gam(NumSpe_mes ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_mes) + s(Depth_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.numsp.depth)
summary(meso.numsp.depth)

meso.numsp.margin <- gam(NumSpe_mes ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_mes) + s(Margin_Sum), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.numsp.margin)
summary(meso.numsp.margin)

meso.numsp.current <- gam(NumSpe_mes ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_mes) + s(CurVel_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.numsp.current)
summary(meso.numsp.current)

meso.numsp.humimp <- gam(NumSpe_mes ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_mes) + s(HumImp_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.numsp.humimp)
summary(meso.numsp.humimp)

#meso.numsp.icecover <- gam(NumSpe_mes ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_mes) + s(IceCov_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(meso.numsp.icecover)
#summary(meso.numsp.icecover)

meso.numsp.nitrate <- gam(NumSpe_mes ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_mes) + s(Nitrate_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.numsp.nitrate)
summary(meso.numsp.nitrate)

#meso.numsp.o2 <- gam(NumSpe_mes ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_mes) + s(O2_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(meso.numsp.o2)
#summary(meso.numsp.o2)

#meso.numsp.PhotoActi <- gam(NumSpe_mes ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_mes) + s(PhotoActi_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(meso.numsp.PhotoActi)
#summary(meso.numsp.PhotoActi)

meso.numsp.PrimProd <- gam(NumSpe_mes ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_mes) + s(PrimProd_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.numsp.PrimProd)
summary(meso.numsp.PrimProd)

meso.numsp.ThemM <- gam(NumSpe_mes ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_mes) + s(ThemM_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.numsp.ThemM)
summary(meso.numsp.ThemM)

meso.numsp.ThemR <- gam(NumSpe_mes ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_mes) + s(ThemR_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.numsp.ThemR)
summary(meso.numsp.ThemR)

meso.numsp.env <- gam(NumSpe_mes ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_mes) + s(Depth_Mean) +  s(Margin_Sum) + s(CurVel_Mean) + s(HumImp_Mean) + s(Nitrate_Mean) + s(PrimProd_Mean) + s(ThemM_mean) + s(ThemR_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.numsp.env)
summary(meso.numsp.env)

meso.numsp.models <- list(Intercept = meso.numsp.intercept,
                          LatLon = meso.numsp.latlon,
                          Depth = meso.numsp.depth,
                          ConMAr = meso.numsp.margin,
                          CurVel = meso.numsp.current,
                          HumImp = meso.numsp.humimp,
                          #IceCov = meso.numsp.icecover,
                          Nitrate = meso.numsp.nitrate,
                          #O2 = meso.numsp.o2,
                          #PhoAct = meso.numsp.PhotoActi,
                          PriPro = meso.numsp.PrimProd,
                          TheMea = meso.numsp.ThemM,
                          TheRan = meso.numsp.ThemR,
                          Environment = meso.numsp.env)
meso.numsp.aic.df <- data.frame(Model = names(meso.numsp.models),
                                   AIC = sapply(meso.numsp.models, function(x) x$aic),
                                   akaike.weights(sapply(meso.numsp.models, function(x) x$aic)))

meso.numsp.aic.df <- meso.numsp.aic.df[order(meso.numsp.aic.df$AIC),]
meso.numsp.aic.df$Cumulative.Weight <- cumsum(meso.numsp.aic.df$weights)

kable(meso.numsp.aic.df, row.names = FALSE)

#write.csv(meso.numsp.aic.df, file = "meso.numsp.aic.GAM.csv")

#Plots for number of species, Meso
Predicted_NumSp_depth_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Depth_Mean, y = predict(meso.numsp.depth, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Depth (m)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_depth_mes_gam.tiff", plot = Predicted_NumSp_depth_mes_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_margin_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Margin_Sum, y = predict(meso.numsp.margin, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Continental Margin (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_margin_mes_gam.tiff", plot = Predicted_NumSp_margin_mes_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_current_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = CurVel_Mean, y = predict(meso.numsp.current, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Current Velocity (m.s-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_current_mes_gam.tiff", plot = Predicted_NumSp_current_mes_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_humimp_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = HumImp_Mean, y = predict(meso.numsp.humimp, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Human Impact",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_humimp_mes_gam.tiff", plot = Predicted_NumSp_humimp_mes_gam, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_icecover_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = IceCov_Mean, y = predict(meso.numsp.icecover, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Ice Cover (fraction)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_icecover_mes_gam.tiff", plot = Predicted_NumSp_icecover_mes_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_nitrate_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Nitrate_Mean, y = predict(meso.numsp.nitrate, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Nitrate (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_nitrate_mes_gam.tiff", plot = Predicted_NumSp_nitrate_mes_gam, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_o2_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = O2_Mean, y = predict(meso.numsp.o2, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "O2 (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_o2_mes_gam.tiff", plot = Predicted_NumSp_o2_mes_gam, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_PhotoActi_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = PhotoActi_Mean, y = predict(meso.numsp.PhotoActi, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Phot. Acti. Radi (E.m-2.day-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_PhotoActi_mes_gam.tiff", plot = Predicted_NumSp_PhotoActi_mes_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_PrimProd_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = PrimProd_Mean, y = predict(meso.numsp.PrimProd, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Primary Productivity (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_PrimProd_mes_gam.tiff", plot = Predicted_NumSp_PrimProd_mes_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_ThemM_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = ThemM_mean, y = predict(meso.numsp.ThemM, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Mean (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_ThemM_mes_gam.tiff", plot = Predicted_NumSp_ThemM_mes_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_ThemR_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = ThemR_mean, y = predict(meso.numsp.ThemR, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Range (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_ThemR_mes_gam.tiff", plot = Predicted_NumSp_ThemR_mes_gam, width = 6, height = 4, dpi = 600)

# ===================================
#                MESO-ES50-GAM
# ===================================

meso.ES50.intercept <- gam(ES50_mes ~ 1, data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.ES50.intercept)
summary(meso.ES50.intercept)

meso.ES50.latlon <- gam(ES50_mes ~ s(Latitude, Longitude, bs = "sos"), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.ES50.latlon)
summary(meso.ES50.latlon)

meso.ES50.depth <- gam(ES50_mes ~ s(Latitude, Longitude, bs = "sos") + s(Depth_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.ES50.depth)
summary(meso.ES50.depth)

meso.ES50.margin <- gam(ES50_mes ~ s(Latitude, Longitude, bs = "sos") + s(Margin_Sum), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.ES50.margin)
summary(meso.ES50.margin)

meso.ES50.current <- gam(ES50_mes ~ s(Latitude, Longitude, bs = "sos") + s(CurVel_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.ES50.current)
summary(meso.ES50.current)

meso.ES50.humimp <- gam(ES50_mes ~ s(Latitude, Longitude, bs = "sos") + s(HumImp_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.ES50.humimp)
summary(meso.ES50.humimp)

#meso.ES50.icecover <- gam(ES50_mes ~ s(Latitude, Longitude, bs = "sos") + s(IceCov_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(meso.ES50.icecover)
#summary(meso.ES50.icecover)

meso.ES50.nitrate <- gam(ES50_mes ~ s(Latitude, Longitude, bs = "sos") + s(Nitrate_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.ES50.nitrate)
summary(meso.ES50.nitrate)

#meso.ES50.o2 <- gam(ES50_mes ~ s(Latitude, Longitude, bs = "sos") + s(O2_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(meso.ES50.o2)
#summary(meso.ES50.o2)

#meso.ES50.PhotoActi <- gam(ES50_mes ~ s(Latitude, Longitude, bs = "sos") + s(PhotoActi_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(meso.ES50.PhotoActi)
#summary(meso.ES50.PhotoActi)

meso.ES50.PrimProd <- gam(ES50_mes ~ s(Latitude, Longitude, bs = "sos") + s(PrimProd_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.ES50.PrimProd)
summary(meso.ES50.PrimProd)

meso.ES50.ThemM <- gam(ES50_mes ~ s(Latitude, Longitude, bs = "sos") + s(ThemM_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.ES50.ThemM)
summary(meso.ES50.ThemM)

meso.ES50.ThemR <- gam(ES50_mes ~ s(Latitude, Longitude, bs = "sos") + s(ThemR_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.ES50.ThemR)
summary(meso.ES50.ThemR)

meso.ES50.env <- gam(ES50_mes ~ s(Latitude, Longitude, bs = "sos") + s(Depth_Mean) +  s(Margin_Sum) + s(CurVel_Mean) + s(HumImp_Mean) + s(Nitrate_Mean) + s(PrimProd_Mean) + s(ThemM_mean) + s(ThemR_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(meso.ES50.env)
summary(meso.ES50.env)

meso.ES50.models <- list(Intercept = meso.ES50.intercept,
                            LatLon = meso.ES50.latlon,
                            Depth = meso.ES50.depth,
                            ConMar = meso.ES50.margin,
                            CurVel = meso.ES50.current,
                            HumImp = meso.ES50.humimp,
                            #IceCov = meso.ES50.icecover,
                            Nitrate = meso.ES50.nitrate,
                            #O2 = meso.ES50.o2,
                            #PhoAct = meso.ES50.PhotoActi,
                            PriPro = meso.ES50.PrimProd,
                            TheMea = meso.ES50.ThemM,
                            TheRan = meso.ES50.ThemR,
                            Environment = meso.ES50.env)
meso.ES50.aic.df <- data.frame(Model = names(meso.ES50.models),
                                  AIC = sapply(meso.ES50.models, function(x) x$aic),
                                  akaike.weights(sapply(meso.ES50.models, function(x) x$aic)))

meso.ES50.aic.df <- meso.ES50.aic.df[order(meso.ES50.aic.df$AIC),]
meso.ES50.aic.df$Cumulative.Weight <- cumsum(meso.ES50.aic.df$weights)

kable(meso.ES50.aic.df, row.names = FALSE)

#write.csv(meso.ES50.aic.df, file = "meso.ES50.aic.GAM.csv")

#Plots for Es50, Meso
Predicted_ES50_depth_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Depth_Mean, y = predict(meso.ES50.depth, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Depth (m)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_depth_mes_gam.tiff", plot = Predicted_ES50_depth_mes_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_margin_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Margin_Sum, y = predict(meso.ES50.margin, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Continental Margin (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_margin_mes_gam.tiff", plot = Predicted_ES50_margin_mes_gam, width = 6, height = 4, dpi = 600)


Predicted_ES50_current_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = CurVel_Mean, y = predict(meso.ES50.current, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Current Velocity (m.s-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_current_mes_gam.tiff", plot = Predicted_ES50_current_mes_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_humimp_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = HumImp_Mean, y = predict(meso.ES50.humimp, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Human Impact",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_humimp_mes_gam.tiff", plot = Predicted_ES50_humimp_mes_gam, width = 6, height = 4, dpi = 600)

#Predicted_ES50_icecover_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = IceCov_Mean, y = predict(meso.ES50.icecover, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Ice Cover (fraction)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_icecover_mes_gam.tiff", plot = Predicted_ES50_icecover_mes_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_nitrate_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Nitrate_Mean, y = predict(meso.ES50.nitrate, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Nitrate (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_nitrate_mes_gam.tiff", plot = Predicted_ES50_nitrate_mes_gam, width = 6, height = 4, dpi = 600)

#Predicted_ES50_o2_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = O2_Mean, y = predict(meso.ES50.o2, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "O2 (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_o2_mes_gam.tiff", plot = Predicted_ES50_o2_mes_gam, width = 6, height = 4, dpi = 600)

#Predicted_ES50_PhotoActi_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = PhotoActi_Mean, y = predict(meso.ES50.PhotoActi, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Phot. Acti. Radi (E.m-2.day-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_PhotoActi_mes_gam.tiff", plot = Predicted_ES50_PhotoActi_mes_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_PrimProd_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = PrimProd_Mean, y = predict(meso.ES50.PrimProd, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Primary Productivity (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_PrimProd_mes_gam.tiff", plot = Predicted_ES50_PrimProd_mes_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_ThemM_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = ThemM_mean, y = predict(meso.ES50.ThemM, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Mean (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_ThemM_mes_gam.tiff", plot = Predicted_ES50_ThemM_mes_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_ThemR_mes_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = ThemR_mean, y = predict(meso.ES50.ThemR, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Range (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_ThemR_mes_gam.tiff", plot = Predicted_ES50_ThemR_mes_gam, width = 6, height = 4, dpi = 600)

# ===================================
#                DEEP-NUMSPe-GAM
# ===================================

#Species Counts and Environment, hexagons
Ecological_Data_Global_Hex_sp <- read.csv("Ecological_Data_Global_Hex_sp.csv", sep = ";")

summary(Ecological_Data_Global_Hex_sp)
#Ecological_Data_Global_Hex_sp <- write.csv(Ecological_Data_Global_Hex_sp, file="Ecological_Data_Global_Hex_sp.csv")

#Correlation analyses for variables
#select only numeric values
numeric_data <- Ecological_Data_Global_Hex_sp %>% select_if(is.numeric)

#select for all columns
correlation_matrix <- cor(numeric_data, use = "complete.obs", method = "spearman")

# Function to calculate p-values for each pairwise correlation
cor_test_pval <- function(x) {
  mat <- numeric(ncol(x) * ncol(x))
  dim(mat) <- c(ncol(x), ncol(x))
  
  for (i in 1:ncol(x)) {
    for (j in i:ncol(x)) {
      test <- cor.test(x[, i], x[, j], method = "spearman", use = "complete.obs")
      mat[i, j] <- test$p.value
      mat[j, i] <- mat[i, j]
    }
  }
  colnames(mat) <- colnames(x)
  rownames(mat) <- colnames(x)
  mat
}

# Calculate p-values
p_values <- cor_test_pval(numeric_data)

# Set correlations to NA where p > 0.05 (not significant)
correlation_matrix[p_values > 0.05 | abs(correlation_matrix) < 0.7] <- NA

# Plot only significant correlations
corrplot(correlation_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, na.label = " ", tl.cex = 0.4, # Adjust font size with tl.cex
         title = "Significant Correlations (Spearman, p < 0.05)", mar = c(0, 0, 1, 0))

#First we’re going to load in our data and then trim the data frame down to just the columns we need.
analysis.cols <- c("Latitude", "Longitude", "Depth_Mean", "Margin_Sum", "CurVel_Bot_Mean", "HumImp_Mean", 
                   "Nitrate_Bot_Mean", "O2_Bot_Mean", "PrimProd_Bot_Mean", "ThemM_Bot_mean", "NumRec_dee", 
                   "NumSpe_dee", "ES50_dee")
Ecological_Data_Global_Hex_sp <- Ecological_Data_Global_Hex_sp [,analysis.cols]
Ecological_Data_Global_Hex_sp  <- Ecological_Data_Global_Hex_sp [complete.cases(Ecological_Data_Global_Hex_sp),]

# Calculate the correlation matrix
corr_matrix <- cor(Ecological_Data_Global_Hex_sp)

# Create the correlation plot with black font for text
corrplot(corr_matrix, tl.col = "black")

# Function to calculate significance (p-values) matrix
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Calculate the p-value matrix
p_matrix <- cor.mtest(Ecological_Data_Global_Hex_sp)

# Set correlation values less than 0.6 to NA
corr_matrix[abs(corr_matrix) < 0.4] <- NA
p_matrix[abs(corr_matrix) < 0.4] <- NA

# Set the tiff device to save the plot
tiff("Correlation_Plot_Global_GAM_deep.tiff", width = 6, height = 4, units = "in", res = 600)

# Create the correlation plot with significant values and formatted coefficients
corrplot(corr_matrix, 
         tl.col = "black",        # Set font color to black
         addCoef.col = "white",   # Add correlation coefficients in white
         number.cex = 0.5,        # Set size of the correlation coefficients
         p.mat = p_matrix,        # Use p-value matrix to highlight significance
         sig.level = 0.05,        # Only show significant correlations (p < 0.05)
         insig = "blank",         # Leave insignificant correlations blank
         na.label = " ",          # Leave NA cells blank (no text or label)
         tl.cex = 0.7,            # Set smaller font size for row/column labels
         number.digits = 1)       # Display correlation coefficients with 1 digit after decimal

# Close the tiff device
dev.off()

#GAMs for number of species, deep 
deep.numsp.intercept <- gam(NumSpe_dee ~ 1, data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.numsp.intercept)
summary(deep.numsp.intercept)

deep.numsp.latlon <- gam(NumSpe_dee ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_dee), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.numsp.latlon)
summary(deep.numsp.latlon)

deep.numsp.depth <- gam(NumSpe_dee ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_dee) + s(Depth_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.numsp.depth)
summary(deep.numsp.depth)

deep.numsp.margin <- gam(NumSpe_dee ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_dee) + s(Margin_Sum), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.numsp.margin)
summary(deep.numsp.margin)

deep.numsp.current <- gam(NumSpe_dee ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_dee) + s(CurVel_Bot_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.numsp.current)
summary(deep.numsp.current)

deep.numsp.humimp <- gam(NumSpe_dee ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_dee) + s(HumImp_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.numsp.humimp)
summary(deep.numsp.humimp)

#deep.numsp.icecover <- gam(NumSpe_dee ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_dee) + s(IceCov_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(deep.numsp.icecover)
#summary(deep.numsp.icecover)

deep.numsp.nitrate <- gam(NumSpe_dee ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_dee) + s(Nitrate_Bot_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.numsp.nitrate)
summary(deep.numsp.nitrate)

deep.numsp.o2 <- gam(NumSpe_dee ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_dee) + s(O2_Bot_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.numsp.o2)
summary(deep.numsp.o2)

deep.numsp.PrimProd <- gam(NumSpe_dee ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_dee) + s(PrimProd_Bot_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.numsp.PrimProd)
summary(deep.numsp.PrimProd)

deep.numsp.ThemM <- gam(NumSpe_dee ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_dee) + s(ThemM_Bot_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.numsp.ThemM)
summary(deep.numsp.ThemM)

#deep.numsp.ThemR <- gam(NumSpe_dee ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_dee) + s(ThemR_Bot_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(deep.numsp.ThemR)
#summary(deep.numsp.ThemR)

deep.numsp.env <- gam(NumSpe_dee ~ s(Latitude, Longitude, bs = "sos") + s(NumRec_dee) + s(Depth_Mean) +  s(Margin_Sum) + s(CurVel_Bot_Mean) + s(HumImp_Mean) + s(Nitrate_Bot_Mean) + s(O2_Bot_Mean) + s(PrimProd_Bot_Mean) + s(ThemM_Bot_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.numsp.env)
summary(deep.numsp.env)

deep.numsp.models <- list(Intercept = deep.numsp.intercept,
                          LatLon = deep.numsp.latlon,
                          Depth = deep.numsp.depth,
                          ConMar = deep.numsp.margin,
                          CurVel = deep.numsp.current,
                          HumImp = deep.numsp.humimp,
                          #IceCov = deep.numsp.icecover,
                          Nitrate = deep.numsp.nitrate,
                          O2 = deep.numsp.o2,
                          PriPro = deep.numsp.PrimProd,
                          ThemMea = deep.numsp.ThemM,
                          #TheRan = deep.numsp.ThemR,
                          Environment = deep.numsp.env)
deep.numsp.aic.df <- data.frame(Model = names(deep.numsp.models),
                                   AIC = sapply(deep.numsp.models, function(x) x$aic),
                                   akaike.weights(sapply(deep.numsp.models, function(x) x$aic)))

deep.numsp.aic.df <- deep.numsp.aic.df[order(deep.numsp.aic.df$AIC),]
deep.numsp.aic.df$Cumulative.Weight <- cumsum(deep.numsp.aic.df$weights)

kable(deep.numsp.aic.df, row.names = FALSE)

#write.csv(deep.numsp.aic.df, file = "deep.numsp.aic.GAM.csv")

#Plots for number of species, deep water
Predicted_NumSp_depth_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Depth_Mean, y = predict(deep.numsp.depth, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Depth (m)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_depth_dee_gam.tiff", plot = Predicted_NumSp_depth_dee_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_margin_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Margin_Sum, y = predict(deep.numsp.margin, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Continental Margin (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_margin_dee_gam.tiff", plot = Predicted_NumSp_margin_dee_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_current_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = CurVel_Bot_Mean, y = predict(deep.numsp.current, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Current Velocity (m.s-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_current_dee_gam.tiff", plot = Predicted_NumSp_current_dee_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_humimp_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = HumImp_Mean, y = predict(deep.numsp.humimp, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Human Impact",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_humimp_dee_gam.tiff", plot = Predicted_NumSp_humimp_dee_gam, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_icecover_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = IceCov_Mean, y = predict(deep.numsp.icecover, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Ice Cover (fraction)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_icecover_dee_gam.tiff", plot = Predicted_NumSp_icecover_dee_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_nitrate_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Nitrate_Bot_Mean, y = predict(deep.numsp.nitrate, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Nitrate (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_nitrate_dee_gam.tiff", plot = Predicted_NumSp_nitrate_dee_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_o2_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = O2_Bot_Mean, y = predict(deep.numsp.o2, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "O2 (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_o2_dee_gam.tiff", plot = Predicted_NumSp_o2_dee_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_PrimProd_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = PrimProd_Bot_Mean, y = predict(deep.numsp.PrimProd, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Primary Productivity (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_PrimProd_dee_gam.tiff", plot = Predicted_NumSp_PrimProd_dee_gam, width = 6, height = 4, dpi = 600)

Predicted_NumSp_ThemM_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = ThemM_Bot_mean, y = predict(deep.numsp.ThemM, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Mean (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_NumSp_ThemM_dee_gam.tiff", plot = Predicted_NumSp_ThemM_dee_gam, width = 6, height = 4, dpi = 600)

#Predicted_NumSp_ThemR_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = ThemR_Bot_mean, y = predict(deep.numsp.ThemR, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Range (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_NumSp_ThemR_dee_gam.tiff", plot = Predicted_NumSp_ThemR_dee_gam, width = 6, height = 4, dpi = 600)

# ===================================
#                DEEP-ES50-GAM
# ===================================

#GAMs for ES50, deep water
deep.ES50.intercept <- gam(ES50_dee ~ 1, data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.ES50.intercept)
summary(deep.ES50.intercept)

deep.ES50.latlon <- gam(ES50_dee ~ s(Latitude, Longitude, bs = "sos"), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.ES50.latlon)
summary(deep.ES50.latlon)

deep.ES50.depth <- gam(ES50_dee ~ s(Latitude, Longitude, bs = "sos") + s(Depth_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.ES50.depth)
summary(deep.ES50.depth)

deep.ES50.margin <- gam(ES50_dee ~ s(Latitude, Longitude, bs = "sos") + s(Margin_Sum), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.ES50.margin)
summary(deep.ES50.margin)

deep.ES50.current <- gam(ES50_dee ~ s(Latitude, Longitude, bs = "sos") + s(CurVel_Bot_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.ES50.current)
summary(deep.ES50.current)

deep.ES50.humimp <- gam(ES50_dee ~ s(Latitude, Longitude, bs = "sos") + s(HumImp_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.ES50.humimp)
summary(deep.ES50.humimp)

#deep.ES50.icecover <- gam(ES50_dee ~ s(Latitude, Longitude, bs = "sos") + s(IceCov_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(deep.ES50.icecover)
#summary(deep.ES50.icecover)

deep.ES50.nitrate <- gam(ES50_dee ~ s(Latitude, Longitude, bs = "sos") + s(Nitrate_Bot_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.ES50.nitrate)
summary(deep.ES50.nitrate)

deep.ES50.o2 <- gam(ES50_dee ~ s(Latitude, Longitude, bs = "sos") + s(O2_Bot_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.ES50.o2)
summary(deep.ES50.o2)

deep.ES50.PrimProd <- gam(ES50_dee ~ s(Latitude, Longitude, bs = "sos") + s(PrimProd_Bot_Mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.ES50.PrimProd)
summary(deep.ES50.PrimProd)

deep.ES50.ThemM <- gam(ES50_dee ~ s(Latitude, Longitude, bs = "sos") + s(ThemM_Bot_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.ES50.ThemM)
summary(deep.ES50.ThemM)

#deep.ES50.ThemR <- gam(ES50_dee ~ s(Latitude, Longitude, bs = "sos") + s(ThemR_Bot_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
#gam.check(deep.ES50.ThemR)
#summary(deep.ES50.ThemR)

deep.ES50.env <- gam(ES50_dee ~ s(Latitude, Longitude, bs = "sos") + s(Depth_Mean) +  s(Margin_Sum) + s(CurVel_Bot_Mean) + s(HumImp_Mean) + s(Nitrate_Bot_Mean) + s(O2_Bot_Mean) + s(PrimProd_Bot_Mean) + s(ThemM_Bot_mean), data = Ecological_Data_Global_Hex_sp, family = "nb", method = "REML", select = TRUE)
gam.check(deep.ES50.env)
summary(deep.ES50.env)

deep.ES50.models <- list(Intercept = deep.ES50.intercept,
                          LatLon = deep.ES50.latlon,
                          Depth = deep.ES50.depth,
                          ConMAr = deep.ES50.margin,
                          CUrVel = deep.ES50.current,
                          HumImp = deep.ES50.humimp,
                          #IceCov = deep.ES50.icecover,
                          Nitrate = deep.ES50.nitrate,
                          O2 = deep.ES50.o2,
                          PriPro = deep.ES50.PrimProd,
                          TheMea = deep.ES50.ThemM,
                          #TheRan = deep.ES50.ThemR,
                          Environment = deep.ES50.env)
deep.ES50.aic.df <- data.frame(Model = names(deep.ES50.models),
                                AIC = sapply(deep.ES50.models, function(x) x$aic),
                                akaike.weights(sapply(deep.ES50.models, function(x) x$aic)))

deep.ES50.aic.df <- deep.ES50.aic.df[order(deep.ES50.aic.df$AIC),]
deep.ES50.aic.df$Cumulative.Weight <- cumsum(deep.ES50.aic.df$weights)

kable(deep.ES50.aic.df, row.names = FALSE)

#write.csv(deep.ES50.aic.df, file = "deep.ES50.aic.GAM.csv")

#Plots for ES50, deep 
Predicted_ES50_depth_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Depth_Mean, y = predict(deep.ES50.depth, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Depth (m)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_depth_dee_gam.tiff", plot = Predicted_ES50_depth_dee_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_margin_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Margin_Sum, y = predict(deep.ES50.margin, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Continental Margin (km2)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_margin_dee_gam.tiff", plot = Predicted_ES50_margin_dee_gam, width = 6, height = 4, dpi = 600)


Predicted_ES50_current_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = CurVel_Bot_Mean, y = predict(deep.ES50.current, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Current Velocity (m.s-1)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_current_dee_gam.tiff", plot = Predicted_ES50_current_dee_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_humimp_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = HumImp_Mean, y = predict(deep.ES50.humimp, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Human Impact",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_humimp_dee_gam.tiff", plot = Predicted_ES50_humimp_dee_gam, width = 6, height = 4, dpi = 600)

#Predicted_ES50_icecover_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = IceCov_Mean, y = predict(deep.ES50.icecover, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Ice Cover (fraction)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_icecover_dee_gam.tiff", plot = Predicted_ES50_icecover_dee_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_nitrate_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = Nitrate_Bot_Mean, y = predict(deep.ES50.nitrate, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Nitrate (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_nitrate_dee_gam.tiff", plot = Predicted_ES50_nitrate_dee_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_o2_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = O2_Bot_Mean, y = predict(deep.ES50.o2, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "O2 (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_o2_dee_gam.tiff", plot = Predicted_ES50_o2_dee_gam, width = 6, height = 4, dpi = 600)


Predicted_ES50_PrimProd_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = PrimProd_Bot_Mean, y = predict(deep.ES50.PrimProd, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Primary Productivity (mmol . m-3)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_PrimProd_dee_gam.tiff", plot = Predicted_ES50_PrimProd_dee_gam, width = 6, height = 4, dpi = 600)

Predicted_ES50_ThemM_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = ThemM_Bot_mean, y = predict(deep.ES50.ThemM, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Mean (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
ggsave("Predicted_ES50_ThemM_dee_gam.tiff", plot = Predicted_ES50_ThemM_dee_gam, width = 6, height = 4, dpi = 600)

#Predicted_ES50_ThemR_dee_gam <- ggplot(Ecological_Data_Global_Hex_sp, aes(x = ThemR_Bot_mean, y = predict(deep.ES50.ThemR, Ecological_Data_Global_Hex_sp))) +
  geom_smooth(method = "gam", formula = y ~ x, color = "#1a80bb", fill = "#85bede") +  # Add a smooth dark blue line with light blue shadow
  geom_point(size = 3) +  # Add scatter plot points
  theme_bw() +  # Use the black and white theme
  labs(
    x = "Temperature Range (ºC)",  # Shorten the x-axis title
    y = "Predicted Value"  # Shorten the y-axis title
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size =20),  # Increase x-axis text size
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22)  # Increase y-axis title size
  )
#ggsave("Predicted_ES50_ThemR_dee_gam.tiff", plot = Predicted_ES50_ThemR_dee_gam, width = 6, height = 4, dpi = 600)

