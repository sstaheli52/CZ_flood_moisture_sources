library(terra)
library(dplyr)
library(ggplot2)
library(sf)
library(grid)
library(ggpubr)
library(writexl)

#------------------------------
# 1. File Paths and Basin Setup
#------------------------------

#Data paths
yearly_files_path <- "data/yearly_means/"
years <- 1980:2020
yearly_precip_files_14 <- paste0(yearly_files_path, "precip_14days_mean_", years, ".nc")
yearly_precip_files_31 <- paste0(yearly_files_path, "precip_31days_mean_", years, ".nc")

#Basin shapefile
hydrobasin_path <- "data/dissolved_eubasins/dissolved_hydrobasins_eu.shp"
hydrobasins <- st_read(hydrobasin_path)

#Selected key basins
selected_basins <- c("Dnieper", "Danube", "Sweden", "Wisla", "Baltic Sea Coast",
                     "Daugava", "Narva", "Neva", "Neman", "Oder", "Dniester",
                     "Volga", "Scandinavia, North Coast", "Poland Coast", "Elbe")

hydrobasins_filtered <- hydrobasins %>%
  filter(MAJ_NAME %in% selected_basins)

#Convert shapefile to terra format
hydrobasins_vect_filtered <- vect(hydrobasins_filtered)

#------------------------------
# 2. Exceedance Calculation Function
#------------------------------

#Calculates exceedance probability of 1997 based on historical ranks
exceedance_probability <- function(yearly_precip_files) {
  climatology <- rast(yearly_precip_files)
  
  #Rank values in descending order, 1997 is the 18th year
  rank_1997 <- app(climatology, fun = function(x) rank(-x, ties.method = "min"))[[18]]  
  
  #Calculate exceedance %: (rank - 1) / total years
  exceedance_1997 <- (rank_1997 - 1) / 41 * 100  
  
  #Mask to basins
  exceedance_clipped <- mask(exceedance_1997, hydrobasins_vect_filtered)
  
  #Convert to dataframe
  exceedance_df <- as.data.frame(exceedance_clipped, xy = TRUE)
  colnames(exceedance_df) <- c("longitude", "latitude", "exceedance")
  exceedance_df <- na.omit(exceedance_df)
  
  #Bin into exceedance categories
  exceedance_df$category <- cut(exceedance_df$exceedance,
                                breaks = c(0, 2, 5, 10, 20, 30, 40, 50, 100),
                                labels = c("<2%", "2-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", ">50%"),
                                include.lowest = TRUE)
  return(exceedance_df)
}

#Run for 14 and 31-days
exceedance_df_14 <- exceedance_probability(yearly_precip_files_14)
exceedance_df_31 <- exceedance_probability(yearly_precip_files_31)

#------------------------------
# 3. Areal Fractions (Whole Domain)
#------------------------------

#Compute fraction of area falling into each exceedance category
calc_fractions <- function(df) {
  df %>%
    count(category) %>%
    mutate(fraction = n / sum(n)) %>%
    mutate(category = factor(category,
                             levels = c("<2%", "2-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", ">50%")))
}

category_arealFractions_14 <- calc_fractions(exceedance_df_14)
category_arealFractions_31 <- calc_fractions(exceedance_df_31)

#Save to Excel
write_xlsx(
  list(
    "14-day Global Fractions" = category_arealFractions_14,
    "31-day Global Fractions" = category_arealFractions_31
  ),
  path = "outputs/global_exceedance_fractions.xlsx"
)

#------------------------------
# 4. Exceedance Map Plotting
#------------------------------

exceedance_colors <- c("<2%" = "#6a0301", "2-5%" = "#E60000", "5-10%" = "#fda52e", "10-20%" = "#45a05c",
                       "20-30%" = "#126a4a", "30-40%" = "#f9d869", "40-50%" = "#4b4a4a", ">50%" = "#D3D3D3")  

#Generate map
exceedance_map <- function(df) {
  ggplot() +
    geom_tile(data = df, aes(x = longitude, y = latitude, fill = category)) +
    scale_fill_manual(values = exceedance_colors, name = "Exceedance\nProbability (%)") +
    geom_sf(data = hydrobasins, fill = NA, color = "black", linewidth = 1) +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none") +
    labs(x = NULL, y = NULL)
}

#Run function and save maps
plot_14days <- exceedance_map(exceedance_df_14)
plot_31days <- exceedance_map(exceedance_df_31)

ggsave("outputs/Exceedance_plots/exceedance_14days_map.pdf", plot_14days, width = 8, height = 6)
ggsave("outputs/Exceedance_plots/exceedance_31days_map.pdf", plot_31days, width = 8, height = 6)

#------------------------------
# 5. Legend Export
#------------------------------

#Create dummy plot just for legend
legend_plot <- ggplot(data.frame(x = 1, y = 1, category = factor(names(exceedance_colors), levels = names(exceedance_colors))),
                      aes(x = x, y = y, fill = category)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = exceedance_colors, name = "Exceedance \nProbability (%)", guide = guide_legend(nrow = 1)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.width = unit(.5, "cm"),
        legend.key.height = unit(.5, "cm"),
        legend.spacing.x = unit(0.2, "cm"))

legend_horizontal <- get_legend(legend_plot)
ggsave("outputs/Exceedance_plots/legend_exceedance_horizontal.pdf", plot = as_ggplot(legend_horizontal), width = 8, height = 1)

#------------------------------
# 6. Areal Fractions by Basin
#------------------------------

#Match shapefile with raster resolution
hydrobasins_rast <- rasterize(vect(hydrobasins_filtered), rast(yearly_precip_files_14[1]), field = "MAJ_NAME")

#Get basin name for each raster cell in exceedance dataframe
exceedance_df_14$MAJ_NAME <- terra::extract(hydrobasins_rast, exceedance_df_14[, c("longitude", "latitude")], ID = FALSE)
exceedance_df_31$MAJ_NAME <- terra::extract(hydrobasins_rast, exceedance_df_31[, c("longitude", "latitude")], ID = FALSE)

#Function to calculate exceedance category fractions per basin
calc_basin_fractions <- function(df) {
  df %>%
    filter(!is.na(MAJ_NAME)) %>%
    group_by(MAJ_NAME, category) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(MAJ_NAME) %>%
    mutate(fraction = count / sum(count)) %>%
    mutate(MAJ_NAME = as.character(MAJ_NAME$MAJ_NAME))
}

areal_fractions_14 <- calc_basin_fractions(exceedance_df_14)
areal_fractions_31 <- calc_basin_fractions(exceedance_df_31)

#Save as Excel
write_xlsx(
  list("Areal Fractions 14 days" = areal_fractions_14, 
       "Areal Fractions 31 days" = areal_fractions_31),
  path = "outputs/Exceedance_plots/areal_fractions_basins_exceedance_table.xlsx"
)
