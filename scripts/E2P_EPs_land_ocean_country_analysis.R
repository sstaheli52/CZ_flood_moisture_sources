library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(openxlsx)
library(exactextractr)

#Load the evaporation NetCDF file
et_data <- rast("thesis_r_code/data/E2P_EPs_only.nc")

#Load the countries shapefile
countries <- st_read("thesis_r_code/data/countries_global_OR/CNTR_RG_20M_2016_4326_latlon.shp")

#Ensure CRS consistency between raster and shapefile
crs_et <- crs(et_data)
crs_countries <- st_crs(countries)

if (crs_countries != crs_et) {
  countries <- st_transform(countries, crs_et)
}

#Rasterize the countries shapefile (land = 1, ocean = 0)
land_ocean_raster <- rasterize(countries, et_data, field = 1, background = 0)

#Reclassify the raster: land = 1, ocean = NA
reclass_matrix <- matrix(c(0, NA, 1, 1), ncol = 2, byrow = TRUE)
land_ocean_raster <- classify(land_ocean_raster, reclass_matrix)

#Apply the land mask to the evaporation data
et_land <- mask(et_data, land_ocean_raster)

#Mask the ET data to extract ocean values (ocean = NA, land = 1)
land_ocean_raster <- classify(land_ocean_raster, matrix(c(1, NA, NA, 1), ncol = 2, byrow = TRUE))
et_ocean <- mask(et_data, land_ocean_raster)

#Extract daily evaporation values for each country using exact_extract
country_et_values <- exact_extract(et_land, countries, "sum")

#Convert the results to a data frame
country_et_values_df <- data.frame(NAME_ENGL = countries$NAME_ENGL)

#Add each band (day) as a separate column
for (i in 1:nlyr(et_land)) {
  country_et_values_df[[paste0("Day_", i)]] <- exact_extract(et_land[[i]], countries, "sum")
}

#Reshape the country evaporation values to a long format
country_et_long <- country_et_values_df %>%
  pivot_longer(cols = starts_with("Day_"), 
               names_to = "Day", 
               values_to = "ET_sum") %>%
  mutate(Day = as.integer(gsub("Day_", "", Day)))

#Extract total evaporation values for the ocean
ocean_et_values <- sapply(1:nlyr(et_ocean), function(i) as.numeric(global(et_ocean[[i]], fun = "sum", na.rm = TRUE)))
ocean_et_values_df <- data.frame(Day = 1:31, Ocean_Evaporation = ocean_et_values)

#Combine land and ocean evaporation data for summaries
total_evaporation <- country_et_long %>%
  group_by(Day) %>%
  summarize(Land_Evaporation = sum(ET_sum, na.rm = TRUE), .groups = "drop") %>%
  left_join(ocean_et_values_df, by = "Day") %>%
  mutate(Total_Evaporation = Land_Evaporation + Ocean_Evaporation,
         Land_Percentage = ifelse(Total_Evaporation > 0, 100 * Land_Evaporation / Total_Evaporation, 0),
         Ocean_Percentage = ifelse(Total_Evaporation > 0, 100 * Ocean_Evaporation / Total_Evaporation, 0))

#Add percentage contribution of each country to the total evaporation
country_et_long <- country_et_long %>%
  left_join(total_evaporation %>% select(Day, Total_Evaporation), by = "Day") %>%
  mutate(Percentage_Evaporation = 100 * ET_sum / Total_Evaporation)

#Export both results to Excel
write.xlsx(country_et_long %>% filter(ET_sum > 0), "thesis_r_code/outputs/E2P_EPs_country_evaporation_summary.xlsx")

write.xlsx(total_evaporation, "thesis_r_code/outputs/E2P_EPs_land_ocean_evaporation_summary.xlsx")


