#Required Libraries
library(sf)
library(stars)
library(terra)
library(ggplot2)
library(dplyr)
library(writexl)
library(slider)
library(tidyr)
library(purrr)
library(cowplot)
library(patchwork)

#------------------------------
# 1. Load and Prepare Data
#------------------------------

#File paths
basin_shapefile <- "data/dissolved_eubasins/dissolved_hydrobasins_eu.shp"
precip_file <- "data/precip_June_July_1997.nc"
evap_file <- "data/E2P_EPs_only.nc"

#Load shapefile and NetCDF files
basins <- st_read(basin_shapefile)
precip <- read_stars(precip_file)
evap <- read_stars(evap_file)

#Set CRS
st_crs(precip) <- 4326
st_crs(evap) <- 4326

#Convert to terra raster format
precip_terra <- rast(precip)
evap_terra <- rast(evap)

#Convert basins to terra vector
basins_terra <- vect(basins)

#------------------------------
# 2. Mask Rasters and Extract Basin Values
#------------------------------

#Crop and mask to basin area
precip_masked <- mask(crop(precip_terra, basins_terra), basins_terra)
evap_masked <- mask(crop(evap_terra, basins_terra), basins_terra)

#Extract raster sums per basin
precip_basin_sum <- terra::extract(precip_masked, basins_terra, fun = sum, na.rm = TRUE)
evap_basin_sum <- terra::extract(evap_masked, basins_terra, fun = sum, na.rm = TRUE)

#Count valid pixels per basin
precip_pixel_counts <- terra::extract(!is.na(precip_masked), basins_terra, fun = sum)
colnames(precip_pixel_counts)[2] <- "n_pixels_precip"
precip_pixel_counts$basin_id <- basins_terra$LEGEND

evap_pixel_counts <- terra::extract(!is.na(evap_masked), basins_terra, fun = sum)
colnames(evap_pixel_counts)[2] <- "n_pixels_evap"
evap_pixel_counts$basin_id <- basins_terra$LEGEND

#------------------------------
# 3. Normalize and Reshape Data
#------------------------------

#Prepare data frames
precip_df <- precip_basin_sum %>%
  rename_with(~ gsub("lyr.", "Day_", .x)) %>%
  mutate(basin_id = ID) %>%
  left_join(precip_pixel_counts, by = "basin_id") %>%
  mutate(across(starts_with("Day_"), ~ .x / n_pixels_precip))

evap_df <- evap_basin_sum %>%
  rename_with(~ gsub("lyr.", "Day_", .x)) %>%
  mutate(basin_id = ID) %>%
  left_join(evap_pixel_counts[, c("basin_id", "n_pixels_evap")], by = "basin_id") %>%
  mutate(across(starts_with("Day_"), ~ .x / n_pixels_evap))

#Wide to long
precip_long <- precip_df %>%
  select(-starts_with("lyr."), -ID.x, -ID.y) %>%
  pivot_longer(cols = starts_with("Day_"), names_to = "date", values_to = "precip") %>%
  mutate(date = as.Date("1997-06-01") + as.numeric(gsub("Day_", "", date)) - 1)

evap_long <- evap_df %>%
  pivot_longer(cols = starts_with("Day_"), names_to = "date", values_to = "evap") %>%
  mutate(date = as.Date("1997-07-01") + as.numeric(gsub("Day_", "", date)) - 1)

#------------------------------
# 4. Merge and Calculate Rolling Correlations
#------------------------------

#Merge and fill NAs with 0
merged_data <- full_join(precip_long, evap_long, by = c("basin_id", "date")) %>%
  mutate(
    precip = ifelse(is.na(precip), 0, precip),
    evap = ifelse(is.na(evap), 0, evap)
  )

#Add rolling sums
merged_data <- merged_data %>%
  group_by(basin_id) %>%
  arrange(date) %>%
  mutate(
    precip_15day = slide_dbl(precip, ~ sum(.x, na.rm = TRUE), .before = 14, .complete = TRUE),
    precip_30day = slide_dbl(precip, ~ sum(.x, na.rm = TRUE), .before = 29, .complete = TRUE)
  ) %>%
  ungroup()

#Function for rolling Spearman's correlation
calc_rolling_cor_with_p <- function(df, window_size = 7) {
  df <- df %>%
    arrange(date) %>%
    mutate(
      rho = slide2_dbl(
        .x = precip,
        .y = evap,
        .f = ~ if (sum(complete.cases(.x, .y)) >= 3) {
          cor(.x, .y, method = "spearman", use = "complete.obs")
        } else { NA_real_ },
        .before = window_size - 1,
        .complete = TRUE
      ),
      p_value = slide2_dbl(
        .x = precip,
        .y = evap,
        .f = ~ if (sum(complete.cases(.x, .y)) >= 3) {
          suppressWarnings(cor.test(.x, .y, method = "spearman")$p.value)
        } else { NA_real_ },
        .before = window_size - 1,
        .complete = TRUE
      )
    )
  return(df)
}

#Basin name lookup
id_name_lookup <- data.frame(
  ID = basins_terra$LEGEND,
  basin_name = basins_terra$MAJ_NAME
)

#Run for both 15 and 30-day windows
rolling_multiwindow <- map_dfr(c(15, 30), function(w) {
  merged_data %>%
    group_by(basin_id) %>%
    group_modify(~ calc_rolling_cor_with_p(.x, window_size = w)) %>%
    mutate(window = paste0(w, "-day"))
}) %>%
  ungroup() %>%
  left_join(id_name_lookup, by = c("basin_id" = "ID"))

#------------------------------
# 5. Plotting Function
#------------------------------

plot_basin_correlation <- function(basin_name_input) {
  filtered <- rolling_multiwindow %>%
    filter(basin_name == basin_name_input,
           date >= as.Date("1997-07-01") & date <= as.Date("1997-07-31")) %>%
    mutate(
      significant = case_when(
        is.na(p_value) ~ "NA",
        p_value < 0.05 ~ "Significant",
        TRUE ~ "Not Significant"
      )
    )
  
  plots <- lapply(c("15-day", "30-day"), function(w) {
    df_win <- filter(filtered, window == w)
    
    ggplot(df_win, aes(x = date, y = rho)) +
      geom_line(color = "steelblue", na.rm = TRUE) +
      geom_point(aes(shape = significant, color = significant), size = 2, stroke = 1) +
      geom_point(
        data = filter(df_win, is.na(rho)),
        aes(x = date, y = 0),
        shape = 4,
        color = "lightblue",
        size = 3,
        stroke = 1
      ) +
      scale_shape_manual(
        values = c("Significant" = 16, "Not Significant" = 16, "NA" = 4),
        name = "Significance", drop = FALSE
      ) +
      scale_color_manual(
        values = c("Significant" = "red", "Not Significant" = "black", "NA" = "lightblue"),
        name = "Significance", drop = FALSE
      ) +
      scale_x_date(date_breaks = "3 days", date_labels = "%d") +
      labs(title = w, x = "Date", y = expression("Spearman " * rho)) +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  combined_plot <- cowplot::plot_grid(plots[[1]], plots[[2]], ncol = 2)
  
  ggsave(
    paste0("outputs/Correlation_plots/", gsub(" ", "", basin_name_input), "_timeseries_correlation.pdf"),
    combined_plot, width = 10, height = 4
  )
  
  return(list(combined_plot = combined_plot, individual_plots = plots))
}


#Example use:
plot_basin_correlation("Danube")

#Save the legend (Narva's plot includes all important legend information)
result <- plot_basin_correlation("Narva")
print(result)

legend_only <- get_legend(result$individual_plots[[1]] + theme(legend.position = "right"))
legend_plot <- plot_grid(legend_only)
ggsave("outputs/Correlation_plots/rolling_correlation_legend.pdf", legend_plot)

#------------------------------
# 6. Summary Table for Key Basins
#------------------------------

key_basins <- c("Danube", "Elbe", "Sweden", "Wisla", "Oder", "Dnieper", "Poland Coast",
                "Neman", "Dniester", "Volga", "Daugava", "Narva", "Baltic Sea Coast",
                "Neva", "Scandinavia, North Coast", "Rhine")

filtered_corr <- rolling_multiwindow %>%
  filter(
    basin_name %in% key_basins,
    date >= as.Date("1997-07-01") & date <= as.Date("1997-07-31"),
    !is.na(p_value)
  ) %>%
  mutate(
    significant = p_value < 0.05,
    evap_zero = evap == 0,
    precip_zero = precip == 0,
    zero_involved = evap_zero | precip_zero
  )

summary_table <- filtered_corr %>%
  group_by(basin_name, window) %>%
  summarise(
    mean_rho = round(mean(rho, na.rm = TRUE), 2),
    significant_dates = paste(format(date[significant], "%b %d"), collapse = ", "),
    zero_sig_dates = paste(format(date[significant & zero_involved], "%b %d"), collapse = ", ")
  ) %>%
  ungroup()

#Save summary table
write_xlsx(summary_table, "outputs/rolling_correlation_summary_table.xlsx")
