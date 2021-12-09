# notes ----

## estimation of goa trawl survey data with missing stations
## tyler jackson


# library ----

library(tidyverse)
library(lubridate)
library(patchwork)

# custom function to correct timeseries stations
f_timeseries_station_correct <- function(data, select){
  
  # data for function
  std_stations <- read_csv("data/misc/tanner_standard_stations.csv")
  
  # change data object name to avoiding masking
  dat = data
  
  # minor function needed within function
  ## pull mean from reference data 
  f_pmfd <- function(input, data, select) {
    # change names to avoid masking
    input %>%
      rename(stn = station, 
             yr = year) -> input
    
    # fill in data for stations two in a given year
    if(input$stn %in% (data %>% filter(year == input$yr) %>% pull(station))) {
      data %>%
        filter(year == input$yr, station == input$stn) %>%
        pull(select) %>%
        .[1] -> out
    }
    # fill in data for stations that need avg from references in given year
    if(!(input$stn %in% (data %>% filter(year == input$yr) %>% pull(station)))) {
      data %>%
        filter(year == input$yr, station %in% c(input$stn, input$reference_station, input$reference_station2, input$reference_station3)) %>%
        pull(select) %>%
        mean -> res
      out <- ifelse(is.nan(res), NA, res)
      
    }
    
    # return 
    return(out)
  }
  
  expand_grid(year = min(data$year):max(data$year), std_stations) %>%
    # remove the text stations
    dplyr::select(-district, -section, -subsection) %>%
    # create dummy column year/station
    mutate(yrstn = paste(year, station, sep = "_")) %>%
    group_by(yrstn) %>%
    nest %>%
    mutate(select = purrr::map_dbl(data, f_pmfd, data = dat, select = select)) %>%
    unnest() %>%
    ungroup() %>%
    dplyr::select(year, station, select) -> out
  
  names(out)[3] <- select
  
  return(out)
  
}

# ggplot options
theme_set(FNGr::theme_sleek())
yr_axis <- FNGr::tickr(tibble(year = 1988:2021), year, 5)

# data ----

## standard stations
std_stations <- read_csv("data/misc/tanner_standard_stations.csv")

## hauls
hauls_raw <- do.call("rbind", lapply(list.files(path = "data/trawl_data", pattern = "haul", full.names = T),
                                 read_csv))

## specimen data
spec_raw <- do.call("rbind", lapply(list.files(path = "data/trawl_data", pattern = "specimen", full.names = T),
                                 read_csv))

## cpue data estimated with the rolling block
abundance_rb <- read_csv("./data/estimation_comparison/rolling_block_abundance_est.csv")


# data merging ----

hauls_raw %>%
  rename_all(tolower) %>%
  # remove bad tows
  #filter(performance <= 4) %>%
  # add year
  mutate(year = year(mdy(haul_date))) %>%
  # trim data
  dplyr::select(year, tdist, tsect, tsubsect, haul, station, km_towed) %>%
  # fill in missing tow info with mode
  mutate(km_towed = ifelse(is.na(km_towed), 1.852, km_towed)) -> hauls


spec_raw %>%
  rename_all(tolower) %>%
  # filter for tanner crab, remove bad hauls
  filter(race_code == 68560) %>%
  # add year
  # extract haul number from tow
  mutate(year = year(mdy(haul_date)),
         haul = as.numeric(substring(tow, 6, 8))) %>%
  dplyr::select(year, station, haul, tdist, tsect, tsub_sect, distance, sex, nsize, legal, maturity, sampfrac) -> spec

# total area per section ----
std_stations %>%
  group_by(tdist, tsect, tsub_sect) %>%
  summarise(subsect_area_nmi2 = sum(area_nmi2)) %>%
  group_by(tdist, tsect) %>%
  mutate(sect_area_nmi2 = sum(subsect_area_nmi2)) -> section_area


# compute timeseries of mature male abundance ----

## compute cpue by station
spec %>%
  # filter for mature males
  filter(sex == 1, nsize > 114) %>% 
  # sum crab caught per haul
  group_by(year, station, haul, distance) %>%
  summarise(catch = round(sum(sampfrac, na.rm = T))) %>%
  # average catch and distance towed (if there was more than one tow per station)
  group_by(year, station) %>%
  summarise(catch = mean(catch), 
            distance = mean(distance)) %>%
  ungroup() %>%
  # join to haul data and fill in zero catches
  # unit conversions: 1.852 km to 1 nmi
  right_join(hauls, by = c("year", "station")) %>%
  replace_na(list(catch = 0)) %>%
  mutate(distance = ifelse(is.na(distance), km_towed / 1.852, distance)) %>%
  # compute cpue, net width = 40 ft
  # unit conversions: 6076 ft to 1 nmi
  mutate(cpue = catch / (distance * (40 / 6076))) %>%
  f_timeseries_station_correct(., select = "cpue") %>%
  left_join(std_stations %>%
              dplyr::select(-8:-10)) -> mm_cpue

## expand to abundance
### kodiak, east aleutian (stratified by section)
mm_cpue %>%
  filter(!(tdist %in% c(2, 6)),
         !is.na(tdist)) %>%
  group_by(year, tdist, tsect) %>%
  summarise(mean_cpue = weighted.mean(cpue, w = area_nmi2, na.rm = T),
            var_cpue = modi::weighted.var(cpue, w = area_nmi2, na.rm = T),
            n_tows = sum(!is.na(cpue))) %>%
  ungroup() %>%
  left_join(section_area %>%
              group_by(tdist, tsect) %>%
              summarise(sect_area_nmi2 = mean(sect_area_nmi2)), 
            by = c("tdist", "tsect")) %>%
  mutate(abundance = mean_cpue * sect_area_nmi2,
         se = sqrt((var_cpue * sect_area_nmi2^2) / n_tows)) -> mm_abundance

### south peninsula and chignik (stratified by subsection)
mm_cpue %>%
  filter(tdist %in% c(2, 6)) %>%
  group_by(year, tdist, tsect, tsub_sect) %>%
  summarise(mean_cpue = weighted.mean(cpue, w = area_nmi2, na.rm = T),
            var_cpue = modi::weighted.var(cpue, w = area_nmi2, na.rm = T),
            n_tows = sum(!is.na(cpue))) %>%
  ungroup() %>%
  left_join(section_area, 
            by = c("tdist", "tsect", "tsub_sect")) %>%
  mutate(abundance = mean_cpue * subsect_area_nmi2,
         se = sqrt((var_cpue * subsect_area_nmi2^2) / n_tows)) -> mm_abundance_spen

### combined dataset for all districts 
mm_abundance <- bind_rows(mm_abundance, mm_abundance_spen)

 
# plots of mature male abundance ----

## create plot data
abundance_rb %>%
  dplyr::select(-mfa) %>%
  left_join(mm_abundance, by = c("year", "tdist", "tsect", "tsub_sect")) %>%
  # rename fields
  rename(`Rolling Block` = mma,
         `Design Based` = abundance) %>%
  # rearrange data
  dplyr::select(year, tdist, tsect, tsub_sect, `Rolling Block`,`Design Based`) %>%
  pivot_longer(5:6, names_to = "method", values_to = "abundance") %>%
  # add text section and subsection names
  ## do join in two steps to accomodate kodiak, and east aleut which are missing subsection
  left_join(std_stations %>%
              distinct(tdist, tsect, section), by = c("tdist", "tsect")) %>%
  left_join(std_stations %>%
              distinct(tdist, tsect, tsub_sect, subsection), by = c("tdist", "tsect", "tsub_sect")) -> mma_plot_data

## kodiak district
mma_plot_data %>%
  filter(tdist == 1) %>%
  ggplot()+
  geom_point(aes(x = year, y = abundance / 1e6, color = method), alpha = 0.6)+
  geom_line(aes(x = year, y = abundance / 1e6, linetype = method, color = method), alpha = 0.6)+
  facet_grid(rows = vars(section), scales = "free")+
  scale_y_continuous(labels = scales::comma)+
  scale_x_continuous(breaks = yr_axis$breaks, labels = yr_axis$labels)+
  labs(x = NULL, y = "Mature Male Abundance (millions)", color = NULL, linetype = NULL)+
  theme(legend.position = c(0.13, 0.95)) -> mma_kod
ggsave("./figures/estimation/mature_males_kodiak_rb_design.png", plot = mma_kod, height = 6, width = 6, units = "in")

## south peninsula district
### subsection
mma_plot_data %>%
  filter(tdist == 2) %>%
  ggplot()+
  geom_point(aes(x = year, y = abundance / 1e3, color = method), alpha = 0.6)+
  geom_line(aes(x = year, y = abundance / 1e3, linetype = method, color = method), alpha = 0.6)+
  facet_grid(rows = vars(subsection), scales = "free")+
  scale_y_continuous(labels = scales::comma)+
  scale_x_continuous(breaks = yr_axis$breaks, labels = yr_axis$labels)+
  labs(x = NULL, y = "Mature Male Abundance (thousands)", color = NULL, linetype = NULL)+
  theme(legend.position = c(0.13, 0.95)) -> mma_spen_subsect
ggsave("./figures/estimation/mature_males_southpensubsect_rb_design.png", plot = mma_spen_subsect, height = 7, width = 6, units = "in")
### section
mma_plot_data %>%
  filter(tdist == 2) %>%
  group_by(year, tdist, tsect, section, method) %>%
  summarise(abundance = sum(abundance, na.rm = T)) %>%
  ggplot()+
  geom_point(aes(x = year, y = abundance / 1e3, color = method), alpha = 0.6)+
  geom_line(aes(x = year, y = abundance / 1e3, linetype = method, color = method), alpha = 0.6)+
  facet_grid(rows = vars(section), scales = "free")+
  scale_y_continuous(labels = scales::comma)+
  scale_x_continuous(breaks = yr_axis$breaks, labels = yr_axis$labels)+
  labs(x = NULL, y = "Mature Male Abundance (thousands)", color = NULL, linetype = NULL)+
  theme(legend.position = c(0.13, 0.95)) -> mma_spen_section
ggsave("./figures/estimation/mature_males_southpensect_rb_design.png", plot = mma_spen_section, height = 5, width = 6, units = "in")

## eastern aleutian district
mma_plot_data %>%
    filter(tdist == 3) %>%
    ggplot()+
    geom_point(aes(x = year, y = abundance / 1e3, color = method), alpha = 0.6)+
    geom_line(aes(x = year, y = abundance / 1e3, linetype = method, color = method), alpha = 0.6)+
    facet_grid(rows = vars(section), scales = "free")+
    scale_y_continuous(labels = scales::comma)+
    scale_x_continuous(breaks = yr_axis$breaks, labels = yr_axis$labels)+
    labs(x = NULL, y = "Mature Male Abundance (thousands)", color = NULL, linetype = NULL)+
    theme(legend.position = c(0.13, 0.95)) -> mma_ealeut
ggsave("./figures/estimation/mature_males_ealeut_rb_design.png", plot = mma_ealeut, height = 6, width = 6, units = "in")

## chignik district
mma_plot_data %>%
    filter(tdist == 6) %>%
    ggplot()+
    geom_point(aes(x = year, y = abundance / 1e3, color = method), alpha = 0.6)+
    geom_line(aes(x = year, y = abundance / 1e3, linetype = method, color = method), alpha = 0.6)+
    facet_grid(rows = vars(subsection), scales = "free")+
    scale_y_continuous(labels = scales::comma)+
    scale_x_continuous(breaks = yr_axis$breaks, labels = yr_axis$labels)+
    labs(x = NULL, y = "Mature Male Abundance (thousands)", color = NULL, linetype = NULL)+
    theme(legend.position = c(0.13, 0.95)) -> mma_chignik
ggsave("./figures/estimation/mature_males_chignik_rb_design.png", plot = mma_chignik, height = 6, width = 6, units = "in")

## histogram of proportion
mma_plot_data %>%
  pivot_wider(names_from = method, values_from = abundance) %>%
  mutate(prop_diff = `Design Based` / `Rolling Block`) %>%
  ggplot()+
  geom_histogram(aes(x = prop_diff), binwidth = 0.005, fill = "black", color = NA)+
  labs(x = "Design Based / Rolling Block \n (Males)", y = "Count") -> mm_hist


  

# compute timeseries of mature female abundance ----

## compute cpue by station
spec %>%
  # filter for mature males
  filter(sex == 2, maturity == 2) %>% 
  # sum crab caught per haul
  group_by(year, station, haul, distance) %>%
  summarise(catch = round(sum(sampfrac, na.rm = T))) %>%
  # average catch and distance towed (if there was more than one tow per station)
  group_by(year, station) %>%
  summarise(catch = mean(catch), 
            distance = mean(distance)) %>%
  ungroup() %>%
  # join to haul data and fill in zero catches
  # unit conversions: 1.852 km to 1 nmi
  right_join(hauls, by = c("year", "station")) %>%
  replace_na(list(catch = 0)) %>%
  mutate(distance = ifelse(is.na(distance), km_towed / 1.852, distance)) %>%
  # compute cpue, net width = 40 ft
  # unit conversions: 6076 ft to 1 nmi
  mutate(cpue = catch / (distance * (40 / 6076))) %>%
  f_timeseries_station_correct(., select = "cpue") %>%
  left_join(std_stations %>%
              dplyr::select(-8:-10)) -> mf_cpue

## expand to abundance
### kodiak, east aleutian (stratified by section)
mf_cpue %>%
  filter(!(tdist %in% c(2, 6)),
         !is.na(tdist)) %>%
  group_by(year, tdist, tsect) %>%
  summarise(mean_cpue = weighted.mean(cpue, w = area_nmi2, na.rm = T),
            var_cpue = modi::weighted.var(cpue, w = area_nmi2, na.rm = T),
            n_tows = sum(!is.na(cpue))) %>%
  ungroup() %>%
  left_join(section_area %>%
              group_by(tdist, tsect) %>%
              summarise(sect_area_nmi2 = mean(sect_area_nmi2)), 
            by = c("tdist", "tsect")) %>%
  mutate(abundance = mean_cpue * sect_area_nmi2,
         se = sqrt((var_cpue * sect_area_nmi2^2) / n_tows)) -> mf_abundance

### south peninsula and chignik (stratified by subsection)
mf_cpue %>%
  filter(tdist %in% c(2, 6)) %>%
  group_by(year, tdist, tsect, tsub_sect) %>%
  summarise(mean_cpue = weighted.mean(cpue, w = area_nmi2, na.rm = T),
            var_cpue = modi::weighted.var(cpue, w = area_nmi2, na.rm = T),
            n_tows = sum(!is.na(cpue))) %>%
  ungroup() %>%
  left_join(section_area, 
            by = c("tdist", "tsect", "tsub_sect")) %>%
  mutate(abundance = mean_cpue * subsect_area_nmi2,
         se = sqrt((var_cpue * subsect_area_nmi2^2) / n_tows)) -> mf_abundance_spen

### combined dataset for all districts 
mf_abundance <- bind_rows(mf_abundance, mf_abundance_spen)

# plots of mature female abundance ----

## create plot data
abundance_rb %>%
  dplyr::select(-mma) %>%
  left_join(mf_abundance, by = c("year", "tdist", "tsect", "tsub_sect")) %>%
  # rename fields
  rename(`Rolling Block` = mfa,
         `Design Based` = abundance) %>%
  # rearrange data
  dplyr::select(year, tdist, tsect, tsub_sect, `Rolling Block`,`Design Based`) %>%
  pivot_longer(5:6, names_to = "method", values_to = "abundance") %>%
  # add text section and subsection names
  ## do join in two steps to accomodate kodiak, and east aleut which are missing subsection
  left_join(std_stations %>%
              distinct(tdist, tsect, section), by = c("tdist", "tsect")) %>%
  left_join(std_stations %>%
              distinct(tdist, tsect, tsub_sect, subsection), by = c("tdist", "tsect", "tsub_sect")) -> mfa_plot_data

## kodiak district
mfa_plot_data %>%
  filter(tdist == 1) %>%
  ggplot()+
  geom_point(aes(x = year, y = abundance / 1e6, color = method), alpha = 0.6)+
  geom_line(aes(x = year, y = abundance / 1e6, linetype = method, color = method), alpha = 0.6)+
  facet_grid(rows = vars(section), scales = "free")+
  scale_y_continuous(labels = scales::comma)+
  scale_x_continuous(breaks = yr_axis$breaks, labels = yr_axis$labels)+
  labs(x = NULL, y = "Mature Female Abundance (millions)", color = NULL, linetype = NULL)+
  theme(legend.position = c(0.13, 0.95)) -> mfa_kod
ggsave("./figures/estimation/mature_females_kodiak_rb_design.png", plot = mfa_kod, height = 6, width = 6, units = "in")

## south peninsula district
### subsection
mfa_plot_data %>%
  filter(tdist == 2) %>%
  ggplot()+
  geom_point(aes(x = year, y = abundance / 1e3, color = method), alpha = 0.6)+
  geom_line(aes(x = year, y = abundance / 1e3, linetype = method, color = method), alpha = 0.6)+
  facet_grid(rows = vars(subsection), scales = "free")+
  scale_y_continuous(labels = scales::comma)+
  scale_x_continuous(breaks = yr_axis$breaks, labels = yr_axis$labels)+
  labs(x = NULL, y = "Mature Female Abundance (thousands)", color = NULL, linetype = NULL)+
  theme(legend.position = c(0.13, 0.95)) -> mfa_spen_subsect
ggsave("./figures/estimation/mature_females_southpensubsect_rb_design.png", plot = mfa_spen_subsect, height = 7, width = 6, units = "in")
### section
mfa_plot_data %>%
  filter(tdist == 2) %>%
  group_by(year, tdist, tsect, section, method) %>%
  summarise(abundance = sum(abundance, na.rm = T)) %>%
  ggplot()+
  geom_point(aes(x = year, y = abundance / 1e3, color = method), alpha = 0.6)+
  geom_line(aes(x = year, y = abundance / 1e3, linetype = method, color = method), alpha = 0.6)+
  facet_grid(rows = vars(section), scales = "free")+
  scale_y_continuous(labels = scales::comma)+
  scale_x_continuous(breaks = yr_axis$breaks, labels = yr_axis$labels)+
  labs(x = NULL, y = "Mature Female Abundance (thousands)", color = NULL, linetype = NULL)+
  theme(legend.position = c(0.13, 0.95)) -> mfa_spen_section
ggsave("./figures/estimation/mature_females_southpensect_rb_design.png", plot = mfa_spen_section, height = 5, width = 6, units = "in")

## eastern aleutian district
mfa_plot_data %>%
  filter(tdist == 3) %>%
  ggplot()+
  geom_point(aes(x = year, y = abundance / 1e3, color = method), alpha = 0.6)+
  geom_line(aes(x = year, y = abundance / 1e3, linetype = method, color = method), alpha = 0.6)+
  facet_grid(rows = vars(section), scales = "free")+
  scale_y_continuous(labels = scales::comma)+
  scale_x_continuous(breaks = yr_axis$breaks, labels = yr_axis$labels)+
  labs(x = NULL, y = "Mature Female Abundance (thousands)", color = NULL, linetype = NULL)+
  theme(legend.position = c(0.13, 0.95)) -> mfa_ealeut
ggsave("./figures/estimation/mature_females_ealeut_rb_design.png", plot = mfa_ealeut, height = 6, width = 6, units = "in")

## chignik district
mfa_plot_data %>%
  filter(tdist == 6) %>%
  ggplot()+
  geom_point(aes(x = year, y = abundance / 1e3, color = method), alpha = 0.6)+
  geom_line(aes(x = year, y = abundance / 1e3, linetype = method, color = method), alpha = 0.6)+
  facet_grid(rows = vars(subsection), scales = "free")+
  scale_y_continuous(labels = scales::comma)+
  scale_x_continuous(breaks = yr_axis$breaks, labels = yr_axis$labels)+
  labs(x = NULL, y = "Mature Female Abundance (thousands)", color = NULL, linetype = NULL)+
  theme(legend.position = c(0.13, 0.95)) -> mfa_chignik
ggsave("./figures/estimation/mature_females_chignik_rb_design.png", plot = mfa_chignik, height = 6, width = 6, units = "in")

## histogram of proportion
mfa_plot_data %>%
  pivot_wider(names_from = method, values_from = abundance) %>%
  mutate(prop_diff = `Design Based` / `Rolling Block`) %>%
  ggplot()+
  geom_histogram(aes(x = prop_diff), binwidth = 0.005, fill = "black", color = NA)+
  labs(x = "Design Based / Rolling Block \n (Females)", y = "Count") -> mf_hist
ggsave("./figures/estimation/prop_diff_histograms.png", plot = mm_hist+mf_hist, height = 3, width = 6, units = "in")


