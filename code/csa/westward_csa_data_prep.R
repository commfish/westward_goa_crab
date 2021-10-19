
library(tidyverse)
library(lubridate)

# data ----

## survey
### load specimen data
purrr::map(grep("specimenDump", list.files("./data/csa/trawl_data", full.names = T), value = T), read_csv) %>%
  do.call("rbind", .) %>%
  rename_all(tolower) -> spec

### load catch dump data
purrr::map(grep("catchDump", list.files("./data/csa/trawl_data", full.names = T), value = T), read_csv) %>%
  do.call("rbind", .) %>%
  rename_all(tolower) -> catch

### list of standard stations 
std_stations <- read_csv('./data/misc/tanner_standard_stations.csv')


## fishery
### retained catch
fish <- read_csv("./data/csa/retained_catch_est.csv") %>%
            #convert district and section names
            left_join(std_stations %>%
                      distinct(tdist, tsect, district, section))

# compute cpue timeseries per group ----

## extract haul data
catch %>%
  mutate(year = as.numeric(substring(survey, 5, 8))) %>%
  group_by(year, station) %>%
  summarise(distance_km = mean(`distance(km)`, na.rm = T),
            area_swept = distance_km * 0.0122) -> tow_dist

  
## cpue count per station
spec %>%
  # filter for male tanner crab with shell condition 2-3
  filter(race_code %in% c(68560, 68541),
         sex == 1) %>%
  # round nsize to nearest mm
  mutate(cw = round(nsize, 0)) %>%
  # define stages based on size and shell condition
  mutate(group = case_when((cw > 114 & cw <= 139) ~ "pre_r",
                           (cw >= 140 & cw < 165 & shell %in% c(1, 2, 9)) ~ "r",
                           (cw >= 165 | (cw >= 140 & shell %in% 3:4)) ~ "post_r"),
         # extract year
         year = as.numeric(substring(survey, 5, 8))) %>%
  # filter for standard stations
  filter(station %in% std_stations$station,
         group %in% c("pre_r", "r", "post_r")) %>%
  # get counts per station
  group_by(year, station, group, tdist, tsect) %>%
  summarise(cnt = sum(sampfrac, na.rm = T)) %>%
  # join to tow distance and station area
  left_join(tow_dist) %>%
  left_join(std_stations %>% dplyr::select(station, district, section, area_km2), by = "station") %>%
  # compute station abundance
  mutate(abundance_station = cnt / area_swept * area_km2,
         cpue_tmp = cnt / area_swept) %>%
  # compute a variance 
  # compute average cpue per year
  group_by(year, group, district, section) %>%
  # add station weight
  mutate(station_wt = area_km2 / sum(area_km2, na.rm = T)) %>%
  summarise(abundance = sum(abundance_station),
            cpue = weighted.mean((cnt / area_swept), w = station_wt, na.rm = T),
            # weighted std error 
            se_cpue = sqrt(modi::weighted.var((cnt / area_swept), w = station_wt, na.rm = T) / n()),
            ## se_cpue = sqrt(((sum(station_wt*((cnt/area_swept)-cpue)^2))/(((n()-1)/n())*sum(station_wt))) / n()),
            se_abund = sqrt(se_cpue^2 * sum(area_km2)^2)) %>% 
  ungroup -> index

# compute mid dates, and tau ----

## survey mid dates
catch %>%
  mutate(year = as.numeric(substring(survey, 5, 8)),
         haul_date = mdy(haul_date)) %>%
  group_by(year, tdist, tsect) %>%
  summarise(survey_mid_date = mean(haul_date)) %>%
  ungroup -> survey_mid_date

## fishery mid dates
fish %>%
  mutate(open_date = mdy(open_date),
         close_date = mdy(close_date)) %>%
  group_by(year, tdist, tsect) %>%
  mutate(fishery_mid_date = open_date + floor(close_date-open_date)/2) %>%
  left_join(survey_mid_date, by = c("year", "tdist", "tsect")) %>%
  group_by(tdist, tsect) %>%
  nest %>%
  mutate(tau_cs = purrr::map(data, function(x) {as.numeric((yday(x$survey_mid_date) - lag(yday(x$fishery_mid_date), 1)) / 365)})) %>%
  unnest(c(data, tau_cs)) %>%
  dplyr::select(year, tdist, tsect, tau_cs) %>%
  #convert district and section names
  left_join(std_stations %>%
              distinct(tdist, tsect, district, section)) %>%
  ungroup -> tau_cs

## tau s
survey_mid_date %>%
  group_by(tdist, tsect) %>%
  nest %>%
  mutate(tau_s = purrr::map(data, function(x) {as.numeric((x$survey_mid_date - lag(x$survey_mid_date, 1)) / 365)})) %>%
  unnest(c(data, tau_s)) %>%
  dplyr::select(year, tdist, tsect, tau_s) %>%
  #convert district and section names
  left_join(std_stations %>%
              distinct(tdist, tsect, district, section)) %>%
  ungroup -> tau_s


# compute annual phi ----
## extract haul data
catch %>%
  mutate(year = as.numeric(substring(survey, 5, 8))) %>%
  group_by(year, station) %>%
  summarise(distance_km = mean(`distance(km)`, na.rm = T),
            area_swept = distance_km * 0.0122) -> tow_dist


## cpue count per station
spec %>%
  # filter for male tanner crab with shell condition 2-3
  filter(race_code %in% c(68560, 68541),
         sex == 1) %>%
  # round nsize to nearest mm
  mutate(cw = round(nsize, 0)) %>%
  # define stages based on size and shell condition
  mutate(group = case_when((cw > 114 & cw <= 139) ~ "pre_r",
                           (cw >= 140 & cw < 165 & shell %in% c(1, 2, 9)) ~ "r",
                           (cw >= 165 | (cw >= 140 & shell %in% 3:4)) ~ "post_r"),
         # extract year
         year = as.numeric(substring(survey, 5, 8))) %>%
  # filter for standard stations
  filter(station %in% std_stations$station,
         group %in% c("pre_r")) %>%
  # get counts per station per 1mm size bin
  group_by(year, station, group, tdist, tsect, cw) %>%
  summarise(cnt = sum(sampfrac, na.rm = T)) %>%
  # join to tow distance and station area
  left_join(tow_dist) %>%
  left_join(std_stations %>% dplyr::select(station, district, section, area_km2), by = "station") %>%
  # compute station abundance
  mutate(abundance_station = cnt / area_swept * area_km2) %>%
  ungroup %>%
  mutate(theta_i = 1 / (1 + exp(-0.081009 * (cw - (-(-8.843975/0.081009)))))) %>%
  # compute the weighted mean theta per year
  group_by(year, tdist, tsect) %>%
  summarise(phi = 1-weighted.mean(theta_i, w = abundance_station)) %>%
  #convert district and section names
  left_join(std_stations %>%
            distinct(tdist, tsect, district, section)) %>%
  ungroup  -> phi

# extract district, section datasets ----


index %>%
  dplyr::select(year, group, district, section, abundance) %>%
  pivot_wider(names_from = group, values_from = abundance) %>%
  dplyr::select(year, district, section, pre_r, r, post_r) %>%
  # join to se
  left_join(index %>%
              dplyr::select(year, group, district, section, se_abund) %>%
              pivot_wider(names_from = group, values_from = se_abund) %>%
              rename(pre_r_se = pre_r, r_se = r, post_r_se = post_r) %>%
              dplyr::select(year, district, section, pre_r_se, r_se, post_r_se)) %>%
  # join to tau
  left_join(tau_s %>% dplyr::select(-tdist, -tsect), by = c("year", "district", "section")) %>%
  left_join(tau_cs %>% dplyr::select(-tdist, -tsect), by = c("year", "district", "section")) %>%
  # join to phi
  left_join(phi %>% dplyr::select(-tdist, -tsect), by = c("year", "district", "section")) %>%
  # retained catch 
  left_join(fish %>%
              dplyr::select(year, district, section, retained_num), by = c("year", "district", "section")) %>%
  replace_na(list(tau_cs = 1, retained_num = 0)) %>%
  group_by(district, section) %>%
  nest() %>%
  mutate(area = gsub(" |/", "_", paste(district, section, sep = "_")),
         csv = purrr::map2(data, area, function(x, y){write_csv(x, paste0("./data/csa/", y, "_in.csv"))}))
  





# average weights by relevant groups ----
alpha = 0.00027
beta = 3.022134

spec %>%
  mutate(year = as.numeric(substring(survey, 5, 8)),
       # add group
       group = case_when((sex == 1 & nsize > 114 & nsize < 125) ~ "pre_r",
                         (sex == 1 & nsize >= 125 & nsize < 165 & shell == 2) ~ "r",
                         (sex == 1 & (nsize >= 165 | (nsize >= 125 & shell > 2))) ~ "post_r")) %>%
  filter(station %in% std_stations,
         group %in% c("pre_r", "r", "post_r"))