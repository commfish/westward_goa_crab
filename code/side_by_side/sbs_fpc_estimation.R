# notes ----
## resolution modification side by side analysis
## author: Tyler Jackson, continuation from W. Gaueman's analysis
## last updated 4/16/2021


# load ----

library(tidyverse)
library(magrittr)
library(fishmethods)
library(FNGr) # only for ggplot theme, from W. Williams GitHub
library(scales)

# function to rename data fields
f_rename <- function(data) {
  data %>%
    transmute(pair = as.character(rep(1:(dim(.)[1]/2), each = 2)),
              tow = tow,
              vessel = ifelse(vessel_id == 30, "res", "sol"),
              juv_fem = juv_fem,
              adult_fem = adult_fem,
              tot_fem = tot_fem,
              sl70 = `sl_<70`, 
              sl70_91 = `sl_70-91`, 
              sl92_114 = `sl_92-114`, 
              submature = `sl_>114`,
              recruit = recruit, 
              legal = tot_legal, 
              mature = tot_mature, 
              tot_male = tot_male,
              arrowtooth = arrowtooth,
              flathead = flathead, 
              pcod = pacific_cod, 
              pollock = pollock)
}

## global options
### custom color/fill pallete for plotting
cb_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

### set theme (from FNGr)
theme_set(theme_sleek() + theme(legend.position = "bottom"))

# data ----

## cpue data
### pre-sponson
read_csv("./data/side_by_side/sbs2015-17.csv") %>%
  f_rename -> pre
### post-sponson
read_csv("./data/side_by_side/sbs2019.csv") %>%
  f_rename -> post

## numbers of crab per bin
read_csv("./data/side_by_side/SBS-crabdumpallyrs.csv") %>%
  rename_all(tolower) %>%
  ## extract year and vessel
  mutate(vessel = tolower(substring(survey, 4, 6)), 
         year = as.numeric(substring(survey, 7, 11))) %>%
  ## compute area swept
  mutate(km_towed = distance * 1.852) %>%
  ## create size bin, combine bins smaller than 60, greater than 120
  filter(!is.na(isize)) %>%
  
  ## in 10 mm bins
  # mutate(sizebin = ceiling(isize / 10) * 10,
  #        sizebin = ifelse(sizebin < 40, paste0("crab", 40),
  #                  ifelse(sizebin > 130, "crab131_plus", paste0("crab", sizebin)))) %>%
  
  ## bin in 3 size bins
  mutate(sizebin = case_when(isize <= 50 ~ "crab_leq50",
                             isize > 50 & isize <= 100 ~ "crab51_100",
                             isize > 100 ~ "crab_geq101")) %>%

  
  ## get count per size in 1 mm bins
  group_by(vessel, tow, km_towed, sizebin) %>%
  summarise(cpue = round(sum(sampfrac, na.rm = T) / mean(km_towed)),
            .groups = "drop") %>%
  pivot_wider(names_from = sizebin, values_from = cpue) %>%
  dplyr::select(-km_towed) %>%
  replace(is.na(.), values = 0) -> spec

# set spp group order for tables and plots
spp_order <- c("juv_fem", "adult_fem", "tot_fem", "sl70", "sl70_91", "sl92_114", "recruit", 
               "legal", "submature", "mature", "tot_male", names(spec)[3:ncol(spec)], 
               "arrowtooth", "flathead", "pcod", "pollock")

## size comp data
read_csv("./data/side_by_side/SBS-crabdumpallyrs.csv") %>%
  rename_all(tolower) %>%
  ## extract year and vessel
  mutate(vessel = tolower(substring(survey, 4, 6)),
         vessel = ifelse(vessel == "res", "R/V Resolution", "R/V Solstice"),
         year = as.numeric(substring(survey, 7, 11)),
         sponson = ifelse(year < 2019, "Pre", "Post"),
         sponson = factor(sponson, levels = c("Pre", "Post"))) %>%
  group_by(year, vessel, sponson, sex, isize, maturity) %>%
  summarise(ncrab = round(sum(sampfrac)), .groups = "drop") -> size_comp

# fpc functions ----

## ratio estimator
f_ratio_fpc <- function(data, ves1, ves2, vessel_name = "vessel", cpue_name = "cpue") {
  # ensure vessel and cpue columns are properly named
  data %>%
    rename(vessel = as.name(vessel_name),
           cpue = as.name(cpue_name)) %>%
  # pivot data to wide format (i.e., in pairs)
    pivot_wider(names_from = vessel, values_from = cpue) -> data
  
  # extract cpue vectors from data
  y = data %>% pull(as.name(ves1)) # standard vessel
  x = data %>% pull(as.name(ves2))
  
  # get n
  if(length(x) != length(y)) {stop("x and y are different lengths")}
  n = length(x)
  # get means
  x_bar = mean(x)
  y_bar = mean(y)
  # compute ratio
  r_hat = y_bar / x_bar
  # variance in ratio
  var_r_hat = (1 / (n * x_bar^2)) * (sum((y - r_hat * x)^2) / (n - 1))
  
  tibble(n = n,
         # cv of mean cpue (se / mean)
         cv_x = sd(x) / sqrt(n) / mean(x), 
         cv_y = sd(y) / sqrt(n) / mean(y),
         fpc = r_hat,
         l95 = fpc + qnorm(0.025) * sqrt(var_r_hat),
         u95 = fpc + qnorm(0.975) * sqrt(var_r_hat)) -> out
  names(out)[3] = paste0("cv_", ves1, "_cpue")
  names(out)[2] = paste0("cv_", ves2, "_cpue")
  
  return(out)
  }

## random block anova
f_rand_block <- function(data, vessel_name = "vessel", cpue_name = "cpue") {
  
  # ensure vessel and cpue columns are properly named
  data %>%
    rename(vessel = as.name(vessel_name),
           cpue = as.name(cpue_name)) -> data
  
  # use sum contrasts for linear model
  options(contrasts = rep("contr.sum", 2))
  
  mutate(data, pair = as.character(pair)) -> data
  
  # fit linear model and extract coefficients
  fit <- lm(log(cpue + 1) ~ pair + vessel, data = data)
  coefs <- summary(fit)[["coefficients"]]
  # number of hauls
  n = length(unique(data$pair)) 
  # grand mean
  mu = coefs[1, 1] 
  # resolution treatment effect
  nu = coefs[nrow(coefs), 1] 
  nu_se = coefs[nrow(coefs), 2] 
  p_val = coefs[nrow(coefs), 4] 
  # fpc
  fpc = exp(2 * nu * (1 + (0.5 * nu_se^2))) 
  l95 = exp((2 * nu) - (1.96 * 2 * nu_se))
  u95 = exp((2 * nu) + (1.96 * 2 * nu_se))
  
  # output a tibble
  tibble(lmod = list(fit),
         n = n,
         mu = mu,
         vessel_effect = nu,
         se = nu_se,
         p_val = p_val,
         fpc = fpc,
         l95 = l95,
         u95 = u95)
}

## kappenman (1992)
### wrapper for fishmethods::fpc for paired data
f_kapp <- function(data, ves1, ves2, vessel_name = "vessel", cpue_name = "cpue") {
  
  # ensure vessel and cpue columns are properly named
  data %>%
    rename(vessel = as.name(vessel_name),
           cpue = as.name(cpue_name)) -> data
  
  # pivot data to wide format (i.e., in pairs)
  data %>%
    pivot_wider(names_from = vessel, values_from = cpue) -> data
  # extract cpue vectors from data
  cpue1 = data %>% pull(as.name(ves1)) # standard vessel
  cpue2 = data %>% pull(as.name(ves2)) 
  
  # call fishmethods::fpc
  #kapp = try(fishmethods::fpc(cpue1, cpue2, method = 4, kapp_zeros = "ind", boot_type = "unpaired"), silent = T)
  kapp = try(fishmethods::fpc(cpue1, cpue2, method = 4), silent = T)
  
  # organize results
  if(class(kapp) == "try-error"){
    tibble(n = NA,
           cpue1 = NA,
           cpue2 = NA,
           fpc = NA,
           boot_se = NA,
           boot_l95 = NA,
           boot_u95 = NA) -> out
  } else{
    kapp %>%
      as_tibble() %>%
      rename(n = n1,
             cpue1 = `mean cpue1`,
             cpue2 = `mean cpue2`,
             fpc = FPC,
             boot_se = Boot_SE,
             boot_l95 = `Boot_95%_LCI`,
             boot_u95 = `Boot_95%_UCI`) %>%
      dplyr::select(n, cpue1, cpue2, fpc, boot_se, boot_l95, boot_u95) -> out
  }
    
  # fix names
  names(out)[2] = paste0(ves1, "_cpue")
  names(out)[3] = paste0(ves2, "_cpue")
  
  # output
  return(out)
  
}

### mse simulation (Munro 1998; von Szalay and Brown 2001)
f_kapp_mse_sim <- function(data, n, method){
  
  # pull cpue from data
  data %>%
    pull(cpue) -> cpue  
    
  # compute parameters of lognormal dist
  mu = log(mean(cpue)) - 0.5*log(1+(sd(cpue)/mean(cpue))^2)
  var = log(1+(sd(cpue)/mean(cpue))^2)
  
  # check pdf
  #hist(cpue, freq = F)
  #lines(density(rlnorm(1000, meanlog = mu, sdlog = sqrt(var))), col = 2)
  
  cpue_u = NULL
  cpue_c = NULL
  mse_u = NULL
  mse_c = NULL
  fpd = seq(0.1, 2, 0.1)
  for(i in 1:length(fpd)){
    for(j in 1:200){
      # simulate cpue standard vessel
      res_sim = rlnorm(n, meanlog = mu, sdlog = sqrt(var))
      # simulate cpue non-standard vessel
      sol_sim = rlnorm(n, meanlog = mu, sdlog = sqrt(var)) * fpd[i]
      # estimate FPC
      fpc <- fpc(cpue1 = res_sim, cpue2 = sol_sim, method = method, nboot = 0, kapp_zeros = "ind")
      
      # join vessel cpues to one vector
      ## uncorrected
      survey_u = c(res_sim, sol_sim)
      ## corrected
      survey_c = c(res_sim, (sol_sim * fpc$FPC))  # scale to standard vessel
      
      # compute means of uncorrected and corrected data
      cpue_u[j] = mean(survey_u)
      cpue_c[j] = mean(survey_c)
      
    }
    mse_u[i] = mean((cpue_u - mean(cpue))^2)
    mse_c[i] = mean((cpue_c - mean(cpue))^2)
  }
  
  tibble(fpd = fpd,
         mse_u = mse_u,
         mse_c = mse_c)
  
}

### simulation plot function
f_sim_plot <- function(data, spp_name, path_prefix = "./figures/side_by_side/"){
  # plot data
  data %>%
    pivot_longer(c(mse_u, mse_c), names_to = "sim", values_to = "mse") %>%
    mutate(sim = case_when(sim == "mse_u" ~ "Uncorrected",
                           sim == "mse_c" ~ "Corrected")) %>%
    ggplot()+
    geom_point(aes(x = fpd, y = mse, color = sim), alpha = 0.3)+
    geom_smooth(aes(x = fpd, y = mse, color = sim), method = "lm",formula = y ~ poly(x, 2), se = F)+
    scale_color_manual(values = cb_palette[c(1,3)])+
    geom_vline(aes(xintercept = fpc), linetype = 2)+
    labs(x = "Fishing Power Difference", y = "MSE", color = NULL)+
    scale_x_continuous(breaks = seq(0, 2, 0.2), limits = c(0, 2)) -> x
  # save
  ggsave(paste0(path_prefix, spp_name, ".png"), plot = x, width = 5, height = 4, units = "in")
  return(x)
}

# pre ----

## compute fpc
pre %>%
  left_join(spec, by = c("tow", "vessel")) %>%
  # replace NA in tows that did not catch any crab
  replace(is.na(.), 0) %>%
  dplyr::select(-tow) %>%
  pivot_longer(c(3:ncol(.)), names_to = "spp", values_to = "cpue") %>%
  # remove pairs with total zero catch
  group_by(pair, spp) %>%
  filter(sum(cpue) > 0) %>%
  # nest by spp
  group_by(spp) %>%
  nest() %>%
  mutate(spp = factor(spp, levels = spp_order)) %>%
  arrange(spp) %>%
  # ratio fpc
  mutate(ratio = purrr::map(data, f_ratio_fpc, ves1 = "res", ves2 = "sol")) %>%
  # random block anova
  mutate(rand_block = purrr::map(data, f_rand_block)) %>%
  # kappenman
  mutate(kapp = purrr::map(data, f_kapp, ves1 = "res", ves2 = "sol")) -> pre_result

## ratio result
pre_result %>%
  dplyr::select(spp, ratio) %>%
  unnest(ratio) %T>%
  # write csv
  write_csv("./output/side_by_side/pre_ratio_fpc.csv")
  
## random block anova result
pre_result %>%
  dplyr::select(spp, rand_block) %>%
  unnest(rand_block) %>%
  dplyr::select(-lmod) %T>%
  # write csv
  write_csv("./output/side_by_side/pre_random_block_fpc.csv")
  
## kappenman result
pre_result %>%
  dplyr::select(spp, kapp) %>%
  unnest(kapp) %T>%
  # write csv
  write_csv("./output/side_by_side/pre_kapp_fpc_paired.csv")


# post ----

## compute fpc
post %>%
  left_join(spec, by = c("tow", "vessel")) %>%
  # replace NA in tows that did not catch any crab
  replace(is.na(.), 0) %>%
  dplyr::select(-tow) %>%
  replace_na(list(crab60 = 100, crab80 = 0, crab100 = 0, crab120 = 0, crab121_plus = 0)) %>%
  pivot_longer(c(3:ncol(.)), names_to = "spp", values_to = "cpue") %>%
  # remove pairs with total zero catch
  group_by(pair, spp) %>%
  filter(sum(cpue) > 0) %>%
  # nest by spp
  group_by(spp) %>%
  nest() %>%
  mutate(spp = factor(spp, levels = spp_order)) %>%
  arrange(spp) %>%
  # ratio fpc
  mutate(ratio = purrr::map(data, f_ratio_fpc, ves1 = "res", ves2 = "sol")) %>%
  # random block anova
  mutate(rand_block = purrr::map(data, f_rand_block)) %>%
  # kappenman
  mutate(kapp = purrr::map(data, f_kapp, ves1 = "res", ves2 = "sol")) -> post_result

## ratio result
post_result %>%
  dplyr::select(spp, ratio) %>%
  unnest(ratio) %T>%
  # write csv
  write_csv("./output/side_by_side/post_ratio_fpc.csv")

## random block anova result
post_result %>%
  dplyr::select(spp, rand_block) %>%
  unnest(rand_block) %>%
  dplyr::select(-lmod) %T>%
  # write csv
  write_csv("./output/side_by_side/post_random_block_fpc.csv")

## kappenman result
post_result %>%
  dplyr::select(spp, kapp) %>%
  unnest(kapp) %T>%
  # write csv
  write_csv("./output/side_by_side/post_kapp_fpc.csv")


# size comp ---- 

## size comp without sex
size_comp %>%
  # filter pre sponson
  filter(year < 2019) %>%
  #bin to 3mm increments
  mutate(bin = round(isize / 3) * 3) %>%
  # combine sexes
  group_by(vessel, sponson, bin) %>%
  summarise(ncrab = sum(ncrab)) %>%
  # plot
  ggplot()+
  geom_line(aes(x = bin, y = ncrab, color = vessel))+
  scale_color_manual(values = c(1, 2))+
  scale_y_continuous(labels = comma)+
  #facet_wrap(~sponson, ncol = 1, scales = "free")+
  labs(x = "Carapace Width (mm)", y = "Number of crab", color = NULL) -> x

ggsave("./figures/side_by_side/pre_mod_sizecomp.png", height = 3, width = 4, units = "in")

## cumulative length frequency distribution
size_comp %>%
  filter(year < 2019) %>%
  #bin to 3mm increments
  mutate(bin = round(isize / 10) * 10) %>%
  group_by(sponson, vessel) %>%
  mutate(total = sum(ncrab, na.rm = T)) %>%
  group_by(sponson, vessel, total, bin) %>%
  summarise(ncrab = sum(ncrab)) %>%
  mutate(cumsum = cumsum(ncrab),
         cdf = cumsum / total) -> tmp
  
  # compute difference
  ungroup(tmp) %>%
    dplyr::select(vessel, bin, cdf) %>%
    pivot_wider(names_from = vessel, values_from = cdf) %>%
    # remove weird NA row
    filter(!is.na(bin)) %>%
    # fill cdf value at largest size for res (not observed)
    replace(is.na(.), 1) %>%
    rename_at(2:3, ~c("res", "sol")) %>%
    mutate(diff = res - sol,
           rate_diff = lead(diff) - diff) %>%
    arrange(-rate_diff) 
  

# plot cdf
tmp %>%
  ggplot()+
  geom_point(aes(x = bin, y = cdf, shape = vessel), alpha = 0.5)+
  geom_line(aes(x = bin, y = cdf, linetype = vessel))+
  scale_x_continuous(breaks = seq(0, 200, 25))+
  labs(x = "Carapace Width", y = "Cumulative Length Frequency", shape = NULL, linetype = NULL) -> x
ggsave("figures/side_by_side/pre_mod_size_comp_cdf.png", plot = x, height = 3, width = 4, units = "in")


# sim ----
## pre-modification simlulation
pre %>%
  left_join(spec, by = c("tow", "vessel")) %>%
  # use only resolution cpue
  filter(vessel == "res") %>%
  # replace NA in tows that did not catch any crab
  replace(is.na(.), 0) %>%
  # remove unneeded columns
  ungroup() %>%
  dplyr::select(-tow, -vessel, -pair) %>%
  # pivot to long format
  pivot_longer(c(1:ncol(.)), names_to = "spp", values_to = "cpue") %>%
  # nest by species group
  group_by(spp) %>%
  nest() %>%
  # get number of hauls
  # run sim
  mutate(n_tow = purrr::map_dbl(data, nrow),
         sim = purrr::map2(data, n_tow, f_kapp_mse_sim, method = 4)) -> mse_sim_kapp
  saveRDS(mse_sim_kapp, "./output/side_by_side/sim_results_pre.RDS")
  
## plots
  readRDS("./output/side_by_side/sim_results_pre.RDS") %>%
    # join to kappenman estimate
    left_join(read_csv("./output/side_by_side/pre_kapp_fpc.csv"), by = "spp") %>%
    # trim data
    dplyr::select(spp, fpc, sim) %>%
    # add fpc to the nested data
    unnest(sim) %>%
    nest() %>%
    # plot data 
    mutate(plot = purrr::map2(data, spp, f_sim_plot))

## mse simulation using ANOVA FPC
  pre %>%
    left_join(spec, by = c("tow", "vessel")) %>%
    # use only resolution cpue
    filter(vessel == "res") %>%
    # replace NA in tows that did not catch any crab
    replace(is.na(.), 0) %>%
    # remove unneeded columns
    ungroup() %>%
    dplyr::select(-tow, -vessel, -pair) %>%
    # pivot to long format
    pivot_longer(c(1:ncol(.)), names_to = "spp", values_to = "cpue") %>%
    # nest by species group
    group_by(spp) %>%
    nest() %>%
    
    filter(spp %in% c("tot_fem", "tot_male", "sl92_114")) %>%
    
    # get number of hauls
    # run sim
    mutate(n_tow = purrr::map_dbl(data, nrow),
           sim = purrr::map2(data, n_tow, f_kapp_mse_sim, method = 2)) -> mse_sim_anova
  
  ## plots
  mse_sim_anova %>%
    # join to kappenman estimate
    left_join(read_csv("./output/side_by_side/pre_random_block_fpc.csv"), by = "spp") %>%
    # trim data
    dplyr::select(spp, fpc, sim) %>%
    # add fpc to the nested data
    unnest(sim) %>%
    nest() %>%
    # plot data 
    mutate(plot = purrr::map2(data, spp, f_sim_plot, path_prefix = "./figures/side_by_side/rb_anova"))
  