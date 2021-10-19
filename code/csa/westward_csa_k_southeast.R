# notes ----
## Tanner crab CSA (~Collie et al 2005), fit to Westward Tanner crab
## kodiak Eastside only
## Tyler Jackson
## 9/15/2021

# load ----
library(tidyverse)
library(TMB)

## custom R functinos
source("./code/csa/csa_functions.R")


## compile and load models
### model with estimated pre-recruit natural mortality
compile("./code/csa/tanner_csa_a.cpp")
dyn.load(dynlib("./code/csa/tanner_csa_a"))

# read input data ----
k_southeast <- read_csv("./data/csa/Kodiak_Southeast_in.csv")

# fit models ----

## build input data list 
kse_csa_data <- list(
  ## survey years
  survey_yrs = k_southeast$year,
  ## number of stages
  nstage = 3,
  ## retain catch in numbers
  ret_cat_num = k_southeast$retained_num,
  ## survey cpue index (numbers)
  index = k_southeast %>%
    dplyr::select(pre_r, r, post_r) %>%
    replace_na(list(pre_r = 1, r = 1, post_r = 1)) %>%
    as.matrix()/1e6,
  # tau cs and s
  tau_cs = k_southeast$tau_cs,
  tau_s = k_southeast$tau_s[-1],
  # natural mortality
  M = 0.6,
  # pre_recruit molt probability
  molt = rep(1, length(k_southeast$year)-1),
  # data weighting on annual survey index
  wt_survey = c(1, 1, 1)
)

## parameter starting values
kse_par_start <- list(ln_index_init = rep(1, 3), 
                      ln_rec = rep(1, length(k_southeast$year)-1),
                      ln_preM = log(0.3),
                      ln_q = log(1e-15))


## builds and fit model
kse_obj <- MakeADFun(kse_csa_data, kse_par_start, DLL = "tanner_csa_a")
kse_opt <- nlminb(start = kse_obj$par, obj = kse_obj$fn, gr = kse_obj$gr)

## parameter estimates with se
kse_par_est <- summary.sdreport(sdreport(kse_obj))

# par_est %>%
#   as.data.frame(row.names = rownames(.)) %>%
#   rownames_to_column() %>%
#   mutate(par = c("ln_N_init", "ln_R_init", "ln_P_init", paste0("ln_N_", 1989:2020), "ln_N_natM", "ln_q")) %>%
# write_csv("./output/csa/eastside_1988_2020_par_est.csv")

## extract report
kse_obj$report() %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename_all(~c("pre_r_fit","r_fit", "post_r_fit")) %>%
  mutate(year = k_eastside$year) -> ke_rep

## plot fits
f_csa_plots(model = kse_obj, 
            par_est = kse_par_est, 
            input = k_southeast %>%
              mutate_at(2:7, function(x){x/1e6}),
            dir = "./figures/csa/southeast_1988_2020")

# check optimization of M
M_range <- seq(0.1, 0.9, 0.05)
obj <- NULL
for(i in 1:length(M_range)){
  tmp <- f_update_csa("tanner_csa_a", kse_csa_data, kse_par_start, M = M_range[i]) 
  obj[i]<- ifelse(tmp$opt$convergence == 0, as.numeric(tmp$opt$objective), NA)
}
tibble(M = M_range,
       obj = obj) %>%
  ggplot()+
  geom_point(aes(x = M, y = obj))+
  geom_line(aes(x = M, y = obj))+
  labs(x = "Natural Mortality (M)", y = "Objective") -> x
ggsave("./figures/csa/southeast_1988_2020/M_profile.png", plot = x, height = 3, width = 3, units = "in")

## print parameters
kse_par_est %>%
  as.data.frame(row.names = rownames(.)) %>%
  rownames_to_column() %>%
  mutate(par = c("ln_N_init", "ln_R_init", "ln_P_init", paste0("ln_N_", 1989:2020), "ln_N_natM", "ln_q")) %>%
  write_csv("./output/csa/southeast_1988_2020_par_est.csv")

