# notes ----
## Tanner crab CSA (~Collie et al 2005), fit to Westward Tanner crab
## kodiak Eastside only
## Tyler Jackson
## 9/6/2021

# load ----
library(tidyverse)
library(TMB)

## custom R functinos
source("./code/csa/csa_functions.R")


## compile and load models
### model with estimated pre-recruit natural mortality
compile("./code/csa/tanner_csa_a.cpp")
dyn.load(dynlib("./code/csa/tanner_csa_a"))

### model with estimated molt probability
compile("./code/csa/tanner_csa_b.cpp")
dyn.load(dynlib("./code/csa/tanner_csa_b"))

# read input data ----
k_eastside <- read_csv("./data/csa/Kodiak_Eastside_in.csv")

# explore base model fits ----

## build input data list 
csa_data <- list(
  ## survey years
  survey_yrs = k_eastside$year,
  ## number of stages
  nstage = 3,
  ## retain catch in numbers
  ret_cat_num = k_eastside$retained_num,
  ## survey cpue index (numbers)
  index = k_eastside %>%
    dplyr::select(pre_r, r, post_r) %>%
    as.matrix()/1e6,
  # tau cs and s
  tau_cs = k_eastside$tau_cs,
  tau_s = k_eastside$tau_s[-1],
  # natural mortality
  M = 0.3,
  # pre_recruit molt probability
  molt = rep(1, nrow(k_eastside)-1),
  # data weighting on annual survey index
  wt_survey = c(1, 1, 1)
)

## parameter starting values
par_start <- list(ln_index_init = rep(1, 3), 
                  ln_rec = rep(1, length(k_eastside$year)-1),
                  ln_preM = log(0.4),
                  ln_q = log(1e-8))


## builds and fit model
obj <- MakeADFun(csa_data, par_start, DLL = "tanner_csa_a")
opt <- nlminb(start = obj$par, obj = obj$fn, gr = obj$gr)

## parameter estimates with se
par_est <- summary.sdreport(sdreport(obj))

## save parameter estimates
par_est %>%
   as.data.frame(row.names = rownames(.)) %>%
   rownames_to_column() %>%
   mutate(par = c("ln_N_init", "ln_R_init", "ln_P_init", paste0("ln_N_", 1989:2020), "ln_N_natM", "ln_q")) %>%
write_csv("./output/csa/eastside_1988_2020_par_est_model_a.csv")

## extract report
obj$report() %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename_all(~c("pre_r_fit","r_fit", "post_r_fit")) %>%
  mutate(year = k_eastside$year) -> rep

## plot fits
f_csa_plots(model = obj, 
            par_est = par_est, 
            input = k_eastside %>%
                      mutate_at(2:7, function(x){x/1e6}),
            dir = "./figures/csa/eastside_1988_2020_model_a")

## optimize recruit and post recruit natural mortality 
M_range <- seq(0.1, 0.9, 0.05)
obj <- NULL
for(i in 1:length(M_range)){
  tmp <- f_update_csa("tanner_csa_a", csa_data, par_start, M = M_range[i]) 
  obj[i]<- ifelse(tmp$opt$convergence == 0, as.numeric(tmp$opt$objective), NA)
}
tibble(M = M_range,
       obj = obj) %>%
  ggplot()+
  geom_point(aes(x = M, y = obj))+
  geom_line(aes(x = M, y = obj))+
  labs(x = "Natural Mortality (M)", y = "Objective") -> x
ggsave("./figures/csa/eastside_1988_2020_model_a/M_profile.png", plot = x, height = 3, width = 3, units = "in")


## optimize pre-r molt probability
molt_range <- seq(0.05, 1, 0.05)
molt <- list()
for(i in 1:length(molt_range)){
molt[i]<- list(rep(molt_range[i], nrow(k_eastside-1)))
}
obj <- NULL
preM <- NULL
q <- NULL
for(i in 1:length(molt)){
  tmp <- f_update_csa("tanner_csa_a", csa_data, par_start, M = 0.6, molt = molt[[i]]) 
  obj[i]<- ifelse(tmp$opt$convergence == 0, as.numeric(tmp$opt$objective), NA)
  preM[i]<- ifelse(tmp$opt$convergence == 0, tmp$opt$par[length(tmp$opt$par)-1], NA)
  q[i]<- ifelse(tmp$opt$convergence == 0, tmp$opt$par[length(tmp$opt$par)], NA)
}

tibble(molt = molt_range,
       obj = obj) %>%
  ggplot()+
  geom_point(aes(x = molt, y = obj))+
  geom_line(aes(x = molt, y = obj))+
  labs(x = "Pre-R Molt Probability", y = "Objective") -> x1

tibble(molt = molt_range,
       preM = preM) %>%
  ggplot()+
  geom_point(aes(x = molt, y = exp(preM)))+
  geom_line(aes(x = molt, y = exp(preM)))+
  labs(x = "Pre-R Molt Probability", y = "Pre-R Nat Mortality (M)") -> x2

tibble(molt = molt_range,
       q = q) %>%
  ggplot()+
  geom_point(aes(x = molt, y = exp(q)))+
  geom_line(aes(x = molt, y = exp(q)))+
  labs(x = "Pre-R Molt Probability", y = "Catchability (q)") -> x3

## effect of molt prob on fits
axis <- tickr(k_eastside, year, 5)

ggsave("./figures/csa/eastside_1988_2020_model_a/molt_prob_profile.png", plot = cowplot::plot_grid(x1, x2, x3, nrow = 1),
       height = 3, width = 7, units = "in")




# fit final models ----
## models that estimate phi outide the model and preM within
modA.1 <- f_update_csa("tanner_csa_a", csa_data, par_start)
modA.2 <- f_update_csa("tanner_csa_a", csa_data, par_start, M = 0.6) # high!
modA.3 <- f_update_csa("tanner_csa_a", csa_data, par_start, M = 0.6, molt = k_eastside$phi[-nrow(k_eastside)]) 
## models that estimate phi inside the model and fix preM = M
csa_data_b <- csa_data
csa_data_b$molt <- NULL
par_start_b <- par_start
par_start_b$ln_preM <- NULL
par_start_b$molt <- 0.2
modB.1 <- f_update_csa("tanner_csa_b", csa_data_b, par_start_b)
modB.2 <- f_update_csa("tanner_csa_b", csa_data_b, par_start_b, M = 0.6)

tibble(model = c("modA.1", "modA.2", "modA.3", "modB.1", "modB.2"),
       out = list(modA.1, modA.2, modA.3, modB.1, modB.2)) %>%
  mutate(par_est = purrr::map(out, function(x){x$par_est}),
         rep = purrr::map(out, function(x){x$rep})) -> mods

## table of parameter estimate by model
mods %>%
  filter(grepl("modA", model)) %>%
  mutate(par = purrr::map(par_est, function(x){x %>% as.data.frame(row.names = rownames(.)) %>%
                            rownames_to_column() %>%
                            mutate(par = c("ln_N_init", "ln_R_init", "ln_P_init", paste0("ln_N_", 1989:2020), "ln_N_natM", "ln_q"))})) %>%
  bind_rows(mods %>%
              filter(grepl("modB", model)) %>%
              mutate(par = purrr::map(par_est, function(x){x %>% as.data.frame(row.names = rownames(.)) %>%
                  rownames_to_column() %>%
                  mutate(par = c("ln_N_init", "ln_R_init", "ln_P_init", paste0("ln_N_", 1989:2020), "phi", "ln_q"))}))) %>%
  dplyr::select(model, par) %>%
  unnest(par) %>%
  # select only estimate, names and model
  dplyr::select(model, par, Estimate) %>%
  # reformat
  pivot_wider(names_from = model, values_from = Estimate) %>%
  # write output
  write_csv("./output/csa/eastside_1988_2020_par_est_models.csv")

## plot fits by model
mods %>%
  unnest(rep) %>%
  dplyr::select(model, year, pre_r_fit, r_fit, post_r_fit) %>%
  mutate(mature_fit = pre_r_fit + r_fit + post_r_fit,
         legal_fit = r_fit + post_r_fit) -> rep

## pre_recruit fit
k_eastside %>%
  ggplot()+
  geom_point(aes(x = year, y = pre_r/1e6), color = "grey40")+
  geom_errorbar(aes(x = year, ymax = pre_r/1e6 + 1.96*pre_r_se/1e6, ymin = pre_r/1e6 - 1.96*pre_r_se/1e6), 
                color = "grey40", width = 0)+
  geom_line(data = rep, aes(x = year, y = pre_r_fit, color = model))+
  scale_x_continuous(labels = axis$labels, breaks = axis$breaks)+
  coord_cartesian(ylim = c(0, NA))+
  labs(x = NULL, y = "Pre-Recruit Index", color = NULL)+
  theme(legend.position = c(0.1, 0.8)) -> x1
ggsave("./figures/csa/eastside_1988_2020_pre_r_all_models.png", plot = x1,
       height = 3, width = 5, units = "in")
## recruit fit
k_eastside %>%
  ggplot()+
  geom_point(aes(x = year, y = r/1e6), color = "grey40")+
  geom_errorbar(aes(x = year, ymax = r/1e6 + 1.96*r_se/1e6, ymin = r/1e6 - 1.96*r_se/1e6), 
                color = "grey40", width = 0)+
  geom_line(data = rep, aes(x = year, y = r_fit, color = model))+
  scale_x_continuous(labels = axis$labels, breaks = axis$breaks)+
  coord_cartesian(ylim = c(0, NA))+
  labs(x = NULL, y = "Recruit Index", color = NULL)+
  theme(legend.position = c(0.1, 0.8)) -> x2
ggsave("./figures/csa/eastside_1988_2020_r_all_models.png", plot = x2,
       height = 3, width = 5, units = "in")
## post recruit fit
k_eastside %>%
  ggplot()+
  geom_point(aes(x = year, y = post_r/1e6), color = "grey40")+
  geom_errorbar(aes(x = year, ymax = post_r/1e6 + 1.96*post_r_se/1e6, ymin = post_r/1e6 - 1.96*post_r_se/1e6), 
                color = "grey40", width = 0)+
  geom_line(data = rep, aes(x = year, y = post_r_fit, color = model))+
  scale_x_continuous(labels = axis$labels, breaks = axis$breaks)+
  coord_cartesian(ylim = c(0, NA))+
  labs(x = NULL, y = "Pre-Recruit Index", color = NULL)+
  theme(legend.position = c(0.1, 0.8)) -> x3
ggsave("./figures/csa/eastside_1988_2020_post_r_all_models.png", plot = x3,
       height = 3, width = 5, units = "in")
# combined plot
ggsave("./figures/csa/eastside_1988_2020_post_r_all_models.png", 
       plot = cowplot::plot_grid(x1, x2, x3, ncol = 1),
       height = 9, width = 5, units = "in")


## legal fit
k_eastside %>%
  mutate(legal = r + post_r,
         legal_se = sqrt(r_se^2 + post_r_se^2)) %>%
  ggplot()+
  geom_point(aes(x = year, y = legal/1e6), color = "grey40")+
  geom_errorbar(aes(x = year, ymax = legal/1e6 + 1.96*legal_se/1e6, ymin = legal/1e6 - 1.96*legal_se/1e6), 
                color = "grey40", width = 0)+
  geom_line(data = rep, aes(x = year, y = legal_fit, color = model))+
  scale_x_continuous(labels = axis$labels, breaks = axis$breaks)+
  coord_cartesian(ylim = c(0, NA))+
  labs(x = NULL, y = "Legal Index", color = NULL)+
  theme(legend.position = c(0.1, 0.8)) -> x
ggsave("./figures/csa/eastside_1988_2020_legal_all_models.png", plot = x,
       height = 3, width = 5, units = "in")
## mature fit
k_eastside %>%
  mutate(mature = pre_r + r + post_r,
         mature_se = sqrt(pre_r_se^2 + r_se^2 + post_r_se^2)) %>%
  ggplot()+
  geom_point(aes(x = year, y = mature/1e6), color = "grey40")+
  geom_errorbar(aes(x = year, ymax = mature/1e6 + 1.96*mature_se/1e6, ymin = mature/1e6 - 1.96*mature_se/1e6), 
                color = "grey40", width = 0)+
  geom_line(data = rep, aes(x = year, y = mature_fit, color = model))+
  scale_x_continuous(labels = axis$labels, breaks = axis$breaks)+
  coord_cartesian(ylim = c(0, NA))+
  labs(x = NULL, y = "Mature Index", color = NULL)+
  theme(legend.position = c(0.1, 0.8)) -> x
ggsave("./figures/csa/eastside_1988_2020_mature_all_models.png", plot = x,
       height = 3, width = 5, units = "in")

## objective function table
mods %>%
  # extract stuff
  mutate(obj = purrr::map_dbl(out, function(x){x$opt$objective}),
         convergence = purrr::map_dbl(out, function(x){x$opt$convergence})) %>%
  dplyr::select(model, obj, convergence) %>%
  write_csv("./output/csa/eastside_1988_2020_objective.csv")
