# notes ----
## functions for tanner crab csa models
## Tyler Jackson

# functions ---- needs documentation

f_extract_rec <- function(par_est, survey_yrs){
  par_est %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(grepl("rec", rowname)) %>%
    mutate(year = survey_yrs[-1]) %>%
    rename(ln_rec = Estimate, se = `Std. Error`) %>%
    dplyr::select(year, ln_rec, se) -> out
  return(out)
}
f_csa_plots <- function(model, par_est, input, dir) {
  
  # load plot style library from B. Williams (GitHub)
  library(FNGr)
  theme_set(theme_sleek())
  axis <- tickr(input, year, 5)
  
  # create directory
  if(!exists(dir)){dir.create(dir)}
  
  # reformat report and par----
  model$report() %>%
    as.data.frame() %>%
    as_tibble() %>%
    rename_all(~c("pre_r_fit","r_fit", "post_r_fit")) %>%
    mutate(year = input$year) -> rep
  
  par_est %>% 
    as_tibble() %>%
    rename_all(~c("est", "se")) %>%
    mutate(par = rownames(par_est)) -> par
  
  # index fit plots ----
  ## pre recruits
  ggplot()+
    geom_point(data = input, aes(x = year, y = pre_r), color = "grey40")+
    geom_errorbar(data = input, aes(x = year, ymin = pre_r - 1.96*pre_r_se, ymax = pre_r + 1.96*pre_r_se), 
                  width = 0, color = "grey40")+
    geom_line(data = rep, aes(x = year, y = pre_r_fit), color = "blue")+
    scale_x_continuous(labels = axis$labels, breaks = axis$breaks)+
    coord_cartesian(ylim = c(0, NA))+
    labs(x = NULL, y = "Pre-Recruit Index") -> p1
  ggsave(file.path(dir, "pre_recruit_fit.png"), plot = p1, width = 5, height = 3, units = "in")
  ## recruits
  ggplot()+
    geom_point(data = input, aes(x = year, y = r), color = "grey40")+
    geom_errorbar(data = input, aes(x = year, ymin = r - 1.96*r_se, ymax = r + 1.96*r_se), 
                  width = 0, color = "grey40")+
    geom_line(data = rep, aes(x = year, y = r_fit), color = "blue")+
    scale_x_continuous(labels = axis$labels, breaks = axis$breaks)+
    coord_cartesian(ylim = c(0, NA))+
    labs(x = NULL, y = "Recruit Index") -> p2
  ggsave(file.path(dir, "recruit_fit.png"), plot = p2, width = 5, height = 3, units = "in")
  ## post recruits
  ggplot()+
    geom_point(data = input, aes(x = year, y = post_r), color = "grey40")+
    geom_errorbar(data = input, aes(x = year, ymin = post_r - 1.96*post_r_se, ymax = post_r + 1.96*post_r_se), 
                  width = 0, color = "grey40")+
    geom_line(data = rep, aes(x = year, y = post_r_fit), color = "blue")+
    scale_x_continuous(labels = axis$labels, breaks = axis$breaks)+
    coord_cartesian(ylim = c(0, NA))+
    labs(x = NULL, y = "Post-Recruit Index") -> p3
  ggsave(file.path(dir, "post_recruit_fit.png"), plot = p3, width = 5, height = 3, units = "in")
  
  ## all size classes, one plot
  ggsave(file.path(dir, "all_index_fit.png"), plot = cowplot::plot_grid(p1, p2, p3, nrow = 3), 
         width = 5, height = 9, units = "in")
  
  # mature male abundance ----
  input %>%
    mutate(mature = pre_r + r + post_r,
           mature_se = sqrt(pre_r_se^2 + r_se^2 + post_r_se^2)) %>%
    ggplot()+
    geom_point(aes(x = year, y = mature), color = "grey40")+
    geom_errorbar(aes(x = year, ymin = mature - 1.96*mature_se, ymax = mature + 1.96*mature_se), 
                  width = 0, color = "grey40")+
    geom_line(data = rep, aes(x = year, y = pre_r_fit + r_fit + post_r_fit), color = "blue")+
    scale_x_continuous(labels = axis$labels, breaks = axis$breaks)+
    coord_cartesian(ylim = c(0, NA))+
    labs(x = NULL, y = "Mature Males (mil)") -> x
  ggsave(file.path(dir, "mature_males_fit.png"), plot = x, width = 5, height = 3, units = "in")
  
  # legal male abundance ----
  input %>%
    mutate(legal = r + post_r,
           legal_se = sqrt(r_se^2 + post_r_se^2)) %>%
    ggplot()+
    geom_point(aes(x = year, y = legal), color = "grey40")+
    geom_errorbar(aes(x = year, ymin = legal - 1.96*legal_se, ymax = legal + 1.96*legal_se), 
                  width = 0, color = "grey40")+
    geom_line(data = rep, aes(x = year, y = r_fit + post_r_fit), color = "blue")+
    scale_x_continuous(labels = axis$labels, breaks = axis$breaks)+
    coord_cartesian(ylim = c(0, NA))+
    labs(x = NULL, y = "Legal Males (mil)") -> x
  ggsave(file.path(dir, "legal_males_fit.png"), plot = x, width = 5, height = 3, units = "in")
  
  # recruitment ----
  f_extract_rec(par_est = par_est, survey_yrs = input$year) %>%
    ggplot()+
    geom_point(aes(x = year, y = exp(ln_rec)), color = "blue")+
    geom_line(aes(x = year, y = exp(ln_rec)), color = "blue")+
    scale_x_continuous(labels = axis$labels, breaks = axis$breaks)+
    coord_cartesian(ylim = c(0, NA))+
    labs(x = NULL, y = "Recruitment (mil)") -> x
  ggsave(file.path(dir, "annual_rec_est.png"), plot = x, width = 5, height = 3, units = "in")
  
}
f_update_csa <- function(model_name, input_data, par_start, survey_yrs = NULL, M = NULL, molt = NULL, wt_survey = NULL){
  
  # update input data ----
  ## timeseries
  if(!is.null(survey_yrs)){
    ## get the original timeseries yrs
    og_yrs <- input_data$survey_yrs
    ## update timseries yrs
    input_data$survey_yrs <- survey_yrs
    ## update retained catch timeseries
    tibble(yr = og_yrs,
           catch = input_data$ret_cat_num) %>%
      filter(yr %in% survey_yrs) %>%
      pull(catch) -> input_data$ret_cat_num
    ## update index timseries
    input_data$index %>%
      as_tibble() %>%
      mutate(yr = og_yrs) %>%
      filter(yr %in% survey_yrs) %>%
      dplyr::select(-yr) %>%
      as.matrix() -> input_data$index
    ## update taus
    tibble(yr = og_yrs,
           tau_cs = input_data$tau_cs) %>%
      filter(yr %in% survey_yrs) %>%
      pull(tau_cs) -> input_data$tau_cs
    tibble(yr = og_yrs,
           tau_s = c(NA, input_data$tau_s)) %>%
      filter(yr %in% survey_yrs) %>%
      pull(tau_s) %>% 
      na.omit() %>%
      as.numeric() -> input_data$tau_s
    ## update recruitment starting values
    tibble(yr = og_yrs,
           ln_rec = c(NA, par_start$ln_rec)) %>%
      filter(yr %in% survey_yrs) %>%
      pull(ln_rec) -> tmp
    tmp[1] <- NA  
    na.omit(tmp) %>%
      as.numeric() -> par_start$ln_rec
  }
  ## recruit natural mortality
  if(!is.null(M)){input_data$M <- M}
  ## recruit natural mortality
  if(!is.null(molt)){input_data$molt <- molt}
  ## objective function weights
  if(!is.null(wt_survey)){input_data$wt_survey <- wt_survey}
  
  # fit model ----
  obj_update <- MakeADFun(input_data, par_start, DLL = model_name)
  opt_update <- nlminb(start = obj_update $par, obj = obj_update$fn, gr = obj_update $gr)
  
  out = list(
    obj = obj_update,
    opt = opt_update,
    par_est = summary.sdreport(sdreport(obj_update)),
    rep = obj_update$report() %>%
      as.data.frame() %>%
      as_tibble() %>%
      rename_all(~c("pre_r_fit","r_fit", "post_r_fit")) %>%
      mutate(year = input_data$survey_yrs)
  )
  
  return(out)
}
