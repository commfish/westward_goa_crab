# This script analyzes Resolution/Solstice comparison data based on the randomized block
# design methodology described in the ROP. The first part produces the statistics of interest
# and writes them to a csv file and the second part generates normal qq plots for the residuals. 
# Datasets sbs2015-17.csv and sbs2019.csv must be treated separately.
# w.gaeuman 2-12-2021

# load required packages
library(tidyverse)

# want sum contrasts for linear model (R default is treatment contrasts)
options(contrasts = rep("contr.sum", 2))

# CHANGE INPUT DATA NAME IN THIS PART AS NEEDED
# ---------------------------------------------
# import and restructure data; apply log(x + 1) transform to all CPUEs; assumes tow pairs are 
# in adjacent rows
data0 <- read.csv("sbs2019.csv", as.is = TRUE) %>%
	transmute(pair = as.character(rep(1:(dim(.)[1]/2), each = 2)),
		vessel = ifelse(vessel_id == 30, "R", "S"),
		juv.fem = juv_fem, adult.fem = adult_fem, tot.fem = tot_fem, 
		sub70 = sl_.70, sub70.91 = sl_70.91, sub92.114 = sl_92.114, submature = sl_.114,
		recruit = recruit, legal = tot_legal, mature = tot_mature, tot.male = tot_male,
		arrowtooth = arrowtooth, flathead = flathead, pcod = pacific_cod, pollock = pollock) %>%
	mutate_if(is.numeric, function(x) log(x + 1))

# A. STATISTICS OF INTEREST
# -------------------------
# define analysis funtion for variable of interest (expects quoted variable name, e.g. "tot.fem")
f <- function(data, quoted.var){
	# use only pairs with at least one nonzero catch
	good.pairs <- group_by(data, pair) %>%
		summarize(x = sum(eval(as.name(quoted.var)))) %>%
		filter(x > 0) %>%
		pull(pair)
	data <- filter(data, pair %in% good.pairs)
	# fit linear model and obtain quantities of interest
	fit <- lm(eval(as.name(quoted.var)) ~ pair + vessel, data = data)
	coefs <- summary(fit)[["coefficients"]]
	n = length(good.pairs) # number of hauls
	mu = coefs[1, 1] # grand mean
	nu = coefs[nrow(coefs), 1] # Resolution treatment effect
	nu.se = coefs[nrow(coefs), 2] # treatment effect standard error
	nu.p = coefs[nrow(coefs), 4] # treatment effect p-value
	fpc = exp(2*nu*(1 + 0.5*nu.se^2)) # FPC factor (Resolution/Soltice)
	low.95 = exp(2*nu - 1.96*2*nu.se); hi.95 = exp(2*nu + 1.96*2*nu.se) # FPC CI
	# return everything in a single vector with named components
	return(c(n = n, mu = mu, nu = nu, nu.se = nu.se, nu.p = nu.p, 
		fpc = fpc, low.95 = low.95, hi.95 = hi.95))
}

# apply function in a loop to all variables of interest, collecting results in matrix 'results'
results <- numeric()
for(quoted.var in names(data0)[3:17])
	results <- rbind(results, f(data0, quoted.var))
# give rows analysis variable names
row.names(results) <- names(data0)[3:17]

# write results to appropriately named csv file
 write.csv(results, "results2019.csv")

# B. RESIDUAL PLOTS 
# -----------------
# (this part is independent of A above)

# define function to obtain residuals (reproduces a lot of function f above)
g <- function(data, quoted.var){
	# use only pairs with at least one nonzero catch
	good.pairs <- group_by(data, pair) %>%
		summarize(x = sum(eval(as.name(quoted.var)))) %>%
		filter(x > 0) %>%
		pull(pair)
	data <- filter(data, pair %in% good.pairs)
	# fit linear model and obtain quantities of interest
	fit <- lm(eval(as.name(quoted.var)) ~ pair + vessel, data = data)
	resids <- fit$res
	return(resids)
}

# loop through variables of interest and collect residuals in list 'all.resids'
all.resids <- list()
for(quoted.var in names(data0)[3:17])
	all.resids[[quoted.var]] <- g(data0, quoted.var)

# convert list to a dataframe for plotting
all.resids <- stack(all.resids)
names(all.resids) <- c("resids", "variable")

# plot using separate facets for each variable
ggplot(all.resids, aes(sample = resids)) +
	stat_qq() + stat_qq_line(color = "red") +
	theme_bw() +
	facet_wrap( ~ variable, ncol = 3, scales = "free_y") 

ggsave("resids2019.png")











