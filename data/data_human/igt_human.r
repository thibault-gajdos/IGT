rm(list=ls(all=TRUE))  ## efface les données
source('~/thib/projects/tools/R_lib.r')
setwd('~/thib/projects/IGT/data/data_human/')

## * prepare data
## load data
gain <-  read_excel('win_human_subjects_x_trials.xlsx', col_names = FALSE)
loss <-  read_excel('loss_human_subjects_x_trials.xlsx', col_names = FALSE)
choice <-  read_excel('choice_human_subjects_x_trials.xlsx', col_names = FALSE)

N <- nrow(choice) ## number of subjects
T <- ncol(choice) ## max number of trials by subject
gain <- as.data.frame(gain) %>%
    pivot_longer(cols = everything(), names_to = "trial", values_to = "gain") %>%
    mutate(gain = gain/100)
loss <- as.data.frame(loss)  %>%
    pivot_longer(cols = everything(), names_to = "trial", values_to = "loss") %>%
    mutate(loss = loss/100)
choice <- as.data.frame(choice) %>%
    pivot_longer(cols = everything(), names_to = "trial", values_to = "choice")
subjID <- sort(rep(1:N,T))
data <- data.frame(subjID = subjID, trial = gain$trial, outcome = gain$gain - abs(loss$loss), choice = choice$choice)

## Missing data
data <- data %>%
    filter(choice > 0)

## build data list
N <- length(unique(data$subjID))
T <- length(unique(data$trial))
## compute trials by subject
d <- data %>%
    group_by(subjID) %>%
    summarise(t_subjs = n())
t_subjs <- d$t_subjs

## Initialize data arrays
Ydata    <- array(-1, c(N, T))
RLmatrix <- array( 0, c(N, T))

## Write from raw_data to the data arrays
for (i in 1:N) {
    t <- t_subjs[i]
    data_subj <- data %>% filter(subjID == i)
    Ydata[i, 1:t]    <- data_subj$choice
    RLmatrix[i, 1:t] <- data_subj$outcome
}

## Wrap into a list for Stan
data_list <- list(
    N        = N,
    T        = T,
    Tsubj    = t_subjs,
    choice   = Ydata,
    outcome  = RLmatrix
    )

## * fit pvld noise
modelFile <- '~/thib/projects/IGT/data/igt_plv_delta_noise.stan'
fit.human <- stan(modelFile,
            data = data_list,
            iter = 10000,
            warmup = 5000,
            chains = 4,
            init = "random",
            control = list(adapt_delta = 0.90,  max_treedepth = 13)
            )
save(fit.human, file = 'fit_human_noise.rdata')

## list_of_draws <- extract(fit)
## print(names(list_of_draws))
## posterior <- as.array(fit)
## np <- nuts_params(fit)

## mcmc_parcoord(posterior, np = np)

## color_scheme_set("mix-brightblue-gray")
## mcmc_trace(posterior, pars = params, np = np) +
##   xlab("Post-warmup iteration")

## mcmc_pairs(posterior, np = np, pars = params,
##            off_diag_args = list(size = 0.75))

## params_indiv <- c('alpha', 'A', 'cons', 'lambda', 'zeta')
## params_group <- c('mu_alpha', 'mu_A', 'mu_cons', 'mu_lambda', 'mu_zeta')
## fit_summary_group <- summary(fit, par = params_group, probs = c(0.05, 0.95))$summary
## fit_summary_indiv <- summary(fit, par = params_indiv, probs = c(0.05, 0.95))$summary
## fit_summary_group
## fit_summary_indiv
## pairs(fit, par = params, res=150, width = 1000, height = 1000)
