rm(list=ls(all=TRUE))  ## efface les données


source('~/thib/projects/tools/R_lib.r')
setwd('~/thib/projects/IGT/data/data_mice/')
## * Prepare data
## load data
gain <-  read_excel('cons_mice_subjects_x_trials.xlsx', col_names = FALSE)
loss <-  read_excel('loss_mice_subjects_x_trials.xlsx', col_names = FALSE)
choice <-  read_excel('choice_mice_subjects_x_trials.xlsx', col_names = FALSE)

N <- nrow(choice) ## number of subjects
T <- ncol(choice) ## max number of trials by subject
gain <- as.data.frame(gain) %>%
    pivot_longer(cols = everything(), names_to = "trial", values_to = "gain")
loss <- as.data.frame(loss)  %>%
    pivot_longer(cols = everything(), names_to = "trial", values_to = "loss") %>%
    mutate(loss = min(abs(loss), 1)) ## replace loss by 1
choice <- as.data.frame(choice) %>%
    pivot_longer(cols = everything(), names_to = "trial", values_to = "choice")
subjID <- sort(rep(1:N,T))
data <- data.frame(subjID = subjID, trial = gain$trial, outcome = gain$gain - abs(loss$loss), choice = choice$choice)

## Missing data
which(is.na(data), arr.ind=TRUE) ##check for NA
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
#noise    <- array(-1, c(N, T))

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

## ** fit pvld
modelFile <- '~/thib/projects/IGT/data/igt_plv_delta_noise.stan'
fit.mice <- stan(modelFile,
            data = data_list,
            iter = 10000,
            warmup = 4000,
            chains = 4,
            init = "random",
            control = list(adapt_delta = 0.95,  max_treedepth = 12)
            )

save(fit.mice, file = 'fit_mice_noise.rdata')
