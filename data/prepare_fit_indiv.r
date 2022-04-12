## Create the following data
## summary.all (summary_indiv.rdata): parametres des fits par sujet et modèle
## pred.all (pred_indiv.rdata): predictions des fits par essai, sujet et modèle (choix le plus fréquent)
## obs.all (obs.rdata): observations des choix par essai et sujet
## accuracy (accuracy.rdata) (compute prediction accuracy (1 = correct, 0 = incorrect)
## by trial for each subject and model)

rm(list=ls(all=TRUE))  ## efface les données

source('~/thib/projects/tools/R_lib.r')
setwd('~/thib/projects/IGT/data/')

load('./cluster/data_igt.rdata')
models <- c('pvl_delta_lambda', 'pvl_delta_lambda_noise', 'vse', 'vse_noise')
range.d <- c(1:12)
names <- data.frame(study = unique(data.igt$study), d = unique(data.igt$d))

## predictions and fitted parameters
params.pvl = c('A', 'alpha', 'cons', 'persev')
params.pvl.noise = c('A', 'alpha', 'cons', 'zeta','persev')
params.vse = c('alpha', 'cons', 'gamma', 'delta', 'phi')
params.vse.noise = c('alpha', 'cons', 'gamma', 'delta', 'phi','zeta')
mat = matrix(ncol = 0, nrow = 0)
summary.all = data.frame(mat)
pred.all = data.frame(mat)

##for (m in c('vse','vse_noise','pvl_delta_lambda', 'pvl_delta_lambda_noise')){
for (m in models){
    if (m == 'vse'){
            params = params.vse
        }else if (m == 'vse_noise'){
            params = params.vse.noise
        }else if (m == 'pvl_delta_lambda'){
            params = params.pvl
        }else{
            params = params.pvl.noise
        }
    for (s in range.d){
        nsubjs <- length(unique(data.igt[data.igt$d == s,]$subjID))
        for (i in c(1:nsubjs)){
            data.name <- paste('./cluster/output/',m,'_',s,'_',i,'.rdata', sep='')
            load(data.name)
            print(data.name)

            ## prediction
            pred <- extract(fit)$y_pred %>%
                               as.data.frame() %>%
                               summarise_all(mlv, method = 'mfv') %>%
                               mutate(d = s, study = names[names$d == s,]$study, model = m, subjID = i)
            pred.all <-  bind_rows(pred.all, pred[1,])

            ## log_lik and parameters
            log_lik = loo::extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)
            l <- loo::elpd(log_lik)$estimates[[1]]
            sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
            div = sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
            su <- as.data.frame(summary(fit, pars = params)$summary) %>%
                rownames_to_column() %>%
                rename(param = rowname) %>%
                mutate(l = l, d =s, study = names[names$d == s,]$study, model = m, subjID = i, divergent = div)
            summary.all <-  bind_rows(summary.all, su)
        }
    }
}

summary.all <- summary.all %>%
    mutate(model = case_when(model == 'pvl_delta_lambda' ~ 'pvl',
                             model == 'pvl_delta_lambda_noise' ~ 'pvl_noise',
                             model == 'vse' ~ 'vse',
                             model == 'vse_noise' ~ 'vse_noise'))
pred.all <- pred.all %>%
    mutate(model = case_when(model == 'pvl_delta_lambda' ~ 'pvl',
                             model == 'pvl_delta_lambda_noise' ~ 'pvl_noise',
                             model == 'vse' ~ 'vse',
                             model == 'vse_noise' ~ 'vse_noise'))

save(summary.all, file = 'summary_indiv.rdata')
save(pred.all, file = 'pred_indiv.rdata')


## observed choices
mat = matrix(ncol = 0, nrow = 0)
obs.all = data.frame(mat)

load('./cluster/data_igt.rdata')
for (s in range.d){
    data <- data.igt %>% filter(d == s)
    nsubjs <- length(unique(data.igt[data.igt$d == s,]$subjID))

    for (i in c(1:nsubjs)){
        obs <- data %>%
            filter(subjID == i)
        obs_vector  <- t(as.vector(obs$choice))
        obs <- as.data.frame(obs_vector) %>%
            mutate(study =  names[names$d == s,]$study, subjID = i)
        if (nrow(obs.all) == 0){
            obs.all <- obs
        }else{
            obs.all <-  bind_rows(obs.all, obs)
        }
    }
}

save(obs.all, file = 'obs.rdata')

## * Accuracy
pred.all <- pred.all %>%
    select(-d) %>%
    relocate(c(model, study,  subjID), .before = 1)


## compute prediction accuracy (1 = correct, 0 = incorrect)
## by trial for each subject and model
accuracy <- pred.all %>% filter(V1 ==  -1) ## create empty dataframe
obs <- obs.all %>% select(starts_with('V'))
for (m in unique(pred.all$model)){
    pred <- pred.all %>%
        filter(model == m) %>%
        select(starts_with('V'))
    diff <- as.data.frame(pred == obs)
    diff <- add_column(diff, model = m, study = obs.all$study, subjID = obs.all$subjID, .before = 1)
    accuracy <- bind_rows(accuracy, diff)
}
save(accuracy, file = 'accuracy.rdata' )
