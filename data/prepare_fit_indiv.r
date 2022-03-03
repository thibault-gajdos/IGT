## Create the following data
## summary.all (summary_indiv.rdata): parametres des fits par sujet et modèle
## pred.all (pred_indiv.rdata): predictions des fits par essai, sujet et modèle (choix le plus fréquent)
## obs.all (obs.rdata): observations des choix par essai et sujet
## accuracy (accuracy.rdata) (compute prediction accuracy (1 = correct, 0 = incorrect)
## by trial for each subject and model)

rm(list=ls(all=TRUE))  ## efface les données

source('~/thib/projects/tools/R_lib.r')
setwd('~/thib/projects/IGT/data/')

load('./cluster/data_all/data_all.rdata')
names <- c(unique(data.all$study),'mice', 'human')


## predictions and fitted parameters
params.pvl = c('A', 'alpha', 'cons', 'lambda')
params.pvl.noise = c('A', 'alpha', 'cons', 'lambda', 'zeta')
params.vse = c('alpha', 'cons', 'gamma', 'delta', 'phi')
params.vse.noise = c('alpha', 'cons', 'gamma', 'delta', 'phi','zeta')
mat = matrix(ncol = 0, nrow = 0)
summary.all = data.frame(mat)
pred.all = data.frame(mat)


for (m in c('vse','vse_noise','pvl_delta_lambda', 'pvl_delta_lambda_noise')){
    if (m == 'vse'){
            params = params.vse
        }else if (m == 'vse_noise'){
            params = params.vse.noise
        }else if (m == 'pvl_delta_lambda'){
            params = params.pvl
        }else{
            params = params.pvl.noise
        }
    for (name in names){
        if (name %in% c('human','mice')){
            nsubjs <- 40
        }else{
            nsubjs <- length(unique(data.all[data.all$study == name,]$subjID))
        }
        for (i in c(1:nsubjs)){
            data.name <- paste('./cluster/output/',m,'_',name,'_',i,'.rdata', sep='')
            load(data.name)
            print(data.name)

            ## prediction
            pred <- extract(fit)$y_pred %>%
                               as.data.frame() %>%
                               summarise_all(mlv, method = 'mfv') %>%
                               mutate(study = name, model = m, subjID = i)
            pred.all <-  bind_rows(pred.all, pred[1,])

            ## log_lik and parameters
            log_lik = loo::extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)
            l <- loo::elpd(log_lik)$estimates[[1]]
            s <- as.data.frame(summary(fit, pars = params)$summary) %>%
                rownames_to_column() %>%
                rename(param = rowname) %>%
                mutate(l = l, study = name, model = m, subjID = i)
            summary.all <-  bind_rows(summary.all, s)
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

for (name in names){
    ## load data
    if (name == 'human') {
        load('./cluster/data_human/data_human.rdata')
        data <- data.human
    }else if (name == 'mice'){
        load('./cluster/data_mice/data_mice.rdata')
        data <- data.mice
    }else{
        load('./cluster/data_all/data_all.rdata')
        data <- data.all %>% filter(study == name)
    }
    if (name %in% c('human','mice')){
        nsubjs <- 40
    }else{
        nsubjs <- length(unique(data.all[data.all$study == name,]$subjID))
    }
    for (i in c(1:nsubjs)){
        obs <- data %>%
            filter(subjID == i)
        obs_vector  <- t(as.vector(obs$choice))
        obs <- as.data.frame(obs_vector) %>%
            mutate(study = name, subjID = i)
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
    relocate(c(model, study, subjID), .before = 1)
## create an outlier variable
## = 1 if max(Rhat)>1.2 (fit did not converge), 0 otherwise
outlier <- summary.all %>%
    group_by(model, study, subjID) %>%
    mutate(Rmax = max(Rhat, na.rm = T)) %>%
    mutate(outlier = ifelse(Rmax>1.2, 1,0)) %>%
    filter(param == 'alpha') %>%
    select(model, study, subjID, outlier) %>%
    ungroup()
out  <- outlier %>%
    group_by(model, study) %>%
    summarise(outlier = sum(outlier, na.rm = TRUE)) %>%
    pivot_wider(names_from = model, values_from = outlier)



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
