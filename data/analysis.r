rm(list=ls(all=TRUE))## efface les données
source('~/thib/projects/tools/R_lib.r')
setwd('~/thib/projects/IGT/data/')
load('summary_indiv.rdata') ## summary.all
load( 'pred_indiv.rdata') ## pred.all
load('obs.rdata') ## obs.al
load('accuracy.rdata') ## obs.all

## * fit
setwd('~/thib/projects/IGT/data/cluster/')
source('fit_util.r')
load('data_igt.rdata')
data <- data.igt %>%
    filter(study == 'mice', subjID == 17)
fit_indiv(data = data, model = 'pvl_delta_lambda_noise', delta = .99,  treedepth = 12, name = 'bis')
load('~/thib/projects/IGT/data/cluster/output/pvl_delta_lambda_noise_6_17.rdata')
launch_shinystan(fit)
p <- pairs(fit, pars = c('A','alpha','cons','persev', 'zeta'))
p
setwd('~/thib/projects/IGT/data/')

## * fit quality

fit.quality <- summary.all %>%
      group_by(model, study, subjID) %>%
      mutate(Rmax = max(Rhat, na.rm = T)) %>%
      filter(param == 'alpha') %>%
      select(model, study, subjID, Rmax, divergent) %>%
      mutate(outlier = ifelse(Rmax>1.2 | divergent>0, 1,0)) %>%
      ungroup()

dd <- fit.quality %>%
    group_by(model, study) %>%
    summarise(outlier = sum(outlier), Rmax = sum(Rmax-1))
dd
print(dd, n = 200)
## * Accuracy

a <- accuracy %>%
    rowwise() %>%
    mutate(acc = mean(c_across(starts_with("V")), na.rm = TRUE), .keep = "unused") %>%
    ungroup()
accuracy <- merge(a, fit.quality)

accuracy.summary <- accuracy %>%
    filter(outlier == 0) %>%
    group_by(model,study) %>%
    summarise(acc = mean(acc, na.rm = TRUE)*100)
a.summary <- accuracy.summary %>%
    pivot_wider(names_from = model, values_from = acc)
print(kable(a.summary, digits = 2))

## * Parmeters
d <- merge(summary.all, fit.quality)
d.param <-  d %>%
    mutate(species = ifelse(study == 'mice', 'mice', 'human')) %>%
    filter(outlier == 0) %>%
    select(model, param, mean, study, species) %>%
    group_by(param, model, species) %>%
    summarise(mean = mean(mean, na.rm = TRUE)) %>%
    pivot_wider(names_from = param, values_from = mean)
print(kable(d.param, digit = 2))
dd <- d %>% filter(model == 'vse', param == 'cons', study == 'mice')
hist(dd$mean)



