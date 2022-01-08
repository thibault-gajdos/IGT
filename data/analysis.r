
start = Sys.time()
fit <- stan(modelFile,
            data    = dataList,
            iter    = nIter,
            chains  = 1,
            warmup  = nWarmup,
            init    =  "random",
            seed    = 123,
            control = list(adapt_delta = 0.8,  max_treedepth = 12))
end = Sys.time()
t = end - start
t

x <- extract(fit)
glimpse(x)
x <- x$weight
z <- x[500,,,]
w <- z[1,,]
tail(w)
glimpse(x)
save(fit, file = 'fit_hierarchical_pvld.rda'),
list_of_draws <- extract(fit)
print(names(list_of_draws))
params <- c('mu_alpha',  'mu_c',  'mu_theta')
posterior <- as.array(fit)
np <- nuts_params(fit)

mcmc_parcoord(posterior, np = np)

color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior, pars = params, np = np) +
  xlab("Post-warmup iteration")

mcmc_pairs(posterior, np = np, pars = params,
           off_diag_args = list(size = 0.75))
dev.off()
x11()
res=150, width = 1000, height = 1000)
par(mar = rep(0, 4))
pairs(fit, par = params, res=150, width = 1000, height = 1000)
fit_summary <- summary(fit, par = params, probs = c(0.05, 0.95))$summary


## ** fit hierarch
n_trials <- ncol(choice)
n_subjects <- nrow(choice)

dataList <- list(n_trials = n_trials,
                 n_subjects = n_subjects,
                 choice = choice,
                 gain = gain,
                 loss = loss)

modelFile <- '../pvld_hierarchical.stan'
nIter     <- 10000
nChains   <- 4
nWarmup   <- 5000


##start = Sys.time()
fit <- stan(modelFile,
            data    = dataList,
            chains  = nChains,
            iter    = nIter,
            warmup  = nWarmup,
            init    =  "random",
            seed    = 123,
            control = list(adapt_delta = 0.99,  max_treedepth = 12))

##end = Sys.time()
##t = end - start


x <- extract(fit)
glimpse(x)
x <- x$weight
z <- x[500,,,]
w <- z[1,,]
tail(w)
glimpse(x)
save(fit, file = 'fit_hierarchical_pvld.rda'),
list_of_draws <- extract(fit)
print(names(list_of_draws))
params <- c('mu_alpha',  'mu_c',  'mu_theta')
posterior <- as.array(fit)
np <- nuts_params(fit)

mcmc_parcoord(posterior, np = np)

color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior, pars = params, np = np) +
  xlab("Post-warmup iteration")

mcmc_pairs(posterior, np = np, pars = params,
           off_diag_args = list(size = 0.75))
dev.off()
x11()
res=150, width = 1000, height = 1000)
par(mar = rep(0, 4))
pairs(fit, par = params, res=150, width = 1000, height = 1000)
fit_summary <- summary(fit, par = params, probs = c(0.05, 0.95))$summary

## * fit pvld (Steingroever et al.)
## data : subjects x trials
## ** Prepare data
gain <-  read_excel('win_mice_subjects_x_trials.xlsx', col_names = FALSE)
loss <-  read_excel('loss_mice_subjects_x_trials.xlsx', col_names = FALSE)
conso <-  read_excel('cons_mice_subjects_x_trials.xlsx', col_names = FALSE)   ##subj 7, trial 57: missing, replaced by gain
choice <-  read_excel('choice_mice_subjects_x_trials.xlsx', col_names = FALSE)
unit <- matrix(1, 40, 100)

## replace loss by min(loss,1)
for (i in c(1:nrow(loss))){
    for (j in c(1:ncol(loss))){
        loss[i,j] = min(loss[i,j], 1)
    }
}
net <- conso-loss

## choice = 0
nul <- as.data.frame(which(choice == 0, arr.ind=TRUE)) %>%
    mutate(gain = -1, loss = -1, conso = -1)
n = nrow(nul)
for (i in c(1:n)){
    nul$gain[i] = gain[nul$row[i], nul$col[i]]
    nul$conso[i] = conso[nul$row[i], nul$col[i]]
    nul$loss[i] = loss[nul$row[i], nul$col[i]]
}
nul

## ** Fit
n_s <- nrow(choice)
n_t <- ncol(choice)

dataList <- list(n_t = n_trials,
                 n_s = n_subjects,
                 choice = choice,
                 net = net)

modelFile <- '../pvl_d.stan'
nIter     <- 15000
nChains   <- 4
nWarmup   <- 10000
nThin     <- 1

start = Sys.time()
fit <- stan(modelFile,
            data    = dataList,
            chains  = nChains,
            iter    = nIter,
            warmup  = nWarmup,
            thin    = nThin,
            init    = "random",
            seed    = 123,
            control = list(adapt_delta = 0.95,  max_treedepth = 12))
t = end - start
t
save(fit, file = 'fit.rda'),

## * brouillon
load('fit.dta')
list_of_draws <- extract(fit)
print(names(list_of_draws))
params <- c('mu_alpha',  'mu_c', 'mu_w', 'mu_theta')
fit_summary <- summary(fit, par = params, probs = c(0.05, 0.95))$summary
print(fit_summary)
warnings(fit)
plot.summary <- stan_plot(fit, pars = params, show_outer_line = FALSE, ci_level = .95)
plot.summary
summary(fit)
fit_main <- extract(fit, pars = params)
sample <- posterior_samples(fit_main)

plot_alpha  <- post_plot(sample$mu_alpha) + xlab('Learning (alpha)')
plot_beta  <- post_plot(sample$mu_beta) + xlab('inverse temperature (beta)')
plot_theta  <- post_plot(sample$mu_theta) + xlab('utility (theta)')
plot_delta  <- post_plot(sample$mu_delta) + xlab('discounting (delta)')
plot_phi  <- post_plot(sample$mu_phi) + xlab('exploration (phi)')

plot_hist  <- ggarrange(plot_alpha , plot_beta, plot_theta, plot_delta, plot_phi,
              labels = 'AUTO', ncol = 2, nrow = 3)
plot_hist

mcmc_areas(
  sample,
  pars = params,
  prob = 0.95, # 80% intervals
  prob_outer = 0.95, # 99%
  point_est = "mean"
)



##
y_pred <- extract(fit, 'y_pred')
y_pred <- unlist(y_pred, use.names=FALSE)
str(y_pred)
plot(density(y_pred))

## * individual

## data
gain <- read.csv2('win_Steingroever2011.csv', header = F) %>%  t()
loss <- read.csv2('loss_Steingroever2011.csv', header = F) %>% t()
choice <- read.csv2('choice_Steingroever2011.csv', header = F) %>% t()

n_trials  <- nrow(gain)
n <- 22
dataList <- list(n_trials = n_trials,
                 choice = choice[,n],
                 gain = gain[,n]/100,
                 loss = abs(loss[,n])/100)

modelFile <- 'IGT_indiv.stan'
modelFile <- 'PVL_noise.stan'
modelFile <- 'PVL_nonoise.stan'
## modelFile <- 'pvl.stan'

nIter     <- 14000
nChains   <- 4
nWarmup   <- 10000
nThin     <- 1

start = Sys.time()
fit <- stan('PVL_noise.stan',
            data    = dataList,
            chains  =  nChains,
            iter    = nIter,
            warmup  = nWarmup,
            thin    = nThin,
            init    = "random",
            seed    = 123,
            control = list(adapt_delta = 0.95,  max_treedepth = 12))
end = Sys.time()
t = end - start

params <- c('alpha', 'c',  'theta', 'w', 'zeta', 'log_lik')
fit_summary <- summary(fit, par = params, probs = c(0.05, 0.95))$summary
print(fit_summary)

fit.nonoise <- stan('PVL_nonoise.stan',
            data    = dataList,
            chains  =  nChains,
            iter    = nIter,
            warmup  = nWarmup,
            thin    = nThin,
            init    = "random",
            seed    = 123,
            control = list(adapt_delta = 0.95,  max_treedepth = 12))
params.nonoise <- c('alpha', 'c', 'w', 'theta', 'log_lik')
fit.nonoise_summary <- summary(fit.nonoise, par = params.nonoise, probs = c(0.05, 0.95))$summary

print(fit.nonoise_summary)
print(fit_summary)

posterior <- as.matrix(fit.nonoise)
mcmc_areas(posterior,
           pars = c("zeta"),
           prob = 0.95)
mcmc_areas(posterior,
           pars = c("c"),
           prob = 0.95)
mcmc_areas(posterior,
           pars = c("w"),
           prob = 0.95)

predicted <- extract(fit,  permuted = TRUE, inc_warmup = FALSE, include = TRUE)
sim <- predicted$choice_sim[1,]
alpha <- predicted$alpha
describe(alpha)
c <- predicted$c
describe(c)
delta <- predicted$delta
describe(delta)
theta <- predicted$theta
describe(theta)
Ev <- predicted$Ev[1,]
head(Ev)
Ev_sim <- predicted$Ev_sim[1,]
Ev_sim
u <- predicted$u
u
sim
str(sim)
color_scheme_set("darkgray")
posterior_fit <- as.array(fit)
np_fit <- nuts_params(fit)
mcmc_parcoord(posterior_fit, np = np_fit, pars = c("alpha", "beta", "theta"))
pair <- mcmc_pairs(
  fit,
  pars = c("alpha", "c", "theta"),
  off_diag_args = list(size = 3/4, alpha = 1/3), # size and transparency of scatterplot points
  np_style = pairs_style_np(div_color = "black", div_shape = 2) # color and shape of the divergences
)
pair

save(fit, file = 'fit_pvl.rda')
save(fit, file = 'fit_pvl_noise.rda')
pairs(fit)
pairs(fit, pars = c("alpha", "beta"))

mcmc_pairs(posterior_fit, np = np_fit, pars = c("alpha","c","theta"),
           off_diag_args = list(size = 0.75))

list_of_draws <- extract(fit)
print(names(list_of_draws))
pairs(fit, pars = c("alpha", "c", 'theta'))

## params <- c('alpha', 'beta', 'delta', 'w', 'theta', 'zeta', 'log_lik')

predicted <- extract(fit,  permuted = TRUE, inc_warmup = FALSE, include = TRUE)
glimpse(predicted)


x <- as.matrix(fit, pars = 'log_lik')


plot.summary <- stan_plot(fit, pars = params, show_outer_line = FALSE, ci_level = .95)
plot.summary


dataList <- list(n_trials = n_trials,
                 #n_subjects = n_subjects,
                 choice = choice[,1],
                 gain = conso[,1],
                 loss = abs(loss[,1]))
