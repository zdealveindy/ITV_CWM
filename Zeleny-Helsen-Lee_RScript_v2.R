# Extending the CWM approach to intraspecific trait variation: how to deal with overly optimistic standard tests?
# Oecologia, Methods paper
# Author of the script: David Zeleny

# Install and upload libraries ----
#devtools::install_github ('zdealveindy/simcom')
#devtools::install_github ('zdealveindy/weimea')

library (weimea)  # https://github.com/zdealveindy/weimea
library (simcom)  # https://github.com/zdealveindy/simcom
library (vegan)
library (parallel)
library (tidyverse)
library (abind)

# Definition of functions ----

# cwm_ss - calculates site-specific CWM from "traits" (matrix with rows = sites, cols = species)
cwm_ss <- function (com, traits_ss, wstand = FALSE)
{
  com <- as.matrix (com)
  traits_ss <- as.matrix (traits_ss)
  pij <- com/rowSums (com, na.rm = TRUE)
  ci <- rowSums (traits_ss*pij, na.rm = TRUE)
  return (ci)
}

# cwm_f - calculates fixed CWM from "traits" (matrix with rows = sites, cols = species)
cwm_f <- function (com, traits_ss, wstand = FALSE) # calculates fixed CWM values from traits; wstand not implemented yet
{
  com <- as.matrix(com)
  traits_f <- colMeans (traits_ss, na.rm = TRUE)
  cwm_temp <- rowSums (t(apply (com, 1, FUN = function (x) x/sum  (x)) * traits_f), na.rm = T)
  return (cwm_temp)
}

# cwm_itv - calculates intraspecific CWM from "traits" (matrix with rows = sites, cols = species)
cwm_itv <- function (com, traits_ss, wstand = FALSE) # calculates ITV CWM values from traits; wstand not implemented yet
{
  com <- as.matrix(com)
  mean_t <- colMeans (traits_ss, na.rm = TRUE)
  delta_t <- sweep (traits_ss, 2, mean_t, FUN = '-')
  cwm_temp <- rowSums (t(apply (com, 1, FUN = function (x) x/sum  (x))) * delta_t, na.rm = T)
  return(cwm_temp)
}

# sample_traits_ss - permutes t_ij matrix (intraspecific trait variation matrix with site specific trait values)
# Arguments:
# traits_ss - matrix of site-specific trait values (T); sites in rows, species in columns
# intra - logical (default FALSE), should intraspecific trait values be permutated? (column-wise permutations in the matrix with intraspecific trait values, ??T)
# inter - logical (default FALSE), should interspecific trait values be permuted (elements in the vector mean_t)
# rand - logical (default FALSE), should all (non-missing) values in traits_ss be randomly permuted?
# replace_by_norm - logical (default FALSE), should all (non-missing) values in traits_ss be replaced by those randomly drawn from normal distribution (mean = 0, sd = 1)?
sample_traits_ss <- function (traits_ss, intra = FALSE, inter = FALSE, rand = FALSE, replace_by_norm = FALSE) # permutes traits in t_ij
{
  t_ij_perm <- traits_ss  # in case that both inter and intra stay FALSE
  
  if (rand) {
    if (replace_by_norm) t_ij_perm[!is.na (t_ij_perm)] <- rnorm (n = sum (!is.na (t_ij_perm)), mean = 0, sd = 1) else
      t_ij_perm[!is.na (t_ij_perm)] <- sample (t_ij_perm[!is.na (t_ij_perm)])
  }else {
    
    if (inter & !intra) {
      t_mean <- colMeans (traits_ss, na.rm = TRUE)
      dt_ij <- sweep (traits_ss, 2, FUN = '-', t_mean)
      t_mean_perm <- sample (t_mean)
      t_ij_perm <- sweep (dt_ij, 2, FUN = '+', t_mean_perm)
    }
    if (intra & !inter) t_ij_perm <- apply (traits_ss, 2, FUN = function (x) {tr_r <- x; tr_r[!is.na(x)] <- x[!is.na(x)][sample (length (x[!is.na(x)]))]; tr_r}) # alt: sample (x[!is.na(x)], size = length (x[!is.na(x)]), replace = FALSE)
    
    if (intra & inter) {
      t_mean <- colMeans (traits_ss, na.rm = TRUE)
      dt_ij <- sweep (traits_ss, 2, FUN = '-', t_mean)
      t_mean_perm <- sample (t_mean)
      t_ij_perm <- sweep (dt_ij, 2, FUN = '+', t_mean_perm)
      t_ij_perm <- apply (t_ij_perm, 2, FUN = function (x) {tr_r <- x; tr_r[!is.na(x)] <- x[!is.na(x)][sample (length (x[!is.na(x)]))]; tr_r}) # alt: sample (x[!is.na(x)], size = length (x[!is.na(x)]), replace = FALSE)
    }
  }
  return (t_ij_perm)
}
  
# gen_traits_ss - generates t_ij matrix from com and env simulated data
gen_traits_ss <- function (com, env, scale = TRUE, m = 1, jitter = NULL, amount = NULL)
{
  t_ij_orig <- sweep (ifelse (com>0,1,NA), 1, env, FUN = '*')
  mean_tij <- colMeans (t_ij_orig, na.rm = T)
  dt_ij <- sweep (t_ij_orig, 2, mean_tij, FUN = '-')
  t_ij <- sweep (dt_ij*m, 2, mean_tij, FUN = '+')
  if (!is.null (jitter)) t_ij <- jitter (t_ij, jitter, amount)
  if (scale) t_ij <- (t_ij - min (t_ij, na.rm = T))/(max (t_ij, na.rm = T) - min (t_ij, na.rm = T))
  t_ij
}

# Calculates intraspecific vs interspecific trait variation ratio from traits_ss
intra_inter_ratio <- function (traits_ss) {
  sd_intra <- apply (traits_ss, 2, sd, na.rm = TRUE)
  t_mean <- colMeans (traits_ss, na.rm = TRUE)
  sd_inter <- sd (t_mean, na.rm = T)
  return (mean (sd_intra/sd_inter, na.rm = TRUE))
}

# plot_ss and plot_reg - plots trait-env relationship with intraspecific trait variation (plot_reg is just helper function)
plot_reg <- function (x, y){
  LM <- lm (y ~ x)
  pred <- predict (LM, newdata = list (x = c(min (x), max (x))))
  lines (x = c(min (x), max (x)), y = pred)
}

# plot_ss was previously called plot_itv
plot_ss <- function (com, traits_ss, env, ...)
{
  t_e <- lapply (1:ncol (com), FUN = function (x) {t <- traits_ss[com[,x]>0, x]; e <- env[com[,x]>0]; data.frame (e = e, t = t)})
  t_range <- range (traits_ss, na.rm = TRUE)
  plot.new ()
  plot.window (xlim = c(min (env), max (env)), ylim = c(t_range[1], t_range[2]))
  axis (1)
  axis (2)
  title (xlab = 'env', ylab = 'trait')
  box (bty = 'l')
  #lapply (t_e, FUN = function (sp) points (t ~ e, sp))
  lapply (t_e, FUN = function (sp) if (nrow (sp)>0) plot_reg (sp$e, sp$t))
}

plot_itv <- function (com, traits_ss, env)
{
  t_e <- lapply (1:ncol (com), FUN = function (x) {t <- traits_ss[com[,x]>0, x]; t_mean <- mean (t, na.rm = TRUE); e <- env[com[,x]>0]; data.frame (e = e, t_itv = t-t_mean)})
  t_itv_range <- range (unlist (lapply (t_e, FUN = function (x) x$t_itv)), na.rm = TRUE)
  plot.new ()
  plot.window (xlim = c(min (env), max (env)), ylim = c(t_itv_range[1], t_itv_range[2]))
  axis (1)
  axis (2)
  title (xlab = 'env', ylab = 'trait')
  box (bty = 'l')
  #lapply (t_e, FUN = function (sp) points (t ~ e, sp))
  lapply (t_e, FUN = function (sp) if (nrow (sp)>0) plot_reg (sp$e, sp$t_itv))
}

plot_f <- function (com, traits_ss, env)
{
  t_e <- lapply (1:ncol (com), FUN = function (x) {t <- traits_ss[com[,x]>0, x]; t_mean <- mean (t, na.rm = TRUE); t_fix <- rep (t_mean, length (t)); e <- env[com[,x]>0]; data.frame (e = e, t_fix = t_fix)})
  t_range <- range (traits_ss, na.rm = TRUE)
  plot.new ()
  plot.window (xlim = c(min (env), max (env)), ylim = c(t_range[1], t_range[2]))
  axis (1)
  axis (2)
  title (xlab = 'env', ylab = 'trait')
  box (bty = 'l')
  #lapply (t_e, FUN = function (sp) points (t ~ e, sp))
  lapply (t_e, FUN = function (sp) if (nrow (sp)>0) plot_reg (sp$e, sp$t_fix))
}

# test_cwm_ss - max permutation test of site-specific CWM ~ env regression
test_cwm_ss <- function (com, traits_ss, env, perm = 199)
{
  CWM_ss_obs <- cwm_ss (com = com, traits_ss = traits_ss)
  LM <- lm (scale (CWM_ss_obs) ~ scale (env))
  an <- anova (LM)
  su <- summary (LM)
  F_obs <- su$fstatistic[1]

  mean_t <- colMeans (traits_ss, na.rm = TRUE)
  delta_T <- apply (traits_ss, 2, FUN = function (x) x - mean (x, na.rm = T))
  
  F_col_rand <- replicate (perm, expr = {
    traits_ss_rand <- sample_traits_ss (traits_ss = traits_ss, intra = TRUE, inter = TRUE)
    CWM_ss_rand <- cwm_ss (com, traits_ss_rand)
    summary (lm (scale (CWM_ss_rand) ~ scale (env)))$fstatistic[1]
  })
  F_col_all <- c(F_col_rand, F_obs)
  P_col_i <- sum (F_col_all >= F_obs)/(perm+1)
  
  F_row_rand <- replicate (perm, expr = {
    summary (lm (scale (CWM_ss_obs) ~ scale (sample (env))))$fstatistic[1]
  })
  F_row_all <- c(F_row_rand, F_obs)
  P_row <- sum (F_row_all >= F_obs)/(perm+1)
  P_max <- max (P_row, P_col_i)
  
  return (list (SSR = an$`Sum Sq`[1], SSE = an$`Sum Sq`[2], b1 = su$coefficients[2,1], b1_se = su$coefficients[2,2], F_obs = F_obs, r2 = su$r.squared, P_par = an$`Pr(>F)`[1], P_row = P_row, P_col_i = P_col_i, P_max = P_max, perm = perm))
}

# test_cwm_f - max permutation test of fixed CWM ~ env regression (using weimea package)
test_cwm_f <- function (com, traits_ss, env, perm = 199) 
{
  require (weimea)  
  traits_f <- colMeans (traits_ss, na.rm = TRUE)
  CWM_f <- cwm (com = com, traits = traits_f)
  LM <- lm (scale (CWM_f$traits) ~ scale (env))  # calculates standardized regression
  an <- anova (LM)
  su <- summary (LM)
  res_test_cwm <- test_cwm (cwm = CWM_f, env = env, method = 'lm', test = c('par', 'max'), perm = perm)
  return (list (SSR = an$`Sum Sq`[1], SSE = an$`Sum Sq`[2], b1 = su$coefficients[2,1], b1_se = su$coefficients[2,2], F_obs = res_test_cwm$out$F, r2 = res_test_cwm$out$r2, P_par = res_test_cwm$out$P_par, P_max = res_test_cwm$out$P_max, perm = perm))
}

# test_cwm_itv - parametric test of site-specific CWM ~ env regression  
test_cwm_itv <- function (com, traits_ss, env, perm = 199) 
{
  CWM_itv_obs <- cwm_itv (com = com, traits_ss = traits_ss)
  LM <- lm (scale (CWM_itv_obs) ~ scale (env))
  an <- anova (LM)
  su <- summary (LM)

  return (list (SSR = an$`Sum Sq`[1], SSE = an$`Sum Sq`[2], b1 = su$coefficients[2,1], b1_se = su$coefficients[2,2], F_obs = su$fstatistic[1], r2 = su$r.squared, P_par = an$`Pr(>F)`[1], perm = perm))
}
  
# Helper function: permutes elements in vector of any length (incl 1), from package gtools
permute <- function (x) {
  sample(x, size = length(x), replace = FALSE)
}

# Creates simple simulated data using the simcom package and plots them ----
sim <- simul.comm (totS = 50, min.niche.breath = 500, max.niche.breath = 5000) # 50 species
sam <- sample.comm (sim, Np = 25) # 25 sites

com <- sam$a.mat
env <- sam$sample.x

traits_ss <- gen_traits_ss (com = com, env = env, m = .2, scale = FALSE, jitter = 0, amount = NULL)
traits_fixed <- colMeans (traits_ss, na.rm = T)

CWM_ITV <- cwm_ss (com = com, traits_ss = traits_ss)
plot (CWM_ITV ~ env)

CWM_SS_noITV <- cwm_ss (com = com, traits_ss = sample_traits_ss (traits_ss, intra = T, inter = F))
plot (CWM_SS_noITV ~ env)

CWM_fixed <- cwm (com, traits_fixed)
plot (CWM_fixed$traits ~ env)

plot_ss (com = com, traits_ss = traits_ss, env = env)
plot_itv (com = com, traits_ss = traits_ss, env = env)
plot_f (com = com, traits_ss = traits_ss, env = env)
plot_ss (com = com, traits_ss = sample_traits_ss (traits_ss, intra = TRUE), env = env)  # plots intraspecific trait variation

# CWM site specific tested by combined permutation test + all other tests ----
perm <- 1000 #1000  # how many standard test per community is done
no_rep <- 50 #50  # replicate community data (with the same parameters)

set.seed (35974)

cl <- makeCluster(16)
clusterEvalQ (cl, {library (weimea); library (simcom)})
clusterExport (cl, varlist = c('cwm_ss', 'cwm_itv', 'gen_traits_ss', 'sample_traits_ss', 'perm', 'intra_inter_ratio', 'test_cwm_ss'))

P_ss_m_all_max <- parLapply (cl, 1:no_rep, fun = function (x) {
  
  sim <- simul.comm (totS = 100, min.niche.breath = 2500, max.niche.breath = 5000) # orign 50 species
  sam <- sample.comm (sim, Np = 100) # orig 25 sites
  com <- sam$a.mat
  env <- sam$sample.x
  
  P_ss_m <- lapply (seq (0, 5, by = 0.5), FUN = function (m) { 
    traits_ss <- gen_traits_ss (com = com, env = env, m = m, scale = TRUE, jitter = 1, amount = ifelse (m == 0, 0, 1000))
    iir <- intra_inter_ratio (traits_ss)
    P_ss <- sapply (1:perm, FUN = function (x) {
      traits_ss_rand <- sample_traits_ss (traits_ss, intra = T, inter = T)
      CWM_intra_rand <- cwm_ss (com = com, traits_ss = traits_ss_rand)
      result_test_cwm_ss <- test_cwm_ss (com = com, traits_ss = traits_ss_rand, env = env, perm = 199)
      CWM_itv_rand <- cwm_itv (com = com, traits_ss = traits_ss_rand)
      result_test_cwm_itv_par <- summary (lm (CWM_itv_rand ~ env))$coefficients[2,4]  
      unlist (c (m = m, perm_id = x, iir = iir, result_test_cwm_ss, P_par_itv = result_test_cwm_itv_par))
    })
  })
  
  # add completely randomized t_ij
  P_ss_rand <- sapply (1:perm, FUN = function (x) {
    traits_ss <- gen_traits_ss (com = com, env = env, m = 1, scale = T, amount = .1)
    traits_ss_all_rand <- sample_traits_ss (traits_ss, intra = F, inter = F, rand = T, replace_by_norm = T)
    CWM_intra_rand <- cwm_ss (com = com, traits_ss = traits_ss_all_rand)
    result_test_cwm_ss_all_rand <- test_cwm_ss (com = com, traits_ss = traits_ss_all_rand, env = env, perm = 199)
    CWM_itv_rand <- cwm_itv (com = com, traits_ss = traits_ss_all_rand)
    result_test_cwm_itv_par_all_rand <- summary (lm (CWM_itv_rand ~ env))$coefficients[2,4]  
    unlist (c (m = 100, perm_id = x, irr = 100, result_test_cwm_ss_all_rand, P_par_itv = result_test_cwm_itv_par_all_rand))
    })
  
  P_ss_all <- c (P_ss_m, list (P_ss_rand))  # adds results to one list
})

stopCluster (cl)

# Save the data which are result of costly calculation
# save (P_ss_m_all_max, file = 'P_ss_m_all_max_20240330.RData')
# load (file = 'P_ss_m_all_max_20240330.RData')

summary_P_ss_m_all_max <- sapply (P_ss_m_all_max, simplify = 'array', FUN = function (repl)
  sapply (repl, FUN = function (x)
    c(m = unique (x[1,]), # m
      iir = mean (x[3,]), # iir
      SSR_mean = mean (x["SSR",]),
      SSE_mean = mean (x["SSE",]),
      b1_mean = mean (x["b1",]),
      b1_se_mean = mean (x["b1_se",]),
      R2_ss_mean = mean (x["r2",]),
      P_par_IF = sum(x["P_par",]<0.05)/(ncol (x)*0.05), # P_par
      P_row_IF = sum(x["P_row",]<0.05)/(ncol (x)*0.05), # P_row
      P_col_i_IF = sum(x["P_col_i",]<0.05)/(ncol (x)*0.05), # P_col_i
      P_max_IF = sum(x["P_max",]<0.05)/(ncol (x)*0.05),
      P_par_itv_IF = sum (x["P_par_itv",]<0.05)/(ncol (x)*0.05))   # P_max
  ))

mean_P_ss_m_all_max <- apply (summary_P_ss_m_all_max, MARGIN = c(2,1), mean)
sd_P_ss_m_all_max <- apply (summary_P_ss_m_all_max, MARGIN = c(2,1), sd)


mean_P_ss_m_all_max[12,'iir'] <- 3.5  # arbitrary value for IIR of completely random scenario

# Plot the dependence of inflation index on magnitude of intraspecific trait variation ----
jpeg ('standard_max_test_cwm_ss.jpg', width = 12, height = 12, units = 'cm', res = 600, pointsize = 8)
par (mfrow = c(2,2))

# Plot results for CWM_SS and standard parametric test
plot (P_par_IF ~ iir, mean_P_ss_m_all_max, xlim = c(0, 3.5), ylim = c(0, 16), type = 'n', las = 1, xaxt = 'n', bty = 'l', ann = F)
title (xlab = list ('Relative ITV index'), ylab = expression (Inflation~index~(alpha==0.05)), line = 2.5)
title (main = list ("Site-specific CWM ~ env\nStandard parametric test", cex = 1, font = 1))
axis (1, at = seq (0, 3, by = 0.5), labels = format (seq (0, 3, by = 0.5)))
axis (1, at = 3.5, labels = expression (infinity))
abline (h = 1, col = 'grey', lty = 'dashed')
for (i in seq (1, nrow (mean_P_ss_m_all_max)))
{
  lines (x = c(mean_P_ss_m_all_max[i, 'iir'], mean_P_ss_m_all_max[i, 'iir']), y = c(mean_P_ss_m_all_max[i, 'P_par_IF'] - sd_P_ss_m_all_max[i, 'P_par_IF'], c(mean_P_ss_m_all_max[i, 'P_par_IF'] + sd_P_ss_m_all_max[i, 'P_par_IF'])), col = 'darkgrey')
}
lines (P_par_IF ~ iir, mean_P_ss_m_all_max[-12,], lty = c('solid'))
lines (P_par_IF ~ iir, mean_P_ss_m_all_max[-1:-10,], lty = c('dashed'))
points (P_par_IF ~ iir, mean_P_ss_m_all_max, pch = c(22, rep (21, 11)), bg = c('black', rep ('gray', 10), 'white'))
mtext (side = 3, text = '(a)', adj = -.25, line = 2, cex = 1.1, font = 2)

# Plot results for CWM_ITV and standard parametric test
plot (P_par_itv_IF ~ iir, mean_P_ss_m_all_max, xlim = c(0, 3.5), ylim = c(0, 1.5), type = 'n', las = 1, xaxt = 'n', bty = 'l', ann = FALSE)
title (xlab = list ('Relative ITV index'), ylab = expression (Inflation~index~(alpha==0.05)), line = 2.5)
title (main = list ("Intraspecific CWM ~ env\nStandard parametric test", cex = 1, font = 1))
axis (1, at = seq (0, 3, by = 0.5), labels = format (seq (0, 3, by = 0.5)))
axis (1, at = 3.5, labels = expression (infinity))
abline (h = 1, col = 'grey', lty = 'dashed')
for (i in seq (1, nrow (mean_P_ss_m_all_max)))
{
  lines (x = c(mean_P_ss_m_all_max[i, 'iir'], mean_P_ss_m_all_max[i, 'iir']), y = c(mean_P_ss_m_all_max[i, 'P_par_itv_IF'] - sd_P_ss_m_all_max[i, 'P_par_itv_IF'], c(mean_P_ss_m_all_max[i, 'P_par_itv_IF'] + sd_P_ss_m_all_max[i, 'P_par_itv_IF'])), col = 'darkgrey')
}
lines (P_par_itv_IF ~ iir, mean_P_ss_m_all_max[c(-1,-12),], lty = c('solid'))
lines (P_par_itv_IF ~ iir, mean_P_ss_m_all_max[-1:-10,], lty = c('dashed'))
points (P_par_itv_IF ~ iir, mean_P_ss_m_all_max[-1,], pch = rep (21, 11), bg = c(rep ('gray', 10), 'white'))
mtext (side = 3, text = '(b)', adj = -.25, line = 2, cex = 1.1, font = 2)


# Plot results for CWM_ss and max test of row- and combined col-based permutation tests
plot (P_max_IF ~ iir, mean_P_ss_m_all_max, xlim = c(0, 3.5), ylim = c(0, 16), type = 'n', las = 1, xaxt = 'n', bty = 'l', ann = FALSE)
title (xlab = list ('Relative ITV index'), ylab = expression (Inflation~index~(alpha==0.05)), line = 2.5)
title (main = list ("Site-specific CWM ~ env\n'Max' permutation test", cex = 1, font = 1))
axis (1, at = seq (0, 3, by = 0.5), labels = format (seq (0, 3, by = 0.5)))
axis (1, at = seq (0, 3, by = 0.5), labels = format (seq (0, 3, by = 0.5)))
axis (1, at = 3.5, labels = expression (infinity))
abline (h = 1, col = 'grey', lty = 'dashed')
for (i in seq (1, nrow (mean_P_ss_m_all_max)))
{
  lines (x = c(mean_P_ss_m_all_max[i, 'iir'], mean_P_ss_m_all_max[i, 'iir']), y = c(mean_P_ss_m_all_max[i, 'P_max_IF'] - sd_P_ss_m_all_max[i, 'P_max_IF'], c(mean_P_ss_m_all_max[i, 'P_max_IF'] + sd_P_ss_m_all_max[i, 'P_par_IF'])), col = 'darkgrey')
}

lines (P_max_IF ~ iir, mean_P_ss_m_all_max[-12,], lty = c('solid'))
lines (P_max_IF ~ iir, mean_P_ss_m_all_max[-1:-10,], lty = c('dashed'))
points (P_max_IF ~ iir, mean_P_ss_m_all_max, pch = c(22, rep (21, 11)), bg = c('black', rep ('gray', 10), 'white'))
mtext (side = 3, text = '(c)', adj = -.25, line = 2, cex = 1.1, font = 2)

dev.off ()

# Effect sizes and ITV ----
plot (R2_ss_mean ~ iir, mean_P_ss_m_all_max, xlim = c(0, 3.5), ylim = c(0, 0.5), type = 'n', las = 1, xaxt = 'n', bty = 'l', ann = F)
title (xlab = list ('Relative ITV index'), ylab = expression (r^2), line = 2.5)
title (main = list ("Site-specific CWM ~ env\nExplained variation", cex = 1))
axis (1, at = seq (0, 3, by = 0.5), labels = format (seq (0, 3, by = 0.5)))
axis (1, at = 3.5, labels = expression (infinity))
abline (h = 1/(25-1), col = 'grey', lty = 'dashed')
for (i in seq (1, nrow (mean_P_ss_m_all_max)))
{
  lines (x = c(mean_P_ss_m_all_max[i, 'iir'], mean_P_ss_m_all_max[i, 'iir']), y = c(mean_P_ss_m_all_max[i, 'R2_ss_mean'] - sd_P_ss_m_all_max[i, 'R2_ss_mean'], c(mean_P_ss_m_all_max[i, 'R2_ss_mean'] + sd_P_ss_m_all_max[i, 'R2_ss_mean'])), col = 'darkgrey')
}
lines (R2_ss_mean ~ iir, mean_P_ss_m_all_max[-12,], lty = c('solid'))
lines (R2_ss_mean ~ iir, mean_P_ss_m_all_max[-1:-10,], lty = c('dashed'))
points (R2_ss_mean ~ iir, mean_P_ss_m_all_max, pch = c(22, rep (21, 11)), bg = c('black', rep ('gray', 10), 'white'))
mtext (side = 3, text = '(d)', adj = -.25, line = 2, cex = 1.1, font = 2)

plot (SSR_mean ~ iir, mean_P_ss_m_all_max, xlim = c(0, 3.5), ylim = c(0, 10), type = 'n', las = 1, xaxt = 'n', bty = 'l', ann = F)
title (xlab = list ('Relative ITV index'), ylab = 'SSR', line = 2.5)
title (main = list ("Site-specific CWM ~ env\nExplained sum of squares", cex = 1))
axis (1, at = seq (0, 3, by = 0.5), labels = format (seq (0, 3, by = 0.5)))
axis (1, at = 3.5, labels = expression (infinity))
for (i in seq (1, nrow (mean_P_ss_m_all_max)))
{
  lines (x = c(mean_P_ss_m_all_max[i, 'iir'], mean_P_ss_m_all_max[i, 'iir']), y = c(mean_P_ss_m_all_max[i, 'SSR_mean'] - sd_P_ss_m_all_max[i, 'SSR_mean'], c(mean_P_ss_m_all_max[i, 'SSR_mean'] + sd_P_ss_m_all_max[i, 'SSR_mean'])), col = 'darkgrey')
}
lines (SSR_mean ~ iir, mean_P_ss_m_all_max[-12,], lty = c('solid'))
lines (SSR_mean ~ iir, mean_P_ss_m_all_max[-1:-10,], lty = c('dashed'))
points (SSR_mean ~ iir, mean_P_ss_m_all_max, pch = c(22, rep (21, 11)), bg = c('black', rep ('gray', 10), 'white'))
mtext (side = 3, text = '(d)', adj = -.25, line = 2, cex = 1.1, font = 2)

plot (b1_mean ~ iir, mean_P_ss_m_all_max, xlim = c(0, 3.5), ylim = c(-.1, .1), type = 'n', las = 1, xaxt = 'n', bty = 'l', ann = F)
title (xlab = list ('Relative ITV index'), ylab = 'SSR', line = 2.5)
title (main = list ("Site-specific CWM ~ env\nRegression slope estimate", cex = 1))
axis (1, at = seq (0, 3, by = 0.5), labels = format (seq (0, 3, by = 0.5)))
axis (1, at = 3.5, labels = expression (infinity))
for (i in seq (1, nrow (mean_P_ss_m_all_max)))
{
  lines (x = c(mean_P_ss_m_all_max[i, 'iir'], mean_P_ss_m_all_max[i, 'iir']), y = c(mean_P_ss_m_all_max[i, 'b1_mean'] - sd_P_ss_m_all_max[i, 'b1_mean'], c(mean_P_ss_m_all_max[i, 'b1_mean'] + sd_P_ss_m_all_max[i, 'b1_mean'])), col = 'darkgrey')
}
lines (b1_mean ~ iir, mean_P_ss_m_all_max[-12,], lty = c('solid'))
lines (b1_mean ~ iir, mean_P_ss_m_all_max[-1:-10,], lty = c('dashed'))
points (b1_mean ~ iir, mean_P_ss_m_all_max, pch = c(22, rep (21, 11)), bg = c('black', rep ('gray', 10), 'white'))
mtext (side = 3, text = '(d)', adj = -.25, line = 2, cex = 1.1, font = 2)

plot (b1_se_mean ~ iir, mean_P_ss_m_all_max, xlim = c(0, 3.5), ylim = c(0, .2), type = 'n', las = 1, xaxt = 'n', bty = 'l', ann = F)
title (xlab = list ('Relative ITV index'), ylab = 'SSR', line = 2.5)
title (main = list ("Site-specific CWM ~ env\nError of the regression slope", cex = 1))
axis (1, at = seq (0, 3, by = 0.5), labels = format (seq (0, 3, by = 0.5)))
axis (1, at = 3.5, labels = expression (infinity))
for (i in seq (1, nrow (mean_P_ss_m_all_max)))
{
  lines (x = c(mean_P_ss_m_all_max[i, 'iir'], mean_P_ss_m_all_max[i, 'iir']), y = c(mean_P_ss_m_all_max[i, 'b1_se_mean'] - sd_P_ss_m_all_max[i, 'b1_se_mean'], c(mean_P_ss_m_all_max[i, 'b1_se_mean'] + sd_P_ss_m_all_max[i, 'b1_se_mean'])), col = 'darkgrey')
}
lines (b1_se_mean ~ iir, mean_P_ss_m_all_max[-12,], lty = c('solid'))
lines (b1_se_mean ~ iir, mean_P_ss_m_all_max[-1:-10,], lty = c('dashed'))
points (b1_se_mean ~ iir, mean_P_ss_m_all_max, pch = c(22, rep (21, 11)), bg = c('black', rep ('gray', 10), 'white'))
mtext (side = 3, text = '(d)', adj = -.25, line = 2, cex = 1.1, font = 2)

# Case study: Lalashan Forest Dynamics Plot data ----

## Import data ----

#setwd ("c:\\Users\\Zeleny\\Dropbox\\CLANKY\\CWM and intraspecific trait variation\\scripts\\")
#setwd ('c:\\Users\\zeleny\\Documents\\David clanky\\CWM ITV paper\\')
setwd ('h:\\Dropbox_select\\CLANKY\\CWM and intraspecific trait variation\\scripts\\')

com <- read.delim ("LFDP_com.txt")
traits <- read.delim ("LFDP_traits.txt")
topo <- read.delim ("LFDP_topo.txt")

inflation_ss <- function (com, traits_ss, env, alpha = 0.05, nrep = 1000) 
{
  p_null <- replicate (nrep, expr = {
    traits_ss_null <- sample_traits_ss (traits_ss, inter = T, intra = T)
    CWM_ss_null <- cwm_ss (com, traits_ss_null)
    summary (lm (CWM_ss_null ~ env))$coefficient [2, 4]
  })
  inflation <- sum (p_null < alpha) / (alpha * nrep)
  return (inflation)
}

inflation_itv <- function (com, traits_ss, env, alpha = 0.05, nrep = 1000) 
{
  p_null <- replicate (nrep, expr = {
    traits_itv_null <- sample_traits_ss (traits_ss, inter = F, intra = T)
    CWM_itv_null <- cwm_itv (com, traits_itv_null)
    summary (lm (CWM_itv_null ~ env))$coefficient [2, 4]
  })
  inflation <- sum (p_null < alpha) / (alpha * nrep)
  return (inflation)
}    

inflation_ss_max <- function (com, traits_ss, env, alpha = 0.05, nrep = 1000, perm = 199) 
{
  p_null <- replicate (nrep, expr = {
    traits_ss_rand <- sample_traits_ss (traits_ss = traits_ss, intra = TRUE, inter = TRUE)
    test_cwm_ss (com = com, traits_ss = traits_ss_rand, env = env, perm = perm)$P_max
  })
  inflation <- sum (p_null < alpha) / (alpha * nrep)
  return (inflation)
}

inflation_ss_max_parallel <- function (com, traits_ss, env, alpha = 0.05, nrep = 1000, perm = 199, nclust = 20) 
{
  cl <- makeCluster (nclust)
  clusterExport (cl, varlist = c("test_cwm_ss", "cwm_ss", "sample_traits_ss", "com", "traits_ss", "env", "perm"), env = environment ())
  p_null <- parLapply (cl, 1:nrep, fun = function (x) {
    traits_ss_rand <- sample_traits_ss (traits_ss = traits_ss, intra = TRUE, inter = TRUE)
    test_cwm_ss (com = com, traits_ss = traits_ss_rand, env = env, perm = perm)$P_max
  })
  stopCluster (cl)
  inflation <- sum (unlist (p_null) < alpha) / (alpha * nrep)
  return (inflation)
}

# Calculate CA1 of each site ----
CA1 <- scores (cca (com), choices = 1, display = "site")

# Calculate site-specific traits ----
trait_types <- factor (c("LA", "Lth", "SLA", "LDMC"), levels = c("LA", "Lth", "SLA", "LDMC"), ordered = T)

# traits_ss <- traits %>%
#   mutate (LA = log10 (LA),
#           SLA = log10 (SLA)) %>%
#   group_by (subplot, species) %>%
#   summarise_at (as.character (trait_types), mean)

# First summarise and average, then transform
traits_ss <- traits %>%
  group_by (subplot, species) %>%
  summarise_at (as.character (trait_types), mean) %>%
  mutate (LA = log10 (LA), SLA = log10 (SLA))

traits_ss_array <- traits_ss %>%
  gather (key = traits, value = value, all_of (trait_types)) %>%
  mutate (traits = factor (traits, levels = trait_types, ordered = T)) %>%
  split (.$traits) %>%
  lapply (function (x) select (x, -traits)) %>%
  lapply (function (x) spread (x, key = species, value = value)) %>%
  lapply (function (x) column_to_rownames (x, "subplot")) %>%
  abind (along = 3)

# Calculates inflation factors ----
var_ratio <- apply (traits_ss_array, 3, FUN = function (x) intra_inter_ratio (x))

set.seed (48468) # set for reproducibility
# inflation of site-specific CWM 
IF_ss_CA1 <- apply (traits_ss_array, 3, FUN = function (x) inflation_ss (com, x, CA1, nrep = 10000))
IF_ss_windwardness <- apply (traits_ss_array, 3, FUN = function (x) inflation_ss (com, x, topo$windwardness, nrep = 10000))


set.seed (62654) # set for reproducibility
# inflation of intraspecific CWM
IF_itv_CA1 <- apply (traits_ss_array, 3, FUN = function (x) inflation_itv (com = com, traits_ss = x, env = CA1, nrep = 10000))
IF_itv_windwardness <- apply (traits_ss_array, 3, FUN = function (x) inflation_itv (com, x, topo$windwardness, nrep = 10000))

set.seed (98461) # set for reproducibility
# inflation of site-specific CWM tested by newly introduced max test
IF_ss_max_CA1 <- apply (traits_ss_array, 3, FUN = function (x) inflation_ss_max_parallel (com = com, traits_ss = x, env = CA1, nrep = 10000, perm = 199, nclust = 10))

IF_ss_max_windwardness <- apply (traits_ss_array, 3, FUN = function (x) inflation_ss_max_parallel (com = com, traits_ss = x, env = topo$windwardness, nrep = 10000, perm = 199, nclust = 10)) 

# save (list = c('IF_ss_max_CA1','IF_ss_max_windwardness'), file = 'IF_ss_max.RData')
# load ('IF_ss_max.RData')

# Drawing figures with dependence of inflation factor on intraspecific trait variation ----
jpeg (filename = "inflation.factors2.jpg", width = 12, height = 12, units = 'cm', res = 600, pointsize = 8)
par (mfrow = c(2,2))
par (mar = c(5, 4, 2, 2), xpd = F)

plot (IF_ss_CA1 ~ var_ratio, type = "n", las = 1, bty = 'l', ann = F, xlim = c(0.2, 0.7), ylim = c(0, 6.5))
title (xlab = list ('Relative ITV index'), ylab = expression (Inflation~ index~(alpha==0.05)), line = 2.5)
title (main = list ("Site-specific CWM ~ env\nStandard parametric test", cex = 1, font = 1))
text (IF_ss_CA1 ~ var_ratio, labels = trait_types, cex = 0.7, font = 1, col = "black")
text (IF_ss_windwardness ~ var_ratio, labels = trait_types, cex = 0.7, font = 3, col = "grey50")
abline (h = 1, col = 'grey', lty = 'dashed')
mtext (side = 3, text = '(a)', adj = -.25, line = 0, cex = 1.1, font = 2)

plot (IF_ss_max_CA1 ~ var_ratio, type = "n", las = 1, bty = 'l', ann = F, xlim = c(0.2, 0.7), ylim = c(0, 6.5))
title (xlab = list ('Relative ITV index'), ylab = expression (Inflation~ index~(alpha==0.05)), line = 2.5)
title (main = list ("Intraspecific CWM ~ env\nStandard parametric test", cex = 1, font = 1))
text (IF_ss_max_CA1 ~ var_ratio, labels = trait_types,  cex = 0.7, font = 1, col = "black")
text (IF_ss_max_windwardness ~ var_ratio, labels = trait_types, cex = 0.7, font = 3, col = "grey50")
abline (h = 1, col = 'grey', lty = 'dashed')
mtext (side = 3, text = '(b)', adj = -.25, line = 0, cex = 1.1, font = 2)

plot (IF_itv_CA1 ~ var_ratio, type = "n", las = 1, bty = 'l', ann = F, xlim = c(0.2, 0.7), ylim = c(0, 6.5))
title (xlab = list ('Relative ITV index'), ylab = expression (Inflation~ index~(alpha==0.05)), line = 2.5)
title (main = list ("Site-specific CWM ~ env\n'Max' permutation test", cex = 1, font = 1))
text (IF_itv_CA1 ~ var_ratio, labels = trait_types, cex = 0.7, font = 1, col = "black")
text (IF_itv_windwardness ~ var_ratio, labels = trait_types, cex = 0.7, font = 3, col = "grey50")
abline (h = 1, col = 'grey', lty = 'dashed')
mtext (side = 3, text = '(c)', adj = -.25, line = 0, cex = 1.1, font = 2)

dev.off ()


# Power tests ----
#library (simcom)
#library (parallel)
#library (tidyverse)
complete_test <- function (totS, Np, m){
  sim <- simul.comm (totS = totS, min.niche.breath = 500, max.niche.breath = 5000) # 50 species
  sam <- sample.comm (sim, Np = Np) # 25 sites
  com <- sam$a.mat
  env <- sam$sample.x
  traits_ss <- gen_traits_ss (com = com, env = env, m = m, scale = TRUE, jitter = 1, amount = 500)
  traits_fixed <- colMeans (traits_ss, na.rm = T)
  CWM_SS <- cwm_ss (com = com, traits_ss = traits_ss)
  CWM_fixed <- cwm (com, traits_fixed)
  CWM_ITV <- cwm_itv (com = com, traits_ss = sample_traits_ss (traits_ss))
  TEST_ITV <- test_cwm_itv (com = com, traits_ss = traits_ss, env = env)
  TEST_SS <- test_cwm_ss (com = com, traits_ss = traits_ss, env = env)
  P <- c (TEST_SS$P_par, TEST_SS$P_max, TEST_ITV$P_par)
  return (P < 0.05)
}

res <- expand.grid (no_spe = c(10, 25, 50, 100), no_sit = c(10, 25, 50, 100), m = c(-1, -0.2, 0, 0.2, 1), repl = 1:1000)
res <- cbind (res, SS_P_par = 0, SS_P_max = 0, ITV_P_par = 0)

cl <- makeCluster(16)
clusterExport (cl, list ('complete_test', 'gen_traits_ss', 'cwm_ss', 'cwm', 'cwm_itv', 'sample_traits_ss', 'test_cwm_itv', 'test_cwm_ss'))
clusterEvalQ (cl, library (simcom))
res_apply <- parApply (cl, res, 1, FUN = function (x) try (expr = complete_test (totS = x[1], Np = x[2], m = x[3]), silent = TRUE))
stopCluster (cl)

res_apply_backup <- res_apply
res_apply_orig <- res_apply

lapply (which (sapply (res_apply_orig, FUN = function (x) !is.logical (x))), FUN = function (i) res_apply[[i]] <- as.logical (c(NA, NA, NA)))

res[,5:7] <- do.call ('rbind.data.frame',  (res_apply))
res2 <- res[res[,5] %in% c('TRUE', 'FALSE'),]

res_sum <- res %>% group_by (no_spe, no_sit, m) %>% summarise (SS_P_par = sum (as.logical (SS_P_par), na.rm = TRUE)/n (), SS_P_max = sum (as.logical (SS_P_max), na.rm = TRUE)/n(), ITV_P_par = sum (as.logical (ITV_P_par), na.rm = TRUE)/n())

p <- res_sum %>% rename (c('No. species' = 'no_spe')) %>%
  gather (key = P_type, value = P_val, SS_P_par:ITV_P_par) %>%
  ggplot (aes (x = no_sit, y = P_val, colour = P_type )) + 
  geom_line () + 
  facet_grid (rows = vars (`No. species`), cols = vars (m), as.table = FALSE, labeller = label_both) + 
  theme_bw () +
  xlab ('No. samples') + 
  ylab ('Proportion of significant results (P < 0.05)') +
  guides(colour=guide_legend(title="Test type"))

ggsave (filename = 'power_test.jpeg', plot = p, width = 6.25, height = 5, units = 'in')

# Testing pairwise significance of traits and env. variables from LFDP  ----

topo_types <- factor (colnames (topo), levels = colnames (topo), ordered = T)
CWM_types <- c("sitespecific", "fixed", "intraspecific")
test_types <- c("R_squared", "F", "parametric", "max")

row_names_pvalue <- paste (sort (rep (topo_types, 3)), rep (CWM_types, 3), sep = "_")
col_names_pvalue <- paste (sort (rep (trait_types, 4)),  rep (test_types, 4), sep = "_")

p_values <- matrix (nrow = length (row_names_pvalue), 
                    ncol = length (col_names_pvalue), 
                    dimnames = list (row_names_pvalue, col_names_pvalue))

set.seed (46465)
p_values [1, ] <- apply (traits_ss_array, 3, function (x) 
  {
  res_test_cwm_ss <- test_cwm_ss (com, x, topo$elevation, perm = 999)
  c(r2 = res_test_cwm_ss$r2, F_obs = res_test_cwm_ss$F_obs, P_par = res_test_cwm_ss$P_par, P_max = res_test_cwm_ss$P_max)
  })
p_values [2, ] <- apply (traits_ss_array, 3, function (x)
  {
  res_test_cwm_f <- test_cwm_f (com, x, topo$elevation, perm = 999)
  c(r2 = res_test_cwm_f$r2, F_obs = res_test_cwm_f$F_obs, P_par = res_test_cwm_f$P_par, P_max = res_test_cwm_f$P_max)
  })
p_values [3, ] <- apply (traits_ss_array, 3, function (x)
  {
  res_test_cwm_itv <- test_cwm_itv (com, x, topo$elevation)
  c(r2 = res_test_cwm_itv$r2, F_obs = res_test_cwm_itv$F_obs, P_par = res_test_cwm_itv$P_par, P_max = NA)
  })

p_values [4, ] <- apply (traits_ss_array, 3, function (x) 
{
  res_test_cwm_ss <- test_cwm_ss (com, x, topo$convexity, perm = 999)
  c(r2 = res_test_cwm_ss$r2, F_obs = res_test_cwm_ss$F_obs, P_par = res_test_cwm_ss$P_par, P_max = res_test_cwm_ss$P_max)
})
p_values [5, ] <- apply (traits_ss_array, 3, function (x)
{
  res_test_cwm_f <- test_cwm_f (com, x, topo$convexity, perm = 999)
  c(r2 = res_test_cwm_f$r2, F_obs = res_test_cwm_f$F_obs, P_par = res_test_cwm_f$P_par, P_max = res_test_cwm_f$P_max)
})
p_values [6, ] <- apply (traits_ss_array, 3, function (x)
{
  res_test_cwm_itv <- test_cwm_itv (com, x, topo$convexity)
  c(r2 = res_test_cwm_itv$r2, F_obs = res_test_cwm_itv$F_obs, P_par = res_test_cwm_itv$P_par, P_max = NA)
})

p_values [7, ] <- apply (traits_ss_array, 3, function (x) 
{
  res_test_cwm_ss <- test_cwm_ss (com, x, topo$windwardness, perm = 999)
  c(r2 = res_test_cwm_ss$r2, F_obs = res_test_cwm_ss$F_obs, P_par = res_test_cwm_ss$P_par, P_max = res_test_cwm_ss$P_max)
})
p_values [8, ] <- apply (traits_ss_array, 3, function (x)
{
  res_test_cwm_f <- test_cwm_f (com, x, topo$windwardness, perm = 999)
  c(r2 = res_test_cwm_f$r2, F_obs = res_test_cwm_f$F_obs, P_par = res_test_cwm_f$P_par, P_max = res_test_cwm_f$P_max)
})
p_values [9, ] <- apply (traits_ss_array, 3, function (x)
{
  res_test_cwm_itv <- test_cwm_itv (com, x, topo$windwardness)
  c(r2 = res_test_cwm_itv$r2, F_obs = res_test_cwm_itv$F_obs, P_par = res_test_cwm_itv$P_par, P_max = NA)
})

write.csv (p_values, file = "p.values.csv")

# Drawing regressions of CWM and env from LFDP data ----
plot_cwm_reg <- function (CWM_ss, CWM_f, CWM_itv, env, P_ss = NULL, P_f = NULL, P_itv = NULL, xlab = 'Env', ylab = 'CWM', plot.legend = TRUE, pos.legend = 'topright', ...)
{
  CWM_ss_scaled <- scale (CWM_ss)
  CWM_f_scaled <- scale (CWM_f)
  CWM_itv_scaled <- scale (CWM_itv)
  
  cols <- grey.colors(3, start = 0.2, end = 0.8)
  col_CWM_ss <- cols[1]
  col_CWM_f <- cols[2]
  col_CWM_itv <- cols[3]
  
  pch_CWM_ss <- 15
  pch_CWM_f <- 16
  pch_CWM_itv <- 17
  
  ylim_all <- range (c(CWM_ss_scaled, CWM_f_scaled, CWM_itv_scaled))
  
  plot (scale (CWM_ss) ~ env, type = "n", ylim = ylim_all, xlab = xlab, ylab = ylab, bty = 'l', las = 1, ...)
  
  points (CWM_ss_scaled ~ env, pch = pch_CWM_ss, cex = 0.9, col = col_CWM_ss)
  points (CWM_f_scaled ~ env, pch = pch_CWM_f, cex = 0.9, col = col_CWM_f)
  points (CWM_itv_scaled ~ env, pch = pch_CWM_itv, cex = 0.9, col = col_CWM_itv)
  
  abline (lm (CWM_ss_scaled ~ env), lty = ifelse (P_ss < 0.05, "solid", "dotted"), lwd = 1, col = col_CWM_ss)
  abline (lm (CWM_f_scaled ~ env), lty = ifelse (P_f < 0.05, "solid", "dotted"), lwd = 1, col = col_CWM_f)
  abline (lm (CWM_itv_scaled ~ env), lty = ifelse (P_itv < 0.05, "solid", "dotted"), lwd = 1, col = col_CWM_itv)  
  
  if (plot.legend) legend (pos.legend, pch = 15:17, col = cols, cex = 0.9, legend = c('ss CWM', 'f CWM', 'is CWM'), bty = 'n', lwd = 2, lty = c(ifelse (P_ss < 0.05, "solid", "dotted"), ifelse (P_f < 0.05, "solid", "dotted"), ifelse (P_itv < 0.05, "solid", "dotted")))
  
}

jpeg (filename = "CWM-env_regressions.jpg", width = 12, height = 12, units = 'cm', res = 600, pointsize = 8)
par (mfrow = c(2,2))
par (mar = c(5, 4, 2, 2), xpd = F)
plot_cwm_reg (CWM_ss = CWM_ss_LA, CWM_f = CWM_f_LA, CWM_itv = CWM_itv_LA, env = topo$windwardness, P_ss = 0.104, P_f = 0.488, P_itv = 0.039, xlab = 'Windwardness', ylab = 'CWM of log (LA)', main = list ('CWM of log (LA) ~ Windwardness', cex = 1), plot.legend = F)
mtext (side = 3, text = '(a)', adj = -.25, line = 0, cex = 1.1, font = 2)

plot_cwm_reg (CWM_ss = CWM_ss_Lth, CWM_f = CWM_f_Lth, CWM_itv = CWM_itv_Lth, env = topo$windwardness, P_ss = 0.007, P_f = 0.364, P_itv = 0.001, xlab = 'Windwardness', ylab = 'CWM of Lth', main = list ('CWM of Lth ~ Windwardness', cex = 1), plot.legend = F)
mtext (side = 3, text = '(b)', adj = -.25, line = 0, cex = 1.1, font = 2)

plot_cwm_reg (CWM_ss = CWM_ss_SLA, CWM_f = CWM_f_SLA, CWM_itv = CWM_itv_SLA, env = topo$windwardness, P_ss = 0.003, P_f = 0.311, P_itv = 0.000, xlab = 'Windwardness', ylab = 'CWM of log (SLA)', main = list ('CWM of log (SLA) ~ Windwardness', cex = 1), plot.legend = F)
mtext (side = 3, text = '(c)', adj = -.25, line = 0, cex = 1.1, font = 2)

plot.new ()
plot.window (xlim = 0:1, ylim = 0:1)
legend ('topleft', pch = c(15:17, NA, NA, NA), col = c(grey.colors(3, start = 0.2, end = 0.8), 'black', 'black', 'black'), cex = 1.2, legend = c('site-specific CWM', 'fixed CWM', 'intraspecific CWM', '', 'significant (P < 0.05)', 'not significant'), bty = 'n', lwd = 1, lty = c('solid', 'solid', 'solid', 'blank', 'solid', 'dotted'))

dev.off ()

# Variation partitioning into effect of turnover, ITV and covariation (sensu Leps et al. 2011) ----
relvar_total <- function (com, traits_ss) {
  CWM_ss <- cwm_ss (com, traits_ss)
  CWM_f <- cwm_f (com, traits_ss)
  CWM_itv <- cwm_itv (com, traits_ss)
  
  null_ss <- sum (resid (lm (CWM_ss ~ 1))^2)
  null_f <- sum (resid (lm (CWM_f ~ 1))^2)
  null_itv <- sum (resid (lm (CWM_itv ~ 1))^2)
  null_cov <- null_ss - null_f - null_itv
  
  rel_null_f <- (null_f / null_ss) * 100
  rel_null_itv <- (null_itv / null_ss) * 100
  rel_null_cov <- (null_cov / null_ss) * 100
  return (c(rel_null_f, rel_null_itv, rel_null_cov))
} 

relvar_part <- function (com, traits_ss, env) {
  CWM_ss <- cwm_ss (com, traits_ss)
  CWM_f <- cwm_f (com, traits_ss)
  CWM_itv <- cwm_itv (com, traits_ss)
  
  null_ss <- sum (resid (lm (CWM_ss ~ 1))^2)
  ESS_ss <- sum ((predict (lm (CWM_ss ~ env)) - mean (CWM_ss))^2)
  ESS_f <- sum ((predict (lm (CWM_f ~ env)) - mean (CWM_f))^2)
  ESS_itv <- sum ((predict (lm (CWM_itv ~ env)) - mean (CWM_itv))^2)
  cov <- ESS_ss - ESS_f - ESS_itv
  
  rel_f <- (ESS_f / null_ss) * 100
  rel_itv <- (ESS_itv / null_ss) * 100
  rel_cov <- (cov / null_ss) * 100
  return (c(rel_f, rel_itv, rel_cov))
}


env_types <- c("total", colnames (topo))

var_resource <- c("Fixed", "ITV", "Cov")

rel_LA <- matrix (nrow = length (env_types), ncol = length (var_resource),
                  dimnames = list (env_types, var_resource))

rel_LA [1, ] <- relvar_total (com, traits_ss_array [, , "LA"])
rel_LA [2:4, ] <- t (apply (topo, 2, function (x) 
  relvar_part (com, traits_ss_array [, , "LA"], x)))
rel_LA <- rel_LA %>%
  as.data.frame () %>%
  rownames_to_column ("env") %>%
  mutate (trait = "LA") %>%
  gather (key = variation, value = value, all_of (as.character (var_resource)))

rel_Lth <- matrix (nrow = length (env_types), ncol = length (var_resource),
                   dimnames = list (env_types, var_resource))

rel_Lth [1, ] <- relvar_total (com, traits_ss_array [, , "Lth"])
rel_Lth [2:4, ] <- t (apply (topo, 2, function (x) 
  relvar_part (com, traits_ss_array [, , "Lth"], x)))
rel_Lth <- rel_Lth %>%
  as.data.frame () %>%
  rownames_to_column ("env") %>%
  mutate (trait = "Lth") %>%
  gather (key = variation, value = value, all_of (as.character (var_resource)))

rel_SLA <- matrix (nrow = length (env_types), ncol = length (var_resource),
                   dimnames = list (env_types, var_resource))

rel_SLA [1, ] <- relvar_total (com, traits_ss_array [, , "SLA"])
rel_SLA [2:4, ] <- t (apply (topo, 2, function (x) 
  relvar_part (com, traits_ss_array [, , "SLA"], x)))
rel_SLA <- rel_SLA %>%
  as.data.frame () %>%
  rownames_to_column ("env") %>%
  mutate (trait = "SLA") %>%
  gather (key = variation, value = value, all_of (as.character (var_resource)))


rel_LDMC <- matrix (nrow = length (env_types), ncol = length (var_resource),
                    dimnames = list (env_types, var_resource))

rel_LDMC [1, ] <- relvar_total (com, traits_ss_array [, , "LDMC"])
rel_LDMC [2:4, ] <- t (apply (topo, 2, function (x) 
  relvar_part (com, traits_ss_array [, , "LDMC"], x)))
rel_LDMC <- rel_LDMC %>%
  as.data.frame () %>%
  rownames_to_column ("env") %>%
  mutate (trait = "LDMC") %>%
  gather (key = variation, value = value, all_of (as.character (var_resource)))

plot_table <- rbind (rel_LA, rel_Lth, rel_SLA, rel_LDMC) %>%
  mutate (env = factor (env, levels = rev (env_types), ordered = T),
          trait = factor (trait, levels = trait_types, ordered = T),
          variation = factor (variation, levels = var_resource, ordered = T)) 

library (lattice)
jpeg (filename = "variation.partitioning.jpg", width = 8, height = 6,
      units = "in", res = 600, quality = 100)
barchart (value ~ variation | trait + env, data = plot_table,
          origin = 0, strip = F,
          col = c("grey55", "grey75", "white"),
          xlab = "Source of trait variations",
          ylab = "Relative variation (%)",
          xlab.top = c("log (LA)", "Lth", "log (SLA)", "LDMC"),
          ylab.right = c("Windwardness", "Convexity", "Elevation", "Total"),
          scales = list (alternating = 
                           rep (1, length (trait_types) * length (env_types)),
                         tck = c(0.8, 0)),
          key = list (
            space = "right",
            points = list (col = rep ("transparent", 5)),
            rectangles = list (size = rep (1.5, 5),
                               height = rep (0.6, 5),
                               alpha = c(1, 0, 1, 0, 1),
                               col = c("grey55", "transparent","grey75",
                                       "transparent", "white")),
            text = list (c("Fixed", "", "ITV", "", "Covariation"),
                         cex = c(1, 0.1, 1, 0.1, 1))
          ))
dev.off ()

