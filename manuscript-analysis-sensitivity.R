

library(tidyverse)
library(multcomp)

my_fn <- function(data, mapping, method="loess", ...){
  p <- ggplot(data = data, mapping = mapping) +
    geom_point() +
    geom_smooth(method=method, ...)
  p
}
stargazer2 <- function(model, odd.ratio = T, ...) {
  if(!("list" %in% class(model))) model <- list(model)
  if (odd.ratio) {
    coefOR2 <- lapply(model, function(x) exp(coef(x)))
    seOR2 <- lapply(model, function(x) exp(coef(x)) * summary(x)$coef[, 2])
    p2 <- lapply(model, function(x) summary(x)$coefficients[, 4])
    stargazer::stargazer(model, coef = coefOR2, se = seOR2, p = p2, ...)
  } else {
    stargazer::stargazer(model, ...)
  }
}
CalcZoneLength <- function(x_logical){
  # Calculates the length of each zone exposure. 
  # @param x_logical A logical variable where TRUE indicates in-zone.
  if(!is.logical(x_logical)) stop("The input should be a logical variable.")
  xrle <- rle(x_logical)
  if(all(xrle$values == FALSE)){
    res <- data.frame(
      episode = 1,
      zlength = 0
    )
  } else {
    zone.length <- xrle$length[xrle$values]
    res <- data.frame(
      episode = 1:length(zone.length),
      zlength = zone.length
    )
  }
  return(res)
}

datapath <- "~/path-to-data/"
dcovariates <- read_csv(paste0(datapath, "analytic/covariates_post_exclusion.csv")) %>%
  mutate(across(starts_with("cat"), as.factor)) %>%
  filter(cat_surgcat4 == "1_CABG") %>%
  dplyr::select(
    id, bin_aki48h, #bin_akicomposite, cat_aki7d, 
    val_age, cat_gender, val_bmi, bin_hypertn, bin_dm,    
    bin_copd, bin_chf, bin_priormi, bin_cva, bin_pvd,       
    bin_emergent, bin_redo, val_hct, val_hdef, bin_statin,  
    bin_acearb, bin_betablocker, bin_iabp, val_perfustm, val_visscore,   
    val_creatlst, cat_transbc, val_predmort, val_crystalloid, val_proctime
  )
dat <- read_rds(paste0(datapath, "analytic/hemo_timeseries_interp_post_exclusion.rds")) %>%
  dplyr::select(id, time, cat_cpb, val_MAP, val_CVP) %>%
  rename(map = val_MAP, cvp = val_CVP) %>%
  filter(!(is.na(map) | is.na(cvp))) %>%
  mutate(
    hypotension = map < 65,
    cond1 = map <= 55 & (cvp >=0 & cvp <= 20), # 
    cond2 = (map > 55 & map <= 65) & (cvp > 8 & cvp <= 20), # 
    cond3 = (map > 65 & map <= 115) & (cvp > 12 & cvp <= 20), # don't aggresively make moves on
    zone45 = cond1 | cond2 | cond3,
    cond4 = (map > 55 & map <= 75) & (cvp > 0 & cvp <= 8), 
    cond5 = (map > 65 & map <= 85) & (cvp > 8 & cvp <= 10), 
    cond6 = (map > 65 & map <= 115) & (cvp > 10 & cvp <= 12), 
    zone3 = cond4 | cond5 | cond6,
    cond7 = (map > 75 & map <= 115) & (cvp > 0 & cvp <= 8), 
    cond8 = (map > 85 & map <= 115) & (cvp > 8 & cvp <= 10), 
    zone12 = cond7 | cond8
  ) %>%
  rename(phase = cat_cpb)
id.no.phase <- unique(filter(dat, is.na(phase))$id) # remove 5 subjects, only 3 of them are CABG
dat <- dat %>% filter(!(id %in% id.no.phase))

depisodes.45 <- dat %>% 
  group_by(id, phase) %>%
  reframe(CalcZoneLength(zone45)) %>%
  arrange(id)
depisodes.summary.45 <- depisodes.45 %>% 
  group_by(id, phase) %>%
  summarise(
    ztotal_45 = sum(zlength)
  ) %>%
  arrange(desc(ztotal_45)) %>%
  ungroup() %>%
  arrange(id) %>%
  pivot_wider(names_from = "phase", values_from = c("ztotal_45"), names_prefix = "z45_")
depisodes.3 <- dat %>% 
  group_by(id, phase) %>%
  reframe(CalcZoneLength(zone3)) %>%
  arrange(id)
depisodes.summary.3 <- depisodes.3 %>% 
  group_by(id, phase) %>%
  summarise(
    # nepisodes = n(),
    ztotal_3 = sum(zlength)
  ) %>%
  arrange(desc(ztotal_3)) %>%
  ungroup() %>%
  arrange(id) %>%
  pivot_wider(names_from = "phase", values_from = c("ztotal_3"), names_prefix = "z3_")
depisodes.12 <- dat %>% 
  group_by(id, phase) %>%
  reframe(CalcZoneLength(zone12)) %>%
  arrange(id)
depisodes.summary.12 <- depisodes.12 %>% 
  group_by(id, phase) %>%
  summarise(
    # nepisodes = n(),
    ztotal_12 = sum(zlength)
  ) %>%
  arrange(desc(ztotal_12)) %>%
  ungroup() %>%
  arrange(id) %>%
  pivot_wider(names_from = "phase", values_from = c("ztotal_12"), names_prefix = "z12_")
depisodes.hypotension <- dat %>% 
  group_by(id, phase) %>%
  reframe(CalcZoneLength(hypotension)) %>%
  arrange(id)
depisodes.summary.hypotension <- depisodes.hypotension %>% 
  group_by(id, phase) %>%
  summarise(
    ztotal_hypo = sum(zlength)
  ) %>%
  arrange(desc(ztotal_hypo)) %>%
  ungroup() %>%
  arrange(id) %>%
  pivot_wider(names_from = "phase", values_from = c("ztotal_hypo"), names_prefix = "zhypo_")
depisodes.summary <- depisodes.summary.12 %>% 
  left_join(depisodes.summary.3) %>% 
  left_join(depisodes.summary.45) %>% 
  left_join(depisodes.summary.hypotension) %>% 
  inner_join(dcovariates) %>%
  mutate(
    zhypo = zhypo_pre + zhypo_intra + zhypo_post,
    z12 = z12_pre + z12_intra + z12_post,
    z3 = z3_pre + z3_intra + z3_post,
    z45 = z45_pre + z45_intra + z45_post
  ) %>% 
  mutate(across(starts_with("cat"), as.factor)) %>%
  mutate(across(starts_with("val"), function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T)))



# Model Sensitivity ---------------------------------------------
m.phase.3zone.nocpb.nototal <- glm(
  as.integer(bin_aki48h) ~
    I(z12_pre/5) + I(z3_pre/5) + I(z45_pre/5) +
    I(z12_intra/5) + I(z3_intra/5) + I(z45_intra/5) +
    I(z12_post/5) + I(z3_post/5) + I(z45_post/5) +
    val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef +
    bin_statin + bin_acearb + bin_betablocker + bin_iabp +
    val_creatlst + val_predmort + cat_transbc + val_visscore + val_crystalloid,
  data = depisodes.summary,
  family = "binomial"
)
summary(m.phase.3zone.nocpb.nototal)


# Model with zone and hypo -------------------------------------
m.phase.hypo.phase.3zone.nocpb.nototal <- glm(
  as.integer(bin_aki48h) ~ 
    I(zhypo_pre/5) + I(z12_pre/5) + I(z3_pre/5) + I(z45_pre/5) + 
    I(zhypo_intra/5) + I(z12_intra/5) + I(z3_intra/5) + I(z45_intra/5) + 
    I(zhypo_post/5) + I(z12_post/5) + I(z3_post/5) + I(z45_post/5) + 
    val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + 
    bin_statin + bin_acearb + bin_betablocker + bin_iabp + 
    val_creatlst + val_predmort + cat_transbc + val_visscore + val_crystalloid, 
  data = depisodes.summary, 
  family = "binomial"
)
summary(m.phase.hypo.phase.3zone.nocpb.nototal)
coef.m <- exp(coef(m.phase.hypo.phase.3zone.nocpb.nototal))
ci.m <- exp(confint(m.phase.hypo.phase.3zone.nocpb.nototal))
p.val <- cut(
  summary(m.phase.hypo.phase.3zone.nocpb.nototal)$coefficients[, "Pr(>|z|)"],
  breaks = c(0, 0.001, 0.01, 0.05, 1),
  labels = c("***", "**", "*", "")
)
var.names <- c(
  "Intercept", 
  "Hypo Pre-CPB", "Zone 1 Pre-CPB", "Zone 2 Pre-CPB", "Zone 3 Pre-CPB",
  "Hypo Intra-CPB", "Zone 1 Intra-CPB", "Zone 2 Intra-CPB", "Zone 3 Intra-CPB",
  "Hypo Post-CPB", "Zone 1 Post-CPB", "Zone 2 Post-CPB", "Zone 3 Post-CPB",
  "Age", "Gender:Male", "BMI", "Hypertension", "Diabetes",
  "COPD", "CHF", "Prior MI", "CVA", "PVD", "Emergent", "Redo",
  "HCT", "HDEF", "Statin", "Acearb", "Betablocker", "IADP",
  "Baseline Creatinine", "Predicted Mortality", 
  "TrasBC = 0", "TransBC = 1", "VIS Score", "Crystalloid"
)
tab.m <- data.frame(
    Variable = var.names,
    OR = round(coef.m, 3),
    CI = paste0(
      "(", round(ci.m[, 1], 3), ", ", round(ci.m[, 2], 3), ")", 
      p.val
    )
)
row.names(tab.m) <- NULL



varnames <- c(
  "I(zhypo_pre/5)", "I(z12_pre/5)", "I(z3_pre/5)", "I(z45_pre/5)", 
  "I(zhypo_intra/5)", "I(z12_intra/5)", "I(z3_intra/5)", "I(z45_intra/5)", 
  "I(zhypo_post/5)", "I(z12_post/5)", "I(z3_post/5)", "I(z45_post/5)"
)
var.levels <- c(
  "Hypo Pre", "Zone 1 Pre", "Zone 2 Pre", "Zone 3 Pre",
  "Hypo Intra", "Zone 1 Intra", "Zone 2 Intra", "Zone 3 Intra",
  "Hypo Post", "Zone 1 Post", "Zone 2 Post", "Zone 3 Post"
)
var.labels <- c(
  "Hypo \n Pre", "Zone 1 \n  Pre", "Zone 2 \n  Pre", "Zone 3 \n  Pre",
  "Hypo  \n Intra", "Zone 1  \n Intra", "Zone 2  \n Intra", "Zone 3  \n Intra",
  "Hypo \n Post", "Zone 1  \n Post", "Zone 2  \n Post", "Zone 3  \n Post"
)
dres.m <- data.frame(
  var = factor(var.levels, levels = var.levels, labels = var.labels),
  or = exp(coef(m.phase.hypo.phase.3zone.nocpb.nototal)[varnames]),
  lower = ci.m[varnames, 1],
  upper = ci.m[varnames, 2],
  Exposure = rep(c("Hypotension", "Zone 1", "Zone 2", "Zone 3"), times = 3)
) 
dres.m <- data.frame(
  var = factor(
    c(var.levels[-c(1,5,9)], var.levels), levels = var.levels, 
    labels = var.labels
  ),
  or = c(
    exp(coef(m.phase.3zone.nocpb.nototal)[varnames0]),
    exp(coef(m.phase.hypo.phase.3zone.nocpb.nototal)[varnames])
  ),
  lower = c(ci.m0[varnames0, 1], ci.m[varnames, 1]),
  upper = c(ci.m0[varnames0, 2], ci.m[varnames, 2]),
  Exposure = c(
    rep(c("Zone 1", "Zone 2", "Zone 3"), times = 3),
    rep(c("Hypotension", "Zone 1", "Zone 2", "Zone 3"), times = 3)
  ),
  model = c(rep("No Hypo", 9), rep("With Hypo", 12))
) 
custom_colors <- c("Hypotension" = "black", "Zone 1" = "blue", "Zone 2" = "darkgrey", "Zone 3" = "red")
ggplot(dres.m) + 
  geom_point(
    aes(x = var, y = or, color = Exposure, shape = model),
    position = position_dodge(width = 0.3)
  ) + 
  geom_linerange(aes(
    x = var, y = or, ymin = lower, ymax = upper, color = Exposure, linetype = model
  ), position = position_dodge(width = 0.3)) + 
  geom_hline(yintercept = 1) + 
  scale_color_manual(values = custom_colors) + 
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(x = "", y = "Odds Ratios (Per 5 Minutes)")

stargazer2(
  list(m.phase.hypo.phase.3zone.nocpb.nototal),
  digits = 3, type = "text", star.cutoffs = c(0.05, 0.01, 0.001)
)
GGally::ggpairs(
  depisodes.summary %>% dplyr::select(
    zhypo_intra, z45_intra, z3_intra, z12_intra
  ),
  lower = list(continuous = GGally::wrap(my_fn, method="lm"))
)


# Plot AKI vs Intra zone minutes by levels of hypo
q.hypo.intra <- quantile(
  depisodes.summary$zhypo_intra, probs = c(0, 0.33, 0.66, 1)
)
    
depisodes.summary <- depisodes.summary %>%
  mutate(
    cat_hypo_intra = cut(
      x = zhypo_intra, 
      labels = c("hypo intra low", "hypo intra mid", "hypo intra high"),
      breaks = q.hypo.intra, include.lowest = T, 
      right = T, ordered_result = T
    )
  )


# filter(depisodes.summary, z12_intra < 75) %>%
depisodes.summary %>%
ggplot() +
  geom_point(
    aes(y = bin_aki48h, x = z12_intra), 
    position = position_jitter(height = 0.07)
  ) + 
  geom_smooth(
    aes(y = bin_aki48h, x = z12_intra), 
    method = "glm", method.args = list(family = "binomial")
  ) + 
  facet_wrap(~ cat_hypo_intra, scales = "free") +
  theme_bw() + labs(x = "Zone 1 Intra")

depisodes.summary %>%
  ggplot() +
  geom_point(
    aes(y = bin_aki48h, x = z3_intra), 
    position = position_jitter(height = 0.07)
  ) + 
  geom_smooth(
    aes(y = bin_aki48h, x = z3_intra), 
    method = "glm", method.args = list(family = "binomial")
  ) + 
  facet_wrap(~ cat_hypo_intra, scales = "free") +
  theme_bw()  + labs(x = "Zone 2 Intra")


depisodes.summary %>%
  ggplot() +
  geom_point(
    aes(y = bin_aki48h, x = z45_intra), 
    position = position_jitter(height = 0.07)
  ) + 
  geom_smooth(
    aes(y = bin_aki48h, x = z45_intra), 
    method = "glm", method.args = list(family = "binomial")
  ) + 
  facet_wrap(~ cat_hypo_intra, scales = "free") +
  theme_bw()  + labs(x = "Zone 3 Intra")




gam.phase.hypo.phase.3zone.nocpb.nototal <- gam(
  as.integer(bin_aki48h) ~ 
    I(zhypo_pre/5) + I(zhypo_intra/5) + I(zhypo_post/5) + 
    I(z12_pre/5) + I(z3_pre/5) + I(z45_pre/5) + 
    I(z12_intra/5) + I(z3_intra/5) + s(I(z45_intra/5)) + 
    I(z12_post/5) + I(z3_post/5) + I(z45_post/5) + 
    val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + 
    bin_statin + bin_acearb + bin_betablocker + bin_iabp + 
    val_creatlst + val_predmort + cat_transbc + val_visscore + val_crystalloid, 
  data = depisodes.summary, 
  family = "binomial"
)

plot(gam.phase.hypo.phase.3zone.nocpb.nototal)






# Bootstrap ----------------------------------------------------
set.seed(123)          # for reproducibility
B <- 2000              # number of bootstrap samples
N <- nrow(depisodes.summary)
m1.boot <- data.frame() # m.phase.hypo.phase.3zone
m2.boot <- data.frame() # m.phase.hypo.phase.3zone.nocpb.nototal

for (b in 1:B) {
  # 1. Resample indices with replacement
  idx <- sample(seq_len(N), size = N, replace = TRUE)
  boot_data <- depisodes.summary[idx, ]
  
  # 2. Fit logistic regression on bootstrap sample
  m.phase.hypo.phase.3zone.boot <- glm(
    as.integer(bin_aki48h) ~ 
      I(zhypo_pre/5) + I(zhypo_intra/5) + I(zhypo_post/5) + 
      I(z12_pre/5) + I(z3_pre/5) + I(z45_pre/5) + 
      I(z12_intra/5) + I(z3_intra/5) + I(z45_intra/5) + 
      I(z12_post/5) + I(z3_post/5) + I(z45_post/5) + 
      val_perfustm + val_proctime + 
      val_age + cat_gender + val_bmi +
      bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
      bin_emergent + bin_redo + val_hct + val_hdef + 
      bin_statin + bin_acearb + bin_betablocker + bin_iabp + 
      val_creatlst + val_predmort + cat_transbc + val_visscore + val_crystalloid, 
    data = boot_data, 
    family = "binomial"
  )
  m.phase.hypo.phase.3zone.nocpb.nototal.boot <- glm(
    as.integer(bin_aki48h) ~ 
      I(zhypo_pre/5) + I(zhypo_intra/5) + I(zhypo_post/5) + 
      I(z12_pre/5) + I(z3_pre/5) + I(z45_pre/5) + 
      I(z12_intra/5) + I(z3_intra/5) + I(z45_intra/5) + 
      I(z12_post/5) + I(z3_post/5) + I(z45_post/5) + 
      val_age + cat_gender + val_bmi +
      bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
      bin_emergent + bin_redo + val_hct + val_hdef + 
      bin_statin + bin_acearb + bin_betablocker + bin_iabp + 
      val_creatlst + val_predmort + cat_transbc + val_visscore + val_crystalloid, 
    data = boot_data, 
    family = "binomial"
  )
  
  # 3. Extract coefficients and their summary info
  m1.boot.summary <- summary(m.phase.hypo.phase.3zone.boot)
  m2.boot.summary <- summary(m.phase.hypo.phase.3zone.nocpb.nototal.boot)
  m1.boot.coefs <- m1.boot.summary$coefficients
  m2.boot.coefs <- m2.boot.summary$coefficients
  m1.boot <- rbind(
    m1.boot,
    data.frame(
      b = b,
      var = rownames(m1.boot.coefs),
      est = m1.boot.coefs[, "Estimate"],
      pval = m1.boot.coefs[, "Pr(>|z|)"]
    )
  )
  m2.boot <- rbind(
    m2.boot,
    data.frame(
      b = b,
      var = rownames(m2.boot.coefs),
      est = m2.boot.coefs[, "Estimate"],
      pval = m2.boot.coefs[, "Pr(>|z|)"]
    )
  )
}

# Calculate percentile 95% CI
m1.boot.percentile <- m1.boot %>%
  group_by(var) %>%
  summarise(
    lower = quantile(est, probs = 0.02500),
    upper = quantile(est, probs = 0.97500)
  ) %>%
  mutate(
    lower = exp(lower), upper = exp(upper)
  )
m2.boot.percentile <- m2.boot %>%
  group_by(var) %>%
  summarise(
    lower = quantile(exp(est), probs = 0.02500),
    upper = quantile(exp(est), probs = 0.97500)
  ) 
  # mutate(
  #   lower = exp(lower), upper = exp(upper)
  # )

exp(mean(m2.boot$est[m2.boot$var == "I(z45_pre/5)"]))
exp(mean(m2.boot$est[m2.boot$var == "I(z3_pre/5)"]))
exp(mean(m2.boot$est[m2.boot$var == "I(z45_post/5)"]))
exp(mean(m2.boot$est[m2.boot$var == "I(zhypo_intra/5)"]))


exp(quantile(m2.boot$est[m2.boot$var == "I(z45_pre/5)"], probs = c(0.025, 0.975)))
exp(quantile(m2.boot$est[m2.boot$var == "I(z3_pre/5)"], probs = c(0.025, 0.975)))
exp(quantile(m2.boot$est[m2.boot$var == "I(z45_post/5)"], probs = c(0.025, 0.975)))
exp(quantile(m2.boot$est[m2.boot$var == "I(zhypo_intra/5)"], probs = c(0.025, 0.975)))


# Get full data estimate
m1.summary <- summary(m.phase.hypo.phase.3zone)
m2.summary <- summary(m.phase.hypo.phase.3zone.nocpb.nototal)
m1.coefs <- m1.summary$coefficients
m2.coefs <- m2.summary$coefficients
m1.df <- data.frame(
  var = rownames(m1.coefs),
  est = exp(m1.coefs[, "Estimate"]),
  pval = m1.coefs[, "Pr(>|z|)"]
)
m2.df <- data.frame(
  var = rownames(m2.coefs),
  est = exp(m2.coefs[, "Estimate"]),
  pval = m2.coefs[, "Pr(>|z|)"]
)

exp(confint(m.phase.hypo.phase.3zone.nocpb.nototal))

var.plot <- c(
  "I(zhypo_intra/5)", "I(z3_pre/5)", "I(z45_pre/5)"
)
m1.boot %>%
  filter(var %in% var.plot) %>% 
  mutate(var = factor(var, levels = var.plot)) %>%
  mutate(est = exp(est)) %>%
  ggplot() +
  geom_density(aes(x = est), alpha = 0.6) +
  geom_vline(
    data = m1.df %>% filter(var %in% var.plot) %>% 
      mutate(var = factor(var, levels = var.plot)),
    aes(xintercept = est),
    color = "red", 
    linetype = "dashed", 
    size = 1
  ) +
  geom_rect(
    data = m1.boot.percentile %>% filter(var %in% var.plot) %>% 
      mutate(var = factor(var, levels = var.plot)),
    aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
    color = "grey", alpha = 0.1
  ) + 
  facet_wrap(~ var, scales = "free", nrow = 3) +
  labs(title = "Bootstrap Distribution of Odds Ratios",
       x = "Odds Ratio",
       y = "Density") +
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 0.01) # Dynamic breaks for each facet
  ) + 
  theme_bw()
m2.boot %>%
  filter(var %in% var.plot) %>% 
  mutate(var = factor(var, levels = var.plot)) %>%
  mutate(est = exp(est)) %>%
  ggplot() +
  geom_density(aes(x = est), alpha = 0.6) +
  geom_vline(
    data = m2.df %>% filter(var %in% var.plot) %>% 
      mutate(var = factor(var, levels = var.plot)),
    aes(xintercept = est),
    color = "red", 
    linetype = "dashed", 
    size = 1
  ) +
  geom_rect(
    data = m2.boot.percentile %>% filter(var %in% var.plot) %>% 
      mutate(var = factor(var, levels = var.plot)),
    aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
    color = "grey", alpha = 0.1
  ) + 
  facet_wrap(~ var, scales = "free", nrow = 3) +
  labs(title = "Bootstrap Distribution of Odds Ratios (No CPB or Surgery Time)",
       x = "Odds Ratio",
       y = "Density") +
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 0.01) # Dynamic breaks for each facet
  ) + 
  theme_bw()








