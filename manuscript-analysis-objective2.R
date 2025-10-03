
library(tidyverse)
library(multcomp)

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
    id, bin_aki48h, 
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
id.no.phase <- unique(filter(dat, is.na(phase))$id) 
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
    ztotal_12 = sum(zlength)
  ) %>%
  arrange(desc(ztotal_12)) %>%
  ungroup() %>%
  arrange(id) %>%
  pivot_wider(names_from = "phase", values_from = c("ztotal_12"), names_prefix = "z12_")

depisodes.summary <- depisodes.summary.12 %>% 
  left_join(depisodes.summary.3) %>% 
  left_join(depisodes.summary.45) %>% 
  inner_join(dcovariates)





# Pre-CPB --------------
mglm.pre <- glm(
  as.integer(bin_aki48h) ~ I(z12_pre/5) + I(z3_pre/5) + I(z45_pre/5) + 
    val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + bin_statin +
    bin_acearb + bin_betablocker + bin_iabp + val_perfustm  +
    val_creatlst + val_predmort + val_proctime,
  data = depisodes.summary,
  family = "binomial"
)
summary(mglm.pre)
coef.pre <- coef(mglm.pre)
ci.pre <- confint(mglm.pre)
dplot.pre <- data.frame(
  phase = "Pre-CPB",
  Zone = c("Zone 1", "Zone 2", "Zone 3"),
  coef = exp(coef.pre[c("I(z12_pre/5)", "I(z3_pre/5)", "I(z45_pre/5)")]),
  ci_lower = exp(ci.pre[c("I(z12_pre/5)", "I(z3_pre/5)", "I(z45_pre/5)"), 1]),
  ci_upper = exp(ci.pre[c("I(z12_pre/5)", "I(z3_pre/5)", "I(z45_pre/5)"), 2])
)
paste0(
  round(exp(coef.pre["I(z12_pre/5)"]), 3), 
  " (", round(exp(ci.pre["I(z12_pre/5)", 1]), 3), ", ", 
  round(exp(ci.pre["I(z12_pre/5)", 2]), 3), ")"
)
paste0(
  round(exp(coef.pre["I(z3_pre/5)"]), 3), 
  " (", round(exp(ci.pre["I(z3_pre/5)", 1]), 3), ", ", 
  round(exp(ci.pre["I(z3_pre/5)", 2]), 3), ")"
)
paste0(
  round(exp(coef.pre["I(z45_pre/5)"]), 3), 
  " (", round(exp(ci.pre["I(z45_pre/5)", 1]), 3), ", ", 
  round(exp(ci.pre["I(z45_pre/5)", 2]), 3), ")"
)

k_3_12 <- matrix(c(0, -1, 1, 0, rep(0, length(coef(mglm.pre)) - 4)), 1)
t_3_12 <- glht(mglm.pre, linfct = k_3_12)
summary(t_3_12)

k_45_3 <- matrix(c(0, 0, -1, 1, rep(0, length(coef(mglm.pre)) - 4)), 1)
t_45_3 <- glht(mglm.pre, linfct = k_45_3)
summary(t_45_3)

k_45_12 <- matrix(c(0, -1, 0, 1, rep(0, length(coef(mglm.pre)) - 4)), 1)
t_45_12 <- glht(mglm.pre, linfct = k_45_12)
summary(t_45_12)
confint(t_45_12)

dtest.pre <- data.frame(
  phase = "Pre-CPB",
  y = c(1.08, 1.11, 1.14),
  xmin = c("Zone 1", "Zone 2", "Zone 1"),
  xmax = c("Zone 2", "Zone 3", "Zone 3"),
  or = c(
    summary(t_3_12)$test$coefficients, 
    summary(t_45_3)$test$coefficients, 
    summary(t_45_12)$test$coefficients
  ),
  pval = p.adjust(c(
    summary(t_3_12)$test$pvalues, 
    summary(t_45_3)$test$pvalues, 
    summary(t_45_12)$test$pvalues
  ), method = "bonferroni")
) %>% mutate(
  or_pval = paste0(
    round(exp(or), 3), 
    " (p ", 
    ifelse(
      pval < 0.001, "< 0.001)", paste0("= ", round(pval, 3), ")")
    )
  )
)




# Intra-CPB --------------
mglm.intra <- glm(
  as.integer(bin_aki48h) ~ I(z12_intra/5) + I(z3_intra/5) + I(z45_intra/5) + 
    val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + bin_statin +
    bin_acearb + bin_betablocker + bin_iabp +
    val_creatlst + val_predmort + val_proctime, # + val_perfustm,
  data = depisodes.summary,
  family = "binomial"
)
summary(mglm.intra)
ci.intra <- confint(mglm.intra)
coef.intra <- coef(mglm.intra)
dplot.intra <- data.frame(
  phase = "Intra-CPB",
  Zone = c("Zone 1", "Zone 2", "Zone 3"),
  coef = exp(coef.intra[c("I(z12_intra/5)", "I(z3_intra/5)", "I(z45_intra/5)")]),
  ci_lower = exp(ci.intra[c("I(z12_intra/5)", "I(z3_intra/5)", "I(z45_intra/5)"), 1]),
  ci_upper = exp(ci.intra[c("I(z12_intra/5)", "I(z3_intra/5)", "I(z45_intra/5)"), 2])
)
paste0(
  round(exp(coef.intra["I(z12_intra/5)"]), 3), 
  " (", round(exp(ci.intra["I(z12_intra/5)", 1]), 3), ", ", 
  round(exp(ci.intra["I(z12_intra/5)", 2]), 3), ")"
)
paste0(
  round(exp(coef.intra["I(z3_intra/5)"]), 3), 
  " (", round(exp(ci.intra["I(z3_intra/5)", 1]), 3), ", ", 
  round(exp(ci.intra["I(z3_intra/5)", 2]), 3), ")"
)
paste0(
  round(exp(coef.intra["I(z45_intra/5)"]), 3), 
  " (", round(exp(ci.intra["I(z45_intra/5)", 1]), 3), ", ", 
  round(exp(ci.intra["I(z45_intra/5)", 2]), 3), ")"
)

k_3_12 <- matrix(c(0, -1, 1, 0, rep(0, length(coef(mglm.intra)) - 4)), 1)
t_3_12 <- glht(mglm.intra, linfct = k_3_12)
summary(t_3_12)

k_45_3 <- matrix(c(0, 0, -1, 1, rep(0, length(coef(mglm.intra)) - 4)), 1)
t_45_3 <- glht(mglm.intra, linfct = k_45_3)
summary(t_45_3)

k_45_12 <- matrix(c(0, -1, 0, 1, rep(0, length(coef(mglm.intra)) - 4)), 1)
t_45_12 <- glht(mglm.intra, linfct = k_45_12)
summary(t_45_12)

dtest.intra <- data.frame(
  phase = "Intra-CPB",
  y = c(1.08, 1.11, 1.14),
  xmin = c("Zone 1", "Zone 2", "Zone 1"),
  xmax = c("Zone 2", "Zone 3", "Zone 3"),
  or = c(
    summary(t_3_12)$test$coefficients, 
    summary(t_45_3)$test$coefficients, 
    summary(t_45_12)$test$coefficients
  ),
  pval = p.adjust(c(
    summary(t_3_12)$test$pvalues, 
    summary(t_45_3)$test$pvalues, 
    summary(t_45_12)$test$pvalues
  ), method = "bonferroni")
) %>% mutate(
  or_pval = paste0(
    round(exp(or), 3), 
    " (p ", 
    ifelse(
      pval < 0.001, "< 0.001)", paste0("= ", round(pval, 3), ")")
    )
  )
)




# Post-CPB --------------
mglm.post <- glm(
  as.integer(bin_aki48h) ~ I(z12_post/5) + I(z3_post/5) + I(z45_post/5) + 
    val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + bin_statin +
    bin_acearb + bin_betablocker + bin_iabp + val_perfustm  +
    val_creatlst + val_predmort + val_proctime,
  data = depisodes.summary,
  family = "binomial"
)
summary(mglm.post)
ci.post <- confint(mglm.post)
coef.post <- coef(mglm.post)
dplot.post <- data.frame(
  phase = "Post-CPB",
  Zone = c("Zone 1", "Zone 2", "Zone 3"),
  coef = exp(coef.post[c("I(z12_post/5)", "I(z3_post/5)", "I(z45_post/5)")]),
  ci_lower = exp(ci.post[c("I(z12_post/5)", "I(z3_post/5)", "I(z45_post/5)"), 1]),
  ci_upper = exp(ci.post[c("I(z12_post/5)", "I(z3_post/5)", "I(z45_post/5)"), 2])
)
paste0(
  round(exp(coef.post["I(z12_post/5)"]), 3), 
  " (", round(exp(ci.post["I(z12_post/5)", 1]), 3), ", ", 
  round(exp(ci.post["I(z12_post/5)", 2]), 3), ")"
)
paste0(
  round(exp(coef.post["I(z3_post/5)"]), 3), 
  " (", round(exp(ci.post["I(z3_post/5)", 1]), 3), ", ", 
  round(exp(ci.post["I(z3_post/5)", 2]), 3), ")"
)
paste0(
  round(exp(coef.post["I(z45_post/5)"]), 3), 
  " (", round(exp(ci.post["I(z45_post/5)", 1]), 3), ", ", 
  round(exp(ci.post["I(z45_post/5)", 2]), 3), ")"
)

k_3_12 <- matrix(c(0, -1, 1, 0, rep(0, length(coef(mglm.post)) - 4)), 1)
t_3_12 <- glht(mglm.post, linfct = k_3_12)
summary(t_3_12)

k_45_3 <- matrix(c(0, 0, -1, 1, rep(0, length(coef(mglm.post)) - 4)), 1)
t_45_3 <- glht(mglm.post, linfct = k_45_3)
summary(t_45_3)

k_45_12 <- matrix(c(0, -1, 0, 1, rep(0, length(coef(mglm.post)) - 4)), 1)
t_45_12 <- glht(mglm.post, linfct = k_45_12)
summary(t_45_12)
confint(t_45_12)

dtest.post <- data.frame(
  phase = "Post-CPB",
  y = c(1.08, 1.11, 1.14),
  xmin = c("Zone 1", "Zone 2", "Zone 1"),
  xmax = c("Zone 2", "Zone 3", "Zone 3"),
  or = c(
    summary(t_3_12)$test$coefficients, 
    summary(t_45_3)$test$coefficients, 
    summary(t_45_12)$test$coefficients
  ),
  pval = p.adjust(c(
    summary(t_3_12)$test$pvalues, 
    summary(t_45_3)$test$pvalues, 
    summary(t_45_12)$test$pvalues
  ), method = "bonferroni")
) %>% mutate(
  or_pval = paste0(
    round(exp(or), 3), 
    " (p ", 
    ifelse(
      pval < 0.001, "< 0.001)", paste0("= ", round(pval, 3), ")")
    )
  )
)

# a <- exp(c(0.02863, 0.0395, 0.06813))
# b <- p.adjust(c(0.361, 0.0273, 0.00533), method = "bonferroni")
# paste0(round(a[1], 3), " (p = ", round(b[1], 3), ")")
# paste0(round(a[2], 3), " (p = ", round(b[2], 3), ")")
# paste0(round(a[3], 3), " (p = ", round(b[3], 3), ")")


# Make plot --------------
dtest <- rbind(dtest.pre, dtest.intra, dtest.post) %>% 
  mutate(phase = factor(phase, levels = c("Pre-CPB", "Intra-CPB", "Post-CPB"))) %>%
  mutate(
    xmin = factor(xmin, levels = c("Zone 1", "Zone 2", "Zone 3")),
    xmax = factor(xmax, levels = c("Zone 1", "Zone 2", "Zone 3"))
  )
dplot <- rbind(dplot.pre, dplot.intra, dplot.post) %>%
  mutate(phase = factor(phase, levels = c("Pre-CPB", "Intra-CPB", "Post-CPB"))) %>%
  mutate(Zone = factor(Zone, levels = c("Zone 1", "Zone 2", "Zone 3")))
custom_colors <- c("Zone 1" = "blue", "Zone 2" = "darkgrey", "Zone 3" = "red")
ggplot(dplot) + 
  geom_point(aes(x = Zone, y = coef, color = Zone, shape = Zone), size = 5) + 
  geom_errorbar(
    aes(x = Zone, ymin = ci_lower, ymax = ci_upper, color = Zone),
    linewidth = 1.2, width = 0.6
  ) + 
  geom_hline(
    yintercept = 1, alpha = 0.7, linetype = 2, linewidth = 1, color = "black"
  ) + 
  ggpubr::geom_bracket(
    data = dtest, 
    aes(xmin = xmin, xmax = xmax, y.position = y, label = or_pval)
  ) + 
  facet_wrap(~ phase) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0.9, 1.15), breaks = seq(0.9, 1.1, by = 0.05)) + 
  scale_color_manual(values = custom_colors) + 
  theme(legend.position = "top", legend.title = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),    
        axis.title = element_text(size = 14),
        legend.key.size = unit(0.35, units = "in"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  labs(x = "", y = "Odds Ratio (Per 5 Minutes)")







# ------------ Not for this manuscript ------------------



# Analysis 2: All zones and all phases
mglm.all <- glm(
  as.integer(bin_aki48h) ~
    I(z12_pre/5) + I(z3_pre/5) + I(z45_pre/5) +
    I(z12_intra/5) + I(z3_intra/5) + I(z45_intra/5) +
    I(z12_post/5) + I(z3_post/5) + I(z45_post/5) +
    val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + bin_statin +
    bin_acearb + bin_betablocker + bin_iabp  +
    val_creatlst + val_predmort,
  data = depisodes.summary,
  family = "binomial"
)
summary(mglm.all)
coef.all <- coef(mglm.all)
ci.all <- confint(mglm.all)
ci.pasted <- apply(round(exp(ci.all), 3), 1, function(x) paste(x[1], x[2], sep = ", "))
dtab <- data.frame(
  coef = exp(coef.all),
  ci = ci.pasted
)

paste(ci.all[1,1], ci.all[1,2], collapse = ", ")

k_2_1_pre <- matrix(c(
  0, -1, 1, 0, rep(0, 6), rep(0, length(coef.all) - 10)
), 1)
t_2_1_pre <- glht(mglm.all, linfct = k_2_1_pre)
summary(t_2_1_pre)
k_3_2_pre <- matrix(c(
  0, 0, -1, 1, rep(0, 6), rep(0, length(coef.all) - 10)
), 1)
t_3_2_pre <- glht(mglm.all, linfct = k_3_2_pre)
summary(t_3_2_pre)
k_3_1_pre <- matrix(c(
  0, -1, 0, 1, rep(0, 6), rep(0, length(coef.all) - 10)
), 1)
t_3_1_pre <- glht(mglm.all, linfct = k_3_1_pre)
summary(t_3_1_pre)

k_2_1_intra <- matrix(c(
  0, rep(0, 3), -1, 1, 0, rep(0, 3), rep(0, length(coef.all) - 10)
), 1)
t_2_1_intra <- glht(mglm.all, linfct = k_2_1_intra)
summary(t_2_1_intra)
k_3_2_intra <- matrix(c(
  0, rep(0, 3), 0, -1, 1, rep(0, 3), rep(0, length(coef.all) - 10)
), 1)
t_3_2_intra <- glht(mglm.all, linfct = k_3_2_intra)
summary(t_3_2_intra)
k_3_1_intra <- matrix(c(
  0, rep(0, 3), -1, 0, 1, rep(0, 3), rep(0, length(coef.all) - 10)
), 1)
t_3_1_intra <- glht(mglm.all, linfct = k_3_1_intra)
summary(t_3_1_intra)

k_2_1_post <- matrix(c(
  0, rep(0, 6), -1, 1, 0, rep(0, length(coef.all) - 10)
), 1)
t_2_1_post <- glht(mglm.all, linfct = k_2_1_post)
summary(t_2_1_post)
k_3_2_post <- matrix(c(
  0, rep(0, 6), 0, -1, 1, rep(0, length(coef.all) - 10)
), 1)
t_3_2_post <- glht(mglm.all, linfct = k_3_2_post)
summary(t_3_2_post)
k_3_1_post <- matrix(c(
  0, rep(0, 6), -1, 0, 1, rep(0, length(coef.all) - 10)
), 1)
t_3_1_post <- glht(mglm.all, linfct = k_3_1_post)
summary(t_3_1_post)

dplot.all <- data.frame(
  zone = rep(c("Zone 1", "Zone 2", "Zone 3"), 3),
  phase = rep(c("Pre-CPB", "Intra-CPB", "Post-CPB"), each  = 3),
  coef = exp(coef.all[c(
    "I(z12_pre/5)", "I(z3_pre/5)", "I(z45_pre/5)",
    "I(z12_intra/5)", "I(z3_intra/5)", "I(z45_intra/5)",
    "I(z12_post/5)", "I(z3_post/5)", "I(z45_post/5)"
  )]),
  ci_lower = exp(ci.all[c(
    "I(z12_pre/5)", "I(z3_pre/5)", "I(z45_pre/5)",
    "I(z12_intra/5)", "I(z3_intra/5)", "I(z45_intra/5)",
    "I(z12_post/5)", "I(z3_post/5)", "I(z45_post/5)"
  ), 1]),
  ci_upper = exp(ci.all[c(
    "I(z12_pre/5)", "I(z3_pre/5)", "I(z45_pre/5)",
    "I(z12_intra/5)", "I(z3_intra/5)", "I(z45_intra/5)",
    "I(z12_post/5)", "I(z3_post/5)", "I(z45_post/5)"
  ), 2])
) %>% mutate(
  phase = factor(phase, levels = c("Pre-CPB", "Intra-CPB", "Post-CPB"))
) %>% mutate(
  zone = factor(zone, levels = c("Zone 1", "Zone 2", "Zone 3"))
)
dtest.all <- data.frame(
  phase = rep(c("Pre-CPB", "Intra-CPB", "Post-CPB"), each  = 3),
  y = rep(c(1.08, 1.11, 1.14), 3),
  xmin = rep(c("Zone 1", "Zone 2", "Zone 1"), 3),
  xmax = rep(c("Zone 2", "Zone 3", "Zone 3"), 3),
  or = c(
    summary(t_2_1_pre)$test$coefficients, 
    summary(t_3_2_pre)$test$coefficients, 
    summary(t_3_1_pre)$test$coefficients,
    summary(t_2_1_intra)$test$coefficients, 
    summary(t_3_2_intra)$test$coefficients, 
    summary(t_3_1_intra)$test$coefficients,
    summary(t_2_1_post)$test$coefficients, 
    summary(t_3_2_post)$test$coefficients, 
    summary(t_3_1_post)$test$coefficients
  ),
  pval = c(
    summary(t_2_1_pre)$test$pvalues, 
    summary(t_3_2_pre)$test$pvalues, 
    summary(t_3_1_pre)$test$pvalues,
    summary(t_2_1_intra)$test$pvalues, 
    summary(t_3_2_intra)$test$pvalues, 
    summary(t_3_1_intra)$test$pvalues,
    summary(t_2_1_post)$test$pvalues, 
    summary(t_3_2_post)$test$pvalues, 
    summary(t_3_1_post)$test$pvalues
  )
) %>% mutate(
  or_pval = paste0(
    round(exp(or), 3), 
    " (p ", 
    ifelse(
      pval < 0.001, "< 0.001)", paste0("= ", round(pval, 3), ")")
    )
  )
) %>% 
  mutate(phase = factor(phase, levels = c("Pre-CPB", "Intra-CPB", "Post-CPB"))) %>%
  mutate(
    xmin = factor(xmin, levels = c("Zone 1", "Zone 2", "Zone 3")),
    xmax = factor(xmax, levels = c("Zone 1", "Zone 2", "Zone 3"))
  )
custom_colors <- c("Zone 1" = "blue", "Zone 2" = "darkgrey", "Zone 3" = "red")


# Without pairwise testing
pz <- c(
  "Pre-CPB 
      Zone 1", 
  "Pre-CPB
      Zone 2",
  "Pre-CPB
      Zone 3",
  "Intra-CPB 
      Zone 1", 
  "Intra-CPB
      Zone 2",
  "Intra-CPB
      Zone 3",
  "Post-CPB 
      Zone 1", 
  "Post-CPB
      Zone 2",
  "Post-CPB
      Zone 3"
)
dplot.all <- dplot.all %>%
  mutate(phase_zone = factor(pz, levels = pz))
ggplot(dplot.all) +
  geom_point(aes(x = phase_zone, y = coef, color = zone)) +
  geom_errorbar(aes(x = phase_zone, ymin = ci_lower, ymax = ci_upper, color = zone)) +
  geom_hline(yintercept = 1, alpha = 0.15) +
  theme_bw() +
  # facet_wrap(~ phase) + 
  scale_y_continuous(limits = c(0.95, 1.1), breaks = seq(0.95, 1.1, by = 0.05)) +
  scale_color_manual(values = custom_colors) +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),    
        axis.title = element_text(size = 14),
        legend.key.size = unit(0.35, units = "in"),
        legend.text = element_text(size = 12)) +
  labs(x = "", y = "Odds Ratio (Per 5 Minutes)")

