# - Objective 1: risk gradient.

# - Hypothesis 1. Figure 1: Primary Analysis: Three lines by zone, y: prob(AKI). 
#                 Fully adjusted for all covariates. 
# - Three marginal models. 

# - Hypothesis 2. Table 2: logistic regression coef + test for gradients. 
# - One conditional model. 


library(tidyverse)
library(mgcv)
library(tableone)
library(multcomp)

CalcMode <- function(x) {
  # Create a frequency table
  freq_table <- table(x)
  # Find the value(s) with the highest frequency
  mode_values <- names(freq_table[freq_table == max(freq_table)])
  # Return the mode(s)
  return(mode_values)
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
    id, bin_aki48h, 
    val_age, cat_gender, val_bmi, bin_hypertn, bin_dm,    
    bin_copd, bin_chf, bin_priormi, bin_cva, bin_pvd,       
    bin_emergent, bin_redo, val_hct, val_hdef, bin_statin,  
    bin_acearb, bin_betablocker, bin_iabp, val_perfustm, val_visscore,   
    val_creatlst, cat_transbc, val_predmort, val_crystalloid, val_proctime, val_predrenf
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
  )  %>%
  rename(phase = cat_cpb)
id.no.phase <- unique(filter(dat, is.na(phase))$id) 
dat <- dat %>% filter(!(id %in% id.no.phase))
  
depisodes.45 <- dat %>% 
  group_by(id) %>%
  reframe(CalcZoneLength(zone45)) 
depisodes.summary.45 <- depisodes.45 %>% 
  group_by(id) %>%
  summarise(
    nepisodes = n(),
    ztotal = sum(zlength)
  ) %>%
  arrange(desc(ztotal), desc(nepisodes)) %>%
  ungroup()
depisodes.3 <- dat %>% 
  group_by(id) %>%
  reframe(CalcZoneLength(zone3)) 
depisodes.summary.3 <- depisodes.3 %>% 
  group_by(id) %>%
  summarise(
    nepisodes = n(),
    ztotal = sum(zlength)
  ) %>%
  arrange(desc(ztotal), desc(nepisodes)) %>%
  ungroup()
depisodes.12 <- dat %>% 
  group_by(id) %>%
  reframe(CalcZoneLength(zone12)) 
depisodes.summary.12 <- depisodes.12 %>% 
  group_by(id) %>%
  summarise(
    nepisodes = n(),
    ztotal = sum(zlength)
  ) %>%
  arrange(desc(ztotal), desc(nepisodes)) %>%
  ungroup()




# Fit models for plots --------------------------------
# Hypothesis 1, Marginal risk of each zone
dgam.45 <- depisodes.summary.45 %>% 
  inner_join(dcovariates) %>% 
  mutate(across(starts_with("bin"), as.factor)) %>%
  mutate(across(starts_with("cat"), as.factor)) %>%
  mutate(across(starts_with("val"), function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T)))
mglm.45 <- glm(
  bin_aki48h ~ ztotal + val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + bin_statin +
    bin_acearb + bin_betablocker + bin_iabp + val_perfustm + val_visscore +
    val_creatlst + cat_transbc + val_predmort + val_crystalloid + val_proctime,
  data = dgam.45,
  family = "binomial"
)

dpred.45 <- dgam.45 %>% 
  dplyr::select(-c(id, nepisodes, ztotal, bin_aki48h)) %>% 
  summarise(across(everything(), ~ {
    if (is.numeric(.)) {
      mean(., na.rm = TRUE)  # Mean for numeric columns
    } else {
      CalcMode(.)  # Mode for categorical/binary columns
    }
  })) 
dpred.45.ztotal <- data.frame(ztotal = seq(0, max(dgam.45$ztotal), by = 1))
dpred.45 <- dpred.45[rep(1, nrow(dpred.45.ztotal)), ]
dpred.45 <- cbind(dpred.45, dpred.45.ztotal)
pglm.45.temp <- predict(
  mglm.45, newdata = dpred.45, 
  type = "link", se.fit = T
)
pglm.45 <- data.frame(
  x =  seq(0, max(dgam.45$ztotal), by = 1),
  fit_prob = plogis(pglm.45.temp$fit),
  ci_upper = plogis(pglm.45.temp$fit + 1.96 * pglm.45.temp$se.fit),
  ci_lower = plogis(pglm.45.temp$fit - 1.96 * pglm.45.temp$se.fit),
  zone = "Zone 3",
  method = "GLM"
)

dgam.3 <- depisodes.summary.3 %>% 
  inner_join(dcovariates) %>% 
  mutate(across(starts_with("bin"), as.factor)) %>%
  mutate(across(starts_with("cat"), as.factor)) %>%
  mutate(across(starts_with("val"), function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T)))
mglm.3 <- glm(
  bin_aki48h ~ ztotal + val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + bin_statin +
    bin_acearb + bin_betablocker + bin_iabp + val_perfustm + val_visscore +
    val_creatlst + cat_transbc + val_predmort + val_crystalloid + val_proctime,
  data = dgam.3,
  family = "binomial"
)

dpred.3 <- dgam.3 %>% 
  dplyr::select(-c(id, nepisodes, ztotal, bin_aki48h)) %>% 
  summarise(across(everything(), ~ {
    if (is.numeric(.)) {
      mean(., na.rm = TRUE)  # Mean for numeric columns
    } else {
      CalcMode(.)  # Mode for categorical/binary columns
    }
  })) 
dpred.3.ztotal <- data.frame(ztotal = seq(0, max(dgam.3$ztotal), by = 1))
dpred.3 <- dpred.3[rep(1, nrow(dpred.3.ztotal)), ]
dpred.3 <- cbind(dpred.3, dpred.3.ztotal)
pglm.3.temp <- predict(
  mglm.3, newdata = dpred.3, 
  type = "link", se.fit = T
)
pglm.3 <- data.frame(
  x =  seq(0, max(dgam.3$ztotal), by = 1),
  fit_prob = plogis(pglm.3.temp$fit),
  ci_upper = plogis(pglm.3.temp$fit + 1.96 * pglm.3.temp$se.fit),
  ci_lower = plogis(pglm.3.temp$fit - 1.96 * pglm.3.temp$se.fit),
  zone = "Zone 2",
  method = "GLM"
)

dgam.12 <- depisodes.summary.12 %>% 
  inner_join(dcovariates) %>% 
  mutate(across(starts_with("bin"), as.factor)) %>%
  mutate(across(starts_with("cat"), as.factor)) %>%
  mutate(across(starts_with("val"), function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T)))
mglm.12 <- glm(
  bin_aki48h ~ ztotal + val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + bin_statin +
    bin_acearb + bin_betablocker + bin_iabp + val_perfustm + val_visscore +
    val_creatlst + cat_transbc + val_predmort + val_crystalloid + val_proctime,
  data = dgam.12,
  family = "binomial"
)
dpred.12 <- dgam.12 %>% 
  dplyr::select(-c(id, nepisodes, ztotal, bin_aki48h)) %>% 
  summarise(across(everything(), ~ {
    if (is.numeric(.)) {
      mean(., na.rm = TRUE)  # Mean for numeric columns
    } else {
      CalcMode(.)  # Mode for categorical/binary columns
    }
  })) 
dpred.12.ztotal <- data.frame(ztotal = seq(0, max(dgam.12$ztotal), by = 1))
dpred.12 <- dpred.12[rep(1, nrow(dpred.12.ztotal)), ]
dpred.12 <- cbind(dpred.12, dpred.12.ztotal)
pglm.12.temp <- predict(
  mglm.12, newdata = dpred.12, 
  type = "link", se.fit = T
)
pglm.12 <- data.frame(
  x =  seq(0, max(dgam.12$ztotal), by = 1),
  fit_prob = plogis(pglm.12.temp$fit),
  ci_upper = plogis(pglm.12.temp$fit + 1.96 * pglm.12.temp$se.fit),
  ci_lower = plogis(pglm.12.temp$fit - 1.96 * pglm.12.temp$se.fit),
  zone = "Zone 1",
  method = "GLM"
)



# Fit models for marginal odds ratios -------------------------------------------
# For the table, present odds ratios per five minutes. 
mglm.12 <- glm(
  bin_aki48h ~ I(ztotal/5) + val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + bin_statin +
    bin_acearb + bin_betablocker + bin_iabp + val_perfustm + val_visscore +
    val_creatlst + cat_transbc + val_predmort + val_crystalloid + val_proctime,
  data = dgam.12,
  family = "binomial"
)
summary(mglm.12)$coefficients["I(ztotal/5)", ]
mglm.3 <- glm(
  bin_aki48h ~ I(ztotal/5) + val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + bin_statin +
    bin_acearb + bin_betablocker + bin_iabp + val_perfustm + val_visscore +
    val_creatlst + cat_transbc + val_predmort + val_crystalloid + val_proctime,
  data = dgam.3,
  family = "binomial"
)
summary(mglm.3)$coefficients["I(ztotal/5)", ]
mglm.45 <- glm(
  bin_aki48h ~ I(ztotal/5) + val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + bin_statin +
    bin_acearb + bin_betablocker + bin_iabp + val_perfustm + val_visscore +
    val_creatlst + cat_transbc + val_predmort + val_crystalloid + val_proctime,
  data = dgam.45,
  family = "binomial"
)
summary(mglm.45)$coefficients["I(ztotal/5)", ]

citemp.12 <- exp(confint(mglm.12)["I(ztotal/5)", ])
citemp.3 <- exp(confint(mglm.3)["I(ztotal/5)", ])
citemp.45 <- exp(confint(mglm.45)["I(ztotal/5)", ])
dtext <- data.frame(
  x = rep(3.5, 3),
  y = rep(0.75, 3),
  zone = c("Zone 1", "Zone 2", "Zone 3"),
  text = c(
    paste0(
      round(exp(coef(mglm.12)["I(ztotal/5)"]), 3), 
      " (", round(citemp.12[1], 3),", ", round(citemp.12[2], 3), ")"
    ),
    paste0(
      round(exp(coef(mglm.3)["I(ztotal/5)"]), 3), 
      " (", round(citemp.3[1], 3),", ", round(citemp.3[2], 3), ")"
    ),
    paste0(
      round(exp(coef(mglm.45)["I(ztotal/5)"]), 3), 
      " (", round(citemp.45[1], 3),", ", round(citemp.45[2], 3), ")"
    )
  ),
  pval = c(
    round(summary(mglm.12)$coefficients["I(ztotal/5)", "Pr(>|z|)"], 4),
    round(summary(mglm.3)$coefficients["I(ztotal/5)", "Pr(>|z|)"], 4),
    round(summary(mglm.45)$coefficients["I(ztotal/5)", "Pr(>|z|)"], 4)
  )
)

# Marginal prob plot ------------------------------------------------------------
pglm <- rbind(pglm.45, pglm.3, pglm.12)
custom_colors <- c("Zone 1" = "blue", "Zone 2" = "darkgrey", "Zone 3" = "red")
ggplot(pglm) +
  # geom_point(aes(x = ztotal, y = bin_aki48h), color = "blue", shape = 1) + 
  geom_ribbon(aes(x=x/60, ymin = ci_lower, ymax = ci_upper, fill = zone), alpha = 0.5) +
  geom_line(aes(x=x/60, y = fit_prob, color = zone), size = 1) +
  # geom_text(data = dtext, aes(x = x, y = y, label = text), size = 5) + 
  facet_wrap( ~ zone) + 
  labs(x = "Cumulative Time in Zone (Hours)", y = "Estimated AKI Probability") +
  theme_bw() + 
  theme(legend.position = "") + 
  scale_x_continuous(limits = c(0, 7.5), breaks = 0:8) + 
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.85, by = 0.1)) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors)






# Hypothesis 2: Risk Gradient -------------------------------------------------------
dglm <- rename(dplyr::select(depisodes.summary.45, -nepisodes), z45 = ztotal) %>% 
  left_join(rename(dplyr::select(depisodes.summary.3, -nepisodes), z3 = ztotal), by = "id") %>%
  left_join(rename(dplyr::select(depisodes.summary.12, -nepisodes), z12 = ztotal), by = "id") %>%
  inner_join(dcovariates, by = "id") %>%
  dplyr::select(-id) %>%
  mutate(across(starts_with("bin"), as.factor)) %>%
  mutate(across(starts_with("cat"), as.factor)) %>%
  mutate(across(starts_with("val"), function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T)))
  
dsummary <- dglm %>%
  mutate(id = 1:1196) %>%
  dplyr::select(id, bin_aki48h, z12, z3, z45) %>%
  rename(aki = bin_aki48h, zone1 = z12, zone2 = z3, zone3 = z45)

dsummary <- dsummary[sample(1:1196, size = 1196, replace = F), ] %>%
  mutate(id = 1:1196)

m2 <- glm(
  bin_aki48h ~ val_age + cat_gender + val_bmi + 
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd + 
    bin_emergent + bin_redo + val_hct + val_hdef + 
    bin_statin + bin_acearb + bin_betablocker + bin_iabp + val_perfustm + val_visscore + 
    val_creatlst + cat_transbc + val_crystalloid + val_predmort + 
    I(z45/5) + I(z3/5) + I(z12/5),
  data = dglm, 
  family = "binomial"
)
summary(m2)
ci.whole <- confint(m2)
coef.whole <- coef(m2)
dplot.whole <- data.frame(
  Zone = c("Zone 1", "Zone 2", "Zone 3"),
  coef = exp(coef.whole[c("I(z12/5)", "I(z3/5)", "I(z45/5)")]),
  ci_lower = exp(ci.whole[c("I(z12/5)", "I(z3/5)", "I(z45/5)"), 1]),
  ci_upper = exp(ci.whole[c("I(z12/5)", "I(z3/5)", "I(z45/5)"), 2])
)  %>% mutate(Zone = factor(Zone, levels = c("Zone 1", "Zone 2", "Zone 3")))

# Test for gradient
k_3_12 <- matrix(c(rep(0, length(coef(m2)) - 3), 0, 1, -1), 1)
t_3_12 <- glht(m2, linfct = k_3_12)
summary(t_3_12)
k_45_3 <- matrix(c(rep(0, length(coef(m2)) - 3), 1, -1, 0), 1)
t_45_3 <- glht(m2, linfct = k_45_3)
summary(t_45_3)
k_45_12 <- matrix(c(rep(0, length(coef(m2)) - 3), 1, 0, -1), 1)
t_45_12 <- glht(m2, linfct = k_45_12)
summary(t_45_12)

dtest.whole <- data.frame(
  y = c(1.05, 1.07, 1.09),
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
) %>% mutate(
    xmin = factor(xmin, levels = c("Zone 1", "Zone 2", "Zone 3")),
    xmax = factor(xmax, levels = c("Zone 1", "Zone 2", "Zone 3"))
)
dtest.whole

custom_colors <- c("Zone 1" = "blue", "Zone 2" = "darkgrey", "Zone 3" = "red")
ggplot(dplot.whole) + 
  geom_point(aes(x = Zone, y = coef, color = Zone)) + 
  geom_errorbar(aes(x = Zone, ymin = ci_lower, ymax = ci_upper, color = Zone)) + 
  geom_hline(yintercept = 1, alpha = 0.15) + 
  ggpubr::geom_bracket(
    data = dtest.whole, 
    aes(xmin = xmin, xmax = xmax, y.position = y, label = or_pval)
  ) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0.95, 1.1), breaks = seq(0.95, 1.1, by = 0.025)) + 
  scale_color_manual(values = custom_colors) + 
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(x = "", y = "Odds Ratio (Per 5 Minutes)")








