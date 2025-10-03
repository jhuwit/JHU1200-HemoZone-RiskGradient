

library(tidyverse)
library(multcomp)
library(stargazer)

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
    cond1 = map <= 55 & (cvp >=0 & cvp <= 20), 
    cond2 = (map > 55 & map <= 65) & (cvp > 8 & cvp <= 20), 
    cond3 = (map > 65 & map <= 115) & (cvp > 12 & cvp <= 20),
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
  inner_join(dcovariates) %>%
  mutate(across(starts_with("cat"), as.factor)) %>%
  mutate(across(starts_with("val"), function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T)))



# Model for extracting OR for every 5 minutes ------------------------------------------
mglm.all <- glm(
  as.integer(bin_aki48h) ~ 
    I(z12_pre/5) + I(z3_pre/5) + I(z45_pre/5) + 
    I(z12_intra/5) + I(z3_intra/5) + I(z45_intra/5) + 
    I(z12_post/5) + I(z3_post/5) + I(z45_post/5) +
    val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + 
    bin_statin + bin_acearb + bin_betablocker + bin_iabp + 
    val_creatlst + val_predmort +
    cat_transbc + val_visscore + val_crystalloid,
  data = depisodes.summary,
  family = "binomial"
)
summary(mglm.all)
coef.all <- coef(mglm.all)
ci.all <- confint(mglm.all)
vars = c(
  "I(z12_pre/5)", "I(z3_pre/5)", "I(z45_pre/5)",
  "I(z12_intra/5)", "I(z3_intra/5)", "I(z45_intra/5)",
  "I(z12_post/5)", "I(z3_post/5)", "I(z45_post/5)"
)
custom_colors <- c("Zone 1" = "blue", "Zone 2" = "darkgrey", "Zone 3" = "red")
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
pz = factor(pz, levels = pz)
dplot <- data.frame(
  zone = rep(paste("Zone", 1:3), 3),
  phase_zone = pz,
  coef = exp(coef.all[vars]),
  ci_lower = exp(ci.all[vars, 1]),
  ci_upper = exp(ci.all[vars, 2])
) %>% mutate(phase_zone)
ggplot(dplot) +
  geom_point(aes(x = phase_zone, y = coef, color = zone, shape = zone), size = 5) +
  geom_errorbar(
    aes(x = phase_zone, ymin = ci_lower, ymax = ci_upper, color = zone), linewidth = 1.2, width = 0.6
  ) +
  geom_hline(
    yintercept = 1, alpha = 0.7, linetype = 2, linewidth = 1, color = "black"
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0.93, 1.1), breaks = seq(0.95, 1.1, by = 0.05)) +
  scale_color_manual(values = custom_colors) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 11),      
        axis.title = element_text(size = 14)) +
  labs(x = "", y = "Odds Ratio (Per 5 Minutes)")
p.val <- cut(
  summary(mglm.all)$coefficients[, "Pr(>|z|)"],
  breaks = c(0, 0.001, 0.01, 0.05, 1),
  labels = c("***", "**", "*", "")
)
var.names <- c(
  "Intercept", 
  "Zone 1 Pre-CPB", "Zone 2 Pre-CPB", "Zone 3 Pre-CPB",
  "Zone 1 Intra-CPB", "Zone 2 Intra-CPB", "Zone 3 Intra-CPB",
  "Zone 1 Post-CPB", "Zone 2 Post-CPB", "Zone 3 Post-CPB",
  "Age", "Gender:Male", "BMI", "Hypertension", "Diabetes",
  "COPD", "CHF", "Prior MI", "CVA", "PVD", "Emergent", "Redo",
  "HCT", "HDEF", "Statin", "Acearb", "Betablocker", "IADP",
  "Baseline Creatinine", "Predicted Mortality", 
  "TrasBC = 0", "TransBC = 1", "VIS Score", "Crystalloid"
)
coef.all <- exp(coef.all)
ci.all <- exp(ci.all)
tab.m <- data.frame(
  Variable = var.names,
  OR = round(coef.all, 3),
  CI = paste0(
    "(", round(ci.all[, 1], 3), ", ", round(ci.all[, 2], 3), ")", 
    p.val
  )
)
row.names(tab.m) <- NULL








# Model for calculating probability -------------------------------------
mglm.all <- glm(
  bin_aki48h ~ 
    z12_pre + z3_pre + z45_pre + 
    z12_intra + z3_intra + z45_intra + 
    z12_post + z3_post + z45_post + 
    val_age + cat_gender + val_bmi +
    bin_hypertn + bin_dm + bin_copd + bin_chf + bin_priormi + bin_cva + bin_pvd +
    bin_emergent + bin_redo + val_hct + val_hdef + 
    bin_statin + bin_acearb + bin_betablocker + bin_iabp + 
    val_creatlst + val_predmort +
    cat_transbc + val_visscore + val_crystalloid,
  data = depisodes.summary,
  family = "binomial"
)

# Impute mean for exposures that have not occurred ------------
dpred.baseline <- depisodes.summary %>%
  mutate(across(ends_with(c("pre", "intra", "post")), ~ mean(., na.rm = TRUE)))
ppred.baseline <- predict(
  mglm.all, newdata = dpred.baseline, 
  type = "link", se.fit = F
)
dpred.precpb <- depisodes.summary %>%
  mutate(across(ends_with(c("intra", "post")), ~ mean(., na.rm = TRUE)))
ppred.precpb <- predict(
  mglm.all, newdata = dpred.precpb, 
  type = "link", se.fit = F
)
dpred.intracpb <- depisodes.summary %>%
  mutate(across(ends_with(c("post")), ~ mean(., na.rm = TRUE)))
ppred.intracpb <- predict(
  mglm.all, newdata = dpred.intracpb, 
  type = "link", se.fit = F
)
ppred.postcpb <- predict(
  mglm.all, newdata = depisodes.summary, 
  type = "link", se.fit = F
)
dres <- cbind(
  depisodes.summary,
  data.frame(
    prob_Baseline = plogis(ppred.baseline),
    prob_PreCPB = plogis(ppred.precpb),
    prob_IntraCPB = plogis(ppred.intracpb),
    prob_PostCPB = plogis(ppred.postcpb)
  )
)
q.prob.baseline <- quantile(dres$prob_Baseline, probs = seq(0, 1, by = 0.1))
dres.catbaseline <- dres %>%
  mutate(
    cat_prob_baseline = cut(
      x = prob_Baseline, labels = paste("Decile", 1:10),
      breaks = q.prob.baseline, include.lowest = T, right = T, ordered_result = T
    )
  ) %>% 
  mutate(
    prob_change = prob_PostCPB - prob_Baseline
  ) %>%
  mutate(
    cat_prob_change = cut(
      x = prob_change, labels = c("<= -0.05", "(-0.05, 0.05]", "> 0.05"),
      breaks = c(-1, -0.05, 0.05, 1), include.lowest = T, right = T, ordered_result = T
    )
  )

dres.catbaseline.long <- dres.catbaseline %>% 
  dplyr::select(
    id, bin_aki48h, prob_Baseline, prob_PreCPB, prob_IntraCPB, prob_PostCPB,
    cat_prob_baseline, cat_prob_change
  ) %>%
  pivot_longer(
    cols = -c(id, bin_aki48h, cat_prob_baseline, cat_prob_change), 
    names_to = "time", values_to = "prob"
  ) %>%
  mutate(time = gsub("prob_", "", time)) %>%
  mutate(
    time = factor(
      time, levels = c("Baseline", "PreCPB", "IntraCPB", "PostCPB"),
      labels = c("Base", "Pre", "Intra", "Post")
    )
  )
custom_colors <- c("<= -0.05" = "#1A85FF", "(-0.05, 0.05]" = "grey", "> 0.05" = "#D41159")
pa <- ggplot(dres.catbaseline.long) + 
  geom_line(aes(x = time, y = prob, group = id, color = cat_prob_change)) + 
  facet_wrap(~ cat_prob_baseline, nrow = 2) + 
  labs(x = "Time", y = "Estimated AKI Probability") +
  theme_bw() + 
  theme(legend.position = "top", 
        legend.key.size = unit(0.35, units = "in"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust = 0.25),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 13)) +
  scale_x_discrete(expand = c(0.05, 0.05)) + 
  guides(color = guide_legend(override.aes = list(size = 2))) +
  labs(color = "Panel A.        Change in AKI Probability") + 
  scale_color_manual(values = custom_colors)




# Take one subject, modify Pre exposure, show effect -------------------
# One subject from Decile 4
eg.id <- dres.catbaseline$id[
  with(dres.catbaseline, which(cat_prob_baseline == "Decile 4" & prob_PreCPB > 0.5))
] # "30514_377019"
eg.id <- "30514_377019"
dpred.eg <- filter(dres.catbaseline, id == eg.id)
# calculate baseline prob. condition on the population mean for exposure. 
dpred.eg.baseline <- depisodes.summary %>%
  mutate(across(ends_with(c("pre", "intra", "post")), ~ mean(., na.rm = TRUE))) %>%
  filter(id == eg.id)
ppred.eg.baseline <- predict(
  mglm.all, newdata = dpred.eg.baseline, 
  type = "link", se.fit = F
)
plogis(ppred.eg.baseline)

# calculate preCPB prob. 
dpred.eg.precpb <- depisodes.summary %>%
  mutate(across(ends_with(c("intra", "post")), ~ mean(., na.rm = TRUE))) %>%
  filter(id == eg.id)
dpred.eg.precpb <- rbind(
  dpred.eg.precpb,
  dpred.eg.precpb %>% mutate(z12_pre = z12_pre + 30, z3_pre = z3_pre, z45_pre = z45_pre - 30),
  dpred.eg.precpb %>% mutate(z12_pre = z12_pre + 60, z3_pre = z3_pre, z45_pre = z45_pre - 60),
  dpred.eg.precpb %>% mutate(z12_pre = z12_pre + z3_pre + z45_pre, z3_pre = z3_pre - z3_pre, z45_pre = z45_pre - z45_pre)
) 
ppred.eg.precpb <- predict(
  mglm.all, newdata = dpred.eg.precpb, 
  type = "link", se.fit = F
) 
plogis(ppred.eg.precpb)

# calculate intraCPB prob. 
dpred.eg.intracpb <- depisodes.summary %>%
  mutate(across(ends_with(c("post")), ~ mean(., na.rm = TRUE))) %>%
  filter(id == eg.id)
dpred.eg.intracpb <- rbind(
  dpred.eg.intracpb,
  dpred.eg.intracpb %>% mutate(z12_pre = z12_pre + 30, z3_pre = z3_pre, z45_pre = z45_pre - 30),
  dpred.eg.intracpb %>% mutate(z12_pre = z12_pre + 60, z3_pre = z3_pre, z45_pre = z45_pre - 60),
  dpred.eg.intracpb %>% mutate(z12_pre = z12_pre + z3_pre + z45_pre, z3_pre = z3_pre - z3_pre, z45_pre = z45_pre - z45_pre)
) 
ppred.eg.intracpb <- predict(
  mglm.all, newdata = dpred.eg.intracpb, 
  type = "link", se.fit = F
) 
plogis(ppred.eg.intracpb)

# calculate postCPB prob. 
dpred.eg.postcpb <- depisodes.summary %>%
  filter(id == eg.id)
dpred.eg.postcpb <- rbind(
  dpred.eg.postcpb,
  dpred.eg.postcpb %>% mutate(z12_pre = z12_pre + 30, z3_pre = z3_pre, z45_pre = z45_pre - 30),
  dpred.eg.postcpb %>% mutate(z12_pre = z12_pre + 60, z3_pre = z3_pre, z45_pre = z45_pre - 60),
  dpred.eg.postcpb %>% mutate(z12_pre = z12_pre + z3_pre + z45_pre, z3_pre = z3_pre - z3_pre, z45_pre = z45_pre - z45_pre)
) 
ppred.eg.postcpb <- predict(
  mglm.all, newdata = dpred.eg.postcpb, 
  type = "link", se.fit = F
) 
plogis(ppred.eg.postcpb)

dprob.eg <- 
  data.frame(
    scenario = c(
      "Original", 
      "Scenario 1",
      "Scenario 2",
      "Scenario 3"
    ),
    prob_Baseline = rep(plogis(ppred.eg.baseline), 4),
    prob_PreCPB = plogis(ppred.eg.precpb),
    prob_IntraCPB = plogis(ppred.eg.intracpb),
    prob_PostCPB = plogis(ppred.eg.postcpb)
  )
dprob.eg.long <- dprob.eg %>% 
  pivot_longer(
    cols = -scenario, 
    names_to = "time", values_to = "prob"
  ) %>%
  mutate(time = gsub("prob_", "", time)) %>%
  mutate(
    time = factor(
      time, levels = c("Baseline", "PreCPB", "IntraCPB", "PostCPB"),
      labels = c("Baseline", "Pre-CPB", "Intra-CPB", "Post-CPB")
    )
  ) %>% mutate(
    scenario = factor(scenario, levels = c("Original", paste("Scenario", 1:3)))
  )

pb <- ggplot(dprob.eg.long) + 
  geom_line(aes(x = time, y = prob, group = scenario, color = scenario)) + 
  geom_rect(
    data = data.frame(xfrom = "Baseline", xto = "Pre-CPB", yfrom = 0.19, yto = 0.75), 
    aes(xmin = xfrom, xmax = xto, ymin = yfrom, ymax = yto), alpha = 0.08
  ) + 
  geom_text(
    data = data.frame(
      x = rep("Post-CPB", 4),
      y = filter(dprob.eg.long, time == "Post-CPB")$prob,
      scenario = c("Original", paste("Scenario", 1:3))
    ),
    aes(x = x, y = y, label = scenario, color = scenario),
    hjust = -0.1
  ) + 
  scale_y_continuous(
    limits = c(0.19, 0.75), 
    breaks = seq(0.15, 0.75, by = 0.05), 
    minor_breaks = NULL
  ) +
  scale_x_discrete(expand = c(0.02, 0.01, 0.2, 0.1)) + 
  theme_bw() + 
  labs(x = "Time", y = "Estimated AKI Probability") +
  theme(
    legend.position = "none", 
    legend.key.size = unit(0.35, units = "in"),
    plot.title = element_text(hjust=0.5), 
    axis.text = element_text(size = 12),      
    axis.title = element_text(size = 14)
  ) +
  ggtitle("Panel B.      Modifying Pre-CPB Exposure")







# Take one subject, modify Post exposure, show effect -------------------

dres.catbaseline.subset <- dres.catbaseline %>%
  mutate(
    zone_total = z12_pre + z3_pre + z45_pre + 
      z12_intra + z3_intra + z45_intra + 
      z12_post + z3_post + z45_post,
    post_total = z12_post + z3_post + z45_post
  ) %>%
  filter(val_proctime - zone_total <= 60) %>%
  arrange(desc(post_total))
eg.id <- "85645_505323"
dpred.eg <- filter(dres.catbaseline, id == eg.id)
# calculate baseline prob. condition on the population mean for exposure. 
dpred.eg.baseline <- depisodes.summary %>%
  mutate(across(ends_with(c("pre", "intra", "post")), ~ mean(., na.rm = TRUE))) %>%
  filter(id == eg.id)
ppred.eg.baseline <- predict(
  mglm.all, newdata = dpred.eg.baseline, 
  type = "link", se.fit = F
)
plogis(ppred.eg.baseline)

# calculate preCPB prob. 
dpred.eg.precpb <- depisodes.summary %>%
  mutate(across(ends_with(c("intra", "post")), ~ mean(., na.rm = TRUE))) %>%
  filter(id == eg.id)
ppred.eg.precpb <- predict(
  mglm.all, newdata = dpred.eg.precpb, 
  type = "link", se.fit = F
) 
plogis(ppred.eg.precpb)

# calculate intraCPB prob. 
dpred.eg.intracpb <- depisodes.summary %>%
  mutate(across(ends_with(c("post")), ~ mean(., na.rm = TRUE))) %>%
  filter(id == eg.id)
ppred.eg.intracpb <- predict(
  mglm.all, newdata = dpred.eg.intracpb, 
  type = "link", se.fit = F
) 
plogis(ppred.eg.intracpb)

# calculate postCPB prob. 
dpred.eg.postcpb <- depisodes.summary %>%
  filter(id == eg.id)
dpred.eg.postcpb %>% dplyr::select(z12_post, z3_post, z45_post)
dpred.eg.postcpb <- rbind(
  dpred.eg.postcpb,
  dpred.eg.postcpb %>% mutate(
    z12_post = z12_post - 30, z3_post = z3_post, 
    z45_post = z45_post + 30
  ),
  dpred.eg.postcpb %>% mutate(
    z12_post = z12_post - 60, z3_post = z3_post + 30, 
    z45_post = z45_post + 30
  ),
  dpred.eg.postcpb %>% mutate(
    z45_post = z12_post + z3_post + z45_post,
    z12_post = z12_post - z12_post, 
    z3_post = z3_post - z3_post
  )
) 
ppred.eg.postcpb <- predict(
  mglm.all, newdata = dpred.eg.postcpb, 
  type = "link", se.fit = F
) 
plogis(ppred.eg.postcpb)

# depisodes.summary %>% filter(id == eg.id)
dprob.eg <- 
  data.frame(
    scenario = c(
      "Original",
      "Scenario 1",
      "Scenario 2",
      "Scenario 3"
    ),
    prob_Baseline = c(plogis(ppred.eg.baseline), rep(NA, 3)),
    prob_PreCPB = c(plogis(ppred.eg.precpb), rep(NA, 3)),
    prob_IntraCPB = rep(plogis(ppred.eg.intracpb), 4),
    prob_PostCPB = plogis(ppred.eg.postcpb)
  )

dprob.eg.long <- dprob.eg %>% 
  pivot_longer(
    cols = -scenario, 
    names_to = "time", values_to = "prob"
  ) %>%
  mutate(time = gsub("prob_", "", time)) %>%
  mutate(
    time = factor(
      time, levels = c("Baseline", "PreCPB", "IntraCPB", "PostCPB"),
      labels = c("Baseline", "Pre-CPB", "Intra-CPB", "Post-CPB")
    )
  ) %>%
  mutate(
    scenario = factor(scenario, levels = c("Original", paste("Scenario", 1:3)))
  )
pc <- ggplot(dprob.eg.long) + 
  geom_line(aes(x = time, y = prob, group = scenario, color = scenario)) + 
  geom_rect(
    data = data.frame(xfrom = "Intra-CPB", xto = "Post-CPB", yfrom = 0.1, yto = 0.25), 
    aes(xmin = xfrom, xmax = xto, ymin = yfrom, ymax = yto), alpha = 0.08
  ) +
  geom_text(
    data = data.frame(
      x = rep("Post-CPB", 4),
      y = filter(dprob.eg.long, time == "Post-CPB")$prob,
      scenario = c("Original", paste("Scenario", 1:3))
    ),
    aes(x = x, y = y, label = scenario, color = scenario),
    hjust = -0.1
  ) + 
  ggtitle("Panel C.      Modifying Post-CPB Exposure") +
  scale_y_continuous(
    limits = c(0.1, 0.25), breaks = seq(0.1, 0.25, by = 0.025), 
    minor_breaks = NULL
  ) +
  scale_x_discrete(expand = c(0.02, 0.01, 0.2, 0.1)) + 
  theme_bw() + 
  labs(x = "Time", y = "Estimated AKI Probability") +
  theme(
    legend.position = "none", 
    legend.key.size = unit(0.35, units = "in"),
    plot.title = element_text(hjust=0.5), 
    axis.text = element_text(size = 12),      
    axis.title = element_text(size = 14)
  ) 



library(patchwork)
pa/(pb/pc + plot_layout(axis_titles = "collect")) + 
  plot_layout(heights = c(5, 5))

