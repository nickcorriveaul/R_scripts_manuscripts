### Code for the statistical analyses of the manuscript "Clinical criteria for a limbic-predominant amnestic neurodegenerative syndrome" by Corriveau-Lecavalier et al. 

# Clear workspace --------------------------------------------------------------
rm(list=ls())

# Parameters  ------------------------------------------------------------------
options(max.print=1000000000)
options(scipen = 999, "digits"=3)

# Set working directory --------------------------------------------------------
setwd("/yourdirectory/")

# Loading relevant libraries ---------------------------------------------------
library(pROC); library(ggplot2); library(hrbrthemes); library(RVAideMemoire); 
library(rcompanion); library(caret); library(dplyr); library(lme4); 
library(patchwork); library(lm.beta); library(emmeans)

# Assessment of demographic variable differences -------------------------------
df = read.delim('XXX.csv', header = T, sep = ',')

df_baseline = df %>%
  group_by(clinic) %>%
  slice_min(visitage) %>%
  ungroup()

table(df_baseline$Path_dx)

t_sex = table(df_baseline$Path_dx, df_baseline$gender)
chisq.test(t_sex)
chisq.multcomp(t_sex, p.method = "fdr")

t_diag = table(df_baseline$Path_dx, df_baseline$Diagnosis)
chisq.test(t_diag)
chisq.multcomp(t_diag, p.method = "fdr")

t_apoe = table(df_baseline$Path_dx, df_baseline$apoe)
chisq.test(t_apoe)
chisq.multcomp(t_apoe, p.method = "fdr")

df_baseline %>%                               
  group_by(Path_dx) %>% 
  summarize(median = median(visitage, na.rm = T),
            q1 = quantile(visitage, 0.25, na.rm = T),
            q3 = quantile(visitage, 0.75, na.rm = T))

anova_age = aov(data=df_baseline, visitage ~ Path_dx)
summary(anova_age)
TukeyHSD(anova_age)

df_baseline %>%                               
  group_by(Path_dx) %>% 
  summarize(median = median(deathage, na.rm = T),
            q1 = quantile(deathage, 0.25, na.rm = T),
            q3 = quantile(deathage, 0.75, na.rm = T))

anova_agedeath = aov(data=df_baseline, deathage ~ Path_dx)
summary(anova_agedeath)
TukeyHSD(anova_agedeath)

df_baseline %>%                               
  group_by(Path_dx) %>% 
  summarize(median = median(education, na.rm = T),
            q1 = quantile(education, 0.25, na.rm = T),
            q3 = quantile(education, 0.75, na.rm = T))

anova_education = aov(data=df_baseline, education ~ Path_dx)
summary(anova_education)
TukeyHSD(anova_education)

df_baseline %>%                               
  group_by(Path_dx) %>% 
  summarize(median = median(calc.cdrsum, na.rm = T),
            q1 = quantile(calc.cdrsum, 0.25, na.rm = T),
            q3 = quantile(calc.cdrsum, 0.75, na.rm = T))

anova_cdrsob = aov(data=df_baseline, calc.cdrsum ~ Path_dx)
summary(anova_cdrsob)
TukeyHSD(anova_cdrsob)

df_lastpib = df %>%
  group_by(clinic) %>%
  slice_max(pibdate) %>%
  ungroup()

df_lastpib %>%                               
  group_by(Path_dx) %>% 
  summarize(median = median(amyloid_centiloid, na.rm = T),
            q1 = quantile(amyloid_centiloid, 0.25, na.rm = T),
            q3 = quantile(amyloid_centiloid, 0.75, na.rm = T))

anova_pib = aov(data=df_lastpib, amyloid_centiloid ~ Path_dx)
summary(anova_pib)
TukeyHSD(anova_pib)

df_lasttau = df %>%
  group_by(clinic) %>%
  slice_max(taudate) %>%
  ungroup()

df_lasttau %>%                               
  group_by(Path_dx) %>% 
  summarize(median = median(spm12.tau.ratio, na.rm = T),
            q1 = quantile(spm12.tau.ratio, 0.25, na.rm = T),
            q3 = quantile(spm12.tau.ratio, 0.75, na.rm = T))

anova_av1451 = aov(data=df_lasttau, spm12.tau.ratio ~ Path_dx)
summary(anova_av1451)
TukeyHSD(anova_av1451)

#Definition of MRI cut-point --------------------------------------------------
df = read.delim('XXX.csv', header = T, sep = ',')

# fit linear model
lm_hippvol = lm(spm12.hva ~ calc.cdrsum + (1 | clinic), data = df)
summary(lm_hippvol)

# extract and scale residuals
fitted_spm12.hva = predict(lm_hippvol, df)
df_fitted_hippvol = cbind(df, fitted_spm12.hva)
df_fitted_hippvol$spm12.hva_residuals = df_fitted_hippvol$spm12.hva - df_fitted_hippvol$fitted_spm12.hva
df_fitted_hippvol$spm12.hva_residuals_z = scale(df_fitted_hippvol$spm12.hva_residuals)
df = cbind(df, df_fitted_hippvol$spm12.hva_residuals_z)

# determine threshold using a logistic regression and ROC analysis on the last MRI
df$Path_dx[df$Path_dx == "ADNC"] = 0
df$Path_dx[df$Path_dx == "LBD"] = 0
df$Path_dx[df$Path_dx == "ADNC/LBD"] = 0
df$Path_dx[df$Path_dx == "PART"] = 0
df$Path_dx[df$Path_dx == "LATE-NC"] = 1
df$Path_dx[df$Path_dx == "ADNC/LATE-NC"] = 1
df$Path_dx = as.numeric(df$Path_dx)

df_lastmri = df %>%
  group_by(clinic) %>%
  slice_max(mridate) %>%
  ungroup()

model <- glm(Path_dx ~ spm12.hva_residuals_z, data = df_lastmri, family = "binomial")
summary(model)
roc_obj <- roc(df_lastmri$Path_dx, df_lastmri$spm12.hva_residuals_z)
optimal_threshold <- coords(roc_obj, "best", ret="threshold")$threshold

#extracting residuals from the Mayo cohort using the ADNI model
df2 = read.delim('XXX.csv', header = T, sep = ',')

fitted_spm12.hva = predict(lm_hippvol, df2)
df_fitted_hippvol = cbind(df2, fitted_spm12.hva)
df_fitted_hippvol$spm12.hva_residuals = df_fitted_hippvol$spm12.hva - df_fitted_hippvol$fitted_spm12.hva
df_fitted_hippvol$spm12.hva_residuals_z = scale(df_fitted_hippvol$spm12.hva_residuals)
df = cbind(df, df_fitted_hippvol$spm12.hva_residuals_z)

#Defining FDG-PET IMT cut-point ------------------------------------------------
df = read.delim('XXX.csv', header = T, sep = ',')
df$Path_dx[df$Path_dx == "ADNC"] = 0
df$Path_dx[df$Path_dx == "LBD"] = 0
df$Path_dx[df$Path_dx == "ADNC/LBD"] = 0
df$Path_dx[df$Path_dx == "PART"] = 0
df$Path_dx[df$Path_dx == "LATE-NC"] = 1
df$Path_dx[df$Path_dx == "ADNC/LATE-NC"] = 1
df$Path_dx = as.numeric(df$Path_dx)

df_lastfdg = df %>%
  group_by(clinic) %>%
  slice_max(fdgdate) %>%
  ungroup()

roc_obj <- roc(df_lastfdg$Path_dx, df_lastfdg$imt_ratio)
optimal_threshold <- coords(roc_obj, "best", ret="threshold")$threshold 

# Assessment of pathological diagnoses across LANS likelihoods -----------------
df = read.delim('XXX.csv', header = T, sep = ',')

df_baseline = df %>%
  group_by(clinic) %>%
  slice_min(visitage) %>%
  ungroup()

df_baseline$likelihood[df_baseline$likelihood == "Highest"] <- "High"

t_likelihood = table(df_baseline$Path_dx, df_baseline$likelihood)
t_likelihood = t_likelihood + 1
chisq.test(t_likelihood)
chisq.multcomp(t_likelihood, p.method = "fdr")

#GLM model for the binary classification of LATE -------------------------------
#Fitting the model in the Mayo cohort
df = read.delim('XXX.csv', header = T, sep = ',')

df_baseline = df %>%
  group_by(clinic) %>%
  slice_min(visitage) %>%
  ungroup()

df_m0 = df_baseline[,c("clinic", "Path_dx", "visitage", "calc.cdrsum", "imt_ratio", "spm12.fdg.ad.ratio", "spm12.hva_residuals_z", "tau_status")]
df_m0$Path_dx[df_m0$Path_dx == "AD"] <- 0
df_m0$Path_dx[df_m0$Path_dx == "AD-LATE"] <- 1
df_m0$Path_dx[df_m0$Path_dx == "Pure LATE"] <- 1
df_m0$Path_dx = as.numeric(df_m0$Path_dx)

glm.fit = glm(Path_dx ~ visitage + calc.cdrsum + imt_ratio + spm12.fdg.ad.ratio + spm12.hva_residuals_z + tau_status, family=binomial, data=df_m0)
summary(glm.fit)

par(pty="s")

predicted_probs <- predict(glm.fit, newdata = df_m0, type = "response")
predicted_class <- ifelse(predicted_probs >= 0.5, 1, 0)
mat = cbind(df_m0$Path_dx, predicted_class)
colnames(mat) = c("True", "Predicted")
mat = as.data.frame(mat)
mat = table(mat$True, mat$Predicted)
confusionMatrix(mat)

#Out of sample predictions in the ADNI cohort
df = read.delim('XXX.csv', header = T, sep = ',')

df_baseline = df %>%
  group_by(clinic) %>%
  slice_min(visitage) %>%
  ungroup()

df_m1 = df_baseline[,c("clinic", "Path_dx", "visitage", "calc.cdrsum", "imt_ratio", "spm12.fdg.ad.ratio", "spm12.hva_residuals_z", "tau_status")]
df_m1$Path_dx[df_m1$Path_dx == "AD"] <- 0
df_m1$Path_dx[df_m1$Path_dx == "AD-LATE"] <- 1
df_m1$Path_dx[df_m1$Path_dx == "Pure LATE"] <- 1
df_m1$Path_dx = as.numeric(df_m1$Path_dx)

predicted_probs <- predict(glm.fit, newdata = df_m1, type = "response")
predicted_class <- ifelse(predicted_probs >= 0.5, 1, 0)
mat = cbind(df_m1$Path_dx, predicted_class)
colnames(mat) = c("True", "Predicted")
mat = as.data.frame(mat)
mat = table(mat$True, mat$Predicted)
confusionMatrix(mat)

### Visualization and analysis of CDR-SB trajectories by LANS likelihoods ------
df = read.delim('XXX.csv', header = T, sep = ',')
df$likelihood = factor(df$likelihood, levels = c("Highest", "High", "Moderate", "Low"))

cleanup = theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "white", colour = "white"),
                axis.text =element_text(size=15, face= "bold", color = "black"), 
                axis.title = element_text(size = 18, face = "bold", color = "black"),
                plot.title = element_text(size = 20, face = "bold", color = "black"),
                legend.text = element_text(size=12, face="bold", color = "black"),
                axis.line = element_line(colour = "black"),
                legend.title = element_blank())

p_calc.cdrsum_long = ggplot(data=df, aes(x = time_since_first, y = calc.cdrsum, color = likelihood)) +
  geom_line(aes(group = clinic), alpha = 0.3) +
  geom_point(aes(fill = likelihood), pch = 21, size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", aes(group = likelihood, fill = likelihood), size = 2, alpha = 0.3) +
  scale_fill_manual(values=c('#997700', '#882255', '#0077BB', '#117733')) +
  scale_color_manual(values=c('#997700', '#882255', '#0077BB', '#117733')) +
  labs(y = "CDR-SB", x = "Time since baseline (years)", fill = "Diagnosis") + 
  cleanup  +
  theme(legend.position = c(0.5, 0.95), legend.direction = "horizontal") +
  guides(fill = guide_legend(override.aes = list(color = NULL)))

lm_calc.cdrsum = lmer(calc.cdrsum ~ likelihood*time_since_first + (1 | clinic), data = df)
summary(lm_calc.cdrsum)

trends.emm = emtrends(lm_calc.cdrsum, "likelihood", var= "time_since_first")
trends.emm.diff=pairs(trends.emm,data=df)

emmeans_table <- emmeans(
  lm_calc.cdrsum,
  ~ likelihood,
  at = list(time_since_first = 0)
)

pairwise_comparisons <- pairs(emmeans_table)

### Analysis of CDR-SB trajectories by LANS X Path diagnosis -------------------
df = read.delim('XXX.csv', header = T, sep = ',')
df$likelihood_path[df$likelihood_path == "ADNC (Highest)"] <- "ADNC/LATE-NC (Highest/high)"
df$likelihood_path[df$likelihood_path == "ADNC (High)"] <- "ADNC/LATE-NC (Highest/high)"
df$likelihood_path = factor(df$likelihood_path, levels = c("ADNC", "ADNC/LATE-NC (Highest/high)", "ADNC/LATE-NC (Moderate)", "ADNC/LATE-NC (Low)", "LATE-NC"))

cleanup = theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "white", colour = "white"),
                axis.text =element_text(size=15, face= "bold", color = "black"), 
                axis.title = element_text(size = 18, face = "bold", color = "black"),
                plot.title = element_text(size = 20, face = "bold", color = "black"),
                legend.text = element_text(size=12, face="bold", color = "black"),
                axis.line = element_line(colour = "black"),
                legend.title = element_blank())

p_calc.cdrsum_long = ggplot(data=df, aes(x = time_since_first, y = calc.cdrsum, color = likelihood_path)) +
  geom_line(aes(group = clinic), alpha = 0.3) +
  geom_point(aes(fill = likelihood_path), pch = 21, size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", aes(group = likelihood_path, fill = likelihood_path), size = 2, alpha = 0.3) +
  scale_fill_manual(values=c('#997700', '#882255', '#0077BB', '#117733', '#332288')) +
  scale_color_manual(values=c('#997700', '#882255', '#0077BB', '#117733', '#332288')) +
  labs(y = "CDR-SB", x = "Time since baseline (years)") + 
  cleanup  +   theme(legend.position = "top") +
  xlim(0,13) +
  theme(legend.position = c(0.5, 0.95), legend.direction = "horizontal") +
  guides(fill = guide_legend(override.aes = list(color = NULL), nrow = 2)) # Remove color legend

lm_calc.cdrsum = lmer(calc.cdrsum ~ likelihood_path*time_since_first + (1 | clinic), data = df)
summary(lm_calc.cdrsum)

trends.emm = emtrends(lm_calc.cdrsum, "likelihood_path", var= "time_since_first")
trends.emm.diff=pairs(trends.emm,data=df)

emmeans_table <- emmeans(
  lm_calc.cdrsum,
  ~ likelihood_path,
  at = list(time_since_first = 0)
)

pairwise_comparisons <- pairs(emmeans_table)

#Pie chart of pathological diagnoses (Figure 1) --------------------------------
df = read.delim('XXX.csv', header = T, sep = ',')

df_baseline = df %>%
  group_by(clinic) %>%
  slice_min(visitage) %>%
  ungroup()

cleanup = theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "white", colour = "white"),
                axis.text.y =element_text(size=12, face= "bold", color = "black"), 
                axis.text.x =element_text(size=10, face= "bold", color = "black"), 
                axis.title = element_blank(),
                plot.title = element_text(size = 15, face = "bold", color = "black"),
                legend.text = element_text(size=10, face="bold", color = "black"),
                legend.title = element_blank())

path_dx_freq <- table(df_baseline$Path_dx)
path_dx_data <- data.frame(Path_dx = names(path_dx_freq), Frequency = as.vector(path_dx_freq))
path_dx_data <- path_dx_data[order(path_dx_data$Frequency, decreasing = TRUE), ]
path_dx_data$Percentage <- path_dx_data$Frequency / sum(path_dx_data$Frequency) * 100

pie_chart = ggplot(path_dx_data, aes(x="", y=Frequency, fill=Path_dx)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values = c("#1f77b4", '#ff7f0e', '#2ca02c', '#9467bd',
                               '#8c564b',  '#e377c2', '#d62728','#bcbd22', 
                               '#17becf', '#aec7e8',  '#ffbb78', '#98df8a',
                               '#7f7f7f')) +
  theme(legend.position = "none")